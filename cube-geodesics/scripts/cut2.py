import sympy as sp, mpmath as mp, pickle, itertools, sys
from fractions import Fraction as Fr
x, y, u, v = sp.symbols('x y u v')
ar = pickle.load(open('areas.pkl', 'rb'))
pos = {s.name: s for e in ar.values() for s in e.free_symbols}
T = {k: e.subs({pos['x']: x, pos['y']: y}) for k, e in ar.items()}
# triangle above v=1 when apex is outside the box
Li = (1+x+y)*u-(1+x-y)*v; Lii = -(y-x)*u+(x+y)*v-(y-x)
A = sp.solve([Li, Lii], [u, v], dict=True)[0]
p1 = sp.solve(Li.subs(v, 1), u)[0]; p2 = sp.solve(Lii.subs(v, 1), u)[0]
T['Tb'] = sp.factor(sp.Rational(1, 2)*(p1-p2)*(A[v]-1))
print("Tbox =", T['Tb'])
N3 = x**3 - x**2*y - 2*x**2 + x*y**2 - 2*x*y - 3*x - y**3 + y
N5 = 2*x**3 - x**2 + 2*x*y**2 - 3*x - y**2 + y
Nc = x**2 + y**2 + 2*y - 1
N15 = 3*x**3 + x**2*y - 6*x**2 + 3*x*y**2 + 2*x*y + 3*x + y**3 - y
Nb = y - y**2 - 3*x**2 - 3*x
polys = [N3, N5, Nc, N15, Nb]
fT = {k: sp.lambdify((x, y), e, 'mpmath') for k, e in T.items()}
fN = [sp.lambdify((x, y), e, 'mpmath') for e in polys]
def cut(xv, yv):
    n3, n5, nc, n15, nb = (f(xv, yv) for f in fN)
    if nc > 0 and n15 < 0: c = fT['T1'](xv, yv)
    else:
        c = 0
        if n3 > 0: c = max(c, fT['T3'](xv, yv))
        if n5 > 0: c = max(c, fT['T5'](xv, yv))
    if nb > 0: c -= fT['Tb'](xv, yv)
    return c
def real_roots_01(P, var, dps):
    out = []
    P = sp.Poly(P, var)
    if P.degree() < 1: return out
    for fac, _ in sp.factor_list(P.as_expr())[1]:
        fp = sp.Poly(fac, var)
        if fp.degree() < 1: continue
        for r in fp.nroots(n=dps+10, maxsteps=2000):
            if abs(sp.im(r)) < sp.Float(10)**(-dps) and 0 < sp.re(r) < 1: out.append(mp.mpf(str(sp.re(r))))
    return out
if __name__ == "__main__":
    from chambers import cons, clip, area
    def exact_cut(xf, yf):
        C = cons(xf, yf)
        p = [((Fr(0), Fr(0)), 'v0'), ((Fr(1), Fr(0)), 'u1'), ((Fr(1), Fr(1)), 'v1'), ((Fr(0), Fr(1)), 'u0')]
        for nm in ['i', 'ii']: p = clip(p, (C[nm], nm))
        a0 = area(p) if p else 0
        for nm in ['iii', 'v']:
            if p: p = clip(p, (C[nm], nm))
        return a0 - (area(p) if p else 0)
    mp.mp.dps = 30; bad = 0; M = 90
    for i in range(M):
        for j in range(i+1, M):
            xf, yf = Fr(2*i+1, 2*M), Fr(2*j+1, 2*M)
            e = exact_cut(xf, yf); c = cut(mp.mpf(xf.numerator)/xf.denominator, mp.mpf(yf.numerator)/yf.denominator)
            if abs(mp.mpf(e.numerator)/e.denominator - c) > mp.mpf(10)**-20: bad += 1
    print("piecewise cut formula mismatches on grid:", bad)
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    mp.mp.dps = dps
    cand = []
    for P in polys:
        cand += [sp.discriminant(P, x), P.subs(x, 0), P.subs(x, y)]
    for P, Q in itertools.combinations(polys, 2): cand.append(sp.resultant(P, Q, x))
    ycrit = [mp.mpf(0), mp.mpf(1)]
    for c in cand: ycrit += real_roots_01(c, y, dps)
    ycrit = sorted(ycrit); yc = [ycrit[0]]
    for c in ycrit[1:]:
        if c - yc[-1] > mp.mpf(10)**(-dps+5): yc.append(c)
    print("critical y:", [mp.nstr(c, 12) for c in yc])
    cx = [[sp.lambdify(y, co, 'mpmath') for co in sp.Poly(P, x).all_coeffs()] for P in polys]
    def inner(yv):
        br = [mp.mpf(0), yv]
        for cs in cx:
            for r in mp.polyroots([c(yv) for c in cs], maxsteps=500, extraprec=300):
                if abs(mp.im(r)) < mp.mpf(10)**(-dps+8) and 0 < mp.re(r) < yv: br.append(mp.re(r))
        br = sorted(br); tot = mp.mpf(0)
        for a, b in zip(br, br[1:]):
            if b - a < mp.mpf(10)**(-dps+6): continue
            if cut((a+b)/2, yv) == 0 and cut(a+(b-a)/7, yv) == 0: continue
            tot += mp.quad(lambda xv: cut(xv, yv), [a, b])
        return tot
    Vcut = mp.quad(inner, yc)
    print("Vcut =", Vcut)
    V0x8 = mp.mpf('0.159804142167423584938385007215218223528213281995088507800313')
    print("P4 = 8V0 - 16 Vcut =", V0x8 - 16*Vcut)
