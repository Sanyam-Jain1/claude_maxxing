# Mean geodesic distance between random points on two ADJACENT faces of the unit cube, to high precision.
import mpmath as mp
mp.mp.dps = 40
def g(m, c):
    # int_0^m lam*sqrt(1+lam^2-2*lam*c) dlam (closed form)
    k = 1-c*c
    def F(l):
        q = 1+l*l-2*l*c
        if q < 0: q = mp.mpf(0)
        a = (k/2*mp.asinh((l-c)/mp.sqrt(k))) if k > mp.mpf(10)**(-mp.mp.dps+2) else 0
        return q**mp.mpf(1.5)/3 + c*((l-c)*mp.sqrt(q)/2 + a)
    return F(m)-F(mp.mpf(0))
def J(c, R1, R2): return R1**5/5*g(R2/R1, c) + R2**5/5*g(R1/R2, c)
# sanity check of J against direct quadrature
c0, R1, R2 = mp.mpf('-0.3'), mp.mpf('1.2'), mp.mpf('0.7')
chk = mp.quad(lambda r1: mp.quad(lambda r2: r1*r2*mp.sqrt(r1*r1+r2*r2-2*r1*r2*c0), [0, R2]), [0, R1])
print("J check:", J(c0, R1, R2), chk)
def integrand(t1, t2):
    D = t1+t2; R1 = 1/mp.sin(t1); R2 = 1/mp.sin(t2)
    return J(mp.cos(D), R1, R2) - J(mp.cos(3*mp.pi/2-D), R1, R2)
corr = mp.quad(lambda t1: mp.quad(lambda t2: integrand(t1, t2), [3*mp.pi/4-t1, mp.pi/2]), [mp.pi/4, mp.pi/2])
print("E[(d0-d1)^+] =", corr)
def rect(a, b):
    d = mp.sqrt(a*a+b*b)
    return (a**3/b**2 + b**3/a**2 + d*(3-a*a/(b*b)-b*b/(a*a)) + mp.mpf(5)/2*(b*b/a*mp.log((a+d)/b) + a*a/b*mp.log((b+d)/a)))/15
Esq = rect(mp.mpf(1), mp.mpf(1)); E12 = rect(mp.mpf(1), mp.mpf(2)); Ed0 = 2*E12-Esq
print("E same face =", Esq); print("E[d0] (flat adjacent squares) =", Ed0)
Dadj = Ed0 - 2*corr
print("D_adj =", Dadj)
open('adjmean.txt', 'w').write(mp.nstr(corr, 38)+"\n"+mp.nstr(Dadj, 38)+"\n")
s2, s5 = mp.sqrt(2), mp.sqrt(5)
basis = {'1': mp.mpf(1), 'sqrt2': s2, 'log(1+sqrt2)': mp.log(1+s2), 'log2': mp.log(2), 'pi': mp.pi, 'sqrt2*log2': s2*mp.log(2),
         'sqrt2*pi': s2*mp.pi, 'sqrt2*log(1+sqrt2)': s2*mp.log(1+s2), 'pi*log2': mp.pi*mp.log(2), 'Catalan': mp.catalan,
         'log2^2': mp.log(2)**2, 'pi^2': mp.pi**2, 'log(1+sqrt2)^2': mp.log(1+s2)**2, 'log2*log(1+sqrt2)': mp.log(2)*mp.log(1+s2), 'pi*log(1+sqrt2)': mp.pi*mp.log(1+s2)}
names = list(basis)
for sub in (names[:5], names[:8], names[:11], names):
    rel = mp.pslq([corr]+[basis[n] for n in sub], maxcoeff=10**7, maxsteps=10**6)
    print(len(sub), "PSLQ:", None if not rel else {('corr' if i == 0 else sub[i-1]): r for i, r in enumerate(rel) if r})
