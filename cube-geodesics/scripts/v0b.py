import sympy as sp, mpmath as mp
mp.mp.dps = 60
s, t = sp.symbols('s t', positive=True)
p, q = 1-s, 1+s
D1 = p+q*t; D2 = q-p*t
f = p**2*(1-t)**2/(D1**2*D2**2)
# manual partial fractions in t: A/D1^2 + B/D1 + C/D2^2 + E/D2
A, B, C, E = sp.symbols('A B C E')
num = sp.expand(A*D2**2 + B*D1*D2**2 + C*D1**2 + E*D2*D1**2 - p**2*(1-t)**2)
sol = sp.solve(sp.Poly(num, t).coeffs(), [A, B, C, E], dict=True)[0]
sol = {k: sp.factor(v) for k, v in sol.items()}
print(sol)
# antiderivative: -A/(q D1) + (B/q) log D1 + C/(p D2) - (E/p) log D2
F = -sol[A]/(q*D1) + sol[B]/q*sp.log(D1) + sol[C]/(p*D2) - sol[E]/p*sp.log(D2)
assert sp.simplify(sp.diff(F, t) - f) == 0
lo = s/(2+s); hi1 = 2*s/(1-s)
G1 = sp.simplify(F.subs(t, hi1) - F.subs(t, lo))   # s < 1/3
G2 = sp.simplify(F.subs(t, 1) - F.subs(t, lo))     # s > 1/3
print("G1 =", G1); print("G2 =", G2)
g1 = sp.lambdify(s, G1, 'mpmath'); g2 = sp.lambdify(s, G2, 'mpmath')
IB = mp.quad(g1, [0, mp.mpf(1)/3]) + mp.quad(g2, [mp.mpf(1)/3, 1])
print("I_B =", IB)
IA = mp.mpf(1)/2 - 2*mp.log(2) + mp.pi*mp.sqrt(3)/6
print("8V0 =", 2*(2*IA+IB))
G = mp.catalan
basis = {'1': mp.mpf(1), 'pi*sqrt3': mp.pi*mp.sqrt(3), 'log2': mp.log(2), 'log3': mp.log(3), 'pi': mp.pi,
         'pi^2': mp.pi**2, 'log2^2': mp.log(2)**2, 'log3^2': mp.log(3)**2, 'log2log3': mp.log(2)*mp.log(3),
         'Catalan': G, 'pi*log2': mp.pi*mp.log(2), 'pi*log3': mp.pi*mp.log(3), 'Li2(1/3)': mp.polylog(2, mp.mpf(1)/3),
         'Li2(1/4)': mp.polylog(2, mp.mpf(1)/4), 'log5': mp.log(5), 'atan(1/2)': mp.atan(mp.mpf(1)/2), 'atan(1/3)': mp.atan(mp.mpf(1)/3)}
names = list(basis)
rel = mp.pslq([IB] + [basis[n] for n in names], maxcoeff=10**8, maxsteps=10**6)
print("PSLQ:", rel)
if rel: print({('I_B' if i == 0 else names[i-1]): r for i, r in enumerate(rel) if r})
