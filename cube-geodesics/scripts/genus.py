import sympy as sp
x, y, z = sp.symbols('x y z')
curves = {'c_iii': x**3 - x**2*y - 2*x**2 + x*y**2 - 2*x*y - 3*x - y**3 + y,
          'c_v': 2*x**3 - x**2 + 2*x*y**2 - 3*x - y**2 + y,
          'c_15': 3*x**3 + x**2*y - 6*x**2 + 3*x*y**2 + 2*x*y + 3*x + y**3 - y,
          'c_iv': x**3 - 3*x**2*y + 3*x**2 + x*y**2 - 4*x*y + 2*x - 3*y**3 + 3*y**2,
          'c_vi': 2*x**2*y - 3*x**2 + 4*x*y - 3*x + 2*y**3 - 3*y**2 + y}
for nm, F in curves.items():
    Fh = sp.Poly(F, x, y).homogenize(z).as_expr()
    G = sp.groebner([Fh, sp.diff(Fh, x), sp.diff(Fh, y), sp.diff(Fh, z)], x, y, z, order='grevlex')
    # singular points: solve in affine charts
    sing = []
    for chart in ({z: 1}, {z: 0, y: 1}, {z: 0, y: 0, x: 1}):
        eqs = [e.subs(chart) for e in (Fh, sp.diff(Fh, x), sp.diff(Fh, y), sp.diff(Fh, z))]
        vars_ = [s for s in (x, y, z) if s not in chart]
        if vars_:
            sol = sp.solve(eqs, vars_, dict=True)
            sing += [(chart, s) for s in sol]
        else:
            if all(sp.simplify(e) == 0 for e in eqs): sing.append((chart, {}))
    print(nm, "factor:", sp.factor(F), "| singular points:", sing)
