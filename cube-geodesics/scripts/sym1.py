import sympy as sp
x, y, u, v = sp.symbols('x y u v', positive=True)
L = {'i': (1+x+y)*u-(1+x-y)*v, 'ii': -(y-x)*u+(x+y)*v-(y-x), 'iii': x*u+(1+y)*v+x-(1+y),
     'iv': (1+x)*u-(1-y)*v-(1-y), 'v': (1+x+y)*u+(3-x+y)*v+4*x-4, 'vi': (4+x-y)*u+(x+y)*v+x+3*y-4}
def inter(a, b):
    s = sp.solve([L[a], L[b]], [u, v], dict=True)[0]
    return (sp.factor(s[u]), sp.factor(s[v]))
def ev(name, P): return sp.factor(L[name].subs({u: P[0], v: P[1]}))
def area(pts):
    s = 0
    for k in range(len(pts)):
        p, q = pts[k], pts[(k+1) % len(pts)]
        s += p[0]*q[1]-p[1]*q[0]
    return sp.factor(s/2)
A = inter('i', 'ii'); print("apex A =", A)
for nm in ['iii', 'v', 'iv', 'vi']: print(f"{nm}(A) =", ev(nm, A))
P_i_iii, P_ii_iii, P_i_v, P_ii_v, P_v_iii = inter('i','iii'), inter('ii','iii'), inter('i','v'), inter('ii','v'), inter('v','iii')
print("v∩iii =", P_v_iii)
print("i(v∩iii) =", ev('i', P_v_iii)); print("ii(v∩iii) =", ev('ii', P_v_iii))
T3 = area([A, P_ii_iii, P_i_iii]); T5 = area([A, P_ii_v, P_i_v]); T1 = area([A, P_ii_iii, P_v_iii, P_i_v])
print("T3 =", T3); print("T5 =", T5); print("T1 =", T1)
T0 = area([(0,0), A, (0, (y-x)/(y+x))]); print("T0 (full triangle) =", T0)
import pickle; pickle.dump(dict(T0=T0, T1=T1, T3=T3, T5=T5), open('areas.pkl', 'wb'))
