# Exact check (rational arithmetic): for source points S on a rational grid, the polygon cut out by the six
# inequalities (i)-(vi) is exactly the set of targets for which the four-face unfolding D1 is the unique shortest.
# D1 - D_k is affine in T for fixed S, so checking polygon vertices suffices.
from fractions import Fraction as Fr
from chambers import polygon
def others(x, y, u, v):
    a, b, c, d = 1-x, 1-y, 1-v, 1-u
    D1 = (3-d-a)**2+(2-c-b)**2
    rest = [(b-d)**2+(3-a-c)**2, (b-d)**2+(1+a+c)**2, (a-c)**2+(3-b-d)**2, (a-c)**2+(1+b+d)**2,
            (3-c-b)**2+(2-d-a)**2, (2-d+a)**2+(1+c-b)**2, (2+c-b)**2+(1-d+a)**2,
            (2+d-a)**2+(1-c+b)**2, (2-c+b)**2+(1+d-a)**2, (1+d+a)**2+(c+b)**2, (1+c+b)**2+(d+a)**2]
    return D1, rest
M = 150; bad = 0; npoly = 0
for i in range(M):
    for j in range(M):
        x, y = Fr(2*i+1, 2*M), Fr(2*j+1, 2*M)
        p = polygon(x, y)
        if not p: continue
        npoly += 1
        for (u, v), _ in p:
            D1, rest = others(x, y, u, v)
            if any(r < D1 for r in rest): bad += 1
print(f"grid {M}x{M}: {npoly} non-empty polygons, vertices where some other unfolding is strictly shorter than D1: {bad}")
