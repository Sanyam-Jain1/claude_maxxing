from fractions import Fraction as Fr
import sys, collections
def cons(x, y):
    # each (alpha, beta, gamma): alpha*u + beta*v + gamma < 0
    return {
     'i':  (1+x+y, -(1+x-y), 0),
     'ii': (-(y-x), x+y, -(y-x)),
     'iii':(x, 1+y, x-(1+y)),
     'iv': (1+x, -(1-y), -(1-y)),
     'v':  (1+x+y, 3-x+y, 4*x-4),
     'vi': (4+x-y, x+y, x+3*y-4),
     'u0': (-1, 0, 0), 'u1': (1, 0, -1), 'v0': (0, -1, 0), 'v1': (0, 1, -1)}
def clip(poly, line):
    a, b, c = line[0]; name = line[1]
    out = []
    n = len(poly)
    for k in range(n):
        (p, ep), (q, eq) = poly[k], poly[(k+1) % n]   # ep = name of edge from p to next
        fp = a*p[0]+b*p[1]+c; fq = a*q[0]+b*q[1]+c
        if fp < 0:
            if fq < 0: out.append((p, ep))
            elif fq == 0: out.append((p, ep))
            else:
                t = fp/(fp-fq); r = (p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1]))
                out.append((p, ep)); out.append((r, name))
        elif fp == 0:
            if fq < 0: out.append((p, ep))
            else: out.append((p, name))
        else:
            if fq < 0:
                t = fp/(fp-fq); r = (p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1]))
                out.append((r, ep))
    # dedupe consecutive identical points
    res = []
    for pt in out:
        if res and res[-1][0] == pt[0]: res[-1] = pt
        else: res.append(pt)
    if len(res) > 1 and res[0][0] == res[-1][0]: res.pop(0)
    return res if len(res) >= 3 else []
def polygon(x, y):
    C = cons(x, y)
    poly = [((Fr(0), Fr(0)), 'v0'), ((Fr(1), Fr(0)), 'u1'), ((Fr(1), Fr(1)), 'v1'), ((Fr(0), Fr(1)), 'u0')]
    for name in ['i', 'ii', 'iii', 'iv', 'v', 'vi']:
        poly = clip(poly, (C[name], name))
        if not poly: return []
    return poly
def area(poly):
    s = 0
    for k in range(len(poly)):
        p, q = poly[k][0], poly[(k+1) % len(poly)][0]
        s += p[0]*q[1]-p[1]*q[0]
    return s/2
def label(poly):
    names = [e for _, e in poly]
    k = min(range(len(names)), key=lambda i: names[i:]+names[:i])
    return tuple(names[k:]+names[:k])
if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    labs = collections.Counter(); tot = 0; grid = {}
    for i in range(M):
        for j in range(M):
            x = Fr(2*i+1, 2*M); y = Fr(2*j+1, 2*M)
            p = polygon(x, y)
            if p:
                l = label(p); labs[l] += 1; tot += area(p); grid[(i, j)] = l
    print("midpoint-rule 8*V =", float(8*tot/M/M))
    ids = {}
    for l, cnt in labs.most_common(): ids[l] = len(ids); print(ids[l], cnt, l)
    # ascii map
    step = max(1, M//60)
    for j in range(M-1, -1, -step):
        print(''.join(('%x' % ids[grid[(i, j)]]) if (i, j) in grid else '.' for i in range(0, M, step)))
