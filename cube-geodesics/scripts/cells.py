# Semi-analytic integration: for a fixed source S, the target face is partitioned into cells
# (nearest unfolded image of S). Cell areas and integrals of distance are computed in closed form;
# the outer integral over S uses composite Gauss-Legendre.
import math, sys
from multiprocessing import Pool
import numpy as np
def clip(poly, a, b, c):
    # keep a*x+b*y+c <= 0
    out = []; n = len(poly)
    for k in range(n):
        p = poly[k]; q = poly[(k+1) % n]
        fp = a*p[0]+b*p[1]+c; fq = a*q[0]+b*q[1]+c
        if fp <= 0: out.append(p)
        if (fp < 0 and fq > 0) or (fp > 0 and fq < 0):
            t = fp/(fp-fq); out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    return out
SQ = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
def cell_stats(sites):
    # returns list of (area, integral of distance to own site) per site
    res = []
    for i, (px, py) in enumerate(sites):
        poly = SQ
        for j, (qx, qy) in enumerate(sites):
            if j == i: continue
            # closer to p than q: 2(q-p).T <= |q|^2-|p|^2
            poly = clip(poly, 2*(qx-px), 2*(qy-py), px*px+py*py-qx*qx-qy*qy)
            if len(poly) < 3: break
        if len(poly) < 3: res.append((0., 0.)); continue
        A = 0.; I = 0.; n = len(poly)
        for k in range(n):
            ax, ay = poly[k]; bx, by = poly[(k+1) % n]
            A += ax*by-ay*bx
            ex, ey = bx-ax, by-ay; L = math.hypot(ex, ey)
            if L < 1e-300: continue
            ex /= L; ey /= L
            ta = (ax-px)*ex+(ay-py)*ey; tb = ta+L
            h = (ax-px)*ey-(ay-py)*ex          # signed distance (CCW polygon, p inside -> h>0)
            h2 = h*h
            if h2 < 1e-300: continue
            ah = abs(h)
            F = lambda t: 0.5*(t*math.sqrt(h2+t*t)+h2*math.asinh(t/ah))
            I += h/3.0*(F(tb)-F(ta))
        res.append((A/2, I))
    return res
def sites_opp(a, b):
    # images of S=(a,b) (bottom face) in the frame of the top face (c,d); first 4 = three-face paths
    return [(3-a, b), (-1-a, b), (a, 3-b), (a, -1-b),
            (2-b, 3-a), (3-b, 2-a), (b-1, 2+a), (b-2, 1+a), (1+b, a-2), (2+b, a-1), (-b, -1-a), (-1-b, -a)]
def sites_adj(a, b):
    # S=(a,b) on bottom, T=(z,y) on face x=1: direct, via y=1 face, via y=0 face
    return [(a-1, b), (b-1, 2-a), (-b, a-1)]
gl = {}
def panel_job(args):
    kind, x0, x1, y0, y1, g, tri = args
    xs, ws = np.polynomial.legendre.leggauss(g)
    out = np.zeros(4)
    for i in range(g):
        for j in range(g):
            if tri:   # Duffy map of the triangle (x0,y0),(x1,y0),(x1,y1): a in [x0,x1], b in [y0, y0+(a-x0)*(y1-y0)/(x1-x0)]
                a = 0.5*(x0+x1)+0.5*(x1-x0)*xs[i]; hgt = (a-x0)*(y1-y0)/(x1-x0)
                b = y0+hgt*(0.5+0.5*xs[j]); w = ws[i]*ws[j]*0.25*(x1-x0)*hgt
            else:
                a = 0.5*(x0+x1)+0.5*(x1-x0)*xs[i]; b = 0.5*(y0+y1)+0.5*(y1-y0)*xs[j]; w = ws[i]*ws[j]*0.25*(x1-x0)*(y1-y0)
            if kind == 'opp':
                st = cell_stats(sites_opp(a, b))
                out += w*np.array([sum(s[0] for s in st[:4]), sum(s[0] for s in st[4:]), sum(s[1] for s in st[:4]), sum(s[1] for s in st[4:])])
            else:
                st = cell_stats(sites_adj(a, b))
                out += w*np.array([st[0][0], st[1][0]+st[2][0], st[0][1], st[1][1]+st[2][1]])
    return out
def run(kind, K, g):
    jobs = []
    if kind == 'opp':   # fundamental domain 0<=b<=a<=1/2 (x8)
        h = 0.5/K
        for i in range(K):
            for j in range(i+1):
                if j < i: jobs.append((kind, i*h, (i+1)*h, j*h, (j+1)*h, g, False))
                else: jobs.append((kind, i*h, (i+1)*h, j*h, (j+1)*h, g, True))
        mult = 8
    else:               # symmetry b -> 1-b: integrate b in [0,1/2] (x2)
        h = 1.0/K
        for i in range(K):
            for j in range(K//2): jobs.append((kind, i*h, (i+1)*h, j*h, (j+1)*h, g, False))
        mult = 2
    with Pool(4) as p: res = p.map(panel_job, jobs, chunksize=64)
    return mult*np.sum(res, axis=0)
if __name__ == "__main__":
    kind = sys.argv[1]; K = int(sys.argv[2]); g = int(sys.argv[3])
    r = run(kind, K, g)
    np.set_printoptions(precision=15)
    print(kind, "K=%d g=%d" % (K, g), "area_main, area_other, Idist_main, Idist_other =", r)
    print("  total area check:", r[0]+r[1], " P(other) =", r[1], " mean distance =", r[2]+r[3])
