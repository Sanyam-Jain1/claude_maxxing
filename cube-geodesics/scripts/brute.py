# Independent brute-force geodesic distance on the unit cube surface by enumerating
# all face-sequence unfoldings (with validity check), to validate closed-form candidates.
import numpy as np, itertools, sys
rng = np.random.default_rng(1)
# faces: (axis, value)
faces = [(ax, v) for ax in range(3) for v in (0, 1)]
def adjacent(f, g): return f[0] != g[0]
def edge(f, g):
    # shared edge endpoints of adjacent faces f,g
    ax = 3 - f[0] - g[0]
    p = np.zeros(3); p[f[0]] = f[1]; p[g[0]] = g[1]
    q = p.copy(); q[ax] = 1.0
    return p, q
def center(f):
    c = np.full(3, .5); c[f[0]] = f[1]; return c
def rot_about(p, q, ang):
    # returns function x -> rotation of x about line pq by ang
    k = (q - p) / np.linalg.norm(q - p)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
    return R, p - R @ p   # x -> R x + t
def unfold_maps(seq):
    # returns list of (R,t) mapping face seq[i] into plane of seq[0], plus edge images
    maps = [(np.eye(3), np.zeros(3))]
    edges = []
    for i in range(len(seq) - 1):
        f, g = seq[i], seq[i + 1]
        p, q = edge(f, g)
        best = None
        for ang in (np.pi / 2, -np.pi / 2):
            R, t = rot_about(p, q, ang)
            c = R @ center(g) + t
            # must lie in plane of f and not coincide with center(f)
            if abs(c[f[0]] - f[1]) < 1e-9 and np.linalg.norm(c - center(f)) > 1e-6:
                best = (R, t)
        R, t = best
        Rp, tp = maps[-1]
        maps.append((Rp @ R, Rp @ t + tp))
        edges.append((Rp @ p + tp, Rp @ q + tp))
    return maps, edges
def seg_cross(a, b, p, q, drop):
    # do segments ab and pq intersect (2D after dropping coordinate 'drop')?
    idx = [i for i in range(3) if i != drop]
    a, b, p, q = a[idx], b[idx], p[idx], q[idx]
    def cr(u, v): return u[0] * v[1] - u[1] * v[0]
    d1 = cr(b - a, p - a); d2 = cr(b - a, q - a); d3 = cr(q - p, a - p); d4 = cr(q - p, b - p)
    return d1 * d2 <= 1e-12 and d3 * d4 <= 1e-12
def sequences(f0, fk, maxlen):
    out = []
    def rec(seq):
        if seq[-1] == fk and len(seq) > 1 or (len(seq) == 1 and f0 == fk):
            out.append(tuple(seq))
        if len(seq) >= maxlen: return
        for g in faces:
            if adjacent(seq[-1], g) and g not in seq:
                rec(seq + [g])
    rec([f0])
    return out
cache = {}
def geodesic(S, f0, T, fk, maxlen=6, need_valid=True):
    key = (f0, fk, maxlen)
    if key not in cache:
        cache[key] = [(s,) + unfold_maps(s) for s in sequences(f0, fk, maxlen)]
    best = (1e9, None)
    for s, maps, edges in cache[key]:
        R, t = maps[-1]
        Ti = R @ T + t
        d = np.linalg.norm(Ti - S)
        if d < best[0]:
            ok = all(seg_cross(S, Ti, p, q, f0[0]) for p, q in edges) if need_valid else True
            if ok: best = (d, s)
    return best
def pt(f, u, v):
    x = np.zeros(3); x[f[0]] = f[1]
    o = [i for i in range(3) if i != f[0]]
    x[o[0]] = u; x[o[1]] = v
    return x
if __name__ == "__main__":
    bot, top, right = (2, 0), (2, 1), (0, 1)
    # opposite faces: closed-form 12 candidates
    def opp(a, b, c, d):
        c3 = [(b-d)**2+(3-a-c)**2, (b-d)**2+(1+a+c)**2, (a-c)**2+(3-b-d)**2, (a-c)**2+(1+b+d)**2]
        c4 = [(3-d-a)**2+(2-c-b)**2, (3-c-b)**2+(2-d-a)**2, (2-d+a)**2+(1+c-b)**2, (2+c-b)**2+(1-d+a)**2,
              (2+d-a)**2+(1-c+b)**2, (2-c+b)**2+(1+d-a)**2, (1+d+a)**2+(c+b)**2, (1+c+b)**2+(d+a)**2]
        return min(c3), min(c4)
    def adj(a, b, y, z):
        return [(1+z-a)**2+(y-b)**2, (2-y-a)**2+(1+z-b)**2, (1+y-a)**2+(z+b)**2]
    mx = 0; n4 = 0; nlen = {}
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    for _ in range(N):
        a, b, c, d = rng.random(4)
        g, s = geodesic(pt(bot, a, b), bot, pt(top, c, d), top)
        g2, _ = geodesic(pt(bot, a, b), bot, pt(top, c, d), top, need_valid=False)
        m3, m4 = opp(a, b, c, d)
        mx = max(mx, abs(g - np.sqrt(min(m3, m4))), abs(g - g2)); n4 += m4 < m3
        nlen[len(s)] = nlen.get(len(s), 0) + 1
    print("opposite: max |brute - formula| =", mx, " frac 4-face:", n4 / N, nlen)
    mx = 0; n3 = 0; nlen = {}
    for _ in range(N):
        a, b, y, z = rng.random(4)
        g, s = geodesic(pt(bot, a, b), bot, pt(right, y, z), right)
        g2, _ = geodesic(pt(bot, a, b), bot, pt(right, y, z), right, need_valid=False)
        cand = adj(a, b, y, z)
        mx = max(mx, abs(g - np.sqrt(min(cand))), abs(g - g2)); n3 += min(cand[1:]) < cand[0]
        nlen[len(s)] = nlen.get(len(s), 0) + 1
    print("adjacent: max |brute - formula| =", mx, " frac via third face:", n3 / N, nlen)
