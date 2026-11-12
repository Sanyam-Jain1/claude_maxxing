# Rigorous finite check: every unfolding (simple face sequence) other than the listed candidates is dominated
# everywhere by a single candidate. D_k - D_j is bi-affine in (S,T), so it suffices to check the 16 vertex pairs.
import numpy as np, itertools
from fractions import Fraction as Fr
from brute import sequences, unfold_maps, pt
def sqdist_fn(seq):
    maps, _ = unfold_maps(list(seq))
    R, t = maps[-1]; R = np.rint(R).astype(int); t = np.rint(t).astype(int)
    f0, fk = seq[0], seq[-1]
    def D(su, sv, tu, tv):
        S = [Fr(0)]*3; S[f0[0]] = Fr(f0[1]); o = [i for i in range(3) if i != f0[0]]; S[o[0]] = su; S[o[1]] = sv
        T = [Fr(0)]*3; T[fk[0]] = Fr(fk[1]); o = [i for i in range(3) if i != fk[0]]; T[o[0]] = tu; T[o[1]] = tv
        Ti = [sum(int(R[i][j])*T[j] for j in range(3))+int(t[i]) for i in range(3)]
        return sum((S[i]-Ti[i])**2 for i in range(3))
    return D
corners = [tuple(Fr(v) for v in c) for c in itertools.product((0, 1), repeat=4)]
def analyse(f0, fk, ncand_len):
    seqs = sequences(f0, fk, 6)
    fns = {s: sqdist_fn(s) for s in seqs}
    cands = [s for s in seqs if len(s) <= ncand_len]
    others = [s for s in seqs if len(s) > ncand_len]
    print(f"faces {f0}->{fk}: {len(seqs)} simple face sequences; {len(cands)} candidates (<= {ncand_len} faces), {len(others)} longer")
    vals = {s: [fns[s](*c) for c in corners] for s in seqs}
    undominated = []
    for s in others:
        if not any(all(a >= b for a, b in zip(vals[s], vals[c])) for c in cands): undominated.append(s)
    print("   longer sequences NOT dominated by a single candidate:", len(undominated))
    return undominated
bot, top, right = (2, 0), (2, 1), (0, 1)
u1 = analyse(bot, right, 3)
u2 = analyse(bot, top, 4)
for s in (u1+u2)[:10]: print("   ", s)
