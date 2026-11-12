import sys; sys.path.insert(0, 'scripts')
from load import *
d = load()
P = {m: wide(d, m, 'p') for m in ['DAM','GDAM','RTM']}
idx = P['DAM'].index.intersection(P['RTM'].dropna().index)
np.set_printoptions(linewidth=250, precision=3, suppress=True)
for m in P:
    mp = P[m].loc[idx].mean()
    st = mp.diff().values
    print(m, 'mean signed step into block (rows=hour 0..23, cols = into :00 :15 :30 :45)')
    print(np.round(np.r_[np.nan, st[1:]].reshape(24,4), 3))
sp = (P['DAM'].loc[idx] - P['RTM'].loc[idx])
print('mean DAM-RTM spread by block (rows=hour; cols :00 :15 :30 :45)')
print(sp.mean().values.reshape(24,4))
