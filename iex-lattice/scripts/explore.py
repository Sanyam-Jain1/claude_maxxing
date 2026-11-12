import sys; sys.path.insert(0, 'scripts')
from load import *
d = load()
print(d.groupby('mkt').agg(n=('p','size'), d0=('date','min'), d1=('date','max'), days=('date','nunique'), nanp=('p', lambda s: s.isna().sum()), pmin=('p','min'), pmax=('p','max')))
print(d.groupby(['mkt','date']).size().value_counts().head())
for m in ['DAM','GDAM','RTM']:
    for col in ['p','buy','sell','v']:
        w = wide(d, m, col)
        if col != 'p': w = np.log(w)
        dif = w.diff(axis=1).iloc[:, 1:]            # dif[b] = x[b]-x[b-1]; transition INTO block b
        pos = (np.arange(1,96) % 4)                 # position of destination block in hour: 0 = hour boundary
        out = []
        for q in range(4):
            x = dif.loc[:, pos == q].values.ravel(); x = x[~np.isnan(x)]
            out.append((np.mean(np.abs(x)), np.sqrt(np.mean(x**2)), np.mean(x == 0)))
        print(m, col, ' into-q0(hr bdry) q1 q2(half-hr bdry) q3  mean|d|:', ' '.join(f'{o[0]:.4f}' for o in out), '| rms:', ' '.join(f'{o[1]:.4f}' for o in out), '| zero-frac:', ' '.join(f'{o[2]:.3f}' for o in out))
