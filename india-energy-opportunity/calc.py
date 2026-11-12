"""Kill-tests with public IEX 15-min data (Jan 2023 - Sep 2025).
A) How much of a 2h battery's day-ahead arbitrage value does a dumb fixed-clock rule capture vs perfect foresight?
B) G-DAM premium over DAM vs the REC/buyout value."""
import sys; sys.path.insert(0, '../iex-lattice/scripts')
from load import *
from scipy.optimize import linprog
import scipy.sparse as sp
d = load()
P = wide(d,'DAM','p').dropna(); G = wide(d,'GDAM','p').reindex(P.index)
eta = 0.85**0.5; E = 2.0; dt = 0.25; maxcyc = 2.0
def opt(p):
    n = 96  # vars: c[96], dch[96], soc[96]
    cost = np.r_[p*dt, -p*dt, np.zeros(n)]
    # soc_t = soc_{t-1} + eta*c*dt - dch*dt/eta ; soc_{-1}=0, end free>=0
    A = sp.lil_matrix((n, 3*n)); b = np.zeros(n)
    for t in range(n):
        A[t, 2*n+t] = 1; A[t, t] = -eta*dt; A[t, n+t] = dt/eta
        if t: A[t, 2*n+t-1] = -1
    Aub = sp.lil_matrix((1, 3*n)); Aub[0, n:2*n] = dt     # discharge throughput <= maxcyc*E
    r = linprog(cost, A_eq=A.tocsr(), b_eq=b, A_ub=Aub.tocsr(), b_ub=[maxcyc*E*eta],
                bounds=[(0,1)]*(2*n)+[(0,E)]*n, method='highs')
    return -r.fun
def rule(p, cw, dw):
    """charge full power in windows cw, discharge in dw (lists of (start_blk, nblk)); 1MW, 2MWh usable at soc side"""
    rev = 0
    for (cs, cn), (ds, dn) in zip(cw, dw):
        ein = p[cs:cs+cn].sum()*dt                  # 2h charging at 1 MW -> 2 MWh from grid, stores 2*eta
        stored = cn*dt*eta; out = stored*eta        # delivered MWh, spread evenly over discharge window
        rev += out*p[ds:ds+dn].mean() - ein
    return rev
res = []
for day, row in P.iterrows():
    p = row.values*1000   # Rs/MWh
    o = opt(p)
    r1 = rule(p, [(44,8)], [(76,8)])                          # 1 cycle: charge 11-13h, discharge 19-21h
    r2 = r1 + max(0, rule(p, [(8,8)], [(24,8)]))              # + night->morning cycle 02-04h -> 06-08h, skipped if loss (known day-ahead)
    res.append((day, o, r1, r2))
R = pd.DataFrame(res, columns=['day','opt','rule1','rule2']).set_index('day')
R['yr'] = R.index.year
g = R.groupby('yr')[['opt','rule1','rule2']].mean()*30.4/1e5     # Rs lakh per MW per month
g['rule2/opt'] = g.rule2/g.opt
print('A) 1 MW / 2 MWh, 85% RTE, <=2 cycles/day, DAM prices; Rs lakh/MW/month'); print(g.round(3))
print('   whole sample: opt %.3f, fixed-clock 2-cycle %.3f (%.0f%%)' % (R.opt.mean()*30.4/1e5, R.rule2.mean()*30.4/1e5, 100*R.rule2.sum()/R.opt.sum()))
# adaptive rule: use yesterday's optimal windows? simpler: best fixed windows chosen on previous 30 days
sprd = (G - P)
solar = list(range(36,64))
print('B) G-DAM minus DAM, Rs/kWh: all blocks / solar hours 09-16h, by year')
for y in sorted(set(P.index.year)):
    s = sprd[sprd.index.year==y]; print('  ', y, round(float(np.nanmean(s.values)),3), round(float(np.nanmean(s[solar].values)),3), ' share of blocks GDAM>DAM:', round(float(np.nanmean(s.values>0.001)),2))
