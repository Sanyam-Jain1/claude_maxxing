"""More realistic dispatch comparisons for a discom-contracted 1 MW / 2 MWh BESS (85% RTE) on IEX DAM prices.
 - opt2 / opt1: perfect knowledge of that day's DAM prices (available ~D-1 afternoon), 2 or 1 cycles/day
 - clock: fixed windows all year (charge 11-13h, discharge 19-21h; second cycle 02-04h -> 06-08h)
 - seasonal clock: best fixed windows re-chosen each month from the previous month's average price shape
 - persistence: optimise on yesterday's prices, execute on today's"""
import sys; sys.path.insert(0, '../iex-lattice/scripts')
from load import *
from scipy.optimize import linprog
import scipy.sparse as sp
d = load(); P = wide(d,'DAM','p').dropna()*1000
eta = 0.85**0.5; E = 2.0; dt = 0.25; n = 96
A = sp.lil_matrix((n, 3*n))
for t in range(n):
    A[t, 2*n+t] = 1; A[t, t] = -eta*dt; A[t, n+t] = dt/eta
    if t: A[t, 2*n+t-1] = -1
A = A.tocsr(); Aub = sp.lil_matrix((1, 3*n)); Aub[0, n:2*n] = dt; Aub = Aub.tocsr()
def solve(p, cyc):
    r = linprog(np.r_[p*dt, -p*dt, np.zeros(n)], A_eq=A, b_eq=np.zeros(n), A_ub=Aub, b_ub=[cyc*E*eta],
                bounds=[(0,1)]*(2*n)+[(0,E)]*n, method='highs')
    return r.x[:n], r.x[n:2*n]
val = lambda p, c, dch: float(((dch - c)*p).sum()*dt)
days = P.index; X = P.values
out = {k: [] for k in ['opt2','opt1','clock','season','persist']}
clk_c = np.zeros(n); clk_d = np.zeros(n); clk_c[44:52] = 1; clk_d[76:84] = eta*eta; clk_c[8:16] = 1; clk_d[24:32] = eta*eta
sched_prev = None; season = {}
for i, p in enumerate(X):
    c2, d2 = solve(p, 2); c1, d1 = solve(p, 1)
    out['opt2'].append(val(p, c2, d2)); out['opt1'].append(val(p, c1, d1)); out['clock'].append(val(p, clk_c, clk_d))
    out['persist'].append(val(p, *sched_prev) if sched_prev else np.nan); sched_prev = (c2, d2)
    key = (days[i].year, days[i].month)
    if key not in season:
        prev = P[(P.index < pd.Timestamp(key[0], key[1], 1)) & (P.index >= pd.Timestamp(key[0], key[1], 1) - pd.Timedelta(days=30))]
        season[key] = solve(prev.mean().values, 2) if len(prev) > 10 else None
    out['season'].append(val(p, *season[key]) if season[key] else np.nan)
R = pd.DataFrame(out, index=days); k = 30.4/1e5
m = R.dropna()
print('Rs lakh / MW / month (mean over %d days, Feb 2023 - Sep 2025):' % len(m)); print((m.mean()*k).round(2).to_string())
print('share of opt2:', (m.sum()/m.opt2.sum()).round(2).to_dict())
print('by year:'); print((m.groupby(m.index.year).mean()*k).round(2))
print('gap opt2 - seasonal clock, Rs lakh/MW/yr:', round(float((m.opt2-m.season).mean()*365/1e5),1), '; opt2 - persistence:', round(float((m.opt2-m.persist).mean()*365/1e5),1))
