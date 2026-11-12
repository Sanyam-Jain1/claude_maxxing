"""Detrended within-hour shape: for each (day,hour) with hours h-1,h,h+1 available, the 4 blocks of hour h minus
a centred 2-hour moving average (weights 1/8,1/4,1/4,1/4... i.e. the 2x4 MA that kills any hour-periodic pattern
and passes linear trends). Residual pattern r_q, q=0..3, then split by sign of local trend (mean(h+1)-mean(h-1))."""
import sys; sys.path.insert(0, 'scripts')
from load import *
d = load()
np.set_printoptions(linewidth=250, precision=3, suppress=True)
def resid(w):
    x = w.values.astype(float)                       # days x 96
    k = np.array([.5,1,1,1,.5])/4                    # centred 1-hour (2x4) MA
    ma = np.full_like(x, np.nan)
    for j in range(2, 94): ma[:, j] = x[:, j-2:j+3] @ k
    return x - ma, x
def report(m, col, log=False):
    w = wide(d, m, col)
    if log: w = np.log(w.where(w > 0))
    r, x = resid(w)
    hm = np.nanmean(x.reshape(len(x), 24, 4), axis=2)             # hourly means
    tr = np.full_like(hm, np.nan); tr[:, 1:-1] = (hm[:, 2:] - hm[:, :-2]) / 8   # trend per block
    R = r.reshape(len(x), 24, 4)
    print(f'--- {m} {col}{" (log)" if log else ""}: mean detrended residual by quarter position q0..q3')
    for name, mask in [('all', np.isfinite(tr)), ('rising', tr > 0), ('falling', tr < 0)]:
        mm = np.where(mask[:, :, None], R, np.nan)
        mu = np.nanmean(mm, axis=(0,1)); n = np.sum(np.isfinite(mm[:, :, 0]))
        se = np.nanstd(mm, axis=(0,1)) / np.sqrt(n)
        print(f'  {name:8s} n={n:6d}  mean {mu}  se {se}')
    # slope form: regress (r_q3 - r_q0) on trend  -> >0 means sawtooth (overshoot), <0 means staircase (damped)
    y = (R[:, :, 3] - R[:, :, 0]).ravel(); t = tr.ravel(); ok = np.isfinite(y) & np.isfinite(t)
    b = np.sum(y[ok]*t[ok]) / np.sum(t[ok]**2)
    print(f'  slope of (x_q3-x_q0 residual) on per-block trend: {b:+.3f}   (0 = smooth; -3 = perfect staircase; >0 = sawtooth/overshoot)')
    return R, tr
for m in ['DAM','GDAM','RTM']:
    report(m, 'p')
for m in ['DAM','RTM']:
    for c in ['buy','sell','v']: report(m, c, log=True)
