"""Top-of-hour (q0) effect in IEX RTM/DAM/GDAM: robustness across two independent public data sets."""
import sys; sys.path.insert(0, 'scripts')
from load import *
np.set_printoptions(linewidth=250, precision=3, suppress=True)
d = load()
def load2(f='/home/user/nandanrajesh/indian-real-time-market-analysis/RTM_Market Snapshot april2022-june2025.xlsx'):
    x = pd.read_excel(f)
    x.columns = ['date','hour','sess','tb','buy','sell','v','fsv','p']
    x['date'] = pd.to_datetime(x['date'], format='%d-%m-%Y', errors='coerce'); x = x.dropna(subset=['date'])
    hh = x['tb'].str.slice(0,2).astype(int); mm = x['tb'].str.slice(3,5).astype(int)
    x['blk'] = hh*4 + mm//15
    for c in ['buy','sell','v','p']: x[c] = pd.to_numeric(x[c], errors='coerce')
    x['p'] = x['p']/1000; x['mkt'] = 'RTM2'
    for c in ['buy','sell','v']: x[c] = x[c]/1000
    return x.drop_duplicates(['date','blk'])
d2 = load2()
print('source2', d2.date.min(), d2.date.max(), d2.date.nunique(), 'sessions ok:', ((d2.sess-1)==d2.blk//2).mean())
# agreement of two sources on overlap
a = wide(d,'RTM','p'); b = d2.pivot(index='date', columns='blk', values='p').reindex(columns=range(96))
ix = a.index.intersection(b.index); df = (a.loc[ix]-b.loc[ix]).abs()
print('overlap days', len(ix), 'max |diff|', np.nanmax(df.values), 'share >0.001:', np.nanmean(df.values>0.001))

def hourres(w):
    """residual vs centred 2x4 moving average -> array days x 24 x 4 (hours 1..22 valid)"""
    x = w.values.astype(float); k = np.array([.5,1,1,1,.5])/4
    ma = np.full_like(x, np.nan)
    for j in range(2, 94): ma[:, j] = x[:, j-2:j+3] @ k
    return (x-ma).reshape(len(x),24,4)
def boot(daily, B=2000, seed=0):
    """day-block bootstrap SE of mean of daily values"""
    v = daily[np.isfinite(daily)]; rng = np.random.default_rng(seed)
    m = [v[rng.integers(0,len(v),len(v))].mean() for _ in range(B)]
    return v.mean(), np.std(m), np.percentile(m,[2.5,97.5])
W = {'RTM(src1 2023-25)': a, 'RTM(src2 2022-25)': b, 'DAM': wide(d,'DAM','p'), 'GDAM': wide(d,'GDAM','p')}
for name, w in W.items():
    R = hourres(w)
    print(f'\n=== {name}: detrended residual (Rs/kWh) by quarter position; day-bootstrap 95% CI')
    for q in range(4):
        daily = np.nanmean(R[:,:,q], axis=1)
        m, se, ci = boot(daily)
        print(f'  q{q} (:{15*q:02d})  mean {m:+.4f}  se {se:.4f}  CI [{ci[0]:+.4f},{ci[1]:+.4f}]  median-of-days {np.nanmedian(daily):+.4f}  share of days <0: {np.nanmean(daily[np.isfinite(daily)]<0):.3f}')
    # paired test inside the same RTM session: p(:15) - p(:00) and placebo p(:45)-p(:30)
    x = w.values.reshape(len(w),24,4)
    for lab,(i,j) in {'p(:15)-p(:00) [same session]':(0,1), 'p(:45)-p(:30) [same session, placebo]':(2,3)}.items():
        dd = x[:,:,j]-x[:,:,i]; m,se,ci = boot(np.nanmean(dd,axis=1))
        print(f'  {lab}: mean {m:+.4f} CI [{ci[0]:+.4f},{ci[1]:+.4f}]; median block-level {np.nanmedian(dd):+.4f}; share>0 {np.nanmean(dd[np.isfinite(dd)]>0):.3f} share<0 {np.nanmean(dd[np.isfinite(dd)]<0):.3f}')
    print('  q0 residual by hour 1..22:', np.nanmean(R[:,1:23,0],axis=0))
    yrs = w.index.year
    print('  q0 residual by year:', {int(y): round(float(np.nanmean(R[yrs==y][:,:,0])),4) for y in sorted(set(yrs))})
    mo = w.index.month
    print('  q0 residual by month:', np.array([np.nanmean(R[mo==k][:,:,0]) for k in range(1,13)]))
    wd = w.index.dayofweek
    print('  q0 residual Mon..Sun:', np.array([np.nanmean(R[wd==k][:,:,0]) for k in range(7)]))
# quantities (log) for RTM src2
print('\n=== RTM src2 quantities: detrended log residual by quarter (x100 = %)')
for c in ['buy','sell','v']:
    w = d2.pivot(index='date', columns='blk', values=c).reindex(columns=range(96)); w = np.log(w.where(w>0))
    R = hourres(w); print(f'  {c:4s} all-hours q0..q3 %:', 100*np.nanmean(R,axis=(0,1)))
    print(f'       q0 by hour %:', 100*np.nanmean(R[:,1:23,0],axis=0))
