"""Attempts to break / explain the top-of-hour (q0) RTM discount."""
import sys; sys.path.insert(0, 'scripts')
from load import *
np.set_printoptions(linewidth=250, precision=3, suppress=True)
d = load()
def hourres(w):
    x = w.values.astype(float); k = np.array([.5,1,1,1,.5])/4
    ma = np.full_like(x, np.nan)
    for j in range(2, 94): ma[:, j] = x[:, j-2:j+3] @ k
    return (x-ma).reshape(len(x),24,4)
days = wide(d,'RTM','p').dropna(thresh=96).index.intersection(wide(d,'DAM','p').index)
W = {(m,c): wide(d,m,c).loc[days] for m in ['DAM','RTM'] for c in ['p','buy','sell','v']}
Rr = hourres(W['RTM','p']); Rd = hourres(W['DAM','p'])
def bootci(daily, B=2000):
    v = daily[np.isfinite(daily)]; rng = np.random.default_rng(1)
    m = [v[rng.integers(0,len(v),len(v))].mean() for _ in range(B)]; return v.mean(), np.percentile(m,[2.5,97.5])
print('1) Is it tails / the price cap?')
x = W['RTM','p'].values.reshape(len(days),24,4)
ok = (x.max(axis=2) < 9.99) & (x.min(axis=2) > 0.5)
print('   q0 resid, hours with all 4 blocks strictly inside (0.5, 9.99):', np.nanmean(np.where(ok, Rr[:,:,0], np.nan)), ' share of hours kept', ok.mean())
r0 = Rr[:,:,0]; lo, hi = np.nanpercentile(r0,[1,99]); print('   q0 resid winsorised 1/99%:', np.nanmean(np.clip(r0,lo,hi)), ' median:', np.nanmedian(r0), ' share<0:', np.nanmean(r0[np.isfinite(r0)]<0), 'share>0:', np.nanmean(r0[np.isfinite(r0)]>0))
rel = Rr[:,:,0] / np.nanmean(x,axis=2); print('   q0 resid as % of hour-mean price (mean of ratios):', 100*np.nanmean(rel))
print('2) DAM-RTM spread by quarter (detrended), day-bootstrap CI')
Rs = Rd - Rr
for q in range(4):
    m, ci = bootci(np.nanmean(Rs[:,:,q],axis=1)); print(f'   q{q}: {m:+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}]')
print('3) Mechanism regression on (day,hour) cells, hours 1..22: RTM q0 residual ~ a + b1*price trend + b2*DAM volume step at boundary + b3*DAM q0 price resid')
hm = np.nanmean(x,axis=2); tr = np.full_like(hm,np.nan); tr[:,1:-1] = (hm[:,2:]-hm[:,:-2])/8
vd = W['DAM','v'].values; vstep = np.full((len(days),24),np.nan)
vstep[:,1:] = vd[:,4::4] - vd[:,3:-1:4]                  # DAM cleared volume: block h:00 minus block (h-1):45, GW
vtr = np.full((len(days),24),np.nan); vh = vd.reshape(len(days),24,4).mean(axis=2); vtr[:,1:-1] = (vh[:,2:]-vh[:,:-2])/8
exc = vstep - vtr                                          # excess DAM volume step beyond smooth trend (GW)
bb = np.log(W['RTM','buy'].values).reshape(len(days),24,4).mean(axis=2); btr = np.full_like(bb,np.nan); btr[:,1:-1]=(bb[:,2:]-bb[:,:-2])/8
Y = Rr[:,:,0].ravel(); X = np.c_[np.ones(Y.size), tr.ravel(), exc.ravel(), Rd[:,:,0].ravel(), btr.ravel()]
okk = np.isfinite(Y) & np.isfinite(X).all(axis=1)
names = ['intercept','RTM price trend/blk','DAM excess vol step (GW)','DAM q0 price resid','RTM log-buy trend/blk']
for cols in ([0],[0,1],[0,2],[0,3],[0,4],[0,1,2,3,4]):
    A = X[okk][:,cols]; b,*_ = np.linalg.lstsq(A, Y[okk], rcond=None); e = Y[okk]-A@b
    # day-clustered SE
    g = np.repeat(np.arange(len(days)),24)[okk]; bread = np.linalg.inv(A.T@A); S = np.zeros((len(cols),len(cols)))
    sc = A*e[:,None]; dfc = pd.DataFrame(sc).groupby(g).sum().values; S = dfc.T@dfc; se = np.sqrt(np.diag(bread@S@bread))
    print('   ', '; '.join(f'{names[c]}={bb_:+.4f} (se {s:.4f})' for c,bb_,s in zip(cols,b,se)), f'  R2={1-e.var()/Y[okk].var():.4f}')
print('   mean DAM excess volume step at hour boundary (GW):', np.nanmean(exc), ' by hour 1..22:', np.nanmean(exc[:,1:23],axis=0))
print('4) Same-session relative value, split-sample by year: mean p(:15)-p(:00), and p(mean of :15,:30,:45)-p(:00), Rs/kWh')
yr = days.year
for y in sorted(set(yr)):
    xx = x[yr==y]; print(f'   {y}: p15-p00 {np.nanmean(xx[:,:,1]-xx[:,:,0]):+.4f}   rest-of-hour minus p00 {np.nanmean(xx[:,:,1:].mean(axis=2)-xx[:,:,0]):+.4f}   hours with p00<p15: {np.nanmean(xx[:,:,0]<xx[:,:,1]):.3f}, >: {np.nanmean(xx[:,:,0]>xx[:,:,1]):.3f}')
print('5) Hour-of-day stability: correlation of the 22-hour q0 profile between years')
prof = np.array([np.nanmean(Rr[yr==y][:,1:23,0],axis=0) for y in sorted(set(yr))]); print(np.corrcoef(prof)); print(prof)
print('6) Full detrended RTM shape per hour (rows hour 1..22, cols q0..q3)'); print(np.nanmean(Rr[:,1:23,:],axis=0))
print('7) RTM volume-weighted: avg price paid if buying equal energy each block vs only q0, Rs/kWh:', np.nanmean(x), np.nanmean(x[:,:,0]))
