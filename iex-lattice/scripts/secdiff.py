"""Curvature-robust check: second difference SD_j = p_j - (p_{j-1}+p_{j+1})/2 centred on each block.
A pure dip of depth delta at :00 gives SD(:00)=-delta, SD(:15)=SD(:45)=+delta/2, SD(:30)=0; smooth curvature hits all equally."""
import sys, json; sys.path.insert(0, 'scripts')
from load import *
np.set_printoptions(linewidth=250, precision=3, suppress=True)
d = load()
out = {}
for m in ['RTM','DAM','GDAM']:
    w = wide(d,m,'p'); x = w.values
    sd = np.full_like(x, np.nan); sd[:,1:-1] = x[:,1:-1] - (x[:,:-2]+x[:,2:])/2
    pos = np.arange(96) % 4
    print(m, 'mean SD by position :00 :15 :30 :45 ->', [round(float(np.nanmean(sd[:,pos==q])),4) for q in range(4)])
    S = sd.reshape(len(x),24,4)
    dip = S[:,1:,0] - (S[:,1:,2] + S[:,:-1,2])/2          # SD(h:00) minus avg of SD(h:30), SD(h-1:30): hours 1..23
    daily = np.nanmean(dip,axis=1); v = daily[np.isfinite(daily)]; rng = np.random.default_rng(0)
    bs = [v[rng.integers(0,len(v),len(v))].mean() for _ in range(2000)]
    print(f'   curvature-corrected :00 effect: {v.mean():+.4f}  CI [{np.percentile(bs,2.5):+.4f},{np.percentile(bs,97.5):+.4f}]  days<0: {np.mean(v<0):.3f}')
    byh = np.nanmean(dip,axis=0); se = np.nanstd(dip,axis=0)/np.sqrt(np.sum(np.isfinite(dip),axis=0))
    print('   by hour 1..23:', byh); print('   t-stat      :', byh/se)
    yr = w.index.year
    byy = {int(y): round(float(np.nanmean(dip[yr==y])),4) for y in sorted(set(yr))}; print('   by year', byy)
    out[m] = dict(pos=[float(np.nanmean(sd[:,pos==q])) for q in range(4)], mean=float(v.mean()), ci=[float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))],
                  share_neg=float(np.mean(v<0)), byhour=[float(a) for a in byh], se=[float(a) for a in se], byyear=byy)
    if m=='RTM':
        mp = np.nanmean(x,axis=0); out['rtm_profile']=[float(a) for a in mp]
    if m=='DAM': out['dam_profile']=[float(a) for a in np.nanmean(x,axis=0)]
json.dump(out, open('results.json','w'))
