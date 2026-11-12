"""Load the public IEX 15-min dataset (divyeshBhartiya/IEXInsiderMCP, IEX_Market_Data.csv)."""
import sys, pandas as pd, numpy as np
SRC = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1].endswith('.csv') else '/home/user/divyeshbhartiya/iexinsidermcp/IEX_Market_Data.csv'
def load(src=SRC):
    d = pd.read_csv(src, low_memory=False).iloc[:, :8]
    d.columns = ['mkt', 'year', 'date', 'tb', 'buy', 'sell', 'p', 'v']
    for c in ['buy', 'sell', 'p', 'v']:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d['date'] = pd.to_datetime(d['date'])
    t = pd.to_timedelta(d['tb'])
    d['blk'] = (t.dt.total_seconds() // 900).astype(int)      # 0..95, block start
    d = d.drop_duplicates(['mkt', 'date', 'blk']).sort_values(['mkt', 'date', 'blk'])
    d['ts'] = d['date'] + t
    return d.reset_index(drop=True)
def wide(d, mkt, col):
    w = d[d.mkt == mkt].pivot(index='date', columns='blk', values=col)
    return w.reindex(columns=range(96))
