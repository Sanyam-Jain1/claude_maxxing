# IEX top-of-hour discount

Research note: on IEX's Real-Time Market the 15-minute block starting at :00 clears about ₹0.10/kWh
below the local trend (2022–2025), including against the :15 block sold in the same auction.
See `report.html` for the question, evidence, failed hypotheses, prior art and tests.

Data (public, not redistributed here):
- `divyeshBhartiya/IEXInsiderMCP` → `IEX_Market_Data.csv` (DAM/GDAM/RTM, Jan 2023 – Sep 2025)
- `nandanrajesh/Indian-Real-Time-Market-Analysis` → RTM snapshot xlsx (Apr 2022 – Jun 2025)

Scripts (`scripts/`, need pandas, numpy, openpyxl; paths to the data are at the top of `load.py` / `toh.py`):
- `explore.py`  step sizes by clock position
- `profile.py`  mean signed steps per block
- `sawtooth.py` ramp-conditional sawtooth test (rejected)
- `toh.py`      two-source robustness, bootstrap CIs, by hour/year/month/weekday, bid volumes
- `mech.py`     outlier, price-cap, regression and split-sample checks
- `secdiff.py`  curvature-robust second-difference statistic; writes `results.json`
- `build_report.py` renders `report.html` from `template.html`
