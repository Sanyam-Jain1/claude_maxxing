# Shortest paths on a cube: exact and high-precision statistics

Research note and reproducible checks for the question: for two uniform random points on the
surface of a unit cube, how often does the shortest surface path use more faces than necessary,
and how long is it on average?

`report.html` is the full write-up (page body; open the published artifact or wrap it in an HTML
skeleton to view locally).

## Results

| Quantity | Value | Status |
|---|---|---|
| Adjacent faces: P(shortest path avoids the shared edge) | ln 2 − 1/2 = 0.19314718… | proved |
| Opposite faces: P(shortest path crosses four faces) | 0.155880815524256750509896159877… | computed (no closed form found) |
| Upper bound for the above from the two corner-move conditions | 0.1598041421674… | proved bound, numerical value |
| Mean geodesic distance, adjacent faces | 1.065853075348094653616848853 | computed |
| Mean geodesic distance, opposite faces | 1.5930341029 | computed |
| Mean geodesic distance, whole unit cube | 1.0629753062 | computed |

## Scripts (run from `scripts/`; needs `numpy sympy mpmath`)

| Script | What it does |
|---|---|
| `brute.py` | Brute-force geodesic distance by enumerating all face-sequence unfoldings with a validity check; validates the closed-form candidate distances. |
| `dominance.py` | Exact integer check that every longer face sequence is dominated everywhere by one of the 3 (adjacent) or 12 (opposite) candidates. |
| `chambers.py` | Exact rational polygon clipping for the six-inequality region; prints the polygon types over a grid of source points. |
| `six.py` | Exact check that the six inequalities cut out exactly the winning region (vertex test on a 150×150 rational grid). |
| `sym1.py` | Symbolic vertices, area formulas and region-boundary curves; writes `areas.pkl` (needed by `cut2.py`). |
| `genus.py` | Shows the cubic boundary curves have no singular points (genus 1). |
| `v0b.py` | The two-inequality upper bound to 60 digits. |
| `cut2.py N` | Four-face probability at N-digit working precision (run `sym1.py` first). |
| `cells.py adj|opp K g` | Independent method: per-source nearest-image cells with closed-form distance integrals, composite Gauss rule with K panels. |
| `adjmean.py` | Adjacent-face mean distance to 40 digits. |
| `pslq.py` | Integer-relation searches (negative result). |
| `mc.cpp` | Brute-force Monte Carlo (`g++ -O3 -o mc mc.cpp -lpthread && ./mc 40000000000`). |

## Not verified

Novelty rests on a limited literature search. Sufficiency of the six inequalities is checked exactly
on a rational grid and by Monte Carlo, not proved symbolically. Whether the four-face probability
has a closed form is open.
