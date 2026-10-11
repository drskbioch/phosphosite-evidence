# Number checkers as frozen at tag rev6-final

This folder holds the three number checkers (`check_numbers_sourced.py`, `check_bindings.py`, `check_literals.py`) exactly as they were at the tag `rev6-final` (commit `3521d5c`). They are kept for reference only. They are not part of the pipeline and `REPRODUCE.md` does not run them. The checkers to run are the ones in the folder above.

`check_literals.py` in this folder reports **three unclassified literals on the current claims builder** (`build_claims.py`: `2018` at lines 408 and 412, `2017` at line 411, the restriction years of the predictor sensitivity rows). It does so by design: it is the frozen version, and the current claims builder contains rows that were added after the freeze. The current `check_literals.py` classifies those two years as definitions (the plan names them in section 11, set at the freeze) and reports 0 unclassified.

What changed between these copies and the current checkers, and why, is recorded in section 15 of `study/analysis/ANALYSIS_PLAN.md` (the entries headed "Post-freeze extension", `check_bindings.py`: adds a check; `check_numbers_sourced.py`: removes sources, which makes the check stricter; `check_literals.py`: adds an exemption).
