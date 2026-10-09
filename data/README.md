# Input data (not included)

The original analysis reads `data/Diff_runs.txt` using `read.delim()`.

To reproduce the analysis with authorized data, this tab-delimited file should contain:

- 23 rows: one per 400 m distance segment from 0.4 to 9.2 km.
- The first column: segment/distance identifier (excluded by the R script).
- The next 14 columns: heart rate in beats per minute for each run.
- Columns ordered as 6 aerobic, 5 interval, and 3 race workouts, matching `group` in `run_analysis.R`.

The example plotting commands also refer to specific original run column names (`Dec17`, `Oct_25`, and `Jun_Race`). Replace those names if you use a different dataset. All measurements should come from data you are authorized to publish or analyze.
