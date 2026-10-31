# Fed-NSGA2-CH / CHS - federated multi-objective biomarker discovery

Datasets: **KIDNEY** (TCGA KIRC+KIRP+KICH -> ccRCC / PRCC / ChRCC) and **LUNG** (TCGA LUAD+LUSC -> LUAD / LUSC).
Both need ONLY the expression files you already downloaded - no clinical files. Labels are the TCGA cohort
of origin. (BRCA/PAM50 is an optional extra that switches on only if you obtain `BRCA_clinicalMatrix`.)

## Step 1 - get the data (choose ONE)
The `data\` folder of a fresh unzip is EMPTY (the big files are not in the zip). KIDNEY and LUNG need these five
files, unzipped, named exactly: `KIRC_HiSeqV2  KIRP_HiSeqV2  KICH_HiSeqV2  LUAD_HiSeqV2  LUSC_HiSeqV2`.

**A. Copy them from your earlier folder** (you already downloaded them there):
```
Copy-Item C:\Users\hp\Desktop\fed-nsga2-ch-v5\fed-nsga2-ch\data\*_HiSeqV2 data\
```
**B. Download** (internet, ~250 MB): `python data\download_data.py`
Direct links (same ones the script uses; the .gz must be unzipped, e.g. with 7-Zip, and renamed as above):
```
https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz   ->  data\KIRC_HiSeqV2
https://tcga.xenahubs.net/download/TCGA.KIRP.sampleMap/HiSeqV2.gz   ->  data\KIRP_HiSeqV2
https://tcga.xenahubs.net/download/TCGA.KICH.sampleMap/HiSeqV2.gz   ->  data\KICH_HiSeqV2
https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz   ->  data\LUAD_HiSeqV2
https://tcga.xenahubs.net/download/TCGA.LUSC.sampleMap/HiSeqV2.gz   ->  data\LUSC_HiSeqV2
```
**C. No data at all - built-in DEMO cohort:** add `--demo` to any `run_all.py` command. It generates simulated files in the
same format and folder layout, with the real cohorts' sample counts and class imbalance (kidney 530/290/66, lung 515/501),
correlated genes and overlapping classes, and runs the same real-data code path. It is SIMULATED: use it to see that
everything works and what the outputs look like; the report carries a banner and the numbers say nothing about cancer.

## Step 2 - check, then run (PowerShell, from this folder)
```
python data\check_data.py                 # must say READY (add --demo to check the demo data)
python run_all.py --preset quick           # verifies the whole chain, roughly 10-20 min
python run_all.py --preset standard        # the full study, roughly 1 h; resumable
python run_all.py --preset quick --demo    # same thing on the demo data
```
(`pip` itself may be blocked on your PC; use `python -m pip install ...` if a package is ever missing.)
Double-clicking `run_quick.bat` / `run_standard.bat` does the same. If a run is interrupted, run the same command again:
finished runs and finished per-cohort experiments are skipped. A failing step never stops the others; it is listed in the report.

## What you get (`results\<preset>\`)
* `report.html` - open in a browser: findings, tables, statistical tests and all graphs in one page
* `REPORT.md`, and CSV tables: `benchmark\summary.csv|runs.csv`, `stats\paired.csv|descriptive.csv|variants.csv|comm.csv`,
  `sweeps\sweeps_summary.csv`, `main\<COHORT>\summary.csv|per_client.csv|panels_*.csv`
* Graphs (PNG): per cohort `main\<COHORT>\figures\` (Pareto fronts, method comparison, accuracy vs communication,
  per-client, train-vs-test, stability); `benchmark\figures\`; `stats\figures\forest.png`; `sweeps\figures\sweeps.png`

## What is measured
Objectives (maximised): inner-CV balanced accuracy and leanness 1/(n+1). Methods: Fed-NSGA2 (CH or CHS),
Centralized NSGA2, Local-only (no aggregation), Naive FedAvg (ANOVA). Metrics: Fed-CHV (also soft-leanness
version and the original pooled CHV), external balanced accuracy (best, and best with <= 10 genes), measured
communication bytes, per-client accuracy/CHV, Jaccard/Dice panel stability.
Statistics: mean, sd and 95% t-interval per method; paired Wilcoxon signed-rank across seeds with Holm correction,
bootstrap 95% CI of the difference, Cohen's dz, win rate; CHS vs CH; communication ratio.
Sweeps: client skew (IID / Zipf s=1,1.5,2), platform batch effects, DP noise, number of clients, algorithm variant.

## Reading the results - one important caveat
The search budget here (standard: 5 rounds x 5 generations, population 30) is orders of magnitude smaller than the
benchmark paper's (population 500, 1000 generations). With a small budget the genetic search finds only small panels and a
plain ANOVA ranking can match or beat it; that is a budget effect, not a bug. To give the GA a fair chance raise
`--rounds`/`--local-gen`/`--pop-size` in `run_experiment.py` or edit the `standard` preset in `run_all.py`.

## Verified vs not
* Verified by automated tests (38): CHV equals the authors' reference code to 1e-16 on random inputs; GA operators
  (clone index, LASSO init, clone repurposing, symmetric mutation); Fed-CHV properties; separate-OS-process clients give
  identical results to in-process; loaders on Xena-format files; the whole `run_all` chain including resume.
* Your own KIDNEY benchmark ran 20/20 runs without failure on the real files.
* NOT verified: real-data scientific conclusions - read the report's "Limitations" section. In particular: labels are
  coarse (cohort of origin), the "external" set is a TCGA hold-out (same platform), Fed-CHV's federated extension is this
  project's own definition, DP uses a heuristic sensitivity, and the `quick` preset (3 seeds) cannot reach p < 0.05
  (smallest possible Wilcoxon p is 0.25); `standard` uses 6 seeds (smallest 0.031).

## Other commands
```
python run_experiment.py --real --cohort KIDNEY --variant CHS --executor process   # one run, clients as OS processes
python run_benchmark.py --real --seeds 6      python run_sweeps.py --real --cohort KIDNEY --seeds 3
python run_all.py --synthetic --preset quick   # no data needed (wiring test only)
python -m pip install pytest ; python -m pytest -q
```
