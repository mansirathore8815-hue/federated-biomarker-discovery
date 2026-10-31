"""
ONE COMMAND for the whole project: experiments, statistics, graphs and a report.

    python run_all.py --preset quick        # first: checks everything works (minutes)
    python run_all.py --preset standard     # the real thing (longer; resumable - just run it again)
    python run_all.py --synthetic --preset quick   # no data needed; wiring test only

Open  results/<preset>/report.html  in a browser when it finishes. Every step is wrapped: if one
fails, the others still run and the report lists the problem. Interrupted? Run the same command
again - finished runs are not repeated.
"""
from __future__ import annotations
import argparse
import json
import os
import time
import traceback
import warnings

import pandas as pd

from config import Config
from run_experiment import run
from src import stats
from src.benchmark import SWEEPS, run_benchmark, run_sweeps
from src.data_loader import COHORTS, DEFAULT_COHORTS, available
from src.report import build_report

SWEEP_BUDGET = dict(rounds=3, local_gen=3, pop_size=20, front_cap=20, top_genes=150, n_clients=5)
PRESETS = {
    "quick": dict(seeds=3, sweep_seeds=1, sweep_budget=SWEEP_BUDGET,
                  budget=dict(rounds=3, local_gen=3, pop_size=20, front_cap=20, top_genes=150, n_clients=5)),
    "standard": dict(seeds=6, sweep_seeds=3, sweep_budget=SWEEP_BUDGET,
                     budget=dict(rounds=5, local_gen=5, pop_size=30, front_cap=30, top_genes=300, n_clients=5)),
    "smoke": dict(seeds=2, sweep_seeds=1,                      # used by the automated tests only
                  sweep_budget=dict(rounds=1, local_gen=1, pop_size=8, front_cap=8, top_genes=40, n_clients=3),
                  budget=dict(rounds=2, local_gen=2, pop_size=10, front_cap=10, top_genes=40, n_clients=3)),
}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", choices=list(PRESETS), default="quick")
    ap.add_argument("--cohorts", default=",".join(DEFAULT_COHORTS))
    ap.add_argument("--synthetic", action="store_true", help="single synthetic cohort (wiring test only)")
    ap.add_argument("--demo", action="store_true", help="simulated KIDNEY/LUNG files in Xena format (no download needed)")
    ap.add_argument("--data-dir", default=None, help="folder holding the *_HiSeqV2 files (default: data/)")
    ap.add_argument("--seeds", type=int, default=None)
    ap.add_argument("--variants", default="CH,CHS")
    ap.add_argument("--partitions", default="iid,zipfian")
    ap.add_argument("--sweeps", default="all", help="comma list of %s, 'all' or 'none'" % list(SWEEPS))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)
    warnings.filterwarnings("ignore")
    t0 = time.time()
    pre = PRESETS[a.preset]
    seeds = a.seeds or pre["seeds"]
    if a.demo:
        from src.demo import make_demo
        a.data_dir = a.data_dir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_demo")
        if not os.path.exists(os.path.join(a.data_dir, "KIRC_HiSeqV2")):
            print("creating demo data in", a.data_dir, flush=True)
            make_demo(a.data_dir)
    if a.data_dir:
        os.environ["FED_DATA_DIR"] = a.data_dir
    out = a.out_dir or os.path.join("results", a.preset + ("_demo" if a.demo else "_synthetic" if a.synthetic else ""))
    os.makedirs(out, exist_ok=True)
    real = not a.synthetic
    problems = []
    ctx = dict(title="Fed-NSGA2-CH/CHS: federated multi-objective biomarker discovery", preset=a.preset,
                             budget=pre["budget"], seeds=seeds, synthetic=a.synthetic or a.demo, problems=problems,
                             main={}, cohort_notes={}, narrative=[], bench_dir=os.path.join(out, "benchmark"))

    def step(name, fn):
        s = time.time()
        print(f"\n=== {name} ===", flush=True)
        try:
            r = fn()
            print(f"--- {name}: done in {time.time() - s:.0f}s", flush=True)
            return r
        except Exception as e:  # noqa: BLE001
            problems.append(f"{name} FAILED: {type(e).__name__}: {e}")
            print(f"!!! {name} failed:\n{traceback.format_exc()}", flush=True)
            return None

    # ---- 0. which cohorts can run
    cohorts = []
    if real:
        for c in [x.strip().upper() for x in a.cohorts.split(",") if x.strip()]:
            if c not in COHORTS:
                problems.append(f"unknown cohort {c}"); continue
            ok, why = available(c)
            print(f"cohort {c}: {'available' if ok else 'NOT available - ' + why}")
            if ok:
                cohorts.append(c)
            else:
                problems.append(f"cohort {c} skipped: {why}")
        if not cohorts:
            print("\nNo cohort available. Run `python data\\check_data.py` for options, or add --demo to use built-in demo data.")
            raise SystemExit(1)
    else:
        cohorts = ["SYN"]
    ctx["cohorts_run"] = cohorts

    # ---- 1. one full experiment per cohort (figures + per-client tables)
    for c in cohorts:
        def one(c=c):
            d = os.path.join(out, "main", c)
            need = ("summary.csv", "per_client.csv", "run_info.json")
            if all(os.path.exists(os.path.join(d, f)) for f in need):      # resume: reuse a finished experiment
                info = json.load(open(os.path.join(d, "run_info.json")))
                ctx["main"][c] = dict(summary=pd.read_csv(os.path.join(d, "summary.csv")),
                                      per_client=pd.read_csv(os.path.join(d, "per_client.csv")),
                                      rep_genes=info["representative_panel"], stability=info["stability"],
                                      fig_dir=os.path.join(d, "figures"))
                ctx["cohort_notes"].setdefault(c, dict(samples=0, notes=info.get("data", {}).get("notes", [])))
                print("already done - reusing", d)
                return
            res = run(Config(real=real, cohort=c if real else "KIDNEY", out_dir=d, partition="zipfian", variant="CH", **pre["budget"]),
                      verbose=False, save=True)
            ctx["main"][c] = dict(summary=res["summary"], per_client=res["per_client"], rep_genes=res["rep_genes"],
                                  stability=res["stability"], fig_dir=os.path.join(d, "figures"))
            ctx["cohort_notes"].setdefault(c, dict(samples=0, notes=res["data"]["notes"]))
            pd.set_option("display.width", 220, "display.float_format", "{:.3f}".format)
            print(res["summary"].to_string(index=False))
        step(f"single experiment: {c}", one)

    # ---- 2. benchmark over seeds
    bench = step("benchmark", lambda: run_benchmark(cohorts, a.variants.split(","), a.partitions.split(","), seeds, real,
                                                    ctx["bench_dir"], pre["budget"]))
    runs = None
    if bench:
        runs = bench["runs"]
        ctx["bench_agg"] = bench["agg"]
        for k, n in bench["notes"].items():
            ctx["cohort_notes"][k] = n
        for s in bench["skipped"] + bench["failed"]:
            problems.append(f"benchmark: {s}")
        try:
            from src.plots import plot_benchmark
            plot_benchmark(bench["agg"], bench["curves"], ctx["bench_dir"])
        except Exception as e:  # noqa: BLE001
            problems.append(f"benchmark figures skipped ({type(e).__name__})")

    # ---- 3. statistics
    if runs is not None:
        def do_stats():
            sd = os.path.join(out, "stats")
            os.makedirs(os.path.join(sd, "figures"), exist_ok=True)
            ctx["descriptive"] = stats.descriptive(runs)
            ctx["paired"] = stats.paired_vs_fed(runs)
            ctx["variants"] = stats.paired_variants(runs)
            ctx["comm"] = stats.comm_savings(runs)
            for n in ("descriptive", "paired", "variants", "comm"):
                ctx[n].to_csv(os.path.join(sd, f"{n}.csv"), index=False)
            ctx["narrative"] = stats.narrative(ctx["paired"], ctx["variants"], ctx["comm"])
            for line in ctx["narrative"]:
                print(" *", line)
            try:
                from src.plots import plot_forest
                if len(ctx["paired"]):
                    plot_forest(ctx["paired"], os.path.join(sd, "figures", "forest.png"))
                    ctx["forest"] = os.path.join(sd, "figures", "forest.png")
            except Exception as e:  # noqa: BLE001
                problems.append(f"forest plot skipped ({type(e).__name__})")
        step("statistics", do_stats)

    # ---- 4. sweeps on the first cohort
    if a.sweeps != "none":
        def do_sweeps():
            which = None if a.sweeps == "all" else a.sweeps.split(",")
            sd = os.path.join(out, "sweeps")
            r = run_sweeps(cohorts[0], real, pre["sweep_seeds"], sd, pre["sweep_budget"], which)
            if r:
                ctx["sweeps_agg"] = r["agg"]
                try:
                    from src.plots import plot_sweeps
                    os.makedirs(os.path.join(sd, "figures"), exist_ok=True)
                    plot_sweeps(r["runs"], os.path.join(sd, "figures", "sweeps.png"))
                    ctx["sweeps_fig"] = os.path.join(sd, "figures", "sweeps.png")
                except Exception as e:  # noqa: BLE001
                    problems.append(f"sweeps figure skipped ({type(e).__name__})")
        step("sweeps", do_sweeps)

    # ---- 5. report
    ctx["elapsed"] = time.time() - t0
    path = step("report", lambda: build_report(out, ctx))
    print("\n" + "=" * 60)
    print(f"FINISHED in {ctx['elapsed'] / 60:.1f} min.  Problems: {len(problems)}")
    for p in problems:
        print("  -", p)
    if path:
        print(f"OPEN THIS IN YOUR BROWSER:  {os.path.abspath(path)}")
        print(f"Tables (CSV) and figures (PNG) are in: {os.path.abspath(out)}")
    return 0 if path else 1


if __name__ == "__main__":
    raise SystemExit(main())
