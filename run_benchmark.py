"""Benchmark only (run_all.py does everything).
    python run_benchmark.py --real --cohorts KIDNEY,LUNG --seeds 6
    python run_benchmark.py                       # synthetic wiring check
Resumable: run the same command again after an interruption."""
import argparse
import os

from src.benchmark import run_benchmark
from src.data_loader import DEFAULT_COHORTS

ap = argparse.ArgumentParser()
ap.add_argument("--real", action="store_true")
ap.add_argument("--cohorts", default=",".join(DEFAULT_COHORTS))
ap.add_argument("--variants", default="CH,CHS")
ap.add_argument("--partitions", default="iid,zipfian")
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--n-clients", type=int, default=5)
ap.add_argument("--rounds", type=int, default=5)
ap.add_argument("--local-gen", type=int, default=5)
ap.add_argument("--pop-size", type=int, default=30)
ap.add_argument("--top-genes", type=int, default=300)
ap.add_argument("--out-dir", default=os.path.join("results", "benchmark"))
a = ap.parse_args()
if __name__ == "__main__":
    r = run_benchmark(a.cohorts.split(","), a.variants.split(","), a.partitions.split(","), a.seeds, a.real, a.out_dir,
                      dict(rounds=a.rounds, local_gen=a.local_gen, pop_size=a.pop_size, front_cap=a.pop_size,
                           top_genes=a.top_genes, n_clients=a.n_clients))
    if r is None:
        raise SystemExit("No run completed; see skipped_or_failed.csv")
    print(r["agg"].round(3).to_string(index=False))
    from src.plots import plot_benchmark
    plot_benchmark(r["agg"], r["curves"], a.out_dir)
