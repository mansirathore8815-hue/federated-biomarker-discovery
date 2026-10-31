"""Sweeps only (run_all.py does everything).
    python run_sweeps.py --real --cohort KIDNEY --seeds 3
    python run_sweeps.py                          # synthetic wiring check
Resumable."""
import argparse
import os

from src.benchmark import SWEEPS, run_sweeps

ap = argparse.ArgumentParser()
ap.add_argument("--real", action="store_true")
ap.add_argument("--cohort", default="KIDNEY")
ap.add_argument("--seeds", type=int, default=2)
ap.add_argument("--only", default=None, help="comma list of: %s" % ", ".join(SWEEPS))
ap.add_argument("--out-dir", default=os.path.join("results", "sweeps"))
a = ap.parse_args()
if __name__ == "__main__":
    r = run_sweeps(a.cohort, a.real, a.seeds, a.out_dir,
                   dict(rounds=3, local_gen=3, pop_size=20, front_cap=20, top_genes=150, n_clients=5),
                   a.only.split(",") if a.only else None)
    from src.plots import plot_sweeps
    os.makedirs(os.path.join(a.out_dir, "figures"), exist_ok=True)
    plot_sweeps(r["runs"], os.path.join(a.out_dir, "figures", "sweeps.png"))
    print(r["agg"].round(3).to_string(index=False))
