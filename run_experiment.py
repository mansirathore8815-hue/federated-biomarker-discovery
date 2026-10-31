"""
    python run_experiment.py                                  # synthetic dry run
    python run_experiment.py --real --cohort BRCA --variant CHS --n-clients 8 --dp-epsilon 1
    python run_experiment.py --executor process               # clients as separate OS processes
Research proof-of-concept: simulated institutions; see README for what is and is not verified.
"""
from __future__ import annotations
import argparse
import warnings
from dataclasses import fields

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from config import Config
from src import baselines, comm_cost
from src.chv_metric import chv_from_report, fed_chv, fed_chv_soft
from src.client import make_clients, wire_bytes
from src.data_loader import load_dataset, split_external
from src.evaluate import dice, jaccard, pairwise, panel_accuracies, representative
from src.federation import (aggregate_accuracy, federated_prefilter, run_federated,
                            select_front)
from src.ga_problem import Evaluator, holdout_accuracy, to_min
from src.hypervolume import cross_hypervolume, leanness
from src.partition import apply_batch_effect, partition_indices
from src.reporting import save_all


def _split(Xc, yc, frac, seed):
    try:
        return train_test_split(Xc, yc, test_size=frac, random_state=seed, stratify=yc)
    except ValueError:
        return train_test_split(Xc, yc, test_size=frac, random_state=seed)


def _zscore(X):
    sd = X.std(axis=0)
    return (X - X.mean(axis=0)) / np.where(sd > 0, sd, 1.0)


def run(cfg: Config, verbose=True, save=True, data=None):
    cfg.validate()
    warnings.filterwarnings("ignore")
    log = print if verbose else (lambda *a, **k: None)
    rng = np.random.default_rng(cfg.seed)

    # ---- data: external hold-out, partition, batch effects, federated filter + standardisation
    ds = data if data is not None else load_dataset(cfg)
    X, y, genes = ds["X"], ds["y"], ds["genes"]
    n_raw = X.shape[1]
    if ds["ext"] is None:
        Xa, ya, Xe, ye = split_external(X, y, cfg.external_frac, cfg.seed)
        ext_kind, ext_independent = "held-out TCGA split (NOT a different platform)", False
    else:
        Xa, ya = X, y
        Xe, ye = ds["ext"]
        ext_kind, ext_independent = "independent external cohort", True
    parts = partition_indices(len(ya), cfg.n_clients, cfg.partition, cfg.zipf_s,
                              cfg.min_client_size, cfg.seed)
    sd, gm = Xa.std(axis=0), Xa.mean(axis=0)
    raw = [(apply_batch_effect(Xa[i], cfg.batch_strength, rng, sd, gm), ya[i]) for i in parts]
    keep, mu, sdv = federated_prefilter([x for x, _ in raw], cfg.top_genes)
    gene_names = [genes[i] for i in keep]
    G = len(keep)
    Xe = _zscore(Xe[:, keep]) if ext_independent else (Xe[:, keep] - mu) / sdv
    ext = (Xe.astype(np.float32), ye)
    train, test = [], []
    for k, (x, yy) in enumerate(raw):
        xz = ((x[:, keep] - mu) / sdv).astype(np.float32)
        xtr, xte, ytr, yte = _split(xz, yy, cfg.local_test_frac, cfg.seed + k)
        train.append((xtr, ytr)); test.append((xte, yte))
    pooled_tr = (np.vstack([a for a, _ in train]), np.concatenate([b for _, b in train]))
    test_sizes = [len(b) for _, b in test]
    log(f"[data] {ds['name']}: {len(y)} samples x {n_raw} genes | external: {len(ye)} ({ext_kind})")
    for n in ds["notes"]:
        log(f"       {n}")
    log(f"[clients] {cfg.n_clients} x {cfg.partition} train {[len(b) for _, b in train]} "
        f"test {test_sizes} | {G} genes after federated filter | executor={cfg.executor}")

    setup = comm_cost.prefilter_bytes(cfg.n_clients, n_raw, G)
    cen_bytes = comm_cost.centralized_bytes(len(pooled_tr[1]), n_raw)
    clients = make_clients(train, test, cfg.to_dict())
    try:
        log(f"[run] Fed-NSGA2-{cfg.variant} ({cfg.inner_model}): {cfg.rounds} rounds x {cfg.local_gen} "
            f"gens, pop {cfg.pop_size}, DP eps={cfg.dp_epsilon}")
        fed = run_federated(clients, cfg, G, pooled_tr, ext, setup_extra_bytes=setup)

        b = wire_bytes(clients)
        nav_masks = baselines.naive_fedavg(clients, G)
        nav_bytes = wire_bytes(clients) - b + setup
        b = wire_bytes(clients)
        loc_masks = baselines.local_only(clients, cfg, G)
        aggregate_accuracy(clients, loc_masks)           # the one scoring round local-only needs
        loc_bytes = wire_bytes(clients) - b + setup
        cen_masks = baselines.centralized(pooled_tr, cfg)
        methods = {
            f"Fed-NSGA2-{cfg.variant}": (fed["masks"], fed["comm_bytes"]),
            f"Centralized NSGA2-{cfg.variant}": (cen_masks, cen_bytes),
            "Local-only (no aggregation)": (loc_masks, loc_bytes),
            "Naive FedAvg (ANOVA)": (nav_masks, nav_bytes),
        }

        # ---- one common yardstick for every method (noise-free federated scoring)
        pooled_ev = Evaluator(*pooled_tr, cfg.cv_folds, cfg.inner_model, cfg.seed)
        rows, panels, keepers = [], {}, {}
        for name, (m, wire) in methods.items():
            agg = aggregate_accuracy(clients, m)
            m, F = select_front(m, to_min(agg, m.sum(axis=1)), cfg.front_cap)
            acc = 1 - F[:, 0]
            from src.comm_cost import pack_masks
            reports = [c.call("holdout_report", pack_masks(m)) for c in clients]
            fc = fed_chv(m, reports, test_sizes, cfg.min_test_samples)
            # plain CHV with one pooled "client" (the original, centralized definition)
            sizes = m.sum(axis=1)
            pooled_te = np.array([holdout_accuracy(mm, *pooled_tr, *[np.vstack([a for a, _ in test]),
                                  np.concatenate([b for _, b in test])], cfg.inner_model) for mm in m])
            chv_p = cross_hypervolume(np.column_stack([pooled_ev(m), leanness(sizes)]),
                                      np.column_stack([pooled_te, leanness(sizes)]))
            sc = panel_accuracies(m, *pooled_tr, *ext, cfg.inner_model)
            le10 = sizes <= 10
            rows.append(dict(method=name, front_size=len(m), fed_chv=fc["fed_chv"],
                             fed_chv_worst_client=fc["fed_chv_min"],
                             fed_chv_soft=fed_chv_soft(m, reports, test_sizes, cfg.min_test_samples),
                             chv_pooled=chv_p, best_ext_bal_acc=sc.max(), size_at_best=int(sizes[sc.argmax()]),
                             best_ext_bal_acc_le10=sc[le10].max() if le10.any() else np.nan,
                             min_panel=int(sizes.min()), wire_bytes=int(wire)))
            panels[name] = pd.DataFrame(dict(
                size=sizes, agg_cv_acc=acc, ext_bal_acc=sc,
                genes=[";".join(np.array(gene_names)[mm]) for mm in m]))
            keepers[name] = (m, acc, reports, fc)
    finally:
        for c in clients:
            c.close()
    summary = pd.DataFrame(rows)

    # ---- per-client + aggregate report for the federated method
    fm, facc, freports, ffc = keepers[f"Fed-NSGA2-{cfg.variant}"]
    rep = representative(fm, facc)
    pc = []
    for c in range(cfg.n_clients):
        ok = test_sizes[c] >= cfg.min_test_samples
        pc.append(dict(client=c, n_train=len(train[c][1]), n_test=test_sizes[c],
                       panel_size=int(fm[rep].sum()), test_bal_acc=float(freports[c][1][rep]),
                       client_chv=ffc["per_client"][c] if ok else np.nan,
                       local_best_size=int(fed["client_best"][c].sum()),
                       jaccard_global_vs_local=jaccard(fm[rep], fed["client_best"][c])))
    per_client = pd.DataFrame(pc)
    stab = dict(jaccard=pairwise(fed["client_best"], jaccard), dice=pairwise(fed["client_best"], dice))
    res = dict(config=cfg.to_dict(), summary=summary, per_client=per_client,
               history=pd.DataFrame(fed["history"]), panels=panels, stability=stab,
               dp=dict(epsilon_per_release=cfg.dp_epsilon, releases=fed["n_releases"],
                       composed_epsilon=fed["eps_total"]),
               comm=dict(federated_total=fed["comm_bytes"], setup=setup, centralized_bytes=cen_bytes),
               rep_genes=[gene_names[i] for i in np.where(fm[rep])[0]],
               data=dict(name=ds["name"], external=ext_kind, notes=ds["notes"]),
               reports=freports, client_best=fed["client_best"])
    if save:
        save_all(res, cfg.out_dir)
    if verbose:
        pd.set_option("display.width", 220, "display.float_format", "{:.3f}".format)
        print("\n== Method comparison (same yardstick; wire_bytes = search communication) ==\n",
              summary.to_string(index=False))
        print("\n== Per client (representative panel) ==\n", per_client.to_string(index=False))
        print("\n== Aggregate ==\n", per_client.drop(columns="client").mean().to_string())
        print(f"\nclient-panel stability: Jaccard {stab['jaccard']:.3f}  Dice {stab['dice']:.3f}")
        print(f"communication: federated {res['comm']['federated_total']:,} B vs centralized "
              f"{cen_bytes:,} B (measured message sizes)")
        if cfg.dp_epsilon:
            print(f"DP: eps/release {cfg.dp_epsilon}, {fed['n_releases']} releases -> composed eps "
                  f"{fed['eps_total']:.0f} (heuristic sensitivity, not a certified guarantee)")
        print("representative panel:", res["rep_genes"])
    return res


def parse_args() -> Config:
    p = argparse.ArgumentParser()
    for f in fields(Config):
        flag, d = "--" + f.name.replace("_", "-"), f.default
        if isinstance(d, bool):
            p.add_argument(flag, action="store_true")
        elif d is None:
            p.add_argument(flag, type=float, default=None)
        else:
            p.add_argument(flag, type=type(d), default=d)
    return Config(**vars(p.parse_args()))


if __name__ == "__main__":
    run(parse_args())
