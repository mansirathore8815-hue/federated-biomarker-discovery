"""Every parameter named in the proposal lives here, with the section it comes from."""
from __future__ import annotations
from dataclasses import dataclass, asdict


@dataclass
class Config:
    # data (sec 5.1, 7)
    cohort: str = "KIDNEY"          # KIDNEY | LUNG | BRCA(optional)  (see data_loader.COHORTS)
    real: bool = False              # False = synthetic stand-in
    top_genes: int = 300            # candidate genes after federated variance filter
    external_frac: float = 0.2      # globally reserved external validation set (sec 5.1)
    local_test_frac: float = 0.3    # per-client held-out test split
    # federation simulation (sec 5.1)
    n_clients: int = 5              # proposal: 3-10
    partition: str = "zipfian"      # iid | zipfian
    zipf_s: float = 1.5             # proposal: s = 1.0-2.0
    min_client_size: int = 20
    batch_strength: float = 0.0     # simulated platform batch effects (sec 7); 0 = off
    # algorithm (sec 5.2)
    variant: str = "CH"             # NSGA2 | CH | CHS  (Cattelani & Fortino 2022)
    inner_model: str = "nb"         # nb | lr   (inner model inside the GA fitness)
    executor: str = "inprocess"     # inprocess | process (clients as separate OS processes)
    rounds: int = 5
    local_gen: int = 5
    pop_size: int = 30
    front_cap: int = 30
    cv_folds: int = 3
    # privacy (sec 5.4)
    dp_epsilon: float | None = None  # per-release Laplace epsilon on aggregated fitness; None = off
    # Fed-CHV (sec 5.3): see src/chv_metric.py
    min_test_samples: int = 10      # clients with fewer test samples are skipped in Fed-CHV
    seed: int = 0
    out_dir: str = "results"

    def validate(self):
        assert 3 <= self.n_clients <= 10, "proposal specifies 3-10 clients"
        assert self.partition in ("iid", "zipfian")
        assert 1.0 <= self.zipf_s <= 2.0, "proposal specifies s=1.0-2.0"
        assert self.variant in ("NSGA2", "CH", "CHS")
        assert self.inner_model in ("nb", "lr")
        assert self.executor in ("inprocess", "process")
        assert self.front_cap <= self.pop_size
        return self

    def to_dict(self):
        return asdict(self)
