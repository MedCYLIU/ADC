"""Exploratory, patient-level reanalysis of GSE243375; no selected cutpoint."""

from pathlib import Path
import pandas as pd, numpy as np, json, math

r = Path(__file__).parent
p = pd.read_csv(r / "jci_patient_baseline.csv")
s = pd.read_csv(r / "jci_sample_joined.csv")
mu = p.ERBB2.mean()
sd = p.ERBB2.std(ddof=1)


def logistic(X, y, offset=0):
    b = np.zeros(X.shape[1])
    ll = -np.inf
    for _ in range(100):
        eta = offset + X @ b
        prob = np.exp(-np.logaddexp(0, -eta))
        score = X.T @ (y - prob)
        info = X.T @ (X * (prob * (1 - prob))[:, None])
        if np.max(abs(score)) < 1e-9:
            break
        step = np.linalg.solve(info, score)
        curr = np.sum(y * eta - np.logaddexp(0, eta))
        fac = 1
        while fac > 2**-20:
            bn = b + fac * step
            en = offset + X @ bn
            new = np.sum(y * en - np.logaddexp(0, en))
            if new >= curr - 1e-12:
                break
            fac /= 2
        b = bn
    assert np.max(abs(score)) < 1e-6, score
    eta = offset + X @ b
    return b, float(np.sum(y * eta - np.logaddexp(0, eta)))


def model(d, col, label):
    y = d.pcr.to_numpy()
    z = (d[col].to_numpy() - mu) / sd
    base = np.column_stack([np.ones(len(d)), d.hr])
    X = np.column_stack([base, z])
    b, ll = logistic(X, y)

    def lr(v):
        return 2 * (ll - logistic(base, y, v * z)[1])

    ci = []
    for dr in [-1, 1]:
        a = b[-1]
        c = a + dr * 0.5
        while lr(c) < 3.841458820694124:
            c = a + 2 * (c - a)
        for _ in range(55):
            mid = (a + c) / 2
            if lr(mid) < 3.841458820694124:
                a = mid
            else:
                c = mid
        ci.append(float(np.exp((a + c) / 2)))
    lr0 = 2 * (ll - logistic(base, y)[1])
    return {
        "model": label,
        "n": len(d),
        "pcr": int(y.sum()),
        "beta": float(b[-1]),
        "OR": float(np.exp(b[-1])),
        "lo": ci[0],
        "hi": ci[1],
        "p_LR": math.erfc(math.sqrt(max(0, lr0) / 2)),
        "log_likelihood": ll,
    }


q = p[p.ncores == 2].copy()
models = [
    model(p, "ERBB2", "All baseline patients: mean"),
    model(q, "ERBB2", "Two-core subset: mean"),
    model(q, "ERBB2_min", "Two-core subset: lower core"),
    model(q, "ERBB2_max", "Two-core subset: higher core"),
]
pd.DataFrame(models).to_csv(r / "jci_models.csv", index=False)
rng = np.random.default_rng(243375)
B = 5000
y = p.pcr.to_numpy()
base = np.column_stack([np.ones(len(p)), p.hr])
single = []
for b in range(B):
    x = np.where(rng.integers(2, size=len(p)), p.ERBB2_min, p.ERBB2_max)
    z = (x - mu) / sd
    bet, ll = logistic(np.column_stack([base, z]), y)
    single.append(float(np.exp(bet[-1])))
pd.DataFrame({"iteration": np.arange(B) + 1, "OR": single}).to_csv(
    r / "single_core_resampling.csv", index=False
)
# Patient bootstrap of descriptive sampling variation, all repeated-core subjects have equal weight.
a = q["range"].to_numpy()
boot = np.median(a[rng.integers(len(q), size=(B, len(q)))], axis=1)
base_s = s[s.Treatment_Group == "Pre"].copy()
post = s[s.Treatment_Group == "Post"].copy()
pair = (
    post[["Patient ID", "ERBB2"]]
    .rename(columns={"ERBB2": "post_ERBB2"})
    .merge(p[["Patient ID", "ERBB2", "pcr", "hr", "het"]], on="Patient ID")
)
assert (pair.pcr == 0).all()
pair["delta"] = pair.post_ERBB2 - pair.ERBB2
pair.to_csv(r / "jci_paired_change.csv", index=False)
d = pair.delta.to_numpy()
boots = np.median(d[rng.integers(len(d), size=(B, len(d)))], axis=1)
# Descriptive association of core-to-core variation with deposited FISH variation.
base_s["ERBB2_FISH_Amp%"] = pd.to_numeric(base_s["ERBB2_FISH_Amp%"], errors="coerce")
fish = base_s.groupby("Patient ID")["ERBB2_FISH_Amp%"].agg(lambda v: v.max() - v.min())
q["FISH_range"] = q["Patient ID"].map(fish)
rho = float(np.corrcoef(q["range"].rank(), q.FISH_range.rank())[0, 1])
q.to_csv(r / "jci_two_core_patient.csv", index=False)
summary = {
    "cohort": "NCT02326974 / GSE243375",
    "unit": "patient",
    "baseline_n": len(p),
    "pcr": int(p.pcr.sum()),
    "hr_positive": int(p.hr.sum()),
    "heterogeneous": int(p.het.sum()),
    "mean_ERBB2": mu,
    "SD_ERBB2": sd,
    "models": models,
    "same_scale": "All ORs per SD of the 129 patient baseline mean ERBB2; scale reused for all subset/single-core models.",
    "two_core_n": len(q),
    "within_patient_abs_difference_median": float(np.median(a)),
    "within_patient_abs_difference_IQR": np.quantile(a, [0.25, 0.75]).tolist(),
    "within_patient_abs_difference_median_CI": np.quantile(
        boot, [0.025, 0.975]
    ).tolist(),
    "within_patient_abs_difference_gt1_n": int((a > 1).sum()),
    "within_patient_abs_difference_max": float(a.max()),
    "FISH_range_Spearman_rho": rho,
    "single_core_selection_OR_median": float(np.median(single)),
    "single_core_selection_OR_2.5_97.5_percentiles": np.quantile(
        single, [0.025, 0.975]
    ).tolist(),
    "single_core_interval_note": "Conditional sampling-sensitivity interval, NOT confidence interval. Each simulation randomly selects one observed pretreatment core per patient and retains outcome/HR; no patient resampling.",
    "paired_post_n": len(pair),
    "paired_delta_median": float(np.median(d)),
    "paired_delta_median_CI": np.quantile(boots, [0.025, 0.975]).tolist(),
    "paired_note": "Residual disease subset only, all non-pCR; no untreated comparator; bulk composition and tumor cellularity confound longitudinal change.",
    "bootstrap_replicates": B,
    "seed": 243375,
}
(r / "jci_statistics.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
