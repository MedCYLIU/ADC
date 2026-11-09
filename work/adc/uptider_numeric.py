"""Reproduce the unrestricted UPTIDER specimen calculations.
Run in the directory containing source.xlsx. Outputs are CSV and JSON.
"""

import csv, hashlib, itertools, json
from pathlib import Path
import numpy as np
import openpyxl

DATA = Path("source.xlsx")
assert hashlib.md5(DATA.read_bytes()).hexdigest() == "6d7b721677e3ab6e83017cb1915396a5"
wb = openpyxl.load_workbook(DATA, read_only=True, data_only=True)
rows = iter(wb["Fig 4 - sup fig 6 - sup data 5"].values)
header = next(rows)
lesions = []
for i, r in enumerate(rows, 2):
    d = dict(zip(header, r))
    assert d["samplekind"] == "metastatic"
    site = "Lymph nodes" if d["tissue"].startswith("Lymph nodes") else d["tissue"]
    lesions.append(
        dict(
            patient=d["Patient_ID"],
            row=i,
            site=site,
            tissue=d["tissue"],
            HER2=d["HER2"],
            TROP2=d["TROP2"],
        )
    )
assert len(lesions) == 749


def write(name, records):
    with open(name, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=records[0])
        w.writeheader()
        w.writerows(records)


write("lesions.csv", lesions)
ids = sorted(set(r["patient"] for r in lesions))
groups = ["Single lesion", "Same-site pair", "Cross-site pair"]
assert len(ids) == 30
out = []
curves = []
summary = []
comp = []
counts = []
rng = np.random.default_rng(20261107)
boot = rng.integers(0, len(ids), (5000, len(ids)))
for target in ["HER2", "TROP2"]:
    patient_err = {}
    for ident in ids:
        sub = [r for r in lesions if r["patient"] == ident]
        v = np.array([r[target] for r in sub], dtype=float)
        sites = [r["site"] for r in sub]
        assert np.isfinite(v).all()
        idx = np.arange(len(v))
        allmed = np.median(v)
        e = {g: [] for g in groups}
        sens = {g: [] for g in groups}
        for i in idx:
            e[groups[0]].append(abs(v[i] - np.median(v[idx != i])))
            sens[groups[0]].append(abs(v[i] - allmed))
        for i, j in itertools.combinations(idx, 2):
            g = groups[1] if sites[i] == sites[j] else groups[2]
            sample = (v[i] + v[j]) / 2
            e[g].append(abs(sample - np.median(v[(idx != i) & (idx != j)])))
            sens[g].append(abs(sample - allmed))
        for g in groups:
            a = np.array(e[g])
            assert len(a) > 0
            patient_err[ident, g] = a
            out.append(
                dict(
                    patient=ident,
                    target=target,
                    strategy=g,
                    draws=len(a),
                    mae=a.mean(),
                    over_one=(a > 1).mean(),
                    full_median_mae=np.mean(sens[g]),
                )
            )
    for g in groups:
        ps = [patient_err[k, g] for k in ids]
        m = np.array([v.mean() for v in ps])
        h = np.array([(v > 1).mean() for v in ps])
        ci = lambda x: list(np.quantile(x[boot].mean(1), [0.025, 0.975]))
        summary.append(
            dict(
                target=target,
                strategy=g,
                n=len(ps),
                draws=sum(map(len, ps)),
                mae=m.mean(),
                mae_ci=ci(m),
                p=h.mean(),
                p_ci=ci(h),
            )
        )
        # Strict exceedance probability under the empirical distribution at every distinct observed error.
        values = np.concatenate(ps)
        weights = np.concatenate([np.full(len(v), 1 / (len(ids) * len(v))) for v in ps])
        uniq, inv = np.unique(values, return_inverse=True)
        tail = np.maximum(0, 1 - np.cumsum(np.bincount(inv, weights=weights)))
        for xx, yy in zip(uniq, tail):
            curves.append(dict(target=target, strategy=g, error=xx, probability=yy))
    for g in groups[1:]:
        ref = groups[0] if g == groups[1] else groups[1]
        dif = np.array(
            [patient_err[k, g].mean() - patient_err[k, ref].mean() for k in ids]
        )
        comp.append(
            dict(
                target=target,
                contrast=g + " minus " + ref,
                delta=dif.mean(),
                ci=ci(dif),
            )
        )
write("patient_metrics.csv", out)
write("threshold_curves.csv", curves)
Path("summary.json").write_text(
    json.dumps(dict(summary=summary, contrast=comp), indent=2)
)
