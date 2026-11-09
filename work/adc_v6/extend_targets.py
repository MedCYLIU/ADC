"""Exhaustive 72-target extension of the audited specimen-resampling analysis."""

from pathlib import Path
import csv, hashlib, importlib.util, json
import numpy as np
import openpyxl

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT = HERE / "all_targets"
OUT.mkdir(exist_ok=True)
SOURCE = BASE / "adc/source.xlsx"
spec = importlib.util.spec_from_file_location(
    "audit", BASE / "adc_v2/audit/audit_sampling.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
book = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
sheet = book["Fig 4 - sup fig 6 - sup data 5"]
records = iter(sheet.values)
fields = next(records)
targets = list(fields[7:-1])
assert len(targets) == 72
audit.TARGETS = targets
data = audit.read_data(SOURCE)
results = []
manifest = []
donor_rows = []


def high_cellularity(value):
    try:
        return float(value) >= 50
    except (TypeError, ValueError):
        return False


for subset in ["all", "tumor_nuclei_50"]:
    counts = []
    subset_donors = []
    for donor, raw in sorted(data.items()):
        rows = raw if subset == "all" else [r for r in raw if high_cellularity(r[3])]
        if len(rows) < 4:
            continue
        sites = np.array([r[0] for r in rows])
        y = np.array([r[1] for r in rows])
        numbers = np.array([r[2] for r in rows])
        scores = audit.matched_indices(sites, y, numbers)
        if not scores:
            continue
        anchor = [r for r in scores if r["target"] == targets[0]]
        counts.append(
            dict(
                donor=donor,
                records=len(rows),
                indices=len(anchor),
                triples=sum(r["triples"] for r in anchor),
            )
        )
        for target in targets:
            s = [r for r in scores if r["target"] == target]
            v = {
                k: float(np.mean([r[k] for r in s]))
                for k in ["index", "same_pair", "cross_pair"]
            }
            v["same_minus_index"] = v["same_pair"] - v["index"]
            v["cross_minus_index"] = v["cross_pair"] - v["index"]
            v["cross_minus_same"] = v["cross_pair"] - v["same_pair"]
            v["relative_same_reduction"] = (
                1 - v["same_pair"] / v["index"] if v["index"] > 0 else float("nan")
            )
            v["relative_cross_reduction"] = (
                1 - v["cross_pair"] / v["index"] if v["index"] > 0 else float("nan")
            )
            row = dict(subset=subset, donor=donor, target=target, **v)
            subset_donors.append(row)
            donor_rows.append(row)
        print(subset, donor, len(rows), len(anchor), flush=True)
    n = len(counts)
    boot = np.random.default_rng(20261107).integers(0, n, size=(5000, n))
    for target in targets:
        group = [r for r in subset_donors if r["target"] == target]
        for key in [
            "index",
            "same_pair",
            "cross_pair",
            "same_minus_index",
            "cross_minus_index",
            "cross_minus_same",
        ]:
            vals = np.array([r[key] for r in group])
            lo, hi = np.quantile(vals[boot].mean(axis=1), [0.025, 0.975])
            results.append(
                dict(
                    subset=subset,
                    target=target,
                    estimand=key,
                    donors=n,
                    mean=float(vals.mean()),
                    low=float(lo),
                    high=float(hi),
                )
            )
    manifest.append(
        dict(
            subset=subset,
            donors=n,
            records=sum(r["records"] for r in counts),
            eligible_indices=sum(r["indices"] for r in counts),
            triples=sum(r["triples"] for r in counts),
            by_donor=counts,
        )
    )
audit.write_csv(OUT / "target_estimates.csv", results)
audit.write_csv(OUT / "donor_estimates.csv", donor_rows)
(OUT / "manifest.json").write_text(
    json.dumps(
        dict(
            source_url=audit.SOURCE_URL,
            source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            targets=targets,
            bootstrap=5000,
            seed=20261107,
            populations=manifest,
            unit="specimen record; distinct lesion and collection time not verified",
            interval_scope="nominal percentile intervals, not simultaneous",
        ),
        indent=2,
    )
)

# The first two-target results must agree with the previously independent audit.
checks = []
for subset, filename in [
    ("all", "matched_results.csv"),
    ("tumor_nuclei_50", "cellularity_results.csv"),
]:
    old = list(csv.DictReader((BASE / "adc_v2/audit" / filename).open()))
    for row in old:
        r = next(
            v
            for v in results
            if v["subset"] == subset
            and v["target"] == row["target"]
            and v["estimand"] == row["estimand"]
        )
        error = max(abs(r[k] - float(row[k])) for k in ["mean", "low", "high"])
        assert error < 1e-10, (subset, row["target"], row["estimand"], error)
        checks.append(
            dict(
                subset=subset,
                target=row["target"],
                estimand=row["estimand"],
                max_error=error,
            )
        )
(OUT / "validation.json").write_text(json.dumps(checks, indent=2))
print("Saved 72 targets, both subsets; established HER2/TROP2 results reproduced.")
