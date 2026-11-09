"""GSE209998: paired expression in breast primary tumors and matched metastases.
Run after download.py; upper-quartile normalized counts were deposited by the authors.
"""

from pathlib import Path
import csv, gzip, json, hashlib, math, collections, itertools
import numpy as np
import pandas as pd
from download import GEO, MATRIX, META
from collections import Counter

GENES = ("ERBB2", "TACSTD2")


def main():
    with gzip.open(META, "rt", encoding="utf-8") as f:
        meta = {}
        characteristics = {}
        for row in csv.reader(f, delimiter="\t"):
            if not row:
                continue
            key = row[0]
            if key == "!Sample_characteristics_ch1":
                field = row[1].split(": ", 1)[0]
                vals = [x.split(": ", 1) for x in row[1:]]
                assert all(len(v) == 2 and v[0] == field for v in vals)
                characteristics[field] = [v[1] for v in vals]
            elif key == "!Sample_description":
                if row[1].startswith("AUR-"):
                    meta["barcode"] = row[1:]
                else:
                    meta["preservation"] = row[1:]
            elif key == "!Sample_geo_accession":
                meta["gsm"] = row[1:]
            elif key == "!Sample_source_name_ch1":
                meta["source"] = row[1:]
    ann = pd.DataFrame(meta | characteristics)
    ann["patient"] = ann.barcode.str.extract(r"^(AUR-[A-Z0-9]+)-", expand=False)
    assert (
        len(ann) == 129
        and len(set(ann.barcode)) == len(ann)
        and len(set(ann.gsm)) == len(ann)
        and ann.patient.notna().all()
        and (ann.source == ann.disease).all()
    )
    with gzip.open(MATRIX, "rt", encoding="utf-8") as f:
        read = csv.reader(f, delimiter="\t")
        cols = next(read)[1:]
        assert set(cols) == set(ann.barcode) and len(cols) == len(ann)
        expr = {}
        counts = Counter()
        for row in read:
            if row[0] in GENES:
                counts[row[0]] += 1
                assert len(row) == len(cols) + 1
                expr[row[0]] = dict(zip(cols, map(float, row[1:])))
    assert all(counts[g] == 1 for g in GENES)
    for g in GENES:
        ann[g + "_UQN"] = ann.barcode.map(expr[g])
        assert np.isfinite(ann[g + "_UQN"]).all() and (ann[g + "_UQN"] >= 0).all()
        ann[g] = np.log2(ann[g + "_UQN"] + 1)
    ann.to_csv("samples_audit.csv", index=False)
    rows = []
    for pid, group in ann.groupby("patient"):
        primary = group[group.source == "Primary tumor"]
        met = group[group.source == "Metastatic tumor"]
        if primary.empty or met.empty:
            continue
        assert len(primary) == 1
        for g in GENES:
            a = primary[g].iat[0]
            b = met[g].median()
            match = met[met.preservation == primary.preservation.iat[0]]
            rows.append(
                dict(
                    patient=pid,
                    gene=g,
                    primary=a,
                    metastatic_median=b,
                    delta=b - a,
                    primary_gsm=primary.gsm.iat[0],
                    primary_treatment=primary.treatment.iat[0],
                    primary_preservation=primary.preservation.iat[0],
                    n_met=len(met),
                    n_sites=met.tissue.nunique(),
                    met_gsm=";".join(met.gsm),
                    same_preservation_n=len(match),
                    same_preservation_delta=(
                        match[g].median() - a if len(match) else None
                    ),
                )
            )
    data = pd.DataFrame(rows)
    data.to_csv("patient_paired.csv", index=False)
    pids = sorted(set(data.patient))
    n = len(pids)
    # donor bootstraps (paired differences) with same random indices for both targets
    rng = np.random.default_rng(20261107)
    draw = rng.integers(n, size=(10000, n))
    report = {}
    for g in GENES:
        group = data[data.gene == g].sort_values("patient")
        v = group.delta.to_numpy()
        report[g] = {
            "n_donors": len(v),
            "n_lower": int((v < 0).sum()),
            "n_positive": int((v > 0).sum()),
            "n_tied": int((v == 0).sum()),
            "n_abs_delta_ge_1": int((abs(v) >= 1).sum()),
            "median_delta": float(np.median(v)),
            "median_delta_ci": np.quantile(
                np.median(v[draw], axis=1), [0.025, 0.975]
            ).tolist(),
            "median_absolute_delta": float(np.median(abs(v))),
            "median_absolute_delta_ci": np.quantile(
                np.median(abs(v[draw]), axis=1), [0.025, 0.975]
            ).tolist(),
            "iqr": np.quantile(v, [0.25, 0.75]).tolist(),
        }
        for name, col in [
            (
                "pretreatment_primary",
                group.loc[group.primary_treatment == "Pre-treatment", "delta"],
            ),
            ("same_preservation", group["same_preservation_delta"].dropna()),
        ]:
            sv = col.to_numpy()
            srng = np.random.default_rng(20261108)
            sdraw = srng.integers(len(sv), size=(10000, len(sv)))
            report[g][name] = {
                "n": len(col),
                "median_delta": float(col.median()),
                "median_absolute_delta": float(col.abs().median()),
                "median_delta_ci": np.quantile(
                    np.median(sv[sdraw], axis=1), [0.025, 0.975]
                ).tolist(),
            }
    paired = ann[ann.patient.isin(pids)]
    info = {
        "source": GEO,
        "expression": MATRIX,
        "metadata": META,
        "sha256": {
            file: hashlib.sha256(Path(file).read_bytes()).hexdigest()
            for file in (MATRIX, META)
        },
        "n_GEO_samples": len(ann),
        "n_all_patients": ann.patient.nunique(),
        "counts": dict(Counter(ann.source)),
        "n_paired_primary": int((paired.source == "Primary tumor").sum()),
        "n_paired_metastatic": int((paired.source == "Metastatic tumor").sum()),
        "n_paired_donors": n,
        "n_paired_with_multiple_metastases": int(
            (data[data.gene == GENES[0]].n_met > 1).sum()
        ),
        "results": report,
    }
    Path("results.json").write_text(
        json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
