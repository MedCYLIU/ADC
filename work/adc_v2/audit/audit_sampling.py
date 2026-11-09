"""Exact UPTIDER enumerations. Usage: python audit_sampling.py source.xlsx --outdir .
No access to clinical identifiers or expression files beyond the public workbook.
"""

import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from openpyxl import load_workbook

SOURCE_URL = "https://pmc-oa-opendata.s3.amazonaws.com/PMC12852894.1/41467_2025_67840_MOESM13_ESM.xlsx"
SHA256 = "0c82f309085bbc0709d681e1bd169438c2593687c9f0f1d92f7698fe8a6b3ca2"
TARGETS = ["HER2", "TROP2"]
STRATEGIES = ["Index alone", "Same-site second lesion", "Cross-site second lesion"]
SEED = 20261107
N_BOOT = 5000


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_data(workbook):
    assert hashlib.sha256(workbook.read_bytes()).hexdigest() == SHA256
    book = load_workbook(workbook, read_only=True, data_only=True)
    sheet = book["Fig 4 - sup fig 6 - sup data 5"]
    values = iter(sheet.values)
    fields = next(values)
    participants = defaultdict(list)
    for row_number, values in enumerate(values, 2):
        row = dict(zip(fields, values))
        assert row["samplekind"] == "metastatic" and row["tissue"] and row["Patient_ID"]
        tissue = row["tissue"]
        site = "Lymph nodes" if tissue.startswith("Lymph nodes ") else tissue
        expression = [row[target] for target in TARGETS]
        assert np.isfinite(expression).all()
        participants[row["Patient_ID"]].append(
            (site, expression, row_number, row["tumor_cellularity_nuclei"])
        )
    assert len(participants) == 30 and sum(map(len, participants.values())) == 749
    return participants


def primary_errors(sites, y):
    """Independent vectorization; medians formed explicitly from unselected rows."""
    combinations = {
        "Index alone": [(i,) for i in range(len(y))],
        "Same-site second lesion": [],
        "Cross-site second lesion": [],
    }
    for i, j in itertools.combinations(range(len(y)), 2):
        which = (
            "Same-site second lesion"
            if sites[i] == sites[j]
            else "Cross-site second lesion"
        )
        combinations[which].append((i, j))
    output = []
    for which, draws in combinations.items():
        selected = np.asarray(draws)
        mask = (np.arange(len(y))[None, :, None] == selected[:, None, :]).any(axis=2)
        reference = np.nanmedian(
            np.where(mask[:, :, None], np.nan, y[None, :, :]), axis=1
        )
        absolute = np.abs(y[selected].mean(axis=1) - reference)
        for k, target in enumerate(TARGETS):
            output.append(
                dict(
                    target=target,
                    strategy=which,
                    draws=len(draws),
                    mae=absolute[:, k].mean(),
                    exceed_one=(absolute[:, k] > 1).mean(),
                )
            )
    return output


def matched_indices(sites, y, row_numbers):
    """The two candidate additions are always compared to exactly the same reference.
    Excluding the three lesions prevents leakage from either option into the reference.
    """
    indices = np.arange(len(y))
    scores = []
    for i in indices:
        same = indices[(sites == sites[i]) & (indices != i)]
        cross = indices[sites != sites[i]]
        if not len(same) or not len(cross):
            continue
        j, k = np.asarray(list(itertools.product(same, cross))).T
        excluded = (
            (indices[None, :] == i)
            | (indices[None, :] == j[:, None])
            | (indices[None, :] == k[:, None])
        )
        reference = np.nanmedian(
            np.where(excluded[:, :, None], np.nan, y[None, :, :]), axis=1
        )
        scores_by_strategy = [
            np.abs(y[i] - reference),
            np.abs((y[i] + y[j]) / 2 - reference),
            np.abs((y[i] + y[k]) / 2 - reference),
        ]
        for g, target in enumerate(TARGETS):
            scores.append(
                dict(
                    index_row=row_numbers[i],
                    index_site=sites[i],
                    same=len(same),
                    cross=len(cross),
                    triples=len(j),
                    target=target,
                    index=scores_by_strategy[0][:, g].mean(),
                    same_pair=scores_by_strategy[1][:, g].mean(),
                    cross_pair=scores_by_strategy[2][:, g].mean(),
                )
            )
    return scores


def interval_bootstrap(values, bootstrap_draws):
    probs = np.quantile(
        np.asarray(values)[bootstrap_draws].mean(axis=1), [0.025, 0.975]
    )
    return dict(mean=np.mean(values), low=probs[0], high=probs[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("."))
    args = parser.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    data = read_data(args.source)
    base = []
    all_indices = []
    patient_records = []
    match_counts = []
    for ident, lesions in sorted(data.items()):
        sites = np.array([r[0] for r in lesions])
        y = np.array([r[1] for r in lesions])
        rows = np.array([r[2] for r in lesions])
        base.extend(dict(patient=ident, **row) for row in primary_errors(sites, y))
        index_results = matched_indices(sites, y, rows)
        all_indices.extend(dict(patient=ident, **row) for row in index_results)
        unique_anchors = [r for r in index_results if r["target"] == TARGETS[0]]
        match_counts.append(
            dict(
                patient=ident,
                lesions=len(lesions),
                eligible_anchors=len(unique_anchors),
                matched_triples=sum(r["triples"] for r in unique_anchors),
            )
        )
        for target in TARGETS:
            subset = [r for r in index_results if r["target"] == target]
            assert len(subset) > 0
            record = dict(patient=ident, target=target, anchors=len(subset))
            for key in ["index", "same_pair", "cross_pair"]:
                record[key] = np.mean([r[key] for r in subset])
            record["cross_minus_same"] = record["cross_pair"] - record["same_pair"]
            record["same_minus_index"] = record["same_pair"] - record["index"]
            record["cross_minus_index"] = record["cross_pair"] - record["index"]
            patient_records.append(record)
    write_csv(out / "primary_donor_reproduction.csv", base)
    write_csv(out / "matched_index_estimates.csv", all_indices)
    write_csv(out / "matched_donor_estimates.csv", patient_records)
    write_csv(out / "matched_eligibility.csv", match_counts)
    rng = np.random.default_rng(SEED)
    boot = rng.integers(0, len(data), size=(N_BOOT, len(data)))
    primary = []
    matched = []
    for target in TARGETS:
        for strategy in STRATEGIES:
            values = [
                r for r in base if r["target"] == target and r["strategy"] == strategy
            ]
            for metric in ("mae", "exceed_one"):
                primary.append(
                    dict(
                        target=target,
                        strategy=strategy,
                        draws=sum(r["draws"] for r in values),
                        metric=metric,
                        **interval_bootstrap([r[metric] for r in values], boot)
                    )
                )
        donors = [r for r in patient_records if r["target"] == target]
        for key in [
            "index",
            "same_pair",
            "cross_pair",
            "same_minus_index",
            "cross_minus_index",
            "cross_minus_same",
        ]:
            vals = np.array([r[key] for r in donors])
            loo = (vals.sum() - vals) / (len(vals) - 1)
            matched.append(
                dict(
                    target=target,
                    estimand=key,
                    n=len(vals),
                    **interval_bootstrap(vals, boot),
                    loo_low=min(loo),
                    loo_high=max(loo)
                )
            )
    write_csv(out / "primary_reproduction.csv", primary)
    write_csv(out / "matched_results.csv", matched)
    summary = dict(
        source=SOURCE_URL,
        sha256=SHA256,
        donors=len(data),
        metastatic_lesions=sum(r["lesions"] for r in match_counts),
        eligible_anchors=sum(r["eligible_anchors"] for r in match_counts),
        matched_triples=sum(r["matched_triples"] for r in match_counts),
        seed=SEED,
        bootstrap_replicates=N_BOOT,
    )
    (out / "manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(summary)
    for r in matched:
        print(
            r["target"],
            r["estimand"],
            [round(r[k], 4) for k in ("mean", "low", "high")],
        )


if __name__ == "__main__":
    main()
