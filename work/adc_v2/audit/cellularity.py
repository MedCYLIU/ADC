"""Repeat the matched-index experiment with tumor nuclei cellularity >=50%.
Usage: python cellularity.py source.xlsx --outdir .
"""

import argparse
import json
from pathlib import Path
from numbers import Real
import numpy as np
from audit_sampling import (
    read_data,
    matched_indices,
    interval_bootstrap,
    write_csv,
    TARGETS,
    N_BOOT,
    SEED,
)


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("source", type=Path)
    arg.add_argument("--outdir", type=Path, default=Path("."))
    opts = arg.parse_args()
    participants = read_data(opts.source)
    counts = []
    output = []
    excluded_missing = 0
    excluded_low = 0
    for donor, all_lesions in sorted(participants.items()):
        lesions = []
        for lesion in all_lesions:
            purity = lesion[3]
            if not isinstance(purity, Real) or not np.isfinite(purity):
                excluded_missing += 1
            elif purity < 50:
                excluded_low += 1
            else:
                assert purity <= 100
                lesions.append(lesion)
        r = []
        if len(lesions) >= 4:
            r = matched_indices(
                np.array([l[0] for l in lesions]),
                np.array([l[1] for l in lesions]),
                np.array([l[2] for l in lesions]),
            )
        for target in TARGETS:
            selected = [a for a in r if a["target"] == target]
            if not selected:
                continue
            rec = dict(patient=donor, target=target, eligible_indices=len(selected))
            for st in ["index", "same_pair", "cross_pair"]:
                rec[st] = np.mean([a[st] for a in selected])
            rec["same_minus_index"] = rec["same_pair"] - rec["index"]
            rec["cross_minus_index"] = rec["cross_pair"] - rec["index"]
            rec["cross_minus_same"] = rec["cross_pair"] - rec["same_pair"]
            output.append(rec)
        eligible = [a for a in r if a["target"] == TARGETS[0]]
        counts.append(
            dict(
                patient=donor,
                total=len(all_lesions),
                cellularity_pass=len(lesions),
                eligible_indices=len(eligible),
                triples=sum(a["triples"] for a in eligible),
            )
        )
    write_csv(opts.outdir / "cellularity_eligibility.csv", counts)
    write_csv(opts.outdir / "cellularity_donor_estimates.csv", output)
    included = [a for a in counts if a["eligible_indices"]]
    n = len(included)
    rng = np.random.default_rng(SEED)
    boot = rng.integers(0, n, (N_BOOT, n))
    rows = []
    for g in TARGETS:
        donors = [x for x in output if x["target"] == g]
        for z in [
            "index",
            "same_pair",
            "cross_pair",
            "same_minus_index",
            "cross_minus_index",
            "cross_minus_same",
        ]:
            values = [x[z] for x in donors]
            rows.append(
                dict(target=g, estimand=z, n=n, **interval_bootstrap(values, boot))
            )
    write_csv(opts.outdir / "cellularity_results.csv", rows)
    info = dict(
        source_donors=len(participants),
        source_lesions=sum(a["total"] for a in counts),
        cellularity_missing=excluded_missing,
        cellularity_under_50=excluded_low,
        lesions_with_cellularity_at_least_50=sum(a["cellularity_pass"] for a in counts),
        eligible_donors=n,
        lesions_in_eligible_donors=sum(a["cellularity_pass"] for a in included),
        eligible_index_lesions=sum(a["eligible_indices"] for a in included),
        triples=sum(a["triples"] for a in included),
    )
    (opts.outdir / "cellularity_manifest.json").write_text(
        json.dumps(info, indent=2), encoding="utf-8"
    )
    print(info)
    for a in rows:
        print(
            a["target"],
            a["estimand"],
            *[f"{a[k]:.3f}" for k in ("mean", "low", "high")],
        )


if __name__ == "__main__":
    main()
