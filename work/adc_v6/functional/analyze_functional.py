"""Descriptive secondary calculations; no dose-response refitting or inferential tests."""

from pathlib import Path
from openpyxl import load_workbook
import numpy as np, pandas as pd, json, re, sys, hashlib

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parent
files = {
    i: ROOT / f"ccr-25-3321_supplementary_table_{i}_suppst_{i}.xlsx"
    for i in range(1, 6)
}
rows = []
ws = load_workbook(files[1], data_only=True).active
cell = None
for r in range(5, 29):
    cell = ws.cell(r, 2).value or cell
    for protein, col in [("beta-actin", 4), ("HER2", 6), ("TROP2", 8), ("BCRP", 10)]:
        rows.append(
            dict(
                cell_line=cell,
                replicate=ws.cell(r, 3).value,
                protein=protein,
                deposited_ratio=ws.cell(r, col).value,
                source_row=r,
                source_col=col,
            )
        )
df = pd.DataFrame(rows)
df.to_csv(ROOT / "immunoblot_deposited_ratios.csv", index=False)
ag = (
    df.groupby(["cell_line", "protein"])
    .deposited_ratio.agg(["count", "mean", "std"])
    .reset_index()
)
ag.to_csv(ROOT / "immunoblot_descriptive_summary.csv", index=False)
ws = load_workbook(files[3], data_only=True).active


def value(x):
    m = re.match(r"\s*(>?\s*)([0-9.]+)\s*nM\s*\[(.*?)\]", x)
    assert m, x
    return float(m[2]), bool(m[1].strip()), m[3]


res = []
for start, end, models in [
    (4, 9, [("HCC1806", 3, 4, 5), ("MDA-MB-468", 6, 7, 8)]),
    (12, 19, [("HCC1954", 3, 4, 5), ("SUM190", 6, 7, 8)]),
]:
    for r in range(start, end + 1):
        for cell, pc, rc, pv in models:
            p = ws.cell(r, pc).value
            q = ws.cell(r, rc).value
            a, ac, aci = value(p)
            b, bc, bci = value(q)
            assert not ac
            res.append(
                dict(
                    model=cell,
                    selection="SG" if cell in ["HCC1806", "MDA-MB-468"] else "T-DXd",
                    drug=ws.cell(r, 2).value,
                    parental_IC50_nM=a,
                    resistant_IC50_nM=b,
                    resistant_censored_above=bc,
                    parental_reported_CI=aci,
                    resistant_reported_CI=bci,
                    parental_source=p,
                    resistant_source=q,
                    source_p_value=ws.cell(r, pv).value,
                    resistance_ratio=b / a,
                    ratio_is_lower_bound=bc,
                    log2_resistance_ratio=float(np.log2(b / a)),
                    source_row=r,
                )
            )
ic = pd.DataFrame(res)
assert len(ic) == 28
ic.to_csv(ROOT / "ic50_descriptive_ratios.csv", index=False)
comp = []
for cell in ic.model.unique():
    s = ic[ic.model.eq(cell)].set_index("drug")
    for adc, payload in [
        ("Sacituzumab Govitecan", "SN-38"),
        ("T-DXd", "DXd"),
        ("Trastuzumab-MMAF", "MMAF"),
        ("Sacituzumab-MMAF", "MMAF"),
    ]:
        if adc in s.index:
            comp.append(
                dict(
                    model=cell,
                    ADC=adc,
                    payload=payload,
                    ADC_resistance_ratio=float(s.loc[adc, "resistance_ratio"]),
                    free_payload_resistance_ratio=float(
                        s.loc[payload, "resistance_ratio"]
                    ),
                    ratio_of_ratios=float(
                        s.loc[adc, "resistance_ratio"]
                        / s.loc[payload, "resistance_ratio"]
                    ),
                    ratio_is_lower_bound=bool(s.loc[adc, "ratio_is_lower_bound"]),
                    interpretation="Descriptive decomposition only; does not estimate delivery loss or causal mediation.",
                )
            )
pd.DataFrame(comp).to_csv(ROOT / "adc_payload_descriptive_ratios.csv", index=False)
audit = {
    "data_retrieved": "2026-11-09",
    "original_article_DOI": "10.1158/1078-0432.CCR-25-3321",
    "original_article_source_id": 18,
    "study_count": 1,
    "parental_resistant_model_pairs": 4,
    "normalized_immunoblot_rows": 24,
    "protein_values": 96,
    "IC50_paired_model_drug_rows": 28,
    "censored_IC50_ratios": 1,
    "inference_performed": False,
    "raw_dose_response_available_in_tables": False,
    "individual_mouse_data_available_in_tables": False,
    "table2_excluded": "Literature summary, not new clinical observations; overlaps cohorts already cited",
    "tables4_5_excluded": "Published ANOVA/Tukey output only; underlying animal-level longitudinal data absent",
    "immunoblot_limitation": "Rep1-3 are assay repetitions within four selected model pairs, not twelve independent cell models. Beta-actin included as deposited; no undocumented renormalization.",
    "IC50_limitation": "Published fitted values transcribed and ratios calculated; 95% CIs are original curve-fit intervals, frequently NE. No CI for ratio calculated because fits/covariance/raw data unavailable. Extreme ratios may reflect unstable parental fits. No pooled resistance estimate or statistical significance claimed.",
    "source_files_sha256": {
        f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files.values()
    },
}
(ROOT / "functional_analysis_audit.json").write_text(
    json.dumps(audit, indent=2), encoding="utf-8"
)
print(
    ic[["model", "drug", "resistance_ratio", "ratio_is_lower_bound"]].to_string(
        index=False
    )
)
