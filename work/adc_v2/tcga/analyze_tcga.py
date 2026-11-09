"""Analyze public TCGA-BRCA GE and IHC snapshots in this directory. Re-fetch with fetch_tcga.py.
Requires numpy, reportlab and pypdfium2; estimation does not use outcome data or clinical cutoffs.
"""

from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
import pypdfium2 as pdfium

here = Path(__file__).resolve().parent


def read(name):
    return json.loads((here / (name + ".json")).read_text(encoding="utf-8"))


def write_csv(name, data):
    with (here / name).open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=data[0].keys())
        w.writeheader()
        w.writerows(data)


def unique(records, key):
    result = {r[key]: r for r in records}
    assert len(result) == len(records), f"Duplicate {key} detected"
    return result


for src in read("sources"):
    raw = (here / src["file"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == src["sha256"]
samples = unique(read("samples"), "sampleId")
types = unique(read("sample_type"), "sampleId")
scores = unique(read("her2_ihc_score"), "patientId")
assert all(
    r["clinicalAttributeId"] == "HER2_IHC_SCORE" and r["value"] in list("0123")
    for r in scores.values()
)
molecular = {
    gene: unique(read(gene.lower()), "sampleId") for gene in ["ERBB2", "TACSTD2"]
}
for gene, records in molecular.items():
    for s, row in records.items():
        assert row["molecularProfileId"] == "brca_tcga_rna_seq_v2_mrna"
        assert row["entrezGeneId"] == {"ERBB2": 2064, "TACSTD2": 4070}[gene]
        assert samples[s]["patientId"] == row["patientId"]
        assert np.isfinite(row["value"]) and row["value"] >= 0
assert all(
    (r["sampleType"] == "Primary Solid Tumor") == (types[s]["value"] == "Primary")
    for s, r in samples.items()
)
primary = {s: r for s, r in samples.items() if r["sampleType"] == "Primary Solid Tumor"}
assert len(primary) == len({r["patientId"] for r in primary.values()})
complete = sorted(
    set(primary).intersection(*[set(records) for records in molecular.values()])
)
data = []
for s in complete:
    pt = samples[s]["patientId"]
    if pt not in scores:
        continue
    row = {"patient_id": pt, "sample_id": s, "her2_ihc_score": int(scores[pt]["value"])}
    for g, v in molecular.items():
        row[g + "_rsem"] = v[s]["value"]
        row[g + "_log2_rsem_plus1"] = float(np.log2(1 + v[s]["value"]))
    data.append(row)
flow = {
    "all_specimens": len(samples),
    "metastatic_excluded": len(samples) - len(primary),
    "primary_specimens": len(primary),
    "primary_complete_expression": len(complete),
    "missing_score_among_complete_expression": len(complete) - len(data),
    "analysis_patients": len(data),
    "total_recorded_ihc_scores": len(scores),
    "scored_patients_missing_primary_expression": len(scores) - len(data),
    "score_counts": dict(sorted(Counter(r["her2_ihc_score"] for r in data).items())),
}
write_csv("tcga_plot_data.csv", data)
scores_v = np.array([r["her2_ihc_score"] for r in data])


def average_rank(values):
    unique_values, inverse, counts = np.unique(
        values, return_inverse=True, return_counts=True
    )
    midranks = np.cumsum(counts) - (counts - 1) / 2.0
    return midranks[inverse]


def spearman(x, y):
    return float(np.corrcoef(average_rank(x), average_rank(y))[0, 1])


def superiority(a, b):
    # P(a > b) + 0.5 P(a == b), estimated from all between-group pairs.
    bs = np.sort(b)
    return float(
        np.mean(
            (np.searchsorted(bs, a, "left") + np.searchsorted(bs, a, "right"))
            / (2 * len(b))
        )
    )


summary = []
stats = {}
rng = np.random.default_rng(20261107)
# Use identical patient draws for both genes. Do not treat paired measurements as independent.
boots = rng.integers(0, len(data), size=(5000, len(data)))
for g in molecular:
    vals = np.array([r[g + "_log2_rsem_plus1"] for r in data])
    for score in range(4):
        a = vals[scores_v == score]
        q1, med, q3 = np.quantile(a, [0.25, 0.50, 0.75])
        summary.append(
            dict(gene=g, HER2_IHC_score=score, n=len(a), q25=q1, median=med, q75=q3)
        )
    rho = spearman(vals, scores_v)
    seps = superiority(vals[scores_v == 1], vals[scores_v == 0])
    bs_rho = []
    bs_sep = []
    for ix in boots:
        a = vals[ix]
        q = scores_v[ix]
        bs_rho.append(spearman(a, q))
        bs_sep.append(superiority(a[q == 1], a[q == 0]))
    stats[g] = {
        "rho": rho,
        "rho_ci": np.quantile(bs_rho, [0.025, 0.975]).tolist(),
        "p_1_gt_0": seps,
        "p_1_gt_0_ci": np.quantile(bs_sep, [0.025, 0.975]).tolist(),
    }
write_csv("tcga_summary.csv", summary)
(here / "tcga_results.json").write_text(
    json.dumps({"flow": flow, "stats": stats}, indent=2), encoding="utf8"
)
# Display source measurements as points. Horizontal displacement is jitter only.
W, H = 510, 278
c = canvas.Canvas(str(here / "tcga.pdf"), pagesize=(W, H))
c.setTitle("Primary breast-tumor transcripts by recorded HER2 IHC score (TCGA-BRCA)")


def text(x, y, s, pt=8, bold=False, col="#26383B", align="left"):
    c.setFont("Helvetica-Bold" if bold else "Helvetica", pt)
    c.setFillColor(HexColor(col))
    getattr(
        c,
        {
            "left": "drawString",
            "center": "drawCentredString",
            "right": "drawRightString",
        }[align],
    )(x, H - y, s)


def lin(x, y, x1, y1, col, w=0.45):
    c.setStrokeColor(HexColor(col))
    c.setLineWidth(w)
    c.line(x, H - y, x1, H - y1)


jy = np.random.default_rng(31811)
jitter = jy.uniform(-13, 13, len(data))
for j, g in enumerate(["ERBB2", "TACSTD2"]):
    left = 39 + 243 * j
    right = left + 211
    text(left - 16, 17, "AB"[j], 11, True)
    text(left, 17, g + " transcript", 11, True)
    info = stats[g]
    lo, hi = info["rho_ci"]
    text(left, 33, f"Spearman r = {info['rho']:.2f} (95% CI {lo:.2f} to {hi:.2f})", 8)
    top, bottom = 42, 216
    yy = lambda v: bottom - (v - 3) * (bottom - top) / 16
    for tick in [4, 8, 12, 16]:
        lin(left, yy(tick), right, yy(tick), "#DBE1E1")
        text(left - 5, yy(tick) + 2.5, str(tick), 7.5, align="right")
    lin(left, top, left, bottom, "#768486")
    lin(left, bottom, right, bottom, "#768486")
    color = ("#007C77", "#B95325")[j]
    for score in range(4):
        x = left + 24 + score * 53
        data_indexes = np.where(scores_v == score)[0]
        for i in data_indexes:
            v = data[i][g + "_log2_rsem_plus1"]
            c.saveState()
            c.setFillAlpha(0.43)
            c.setFillColor(HexColor(color))
            c.circle(x + jitter[i], H - yy(v), 1.25, stroke=0, fill=1)
            c.restoreState()
        vals = [data[i][g + "_log2_rsem_plus1"] for i in data_indexes]
        q1, med, q3 = np.quantile(vals, [0.25, 0.50, 0.75])
        # Median tick plus IQR at the center. No parametric error bars.
        lin(x, yy(q1), x, yy(q3), "#192D36", 1.8)
        lin(x - 4, yy(q1), x + 4, yy(q1), "#192D36", 0.85)
        lin(x - 4, yy(q3), x + 4, yy(q3), "#192D36", 0.85)
        lin(x - 6, yy(med), x + 6, yy(med), "#192D36", 1.9)
        text(x, 229, ["0", "1+", "2+", "3+"][score], 8.6, True, align="center")
        text(x, 241, "n = " + str(len(vals)), 7.5, col="#5C666A", align="center")
    text((left + right) / 2, 256, "Recorded HER2 IHC score", 8.4, align="center")
c.saveState()
c.translate(10, H - 128)
c.rotate(90)
c.setFont("Helvetica", 8.6)
c.setFillColor(HexColor("#26383B"))
c.drawCentredString(0, 0, "log2(normalized RSEM + 1)")
c.restoreState()
text(
    W / 2,
    273,
    "One primary tumor per patient; points are individual values.  Black bars: median and interquartile range.",
    7.4,
    col="#5C666A",
    align="center",
)
c.showPage()
c.save()
pdf = pdfium.PdfDocument(str(here / "tcga.pdf"))
pdf[0].render(scale=4.5).to_pil().save(here / "tcga.png")
print(json.dumps({"flow": flow, "stats": stats}, indent=2))
