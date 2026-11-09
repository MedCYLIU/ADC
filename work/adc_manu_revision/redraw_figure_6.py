"""Redraw the frozen Figure 6 measurements; no statistical fitting is performed."""

from pathlib import Path
import hashlib, json, math, shutil
import pandas as pd
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
import pypdfium2 as pdfium

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
OUT = BASE / "figures"
OUT.mkdir(exist_ok=True)
SRC = ROOT / "work/adc_v6/public_cohorts"
files = [
    "jci_sample_joined.csv",
    "jci_patient_baseline.csv",
    "jci_paired_change.csv",
    "jci_models.csv",
]
manifest = {
    name: hashlib.sha256((SRC / name).read_bytes()).hexdigest() for name in files
}
for n, f in [("Arial", "arial.ttf"), ("Arial-Bold", "arialbd.ttf")]:
    pdfmetrics.registerFont(TTFont(n, str(Path("C:/Windows/Fonts") / f)))
BLUE = "#2F5D7C"
GREEN = "#4C8C82"
GREY = "#7A7A7A"
INK = "#222222"
W = 180 / 25.4 * 72
H = 516
c = canvas.Canvas(str(OUT / "Figure_6.pdf"), pagesize=(W, H))
c.setTitle("Pretreatment specimen variation and response associations in NCT02326974")
c.setAuthor(
    "Chuyu Liu, Qi Wang, Kaiguo Long, Yujing Zhang, Feng Zeng, Xiaoyu Wang, Yulan Cai"
)


def txt(x, y, t, size=9.5, bold=False, col=INK, align="left"):
    c.setFillColor(HexColor(col))
    c.setFont("Arial-Bold" if bold else "Arial", size)
    {"left": c.drawString, "center": c.drawCentredString, "right": c.drawRightString}[
        align
    ](x, H - y, str(t))


def line(x, y, X, Y, col=INK, width=0.7, dash=()):
    c.setStrokeColor(HexColor(col))
    c.setLineWidth(width)
    c.setDash(dash)
    c.line(x, H - y, X, H - Y)
    c.setDash([])


def dot(x, y, r=2, col=INK, open=False):
    c.setStrokeColor(HexColor(col))
    c.setFillColor(HexColor("#FFFFFF" if open else col))
    c.setLineWidth(0.7)
    c.circle(x, H - y, r, fill=1, stroke=1)


def square(x, y, r=2.8, col=BLUE):
    c.setStrokeColor(HexColor(col))
    c.setFillColor(HexColor(col))
    c.rect(x - r, H - y - r, 2 * r, 2 * r, fill=1, stroke=0)


def rot(x, y, t):
    c.saveState()
    c.translate(x, H - y)
    c.rotate(90)
    c.setFont("Arial", 9.5)
    c.setFillColor(HexColor(INK))
    c.drawCentredString(0, 0, t)
    c.restoreState()


s = pd.read_csv(SRC / "jci_sample_joined.csv")
p = pd.read_csv(SRC / "jci_patient_baseline.csv")
q = p[p.ncores == 2]
pairs = []
for pid in q["Patient ID"]:
    z = s[(s["Patient ID"] == pid) & (s.Treatment_Group == "Pre")].sort_values(
        "expression_id"
    )
    pairs.append((z.ERBB2.iloc[0], z.ERBB2.iloc[1], z.pCR.iloc[0]))
assert len(pairs) == 113
left, top, right, bottom = 48, 53, 231, 235
x = lambda v: left + (v - 5) / 10 * (right - left)
y = lambda v: bottom - (v - 5) / 10 * (bottom - top)
txt(4, 16, "A", 12, True)
txt(22, 16, "Paired baseline cores", 10.5, True)
txt(48, 36, "113 patients", 9.5)
for v in [6, 8, 10, 12, 14]:
    line(x(v), bottom, x(v), bottom + 3)
    txt(x(v), bottom + 17, v, align="center")
    line(left - 3, y(v), left, y(v))
    txt(left - 7, y(v) + 3, v, align="right")
line(left, top, left, bottom)
line(left, bottom, right, bottom)
line(x(5), y(5), x(15), y(15), "#A6A6A6", 0.7, (3, 2))
for a, b, resp in pairs:
    dot(x(a), y(b), 2, BLUE if resp == "Yes" else GREY, resp == "No")
txt((left + right) / 2, 270, "Specimen A ERBB2 (log2 CPM)", 9.5, align="center")
rot(13, (top + bottom) / 2, "Specimen B ERBB2 (log2 CPM)")
dot(52, 288, 2.1, BLUE)
txt(60, 291, "pCR")
dot(111, 288, 2.1, GREY, True)
txt(119, 291, "No pCR")
txt(48, 309, "Median absolute difference: 0.487", 9)

L, T, R, B = 302, 53, 493, 235
yy = lambda v: B - (v - 4) / 12 * (B - T)
xx = [337, 455]
txt(265, 16, "B", 12, True)
txt(283, 16, "Baseline and residual tumors", 10.5, True)
txt(302, 36, "43 patients without pCR", 9.5)
for v in [4, 7, 10, 13, 16]:
    line(L - 3, yy(v), L, yy(v))
    txt(L - 7, yy(v) + 3, v, align="right")
line(L, T, L, B)
line(L, B, R, B)
data = pd.read_csv(SRC / "jci_paired_change.csv")
rng = np.random.default_rng(6)
assert len(data) == 43
for row in data.itertuples():
    jit = rng.uniform(-6, 6)
    line(xx[0] + jit, yy(row.ERBB2), xx[1] + jit, yy(row.post_ERBB2), "#BBBBBB", 0.5)
    dot(xx[0] + jit, yy(row.ERBB2), 1.8, GREY)
    dot(xx[1] + jit, yy(row.post_ERBB2), 1.8, GREEN)
for a, label in zip(xx, ["Pretreatment", "Residual"]):
    txt(a, B + 17, label, 9.5, align="center")
rot(267, (T + B) / 2, "ERBB2 (log2 CPM)")
txt(302, 285, "Median paired change: +1.528", 9)
txt(302, 303, "95% CI 0.383 to 2.404", 9)

m = pd.read_csv(SRC / "jci_models.csv")
txt(4, 341, "C", 12, True)
txt(22, 341, "ERBB2 and pCR", 10.5, True)
txt(493, 341, "Adjusted for hormone receptor status", 9, align="right")
txt(20, 365, "Baseline expression summary", 9.5, True)
txt(222, 365, "n", 9.5, True, align="right")
txt(499, 365, "OR (95% CI)", 9.5, True, align="right")
fl, fr = 244, 375
fx = lambda v: fl + math.log(v / 0.7) / math.log(5 / 0.7) * (fr - fl)
line(fx(1), 371, fx(1), 462, "#AAAAAA", 0.7, (3, 2))
labels = [
    "All baseline patients: mean",
    "Two-core subset: mean",
    "Two-core subset: lower",
    "Two-core subset: higher",
]
for i, row in enumerate(m.itertuples()):
    yv = 384 + 24 * i
    col = BLUE if i == 0 else GREEN
    txt(20, yv + 3, labels[i], 9.5)
    txt(222, yv + 3, row.n, 9.5, align="right")
    line(fx(row.lo), yv, fx(row.hi), yv, col, 1)
    square(fx(row.OR), yv, 2.8, col)
    txt(499, yv + 3, f"{row.OR:.2f} ({row.lo:.2f}–{row.hi:.2f})", 9.2, align="right")
line(fl, 472, fr, 472)
for v in [1, 2, 4]:
    line(fx(v), 472, fx(v), 476)
    txt(fx(v), 490, v, 9.5, align="center")
txt((fl + fr) / 2, 509, "Odds ratio (log scale)", 9.5, align="center")
c.save()
f = pdfium.PdfDocument(OUT / "Figure_6.pdf")
f[0].render(scale=600 / 72).to_pil().save(OUT / "Figure_6.png", dpi=(600, 600))
f[0].render(scale=1.7).to_pil().save(OUT / "Figure_6_preview.png")
assert manifest == {
    name: hashlib.sha256((SRC / name).read_bytes()).hexdigest() for name in files
}
(BASE / "figure6_manifest.json").write_text(
    json.dumps(
        {
            "source_sha256": manifest,
            "panel_a_points": 113,
            "panel_b_pairs": 43,
            "models": m.to_dict(orient="records"),
            "new_model_fits": False,
            "width_mm": 180,
            "height_points": H,
            "font_points": {
                "body": 9.5,
                "minor": 9,
                "panel_title": 10.5,
                "panel_label": 12,
            },
            "jitter_seed": 6,
            "palette": [BLUE, GREEN, GREY],
        },
        indent=2,
    ),
    encoding="utf-8",
)
print(
    "Figure 6 generated from frozen sources; 113 points, 43 pairs, 4 model intervals."
)
