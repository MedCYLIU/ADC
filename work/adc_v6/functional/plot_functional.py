from pathlib import Path
import os, math, json
import numpy as np, pandas as pd, pypdfium2 as pdfium
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, Color, white

ROOT = Path(__file__).resolve().parent
for face, file in [("Arial", "arial.ttf"), ("ArialBold", "arialbd.ttf")]:
    pdfmetrics.registerFont(
        TTFont(face, str(Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / file))
    )
W, H = 180 / 25.4 * 72, 407
out = ROOT / "Supplementary_Figure_2.pdf"
c = canvas.Canvas(str(out), pagesize=(W, H))
c.setTitle("Descriptive reanalysis of deposited ADC-resistant cell-model measurements")


def text(x, y, s, size=8, bold=False, align="left", color="#111111"):
    c.setFillColor(HexColor(color))
    c.setFont("ArialBold" if bold else "Arial", size)
    getattr(
        c,
        {
            "left": "drawString",
            "right": "drawRightString",
            "center": "drawCentredString",
        }[align],
    )(x, H - y, str(s))


def line(x, y, x2, y2, color="#111111", width=0.55, dash=None):
    c.setStrokeColor(HexColor(color))
    c.setLineWidth(width)
    c.setDash(dash or [])
    c.line(x, H - y, x2, H - y2)
    c.setDash([])


def dot(x, y, fill=True, color="#255C88", r=2):
    c.setStrokeColor(HexColor(color))
    c.setFillColor(HexColor(color) if fill else white)
    c.setLineWidth(0.6)
    c.circle(x, H - y, r, stroke=1, fill=1)


def square(x, y):
    c.setStrokeColor(HexColor("#111111"))
    c.setFillColor(white)
    c.rect(x - 2.5, H - y - 2.5, 5, 5, fill=1, stroke=1)


models = ["HCC1806", "MDA-MB-468", "HCC1954", "SUM190"]
labels = [
    "HCC1806 / SG-R",
    "MDA-MB-468 / SG-R",
    "HCC1954 / T-DXd-R",
    "SUM190 / T-DXd-R",
]
raw = pd.read_csv(ROOT / "immunoblot_deposited_ratios.csv")
ic = pd.read_csv(ROOT / "ic50_descriptive_ratios.csv")
text(3, 13, "A", 11, True)
text(18, 13, "Deposited protein ratios", 9, True)
square(270, 10)
text(277, 13, "Parental normalizer", 7)
dot(390, 10)
text(397, 13, "Resistant replicates", 7)
xstart = 133
colw = 117
for j, (protein, maxx, ticks) in enumerate(
    [
        ("HER2", 2, [0, 1, 2]),
        ("TROP2", 1.5, [0, 0.5, 1, 1.5]),
        ("BCRP", 7, [0, 2, 4, 6]),
    ]
):
    x = xstart + j * colw
    pw = 102
    sx = lambda z: x + float(z) / maxx * pw
    text(x + pw / 2, 34, protein, 8.5, True, "center")
    line(sx(1), 41, sx(1), 131, "#BBBBBB", 0.4, [2, 2])
    for i, model in enumerate(models):
        y = 52 + i * 23
        parental = raw[
            raw.cell_line.eq(model) & raw.protein.eq(protein)
        ].deposited_ratio.to_numpy()
        resistant = raw[
            raw.cell_line.str.startswith(model + "-") & raw.protein.eq(protein)
        ].deposited_ratio.to_numpy()
        assert len(parental) == len(resistant) == 3
        if j == 0:
            text(3, y + 2.5, labels[i], 7.3)
        line(sx(parental.mean()), y, sx(resistant.mean()), y, "#A6A6A6", 0.6)
        square(sx(parental.mean()), y)
        for delta, val in zip([-3, 0, 3], resistant):
            dot(sx(val), y + delta)
    line(x, 135, x + pw, 135)
    for t in ticks:
        line(sx(t), 135, sx(t), 138)
        text(sx(t), 148, f"{t:g}", 7, align="center")
text(310, 162, "Ratio relative to deposited parental reference", 8, align="center")
text(3, 184, "B", 11, True)
text(18, 184, "IC50 resistant / parental: matched MMAF comparisons", 9, True)
drugs = ["MMAF", "Sacituzumab-MMAF", "Trastuzumab-MMAF"]
labels2 = ["Free MMAF", "Sacituzumab-MMAF", "Trastuzumab-MMAF"]
cx, cellw = 132, 88.5
top = 225
cellh = 27
for j, model in enumerate(models):
    text(cx + (j + 0.5) * cellw, 204, model, 7.5, True, "center")
    text(
        cx + (j + 0.5) * cellw,
        216,
        "SG-selected" if j < 2 else "T-DXd-selected",
        7,
        align="center",
        color="#555555",
    )


def color(v):
    # A symmetric neutral anchor at zero with fixed -2 and +8 limits.
    base = np.array([247, 247, 247.0])
    end = np.array([43, 93, 132.0]) if v < 0 else np.array([156, 46, 49.0])
    t = min(abs(v) / (2 if v < 0 else 8), 1)
    return Color(*((base * (1 - t) + end * t) / 255))


for i, drug in enumerate(drugs):
    text(3, top + (i + 0.5) * cellh + 2.5, labels2[i], 8)
    for j, model in enumerate(models):
        x, y = cx + j * cellw, top + i * cellh
        a = ic[ic.drug.eq(drug) & ic.model.eq(model)]
        c.setFillColor(
            color(float(a.iloc[0].log2_resistance_ratio))
            if len(a)
            else HexColor("#F5F5F5")
        )
        c.rect(x, H - y - cellh, cellw - 0.8, cellh - 0.8, fill=1, stroke=0)
        if len(a):
            a = a.iloc[0]
            v = a.resistance_ratio
            label = f"{v:.2f}" if v < 10 else f"{v:.1f}"
            if drug == "Trastuzumab-MMAF":
                label += "*"
            text(
                x + cellw / 2,
                y + cellh / 2 + 2.7,
                label,
                8.5,
                True,
                "center",
                color="#FFFFFF" if a.log2_resistance_ratio > 5 else "#111111",
            )
        else:
            text(
                x + cellw / 2,
                y + cellh / 2 + 2.7,
                "Not reported",
                7,
                align="center",
                color="#777777",
            )
text(
    3,
    326,
    "* Parental IC50 interval has a non-estimable bound; fold-change precision is unresolved.",
    7,
)
text(
    3,
    338,
    "All ratios are descriptive calculations from published fits; no new tests or pooled estimates.",
    7,
)
barx, bary, bw, bh = 220, 354, 210, 8
for j in range(210):
    v = -2 + j / 209 * 10
    c.setFillColor(color(v))
    c.rect(barx + j, H - bary - bh, 1.1, bh, fill=1, stroke=0)
for t in [-2, 0, 2, 4, 6, 8]:
    x = barx + (t + 2) / 10 * bw
    line(x, bary + bh, x, bary + bh + 3)
    text(x, bary + bh + 13, str(t), 7, align="center")
text(3, 361, "Cell colour: log2 resistance ratio", 7.5)
text(
    3,
    395,
    "Rampa et al., Clin Cancer Res 2026;32:1454–61. Supplementary Tables 1 and 3 (CC BY 4.0).",
    7,
    color="#555555",
)
c.save()
d = pdfium.PdfDocument(str(out))
p = d[0]
p.render(scale=300 / 72).to_pil().save(out.with_suffix(".png"), dpi=(300, 300))
p.render(scale=1.75).to_pil().save(ROOT / "Supplementary_Figure_2_preview.png")
p.close()
d.close()
# Tab-separated companion source measurements, preserving the full IC50 audit.
raw.to_csv(ROOT / "immunoblot_deposited_ratios.tsv", sep="\t", index=False)
ic.to_csv(ROOT / "ic50_descriptive_ratios.tsv", sep="\t", index=False)
