"""Vector forest display of every deposited target; no selection by estimates."""

from pathlib import Path
import csv
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
import pypdfium2 as pdfium

BASE = Path(__file__).resolve().parent
OUT = BASE / "figures"
OUT.mkdir(exist_ok=True)
for n, f in [("Arial", "arial.ttf"), ("Arial-Bold", "arialbd.ttf")]:
    pdfmetrics.registerFont(TTFont(n, "C:/Windows/Fonts/" + f))
rows = list(csv.DictReader((BASE / "all_targets/target_estimates.csv").open()))
lookup = {
    (r["subset"], r["target"]): r for r in rows if r["estimand"] == "cross_minus_same"
}
targets = sorted({r["target"] for r in rows})
assert len(targets) == 72
W = 180 / 25.4 * 72
H = 563
c = canvas.Canvas(str(OUT / "Supplementary_Figure_1.pdf"), pagesize=(W, H))
c.setTitle("All 72 deposited target transcripts: matched specimen contrast")


def text(x, y, t, size=7, bold=False):
    c.setFillColor(HexColor("#111111"))
    c.setFont("Arial-Bold" if bold else "Arial", size)
    c.drawString(x, H - y, str(t))


def line(x, y, x2, y2, color="#111111", width=0.5, dash=None):
    c.setStrokeColor(HexColor(color))
    c.setLineWidth(width)
    c.setDash(dash or [])
    c.line(x, H - y, x2, H - y2)
    c.setDash([])


text(0, 12, "Matched different-tissue minus same-tissue error", 10, True)
text(0, 28, "All 72 source-listed target transcripts; equal donor weighting", 8)
text(0, 43, "●  All eligible donors (n = 30)       ○  ≥50% tumor nuclei (n = 28)", 8)
for col in range(3):
    left = col * 171
    x = lambda v: left + 47 + (v + 0.2) / 0.34 * 113
    line(x(0), 65, x(0), 499, "#999999", 0.5, [2, 2])
    for i, gene in enumerate(targets[col * 24 : (col + 1) * 24]):
        y = 71 + i * 18
        text(left, y + 2, gene, 7.1, bold=gene in ["HER2", "TROP2"])
        for subset, offset, filled, color in [
            ("all", -2.5, True, "#111111"),
            ("tumor_nuclei_50", 3.5, False, "#888888"),
        ]:
            r = lookup[subset, gene]
            m, lo, hi = [float(r[k]) for k in ["mean", "low", "high"]]
            line(x(lo), y + offset, x(hi), y + offset, color, 0.65)
            c.setStrokeColor(HexColor(color))
            c.setFillColor(HexColor(color if filled else "#FFFFFF"))
            c.circle(x(m), H - y - offset, 1.5, stroke=1, fill=1)
    line(x(-0.2), 507, x(0.14), 507)
    for v in [-0.2, -0.1, 0, 0.1]:
        line(x(v), 507, x(v), 510)
        text(x(v) - 5, 521, f"{v:g}".replace("-", "−"), 7)
    text(left + 45, 536, "Difference in MAE (VST)", 7.5)
text(
    0,
    557,
    "Nominal 95% donor-bootstrap intervals; not simultaneous. Negative values favor different-tissue pairs.",
    7.2,
)
c.save()
p = pdfium.PdfDocument(OUT / "Supplementary_Figure_1.pdf")
p[0].render(scale=300 / 72).to_pil().save(
    OUT / "Supplementary_Figure_1.png", dpi=(300, 300)
)
p[0].render(scale=1.5).to_pil().save(OUT / "Supplementary_Figure_1_preview.png")
print("Saved all-target supplement figure.")
