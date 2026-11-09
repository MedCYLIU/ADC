"""Plot-only typography/layout update of frozen Figures 3-5.

Run with the bundled Python runtime from any directory. No analysis script is
imported or executed; no model, confidence interval, bootstrap, sample or outcome
is changed. Quartiles/medians used as graphical annotations follow the original
plot code's numpy linear quantile definition. Every input is hashed before and
after rendering, and all plotted coordinates/estimates are recorded in an audit.
"""

from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pypdfium2 as pdfium
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor

HERE = Path(__file__).resolve().parent
WORK = HERE.parent
OUT = HERE / "figures"
OUT.mkdir(parents=True, exist_ok=True)
W = 180 / 25.4 * 72
INK, BLUE, TEAL, GRAY = "#202020", "#2F5D7C", "#4C8C82", "#7A7A7A"
LIGHT = "#D2D2D2"
FONT, BOLD = "ADC_Redesign_Arial", "ADC_Redesign_ArialBold"
pdfmetrics.registerFont(TTFont(FONT, "C:/Windows/Fonts/arial.ttf"))
pdfmetrics.registerFont(TTFont(BOLD, "C:/Windows/Fonts/arialbd.ttf"))


def rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Figure:
    def __init__(self, n, height):
        self.n, self.h = n, height
        self.path = OUT / f"Figure_{n}.pdf"
        self.c = canvas.Canvas(str(self.path), pagesize=(W, height), pageCompression=1)
        self.c.setTitle(f"Figure {n}")
        self.c.setSubject(
            "Plot-only revision of unchanged deposited-data analysis outputs"
        )
        self.text_boxes = []

    def text(self, x, y, s, size=10, bold=False, color=INK, align="left"):
        s = str(s)
        font = BOLD if bold else FONT
        width = pdfmetrics.stringWidth(s, font, size)
        left = (
            x if align == "left" else x - width if align == "right" else x - width / 2
        )
        assert left >= -0.1 and left + width <= W + 0.1, (self.n, s, left, width)
        assert 0 <= y <= self.h, (self.n, s, y)
        self.text_boxes.append(
            {
                "text": s,
                "size_pt": size,
                "x_left": left,
                "baseline_y": y,
                "width": width,
            }
        )
        self.c.setFillColor(HexColor(color))
        self.c.setFont(font, size)
        {
            "left": self.c.drawString,
            "right": self.c.drawRightString,
            "center": self.c.drawCentredString,
        }[align](x, self.h - y, s)

    def line(self, x, y, xx, yy, color=INK, width=0.65, dash=()):
        self.c.setStrokeColor(HexColor(color))
        self.c.setLineWidth(width)
        self.c.setDash(dash)
        self.c.line(x, self.h - y, xx, self.h - yy)
        self.c.setDash([])

    def marker(self, x, y, kind="circle", color=BLUE, radius=2.1, opened=False):
        yy = self.h - y
        c = self.c
        r = radius
        c.setStrokeColor(HexColor(color))
        c.setFillColor(HexColor("#FFFFFF" if opened else color))
        c.setLineWidth(0.65)
        if kind == "circle":
            c.circle(x, yy, r, fill=1, stroke=1)
        elif kind == "square":
            c.rect(x - r, yy - r, 2 * r, 2 * r, fill=1, stroke=1)
        elif kind == "triangle":
            p = c.beginPath()
            p.moveTo(x, yy + r * 1.18)
            p.lineTo(x - r * 1.08, yy - r * 0.82)
            p.lineTo(x + r * 1.08, yy - r * 0.82)
            p.close()
            c.drawPath(p, fill=1, stroke=1)
        elif kind == "cross":
            c.line(x - r, yy - r, x + r, yy + r)
            c.line(x - r, yy + r, x + r, yy - r)
        else:
            raise ValueError(kind)

    def ci(self, low, high, y, mid, color=BLUE):
        self.line(low, y, high, y, color, 1.0)
        self.line(low, y - 2.2, low, y + 2.2, color, 0.8)
        self.line(high, y - 2.2, high, y + 2.2, color, 0.8)
        self.marker(mid, y, "square", color, 2.4)

    def panel(self, letter, title, y, x=4):
        self.text(x, y, letter, 12, True)
        self.text(x + 20, y, title, 11, True)

    def axis(self, xp, lo, hi, y, ticks, label):
        self.line(xp(lo), y, xp(hi), y)
        for v in ticks:
            self.line(xp(v), y, xp(v), y + 3)
            self.text(xp(v), y + 15, f"{v:g}", align="center")
        self.text((xp(lo) + xp(hi)) / 2, y + 30, label, align="center")

    def rotated(self, x, y, s):
        self.c.saveState()
        self.c.translate(x, self.h - y)
        self.c.rotate(90)
        self.c.setFillColor(HexColor(INK))
        self.c.setFont(FONT, 10)
        self.c.drawCentredString(0, 0, s)
        self.c.restoreState()

    def finish(self):
        self.c.showPage()
        self.c.save()
        pdf = pdfium.PdfDocument(self.path)
        pdf[0].render(scale=600 / 72).to_pil().save(
            OUT / f"Figure_{self.n}.png", dpi=(600, 600)
        )
        pdf[0].render(scale=1.75).to_pil().save(
            OUT / f"Figure_{self.n}_preview.png", dpi=(126, 126)
        )
        pdf.close()
        return {
            "size_mm": [180, self.h / 72 * 25.4],
            "body_font_pt": 10,
            "png_dpi": 600,
            "text_bounds": self.text_boxes,
        }


def donors():
    data = rows(WORK / "adc/lesions.csv")
    ids = sorted({r["patient"] for r in data})
    counts = Counter(r["patient"] for r in data)
    assert len(data) == 749 and len(ids) == 30
    f = Figure(3, 200 / 25.4 * 72)
    styles = {
        "Liver": ("circle", BLUE, False),
        "Lymph nodes": ("circle", BLUE, True),
        "Lung": ("triangle", TEAL, False),
        "Bones": ("square", GRAY, True),
        "Other": ("cross", GRAY, False),
    }
    for x, site in zip([21, 107, 243, 327, 420], styles):
        kind, color, opened = styles[site]
        f.marker(x, 38, kind, color, 2.2, opened)
        f.text(x + 8, 41, site)
    f.text(4, 62, "Donor (n)", bold=True)
    top, spacing, bottom = 70, 14.95, 522
    low, high = 1, 11.2
    lefts, width = [76, 309], W - 309 - 7
    plotted = {}
    for panel, (target, label, left) in enumerate(
        zip(["HER2", "TROP2"], ["HER2 / ERBB2", "TROP2 / TACSTD2"], lefts)
    ):
        f.panel("AB"[panel], label, 16, left - 17)
        xp = lambda v: left + (v - low) / (high - low) * width
        assert all(low <= float(r[target]) <= high for r in data)
        f.line(left, top - 6, left, bottom)
        f.line(left, bottom, left + width, bottom)
        for tick in [2, 4, 6, 8, 10]:
            f.line(xp(tick), bottom, xp(tick), bottom + 3)
            f.text(xp(tick), bottom + 15, tick, align="center")
        f.text(left + width / 2, bottom + 31, "VST expression", align="center")
        plotted[target] = []
        for index, ident in enumerate(ids):
            y = top + (index + 0.5) * spacing
            group = [r for r in data if r["patient"] == ident]
            values = [float(r[target]) for r in group]
            if panel == 0:
                f.text(
                    left - 9, y + 3.2, f"{ident[2:]} ({counts[ident]})", align="right"
                )
            f.line(xp(min(values)), y, xp(max(values)), y, LIGHT, 0.6)
            for r in group:
                site = r["site"] if r["site"] in styles else "Other"
                kind, color, opened = styles[site]
                jitter = ((int(r["row"]) * 37) % 11 - 5) * 0.78
                f.marker(xp(float(r[target])), y + jitter, kind, color, 1.6, opened)
                plotted[target].append([r["patient"], r["row"], float(r[target])])
            med = xp(float(np.median(values)))
            f.line(med, y - 5.2, med, y + 5.2, INK, 1.1)
    audit = f.finish()
    audit.update(
        {
            "specimens": 749,
            "donors": 30,
            "donor_order": ids,
            "specimens_per_donor": dict(counts),
            "all_plotted_values": plotted,
            "tissue_styles": styles,
            "x_range": [1, 11.2],
            "caption_note": "Point shapes and colours identify tissue categories; pale segments show observed ranges and black ticks show medians. All 749 observations per target retained.",
        }
    )
    return audit


def interval_string(r):
    m, a, b = [float(r[k]) for k in ["mean", "low", "high"]]
    return f"{m:.3f} ({a:.3f} to {b:.3f})"


def matched():
    d = WORK / "adc_v2/audit"
    s = {(r["target"], r["estimand"]): r for r in rows(d / "matched_results.csv")}
    sens = {
        (r["target"], r["estimand"]): r for r in rows(d / "cellularity_results.csv")
    }
    donors_data = rows(d / "matched_donor_estimates.csv")
    f = Figure(4, 532)
    f.panel("A", "Absolute error of expression estimates", 15)
    f.text(20, 36, "Target / sampling", bold=True)
    f.text(180, 36, "Donors", bold=True, align="center")
    f.text(506, 36, "Mean (95% CI)", bold=True, align="right")
    f.line(20, 43, 506, 43, LIGHT)
    xp = lambda v: 205 + 150 * v
    for target, start, color in [("HER2", 58, BLUE), ("TROP2", 128, TEAL)]:
        for j, (key, label) in enumerate(
            [
                ("index", "Index specimen"),
                ("same_pair", "Same-tissue pair"),
                ("cross_pair", "Different-tissue pair"),
            ]
        ):
            y = start + 20 * j
            if j == 0:
                f.text(20, y + 3.3, target, bold=True)
            f.text(65, y + 3.3, label)
            f.text(180, y + 3.3, "30", align="center")
            r = s[target, key]
            m, a, b = [float(r[k]) for k in ["mean", "low", "high"]]
            f.ci(xp(a), xp(b), y, xp(m), color)
            f.text(506, y + 3.3, interval_string(r), align="right")
    f.axis(xp, 0, 1, 187, [0, 0.2, 0.4, 0.6, 0.8, 1], "Mean absolute error (VST units)")
    f.panel("B", "Different-tissue minus same-tissue error", 252)
    f.text(506, 275, "Mean (95% CI)", bold=True, align="right")
    xp = lambda v: 128 + (v + 0.36) / 0.88 * 227
    f.line(xp(0), 282, xp(0), 334, GRAY, 0.7, [2, 2])
    for target, y, color in [("HER2", 297, BLUE), ("TROP2", 326, TEAL)]:
        f.text(20, y + 3.3, target, bold=True)
        vals = sorted(
            float(r["cross_minus_same"]) for r in donors_data if r["target"] == target
        )
        assert len(vals) == 30
        for i, v in enumerate(vals):
            f.marker(xp(v), y - 4 + ((i * 5) % 7 - 3) * 1.05, "circle", GRAY, 1.8, True)
        r = s[target, "cross_minus_same"]
        m, a, b = [float(r[k]) for k in ["mean", "low", "high"]]
        f.ci(xp(a), xp(b), y + 5, xp(m), color)
        f.text(506, y + 3.3, interval_string(r), align="right")
    f.axis(
        xp,
        -0.36,
        0.52,
        345,
        [-0.3, -0.1, 0, 0.1, 0.3, 0.5],
        "Error difference (VST units)",
    )
    f.panel("C", "Sensitivity: at least 50% tumor nuclei", 408)
    f.text(20, 432, "Target", bold=True)
    f.text(180, 432, "Donors", bold=True, align="center")
    f.text(506, 432, "Mean difference (95% CI)", bold=True, align="right")
    f.line(20, 439, 506, 439, LIGHT)
    xp = lambda v: 205 + (v + 0.12) / 0.24 * 150
    f.line(xp(0), 443, xp(0), 483, GRAY, 0.7, [2, 2])
    for target, y, color in [("HER2", 455, BLUE), ("TROP2", 479, TEAL)]:
        f.text(20, y + 3.3, target, bold=True)
        f.text(180, y + 3.3, "28", align="center")
        r = sens[target, "cross_minus_same"]
        m, a, b = [float(r[k]) for k in ["mean", "low", "high"]]
        f.ci(xp(a), xp(b), y, xp(m), color)
        f.text(506, y + 3.3, interval_string(r), align="right")
    f.axis(
        xp,
        -0.12,
        0.12,
        496,
        [-0.1, 0, 0.1],
        "Different-tissue minus same-tissue error (VST)",
    )
    audit = f.finish()
    audit.update(
        {
            "matched_estimates": list(s.values()),
            "sensitivity_estimates": list(sens.values()),
            "donor_contrast_points": donors_data,
            "models_refitted": False,
            "caption_note": "Coloured squares and lines show mean estimates and 95% CIs; open grey circles in B show donor-level differences. Blue denotes HER2 and teal denotes TROP2.",
        }
    )
    return audit


def ispy():
    d = WORK / "adc_v2/main"
    data = rows(d / "ispy_joined.csv")
    models = rows(d / "ispy_models.csv")
    stats = json.loads((d / "ispy_stats.json").read_text(encoding="utf8"))
    assert len(data) == len({r["id"] for r in data}) == 127
    assert np.isclose(
        np.std([float(r["ERBB2"]) for r in data], ddof=1), stats["sd_ERBB2"]
    )
    f = Figure(5, 490)
    f.panel("A", "Pretreatment ERBB2 by pathological response", 15)
    top, bottom = 51, 258
    yp = lambda v: bottom - (v - 8.4) / (13 - 8.4) * (bottom - top)
    f.line(44, top, 44, bottom)
    f.line(44, bottom, 503, bottom)
    for tick in [9, 10, 11, 12, 13]:
        f.line(41, yp(tick), 44, yp(tick))
        f.text(35, yp(tick) + 3.3, tick, align="right")
    f.rotated(12, (top + bottom) / 2, "ERBB2 expression (log2)")
    count = 0
    plotted = []
    graphical_summaries = []
    for center, arm, label, model in zip(
        [121, 274, 427],
        ["T-DM1/P", "THP", "TH Control"],
        ["T-DM1/P", "THP", "TH"],
        models,
    ):
        selected = [r for r in data if r["arm"] == arm]
        assert len(selected) == int(model["n"])
        assert sum(int(r["pcr"]) for r in selected) == int(model["pcr"])
        f.text(center, 37, label, bold=True, align="center")
        for outcome, offset, color in [(0, -33, GRAY), (1, 33, BLUE)]:
            vals = sorted(
                float(r["ERBB2"]) for r in selected if int(r["pcr"]) == outcome
            )
            xc = center + offset
            used = []
            for v in vals:
                y = yp(v)
                candidates = [0] + [
                    sign * j * 4.6 for j in range(1, 10) for sign in (-1, 1)
                ]
                valid = [
                    shift
                    for shift in candidates
                    if all((shift - u) ** 2 + (y - z) ** 2 >= 4.6**2 for u, z in used)
                ]
                assert valid
                shift = valid[0]
                assert abs(shift) <= 23
                used.append((shift, y))
                f.marker(xc - 5 + shift, y, "circle", color, 2.0, outcome == 0)
                count += 1
                plotted.append({"arm": arm, "pcr": outcome, "ERBB2": v})
            q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75], method="linear")
            qx = xc + 22
            f.line(qx, yp(q1), qx, yp(q3), color, 1.3)
            f.line(qx - 3, yp(q1), qx + 3, yp(q1), color, 1.0)
            f.line(qx - 3, yp(q3), qx + 3, yp(q3), color, 1.0)
            f.line(qx - 4.5, yp(med), qx + 4.5, yp(med), color, 1.7)
            f.text(xc, 277, "Residual" if outcome == 0 else "pCR", align="center")
            f.text(xc, 292, f"n = {len(vals)}", align="center")
            graphical_summaries.append(
                {"arm": arm, "pcr": outcome, "q1": q1, "median": med, "q3": q3}
            )
    assert count == 127
    f.panel("B", "Association of ERBB2 with pCR", 330)
    f.text(20, 354, "Regimen", bold=True)
    f.text(112, 354, "pCR / n", bold=True, align="center")
    f.text(506, 354, "Adjusted OR (95% CI)", bold=True, align="right")
    f.line(20, 361, 506, 361, LIGHT)
    xp = (
        lambda v: 168
        + (math.log(v) - math.log(0.4)) / (math.log(32) - math.log(0.4)) * 188
    )
    f.line(xp(1), 368, xp(1), 438, GRAY, 0.7, [2, 2])
    for y, label, m, original in zip(
        [380, 405, 430], ["T-DM1/P", "THP", "TH"], models, stats["models"]
    ):
        assert all(float(m[k]) == original[k] for k in ["OR", "lo", "hi"])
        mid, lo, hi = [float(m[k]) for k in ["OR", "lo", "hi"]]
        f.text(20, y + 3.3, label)
        f.text(112, y + 3.3, f'{m["pcr"]} / {m["n"]}', align="center")
        f.ci(xp(lo), xp(hi), y, xp(mid), TEAL)
        f.text(506, y + 3.3, f"{mid:.2f} ({lo:.2f} to {hi:.2f})", align="right")
    f.axis(xp, 0.4, 32, 449, [0.5, 1, 2, 5, 10, 20], "Odds ratio per 1 SD higher ERBB2")
    audit = f.finish()
    audit.update(
        {
            "points": plotted,
            "graphical_quartiles": graphical_summaries,
            "frozen_models": models,
            "sd_ERBB2": stats["sd_ERBB2"],
            "models_refitted": False,
            "caption_note": "Open grey circles denote residual disease, filled blue circles denote pCR; adjacent lines show median and IQR. Teal squares and lines show the unchanged adjusted ORs and profile-likelihood CIs.",
        }
    )
    return audit


if __name__ == "__main__":
    inputs = [
        WORK / "adc/lesions.csv",
        WORK / "adc_v2/audit/matched_results.csv",
        WORK / "adc_v2/audit/cellularity_results.csv",
        WORK / "adc_v2/audit/matched_donor_estimates.csv",
        WORK / "adc_v2/main/ispy_joined.csv",
        WORK / "adc_v2/main/ispy_models.csv",
        WORK / "adc_v2/main/ispy_stats.json",
    ]
    before = {str(p): digest(p) for p in inputs}
    figs = {"3": donors(), "4": matched(), "5": ispy()}
    assert all(digest(k) == v for k, v in before.items())
    audit = {
        "operation": "plot-only typography, geometry and palette revision",
        "source_sha256": before,
        "inputs_unchanged_after_render": True,
        "inferential_analysis_executed": False,
        "models_refitted": False,
        "figures": figs,
        "embedding_width_mm": 165,
        "font_size_at_165_mm": 10 * 165 / 180,
        "software": {
            "numpy": np.__version__,
            "plot_backend": "ReportLab with embedded Arial",
            "rasterizer": "pypdfium2",
        },
    }
    (HERE / "figures_3_5_render_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf8"
    )
    print(
        "Figures 3-5: 180-mm vector PDFs and 600-dpi PNGs written; frozen inputs unchanged."
    )
