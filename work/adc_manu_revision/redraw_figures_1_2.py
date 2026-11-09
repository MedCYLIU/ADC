"""Render TCGA and AURORA panels from the derived analysis tables."""

from pathlib import Path
import hashlib
import json
import os

import numpy as np
import pandas as pd
import pypdfium2 as pdfium
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

BASE = Path(__file__).resolve().parent
DATA = BASE / "figure12_source"
OUT = BASE / "figures"
WIDTH = 180 / 25.4 * 72
INK = "#202020"
GRAY = "#7A7A7A"
BLUE = "#2F5D7C"
TEAL = "#4C8C82"
GENES = ("ERBB2", "TACSTD2")
PALETTE = (BLUE, TEAL)


def register_fonts():
    windows = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    pairs = [
        (windows / "arial.ttf", windows / "arialbd.ttf", "Arial"),
        (
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
            "DejaVu Sans",
        ),
    ]
    for regular, bold, family in pairs:
        if regular.is_file() and bold.is_file():
            pdfmetrics.registerFont(TTFont("PlotSans", str(regular)))
            pdfmetrics.registerFont(TTFont("PlotSansBold", str(bold)))
            return family
    raise FileNotFoundError(
        "Install Arial or DejaVu Sans; no silent font substitution."
    )


class Figure:
    def __init__(self, number, height, title):
        self.number, self.height = number, height
        self.path = OUT / f"Figure_{number}.pdf"
        self.c = canvas.Canvas(str(self.path), pagesize=(WIDTH, height))
        self.c.setTitle(title)  # PDF metadata only; no figure-level heading.
        self.c.setAuthor("Liu et al.")
        self.c.setSubject("Display-only revision of frozen public-data estimates")

    def text(self, x, y, text, size=9.5, bold=False, align="left", color=INK):
        self.c.setFont("PlotSansBold" if bold else "PlotSans", size)
        self.c.setFillColor(HexColor(color))
        method = {
            "left": "drawString",
            "center": "drawCentredString",
            "right": "drawRightString",
        }[align]
        getattr(self.c, method)(x, self.height - y, str(text))

    def line(self, x1, y1, x2, y2, width=0.65, color=GRAY, dash=None):
        self.c.setLineWidth(width)
        self.c.setStrokeColor(HexColor(color))
        self.c.setDash(dash or [])
        self.c.line(x1, self.height - y1, x2, self.height - y2)
        self.c.setDash([])

    def point(self, x, y, radius=2.65, filled=True, color=BLUE, alpha=1):
        self.c.saveState()
        self.c.setFillColor(HexColor(color) if filled else white)
        self.c.setStrokeColor(HexColor(color))
        self.c.setLineWidth(0.85)
        self.c.setFillAlpha(alpha)
        self.c.setStrokeAlpha(alpha)
        self.c.circle(x, self.height - y, radius, stroke=1, fill=1)
        self.c.restoreState()

    def square(self, x, y, color, side=5):
        self.c.setFillColor(HexColor(color))
        self.c.rect(
            x - side / 2, self.height - y - side / 2, side, side, stroke=0, fill=1
        )

    def vertical_text(self, x, y, text, size=10):
        self.c.saveState()
        self.c.translate(x, self.height - y)
        self.c.rotate(90)
        self.c.setFont("PlotSans", size)
        self.c.setFillColor(HexColor(INK))
        self.c.drawCentredString(0, 0, text)
        self.c.restoreState()

    def save(self):
        self.c.showPage()
        self.c.save()
        pdf = pdfium.PdfDocument(self.path)
        pdf[0].render(scale=600 / 72).to_pil().save(
            self.path.with_suffix(".png"), dpi=(600, 600)
        )
        pdf[0].render(scale=2).to_pil().save(OUT / f"Figure_{self.number}_preview.png")
        pdf.close()


def load_frozen():
    """Load the TCGA and AURORA tables produced by run_analysis.py."""
    data = BASE.parent / "adc_v2"
    paths = {
        "tcga_plot_data": data / "tcga/tcga_plot_data.csv",
        "tcga_summary": data / "tcga/tcga_summary.csv",
        "tcga_results": data / "tcga/tcga_results.json",
        "patient_paired": data / "met/patient_paired.csv",
        "aurora_results": data / "met/results.json",
    }
    manifest = {
        key: {
            "path": str(path.relative_to(BASE.parent.parent)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for key, path in paths.items()
    }
    td = pd.read_csv(paths["tcga_plot_data"])
    ts = pd.read_csv(paths["tcga_summary"])
    tr = json.loads(paths["tcga_results"].read_text(encoding="utf-8"))
    ad = pd.read_csv(paths["patient_paired"])
    ar = json.loads(paths["aurora_results"].read_text(encoding="utf-8"))
    assert len(td) == td.patient_id.nunique() == 620
    assert len(ts) == 8
    assert len(ad) == 78 and ad.patient.nunique() == 39
    assert not ad.duplicated(["patient", "gene"]).any()
    for gene in GENES:
        assert len(ad[ad.gene.eq(gene)]) == ar["results"][gene]["n_donors"] == 39
        same = (
            ad.loc[ad.gene.eq(gene), "same_preservation_n"].to_numpy()
            == ad.loc[ad.gene.eq(gene), "n_met"].to_numpy()
        )
        assert same.sum() == 13 and (~same).sum() == 26
        for score, n in enumerate([61, 270, 199, 90]):
            summary = ts[ts.gene.eq(gene) & ts.HER2_IHC_score.eq(score)]
            assert len(summary) == 1 and int(summary.iloc[0].n) == n
            assert td.her2_ihc_score.eq(score).sum() == n
    return td, ts, tr, ad, ar, manifest


def figure1(td, ts, tr):
    fig = Figure(1, 312, "TCGA-BRCA transcript abundance by recorded HER2 IHC score")
    top, bottom, side = 67, 253, 202
    y = lambda value: bottom - (value - 3) / 16 * (bottom - top)
    # Exact same display-only displacement, row order and seed as adc_v6.
    jitter = np.random.default_rng(31811).uniform(-12.5, 12.5, len(td))
    for panel, (gene, color) in enumerate(zip(GENES, PALETTE)):
        left = [45, 301][panel]
        right = left + side
        fig.text(left - 39, 16, "AB"[panel], size=12, bold=True)
        fig.text(left, 16, gene, size=10.5, bold=True)
        stat = tr["stats"][gene]
        lo, hi = stat["rho_ci"]
        fig.text(left, 37, f"Spearman ρ = {stat['rho']:.2f}", size=10)
        fig.text(left, 53, f"95% CI {lo:.2f} to {hi:.2f}", size=9.5)
        fig.line(left, top, left, bottom)
        fig.line(left, bottom, right, bottom)
        for tick in [4, 8, 12, 16]:
            fig.line(left - 3.5, y(tick), left, y(tick))
            fig.text(left - 6, y(tick) + 3.2, tick, align="right")
        for score in range(4):
            x = left + 25 + score * (side - 50) / 3
            for i in np.flatnonzero(td.her2_ihc_score.to_numpy() == score):
                fig.point(
                    x + jitter[i],
                    y(td.iloc[i][gene + "_log2_rsem_plus1"]),
                    radius=1.2,
                    alpha=0.32,
                    color=color,
                )
            row = ts[ts.gene.eq(gene) & ts.HER2_IHC_score.eq(score)].iloc[0]
            fig.line(x, y(row.q25), x, y(row.q75), width=1.65, color=INK)
            for q in [row.q25, row.q75]:
                fig.line(x - 4, y(q), x + 4, y(q), width=1, color=INK)
            fig.line(
                x - 6.5,
                y(row["median"]),
                x + 6.5,
                y(row["median"]),
                width=1.8,
                color=INK,
            )
            fig.line(x, bottom, x, bottom + 3.5)
            fig.text(x, 269, ["0", "1+", "2+", "3+"][score], align="center")
            fig.text(x, 286, f"n = {int(row.n)}", align="center")
        fig.text(
            (left + right) / 2, 307, "Recorded HER2 IHC score", size=10, align="center"
        )
    fig.vertical_text(
        12, (top + bottom) / 2, "Expression, log2(normalized RSEM + 1)", size=10
    )
    fig.save()
    return hashlib.sha256(jitter.tobytes()).hexdigest()


def figure2(ad, ar):
    fig = Figure(2, 526, "AURORA US paired primary and metastatic transcript abundance")
    top, side, lo, hi = 34, 195, 6, 20
    bottom = top + side
    for panel, (gene, color) in enumerate(zip(GENES, PALETTE)):
        left = [45, 301][panel]
        x = lambda value: left + (value - lo) / (hi - lo) * side
        y = lambda value: bottom - (value - lo) / (hi - lo) * side
        fig.text(left - 39, 16, "AB"[panel], size=12, bold=True)
        fig.text(left, 16, gene, size=10.5, bold=True)
        fig.line(left, top, left, bottom)
        fig.line(left, bottom, left + side, bottom)
        fig.line(x(lo), y(lo), x(hi), y(hi), dash=[3, 3], color=GRAY)
        for tick in range(6, 21, 2):
            fig.line(x(tick), bottom, x(tick), bottom + 3.5)
            fig.line(left - 3.5, y(tick), left, y(tick))
            fig.text(x(tick), bottom + 16, tick, align="center")
            fig.text(left - 6, y(tick) + 3.2, tick, align="right")
        for row in ad[ad.gene.eq(gene)].itertuples():
            assert lo <= row.primary <= hi and lo <= row.metastatic_median <= hi
            fig.point(
                x(row.primary),
                y(row.metastatic_median),
                color=color,
                filled=(row.same_preservation_n == row.n_met),
            )
        fig.text(
            left + side / 2,
            265,
            "Primary expression, log2(UQN + 1)",
            size=9.5,
            align="center",
        )
    fig.vertical_text(12, top + side / 2, "Metastatic median, log2(UQN + 1)", size=10)
    fig.point(13, 286, color=INK)
    fig.text(23, 289.2, "All paired samples share preservation", size=9.3)
    fig.point(273, 286, color=INK, filled=False)
    fig.text(283, 289.2, "At least one preservation difference", size=9.3)
    fig.text(6, 314, "C", size=12, bold=True)
    fig.text(29, 314, "Paired changes and sensitivity analyses", size=10.5, bold=True)
    table_left, table_right = 6, WIDTH - 7
    forest_left, forest_right = 232, 367
    low, high = -1, 1.3
    xp = lambda value: forest_left + (value - low) / (high - low) * (
        forest_right - forest_left
    )
    fig.text(6, 336, "Target", bold=True)
    fig.text(64, 336, "Analysis", bold=True)
    fig.text(217, 336, "n", bold=True, align="right")
    fig.text(table_right, 336, "Median change (95% CI)", bold=True, align="right")
    fig.line(table_left, 343, table_right, 343)
    fig.line(xp(0), 348, xp(0), 470, color=GRAY, dash=[2, 2])
    rows = [
        (g, label, key)
        for g in GENES
        for label, key in [
            ("All paired individuals", None),
            ("Pre-treatment primary", "pretreatment_primary"),
            ("Matched preservation", "same_preservation"),
        ]
    ]
    estimates = []
    for i, ((gene, label, key), row_y) in enumerate(
        zip(rows, [360, 379, 398, 424, 443, 462])
    ):
        result = ar["results"][gene] if key is None else ar["results"][gene][key]
        n = result["n_donors"] if key is None else result["n"]
        value, (lower, upper) = result["median_delta"], result["median_delta_ci"]
        assert low <= lower <= value <= upper <= high
        color = PALETTE[GENES.index(gene)]
        if i in [0, 3]:
            fig.text(6, row_y + 3.2, gene, bold=True)
        fig.text(64, row_y + 3.2, label)
        fig.text(217, row_y + 3.2, n, align="right")
        fig.line(xp(lower), row_y, xp(upper), row_y, width=1, color=color)
        for bound in [lower, upper]:
            fig.line(
                xp(bound), row_y - 2.5, xp(bound), row_y + 2.5, width=1, color=color
            )
        fig.square(xp(value), row_y, color=color)
        fig.text(
            table_right,
            row_y + 3.2,
            f"{value:.2f} ({lower:.2f} to {upper:.2f})",
            align="right",
        )
        estimates.append(
            {
                "gene": gene,
                "analysis": key or "all",
                "n": n,
                "estimate": value,
                "lower": lower,
                "upper": upper,
            }
        )
    fig.line(table_left, 475, table_right, 475)
    fig.line(forest_left, 482, forest_right, 482)
    for tick in [-1, 0, 1]:
        fig.line(xp(tick), 482, xp(tick), 486)
        fig.text(xp(tick), 499, tick, align="center")
    fig.text(
        (forest_left + forest_right) / 2,
        520,
        "Metastatic − primary, log2(UQN + 1)",
        size=9.5,
        align="center",
    )
    fig.save()
    return estimates


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    font = register_fonts()
    td, ts, tr, ad, ar, sources = load_frozen()
    jitter_hash = figure1(td, ts, tr)
    estimates = figure2(ad, ar)
    manifest = {
        "purpose": "Display changes only; source estimates and intervals unchanged",
        "width_mm": 180,
        "png_dpi": 600,
        "font_family": font,
        "font_sizes_pt": {
            "body": 9.5,
            "axis_title": 10,
            "panel_title": 10.5,
            "panel_letter": 12,
            "legend": 9.3,
        },
        "palette": {"ERBB2": BLUE, "TACSTD2": TEAL, "axes": GRAY, "summaries": INK},
        "source_files": sources,
        "Figure_1": {
            "height_pt": 312,
            "points_per_panel": 620,
            "jitter_seed": 31811,
            "jitter_sha256": jitter_hash,
            "same_displacement_in_both_panels": True,
            "stats": tr["stats"],
            "summaries": ts.to_dict("records"),
        },
        "Figure_2": {
            "height_pt": 526,
            "points_per_panel": 39,
            "filled": 13,
            "open": 26,
            "estimates": estimates,
            "jitter": "None; source coordinates retained",
        },
    }
    (OUT / "Figure_1_2_render_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "outputs": [
                    "Figure_1.pdf",
                    "Figure_1.png",
                    "Figure_2.pdf",
                    "Figure_2.png",
                ],
                "width_mm": 180,
                "dpi": 600,
                "input_checks_passed": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
