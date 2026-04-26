"""Matplotlib-based pipeline visualization tool.

Generates a single publication-quality multi-panel figure from the
pipeline summary dict, plus individual panel PNGs, and saves them
into ``output/visualizations/``.

Style: white backgrounds, serif type, Tol-bright colorblind-safe palette.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import matplotlib.gridspec as gridspec
import numpy as np


# ── Scientific palette (Tol bright – colorblind-safe) ────────────────
_BLUE     = "#4477AA"
_CYAN     = "#66CCEE"
_GREEN    = "#228833"
_YELLOW   = "#CCBB44"
_RED      = "#EE6677"
_PURPLE   = "#AA3377"
_GREY     = "#BBBBBB"
_BLACK    = "#222222"
_DKGREY   = "#666666"
_LTGREY   = "#E5E5E5"
_PALETTE  = [_BLUE, _RED, _GREEN, _PURPLE, _CYAN, _YELLOW]

# Evo2 log-prob quality tiers (mean log-prob per nucleotide)
_EVO2_GOOD  = "#228833"  # > -0.3
_EVO2_OK    = "#CCBB44"  # -0.3 to -0.5
_EVO2_POOR  = "#EE6677"  # < -0.5

_DPI = 300

# ── Readable labels for check names ──────────────────────────────────
_LABEL_MAP = {
    "gc_content":             "GC Content",
    "cai_score":              "CAI Score",
    "restriction_sites":      "Restriction Sites",
    "rna_secondary_structure": "RNA MFE (kcal/mol)",
    "back_translation":       "Back-translation",
    "rare_codons":            "Rare Codons",
    "direct_repeats":         "Direct Repeats",
}


def _apply_theme() -> None:
    """Set global rcParams for clean scientific figures."""
    plt.rcParams.update({
        "font.family":        "serif",
        "font.size":          10,
        "axes.titlesize":     11,
        "axes.labelsize":     10,
        "xtick.labelsize":    9,
        "ytick.labelsize":    9,
        "legend.fontsize":    8,
        "text.color":         _BLACK,
        "axes.labelcolor":    _BLACK,
        "axes.edgecolor":     _DKGREY,
        "xtick.color":        _DKGREY,
        "ytick.color":        _DKGREY,
        "figure.facecolor":   "white",
        "axes.facecolor":     "white",
        "savefig.facecolor":  "white",
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.linewidth":     0.8,
        "axes.grid":          False,
        "grid.color":         _LTGREY,
        "grid.alpha":         0.6,
        "legend.frameon":     True,
        "legend.edgecolor":   _LTGREY,
        "legend.facecolor":   "white",
        "legend.framealpha":  1.0,
        "patch.linewidth":    0.5,
    })


_apply_theme()


# ── Helpers ──────────────────────────────────────────────────────────


# ── Panel C: Expression cassette architecture ────────────────────────

def _panel_cassette(ax: plt.Axes, summary: dict) -> None:
    """Horizontal stacked bar showing expression cassette elements."""
    chains = summary.get("chains", [])
    plasmids = summary.get("plasmids", [])
    if not chains or not plasmids:
        ax.set_visible(False)
        return

    # Calculate cassette element sizes from known constants
    # ATG (3 bp) + 6xHis (18 bp) + EK site (15 bp) + gene + stop (6 bp)
    START_LEN = 3   # ATG
    TAG_LEN = 18    # 6xHis = CAC×6
    EK_LEN = 15     # DDDDK = 5 codons
    STOP_LEN = 6    # TAATAA

    ids = []
    elements_list = []
    for pm in plasmids:
        chain_id = pm.get("chain_id", "?")
        insert_size = pm.get("insert_size", 0)
        gene_len = max(0, insert_size - START_LEN - TAG_LEN - EK_LEN - STOP_LEN)

        ids.append(chain_id)
        elements_list.append([
            ("ATG", START_LEN),
            ("6×His", TAG_LEN),
            ("EK site", EK_LEN),
            ("Gene", gene_len),
            ("Stop", STOP_LEN),
        ])

    if not ids:
        ax.set_visible(False)
        return

    y = np.arange(len(ids))
    element_colors = {
        "ATG": _RED,
        "6×His": _BLUE,
        "EK site": _CYAN,
        "Gene": _GREEN,
        "Stop": _PURPLE,
    }

    for i, elements in enumerate(elements_list):
        left = 0
        for name, size in elements:
            ax.barh(i, size, left=left, height=0.5,
                    color=element_colors.get(name, _GREY),
                    edgecolor="white", linewidth=0.8, zorder=3)
            # Label inside bar if wide enough
            if size > 20:
                ax.text(left + size / 2, i, f"{name}\n{size} bp",
                        ha="center", va="center", fontsize=6.5,
                        color="white", fontweight="bold")
            left += size

        # Total annotation
        total = sum(s for _, s in elements)
        ax.text(left + 3, i, f"{total} bp",
                va="center", fontsize=8, fontweight="bold", color=_BLACK)

    ax.set_yticks(y)
    ax.set_yticklabels(ids, fontfamily="monospace", fontsize=8)
    ax.set_xlabel("Cassette length (bp)")
    
    # Invert and squash y-axis (min height to prevent large stretching for single chains)
    min_rows_to_render = 2
    ax.set_ylim(max(len(ids), min_rows_to_render) - 0.5, -0.5)
    
    ax.set_title("(a) Expression Cassette Architecture", loc="left",
                 fontweight="bold", pad=8)

    # Legend
    legend_patches = [
        mpatches.Patch(color=element_colors["ATG"], label="Start codon"),
        mpatches.Patch(color=element_colors["6×His"], label="6×His tag"),
        mpatches.Patch(color=element_colors["EK site"], label="Protease site"),
        mpatches.Patch(color=element_colors["Gene"], label="Gene"),
        mpatches.Patch(color=element_colors["Stop"], label="Double stop"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=6, ncol=2)


# ── Panel D: Plasmid construct (donut) ───────────────────────────────

def _panel_plasmid(ax: plt.Axes, summary: dict) -> None:
    """Donut chart – first plasmid construct."""
    plasmids = summary.get("plasmids", [])
    if not plasmids:
        ax.set_visible(False)
        return

    pm = plasmids[0]
    insert = pm.get("insert_size", 0)
    backbone = pm.get("backbone_size", 0)
    total = pm.get("total_size", insert + backbone)
    vector = pm.get("vector", "pET-28a(+)")
    sites = pm.get("cloning_sites", [])
    chain_id = pm.get("chain_id", "")

    if total == 0:
        ax.set_visible(False)
        return

    wedges, _ = ax.pie(
        [backbone, insert],
        colors=[_BLUE, _GREEN],
        startangle=90, counterclock=False,
        wedgeprops=dict(width=0.30, edgecolor="white", linewidth=1.2),
    )

    # Centre text
    ax.text(0, 0.05, f"{total:,}", fontsize=16, fontweight="bold",
            ha="center", va="center", color=_BLACK, fontfamily="serif")
    ax.text(0, -0.12, "bp", fontsize=9, ha="center", va="center",
            color=_DKGREY)

    ins_pct = insert / total * 100
    bb_pct = backbone / total * 100
    ring_labels = [
        f"Backbone  {backbone:,} bp ({bb_pct:.1f}%)",
        f"Insert  {insert:,} bp ({ins_pct:.1f}%)",
    ]
    ax.legend(wedges, ring_labels, loc="lower center", ncol=1,
              bbox_to_anchor=(0.5, -0.12), fontsize=7)

    title = f"(b) Plasmid: {vector}"
    if chain_id:
        title += f" :: {chain_id}"
    ax.set_title(title, loc="left", fontweight="bold", pad=10, fontsize=10)

    if sites:
        ax.text(0.5, -0.02, f"Cloning sites: {' / '.join(sites)}",
                ha="center", fontsize=7.5, color=_DKGREY, style="italic",
                transform=ax.transAxes)


# ── Panel E: Pipeline overview gauges ────────────────────────────────

def _panel_overview(ax: plt.Axes, summary: dict) -> None:
    """Normalised 0-100% gauge bars for key metrics.

    All metrics are rescaled to a common 0-100 percentage axis so they
    are visually comparable. The original values are shown in annotations.
    """
    chains = summary.get("chains", [])
    if not chains:
        ax.set_visible(False)
        return

    # (label, raw_value, raw_range, raw_threshold, color)
    raw: list[tuple[str, float, tuple, float | None, str]] = []

    ch = chains[0]
    for ck in ch.get("checks", []):
        name = ck["name"]
        val = ck.get("value")
        if name == "gc_content" and isinstance(val, (int, float)):
            raw.append(("GC Content", val, (0.3, 0.7), None, _BLUE))
        elif name == "cai_score" and isinstance(val, (int, float)):
            raw.append(("CAI", val, (0.0, 1.0), 0.8, _GREEN))

    sol = ch.get("solubility_score")
    if sol is not None:
        raw.append(("Solubility", sol, (0.0, 1.0), 0.45, _CYAN))

    # Evo2 mean log-prob (replaces old pLDDT)
    struct = summary.get("structure_report", [])
    if struct:
        rec = struct[0]
        mlp = rec.get("mean_log_prob")
        if mlp is not None:
            # Scale: log-prob range roughly -1.0 to 0.0, threshold at -0.3
            raw.append(("Evo2 (log-prob/nt)", mlp, (-1.0, 0.0), -0.3, _PURPLE))

    if not raw:
        ax.set_visible(False)
        return

    # Normalise everything to 0-100%
    labels: list[str] = []
    pct_vals: list[float] = []
    pct_thresholds: list[float | None] = []
    raw_strs: list[str] = []
    colors: list[str] = []

    for label, val, (lo, hi), thr, color in raw:
        span = hi - lo if hi != lo else 1.0
        pct = (val - lo) / span * 100.0
        pct = max(0.0, min(100.0, pct))
        pct_thr = ((thr - lo) / span * 100.0) if thr is not None else None

        labels.append(label)
        pct_vals.append(pct)
        pct_thresholds.append(pct_thr)
        raw_strs.append(f"{val:.4f}" if abs(val) < 1 else f"{val:.3f}" if abs(val) < 10 else f"{val:.1f}")
        colors.append(color)

    y = np.arange(len(labels))

    # Background bar (full 0-100)
    for i in range(len(labels)):
        ax.barh(i, 100, left=0, height=0.55, color=_LTGREY,
                edgecolor="none", zorder=1)

    # Value bars
    for i, (pct, color) in enumerate(zip(pct_vals, colors)):
        ax.barh(i, pct, left=0, height=0.38, color=color,
                edgecolor="white", linewidth=0.5, zorder=3)

    # Threshold markers
    for i, pct_thr in enumerate(pct_thresholds):
        if pct_thr is not None:
            ax.plot(pct_thr, i, marker="|", color=_RED, markersize=14,
                    markeredgewidth=1.5, zorder=4)

    # Raw-value annotations
    for i, raw_s in enumerate(raw_strs):
        ax.text(102, i, raw_s, va="center", ha="left", fontsize=8,
                fontweight="bold", color=_BLACK)

    ax.set_xlim(0, 110)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontfamily="serif", fontsize=10)
    ax.invert_yaxis()
    ax.set_title("(c) Key Quality Metrics", loc="left", fontweight="bold", pad=12)

    handles = [
        mpatches.Patch(color=_LTGREY, label="Full range"),
        plt.Line2D([0], [0], marker="|", color=_RED, linestyle="None",
                   markersize=8, markeredgewidth=1.5, label="Threshold"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6)


# ── Combined figure ──────────────────────────────────────────────────

def _build_combined_figure(summary: dict, out: Path) -> Path:
    """Create a single multi-panel figure combining all visualisations.

    Layout (2 rows × 2 cols):
        Row 0:  (a) Cassette architecture  |  (b) Plasmid construct
        Row 1:  (c) Key metrics overview   |  (Empty)
    """
    _apply_theme()

    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.30)

    # (a) Cassette architecture
    ax_cass = fig.add_subplot(gs[0, 0])
    _panel_cassette(ax_cass, summary)

    # (b) Plasmid
    ax_plasm = fig.add_subplot(gs[0, 1])
    _panel_plasmid(ax_plasm, summary)

    # (c) Overview gauges
    ax_over = fig.add_subplot(gs[1, 0])
    _panel_overview(ax_over, summary)

    # Suptitle
    protein = summary.get("protein", {})
    name = protein.get("name", summary.get("compound_id", ""))
    uid = protein.get("uniprot_id", summary.get("compound_id", ""))
    host = summary.get("host_recommendation", {}).get("primary_host", "")

    suptitle = f"{name} ({uid})"
    if host:
        suptitle += f" — {host}"
    fig.suptitle(suptitle, fontsize=15, fontweight="bold", y=0.995)

    path = out / "pipeline_report.png"
    fig.savefig(path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


# ── Individual panel exports (supplementary) ─────────────────────────

def _save_individual(summary: dict, out: Path) -> list[Path]:
    """Save each panel as a standalone figure."""
    _apply_theme()
    paths: list[Path] = []

    # Cassette architecture
    plasmids = summary.get("plasmids", [])
    if plasmids:
        fig, ax = plt.subplots(figsize=(7, max(2.5, len(plasmids) * 0.8 + 1)))
        _panel_cassette(ax, summary)
        fig.tight_layout()
        p = out / "cassette_architecture.png"
        fig.savefig(p, dpi=_DPI)
        plt.close(fig)
        paths.append(p)

    # Plasmid
    for i, pm in enumerate(plasmids):
        fig, ax = plt.subplots(figsize=(5, 5))
        # Build a mini-summary with just this plasmid
        _panel_plasmid(ax, {"plasmids": [pm]})
        fig.tight_layout()
        safe = pm.get("chain_id", f"plasmid_{i}").replace(" ", "_").replace("/", "_")
        p = out / f"plasmid_{safe}.png"
        fig.savefig(p, dpi=_DPI)
        plt.close(fig)
        paths.append(p)

    # Overview gauges
    fig, ax = plt.subplots(figsize=(6, 3))
    _panel_overview(ax, summary)
    if ax.get_visible():
        fig.tight_layout()
        p = out / "quality_metrics.png"
        fig.savefig(p, dpi=_DPI)
        plt.close(fig)
        paths.append(p)
    else:
        plt.close(fig)

    return paths


# ── Public API ───────────────────────────────────────────────────────

def generate_visualizations(summary: dict, output_dir: str | Path) -> list[Path]:
    """Generate all matplotlib visualisations.

    Returns list of written file paths. The primary output is a single
    combined multi-panel figure (``pipeline_report.png``), plus
    individual panel PNGs for selective inclusion.
    """
    _apply_theme()

    out = Path(output_dir) / "visualizations"
    out.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []

    # Primary: single combined figure
    combined = _build_combined_figure(summary, out)
    paths.append(combined)

    # Secondary: individual panels
    paths.extend(_save_individual(summary, out))

    return paths
