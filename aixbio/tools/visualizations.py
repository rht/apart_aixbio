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
import matplotlib.ticker as mticker
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

# AlphaFold canonical pLDDT colour bands
_AF_VHIGH = "#0053D6"
_AF_CONF  = "#65CBF3"
_AF_LOW   = "#FFDB13"
_AF_VLOW  = "#FF7D45"

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
        "grid.linewidth":     0.4,
        "grid.alpha":         0.6,
        "legend.frameon":     True,
        "legend.edgecolor":   _LTGREY,
        "legend.facecolor":   "white",
        "legend.framealpha":  1.0,
        "patch.linewidth":    0.5,
    })


_apply_theme()


# ── Helpers ──────────────────────────────────────────────────────────

def _normalise_plddt(v: float) -> float:
    """Ensure pLDDT is on a 0-100 scale (some pipelines return 0-1)."""
    return v * 100 if v <= 1.0 else v


def _plddt_color(v: float) -> str:
    """Return AlphaFold confidence-band colour for a pLDDT value (0-100)."""
    if v >= 90:
        return _AF_VHIGH
    if v >= 70:
        return _AF_CONF
    if v >= 50:
        return _AF_LOW
    return _AF_VLOW


def _extract_numeric_checks(chains: list[dict]) -> list[tuple[str, str, float, bool, str]]:
    """Return (chain_id, check_name, value, passed, threshold) for numeric checks."""
    rows = []
    for ch in chains:
        cid = ch.get("chain_id", "?")
        for ck in ch.get("checks", []):
            v = ck.get("value")
            if isinstance(v, (int, float)):
                rows.append((cid, ck["name"], float(v), ck["passed"],
                             ck.get("threshold", "")))
    return rows


# ── Panel A: Sequence quality metrics (normalised 0-1 radar) ────────

def _panel_validation(ax: plt.Axes, summary: dict) -> None:
    """Horizontal bar chart with paired value + threshold annotations."""
    chains = summary.get("chains", [])
    rows = _extract_numeric_checks(chains)
    if not rows:
        ax.set_visible(False)
        return

    labels = [_LABEL_MAP.get(r[1], r[1]) for r in rows]
    values = [r[2] for r in rows]
    passed = [r[3] for r in rows]
    thresholds = [r[4] for r in rows]
    colors = [_GREEN if p else _RED for p in passed]

    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=colors, height=0.5,
                   edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Measured value")
    ax.axvline(x=0, color=_DKGREY, linewidth=0.5, zorder=1)
    ax.set_title("(a) Sequence Validation Metrics", loc="left",
                 fontweight="bold", pad=8)

    # Annotations: value + threshold placed to the right of each bar
    max_abs = max(abs(v) for v in values) if values else 1
    pad = max_abs * 0.04
    for bar, val, thr in zip(bars, values, thresholds):
        # For positive values: right of bar.  For negative: right of zero line.
        if val >= 0:
            x_txt = bar.get_width() + pad
        else:
            x_txt = pad  # just right of the zero line
        ax.text(x_txt, bar.get_y() + bar.get_height() / 2,
                f"{val:.2f}" if abs(val) < 10 else f"{val:.0f}",
                va="center", ha="left", fontsize=8, color=_BLACK,
                fontweight="bold")
        if thr:
            ax.text(x_txt, bar.get_y() + bar.get_height() / 2 + 0.22,
                    f"thr: {thr}", va="center", ha="left",
                    fontsize=6.5, color=_DKGREY, style="italic")

    handles = [mpatches.Patch(facecolor=_GREEN, edgecolor="white", label="Pass"),
               mpatches.Patch(facecolor=_RED, edgecolor="white", label="Fail")]
    ax.legend(handles=handles, loc="lower right", fontsize=7)


# ── Panel B: Structure confidence (pLDDT + RMSD) ────────────────────

def _panel_structure(ax_plddt: plt.Axes, ax_rmsd: plt.Axes | None,
                     summary: dict) -> None:
    """pLDDT bar with AlphaFold confidence bands, plus RMSD panel."""
    recs = summary.get("structure_report", [])
    if not recs:
        ax_plddt.set_visible(False)
        if ax_rmsd is not None:
            ax_rmsd.set_visible(False)
        return

    ids = [r.get("chain_id", "?") for r in recs]
    plddts_raw = [r.get("plddt_mean", 0) for r in recs]
    plddts = [_normalise_plddt(p) for p in plddts_raw]
    rmsds = [r.get("rmsd_to_ref") for r in recs]
    methods = [r.get("method", "?") for r in recs]

    y = np.arange(len(ids))

    # ── pLDDT ──
    colors = [_plddt_color(p) for p in plddts]
    bars = ax_plddt.barh(y, plddts, color=colors, height=0.45,
                         edgecolor="white", linewidth=0.5, zorder=3)
    ax_plddt.set_yticks(y)
    ax_plddt.set_yticklabels(ids, fontfamily="monospace", fontsize=8)
    ax_plddt.set_xlim(0, 100)
    ax_plddt.set_xlabel("pLDDT score")
    ax_plddt.invert_yaxis()

    # Shaded confidence bands
    for lo, hi, c in [(0, 50, _AF_VLOW), (50, 70, _AF_LOW),
                      (70, 90, _AF_CONF), (90, 100, _AF_VHIGH)]:
        ax_plddt.axvspan(lo, hi, color=c, alpha=0.07, zorder=0)

    for bar, plddt, method in zip(bars, plddts, methods):
        label = {"esmfold": "ESMFold", "afdb": "AFDB",
                 "afdb_fallback": "AFDB*"}.get(method, method)
        ax_plddt.text(
            min(plddt + 2, 92), bar.get_y() + bar.get_height() / 2,
            f"{plddt:.1f}  ({label})", va="center", fontsize=8, color=_BLACK)

    ax_plddt.set_title("(b) Structure Confidence (pLDDT)", loc="left",
                       fontweight="bold", pad=8)

    band_patches = [
        mpatches.Patch(color=_AF_VHIGH, label="Very high (≥ 90)"),
        mpatches.Patch(color=_AF_CONF,  label="Confident (≥ 70)"),
        mpatches.Patch(color=_AF_LOW,   label="Low (≥ 50)"),
        mpatches.Patch(color=_AF_VLOW,  label="Very low (< 50)"),
    ]
    ax_plddt.legend(handles=band_patches, loc="lower right", fontsize=6,
                    title="AlphaFold confidence", title_fontsize=6)

    # ── RMSD ──
    if ax_rmsd is None:
        return
    has_rmsd = any(r is not None for r in rmsds)
    if not has_rmsd:
        ax_rmsd.set_visible(False)
        return

    rmsd_vals = [r if r is not None else 0 for r in rmsds]
    rmsd_colors = []
    for r in rmsds:
        if r is None:
            rmsd_colors.append(_GREY)
        elif r < 2:
            rmsd_colors.append(_GREEN)
        elif r < 5:
            rmsd_colors.append(_YELLOW)
        else:
            rmsd_colors.append(_RED)

    bars2 = ax_rmsd.barh(y, rmsd_vals, color=rmsd_colors, height=0.45,
                         edgecolor="white", linewidth=0.5, zorder=3)
    ax_rmsd.set_yticks(y)
    ax_rmsd.set_yticklabels(ids, fontfamily="monospace", fontsize=8)
    ax_rmsd.set_xlabel("RMSD (Å)")
    ax_rmsd.invert_yaxis()
    ax_rmsd.set_title("(c) RMSD to Reference", loc="left",
                      fontweight="bold", pad=8)

    # Quality thresholds
    ax_rmsd.axvline(x=2, color=_GREEN, linestyle=":", linewidth=0.7, alpha=0.8)
    ax_rmsd.axvline(x=5, color=_RED, linestyle=":", linewidth=0.7, alpha=0.8)

    for bar, val in zip(bars2, rmsd_vals):
        if val > 0:
            ax_rmsd.text(bar.get_width() + 0.1,
                         bar.get_y() + bar.get_height() / 2,
                         f"{val:.2f} Å", va="center", fontsize=8, color=_BLACK)

    rmsd_patches = [
        mpatches.Patch(color=_GREEN, label="< 2 Å (excellent)"),
        mpatches.Patch(color=_YELLOW, label="2–5 Å (acceptable)"),
        mpatches.Patch(color=_RED, label="> 5 Å (poor)"),
    ]
    ax_rmsd.legend(handles=rmsd_patches, loc="lower right", fontsize=6)


# ── Panel C: Plasmid construct (donut) ───────────────────────────────

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

    title = f"(d) Plasmid: {vector}"
    if chain_id:
        title += f" :: {chain_id}"
    ax.set_title(title, loc="left", fontweight="bold", pad=10, fontsize=10)

    if sites:
        ax.text(0.5, -0.02, f"Cloning sites: {' / '.join(sites)}",
                ha="center", fontsize=7.5, color=_DKGREY, style="italic",
                transform=ax.transAxes)


# ── Panel D: Synthesis vendor comparison ─────────────────────────────

def _panel_synthesis(ax: plt.Axes, summary: dict) -> None:
    """Grouped bar chart with cost + feasibility + turnaround."""
    synthesis = summary.get("synthesis_quotes", [])
    if not synthesis:
        ax.set_visible(False)
        return

    chains: list[str] = []
    vendor_data: dict[str, list[dict]] = {}

    for sq in synthesis:
        chain_id = sq.get("chain_id", "?")
        chains.append(chain_id)
        for v in sq.get("vendors", []):
            vname = v.get("vendor", "?")
            vendor_data.setdefault(vname, []).append(v)

    vendors = list(vendor_data.keys())
    n_chains = len(chains)
    n_vendors = len(vendors)
    x = np.arange(n_chains)
    width = 0.30

    for i, vendor in enumerate(vendors):
        recs = vendor_data[vendor]
        costs = [r.get("estimated_cost_usd", 0) or 0 for r in recs]
        feas = [r.get("feasible", False) for r in recs]
        turnaround = [r.get("estimated_turnaround", "") for r in recs]
        offset = (i - (n_vendors - 1) / 2) * width
        color = _PALETTE[i % len(_PALETTE)]

        bars = ax.bar(x + offset, costs, width * 0.85, label=vendor,
                      color=color, edgecolor="white", linewidth=0.5, zorder=3)

        for bar, cost, ok, ta in zip(bars, costs, feas, turnaround):
            # Cost + feasibility label
            status = "feasible" if ok else "infeasible"
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 1,
                    f"${cost:.0f}\n({status})",
                    ha="center", va="bottom", fontsize=7,
                    color=_GREEN if ok else _RED, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(chains, fontfamily="monospace", fontsize=8)
    ax.set_ylabel("Estimated cost (USD)")
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("${x:,.0f}"))
    ax.set_title("(e) Gene Synthesis Vendor Comparison", loc="left",
                 fontweight="bold", pad=8)
    ax.legend(fontsize=7)

    # Add headroom for annotations
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.25)


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

    struct = summary.get("structure_report", [])
    if struct:
        plddt = _normalise_plddt(struct[0].get("plddt_mean", 0))
        raw.append(("pLDDT", plddt, (0.0, 100.0), 70.0, _PURPLE))

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
        raw_strs.append(f"{val:.3f}" if val < 10 else f"{val:.1f}")
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

    ax.set_xlim(0, 115)  # give room for annotations
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Normalised score (%)")
    ax.invert_yaxis()
    ax.set_title("(f) Key Quality Metrics", loc="left",
                 fontweight="bold", pad=8)

    handles = [
        mpatches.Patch(color=_LTGREY, label="Full range"),
        plt.Line2D([0], [0], marker="|", color=_RED, linestyle="None",
                   markersize=8, markeredgewidth=1.5, label="Threshold"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6)


# ── Combined figure ──────────────────────────────────────────────────

def _build_combined_figure(summary: dict, out: Path) -> Path:
    """Create a single multi-panel figure combining all visualisations.

    Layout (3 rows × 2 cols):
        Row 0:  (a) Validation checks  |  (b) pLDDT confidence
        Row 1:  (c) RMSD to ref        |  (d) Plasmid construct
        Row 2:  (e) Synthesis costs     |  (f) Key metrics overview
    """
    _apply_theme()

    fig = plt.figure(figsize=(14, 16))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.38, wspace=0.30)

    # (a) Validation checks
    ax_val = fig.add_subplot(gs[0, 0])
    _panel_validation(ax_val, summary)

    # (b) pLDDT
    ax_plddt = fig.add_subplot(gs[0, 1])

    # (c) RMSD
    ax_rmsd = fig.add_subplot(gs[1, 0])
    _panel_structure(ax_plddt, ax_rmsd, summary)

    # (d) Plasmid
    ax_plasm = fig.add_subplot(gs[1, 1])
    _panel_plasmid(ax_plasm, summary)

    # (e) Synthesis
    ax_synth = fig.add_subplot(gs[2, 0])
    _panel_synthesis(ax_synth, summary)

    # (f) Overview gauges
    ax_over = fig.add_subplot(gs[2, 1])
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

    # Validation
    chains = summary.get("chains", [])
    rows = _extract_numeric_checks(chains)
    if rows:
        fig, ax = plt.subplots(figsize=(7, max(2.5, len(rows) * 0.6)))
        _panel_validation(ax, summary)
        fig.tight_layout()
        p = out / "validation_checks.png"
        fig.savefig(p, dpi=_DPI)
        plt.close(fig)
        paths.append(p)

    # Structure quality
    struct = summary.get("structure_report", [])
    if struct:
        has_rmsd = any(s.get("rmsd_to_ref") is not None for s in struct)
        ncols = 2 if has_rmsd else 1
        fig, axes = plt.subplots(1, ncols,
                                 figsize=(5 * ncols, max(3, len(struct) * 0.8 + 1)))
        if ncols == 1:
            axes = [axes]
        _panel_structure(axes[0], axes[1] if has_rmsd else None, summary)
        fig.tight_layout()
        p = out / "structure_quality.png"
        fig.savefig(p, dpi=_DPI)
        plt.close(fig)
        paths.append(p)

    # Plasmid
    plasmids = summary.get("plasmids", [])
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

    # Synthesis
    if summary.get("synthesis_quotes"):
        fig, ax = plt.subplots(figsize=(max(5, 2.2), 4.5))
        _panel_synthesis(ax, summary)
        fig.tight_layout()
        p = out / "synthesis_comparison.png"
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
