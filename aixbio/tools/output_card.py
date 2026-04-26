from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from pathlib import Path


def generate_output_card(summary: dict, output_dir: str | Path) -> Path:
    """Generate a self-contained HTML visualization card from a pipeline summary."""
    out = Path(output_dir)
    compound_id = summary.get("compound_id", "unknown")
    protein = summary.get("protein", {})
    protein_name = protein.get("name", compound_id)
    status = summary.get("status", "unknown")
    host_rec = summary.get("host_recommendation", {})
    chains = summary.get("chains", [])
    warnings = summary.get("warnings", [])
    synthesis = summary.get("synthesis_quotes", [])

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    status_color = {"completed": "#16a34a", "failed": "#dc2626"}.get(status, "#a16207")
    status_label = status.upper()

    chains_html = ""
    for ch in chains:
        chain_id = html.escape(ch.get("chain_id", ""))
        ch_status = ch.get("status", "unknown")
        ch_color = {"passed": "#16a34a", "failed": "#dc2626", "max_retries_exceeded": "#dc2626"}.get(ch_status, "#a16207")
        insert_size = ch.get("insert_size", "?")
        remediation = ch.get("remediation_rounds", 0)
        sol_score = ch.get("solubility_score")
        disulfide = ch.get("disulfide_risk", False)

        sol_html = ""
        if sol_score is not None:
            sol_pct = round(sol_score * 100)
            sol_label = "Soluble" if sol_score >= 0.45 else "Inclusion body risk"
            sol_bar_color = "#16a34a" if sol_score >= 0.45 else "#f59e0b"
            sol_html = f"""
            <div class="metric">
              <span class="metric-label">Solubility</span>
              <div class="bar-track"><div class="bar-fill" style="width:{sol_pct}%;background:{sol_bar_color}"></div></div>
              <span class="metric-value">{sol_score:.2f} &mdash; {sol_label}</span>
            </div>"""

        disulfide_html = ""
        if disulfide:
            disulfide_html = '<span class="tag tag-warn">Disulfide risk &mdash; refolding required</span>'

        checks_rows = ""
        for ck in ch.get("checks", []):
            ck_name = html.escape(str(ck.get("name", "")))
            ck_passed = ck.get("passed", False)
            ck_value = html.escape(str(ck.get("value", "")))
            ck_thresh = html.escape(str(ck.get("threshold", "")))
            icon = '<span class="check-pass">&#10003;</span>' if ck_passed else '<span class="check-fail">&#10007;</span>'
            checks_rows += f"<tr><td>{icon}</td><td>{ck_name}</td><td>{ck_value}</td><td class='thresh'>{ck_thresh}</td></tr>\n"

        chains_html += f"""
        <div class="chain-card">
          <div class="chain-header">
            <span class="chain-id">{chain_id}</span>
            <span class="status-badge" style="background:{ch_color}">{ch_status.upper()}</span>
          </div>
          <div class="chain-stats">
            <div class="stat"><span class="stat-val">{insert_size}</span><span class="stat-label">bp insert</span></div>
            <div class="stat"><span class="stat-val">{remediation}</span><span class="stat-label">remediation rounds</span></div>
          </div>
          {sol_html}
          {disulfide_html}
          <table class="checks-table">
            <thead><tr><th></th><th>Check</th><th>Value</th><th>Threshold</th></tr></thead>
            <tbody>{checks_rows}</tbody>
          </table>
        </div>"""

    host_html = ""
    if host_rec:
        h_name = html.escape(host_rec.get("primary_host", ""))
        h_conf = html.escape(host_rec.get("confidence", ""))
        h_reason = html.escape(host_rec.get("reasoning", ""))
        alts = ", ".join(html.escape(a) for a in host_rec.get("alternative_hosts", []))
        features = host_rec.get("features", {})

        caveats_html = ""
        for c in host_rec.get("caveats", []):
            caveats_html += f'<li>{html.escape(c)}</li>'

        host_html = f"""
        <div class="section">
          <h2>Host Recommendation</h2>
          <div class="host-grid">
            <div class="host-main">
              <span class="host-name">{h_name}</span>
              <span class="status-badge" style="background:#2563eb">{h_conf.upper()} CONFIDENCE</span>
            </div>
            <p class="host-reasoning">{h_reason}</p>
            {"<p class='host-alt'>Alternatives: " + alts + "</p>" if alts else ""}
            <div class="feature-chips">
              <span class="chip">{features.get("total_chain_length", "?")} aa</span>
              <span class="chip">{features.get("total_cysteine_count", "?")} Cys</span>
              <span class="chip">{features.get("n_glycosylation_sites", "?")} N-glyc sites</span>
              <span class="chip">GRAVY {features.get("max_gravy", "?")}</span>
            </div>
            {"<ul class='caveats'>" + caveats_html + "</ul>" if caveats_html else ""}
          </div>
        </div>"""

    synthesis_html = ""
    if synthesis:
        vendor_rows = ""
        for sq in synthesis:
            sq_chain = html.escape(sq.get("chain_id", ""))
            sq_len = sq.get("insert_length", "?")
            sq_gc = sq.get("gc_content", 0)
            for v in sq.get("vendors", []):
                v_name = html.escape(v.get("vendor", ""))
                v_ok = v.get("feasible", False)
                v_color = "#16a34a" if v_ok else "#dc2626"
                v_label = "FEASIBLE" if v_ok else "FLAGGED"
                v_cost = f"~${v['estimated_cost_usd']:.0f}" if v.get("estimated_cost_usd") else "N/A"
                v_eta = html.escape(v.get("estimated_turnaround") or "N/A")
                notes = "; ".join(html.escape(n) for n in v.get("notes", []))
                vendor_rows += f"""<tr>
                  <td>{sq_chain}</td>
                  <td>{v_name}</td>
                  <td><span class="status-badge-sm" style="background:{v_color}">{v_label}</span></td>
                  <td>{v_cost}</td>
                  <td>{v_eta}</td>
                  <td class="notes">{notes}</td>
                </tr>"""

        synthesis_html = f"""
        <div class="section">
          <h2>Synthesis Feasibility</h2>
          <table class="synth-table">
            <thead><tr><th>Chain</th><th>Vendor</th><th>Status</th><th>Est. Cost</th><th>Turnaround</th><th>Notes</th></tr></thead>
            <tbody>{vendor_rows}</tbody>
          </table>
        </div>"""

    warnings_html = ""
    if warnings:
        items = "".join(f"<li>{html.escape(w)}</li>" for w in warnings)
        warnings_html = f"""
        <div class="section warnings-section">
          <h2>Warnings</h2>
          <ul class="warnings-list">{items}</ul>
        </div>"""

    files_html = ""
    artifact_files = sorted(out.glob("*"))
    if artifact_files:
        file_items = ""
        for f in artifact_files:
            if f.name == "output_card.html":
                continue
            icon = _file_icon(f.suffix)
            file_items += f'<div class="file-item"><span class="file-icon">{icon}</span><span class="file-name">{html.escape(f.name)}</span><span class="file-size">{_human_size(f.stat().st_size)}</span></div>'
        files_html = f"""
        <div class="section">
          <h2>Output Artifacts</h2>
          <div class="file-grid">{file_items}</div>
        </div>"""

    card_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AIXBio &mdash; {html.escape(protein_name)} ({html.escape(compound_id)})</title>
<style>
  :root {{
    --bg: #0f172a; --surface: #1e293b; --surface2: #334155;
    --text: #f1f5f9; --text2: #94a3b8; --accent: #38bdf8;
    --green: #16a34a; --red: #dc2626; --amber: #f59e0b;
    --radius: 12px; --font: 'Inter', system-ui, -apple-system, sans-serif;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: var(--font); padding: 24px; min-height: 100vh; }}
  .card {{ max-width: 900px; margin: 0 auto; }}
  .header {{ background: linear-gradient(135deg, #1e3a5f 0%, #0f172a 100%); border: 1px solid var(--surface2); border-radius: var(--radius); padding: 32px; margin-bottom: 20px; position: relative; overflow: hidden; }}
  .header::before {{ content: ''; position: absolute; top: -40px; right: -40px; width: 180px; height: 180px; background: radial-gradient(circle, rgba(56,189,248,0.12) 0%, transparent 70%); }}
  .header-top {{ display: flex; justify-content: space-between; align-items: flex-start; flex-wrap: wrap; gap: 12px; }}
  .brand {{ font-size: 13px; text-transform: uppercase; letter-spacing: 2px; color: var(--accent); margin-bottom: 8px; }}
  .protein-name {{ font-size: 28px; font-weight: 700; margin-bottom: 4px; }}
  .uniprot {{ font-size: 14px; color: var(--text2); font-family: monospace; }}
  .status-pill {{ display: inline-block; padding: 6px 16px; border-radius: 20px; font-size: 13px; font-weight: 600; color: white; letter-spacing: 0.5px; }}
  .timestamp {{ font-size: 12px; color: var(--text2); margin-top: 12px; }}

  .section {{ background: var(--surface); border: 1px solid var(--surface2); border-radius: var(--radius); padding: 24px; margin-bottom: 16px; }}
  h2 {{ font-size: 16px; font-weight: 600; margin-bottom: 16px; color: var(--accent); }}

  .chain-card {{ background: var(--surface); border: 1px solid var(--surface2); border-radius: var(--radius); padding: 24px; margin-bottom: 16px; }}
  .chain-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px; }}
  .chain-id {{ font-size: 18px; font-weight: 600; font-family: monospace; }}
  .status-badge {{ display: inline-block; padding: 3px 10px; border-radius: 6px; font-size: 11px; font-weight: 700; color: white; letter-spacing: 0.5px; }}
  .status-badge-sm {{ display: inline-block; padding: 2px 8px; border-radius: 4px; font-size: 10px; font-weight: 700; color: white; letter-spacing: 0.5px; }}

  .chain-stats {{ display: flex; gap: 32px; margin-bottom: 16px; }}
  .stat {{ display: flex; flex-direction: column; }}
  .stat-val {{ font-size: 24px; font-weight: 700; color: var(--text); }}
  .stat-label {{ font-size: 12px; color: var(--text2); }}

  .metric {{ margin-bottom: 12px; }}
  .metric-label {{ font-size: 12px; color: var(--text2); display: block; margin-bottom: 4px; }}
  .bar-track {{ height: 8px; background: var(--surface2); border-radius: 4px; overflow: hidden; margin-bottom: 4px; }}
  .bar-fill {{ height: 100%; border-radius: 4px; transition: width 0.4s; }}
  .metric-value {{ font-size: 13px; color: var(--text2); }}

  .tag {{ display: inline-block; padding: 4px 10px; border-radius: 6px; font-size: 12px; margin-bottom: 12px; }}
  .tag-warn {{ background: rgba(245,158,11,0.15); color: var(--amber); border: 1px solid rgba(245,158,11,0.3); }}

  .checks-table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
  .checks-table th {{ text-align: left; color: var(--text2); font-weight: 500; padding: 8px 12px; border-bottom: 1px solid var(--surface2); }}
  .checks-table td {{ padding: 8px 12px; border-bottom: 1px solid rgba(51,65,85,0.5); }}
  .checks-table .thresh {{ color: var(--text2); }}
  .check-pass {{ color: var(--green); font-weight: bold; }}
  .check-fail {{ color: var(--red); font-weight: bold; }}

  .host-grid {{ display: flex; flex-direction: column; gap: 8px; }}
  .host-main {{ display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }}
  .host-name {{ font-size: 20px; font-weight: 600; }}
  .host-reasoning {{ font-size: 14px; color: var(--text2); line-height: 1.5; }}
  .host-alt {{ font-size: 13px; color: var(--text2); }}
  .feature-chips {{ display: flex; gap: 8px; flex-wrap: wrap; margin-top: 4px; }}
  .chip {{ background: var(--surface2); padding: 4px 10px; border-radius: 6px; font-size: 12px; font-family: monospace; }}
  .caveats {{ list-style: none; margin-top: 8px; }}
  .caveats li {{ font-size: 13px; color: var(--amber); padding: 4px 0; }}
  .caveats li::before {{ content: '!'; display: inline-block; width: 18px; height: 18px; line-height: 18px; text-align: center; background: rgba(245,158,11,0.15); border-radius: 4px; font-size: 11px; font-weight: 700; margin-right: 8px; }}

  .synth-table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
  .synth-table th {{ text-align: left; color: var(--text2); font-weight: 500; padding: 8px 12px; border-bottom: 1px solid var(--surface2); }}
  .synth-table td {{ padding: 8px 12px; border-bottom: 1px solid rgba(51,65,85,0.5); }}
  .synth-table .notes {{ color: var(--text2); font-size: 12px; }}

  .warnings-section {{ border-color: rgba(245,158,11,0.3); }}
  .warnings-list {{ list-style: none; }}
  .warnings-list li {{ font-size: 13px; color: var(--amber); padding: 6px 0; border-bottom: 1px solid rgba(245,158,11,0.1); }}
  .warnings-list li:last-child {{ border-bottom: none; }}
  .warnings-list li::before {{ content: '\\26A0'; margin-right: 8px; }}

  .file-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 8px; }}
  .file-item {{ display: flex; align-items: center; gap: 8px; background: var(--bg); border: 1px solid var(--surface2); border-radius: 8px; padding: 10px 12px; }}
  .file-icon {{ font-size: 18px; }}
  .file-name {{ font-size: 13px; font-family: monospace; flex: 1; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
  .file-size {{ font-size: 11px; color: var(--text2); white-space: nowrap; }}

  .footer {{ text-align: center; padding: 20px; font-size: 12px; color: var(--text2); }}

  @media (max-width: 640px) {{
    body {{ padding: 12px; }}
    .header {{ padding: 20px; }}
    .protein-name {{ font-size: 22px; }}
    .chain-stats {{ gap: 16px; }}
    .stat-val {{ font-size: 20px; }}
  }}
</style>
</head>
<body>
<div class="card">
  <div class="header">
    <div class="header-top">
      <div>
        <div class="brand">AIXBio Pipeline</div>
        <div class="protein-name">{html.escape(protein_name)}</div>
        <div class="uniprot">{html.escape(compound_id)}</div>
      </div>
      <span class="status-pill" style="background:{status_color}">{status_label}</span>
    </div>
    <div class="timestamp">{timestamp}</div>
  </div>

  {host_html}
  {chains_html}
  {synthesis_html}
  {warnings_html}
  {files_html}

  <div class="footer">
    Generated by AIXBio &middot; Digital-to-Biological Pipeline
  </div>
</div>
</body>
</html>"""

    card_path = out / "output_card.html"
    card_path.write_text(card_html)
    return card_path


def _file_icon(suffix: str) -> str:
    return {
        ".json": "&#x1F4CA;",
        ".fasta": "&#x1F9EC;",
        ".gb": "&#x1F4C4;",
        ".tsv": "&#x1F4CB;",
        ".md": "&#x1F4DD;",
        ".txt": "&#x1F4C3;",
        ".html": "&#x1F310;",
        ".cif": "&#x1F9CA;",
    }.get(suffix, "&#x1F4C1;")


def _human_size(size: int) -> str:
    for unit in ("B", "KB", "MB"):
        if size < 1024:
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"
