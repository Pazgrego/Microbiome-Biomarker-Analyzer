"""
Microbiome OTU Analysis Pipeline
=================================
Reads a CSV (or .numbers) OTU table, calculates relative abundances,
extracts 5 key biomarkers with biological thresholds, computes alpha diversity,
generates a polished interactive Plotly donut chart, and writes a standalone HTML report
with interactive dynamic status tooltips (info icons).

Usage:
    python microbiome_pipeline.py --input data.csv --output report.html
"""

import argparse
from html import escape as html_escape
import os
import sys
import warnings
from datetime import datetime

import pandas as pd
import numpy as np
import plotly.express as px  # Advanced interactive visualization engine
import plotly.io as pio

warnings.filterwarnings("ignore")

# ── Configuration ──────────────────────────────────────────────────────────────

BIOMARKERS = [
    "Akkermansia muciniphila",
    "Ruminococcus bromii",
    "Eubacterium sp",
    "Faecalibacterium prausnitzii",
    "Roseburia sp",
]

BIOMARKER_SEARCH = {
    "Akkermansia muciniphila": ["akkermansia"],
    "Ruminococcus bromii":     ["ruminococcus bromii", "ruminococcus 1", "ruminococcus 2",
                                "ruminococcus sp"],
    "Eubacterium sp":          ["eubacterium"],
    "Faecalibacterium prausnitzii": ["faecalibacterium", "prausnitzii"],
    "Roseburia sp":            ["roseburia"],
}

BIOMARKER_THRESHOLDS = {
    "Faecalibacterium prausnitzii": 0.010,  # 1.0%
    "Eubacterium sp":              0.010,  # 1.0%
    "Akkermansia muciniphila": 0.003,      # 0.3%
    "Roseburia sp":            0.005,      # 0.5%
    "Ruminococcus bromii":     0.005       # 0.5%
}

DEFAULT_BIOMARKER_THRESHOLD = 0.001        # 0.1%

# ── Helper functions ───────────────────────────────────────────────────────────

def load_data(file_path: str, abundance_column: str = None) -> pd.DataFrame:
    """Loads CSV or Apple .numbers file into a pandas DataFrame."""
    if not os.path.exists(file_path):
        print(f"Error: file not found at '{file_path}'", file=sys.stderr)
        sys.exit(1)

    _, ext = os.path.splitext(file_path.lower())

    if ext == ".csv":
        try:
            df = pd.read_csv(file_path)
        except Exception as e:
            print(f"Error reading CSV '{file_path}': {e}", file=sys.stderr)
            sys.exit(1)

    elif ext == ".numbers":
        try:
            from numbers_parser import Document
        except ImportError:
            print("Error: 'numbers-parser' package is required to read .numbers files.", file=sys.stderr)
            print("Install it with: pip install numbers-parser", file=sys.stderr)
            sys.exit(1)

        try:
            doc = Document(file_path)
            sheets = doc.sheets
            if not sheets or not sheets[0].tables:
                print(f"Error: no tables found in .numbers file '{file_path}'", file=sys.stderr)
                sys.exit(1)
            table = sheets[0].tables[0]
            data = table.rows(as_list=True)
            if not data:
                print(f"Error: table is empty in '{file_path}'", file=sys.stderr)
                sys.exit(1)

            headers = [str(h) if h is not None else f"Col{i}" for i, h in enumerate(data[0])]
            rows = data[1:]
            df = pd.DataFrame(rows, columns=headers)
        except Exception as e:
            print(f"Error parsing .numbers file '{file_path}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Error: Unsupported file extension '{ext}'. Must be .csv or .numbers", file=sys.stderr)
        sys.exit(1)

    df.columns = [c.strip() for c in df.columns]

    required_base = {"Name", "Taxonomy"}
    missing_base = required_base - set(df.columns)
    if missing_base:
        print(f"Error: Missing required columns {missing_base}.", file=sys.stderr)
        print(f"Found columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    all_cols = list(df.columns)
    abundance_candidates = [c for c in all_cols if "abundance" in c.lower() or c == "Abundance"]

    if abundance_column:
        if abundance_column not in df.columns:
            print(f"Error: Specified abundance column '{abundance_column}' not found in the dataset.", file=sys.stderr)
            print(f"Available columns: {all_cols}", file=sys.stderr)
            sys.exit(1)
        chosen_abundance = abundance_column
    else:
        if "Combined Abundance" in df.columns:
            chosen_abundance = "Combined Abundance"
        elif "Abundance" in df.columns:
            chosen_abundance = "Abundance"
        elif len(abundance_candidates) == 1:
            chosen_abundance = abundance_candidates[0]
        elif len(abundance_candidates) > 1:
            chosen_abundance = abundance_candidates[0]
            print(f"[load_data] Multiple abundance columns detected. Defaulting to '{chosen_abundance}'.")
        else:
            print("Error: Could not identify an abundance column.", file=sys.stderr)
            print(f"Columns found: {all_cols}", file=sys.stderr)
            sys.exit(1)

    print(f"[load_data] Using '{chosen_abundance}' as the primary analytical abundance metric.")
    df = df.rename(columns={chosen_abundance: "Abundance"})
    df["Abundance"] = pd.to_numeric(df["Abundance"], errors="coerce").fillna(0).astype(int)
    return df[["Name", "Taxonomy", "Abundance"]]


def _match_rows(df: pd.DataFrame, terms: list) -> pd.DataFrame:
    """Returns rows where ALL terms are found in the 'Taxonomy' string (case-insensitive)."""
    mask = pd.Series(True, index=df.index)
    tax_lower = df["Taxonomy"].astype(str).str.lower()
    for t in terms:
        mask = mask & tax_lower.str.contains(t.lower(), regex=False)
    return df[mask]


def _pretty_label(row: pd.Series) -> str:
    """Extracts a shorter, human-readable name from raw taxonomy."""
    tax = str(row["Taxonomy"])
    parts = [p.strip() for p in tax.split(";")]
    if len(parts) >= 7 and parts[6] and parts[6].lower() != "unclassified":
        return parts[6]
    if len(parts) >= 6 and parts[5] and parts[5].lower() != "unclassified":
        return f"{parts[5]} sp."
    return str(row["Name"])


def analyze_microbiome(df: pd.DataFrame) -> dict:
    """
    Executes the core analytics pipeline, including computational alpha diversity,
    biomarker verification, and high-performance HTML graph rendering using Plotly.
    """
    try:
        total_reads = df["Abundance"].sum()
        
        if total_reads <= 0:
            raise ValueError("Total Abundance sum is 0 or negative. Cannot compute relative abundances or diversity indices.")

        df["RelativeAbundance"] = df["Abundance"] / total_reads
        df_enriched = df[df["Abundance"] > 0].copy()
        print(f"[load_data] Loaded {len(df_enriched)} active OTU rows for analysis.")

        # ── 1. Computational Biology: Alpha Diversity Indices ──
        p_i = df_enriched["RelativeAbundance"].to_numpy()
        shannon_index = float(-np.sum(p_i * np.log(p_i)))
        simpson_index = float(1.0 - np.sum(p_i ** 2))

        # ── 2. Biomarker extraction with custom biological thresholds ──
        biomarker_results = []
        for display_name, search_terms in BIOMARKER_SEARCH.items():
            matched = _match_rows(df_enriched, search_terms)
            agg_abundance = int(matched["Abundance"].sum())
            agg_pct = float(matched["RelativeAbundance"].sum().round(4))
            
            current_threshold = BIOMARKER_THRESHOLDS.get(display_name, DEFAULT_BIOMARKER_THRESHOLD)
            status = "Sufficient" if agg_pct >= current_threshold else "Low / Absent"
            
            biomarker_results.append(
                {
                    "name": display_name,
                    "abundance": agg_abundance,
                    "relative_pct": agg_pct,
                    "status": status,
                    "threshold_pct": current_threshold,
                    "match_count": len(matched),
                }
            )

        # ── 3. Top 10 taxonomic profiling ──
        top10 = df_enriched.nlargest(10, "Abundance").copy()
        if not top10.empty:
            top10["Label"] = top10.apply(_pretty_label, axis=1)
            top10["Percentage"] = (top10["RelativeAbundance"] * 100).round(2)
        else:
            top10["Label"] = []
            top10["Percentage"] = []

        # ── 4. Interactive Data Visualization: Polished Plotly Donut Chart ──
        plotly_html_snippet = ""
        if not top10.empty:
            fig = px.pie(
                top10, 
                values='RelativeAbundance', 
                names='Label', 
                hole=0.5,
                color_discrete_sequence=px.colors.qualitative.Pastel,
                custom_data=['Percentage']
            )
            
            # FIXED: Removed 'insidetextanchor' and used 'texttemplate' to cleanly align inner slice texts
            fig.update_traces(
                textposition='inside',
                texttemplate='%{percent:.1%}',
                insidetextorientation='radial',
                hovertemplate="<b>Taxon:</b> %{label}<br><b>Proportion:</b> %{customdata[0]}%<extra></extra>"
            )
            
            fig.update_layout(
                height=450,
                margin=dict(t=20, b=20, l=20, r=20),
                showlegend=True,
                legend=dict(
                    orientation="h", 
                    yanchor="top", 
                    y=-0.1, 
                    xanchor="center", 
                    x=0.5,
                    font=dict(size=11)
                ),
                annotations=[
                    dict(
                        text='<b>Top 10<br>Species</b>',
                        x=0.5,
                        y=0.5,
                        font_size=15,
                        font_color='#1e1b4b',
                        showarrow=False,
                        align='center'
                    )
                ],
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)'
            )
            
            plotly_html_snippet = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')

        return {
            "total_reads": int(total_reads),
            "distinct_otus": int(len(df_enriched)),
            "shannon_index": round(shannon_index, 3),
            "simpson_index": round(simpson_index, 3),
            "biomarkers": biomarker_results,
            "top10": top10.to_dict(orient="records"),
            "plotly_html": plotly_html_snippet,
        }

    except Exception as e:
        print(f"\n[CRITICAL ERROR] Pipeline processing failed during data computation: {e}", file=sys.stderr)
        sys.exit(1)


def generate_html_report(analysis: dict, output_path: str):
    """Generates a highly styled, standalone HTML report with advanced layout structures."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Format Health Biomarkers rows with dynamic tooltip feature
    biomarker_rows_html = ""
    for b in analysis["biomarkers"]:
        status_class = "status-sufficient" if b["status"] == "Sufficient" else "status-low"
        pct_display = f"{b['relative_pct'] * 100:.2f}%"
        threshold_display = f"{b['threshold_pct'] * 100:.2f}%"
        
        biomarker_rows_html += f"""
        <tr>
          <td><strong>{html_escape(b['name'])}</strong></td>
          <td style="text-align:right; font-weight:600;">{pct_display}</td>
          <td>
            <div style="display: flex; align-items: center; gap: 8px;">
              <span class="status-badge {status_class}">{b['status']}</span>
              <div class="tooltip-container">
                <span class="info-icon">ⓘ</span>
                <span class="tooltip-text">Threshold: {threshold_display}</span>
              </div>
            </div>
          </td>
        </tr>"""

    # Format Top 10 profiling rows
    top10_rows_html = ""
    for idx, r in enumerate(analysis["top10"]):
        pct_display = f"{r['RelativeAbundance'] * 100:.2f}%"
        short_tax = str(r.get('Taxonomy', ''))
        if len(short_tax) > 75:
            short_tax = short_tax[:72] + "..."
        
        top10_rows_html += f"""
        <tr>
          <td>
            <div style="font-weight:600">{html_escape(r.get('Label', r['Name']))}</div>
            <div style="font-size:11px; color:#666;" title="{html_escape(str(r.get('Taxonomy','')))}">{html_escape(short_tax)}</div>
          </td>
          <td style="text-align:right; font-weight:600; color:#4f46e5; font-size:14px;">{pct_display}</td>
        </tr>"""

    chart_content = analysis["plotly_html"] if analysis["plotly_html"] else "<p>No visualization data available</p>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Microbiome Quality Report</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background-color: #f9fafb;
    color: #111827;
    margin: 0;
    padding: 40px 20px;
  }}
  header {{
    max-width: 1000px;
    margin: 0 auto 30px auto;
    border-bottom: 1px solid #e5e7eb;
    padding-bottom: 20px;
  }}
  h1 {{ margin: 0 0 5px 0; font-size: 28px; color: #1e1b4b; }}
  .subtitle {{ color: #4b5563; font-size: 14px; margin: 0; }}
  main {{
    max-width: 1000px;
    margin: 0 auto;
    display: grid;
    grid-template-columns: 1fr;
    gap: 30px;
  }}
  @media(min-width: 768px) {{
    main {{ grid-template-columns: 1fr; }}
  }}
  section h2 {{
    font-size: 18px;
    color: #374151;
    margin: 0 0 12px 0;
    border-left: 4px solid #4f46e5;
    padding-left: 10px;
  }}
  .card {{
    background: white;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    padding: 20px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.02);
    box-sizing: border-box;
  }}
  .metrics-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 15px;
  }}
  .metric-box {{
    background: #f3f4f6;
    padding: 15px;
    border-radius: 6px;
    text-align: center;
  }}
  .metric-val {{ font-size: 24px; font-weight: bold; color: #4f46e5; margin-bottom: 4px; }}
  .metric-lbl {{ font-size: 11px; color: #6b7280; text-transform: uppercase; letter-spacing: 0.5px; }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 13.5px;
  }}
  th {{
    background: #f9fafb;
    color: #4b5563;
    text-align: left;
    font-weight: 600;
    padding: 10px 12px;
    border-bottom: 1px solid #e5e7eb;
  }}
  td {{
    padding: 12px;
    border-bottom: 1px solid #f3f4f6;
    vertical-align: middle;
  }}
  tr:last-child td {{ border-bottom: none; }}
  .status-badge {{
    display: inline-block;
    padding: 3px 8px;
    border-radius: 12px;
    font-size: 11px;
    font-weight: 600;
    text-transform: uppercase;
  }}
  .status-sufficient {{ background: #dcfce7; color: #166534; }}
  .status-low {{ background: #fee2e2; color: #991b1b; }}
  .chart-wrap {{ 
    padding: 20px;
  }}
  
  /* ── CSS Tooltip Styling ── */
  .tooltip-container {{
    position: relative;
    display: inline-block;
    cursor: pointer;
  }}
  .info-icon {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 16px;
    height: 16px;
    border-radius: 50%;
    background-color: #9ca3af;
    color: white;
    font-size: 11px;
    font-weight: bold;
    text-align: center;
  }}
  .tooltip-container:hover .info-icon {{
    background-color: #4f46e5;
  }}
  .tooltip-text {{
    visibility: hidden;
    width: 110px;
    background-color: #1f2937;
    color: #fff;
    text-align: center;
    border-radius: 6px;
    padding: 5px 8px;
    font-size: 11px;
    position: absolute;
    z-index: 10;
    bottom: 125%;
    left: 50%;
    margin-left: -55px;
    opacity: 0;
    transition: opacity 0.2s;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  }}
  .tooltip-text::after {{
    content: "";
    position: absolute;
    top: 100%;
    left: 50%;
    margin-left: -5px;
    border-width: 5px;
    border-style: solid;
    border-color: #1f2937 transparent transparent transparent;
  }}
  .tooltip-container:hover .tooltip-text {{
    visibility: visible;
    opacity: 1;
  }}
  
  footer {{
    max-width: 1000px;
    margin: 40px auto 0 auto;
    text-align: center;
    font-size: 12px;
    color: #9ca3af;
    border-top: 1px solid #e5e7eb;
    padding-top: 20px;
  }}
</style>
</head>
<body>

<header>
  <h1>🔬 Microbiome Analysis Report</h1>
  <p class="subtitle">Automated high-throughput OTU profiling pipeline</p>
</header>

<main>

  <section>
    <h2>📈 Sequencing Overview & Alpha Diversity</h2>
    <div class="card metrics-grid">
      <div class="metric-box">
        <div class="metric-val">{analysis['total_reads']}</div>
        <div class="metric-lbl">Total Sample Depth</div>
      </div>
      <div class="metric-box">
        <div class="metric-val">{analysis['distinct_otus']}</div>
        <div class="metric-lbl">Observed Richness (OTUs)</div>
      </div>
      <div class="metric-box">
        <div class="metric-val">{analysis['shannon_index']}</div>
        <div class="metric-lbl">Shannon Index (H')</div>
      </div>
      <div class="metric-box">
        <div class="metric-val">{analysis['simpson_index']}</div>
        <div class="metric-lbl">Simpson's Diversity (1-D)</div>
      </div>
    </div>
  </section>

  <section>
    <h2>🛡️ Health Biomarkers</h2>
    <div class="card" style="padding:0; overflow:hidden;">
      <table>
        <thead>
          <tr>
            <th>Target Taxon</th>
            <th style="text-align:right; width:30%;">Relative Abundance</th>
            <th style="width:25%;">Status</th>
          </tr>
        </thead>
        <tbody>{biomarker_rows_html}</tbody>
      </table>
    </div>
  </section>

  <section>
    <h2>🎨 Composition Profile (Interactive)</h2>
    <div class="card chart-wrap">
      <div style="width: 100%;">
        {chart_content}
      </div>
    </div>
  </section>

  <section>
    <h2>📊 Top 10 Most Abundant Species</h2>
    <div class="card" style="padding:0; overflow:hidden;">
      <table>
        <thead>
          <tr>
            <th>Species / Taxon</th>
            <th style="text-align:right; width:30%;">Relative Abundance</th>
          </tr>
        </thead>
        <tbody>{top10_rows_html}</tbody>
      </table>
    </div>
  </section>

</main>

<footer>
  <p>Microbiome Pipeline · generated {ts} · for research purposes only</p>
</footer>

</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"[generate_html_report] Report saved → {output_path}")


# ── CLI entry point ────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Microbiome OTU analysis pipeline")
    parser.add_argument("--input",  "-i", required=True, help="Path to OTU table (.csv or .numbers)")
    parser.add_argument("--output", "-o", default="report.html", help="Output HTML file (default: report.html)")
    parser.add_argument(
        "--abundance", "-a",
        default=None,
        help="Abundance column to use (default: Combined Abundance, or sole abundance column)",
    )
    args = parser.parse_args()

    df       = load_data(args.input, abundance_column=args.abundance)
    analysis = analyze_microbiome(df)
    generate_html_report(
        analysis,
        output_path=args.output
    )


if __name__ == "__main__":
    main()