"""
Microbiome OTU Analysis Pipeline
=================================
Reads a CSV (or .numbers) OTU table, calculates relative abundances,
extracts 5 key biomarkers, generates a top-10 pie chart, and writes
a standalone HTML report.

Usage:
    python microbiome_pipeline.py --input data.csv --output report.html
    python microbiome_pipeline.py --input data.numbers --output report.html

CSV format expected:
    Name, Taxonomy, Abundance
    OTU_001, "Bacteria; Firmicutes; ...; Faecalibacterium prausnitzii", 5400
    ...
"""

import argparse
import base64
import io
import os
import sys
import warnings
from datetime import datetime

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

# ── Configuration ──────────────────────────────────────────────────────────────

BIOMARKERS = [
    "Akkermansia muciniphila",
    "Ruminococcus bromii",
    "Eubacterium sp",
    "Faecalibacterium prausnitzii",
    "Roseburia sp",
]

# Keyword(s) used for partial / case-insensitive taxonomy matching.
# Each entry maps display name → list of search terms (OR logic).
BIOMARKER_SEARCH = {
    "Akkermansia muciniphila": ["akkermansia"],
    "Ruminococcus bromii":     ["ruminococcus bromii", "ruminococcus 1", "ruminococcus 2",
                                 "ruminococcus sp"],
    # Only match rows where "Eubacterium" appears at the genus level (not as a
    # bracketed group prefix like "[Eubacterium] eligens group").
    # We require the Taxonomy field to contain "; Eubacterium;" or end with it.
    "Eubacterium sp":          ["eubacteriaceae; eubacterium"],
    "Faecalibacterium prausnitzii": ["faecalibacterium"],
    "Roseburia sp":            ["roseburia"],
}

STATUS_THRESHOLD = 0.1  # % – above this → "Sufficient" (green), else "Low" (red)

PIE_COLORS = [
    "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", "#264653",
    "#457B9D", "#A8DADC", "#F1FAEE", "#E63946", "#8ECAE6",
]


# ── 1. Data Loading ────────────────────────────────────────────────────────────

def load_data(filepath: str) -> pd.DataFrame:
    """
    Load OTU data from a CSV file or an Apple Numbers (.numbers) file.

    Returns a DataFrame with columns:
        Name       – OTU identifier or organism label
        Taxonomy   – full taxonomy string
        Abundance  – raw read count (numeric)
    """
    ext = os.path.splitext(filepath)[1].lower()

    if ext == ".numbers":
        df = _load_numbers(filepath)
    elif ext in (".csv", ".tsv", ".txt"):
        sep = "\t" if ext == ".tsv" else ","
        df = pd.read_csv(filepath, sep=sep)
        df = _normalise_columns(df)
    else:
        # Attempt CSV anyway
        df = pd.read_csv(filepath)
        df = _normalise_columns(df)

    # Coerce abundance to numeric; drop zero/missing rows
    df["Abundance"] = pd.to_numeric(df["Abundance"], errors="coerce").fillna(0)
    df = df[df["Abundance"] > 0].copy()
    df.reset_index(drop=True, inplace=True)

    print(f"[load_data] Loaded {len(df):,} OTU rows from '{filepath}'.")
    return df


def _load_numbers(filepath: str) -> pd.DataFrame:
    """Parse an Apple Numbers workbook and return the first table."""
    try:
        from numbers_parser import Document
    except ImportError:
        sys.exit(
            "numbers-parser is required for .numbers files.\n"
            "Install it with:  pip install numbers-parser"
        )
    doc = Document(filepath)
    table = doc.sheets[0].tables[0]
    rows = [
        [str(c.value) if c.value is not None else "" for c in row]
        for row in table.iter_rows()
    ]
    # Row 0 is a merged header; row 1 contains real column names
    header = rows[1]
    data = rows[2:]
    df = pd.DataFrame(data, columns=header)
    return _normalise_columns(df)


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Map varied column names to canonical: Name, Taxonomy, Abundance.
    Falls back gracefully when columns are missing.
    """
    col_map: dict[str, str] = {}
    for col in df.columns:
        low = col.lower().strip()
        if low in ("name", "otu", "otu id", "otu_id", "id"):
            col_map[col] = "Name"
        elif low in ("taxonomy", "taxon", "lineage", "microorganism assingee",
                     "microorganism assignee"):
            col_map[col] = "Taxonomy"
        elif "combined abundance" in low or low in ("abundance", "count",
                                                     "reads", "counts"):
            col_map[col] = "Abundance"

    df = df.rename(columns=col_map)

    # Ensure all three columns exist
    for required in ("Name", "Taxonomy", "Abundance"):
        if required not in df.columns:
            df[required] = ""

    return df[["Name", "Taxonomy", "Abundance"]]


# ── 2. Analysis ────────────────────────────────────────────────────────────────

def analyze_microbiome(df: pd.DataFrame) -> dict:
    """
    Perform relative-abundance calculations and extract biomarker data.

    Returns a dict:
        total_reads     – int
        df_enriched     – DataFrame with added RelativeAbundance (%) column
        top10           – DataFrame of the 10 most abundant species
        biomarkers      – list of dicts per biomarker with keys:
                          name, abundance, relative_pct, status
        summary_stats   – dict with high-level stats
    """
    total = df["Abundance"].sum()

    df_enriched = df.copy()
    df_enriched["RelativeAbundance"] = (df_enriched["Abundance"] / total * 100).round(4)

    # Top-10 by relative abundance
    top10 = (
        df_enriched
        .sort_values("RelativeAbundance", ascending=False)
        .head(10)
        .reset_index(drop=True)
    )
    # Pretty label for chart: last meaningful part of taxonomy or Name
    top10["Label"] = top10.apply(_pretty_label, axis=1)

    # Biomarker extraction
    biomarker_results = []
    for display_name, search_terms in BIOMARKER_SEARCH.items():
        matched = _match_rows(df_enriched, search_terms)
        agg_abundance = matched["Abundance"].sum()
        agg_pct = matched["RelativeAbundance"].sum().round(4)
        status = "Sufficient" if agg_pct >= STATUS_THRESHOLD else "Low / Absent"
        biomarker_results.append(
            {
                "name": display_name,
                "abundance": agg_abundance,
                "relative_pct": agg_pct,
                "status": status,
                "match_count": len(matched),
            }
        )

    summary_stats = {
        "total_reads": int(total),
        "unique_otus": len(df_enriched),
        "top_species_label": top10.iloc[0]["Label"] if len(top10) else "N/A",
        "top_species_pct": float(top10.iloc[0]["RelativeAbundance"]) if len(top10) else 0,
        "biomarkers_sufficient": sum(
            1 for b in biomarker_results if b["status"] == "Sufficient"
        ),
    }

    print(f"[analyze_microbiome] Total reads: {total:,.0f} | OTUs: {len(df_enriched):,}")
    for b in biomarker_results:
        print(
            f"  {b['name']:<35} {b['relative_pct']:6.3f}%  "
            f"({b['match_count']} OTUs)  → {b['status']}"
        )

    return {
        "total_reads": int(total),
        "df_enriched": df_enriched,
        "top10": top10,
        "biomarkers": biomarker_results,
        "summary_stats": summary_stats,
    }


def _match_rows(df: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    """Return rows whose Taxonomy (or Name) contains any of the search terms."""
    mask = pd.Series(False, index=df.index)
    for term in terms:
        mask |= df["Taxonomy"].str.contains(term, case=False, na=False)
        mask |= df["Name"].str.contains(term, case=False, na=False)
    return df[mask]


def _pretty_label(row: pd.Series) -> str:
    """Extract the last non-ambiguous part of a taxonomy string."""
    tax = row["Taxonomy"]
    parts = [p.strip() for p in tax.split(";") if p.strip()]
    # Walk backwards; skip generic tokens
    for part in reversed(parts):
        low = part.lower()
        if low not in ("ambiguous_taxa", "uncultured bacterium", "metagenome",
                        "unidentified", "uncultured", "ambiguous taxa"):
            return part
    return row["Name"]


# ── 3. Visualisation ──────────────────────────────────────────────────────────

def _build_pie_chart(top10: pd.DataFrame) -> str:
    """Render a matplotlib donut chart and return it as a base-64 PNG data-URI."""
    labels = top10["Label"].tolist()
    values = top10["RelativeAbundance"].tolist()
    colors = PIE_COLORS[: len(labels)]

    fig, ax = plt.subplots(figsize=(10, 6), facecolor="none")

    wedges, texts, autotexts = ax.pie(
        values,
        labels=None,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%" if pct > 2 else "",
        pctdistance=0.78,
        startangle=140,
        wedgeprops=dict(width=0.62, edgecolor="white", linewidth=1.5),
        explode=[0.04] + [0] * (len(values) - 1),
    )
    for at in autotexts:
        at.set_fontsize(8.5)
        at.set_color("#1a1a2e")
        at.set_fontweight("bold")

    # Donut centre label
    ax.text(0, 0, "Top 10\nSpecies", ha="center", va="center",
            fontsize=11, color="#1a1a2e", fontweight="bold", linespacing=1.4)

    # Legend on the right
    legend_labels = [
        f"{lbl[:38]}…  ({val:.2f}%)" if len(lbl) > 38 else f"{lbl}  ({val:.2f}%)"
        for lbl, val in zip(labels, values)
    ]
    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, legend_labels)]
    ax.legend(
        handles=patches,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=8.5,
        frameon=False,
        labelspacing=0.75,
    )

    ax.set_title("Top 10 Most Abundant Species", fontsize=14,
                 fontweight="bold", color="#1a1a2e", pad=16)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight",
                facecolor="none", transparent=True)
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    return f"data:image/png;base64,{b64}"


# ── 4. HTML Report ────────────────────────────────────────────────────────────

def generate_html_report(
    analysis: dict,
    output_path: str = "report.html",
    input_filename: str = "",
) -> None:
    """Write a self-contained HTML report to *output_path*."""
    top10 = analysis["top10"]
    biomarkers = analysis["biomarkers"]
    stats = analysis["summary_stats"]

    pie_src = _build_pie_chart(top10)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")

    # ── Biomarker table rows
    bm_rows_html = ""
    for bm in biomarkers:
        is_ok = bm["status"] == "Sufficient"
        badge_cls = "badge-ok" if is_ok else "badge-low"
        icon = "✔" if is_ok else "✘"
        bm_rows_html += f"""
        <tr>
          <td class="bm-name">{bm['name']}</td>
          <td class="num">{bm['abundance']:,.0f}</td>
          <td class="num">{bm['relative_pct']:.4f}%</td>
          <td class="num muted">{bm['match_count']}</td>
          <td><span class="badge {badge_cls}">{icon} {bm['status']}</span></td>
        </tr>"""

    # ── Top-10 table rows
    top10_rows_html = ""
    for _, row in top10.iterrows():
        bar_w = min(int(row["RelativeAbundance"] / top10["RelativeAbundance"].max() * 100), 100)
        top10_rows_html += f"""
        <tr>
          <td class="tax-label">{row['Label']}</td>
          <td class="num">{row['Abundance']:,.0f}</td>
          <td class="num">{row['RelativeAbundance']:.4f}%</td>
          <td class="bar-cell">
            <div class="bar" style="width:{bar_w}%"></div>
          </td>
        </tr>"""

    sufficient_count = stats["biomarkers_sufficient"]
    total_biomarkers = len(biomarkers)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Microbiome Analysis Report</title>
<style>
  /* ── Reset & Base ── */
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  :root {{
    --bg:       #f5f7fa;
    --surface:  #ffffff;
    --ink:      #1a1a2e;
    --muted:    #6b7280;
    --accent:   #2A9D8F;
    --accent2:  #E76F51;
    --ok:       #1a7a5e;
    --ok-bg:    #d1fae5;
    --low:      #b91c1c;
    --low-bg:   #fee2e2;
    --radius:   10px;
    --shadow:   0 2px 12px rgba(0,0,0,.07);
    --font:     'Segoe UI', system-ui, sans-serif;
  }}
  body {{
    font-family: var(--font);
    background: var(--bg);
    color: var(--ink);
    min-height: 100vh;
  }}

  /* ── Header ── */
  header {{
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
    color: #fff;
    padding: 2.5rem 2rem 2rem;
    position: relative;
    overflow: hidden;
  }}
  header::after {{
    content: '';
    position: absolute;
    inset: 0;
    background: url("data:image/svg+xml,%3Csvg width='60' height='60' viewBox='0 0 60 60' xmlns='http://www.w3.org/2000/svg'%3E%3Ccircle cx='30' cy='30' r='1.5' fill='%23ffffff08'/%3E%3C/svg%3E");
  }}
  header .inner {{ position: relative; max-width: 960px; margin: auto; }}
  header h1 {{ font-size: 1.9rem; font-weight: 700; letter-spacing: -.02em; }}
  header .meta {{ margin-top: .4rem; font-size: .85rem; opacity: .65; }}
  header .meta span {{ margin-right: 1.4rem; }}

  /* ── Main container ── */
  main {{ max-width: 960px; margin: 0 auto; padding: 2rem 1.5rem 4rem; }}

  /* ── Stat cards ── */
  .stat-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
    margin-bottom: 2rem;
  }}
  .stat-card {{
    background: var(--surface);
    border-radius: var(--radius);
    padding: 1.2rem 1.4rem;
    box-shadow: var(--shadow);
    border-top: 4px solid var(--accent);
  }}
  .stat-card.accent2 {{ border-color: var(--accent2); }}
  .stat-card .label {{ font-size: .72rem; text-transform: uppercase; letter-spacing: .08em; color: var(--muted); margin-bottom: .3rem; }}
  .stat-card .value {{ font-size: 1.8rem; font-weight: 700; color: var(--ink); line-height: 1; }}
  .stat-card .sub   {{ font-size: .78rem; color: var(--muted); margin-top: .25rem; }}

  /* ── Section ── */
  section {{ margin-bottom: 2.5rem; }}
  h2 {{
    font-size: 1.15rem; font-weight: 700; color: var(--ink);
    padding-bottom: .5rem; margin-bottom: 1rem;
    border-bottom: 2px solid #e5e7eb;
    display: flex; align-items: center; gap: .5rem;
  }}
  h2 .tag {{
    font-size: .7rem; background: var(--accent); color: #fff;
    padding: .15rem .5rem; border-radius: 99px; font-weight: 600;
    letter-spacing: .05em;
  }}

  /* ── Tables ── */
  .card {{
    background: var(--surface);
    border-radius: var(--radius);
    box-shadow: var(--shadow);
    overflow: hidden;
  }}
  table {{ width: 100%; border-collapse: collapse; font-size: .88rem; }}
  thead th {{
    background: #f9fafb; font-size: .72rem; text-transform: uppercase;
    letter-spacing: .07em; color: var(--muted);
    padding: .75rem 1rem; text-align: left; border-bottom: 1px solid #e5e7eb;
  }}
  tbody tr {{ border-bottom: 1px solid #f3f4f6; transition: background .15s; }}
  tbody tr:last-child {{ border-bottom: none; }}
  tbody tr:hover {{ background: #f9fafb; }}
  td {{ padding: .7rem 1rem; vertical-align: middle; }}
  .num   {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .muted {{ color: var(--muted); }}
  .bm-name   {{ font-weight: 600; }}
  .tax-label {{ max-width: 340px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}

  /* ── Badges ── */
  .badge {{
    display: inline-flex; align-items: center; gap: .3rem;
    padding: .25rem .75rem; border-radius: 99px;
    font-size: .78rem; font-weight: 600;
  }}
  .badge-ok  {{ background: var(--ok-bg);  color: var(--ok);  }}
  .badge-low {{ background: var(--low-bg); color: var(--low); }}

  /* ── Bar chart column ── */
  .bar-cell {{ width: 180px; }}
  .bar {{
    height: 8px; border-radius: 4px;
    background: linear-gradient(90deg, var(--accent), var(--accent2));
    min-width: 3px;
  }}

  /* ── Pie chart ── */
  .pie-wrap {{ background: var(--surface); border-radius: var(--radius); box-shadow: var(--shadow); padding: 1.5rem; text-align: center; }}
  .pie-wrap img {{ max-width: 100%; height: auto; }}

  /* ── Footer ── */
  footer {{ text-align: center; font-size: .78rem; color: var(--muted); margin-top: 3rem; }}
</style>
</head>
<body>

<header>
  <div class="inner">
    <h1>🧬 Microbiome Analysis Report</h1>
    <div class="meta">
      <span>📁 {input_filename or 'OTU Table'}</span>
      <span>🕐 Generated {ts}</span>
    </div>
  </div>
</header>

<main>

  <!-- ── Summary cards ── -->
  <div class="stat-grid">
    <div class="stat-card">
      <div class="label">Total Reads</div>
      <div class="value">{stats['total_reads']:,}</div>
      <div class="sub">raw abundance counts</div>
    </div>
    <div class="stat-card">
      <div class="label">Unique OTUs</div>
      <div class="value">{stats['unique_otus']:,}</div>
      <div class="sub">operational taxonomic units</div>
    </div>
    <div class="stat-card">
      <div class="label">Top Species</div>
      <div class="value" style="font-size:1.05rem;line-height:1.3">{stats['top_species_label']}</div>
      <div class="sub">{stats['top_species_pct']:.2f}% relative abundance</div>
    </div>
    <div class="stat-card accent2">
      <div class="label">Biomarkers Sufficient</div>
      <div class="value">{sufficient_count}/{total_biomarkers}</div>
      <div class="sub">above {STATUS_THRESHOLD}% threshold</div>
    </div>
  </div>

  <!-- ── Biomarker table ── -->
  <section>
    <h2>🔬 Key Biomarker Summary <span class="tag">5 targets</span></h2>
    <div class="card">
      <table>
        <thead>
          <tr>
            <th>Biomarker</th>
            <th style="text-align:right">Raw Counts</th>
            <th style="text-align:right">Relative Abundance</th>
            <th style="text-align:right">Matching OTUs</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody>{bm_rows_html}</tbody>
      </table>
    </div>
    <p style="font-size:.78rem;color:#6b7280;margin-top:.6rem">
      ℹ️ Status threshold: &gt; {STATUS_THRESHOLD}% relative abundance = <strong>Sufficient</strong>.
      Matching uses case-insensitive taxonomy search; multiple OTUs may map to one biomarker.
    </p>
  </section>

  <!-- ── Pie chart ── -->
  <section>
    <h2>🥧 Top 10 Species Distribution <span class="tag">relative abundance</span></h2>
    <div class="pie-wrap">
      <img src="{pie_src}" alt="Top 10 species pie chart">
    </div>
  </section>

  <!-- ── Top-10 table ── -->
  <section>
    <h2>📊 Top 10 Most Abundant Species</h2>
    <div class="card">
      <table>
        <thead>
          <tr>
            <th>Species / Taxon</th>
            <th style="text-align:right">Raw Counts</th>
            <th style="text-align:right">Relative Abundance</th>
            <th>Proportion</th>
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
    args = parser.parse_args()

    df       = load_data(args.input)
    analysis = analyze_microbiome(df)
    generate_html_report(
        analysis,
        output_path=args.output,
        input_filename=os.path.basename(args.input),
    )


if __name__ == "__main__":
    main()
