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
from html import escape as html_escape
import io
import json
import math
import os
import sys
import warnings
from datetime import datetime

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    # Only match rows explicitly indicating Eubacterium, avoiding 'Pseudoeubacterium' or unrelated matches
    "Eubacterium sp":          ["eubacterium"],
    "Faecalibacterium prausnitzii": ["faecalibacterium", "prausnitzii"],
    "Roseburia sp":            ["roseburia"],
}

# Custom thresholds mapping based on biological relevance (values in proportion: 0.01 = 1.0%)
BIOMARKER_THRESHOLDS = {
    "Faecalibacterium prausnitzii": 0.010,  # High threshold for a highly dominant species (1.0%)
    "Eubacterium sp":              0.010,  # High threshold for a dominant core genus (1.0%)
    "Akkermansia muciniphila": 0.003,      # Lower threshold for a specialized mucosal niche species (0.3%)
    "Roseburia sp":            0.005,      # Medium threshold for butyrate producers (0.5%)
    "Ruminococcus bromii":     0.005       # Medium threshold for key starch degraders (0.5%)
}

# Default threshold for any other biomarkers not explicitly defined in BIOMARKER_THRESHOLDS
DEFAULT_BIOMARKER_THRESHOLD = 0.001        # (0.1%)

# Educational copy for the report (key = biomarker display name in BIOMARKER_SEARCH)
BIOMARKER_ROLES = {
    "Akkermansia muciniphila": {
        "title": "Akkermansia muciniphila",
        "body": (
            "This specialized bacterium is a crucial gatekeeper of gut barrier integrity. "
            "It uniquely resides in and feeds on the mucus layer lining the intestinal walls. "
            "By degrading this mucus, it stimulates the body to constantly produce fresh, healthy mucus, "
            "effectively strengthening the gut barrier and preventing toxins from leaking into the bloodstream "
            '(a condition often referred to as "leaky gut"). A robust population of Akkermansia muciniphila '
            "is strongly associated with healthy metabolic functions, reduced systemic inflammation, and a lower "
            "risk of developing obesity and type 2 diabetes."
        ),
    },
    "Ruminococcus bromii": {
        "title": "Ruminococcus bromii",
        "body": (
            'Commonly described as a "keystone species" in the human colon, this bacterium possesses a highly '
            "specialized enzymatic machinery required to break down resistant starch and complex dietary fibers "
            "that our bodies cannot digest on their own. By initiating the degradation of these tough carbohydrates, "
            "Ruminococcus bromii acts as a primary recycler in the gut ecosystem, breaking food down into simpler "
            "sugars. This process releases vital nutrients and paves the way for other beneficial bacteria to "
            "consume those leftovers and convert them into protective compounds."
        ),
    },
    "Eubacterium sp": {
        "title": "Eubacterium group (Eubacterium sp.)",
        "body": (
            "This core group of core intestinal bacteria plays a massive role in maintaining a balanced gut microbiome "
            "through the fermentation of complex plant polysaccharides and prebiotic fibers. As they ferment these "
            "fibers, they generate large quantities of short-chain fatty acids (SCFAs), which lower the pH of the "
            "colon, making it an unwelcoming environment for harmful pathogens. Furthermore, specific species within "
            "this genus have been shown to actively participate in anti-inflammatory pathways, support metabolic "
            "homeostasis, and influence the gut-brain axis, thereby contributing to overall cognitive health."
        ),
    },
    "Faecalibacterium prausnitzii": {
        "title": "Faecalibacterium prausnitzii",
        "body": (
            "Widely recognized as one of the most abundant and vital pillars of a healthy adult gut microbiome, "
            "this bacterium serves as a primary factory for butyrate, a critical short-chain fatty acid. Butyrate "
            "acts as the main fuel source for the cells lining the colon, keeping them healthy and functional. "
            "Beyond nourishment, Faecalibacterium prausnitzii exerts powerful, direct anti-inflammatory effects by "
            "modulating the immune system and blocking inflammatory signaling pathways. Low levels of this bacterium "
            "are consistently documented in individuals suffering from inflammatory bowel diseases (IBD), metabolic "
            "disorders, and depression."
        ),
    },
    "Roseburia sp": {
        "title": "Roseburia group (Roseburia sp.)",
        "body": (
            "This genus consists of highly active carbohydrate-fermenting bacteria that work in close harmony with "
            "other fiber degraders to produce butyrate. By generating this fatty acid, Roseburia species help maintain "
            "the mucosal lining, enhance intestinal motility, and provide the energy necessary to sustain a robust "
            "immune defense within the gut. Their prevalence is strongly tied to a fiber-rich, Mediterranean-style diet, "
            "and maintaining an optimal abundance of Roseburia is essential for preventing low-grade chronic "
            "inflammation, metabolic syndrome, and irritable bowel symptoms."
        ),
    },
}

# ── Helper functions ───────────────────────────────────────────────────────────

def _normalise_columns(df: pd.DataFrame, abundance_column=None) -> pd.DataFrame:
    """Map varied column names to Name, Taxonomy, Abundance."""
    col_map: dict[str, str] = {}
    for col in df.columns:
        low = col.lower().strip()
        if low in ("name", "otu", "otu id", "otu_id", "id"):
            col_map[col] = "Name"
        elif low in ("taxonomy", "taxon", "lineage"):
            col_map[col] = "Taxonomy"

    df = df.rename(columns=col_map)

    if abundance_column:
        if abundance_column not in df.columns:
            print(
                f"Error: abundance column '{abundance_column}' not found.",
                file=sys.stderr,
            )
            print(f"Found columns: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)
        df["Abundance"] = df[abundance_column]
    elif "Abundance" not in df.columns:
        combined = next(
            (c for c in df.columns if c.lower().strip() == "combined abundance"),
            None,
        )
        if combined:
            df["Abundance"] = df[combined]
        else:
            abundance_cols = [
                c for c in df.columns
                if "abundance" in c.lower() or c.lower() in ("count", "counts", "reads")
            ]
            if len(abundance_cols) == 1:
                df["Abundance"] = df[abundance_cols[0]]
            elif abundance_cols:
                print(
                    "Error: multiple abundance columns found; use --abundance to pick one.",
                    file=sys.stderr,
                )
                print(f"Options: {abundance_cols}", file=sys.stderr)
                sys.exit(1)

    for required in ("Name", "Taxonomy", "Abundance"):
        if required not in df.columns:
            df[required] = ""

    return df[["Name", "Taxonomy", "Abundance"]]


def load_data(file_path: str, abundance_column=None) -> pd.DataFrame:
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
    df = _normalise_columns(df, abundance_column=abundance_column)

    missing = {"Name", "Taxonomy", "Abundance"} - set(df.columns)
    if missing:
        print(f"Error: Missing required columns {missing}.", file=sys.stderr)
        print(f"Found columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    df["Abundance"] = pd.to_numeric(df["Abundance"], errors="coerce").fillna(0)
    df = df[df["Abundance"] > 0].copy()
    df.reset_index(drop=True, inplace=True)
    print(f"[load_data] Loaded {len(df):,} OTU rows from '{file_path}'.")
    return df


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
    # try species first
    if len(parts) >= 7 and parts[6] and parts[6].lower() != "unclassified":
        return parts[6]
    if len(parts) >= 6 and parts[5] and parts[5].lower() != "unclassified":
        return f"{parts[5]} sp."
    return str(row["Name"])


def analyze_microbiome(df: pd.DataFrame) -> dict:
    """Executes the core analytics pipeline."""
    total_reads = df["Abundance"].sum()
    if total_reads <= 0:
        print("[Warning] Total Abundance sum is 0. Check your data input.")
        df["RelativeAbundance"] = 0.0
    else:
        df["RelativeAbundance"] = df["Abundance"] / total_reads

    # Remove 0-abundance rows for internal statistics
    df_enriched = df[df["Abundance"] > 0].copy()
    print(f"[load_data] Loaded {len(df_enriched)} OTU rows.")

    # ── Biomarker extraction with custom thresholds ──
    biomarker_results = []
    for display_name, search_terms in BIOMARKER_SEARCH.items():
        matched = _match_rows(df_enriched, search_terms)
        agg_abundance = matched["Abundance"].sum()
        agg_pct = matched["RelativeAbundance"].sum().round(4)
        
        # Retrieve the custom threshold for the current biomarker, fallback to default if not listed
        current_threshold = BIOMARKER_THRESHOLDS.get(display_name, DEFAULT_BIOMARKER_THRESHOLD)
        
        # Evaluate the status dynamically using the customized threshold
        status = "Sufficient" if agg_pct >= current_threshold else "Low / Absent"
        
        biomarker_results.append(
            {
                "name": display_name,
                "abundance": agg_abundance,
                "relative_pct": agg_pct,
                "status": status,
                "threshold": current_threshold,
                "match_count": len(matched),
            }
        )

    # ── Top 10 profiling ──
    top10 = df_enriched.nlargest(10, "Abundance").copy()
    if not top10.empty:
        top10["Label"] = top10.apply(_pretty_label, axis=1)
    else:
        top10["Label"] = []

    pie_base64 = ""
    if not top10.empty:
        colors = plt.cm.get_cmap("tab20c")(range(len(top10)))
        labels = top10["Label"].tolist() if "Label" in top10.columns else top10["Name"].tolist()

        legend_labels = []
        for label in labels:
            short = str(label)
            if len(short) > 32:
                short = short[:29] + "..."
            legend_labels.append(short)

        fig = plt.figure(figsize=(10, 5.2))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.12], wspace=0.04)
        ax_pie = fig.add_subplot(gs[0, 0])
        ax_leg = fig.add_subplot(gs[0, 1])
        ax_leg.axis("off")

        ring_width = 0.55
        ring_mid = 1 - ring_width / 2
        min_slice_label_pct = 9.0

        def _slice_pct_label(pct):
            return f"{pct:.1f}%" if pct >= min_slice_label_pct else ""

        wedges, _, autotexts = ax_pie.pie(
            top10["Abundance"],
            colors=colors,
            startangle=140,
            wedgeprops=dict(width=ring_width, edgecolor="w"),
            autopct=_slice_pct_label,
            textprops={"fontsize": 10, "fontweight": "bold", "color": "white"},
        )
        ax_pie.set_aspect("equal")
        for wedge, autotext in zip(wedges, autotexts):
            if not autotext.get_text():
                autotext.set_visible(False)
                continue
            autotext.set_color("white")
            angle = math.radians((wedge.theta1 + wedge.theta2) / 2)
            autotext.set_position((
                ring_mid * math.cos(angle),
                ring_mid * math.sin(angle),
            ))
            autotext.set_horizontalalignment("center")
            autotext.set_verticalalignment("center")

        ax_pie.text(
            0, 0, "Top 10\nSpecies",
            ha="center", va="center",
            fontsize=14, fontweight="bold", color="#1a1a2e", linespacing=1.15,
        )

        ax_leg.legend(
            wedges,
            legend_labels,
            title="Taxa",
            loc="center left",
            fontsize=11,
            title_fontsize=12,
            frameon=False,
            borderaxespad=0,
        )

        fig.subplots_adjust(left=0.02, right=0.99, top=0.98, bottom=0.02)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, facecolor="white", pad_inches=0.05)
        plt.close(fig)
        buf.seek(0)
        pie_base64 = base64.b64encode(buf.read()).decode("utf-8")

    return {
        "total_reads": total_reads,
        "distinct_otus": len(df_enriched),
        "biomarkers": biomarker_results,
        "top10": top10.to_dict(orient="records"),
        "pie_base64": pie_base64,
    }


def generate_html_report(analysis: dict, output_path: str):
    """Generates a highly styled, standalone HTML report with embedded visualizations."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Format Biomarker rows
    biomarker_rows_html = ""
    for b in analysis["biomarkers"]:
        status_class = "status-sufficient" if b["status"] == "Sufficient" else "status-low"
        pct_display = f"{b['relative_pct'] * 100:.2f}%"
        thresh_pct = b["threshold"] * 100
        tip = f"Threshold: ≥ {thresh_pct:.2f}% relative abundance"
        tip_attr = html_escape(tip, quote=True)
        role_key = html_escape(b["name"], quote=True)
        biomarker_rows_html += f"""
        <tr class="biomarker-row" data-role-key="{role_key}" tabindex="0" role="button"
            aria-expanded="false" aria-controls="biomarker-detail">
          <td class="biomarker-taxon-cell"><strong>{html_escape(b['name'])}</strong></td>
          <td style="text-align:right">{b['abundance']:,}</td>
          <td style="text-align:right; font-weight:600;">{pct_display}</td>
          <td class="status-cell">
            <span class="status-badge has-threshold {status_class}"
                  data-tip="{tip_attr}" title="{tip_attr}" tabindex="0"
                  aria-label="{html_escape(b['status'] + '. ' + tip, quote=True)}">
              {b['status']}<span class="status-info-mark" aria-hidden="true"></span>
            </span>
          </td>
        </tr>"""

    biomarker_roles_json = json.dumps(BIOMARKER_ROLES)

    # Format Top 10 rows
    top10_rows_html = ""
    for idx, r in enumerate(analysis["top10"]):
        pct_display = f"{r['RelativeAbundance'] * 100:.2f}%"
        short_tax = str(r.get('Taxonomy', ''))
        if len(short_tax) > 65:
            short_tax = short_tax[:62] + "..."
        
        top10_rows_html += f"""
        <tr>
          <td>
            <div style="font-weight:600">{r.get('Label', r['Name'])}</div>
            <div style="font-size:11px; color:#666;" title="{r.get('Taxonomy','')}">{short_tax}</div>
          </td>
          <td style="text-align:right">{r['Abundance']:,}</td>
          <td style="text-align:right; font-weight:600;">{pct_display}</td>
          <td>
            <div style="background:#eee; border-radius:3px; height:8px; width:100px; overflow:hidden;">
              <div style="background:#4f46e5; height:100%; width:{r['RelativeAbundance']*100}%"></div>
            </div>
          </td>
        </tr>"""

    pie_src = f"data:image/png;base64,{analysis['pie_base64']}" if analysis["pie_base64"] else ""

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
    max-width: 900px;
    margin: 0 auto 30px auto;
    border-bottom: 1px solid #e5e7eb;
    padding-bottom: 20px;
  }}
  h1 {{ margin: 0 0 5px 0; font-size: 28px; color: #1e1b4b; }}
  .subtitle {{ color: #4b5563; font-size: 14px; margin: 0; }}
  main {{
    max-width: 900px;
    margin: 0 auto;
    display: grid;
    grid-template-columns: 1fr;
    gap: 30px;
  }}
  @media(min-width: 768px) {{
    main {{ grid-template-columns: 1fr 1fr; }}
    .full-width {{ grid-column: span 2; }}
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
  }}
  .metrics-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 15px;
  }}
  .metric-box {{
    background: #f3f4f6;
    padding: 15px;
    border-radius: 6px;
    text-align: center;
  }}
  .metric-val {{ font-size: 24px; font-weight: bold; color: #4f46e5; margin-bottom: 4px; }}
  .metric-lbl {{ font-size: 12px; color: #6b7280; text-transform: uppercase; letter-spacing: 0.5px; }}
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
  .status-cell {{ position: relative; overflow: visible; }}
  .status-badge {{
    display: inline-block;
    padding: 3px 8px;
    border-radius: 12px;
    font-size: 11px;
    font-weight: 600;
    text-transform: uppercase;
  }}
  .status-badge.has-threshold {{
    position: relative;
    display: inline-flex;
    align-items: center;
    gap: 5px;
    cursor: default;
  }}
  .status-info-mark {{
    font-size: 11px;
    line-height: 1;
    flex-shrink: 0;
    opacity: 0.5;
    transition: opacity 0.15s ease;
  }}
  .status-info-mark::before {{
    content: "i";
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 13px;
    height: 13px;
    border-radius: 50%;
    border: 1px solid currentColor;
    font-size: 9px;
    font-weight: 700;
    font-style: italic;
  }}
  .status-badge.has-threshold:hover .status-info-mark,
  .status-badge.has-threshold:focus-visible .status-info-mark {{
    opacity: 0.85;
  }}
  .status-badge.has-threshold:hover,
  .status-badge.has-threshold:focus-visible {{
    outline: none;
  }}
  .status-badge.has-threshold::after {{
    content: attr(data-tip);
    position: absolute;
    left: 50%;
    bottom: calc(100% + 10px);
    transform: translateX(-50%);
    max-width: 280px;
    white-space: normal;
    text-align: center;
    text-transform: none;
    font-size: 11px;
    font-weight: 500;
    line-height: 1.4;
    letter-spacing: 0;
    padding: 8px 10px;
    border-radius: 6px;
    background: #1f2937;
    color: #f9fafb;
    box-shadow: 0 4px 14px rgba(0, 0, 0, 0.15);
    pointer-events: none;
    opacity: 0;
    visibility: hidden;
    transition: opacity 0.15s ease, visibility 0.15s ease;
    z-index: 20;
  }}
  .status-badge.has-threshold:hover::after,
  .status-badge.has-threshold:focus-visible::after {{
    opacity: 1;
    visibility: visible;
  }}
  .status-sufficient {{ background: #dcfce7; color: #166534; }}
  .status-low {{ background: #fee2e2; color: #991b1b; }}
  .biomarker-card {{ overflow: visible; }}
  .biomarker-hint {{
    margin: -4px 0 12px 0;
    font-size: 13px;
    color: #6b7280;
  }}
  .biomarker-row {{
    cursor: pointer;
    transition: background 0.12s ease;
  }}
  .biomarker-row:hover {{
    background: #f3f4f6;
  }}
  .biomarker-row.is-selected {{
    background: #eef2ff;
  }}
  .biomarker-row.is-selected .biomarker-taxon-cell strong {{
    color: #4f46e5;
  }}
  .biomarker-taxon-cell strong::after {{
    content: " ▸";
    font-size: 11px;
    color: #9ca3af;
    font-weight: 400;
  }}
  .biomarker-row.is-selected .biomarker-taxon-cell strong::after {{
    content: " ▾";
    color: #4f46e5;
  }}
  .biomarker-detail {{
    display: none;
    margin: 0;
    padding: 16px 18px;
    border-top: 1px solid #e5e7eb;
    background: #f8fafc;
  }}
  .biomarker-detail.is-open {{
    display: block;
  }}
  .biomarker-detail h4 {{
    margin: 0 0 10px 0;
    font-size: 15px;
    font-weight: 600;
    color: #4f46e5;
  }}
  .biomarker-detail p {{
    margin: 0;
    font-size: 13px;
    line-height: 1.6;
    color: #374151;
  }}
  .biomarker-detail-placeholder {{
    margin: 0;
    padding: 14px 18px;
    border-top: 1px solid #e5e7eb;
    font-size: 13px;
    color: #9ca3af;
    font-style: italic;
    background: #fafafa;
  }}
  .composition-section h2 {{
    font-size: 22px;
    margin-bottom: 10px;
  }}
  .composition-card {{
    padding: 8px 12px 12px;
    overflow: hidden;
  }}
  .composition-card .img-wrap {{
    width: 100%;
    line-height: 0;
  }}
  .composition-card .img-wrap img {{
    display: block;
    width: 100%;
    height: auto;
    max-height: 520px;
    object-fit: contain;
    object-position: center;
  }}
  footer {{
    max-width: 900px;
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
  <h1>Microbiome Analysis Report</h1>
  <p class="subtitle">Automated high-throughput OTU profiling pipeline</p>
</header>

<main>

  <section class="full-width">
    <h2>Sequencing Overview</h2>
    <div class="card metrics-grid">
      <div class="metric-box">
        <div class="metric-val">{analysis['total_reads']:,}</div>
        <div class="metric-lbl">Total Sample Counts (Depth)</div>
      </div>
      <div class="metric-box">
        <div class="metric-val">{analysis['distinct_otus']:,}</div>
        <div class="metric-lbl">Observed Richness (Active OTUs)</div>
      </div>
    </div>
  </section>

  <section class="full-width">
    <h2>Health Biomarkers</h2>
    <p class="biomarker-hint">Click a target taxon to learn about its role in gut health.</p>
    <div class="card biomarker-card" style="padding:0;">
      <table class="biomarker-table">
        <thead>
          <tr>
            <th>Target Taxon</th>
            <th style="text-align:right">Raw Abundance</th>
            <th style="text-align:right">Relative Abundance</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody>{biomarker_rows_html}</tbody>
      </table>
      <p id="biomarker-detail-placeholder" class="biomarker-detail-placeholder">
        Select a target taxon above to read about its role.
      </p>
      <div id="biomarker-detail" class="biomarker-detail" aria-live="polite">
        <h4 id="biomarker-detail-title"></h4>
        <p id="biomarker-detail-body"></p>
      </div>
    </div>
    <script type="application/json" id="biomarker-roles-json">{biomarker_roles_json}</script>
    <script>
    (function () {{
      const roles = JSON.parse(document.getElementById("biomarker-roles-json").textContent);
      const panel = document.getElementById("biomarker-detail");
      const placeholder = document.getElementById("biomarker-detail-placeholder");
      const titleEl = document.getElementById("biomarker-detail-title");
      const bodyEl = document.getElementById("biomarker-detail-body");
      let selectedRow = null;

      function closeDetail() {{
        selectedRow = null;
        document.querySelectorAll(".biomarker-row").forEach((r) => {{
          r.classList.remove("is-selected");
          r.setAttribute("aria-expanded", "false");
        }});
        panel.classList.remove("is-open");
        placeholder.style.display = "block";
      }}

      function openDetail(row) {{
        const key = row.getAttribute("data-role-key");
        const role = roles[key];
        if (!role) return;

        if (selectedRow === row) {{
          closeDetail();
          return;
        }}

        document.querySelectorAll(".biomarker-row").forEach((r) => {{
          r.classList.remove("is-selected");
          r.setAttribute("aria-expanded", "false");
        }});
        row.classList.add("is-selected");
        row.setAttribute("aria-expanded", "true");
        selectedRow = row;

        titleEl.textContent = role.title;
        bodyEl.textContent = role.body;
        placeholder.style.display = "none";
        panel.classList.add("is-open");
        panel.scrollIntoView({{ behavior: "smooth", block: "nearest" }});
      }}

      document.querySelectorAll(".biomarker-row").forEach((row) => {{
        row.addEventListener("click", (e) => {{
          if (e.target.closest(".status-cell")) return;
          openDetail(row);
        }});
        row.addEventListener("keydown", (e) => {{
          if (e.target.closest(".status-cell")) return;
          if (e.key === "Enter" || e.key === " ") {{
            e.preventDefault();
            openDetail(row);
          }}
        }});
      }});
    }})();
    </script>
  </section>

  <section class="full-width composition-section">
    <h2>Composition Profile</h2>
    <div class="card composition-card">
      <div class="img-wrap">
        <img src="{pie_src}" alt="Top 10 species pie chart — relative abundance">
      </div>
    </div>
  </section>

  <section class="full-width">
    <h2>Top 10 Most Abundant Species</h2>
    <div class="card" style="padding:0; overflow:hidden;">
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