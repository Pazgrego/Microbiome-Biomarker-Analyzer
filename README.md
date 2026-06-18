# Microbiome Biomarker Analyzer

## Project Overview
Most people who take a gut microbiome test receive a confusing, massive spreadsheet (an OTU table) filled with complex bacterial names and raw numbers. To a regular person, this data is completely unreadable. 

This project is an automated data pipeline designed to bridge that gap. It takes raw, technical microbiome data and transforms it into a clean, easy-to-understand **HTML visual report designed for everyday users (patients, biohackers, or health enthusiasts)**. The tool simplifies the science, providing clear visual insights into their gut health without requiring a background in biology.

## How to run

1. **Install dependencies** (Python 3):

   ```bash
   pip install pandas matplotlib
   ```

2. **Generate the report** from an OTU table (CSV):

   ```bash
   python microbiome_pipeline.py --input docs/OTU_Table_P37.csv --output report.html
   ```

   Use your own file with `--input` / `-i`. Optional: `--output` / `-o` (default: `report.html`).

3. **Open the HTML** in your browser:

   ```bash
   open report.html
   ```


## Core Features

1. **OTU table ingestion and normalization** — Reads CSV files, requires `Name` and `Taxonomy` columns, auto-detects an abundance column (`Combined Abundance`, `Abundance`, or any column containing `"abundance"`), and supports an optional `--abundance` override.

2. **Relative abundance calculation** — Sums total reads and converts raw counts into per-OTU relative abundances (proportions of the sample). Rows with zero abundance are filtered out before analysis.

3. **Alpha diversity metrics** — Computes Shannon index (species richness and evenness) and Simpson diversity (1 − D), plus total sample depth (total reads) and observed richness (count of non-zero OTUs).

4. **Biomarker screening against biological thresholds** — Searches five gut-health taxa via fuzzy taxonomy matching and classifies each as **Sufficient** or **Low / Absent** based on combined relative abundance:

   | Biomarker | Threshold |
   |---|---|
   | *Faecalibacterium prausnitzii* | 1.0% |
   | *Eubacterium sp* | 1.0% |
   | *Akkermansia muciniphila* | 0.3% |
   | *Roseburia sp* | 0.5% |
   | *Ruminococcus bromii* | 0.5% |

5. **Top-10 composition profiling** — Ranks OTUs by abundance, extracts human-readable labels from taxonomy strings (species name, or genus + `sp.`), and builds a top-10 species table with percentages.

6. **Interactive visualization** — Generates a Plotly donut chart of the top 10 taxa with hover tooltips showing proportions, embedded in the output HTML via Plotly's CDN.

7. **Standalone HTML report generation** — Assembles everything into a single self-contained `report.html` with sequencing overview and diversity metrics (with explanatory tooltips), a biomarker status table (color-coded badges and threshold tooltips), a clinical encyclopedia with descriptions and PubMed reference links, the interactive composition chart, a top-10 abundance table, and introductory educational text about the microbiome.


---

## Technical Architecture & Dependencies

The project is a single Python script (`microbiome_pipeline.py`) with three main stages:

1. **`load_data()`** — Ingests a CSV file, validates required columns, and normalizes abundance data.
2. **`analyze_microbiome()`** — Computes relative abundances, alpha diversity, biomarker status, top-10 taxa, and a Plotly chart snippet.
3. **`generate_html_report()`** — Renders a standalone HTML report using Python f-strings and embedded CSS.

### Dependencies

**Required** (Python 3):

```bash
pip install pandas numpy plotly
```

| Package | Role |
|---|---|
| `pandas` | Load OTU tables, filter and aggregate taxa, rank abundances |
| `numpy` | Shannon and Simpson diversity calculations |
| `plotly` | Interactive top-10 donut chart (served via Plotly CDN in the HTML output) |


### Expected Input Format

The pipeline expects a **CSV** file with at least `Name`, `Taxonomy`, and an abundance column (`Abundance`, `Combined Abundance`, or any column whose name contains `"abundance"`):

```csv
Name,Taxonomy,Abundance
OTU_001,"Bacteria; Firmicutes; Clostridia; Clostridiales; Ruminococcaceae; Faecalibacterium prausnitzii",32781
OTU_005,"Bacteria; Verrucomicrobia; Verrucomicrobiae; Verrucomicrobiales; Akkermansiaceae; Akkermansia muciniphila",6100
...
```

Use `--abundance <column>` if your file has multiple abundance columns and auto-detection should not be used.