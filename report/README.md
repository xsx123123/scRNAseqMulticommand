# scRNAseqMulticommand Quarto Report

This Quarto project generates an interactive HTML report from scRNAseqMulticommand pipeline outputs.

## Quick Start

1. Link or copy your pipeline output directory to `report/data/`:

```bash
cd report
ln -s /path/to/pipeline/output data
```

2. Render the report:

```bash
quarto render
```

3. Open `_site/index.html` in your browser.

## Project Structure

- `index.qmd` — Project overview dashboard
- `01-qc.qmd` — Quality control results
- `02-integration.qmd` — Batch effect and integration
- `03-clustering.qmd` — Clustering and marker genes
- `04-annotation.qmd` — Cell type annotation
- `05-deg.qmd` — Differential expression
- `06-methods.qmd` — Methods and parameters
- `R/` — Helper R functions
- `styles/` — Custom CSS/SCSS themes
- `_quarto.yml` — Quarto project configuration

## Data Contract

The report expects the pipeline output directory to contain:

- `manifest.json` — Global summary of the analysis run
- `summary.json` files in each step directory (QC, doublet, ambient_rna, etc.)
- Plots and CSV files referenced by the summaries

See `docs/QUARTO_REPORT_ARCHITECTURE.md` for the full data contract.
