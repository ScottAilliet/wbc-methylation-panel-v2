# WBC Methylation Panel v2 — Primer Design Pipeline

A complete pipeline for designing bisulfite-PCR primers for WBC methylation markers,
based on the Loyfer et al. (Nature 2023) human methylation atlas (GSE186458).

## Features

- **DMR discovery**: Top 300 DMR regions per immune cell type, using the full 82-cell-type human atlas as background
- **Per-CpG resolution**: Individual CpG labels, positions, and beta values for all 207 samples
- **Primer3Plus-compatible**: All 161 Primer3Plus parameters configurable
- **6-strand bowtie specificity**: Screens against unconverted, converted-methylated, and converted-unmethylated genomes (sense + antisense)
- **Secondary structure screening**: Hairpin MFE with -1.5 kcal/mol cutoff
- **Common SNP screening**: dbSNP-based screening with 3' end priority
- **DimerDetective**: Primer-dimer prediction with -1.0 kcal/mol conservative cutoff
- **U-assays output format**: XLSX (27 columns) + PDF (one page per primer pair)
- **Modular execution**: Run any combination of pipeline steps
- **macOS compatible**: Runs locally without biomni

## Quick Start

```bash
# 1. Install dependencies
chmod +x install_macos.sh
./install_macos.sh

# 2. Download data (~15 GB)
chmod +x download_data.sh
./download_data.sh

# 3. Run the full pipeline
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/

# Or run specific steps only
python -m methyl_panel.pipeline --steps 1,2,3 \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62
```

## Pipeline Steps

| Step | Module | Description |
|------|--------|-------------|
| 1 | `phase1_dmr_loader` | Load DMR blocks from Excel |
| 2 | `phase2_bisulfite_convert` | Bisulfite convert to 6 strands |
| 3 | `phase3_primer3_design` | Design primers with Primer3 |
| 4 | `phase4_bowtie_specificity` | 6-strand bowtie2 specificity |
| 5 | `phase5_structure` | Secondary structure screening |
| 6 | `phase6_snp` | Common SNP screening |
| 7 | `phase7_dimer` | DimerDetective primer-dimer prediction |
| 8 | `output_xlsx` | U-assays-style XLSX output |
| 9 | `output_pdf` | U-assays-style PDF output |

## Configuration

The pipeline uses `Primer3PlusConfig` (in `config.py`) which supports all 161 Primer3Plus parameters.
Load settings from a Primer3Plus settings file:

```python
from methyl_panel.config import Primer3PlusConfig
config = Primer3PlusConfig.from_primer3plus_file("Primer3Plus_Settings.txt")
```

Or configure programmatically:

```python
config = Primer3PlusConfig()
config.primer_min_tm = 58.0  # REQUIRED — no defaults
config.primer_opt_tm = 60.0
config.primer_max_tm = 62.0
config.validate_tm()  # Raises if Tm not set
```

### Salt Presets

- **primer3plus**: 50mM Na, 3mM Mg, 1.2mM dNTP, 250nM primer (SantaLucia correction)
- **roche_dlc**: 100mM Na, 4.5mM Mg, 1.2mM dNTP, 250nM primer (Owczarzy correction)

## Cell Types

| Cell type | ID | Markers | Background |
|---|---|---|---|
| Monocytes | MONO | 290 | All 81 other cell types |
| B cells | BCELL | 300 | All 81 other cell types |
| NK | NK | 300 | All 81 other cell types |
| Granulocytes | GRAN | 300 | All 81 other cell types |
| CD3 T | CD3T | 300 | 73 non-T cell types |
| CD8 T | CD8T | 300 | Non-T + CD3 + CD4 lineage |
| CD4 T | CD4T | 300 | Non-T + CD3 + CD8 lineage |

## Data Source

Loyfer et al., "A human methylation atlas for normal and cancer cell types",
Nature 617, 690–697 (2023). GEO accession: GSE186458.

## License

MIT
