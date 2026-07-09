# WBC Methylation Panel v2 — Operator Handbook

**Version:** Pipeline v2.0.0
**Repository:** `wbc-methylation-panel-v2`
**Last updated:** 2026

---

## PREFACE

### What This Handbook Is

This handbook is the complete guide for running the WBC Methylation Panel v2 primer design pipeline on your Mac. It covers installation, data download, running the pipeline, and understanding the output — written for a lab scientist who is not a programmer.

### Who This Handbook Is For

You are a doctor, medical researcher, or lab scientist. You know molecular biology (PCR, DNA, methylation) but you are not a programmer. This handbook explains the computing concepts you need, in plain language with analogies to things you already know.

If you want to understand the science behind the pipeline, see [Materials and Methods](materials_and_methods.md). This handbook focuses on how to run it.

---

## Chapter 0 — Computing Basics for Lab Scientists

### 0.1 The Terminal (Command Line)

**What it is:** The Terminal is a text-based way to control your computer. Instead of clicking icons, you type commands. On Mac, open it from Applications → Utilities → Terminal (or press Cmd+Space and type "Terminal").

**The analogy:** Think of the Terminal as a very literal lab assistant. You tell it exactly what to do, step by step. If you type the wrong instruction, it does the wrong thing — but it never guesses.

**What you see when you open it:**

```
scott@MacBook-Pro ~ %
```

The `%` symbol is the **prompt** — it means "I am ready for your command." You type after it and press Enter.

**Your first commands:**

```bash
pwd
```
"Print working directory" — tells you which folder you are in. Like asking "which lab am I standing in?"

```bash
ls
```
"List" — shows files and folders in your current location. Like looking at the shelves.

```bash
ls -la
```
Options (flags) always start with a dash. `-l` means "long format", `-a` means "all files including hidden ones".

---

### 0.2 Folders, Paths, and Navigation

**Folders are called directories** in terminal language. Same thing.

**A path is the address of a file or folder.** Two kinds:

- **Absolute path:** Starts from the root of your hard drive. Always starts with `/`.
  - `/Users/scottailliet/Documents/data/` — like a full postal address
  - Your home folder: `/Users/scottailliet/` — abbreviated as `~`

- **Relative path:** Starts from wherever you are now.
  - `data/beta_files/` means "from here, go into data, then beta_files"
  - `.` means "this folder"
  - `..` means "the parent folder (one level up)"

**Moving between folders:**

```bash
cd Documents          # Move into the Documents folder
cd ..                 # Go up one level
cd ~                  # Go back to your home folder
cd /Users/scottailliet/wbc-methylation-panel-v2  # Jump to a specific folder
```

**Pro tip:** Type the first few letters of a folder name, then press **Tab** — the terminal auto-completes it. This saves time and prevents typos.

---

### 0.3 Virtual Environments (Instead of Conda)

This pipeline uses a Python **virtual environment** (`.venv`) instead of conda. A virtual environment is an isolated Python setup that keeps this project's packages separate from everything else on your Mac.

**The analogy:** Like a clean bench in the lab — this project's reagents stay on this bench and don't mix with other projects.

**Creating and using it:**

```bash
# Create the virtual environment (done once by install_macos.sh)
python3.11 -m venv .venv

# Activate it (every time you open a new Terminal)
source .venv/bin/activate

# Your prompt changes to show you're in the environment:
(.venv) scott@MacBook-Pro wbc-methylation-panel-v2 %

# Deactivate (exit the environment)
deactivate
```

**Key point:** You must activate the environment every time you open a new Terminal window and want to use the pipeline. The environment does not persist between Terminal sessions.

---

### 0.4 Installing Software

The `install_macos.sh` script handles everything automatically. It installs:
- **Homebrew** (Mac package manager) — if not already installed
- **System tools:** samtools, bowtie2, tabix, bedtools
- **Python packages:** primer3-py, openpyxl, pandas, numpy, reportlab, pysam
- **wgbstools** (from GitHub source — not on PyPI)

```bash
chmod +x install_macos.sh    # Make the script executable (one-time)
./install_macos.sh           # Run it
```

`chmod +x` means "make this file executable" — you only do this once. `./` means "run the file in this folder."

---

### 0.5 Running Commands and Scripts

**Running the pipeline:**

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

The `\` at the end of a line means "this command continues on the next line." You can also type it all on one line.

**Arguments (flags):**
- `--steps all` means "run all 9 steps"
- `--steps 1,2,3` means "run only steps 1, 2, and 3"
- `--min-tm 58` means "set the minimum Tm to 58°C"
- `--help` shows all available arguments

```bash
python -m methyl_panel.pipeline --help
```

---

### 0.6 Files You Will Encounter

| Extension | What it is | How to open |
|-----------|-----------|-------------|
| `.py` | Python script (code) | Text editor or run with `python` |
| `.md` | Markdown document (like this handbook) | Text editor, or rendered in GitHub |
| `.xlsx` | Excel spreadsheet | Microsoft Excel, Google Sheets |
| `.bed.gz` | Compressed genomic coordinates file | `gunzip` to decompress |
| `.beta` | Methylation data file (Loyfer format) | Processed by wgbstools |
| `.csv` | Comma-separated values | Excel or text editor |
| `.json` | Data file (JavaScript Object Notation) | Text editor or Python |
| `.fa.gz` | Compressed genome FASTA | Used by samtools/bowtie2 |
| `.vcf.gz` | Compressed variant file (dbSNP) | Used by tabix |
| `.pdf` | PDF report (pipeline output) | Any PDF viewer |

---

### 0.7 Common Mistakes and How to Avoid Them

| Mistake | What happens | Fix |
|---------|-------------|-----|
| Not in the right folder | "File not found" errors | Run `pwd` to check; `cd` to the right place |
| Forgot to activate virtual environment | "Module not found" errors | Run `source .venv/bin/activate` |
| Typo in a command | "Command not found" | Use Tab-completion; copy-paste from this handbook |
| Wrong file path | "No such file or directory" | Check with `ls`; use absolute paths if unsure |
| Missing Tm values | "Tm values are REQUIRED" error | Always provide `--min-tm`, `--opt-tm`, `--max-tm` |

---

### 0.8 Quick Reference: The 5 Commands You Will Use Most

| Command | What it does | When to use it |
|---------|-------------|----------------|
| `source .venv/bin/activate` | Enter the pipeline environment | Every time you open a new Terminal |
| `./install_macos.sh` | Install all dependencies | First time only |
| `./download_data.sh` | Download all data (~15 GB) | First time only |
| `python -m methyl_panel.pipeline --steps all ...` | Run the pipeline | To design primers |
| `python -m methyl_panel.pipeline --help` | Show all options | When unsure |

---

## Chapter 1 — Installation

### 1.1 Prerequisites

- **macOS** (Intel or Apple Silicon)
- **Xcode Command Line Tools:** Install with `xcode-select --install` (a dialog will appear — click Install)
- **~20 GB free disk space** (15 GB for data, 5 GB for tools and genome)

### 1.2 Step-by-Step Installation

```bash
# 1. Clone the repository from GitHub
git clone https://github.com/ScottAilliet/wbc-methylation-panel-v2.git
cd wbc-methylation-panel-v2

# 2. Run the installation script
chmod +x install_macos.sh
./install_macos.sh
```

The script will:
1. Install Homebrew (if not present)
2. Install system tools: samtools, bowtie2, tabix, bedtools
3. Create a Python virtual environment (`.venv`)
4. Install Python packages: primer3-py, openpyxl, pandas, numpy, reportlab, pysam
5. Install wgbstools from GitHub source
6. Initialize the hg19 reference genome for wgbstools
7. Verify all installations

This takes about 10–15 minutes depending on your internet speed.

### 1.3 Verify the Installation

```bash
# Activate the environment
source .venv/bin/activate

# Check that the pipeline module loads
python -c "from methyl_panel.pipeline import main; print('Pipeline OK')"

# Check primer3
python -c "import primer3; print(f'primer3-py {primer3.__version__}')"

# Check bowtie2
bowtie2 --version | head -1

# Check samtools
samtools --version | head -1
```

If all commands print a version or "Pipeline OK", you're ready.

---

## Chapter 2 — Data Download

### 2.1 What Data You Need

The pipeline requires:

| Data | Size | Purpose |
|------|------|---------|
| GSE186458 beta files (207 samples) | ~2.4 GB | WGBS methylation data for DMR discovery |
| GSE186458 blocks file | ~1 MB | Genomic block definitions (205 CpG blocks) |
| hg19 genome FASTA | ~3 GB | Reference genome for sequence fetching |
| dbSNP common variants VCF | ~10 GB | SNP screening (Phase 6) |
| Bowtie2 indices | ~4 GB | Specificity screening (Phase 4) |

**Total: ~15–20 GB**

### 2.2 Download Everything

```bash
# Make sure you're in the repo folder and the environment is active
cd ~/wbc-methylation-panel-v2
source .venv/bin/activate

# Download all data
chmod +x download_data.sh
./download_data.sh
```

This downloads:
1. Blocks file from GEO FTP
2. All 207 beta files from GEO (using the manifest in `data/full_atlas_manifest.csv`)
3. hg19 genome via wgbstools
4. dbSNP common variants (MAF ≥ 1%) from NCBI
5. Builds bowtie2 indices for the unconverted genome

**This takes 30–60 minutes** depending on your internet speed. The script skips files that are already downloaded, so you can re-run it if the download is interrupted.

### 2.3 If You Already Have Some Data

If you already have beta files or the blocks file from a previous project, you can copy them into the right locations:

```bash
# Beta files go in data/beta_files/
cp /path/to/your/beta_files/*.beta data/beta_files/

# Blocks file goes in data/
cp /path/to/GSE186458_blocks.s205.bed.gz data/
```

Then run `./download_data.sh` — it will skip what you already have and download the rest.

---

## Chapter 3 — Running the Pipeline

### 3.1 The 9 Pipeline Steps

| Step | Module | What it does |
|------|--------|-------------|
| 1 | `phase1_dmr_loader` | Load DMR blocks from the Excel file |
| 2 | `phase2_bisulfite_convert` | Bisulfite-convert genomic sequences to 6 strands |
| 3 | `phase3_primer3_design` | Design primers with Primer3 (all 161 Primer3Plus parameters) |
| 4 | `phase4_bowtie_specificity` | Screen primers against 6 genome states with bowtie2 |
| 5 | `phase5_structure` | Screen for hairpin secondary structures |
| 6 | `phase6_snp` | Screen for common SNPs (dbSNP, MAF ≥ 1%) |
| 7 | `phase7_dimer` | Predict primer-dimer formation (DimerDetective) |
| 8 | `output_xlsx` | Generate U-assays-style Excel output (27 columns) |
| 9 | `output_pdf` | Generate U-assays-style PDF (one page per primer pair) |

### 3.2 Running All Steps

```bash
source .venv/bin/activate

python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

**Required arguments:**
- `--dmr-xlsx`: Path to the DMR per-CpG Excel file
- `--genome`: Path to the hg19 genome FASTA (can be `.gz`)
- `--min-tm`, `--opt-tm`, `--max-tm`: Tm range in °C — **required, no defaults**

**Optional arguments:**
- `--output-dir`: Where to save results (default: `results/`)
- `--cell-type`: Filter to one cell type (e.g. `MONO`, `BCELL`, `NK`, `GRAN`, `CD3T`, `CD8T`, `CD4T`)
- `--settings`: Load Primer3Plus settings from a file
- `--flank`: Flanking bp around each DMR (default: 100)
- `--max-blocks`: Limit number of blocks (for testing)
- `--bowtie-index-dir`: Directory with bowtie2 indices (for Phase 4)
- `--dbsnp`: Path to dbSNP VCF (for Phase 6)

### 3.3 Running Specific Steps

You can run any combination of steps. This is useful for re-running just the QC steps after changing parameters:

```bash
# Run only primer design (steps 1-3), skip QC
python -m methyl_panel.pipeline --steps 1,2,3 \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62

# Re-run only QC steps on existing primers
python -m methyl_panel.pipeline --steps 5,7 \
    --output-dir results/

# Run only output generation
python -m methyl_panel.pipeline --steps 8,9 \
    --output-dir results/
```

Each step reads from and writes to JSON files in the output directory, so steps can be run independently as long as the prerequisite JSON files exist.

### 3.4 Running a Single Cell Type

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --cell-type CD8T \
    --output-dir results/CD8T/
```

Cell type IDs:

| Cell type | ID |
|-----------|-----|
| Monocytes | `MONO` |
| B cells | `BCELL` |
| NK cells | `NK` |
| Granulocytes | `GRAN` |
| CD3 T cells (pan-T) | `CD3T` |
| CD8 T cells | `CD8T` |
| CD4 T cells | `CD4T` |

### 3.5 Choosing Tm Values

The Tm range is **required** — the pipeline refuses to run without it. This forces an explicit decision about the Tm window for each run.

| Platform | Recommended Tm range |
|----------|---------------------|
| Primer3Plus default | `--min-tm 55 --opt-tm 60 --max-tm 65` |
| Roche Digital LightCycler | `--min-tm 58 --opt-tm 60 --max-tm 62` |
| QIAcuity dPCR | `--min-tm 56 --opt-tm 58 --max-tm 60` |

To use Roche DLC salt conditions, load from a Primer3Plus settings file that has the Roche preset, or modify the salt parameters in `config.py`.

### 3.6 Using a Primer3Plus Settings File

If you have a Primer3Plus settings file (exported from the Primer3Plus web interface):

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --settings Primer3Plus_Settings.txt \
    --output-dir results/
```

The settings file overrides all Primer3 parameters, including Tm. If the settings file has Tm values, you don't need `--min-tm` etc.

---

## Chapter 4 — Understanding the Output

### 4.1 Output Files

After running the pipeline, the output directory contains:

| File | Description |
|------|-------------|
| `dmr_blocks.json` | DMR blocks loaded in Step 1 |
| `converted_sequences.json` | Bisulfite-converted sequences from Step 2 |
| `primers.json` | All primer pairs with QC results (updated by each step) |
| `primer_assays.xlsx` | U-assays-style Excel output (27 columns) |
| `primer_assays.pdf` | U-assays-style PDF (one page per primer pair) |

### 4.2 Excel Output Columns (27 columns)

| # | Column | Description |
|---|--------|-------------|
| 1 | `assay_id` | Unique assay identifier (e.g. MONO_0001_M_001) |
| 2 | `seq_id` | Source DMR block ID |
| 3 | `template_used` | Which bisulfite strand was used (SM, AM, SU, AU) |
| 4 | `assay` | Assay type: M (methylated) or U (unmethylated) |
| 5 | `left_primer_display` | Left primer with CpG notation (y = CpG C) |
| 6 | `right_primer_display` | Right primer with CpG notation |
| 7 | `left_tm_C` | Left primer Tm (°C) |
| 8 | `right_tm_C` | Right primer Tm (°C) |
| 9 | `product_size_bp` | Amplicon size (bp) |
| 10 | `c_total_tail` | Total CpGs in both primer 3' tails |
| 11 | `c_total` | Total CpGs in amplicon |
| 12 | `left_c_total` | CpGs in left primer |
| 13 | `right_c_total` | CpGs in right primer |
| 14 | `left_c_tail` | CpGs in left primer 3' tail (last 5 nt) |
| 15 | `right_c_tail` | CpGs in right primer 3' tail |
| 16 | `sense_meth_mismatch_score` | Mismatches vs sense methylated strand |
| 17 | `sense_unmeth_mismatch_score` | Mismatches vs sense unmethylated strand |
| 18 | `anti_meth_mismatch_score` | Mismatches vs antisense methylated strand |
| 19 | `anti_unmeth_mismatch_score` | Mismatches vs antisense unmethylated strand |
| 20 | `bowtie_passes_filter` | TRUE if unique mapping (Phase 4) |
| 21 | `bowtie_intended_genome` | Which genome state the primer targets |
| 22 | `left_structure_mfe_kcal_mol` | Left primer hairpin MFE (Phase 5) |
| 23 | `right_structure_mfe_kcal_mol` | Right primer hairpin MFE |
| 24 | `primer_dimer_prediction` | Dimer risk tier: high, medium, low (Phase 7) |
| 25 | `primer_dimer_end_min_dg` | Minimum end stability ΔG (kcal/mol) |
| 26 | `common_variant_score` | SNP score: 00=clean, 10=body SNP, 01=3' SNP, 11=both |
| 27 | `mapping_error_note` | Notes on bowtie mapping results |

### 4.3 PDF Output

Each page of the PDF contains:
- **Assay header:** assay_id, seq_id, template, cell type
- **6-strand sequence visualization:** S, SM, SU, A, AM, AU with CpG sites highlighted
- **Primer information table:** sequences, Tm, length, GC%, CpG counts
- **Assay characteristics:** product size, CpG counts, mismatch scores
- **QC filter results:** bowtie, structure MFE, dimer prediction, SNP score
- **DimerDetective footer:** thermodynamic parameters and risk tier explanation

### 4.4 Interpreting QC Results

**Bowtie specificity (Phase 4):**
- `bowtie_passes_filter = TRUE`: Primer aligns uniquely to the intended genome state with no off-target hits
- `bowtie_passes_filter = FALSE`: Off-target alignments found — primer may amplify unintended regions
- `NULL`: Bowtie indices not available — screening was skipped

**Secondary structure (Phase 5):**
- MFE ≥ −1.5 kcal/mol: PASS (no stable hairpin)
- MFE < −1.5 kcal/mol: FAIL (stable hairpin may form, reducing PCR efficiency)
- More negative = more stable = worse

**Primer-dimer (Phase 7, DimerDetective):**
- `low`: end_min_dg > −0.18 kcal/mol — no dimer risk
- `medium`: −2.48 < end_min_dg ≤ −0.18 — mixed zone, may form dimers
- `high`: end_min_dg ≤ −2.48 kcal/mol — high dimer risk
- Passes filter if end_min_dg ≥ −1.0 kcal/mol (conservative cutoff)

> **Bisulfite caveat:** DimerDetective was validated on 40–60% GC primers. Bisulfite MSP primers typically have 10–25% GC, outside the validation range. Treat dimer predictions as indicative, not definitive.

**Common SNPs (Phase 6):**
- Score `00`: No common SNPs in primer — clean
- Score `10`: SNP in primer body (not at 3' end) — moderate risk
- Score `01`: SNP in 3' end (last 5 nt) — high risk, may disrupt primer binding
- Score `11`: SNP in both body and 3' end — highest risk
- `NULL`: dbSNP VCF not available — screening was skipped

---

## Chapter 5 — The 6 Bisulfite Strands

The pipeline converts each genomic region into 6 strands for primer design:

| Strand | Name | Description |
|--------|------|-------------|
| S | Sense | Original top strand (unconverted) |
| SM | Sense Methylated | Top strand, C→T except CpG C's preserved |
| SU | Sense Unmethylated | Top strand, all C→T |
| A | Antisense | Original bottom strand (reverse complement) |
| AM | Antisense Methylated | Bottom strand, C→T except CpG C's preserved |
| AU | Antisense Unmethylated | Bottom strand, all C→T |

**MSP primer design:**
- **Methylated primers** (M assays): designed from SM or AM — CpG C's are preserved, so the primer only matches methylated DNA
- **Unmethylated primers** (U assays): designed from SU or AU — all C's converted to T, so the primer only matches unmethylated DNA

The `template_used` column in the output tells you which strand each primer was designed from.

---

## Chapter 6 — Configuration

### 6.1 Primer3PlusConfig

All 161 Primer3Plus parameters are configurable. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `primer_min_tm` | None (required) | Minimum Tm (°C) |
| `primer_opt_tm` | None (required) | Optimal Tm (°C) |
| `primer_max_tm` | None (required) | Maximum Tm (°C) |
| `primer_min_size` | 18 | Minimum primer length (bp) |
| `primer_opt_size` | 20 | Optimal primer length (bp) |
| `primer_max_size` | 23 | Maximum primer length (bp) |
| `primer_product_size_range` | [(60, 150)] | Amplicon size range (bp) |
| `primer_min_gc` | 30.0 | Minimum GC% |
| `primer_max_gc` | 80.0 | Maximum GC% |
| `primer_max_poly_x` | 5 | Max consecutive identical bases |
| `primer_gc_clamp` | 1 | Require GC clamp at 3' end |
| `primer_max_end_stability` | 9.0 | Max end stability (kcal/mol) |
| `primer_max_self_any` | 8.0 | Max self-complementarity |
| `primer_max_self_end` | 3.0 | Max self-end complementarity |
| `primer_pair_max_diff_tm` | 1.0 | Max Tm difference between primers |

### 6.2 Salt Presets

| Preset | Na+ (mM) | Mg2+ (mM) | dNTP (mM) | Primer (nM) | Correction |
|--------|---------|----------|----------|------------|------------|
| `primer3plus` | 50 | 3 | 1.2 | 250 | SantaLucia |
| `roche_dlc` | 100 | 4.5 | 1.2 | 250 | Owczarzy 2008 |

Default: `primer3plus`. To use Roche DLC conditions, set `salt_preset = "roche_dlc"` in config or load from a Primer3Plus settings file with Roche parameters.

### 6.3 PipelineConfig

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_cpg` | 5 | Minimum CpGs per DMR block |
| `preferred_min_cpg` | 7 | Preferred minimum CpGs |
| `min_cpg_per_primer` | 2 | Minimum CpGs per primer |
| `min_cpg_pair_total` | 4 | Minimum total CpGs per primer pair |
| `max_primers_per_block` | 10 | Max primer pairs per block |
| `delta_means` | 0.3 | DMR delta_means threshold |
| `only_hypo` | True | Only hypomethylated markers |

---

## Chapter 7 — Troubleshooting

### 7.1 Installation Problems

**Problem:** `brew: command not found`
- **Cause:** Homebrew not installed
- **Solution:** The install script should handle this. If not, install manually from https://brew.sh

**Problem:** `python3.11: command not found`
- **Cause:** Python 3.11 not installed
- **Solution:** `brew install python@3.11`

**Problem:** `ModuleNotFoundError: No module named 'methyl_panel'`
- **Cause:** Virtual environment not activated, or pipeline not installed
- **Solution:** Run `source .venv/bin/activate` from the repo folder

**Problem:** `ModuleNotFoundError: No module named 'primer3'`
- **Cause:** primer3-py not installed in the virtual environment
- **Solution:** `pip install primer3-py==2.3.0`

### 7.2 Data Problems

**Problem:** `FileNotFoundError: data/DMR_percpg_full_atlas_all_cell_types.xlsx`
- **Cause:** The DMR Excel file is not in the `data/` folder
- **Solution:** This file is generated by the DMR discovery step (wgbstools find_markers). If you don't have it, you need to run DMR discovery first, or obtain it from a colleague.

**Problem:** `samtools faidx failed` or `No sequence returned`
- **Cause:** hg19 genome FASTA not found or not indexed
- **Solution:** Run `./download_data.sh` to download and index the genome. Or run `samtools faidx data/hg19/hg19.fa.gz` manually.

**Problem:** Beta file download fails for some samples
- **Cause:** NCBI FTP server may be slow or unavailable
- **Solution:** Re-run `./download_data.sh` — it skips files that are already downloaded

### 7.3 Pipeline Problems

**Problem:** `Tm values are REQUIRED`
- **Cause:** You didn't provide `--min-tm`, `--opt-tm`, `--max-tm`
- **Solution:** Always provide all three Tm values, or use `--settings` with a file that contains them

**Problem:** `No primers found for block X`
- **Cause:** The DMR block doesn't have enough CpGs, or the sequence is too short, or no primers passed Primer3's filters
- **Solution:** Try widening the Tm range, or increasing `--flank` (default 100), or lowering `--min-cpg`

**Problem:** Bowtie screening shows `NULL` for all primers
- **Cause:** Bowtie2 indices not found in the specified directory
- **Solution:** Provide `--bowtie-index-dir` pointing to the directory with indices, or run `./download_data.sh` to build them

**Problem:** SNP screening shows `NULL` for all primers
- **Cause:** dbSNP VCF not provided
- **Solution:** Provide `--dbsnp data/dbsnp_common.vcf.gz` (downloaded by `download_data.sh`)

**Problem:** Pipeline is slow
- **Cause:** Processing all 207 samples × 300 markers × 4 templates is computationally intensive
- **Solution:** Use `--cell-type` to process one cell type at a time, or `--max-blocks 10` for testing

---

## Chapter 8 — Repository Structure

```
wbc-methylation-panel-v2/
├── methyl_panel/                    # Python package
│   ├── __init__.py                  # Package init (version 2.0.0)
│   ├── config.py                    # Primer3PlusConfig (161 params) + PipelineConfig
│   ├── pipeline.py                  # CLI entry point (9 steps)
│   ├── phase1_dmr_loader.py         # Step 1: Load DMR blocks from Excel
│   ├── phase2_bisulfite_convert.py  # Step 2: Bisulfite convert to 6 strands
│   ├── phase3_primer3_design.py     # Step 3: Primer3 primer design
│   ├── phase4_bowtie_specificity.py # Step 4: 6-strand bowtie2 screening
│   ├── phase5_structure.py          # Step 5: Hairpin MFE screening
│   ├── phase6_snp.py                # Step 6: dbSNP common variant screening
│   ├── phase7_dimer.py              # Step 7: DimerDetective primer-dimer prediction
│   ├── output_xlsx.py               # Step 8: U-assays-style Excel output
│   └── output_pdf.py                # Step 9: U-assays-style PDF output
├── data/                            # Data files
│   ├── full_atlas_manifest.csv      # 207 sample download manifest
│   ├── full_atlas_groups.csv        # Full atlas cell type groups (82 types)
│   ├── immune_groups.csv            # Immune cell type groups (7 types)
│   └── download_list_nonimmune.csv  # Non-immune sample download list
├── docs/
│   ├── handbook.md                  # This document
│   └── materials_and_methods.md     # Scientific methods reference
├── install_macos.sh                 # macOS installation script
├── download_data.sh                 # Data download script (~15 GB)
├── README.md                        # Project overview
└── .gitignore
```

---

## Chapter 9 — Data Source

Loyfer N et al. (2023). A human methylation atlas for normal and cancer cell types. *Nature* 617, 690–697.

Dataset: [GSE186458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458) — 207 samples, 82 cell types, hg19 (GRCh37).

---

## Appendix A — Glossary

| Term | Definition |
|------|-----------|
| **DMR** | Differentially Methylated Region — a genomic block where methylation differs between cell types |
| **MSP** | Methylation-Specific PCR — PCR using primers designed on bisulfite-converted DNA |
| **CpG** | A cytosine followed by a guanine (5'→3') — the site of DNA methylation in mammals |
| **MFE** | Minimum Free Energy — the stability of a secondary structure (more negative = more stable) |
| **ΔG** | Gibbs free energy — thermodynamic stability (negative = stable, positive = unstable) |
| **MAF** | Minor Allele Frequency — how common a SNP is in the population |
| **Tm** | Melting temperature — the temperature at which half of DNA duplexes dissociate |
| **SM/AM** | Sense/Antisense Methylated — bisulfite-converted strand with CpG C's preserved |
| **SU/AU** | Sense/Antisense Unmethylated — bisulfite-converted strand with all C→T |
| **DimerDetective** | A primer-dimer prediction method based on 3' end thermodynamic stability |
| **U-assays** | The output format (XLSX + PDF) matching the reference assay tracking format |
| **wgbstools** | Tools for analyzing whole-genome bisulfite sequencing data |
| **Beta file** | Binary methylation data file (2 bytes per CpG: methylated count + total count) |

---

## Appendix B — Quick Start (Copy-Paste)

```bash
# 1. Clone and install
git clone https://github.com/ScottAilliet/wbc-methylation-panel-v2.git
cd wbc-methylation-panel-v2
chmod +x install_macos.sh
./install_macos.sh

# 2. Download data (~15 GB, 30-60 min)
chmod +x download_data.sh
./download_data.sh

# 3. Activate environment
source .venv/bin/activate

# 4. Run the pipeline (all steps, all cell types)
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/

# 5. Or run a single cell type
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/DMR_percpg_full_atlas_all_cell_types.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --cell-type CD8T \
    --output-dir results/CD8T/

# 6. Check results
open results/primer_assays.xlsx
open results/primer_assays.pdf
```
