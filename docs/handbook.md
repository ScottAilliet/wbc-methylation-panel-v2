# WBC Methylation Panel v2 — Operator Handbook

**Version:** Pipeline v2.1.2
**Repository:** `wbc-methylation-panel-v2`
**Last updated:** 2026-07-09

---

## PREFACE

### What This Handbook Is

This handbook is the complete guide for running the WBC Methylation Panel v2 primer design pipeline on your Mac. It covers installation, data download, running the pipeline, and understanding the output — written for a lab scientist who is not a programmer.

### What Changed from v1

This handbook (v2.1.0) replaces the previous version (saved as `handbook_v1.md`). Key changes:

1. **DMR input file:** The pipeline now uses `WBC_Panel_Top200_v7.9.xlsx` (the output of the DMR discovery pipeline) instead of a non-existent `DMR_percpg_full_atlas_all_cell_types.xlsx`. The loader auto-detects the format.
2. **Cell type IDs:** Both short IDs (`MONO`, `BCELL`) and sheet names (`Mono`, `B`) are accepted.
3. **Genome path:** Default changed to `data/hg19/hg19.fa.gz`.
4. **Bug fixes:** `PipelineConfig` missing fields fixed, `download_data.sh` NCBI truncation fixed, `install_macos.sh` tabix→htslib fixed.
5. **Error handling:** Blocks on chromosomes not in the genome are now skipped gracefully instead of crashing.

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
cd Documents
cd ..
cd ~
cd /Users/scottailliet/wbc-methylation-panel-v2
```

**Pro tip:** Type the first few letters of a folder name, then press **Tab** — the terminal auto-completes it. This saves time and prevents typos.

---

### 0.3 Virtual Environments (Instead of Conda)

This pipeline uses a Python **virtual environment** (`.venv`) instead of conda. A virtual environment is an isolated Python setup that keeps this project's packages separate from everything else on your Mac.

**The analogy:** Like a clean bench in the lab — this project's reagents stay on this bench and don't mix with other projects.

**Creating and using it:**

```bash
python3.11 -m venv .venv
source .venv/bin/activate
deactivate
```

**Key point:** You must activate the environment every time you open a new Terminal window and want to use the pipeline. The environment does not persist between Terminal sessions.

---

### 0.4 Installing Software

The `install_macos.sh` script handles everything automatically. It installs:
- **Homebrew** (Mac package manager) — if not already installed
- **System tools:** samtools, bowtie2, htslib (includes tabix), bedtools
- **Python packages:** primer3-py, openpyxl, pandas, numpy, reportlab, pysam
- **wgbstools** (from GitHub source — not on PyPI)

```bash
chmod +x install_macos.sh
./install_macos.sh
```

`chmod +x` means "make this file executable" — you only do this once. `./` means "run the file in this folder."

> **If you get "permission denied":** You forgot `chmod +x`. Run `chmod +x install_macos.sh` first, then `./install_macos.sh`.

> **If you get "command not found: #":** You pasted a `#` comment line from this handbook. Lines starting with `#` are explanations, not commands. Only paste lines that do NOT start with `#`.

---

### 0.5 Running Commands and Scripts

**Running the pipeline:**

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
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
| Pasting `#` comment lines | `zsh: command not found: #` | Only paste lines without `#` |
| Forgot `chmod +x` | "permission denied" | Run `chmod +x` on the script first |
| Stuck at `quote>` prompt | Pasted a line with an apostrophe | Press Ctrl-C to cancel, re-run the command |

---

### 0.8 Quick Reference: The 5 Commands You Will Use Most

| Command | What it does | When to use it |
|---------|----------------|----------------|
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
git clone https://github.com/ScottAilliet/wbc-methylation-panel-v2.git
cd wbc-methylation-panel-v2
chmod +x install_macos.sh
./install_macos.sh
```

> **Do not paste the `#` comment lines** — they are explanations, not commands. Only paste the lines without `#`. If you paste a `#` line, you will see `zsh: command not found: #` which is harmless but confusing.

The script will:
1. Install Homebrew (if not present)
2. Install system tools: samtools, bowtie2, htslib (includes tabix), bedtools
3. Create a Python virtual environment (`.venv`)
4. Install Python packages: primer3-py, openpyxl, pandas, numpy, reportlab, pysam
5. Install wgbstools from GitHub source
6. Initialize the hg19 reference genome for wgbstools
7. Verify all installations

This takes about 10–15 minutes depending on your internet speed.

### 1.3 Verify the Installation

```bash
source .venv/bin/activate

python -c "from methyl_panel.pipeline import main; print('Pipeline OK')"

python -c "import primer3; print(f'primer3-py {primer3.__version__}')"

bowtie2 --version | head -1

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
| dbSNP common variants VCF | ~10 GB | SNP screening (Step 6) |
| Bowtie2 indices | ~4 GB | Specificity screening (Step 4) |

**Total: ~15–20 GB**

### 2.2 Download Everything

```bash
cd ~/wbc-methylation-panel-v2
source .venv/bin/activate

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

> **Note:** The download script uses `wget -q` (quiet mode) for beta files, so there are no progress bars. You will see `[$COUNT/$TOTAL] OK:` lines as each file downloads. The dbSNP download (~10 GB) shows a progress bar.

### 2.3 The DMR Input File

The pipeline reads DMR blocks from an Excel file. The default is:

```
data/WBC_Panel_Top200_v7.9.xlsx
```

This file contains 7 cell-type sheets (CD4T, CD8T, B, NK, Mono, Gran, Blood-T-CD3), each with 105–196 DMR blocks including per-CpG methylation data, genomic coordinates, gene annotations, and cleanliness scores.

**If you don't have this file:** It is the output of the DMR discovery pipeline (repo 1, `wbc-methylation-panel`). You can also generate it by running `wgbstools find_markers` with the downloaded beta files and blocks file. See the [Materials and Methods](materials_and_methods.md) for details.

### 2.4 If You Already Have Some Data

If you already have beta files or the blocks file from a previous project, you can copy them into the right locations:

```bash
cp /path/to/your/beta_files/*.beta data/beta_files/
cp /path/to/GSE186458_blocks.s205.bed.gz data/
```

Then run `./download_data.sh` — it will skip what you already have and download the rest.

---

## Chapter 3 — Running the Pipeline

### 3.1 The 9 Pipeline Steps

| Step | Module | What it does | Needs full genome? |
|------|--------|-------------|-------------------|
| 1 | `phase1_dmr_loader` | Load DMR blocks from the Excel file | No |
| 2 | `phase2_bisulfite_convert` | Bisulfite-convert genomic sequences to 6 strands | Yes (samtools faidx) |
| 3 | `phase3_primer3_design` | Design primers with Primer3 | No (uses step 2 output) |
| 4 | `phase4_bowtie_specificity` | Screen primers against 6 genome states with bowtie2 | Yes (bowtie2 indices) |
| 5 | `phase5_structure` | Screen for hairpin secondary structures | No |
| 6 | `phase6_snp` | Screen for common SNPs (dbSNP, MAF ≥ 1%) | Yes (dbSNP VCF) |
| 7 | `phase7_dimer` | Predict primer-dimer formation (DimerDetective) | No |
| 8 | `output_xlsx` | Generate U-assays-style Excel output (27 columns) | No |
| 9 | `output_pdf` | Generate U-assays-style PDF (one page per primer pair) | No |

**Steps that work without the full data download:** 1, 2, 3, 5, 7, 8, 9 (need only the DMR Excel + hg19 genome)

**Steps that need the full data download:** 4 (bowtie2 indices), 6 (dbSNP VCF)

### 3.2 Running All Steps

```bash
source .venv/bin/activate

python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

**Required arguments:**
- `--dmr-xlsx`: Path to the DMR Excel file (default: `data/WBC_Panel_Top200_v7.9.xlsx`)
- `--genome`: Path to the hg19 genome FASTA (can be `.gz`; default: `data/hg19/hg19.fa.gz`)
- `--min-tm`, `--opt-tm`, `--max-tm`: Tm range in °C — **required, no defaults**

**Optional arguments:**
- `--output-dir`: Where to save results (default: `results/`)
- `--cell-type`: Filter to one cell type (e.g. `MONO`, `BCELL`, `NK`, `GRAN`, `CD3T`, `CD8T`, `CD4T`)
- `--settings`: Load Primer3Plus settings from a file
- `--flank`: Flanking bp around each DMR (default: 100)
- `--max-blocks`: Limit number of blocks (for testing)
- `--min-cpg`: Minimum CpGs per block (default: 5)
- `--bowtie-index-dir`: Directory with bowtie2 indices (for Step 4)
- `--dbsnp`: Path to dbSNP VCF (for Step 6)

### 3.3 Running Without Bowtie2 and dbSNP

If you haven't downloaded the bowtie2 indices or dbSNP VCF yet, you can still run steps 1, 2, 3, 5, 7, 8, 9:

```bash
python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

The bowtie_passes_filter and common_variant_score columns will be NULL in the output. You can run steps 4 and 6 later once the data is downloaded.

### 3.4 Running Specific Steps

You can run any combination of steps. This is useful for re-running just the QC steps after changing parameters:

```bash
python -m methyl_panel.pipeline --steps 1,2,3 \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62

python -m methyl_panel.pipeline --steps 5,7 \
    --output-dir results/

python -m methyl_panel.pipeline --steps 8,9 \
    --output-dir results/
```

Each step reads from and writes to JSON files in the output directory, so steps can be run independently as long as the prerequisite JSON files exist.

### 3.5 Running a Single Cell Type

```bash
python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --cell-type CD8T \
    --output-dir results/CD8T/
```

Cell type IDs:

| Cell type | ID | Sheet name in v7.9 Excel |
|-----------|-----|--------------------------|
| Monocytes | `MONO` | Mono |
| B cells | `BCELL` | B |
| NK cells | `NK` | NK |
| Granulocytes | `GRAN` | Gran |
| CD3 T cells (pan-T) | `CD3T` | Blood-T-CD3 |
| CD8 T cells | `CD8T` | CD8T |
| CD4 T cells | `CD4T` | CD4T |

Both the short ID and the sheet name are accepted (e.g. `--cell-type MONO` and `--cell-type Mono` both work).

### 3.6 Running All 7 Cell Types

To run all 7 immune cell types, run the pipeline 7 times with a different `--cell-type` each time:

```bash
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
        --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --cell-type $CT \
        --output-dir results/$CT/
done
```

Each cell type's results are saved in its own subdirectory (`results/MONO/`, `results/BCELL/`, etc.).

### 3.7 Choosing Tm Values

The Tm range is **required** — the pipeline refuses to run without it. This forces an explicit decision about the Tm window for each run.

| Platform | Recommended Tm range |
|----------|---------------------|
| Primer3Plus default | `--min-tm 55 --opt-tm 60 --max-tm 65` |
| Roche Digital LightCycler | `--min-tm 58 --opt-tm 60 --max-tm 62` |
| QIAcuity dPCR | `--min-tm 56 --opt-tm 58 --max-tm 60` |

To use Roche DLC salt conditions, load from a Primer3Plus settings file that has the Roche preset, or modify the salt parameters in `config.py`.

### 3.8 Using a Primer3Plus Settings File

If you have a Primer3Plus settings file (exported from the Primer3Plus web interface):

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --settings Primer3Plus_Settings.txt \
    --output-dir results/
```

The settings file overrides all Primer3 parameters, including Tm. If the settings file has Tm values, you don't need `--min-tm` etc.

### 3.9 Testing with a Small Subset

To test the pipeline quickly without processing all blocks:

```bash
python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --cell-type MONO \
    --max-blocks 5 \
    --output-dir results/test/
```

`--max-blocks 5` processes only the first 5 blocks that have sequences available in the genome. This is useful for verifying the pipeline works before running the full set.

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
| 1 | `assay_id` | Unique assay identifier (e.g. MONO_0001_001) |
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
| 20 | `bowtie_passes_filter` | TRUE if unique mapping (Step 4) |
| 21 | `bowtie_intended_genome` | Which genome state the primer targets |
| 22 | `left_structure_mfe_kcal_mol` | Left primer hairpin MFE (Step 5) |
| 23 | `right_structure_mfe_kcal_mol` | Right primer hairpin MFE |
| 24 | `primer_dimer_prediction` | Dimer risk tier: high, medium, low (Step 7) |
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

**Bowtie specificity (Step 4):**
- `bowtie_passes_filter = TRUE`: Primer aligns uniquely to the intended genome state with no off-target hits
- `bowtie_passes_filter = FALSE`: Off-target alignments found — primer may amplify unintended regions
- `NULL`: Bowtie indices not available — screening was skipped

**Secondary structure (Step 5):**
- MFE ≥ −1.5 kcal/mol: PASS (no stable hairpin)
- MFE < −1.5 kcal/mol: FAIL (stable hairpin may form, reducing PCR efficiency)
- More negative = more stable = worse

**Primer-dimer (Step 7, DimerDetective):**
- `low`: end_min_dg > −0.18 kcal/mol — no dimer risk
- `medium`: −2.48 < end_min_dg ≤ −0.18 — mixed zone, may form dimers
- `high`: end_min_dg ≤ −2.48 kcal/mol — high dimer risk
- Passes filter if end_min_dg ≥ −1.0 kcal/mol (conservative cutoff)

> **Bisulfite caveat:** DimerDetective was validated on 40–60% GC primers. Bisulfite MSP primers typically have 10–25% GC, outside the validation range. Treat dimer predictions as indicative, not definitive.

**Common SNPs (Step 6):**
- Score `00`: No common SNPs in primer — clean
- Score `10`: SNP in primer body (not at 3' end) — moderate risk
- Score `01`: SNP in 3' end (last 5 nt) — high risk, may disrupt primer binding
- Score `11`: SNP in both body and 3' end — highest risk
- `NULL`: dbSNP VCF not available — screening was skipped

### 4.5 M-assays vs U-assays

The pipeline designs primers on 4 bisulfite-converted templates:
- **Methylated (M assays):** from SM (sense methylated) or AM (antisense methylated) — CpG C's are preserved
- **Unmethylated (U assays):** from SU (sense unmethylated) or AU (antisense unmethylated) — all C's converted to T

For hypomethylated DMRs (the target cell type has low methylation), most primers will be M-assays. This is expected: the methylated templates retain CpG sites that provide the sequence variation needed for Primer3 to design primers within the Tm window. The unmethylated templates become very AT-rich (all C→T), making it harder for Primer3 to find suitable primers.

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
- **Cause:** Virtual environment not activated, or not in the repo folder
- **Solution:** Run `source .venv/bin/activate` from the repo folder

**Problem:** `ModuleNotFoundError: No module named 'primer3'`
- **Cause:** primer3-py not installed in the virtual environment
- **Solution:** `pip install primer3-py==2.3.0`

**Problem:** `brew install tabix` fails — "No available formula with the name 'tabix'"`
- **Cause:** Homebrew removed the standalone tabix formula; it's now bundled with htslib
- **Solution:** Use `brew install htslib` instead. The install script has been fixed to use htslib.

**Problem:** `zsh: permission denied: ./install_macos.sh` or `./download_data.sh`
- **Cause:** The script file is not executable. This happens on first clone, and also after `git checkout --` or `git pull` resets file permissions.
- **Solution:** Run `chmod +x install_macos.sh` (or `chmod +x download_data.sh`), then run the script again.

**Problem:** `git pull` fails with "Your local changes to the following files would be overwritten by merge"
- **Cause:** You modified a file locally (e.g. `chmod +x` changed the file, or you edited it) and the remote version also changed.
- **Solution:** Discard your local copy and take the remote version:
  ```bash
  git checkout -- download_data.sh
  git pull
  chmod +x download_data.sh
  ```
  `git checkout --` replaces your local file with the last committed version. Then `git pull` gets the fix. Then `chmod +x` because `git checkout --` resets permissions.

### 7.2 Data Problems

**Problem:** `FileNotFoundError: data/WBC_Panel_Top200_v7.9.xlsx`
- **Cause:** The DMR Excel file is not in the `data/` folder
- **Solution:** This file is the output of the DMR discovery pipeline (repo 1). Copy it to `data/` or generate it by running `wgbstools find_markers`.

**Problem:** `samtools faidx failed` or `No sequence returned`
- **Cause:** hg19 genome FASTA not found, not indexed, or the chromosome is not in the genome
- **Solution:** Run `./download_data.sh` to download and index the genome. If only some chromosomes are available, blocks on missing chromosomes are automatically skipped.

**Problem:** hg19 genome download fails or `wgbstools init_genome` fails
- **Cause:** wgbstools `init_genome` may fail on some macOS configurations. The download script now downloads hg19 directly from UCSC instead.
- **Solution:** Re-run `./download_data.sh` after `git pull`. If it still fails, download manually:
  ```bash
  wget -O data/hg19/hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  samtools faidx data/hg19/hg19.fa.gz
  ```

**Problem:** dbSNP download fails with 404 Not Found
- **Cause:** NCBI reorganized the dbSNP FTP paths. The file is now `00-common_all.vcf.gz` (with `00-` prefix), not `common_all.vcf.gz`.
- **Solution:** Re-run `./download_data.sh` after `git pull`. The URL has been corrected.

**Problem:** `File truncated at line 1` when indexing genome
- **Cause:** The FASTA file was compressed with gzip instead of bgzip
- **Solution:** `gunzip -c file.fa.gz | bgzip -c > file.fa.bgz && mv file.fa.bgz file.fa.gz && samtools faidx file.fa.gz`

**Problem:** Beta file download fails for some samples
- **Cause:** NCBI FTP server may be slow or unavailable, or filenames are truncated
- **Solution:** Re-run `./download_data.sh` — it skips files that are already downloaded and uses directory listing to find truncated filenames

**Problem:** `zsh: command not found: #`
- **Cause:** You pasted a `#` comment line from this handbook
- **Solution:** Only paste lines that do NOT start with `#`. The `#` lines are explanations.

**Problem:** Terminal shows `quote>` and seems stuck
- **Cause:** You pasted a line containing an apostrophe (like "you're"), and zsh is waiting for a closing quote
- **Solution:** Press Ctrl-C to cancel, then re-run the command without the comment line

### 7.3 Pipeline Problems

**Problem:** `Tm values are REQUIRED`
- **Cause:** You didn't provide `--min-tm`, `--opt-tm`, `--max-tm`
- **Solution:** Always provide all three Tm values, or use `--settings` with a file that contains them

**Problem:** `No primers found for block X`
- **Cause:** The DMR block doesn't have enough CpGs, or the sequence is too short, or no primers passed Primer3's filters
- **Solution:** Try widening the Tm range, or increasing `--flank` (default 100), or lowering `--min-cpg`

**Problem:** `Skipping X (chrN): samtools faidx failed`
- **Cause:** The chromosome is not in your genome FASTA file
- **Solution:** This is expected if you don't have the full hg19 genome. The pipeline skips these blocks and continues with the ones it can fetch. Download the full genome with `./download_data.sh` to process all blocks.

**Problem:** Bowtie screening shows `NULL` for all primers
- **Cause:** Bowtie2 indices not found in the specified directory
- **Solution:** Provide `--bowtie-index-dir` pointing to the directory with indices, or run `./download_data.sh` to build them

**Problem:** SNP screening shows `NULL` for all primers
- **Cause:** dbSNP VCF not provided
- **Solution:** Provide `--dbsnp data/dbsnp_common.vcf.gz` (downloaded by `download_data.sh`)

**Problem:** Pipeline is slow
- **Cause:** Processing all blocks × 4 templates is computationally intensive
- **Solution:** Use `--cell-type` to process one cell type at a time, or `--max-blocks 10` for testing

**Problem:** All primers are M-assays (no U-assays)
- **Cause:** This is expected for hypomethylated DMRs. The unmethylated templates (SU/AU) are very AT-rich after bisulfite conversion, making it difficult for Primer3 to design primers within the Tm window.
- **Solution:** This is normal behavior, not an error. If you need U-assays, try widening the Tm range or increasing `--flank`.

---

## Chapter 8 — Repository Structure

```
wbc-methylation-panel-v2/
├── methyl_panel/                    # Python package
│   ├── __init__.py                  # Package init (version 2.0.0)
│   ├── config.py                    # Primer3PlusConfig (181 params) + PipelineConfig
│   ├── pipeline.py                  # CLI entry point (9 steps)
│   ├── phase1_dmr_loader.py         # Step 1: Load DMR blocks from Excel (v7.9 + Block_Summary)
│   ├── phase2_bisulfite_convert.py  # Step 2: Bisulfite convert to 6 strands
│   ├── phase3_primer3_design.py     # Step 3: Primer3 primer design
│   ├── phase4_bowtie_specificity.py # Step 4: 6-strand bowtie2 screening
│   ├── phase5_structure.py          # Step 5: Hairpin MFE screening
│   ├── phase6_snp.py                # Step 6: dbSNP common variant screening
│   ├── phase7_dimer.py              # Step 7: DimerDetective primer-dimer prediction
│   ├── output_xlsx.py               # Step 8: U-assays-style Excel output
│   └── output_pdf.py                # Step 9: U-assays-style PDF output
├── data/                            # Data files
│   ├── WBC_Panel_Top200_v7.9.xlsx   # DMR blocks (7 cell types, per-CpG methylation)
│   ├── full_atlas_manifest.csv      # 207 sample download manifest
│   ├── full_atlas_groups.csv        # Full atlas cell type groups (82 types)
│   ├── immune_groups.csv            # Immune cell type groups (7 types)
│   └── download_list_nonimmune.csv  # Non-immune sample download list
├── docs/
│   ├── handbook.md                  # This document (v2.1.0)
│   ├── handbook_v1.md               # Previous handbook version (preserved)
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
| **bgzip** | Block-compressed gzip — required by samtools faidx (not regular gzip) |

---

## Appendix B — Quick Start (Copy-Paste)

```bash
git clone https://github.com/ScottAilliet/wbc-methylation-panel-v2.git
cd wbc-methylation-panel-v2
chmod +x install_macos.sh
./install_macos.sh
```

```bash
chmod +x download_data.sh
./download_data.sh
```

```bash
source .venv/bin/activate
```

```bash
python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

```bash
python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --cell-type CD8T \
    --output-dir results/CD8T/
```

```bash
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps 1,2,3,5,7,8,9 \
        --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --cell-type $CT \
        --output-dir results/$CT/
done
```

```bash
open results/primer_assays.xlsx
open results/primer_assays.pdf
```

---

## Appendix C — Change Log

### v2.1.2 (2026-07-09)

- `download_data.sh`: Fixed dbSNP URL — NCBI renamed `common_all.vcf.gz` to `00-common_all.vcf.gz`. Now downloads pre-built `.tbi` tabix index from NCBI instead of building it locally.
- `download_data.sh`: Fixed hg19 genome download — replaced wgbstools `init_genome` (which fails silently on some macOS configs) with direct download from UCSC (`hg19.fa.gz`, ~900 MB). Includes automatic bgzip recompression if needed.
- Handbook: Added troubleshooting for dbSNP 404 and hg19/wgbstools init_genome failure.

### v2.1.1 (2026-07-09)

- `download_data.sh`: Fixed wget consuming stdin from the manifest file, which killed the download loop after the first iteration. The manifest is now read on file descriptor 3 and wget stdin is redirected from /dev/null.
- Handbook: Added troubleshooting for `git pull` conflict ("Your local changes would be overwritten") and `chmod +x` needed after `git checkout --`.

### v2.1.0 (2026-07-09)

**Code fixes:**
- `phase1_dmr_loader.py`: Rewritten to auto-detect and load v7.9 Excel format (one sheet per cell type). Added cell type ID mapping (Mono→MONO, B→BCELL, etc.). Fallback CpG position distribution when `block_cpg_coords` is missing.
- `config.py`: Added missing `PipelineConfig` fields (`max_primers_per_block`, `min_cpg_per_primer`, `min_cpg_pair_total`) that caused step 3 to crash.
- `pipeline.py`: Added graceful skip for blocks on chromosomes not in the genome. Updated default DMR and genome paths. Updated `--cell-type` help text.
- `download_data.sh`: Removed `set -e` (was killing script on first download failure). Use FTP directory listing to handle NCBI filename truncation (74/207 files affected). Added progress counts.
- `install_macos.sh`: Changed `tabix` to `htslib` (Homebrew removed standalone tabix formula).

**Handbook fixes:**
- Corrected DMR file path to `WBC_Panel_Top200_v7.9.xlsx`
- Added cell type ID ↔ sheet name mapping table
- Added "Running Without Bowtie2 and dbSNP" section
- Added "Running All 7 Cell Types" section
- Added "Testing with a Small Subset" section
- Added troubleshooting entries for: permission denied, `#` comment lines, `quote>` stuck prompt, bgzip vs gzip, NCBI filename truncation, missing chromosomes, M-assays only
- Added "What Changed from v1" section
- Added "M-assays vs U-assays" explanation
- Preserved old handbook as `handbook_v1.md`

**Verified by stress testing:** All 7 cell types (MONO, BCELL, NK, GRAN, CD3T, CD8T, CD4T) run end-to-end through steps 1,2,3,5,7,8,9 with 20 blocks each. Total: 123 primer pairs, 7 XLSX files (27 columns each), 7 PDF files.
