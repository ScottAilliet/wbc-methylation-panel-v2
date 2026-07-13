# WBC Methylation Panel v2 — Operator Handbook

**Version:** Pipeline v2.2.1
**Repository:** `wbc-methylation-panel-v2`
**Last updated:** 2026-07-10

---

## PREFACE

### What This Handbook Is

This handbook is the complete guide for running the WBC Methylation Panel v2 primer design pipeline on your Mac. It covers installation, data download, running the pipeline, and understanding the output — written for a lab scientist who is not a programmer.

### What Changed in v2.2.1

This version restricts the DMR discovery background to **blood cells only**. Previously, find_markers used all 207 atlas samples (including liver, brain, colon, etc.) as background. For a WBC assay, this is biologically incorrect — the DMRs must distinguish the target cell type from other blood cells, not from non-blood tissues. Now only the 36 blood samples (14 blood cell types) are used as background.

Key change:

1. **Blood-only background:** `generate_groups_file()` now filters to `Blood-*` groups only. Each cell type's background is the other blood cell types (e.g. MONO background = B, NK, granulocytes, T cells = 33 samples; CD3T background = B, NK, granulocytes, monocytes = 14 samples).

### What Changed in v2.2.0

This version adds **DMR discovery from scratch** — the pipeline is now fully self-contained. You no longer need a pre-computed DMR Excel file. The new Step 0 uses `wgbstools find_markers` to discover differentially methylated regions directly from the raw 207 beta files downloaded by `download_data.sh`.

Key changes:

1. **New Step 0 (DMR Discovery):** `--discover-dmrs` flag runs `wgbstools find_markers` on raw beta files to find DMRs from scratch. No external DMR file needed.
2. **New module:** `phase0_dmr_discovery.py` — generates groups files, runs find_markers, extracts per-CpG methylation, computes cleanliness scores.
3. **`install_macos.sh` fixed:** wgbstools installed via symlink to `.venv/bin/wgbstools` (not `pip install -e .`, which fails because the package has no `__init__.py`). Added scipy dependency. Removed unnecessary `init_genome` call.
4. **Cell type target mapping:** Each of the 7 cell types maps to specific atlas groups (e.g. CD3T = all T cell subtypes, BCELL = B + B-Memory).
5. **Backward compatible:** Step 1 (load from Excel) still works for users who have a pre-computed DMR file.

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
4. Install Python packages: primer3-py, openpyxl, pandas, numpy, reportlab, pysam, scipy
5. Install wgbstools from GitHub source (compiles C++ tools + pip install)
6. Verify all installations

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

Also verify wgbstools:

```bash
wgbstools --version
```

If `wgbstools: command not found`, the install script should have created a symlink. If not:

```bash
ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
wgbstools --version
```

### 1.2 Updating the Pipeline

When a new version is released (bug fixes, new features), update your local copy:

```bash
cd ~/wbc-methylation-panel-v2
git pull
```

If `git pull` fails with "Your local changes would be overwritten", see the troubleshooting section (Chapter 7). After updating, you may need to re-run `./install_macos.sh` if dependencies changed.

To verify you have the latest version:

```bash
git log --oneline -1
```

Check the [change log](#appendix-c--change-log) at the end of this handbook to see what changed in each version.

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

### 2.3 The DMR Input — Two Options

The pipeline needs DMR blocks (differentially methylated regions) as input. There are two ways to provide them:

**Option A — Discover DMRs from scratch (recommended, fully self-contained):**

Use `--discover-dmrs` to run Step 0, which uses `wgbstools find_markers` to discover DMRs directly from the 207 beta files you downloaded. No external file needed — this is the default workflow for v2.2.0.

**Option B — Load a pre-computed DMR Excel file:**

Use `--dmr-xlsx` to load a pre-computed Excel file (Step 1). This is useful if you already have DMR blocks from a previous analysis. The file must be in v7.9 format (one sheet per cell type) or Block_Summary format.

### 2.4 If You Already Have Some Data

If you already have beta files or the blocks file from a previous project, you can copy them into the right locations:

```bash
cp /path/to/your/beta_files/*.beta data/beta_files/
cp /path/to/GSE186458_blocks.s205.bed.gz data/
```

Then run `./download_data.sh` — it will skip what you already have and download the rest.

---

## Chapter 3 — Running the Pipeline

### 3.1 The Pipeline Steps

| Step | Module | What it does | Needs full genome? |
|------|--------|-------------|-------------------|
| 0 | `phase0_dmr_discovery` | Discover DMRs from raw beta files (wgbstools find_markers) | No |
| 1 | `phase1_dmr_loader` | Load DMR blocks from a pre-computed Excel file | No |
| 2 | `phase2_bisulfite_convert` | Bisulfite-convert genomic sequences to 6 strands | Yes (samtools faidx) |
| 3 | `phase3_primer3_design` | Design primers with Primer3 | No (uses step 2 output) |
| 4 | `phase4_bowtie_specificity` | Screen primers against 6 genome states with bowtie2 | Yes (bowtie2 indices) |
| 5 | `phase5_structure` | Screen for hairpin secondary structures | No |
| 6 | `phase6_snp` | Screen for common SNPs (dbSNP, MAF ≥ 1%) | Yes (dbSNP VCF) |
| 7 | `phase7_dimer` | Predict primer-dimer formation (DimerDetective) | No |
| 8 | `output_xlsx` | Generate U-assays-style Excel output (27 columns) | No |
| 9 | `output_pdf` | Generate U-assays-style PDF (one page per primer pair) | No |

**Step 0 vs Step 1:** You use either Step 0 (`--discover-dmrs`) or Step 1 (default), not both. Step 0 discovers DMRs from raw beta files. Step 1 loads them from a pre-computed Excel file. Both produce the same `dmr_blocks.json` output that steps 2–9 consume.

**Steps that work without bowtie2/dbSNP:** 0, 1, 2, 3, 5, 7, 8, 9 (need only beta files + blocks file + hg19 genome)

**Steps that need the full data download:** 4 (bowtie2 indices), 6 (dbSNP VCF)

### 3.2 Running with DMR Discovery (Recommended)

This is the default workflow for v2.2.0 — fully self-contained, no external DMR file needed:

```bash
source .venv/bin/activate

python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/
```

**Required arguments for `--discover-dmrs`:**
- `--cell-type`: One of `MONO`, `BCELL`, `NK`, `GRAN`, `CD3T`, `CD8T`, `CD4T` (required)
- `--genome`: Path to the hg19 genome FASTA (can be `.gz`; default: `data/hg19/hg19.fa.gz`)
- `--min-tm`, `--opt-tm`, `--max-tm`: Tm range in °C — **required, no defaults**

**DMR discovery arguments (optional):**
- `--beta-dir`: Directory with .beta files (default: `data/beta_files/`)
- `--blocks-file`: Path to blocks BED file (default: `data/GSE186458_blocks.s205.bed.gz`)
- `--groups-csv`: Path to full atlas groups CSV (default: `data/full_atlas_groups.csv`)
- `--threads`: Number of threads for find_markers (default: 2)
- `--top-markers`: Number of top markers per cell type (default: 200)
- `--max-bg-samples`: Max background samples for per-CpG extraction (default: 30)

**What Step 0 does:**
1. Generates a groups CSV file — target samples (your cell type) vs background (other blood cell types only)
2. Generates a beta list file — paths to all .beta files
3. Runs `wgbstools find_markers` — discovers hypomethylated DMRs (delta_means ≥ 0.3, min 3 CpGs, top 200)
4. Parses the output BED file
5. Extracts per-CpG methylation from beta files for each DMR
6. Computes cleanliness scores (target near-zero, background near-one, consistency, coverage)
7. Saves `dmr_blocks.json` — same format as Step 1, so steps 2–9 work unchanged

**Runtime:** ~5 minutes per cell type with 2 threads on a Mac (find_markers ~1 min, per-CpG extraction ~4 min).

### 3.3 Running with a Pre-Computed Excel File

If you already have a DMR Excel file (v7.9 format):

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

### 3.4 Running Without Bowtie2 and dbSNP

If you haven't downloaded the bowtie2 indices or dbSNP VCF yet, you can still run steps 0, 2, 3, 5, 7, 8, 9:

```bash
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/
```

The bowtie_passes_filter and common_variant_score columns will be NULL in the output. You can run steps 4 and 6 later once the data is downloaded.

### 3.5 Running Specific Steps

You can run any combination of steps. This is useful for re-running just the QC steps after changing parameters:

```bash
# DMR discovery only:
python -m methyl_panel.pipeline --steps 0 --discover-dmrs \
    --cell-type MONO --output-dir results/MONO/

# Primer design only (after DMR discovery):
python -m methyl_panel.pipeline --steps 2,3 \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/

# QC steps only:
python -m methyl_panel.pipeline --steps 5,7 \
    --output-dir results/MONO/

# Output only:
python -m methyl_panel.pipeline --steps 8,9 \
    --output-dir results/MONO/
```

Each step reads from and writes to JSON files in the output directory, so steps can be run independently as long as the prerequisite JSON files exist.

### 3.6 Running a Single Cell Type

```bash
python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type CD8T \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/CD8T/
```

Cell type IDs and their target samples in the atlas. Background is always the other blood cell types (36 total blood samples minus target samples):

| Cell type | ID | Atlas groups (target) | # Target | # Background |
|-----------|-----|------------------------|----------|--------------|
| Monocytes | `MONO` | Blood-Monocytes | 3 | 33 |
| B cells | `BCELL` | Blood-B, Blood-B-Mem | 5 | 31 |
| NK cells | `NK` | Blood-NK | 3 | 33 |
| Granulocytes | `GRAN` | Blood-Granulocytes | 3 | 33 |
| CD3 T cells (pan-T) | `CD3T` | All Blood-T-* groups | 22 | 14 |
| CD8 T cells | `CD8T` | Blood-T-CD8, Blood-T-Eff-CD8, Blood-T-EffMem-CD8, Blood-T-Naive-CD8 | 10 | 26 |
| CD4 T cells | `CD4T` | Blood-T-CD4, Blood-T-CenMem-CD4, Blood-T-EffMem-CD4, Blood-T-Naive-CD4 | 10 | 26 |

When using `--discover-dmrs`, the pipeline automatically merges these atlas groups into one target population. The background is restricted to other blood cell types only (36 blood samples total) — not the full atlas — because this is a WBC assay that needs to distinguish the target from other blood cells.

### 3.7 Running All 7 Cell Types

To run all 7 immune cell types with DMR discovery:

```bash
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps all --discover-dmrs \
        --cell-type $CT \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --output-dir results/$CT/
done
```

Each cell type's results are saved in its own subdirectory (`results/MONO/`, `results/BCELL/`, etc.). The find_markers step runs once per cell type (~1 min each), and per-CpG extraction takes ~4 min per cell type. Total runtime: ~35 minutes for all 7.

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
python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --max-blocks 50 \
    --output-dir results/test/
```

`--max-blocks 50` processes only the first 50 blocks that pass the CpG count filter. Note: the first blocks by delta_means may have only 3–4 CpGs (below the `--min-cpg 5` threshold), so use at least `--max-blocks 50` to ensure some blocks pass the filter and produce primers.

To skip re-running find_markers on subsequent runs (reuse the existing BED output):

```bash
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --max-blocks 50 \
    --skip-find-markers \
    --output-dir results/test/
```

---

## Chapter 4 — Understanding the Output

### 4.1 Output Files

After running the pipeline, the output directory contains:

| File | Description |
|------|-------------|
| `groups_MONO.csv` | Groups file for find_markers (target vs background) — Step 0 |
| `beta_list.txt` | List of beta file paths — Step 0 |
| `find_markers_output/` | wgbstools find_markers output (BED files, params) — Step 0 |
| `dmr_blocks.json` | DMR blocks (from Step 0 or Step 1) |
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
- **Cause:** You modified a file locally (e.g. `chmod +x` changed the file, or you edited it) and the remote version also changed. The error message names the conflicting file(s).
- **Solution:** Discard your local copy of the conflicting file(s) and take the remote version. For each file listed in the error message:
  ```bash
  git checkout -- <filename>
  ```
  Then pull and restore permissions:
  ```bash
  git pull
  chmod +x install_macos.sh download_data.sh
  ```
  Example: if the error says `install_macos.sh` would be overwritten:
  ```bash
  git checkout -- install_macos.sh
  git pull
  chmod +x install_macos.sh download_data.sh
  ```
  `git checkout -- <filename>` replaces your local file with the last committed version. Then `git pull` gets the latest changes. Then `chmod +x` because `git checkout --` resets file permissions.

**Problem:** `pip install -e .` fails with `ModuleOrPackageNotFoundError: No file/folder found for package wgbs-tools`
- **Cause:** wgbstools is not a standard Python package — it has no `__init__.py` and its `pyproject.toml` has the `packages` line commented out. Poetry cannot find a package to install.
- **Solution:** This is expected. wgbstools is installed via a symlink, not pip. Run:
  ```bash
  cd wgbs_tools && python3 setup.py && cd ..
  ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
  wgbstools --version
  ```
  The `wgbstools` executable is a symlink to `src/python/wgbs_tools.py`. When Python runs it, it resolves the symlink and adds `src/python/` to `sys.path`, which lets it find `find_markers.py` and other modules. The install script has been fixed to use this approach.

**Problem:** `pipeline.py: error: unrecognized arguments: --discover-dmrs`
- **Cause:** Your local code is out of date. The `--discover-dmrs` flag was added in v2.2.0. You need to `git pull` to get the latest code.
- **Solution:**
  ```bash
  git pull
  git log --oneline -1
  ```
  You should see a recent commit (not `6c47ed1` or earlier). If `git pull` says "already up to date" but `--discover-dmrs` is still not recognized, you may have a **nested directory** problem — the repo was cloned inside a folder with the same name. Check:
  ```bash
  pwd
  git remote -v
  grep -c "discover-dmrs" methyl_panel/pipeline.py
  ```
  If the `grep` returns `0`, you are in the wrong copy of the repo. Look for another `wbc-methylation-panel-v2/` directory (possibly inside the current one) that has the updated code. The correct directory is the one where `grep -c "discover-dmrs" methyl_panel/pipeline.py` returns a number greater than `0`.

### 7.2 Data Problems

**Problem:** `FileNotFoundError: data/WBC_Panel_Top200_v7.9.xlsx`
- **Cause:** You're using Step 1 (load from Excel) but the DMR Excel file is not in the `data/` folder
- **Solution:** Use `--discover-dmrs` instead to discover DMRs from scratch from the beta files. Or, if you have a pre-computed Excel file, copy it to `data/`.

**Problem:** `wgbstools: command not found` (when using --discover-dmrs)
- **Cause:** wgbstools not installed or not on PATH
- **Solution:** Re-run `./install_macos.sh`. If that doesn't fix it, manually install:
  ```bash
  cd wgbs_tools && python3 setup.py && cd ..
  ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
  ```
  Note: `pip install -e .` does NOT work for wgbstools — the package has no `__init__.py` and `pyproject.toml` has the packages line commented out. The symlink approach is the correct installation method.

**Problem:** `wgbstools find_markers failed with return code 1` with `ModuleNotFoundError: No module named 'scipy'`
- **Cause:** scipy is not installed in the venv. This can happen if the wgbstools `pip install -e .` step failed and rolled back other package installations.
- **Solution:** Install scipy directly in the venv:
  ```bash
  source .venv/bin/activate
  pip install scipy
  python -c "import scipy; print(scipy.__version__)"
  ```
  Then re-run the pipeline. The pipeline also checks for scipy before running find_markers and will tell you the exact command to run if it's missing.

**Problem:** `wgbstools find_markers failed with return code 1` (other causes)
- **Cause:** Missing beta files, blocks file, or groups file. Or insufficient memory.
- **Solution:** Ensure `data/beta_files/` has all 207 .beta files (run `./download_data.sh`). Ensure `data/GSE186458_blocks.s205.bed.gz` exists. Check the error message printed by find_markers.

**Problem:** `No markers found` or `0 blocks` from find_markers
- **Cause:** The target cell type has too few samples, or the delta_means threshold is too high, or beta files are missing for the target samples
- **Solution:** Check that the target cell type's beta files are present in `data/beta_files/`. The groups file is automatically filtered to only include samples with available beta files.

**Problem:** `No primers found` with `--max-blocks 10`
- **Cause:** The first 10 blocks by delta_means may have only 3–4 CpGs, below the `--min-cpg 5` threshold
- **Solution:** Use `--max-blocks 50` or more. The first blocks ranked by delta_means tend to have fewer CpGs; blocks with 5+ CpGs appear later in the ranking.

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
│   ├── pipeline.py                  # CLI entry point (steps 0-9)
│   ├── phase0_dmr_discovery.py      # Step 0: DMR discovery via wgbstools find_markers
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
│   ├── full_atlas_manifest.csv      # 207 sample download manifest
│   ├── full_atlas_groups.csv        # Full atlas cell type groups (82 types)
│   ├── immune_groups.csv            # Immune cell type groups (7 types)
│   └── download_list_nonimmune.csv  # Non-immune sample download list
├── docs/
│   ├── handbook.md                  # This document (v2.2.1)
│   ├── handbook_v2.2.md             # Previous handbook version (preserved)
│   ├── handbook_v2.1.md             # Previous handbook version (preserved)
│   ├── handbook_v1.md               # Original handbook version (preserved)
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
# Discover DMRs from scratch and design primers for monocytes
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
    --discover-dmrs --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

```bash
# Same for CD8 T cells
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
    --discover-dmrs --cell-type CD8T \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/CD8T/
```

```bash
# All 7 immune cell types
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
        --discover-dmrs --cell-type $CT \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --output-dir results/$CT/
done
```

```bash
open results/MONO/primer_assays.xlsx
open results/MONO/primer_assays.pdf
```

---

## Appendix C — Change Log

### v2.2.1 (2026-07-10)

- `phase0_dmr_discovery.py`: Background restricted to blood cell types only. Previously, `generate_groups_file()` used all 207 atlas samples as background (including liver, brain, colon, etc.). Now only the 36 blood samples (14 `Blood-*` groups) are used. This is biologically correct for a WBC assay — DMRs must distinguish the target from other blood cells, not from non-blood tissues. Added `BLOOD_GROUPS` constant listing all 14 blood cell type group names.
- Removed multiplex Step 10 (`phase10_multiplex.py`) and all associated CLI arguments (`--multiplex-dirs`, `--cross-dimer-cutoff`, `--tm-tolerance`, `--min-amplicon-diff`).

### v2.2.0 (2026-07-10)

**New feature — DMR discovery from scratch:**
- `phase0_dmr_discovery.py`: New module implementing Step 0. Generates per-cell-type groups files (target vs background), runs `wgbstools find_markers` on raw beta files, parses BED output, extracts per-CpG methylation from beta files, computes 5-component cleanliness scores, and saves `dmr_blocks.json` in the same format as Step 1.
- `pipeline.py`: Added `step0_discover_dmrs()` and new CLI arguments (`--discover-dmrs`, `--beta-dir`, `--blocks-file`, `--groups-csv`, `--wgbstools-path`, `--threads`, `--top-markers`, `--max-bg-samples`, `--skip-find-markers`). When `--discover-dmrs` is set, Step 0 replaces Step 1.
- `install_macos.sh`: Fixed wgbstools installation — uses symlink to `.venv/bin/wgbstools` (not `pip install -e .`, which fails because the package has no `__init__.py`). Added `scipy` dependency. Removed unnecessary `init_genome` call.
- Cell type target mapping: 7 immune cell types mapped to atlas groups (MONO→Blood-Monocytes, BCELL→Blood-B+Blood-B-Mem, NK→Blood-NK, GRAN→Blood-Granulocytes, CD3T→9 T cell subtypes, CD4T→4 CD4 subtypes, CD8T→4 CD8 subtypes).
- Backward compatible: Step 1 (load from Excel) still works via `--dmr-xlsx`.

**Verified by sandbox testing:** Step 0 runs end-to-end on chr22 with 10 beta files (3 Monocytes + 7 background). find_markers finds 1,873 markers, outputs top 200. Per-CpG methylation extracted (target ~0.01, background ~0.93). 197/200 blocks have cleanliness score ≥ 0.6. Full pipeline (steps 0,2,3,5,7,8,9, --max-blocks 50) produces 30 primer pairs, XLSX (30 rows, 27 cols), PDF (76 KB).

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
