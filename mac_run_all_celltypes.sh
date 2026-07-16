#!/bin/bash
# ==============================================================================
# WBC Methylation Panel — Full Pipeline Run for All 7 Cell Types
# v2.2.11 (position-based CpG counting + Roche DLC + v2.2.8 loader)
#
# PREREQUISITES:
#   1. git pull origin main  (get v2.2.11, commit e82fbcb)
#   2. data/hg19/hg19.fa + hg19.fa.fai  (download from UCSC, see mac_setup_and_index_build.md)
#   3. data/bowtie2_indices/ with 3 indices  (see mac_setup_and_index_build.md)
#   4. dmr_regions_v2.2.8_260714.xlsx in the repo root or adjust path below
#   5. For relaxed CD4T/CD8T: data/beta_files/, data/GSE186458_blocks.s205.bed.gz,
#      data/full_atlas_groups.csv, wgbstools installed
#
# USAGE:
#   chmod +x mac_run_all_celltypes.sh
#   ./mac_run_all_celltypes.sh           # all cell types, no bowtie2
#   ./mac_run_all_celltypes.sh --bowtie  # include bowtie2 specificity screening
# ==============================================================================

set -e
cd "$(dirname "$0")"

# --- Configuration ---
V228_XLSX="dmr_regions_v2.2.8_260714.xlsx"
ATLAS_XLSX="dmr_lists_DMR_percpg_full_atlas_all_cell_types.xlsx"
GENOME="data/hg19/hg19.fa"
BOWTIE_DIR="data/bowtie2_indices"
MIN_TM=55
OPT_TM=58
MAX_TM=61
PRODUCT_MIN=40
PRODUCT_MAX=150
MIN_CPG=5
# Primer CpG filter: 3 total / 1 per primer (relaxed for sparse-CpG DMRs like CD4/CD8)
MIN_CPG_PAIR_TOTAL=3
MIN_CPG_PER_PRIMER=1

# Check if --bowtie flag was passed
STEPS="1,2,3,5,7,8,9"
BOWTIE_ARG=""
if [ "$1" = "--bowtie" ]; then
    STEPS="1,2,3,4,5,7,8,9"
    BOWTIE_ARG="--bowtie-index-dir $BOWTIE_DIR"
    echo "Bowtie2 screening ENABLED (step 4 included)"
else
    echo "Bowtie2 screening SKIPPED (pass --bowtie to enable)"
fi

COMMON_ARGS="--steps $STEPS \
  --genome $GENOME \
  --min-tm $MIN_TM --opt-tm $OPT_TM --max-tm $MAX_TM \
  --salt-preset roche_dlc \
  --product-size-min $PRODUCT_MIN --product-size-max $PRODUCT_MAX \
  --min-cpg $MIN_CPG \
  --min-cpg-pair-total $MIN_CPG_PAIR_TOTAL \
  --min-cpg-per-primer $MIN_CPG_PER_PRIMER \
  $BOWTIE_ARG"

echo ""
echo "============================================"
echo "WBC Methylation Panel — Full Pipeline Run"
echo "v2.2.11 | Roche DLC | Tm ${MIN_TM}-${MAX_TM}C | Product ${PRODUCT_MIN}-${PRODUCT_MAX}bp"
echo "CpG filter: min_pair_total=${MIN_CPG_PAIR_TOTAL} min_per_primer=${MIN_CPG_PER_PRIMER}"
echo "============================================"
echo ""

# --- Part 1: v2.2.8 blocks (strict 0.70 filter) for all 7 cell types ---
for CT in BCELL CD3T CD4T CD8T GRAN MONO NK; do
    echo ">>> $CT (v2.2.8 blocks)"
    rm -rf results/$CT
    mkdir -p results/$CT
    python3 -m methyl_panel.pipeline \
        $COMMON_ARGS \
        --dmr-xlsx "$V228_XLSX" \
        --cell-type $CT \
        --output-dir results/$CT/ 2>&1 | tail -5
    echo ""
done

# --- Part 2: CD4T/CD8T relaxed (min_bg_subgroup_meth=0.40) ---
# These re-run find_markers from scratch with a relaxed per-subgroup filter
# to recover markers that fail the strict 0.70 filter due to T-cell cross-reactivity.
# Requires: wgbstools, beta files, blocks file, groups CSV
for CT in CD4T CD8T; do
    echo ">>> $CT (relaxed filter, min_bg_subgroup_meth=0.40)"
    rm -rf results/${CT}_relaxed
    mkdir -p results/${CT}_relaxed
    python3 -m methyl_panel.pipeline \
        --steps 0,2,3,4,5,7,8,9 \
        --discover-dmrs \
        --cell-type $CT \
        --beta-dir data/beta_files/ \
        --blocks-file data/GSE186458_blocks.s205.bed.gz \
        --groups-csv data/full_atlas_groups.csv \
        --genome $GENOME \
        $BOWTIE_ARG \
        --min-tm $MIN_TM --opt-tm $OPT_TM --max-tm $MAX_TM \
        --salt-preset roche_dlc \
        --product-size-min $PRODUCT_MIN --product-size-max $PRODUCT_MAX \
        --min-cpg $MIN_CPG \
        --min-cpg-pair-total $MIN_CPG_PAIR_TOTAL \
        --min-cpg-per-primer $MIN_CPG_PER_PRIMER \
        --min-bg-subgroup-meth 0.40 \
        --output-dir results/${CT}_relaxed/ 2>&1 | tail -10
    echo ""
done

# --- Part 3: CD4T/CD8T from full atlas (300 blocks each, unfiltered) ---
# These use the full atlas Excel (Block_Summary format) which has 300 markers
# per cell type discovered with full atlas background (81 groups).
# NOT blood-only filtered — user must check methylation specificity manually.
for CT in CD4T CD8T; do
    echo ">>> $CT (full atlas, 300 blocks)"
    rm -rf results/${CT}_atlas
    mkdir -p results/${CT}_atlas
    python3 -m methyl_panel.pipeline \
        $COMMON_ARGS \
        --dmr-xlsx "$ATLAS_XLSX" \
        --cell-type $CT \
        --output-dir results/${CT}_atlas/ 2>&1 | tail -5
    echo ""
done

# --- Summary ---
echo "============================================"
echo "SUMMARY"
echo "============================================"
for dir in results/*/; do
    ct=$(basename "$dir")
    if [ -f "$dir/primers.json" ]; then
        n=$(python3 -c "import json; print(len(json.load(open('$dir/primers.json'))))")
        echo "  $ct: $n primers"
    fi
done
echo ""
echo "Done! Check results/<cell_type>/primer_assays.xlsx for each cell type."
