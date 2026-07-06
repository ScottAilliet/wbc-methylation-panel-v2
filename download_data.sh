#!/bin/bash
# =============================================================================
# WBC Methylation Panel — Data Download Script
# =============================================================================
# Downloads all data needed to run the pipeline:
#   1. GSE186458 beta files (207 samples, full Loyfer atlas)
#   2. GSE186458 blocks file
#   3. hg19 genome FASTA (for wgbstools)
#   4. dbSNP common variants VCF (for SNP screening)
#   5. Bowtie2 indices (for specificity screening)
#
# Usage:
#   chmod +x download_data.sh
#   ./download_data.sh
#
# Total download size: ~15 GB
# =============================================================================

set -e

DATA_DIR="data"
mkdir -p "$DATA_DIR"

echo "=========================================="
echo "WBC Methylation Panel — Data Download"
echo "=========================================="
echo ""
echo "Data will be downloaded to: $DATA_DIR/"
echo "Total size: ~15 GB"
echo ""

# --- 1. Blocks file ---
echo "=== 1. Downloading blocks file ==="
BLOCKS="$DATA_DIR/GSE186458_blocks.s205.bed.gz"
if [ ! -f "$BLOCKS" ]; then
    wget -q -O "$BLOCKS" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/suppl/GSE186458_blocks.s205.bed.gz"
    tabix -p bed "$BLOCKS"
    echo "✓ Blocks file downloaded"
else
    echo "✓ Blocks file already exists"
fi

# --- 2. Beta files (207 samples) ---
echo ""
echo "=== 2. Downloading beta files (207 samples) ==="
BETA_DIR="$DATA_DIR/beta_files"
mkdir -p "$BETA_DIR"

# Download each beta file from GEO
# The manifest is in data/full_atlas_manifest.csv
MANIFEST="data/full_atlas_manifest.csv"
if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: $MANIFEST not found. Please copy it from the repository."
    exit 1
fi

TOTAL=$(tail -n +2 "$MANIFEST" | wc -l)
COUNT=0
FAILED=0

# Read manifest and download each file
while IFS=, read -r gsm title cell_type fname need_download; do
    [ "$gsm" = "gsm" ] && continue  # Skip header
    COUNT=$((COUNT + 1))

    # Try to find the actual filename (NCBI sometimes truncates names)
    OUTFILE="$BETA_DIR/${gsm}_${fname}.beta"

    if [ -f "$OUTFILE" ] && [ -s "$OUTFILE" ]; then
        continue  # Already downloaded
    fi

    # Try download with GSM prefix
    URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/${gsm}/suppl/${gsm}_${fname}.beta"
    wget -q -O "$OUTFILE" "$URL" 2>/dev/null

    if [ ! -s "$OUTFILE" ]; then
        # Try without GSM prefix (NCBI sometimes uses different naming)
        rm -f "$OUTFILE"
        # Fetch directory listing to find actual filename
        LISTING=$(wget -q -O - "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/${gsm}/suppl/" 2>/dev/null)
        ACTUAL_FNAME=$(echo "$LISTING" | grep -o 'href="[^"]*\.beta"' | grep -v hg38 | head -1 | sed 's/href="//;s/"//')

        if [ -n "$ACTUAL_FNAME" ]; then
            wget -q -O "$BETA_DIR/$ACTUAL_FNAME" \
                "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/${gsm}/suppl/$ACTUAL_FNAME" 2>/dev/null
            if [ -s "$BETA_DIR/$ACTUAL_FNAME" ]; then
                echo "[$COUNT/$TOTAL] OK: $ACTUAL_FNAME"
            else
                echo "[$COUNT/$TOTAL] FAIL: $gsm ($cell_type)"
                rm -f "$BETA_DIR/$ACTUAL_FNAME"
                FAILED=$((FAILED + 1))
            fi
        else
            echo "[$COUNT/$TOTAL] FAIL: $gsm ($cell_type) — no .beta found"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "[$COUNT/$TOTAL] OK: ${gsm}_${fname}.beta"
    fi
done < "$MANIFEST"

echo "✓ Beta files downloaded: $((TOTAL - FAILED))/$TOTAL ($FAILED failed)"

# --- 3. hg19 genome ---
echo ""
echo "=== 3. Downloading hg19 genome ==="
GENOME_DIR="$DATA_DIR/hg19"
mkdir -p "$GENOME_DIR"

if [ ! -f "$GENOME_DIR/hg19.fa.gz" ]; then
    echo "  Initializing hg19 via wgbstools..."
    python3 -c "
import sys
sys.path.insert(0, 'wgbs_tools/src/python')
from init_genome import init_genome
init_genome('hg19')
" 2>/dev/null || echo "  (Run 'wgbstools init_genome hg19' manually if needed)"
    echo "✓ hg19 genome ready"
else
    echo "✓ hg19 genome already exists"
fi

# --- 4. dbSNP common variants ---
echo ""
echo "=== 4. Downloading dbSNP common variants ==="
DBSNP="$DATA_DIR/dbsnp_common.vcf.gz"
if [ ! -f "$DBSNP" ]; then
    echo "  Downloading dbSNP common variants (MAF ≥ 1%)..."
    # Use NCBI dbSNP — this is a large file
    # For hg19: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all.vcf.gz
    wget -q -O "$DBSNP" \
        "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all.vcf.gz"
    tabix -p vcf "$DBSNP"
    echo "✓ dbSNP downloaded"
else
    echo "✓ dbSNP already exists"
fi

# --- 5. Bowtie2 indices ---
echo ""
echo "=== 5. Building Bowtie2 indices ==="
BOWTIE_DIR="$DATA_DIR/bowtie_indices"
mkdir -p "$BOWTIE_DIR"

# Unconverted genome index
if [ ! -f "$BOWTIE_DIR/unconverted.1.bt2" ]; then
    echo "  Building unconverted genome index..."
    # Use the hg19 FASTA from wgbstools
    HG19_FA="$GENOME_DIR/hg19.fa.gz"
    if [ -f "$HG19_FA" ]; then
        bowtie2-build "$HG19_FA" "$BOWTIE_DIR/unconverted"
        echo "✓ Unconverted index built"
    else
        echo "  WARNING: hg19 FASTA not found, skipping bowtie index"
    fi
else
    echo "✓ Unconverted index exists"
fi

# Note: Converted genome indices are built by the pipeline on first run
echo "  (Converted genome indices will be built by the pipeline on first run)"

echo ""
echo "=========================================="
echo "Data download complete!"
echo "=========================================="
echo ""
echo "Data directory: $DATA_DIR/"
echo "  - Beta files: $BETA_DIR/ ($(ls "$BETA_DIR" | wc -l) files)"
echo "  - Blocks: $BLOCKS"
echo "  - Genome: $GENOME_DIR/hg19.fa.gz"
echo "  - dbSNP: $DBSNP"
echo "  - Bowtie indices: $BOWTIE_DIR/"
