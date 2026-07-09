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

# NOTE: No "set -e" — we handle errors per-file so one failed download
# doesn't kill the entire run.

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
    if [ -s "$BLOCKS" ]; then
        tabix -p bed "$BLOCKS" 2>/dev/null || true
        echo "✓ Blocks file downloaded"
    else
        echo "✗ ERROR: Blocks file download failed"
        rm -f "$BLOCKS"
    fi
else
    echo "✓ Blocks file already exists"
fi

# --- 2. Beta files (207 samples) ---
echo ""
echo "=== 2. Downloading beta files (207 samples) ==="
BETA_DIR="$DATA_DIR/beta_files"
mkdir -p "$BETA_DIR"

MANIFEST="data/full_atlas_manifest.csv"
if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: $MANIFEST not found. Please copy it from the repository."
    exit 1
fi

TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
COUNT=0
OK=0
FAILED=0
SKIPPED=0

# Read manifest and download each file.
# NCBI GEO truncates supplementary filenames to ~30 characters, so we
# can't always construct the URL from the manifest. Instead, we fetch
# the directory listing for each GSM and download the first .beta file
# that isn't an hg38 version.
while IFS=, read -r gsm title cell_type fname need_download <&3; do
    [ "$gsm" = "gsm" ] && continue  # Skip header
    COUNT=$((COUNT + 1))

    # Check if any .beta file for this GSM already exists
    EXISTING=$(ls "$BETA_DIR/${gsm}_"*.beta 2>/dev/null | head -1)
    if [ -n "$EXISTING" ] && [ -s "$EXISTING" ]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Fetch directory listing for this GSM
    GSM_PREFIX=$(echo "$gsm" | sed 's/...$/nnn/')
    LISTING=$(wget -q -O - "https://ftp.ncbi.nlm.nih.gov/geo/samples/${GSM_PREFIX}/${gsm}/suppl/" 2>/dev/null </dev/null)

    if [ -z "$LISTING" ]; then
        echo "[$COUNT/$TOTAL] FAIL: $gsm ($cell_type) — could not reach FTP"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Find the .beta file (exclude hg38 versions)
    ACTUAL_FNAME=$(echo "$LISTING" | grep -o 'href="[^"]*\.beta"' | grep -v hg38 | head -1 | sed 's/href="//;s/"//')

    if [ -z "$ACTUAL_FNAME" ]; then
        echo "[$COUNT/$TOTAL] FAIL: $gsm ($cell_type) — no .beta file found"
        FAILED=$((FAILED + 1))
        continue
    fi

    OUTFILE="$BETA_DIR/$ACTUAL_FNAME"
    wget -q -O "$OUTFILE" "https://ftp.ncbi.nlm.nih.gov/geo/samples/${GSM_PREFIX}/${gsm}/suppl/$ACTUAL_FNAME" 2>/dev/null </dev/null

    if [ -s "$OUTFILE" ]; then
        echo "[$COUNT/$TOTAL] OK: $ACTUAL_FNAME"
        OK=$((OK + 1))
    else
        echo "[$COUNT/$TOTAL] FAIL: $gsm ($cell_type) — download failed"
        rm -f "$OUTFILE"
        FAILED=$((FAILED + 1))
    fi

    # Brief pause every 20 files to be gentle on NCBI's FTP server
    if [ $((COUNT % 20)) -eq 0 ]; then
        sleep 1
    fi

done 3< "$MANIFEST"

echo ""
echo "✓ Beta files: $OK downloaded, $SKIPPED already present, $FAILED failed (out of $TOTAL)"

# --- 3. hg19 genome ---
echo ""
echo "=== 3. Downloading hg19 genome ==="
GENOME_DIR="$DATA_DIR/hg19"
mkdir -p "$GENOME_DIR"

if [ ! -f "$GENOME_DIR/hg19.fa.gz" ] || [ ! -f "$GENOME_DIR/hg19.fa.gz.fai" ]; then
    echo "  Downloading hg19 from UCSC (~900 MB compressed)..."
    wget -O "$GENOME_DIR/hg19.fa.gz" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    if [ -s "$GENOME_DIR/hg19.fa.gz" ]; then
        echo "  Indexing genome with samtools faidx..."
        samtools faidx "$GENOME_DIR/hg19.fa.gz" 2>/dev/null || {
            echo "  samtools faidx failed (file may be gzip not bgzip), recompressing..."
            gunzip -c "$GENOME_DIR/hg19.fa.gz" | bgzip -c > "$GENOME_DIR/hg19.fa.bgz"
            mv "$GENOME_DIR/hg19.fa.bgz" "$GENOME_DIR/hg19.fa.gz"
            samtools faidx "$GENOME_DIR/hg19.fa.gz"
        }
        echo "✓ hg19 genome downloaded and indexed"
    else
        echo "✗ ERROR: hg19 download failed"
        rm -f "$GENOME_DIR/hg19.fa.gz"
    fi
else
    echo "✓ hg19 genome already exists"
fi

# --- 4. dbSNP common variants ---
echo ""
echo "=== 4. Downloading dbSNP common variants ==="
DBSNP="$DATA_DIR/dbsnp_common.vcf.gz"
if [ ! -f "$DBSNP" ] || [ ! -f "$DBSNP.tbi" ]; then
    echo "  Downloading dbSNP common variants (MAF ≥ 1%)..."
    echo "  This is a large file (~10 GB) — please be patient..."
    wget -O "$DBSNP" \
        "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz"
    if [ -s "$DBSNP" ]; then
        # Download pre-built tabix index if available
        wget -q -O "$DBSNP.tbi" \
            "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz.tbi"
        if [ ! -s "$DBSNP.tbi" ]; then
            echo "  Building tabix index (this takes a few minutes)..."
            tabix -p vcf "$DBSNP" 2>/dev/null || echo "  (tabix indexing failed — SNP screening will not work)"
        fi
        echo "✓ dbSNP downloaded"
    else
        echo "✗ ERROR: dbSNP download failed"
        rm -f "$DBSNP"
    fi
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

echo "  (Converted genome indices will be built by the pipeline on first run)"

echo ""
echo "=========================================="
echo "Data download complete!"
echo "=========================================="
echo ""
echo "Data directory: $DATA_DIR/"
echo "  - Beta files: $BETA_DIR/ ($(ls "$BETA_DIR" 2>/dev/null | wc -l | tr -d ' ') files)"
echo "  - Blocks: $BLOCKS"
echo "  - Genome: $GENOME_DIR/hg19.fa.gz"
echo "  - dbSNP: $DBSNP"
echo "  - Bowtie indices: $BOWTIE_DIR/"
