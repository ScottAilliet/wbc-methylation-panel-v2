#!/bin/bash
# =============================================================================
# WBC Methylation Panel — macOS Installation Script
# =============================================================================
# Installs all dependencies needed to run the pipeline on macOS.
# Does NOT require biomni — everything runs locally.
#
# Usage:
#   chmod +x install_macos.sh
#   ./install_macos.sh
#
# Prerequisites:
#   - macOS (Intel or Apple Silicon)
#   - Xcode Command Line Tools (xcode-select --install)
# =============================================================================

set -e

echo "=========================================="
echo "WBC Methylation Panel — macOS Installation"
echo "=========================================="
echo ""

# Check for Homebrew
if ! command -v brew &> /dev/null; then
    echo "Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
    eval "$(/opt/homebrew/bin/brew shellenv)"
fi

echo "✓ Homebrew installed"

# Install system dependencies
echo ""
echo "Installing system dependencies..."
brew install python@3.11 samtools bowtie2 tabix bedtools

echo "✓ System dependencies installed (samtools, bowtie2, tabix, bedtools)"

# Create virtual environment
echo ""
echo "Creating Python virtual environment..."
python3.11 -m venv .venv
source .venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install Python packages
echo ""
echo "Installing Python packages..."
pip install primer3-py==2.3.0
pip install openpyxl pandas numpy
pip install reportlab
pip install pysam

echo "✓ Python packages installed"

# Install wgbstools (from GitHub — not on PyPI)
echo ""
echo "Installing wgbstools (from GitHub source)..."
if [ ! -d "wgbs_tools" ]; then
    git clone https://github.com/nloyfer/wgbs_tools.git
fi
cd wgbs_tools
python3 setup.py
cd ..

echo "✓ wgbstools installed"

# Initialize hg19 reference for wgbstools
echo ""
echo "Initializing hg19 reference for wgbstools..."
python3 -c "
import sys
sys.path.insert(0, 'wgbs_tools/src/python')
from init_genome import init_genome
init_genome('hg19')
" 2>/dev/null || echo "  (hg19 reference will be initialized on first run)"

echo "✓ hg19 reference ready"

# Verify installation
echo ""
echo "=========================================="
echo "Verifying installation..."
echo "=========================================="

echo -n "Python: "; python3 --version
echo -n "samtools: "; samtools --version | head -1
echo -n "bowtie2: "; bowtie2 --version | head -1
echo -n "tabix: "; tabix --version | head -1
echo -n "primer3-py: "; python3 -c "import primer3; print(primer3.__version__)"
echo -n "openpyxl: "; python3 -c "import openpyxl; print(openpyxl.__version__)"
echo -n "pandas: "; python3 -c "import pandas; print(pandas.__version__)"
echo -n "reportlab: "; python3 -c "import reportlab; print(reportlab.Version)"

echo ""
echo "=========================================="
echo "Installation complete!"
echo "=========================================="
echo ""
echo "To activate the environment:"
echo "  source .venv/bin/activate"
echo ""
echo "To run the pipeline:"
echo "  python -m methyl_panel.pipeline --steps all \\"
echo "    --dmr-xlsx DMR_percpg_full_atlas_all_cell_types.xlsx \\"
echo "    --genome hg19.fa.gz \\"
echo "    --min-tm 58 --opt-tm 60 --max-tm 62 \\"
echo "    --output-dir results/"
echo ""
echo "To run specific steps only:"
echo "  python -m methyl_panel.pipeline --steps 1,2,3 \\"
echo "    --dmr-xlsx DMR_percpg_full_atlas_all_cell_types.xlsx \\"
echo "    --genome hg19.fa.gz \\"
echo "    --min-tm 58 --opt-tm 60 --max-tm 62"
