#!/usr/bin/env bash
# Ubuntu / WSL setup for complete large-VCF analysis
# Usage: bash scripts/ubuntu_setup.sh

set -euo pipefail

echo "=== CML Variant Interpreter — Ubuntu dependencies ==="

sudo apt-get update
sudo apt-get install -y bcftools tabix bgzip samtools

echo ""
echo "bcftools version:"
bcftools --version | head -1

echo ""
echo "=== R CLI packages (data.table, readr, jsonlite) ==="
if command -v Rscript >/dev/null 2>&1; then
  Rscript scripts/install_r_cli_deps.R
else
  echo "Rscript not found. Install R first:"
  echo "  sudo apt-get install -y r-base-core"
  echo "  Rscript scripts/install_r_cli_deps.R"
fi

echo ""
echo "Done. Launch R app from this machine or WSL:"
echo "  R -e \"shiny::runApp('$(pwd)')\""
echo ""
echo "In the app, enable 'Use bcftools (Ubuntu/WSL)' for fastest full-VCF analysis."
