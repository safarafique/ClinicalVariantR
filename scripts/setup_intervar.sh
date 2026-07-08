#!/usr/bin/env bash
# Install InterVar from GitHub (WSL/Linux). LF line endings only.
# Usage: bash scripts/setup_intervar.sh [install_dir]

set -euo pipefail

INSTALL_DIR="${1:-$HOME/tools}"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

if [[ ! -d InterVar ]]; then
  echo "Cloning InterVar from https://github.com/WGLab/InterVar ..."
  git clone --depth 1 https://github.com/WGLab/InterVar.git
else
  echo "InterVar folder exists: $INSTALL_DIR/InterVar"
fi

cd InterVar

cat <<'EOF'

=== InterVar installed ===
GitHub:  https://github.com/WGLab/InterVar
Manual:  https://github.com/WGLab/InterVar/blob/master/docs/user-guide/manual.md
Web UI:  https://wintervar.wglab.org/

NEXT STEPS (required):
1. Download ANNOVAR and copy into this folder:
     annotate_variation.pl  table_annovar.pl  convert2annovar.pl
   Register: http://annovar.openbioinformatics.org/annovar_download_form.php

2. Copy mim2gene.txt into intervardb/

3. Edit config.ini:
     buildver = hg38
     database_locat = humandb

4. For Sample 3 (GRCh38, ~655k variants), use pilot subset first:
     bash scripts/subset_sample3_vcf.sh
     bash scripts/run_intervar_sample3.sh

5. Compare with ACMGamp:
     Rscript standalone/compare_acmgamp_intervar.R \
       --acmgamp report.csv \
       --intervar OUT/prefix.hg38_multianno.txt.intervar

FAST PATH (no InterVar install):
     Rscript standalone/compare_acmgamp_intervar.R \
       --acmgamp report.csv \
       --reference ../testig/testig/260324100042.acmg.tsv

EOF
