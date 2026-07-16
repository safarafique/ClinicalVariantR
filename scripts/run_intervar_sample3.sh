#!/usr/bin/env bash
# Run InterVar on Sample 3 (GRCh38). Requires ANNOVAR + humandb in WSL/Linux.
#
# Usage:
#   bash scripts/run_intervar_sample3.sh
#   bash scripts/run_intervar_sample3.sh /path/to/subset.vcf /path/to/out_prefix
#
# Before first run:
#   1. Register and download ANNOVAR: http://annovar.openbioinformatics.org/annovar_download_form.php
#   2. Copy annotate_variation.pl table_annovar.pl convert2annovar.pl into $INTERVAR_DIR
#   3. export ANNOVAR_DIR=~/tools/annovar  (optional, if scripts are not in InterVar folder)
#   4. Place mim2gene.txt in $INTERVAR_DIR/intervardb/ (from OMIM)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INTERVAR_DIR="${INTERVAR_DIR:-$HOME/tools/InterVar}"
ANNOVAR_DIR="${ANNOVAR_DIR:-$INTERVAR_DIR}"
INPUT_VCF="${1:-/mnt/e/ACGM/results/intervar_compare/Sample3.pilot.vcf}"
OUT_PREFIX="${2:-/mnt/e/ACGM/results/intervar_compare/Sample3.pilot}"
BUILDVER="${BUILDVER:-hg38}"

mkdir -p "$(dirname "$OUT_PREFIX")"

for tool in annotate_variation.pl table_annovar.pl convert2annovar.pl; do
  if [[ ! -f "$ANNOVAR_DIR/$tool" ]]; then
    echo "ERROR: Missing ANNOVAR script: $ANNOVAR_DIR/$tool" >&2
    echo "Download ANNOVAR and copy the three .pl files into ANNOVAR_DIR." >&2
    exit 1
  fi
done

if [[ ! -f "$INTERVAR_DIR/Intervar.py" ]]; then
  echo "ERROR: InterVar not found at $INTERVAR_DIR" >&2
  echo "Run: bash scripts/setup_intervar.sh" >&2
  exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "Input VCF not found: $INPUT_VCF" >&2
  echo "Creating pilot subset (2000 PASS variants)..." >&2
  bash "$SCRIPT_DIR/subset_sample3_vcf.sh"
fi

if [[ ! -f "$INTERVAR_DIR/intervardb/mim2gene.txt" ]]; then
  echo "WARNING: intervardb/mim2gene.txt missing (download from OMIM)." >&2
fi

echo "InterVar input : $INPUT_VCF"
echo "InterVar output: ${OUT_PREFIX}.hg38_multianno.txt.intervar"
echo "Build          : $BUILDVER"
echo "This may take a long time on first run (ANNOVAR downloads hg38 humandb)."

cd "$INTERVAR_DIR"

python3 Intervar.py \
  -b "$BUILDVER" \
  -i "$INPUT_VCF" \
  --input_type=VCF \
  -o "$OUT_PREFIX" \
  --table_annovar="$ANNOVAR_DIR/table_annovar.pl" \
  --convert2annovar="$ANNOVAR_DIR/convert2annovar.pl" \
  --annotate_variation="$ANNOVAR_DIR/annotate_variation.pl"

INTERVAR_OUT="${OUT_PREFIX}.${BUILDVER}_multianno.txt.intervar"
if [[ -f "$INTERVAR_OUT" ]]; then
  echo ""
  echo "=== InterVar finished ==="
  echo "Output: $INTERVAR_OUT"
  echo "Rows  : $(($(wc -l < "$INTERVAR_OUT") - 1))"
  echo ""
  echo "Classification summary (column 'InterVar: InterVar and Evidence'):"
  python3 - <<'PY' "$INTERVAR_OUT"
import sys
from collections import Counter

path = sys.argv[1]
with open(path, encoding="utf-8", errors="replace") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    idx = None
    for i, col in enumerate(header):
        if "InterVar" in col and "Evidence" in col:
            idx = i
            break
    if idx is None:
        for i, col in enumerate(header):
            if col.strip() == "InterVar":
                idx = i
                break
    if idx is None:
        print("Could not find InterVar column. Headers:", header[:8], "...")
        sys.exit(0)
    counts = Counter()
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= idx:
            continue
        field = parts[idx]
        cls = field.split("(")[0].strip() if field else "NA"
        counts[cls] += 1
    for cls, n in counts.most_common():
        print(f"  {cls}: {n}")
PY
  echo ""
  echo "Compare with ClinicalVariantR:"
  RSCRIPT="${RSCRIPT:-/mnt/c/Program Files/R/R-4.6.0/bin/Rscript.exe}"
  if [[ -x "$RSCRIPT" ]] || command -v Rscript >/dev/null 2>&1; then
    RS="${RSCRIPT:-Rscript}"
    "$RS" "/mnt/e/ACGM/cml_variant_interpreter/scripts/compare_sample3_acmg_intervar.R" \
      --intervar "$INTERVAR_OUT" \
      --out-dir "$(dirname "$OUT_PREFIX")"
  else
    echo '  Rscript scripts/compare_sample3_acmg_intervar.R --intervar '"$INTERVAR_OUT"
  fi
else
  echo "WARNING: Expected output not found: $INTERVAR_OUT" >&2
  ls -la "${OUT_PREFIX}"* 2>/dev/null || true
  exit 1
fi
