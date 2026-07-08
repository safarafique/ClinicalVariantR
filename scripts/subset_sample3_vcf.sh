#!/usr/bin/env bash
# Create a small Sample 3 VCF subset for InterVar pilot runs (GRCh38).
# Usage:
#   bash scripts/subset_sample3_vcf.sh [input_vcf] [output_vcf] [mode] [limit]
# Modes: pass (default), chr (requires 5th arg = chromosome), coding

set -euo pipefail

INPUT="${1:-/mnt/e/ACGM/Sample3.haplotypecaller.filtered_VEP.ann (2).vcf}"
OUTPUT="${2:-/mnt/e/ACGM/results/intervar_compare/Sample3.pilot.vcf}"
MODE="${3:-pass}"
LIMIT="${4:-2000}"
CHR="${5:-chr22}"

mkdir -p "$(dirname "$OUTPUT")"

{
  grep '^##' "$INPUT" || true
  grep '^#CHROM' "$INPUT"
  case "$MODE" in
    pass)
      awk -F'\t' -v limit="$LIMIT" '
        !/^#/ && $7 == "PASS" { print; n++; if (n >= limit) exit }
      ' "$INPUT"
      ;;
    chr)
      awk -F'\t' -v chr="$CHR" -v limit="$LIMIT" '
        !/^#/ && $1 == chr && $7 == "PASS" { print; n++; if (n >= limit) exit }
      ' "$INPUT"
      ;;
    coding)
      awk -F'\t' -v limit="$LIMIT" '
        !/^#/ && $7 == "PASS" && $8 ~ /CSQ=/ && tolower($8) ~ /missense|synonymous|stop_gained|frameshift|splice|inframe|start_lost|stop_lost/ {
          print; n++; if (n >= limit) exit
        }
      ' "$INPUT"
      ;;
    *)
      echo "Unknown mode: $MODE (use pass, chr, or coding)" >&2
      exit 1
      ;;
  esac
} > "$OUTPUT"

COUNT="$(grep -vc '^#' "$OUTPUT" || true)"
echo "Wrote $OUTPUT ($COUNT variants, mode=$MODE, limit=$LIMIT)"
