#!/usr/bin/env python3
"""
Prepare an exome VCF (CP-96 and similar) for ClinicalVariantR.

What this script does BEFORE you open ClinicalVariantR:
  1. Read input VCF (SnpEff ANN, hg19, chr-prefixed)
  2. Optionally keep PASS variants only
  3. Optionally keep coding / splice consequence variants only
  4. Join REVEL + CADD_PHRED from local reference files (if provided)
  5. Write a cleaned VCF with REVEL and CADD_PHRED in the INFO column

ClinicalVariantR reads REVEL / CADD from INFO or VEP CSQ — it does NOT compute them.
This script adds scores via lookup tables you download once.

Usage (CP-96 example):
  python prepare_vcf_for_acmg.py \\
    --input "../../Samples folder/CP-96..vcf" \\
    --output "../../results/annotated_vcfs/CP-96/CP-96.prepared.vcf" \\
    --revel ~/vep_refs/revel_grch37.tsv \\
    --cadd ~/vep_refs/whole_genome_SNVs.tsv.gz \\
    --pass-only \\
    --coding-only

Without REVEL/CADD files (cleanup only):
  python prepare_vcf_for_acmg.py \\
    --input "../../Samples folder/CP-96..vcf" \\
    --output "../../results/annotated_vcfs/CP-96/CP-96.prepared.vcf" \\
    --pass-only

Optional: bgzip + tabix after writing (requires bcftools on PATH):
  python prepare_vcf_for_acmg.py ... --bgzip

Dependencies:
  pip install pysam   # recommended for tabix CADD lookup
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional

# Consequences ClinicalVariantR treats as "coding" for rare-variant screens
CODING_CONSEQUENCE_KEYWORDS = (
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "protein_altering_variant",
    "disruptive_inframe_deletion",
    "disruptive_inframe_insertion",
    "conservative_inframe_deletion",
    "conservative_inframe_insertion",
    "initiator_codon_variant",
    "coding_sequence_variant",
)

REVEL_HEADER = '##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL score (dbNSFP/VEP)">'
CADD_HEADER = '##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="CADD Phred score">'


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare CP/AP/BC hg19 VCF for ClinicalVariantR (add REVEL/CADD, filter)."
    )
    p.add_argument(
        "--input", "-i", required=True,
        help="Input VCF path, e.g. '../../Samples folder/CP-96..vcf'",
    )
    p.add_argument(
        "--output", "-o", required=True,
        help="Output prepared VCF path",
    )
    p.add_argument(
        "--revel",
        help="REVEL TSV with columns: chrom,pos,ref,alt,REVEL (hg19/GRCh37)",
    )
    p.add_argument(
        "--cadd",
        help="CADD whole_genome_SNVs.tsv.gz (GRCh37); tabix index .tbi required for speed",
    )
    p.add_argument(
        "--pass-only", action="store_true",
        help="Keep variants with FILTER=PASS only",
    )
    p.add_argument(
        "--coding-only", action="store_true",
        help="Keep variants with coding/splice ANN or CSQ consequence",
    )
    p.add_argument(
        "--bgzip", action="store_true",
        help="Compress output with bgzip and index with tabix (needs bcftools)",
    )
    p.add_argument(
        "--max-variants", type=int, default=0,
        help="Stop after N body lines (0 = no limit; useful for testing)",
    )
    return p.parse_args()


def open_text(path: str):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def variant_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        return chrom if chrom.startswith("chr") else "chr" + chrom[3:]
    return "chr" + chrom


def load_revel_table(path: str) -> Dict[str, float]:
    """Load REVEL scores keyed by chrom:pos:ref:alt."""
    scores: Dict[str, float] = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name.lower(): i for i, name in enumerate(header)}
        for req in ("chrom", "pos", "ref", "alt", "revel"):
            if req not in col:
                raise ValueError(f"REVEL file missing column '{req}'. Header: {header}")

        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = normalize_chrom(parts[col["chrom"]])
            pos = parts[col["pos"]]
            ref = parts[col["ref"]]
            alt = parts[col["alt"]]
            try:
                revel = float(parts[col["revel"]])
            except ValueError:
                continue
            scores[variant_key(chrom, pos, ref, alt)] = revel
            # also store without chr for matching
            scores[variant_key(chrom.replace("chr", ""), pos, ref, alt)] = revel
    print(f"Loaded REVEL: {len(scores):,} keys from {path}", file=sys.stderr)
    return scores


def load_cadd_tabix(path: str) -> Optional[object]:
    """Return pysam TabixFile or None if unavailable."""
    try:
        import pysam  # type: ignore
    except ImportError:
        print(
            "WARN: pysam not installed — CADD lookup skipped.\n"
            "      pip install pysam",
            file=sys.stderr,
        )
        return None

    tbi = path + ".tbi"
    if not os.path.exists(tbi):
        print(
            f"WARN: CADD tabix index missing ({tbi}).\n"
            "      Run: tabix -s 1 -b 2 -e 2 whole_genome_SNVs.tsv.gz",
            file=sys.stderr,
        )
        return None

    return pysam.TabixFile(path)


def lookup_cadd(tabix, chrom: str, pos: str, ref: str, alt: str) -> Optional[float]:
    if tabix is None:
        return None
    import pysam  # type: ignore

    chrom_variants = [chrom, chrom.replace("chr", ""), normalize_chrom(chrom)]
    for c in chrom_variants:
        try:
            rows = tabix.fetch(c, int(pos) - 1, int(pos))
        except (ValueError, OSError):
            continue
        for row in rows:
            parts = row.split("\t")
            if len(parts) < 6:
                continue
            # CADD format: Chrom Pos Ref Alt RawScore Phred
            if parts[1] != str(pos):
                continue
            if parts[2] != ref or parts[3] != alt:
                continue
            try:
                return float(parts[5])
            except (IndexError, ValueError):
                return None
    return None


def parse_info_field(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = ""
    return out


def format_info_field(info: Dict[str, str]) -> str:
    parts = []
    for k, v in info.items():
        if v == "":
            parts.append(k)
        else:
            parts.append(f"{k}={v}")
    return ";".join(parts) if parts else "."


def has_coding_consequence(info: str) -> bool:
    blob = info.upper()
    if "ANN=" in blob or "CSQ=" in blob:
        lower = info.lower()
        return any(kw in lower for kw in CODING_CONSEQUENCE_KEYWORDS)
    return False


def inject_info_scores(
    info: str,
    revel: Optional[float],
    cadd: Optional[float],
) -> str:
    fields = parse_info_field(info)
    if revel is not None:
        fields["REVEL"] = f"{revel:.4f}".rstrip("0").rstrip(".")
    if cadd is not None:
        fields["CADD_PHRED"] = f"{cadd:.2f}".rstrip("0").rstrip(".")
    return format_info_field(fields)


def ensure_info_headers(header_lines: list[str]) -> list[str]:
    has_revel = any("ID=REVEL" in line for line in header_lines)
    has_cadd = any("ID=CADD_PHRED" in line for line in header_lines)
    if has_revel and has_cadd:
        return header_lines

    insert_at = len(header_lines)
    for i, line in enumerate(header_lines):
        if line.startswith("#CHROM"):
            insert_at = i
            break

    new_headers = []
    if not has_revel:
        new_headers.append(REVEL_HEADER)
    if not has_cadd:
        new_headers.append(CADD_HEADER)
    return header_lines[:insert_at] + new_headers + header_lines[insert_at:]


def prepare_vcf(args: argparse.Namespace) -> dict:
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_path}")

    revel_scores = load_revel_table(args.revel) if args.revel else {}
    cadd_tabix = load_cadd_tabix(args.cadd) if args.cadd else None

    stats = {
        "input": str(input_path),
        "output": str(output_path),
        "total": 0,
        "written": 0,
        "skipped_filter": 0,
        "skipped_noncoding": 0,
        "revel_added": 0,
        "cadd_added": 0,
    }

    header: list[str] = []
    with open_text(str(input_path)) as fh_in, open(output_path, "w", encoding="utf-8") as fh_out:
        for line in fh_in:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue

            # write header before first variant
            if stats["total"] == 0:
                header = ensure_info_headers(header)
                for h in header:
                    fh_out.write(h + "\n")

            stats["total"] += 1
            if args.max_variants and stats["written"] >= args.max_variants:
                break

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _id, ref, alt, qual, filt, info = parts[:8]
            rest = parts[8:]

            if args.pass_only and filt != "PASS":
                stats["skipped_filter"] += 1
                continue

            if args.coding_only and not has_coding_consequence(info):
                stats["skipped_noncoding"] += 1
                continue

            key = variant_key(normalize_chrom(chrom), pos, ref, alt)
            revel = revel_scores.get(key)
            cadd = lookup_cadd(cadd_tabix, chrom, pos, ref, alt)

            if revel is not None:
                stats["revel_added"] += 1
            if cadd is not None:
                stats["cadd_added"] += 1

            new_info = inject_info_scores(info, revel, cadd)
            out_parts = [chrom, pos, _id, ref, alt, qual, filt, new_info] + rest
            fh_out.write("\t".join(out_parts) + "\n")
            stats["written"] += 1

            if stats["written"] % 10000 == 0:
                print(f"  ... {stats['written']:,} variants written", file=sys.stderr)

    # edge case: empty VCF (header only)
    if stats["total"] == 0:
        header = ensure_info_headers(header)
        with open(output_path, "w", encoding="utf-8") as fh_out:
            for h in header:
                fh_out.write(h + "\n")

    if args.bgzip:
        bgzip_path(output_path)

    return stats


def bgzip_path(vcf_path: Path) -> None:
    import subprocess

    gz_path = Path(str(vcf_path) + ".gz")
    print(f"Running bgzip -> {gz_path}", file=sys.stderr)
    subprocess.run(["bgzip", "-f", str(vcf_path)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(gz_path)], check=True)
    print(f"Indexed: {gz_path}.tbi", file=sys.stderr)


def print_summary(stats: dict) -> None:
    print("\n=== CP-96 (or sample) VCF preparation complete ===")
    print(f"Input:              {stats['input']}")
    print(f"Output:             {stats['output']}")
    print(f"Variants read:      {stats['total']:,}")
    print(f"Variants written:   {stats['written']:,}")
    print(f"Skipped (filter):   {stats['skipped_filter']:,}")
    print(f"Skipped (noncoding):{stats['skipped_noncoding']:,}")
    print(f"REVEL scores added: {stats['revel_added']:,}")
    print(f"CADD scores added:  {stats['cadd_added']:,}")
    print("\nNext step — run ClinicalVariantR CLI:")
    print("  cd ClinicalVariantR")
    print(f"  Rscript scripts/run_acgm_cli.R \\")
    print(f"    \"{stats['output']}\" \\")
    print(f"    \"../../results/batch/CP-96/cml_panel.csv\" \\")
    print(f"    hematologic_predisposition \\")
    print(f"    \"ABL1,BCR,RUNX1,GATA2,TP53\"")
    print("\nOr open ClinicalVariantR Shiny > Group C > upload prepared VCF.")


def main() -> None:
    args = parse_args()
    stats = prepare_vcf(args)
    print_summary(stats)


if __name__ == "__main__":
    main()
