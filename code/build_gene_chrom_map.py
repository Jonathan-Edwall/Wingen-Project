#!/usr/bin/env python3
"""
build_gene_chrom_map.py - derive chromosome/gene pairs from a
SnpEff-annotated VCF's INFO/ANN field, once, upfront - so the Snakefile can
read it at parse time instead of any per-gene rule discovering it at runtime.

Assumes the standard SnpEff ANN layout, Gene_ID as the 5th '|'-delimited
subfield (Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|...). If your
existing snpeff_vcf_detect_unique_genes.sh uses Gene_Name (4th field) or
something else instead, change GENE_ID_FIELD_INDEX below to match - the
--check-against option (see below) is exactly how to tell whether the index
is wrong for your actual data, without needing to guess.

Usage:
    python3 build_gene_chrom_map.py INPUT_VCF OUTPUT_TSV \
        [--check-against unique_genes.txt]

Output: a two-column TSV, "chrom<TAB>gene_id", one row per chromosome/gene pair.
Chromosome names retain their spelling in the input VCF. This is
what config["gene_chrom_file"] should point to; the Snakefile reads it directly.
"""
import argparse
import subprocess
import sys
import time
from collections import defaultdict

GENE_ID_FIELD_INDEX = 4
EFFECT_FIELD_INDEX = 1
OUTSIDE_GENE_EFFECTS = {"upstream_gene_variant", "downstream_gene_variant", "intergenic_region",
                    }


def extract_gene_chroms(vcf_path):
    proc = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%INFO/ANN\n", vcf_path],
        capture_output=True, text=True, check=True,
    )
    gene_chroms = defaultdict(set)
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        chrom, ann_field = line.split("\t", 1)
        for entry in ann_field.split(","):
            fields = entry.split("|")
            if len(fields) <= max(GENE_ID_FIELD_INDEX, EFFECT_FIELD_INDEX):
                continue
            effect = fields[EFFECT_FIELD_INDEX].strip()
            if effect in OUTSIDE_GENE_EFFECTS:
                continue  # these carry a "GeneA-GeneB" flanking pair, not a real single gene
            gene_id = fields[GENE_ID_FIELD_INDEX].strip()
            if gene_id:
                gene_chroms[gene_id].add(chrom.strip())
    return gene_chroms

def main():
    script_start = time.perf_counter() # Start timer 
    try:
        ap = argparse.ArgumentParser(description=__doc__)
        ap.add_argument("vcf")
        ap.add_argument("out_tsv")
        ap.add_argument("--check-against", help="existing gene list (one gene_id per line) to cross-check against")
        args = ap.parse_args()

        gene_chroms = extract_gene_chroms(args.vcf)
        if not gene_chroms:
            sys.exit(f"ERROR: no ANN/Gene_ID data extracted from {args.vcf} at all - "
                    f"check the VCF actually has an INFO/ANN field, and that "
                    f"GENE_ID_FIELD_INDEX ({GENE_ID_FIELD_INDEX}) is right for its layout.")

        with open(args.out_tsv, "w") as f:
            for gene_id, chroms in sorted(gene_chroms.items()):
                for chrom in sorted(chroms):
                    f.write(f"{chrom}\t{gene_id}\n")
        n_pairs = sum(len(chroms) for chroms in gene_chroms.values())
        print(f"Wrote {n_pairs} chromosome/gene pairs to {args.out_tsv}", file=sys.stderr)

        if args.check_against:
            with open(args.check_against) as f:
                existing = {line.strip() for line in f if line.strip()}
            extracted = set(gene_chroms)
            missing_from_extraction = existing - extracted
            missing_from_existing_list = extracted - existing
            if missing_from_extraction:
                print(f"\nWARNING: {len(missing_from_extraction)} gene(s) in {args.check_against} had no "
                    f"ANN/Gene_ID match at all - likely means GENE_ID_FIELD_INDEX doesn't match how "
                    f"{args.check_against} was built (e.g. it may use Gene_Name, field index 3, instead):",
                    file=sys.stderr)
                for g in sorted(missing_from_extraction)[:10]:
                    print(f"  {g}", file=sys.stderr)
            if missing_from_existing_list:
                print(f"\nNOTE: {len(missing_from_existing_list)} gene(s) were found via ANN extraction "
                    f"but are not in {args.check_against} (e.g. filtered out there for other reasons):",
                    file=sys.stderr)
                for g in sorted(missing_from_existing_list)[:10]:
                    print(f"  {g}", file=sys.stderr)
            if not missing_from_extraction and not missing_from_existing_list:
                print(f"\nOK: gene set matches {args.check_against} exactly.", file=sys.stderr)
    finally:        
        script_runtime = time.perf_counter() - script_start # End timer 
        print(f"\nRuntime: {script_runtime:.2f} seconds" f"\nOutputfile: {args.out_tsv if 'args' in locals() else 'N/A'}", file=sys.stderr,)

if __name__ == "__main__":
    main()
