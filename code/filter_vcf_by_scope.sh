#!/bin/bash
# Start the timer
script_start=$(date +%s)

# filter_vcf_by_scope.sh
#
# Narrows a VCF down to variants where a specific gene's OWN ANN annotation matches a
# specific scope. Meant to run on a VCF that's already been filtered to one gene (e.g.
# the output of filter_vcf_by_chr_or_gene.sh in "gene" mode) - not because it requires
# that, but because that's the efficient order: filter to one gene first (the full VCF,
# once), then narrow by scope (a much smaller file, per gene).
#
# Why this exists, and why a BED coordinate span (like snpeff_vcf_detect_unique_genes.sh
# produces) can't do this job for "exonic"/"cds": a gene's exons and introns interleave,
# so "just the exons" or "just the CDS" is a set of disjoint intervals with gaps between
# them - but a BED row is one contiguous [start, end) span. Any single span that covers
# a gene's first and last qualifying position will also cover every intron sitting
# between them. Coordinate filtering (bcftools -r/-R) can only ever give you a superset
# for exonic/cds - the only way to get the correct set is to look at each variant's own
# ANN annotation, which is what this script does. (For "body", a coordinate span IS
# correct, since introns are legitimately part of the body - there's no gap to worry
# about. But filter_vcf_by_chr_or_gene.sh's "gene" mode currently matches ANY ANN mention
# of a gene regardless of effect - which is scope=all, not scope=body - so "body" needs
# this script too, at least until gene mode itself becomes scope-aware.)
#
# Usage:
#   ./filter_vcf_by_scope.sh <input.vcf.gz|input.vcf> <gene_id> <scope> <output.vcf> [compress:true|false]
#   VCF_RECORD_COUNT_FILE (optional environment variable): atomically write the
#     successfully decoded output record count, including 0 for an empty VCF.
#     Snakemake uses this status to skip empty genes before analysis.
#
#   scope : all | body | exonic | cds   (see snpeff_vcf_detect_unique_genes.sh for the
#           exact definition of each - this uses the identical classification logic)
#     all    - no filtering at all: whatever came in goes out unchanged. Included for
#              symmetry/scriptability, not because it does anything.
#     body   - keep rows where gene_id's own ANN block is NOT upstream/downstream/intergenic
#     exonic - body, minus rows where gene_id's own block is intron_variant
#     cds    - keep rows where gene_id's own ANN block is a coding-change effect
#              (missense/synonymous/stop_gained/stop_lost/start_lost/frameshift/
#              inframe_insertion/inframe_deletion/stop_retained/initiator_codon)
#
#   Multi-gene rows: if a variant's ANN lists more than one gene (e.g. it's intronic in
#   GeneA and simultaneously exonic in overlapping GeneB), only gene_id's OWN block is
#   checked - a neighboring gene's classification at the same position never leaks in.
#   Tested directly: a row with GeneA=missense_variant and GeneB=intron_variant at the
#   same position is kept for "GeneA cds" and dropped for "GeneB cds", correctly.
#
# Examples:
#   ./filter_vcf_by_scope.sh ENSCHAG00000000004_filtered.vcf.gz ENSCHAG00000000004 cds ENSCHAG00000000004_cds.vcf
#   ./filter_vcf_by_scope.sh ENSCHAG00000000004_filtered.vcf.gz ENSCHAG00000000004 body ENSCHAG00000000004_body.vcf true

set -euo pipefail

if [[ $# -ne 4 && $# -ne 5 ]]; then
    echo "Usage: $0 <input.vcf.gz|input.vcf> <gene_id> <scope> <output.vcf> [compress:true|false]" >&2
    exit 1
fi

INPUT="$1"
GENE="$2"
SCOPE="$3"
OUTPUT="$4"

COMPRESS="${5:-false}"
COMPRESS="${COMPRESS,,}"

if [[ "$SCOPE" != "all" && "$SCOPE" != "body" && "$SCOPE" != "exonic" && "$SCOPE" != "cds" ]]; then
    echo "ERROR: scope must be 'all', 'body', 'exonic', or 'cds' (got '$SCOPE')" >&2
    exit 1
fi
if [[ "$COMPRESS" != "true" && "$COMPRESS" != "false" ]]; then
    echo "ERROR: fifth argument (compress) must be 'true' or 'false' (got '$5')" >&2
    exit 1
fi

####################################
# Deciding compressed vs uncompressed output (same convention as filter_vcf_by_chr_or_gene.sh)
####################################
OUTPUT="${OUTPUT%.gz}"
if [[ "$COMPRESS" == "true" ]]; then
    OUTPUT_TYPE="z"
    OUTPUT="${OUTPUT}.gz"
else
    OUTPUT_TYPE="v"
fi

VCF_RECORD_COUNT_FILE="${VCF_RECORD_COUNT_FILE:-}"
if [[ -n "$VCF_RECORD_COUNT_FILE" ]] && [[
    "$VCF_RECORD_COUNT_FILE" == "$INPUT" || "$VCF_RECORD_COUNT_FILE" == "$OUTPUT" ||
    "$VCF_RECORD_COUNT_FILE" -ef "$INPUT" || "$VCF_RECORD_COUNT_FILE" -ef "$OUTPUT"
]]; then
    echo "ERROR: VCF_RECORD_COUNT_FILE must differ from the input and output VCF paths." >&2
    exit 1
fi

####################################
# Confirm this is a SnpEff-annotated VCF - scope is meaningless without ANN to read
####################################
VCF_HEADER=$(bcftools view -h "$INPUT")
if ! grep -q '^##INFO=<ID=ANN,' <<< "$VCF_HEADER"; then
    echo "ERROR: '$INPUT' has no INFO/ANN field (not a SnpEff-annotated VCF) - scope filtering needs it." >&2
    exit 1
fi

####################################
# Running the filter
####################################
if [[ "$SCOPE" == "all" ]]; then
    # No per-row decision needed - just pass through bcftools for consistent output
    # handling (compression, etc.), same as every other bcftools call in this pipeline.
    bcftools view -O "$OUTPUT_TYPE" -o "$OUTPUT" "$INPUT"

else
    # 'bcftools view' with no filter here is just a decompress-and-normalize pass (reads
    # .gz or plain input identically), piped into awk for the actual per-row decision,
    # then back through bcftools for output formatting/compression. Same classification
    # logic (outside_gene / intronic / is_cds) as snpeff_vcf_detect_unique_genes.sh,
    # just applied as a keep/drop test per row instead of a catalog aggregation.
    bcftools view "$INPUT" | awk -v gene="$GENE" -v scope="$SCOPE" -F'\t' '
        /^#/ { print; next }
        {
            n = split($8, info, ";")
            ann = ""
            for (i = 1; i <= n; i++) if (info[i] ~ /^ANN=/) ann = substr(info[i], 5)
            if (ann == "") next

            keep = 0
            m = split(ann, blocks, ",")
            for (j = 1; j <= m; j++) {
                split(blocks[j], f, "|")
                if (f[5] != gene) continue   # only this gene'\''s own block counts

                # SnpEff joins simultaneous effects at one position with "&", e.g.
                # "missense_variant&splice_region_variant" or
                # "splice_region_variant&intron_variant" - a plain effect=="intron_variant"
                # check misses that second example entirely (false exonic match), and
                # effect=="missense_variant" misses the first (false cds non-match). So
                # each "&"-joined part is checked separately against the same term sets,
                # rather than matching the whole effect string at once.
                outside_gene = 0; intronic = 0; is_cds = 0
                ne = split(f[2], eparts, "&")
                for (ei = 1; ei <= ne; ei++) {
                    ep = eparts[ei]
                    if (ep == "upstream_gene_variant" || ep == "downstream_gene_variant" || \
                        ep == "intergenic_region") outside_gene = 1
                    if (ep == "intron_variant") intronic = 1
                    if (ep == "missense_variant"        || ep == "synonymous_variant" || \
                        ep == "stop_gained"              || ep == "stop_lost"          || \
                        ep == "start_lost"               || ep == "frameshift_variant" || \
                        ep == "inframe_insertion"         || ep == "inframe_deletion"   || \
                        ep == "stop_retained_variant"     || ep == "initiator_codon_variant") is_cds = 1
                }

                if      (scope == "body"   && !outside_gene)              keep = 1
                else if (scope == "exonic" && !outside_gene && !intronic) keep = 1
                else if (scope == "cds"    && is_cds)                     keep = 1
            }
            if (keep) print
        }
    ' | bcftools view -O "$OUTPUT_TYPE" -o "$OUTPUT" -
fi

# A failed decode remains fatal under pipefail, rather than becoming count=0.
n_variants=$(bcftools view -H -G "$OUTPUT" | awk 'END { print NR + 0 }')
if [[ "$n_variants" -eq 0 ]]; then
    echo "WARNING: no VCF records remain for gene '$GENE' after scope '$SCOPE'; this empty gene will be skipped by the catalog workflow." >&2
fi

# Parallel jobs publish separate integer files; no concurrent report appends.
if [[ -n "$VCF_RECORD_COUNT_FILE" ]]; then
    mkdir -p "$(dirname "$VCF_RECORD_COUNT_FILE")"
    record_count_tmp=$(mktemp "${VCF_RECORD_COUNT_FILE}.tmp.XXXXXX")
    trap 'rm -f -- "$record_count_tmp"' EXIT
    printf '%s\n' "$n_variants" > "$record_count_tmp"
    mv -f -- "$record_count_tmp" "$VCF_RECORD_COUNT_FILE"
    trap - EXIT
fi

# Ending the timer
script_end=$(date +%s)
script_runtime=$((script_end-script_start))

echo "Scope filter applied: gene=$GENE scope=$SCOPE -> $n_variants variants"
echo "The outputfile is stored as: $OUTPUT"
echo "Runtime: $script_runtime seconds"
