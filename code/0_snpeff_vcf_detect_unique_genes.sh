#!/bin/bash -l

# Start the timer 
script_start=$(date +%s)

# -e/-o pipefail from the very top: any failure from here on stops the
# script with a clear error, instead of continuing on and (as happened
# before) reporting a fake "0 unique genes" success. -u (nounset) is turned
# on further down instead of here, since conda's own activation scripts
# aren't always safe to `source` under nounset.
set -eo pipefail

####################################  
# Defining the working directory
#################################### 

# ROOT=/ROOT/jonathan
# Only activate the local wingen_env conda environment when NOT already running
# inside the Apptainer container. Per the README, the container and wingen_env
# are alternatives, not layered on top of each other - wingen_env (via
# load_conda_env.sh) is only meant for a local, non-containerized run; the
# container already has everything (bcftools included) built in.
#
# APPTAINER_CONTAINER is set by Apptainer itself inside anything it starts -
# confirmed in Apptainer's own docs specifically for `apptainer exec`, which
# is how batch_job_gene_catalog_snakemake.sh runs this script. SINGULARITY_
# CONTAINER is checked too as a free extra safety net - Apptainer's
# Singularity-compatibility layer is documented for *input* variables like
# SINGULARITYENV_*, not confirmed to also mirror this specific one, but
# there's no downside to checking it as well.
if [[ -n "${APPTAINER_CONTAINER:-}" || -n "${SINGULARITY_CONTAINER:-}" ]]; then
    echo "Running inside the Apptainer container - skipping host conda activation (wingen_env is a local-run-only alternative, see README)." >&2
else
    echo "Not running inside the Apptainer container - activating wingen_env for a local run." >&2
    conda_setup_script_path="/proj/uppmax2025-2-125/nobackup/jonathan/anaconda3/etc/profile.d/conda.sh"  # Variable defined in config.sh
    source "$conda_setup_script_path"  # Source Conda initialization script
    # Activate the conda environment
    conda activate wingen_env
fi


export ROOT="/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project" # Variable defined in the main pipeline script!

cd "$ROOT"

# pipeline_script_dir=$ROOT/code/pipeline_scripts

export script_dir="${ROOT}/code" # Variable defined in the main pipeline script!
export pipeline_scripts_dir="${script_dir}" # Variable defined in the main pipeline script!

export results_dir="${ROOT}/results" # Variable defined in the main pipeline script!
export raw_data_dir="${ROOT}/data/raw" # Variable defined in the main pipeline script!

set -u   # nounset from here on - everything below only reads variables this
         # script itself defines, plus "$1".."$5", so this is safe now

####################################
# Options
####################################
# Usage: 0_snpeff_vcf_detect_unique_genes.sh [scope] [format] [input_vcf] [output]
#
#   scope  : all | body | exonic | cds     (default: body)
#     all    - every ANN mention, including upstream_gene_variant / downstream_gene_variant
#              / intergenic_region. NOTE: intergenic_region is not an "umbrella" over
#              upstream/downstream - it's a separate SnpEff category for positions that
#              fall beyond BOTH neighboring genes' up/downstream windows (default 5kb each,
#              SnpEff's -ud option), reported against the flanking gene PAIR
#              ("GeneA-GeneB"), which is exactly the false-match risk the [|]...[|]
#              anchoring in filter_vcf_by_chr_or_gene.sh guards against.
#     body   - the gene region proper: exons + introns, i.e. everywhere BETWEEN a gene's
#              own start and end. Excludes upstream/downstream/intergenic.
#     exonic - exons only: body minus intron_variant. Still includes 5'/3' UTR variants,
#              since UTRs are part of the exon, just not translated.
#     cds    - coding sequence only: missense/synonymous/stop_gained/stop_lost/start_lost/
#              frameshift/inframe_insertion/inframe_deletion/stop_retained/initiator_codon.
#              This is the one that's "just the CDS" - exonic minus the UTRs.
#     (all four are matched against ANN's field 2, the effect/consequence term)
#
#   format : bed | txt                     (default: bed)
#     bed - chrom, start, end, gene_id (0-based, BED-style), always sorted by chromosome
#           then start position (natural order: chr1, chr2, ... chr26, not chr1, chr10,
#           chr11...).
#           Only scope=all or scope=body are allowed here - see the scope caveat above:
#           exonic/cds can't be represented as a single coordinate span without also
#           covering the introns between them, so bed refuses those two rather than
#           silently writing a misleading span. Use format=txt for an exonic/cds gene
#           list, or filter_vcf_by_scope.sh (in 2_wingen_create_gene_catalog.sh) to
#           narrow a specific gene's actual variants to exonic/cds correctly.
#     txt - just the gene IDs, one per line, in the order each gene is first seen in the
#           VCF (deduplicated, not re-sorted - no alphabetical option any more, this is
#           the only order it writes). All four scope values are meaningful here, since a
#           gene list carries no coordinates to mislead - scope just decides which genes
#           qualify for inclusion at all.
#
#   output : output file path, OR a directory to write the default-named file
#            into                          (default: same directory as
#            input_vcf, named unique_genes.<bed|txt>)
#     If this is - or looks like - a directory (it already exists as one, or the path
#     ends in "/"), the file unique_genes.<bed|txt> is written INSIDE it, rather than
#     treating the directory itself as the output file path. That silent mismatch is
#     exactly what broke the first two runs of this script: passing a raw_data_dir-style
#     folder here (out of habit from filter_vcf_by_chr_or_gene.sh's raw_data_dir
#     convention, which this script does NOT use) set OUTPUT to the directory itself, so
#     the final sort/wc redirected straight into it and failed with "... : Is a
#     directory" - and because this script didn't `set -e` at the time, that failure
#     didn't stop it: it went on to print "0 unique genes" and "Runtime: N seconds" as if
#     it had actually succeeded. Both problems are fixed now: this directory-handling,
#     and `set -eo pipefail` above so a failure like it stops the script instead of
#     reporting a fake result.
#
#   input_vcf : path to the SnpEff-annotated input VCF     (default: the full herring
#               cohort VCF hardcoded below)
#     Previously only changeable by hand-editing the hardcoded input_vcf_file line
#     below; now overridable from the command line too, with that same hardcoded path
#     kept as the fallback when this is omitted.
#
# IMPORTANT caveat on the BED output: SnpEff's ANN field does not carry a gene's real
# start/end - only variant positions and (for upstream/downstream/intergenic effects)
# a Distance value. The start/end written here are the leftmost/rightmost variant THIS
# vcf happens to have annotated to that gene at the chosen scope - a lower bound on the
# true gene span, not the span itself. On a subset like head_5000 this can be well short
# of reality, since most of a gene's real length may simply have no variant called in it.
# For the real gene boundaries, use the GTF/GFF3 the SnpEff database was built from.

SCOPE="${1:-body}"
FORMAT="${2:-bed}"

if [[ "$SCOPE" != "all" && "$SCOPE" != "body" && "$SCOPE" != "exonic" && "$SCOPE" != "cds" ]]; then
    echo "ERROR: scope must be 'all', 'body', 'exonic', or 'cds' (got '$SCOPE')" >&2
    exit 1
fi
if [[ "$FORMAT" != "bed" && "$FORMAT" != "txt" ]]; then
    echo "ERROR: format must be 'bed' or 'txt' (got '$FORMAT')" >&2
    exit 1
fi
if [[ "$FORMAT" == "bed" && ( "$SCOPE" == "exonic" || "$SCOPE" == "cds" ) ]]; then
    echo "ERROR: format='bed' only supports scope='all' or 'body' (got scope='$SCOPE')." >&2
    echo "       A single [start,end) span can't represent 'exonic'/'cds' correctly -" >&2
    echo "       exons/CDS are disjoint intervals with introns as gaps between them, so" >&2
    echo "       the span would still include those gaps. Use format='txt' for an" >&2
    echo "       exonic/cds gene list instead, or narrow a specific gene's variants" >&2
    echo "       correctly with filter_vcf_by_scope.sh (called from" >&2
    echo "       2_wingen_create_gene_catalog.sh), which checks each variant's own ANN" >&2
    echo "       rather than relying on a coordinate range." >&2
    exit 1
fi

# input_vcf (3rd)
# input_vcf_file="${raw_data_dir}/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
# input_vcf_file="${raw_data_dir}/head_5000_herring_sentieon_125ind_231031ilter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
input_vcf_file="${3:-${raw_data_dir}/genetic_variation_herring/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz}"

if [[ ! -f "$input_vcf_file" ]]; then
    echo "ERROR: input VCF '$input_vcf_file' does not exist or is not a regular file." >&2
    echo "       (Given as the 3rd argument, or - if omitted - the hardcoded default.)" >&2
    exit 1
fi
# Follows whichever input_vcf_file is actually in play (default or overridden), instead
# of being hardcoded separately and left to drift out of sync with it, as it was before.
vcf_file_dir_gene_catalog="$(dirname "$input_vcf_file")"

# Output gene list file (4th)
OUTPUT="${4:-${vcf_file_dir_gene_catalog}/unique_genes.${FORMAT}}"
# Directory-aware; see the 'output' entry in the usage comment above.
if [[ -d "$OUTPUT" || "$OUTPUT" == */ ]]; then
    OUTPUT="${OUTPUT%/}/unique_genes.${FORMAT}"
fi
mkdir -p "$(dirname "$OUTPUT")"

####################################
# Confirm this really is a SnpEff-annotated VCF
####################################
# Checked from the header (not the filename), so it can't be fooled by a rename, and it
# works the same whether input_vcf_file is plain-text or bgzip-compressed.
if ! bcftools view -h "$input_vcf_file" 2>/dev/null | grep -q '^##INFO=<ID=ANN,'; then
    echo "ERROR: '$input_vcf_file' has no INFO/ANN field (not a SnpEff-annotated VCF)." >&2
    exit 1
fi

####################################
# Explode ANN, classify each mention by scope, write the requested format
####################################
# Extraction goes through 'bcftools query' rather than 'zcat -f | awk' over the raw file.
# Two side effects worth knowing:
#  - No more explicit decompression step: bcftools reads .gz (bgzip or plain gzip) and
#    plain VCF directly through the same code path, confirmed against both.
#  - bcftools has already isolated the ANN value, so the awk below reads it straight from
#    field 3 rather than splitting INFO on ";" and searching for "ANN=" itself. A missing
#    ANN comes through as "." (VCF's missing-value placeholder), checked for accordingly.
# On whether this is actually faster on your ~125-sample file: it should still do less
# per-line work than a generic text pipeline, but I found a real bcftools issue thread
# reporting that bcftools query does slow down with sample count even when FORMAT fields
# aren't requested - so "skips the sample columns entirely" was overstated on my part
# last time. Worth timing directly against the old version on your real file rather than
# trusting my reasoning here. If it's still slow, 'bcftools view -G' (drop genotypes)
# piped in before 'query' is worth trying too - I confirmed the syntax works, but have no
# evidence either way on whether it actually helps versus just moving the cost around.

if [[ "$FORMAT" == "txt" ]]; then

    bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$input_vcf_file" | awk -v scope="$SCOPE" '
        BEGIN { FS="\t" }
        {
            ann = $3
            if (ann == "." || ann == "") next

            # One variant can overlap several transcripts/genes; SnpEff writes one
            # full annotation block per overlap, comma-separated - so each block
            # needs handling separately.
            m = split(ann, blocks, ",")
            for (j = 1; j <= m; j++) {
                split(blocks[j], f, "|")
                gene = f[5]     # field 5 = Gene_ID (stable key - see chat notes on
                                # field 4 Gene_Name vs field 5 Gene_ID)

                # SnpEff joins simultaneous effects at one position with "&", e.g.
                # "missense_variant&splice_region_variant" or
                # "splice_region_variant&intron_variant" - matching the whole effect
                # string against one term misses these entirely, so each "&"-joined
                # part is checked separately against the same term sets instead.
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

                keep = (scope == "all") \
                    || (scope == "body"   && !outside_gene) \
                    || (scope == "exonic" && !outside_gene && !intronic) \
                    || (scope == "cds"    && is_cds)
                if (!keep) continue

                # Dedupe on first sight and print immediately - genes come out in
                # exactly the order each one is first encountered in the VCF, no sorting.
                if (!(gene in seen)) { seen[gene] = 1; print gene }
            }
        }
    ' > "$OUTPUT"

else
    # bed format: track the min/max variant POS seen per gene at the chosen scope.
    # Same classification logic as above, just accumulating into two arrays instead
    # of printing directly, since a gene's span isn't known until every line is seen.
    #
    # Always sorted by chromosome then start position (there's no equivalent of txt's
    # appearance-order here - coordinates are the whole point of the bed output)
    # (-k1,1V uses natural/version order, so chr2 < chr10 the way you'd read it, not the
    # plain lexicographic chr1 < chr10 < chr11 < ... < chr2 that a plain -k1,1 would give
    # on chromosome names like these).
    bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$input_vcf_file" | awk -v scope="$SCOPE" '
        BEGIN { FS="\t"; OFS="\t" }
        {
            chrom = $1; pos = $2; ann = $3
            if (ann == "." || ann == "") next

            m = split(ann, blocks, ",")
            for (j = 1; j <= m; j++) {
                split(blocks[j], f, "|")
                gene = f[5]

                # SnpEff joins simultaneous effects at one position with "&" - see the
                # comment in the txt branch above for why each part is checked separately.
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

                keep = (scope == "all") \
                    || (scope == "body"   && !outside_gene) \
                    || (scope == "exonic" && !outside_gene && !intronic) \
                    || (scope == "cds"    && is_cds)
                if (!keep) continue

                if (!(gene in minpos) || pos < minpos[gene]) minpos[gene] = pos
                if (!(gene in maxpos) || pos > maxpos[gene]) maxpos[gene] = pos
                chrom_of[gene] = chrom
            }
        }
        END {
            for (g in minpos) print chrom_of[g], minpos[g] - 1, maxpos[g], g
        }
    ' | sort -k1,1V -k2,2n > "$OUTPUT"
fi

n_genes=$(wc -l < "$OUTPUT")

# Ending the timer 
script_end=$(date +%s)
# Calculating the script_runtime of the script
script_runtime=$((script_end-script_start))

echo "Gene Catalog Created"
echo "Scope: $SCOPE | Format: $FORMAT | $n_genes unique genes"
echo "The outputfile is stored as: $OUTPUT"
echo "Runtime: $script_runtime seconds"
