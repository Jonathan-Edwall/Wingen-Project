#!/bin/bash
# Start the timer 
script_start=$(date +%s)


# ¤¤¤¤
# Note: input VCF-file needs to be compressed (.gz) and have index file (.tbi or .csi)
#       (In reality the index is only actually required for chr mode - bcftools uses -r
#       for random access, while gene mode streams the file, including when restricted
#       to a chromosome.)
#
# filter_vcf_by_chr_or_gene.sh
# Filters a VCF to one chromosome/region (chr mode), or one gene (gene mode).
#
# Usage:
#   ./filter_vcf_by_chr_or_gene.sh <input.vcf.gz|input.bcf> <chr|gene> <value> <output.vcf> [compress:true|false]
#
#   compress (optional, default: false)
#     true  -> write bgzip-compressed output (-O z)
#     false -> write plain-text output (-O v)
#     Either way, the <output.vcf> filename you pass in gets its extension fixed
#     up automatically to match (.gz appended/stripped as needed), so the file
#     on disk never lies about whether it's actually compressed.
#   FILTER_GENE_CHROM (optional environment variable): limit gene mode to this chromosome.
#   VCF_RECORD_COUNT_FILE (optional environment variable): atomically write the
#     successfully decoded output record count, including 0 for an empty VCF.
#     Snakemake uses this status to skip empty genes before analysis.
#
#   gene mode: <value> is always a single gene ID, e.g. ENSCHAG00000012345. What it's
#   matched against depends on what $INPUT actually is (checked from its header, not the
#   filename, so a rename or symlink can't fool it):
#     - SnpEff-annotated (has INFO/ANN) -> matched against ANN specifically, anchored to
#       pipe boundaries so it can't match as half of an intergenic "GeneA-GeneB" pair.
#     - anything else -> falls back to a plain substring search across the whole INFO
#       column, on the assumption that some other tool may have tagged gene membership
#       there under a different key. This has no known field structure to anchor to, so
#       it's a weaker guarantee than the SnpEff path - see the comments in the gene
#       branch below before relying on it.
#
#   Filtering several genes, or filtering from a BED file of gene coordinates instead of
#   gene-ID text, are both one level up from here, not this script's job:
#     - several genes  -> loop over them, one call per gene (see 2_wingen_create_gene_catalog.sh)
#     - BED coordinates -> resolve the gene's region first (get_gene_region_from_bed() in
#       2_wingen_create_gene_catalog.sh) and call this script in chr mode with that region -
#       that's exactly what 2_wingen_create_gene_catalog.sh already does for
#       vcf_annotation_format="bed".
#
# Examples:
#   ./filter_vcf_by_chr_or_gene.sh herring.vcf.gz chr  chr12                   chr12.vcf
#   ./filter_vcf_by_chr_or_gene.sh herring.vcf.gz chr  chr8:1000000-5000000   inversion.vcf
#   ./filter_vcf_by_chr_or_gene.sh herring.vcf.gz gene ENSCHAG00000012345     gene.vcf
#   ./filter_vcf_by_chr_or_gene.sh herring.vcf.gz chr  chr12                   chr12.vcf   true

set -euo pipefail

# Stop here if the required arguments are missing or there are extra arguments.
if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 <input.vcf.gz|input.bcf> <chr|gene> <value> <output.vcf> [compress:true|false]" >&2
    exit 1
fi

# ROOT="/mnt/c/Users/jonat/pipeline/GitHub/Wingen-Project"
# dataset_folder_name="genetic_variation_herring"
# raw_data_dir="$ROOT/data/raw/$dataset_folder_name"
# input_vcf_file_name="head_5000_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
# INPUT="$raw_data_dir/$input_vcf_file_name"
# preprocessed_dir="$ROOT/data/preprocessed"
# mkdir -p $preprocessed_dir

# FILTER_MODE="gene"
# FILTER_VALUE="ENSCHAG00000000001"
# OUTPUT="$preprocessed_dir/gene_test.vcf"
# bash ./code/filter_vcf_by_chr_or_gene.sh $INPUT $FILTER_MODE $FILTER_VALUE $OUTPUT

# ROOT="/mnt/c/Users/jonat/pipeline/GitHub/Wingen-Project"
# dataset_folder_name="genetic_variation_herring"
# raw_data_dir="$ROOT/data/raw/$dataset_folder_name"
# input_vcf_file_name="head_1488721_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
# INPUT="$raw_data_dir/$input_vcf_file_name"
# preprocessed_dir="$ROOT/data/preprocessed"
# mkdir -p $preprocessed_dir
# FILTER_MODE="chr"
# FILTER_VALUE="chr2"
# OUTPUT="$preprocessed_dir/chr2_test.vcf"
# bash ./code/filter_vcf_by_chr_or_gene.sh $INPUT $FILTER_MODE $FILTER_VALUE $OUTPUT

####################################
# Reading the input arguments
####################################

INPUT="$1" # Input path for the vcf file
echo "$INPUT"
MODE="${2,,}"     # lowercased, so "Chr"/"CHR"/"chr" all work
VALUE="$3"
OUTPUT="$4" # Output path for the filtered vcf file

COMPRESS="${5:-false}"   # optional 5th arg; not given -> "false" (same as the old, only behavior)
COMPRESS="${COMPRESS,,}" # lowercased, so "True"/"TRUE"/"true" all work
GENE_CHROM="${FILTER_GENE_CHROM:-}"

if [[ "$COMPRESS" != "true" && "$COMPRESS" != "false" ]]; then
    echo "ERROR: fifth argument (compress) must be 'true' or 'false' (got '$5')" >&2
    exit 1
fi

####################################
# Deciding compressed vs uncompressed output
####################################
# This is the one place that turns $COMPRESS into bcftools' -O value. Both the
# chr and gene branches below just use $OUTPUT_TYPE, so toggling $COMPRESS
# changes both bcftools calls at once - there's nothing to keep in sync by hand.

OUTPUT="${OUTPUT%.gz}" # strip any .gz the caller already put on, so this is safe to re-run

if [[ "$COMPRESS" == "true" ]]; then
    OUTPUT_TYPE="z"        # bcftools -O z = bgzip-compressed VCF
    OUTPUT="${OUTPUT}.gz"  # keep the filename honest about what's actually inside it
else
    OUTPUT_TYPE="v"        # bcftools -O v = plain-text VCF
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
# Running the filter
####################################

if [[ "$MODE" != "chr" && "$MODE" != "gene" ]]; then
    echo "ERROR: second argument must be 'chr' or 'gene' (got '$2')" >&2
    exit 1
fi

# Read chromosome names from the header instead of scanning the cohort merely
# to inspect its first marker. Exact matches win; only a single chr-prefix
# alias is accepted. A missing chromosome is a valid zero-record selection,
# whereas unreadable input and ambiguous aliases remain errors.
VCF_HEADER=$(bcftools view -h "$INPUT")
VCF_CONTIG_TEXT=$(printf '%s\n' "$VCF_HEADER" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p')
if [[ -z "$VCF_CONTIG_TEXT" ]]; then
    # Older VCFs can omit ##contig definitions. This compatibility fallback
    # reads their records, and propagates query errors rather than hiding them
    # inside process substitution. A valid header-only VCF yields no contigs.
    VCF_CONTIG_TEXT=$(bcftools query -f '%CHROM\n' "$INPUT" | sort -u)
fi
VCF_CONTIGS=()
if [[ -n "$VCF_CONTIG_TEXT" ]]; then
    mapfile -t VCF_CONTIGS <<< "$VCF_CONTIG_TEXT"
fi

resolve_contig() {
    local requested="$1" canonical contig contig_canonical
    local -a aliases=()
    CHROM_MATCH=""
    CONTIG_MISSING=false
    for contig in "${VCF_CONTIGS[@]}"; do
        if [[ "$contig" == "$requested" ]]; then
            CHROM_MATCH="$contig"
            return
        fi
    done
    canonical="$requested"
    [[ "${canonical,,}" != chr* ]] || canonical="${canonical:3}"
    for contig in "${VCF_CONTIGS[@]}"; do
        contig_canonical="$contig"
        [[ "${contig_canonical,,}" != chr* ]] || contig_canonical="${contig_canonical:3}"
        [[ "$contig_canonical" != "$canonical" ]] || aliases+=("$contig")
    done
    if [[ ${#aliases[@]} -gt 1 ]]; then
        echo "ERROR: chromosome '$requested' has ambiguous matches in '$INPUT': ${aliases[*]}" >&2
        exit 1
    elif [[ ${#aliases[@]} -eq 1 ]]; then
        CHROM_MATCH="${aliases[0]}"
    else
        CONTIG_MISSING=true
    fi
}

write_empty_vcf() {
    echo "No matching chromosome '$1' in '$INPUT'; writing a header-only VCF." >&2
    bcftools view -h -O "$OUTPUT_TYPE" -o "$OUTPUT" "$INPUT"
}

if [[ "$MODE" == "chr" ]]; then

    # VALUE can be either:
    #   chr19
    #   19
    #   chr19:20476428-20622999
    #   19:20476428-20622999
    #
    # If VALUE contains ':', separate chromosome and coordinates.
    # If it does not contain ':', treat the whole VALUE as the chromosome.

    if [[ "$VALUE" == *:* ]]; then

        VALUE_CHROM="${VALUE%%:*}"
        VALUE_COORDS="${VALUE#*:}"

        if [[ ! "$VALUE_COORDS" =~ ^[0-9]+(-[0-9]+)?$ ]]; then
            echo "ERROR: invalid region coordinates '$VALUE_COORDS'." >&2
            exit 1
        fi
    else
        VALUE_CHROM="$VALUE"
        VALUE_COORDS=""
    fi
    resolve_contig "$VALUE_CHROM"

    # The compressed VCF needs an index file (.tbi or .csi).
    # Rebuild the index if it is older than the VCF.

    if [[ -f "${INPUT}.csi" ]]; then

        if [[ "${INPUT}" -nt "${INPUT}.csi" ]]; then
            echo "Index is older than VCF. Rebuilding..."
            bcftools index -f -c "$INPUT" || exit 1
        fi

    elif [[ -f "${INPUT}.tbi" ]]; then

        if [[ "${INPUT}" -nt "${INPUT}.tbi" ]]; then
            echo "Index is older than VCF. Rebuilding..."
            bcftools index -f -c "$INPUT" || exit 1
        fi

    else

        echo "No index found. Creating CSI index..."
        bcftools index -f -c "$INPUT" || exit 1

    fi

    # bcftools view = reads the file and writes out a filtered copy.
    #
    # -r "$VALUE"        keep only this chromosome/region
    # -O "$OUTPUT_TYPE"  z = bgzip-compressed VCF,
    #                    v = plain-text VCF
    # -o "$OUTPUT"      path for the saved result
    # "$INPUT"           path for the input file
    if [[ "$CONTIG_MISSING" == true ]]; then
        write_empty_vcf "$VALUE_CHROM"
    else
        VALUE="$CHROM_MATCH${VALUE_COORDS:+:$VALUE_COORDS}"
        echo "Filtering in $MODE with -r (VALUE): $VALUE"
        bcftools view -r "$VALUE" -O "$OUTPUT_TYPE" -o "$OUTPUT" "$INPUT"
    fi

elif [[ "$MODE" == "gene" ]]; then

    # gene mode always takes a single gene ID - filtering several genes, or filtering
    # from a BED file of coordinates, are both handled a level up (see the header
    # comment above and 2_wingen_create_gene_catalog.sh).

    CHROM_ARGS=()
    CONTIG_MISSING=false
    if [[ -n "$GENE_CHROM" ]]; then
        resolve_contig "$GENE_CHROM"
        if [[ "$CONTIG_MISSING" == false ]]; then
            CHROM_ARGS=(-t "$CHROM_MATCH")
        fi
    fi

    # Confirmed against the VCF's own header, whether plain or bgzip-compressed.
    if [[ "$CONTIG_MISSING" == true ]]; then
        write_empty_vcf "$GENE_CHROM"
    elif grep -q '^##INFO=<ID=ANN,' <<< "$VCF_HEADER"; then

        ####################################
        # SnpEff-annotated VCF.
        # bcftools view = reads the file and writes out a filtered copy
        # -i                        keep only rows where this condition is true
        # INFO/ANN ~ "[|]ID[|]"     true if ANN (inside INFO) contains this gene ID as a
        #                           standalone field value - anchored to pipe characters
        #                           on both sides so it can't match as half of an
        #                           intergenic "GeneA-GeneB" neighbor pair, which a plain
        #                           substring search would.
        # -O "$OUTPUT_TYPE"         z = bgzip-compressed VCF, v = plain-text VCF
        # -o "$OUTPUT"              path for the saved result
        # "$INPUT"                  path for the input file
        ####################################
        ANN_GENE_PATTERN="${VALUE//./[.]}"
        bcftools view "${CHROM_ARGS[@]}" -i "INFO/ANN ~ \"[|]${ANN_GENE_PATTERN}[|]\"" -O "$OUTPUT_TYPE" -o "$OUTPUT" "$INPUT"

    else

        ####################################
        # Not SnpEff. bcftools' -i needs a specific INFO/TAG - there is no "search the
        # whole INFO string regardless of tag" expression in its filter language (tried:
        # `INFO ~ "..."` -> bcftools rejects it outright, "tag INFO is not defined").
        # So this falls back to a plain-text substring search over the INFO column
        # itself (awk, not bcftools -i), on the assumption that some OTHER annotation
        # tool may have tagged gene membership somewhere in INFO, just under a different
        # key than SnpEff's ANN.
        #
        # This is a much weaker guarantee than the ANN path above: there's no known
        # field structure to anchor to, so it can false-match (e.g. as a fragment of an
        # unrelated value, or of a different gene's ID that happens to contain this one
        # as a substring) and a zero-variant result is ambiguous - it could mean the
        # gene genuinely has none, or that this VCF has no gene-tagging info in INFO at
        # all. Spot-check a known line's INFO field by eye before trusting this path.
        ####################################
        echo "WARNING: '$INPUT' has no INFO/ANN field (not SnpEff-annotated)." >&2
        echo "         Falling back to a plain, unanchored substring search of the" >&2
        echo "         whole INFO column for '$VALUE'. See script comments for why" >&2
        echo "         this is weaker than the SnpEff path." >&2

        bcftools view "${CHROM_ARGS[@]}" -Ov "$INPUT" \
            | awk -F'\t' -v val="$VALUE" '/^#/ {print; next} index($8, val) > 0 {print}' \
            | bcftools view -O "$OUTPUT_TYPE" -o "$OUTPUT" -

    fi

else
    # MODE was neither "chr" nor "gene" - stop instead of guessing what you meant
    echo "ERROR: second argument must be 'chr' or 'gene' (got '$2')" >&2
    exit 1
fi

# Decide whether extraction was empty here, while creating this gene's VCF.
# pipefail keeps a failed decode from masquerading as a successful zero count.
# Drop sample genotypes during counting to keep this validation inexpensive.
n_variants=$(bcftools view -H -G "$OUTPUT" | awk 'END { print NR + 0 }')
if [[ "$n_variants" -eq 0 ]]; then
    echo "WARNING: no VCF records matched '$3'; this empty gene will be skipped by the catalog workflow." >&2
fi

# Each parallel extraction writes only its own integer sidecar. The workflow
# merges them into one skip report after all extraction jobs have finished.
# Publish only after successful filtering AND successful output decoding.
if [[ -n "$VCF_RECORD_COUNT_FILE" ]]; then
    mkdir -p "$(dirname "$VCF_RECORD_COUNT_FILE")"
    record_count_tmp=$(mktemp "${VCF_RECORD_COUNT_FILE}.tmp.XXXXXX")
    trap 'rm -f -- "$record_count_tmp"' EXIT
    printf '%s\n' "$n_variants" > "$record_count_tmp"
    mv -f -- "$record_count_tmp" "$VCF_RECORD_COUNT_FILE"
    trap - EXIT
fi

echo "Done -> $OUTPUT" >&2

# Ending the timer 
script_end=$(date +%s)
# Calculating the script_runtime of the script
script_runtime=$((script_end-script_start))

echo "Creating filtered VCF-file in BCFtools completed"
echo "The outputfile is stored as: $OUTPUT"
echo "Retained VCF records: $n_variants"
echo "Runtime: $script_runtime seconds"
