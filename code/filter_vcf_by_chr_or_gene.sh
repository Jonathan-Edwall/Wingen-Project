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

####################################
# Running the filter
####################################

if [[ "$MODE" == "chr" ]]; then

    # Check chromosome naming convention in the input VCF.
    # Get the CHROM value from the first actual variant row.
    VCF_FIRST_CHROM=$(bcftools query -f '%CHROM\n' "$INPUT" | sed -n '1p')

    if [[ -z "$VCF_FIRST_CHROM" ]]; then
        echo "ERROR: Could not determine chromosome naming from '$INPUT'." >&2
        exit 1
    fi

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

        # Match VALUE chromosome naming to the VCF chromosome naming.
        if [[ "$VCF_FIRST_CHROM" != chr* && "$VALUE_CHROM" == chr* ]]; then
            VALUE_CHROM="${VALUE_CHROM#chr}"
            echo "VCF uses chromosome names without 'chr'."
            echo "Removed 'chr' prefix from VALUE."

        elif [[ "$VCF_FIRST_CHROM" == chr* && "$VALUE_CHROM" != chr* ]]; then
            VALUE_CHROM="chr${VALUE_CHROM}"
            echo "VCF uses 'chr' chromosome names."
            echo "Added 'chr' prefix to VALUE."
        fi

        # Reconstruct the region.
        VALUE="${VALUE_CHROM}:${VALUE_COORDS}"

    else

        # VALUE is just a chromosome, e.g. chr19 or 19.

        if [[ "$VCF_FIRST_CHROM" != chr* && "$VALUE" == chr* ]]; then
            VALUE="${VALUE#chr}"
            echo "VCF uses chromosome names without 'chr'."
            echo "Removed 'chr' prefix from VALUE."

        elif [[ "$VCF_FIRST_CHROM" == chr* && "$VALUE" != chr* ]]; then
            VALUE="chr${VALUE}"
            echo "VCF uses 'chr' chromosome names."
            echo "Added 'chr' prefix to VALUE."
        fi

    fi

    echo "VCF chromosome naming detected from first variant: $VCF_FIRST_CHROM"
    echo "Filtering in $MODE with -r (VALUE): $VALUE"

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
    bcftools view -r "$VALUE" -O "$OUTPUT_TYPE" -o "$OUTPUT" "$INPUT"

elif [[ "$MODE" == "gene" ]]; then

    # gene mode always takes a single gene ID - filtering several genes, or filtering
    # from a BED file of coordinates, are both handled a level up (see the header
    # comment above and 2_wingen_create_gene_catalog.sh).

    VCF_HEADER=$(bcftools view -h "$INPUT")
    CHROM_ARGS=()
    if [[ -n "$GENE_CHROM" ]]; then
        # Keep old chr-prefixed mapping files compatible with VCFs using bare names.
        mapfile -t VCF_CONTIGS < <(printf '%s\n' "$VCF_HEADER" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p')
        if [[ ${#VCF_CONTIGS[@]} -eq 0 ]]; then
            mapfile -t VCF_CONTIGS < <(bcftools query -f '%CHROM\n' "$INPUT" | sort -u)
        fi
        REQUESTED_CANONICAL="$GENE_CHROM"
        [[ "${REQUESTED_CANONICAL,,}" != chr* ]] || REQUESTED_CANONICAL="${REQUESTED_CANONICAL:3}"
        CHROM_MATCH=""
        CHROM_ALIASES=()
        for CONTIG in "${VCF_CONTIGS[@]}"; do
            if [[ "$CONTIG" == "$GENE_CHROM" ]]; then
                CHROM_MATCH="$CONTIG"
                break
            fi
            CONTIG_CANONICAL="$CONTIG"
            [[ "${CONTIG_CANONICAL,,}" != chr* ]] || CONTIG_CANONICAL="${CONTIG_CANONICAL:3}"
            [[ "$CONTIG_CANONICAL" != "$REQUESTED_CANONICAL" ]] || CHROM_ALIASES+=("$CONTIG")
        done
        if [[ -z "$CHROM_MATCH" && ${#CHROM_ALIASES[@]} -eq 1 ]]; then
            CHROM_MATCH="${CHROM_ALIASES[0]}"
        fi
        if [[ -z "$CHROM_MATCH" ]]; then
            echo "ERROR: chromosome '$GENE_CHROM' has no unambiguous match in '$INPUT'." >&2
            exit 1
        fi
        CHROM_ARGS=(-t "$CHROM_MATCH")
    fi

    # Confirmed against the VCF's own header, whether plain or bgzip-compressed.
    if grep -q '^##INFO=<ID=ANN,' <<< "$VCF_HEADER"; then

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

        n_matches=$(bcftools view -H "$OUTPUT" 2>/dev/null | wc -l)
        if [ "$n_matches" -eq 0 ]; then
            echo "WARNING: 0 variants matched '$VALUE'. With no known field structure to" >&2
            echo "         confirm against, this could mean the gene has none, OR that" >&2
            echo "         this VCF has no gene-tagging info in INFO at all - check a" >&2
            echo "         known annotated line's INFO field by eye before trusting this." >&2
        fi
    fi

else
    # MODE was neither "chr" nor "gene" - stop instead of guessing what you meant
    echo "ERROR: second argument must be 'chr' or 'gene' (got '$2')" >&2
    exit 1
fi

echo "Done -> $OUTPUT" >&2

# Ending the timer 
script_end=$(date +%s)
# Calculating the script_runtime of the script
script_runtime=$((script_end-script_start))

echo "Creating filtered VCF-file in BCFtools completed"
echo "The outputfile is stored as: $OUTPUT"
echo "Runtime: $script_runtime seconds"
