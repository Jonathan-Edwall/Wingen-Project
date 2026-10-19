#!/bin/bash -l
# Create neutral marker coordinates using the bcftools query / inline AWK
# approach of snpeff_vcf_detect_unique_genes.sh. This script is independent
# of the gene-list script and has no external parser dependency.
#
# Usage: bash snpeff_vcf_create_neutral_marker_file.sh [INPUT_VCF] [OUTPUT]
# INPUT_VCF: SnpEff-annotated plain/gzip/BGZF VCF (no index required).
# OUTPUT: file or directory; default neutral_markers.tsv beside INPUT_VCF.
# Output is a headered CHROM/POS TSV, with original 1-based SNP positions.
# No intervals are added; extract_neutral_pool reads this as normalized.
#
# Retain biallelic SNPs only when EVERY comma-separated ANN entry and EVERY
# '&'-joined effect is intergenic_region. Flanking gene names/IDs are allowed.
# Upstream/downstream, gene-body, regulatory and conserved-intergenic effects,
# missing/malformed ANN and annotation errors/warnings are excluded.
set -eo pipefail
script_start=$(date +%s)

####################################
# Environment setup
####################################
# As in the gene-list script, containers use their own dependencies. For a
# local run, use bcftools already on PATH or fall back to the same wingen_env.
# Activate before set -u because conda activation may read unset variables.
if [[ -n "${APPTAINER_CONTAINER:-}" || -n "${SINGULARITY_CONTAINER:-}" ]]; then
    echo "Running inside the Apptainer container - skipping host conda activation." >&2
elif ! command -v bcftools >/dev/null; then
    conda_setup_script_path="/proj/uppmax2025-2-125/nobackup/jonathan/anaconda3/etc/profile.d/conda.sh"
    if [[ ! -f "$conda_setup_script_path" ]]; then
        echo "ERROR: bcftools is required on PATH; the local wingen_env setup was not found." >&2
        exit 1
    fi
    source "$conda_setup_script_path"
    conda activate wingen_env
fi
set -u
command -v bcftools >/dev/null || { echo "ERROR: bcftools is required on PATH." >&2; exit 1; }

####################################
# Input and output
####################################
# Keep the gene script's cohort fallback, but resolve explicit paths from the
# caller's working directory so relative Snakemake paths remain valid.
project_root="${ROOT:-/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project}"
raw_data_dir="${project_root}/data/raw"
input_vcf_file="${1:-${raw_data_dir}/genetic_variation_herring/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz}"
if [[ ! -f "$input_vcf_file" ]]; then
    echo "ERROR: input VCF '$input_vcf_file' does not exist or is not a regular file." >&2
    exit 1
fi
OUTPUT="${2:-$(dirname -- "$input_vcf_file")/neutral_markers.tsv}"
if [[ -d "$OUTPUT" || "$OUTPUT" == */ ]]; then
    OUTPUT="${OUTPUT%/}/neutral_markers.tsv"
fi
if [[ "$input_vcf_file" -ef "$OUTPUT" ]]; then
    echo "ERROR: output must not overwrite the input VCF." >&2
    exit 1
fi
mkdir -p -- "$(dirname -- "$OUTPUT")"

# Stage output beside its destination so errors cannot publish partial marker
# lists or overwrite a previous successful file. Cleanup touches only these
# known files in the freshly allocated working directory.
work_dir="$(mktemp -d "$(dirname -- "$OUTPUT")/.neutral-snpeff.XXXXXX")"
cleanup() {
    rm -f -- "$work_dir/header.vcf" "$work_dir/markers.tsv" "$work_dir/summary.txt"
    rmdir -- "$work_dir"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

####################################
# Confirm SnpEff ANN and extract neutral positions
####################################
bcftools view -h "$input_vcf_file" > "$work_dir/header.vcf"
if ! grep '^##INFO=<ID=ANN,' "$work_dir/header.vcf" >/dev/null; then
    echo "ERROR: '$input_vcf_file' has no INFO/ANN declaration; a SnpEff-annotated VCF is required." >&2
    exit 1
fi

# The gene-list script splits ANN into comma-separated blocks, then splits
# each block's effect field on '&'. Here we apply that same parsing pattern,
# but reject the entire SNP if any block/effect is not intergenic_region.
# REF/ALT are queried too so indels and multiallelic sites are excluded.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "$input_vcf_file" |
    awk -v summary_file="$work_dir/summary.txt" '
        BEGIN { FS = OFS = "\t"; print "CHROM", "POS" }
        {
            chrom = $1; pos = $2; ref = $3; alt = $4; ann = $5
            if (NF != 5 || chrom == "" || chrom == ".") {
                print "ERROR: malformed VCF query record at row " NR > "/dev/stderr"
                fatal = 1; exit 1
            }
            if (pos !~ /^[0-9]+$/ || pos + 0 < 1) {
                print "ERROR: POS must be a positive integer at row " NR > "/dev/stderr"
                fatal = 1; exit 1
            }
            records++
            if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/ || ref == alt) { non_snp++; next }
            if (ann == "." || ann == "") { missing_ann++; next }

            invalid = non_intergenic_site = 0
            m = split(ann, blocks, ",")
            for (j = 1; j <= m; j++) {
                nf = split(blocks[j], f, /[|]/)
                # Full ANN has 16 fields. Gene fields 4/5 are not filters:
                # SnpEff may fill them with IDs of the flanking genes.
                if (nf != 16 || f[1] != alt || f[2] == "" ||
                    f[3] !~ /^(HIGH|MODERATE|LOW|MODIFIER)$/) invalid = 1
                nmessages = split(f[16], messages, "&")
                for (i = 1; i <= nmessages; i++)
                    if (messages[i] ~ /^(ERROR|WARNING)/) invalid = 1

                ne = split(f[2], eparts, "&")
                for (ei = 1; ei <= ne; ei++)
                    if (eparts[ei] != "intergenic_region") non_intergenic_site = 1
            }
            if (invalid) invalid_ann++
            else if (non_intergenic_site) non_intergenic++
            else { selected++; print chrom, pos }
        }
        END {
            if (fatal) exit 1
            summary = sprintf("Scanned %d VCF records; selected %d intergenic biallelic SNP records. " \
                "Excluded: non_snp=%d, non_intergenic=%d, missing_ann=%d, invalid_ann=%d",
                records, selected, non_snp, non_intergenic, missing_ann, invalid_ann)
            if (!selected) {
                print summary > "/dev/stderr"
                print "ERROR: No eligible neutral SNPs; supply markers_file or neutral_vcf for a curated reference." > "/dev/stderr"
                exit 1
            }
            print summary > summary_file
        }
    ' > "$work_dir/markers.tsv"

mv -- "$work_dir/markers.tsv" "$OUTPUT"
cat "$work_dir/summary.txt"
echo "The neutral marker file is stored as: $OUTPUT"
echo "Runtime: $(($(date +%s) - script_start)) seconds"
