#!/usr/bin/env bash
# Standalone, read-only audit of the existing neutral-marker ANN filter.
# Usage: bash audit_snpeff_annotations.sh INPUT.vcf.gz NEW_OUTPUT_DIR [REGION]
# REGION is optional, uses indexed access, and must use the VCF's contig names.
set -euo pipefail
export LC_ALL=C
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: bash $0 INPUT.vcf.gz NEW_OUTPUT_DIR [REGION]" >&2
    exit 2
fi
vcf="$1"
output_dir="$2"
command -v bcftools >/dev/null || { echo "bcftools is required on PATH." >&2; exit 2; }
[[ -f "$vcf" ]] || { echo "Input VCF not found: $vcf" >&2; exit 2; }
mkdir -p -- "$(dirname -- "$output_dir")"
mkdir -- "$output_dir"  # Require a new directory; never replace an earlier audit.
query_args=()
if [[ $# -eq 3 ]]; then query_args=(-r "$3"); fi
{
    printf 'input\t%s\nregion\t%s\nstarted\t%s\n' "$vcf" "${3:-ALL}" "$(date -u +%FT%TZ)"
    bcftools --version | head -n 1
} > "$output_dir/run.txt"

# No sample/genotype fields are requested. Only counters and three examples
# per rejection reason are retained; the full annotation stream is not saved.
bcftools query "${query_args[@]}" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "$vcf" |
awk -v examples="$output_dir/examples.tsv" '
    BEGIN {
        FS = OFS = "\t"
        print "reason", "CHROM", "POS", "REF", "ALT", "ANN" > examples
    }
    function report(section, item, count) {
        printf "%s\t%s\t%.0f\t%.6f\n", section, item, count, (total ? 100 * count / total : 0)
    }
    {
        if (NF != 5 || $1 == "" || $1 == "." || $2 !~ /^[0-9]+$/ || $2 + 0 < 1) {
            print "Malformed query record at row " NR > "/dev/stderr"
            fatal = 1; exit 1
        }
        total++
        if (total % 1000000 == 0) printf "Scanned %.0f records\n", total > "/dev/stderr"
        if ($3 !~ /^[ACGT]$/ || $4 !~ /^[ACGT]$/ || $3 == $4) { non_snp++; next }
        if ($5 == "." || $5 == "") { missing_ann++; next }
        delete reasons
        delete effects
        non_intergenic_site = has_intergenic = 0
        m = split($5, blocks, ",")
        for (j = 1; j <= m; j++) {
            nf = split(blocks[j], f, /[|]/)
            if (nf != 16) reasons["field_count_" nf] = 1
            if (f[1] != $4) reasons["allele_mismatch"] = 1
            if (f[2] == "") reasons["empty_effect"] = 1
            if (f[3] !~ /^(HIGH|MODERATE|LOW|MODIFIER)$/) reasons["invalid_impact"] = 1
            nmessages = split(f[16], messages, "&")
            for (i = 1; i <= nmessages; i++)
                if (messages[i] ~ /^(ERROR|WARNING)/) reasons["message:" messages[i]] = 1
            ne = split(f[2], eparts, "&")
            for (ei = 1; ei <= ne; ei++) {
                effects[eparts[ei]] = 1
                if (eparts[ei] == "intergenic_region") has_intergenic = 1
                else non_intergenic_site = 1
            }
        }
        invalid = 0
        for (reason in reasons) invalid = 1
        if (invalid) {
            invalid_ann++
            for (reason in reasons) {
                reason_counts[reason]++
                if (reason_counts[reason] <= 3) print reason, $1, $2, $3, $4, $5 > examples
            }
            for (effect in effects) effect_counts[effect]++
            group = has_intergenic ? (non_intergenic_site ? "intergenic_plus_other_effects" : "only_intergenic_effects") : "no_intergenic_effect_seen"
            group_counts[group]++
        } else if (non_intergenic_site) non_intergenic++
        else selected++
    }
    END {
        if (fatal) exit 1
        print "section", "item", "sites", "percent_of_all_records"
        report("filter_outcome", "total", total)
        report("filter_outcome", "selected", selected)
        report("filter_outcome", "non_intergenic", non_intergenic)
        report("filter_outcome", "invalid_ann", invalid_ann)
        report("filter_outcome", "missing_ann", missing_ann)
        report("filter_outcome", "non_snp", non_snp)
        for (reason in reason_counts) report("invalid_reason", reason, reason_counts[reason])
        for (effect in effect_counts) report("invalid_effect", effect, effect_counts[effect])
        for (group in group_counts) report("invalid_effect_group", group, group_counts[group])
    }
' > "$output_dir/summary.tsv.partial"

mv -- "$output_dir/summary.tsv.partial" "$output_dir/summary.tsv"
printf 'completed\t%s\n' "$(date -u +%FT%TZ)" >> "$output_dir/run.txt"
cat "$output_dir/summary.tsv"
echo "Audit complete: $output_dir" >&2
