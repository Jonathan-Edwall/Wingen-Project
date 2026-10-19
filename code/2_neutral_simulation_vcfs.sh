#!/usr/bin/env bash
# Neutral simulations, step 2: prepare the neutral pool and selected subset VCFs.
#
# bcftools extracts and indexes variants here. It does not draw random SNPs.
# 3_neutral_simulation_subsets.py supplies seeded subset plans and checks the
# records, sample order and saved provenance. Only Python's standard library,
# Bash and the existing bcftools installation are required.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: bash 2_neutral_simulation_vcfs.sh MODE [options]

extract-pool: extract eligible neutral SNPs directly to a compressed VCF.
  --vcf SOURCE --markers FILE --output POOL.vcf.gz
  [--markers-format normalized|positions|bed]
  [--target-markers-output TARGET_MARKERS.tsv]

  normalized: CHROM/POS table from step 1.
  positions:  headerless CHROM/POS list; positions are 1-based.
  bed:        BED3 containing one SNP per interval; start is 0-based, end exclusive.
  Ready-made positions/BED files bypass the optional step 1 table converter.
  Writes the pool and mapped positions. Snakemake tracks these completed files.

validate-pool: validate an extracted pool and create its CSI index.
  --vcf POOL.vcf.gz --markers TARGET_MARKERS.tsv --index-output POOL.vcf.gz.csi
  Reads the existing compressed pool; it never repeats source extraction.

pool: standalone convenience mode: extract, validate and index in one call.
  Accepts the same options as extract-pool. The workflow uses the two separate
  modes above so Snakemake can reuse completed extraction after a later failure.

existing-pool: validate an already-neutral BGZF VCF and create its cache index.
  --vcf NEUTRAL.vcf.gz --index-output CACHE.csi
  The supplied VCF is read directly; it is not copied or rewritten.

materialize: extract and index the SNP subsets selected by step 3.
  --pool POOL.vcf.gz --gene-vcf GENE.vcf --plan SUBSETS.tsv
  --metadata PLAN.json --output-dir DIRECTORY
  [--pool-index CACHE.csi]
  Writes one VCF per technical replicate.

All modes accept:
  --threads INTEGER      Total CPU budget, including the main thread. Default: 1
                         Uses remaining threads for bcftools compression/indexing.
  --bcftools EXECUTABLE   Default: bcftools
  --python EXECUTABLE     Default: python
  --helper PATH           Default: sibling 3_neutral_simulation_subsets.py

Use --help after the mode, or with no mode, to display this help.
USAGE
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

log() { printf '[%s] [%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "${rule_name:-neutral_preparation}" "$*" >&2; }
run_timed() {
    local phase="$1"; shift
    local started=$SECONDS result
    log "START $phase"
    if "$@"; then
        log "DONE $phase elapsed=$((SECONDS-started))s"
    else
        result=$?
        log "FAILED $phase elapsed=$((SECONDS-started))s exit=$result"
        return "$result"
    fi
}

mode="${1:---help}"
case "$mode" in
    --help|-h) usage; exit 0 ;;
    pool|extract-pool|validate-pool|existing-pool|materialize) shift ;;
    *) die "Unknown mode '$mode'; use --help" ;;
esac
case "$mode" in
    pool|extract-pool) rule_name="extract_neutral_pool" ;;
    validate-pool) rule_name="validate_neutral_pool" ;;
    existing-pool) rule_name="prepare_existing_neutral_pool" ;;
    materialize) rule_name="materialize_neutral_vcfs" ;;
esac

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
helper="$script_dir/3_neutral_simulation_subsets.py"
python="python"
bcftools="bcftools"
threads="1"
vcf=""; markers=""; output=""; metadata=""
target_markers_output=""; index_output=""
markers_format="normalized"
pool=""; gene_vcf=""; plan=""; output_dir=""; pool_index=""

while (($#)); do
    option="$1"
    case "$option" in
        --help|-h) usage; exit 0 ;;
    esac
    # Accept both --option value and --option=value.
    if [[ "$option" == --*=* ]]; then
        value="${option#*=}"
        option="${option%%=*}"
        shift
    else
        (($# >= 2)) || die "Missing value for '$option'"
        value="$2"
        shift 2
    fi
    case "$mode:$option" in
        *:--bcftools) bcftools="$value" ;;
        *:--python) python="$value" ;;
        *:--helper) helper="$value" ;;
        *:--threads) threads="$value" ;;
        materialize:--metadata) metadata="$value" ;;
        pool:--vcf|extract-pool:--vcf|validate-pool:--vcf|existing-pool:--vcf) vcf="$value" ;;
        pool:--markers|extract-pool:--markers|validate-pool:--markers) markers="$value" ;;
        pool:--markers-format|extract-pool:--markers-format) markers_format="$value" ;;
        pool:--output|extract-pool:--output) output="$value" ;;
        pool:--target-markers-output|extract-pool:--target-markers-output) target_markers_output="$value" ;;
        existing-pool:--index-output|validate-pool:--index-output) index_output="$value" ;;
        materialize:--pool) pool="$value" ;;
        materialize:--gene-vcf) gene_vcf="$value" ;;
        materialize:--plan) plan="$value" ;;
        materialize:--output-dir) output_dir="$value" ;;
        materialize:--pool-index) pool_index="$value" ;;
        *) die "Unknown option '$option' for mode '$mode'" ;;
    esac
done

# bcftools requests worker threads in addition to its main processing thread.
# Validate before creating files or calling tools, including integer overflow.
[[ "$threads" =~ ^[1-9][0-9]*$ && ${#threads} -le 10 ]] &&
    ((threads <= 2147483647)) || die "--threads must be a positive integer no larger than 2147483647"
bcftools_threads=$((threads - 1))

case "$mode" in
    pool|extract-pool)
        [[ -n "$vcf" && -n "$markers" && -n "$output" ]] ||
            die "pool requires --vcf, --markers and --output"
        case "$markers_format" in
            normalized|positions|bed) ;;
            *) die "--markers-format must be normalized, positions or bed" ;;
        esac
        target_markers_output="${target_markers_output:-${output%.*}.target_markers.tsv}"
        work_parent="$(dirname -- "$output")"
        ;;
    validate-pool)
        [[ -n "$vcf" && -n "$markers" && -n "$index_output" ]] ||
            die "validate-pool requires --vcf, --markers and --index-output"
        work_parent="$(dirname -- "$index_output")"
        ;;
    existing-pool)
        [[ -n "$vcf" && -n "$index_output" ]] ||
            die "existing-pool requires --vcf and --index-output"
        work_parent="$(dirname -- "$index_output")"
        ;;
    materialize)
        [[ -n "$pool" && -n "$gene_vcf" && -n "$plan" && -n "$metadata" && -n "$output_dir" ]] ||
            die "materialize requires --pool, --gene-vcf, --plan, --metadata and --output-dir"
        work_parent="$output_dir"
        ;;
esac
[[ -f "$helper" ]] || die "Neutral subset helper not found: $helper"
command -v "$python" >/dev/null || die "Python executable not found: $python"
command -v "$bcftools" >/dev/null || die "bcftools executable not found: $bcftools"
bcftools_version="$("$bcftools" --version-only)"

# Every temporary file is inside this newly allocated directory. Cleanup never
# deletes output_dir or a user-supplied path; outputs publish after their phase succeeds.
mkdir -p -- "$work_parent"
work_dir="$(mktemp -d "$work_parent/.neutral-${mode}.XXXXXX")"
cleanup() { rm -rf -- "$work_dir"; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

extract_pool() {
    local destination_pool="$1" destination_markers="$2"
    run_timed "read source VCF header" "$bcftools" view --no-version --header-only -Ov \
        -o "$work_dir/header.vcf" "$vcf"
    local contig_args=()
    if ! grep -q '^##contig=<ID=' "$work_dir/header.vcf"; then
        run_timed "discover source contigs (header has no contig declarations)" \
            "$bcftools" query -f '%CHROM\n' "$vcf" > "$work_dir/contigs.tsv"
        contig_args=(--contigs "$work_dir/contigs.tsv")
    fi
    run_timed "map neutral positions to source contig names" "$python" "$helper" pool-target-markers \
        --vcf-header "$work_dir/header.vcf" "${contig_args[@]}" \
        --markers "$markers" --markers-format "$markers_format" \
        --target-markers-output "$work_dir/mapped_positions.tsv"
    # Filter in bcftools before compression, preserving original genotypes,
    # annotations and sample order. No full plain-text VCF is written.
    run_timed "extract and BGZF-compress neutral SNPs (CPU budget=$threads)" \
        "$bcftools" view --no-version --threads "$bcftools_threads" \
        -T "$work_dir/mapped_positions.tsv" --targets-overlap pos \
        -m2 -M2 -v snps -i 'REF~"^[ACGT]$" && ALT~"^[ACGT]$" && REF!=ALT' \
        -Oz -o "$work_dir/extracted.vcf.gz" "$vcf"
    mkdir -p -- "$(dirname -- "$destination_pool")" "$(dirname -- "$destination_markers")"
    mv -f -- "$work_dir/mapped_positions.tsv" "$destination_markers"
    mv -f -- "$work_dir/extracted.vcf.gz" "$destination_pool"
}

validate_pool() {
    local source_pool="$1" target_positions="$2" destination_index="$3"
    run_timed "validate compressed neutral SNPs (no VCF rewrite)" "$python" "$helper" check-pool-records \
        --candidates "$source_pool" --target-markers "$target_positions" --validate-only
    run_timed "index compressed neutral pool" "$bcftools" index --force --csi \
        --threads "$bcftools_threads" --output "$work_dir/validated_pool.csi" "$source_pool"
    mkdir -p -- "$(dirname -- "$destination_index")"
    mv -f -- "$work_dir/validated_pool.csi" "$destination_index"
}

case "$mode" in
    pool|extract-pool)
        "$python" "$helper" guard-paths \
            --input "$vcf" --input "$markers" \
            --output "$output" --output "$output.csi" --output "$target_markers_output"
        if [[ "$mode" == extract-pool ]]; then
            extract_pool "$output" "$target_markers_output"
        else
            # Standalone convenience mode retains all-or-nothing publication.
            extract_pool "$work_dir/staged_pool.vcf.gz" "$work_dir/staged_markers.tsv"
            validate_pool "$work_dir/staged_pool.vcf.gz" "$work_dir/staged_markers.tsv" "$work_dir/staged_pool.csi"
            mkdir -p -- "$(dirname -- "$output")" "$(dirname -- "$target_markers_output")"
            mv -f -- "$work_dir/staged_markers.tsv" "$target_markers_output"
            mv -f -- "$work_dir/staged_pool.vcf.gz" "$output"
            mv -f -- "$work_dir/staged_pool.csi" "$output.csi"
        fi
        ;;
    validate-pool)
        "$python" "$helper" guard-paths --input "$vcf" --input "$markers" --output "$index_output"
        validate_pool "$vcf" "$markers" "$index_output"
        ;;
    existing-pool)
        "$python" "$helper" guard-paths --input "$vcf" \
            --output "$index_output"
        "$python" "$helper" check-pool --vcf "$vcf"
        # The index lives in the writable cache; the source may be read-only.
        "$bcftools" index --force --csi --threads "$bcftools_threads" --output "$index_output" "$vcf"
        ;;
    materialize)
        index_args=()
        indexed_pool="$pool"
        if [[ -n "$pool_index" ]]; then
            index_args=(--pool-index "$pool_index")
            indexed_pool="${pool}##idx##${pool_index}"
        fi
        run_timed "prepare technical-replicate position lists" "$python" "$helper" prepare-subsets --pool "$pool" --gene-vcf "$gene_vcf" \
            --plan "$plan" --metadata "$metadata" --output-dir "$output_dir" \
            --work-dir "$work_dir" "${index_args[@]}"

        prepared_count=0
        replicate_count="$(wc -l < "$work_dir/groups.tsv")"
        log "START create $replicate_count technical-replicate VCFs"
        while IFS= read -r filename || [[ -n "$filename" ]]; do
            [[ "$filename" =~ ^replicate[0-9]+\.vcf\.gz$ ]] ||
                die "Unexpected neutral subset filename: $filename"
            "$bcftools" view --no-version -R "$work_dir/$filename.positions.tsv" \
                --regions-overlap pos --threads "$bcftools_threads" -Oz -o "$work_dir/$filename" "$indexed_pool"
            "$python" "$helper" check-subset --vcf "$work_dir/$filename" \
                --manifest "$work_dir/manifest.json" --filename "$filename"
            "$bcftools" index --force --csi --threads "$bcftools_threads" "$work_dir/$filename"
            prepared_count=$((prepared_count + 1))
            if ((prepared_count % 10 == 0 || prepared_count == replicate_count)); then
                log "Created and checked $prepared_count/$replicate_count technical-replicate VCFs"
            fi
        done < "$work_dir/groups.tsv"

        # No output is published until every requested subset has passed checks.
        "$python" "$helper" finish-subsets --manifest "$work_dir/manifest.json" \
            --output "$work_dir/final_manifest.json" --bcftools-version "$bcftools_version"
        while IFS= read -r filename || [[ -n "$filename" ]]; do
            mv -f -- "$work_dir/$filename" "$output_dir/$filename"
            mv -f -- "$work_dir/$filename.csi" "$output_dir/$filename.csi"
        done < "$work_dir/groups.tsv"
        mv -f -- "$work_dir/final_manifest.json" "$output_dir/manifest.json"
        log "DONE published $prepared_count technical-replicate VCFs and manifest"
        ;;
esac
