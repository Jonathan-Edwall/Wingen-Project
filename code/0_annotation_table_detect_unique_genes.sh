#!/bin/bash

# Start the timer, as in 0_snpeff_vcf_detect_unique_genes.sh.
script_start=$(date +%s)
set -euo pipefail

####################################
# Options
####################################
# Usage: 0_annotation_table_detect_unique_genes.sh bed|txt INPUT_TABLE [OUTPUT] [options]
#
# Same format/input/output ordering as the SnpEff script after its scope arg.
# There is no scope arg: supplied gene_start/gene_end define the whole interval.
# Required: a headered TSV with chromosome, gene start, gene end and gene ID.
# Defaults match ALL_SNP_Gene_Annotations.txt. Repeated SNP-gene rows become one
# gene; a gene's chromosome/bounds/biotype must agree on every genic row.
# BED uses real gene bounds (start-1, end), not the min/max observed SNP POS.
# txt keeps IDs in first-seen order; BED is naturally chromosome/start sorted.
# Output can be a file or directory; an existing file is never overwritten.
# Only Bash, awk and ordinary GNU coreutils are needed, plus gzip for .gz input.
# No hardcoded project root or conda activation is needed for these tools.

usage() {
    cat <<'HELP'
Usage: 0_annotation_table_detect_unique_genes.sh bed|txt INPUT_TABLE [OUTPUT] [options]

  --chrom-column NAME          default: CHROM
  --start-column NAME          default: gene_start
  --end-column NAME            default: gene_end
  --gene-id-column NAME        default: gene_id
  --biotype-column NAME        default: auto (use gene_biotype if present)
  --annotation-column NAME     default: auto (use annotation if present)
  --pos-column NAME            default: auto (use SNP POS if present)
  --genic-label VALUE          default: genic
  --include-biotypes A,B       keep these exact biotypes; default: all
  --exclude-biotypes A,B       discard these exact biotypes; exclusion wins
  --input-coordinates VALUE    one-based-inclusive (default) or bed
  --header                     BED only: add #chromosome, pos1, pos2, gene_id

Optional columns can be disabled with NAME=none. Explicit names must exist.
With a status column, only --genic-label rows qualify. Without one, valid
gene IDs and bounds qualify. Missing IDs are skipped. Optional SNP positions
are always 1-based and are checked against the supplied gene interval.
Gene scope is implicit in the input bounds; exon/CDS scope cannot be inferred.
OUTPUT defaults to unique_genes.bed/txt beside INPUT_TABLE. Output is BED4
(CHROM, 0-based start, exclusive end, gene ID), or one unique ID per text line.
BED is headerless by default; --header adds a pipeline-compatible comment.
HELP
}

fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then usage; exit 0; fi
[[ $# -ge 2 ]] || { usage >&2; exit 1; }
FORMAT="$1"
INPUT_TABLE="$2"
shift 2
[[ "$FORMAT" == bed || "$FORMAT" == txt ]] || fail "format must be 'bed' or 'txt' (got '$FORMAT')"
OUTPUT=""
if [[ $# -gt 0 && "$1" != --* ]]; then OUTPUT="$1"; shift; fi

CHROM_COLUMN=CHROM
START_COLUMN=gene_start
END_COLUMN=gene_end
GENE_ID_COLUMN=gene_id
BIOTYPE_COLUMN=auto
ANNOTATION_COLUMN=auto
POS_COLUMN=auto
GENIC_LABEL=genic
INCLUDE_BIOTYPES=""
EXCLUDE_BIOTYPES=""
INPUT_COORDINATES=one-based-inclusive
WRITE_HEADER=false
while [[ $# -gt 0 ]]; do
    if [[ "$1" == --header ]]; then WRITE_HEADER=true; shift; continue; fi
    [[ $# -ge 2 ]] || fail "option '$1' requires a value"
    [[ -n "$2" && "$2" != *$'\n'* && "$2" != *$'\t'* && "$2" != *$'\r'* ]] || fail "invalid empty/multiline value for '$1'"
    case "$1" in
        --chrom-column) CHROM_COLUMN="$2" ;;
        --start-column) START_COLUMN="$2" ;;
        --end-column) END_COLUMN="$2" ;;
        --gene-id-column) GENE_ID_COLUMN="$2" ;;
        --biotype-column) BIOTYPE_COLUMN="$2" ;;
        --annotation-column) ANNOTATION_COLUMN="$2" ;;
        --pos-column) POS_COLUMN="$2" ;;
        --genic-label) GENIC_LABEL="$2" ;;
        --include-biotypes) INCLUDE_BIOTYPES="$2" ;;
        --exclude-biotypes) EXCLUDE_BIOTYPES="$2" ;;
        --input-coordinates) INPUT_COORDINATES="$2" ;;
        *) fail "unknown option '$1'; use --help" ;;
    esac
    shift 2
done
[[ "$WRITE_HEADER" == false || "$FORMAT" == bed ]] || fail "--header is supported only with bed output"
[[ "$INPUT_COORDINATES" == one-based-inclusive || "$INPUT_COORDINATES" == bed ]] || fail "input coordinates must be 'one-based-inclusive' or 'bed'"
[[ -f "$INPUT_TABLE" ]] || fail "input table does not exist: $INPUT_TABLE"
OUTPUT="${OUTPUT:-$(dirname -- "$INPUT_TABLE")/unique_genes.${FORMAT}}"
if [[ -d "$OUTPUT" || "$OUTPUT" == */ ]]; then OUTPUT="${OUTPUT%/}/unique_genes.${FORMAT}"; fi
if [[ "$FORMAT" == bed && "${OUTPUT,,}" != *.bed ]]; then fail "BED output must end in .bed so the pipeline recognizes it"; fi
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]] || fail "refusing to overwrite existing output: $OUTPUT"
mkdir -p -- "$(dirname -- "$OUTPUT")"
work_dir=$(mktemp -d -- "${OUTPUT}.tmp.XXXXXX")
raw_output="$work_dir/genes.raw"
sorted_output="$work_dir/genes.sorted"
cleanup() {
    rm -f -- "$raw_output" "$sorted_output"
    rmdir -- "$work_dir"
}
trap cleanup EXIT

####################################
# Read the annotation table
####################################
read_annotation_table() {
    if [[ "${INPUT_TABLE,,}" == *.gz ]]; then
        gzip -cd -- "$INPUT_TABLE"
    else
        cat -- "$INPUT_TABLE"
    fi
}

####################################
# Classify each row, deduplicate genes, write the requested format
####################################
# Like the SnpEff awk loop, classify each association before adding its gene.
# A SNP can overlap several genes; excluding one association must not discard
# the other gene. Unlike ANN, this table already supplies the full gene bounds.
# Stage all output so a bad row or a broken gzip stream never publishes a list.
read_annotation_table | LC_ALL=C awk \
    -v format="$FORMAT" -v chrom_col="$CHROM_COLUMN" -v start_col="$START_COLUMN" \
    -v end_col="$END_COLUMN" -v gene_col="$GENE_ID_COLUMN" \
    -v biotype_col="$BIOTYPE_COLUMN" -v annotation_col="$ANNOTATION_COLUMN" \
    -v pos_col="$POS_COLUMN" -v genic_label="$GENIC_LABEL" \
    -v include_types="$INCLUDE_BIOTYPES" -v exclude_types="$EXCLUDE_BIOTYPES" \
    -v input_coordinates="$INPUT_COORDINATES" '
    function error(message) {
        failed = 1
        printf "ERROR: annotation table line %d: %s\n", NR, message > "/dev/stderr"
        exit 1
    }
    function missing(value) { return value == "" || value == "NA" || value == "." }
    function column(name, fallback, optional) {
        if (optional && name == "none") return 0
        if (optional && name == "auto") return fallback in header ? header[fallback] : 0
        if (!(name in header)) error("missing header column: " name)
        return header[name]
    }
    function parse_types(text, target, parts, n, i, label) {
        if (text == "") return
        n = split(text, parts, ",")
        for (i = 1; i <= n; i++) {
            label = parts[i]
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", label)
            if (label == "") error("empty biotype label in filter")
            target[label] = 1
        }
    }
    function integer(text, minimum, label, value) {
        if (text !~ /^[0-9]+$/) error(label " must be an integer, got: " text)
        value = text + 0
        if (value < minimum || value > 9007199254740991) error("invalid " label ": " text)
        return value
    }
    BEGIN { FS = "\t"; OFS = "\t" }
    {
        sub(/\r$/, "", $0)
        if ($0 == "") next
        if (!header_seen) {
            sub(/^\357\273\277/, "", $1)
            width = NF
            for (i = 1; i <= NF; i++) {
                if ($i == "" || $i in header) error("empty or duplicate header: " $i)
                header[$i] = i
            }
            c = column(chrom_col, "", 0)
            s = column(start_col, "", 0)
            e = column(end_col, "", 0)
            g = column(gene_col, "", 0)
            if (c == s || c == e || c == g || s == e || s == g || e == g)
                error("chromosome, start, end and gene ID require distinct columns")
            b = column(biotype_col, "gene_biotype", 1)
            a = column(annotation_col, "annotation", 1)
            p = column(pos_col, "POS", 1)
            if (pos_col == "auto" && (p == c || p == s || p == e || p == g)) p = 0
            if ((include_types != "" || exclude_types != "") && !b)
                error("biotype filtering requires a biotype column")
            parse_types(include_types, include)
            parse_types(exclude_types, exclude)
            header_seen = 1
            next
        }
        if (NF != width) error("expected " width " tab-separated fields, got " NF)
        total_rows++
        if (a && ($a "") != (genic_label "")) next
        gene = $g ""
        if (missing(gene)) next
        if (gene !~ /^[A-Za-z0-9_][A-Za-z0-9_.:-]*$/)
            error("gene ID contains characters incompatible with pipeline filenames: " gene)
        chrom = $c ""
        if (missing(chrom) || chrom ~ /[[:space:]]/) error("missing or invalid chromosome for " gene)
        start = integer($s, input_coordinates == "bed" ? 0 : 1, "gene start")
        end = integer($e, 1, "gene end")
        start0 = input_coordinates == "bed" ? start : start - 1
        if (end <= start0) error("empty or reversed gene interval for " gene)
        if (p) {
            pos = integer($p, 1, "SNP position")
            if (pos <= start0 || pos > end) error("SNP position outside gene interval for " gene)
        }
        biotype = b ? ($b "") : "NA"
        if (missing(biotype)) biotype = "NA"
        if (gene in first_line) {
            if ((chrom_of[gene] "") != (chrom "") || start_of[gene] != start0 || end_of[gene] != end || (biotype_of[gene] "") != (biotype ""))
                error("conflicting chromosome, boundaries or biotype for " gene " (first seen at line " first_line[gene] ")")
        } else {
            first_line[gene] = NR
            chrom_of[gene] = chrom
            start_of[gene] = start0
            end_of[gene] = end
            biotype_of[gene] = biotype
        }
        if (include_types != "" && !(biotype in include)) next
        if (biotype in exclude) next
        selected_rows++
        if (!(gene in seen)) { seen[gene] = 1; order[++n_genes] = gene }
    }
    END {
        if (failed) exit 1
        if (!header_seen) { print "ERROR: empty input; expected a TSV header" > "/dev/stderr"; exit 1 }
        if (!n_genes) { print "ERROR: no retained genes; check column names, annotation label and biotype filters" > "/dev/stderr"; exit 1 }
        for (i = 1; i <= n_genes; i++) {
            gene = order[i]
            if (format == "txt") print gene
            else printf "%s\t%.0f\t%.0f\t%s\n", chrom_of[gene], start_of[gene], end_of[gene], gene
        }
        printf "Selected %d of %d annotation rows; %d unique genes\n", selected_rows, total_rows, n_genes > "/dev/stderr"
    }
' > "$raw_output"

if [[ "$FORMAT" == bed ]]; then
    # Keep a requested header before the sorted BED4 data, never among the genes.
    if [[ "$WRITE_HEADER" == true ]]; then
        printf '#chromosome\tpos1\tpos2\tgene_id\n' > "$sorted_output"
    else
        : > "$sorted_output"
    fi
    LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -k3,3n -k4,4 -- "$raw_output" >> "$sorted_output"
else
    cp -- "$raw_output" "$sorted_output"
fi
[[ -s "$sorted_output" ]] || fail "no gene list was produced"
# Same filesystem: publishing a completed hard link is atomic and refuses an
# existing target, including one created since the initial existence check.
ln -T -- "$sorted_output" "$OUTPUT"
n_genes=$(wc -l < "$raw_output")
script_end=$(date +%s)
script_runtime=$((script_end-script_start))
echo "Gene Catalog Created"
echo "Interval: supplied gene boundaries | Format: $FORMAT | $n_genes unique genes"
echo "The outputfile is stored as: $OUTPUT"
echo "Runtime: $script_runtime seconds"
