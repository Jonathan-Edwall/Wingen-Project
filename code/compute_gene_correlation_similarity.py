#!/usr/bin/env python3
"""
compute_gene_correlation_similarity.py - "Gene Linkage Screening Test":
which genes have a SIMILAR pattern of genetic-diversity-vs-environment
correlation, as candidates to check for real linkage disequilibrium later.

The N^2 pairwise comparison itself is unavoidable for an EXACT top-K
per gene (there is no way to know which partners are most similar without
comparing every candidate), but nothing below holds, writes or sorts all
N*(N-1)/2 pairs at once. At the project's ceiling of 27,000 genes that is
~364 million pairs - too many to keep in memory, write to a TSV, or fully
SORT for an exact rank. So the genes are handled in blocks
(--block-genes): each block of genes is compared with every gene, its pairs
are reduced to what is kept, and the block is discarded. The memory needed
depends on the block size and not on the number of genes.
  - Every gene's top-K most similar partners (K from --top-k-per-gene,
    default 20) is extracted with a window function over that gene's block
    of pairs - only the N*K-sized result is kept.
  - "Where does this rank overall" uses a histogram of |similarity| over all
    pairs, gathered block by block (HISTOGRAM_BINS bins, so the cutoffs are
    exact to 1e-5 in |similarity|; the smallest and largest values are
    exact), to get a handful of reference cutoffs (50th/90th/95th/99th/
    99.9th percentile of |similarity|). Every KEPT pair's percentile is then
    interpolated against those cutoffs in Python. This is an approximation,
    not an exact percent_rank(); for deciding "is this pair in the
    interesting tail or not", that's the right trade - exactness is not
    reliable to two decimal places anyway, with gene catalogs ranging from a
    dozen genes to 27,000. The result does not change from run to run.
  - The full all-pairs table is only written when the catalog is small
    enough that eyeballing it is actually useful for calibrating what
    similarity level is "interesting" - see --max-full-export-genes. Above
    that it is skipped with a log line, not silently dropped.

Every gene pair is still scored, not just ones clearing a pre-set bar.
The question is "which genes have similar-ish correlations"; what degree
of similarity is interesting enough is deliberately left open, so the
percentile is there to let that threshold be read off the data rather
than guessed in advance.

Similarity metric: Pearson correlation between the two genes' vectors of
LOCATION-level Pearson r (one per environmental variable + latitude).
This is deliberately not a raw distance (Euclidean/Minkowski): those are
sensitive to one gene's correlations being uniformly stronger or weaker
than another's, which is not what "similar pattern" should penalize.
Pearson correlation is exactly cosine similarity computed on each
profile after centering it on its own mean - which is the right notion
of "do these two genes move the same way across the environmental
axes", independent of each gene's own average correlation strength.
Jaccard/simple-matching don't apply here at all - those compare
sets/binary vectors, and binarizing a continuous correlation profile
would throw away the magnitude information the whole test depends on.

No normalization/standardization step is applied before computing
similarity, deliberately: Pearson correlation is invariant to any
affine rescaling of either input vector (corr(aX+b, Y) = corr(X,Y) for
a>0), so z-scoring the profile vectors first would not change a single
similarity value - it would be dead code. (This would NOT be true for a
Euclidean/Minkowski-based similarity, or for feeding these profiles into
PCA/clustering.)

Gene length / variant count is NOT a similarity factor, by design. The
only per-pair "sample size" is n_shared_variables - how many of the ~7-8
environmental/latitude axes BOTH genes have a finite correlation for -
gated by --min-shared-variables, used only as a reliability floor, never
as a weight.

Usage:
    python3 compute_gene_correlation_similarity.py \
        --duckdb gene_catalog.duckdb \
        --min-shared-variables 4 --top-k-per-gene 20 \
        --max-full-export-genes 2000 \
        --output-tsv gene_correlation_similarity_top_matches.tsv
"""
import argparse
import math
import re
import sys
import time

import duckdb

# Quantile levels of the reference curve that kept pairs are placed on, and the resolution of the |similarity| histogram that
# the curve is read from (bin width 1/HISTOGRAM_BINS).
CUTOFF_POINTS = [0.0, 0.5, 0.9, 0.95, 0.99, 0.999, 1.0]
HISTOGRAM_BINS = 100000
# Genes per block. A block holds block_genes x (number of genes) pair groups at a time, roughly 100 bytes each.
DEFAULT_BLOCK_GENES = 256

CANDIDATE_VARIABLES = [
    ("latitude", "latitude", "Latitude"),
    ("salinity", "salinity", "Salinity"),
    ("water_temp_summer", "temp_summer", "Water temp (summer)"),
    ("water_temp_winter", "temp_winter", "Water temp (winter)"),
    ("chlorophyll_a_log", "chl_log", "log(chlorophyll-a)"),
    ("oxygen", "oxygen", "Oxygen"),
    ("current_speed", "current_speed", "Current speed"),
]


def log(message):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] [gene_correlation_similarity] {message}", file=sys.stderr)


def struct_fields(column_type):
    m = re.match(r"STRUCT\((.*)\)$", column_type)
    if not m:
        return set()
    return {f.strip().split(" ")[0].strip('"') for f in m.group(1).split(",")}


def available_variables(con, table, suffix_field):
    """suffix_field: 'r' for Pearson (location_<suffix>_r) or 'rho' for Spearman (location_<suffix>_rho)."""
    types = {name: col_type for name, col_type, *_ in con.execute(f"DESCRIBE {table}").fetchall()}
    available = []
    for group, suffix, label in CANDIDATE_VARIABLES:
        field = f"location_{suffix}_{suffix_field}"
        if group in types and field in struct_fields(types[group]):
            available.append((group, suffix, label))
    return available


def profile_long_sql(table, variables, suffix_field):
    # A gene with no finite value for a variable gets no row for it, so the join of two genes (and
    # n_shared_variables, --min-shared-variables) counts only variables where BOTH have a finite value.
    return " UNION ALL ".join(
        f'SELECT gene_id, \'{suffix}\' AS variable, "{group}".location_{suffix}_{suffix_field} AS val FROM {table} '
        f'WHERE "{group}".location_{suffix}_{suffix_field} IS NOT NULL'
        for group, suffix, _ in variables
    )


def screened_genes_sql(con):
    """The rows of `genes` that take part in the screen: all of them except the pooled neutral reference
    (target_type 'neutral_pool'), which is a reference sample and not a gene."""
    columns = {name for name, *_ in con.execute("DESCRIBE genes").fetchall()}
    if "target_type" not in columns:
        return "SELECT * FROM genes"
    return "SELECT * FROM genes WHERE coalesce(target_type, 'gene') <> 'neutral_pool'"


def histogram_cutoffs(counts, smallest, largest):
    """|similarity| at each of CUTOFF_POINTS, read from a histogram of |similarity| over all pairs (counts[i] pairs with
    |similarity| in [i/HISTOGRAM_BINS, (i+1)/HISTOGRAM_BINS)). The levels 0 and 1 are the exact smallest and largest
    value; in between the value is placed linearly inside the bin that holds the level, and kept within the observed
    range."""
    total = sum(counts)
    cutoffs, passed, i = [], 0, 0
    for level in CUTOFF_POINTS:
        if level <= 0:
            cutoffs.append(smallest)
        elif level >= 1:
            cutoffs.append(largest)
        else:
            target = level * total
            while passed + counts[i] < target:
                passed += counts[i]
                i += 1
            value = (i + (target - passed) / counts[i]) / HISTOGRAM_BINS
            cutoffs.append(min(max(value, smallest), largest))
    return cutoffs


def screen_in_blocks(con, n_genes, min_shared_variables, top_k_per_gene, block_genes):
    """Compare every screened gene with every other one, block_genes genes at a time, from the TEMP view
    pearson_long (gene_id, variable, val). Fills the table top_matches with each gene's top_k_per_gene partners
    (largest |similarity| first, ties by partner id) and returns the |similarity| cutoffs at CUTOFF_POINTS and the
    number of (gene, partner) comparisons, or (None, 0) when no pair qualifies. A pair is compared once from each
    side; both sides give the same value, and the histogram counts both, which leaves every quantile where it is."""
    con.execute("CREATE OR REPLACE TEMP TABLE screen_genes AS "
                "SELECT gene_id, (row_number() OVER (ORDER BY gene_id) - 1)::INTEGER AS pos FROM (SELECT DISTINCT gene_id FROM pearson_long)")
    con.execute("CREATE OR REPLACE TEMP TABLE screen_profile AS "
                "SELECT g.pos, p.variable, p.val FROM pearson_long p JOIN screen_genes g USING (gene_id) ORDER BY g.pos")
    con.execute("CREATE OR REPLACE TABLE top_matches "
                "(gene_id VARCHAR, partner VARCHAR, similarity_r DOUBLE, n_shared_variables BIGINT)")
    counts = [0] * HISTOGRAM_BINS
    smallest, largest, n_pairs = math.inf, -math.inf, 0
    n_blocks = -(-n_genes // block_genes)
    report_every = max(1, n_blocks // 10)
    for block, first in enumerate(range(0, n_genes, block_genes), start=1):
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE block_pairs AS
            SELECT a.pos AS pos, b.pos AS partner_pos, corr(a.val, b.val) AS similarity_r, count(*) AS n_shared_variables
            FROM screen_profile a JOIN screen_profile b ON a.variable = b.variable AND a.pos <> b.pos
            WHERE a.pos >= {first} AND a.pos < {first + block_genes}
            GROUP BY a.pos, b.pos
            HAVING count(*) >= {int(min_shared_variables)} AND isfinite(corr(a.val, b.val))
        """)
        con.execute(f"""
            INSERT INTO top_matches
            SELECT g.gene_id, p.gene_id, k.similarity_r, k.n_shared_variables
            FROM (SELECT pos, partner_pos, similarity_r, n_shared_variables FROM (
                      SELECT *, row_number() OVER (PARTITION BY pos ORDER BY abs(similarity_r) DESC, partner_pos) AS rnk
                      FROM block_pairs) WHERE rnk <= {int(top_k_per_gene)}) k
            JOIN screen_genes g ON g.pos = k.pos JOIN screen_genes p ON p.pos = k.partner_pos
        """)
        for bin_index, n_in_bin, low, high in con.execute(f"""
                SELECT least(floor(abs(similarity_r) * {HISTOGRAM_BINS}), {HISTOGRAM_BINS - 1})::INTEGER AS bin_index,
                       count(*), min(abs(similarity_r)), max(abs(similarity_r))
                FROM block_pairs GROUP BY bin_index""").fetchall():
            counts[bin_index] += n_in_bin
            n_pairs += n_in_bin
            smallest, largest = min(smallest, low), max(largest, high)
        if block % report_every == 0 or block == n_blocks:
            log(f"block {block}/{n_blocks} compared")
    con.execute("DROP TABLE IF EXISTS block_pairs")
    if n_pairs == 0:
        return None, 0
    return histogram_cutoffs(counts, smallest, largest), n_pairs


def interpolate_percentile(value, cutoff_points, cutoff_values):
    """Linear interpolation of `value`'s percentile against the reference
    curve (cutoff_points are the quantile levels, e.g.
    [0, .5, .9, .95, .99, .999, 1], cutoff_values their |similarity| cutoffs).
    An approximation by construction - see the module docstring."""
    if value <= cutoff_values[0]:
        return cutoff_points[0]
    for i in range(1, len(cutoff_values)):
        if value <= cutoff_values[i]:
            lo_v, hi_v = cutoff_values[i - 1], cutoff_values[i]
            lo_p, hi_p = cutoff_points[i - 1], cutoff_points[i]
            if hi_v == lo_v:
                return hi_p
            frac = (value - lo_v) / (hi_v - lo_v)
            return lo_p + frac * (hi_p - lo_p)
    return cutoff_points[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--duckdb", required=True, help="gene_catalog.duckdb from aggregate_gene_summaries.py")
    ap.add_argument("--min-shared-variables", type=int, default=4,
                    help="minimum environmental/latitude variables both genes must have a finite "
                         "Pearson correlation for before a pair is scored at all")
    ap.add_argument("--top-k-per-gene", type=int, default=20,
                    help="a SAFETY CAP on how many of each gene's most similar partners are kept, "
                         "not the primary selection - most genes will show far fewer than this once "
                         "build_gene_catalog_html.py filters down to the ones actually worth "
                         "surfacing (see its DISPLAY_PERCENTILE_FLOOR). Exists only so a gene that "
                         "happens to sit in a genuinely large correlated block doesn't produce an "
                         "unbounded row count; raise it if you ever see it being the limiting factor")
    ap.add_argument("--max-full-export-genes", type=int, default=2000,
                    help="also write the complete all-pairs table, for manual calibration, but only "
                         "when the catalog has at most this many genes (default 2000, ~2M pairs) - "
                         "above that only the top-K-per-gene output and the percentile cutoffs are "
                         "written, so output size stays independent of catalog size")
    ap.add_argument("--output-tsv", required=True, help="per-gene top-K matches")
    ap.add_argument("--output-summary-tsv", help="per-gene single best match (default: alongside --output-tsv)")
    ap.add_argument("--output-full-pairs-tsv", help="complete all-pairs table, when --max-full-export-genes allows it")
    ap.add_argument("--threads", type=int,
                    help="DuckDB thread count. Pass Snakemake's {threads} here explicitly - without "
                         "this, DuckDB auto-detects and may use MORE cores than the rule's own "
                         "threads: declaration promised Snakemake's scheduler, which can oversubscribe "
                         "if other jobs are running concurrently. Omit to let DuckDB auto-detect "
                         "(fine for an interactive/standalone run with nothing else competing).")
    ap.add_argument("--block-genes", type=int, default=DEFAULT_BLOCK_GENES,
                    help="genes compared with all others at a time (default %(default)s). The memory needed is about "
                         "this number x the number of genes x 100 bytes; the result does not depend on it.")
    args = ap.parse_args()
    if args.block_genes < 1:
        sys.exit("ERROR: --block-genes must be at least 1.")
    script_start = time.perf_counter()

    con = duckdb.connect(args.duckdb, read_only=False)
    if args.threads:
        con.execute(f"SET threads={int(args.threads)}")
        log(f"DuckDB threads explicitly set to {args.threads}")
    pearson_vars = available_variables(con, "genes", "r")
    if len(pearson_vars) < 3:
        sys.exit(f"ERROR: only {len(pearson_vars)} environmental/latitude correlation field(s) are "
                 f"present in {args.duckdb} - need at least 3. Check WRITE_LOCATION_STATS: TRUE.")
    spearman_vars = available_variables(con, "genes", "rho")
    log(f"Pearson axes: {', '.join(s for _, s, _ in pearson_vars)}"
        + (f" | Spearman axes: {', '.join(s for _, s, _ in spearman_vars)}" if spearman_vars else " | no Spearman (rho) fields found"))

    con.execute(f"CREATE OR REPLACE TEMP VIEW screened_genes AS {screened_genes_sql(con)}")
    n_genes = con.execute("SELECT count(*) FROM screened_genes").fetchone()[0]
    n_reference = con.execute("SELECT count(*) FROM genes").fetchone()[0] - n_genes
    if n_reference:
        log(f"{n_reference} pooled neutral reference row(s) (target_type 'neutral_pool') left out of the screen")
    if n_genes < 3:
        sys.exit(f"ERROR: only {n_genes} gene(s) in {args.duckdb} - need at least 3 to screen any pairs.")
    n_possible_pairs = n_genes * (n_genes - 1) // 2
    log(f"{n_genes} genes -> {n_possible_pairs:,} possible pairs")

    con.execute(f"CREATE OR REPLACE TEMP VIEW pearson_long AS {profile_long_sql('screened_genes', pearson_vars, 'r')}")

    # Only top_matches (N*K rows) and cutoffs (one row) are stored; the pairs themselves are compared block by block and dropped.
    t_compute = time.perf_counter()
    cutoff_points = CUTOFF_POINTS
    cutoff_values, _n_compared = screen_in_blocks(con, n_genes, args.min_shared_variables, args.top_k_per_gene, args.block_genes)
    if cutoff_values is None:
        sys.exit(f"ERROR: no gene pair had >= {args.min_shared_variables} shared, finite correlation "
                 f"values out of {n_genes} genes - nothing to screen.")
    con.execute("CREATE OR REPLACE TABLE cutoffs AS SELECT ?::DOUBLE[] AS q", [cutoff_values])
    n_matches = con.execute("SELECT count(*) FROM top_matches").fetchone()[0]
    log(f"pairwise computation: {time.perf_counter() - t_compute:.1f}s "
        f"({n_matches:,} kept matches across {n_genes} genes)")
    log(f"reference |similarity| cutoffs (histogram, bin width {1 / HISTOGRAM_BINS:g}): " +
        ", ".join(f"p{p*100:g}={v:.3f}" for p, v in zip(cutoff_points, cutoff_values)))

    # Spearman companion value, computed only for the N*K pairs already kept
    # above (cheap) as a rank-based consistency check alongside the
    # Pearson-based ranking.
    spearman_lookup = {}
    if spearman_vars:
        con.execute(f"CREATE OR REPLACE TEMP VIEW spearman_long AS {profile_long_sql('screened_genes', spearman_vars, 'rho')}")
        rows = con.execute("SELECT gene_id, variable, val FROM spearman_long WHERE val IS NOT NULL").fetchall()
        by_gene = {}
        for gene_id, variable, val in rows:
            by_gene.setdefault(gene_id, {})[variable] = val
        kept_pairs = con.execute("SELECT DISTINCT gene_id, partner FROM top_matches").fetchall()
        for gene_id, partner in kept_pairs:
            va, vb = by_gene.get(gene_id, {}), by_gene.get(partner, {})
            shared = [k for k in va if k in vb]
            if len(shared) < 3:
                continue
            xs = [va[k] for k in shared]
            ys = [vb[k] for k in shared]
            n = len(xs)
            mx, my = sum(xs) / n, sum(ys) / n
            cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            vx = sum((x - mx) ** 2 for x in xs)
            vy = sum((y - my) ** 2 for y in ys)
            if vx > 0 and vy > 0:
                spearman_lookup[(gene_id, partner)] = cov / (vx * vy) ** 0.5

    rows = con.execute("SELECT gene_id, partner, similarity_r, n_shared_variables FROM top_matches "
                        "ORDER BY gene_id, abs(similarity_r) DESC").fetchall()
    import csv
    with open(args.output_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "partner", "similarity_r", "spearman_similarity_r", "percentile_vs_all_pairs", "n_shared_variables"])
        for gene_id, partner, similarity_r, n_shared in rows:
            pct = interpolate_percentile(abs(similarity_r), cutoff_points, cutoff_values)
            spear = spearman_lookup.get((gene_id, partner), spearman_lookup.get((partner, gene_id)))
            w.writerow([gene_id, partner, f"{similarity_r:.6g}",
                       f"{spear:.6g}" if spear is not None else "", f"{pct:.4f}", n_shared])

    summary_path = args.output_summary_tsv or args.output_tsv.replace(".tsv", "_best_per_gene.tsv")
    seen = set()
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "best_match_gene", "best_match_r", "best_match_percentile"])
        for gene_id, partner, similarity_r, n_shared in rows:
            if gene_id in seen:
                continue
            seen.add(gene_id)
            pct = interpolate_percentile(abs(similarity_r), cutoff_points, cutoff_values)
            w.writerow([gene_id, partner, f"{similarity_r:.6g}", f"{pct:.4f}"])

    # Gene IDs stay text even when every ID looks like a number (type inference would make them integers,
    # drop leading zeros and no longer match genes.gene_id).
    con.execute("CREATE OR REPLACE TABLE gene_correlation_similarity_pairs AS "
                f"SELECT * FROM read_csv_auto('{args.output_tsv}', delim='\t', "
                "types={'gene_id': 'VARCHAR', 'partner': 'VARCHAR'})")
    con.execute("CREATE OR REPLACE TABLE gene_correlation_similarity_summary AS "
                f"SELECT * FROM read_csv_auto('{summary_path}', delim='\t', "
                "types={'gene_id': 'VARCHAR', 'best_match_gene': 'VARCHAR'})")

    if n_genes <= args.max_full_export_genes:
        full_path = args.output_full_pairs_tsv or args.output_tsv.replace(".tsv", "_full_all_pairs.tsv")
        con.execute(f"""
            COPY (SELECT a.gene_id AS gene_a, b.gene_id AS gene_b, corr(a.val, b.val) AS similarity_r,
                         count(*) AS n_shared_variables
                  FROM pearson_long a JOIN pearson_long b ON a.variable = b.variable AND a.gene_id < b.gene_id
                  GROUP BY a.gene_id, b.gene_id HAVING count(*) >= {int(args.min_shared_variables)}
                  ORDER BY abs(similarity_r) DESC)
            TO '{full_path}' (HEADER, DELIMITER '\t')
        """)
        log(f"catalog has {n_genes} genes (<= --max-full-export-genes {args.max_full_export_genes}) "
            f"-> also wrote the complete all-pairs table to {full_path}")
    else:
        log(f"catalog has {n_genes} genes (> --max-full-export-genes {args.max_full_export_genes}) "
            f"-> skipping the full all-pairs export; use the percentile cutoffs above and the "
            f"top-{args.top_k_per_gene}-per-gene output instead")

    con.close()
    log(f"wrote {n_matches:,} top-match rows -> {args.output_tsv}")
    log(f"wrote {len(seen)} per-gene summary rows -> {summary_path}")
    log(f"total runtime: {time.perf_counter() - script_start:.1f}s")


if __name__ == "__main__":
    main()
