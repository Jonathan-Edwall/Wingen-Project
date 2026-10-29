#!/usr/bin/env python3
"""
compute_environmental_correlation_matrix.py - how correlated are the
environmental predictor variables with EACH OTHER, across the empirical
sample locations?

This is unrelated to the Gene Linkage Screening Test - no gene, no genetic
diversity value, no gene-pair enters into it anywhere. It is a property of
the environmental inputs alone: is latitude basically redundant with
winter temperature at these sample sites? One matrix per catalog run, not
one per gene.

Source data: every gene's {gene}_per_location.tsv already carries each
location's environmental values - written once per gene by
2_Herring_Wingen.Rmd's WRITE_PER_LOCATION_TSV block. Since these values
don't depend on genetics, every gene's file has (subject to that gene's
own missing-data pattern) the same numbers for the same cluster_id. This
script globs ALL of them, groups by cluster_id ("sample locations" here
means the same 10 km complete-linkage location CLUSTERS used everywhere
else in this pipeline, not raw individual coordinates), and takes one
representative value per column - null-tolerant, so one gene's missing
value at a location doesn't block using another gene's value at the same
location. Deliberately reads only the per-location values already
extracted at the empirical sample sites, never the full downloaded
CopernicusMarine grid - correlating against the whole raster would answer
a different (and far more expensive) question than "are my predictors
redundant at the places I actually sampled".

Which columns get compared is discovered, not hardcoded: EXCLUDE_COLUMNS
below is deliberately an EXCLUDE list (gene/location bookkeeping fields),
not an include list of named marine variables - every other numeric
column in the per_location.tsv files is treated as an environmental
predictor automatically. Add a new variable (e.g. the water-route distance
to a reference point), and it is picked up the next time this runs with NO
change here - same for a future domain_type=land run with a completely
different, non-overlapping variable set. KNOWN_LABELS only improves the
DISPLAY name for variables this script happens to recognize; an unknown
column still gets compared and shown, just under a plainer auto-generated
label (see prettify_label()) until a nicer one is worth adding.

Similarity metric and standardization: the same reasoning as
compute_gene_correlation_similarity.py applies identically here - plain
Pearson r, no normalization needed first (r is invariant to each
variable's own scale/units; applies the same way whether the two things
being correlated are two genes' profiles or two raw environmental
variables).

Usage:
    python3 compute_environmental_correlation_matrix.py \
        --results-dir results/gene_catalog_herring_secd/results \
        --duckdb results/gene_catalog_herring_secd/results/gene_catalog.duckdb \
        --output-tsv env_correlation_matrix.tsv
"""
import argparse
import glob
import os
import re
import sys
import time

import duckdb

# Bookkeeping/identifier fields written into every {gene}_per_location.tsv
# that are never themselves environmental predictors. Everything else
# numeric is treated as one. lon_mean is excluded deliberately, matching
# how the rest of this pipeline already treats longitude (a coordinate for
# mapping, never used as a gd-correlation predictor the way latitude is).
EXCLUDE_COLUMNS = {
    "gene_id", "cluster_id", "n_individuals", "n_individuals_for_gd", "lon_mean", "gd_stat", "gd",
    "water_dist_range_km", "water_dist_span_flag",  # audit columns of water_dist_km: spread within a cluster, straddle flag
}
# Per-variable data-quality counts, not predictor values.
EXCLUDE_PATTERNS = [re.compile(r"^n_individuals_filled_"), re.compile(r"^n_individuals_no_")]

# Nicer display names for variables this script happens to recognize.
# Purely cosmetic - an unlisted column still participates fully, see
# prettify_label(). lat_mean is the one renamed key (-> "latitude"), to
# match the name the rest of the pipeline already uses for it elsewhere.
KNOWN_LABELS = {
    "lat_mean": "Latitude", "salinity": "Salinity", "temp_summer": "Water temp (summer)",
    "temp_winter": "Water temp (winter)", "chl_log": "log(chlorophyll-a)",
    "oxygen": "Oxygen", "current_speed": "Current speed", "water_dist_km": "Water-route distance (km)",
}


def log(message):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] [env_correlation_matrix] {message}", file=sys.stderr)


def prettify_label(column_name):
    """Fallback label for a column this script doesn't specifically know
    about: 'oresund_distance_km' -> 'Oresund distance km'. Not fancy -
    good enough to be readable until it earns an entry in KNOWN_LABELS."""
    words = column_name.replace("_", " ").split()
    return " ".join([words[0].capitalize()] + words[1:]) if words else column_name


def normalized_key(column_name):
    return "latitude" if column_name == "lat_mean" else column_name


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True, help="{results_dir}/results - globs */*_per_location.tsv beneath it")
    ap.add_argument("--duckdb", required=True, help="gene_catalog.duckdb - written into, same file build_gene_catalog_html.py reads")
    ap.add_argument("--output-tsv", required=True, help="long-format pairs: variable_a, variable_b, r, n_locations")
    args = ap.parse_args()
    script_start = time.perf_counter()

    glob_pattern = os.path.join(args.results_dir, "*", "*_per_location.tsv")
    files = sorted(glob.glob(glob_pattern))
    if not files:
        sys.exit(f"ERROR: no files matched {glob_pattern!r} - check --results-dir, and that "
                 f"WRITE_PER_LOCATION_TSV: TRUE and USE_COPERNICUS_MARINE: TRUE were set for this run.")
    log(f"found {len(files)} per_location.tsv files under {args.results_dir!r}")

    con = duckdb.connect(args.duckdb, read_only=False)
    # R writes a missing value as NA. Without nullstr a column holding one NA in any file is read as text and
    # would be left out of the matrix.
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE all_locations AS
        SELECT * FROM read_csv_auto({files!r}, union_by_name=true, nullstr=['', 'NA'])
    """)
    schema = con.execute("DESCRIBE all_locations").fetchall()
    if "cluster_id" not in {r[0] for r in schema}:
        sys.exit("ERROR: none of the per_location.tsv files have a cluster_id column - "
                 "unexpected schema, check WRITE_PER_LOCATION_TSV output by hand.")

    numeric_types = {"DOUBLE", "FLOAT", "REAL", "BIGINT", "INTEGER", "HUGEINT", "DECIMAL"}
    variables = []  # (source_col, key, label)
    for col_name, col_type, *_ in schema:
        if col_name in EXCLUDE_COLUMNS or any(p.match(col_name) for p in EXCLUDE_PATTERNS):
            continue
        if not any(col_type.upper().startswith(t) for t in numeric_types):
            continue
        key = normalized_key(col_name)
        label = KNOWN_LABELS.get(col_name, prettify_label(key))
        variables.append((col_name, key, label))
    if len(variables) < 2:
        sys.exit(f"ERROR: only {len(variables)} environmental column(s) found (auto-discovered, "
                 f"excluding {sorted(EXCLUDE_COLUMNS)}) across {len(files)} per_location.tsv files - "
                 f"need at least 2 to correlate anything.")
    log(f"auto-discovered {len(variables)} environmental variables: {', '.join(key for _, key, _ in variables)}")

    select_cols = ", ".join(f'max("{src}") AS "{key}"' for src, key, _ in variables)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE locations AS
        SELECT cluster_id, {select_cols}
        FROM all_locations
        GROUP BY cluster_id
    """)
    n_locations = con.execute("SELECT count(*) FROM locations").fetchone()[0]
    log(f"{n_locations} distinct sample-location clusters after deduplicating across genes")

    union_sql = " UNION ALL ".join(
        f'SELECT cluster_id, \'{key}\' AS variable, "{key}" AS val FROM locations' for _, key, _ in variables
    )
    con.execute(f"CREATE OR REPLACE TEMP VIEW env_long AS {union_sql}")

    con.execute("""
        CREATE OR REPLACE TABLE env_correlation_matrix AS
        SELECT a.variable AS variable_a, b.variable AS variable_b,
               corr(a.val, b.val) AS r, count(*) AS n_locations
        FROM env_long a JOIN env_long b ON a.cluster_id = b.cluster_id
        WHERE a.variable <= b.variable
        GROUP BY a.variable, b.variable
        ORDER BY a.variable, b.variable
    """)
    con.execute(
        "COMMENT ON TABLE env_correlation_matrix IS "
        "'Pairwise Pearson r among the environmental predictor variables themselves, across sample "
        "locations - NOT related to gene diversity or the Gene Linkage Screening Test. One row per "
        "unordered variable pair (variable_a <= variable_b, so the diagonal appears once). Variable "
        "set is auto-discovered per run, not hardcoded. '"
        "'See compute_environmental_correlation_matrix.py.'"
    )
    # Labels travel with the data (not just this run's log) so
    # build_gene_catalog_html.py can display a nice name without needing
    # its own copy of KNOWN_LABELS/prettify_label - one source of truth.
    con.execute("CREATE OR REPLACE TABLE env_correlation_matrix_labels AS SELECT * FROM (VALUES "
                + ", ".join(f"('{key}', '{label}')" for _, key, label in variables)
                + ") AS t(variable, label)")
    con.execute(f"COPY env_correlation_matrix TO '{args.output_tsv}' (HEADER, DELIMITER '\t')")

    rows = con.execute("SELECT variable_a, variable_b, r FROM env_correlation_matrix "
                        "WHERE variable_a != variable_b ORDER BY abs(r) DESC").fetchall()
    con.close()

    log(f"wrote {len(variables) * (len(variables) + 1) // 2} matrix entries -> {args.output_tsv}")
    if rows:
        top = rows[0]
        log(f"most correlated pair: {top[0]} vs {top[1]}, r={top[2]:.3f} "
            f"({'candidate for redundancy' if abs(top[2]) >= 0.8 else 'not obviously redundant'})")
    log(f"runtime: {time.perf_counter() - script_start:.2f} seconds")


if __name__ == "__main__":
    main()
