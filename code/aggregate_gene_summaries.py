"""
Snakemake `script:` for the `aggregate_gene_summaries` rule.

Runs exactly once, after every gene's .done marker exists - that is what
makes this safe: Snakemake's own DAG already guarantees this rule doesn't
start until every gene has finished, so there's no concurrent-reader-vs-
writer problem to solve here either.

Each gene writes its OWN file, "{gene}_gene_summary.json", into a folder
SHARED by every gene on that chromosome ({results_dir}/results/{chrom}/,
where results_dir is config.yaml's results_dir - see the Snakefile) -
the gene-ID prefix is what keeps genes on the same chromosome from
overwriting each other, not a per-gene folder. The glob below matches
that: one level of {chrom} wildcard, then any *_gene_summary.json file
directly inside it (not gene_summary.json inside a further {gene}/
subfolder). It's built from --duckdb's own directory (see below) rather
than hardcoding "results/gene_catalog", so this stays correct for
whichever results_dir a given catalog's config.yaml actually names.

DuckDB's read_json_auto() takes the glob directly and infers a schema from
the union of fields it finds, so genes that errored out before computing
every field (e.g. a gene with 0 variants that never reached the
gene-summary chunk) just get NULLs there rather than breaking the load -
worth querying for after a big run (see the printed null-count summary
below) rather than assuming every row is complete. A nested field (e.g.
gene_summary.json's "windows": {"variance": ..., "mean": ..., "median":
...}) is read as a genuine STRUCT column - the .duckdb file and the html
catalog (build_gene_catalog_html.py) both use that structure directly, but
struct cells don't make a great flat-file TSV, so the TSV export below
flattens any top-level struct into "structname_field" columns instead.

Which genetic diversity statistic a catalog holds (config.yaml's gd_stat: pi,
Ho, biallelic_richness, tajimas_d or watterson_theta) travels with the data:
the Rmd writes gd_stat, gd_label and gd_label_short into every gene's JSON, so they show up here as
plain columns - first in the TSV, right next to gene_id - and
build_gene_catalog_html.py takes its column labels from them rather than from
config.yaml. Columns such as windows_mean or location_mean are just numbers
and do not say which statistic they are a mean of, so a catalog must never mix
statistics: that is checked right after the load below.

Per-gene statistics are passed through without multiple-testing adjustment.

Gene names: with --gene-list (the gene-list BED), the optional columns 5 (short gene name) and
6 (long description) are added to the catalog as `gene_name` and `gene_description`, right after
`gene_id`; a missing value (".", NA or empty) becomes NULL. They are joined in HERE, at the very
end, rather than passed into each gene's analysis, so a name can change without any gene being
re-run. A BED with only four columns, or no --gene-list at all, leaves the catalog's columns
exactly as they were.
"""
import argparse
import json
import os

import duckdb


def quote_identifier(name):
    """Quote one actual field name, including any embedded double quotes."""
    return '"' + name.replace('"', '""') + '"'


MISSING_NAME = {"", ".", "NA"}


def read_gene_names(path):
    """gene_id -> [short name or None, description or None] from a gene-list BED.

    Columns 1-4 are chrom, start, end, gene_id; the optional columns 5 and 6 are the short name
    and the long description. Lines are split the same way the Snakefile's read_gene_bed does
    (on tabs when there is one). A value on any line of a gene fills a gap left by its other
    lines; a gene with neither value is left out.
    """
    names = {}
    with open(path, encoding="utf-8-sig") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 4:
                raise RuntimeError(f"{path}, line {number}: expected at least 4 columns "
                                   f"(chrom, start, end, gene_id), got {len(fields)}")
            gene = fields[3].strip()
            values = [fields[i].strip() if len(fields) > i else "" for i in (4, 5)]
            values = [None if v in MISSING_NAME else v for v in values]
            if gene and any(values):
                known = names.setdefault(gene, [None, None])
                for i, value in enumerate(values):
                    if known[i] is None:
                        known[i] = value
    return names


parser = argparse.ArgumentParser()
parser.add_argument("--duckdb", required=True)
parser.add_argument("--tsv", required=True)
parser.add_argument("--availability-manifest", help="Read only retained summaries from this run.")
parser.add_argument("--gene-list", help="Gene-list BED; its optional columns 5 and 6 (gene name, description) "
                                        "become the catalog's gene_name and gene_description columns.")
args = parser.parse_args()

duckdb_path = args.duckdb
tsv_path = args.tsv

# The gene-catalog directory (parent of every {chrom}/ folder) is wherever
# --duckdb itself is being written - by the Snakefile's own convention
# that's always "{results_dir}/results", so deriving it from
# duckdb_path here means this script never needs a second --results-dir/
# --gene-catalog-dir flag that could drift out of sync with whatever path
# the Snakefile actually passed for --duckdb. `or "."` only matters if
# duckdb_path is a bare filename with no directory part at all (e.g. run
# by hand as `--duckdb gene_catalog.duckdb` from inside that folder).
gene_catalog_dir = os.path.dirname(duckdb_path) or "."
SUMMARY_GLOB = os.path.join(gene_catalog_dir, "*", "*_gene_summary.json")

# Restrict a workflow run to its retained genes and explicitly enabled pooled
# reference. Stale skipped/removed/disabled summaries are not loaded. Validate
# before opening DuckDB so an all-skipped run cannot create catalog outputs.
summary_source = SUMMARY_GLOB
if args.availability_manifest:
    with open(args.availability_manifest, encoding="utf-8") as handle:
        manifest = json.load(handle)
    retained = manifest.get("retained")
    if manifest.get("schema_version") != 1 or not isinstance(retained, list):
        raise RuntimeError("Invalid gene availability manifest: expected schema_version 1 and retained list.")
    reference = manifest.get("neutral_reference")
    if reference is not None:
        if (not isinstance(reference, dict) or reference.get("target_type") != "neutral_pool"
                or not isinstance(reference.get("gene_id"), str) or not reference["gene_id"].strip()
                or not isinstance(reference.get("summary_path"), str) or not reference["summary_path"].strip()):
            raise RuntimeError("Invalid neutral reference record in gene availability manifest.")
        if any(record.get("gene_id") == reference["gene_id"] for record in retained):
            raise RuntimeError("Neutral reference gene_id collides with a retained gene.")
    # A reference-only run is explicit (gene_list_limit=0 / requested_count=0).
    # Enabling the reference must not hide an all-skipped real-gene run.
    if not retained and not (reference is not None and manifest.get("requested_count") == 0):
        raise RuntimeError("No genes have retained VCF records; see skipped_genes.tsv. No catalog was created.")
    summary_source = []
    for record in retained + ([reference] if reference is not None else []):
        path = record["summary_path"]
        if not os.path.isfile(path):
            raise RuntimeError(f"Missing summary for retained gene {record['gene_id']!r}: {path}")
        with open(path, encoding="utf-8") as handle:
            summary = json.load(handle)
        if not isinstance(summary, dict) or summary.get("gene_id") != record["gene_id"]:
            raise RuntimeError(f"Summary gene_id does not match retained gene {record['gene_id']!r}: {path}")
        summary_source.append(path)

con = duckdb.connect(duckdb_path)
# import os
# import re

# import duckdb

# SUMMARY_GLOB = "results/gene_catalog/*/*_gene_summary.json"

# con = duckdb.connect(snakemake.output.duckdb)

con.execute(
    """
    CREATE OR REPLACE TABLE genes AS
    SELECT *
    FROM read_json_auto(?, filename = true, union_by_name = true)
    """,
    [summary_source],
)
# con.execute(
#     f"""
#     CREATE OR REPLACE TABLE genes AS
#     SELECT *
#     FROM read_json_auto('{SUMMARY_GLOB}', filename = true, union_by_name = true)
#     """
# )

# filename -> summary_path (clearer name) + gene_dir (the folder
# build_gene_catalog_html.py will look in for this gene's plots - note
# that's the SHARED per-chromosome folder now, not a per-gene one, so
# that script has to filter what it finds there by gene_id, not just
# list the whole folder) - both computed once, here, so neither the Rmd
# nor the HTML builder need to know the glob pattern above.
con.execute("ALTER TABLE genes RENAME filename TO summary_path")
con.execute("ALTER TABLE genes ADD COLUMN gene_dir VARCHAR")
con.execute(
    r"UPDATE genes SET gene_dir = regexp_replace(summary_path, '[/\\][^/\\]*_gene_summary\.json$', '')"
)

# chr: gene_dir's own basename (e.g. {results_dir}/results/chr1 -> chr1). No
# per-gene JSON carries a "chr" field on its own - only the shared
# per-chromosome folder does - so without this, build_gene_catalog_html.py's
# COLUMN_LAYOUT entry for "chr" has nothing to bind to and available_columns()
# silently drops the whole column from the rendered table.
con.execute("ALTER TABLE genes ADD COLUMN chr VARCHAR")
con.execute(r"UPDATE genes SET chr = regexp_replace(gene_dir, '^.*[/\\]', '')")

# Gene names, from the gene-list BED's optional columns. A LEFT JOIN rebuild (not UPDATE) keeps every
# STRUCT column intact and puts the two new columns right after gene_id in the TSV.
# ORDER BY g.rowid: the join may be planned with `genes` on the build side, which emits the unmatched
# (unnamed) genes LAST, so on a large catalog the TSV / DuckDB row order would silently depend on whether
# names are attached. rowid is the insertion order (= the manifest order), i.e. the order this script
# always produced, so attaching names changes only the two new columns.
gene_names = read_gene_names(args.gene_list) if args.gene_list else {}
if args.gene_list and not gene_names:
    print(f"gene names: {args.gene_list} has no name or description columns - the catalog gets no name columns")
if gene_names:
    con.execute("CREATE TEMP TABLE gene_names (gene_id VARCHAR, gene_name VARCHAR, gene_description VARCHAR)")
    con.executemany("INSERT INTO gene_names VALUES (?, ?, ?)",
                    [(gene, values[0], values[1]) for gene, values in gene_names.items()])
    present = {r[0] for r in con.execute("DESCRIBE genes").fetchall()}
    pick = lambda col: f"COALESCE(n.{col}, g.{col})" if col in present else f"n.{col}"
    exclude = ", ".join(quote_identifier(c) for c in ["gene_id"] + [c for c in ("gene_name", "gene_description") if c in present])
    con.execute(f"""
        CREATE TABLE genes_named AS
        SELECT g.gene_id, {pick('gene_name')} AS gene_name, {pick('gene_description')} AS gene_description,
               g.* EXCLUDE ({exclude})
        FROM genes g LEFT JOIN gene_names n ON g.gene_id = n.gene_id
        ORDER BY g.rowid
    """)
    con.execute("DROP TABLE genes")
    con.execute("ALTER TABLE genes_named RENAME TO genes")
    n_short, n_desc = con.execute("SELECT count(gene_name), count(gene_description) FROM genes").fetchone()
    print(f"gene names: {n_short} gene(s) have a short name and {n_desc} a description (from {args.gene_list})")

n_genes = con.execute("SELECT count(*) FROM genes").fetchone()[0]
if n_genes == 0:
    raise RuntimeError(
        f"No files matched {SUMMARY_GLOB!r} - is calc_gen_diversity_stats "
        f"TRUE in config.yaml? (gene_summary.json is only written when it is)"
    )

# One catalog = one genetic diversity statistic. A JSON from before gd_stat existed has no
# such field (it was always pi), so a catalog where NONE of the genes carry it is an old but
# consistent one - while a mix of "has gd_stat" and "hasn't" (or of two different values) is
# a half-rebuilt results folder whose numbers can't be compared column by column.
genes_columns = [r[0] for r in con.execute("DESCRIBE genes").fetchall()]
if "gd_stat" in genes_columns:
    by_stat = con.execute(
        "SELECT gd_stat, count(*) FROM genes GROUP BY gd_stat ORDER BY count(*) DESC, gd_stat"
    ).fetchall()
    if len(by_stat) > 1:
        detail = ", ".join(
            f"{stat if stat is not None else 'no gd_stat (older run, = pi)'}: {n} gene(s)"
            for stat, n in by_stat
        )
        raise RuntimeError(
            f"The gene summaries under {gene_catalog_dir} were written with different genetic "
            f"diversity statistics ({detail}). One catalog can only hold one - use a separate "
            f"results_dir per gd_stat, or delete the stale *_gene_summary.json files and re-run "
            f"those genes."
        )
    print(f"genetic diversity statistic: {by_stat[0][0]}")
else:
    print("genetic diversity statistic: pi (these summaries predate the gd_stat field)")

# One-line-per-column null-count so a bad batch is visible immediately
# rather than discovered later inside the html catalog. For a struct
# column (e.g. "windows") this counts genes missing the struct entirely
# (old-schema genes without that field at all), not genes with just one
# null sub-field - a coarser check, but still the first thing worth
# knowing about a bad batch.
cols = [r[0] for r in con.execute("DESCRIBE genes").fetchall() if r[0] not in ("gene_dir", "gene_name", "gene_description")]
null_counts = con.execute(
    "SELECT " + ", ".join(
        f"count(*) - count({quote_identifier(c)}) AS {quote_identifier(c)}"
        for c in cols
    ) + " FROM genes"
).fetchone()
print(f"aggregated {n_genes} gene summaries into {duckdb_path}")
# print(f"aggregated {n_genes} gene summaries into {snakemake.output.duckdb}")
for col, n_null in zip(cols, null_counts):
    if n_null:
        print(f"  NOTE: {col} is NULL for {n_null}/{n_genes} genes")

# TSV export: expand any top-level STRUCT column (e.g. windows -> struct
# with variance/mean/median fields) into separate windows_variance /
# windows_mean / windows_median columns, rather than one cell per gene
# holding a stringified {'variance': ..., 'mean': ..., 'median': ...} -
# a TSV's one job is to open cleanly in a spreadsheet. Only handles a
# struct of plain scalar fields (which is what this pipeline produces) -
# not structs nested inside structs or inside lists.
select_parts = []
for name, col_type, *_ in con.execute("DESCRIBE genes").fetchall():
    quoted_name = quote_identifier(name)
    if col_type.startswith("STRUCT("):
        # Bind the struct to get actual field names. Parsing its type string
        # double-quotes reserved names such as "scope" and splits comma types.
        fields = con.execute(f"SELECT {quoted_name}.* FROM genes LIMIT 0").description
        select_parts += [
            f"{quoted_name}.{quote_identifier(field[0])} "
            f"AS {quote_identifier(name + '_' + field[0])}"
            for field in fields
        ]
    else:
        select_parts.append(quoted_name)

con.execute(
    f"COPY (SELECT {', '.join(select_parts)} FROM genes) "
    "TO ? (HEADER, DELIMITER '\\t')",
    [tsv_path],
)
# con.execute(
#     f"COPY (SELECT {', '.join(select_parts)} FROM genes) "
#     f"TO '{snakemake.output.tsv}' (HEADER, DELIMITER '\t')"
# )

con.close()
