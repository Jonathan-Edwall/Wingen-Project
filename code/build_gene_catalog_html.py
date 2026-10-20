"""
Snakemake `script:` for the `build_gene_catalog_html` rule.

Reads the catalog DuckDB built by aggregate_gene_summaries, finds each
gene's own plots within its (shared-per-chromosome) folder by filename,
and renders templates/gene_catalog.html through Jinja2.

Column layout supports nested grouping, to any depth (COLUMN_LAYOUT
below) so gene_summary.json can nest related stats under one key - e.g. a
gene's "windows" field as {"variance": ..., "mean": ..., "median": ...}
instead of three separate top-level pi_variance_raw_windows /
pi_mean_windows / pi_median_windows keys, or a covariate's Window/Location
halves nested a level deeper still (see _covariate_group below).
DuckDB's read_json_auto() turns that nested object
into a genuine STRUCT column on load (confirmed: DESCRIBE genes shows
`windows STRUCT(variance DOUBLE, mean DOUBLE, median DOUBLE)`), which
comes back from a Python fetch as an ordinary nested dict - so this script
doesn't need to do anything special to produce nested JSON, `dict(zip(...))`
already gives one. The grouped header itself (the multi-row <thead>,
colspan/rowspan etc.) is built client-side in the template from the same
COLUMN_LAYOUT list, sent through unchanged as JSON - one source of truth
for which columns exist and how they're grouped, rather than the header
markup and the DataTables column config having to be kept in sync by hand.

Which genetic diversity statistic the catalog is about (gd_stat in config.yaml) is read back
from the data itself - see read_gd_metric() - and only changes the column labels and two
field-name prefixes; the layout is otherwise the same for every statistic.

Run once, at the end of the DAG, on purpose (see the rule's own comment in
the Snakefile) - re-run any time via `snakemake --forcerun build_gene_catalog_html`.
"""
import argparse
import csv
import io
import json
import math
import os
import re

import duckdb
from jinja2 import Environment, FileSystemLoader
# import json
# import os

# import duckdb
# from jinja2 import Environment, FileSystemLoader

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")

# One entry per top-level field. A plain column: {"key", "label", "numeric"}.
# A grouped column: {"group", "columns": [...]} - "key" inside a group is
# dot-notation into the struct (DataTables resolves "windows.variance"
# against row.windows.variance natively, no custom lookup code needed).
# A group's own "columns" can themselves be groups (see _covariate_group
# below) - nesting is not limited to one level; the client-side header
# builder and available_columns() below both recurse to whatever depth
# the layout actually uses.


# Extra "Location" columns per covariate: signed correlations, conventional Pearson p and
# the neutral reference comparison written by the Rmd. Set SHOW_LOCATION_STATS = False
# to hide them again; a column whose data does not exist in a given run is dropped automatically either way.
SHOW_LOCATION_STATS = True
LEGACY_EXTREME_VALUE_PERCENTILES = (95.0, 99.0)


def _percentile_text(percentile):
    """Match R's JSON field suffixes without scientific notation."""
    return f"{percentile:.6f}".rstrip("0").rstrip(".")


def _percentile_ordinal(percentile):
    label = _percentile_text(percentile)
    if not percentile.is_integer():
        return f"{label}th"
    n = int(percentile)
    ending = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{label}{ending}"


def _alpha_text(percentile):
    """The reference upper-tail fraction, with the percentile's precision."""
    return f"{(100 - percentile) / 100:.8f}".rstrip("0").rstrip(".")


def add_window_correlations(records, gd_key):
    """Fill display r from older simple-intercept fits without changing saved data.

    New summaries provide r directly. Older Rmd fits used the same observation
    pairs and an intercept, so r = sign(slope) * sqrt(R²). A missing or invalid
    fit stays unavailable. Counts let layout pruning see these in-memory fields.
    """
    counts = {}
    for group, suffix in (("latitude", "latitude"), ("salinity", "salinity"),
                          ("water_temp_summer", "temp_summer"),
                          ("water_temp_winter", "temp_winter"),
                          ("chlorophyll_a_log", "chl_log"),
                          ("oxygen", "oxygen"), ("current_speed", "current_speed")):
        prefix = f"window_{gd_key}_{suffix}"
        found = False
        count = 0
        for record in records:
            values = record.get(group)
            if not isinstance(values, dict):
                continue
            key = prefix + "_r"
            if key not in values and not all(prefix + tail in values for tail in ("_r2", "_slope")):
                continue
            found = True
            if values.get(key) is None:
                r2, slope = (values.get(prefix + tail) for tail in ("_r2", "_slope"))
                windows = record.get("windows") or {}
                constant_windows = isinstance(windows, dict) and windows.get("variance") == 0
                valid = all(isinstance(value, (int, float)) and not isinstance(value, bool)
                            and math.isfinite(value) for value in (r2, slope)) and not constant_windows
                values[key] = (
                    math.copysign(math.sqrt(min(r2, 1.0)), slope)
                    if valid and 0 <= r2 <= 1 + 1e-12 and (slope != 0 or r2 == 0)
                    else None
                )
            value = values[key]
            count += isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
        if found:
            counts[f"{group}.{prefix}_r"] = count
    return counts


def read_extreme_value_percentiles(rows, existing_cols):
    """Read the cutoffs actually computed, with 95/99 for older summaries.

    R stores this setting as a scalar CSV value so DuckDB JSON inference
    cannot collapse a one-element vector into a scalar. Mixed settings
    require rebuilding the affected summaries before rendering a catalog.
    """
    settings = set()
    for row in rows:
        metadata = dict(zip(existing_cols, row)).get("neutral_model") or {}
        raw = metadata.get("extreme_value_percentiles")
        if "extreme_value_percentiles" not in metadata:
            values = LEGACY_EXTREME_VALUE_PERCENTILES
        else:
            values = (() if not raw.strip() else raw.split(",")) if isinstance(raw, str) else raw
            try:
                if any(isinstance(value, bool) for value in values):
                    raise ValueError("Boolean percentile")
                values = tuple(round(float(value), 6) for value in values)
            except (TypeError, ValueError) as exc:
                raise RuntimeError("Invalid neutral_model.extreme_value_percentiles metadata.") from exc
            if any(not math.isfinite(p) or not 0 < p < 100 for p in values):
                raise RuntimeError("Neutral reference percentiles must be numbers strictly between 0 and 100.")
            if len(set(values)) != len(values):
                raise RuntimeError("Neutral reference percentiles must be distinct after rounding to 6 decimals.")
        settings.add(tuple(values))
    if len(settings) > 1:
        raise RuntimeError(
            "The catalog mixes neutral reference percentile settings. "
            "Rebuild the gene summaries with one extreme_value_percentiles setting."
        )
    return next(iter(settings), LEGACY_EXTREME_VALUE_PERCENTILES)


def _neutral_cutoff_columns(key_prefix, suffix, extreme_value_percentiles, coefficient):
    symbol = "r" if coefficient == "r" else "\u03c1"
    method = "Pearson" if coefficient == "r" else "Spearman"
    return [
        {
            "key": f"{key_prefix}.location_{suffix}_neutral_abs_{coefficient}_q{_percentile_text(p).replace('.', '_')}",
            "label": f"Neutral |{symbol}| (\u03b1={_alpha_text(p)})", "numeric": True,
            "tooltip": f"c_alpha(|{symbol}|): estimated {_percentile_ordinal(float(p))} percentile of absolute {method} correlations across valid neutral subsets. Alpha={_alpha_text(p)} is the chosen nominal significance level (reference upper-tail fraction); this cell contains a correlation magnitude, not a p-value. Compare the magnitude of the gene's {method} coefficient with its own cutoff. Quantile interpolation can differ from the discrete add-one neutral p-value comparison when the simulation has finitely many subsets.",
        }
        for p in extreme_value_percentiles
    ]


def _location_stat_columns(key_prefix, suffix, extreme_value_percentiles=LEGACY_EXTREME_VALUE_PERCENTILES,
                           gd_short=None):
    if not SHOW_LOCATION_STATS:
        return []
    groups = []
    for method, coefficient, p_suffix, valid_suffix, status_suffix in (
        ("Pearson", "r", "p_neutral", "neutral_n_valid", "neutral_status"),
        ("Spearman", "rho", "p_neutral_rho", "neutral_n_valid_rho", "neutral_rho_status"),
    ):
        symbol = "r" if coefficient == "r" else "\u03c1"
        columns = [
            {"key": f"{key_prefix}.location_{suffix}_{p_suffix}", "label": "p neutral", "numeric": True,
             "tooltip": f"Monte Carlo p-value comparing the magnitude of the observed {method} correlation with the neutral subsets: (1 + number at least as extreme) / (1 + valid subsets).",
             "neutral_test": {
                 "predictor": key_prefix,
                 "method": method.lower(),
                 "p_key": f"{key_prefix}.location_{suffix}_{p_suffix}",
                 "status_key": f"{key_prefix}.location_{suffix}_{status_suffix}",
                 "coefficient_key": f"{key_prefix}.location_{suffix}_{coefficient}",
                 "absolute_key": f"{key_prefix}.location_{suffix}_abs_{coefficient}",
                 "alphas": [round((100 - float(p)) / 100, 8) for p in extreme_value_percentiles],
             }},
            {"key": f"{key_prefix}.location_{suffix}_{coefficient}", "label": symbol, "numeric": True,
             "tooltip": "Pearson correlation of diversity and predictor values." if coefficient == "r" else "Spearman rank correlation of diversity and predictor values."},
            {"key": f"{key_prefix}.location_{suffix}_abs_{coefficient}",
             "source_key": f"{key_prefix}.location_{suffix}_{coefficient}", "transform": "abs",
             "label": f"|{symbol}|", "numeric": True,
             "tooltip": f"Absolute observed {method} correlation for filtering by strength in either direction. Compare with this method's neutral threshold, or filter p neutral to assess the neutral comparison."},
        ] + _neutral_cutoff_columns(key_prefix, suffix, extreme_value_percentiles, coefficient)
        if method == "Pearson":
            columns.append({"key": f"{key_prefix}.location_{suffix}_p", "label": "p Pearson", "numeric": True,
                            "location_extra": True,
                            "tooltip": "Conventional Pearson test against zero correlation; separate from the neutral simulation comparison."})
            if gd_short is not None:
                columns.extend([
                    {"key": f"{key_prefix}.location_{suffix}_r2", "label": f"{gd_short} R\u00b2", "numeric": True, "location_extra": True},
                    {"key": f"{key_prefix}.location_{suffix}_slope", "label": "Slope", "numeric": True, "location_extra": True},
                ])
        columns.extend([
            {"key": f"{key_prefix}.location_{suffix}_{valid_suffix}", "label": "N sim", "numeric": True,
             "tooltip": f"Number of valid neutral subsets used for the {method} comparison."},
            {"key": f"{key_prefix}.location_{suffix}_{status_suffix}", "label": "Sim status", "numeric": False,
             "tooltip": f"Availability of the {method} neutral comparison. 'ok' means the observed correlation is defined and enough valid simulations were available. It does not indicate rejection; use p neutral for that assessment."},
        ])
        groups.append({"group": method, "columns": columns})
    return groups


def _covariate_group(display_name, key_prefix, suffix, gd_key, gd_short, extreme_value_percentiles,
                     tooltip=""):
    """One window Pearson reference and grouped location comparisons per predictor."""
    return {
        "group": display_name,
        "tooltip": tooltip,
        "columns": [
            {
                "group": "Window",
                "columns": [
                    {"key": f"{key_prefix}.window_{gd_key}_{suffix}_r", "label": "r", "numeric": True,
                     "tooltip": "Pearson correlation between reliable window diversity values and the environmental values at their cell centres. Descriptive reference; overlapping windows are not independent observations."},
                ],
            },
            {
                "group": "Location",
                "columns": _location_stat_columns(key_prefix, suffix, extreme_value_percentiles,
                                                  gd_short=gd_short),
            },
        ],
    }


def build_column_layout(gd_key, gd_short, extreme_value_percentiles=LEGACY_EXTREME_VALUE_PERCENTILES):
    """The catalog's column layout for one genetic diversity statistic.

    Only two things depend on the statistic, both taken from the catalog's own
    data (see read_gd_metric below), never hardcoded here:
      gd_key   - lower-cased gd_stat; part of the Window fields' names in
                 gene_summary.json (window_<gd_key>_latitude_r, ...)
      gd_short - short display label ("\u03c0", "Ho", "AR", ...) for the column
                 headers that say which statistic a regression is about
    Everything else (windows.*, location.*, n_individuals_vcf.*, ...) is
    statistic-neutral in the JSON already and identical for every statistic.
    """
    return [
        {"key": "gene_id", "label": "Gene", "numeric": False},
        {"key": "chr", "label": "Chr", "numeric": False},
        {"key": "n_variants", "label": "Variants", "numeric": True},
        {"key": "n_locations_used", "label": "N sites", "numeric": True},
        {
            "group": "N ind",
            "columns": [
                {"key": "n_individuals_vcf.retained", "label": "Retained", "numeric": True},
                {"key": "n_individuals_vcf.initial", "label": "Initial", "numeric": True},
                {"key": "n_individuals_vcf.vcf_error", "label": "VCF Error", "numeric": True},
                {"key": "n_individuals_vcf.coords_mismatch", "label": "Coords Error", "numeric": True},
            ],
        },

        {
            "group": "Windows",
            "columns": [
                {"key": "windows.variance", "label": "Var", "numeric": True},
                {"key": "windows.mean", "label": "Mean", "numeric": True},
                {"key": "windows.abs_mean", "label": "|Mean|", "numeric": True,
                 "tooltip": "Absolute value of the mean across eligible raw windows, abs(mean). Stored only when this gene's observed window min < 0 and max > 0; otherwise unavailable."},
                {"key": "windows.median", "label": "Median", "numeric": True},
                {"key": "windows.min", "label": "Min", "numeric": True},
                {"key": "windows.max", "label": "Max", "numeric": True},

            ],
        },
        {
            "group": "Location",
            "columns": [
                {"key": "location.variance", "label": "Var", "numeric": True},
                {"key": "location.mean", "label": "Mean", "numeric": True},
                {"key": "location.abs_mean", "label": "|Mean|", "numeric": True,
                 "tooltip": "Absolute value of the mean across defined locations, abs(mean). Stored only when this gene's observed location min < 0 and max > 0; otherwise unavailable."},
                {"key": "location.median", "label": "Median", "numeric": True},
                {"key": "location.min", "label": "Min", "numeric": True},
                {"key": "location.max", "label": "Max", "numeric": True},

            ],
        },
        #     {
        #     "group": "Centroid",
        #     "columns": [
        #         {"key": "centroid.variance", "label": "Var", "numeric": True},
        #         {"key": "centroid.mean", "label": "Mean", "numeric": True},
        #         {"key": "centroid.median", "label": "Median", "numeric": True},
        #         {"key": "centroid.min", "label": "Min", "numeric": True},
        #         {"key": "centroid.max", "label": "Max", "numeric": True},

        #     ],
        # },
        # {"key": "pi_variance_per_sample", "label": "Var(\u03c0) sample", "numeric": True},
        #  [°C]
        {
            "group": "Latitude",
            "columns": [
                {"group": "Window", "columns": [
                    {"key": f"latitude.window_{gd_key}_latitude_r", "label": "r", "numeric": True,
                     "tooltip": "Pearson correlation between reliable window diversity values and the latitude of their cell centres. Descriptive reference; overlapping windows are not independent observations."},
                ]},
                {"group": "Location", "columns": _location_stat_columns("latitude", "latitude", extreme_value_percentiles)},
            ],
        },

        _covariate_group("Salinity", "salinity", "salinity", gd_key, gd_short, extreme_value_percentiles),
        _covariate_group("Water Temp Summer", "water_temp_summer", "temp_summer", gd_key, gd_short, extreme_value_percentiles),
        _covariate_group("Water Temp Winter", "water_temp_winter", "temp_winter", gd_key, gd_short, extreme_value_percentiles),
        _covariate_group("log(chlorophyll-a)", "chlorophyll_a_log", "chl_log", gd_key, gd_short, extreme_value_percentiles),
        _covariate_group("Oxygen", "oxygen", "oxygen", gd_key, gd_short, extreme_value_percentiles,
                         "Dissolved oxygen concentration in mmol/m³, averaged over the configured baseline at the configured depth."),
        _covariate_group("Current Speed", "current_speed", "current_speed", gd_key, gd_short, extreme_value_percentiles,
                         "Horizontal current speed in m/s: sqrt(uo² + vo²) is calculated for each monthly field, then averaged over the physics baseline."),

    ]


def _leaf_keys(node):
    """Every leaf column key under a layout node, recursively - a plain
    column ({"key": ...}) is its own single-item list; a group ({"group":
    ..., "columns": [...]}) is the concatenation of its children's. Used
    to find a representative key (leaf_keys(node)[0]) for whatever struct
    a group ultimately belongs to, however deep it's nested."""
    if "group" not in node:
        return [node["key"]]
    keys = []
    for child in node["columns"]:
        keys.extend(_leaf_keys(child))
    return keys


def _nonnull_count(con, table, key):
    """COUNT() of non-null values for one leaf key. A dotted struct path
    (e.g. "salinity.window_pi_salinity_r2") is resolved as quoted
    "field"."subfield" access, which DuckDB resolves natively against a
    STRUCT column - replaces the old _struct_field_names() + one combined
    SUM-of-counts query with a plain per-leaf count, simpler now that
    emptiness is checked per-leaf rather than once for a whole flat group."""
    quoted = ".".join(f'"{part}"' for part in key.split("."))
    try:
        return con.execute(f"SELECT count({quoted}) FROM {table}").fetchone()[0]
    except duckdb.BinderException:
        # The struct exists but this FIELD does not: a catalog whose gene summaries were all written before the field
        # was added to the pipeline (e.g. the r / p columns of WRITE_LOCATION_STATS). Callers treat None like a column
        # missing from the schema - dropped silently, nothing to disclose - instead of crashing.
        return None


def _prune_empty(con, table, node, existing_cols, path=(), computed_counts=None):
    """Recursively drop a column or group with no data anywhere in this
    run, at any nesting depth - the generalization of available_columns's
    two-case logic (see its docstring) to a layout tree of arbitrary depth
    rather than just one level of grouping.

    Returns (kept_node_or_None, dropped_labels):
      kept_node_or_None - `node`, with any empty descendants pruned out of
                           its "columns" list, or None if `node` itself
                           should disappear entirely (schema-absent, or
                           nothing under it survived).
      dropped_labels     - display labels for case 2 drops (see
                           available_columns) found anywhere under `node`,
                           " -> "-joined with `path` once below the top
                           level - e.g. "Salinity -> Location" for that one
                           nested half disappearing, vs. a bare "Latitude"
                           when a whole top-level group does. Empty
                           wherever nothing was dropped, and always empty
                           for case 1 (nothing to disclose).
    """
    if "group" not in node:
        struct_key = node["key"].split(".")[0]
        if struct_key not in existing_cols:
            return None, []  # case 1
        source_key = node.get("source_key", node["key"])
        n_nonnull = (computed_counts[source_key] if computed_counts and source_key in computed_counts
                     else _nonnull_count(con, table, source_key))
        if n_nonnull is None:
            return None, []  # case 1 - field absent from this schema
        if n_nonnull == 0:
            return None, [" \u2192 ".join((*path, node["label"]))]  # case 2
        return node, []

    struct_key = _leaf_keys(node)[0].split(".")[0]
    if struct_key not in existing_cols:
        return None, []  # case 1 - whole struct absent from this schema

    kept_children, dropped = [], []
    for child in node["columns"]:
        kept_child, child_dropped = _prune_empty(
            con, table, child, existing_cols, path=(*path, node["group"]),
            computed_counts=computed_counts
        )
        if kept_child is not None:
            kept_children.append(kept_child)
        dropped.extend(child_dropped)

    if not kept_children:
        # Nothing under this group survived - one combined label for the
        # group itself reads better than repeating every leaf already
        # implied by it (one "Salinity", not "Salinity -> Window -> pi
        # r2", "Salinity -> Window -> slope", ...).
        return None, [" \u2192 ".join((*path, node["group"]))]

    return {**node, "columns": kept_children}, dropped


def available_columns(con, table, layout, existing_cols, computed_counts=None):
    """Drop a plain column or group two different ways, for two different
    reasons:

    1. The key (or, for a group, its whole underlying struct) never
       appears in this run's schema at all (e.g. Centroid, commented out
       above - old-schema genes without "windows" at all). Nothing to
       disclose here: the column was never a real possibility for this
       run's data in the first place.
    2. The key/struct IS in the schema, but every single row is null for
       it - the case a domain_type-gated group (Salinity, Water Temp,
       ...) hits on a run where that covariate genuinely doesn't apply
       (e.g. a land species run), whether the Rmd wrote an explicit
       null-filled struct per gene or just never wrote the key at all -
       both look identical here: zero non-null values across every row.
       This case IS disclosed (see dropped_empty below), since a column
       disappearing because it has real data withheld is a meaningfully
       different situation from case 1, worth a reader knowing about.

    A nested group (Salinity's Window/Location sub-groups, say) is pruned
    per sub-group rather than as one all-or-nothing unit: if Window has
    real data but Location doesn't (Track A/B in context.md run on
    genuinely separate code paths, so this can happen for real, not just
    hypothetically), Salinity is kept with only its Window half, and
    "Salinity -> Location" is disclosed on its own - a strictly more
    informative version of the old flat behaviour, which could only keep
    or drop each group as a whole.

    Returns (kept_layout, dropped_empty, strip_keys):
      kept_layout   - same shape as `layout`, minus anything dropped either
                      way, and with any nested groups themselves pruned
      dropped_empty - display labels for case 2 only, for the page's own
                      one-line disclosure
      strip_keys    - top-level dict keys to delete from every record
                      before they're serialized to JSON, for case 2 only,
                      and only when an ENTIRE top-level group is empty - a
                      partially-empty nested sub-group is dropped from the
                      rendered layout only, without touching stored
                      per-record JSON. A few unused bytes staying embedded
                      in that case isn't worth stripping at a nested path
                      for. Case 1 needs no stripping either: DuckDB never
                      returned that column to begin with, so it was never
                      in any record.
    """
    kept, dropped_empty, strip_keys = [], [], []
    for col in layout:
        if "group" in col:
            kept_node, dropped = _prune_empty(con, table, col, existing_cols,
                                             computed_counts=computed_counts)
            if kept_node is None:
                if dropped:  # case 2: whole top-level group emptied out
                    strip_keys.append(_leaf_keys(col)[0].split(".")[0])
                dropped_empty.extend(dropped)  # stays [] for case 1
                continue
            kept.append(kept_node)
            dropped_empty.extend(dropped)
        else:
            if col["key"] not in existing_cols:
                continue  # case 1
            n_nonnull = _nonnull_count(con, table, col["key"])
            if n_nonnull is None:
                continue  # case 1 - field absent from this schema
            if n_nonnull == 0:
                dropped_empty.append(col["label"])
                strip_keys.append(col["key"])
                continue  # case 2
            kept.append(col)
    return kept, dropped_empty, strip_keys


def embeddable_json(value):
    """json.dumps, but safe to sit inside a <script> tag: escapes '<' so a
    string value (a file path, say) can never contain a literal
    "</script>" that would close the tag early."""
    return json.dumps(value).replace("<", "\\u003c")


def read_skipped_genes(path):
    """Read the optional report; absent or row-free reports add no HTML control."""
    if not path:
        return []
    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            text = handle.read()
    except FileNotFoundError:
        return []
    if not text.strip():
        return []
    reader = csv.reader(io.StringIO(text.lstrip("\r\n")), delimiter="\t")
    columns = (
        "gene_id", "chromosome", "requested_region", "gene_scope", "reason",
        "raw_vcf_records", "retained_vcf_records",
    )
    header = next(reader, [])
    if not set(columns).issubset(header) or len(header) != len(set(header)):
        raise ValueError(f"Invalid skipped-gene report header: {path}")
    records = []
    for values in reader:
        if not any(value.strip() for value in values):
            continue
        if len(values) != len(header):
            raise ValueError(f"Invalid skipped-gene report row {reader.line_num}: {path}")
        row = dict(zip(header, values))
        if not row["gene_id"].strip():
            raise ValueError(f"Missing gene_id in skipped-gene report row {reader.line_num}: {path}")
        records.append(row)
    return records


def gene_plot_paths(gene_id, gene_dir, html_output_path):
    """PNGs for this specific gene, as paths relative to wherever the HTML
    file itself ends up (so the page keeps working if results/ is zipped
    up and shared or moved as a whole).

    gene_dir is the SHARED per-chromosome folder ({results_dir}/results/
    {chrom}/), not a private per-gene one - every gene on that chromosome
    has its plots sitting there together, distinguished only by filename
    (the Rmd's Sys.Date()-stamped, settings-encoded names already embed
    the gene ID). So this can't just list the folder; it has to keep only
    the filenames that mention this particular gene_id. Ensembl-style IDs
    here are fixed-length and zero-padded, so one gene's ID can't be a
    proper substring of another's - a plain substring check is safe.
    """
    if not gene_dir or not os.path.isdir(gene_dir):
        return []
    rel_dir = os.path.relpath(gene_dir, start=os.path.dirname(html_output_path))
    return [
        f"{rel_dir}/{fn}"
        for fn in sorted(os.listdir(gene_dir))
        if fn.lower().endswith(".png") and gene_id in fn
    ]

# What the Rmd falls back to for gd_stat = "pi". Also what a catalog built BEFORE gd_stat existed
# gets, since every such catalog is nucleotide diversity.
LEGACY_GD_METRIC = ("pi", "Nucleotide diversity (\u03c0)", "\u03c0")


def read_gd_metric(con, table, existing_cols):
    """(gd_stat, gd_label, gd_label_short) of the statistic this catalog was built with, read from
    the catalog's own data (the Rmd writes all three into every gene's JSON) instead of from
    config.yaml - so the page always describes what was actually computed, even if config.yaml
    was edited after the run, and this script needs no config plumbing of its own.

    Two genes disagreeing on gd_stat itself is an error (aggregate_gene_summaries.py refuses to
    build such a catalog in the first place, this is the same check for a DuckDB file built some
    other way). Genes that only differ in label wording - someone reworded gd_label in
    config.yaml and re-ran part of the catalog - are not worth stopping for: the most common
    wording wins and a note says so.
    """
    if "gd_stat" not in existing_cols:
        return LEGACY_GD_METRIC
    label_col = "gd_label" if "gd_label" in existing_cols else "NULL"
    short_col = "gd_label_short" if "gd_label_short" in existing_cols else "NULL"
    rows = con.execute(
        f"SELECT gd_stat, {label_col}, {short_col}, count(*) FROM {table} "
        f"GROUP BY 1, 2, 3 ORDER BY count(*) DESC"
    ).fetchall()
    stats = {r[0] for r in rows}
    if len(stats) > 1:
        raise RuntimeError(
            f"{table} mixes genetic diversity statistics ({sorted(map(str, stats))}) - "
            f"one catalog can only hold one."
        )
    stat, label, short, _ = rows[0]
    if len(rows) > 1:
        print(f"  NOTE: {sum(r[3] for r in rows[1:])} gene(s) carry a different gd_label/gd_label_short "
              f"wording than the {rows[0][3]} that win; using {label!r} / {short!r}.")
    # gd_stat is never NULL here; if a label is, fall back to the wingen name rather than print "None"
    return stat, (label or stat), (short or stat)


parser = argparse.ArgumentParser()
parser.add_argument("--duckdb", required=True)
parser.add_argument("--html", required=True)
parser.add_argument("--skipped-genes", help="Optional skipped_genes.tsv to embed in a collapsed report")
args = parser.parse_args()

duckdb_path = args.duckdb
html_path = args.html
skipped_genes = read_skipped_genes(args.skipped_genes)

con = duckdb.connect(duckdb_path, read_only=True)
# con = duckdb.connect(snakemake.input.duckdb, read_only=True)
cols = [r[0] for r in con.execute("DESCRIBE genes").fetchall()]
rows = con.execute("SELECT * FROM genes ORDER BY gene_id").fetchall()

# Needs `con` still open - the emptiness check in available_columns runs
# its own COUNT queries against "genes" rather than working from `rows`
# in Python, so it stays correct even if this script's own row-fetch
# above ever changes shape.
gd_stat, gd_label, gd_label_short = read_gd_metric(con, "genes", cols)
# gd_key has to match how the Rmd builds the Window fields' names: tolower(gd_stat).
extreme_value_percentiles = read_extreme_value_percentiles(rows, cols)
COLUMN_LAYOUT = build_column_layout(gd_stat.lower(), gd_label_short, extreme_value_percentiles)
records = [dict(zip(cols, row)) for row in rows]
window_counts = add_window_correlations(records, gd_stat.lower())
column_layout, dropped_empty, strip_keys = available_columns(con, "genes", COLUMN_LAYOUT, cols,
                                                             computed_counts=window_counts)
con.close()

for rec in records:
    rec["plots"] = gene_plot_paths(
        rec["gene_id"],
        rec.get("gene_dir"),
        html_path,
    )
    # rec["plots"] = gene_plot_paths(rec["gene_id"], rec.get("gene_dir"), snakemake.output.html)
    for key in strip_keys + ["gd_stat", "gd_label", "gd_label_short"]:
        rec.pop(key, None)

# One line disclosing what got hidden and why, only when something did -
# e.g. "4 columns hidden - no data for any gene in this run: Salinity,
# Water Temp Summer, Water Temp Winter, log(chlorophyll-a)." Deliberately
# doesn't try to explain WHY in domain_type terms (that's config.yaml's
# concern, not this script's) - just states the fact plainly enough that
# it reads as "expected, not broken" to whoever's looking at the page.
empty_columns_note = (
    f"{len(dropped_empty)} column{'s' if len(dropped_empty) != 1 else ''} hidden "
    f"- no data for any gene in this run: {', '.join(dropped_empty)}."
    if dropped_empty else ""
)

generated_at = duckdb.execute("SELECT current_localtimestamp()").fetchone()[0]
# str(generated_at) would include microseconds (e.g. "...22:41.905728") -
# nobody reading the page needs sub-second precision on a build timestamp.
generated_at_display = generated_at.strftime("%Y-%m-%d %H:%M")

env = Environment(loader=FileSystemLoader(TEMPLATE_DIR), autoescape=True)
html = env.get_template("gene_catalog_template.html").render(
    generated_at=generated_at_display,
    gd_stat=gd_stat,
    gd_label=gd_label,
    skipped_genes=skipped_genes,
    empty_columns_note=empty_columns_note,
    neutral_thresholds_present=any(
        "_neutral_abs_r_q" in key or "_neutral_abs_rho_q" in key
        for node in column_layout for key in _leaf_keys(node)
    ),
    data_json=embeddable_json(records),
    columns_json=embeddable_json(column_layout),
)

with open(html_path, "w", encoding="utf-8") as f:
    f.write(html)

print(f"wrote {len(records)} rows to {html_path}")

# with open(snakemake.output.html, "w") as f:
#     f.write(html)

# print(f"wrote {len(records)} rows to {snakemake.output.html}")
