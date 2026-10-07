# Config
`config/config.yaml` is the configuration file for the gene-catalog run and includes:
- Definition of the input and output paths

- Wingen specific settings

## Definition of the input and output paths
| Variable | Purpose |
|---|---|
| `results_dir` | Output folder for the gene catalog run. Defaults to `results/gene_catalog`. To avoid two gene catalogs colliding, simply assign a unique name for `results_dir`.|
| `input_vcf`, `sample_coordinates_file` | Paths to the cohort VCF and the per-individual coordinates file |
# Gene selection
Two independent settings decides what the pipeline does with the genes:
* **which genes run:** `gene_list_file` or `genes`
* **how each gene's variants are found in the VCF:** `vcf_annotation_format`.
| Variable | Purpose |
|---|---|


| `gene_list_file`, `genes` | Defines which genes to run. The file type is read from the extension: `.bed` is a BED file (chrom, start, end, gene ID; one row per gene), any other extension is a text file with one gene ID per line. If the path is set but the file doesn't exist yet, `batch_job_gene_catalog_snakemake.sh` creates it from the SnpEff annotation of `input_vcf` (needs a SnpEff-annotated VCF). |

| `genes` | a hardcoded array of genes to run, used when `gene_list_file` is not set. |

| `gene_list_limit` | Run only the first N genes of the list, useful for a quick test run. |

| `vcf_annotation_format` | `snpeff` \| `bed`: How a gene's variants are found in the VCF. `snpeff` matches the gene ID in the `INFO/ANN` field. `bed` extracts the gene's `chrom:start-end` region using the coordinates in the BED gene list, so it requires `gene_list_file` to be a `.bed` file, but it works for VCF files without SnpEff annotation. |

| `gene_chrom_file` | Chromosome of each gene (`chrom<TAB>gene_id` per line, built by `code/build_gene_chrom_map.py`). Only needed for a text gene list or `genes` array combined with `vcf_annotation_format: snpeff`. With a BED gene list the chromosomes are read from the BED and this setting is not used. If the path is set but the file doesn't exist yet, the batch script builds it from the SnpEff annotation. |

| `gene_scope` | `all` \| `body` \| `exonic` \| `cds`: Used if `vcf_annotation_format`:`snpeff`. Runs filter_vcf_by_scope.sh to narrow down the gene-filtered VCF down (based off the ANN block) to filter the scope of the genes down to: `body` (gene body), `exonic` (exonic region), `cds` (coding sequence). Default option is `all`, which includes the whole gene region including upstream and downstream region annotations. |

## Gene list file formats
**BED** (`.bed`, case-insensitive): the genes are listed with stamdard bed format with one gene per row and the following columns:
- column 1-3: chrom, start, end
- column 4: gene ID;
any further columns are ignored currently. 

**Text** (any other extension): one gene ID per line. A line with several columns stops the run.

## Wingen script parameters
These parameter settings are environment variables used in `code/2_Herring_Wingen.Rmd` and can optionally be set, otherwise fallback values are used. More information about the wingen function arguments can be found from the official documentation: https://cran.r-project.org/web/packages/wingen/refman/wingen.html
### Domain
`domain_type`: `"ocean"` (default) or `"land"`: intended to switch which side of the coastline is kept (marine vs. terrestrial species).
**Note**: Not fully implemented into `2_Herring_Wingen.Rmd` yet.
### Booleans
`calc_gen_diversity_stats`: decides whether a statistical summary file is created, `{gene}_gene_summary.json` (`TRUE`). Required to be set `TRUE` for creating the HTML gene catalog.

`use_copernicus_marine`: decides whether the whole environmental-covariates section (salinity, water temperature, chlorophyll, oxygen and current speed) is executed (default `FALSE`). See its own section below.
### Naming notation
`study_species`: free-text label used in the output filenames and plot titles.
### Genetic diversity statistic
`gd_stat`: `"pi"` (default) \| `"Ho"` \| `"biallelic_richness"`. Determines which statistic `window_gd()` maps, carried through everything downstream: window stats, per-location stats, plot titles, `gene_catalog.html`'s column headers, and output filenames.

`gd_label` / `gd_label_short`: optional text overrides for how the chosen statistic is displayed. `gd_label` is the long form (plot titles, the `gene_catalog.html` heading); `gd_label_short` is the short form (Legends, column headers). If not set, following defaults are used:

| `gd_stat` | `gd_label` (default) | `gd_label_short` (default) |
|---|---|---|
| `pi` | Nucleotide diversity (π) | π |
| `Ho` | Observed heterozygosity (Ho) | Ho |
| `biallelic_richness` | Allelic richness (biallelic) | AR |

**Worth knowing:**

- **One `results_dir` can only hold one statistic.** A re-run under a different statistic would overwrite the JSONs but leave the old statistic's PNGs sitting next to the new ones in the same per-chromosome folders. Use a separate `results_dir` per statistic instead.

### Geographical Area Settings
To subset your dataset into a "boxed" region to restrict which *individuals* are included, the allowed interval of latitude and longitude coordinates can be defined from `lat_min`/`lat_max`/`lon_min`/`lon_max` (WGS84 degrees). All four must be set together for this filter to activate. 
#### Raster construction
`RESOLUTION_KM`: cell size of the working raster, in km (pipeline default `2.5`). Tune this against your sample density via the Rmd's own `preview_gd()` step.

`CTR_DISAGG`: disaggregation factor applied to the raster at build time (subdivides each cell). Leave unset to disable (default).

## Rolling Window Computation - `window_gd()` (or `window_general()`)
Which statistic this section's settings tune is controlled by `gd_stat` above (`pi` by default). (`window_general()`, wingen's dataframe-based alternative, isn't wired into this pipeline at all yet — see `README.md`.)

`WDIM` (`wdim` in `window_gd()`): window dimensions in aggregated cells.
`FACT` (`fact` in `window_gd()`): aggregation factor applied to the raster before windowing.
`MIN_N_SAMPLES_PER_CALC` (`min_n` in `window_gd()`): minimum number of empirical samples a window needs; fewer than this and the focal cell gets `NA`. Default: `2`.

`use_rarify`: equivalent to rarify in window_gd()
`rarify_n`: equivalent to rarify_n in window_gd(), only used if rarify = TRUE. Defines the number of points to use for rarefaction (defaults to min_n). [Default: 2]
`RARIFY_N_ITERATIONS` (`rarify_nit` in `window_gd()`, only used if `rarify = TRUE`): number of iterations to use for rarefaction (defaults to 5). Can also be set to `"all"` to use every combination of samples of size `RARIFY_N` within the window.

`LOCATION_RARIFY`: `auto` (default) | `TRUE` | `FALSE`. Whether location diversity uses repeated draws of `rarify_n` individuals. `auto` uses individual subsampling only for `biallelic_richness` when `USE_RARIFY` is enabled, while Pi and Ho instead uses all eligible individuals in a cluster. 

Location statistics first group distinct sampling sites by complete-linkage
clustering with `CLUSTER_VARIANCE_POINTS_KM` (default 10 km), then recompute
diversity from each cluster's pooled genotypes. Environmental values are then averaged over the members of each cluster and each cluster contributes one equally weighted diversity/environment pair to the correlation, irrespective of its individual count.

## Kriging Settings (Interpolation) - `wkrig_gd()`
`GRD_DISAGG_FACTOR`: parameter previously existing in the old `krig_gd(disagg_grd)` to disaggregate the prediction grid (smoothing) prior to kriging. Here, `terra`'s `disagg` function is used instead since `wkrig_gd()` lacks this functionality. This parameter subdivides each lyr cell into `grd_disagg_factor^2` sub-cells for smoothing (finer pixel spacing) in the kriging prediction grid. Setting this parameter to 1 keeps grd stays at lyr's native resolution.
`USE_WEIGHTED_KRIGING`: Boolean determining whether weighted kriging should be performed (currently weights are assigned based on sample count: `wkrig_gd(weight_r = wgd_confident[["sample_count"]])`).

`NMAX_NEIGHBOURS`(`nmax` in `wkrig_gd()`): decides the maximum number of neighboring observations to use for kriging at each prediction location.
`MAXDIST` (`maxdist` in `wkrig_gd()`): maximum distance to a neighboring point to still consider it. Left unset by default and resolved automatically from the computed `MAX_INTERP_DIST` below (the smallest value that won't exclude a centroid a surviving prediction could need). Set explicitly only to override that.

`MAX_RANGE_FRAC`: equivalent to `max_range_frac` in `wkrig_gd()` — caps how much of the fitted variogram's range parameter is used when searching for neighbors. Default `0.5`. Live, not overridden.

## Masking
### Sample-count masking
`MIN_SAMPLE_COUNT_MODE`: controls whether automatic masking based on empirical sample count is applied. With `"hardcoded"`, interpolated cells are masked when the corresponding spatial region contains fewer than `MIN_SAMPLE_COUNT` empirical samples. Setting this to `"none"` disables sample-count-based masking.

`MIN_SAMPLE_COUNT`: minimum number of empirical samples required for a spatial region to remain unmasked when `MIN_SAMPLE_COUNT_MODE="hardcoded"`.

### Interpolation-distance masking
**Note:** separate from sample-count masking above. This masking is purely based on the distance from empirical data, regardless of how many samples are nearby.
`MAX_INTERP_DIST_MODE`:
- `"computed"` (default): a window-footprint-based reach, ignoring real sample positions.
- `"clustered"`: median nearest-neighbor distance between clustered real sampling sites; more conservative.
- `"hardcoded"`: use the numerical value computed from `MAX_INTERP_DIST_KM`.
- `"none"`: no distance cutoff at all.

`MAX_INTERP_DIST_KM`: required only if `MAX_INTERP_DIST_MODE: "hardcoded"`.

`MAX_INTERP_DIST_SCALE`: A factor which controls how far interpolation may extend from observed data whichever `MAX_INTERP_DIST` mode resolves to. Controls how far the interpolation may extend beyond observed data (default `1`).

(`MAX_INTERP_DIST_SCALE`), which controls how far interpolation may extend from observed data.
`MIN_SAMPLE_COUNT_VALUE`: Minimum number of empirical samples required for a spatial region to remain unmasked when `MIN_SAMPLE_COUNT_MODE="hardcoded"`.

### Species-range polygon masking
`MASKING_POLYGON`: path to a raster mask — typically a species distribution or habitat layer (e.g. a spawning-area or range `.tif`) — reprojected and aligned onto the working grid, then applied via `mask_gd()`. The practical use: cropping the interpolation to where the species actually occurs, rather than letting it extend across open water (or land) the species doesn't naturally reach. If set, the file must exist or the render stops immediately; if unset, this step is skipped entirely. Two real example layers (both HELCOM/HOLAS products — a herring spawning-area layer and a *Zostera marina* distribution layer) sit commented out right next to the setting in the Rmd already.


## Plot Display Settings
`SHOW_SAMPLE_POINTS`: controls whether the empirical sample locations are overlaid on the final masked π map. Set to `TRUE` to display the sample points; set to `FALSE` to omit them.

`SAMPLE_POINT_SIZE_BY`: Only used when `SHOW_SAMPLE_POINTS=TRUE` and determines which statistic controls the size of the overlaid sample points. `"count"` scales point size by the number of individuals sharing each location (or sample-point cluster), while `"value"` scales point size by the empirical genetic diversity metric value at that location. `"All"` (default) generates both versions as separate output files.

`HEATMAP_COLOR_PALETTE`: selects the color palette used for the heatmaps. The available options correspond to the palettes supported by `scale_fill_viridis_c()`: `"inferno"` (wingens default), `"mako"`, `"rocket"`, `"turbo"`, `"cividis"`, `"plasma"`, `"viridis"`, and `"magma"`. This setting only affects the genetic diversity layer heatmaps; the sample-count layer uses its own fixed `"viridis"` palette.

`SPHERICAL_GRID`: draws the graticule with a curved (Lambert Azimuthal Equal-Area) projection instead of plain lat/long, to visually emphasize the map's 2D-projection-of-a-3D-surface nature. Display-only. Default `FALSE`.

`CLUSTER_SAMPLE_POINTS` / `CLUSTER_SAMPLE_POINTS_KM`: groups sample dots within `CLUSTER_SAMPLE_POINTS_KM` of each other into one plotted dot, purely for map readability (default `FALSE` / `2`). A separate computation from the sample-clustering used to derive `MAX_INTERP_DIST`'s `"clustered"` mode above; changing one never silently changes the other.

## Other
`CLUSTER_VARIANCE_POINTS_KM`: merge distance (km) for the gene-summary chunk's per-location clustering specifically. A separate "knob" from `CLUSTER_SAMPLE_POINTS_KM` above, which only affects the map plot. Default `10`.

## Environmental covariates (CopernicusMarine)
The following settings are only relevant if `use_copernicus_marine: TRUE`.

`ENV_PRESET`: which Copernicus Marine datasets the environmental values (salinity, winter and summer temperature, log chlorophyll-a, dissolved oxygen and current speed) come from. `global` (GLORYS12V1 reanalysis, ~8 km + GlobColour ocean colour, 4 km; the prototype's data) or `baltic` (Baltic reanalysis, ~1.9 km + Baltic ocean colour with a regional chlorophyll algorithm, 1 km). [Default if unset: `global`]

The presets themselves (product id, dataset/layer id, variable names, the years each dataset covers, and notes on when to use which) live in **`config/env_presets.yaml`**. To use another dataset, add an entry there and set `ENV_PRESET` to its name. The first run with a new preset downloads the layers (cached in `{results_dir}/temp/env_cache`); test it on one gene first (`gene_list_limit: 1`).

`ENV_OXYGEN_ENABLED`: include dissolved oxygen in the downloads, per-location values, window Pearson correlation, location Pearson/Spearman statistics and neutral comparisons. Set `FALSE` to skip it. [Default if unset: FALSE; `config.yaml` uses TRUE]

`ENV_CURRENT_SPEED_ENABLED`: include current speed in the same outputs and analyses. Set `FALSE` to skip it. [Default if unset: FALSE; `config.yaml` uses TRUE]

For `baltic`, salinity, temperature and current speed components come from `BALTICSEA_MULTIYEAR_PHY_003_011`, dataset `cmems_mod_bal_phy_my_P1M-m`. Oxygen comes from the separate `BALTICSEA_MULTIYEAR_BGC_003_012` product, dataset `cmems_mod_bal_bgc_my_P1M-m`. Both are documented in the [Baltic product manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-BAL-PUM-003-011-012.pdf). The `global` oxygen dataset is documented in the [global biogeochemical product manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-GLO-PUM-001-029.pdf).

`ENV_BASELINE_START_YEAR`, `ENV_BASELINE_END_YEAR`: the years the environmental values are averaged over. Trimmed automatically to what each dataset covers (announced in the log); the years actually used are written to the gene's JSON (`run_settings`) and are part of the cache file names. [Default if unset: 2019 / 2022; `config.yaml` uses 1998 / 2022]

`ENV_DEPTH_M`: depth in metres below the surface for salinity, temperature and current speed. Oxygen also uses this depth unless overridden. `0` requests the top model layer. [Default: 0]

`ENV_OXYGEN_DEPTH_M`: optional oxygen-specific depth in metres below the surface. Omit it to use `ENV_DEPTH_M`; it selects a model layer at that depth, not the sea floor. [Default: `ENV_DEPTH_M`]

`ENV_CACHE_DIR`: where downloaded/derived environmental rasters are cached, shared across every gene run.

## Reproducibility and optional outputs
`SEED`: seed for the random individual draws of rarefaction (`window_gd()`, and the location values when `LOCATION_RARIFY` rarefies them). With a seed the window values are identical on every re-run. [Default: 42]

`WRITE_PER_LOCATION_TSV`: also write `{gene}_per_location.tsv` next to the plots: one row per location with its individuals count, centroid, the diversity statistic (`gd`) and the enabled environmental values - the exact numbers the location-level fits were computed from. [Default: TRUE]

`WRITE_LOCATION_STATS`: add `n`, signed `r`, an analytic 95% confidence interval of `r` (Fisher z-transform, no bootstrapping), conventional Pearson `p` and Spearman `rho` to stored location statistics. CI endpoints are stored in JSON/TSV but are not displayed in the HTML catalog. [Default: TRUE]

`WRITE_RUN_SETTINGS`: write the settings that produced the gene's numbers (seed, rarefaction, window and kriging settings, environmental preset and years) to the JSON as `run_settings`. [Default: TRUE]

The HTML catalog's **Window, Location & run settings** info panel reads these
saved values, including environmental product IDs and requested depths. 

## Neutral-marker comparison

Set `run_neutral_sim: true` and provide `markers_file` or `neutral_vcf`.
`markers_file` can be omitted when `neutral_vcf` is supplied.
This toggle controls both the neutral Snakemake rules and the Rmd's
`RUN_NEUTRAL_SIM` boolean. All settings below are top-level `config.yaml`
entries. 
### Input and sampling parameters

| Parameter | Default | Purpose |
|---|---|---|
| `run_neutral_sim` | `false` | Enable neutral preparation and comparison. The observed gene analysis runs with either setting. Requires `calc_gen_diversity_stats: TRUE`. |
| `markers_file` | | File path to a XLSX/CSV/TSV table or bed file containing neutral marker annotation. Required only when simulations are enabled and `neutral_vcf` is absent. |
| `neutral_vcf` | `null` | Path to an already neutral marker filtered vcf-file. If defined, `markers_file` is not considered. |
| `input_format` | `auto` | `auto`, `table`, `positions` or `bed`. `auto` uses BED for `.bed`/`.bed.gz`, otherwise a table. Only tables run the table converter; direct position/BED inputs bypass it. |
| `technical_replicates` | `100` | Number of random neutral sets per gene. Each set has the target gene's final SNP count. A gene with 40 SNPs produces 100 VCFs of 40 SNPs each. |
| `seed` | `42` | Base seed for reproducible marker sampling. This is separate from `wingen_params.SEED`, which controls individual rarefaction. |
| `min_valid_fraction` | `0.9` | Minimum fraction of requested sets yielding finite correlations, checked separately for Pearson and Spearman for each predictor. With 100 sets, at least 90 are required; otherwise that empirical comparison is unavailable. Statistical failures and job errors differ; see below. |
| `extreme_value_percentiles` | `[95, 99]` | Cutoff percentiles in the neutral model correlation distributionions (Pearson and Spearman). The values must be between 0 and 100. Examples: `[95]`, `[90, 97.5]`. |

### Table parameters

These parameters apply to XLSX/CSV/TSV tables; they are ignored for ready
`positions` or `bed` inputs and when `neutral_vcf` is supplied.

| Parameter | Default | Purpose |
|---|---|---|
| `sheet` | `Sheet1` | XLSX worksheet; ignored for TSV/CSV. |
| `header_row` | `auto` | Definition of the header row. Can either be detected automatically or supplied as a value i.e: `3`. |
| `chrom_column` | `1` | Chromosome column, given as either a column index, or a case-insensitive header name. |
| `pos_column` | `2` | Position column, given as either a column index, or a case-insensitive header name. |
| `marker_label_column` | `TYPE` | The column containing marker annotations, given as either a column index, or a case-insensitive header name. Ignored when `neutral_only: true`. |
| `neutral_marker_label` | `NEUTRAL` | The annotation text in `marker_label_column` that identifies rows to retain, matched case-insensitively. For example, use `"Neutral"` for cells labelled Neutral. Ignored when `neutral_only: true`. |
| `neutral_only` | `false` | `true` declares that every data row in this table is already neutral. The converter keeps all coordinates and skips annotation filtering, so no marker-label column is needed. Ready `positions` and `bed` files already have this requirement and ignore this setting. |

Automatic header detection accepts `Chr`, `CHROM` and `Chromosome` for the
chromosome column and `Pos`, `Position`, `Pos 1`, `pos1` and `pos_1` for position,
case-insensitively, including a leading `#`. It finds the header in the configured
coordinate columns. If no header or multiple possible headers are found, set
`header_row` explicitly. 

Optional:
A Copernicus Marine Service account, with `COPERNICUSMARINE_SERVICE_USERNAME`/`COPERNICUSMARINE_SERVICE_PASSWORD` set as OS-level environment variables.