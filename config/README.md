# Config
`config/config.yaml` is the configuration file for the gene-catalog run and includes:
- Definition of the input and output paths

- Wingen specific settings

## Definition of the input and output paths
| Variable | Purpose |
|---|---|
| `results_dir` | Output folder for the gene catalog run. Defaults to `results/gene_catalog`. To avoid two gene catalogs colliding, simply assign a unique name for `results_dir`.|
| `input_vcf`, `sample_coordinates_file` | Paths to the cohort VCF and the per-individual coordinates file |
## Gene selection
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
**BED** (`.bed`, case-insensitive): the genes are listed with standard bed format with one gene per row and the following columns:
- column 1-3: chrom, start, end
- column 4: gene ID;
any further columns are ignored currently. 

**Text** (any other extension): one gene ID per line. A line with several columns stops the run.

## First-gene warm-up
Every report reads and fills the same caches under `{results_dir}/temp` (the environmental layers in `env_cache`, the Natural Earth maps in `ne_cache`). The first report on an empty cache downloads and builds them, and it is by far the heaviest and slowest report of a run (measured on herring, 125 individuals: about 43 GiB peak and 44 min, against about 4 GiB and 11 min once the layers are cached). 

With `first_gene_warmup: true` (the default, also when the key is missing) the workflow runs **the first retained gene of the list alone, with the whole allocation, and starts every other gene after it**: that gene's extraction, neutral preparation and report finish first; only then are the other genes planned and run together. The wait is on a finished report and not on a particular download, so it covers any environmental data source and any `domain_type` that caches under `results_dir/temp`; the Rmd is not involved. The reference row (`run_neutral_reference`) is not affected: it already waits for every ordinary report.

The mechanism is a second checkpoint, `first_gene_warmup_done`, whose only input is the first gene's `.done` marker. The final aggregation (and the reference row) list all reports only after that checkpoint has run, so Snakemake creates the other genes' jobs once the first report is finished, and with exactly the inputs they have without the warm-up: which reports are out of date is decided as it always was. The checkpoint writes `{results_dir}/temp/first_gene_warmup.json` (the first report's path and finish time). `rule all` asks for it (when at least two genes are retained) and only the aggregation reads it, so re-running the first gene alone re-aggregates the catalog.

| Setting | Behavior |
|---|---|
| `first_gene_warmup: true` | `true`/`false` only. `false`: all genes start together, with no checkpoint and no waiting for the first gene. On the command line Snakemake reads only `True`/`False`: `--config first_gene_warmup=False`. |
| `first_gene_threads: auto` | Snakemake's `--cores` budget, which is also the number of R workers (`future::multisession`) of that report. A positive integer is accepted and limited by `--cores`. R workers are separate processes, so more threads can mean more memory; the peak of the first report with many workers has not been measured yet, so lower the number if it runs out of memory. With an executor that submits one Slurm job per rule, `--cores`/`--jobs` is not the size of a node: give a number. |
| `first_gene_mem_mb: auto` | Smaller declared `mem_mb` (`--resources`) / `SLURM_MEM_PER_NODE` budget minus `first_gene_headroom_mb`; 40000 MB when no budget is declared. An explicit number is limited to the same value (budget minus headroom). There is no retry multiplier. |
| `first_gene_headroom_mb` | Memory kept back from the declared budget for the scheduler and everything else in the job (default 4096); it applies to `auto` and to the limit of an explicit number. |
| `first_gene_runtime_minutes` | Scheduling metadata, default 240. The enclosing Slurm allocation sets the real time limit. |

## Wingen script parameters
These parameter settings are environment variables used in `code/2_Herring_Wingen.Rmd` and can optionally be set, otherwise fallback values are used. More information about the wingen function arguments can be found from the official documentation: https://cran.r-project.org/web/packages/wingen/refman/wingen.html
### Domain
`domain_type`: `"ocean"` (default) or `"land"`: intended to switch which side of the coastline is kept (marine vs. terrestrial species).
The selected domain determines which side of the coastline is retained for prediction eligibility and the final domain mask. The water-route distance (`WATER_DIST_*` below) uses it too: `"land"` routes over land instead of through water.
### Booleans
`calc_gen_diversity_stats`: decides whether a statistical summary file is created, `{gene}_gene_summary.json` (`TRUE`). Required to be set `TRUE` for creating the HTML gene catalog.

`use_copernicus_marine`: decides whether the whole environmental-covariates section (salinity, water temperature, chlorophyll, oxygen and current speed) is executed (default `FALSE`). See its own section below.
### Naming notation
`study_species`: free-text label used in the output filenames and plot titles.
### Genetic diversity statistic
`gd_stat`: `"pi"` (default) \| `"Ho"` \| `"biallelic_richness"`\|`"theta"` (hierfstat::theta.Watt.dosage)\| `"tajimas_d"` (hierfstat::TajimaD.dosage). Determines which statistic `window_gd()` maps, carried through everything downstream: window stats, per-location stats, plot titles, `gene_catalog.html`'s column headers, and output filenames.

additional support has been added for 
`gd_label` / `gd_label_short`: optional text overrides for how the chosen statistic is displayed. `gd_label` is the long form (plot titles, the `gene_catalog.html` heading); `gd_label_short` is the short form (Legends, column headers). If not set, following defaults are used:

| `gd_stat` | `gd_label` (default) | `gd_label_short` (default) |
|---|---|---|
| `pi` | Nucleotide diversity (π) | π |
| `Ho` | Observed heterozygosity (Ho) | Ho |
| `biallelic_richness` | Allelic richness (biallelic) | AR |
| `theta` | Wattersons Theta (θ) | θ |
| `tajimas_d` | Tajimas D (D) | D |

**Worth knowing:**

- **One `results_dir` can only hold one statistic.** A re-run under a different statistic would overwrite the JSONs but leave the old statistic's PNGs sitting next to the new ones in the same per-chromosome folders. Use a separate `results_dir` per statistic instead.

### Geographical Area Settings
To subset your dataset into a "boxed" region to restrict which *individuals* are included, the allowed interval of latitude and longitude coordinates can be defined from `lat_min`/`lat_max`/`lon_min`/`lon_max` (WGS84 degrees). All four must be set together for this filter to activate. 
#### Raster construction
`RESOLUTION_KM`: cell size of the working raster, in km (Rmd fallback `2.5`; supplied global test config `10`). Tune this against your sample density via the Rmd's own `preview_gd()` step.

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

`MAXDIST` (`maxdist` in `wkrig_gd()`): maximum distance to a kriging observation. When omitted, the inherited policy uses `MAX_INTERP_DIST + half_window_width`. This is an approximate search policy: square-window corners extend farther and Mollweide planar distances differ from geographic distances. Finite `MAXDIST` or `NMAX_NEIGHBOURS` can change predictions; neither is a guaranteed no-loss speed control.

`MAX_RANGE_FRAC`: equivalent to `max_range_frac` in `wkrig_gd()` — caps how much of the fitted variogram's range parameter is used when searching for neighbors. Default `0.5`. Live, not overridden.

## Masking
### Sample-count masking
`MIN_SAMPLE_COUNT_MODE`: controls whether automatic masking based on empirical sample count is applied. With `"hardcoded"`, interpolated cells are masked when the corresponding spatial region contains fewer than `MIN_SAMPLE_COUNT` empirical samples. Setting this to `"none"` disables sample-count-based masking.

`MIN_SAMPLE_COUNT`: minimum number of empirical samples required for a spatial region to remain unmasked when `MIN_SAMPLE_COUNT_MODE="hardcoded"`.

### Interpolation-distance masking
**Note:** separate from sample-count masking above. This masking is purely based on the distance from empirical data, regardless of how many samples are nearby.
`MAX_INTERP_DIST_MODE`:
- `"computed"` (default): a window-footprint-based reach, ignoring real sample positions.
- `"clustered"`: median of each spacing-cluster centroid's distance to its nearest other centroid. Exact duplicate sites are removed before the separate `SPACING_CLUSTER_KM` cut (default 2 km). This reach can be smaller or larger than computed mode and does not validate interpolation across unsampled gaps.
- `"hardcoded"`: use the numerical value computed from `MAX_INTERP_DIST_KM`.
- `"none"`: no distance cutoff at all.

`MAX_INTERP_DIST_KM`: required only if `MAX_INTERP_DIST_MODE: "hardcoded"`.

`MAX_INTERP_DIST_SCALE`: A factor which controls how far interpolation may extend from observed data whichever `MAX_INTERP_DIST` mode resolves to. Controls how far the interpolation may extend beyond observed data (default `1`).


### Species-range polygon masking
`MASKING_POLYGON`: path to a raster mask — typically a species distribution or habitat layer (e.g. a spawning-area or range `.tif`) — reprojected and aligned onto the working grid, then applied via `mask_gd()`. The practical use: cropping the interpolation to where the species actually occurs, rather than letting it extend across open water (or land) the species doesn't naturally reach. If set, the file must exist or the render stops immediately; if unset, this step is skipped entirely. Two real example layers (both HELCOM/HOLAS products — a herring spawning-area layer and a *Zostera marina* distribution layer) sit commented out right next to the setting in the Rmd already.


## Plot Display Settings
`SHOW_SAMPLE_POINTS`: controls whether the empirical sample locations are overlaid on the final masked π map. Set to `TRUE` to display the sample points; set to `FALSE` to omit them.

`SAMPLE_POINT_SIZE_BY`: Only used when `SHOW_SAMPLE_POINTS=TRUE` and determines which statistic controls the size of the overlaid sample points. `"count"` scales point size by the number of individuals sharing each location (or sample-point cluster), while `"value"` scales point size by the empirical genetic diversity metric value at that location. `"All"` (default) generates both versions as separate output files.

`HEATMAP_COLOR_PALETTE`: selects the color palette used for the heatmaps. The available options correspond to the palettes supported by `scale_fill_viridis_c()`: `"inferno"` (wingens default), `"mako"`, `"rocket"`, `"turbo"`, `"cividis"`, `"plasma"`, `"viridis"`, and `"magma"`. This setting only affects the genetic diversity layer heatmaps; the sample-count layer uses its own fixed `"viridis"` palette.

`SPHERICAL_GRID`: draws the graticule with a curved (Lambert Azimuthal Equal-Area) projection instead of plain lat/long, to visually emphasize the map's 2D-projection-of-a-3D-surface nature. Display-only. Default `FALSE`.

`CLUSTER_SAMPLE_POINTS_KM`: display-only cut of the shared complete-linkage site tree. Omitted/null inherits the resolved `CLUSTER_VARIANCE_POINTS_KM` threshold; zero uses exact numeric coordinate pairs; a positive value cuts the same tree at that distance in km. The superseded `CLUSTER_SAMPLE_POINTS` Boolean is not used. This does not change genotype groups, exact empirical mask sites or the separate `SPACING_CLUSTER_KM` grouping.

## Other
`CLUSTER_VARIANCE_POINTS_KM`: complete-linkage cluster distance for genotype-based location statistics (default 10 km). Display dots inherit this resolved threshold unless `CLUSTER_SAMPLE_POINTS_KM` overrides it. Dot positions average unique member sites; counts retain all individuals, and dot values average extracted window-raster values rather than the pooled-genotype location statistic. Nonfinite mean window values are omitted from both count and value dot plots.

## Environmental covariates (CopernicusMarine)
The following settings are only relevant if `use_copernicus_marine: TRUE`.

`ENV_PRESET`: which Copernicus Marine datasets the environmental values (salinity, winter and summer temperature, log chlorophyll-a, dissolved oxygen and current speed) come from. `global` (GLORYS12V1 reanalysis, ~8 km + GlobColour ocean colour, 4 km; the prototype's data) or `baltic` (Baltic reanalysis, ~1.9 km + Baltic ocean colour with a regional chlorophyll algorithm, 1 km). [Default if unset: `global`]

The presets themselves (product id, dataset/layer id, variable names, the years each dataset covers, and notes on when to use which) live in **`config/env_presets.yaml`**. To use another dataset, add an entry there and set `ENV_PRESET` to its name. The first run with a new preset downloads the layers (cached in `{results_dir}/temp/env_cache`); test it on one gene first (`gene_list_limit: 1`). With the default `first_gene_warmup: true` the workflow does this by itself: the first gene runs alone and fills the cache before the others start.

`ENV_OXYGEN_ENABLED`: include dissolved oxygen in the downloads, per-location values, window Pearson correlation, location Pearson/Spearman statistics and neutral comparisons. Set `FALSE` to skip it. [Default if unset: FALSE; `config.yaml` uses TRUE]

`ENV_CURRENT_SPEED_ENABLED`: include current speed in the same outputs and analyses. Set `FALSE` to skip it. [Default if unset: FALSE; `config.yaml` uses TRUE]

For `baltic`, salinity, temperature and current speed components come from `BALTICSEA_MULTIYEAR_PHY_003_011`, dataset `cmems_mod_bal_phy_my_P1M-m`. Oxygen comes from the separate `BALTICSEA_MULTIYEAR_BGC_003_012` product, dataset `cmems_mod_bal_bgc_my_P1M-m`. Both are documented in the [Baltic product manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-BAL-PUM-003-011-012.pdf). The `global` oxygen dataset is documented in the [global biogeochemical product manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-GLO-PUM-001-029.pdf).

`ENV_BASELINE_START_YEAR`, `ENV_BASELINE_END_YEAR`: the years the environmental values are averaged over. Trimmed automatically to what each dataset covers (announced in the log); the years actually used are written to the gene's JSON (`run_settings`) and are part of the cache file names. [Default if unset: 2019 / 2022; `config.yaml` uses 1998 / 2022]

`ENV_DEPTH_M`: depth in metres below the surface for salinity, temperature and current speed. Oxygen also uses this depth unless overridden. `0` requests the top model layer. [Default: 0]

`ENV_OXYGEN_DEPTH_M`: optional oxygen-specific depth in metres below the surface. Omit it to use `ENV_DEPTH_M`; it selects a model layer at that depth, not the sea floor. [Default: `ENV_DEPTH_M`]

`ENV_CACHE_DIR`: where downloaded/derived environmental rasters are cached, shared across every gene run.

## Water-route distance (optional)
The distance from each sampling location to the nearest start point or start box, travelling **only through water** (`domain_type: "ocean"`) or **only over land** (`"land"`), **never** *as the crow flies*. It is treated like the other environmental variables at the location level (cluster mean, location Pearson/Spearman statistics, neutral comparison, catalog column "Water Route Distance"); there is no window-level correlation. Distances are computed once per dataset and cached in `ENV_CACHE_DIR` (`water_distance_<hash>.tsv`, `.json` and `_routes.png`). Which start to use is a methods choice, e.g. the HELCOM entrance line vs. the three Danish straits for the Baltic. Method, cache columns and test results: `docs/water_distance.md`.

| Setting | Default | Purpose |
|---|---|---|
| `WATER_DIST_ENABLED` | `FALSE` | `TRUE` switches it on and requires at least one start. The Herring SECD block in `config.yaml` is a working example. |
| `WATER_DIST_REFERENCE_POINTS` | none | Start points as `"lon,lat"` pairs separated by semicolons, e.g. `"12.80,55.75; 10.95,55.45"` (Øresund; Great Belt). A single point can also be given as `WATER_DIST_REF_LON` and `WATER_DIST_REF_LAT`. |
| `WATER_DIST_REFERENCE_BOXES` | none | Start boxes separated by semicolons, written as `"latitude = [57.7395, 57.7415] longitude = [10.60, 12.00]"` (the HELCOM entrance to the Baltic) or `"lon_min, lon_max, lat_min, lat_max"`. Every passable cell a box touches is a start at distance 0. |
| `WATER_DIST_RES_KM` | `1` | Cell size of the route search. Straits narrower than a cell close: 1 km keeps Øresund open, 2 km or coarser can close it. |
| `WATER_DIST_MIN_COVER` | `0.99` | Fraction of a cell that must lie in the medium for it to be passable. `0.99` lets any land block a cell; lower values keep narrow straits open but can leak across thin land. At 2 km, `0.9` was the reliable value. |
| `WATER_DIST_SNAP_TOL_KM` | `max(3 × WATER_DIST_RES_KM, 3)` | How far a site may be moved to the nearest passable cell; further than this it fails as `no_passable_cell_near`. |
| `WATER_DIST_NEGATIVE_BOX` | none | One box whose sites get negative distances (e.g. the Kattegat/Skagerrak side of Øresund), so the gradient runs through zero. Leave unset for a study inside one sea. |
| `WATER_DIST_ON_FAILED_SITES` | `stop` | `stop` lists every site that could not be routed and stops the run; `drop` leaves those sites out and continues. |
| `WATER_DIST_QA_PLOT` | `TRUE` | Writes a map of the routes next to the cache file; check that the routes stay in the water. |

`{gene}_per_location.tsv` also gets `water_dist_km` (the cluster's mean distance), `water_dist_range_km` and `water_dist_span_flag`, which flags a cluster whose sites' distances differ by more than 1.5 × `CLUSTER_VARIANCE_POINTS_KM`.

## Reproducibility and optional outputs
`SEED`: seed for the random individual draws of rarefaction (`window_gd()`, and the location values when `LOCATION_RARIFY` rarefies them). With a seed the window values are identical on every re-run. [Default: 42]

`WRITE_PER_LOCATION_TSV`: also write `{gene}_per_location.tsv` next to the plots: one row per location with its individuals count, centroid, the diversity statistic (`gd`) and the enabled environmental values - the exact numbers the location-level fits were computed from. [Default: TRUE]

`WRITE_LOCATION_STATS`: add `n`, signed `r`, an analytic 95% confidence interval of `r` (Fisher z-transform, no bootstrapping), conventional Pearson `p` and Spearman `rho` to stored location statistics. CI endpoints are stored in JSON/TSV but are not displayed in the HTML catalog. [Default: TRUE]

`WRITE_RUN_SETTINGS`: write the settings that produced the gene's numbers (seed, rarefaction, window and kriging settings, environmental preset and years) to the JSON as `run_settings`. [Default: TRUE]

The HTML catalog's **Window, Location & run settings** info panel reads these
saved values, including environmental product IDs and requested depths. 

## Neutral-marker comparison

Set `run_neutral_sim: true` and provide `markers_file` or `neutral_vcf`. A third option, if the user has a SnpEff-annotated VCF-file, is to automatically create a neutral marker file by setting `vcf_annotation_format: snpeff` and leave both neutral inputs null/unset to extract intergenic SNP coordinates automatically from `input_vcf`.
Input precedence is `neutral_vcf`, then `markers_file`, then automatic extraction.
This toggle controls both the neutral Snakemake rules and the Rmd's `RUN_NEUTRAL_SIM` boolean. All settings below are top-level `config.yaml` entries. 
### Input and sampling parameters

| Parameter | Default | Purpose |
|---|---|---|
| `run_neutral_sim` | `false` | Enable neutral preparation and comparison. The observed gene analysis runs with either setting. Requires `calc_gen_diversity_stats: TRUE`. |
| `markers_file` | | File path to an XLSX/CSV/TSV marker table, ready CHROM/POS list or single-base BED, containing neutral marker annotation and positions. Required in `bed` mode when simulations are enabled and `neutral_vcf` is absent. In `snpeff` mode, leave null/unset for automatic intergenic extraction. A supplied file overrides automatic extraction. |
| `neutral_vcf` | `null` | Path to an already neutral marker filtered vcf-file. If defined, `markers_file` is not considered. |
| `input_format` | `auto` | `auto`, `table`, `positions` or `bed`. `auto` uses BED for `.bed`/`.bed.gz`, otherwise a table. Only tables run the table converter; direct position/BED inputs bypass it. Ignored for automatic SnpEff extraction and supplied `neutral_vcf`. |
| `technical_replicates` | `100` | Number of random neutral sets per gene. Each set has the target gene's final SNP count. A gene with 40 SNPs produces 100 VCFs of 40 SNPs each. |
| `seed` | `42` | Base seed for reproducible marker sampling. This is separate from `wingen_params.SEED`, which controls individual rarefaction. |
| `min_valid_fraction` | `0.9` | Minimum fraction of requested sets yielding finite correlations, checked separately for Pearson and Spearman for each predictor. With 100 sets, at least 90 are required; otherwise that empirical comparison is unavailable. Statistical failures and job errors differ; see below. |
| `extreme_value_percentiles` | `[95, 99]` | Cutoff percentiles in the neutral model correlation distributions (Pearson and Spearman). The values must be between 0 and 100. Examples: `[95]`, `[90, 97.5]`. |

### Table parameters

These parameters apply to XLSX/CSV/TSV tables; they are ignored for ready `positions` or `bed` inputs, automatic SnpEff extraction, and when `neutral_vcf` is supplied.

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


## Neutral source, simulation pool and optional catalog reference

`run_neutral_sim` controls matched Monte Carlo replicates. Independently,
`run_neutral_reference` (true/false only) adds one pooled catalog reference (default
false; supplied config true) named by `neutral_reference_id` (supplied
`neutral_marker_reference`, omitted fallback `All_neutral_markers`). It is an extra
`neutral_pool` category, not a genomic interval or a matched replicate. Either
switch requires a valid neutral source; both false permit no neutral inputs.

`NEUTRAL_SOURCE_POOL` is the complete validated source. `neutral_pool_max_markers`
(omitted/null = all; supplied null) optionally makes a separate DEBUG simulation
pool by seeded uniform sampling without replacement. This changes the sampling
population and may change neutral p-values. `neutral_reference_max_markers`
(omitted/null = the whole effective simulation pool; supplied 100000) independently
caps only the extra reference (and is ignored while `run_neutral_reference` is false).


| Resource setting | Behavior |
|---|---|
| `neutral_preparation_threads: auto` | Snakemake's declared `workflow.cores` budget; explicit positive integers are supported. |
| `neutral_preparation_mem_mb: auto` | Smaller declared global `mem_mb` / `SLURM_MEM_PER_NODE` budget minus `neutral_preparation_headroom_mb` (4096); fallback 8000 MB. |
| `neutral_subset_threads`, `neutral_subset_mem_mb` | Supplied 2 threads / 8000 MB per gene; replicate extraction within a gene remains serial. |
| `neutral_reference_mem_mb: auto` | Same allocation calculation with `neutral_reference_headroom_mb` (4096); fallback 40000 MB. |
| `neutral_reference_threads`, `neutral_reference_runtime_minutes` | Supplied one R worker and 720 minutes for the extra reference. |


## Efficient execution, preview and supplied test settings

`GRID_MIN_CELLS` compares the total rectangular aggregated-grid cell count using
`>=`: zero selects the efficient adapters, `"Inf"` native execution.  `WINGEN_DIAGNOSTICS` enables profiling; expensive native comparisons also
require `GRID_VALIDATE`. `PREVIEW_WINDOWS` defaults true independently. By default the efficient execution is used, diagnostics/validation false and preview true.
These execution controls preserve the configured grid, estimator, rarefaction and
kriging inputs. The final geographic distance mask uses nearest unique empirical
sites, with nearest-neighbor reprojection to the output grid (Mollweide preserves
area, not distance).

## Optional analyses over the finished catalog
Two analyses can be added after the final aggregation. Each reads the finished catalog, stores its result as tables in `gene_catalog.duckdb` and as a TSV next to it, and the catalog page shows it when the tables exist. Both are off by default; with both off, neither job is part of the plan and the DuckDB file and the page hold nothing of them.

| Setting | Behavior |
|---|---|
| `run_correlation_similarity_screen: false` | `true`/`false` only. Adds the gene linkage screening below. On the command line Snakemake reads only `True`/`False`: `--config run_correlation_similarity_screen=True`. |
| `correlation_similarity_min_shared_variables: 4` | A pair of genes is scored only if both have a location-level `r` for at least this many of the environmental variables (latitude included). A whole number of 3 or more: over two variables a correlation is always +1 or -1. |
| `correlation_similarity_top_k_per_gene: 20` | Partners stored per gene, most similar first. The page lists only the partners at or above the 0.90 catalog percentile, so this is a cap. |
| `correlation_similarity_max_full_export_genes: 2000` | The table of every scored pair is also written when the catalog has at most this many genes (2000 genes are about 2 million pairs); `0` never writes it. |
| `run_env_correlation_matrix: false` | `true`/`false` only. Adds the environmental correlation matrix below. |

### Gene linkage screening
Which genes have a similar pattern of location-level correlations with the environment? Each gene has a vector of its location-level Pearson `r` with every environmental variable that the catalog holds (latitude, salinity, summer and winter water temperature, log chlorophyll-a, oxygen, current speed). The similarity of two genes is the Pearson correlation of their two vectors over the variables both have (`n_shared_variables`, at least `correlation_similarity_min_shared_variables`); the Spearman correlation of their `rho` vectors is stored beside it for the partners that are kept, and the percentile says where the pair's |similarity| lies among all pairs of the catalog. The result lists candidate pairs for a later linkage study. It is not a measurement of linkage disequilibrium.

It needs `WRITE_LOCATION_STATS: TRUE` (the location-level `r` and `rho` come from the gene summaries), at least three variables with such values, and at least three genes; without them the job stops with a message that names the missing piece, and the result tables and TSVs are not written. The pooled neutral reference row (`run_neutral_reference`) is a reference sample and not a gene: it stays in the catalog and on the page, and is not screened.

* **What is stored.** In `gene_catalog.duckdb`: `gene_correlation_similarity_pairs` (`gene_id`, `partner`, `similarity_r`, `spearman_similarity_r`, `percentile_vs_all_pairs`, `n_shared_variables`; up to `correlation_similarity_top_k_per_gene` rows per gene) and `gene_correlation_similarity_summary` (`gene_id`, `best_match_gene`, `best_match_r`, `best_match_percentile`: the best partner of each gene), plus the working tables `top_matches` and `cutoffs`. As files in `{results_dir}/results`: `gene_correlation_similarity_pairs.tsv`, `gene_correlation_similarity_summary.tsv` and, for catalogs of at most `correlation_similarity_max_full_export_genes` genes, `gene_correlation_similarity_all_pairs.tsv` (every scored pair, to see what a high similarity looks like in this catalog).
* **On the page.** A "Gene Linkage Screening" group of the column layout (best partner, its similarity and percentile), a "Similar genes" column (the partners at or above the 0.90 percentile, `DISPLAY_PERCENTILE_FLOOR` in `code/build_gene_catalog_html.py`, at most 20, `MAX_SIMILAR_PER_GENE_DISPLAY`, most similar first; a click opens the Pearson and Spearman values, the percentile, the number of shared variables and the breakdown by environmental variable), the chip "Color linkage screening" and the filter "Only genes with a linkage match".  Lowering `MAX_SIMILAR_PER_GENE_DISPLAY` or raising `DISPLAY_PERCENTILE_FLOOR` in `code/build_gene_catalog_html.py` shortens the lists and the file.
* **How it is computed.** Genes are compared with all other genes one block at a time (`--block-genes` of `code/compute_gene_correlation_similarity.py`, 256 genes by default), and only each gene's best partners and a histogram of |similarity| are kept, so the memory depends on the block size and not on the number of pairs. The percentile is read from seven cutoffs of that histogram (0, 50, 90, 95, 99, 99.9 and 100 % of |similarity|; exact to 1e-5) and the pair is placed on that curve: an approximation by design, good for "is this pair in the upper tail", not for the second decimal. The result is the same on every run, with any thread count and any block size. A gene whose values are all equal has no correlation with anything and is not paired.
* **Reading the result.** The percentile ranks a pair among the gene pairs of this catalog.
* **Order and re-runs.** The aggregation writes `gene_catalog.duckdb` first; then the screen, then the matrix, then the page, one after the other because DuckDB lets only one process write to a file. The two jobs follow the catalog TSV, so they run again whenever the aggregation does; a changed setting of the screen re-runs the screen, the matrix (which starts after it) and the page, and nothing else. Switching the screen on for a catalog that is already finished runs only the new job and the page; no gene is touched. Switching it off again leaves its tables in the DuckDB file, and on the page, until the aggregation next runs (for example `--forcerun aggregate_gene_summaries`), because the aggregation writes a new file.
* **Cores.** The screen asks for all of `--cores`; DuckDB uses what it can of them inside a block (on the 2-core test machine two threads were about 15 % faster than one; more cores have not been measured). With an executor that submits one Slurm job per rule, `--cores`/`--jobs` is not the size of a node.

### Environmental correlation matrix
How correlated are the environmental variables with each other across the sample locations? This is a redundancy check on the predictors (the environmental metrics) and involves no gene. It reads the `{gene}_per_location.tsv` files of the catalog (so it needs `WRITE_PER_LOCATION_TSV: TRUE`), takes one value per variable for each of the 10 km location clusters, and correlates every pair of variables over the clusters. Every numeric column of those files except the bookkeeping columns (gene, cluster, counts, longitude, the statistic itself and the `n_individuals_filled_*` flags) is a variable, so a variable added to the Rmd is picked up without a change here. It adds the tables `env_correlation_matrix` (`variable_a`, `variable_b`, `r`, `n_locations`, one row per unordered pair) and `env_correlation_matrix_labels` to `gene_catalog.duckdb`, writes `{results_dir}/results/env_correlation_matrix.tsv`, and the page shows the matrix in a section of its own.

<!-- BEGIN OLS_MODEL_AGREEMENT_DOCUMENTATION -->
## Per-gene environmental model gallery

The main Rmd fits the straight OLS and degree-two quadratic models in R and saves their omitted-cluster predictions, training-mean baseline, curve points and model definitions alongside observed cluster pairs, labelled map PNGs and marker coordinates.
The opened gene/metric page displays those saved values; changing categories or selecting the quadratic comparison does not fit another model.
Existing pipeline genetic calculations and ordinary correlation statistics are unchanged. No server, R installation or internet access is required while viewing results. HTML-only refreshes require these saved R outputs; older raw browser-fit exports must first be regenerated by the R report.

| Setting under `wingen_params` | Supplied value | Purpose |
|---|---|---|
| `LOCATION_SCREENING_ENABLED` | `TRUE` | Enable saved R model, cluster-pair and map exports when genetic statistics and marine environmental data are enabled. |
| `LINKED_PERCENTILE_ENABLED` | `TRUE` | Export the interactive model-agreement page. The master switch above must also be enabled. |
| `LOCATION_SCREENING_ENV_WEIGHTING` | `individual` | Retain inherited means over cluster member individuals. This integration accepts only `individual`. |
| `LINKED_PERCENTILE_GD_PROBS` | `"0.3333333333333333,0.6666666666666667"` | Observed-value thirds used as initial genetic category boundaries. |
| `LINKED_PERCENTILE_MATRIX_ENV_BREAKS` | `"0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"` | Initial environmental deciles; strata can be changed in the browser. |

All paired clusters are shown initially. Disabling either new export switch omits the added pages; inherited catalog plots/statistics remain.
The workflow tracks the main Rmd, required efficiency helper, optional diagnostics when present, and the prediction-agreement HTML template when enabled.
The catalog embeds pinned JavaScript, CSS and licenses for offline use.

The environmental backgrounds use the derived full environmental raster field (not a surface interpolated from sampling-cluster values). Identical render inputs
share a PNG under delivered `results/environment_layers/`. The cache key includes
raster values/grid/CRS, display extent, palette/style, context and render settings;
files are checked and published under a lock. Each page retains its own correctly
transformed marker positions. Genetic backgrounds remain gene-specific.

Neutral simulations are independent of these model fits. With both `run_neutral_sim: false` and `run_neutral_reference: false`, no neutral-marker input is required.
Ordinary environmental correlations and model pages remain available. If either neutral option is enabled, the existing workflow still requires a supplied neutral VCF/markers or configured automatic SnpEff extraction.
Missing requested neutral inputs produce the existing configuration error.

The metric dropdown groups the smaller available Pearson or Spearman neutral
p-value into `Neutral p ≤ 0.01`, `Neutral 0.01 < p ≤ 0.05`, `Neutral p > 0.05`, or
`Neutral p unavailable`. An unavailable neutral comparison never removes an
otherwise usable environmental analysis. Saved failed-repeat statuses are shown
as unavailable; their raw values remain in the tooltip/audit metadata. The numeric
thresholds are applied separately to each gene/metric, without adjusting for the
number of comparisons across the catalog.

<!-- END OLS_MODEL_AGREEMENT_DOCUMENTATION -->
