#!/bin/bash -l
#SBATCH -A uppmax2026-1-187 # Specifying the project ID for compute allocation
#SBATCH --cpus-per-task=8 # 2 genes concurrently x threads=2 each (Snakefile's run_wingen_report) - matches max_parallel_jobs_gene_catalog_creation=2 in the non-Snakemake pipeline
#SBATCH --mem=160G # 2 genes concurrently x ~28.3GB observed peak each (sacct, job 6675358) - same reasoning as the 80G batch_job_gene_catalog.sh now requests
#SBATCH -t 10:00:00
#SBATCH -J gene_catalog_snakemake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edwalljonathan@gmail.com

# Submit from the project root with sbatch; direct Bash runs locate it via workflow/.
export ROOT="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}"
cat "$0"
cd "$ROOT"
# module load snakemake/9.21.0 # didn't work
# module load snakemake 
# Loading pixi environment including snakemake and apptainer
eval "$(pixi shell-hook --manifest-path "$ROOT/pixi.toml")"

set -euo pipefail   # exit on error / unset var / failed pipe stage
trap 'echo "ERROR: command [$BASH_COMMAND] failed at line $LINENO"' ERR
echo "Start time: $(date)"
echo "USER = $USER"
echo "QOS = ${SLURM_JOB_QOS:-N/A}"  # N/A when run directly with bash instead of sbatch

if [ -n "${SLURM_JOBID:-}" ]; then
    echo "--- Allocated resources (scontrol) ---"
    scontrol show job "$SLURM_JOBID"
fi

# conda_setup_script_path="/proj/uppmax2025-2-125/nobackup/jonathan/anaconda3/etc/profile.d/conda.sh"
# source "$conda_setup_script_path"  # Source Conda initialization script
# # Activate the conda environment - snakemake itself comes from the module
# # loaded above; wingen_env here is what the actual rule shell commands
# # (bcftools, R) run against directly now, since neither rule uses a
# # container: any more.
# conda activate wingen_env
# # wait

# which python
# python -c "import sys; print(sys.executable)"
# python -c "import duckdb; print(duckdb.__version__)"
# export ROOT="${ROOT:-/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project}"
# cd "$ROOT"

# # PROJ can't find its own database inside the container unless told where it is
# export APPTAINERENV_PROJ_DATA=/opt/conda/share/proj
# export APPTAINERENV_PROJ_LIB=/opt/conda/share/proj

# Note what's gone from here vs. the old batch_job_gene_catalog.sh: no
# script_dir/pipeline_scripts_dir/wingen_output_dir/LAT_MIN/etc exports -
# all of that lives in config/config.yaml now, read once by the Snakefile
# instead of threaded through the environment by hand.

echo "Pipeline code version (git commit): $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo, or git unavailable')"
echo "Starting snakemake at: $(date)"

# --cores 4 matches -c above: local parallelism within this one
# allocation, same overall shape as the old script - Snakemake's DAG
# handles ordering/dependencies instead of the wait -n loop.
# No --snakefile flag needed: workflow/Snakefile is Snakemake's own
# standard project layout, auto-discovered from the repo root the same
# way a bare ./Snakefile would be - this is the payoff of following that
# convention rather than the code/Snakefile placement from before.
# --resources mem_mb=80000 matches --mem above - without this, Snakemake's
# local scheduler only throttles concurrency by cores/threads, not
# memory, and could try to start more concurrent renders than this
# allocation can actually hold; this makes it respect the same ceiling
# SLURM already gave the job.

# No --sdm apptainer any more - all three rules (including
# narrow_vcf_by_scope, easy to miss since it's not always in the DAG)
# dropped their container: directives and now run straight against the
# already-activated wingen_env, same as the non-Snakemake pipeline does.
# The Dockerfile is still in the repo if this ever needs to change back -
# just not part of the active path right now.

# input_vcf: "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
# gene_list_file: "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/unique_genes.txt"  # required if gene_list_source is txt or bed
# bash ./code/0_snpeff_vcf_detect_unique_genes.sh body txt /proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/
# # bash ./code/0_snpeff_vcf_detect_unique_genes.sh body txt "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
# python3 ./code/build_gene_chrom_map.py /proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/gene_chrom.tsv" --check-against "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/unique_genes.txt"

# CONTAINER="wingen-gene-catalog_v1.sif"
# CONTAINER="wingen-gene-catalog_v2.sif"
# CONTAINER="wingen-gene-catalog_v3.sif"
CONTAINER="wingen-gene-catalog.sif"
# Read input_vcf - and the two gene-list paths build_gene_chrom_map.py needs - straight
# from config/config.yaml instead of hardcoding them a second time here. PyYAML rides
# along with the snakemake install already in this pixi env, so nothing extra to add.
# Optional annotation-table settings are read from that same config below.
#
# NOTE: input_vcf in config.yaml selects the dataset for both this preparation step
# and the subsequent Snakemake run. Update it there when switching datasets.
CONFIG_FILE="config/config.yaml"
get_config_value_or_default () {
    python3 -c '
import sys, yaml
value = yaml.safe_load(open(sys.argv[1], encoding="utf-8")).get(sys.argv[2])
if isinstance(value, list):
    value = ",".join(value)
print(sys.argv[3] if value is None else value)
' "$CONFIG_FILE" "$1" "${2-}"
}
INPUT_VCF="$(get_config_value_or_default input_vcf '')"
if [[ -z "$INPUT_VCF" ]]; then
    echo "ERROR: input_vcf must be set in ${CONFIG_FILE}." >&2
    exit 1
fi
# gene_chrom_file is optional: it is only needed for a txt gene list (or the 'genes' list array) combined with vcf_annotation_format snpeff. With a BED gene list the
# Snakefile reads each gene's chromosome from the BED instead (see the chromosome-map step below).
GENE_CHROM_FILE="$(get_config_value_or_default gene_chrom_file '')"

GENE_LIST_FILE="$(get_config_value_or_default gene_list_file '')"
VCF_ANNOTATION_FORMAT="$(get_config_value_or_default vcf_annotation_format '')"
if [[ "$VCF_ANNOTATION_FORMAT" != "snpeff" && "$VCF_ANNOTATION_FORMAT" != "bed" ]]; then
    echo "ERROR: vcf_annotation_format must be set to snpeff or bed in ${CONFIG_FILE}." >&2
    exit 1
fi
# A configured gene file takes precedence over the inline gene list array.
if [[ -z "$GENE_LIST_FILE" ]]; then
    GENE_LIST_SOURCE="array"
elif [[ "${GENE_LIST_FILE,,}" == *.bed ]]; then
    GENE_LIST_SOURCE="bed"
else
    GENE_LIST_SOURCE="txt"
fi
# Root of THIS catalog's whole output tree (filtered VCFs, .done markers,
# benchmarks, per-chrom PNGs/JSONs, the aggregated duckdb/tsv/html) - see
# the comment on results_dir in config.yaml itself. Optional: old configs
# that don't set it keep writing to "results", exactly as before.
RESULTS_DIR="$(get_config_value_or_default results_dir results)"

echo "input_vcf (from ${CONFIG_FILE}): $INPUT_VCF"
echo "gene_chrom_file (from ${CONFIG_FILE}): $GENE_CHROM_FILE"
echo "gene list type (inferred from gene_list_file / genes): $GENE_LIST_SOURCE"
echo "gene_list_file (from ${CONFIG_FILE}): $GENE_LIST_FILE"
echo "vcf_annotation_format (from ${CONFIG_FILE}): $VCF_ANNOTATION_FORMAT"
echo "results_dir (from ${CONFIG_FILE}): $RESULTS_DIR"
# Gene-list filenames without a directory are placed beside input_vcf.
# Paths containing a directory are kept as configured, as in the Snakefile.
#
# --bind /proj/uppmax2025-2-125 added here to match the snakemake --apptainer-args call
# further down; drop it if these two apptainer exec calls already see that path fine
# without it (untested - the failed run being fixed here was run outside the container).
vcf_file_dir_gene_catalog="$(dirname "$INPUT_VCF")"

# If no gene list file is defined and the VCF-file annotation format is set to "snpeff", then an arbitrarily named gene list file will be created
if [[ -z "$GENE_LIST_FILE" && "$VCF_ANNOTATION_FORMAT" == "snpeff" ]]; then
    GENE_LIST_FILE="${vcf_file_dir_gene_catalog}/unique_genes.txt"
    GENE_LIST_SOURCE="txt"
fi

# BED mode uses an existing annotation table to create a missing gene list.
if [[ -z "$GENE_LIST_FILE" && "$VCF_ANNOTATION_FORMAT" == "bed" ]]; then
    GENE_ANNOTATION_FILE="$(get_config_value_or_default gene_annotation_file '')"
    if [[ -n "$GENE_ANNOTATION_FILE" && -f "$GENE_ANNOTATION_FILE" ]]; then
        GENE_LIST_FILE="${vcf_file_dir_gene_catalog}/unique_genes.bed"
        GENE_LIST_SOURCE="bed"
    fi
fi

if [[ -n "$GENE_LIST_FILE" && "$GENE_LIST_FILE" != */* ]]; then
    GENE_LIST_FILE="${vcf_file_dir_gene_catalog}/${GENE_LIST_FILE}"
fi

# ---------------------------------------------------------
# Creating the gene file if gene_list_file is set but the file does not exist yet
# ---------------------------------------------------------
# Without a configured or automatically named file, the inline 'genes' list is used.
# Otherwise the file is written exactly where config.yaml says the Snakefile will read it
# from. The selected vcf_annotation_format determines which converter creates it.
if [[ "$GENE_LIST_SOURCE" == "array" ]]; then
    echo "No gene_list_file set; using the 'genes' list from ${CONFIG_FILE}."
elif [[ -d "$GENE_LIST_FILE" ]]; then
    echo "ERROR: gene_list_file must name a file, not an existing directory: $GENE_LIST_FILE" >&2
    exit 1
elif [[ ! -f "$GENE_LIST_FILE" ]]; then
    echo "GENE_LIST_FILE does not exist: $GENE_LIST_FILE"
    # snpeff VCF-file: Creating a gene list from a SnpEff-processed VCF-file using 0_snpeff_vcf_detect_unique_genes.sh
    if [[ "$VCF_ANNOTATION_FORMAT" == "snpeff" ]]; then
        echo "Creating gene file from SnpEff annotations ($GENE_LIST_SOURCE)..."
        apptainer exec --bind /proj/uppmax2025-2-125 "$CONTAINER" \
            bash ./code/0_snpeff_vcf_detect_unique_genes.sh body "${GENE_LIST_SOURCE}" "$INPUT_VCF" "$GENE_LIST_FILE"
    else
        GENE_ANNOTATION_FILE="$(get_config_value_or_default gene_annotation_file '')"
        if [[ -z "$GENE_ANNOTATION_FILE" || ! -f "$GENE_ANNOTATION_FILE" ]]; then
            echo "ERROR: missing BED gene list; set gene_annotation_file to an existing gene-coordinate table, or supply the BED file yourself." >&2
            exit 1
        fi
        if [[ "$GENE_LIST_SOURCE" != "bed" ]]; then
            echo "ERROR: BED mode requires gene_list_file ending in .bed." >&2
            exit 1
        fi
        if [[ "$(get_config_value_or_default gene_scope all)" != "all" ]]; then
            echo "ERROR: table-derived BED intervals require gene_scope: all; biotype filters select genes, not exons/CDS." >&2
            exit 1
        fi
        # Creating a .bed-format gene list from a supplied annotation table using 0_annotation_table_detect_unique_genes.sh
        # The converter supplies column names, gene bounds and BED formatting.
        annotation_args=(bed "$GENE_ANNOTATION_FILE" "$GENE_LIST_FILE" --header)
        include_biotypes="$(get_config_value_or_default gene_annotation_include_biotypes '')"
        [[ -z "$include_biotypes" ]] || annotation_args+=(--include-biotypes "$include_biotypes")

        echo "Creating BED gene file from annotation table: $GENE_ANNOTATION_FILE"
        apptainer exec --bind /proj/uppmax2025-2-125 "$CONTAINER" \
            bash ./code/0_annotation_table_detect_unique_genes.sh "${annotation_args[@]}"
    fi
else
    echo "Using existing gene list file: ${GENE_LIST_FILE}"
fi

# ---------------------------------------------------------
# Extracting chromosomes if GENE_CHROM_FILE is empty or missing
# ---------------------------------------------------------
if [[ "$GENE_LIST_SOURCE" != "bed" && "$VCF_ANNOTATION_FORMAT" != "bed" ]]; then
    if [[ -z "$GENE_CHROM_FILE" || ! -f "$GENE_CHROM_FILE" ]]; then
        # If gene_chrome_file is not defined, then an arbitrarily named gene_chrome_file will be created
        if [[ -z "$GENE_CHROM_FILE" ]]; then
            GENE_CHROM_FILE="${vcf_file_dir_gene_catalog}/gene_chrom.tsv"
        else
            if [[ "$GENE_CHROM_FILE" != */* ]]; then
                GENE_CHROM_FILE="${vcf_file_dir_gene_catalog}/${GENE_CHROM_FILE}"
            fi
        fi

        echo "GENE_CHROM_FILE is empty or does not exist.\nCreating gene chromosome map..."
        mkdir -p "$(dirname "$GENE_CHROM_FILE")"
        if [[ -n "$GENE_LIST_FILE" ]]; then
            apptainer exec --bind /proj/uppmax2025-2-125 "$CONTAINER" \
                python3 ./code/build_gene_chrom_map.py "$INPUT_VCF" "$GENE_CHROM_FILE" --check-against "$GENE_LIST_FILE"
        else
            apptainer exec --bind /proj/uppmax2025-2-125 "$CONTAINER" \
                python3 ./code/build_gene_chrom_map.py "$INPUT_VCF" "$GENE_CHROM_FILE"
        fi
        # script_end=$(date +%s) # Ending the timer 
        # script_runtime=$((script_end-script_start)) # Calculating the script_runtime of the script
        # echo "./code/build_gene_chrom_map.py Runtime: $script_runtime seconds. \n Outputfile: ${GENE_CHROM_FILE}"
    else
        echo "Using existing gene chromosome file: ${GENE_CHROM_FILE}"
    fi
else
    echo "GENE_LIST_SOURCE is '$GENE_LIST_SOURCE' and vcf_annotation_format is '$VCF_ANNOTATION_FORMAT'; a BED file supplies the chromosomes, so no GENE_CHROM_FILE is required."
fi

# apptainer exec "$CONTAINER" bash ./code/0_snpeff_vcf_detect_unique_genes.sh body txt data/raw/genetic_variation_herring/
# python3 code/build_gene_chrom_map.py \
#     data/raw/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz \
#     data/raw/genetic_variation_herring/gene_chrom.tsv \
#     --check-against data/raw/genetic_variation_herring/unique_genes.txt

# snakemake \
#     --cores 4 \
#     --resources mem_mb=80000 \
#     --rerun-incomplete \
#     --notemp \
#     --sdm conda
# snakemake --unlock

# snakemake --cores 4 --resources mem_mb=80000 --rerun-incomplete --notemp \
#           --sdm conda apptainer \
#           --apptainer-args "--bind /proj/uppmax2025-2-125"

# snakemake --cores 4 --resources mem_mb=80000 --rerun-incomplete --notemp \
#           --sdm apptainer \
#           --apptainer-args "--bind /proj/uppmax2025-2-125"apptainer exec --bind /proj/uppmax2025-2-125 "$CONTAINER" \
snakemake --unlock --cores "${SLURM_CPUS_PER_TASK:-4}"
# snakemake --unlock
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-<unset>}"
echo "SLURM_MEM_PER_NODE  = ${SLURM_MEM_PER_NODE:-<unset>}"
snakemake --cores "${SLURM_CPUS_PER_TASK:-4}" \
          --resources mem_mb="${SLURM_MEM_PER_NODE:-80000}" \
          --rerun-incomplete --notemp \
          --sdm apptainer \
          --apptainer-args "--bind /proj/uppmax2025-2-125" \
          --report "$RESULTS_DIR/results/workflow_timing.html" \
          --report-after-run

# snakemake --cores 2 --resources mem_mb=40000 --rerun-incomplete --notemp \
#           --sdm apptainer \
#           --apptainer-args "--bind /proj/uppmax2025-2-125"

# snakemake --cores 2 --resources mem_mb=40000 --rerun-incomplete --notemp \
#           --sdm apptainer conda \
#           --apptainer-args "--bind /proj/uppmax2025-2-125"

# snakemake --cores 4 --resources mem_mb=80000 --rerun-incomplete --notemp

# snakemake --unlock
# snakemake --cores 2 \
#           --resources mem_mb=40000 \
#           --rerun-incomplete \
#           --notemp \
#           --sdm apptainer

# snakemake --cores 6 --resources mem_mb=120000

echo "Finished snakemake at: $(date)"

# ls -lR "$RESULTS_DIR" # Checking if outputs was generated (.snakemake_gene_catalog/, benchmarks/, temp/, results/)
echo "Results directory: $RESULTS_DIR"
module list # Useful for keeping track of the versions used
echo "JOB = ${SLURM_JOBID:-N/A}" # N/A when run directly with bash instead of sbatch
echo "End time: $(date)"

if [ -n "${SLURM_JOBID:-}" ]; then
    echo "Queuing a sacct check to run once this job actually completes (not from inside it - that only gives partial/RUNNING-state data)..."
    sbatch --dependency=afterany:"$SLURM_JOBID" -A uppmax2026-1-187 -c 1 -t 00:05:00 \
           -J "sacct_${SLURM_JOBID}" --chdir="$ROOT" \
           --wrap="sacct -j $SLURM_JOBID --format=JobID,Elapsed,AllocCPUS,TotalCPU,CPUTime,MaxRSS,ReqMem,State"
fi
