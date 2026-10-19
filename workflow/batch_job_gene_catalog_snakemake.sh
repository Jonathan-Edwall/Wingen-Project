#!/bin/bash -l
#SBATCH -A uppmax2026-1-187 # Specifying the project ID for compute allocation
#SBATCH --cpus-per-task=8 # 2 genes concurrently x threads=2 each (Snakefile's run_wingen_report) - matches max_parallel_jobs_gene_catalog_creation=2 in the non-Snakemake pipeline
#SBATCH --mem=160G # 2 genes concurrently x ~28.3GB observed peak each (sacct, job 6675358) - same reasoning as the 80G batch_job_gene_catalog.sh now requests
#SBATCH -t 10:00:00
#SBATCH -J gene_catalog_snakemake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edwalljonathan@gmail.com

# module load snakemake/9.21.0 # didn't work
# module load snakemake 
# Loading pixi environment including snakemake and apptainer
eval "$(pixi shell-hook --manifest-path /proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/pixi.toml)"


set -euo pipefail   # exit on error / unset var / failed pipe stage
trap 'echo "ERROR: command [$BASH_COMMAND] failed at line $LINENO"' ERR
echo "Start time: $(date)"
cat "$0"
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
export ROOT="/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project"
cd "$ROOT"

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
# bash ./code/snpeff_vcf_detect_unique_genes.sh body txt /proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/
# # bash ./snpeff_vcf_detect_unique_genes.sh body txt "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"
# python3 ./code/build_gene_chrom_map.py /proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/gene_chrom.tsv" --check-against "/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project/data/raw/genetic_variation_herring/unique_genes.txt"

# CONTAINER="wingen-gene-catalog_v1.sif"
# CONTAINER="wingen-gene-catalog_v2.sif"
# CONTAINER="wingen-gene-catalog_v3.sif"
CONTAINER="wingen-gene-catalog.sif"
# Read input_vcf - and the two gene-list paths build_gene_chrom_map.py needs - straight
# from config/config.yaml instead of hardcoding them a second time here. PyYAML rides
# along with the snakemake install already in this pixi env, so nothing extra to add.
#
# NOTE: config/config.yaml currently points input_vcf at the small head_chr1_1-1e6 TEST
# SUBSET, not the full cohort VCF you ran snpeff_vcf_detect_unique_genes.sh against by
# hand. Update input_vcf in config.yaml to the full file first if you want this step -
# and the snakemake run after it - to build the catalog from the real, full VCF rather
# than the subset.
CONFIG_FILE="config/config.yaml"
get_config_value () {
    python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['$1'])"
}
# Like get_config_value, but for genuinely OPTIONAL keys: prints $2 instead
# of raising a KeyError/traceback when $1 is absent or explicitly null in
# config.yaml. Kept separate from get_config_value above rather than
# changing that one's behavior - the four keys already read through it
# (input_vcf etc.) are all required, and should keep failing loudly if
# they're ever missing, not silently fall back to something plausible-
# looking.
get_config_value_or () {
    python3 -c "
import yaml
val = yaml.safe_load(open('$CONFIG_FILE')).get('$1')
print('$2' if val is None else val)
"
}
INPUT_VCF="$(get_config_value input_vcf)"
# gene_chrom_file is optional: it is only needed for a txt gene list (or the 'genes' list array) combined with vcf_annotation_format snpeff. With a BED gene list the
# Snakefile reads each gene's chromosome from the BED instead (see the chromosome-map step below).
GENE_CHROM_FILE="$(get_config_value_or gene_chrom_file '')"

GENE_LIST_FILE="$(get_config_value_or gene_list_file '')"
VCF_ANNOTATION_FORMAT="$(get_config_value vcf_annotation_format)"
# Determining the source of the gene list file.
# If no file is defined, the 'genes' list ("array") is used as the gene list.
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
RESULTS_DIR="$(get_config_value_or results_dir results)"

echo "input_vcf (from ${CONFIG_FILE}): $INPUT_VCF"
echo "gene_chrom_file (from ${CONFIG_FILE}): $GENE_CHROM_FILE"
echo "gene list type (inferred from gene_list_file / genes): $GENE_LIST_SOURCE"
echo "gene_list_file (from ${CONFIG_FILE}): $GENE_LIST_FILE"
echo "vcf_annotation_format (from ${CONFIG_FILE}): $VCF_ANNOTATION_FORMAT"
echo "results_dir (from ${CONFIG_FILE}): $RESULTS_DIR"
# Both scripts below also fall back to their own hardcoded default if these are left
# out - but the actual pipeline run should always match whatever config.yaml says, not
# silently drift from it, so INPUT_VCF is always passed explicitly here instead of
# relying on that fallback. snpeff_vcf_detect_unique_genes.sh's own 3rd arg (output) is
# passed as "" so it keeps writing next to input_vcf by default, same as before.
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

if [[ "$GENE_LIST_FILE" != */* ]]; then
    GENE_LIST_FILE="${vcf_file_dir_gene_catalog}/${GENE_LIST_FILE}"
fi


# ---------------------------------------------------------
# Creating the gene file if gene_list_file is set but the file does not exist yet
# ---------------------------------------------------------
# No gene_list_file means the inline 'genes' list is used, so there is nothing to create.
# Otherwise the file is written exactly where config.yaml says the Snakefile will read it
# from, in the format its extension asks for (txt or bed - needs a SnpEff-annotated VCF).
if [[ "$GENE_LIST_SOURCE" == "array" ]]; then
    echo "No gene_list_file set; using the 'genes' list from ${CONFIG_FILE}."
elif [[ ! -f "$GENE_LIST_FILE" ]]; then
    echo "GENE_LIST_FILE does not exist: $GENE_LIST_FILE \nCreating gene file ($GENE_LIST_SOURCE)..."
    echo "Creating gene file ($GENE_LIST_SOURCE)..."
    apptainer exec --bind /proj/uppmax2025-2-125 "$CONTAINER" \
        bash ./code/snpeff_vcf_detect_unique_genes.sh body "${GENE_LIST_SOURCE}" "$INPUT_VCF" "$GENE_LIST_FILE"
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

# apptainer exec "$CONTAINER" bash ./code/snpeff_vcf_detect_unique_genes.sh body txt data/raw/genetic_variation_herring/
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

snakemake --unlock
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
