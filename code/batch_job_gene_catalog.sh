#!/bin/bash -l
#SBATCH -A uppmax2026-1-187 # Specifying the project ID for compute allocation
#SBATCH -c 4 # 2 genes running concurrently now, ~1 core each observed via sacct -> 2 needed + margin
#SBATCH --mem=80G # ~28GB/gene peak observed; additive assumption for 2 concurrent (~56GB) + margin
#SBATCH -t 01:00:00
#SBATCH -J gene_catalog_job # Job_name
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edwalljonathan@gmail.com

####################################  
# Defining the working directory
#################################### 
export ROOT="/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project" # Variable defined in the main pipeline script!
cd $ROOT
# pipeline_script_dir=$ROOT/code/pipeline_scripts
export script_dir="${ROOT}/code" # Variable defined in the main pipeline script!
export pipeline_scripts_dir="${script_dir}" # Variable defined in the main pipeline script!

export results_dir="${ROOT}/results" # Variable defined in the main pipeline script!
export raw_data_dir="${ROOT}/data/raw" # Variable defined in the main pipeline script!




######################################  
####### Defining parameter values #######
######################################
export gene_list_source="txt" # temp 2026-09-14
export gene_list_file="${raw_data_dir}/genetic_variation_herring/unique_genes.txt" # temp 2026-09-14

# export gene_scope="cds" # "all" (default) | "body" | "exonic" | "cds"
export gene_scope="all" # "all" (default) | "body" | "exonic" | "cds"
export gene_list_limit="2" # temp 2026-09-14


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

conda_setup_script_path="/proj/uppmax2025-2-125/nobackup/jonathan/anaconda3/etc/profile.d/conda.sh"  # Variable defined in config.sh
source "$conda_setup_script_path"  # Source Conda initialization script
# Activate the conda environment
conda activate krill_env
export ROOT="/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project" # Variable defined in the main pipeline script!

cd "$ROOT"
export script_dir="${ROOT}/code" # Variable defined in the main pipeline script!
export pipeline_scripts_dir="${script_dir}" # Variable defined in the main pipeline script!
script_name="2_wingen_create_gene_catalog.sh"

echo "Pipeline code version (git commit): $(git -C "$script_dir" rev-parse HEAD 2>/dev/null || echo 'not a git repo, or git unavailable')"
echo "Starting $script_name at: $(date)"
source "$pipeline_scripts_dir/$script_name"
echo "Finished $script_name at: $(date)"

ls -l # Checking if outputs was generated

module list # Useful for keeping track of the versions used
echo "JOB = ${SLURM_JOBID:-N/A}" # N/A when run directly with bash instead of sbatch
echo "End time: $(date)"

# Resource sizing: single-gene baseline was -c 2 --mem=40G (see git history of this
# comment). Now running 2 genes concurrently (max_parallel_jobs_gene_catalog_creation=2 in 2_wingen_create_gene_catalog.sh), scaled accordingly.
# re-check with sacct below once this actually runs, since additive scaling is an assumption, not yet confirmed:
#   TotalCPU/Elapsed should show ~2 (vs ~1 single-gene) if both truly run concurrently
#   MaxRSS reflects the combined peak across both concurrent R processes (same cgroup)

if [ -n "${SLURM_JOBID:-}" ]; then
    echo "Queuing a sacct check to run once this job actually completes (not from inside it - that only gives partial/RUNNING-state data)..."
    sbatch --dependency=afterany:"$SLURM_JOBID" -A uppmax2026-1-187 -c 1 -t 00:05:00 \
           -J "sacct_${SLURM_JOBID}" --chdir="$ROOT" \
           --wrap="sacct -j $SLURM_JOBID --format=JobID,Elapsed,AllocCPUS,TotalCPU,CPUTime,MaxRSS,ReqMem,State"
fi

# if [ -n "${SLURM_JOBID:-}" ]; then
#     echo "Queuing a jobstats check to run once this job completes..."
#     sbatch --dependency=afterany:"$SLURM_JOBID" -A uppmax2026-1-187 -c 1 -t 00:10:00 \
#            -J "jobstats_${SLURM_JOBID}" --chdir="$ROOT" --wrap="bash -l -c 'jobstats -p $SLURM_JOBID'"
# fi