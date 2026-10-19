#!/bin/bash -l
#SBATCH -A uppmax2026-1-187
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH -t 00:45:00
#SBATCH -J wingen_environment_debug

set -euo pipefail

# sbatch executes a copy of this file from its spool directory. Submit from
# the project root so SLURM_SUBMIT_DIR identifies the installed code instead.
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/code/debug_environment.R" ]]; then
    ROOT="$(cd -- "$SLURM_SUBMIT_DIR" && pwd)"
else
    ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export ROOT
cd -- "$ROOT"

if [[ ! -f "$ROOT/code/debug_environment.R" ]]; then
    echo "Missing $ROOT/code/debug_environment.R. Copy the diagnostic files into the project and submit from its root." >&2
    exit 2
fi

# Use the project's existing environment and container. Nothing is installed.
if command -v pixi >/dev/null 2>&1 && [[ -f "$ROOT/pixi.toml" ]]; then
    pixi_hook="$(pixi shell-hook --manifest-path "$ROOT/pixi.toml")"
    eval "$pixi_hook"
    unset pixi_hook
fi
if ! command -v apptainer >/dev/null 2>&1; then
    echo "apptainer is unavailable. Activate the same environment used by the pipeline before submission." >&2
    exit 2
fi
container="$ROOT/wingen-gene-catalog.sif"
if [[ ! -f "$container" ]]; then
    echo "Missing existing container: $container" >&2
    exit 2
fi

mkdir -p -- "$ROOT/debug"
output_dir="$(mktemp -d "$ROOT/debug/environment-${SLURM_JOB_ID:-manual}-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")"
log_file="$output_dir/debug.log"
{
    echo "Start: $(date -u +%FT%TZ)"
    echo "Project: $ROOT"
    echo "Job: ${SLURM_JOB_ID:-manual}"
    echo "Diagnostics directory: $output_dir"
    echo "TMPDIR is left unchanged; R will report its actual temporary directory."
} | tee "$log_file"

# Arguments supplied to this wrapper are passed to the R diagnostic. Later
# --config, --rmd or --out arguments override these defaults. The wrapper's
# debug.log remains in the unique directory printed above if --out is changed.
if apptainer exec --home "$ROOT" --bind /proj/uppmax2025-2-125 --pwd "$ROOT" "$container" \
    Rscript "$ROOT/code/debug_environment.R" \
    --config "$ROOT/config/config.yaml" \
    --rmd "$ROOT/code/2_Herring_Wingen.Rmd" \
    --out "$output_dir" "$@" 2>&1 | tee -a "$log_file"; then
    status=0
else
    status=$?
fi
{
    echo "Finish: $(date -u +%FT%TZ)"
    echo "Exit status: $status"
    echo "Diagnostic log: $log_file"
} | tee -a "$log_file"
exit "$status"
