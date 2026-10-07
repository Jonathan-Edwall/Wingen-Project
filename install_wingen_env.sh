#!/usr/bin/env bash
set -euo pipefail

LOCKFILE="environment_wingen_lock.conda-lock.yml"
ENV_NAME="wingen_env"

conda-lock install --name "$ENV_NAME" "$LOCKFILE"

conda run -n "$ENV_NAME" Rscript -e '
remotes::install_version("wingen", version = "2.2.0", repos = "https://cran.r-project.org")
remotes::install_version("CopernicusMarine", version = "0.4.9", repos = "https://cran.r-project.org")
remotes::install_version("SpatialKDE", version = "0.8.2", repos = "https://cran.r-project.org")
'
