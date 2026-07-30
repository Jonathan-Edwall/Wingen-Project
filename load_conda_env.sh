#!/usr/bin/env bash
# load_conda_env.sh
set -e

source "$(conda info --base)/etc/profile.d/conda.sh"


conda activate base
# This step requires that conda-lock is installed
conda-lock install \
    --name krill_env \
    environment_krill_lock.yml
    
# Activating the environment
conda activate krill_env

# Installing wingen for R
Rscript install_wingen.R
