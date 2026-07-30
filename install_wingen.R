# install_wingen.R
#
# wingen (https://CRAN.R-project.org/package=wingen) is not distributed as a
# conda/mamba package on any channel - only CRAN and GitHub host it. Every
# dependency IT needs, however, has already been installed via conda from
# environment_krill.yml. So this step just fetches the wingen sources
# themselves; dependencies = FALSE tells install.packages() not to reach
# out to CRAN for anything else, so if something is missing you'll get an
# immediate, clear "there is no package called X" error instead of a silent
# source build that could fight with the conda-installed binaries.
#
# Usage:
#   mamba env create -f environment_krill.yml
#   conda activate krill_env
#   Rscript install_wingen.R

install.packages("wingen", repos = "https://cloud.r-project.org", dependencies = FALSE)

library(wingen)
cat("wingen", as.character(packageVersion("wingen")), "installed OK\n")
