# ============================================================================
# EZLineagePlotter - container image for running the FULL app locally.
#
# Based on the Bioconductor image because ggtree/treeio are Bioconductor
# packages that install cleanly there (no slow, fragile source compile).
#
# Build (normally done by GitHub Actions, not by hand):
#   docker build -t ezlineageplotter .
# Run:
#   docker run --rm -p 3838:3838 ezlineageplotter
# then open http://localhost:3838
# ============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_19

# System libraries some R packages need (fonts for showtext, units for ggforce,
# image libs for jpeg/png). Most are already present in the base image; this
# makes the build self-contained.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libudunits2-dev \
      libfontconfig1-dev \
      libfreetype6-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libpng-dev \
      libjpeg-dev \
      libtiff5-dev \
      libxml2-dev \
      fonts-dejavu-core \
      fonts-dejavu \
      fonts-liberation \
 && rm -rf /var/lib/apt/lists/*

# R packages (CRAN + Bioconductor). BiocManager::install handles both and picks
# versions matching the Bioconductor release. The final check fails the build
# loudly if any package did not install.
RUN R -q -e 'options(Ncpus = parallel::detectCores()); \
  pkgs <- c("shiny","shinydashboard","shinyjs","shinyWidgets","shinyBS","colourpicker","DT", \
            "ape","dplyr","ggplot2","ggforce","yaml","stringr","scales","tidyverse", \
            "ggnewscale","gridExtra","cowplot","jpeg","png","combinat","infotheo","aricode", \
            "showtext","sysfonts","svglite","data.table","RColorBrewer","htmltools","later", \
            "ggtree","treeio"); \
  BiocManager::install(pkgs, update = FALSE, ask = FALSE); \
  miss <- pkgs[!pkgs %in% rownames(installed.packages())]; \
  if (length(miss)) stop("Failed to install: ", paste(miss, collapse = ", "))'

# App files: single-tree app, multi-tree module, combined launcher, entry point.
WORKDIR /app
COPY EZlineagePlotter56.R EZlineagePlotter56_mt.R EZlineagePlotter56_combined.R app.R /app/

EXPOSE 3838

# Serve the full app (single + multiple-trees modes) on all interfaces so the
# host browser can reach it. Single Tree is the default; ?mode=multi for multi.
CMD ["R", "-q", "-e", "options(shiny.devmode = FALSE, shiny.autoreload = FALSE); shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
