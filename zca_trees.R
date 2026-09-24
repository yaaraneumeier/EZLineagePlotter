# ============================================================================
# zca_trees.R - one figure per tree (10 individuals) with the ZCA cluster row,
# from each individual's latest saved config, rendered with masterlist v84.
#
#   Rscript zca_trees.R            # all 10
#   Rscript zca_trees.R MM-412     # one
#
# BRCA-795 / MM-127: their 07-2026 configs already start with a "Cluster" (ZCA)
# heatmap; used as-is (+ patches.yaml). The 04-2026 configs have none, so a
# "Cluster" heatmap is inserted on top (Set1 over that tree's cluster names).
#
# Output: out/ZCA_trees/Tree___<individual>_ZCA_v84.pdf
#         out/ZCA_trees/<individual>.rds   (rendered tip order + ggplot, for zca_check.R)
# ============================================================================
source("individuals.R")

LATEST <- list(
  "BRCA-775"   = "Trees 04-2026/BRCA-775_v72_config.yaml",
  "BRCA-795"   = "Trees 07-2026/BRCA-795_v82_v10 with CNV_with node81 rotation_colorAdjust_config.yaml",
  "BRCA-841"   = "Trees 04-2026/BRCA-841_v72_config.yaml",
  "MM-127"     = "Trees 07-2026/MM-127_v82_v7 with CNV_config.yaml",
  "MM-412"     = "Trees 04-2026/MM-412_v72_config (2).yaml",
  "MM-423"     = "Trees 04-2026/MM-423_v72_config (1).yaml",
  "NSCLC-0267" = "Trees 04-2026/NSCLC-0267_v72_config.yaml",
  "NSCLC-0401" = "Trees 04-2026/NSCLC-0401_v72_config.yaml",
  "NSCLC-0545" = "Trees 04-2026/NSCLC-0545_v72_V2_config.yaml",
  "NSCLC-2680" = "Trees 04-2026/NSCLC-2680_v72_config.yaml"
)

# Cluster heatmap to insert: the 07-2026 MM-127 entry, restyled to the host config.
cluster_template <- yaml::read_yaml(file.path(KLEIN, LATEST[["MM-127"]]))[["visual definitions"]]$heatmaps[[1]]
stopifnot(identical(unlist(cluster_template$columns), "ZCA"))

cluster_patch <- function(cfg_path) {
  hms <- yaml::read_yaml(cfg_path)[["visual definitions"]]$heatmaps
  if (any(vapply(hms, function(h) "ZCA" %in% unlist(h$columns), logical(1)))) return(NULL)
  h1 <- hms[[1]]
  hm <- cluster_template
  hm$title <- "Cluster"
  hm$distance <- h1$distance        # new row takes the tree-to-heatmap gap ...
  hm$height <- h1$height            # ... and the host's row height
  hm$row_label_font_size <- h1$row_label_font_size
  hm$colnames_angle <- h1$colnames_angle
  hm$row_label_source <- "mapping"
  hm$label_mapping <- list(ZCA = "Cluster")
  hm$custom_discrete <- "no"        # Set1 over this tree's sorted cluster names
  hm$custom_colors <- list()
  list(heatmaps = list(`1` = list(distance = 0)),   # old first row now sits under Cluster
       insert_heatmaps = list(list(at = 1, heatmap = hm)))
}

patches <- yaml::read_yaml("patches.yaml")
args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) names(LATEST)[names(LATEST) %in% args] else names(LATEST)
dir.create("out/ZCA_trees", showWarnings = FALSE, recursive = TRUE)

ez_load_app()
for (ind in todo) {
  cfg <- file.path(KLEIN, LATEST[[ind]])
  spec <- INDIVIDUALS[[ind]]
  y <- yaml::read_yaml(cfg)
  uses_rdata <- any(vapply(y[["visual definitions"]]$heatmaps, function(h) identical(h$data_source, "rdata"), logical(1)))
  patch <- patches[[basename(cfg)]]
  if (is.null(patch)) patch <- cluster_patch(cfg)
  patch <- mutation_label_patch(cfg, patch)          # "ESR1 p.S463P" row labels
  if (ind == "BRCA-795") {
    # The 07-2026 BRCA-795 configs colour Cluster 1B and Cluster 2 the same blue
    # (#377EB8), so 4 Cluster-2 leaves look like 1B. Give Cluster 2 the unused
    # Set1 orange; the other cluster colours are kept as configured.
    cc <- y[["visual definitions"]]$heatmaps[[1]]$custom_colors
    cc[["Cluster 2"]] <- "#FF7F00"
    patch$heatmaps[["1"]]$custom_colors <- cc
  }
  out_file <- file.path("out/ZCA_trees", sprintf("Tree___%s_ZCA_v84.pdf", ind))
  message("\n[zca_trees] ", ind, " <- ", LATEST[[ind]])
  t0 <- Sys.time()
  res <- ez_render(cfg, file.path(TREES, spec$tree), CSV, out_file,
                   rdata = if (uses_rdata) file.path(RDATA, spec$rdata) else NULL,
                   annot = if (uses_rdata) file.path(RDATA, spec$annot) else NULL,
                   patch = patch)
  saveRDS(list(individual = ind, config = LATEST[[ind]], tip_order = res$tip_order, plot = res$plot,
               unmatched_classes = res$unmatched_classes),
          file.path("out/ZCA_trees", paste0(ind, ".rds")))
  message(sprintf("[zca_trees] %s ok=%s %.0fs", ind, isTRUE(res$ok), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
