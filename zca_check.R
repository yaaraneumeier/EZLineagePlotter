# ============================================================================
# zca_check.R - check the ZCA "Cluster" row as drawn in the zca_trees.R figures
# against the leaf -> cluster map read off the slides (ZCA_CLUSTER_HANDOFF.md).
#
# Reads what is drawn, not the masterlist: from each rendered ggplot it takes
# the Cluster heatmap's tiles (tip position + fill colour), decodes each colour
# through that layer's fill scale back to a cluster name, and pairs tiles with
# the tip labels at the same position.
#
# Checks per tree:
#   drawn cluster == slide cluster, per leaf
#   one colour per cluster name (legend decodable)
#   left-to-right leaf order in the figure == slide plot_index order
#
# Output: out/ZCA_trees/zca_check_summary.tsv, out/ZCA_trees/zca_check_leaves.tsv
# ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggtree) })
KLEIN <- normalizePath(Sys.getenv("EZ_DATA", file.path(Sys.getenv("EZ_ROOT", "."), "data")), mustWork = TRUE)   # run from the repo root
slides <- read.csv(file.path(KLEIN, "zca_leaf_cluster_map_20260909.csv"), colClasses = "character")
IND2SLIDE <- c("BRCA-775" = "775-13", "BRCA-795" = "795-09", "BRCA-841" = "841-12", "MM-127" = "MM15-127",
               "MM-412" = "MM16-412", "MM-423" = "MM16-423", "NSCLC-0267" = "BC15-0267",
               "NSCLC-0401" = "BC16-0401", "NSCLC-0545" = "BC16-0545", "NSCLC-2680" = "BC14-2680")

drawn_clusters <- function(p) {
  b <- suppressWarnings(ggplot_build(p))
  # the heatmap layer whose source column is ZCA
  k <- which(vapply(p$layers, function(l) is.data.frame(l$data) && "column" %in% names(l$data) &&
                      all(l$data$column == "ZCA"), logical(1)))
  stopifnot(length(k) == 1)
  d <- b$data[[k]]
  fill_col <- grep("^fill", names(d), value = TRUE)[1]
  aes_name <- fill_col
  sc <- b$plot$scales$get_scales(aes_name)
  brk <- sc$get_limits()
  key <- setNames(as.character(brk), toupper(sc$map(brk)))           # colour -> cluster name
  # tip labels at each tile row (tile y and label y share the tip axis)
  lab_layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomTextGGtree"), logical(1)))[1]
  labs <- b$data[[lab_layer]][, c("y", "label")]
  tiles <- data.frame(y = d$y, colour = toupper(d[[fill_col]]))
  tiles$drawn <- unname(key[tiles$colour])
  m <- merge(tiles, labs, by = "y")
  m <- m[order(m$y), ]                                               # left -> right in the figure (tip axis runs -N..-1)
  list(leaves = data.frame(leaf_id = m$label, fig_index = seq_len(nrow(m)) - 1, drawn = m$drawn, colour = m$colour),
       key = key, n_tiles = nrow(d))
}

summ <- list(); leaves_all <- list()
for (ind in names(IND2SLIDE)) {
  f <- file.path("out/ZCA_trees", paste0(ind, ".rds"))
  if (!file.exists(f)) { message("missing ", f); next }
  r <- readRDS(f)
  dc <- drawn_clusters(r$plot)
  s <- slides[slides$individual == IND2SLIDE[[ind]], c("leaf_id", "plot_index", "cluster")]
  s$plot_index <- as.integer(s$plot_index)
  j <- merge(dc$leaves, s, by = "leaf_id", all = TRUE)
  j <- j[order(j$fig_index), ]
  j$match <- !is.na(j$drawn) & !is.na(j$cluster) & j$drawn == j$cluster
  same_order <- isTRUE(all(j$fig_index == j$plot_index))
  first_off <- if (same_order) NA else j$leaf_id[which(j$fig_index != j$plot_index)[1]]
  summ[[ind]] <- data.frame(
    individual = ind, slide_individual = IND2SLIDE[[ind]], config = r$config,
    leaves_fig = nrow(dc$leaves), leaves_slide = nrow(s), tiles = dc$n_tiles,
    cluster_match = sum(j$match), cluster_mismatch = sum(!j$match),
    colours_unique = !anyDuplicated(names(dc$key)),   # one colour per cluster name
    order_equals_slides = same_order,
    order_positions_differ = sum(j$fig_index != j$plot_index, na.rm = TRUE),
    first_leaf_out_of_order = first_off,
    tip_order_equals_app = identical(as.character(dc$leaves$leaf_id), as.character(r$tip_order)),
    legend = paste(sprintf("%s=%s", dc$key, names(dc$key)), collapse = "; "),
    stringsAsFactors = FALSE)
  leaves_all[[ind]] <- cbind(individual = ind, j)
}
S <- do.call(rbind, summ); L <- do.call(rbind, leaves_all)
write.table(S, "out/ZCA_trees/zca_check_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(L, "out/ZCA_trees/zca_check_leaves.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(S[, c("individual", "leaves_fig", "leaves_slide", "cluster_match", "cluster_mismatch",
            "colours_unique", "order_equals_slides", "order_positions_differ", "tip_order_equals_app")], row.names = FALSE)
