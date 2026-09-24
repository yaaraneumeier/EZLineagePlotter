# ============================================================================
# render_paper.R - app render of one tree with the paper's settings, for
# paper_fig.R to lay out.
#
#   Rscript paper/render_paper.R BRCA-775      (run from ez_headless/, ~2-3 min)
#
# Heatmap rows: one template for all patients - the BRCA-795 07-2026 config +
# patches.yaml (row labels, CNV chromosome settings) - with PAPER_ROWS edits
# (rows dropped, one palette per row). Per patient only the mutation columns
# change (taken from that patient's own config); rows whose column is empty
# for the patient, and the CNV row when there is no RData, are dropped.
# Tree settings (rotation, stretch, bootstrap, ...) come from the patient's own
# latest config. Output: out/paper/cache/<ind>_plot.rds (the app's ggplot).
# ============================================================================
source("individuals.R")

TEMPLATE <- "Trees 07-2026/BRCA-795_v82_v10 with CNV_with node81 rotation_colorAdjust_config.yaml"
CONFIGS <- list(            # latest config per patient (tree settings + mutation columns)
  "BRCA-795" = TEMPLATE,
  "BRCA-775" = "Trees 04-2026/BRCA-775_v72_config.yaml",
  "BRCA-841" = "Trees 04-2026/BRCA-841_v72_config.yaml",
  "MM-127"   = "Trees 07-2026/MM-127_v82_v7 with CNV_config.yaml",
  "MM-412"   = "Trees 04-2026/MM-412_v72_config (2).yaml",
  "MM-423"   = "Trees 04-2026/MM-423_v72_config (1).yaml",
  "NSCLC-0267" = "Trees 04-2026/NSCLC-0267_v72_config.yaml",
  "NSCLC-0401" = "Trees 04-2026/NSCLC-0401_v72_config.yaml",
  "NSCLC-0545" = "Trees 04-2026/NSCLC-0545_v72_V2_config.yaml",
  "NSCLC-2680" = "Trees 04-2026/NSCLC-2680_v72_config.yaml"
)
# patients sampled at a single time point: no Sampling date row (NSCLC-0267 has
# one 2009_07 cell among 63 from 2015_11 in masterlist v84 - dropped all the same)
# per-patient tree edits (app "manual rotation": node -> children left to right;
# ape node numbers of the trees2025_ordering newick, as the app numbers them)
TREE_EDITS <- list(
  # Cluster 1 in one block: Cluster 2+3 clade (node 93) to the left edge
  "NSCLC-2680" = list(manual_rotation = list(display = "yes", nodes = list(
    list(node = 57, order = c(92, 108, 58)),
    list(node = 92, order = c(93, 104)))))
)
SINGLE_TIMEPOINT <- c("NSCLC-0267", "NSCLC-0401", "NSCLC-0545", "NSCLC-2680")

PAPER_ROWS <- list(
  drop_columns = c("Cell.type.2"),   # Cell type row: redundant with the tree's branch colours
  # one palette per annotation row (column -> value -> colour)
  colors = list(
    Cell.type.1   = c("tumor" = "#8C510A", "non-tumor" = "#DFC27D"),            # browns
    Tissue.type.1 = c("Blood" = "#D55E00", "BM" = "#E69F00", "LN" = "#56B4E9",   # Okabe-Ito
                      "Met" = "#CC79A7", "Oral mucosa" = "#009E73", "PT" = "#F0E442"),
    Cell.type.3   = c("CD3+" = "#1B9E77", "CD68+" = "#7570B3", "CK+" = "#E7298A", # Dark2
                      "OEC" = "#E6AB02", "CD31+" = "#D95F02", "Cellline" = "#66A61E",
                      "EpCAM+" = "#A6761D", "gp100+" = "#666666",
                      "CD19+" = "#A6CEE3", "gp100-MCSP-" = "#BDBDBD"),
    Tissue.type.2 = c("SLN" = "#6A51A3", "NSLN" = "#CBC9E2"),                     # Purples
    NRAS_pQ61K_new3 = c("M" = "#B2182B", "WT" = "#FDDBC7"),                       # RdBu reds
    .mutation     = c("M" = "#1A1A1A", "WT" = "#D9D9D9"),                        # all mutation rows
    WGD           = c("1" = "#8E0152", "0" = "#F1B6DA")                           # PiYG pinks
  ),
  # rows not in the BRCA-795 template: cloned from its Tissue source row and
  # inserted after the row holding `after` (".mutation" = the mutation block)
  extra = list(
    list(column = "Tissue.type.2", title = "LN type", after = "Tissue.type.1"),
    list(column = "NRAS_pQ61K_new3", title = "NRAS p.Q61R survey", after = ".mutation", only = "MM-127")
  ),
  # per-patient palettes: sampling dates light -> dark in time order, evenly
  # spaced by rank (spacing by time crushed 2016-2018 into one dark blue next to
  # a 2009 sample - legend says the colours are ordered); clusters only need
  # distinct colours (paper_fig.R draws the row uncoloured)
  ramp = list(Date = c("#C7E9B4", "#7FCDBB", "#41B6C4", "#225EA8", "#081D58"))
)
# branch colours: the app's lab PAPER_PALETTE, with paper-only edits (keys are
# func.norm.cat() names). Cell line yellow is invisible as a thin branch on white;
# PT-derived epithelial cell has no palette entry (would get a rainbow default).
PALETTE_EDITS <- c("cellline" = "#B59B00", "ptderivedepithelialcell" = "#8FD19E")
# colour for each value of a ramp row, evenly spaced in sort order
ramp_colours <- function(ramp, vals) setNames(colorRampPalette(ramp)(length(vals)), sort(vals))

args <- commandArgs(trailingOnly = TRUE)
IND <- if (length(args)) args[1] else "BRCA-795"
inp <- INDIVIDUALS[[IND]]
own_cfg <- file.path(KLEIN, CONFIGS[[IND]])
dir.create("out/paper/cache", showWarnings = FALSE, recursive = TRUE)
is_mut <- function(cols) length(cols) && all(grepl("^[A-Za-z0-9.]+___[0-9XY]+___[0-9]+___", cols))

# 1. template rows: BRCA-795 config + its patches.yaml entry
tpl_path <- file.path(KLEIN, TEMPLATE)
tpl <- yaml::read_yaml(ez_patch_yaml(tpl_path, yaml::read_yaml("patches.yaml")[[basename(tpl_path)]],
                                     tempfile(fileext = ".yaml")))
hms <- tpl[["visual definitions"]]$heatmaps

row_key <- function(h) { cols <- unlist(h$columns); if (is_mut(cols)) ".mutation" else if (length(cols)) cols[1] else "" }
base_row <- hms[[which(vapply(hms, row_key, "") == "Tissue.type.1")]]
for (ex in PAPER_ROWS$extra) {
  if (!is.null(ex$only) && !IND %in% ex$only) next
  h <- base_row; h$title <- ex$title; h$columns <- list(ex$column)
  h$label_mapping <- setNames(list(ex$title), ex$column)
  at <- which(vapply(hms, row_key, "") == ex$after)
  hms <- append(hms, list(h), after = at)
}

# 2. this patient's mutation columns
own <- yaml::read_yaml(own_cfg)
own_mut <- Filter(function(h) is_mut(unlist(h$columns)), own[["visual definitions"]]$heatmaps)
mut_cols <- unique(unlist(lapply(own_mut, function(h) unlist(h$columns))))
# ... that have calls for this patient (masterlist values other than #N/A / empty)
csv0 <- read.csv(CSV, check.names = TRUE, fileEncoding = "latin1")
csv0 <- csv0[as.character(csv0$new_chosen_sr) %in% ape::read.tree(file.path(TREES, inp$tree))$tip.label, ]
no_calls <- mut_cols[!vapply(mut_cols, function(cn) any(!as.character(csv0[[cn]]) %in% c("", "#N/A", "N/A", "NA", NA)), TRUE)]
if (length(no_calls)) message("[render_paper] ", IND, ": no calls in masterlist for ", paste(no_calls, collapse = ", "))
mut_cols <- setdiff(mut_cols, no_calls)
message("[render_paper] ", IND, " mutation columns: ", paste(mut_cols, collapse = ", "))

# 3. per-patient data: values present for each single-column row
csv <- read.csv(CSV, check.names = TRUE, fileEncoding = "latin1")
tips <- ape::read.tree(file.path(TREES, inp$tree))$tip.label
csv <- csv[as.character(csv$new_chosen_sr) %in% tips, ]
na_like <- c("", "#N/A", "N/A", "NA")
present <- function(col) { v <- as.character(csv[[col]]); sort(unique(v[!is.na(v) & !v %in% na_like])) }

keep <- rep(TRUE, length(hms))
for (i in seq_along(hms)) {
  h <- hms[[i]]; cols <- unlist(h$columns)
  if (identical(h$data_source, "rdata")) { keep[i] <- !is.null(inp$rdata); next }
  if (is_mut(cols)) {
    if (!length(mut_cols)) { keep[i] <- FALSE; next }
    hms[[i]]$columns <- as.list(mut_cols)
    hms[[i]]$label_mapping <- list()          # re-filled by mutation_label_patch()
    hms[[i]]$custom_discrete <- "yes"; hms[[i]]$custom_colors <- as.list(PAPER_ROWS$colors$.mutation)
    next
  }
  if (any(cols %in% PAPER_ROWS$drop_columns)) { keep[i] <- FALSE; next }
  if (IND %in% SINGLE_TIMEPOINT && cols[1] == "Date") { keep[i] <- FALSE; next }
  vals <- present(cols[1])
  if (!length(vals)) { keep[i] <- FALSE; next }
  pal <- if (cols[1] %in% names(PAPER_ROWS$colors)) PAPER_ROWS$colors[[cols[1]]]
         else if (cols[1] %in% names(PAPER_ROWS$ramp)) ramp_colours(PAPER_ROWS$ramp[[cols[1]]], vals)
         else if (cols[1] == "ZCA") setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(max(9, length(vals)))[seq_along(vals)], vals)
         else NULL
  if (!is.null(pal)) {
    miss <- setdiff(vals, names(pal))
    if (length(miss)) stop(IND, ": no colour for ", cols[1], " = ", paste(miss, collapse = ", "))
    hms[[i]]$custom_discrete <- "yes"; hms[[i]]$custom_colors <- as.list(pal[vals])
  }
}
message("[render_paper] rows dropped: ", paste(vapply(hms[!keep], `[[`, "", "title"), collapse = ", "))
hms <- hms[keep]
if (length(hms)) hms[[1]]$distance <- tpl[["visual definitions"]]$heatmaps[[1]]$distance   # tree gap

# 4. patient config with the paper rows, then mutation labels
y <- own
y[["visual definitions"]]$heatmaps <- hms
y[["visual definitions"]]$legend <- tpl[["visual definitions"]]$legend
for (k in names(TREE_EDITS[[IND]])) y[["visual definitions"]][[k]] <- TREE_EDITS[[IND]][[k]]
paper_yaml <- file.path(normalizePath("out/paper/cache"), paste0(IND, "_paper_config.yaml"))
yaml::write_yaml(y, paper_yaml)
ez_patch_yaml(paper_yaml, mutation_label_patch(paper_yaml), paper_yaml)

# 5. render
uses_rdata <- any(vapply(hms, function(h) identical(h$data_source, "rdata"), logical(1)))
ez_load_app()
PAPER_PALETTE_NORM[names(PALETTE_EDITS)] <- PALETTE_EDITS      # read by ez_render from globalenv
res <- ez_render(paper_yaml, file.path(TREES, inp$tree), CSV,
                 file.path("out/paper/cache", paste0(IND, "_app.pdf")),
                 rdata = if (uses_rdata) file.path(RDATA, inp$rdata) else NULL,
                 annot = if (uses_rdata) file.path(RDATA, inp$annot) else NULL)
stopifnot(isTRUE(res$ok))
saveRDS(list(individual = IND, tip_order = res$tip_order, plot = res$plot),
        file.path("out/paper/cache", paste0(IND, "_plot.rds")))
zca <- setNames(as.character(csv$ZCA), as.character(csv$new_chosen_sr))[res$tip_order]
r <- rle(zca)
message("[render_paper] cluster runs left->right: ", paste(sprintf("%s(%d)", sub("Cluster ", "", r$values), r$lengths), collapse = " "))
message("[render_paper] ", IND, " done")
