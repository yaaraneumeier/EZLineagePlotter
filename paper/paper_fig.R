# ============================================================================
# paper_fig.R - production layout of one tree figure for the paper.
#
#   Rscript paper/paper_fig.R <IND> <tag> [--pdf]     (run from ez_headless/)
#
# Input: out/paper/cache/<IND>_plot.rds - the app's ggplot, from render_paper.R.
# Re-lays it out without re-rendering (~1 min, mostly loading the app):
#   - tree depth on a pseudo-log scale, thinner branches;
#   - cuts the plot into bands along the depth axis (tree + tip labels,
#     annotation rows, CNV) and stacks them at fixed mm heights;
#   - Cluster row drawn uncoloured with one label per cluster, dashed
#     separators from it up through the tree;
#   - legend rebuilt as two rows aligned to the heatmap edges;
#   - banner (crop of patient history-5.pdf) on top, or the patient name.
# Output: out/paper/<IND>_<tag>.png (300 dpi) [+ .pdf, vector except the banner]
# All sizes/fonts are the constants under "layout parameters".
# ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid); library(gtable) })

args <- commandArgs(trailingOnly = TRUE)
IND <- if (length(args)) args[1] else "BRCA-795"
TAG <- if (length(args) > 1) args[2] else "iter"
PDF <- "--pdf" %in% args
dir.create("out/paper/cache", showWarnings = FALSE, recursive = TRUE)

# ---- per-figure settings ---------------------------------------------------
# branch widths: hypergeometric test of the leaves of the branch's colour under
# that branch vs randomly sampled leaf colours (the app's simulated test)
BRANCH_P_TITLE <- "Branch p-value\n(hypergeometric)"
TITLES <- c("scIMPACT Mutation Calls" = "Mutation (scIMPACT)",
            "CNV data" = "Copy number", "p_val_new" = BRANCH_P_TITLE,
            "Bootstrap" = "Bootstrap", "Tumor / Non-tumor" = "Tumor / non-tumor",
            "Cell Type" = "Cell type", "Tissue Source" = "Tissue source",
            "Sampling Date" = "Sampling date", "NRAS p.Q61R survey" = "NRAS p.Q61R survey")
# banner: crop of patient history-5.pdf, px at 600 dpi
FIG <- list(
  "BRCA-795" = list(banner = list(page = 1, x = 250, y = 1150, w = 4800, h = 1100)),
  "BRCA-775" = list(banner = list(page = 1, x = 250, y = 2330, w = 4800, h = 1420)),
  "MM-127"   = list(banner = list(page = 2, x = 100, y = 110,  w = 4520, h = 1250)),
  "MM-412"   = list(banner = list(page = 2, x = 100, y = 1500, w = 4520, h = 1240)),
  "MM-423"   = list(banner = list(page = 2, x = 100, y = 2890, w = 4520, h = 1240))
)
cfg <- FIG[[IND]]; cfg$titles <- TITLES

# ---- layout parameters (mm, pt) ---------------------------------------------
W_MM     <- 180                  # double-column width
OUTER_MM <- 2                    # outer margin
H_TREE   <- 62                   # tree + tip labels
ROW_MM   <- 2.4                  # annotation row thickness (same for every patient)
FONT     <- "Nimbus Sans"        # Helvetica metrics, has a bold face (the default "sans" here has none)
H_CNV    <- 34
H_TITLE  <- 10; PT_TITLE <- 11; TITLE_COL <- "#C03830"   # title when there is no banner (banner red)
H_LEG    <- 34
TREE_LOG_C  <- 0.05             # pseudo-log depth scale (0 = linear)
BRANCH_W_SCALE <- 0.45           # branch width multiplier
BOOT_SCALE <- 0.6                # bootstrap triangles (tree + legend keys) size multiplier
RIGHT_MM <- 27                   # mm right of the last tip, for the row labels
TIP_PAD  <- 0.075                # depth units below the deepest tip for its label
SHOW_TIPS <- FALSE               # tip labels hidden in every paper figure (< 6 pt on all but the smallest trees)
PT_ROW   <- 6                    # heatmap row labels (all text >= 6 pt at 180 mm)
PT_MUT   <- 6                    # mutation row labels
PT_CHR   <- 6                    # chromosome numbers
CHR_MIN_MM <- 2.4                # mm between kept chromosome labels (6 pt ~ 2.1 mm)
DENSE_TIPS <- 100                # above this many cells: no per-cell grid lines in the rows
LEG_ROW1 <- c("Cluster", "Cell type", "Tumor / non-tumor", "Tissue source", "Phenotype")
KEY_MM   <- 2.6
LEG_GAP_MAX <- 30                # mm; short legend rows stay left-aligned
LW_GRID  <- 0.2                  # mm; heatmap grid + chromosome lines
PT_CLUSTER <- 6; CLUSTER_PREFIX <- "^Cluster "      # Cluster row labels
SEP_LTY <- "22"; SEP_LW <- 0.3; SEP_COL <- "grey45"   # cluster separators
PT_LEG_T <- 7; PT_LEG_X <- 6
# copy number. Absolute (sunshine: integer copies, 2 = diploid = white): one
# fixed colour per state, shared by every patient, top state open-ended.
# Relative (centred on 0 = white): symmetric blue-white-red, limits clipped at
# the 99th percentile of |value| and disclosed in the legend.
CN_ABS <- c("0" = "#08519C", "1" = "#9ECAE1", "2" = "#FFFFFF", "3" = "#FCBBA1",
            "4" = "#FB6A4A", "5" = "#CB181D", "6+" = "#67000D")
CN_REL <- c("darkblue", "white", "darkred")
CN_NA  <- "#8C8C8C"              # missing copy number: a grey the blue-white-red ramp never reaches
# legend label for the NA key, per legend (title after renaming); default "no data"
NA_LABELS <- c("Mutation (scIMPACT)" = "not called", "Cell type" = "unassigned",
               "LN type" = "not LN", "NRAS p.Q61R survey" = "not tested")
STAMP_COL <- "grey40"; PT_STAMP <- 6   # provenance stamp on non-final tags

# ---- load the plot (cached copy with the app's helper functions resolved) ---
# Drawing the legend keys needs the app's functions (debug_cat, draw_key_*).
source("individuals.R"); ez_load_app()
# app render with the paper settings (render_paper.R)
p <- readRDS(file.path("out/paper/cache", paste0(IND, "_plot.rds")))$plot

# ---- tree depth: pseudo-log so the short branches near the root get room ----
# p$data$x = depth - D (root at -D, deepest tip at 0). New depth
# D * log1p(d / c) / log1p(D / c) with c = TREE_LOG_C * D: root and deepest tip
# stay put, everything in between moves down. (Non-linear: note in the legend.)
if (TREE_LOG_C > 0) {
  D <- -min(p$data$x); d <- p$data$x + D; cc <- TREE_LOG_C * D
  p$data$x <- D * log1p(d / cc) / log1p(D / cc) - D
}
# ---- branch widths: thinner, same three p-value classes -----------------------
b0 <- ggplot_build(p)
ss <- b0$plot$scales$get_scales("size")
lv <- ss$get_limits(); old_w <- ss$map(lv)
p <- p + scale_size_manual(values = setNames(old_w * BRANCH_W_SCALE, lv), limits = lv, drop = FALSE)
message("[paper] branch widths ", paste(old_w, collapse = "/"), " -> ", paste(old_w * BRANCH_W_SCALE, collapse = "/"))

b <- ggplot_build(p)
kind <- vapply(p$layers, function(l) class(l$geom)[1], "")

# band limits along the depth axis (x before coord_flip)
xr <- function(i) range(unlist(lapply(i, function(k) {
  d <- b$data[[k]]; c(d$x, d$xmin, d$xmax, d$xend)
})), na.rm = TRUE)
# heatmap tile layers (ggnewscale renames all but the last heatmap's: GeomTile)
tiles <- which(kind %in% c("NewGeomTile", "GeomTile"))
# CNV heatmap = tile layer with many tiles per tip (mutation blocks have <= ~10) (absent when there is no RData)
n_tip <- sum(p$data$isTip)
cnv <- tiles[vapply(tiles, function(k) nrow(b$data[[k]]) > 25 * n_tip, TRUE)]
ann <- setdiff(tiles, cnv)
tipd <- b$data[[ann[1]]]                    # one tile per tip: tip positions
tree_r <- xr(1:4)
ann_r  <- range(unlist(lapply(ann, function(k) b$data[[k]][c("xmin", "xmax")])))
cnv_r  <- if (length(cnv)) range(unlist(b$data[[cnv]][c("xmin", "xmax")])) else c(NA, NA)
# coord_flip: panel x.range is the tip axis. Both scales are reversed, so coord
# limits (given in data units) are the negated transformed ranges; check once.
yl <- -rev(b$layout$panel_params[[1]]$x.range)
chk <- ggplot_build(p + coord_flip(xlim = c(-1, 0), ylim = yl, expand = FALSE))$layout$panel_params[[1]]
stopifnot(isTRUE(all.equal(chk$x.range, b$layout$panel_params[[1]]$x.range)), isTRUE(all.equal(chk$y.range, c(0, 1))))
# room for the row labels right of the heatmaps: widen the tip axis so that
# (pad + extra) tip units = RIGHT_MM, with mm per unit = Wp / (span + extra)
Wp <- W_MM - 2 * OUTER_MM
pad0 <- min(-tipd$y) - 0.5 - yl[1]     # tip units right of the first tip's edge now
extra <- max(0, (RIGHT_MM * diff(yl) - pad0 * Wp) / (Wp - RIGHT_MM))
yl[1] <- yl[1] - extra
row_w <- median(unlist(lapply(ann, function(k) b$data[[k]]$xmax - b$data[[k]]$xmin)))
H_ROWS <- ROW_MM * (diff(ann_r) + 0.008) / row_w
mm_per <- (W_MM - 2 * OUTER_MM) / diff(yl)      # mm per tip
if (!SHOW_TIPS) TIP_PAD <- 0.01                 # no label space under the tips
message(sprintf("[paper] tree %.3f..%.3f  rows %.3f..%.3f  cnv %.3f..%.3f", tree_r[1], tree_r[2],
                ann_r[1], ann_r[2], cnv_r[1], cnv_r[2]))

# ---- copy number: recolour the tiles from their values -----------------------
# The app pre-computes fill_color per tile (blue-white-red split at 2, each side
# stretched to the patient's own extreme), so gains and losses of one copy got
# different saturation and the colours differed between patients. Replaced here.
cn_leg_df <- NULL
if (length(cnv)) {
  d <- p$layers[[cnv]]$data; v <- d$value
  CN_IS_ABS <- all(abs(v - round(v)) < 1e-6, na.rm = TRUE) && min(v, na.rm = TRUE) >= 0
  if (CN_IS_ABS) {
    st <- ifelse(v >= length(CN_ABS) - 1, names(CN_ABS)[length(CN_ABS)], as.character(round(v)))
    d$fill_color <- ifelse(is.na(v), CN_NA, unname(CN_ABS[st]))
    cn_title <- "Copy number\n(absolute)"
    message(sprintf("[paper] copy number absolute, max %g, %.2f%% of bins >= %s", max(v, na.rm = TRUE),
                    100 * mean(v >= length(CN_ABS) - 1, na.rm = TRUE), names(CN_ABS)[length(CN_ABS)]))
  } else {
    lim <- quantile(abs(v), 0.99, na.rm = TRUE); mx <- max(abs(v), na.rm = TRUE)
    shades <- colorRampPalette(CN_REL)(1000)
    d$fill_color <- ifelse(is.na(v), CN_NA, shades[1 + round((pmin(pmax(v, -lim), lim) + lim) / (2 * lim) * 999)])
    cn_title <- if (mx > lim) sprintf("Relative copy number\nmax |%.3g| > %.3g CLIPPED", mx, lim) else "Relative copy number"
    message(sprintf("[paper] copy number relative, |max| %.3g, colour limit +-%.3g", mx, lim))
  }
  p$layers[[cnv]]$data <- d
  cn_has_na <- anyNA(v)
}

# ---- dashed cluster separators: Cluster row up through the tree ------------
# Cluster row = the annotation row nearest the tree. Boundaries where the tile
# colour changes between neighbouring tips. Built coordinates are the negated
# data coordinates (both scales reversed), so negate back for the new layer.
cl <- ann[which.max(vapply(ann, function(k) mean(b$data[[k]]$x), 1))]
cd <- b$data[[cl]]; cd <- cd[order(cd$y), ]
fcol <- names(cd)[grepl("^fill", names(cd))][1]
chg <- which(cd[[fcol]][-1] != cd[[fcol]][-nrow(cd)])
sep <- data.frame(y = -(cd$y[chg] + cd$y[chg + 1]) / 2)
if (nrow(sep)) { sep$x <- -min(cd$xmin); sep$xend <- -(tree_r[2] + 0.005) }
message("[paper] cluster separators at ", length(chg), " boundaries")
sep_layer <- if (nrow(sep)) geom_segment(data = sep, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE,
                                         linetype = SEP_LTY, linewidth = SEP_LW, colour = SEP_COL)

# ---- text sizes (ggplot size is mm) ----------------------------------------
for (k in seq_along(p$layers)) {
  l <- p$layers[[k]]
  if (kind[k] %in% c("GeomTextGGtree", "GeomText")) l$aes_params$family <- FONT
  if (kind[k] == "GeomTextGGtree" && !SHOW_TIPS) l$aes_params$alpha <- 0
  if (kind[k] == "GeomText") {
    lab <- b$data[[k]]$label
    l$aes_params$size <- if (all(grepl("^([0-9]+|X|Y)$", lab))) PT_CHR / .pt
                         else if (length(lab) > 1) PT_MUT / .pt else PT_ROW / .pt
  }
}

# The app's (invisible) Bootstrap-legend point layer has size 5 and joins every
# legend, which inflates all keys to 5 mm; its triangle keys use fixed pt sizes.
for (k in which(kind == "GeomPoint")) p$layers[[k]]$aes_params$size <- 0.5

# bootstrap triangles smaller: the tree's three point layers (app sizes 4.5/3.5/2.5)
# and the legend keys, which the app draws as its own triangles at 5/4/3 pt
for (k in which(kind == "GeomPointGGtree")) {
  l <- p$layers[[k]]
  if (identical(l$aes_params$shape, 24)) l$aes_params$size <- l$aes_params$size * BOOT_SCALE
}
boot_key <- function(data, params, size) {
  sz <- c("1" = 5, "2" = 4, "3" = 3)[as.character(data$shape)] * BOOT_SCALE
  if (!length(sz) || is.na(sz)) return(grid::nullGrob())
  grid::polygonGrob(x = unit(0.5, "npc") + unit(c(-sz, sz, 0), "pt"),
                    y = unit(0.5, "npc") + unit(c(-sz * 0.6, -sz * 0.6, sz * 0.8), "pt"),
                    gp = gpar(fill = "grey36", col = "grey20", alpha = 0.5))
}
for (k in which(kind == "GeomPoint")) if ("bootstrap_shape" %in% vapply(p$layers[[k]]$mapping, rlang::as_label, ""))
  p$layers[[k]]$geom$draw_key <- boot_key    # this layer's own geom (the app made it via key_glyph)

# heatmap grid lines and chromosome separators: thinner
hm_first <- min(which(kind %in% c("NewGeomTile", "GeomTile")))
for (k in which(kind == "GeomSegment")) if (k > hm_first) p$layers[[k]]$aes_params$linewidth <- LW_GRID
# dense trees: per-cell grid lines (segments along the depth axis) turn the rows
# into grey stripes - keep only each block's outer two
if (n_tip > DENSE_TIPS) for (k in which(kind == "GeomSegment")) {
  d <- p$layers[[k]]$data
  if (k <= hm_first || !is.data.frame(d) || !all(c("x", "xend", "y", "yend") %in% names(d))) next
  along <- abs(d$y - d$yend) < 1e-9 & abs(d$x - d$xend) > 1e-9
  if (sum(along) < 3) next
  edge <- along & (abs(d$y - min(d$y[along])) < 1e-9 | abs(d$y - max(d$y[along])) < 1e-9)
  p$layers[[k]]$data <- d[!along | edge, ]
}
# chromosome numbers: drop a label closer than CHR_MIN_MM to the previous kept one
# (depth units per mm in the CNV band: its range + the band padding over H_CNV)
row_lab_y <- p$layers[[which(kind == "GeomText")[1]]]$data$y[1]    # data y of the row labels
CHR_MIN_GAP <- if (length(cnv)) CHR_MIN_MM * (diff(cnv_r) + 0.024) / H_CNV else 0
for (k in which(kind == "GeomText")) {
  d <- p$layers[[k]]$data
  if (!is.data.frame(d) || !all(grepl("^([0-9]+|X|Y)$", d$label))) next
  pos <- d$x; o <- order(-pos); keep_i <- o[1]; last <- pos[o[1]]
  for (i in o[-1]) if (abs(pos[i] - last) >= CHR_MIN_GAP) { keep_i <- c(keep_i, i); last <- pos[i] }
  message("[paper] chromosome labels dropped: ", paste(d$label[-keep_i], collapse = " "))
  d <- d[sort(keep_i), ]
  # the app right-aligns them 3 tips past the heatmap - at 6 pt they reach into
  # it on dense trees; put them in the row-label column, left-aligned
  d$y <- row_lab_y; p$layers[[k]]$aes_params$hjust <- 0
  p$layers[[k]]$data <- d
}

# ---- legend titles ------------------------------------------------------------
for (s in p$scales$scales) {
  nm <- tryCatch(as.character(s$name), error = function(e) "")
  if (length(nm) == 1 && nm %in% names(cfg$titles)) s$name <- cfg$titles[[nm]]
}
# NA keys say what missing means ("not called", "unassigned", ...); dates as YYYY-MM
relabel <- function(s, na_lab, fix = identity) {
  old <- s$labels; force(na_lab); force(fix)      # evaluated now, not when the legend is drawn
  s$labels <- function(x) {
    out <- if (is.function(old)) old(x) else if (inherits(old, "waiver")) as.character(x)
           else if (!is.null(names(old))) { r <- unname(old[as.character(x)]); ifelse(is.na(r), as.character(x), r) }
           else if (length(old) == length(x)) as.character(old) else as.character(x)
    out <- fix(out)
    out[is.na(x) | x %in% "NA" | out %in% "NA"] <- na_lab
    out
  }
}
for (s in p$scales$scales) {
  if (!inherits(s, "ScaleDiscrete")) next
  nm <- tryCatch(as.character(s$name), error = function(e) "")
  if ("colour" %in% s$aesthetics) nm <- "Cell type"      # branch colours; title pinned in p$guides
  if (length(nm) != 1 || is.na(nm) || !nzchar(nm) || nm == "Cluster") next
  relabel(s, if (nm %in% names(NA_LABELS)) NA_LABELS[[nm]] else "no data",
          if (nm == "Sampling date") function(v) sub("^([0-9]{4})_([0-9]{2})$", "\\1-\\2", v) else identity)
}

# tree-branch colour legend duplicates the Cell type row legend
p <- p + labs(size = BRANCH_P_TITLE)
# row labels drawn next to the heatmaps use the same wording as the legend titles
for (k in which(kind == "GeomText")) {
  d <- p$layers[[k]]$data
  if (is.data.frame(d) && "label" %in% names(d)) {
    hit <- d$label %in% names(cfg$titles)
    d$label[hit] <- unname(cfg$titles[d$label[hit]]); p$layers[[k]]$data <- d
  }
}

# Cluster row label carries the total number of cells
for (k in which(kind == "GeomText")) {
  d <- p$layers[[k]]$data
  if (is.data.frame(d) && "label" %in% names(d) && any(d$label == "Cluster"))
    { d$label[d$label == "Cluster"] <- sprintf("Cluster (%d cells)", n_tip); p$layers[[k]]$data <- d }
}

# ---- Cluster row: no colours, an outlined row with one label per cluster ------
# The separators run through it. Cluster names come from the row's fill scale
# (colour -> value); a cluster split into several runs is labelled on its longest.
fs <- b$plot$scales$get_scales(fcol)
lv <- fs$get_limits(); col2cl <- setNames(lv, fs$map(lv))
run <- cumsum(c(TRUE, cd[[fcol]][-1] != cd[[fcol]][-nrow(cd)]))
runs <- do.call(rbind, lapply(split(seq_len(nrow(cd)), run), function(i)
  data.frame(cluster = col2cl[[cd[[fcol]][i[1]]]], n = length(i), y = -mean(cd$y[i]))))
n_cl <- tapply(runs$n, runs$cluster, sum)
runs <- runs[order(-runs$n), ]; runs <- runs[!duplicated(runs$cluster), ]
runs$label <- sub(CLUSTER_PREFIX, "", runs$cluster)
# cell count in the label where it fits in the run (bold 6 pt ~ 0.62 * pt * 0.353 mm per char)
fits <- function(txt, n) nchar(txt) * 0.62 * PT_CLUSTER * 0.353 + 1 <= n * mm_per
with_n <- sprintf("%s (%d)", runs$label, n_cl[runs$cluster])
runs$label <- ifelse(fits(with_n, runs$n), with_n, runs$label)
message("[paper] cluster labels: ", paste(runs$label, collapse = " "))
cl_x <- -mean(range(c(cd$xmin, cd$xmax)))
box <- data.frame(xmin = -max(cd$xmax), xmax = -min(cd$xmin),
                  ymin = -max(cd$y) - 0.5, ymax = -min(cd$y) + 0.5)
cl_layers <- list(
  geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = LW_GRID),
  geom_text(data = runs, aes(x = cl_x, y = y, label = label), inherit.aes = FALSE,
            size = PT_CLUSTER / .pt, fontface = "bold", family = FONT))
# drop the Cluster tiles and their grid (the two segment layers drawn with them)
stopifnot(all(kind[cl + 1:2] == "GeomSegment"))
p$layers <- p$layers[-(cl + 0:2)]; kind <- kind[-(cl + 0:2)]

p$layers <- c(if (!is.null(sep_layer)) list(sep_layer), p$layers, cl_layers)   # separators under the branches and tiles
kind <- c(if (!is.null(sep_layer)) "sep", kind, "cl_box", "cl_text")

# ---- bands ------------------------------------------------------------------
band <- function(lim, pad = c(0, 0)) {
  p + coord_flip(xlim = -c(lim[2] + pad[2], lim[1] - pad[1]), ylim = yl, expand = FALSE, clip = "on") +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), text = element_text(family = FONT))
}
pad <- 0.004
p_tree <- band(c(-TIP_PAD, tree_r[2]), c(0, 0.01))
p_rows <- band(ann_r, c(pad, pad))
p_cnv  <- if (length(cnv)) band(cnv_r, c(0.02, pad)) else NULL    # extra below for the last chromosome label

# ---- legends: pull each guide out of the app's guide box, lay out as a grid ---
leg_theme <- theme(legend.title = element_text(size = PT_LEG_T, face = "bold", hjust = 0, family = FONT, margin = margin(b = 1, unit = "mm")),
                   legend.text = element_text(size = PT_LEG_X, family = FONT, margin = margin(l = 0.8, unit = "mm")),
                   legend.key.size = unit(KEY_MM, "mm"), legend.key.width = unit(KEY_MM, "mm"),
                   legend.key.height = unit(KEY_MM, "mm"), legend.key.spacing.y = unit(0.5, "mm"),
                   legend.key.spacing.x = unit(1, "mm"),
                   legend.margin = margin(0, 0, 0, 0), legend.position = "bottom",
                   legend.direction = "vertical")
# longer colour bar for the copy-number legend: the continuous fill scale (the app
# draws that legend from its own scale; the CNV tile layer's guide is hidden).
# Only that scale - a blanket guides(fill = ...) would drop a discrete legend on "fill".
cont_fill <- Filter(function(sc) inherits(sc, "ScaleContinuous") && any(grepl("^fill", sc$aesthetics)),
                    b$plot$scales$scales)
# The app's copy-number colourbar is dropped and replaced by a legend of the
# recoloured tiles (built from a stand-alone plot, same theme).
cnv_guide <- if (length(cont_fill)) do.call(guides, setNames(list("none"), cont_fill[[1]]$aesthetics[1]))
g <- ggplotGrob(p + leg_theme + cnv_guide)
box <- g$grobs[[grep("guide-box-bottom", g$layout$name)]]
guides <- box$grobs[box$layout$name == "guides"]
if (length(cnv)) {
  cn_plot <- if (CN_IS_ABS) {
    st <- c(names(CN_ABS), if (cn_has_na) NA)
    ggplot(data.frame(s = factor(st, levels = names(CN_ABS))), aes(0, 0, fill = s)) +
      geom_tile(colour = "grey60", linewidth = 0.2) +
      scale_fill_manual(values = CN_ABS, na.value = CN_NA, drop = FALSE, name = cn_title,
                        labels = function(x) ifelse(is.na(x), "no data", x)) +
      guides(fill = guide_legend(ncol = 2))
  } else {
    ggplot(data.frame(v = c(-lim, lim)), aes(0, 0, fill = v)) + geom_tile() +
      scale_fill_gradientn(colours = CN_REL, limits = c(-lim, lim), oob = scales::squish, name = cn_title) +
      guides(fill = guide_colourbar(theme = theme(legend.key.height = unit(14, "mm"), legend.key.width = unit(2.2, "mm"))))
  }
  cg <- ggplotGrob(cn_plot + theme_void() + leg_theme)
  cb <- cg$grobs[[grep("guide-box-bottom", cg$layout$name)]]
  guides <- c(guides, cb$grobs[cb$layout$name == "guides"])
}
guide_title <- function(gt) {
  i <- grep("^title", gt$layout$name)
  if (!length(i)) return("")
  tx <- tryCatch(gt$grobs[[i[1]]]$children[[1]]$children[[1]]$label, error = function(e) NULL)
  if (is.null(tx)) tx <- tryCatch(gt$grobs[[i[1]]]$children[[1]]$label, error = function(e) "")
  paste(tx, collapse = " ")
}
justify_top <- function(x) {
  x$vp <- viewport(x = 0, y = 1, just = c("left", "top"), width = sum(x$widths), height = sum(x$heights))
  x
}
titles <- vapply(guides, guide_title, "")
# thin outline on every key of the tile legends, as the rows have one - a white
# key ("not called", "no data") is otherwise invisible
outline_keys <- function(gt) {
  i <- grep("^key-[0-9]+-[0-9]+-bg$", gt$layout$name)
  for (j in i) gt <- gtable::gtable_add_grob(gt, rectGrob(gp = gpar(fill = NA, col = "grey55", lwd = 0.5)),
                                             t = gt$layout$t[j], l = gt$layout$l[j], name = paste0("outline-", j), z = Inf)
  gt
}
tile_leg <- !sub("\n.*", "", titles) %in% c("Cell type", "Bootstrap", "Branch p-value", "") & !grepl("^Copy number|^Relative copy", titles)
guides[tile_leg] <- lapply(guides[tile_leg], outline_keys)
message("[paper] guides: ", paste(titles, collapse = " | "))
keep <- nzchar(titles)
order_pref <- c("Cluster", "Cell type", "Tumor / non-tumor", "Tissue source", "Phenotype",
                "Mutation (scIMPACT)", "Sampling date", "WGD", "Copy number", "Relative copy number",
                "Bootstrap", "Branch p-value")
ord <- order(match(sub("\n.*", "", titles[keep]), order_pref, nomatch = 99))
guides <- guides[keep][ord]
titles <- titles[keep][ord]
# heatmap left/right edges in mm from the figure's left edge (bands span W_MM)
tip_rng <- range(-tipd$y)                      # data units of first / last tip
x_left  <- (yl[2] - (tip_rng[2] + 0.5)) * mm_per
x_right <- (yl[2] - (tip_rng[1] - 0.5)) * mm_per
message(sprintf("[paper] heatmap spans %.1f..%.1f mm", x_left, x_right))
# one legend row: first guide at the heatmap's left edge, last one ending at its right edge
legend_row <- function(gs) {
  w <- vapply(gs, function(x) convertWidth(sum(x$widths), "mm", valueOnly = TRUE), 1)
  gap <- min(LEG_GAP_MAX, max(2, (x_right - x_left - sum(w)) / max(1, length(gs) - 1)))
  ws <- c(x_left, as.vector(rbind(w, gap))[-2 * length(gs)])
  grobs <- c(list(nullGrob()), unlist(lapply(seq_along(gs), function(i)
    if (i < length(gs)) list(justify_top(gs[[i]]), nullGrob()) else list(justify_top(gs[[i]]))), recursive = FALSE))
  gridExtra::arrangeGrob(grobs = c(grobs, list(nullGrob())), nrow = 1,
                         widths = unit.c(unit(ws, "mm"), unit(1, "null")))
}
row1 <- titles %in% LEG_ROW1
lh <- function(gs) max(do.call(unit.c, lapply(gs, function(x) convertHeight(sum(x$heights), "mm"))))
legend_grob <- gridExtra::arrangeGrob(legend_row(guides[row1]), legend_row(guides[!row1]), ncol = 1,
                                      heights = unit.c(lh(guides[row1]) + unit(2, "mm"), lh(guides[!row1])))

H_LEG <- convertHeight(sum(legend_grob$heights), "mm", valueOnly = TRUE) + 3
message(sprintf("[paper] legend height %.1f mm", H_LEG))

# ---- banner (patient-history timeline), or just the patient name ----------------
no_margin <- theme(plot.margin = margin(0, 0, 0, 0))
if (is.null(cfg$banner)) {
  # no timeline (single sampling time point): the name as the banners print it
  H_BANNER <- H_TITLE
  p_ban <- wrap_elements(full = textGrob(IND, x = unit(0, "npc"), hjust = 0,
                                         gp = gpar(fontsize = PT_TITLE, fontface = 2, fontfamily = FONT, col = TITLE_COL))) + no_margin
} else {
bn <- cfg$banner
ban_png <- file.path("out/paper/cache", paste0(IND, "_banner.png"))
if (!file.exists(ban_png)) {
  system2("pdftoppm", c("-r 600", "-f", bn$page, "-l", bn$page, "-x", bn$x, "-y", bn$y, "-W", bn$w, "-H", bn$h,
                        "-png", "-singlefile", shQuote(file.path(KLEIN, "patient history-5.pdf")),
                        shQuote(sub("\\.png$", "", ban_png))))
}
ban <- png::readPNG(ban_png)
H_BANNER <- (W_MM - 2 * OUTER_MM) * nrow(ban) / ncol(ban)     # full width, aspect kept
p_ban <- wrap_elements(full = rasterGrob(ban, interpolate = TRUE)) + no_margin
}

# ---- compose ------------------------------------------------------------------
if (is.null(p_cnv)) H_CNV <- 0
parts <- Filter(Negate(is.null), list(p_ban, p_tree, p_rows, p_cnv, wrap_elements(full = legend_grob) + no_margin))
hts <- c(H_BANNER, H_TREE, H_ROWS, if (H_CNV > 0) H_CNV, H_LEG)
fig <- wrap_plots(parts, ncol = 1) + plot_layout(heights = unit(hts, "mm"))
# no per-part margins (keeps legend / banner offsets in the bands' coordinates), one outer margin
# provenance stamp (masterlist, tag, date) on every tag except "final"
stamp <- if (TAG != "final") sprintf("%s | %s | %s | %s", IND, TAG, basename(CSV), format(Sys.Date()))
fig <- fig + plot_annotation(caption = stamp,
                             theme = theme(plot.margin = margin(OUTER_MM, OUTER_MM, OUTER_MM, OUTER_MM, "mm"),
                                           plot.caption = element_text(size = PT_STAMP, colour = STAMP_COL, family = FONT)))
H_MM <- H_BANNER + H_TREE + H_ROWS + H_CNV + H_LEG + 2 * OUTER_MM + if (is.null(stamp)) 0 else 3
out <- file.path("out/paper", sprintf("%s_%s.png", IND, TAG))
ggsave(out, fig, width = W_MM, height = H_MM, units = "mm", dpi = 300, bg = "white")
if (PDF) ggsave(sub("\\.png$", ".pdf", out), fig, width = W_MM, height = H_MM, units = "mm", device = cairo_pdf)
message("[paper] wrote ", out, sprintf(" (%d x %.0f mm)", W_MM, H_MM))
