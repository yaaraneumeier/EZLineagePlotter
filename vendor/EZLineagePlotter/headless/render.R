# ============================================================================
# render.R - command-line entry point for ez_headless.R
#
#   Rscript headless/render.R --yaml cfg.yaml --tree tree.newick --csv data.csv \
#       --out Tree___name.pdf [--rdata cnv.RData] [--annot annot.RData] \
#       [--patch patch.yaml] [--class-column Cell.type.2] [--class-title "Cell type"] \
#       [--seed 1]
#
# --patch: YAML with recovered settings the app's config export does not save,
#   e.g.  heatmaps: {1: {label_mapping: {ZCA: Cluster}}, 9: {cnv_chr_lines: 'yes'}}
#   (1-based heatmap index -> fields, same names as the app's YAML import).
# The output format follows the --out extension (pdf, png, svg).
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- list(`class-column` = "Cell.type.2", `class-title` = "Cell type", seed = "1")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || i == length(args)) stop("bad argument: ", args[i])
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
for (req in c("yaml", "tree", "csv", "out")) if (is.null(opt[[req]])) stop("missing --", req)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
if (!nzchar(Sys.getenv("EZ_APP_DIR"))) Sys.setenv(EZ_APP_DIR = dirname(script_dir))
source(file.path(script_dir, "ez_headless.R"))

patch <- if (!is.null(opt$patch)) yaml::read_yaml(opt$patch) else NULL  # before ez_render() changes dir
res <- ez_render(opt$yaml, opt$tree, opt$csv, opt$out,
                 rdata = opt$rdata, annot = opt$annot,
                 patch = patch,
                 class_column = opt$`class-column`, class_title = opt$`class-title`,
                 seed = as.integer(opt$seed))
if (length(res$unmatched_classes))
  message("classes without a paper-palette colour (rainbow default): ", paste(res$unmatched_classes, collapse = ", "))
if (!isTRUE(res$ok)) quit(status = 1)
message("wrote ", res$out_file)
