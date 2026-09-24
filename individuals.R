# Shared inputs for run_all.R / zca_trees.R: data locations and, per individual
# prefix used in the config names, its newick and CNV RData (+ Annot) files.

# Self-contained layout (run every script from the repo root):
#   vendor/EZLineagePlotter/   the app (branch headless-render @ 3a1b0c5), headless driver in headless/
#   data/                      inputs, same relative layout as the klein/ project folder
# Override with EZ_ROOT (repo root), EZ_APP_DIR (another app checkout) or EZ_DATA (another data folder).
EZ_ROOT <- normalizePath(Sys.getenv("EZ_ROOT", "."), mustWork = TRUE)
if (!nzchar(Sys.getenv("EZ_APP_DIR"))) Sys.setenv(EZ_APP_DIR = file.path(EZ_ROOT, "vendor", "EZLineagePlotter"))
Sys.setenv(EZ_APP_DIR = normalizePath(Sys.getenv("EZ_APP_DIR"), mustWork = TRUE))
source(file.path(Sys.getenv("EZ_APP_DIR"), "headless", "ez_headless.R"))

KLEIN <- normalizePath(Sys.getenv("EZ_DATA", file.path(EZ_ROOT, "data")), mustWork = TRUE)
CSV   <- file.path(KLEIN, "masterlist_lineage_tree_v84_20260923.csv")
TREES <- file.path(KLEIN, "trees2025_ordering")
RDATA <- file.path(KLEIN, "patients_samples_by_name")

# individual prefix in the config name -> newick, CNV RData
INDIVIDUALS <- list(
  "BRCA-775"   = list(tree = "775_13_bootstrapped_rerooted.newick",
                      rdata = "775-13/RData/results_CNV_tool_final.RData",
                      annot = "775-13/RData/Data_Seg_Copy.RData"),
  "BRCA-795"   = list(tree = "795_09_bootstrapped_rerooted.newick",
                      rdata = "795-09/RData/results_CNV_tool_final.RData",
                      annot = "795-09/RData/Data_Seg_Copy.RData"),
  "BRCA-841"   = list(tree = "841_12_bootstrapped_rerooted.newick"),
  "MM-127"     = list(tree = "mm15_127_bootstrapped_rerooted.newick",
                      rdata = "MM15-127/RData/results_CNV_tool_final.RData",
                      annot = "MM15-127/RData/Data_Seg_Copy.RData"),
  "MM-412"     = list(tree = "mm16_412_bootstrapped_rerooted.newick",
                      rdata = "MM16-412/RData/results_CNV_tool_final.RData",
                      annot = "MM16-412/RData/Data_Seg_Copy.RData"),
  "MM-423"     = list(tree = "mm16_423_bootstrapped_rerooted.newick",
                      rdata = "MM16-423/RData/results_CNV_tool_final.RData",
                      annot = "MM16-423/RData/Data_Seg_Copy.RData"),
  "NSCLC-0267" = list(tree = "bc15_0267_bootstrapped_rerooted.newick",
                      rdata = "BC15-0267/RData/results_CNV_tool_final.RData",
                      annot = "BC15-0267/RData/Data_Seg_Copy.RData"),
  "NSCLC-0401" = list(tree = "bc16_0401_bootstrapped_rerooted.newick"),
  "NSCLC-0545" = list(tree = "bc16_0545_bootstrapped_rerooted.newick",
                      rdata = "BC16-0545/RData/results_CNV_tool_final.RData",
                      annot = "BC16-0545/RData/Data_Seg_Copy.RData"),
  "NSCLC-2680" = list(tree = "bc14_2680_bootstrapped_rerooted.newick",
                      rdata = "BC14-2680/RData/results_CNV_tool_final.RData",
                      annot = "BC14-2680/RData/Data_Seg_Copy.RData")
)


# ---- mutation row labels ("ESR1 p.S463P") -----------------------------------
# From the scIMPACT classifier workbook via mutation_labels.py -> mutation_labels.tsv.
# Columns named SYMBOL___chr___pos___ref___alt[___id] are matched on the first 5 parts.
MUT_LABELS <- local({
  t <- read.delim(file.path(KLEIN, "mutation_labels.tsv"), colClasses = "character")
  setNames(t$label, t$key5)
})
# Add label_mapping (row_label_source: mapping) to every heatmap whose columns are
# all mutation keys; indices refer to the config's own heatmap list, so this can
# be combined with other patch entries.
mutation_label_patch <- function(cfg_path, patch = NULL) {
  hms <- yaml::read_yaml(cfg_path)[["visual definitions"]]$heatmaps
  if (is.null(patch)) patch <- list()
  for (i in seq_along(hms)) {
    cols <- unlist(hms[[i]]$columns)
    if (!length(cols) || !all(grepl("^[A-Za-z0-9.]+___[0-9XY]+___[0-9]+___", cols))) next
    key <- vapply(strsplit(cols, "___", fixed = TRUE), function(p) paste(p[1:5], collapse = "___"), "")
    lab <- unname(MUT_LABELS[key])
    if (anyNA(lab)) warning(basename(cfg_path), ": no label for ", paste(cols[is.na(lab)], collapse = ", "))
    lab[is.na(lab)] <- cols[is.na(lab)]
    e <- patch$heatmaps[[as.character(i)]]; if (is.null(e)) e <- list()
    e$row_label_source <- "mapping"
    e$show_row_labels <- "yes"
    e$label_mapping <- as.list(setNames(lab, cols))
    patch$heatmaps[[as.character(i)]] <- e
  }
  patch
}
