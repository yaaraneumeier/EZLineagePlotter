# Preflight: packages at the expected versions, pdftoppm, Nimbus Sans with a
# real bold face, vendored app and data present. Run from the repo root.
#   Rscript env/check_env.R
ok <- TRUE
bad <- function(...) { message("FAIL: ", ...); ok <<- FALSE }
want <- c(ggplot2 = "4.0.2", S7 = "0.2.1", scales = "1.4.0", gtable = "0.3.6", ggtree = "3.10.0",
          treeio = "1.26.0", ggnewscale = "0.5.0", patchwork = "1.3.0", ape = "5.8.1", shiny = "1.8.1.1",
          yaml = "2.3.9", gridExtra = "2.3", png = "0.1.8", systemfonts = "1.1.0")
for (p in names(want)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  if (is.na(v)) bad(p, " not installed")
  else if (v != want[[p]]) message("warn: ", p, " ", v, " (figures made with ", want[[p]], ")")
}
other <- c("shinyjs", "shinyWidgets", "shinydashboard", "shinyBS", "colourpicker", "DT", "data.table", "digest",
           "RColorBrewer", "viridis", "cowplot", "ggforce", "jpeg", "xml2", "stringr", "dplyr", "rlang",
           "infotheo", "combinat", "later", "tidyverse")
for (p in other) if (!requireNamespace(p, quietly = TRUE)) bad(p, " not installed")
if (!nzchar(Sys.which("pdftoppm"))) bad("pdftoppm not on PATH (conda: poppler)")
if (R.version$major != "4" || R.version$minor != "3.3") message("warn: R ", R.version.string, " (figures made with 4.3.3)")
f <- systemfonts::match_fonts("Nimbus Sans", weight = "bold")$path
if (!grepl("NimbusSans-Bold", f)) bad("'Nimbus Sans' bold resolves to ", f, " - set FONTCONFIG_FILE=fonts/fonts.conf (source env/activate.sh)")
for (x in c("vendor/EZLineagePlotter/headless/ez_headless.R", "data/masterlist_lineage_tree_v84_20260923.csv",
            "data/patient history-5.pdf", "data/trees2025_ordering", "data/mutation_labels.tsv"))
  if (!file.exists(x)) bad(x, " missing (run from the repo root; data/ comes from the data archive, see README)")
if (ok) message("check_env: all good") else quit(status = 1)
