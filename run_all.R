# ============================================================================
# run_all.R - re-render every figure under "Trees 04-2026" / "Trees 07-2026"
# from its saved config with ez_headless.R.
#
#   Rscript run_all.R            # all 19 configs
#   Rscript run_all.R MM-127     # only configs whose file name matches
#
# Output: out/<Trees folder>/<same file name as the reference figure>
#         out/run_log.tsv (out/run_log_<filter>.tsv when a filter is given)
# ============================================================================
source("individuals.R")

# config -> reference figure(s) it produced (same folder)
reference_for <- function(cfg) {
  stem <- sub("_config( \\([0-9]+\\))?\\.yaml$", "", basename(cfg))
  refs <- list.files(dirname(cfg), pattern = "^Tree___.*\\.(pdf|png)$")
  norm <- function(x) gsub("[ _]+|\\([0-9]+\\)", "", tools::file_path_sans_ext(sub("^Tree___", "", x)))
  refs[norm(refs) == norm(stem)]
}

patches <- yaml::read_yaml("patches.yaml")
args <- commandArgs(trailingOnly = TRUE)
configs <- list.files(file.path(KLEIN, c("Trees 04-2026", "Trees 07-2026")), pattern = "_config.*\\.yaml$", full.names = TRUE)
if (length(args)) configs <- configs[grepl(args[1], basename(configs), fixed = TRUE)]

ez_load_app()
log_rows <- list()
for (cfg in configs) {
  ind <- names(INDIVIDUALS)[vapply(names(INDIVIDUALS), function(p) startsWith(basename(cfg), p), logical(1))]
  if (length(ind) != 1) { message("[run_all] no individual for ", basename(cfg)); next }
  spec <- INDIVIDUALS[[ind]]
  y <- yaml::read_yaml(cfg)
  uses_rdata <- any(vapply(y[["visual definitions"]]$heatmaps, function(h) identical(h$data_source, "rdata"), logical(1)))
  refs <- reference_for(cfg)
  fmt <- y[["Individual general definitions"]]$out_file$file_type
  out_name <- if (length(refs)) sub("\\.(pdf|png)$", paste0(".", fmt), refs[grepl(paste0("\\.", fmt, "$"), refs)][1]) else paste0("Tree___", basename(cfg), ".", fmt)
  out_file <- file.path("out", basename(dirname(cfg)), out_name)

  message("\n[run_all] ", basename(cfg), " -> ", out_file)
  t0 <- Sys.time()
  res <- tryCatch(
    ez_render(cfg, file.path(TREES, spec$tree), CSV, out_file,
              rdata = if (uses_rdata) file.path(RDATA, spec$rdata) else NULL,
              annot = if (uses_rdata) file.path(RDATA, spec$annot) else NULL,
              patch = mutation_label_patch(cfg, patches[[basename(cfg)]])),
    error = function(e) list(ok = FALSE, error = conditionMessage(e), unmatched_classes = character()))
  log_rows[[length(log_rows) + 1]] <- data.frame(
    config = basename(cfg), folder = basename(dirname(cfg)), individual = ind,
    reference = paste(refs, collapse = ";"), output = out_file, ok = isTRUE(res$ok),
    seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
    patched = !is.null(patches[[basename(cfg)]]), rdata = if (uses_rdata) spec$rdata else "",
    unmatched_classes = paste(res$unmatched_classes, collapse = ";"),
    error = if (is.null(res$error)) "" else res$error, stringsAsFactors = FALSE)
}
log_df <- do.call(rbind, log_rows)
log_file <- if (length(args)) sprintf("out/run_log_%s.tsv", gsub("[^A-Za-z0-9_-]", "_", args[1])) else "out/run_log.tsv"
write.table(log_df, log_file, sep = "\t", quote = FALSE, row.names = FALSE)
print(log_df[, c("config", "ok", "seconds", "unmatched_classes", "error")])
