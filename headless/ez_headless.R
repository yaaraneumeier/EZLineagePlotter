# ============================================================================
# ez_headless.R - render an EZLineagePlotter figure from a saved config YAML
# without opening the Shiny UI.
#
# Drives the real app server (EZlineagePlotter56.R, v.16) through
# shiny::testServer and replays the upload / import / classify / download
# steps a user performs in the browser. No plotting code is copied.
#
# Gaps in the saved YAMLs that this fills in:
#   - tree / CSV / RData paths (YAML points at deleted /tmp uploads)
#   - the tree classification (never written to the YAML): the column given in
#     `class_column`, coloured with the app's PAPER_PALETTE
#
# Usage (see README.md, render.R):
#   source("headless/ez_headless.R")
#   ez_render(yaml, tree, csv, out_file, rdata = NULL, annot = NULL, patch = NULL,
#             class_column = "Cell.type.2")
# ============================================================================

# App directory: $EZ_APP_DIR, else the repository root (parent of headless/).
APP_DIR <- local({
  env <- Sys.getenv("EZ_APP_DIR")
  if (nzchar(env)) return(normalizePath(env))
  this <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (!is.na(this)) return(dirname(dirname(this)))
  normalizePath("..")
})

suppressPackageStartupMessages({
  library(shiny)
  library(xml2)
})

# ---- load the app (everything except shinyApp(); aricode is unused) --------
ez_load_app <- function(app_dir = APP_DIR) {
  lines <- readLines(file.path(app_dir, "EZlineagePlotter56.R"), warn = FALSE)
  lines[grep("^shinyApp\\(|^library\\(aricode\\)", lines)] <- ""
  old <- setwd(app_dir); on.exit(setwd(old))
  suppressPackageStartupMessages(eval(parse(text = lines), envir = globalenv()))
  invisible(TRUE)
}

# ---- initial input values, read from the static UI ---------------------------
# testServer starts every input as NULL; the browser would start them at the
# UI defaults. Walk the rendered HTML and collect id -> default value.
ez_ui_defaults <- function(ui_obj) {
  doc <- read_html(as.character(ui_obj))
  out <- list()
  num_or_chr <- function(x) { n <- suppressWarnings(as.numeric(x)); if (!is.na(n)) n else x }

  for (n in xml_find_all(doc, "//input[@id]")) {
    id <- xml_attr(n, "id"); type <- tolower(xml_attr(n, "type") %||% "text")
    cls <- xml_attr(n, "class") %||% ""
    if (type == "file" || grepl("js-range-slider", cls)) next
    if (type == "checkbox") {
      out[[id]] <- !is.na(xml_attr(n, "checked"))
    } else if (type == "number") {
      out[[id]] <- num_or_chr(xml_attr(n, "value") %||% NA)
    } else if (type %in% c("text", "password", "")) {
      out[[id]] <- xml_attr(n, "value") %||% ""
    }
  }
  for (n in xml_find_all(doc, "//input[contains(@class,'js-range-slider')]")) {
    from <- as.numeric(xml_attr(n, "data-from")); to <- xml_attr(n, "data-to")
    out[[xml_attr(n, "id")]] <- if (!is.na(to) && (xml_attr(n, "data-type") %||% "") == "double") c(from, as.numeric(to)) else from
  }
  for (n in xml_find_all(doc, "//div[contains(@class,'shiny-input-radiogroup')][@id]")) {
    chk <- xml_find_first(n, ".//input[@checked]")
    if (!is.na(chk)) out[[xml_attr(n, "id")]] <- xml_attr(chk, "value")
  }
  for (n in xml_find_all(doc, "//select[@id]")) {
    opts <- xml_find_all(n, ".//option")
    if (length(opts) == 0) next
    sel <- opts[!is.na(xml_attr(opts, "selected"))]
    out[[xml_attr(n, "id")]] <- xml_attr(if (length(sel)) sel[[1]] else opts[[1]], "value")
  }
  for (n in xml_find_all(doc, "//textarea[@id]")) out[[xml_attr(n, "id")]] <- xml_text(n)
  out
}

ez_file_input <- function(path) {
  data.frame(name = basename(path), size = file.size(path), type = "",
             datapath = normalizePath(path), stringsAsFactors = FALSE)
}

# ---- main entry -------------------------------------------------------------
# Merge recovered settings (see patches.yaml) into a temp copy of the config.
ez_patch_yaml <- function(yaml_path, patch, dest) {
  y <- yaml::read_yaml(yaml_path)
  for (k in names(patch$heatmaps)) {
    i <- as.integer(k)
    stopifnot(i <= length(y[["visual definitions"]]$heatmaps))
    for (f in names(patch$heatmaps[[k]])) y[["visual definitions"]]$heatmaps[[i]][[f]] <- patch$heatmaps[[k]][[f]]
  }
  yaml::write_yaml(y, dest)
  dest
}

ez_render <- function(yaml, tree, csv, out_file, rdata = NULL, annot = NULL, patch = NULL,
                      class_column = "Cell.type.2", class_title = "Cell type",
                      seed = 1) {
  stopifnot(file.exists(yaml), file.exists(tree), file.exists(csv), is.null(rdata) || file.exists(rdata),
            is.null(annot) || file.exists(annot))
  tree <- normalizePath(tree); csv <- normalizePath(csv)
  if (!is.null(rdata)) rdata <- normalizePath(rdata)
  if (!is.null(annot)) annot <- normalizePath(annot)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(normalizePath(dirname(out_file)), basename(out_file))  # before setwd()
  if (!exists("server", envir = globalenv()) || !exists("PAPER_PALETTE_NORM", envir = globalenv())) ez_load_app()
  defaults <- ez_ui_defaults(get("ui", envir = globalenv()))
  result <- list(out_file = NA_character_, unmatched_classes = character(), ok = FALSE)

  # The engine also ggsave()s its own copy to "./"; keep that out of the app dir.
  work <- tempfile("ez_work_"); dir.create(work)
  old <- setwd(work); on.exit({ setwd(old); unlink(work, recursive = TRUE) }, add = TRUE)
  yaml <- normalizePath(yaml)
  if (!is.null(patch)) yaml <- ez_patch_yaml(yaml, patch, file.path(work, basename(yaml)))
  set.seed(seed)  # fisher.test(simulate.p.value = TRUE) inside the engine

  testServer(get("server", envir = globalenv()), {
    # Route update*Input() messages back into `input`, as the browser would.
    pending <- list()
    session$sendInputMessage <- function(inputId, message) {
      # individual_value: with "use all data" the browser's selectize would pick
      # the first individual in the CSV and rename the plot; keep the YAML's name.
      if (!is.null(message$value) && inputId != "individual_value") pending[[inputId]] <<- message$value
    }
    # Surface showNotification() text (the app reports early exits this way).
    session$sendNotification <- function(type, message) {
      if (identical(type, "show")) message("[ez-notify] ", gsub("<[^>]+>", "", paste(unlist(message$html), collapse = " ")))
    }
    pump <- function(max_rounds = 20) {
      for (i in seq_len(max_rounds)) {
        if (length(pending) == 0) return(invisible())
        batch <- pending; pending <<- list()
        batch <- lapply(batch, function(v) {
          if (is.character(v) && length(v) == 1) { n <- suppressWarnings(as.numeric(v)); if (!is.na(n) && grepl("^-?[0-9.eE+-]+$", v)) return(n) }
          v
        })
        do.call(session$setInputs, batch)
      }
    }
    set <- function(...) { session$setInputs(...); pump() }

    do.call(session$setInputs, defaults); pump()

    set(tree_file = ez_file_input(tree))
    set(csv_file  = ez_file_input(csv))
    if (!is.null(rdata)) set(rdata_file = ez_file_input(rdata))
    if (!is.null(annot)) set(annot_rdata_file = ez_file_input(annot))  # chromosome mapping (Annot)
    set(yaml_config = ez_file_input(yaml))

    # The YAML import copies heatmap_configs into the render list
    # values$heatmaps but drops the cnv_chr_* fields; in the browser they come
    # back when the heatmap tab is applied from its widgets. Copy them across.
    rd_cfg <- Filter(function(h) identical(h$data_source, "rdata"), values$heatmap_configs)
    if (length(rd_cfg)) {
      hm <- values$heatmaps
      rd_idx <- which(vapply(hm, function(h) identical(h$data_source, "rdata"), logical(1)))
      stopifnot(length(rd_idx) == length(rd_cfg))
      for (k in seq_along(rd_idx)) {
        chr_fields <- grep("^cnv_chr_", names(rd_cfg[[k]]), value = TRUE)
        hm[[rd_idx[k]]][chr_fields] <- rd_cfg[[k]][chr_fields]
      }
      values$heatmaps <- hm
    }

    # Process Data & Match IDs (tree tips vs CSV rows); use all rows, the
    # tree only contains one individual.
    set(use_all_data = TRUE)
    # (the process_data button observer does not fire under MockShinySession;
    # it only calls match_tree_with_csv(), so call that directly)
    match_tree_with_csv(); pump()
    message("[ez] after process_data: temp_csv=", format(values$temp_csv_path),
            " matched=", format(values$id_match$summary$matched), "/", format(values$id_match$summary$total_tree_labels),
            " filtered_rows=", format(nrow(values$filtered_csv)), " id_col=", format(input$id_column),
            " use_all=", format(input$use_all_data), " process_data=", format(input$process_data))

    # Discrete heatmaps with "custom colours": the saved YAML holds only the
    # colours the user changed (app S2.7 skips saving palette defaults). The
    # browser's colour pickers fill the rest from the palette by sorted value
    # position before Apply; do the same so unlisted values are not left NA.
    na_like <- c("#N/A", "#N/A!", "N/A", "#NA", "#VALUE!", "#REF!", "#DIV/0!", "#NULL!", "#NAME?", "#NUM!")
    hm <- values$heatmaps
    for (k in seq_along(hm)) {
      h <- hm[[k]]
      if (!isTRUE(h$is_discrete) || !isTRUE(h$man_define_colors) || identical(h$data_source, "rdata")) next
      cols_k <- intersect(unlist(h$columns), names(values$filtered_csv))
      vals <- unlist(lapply(cols_k, function(cn) {
        v <- as.character(values$filtered_csv[[cn]]); v[v %in% na_like | v == ""] <- NA
        nn <- v[!is.na(v)]; isnum <- !is.na(suppressWarnings(as.numeric(nn)))
        if (length(nn) && sum(isnum) > 0 && sum(isnum) / length(nn) < 0.5) v[!is.na(v) & !is.na(suppressWarnings(as.numeric(v)))] <- NA
        v
      }))
      uv_k <- sort(unique(na.omit(vals)))
      if (!length(uv_k) || length(uv_k) > 30) next
      pal <- if (!is.null(h$color_scheme)) h$color_scheme else "Set1"
      mx <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
      defaults_k <- if (length(uv_k) <= mx) RColorBrewer::brewer.pal(max(3, length(uv_k)), pal)[seq_along(uv_k)]
                    else colorRampPalette(RColorBrewer::brewer.pal(mx, pal))(length(uv_k))
      full <- setNames(as.list(defaults_k), uv_k)
      stored <- h$custom_colors
      for (nm in intersect(names(stored), uv_k)) full[[nm]] <- stored[[nm]]
      filled <- setdiff(uv_k, names(stored))
      if (length(filled)) message("[ez] heatmap '", h$title, "': palette default for ", paste(filled, collapse = ", "))
      hm[[k]]$custom_colors <- full
    }
    values$heatmaps <- hm

    # Classification: column coloured with the paper palette (what the
    # Classification tab's "Paper palette" radio button fills in).
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) values$filtered_csv else values$csv_data
    uv <- unique(csv_to_use[[class_column]]); uv <- uv[!is.na(uv)]
    cols <- unname(PAPER_PALETTE_NORM[func.norm.cat(uv)])
    result$unmatched_classes <<- as.character(uv[is.na(cols)])
    cols[is.na(cols)] <- rainbow(length(uv))[is.na(cols)]  # app default for a manual colour
    set(classification_column = class_column, classification_title = class_title,
        no_cluster_color = unname(PAPER_PALETTE_NORM[["na"]]))
    # Same object observeEvent(input$add_classification) builds (button
    # observers do not fire under MockShinySession, so store it directly).
    values$classifications <- list(list(
      title = input$classification_title,
      column = class_column,
      classes = lapply(seq_along(uv), function(i) list(column = class_column, value = uv[i],
                                                      display_name = as.character(uv[i]), color = cols[i])),
      fdr = input$fdr_perc,
      no_cluster_title = "No cluster",
      no_cluster_color = input$no_cluster_color,
      highlight = list(enabled = FALSE)
    ))
    values$active_classification_index <- 1
    values$temp_classification_preview <- NULL

    # Generate (bypass the 500 ms rapid-call guard) and download.
    Sys.sleep(0.6)
    update_yaml()
    generate_plot()
    if (is.null(values$current_plot)) stop("generate_plot() produced no plot")

    fmt <- tolower(tools::file_ext(out_file))
    set(output_format = fmt)
    tmp <- output$download_plot
    file.copy(tmp, out_file, overwrite = TRUE)
    result$out_file <<- out_file
    result$ok <<- file.exists(out_file)
    result$tip_order <<- values$rendered_tip_order
  })
  result
}
