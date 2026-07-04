# ================================================================
# COMPARISON MODE MODULE
# Compares TWO trees as a tanglegram (facing trees + tip-connecting
# lines), with an entanglement score and (later) branch untangling.
# Sourced from EZlineagePlotter56_combined.R — does NOT touch single-
# or multi-mode code. All Shiny input IDs are prefixed with "cmp_".
#
# Milestones done: 1 (skeleton), 2 (upload + tip match + prune +
#   gray "Original" tanglegram + baseline crossing score).
# Engine helpers are namespaced "cmp." so they never clobber the
# single-mode func.* definitions. See docs/comparison_mode_plan.md.
# ================================================================

# ---------------------------------------------------------------
# Engine helpers (namespaced cmp.*)
# ---------------------------------------------------------------

# Read one tree; if a file holds several, take the first.
cmp.read.tree <- function(path) {
  t <- ape::read.tree(path)
  if (inherits(t, "multiPhylo")) t <- t[[1]]
  t
}

# Map each tip label to a canonical CSV ID using the app's general matcher
# (exact / numeric / prefix-suffix). Unmatched tips -> NA.
cmp.canonical.ids <- function(tip_labels, csv_ids) {
  tip_labels <- as.character(tip_labels)
  res <- match_tree_ids_with_csv(tip_labels, csv_ids)   # single-mode global helper
  m <- res$mapping
  out <- vapply(tip_labels, function(l) {
    v <- m[[l]]
    if (is.null(v) || length(v) == 0) NA_character_ else as.character(v[1])
  }, character(1))
  names(out) <- tip_labels
  out
}

# Relabel a tree's tips to their canonical CSV ID (unmatched tips keep their
# original label, so they still draw but won't connect to anything).
cmp.relabel.tree <- function(tree, canon) {
  new <- unname(canon[tree$tip.label])
  na_idx <- is.na(new)
  new[na_idx] <- tree$tip.label[na_idx]
  tree$tip.label <- new
  tree
}

# Keep only tips whose (relabeled) id is in shared_ids.
cmp.prune.to.shared <- function(tree, shared_ids) {
  drop <- tree$tip.label[!(tree$tip.label %in% shared_ids)]
  if (length(drop) > 0 && length(drop) < length(tree$tip.label)) {
    tree <- ape::drop.tip(tree, drop)
  }
  tree
}

# Do the two connecting segments (A: x=-1..1 at yA1..yB1) and (B: yA2..yB2) cross?
cmp.check.lines.intersect <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
  d1 <- c(x2 - x1, y2 - y1)
  d2 <- c(x4 - x3, y4 - y3)
  det <- d1[1] * d2[2] - d1[2] * d2[1]
  if (det == 0) return(FALSE)
  A <- matrix(c(d1[1], -d2[1], d1[2], -d2[2]), nrow = 2)
  b <- c(x3 - x1, y3 - y1)
  sol <- tryCatch(solve(A, b), error = function(e) NULL)
  if (is.null(sol)) return(FALSE)
  t <- sol[1]; s <- sol[2]
  (t >= 0 && t <= 1 && s >= 0 && s <= 1)
}

# Entanglement score: number of crossing tip-connecting lines, over the tips
# shared by both trees (by canonical id). Uses tip y-positions from the two
# ggtree data frames. Counts each unordered pair once.
cmp.make.X.score <- function(d1, d2) {
  tA <- d1[d1$isTip == TRUE & !is.na(d1$isTip), c("label", "y")]
  tB <- d2[d2$isTip == TRUE & !is.na(d2$isTip), c("label", "y")]
  shared <- intersect(as.character(tA$label), as.character(tB$label))
  n <- length(shared)
  if (n < 2) return(0L)
  yA <- setNames(tA$y, as.character(tA$label))[shared]
  yB <- setNames(tB$y, as.character(tB$label))[shared]
  cross <- 0L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (cmp.check.lines.intersect(-1, yA[i], 1, yB[i], -1, yA[j], 1, yB[j])) {
        cross <- cross + 1L
      }
    }
  }
  cross
}

# Compute the rectangular tree layout via ggtree (only to get x/y coordinates),
# returned as a plain data frame. ggtree may emit harmless fortify() warnings
# with newer ggplot2 — suppress them; the coordinates are still valid.
cmp.layout <- function(tree) {
  d <- suppressWarnings(ggtree::ggtree(tree)$data)
  as.data.frame(d)[, c("parent", "node", "x", "y", "isTip", "label")]
}

# Reconstruct the rectangular ("elbow") edge segments of a tree from its layout:
# a horizontal segment at the child's y from the parent's x to the child's x,
# and a vertical segment at the parent's x from the parent's y to the child's y.
cmp.edges <- function(d) {
  pc <- d[, c("node", "x", "y")]
  names(pc) <- c("node", "px", "py")
  m <- merge(d, pc, by.x = "parent", by.y = "node", all.x = TRUE)
  m <- m[!is.na(m$px) & m$node != m$parent, ]   # drop the root (parent == node)
  h <- data.frame(x = m$px, xend = m$x,  y = m$y,  yend = m$y)   # horizontal
  v <- data.frame(x = m$px, xend = m$px, y = m$py, yend = m$y)   # vertical
  rbind(h, v)
}

# Build the tanglegram ggplot: tree A on the left, tree B mirrored on the right,
# grey lines connecting tips that share a canonical id. Drawn entirely with
# geom_segment (version-robust; no geom_tree layer). Gray for now (classification
# coloring arrives later). Returns list($plot, $score, $shared).
cmp.draw.tanglegram <- function(treeA, treeB) {
  dA <- cmp.layout(treeA)
  dB <- cmp.layout(treeB)

  # Depth-normalize x so both trees reach the same maximum depth.
  mxA <- max(dA$x, na.rm = TRUE); mxB <- max(dB$x, na.rm = TRUE)
  desired <- max(mxA, mxB)
  if (is.finite(desired) && mxA > 0) dA$x <- dA$x * desired / mxA
  if (is.finite(desired) && mxB > 0) dB$x <- dB$x * desired / mxB

  # Put both trees on a common vertical span so they line up.
  yr <- range(c(dA$y, dB$y), na.rm = TRUE)
  if (diff(yr) > 0) {
    dA$y <- scales::rescale(dA$y, to = yr)
    dB$y <- scales::rescale(dB$y, to = yr)
  }

  score <- cmp.make.X.score(dA, dB)   # computed before mirroring (uses y only)

  # Mirror tree B and push it to the right of tree A (with a small gap).
  gap <- 0.15 * max(dA$x, na.rm = TRUE)
  if (!is.finite(gap) || gap <= 0) gap <- 1
  dB$x <- max(dB$x, na.rm = TRUE) - dB$x + max(dA$x, na.rm = TRUE) + gap

  eA <- cmp.edges(dA)
  eB <- cmp.edges(dB)

  tA <- dA[dA$isTip %in% TRUE, c("x", "y", "label")]
  tB <- dB[dB$isTip %in% TRUE, c("x", "y", "label")]
  shared_labels <- intersect(as.character(tA$label), as.character(tB$label))
  conn <- NULL
  if (length(shared_labels) > 0) {
    conn <- do.call(rbind, lapply(shared_labels, function(lb) {
      a <- tA[as.character(tA$label) == lb, ][1, ]
      b <- tB[as.character(tB$label) == lb, ][1, ]
      data.frame(x = a$x, xend = b$x, y = a$y, yend = b$y)
    }))
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = eA, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                          color = "grey20", linewidth = 0.5) +
    ggplot2::geom_segment(data = eB, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                          color = "grey20", linewidth = 0.5)
  if (!is.null(conn) && nrow(conn) > 0) {
    p <- p + ggplot2::geom_segment(data = conn, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   color = "grey60", linewidth = 0.35)
  }
  p <- p + ggplot2::theme_void()

  list(plot = p, score = score, shared = length(shared_labels))
}

# ---------------------------------------------------------------
# UI: one tabItem per sidebar menuItem
# ---------------------------------------------------------------

cmp_placeholder <- function(title, note) {
  box(title = title, status = "primary", solidHeader = TRUE, width = 12,
      tags$p(style = "color:#555;", note),
      tags$p(style = "color:#999;font-size:12px;",
             "This tab's functionality is being built in stages — see the roadmap on the Upload tab."))
}

cmp_tabItem_upload <- function() {
  tabItem(
    tabName = "cmp_upload",
    fluidRow(
      box(title = "Comparison Mode", status = "success", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          tags$div(style = "background:#d4edda;padding:15px;border-radius:5px;border:2px solid #155724;",
            tags$h4(style = "color:#155724;margin:0;", "Compare two trees (tanglegram)"),
            tags$p(style = "margin:10px 0 0 0;color:#155724;",
              "Draw two trees facing each other with lines connecting shared tips, score how ",
              "'entangled' they are, and (soon) rotate branches to untangle them."))
      )
    ),
    fluidRow(
      box(title = "1. Upload", status = "primary", solidHeader = TRUE, width = 5,
          fileInput("cmp_treeA", "Tree A (Newick):",
                    accept = c(".newick", ".nwk", ".tree", ".txt")),
          fileInput("cmp_treeB", "Tree B (Newick):",
                    accept = c(".newick", ".nwk", ".tree", ".txt")),
          fileInput("cmp_csv", "Classification CSV (optional, for coloring):", accept = c(".csv")),
          selectInput("cmp_id_column", "ID column (for classification):", choices = NULL),
          actionButton("cmp_check_match", "Check tip matching", class = "btn-info"),
          tags$p(style = "color:#999;font-size:12px;margin-top:8px;",
                 "Tip matching compares the two trees' labels directly with the app's general matcher ",
                 "(exact / numeric / prefix) — no CSV needed. The CSV/ID column is only used later for ",
                 "classification coloring.")
      ),
      box(title = "2. Tip matching", status = "primary", solidHeader = TRUE, width = 7,
          uiOutput("cmp_match_report"),
          tags$hr(),
          radioButtons("cmp_prune", "If the two trees don't fully overlap:",
                       choices = c("Keep all tips" = "keep",
                                   "Use only shared tips" = "shared"),
                       selected = "keep"),
          actionButton("cmp_build", "Build comparison", class = "btn-success")
      )
    )
  )
}

cmp_tabItem_classification <- function() {
  tabItem(tabName = "cmp_classification",
    fluidRow(cmp_placeholder("Classification",
      "Pick one CSV column and palette, applied to BOTH trees so the comparison is meaningful. Trees start gray."))
  )
}

cmp_tabItem_compare <- function() {
  tabItem(tabName = "cmp_compare",
    fluidRow(
      box(title = "Comparison preview", status = "primary", solidHeader = TRUE, width = 12,
          textOutput("cmp_score"),
          tags$br(),
          plotOutput("cmp_preview", height = "650px"))
    ),
    fluidRow(cmp_placeholder("Untangle (coming next)",
      "Untangle controls, the versions list, and per-version scores will live here in the next milestone."))
  )
}

cmp_tabItem_download <- function() {
  tabItem(tabName = "cmp_download",
    fluidRow(cmp_placeholder("Download",
      "Export the current version as PDF/PNG, and download the rotated tree(s) as Newick."))
  )
}

cmp_tabItem_config <- function() {
  tabItem(tabName = "cmp_config",
    fluidRow(
      box(title = "Configuration", status = "primary", solidHeader = TRUE, width = 12,
          tags$p("Comparison YAML configuration (save/load) will appear here."),
          verbatimTextOutput("cmp_yaml_output"))
    )
  )
}

# ---------------------------------------------------------------
# Server
# ---------------------------------------------------------------
cmp_install_server <- function(input, output, session) {

  cmp_values <- reactiveValues(
    treeA = NULL, treeB = NULL,          # phylo objects (raw)
    csv_data = NULL,
    match = NULL,                        # list(canonA, canonB, shared, aOnly, bOnly, summary)
    fig = NULL,                          # ggplot of current comparison
    score = NULL, shared = NULL
  )

  # --- File uploads ---
  observeEvent(input$cmp_treeA, {
    cmp_values$treeA <- tryCatch(cmp.read.tree(input$cmp_treeA$datapath),
                                 error = function(e) { showNotification(paste("Tree A:", e$message), type = "error"); NULL })
  })
  observeEvent(input$cmp_treeB, {
    cmp_values$treeB <- tryCatch(cmp.read.tree(input$cmp_treeB$datapath),
                                 error = function(e) { showNotification(paste("Tree B:", e$message), type = "error"); NULL })
  })
  observeEvent(input$cmp_csv, {
    df <- tryCatch(data.table::fread(input$cmp_csv$datapath, data.table = FALSE),
                   error = function(e) { showNotification(paste("CSV:", e$message), type = "error"); NULL })
    cmp_values$csv_data <- df
    if (!is.null(df)) updateSelectInput(session, "cmp_id_column", choices = names(df), selected = names(df)[1])
  })

  # --- Tip matching: match Tree A's tips DIRECTLY against Tree B's tips with the
  # app's general matcher (exact / numeric / prefix). The CSV is NOT needed here;
  # it is only used for classification coloring (later milestone). ---
  observeEvent(input$cmp_check_match, {
    if (is.null(cmp_values$treeA) || is.null(cmp_values$treeB)) {
      showNotification("Upload both trees first.", type = "error"); return()
    }
    labA <- as.character(cmp_values$treeA$tip.label)
    labB <- as.character(cmp_values$treeB$tip.label)

    # Map each B tip label to an A tip label (or NA if it has no counterpart).
    res <- match_tree_ids_with_csv(labB, labA)   # single-mode global matcher
    mapB <- res$mapping
    canonB <- vapply(labB, function(l) {
      v <- mapB[[l]]
      if (is.null(v) || length(v) == 0) NA_character_ else as.character(v[1])
    }, character(1))
    names(canonB) <- labB

    matched_A <- unique(stats::na.omit(canonB))  # A labels that have a B counterpart
    shared <- intersect(labA, matched_A)
    cmp_values$match <- list(
      canonA = setNames(labA, labA),             # A maps to itself
      canonB = canonB,                           # B label -> matched A label (NA = unmatched)
      shared = shared,
      aOnly = setdiff(labA, shared),             # A tips with no B counterpart
      bOnly = labB[is.na(canonB)],               # B tips with no A counterpart
      nA = length(labA), nB = length(labB)
    )
  })

  output$cmp_match_report <- renderUI({
    m <- cmp_values$match
    if (is.null(m)) return(tags$p(style = "color:#888;", "Upload two trees, then click “Check tip matching.”"))
    fmt <- function(v, k = 10) if (length(v) == 0) "(none)" else paste(c(utils::head(v, k), if (length(v) > k) paste0("… (+", length(v) - k, " more)")), collapse = ", ")
    tagList(
      tags$p(sprintf("Tree A: %d tips.", m$nA)),
      tags$p(sprintf("Tree B: %d tips.", m$nB)),
      tags$p(tags$b(sprintf("Shared tips: %d", length(m$shared)))),
      if (length(m$aOnly) > 0) tags$p(style = "color:#a00;", sprintf("Only in A (%d): %s", length(m$aOnly), fmt(m$aOnly))),
      if (length(m$bOnly) > 0) tags$p(style = "color:#a00;", sprintf("Only in B (%d): %s", length(m$bOnly), fmt(m$bOnly))),
      if (length(m$shared) == 0) tags$p(style = "color:#a00;font-weight:bold;", "No shared tips — the two trees' tip labels don't match (check label formats).")
    )
  })

  # --- Build the (gray) Original comparison ---
  observeEvent(input$cmp_build, {
    m <- cmp_values$match
    if (is.null(m)) { showNotification("Run “Check tip matching” first.", type = "error"); return() }
    if (length(m$shared) == 0) { showNotification("No shared tips to compare.", type = "error"); return() }

    res <- tryCatch({
      tA <- cmp.relabel.tree(cmp_values$treeA, m$canonA)
      tB <- cmp.relabel.tree(cmp_values$treeB, m$canonB)
      if (identical(input$cmp_prune, "shared")) {
        tA <- cmp.prune.to.shared(tA, m$shared)
        tB <- cmp.prune.to.shared(tB, m$shared)
      }
      cmp.draw.tanglegram(tA, tB)
    }, error = function(e) { showNotification(paste("Build failed:", e$message), type = "error"); NULL })

    if (!is.null(res)) {
      cmp_values$fig <- res$plot
      cmp_values$score <- res$score
      cmp_values$shared <- res$shared
      showNotification("Comparison built — see the Compare / Untangle tab.", type = "message")
    }
  })

  output$cmp_preview <- renderPlot({
    if (is.null(cmp_values$fig)) {
      plot.new()
      text(0.5, 0.6, "Comparison mode", cex = 1.6, font = 2)
      text(0.5, 0.42, "Upload two trees + CSV, check matching, then Build comparison.", cex = 1.05, col = "#666666")
    } else {
      print(cmp_values$fig)
    }
  })

  output$cmp_score <- renderText({
    if (is.null(cmp_values$score)) "No comparison built yet."
    else sprintf("Original: %d crossings over %d shared tips.", cmp_values$score, cmp_values$shared)
  })

  output$cmp_yaml_output <- renderText({ "No comparison configured yet." })
}
