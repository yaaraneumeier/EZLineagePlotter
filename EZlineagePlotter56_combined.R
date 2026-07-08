# ============================================================================
# EZLineagePlotter - Combined Launcher
# Loads the single-mode app and multi-tree module, switches via URL parameter.
# EZlineagePlotter56.R is NOT modified - run it directly for guaranteed single mode.
#
# Usage:
#   source("EZlineagePlotter56_combined.R")
#   - Default (no parameter): Single Tree mode
#   - Add ?mode=multi to URL: Multiple Trees mode
#   - A floating button in the bottom-right corner switches between modes
# ============================================================================

# Load everything from the single-mode file except the shinyApp() launch.
# This gives us all library() calls, func.* definitions, plus `ui` and `server`.
main_lines <- readLines(file.path(getwd(), "EZlineagePlotter56.R"), warn = FALSE)
main_lines[grep("^shinyApp\\(", main_lines)] <- ""
eval(parse(text = main_lines))

single_ui <- ui
single_server <- server
rm(ui, server)

# Load multi-tree module (function definitions only)
source(file.path(getwd(), "EZlineagePlotter56_mt.R"), local = TRUE)

# --- Multi-tree UI (standalone dashboard) ---
multi_ui <- dashboardPage(
  dashboardHeader(title = "Multi-Tree Lineage Plotter"),
  dashboardSidebar(width = 300, sidebarMenu(
    menuItem("Upload Data", tabName = "mt_upload", icon = icon("upload")),
    menuItem("Tree Display", tabName = "mt_tree_display", icon = icon("tree")),
    menuItem("Classification", tabName = "mt_classification", icon = icon("palette")),
    menuItem("Bootstrap Values", tabName = "mt_bootstrap", icon = icon("percentage")),
    menuItem("Highlighting", tabName = "mt_highlighting", icon = icon("highlighter")),
    menuItem("Legend", tabName = "mt_legend", icon = icon("list")),
    menuItem("Extra", tabName = "mt_extra", icon = icon("plus-circle")),
    menuItem("Download", tabName = "mt_download", icon = icon("download")),
    menuItem("Configuration", tabName = "mt_config", icon = icon("cogs"))
  )),
  dashboardBody(
    shinyjs::useShinyjs(),
    tabItems(
      mt_tabItem_upload(),
      mt_tabItem_tree_display(),
      mt_tabItem_classification(),
      mt_tabItem_bootstrap(),
      mt_tabItem_highlighting(),
      mt_tabItem_legend(),
      mt_tabItem_extra(),
      mt_tabItem_download(),
      mt_tabItem_config()
    )
  )
)

# Mode indicator + switch button (injected via JS into sidebar)
mode_switch_btn <- function(current_mode) {
  if (current_mode == "single") {
    label <- "Switch to Multiple Trees ▶"
    url <- "?mode=multi"
    badge <- "Single Tree Mode"
  } else {
    label <- "◀ Switch to Single Tree"
    url <- "?mode=single"
    badge <- "Multiple Trees Mode"
  }
  tags$script(HTML(sprintf("
    $(document).ready(function() {
      $('.sidebar').append(
        '<div style=\"padding:10px 15px;margin-top:20px;border-top:1px solid #4b646f;\">' +
        '<div style=\"text-align:center;margin-bottom:8px;padding:4px 8px;background:#1a2226;' +
        'border-radius:3px;color:#b8c7ce;font-size:11px;font-weight:bold;\">%s</div>' +
        '<a href=\"%s\" style=\"display:block;padding:8px 12px;background:#3c8dbc;color:white;' +
        'border-radius:4px;text-decoration:none;font-size:12px;font-weight:bold;text-align:center;\">' +
        '%s</a></div>'
      );
    });
  ", badge, url, label)))
}

# --- Dynamic UI: serves the right dashboard based on ?mode= ---
ui <- function(request) {
  query <- parseQueryString(request$QUERY_STRING)
  mode <- if (identical(query$mode, "multi")) "multi" else "single"

  if (mode == "multi") {
    htmltools::tagAppendChildren(multi_ui, mode_switch_btn("multi"))
  } else {
    htmltools::tagAppendChildren(single_ui, mode_switch_btn("single"))
  }
}

# --- Combined server: delegates to the right mode ---
server <- function(input, output, session) {
  query <- isolate(parseQueryString(session$clientData$url_search))
  mode <- if (identical(query$mode, "multi")) "multi" else "single"

  if (mode == "multi") {
    mt_install_server(input, output, session)
  } else {
    single_server(input, output, session)
  }
}

# Dev mode is intentionally OFF for end users: it serves un-versioned, full
# Shiny JS and enables file-autoreload, which leaves the browser holding stale
# assets after a restart (requiring a hard refresh). It also overrode the
# container's shiny.devmode=FALSE. Keep it off here.
shiny::devmode(FALSE)
shinyApp(ui = ui, server = server)
