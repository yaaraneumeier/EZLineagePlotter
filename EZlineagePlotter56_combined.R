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

# Load comparison module (function definitions only). Sourced local = TRUE so its
# namespaced cmp.* engine (added in later milestones) can't clobber single-mode func.*
source(file.path(getwd(), "EZlineagePlotter56_cmp.R"), local = TRUE)

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

# --- Comparison UI (standalone dashboard) ---
compare_ui <- dashboardPage(
  dashboardHeader(title = "Compare Two Trees"),
  dashboardSidebar(width = 300, sidebarMenu(
    menuItem("Upload Data", tabName = "cmp_upload", icon = icon("upload")),
    menuItem("Classification", tabName = "cmp_classification", icon = icon("palette")),
    menuItem("Compare / Untangle", tabName = "cmp_compare", icon = icon("code-branch")),
    menuItem("Download", tabName = "cmp_download", icon = icon("download")),
    menuItem("Configuration", tabName = "cmp_config", icon = icon("cogs"))
  )),
  dashboardBody(
    shinyjs::useShinyjs(),
    tabItems(
      cmp_tabItem_upload(),
      cmp_tabItem_classification(),
      cmp_tabItem_compare(),
      cmp_tabItem_download(),
      cmp_tabItem_config()
    )
  )
)

# Mode indicator + switch buttons (injected via JS into sidebar).
# Three modes: single, multi, compare. Shows the current mode as a badge and a
# link to each of the other two modes.
mode_switch_btn <- function(current_mode) {
  modes <- list(
    single  = list(url = "?mode=single",  label = "Single Tree"),
    multi   = list(url = "?mode=multi",   label = "Multiple Trees"),
    compare = list(url = "?mode=compare", label = "Comparison")
  )
  badge <- paste0(modes[[current_mode]]$label, " Mode")
  others <- setdiff(names(modes), current_mode)
  links <- paste0(sapply(others, function(m) {
    sprintf(paste0("<a href=\"%s\" style=\"display:block;padding:8px 12px;margin-top:6px;",
                   "background:#3c8dbc;color:white;border-radius:4px;text-decoration:none;",
                   "font-size:12px;font-weight:bold;text-align:center;\">Switch to %s</a>"),
            modes[[m]]$url, modes[[m]]$label)
  }), collapse = "")
  tags$script(HTML(sprintf("
    $(document).ready(function() {
      $('.sidebar').append(
        '<div style=\"padding:10px 15px;margin-top:20px;border-top:1px solid #4b646f;\">' +
        '<div style=\"text-align:center;margin-bottom:4px;padding:4px 8px;background:#1a2226;' +
        'border-radius:3px;color:#b8c7ce;font-size:11px;font-weight:bold;\">%s</div>' +
        '%s</div>'
      );
    });
  ", badge, links)))
}

# --- Dynamic UI: serves the right dashboard based on ?mode= ---
resolve_mode <- function(query) {
  if (identical(query$mode, "multi")) "multi"
  else if (identical(query$mode, "compare")) "compare"
  else "single"
}

ui <- function(request) {
  mode <- resolve_mode(parseQueryString(request$QUERY_STRING))

  if (mode == "multi") {
    htmltools::tagAppendChildren(multi_ui, mode_switch_btn("multi"))
  } else if (mode == "compare") {
    htmltools::tagAppendChildren(compare_ui, mode_switch_btn("compare"))
  } else {
    htmltools::tagAppendChildren(single_ui, mode_switch_btn("single"))
  }
}

# --- Combined server: delegates to the right mode ---
server <- function(input, output, session) {
  mode <- resolve_mode(isolate(parseQueryString(session$clientData$url_search)))

  if (mode == "multi") {
    mt_install_server(input, output, session)
  } else if (mode == "compare") {
    cmp_install_server(input, output, session)
  } else {
    single_server(input, output, session)
  }
}

shiny::devmode(TRUE)
shinyApp(ui = ui, server = server)
