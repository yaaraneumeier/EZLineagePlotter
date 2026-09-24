# ============================================================================
# MULTI-TREE LINEAGE PLOTTER - Standalone App
# Sources shared functions from the single-mode file, then defines its own
# UI and server. The single-mode file (EZlineagePlotter56.R) is NOT modified.
# ============================================================================

# Source only the function definitions (lines 1-9168) from the single-mode file.
# This gets us all library() calls and func.* definitions without triggering
# the single-mode UI/server/shinyApp.
main_file <- file.path(getwd(), "EZlineagePlotter56.R")
main_lines <- readLines(main_file)
ui_start <- grep("^ui <- dashboardPage", main_lines)[1]
eval(parse(text = main_lines[1:(ui_start - 1)]))

# Source the multi-tree module (UI builders + server installer)
source(file.path(getwd(), "EZlineagePlotter56_mt.R"), local = TRUE)

# --- UI ---
ui <- dashboardPage(
  dashboardHeader(title = "Multi-Tree Lineage Plotter"),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Upload Data", tabName = "mt_upload", icon = icon("upload")),
      menuItem("Tree Display", tabName = "mt_tree_display", icon = icon("tree")),
      menuItem("Classification", tabName = "mt_classification", icon = icon("palette")),
      menuItem("Bootstrap Values", tabName = "mt_bootstrap", icon = icon("percentage")),
      menuItem("Highlighting", tabName = "mt_highlighting", icon = icon("highlighter")),
      menuItem("Legend", tabName = "mt_legend", icon = icon("list")),
      menuItem("Extra", tabName = "mt_extra", icon = icon("plus-circle")),
      menuItem("Download", tabName = "mt_download", icon = icon("download")),
      menuItem("Configuration", tabName = "mt_config", icon = icon("cogs"))
    )
  ),

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

# --- Server ---
server <- function(input, output, session) {
  mt_install_server(input, output, session)
}

shiny::devmode(TRUE)
shinyApp(ui = ui, server = server)
