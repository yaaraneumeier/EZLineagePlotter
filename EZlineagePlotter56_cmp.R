# ================================================================
# COMPARISON MODE MODULE  (Milestone 1: skeleton)
# Compares TWO trees as a tanglegram (facing trees + tip-connecting
# lines), with an entanglement score and optional branch untangling.
# Sourced from EZlineagePlotter56_combined.R — does NOT touch single-
# mode code. All Shiny input IDs are prefixed with "cmp_".
#
# The comparison ENGINE (ported from the original notebook:
# func.make.compare_fig / func.make.X.score / func_tree_rotation ...)
# will be added, namespaced "cmp.", in later milestones. It must stay
# namespaced so it never clobbers the single-mode func.* definitions.
#
# See docs/comparison_mode_plan.md for the full design.
# ================================================================

# Small helper: a placeholder panel for tabs whose functionality lands
# in a later milestone, so the mode loads and is clickable now.
cmp_placeholder <- function(title, note) {
  box(title = title, status = "primary", solidHeader = TRUE, width = 12,
      tags$p(style = "color:#555;", note),
      tags$p(style = "color:#999;font-size:12px;",
             "This tab's functionality is being built in stages — see the roadmap on the Upload tab."))
}

# --- Individual tabItem functions (one per sidebar menuItem) ---

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
              "'entangled' they are, and optionally rotate branches to untangle them."),
            tags$p(style = "margin:8px 0 0 0;color:#155724;", tags$strong("Build roadmap:")),
            tags$ul(style = "color:#155724;",
              tags$li("1. Mode skeleton + switch (this build)"),
              tags$li("2. Upload two trees + CSV, check tip matching, prune or keep"),
              tags$li("3. Classification: one shared column + palette for both trees"),
              tags$li("4. Untangle: rotate to reduce crossings (off by default)"),
              tags$li("5. Versions: keep and switch between untangle results"),
              tags$li("6. Download: plot (PDF/PNG) + rotated Newick"),
              tags$li("7. Save/Load comparison configuration (YAML)"))
          )
      )
    ),
    fluidRow(
      cmp_placeholder("Upload",
        "Here you will upload Tree A, Tree B and a classification CSV, then see a tip-match report and choose whether to prune non-matching tips.")
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
    fluidRow(cmp_placeholder("Compare / Untangle",
      "The tanglegram preview, untangle controls, the versions list, and the entanglement score will live here.")),
    fluidRow(
      box(title = "Preview", status = "primary", solidHeader = TRUE, width = 12,
          plotOutput("cmp_preview", height = "600px"))
    )
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

# --- Server ---
cmp_install_server <- function(input, output, session) {

  # Reactive state, isolated from single/multi modes (all prefixed cmp_).
  cmp_values <- reactiveValues(
    treeA = NULL, treeB = NULL,          # phylo objects
    nameA = NULL, nameB = NULL,          # original file names
    pathA = NULL, pathB = NULL,          # temp file paths
    csv_data = NULL, csv_path = NULL,
    shared_column = NULL, palette = NULL,
    prune_choice = "keep",               # "shared" | "keep"
    matched = NULL,                      # tip-match report
    versions = list(),                   # id -> version object
    current_version = NULL,              # id of displayed version
    baseline_score = NULL                # crossings of the Original
  )

  # Milestone 1: a simple placeholder preview so the mode renders cleanly.
  output$cmp_preview <- renderPlot({
    plot.new()
    text(0.5, 0.6, "Comparison mode", cex = 1.6, font = 2)
    text(0.5, 0.42, "Upload two trees to begin (coming in Milestone 2).", cex = 1.1, col = "#666666")
  })

  output$cmp_yaml_output <- renderText({
    "No comparison configured yet."
  })
}
