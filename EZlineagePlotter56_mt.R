# ================================================================
# MULTI-TREE MODE MODULE
# Plots multiple trees on one page with shared classifications.
# Sourced from EZlineagePlotter56.R — does NOT touch single-mode code.
# All Shiny input IDs prefixed with "mt_".
# ================================================================

# (Old mt_ui_*_tab and mt_menuItem functions removed — content is now
#  inline in the mt_tabItem_*() functions below)

# --- Individual tabItem functions (one per sidebar menuItem) ---

mt_tabItem_upload <- function() {
  tabItem(
    tabName = "mt_upload",

    # Row 1: Version bar (identical to single mode)
    fluidRow(
      box(
        title = "EZLineagePlotter - Stable Release",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        tags$div(style = "background: #d4edda; padding: 15px; border-radius: 5px; border: 2px solid #155724;",
                 tags$h4(style = "color: #155724; margin: 0;", "Version S3.13 Stable"),
                 tags$p(style = "margin: 10px 0 0 0; color: #155724;",
                        tags$strong("New in S3.13:"),
                        tags$ul(
                          tags$li("Multiple Trees Mode: plot several trees on one page with shared classification coloring"),
                          tags$li("Per-tree rotation, background color, and display titles"),
                          tags$li("YAML import support for multi-tree shared settings")
                        ),
                        tags$strong("From S3.12:"),
                        tags$ul(
                          tags$li("Fix crash when YAML references CSV columns that don't exist"),
                          tags$li("Filter stray numeric values from discrete heatmap legends"),
                          tags$li("Fix legend showing #N/A and other Excel NA-like strings as discrete values"),
                          tags$li("Fix discrete color corruption during column changes and UI rebuilds"),
                          tags$li("Fix CSV heatmaps losing data_source after Apply Heatmaps"),
                          tags$li("Fix discrete heatmap legend mislabeling when columns have different value sets")
                        ),
                        tags$strong("From S3.11:"),
                        tags$ul(
                          tags$li("SNP Analysis Tab: exploratory mutation analysis with VAF heatmaps and classification calls"),
                          tags$li("Chromosome boundary lines and labels for RData CNV heatmaps"),
                          tags$li("Separate chromosome mapping RData file upload"),
                          tags$li("Detailed display mode for RData CNV heatmaps (per-cell resolution)"),
                          tags$li("Rotated tree Newick download"),
                          tags$li("Manual RGB/Hex color input for heatmap colors"),
                          tags$li("Google Fonts support for legend text")
                        ),
                        tags$strong("Core Features:"),
                        tags$ul(
                          tags$li("RData CNV heatmaps from QDNAseq/scIMPACT pipelines"),
                          tags$li("Multiple heatmaps with discrete/continuous color scales"),
                          tags$li("Tree visualization with classification coloring"),
                          tags$li("Highlight regions with customizable ellipses"),
                          tags$li("Bootstrap value display and flexible legend styling"),
                          tags$li("Export to PDF/PNG with custom dimensions")
                        )
                 )
        )
      )
    ),

    # Row 2: Tree Files (left) + Classification Data (right)
    fluidRow(
      box(
        title = "Tree Files",
        status = "primary",
        solidHeader = TRUE,
        width = 6,
        tags$p(class = "text-muted", tags$small(
          "Upload one or more newick files. Click '+ Add another tree' to add a new upload slot.")),
        uiOutput("mt_tree_upload_slots"),
        actionButton("mt_add_upload_slot", "+ Add another tree",
                     class = "btn-success btn-sm", icon = icon("plus")),
        tags$br(), tags$br(),
        verbatimTextOutput("mt_tree_summary")
      ),

      box(
        title = "Classification Data",
        status = "primary",
        solidHeader = TRUE,
        width = 6,
        fileInput("mt_csv_file", "Upload CSV Classification File",
                  accept = c(".csv", ".tsv", ".txt")),
        selectizeInput("mt_id_column", "Select ID Column", choices = NULL,
                       options = list(placeholder = "Select...")),
        selectizeInput("mt_individual_column", "Select Individual Column",
                       choices = NULL,
                       options = list(placeholder = "Select...")),
        checkboxInput("mt_use_all_data",
                      "Use all data (ignore individual filtering)",
                      value = FALSE),
        verbatimTextOutput("mt_csv_summary"),
        br(),
        conditionalPanel(
          condition = "output.mt_files_loaded == 'TRUE'",
          actionButton("mt_process_data", "Process Data & Match IDs",
                       class = "btn-success btn-lg",
                       icon = icon("check-circle"),
                       style = "width: 100%;")
        )
      )
    ),

    # Row 3: Per-tree individual mapping (only when individual_column is selected)
    conditionalPanel(
      condition = "input.mt_individual_column != null && input.mt_individual_column != '' && input.mt_use_all_data == false",
      fluidRow(
        box(
          title = "Tree-to-Individual Mapping",
          status = "warning",
          solidHeader = TRUE,
          width = 12,
          tags$div(
            style = "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            tags$p(style = "margin: 0; color: #856404;",
                   icon("info-circle"),
                   " Each tree must be matched to an individual in the CSV.",
                   tags$br(),
                   tags$small("Select which individual value each uploaded tree corresponds to. Use 'Auto-match by filename' to match tree names against individual names automatically."))
          ),
          uiOutput("mt_per_tree_individual_ui"),
          tags$hr(),
          actionButton("mt_auto_match_individuals", "Auto-match by filename",
                       icon = icon("magic"), class = "btn-info btn-sm")
        )
      )
    ),

    # Row 4: Matching status
    fluidRow(
      box(
        title = "Matching Status",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        verbatimTextOutput("mt_log")
      )
    ),

    # Row 5: Optional YAML import
    fluidRow(
      box(
        title = "Import Saved Settings (Optional)",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        collapsed = TRUE,
        tags$div(
          style = "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
          tags$p(style = "margin: 0; color: #856404;",
                 icon("info-circle"),
                 " Load shared visual settings from a previously saved YAML configuration.",
                 tags$br(),
                 tags$small("This applies classification, bootstrap, highlighting, legend, and extra settings. Heatmap/SNP blocks are ignored."))
        ),
        fileInput("mt_yaml_config", "Choose YAML Configuration File",
                  accept = c(".yaml", ".yml"))
      )
    ),

    # Row 6: Preview
    fluidRow(
      box(
        title = "Multi-Tree Preview",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        actionButton("mt_update_preview", "Update Preview",
                     class = "btn-primary", icon = icon("refresh")),
        tags$hr(),
        imageOutput("mt_combined_plot", height = "auto")
      )
    )
  )
}

mt_tabItem_tree_display <- function() {
  tabItem(tabName = "mt_tree_display",
    fluidRow(
      box(
        title = "Tree Appearance",
        status = "primary",
        solidHeader = TRUE,
        width = 4,
        checkboxInput("mt_trim_tips", "Trim Tips", value = FALSE),
        conditionalPanel(
          condition = "input.mt_trim_tips == true",
          sliderInput("mt_tip_length", "Tip Length", min = 0.01, max = 0.2, value = 0.05, step = 0.01)
        ),
        sliderInput("mt_edge_width", "Edge Width Multiplier", min = 0.5, max = 3, value = 1, step = 0.1),
        sliderInput("mt_tip_font_size", "Tip Label Font Size", min = 1, max = 10, value = 3, step = 0.5),
        checkboxInput("mt_display_node_numbers", "Display Node Numbers", value = FALSE),
        sliderInput("mt_node_number_font_size", "Node Number Font Size", min = 1, max = 8, value = 3.5, step = 0.5),
        checkboxInput("mt_use_pvalues", "Use P-values for Branch Width", value = TRUE),
        conditionalPanel(
          condition = "input.mt_use_pvalues == true",
          sliderInput("mt_fdr_perc", "FDR Percentage", min = 0.01, max = 0.25, value = 0.1, step = 0.01)
        ),
        checkboxInput("mt_ladderize", "Ladderize Tree", value = FALSE),
        tags$hr(),
        actionButton("mt_apply_tree_display", "Apply & Preview",
                     icon = icon("eye"), class = "btn-primary btn-block")
      ),

      box(
        title = "Rotation Options",
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        selectInput("mt_tree_selector", "Select Tree:", choices = NULL),
        checkboxInput("mt_enable_rotation", "Enable Tree Rotation", value = FALSE),
        conditionalPanel(
          condition = "input.mt_enable_rotation == true",
          radioButtons("mt_rotation_type", "Rotation Type:",
                       choices = list(
                         "Primary Classification First" = "primary",
                         "Secondary Classification First" = "secondary",
                         "Manual Node Rotation" = "manual"
                       ), selected = "primary"),
          conditionalPanel(
            condition = "input.mt_rotation_type == 'manual'",
            selectizeInput("mt_nodes_to_rotate", "Select Nodes to Rotate",
                           choices = NULL, multiple = TRUE,
                           options = list(placeholder = "Select nodes to rotate")),
            checkboxInput("mt_highlight_selected_nodes", "Highlight Selected Nodes on Tree", value = FALSE),
            tags$div(style = "margin-left: 20px; margin-top: -10px; margin-bottom: 10px;",
                     tags$small(style = "color: #666;", "Shows red circles on selected nodes to help you verify your selection")
            ),
            br(),
            actionButton("mt_apply_manual_rotation", "Apply Rotation", icon = icon("rotate"), class = "btn-primary"),
            actionButton("mt_clear_manual_rotation", "Clear Manual Configuration", icon = icon("trash"), class = "btn-warning")
          ),
          uiOutput("mt_rotation_status_box"),
          uiOutput("mt_rotation_classes_ui")
        )
      )
    ),

    fluidRow(
      box(
        title = tagList(
          "Tree Preview ",
          span(id = "mt_status_waiting",
            style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
            icon("clock"), " Waiting for data"
          ),
          span(id = "mt_status_processing",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("spinner", class = "fa-spin"), " Processing..."
          ),
          span(id = "mt_status_ready",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("check-circle"), " Ready"
          )
        ),
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        imageOutput("mt_tree_display_preview", height = "auto")
      )
    )
  )
}

mt_tabItem_classification <- function() {
  tabItem(tabName = "mt_classification",
    fluidRow(
      box(
        title = "Classification Settings",
        status = "primary",
        solidHeader = TRUE,
        width = 4,
        selectizeInput("mt_classification_column", "Select Classification Column:",
                      choices = NULL,
                      options = list(placeholder = "Select column...")),
        textInput("mt_classification_title", "Legend Title:", value = "Cell type"),
        selectInput("mt_no_cluster_color", "No Cluster Color:",
                    choices = c("gray", "black", "white", "red"), selected = "gray"),
        tags$hr(),
        radioButtons("mt_classification_scope", "Apply to:",
                     choices = c("All trees" = "all", "Single tree" = "per_tree"),
                     selected = "all", inline = TRUE),
        conditionalPanel(
          condition = "input.mt_classification_scope == 'per_tree'",
          selectInput("mt_classification_tree_selector", "Which tree:",
                      choices = NULL),
          tags$p(class = "text-muted", style = "font-size: 12px;",
                 "This will override the shared classification for the selected tree only.")
        ),
        tags$hr(),
        actionButton("mt_update_classification_preview", "Apply & Preview",
                     icon = icon("eye"), class = "btn-primary btn-block"),
        conditionalPanel(
          condition = "input.mt_classification_scope == 'per_tree'",
          actionButton("mt_save_classification_no_render", "Save for This Tree (no render)",
                       icon = icon("save"), class = "btn-info btn-block btn-sm"),
          tags$p(class = "text-muted", style = "font-size: 11px; margin-top: 5px;",
                 "Tip: Save classifications for each tree first, then click Apply & Preview once to render all.")
        ),
        tags$br(),
        actionButton("mt_remove_classification", "Reset Classification",
                     icon = icon("undo"), class = "btn-warning btn-block btn-sm")
      ),
      box(
        title = "Classification Values",
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        uiOutput("mt_classification_values_ui")
      )
    ),
    fluidRow(
      box(
        title = tagList(
          "Preview ",
          span(id = "mt_class_status_waiting",
            style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
            icon("clock"), " Waiting for data"
          ),
          span(id = "mt_class_status_processing",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("spinner", class = "fa-spin"), " Processing..."
          ),
          span(id = "mt_class_status_ready",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("check-circle"), " Ready"
          )
        ),
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        imageOutput("mt_combined_plot_class", height = "auto")
      )
    )
  )
}

mt_tabItem_bootstrap <- function() {
  tabItem(tabName = "mt_bootstrap",
    fluidRow(
      box(
        title = "Bootstrap Display",
        status = "primary",
        solidHeader = TRUE,
        width = 4,
        checkboxInput("mt_enable_bootstrap", "Show Bootstrap Values", value = TRUE),
        conditionalPanel(
          condition = "input.mt_enable_bootstrap == true",
          radioButtons("mt_bootstrap_format", "Display Format:",
                       choices = list(
                         "Triangles" = "triangles",
                         "Raw Values" = "raw",
                         "Percentage" = "percentage",
                         "Color-coded Numbers" = "numbered_color",
                         "Color-coded Percentage" = "percentage_color"
                       ), selected = "triangles"),
          sliderInput("mt_bootstrap_param", "Bootstrap Precision (decimal places)",
                      min = 1, max = 5, value = 1, step = 1),
          sliderInput("mt_bootstrap_label_size", "Bootstrap Triangle Size:",
                      min = 0, max = 10, value = 1.5, step = 0.5, width = "100%"),
          sliderInput("mt_man_boot_x_offset", "Bootstrap Position (higher/lower):",
                      min = -0.5, max = 0.5, value = 0, step = 0.001, width = "100%"),
          tags$hr(),
          actionButton("mt_apply_bootstrap", "Apply & Preview",
                       icon = icon("eye"), class = "btn-primary btn-block")
        )
      ),

      box(
        title = tagList(
          "Preview ",
          span(id = "mt_boot_status_waiting",
            style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
            icon("clock"), " Waiting for data"
          ),
          span(id = "mt_boot_status_processing",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("spinner", class = "fa-spin"), " Processing..."
          ),
          span(id = "mt_boot_status_ready",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("check-circle"), " Ready"
          )
        ),
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        imageOutput("mt_bootstrap_preview", height = "auto")
      )
    )
  )
}

mt_tabItem_highlighting <- function() {
  tabItem(tabName = "mt_highlighting",
    fluidRow(
      box(
        title = "Highlight Settings",
        status = "primary",
        solidHeader = TRUE,
        width = 4,
        checkboxInput("mt_enable_highlight", "Enable Highlighting", value = FALSE),
        conditionalPanel(
          condition = "input.mt_enable_highlight == true",
          selectInput("mt_highlight_column", "Select Column for Highlighting",
                      choices = NULL, selected = NULL),
          selectizeInput("mt_highlight_values", "Select Values to Highlight",
                         choices = NULL, multiple = TRUE,
                         options = list(maxOptions = 1000,
                                        placeholder = "Select values to highlight")),
          textInput("mt_highlight_title", "Highlight Legend Title", value = "Highlight"),
          hr(),
          h5("Global Positioning:"),
          fluidRow(
            column(12,
                   tags$label("Vertical Offset (up/down)"),
                   fluidRow(
                     column(4, actionButton("mt_offset_down_big", "-1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")),
                     column(4, actionButton("mt_offset_down_small", "-0.1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")),
                     column(4, actionButton("mt_offset_down_tiny", "-0.01", class = "btn-sm btn-outline-secondary", style = "width: 100%;"))
                   ),
                   numericInput("mt_highlight_offset", NULL,
                                value = 0, min = -10, max = 10, step = 0.01, width = "100%"),
                   fluidRow(
                     column(4, actionButton("mt_offset_up_big", "+1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")),
                     column(4, actionButton("mt_offset_up_small", "+0.1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")),
                     column(4, actionButton("mt_offset_up_tiny", "+0.01", class = "btn-sm btn-outline-secondary", style = "width: 100%;"))
                   ),
                   tags$small(class = "text-muted", "Use buttons for quick adjustments, or type exact value")
            )
          ),
          br(),
          sliderInput("mt_highlight_vertical_offset", "Horizontal Offset (left/right)",
                      min = -1, max = 1, value = 0, step = 0.01),
          sliderInput("mt_highlight_adjust_height", "Ellipse Height",
                      min = 0.1, max = 3, value = 1, step = 0.05),
          sliderInput("mt_highlight_adjust_width", "Ellipse Width",
                      min = 0.1, max = 5, value = 1.5, step = 0.1),
          hr(),
          actionButton("mt_apply_highlight", "Apply & Preview",
                       icon = icon("eye"), class = "btn-primary btn-block"),
          actionButton("mt_remove_highlight", "Remove Highlight",
                       icon = icon("minus"), class = "btn-danger btn-block btn-sm"),
          hr(),
          tags$h5("Per-Tree Ellipse Adjustments"),
          uiOutput("mt_per_tree_highlight_ui")
        )
      ),

      box(
        title = "Highlight Values Settings",
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        uiOutput("mt_highlight_values_settings_ui")
      )
    ),

    fluidRow(
      box(
        title = tagList(
          "Preview ",
          span(id = "mt_high_status_waiting",
            style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
            icon("clock"), " Waiting for data"
          ),
          span(id = "mt_high_status_processing",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("spinner", class = "fa-spin"), " Processing..."
          ),
          span(id = "mt_high_status_ready",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("check-circle"), " Ready"
          )
        ),
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        imageOutput("mt_highlight_preview", height = "auto")
      )
    )
  )
}

mt_tabItem_legend <- function() {
  font_choices <- c("Sans-serif (default)" = "sans", "Serif" = "serif", "Monospace" = "mono")
  if (exists("SHOWTEXT_AVAILABLE") && SHOWTEXT_AVAILABLE) {
    font_choices <- c(font_choices,
      "Roboto" = "Roboto", "Open Sans" = "Open Sans", "Lato" = "Lato",
      "Montserrat" = "Montserrat", "PT Sans" = "PT Sans",
      "Playfair Display" = "Playfair Display", "Merriweather" = "Merriweather",
      "Source Code Pro" = "Source Code Pro")
  }

  tabItem(tabName = "mt_legend",
    fluidRow(
      column(4,
        box(title = NULL, status = "primary", solidHeader = FALSE, width = 12,
          tags$h4(icon("arrows-alt"), " Legend Position", style = "margin-top: 0;"),
          tags$p(class = "text-muted", "Where legends appear on the plot:"),
          selectInput("mt_legend_position", NULL,
                      choices = c("Right (default)" = "right", "Left" = "left",
                                  "Top" = "top", "Bottom" = "bottom",
                                  "None (hide all)" = "none"),
                      selected = "right")
        ),
        box(title = NULL, status = "info", solidHeader = FALSE, width = 12,
          tags$h4(icon("eye"), " Legend Visibility", style = "margin-top: 0;"),
          tags$p(class = "text-muted", "Toggle which legends to show:"),
          checkboxInput("mt_legend_show_classification", "Classification Legend", value = TRUE),
          checkboxInput("mt_legend_show_highlight", "Highlight Legend", value = TRUE),
          checkboxInput("mt_legend_show_bootstrap", "Bootstrap Legend", value = TRUE),
          checkboxInput("mt_legend_show_pvalue", "P Value Legend", value = TRUE)
        ),
        box(title = NULL, status = "primary", solidHeader = FALSE, width = 12,
          tags$h4(icon("font"), " Font Settings", style = "margin-top: 0;"),
          sliderInput("mt_legend_title_size", "Title Size:", min = 4, max = 48, value = 12, step = 1),
          sliderInput("mt_legend_text_size", "Text Size:", min = 2, max = 36, value = 10, step = 1),
          selectInput("mt_legend_font_family", "Font Family:", choices = font_choices, selected = "sans"),
          tags$small(class = "text-muted",
            if (exists("SHOWTEXT_AVAILABLE") && SHOWTEXT_AVAILABLE) "Google Fonts available via showtext"
            else "Install 'showtext' package for more fonts")
        ),
        box(title = NULL, status = "primary", solidHeader = FALSE, width = 12,
          tags$h4(icon("th"), " Symbol & Spacing", style = "margin-top: 0;"),
          sliderInput("mt_legend_key_size", "Key Size:", min = 0.1, max = 5, value = 1, step = 0.1),
          sliderInput("mt_legend_key_width", "Key Width:", min = 0.5, max = 3, value = 1, step = 0.1),
          sliderInput("mt_legend_key_height", "Key Height:", min = 0.5, max = 3, value = 1, step = 0.1),
          sliderInput("mt_legend_spacing", "Horizontal Spacing:", min = 0.05, max = 3, value = 0.3, step = 0.05),
          sliderInput("mt_legend_spacing_vertical", "Vertical Spacing:", min = 0.1, max = 5, value = 1, step = 0.1),
          sliderInput("mt_legend_title_key_spacing", "Title to Keys Spacing:", min = 0, max = 2, value = 0.2, step = 0.05),
          sliderInput("mt_legend_key_spacing", "Between Keys Spacing:", min = 0, max = 2, value = 0.1, step = 0.05),
          checkboxInput("mt_legend_reverse_order", "Reverse item order in legends", value = FALSE)
        ),
        box(title = NULL, status = "warning", solidHeader = FALSE, width = 12,
          collapsible = TRUE, collapsed = TRUE,
          tags$h4(icon("palette"), " Legend Background", style = "margin-top: 0;"),
          tags$p(class = "text-muted", tags$small("Legend box appearance:")),
          fluidRow(
            column(6,
              colourpicker::colourInput("mt_legend_box_background", "Box Background",
                                        value = "transparent", showColour = "both",
                                        allowTransparent = TRUE)
            ),
            column(6,
              sliderInput("mt_legend_margin", "Legend Margin",
                          min = 0, max = 2, value = 0.2, step = 0.05)
            )
          )
        ),
        box(title = NULL, status = "primary", solidHeader = FALSE, width = 12,
          tags$h4(icon("sitemap"), " Legend Scope", style = "margin-top: 0;"),
          radioButtons("mt_legend_scope", NULL,
                       choices = c("Same legend for all trees" = "all",
                                   "Legend per individual tree" = "per_tree"),
                       selected = "all"),
          conditionalPanel(
            condition = "input.mt_legend_scope == 'per_tree'",
            selectInput("mt_legend_tree_selector", "Select Tree:",
                        choices = NULL, selected = NULL)
          )
        ),
        box(title = NULL, status = "primary", solidHeader = FALSE, width = 12,
          actionButton("mt_apply_legend", "Apply Legend Settings",
                       class = "btn-success btn-lg btn-block", icon = icon("check"))
        )
      ),
      column(8,
        box(
          title = tagList(
            "Legend Preview ",
            span(id = "mt_legend_status_waiting",
              style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
              icon("clock"), " Waiting for data"
            ),
            span(id = "mt_legend_status_processing",
              style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
              icon("spinner", class = "fa-spin"), " Processing..."
            ),
            span(id = "mt_legend_status_ready",
              style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
              icon("check-circle"), " Ready"
            )
          ),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          imageOutput("mt_legend_preview", height = "auto")
        )
      )
    )
  )
}

mt_tabItem_extra <- function() {
  tabItem(tabName = "mt_extra",
    fluidRow(
      box(
        title = tagList(
          icon("arrows-alt"), " Plot Position ",
          span(id = "mt_extra_status_waiting",
            style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
            icon("clock"), " Waiting for data"),
          span(id = "mt_extra_status_processing",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("spinner", class = "fa-spin"), " Processing..."),
          span(id = "mt_extra_status_ready",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px;",
            icon("check"), " Ready")
        ),
        status = "primary", solidHeader = TRUE, width = 6,
        selectInput("mt_extra_tree_selector", "Apply to Tree:",
                    choices = c("All Trees" = "all"), selected = "all"),
        tags$p(class = "text-muted", "Move the plot on the page. Select a specific tree or all."),
        fluidRow(
          column(6, sliderInput("mt_plot_offset_x", "Horizontal Position:",
                                min = -5, max = 5, value = 0, step = 0.1, post = " (left/right)")),
          column(6, sliderInput("mt_plot_offset_y", "Vertical Position:",
                                min = -10, max = 5, value = 0, step = 0.1, post = " (down/up)"))
        ),
        fluidRow(column(12, actionButton("mt_reset_plot_position", "Reset to Center",
                                          class = "btn-secondary btn-sm", icon = icon("undo")))),
        hr(),
        tags$p(class = "text-muted", "Scale the entire plot up or down (zoom)."),
        sliderInput("mt_plot_scale_percent", "Plot Scale:", min = 25, max = 200, value = 100, step = 5, post = "%"),
        fluidRow(column(12, actionButton("mt_reset_plot_scale", "Reset to 100%",
                                          class = "btn-secondary btn-sm", icon = icon("undo")))),
        hr(),
        tags$p(class = "text-muted", tags$strong("Background:"), " Set page background color."),
        fluidRow(
          column(6, colourpicker::colourInput("mt_background_color", "Background Color:",
                                              value = "#FFFFFF", showColour = "both", allowTransparent = TRUE)),
          column(6, actionButton("mt_reset_background", "Reset to White",
                                  class = "btn-secondary btn-sm", icon = icon("undo"), style = "margin-top: 25px;"))
        ),
        hr(),
        checkboxInput("mt_a4_output", "A4 Page Format (297x210mm)", value = FALSE)
      ),
      box(title = "Preview", status = "primary", solidHeader = TRUE, width = 6,
        imageOutput("mt_extra_preview", height = "auto"),
        actionButton("mt_extra_apply", "Apply to Plot", class = "btn-primary", style = "margin-top: 10px;")
      )
    ),

    fluidRow(
      box(
        title = tagList(icon("heading"), " Titles"),
        status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = TRUE,
        fluidRow(
          column(6,
            tags$h5("Page Title"),
            checkboxInput("mt_enable_page_title", "Enable Page Title", value = FALSE),
            conditionalPanel(
              condition = "input.mt_enable_page_title == true",
              textInput("mt_page_title_text", "Title Text:", value = ""),
              fluidRow(
                column(6, numericInput("mt_page_title_x", "X Position:", value = 0.5, min = 0, max = 1, step = 0.01)),
                column(6, numericInput("mt_page_title_y", "Y Position:", value = 0.95, min = 0, max = 1, step = 0.01))
              ),
              sliderInput("mt_page_title_size", "Font Size:", min = 6, max = 72, value = 18, step = 1),
              colourpicker::colourInput("mt_page_title_color", "Color:", value = "#000000"),
              fluidRow(
                column(6, checkboxInput("mt_page_title_bold", "Bold", value = TRUE)),
                column(6, selectInput("mt_page_title_hjust", "Alignment:",
                                      choices = c("Left" = "0", "Center" = "0.5", "Right" = "1"), selected = "0.5"))
              )
            )
          ),
          column(6,
            tags$h5("Per-Tree Titles"),
            tags$p(class = "text-muted", "Change titles for individual trees."),
            uiOutput("mt_per_tree_extra_ui")
          )
        )
      )
    ),

    fluidRow(
      box(
        title = tagList(icon("font"), " Text Overlay"),
        status = "success", solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = FALSE,
        tags$div(
          style = "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
          tags$p(style = "margin: 0; color: #155724;", icon("info-circle"),
            " Text is placed as an ", tags$b("overlay"), " on top of the plot.")
        ),
        fluidRow(
          column(3, textInput("mt_custom_text_content", "Text:", value = "")),
          column(2, numericInput("mt_custom_text_x", "X Position:", value = 0.5, min = 0, max = 1, step = 0.01)),
          column(2, numericInput("mt_custom_text_y", "Y Position:", value = 0.5, min = 0, max = 1, step = 0.01)),
          column(2, sliderInput("mt_custom_text_size", "Size:", min = 4, max = 36, value = 12, step = 1)),
          column(2, colourpicker::colourInput("mt_custom_text_color", "Color:", value = "#000000"))
        ),
        fluidRow(
          column(3, selectInput("mt_custom_text_fontface", "Font Style:",
                                choices = c("Plain" = "plain", "Bold" = "bold", "Italic" = "italic", "Bold Italic" = "bold.italic"))),
          column(3, selectInput("mt_custom_text_hjust", "H. Align:",
                                choices = c("Left" = "0", "Center" = "0.5", "Right" = "1"), selected = "0.5")),
          column(3, selectInput("mt_custom_text_vjust", "V. Align:",
                                choices = c("Top" = "1", "Middle" = "0.5", "Bottom" = "0"), selected = "0.5")),
          column(3, numericInput("mt_custom_text_angle", "Rotation:", value = 0, min = -180, max = 180, step = 1))
        ),
        fluidRow(
          column(4, actionButton("mt_add_custom_text", "Add Text", class = "btn-success", icon = icon("plus"))),
          column(4, actionButton("mt_clear_custom_texts", "Clear All Texts", class = "btn-warning", icon = icon("trash")))
        ),
        hr(),
        h5("Added Texts:"),
        uiOutput("mt_custom_texts_list")
      )
    )
  )
}

mt_tabItem_config <- function() {
  tabItem(tabName = "mt_config", fluidRow(
    box(title = "Configuration", status = "primary",
        solidHeader = TRUE, width = 12,
        tags$p("Current multi-tree YAML configuration:"),
        verbatimTextOutput("mt_yaml_output"),
        downloadButton("mt_download_yaml_config", "Download Configuration", class = "btn-success")
    )
  ))
}

mt_tabItem_download <- function() {
  tabItem(tabName = "mt_download",
    fluidRow(
      box(title = "Output Options", status = "primary", solidHeader = TRUE, width = 4,
        textInput("mt_individual_name", "Sample/Individual Name:", value = "Sample1"),
        selectInput("mt_output_format", "File Format:",
                    choices = c("pdf", "png", "tiff", "svg", "jpeg"), selected = "pdf"),
        selectInput("mt_page_orientation", "Page Orientation:",
                    choices = c("Landscape" = "landscape", "Portrait" = "portrait"), selected = "landscape"),
        numericInput("mt_output_width", "Width:", value = 29.7),
        numericInput("mt_output_height", "Height:", value = 42),
        selectInput("mt_output_units", "Units:", choices = c("cm", "mm", "in"), selected = "cm"),
        textInput("mt_output_path", "Output Directory:", value = "./"),
        checkboxInput("mt_replace_name", "Use Custom File Name", value = FALSE),
        conditionalPanel(
          condition = "input.mt_replace_name == true",
          textInput("mt_custom_name", "Custom File Name:", value = "multi_tree_plot")
        ),
        textInput("mt_prefix_text", "Optional Text at Beginning:", value = "MultiTree__"),
        textInput("mt_suffix_text", "Optional Text at End:", value = "")
      ),
      box(
        title = tagList(
          "Final Preview ",
          span(id = "mt_download_status_waiting",
            style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
            icon("clock"), " Waiting for data"),
          span(id = "mt_download_status_processing",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("spinner", class = "fa-spin"), " Processing..."),
          span(id = "mt_download_status_ready",
            style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
            icon("check-circle"), " Ready")
        ),
        status = "primary", solidHeader = TRUE, width = 8,
        imageOutput("mt_download_preview", height = "auto"),
        downloadButton("mt_download_plot", "Download Plot", class = "btn-primary"),
        downloadButton("mt_download_yaml", "Download YAML Configuration", class = "btn-success"),
        tags$p(class = "text-muted", style = "margin-top: 5px;",
          tags$small("YAML file captures current settings for reproducibility"))
      )
    )
  )
}

# ================================================================
# HELPER FUNCTIONS
# ================================================================

# Build a yaml_data list matching the single-mode structure.
# No heatmap block — func.print.lineage.tree will skip heatmap rendering.
mt_build_yaml_data <- function(newick_path, csv_path, shared, per_tree) {
  # Classification according list
  according_list <- list()
  if (!is.null(shared$classification_groups) && length(shared$classification_groups) > 0) {
    for (i in seq_along(shared$classification_groups)) {
      grp <- shared$classification_groups[[i]]
      entry <- list()
      entry[[as.character(i)]] <- list(
        title1 = shared$classification_column,
        value1 = grp$values,
        display_name = grp$display_name,
        color = grp$color
      )
      according_list[[i]] <- entry
    }
  }

  # Must have at least one according entry (matches single-mode default from line 12972)
  if (length(according_list) == 0) {
    default_acc <- list()
    default_acc[["1"]] <- list(
      title1 = if (!is.null(shared$id_column) && nzchar(shared$id_column)) shared$id_column else "ID",
      value1 = list("Default"),
      display_name = "Default",
      color = "gray"
    )
    according_list <- list(default_acc)
  }

  # Build classification block
  classification_entry <- list()
  classification_entry[["1"]] <- list(
    according = according_list,
    FDR_perc = shared$fdr_perc,
    non_cluster_title = "No cluster",
    non_cluster_color = if (shared$no_cluster_color) "gray" else "white",
    not_in_legend = list("No cluster", "0"),
    title = shared$classification_title,
    na_name = "na"
  )

  # Add highlight if enabled
  if (isTRUE(shared$enable_highlight) && !is.null(shared$highlight_yaml)) {
    # Adjust ellipse per-tree if available
    hl <- shared$highlight_yaml
    if (!is.null(per_tree$highlight_adjust_height)) {
      hl$adjust_height <- per_tree$highlight_adjust_height
    }
    if (!is.null(per_tree$highlight_adjust_width)) {
      hl$adjust_width <- per_tree$highlight_adjust_width
    }
    classification_entry[["1"]]$highlight <- hl
  }

  # Bootstrap
  bootstrap_block <- list(
    display = if (isTRUE(shared$enable_bootstrap)) "yes" else "no",
    format = shared$bootstrap_format,
    param = as.character(shared$bootstrap_param)
  )

  # Rotation - supports all three types: rotation1 (primary), rotation2 (secondary), manual
  rotation_nodes <- list()
  if (!is.null(per_tree$rotate) && length(per_tree$rotate) > 0) {
    for (k in seq_along(per_tree$rotate)) {
      node_entry <- list()
      node_entry[[as.character(k)]] <- per_tree$rotate[k]
      rotation_nodes[[k]] <- node_entry
    }
  }

  rotation_type <- per_tree$rotation_type
  manual_rotation <- list(
    display = if (length(rotation_nodes) > 0 &&
                  (is.null(rotation_type) || rotation_type == "manual")) "yes" else "no",
    nodes = rotation_nodes
  )

  # rotation1 (primary classification first)
  rot1_according <- list()
  if (!is.null(per_tree$rotation1_config) && length(per_tree$rotation1_config) > 0) {
    for (k in seq_along(per_tree$rotation1_config)) {
      cfg <- per_tree$rotation1_config[[k]]
      if (!is.null(cfg$col) && !is.null(cfg$val) && cfg$col != "" && cfg$val != "") {
        rot1_according[[as.character(k)]] <- list(title = cfg$col, value = cfg$val)
      }
    }
  }
  rotation1_block <- list(
    display = if (length(rot1_according) > 0 && !is.null(rotation_type) &&
                  (rotation_type == "primary" || rotation_type == "manual")) "yes" else "no",
    according = rot1_according
  )

  # rotation2 (secondary classification first)
  rot2_according <- list()
  if (!is.null(per_tree$rotation2_config) && length(per_tree$rotation2_config) > 0) {
    for (k in seq_along(per_tree$rotation2_config)) {
      cfg <- per_tree$rotation2_config[[k]]
      if (!is.null(cfg$col) && !is.null(cfg$val) && cfg$col != "" && cfg$val != "") {
        rot2_according[[as.character(k)]] <- list(title = cfg$col, value = cfg$val)
      }
    }
  }
  rotation2_block <- list(
    display = if (length(rot2_according) > 0 && !is.null(rotation_type) &&
                  rotation_type == "secondary") "yes" else "no",
    according = rot2_according
  )

  # Assemble the full yaml_data
  yaml_data <- list(
    "Individual general definitions" = list(
      Individual = shared$individual_name,
      "tree path" = list(newick_path),
      "mapping csv file" = csv_path,
      "individual column" = shared$individual_column,
      "out_file" = list(
        "base_path" = "./",
        "file_type" = "pdf",
        "optional text at beggining" = "",
        "optional text at end" = "",
        "replace name" = list(flag = "no", name = "tree_plot")
      )
    ),
    "Mapping exl renaming titles" = list(
      "ID column" = shared$id_column
    ),
    "visual definitions" = list(
      classification = list(classification_entry),
      Bootstrap = bootstrap_block,
      rotation1 = rotation1_block,
      rotation2 = rotation2_block,
      manual_rotation = manual_rotation,
      "trim tips" = list(
        display = if (isTRUE(shared$trim_tips)) "yes" else "no",
        length = shared$tip_length
      ),
      edge_width_multiplier = list(size = shared$edge_width),
      font_size = list(
        tips = shared$tip_font_size,
        legend_title = shared$legend_title_size,
        legend_text = shared$legend_text_size,
        legend_box = 15,
        heat_map_title = 25,
        heat_map_legend = 3.8
      ),
      laderize_flag = if (isTRUE(shared$ladderize)) "yes" else "no",
      compare_two_trees = "no",
      simulate_p_value = "yes",
      flag_calc_scores_for_tree = "no",
      flag_make_newick_file = "no",
      debug_mode = "no",
      id_tip_trim_flag = shared$trim_tips,
      id_tip_trim_start = 1,
      id_tip_trim_end = 100,
      id_tip_prefix = ""
    )
  )

  # Add man_* params from YAML import if available
  if (!is.null(shared$man_params)) {
    for (pname in names(shared$man_params)) {
      yaml_data[["visual definitions"]][[pname]] <- shared$man_params[[pname]]
    }
  }

  yaml_data
}

# Render multiple trees and combine into a single gtable.
# Returns a gtable suitable for plotOutput / ggsave. Does NOT save to disk.
func.render.single.tree.in.app <- function(
  newick_path, tree_name, csv_path,
  shared_settings, per_tree, tree_title, tree_bg_color,
  log_fn = function(msg) {}, tree_cache = list()
) {
  tree_shared <- shared_settings
  if (!is.null(per_tree$individual_value) && nzchar(per_tree$individual_value)) {
    tree_shared$individual_name <- per_tree$individual_value
  }
  if (!is.null(per_tree$classification_groups) && length(per_tree$classification_groups) > 0) {
    tree_shared$classification_groups <- per_tree$classification_groups
  }
  if (!is.null(per_tree$classification_title)) {
    tree_shared$classification_title <- per_tree$classification_title
  }
  if (!is.null(per_tree$no_cluster_color)) {
    tree_shared$no_cluster_color <- (per_tree$no_cluster_color == "gray")
  }
  yaml_data <- mt_build_yaml_data(
    newick_path = newick_path,
    csv_path = csv_path,
    shared = tree_shared,
    per_tree = per_tree
  )

  temp_yaml <- tempfile(fileext = ".yaml")
  writeLines(yaml::as.yaml(yaml_data, indent.mapping.sequence = TRUE), temp_yaml)

  treesi <- tryCatch({
    suppressWarnings(func.print.lineage.tree(
      conf_yaml_path = temp_yaml,
      csv_path_bash = csv_path,
      width = shared_settings$output_width,
      height = shared_settings$output_height,
      units_out = shared_settings$output_units,
      debug_mode = FALSE,
      compare_two_trees = FALSE,
      list_nodes_to_rotate = if (!is.null(per_tree$rotate) && length(per_tree$rotate) > 0) per_tree$rotate else NA,
      flag_display_nod_number_on_tree = isTRUE(shared_settings$display_node_numbers),
      node_number_font_size = shared_settings$node_number_font_size,
      bootstrap_label_size = shared_settings$bootstrap_label_size,
      man_boot_x_offset = shared_settings$man_boot_x_offset,
      heatmap_tree_distance = 0.02,
      heatmap_global_gap = 0.05,
      legend_settings = shared_settings$legend_settings,
      cached_p_list_of_pairs = tree_cache$p_list_of_pairs,
      cached_p_list_hash = tree_cache$p_list_hash,
      p_list_cache = if (!is.null(tree_cache$p_list_cache)) tree_cache$p_list_cache else list(),
      heatmap_cache = if (!is.null(tree_cache$heatmap_cache)) tree_cache$heatmap_cache else list()
    ))
  }, error = function(e) {
    tb <- paste(capture.output(traceback(4)), collapse = "\n")
    log_fn(paste0("ERROR rendering tree ", tree_name, ": ", e$message))
    log_fn(paste0("Traceback:\n", tb))
    NULL
  })

  if (is.null(treesi)) return(NULL)

  one_plot <- NULL
  if (!is.null(treesi$plots) && length(treesi$plots) > 0) {
    one_plot <- treesi$plots[[1]]
  } else if (is.list(treesi) && length(treesi) > 0) {
    one_plot <- treesi[[1]]
  }

  updated_cache <- tree_cache
  if (!is.null(treesi$cache_data)) {
    updated_cache <- list(
      p_list_of_pairs = treesi$cache_data$p_list_of_pairs,
      p_list_hash = treesi$cache_data$p_list_hash,
      p_list_cache = if (!is.null(tree_cache$p_list_cache)) tree_cache$p_list_cache else list(),
      heatmap_cache = treesi$cache_data$heatmap_cache
    )
    new_hash <- treesi$cache_data$p_list_hash
    if (!is.null(new_hash) && !is.null(treesi$cache_data$p_list_of_pairs)) {
      updated_cache$p_list_cache[[new_hash]] <- treesi$cache_data$p_list_of_pairs
    }
  }

  if (is.null(one_plot)) return(NULL)

  styled_plot <- one_plot +
    ggtree::theme_tree(bgcolor = tree_bg_color) +
    ggplot2::ggtitle(tree_title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
      plot.background = ggplot2::element_rect(fill = tree_bg_color, color = tree_bg_color)
    ) +
    ggplot2::guides(color = "none", size = "none", shape = "none")

  list(plot = styled_plot, cache = updated_cache)
}

func.multiple.trees.one.page.in.app <- function(
  newick_paths, newick_names, csv_path,
  shared_settings, per_tree_list,
  tree_titles, tree_bg_colors, log_fn = function(msg) {},
  p_cache_per_tree = list(),
  cached_tree_plots = list(),
  dirty_trees = NULL
) {
  list_out <- list()
  updated_cache <- p_cache_per_tree
  updated_plots <- cached_tree_plots
  all_dirty <- is.null(dirty_trees)

  for (i in seq_along(newick_paths)) {
    tn <- newick_names[i]
    per_tree <- if (!is.null(per_tree_list[[tn]])) per_tree_list[[tn]] else list()

    need_render <- all_dirty || (tn %in% dirty_trees) || is.null(cached_tree_plots[[tn]])

    if (!need_render) {
      log_fn(paste0("Using cached plot for tree ", i, "/", length(newick_paths), ": ", tn))
      list_out[[length(list_out) + 1]] <- cached_tree_plots[[tn]]
      next
    }

    log_fn(paste0("Rendering tree ", i, "/", length(newick_paths), ": ", tn))
    tree_cache <- if (!is.null(p_cache_per_tree[[tn]])) p_cache_per_tree[[tn]] else list()

    result_single <- func.render.single.tree.in.app(
      newick_path = newick_paths[i],
      tree_name = tn,
      csv_path = csv_path,
      shared_settings = shared_settings,
      per_tree = per_tree,
      tree_title = tree_titles[i],
      tree_bg_color = tree_bg_colors[i],
      log_fn = log_fn,
      tree_cache = tree_cache
    )

    if (is.null(result_single)) next

    list_out[[length(list_out) + 1]] <- result_single$plot
    updated_plots[[tn]] <- result_single$plot
    updated_cache[[tn]] <- result_single$cache
  }

  if (length(list_out) == 0) return(NULL)

  log_fn(paste0("Combining ", length(list_out), " trees into one page"))
  result <- do.call(gridExtra::grid.arrange, c(list_out, ncol = 1))
  attr(result, "p_cache_per_tree") <- updated_cache
  attr(result, "cached_tree_plots") <- updated_plots
  result
}

# ================================================================
# SERVER LOGIC INSTALLER
# Called once from the main server function.
# ================================================================
mt_install_server <- function(input, output, session) {

  # --- Reactive state (isolated from single-mode) ---
  mt_values <- reactiveValues(
    newick_paths = NULL,        # character vector of temp file paths
    newick_names = NULL,        # character vector of original filenames (sans extension)
    slot_files = list(),        # list keyed by slot number, each = list(name, path, tree_name)
    csv_path = NULL,            # single temp file path
    csv_data = NULL,            # data.frame
    trees = list(),             # list of treedata objects per newick
    matched_per_tree = list(),  # list[[tree_name]] = filtered CSV subset for that tree
    filtered_csv = NULL,        # combined filtered CSV (only matched rows across all trees)
    temp_csv_path = NULL,        # path to filtered CSV on disk (passed to func.print.lineage.tree)
    per_tree = list(),          # list[[tree_name]] = list(rotate, highlight_adjust_height, highlight_adjust_width, title, bg_color)
    man_params = list(),        # man_* params from YAML import
    highlight_yaml = NULL,      # shared highlight YAML block
    classification_groups = list(), # classification groups for YAML builder
    temp_classification_preview = NULL, # temp classification for preview
    last_plot = NULL,           # last rendered gtable
    last_plot_file = NULL,      # path to cached SVG of last render
    plot_counter = 0,           # incremented after each successful render
    plot_generating = FALSE,    # recursion guard
    log_messages = "",          # log text for mt_log output
    p_cache_per_tree = list(),  # per-tree p-value cache: list[[tree_name]] = list(p_list, hash, p_list_cache)
    rotation1_config_per_tree = list(), # per-tree rotation1 config
    rotation2_config_per_tree = list(), # per-tree rotation2 config
    manual_rotation_per_tree = list(),  # per-tree manual rotation nodes
    cached_tree_plots = list(), # per-tree ggplot objects for incremental rendering
    dirty_trees = NULL          # character vector of tree names needing re-render (NULL = all dirty)
  )

  # Helper to append log messages
  mt_log <- function(msg) {
    ts <- format(Sys.time(), "%H:%M:%S")
    mt_values$log_messages <- paste0(mt_values$log_messages, "[", ts, "] ", msg, "\n")
  }

  # --- Dynamic upload slots for newick trees ---
  mt_slot_count <- reactiveVal(1)
  mt_observed_slots <- reactiveVal(integer(0))

  # "+ Add another tree" button adds a new upload widget
  observeEvent(input$mt_add_upload_slot, ignoreInit = TRUE, {
    mt_slot_count(mt_slot_count() + 1)
  })

  # Render the file input slots
  output$mt_tree_upload_slots <- renderUI({
    n <- mt_slot_count()
    slot_uis <- lapply(seq_len(n), function(i) {
      file_id <- paste0("mt_tree_slot_", i)
      remove_id <- paste0("mt_remove_slot_", i)
      # Check if this slot already has a file loaded
      loaded_name <- NULL
      if (!is.null(mt_values$slot_files[[as.character(i)]])) {
        loaded_name <- mt_values$slot_files[[as.character(i)]]$name
      }
      tags$div(
        style = "display: flex; align-items: flex-start; gap: 5px; margin-bottom: 5px;",
        tags$div(style = "flex: 1;",
          if (!is.null(loaded_name)) {
            tags$div(
              style = "padding: 8px 12px; background: #d4edda; border-radius: 4px; border: 1px solid #c3e6cb;",
              tags$span(style = "color: #155724;", icon("check"), paste0(" ", loaded_name))
            )
          } else {
            fileInput(file_id, NULL, multiple = FALSE,
                      accept = c(".nwk", ".newick", ".tree", ".nw", ".txt"),
                      width = "100%")
          }
        ),
        if (!is.null(loaded_name)) {
          actionButton(remove_id, "", icon = icon("times"),
                       class = "btn-xs btn-danger", style = "margin-top: 8px; padding: 2px 7px;")
        }
      )
    })
    tagList(slot_uis)
  })

  # Watch each slot's fileInput for uploads
  observe({
    n <- mt_slot_count()
    already <- mt_observed_slots()
    new_slots <- setdiff(seq_len(n), already)
    for (i in new_slots) {
      local({
        local_i <- i
        file_id <- paste0("mt_tree_slot_", local_i)
        remove_id <- paste0("mt_remove_slot_", local_i)

        # File upload observer
        observeEvent(input[[file_id]], {
          file_info <- input[[file_id]]
          req(file_info)
          tn <- tools::file_path_sans_ext(file_info$name)

          # Copy to persistent temp
          perm_path <- file.path(tempdir(), paste0("mt_slot", local_i, "_", basename(file_info$datapath)))
          file.copy(file_info$datapath, perm_path, overwrite = TRUE)

          # Ensure unique name
          existing_names <- mt_values$newick_names
          if (tn %in% existing_names) {
            suffix <- 2
            while (paste0(tn, "_", suffix) %in% existing_names) suffix <- suffix + 1
            tn <- paste0(tn, "_", suffix)
          }

          # Store in slot_files for UI rendering
          sf <- mt_values$slot_files
          sf[[as.character(local_i)]] <- list(name = file_info$name, path = perm_path, tree_name = tn)
          mt_values$slot_files <- sf

          # Rebuild newick paths/names from all slots
          mt_rebuild_from_slots()

          mt_log(paste0("Uploaded tree '", tn, "' in slot ", local_i,
                        ". Total: ", length(mt_values$newick_names), " trees"))
        }, ignoreInit = TRUE)

        # Remove button observer
        observeEvent(input[[remove_id]], {
          sf <- mt_values$slot_files
          old_tn <- sf[[as.character(local_i)]]$tree_name
          sf[[as.character(local_i)]] <- NULL
          mt_values$slot_files <- sf
          if (!is.null(old_tn)) mt_values$per_tree[[old_tn]] <- NULL
          mt_rebuild_from_slots()
          mt_log(paste0("Removed tree from slot ", local_i,
                        ". Remaining: ", length(mt_values$newick_names)))
        }, ignoreInit = TRUE)
      })
    }
    mt_observed_slots(seq_len(n))
  })

  # Rebuild newick_paths/newick_names from slot_files
  mt_rebuild_from_slots <- function() {
    sf <- mt_values$slot_files
    paths <- c()
    names <- c()
    for (key in sort(as.numeric(names(sf)))) {
      entry <- sf[[as.character(key)]]
      if (!is.null(entry)) {
        paths <- c(paths, entry$path)
        names <- c(names, entry$tree_name)
      }
    }
    mt_values$newick_paths <- if (length(paths) > 0) paths else NULL
    mt_values$newick_names <- if (length(names) > 0) names else NULL

    # Init per-tree state for any new names
    for (tn in names) {
      if (is.null(mt_values$per_tree[[tn]])) {
        mt_values$per_tree[[tn]] <- list(
          rotate = NULL, highlight_adjust_height = 1,
          highlight_adjust_width = 1.5, title = tn, bg_color = "#FFFFFF",
          individual_value = NULL
        )
      }
    }

    # Update tree selector
    if (length(names) > 0) {
      updateSelectInput(session, "mt_tree_selector",
                        choices = names, selected = names[length(names)])
    } else {
      updateSelectInput(session, "mt_tree_selector", choices = NULL)
    }
  }

  # --- CSV upload observer ---
  # Mirrors single-mode CSV read (lines 12346-12373): fread + column filtering
  observeEvent(input$mt_csv_file, ignoreInit = TRUE, {
    req(input$mt_csv_file)
    mt_values$csv_path <- input$mt_csv_file$datapath
    csv_data_raw <- tryCatch({
      dt <- data.table::fread(input$mt_csv_file$datapath, data.table = FALSE)
      names(dt) <- make.names(names(dt), unique = TRUE)
      dt
    }, error = function(e) {
      mt_log(paste0("ERROR reading CSV: ", e$message))
      NULL
    })
    if (!is.null(csv_data_raw)) {
      # Filter out empty/auto-generated columns (matches single-mode lines 12358-12373)
      col_names <- names(csv_data_raw)
      valid_cols <- !grepl("^V\\d+$", col_names) &
                    !grepl("^\\.\\.\\.", col_names) &
                    !grepl("^X(\\.\\d+)?$", col_names) &
                    col_names != "" & !is.na(col_names)
      removed_count <- sum(!valid_cols)
      if (removed_count > 0) {
        mt_log(paste0("Removed ", removed_count, " empty/auto-named columns (keeping ", sum(valid_cols), ")"))
        csv_data_raw <- csv_data_raw[, valid_cols, drop = FALSE]
      }
      mt_values$csv_data <- csv_data_raw

      cols <- names(mt_values$csv_data)
      mt_log(paste0("CSV loaded: ", nrow(mt_values$csv_data), " rows, ",
                    length(cols), " columns"))
      updateSelectizeInput(session, "mt_id_column", choices = cols,
                           selected = cols[1], server = TRUE)
      updateSelectizeInput(session, "mt_individual_column", choices = c("", cols),
                           server = TRUE)
      updateSelectizeInput(session, "mt_classification_column", choices = c("Select column..." = "", cols),
                           server = TRUE)
      updateSelectizeInput(session, "mt_highlight_column", choices = c("Select column..." = "", cols),
                           server = TRUE)
    }
  })

  # --- Individual column observer: populate per-tree individual dropdowns ---
  # Mirrors single-mode logic from lines 12434-12479 of main file.
  mt_individual_choices <- reactive({
    req(mt_values$csv_data, input$mt_individual_column)
    if (input$mt_individual_column == "") return(character(0))
    col <- input$mt_individual_column
    if (!(col %in% names(mt_values$csv_data))) return(character(0))
    vals <- unique(mt_values$csv_data[[col]])
    vals <- vals[!is.na(vals)]
    as.character(sort(vals))
  })

  # Per-tree individual mapping UI — uses selectizeInput like single mode
  output$mt_per_tree_individual_ui <- renderUI({
    req(mt_values$newick_names, input$mt_individual_column)
    if (input$mt_individual_column == "") {
      return(tags$p(class = "text-muted",
                    tags$em("Select an Individual Column first.")))
    }
    names <- mt_values$newick_names
    if (is.null(names) || length(names) == 0) {
      return(tags$p(class = "text-muted",
                    tags$em("Upload trees first.")))
    }
    choices <- mt_individual_choices()
    if (length(choices) == 0) {
      return(tags$p(class = "text-muted",
                    tags$em("No individual values found in CSV.")))
    }
    rows <- lapply(seq_along(names), function(i) {
      tn <- names[i]
      sel_id <- paste0("mt_tree_indiv_", i)
      current_val <- isolate(mt_values$per_tree[[tn]]$individual_value)
      tags$div(
        style = "display: flex; align-items: center; gap: 10px; padding: 5px;",
        tags$span(style = "min-width: 200px; font-weight: bold;", paste0(i, ". ", tn)),
        tags$span(style = "font-size: 16px;", HTML("&rarr;")),
        tags$div(style = "flex: 1;",
          selectizeInput(sel_id, NULL,
                         choices = c("(none)" = "", choices),
                         selected = if (!is.null(current_val) && current_val %in% choices) current_val else "",
                         width = "100%",
                         options = list(placeholder = "Select an individual..."))
        )
      )
    })
    tagList(rows)
  })

  # Observer: watch each per-tree individual dropdown and save selection
  mt_indiv_observers_installed <- reactiveVal(integer(0))

  observe({
    names <- mt_values$newick_names
    if (is.null(names)) return()
    already <- mt_indiv_observers_installed()
    new_indices <- setdiff(seq_along(names), already)
    for (i in new_indices) {
      local({
        local_i <- i
        sel_id <- paste0("mt_tree_indiv_", local_i)
        observeEvent(input[[sel_id]], {
          nn <- mt_values$newick_names
          if (local_i <= length(nn)) {
            tn <- nn[local_i]
            val <- input[[sel_id]]
            if (!is.null(mt_values$per_tree[[tn]])) {
              mt_values$per_tree[[tn]]$individual_value <-
                if (is.null(val) || val == "") NULL else val
            }
          }
        }, ignoreInit = TRUE)
      })
    }
    mt_indiv_observers_installed(union(already, new_indices))
  })

  # Normalize a name for fuzzy comparison: lowercase, replace hyphens/underscores/spaces
  mt_normalize_name <- function(x) {
    x <- tolower(x)
    gsub("[-_ ]+", "", x)
  }

  # Auto-match individuals by filename
  observeEvent(input$mt_auto_match_individuals, ignoreInit = TRUE, {
    req(mt_values$newick_names, mt_values$csv_data, input$mt_individual_column)
    if (input$mt_individual_column == "") return()
    choices <- mt_individual_choices()
    mt_log(paste0("Auto-match: ", length(choices), " individual values in CSV column '",
                  input$mt_individual_column, "'"))
    mt_log(paste0("  CSV values: ", paste(choices, collapse = ", ")))
    mt_log(paste0("  Tree names: ", paste(mt_values$newick_names, collapse = ", ")))

    # Pre-compute normalized CSV values
    choices_norm <- mt_normalize_name(choices)

    matched <- 0
    for (i in seq_along(mt_values$newick_names)) {
      tn <- mt_values$newick_names[i]
      match_val <- NULL

      # Strip common suffixes first
      tn_stripped <- gsub("_bootstrapped.*$|_rerooted.*$|_rooted.*$", "", tn, ignore.case = TRUE)
      tn_stripped <- sub("\\(\\d+\\)$", "", tn_stripped)
      tn_stripped <- trimws(tn_stripped)
      tn_norm <- mt_normalize_name(tn_stripped)

      mt_log(paste0("  Tree '", tn, "' -> stripped '", tn_stripped,
                    "' -> normalized '", tn_norm, "'"))

      # 1. Exact match (original or stripped)
      if (tn %in% choices) {
        match_val <- tn
      } else if (tn_stripped %in% choices) {
        match_val <- tn_stripped
      }

      # 2. Normalized match (handles hyphens vs underscores: 775-13 == 775_13)
      if (is.null(match_val)) {
        idx <- which(choices_norm == tn_norm)
        if (length(idx) > 0) {
          match_val <- choices[idx[1]]
        }
      }

      # 3. Normalized substring: CSV value in tree name or tree in CSV value
      if (is.null(match_val)) {
        for (j in order(nchar(choices_norm), decreasing = TRUE)) {
          cv_norm <- choices_norm[j]
          if (nzchar(cv_norm) && (grepl(cv_norm, tn_norm, fixed = TRUE) ||
                                   grepl(tn_norm, cv_norm, fixed = TRUE))) {
            match_val <- choices[j]
            break
          }
        }
      }

      if (!is.null(match_val)) {
        mt_values$per_tree[[tn]]$individual_value <- match_val
        sel_id <- paste0("mt_tree_indiv_", i)
        updateSelectizeInput(session, sel_id, selected = match_val)
        matched <- matched + 1
        mt_log(paste0("  -> MATCHED '", match_val, "'"))
      } else {
        mt_log(paste0("  -> NO MATCH"))
      }
    }
    mt_log(paste0("Auto-matched ", matched, "/", length(mt_values$newick_names),
                  " trees to individuals"))
    showNotification(paste0("Auto-matched ", matched, "/",
                            length(mt_values$newick_names), " trees"),
                     type = if (matched > 0) "message" else "warning")
  })

  # --- files_loaded output: TRUE when at least one tree + CSV are loaded ---
  output$mt_files_loaded <- reactive({
    !is.null(mt_values$newick_paths) && length(mt_values$newick_paths) > 0 &&
      !is.null(mt_values$csv_data)
  })
  outputOptions(output, "mt_files_loaded", suspendWhenHidden = FALSE)

  # --- Tree summary output ---
  output$mt_tree_summary <- renderText({
    if (is.null(mt_values$newick_names) || length(mt_values$newick_names) == 0) {
      return("No trees uploaded.")
    }
    paste0(length(mt_values$newick_names), " tree(s) uploaded: ",
           paste(mt_values$newick_names, collapse = ", "))
  })

  # --- CSV summary output ---
  output$mt_csv_summary <- renderText({
    if (is.null(mt_values$csv_data)) return("No CSV loaded.")
    paste0("CSV: ", nrow(mt_values$csv_data), " rows, ",
           ncol(mt_values$csv_data), " columns")
  })

  # --- Highlight column observer (matches single mode lines 14225-14270) ---
  observeEvent(input$mt_highlight_column, ignoreInit = TRUE, {
    if (is.null(input$mt_highlight_column) || input$mt_highlight_column == "") {
      updateSelectizeInput(session, "mt_highlight_values", choices = character(0), selected = character(0))
      output$mt_highlight_values_settings_ui <- renderUI({
        tags$p(style = "color: gray; font-style: italic;", "Please select a column first")
      })
      return(NULL)
    }

    req(mt_values$csv_data, input$mt_highlight_column)
    if (!(input$mt_highlight_column %in% names(mt_values$csv_data))) return(NULL)

    csv_to_use <- if (!is.null(mt_values$filtered_csv) && nrow(mt_values$filtered_csv) > 0) {
      mt_values$filtered_csv
    } else {
      mt_values$csv_data
    }

    unique_values <- unique(csv_to_use[[input$mt_highlight_column]])
    unique_values <- unique_values[!is.na(unique_values)]

    if (length(unique_values) > 100) {
      output$mt_highlight_values_settings_ui <- renderUI({
        tags$div(
          style = "padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
          tags$h5("Too many unique values!", style = "color: #856404; margin-top: 0;"),
          tags$p(paste("This column has", length(unique_values), "unique values."))
        )
      })
      updateSelectizeInput(session, "mt_highlight_values", choices = character(0))
      return(NULL)
    }

    updateSelectizeInput(session, "mt_highlight_values", choices = unique_values, selected = character(0), server = TRUE)

    output$mt_highlight_values_settings_ui <- renderUI({
      tags$div(
        tags$p(tags$strong(paste("Available values:", length(unique_values)))),
        tags$p("Select values from the dropdown to highlight them."),
        tags$p(style = "color: #666;", "After selecting, you can customize color and settings for each value below.")
      )
    })
  }, ignoreInit = TRUE)

  # --- Process & Match IDs ---
  observeEvent(input$mt_process_data, ignoreInit = TRUE, {
    req(mt_values$csv_data, mt_values$newick_paths)
    mt_log("Processing and matching IDs for each tree...")

    csv_full <- mt_values$csv_data
    id_col <- input$mt_id_column
    indiv_col <- input$mt_individual_column
    use_all <- isTRUE(input$mt_use_all_data)

    # Read each newick and match against its own CSV subset
    mt_values$trees <- list()
    mt_values$matched_per_tree <- list()
    total_ok <- 0
    total_missing <- 0

    for (i in seq_along(mt_values$newick_paths)) {
      tn <- mt_values$newick_names[i]
      tree <- tryCatch(
        treeio::read.newick(mt_values$newick_paths[i]),
        error = function(e) {
          mt_log(paste0("ERROR reading newick '", tn, "': ", e$message))
          NULL
        }
      )
      if (is.null(tree)) next
      mt_values$trees[[tn]] <- tree

      # Determine CSV subset for this tree
      if (use_all || is.null(indiv_col) || indiv_col == "") {
        csv_sub <- csv_full
        subset_note <- "all rows"
      } else {
        iv <- mt_values$per_tree[[tn]]$individual_value
        if (is.null(iv) || iv == "") {
          mt_log(paste0("Tree '", tn, "': no individual selected, skipping"))
          total_missing <- total_missing + 1
          next
        }
        csv_sub <- csv_full[csv_full[[indiv_col]] == iv, , drop = FALSE]
        subset_note <- paste0("individual '", iv, "' (", nrow(csv_sub), " rows)")
      }

      mt_values$matched_per_tree[[tn]] <- csv_sub

      tip_labels <- tree@phylo$tip.label
      matched <- sum(tip_labels %in% csv_sub[[id_col]])
      mt_log(paste0("Tree '", tn, "': ", length(tip_labels), " tips, ",
                    matched, " matched to ", subset_note))
      total_ok <- total_ok + 1
    }

    # Build filtered_csv: combine all matched per-tree subsets, filtered to matched IDs
    all_matched <- do.call(rbind, mt_values$matched_per_tree)
    if (!is.null(all_matched) && nrow(all_matched) > 0 && !is.null(id_col) && id_col %in% names(all_matched)) {
      all_tip_labels <- c()
      for (tn in names(mt_values$trees)) {
        all_tip_labels <- c(all_tip_labels, mt_values$trees[[tn]]@phylo$tip.label)
      }
      all_tip_labels <- unique(all_tip_labels)
      mt_values$filtered_csv <- all_matched[all_matched[[id_col]] %in% all_tip_labels, , drop = FALSE]
      mt_log(paste0("Filtered CSV: ", nrow(mt_values$filtered_csv), " rows (matched tree tips)"))
    } else {
      mt_values$filtered_csv <- all_matched
    }

    # Write filtered CSV to disk (matches single-mode pattern from line 13830)
    # func.print.lineage.tree reads CSV from disk, so we must give it the filtered version
    if (!is.null(mt_values$filtered_csv) && nrow(mt_values$filtered_csv) > 0) {
      mt_values$temp_csv_path <- tempfile(fileext = ".csv")
      csv_to_write <- mt_values$filtered_csv
      col_names <- names(csv_to_write)
      valid_cols <- !grepl("^\\.\\.\\.", col_names) & col_names != "" & !is.na(col_names)
      if (sum(valid_cols) < length(col_names)) {
        mt_log(paste0("Filtering out ", sum(!valid_cols), " empty/unnamed columns from temp CSV"))
        csv_to_write <- csv_to_write[, valid_cols, drop = FALSE]
      }
      write.csv(csv_to_write, mt_values$temp_csv_path, row.names = FALSE)
      mt_log(paste0("Temp CSV written: ", nrow(csv_to_write), " rows x ", ncol(csv_to_write),
                    " cols -> ", mt_values$temp_csv_path))
    }

    # Auto-generate classification groups from full CSV (shared across trees)
    if (!is.null(input$mt_classification_column) && input$mt_classification_column != "") {
      class_col <- input$mt_classification_column
      unique_vals <- sort(unique(na.omit(csv_full[[class_col]])))
      n_colors <- max(3, length(unique_vals))
      palette <- tryCatch(
        RColorBrewer::brewer.pal(min(n_colors, 12), "Set3"),
        error = function(e) rainbow(n_colors)
      )
      if (length(unique_vals) > length(palette)) {
        palette <- rep_len(palette, length(unique_vals))
      }
      groups <- list()
      for (j in seq_along(unique_vals)) {
        groups[[j]] <- list(
          values = list(unique_vals[j]),
          display_name = as.character(unique_vals[j]),
          color = palette[j]
        )
      }
      mt_values$classification_groups <- groups
      mt_log(paste0("Classification: ", length(unique_vals),
                    " groups from column '", class_col, "'"))
    }

    if (total_missing > 0) {
      showNotification(paste0("Processed ", total_ok, " trees. ",
                              total_missing, " trees missing individual mapping."),
                       type = "warning", duration = 5)
    } else {
      showNotification("Multi-tree processing complete!", type = "message")
    }
  })

  # --- Tree selector: save/load per-tree state ---
  # Store current per-tree values when switching away
  mt_prev_tree <- reactiveVal(NULL)

  # Save current tree's rotation state before switching
  mt_save_current_tree_rotation <- function() {
    prev <- mt_prev_tree()
    if (is.null(prev) || is.null(mt_values$per_tree[[prev]])) return()
    # Save rotation type
    mt_values$per_tree[[prev]]$rotation_type <- input$mt_rotation_type
    mt_values$per_tree[[prev]]$enable_rotation <- isTRUE(input$mt_enable_rotation)
    # Manual rotation nodes saved by apply button
  }

  observeEvent(input$mt_tree_selector, ignoreInit = TRUE, {
    req(input$mt_tree_selector)

    # Save previous tree's rotation state
    mt_save_current_tree_rotation()

    # Load new tree's per-tree values
    current <- input$mt_tree_selector
    mt_prev_tree(current)

    pt <- mt_values$per_tree[[current]]
    if (!is.null(pt)) {
      # Restore rotation settings
      updateCheckboxInput(session, "mt_enable_rotation",
                          value = isTRUE(pt$enable_rotation))
      if (!is.null(pt$rotation_type)) {
        updateRadioButtons(session, "mt_rotation_type",
                           selected = pt$rotation_type)
      }
      # Restore manual nodes
      if (!is.null(pt$rotate) && length(pt$rotate) > 0) {
        updateSelectizeInput(session, "mt_nodes_to_rotate",
                             selected = as.character(pt$rotate))
      } else {
        updateSelectizeInput(session, "mt_nodes_to_rotate",
                             selected = character(0))
      }
    }
  })

  # ================================================================
  # ROTATION SERVER LOGIC (mirrors single-mode rotation)
  # ================================================================

  # Populate node selector with tree node numbers when a tree is selected
  observe({
    req(input$mt_tree_selector, mt_values$trees)
    tn <- input$mt_tree_selector
    tree_obj <- mt_values$trees[[tn]]
    if (!is.null(tree_obj)) {
      all_nodes <- tree_obj$edge[, 1]
      unique_nodes <- sort(unique(all_nodes))
      updateSelectizeInput(session, "mt_nodes_to_rotate",
                           choices = as.character(unique_nodes),
                           server = TRUE)
    }
  })

  # Rotation classes UI (primary/secondary) - mirrors single mode
  output$mt_rotation_classes_ui <- renderUI({
    req(input$mt_enable_rotation, input$mt_rotation_type != "manual")
    if (input$mt_rotation_type == "primary") {
      tagList(
        numericInput("mt_rotation1_num_groups", "Number of rotation groups:", value = 2, min = 2, max = 10),
        uiOutput("mt_rotation1_groups_ui"),
        br(),
        actionButton("mt_apply_rotation1", "Apply Rotation", icon = icon("rotate"), class = "btn-primary"),
        actionButton("mt_clear_rotation1", "Clear Primary Configuration", icon = icon("trash"), class = "btn-warning")
      )
    } else if (input$mt_rotation_type == "secondary") {
      tagList(
        numericInput("mt_rotation2_num_groups", "Number of rotation groups:", value = 2, min = 2, max = 10),
        uiOutput("mt_rotation2_groups_ui"),
        br(),
        actionButton("mt_apply_rotation2", "Apply Rotation", icon = icon("rotate"), class = "btn-primary"),
        actionButton("mt_clear_rotation2", "Clear Secondary Configuration", icon = icon("trash"), class = "btn-warning")
      )
    }
  })

  # Rotation status box
  output$mt_rotation_status_box <- renderUI({
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    pt <- mt_values$per_tree[[tn]]
    if (is.null(pt)) return(NULL)

    status_items <- list()
    if (!is.null(pt$rotation1_config) && length(pt$rotation1_config) > 0) {
      status_items <- c(status_items, list(tags$li("Primary rotation configured")))
    }
    if (!is.null(pt$rotation2_config) && length(pt$rotation2_config) > 0) {
      status_items <- c(status_items, list(tags$li("Secondary rotation configured")))
    }
    if (!is.null(pt$rotate) && length(pt$rotate) > 0) {
      status_items <- c(status_items, list(tags$li(paste0("Manual rotation: nodes ", paste(pt$rotate, collapse = ", ")))))
    }
    if (length(status_items) == 0) return(NULL)

    tags$div(
      style = "margin-top: 10px; padding: 10px; background-color: #d4edda; border-radius: 5px;",
      tags$strong(paste0("Rotation status for: ", tn)),
      tags$ul(status_items)
    )
  })

  # rotation1 groups UI
  output$mt_rotation1_groups_ui <- renderUI({
    req(input$mt_rotation1_num_groups)
    csv_data <- mt_csv_for_classification()
    if (is.null(csv_data)) return(tags$p("Load CSV data first."))
    col_names <- names(csv_data)
    group_uis <- lapply(1:input$mt_rotation1_num_groups, function(i) {
      tagList(
        hr(),
        h5(paste("Rotation Group", i)),
        fluidRow(
          column(6, selectInput(inputId = paste0("mt_rotation1_col_", i), label = "Column:",
                                choices = c("-- Select --" = "", col_names))),
          column(6, uiOutput(paste0("mt_rotation1_val_ui_", i)))
        )
      )
    })
    tagList(group_uis)
  })

  # rotation2 groups UI
  output$mt_rotation2_groups_ui <- renderUI({
    req(input$mt_rotation2_num_groups)
    csv_data <- mt_csv_for_classification()
    if (is.null(csv_data)) return(tags$p("Load CSV data first."))
    col_names <- names(csv_data)
    group_uis <- lapply(1:input$mt_rotation2_num_groups, function(i) {
      tagList(
        hr(),
        h5(paste("Rotation Group", i)),
        fluidRow(
          column(6, selectInput(inputId = paste0("mt_rotation2_col_", i), label = "Column:",
                                choices = c("-- Select --" = "", col_names))),
          column(6, uiOutput(paste0("mt_rotation2_val_ui_", i)))
        )
      )
    })
    tagList(group_uis)
  })

  # Dynamic value dropdowns for rotation1 groups
  observe({
    req(input$mt_rotation1_num_groups)
    csv_data <- isolate(mt_csv_for_classification())
    if (is.null(csv_data)) return()
    for (i in 1:input$mt_rotation1_num_groups) {
      local({
        my_i <- i
        output[[paste0("mt_rotation1_val_ui_", my_i)]] <- renderUI({
          col_input <- input[[paste0("mt_rotation1_col_", my_i)]]
          if (is.null(col_input) || col_input == "") {
            selectInput(inputId = paste0("mt_rotation1_val_", my_i), label = "Value:",
                        choices = c("Select column first" = ""))
          } else {
            csv_to_use <- isolate(mt_csv_for_classification())
            unique_vals <- unique(csv_to_use[[col_input]])
            unique_vals <- unique_vals[!is.na(unique_vals)]
            selectInput(inputId = paste0("mt_rotation1_val_", my_i), label = "Value:",
                        choices = c("-- Select --" = "", unique_vals))
          }
        })
      })
    }
  })

  # Dynamic value dropdowns for rotation2 groups
  observe({
    req(input$mt_rotation2_num_groups)
    csv_data <- isolate(mt_csv_for_classification())
    if (is.null(csv_data)) return()
    for (i in 1:input$mt_rotation2_num_groups) {
      local({
        my_i <- i
        output[[paste0("mt_rotation2_val_ui_", my_i)]] <- renderUI({
          col_input <- input[[paste0("mt_rotation2_col_", my_i)]]
          if (is.null(col_input) || col_input == "") {
            selectInput(inputId = paste0("mt_rotation2_val_", my_i), label = "Value:",
                        choices = c("Select column first" = ""))
          } else {
            csv_to_use <- isolate(mt_csv_for_classification())
            unique_vals <- unique(csv_to_use[[col_input]])
            unique_vals <- unique_vals[!is.na(unique_vals)]
            selectInput(inputId = paste0("mt_rotation2_val_", my_i), label = "Value:",
                        choices = c("-- Select --" = "", unique_vals))
          }
        })
      })
    }
  })

  # Apply rotation1 (primary)
  observeEvent(input$mt_apply_rotation1, ignoreInit = TRUE, {
    req(input$mt_tree_selector, input$mt_rotation1_num_groups)
    tn <- input$mt_tree_selector
    config <- list()
    for (i in 1:input$mt_rotation1_num_groups) {
      col <- input[[paste0("mt_rotation1_col_", i)]]
      val <- input[[paste0("mt_rotation1_val_", i)]]
      if (!is.null(col) && col != "" && !is.null(val) && val != "") {
        config[[i]] <- list(col = col, val = val)
      }
    }
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$rotation1_config <- config
      mt_values$per_tree[[tn]]$rotation_type <- "primary"
      mt_values$per_tree[[tn]]$enable_rotation <- TRUE
    }
    showNotification(paste0("Primary rotation applied to ", tn), type = "message")
    mt_request_render(dirty = tn)
  })

  # Apply rotation2 (secondary)
  observeEvent(input$mt_apply_rotation2, ignoreInit = TRUE, {
    req(input$mt_tree_selector, input$mt_rotation2_num_groups)
    tn <- input$mt_tree_selector
    config <- list()
    for (i in 1:input$mt_rotation2_num_groups) {
      col <- input[[paste0("mt_rotation2_col_", i)]]
      val <- input[[paste0("mt_rotation2_val_", i)]]
      if (!is.null(col) && col != "" && !is.null(val) && val != "") {
        config[[i]] <- list(col = col, val = val)
      }
    }
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$rotation2_config <- config
      mt_values$per_tree[[tn]]$rotation_type <- "secondary"
      mt_values$per_tree[[tn]]$enable_rotation <- TRUE
    }
    showNotification(paste0("Secondary rotation applied to ", tn), type = "message")
    mt_request_render(dirty = tn)
  })

  # Apply manual rotation
  observeEvent(input$mt_apply_manual_rotation, ignoreInit = TRUE, {
    req(input$mt_tree_selector, input$mt_nodes_to_rotate)
    tn <- input$mt_tree_selector
    nodes <- as.numeric(input$mt_nodes_to_rotate)
    nodes <- nodes[!is.na(nodes)]
    if (length(nodes) == 0) {
      showNotification("Please select at least one node to rotate", type = "error", duration = 5)
      return()
    }
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$rotate <- nodes
      mt_values$per_tree[[tn]]$rotation_type <- "manual"
      mt_values$per_tree[[tn]]$enable_rotation <- TRUE
    }
    showNotification(paste0("Manual rotation applied to ", tn), type = "message")
    mt_request_render(dirty = tn)
  })

  # Clear rotation buttons
  observeEvent(input$mt_clear_manual_rotation, ignoreInit = TRUE, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$rotate <- NULL
    }
    updateSelectizeInput(session, "mt_nodes_to_rotate", selected = character(0))
    showNotification(paste0("Manual rotation cleared for ", tn), type = "warning")
  })

  observeEvent(input$mt_clear_rotation1, ignoreInit = TRUE, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$rotation1_config <- list()
    }
    updateNumericInput(session, "mt_rotation1_num_groups", value = 2)
    showNotification(paste0("Primary rotation cleared for ", tn), type = "warning")
  })

  observeEvent(input$mt_clear_rotation2, ignoreInit = TRUE, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$rotation2_config <- list()
    }
    updateNumericInput(session, "mt_rotation2_num_groups", value = 2)
    showNotification(paste0("Secondary rotation cleared for ", tn), type = "warning")
  })

  # --- Per-tree highlight UI (rendered dynamically) ---
  output$mt_per_tree_highlight_ui <- renderUI({
    req(input$mt_tree_selector, mt_values$per_tree)
    tn <- input$mt_tree_selector
    pt <- mt_values$per_tree[[tn]]
    tagList(
      tags$p(tags$strong(paste0("Tree: ", tn))),
      sliderInput("mt_pt_highlight_height", "Ellipse Height (this tree):",
                  min = 0.1, max = 10,
                  value = if (!is.null(pt$highlight_adjust_height)) pt$highlight_adjust_height else 1,
                  step = 0.1),
      sliderInput("mt_pt_highlight_width", "Ellipse Width (this tree):",
                  min = 0.1, max = 10,
                  value = if (!is.null(pt$highlight_adjust_width)) pt$highlight_adjust_width else 1.5,
                  step = 0.1)
    )
  })

  # Save per-tree highlight values
  observeEvent(input$mt_pt_highlight_height, ignoreInit = TRUE, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$highlight_adjust_height <- input$mt_pt_highlight_height
    }
  })

  observeEvent(input$mt_pt_highlight_width, ignoreInit = TRUE, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$highlight_adjust_width <- input$mt_pt_highlight_width
    }
  })

  # --- Highlight values observer: per-value color + transparency settings (matches single mode lines 14621-14689) ---
  observeEvent(input$mt_highlight_values, ignoreInit = TRUE, {
    req(input$mt_highlight_column)

    selected_values <- input$mt_highlight_values

    if (is.null(selected_values) || length(selected_values) == 0) {
      output$mt_highlight_values_settings_ui <- renderUI({
        tags$div(
          tags$p(tags$strong("No values selected")),
          tags$p("Select values from the dropdown to customize their highlight settings.")
        )
      })
      return(NULL)
    }

    output$mt_highlight_values_settings_ui <- renderUI({
      value_settings <- lapply(seq_along(selected_values), function(i) {
        value <- selected_values[i]
        default_color <- rainbow(length(selected_values))[i]
        default_transparency <- 0.5

        tagList(
          tags$hr(style = if (i == 1) "margin-top: 0;" else ""),
          fluidRow(
            column(12, tags$h5(tags$strong(as.character(value))))
          ),
          fluidRow(
            column(4,
                   colourpicker::colourInput(
                     inputId = paste0("mt_highlight_color_", i),
                     label = "Color",
                     value = default_color,
                     allowTransparent = FALSE
                   )
            ),
            column(4,
                   sliderInput(
                     inputId = paste0("mt_highlight_transparency_", i),
                     label = "Transparency",
                     min = 0.1, max = 1, value = default_transparency, step = 0.1
                   )
            ),
            column(4,
                   tags$div(
                     style = "margin-top: 25px;",
                     tags$div(
                       style = sprintf(
                         "width: 100%%; height: 40px; border-radius: 5px; border: 2px solid #999; background-color: %s; opacity: %s;",
                         default_color, default_transparency
                       ),
                       id = paste0("mt_highlight_preview_box_", i)
                     )
                   )
            )
          )
        )
      })

      tagList(
        tags$h4("Individual Value Settings:"),
        tags$p(style = "color: #666;", "Customize color and transparency for each highlighted value."),
        value_settings
      )
    })
  }, ignoreNULL = FALSE)

  # --- Update highlight preview boxes when color/transparency changes ---
  observe({
    req(input$mt_highlight_values)
    lapply(seq_along(input$mt_highlight_values), function(i) {
      color_input <- paste0("mt_highlight_color_", i)
      trans_input <- paste0("mt_highlight_transparency_", i)
      if (!is.null(input[[color_input]]) && !is.null(input[[trans_input]])) {
        shinyjs::runjs(sprintf(
          "$('#mt_highlight_preview_box_%d').css({'background-color': '%s', 'opacity': '%s'});",
          i, input[[color_input]], input[[trans_input]]
        ))
      }
    })
  })

  # --- Highlight offset quick-adjust buttons ---
  observeEvent(input$mt_offset_down_big, ignoreInit = TRUE, {
    cur <- if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0
    updateNumericInput(session, "mt_highlight_offset", value = cur - 1)
  })
  observeEvent(input$mt_offset_down_small, ignoreInit = TRUE, {
    cur <- if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0
    updateNumericInput(session, "mt_highlight_offset", value = cur - 0.1)
  })
  observeEvent(input$mt_offset_down_tiny, ignoreInit = TRUE, {
    cur <- if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0
    updateNumericInput(session, "mt_highlight_offset", value = cur - 0.01)
  })
  observeEvent(input$mt_offset_up_big, ignoreInit = TRUE, {
    cur <- if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0
    updateNumericInput(session, "mt_highlight_offset", value = cur + 1)
  })
  observeEvent(input$mt_offset_up_small, ignoreInit = TRUE, {
    cur <- if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0
    updateNumericInput(session, "mt_highlight_offset", value = cur + 0.1)
  })
  observeEvent(input$mt_offset_up_tiny, ignoreInit = TRUE, {
    cur <- if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0
    updateNumericInput(session, "mt_highlight_offset", value = cur + 0.01)
  })

  # --- Highlight Apply & Preview button ---
  observeEvent(input$mt_apply_highlight, ignoreInit = TRUE, {
    mt_request_render()
  })

  # --- Highlight Remove button ---
  observeEvent(input$mt_remove_highlight, ignoreInit = TRUE, {
    updateCheckboxInput(session, "mt_enable_highlight", value = FALSE)
    updateSelectizeInput(session, "mt_highlight_values", selected = character(0))
    mt_values$highlight_yaml <- NULL
    mt_request_render()
  })

  # --- Per-tree extra UI (title + bg color) ---
  output$mt_per_tree_extra_ui <- renderUI({
    req(input$mt_extra_tree_selector, mt_values$per_tree)
    tn <- input$mt_extra_tree_selector
    if (tn == "all") return(tags$p(class = "text-muted", "Select a specific tree to edit its title and background color."))
    pt <- mt_values$per_tree[[tn]]
    if (is.null(pt)) return(NULL)
    tagList(
      tags$p(tags$strong(paste0("Tree: ", tn))),
      textInput("mt_pt_title", "Tree Title:",
                value = if (!is.null(pt$title)) pt$title else tn),
      colourpicker::colourInput("mt_pt_bg_color", "Tree Background Color:",
                                value = if (!is.null(pt$bg_color)) pt$bg_color else "#FFFFFF")
    )
  })

  # Save per-tree extra values
  observeEvent(input$mt_pt_title, ignoreInit = TRUE, {
    req(input$mt_extra_tree_selector)
    tn <- input$mt_extra_tree_selector
    if (tn == "all") return()
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$title <- input$mt_pt_title
    }
  })

  observeEvent(input$mt_pt_bg_color, ignoreInit = TRUE, {
    req(input$mt_extra_tree_selector)
    tn <- input$mt_extra_tree_selector
    if (tn == "all") return()
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$bg_color <- input$mt_pt_bg_color
    }
  })

  # --- Build shared settings list from current UI inputs ---
  mt_gather_shared <- function() {
    list(
      individual_name = if (!is.null(input$mt_individual_name)) input$mt_individual_name else "Sample1",
      id_column = input$mt_id_column,
      individual_column = if (!is.null(input$mt_individual_column)) input$mt_individual_column else "",
      classification_column = input$mt_classification_column,
      classification_title = if (!is.null(input$mt_classification_title)) input$mt_classification_title else "Cell type",
      fdr_perc = if (!is.null(input$mt_fdr_perc)) input$mt_fdr_perc else 0.1,
      no_cluster_color = if (!is.null(input$mt_no_cluster_color)) (input$mt_no_cluster_color == "gray") else TRUE,
      classification_groups = mt_values$classification_groups,
      enable_bootstrap = isTRUE(input$mt_enable_bootstrap),
      bootstrap_format = if (!is.null(input$mt_bootstrap_format)) input$mt_bootstrap_format else "triangles",
      bootstrap_param = if (!is.null(input$mt_bootstrap_param)) input$mt_bootstrap_param else 1,
      bootstrap_label_size = if (!is.null(input$mt_bootstrap_label_size)) input$mt_bootstrap_label_size else 1.5,
      man_boot_x_offset = if (!is.null(input$mt_man_boot_x_offset)) input$mt_man_boot_x_offset else 0,
      enable_highlight = isTRUE(input$mt_enable_highlight),
      highlight_yaml = mt_values$highlight_yaml,
      tip_font_size = if (!is.null(input$mt_tip_font_size)) input$mt_tip_font_size else 3,
      edge_width = if (!is.null(input$mt_edge_width)) input$mt_edge_width else 1,
      display_node_numbers = isTRUE(input$mt_display_node_numbers),
      node_number_font_size = if (!is.null(input$mt_node_number_font_size)) input$mt_node_number_font_size else 3.5,
      tip_length = if (!is.null(input$mt_tip_length)) input$mt_tip_length else 0.05,
      ladderize = isTRUE(input$mt_ladderize),
      trim_tips = isTRUE(input$mt_trim_tips),
      legend_title_size = if (!is.null(input$mt_legend_title_size)) input$mt_legend_title_size else 12,
      legend_text_size = if (!is.null(input$mt_legend_text_size)) input$mt_legend_text_size else 10,
      output_width = if (!is.null(input$mt_output_width)) input$mt_output_width else 29.7,
      output_height = if (!is.null(input$mt_output_height)) input$mt_output_height else 42,
      output_units = if (!is.null(input$mt_output_units)) input$mt_output_units else "cm",
      output_format = if (!is.null(input$mt_output_format)) input$mt_output_format else "pdf",
      man_params = mt_values$man_params,
      legend_settings = list(
        position = if (!is.null(input$mt_legend_position)) input$mt_legend_position else "right",
        show_classification = isTRUE(input$mt_legend_show_classification),
        show_highlight = isTRUE(input$mt_legend_show_highlight),
        show_bootstrap = isTRUE(input$mt_legend_show_bootstrap),
        show_pvalue = isTRUE(input$mt_legend_show_pvalue),
        title_size = if (!is.null(input$mt_legend_title_size)) input$mt_legend_title_size else 12,
        text_size = if (!is.null(input$mt_legend_text_size)) input$mt_legend_text_size else 10,
        font_family = if (!is.null(input$mt_legend_font_family)) input$mt_legend_font_family else "sans",
        key_size = if (!is.null(input$mt_legend_key_size)) input$mt_legend_key_size else 1,
        key_width = if (!is.null(input$mt_legend_key_width)) input$mt_legend_key_width else 1,
        key_height = if (!is.null(input$mt_legend_key_height)) input$mt_legend_key_height else 1,
        spacing = if (!is.null(input$mt_legend_spacing)) input$mt_legend_spacing else 0.3,
        spacing_vertical = if (!is.null(input$mt_legend_spacing_vertical)) input$mt_legend_spacing_vertical else 1,
        title_key_spacing = if (!is.null(input$mt_legend_title_key_spacing)) input$mt_legend_title_key_spacing else 0.2,
        key_spacing = if (!is.null(input$mt_legend_key_spacing)) input$mt_legend_key_spacing else 0.1,
        reverse_order = isTRUE(input$mt_legend_reverse_order),
        box_background = if (!is.null(input$mt_legend_box_background)) input$mt_legend_box_background else "transparent",
        margin = if (!is.null(input$mt_legend_margin)) input$mt_legend_margin else 0.2,
        scope = if (!is.null(input$mt_legend_scope)) input$mt_legend_scope else "all",
        legend_tree = if (!is.null(input$mt_legend_tree_selector)) input$mt_legend_tree_selector else NULL
      )
    )
  }

  mt_last_render_time <- reactiveVal(0)

  # ==========================================================================
  # DEBOUNCED REACTIVE INPUTS (matches single-mode pattern, prevents rapid UI updates)
  # ==========================================================================
  MT_DEBOUNCE_MS <- 300

  # Tree display sliders
  mt_tip_font_size_d <- debounce(reactive(input$mt_tip_font_size), MT_DEBOUNCE_MS)
  mt_edge_width_d <- debounce(reactive(input$mt_edge_width), MT_DEBOUNCE_MS)
  mt_node_number_font_size_d <- debounce(reactive(input$mt_node_number_font_size), MT_DEBOUNCE_MS)
  mt_tip_length_d <- debounce(reactive(input$mt_tip_length), MT_DEBOUNCE_MS)

  # Bootstrap sliders
  mt_bootstrap_label_size_d <- debounce(reactive(input$mt_bootstrap_label_size), MT_DEBOUNCE_MS)
  mt_man_boot_x_offset_d <- debounce(reactive(input$mt_man_boot_x_offset), MT_DEBOUNCE_MS)

  # Highlight sliders
  mt_highlight_offset_d <- debounce(reactive(input$mt_highlight_offset), MT_DEBOUNCE_MS)
  mt_highlight_vertical_offset_d <- debounce(reactive(input$mt_highlight_vertical_offset), MT_DEBOUNCE_MS)
  mt_highlight_adjust_height_d <- debounce(reactive(input$mt_highlight_adjust_height), MT_DEBOUNCE_MS)
  mt_highlight_adjust_width_d <- debounce(reactive(input$mt_highlight_adjust_width), MT_DEBOUNCE_MS)

  # Legend sliders
  mt_legend_title_size_d <- debounce(reactive(input$mt_legend_title_size), MT_DEBOUNCE_MS)
  mt_legend_text_size_d <- debounce(reactive(input$mt_legend_text_size), MT_DEBOUNCE_MS)
  mt_legend_key_size_d <- debounce(reactive(input$mt_legend_key_size), MT_DEBOUNCE_MS)
  mt_legend_key_width_d <- debounce(reactive(input$mt_legend_key_width), MT_DEBOUNCE_MS)
  mt_legend_key_height_d <- debounce(reactive(input$mt_legend_key_height), MT_DEBOUNCE_MS)
  mt_legend_spacing_d <- debounce(reactive(input$mt_legend_spacing), MT_DEBOUNCE_MS)
  mt_legend_spacing_vertical_d <- debounce(reactive(input$mt_legend_spacing_vertical), MT_DEBOUNCE_MS)
  mt_legend_title_key_spacing_d <- debounce(reactive(input$mt_legend_title_key_spacing), MT_DEBOUNCE_MS)
  mt_legend_key_spacing_d <- debounce(reactive(input$mt_legend_key_spacing), MT_DEBOUNCE_MS)
  mt_legend_margin_d <- debounce(reactive(input$mt_legend_margin), MT_DEBOUNCE_MS)

  # Extra tab sliders
  mt_plot_offset_x_d <- debounce(reactive(input$mt_plot_offset_x), MT_DEBOUNCE_MS)
  mt_plot_offset_y_d <- debounce(reactive(input$mt_plot_offset_y), MT_DEBOUNCE_MS)
  mt_plot_scale_percent_d <- debounce(reactive(input$mt_plot_scale_percent), MT_DEBOUNCE_MS)
  mt_page_title_size_d <- debounce(reactive(input$mt_page_title_size), MT_DEBOUNCE_MS)
  mt_page_title_x_d <- debounce(reactive(input$mt_page_title_x), MT_DEBOUNCE_MS)
  mt_page_title_y_d <- debounce(reactive(input$mt_page_title_y), MT_DEBOUNCE_MS)
  mt_background_color_d <- debounce(reactive(input$mt_background_color), MT_DEBOUNCE_MS)
  # ==========================================================================
  # END DEBOUNCED REACTIVE INPUTS
  # ==========================================================================

  # --- Build highlight YAML (called from mt_do_render, already inside isolate) ---
  mt_build_highlight_yaml <- function() {
    if (!isTRUE(input$mt_enable_highlight) || is.null(input$mt_highlight_values) ||
        length(input$mt_highlight_values) == 0) {
      mt_values$highlight_yaml <- NULL
      return()
    }
    hl_vals <- input$mt_highlight_values
    hl_according <- list()
    for (j in seq_along(hl_vals)) {
      color_id <- paste0("mt_highlight_color_", j)
      trans_id <- paste0("mt_highlight_transparency_", j)
      val_color <- if (!is.null(input[[color_id]])) input[[color_id]] else rainbow(length(hl_vals))[j]
      val_trans <- if (!is.null(input[[trans_id]])) input[[trans_id]] else 0.5
      entry <- list()
      entry[[as.character(j)]] <- list(
        title1 = input$mt_highlight_column,
        value1 = list(hl_vals[j]),
        display_name = hl_vals[j],
        color = val_color,
        transparency = val_trans
      )
      hl_according[[j]] <- entry
    }
    mt_values$highlight_yaml <- list(
      display = "yes",
      title = if (!is.null(input$mt_highlight_title) && nchar(input$mt_highlight_title) > 0) input$mt_highlight_title else "Highlight",
      offset = if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0,
      vertical_offset = if (!is.null(input$mt_highlight_vertical_offset)) input$mt_highlight_vertical_offset else 0,
      adjust_height = if (!is.null(input$mt_highlight_adjust_height)) input$mt_highlight_adjust_height else 1,
      adjust_width = if (!is.null(input$mt_highlight_adjust_width)) input$mt_highlight_adjust_width else 1.5,
      according = hl_according
    )
  }

  # --- Core render function (called via session$onFlushed — NOT a reactive context) ---
  mt_do_render <- function() {
    isolate({
      if (is.null(mt_values$newick_paths) || is.null(mt_values$csv_path)) return()
      if (isTRUE(mt_values$plot_generating)) return()
      if (mt_classification_loading()) return()

      mt_values$plot_generating <- TRUE
      on.exit({
        mt_values$plot_generating <- FALSE
        mt_show_ready()
      })

      # Clean up old plot files (keep last 5)
      old_svgs <- list.files(tempdir(), pattern = paste0("^file.*\\.svg$"), full.names = TRUE)
      if (length(old_svgs) > 5) {
        old_sorted <- old_svgs[order(file.info(old_svgs)$mtime)]
        sapply(old_sorted[1:(length(old_sorted) - 5)], unlink)
      }

      mt_values$log_messages <- ""
      mt_log("Starting multi-tree render...")

      mt_save_current_tree_rotation()
      mt_build_highlight_yaml()

      shared <- mt_gather_shared()

      if (isTRUE(input$mt_a4_output)) {
        shared$output_width <- 297
        shared$output_height <- 210
        shared$output_units <- "mm"
      }

      tree_titles <- sapply(mt_values$newick_names, function(tn) {
        pt <- mt_values$per_tree[[tn]]
        if (!is.null(pt$title) && nchar(pt$title) > 0) pt$title else tn
      })
      tree_bg_colors <- sapply(mt_values$newick_names, function(tn) {
        pt <- mt_values$per_tree[[tn]]
        if (!is.null(pt$bg_color)) pt$bg_color else "#FFFFFF"
      })

      render_csv_path <- if (!is.null(mt_values$temp_csv_path) && file.exists(mt_values$temp_csv_path)) {
        mt_values$temp_csv_path
      } else {
        mt_values$csv_path
      }

      dirty <- mt_values$dirty_trees

      result <- tryCatch(
        func.multiple.trees.one.page.in.app(
          newick_paths = mt_values$newick_paths,
          newick_names = mt_values$newick_names,
          csv_path = render_csv_path,
          shared_settings = shared,
          per_tree_list = mt_values$per_tree,
          tree_titles = tree_titles,
          tree_bg_colors = tree_bg_colors,
          log_fn = mt_log,
          p_cache_per_tree = mt_values$p_cache_per_tree,
          cached_tree_plots = mt_values$cached_tree_plots,
          dirty_trees = dirty
        ),
        error = function(e) {
          mt_log(paste0("RENDER ERROR: ", e$message))
          NULL
        }
      )

      if (!is.null(result)) {
        cached_p <- attr(result, "p_cache_per_tree")
        if (!is.null(cached_p)) mt_values$p_cache_per_tree <- cached_p
        cached_plots <- attr(result, "cached_tree_plots")
        if (!is.null(cached_plots)) mt_values$cached_tree_plots <- cached_plots

        mt_values$last_plot <- result
        mt_values$dirty_trees <- character(0)
        mt_log("Render complete!")

        tmp <- tempfile(fileext = ".svg")
        w <- if (isTRUE(input$mt_a4_output)) 297 else {
          v <- input$mt_output_width; if (!is.null(v)) v else 29.7
        }
        h <- if (isTRUE(input$mt_a4_output)) 210 else {
          v <- input$mt_output_height; if (!is.null(v)) v else 42
        }
        u <- if (isTRUE(input$mt_a4_output)) "mm" else {
          v <- input$mt_output_units; if (!is.null(v)) v else "cm"
        }
        ggplot2::ggsave(tmp, plot = result, width = w, height = h, units = u,
                        limitsize = FALSE)
        mt_values$last_plot_file <- tmp
        mt_values$plot_counter <- mt_values$plot_counter + 1
      }

      mt_last_render_time(as.numeric(Sys.time()) * 1000)

      current_count <- mt_values$plot_counter
      if (current_count > 0 && current_count %% 3 == 0) {
        later::later(function() { gc(verbose = FALSE) }, delay = 0.1)
      }
    })
  }

  # Request a render: shows Processing, flushes to browser, THEN does heavy work
  mt_request_render <- function(dirty = NULL) {
    if (isTRUE(mt_values$plot_generating)) return()
    if (is.null(mt_values$newick_paths) || is.null(mt_values$csv_path)) return()

    if (!is.null(dirty)) {
      current_dirty <- isolate(mt_values$dirty_trees)
      if (is.null(current_dirty)) {
        mt_values$dirty_trees <- NULL
      } else {
        mt_values$dirty_trees <- union(current_dirty, dirty)
      }
    } else {
      mt_values$dirty_trees <- NULL
    }

    mt_show_processing()

    session$onFlushed(function() {
      mt_do_render()
    }, once = TRUE)
  }

  # --- Status indicator helpers ---
  mt_show_processing <- function() {
    for (pfx in c("mt_status", "mt_class_status", "mt_boot_status", "mt_high_status",
                   "mt_legend_status", "mt_extra_status", "mt_download_status")) {
      shinyjs::hide(paste0(pfx, "_waiting"))
      shinyjs::hide(paste0(pfx, "_ready"))
      shinyjs::show(paste0(pfx, "_processing"))
    }
  }

  mt_show_ready <- function() {
    for (pfx in c("mt_status", "mt_class_status", "mt_boot_status", "mt_high_status",
                   "mt_legend_status", "mt_extra_status", "mt_download_status")) {
      shinyjs::hide(paste0(pfx, "_waiting"))
      shinyjs::hide(paste0(pfx, "_processing"))
      shinyjs::show(paste0(pfx, "_ready"))
    }
  }

  # --- Plot outputs: all read cached SVG via plot_counter (no cascading renders) ---
  mt_plot_image <- function(alt_empty) {
    mt_values$plot_counter
    pf <- isolate(mt_values$last_plot_file)
    if (is.null(pf) || !file.exists(pf)) {
      return(list(src = "", alt = alt_empty))
    }
    list(src = pf, contentType = "image/svg+xml", width = "100%", alt = "Multi-tree plot")
  }

  output$mt_combined_plot <- renderImage({
    mt_plot_image("No plot yet. Upload data, process, and click Update Preview.")
  }, deleteFile = FALSE)

  output$mt_tree_display_preview <- renderImage({
    mt_plot_image("Click Apply & Preview to render.")
  }, deleteFile = FALSE)

  output$mt_bootstrap_preview <- renderImage({
    mt_plot_image("Click Apply & Preview to render.")
  }, deleteFile = FALSE)

  output$mt_highlight_preview <- renderImage({
    mt_plot_image("Click Apply & Preview to render.")
  }, deleteFile = FALSE)

  # --- Main Update Preview button (tree display tab) ---
  observeEvent(input$mt_update_preview, ignoreInit = TRUE, {
    mt_request_render()
  })

  # --- Tree display Apply & Preview button ---
  observeEvent(input$mt_apply_tree_display, ignoreInit = TRUE, {
    mt_request_render()
  })

  # --- Bootstrap Apply & Preview button ---
  observeEvent(input$mt_apply_bootstrap, ignoreInit = TRUE, {
    mt_request_render()
  })

  # --- Legend Apply button ---
  observeEvent(input$mt_apply_legend, ignoreInit = TRUE, {
    mt_request_render()
    showNotification("Legend settings applied", type = "message", duration = 2)
  })

  # --- Legend preview output ---
  output$mt_legend_preview <- renderImage({
    mt_plot_image("Click Apply Legend Settings to render.")
  }, deleteFile = FALSE)

  # Keep legend tree selector in sync
  observe({
    nn <- mt_values$newick_names
    if (!is.null(nn) && length(nn) > 0) {
      updateSelectInput(session, "mt_legend_tree_selector", choices = nn, selected = nn[1])
    }
  })

  # --- Log output ---
  output$mt_log <- renderText({
    mt_values$log_messages
  })

  # ================================================================
  # EXTRA TAB SERVER LOGIC
  # ================================================================

  # Keep extra tree selector in sync
  observe({
    nn <- mt_values$newick_names
    if (!is.null(nn) && length(nn) > 0) {
      updateSelectInput(session, "mt_extra_tree_selector",
                        choices = c("All Trees" = "all", setNames(nn, nn)),
                        selected = "all")
    }
  })

  # Extra preview output
  output$mt_extra_preview <- renderImage({
    mt_plot_image("Click Apply to Plot to render.")
  }, deleteFile = FALSE)

  # Extra Apply button
  observeEvent(input$mt_extra_apply, ignoreInit = TRUE, {
    mt_request_render()
  })

  # Reset position
  observeEvent(input$mt_reset_plot_position, ignoreInit = TRUE, {
    updateSliderInput(session, "mt_plot_offset_x", value = 0)
    updateSliderInput(session, "mt_plot_offset_y", value = 0)
  })

  # Reset scale
  observeEvent(input$mt_reset_plot_scale, ignoreInit = TRUE, {
    updateSliderInput(session, "mt_plot_scale_percent", value = 100)
  })

  # Reset background
  observeEvent(input$mt_reset_background, ignoreInit = TRUE, {
    colourpicker::updateColourInput(session, "mt_background_color", value = "#FFFFFF")
  })

  # Custom text annotations storage
  mt_custom_texts <- reactiveVal(list())

  observeEvent(input$mt_add_custom_text, ignoreInit = TRUE, {
    req(input$mt_custom_text_content)
    if (nchar(trimws(input$mt_custom_text_content)) > 0) {
      new_text <- list(
        content = input$mt_custom_text_content,
        x = if (!is.null(input$mt_custom_text_x)) input$mt_custom_text_x else 0.5,
        y = if (!is.null(input$mt_custom_text_y)) input$mt_custom_text_y else 0.5,
        size = if (!is.null(input$mt_custom_text_size)) input$mt_custom_text_size else 12,
        color = if (!is.null(input$mt_custom_text_color)) input$mt_custom_text_color else "#000000",
        fontface = if (!is.null(input$mt_custom_text_fontface)) input$mt_custom_text_fontface else "plain",
        hjust = if (!is.null(input$mt_custom_text_hjust)) as.numeric(input$mt_custom_text_hjust) else 0.5,
        vjust = if (!is.null(input$mt_custom_text_vjust)) as.numeric(input$mt_custom_text_vjust) else 0.5,
        angle = if (!is.null(input$mt_custom_text_angle)) input$mt_custom_text_angle else 0
      )
      mt_custom_texts(c(mt_custom_texts(), list(new_text)))
      updateTextInput(session, "mt_custom_text_content", value = "")
    }
  })

  observeEvent(input$mt_clear_custom_texts, ignoreInit = TRUE, {
    mt_custom_texts(list())
  })

  output$mt_custom_texts_list <- renderUI({
    texts <- mt_custom_texts()
    if (length(texts) == 0) {
      return(tags$p(class = "text-muted", "No custom texts added yet."))
    }
    text_items <- lapply(seq_along(texts), function(i) {
      txt <- texts[[i]]
      tags$div(
        style = "padding: 5px; margin: 5px 0; border: 1px solid #ddd; border-radius: 4px; background-color: #f9f9f9;",
        tags$span(style = paste0("color: ", txt$color, "; font-weight: ",
                                  if (txt$fontface %in% c("bold", "bold.italic")) "bold" else "normal", ";"),
                  paste0(i, ". \"", substr(txt$content, 1, 30), if (nchar(txt$content) > 30) "..." else "", "\"")),
        tags$span(class = "text-muted",
                  paste0(" (x:", round(txt$x, 2), ", y:", round(txt$y, 2), ", size:", txt$size, ")")),
        actionButton(paste0("mt_delete_text_", i), "", icon = icon("times"),
                     class = "btn-xs btn-danger", style = "float: right; padding: 2px 6px;")
      )
    })
    do.call(tagList, text_items)
  })

  observe({
    texts <- mt_custom_texts()
    lapply(seq_along(texts), function(i) {
      observeEvent(input[[paste0("mt_delete_text_", i)]], {
        current <- mt_custom_texts()
        if (i <= length(current)) {
          mt_custom_texts(current[-i])
        }
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # ================================================================
  # DOWNLOAD TAB SERVER LOGIC
  # ================================================================

  # Download preview output
  output$mt_download_preview <- renderImage({
    mt_plot_image("No plot yet.")
  }, deleteFile = FALSE)

  # Download plot handler
  output$mt_download_plot <- downloadHandler(
    filename = function() {
      if (isTRUE(input$mt_replace_name)) {
        paste0(input$mt_custom_name, ".", input$mt_output_format)
      } else {
        paste0(input$mt_prefix_text, input$mt_individual_name, input$mt_suffix_text, ".", input$mt_output_format)
      }
    },
    content = function(file) {
      plot_obj <- mt_values$last_plot
      if (is.null(plot_obj)) {
        showNotification("No plot to download. Click Update Preview first.", type = "error")
        return()
      }
      w <- if (!is.null(input$mt_output_width)) input$mt_output_width else 29.7
      h <- if (!is.null(input$mt_output_height)) input$mt_output_height else 42
      u <- if (!is.null(input$mt_output_units)) input$mt_output_units else "cm"
      tryCatch({
        ggplot2::ggsave(file, plot = plot_obj, width = w, height = h, units = u,
                        device = input$mt_output_format, limitsize = FALSE, dpi = 300)
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )

  # Download YAML handler
  output$mt_download_yaml <- downloadHandler(
    filename = function() {
      paste0(input$mt_individual_name, "_multi_tree_config.yaml")
    },
    content = function(file) {
      shared <- tryCatch(mt_gather_shared(), error = function(e) list())
      yaml_data <- tryCatch(
        mt_build_yaml_data(
          newick_path = if (!is.null(mt_values$newick_paths)) mt_values$newick_paths[1] else "tree.nwk",
          csv_path = if (!is.null(mt_values$csv_path)) mt_values$csv_path else "data.csv",
          shared = shared,
          per_tree = list()
        ),
        error = function(e) list()
      )
      yaml_text <- tryCatch(yaml::as.yaml(yaml_data, indent = 2), error = function(e) "# Error generating YAML")
      writeLines(yaml_text, file)
    }
  )

  # Config tab YAML download (same content, different button ID)
  output$mt_download_yaml_config <- downloadHandler(
    filename = function() { "multi_tree_config.yaml" },
    content = function(file) {
      shared <- tryCatch(mt_gather_shared(), error = function(e) list())
      yaml_data <- tryCatch(
        mt_build_yaml_data(
          newick_path = if (!is.null(mt_values$newick_paths)) mt_values$newick_paths[1] else "tree.nwk",
          csv_path = if (!is.null(mt_values$csv_path)) mt_values$csv_path else "data.csv",
          shared = shared, per_tree = list()
        ),
        error = function(e) list()
      )
      writeLines(tryCatch(yaml::as.yaml(yaml_data, indent = 2), error = function(e) "# Error"), file)
    }
  )

  # ================================================================
  # CLASSIFICATION TAB SERVER LOGIC
  # Mirrors single-mode classification workflow (lines 13907-18780)
  # ================================================================

  mt_classification_loading <- reactiveVal(FALSE)

  # Keep classification tree selector in sync with uploaded trees
  observe({
    nn <- mt_values$newick_names
    if (!is.null(nn) && length(nn) > 0) {
      updateSelectInput(session, "mt_classification_tree_selector",
                        choices = nn, selected = nn[1])
    }
  })

  # Helper: get the appropriate CSV for classification (filtered if available)
  mt_csv_for_classification <- function() {
    if (!is.null(mt_values$filtered_csv) && nrow(mt_values$filtered_csv) > 0) {
      mt_values$filtered_csv
    } else {
      mt_values$csv_data
    }
  }

  # --- Classification values UI: color pickers for each unique value ---
  # Mirrors single-mode output$classification_values_ui (lines 13907-14069)
  output$mt_classification_values_ui <- renderUI({
    req(input$mt_classification_column)
    mt_classification_loading(TRUE)
    on.exit(mt_classification_loading(FALSE))

    csv_data <- isolate(mt_csv_for_classification())
    column <- isolate(input$mt_classification_column)
    if (is.null(csv_data) || is.null(column) || column == "" || !(column %in% names(csv_data))) {
      return(tags$p(class = "text-muted", "Select a classification column. Process Data first for best performance."))
    }
    unique_values <- unique(csv_data[[column]])
    unique_values <- unique_values[!is.na(unique_values)]
    if (length(unique_values) > 100) {
      return(tags$div(
        style = "padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
        tags$h5("Too many unique values!", style = "color: #856404;"),
        tags$p(paste("Column has", length(unique_values), "values. Select a column with fewer (< 100)."))
      ))
    }
    if (length(unique_values) == 0) {
      return(tags$p(class = "text-muted", "No values found in selected column."))
    }
    r_colors <- c("Select color name..." = "",
      "red", "blue", "green", "yellow", "orange", "purple", "pink", "brown",
      "gray", "black", "white", "cyan", "magenta",
      "darkred", "darkblue", "darkgreen", "darkorange", "darkviolet",
      "lightblue", "lightgreen", "lightyellow", "lightpink", "lightgray",
      "steelblue", "skyblue", "navy", "maroon", "olive", "teal",
      "coral", "salmon", "tomato", "firebrick", "indianred",
      "forestgreen", "seagreen", "limegreen", "springgreen",
      "royalblue", "dodgerblue", "cornflowerblue", "midnightblue",
      "orchid", "plum", "violet", "mediumpurple", "darkorchid",
      "gold", "khaki", "goldenrod", "darkgoldenrod",
      "sienna", "chocolate", "peru", "tan", "wheat",
      "slategray", "dimgray", "darkgray", "lightslategray",
      "aquamarine", "turquoise", "cadetblue", "darkturquoise")

    palette_ui <- tags$div(
      style = "background-color: #f8f9fa; padding: 10px; margin-bottom: 15px; border-radius: 5px;",
      tags$h5("Quick Palette Options:", style = "margin-top: 0;"),
      fluidRow(
        column(6, selectInput("mt_color_palette_choice", "Apply Color Palette:",
                              choices = c("Rainbow" = "rainbow", "Heat Colors" = "heat.colors",
                                          "Terrain Colors" = "terrain.colors", "Topo Colors" = "topo.colors",
                                          "CM Colors" = "cm.colors",
                                          "Viridis" = "viridis", "Plasma" = "plasma",
                                          "Inferno" = "inferno", "Magma" = "magma",
                                          "Cividis" = "cividis"),
                              selected = "rainbow")),
        column(6, style = "padding-top: 25px;",
               actionButton("mt_apply_palette", "Apply Palette to All",
                            icon = icon("palette"), class = "btn-primary btn-sm"))
      )
    )
    default_colors <- rainbow(length(unique_values))
    class_values <- lapply(seq_along(unique_values), function(i) {
      fluidRow(
        column(4, tags$b(as.character(unique_values[i]))),
        column(4,
               colourpicker::colourInput(
                 inputId = paste0("mt_class_color_", i),
                 label = NULL,
                 value = default_colors[i],
                 allowTransparent = FALSE
               )
        ),
        column(4,
               selectInput(
                 inputId = paste0("mt_class_color_name_", i),
                 label = NULL,
                 choices = r_colors,
                 selected = ""
               )
        )
      )
    })
    tagList(palette_ui, tags$h5(paste0("Colors (", length(unique_values), " values):")), class_values)
  })

  # Apply palette button
  observeEvent(input$mt_apply_palette, ignoreInit = TRUE, {
    req(input$mt_classification_column, input$mt_color_palette_choice)
    csv_data <- mt_csv_for_classification()
    if (is.null(csv_data)) return()
    column <- input$mt_classification_column
    if (!(column %in% names(csv_data))) return()
    unique_values <- unique(csv_data[[column]])
    unique_values <- unique_values[!is.na(unique_values)]
    n <- length(unique_values)
    colors <- tryCatch(
      switch(input$mt_color_palette_choice,
             "rainbow" = rainbow(n), "heat.colors" = heat.colors(n),
             "terrain.colors" = terrain.colors(n), "topo.colors" = topo.colors(n),
             "cm.colors" = cm.colors(n),
             "viridis" = viridisLite::viridis(n), "plasma" = viridisLite::plasma(n),
             "inferno" = viridisLite::inferno(n), "magma" = viridisLite::magma(n),
             "cividis" = viridisLite::cividis(n),
             rainbow(n)),
      error = function(e) rainbow(n)
    )
    lapply(seq_along(unique_values), function(i) {
      colourpicker::updateColourInput(session, paste0("mt_class_color_", i), value = colors[i])
    })
    showNotification(paste("Applied", input$mt_color_palette_choice, "palette"), type = "message", duration = 2)
  })

  # R color name dropdown → colourInput sync: per-value observers (matches single mode pattern)
  # Single mode creates individual observeEvent per dropdown to avoid N-squared cascade
  observeEvent(input$mt_classification_column, {
    req(input$mt_classification_column)
    csv_data <- isolate(mt_csv_for_classification())
    if (is.null(csv_data)) return()
    column <- isolate(input$mt_classification_column)
    if (!(column %in% names(csv_data))) return()
    unique_values <- unique(csv_data[[column]])
    unique_values <- unique_values[!is.na(unique_values)]
    lapply(seq_along(unique_values), function(i) {
      observeEvent(input[[paste0("mt_class_color_name_", i)]], {
        color_name <- input[[paste0("mt_class_color_name_", i)]]
        if (!is.null(color_name) && color_name != "") {
          colourpicker::updateColourInput(session, paste0("mt_class_color_", i), value = color_name)
        }
      }, ignoreInit = TRUE)
    })
  }, ignoreInit = TRUE)

  # Helper: gather current classification from UI inputs
  mt_gather_classification <- function() {
    req(input$mt_classification_column)
    csv_data <- mt_csv_for_classification()
    if (is.null(csv_data)) return(NULL)
    column <- input$mt_classification_column
    if (!(column %in% names(csv_data))) return(NULL)
    unique_values <- unique(csv_data[[column]])
    unique_values <- unique_values[!is.na(unique_values)]
    classes <- lapply(seq_along(unique_values), function(i) {
      color_val <- input[[paste0("mt_class_color_", i)]]
      if (is.null(color_val)) color_val <- rainbow(length(unique_values))[i]
      list(
        column = column,
        value = unique_values[i],
        display_name = as.character(unique_values[i]),
        color = color_val
      )
    })
    list(
      title = if (!is.null(input$mt_classification_title)) input$mt_classification_title else "Cell type",
      column = column,
      classes = classes,
      fdr = if (!is.null(input$mt_fdr_perc)) input$mt_fdr_perc else 0.1,
      no_cluster_title = "No cluster",
      no_cluster_color = if (!is.null(input$mt_no_cluster_color)) input$mt_no_cluster_color else "gray",
      highlight = list(enabled = FALSE)
    )
  }

  # Reset Classification: clear all classification state
  observeEvent(input$mt_remove_classification, ignoreInit = TRUE, {
    mt_values$classification_groups <- list()
    mt_values$temp_classification_preview <- NULL
    # Clear per-tree overrides
    for (tn in names(mt_values$per_tree)) {
      mt_values$per_tree[[tn]]$classification_groups <- NULL
      mt_values$per_tree[[tn]]$classification_title <- NULL
      mt_values$per_tree[[tn]]$no_cluster_color <- NULL
    }
    showNotification("Classification reset", type = "warning")
  })

  # Save for This Tree (without rendering) - speeds up per-tree workflow
  observeEvent(input$mt_save_classification_no_render, ignoreInit = TRUE, {
    cls <- mt_gather_classification()
    if (is.null(cls)) return()
    tn <- input$mt_classification_tree_selector
    if (is.null(tn) || is.null(mt_values$per_tree[[tn]])) return()
    grps <- lapply(cls$classes, function(entry) {
      list(values = list(entry$value), display_name = entry$display_name, color = entry$color)
    })
    mt_values$per_tree[[tn]]$classification_groups <- grps
    mt_values$per_tree[[tn]]$classification_title <- cls$title
    mt_values$per_tree[[tn]]$no_cluster_color <- cls$no_cluster_color
    showNotification(paste0("Classification saved for: ", tn, " (click Apply & Preview to render)"),
                     type = "message", duration = 3)
  })

  # Apply & Preview: save classification (respecting scope) and trigger main render
  observeEvent(input$mt_update_classification_preview, ignoreInit = TRUE, {
    cls <- mt_gather_classification()
    if (is.null(cls)) return()

    scope <- if (!is.null(input$mt_classification_scope)) input$mt_classification_scope else "all"
    grps <- lapply(cls$classes, function(entry) {
      list(values = list(entry$value), display_name = entry$display_name, color = entry$color)
    })

    dirty <- NULL
    if (scope == "per_tree" && !is.null(input$mt_classification_tree_selector)) {
      tn <- input$mt_classification_tree_selector
      if (!is.null(mt_values$per_tree[[tn]])) {
        mt_values$per_tree[[tn]]$classification_groups <- grps
        mt_values$per_tree[[tn]]$classification_title <- cls$title
        mt_values$per_tree[[tn]]$no_cluster_color <- cls$no_cluster_color
      }
      dirty <- tn
      showNotification(paste0("Classification applied to tree: ", tn), type = "message", duration = 3)
    } else {
      mt_values$classification_groups <- grps
      mt_values$temp_classification_preview <- cls
      showNotification("Classification applied to all trees", type = "message", duration = 3)
    }

    mt_request_render(dirty = dirty)
  })

  output$mt_combined_plot_class <- renderImage({
    mt_plot_image("Click Apply & Preview to render.")
  }, deleteFile = FALSE)

  # --- Configuration YAML output ---
  output$mt_yaml_output <- renderText({
    mt_values$plot_counter
    isolate({
      shared <- tryCatch(mt_gather_shared(), error = function(e) list())
      yaml_data <- tryCatch(
        mt_build_yaml_data(
          newick_path = if (!is.null(mt_values$newick_paths)) mt_values$newick_paths[1] else "tree.nwk",
          csv_path = if (!is.null(mt_values$csv_path)) mt_values$csv_path else "data.csv",
          shared = shared,
          per_tree = list()
        ),
        error = function(e) list(error = e$message)
      )
      yaml::as.yaml(yaml_data, indent.mapping.sequence = TRUE)
    })
  })

  # --- YAML import observer (option ii: imports all shared settings) ---
  observeEvent(input$mt_yaml_config, ignoreInit = TRUE, {
    req(input$mt_yaml_config)
    mt_log("Importing YAML configuration...")

    yaml_data <- tryCatch(
      yaml::read_yaml(input$mt_yaml_config$datapath),
      error = function(e) {
        mt_log(paste0("ERROR reading YAML: ", e$message))
        NULL
      }
    )
    if (is.null(yaml_data)) return()

    # Import shared settings from YAML
    vis <- yaml_data[["visual definitions"]]
    if (!is.null(vis)) {
      # man_* params
      man_names <- grep("^man_", names(vis), value = TRUE)
      for (mn in man_names) {
        mt_values$man_params[[mn]] <- vis[[mn]]
      }
      if (length(man_names) > 0) {
        mt_log(paste0("Imported ", length(man_names), " man_* parameters"))
      }

      # Font sizes
      if (!is.null(vis$font_size)) {
        if (!is.null(vis$font_size$tips)) {
          updateSliderInput(session, "mt_tip_font_size", value = vis$font_size$tips)
        }
        if (!is.null(vis$font_size$legend_title)) {
          updateSliderInput(session, "mt_legend_title_size", value = vis$font_size$legend_title)
        }
        if (!is.null(vis$font_size$legend_text)) {
          updateSliderInput(session, "mt_legend_text_size", value = vis$font_size$legend_text)
        }
      }

      # Bootstrap
      if (!is.null(vis$Bootstrap)) {
        if (vis$Bootstrap$display == "yes") {
          updateCheckboxInput(session, "mt_enable_bootstrap", value = TRUE)
        }
        if (!is.null(vis$Bootstrap$format)) {
          updateSelectInput(session, "mt_bootstrap_format", selected = vis$Bootstrap$format)
        }
      }

      # Edge width
      if (!is.null(vis$edge_width_multiplier$size)) {
        updateSliderInput(session, "mt_edge_width", value = vis$edge_width_multiplier$size)
      }

      # Trim tips
      if (!is.null(vis$`trim tips`)) {
        updateCheckboxInput(session, "mt_trim_tips",
                            value = (vis$`trim tips`$display == "yes"))
        if (!is.null(vis$`trim tips`$length)) {
          updateSliderInput(session, "mt_tip_length", value = vis$`trim tips`$length)
        }
      }

      # Skip heatmap/SNP blocks
      if (!is.null(vis$classification)) {
        for (ci in seq_along(vis$classification)) {
          cls <- vis$classification[[ci]]
          cls_key <- names(cls)[1]
          if (!is.null(cls[[cls_key]]$heatmap_display)) {
            mt_log("NOTE: Heatmap blocks in YAML ignored (not available in multi-tree mode)")
          }
        }
      }
    }

    # Import individual settings
    gen <- yaml_data[["Individual general definitions"]]
    if (!is.null(gen)) {
      if (!is.null(gen$Individual)) {
        updateTextInput(session, "mt_individual_name", value = gen$Individual)
      }
    }

    mt_log("YAML import complete")
    showNotification("YAML configuration imported!", type = "message")
  })

} # End of mt_install_server
