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
        selectInput("mt_id_column", "Select ID Column", choices = NULL),
        selectInput("mt_individual_column", "Select Individual Column",
                    choices = NULL),
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
  tabItem(tabName = "mt_tree_display", fluidRow(
    box(title = "Tree Display Settings", status = "primary",
        solidHeader = TRUE, width = 4,
        tags$p(class = "text-muted", "These settings apply to all trees."),
        sliderInput("mt_tip_font_size", "Tip Font Size:",
                    min = 0.5, max = 20, value = 3, step = 0.5),
        sliderInput("mt_edge_width", "Edge Width:",
                    min = 0.1, max = 5, value = 1, step = 0.1),
        checkboxInput("mt_display_node_numbers", "Display Node Numbers", value = FALSE),
        conditionalPanel(
          condition = "input.mt_display_node_numbers == true",
          sliderInput("mt_node_number_font_size", "Node Number Font Size:",
                      min = 0.5, max = 10, value = 3.5, step = 0.5)
        ),
        sliderInput("mt_tip_length", "Tip Length:",
                    min = 0, max = 1, value = 0.05, step = 0.01),
        checkboxInput("mt_ladderize", "Ladderize Tree", value = TRUE),
        checkboxInput("mt_trim_tips", "Trim Tips", value = TRUE),
        tags$hr(),
        tags$h5("Per-Tree: Rotation"),
        selectInput("mt_tree_selector", "Select Tree:", choices = NULL),
        textInput("mt_nodes_to_rotate", "Nodes to Rotate (comma-separated):",
                  value = "", placeholder = "e.g. 15, 23, 42")
    )
  ))
}

mt_tabItem_classification <- function() {
  tabItem(tabName = "mt_classification", fluidRow(
    box(title = "Classification Settings", status = "primary",
        solidHeader = TRUE, width = 12,
        tags$p(class = "text-muted", "Classification applies identically to all trees."),
        selectInput("mt_classification_column", "Classification Column:", choices = NULL),
        numericInput("mt_fdr_perc", "FDR Percentage:", value = 0.1, min = 0, max = 1, step = 0.01),
        textInput("mt_classification_title", "Classification Title:", value = "Cell type"),
        checkboxInput("mt_no_cluster_color", "Gray for unclassified", value = TRUE),
        tags$hr(),
        uiOutput("mt_classification_preview")
    )
  ))
}

mt_tabItem_bootstrap <- function() {
  tabItem(tabName = "mt_bootstrap", fluidRow(
    box(title = "Bootstrap Settings", status = "primary",
        solidHeader = TRUE, width = 12,
        checkboxInput("mt_enable_bootstrap", "Enable Bootstrap Display", value = TRUE),
        conditionalPanel(
          condition = "input.mt_enable_bootstrap == true",
          selectInput("mt_bootstrap_format", "Bootstrap Format:",
                      choices = c("triangles", "circles"), selected = "triangles"),
          numericInput("mt_bootstrap_param", "Bootstrap Parameter:", value = 1, min = 0, max = 100),
          sliderInput("mt_bootstrap_label_size", "Bootstrap Label Size:",
                      min = 0.5, max = 10, value = 1.5, step = 0.5),
          sliderInput("mt_man_boot_x_offset", "Bootstrap X Offset:",
                      min = -0.1, max = 0.1, value = 0, step = 0.001)
        )
    )
  ))
}

mt_tabItem_highlighting <- function() {
  tabItem(tabName = "mt_highlighting", fluidRow(
    box(title = "Highlighting Settings", status = "primary",
        solidHeader = TRUE, width = 12,
        checkboxInput("mt_enable_highlight", "Enable Highlighting", value = FALSE),
        conditionalPanel(
          condition = "input.mt_enable_highlight == true",
          selectInput("mt_highlight_column", "Highlight Column:", choices = NULL),
          selectizeInput("mt_highlight_values", "Values to Highlight:",
                         choices = NULL, multiple = TRUE),
          textInput("mt_highlight_title", "Highlight Title:", value = "Highlight"),
          numericInput("mt_highlight_offset", "Highlight Offset:", value = 0, step = 0.1),
          sliderInput("mt_highlight_vertical_offset", "Vertical Offset:",
                      min = -5, max = 5, value = 0, step = 0.1),
          sliderInput("mt_highlight_adjust_height", "Ellipse Height:",
                      min = 0.1, max = 10, value = 1, step = 0.1),
          sliderInput("mt_highlight_adjust_width", "Ellipse Width:",
                      min = 0.1, max = 10, value = 1.5, step = 0.1),
          tags$hr(),
          tags$h5("Per-Tree Ellipse Adjustments"),
          uiOutput("mt_per_tree_highlight_ui")
        )
    )
  ))
}

mt_tabItem_legend <- function() {
  tabItem(tabName = "mt_legend", fluidRow(
    box(title = "Legend Settings", status = "primary",
        solidHeader = TRUE, width = 12,
        selectInput("mt_legend_position", "Legend Position:",
                    choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
        checkboxInput("mt_legend_show_classification", "Show Classification Legend", value = TRUE),
        checkboxInput("mt_legend_show_highlight", "Show Highlight Legend", value = TRUE),
        checkboxInput("mt_legend_show_bootstrap", "Show Bootstrap Legend", value = TRUE),
        tags$hr(),
        sliderInput("mt_legend_title_size", "Title Size:", min = 4, max = 30, value = 12, step = 1),
        sliderInput("mt_legend_text_size", "Text Size:", min = 4, max = 30, value = 10, step = 1),
        selectInput("mt_legend_font_family", "Font Family:",
                    choices = c("sans", "serif", "mono"), selected = "sans"),
        sliderInput("mt_legend_key_size", "Key Size:", min = 0.1, max = 3, value = 1, step = 0.1),
        sliderInput("mt_legend_spacing", "Spacing:", min = 0, max = 2, value = 0.3, step = 0.05),
        checkboxInput("mt_legend_reverse_order", "Reverse Legend Order", value = FALSE)
    )
  ))
}

mt_tabItem_extra <- function() {
  tabItem(tabName = "mt_extra", fluidRow(
    box(title = "Extra Settings", status = "primary",
        solidHeader = TRUE, width = 12,
        tags$h5("Page Layout"),
        checkboxInput("mt_a4_output", "A4 Page Format (297x210mm)", value = FALSE),
        colourpicker::colourInput("mt_background_color", "Page Background:", value = "#FFFFFF"),
        tags$hr(),
        tags$h5("Per-Tree Settings"),
        tags$p(class = "text-muted", "Select tree from Tree Display tab's tree selector."),
        uiOutput("mt_per_tree_extra_ui"),
        tags$hr(),
        tags$h5("Page Title"),
        checkboxInput("mt_enable_page_title", "Enable Page Title", value = FALSE),
        conditionalPanel(
          condition = "input.mt_enable_page_title == true",
          textInput("mt_page_title_text", "Title Text:", value = ""),
          sliderInput("mt_page_title_size", "Title Size:", min = 5, max = 30, value = 14, step = 1)
        )
    )
  ))
}

mt_tabItem_config <- function() {
  tabItem(tabName = "mt_config", fluidRow(
    box(title = "Configuration", status = "primary",
        solidHeader = TRUE, width = 12,
        tags$p("Current multi-tree YAML configuration:"),
        verbatimTextOutput("mt_yaml_output"),
        downloadButton("mt_download_yaml", "Download Configuration", class = "btn-success")
    )
  ))
}

mt_tabItem_download <- function() {
  tabItem(tabName = "mt_download", fluidRow(
    box(title = "Download Multi-Tree Plot", status = "primary",
        solidHeader = TRUE, width = 12,
        textInput("mt_individual_name", "Individual Name:", value = "Sample1"),
        selectInput("mt_output_format", "Output Format:",
                    choices = c("pdf", "png", "svg", "tiff"), selected = "pdf"),
        numericInput("mt_output_width", "Width:", value = 29.7, min = 1, max = 200),
        numericInput("mt_output_height", "Height:", value = 42, min = 1, max = 200),
        selectInput("mt_output_units", "Units:",
                    choices = c("cm", "mm", "in"), selected = "cm"),
        tags$hr(),
        checkboxInput("mt_replace_name", "Use Custom Filename", value = FALSE),
        conditionalPanel(
          condition = "input.mt_replace_name == true",
          textInput("mt_custom_name", "Custom Filename:", value = "multi_tree_plot")
        ),
        textInput("mt_prefix_text", "Prefix:", value = "MultiTree__"),
        textInput("mt_suffix_text", "Suffix:", value = ""),
        tags$hr(),
        downloadButton("mt_download_plot", "Download Plot", class = "btn-primary btn-block")
    )
  ))
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

  # Rotation
  rotation_nodes <- list()
  if (!is.null(per_tree$rotate) && length(per_tree$rotate) > 0) {
    for (k in seq_along(per_tree$rotate)) {
      node_entry <- list()
      node_entry[[as.character(k)]] <- per_tree$rotate[k]
      rotation_nodes[[k]] <- node_entry
    }
  }
  manual_rotation <- list(
    display = if (length(rotation_nodes) > 0) "yes" else "no",
    nodes = rotation_nodes
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
        "file_type" = shared$output_format,
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
      rotation1 = list(display = "no", according = list()),
      rotation2 = list(display = "no", according = list()),
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
      compare_two_trees = "no",
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
func.multiple.trees.one.page.in.app <- function(
  newick_paths, newick_names, csv_path,
  shared_settings, per_tree_list,
  tree_titles, tree_bg_colors, log_fn = function(msg) {}
) {
  list_out <- list()

  for (i in seq_along(newick_paths)) {
    tn <- newick_names[i]
    log_fn(paste0("Rendering tree ", i, "/", length(newick_paths), ": ", tn))

    # Build YAML for this tree
    per_tree <- if (!is.null(per_tree_list[[tn]])) per_tree_list[[tn]] else list()
    # Override Individual name per-tree if set (each tree can match a different CSV individual)
    tree_shared <- shared_settings
    if (!is.null(per_tree$individual_value) && nzchar(per_tree$individual_value)) {
      tree_shared$individual_name <- per_tree$individual_value
    }
    yaml_data <- mt_build_yaml_data(
      newick_path = newick_paths[i],
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
        legend_settings = shared_settings$legend_settings
      ))
    }, error = function(e) {
      log_fn(paste0("ERROR rendering tree ", tn, ": ", e$message))
      NULL
    })

    if (is.null(treesi)) next

    # Return-value fix: current func.print.lineage.tree returns
    # list(plots=..., rotated_tree=..., cache_data=...)
    one_plot <- NULL
    if (!is.null(treesi$plots) && length(treesi$plots) > 0) {
      one_plot <- treesi$plots[[1]]
    } else if (is.list(treesi) && length(treesi) > 0) {
      one_plot <- treesi[[1]]
    }
    if (is.null(one_plot)) {
      log_fn(paste0("WARNING: No plot extracted for tree ", tn))
      next
    }

    # Apply per-tree styling
    list_out[[length(list_out) + 1]] <- one_plot +
      ggtree::theme_tree(bgcolor = tree_bg_colors[i]) +
      ggplot2::ggtitle(tree_titles[i]) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
        plot.background = ggplot2::element_rect(
          fill = tree_bg_colors[i], color = tree_bg_colors[i]
        )
      ) +
      ggplot2::guides(color = "none", size = "none", shape = "none")
  }

  if (length(list_out) == 0) return(NULL)

  log_fn(paste0("Combining ", length(list_out), " trees into one page"))
  do.call(gridExtra::grid.arrange, c(list_out, ncol = 1))
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
    per_tree = list(),          # list[[tree_name]] = list(rotate, highlight_adjust_height, highlight_adjust_width, title, bg_color)
    man_params = list(),        # man_* params from YAML import
    highlight_yaml = NULL,      # shared highlight YAML block
    classification_groups = list(), # auto-generated classification groups
    last_plot = NULL,           # last rendered gtable
    log_messages = ""           # log text for mt_log output
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
  observeEvent(input$mt_add_upload_slot, {
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
  observeEvent(input$mt_csv_file, {
    req(input$mt_csv_file)
    mt_values$csv_path <- input$mt_csv_file$datapath
    mt_values$csv_data <- tryCatch(
      read.csv(input$mt_csv_file$datapath, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) {
        mt_log(paste0("ERROR reading CSV: ", e$message))
        NULL
      }
    )
    if (!is.null(mt_values$csv_data)) {
      cols <- names(mt_values$csv_data)
      mt_log(paste0("CSV loaded: ", nrow(mt_values$csv_data), " rows, ",
                    length(cols), " columns"))
      updateSelectInput(session, "mt_id_column", choices = cols, selected = cols[1])
      updateSelectInput(session, "mt_individual_column", choices = c("", cols))
      updateSelectInput(session, "mt_classification_column", choices = cols)
      updateSelectInput(session, "mt_highlight_column", choices = cols)
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
  observeEvent(input$mt_auto_match_individuals, {
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

  # --- Highlight column observer ---
  observeEvent(input$mt_highlight_column, {
    req(input$mt_highlight_column, mt_values$csv_data)
    vals <- sort(unique(na.omit(mt_values$csv_data[[input$mt_highlight_column]])))
    updateSelectizeInput(session, "mt_highlight_values",
                         choices = vals, selected = NULL)
  })

  # --- Process & Match IDs ---
  observeEvent(input$mt_process_data, {
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

  observeEvent(input$mt_tree_selector, {
    req(input$mt_tree_selector)

    # Save previous tree's per-tree values
    prev <- mt_prev_tree()
    if (!is.null(prev) && !is.null(mt_values$per_tree[[prev]])) {
      # Parse rotation nodes
      rot_text <- input$mt_nodes_to_rotate
      rot_nodes <- NULL
      if (!is.null(rot_text) && nchar(trimws(rot_text)) > 0) {
        rot_nodes <- suppressWarnings(
          as.numeric(trimws(unlist(strsplit(rot_text, ","))))
        )
        rot_nodes <- rot_nodes[!is.na(rot_nodes)]
      }
      mt_values$per_tree[[prev]]$rotate <- rot_nodes
    }

    # Load new tree's per-tree values
    current <- input$mt_tree_selector
    mt_prev_tree(current)

    pt <- mt_values$per_tree[[current]]
    if (!is.null(pt)) {
      rot_str <- if (!is.null(pt$rotate) && length(pt$rotate) > 0) {
        paste(pt$rotate, collapse = ", ")
      } else ""
      updateTextInput(session, "mt_nodes_to_rotate", value = rot_str)
    }
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
  observeEvent(input$mt_pt_highlight_height, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$highlight_adjust_height <- input$mt_pt_highlight_height
    }
  })

  observeEvent(input$mt_pt_highlight_width, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$highlight_adjust_width <- input$mt_pt_highlight_width
    }
  })

  # --- Per-tree extra UI (title + bg color) ---
  output$mt_per_tree_extra_ui <- renderUI({
    req(input$mt_tree_selector, mt_values$per_tree)
    tn <- input$mt_tree_selector
    pt <- mt_values$per_tree[[tn]]
    tagList(
      tags$p(tags$strong(paste0("Tree: ", tn))),
      textInput("mt_pt_title", "Tree Title:",
                value = if (!is.null(pt$title)) pt$title else tn),
      colourpicker::colourInput("mt_pt_bg_color", "Tree Background Color:",
                                value = if (!is.null(pt$bg_color)) pt$bg_color else "#FFFFFF")
    )
  })

  # Save per-tree extra values
  observeEvent(input$mt_pt_title, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
    if (!is.null(mt_values$per_tree[[tn]])) {
      mt_values$per_tree[[tn]]$title <- input$mt_pt_title
    }
  })

  observeEvent(input$mt_pt_bg_color, {
    req(input$mt_tree_selector)
    tn <- input$mt_tree_selector
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
      no_cluster_color = isTRUE(input$mt_no_cluster_color),
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
        title_size = if (!is.null(input$mt_legend_title_size)) input$mt_legend_title_size else 12,
        text_size = if (!is.null(input$mt_legend_text_size)) input$mt_legend_text_size else 10,
        font_family = if (!is.null(input$mt_legend_font_family)) input$mt_legend_font_family else "sans",
        key_size = if (!is.null(input$mt_legend_key_size)) input$mt_legend_key_size else 1,
        spacing = if (!is.null(input$mt_legend_spacing)) input$mt_legend_spacing else 0.3,
        reverse_order = isTRUE(input$mt_legend_reverse_order)
      )
    )
  }

  # --- Build highlight YAML when user configures highlighting ---
  observeEvent(list(input$mt_highlight_column, input$mt_highlight_values,
                    input$mt_highlight_title, input$mt_highlight_offset,
                    input$mt_highlight_vertical_offset,
                    input$mt_highlight_adjust_height, input$mt_highlight_adjust_width), {
    if (!isTRUE(input$mt_enable_highlight) || is.null(input$mt_highlight_values) ||
        length(input$mt_highlight_values) == 0) {
      mt_values$highlight_yaml <- NULL
      return()
    }
    # Build according list for highlight
    hl_according <- list()
    for (j in seq_along(input$mt_highlight_values)) {
      entry <- list()
      entry[[as.character(j)]] <- list(
        title1 = input$mt_highlight_column,
        value1 = list(input$mt_highlight_values[j]),
        display_name = input$mt_highlight_values[j],
        color = "yellow",
        transparency = 0.5
      )
      hl_according[[j]] <- entry
    }
    mt_values$highlight_yaml <- list(
      display = "yes",
      offset = if (!is.null(input$mt_highlight_offset)) input$mt_highlight_offset else 0,
      vertical_offset = if (!is.null(input$mt_highlight_vertical_offset)) input$mt_highlight_vertical_offset else 0,
      adjust_height = if (!is.null(input$mt_highlight_adjust_height)) input$mt_highlight_adjust_height else 1,
      adjust_width = if (!is.null(input$mt_highlight_adjust_width)) input$mt_highlight_adjust_width else 1.5,
      according = hl_according
    )
  }, ignoreNULL = FALSE)

  # --- Render: gated by Update Preview button ---
  mt_render_plot <- eventReactive(input$mt_update_preview, {
    req(mt_values$newick_paths, mt_values$csv_path)
    mt_values$log_messages <- ""  # clear log
    mt_log("Starting multi-tree render...")

    # Save current tree selector state before rendering
    if (!is.null(input$mt_tree_selector) && !is.null(mt_values$per_tree[[input$mt_tree_selector]])) {
      rot_text <- input$mt_nodes_to_rotate
      rot_nodes <- NULL
      if (!is.null(rot_text) && nchar(trimws(rot_text)) > 0) {
        rot_nodes <- suppressWarnings(as.numeric(trimws(unlist(strsplit(rot_text, ",")))))
        rot_nodes <- rot_nodes[!is.na(rot_nodes)]
      }
      mt_values$per_tree[[input$mt_tree_selector]]$rotate <- rot_nodes
    }

    shared <- mt_gather_shared()

    # A4 override
    if (isTRUE(input$mt_a4_output)) {
      shared$output_width <- 297
      shared$output_height <- 210
      shared$output_units <- "mm"
    }

    # Collect per-tree titles and bg colors
    tree_titles <- sapply(mt_values$newick_names, function(tn) {
      pt <- mt_values$per_tree[[tn]]
      if (!is.null(pt$title) && nchar(pt$title) > 0) pt$title else tn
    })
    tree_bg_colors <- sapply(mt_values$newick_names, function(tn) {
      pt <- mt_values$per_tree[[tn]]
      if (!is.null(pt$bg_color)) pt$bg_color else "#FFFFFF"
    })

    result <- tryCatch(
      func.multiple.trees.one.page.in.app(
        newick_paths = mt_values$newick_paths,
        newick_names = mt_values$newick_names,
        csv_path = mt_values$csv_path,
        shared_settings = shared,
        per_tree_list = mt_values$per_tree,
        tree_titles = tree_titles,
        tree_bg_colors = tree_bg_colors,
        log_fn = mt_log
      ),
      error = function(e) {
        mt_log(paste0("RENDER ERROR: ", e$message))
        NULL
      }
    )

    if (!is.null(result)) {
      mt_values$last_plot <- result
      mt_log("Render complete!")
    }
    result
  })

  # --- Plot output (renderImage to match imageOutput) ---
  output$mt_combined_plot <- renderImage({
    result <- mt_render_plot()
    if (is.null(result)) {
      return(list(src = "", alt = "No plot yet. Upload data, process, and click Update Preview."))
    }
    tmp <- tempfile(fileext = ".png")
    w <- if (isTRUE(input$mt_a4_output)) 297 else {
      if (!is.null(input$mt_output_width)) input$mt_output_width else 29.7
    }
    h <- if (isTRUE(input$mt_a4_output)) 210 else {
      if (!is.null(input$mt_output_height)) input$mt_output_height else 42
    }
    u <- if (isTRUE(input$mt_a4_output)) "mm" else {
      if (!is.null(input$mt_output_units)) input$mt_output_units else "cm"
    }
    ggplot2::ggsave(tmp, plot = result, width = w, height = h, units = u,
                    limitsize = FALSE, dpi = 150)
    list(src = tmp, contentType = "image/png", width = "100%",
         alt = "Multi-tree combined plot")
  }, deleteFile = TRUE)

  # --- Log output ---
  output$mt_log <- renderText({
    mt_values$log_messages
  })

  # --- Classification preview ---
  output$mt_classification_preview <- renderUI({
    grps <- mt_values$classification_groups
    if (is.null(grps) || length(grps) == 0) {
      return(tags$p(class = "text-muted", "No classification groups yet. Upload CSV and press Process."))
    }
    tags$div(
      tags$h5(paste0(length(grps), " groups:")),
      lapply(grps, function(g) {
        tags$div(
          style = paste0("padding: 2px 8px; margin: 2px; display: inline-block; ",
                         "background-color: ", g$color, "; border-radius: 4px;"),
          g$display_name
        )
      })
    )
  })

  # --- Configuration YAML output ---
  output$mt_yaml_output <- renderText({
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

  # --- Download YAML config ---
  output$mt_download_yaml <- downloadHandler(
    filename = function() { "multi_tree_config.yaml" },
    content = function(file) {
      shared <- mt_gather_shared()
      yaml_data <- mt_build_yaml_data(
        newick_path = if (!is.null(mt_values$newick_paths)) mt_values$newick_paths[1] else "tree.nwk",
        csv_path = if (!is.null(mt_values$csv_path)) mt_values$csv_path else "data.csv",
        shared = shared,
        per_tree = list()
      )
      writeLines(yaml::as.yaml(yaml_data, indent.mapping.sequence = TRUE), file)
    }
  )

  # --- Download plot ---
  output$mt_download_plot <- downloadHandler(
    filename = function() {
      fmt <- if (!is.null(input$mt_output_format)) input$mt_output_format else "pdf"
      if (isTRUE(input$mt_replace_name) && !is.null(input$mt_custom_name) && nchar(input$mt_custom_name) > 0) {
        paste0(input$mt_custom_name, ".", fmt)
      } else {
        prefix <- if (!is.null(input$mt_prefix_text)) input$mt_prefix_text else "MultiTree__"
        suffix <- if (!is.null(input$mt_suffix_text)) input$mt_suffix_text else ""
        name <- if (!is.null(input$mt_individual_name)) input$mt_individual_name else "Sample1"
        paste0(prefix, name, suffix, ".", fmt)
      }
    },
    content = function(file) {
      req(mt_values$last_plot)
      w <- if (isTRUE(input$mt_a4_output)) 297 else input$mt_output_width
      h <- if (isTRUE(input$mt_a4_output)) 210 else input$mt_output_height
      u <- if (isTRUE(input$mt_a4_output)) "mm" else input$mt_output_units
      ggplot2::ggsave(file, plot = mt_values$last_plot,
                      width = w, height = h, units = u, limitsize = FALSE)
    }
  )

  # --- YAML import observer (option ii: imports all shared settings) ---
  observeEvent(input$mt_yaml_config, {
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
