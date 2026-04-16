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
  tabItem(tabName = "mt_upload", fluidRow(
    box(title = "Upload Data (Multiple Trees)", status = "primary",
        solidHeader = TRUE, width = 4,
        fileInput("mt_tree_files", "Upload Newick Files",
                  multiple = TRUE,
                  accept = c(".nwk", ".newick", ".tree", ".nw", ".txt")),
        tags$hr(),
        fileInput("mt_csv_file", "Upload CSV File",
                  accept = c(".csv", ".tsv", ".txt")),
        tags$hr(),
        selectInput("mt_id_column", "ID Column:", choices = NULL),
        selectInput("mt_individual_column", "Individual Column:", choices = NULL),
        selectizeInput("mt_individual_value", "Individual Value:", choices = NULL),
        checkboxInput("mt_use_all_data", "Use all data (ignore individual filter)", value = FALSE),
        tags$hr(),
        actionButton("mt_process_data", "Process & Match IDs",
                     class = "btn-primary btn-block", icon = icon("check")),
        tags$hr(),
        tags$h5("Optional: Import YAML Configuration"),
        fileInput("mt_yaml_config", "Import YAML", accept = c(".yaml", ".yml")),
        tags$p(class = "text-muted",
               tags$small("Import a YAML to pre-fill shared settings. Heatmap/SNP blocks are ignored."))
    ),
    box(title = "Multi-Tree Preview", status = "primary",
        solidHeader = TRUE, width = 8,
        actionButton("mt_update_preview", "Update Preview",
                     class = "btn-primary", icon = icon("refresh")),
        tags$hr(),
        imageOutput("mt_combined_plot", height = "auto"),
        tags$hr(),
        verbatimTextOutput("mt_log")
    )
  ))
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
    yaml_data <- mt_build_yaml_data(
      newick_path = newick_paths[i],
      csv_path = csv_path,
      shared = shared_settings,
      per_tree = per_tree
    )

    temp_yaml <- tempfile(fileext = ".yaml")
    writeLines(yaml::as.yaml(yaml_data, indent.mapping.sequence = TRUE), temp_yaml)

    # Call existing func.print.lineage.tree — UNCHANGED
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
    csv_path = NULL,            # single temp file path
    csv_data = NULL,            # data.frame
    trees = list(),             # list of treedata objects per newick
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

  # --- Newick upload observer ---
  observeEvent(input$mt_tree_files, {
    req(input$mt_tree_files)
    files <- input$mt_tree_files
    mt_values$newick_paths <- files$datapath
    # Strip extension for display names
    mt_values$newick_names <- tools::file_path_sans_ext(files$name)
    mt_log(paste0("Uploaded ", nrow(files), " newick file(s): ",
                  paste(files$name, collapse = ", ")))

    # Initialize per-tree state
    mt_values$per_tree <- list()
    for (tn in mt_values$newick_names) {
      mt_values$per_tree[[tn]] <- list(
        rotate = NULL,
        highlight_adjust_height = 1,
        highlight_adjust_width = 1.5,
        title = tn,
        bg_color = "#FFFFFF"
      )
    }

    # Update tree selector dropdown
    updateSelectInput(session, "mt_tree_selector",
                      choices = mt_values$newick_names,
                      selected = mt_values$newick_names[1])
  })

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

  # --- Individual column observer ---
  observeEvent(input$mt_individual_column, {
    req(input$mt_individual_column, mt_values$csv_data)
    if (input$mt_individual_column != "") {
      vals <- unique(mt_values$csv_data[[input$mt_individual_column]])
      updateSelectizeInput(session, "mt_individual_value",
                           choices = vals, selected = vals[1])
    }
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
    mt_log("Processing and matching IDs...")

    csv <- mt_values$csv_data
    id_col <- input$mt_id_column

    # Filter by individual if needed
    if (!isTRUE(input$mt_use_all_data) && !is.null(input$mt_individual_column) &&
        input$mt_individual_column != "" && !is.null(input$mt_individual_value)) {
      csv <- csv[csv[[input$mt_individual_column]] == input$mt_individual_value, , drop = FALSE]
      mt_log(paste0("Filtered to individual '", input$mt_individual_value,
                    "': ", nrow(csv), " rows"))
    }

    # Read each newick and match tips
    mt_values$trees <- list()
    for (i in seq_along(mt_values$newick_paths)) {
      tn <- mt_values$newick_names[i]
      tree <- tryCatch(
        treeio::read.newick(mt_values$newick_paths[i]),
        error = function(e) {
          mt_log(paste0("ERROR reading newick '", tn, "': ", e$message))
          NULL
        }
      )
      if (!is.null(tree)) {
        mt_values$trees[[tn]] <- tree
        tip_labels <- tree@phylo$tip.label
        matched <- sum(tip_labels %in% csv[[id_col]])
        mt_log(paste0("Tree '", tn, "': ", length(tip_labels), " tips, ",
                      matched, " matched to CSV"))
      }
    }

    # Auto-generate classification groups from the CSV
    if (!is.null(input$mt_classification_column) && input$mt_classification_column != "") {
      class_col <- input$mt_classification_column
      unique_vals <- sort(unique(na.omit(csv[[class_col]])))
      # Generate colors
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
      mt_log(paste0("Classification: ", length(unique_vals), " groups from column '", class_col, "'"))
    }

    showNotification("Multi-tree processing complete!", type = "message")
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
