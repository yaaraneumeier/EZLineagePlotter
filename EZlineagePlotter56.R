############# part 1

# ============================================================================
# LINEAGE PLOTTER SHINY APP
# A tool for visualizing and analyzing phylogenetic trees
# ============================================================================

# Load required packages
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(shinyBS)
library(colourpicker)
library(DT)
library(ape)
library(dplyr) 
library(ggplot2)
library(ggtree)
library(ggforce)
library(yaml)
library(stringr)
library(scales)
library(tidyverse)
library(ggnewscale)
library(stats)
library(gridExtra)
library(cowplot)  # v144: For proper plot positioning using draw_plot
library(jpeg)     # v144: For JPEG image overlay support
# Note: data.table removed due to function masking issues with dplyr
library(combinat)
library(infotheo)
library(aricode)


options(shiny.reactlog = TRUE)
# v146: Increase max upload size for Shiny server deployments (100MB)
options(shiny.maxRequestSize = 100*1024^2)

# v147: Fixed ellipse x0 positioning in heatmap mode (was using incorrect -15.4 multiplier)
#       Added debug output for PLOT ellipse alpha to compare with legend ellipse alpha
#       Improved legend coordinate output with coordinate system explanations
# v148: Fixed ellipse x0 using SAME data source (p$data) for both tree_max and node positions
#       Fixed legend ellipse size to match plot ellipse (removed man_multiply_elipse and +0.3)
#       Changed y_off_base to position highlight/bootstrap legends on RIGHT side of plot
# v149: Fixed legend ellipse transparency using scales::alpha() to embed alpha in fill color
#       Repositioned legends to TOP-RIGHT using -(y_off_base-5) for x-position
#       Reduced bootstrap title size (0.15 multiplier, capped at 3, min 2)
# v167: OPTION C TRIAL - Native ggplot legends approach
# v168: Fixed legend bleeding with show.legend = c(shape = TRUE) and guides() override
# v169: Added SIZE and ALPHA scale overrides to guides() - still only fixed LAST heatmap
# v170: CRITICAL FIX - Custom key_glyph to prevent ALL legend bleeding
# v171: Implement actual HIGHLIGHT and BOOTSTRAP legends using native ggplot
#       BUG: scale_fill_manual and scale_size_manual REPLACED existing heatmap/P value scales
# v172: Tried ggnewscale::new_scale_fill() and new_scale("size") - still broke P value legend
# v173: SHAPE-ONLY approach - use ONLY shape aesthetic for both legends
#       This avoids ALL scale conflicts: fill (heatmaps), size (P value), colour (classification)
# v174: Fixed Highlight ellipse colors/transparency using force() to capture closure
#       Fixed Bootstrap triangle sizes to vary visibly (both width and height)
# v175: CRITICAL FIX - Use data$.label instead of data$shape in key_glyph
#       data$shape returns numeric shape code (15, 17, 19), NOT the text label
# v176: Fix NULL/empty .label handling - check length() before is.na()
#       as.character(NULL) returns character(0), causing "missing value" error
# v177: Use unique shape values per legend item to identify items in key_glyph
#       Problem: all items mapped to shape=15 made them indistinguishable
# v178: Fix legend bleeding - return nullGrob if shape value not in our set
#       Prevents ellipses/triangles from appearing on other legends
# v179: Enhanced Legend, Extra, and Download tabs with new controls
#       Added tree stretch, background color, symbol & spacing settings
# v180: Multiple enhancements:
#       - Fixed tip labels being hidden by highlight ellipses (layer reordering)
#       - Added legend controls: key width/height, box background, margin
#       - Removed byrow, kept reverse order for legend item control
#       - Fixed download preview to show full page with plot proportions preserved
#       - Fixed Extra tab spinner animation
#       - Bootstrap legend only for triangles format (no legend for percentage/raw)
# S1: First stable release
#       - Performance: Disabled verbose debug output
#       - All v180 features included and tested

# ============================================================================
# VERSION S2.10
# ============================================================================
# S2.10: Fix discrete heatmap colors being lost when Apply Heatmaps is clicked
#        after adding a second heatmap (S2.8/S2.9 didn't fully resolve this)
#        - BUG: Even though S2.8/S2.9 protected stored colors during UI rebuild,
#          the Apply Heatmaps code was reading colors from input[[...]] which
#          could contain stale/default values if the UI hadn't fully updated.
#        - EVIDENCE: Debug showed stored colors were correct in heatmap_configs,
#          but the rendered plot used default palette colors instead.
#        - ROOT CAUSE: Apply Heatmaps (lines 16253-16259) read color values
#          directly from input color pickers instead of from heatmap_configs.
#          When UI rebuilds during "Add Heatmap", input values may be stale.
#        - FIX: Modified Apply Heatmaps to prioritize cfg$custom_colors (stored
#          in heatmap_configs by the color observer) over input values.
#          stored_color is now the source of truth; input is only a fallback.
#        - Also fixed NA color to use stored cfg$na_color over input.
#        - Debug logging ([S2.10-DEBUG]) shows color sources during Apply.
#
# ============================================================================
# VERSION S2.9
# ============================================================================
# S2.9: Fix inhibit_color_save timing issue (S2.8 still had timing problems)
#       - BUG: S2.8's 500ms delay was too short. The UI rebuild process takes longer:
#         heatmap_ui_trigger change -> renderUI runs -> colourInput widgets created
#         -> change events fire. By the time these events fired, 500ms had passed
#         and the inhibit flag was already FALSE.
#       - EVIDENCE: Debug output showed inhibit_color_save=FALSE BEFORE renderUI ran:
#         "[S2.8-DEBUG] inhibit_color_save set to FALSE (after delay)"
#         "[S2.6-DEBUG] renderUI for discrete colors heatmap 1"  <- happens AFTER!
#       - FIX: Increased delay from 500ms to 2500ms to cover full UI rebuild cycle
#       - ADDITIONAL FIX: Skip saving if new color equals stored color (redundant save)
#         This provides defense-in-depth when renderUI creates widgets with stored values
#       - Debug logging ([S2.9-DEBUG]) shows when saves complete after longer delay
#
# S2.8: (S2.9 extended this fix) Added inhibit_color_save flag mechanism
#       - BUG: S2.7's default color detection wasn't reliably preventing resets
#         Adding a second discrete heatmap still reset first heatmap's custom colors
#         (e.g., user's custom pink #E08F8F became default red #E41A1C)
#       - ROOT CAUSE: Timing issue during UI rebuild. Color picker observers were
#         firing with default values before the comparison logic could properly
#         detect and skip them. The S2.7 approach of comparing incoming color to
#         default palette colors had race conditions with Shiny's reactivity.
#       - FIX: Added inhibit_color_save flag mechanism:
#         1. When add/remove/move heatmap is triggered, set inhibit_color_save=TRUE
#         2. Color save observer checks this flag FIRST and skips ALL saves if TRUE
#         3. After delay (UI rebuild complete), flag is reset to FALSE
#         4. This completely blocks any spurious saves during the critical window
#       - This approach is more robust than S2.7's color comparison because it
#         doesn't rely on correctly identifying default vs custom colors
#       - Debug logging ([S2.8-DEBUG]) shows when saves are blocked
#
# S2.7: (Superseded by S2.8) Attempted fix using default color detection
#       - Tried to detect default palette colors and skip saving them
#       - Did not fully resolve the issue due to timing/reactivity problems
#
# S2.6: Debug version for investigating discrete heatmap color reset issue
#       - Added detailed debug logging for custom_colors storage and retrieval
#       - Debug points: add_new_heatmap observer, discrete color picker observers,
#         renderUI for discrete colors
#       - This version helps diagnose why colors reset when adding second heatmap
#
# S2.5: Fix heatmap column order changing when adding second heatmap
#       - BUG: Adding a second heatmap caused first heatmap's column order to change
#         (e.g., NRAS,BRAF,MET became BRAF,NRAS,MET) which changed heatmap appearance
#       - ROOT CAUSE: Shiny selectize reorders selected items to match the order
#         they appear in 'choices' when the UI rebuilds. Since choices were alphabetical,
#         the user's selection order was lost after clicking "Add Heatmap".
#       - FIX: Reorder the 'choices' list to put already-selected columns first
#         (in their original order) before remaining columns. This ensures selectize
#         preserves the user's intended column order through UI rebuilds.
#
# S2.4: Fix discrete heatmap colors resetting when adding second heatmap
#       - BUG: Adding a second discrete heatmap caused first heatmap colors to reset
#       - ROOT CAUSE: Colors were only stored in Shiny inputs (not heatmap_configs).
#         When "Add Heatmap" is clicked, the UI rebuilds (heatmap_ui_trigger).
#         During rebuild, old color picker inputs don't exist yet, so
#         isolate(input[[...]]) returns NULL and default palette colors were used.
#       - FIX: Added observer to persist color picker values to heatmap_configs.
#         Modified renderUI to use stored colors as fallback when inputs don't exist.
#
# S2.3: Added debug output for multiple discrete heatmap column mapping investigation
#       - Debug output shows column selections in apply handler
#       - Debug output shows YAML according columns for each heatmap
#       - Debug output shows column extraction during rendering
#
# S2.2: Bug fixes for discrete heatmap colors
#       - Fixed discrete heatmap colors not changing when color pickers changed
#       - Fixed heatmap legend colors not matching heatmap tile colors
#       - Root cause: Color picker values were collected using unfiltered data
#         but color pickers were created using filtered data (tree tips only),
#         causing wrong value-to-color mappings
# S2.1: Performance & stability improvements
#       - Two-tier caching system for faster plot updates
#       - Async garbage collection for responsive UI
#       - SVG preview rendering (faster than PNG)
#       - Server-side selectize for large datasets
#       - Fix cascading plot regeneration that freezes UI
#       - Tab switch optimization with ignoreInit
#       - RData heatmap fix for custom classification path
#       - Improved CNV sample matching with better diagnostics
#       - Dropdown for RData sample name mapping column
#       - CSV loading fixes for classification coloring
#       - Manual rotation and multifurcating node rotation fixes
# S2.0: Major stable release with RData CNV heatmap support
#       - RData CNV file import for heatmaps (from QDNAseq/scIMPACT pipelines)
#       - Automatic sample matching via CSV lookup columns
#       - Red-white-blue default color scheme for CNV data
#       - Vertical column lines option for heatmaps
#       - Horizontal row lines option for heatmaps
#       - Fixed multiple heatmap support (CSV + RData together)
# S1.6: Stable release with performance optimizations and bug fixes
#       - All S1.x improvements consolidated and tested
# S1.5: Fixed Legend Background not working
#       - Added legend.background theme element for individual legend panel backgrounds
#       - Fixed RGBA color handling from colourpicker (convert to RGB)
# S1.4: Performance optimization - Reduce redundant recalculations
#       - Added plot_trigger mechanism to batch rapid updates (100ms debounce)
#       - Converted observers to use request_plot_update() instead of direct generate_plot()
#       - Removed duplicate guards in generate_plot() function
#       - Reduced debug logging overhead
# S1.3: Performance optimization - Layer management
#       - Reduced func.move.tiplabels.to.front() calls from 3+ per render to 1
#       - Removed redundant intermediate calls in func_highlight and func.print.lineage.tree
#       - Optimized function: uses vapply instead of sapply, removed debug I/O overhead
#       - Layer reordering now happens ONCE at the end in generate_plot()
# S1.2: Fixed undefined x_range_min in func_highlight causing "Problem while
#       computing aesthetics" error when adding 2+ highlights with a heatmap.
VERSION <- "S2.292dev"

# Debug output control - set to TRUE to enable verbose console logging
# For production/stable use, keep this FALSE for better performance
DEBUG_VERBOSE <- FALSE

# Helper function for debug output - only prints when DEBUG_VERBOSE is TRUE
debug_cat <- function(...) {
  if (DEBUG_VERBOSE) {
    cat(file = stderr(), ...)
  }
}

###### part 1 a:
# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#read yaml file
# Function to read and parse YAML files
func.read_yaml <- function(yaml_path) {
  tryCatch({
    # Check if file exists
    if (!file.exists(yaml_path)) {
      stop(paste("YAML configuration file not found:", yaml_path))
    }
    
    # Read YAML file
    yaml_content <- yaml::read_yaml(yaml_path)
    
    # Validate basic structure
    if (!is.list(yaml_content) || 
        !"Individual general definitions" %in% names(yaml_content)) {
      stop("Invalid YAML configuration format")
    }
    
    return(yaml_content)
  }, error = function(e) {
    message(paste("Error reading YAML file:", e$message))
    return(NULL)
  })
}


# Function to get filename without extension from a filepath
get_filename <- function(filepath) {
  # Extract the base name (file name with extension) from the full path
  base_name <- basename(filepath)
  
  # Remove the file extension
  file_name <- tools::file_path_sans_ext(base_name)
  
  return(file_name)
}



# Function to check for string/integer mismatches and common prefixes/suffixes
match_tree_ids_with_csv <- function(tree_labels, csv_ids) {
  # S1-PERF: Optimized with vectorized operations instead of nested loops

  # Ensure tree_labels and csv_ids are character vectors
  tree_labels <- as.character(tree_labels)
  csv_ids <- as.character(csv_ids)

  # Initialize matching results
  matched <- list(
    exact_matches = character(0),
    prefix_suffix_matches = list(),
    numeric_matches = list(),
    unmatched = character(0)
  )

  # Create a copy of original tree labels
  tree_labels_original <- tree_labels

  # Try exact matches first (already vectorized - good)
  exact_matches <- intersect(tree_labels, csv_ids)
  matched$exact_matches <- exact_matches

  # Remove exact matches from consideration
  remaining_tree_labels <- setdiff(tree_labels, exact_matches)

  if (length(remaining_tree_labels) > 0) {
    # S1-PERF: Vectorized numeric conversion
    numeric_tree <- suppressWarnings(as.numeric(remaining_tree_labels))
    numeric_csv <- suppressWarnings(as.numeric(csv_ids))

    valid_numeric_tree <- !is.na(numeric_tree)
    valid_numeric_csv <- !is.na(numeric_csv)

    # S1-PERF: Use match() for numeric matching instead of loop
    if (any(valid_numeric_tree) && any(valid_numeric_csv)) {
      numeric_tree_valid <- numeric_tree[valid_numeric_tree]
      numeric_csv_valid <- numeric_csv[valid_numeric_csv]
      csv_ids_valid <- csv_ids[valid_numeric_csv]

      # Use match to find all numeric matches at once
      match_indices <- match(numeric_tree_valid, numeric_csv_valid)
      matched_mask <- !is.na(match_indices)

      if (any(matched_mask)) {
        matched_tree_labels <- remaining_tree_labels[valid_numeric_tree][matched_mask]
        matched_csv_ids <- csv_ids_valid[match_indices[matched_mask]]
        for (i in seq_along(matched_tree_labels)) {
          matched$numeric_matches[[matched_tree_labels[i]]] <- matched_csv_ids[i]
        }
      }
    }

    # S1-PERF: Vectorized prefix/suffix matching
    # Only process labels not already matched
    labels_to_check <- setdiff(remaining_tree_labels, names(matched$numeric_matches))

    if (length(labels_to_check) > 0 && length(csv_ids) > 0) {
      # Pre-filter: remove NA values
      valid_csv <- csv_ids[!is.na(csv_ids)]

      for (tree_label in labels_to_check) {
        if (is.na(tree_label)) {
          matched$unmatched <- c(matched$unmatched, tree_label)
          next
        }

        # S1-PERF: Use startsWith/endsWith which are much faster than grepl
        # Check if tree_label is prefix of any csv_id
        prefix_matches <- valid_csv[startsWith(valid_csv, tree_label)]

        # Check if tree_label is suffix of any csv_id
        suffix_matches <- valid_csv[endsWith(valid_csv, tree_label)]

        # Check if any csv_id is prefix of tree_label
        csv_as_prefix <- valid_csv[vapply(valid_csv, function(x) startsWith(tree_label, x), logical(1))]

        # Check if any csv_id is suffix of tree_label
        csv_as_suffix <- valid_csv[vapply(valid_csv, function(x) endsWith(tree_label, x), logical(1))]

        all_matches <- unique(c(prefix_matches, suffix_matches, csv_as_prefix, csv_as_suffix))

        if (length(all_matches) > 0) {
          matched$prefix_suffix_matches[[tree_label]] <- all_matches
        } else {
          matched$unmatched <- c(matched$unmatched, tree_label)
        }
      }
    }
  }

  # S1-PERF: Vectorized mapping construction
  final_mapping <- list()

  # Add exact matches (vectorized)
  if (length(matched$exact_matches) > 0) {
    for (label in matched$exact_matches) {
      final_mapping[[label]] <- label
    }
  }

  # Add numeric matches (pick first if multiple)
  for (tree_label in names(matched$numeric_matches)) {
    final_mapping[[tree_label]] <- matched$numeric_matches[[tree_label]][1]
  }

  # Add prefix/suffix matches (pick first if multiple)
  for (tree_label in names(matched$prefix_suffix_matches)) {
    final_mapping[[tree_label]] <- matched$prefix_suffix_matches[[tree_label]][1]
  }

  # Infer trimming parameters from the matching results
  id_tip_trim_flag <- FALSE
  id_tip_trim_start <- 1
  id_tip_trim_end <- NA
  id_tip_prefix <- ""

  # Check if we need trimming based on the matches
  if (length(matched$exact_matches) == length(tree_labels_original)) {
    # Perfect exact matches - no trimming needed
    id_tip_trim_flag <- FALSE
  } else if (length(matched$prefix_suffix_matches) > 0) {
    # Analyze prefix/suffix patterns to determine trimming

    # Sample a few matches to determine pattern
    sample_matches <- head(matched$prefix_suffix_matches, 5)

    for (tree_label in names(sample_matches)) {
      csv_match <- matched$prefix_suffix_matches[[tree_label]][1]

      # Check if tree label has a prefix that csv doesn't have
      if (nchar(tree_label) > nchar(csv_match)) {
        # Tree label is longer - likely has a prefix
        # Find where they start matching
        for (start_pos in 1:nchar(tree_label)) {
          substring_tree <- substring(tree_label, start_pos)
          if (grepl(csv_match, substring_tree, fixed = TRUE)) {
            id_tip_trim_flag <- TRUE
            id_tip_trim_start <- start_pos
            id_tip_trim_end <- nchar(tree_label)
            break
          }
        }
        break
      } else if (nchar(csv_match) > nchar(tree_label)) {
        # CSV ID is longer - tree might be missing prefix
        # Find if tree label is contained in csv
        if (grepl(tree_label, csv_match, fixed = TRUE)) {
          # Determine what prefix to add
          prefix_start <- regexpr(tree_label, csv_match, fixed = TRUE)[1]
          if (prefix_start > 1) {
            id_tip_prefix <- substring(csv_match, 1, prefix_start - 1)
          }
        }
        break
      }
    }
  } else if (length(matched$numeric_matches) > 0) {
    # Numeric matches - might need trimming
    id_tip_trim_flag <- TRUE
  }

  # Add trimming parameters to the return
  return(list(
    mapping = final_mapping,
    matched = matched,
    summary = list(
      total_tree_labels = length(tree_labels_original),
      exact_matches = length(matched$exact_matches),
      numeric_matches = length(matched$numeric_matches),
      prefix_suffix_matches = length(matched$prefix_suffix_matches),
      unmatched = length(matched$unmatched)
    ),
    trimming_params = list(
      id_tip_trim_flag = id_tip_trim_flag,
      id_tip_trim_start = id_tip_trim_start,
      id_tip_trim_end = id_tip_trim_end,
      id_tip_prefix = id_tip_prefix
    )
  ))


}

# Function to check if a value is TRUE/YES in various formats
func.check.bin.val.from.conf <- function(val) {
  # v53: print("val is")
  # v53: print(val)
  if (tolower(val) == "yes" || val == TRUE || tolower(val) == "true") {
    out <- TRUE
  } else {
    out <- FALSE
  }
  return(out)
}

# v82: Enhanced function to repair corrupted ggtree/ggplot mapping attribute
# This fixes the error: "@mapping must be <ggplot2::mapping>, not S3<data.frame>"
# which occurs in newer versions of ggplot2 (3.4+) when gheatmap or other operations
# accidentally corrupt the mapping slot
# CRITICAL v82 FIX: Do NOT reset layer mappings - only fix top-level mapping
# Resetting layer mappings (like fill = value for GeomTile) breaks the heatmap
func.repair.ggtree.mapping <- function(p, verbose = FALSE) {
  repaired <- FALSE

  # Check if top-level mapping is valid (should be class "uneval" from aes())
  if (!inherits(p$mapping, "uneval")) {
    if (verbose) {
      debug_cat(paste0("\n=== v82: Repairing corrupted plot mapping ===\n"))
      debug_cat(paste0("  Original mapping class: ", paste(class(p$mapping), collapse=", "), "\n"))
    }

    tryCatch({
      p$mapping <- aes()
      repaired <- TRUE
      if (verbose) {
        debug_cat(paste0("  Fixed mapping class: ", paste(class(p$mapping), collapse=", "), "\n"))
      }
    }, error = function(e) {
      debug_cat(paste0("  v2: Could not repair mapping: ", e$message, "\n"))
    })
  }

  # v82: REMOVED layer mapping repairs - these were breaking the heatmap fill mapping
  # The layer mappings are set correctly by gheatmap and should not be modified.
  # Only the top-level plot mapping sometimes gets corrupted to a data.frame.

  if (repaired && verbose) {
    debug_cat(paste0("================================\n"))
  }

  return(p)
}

# v180: Function to move tip label layers to the end of the layer stack
# This ensures tip labels render ON TOP of other elements (highlight ellipses, heatmaps, etc.)
# Called ONCE after all layers are added and before final rendering (in generate_plot)
# S1.3-PERF: Optimized - removed redundant debug logging, uses vectorized operations
func.move.tiplabels.to.front <- function(p, verbose = FALSE) {
  # S1.3-PERF: Removed per-call debug logging to reduce I/O overhead
  # Only log when verbose=TRUE (disabled by default now)

  # Early exit: no layers to process
  if (is.null(p$layers) || length(p$layers) == 0) {
    return(p)
  }

  # S1.3-PERF: Use vapply for type-safe vectorized operation (faster than sapply)
  layer_types <- vapply(p$layers, function(l) class(l$geom)[1], character(1))

  # Find text-type layers (tip labels: GeomText, GeomLabel, or GeomTiplab)
  tiplab_indices <- which(layer_types %in% c("GeomText", "GeomLabel", "GeomTiplab"))

  # Early exit: no text layers to move
  if (length(tiplab_indices) == 0) {
    return(p)
  }

  # S1.3-PERF: Single-operation layer reordering (vectorized)
  other_indices <- setdiff(seq_along(p$layers), tiplab_indices)
  p$layers <- c(p$layers[other_indices], p$layers[tiplab_indices])

  if (verbose) {
    debug_cat(paste0("[PERF] func.move.tiplabels.to.front: moved ", length(tiplab_indices),
                     " text layer(s) to front\n"))
  }

  return(p)
}

# v71: Function to diagnose which layer is causing ggplot_build to fail
func.diagnose.layer.issues <- function(p, verbose = TRUE) {
  if (verbose) {
    debug_cat(paste0("\n=== v71: DIAGNOSING LAYER ISSUES ===\n"))
  }

  problematic_layers <- c()

  if (!is.null(p$layers) && length(p$layers) > 0) {
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      geom_class <- class(layer$geom)[1]

      # Try to set up this layer's data
      layer_ok <- tryCatch({
        # Get plot data
        plot_data <- if (!is.null(p$data)) p$data else data.frame()

        # Try basic setup
        if (!is.null(layer$data)) {
          if (is.function(layer$data)) {
            test_data <- layer$data(plot_data)
          } else {
            test_data <- layer$data
          }
        }
        TRUE
      }, error = function(e) {
        if (verbose) {
          debug_cat(paste0("  Layer ", i, " (", geom_class, "): FAILED - ", e$message, "\n"))
        }
        FALSE
      })

      if (!layer_ok) {
        problematic_layers <- c(problematic_layers, i)
      } else if (verbose) {
        debug_cat(paste0("  Layer ", i, " (", geom_class, "): OK\n"))
      }
    }
  }

  if (verbose) {
    if (length(problematic_layers) == 0) {
      debug_cat(paste0("  No obvious layer issues detected\n"))
    } else {
      debug_cat(paste0("  Problematic layers: ", paste(problematic_layers, collapse=", "), "\n"))
    }
    debug_cat(paste0("================================\n"))
  }

  return(problematic_layers)
}


# Parse YAML configuration file
parse_yaml_config <- function(file_path) {
  tryCatch({
    config <- yaml::read_yaml(file_path)
    return(config)
  }, error = function(e) {
    return(list(error = paste("Failed to parse YAML file:", e$message)))
  })
}

# S1.62dev: Extract CNV data from RData file
# The RData file should contain 'results_CNV_tool_final' with nested structure:
# results_CNV_tool_final[[sample_name]][[1]]$copy contains the CNV values
# results_CNV_tool_final[[sample_name]][[1]]$chr may contain chromosome info
# results_CNV_tool_final[[sample_name]][[1]]$start/end may contain positions
# Returns a list with:
#   - matrix: rows = genomic positions, columns = samples
#   - sample_names: vector of sample names
#   - chr_info: chromosome information for each position (if available)
#   - position_info: start/end positions (if available)
#   - error: error message if failed (NULL if success)
func.extract.cnv.from.rdata <- function(rdata_path, downsample_factor = 10) {
  tryCatch({
    # Load RData into a new environment to avoid polluting global namespace
    env <- new.env()
    load(rdata_path, envir = env)

    # Check if results_CNV_tool_final exists
    if (!"results_CNV_tool_final" %in% names(env)) {
      return(list(error = "RData file does not contain 'results_CNV_tool_final'"))
    }

    results <- env$results_CNV_tool_final
    sample_names <- names(results)

    if (length(sample_names) == 0) {
      return(list(error = "No samples found in results_CNV_tool_final"))
    }

    cat(file=stderr(), paste0("[RDATA-CNV] Found ", length(sample_names), " samples\n"))

    # Extract CNV data from each sample
    cnv_matrix <- NULL
    chr_info <- NULL
    start_info <- NULL
    end_info <- NULL
    first <- TRUE

    for (sample_name in sample_names) {
      # Navigate to the copy data: results_CNV_tool_final[[sample]][[1]]$copy
      sample_data <- results[[sample_name]]

      if (is.null(sample_data) || length(sample_data) == 0) {
        cat(file=stderr(), paste0("[RDATA-CNV] Warning: No data for sample ", sample_name, "\n"))
        next
      }

      # Get the copy values (first element's $copy)
      copy_data <- NULL
      if (is.list(sample_data) && length(sample_data) >= 1) {
        if (!is.null(sample_data[[1]]$copy)) {
          copy_data <- as.data.frame(sample_data[[1]]$copy)
        }

        # S1.62dev: Extract chromosome and position information (only on first valid sample)
        if (first && !is.null(sample_data[[1]])) {
          # Check available fields in the data structure
          available_fields <- names(sample_data[[1]])
          cat(file=stderr(), paste0("[RDATA-CNV] Available fields in sample data: ", paste(available_fields, collapse=", "), "\n"))

          # Extract chromosome info if available
          if (!is.null(sample_data[[1]]$chr)) {
            chr_info <- sample_data[[1]]$chr
            cat(file=stderr(), paste0("[RDATA-CNV] Found chromosome info: ", length(chr_info), " entries\n"))
            cat(file=stderr(), paste0("[RDATA-CNV] Unique chromosomes: ", paste(unique(chr_info), collapse=", "), "\n"))
          } else if ("chr" %in% available_fields) {
            cat(file=stderr(), "[RDATA-CNV] chr field exists but is NULL\n")
          }

          # Extract start positions if available
          if (!is.null(sample_data[[1]]$start)) {
            start_info <- sample_data[[1]]$start
            cat(file=stderr(), paste0("[RDATA-CNV] Found start positions: ", length(start_info), " entries\n"))
          }

          # Extract end positions if available
          if (!is.null(sample_data[[1]]$end)) {
            end_info <- sample_data[[1]]$end
            cat(file=stderr(), paste0("[RDATA-CNV] Found end positions: ", length(end_info), " entries\n"))
          }
        }
      }

      if (is.null(copy_data)) {
        cat(file=stderr(), paste0("[RDATA-CNV] Warning: No copy data for sample ", sample_name, "\n"))
        next
      }

      colnames(copy_data) <- sample_name

      if (first) {
        cnv_matrix <- copy_data
        first <- FALSE
      } else {
        # Ensure same number of rows before binding
        if (nrow(copy_data) == nrow(cnv_matrix)) {
          cnv_matrix <- cbind(cnv_matrix, copy_data)
        } else {
          cat(file=stderr(), paste0("[RDATA-CNV] Warning: Row mismatch for sample ", sample_name,
                                    " (", nrow(copy_data), " vs ", nrow(cnv_matrix), ")\n"))
        }
      }
    }

    if (is.null(cnv_matrix) || ncol(cnv_matrix) == 0) {
      return(list(error = "No valid CNV data could be extracted"))
    }

    cat(file=stderr(), paste0("[RDATA-CNV] Extracted matrix: ", nrow(cnv_matrix), " rows x ", ncol(cnv_matrix), " columns\n"))

    # Replace dots with dashes in column names (sample names)
    colnames(cnv_matrix) <- gsub("\\.", "-", colnames(cnv_matrix))

    # Downsample rows if requested - also downsample chr/position info
    if (downsample_factor > 1) {
      row_indices <- unique((1:nrow(cnv_matrix) %/% downsample_factor) * downsample_factor)
      row_indices <- row_indices[row_indices > 0 & row_indices <= nrow(cnv_matrix)]
      cnv_matrix <- cnv_matrix[row_indices, , drop = FALSE]

      # Also downsample chr and position info if available
      if (!is.null(chr_info) && length(chr_info) >= max(row_indices)) {
        chr_info <- chr_info[row_indices]
      }
      if (!is.null(start_info) && length(start_info) >= max(row_indices)) {
        start_info <- start_info[row_indices]
      }
      if (!is.null(end_info) && length(end_info) >= max(row_indices)) {
        end_info <- end_info[row_indices]
      }

      cat(file=stderr(), paste0("[RDATA-CNV] After downsampling (factor ", downsample_factor, "): ",
                                nrow(cnv_matrix), " rows\n"))
    }

    # S1.62dev: Transpose matrix so rows = samples and columns = genomic positions
    # gheatmap() expects row names to match tree tip labels
    cnv_matrix <- t(cnv_matrix)
    cat(file=stderr(), paste0("[RDATA-CNV] After transpose: ", nrow(cnv_matrix), " samples x ", ncol(cnv_matrix), " positions\n"))

    # S1.62dev: Build chromosome summary if available
    chr_summary <- NULL
    if (!is.null(chr_info)) {
      chr_summary <- paste0(
        "Chromosomes found: ", paste(unique(chr_info), collapse=", "), "\n",
        "Positions per chromosome:\n"
      )
      chr_counts <- table(chr_info)
      for (chr_name in names(chr_counts)) {
        chr_summary <- paste0(chr_summary, "  ", chr_name, ": ", chr_counts[[chr_name]], " positions\n")
      }
      cat(file=stderr(), paste0("[RDATA-CNV] ", chr_summary))
    }

    return(list(
      matrix = as.matrix(cnv_matrix),
      sample_names = rownames(cnv_matrix),  # Now samples are rows
      n_positions = ncol(cnv_matrix),       # Positions are now columns
      n_samples = nrow(cnv_matrix),
      chr_info = chr_info,                  # S1.62dev: Chromosome info (position-wise)
      start_info = start_info,              # S1.62dev: Start positions
      end_info = end_info,                  # S1.62dev: End positions
      chr_summary = chr_summary,            # S1.62dev: Summary string
      error = NULL
    ))

  }, error = function(e) {
    return(list(error = paste("Failed to extract CNV data:", e$message)))
  })
}


# Convert to YAML
#yaml_text <- as.yaml(yaml_structure, indent.mapping.sequence = TRUE)
#return(yaml_text)
#}

# Convert settings to YAML format
#settings_to_yaml <- function(settings) {
#  # Create YAML from settings
#  yaml_text <- yaml::as.yaml(settings, indent.mapping.sequence = TRUE)
#  return(yaml_text)
#}


# Check if two tree tips should be connected
check_lines_intersect_2d_coords <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
  # Calculate direction vectors for the two lines
  d1 <- c(x2 - x1, y2 - y1)  # Direction vector of the first line
  d2 <- c(x4 - x3, y4 - y3)  # Direction vector of the second line
  
  # Check if the lines are parallel using the determinant (cross product)
  det <- d1[1] * d2[2] - d1[2] * d2[1]
  
  if (det == 0) {
    return(FALSE)  # Lines are parallel and do not intersect
  }
  
  # Solve for the intersection point using parametric equations
  A <- matrix(c(d1[1], -d2[1], d1[2], -d2[2]), nrow = 2)
  b <- c(x3 - x1, y3 - y1)
  
  sol <- solve(A, b)
  t <- sol[1]
  s <- sol[2]
  
  # Check if the intersection point lies within the segments
  if (t >= 0 && t <= 1 && s >= 0 && s <= 1) {
    return(TRUE)  # The lines intersect within the segments
  } else {
    return(FALSE)  # The lines do not intersect within the segments
  }
}

##### part 2

# ============================================================================
# CORE TREE PROCESSING FUNCTIONS
# ============================================================================

func.make.list_id_by_class <- function(cls_num, cls_renaming_list,yaml_file,
                                       title.id,leaves_id_from_tree,readfile440,acc,debug_mode,
                                       id_tip_trim_flag= TRUE,id_tip_trim_start=3,id_tip_trim_end=7) {
  
  #print("In func.make.list_id_by_class")
  #print(readfile440)
  #print("_____")
  #print("leaves_id_from_tree is")
  #print(leaves_id_from_tree)
  
  
  list_id_by_class <- c()
  for (cls in cls_renaming_list) {
    
    list_id_by_class[[cls]] <- ""
  }
  
  
  #print("leaves_id_from_tree is")
  #print(leaves_id_from_tree)
  
  
  list_id_is_in_clss_flags <- c()
  for (id_tip in leaves_id_from_tree) {
    list_id_is_in_clss_flags[[id_tip]] <- 0
  }
  
  #print("title.id is")
  #print(title.id)
  #print("readfile440[[title.id] is")
  #print(readfile440[title.id,])
  #print("___________")
  #print('readfile440 is')
  #print(readfile440)
  #print('________')
  for (id_tip in leaves_id_from_tree) {
    
    #question
    
    if (id_tip_trim_flag == TRUE) {
      # Trimming ENABLED: apply substr to extract portion of tip ID
      #print("In TRUE tip - applying trimming")
      #print(substr(id_tip,id_tip_trim_start,id_tip_trim_end))
      row <- readfile440[readfile440[[title.id]] == as.integer(substr(id_tip,id_tip_trim_start,id_tip_trim_end)),]
      
    } else {
      # Trimming DISABLED: use full tip ID as-is
      #print("In FALSE tip - using full ID")
      row <- readfile440[readfile440[[title.id]] == id_tip,]
      
    }
    #print("row is")
    #print(row)
    
    flag_is_in <-TRUE
    if (nrow(row)==0) {
      flag_is_in <- FALSE
    }
    
    
    
    
    
    for (ac in acc) {
      
      flag_op <-TRUE
      
      is_in_class <- c()
      index_op <- 1
      while(flag_op== TRUE) {
        
        
        title_i <- paste0('title',index_op)
        value_i <- paste0('value',index_op)
        
        index_op <- index_op+1
        
        names_ac <- names(ac)
        
        
        
        title_i_name <- ac[[names_ac]][[title_i]]
        value_i_name <- ac[[names_ac]][[value_i]]
        
        
        
        if ( (paste0('title',index_op)) %in% names(ac[[names_ac]])) {
          
        } else {
          
          flag_op <-FALSE
        }
        
        if (flag_is_in== TRUE) {
          val <- row[[title_i_name]]
          is_in <- func.check.if.id.in.sub.class(val, title_i_name, value_i_name)
        } else {
          val <- NA
          is_in <- FALSE
        }
        
        is_in_class <- c(is_in_class,is_in)
        
        
      }
      
      
      
      
      if (sum(is_in_class)== length(is_in_class)) {
        is_in_class_res<- TRUE
      } else {
        is_in_class_res<- FALSE
      }
      
      
      
      if (is_in_class_res == TRUE) {
        
        
        if (list_id_is_in_clss_flags[[id_tip]]==0) {
          
          cls <- ac[[names_ac]][['display_name']]
          
          
          list_id_is_in_clss_flags[[id_tip]] <-1
          
          
          list_id_by_class[[cls]] <- c(list_id_by_class[[cls]], id_tip)
          
        } else {
          # v53: print(paste0("Warining: double mapping for ", id_tip))
        }
      }
    }
  }
  #print("length(list_id_is_in_clss_flags) is")
  #print(length(list_id_is_in_clss_flags))
  #print("length(list_id_is_in_clss_flags[list_id_is_in_clss_flags==0]) is")
  #print(length(list_id_is_in_clss_flags[list_id_is_in_clss_flags==0]))
  #if (length(list_id_is_in_clss_flags[list_id_is_in_clss_flags==0])>0) {
  #  print("Warning: the following cells are not classified")
  #  print(list_id_is_in_clss_flags[list_id_is_in_clss_flags==0])
  #}
  
  #print("returned list_id_by_class is")
  # print(list_id_by_class)
  
  
  list_id_by_class<- list_id_by_class
}

# Function to create lists of tips
func.create.cc_tipss <- function(d440){
  subframe_of_tips <- d440[d440$isTip=="TRUE",2]
  list_tips <- list(subframe_of_tips$node)
  cc_tips <- c(list_tips)
  cc_tipss <- cc_tips[[1]]
  return(cc_tipss)
}

# Function to create lists of intermediate nodes
func.create.cc_nodss <- function(d440){
  subframe_of_nodes <- d440[d440$isTip=="FALSE",2]
  list_nodes <- list(subframe_of_nodes$node)
  cc_nodes <- c(list_nodes)
  cc_nodss <- cc_nodes[[1]]
  return(cc_nodss)
}

# Function to create lists of all nodes
func.create.cc_totss <- function(d440){
  list_all_tree <- list(d440$node)
  cc_tot <- c(list_all_tree)
  cc_totss <- cc_tot[[1]]
  return(cc_totss)
}


# Wrapper to get subtree info (handles different ggtree versions)
func.get.subtree.wrapper.return.list.full <- function(tree, nod) {
  su <- ggtree:::getSubtree(tree, nod)  # kids contain all the nods in the subtree including nod (its root)
  check_df <- inherits(su, "data.frame")
  if (check_df == TRUE) {
    out <- su$node
  } else {
    out <- su
  }
  out <- as.numeric(out)
  return(out)
}


# Wrapper to get subtree tips (handles different ggtree versions)
func.get.subtree.wrapper.return.list.tips <- function(tree, nod) {
  su <- ggtree:::getSubtree(tree, nod)  # kids contain all the nods in the subtree including nod (its root)
  check_df <- inherits(su, "data.frame")
  if (check_df == TRUE) {
    sub <- subset(su, isTip == TRUE)
    out <- sub$node
  } else {
    # v56b: Suppress harmless fortify warnings
    subb <- subset(suppressWarnings(ggtree(tree))$data, node %in% c(su))
    out <- subset(subb, isTip == TRUE)
    out <- out$node
  }
  return(out)
}


# Fix missing leaves in readfile
fix.readfile440.with.missing.leaves <- function(readfile440, title.id, tree440, ids_list, no_name,
                                                id_tip_trim_flag, id_tip_trim_start, id_tip_trim_end, id_tip_prefix, debug_mode = FALSE) {
  readfile440n <- readfile440
  #print("in fix.readfile440.with.missing.leaves")

  # v56b: Suppress harmless fortify warnings
  tree_data <- suppressWarnings(ggtree(tree440))$data
  leaves_id_from_tree1 <- tree_data[tree_data$isTip == TRUE, 'label']
  leaves_id_from_tree <- as.list(leaves_id_from_tree1['label'])$label
  #print(leaves_id_from_tree)
  #print("id_tip_trim_flag is")
  #print(id_tip_trim_flag)
  
  
  if (debug_mode == TRUE) {
    # v53: print("leaves_id_from_tree is")
    # v53: print(leaves_id_from_tree)
  }
  
  for (id in leaves_id_from_tree) {
    if (id == "root") {
      idd <- "root"
    } else {
      if (id_tip_trim_flag == TRUE) {
        idd <- substring(id, id_tip_trim_start, nchar(id))
      } else {
        idd <- id
      }
    }
    
    if (!(idd %in% ids_list)) {
      place <- nrow(readfile440n) + 1
      readfile440n[place,] = rep(NA, length(colnames(readfile440n)))
      readfile440n[place, title.id] <- idd
    }
  } 
  
  if (id_tip_trim_flag == TRUE) {
    if (debug_mode == TRUE) {
      # v53: print("cell names originally from")
      # v53: print(ids_list)
      # v53: print("now have the prefix from id_tip_prefix")
      # v53: print(id_tip_prefix)
    }
  }
  
  ids_list <- as.character(ids_list)
  ineter <- intersect(leaves_id_from_tree, ids_list)
  
  if (length(ineter) == 0) {
    # v53: print("ERROR id tips from tree and csv don't match")
  } else {
    not_in_csv <- setdiff(leaves_id_from_tree, ids_list)
    if (length(not_in_csv) > 0) {
      # v53: print(paste0(paste0("WARNING: cant find ", not_in_csv), " in csv"))
    }
    
    readfile440n <- readfile440n[readfile440n[[title.id]] %in% ineter, ]
  }
  
  return(readfile440n)
}



# Check if an ID is in a subclass
func.check.if.id.in.sub.class <- function(val, title_i_name, value_i_name) {
  
  is_in <- FALSE
  if (is.null(val)) {
    return(is_in)
  }
  
  if (length(value_i_name) > 1) {
    if (length(val) == 1 && val %in% value_i_name) {
      is_in <- TRUE
    }
    if (length(val) == 1 && is.na(val)) {
      if ("na" %in% value_i_name || "NA" %in% value_i_name) {
        is_in <- TRUE
      }
    }
  } else {
    if (substr(value_i_name, 1, 1) == '(') {
      be <- strsplit(value_i_name, "-")
      begin <- be[[1]]
      end <- begin[2]
      begin <- begin[1]
      begin <- as.numeric(substr(begin, 2, nchar(begin)))
      begin <- begin[1]
      end <- as.numeric(substr(end, 1, nchar(end) - 1))
      
      if (is.na(begin) | is.na(end)) {
        # v53: print("Error")
        # v53: print(paste0(value_i_name, " is an illegal value inside brackets"))
        # v53: print("Inside brackets numerical range should be placed")
        # v53: print("If the brackets are a part of the desired value, put them inside a list using []")
        stop(" ")
      }
      
      if (begin >= end) {
        stop("illegal range")
      }
      
      if (length(val) == 1 && !is.na(val) && (begin <= val) && (val <= end)) {
        is_in <- TRUE
      }
    } else {
      if (length(val) == 1 && is.na(val)) {
        if (tolower(value_i_name) == "na") {
          is_in <- TRUE
        }
      } else {
        if (length(val) == 1 && val == value_i_name) {
          is_in <- TRUE
        }
      }
    }
  }
  
  return(is_in)
}


# Create list of classes
function.create.cls.list <- function(ids_list, dx_rx_types1_short, list_id_by_class, na_name) {
  cls <- c()
  for (j in ids_list) {
    l <- ""
    for (n in dx_rx_types1_short) {
      b <- as.character(j)
      if (b %in% list_id_by_class[[n]]) {
        l <- n
      }
    } 
    
    if (nchar(l) == 0) {
      l <- na_name
    }
    
    cls <- c(cls, l)
  }
  
  return(cls)
}



# Create a list of nodes by class
func.create.list_node_by_class <- function(dx_rx_types1_short, list_id_by_class,
                                           dxdf440_dataf, title.id, tree_size, d440, cc_totss, debug_mode = FALSE,
                                           id_tip_trim_flag, id_tip_prefix) { 
  
  classication <- "Mapping"
  if (debug_mode == TRUE) {
    # v53: print("In func.create.list_node_by_class")
  }
  
  # Mapping nodes to classes
  convert_id_to_node <- func.create.convert_id_to_node(tree_size, d440, cc_totss)
  
  if (debug_mode == TRUE) {
    #print("convert_id_to_node is")
    #print(convert_id_to_node)
    #print("dx_rx_types1_short is")
    #print(dx_rx_types1_short)
  }
  
  list_node_by_class <- c()
  for (i in 1:length(dx_rx_types1_short)) {
    name <- dx_rx_types1_short[i]
    subframe <- dxdf440_dataf[is.element(dxdf440_dataf[[classication]], dx_rx_types1_short[i]), ]
    
    if (nrow(subframe) == 0) {
      vec <- NA
      vec2 <- NA
    } else {
      if (id_tip_trim_flag == FALSE) {
        vec <- c(paste0(id_tip_prefix, subframe[[title.id]]))
      } else {
        vec <- c(paste0("", subframe[[title.id]]))
      }
      
      vec2 <- c()
      for (v in vec) {
        if (v == "IDroot") {
          vec2 <- c(vec2, "root")
        } else {
          vec2 <- c(vec2, v)
        }
      }
    }
    
    vec <- vec2
    vec1 <- vec
    
    if (is.na(vec[1]) && length(vec) == 1) {
      vec1 <- c(" ")
      list_node_by_class[[name]] <- as.numeric(vec1)
    } else {
      for (i in 1:length(vec)) {
        vec1[i] <- as.numeric(convert_id_to_node[vec[i]])
      }
      
      vec1 <- c(" ", vec1)
      list_node_by_class[[name]] <- as.numeric(sort(as.numeric(vec1)))
    }
  }
  
  if (debug_mode == TRUE) {
    #print("list_node_by_class is")
    #print(list_node_by_class)
  }
  
  return(list_node_by_class)
}



# Convert node ID to node number
func.create.convert_id_to_node <- function(tree_size, d440, cc_totss) {
  convert_node_to_id <- vector(mode = "list", length = tree_size)
  convert_id_to_node <- vector(mode = "list", length = tree_size)
  names(convert_node_to_id) <- c(d440$node)
  names(convert_id_to_node) <- c(d440$label)
  
  for (i in cc_totss) {
    convert_node_to_id[[i]] <- d440$label[i]
    convert_id_to_node[[i]] <- d440$node[i]
  }
  
  return(convert_id_to_node)
}


# Create list of nodes with bootstrap greater than 0.9
func.create.cc_nodss90 <- function(subframe_of_nodes) {
  subframe_of_nodes90 <- subframe_of_nodes[subframe_of_nodes$label >= 0.9, ]
  list_nodes90 <- list(subframe_of_nodes90$node)
  cc_nods90 <- c(list_nodes90)
  cc_nodss90 <- cc_nods90[[1]]
  
  temp <- c()
  for (el in cc_nodss90) {
    temp <- c(temp, el)
  }
  
  return(temp)
}


# Create list of nodes with bootstrap between 0.8 and 0.9
func.create.cc_nodss80 <- function(subframe_of_nodes) {
  subframe_of_nodes80 <- subframe_of_nodes[subframe_of_nodes$label >= 0.8 & subframe_of_nodes$label < 0.9, ]
  list_nodes80 <- list(subframe_of_nodes80$node)
  cc_nods80 <- c(list_nodes80)
  cc_nodss80 <- cc_nods80[[1]]
  
  return(cc_nodss80)
}

# Create list of nodes with bootstrap between 0.7 and 0.8
func.create.cc_nodss70 <- function(subframe_of_nodes) {
  subframe_of_nodes70 <- subframe_of_nodes[subframe_of_nodes$label >= 0.7 & subframe_of_nodes$label < 0.8, ]
  list_nodes70 <- list(subframe_of_nodes70$node)
  cc_nods70 <- c(list_nodes70)
  cc_nodss70 <- cc_nods70[[1]]
  
  return(cc_nodss70)
}



# Create new colors list for each node
func.create.new_colors_list <- function(FDR_perc, tree_TRY, tree_with_group, no_name, tree_size) {
  list_uncertain <- c()
  list_certain <- c() 
  list_equals <- c()
  new_colors_list <- rep(no_name, tree_size)
  
  for (nod in tree_TRY$data$node) {
    isT <- tree_TRY$data[[nod, "isTip"]]
    count_list <- c()
    
    if (isT == TRUE) {
      new_colors_list[nod] <- as.character.factor(tree_TRY$data[[nod, "group"]])
    } else {
      kids <- ggtree:::getSubtree(tree_with_group, nod)
      tips_list <- c()
      
      for (j in kids) {
        check <- tree_TRY$data[[j, "isTip"]]
        if (check == TRUE) {
          tips_list <- append(tips_list, j)
        }
      }
      
      for (j in tips_list) {
        color <- as.character.factor(tree_TRY$data[[j, "group"]])
        if (color %in% names(count_list)) {
          count_list[color] <- count_list[color] + 1
        } else {
          count_list[color] <- 1
        }
      }
      
      identical_test <- length(which(count_list == max(count_list)))
      if (identical_test == 1) { 
        color <- names(which(count_list == max(count_list)))
        new_colors_list[nod] <- color
        list_certain <- append(list_certain, nod)
      } else {
        list_uncertain <- append(list_uncertain, nod)
        vec <- names(which(count_list == max(count_list)))
        list_equals[[nod]] <- vec
      }
    }
  }
  
  for (nod in list_uncertain) {
    flag_end <- FALSE
    flag <- FALSE
    prev_par <- nod
    color <- no_name
    color_temp <- no_name
    vec <- list_equals[[nod]]
    
    while (flag_end == FALSE & flag == FALSE) {
      par <- tree_TRY$data[[prev_par, "parent"]]
      if (par == prev_par) {
        flag_end <- TRUE
      }
      if (par %in% list_certain) {
        flag <- TRUE
        color_temp <- new_colors_list[par]
        if (color_temp %in% vec) {
          color <- color_temp
        }
      } 
      prev_par <- par
    }
    
    new_colors_list[nod] <- color
  }
  
  return(new_colors_list)
}


# S2.0-PERF: Create hash for p_list_of_pairs caching (Option 3A)
# This hash is used to determine if expensive p-value calculations can be skipped.
# Inputs that affect p-values: tree structure, classification mapping, FDR, simulate.p.value
# IMPORTANT: Classification column changes MUST invalidate cache (user requirement)
func.create.p_list_cache_hash <- function(tree440, list_id_by_class, dx_rx_types1_short,
                                           FDR_perc, simulate.p.value) {
  # Create a hash from the key inputs that affect p_list_of_pairs calculation
  # Using digest package for consistent hashing

  # Build a list of all cache-relevant data
  cache_data <- list(
    # Tree structure - use tip labels and edge structure
    tree_tips = sort(tree440$tip.label),
    tree_Nnode = tree440$Nnode,
    tree_edge = tree440$edge,
    # Classification data - this is CRITICAL for cache invalidation
    classification_types = sort(dx_rx_types1_short),
    classification_mapping = lapply(list_id_by_class, function(x) sort(x)),
    # Statistical parameters
    FDR_perc = FDR_perc,
    simulate_p_value = simulate.p.value
  )

  # Use digest for reliable hashing
  hash_value <- digest::digest(cache_data, algo = "md5")

  return(hash_value)
}

# S2.9-PERF: Create hash for heatmap caching
# This hash determines if a heatmap can be reused from cache.
# Inputs that affect heatmap rendering: data, colors, dimensions, display settings
func.create.heatmap_cache_hash <- function(heat_param, heat_data = NULL, tree_tips = NULL) {
  # Build a list of all cache-relevant data for this heatmap
  cache_data <- list(
    # Data source and content
    data_source = heat_param[['data_source']],
    # For RData, include a hash of the actual data matrix
    data_hash = if (!is.null(heat_data)) digest::digest(heat_data, algo = "xxhash32") else NULL,
    # Tree tips affect y-positioning
    tree_tips_hash = if (!is.null(tree_tips)) digest::digest(sort(tree_tips), algo = "xxhash32") else NULL,
    # Color parameters
    low_color = heat_param[['low']],
    mid_color = heat_param[['mid']],
    high_color = heat_param[['high']],
    midpoint = heat_param[['midpoint']],
    na_color = heat_param[['na_color']],
    # Display settings
    cnv_display_mode = heat_param[['cnv_display_mode']],
    cnv_height_scale = heat_param[['cnv_height_scale']],
    cnv_render_downsample = heat_param[['cnv_render_downsample']],
    # Dimensions
    distance = heat_param[['distance']],
    height = heat_param[['height']],
    # Other visual settings
    show_col_lines = heat_param[['show_col_lines']],
    col_line_color = heat_param[['col_line_color']],
    col_line_size = heat_param[['col_line_size']],
    show_row_lines = heat_param[['show_row_lines']],
    row_line_color = heat_param[['row_line_color']],
    row_line_size = heat_param[['row_line_size']]
  )

  # Use digest for reliable hashing (xxhash64 is fast)
  hash_value <- digest::digest(cache_data, algo = "xxhash64")

  return(hash_value)
}


# Create list of p-values for each pair of nodes
func.create.p_list_of_pairs <- function(list_node_by_class, d440, dx_rx_types1_short,
                                        cc_nodss, tree_with_group, FDR_perc, tree, cc_tipss,
                                        tree_TRY, tree_size, no_name, simulate.p.value) {
  
  in_sub_tree <- rep(0, 2)
  Not_in_sub_tree <- rep(0, 2)
  p_list_of_pairs <- rep(1, tree_size)
  
  for (nod in cc_nodss) {   
    kids <- ggtree:::getSubtree(tree_with_group, nod)  # kids contain all the nodes in the subtree including nod (its root)
    sub_tree_size <- length(kids)
    tips_list <- c()
    
    in_sub_tree <- rep(0, 2)
    Not_in_sub_tree <- rep(0, 2)
    
    hyp <- tree_TRY$data[nod, "new_class"]
    if (hyp != no_name) {
      for (j in kids) {
        check <- tree_TRY$data[[j, "isTip"]]
        if (check == TRUE) {
          tips_list <- append(tips_list, j)
        }
      }
      
      for (j in cc_tipss) {
        j_typ <- tree_TRY$data[j, "new_class"]
        if (j %in% tips_list) {
          if (j_typ == hyp) {
            in_sub_tree[1] <- in_sub_tree[1] + 1
          } else {
            in_sub_tree[2] <- in_sub_tree[2] + 1
          }
        } else {
          if (j_typ == hyp) {
            Not_in_sub_tree[1] <- Not_in_sub_tree[1] + 1
          } else {
            Not_in_sub_tree[2] <- Not_in_sub_tree[2] + 1
          }
        }
      }
      
      Not_in_sub_tree1 <- abs(Not_in_sub_tree)
      in_sub_tree1 <- abs(in_sub_tree)
      
      xtab <- as.table(rbind(
        in_sub_tree1,
        Not_in_sub_tree1
      ))
      
      dimnames(xtab) <- list(
        In_sub_tree = c("Yes", "No"),
        Groups = c("In_type", "Not_in_type")
      )
      
      conf_comp <- 1 - FDR_perc
      
      test <- fisher.test(xtab, conf.level = conf_comp, simulate.p.value = simulate.p.value)
      p_val <- test$p.value
      p_list_of_pairs[nod] <- p_val  
    }
  }
  
  return(p_list_of_pairs)
}


# Create p-value list from list of pairs
func.create.p_val_list_FROM_LIST <- function(FDR_perc, tree_TRY, p_list_of_pairs, op_list) {
  p_val_list <- c()
  
  for (nod in tree_TRY$data$node) {
    p_val_op <- 0
    p_val_precise <- p_list_of_pairs[nod]
    
    if (FDR_perc < p_val_precise) {
      p_val_op <- op_list[1]
    } 
    if (FDR_perc >= p_val_precise && p_val_precise > 0.5 * FDR_perc) {
      p_val_op <- op_list[2]
    }
    if (0.5 * FDR_perc >= p_val_precise) {
      p_val_op <- op_list[3]
    }
    
    p_val_list <- append(p_val_list, p_val_op)
  }
  
  return(p_val_list)
}



# Create rotation parameters
func.make.rot.params <- function(rot_index, yaml_file, title.id, ids_list, tree440, readfile,
                                 id_tip_trim_flag, id_tip_prefix) {
  # Rotation function
  rotation_name <- paste0('rotation', rot_index)
  if (rot_index == 1) {
    rot_class_list <- yaml_file[['visual definitions']]$rotation1$'according'
  } else {
    rot_class_list <- yaml_file[['visual definitions']]$rotation2$'according'
  }
  
  how_many_rot <- length(rot_class_list)
  column_base <- 3
  
  list_id_by_rot <- func.make.mapping.list(how_many_rot, rot_class_list, yaml_file, title.id,
                                           ids_list, readfile, FALSE, NA,
                                           id_tip_trim_flag, id_tip_prefix)
  
  TREE_OTU_dx.rx <- groupOTU(tree440, list_id_by_rot)
  # v56b: Suppress harmless fortify warnings
  iis <- subset(suppressWarnings(ggtree(TREE_OTU_dx.rx))$data, group == 0)
  types_list_dx.rx <- names(list_id_by_rot)
  list_weight_dx.rx <- func.set.list_weight(types_list_dx.rx)
  
  ret_rot <- c()
  ret_rot[['list_id_by_rot']] <- list_id_by_rot
  ret_rot[['TREE_OTU_dx.rx']] <- TREE_OTU_dx.rx
  ret_rot[['types_list_dx.rx']] <- types_list_dx.rx
  ret_rot[['list_weight_dx.rx']] <- list_weight_dx.rx
  
  return(ret_rot)
}


# Create mapping list
func.make.mapping.list <- function(cls_num, cls_list, yaml_file, title.id, ids_list, readfile,
                                   debug_mode = FALSE, no_nam = NA,
                                   id_tip_trim_flag = FALSE, id_tip_prefix = "AAA") {
  if (debug_mode == TRUE) {
    # v53: print("In func.make.mapping.list")
  }
  
  nams <- c()
  for (i in names(cls_list)) {
    value <- cls_list[[i]]$'value'
    nams <- c(nams, value)
  }
  
  class1 <- c()
  for (ty in 1:length(nams)) {
    name <- nams[ty]
    
    # Safety check: skip if name is NULL, NA, or empty
    if (is.null(name) || is.na(name) || length(name) == 0 || name == "") {
      if (debug_mode == TRUE) {
        # v53: print(paste0("Skipping empty name at index ", ty))
      }
      next
    }
    
    vec_num <- c("")
    vec <- c("")
    
    for (ii in ids_list) {
      if (ii == "root") {
        i <- "root"
      } else {
        if (nchar(ii) > 5) {
          i <- substring(ii, 3, nchar(ii))
        } else {
          i <- ii
        }
      }
      
      check <- func.check.if.id.in.class.NEW(name, i, yaml_file, readfile, title.id, cls_list) 
      
      if (check == TRUE) {
        if (i == "root") {
          vec <- c(vec, i)
        } else {
          if (id_tip_trim_flag == FALSE) {
            vec <- c(vec, paste0(id_tip_prefix, i))
          } else {
            vec <- c(vec, paste0("", i))
          }
        }
      }
    }
    
    # Only assign if name is valid and not empty
    if (!is.null(name) && !is.na(name) && length(name) > 0 && name != "") {
      class1[[name]] <- vec
    } else {
      if (debug_mode == TRUE) {
        # v53: print(paste0("Warning: Attempted to assign with invalid name: ", name))
      }
    }
  }
  
  return(class1)
}


# Check if ID is in class
func.check.if.id.in.class.NEW <- function(class_name, id, yaml_file, readfile, title.id, cls_list) {
  lin <- readfile[which(readfile[title.id] == id), ]
  
  ans <- FALSE
  for (i in names(cls_list)) {
    if ('title' %in% names(cls_list[[i]])) {
      title <- cls_list[[i]]$'title'
      value <- cls_list[[i]]$'value'
      
      if (title %in% colnames(lin)) {
        val <- lin[[title]]
        if (length(val) > 0) {
          if (val == class_name) {
            ans <- TRUE
          }
        }
      }
    } else {
      ix <- 1
      flag <- TRUE
      ans <- TRUE
      
      while (flag == TRUE) {
        ti <- paste0("title", ix)
        vi <- paste0("value", ix)
        
        if (ti %in% names(cls_list[[i]])) {
          title <- cls_list[[i]][[ti]]
          value <- cls_list[[i]][[vi]]
          
          if (title %in% colnames(lin)) {
            val <- lin[[title]]
            if (length(val) > 0) {
              if (val %in% value & ans == TRUE) {
                ans <- TRUE
              }
            }
          }
          
          ix <- ix + 1
        } else {
          flag <- FALSE
          if (ix == 1) {
            ans <- FALSE
          }
        }
      }
    }
  }
  
  return(ans)
}



# Create weight list for nodes
func.set.list_weight <- function(types_list) {
  list.weight <- c()
  start <- 1
  jmp <- length(types_list) + 3
  pow <- 1
  weight <- start
  
  for (type in types_list) {
    list.weight[type] <- weight
    weight <- weight + 2 * ((jmp^pow))
    pow <- pow + 1
  }
  
  return(list.weight)
}


# Create weight list for tree
func.create.weight_list <- function(tree_TRY, weights_list, TREE_OTU_class, tree_size) {
  # Define weights to tips for rotation
  # v56b: Suppress harmless fortify warnings
  g <- suppressWarnings(ggtree(TREE_OTU_class))
  weight_list <- rep(0, tree_size)
  
  for (nod in g$data$node) {
    kids <- ggtree:::getSubtree(TREE_OTU_class, nod)
    tips_list <- c()
    weight <- 0
    s <- 0
    
    for (j in kids) {
      check <- tree_TRY$data[[j, "isTip"]]
      if (check == TRUE) {
        tips_list <- append(tips_list, j)
      }
    }
    
    sl <- c()
    for (j in tips_list) {
      typ <- as.character.factor(g$data[[j, "group"]])
      weight_tip <- 0
      
      if (typ %in% names(weights_list)) {
        weight_tip <- weights_list[[typ]]
      }
      
      s <- s + weight_tip
      sl <- c(sl, weight_tip)
    }
    
    weight_op1 <- min(sl)
    weight_op2 <- func.get.most.common.element.in.list(sl)
    max_weight <- max(weights_list)
    weight_op3 <- func.get.normalized.price(sl, max_weight)
    weight <- weight_op3
    weight_list[nod] <- weight * 100000
  }
  
  return(weight_list)
}



# Get most common element in list
func.get.most.common.element.in.list <- function(l) {
  return(as.numeric(names(sort(table(l), decreasing = TRUE)[1])))
}

# Get normalized price
func.get.normalized.price <- function(l, max_weight) {
  val <- sum(l) / (max_weight * length(l))
  return(val)
}

# Find plot boundaries
func.find.plot.boundaries <- function(tt, debug_mode = FALSE) {
  buildt <- ggplot_build(tt)
  
  mmin <- 0
  mmax <- 0
  build_length <- length(buildt$data)
  
  for (index_build in 1:build_length) {
    build_frame <- buildt$data[[index_build]]
    if ('y' %in% colnames(buildt$data[[index_build]])) {
      mmin <- min(mmin, min(build_frame$y))
      mmax <- max(mmax, max(build_frame$y))
    }
    if ('yend' %in% colnames(buildt$data[[index_build]])) {
      mmin <- min(mmin, min(build_frame$yend))
      mmax <- max(mmax, max(build_frame$yend))
    }
  }
  
  xmmin <- 0
  xmmax <- 0
  
  for (index_build in 1:build_length) {
    build_frame <- buildt$data[[index_build]]
    if ('x' %in% colnames(buildt$data[[index_build]])) {
      xmmin <- min(xmmin, min(build_frame$x))
      xmmax <- max(xmmax, max(build_frame$x))
    }
    if ('yend' %in% colnames(buildt$data[[index_build]])) {
      xmmin <- min(xmmin, min(build_frame$xend))
      xmmax <- max(xmmax, max(build_frame$xend))
    }
  }
  
  if (debug_mode == TRUE) {
    # v53: print("xmmin is")
    # v53: print(xmmin)
    # v53: print("xmmax is")
    # v53: print(xmmax)
  }
  
  plot_bounds <- c()
  plot_bounds$xmin <- xmmin
  plot_bounds$xmax <- xmmax
  plot_bounds$ymin <- mmin
  plot_bounds$ymax <- mmax
  
  if (debug_mode == TRUE) {
    # v53: print("plot_bounds is")
    # v53: print(plot_bounds)
  }
  
  return(plot_bounds)
}


##### part 3

# ============================================================================
# TREE VISUALIZATION FUNCTIONS
# ============================================================================

# Function to prepare leaves ID ordered for dataframe
func.make.leaves_id_ordered_for_df440 <- function(leaves_id_from_tree1, dxdf440_dataf, title.id, id_tip_trim_flag, id_tip_prefix) {
  leaves_id_ordered_for_df440 <- c()
  base.list <- dxdf440_dataf[[title.id]]
  tips_list <- as.list(leaves_id_from_tree1)
  
  for (id in base.list) {
    if (id_tip_trim_flag == "FALSE") {
      temp <- paste0(id_tip_prefix, id)
    } else {
      temp <- id
    }
    
    if (temp %in% tips_list$label) {
      leaves_id_ordered_for_df440 <- c(leaves_id_ordered_for_df440, temp)
    }
  }
  
  return(leaves_id_ordered_for_df440)
}


# Function to handle highlighting
# v139: Added high_alpha_list parameter for transparency control
func_highlight <- function(p, how_many_hi, heat_flag, high_color_list, a, b, man_adjust_elipse, pr440_short_tips_TRY,
                           boudariestt, debug_mode = FALSE, high_offset = 0, high_vertical_offset = 0,
                           high_alpha_list = NULL) {
  # DEBUG-2ND-HIGHLIGHT: Track entry with details
  cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] ENTER func_highlight at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
  cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   how_many_hi=", how_many_hi, ", heat_flag=", heat_flag, "\n"))
  cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   p$layers count=", length(p$layers), "\n"))

  up_offset <- -1 # -3
  y_off_base <- -8

  # v139: Default alpha list if not provided
  if (is.null(high_alpha_list)) {
    high_alpha_list <- rep(0.5, how_many_hi)
  }

  # Debug output for Bug #11
  # v53: cat(file=stderr(), "🔵 func_highlight ENTRY:\n")
  # v53: cat(file=stderr(), paste0("🔵   how_many_hi: ", how_many_hi, "\n"))
  # v53: cat(file=stderr(), paste0("🔵   heat_flag: ", heat_flag, "\n"))
  # v53: cat(file=stderr(), paste0("🔵   a (ellipse height): ", a, "\n"))
  # v53: cat(file=stderr(), paste0("🔵   b (ellipse width): ", b, "\n"))
  # v53: cat(file=stderr(), paste0("🔵   high_offset (horizontal): ", high_offset, "\n"))
  # v53: cat(file=stderr(), paste0("🔵   high_vertical_offset (vertical): ", high_vertical_offset, "\n"))
  # v53: cat(file=stderr(), paste0("🔵   high_color_list: ", paste(unlist(high_color_list), collapse=", "), "\n"))
  # v53: cat(file=stderr(), paste0("🔵   Columns in pr440_short_tips_TRY$data: ", paste(names(pr440_short_tips_TRY$data), collapse=", "), "\n"))
  # v53: cat(file=stderr(), paste0("🔵   'high1' column exists: ", "high1" %in% names(pr440_short_tips_TRY$data), "\n"))
  
  for (index_high in 1:how_many_hi) {
    if (index_high == 1) {
      # Check if high1 column exists before trying to filter
      if (!"high1" %in% names(pr440_short_tips_TRY$data)) {
        # v53: cat(file=stderr(), "🔴 ERROR: 'high1' column does NOT exist in pr440_short_tips_TRY$data!\n")
        # v53: cat(file=stderr(), paste0("🔴 Available columns: ", paste(names(pr440_short_tips_TRY$data), collapse=", "), "\n"))
        next
      }

      # v146: When heat_flag is TRUE, get node positions from p$data (after heatmap transform)
      # to ensure ellipses align with the transformed tree coordinates
      if (heat_flag == TRUE && "high1" %in% names(p$data)) {
        high_nodes_table1 <- p$data[p$data$high1 == TRUE,]
        debug_cat(paste0("  v46: Using p$data for ellipse positioning (heat_flag=TRUE)\n"))
        debug_cat(paste0("  v46: high_nodes_table1 rows: ", nrow(high_nodes_table1), "\n"))
        if (nrow(high_nodes_table1) > 0) {
          debug_cat(paste0("  v46: y-range: ", round(min(high_nodes_table1$y, na.rm=TRUE), 2),
                                    " to ", round(max(high_nodes_table1$y, na.rm=TRUE), 2), "\n"))
        }
      } else {
        high_nodes_table1 <- pr440_short_tips_TRY$data[pr440_short_tips_TRY$data$high1 == TRUE,]
      }
      # v53: cat(file=stderr(), paste0("🔵   high_nodes_table1 rows: ", nrow(high_nodes_table1), "\n"))
      
      if (debug_mode == TRUE) {
        # v53: print("high_nodes_table1 is")
        # v53: print(high_nodes_table1)
        # v53: print("high_nodes_table1$x is")
        # v53: print(high_nodes_table1$x)
        # v53: print("max high_nodes_table1$x is")
        # v53: print(max(pr440_short_tips_TRY$data[,'x']))
      }
      
      if (length(high_nodes_table1$x) == 0) {
        # v53: print(paste0(paste0("Highlist number ", index_high), "- None of the nodes is highlighted, skipped"))
        next
      }
      
      man_adjust_elipse <- man_adjust_elipse + high_offset
      bas <- max(pr440_short_tips_TRY$data[,'x']) - as.integer(high_nodes_table1$x)
      
      # Bug #12 fix: Removed hardcoded -15 that pushed ellipses too far left
      # Bug #14 fix: Inverted offset sign so positive moves RIGHT (intuitive)
      # v53: cat(file=stderr(), paste0("🔵 Ellipse positioning: man_adjust_elipse=", man_adjust_elipse, 
      #                            " (inverted), max_x=", max(pr440_short_tips_TRY$data[,'x']), "\n"))
      
      # v139: Use high_alpha_list for transparency instead of hardcoded 0.5
      alpha_val <- if (length(high_alpha_list) >= 1 && !is.null(high_alpha_list[[1]])) high_alpha_list[[1]] else 0.5

      # v148: Debug output for PLOT ellipse alpha and dimensions
      debug_cat(paste0("\n=== v148: PLOT ELLIPSE (high1) ===\n"))
      debug_cat(paste0("  Alpha: ", alpha_val, "\n"))
      debug_cat(paste0("  Dimensions: a=", round(a, 4), ", b=", round(b, 4), "\n"))
      debug_cat(paste0("==================================\n"))

      if (heat_flag == FALSE) {
        p <- p +
          geom_ellipse(data = high_nodes_table1,
                       aes(x0 = ((max(pr440_short_tips_TRY$data[,'x']) - x) * (-1) - man_adjust_elipse),
                           y0 = y + high_vertical_offset, a = a, b = b, angle = 0),
                       fill = high_color_list[[1]], alpha = alpha_val, linetype = "blank", show.legend = FALSE)
      } else {
        # v148: Fixed ellipse x0 positioning for heatmap mode
        # CRITICAL: Must use the SAME data source for both tree_max_x and node x-values
        # When heat_flag=TRUE, we get nodes from p$data, so we must also get tree_max from p$data
        # Using different data sources causes coordinate mismatch and wrong positioning

        # Get tree_max_x from the SAME source as node positions (p$data)
        tree_max_x_from_p <- max(p$data[p$data$isTip == TRUE, 'x'], na.rm = TRUE)

        # Also get the original tree_max for comparison
        tree_max_x_original <- max(pr440_short_tips_TRY$data[,'x'], na.rm = TRUE)

        debug_cat(paste0("\n=== v148: ELLIPSE X0 POSITIONING (heat mode) ===\n"))
        debug_cat(paste0("  tree_max_x (from p$data tips): ", round(tree_max_x_from_p, 4), "\n"))
        debug_cat(paste0("  tree_max_x (from pr440 - for reference): ", round(tree_max_x_original, 4), "\n"))
        debug_cat(paste0("  man_adjust_elipse: ", man_adjust_elipse, "\n"))
        debug_cat(paste0("  Sample node x values: ", paste(round(head(high_nodes_table1$x, 3), 4), collapse=", "), "\n"))
        debug_cat(paste0("  Sample calculated x0: ", paste(round((tree_max_x_from_p - head(high_nodes_table1$x, 3)) * (-1) - man_adjust_elipse, 4), collapse=", "), "\n"))
        debug_cat(paste0("=============================================\n"))

        p <- p +
          geom_ellipse(data = high_nodes_table1,
                       aes(x0 = ((tree_max_x_from_p - x) * (-1) - man_adjust_elipse),
                           y0 = y + high_vertical_offset, a = a, b = b, angle = 0),
                       fill = high_color_list[[1]], alpha = alpha_val, linetype = "blank", show.legend = FALSE)
      }
    } else if (index_high == 2) {
      high_nodes_table2 <- p$data[p$data$high2 == TRUE,]
      # v139: Use high_alpha_list for transparency
      alpha_val2 <- if (length(high_alpha_list) >= 2 && !is.null(high_alpha_list[[2]])) high_alpha_list[[2]] else 0.5

      # S1.2: Fixed - use same positioning logic as highlight 1 (was using undefined x_range_min)
      # Get tree_max_x from p$data for consistent positioning
      tree_max_x <- max(p$data[p$data$isTip == TRUE, 'x'], na.rm = TRUE)

      p <- p +
        geom_ellipse(data = high_nodes_table2,
                     aes(x0 = ((tree_max_x - x) * (-1) - man_adjust_elipse),
                         y0 = y + high_vertical_offset, a = a, b = b, angle = 0),
                     fill = high_color_list[[2]], alpha = alpha_val2, linetype = "blank", show.legend = FALSE)
    } else if (index_high == 3) {
      high_nodes_table3 <- p$data[p$data$high3 == TRUE,]
      # v139: Use high_alpha_list for transparency
      alpha_val3 <- if (length(high_alpha_list) >= 3 && !is.null(high_alpha_list[[3]])) high_alpha_list[[3]] else 0.5

      # S1.2: Fixed - use same positioning logic as highlight 1 (was using undefined x_range_min)
      # Get tree_max_x from p$data for consistent positioning
      tree_max_x <- max(p$data[p$data$isTip == TRUE, 'x'], na.rm = TRUE)

      p <- p +
        geom_ellipse(data = high_nodes_table3,
                     aes(x0 = ((tree_max_x - x) * (-1) - man_adjust_elipse),
                         y0 = y + high_vertical_offset, a = a, b = b, angle = 0),
                     fill = high_color_list[[3]], alpha = alpha_val3, linetype = "blank", show.legend = FALSE)
    }
  }

  # v180: REMOVED legacy layer reordering code that was causing heatmap corruption
  # The old code (Bug #11 fix) assumed a specific layer order and would scramble layers
  # when there were 8+ layers. With multiple heatmaps (3+), this caused the first heatmap
  # to disappear. Layer ordering is now handled by func.move.tiplabels.to.front() at the
  # END of func.print.lineage.tree (line ~7065), not here after each ellipse addition.

  # S1.3-PERF: Removed intermediate call to func.move.tiplabels.to.front() here
  # Layer reordering is now done ONCE at the end in generate_plot()
  # The "second highlight stuck" bug was caused by undefined x_range_min (fixed in S1.2),
  # not by missing layer reordering calls. This saves ~2-3 redundant function calls per render.

  cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] EXIT func_highlight at ", format(Sys.time(), "%H:%M:%OS3"), ", layers=", length(p$layers), "\n"))
  return(p)
}
# v160: Helper function to add custom legends to gtable's legend area
# NOTE: This function is NOT CURRENTLY USED (v161 reverted to annotation-based approach)
# Kept for potential future use if gtable approach is revisited
# This places highlight and bootstrap legends alongside ggplot legends (Classification, Heatmap)
func.add.custom.legends.to.gtable <- function(gt, legend_info) {
  if (is.null(legend_info)) {
    return(gt)
  }

  debug_cat(paste0("\n=== v160: ADDING CUSTOM LEGENDS TO GTABLE ===\n"))

  # Find the guide-box that has actual content
  guide_box_names <- c("guide-box-right", "guide-box-bottom", "guide-box-left", "guide-box-top", "guide-box")
  guide_box_idx <- NULL
  guide_box_name <- NULL

  for (gb_name in guide_box_names) {
    idx <- which(gt$layout$name == gb_name)
    if (length(idx) > 0) {
      gb_row <- gt$layout$t[idx]
      gb_col <- gt$layout$l[idx]
      cell_width <- gt$widths[gb_col]
      width_cm <- tryCatch(grid::convertWidth(cell_width, "cm", valueOnly = TRUE), error = function(e) 0)

      debug_cat(paste0("  Checking ", gb_name, " at row ", gb_row, ", col ", gb_col,
                                 " (width=", round(width_cm, 2), "cm)\n"))

      if (width_cm > 0.1) {
        guide_box_idx <- idx
        guide_box_name <- gb_name
        break
      }
    }
  }

  if (is.null(guide_box_idx)) {
    guide_box_idx <- which(gt$layout$name == "guide-box-right")[1]
    if (length(guide_box_idx) == 0) {
      debug_cat(paste0("  ERROR: No guide-box found, skipping custom legends\n"))
      return(gt)
    }
    guide_box_name <- "guide-box-right"
  }

  guide_row <- gt$layout$t[guide_box_idx]
  guide_col <- gt$layout$l[guide_box_idx]
  debug_cat(paste0("  Using ", guide_box_name, " at row ", guide_row, ", col ", guide_col, "\n"))

  # Calculate legend content
  n_highlight <- 0
  n_bootstrap <- 0
  if (!is.null(legend_info$highlight) && legend_info$show_highlight) {
    n_highlight <- length(legend_info$highlight$labels)
  }
  if (!is.null(legend_info$bootstrap) && legend_info$show_bootstrap) {
    n_bootstrap <- 3
  }

  if (n_highlight == 0 && n_bootstrap == 0) {
    debug_cat(paste0("  No custom legends to add\n"))
    return(gt)
  }

  debug_cat(paste0("  Highlight: ", n_highlight, " items, Bootstrap: ", n_bootstrap, " items\n"))

  # v160: Build legend as a proper nested gtable with rows for each item
  title_fontsize <- 12
  text_fontsize <- 10
  row_height <- grid::unit(0.5, "cm")
  key_width <- grid::unit(0.8, "cm")
  label_width <- grid::unit(2.2, "cm")

  # Calculate rows needed
  n_rows <- 0
  if (n_highlight > 0) n_rows <- n_rows + 1 + n_highlight  # title + items
  if (n_bootstrap > 0) n_rows <- n_rows + 1 + 3  # title + 3 bootstrap levels
  if (n_highlight > 0 && n_bootstrap > 0) n_rows <- n_rows + 1  # gap

  # Create heights using unit.c to properly combine units
  heights_list <- lapply(1:n_rows, function(i) row_height)
  heights <- do.call(grid::unit.c, heights_list)
  widths <- grid::unit.c(key_width, label_width)

  # Create the legend gtable
  legend_gt <- gtable::gtable(widths = widths, heights = heights)

  current_row <- 1

  # Add highlight legend
  if (!is.null(legend_info$highlight) && legend_info$show_highlight) {
    hi <- legend_info$highlight

    # Title (spans both columns)
    legend_gt <- gtable::gtable_add_grob(
      legend_gt,
      grobs = grid::textGrob(hi$title, x = 0, hjust = 0,
                             gp = grid::gpar(fontsize = title_fontsize, fontface = "bold")),
      t = current_row, l = 1, r = 2
    )
    current_row <- current_row + 1

    # Each highlight item
    for (i in seq_along(hi$labels)) {
      # Key (ellipse)
      theta <- seq(0, 2*pi, length.out = 30)
      legend_gt <- gtable::gtable_add_grob(
        legend_gt,
        grobs = grid::polygonGrob(
          x = grid::unit(0.5, "npc") + grid::unit(6 * cos(theta), "pt"),
          y = grid::unit(0.5, "npc") + grid::unit(4 * sin(theta), "pt"),
          gp = grid::gpar(fill = hi$colors[[i]], alpha = hi$alphas[[i]], col = NA)
        ),
        t = current_row, l = 1
      )

      # Label
      legend_gt <- gtable::gtable_add_grob(
        legend_gt,
        grobs = grid::textGrob(hi$labels[[i]], x = 0, hjust = 0,
                               gp = grid::gpar(fontsize = text_fontsize)),
        t = current_row, l = 2
      )
      current_row <- current_row + 1
    }

    # Gap before bootstrap
    if (n_bootstrap > 0) {
      current_row <- current_row + 1
    }
  }

  # Add bootstrap legend
  if (!is.null(legend_info$bootstrap) && legend_info$show_bootstrap) {
    # Title
    legend_gt <- gtable::gtable_add_grob(
      legend_gt,
      grobs = grid::textGrob("Bootstrap", x = 0, hjust = 0,
                             gp = grid::gpar(fontsize = title_fontsize, fontface = "bold")),
      t = current_row, l = 1, r = 2
    )
    current_row <- current_row + 1

    tri_labels <- c(">90%", ">80%", ">70%")
    tri_sizes <- c(5, 4, 3)

    for (i in 1:3) {
      sz <- tri_sizes[i]
      # Triangle key
      legend_gt <- gtable::gtable_add_grob(
        legend_gt,
        grobs = grid::polygonGrob(
          x = grid::unit(0.5, "npc") + grid::unit(c(-sz, sz, 0), "pt"),
          y = grid::unit(0.5, "npc") + grid::unit(c(-sz*0.6, -sz*0.6, sz*0.8), "pt"),
          gp = grid::gpar(fill = "grey36", col = "grey20", alpha = 0.5)
        ),
        t = current_row, l = 1
      )

      # Label
      legend_gt <- gtable::gtable_add_grob(
        legend_gt,
        grobs = grid::textGrob(tri_labels[i], x = 0, hjust = 0,
                               gp = grid::gpar(fontsize = text_fontsize)),
        t = current_row, l = 2
      )
      current_row <- current_row + 1
    }
  }

  # Calculate total height
  total_height <- grid::unit(n_rows * 0.5, "cm")
  debug_cat(paste0("  Legend has ", n_rows, " rows, height=", n_rows * 0.5, "cm\n"))

  # v160: NEW APPROACH - Add legend row to MAIN gtable directly above guide-box
  # This ensures visibility since it's at the top-level gtable, not nested inside

  # Add a new row to the main gtable ABOVE the guide-box row
  gt <- gtable::gtable_add_rows(gt, heights = total_height, pos = guide_row - 1)

  # All grobs below the insertion point shift down by 1 row, so update guide_row
  # The new row is now at position guide_row, and old content is at guide_row + 1

  # Add our legend to the new row, in the same column as the guide-box
  gt <- gtable::gtable_add_grob(
    gt,
    grobs = legend_gt,
    t = guide_row,
    l = guide_col,
    name = "custom-legend"
  )

  debug_cat(paste0("  v60: Added legend as new row at position ", guide_row, " in main gtable\n"))
  debug_cat(paste0("  v60: Legend placed in column ", guide_col, " (same as guide-box)\n"))
  debug_cat(paste0("===============================================\n"))

  return(gt)
}


# Function to create the second legend
# v145: Added high_alpha_list parameter for transparency control
# v167: OPTION C - Using native ggplot legends instead of gtable manipulation
func.make.second.legend <- function(p, FLAG_BULK_DISPLAY, how_many_hi, heat_flag, how_many_boxes,
                                    how_mant_rows, boudariestt, y_off_base, high_title_list,
                                    size_font_legend_title, high_label_list, size_font_legend_text,
                                    high_color_list, a, b, x_range_min, show_boot_flag, size_90,
                                    size_80, size_70, man_adjust_image_of_second_legend,
                                    man_multiply_second_legend, man_multiply_second_legend_text,
                                    man_multiply_elipse, man_space_second_legend,
                                    man_space_second_legend_multiplier, man_offset_for_highlight_legend_x,
                                    debug_mode = FALSE, boot_values = NA, man_offset_second_legend = 0, width,
                                    bootstrap_label_size = 1.5,
                                    highlight_x_offset = 0, highlight_y_offset = 0,
                                    highlight_title_size = NULL, highlight_text_size = NULL,
                                    highlight_title_gap = 1, highlight_label_gap = 0.5,
                                    bootstrap_x_offset = 0, bootstrap_y_offset = 0,
                                    bootstrap_title_x_offset = 2,
                                    bootstrap_title_size_mult = NULL, bootstrap_text_size_mult = NULL,
                                    bootstrap_title_gap = 2, bootstrap_label_gap = 2,
                                    show_highlight_legend = TRUE, show_bootstrap_legend = TRUE,
                                    high_alpha_list = NULL) {

  # v178: OPTION C - NATIVE GGPLOT LEGENDS (SHAPE-ONLY approach)
  # Use ONLY shape aesthetic to avoid conflicts with existing scales:
  # - fill: used by heatmaps
  # - size: used by P value legend
  # - colour: used by classification
  # By using only shape (not used elsewhere), we avoid ALL scale conflicts
  # v178: Return nullGrob if shape value not in our set - prevents legend bleeding
  #       This stops ellipses/triangles from appearing on other legends

  debug_cat(paste0("\n=== v178: NATIVE GGPLOT LEGENDS - SHAPE-ONLY APPROACH ===\n"))
  debug_cat(paste0("  Using ONLY shape aesthetic to avoid scale conflicts\n"))
  debug_cat(paste0("  v78: key_glyph returns nullGrob for shapes not in valid set\n"))

  # Initialize high_alpha_list if NULL
  if (is.null(high_alpha_list) || length(high_alpha_list) == 0) {
    high_alpha_list <- rep(0.5, how_many_hi)
  }

  # Title fontsize should match ggplot legend titles
  title_fontsize <- if (!is.null(highlight_title_size)) highlight_title_size else size_font_legend_title
  text_fontsize <- if (!is.null(highlight_text_size)) highlight_text_size else size_font_legend_text
  boot_title_fontsize <- if (!is.null(bootstrap_title_size_mult)) bootstrap_title_size_mult else size_font_legend_title
  boot_text_fontsize <- if (!is.null(bootstrap_text_size_mult)) bootstrap_text_size_mult else size_font_legend_text

  debug_cat(paste0("  v78: Highlight title fontsize: ", title_fontsize, "\n"))
  debug_cat(paste0("  v78: Bootstrap title fontsize: ", boot_title_fontsize, "\n"))

  # ============================================
  # v178: HIGHLIGHT LEGEND using shape aesthetic
  # Custom key_glyph draws ellipses with correct colors/transparency
  # v178: Returns nullGrob for shapes not in valid set (prevents bleeding)
  # ============================================
  if (FLAG_BULK_DISPLAY == TRUE && show_highlight_legend == TRUE && how_many_hi > 0 &&
      !is.null(high_label_list) && length(high_label_list) > 0) {

    # Get highlight title
    highlight_title <- if (!is.null(high_title_list) && length(high_title_list) > 0) {
      high_title_list[[1]]
    } else {
      "Highlight"
    }

    debug_cat(paste0("\n  v178: Creating HIGHLIGHT legend (shape-only)\n"))
    debug_cat(paste0("    Title: '", highlight_title, "'\n"))
    debug_cat(paste0("    Items: ", length(high_label_list), "\n"))

    # v177: Create unique shape values for each item (1, 2, 3, ...)
    n_highlights <- length(high_label_list)
    highlight_shape_values <- seq_len(n_highlights)
    highlight_labels <- unlist(high_label_list)

    # Create data frame for highlight legend
    highlight_legend_data <- data.frame(
      x = rep(NA_real_, n_highlights),
      y = rep(NA_real_, n_highlights),
      highlight_shape = factor(highlight_labels, levels = highlight_labels),
      stringsAsFactors = FALSE
    )

    # v177: Create mappings from shape value to label, color, alpha
    # shape_to_label: 1 -> "tumor", 2 -> "normal", etc.
    shape_to_label <- setNames(highlight_labels, as.character(highlight_shape_values))
    shape_to_color <- setNames(unlist(high_color_list), as.character(highlight_shape_values))
    shape_to_alpha <- setNames(unlist(high_alpha_list), as.character(highlight_shape_values))

    for (i in seq_along(high_label_list)) {
      debug_cat(paste0("    Item ", i, ": shape=", i, " label='", high_label_list[[i]],
                                 "' color=", high_color_list[[i]],
                                 " alpha=", high_alpha_list[[i]], "\n"))
    }

    # v178: Custom key_glyph for ELLIPSE - draws ellipse with stored color/alpha
    # Uses shape value (1, 2, 3...) to look up the correct color and alpha
    # v178: Returns nullGrob if shape value not in our set (prevents legend bleeding)
    draw_key_highlight_ellipse <- local({
      shape_to_color_local <- shape_to_color
      shape_to_alpha_local <- shape_to_alpha
      valid_shapes <- names(shape_to_color)  # v178: Track valid shape values
      function(data, params, size) {
        # v178: Get the shape value and use it to look up color/alpha
        shape_val <- data$shape
        debug_cat(paste0("    v178: draw_key_highlight_ellipse called, shape=", shape_val, "\n"))

        if (is.null(shape_val) || length(shape_val) == 0) {
          debug_cat(paste0("    v178: shape is NULL/empty, returning nullGrob\n"))
          return(grid::nullGrob())
        }

        shape_key <- as.character(shape_val)

        # v178: CRITICAL - Return nullGrob if this shape is not one of ours
        # This prevents ellipses from appearing on other legends (classification, heatmap, etc.)
        if (!(shape_key %in% valid_shapes)) {
          debug_cat(paste0("    v178: shape ", shape_key, " not in valid set (", paste(valid_shapes, collapse=","), "), returning nullGrob\n"))
          return(grid::nullGrob())
        }

        fill_color <- shape_to_color_local[shape_key]
        fill_alpha <- shape_to_alpha_local[shape_key]

        debug_cat(paste0("    v178: Looked up color='", fill_color, "', alpha='", fill_alpha, "'\n"))

        # v178: Extra safety - if lookup still returns NA, don't draw
        if (is.na(fill_color) || is.na(fill_alpha)) {
          return(grid::nullGrob())
        }

        # Draw an ellipse matching old legend style (6pt x 4pt radius)
        theta <- seq(0, 2*pi, length.out = 30)
        grid::polygonGrob(
          x = grid::unit(0.5, "npc") + grid::unit(6 * cos(theta), "pt"),
          y = grid::unit(0.5, "npc") + grid::unit(4 * sin(theta), "pt"),
          gp = grid::gpar(fill = fill_color, col = NA, alpha = fill_alpha)
        )
      }
    })

    # Add highlight legend using SHAPE aesthetic only
    # v177: Use unique shape values (1, 2, 3...) for each label
    p <- p +
      geom_point(
        data = highlight_legend_data,
        aes(x = x, y = y, shape = highlight_shape),
        size = 5,
        na.rm = TRUE,
        inherit.aes = FALSE,
        show.legend = TRUE,
        key_glyph = draw_key_highlight_ellipse
      ) +
      scale_shape_manual(
        name = highlight_title,
        values = setNames(highlight_shape_values, highlight_labels),
        guide = guide_legend(order = 97)
      )

    debug_cat(paste0("  v78: Highlight legend added (shape values 1-", n_highlights, ")\n"))
  }

  # ============================================
  # v178: BOOTSTRAP LEGEND using shape aesthetic
  # Custom key_glyph draws triangles with varying sizes
  # v178: Returns nullGrob for shapes not in valid set (prevents bleeding)
  # v178: Match tree's alpha=0.5, fill="grey36", col="grey20"
  # ============================================
  if (show_boot_flag == TRUE && show_bootstrap_legend == TRUE) {
    # Check if bootstrap is in triangles format
    boot_format <- if (!is.null(boot_values) && is.list(boot_values) && !is.null(boot_values$'format')) {
      boot_values$'format'
    } else {
      "triangles"  # default
    }

    debug_cat(paste0("\n  v178: Bootstrap format: '", boot_format, "'\n"))

    if (boot_format == "triangles") {
      debug_cat(paste0("  v78: Creating BOOTSTRAP legend (shape-only)\n"))

      # Bootstrap legend items - match old legend sizes (5, 4, 3 pt)
      bootstrap_labels <- c(">90%", ">80%", ">70%")
      bootstrap_shape_values <- c(1, 2, 3)  # v178: Unique shape values
      bootstrap_pt_sizes <- c(5, 4, 3)  # Match old legend: tri_sizes <- c(5, 4, 3)

      bootstrap_legend_data <- data.frame(
        x = rep(NA_real_, 3),
        y = rep(NA_real_, 3),
        bootstrap_shape = factor(bootstrap_labels, levels = bootstrap_labels),
        stringsAsFactors = FALSE
      )

      # v177: Create mapping from shape value to size (in pt)
      shape_to_size <- setNames(bootstrap_pt_sizes, as.character(bootstrap_shape_values))

      debug_cat(paste0("  v78: Bootstrap shape values: 1=5pt, 2=4pt, 3=3pt\n"))

      # v178: Custom key_glyph for TRIANGLE - draws triangle with varying size
      # Uses shape value (1, 2, 3) to look up the correct size
      # v178: Returns nullGrob if shape value not in our set (prevents legend bleeding)
      draw_key_bootstrap_triangle <- local({
        shape_to_size_local <- shape_to_size
        valid_shapes <- names(shape_to_size)  # v178: Track valid shape values ("1", "2", "3")
        function(data, params, size) {
          # v178: Get the shape value and use it to look up size
          shape_val <- data$shape
          debug_cat(paste0("    v178: draw_key_bootstrap_triangle called, shape=", shape_val, "\n"))

          if (is.null(shape_val) || length(shape_val) == 0) {
            debug_cat(paste0("    v178: shape is NULL/empty, returning nullGrob\n"))
            return(grid::nullGrob())
          }

          shape_key <- as.character(shape_val)

          # v178: CRITICAL - Return nullGrob if this shape is not one of ours
          # This prevents triangles from appearing on other legends (classification, heatmap, highlight, etc.)
          if (!(shape_key %in% valid_shapes)) {
            debug_cat(paste0("    v178: shape ", shape_key, " not in valid set (", paste(valid_shapes, collapse=","), "), returning nullGrob\n"))
            return(grid::nullGrob())
          }

          sz <- shape_to_size_local[shape_key]

          debug_cat(paste0("    v178: Triangle size='", sz, "pt'\n"))

          # v178: Extra safety - if lookup returns NA, don't draw
          if (is.na(sz)) {
            return(grid::nullGrob())
          }

          # Draw triangle matching old legend style
          # Old code: x = c(-sz, sz, 0), y = c(-sz*0.6, -sz*0.6, sz*0.8)
          grid::polygonGrob(
            x = grid::unit(0.5, "npc") + grid::unit(c(-sz, sz, 0), "pt"),
            y = grid::unit(0.5, "npc") + grid::unit(c(-sz*0.6, -sz*0.6, sz*0.8), "pt"),
            gp = grid::gpar(fill = "grey36", col = "grey20", alpha = 0.5)
          )
        }
      })

      # Add bootstrap legend using SHAPE aesthetic only
      # Use new_scale("shape") to create a separate shape scale for bootstrap
      p <- p + ggnewscale::new_scale("shape")

      p <- p +
        geom_point(
          data = bootstrap_legend_data,
          aes(x = x, y = y, shape = bootstrap_shape),
          size = 5,
          color = "grey36",
          na.rm = TRUE,
          inherit.aes = FALSE,
          show.legend = TRUE,
          key_glyph = draw_key_bootstrap_triangle
        ) +
        scale_shape_manual(
          name = "Bootstrap",
          values = setNames(bootstrap_shape_values, bootstrap_labels),  # v177: unique values 1, 2, 3
          guide = guide_legend(order = 98)
        )

      debug_cat(paste0("  v78: Bootstrap legend added (shape values 1-3)\n"))
    } else {
      # v180: For non-triangle formats (percentage, raw, color-coded), show no legend
      debug_cat(paste0("  v80: Bootstrap format '", boot_format, "' - no legend (only triangles have legend)\n"))
    }
  }

  debug_cat(paste0("\n  v180: Highlight and Bootstrap legends complete\n"))
  debug_cat(paste0("  v78: Existing scales preserved: fill (heatmaps), size (P value), colour (classification)\n"))
  debug_cat(paste0("=================================================\n"))

  return(p)
}


# Function to calculate the minimum column x of frame of previous heat
func.calc.min_col_x_of_frame_of_prev_heat <- function(pr440_short_tips_TRY_heat, colnames_sub_df_heat_j) {
  buildd <- ggplot_build(pr440_short_tips_TRY_heat)
  
  build_length <- length(buildd$data)
  heat_build_idx <- build_length
  
  for (index_build in 1:build_length) {
    if ('label' %in% colnames(buildd$data[[index_build]])) {
      build_frame <- buildd$data[[index_build]]
      col_label <- build_frame$label
      if (col_label[1] == colnames_sub_df_heat_j[1]) {
        heat_build_idx <- index_build
      }
    }
  }
  
  frame_of_prev_heat <- buildd$data[[heat_build_idx - 1]]
  col_x_of_frame_of_prev_heat <- frame_of_prev_heat$x
  min_col_x_of_frame_of_prev_heat <- min(frame_of_prev_heat$x)
  min_col_x_of_frame_of_prev_heat <- min_col_x_of_frame_of_prev_heat * (-1)
  
  return(min_col_x_of_frame_of_prev_heat)
}


# Function to make highlight parameters
func.make.highlight.params.NEW <- function(yaml_file, title.id, ids_list, tree440, readfile440, hi_def, debug_mode,
                                           id_tip_trim_start, id_tip_trim_end, id_tip_trim_flag) {
  
  # === COMPREHENSIVE DEBUG OUTPUT ===
  # v53: cat(file=stderr(), "\nÃ°Å¸â€ÂÃ°Å¸â€ÂÃ°Å¸â€Â ENTERING func.make.highlight.params.NEW Ã°Å¸â€ÂÃ°Å¸â€ÂÃ°Å¸â€Â\n")
  # v53: debug_cat("==========================================\n")
  # v53: cat(file=stderr(), "hi_def structure:\n")
  # v54: str(hi_def)
  # v53: cat(file=stderr(), "\n")
  # v53: cat(file=stderr(), "readfile440 dimensions:", nrow(readfile440), "x", ncol(readfile440), "\n")
  # v53: cat(file=stderr(), "readfile440 column names:", paste(colnames(readfile440), collapse=", "), "\n")
  # v53: cat(file=stderr(), "title.id:", title.id, "\n")
  # v53: cat(file=stderr(), "id_tip_trim_flag:", id_tip_trim_flag, "\n")
  # v53: cat(file=stderr(), "id_tip_trim_start:", id_tip_trim_start, "\n")
  # v53: cat(file=stderr(), "id_tip_trim_end:", id_tip_trim_end, "\n")
  
  # Check first few rows of readfile440
  # v53: cat(file=stderr(), "\nFirst 3 rows of readfile440:\n")
  # v53: print(head(readfile440, 3))
  # v53: cat(file=stderr(), "\n")
  # v53: debug_cat("==========================================\n\n")
  # === END DEBUG ===
  
  len_hi <- length(hi_def$according)
  indexes_hi <- 1:len_hi
  
  # v53: print("in func.make.highlight.params.NEW")
  # v53: print("indexes_hi is")
  # v53: print(indexes_hi)
  # v53: print("hi_def is")
  # v53: print(hi_def)
  # v53: print("id_tip_trim_start is")
  # v53: print(id_tip_trim_start)
  # v53: print("id_tip_trim_end is")
  # v53: print(id_tip_trim_end)
  # v53: print("id_tip_trim_flag")
  
  how_many_hi <- length(indexes_hi)
  high_label_list <- c()
  
  for (in_hi in indexes_hi) {
    high_label_list[[in_hi]] <- hi_def$according[[in_hi]][[as.character(in_hi)]]$display_name
  }
  
  high_color_list <- c()
  for (in_hi in indexes_hi) {
    high_color_list[[in_hi]] <- hi_def$according[[in_hi]][[as.character(in_hi)]]$color
  }

  # v139: Extract transparency (alpha) for each highlight
  # v146: Added debug output to trace transparency values
  debug_cat(paste0("\n=== v146: EXTRACTING HIGHLIGHT TRANSPARENCY ===\n"))
  high_alpha_list <- c()
  for (in_hi in indexes_hi) {
    alpha_val <- hi_def$according[[in_hi]][[as.character(in_hi)]]$transparency
    debug_cat(paste0("  Highlight ", in_hi, ":\n"))
    debug_cat(paste0("    hi_def$according[[", in_hi, "]][[\"", in_hi, "\"]]$transparency = ",
                              if(is.null(alpha_val)) "NULL" else alpha_val, "\n"))
    # Default to 0.5 if transparency not specified
    high_alpha_list[[in_hi]] <- if (!is.null(alpha_val)) alpha_val else 0.5
    debug_cat(paste0("    high_alpha_list[[", in_hi, "]] = ", high_alpha_list[[in_hi]], "\n"))
  }
  debug_cat(paste0("  Final high_alpha_list: ", paste(high_alpha_list, collapse=", "), "\n"))
  debug_cat(paste0("==============================================\n"))

  high_title_list <- c()
  for (in_hi in indexes_hi) {
    high_title_list[[in_hi]] <- hi_def$according[[in_hi]][[as.character(in_hi)]]$display_title
  }
  
  offset_hi <- 0
  if ('offset' %in% names(hi_def)) {
    offset_hi <- hi_def$offset
  }
  
  vertical_offset_hi <- 0
  if ('vertical_offset' %in% names(hi_def)) {
    vertical_offset_hi <- hi_def$vertical_offset
  }
  
  # v50: Default to 1.0 for multiplicative scaling (1.0 = no change)
  adjust_height_ecliplse <- 1.0
  if ("adjust_height" %in% names(hi_def)) {
    adjust_height_ecliplse <- hi_def$adjust_height
  }
  
  # v50: Default to 1.0 for multiplicative scaling
  adjust_width_eclipse <- 1.0
  if ("adjust_width" %in% names(hi_def)) {
    adjust_width_eclipse <- hi_def$adjust_width
  }
  
  # v56d: Suppress harmless fortify warnings
  tab <- suppressWarnings(ggtree(tree440))$data
  
  lists_list_hi <- c()
  for (in_hi in indexes_hi) {
    # v53: print("in_hi is")
    # v53: print(in_hi)
    temp_list <- c()
    node_num_list <- tab$node
    # v53: print("node_num_list is")
    # v53: print(node_num_list)
    title_i <- paste0('title', as.character(1))
    value_i <- paste0('value', as.character(1))
    # v53: print("title_i is")
    # v53: print(title_i)
    # v53: print("value_i is")
    # v53: print(value_i)
    
    title_i_name <- hi_def$according[[in_hi]][[as.character(in_hi)]][[title_i]]
    value_i_name <- hi_def$according[[in_hi]][[as.character(in_hi)]][[value_i]]
    
    # === MORE DEBUG ===
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â EXTRACTED FROM YAML:\n")
    # v53: cat(file=stderr(), "title_i_name (column):", title_i_name, "\n")
    # v53: cat(file=stderr(), "value_i_name (value to find):", value_i_name, "\n")
    
    # Check if column exists
    if (title_i_name %in% colnames(readfile440)) {
      # v53: cat(file=stderr(), "Ã¢Å“â€œ Column '", title_i_name, "' EXISTS in readfile440\n", sep="")
      unique_vals <- unique(readfile440[[title_i_name]])
      # v53: cat(file=stderr(), "  Unique values in column:", paste(unique_vals, collapse=", "), "\n")
      # v53: cat(file=stderr(), "  value_i_name ('", value_i_name, "') in unique values:", 
      #     value_i_name %in% unique_vals, "\n", sep="")
    } else {
      # v53: cat(file=stderr(), "Ã¢Å“â€” Column '", title_i_name, "' DOES NOT EXIST in readfile440\n", sep="")
    }
    # v53: cat(file=stderr(), "\n")
    # === END DEBUG ===
    
    # v53: print("title_i_name is")
    # v53: print(title_i_name)
    # v53: print("value_i_name is")
    # v53: print(value_i_name)
    
    if (!title_i_name %in% colnames(readfile440)) {
      # v53: print("wrong highlight parameter")
      stop(paste(title_i_name, "is not a title in classification file"))
    }
    
    match_count <- 0  # Track how many matches we find
    
    for (node_num in node_num_list) {
      row <- subset(tab, node == node_num)
      is_high <- FALSE
      
      if (row$isTip == TRUE) {
        if (id_tip_trim_flag == TRUE) {
          id <- substr(row$label, id_tip_trim_start, id_tip_trim_end)
        } else {
          id <- row$label
        }
        
        # === DEBUG FOR EACH NODE ===
        if (node_num <= 5) {  # Only print first 5 to avoid spam
          # v53: cat(file=stderr(), "  Node", node_num, "- label:", row$label, "- extracted id:", id, "\n")
        }
        # === END DEBUG ===
        
        row_file <- readfile440[readfile440[[title.id]] == id, ]
        
        # === DEBUG ===
        if (node_num <= 5) {
          # v53: cat(file=stderr(), "    Matching in CSV with title.id='", title.id, "' == '", id, "'\n", sep="")
          # v53: cat(file=stderr(), "    Found", nrow(row_file), "matching row(s)\n")
          if (nrow(row_file) > 0) {
            # v53: cat(file=stderr(), "    Value in column:", row_file[[title_i_name]], "\n")
          }
        }
        # === END DEBUG ===
        
        if (nrow(row_file) == 1) {
          val <- row_file[[title_i_name]]
          
          if (is.na(val)) {
            # Skip
          } else {
            if (val == value_i_name) {
              is_high <- TRUE
              match_count <- match_count + 1
              if (match_count <= 5) {  # Print first 5 matches
                # v53: cat(file=stderr(), "    Ã¢Å“â€œ MATCH FOUND! Node", node_num, "has value:", val, "\n")
              }
            }
          }
        }
      }
      
      temp_list <- c(temp_list, is_high)
    }
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â SUMMARY: Found", match_count, "matching nodes out of", length(node_num_list), "total nodes\n\n")
    
    lists_list_hi[[in_hi]] <- temp_list
  }
  
  highlight.params.NEW <- c()
  highlight.params.NEW$how_many_hi <- how_many_hi
  highlight.params.NEW$high_label_list <- high_label_list
  highlight.params.NEW$high_color_list <- high_color_list
  highlight.params.NEW$high_alpha_list <- high_alpha_list  # v139: Add alpha list
  highlight.params.NEW$high_title_list <- high_title_list
  highlight.params.NEW$lists_list_hi <- lists_list_hi
  highlight.params.NEW$offset_hi <- offset_hi
  highlight.params.NEW$vertical_offset_hi <- vertical_offset_hi
  highlight.params.NEW$adjust_height_ecliplse <- adjust_height_ecliplse
  highlight.params.NEW$adjust_width_eclipse <- adjust_width_eclipse

  return(highlight.params.NEW)
}




# Function to make rotation calculations
func.calc.rotation.definitions <- function(rot1, rot2, yaml_file, title.id, ids_list, tree440,
                                           readfile440, debug_mode, id_tip_trim_flag, id_tip_prefix) {
  rotate_flag_for_title <- "no"
  rotate_str <- ""
  rotation_params1 <- c()
  rotation_params2 <- c()
  rotate1_types_list_dx.rx <- NA
  list_weight_frac1 <- NA
  
  if (func.check.bin.val.from.conf(rot1) == FALSE) {
    rotate_flag <- "no"
    rotate_flag_for_title <- "no"
    rotate_str <- "no_rotate"
  } else {
    rotate_flag <- 'RX_first'
    rotate_flag_for_title <- "yes"
    
    rotation_params1 <- func.make.rot.params(1, yaml_file, title.id, ids_list, tree440,
                                             readfile440, id_tip_trim_flag, id_tip_prefix) 
    
    rotate1_types_list_dx.rx <- rotation_params1['types_list_dx.rx']
    
    if (debug_mode == TRUE) {
      # v53: print("##DEBUG##")
      # v53: print("rotation_params1 is")
      # v53: print(rotation_params1)
    }
    
    st <- c()
    for (ind_s in rotate1_types_list_dx.rx) {
      st <- paste0(st, ind_s, collapse = "_")
    }
    
    if (nchar(st) > 30) {
      st <- substr(st, 1, 30)
      st <- paste0(st, "_ETC")
    }
    
    rotate_str <- paste0('rotate_by_', st)
    
    if (func.check.bin.val.from.conf(rot2) == FALSE) {
      rotation_params2 <- rotation_params1
    } else {
      rotation_params2 <- func.make.rot.params(2, yaml_file, title.id, ids_list, tree440,
                                               readfile440, id_tip_trim_flag, id_tip_prefix)  
      list_weight_frac1 <- rotation_params2['list_weight_dx.rx']
      
      if (debug_mode == TRUE) {
        # v53: print("##DEBUG##")
        # v53: print("rotation_params2 is")
        # v53: print(rotation_params2)
        # v53: print("list_weight_frac1=")
        # v53: print(list_weight_frac1)
      }
    }
  }
  
  l <- c()
  l$rotate_flag_for_title <- rotate_flag_for_title
  l$rotate_str <- rotate_str
  l$rotation_params1 <- rotation_params1
  l$rotation_params2 <- rotation_params2
  l$rotate_flag <- rotate_flag
  l$rotate1_types_list_dx.rx <- rotate1_types_list_dx.rx
  l$list_weight_frac1 <- list_weight_frac1
  
  return(l)
}



# Function to rotate tree based on weights
func.rotate.tree.based.on.weights <- function(tree_TRY, list_weight_dx.rx, list_weight_frac,
                                              TREE_OTU_dx.rx, TREE_OTU_frac, tree_size) {
  # Rotating tree
  tree_return <- tree_TRY
  list_weights_for_nodes_dx.rx <- func.create.weight_list(tree_return, list_weight_dx.rx, TREE_OTU_dx.rx, tree_size)
  list_weights_for_nodes_frac <- func.create.weight_list(tree_return, list_weight_frac, TREE_OTU_frac, tree_size)
  
  if (NA %in% list_weights_for_nodes_dx.rx) {
    # v53: print("Error: rotation classes do not cover all tips")
    # v53: print("weight list 1 is")
    # v53: print(list_weights_for_nodes_dx.rx)
    # v53: print("cells without weights are")
    u <- which(is.na(list_weights_for_nodes_dx.rx))
    # v53: print("which are")
    # v53: print(tree_TRY$data$label[u])
    # v53: print(tree_TRY$data$group[u])
    # v53: print("NA weights are rewritten to 0, please verify yml rotation definitions")
    
    idxx <- 1
    for (ww in list_weights_for_nodes_dx.rx) {
      if (is.na(ww)) {
        list_weights_for_nodes_dx.rx[idxx] <- 0
      }
      idxx <- idxx + 1
    }
  }
  
  tree_TRY$data$weight1 <- list_weights_for_nodes_dx.rx
  tree_TRY$data$weight2 <- list_weights_for_nodes_frac
  
  for (nod in tree_TRY$data$node) {
    isT <- tree_TRY$data[[nod, "isTip"]]
    if (isT == FALSE) {
      children <- which(tree_TRY$data$parent == nod)
      children_weights <- rep(0, 2)
      children_weights_SECOND <- rep(0, 2)
      children_side_orig <- rep("l", 2)
      children_side_new <- rep("l", 2)
      
      if (length(children) == 2) {
        children_weights[1] <- list_weights_for_nodes_dx.rx[children[1]]
        children_weights[2] <- list_weights_for_nodes_dx.rx[children[2]]
        children_weights_SECOND[1] <- list_weights_for_nodes_frac[children[1]]
        children_weights_SECOND[2] <- list_weights_for_nodes_frac[children[2]]
        
        if (tree_TRY$data$y[children[1]] > tree_TRY$data$y[children[2]]) {
          children_side_orig[1] <- "l"
          children_side_orig[2] <- "r"
        } else {
          children_side_orig[1] <- "r"
          children_side_orig[2] <- "l"
        }
        
        if (children_weights[1] < children_weights[2]) {
          children_side_new[1] <- "r"
          children_side_new[2] <- "l"    
        } else if (children_weights[1] > children_weights[2]) {
          children_side_new[1] <- "l"
          children_side_new[2] <- "r"   
        } else {
          if (is.na(children_weights_SECOND[1]) | is.na(children_weights_SECOND[2])) {
            children_side_new[1] <- children_side_orig[1]
            children_side_new[2] <- children_side_orig[2]
          } else {
            if (children_weights_SECOND[1] < children_weights_SECOND[2]) {
              children_side_new[1] <- "r"
              children_side_new[2] <- "l"    
            } else if (children_weights_SECOND[1] > children_weights_SECOND[2]) {
              children_side_new[1] <- "l"
              children_side_new[2] <- "r"  
            } else {
              children_side_new[1] <- children_side_orig[1]
              children_side_new[2] <- children_side_orig[2]
            }
          }
        }
        
        if (children_side_new[1] != children_side_orig[1]) {
          tree_return <- flip(tree_return, children[1], children[2])
        }
      } else {
        # Multifurcating node (more than 2 children)
        cat(file=stderr(), paste0("\n[DEBUG-WEIGHT-ROT] Multifurcating node ", nod, " with ", length(children), " children\n"))
        cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Children indices: ", paste(children, collapse=", "), "\n"))

        children_weights <- func.make.children.weight.list(children, nod, list_weights_for_nodes_dx.rx)
        children_weights_SECOND <- func.make.children.weight.list(children, nod, list_weights_for_nodes_frac)

        cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Children weights: ", paste(children_weights, collapse=", "), "\n"))

        # Sort DESCENDING to match binary node behavior (higher weight = left/top)
        children_weights_ordered <- sort(children_weights, decreasing = TRUE)
        cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Sorted weights (DESCENDING - high to low): ", paste(children_weights_ordered, collapse=", "), "\n"))

        # Build destination mapping: dest[i] = which original position should go to position i
        dest <- c()
        for (i in 1:length(children)) {
          wh <- which(children_weights == children_weights_ordered[i])

          if (length(wh) > 1) {
            w <- wh[which(!wh %in% dest)[1]]
          } else {
            w <- wh
          }
          dest <- c(dest, w)
        }

        cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Destination mapping (dest): ", paste(dest, collapse=", "), "\n"))
        cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] This means: position 1 should have child from orig pos ", dest[1], ", pos 2 from orig pos ", dest[2], ", etc.\n"))

        # Track which positions have been processed to avoid undoing swaps
        processed <- rep(FALSE, length(children))

        for (i in 1:length(children)) {
          # Skip if this position was already processed (as part of a previous swap)
          if (processed[i]) {
            cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Loop i=", i, ": Already processed, skipping\n"))
            next
          }

          target_pos <- dest[i]

          if (i == target_pos) {
            # Child is already in correct position
            cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Loop i=", i, ": Child ", children[i], " already in correct position\n"))
            processed[i] <- TRUE
          } else {
            # Need to swap position i with position target_pos
            cat(file=stderr(), paste0("[DEBUG-WEIGHT-ROT] Loop i=", i, ": Swapping children[", i, "]=", children[i], " <-> children[", target_pos, "]=", children[target_pos], "\n"))
            tree_return <- flip(tree_return, children[i], children[target_pos])

            # Mark both positions as processed
            processed[i] <- TRUE
            processed[target_pos] <- TRUE

            # Update dest to reflect the swap (in case there are chained swaps)
            # Find where target_pos was supposed to go and update it
            for (j in (i+1):length(children)) {
              if (j <= length(dest) && dest[j] == i) {
                dest[j] <- target_pos
              }
            }
          }
        }
      }
    } 
  }
  
  return(tree_return)
}


# Function to make children weight list
func.make.children.weight.list <- function(children, nod, tips_weight_list) {
  children_weights <- rep(0, length(children))
  for (idx in 1:length(children)) {
    if (nod == children[idx]) {
      children_weights[idx] <- 0
    } else {
      children_weights[idx] <- tips_weight_list[children[idx]]
    }
  }
  return(children_weights)
}


# Function to rotate specific nodes
func.rotate.specific.nodes <- function(tree_TRY1, list_nodes_to_rotate) {
  cat(file=stderr(), paste0("[DEBUG-ROTATION] func.rotate.specific.nodes called\n"))
  cat(file=stderr(), paste0("[DEBUG-ROTATION] list_nodes_to_rotate = ", paste(list_nodes_to_rotate, collapse=", "), "\n"))
  cat(file=stderr(), paste0("[DEBUG-ROTATION] is.null = ", is.null(list_nodes_to_rotate), ", length = ", length(list_nodes_to_rotate), "\n"))

  # Check if list_nodes_to_rotate is valid (not NA and has length > 0)
  if (!is.null(list_nodes_to_rotate) && length(list_nodes_to_rotate) > 0 && !all(is.na(list_nodes_to_rotate))) {
    cat(file=stderr(), paste0("[DEBUG-ROTATION] Entering rotation loop for nodes: ", paste(list_nodes_to_rotate, collapse=", "), "\n"))
    for (nod in list_nodes_to_rotate) {
      children <- which(tree_TRY1$data$parent == nod & tree_TRY1$data$node != nod)
      cat(file=stderr(), paste0("[DEBUG-ROTATION] Node ", nod, " has ", length(children), " children: ", paste(children, collapse=", "), "\n"))

      if (length(children) < 2) {
        cat(file=stderr(), paste0("[DEBUG-ROTATION] Skipping node ", nod, " - less than 2 children\n"))
        # Cannot flip with less than 2 children - skip this node
        next
      } else if (length(children) == 2) {
        cat(file=stderr(), paste0("[DEBUG-ROTATION] Flipping node ", nod, " children: ", children[1], " <-> ", children[2], "\n"))
        tree_TRY1 <- flip(tree_TRY1, children[1], children[2])
      } else {
        # Multifurcating node - flip first and last children
        cat(file=stderr(), paste0("[DEBUG-ROTATION] Multifurcating node ", nod, " - flipping: ", children[1], " <-> ", children[length(children)], "\n"))
        tree_TRY1 <- flip(tree_TRY1, children[1], children[length(children)])
      }
    }
  } else {
    cat(file=stderr(), paste0("[DEBUG-ROTATION] No valid nodes to rotate\n"))
  }

  return(tree_TRY1)
}




func.print.lineage.tree <- function(conf_yaml_path,
                                    width=170,height= 60,laderize_flag=FALSE,simulate.p.value=TRUE,
                                    man_adjust_elipse=0.005, man_multiply_elipse= 1.3,
                                    man_adj_second_legend= 0,man_space_second_legend = -0.02,
                                    add_date_to_text_flag=TRUE,
                                    man_adjust_image_of_second_legend=0, man_adj_heat_loc=0,man_boot_x_offset=0,
                                    man_adj_heat_loc2=0, man_adj_heat_loc3=0,
                                    id_tip_trim_flag= TRUE, 
                                    id_tip_trim_start= 3,
                                    id_tip_trim_end=7,
                                    id_tip_prefix='',
                                    debug_mode=FALSE, 
                                    debug_print_data_tree= TRUE,
                                    man_multiply_second_legend= 1,
                                    man_multiply_second_legend_text=1.7,
                                    man_multiply_first_legend_text= 1,
                                    man_multiply_first_legend_title_size=1,
                                    man_space_second_legend_multiplier=1,
                                    man_offset_for_highlight_legend_x=0,
                                    a_4_output=FALSE,
                                    man_adjust_elipse_a=0,
                                    man_adjust_elipse_b=0,
                                    tree_path_list_bash= NA,
                                    csv_path_bash= NA,
                                    out_path_bash= NA,
                                    out_file_name_bash=NA,
                                    out_file_type= NA,
                                    heat_maps_titles_angles_vector= NA,
                                    man_offset_second_legend=0,
                                    units_out="cm",
                                    compare_two_trees= NA,
                                    trees_to_compare= NA,
                                    classification_in_compare=NA,
                                    flag_print_tree_data= FALSE,
                                    list_nodes_to_rotate= NA,
                                    flag_display_nod_number_on_tree= FALSE,
                                    node_number_font_size= 3.5,
                                    highlight_manual_nodes= FALSE,
                                    manual_nodes_to_highlight= NA,
                                    flag_calc_scores_for_tree= FALSE,
                                    flag_make_newick_file= FALSE,
                                    path_out_newick="",
                                    width_heatmap=NA,
                                    flag_colnames=NA,
                                    viridis_option_list=viridis_option_list,
                                    rotate_compared_tree= FALSE, 
                                    PERM_BOUND=7,
                                    how_many_rotation_vec=c(1,2),
                                    flag_csv_read_func="def",
                                    rowname_param="",
                                    heat_legend_replace=NA,
                                    tip_name_display_flag=TRUE,
                                    bootstrap_label_size= 1.5,  # v129: Reduced from 3.5 for smaller default legend
                                    heatmap_tree_distance= 0.02,
                                    heatmap_global_gap = 0.05,  # v125: Gap between multiple heatmaps
                                    legend_settings = NULL,  # v135: Legend settings for highlight/bootstrap legends
                                    rdata_cnv_matrix = NULL,  # S1.62dev: CNV matrix from RData file
                                    cached_p_list_of_pairs = NULL,  # S2.0-PERF: Cached p-values (Option 3A)
                                    cached_p_list_hash = NULL,      # S2.0-PERF: Hash for cache validation
                                    p_list_cache = list(),          # S2.7-PERF: Multi-entry cache (hash -> p_list_of_pairs)
                                    heatmap_cache = list()) {       # S2.9-PERF: Heatmap cache (index -> tile_df + hash)

  # === DEBUG CHECKPOINT 2: FUNCTION ENTRY ===
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 2: func.print.lineage.tree ENTRY\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â node_number_font_size received:", node_number_font_size, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â flag_display_nod_number_on_tree:", flag_display_nod_number_on_tree, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â highlight_manual_nodes received:", highlight_manual_nodes, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â manual_nodes_to_highlight received:", paste(manual_nodes_to_highlight, collapse=", "), "\n")
  # v53: debug_cat("================================================\n\n")

  # === PROFILING: func.print.lineage.tree ===
  .prof_func_start <- Sys.time()
  .prof_section_start <- Sys.time()

  yaml_file<- func.read_yaml(conf_yaml_path)
  cat(file=stderr(), sprintf("[PROF-TREE] YAML read: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))

  #get configuration from yaml file
  
  # v53: print(paste0("Configuration file: ",conf_yaml_path))
  
  
  out_list<-c()
  out_index <-1
  
  
  #units_out <- 'cm'
  if (a_4_output== TRUE) {
    # v53: print("output in A4 format")
    width <- 297
    height <- 210
    units_out <-  "mm"
    man_multiply_second_legend_text <- 0.34
    man_multiply_first_legend_title_size <- 0.25
    
  }
  
  
  individual<- yaml_file[['Individual general definitions']]$Individual
  
  if (is.na(individual)) {
    stop("Missing individual name")
  }
  
  if (is.na(csv_path_bash)) {
    # v53: print("Use csv_path from yaml")
    csv_path<- yaml_file[['Individual general definitions']]$'mapping csv file'
  } else {
    csv_path<- csv_path_bash
    # v53: print("Use csv_path_bash")
  }
  
  
  
  if(is.na(csv_path)) {
    stop("Missing csv file")
  }
  
  #get csv file for mapping subgroups
  # v53: print(paste0("Get mapping csv from: ",csv_path))
  .prof_section_start <- Sys.time()
  if (flag_csv_read_func=="fread"){
    fread_rownames(csv_path, row.var = rowname_param)
  } else {
    readfile <- read.csv(csv_path)
  }
  cat(file=stderr(), sprintf("[PROF-TREE] CSV read (in func): %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  #print("readfile is")
  #print(readfile)
  #print("####")
  
  
  
  
  
  tree_path_list<- yaml_file[['Individual general definitions']]$'tree path'
  
  title.id <- yaml_file[['Mapping exl renaming titles']]$'ID column'
  
  
  if (debug_mode== TRUE) {
    # v53: print("title.id")
    # v53: print(title.id)
    
  }
  
  
  if (is.na(tree_path_list_bash)) {
    # v53: print("Use tree_path from yaml")
  } else {
    # v53: print("Use tree_path_list_bash")
    tree_path_list<- tree_path_list_bash
  }
  
  
  trees_number <- length(tree_path_list)
  if (trees_number==0) {
    stop("No tree path was given")
  }
  
  flag_short_tips_yaml <- yaml_file[['visual definitions']]$'trim tips'$'display'
  if (!is.null(flag_short_tips_yaml)) {
    # Convert "yes"/"no" to TRUE/FALSE
    if (is.character(flag_short_tips_yaml)) {
      flag_short_tips <- (tolower(flag_short_tips_yaml) == "yes")
    } else {
      flag_short_tips <- as.logical(flag_short_tips_yaml)
    }
    # v53: print(paste0("flag_short_tips from YAML: ", flag_short_tips))
  } else {
    flag_short_tips <- TRUE  # Default
    # v53: print("flag_short_tips not in YAML, using default: TRUE")
  }
  
  tips_length <- as.numeric(yaml_file[['visual definitions']]$'trim tips'$'length')
  if (is.na(tips_length) || is.null(tips_length)) {
    tips_length <- 0.05  # Default
    # v53: print("tips_length not in YAML, using default: 0.05")
  } else {
    # v53: print(paste0("tips_length from YAML: ", tips_length))
  }
  
  path_base <- yaml_file[['Individual general definitions']]$out_file$base_path
  
  if (is.na(out_path_bash)) {
    # v53: print("Use path_base from yaml for output path")
    
  } else {
    path_base<- out_path_bash
    # v53: print("Use out_path_bash for output path")
  }
  
  
  edge_width_multiplier<-yaml_file[['visual definitions']]$'edge_width_multiplier'$'size'
  if (is.na(edge_width_multiplier)) {
    # v53: print("Missing edge_width_multiplier in yaml file")
    # v53: print("setting to default value")
    edge_width_multiplier<-1
  }
  size_tip_text<-yaml_file[['visual definitions']]$'font_size'$'tips'
  if (is.na(size_tip_text)) {
    # v53: print("Missing size_tip_text in yaml file")
    # v53: print("setting to default value")
    size_tip_text <- 10
  }
  
  size_font_legend_title <- yaml_file[['visual definitions']]$'font_size'$'legend_title'
  if(is.na(size_font_legend_title)) {
    # v53: print("Missig font size for legend title")
    # v53: print("setting to default value")
    size_font_legend_title <- 50
  }
  size_font_legend_text <- yaml_file[['visual definitions']]$'font_size'$'legend_text' 
  if (is.na(size_font_legend_text)) {
    # v53: print("Missing font size for legend text")
    # v53: print("setting to default value")
    size_font_legend_text <- 35
  }
  size_font_legend_box <- yaml_file[['visual definitions']]$'font_size'$'legend_box' 
  if (is.na(size_font_legend_box)) {
    # v53: print("Missing font size for lengend box")
    # v53: print("setting to ddefault value")
    size_font_legend_box <- 30
  }
  
  size_font_heat_map_text <- yaml_file[['visual definitions']]$'font_size'$'heat_map_title' 
  if (is.na(size_font_legend_box)) {
    size_font_heat_map_text <- 25
  }
  
  size_font_heat_map_legend <- yaml_file[['visual definitions']]$'font_size'$'heat_map_legend' 
  if (length(size_font_heat_map_legend)==0){
    size_font_heat_map_legend <- size_tip_text *1.6
  }
  if (is.na(size_font_heat_map_legend)) {
    size_font_heat_map_legend <- size_tip_text *1.6
  }
  
  # Read laderize_flag from YAML (if it exists)
  laderize_flag_yaml <- yaml_file[['visual definitions']]$'laderize_flag'
  if (!is.null(laderize_flag_yaml)) {
    # Convert "yes"/"no" to TRUE/FALSE
    if (is.character(laderize_flag_yaml)) {
      laderize_flag <- (tolower(laderize_flag_yaml) == "yes")
    } else {
      laderize_flag <- as.logical(laderize_flag_yaml)
    }
    # v53: print(paste0("laderize_flag from YAML: ", laderize_flag))
  }
  
  # Read flag_display_nod_number_on_tree from YAML (if it exists)
  display_node_flag_yaml <- yaml_file[['visual definitions']]$'flag_display_nod_number_on_tree'
  if (!is.null(display_node_flag_yaml)) {
    # Convert "yes"/"no" to TRUE/FALSE
    if (is.character(display_node_flag_yaml)) {
      flag_display_nod_number_on_tree <- (tolower(display_node_flag_yaml) == "yes")
    } else {
      flag_display_nod_number_on_tree <- as.logical(display_node_flag_yaml)
    }
    # v53: print(paste0("flag_display_nod_number_on_tree from YAML: ", flag_display_nod_number_on_tree))
  }
  
  # CRITICAL FIX: Read trimming parameters from YAML (if they exist)
  # This overrides the function defaults with the correctly inferred values
  id_tip_trim_flag_yaml <- yaml_file[['visual definitions']]$'id_tip_trim_flag'
  if (!is.null(id_tip_trim_flag_yaml)) {
    if (is.character(id_tip_trim_flag_yaml)) {
      id_tip_trim_flag <- (tolower(id_tip_trim_flag_yaml) == "yes")
    } else {
      id_tip_trim_flag <- as.logical(id_tip_trim_flag_yaml)
    }
    # v53: print(paste0("id_tip_trim_flag from YAML: ", id_tip_trim_flag))
  }
  
  id_tip_trim_start_yaml <- yaml_file[['visual definitions']]$'id_tip_trim_start'
  if (!is.null(id_tip_trim_start_yaml)) {
    id_tip_trim_start <- as.numeric(id_tip_trim_start_yaml)
    # v53: print(paste0("id_tip_trim_start from YAML: ", id_tip_trim_start))
  }
  
  id_tip_trim_end_yaml <- yaml_file[['visual definitions']]$'id_tip_trim_end'
  if (!is.null(id_tip_trim_end_yaml)) {
    id_tip_trim_end <- if (is.na(id_tip_trim_end_yaml)) NA else as.numeric(id_tip_trim_end_yaml)
    # v53: print(paste0("id_tip_trim_end from YAML: ", id_tip_trim_end))
  }
  
  id_tip_prefix_yaml <- yaml_file[['visual definitions']]$'id_tip_prefix'
  if (!is.null(id_tip_prefix_yaml)) {
    id_tip_prefix <- as.character(id_tip_prefix_yaml)
    # v53: print(paste0("id_tip_prefix from YAML: ", id_tip_prefix))
  }
  
  out_trees <- c()   
  
  for (tree_index in 1:trees_number) {
    # v53: print(paste0('Tree number ',tree_index))
    tree_path <- tree_path_list[tree_index]

    ###get newick file and create tree
    # v53: print(paste0("Get tree from: ",tree_path))
    .prof_section_start <- Sys.time()
    tree440 <- read.tree(tree_path)
    cat(file=stderr(), sprintf("[PROF-TREE] read.tree(): %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))

    # Apply ladderize if flag is TRUE
    if (laderize_flag == TRUE) {
      # v53: print("Ladderizing tree...")
      .prof_section_start <- Sys.time()
      tree440 <- ape::ladderize(tree440, right = TRUE)
      cat(file=stderr(), sprintf("[PROF-TREE] ladderize(): %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
      # v53: print("Tree ladderized")
    }
    
    # v53: print(tree440)
    
    # v53: print(paste0("Individual name is ",individual))
    
    if (debug_mode==TRUE) {
      # v53: print("individual is")
      # v53: print(individual)
    }
    
    
    
    
    if ("individual" %in% names(readfile)) {
      idv_title <- "individual"
    } else {
      idv_title <- "Individual"
    }
    
    
    
    # v53: print("individuals list in csv:")
    list_individuals= unique(readfile[[idv_title]])
    # v53: print(list_individuals)
    
    str= paste0("choosing only columns with individual= ",individual)
    
    #make subframe for the data of specific individual
    if (idv_title %in% names(readfile)) {   
      readfile440 <- readfile[readfile[[idv_title]]==individual,]
    } else {
      readfile440 <- readfile
    }
    droplist <- c("?","??","???","????","?????")
    #readfile440<- subset(readfile440,!title.id %in% droplist)
    #print("new readfile440 is")
    #print("readfile440 is")
    #print(readfile440)
    #print("_______")
    
    
    
    
    ids_list<-readfile440[[title.id]]
    
    ids_list <- ids_list[! ids_list %in% droplist]
    readfile440<- readfile[readfile[[title.id]]%in% ids_list,]
    #print("Sample.Reads.ID col is")
    #print(readfile440[['Sample.Reads.ID']])
    #        print("_______")
    
    
    
    if (debug_mode== TRUE) {
      
      # v53: print("title.id is")
      # v53: print(title.id)
      # v53: print("ids_list from csv is")
      # v53: print(ids_list)
    }
    
    # v53: print("title.id is")
    #print(title.id)
    #print("ids_list from csv is")
    #print(ids_list)
    
    
    
    readfile440 <- fix.readfile440.with.missing.leaves(readfile440,title.id,tree440,ids_list,"na",
                                                       id_tip_trim_flag,id_tip_trim_start,id_tip_trim_end,id_tip_prefix,debug_mode)
    
    
    #print("AFTER")
    #print(readfile440)
    #print("###")
    
    #                print("Sample.Reads.ID col is")
    #print(readfile440[['Sample.Reads.ID']])
    #        print("_______")
    
    if (debug_mode==TRUE){
      # v53: print("readfile440 is ")
      # v53: print(readfile440)
    }
    
    
    
    
    
    keys_visual <- names(yaml_file[['visual definitions']])
    disp_opt_num <- length(yaml_file[['visual definitions']]$'classification')
    if (is.na(disp_opt_num)){
      stop("Missing classification in visual definitions")
    }
    if (disp_opt_num==0) {
      stop("Missing classification in visual definitions")
    }
    
    #print("B1")
    rot1<-yaml_file[['visual definitions']]$'rotation1'$'display'
    rot2<-yaml_file[['visual definitions']]$'rotation2'$'display'
    
    # Convert "yes"/"no" to logical if needed
    if (!is.null(rot1) && is.character(rot1)) {
      rot1 <- (tolower(rot1) == "yes")
    }
    if (!is.null(rot2) && is.character(rot2)) {
      rot2 <- (tolower(rot2) == "yes")
    }
    
    if (debug_mode== TRUE) {
      # v53: print("##DEBUG##")
      # v53: print("rot1 is")
      # v53: print(rot1)
      # v53: print("rot2 is")
      # v53: print(rot2)
      
    }
    
    # Check if rotation is requested and validate configuration
    if ((rot1 == TRUE || rot2 == TRUE) && is.null(yaml_file[['visual definitions']]$rotation1$'according') && 
        is.null(yaml_file[['visual definitions']]$rotation2$'according')) {
      # v53: print("WARNING: Rotation enabled but no 'according' classification specified. Disabling rotation.")
      rot1 <- FALSE
      rot2 <- FALSE
    }
    
    
    
    
    list_rotate <- func.calc.rotation.definitions(rot1,rot2,yaml_file,title.id,
                                                  ids_list,tree440,readfile440,debug_mode,
                                                  id_tip_trim_flag,id_tip_prefix)
    
    
    
    #print("list_rotate is")
    #print(list_rotate)
    #print("B2")
    
    rotate_flag_for_title <-list_rotate$rotate_flag_for_title
    rotate_str <- list_rotate$rotate_str
    rotation_params1 <- list_rotate$rotation_params1
    rotation_params2 <- list_rotate$rotation_params2
    rotate_flag <- list_rotate$rotate_flag
    rotate1_types_list_dx.rx <- list_rotate$rotate1_types_list_dx.rx
    list_weight_frac1 <- list_rotate$list_weight_frac1
    
    
    #("B3")
    
    #highlight used to be here, moved inside
    
    show_boot <- yaml_file[['visual definitions']]$'Bootstrap'$'display'
    if (func.check.bin.val.from.conf(show_boot)== TRUE) {
      show_boot_flag <- TRUE
    } else {
      show_boot_flag <- FALSE
    }
    boot_values <- yaml_file[['visual definitions']]$'Bootstrap'
    if (! 'format' %in% names(boot_values)) {
      boot_values$'format' <- 'triangles'
      yaml_file[['visual definitions']]$'Bootstrap'$'format' <- 'triangles'
    }
    
    if (! 'classification' %in% keys_visual) {
      stop("Missing classification in yaml")
    }
    
    #print("B4")
    
    for (disp_index in 1:disp_opt_num) {
      #print("disp_index is")
      #print(disp_index)
      disp_index_str <- as.character(disp_index)
      temp <- yaml_file[['visual definitions']]$'classification'[[disp_index]]
      
      att<- names(yaml_file[['visual definitions']]$'classification'[[disp_index]])
      
      
      ttt <- yaml_file[['visual definitions']]$'classification'[[disp_index]]
      ind1<- as.character(disp_index)
      att1 <- names(ttt[[ind1]])
      heat_map_title <- ""
      
      
      
      disp_indx_ch <- as.character(disp_index)
      
      #print("B5")       
      
      FDR_perc <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'FDR_perc'
      
      
      no_name <-  yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'non_cluster_title'
      
      no_name_color <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'non_cluster_color'
      
      
      labels_not_in_legend <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'not_in_legend'
      
      title_replace <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'title'
      
      
      
      dxdf440_for_heat <- NA
      heat_flag <- FALSE
      heat_map_title_list <-c()
      
      
      heat_display_vec=c()
      heat_display_params_list <- c()

      # v56c: DEBUG - show what attributes are in the classification
      debug_cat(paste0("\n=== v57: Classification attributes (att1): ", paste(att1, collapse=", "), " ===\n"))

      if ('heatmap_display' %in% att1) {

        # v56c: DEBUG
        debug_cat("\n=== v57: FOUND heatmap_display in classification ===\n")

        heat_definitions <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$heatmap_display
        heat_list_len <- length(heat_definitions)
        debug_cat(paste0("  Number of heatmaps found: ", heat_list_len, "\n"))
        debug_cat(paste0("  heat_definitions structure: ", class(heat_definitions), "\n"))
        if (heat_list_len > 0) {
          debug_cat(paste0("  First heatmap names: ", paste(names(heat_definitions[[1]]), collapse=", "), "\n"))
        }
        #print("heat_list_len is")
        #print(heat_list_len)
        dxdf440_for_heat <- c()
        #print("B8")
        #print("heat_list_len is")
        #print(heat_list_len)
        
        indx_for_sav <-0
        for (inx in 1:heat_list_len) {
          #print("inx is")
          #print(inx)
          
          
          ind <- as.character(inx)
          heat_map_i_def <- heat_definitions[[inx]][[ind]]
          # S1.62dev: Debug output to trace heatmap definition structure
          debug_cat(paste0("\n=== S1.62dev DEBUG: Heatmap ", inx, " definition ===\n"))
          debug_cat(paste0("  Available fields: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
          debug_cat(paste0("  data_source field exists: ", 'data_source' %in% names(heat_map_i_def), "\n"))
          if ('data_source' %in% names(heat_map_i_def)) {
            debug_cat(paste0("  data_source value: '", heat_map_i_def$data_source, "'\n"))
          }
          #print("B9")
          
          
          att2 <- heat_map_i_def$display
          
          
          if (func.check.bin.val.from.conf(att2)== TRUE) {
            indx_for_sav <- indx_for_sav+1
          } else {
            next
          }
          # print("B10")
          
          param<- c()
          if (func.check.bin.val.from.conf(att2)== TRUE) {
            
            #print("B11")
            
            if ('is_discrete' %in% names(heat_map_i_def)) {
              
              
              param[['is_discrete']] <- func.check.bin.val.from.conf(heat_map_i_def$is_discrete)
              
            } else {
              
              param['is_discrete'] <- FALSE
            } 
            
            
            if (param[['is_discrete']]== TRUE) {
              
              if ('man_define_colors' %in% names(heat_map_i_def)) {
                param['man_define_colors'] <- func.check.bin.val.from.conf(heat_map_i_def$man_define_colors)
              } else {
                param['man_define_colors'] <- FALSE
              }
              
              
              if ('color_scale_option' %in% names(heat_map_i_def)) {
                # v67: FIX - use double brackets to get the value directly, not wrapped in a list
                # Single brackets heat_map_i_def['color_scale_option'] returns list(color_scale_option = "Set1")
                # Double brackets heat_map_i_def[['color_scale_option']] returns just "Set1"
                param[['color_scale_option']] <- heat_map_i_def[['color_scale_option']]
              } else {
                param[['color_scale_option']] <- NULL
              }
              
              
              
              if ('color_scale_range_start' %in% names(heat_map_i_def)) {
                param['man'] <- TRUE
                param['color_scale_range_start'] <- as.numeric(heat_map_i_def$color_scale_range_start)
              } else {
                param['man'] <- FALSE
              }
              if ('color_scale_range_end' %in% names(heat_map_i_def)) {
                param['color_scale_range_end'] <- as.numeric(heat_map_i_def$color_scale_range_end)
              } else {
                param['color_scale_range_end'] <- 300
              }

              # v70: Get NA color (default white)
              if ('na_color' %in% names(heat_map_i_def)) {
                param[['na_color']] <- heat_map_i_def[['na_color']]
              } else {
                param[['na_color']] <- "white"
              }

              # v104: Get per-heatmap distance (default 0.02)
              if ('distance' %in% names(heat_map_i_def)) {
                param[['distance']] <- as.numeric(heat_map_i_def[['distance']])
              } else {
                param[['distance']] <- 0.02
              }

              # v105: Get per-heatmap height (default 0.8)
              if ('height' %in% names(heat_map_i_def)) {
                param[['height']] <- as.numeric(heat_map_i_def[['height']])
              } else {
                param[['height']] <- 0.8
              }

              # v111: Get per-heatmap row_height (default 1.0)
              if ('row_height' %in% names(heat_map_i_def)) {
                param[['row_height']] <- as.numeric(heat_map_i_def[['row_height']])
              } else {
                param[['row_height']] <- 1.0
              }

              # v111: Get grid settings
              if ('show_grid' %in% names(heat_map_i_def)) {
                param[['show_grid']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_grid']])
              } else {
                param[['show_grid']] <- FALSE
              }
              if ('grid_color' %in% names(heat_map_i_def)) {
                param[['grid_color']] <- heat_map_i_def[['grid_color']]
              } else {
                param[['grid_color']] <- "#000000"
              }
              if ('grid_size' %in% names(heat_map_i_def)) {
                param[['grid_size']] <- as.numeric(heat_map_i_def[['grid_size']])
              } else {
                param[['grid_size']] <- 0.5
              }

              # S1.62dev: Get row line settings (horizontal lines only)
              if ('show_row_lines' %in% names(heat_map_i_def)) {
                param[['show_row_lines']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_row_lines']])
              } else {
                param[['show_row_lines']] <- FALSE
              }
              if ('row_line_color' %in% names(heat_map_i_def)) {
                param[['row_line_color']] <- heat_map_i_def[['row_line_color']]
              } else {
                param[['row_line_color']] <- "#000000"
              }
              if ('row_line_size' %in% names(heat_map_i_def)) {
                param[['row_line_size']] <- as.numeric(heat_map_i_def[['row_line_size']])
              } else {
                param[['row_line_size']] <- 0.5
              }
              # S1.62dev: Get column line settings (vertical lines)
              if ('show_col_lines' %in% names(heat_map_i_def)) {
                param[['show_col_lines']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_col_lines']])
              } else {
                param[['show_col_lines']] <- FALSE
              }
              if ('col_line_color' %in% names(heat_map_i_def)) {
                param[['col_line_color']] <- heat_map_i_def[['col_line_color']]
              } else {
                param[['col_line_color']] <- "#000000"
              }
              if ('col_line_size' %in% names(heat_map_i_def)) {
                param[['col_line_size']] <- as.numeric(heat_map_i_def[['col_line_size']])
              } else {
                param[['col_line_size']] <- 0.5
              }

              # S1.62dev: Get vertical text labels settings
              if ('show_vertical_text' %in% names(heat_map_i_def)) {
                param[['show_vertical_text']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_vertical_text']])
              } else {
                param[['show_vertical_text']] <- FALSE
              }
              if ('vertical_text_column' %in% names(heat_map_i_def)) {
                param[['vertical_text_column']] <- heat_map_i_def[['vertical_text_column']]
              } else {
                param[['vertical_text_column']] <- ""
              }
              if ('vertical_text_size' %in% names(heat_map_i_def)) {
                param[['vertical_text_size']] <- as.numeric(heat_map_i_def[['vertical_text_size']])
              } else {
                param[['vertical_text_size']] <- 3
              }
              if ('vertical_text_offset' %in% names(heat_map_i_def)) {
                param[['vertical_text_offset']] <- as.numeric(heat_map_i_def[['vertical_text_offset']])
              } else {
                param[['vertical_text_offset']] <- 0.5
              }
              if ('vertical_text_color' %in% names(heat_map_i_def)) {
                param[['vertical_text_color']] <- heat_map_i_def[['vertical_text_color']]
              } else {
                param[['vertical_text_color']] <- "#000000"
              }

              # v109: Get colnames_angle (default 45)
              if ('colnames_angle' %in% names(heat_map_i_def)) {
                param[['colnames_angle']] <- as.numeric(heat_map_i_def[['colnames_angle']])
              } else {
                param[['colnames_angle']] <- 45
              }

              # S1.62dev: Get show_colnames flag (default TRUE)
              if ('show_colnames' %in% names(heat_map_i_def)) {
                param[['show_colnames']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_colnames']])
              } else {
                param[['show_colnames']] <- TRUE
              }

              # v105: Get row labels settings
              if ('show_row_labels' %in% names(heat_map_i_def)) {
                param[['show_row_labels']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_row_labels']])
              } else {
                param[['show_row_labels']] <- FALSE
              }
              if ('row_label_source' %in% names(heat_map_i_def)) {
                param[['row_label_source']] <- heat_map_i_def[['row_label_source']]
              } else {
                param[['row_label_source']] <- "colnames"
              }
              if ('row_label_font_size' %in% names(heat_map_i_def)) {
                param[['row_label_font_size']] <- as.numeric(heat_map_i_def[['row_label_font_size']])
              } else {
                param[['row_label_font_size']] <- 2.5
              }
              if ('custom_row_labels' %in% names(heat_map_i_def)) {
                param[['custom_row_labels']] <- heat_map_i_def[['custom_row_labels']]
              } else {
                param[['custom_row_labels']] <- ""
              }
              # v111: Get row label offset and alignment
              if ('row_label_offset' %in% names(heat_map_i_def)) {
                param[['row_label_offset']] <- as.numeric(heat_map_i_def[['row_label_offset']])
              } else {
                param[['row_label_offset']] <- 1.0
              }
              if ('row_label_align' %in% names(heat_map_i_def)) {
                param[['row_label_align']] <- heat_map_i_def[['row_label_align']]
              } else {
                param[['row_label_align']] <- "left"
              }
              # v108: Get label mapping
              if ('label_mapping' %in% names(heat_map_i_def)) {
                param[['label_mapping']] <- heat_map_i_def[['label_mapping']]
              } else {
                param[['label_mapping']] <- list()
              }

              # v117/v121: Get tip guide line settings (for discrete heatmaps)
              # v121: Added comprehensive debug logging
              debug_cat(paste0("\n=== v121: TIP GUIDE SETTINGS DEBUG (discrete) ===\n"))
              debug_cat(paste0("  heat_map_i_def names: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
              if ('show_guides' %in% names(heat_map_i_def)) {
                raw_val <- heat_map_i_def[['show_guides']]
                debug_cat(paste0("  show_guides raw value: ", raw_val, " (class: ", class(raw_val), ")\n"))
                param[['show_guides']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_guides']])
                debug_cat(paste0("  show_guides after conversion: ", param[['show_guides']], "\n"))
              } else {
                debug_cat(paste0("  show_guides NOT FOUND in heat_map_i_def\n"))
                param[['show_guides']] <- FALSE
              }
              if ('guide_color1' %in% names(heat_map_i_def)) {
                param[['guide_color1']] <- heat_map_i_def[['guide_color1']]
              } else {
                param[['guide_color1']] <- "#CCCCCC"
              }
              if ('guide_color2' %in% names(heat_map_i_def)) {
                param[['guide_color2']] <- heat_map_i_def[['guide_color2']]
              } else {
                param[['guide_color2']] <- "#EEEEEE"
              }
              if ('guide_alpha' %in% names(heat_map_i_def)) {
                param[['guide_alpha']] <- as.numeric(heat_map_i_def[['guide_alpha']])
              } else {
                param[['guide_alpha']] <- 0.3
              }
              if ('guide_width' %in% names(heat_map_i_def)) {
                param[['guide_width']] <- as.numeric(heat_map_i_def[['guide_width']])
              } else {
                param[['guide_width']] <- 0.5
              }
              if ('guide_linetype' %in% names(heat_map_i_def)) {
                param[['guide_linetype']] <- heat_map_i_def[['guide_linetype']]
              } else {
                param[['guide_linetype']] <- "solid"
              }

            } else {
              #print("AAAAAAAAAAAA")
              # print("is discrete false")
              ll <- "beige"
              mm <- "seashell2"
              hh <- "firebrick4"
              if ('man_define_colors' %in% names(heat_map_i_def)) {
                #print("IN")
                
                param['man_define_colors'] <- func.check.bin.val.from.conf(heat_map_i_def$man_define_colors)
              } else {
                
                param['man_define_colors'] <- FALSE
              }
              if (param['man_define_colors']== TRUE) {
                #print("innn")
                
                if ('color_scale_option' %in% names(heat_map_i_def)) {
                  #print("INNNN")
                  colors_man_list_from_yml <- heat_map_i_def[['color_scale_option']]
                  #print(colors_man_list_from_yml)
                  #print(length(colors_man_list_from_yml['color_scale_option']))
                  #print(length(colors_man_list_from_yml[['color_scale_option']]))
                  if (length(colors_man_list_from_yml) >3) {
                    colors_man_list_from_yml <- colors_man_list_from_yml[1:3]
                  }
                  if (length(colors_man_list_from_yml) ==3) {
                    ll <- colors_man_list_from_yml[1]
                    mm <- colors_man_list_from_yml[2]
                    hh <- colors_man_list_from_yml[3]
                  }
                }
                
                if ("limits" %in% names(heat_map_i_def)) {
                  
                  param[['limits']]<- heat_map_i_def[['limits']]
                } else {
                  
                  param[['limits']] <- NA
                }
              }
              param['low'] <- ll
              param['mid'] <- mm
              param['high'] <- hh
              param['midpoint']<- .02
              
              #print("param is")
              #print(param)
              # S1.62dev: Debug - show heat_map_i_def color values
              cat(file=stderr(), paste0("[DEBUG-COLOR-PARAM] heat_map_i_def names: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
              cat(file=stderr(), paste0("[DEBUG-COLOR-PARAM] heat_map_i_def$low = ", ifelse('low' %in% names(heat_map_i_def), heat_map_i_def$low, "NOT FOUND"), "\n"))
              cat(file=stderr(), paste0("[DEBUG-COLOR-PARAM] heat_map_i_def$mid = ", ifelse('mid' %in% names(heat_map_i_def), heat_map_i_def$mid, "NOT FOUND"), "\n"))
              cat(file=stderr(), paste0("[DEBUG-COLOR-PARAM] heat_map_i_def$high = ", ifelse('high' %in% names(heat_map_i_def), heat_map_i_def$high, "NOT FOUND"), "\n"))
              cat(file=stderr(), paste0("[DEBUG-COLOR-PARAM] heat_map_i_def$show_col_lines = ", ifelse('show_col_lines' %in% names(heat_map_i_def), heat_map_i_def$show_col_lines, "NOT FOUND"), "\n"))
              if ('low' %in% names(heat_map_i_def)) {
                param['low'] <-heat_map_i_def$low
              }
              if ('mid' %in% names(heat_map_i_def)) {
                param['mid'] <-heat_map_i_def$mid
              }
              if ('high' %in% names(heat_map_i_def)) {
                param['high'] <-heat_map_i_def$high
              }
              if ('midpoint' %in% names(heat_map_i_def)) {
                param['midpoint'] <-as.numeric(heat_map_i_def$midpoint)
              }
              # S1.62dev: Get use_midpoint flag for continuous heatmaps
              if ('use_midpoint' %in% names(heat_map_i_def)) {
                param[['use_midpoint']] <- func.check.bin.val.from.conf(heat_map_i_def[['use_midpoint']])
              } else {
                param[['use_midpoint']] <- FALSE
              }

              # v112: Get NA color for continuous heatmaps (default grey90)
              if ('na_color' %in% names(heat_map_i_def)) {
                param[['na_color']] <- heat_map_i_def[['na_color']]
              } else {
                param[['na_color']] <- "grey90"
              }

              # v104: Get per-heatmap distance for continuous heatmaps too
              if ('distance' %in% names(heat_map_i_def)) {
                param[['distance']] <- as.numeric(heat_map_i_def[['distance']])
              } else {
                param[['distance']] <- 0.02
              }

              # v105: Get per-heatmap height (default 0.8) for continuous heatmaps too
              if ('height' %in% names(heat_map_i_def)) {
                param[['height']] <- as.numeric(heat_map_i_def[['height']])
              } else {
                param[['height']] <- 0.8
              }

              # v111: Get per-heatmap row_height (default 1.0) for continuous heatmaps too
              if ('row_height' %in% names(heat_map_i_def)) {
                param[['row_height']] <- as.numeric(heat_map_i_def[['row_height']])
              } else {
                param[['row_height']] <- 1.0
              }

              # v111: Get grid settings for continuous heatmaps too
              if ('show_grid' %in% names(heat_map_i_def)) {
                param[['show_grid']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_grid']])
              } else {
                param[['show_grid']] <- FALSE
              }
              if ('grid_color' %in% names(heat_map_i_def)) {
                param[['grid_color']] <- heat_map_i_def[['grid_color']]
              } else {
                param[['grid_color']] <- "#000000"
              }
              if ('grid_size' %in% names(heat_map_i_def)) {
                param[['grid_size']] <- as.numeric(heat_map_i_def[['grid_size']])
              } else {
                param[['grid_size']] <- 0.5
              }

              # S1.62dev: Get row line settings for continuous heatmaps too
              if ('show_row_lines' %in% names(heat_map_i_def)) {
                param[['show_row_lines']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_row_lines']])
              } else {
                param[['show_row_lines']] <- FALSE
              }
              if ('row_line_color' %in% names(heat_map_i_def)) {
                param[['row_line_color']] <- heat_map_i_def[['row_line_color']]
              } else {
                param[['row_line_color']] <- "#000000"
              }
              if ('row_line_size' %in% names(heat_map_i_def)) {
                param[['row_line_size']] <- as.numeric(heat_map_i_def[['row_line_size']])
              } else {
                param[['row_line_size']] <- 0.5
              }
              # S1.62dev: Get column line settings for continuous heatmaps too
              if ('show_col_lines' %in% names(heat_map_i_def)) {
                param[['show_col_lines']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_col_lines']])
              } else {
                param[['show_col_lines']] <- FALSE
              }
              if ('col_line_color' %in% names(heat_map_i_def)) {
                param[['col_line_color']] <- heat_map_i_def[['col_line_color']]
              } else {
                param[['col_line_color']] <- "#000000"
              }
              if ('col_line_size' %in% names(heat_map_i_def)) {
                param[['col_line_size']] <- as.numeric(heat_map_i_def[['col_line_size']])
              } else {
                param[['col_line_size']] <- 0.5
              }

              # S1.62dev: Get vertical text labels settings for continuous heatmaps too
              if ('show_vertical_text' %in% names(heat_map_i_def)) {
                param[['show_vertical_text']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_vertical_text']])
              } else {
                param[['show_vertical_text']] <- FALSE
              }
              if ('vertical_text_column' %in% names(heat_map_i_def)) {
                param[['vertical_text_column']] <- heat_map_i_def[['vertical_text_column']]
              } else {
                param[['vertical_text_column']] <- ""
              }
              if ('vertical_text_size' %in% names(heat_map_i_def)) {
                param[['vertical_text_size']] <- as.numeric(heat_map_i_def[['vertical_text_size']])
              } else {
                param[['vertical_text_size']] <- 3
              }
              if ('vertical_text_offset' %in% names(heat_map_i_def)) {
                param[['vertical_text_offset']] <- as.numeric(heat_map_i_def[['vertical_text_offset']])
              } else {
                param[['vertical_text_offset']] <- 0.5
              }
              if ('vertical_text_color' %in% names(heat_map_i_def)) {
                param[['vertical_text_color']] <- heat_map_i_def[['vertical_text_color']]
              } else {
                param[['vertical_text_color']] <- "#000000"
              }

              # v109: Get colnames_angle (default 45) for continuous heatmaps too
              if ('colnames_angle' %in% names(heat_map_i_def)) {
                param[['colnames_angle']] <- as.numeric(heat_map_i_def[['colnames_angle']])
              } else {
                param[['colnames_angle']] <- 45
              }

              # S1.62dev: Get show_colnames flag (default TRUE) for continuous heatmaps too
              if ('show_colnames' %in% names(heat_map_i_def)) {
                param[['show_colnames']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_colnames']])
              } else {
                param[['show_colnames']] <- TRUE
              }

              # v105: Get row labels settings for continuous heatmaps too
              if ('show_row_labels' %in% names(heat_map_i_def)) {
                param[['show_row_labels']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_row_labels']])
              } else {
                param[['show_row_labels']] <- FALSE
              }
              if ('row_label_source' %in% names(heat_map_i_def)) {
                param[['row_label_source']] <- heat_map_i_def[['row_label_source']]
              } else {
                param[['row_label_source']] <- "colnames"
              }
              if ('row_label_font_size' %in% names(heat_map_i_def)) {
                param[['row_label_font_size']] <- as.numeric(heat_map_i_def[['row_label_font_size']])
              } else {
                param[['row_label_font_size']] <- 2.5
              }
              if ('custom_row_labels' %in% names(heat_map_i_def)) {
                param[['custom_row_labels']] <- heat_map_i_def[['custom_row_labels']]
              } else {
                param[['custom_row_labels']] <- ""
              }
              # v111: Get row label offset and alignment for continuous heatmaps too
              if ('row_label_offset' %in% names(heat_map_i_def)) {
                param[['row_label_offset']] <- as.numeric(heat_map_i_def[['row_label_offset']])
              } else {
                param[['row_label_offset']] <- 1.0
              }
              if ('row_label_align' %in% names(heat_map_i_def)) {
                param[['row_label_align']] <- heat_map_i_def[['row_label_align']]
              } else {
                param[['row_label_align']] <- "left"
              }
              # v108: Get label mapping for continuous heatmaps too
              if ('label_mapping' %in% names(heat_map_i_def)) {
                param[['label_mapping']] <- heat_map_i_def[['label_mapping']]
              } else {
                param[['label_mapping']] <- list()
              }

              # v117/v121: Get tip guide line settings (for continuous heatmaps)
              # v121: Added comprehensive debug logging
              debug_cat(paste0("\n=== v121: TIP GUIDE SETTINGS DEBUG (continuous) ===\n"))
              debug_cat(paste0("  heat_map_i_def names: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
              if ('show_guides' %in% names(heat_map_i_def)) {
                raw_val <- heat_map_i_def[['show_guides']]
                debug_cat(paste0("  show_guides raw value: ", raw_val, " (class: ", class(raw_val), ")\n"))
                param[['show_guides']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_guides']])
                debug_cat(paste0("  show_guides after conversion: ", param[['show_guides']], "\n"))
              } else {
                debug_cat(paste0("  show_guides NOT FOUND in heat_map_i_def\n"))
                param[['show_guides']] <- FALSE
              }
              if ('guide_color1' %in% names(heat_map_i_def)) {
                param[['guide_color1']] <- heat_map_i_def[['guide_color1']]
              } else {
                param[['guide_color1']] <- "#CCCCCC"
              }
              if ('guide_color2' %in% names(heat_map_i_def)) {
                param[['guide_color2']] <- heat_map_i_def[['guide_color2']]
              } else {
                param[['guide_color2']] <- "#EEEEEE"
              }
              if ('guide_alpha' %in% names(heat_map_i_def)) {
                param[['guide_alpha']] <- as.numeric(heat_map_i_def[['guide_alpha']])
              } else {
                param[['guide_alpha']] <- 0.3
              }
              if ('guide_width' %in% names(heat_map_i_def)) {
                param[['guide_width']] <- as.numeric(heat_map_i_def[['guide_width']])
              } else {
                param[['guide_width']] <- 0.5
              }
              if ('guide_linetype' %in% names(heat_map_i_def)) {
                param[['guide_linetype']] <- heat_map_i_def[['guide_linetype']]
              } else {
                param[['guide_linetype']] <- "solid"
              }
            }

            # S2.8: Add data_source and cnv_display_mode for RData detailed rendering
            if ('data_source' %in% names(heat_map_i_def)) {
              param[['data_source']] <- heat_map_i_def[['data_source']]
              cat(file=stderr(), paste0("[S2.8-PARAM] Setting data_source to: ", heat_map_i_def[['data_source']], "\n"))
            }
            if ('cnv_display_mode' %in% names(heat_map_i_def)) {
              param[['cnv_display_mode']] <- heat_map_i_def[['cnv_display_mode']]
              cat(file=stderr(), paste0("[S2.8-PARAM] Setting cnv_display_mode to: ", heat_map_i_def[['cnv_display_mode']], "\n"))
            }
            if ('cnv_height_scale' %in% names(heat_map_i_def)) {
              param[['cnv_height_scale']] <- as.numeric(heat_map_i_def[['cnv_height_scale']])
              cat(file=stderr(), paste0("[S2.8-PARAM] Setting cnv_height_scale to: ", heat_map_i_def[['cnv_height_scale']], "\n"))
            }

            heat_display_params_list[[indx_for_sav]] <- param
            # print("B12")
            
            
            
            heat_display_vec <- c(heat_display_vec, TRUE)
            heat_flag <- TRUE
            heat_map_title <- heat_map_i_def$title
            heat_map_title_list <- c(heat_map_title_list,heat_map_title)

            # S1.62dev: Check if this is an RData CNV heatmap
            if ('data_source' %in% names(heat_map_i_def) && heat_map_i_def$data_source == "rdata") {
              cat(file=stderr(), "\n[HEATMAP-RENDER] Processing RData CNV heatmap\n")
              debug_cat(paste0("\n=== S1.62dev: Processing RData CNV heatmap ===\n"))

              # Get the CNV matrix from the function parameter (passed from Shiny app)
              # Note: cnv_matrix doesn't serialize properly to YAML, so we pass it as a parameter
              if (!is.null(rdata_cnv_matrix)) {
                # Apply downsampling if specified in heatmap config
                # Default is 0 (no additional downsampling during render)
                # Check for new field name first, fall back to old name for backwards compatibility
                cnv_downsample <- if ('cnv_render_downsample' %in% names(heat_map_i_def)) {
                  as.numeric(heat_map_i_def$cnv_render_downsample)
                } else if ('cnv_downsample' %in% names(heat_map_i_def)) {
                  as.numeric(heat_map_i_def$cnv_downsample)
                } else {
                  0
                }
                cat(file=stderr(), paste0("[RDATA-CNV] Using render downsample factor: ", cnv_downsample, "\n"))
                cat(file=stderr(), paste0("[RDATA-CNV] heat_map_i_def fields: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
                cat(file=stderr(), paste0("[RDATA-CNV] cnv_render_downsample in heat_map_i_def: ",
                    ifelse('cnv_render_downsample' %in% names(heat_map_i_def),
                           as.character(heat_map_i_def$cnv_render_downsample), "NOT FOUND"), "\n"))

                # Apply downsampling
                cnv_data <- rdata_cnv_matrix
                if (cnv_downsample > 1 && ncol(cnv_data) > cnv_downsample) {
                  keep_cols <- seq(1, ncol(cnv_data), by = cnv_downsample)
                  cnv_data <- cnv_data[, keep_cols, drop = FALSE]
                  debug_cat(paste0("  Applied downsampling (factor ", cnv_downsample, "): ", ncol(cnv_data), " columns\n"))
                }

                # Apply WGD normalization if specified
                cnv_wgd_norm <- if ('cnv_wgd_norm' %in% names(heat_map_i_def)) {
                  func.check.bin.val.from.conf(heat_map_i_def$cnv_wgd_norm)
                } else {
                  FALSE
                }
                if (cnv_wgd_norm) {
                  cnv_data <- cnv_data / 2
                  debug_cat("  Applied WGD normalization (values / 2)\n")
                }

                debug_cat(paste0("  CNV matrix: ", nrow(cnv_data), " samples x ", ncol(cnv_data), " positions\n"))
                debug_cat(paste0("  Sample names (first 5): ", paste(head(rownames(cnv_data), 5), collapse=", "), "\n"))

                # The CNV matrix rownames should be sample IDs
                # We need to match them to tree tips
                # Get tree tip labels
                g_check <- suppressWarnings(ggtree(tree440))$data
                g_check_tip <- subset(g_check, isTip == TRUE)
                tree_tips <- g_check_tip$label

                # Handle NA labels in ggtree
                if (any(is.na(tree_tips))) {
                  tip_node_ids <- g_check_tip$node
                  ntips <- length(tree440$tip.label)
                  if (all(tip_node_ids >= 1 & tip_node_ids <= ntips)) {
                    tree_tips <- tree440$tip.label[tip_node_ids]
                  }
                }

                cat(file=stderr(), paste0("[HEATMAP-RENDER] Tree tips (first 5): ", paste(head(tree_tips, 5), collapse=", "), "\n"))
                debug_cat(paste0("  Tree tips (first 5): ", paste(head(tree_tips, 5), collapse=", "), "\n"))

                # Match CNV sample names to tree tips
                cnv_samples <- rownames(cnv_data)
                cat(file=stderr(), paste0("[HEATMAP-RENDER] CNV samples (first 5): ", paste(head(cnv_samples, 5), collapse=", "), "\n"))

                # S2.0-RDATA: Enhanced matching with multiple strategies
                # Strategy 0: Use user-selected mapping column (if specified)
                # Strategy 1: Use CSV data to create a mapping (tree tip ID -> sample name column)
                # Strategy 2: Direct match (tree tip == CNV sample name)
                # Strategy 3: Partial matching (CNV sample contains tree tip or vice versa)
                # Strategy 4: Numeric extraction - extract numbers from CNV samples and try to match
                # Strategy 5: Flexible partial matching with cleaned identifiers

                # Create mapping: find which CNV samples match which tree tips
                matched_cnv <- data.frame(matrix(NA, nrow = length(tree_tips), ncol = ncol(cnv_data)))
                rownames(matched_cnv) <- tree_tips
                colnames(matched_cnv) <- colnames(cnv_data)

                matches_found <- 0

                # S2.0-RDATA: Smart normalization function for comparing sample names
                # Handles common differences like "-" vs ".", leading "X", case differences
                smart_normalize <- function(x) {
                  x <- as.character(x)
                  x <- gsub("[._-]", "", x)   # Remove common separators
                  x <- gsub("^X", "", x)       # Remove leading X (R adds this to numeric column names)
                  tolower(x)
                }

                # S2.0-RDATA: Check if user specified a mapping column
                user_mapping_column <- NULL
                if ('rdata_mapping_column' %in% names(heat_map_i_def) &&
                    !is.null(heat_map_i_def$rdata_mapping_column) &&
                    heat_map_i_def$rdata_mapping_column != "") {
                  user_mapping_column <- heat_map_i_def$rdata_mapping_column
                  cat(file=stderr(), paste0("[HEATMAP-RENDER] User specified mapping column: '", user_mapping_column, "'\n"))
                }

                # First, try user-specified mapping column with smart matching
                csv_mapping <- NULL
                if (!is.null(user_mapping_column) && exists("readfile440") && !is.null(readfile440) && nrow(readfile440) > 0) {
                  cat(file=stderr(), paste0("[HEATMAP-RENDER] Using user-specified mapping column '", user_mapping_column, "' with smart matching\n"))

                  if (user_mapping_column %in% names(readfile440)) {
                    id_col_vals <- as.character(readfile440[[title.id]])
                    sample_col_vals <- as.character(readfile440[[user_mapping_column]])

                    cat(file=stderr(), paste0("[HEATMAP-RENDER] ID column '", title.id, "' values (first 5): ", paste(head(id_col_vals, 5), collapse=", "), "\n"))
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] Mapping column '", user_mapping_column, "' values (first 5): ", paste(head(sample_col_vals, 5), collapse=", "), "\n"))

                    # Pre-normalize CNV sample names for smart matching
                    cnv_samples_normalized <- smart_normalize(cnv_samples)
                    names(cnv_samples_normalized) <- cnv_samples  # Keep original names as names

                    for (i in seq_along(tree_tips)) {
                      tip <- tree_tips[i]
                      # Find this tree tip in the CSV ID column
                      csv_row <- which(id_col_vals == tip)
                      if (length(csv_row) > 0) {
                        # Get the corresponding sample name from the mapping column
                        sample_val <- sample_col_vals[csv_row[1]]
                        if (!is.na(sample_val) && sample_val != "") {
                          sample_val_normalized <- smart_normalize(sample_val)

                          # Try exact match first
                          if (sample_val %in% cnv_samples) {
                            matched_cnv[tip, ] <- as.numeric(cnv_data[sample_val, ])
                            matches_found <- matches_found + 1
                          } else {
                            # Try smart normalized match
                            match_idx <- which(cnv_samples_normalized == sample_val_normalized)
                            if (length(match_idx) > 0) {
                              matched_sample <- cnv_samples[match_idx[1]]
                              matched_cnv[tip, ] <- as.numeric(cnv_data[matched_sample, ])
                              matches_found <- matches_found + 1
                            } else {
                              # Try partial match - CNV sample contains the mapping value or vice versa
                              for (j in seq_along(cnv_samples)) {
                                cnv_s <- cnv_samples[j]
                                cnv_s_norm <- cnv_samples_normalized[j]
                                if (grepl(sample_val_normalized, cnv_s_norm, fixed = TRUE) ||
                                    grepl(cnv_s_norm, sample_val_normalized, fixed = TRUE)) {
                                  matched_cnv[tip, ] <- as.numeric(cnv_data[cnv_s, ])
                                  matches_found <- matches_found + 1
                                  break
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] After user-specified mapping (smart match): ", matches_found, " matches\n"))
                    csv_mapping <- list(column = user_mapping_column, count = matches_found, user_specified = TRUE)
                  } else {
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] WARNING: User-specified column '", user_mapping_column, "' not found in CSV\n"))
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] Available columns: ", paste(names(readfile440), collapse=", "), "\n"))
                  }
                }

                # If no user-specified column or it didn't work, try auto-detection
                if (is.null(csv_mapping) && exists("readfile440") && !is.null(readfile440) && nrow(readfile440) > 0) {
                  cat(file=stderr(), "[HEATMAP-RENDER] Auto-detecting mapping column in CSV...\n")
                  cat(file=stderr(), paste0("[HEATMAP-RENDER] CSV columns (all): ", paste(names(readfile440), collapse=", "), "\n"))
                  cat(file=stderr(), paste0("[HEATMAP-RENDER] ID column '", title.id, "' values (first 5): ", paste(head(as.character(readfile440[[title.id]]), 5), collapse=", "), "\n"))

                  # Look for a column that contains values matching CNV sample names
                  for (col_name in names(readfile440)) {
                    col_vals <- as.character(readfile440[[col_name]])
                    # Check if any CNV sample names appear in this column (exact match)
                    matches_in_col <- sum(cnv_samples %in% col_vals)
                    if (matches_in_col > 0) {
                      cat(file=stderr(), paste0("[HEATMAP-RENDER] Found potential mapping column '", col_name, "' with ", matches_in_col, " exact matches\n"))
                      if (is.null(csv_mapping) || matches_in_col > csv_mapping$count) {
                        csv_mapping <- list(column = col_name, count = matches_in_col)
                      }
                    }
                    # S2.0-RDATA: Also check for partial matches (column values contained in CNV sample names)
                    if (matches_in_col == 0 && col_name != title.id) {
                      partial_matches <- 0
                      unique_vals <- unique(col_vals[!is.na(col_vals) & col_vals != ""])
                      if (length(unique_vals) > 0 && length(unique_vals) < 500) {  # Only check if reasonable number of unique values
                        for (val in unique_vals) {
                          if (nchar(val) >= 3) {  # Only match if value is at least 3 chars
                            for (cnv_s in cnv_samples) {
                              if (grepl(val, cnv_s, fixed = TRUE)) {
                                partial_matches <- partial_matches + 1
                                break
                              }
                            }
                          }
                        }
                      }
                      if (partial_matches > 0) {
                        cat(file=stderr(), paste0("[HEATMAP-RENDER] Column '", col_name, "' has ", partial_matches, " values that appear as substrings in CNV sample names\n"))
                        cat(file=stderr(), paste0("[HEATMAP-RENDER]   Column sample values: ", paste(head(unique_vals, 5), collapse=", "), "\n"))
                        if (partial_matches > 5 && (is.null(csv_mapping) || partial_matches > csv_mapping$count)) {
                          csv_mapping <- list(column = col_name, count = partial_matches, partial = TRUE)
                        }
                      }
                    }
                  }

                  # Use the best auto-detected mapping column if found
                  if (!is.null(csv_mapping) && csv_mapping$count > 0) {
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] Using auto-detected CSV column '", csv_mapping$column, "' for mapping\n"))

                    # Create a mapping from tree tip (ID column) to sample name
                    id_col_vals <- as.character(readfile440[[title.id]])
                    sample_col_vals <- as.character(readfile440[[csv_mapping$column]])

                    cat(file=stderr(), paste0("[HEATMAP-RENDER] ID column '", title.id, "' values (first 5): ", paste(head(id_col_vals, 5), collapse=", "), "\n"))
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] Sample column '", csv_mapping$column, "' values (first 5): ", paste(head(sample_col_vals, 5), collapse=", "), "\n"))

                    is_partial <- !is.null(csv_mapping$partial) && csv_mapping$partial

                    for (i in seq_along(tree_tips)) {
                      tip <- tree_tips[i]
                      # Find this tree tip in the CSV ID column
                      csv_row <- which(id_col_vals == tip)
                      if (length(csv_row) > 0) {
                        # Get the corresponding sample name from the mapping column
                        sample_val <- sample_col_vals[csv_row[1]]

                        if (is_partial) {
                          # For partial matches, find the CNV sample that contains this value
                          for (cnv_s in cnv_samples) {
                            if (grepl(sample_val, cnv_s, fixed = TRUE)) {
                              matched_cnv[tip, ] <- as.numeric(cnv_data[cnv_s, ])
                              matches_found <- matches_found + 1
                              break
                            }
                          }
                        } else {
                          # Exact match
                          if (sample_val %in% cnv_samples) {
                            matched_cnv[tip, ] <- as.numeric(cnv_data[sample_val, ])
                            matches_found <- matches_found + 1
                          }
                        }
                      }
                    }
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] After auto-detected CSV mapping: ", matches_found, " matches\n"))
                  }
                }

                # If CSV mapping didn't find everything, try direct matching
                if (matches_found < length(tree_tips)) {
                  cat(file=stderr(), "[HEATMAP-RENDER] Trying direct and partial matching...\n")
                  for (tip in tree_tips) {
                    if (all(is.na(matched_cnv[tip, ]))) {  # Only if not already matched
                      # Try exact match
                      if (tip %in% cnv_samples) {
                        matched_cnv[tip, ] <- as.numeric(cnv_data[tip, ])
                        matches_found <- matches_found + 1
                      } else {
                        # Try partial matching - CNV sample contains tree tip or vice versa
                        for (cnv_s in cnv_samples) {
                          if (grepl(tip, cnv_s, fixed = TRUE) || grepl(cnv_s, tip, fixed = TRUE)) {
                            matched_cnv[tip, ] <- as.numeric(cnv_data[cnv_s, ])
                            matches_found <- matches_found + 1
                            break
                          }
                        }
                      }
                    }
                  }
                }

                # S2.0-RDATA: If still no matches, show diagnostic info to help user
                if (matches_found == 0) {
                  cat(file=stderr(), "\n[HEATMAP-RENDER] *** WARNING: NO MATCHES FOUND ***\n")
                  cat(file=stderr(), "[HEATMAP-RENDER] This means the CNV sample names don't match tree tip labels.\n")
                  cat(file=stderr(), "[HEATMAP-RENDER] To fix this, your CSV needs a column that maps tree tip IDs to CNV sample names.\n")
                  cat(file=stderr(), "[HEATMAP-RENDER] Tree tip labels (samples): ", paste(head(tree_tips, 10), collapse=", "), "...\n")
                  cat(file=stderr(), "[HEATMAP-RENDER] CNV sample names: ", paste(head(cnv_samples, 10), collapse=", "), "...\n")
                  cat(file=stderr(), "[HEATMAP-RENDER] Add a column to your CSV that contains the CNV sample names, with rows matching your ID column.\n\n")
                }

                cat(file=stderr(), paste0("[HEATMAP-RENDER] Matched ", matches_found, " out of ", length(tree_tips), " tree tips to CNV data\n"))
                debug_cat(paste0("  Matched ", matches_found, " out of ", length(tree_tips), " tree tips to CNV data\n"))

                # Use the matched CNV data as the heatmap dataframe
                df_heat_temp <- matched_cnv

                # S2.12: Apply per-cell WGD normalization if specified
                cnv_wgd_per_cell <- if ('cnv_wgd_per_cell' %in% names(heat_map_i_def)) {
                  func.check.bin.val.from.conf(heat_map_i_def$cnv_wgd_per_cell)
                } else {
                  FALSE
                }

                if (cnv_wgd_per_cell) {
                  cnv_wgd_column <- heat_map_i_def$cnv_wgd_column
                  debug_cat(paste0("  S2.12: Per-cell WGD enabled, column: ", ifelse(is.null(cnv_wgd_column), "NULL", cnv_wgd_column), "\n"))

                  if (!is.null(cnv_wgd_column) && cnv_wgd_column %in% names(readfile440)) {
                    # Get the WGD values from CSV
                    id_col_vals <- as.character(readfile440[[title.id]])
                    wgd_col_vals <- readfile440[[cnv_wgd_column]]

                    # Create a lookup: tree tip -> WGD value
                    wgd_lookup <- setNames(wgd_col_vals, id_col_vals)

                    # Apply per-cell normalization
                    cells_normalized <- 0
                    for (tip in rownames(df_heat_temp)) {
                      if (tip %in% names(wgd_lookup)) {
                        wgd_val <- wgd_lookup[[tip]]
                        if (!is.na(wgd_val) && wgd_val == 1) {
                          df_heat_temp[tip, ] <- df_heat_temp[tip, ] / 2
                          cells_normalized <- cells_normalized + 1
                        }
                      }
                    }
                    debug_cat(paste0("  S2.12: Applied per-cell WGD normalization to ", cells_normalized, " cells (where ", cnv_wgd_column, " = 1)\n"))
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] S2.12: Per-cell WGD normalization applied to ", cells_normalized, " cells\n"))
                  } else {
                    cat(file=stderr(), paste0("[HEATMAP-RENDER] S2.12: WARNING - WGD column '", cnv_wgd_column, "' not found in CSV. Available columns: ", paste(names(readfile440), collapse=", "), "\n"))
                  }
                }

                dxdf440_for_heat[[indx_for_sav]] <- df_heat_temp

                cat(file=stderr(), paste0("[HEATMAP-RENDER] CNV heatmap data: ", nrow(df_heat_temp), " x ", ncol(df_heat_temp), "\n"))
                debug_cat(paste0("  Created CNV heatmap data: ", nrow(df_heat_temp), " x ", ncol(df_heat_temp), "\n"))
              } else {
                debug_cat("  ERROR: rdata_cnv_matrix parameter is NULL - no RData CNV file loaded\n")
                heat_display_vec <- c(heat_display_vec, FALSE)
                next
              }
            } else {
              # Standard CSV-based heatmap - extract columns from readfile440
              acc_heat_list <- heat_map_i_def$according
              l_titles_for_heat <- c()
              according_from_script <- FALSE

              # print("B13")

              if ('according_from_script' %in% names(heat_map_i_def)) {

                according_from_script<- heat_map_i_def$according_from_script
              }

              if ('according_column_range' %in% names(heat_map_i_def)){
              start1 <- heat_map_i_def$according_column_range[1]
              end1 <- heat_map_i_def$according_column_range[2]
              l_titles_for_heat<- names(readfile440)[start1:end1]
              
              
            } else {
              # v57: DEBUG - trace column extraction
              debug_cat(paste0("\n=== v57: Extracting heatmap columns from 'according' ===\n"))
              debug_cat(paste0("  acc_heat_list length: ", length(acc_heat_list), "\n"))
              ind <-1
              for (j in acc_heat_list) {

                j1 <- names(j)
                debug_cat(paste0("  Column ", ind, ": j1=", j1, "\n"))

                ind<- ind+1
                j2<- j[[j1]]
                debug_cat(paste0("    j2 (column name)=", j2, "\n"))

                l_titles_for_heat <- c(l_titles_for_heat,j2)
              }
              debug_cat(paste0("  Final l_titles_for_heat: ", paste(l_titles_for_heat, collapse=", "), "\n"))
              # S2.3-DEBUG: Always show column extraction for verification
              cat(file=stderr(), paste0("[RENDER-DEBUG] Heatmap inx=", inx, " (indx_for_sav=", indx_for_sav,
                                        ") extracting columns: ", paste(l_titles_for_heat, collapse=", "), "\n"))
            }


            #print("B14")
            #names(readfile440)
            
            
            
            if (according_from_script== TRUE){
              # v53: print("In according_from_script")
              l_titles_for_heat<- heat_cnv_colnames[inx]
            }
            
            #print("l_titles_for_heat is")
            #print(list(l_titles_for_heat))
            #print("index 1 is")
            #idstart= which(l_titles_for_heat=='chrome_1total_length_abbr')
            
            #print(idstart)
            #print("index last is")
            #idend=  which(l_titles_for_heat=='chrome_22total_length_abbr')
            
            #print(idend)
            #print("CHECKKKKKKKKKKK")
            #print("readfile440 is")
            #print(head(readfile440))
            #print("colnames of readfile440 are")
            #print(names(readfile440))
            
            l_titles_for_heat <- as.character(l_titles_for_heat)

            # v57: DEBUG - show column validation
            debug_cat(paste0("\n=== v57: Validating heatmap columns ===\n"))
            debug_cat(paste0("  title.id: ", title.id, "\n"))
            debug_cat(paste0("  l_titles_for_heat: ", paste(l_titles_for_heat, collapse=", "), "\n"))
            debug_cat(paste0("  Available CSV columns: ", paste(head(names(readfile440), 10), collapse=", "), "...\n"))

            valid_columns <- c(title.id, l_titles_for_heat)
            valid_columns <- valid_columns[valid_columns %in% names(readfile440)]
            debug_cat(paste0("  Valid columns (after filtering): ", paste(valid_columns, collapse=", "), "\n"))
            debug_cat(paste0("  Number of valid columns: ", length(valid_columns), "\n"))

            # Select only the valid columns
            df_heat_temp <- readfile440[, valid_columns, drop = FALSE]
            debug_cat(paste0("  df_heat_temp dimensions: ", nrow(df_heat_temp), " rows x ", ncol(df_heat_temp), " cols\n"))
            
            #print("df_heat_temp is")
            #print(df_heat_temp)
            
            
            #df_heat_temp <-readfile440[, c(title.id,l_titles_for_heat)]
            
            
            # print("B15")       
            #print("df_heat_temp is")
            #print(head(df_heat_temp))
            
            duplicated_cols <- names(df_heat_temp)[duplicated(names(df_heat_temp))]
            
            # Print the duplicates (if any)
            if (length(duplicated_cols) > 0) {
              # v53: print(paste("Duplicate column names:", paste(duplicated_cols, collapse = ", ")))
            } else {
              # v53: print("No duplicate column names found.")
            }



            # v56b: Suppress harmless fortify warnings
            g_check <- suppressWarnings(ggtree(tree440))$data
            #print(print("B15a"))

            if (flag_print_tree_data== TRUE) {
              # v53: print("##tree data##")

              tib <- as_tibble(g_check)
              rows <- nrow(tib)
              # v53: print(tib,n=rows)

            }

            g_check_tip <- subset(g_check,isTip==TRUE)

            # v66: ROBUST FIX - Always use tree440$tip.label as source of truth
            # ggtree data frame can have NA labels in many scenarios (ladderizing, reordering, etc.)
            # The tree440$tip.label is indexed 1:Ntip corresponding to tip node IDs 1:Ntip
            ggtree_labels <- g_check_tip$label

            # v66: DEBUG - show initial state
            debug_cat(paste0("\n=== v66: TIP LABEL EXTRACTION DEBUG ===\n"))
            debug_cat(paste0("  ggtree_labels NA count: ", sum(is.na(ggtree_labels)), " out of ", length(ggtree_labels), "\n"))
            debug_cat(paste0("  tree440$tip.label sample: ", paste(head(tree440$tip.label, 5), collapse=", "), "\n"))
            debug_cat(paste0("  g_check_tip$node sample: ", paste(head(g_check_tip$node, 5), collapse=", "), "\n"))

            # v66: Check if ANY labels are NA or empty - if so, use tree440$tip.label
            # This is more robust than only checking if ALL are NA
            if (any(is.na(ggtree_labels)) || any(ggtree_labels == "")) {
              debug_cat(paste0("  v6: Found NA/empty labels in ggtree, using tree440$tip.label\n"))

              # For tip nodes, node IDs 1 to Ntip correspond directly to tree$tip.label indices
              tip_node_ids <- g_check_tip$node
              ntips <- length(tree440$tip.label)

              # Validate that node IDs are in valid range
              if (all(tip_node_ids >= 1 & tip_node_ids <= ntips)) {
                ggtree_labels <- tree440$tip.label[tip_node_ids]
                debug_cat(paste0("  v6: Successfully extracted labels from tree440$tip.label\n"))
              } else {
                # Fallback: use tip labels in their original order from tree440
                debug_cat(paste0("  v6: WARNING - node IDs out of range, using tree tip order\n"))
                # Order g_check_tip by y coordinate (visual order) and assign labels
                tip_order <- order(g_check_tip$y)
                ggtree_labels <- tree440$tip.label[tip_order]
              }

              debug_cat(paste0("  v6: After fix, ggtree_labels sample: ", paste(head(ggtree_labels, 5), collapse=", "), "\n"))
            }

            # v60: FIX - Logic was inverted! When id_tip_trim_flag == TRUE, apply trimming
            # When id_tip_trim_flag == FALSE (default), use full labels as-is
            if (id_tip_trim_flag == TRUE) {
              # Trimming ENABLED: extract substring from tip labels
              tip_list <- substr(ggtree_labels, id_tip_trim_start, id_tip_trim_end)
            } else {
              # Trimming DISABLED: use full tip labels
              tip_list <- ggtree_labels
            }
            # v53: print("TIP LIST CHECK")
            # v53: print("tip_list is")
            # v53: print(tip_list)
            # v53: print("df_heat_temp[[title.id]]) is")
            # v53: print(df_heat_temp[[title.id]])

            # S1.62dev: Check if this is an RData heatmap - skip CSV-specific matching code
            is_rdata_heatmap <- 'data_source' %in% names(heat_map_i_def) && heat_map_i_def$data_source == "rdata"

            if (!is_rdata_heatmap) {
            # v58: DEBUG - Show matching details (CSV only)
            debug_cat(paste0("\n=== v58: Matching heatmap data to tree tips ===\n"))
            debug_cat(paste0("  tip_list length: ", length(tip_list), "\n"))
            debug_cat(paste0("  tip_list sample: ", paste(head(tip_list, 5), collapse=", "), "\n"))
            debug_cat(paste0("  CSV ID column (", title.id, ") sample: ", paste(head(df_heat_temp[[title.id]], 5), collapse=", "), "\n"))

            # Check for matches before applying
            matches <- match(tip_list, df_heat_temp[[title.id]])
            num_matches <- sum(!is.na(matches))
            debug_cat(paste0("  Number of matches found: ", num_matches, " out of ", length(tip_list), " tips\n"))

            df_heat_temp <- df_heat_temp[matches,]
            debug_cat(paste0("  After match(): ", nrow(df_heat_temp), " rows\n"))

            ro= na.omit(df_heat_temp[[title.id]])
            df_heat_temp_filtered<- df_heat_temp[df_heat_temp[[title.id]] %in% (ro), ]

            df_heat_temp<- df_heat_temp_filtered
            debug_cat(paste0("  After filtering NAs: ", nrow(df_heat_temp), " rows\n"))
            # v53: print(df_heat_temp)
            # v53: print("df_heat_temp[[title.id]]) is")
            # v53: print(df_heat_temp[[title.id]])
            
            
            #print("B16")
            
            # v53: print("df_heat_temp is")
            # v53: print(head(df_heat_temp))
            
            
            # v53: print("id_tip_prefix is")
            #print(id_tip_prefix)
            #print("df_heat_temp[[title.id]] is")
            #print(df_heat_temp[[title.id]])
            #print("print(df_heat_temp[title.id]) is")
            #print(df_heat_temp[title.id])
            #    print("title.id is")
            #    print(title.id)
            #    print("names is")
            #    print(names(df_heat_temp))
            #print("len is")
            #print(length(df_heat_temp[[title.id]]))
            #print("len set is")
            #print(length(unique(df_heat_temp[[title.id]])))
            
            
            #question????YO YO
            
            # v58: FIXED - Validate data and ensure unique rownames before setting
            if (nrow(df_heat_temp) == 0) {
              # Skip this heatmap if no data
              debug_cat(paste0("  WARNING: No matching data for heatmap - skipping\n"))
              heat_display_vec <- c(heat_display_vec, FALSE)
              # v58: Reset heat_flag if no heatmaps have data
              if (indx_for_sav == 1) {
                heat_flag <- FALSE
                debug_cat(paste0("  Resetting heat_flag to FALSE\n"))
              }
              next
            }

            if (id_tip_trim_flag== FALSE) {
              #print("inn1")

              # v56b: Create row names and ensure they are unique
              new_rownames <- paste0(id_tip_prefix, df_heat_temp[[title.id]])
              # Handle duplicates by making unique
              if (anyDuplicated(new_rownames) > 0) {
                new_rownames <- make.unique(as.character(new_rownames))
              }
              # v56b: Validate length matches before setting
              if (length(new_rownames) == nrow(df_heat_temp)) {
                rownames(df_heat_temp) <- new_rownames
              }
              # v53: print("rownames(df_heat_temp) is")
              # v53: print(rownames(df_heat_temp))
              # v53: print("inn2")
            } else {
              # v53: print("else21")

              # v56b: Create row names and ensure they are unique
              new_rownames <- paste0("", df_heat_temp[[title.id]])
              # Handle duplicates by making unique
              if (anyDuplicated(new_rownames) > 0) {
                new_rownames <- make.unique(as.character(new_rownames))
              }
              # v56b: Validate length matches before setting
              if (length(new_rownames) == nrow(df_heat_temp)) {
                rownames(df_heat_temp) <- new_rownames
              }

              # v53: print("else22")
            }
            # v53: print("check l_titles_for_heat is")
            # v53: print(l_titles_for_heat)
            # v53: print("check names is")
            # v53: print(names(df_heat_temp))
            
            if (length(l_titles_for_heat)>1) {
              #print("Inn if length")
              
              df_heat_temp <- df_heat_temp[, c(l_titles_for_heat)]
            } else {
              
              df_heat_temp <- df_heat_temp[, c(l_titles_for_heat,l_titles_for_heat)]
              
              df_heat_temp <- df_heat_temp[ -c(2) ]
              
            }


            dxdf440_for_heat[[indx_for_sav]] <- df_heat_temp
            } # S1.62dev: End of if (!is_rdata_heatmap) block - CSV-specific matching code

            } # End of else block (CSV path) - S1.62dev
          } else {

            heat_display_vec <- c(heat_display_vec, FALSE)
          }
          temp <- dxdf440_for_heat[[indx_for_sav]]
          
          #print("B17")                 
          
          if ("with" %in% names(heat_map_i_def)) {
            
            
            
            with_title <- heat_map_i_def$with$title
            with_value <- heat_map_i_def$with$value
            
            
            if (( with_title %in% colnames(readfile440))) {
              
            } else {
              # v53: print(paste0("Error: missing  column in csv file ", with_title))
            }
            
            col_with <- readfile440[, with_title]
            
            
            
            
            for (col_index in colnames(temp)) {
              
              
              
              old_col <- temp[,col_index]
              
              if (length(col_with) < length(old_col) ) {
                gap <- length(old_col) - length(col_with)
                col_with <- c(col_with, rep(NA,gap))
              }
              
              new_col <-c()
              for (row_index in 1:length(old_col)) {
                val <- col_with[[row_index]]
                if (func.check.if.id.in.sub.class(val, with_title, with_value) == TRUE) {
                  new_col <- c(new_col, old_col[[row_index]])
                } else {
                  new_col <- c(new_col,NA)
                }
              }
              temp[,col_index] <- new_col
              
            }
            dxdf440_for_heat[[indx_for_sav]] <- temp
            
          }
          
        }
        
        
      }
      
      
      #print("B18")
      how_many_hi <- 0
      high_alpha_list <- NULL  # v140: Initialize before highlight block to prevent 'object not found' error


      FLAG_BULK_DISPLAY <- FALSE
      
      if ('highlight' %in% att1) {
        # v53: print("check high")
        hi_def <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$highlight
        hi<- hi_def$display
        # v53: print("hi_def=")
        # v53: print(hi_def)
        # v53: print("hi is")
        # v53: print(hi)
        
        if (func.check.bin.val.from.conf(hi)== FALSE) {
          # v53: print("no highlight")
          #no highlight
        } else {
          FLAG_BULK_DISPLAY <- TRUE
          # v53: print("in highlight")
          highlight.params.NEW <- func.make.highlight.params.NEW(yaml_file,title.id,ids_list,tree440,
                                                                 readfile440,hi_def,debug_mode,
                                                                 id_tip_trim_start,id_tip_trim_end,
                                                                 id_tip_trim_flag)
          
          
          # v53: print("highlight.params.NEW is")
          # v53: print(highlight.params.NEW)
          how_many_hi <- highlight.params.NEW$how_many_hi
          high_label_list<- highlight.params.NEW$high_label_list
          high_color_list<- highlight.params.NEW$high_color_list
          high_alpha_list<- highlight.params.NEW$high_alpha_list  # v139: Extract alpha list
          high_title_list<- highlight.params.NEW$high_title_list
          lists_list_hi<- highlight.params.NEW$lists_list_hi
          high_offset<- highlight.params.NEW$offset_hi
          high_vertical_offset<- highlight.params.NEW$vertical_offset_hi
          adjust_height_ecliplse<- highlight.params.NEW$adjust_height_ecliplse
          adjust_width_eclipse<- highlight.params.NEW$adjust_width_eclipse
          
          
        }
        
      }
      
      
      
      ######
      according_list <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'according'
      cls_list <- c()
      cls_renaming_list <- c()
      colors_scale1 <- c()
      
      for (j in 1:length(according_list)) {
        acc <- according_list[[j]]
        
        acc1 <- acc[[as.character(j)]]
        
        j_ch <- as.character(j)
        cls_list <- c(cls_list, j_ch)
        
        nams <- names(acc[[j_ch]])
        for (cri in c("color","display_name","title1","value1")) {
          if (cri %in% nams) {
            
          } else {
            # v53: print(paste0(paste0(paste0("Missing field in according ",acc),": "),cri))
          }
        }
        
        
        color <- acc[[j_ch]]$color
        
        
        display_name <- acc[[j_ch]]$display_name
        cls_renaming_list <- c(cls_renaming_list, display_name)
        
        
      }
      
      
      
      
      cls_num <- length(cls_list)
      
      
      if (debug_mode== TRUE) {
        # v53: print("##DEBUG##")
        # v53: print("disp_index is")
        # v53: print(disp_index)
        
        
        
      }

      # v56b: Suppress harmless fortify warnings
      tree_data <- suppressWarnings(ggtree(tree440))$data




      leaves_id_from_tree1 <- tree_data[tree_data$isTip==TRUE,'label']
      leaves_id_from_tree <- as.list(leaves_id_from_tree1['label'])$label
      
      
      
      
      list_id_by_class<- func.make.list_id_by_class(
        cls_num, cls_renaming_list,yaml_file,title.id,leaves_id_from_tree,readfile440,
        according_list,debug_mode,id_tip_trim_flag,id_tip_trim_start,id_tip_trim_end)
      
      #print("list_id_by_class is")
      #print(list_id_by_class)
      
      
      dx_rx_types1_short <- names(list_id_by_class)
      
      
      
      for (j in 1:length(according_list)) {
        acc <- according_list[[j]]
        
        
        j_ch <- as.character(j)
        
        
        display_name <- as.character(acc[[j_ch]]$display_name)
        
        
        
        
        color <- acc[[j_ch]]$color   
        
        # DEBUG: Print class info
        class_tip_count <- length(list_id_by_class[[display_name]])
        # v53: cat(file=stderr(), "[DEBUG] Building colors_scale1 - Class '", display_name, "': color=", color, ", tips=", class_tip_count, 
        #     ", adding color=", (class_tip_count > 1), "\n", sep="")
        
        if (length(list_id_by_class[[display_name]])>1 ) {
          colors_scale1 <- c(colors_scale1, color)
          # v53: cat(file=stderr(), "[DEBUG]   -> Added to colors_scale1. Current scale:", paste(colors_scale1, collapse=", "), "\n")
        }
        
      }
      
      
      
      
      non_cluster_color <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'non_cluster_color'
      colors_scale1 <- c(colors_scale1,non_cluster_color)
      
      
      na_name <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'na_name'
      if (is.null(na_name)) {
        na_name <-"NA"
      }
      
      #print("first_call function.create.cls.list")
      cls <- function.create.cls.list(leaves_id_from_tree,dx_rx_types1_short,list_id_by_class,na_name)
      
      #print("cls is")
      #print(cls)
      
      
      if (na_name %in% dx_rx_types1_short) {
        
      } else {
        dx_rx_types1_short <- c(dx_rx_types1_short,na_name)
        
        colors_scale1 <- c(colors_scale1,no_name_color)
      }
      
      # DEBUG: Print classification color mapping
      # v53: cat(file=stderr(), "\n🎨 === CLASSIFICATION COLOR DEBUG ===\n")
      # v53: cat(file=stderr(), "🎨 dx_rx_types1_short (class names):", paste(dx_rx_types1_short, collapse=", "), "\n")
      # v53: cat(file=stderr(), "🎨 colors_scale1:", paste(colors_scale1, collapse=", "), "\n")
      # v53: cat(file=stderr(), "🎨 Number of classes:", length(dx_rx_types1_short), "\n")
      # v53: cat(file=stderr(), "🎨 Number of colors:", length(colors_scale1), "\n")
      # v53: cat(file=stderr(), "🎨 list_id_by_class entries:\n")
      for (cls_name in names(list_id_by_class)) {
        # v53: cat(file=stderr(), "🎨   ", cls_name, ": ", length(list_id_by_class[[cls_name]]), " tips\n", sep="")
      }
      # v53: cat(file=stderr(), "🎨 =====================================\n\n")
      
      #question???
      #print("befff leaves_id_from_tree is")
      #print(leaves_id_from_tree)
      
      if (id_tip_trim_flag == TRUE) {
        # Trimming ENABLED: apply substring
        #print("INAAA - applying trimming")
        leaves_id_from_tree_num <- as.numeric(substring(leaves_id_from_tree,id_tip_trim_start))
        
      } else {
        # Trimming DISABLED: use full ID
        #print("INBBB - using full ID")
        #print(leaves_id_from_tree)
        leaves_id_from_tree_num <-  as.numeric(leaves_id_from_tree)
        
      }
      
      #print("leaves_id_from_tree_num is")
      #print(leaves_id_from_tree_num)
      #print("BEFORE check")
      #print("readfile440 is")
      #print(head(readfile440))
      #print("title.id is")
      #print(title.id)
      
      
      readfile440 <- readfile440[readfile440[[title.id]] %in% leaves_id_from_tree_num, ]
      
      #print("AFTER check")
      
      #print("readfile440 is")
      #print(head(readfile440))
      #print("______")
      
      co_nam <- colnames(readfile440)
      
      
      index_id <- which(co_nam==title.id)
      
      
      colnames(readfile440)[index_id] <- "Sample.Reads.ID"
      
      
      readfile440<- distinct(readfile440, Sample.Reads.ID, .keep_all = TRUE)
      colnames(readfile440)[index_id] <- title.id
      
      
      dxdf440_dataf<-readfile440[, c(title.id,title.id)]
      #print("dxdf440_dataf is")
      #print(dxdf440_dataf)
      
      
      leaves_id_ordered_for_df440 <- func.make.leaves_id_ordered_for_df440(leaves_id_from_tree1,dxdf440_dataf,title.id,
                                                                           id_tip_trim_flag,id_tip_prefix)
      
      #print("list_id_by_class is")
      #print(list_id_by_class)
      #print("before second function.create.cls.list")
      #print("leaves_id_ordered_for_df440 is")
      #print(leaves_id_ordered_for_df440)
      
      cls2 <- function.create.cls.list(leaves_id_ordered_for_df440,dx_rx_types1_short,list_id_by_class,na_name)
      
      #print("cls2 is")
      # print(cls2)
      
      dxdf440_dataf['Mapping'] <- cls2 #cls
      #print("dxdf440_dataf is")
      #print(dxdf440_dataf)
      
      #print("BEF")
      #print("colnames(dxdf440_dataf) is")
      #print(colnames(dxdf440_dataf))
      dxdf440_dataf = subset(dxdf440_dataf, select = c(title.id,'Mapping') )
      #print("AFTE")
      
      # v53: print("laderize_flag is")
      # v53: print(laderize_flag)
      
      if (laderize_flag == TRUE){
        rotate_flag_str <- paste0(rotate_flag_for_title,"_Ladderized")
      } else {
        rotate_flag_str <- rotate_flag_for_title
      }
      
      # v53: print("rotate_flag_str is")
      # v53: print(rotate_flag_str)
      
      out_start <- yaml_file[['Individual general definitions']]$out_file$'optional text at beggining'
      out_end <- yaml_file[['Individual general definitions']]$out_file$'optional text at end'
      file_type <- yaml_file[['Individual general definitions']]$out_file$'file_type'
      
      if (is.na(out_file_type)) {
        
      } else {
        file_type <- out_file_type
        # v53: print("use file type from out_file_type parameter: ", file_type)
      }
      
      if (flag_short_tips == TRUE) {
        tips_string <- "tips_normalized"
      } else {
        tips_string <- ""
      }
      
      classication2<- ""
      for (a in dx_rx_types1_short) {
        classication2 <- paste0(paste0(classication2,"_"),a)
      }
      
      
      file_name_end <- paste(out_start,classication2,tips_string,rotate_str, sep="_")
      
      
      
      if (add_date_to_text_flag ==TRUE) {
        str_date <- gsub(" " , "_", Sys.time()) 
        str_date <- gsub(":" , "_", str_date) 
        file_name_end <- paste(file_name_end,'_',str_date, sep = '')
      }
      
      
      out_file_path <- paste0(path_base,individual)
      
      
      path_split <- str_split(tree_path,'/')
      path_split <- path_split[[1]]
      
      
      s_nam1 <- path_split[length(path_split)]
      
      s_nam<- substr(s_nam1,1,nchar(s_nam1)-7)
      
      
      out_file_path <- paste0(out_file_path,paste0(paste0('_for_',s_nam),'__'))
      
      flag_replace_name <- yaml_file[['Individual general definitions']]$out_file$'replace name'$flag
      
      if (func.check.bin.val.from.conf(flag_replace_name) == TRUE) {
        out_file_path <- yaml_file[['Individual general definitions']]$out_file$'replace name'$name
        out_file_path<- paste0(path_base,out_file_path)
        out_file_path <- paste0(out_file_path,out_index)
        #out_file_path <- paste(out_file_path,'.',file_type, sep = '')
        
      }
      
      if (is.na(out_file_name_bash)) {
        # v53: print("use out_file_path from yaml")
      } else {
        # v53: print("use out_file_name_bash for output file name")
        out_file_path <- paste0(path_base,out_file_name_bash) 
        out_file_path <- paste0(out_file_path,"_")
        out_file_path <- paste0(out_file_path,disp_index)
        out_file_path <- paste0(out_file_path,"_")
        out_file_path <- paste0(out_file_path,out_index)
      }
      
      
      
      
      ooo <- str_split(out_file_path," ")
      
      
      
      oooo <- ooo[[1]]
      ind <-1
      for (letter in oooo) {
        if (letter %in% c(".","+","-","#")) {
          oooo[ind] <- ''
        }
        ind <- ind+1
      }
      
      
      ooo1 <- ""
      for (ind in 1:length(oooo)) {
        
        part <- oooo[ind]
        
        
        ooo1 <- paste0(ooo1,part)
        
      }
      out_file_path <- ooo1
      
      
      
      
      
      out_file_path <- paste(out_file_path,'.',file_type, sep = '')
      
      
      # v53: print("Output file path is")
      # v53: print(out_file_path)    
      
      
      
      
      ou<-0
      #print("before")
      
      #print(ou)
      #print(flag_calc_scores_for_tree)
      
      # === DEBUG CHECKPOINT 3: BEFORE INNER FUNCTION ===
      # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 3: BEFORE func.make.plot.tree.heat.NEW\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Passing node_number_font_size:", node_number_font_size, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Passing flag_display_nod_number_on_tree:", flag_display_nod_number_on_tree, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Passing highlight_manual_nodes:", highlight_manual_nodes, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Passing manual_nodes_to_highlight:", paste(manual_nodes_to_highlight, collapse=", "), "\n")
      # v53: debug_cat("================================================\n\n")
      
      ou <-     func.make.plot.tree.heat.NEW(
        tree440 = tree440,
        dx_rx_types1_short = dx_rx_types1_short,
        #classication =classication,
        list_id_by_class= list_id_by_class,
        dxdf440_dataf= dxdf440_dataf,
        title.id=title.id,
        FDR_perc=FDR_perc,
        no_name= no_name,
        rotate_flag = rotate_flag,
        rotation_params1= rotation_params1,
        rotation_params2= rotation_params2,
        flag_short_tips = flag_short_tips, 
        tips_length = tips_length,
        show_boot_flag = show_boot_flag,
        FLAG_BULK_DISPLAY = FLAG_BULK_DISPLAY,
        how_many_hi = how_many_hi,
        high_label_list = high_label_list,
        high_color_list = high_color_list,
        high_title_list = high_title_list,
        lists_list_hi = lists_list_hi,
        simulate.p.value= simulate.p.value,
        width= width,
        height= height,
        colors_scale1 = colors_scale1,
        out_file_path = out_file_path,
        edge_width_multiplier = edge_width_multiplier,
        size_tip_text= size_tip_text,
        size_font_legend_title= size_font_legend_title,
        size_font_legend_text = size_font_legend_text,
        size_font_legend_box= size_font_legend_box,
        labels_not_in_legend= labels_not_in_legend,
        no_name_color="black",
        debug_mode= debug_mode,
        man_adjust_elipse= man_adjust_elipse,
        man_multiply_elipse= man_multiply_elipse,
        man_adj_second_legend= man_adj_second_legend,
        man_space_second_legend= -0.02,
        laderize_flag=FALSE, 
        cls_renaming_list <- cls_renaming_list,
        title_flag = TRUE,
        title_replace = title_replace,
        flag_classification_format=FALSE,
        heat_flag= heat_flag,
        dxdf440_for_heat = dxdf440_for_heat,
        heat_map_title_list= heat_map_title_list,
        man_adjust_image_of_second_legend = man_adjust_image_of_second_legend,
        man_adj_heat_loc= man_adj_heat_loc,
        man_boot_x_offset= man_boot_x_offset,
        man_adj_heat_loc2= man_adj_heat_loc2,
        man_adj_heat_loc3= man_adj_heat_loc3,
        id_tip_trim_flag, 
        id_tip_trim_start,
        id_tip_trim_end, id_tip_prefix, 
        debug_print_data_tree,
        heat_display_vec,
        heat_display_params_list,
        man_multiply_second_legend,
        man_multiply_second_legend_text,
        size_font_heat_map_text,
        man_multiply_first_legend_text,
        man_multiply_first_legend_title_size,
        man_space_second_legend_multiplier,
        man_offset_for_highlight_legend_x,
        units_out,
        size_font_heat_map_legend= size_font_heat_map_legend,
        man_adjust_elipse_a,
        man_adjust_elipse_b,
        heat_maps_titles_angles_vector,
        boot_values,
        man_offset_second_legend,
        disp_index,
        list_nodes_to_rotate,
        flag_display_nod_number_on_tree,
        node_number_font_size,
        highlight_manual_nodes,
        manual_nodes_to_highlight,
        flag_calc_scores_for_tree,
        individual,
        width_heatmap,
        high_offset,
        high_vertical_offset,
        adjust_height_ecliplse,
        adjust_width_eclipse,
        high_alpha_list = high_alpha_list,  # v140: Pass transparency list
        flag_colnames,
        viridis_option_list,
        heat_legend_replace,
        tip_name_display_flag=tip_name_display_flag,
        flag_make_newick_file=flag_make_newick_file,
        bootstrap_label_size =bootstrap_label_size,
        heatmap_tree_distance = heatmap_tree_distance,
        heatmap_global_gap = heatmap_global_gap,  # v125: Gap between multiple heatmaps
        legend_settings = legend_settings,  # v136: Pass legend settings for highlight/bootstrap legends
        cached_p_list_of_pairs = cached_p_list_of_pairs,  # S2.0-PERF: Pass cached p-values
        cached_p_list_hash = cached_p_list_hash,  # S2.0-PERF: Pass cache hash for validation
        p_list_cache = p_list_cache,  # S2.7-PERF: Multi-entry cache
        heatmap_cache = heatmap_cache  # S2.9-PERF: Heatmap cache
      )
      # }

      # S2.0-PERF: Extract plot and cache data from new return structure (Option 3A)
      # func.make.plot.tree.heat.NEW now returns list(plot=..., cache_data=...)
      ou_result <- ou
      ou <- ou_result$plot  # Extract the plot object
      cache_data <- ou_result$cache_data  # Store cache data to return

      #print("ou is")
      #print(class(ou))
      #print(ou)


      out_list[[out_index]] <- out_file_path

      #a<- grid.arrange(ou,ou,ou)

      #ggsave("/home/dcsoft/s/yaara/chec/compare_trees_klein.pdf", plot=a, width = width, height = height*3, units = units_out, limitsize = FALSE)

      # v53: print("ou is")
      # v53: print(levels(ou$data$new_class))
      levels_base= levels(ou$data$new_class)
      #print(ou)
      #print("as.character(out_index) is")
      #print(as.character(out_index))
      #out_trees <- c(out_trees,ou)
      out_trees[[as.character(out_index)]] <-ou 
      out_index <- out_index+1    
      #close for of display options   
    }
    
    #a<- grid.arrange(out_trees[[as.character(2)]],out_trees[[as.character(2)]])
    
    #ggsave("/home/dcsoft/s/yaara/chec/compare_trees_klein.pdf", plot=a, width = width, height = height*3, units = units_out, limitsize = FALSE)
    
    
    
    
    # close for on trees structure     
  }
  
  
  
  #print(out_trees)
  # v53: print("compare_two_trees is")
  # v53: print(compare_two_trees)
  if(compare_two_trees== TRUE) {
    # v53: print("check")
    # v53: print("trees_to_compare is")
    # v53: print(trees_to_compare)
    # v53: print("out_trees is")
    # v53: print(out_trees)
    #print("classification_in_compare is")
    #print(classification_in_compare)
    # v53: print("out_file_path is")
    # v53: print(out_file_path)
    #print("tree_path_list is")
    #print(tree_path_list)
    #print("width is")
    #print(width)
    #print("height is")
    #print(height)
    #print("units_out is")
    #print(units_out)
    #print("path_base is")
    #print(path_base)
    #print("file_type is")
    #print(file_type)
    #print("flag_print_tree_data is")
    #print(flag_print_tree_data)
    #print("flag_display_nod_number_on_tree is")
    #print(flag_display_nod_number_on_tree)
    #print("rotate_compared_tree is")
    #print(rotate_compared_tree)
    #print("how_many_rotation_vec is")
    #print(how_many_rotation_vec)
    #print("trees_to_compare is")
    #print(trees_to_compare)
    #out_tree_A= out_tree[as.character(trees_to_compare[1])]
    #print("out_tree_A is")
    #print(out_tree_A)
    #out_tree_B= out_tree[as.character(trees_to_compare[2])]
    #print("out_tree_B is")
    #print(out_tree_B)
    func.make.compare_fig(trees_to_compare=trees_to_compare,
                          out_trees=out_trees,
                          classification_in_compare=classification_in_compare,
                          out_file_path=out_file_path,
                          tree_path_list=tree_path_list,
                          width=width,
                          height=height,
                          units_out=units_out,
                          path_base=path_base,
                          file_type=file_type,
                          flag_print_tree_data= flag_print_tree_data,
                          flag_display_nod_number_on_tree=flag_display_nod_number_on_tree,
                          rotate_compared_tree= rotate_compared_tree,
                          PERM_BOUND=7,
                          how_many_rotation_vec=how_many_rotation_vec,
                          flag_handle_big_branch=2,
                          colors_scale1=colors_scale1,
                          levels_base=levels_base)
  }
  
  
  out_file_path <- out_trees

  # S2.0-PERF: Return both plots and cache data for two-tier caching (Option 3A)
  # cache_data is set when func.make.plot.tree.heat.NEW is called above
  # If cache_data was never set (e.g., error path), use NULL
  if (!exists("cache_data") || is.null(cache_data)) {
    cache_data <- NULL
  }

  return(list(
    plots = out_trees,
    cache_data = cache_data
  ))
  #close func
}




##### part 4


# ============================================================================
# HEATMAP INTEGRATION FUNCTIONS
# ============================================================================



# Main function for creating the plot tree with heatmap
func.make.plot.tree.heat.NEW <- function(tree440, dx_rx_types1_short, list_id_by_class, dxdf440_dataf, 
                                         title.id, FDR_perc, no_name, rotate_flag, rotation_params1,
                                         rotation_params2, flag_short_tips, tips_length, show_boot_flag,
                                         FLAG_BULK_DISPLAY, how_many_hi = 0, high_label_list,
                                         high_color_list, high_title_list, lists_list_hi,
                                         simulate.p.value, width, height, colors_scale1,
                                         out_file_path, edge_width_multiplier = 1,
                                         size_tip_text = 3, size_font_legend_title = 30,
                                         size_font_legend_text = 20, size_font_legend_box = 15,
                                         labels_not_in_legend, no_name_color, debug_mode,
                                         man_adjust_elipse = 0, man_multiply_elipse,
                                         man_adj_second_legend = 0.15, man_space_second_legend = -0.02,
                                         laderize_flag = FALSE, cls_renaming_list, title_flag,
                                         title_replace, flag_classification_format, heat_flag = FALSE,
                                         dxdf440_for_heat, heat_map_title_list = NA,
                                         man_adjust_image_of_second_legend, man_adj_heat_loc,
                                         man_boot_x_offset, man_adj_heat_loc2, man_adj_heat_loc3,
                                         id_tip_trim_flag, id_tip_trim_start, id_tip_trim_end, id_tip_prefix,
                                         debug_print_data_tree, heat_display_vec, heat_display_params_list,
                                         man_multiply_second_legend, man_multiply_second_legend_text,
                                         size_font_heat_map_text, man_multiply_first_legend_text,
                                         man_multiply_first_legend_title_size, man_space_second_legend_multiplier,
                                         man_offset_for_highlight_legend_x, units_out, size_font_heat_map_legend,
                                         man_adjust_elipse_a, man_adjust_elipse_b, heat_maps_titles_angles_vector,
                                         boot_values, man_offset_second_legend, disp_index,
                                         list_nodes_to_rotate, flag_display_nod_number_on_tree,
                                         node_number_font_size = 3.5,
                                         highlight_manual_nodes = FALSE,
                                         manual_nodes_to_highlight = NA,
                                         flag_calc_scores_for_tree, individual, width_heatmap = NA,
                                         high_offset = 0, high_vertical_offset = 0, adjust_height_ecliplse, adjust_width_eclipse,
                                         high_alpha_list = NULL,  # v140: Added for transparency control
                                         flag_colnames, viridis_option_list, heat_legend_replace = NA,
                                         tip_name_display_flag = TRUE,
                                         flag_make_newick_file=FALSE,
                                         bootstrap_label_size = 1.5,  # v129: Reduced from 3.5 for smaller default legend
                                         heatmap_tree_distance = 0.02,
                                         heatmap_global_gap = 0.05,  # v125: Gap between multiple heatmaps
                                         legend_settings = NULL,  # v136: Legend settings for highlight/bootstrap legends
                                         cached_p_list_of_pairs = NULL,  # S2.0-PERF: Cached p-values (Option 3A)
                                         cached_p_list_hash = NULL,      # S2.0-PERF: Hash for cache validation
                                         p_list_cache = list(),          # S2.7-PERF: Multi-entry cache
                                         heatmap_cache = list()) {       # S2.9-PERF: Heatmap cache

  # === DEBUG CHECKPOINT 4: INNER FUNCTION ENTRY ===
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 4: func.make.plot.tree.heat.NEW ENTRY\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â node_number_font_size:", node_number_font_size, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â flag_display_nod_number_on_tree:", flag_display_nod_number_on_tree, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â highlight_manual_nodes:", highlight_manual_nodes, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â manual_nodes_to_highlight:", paste(manual_nodes_to_highlight, collapse=", "), "\n")
  # v53: debug_cat("================================================\n\n")
  
  # === PROFILING: func.make.plot.tree.heat.NEW ===
  .prof_func_start <- Sys.time()

  if (debug_mode == TRUE) {
    # v53: print("In func.make.plot.tree.HEAT")
  }

  # Prepare data structures
  dx_rx_types1 <- dx_rx_types1_short
  # v53: print("A")
  # v53: print(dx_rx_types1)
  
  # v56b: Wrap in suppressWarnings to suppress harmless ggtree/ggplot2 fortify warnings
  .prof_section_start <- Sys.time()
  pr440 <- suppressWarnings(ggtree(tree440))
  cat(file=stderr(), sprintf("[PROF-TREE] ggtree(tree440) #1: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  d440 <- pr440$data
  cc_tipss <- func.create.cc_tipss(d440)
  cc_nodss <- func.create.cc_nodss(d440)
  cc_totss <- func.create.cc_totss(d440)
  
  # v53: print(cc_tipss)
  
  nods_num <- length(cc_nodss)
  tips_num <- length(cc_tipss)
  tree_size <- nods_num + tips_num 
  #print("tree_size is")
  #print(tree_size)
  
  if (debug_mode == TRUE) {
    #print("tree data is")
    #print(d440)
  }
  
  if (debug_print_data_tree == TRUE) {
    # v53: print(tree_size)
    # v53: print(d440)
  }
  
  s <- subset(d440, isTip == "TRUE")
  
  # Create node classes
  list_node_by_class <- func.create.list_node_by_class(dx_rx_types1_short, list_id_by_class, dxdf440_dataf,
                                                       title.id, tree_size, d440, cc_totss, debug_mode,
                                                       id_tip_trim_flag, id_tip_prefix)
  
  
  # v53: print("C")
  list_rename_by_class <- list_node_by_class
  
  y_off_base <- -8
  yet_another_multiplier <- 0 
  
  # Create tree with group information
  .prof_section_start <- Sys.time()
  tree_with_group <- ggtree::groupOTU(tree440, list_node_by_class)
  cat(file=stderr(), sprintf("[PROF-TREE] groupOTU #1: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  
  # Prepare subframes
  subframe_of_nodes <- d440[d440$isTip == "FALSE", ]
  cc_nodss90 <- func.create.cc_nodss90(subframe_of_nodes)
  cc_nodss80 <- func.create.cc_nodss80(subframe_of_nodes)
  cc_nodss70 <- func.create.cc_nodss70(subframe_of_nodes)
  cc_tips <- subframe_of_nodes$node
  
  # Create array of node counts by type
  df_count_FULL_tree_populations <- data.frame(
    idx = 1:length(dx_rx_types1_short),
    type = dx_rx_types1_short,
    count = rep(0, length(dx_rx_types1_short))
  )
  
  # v53: print("D")
  
  for (opt in dx_rx_types1) {
    indx <- which(dx_rx_types1_short == opt)
    df_count_FULL_tree_populations[indx, 'count'] <- length(list_id_by_class[[opt]]) - 1
  }
  
  # Group tree by class
  .prof_section_start <- Sys.time()
  tree_with_group_CPY <- ggtree::groupOTU(tree440, list_rename_by_class)
  cat(file=stderr(), sprintf("[PROF-TREE] groupOTU #2: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  # v56b: Wrap in suppressWarnings to suppress harmless fortify warnings
  .prof_section_start <- Sys.time()
  levels_groups <- levels(suppressWarnings(ggtree(tree_with_group_CPY))$data$group)
  cat(file=stderr(), sprintf("[PROF-TREE] ggtree for levels: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  # v53: print("E")

  # Create tree with coloring
  # v56b: Wrap in suppressWarnings to suppress harmless fortify warnings
  .prof_section_start <- Sys.time()
  tree_TRY <- suppressWarnings(ggtree(tree_with_group_CPY, aes(color = new_class, size = p_val_new),
                     ladderize = laderize_flag))
  cat(file=stderr(), sprintf("[PROF-TREE] ggtree with aes: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  
  test_fig <- tree_TRY
  
  # Assign colors to nodes
  new_colors_list <- func.create.new_colors_list(FDR_perc, tree_TRY, tree_with_group, no_name, tree_size)
  
  if (no_name %in% new_colors_list) {
    flag_no_name <- 1
  } else {
    flag_no_name <- 0
  }
  # v53: print("F")
  
  tree_TRY$data$new_class <- new_colors_list
  
  # Define p-value categories
  op_list <- c(
    paste0("p>", FDR_perc), 
    paste0(paste0(0.5 * FDR_perc, "<p<="), FDR_perc),
    paste0("p<=", 0.5 * FDR_perc)
  )
  
  #print("G")

  # S2.0-PERF: Two-tier caching for p_list_of_pairs (Option 3A)
  # S2.7-PERF: Multi-entry cache - check p_list_cache first for instant switching
  # Compute current hash to check if cache is valid
  .prof_section_start <- Sys.time()
  current_p_list_hash <- func.create.p_list_cache_hash(
    tree440, list_id_by_class, dx_rx_types1_short, FDR_perc, simulate.p.value
  )

  # S2.7-PERF: Check multi-entry cache FIRST (allows instant switching between classifications)
  multi_cache_hit <- !is.null(p_list_cache) &&
                     length(p_list_cache) > 0 &&
                     (current_p_list_hash %in% names(p_list_cache))

  # Fall back to legacy single-entry cache
  legacy_cache_hit <- !multi_cache_hit &&
                      !is.null(cached_p_list_hash) &&
                      !is.null(cached_p_list_of_pairs) &&
                      (cached_p_list_hash == current_p_list_hash)

  if (multi_cache_hit) {
    # S2.7-PERF: Multi-cache hit - instant retrieval from previously computed classification
    p_list_of_pairs <- p_list_cache[[current_p_list_hash]]
    cat(file=stderr(), sprintf("[PERF-CACHE] MULTI-CACHE HIT (hash: %s, %d entries cached): %.3f sec\n",
                               substr(current_p_list_hash, 1, 8), length(p_list_cache),
                               as.numeric(Sys.time() - .prof_section_start)))
  } else if (legacy_cache_hit) {
    # S2.0-PERF: Legacy cache hit - use single cached p-values
    p_list_of_pairs <- cached_p_list_of_pairs
    cat(file=stderr(), sprintf("[PERF-CACHE] Using cached p_list_of_pairs (hash: %s): %.3f sec\n",
                               substr(current_p_list_hash, 1, 8), as.numeric(Sys.time() - .prof_section_start)))
  } else {
    # Cache miss - need to recalculate
    if (!is.null(cached_p_list_hash) && cached_p_list_hash != current_p_list_hash) {
      cat(file=stderr(), sprintf("[PERF-CACHE] Cache miss - hash %s not in cache (%d entries)\n",
                                 substr(current_p_list_hash, 1, 8), length(p_list_cache)))
    } else {
      cat(file=stderr(), "[PERF-CACHE] No cache available - computing p_list_of_pairs\n")
    }
    p_list_of_pairs <- func.create.p_list_of_pairs(
      list_node_by_class, d440, dx_rx_types1_short,
      cc_nodss, tree_with_group, FDR_perc, tree, cc_tipss,
      tree_TRY, tree_size, no_name, simulate.p.value
    )
    cat(file=stderr(), sprintf("[PERF-CACHE] Computed p_list_of_pairs (new hash: %s): %.3f sec\n",
                               substr(current_p_list_hash, 1, 8), as.numeric(Sys.time() - .prof_section_start)))
  }

  .prof_section_start <- Sys.time()
  p_PAIRS_pval_list <- func.create.p_val_list_FROM_LIST(FDR_perc, tree_TRY, p_list_of_pairs, op_list)
  cat(file=stderr(), sprintf("[PROF-TREE] p_val_list creation: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  
  # Calculate scores for tree if requested
  if (flag_calc_scores_for_tree == TRUE) {
    list_scores <- func.calc.p_val_for_root(
      list_node_by_class, d440, dx_rx_types1_short,
      cc_nodss, tree_with_group, FDR_perc, tree, cc_tipss,
      tree_TRY, tree_size, no_name, simulate.p.value, cc_tips
    ) 
    
    p_score_root <- list_scores[1]
    ari_score <- list_scores[2]
  }
  #print("H")
  
  # Assign p-values to tree
  tree_TRY$data$p_val_new <- p_PAIRS_pval_list 
  tree_TRY$data$p_val_new <- factor(tree_TRY$data$p_val_new, levels = op_list)
  
  id_root <- which(tree_TRY$data$node == tree_TRY$data$parent)
  p_score <- tree_TRY$data$p_val_new[id_root]
  
  # Handle highlighting if requested
  if (FLAG_BULK_DISPLAY == TRUE) {
    for (index_high in 1:how_many_hi) {
      list_high_for <- lists_list_hi[[index_high]]
      
      if (index_high == 1) {
        tree_TRY$data$'high1' <- list_high_for 
      } else if (index_high == 2) {
        tree_TRY$data$'high2' <- list_high_for 
      } else if (index_high == 3) {
        tree_TRY$data$'high3' <- list_high_for 
      } else {
        stop("Too many highlight options. Only 3 are supported")
      }
    }
  }
  
  # v53: print("I")
  
  tree_TRY1 <- tree_TRY
  
  # Handle tree rotation if requested
  tree_newick <- tree_TRY1
  
  if (rotate_flag %in% c("RX_first", "FRAC_first")) {
    list_weight_dx.rx <- rotation_params1[['list_weight_dx.rx']]
    list_weight_frac <- rotation_params2[['list_weight_dx.rx']]
    TREE_OTU_dx.rx <- rotation_params1[['TREE_OTU_dx.rx']]
    TREE_OTU_frac <- rotation_params2[['TREE_OTU_dx.rx']]
  }
  
  # Rotate specific nodes if requested
  tree_TRY1 <- func.rotate.specific.nodes(tree_TRY, list_nodes_to_rotate)
  
  # Apply rotation based on weights if requested
  if (rotate_flag == "RX_first") {
    tree_TRY2 <- func.rotate.tree.based.on.weights(
      tree_TRY1, list_weight_dx.rx, list_weight_frac,
      TREE_OTU_dx.rx, TREE_OTU_frac, tree_size
    )
  } else if (rotate_flag == "FRAC_first") {
    tree_TRY2 <- func.rotate.tree.based.on.weights(
      tree_TRY1, list_weight_frac, list_weight_dx.rx,
      TREE_OTU_frac, TREE_OTU_dx.rx, tree_size
    )
  } else {
    tree_TRY2 <- tree_TRY1
  }
  
  #print("J")
  
  tree_newick <- tree_TRY2
  
  # Handle short tips if requested
  temp <- tree_TRY
  pr440_short_tips_TRY <- tree_TRY2
  
  if (flag_short_tips == TRUE) {
    for (i in cc_tipss) {
      par <- pr440_short_tips_TRY$data$parent[i]
      parX <- pr440_short_tips_TRY$data$x[par]
      pr440_short_tips_TRY$data[pr440_short_tips_TRY$data$node[i], "x"] <- parX + tips_length
    }
  } 
  
  # Save rotated tree to newick file if requested
  if (flag_make_newick_file == TRUE) {
    phylo_tree <- as.phylo(tree_newick)
    pa <- paste0(paste0(path_out_newick, individual), "_rotated.newick")
    write.tree(phylo_tree, file = pa)
  }
  
  
  # v53: print("K")
  # Prepare for bootstrap values
  zs <- rep(0, tree_size)
  base <- rep("black", tree_size)
  pr440_short_tips_TRY$data$boot_val <- zs
  pr440_short_tips_TRY$data$boot_val_rounded <- zs#
  pr440_short_tips_TRY$data$color <- base
  space_vec <- rep("", tree_size)
  
  round_to <- 5
  list_boot_display <- c('raw', 'numbered_colored', 'numbered_color', 'percentage')
  
  if (boot_values$'format' %in% list_boot_display) {
    if ('param' %in% names(boot_values)) {
      round_to <- strtoi(boot_values$'param')
    } 
  }
  
  if (boot_values$'format' %in% c("numbered_colored", "numbered_color", 'percentage', 'percentage_color')) {
    round_to <- 1
  }
  #print("L")
  
  # Prepare for color-coded bootstrap values if requested
  if (boot_values$'format' %in% c("numbered_color", "percentage_color")) { 
    pr440_short_tips_TRY$data$boot_val_rounded0 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded1 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded2 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded3 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded4 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded5 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded6 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded7 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded8 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded9 <- space_vec
    pr440_short_tips_TRY$data$boot_val_rounded10 <- space_vec
  }
  
  jmp <- 1 / 10^round_to
  jmp_vector <- seq(0, 1, by = jmp)
  
  # Initialize bootstrap values
  for (i in 1:tree_size) {
    pr440_short_tips_TRY$data$boot_val_rounded[i] <- ""
  }
  
  #print("M")
  
  # Calculate bootstrap values if requested
  if (show_boot_flag == TRUE) {
    for (i in 1:tree_size) {
      if (i %in% cc_tipss) {
        pr440_short_tips_TRY$data$boot_val[i] <- 0
        pr440_short_tips_TRY$data$boot_val_rounded[i] <- ""
      } else {
        pr440_short_tips_TRY$data$boot_val[i] <- pr440_short_tips_TRY$data$label[i]
        if (pr440_short_tips_TRY$data$boot_val[i] == "") {
          pr440_short_tips_TRY$data$boot_val_rounded[i] <- 1
        } else {
          pr440_short_tips_TRY$data$boot_val_rounded[i] <- round(as.numeric(pr440_short_tips_TRY$data$boot_val[i]), digits = round_to)
        }
      }
    }
  }
  
  # Prepare bootstrap display values
  if (show_boot_flag == TRUE) {
    if (boot_values$'format' == "numbered_color") {
      for (i in 1:tree_size) {
        if (!i %in% cc_tipss) {
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0) {
            pr440_short_tips_TRY$data$boot_val_rounded0[i] <- 0
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.1) {
            pr440_short_tips_TRY$data$boot_val_rounded1[i] <- 0.1
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.2) {
            pr440_short_tips_TRY$data$boot_val_rounded2[i] <- 0.2
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.3) {
            pr440_short_tips_TRY$data$boot_val_rounded3[i] <- 0.3
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.4) {
            pr440_short_tips_TRY$data$boot_val_rounded4[i] <- 0.4
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.5) {
            pr440_short_tips_TRY$data$boot_val_rounded5[i] <- 0.5
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.6) {
            pr440_short_tips_TRY$data$boot_val_rounded6[i] <- 0.6
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.7) {
            pr440_short_tips_TRY$data$boot_val_rounded7[i] <- 0.7
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.8) {
            pr440_short_tips_TRY$data$boot_val_rounded8[i] <- 0.8
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.9) {
            pr440_short_tips_TRY$data$boot_val_rounded9[i] <- 0.9
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 1) {
            pr440_short_tips_TRY$data$boot_val_rounded10[i] <- 1
          }
        } 
      }  
    }
    
    # v53: print("N")
    
    if (boot_values$'format' == "percentage_color") {
      for (i in 1:tree_size) {
        if (!i %in% cc_tipss) {
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0) {
            pr440_short_tips_TRY$data$boot_val_rounded0[i] <- 0
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.1) {
            pr440_short_tips_TRY$data$boot_val_rounded1[i] <- 10
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.2) {
            pr440_short_tips_TRY$data$boot_val_rounded2[i] <- 20
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.3) {
            pr440_short_tips_TRY$data$boot_val_rounded3[i] <- 30
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.4) {
            pr440_short_tips_TRY$data$boot_val_rounded4[i] <- 40
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.5) {
            pr440_short_tips_TRY$data$boot_val_rounded5[i] <- 50
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.6) {
            pr440_short_tips_TRY$data$boot_val_rounded6[i] <- 60
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.7) {
            pr440_short_tips_TRY$data$boot_val_rounded7[i] <- 70
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.8) {
            pr440_short_tips_TRY$data$boot_val_rounded8[i] <- 80
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 0.9) {
            pr440_short_tips_TRY$data$boot_val_rounded9[i] <- 90
          }
          if (pr440_short_tips_TRY$data$boot_val_rounded[i] == 1) {
            pr440_short_tips_TRY$data$boot_val_rounded10[i] <- 100
          }
        } 
      }  
    }
  }
  
  # Prepare class list
  cls_renaming_list_with_no_name <- c(cls_renaming_list, no_name)
  
  if (missing(how_many_hi)) {
    how_many_hi <- 0
  }
  
  if (flag_no_name == 1) {
    fin_color_list <- cls_renaming_list_with_no_name
  } else {
    fin_color_list <- cls_renaming_list
  }
  
  # Set color order
  pr440_short_tips_TRY$data$new_class <- factor(pr440_short_tips_TRY$data$new_class, 
                                                levels = unique(dx_rx_types1_short))
  
  # Set edge width based on p-values
  list_of_sizes <- c(1, 1.6, 2.2) * edge_width_multiplier
  
  #print("O")
  
  # Handle special case for gray color
  lll <- length(colors_scale1)
  colors_scale2 <- colors_scale1
  
  if (colors_scale1[lll] == "gray") {
    if (colors_scale1[lll - 1] == "gray") {
      # Remove the last element if both are "gray"
      colors_scale2 <- colors_scale1[1:(lll - 2)]
    } 
  }
  
  #print("P")
  
  colors_scale2 <- c(colors_scale2, 'black')
  
  # DEBUG: Print factor levels and colors before applying scale
  # v53: cat(file=stderr(), "\n[DEBUG] === SCALE_COLOR_MANUAL MAPPING ===\n")
  # v53: cat(file=stderr(), "[DEBUG] Factor levels (dx_rx_types1_short):", paste(unique(dx_rx_types1_short), collapse=", "), "\n")
  # v53: cat(file=stderr(), "[DEBUG] colors_scale2 (what gets applied):", paste(colors_scale2, collapse=", "), "\n")
  # v53: cat(file=stderr(), "[DEBUG] MAPPING: Each factor level gets the color at the same index:\n")
  for (idx in seq_along(unique(dx_rx_types1_short))) {
    level_name <- unique(dx_rx_types1_short)[idx]
    color_assigned <- if(idx <= length(colors_scale2)) colors_scale2[idx] else "NO COLOR"
    # v53: cat(file=stderr(), "[DEBUG]   Level", idx, ":", level_name, "->", color_assigned, "\n")
  }
  # v53: cat(file=stderr(), "[DEBUG] =====================================\n\n")
  
  # Apply colors and sizes to tree
  .prof_section_start <- Sys.time()
  pr440_short_tips_TRY_new <- pr440_short_tips_TRY +
    scale_color_manual(values = colors_scale2) +
    scale_size_manual(values = list_of_sizes)
  cat(file=stderr(), sprintf("[PROF-TREE] scale_color + scale_size: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))

  levels_base <- levels(pr440_short_tips_TRY$data$new_class)
  
  # Store reference to original tree
  t <- pr440_short_tips_TRY_new
  
  # Set up bootstrap display parameters
  square_size_legend <- size_font_legend_box
  text_legend_title_size <- size_font_legend_title
  text_legend_size <- size_font_legend_text
  
  size_70 <- edge_width_multiplier
  size_80 <- size_70 + 1
  size_90 <- size_70 + 2
  
  pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new
  
  # v53: print("Q")
  
  # Add bootstrap visualization if requested
  # v90: Bootstrap triangle layers are stored separately and added AFTER gheatmap
  # to avoid "missing aesthetics x and y" error when gheatmap transforms the data.
  # See bootstrap_triangles_* variables used later after gheatmap section.
  bootstrap_triangles_enabled <- FALSE
  bootstrap_triangles_params <- NULL

  if (show_boot_flag == TRUE) {
    if (boot_values$'format' == 'triangles') {
      # v90: Store parameters for later - DO NOT add layers here
      # These will be added after gheatmap to avoid breaking the layers
      bootstrap_triangles_enabled <- TRUE
      bootstrap_triangles_params <- list(
        man_boot_x_offset = man_boot_x_offset,
        size_90 = size_90 + bootstrap_label_size,
        size_80 = size_80 + bootstrap_label_size,
        size_70 = size_70 + bootstrap_label_size
      )
      # Set to base plot - triangles will be added later
      pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new
    }
    
    if (boot_values$'format' == 'raw') {
      pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new + 
        geom_text(
          aes(label = boot_val_rounded, angle = 90, colour = "black"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "black"
        )
    }
    
    if (boot_values$'format' == 'percentage') {
      pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new + 
        geom_text(
          aes(label = boot_val_rounded, angle = 90, colour = "black"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "black"
        )
    }
    
    if (boot_values$'format' == 'numbered_colored') {
      pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new + 
        geom_text(
          aes(label = boot_val_rounded, angle = 90, colour = "black"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "black"
        )
    }
    
    if (boot_values$'format' %in% c('numbered_color', 'percentage_color')) {
      pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new + 
        geom_text(
          aes(label = boot_val_rounded0, angle = 90, colour = "gray98"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray98"
        ) +
        geom_text(
          aes(label = boot_val_rounded1, angle = 90, color = "gray90"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray90"
        ) +
        geom_text(
          aes(label = boot_val_rounded2, angle = 90, color = "gray80"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray80"
        ) +
        geom_text(
          aes(label = boot_val_rounded3, angle = 90, color = "gray70"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray70"
        ) +
        geom_text(
          aes(label = boot_val_rounded4, angle = 90, color = "gray60"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray60"
        ) +
        geom_text(
          aes(label = boot_val_rounded5, angle = 90, color = "gray50"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray50"
        ) +
        geom_text(
          aes(label = boot_val_rounded6, angle = 90, color = "gray40"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray40"
        ) +
        geom_text(
          aes(label = boot_val_rounded7, angle = 90, color = "gray30"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray30"
        ) +
        geom_text(
          aes(label = boot_val_rounded8, angle = 90, color = "gray20"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray20"
        ) +
        geom_text(
          aes(label = boot_val_rounded9, angle = 90, color = "gray10"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray10"
        ) +
        geom_text(
          aes(label = boot_val_rounded10, angle = 90, color = "gray1"),
          hjust = -0.5, vjust = -0.4, size = bootstrap_label_size,
          show.legend = FALSE, colour = "gray1"
        )
    }
  }
  
  # Set up the base tree visualization
  group_display_title <- "Cell type"
  if (exists('title_flag')) {
    if (title_flag == TRUE) {
      group_display_title <- title_replace
    }
  } 
  
  min_col_x_of_frame_of_prev_heat <- 0
  new_base <- 0
  
  pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new_with_boot
  
  # Add legends and styling
  text_mu <- round(width / 100) - 2
  if (units_out == "cm") {
    text_mu <- text_mu * 10
  }
  mar <- round(width / 50)
  mar1 <- mar + 5

  .prof_section_start <- Sys.time()
  pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot +
    layout_dendrogram() +
    guides(
      colour = guide_legend(
        title = group_display_title, 
        byrow = TRUE,
        override.aes = list(
          size = square_size_legend * text_mu, 
          alpha = 1,
          linetype = "solid", 
          linewidth = 1
        ),
        breaks = levels(pr440_short_tips_TRY_new_with_boot_more1$data$group)
      ),
      size = guide_legend(title = "p value"),
      fill = guide_legend(byrow = TRUE)
    ) + 
    theme(
      legend.text = element_text(
        size = text_legend_size * text_mu,
        margin = margin(t = mar, b = mar, unit = "pt")
      ),
      legend.title = element_text(
        size = text_legend_title_size,
        face = "bold"
      ),
      legend.spacing.y = unit(text_mu / 6, 'cm')
    ) +
    scale_y_reverse() +
    geom_rootedge()
  cat(file=stderr(), sprintf("[PROF-TREE] guides + theme + scale_y_reverse: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))

  # Add node numbers if requested
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 5: NODE NUMBER RENDERING\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â flag_display_nod_number_on_tree is:", flag_display_nod_number_on_tree, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â node_number_font_size is:", node_number_font_size, "\n")
  
  if (flag_display_nod_number_on_tree == TRUE) {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€œ ADDING NODE NUMBERS with size:", node_number_font_size, "\n")
    # v53: debug_cat("================================================\n\n")
    
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
      geom_text(
        aes(label = node, angle = 90, colour = "black"), 
        hjust = -0.5, vjust = -0.4, size = node_number_font_size,
        show.legend = FALSE, colour = "black"
      )
  } else {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€” NODE NUMBERS NOT ENABLED\n")
    # v53: debug_cat("================================================\n\n")
  }
  
  # Highlight manually selected nodes
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 6: HIGHLIGHTING RENDERING\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â highlight_manual_nodes is:", highlight_manual_nodes, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â manual_nodes_to_highlight is:", paste(manual_nodes_to_highlight, collapse=", "), "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â is.na(manual_nodes_to_highlight[1]):", is.na(manual_nodes_to_highlight[1]), "\n")
  
  if (highlight_manual_nodes == TRUE && !is.na(manual_nodes_to_highlight[1])) {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€œ HIGHLIGHTING ENABLED - Creating highlight data\n")
    
    # Create a data frame for highlighted nodes
    highlight_data <- pr440_short_tips_TRY_new_with_boot_more1$data[
      pr440_short_tips_TRY_new_with_boot_more1$data$node %in% manual_nodes_to_highlight, 
    ]
    
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Highlight data rows:", nrow(highlight_data), "\n")
    if (nrow(highlight_data) > 0) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Nodes found in data:", paste(highlight_data$node, collapse=", "), "\n")
    }
    
    if (nrow(highlight_data) > 0) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€œÃ¢Å“â€œ ADDING RED CIRCLES to", nrow(highlight_data), "nodes\n")
      # v53: debug_cat("================================================\n\n")
      
      pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
        geom_point(data = highlight_data, 
                   aes(x = x, y = y),
                   color = "red", 
                   size = 5, 
                   alpha = 0.6,
                   shape = 21,
                   fill = "red",
                   stroke = 2,
                   show.legend = FALSE)
    } else {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€”Ã¢Å“â€” NO MATCHING NODES FOUND IN DATA\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â Available nodes in data (first 20):", 
      #     paste(head(unique(pr440_short_tips_TRY_new_with_boot_more1$data$node), 20), collapse=", "), "\n")
      # v53: debug_cat("================================================\n\n")
    }
  } else {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€” HIGHLIGHTING NOT ENABLED\n")
    if (highlight_manual_nodes != TRUE) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â   Reason: highlight_manual_nodes is FALSE\n")
    }
    if (is.na(manual_nodes_to_highlight[1])) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â   Reason: manual_nodes_to_highlight is NA\n")
    }
    # v53: debug_cat("================================================\n\n")
  }
  
  # Get plot boundaries
  x_range_min <- min(pr440_short_tips_TRY_new_with_boot_more1$data$x)
  x_range_max <- max(pr440_short_tips_TRY_new_with_boot_more1$data$x)
  
  x_off_base <- x_range_min * 0.4 - man_adj_second_legend
  
  # Add tip labels if requested
  if (tip_name_display_flag == TRUE) {
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
      geom_tiplab(size = size_tip_text, angle = 90, hjust = 1, show.legend = FALSE)
  }   
  
  # Set up for legends
  how_mant_rows <- 0
  how_many_boxes <- 0
  
  if (FLAG_BULK_DISPLAY == TRUE) {
    how_mant_rows <- how_mant_rows + how_many_hi
    how_many_boxes <- how_many_boxes + 1
  }
  
  if (show_boot_flag == TRUE) {
    how_mant_rows <- how_mant_rows + 3
    how_many_boxes <- how_many_boxes + 1
  }
  
  p <- pr440_short_tips_TRY_new_with_boot_more1
  b<-0

  # v132: Initialize boudariestt BEFORE highlight code - it's needed for ellipse sizing when heat_flag == TRUE
  # Previously this was initialized at line ~4795 (after heatmap rendering), causing "object 'boudariestt' not found" error
  tt_for_boundaries <- p
  boudariestt <- tryCatch({
    func.find.plot.boundaries(tt_for_boundaries, debug_mode)
  }, error = function(e) {
    debug_cat(paste0("  v32: Error computing boudariestt: ", e$message, "\n"))
    # Return default values if computation fails
    list(xmin = min(p$data$x, na.rm = TRUE), xmax = max(p$data$x, na.rm = TRUE))
  })
  debug_cat(paste0("  v32: boudariestt initialized early: xmin=", boudariestt$xmin, ", xmax=", boudariestt$xmax, "\n"))

  # Handle highlighting if requested
  # v141: Only apply highlight here when NO heatmap, to prevent double highlighting
  # When heat_flag == TRUE, highlight is applied AFTER heatmap processing (around line 6421)
  if (FLAG_BULK_DISPLAY == TRUE && heat_flag == FALSE) {
    x_adj_hi <- 0

    # Calculate tip length for ellipse sizing
    # If trimming is disabled (id_tip_trim_end is NA), use the max tip label length
    tip_length <- if (is.na(id_tip_trim_end)) {
      max(nchar(as.character(tree440$tip.label)), na.rm = TRUE)
    } else {
      id_tip_trim_end
    }

    # v50: Use multiplicative scaling for height and width sliders
    base_a <- tip_length * size_tip_text / 800 + man_adjust_elipse_a
    base_b <- 0.12 + man_adjust_elipse_b
    a <- base_a * adjust_height_ecliplse  # Slider value is now a multiplier (1.0 = no change)
    b <- base_b * adjust_width_eclipse    # Slider value is now a multiplier

    # v139: Pass high_alpha_list for transparency
    .prof_section_start <- Sys.time()
    p <- func_highlight(
      p, how_many_hi, heat_flag, high_color_list, a, b, man_adjust_elipse,
      pr440_short_tips_TRY, boudariestt, debug_mode, high_offset, high_vertical_offset,
      high_alpha_list
    )
    cat(file=stderr(), sprintf("[PROF-TREE] func_highlight #1 (pre-heatmap): %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  }

  if (length(b) == 0) {
    b <- 0.2
  }

  # Add heatmap if requested
  # v99: MANUAL HEATMAP - Replaced gheatmap() with manual geom_tile() approach
  # because gheatmap was corrupting the plot's @mapping property
  debug_cat(paste0("\n=== v99: HEATMAP RENDERING ===\n"))
  debug_cat(paste0("  heat_flag: ", heat_flag, "\n"))
  if (heat_flag == TRUE) {
    debug_cat(paste0("  heat_map_title_list length: ", length(heat_map_title_list), "\n"))
    debug_cat(paste0("  dxdf440_for_heat length: ", length(dxdf440_for_heat), "\n"))
  }

  # v99: Save a backup of p before heatmap for fallback recovery
  tt <- p

  # v132: Refresh boudariestt after highlighting changes (was already initialized at v132 block above)
  # This ensures func.make.second.legend has current plot boundaries
  boudariestt <- tryCatch({
    func.find.plot.boundaries(tt, debug_mode)
  }, error = function(e) {
    debug_cat(paste0("  v32: Error refreshing boudariestt: ", e$message, "\n"))
    # Return default values if computation fails
    list(xmin = min(p$data$x, na.rm = TRUE), xmax = max(p$data$x, na.rm = TRUE))
  })
  debug_cat(paste0("  v32: boudariestt refreshed: xmin=", boudariestt$xmin, ", xmax=", boudariestt$xmax, "\n"))

  # v122: MULTIPLE HEATMAPS IMPLEMENTATION
  # Refactored from v99 to support multiple heatmaps with spacing control
  # Each heatmap can have its own colors, type (discrete/continuous), and parameters
  if (heat_flag == TRUE && length(dxdf440_for_heat) > 0) {
    debug_cat(paste0("\n=== v122: ENTERING MULTIPLE HEATMAPS CODE ===\n"))
    debug_cat(paste0("  Number of heatmaps to render: ", length(dxdf440_for_heat), "\n"))

    # v122: Get tree info once (shared across all heatmaps)
    tip_data <- subset(p$data, isTip == TRUE)
    tree_tips <- tip_data$label
    tree_xmin <- min(p$data$x, na.rm = TRUE)
    tree_xmax <- max(p$data$x, na.rm = TRUE)
    tree_width <- abs(tree_xmax - tree_xmin)

    # v122: Track x position for placing heatmaps side by side
    # current_heatmap_x_end will be updated after each heatmap
    current_heatmap_x_start <- tree_xmax

    # v125: Use global gap setting for spacing between multiple heatmaps
    heatmap_spacing <- heatmap_global_gap * tree_width
    debug_cat(paste0("  v25: heatmap_global_gap=", heatmap_global_gap, ", heatmap_spacing=", heatmap_spacing, "\n"))

    # v122: Loop through all heatmaps
    for (heat_idx in 1:length(dxdf440_for_heat)) {
      debug_cat(paste0("\n=== v122: PROCESSING HEATMAP ", heat_idx, " of ", length(dxdf440_for_heat), " ===\n"))

      # Get this heatmap's data
      heat_data <- dxdf440_for_heat[[heat_idx]]

      # v122: Validate heatmap data
      debug_cat(paste0("  Initial heat_data dimensions: ", nrow(heat_data), " x ", ncol(heat_data), "\n"))

      # v122: Check if heat_data is valid - skip this heatmap if invalid
      if (is.null(heat_data) || !is.data.frame(heat_data) || nrow(heat_data) == 0) {
        debug_cat(paste0("  ERROR: Invalid heatmap data - skipping heatmap ", heat_idx, "\n"))
        next  # v122: Continue to next heatmap instead of stopping all heatmaps
      }

      debug_cat(paste0("  Tree has ", length(tree_tips), " tips\n"))

      # v122: Check current row names
      current_rownames <- rownames(heat_data)
      debug_cat(paste0("  Current rownames: ", paste(head(current_rownames, 5), collapse=", "), "\n"))

      # v122: Verify row names are valid and match tree tips
      if (is.null(current_rownames) || all(current_rownames == as.character(1:nrow(heat_data)))) {
        debug_cat(paste0("  WARNING: Heat data has default numeric row names - skipping heatmap ", heat_idx, "\n"))
        next
      }

      # v122: Check how many row names match tree tips
      matching_tips <- sum(current_rownames %in% tree_tips)
      debug_cat(paste0("  Row names matching tree tips: ", matching_tips, " / ", nrow(heat_data), "\n"))

      if (matching_tips == 0) {
        debug_cat(paste0("  ERROR: No row names match tree tips - skipping heatmap ", heat_idx, "\n"))
        next
      } else if (matching_tips < nrow(heat_data)) {
        debug_cat(paste0("  WARNING: Only ", matching_tips, " row names match - filtering data\n"))
        heat_data <- heat_data[current_rownames %in% tree_tips, , drop = FALSE]
        debug_cat(paste0("  After filtering: ", nrow(heat_data), " rows\n"))
      }

      debug_cat(paste0("  heat_data dimensions: ", nrow(heat_data), " x ", ncol(heat_data), "\n"))
      debug_cat(paste0("  heat_data rownames sample: ", paste(head(rownames(heat_data), 5), collapse=", "), "\n"))
      debug_cat(paste0("  heat_data columns: ", paste(colnames(heat_data), collapse=", "), "\n"))

      # v122: Get heatmap parameters for THIS heatmap (use heat_idx, not 1)
      heat_param <- if (heat_idx <= length(heat_display_params_list)) heat_display_params_list[[heat_idx]] else list()
      is_discrete <- ifelse(!is.null(heat_param) && !is.na(heat_param['is_discrete']),
                            heat_param['is_discrete'] == TRUE, FALSE)
      debug_cat(paste0("  is_discrete: ", is_discrete, "\n"))

      # v122: Calculate heatmap positioning - use current_heatmap_x_start from previous heatmap
      # For first heatmap, this is tree_xmax; for subsequent, it's the end of the previous heatmap
      per_heatmap_distance <- if (!is.null(heat_param[['distance']])) heat_param[['distance']] else heatmap_tree_distance
      heatmap_offset <- tree_width * per_heatmap_distance
      debug_cat(paste0("  per_heatmap_distance: ", per_heatmap_distance, "\n"))
      debug_cat(paste0("  current_heatmap_x_start: ", current_heatmap_x_start, "\n"))
      # v113: Fixed row height and column width to be INDEPENDENT
      # In coord_flip context with ggtree:
      # - Data y (tip indices 1,2,3...) becomes visual x after flip
      # - Data x (column positions) becomes visual y after flip
      # So:
      # - geom_tile 'height' (data y extent) -> controls visual row height
      # - geom_tile 'width' (data x extent) -> controls visual column width

      # Column Width slider controls how wide each heatmap column appears
      column_width_value <- if (!is.null(heat_param[['height']])) heat_param[['height']] else 0.8

      # Row Height slider controls how tall each row (tip strip) appears
      # Values: 0.5 = half height (gaps between rows), 1.0 = full (touching), 2.0 = overlapping
      row_height_value <- if (!is.null(heat_param[['row_height']])) heat_param[['row_height']] else 1.0

      # tile_width: controls column spacing (data x units)
      # Base is 3% of tree width, scaled by Column Width slider
      tile_width <- tree_width * 0.03 * column_width_value

      # v116: Calculate actual tip spacing for row height
      # Tip y-positions are NOT always at 1,2,3... - they depend on tree topology
      # We calculate the MINIMUM spacing to prevent overlap at row_height_value=1.0
      tip_y_positions <- sort(unique(tip_data$y))
      if (length(tip_y_positions) > 1) {
        tip_spacings <- diff(tip_y_positions)
        base_tip_spacing <- min(tip_spacings)  # v116: Use minimum to prevent default overlap
        mean_spacing <- mean(tip_spacings)
        median_spacing <- median(tip_spacings)
        debug_cat(paste0("  v16: tip spacing stats: min=", base_tip_spacing,
                                   ", median=", median_spacing,
                                   ", mean=", mean_spacing, "\n"))
      } else {
        base_tip_spacing <- 1.0  # Fallback for single tip
      }

      # tile_height: controls row height (data y units)
      # row_height_value=1.0 means tiles touch at minimum spacing (no overlap)
      # <1 means gaps, >1 means overlap
      tile_height <- base_tip_spacing * row_height_value

      debug_cat(paste0("  v16: tip_y range: [", min(tip_y_positions), ", ", max(tip_y_positions), "]\n"))
      debug_cat(paste0("  v16: n_tips: ", length(tip_y_positions), "\n"))
      debug_cat(paste0("  v16: base_tip_spacing (min): ", base_tip_spacing, "\n"))

      debug_cat(paste0("  column_width_value: ", column_width_value, "\n"))
      debug_cat(paste0("  row_height_value: ", row_height_value, "\n"))
      debug_cat(paste0("  tile_width (column spacing): ", tile_width, "\n"))
      debug_cat(paste0("  tile_height (row height): ", tile_height, "\n"))

      debug_cat(paste0("  tree_xmin: ", tree_xmin, "\n"))
      debug_cat(paste0("  tree_xmax: ", tree_xmax, "\n"))
      debug_cat(paste0("  tree_width: ", tree_width, "\n"))
      debug_cat(paste0("  heatmap_offset: ", heatmap_offset, "\n"))
      debug_cat(paste0("  tile_width: ", tile_width, "\n"))

      # v99: Build heatmap data frame for geom_tile
      # We need: x (column position), y (tip position), fill (value)
      debug_cat(paste0("\n=== v99: BUILDING HEATMAP TILE DATA ===\n"))

      tile_data_list <- list()
      for (col_idx in 1:ncol(heat_data)) {
        col_name <- colnames(heat_data)[col_idx]

        for (row_idx in 1:nrow(heat_data)) {
          tip_label <- rownames(heat_data)[row_idx]
          value <- heat_data[row_idx, col_idx]

          # Find y position from tree tip data
          tip_row <- tip_data[tip_data$label == tip_label, ]
          if (nrow(tip_row) > 0) {
            y_pos <- tip_row$y[1]
            # v122: Use current_heatmap_x_start (which tracks position after previous heatmaps)
            x_pos <- current_heatmap_x_start + heatmap_offset + (col_idx - 0.5) * tile_width

            # v111: For continuous data, keep value as-is (numeric)
            # For discrete data, convert to character
            if (is_discrete) {
              tile_value <- as.character(value)
            } else {
              # Try to convert to numeric, keep as character if fails
              tile_value <- suppressWarnings(as.numeric(value))
              if (is.na(tile_value) && !is.na(value)) {
                # If conversion to numeric failed but original wasn't NA,
                # this is likely a non-numeric value - treat as NA for continuous
                tile_value <- NA_real_
              }
            }

            tile_data_list[[length(tile_data_list) + 1]] <- data.frame(
              x = x_pos,
              y = y_pos,
              value = if (is_discrete) tile_value else tile_value,
              column = col_name,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(tile_data_list) > 0) {
        tile_df <- do.call(rbind, tile_data_list)
        debug_cat(paste0("  Created tile_df with ", nrow(tile_df), " tiles\n"))
        debug_cat(paste0("  x range: [", min(tile_df$x), ", ", max(tile_df$x), "]\n"))
        debug_cat(paste0("  y range: [", min(tile_df$y), ", ", max(tile_df$y), "]\n"))
        debug_cat(paste0("  Unique values: ", paste(unique(tile_df$value), collapse=", "), "\n"))

        # v99: Add heatmap tiles using geom_tile
        debug_cat(paste0("\n=== v99: ADDING GEOM_TILE LAYER ===\n"))

        p <- tryCatch({
          # v112: Add tile layer with explicit aesthetics
          # In coord_flip context:
          # - geom_tile 'width' (data x) -> visual y (column spacing/width)
          # - geom_tile 'height' (data y) -> visual x (row height/span)
          # So:
          # - width = tile_width (fixed column spacing)
          # - height = tile_height (Row Height slider * base)

          # v113: Get grid settings with debug output
          show_grid <- if (!is.null(heat_param[['show_grid']])) heat_param[['show_grid']] else FALSE
          grid_color <- if (!is.null(heat_param[['grid_color']])) heat_param[['grid_color']] else "#000000"
          grid_size <- if (!is.null(heat_param[['grid_size']])) as.numeric(heat_param[['grid_size']]) else 0.5

          debug_cat(paste0("  v13: show_grid=", show_grid, " (class: ", class(show_grid), ")\n"))
          debug_cat(paste0("  v13: grid_color=", grid_color, ", grid_size=", grid_size, "\n"))

          # v113: Ensure show_grid is properly evaluated as boolean
          show_grid_bool <- isTRUE(show_grid) || identical(show_grid, TRUE) || identical(show_grid, "yes") || identical(show_grid, "TRUE")
          debug_cat(paste0("  v13: show_grid_bool=", show_grid_bool, "\n"))

          # v123: For heatmaps after the first, add new_scale_fill() BEFORE adding geom_tile
          # This is critical - the scale must be reset before the new layer that uses fill
          if (heat_idx > 1) {
            debug_cat(paste0("  v23: Adding new_scale_fill() BEFORE geom_tile for heatmap ", heat_idx, "\n"))
            p <- p + ggnewscale::new_scale_fill()
          }

          # S2.8: Check for RData detailed display mode
          cat(file=stderr(), paste0("[S2.8-DEBUG] heat_param data_source: ",
              ifelse(!is.null(heat_param[['data_source']]), heat_param[['data_source']], "NULL"), "\n"))
          cat(file=stderr(), paste0("[S2.8-DEBUG] heat_param cnv_display_mode: ",
              ifelse(!is.null(heat_param[['cnv_display_mode']]), heat_param[['cnv_display_mode']], "NULL"), "\n"))

          is_rdata_detailed <- !is.null(heat_param[['data_source']]) &&
                               heat_param[['data_source']] == "rdata" &&
                               !is.null(heat_param[['cnv_display_mode']]) &&
                               heat_param[['cnv_display_mode']] == "detailed"
          cat(file=stderr(), paste0("[S2.8-DEBUG] is_rdata_detailed: ", is_rdata_detailed, "\n"))

          if (is_rdata_detailed) {
            debug_cat(paste0("\n=== S2.8: DETAILED MODE (geom_tile with pre-computed colors like pheatmap) ===\n"))
            cat(file=stderr(), "[HEATMAP-RENDER] Using DETAILED mode with geom_tile + pre-computed colors (pheatmap-style)\n")

            # S2.9-PERF: Check heatmap cache before expensive color computation
            heatmap_cache_key <- as.character(heat_idx)
            current_heatmap_hash <- func.create.heatmap_cache_hash(heat_param, heat_data, tip_data$label)

            use_cached_colors <- FALSE
            cached_tile_df <- NULL
            cached_value_range <- NULL

            if (!is.null(heatmap_cache) && heatmap_cache_key %in% names(heatmap_cache)) {
              cached_entry <- heatmap_cache[[heatmap_cache_key]]
              if (!is.null(cached_entry$hash) && cached_entry$hash == current_heatmap_hash) {
                cat(file=stderr(), paste0("[S2.9-CACHE] Heatmap ", heat_idx, " CACHE HIT - reusing pre-computed colors\n"))
                use_cached_colors <- TRUE
                cached_tile_df <- cached_entry$tile_df
                cached_value_range <- cached_entry$value_range
              } else {
                cat(file=stderr(), paste0("[S2.9-CACHE] Heatmap ", heat_idx, " cache MISS (hash mismatch)\n"))
              }
            } else {
              cat(file=stderr(), paste0("[S2.9-CACHE] Heatmap ", heat_idx, " cache MISS (no cached entry)\n"))
            }

            # For detailed mode, we use geom_tile with pre-computed colors like pheatmap
            # This ensures each tile is at the exact tip y-coordinate while having pheatmap-style coloring

            # Get height scale for detailed mode (allows making the heatmap shorter)
            # Note: Due to 90-degree rotation, "height" on screen = tile_width in code
            height_scale <- if (!is.null(heat_param[['cnv_height_scale']]) && !is.na(heat_param[['cnv_height_scale']])) {
              as.numeric(heat_param[['cnv_height_scale']])
            } else {
              1.0
            }
            cat(file=stderr(), paste0("[S2.8-DETAILED] Height scale: ", height_scale, "\n"))

            # Apply height scale to tile_width for detailed mode (controls heatmap "height" on screen)
            tile_width_detailed <- tile_width * height_scale
            cat(file=stderr(), paste0("[S2.8-DETAILED] Adjusted tile_width: ", tile_width_detailed, " (original: ", tile_width, ")\n"))

            # Get color parameters
            low_color <- if (!is.null(heat_param[['low']]) && !is.na(heat_param[['low']])) heat_param[['low']] else "#FF0000"
            mid_color <- if (!is.null(heat_param[['mid']]) && !is.na(heat_param[['mid']])) heat_param[['mid']] else "#FFFFFF"
            high_color <- if (!is.null(heat_param[['high']]) && !is.na(heat_param[['high']])) heat_param[['high']] else "#0000FF"
            midpoint <- if (!is.null(heat_param[['midpoint']]) && !is.na(heat_param[['midpoint']])) as.numeric(heat_param[['midpoint']]) else 2
            na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "grey90"

            debug_cat(paste0("  S2.8 Detailed colors: low=", low_color, ", mid=", mid_color, ", high=", high_color, "\n"))
            debug_cat(paste0("  S2.8 Detailed midpoint: ", midpoint, "\n"))

            # S2.9-PERF: Use cached colors if available, otherwise compute
            if (use_cached_colors && !is.null(cached_tile_df)) {
              # Use cached tile_df with pre-computed fill_color
              tile_df <- cached_tile_df
              value_min <- cached_value_range$min
              value_max <- cached_value_range$max
              cat(file=stderr(), paste0("[S2.9-CACHE] Using cached tile_df with ", nrow(tile_df), " tiles\n"))

              # Still need to create palette for legend
              colors_below <- colorRampPalette(c(low_color, mid_color))(999)
              colors_above <- colorRampPalette(c(mid_color, high_color))(1000)[-1]
              detailed_palette <- c(colors_below, colors_above)
            } else {
              # Create fine-grained color palette like pheatmap (1998 colors)
              colors_below <- colorRampPalette(c(low_color, mid_color))(999)
              colors_above <- colorRampPalette(c(mid_color, high_color))(1000)[-1]
              detailed_palette <- c(colors_below, colors_above)

              # S2.8: Enhanced debug logging for matrix dimensions
              cat(file=stderr(), paste0("[S2.8-MATRIX] heat_data dimensions: ", nrow(heat_data), " rows x ", ncol(heat_data), " cols\n"))
              cat(file=stderr(), paste0("[S2.8-MATRIX] Number of tree tips in tip_data: ", nrow(tip_data), "\n"))
              cat(file=stderr(), paste0("[S2.8-MATRIX] tile_df rows: ", nrow(tile_df), "\n"))

              # Get data range for color mapping
              data_values <- tile_df$value[!is.na(tile_df$value)]
              if (length(data_values) > 0) {
                data_min <- min(data_values, na.rm = TRUE)
                data_max <- max(data_values, na.rm = TRUE)

                range_below <- midpoint - min(data_min, 0)
                range_above <- max(data_max, midpoint + 2) - midpoint
                max_range <- max(range_below, range_above, 2)

                value_min <- midpoint - max_range
                value_max <- midpoint + max_range

                debug_cat(paste0("  S2.8 Data range: [", round(data_min, 2), ", ", round(data_max, 2), "]\n"))
                debug_cat(paste0("  S2.8 Color range: [", round(value_min, 2), ", ", round(value_max, 2), "]\n"))

                # Pre-compute colors for each tile (like pheatmap does)
                # Normalize values to 0-1 range, then map to palette indices
                normalized_values <- (tile_df$value - value_min) / (value_max - value_min)
                normalized_values[normalized_values < 0] <- 0
                normalized_values[normalized_values > 1] <- 1

                # Convert to color indices (1 to 1998)
                color_indices <- round(normalized_values * (length(detailed_palette) - 1)) + 1

                # Assign pre-computed colors
                tile_df$fill_color <- detailed_palette[color_indices]
                tile_df$fill_color[is.na(tile_df$value)] <- na_color

                cat(file=stderr(), paste0("[S2.8-DETAILED] Pre-computed colors for ", nrow(tile_df), " tiles\n"))
                cat(file=stderr(), paste0("[S2.8-DETAILED] Unique colors: ", length(unique(tile_df$fill_color)), "\n"))

                # S2.9-PERF: Store in cache for future use
                if (!exists("updated_heatmap_cache")) {
                  updated_heatmap_cache <- heatmap_cache
                }
                updated_heatmap_cache[[heatmap_cache_key]] <- list(
                  hash = current_heatmap_hash,
                  tile_df = tile_df,
                  value_range = list(min = value_min, max = value_max)
                )
                cat(file=stderr(), paste0("[S2.9-CACHE] Cached heatmap ", heat_idx, " for future use\n"))
              }
            }

            # Continue with rendering (both cached and fresh paths converge here)
            data_values <- tile_df$value[!is.na(tile_df$value)]
            if (length(data_values) > 0) {

              # Use geom_tile with pre-computed colors and scale_fill_identity
              # This positions each tile at exact tip y-coordinates
              # Use tile_width_detailed (scaled by height_scale) to control heatmap "height" on screen
              p_with_tiles <- p + geom_tile(
                data = tile_df,
                aes(x = x, y = y, fill = fill_color),
                width = tile_width_detailed,
                height = tile_height,
                inherit.aes = FALSE
              ) + scale_fill_identity(guide = "none")  # No legend for identity scale

              # Add a proper gradient legend using invisible points
              # This creates a continuous colorbar legend that matches the pre-computed colors
              legend_df <- data.frame(
                x = rep(min(tile_df$x), 5),
                y = rep(min(tile_df$y), 5),
                value = seq(value_min, value_max, length.out = 5)
              )

              p_with_tiles <- p_with_tiles +
                ggnewscale::new_scale_fill() +
                geom_point(data = legend_df, aes(x = x, y = y, fill = value),
                           alpha = 0, shape = 22, inherit.aes = FALSE) +
                scale_fill_gradientn(
                  colors = detailed_palette,
                  limits = c(value_min, value_max),
                  name = heatmap_title,
                  na.value = na_color
                )

              debug_cat(paste0("  S2.8 Added geom_tile with ", length(detailed_palette), " pre-computed colors\n"))
            } else {
              # Fallback to basic geom_tile if no data
              debug_cat("  S2.8 WARNING: No valid data values, falling back to basic mode\n")
              p_with_tiles <- p + geom_tile(
                data = tile_df,
                aes(x = x, y = y, fill = value),
                width = tile_width_detailed,
                height = tile_height,
                inherit.aes = FALSE
              ) + scale_fill_gradient2(
                low = low_color, mid = mid_color, high = high_color,
                midpoint = midpoint,
                name = heatmap_title,
                na.value = na_color
              )
            }

          } else if (show_grid_bool) {
            # v115: Draw tiles WITHOUT borders first (to avoid overlapping border issues)
            p_with_tiles <- p + geom_tile(
              data = tile_df,
              aes(x = x, y = y, fill = value),
              width = tile_width,     # v112: Column spacing (fixed)
              height = tile_height,   # v112: Row height (slider-controlled)
              inherit.aes = FALSE
            )

            # v115: Draw explicit grid lines to ensure all borders appear
            # Calculate grid line positions
            n_cols <- ncol(heat_data)
            col_x_positions <- unique(tile_df$x)
            col_x_positions <- sort(col_x_positions)

            # Calculate y-range for vertical lines (spanning all tips)
            y_min <- min(tile_df$y) - tile_height / 2
            y_max <- max(tile_df$y) + tile_height / 2

            # Build vertical grid lines (n_cols + 1 lines: left edge, between columns, right edge)
            v_lines_x <- c()
            for (i in seq_along(col_x_positions)) {
              # Left edge of this column
              v_lines_x <- c(v_lines_x, col_x_positions[i] - tile_width / 2)
            }
            # Right edge of last column
            if (length(col_x_positions) > 0) {
              v_lines_x <- c(v_lines_x, col_x_positions[length(col_x_positions)] + tile_width / 2)
            }
            v_lines_x <- unique(v_lines_x)

            # Create vertical line data
            v_lines_df <- data.frame(
              x = rep(v_lines_x, each = 1),
              xend = rep(v_lines_x, each = 1),
              y = y_min,
              yend = y_max
            )

            # Build horizontal grid lines (for each row boundary)
            tip_y_values <- sort(unique(tile_df$y))
            h_lines_y <- c()
            for (i in seq_along(tip_y_values)) {
              # Top and bottom edge of each row
              h_lines_y <- c(h_lines_y, tip_y_values[i] - tile_height / 2)
              h_lines_y <- c(h_lines_y, tip_y_values[i] + tile_height / 2)
            }
            h_lines_y <- unique(h_lines_y)

            # Calculate x-range for horizontal lines
            x_min <- min(tile_df$x) - tile_width / 2
            x_max <- max(tile_df$x) + tile_width / 2

            # Create horizontal line data
            h_lines_df <- data.frame(
              x = x_min,
              xend = x_max,
              y = rep(h_lines_y, each = 1),
              yend = rep(h_lines_y, each = 1)
            )

            # Add grid lines as separate geom_segment layers
            p_with_tiles <- p_with_tiles +
              geom_segment(data = v_lines_df, aes(x = x, xend = xend, y = y, yend = yend),
                           color = grid_color, linewidth = grid_size, inherit.aes = FALSE) +
              geom_segment(data = h_lines_df, aes(x = x, xend = xend, y = y, yend = yend),
                           color = grid_color, linewidth = grid_size, inherit.aes = FALSE)

            debug_cat(paste0("  v15: Added explicit grid lines: ", nrow(v_lines_df), " vertical, ", nrow(h_lines_df), " horizontal\n"))
          } else {
            p_with_tiles <- p + geom_tile(
              data = tile_df,
              aes(x = x, y = y, fill = value),
              width = tile_width,     # v112: Column spacing (fixed)
              height = tile_height,   # v112: Row height (slider-controlled)
              inherit.aes = FALSE
            )
          }

          # S1.62dev: Add horizontal row lines (separate from grid - only horizontal lines)
          show_row_lines_bool <- !is.null(heat_param[['show_row_lines']]) &&
                                  isTRUE(heat_param[['show_row_lines']])
          if (show_row_lines_bool) {
            row_line_color <- if (!is.null(heat_param[['row_line_color']])) heat_param[['row_line_color']] else "#000000"
            row_line_size <- if (!is.null(heat_param[['row_line_size']])) as.numeric(heat_param[['row_line_size']]) else 0.5

            debug_cat(paste0("  S1.62dev: Adding horizontal row lines (color: ", row_line_color, ", size: ", row_line_size, ")\n"))

            # Build horizontal row lines (above and below each row)
            tip_y_values <- sort(unique(tile_df$y))
            row_lines_y <- c()
            for (i in seq_along(tip_y_values)) {
              # Top and bottom edge of each row
              row_lines_y <- c(row_lines_y, tip_y_values[i] - tile_height / 2)
              row_lines_y <- c(row_lines_y, tip_y_values[i] + tile_height / 2)
            }
            row_lines_y <- unique(row_lines_y)

            # Calculate x-range for horizontal lines
            row_x_min <- min(tile_df$x) - tile_width / 2
            row_x_max <- max(tile_df$x) + tile_width / 2

            # Create horizontal line data (NO vertical lines)
            row_lines_df <- data.frame(
              x = row_x_min,
              xend = row_x_max,
              y = rep(row_lines_y, each = 1),
              yend = rep(row_lines_y, each = 1)
            )

            # Add row lines as geom_segment layer
            p_with_tiles <- p_with_tiles +
              geom_segment(data = row_lines_df, aes(x = x, xend = xend, y = y, yend = yend),
                           color = row_line_color, linewidth = row_line_size, inherit.aes = FALSE)

            debug_cat(paste0("  S1.62dev: Added ", nrow(row_lines_df), " horizontal row lines\n"))
          }

          # S1.62dev: Add vertical column lines (separate from grid - only vertical lines)
          cat(file=stderr(), paste0("[DEBUG-COLLINES-RENDER] heat_param[['show_col_lines']] = ", ifelse(is.null(heat_param[['show_col_lines']]), "NULL", heat_param[['show_col_lines']]), "\n"))
          cat(file=stderr(), paste0("[DEBUG-COLLINES-RENDER] 'show_col_lines' in names: ", 'show_col_lines' %in% names(heat_param), "\n"))
          show_col_lines_bool <- !is.null(heat_param[['show_col_lines']]) &&
                                  isTRUE(heat_param[['show_col_lines']])
          cat(file=stderr(), paste0("[DEBUG-COLLINES-RENDER] show_col_lines_bool = ", show_col_lines_bool, "\n"))
          if (show_col_lines_bool) {
            col_line_color <- if (!is.null(heat_param[['col_line_color']])) heat_param[['col_line_color']] else "#000000"
            col_line_size <- if (!is.null(heat_param[['col_line_size']])) as.numeric(heat_param[['col_line_size']]) else 0.5

            debug_cat(paste0("  S1.62dev: Adding vertical column lines (color: ", col_line_color, ", size: ", col_line_size, ")\n"))

            # Build vertical column lines (left and right of each column)
            col_x_values <- sort(unique(tile_df$x))
            col_lines_x <- c()
            for (j in seq_along(col_x_values)) {
              # Left and right edge of each column
              col_lines_x <- c(col_lines_x, col_x_values[j] - tile_width / 2)
              col_lines_x <- c(col_lines_x, col_x_values[j] + tile_width / 2)
            }
            col_lines_x <- unique(col_lines_x)

            # Calculate y-range for vertical lines
            col_y_min <- min(tile_df$y) - tile_height / 2
            col_y_max <- max(tile_df$y) + tile_height / 2

            # Create vertical line data (NO horizontal lines)
            col_lines_df <- data.frame(
              x = rep(col_lines_x, each = 1),
              xend = rep(col_lines_x, each = 1),
              y = col_y_min,
              yend = col_y_max
            )

            # Add column lines as geom_segment layer
            p_with_tiles <- p_with_tiles +
              geom_segment(data = col_lines_df, aes(x = x, xend = xend, y = y, yend = yend),
                           color = col_line_color, linewidth = col_line_size, inherit.aes = FALSE)

            debug_cat(paste0("  S1.62dev: Added ", nrow(col_lines_df), " vertical column lines\n"))
          }

          # S1.62dev: Add vertical text labels below heatmap
          show_vertical_text_bool <- !is.null(heat_param[['show_vertical_text']]) &&
                                      isTRUE(heat_param[['show_vertical_text']])
          if (show_vertical_text_bool) {
            vertical_text_column <- if (!is.null(heat_param[['vertical_text_column']])) heat_param[['vertical_text_column']] else ""
            vertical_text_size <- if (!is.null(heat_param[['vertical_text_size']])) as.numeric(heat_param[['vertical_text_size']]) else 3
            vertical_text_offset <- if (!is.null(heat_param[['vertical_text_offset']])) as.numeric(heat_param[['vertical_text_offset']]) else 0.5
            vertical_text_color <- if (!is.null(heat_param[['vertical_text_color']])) heat_param[['vertical_text_color']] else "#000000"

            debug_cat(paste0("  S1.62dev: Adding vertical text labels (column: '", vertical_text_column,
                             "', size: ", vertical_text_size, ", offset: ", vertical_text_offset, ")\n"))

            # Get column x positions and names
            col_x_positions <- sort(unique(tile_df$x))
            col_names <- colnames(heat_data)

            # Y position below the heatmap
            y_bottom <- min(tile_df$y) - tile_height / 2 - vertical_text_offset

            # Get labels - either from CSV column or use heatmap column names
            vertical_labels <- col_names
            if (!is.null(vertical_text_column) && vertical_text_column != "" && exists("data_table") && !is.null(data_table)) {
              # Try to get labels from CSV column
              # The vertical_text_column contains labels, and we need to map heatmap column names to those labels
              # If CSV has a column with the same values as heatmap column names, use corresponding vertical_text_column values
              if (vertical_text_column %in% names(data_table)) {
                debug_cat(paste0("  S1.62dev: Found column '", vertical_text_column, "' in CSV\n"))
                # For RData CNV heatmaps, column names might be genomic positions
                # Try to find a mapping in the CSV (column names as row values)
                csv_col_values <- data_table[[vertical_text_column]]
                debug_cat(paste0("  S1.62dev: CSV column has ", length(csv_col_values), " values\n"))

                # If the number of CSV rows matches the number of heatmap columns, use them directly
                if (length(csv_col_values) >= length(col_names)) {
                  vertical_labels <- as.character(csv_col_values[1:length(col_names)])
                  debug_cat(paste0("  S1.62dev: Using first ", length(col_names), " values from CSV column as labels\n"))
                }
              } else {
                debug_cat(paste0("  S1.62dev: Column '", vertical_text_column, "' not found in CSV, using column names\n"))
              }
            }

            # Create text label data frame
            if (length(col_x_positions) == length(vertical_labels)) {
              text_df <- data.frame(
                x = col_x_positions,
                y = y_bottom,
                label = vertical_labels,
                stringsAsFactors = FALSE
              )

              # Add vertical text labels using geom_text with 90-degree rotation
              p_with_tiles <- p_with_tiles +
                geom_text(data = text_df, aes(x = x, y = y, label = label),
                          angle = 90, hjust = 1, vjust = 0.5,
                          size = vertical_text_size, color = vertical_text_color,
                          inherit.aes = FALSE)

              debug_cat(paste0("  S1.62dev: Added ", nrow(text_df), " vertical text labels\n"))
            } else {
              debug_cat(paste0("  S1.62dev: Warning - column count mismatch: ", length(col_x_positions),
                               " positions vs ", length(vertical_labels), " labels\n"))
            }
          }

          debug_cat(paste0("  geom_tile added successfully\n"))

          # v100: Add color scale using user-selected colors
          na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "grey90"
          # v122: Use heat_idx to get the correct heatmap title
          heatmap_title <- if (heat_idx <= length(heat_map_title_list)) heat_map_title_list[[heat_idx]] else paste0("Heatmap ", heat_idx)

          # v123: new_scale_fill() is now added BEFORE geom_tile (see line ~4993)

          if (is_discrete) {
            debug_cat(paste0("  Adding discrete color scale\n"))

            # v100: Check for custom colors from user
            man_define_colors <- !is.null(heat_param['man_define_colors']) &&
                                 !is.na(heat_param['man_define_colors']) &&
                                 heat_param['man_define_colors'] == TRUE
            custom_colors_raw <- heat_param[['color_scale_option']]

            # v104: Convert from list to named vector if needed (list preserves names through YAML)
            if (is.list(custom_colors_raw) && !is.null(names(custom_colors_raw))) {
              # Convert named list to named character vector
              custom_colors <- unlist(custom_colors_raw)
              debug_cat(paste0("  v04: Converted custom_colors from named list to named vector\n"))
            } else {
              custom_colors <- custom_colors_raw
            }

            debug_cat(paste0("  man_define_colors: ", man_define_colors, "\n"))
            debug_cat(paste0("  custom_colors: ", paste(custom_colors, collapse=", "), "\n"))
            debug_cat(paste0("  na_color: ", na_color, "\n"))

            if (man_define_colors && !is.null(custom_colors) && length(custom_colors) > 0) {
              # v101: Use custom colors provided by user
              debug_cat(paste0("  Using ", length(custom_colors), " custom colors from user\n"))

              # v104: custom_colors should now be a named vector where names are the value labels
              # Debug: show what we received
              debug_cat(paste0("  custom_colors names: ", paste(names(custom_colors), collapse=", "), "\n"))
              debug_cat(paste0("  custom_colors values: ", paste(custom_colors, collapse=", "), "\n"))

              # Get unique values from data (excluding NA) to match with colors
              unique_vals <- unique(tile_df$value)
              unique_vals <- unique_vals[!is.na(unique_vals)]
              debug_cat(paste0("  Unique values in data: ", paste(unique_vals, collapse=", "), "\n"))

              # v101: Fix color mapping - custom_colors is already a named vector
              # Match colors by value name, not by position
              if (!is.null(names(custom_colors)) && length(names(custom_colors)) > 0) {
                # Use the named vector directly - lookup by value name
                color_vec <- custom_colors[unique_vals]
                # For any values not in custom_colors, use a fallback color
                missing_vals <- unique_vals[is.na(color_vec)]
                if (length(missing_vals) > 0) {
                  debug_cat(paste0("  WARNING: Missing colors for: ", paste(missing_vals, collapse=", "), "\n"))
                  # Use grey for missing values
                  color_vec[is.na(color_vec)] <- "grey50"
                }
              } else {
                # Fallback: custom_colors is unnamed, use positional assignment (sorted order)
                debug_cat(paste0("  custom_colors is unnamed, using positional assignment\n"))
                sorted_vals <- sort(unique_vals)
                if (length(custom_colors) >= length(sorted_vals)) {
                  color_vec <- setNames(custom_colors[1:length(sorted_vals)], sorted_vals)
                } else {
                  color_vec <- setNames(rep(custom_colors, length.out = length(sorted_vals)), sorted_vals)
                }
              }
              debug_cat(paste0("  Color mapping: ", paste(names(color_vec), "=", color_vec, collapse=", "), "\n"))

              p_with_tiles <- p_with_tiles + scale_fill_manual(
                values = color_vec,
                name = heatmap_title,
                na.value = na_color
              )
            } else if (!is.null(custom_colors) && length(custom_colors) == 1 &&
                       custom_colors %in% rownames(RColorBrewer::brewer.pal.info)) {
              # v100: Use RColorBrewer palette
              debug_cat(paste0("  Using RColorBrewer palette: ", custom_colors, "\n"))
              p_with_tiles <- p_with_tiles + scale_fill_brewer(
                palette = custom_colors,
                name = heatmap_title,
                na.value = na_color
              )
            } else {
              # v100: Fallback to viridis
              debug_cat(paste0("  Using default viridis palette\n"))
              p_with_tiles <- p_with_tiles + scale_fill_viridis_d(
                name = heatmap_title,
                na.value = na_color
              )
            }
          } else {
            debug_cat(paste0("  Adding continuous color scale\n"))

            # S1.62dev: Debug - show what heat_param contains for colors
            cat(file=stderr(), paste0("[DEBUG-COLOR-RENDER] heat_param['low'] = ", ifelse(is.null(heat_param['low']), "NULL", paste0("'", heat_param['low'], "'")), "\n"))
            cat(file=stderr(), paste0("[DEBUG-COLOR-RENDER] heat_param[['low']] = ", ifelse(is.null(heat_param[['low']]), "NULL", paste0("'", heat_param[['low']], "'")), "\n"))
            cat(file=stderr(), paste0("[DEBUG-COLOR-RENDER] heat_param names: ", paste(names(heat_param), collapse=", "), "\n"))
            cat(file=stderr(), paste0("[DEBUG-COLOR-RENDER] 'low' in names: ", 'low' %in% names(heat_param), "\n"))

            # v100: Get continuous scale colors from parameters
            # S1.62dev: Use double brackets for proper value extraction
            low_color <- if (!is.null(heat_param[['low']]) && !is.na(heat_param[['low']])) heat_param[['low']] else "beige"
            mid_color <- if (!is.null(heat_param[['mid']]) && !is.na(heat_param[['mid']])) heat_param[['mid']] else "seashell2"
            high_color <- if (!is.null(heat_param[['high']]) && !is.na(heat_param[['high']])) heat_param[['high']] else "firebrick4"
            midpoint <- if (!is.null(heat_param[['midpoint']]) && !is.na(heat_param[['midpoint']])) as.numeric(heat_param[['midpoint']]) else 0.02
            limits <- heat_param[['limits']]

            # v113: Debug output for continuous scale colors including NA color
            cat(file=stderr(), paste0("[DEBUG-COLOR-RENDER] Final colors: low=", low_color, ", mid=", mid_color, ", high=", high_color, ", midpoint=", midpoint, "\n"))
            debug_cat(paste0("  Colors: low=", low_color, ", mid=", mid_color, ", high=", high_color, "\n"))
            debug_cat(paste0("  Midpoint: ", midpoint, "\n"))
            debug_cat(paste0("  v13: na_color for continuous: ", na_color, "\n"))
            debug_cat(paste0("  v13: heat_param na_color value: ", ifelse(is.null(heat_param[['na_color']]), "NULL", heat_param[['na_color']]), "\n"))

            # v111: Values should already be numeric from tile building above
            # No need to convert here as the geom already has the numeric data

            if (!is.null(limits) && !all(is.na(limits))) {
              p_with_tiles <- p_with_tiles + scale_fill_gradient2(
                low = low_color, mid = mid_color, high = high_color,
                midpoint = midpoint,
                name = heatmap_title,
                na.value = na_color,
                limits = limits
              )
            } else {
              p_with_tiles <- p_with_tiles + scale_fill_gradient2(
                low = low_color, mid = mid_color, high = high_color,
                midpoint = midpoint,
                name = heatmap_title,
                na.value = na_color
              )
            }
          }

          debug_cat(paste0("  Color scale added successfully\n"))

          # v105: Add row labels if enabled
          show_row_labels <- if (!is.null(heat_param[['show_row_labels']])) heat_param[['show_row_labels']] else FALSE
          if (show_row_labels) {
            debug_cat(paste0("  Adding row labels...\n"))

            row_label_source <- if (!is.null(heat_param[['row_label_source']])) heat_param[['row_label_source']] else "colnames"
            row_label_font_size <- if (!is.null(heat_param[['row_label_font_size']])) heat_param[['row_label_font_size']] else 2.5
            custom_row_labels <- if (!is.null(heat_param[['custom_row_labels']])) heat_param[['custom_row_labels']] else ""
            label_mapping <- if (!is.null(heat_param[['label_mapping']])) heat_param[['label_mapping']] else list()

            # Determine labels to use
            if (row_label_source == "mapping" && length(label_mapping) > 0) {
              # v108: Use per-column label mapping
              labels_to_use <- sapply(colnames(heat_data), function(col_name) {
                if (!is.null(label_mapping[[col_name]]) && nchar(label_mapping[[col_name]]) > 0) {
                  label_mapping[[col_name]]
                } else {
                  col_name  # Default to column name if not mapped
                }
              })
              debug_cat(paste0("  Using label mapping for ", length(labels_to_use), " columns\n"))
            } else if (row_label_source == "custom" && nchar(custom_row_labels) > 0) {
              labels_to_use <- trimws(strsplit(custom_row_labels, ",")[[1]])
              # Pad or truncate to match number of columns
              if (length(labels_to_use) < ncol(heat_data)) {
                labels_to_use <- c(labels_to_use, rep("", ncol(heat_data) - length(labels_to_use)))
              } else if (length(labels_to_use) > ncol(heat_data)) {
                labels_to_use <- labels_to_use[1:ncol(heat_data)]
              }
            } else {
              # Use column names
              labels_to_use <- colnames(heat_data)
            }

            debug_cat(paste0("  Row labels: ", paste(labels_to_use, collapse=", "), "\n"))

            # v113: Improved row labels positioning
            # Labels appear below the heatmap (at lower y values than the tips)
            # With scale_y_reverse + coord_flip, this places them visually to the right

            # v113: Get row label offset and alignment from heat_param
            row_label_offset <- if (!is.null(heat_param[['row_label_offset']])) as.numeric(heat_param[['row_label_offset']]) else 1.0
            row_label_align <- if (!is.null(heat_param[['row_label_align']])) heat_param[['row_label_align']] else "left"

            debug_cat(paste0("  v13: row_label_offset from heat_param: ", row_label_offset, "\n"))
            debug_cat(paste0("  v13: row_label_align from heat_param: ", row_label_align, "\n"))

            # v113: Calculate label y position based on offset
            # Labels go below minimum tip y (which is typically 1)
            # Offset controls how far below: 0 = at the edge of tiles, 5 = far from tiles
            label_y_pos <- min(tile_df$y) - (tile_height / 2) - row_label_offset

            # v109: Get colnames angle from heat_param if available
            colnames_angle <- if (!is.null(heat_param[['colnames_angle']])) as.numeric(heat_param[['colnames_angle']]) else 0

            # Create label data - one label per column (visual "row")
            label_df <- data.frame(
              x = numeric(ncol(heat_data)),  # Will be set per column
              y = label_y_pos,               # All labels at same y position
              label = labels_to_use,
              stringsAsFactors = FALSE
            )

            # Set x position for each label to match its column position
            for (col_idx in 1:ncol(heat_data)) {
              col_tiles <- tile_df[tile_df$column == colnames(heat_data)[col_idx], ]
              if (nrow(col_tiles) > 0) {
                # Use the column's x position (all tiles in a column have same x)
                label_df$x[col_idx] <- unique(col_tiles$x)[1]
              }
            }

            # v113: Determine hjust based on alignment setting
            # In coord_flip context with angle=0:
            # - hjust controls visual horizontal alignment
            # - left=0 (text starts at anchor), center=0.5, right=1 (text ends at anchor)
            # In coord_flip context with angle=90 or angle=45:
            # - vjust becomes more relevant for visual positioning
            hjust_val <- switch(row_label_align,
                                "left" = 0,
                                "center" = 0.5,
                                "right" = 1,
                                0)  # default to left

            # v113: For angled text, also adjust vjust for better alignment
            vjust_val <- if (colnames_angle == 0) 0.5 else if (colnames_angle > 0) 1 else 0

            # Add text labels
            p_with_tiles <- p_with_tiles + geom_text(
              data = label_df,
              aes(x = x, y = y, label = label),
              size = row_label_font_size,
              hjust = hjust_val,
              vjust = vjust_val,
              angle = colnames_angle,
              inherit.aes = FALSE
            )
            debug_cat(paste0("  Row labels added successfully\n"))
            debug_cat(paste0("  Label position: label_y_pos=", label_y_pos, ", angle=", colnames_angle, ", hjust=", hjust_val, ", vjust=", vjust_val, "\n"))
          }

          # v116/v119: Add tip guide lines (vertical lines from tips through heatmap)
          show_guides <- if (!is.null(heat_param[['show_guides']])) heat_param[['show_guides']] else FALSE
          show_guides_bool <- isTRUE(show_guides) || identical(show_guides, TRUE) || identical(show_guides, "yes") || identical(show_guides, "TRUE")

          # v119: Debug output to trace guide line settings
          debug_cat(paste0("\n=== v119: TIP GUIDE LINES CHECK ===\n"))
          debug_cat(paste0("  show_guides raw value: ", show_guides, " (class: ", class(show_guides), ")\n"))
          debug_cat(paste0("  show_guides_bool: ", show_guides_bool, "\n"))

          if (show_guides_bool) {
            debug_cat(paste0("\n=== v116: ADDING TIP GUIDE LINES ===\n"))

            # Get guide line settings
            guide_color1 <- if (!is.null(heat_param[['guide_color1']])) heat_param[['guide_color1']] else "#CCCCCC"
            guide_color2 <- if (!is.null(heat_param[['guide_color2']])) heat_param[['guide_color2']] else "#EEEEEE"
            guide_alpha <- if (!is.null(heat_param[['guide_alpha']])) as.numeric(heat_param[['guide_alpha']]) else 0.3
            guide_width <- if (!is.null(heat_param[['guide_width']])) as.numeric(heat_param[['guide_width']]) else 0.5
            guide_linetype <- if (!is.null(heat_param[['guide_linetype']])) heat_param[['guide_linetype']] else "solid"

            debug_cat(paste0("  guide_color1: ", guide_color1, "\n"))
            debug_cat(paste0("  guide_color2: ", guide_color2, "\n"))
            debug_cat(paste0("  guide_alpha: ", guide_alpha, "\n"))
            debug_cat(paste0("  guide_width: ", guide_width, "\n"))
            debug_cat(paste0("  guide_linetype: ", guide_linetype, "\n"))

            # v125: Get tip positions from tip_data to start guide lines at actual tip locations
            # Each tip may have a different x position based on branch lengths
            n_tips <- nrow(tip_data)
            debug_cat(paste0("  Number of tips: ", n_tips, "\n"))

            # v125: Calculate x-end (right edge of heatmap)
            x_max <- max(tile_df$x) + tile_width / 2  # End at right edge of heatmap

            # v125: Build guide lines data using each tip's actual x position
            # This connects each guide line directly to its tree tip
            guide_lines_list <- lapply(seq_len(nrow(tip_data)), function(tip_idx) {
              tip_label <- tip_data$label[tip_idx]
              tip_x <- tip_data$x[tip_idx]  # Actual x position of this tip
              tip_y <- tip_data$y[tip_idx]  # Actual y position of this tip

              data.frame(
                x = tip_x,  # v125: Start at actual tip x position
                xend = x_max,
                y = tip_y,
                yend = tip_y,
                color_idx = tip_idx %% 2  # 0 for even, 1 for odd
              )
            })
            guide_lines_df <- do.call(rbind, guide_lines_list)

            debug_cat(paste0("  v25: Guide lines from individual tip x positions to x_max=", x_max, "\n"))
            debug_cat(paste0("  v25: Tip x range: [", min(tip_data$x), ", ", max(tip_data$x), "]\n"))

            # Apply transparency to the colors
            guide_color1_alpha <- adjustcolor(guide_color1, alpha.f = guide_alpha)
            guide_color2_alpha <- adjustcolor(guide_color2, alpha.f = guide_alpha)

            # Add the guide lines as a layer
            p_with_tiles <- p_with_tiles +
              geom_segment(
                data = guide_lines_df[guide_lines_df$color_idx == 0, ],
                aes(x = x, xend = xend, y = y, yend = yend),
                color = guide_color1_alpha,
                linewidth = guide_width,
                linetype = guide_linetype,
                inherit.aes = FALSE
              ) +
              geom_segment(
                data = guide_lines_df[guide_lines_df$color_idx == 1, ],
                aes(x = x, xend = xend, y = y, yend = yend),
                color = guide_color2_alpha,
                linewidth = guide_width,
                linetype = guide_linetype,
                inherit.aes = FALSE
              )

            debug_cat(paste0("  v16: Added ", n_tips, " tip guide lines\n"))
          }

          debug_cat(paste0("  Final layers: ", length(p_with_tiles$layers), "\n"))
          p_with_tiles

        }, error = function(e) {
          debug_cat(paste0("  ERROR adding heatmap: ", e$message, "\n"))
          debug_cat(paste0("  Returning tree without heatmap\n"))
          p
        })

        debug_cat(paste0("=== v122: HEATMAP ", heat_idx, " COMPLETE ===\n"))
        debug_cat(paste0("  p layers after heatmap: ", length(p$layers), "\n"))
        debug_cat(paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))

        # v125: Update current_heatmap_x_start for next heatmap
        # Calculate the rightmost x position of this heatmap + global gap
        this_heatmap_x_end <- max(tile_df$x) + tile_width / 2
        current_heatmap_x_start <- this_heatmap_x_end + heatmap_spacing  # v125: Add gap between heatmaps
        debug_cat(paste0("  v25: Updated current_heatmap_x_start to ", current_heatmap_x_start, " (added gap=", heatmap_spacing, ")\n"))

      } else {
        debug_cat(paste0("  WARNING: No tile data created - skipping heatmap ", heat_idx, "\n"))
      }
    } # End of v122 for loop for this heatmap
  }
  # END v122 MULTIPLE HEATMAPS

  # v94: Track p right after heatmap block
  debug_cat(paste0("\n=== v94: Immediately after simplified heatmap block ===\n"))
  debug_cat(paste0("  heat_flag: ", heat_flag, "\n"))
  debug_cat(paste0("  p layers: ", length(p$layers), "\n"))
  debug_cat(paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))

  # ========================================================================
  # v91: ORIGINAL COMPLEX HEATMAP CODE COMMENTED OUT BELOW

  # The following 800+ lines of complex heatmap logic have been commented out
  # to simplify debugging. This will be rebuilt step by step.
  # ========================================================================

  if (FALSE) {
  # ORIGINAL CODE START (commented out for v91)
  # v58: FIX - Also check that dxdf440_for_heat has data, not just heat_flag
  # if (heat_flag == TRUE && length(dxdf440_for_heat) > 0) {
    tt_DISABLED <- p

    # v63: DEBUG - show x range BEFORE scaling
    x_range_before <- range(tt_DISABLED$data$x, na.rm = TRUE)
    debug_cat(paste0("\n=== v63: Pre-scaling x range: [", x_range_before[1], ", ", x_range_before[2], "] ===\n"))

    for (i in cc_totss) {
      par <- tt$data$parent[i]
      parX <- tt$data$x[par]
      tt$data[tt$data$node[i], "x"] <- tt$data[tt$data$node[i], "x"] * 15
    }

    # v63: DEBUG - show x range AFTER scaling
    x_range_after <- range(tt$data$x, na.rm = TRUE)
    debug_cat(paste0("=== v63: Post-scaling x range: [", x_range_after[1], ", ", x_range_after[2], "] ===\n"))
    
    how_many_tips <- length(tt$data$isTip)
    
    if (debug_mode == TRUE) {
      # v53: print("how_many_tips is")
      # v53: print(how_many_tips)
    }
    
    # Calculate max column name length for formatting
    # v86: CRITICAL FIX - renamed loop variable from 'j' to 'h_idx' to avoid conflict with
    # the heatmap counter 'j' that's used later in the code for tracking enabled heatmaps.
    # Using 'j' here was overwriting the counter and causing j=2 instead of j=1 for the
    # first heatmap, which led to incorrect scale application and "Problem while setting up geom" errors.
    max_len_col_name_heat <- 0
    for (h_idx in 1:length(heat_map_title_list)) {
      colnames_heat <- colnames(dxdf440_for_heat[[h_idx]])
      max_len_col_name_heat <- max(max_len_col_name_heat, max(nchar(colnames_heat)))
    }
    
    # Get plot boundaries for positioning
    boudariestt <- func.find.plot.boundaries(tt, debug_mode)
    
    if (debug_mode == TRUE) {
      # v53: print("boudariestt is")
      # v53: print(boudariestt)
    }
    
    p <- tt 
    off_base <- 0.7
    param <- 0
    j <- 0
    
    # Process each heatmap
    for (j1 in 1:length(heat_map_title_list)) {
      if (debug_mode == TRUE) {
        # v53: print("heatmap number")
        # v53: print(j)
        # v53: print("off_base is")
        # v53: print(off_base)
        # v53: print("param is")
        # v53: print(param)
        # v53: print("man_adj_heat_loc is")
        # v53: print(man_adj_heat_loc)
      }
      
      # Skip if this heatmap is disabled
      if (heat_display_vec[j1] == FALSE) {
        pr440_short_tips_TRY_heat <- tt
        next
      } else {
        j <- j + 1
      }
      
      # Calculate offset for this heatmap
      off <- off_base + param + man_adj_heat_loc
      
      if (j > 1) {
        off <- off + man_adj_heat_loc2
        tt <- pr440_short_tips_TRY_heat + new_scale_fill()
        # v69: Repair mapping after new_scale_fill (can cause corruption)
        tt <- func.repair.ggtree.mapping(tt)
      }

      if (j > 2) {
        off <- off + man_adj_heat_loc3
      }
      
      # Calculate heatmap width
      if (is.na(width_heatmap[1])) {
        width_heat_row <- 1.4
        wi <- width_heat_row * length(colnames(dxdf440_for_heat[[j1]])) / 23
      } else {
        wi <- width_heatmap[j1]
      }
      
      if (debug_mode == TRUE) {
        # v53: print("off is")
        # v53: print(off)
        # v53: print("length(colnames(dxdf440_for_heat[[j1]])) is")
        # v53: print(length(colnames(dxdf440_for_heat[[j1]])))
        # v53: print("wi is")
        # v53: print(wi)
        # v53: print("range is is")
        # v53: print(tt$data)
      }
      
      # Skip if no data for this heatmap
      if (is.null(dxdf440_for_heat[[j1]]) == TRUE) {
        # v53: print(paste0("No columns for heatmap display ", j1))
        next
      }
      
      # Set up positioning for this heatmap
      sub_df_heat_j <- dxdf440_for_heat[[j1]]
      colnames_sub_df_heat_j <- colnames(sub_df_heat_j)
      
      if (j > 1) {
        new_heat_x <- min_col_x_of_frame_of_prev_heat + 1.5
      } else {
        new_heat_x <- off
      }
      
      boudariestt <- func.find.plot.boundaries(tt, debug_mode)
      
      # Calculate offset for column names
      heat_names_offset <- how_many_tips * max_len_col_name_heat / 380 * (size_font_heat_map_legend / 0.8) - 0.3
      
      # Get heatmap parameters
      heat_param_j <- heat_display_params_list[[j1]]
      heat_param <- heat_display_params_list[[j1]]
      
      if (debug_mode == TRUE) {
        # v53: print("heat_names_offset is")
        # v53: print(heat_names_offset)
        # v53: print("heat_param_j is")
        # v53: print(heat_param_j)
        # v53: print("new_heat_x is")
        # v53: print(new_heat_x)
      }
      
      # Set angle for column labels
      if (is.na(heat_maps_titles_angles_vector[1])) {
        colnames_angle <- 0
      } else {
        if (length(heat_maps_titles_angles_vector) < j1) {
          colnames_angle <- 0
        } else {
          colnames_angle <- heat_maps_titles_angles_vector[j1]
        }
      }
      
      p <- tt
      
      # Handle column names
      # S1.62dev: Use per-heatmap show_colnames from heat_param if available, else fall back to flag_colnames
      show_colnames_for_heatmap <- TRUE  # default to showing
      if (!is.null(heat_param) && 'show_colnames' %in% names(heat_param)) {
        show_colnames_for_heatmap <- heat_param[['show_colnames']]
      } else if (!is.na(flag_colnames[1]) && j1 <= length(flag_colnames)) {
        show_colnames_for_heatmap <- flag_colnames[j1]
      }

      if (show_colnames_for_heatmap == FALSE) {
        dt <- dxdf440_for_heat[[j1]]
        colnames_len <- length(colnames(dt))
        custom_column_labels <- rep("", colnames_len)
      } else {
        dt <- dxdf440_for_heat[[j1]]
        custom_column_labels <- colnames(dt)
      }
      
      # Replace heat legend labels if requested
      # v87: CRITICAL FIX - renamed loop variable from 'j' to 'legend_key' to avoid overwriting heatmap counter
      if (!is.na(heat_legend_replace[1])) {
        for (legend_key in names(heat_legend_replace)) {
          if (legend_key %in% custom_column_labels) {
            custom_column_labels[custom_column_labels == legend_key] <- heat_legend_replace[legend_key]
          }
        }
      }

      # Calculate max column name length again for formatting
      # v87: CRITICAL FIX - renamed loop variable from 'j' to 'col_idx2' to avoid overwriting heatmap counter
      max_len_col_name_heat <- 0
      for (col_idx2 in 1:length(custom_column_labels)) {
        max_len_col_name_heat <- max(max_len_col_name_heat, max(nchar(custom_column_labels)))
      }

      # v89: RESTORED factor conversion for discrete heatmap data
      # The v88 change to skip factor conversion was WRONG - gheatmap needs factors
      # for discrete data to work properly with ggplot2's discrete color scales.
      # Without factors, "Problem while setting up geom" errors occur during rendering.
      if (!is.null(heat_param) && heat_param['is_discrete'] == TRUE) {
        debug_cat(paste0("\n=== v89: Converting discrete heatmap to factors ===\n"))
        for (col_idx in 1:ncol(dxdf440_for_heat[[j1]])) {
          col_name <- colnames(dxdf440_for_heat[[j1]])[col_idx]
          col_vals <- dxdf440_for_heat[[j1]][, col_idx]
          unique_vals <- sort(unique(na.omit(col_vals)))
          debug_cat(paste0("  Column '", col_name, "': ", length(unique_vals), " unique values\n"))
          debug_cat(paste0("  Unique values: ", paste(head(unique_vals, 10), collapse=", "),
                                    if(length(unique_vals) > 10) "..." else "", "\n"))
          # v89: RESTORED - Convert to factor with sorted levels (REQUIRED for discrete heatmaps)
          dxdf440_for_heat[[j1]][, col_idx] <- factor(col_vals, levels = unique_vals)
          debug_cat(paste0("  Converted to factor with ", length(unique_vals), " levels\n"))
        }
        debug_cat(paste0("================================\n"))
      }

      # Create the heatmap
      # v61: DEBUG - show data structure before gheatmap call
      debug_cat(paste0("\n=== v61: GHEATMAP DATA DEBUG ===\n"))
      debug_cat(paste0("  Heatmap index: j1=", j1, ", j=", j, "\n"))
      heat_data <- dxdf440_for_heat[[j1]]
      debug_cat(paste0("  heat_data dimensions: ", nrow(heat_data), " x ", ncol(heat_data), "\n"))
      debug_cat(paste0("  heat_data rownames sample: ", paste(head(rownames(heat_data), 5), collapse=", "), "\n"))
      debug_cat(paste0("  heat_data columns: ", paste(colnames(heat_data), collapse=", "), "\n"))
      tt_tips <- subset(tt$data, isTip == TRUE)
      debug_cat(paste0("  Tree tip labels sample: ", paste(head(tt_tips$label, 5), collapse=", "), "\n"))
      # Check if rownames match tree tip labels
      matching_tips <- sum(rownames(heat_data) %in% tt_tips$label)
      debug_cat(paste0("  Rownames matching tree tips: ", matching_tips, " / ", nrow(heat_data), "\n"))
      debug_cat(paste0("  offset (new_heat_x): ", new_heat_x, "\n"))
      debug_cat(paste0("  width (wi): ", wi, "\n"))
      debug_cat(paste0("================================\n"))

      # v83: RESTORED duplicate gheatmap pattern from v61 - THIS IS INTENTIONAL
      # The user confirmed this pattern was in the original lineage plotter code and is required.
      # DO NOT REMOVE THIS DUPLICATE CALL - it is necessary for gheatmap to work correctly.
      # First gheatmap call creates the initial structure
      .prof_gheatmap_start <- Sys.time()
      pr440_short_tips_TRY_heat <- gheatmap(
        tt,
        data = dxdf440_for_heat[[j1]],
        colnames_angle = colnames_angle,
        offset = new_heat_x,
        width = wi,
        font.size = size_font_heat_map_legend,
        colnames_offset_x = 0,
        colnames_offset_y = heat_names_offset,
        legend_title = heat_map_title_list[[j1]],
        colnames = TRUE,
        custom_column_labels = custom_column_labels,
        color = NA
      )

      # v83: IMMEDIATELY repair mapping after first gheatmap call
      # This MUST happen before any other operations on the plot
      pr440_short_tips_TRY_heat <- func.repair.ggtree.mapping(pr440_short_tips_TRY_heat, verbose = FALSE)

      # v83: REQUIRED second gheatmap call - DO NOT REMOVE
      # For j==1 (first heatmap):
      #   - Continuous: call gheatmap on result (pr440_short_tips_TRY_heat)
      #   - Discrete: call gheatmap on original tree (tt)
      # This duplicate call is intentional and required for proper rendering.
      if (j == 1) {
        if (heat_param['is_discrete'] == FALSE) {
          # Continuous heatmaps: apply gheatmap on the result
          pr440_short_tips_TRY_heat <- gheatmap(
            pr440_short_tips_TRY_heat,
            data = dxdf440_for_heat[[j1]],
            colnames_angle = colnames_angle,
            offset = new_heat_x,
            width = wi,
            font.size = size_font_heat_map_legend,
            colnames_offset_x = 0,
            colnames_offset_y = heat_names_offset,
            legend_title = heat_map_title_list[[j1]],
            colnames = TRUE,
            custom_column_labels = custom_column_labels,
            color = NA
          )
        } else {
          # Discrete heatmaps: apply gheatmap on original tree again
          pr440_short_tips_TRY_heat <- gheatmap(
            tt,
            data = dxdf440_for_heat[[j1]],
            colnames_angle = colnames_angle,
            offset = new_heat_x,
            width = wi,
            font.size = size_font_heat_map_legend,
            colnames_offset_x = 0,
            colnames_offset_y = heat_names_offset,
            legend_title = heat_map_title_list[[j1]],
            colnames = TRUE,
            custom_column_labels = custom_column_labels,
            color = NA
          )
        }
        # v83: IMMEDIATELY repair mapping after second gheatmap call
        pr440_short_tips_TRY_heat <- func.repair.ggtree.mapping(pr440_short_tips_TRY_heat, verbose = FALSE)
      }
      cat(file=stderr(), sprintf("[PROF-TREE] gheatmap (heatmap %d): %.3f sec\n", j, as.numeric(Sys.time() - .prof_gheatmap_start)))

      # v85: REMOVED direct layer data modification - it corrupts ggplot2 internal state
      # gheatmap creates tile layers with character values, and scale_fill_manual can handle
      # both character and factor values. Direct modification of layer$data causes
      # "Problem while setting up geom" errors during ggplot_build.
      #
      # Instead, we collect the unique values for scale setup without modifying layer data.
      tile_values <- c()
      for (layer_idx in seq_along(pr440_short_tips_TRY_heat$layers)) {
        layer <- pr440_short_tips_TRY_heat$layers[[layer_idx]]
        if (inherits(layer$geom, "GeomTile") && !is.null(layer$data) && is.data.frame(layer$data)) {
          if ("value" %in% names(layer$data)) {
            layer_vals <- unique(na.omit(layer$data$value))
            tile_values <- c(tile_values, layer_vals)
            debug_cat(paste0("  v5: Found tile layer ", layer_idx, " with values: ",
                                      paste(head(layer_vals, 5), collapse=", "), "\n"))
          }
        }
      }
      tile_values <- unique(tile_values)
      debug_cat(paste0("  v5: All tile values: ", paste(tile_values, collapse=", "), "\n"))

      # v82: DEBUG - verify gheatmap result and tile layer data
      debug_cat(paste0("\n=== v71: POST-GHEATMAP DEBUG ===\n"))
      debug_cat(paste0("  Number of layers in plot: ", length(pr440_short_tips_TRY_heat$layers), "\n"))
      gheatmap_xrange <- range(pr440_short_tips_TRY_heat$data$x, na.rm = TRUE)
      debug_cat(paste0("  Plot data x range: [", gheatmap_xrange[1], ", ", gheatmap_xrange[2], "]\n"))

      # Check if there's rect/tile data (heatmap)
      layer_types <- sapply(pr440_short_tips_TRY_heat$layers, function(l) class(l$geom)[1])
      debug_cat(paste0("  Layer geom types: ", paste(layer_types, collapse=", "), "\n"))

      # v72: Enhanced debugging to find and inspect the GeomTile layer
      for (layer_idx in seq_along(pr440_short_tips_TRY_heat$layers)) {
        layer <- pr440_short_tips_TRY_heat$layers[[layer_idx]]
        if (inherits(layer$geom, "GeomTile")) {
          debug_cat(paste0("  GeomTile found at layer ", layer_idx, "\n"))
          tryCatch({
            layer_data <- layer$data
            if (is.function(layer_data)) {
              debug_cat(paste0("    Layer data is a function, evaluating...\n"))
              layer_data <- layer_data(pr440_short_tips_TRY_heat$data)
            }
            if (!is.null(layer_data) && is.data.frame(layer_data)) {
              debug_cat(paste0("    Layer data rows: ", nrow(layer_data), "\n"))
              debug_cat(paste0("    Layer data columns: ", paste(names(layer_data), collapse=", "), "\n"))
              if ("x" %in% names(layer_data)) {
                x_vals <- layer_data$x
                debug_cat(paste0("    Tile x range: [", min(x_vals, na.rm=TRUE), ", ", max(x_vals, na.rm=TRUE), "]\n"))
              }
              if ("y" %in% names(layer_data)) {
                y_vals <- layer_data$y
                debug_cat(paste0("    Tile y range: [", min(y_vals, na.rm=TRUE), ", ", max(y_vals, na.rm=TRUE), "]\n"))
              }
              # v72: Check value column (this is what fill maps to in gheatmap)
              if ("value" %in% names(layer_data)) {
                debug_cat(paste0("    Value column (unique): ", paste(unique(layer_data$value), collapse=", "), "\n"))
                debug_cat(paste0("    Value column class: ", class(layer_data$value)[1], "\n"))
              }
              # v72: Check width column (tile width)
              if ("width" %in% names(layer_data)) {
                width_vals <- layer_data$width
                debug_cat(paste0("    Width values: [", min(width_vals, na.rm=TRUE), ", ", max(width_vals, na.rm=TRUE), "]\n"))
              }
              if ("fill" %in% names(layer_data)) {
                debug_cat(paste0("    Fill values (first 5): ", paste(head(layer_data$fill, 5), collapse=", "), "\n"))
              }
              # v72: Check the layer's aesthetic mapping
              if (!is.null(layer$mapping)) {
                debug_cat(paste0("    Layer mapping: ", paste(names(layer$mapping), collapse=", "), "\n"))
              }
            } else {
              debug_cat(paste0("    WARNING: Layer data is not a data.frame: ", class(layer_data)[1], "\n"))
            }
          }, error = function(e) {
            debug_cat(paste0("    ERROR accessing layer data: ", e$message, "\n"))
          })
        }
      }

      # v82: Removed aggressive ggplot_build debugging that was causing premature errors
      # The build will be done at the end when saving the plot

      # Apply correct coloring scale based on heatmap type
      if (heat_param['is_discrete'] == FALSE) {
        limits <- heat_param[['limits']]

        if (is.na(limits[1]) == TRUE) {
          pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
            scale_fill_gradient2(
              low = heat_param['low'],
              mid = heat_param['mid'],
              high = heat_param['high'],
              midpoint = .02,
              name = heat_map_title_list[[j1]],
              na.value = "white"
            )
        } else {
          pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
            scale_fill_gradient2(
              low = heat_param['low'],
              mid = heat_param['mid'],
              high = heat_param['high'],
              midpoint = .02,
              name = heat_map_title_list[[j1]],
              limits = limits,
              na.value = "white"
            )
        }
        # v82: Removed aggressive mapping repair - will do single repair at end of heatmap loop
      } else {
        # v65: DEBUG - trace discrete heatmap color path
        debug_cat(paste0("\n=== v65: DISCRETE HEATMAP COLOR DEBUG ===\n"))
        debug_cat(paste0("  heat_param['is_discrete']: ", heat_param['is_discrete'], "\n"))
        debug_cat(paste0("  heat_param['man']: ", heat_param['man'], "\n"))
        debug_cat(paste0("  heat_param['man_define_colors']: ", heat_param['man_define_colors'], "\n"))
        debug_cat(paste0("  heat_param[['color_scale_option']]: ",
                                  if(is.null(heat_param[['color_scale_option']])) "NULL"
                                  else paste(class(heat_param[['color_scale_option']]), collapse=", "), "\n"))
        if (!is.null(heat_param[['color_scale_option']])) {
          debug_cat(paste0("  color_scale_option value: ",
                                    paste(heat_param[['color_scale_option']], collapse=", "), "\n"))
          if (is.list(heat_param[['color_scale_option']])) {
            debug_cat(paste0("  color_scale_option$color_scale_option: ",
                                      heat_param[['color_scale_option']]$color_scale_option, "\n"))
          }
        }
        debug_cat(paste0("========================================\n"))

        if (heat_param['man'] == FALSE) {
          if (heat_param['man_define_colors'] == FALSE) {
            # v67: SIMPLIFIED - color_scale_option is now directly a palette name string (e.g., "Set1")
            # No more nested list handling needed after the fix at line 2616
            palette_name <- heat_param[['color_scale_option']]

            debug_cat(paste0("\n=== v67: Discrete palette check ===\n"))
            debug_cat(paste0("  palette_name: ", if(is.null(palette_name)) "NULL" else palette_name, "\n"))
            if (!is.null(palette_name)) {
              debug_cat(paste0("  Is valid RColorBrewer palette: ",
                                        palette_name %in% rownames(RColorBrewer::brewer.pal.info), "\n"))
            }
            debug_cat(paste0("================================\n"))

            # v67: Use RColorBrewer palette if available, otherwise use default hue
            if (!is.null(palette_name) && is.character(palette_name) &&
                palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
              # Get unique values to determine number of colors needed
              heat_data_vals <- unique(na.omit(dxdf440_for_heat[[j1]][,1]))
              n_vals <- length(heat_data_vals)
              max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
              n_colors <- min(n_vals, max_colors)

              debug_cat(paste0("\n=== v67: Applying discrete palette ===\n"))
              debug_cat(paste0("  Palette: ", palette_name, "\n"))
              debug_cat(paste0("  Number of unique values: ", n_vals, "\n"))
              debug_cat(paste0("  Colors to use: ", n_colors, "\n"))
              debug_cat(paste0("================================\n"))

              # v70: Get NA color from heat_param (default white)
              na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "white"
              debug_cat(paste0("  NA color: ", na_color, "\n"))

              pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
                scale_fill_brewer(palette = palette_name, name = heat_map_title_list[[j1]], na.value = na_color)
              # v82: Removed aggressive mapping repair - will do single repair at end of heatmap loop
            } else {
              debug_cat(paste0("  v7: Using default hue scale (no valid palette specified)\n"))
              # v70: Get NA color from heat_param (default white)
              na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "white"
              pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
                scale_fill_hue(name = heat_map_title_list[[j1]], na.value = na_color)
              # v82: Removed aggressive mapping repair - will do single repair at end of heatmap loop
            }
          } else {
            # v70: man_define_colors is TRUE - use custom color values
            # color_scale_option should be a named vector of colors
            custom_colors <- heat_param[['color_scale_option']]
            debug_cat(paste0("\n=== v70: Applying custom discrete colors ===\n"))
            debug_cat(paste0("  custom_colors class: ", paste(class(custom_colors), collapse=", "), "\n"))
            debug_cat(paste0("  custom_colors length: ", length(custom_colors), "\n"))
            debug_cat(paste0("  custom_colors names: ", paste(head(names(custom_colors), 10), collapse=", "), "\n"))
            debug_cat(paste0("  custom_colors values: ", paste(head(custom_colors, 10), collapse=", "), "\n"))

            # v85: Use tile_values collected from actual tile layer data (not dxdf440_for_heat)
            # This ensures we match exactly what gheatmap created
            heat_data_vals <- tile_values
            debug_cat(paste0("  v5: Using tile_values from layer: ", paste(head(heat_data_vals, 10), collapse=", "), "\n"))

            # v73: Fix - properly subset custom_colors to match tile values
            # This prevents NA names which cause scale_fill_manual to fail
            n_levels <- length(heat_data_vals)

            if (n_levels == 0) {
              # v85: If no tile values found, fall back to using dxdf440_for_heat
              debug_cat(paste0("  v5: WARNING - no tile values found, using source data\n"))
              for (col_idx in 1:ncol(dxdf440_for_heat[[j1]])) {
                col_vals <- dxdf440_for_heat[[j1]][, col_idx]
                if (is.factor(col_vals)) {
                  col_levels <- levels(col_vals)
                  if (!is.null(col_levels) && length(col_levels) > 0) {
                    heat_data_vals <- c(heat_data_vals, col_levels)
                  }
                } else {
                  col_unique <- unique(na.omit(col_vals))
                  if (length(col_unique) > 0) {
                    heat_data_vals <- c(heat_data_vals, col_unique)
                  }
                }
              }
              heat_data_vals <- unique(heat_data_vals)
              n_levels <- length(heat_data_vals)
            }

            if (is.null(names(custom_colors)) || length(names(custom_colors)) == 0) {
              debug_cat(paste0("  v3: WARNING - custom_colors has no names, subsetting and assigning by position\n"))
              debug_cat(paste0("  v3: Number of tile values: ", n_levels, "\n"))
              debug_cat(paste0("  v3: Number of custom colors: ", length(custom_colors), "\n"))

              # CRITICAL: Only use as many colors as there are values
              if (n_levels > 0 && length(custom_colors) >= n_levels) {
                # Subset colors to match values exactly
                colors_to_use <- custom_colors[1:n_levels]
                names(colors_to_use) <- as.character(heat_data_vals)
                custom_colors <- colors_to_use
                debug_cat(paste0("  v3: Subsetted to ", n_levels, " colors\n"))
                debug_cat(paste0("  v3: Final names: ", paste(names(custom_colors), collapse=", "), "\n"))
                debug_cat(paste0("  v3: Final values: ", paste(custom_colors, collapse=", "), "\n"))
              } else if (n_levels > 0) {
                # Not enough colors, recycle
                debug_cat(paste0("  v3: WARNING - not enough colors, recycling\n"))
                colors_to_use <- rep(custom_colors, length.out = n_levels)
                names(colors_to_use) <- as.character(heat_data_vals)
                custom_colors <- colors_to_use
              }
            } else {
              # custom_colors already has names - ensure they match values
              debug_cat(paste0("  v3: custom_colors already named: ", paste(names(custom_colors), collapse=", "), "\n"))
              # Only keep colors whose names are in heat_data_vals
              valid_names <- names(custom_colors) %in% as.character(heat_data_vals)
              if (any(valid_names)) {
                custom_colors <- custom_colors[valid_names]
                debug_cat(paste0("  v3: After filtering: ", paste(names(custom_colors), collapse=", "), "\n"))
              }
            }
            # v70: Get NA color from heat_param (default white)
            na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "white"
            debug_cat(paste0("  NA color: ", na_color, "\n"))
            debug_cat(paste0("================================\n"))

            # v88: CRITICAL FIX - Test ggplot_build BEFORE applying scale to diagnose root cause
            debug_cat(paste0("  v8: Testing ggplot_build BEFORE scale application...\n"))
            pre_scale_ok <- tryCatch({
              test_build_pre <- ggplot2::ggplot_build(pr440_short_tips_TRY_heat)
              debug_cat(paste0("  v8: Pre-scale ggplot_build: SUCCESS\n"))
              TRUE
            }, error = function(e) {
              debug_cat(paste0("  v8: Pre-scale ggplot_build: FAILED - ", e$message, "\n"))
              debug_cat(paste0("  v8: ERROR IS IN GHEATMAP OUTPUT, NOT SCALE\n"))
              # Get more info about the error
              debug_cat(paste0("  v8: Checking individual layers...\n"))
              for (li in seq_along(pr440_short_tips_TRY_heat$layers)) {
                layer_test <- tryCatch({
                  # Try to compute layer data
                  layer <- pr440_short_tips_TRY_heat$layers[[li]]
                  if (is.function(layer$data)) {
                    test_d <- layer$data(pr440_short_tips_TRY_heat$data)
                  } else {
                    test_d <- layer$data
                  }
                  # Check if layer has valid mapping
                  if (!is.null(layer$mapping)) {
                    mapping_names <- names(layer$mapping)
                    debug_cat(paste0("    Layer ", li, " (", class(layer$geom)[1], "): mapping=", paste(mapping_names, collapse=","), "\n"))
                  } else {
                    debug_cat(paste0("    Layer ", li, " (", class(layer$geom)[1], "): no mapping\n"))
                  }
                  TRUE
                }, error = function(e2) {
                  debug_cat(paste0("    Layer ", li, " (", class(pr440_short_tips_TRY_heat$layers[[li]]$geom)[1], "): FAILED - ", e2$message, "\n"))
                  FALSE
                })
              }
              FALSE
            })

            # v88: Apply scale regardless of pre-test result (the scale itself isn't the problem)
            debug_cat(paste0("  v8: Applying scale_fill_manual...\n"))
            tryCatch({
              pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
                scale_fill_manual(
                  values = custom_colors,
                  name = heat_map_title_list[[j1]],
                  na.value = na_color
                )
              debug_cat(paste0("  v8: scale_fill_manual applied\n"))

              # v88: Test after scale
              post_scale_ok <- tryCatch({
                test_build_post <- ggplot2::ggplot_build(pr440_short_tips_TRY_heat)
                debug_cat(paste0("  v8: Post-scale ggplot_build: SUCCESS\n"))
                TRUE
              }, error = function(e) {
                debug_cat(paste0("  v8: Post-scale ggplot_build: FAILED - ", e$message, "\n"))
                FALSE
              })

            }, error = function(e) {
              debug_cat(paste0("  v8: scale_fill_manual failed: ", e$message, "\n"))
              debug_cat(paste0("  v8: Trying scale_fill_discrete as fallback\n"))
              tryCatch({
                pr440_short_tips_TRY_heat <<- pr440_short_tips_TRY_heat +
                  scale_fill_discrete(name = heat_map_title_list[[j1]], na.value = na_color)
              }, error = function(e2) {
                debug_cat(paste0("  v8: scale_fill_discrete also failed: ", e2$message, "\n"))
              })
            })
          }
        } else {
          pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
            scale_fill_hue(
              h = c(heat_param$color_scale_range_start, heat_param$color_scale_range_end),
              name = heat_map_title_list[[j1]]
            )
        }
      }

      # v68: DEBUG - after scale application
      debug_cat(paste0("\n=== v68: After scale application ===\n"))
      debug_cat(paste0("  j=", j, ", j1=", j1, "\n"))
      debug_cat(paste0("================================\n"))

      # Apply theme for first heatmap
      if (j == 1) {
        pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat + 
          theme(
            legend.text = element_text(
              size = text_legend_size * text_mu,
              margin = margin(t = mar/2, b = mar/2, unit = "pt")
            ),
            legend.title = element_text(size = text_legend_title_size, face = "bold"),
            legend.spacing.y = unit(text_mu/12, 'cm'),
            legend.spacing.x = unit(0.7, 'cm'),
            legend.spacing = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.height = unit(0.8, "line")
          )
        
        off_base <- off
      } else {
        # Apply coloring for additional heatmaps
        if (heat_param['is_discrete'] == FALSE) {
          pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
            scale_fill_gradient2(
              low = heat_param['low'], 
              mid = heat_param['mid'], 
              high = heat_param['high'], 
              midpoint = .02,
              name = heat_map_title_list[[j1]]       
            )
        } else {
          if (j == 2) {
            if (heat_param['is_discrete'] == FALSE) {
              pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
                scale_fill_viridis_d(
                  option = viridis_option_list[j], 
                  name = heat_map_title_list[[j1]], 
                  direction = -1, 
                  na.value = "gray"
                ) +
                scale_fill_hue(name = heat_map_title_list[[j1]])
            }
          } else {
            pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
              scale_fill_viridis_d(
                option = viridis_option_list[j], 
                name = heat_map_title_list[[j1]], 
                direction = 1, 
                na.value = "gray"
              ) +
              scale_fill_hue(name = heat_map_title_list[[j1]])
          }
        }
      }
      
      # Calculate position for next heatmap
      param <- param + 0.25 * length(colnames(dxdf440_for_heat[[j1]])) + 0.25
      
      if (debug_mode == TRUE) {
        # v53: print("param is")
        # v53: print(param)
      }
      
      # v68: DEBUG - before func.calc call
      debug_cat(paste0("\n=== v68: Before func.calc.min_col_x ===\n"))

      # Store position for next heatmap (wrapped in tryCatch to prevent errors from breaking the loop)
      min_col_x_of_frame_of_prev_heat <- tryCatch({
        func.calc.min_col_x_of_frame_of_prev_heat(pr440_short_tips_TRY_heat, colnames_sub_df_heat_j)
      }, error = function(e) {
        debug_cat(paste0("  v8: WARNING - func.calc.min_col_x_of_frame_of_prev_heat failed: ", e$message, "\n"))
        0  # Return default value
      })

      debug_cat(paste0("  v8: min_col_x_of_frame_of_prev_heat = ", min_col_x_of_frame_of_prev_heat, "\n"))
      debug_cat(paste0("================================\n"))
    }

    # v68: DEBUG - after for loop
    debug_cat(paste0("\n=== v68: After heatmap for loop ===\n"))

    # Update plot with heatmap
    p <- pr440_short_tips_TRY_heat

    # v71: Repair mapping before final state check
    p <- func.repair.ggtree.mapping(p, verbose = TRUE)

    # v71: Diagnose any layer issues
    problematic_layers <- func.diagnose.layer.issues(p, verbose = TRUE)

    # v71: Final heatmap state with layer analysis
    debug_cat(paste0("\n=== v71: FINAL HEATMAP STATE ===\n"))
    final_xrange <- range(p$data$x, na.rm = TRUE)
    debug_cat(paste0("  Tree data x range: [", final_xrange[1], ", ", final_xrange[2], "]\n"))
    debug_cat(paste0("  Number of layers: ", length(p$layers), "\n"))

    # v71: Check for heatmap layer (GeomTile) and get its x coordinates
    layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
    debug_cat(paste0("  Layer types: ", paste(layer_types, collapse=", "), "\n"))

    # v71: Find the GeomTile layer (heatmap) and get its data directly
    heatmap_xmax <- NULL
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      if (inherits(layer$geom, "GeomTile")) {
        debug_cat(paste0("  Found GeomTile at layer ", i, "\n"))
        # Try to access the layer data
        tryCatch({
          layer_data <- layer$data
          if (is.function(layer_data)) {
            layer_data <- layer_data(p$data)
          }
          if (!is.null(layer_data) && "x" %in% names(layer_data)) {
            tile_x <- layer_data$x
            if (length(tile_x) > 0 && !all(is.na(tile_x))) {
              tile_xmax <- max(tile_x, na.rm = TRUE)
              debug_cat(paste0("    GeomTile x range: [", min(tile_x, na.rm = TRUE), ", ", tile_xmax, "]\n"))
              heatmap_xmax <- tile_xmax
            }
          } else {
            debug_cat(paste0("    GeomTile layer data does not have x column\n"))
          }
        }, error = function(e) {
          debug_cat(paste0("    Could not access GeomTile data: ", e$message, "\n"))
        })
      }
    }

    # v71: Also try ggplot_build but don't fail if it doesn't work
    if (is.null(heatmap_xmax)) {
      tryCatch({
        # Repair mapping one more time before building
        p <- func.repair.ggtree.mapping(p)
        built <- ggplot2::ggplot_build(p)
        for (i in seq_along(built$data)) {
          if ("x" %in% names(built$data[[i]])) {
            layer_x <- built$data[[i]]$x
            if (length(layer_x) > 0 && !all(is.na(layer_x))) {
              layer_xmax <- max(layer_x, na.rm = TRUE)
              debug_cat(paste0("    Built layer ", i, " x range: [",
                                        min(layer_x, na.rm = TRUE), ", ", layer_xmax, "]\n"))
              if (is.null(heatmap_xmax) || layer_xmax > heatmap_xmax) {
                heatmap_xmax <- layer_xmax
              }
            }
          }
        }
      }, error = function(e) {
        debug_cat(paste0("  v1: ggplot_build failed (will use fallback): ", e$message, "\n"))
      })
    }
    debug_cat(paste0("================================\n"))

    # v71: Calculate expected x range using multiple methods
    # gheatmap places tiles at x positions based on tree width and offset
    tree_width <- abs(final_xrange[2] - final_xrange[1])

    # Method 1: Based on offset and width parameters (gheatmap uses offset relative to tips at x=0)
    # The heatmap should span from x=offset to x=offset+width (approximately)
    calculated_xmax <- new_heat_x + wi + 0.3  # Add margin

    # Method 2: Use a proportion of tree width as margin (safer fallback)
    proportional_xmax <- tree_width * 0.3 + wi + 0.5

    # v71: Method 3: Use a larger fixed margin to ensure heatmap is visible
    fixed_margin_xmax <- 2.5  # Fixed generous margin

    # Use the largest of: calculated, proportional, fixed, or detected from built data
    expected_xmax <- max(
      calculated_xmax,
      proportional_xmax,
      fixed_margin_xmax,
      ifelse(is.null(heatmap_xmax), calculated_xmax, heatmap_xmax + 0.5)
    )

    debug_cat(paste0("\n=== v71: EXPANDING X-AXIS FOR HEATMAP ===\n"))
    debug_cat(paste0("  Tree x range: [", final_xrange[1], ", ", final_xrange[2], "]\n"))
    debug_cat(paste0("  Tree width: ", tree_width, "\n"))
    debug_cat(paste0("  Heatmap offset (new_heat_x): ", new_heat_x, "\n"))
    debug_cat(paste0("  Heatmap width (wi): ", wi, "\n"))
    debug_cat(paste0("  Calculated x max: ", calculated_xmax, "\n"))
    debug_cat(paste0("  Proportional x max: ", proportional_xmax, "\n"))
    debug_cat(paste0("  Fixed margin x max: ", fixed_margin_xmax, "\n"))
    debug_cat(paste0("  Detected from build: ", ifelse(is.null(heatmap_xmax), "NULL", heatmap_xmax), "\n"))
    debug_cat(paste0("  Final expected max x: ", expected_xmax, "\n"))
    debug_cat(paste0("  Setting coord_flip xlim to: [", final_xrange[1], ", ", expected_xmax, "]\n"))
    debug_cat(paste0("========================================\n"))

    # v88: SIMPLIFIED EXPANSION - Skip expansion functions if plot can't be built
    # Test if plot is buildable before trying expansion
    plot_buildable <- tryCatch({
      ggplot2::ggplot_build(p)
      debug_cat(paste0("  v8: Plot is buildable, proceeding with expansion\n"))
      TRUE
    }, error = function(e) {
      debug_cat(paste0("  v8: Plot is NOT buildable: ", e$message, "\n"))
      debug_cat(paste0("  v8: Skipping expansion - will render plot as-is\n"))
      FALSE
    })

    if (plot_buildable) {
      # Calculate how much expansion is needed on the right side (positive x direction)
      expansion_ratio <- (expected_xmax - final_xrange[2]) / abs(final_xrange[1] - final_xrange[2])
      expansion_ratio <- max(0.3, expansion_ratio)

      debug_cat(paste0("  v8: Using hexpand() with ratio: ", expansion_ratio, "\n"))

      # v88: Try hexpand first
      expansion_success <- FALSE
      tryCatch({
        p <- p + ggtree::hexpand(ratio = expansion_ratio, direction = 1)
        expansion_success <- TRUE
        debug_cat(paste0("  v8: hexpand applied successfully\n"))
      }, error = function(e) {
        debug_cat(paste0("  v8: hexpand failed: ", e$message, "\n"))
      })

      # v88: Fallback to xlim_expand
      if (!expansion_success) {
        debug_cat(paste0("  v8: Trying xlim_expand fallback\n"))
        tryCatch({
          p <- p + ggtree::xlim_expand(c(0, expected_xmax), "right")
          expansion_success <- TRUE
          debug_cat(paste0("  v8: xlim_expand applied\n"))
        }, error = function(e2) {
          debug_cat(paste0("  v8: xlim_expand also failed: ", e2$message, "\n"))
        })
      }

      # v88: Final fallback - just add geom_blank with expanded limits
      if (!expansion_success) {
        tryCatch({
          # Use geom_blank to expand plot limits without modifying coordinate system
          p <- p + geom_blank(data = data.frame(x = c(final_xrange[1], expected_xmax), y = c(1, 1)),
                              aes(x = x, y = y))
          debug_cat(paste0("  v8: geom_blank expansion applied\n"))
        }, error = function(e3) {
          debug_cat(paste0("  v8: All expansion methods failed\n"))
        })
      }

      # Repair mapping after changes
      p <- func.repair.ggtree.mapping(p)
    } else {
      # v88: Plot is not buildable - don't add any expansion, just try to render
      debug_cat(paste0("  v8: WARNING - Plot cannot be built, skipping all expansion\n"))
    }
  }
  # ========================================================================

  # Default ellipse parameters if not set
  a <- 1
  b <- 1

  if (!exists("high_title_list")) {
    high_title_list <- ""
  }

  # v132: Initialize/refresh boudariestt before highlight code - required for ellipse sizing when heat_flag == TRUE
  if (!exists("boudariestt") || is.null(boudariestt)) {
    boudariestt <- tryCatch({
      func.find.plot.boundaries(p, debug_mode)
    }, error = function(e) {
      debug_cat(paste0("  v32: Error computing boudariestt (2nd block): ", e$message, "\n"))
      list(xmin = min(p$data$x, na.rm = TRUE), xmax = max(p$data$x, na.rm = TRUE))
    })
    debug_cat(paste0("  v32: boudariestt initialized (2nd block): xmin=", boudariestt$xmin, ", xmax=", boudariestt$xmax, "\n"))
  }

  # v139: Apply highlighting ONLY after heatmap (when heat_flag == TRUE)
  # The heatmap changes the coordinate system, so highlights need to be recalculated
  # When there's no heatmap, the first func_highlight call (line ~4835) is sufficient
  if (FLAG_BULK_DISPLAY == TRUE && heat_flag == TRUE) {
    x_adj_hi <- 0

    # Calculate tip length for ellipse sizing
    # If trimming is disabled (id_tip_trim_end is NA), use the max tip label length
    tip_length <- if (is.na(id_tip_trim_end)) {
      max(nchar(as.character(tree440$tip.label)), na.rm = TRUE)
    } else {
      id_tip_trim_end
    }

    # v146: Scale both 'a' (height) and 'b' (width) based on plot dimensions
    # When heatmap expands the x-range, we need to scale ellipses proportionally
    # to maintain the same visual relationship with the tree

    # Calculate the tree width (from 0 to max tree x)
    tree_width <- max(pr440_short_tips_TRY$data$x, na.rm = TRUE)

    # Calculate the total plot width (tree + heatmaps)
    total_width <- boudariestt$xmax - boudariestt$xmin

    # Scale factor: ratio of tree to total plot (clamped to reasonable range)
    # When no heatmap, this is close to 1; with heatmaps, this gets smaller
    tree_ratio <- tree_width / max(total_width, tree_width)
    tree_ratio <- min(max(tree_ratio, 0.1), 1.0)

    debug_cat(paste0("\n=== v146: ELLIPSE SCALING FOR HEATMAP ===\n"))
    debug_cat(paste0("  Tree width: ", round(tree_width, 2), "\n"))
    debug_cat(paste0("  Total plot width: ", round(total_width, 2), "\n"))
    debug_cat(paste0("  Tree ratio (tree/total): ", round(tree_ratio, 3), "\n"))

    # v146: Scale both height (a) and width (b) based on tree ratio
    # Use the same base calculation as non-heatmap case, then scale by tree_ratio
    base_a <- (tip_length * size_tip_text / 800 + man_adjust_elipse_a) * tree_ratio
    base_b <- (0.12 + man_adjust_elipse_b) * tree_ratio

    a <- base_a * adjust_height_ecliplse
    b <- base_b * adjust_width_eclipse

    debug_cat(paste0("  Base a (height): ", round(base_a, 4), "\n"))
    debug_cat(paste0("  Base b (width): ", round(base_b, 4), "\n"))
    debug_cat(paste0("  Final a (after user adjust): ", round(a, 4), "\n"))
    debug_cat(paste0("  Final b (after user adjust): ", round(b, 4), "\n"))
    debug_cat(paste0("=============================================\n"))

    # v139: Pass high_alpha_list for transparency
    .prof_section_start <- Sys.time()
    p <- func_highlight(
      p, how_many_hi, heat_flag, high_color_list, a, b, man_adjust_elipse,
      pr440_short_tips_TRY, boudariestt, debug_mode, high_offset, high_vertical_offset,
      high_alpha_list
    )
    cat(file=stderr(), sprintf("[PROF-TREE] func_highlight #2 (post-heatmap): %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))
  }

  if (length(b) == 0) {
    b <- 0.2
  }

  # v94: DEBUG - track layers after if(FALSE) block
  debug_cat(paste0("\n=== v94: AFTER if(FALSE) block ===\n"))
  debug_cat(paste0("  p layers: ", length(p$layers), "\n"))
  debug_cat(paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))
  debug_cat(paste0("================================\n"))

  # Add second legend if heatmap exists
  if (length(heat_map_title_list) > 0) {
    debug_cat(paste0("\n=== v95: Before func.make.second.legend ===\n"))
    debug_cat(paste0("  p layers: ", length(p$layers), "\n"))
    debug_cat(paste0("  heat_flag: ", heat_flag, "\n"))
    debug_cat(paste0("  FLAG_BULK_DISPLAY: ", FLAG_BULK_DISPLAY, "\n"))

    # v95: Debug boudariestt
    if (exists("boudariestt") && !is.null(boudariestt)) {
      debug_cat(paste0("  boudariestt$xmax: ", boudariestt$xmax, "\n"))
      debug_cat(paste0("  boudariestt$xmin: ", boudariestt$xmin, "\n"))
    } else {
      debug_cat(paste0("  WARNING: boudariestt is NULL or doesn't exist!\n"))
    }

    # v133: Get legend settings for highlight and bootstrap legends
    # v135: Use legend_settings parameter instead of values$legend_settings
    legend_settings_local <- legend_settings
    highlight_x_off <- if (!is.null(legend_settings_local$highlight_x_offset)) legend_settings_local$highlight_x_offset else 0
    highlight_y_off <- if (!is.null(legend_settings_local$highlight_y_offset)) legend_settings_local$highlight_y_offset else 0
    highlight_title_sz <- legend_settings_local$highlight_title_size  # NULL is ok, will use default
    highlight_text_sz <- legend_settings_local$highlight_text_size    # NULL is ok, will use default
    highlight_title_g <- if (!is.null(legend_settings_local$highlight_title_gap)) legend_settings_local$highlight_title_gap else 1
    highlight_label_g <- if (!is.null(legend_settings_local$highlight_label_gap)) legend_settings_local$highlight_label_gap else 0.5
    bootstrap_x_off <- if (!is.null(legend_settings_local$bootstrap_x_offset)) legend_settings_local$bootstrap_x_offset else 0
    bootstrap_y_off <- if (!is.null(legend_settings_local$bootstrap_y_offset)) legend_settings_local$bootstrap_y_offset else 0
    bootstrap_title_x_off <- if (!is.null(legend_settings_local$bootstrap_title_x_offset)) legend_settings_local$bootstrap_title_x_offset else 2  # v143
    bootstrap_title_sz <- legend_settings_local$bootstrap_title_size  # NULL is ok, will use default
    bootstrap_text_sz <- legend_settings_local$bootstrap_text_size    # NULL is ok, will use default
    bootstrap_title_g <- if (!is.null(legend_settings_local$bootstrap_title_gap)) legend_settings_local$bootstrap_title_gap else 2
    bootstrap_label_g <- if (!is.null(legend_settings_local$bootstrap_label_gap)) legend_settings_local$bootstrap_label_gap else 2
    # v138: Get show/hide settings for highlight and bootstrap legends
    show_highlight_leg <- if (!is.null(legend_settings_local$show_highlight)) legend_settings_local$show_highlight else TRUE
    show_bootstrap_leg <- if (!is.null(legend_settings_local$show_bootstrap)) legend_settings_local$show_bootstrap else TRUE

    debug_cat(paste0("  v33: highlight offsets - x:", highlight_x_off, ", y:", highlight_y_off, "\n"))
    debug_cat(paste0("  v33: bootstrap offsets - x:", bootstrap_x_off, ", y:", bootstrap_y_off, "\n"))
    debug_cat(paste0("  v38: show_highlight_legend:", show_highlight_leg, ", show_bootstrap_legend:", show_bootstrap_leg, "\n"))

    # v151: Calculate y_off_base to position legends on the RIGHT side
    # With coord_flip + scale_y_reverse: y values become horizontal positions
    # y = 1 is at visual LEFT, y = max_tips is at visual RIGHT
    # CRITICAL: Using y > max_tips EXPANDS the plot range and shrinks the tree!
    # Solution: Use y = max_tips (not max_tips + 5) to minimize expansion
    n_tips <- sum(p$data$isTip == TRUE, na.rm = TRUE)
    max_y_in_data <- max(p$data$y, na.rm = TRUE)
    y_off_base <- max(n_tips, max_y_in_data)  # Stay within existing data range
    debug_cat(paste0("  v51: y_off_base=", y_off_base, " (n_tips=", n_tips, ", max_y=", round(max_y_in_data, 2), ")\n"))

    # v95: Wrap in tryCatch to catch any errors
    .prof_section_start <- Sys.time()
    p <- tryCatch({
      result <- func.make.second.legend(
        p,
        FLAG_BULK_DISPLAY,
        how_many_hi,
        heat_flag,
        how_many_boxes,
        how_mant_rows,
        boudariestt,
        y_off_base,
        high_title_list,
        size_font_legend_title,
        high_label_list,
        size_font_legend_text,
        high_color_list,
        a,
        b,
        x_range_min,
        show_boot_flag,
        size_90,
        size_80,
        size_70,
        man_adjust_image_of_second_legend,
        man_multiply_second_legend,
        man_multiply_second_legend_text,
        man_multiply_elipse,
        man_space_second_legend,
        man_space_second_legend_multiplier,
        man_offset_for_highlight_legend_x,
        debug_mode,
        boot_values,
        man_offset_second_legend,
        width,
        bootstrap_label_size,
        # v133: New highlight and bootstrap legend settings
        highlight_x_offset = highlight_x_off,
        highlight_y_offset = highlight_y_off,
        highlight_title_size = highlight_title_sz,
        highlight_text_size = highlight_text_sz,
        highlight_title_gap = highlight_title_g,
        highlight_label_gap = highlight_label_g,
        bootstrap_x_offset = bootstrap_x_off,
        bootstrap_y_offset = bootstrap_y_off,
        bootstrap_title_x_offset = bootstrap_title_x_off,  # v143
        bootstrap_title_size_mult = bootstrap_title_sz,
        bootstrap_text_size_mult = bootstrap_text_sz,
        bootstrap_title_gap = bootstrap_title_g,
        bootstrap_label_gap = bootstrap_label_g,
        # v138: Show/hide legend controls
        show_highlight_legend = show_highlight_leg,
        show_bootstrap_legend = show_bootstrap_leg,
        # v145: Pass transparency list for legend ellipses
        high_alpha_list = high_alpha_list
      )
      debug_cat(paste0("  func.make.second.legend: SUCCESS\n"))
      result
    }, error = function(e) {
      debug_cat(paste0("  func.make.second.legend ERROR: ", e$message, "\n"))
      debug_cat(paste0("  Returning plot without second legend modifications\n"))
      p  # Return original plot on error
    })
    cat(file=stderr(), sprintf("[PROF-TREE] func.make.second.legend: %.3f sec\n", as.numeric(Sys.time() - .prof_section_start)))

    debug_cat(paste0("=== v95: After func.make.second.legend ===\n"))
    debug_cat(paste0("  p layers: ", length(p$layers), "\n"))
    debug_cat(paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))
  }

  # v94: DEBUG - before bootstrap triangles
  debug_cat(paste0("\n=== v94: Before bootstrap triangles ===\n"))
  debug_cat(paste0("  p layers: ", length(p$layers), "\n"))
  debug_cat(paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))

  # v96: FIX - Repair corrupted mapping BEFORE adding bootstrap triangles
  # gheatmap() corrupts the plot's @mapping attribute (changes it from ggplot2::mapping to data.frame)
  # This causes geom_nodepoint() calls below to crash silently
  # We must repair the mapping BEFORE adding any new layers
  debug_cat(paste0("\n=== v96: Repairing mapping before bootstrap triangles ===\n"))
  debug_cat(paste0("  Mapping class before repair: ", paste(class(p$mapping), collapse=", "), "\n"))
  p <- func.repair.ggtree.mapping(p, verbose = TRUE)
  debug_cat(paste0("  Mapping class after repair: ", paste(class(p$mapping), collapse=", "), "\n"))

  # v90: Add bootstrap triangles AFTER heatmap processing to avoid "missing x and y" error
  # When gheatmap transforms the plot data, any geom_nodepoint layers added beforehand
  # lose their x and y aesthetic mappings. By adding them here (after gheatmap),
  # the layers are created with the correct data structure.
  if (bootstrap_triangles_enabled && !is.null(bootstrap_triangles_params)) {
    debug_cat(paste0("\n=== v90: Adding bootstrap triangles after heatmap ===\n"))
    tryCatch({
      # v126: Reverted to v124 approach - separate geom_nodepoint for each bootstrap level
      # Using scale_size_manual with mapped size was conflicting with tree edge width scale
      p <- p +
        geom_nodepoint(
          position = position_nudge(x = bootstrap_triangles_params$man_boot_x_offset, y = 0),
          aes(subset = boot_val >= 0.9),
          size = bootstrap_triangles_params$size_90, shape = 24, fill = "grey36", colour = "grey20",
          show.legend = FALSE, alpha = 1/2
        ) +
        geom_nodepoint(
          position = position_nudge(x = bootstrap_triangles_params$man_boot_x_offset, y = 0),
          aes(subset = boot_val >= 0.8 & boot_val < 0.9),
          size = bootstrap_triangles_params$size_80, shape = 24, fill = "grey36", colour = "grey20",
          show.legend = FALSE, alpha = 1/2
        ) +
        geom_nodepoint(
          position = position_nudge(x = bootstrap_triangles_params$man_boot_x_offset, y = 0),
          aes(subset = boot_val >= 0.7 & boot_val < 0.8),
          size = bootstrap_triangles_params$size_70, shape = 24, fill = "grey36", colour = "grey20",
          show.legend = FALSE, alpha = 1/2
        )
      debug_cat(paste0("  Bootstrap triangles added successfully\n"))
    }, error = function(e) {
      debug_cat(paste0("  v6 ERROR adding bootstrap triangles: ", e$message, "\n"))
      debug_cat(paste0("  Continuing without bootstrap triangles\n"))
    })
    debug_cat(paste0("================================\n"))
  }

  # Add score information if requested
  if (flag_calc_scores_for_tree == TRUE) {
    # v53: print("SCORE")
    # v53: print(p_score_root)
    p_score_rounded <- round(p_score_root, digits = 5)
    ari_score_rounded <- round(ari_score, digits = 5)
    ps_score_text <- paste0("p score: ", p_score_rounded)
    ari_score_text <- paste0(" ari score: ", ari_score_rounded)
    p <- p + ggtitle(paste0(ps_score_text, ari_score_text))
  }

  # v69: Repair any corrupted mapping before saving/returning
  # This fixes the "@mapping must be <ggplot2::mapping>" error in newer ggplot2 versions
  p <- func.repair.ggtree.mapping(p)

  # S1.3-PERF: Removed intermediate call to func.move.tiplabels.to.front() here
  # Layer reordering is now done ONCE at the end in generate_plot() to avoid redundant calls.
  # v180 original comment: "Move tip labels to front so they render on top of highlight ellipses"

  # v72: FINAL DEBUG - verify heatmap layer exists and inspect plot state
  debug_cat(paste0("\n=== v72: FINAL PLOT STATE BEFORE GGSAVE ===\n"))
  debug_cat(paste0("  Number of layers: ", length(p$layers), "\n"))
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  debug_cat(paste0("  Layer types: ", paste(layer_types, collapse=", "), "\n"))

  # Check for GeomTile (heatmap)
  geomtile_idx <- which(layer_types == "GeomTile")
  if (length(geomtile_idx) > 0) {
    debug_cat(paste0("  GeomTile found at layers: ", paste(geomtile_idx, collapse=", "), "\n"))
    for (idx in geomtile_idx) {
      tryCatch({
        tile_layer <- p$layers[[idx]]
        tile_data <- tile_layer$data
        if (is.function(tile_data)) {
          tile_data <- tile_data(p$data)
        }
        if (!is.null(tile_data) && is.data.frame(tile_data)) {
          debug_cat(paste0("    Layer ", idx, " data rows: ", nrow(tile_data), "\n"))
          if ("x" %in% names(tile_data)) {
            debug_cat(paste0("    Layer ", idx, " x range: [", min(tile_data$x, na.rm=TRUE), ", ", max(tile_data$x, na.rm=TRUE), "]\n"))
          }
          if ("value" %in% names(tile_data)) {
            debug_cat(paste0("    Layer ", idx, " values: ", paste(unique(tile_data$value), collapse=", "), "\n"))
          }
        }
      }, error = function(e) {
        debug_cat(paste0("    ERROR: ", e$message, "\n"))
      })
    }
  } else {
    debug_cat(paste0("  WARNING: No GeomTile layers found! Heatmap may not be displayed.\n"))
  }

  # Check coordinate system
  if (!is.null(p$coordinates)) {
    debug_cat(paste0("  Coordinate system: ", class(p$coordinates)[1], "\n"))
    if (inherits(p$coordinates, "CoordCartesian")) {
      if (!is.null(p$coordinates$limits$x)) {
        debug_cat(paste0("  X limits: [", p$coordinates$limits$x[1], ", ", p$coordinates$limits$x[2], "]\n"))
      }
    }
  }

  # Try a final ggplot_build to get computed values
  tryCatch({
    final_built <- ggplot2::ggplot_build(p)
    debug_cat(paste0("  ggplot_build successful\n"))

    # Find tile layer in built data
    for (i in seq_along(final_built$data)) {
      if ("fill" %in% names(final_built$data[[i]]) && "width" %in% names(final_built$data[[i]])) {
        bd <- final_built$data[[i]]
        debug_cat(paste0("  Built layer ", i, ": ", nrow(bd), " rows, x=[",
                                  min(bd$x, na.rm=TRUE), ", ", max(bd$x, na.rm=TRUE),
                                  "], fill=", paste(unique(bd$fill), collapse=","), "\n"))
      }
    }

    # Check panel ranges
    if (!is.null(final_built$layout$panel_params)) {
      for (panel_idx in seq_along(final_built$layout$panel_params)) {
        pp <- final_built$layout$panel_params[[panel_idx]]
        if (!is.null(pp$x.range)) {
          debug_cat(paste0("  Panel ", panel_idx, " x.range: [", pp$x.range[1], ", ", pp$x.range[2], "]\n"))
        }
      }
    }
  }, error = function(e) {
    debug_cat(paste0("  ggplot_build error: ", e$message, "\n"))
  })
  debug_cat(paste0("============================================\n"))

  # v88: Save the plot with comprehensive error handling and multiple fallbacks
  save_success <- FALSE

  # v167: OPTION C - Use standard ggsave (legends are native ggplot layers now)
  # No gtable manipulation needed - legends are part of the plot object
  tryCatch({
    debug_cat(paste0("\n=== v167: Saving plot with native ggplot legends ===\n"))
    ggsave(out_file_path, plot = p, width = width, height = height, units = units_out, limitsize = FALSE)
    save_success <- TRUE
    debug_cat(paste0("=== v167: Plot saved successfully ===\n"))
    debug_cat(paste0("  LOOK FOR: 'v167 Test Legend' with red square alongside other legends\n"))
  }, error = function(e) {
    debug_cat(paste0("\n=== v167: GGSAVE ERROR ===\n"))
    debug_cat(paste0("  Primary error: ", e$message, "\n"))
  })

  # v88: Fallback 1 - repair mapping and try again
  if (!save_success) {
    debug_cat(paste0("  v8: Trying fallback 1 - repair mapping\n"))
    tryCatch({
      p_repaired <- func.repair.ggtree.mapping(p, verbose = TRUE)
      ggsave(out_file_path, plot = p_repaired, width = width, height = height,
             units = units_out, limitsize = FALSE)
      save_success <- TRUE
      debug_cat(paste0("  v8: Fallback 1 succeeded\n"))
    }, error = function(e2) {
      debug_cat(paste0("  v8: Fallback 1 failed: ", e2$message, "\n"))
    })
  }

  # v88: Fallback 2 - remove heatmap scale and try with defaults
  if (!save_success && heat_flag) {
    debug_cat(paste0("  v8: Trying fallback 2 - reset fill scale to defaults\n"))
    tryCatch({
      # Create fresh plot with default scale
      p_default <- p + scale_fill_discrete(na.value = "white")
      p_default <- func.repair.ggtree.mapping(p_default)
      ggsave(out_file_path, plot = p_default, width = width, height = height,
             units = units_out, limitsize = FALSE)
      save_success <- TRUE
      debug_cat(paste0("  v8: Fallback 2 succeeded (using default colors)\n"))
    }, error = function(e3) {
      debug_cat(paste0("  v8: Fallback 2 failed: ", e3$message, "\n"))
    })
  }

  # v88: Fallback 3 - use print() to render instead of ggsave
  if (!save_success) {
    debug_cat(paste0("  v8: Trying fallback 3 - direct PNG rendering\n"))
    tryCatch({
      # Determine file extension
      file_ext <- tolower(tools::file_ext(out_file_path))
      if (file_ext == "png") {
        png(out_file_path, width = width, height = height, units = units_out, res = 300)
      } else if (file_ext == "pdf") {
        pdf(out_file_path, width = width, height = height)
      } else {
        png(out_file_path, width = width, height = height, units = units_out, res = 300)
      }
      print(p)
      dev.off()
      save_success <- TRUE
      debug_cat(paste0("  v8: Fallback 3 succeeded\n"))
    }, error = function(e4) {
      debug_cat(paste0("  v8: Fallback 3 failed: ", e4$message, "\n"))
      tryCatch(dev.off(), error = function(x) {})  # Clean up device
    })
  }

  # v88: Final fallback - save tree without heatmap
  if (!save_success && heat_flag) {
    debug_cat(paste0("  v8: Trying fallback 4 - save tree without heatmap\n"))
    tryCatch({
      # Use the original tree plot (tt) without heatmap
      ggsave(out_file_path, plot = tt, width = width, height = height,
             units = units_out, limitsize = FALSE)
      save_success <- TRUE
      debug_cat(paste0("  v8: Fallback 4 succeeded (saved tree only, no heatmap)\n"))
      debug_cat(paste0("  v8: WARNING - Heatmap could not be rendered\n"))
    }, error = function(e5) {
      debug_cat(paste0("  v8: Fallback 4 failed: ", e5$message, "\n"))
      stop(e5)  # Re-throw if all fallbacks fail
    })
  }

  if (!save_success) {
    debug_cat(paste0("  v8: All save attempts failed\n"))
    stop("Could not save plot after multiple attempts")
  }

  cat(file=stderr(), sprintf("[PROF-TREE] === TOTAL func.make.plot.tree.heat.NEW: %.3f sec ===\n", as.numeric(Sys.time() - .prof_func_start)))

  # S2.0-PERF: Return both the plot and cache data for two-tier caching (Option 3A)
  # The cache data allows generate_plot() to store and reuse p_list_of_pairs
  # S2.9-PERF: Also return updated heatmap cache
  return(list(
    plot = p,
    cache_data = list(
      p_list_of_pairs = p_list_of_pairs,
      p_list_hash = current_p_list_hash,
      heatmap_cache = if (exists("updated_heatmap_cache")) updated_heatmap_cache else heatmap_cache
    )
  ))
}




##### part 5


#### part 5 a ui


# ============================================================================
# SHINY APPLICATION - PART 6 CORRECTIONS
# ============================================================================



# ============================================================================
# SHINY APPLICATION
# ============================================================================


# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Lineage Tree Plotter v159"),
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Upload Data", tabName = "data_upload", icon = icon("upload")),
      menuItem("Tree Display", tabName = "tree_display", icon = icon("tree")),
      menuItem("Classification", tabName = "classification", icon = icon("palette")),
      menuItem("Bootstrap Values", tabName = "bootstrap", icon = icon("percentage")),
      menuItem("Highlighting", tabName = "highlighting", icon = icon("highlighter")),
      menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
      menuItem("Legend", tabName = "legend", icon = icon("list")),
      menuItem("Extra", tabName = "extra", icon = icon("plus-circle")),  # v130: New tab for title, text, images
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("Configuration", tabName = "config", icon = icon("cogs"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    # Global Progress Bar
    tags$div(
      id = "global_progress",
      style = "position: fixed; top: 50px; left: 50%; transform: translateX(-50%); 
               width: 500px; max-width: 90%;
               z-index: 99999; 
               background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
               box-shadow: 0 4px 15px rgba(0,0,0,0.3);
               border-radius: 10px;
               padding: 20px 30px;
               display: none;
               animation: slideDown 0.3s ease-out;",
      tags$style("
        @keyframes slideDown {
          from { opacity: 0; transform: translateX(-50%) translateY(-20px); }
          to { opacity: 1; transform: translateX(-50%) translateY(0); }
        }
      "),
      uiOutput("progress_bar_ui")
    ),
    tabItems(
      # Data Upload Tab
      tabItem(
        tabName = "data_upload",
        fluidRow(
          box(
            title = "EZLineagePlotter - Development Version",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tags$div(style = "background: #fff3cd; padding: 15px; border-radius: 5px; border: 2px solid #856404;",
                     tags$h4(style = "color: #856404; margin: 0;", "Version S2.292dev"),
                     tags$p(style = "margin: 10px 0 0 0; color: #856404;",
                            tags$strong("New in S2.292dev:"),
                            tags$ul(
                              tags$li("Manual RGB/Hex color input for heatmap colors"),
                              tags$li("Font type selection for legend text"),
                              tags$li("Sync between color picker and hex text input")
                            ),
                            tags$strong("From S2.9 (stable):"),
                            tags$ul(
                              tags$li("Configurable two-stage CNV downsampling (import & render)"),
                              tags$li("Height Scale control for detailed RData heatmaps"),
                              tags$li("Heatmap caching to avoid regenerating unchanged heatmaps"),
                              tags$li("RData sample mapping column now properly saved to YAML")
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
        fluidRow(
          box(
            title = "Tree File",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            fileInput("tree_file", "Upload Newick Tree File",
                      accept = c(".nwk", ".tree", ".newick")),
            verbatimTextOutput("tree_summary")
          ),
          
          box(
            title = "Classification Data",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            fileInput("csv_file", "Upload CSV Classification File",
                      accept = c(".csv")),
            selectInput("id_column", "Select ID Column", choices = NULL),
            selectInput("individual_column", "Select Individual Column", 
                        choices = NULL),
            conditionalPanel(
              condition = "input.individual_column != null && input.individual_column != ''",
              # S2.0-PERF: Use selectizeInput for server-side rendering with large option lists
              selectizeInput("individual_value", "Select Individual",
                          choices = NULL,
                          selected = NULL,
                          options = list(placeholder = "Select an individual..."))
            ),
            
            # Checkbox to ignore individual filtering
            checkboxInput("use_all_data", "Use all data (ignore individual filtering)", 
                          value = FALSE),
            verbatimTextOutput("csv_summary"),
            
            # Add spacing
            br(),
            
            # Process button - only enabled when tree and CSV are loaded
            conditionalPanel(
              condition = "output.files_loaded == 'TRUE'",
              actionButton("process_data", "Process Data & Match IDs", 
                           class = "btn-success btn-lg", 
                           icon = icon("check-circle"),
                           style = "width: 100%;")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Matching Status",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            conditionalPanel(
              condition = "output.match_success == 'TRUE' && output.match_warning != 'TRUE'",
              div(style = "color: green;", 
                  tags$b("Success: All tree tips were matched to CSV data"))
            ),
            conditionalPanel(
              condition = "output.match_warning == 'TRUE'",
              div(style = "color: orange;", 
                  tags$b("Warning: Some tree tips couldn't be matched:"), 
                  verbatimTextOutput("unmatched_tips"))
            ),
            conditionalPanel(
              condition = "output.match_error == 'TRUE'",
              div(style = "color: red;", 
                  tags$b("Error: No tree tips could be matched to CSV data"))
            ),
            DTOutput("matching_table")
          )
        ),

        # S1.62: Optional YAML Settings Import
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
                " Load visual settings from a previously saved YAML configuration.",
                tags$br(),
                tags$small("This will apply classification colors, heatmap settings, bootstrap options, etc. to your current tree and CSV data.")
              )
            ),
            conditionalPanel(
              condition = "output.files_loaded == 'TRUE'",
              fileInput("yaml_config", "Choose YAML Configuration File",
                        accept = c(".yaml", ".yml")),
              verbatimTextOutput("yaml_import_status")
            ),
            conditionalPanel(
              condition = "output.files_loaded != 'TRUE'",
              tags$div(
                style = "color: #856404; padding: 10px; text-align: center;",
                icon("exclamation-triangle"),
                " Please upload tree and CSV files first before importing settings."
              )
            )
          )
        ),

        # S1.62dev: Optional RData CNV File Import
        fluidRow(
          box(
            title = "Import CNV Data from RData (Optional)",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            tags$div(
              style = "background: #d1ecf1; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
              tags$p(style = "margin: 0; color: #0c5460;",
                icon("dna"),
                " Load CNV (Copy Number Variation) data from an RData file.",
                tags$br(),
                tags$small("The RData file should contain 'results_CNV_tool_final' with CNV copy number data. This can be displayed as a heatmap aligned with the tree.")
              )
            ),
            conditionalPanel(
              condition = "output.files_loaded == 'TRUE'",
              fileInput("rdata_file", "Choose RData CNV File",
                        accept = c(".RData", ".rdata", ".Rdata")),
              numericInput("rdata_import_downsample",
                          "Import Downsample Factor",
                          value = 10, min = 1, max = 100, step = 1),
              tags$small(style = "color: #666; margin-top: -10px; display: block; margin-bottom: 10px;",
                        "Reduces genomic positions during import. 1 = keep all, 10 = keep every 10th position."),
              verbatimTextOutput("rdata_import_status")
            ),
            conditionalPanel(
              condition = "output.files_loaded != 'TRUE'",
              tags$div(
                style = "color: #0c5460; padding: 10px; text-align: center;",
                icon("exclamation-triangle"),
                " Please upload tree and CSV files first before importing CNV data."
              )
            )
          )
        )
      ),

      # Tree Display Tab
      tabItem(
        tabName = "tree_display",
        fluidRow(
          box(
            title = "Tree Appearance",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            checkboxInput("trim_tips", "Trim Tips", value = FALSE),
            conditionalPanel(
              condition = "input.trim_tips == true",
              sliderInput("tip_length", "Tip Length", min = 0.01, max = 0.2, value = 0.05, step = 0.01)
            ),
            sliderInput("edge_width", "Edge Width Multiplier", min = 0.5, max = 3, value = 1, step = 0.1),
            sliderInput("tip_font_size", "Tip Label Font Size", min = 1, max = 10, value = 3, step = 0.5),
            checkboxInput("display_node_numbers", "Display Node Numbers", value = FALSE),
            sliderInput("node_number_font_size", "Node Number Font Size", min = 1, max = 8, value = 3.5, step = 0.5),
            checkboxInput("use_pvalues", "Use P-values for Branch Width", value = TRUE),
            conditionalPanel(
              condition = "input.use_pvalues == true",
              sliderInput("fdr_perc", "FDR Percentage", min = 0.01, max = 0.25, value = 0.1, step = 0.01)
            ),
            checkboxInput("ladderize", "Ladderize Tree", value = FALSE)
          ),
          
          box(
            title = "Rotation Options",
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            checkboxInput("enable_rotation", "Enable Tree Rotation", value = FALSE),
            conditionalPanel(
              condition = "input.enable_rotation == true",
              radioButtons("rotation_type", "Rotation Type:",
                           choices = list(
                             "Primary Classification First" = "primary",
                             "Secondary Classification First" = "secondary",
                             "Manual Node Rotation" = "manual"
                           ), selected = "primary"),
              conditionalPanel(
                condition = "input.rotation_type == 'manual'",
                selectizeInput("nodes_to_rotate", "Select Nodes to Rotate", 
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "Select nodes to rotate")),
                checkboxInput("highlight_selected_nodes", "Highlight Selected Nodes on Tree", value = FALSE),
                tags$div(style = "margin-left: 20px; margin-top: -10px; margin-bottom: 10px;",
                         tags$small(style = "color: #666;", "Shows red circles on selected nodes to help you verify your selection")
                ),
                br(),
                actionButton("apply_manual_rotation", "Apply Rotation", icon = icon("rotate"), class = "btn-primary"),
                actionButton("clear_manual_rotation", "Clear Manual Configuration", icon = icon("trash"), class = "btn-warning")
              ),
              # Rotation Status Display
              uiOutput("rotation_status_box"),
              # Dynamic UI for rotation classes will be added via server logic
              uiOutput("rotation_classes_ui")
            )
          )
        ),
        
        fluidRow(
          box(
            title = tagList(
              "Tree Preview ",
              # v57: Static status indicator elements for shinyjs toggling (immediate UI updates)
              span(id = "status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("check-circle"), " Ready"
              ),
              span(id = "status_click_to_generate",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #e9ecef; color: #6c757d; font-size: 12px;",
                icon("hourglass-half"), " Click to generate"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            # v146: Changed to imageOutput with auto height for consistent aspect ratio
            imageOutput("tree_preview", height = "auto")
          )
        )
      ),
      
      # Classification Tab
      tabItem(
        tabName = "classification",
        fluidRow(
          box(
            title = "Classification Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            selectInput("classification_column", "Select Classification Column", choices = NULL, selected = character(0)),
            textInput("classification_title", "Legend Title", value = "Cell type"),
            selectInput("no_cluster_color", "No Cluster Color", choices = c("gray", "black", "white", "red"), selected = "gray"),
            actionButton("update_classification_preview", "Update Preview", icon = icon("eye"), class = "btn-info"),
            actionButton("add_classification", "Save Classification", icon = icon("save"), class = "btn-success"),
            actionButton("remove_classification", "Remove Selected", icon = icon("minus"), class = "btn-danger"),
            hr(),
            h4("Added Classifications:"),
            uiOutput("classifications_list_ui")
          ),
          
          box(
            title = "Classification Values",
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            uiOutput("classification_values_ui")
          )
        ),
        
        fluidRow(
          box(
            title = tagList(
              "Preview ",
              # v59: Static status indicators for immediate shinyjs updates (same as Tree Display)
              span(id = "class_status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "class_status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "class_status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("check-circle"), " Ready"
              ),
              span(id = "class_status_click_to_generate",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #e9ecef; color: #6c757d; font-size: 12px;",
                icon("hourglass-half"), " Click to generate"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            imageOutput("classification_preview", height = "auto")
          )
        )
      ),

      # Bootstrap Values Tab
      tabItem(
        tabName = "bootstrap",
        fluidRow(
          box(
            title = "Bootstrap Display",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            checkboxInput("show_bootstrap", "Show Bootstrap Values", value = TRUE),
            conditionalPanel(
              condition = "input.show_bootstrap == true",
              radioButtons("bootstrap_format", "Display Format:",
                           choices = list(
                             "Triangles" = "triangles",
                             "Raw Values" = "raw",
                             "Percentage" = "percentage",
                             "Color-coded Numbers" = "numbered_color",
                             "Color-coded Percentage" = "percentage_color"
                           ), selected = "triangles"),
              sliderInput("bootstrap_param", "Bootstrap Precision (decimal places)", 
                          min = 1, max = 5, value = 1, step = 1),
              # v130: Reduced default from 3 to 1.5 for smaller bootstrap legend by default
              sliderInput("bootstrap_label_size",
                          "Bootstrap Triangle Size:",
                          min = 0,
                          max = 10,
                          value = 1.5,
                          step = 0.5,
                          width = "100%"),
              # v114: Bootstrap position adjustment slider - finest precision with 0.001 step
              sliderInput("man_boot_x_offset",
                          "Bootstrap Position (higher/lower):",
                          min = -0.5,
                          max = 0.5,
                          value = 0,
                          step = 0.001,
                          width = "100%")
            )
          ),
          
          box(
            title = tagList(
              "Preview ",
              # v59: Static status indicators for immediate shinyjs updates
              span(id = "boot_status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "boot_status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "boot_status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("check-circle"), " Ready"
              ),
              span(id = "boot_status_click_to_generate",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #e9ecef; color: #6c757d; font-size: 12px;",
                icon("hourglass-half"), " Click to generate"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            imageOutput("bootstrap_preview", height = "auto")
          )
        )
      ),

      # Highlighting Tab
      # ============================================================================
      # HIGHLIGHTING TAB - COMPLETE REPLACEMENT
      # ============================================================================
      
      # Highlighting Tab UI (in the dashboardBody > tabItems section)
      tabItem(
        tabName = "highlighting",
        fluidRow(
          box(
            title = "Highlight Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            checkboxInput("enable_highlight", "Enable Highlighting", value = FALSE),
            conditionalPanel(
              condition = "input.enable_highlight == true",
              selectInput("highlight_column", "Select Column for Highlighting", 
                          choices = NULL, selected = NULL),
              selectizeInput("highlight_values", "Select Values to Highlight", 
                             choices = NULL, multiple = TRUE, 
                             options = list(maxOptions = 1000, 
                                            placeholder = "Select values to highlight")),
              textInput("highlight_title", "Highlight Legend Title", value = "Highlight"),
              hr(),
              h5("Global Positioning:"),
              # v54: Better precision control with numeric input + quick adjust buttons
              fluidRow(
                column(12,
                       tags$label("Vertical Offset (up/down)"),
                       fluidRow(
                         column(4,
                                actionButton("offset_down_big", "-1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")
                         ),
                         column(4,
                                actionButton("offset_down_small", "-0.1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")
                         ),
                         column(4,
                                actionButton("offset_down_tiny", "-0.01", class = "btn-sm btn-outline-secondary", style = "width: 100%;")
                         )
                       ),
                       numericInput("highlight_offset", NULL, 
                                    value = 0, min = -10, max = 10, step = 0.01, width = "100%"),
                       fluidRow(
                         column(4,
                                actionButton("offset_up_big", "+1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")
                         ),
                         column(4,
                                actionButton("offset_up_small", "+0.1", class = "btn-sm btn-outline-secondary", style = "width: 100%;")
                         ),
                         column(4,
                                actionButton("offset_up_tiny", "+0.01", class = "btn-sm btn-outline-secondary", style = "width: 100%;")
                         )
                       ),
                       tags$small(class = "text-muted", "Use buttons for quick adjustments, or type exact value")
                )
              ),
              br(),
              sliderInput("highlight_vertical_offset", "Horizontal Offset (left/right)", 
                          min = -1, max = 1, value = 0, step = 0.01),
              sliderInput("highlight_adjust_height", "Ellipse Height", 
                          min = 0.1, max = 3, value = 1, step = 0.05),
              sliderInput("highlight_adjust_width", "Ellipse Width", 
                          min = 0.1, max = 5, value = 1.5, step = 0.1),
              hr(),
              checkboxInput("auto_preview_highlight", "Auto-update preview on changes", value = FALSE),
              actionButton("update_highlight_preview", "Update Preview", 
                           icon = icon("eye"), class = "btn-info"),
              actionButton("save_highlight", "Save Highlight", 
                           icon = icon("save"), class = "btn-success"),
              actionButton("remove_highlight", "Remove Selected", 
                           icon = icon("minus"), class = "btn-danger"),
              hr(),
              h4("Saved Highlights:"),
              uiOutput("highlights_list_ui")
            )
          ),
          
          box(
            title = "Highlight Values Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            uiOutput("highlight_values_settings_ui")
          )
        ),
        
        fluidRow(
          box(
            title = tagList(
              "Tree Preview ",
              # v59: Static status indicators for immediate shinyjs updates
              span(id = "high_status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "high_status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "high_status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("check-circle"), " Ready"
              ),
              span(id = "high_status_click_to_generate",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #e9ecef; color: #6c757d; font-size: 12px;",
                icon("hourglass-half"), " Click to generate"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            imageOutput("highlight_preview", height = "auto")
          )
        )
      ),

      # Heatmap Tab
      tabItem(
        tabName = "heatmap",
        fluidRow(
          box(
            title = "Heatmap Configuration",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            checkboxInput("enable_heatmap", "Enable Heatmap Display", value = FALSE),
            conditionalPanel(
              condition = "input.enable_heatmap == true",
              hr(),
              
              # Global settings
              # v105: Removed global Distance from Tree slider - now per-heatmap only
              fluidRow(
                column(6,
                       sliderInput("heatmap_global_gap", "Gap Between Heatmaps",
                                   min = 0, max = 1, value = 0.05, step = 0.01)
                ),
                column(6,
                       tags$p(class = "text-muted", style = "padding-top: 20px;",
                              "Legend font size is controlled in the Legend tab. Per-heatmap settings are in each heatmap box below.")
                )
              ),
              
              hr(),
              
              # Add new heatmap button
              fluidRow(
                column(12,
                       actionButton("add_new_heatmap", "Add New Heatmap",
                                    icon = icon("plus"), class = "btn-success",
                                    style = "margin-bottom: 15px;"),
                       tags$span(class = "text-muted", style = "margin-left: 10px;",
                                 "Maximum 10 heatmaps allowed")  # v141: Increased from 6 to 10
                )
              ),
              
              # Container for heatmap cards
              uiOutput("heatmap_cards_ui"),
              
              hr(),
              
              # Action buttons
              fluidRow(
                column(12,
                       actionButton("apply_heatmaps", "Apply Heatmaps to Plot", 
                                    icon = icon("check"), class = "btn-primary btn-lg")
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = tagList(
              "Preview ",
              # v59: Static status indicators for immediate shinyjs updates
              span(id = "heat_status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "heat_status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "heat_status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("check-circle"), " Ready"
              ),
              span(id = "heat_status_click_to_generate",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #e9ecef; color: #6c757d; font-size: 12px;",
                icon("hourglass-half"), " Click to generate"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            imageOutput("heatmap_preview", height = "auto")
          )
        )
      ),

      # Legend Tab (v122)
      tabItem(
        tabName = "legend",
        fluidRow(
          # v122: Settings on left (4 columns), preview on right (8 columns)
          column(4,
            # Legend Position box
            box(
              title = NULL,  # v122: No side-by-side title
              status = "primary",
              solidHeader = FALSE,
              width = 12,
              tags$h4(icon("arrows-alt"), " Legend Position", style = "margin-top: 0;"),
              tags$p(class = "text-muted", "Where legends appear on the plot:"),
              selectInput("legend_position", NULL,  # v122: Title is above, not in label
                          choices = c("Right (default)" = "right",
                                      "Left" = "left",
                                      "Top" = "top",
                                      "Bottom" = "bottom",
                                      "None (hide all)" = "none"),
                          selected = "right")
            ),

            # v180: Legend Visibility box - simplified heatmap control
            box(
              title = NULL,
              status = "info",
              solidHeader = FALSE,
              width = 12,
              tags$h4(icon("eye"), " Legend Visibility", style = "margin-top: 0;"),
              tags$p(class = "text-muted", "Toggle which legends to show:"),
              checkboxInput("legend_show_classification", "Classification Legend", value = TRUE),
              checkboxInput("legend_show_highlight", "Highlight Legend", value = TRUE),
              checkboxInput("legend_show_bootstrap", "Bootstrap Legend", value = TRUE),
              checkboxInput("legend_show_pvalue", "P Value Legend", value = TRUE),
              checkboxInput("legend_show_heatmap", "Heatmap Legends", value = TRUE)  # v180: Simplified to single checkbox
            ),

            # Font Sizes box - S2.292dev: Added font type selection
            box(
              title = NULL,
              status = "warning",
              solidHeader = FALSE,
              width = 12,
              tags$h4(icon("text-height"), " Font Settings", style = "margin-top: 0;"),
              sliderInput("legend_title_size", "Legend Title Size",
                          min = 4, max = 48, value = 12, step = 1),
              sliderInput("legend_text_size", "Legend Text Size",
                          min = 2, max = 36, value = 10, step = 1),
              tags$hr(style = "margin: 10px 0;"),
              selectInput("legend_font_family", "Font Type",
                          choices = c("Sans-serif (default)" = "sans",
                                      "Serif (Times-like)" = "serif",
                                      "Monospace" = "mono",
                                      "Helvetica" = "Helvetica",
                                      "Arial" = "Arial",
                                      "Times New Roman" = "Times",
                                      "Courier" = "Courier",
                                      "Palatino" = "Palatino"),
                          selected = "sans"),
              tags$p(class = "text-muted", tags$small(
                "Note: Some fonts may not be available on all systems."
              ))
            ),

            # v180: Symbol & Spacing Settings box - enhanced with more controls
            box(
              title = NULL,
              status = "success",
              solidHeader = FALSE,
              width = 12,
              tags$h4(icon("square"), " Symbol & Spacing Settings", style = "margin-top: 0;"),
              sliderInput("legend_key_size", "Legend Key Size (symbols)",
                          min = 0.1, max = 5, value = 1, step = 0.1),
              tags$hr(style = "margin: 10px 0;"),
              tags$p(class = "text-muted", tags$small("Key dimensions:")),
              fluidRow(
                column(6,
                  sliderInput("legend_key_width", "Key Width",
                              min = 0.5, max = 3, value = 1, step = 0.1)
                ),
                column(6,
                  sliderInput("legend_key_height", "Key Height",
                              min = 0.5, max = 3, value = 1, step = 0.1)
                )
              ),
              tags$hr(style = "margin: 10px 0;"),
              tags$p(class = "text-muted", tags$small("Spacing between legends:")),
              sliderInput("legend_spacing", "Horizontal Spacing (top/bottom)",
                          min = 0.05, max = 3, value = 0.3, step = 0.05),
              sliderInput("legend_spacing_vertical", "Vertical Spacing (left/right)",
                          min = 0.1, max = 5, value = 1, step = 0.1),
              tags$hr(style = "margin: 10px 0;"),
              tags$p(class = "text-muted", tags$small("Spacing within each legend:")),
              sliderInput("legend_title_key_spacing", "Title to Keys Spacing",
                          min = 0, max = 2, value = 0.2, step = 0.05),  # v180: Renamed for clarity
              sliderInput("legend_key_spacing", "Between Keys Spacing",
                          min = 0, max = 2, value = 0.1, step = 0.05),  # v180: Actual key spacing
              tags$hr(style = "margin: 10px 0;"),
              tags$p(class = "text-muted", tags$small("Legend item order:")),
              checkboxInput("legend_reverse_order", "Reverse item order in legends", value = FALSE),
              tags$p(class = "text-muted", tags$small(
                "Note: Legend item order is determined by the data order. ",
                "Use this checkbox to reverse it."
              ))
            ),

            # v180: Legend Background box
            box(
              title = NULL,
              status = "warning",
              solidHeader = FALSE,
              width = 12,
              collapsible = TRUE,
              collapsed = TRUE,
              tags$h4(icon("palette"), " Legend Background", style = "margin-top: 0;"),
              tags$p(class = "text-muted", tags$small("Legend box appearance:")),
              fluidRow(
                column(6,
                  colourpicker::colourInput("legend_box_background", "Box Background",
                                            value = "transparent", showColour = "both",
                                            allowTransparent = TRUE)
                ),
                column(6,
                  sliderInput("legend_margin", "Legend Margin",
                              min = 0, max = 2, value = 0.2, step = 0.05)
                )
              )
            ),

            # v179: Removed Highlight Legend Settings and Bootstrap Legend Settings
            # These are no longer needed since legends now use native ggplot positioning

            # Apply button
            box(
              title = NULL,
              status = "primary",
              solidHeader = FALSE,
              width = 12,
              actionButton("apply_legend_settings", "Apply Legend Settings",
                           class = "btn-success btn-lg btn-block",
                           icon = icon("check"))
            )
          ),

          # v127: Plot preview on right (8 columns) - using imageOutput like other tabs
          # Added processing/ready status indicators like other tabs have
          column(8,
            box(
              title = tagList(
                "Legend Preview ",
                # v127: Static status indicators for immediate shinyjs updates
                span(id = "legend_status_waiting",
                  style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                  icon("clock"), " Waiting for data"
                ),
                span(id = "legend_status_processing",
                  style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                  icon("spinner", class = "fa-spin"), " Processing..."
                ),
                span(id = "legend_status_ready",
                  style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                  icon("check-circle"), " Ready"
                )
              ),
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              imageOutput("legend_preview", height = "auto")
            )
          )
        )
      ),

      # v130: Extra Tab - Page title, custom text annotations, and images
      # v141: Redesigned to add plot positioning and clarify text overlay behavior
      tabItem(
        tabName = "extra",
        fluidRow(
          # Plot Position Section (v141: NEW)
          box(
            title = tagList(
              icon("arrows-alt"), " Plot Position ",
              span(id = "extra_status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "extra_status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "extra_status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px;",
                icon("check"), " Ready"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            collapsible = FALSE,
            tags$p(class = "text-muted",
              "Move the entire plot (tree + heatmap + legend) on the page. ",
              "Text annotations are overlaid on top."
            ),
            fluidRow(
              column(6,
                sliderInput("plot_offset_x", "Horizontal Position:",
                           min = -5, max = 5, value = 0, step = 0.1,
                           post = " (left/right)")
              ),
              column(6,
                sliderInput("plot_offset_y", "Vertical Position:",
                           min = -10, max = 5, value = 0, step = 0.1,
                           post = " (down/up)")
              )
            ),
            fluidRow(
              column(12,
                actionButton("reset_plot_position", "Reset to Center",
                            class = "btn-secondary btn-sm", icon = icon("undo"))
              )
            ),
            hr(),
            # v146: Scale slider to zoom plot without distorting proportions
            tags$p(class = "text-muted",
              "Scale the entire plot up or down (zoom) without changing proportions."
            ),
            fluidRow(
              column(12,
                sliderInput("plot_scale_percent", "Plot Scale:",
                           min = 25, max = 200, value = 100, step = 5,
                           post = "%")
              )
            ),
            fluidRow(
              column(12,
                actionButton("reset_plot_scale", "Reset to 100%",
                            class = "btn-secondary btn-sm", icon = icon("undo"))
              )
            ),
            hr(),
            # v179: Tree stretch controls - make tree longer/shorter and wider/narrower
            tags$p(class = "text-muted",
              tags$strong("Tree Stretch:"), " Stretch the tree independently (may distort proportions)."
            ),
            fluidRow(
              column(6,
                sliderInput("tree_stretch_x", "Tree Width (vertical):",
                           min = 0.5, max = 3, value = 1, step = 0.1,
                           post = "x")
              ),
              column(6,
                sliderInput("tree_stretch_y", "Tree Length (horizontal):",
                           min = 0.5, max = 3, value = 1, step = 0.1,
                           post = "x")
              )
            ),
            fluidRow(
              column(12,
                actionButton("reset_tree_stretch", "Reset to 1x",
                            class = "btn-secondary btn-sm", icon = icon("undo"))
              )
            ),
            hr(),
            # v179: Background color control
            tags$p(class = "text-muted",
              tags$strong("Background:"), " Set the plot background color."
            ),
            fluidRow(
              column(6,
                colourpicker::colourInput("background_color", "Background Color:",
                                          value = "#FFFFFF", showColour = "both",
                                          allowTransparent = TRUE)
              ),
              column(6,
                actionButton("reset_background", "Reset to White",
                            class = "btn-secondary btn-sm", icon = icon("undo"),
                            style = "margin-top: 25px;")
              )
            )
          ),

          # Preview Box
          box(
            title = "Preview",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            # v146: Changed to auto height for consistent aspect ratio with other tabs
            imageOutput("extra_preview", height = "auto"),
            actionButton("extra_apply", "Apply to Plot", class = "btn-primary", style = "margin-top: 10px;")
          )
        ),

        # Page Title Section
        fluidRow(
          box(
            title = tagList(icon("heading"), " Page Title"),
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            checkboxInput("enable_page_title", "Enable Page Title", value = FALSE),
            conditionalPanel(
              condition = "input.enable_page_title == true",
              textInput("page_title_text", "Title Text:", value = ""),
              fluidRow(
                column(6,
                  numericInput("page_title_x", "X Position:", value = 0.5, min = 0, max = 1, step = 0.01)
                ),
                column(6,
                  numericInput("page_title_y", "Y Position:", value = 0.95, min = 0, max = 1, step = 0.01)
                )
              ),
              sliderInput("page_title_size", "Font Size:", min = 6, max = 72, value = 18, step = 1),
              colourpicker::colourInput("page_title_color", "Color:", value = "#000000"),
              fluidRow(
                column(6,
                  checkboxInput("page_title_bold", "Bold", value = TRUE)
                ),
                column(6,
                  checkboxInput("page_title_underline", "Underline", value = FALSE)
                )
              ),
              selectInput("page_title_hjust", "Horizontal Alignment:",
                          choices = c("Left" = "0", "Center" = "0.5", "Right" = "1"), selected = "0.5")
            )
          )
        ),

        # Custom Text Annotations Section (Overlay)
        fluidRow(
          box(
            title = tagList(icon("font"), " Text Overlay"),
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            tags$div(
              style = "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
              tags$p(style = "margin: 0; color: #155724;",
                icon("info-circle"),
                " Text is placed as an ", tags$b("overlay"), " on top of the plot. ",
                "Position the plot using the controls above, then add text wherever you want. ",
                "If text covers the plot, it will appear on top."
              )
            ),

            fluidRow(
              column(3,
                textInput("custom_text_content", "Text:", value = "")
              ),
              column(2,
                numericInput("custom_text_x", "X Position:", value = 0.5, min = 0, max = 1, step = 0.01)
              ),
              column(2,
                numericInput("custom_text_y", "Y Position:", value = 0.5, min = 0, max = 1, step = 0.01)
              ),
              column(2,
                sliderInput("custom_text_size", "Size:", min = 4, max = 36, value = 12, step = 1)
              ),
              column(2,
                colourpicker::colourInput("custom_text_color", "Color:", value = "#000000")
              )
            ),
            fluidRow(
              column(3,
                selectInput("custom_text_fontface", "Font Style:",
                            choices = c("Plain" = "plain", "Bold" = "bold", "Italic" = "italic", "Bold Italic" = "bold.italic"))
              ),
              column(3,
                selectInput("custom_text_hjust", "H. Align:",
                            choices = c("Left" = "0", "Center" = "0.5", "Right" = "1"), selected = "0.5")
              ),
              column(3,
                selectInput("custom_text_vjust", "V. Align:",
                            choices = c("Top" = "1", "Middle" = "0.5", "Bottom" = "0"), selected = "0.5")
              ),
              column(3,
                numericInput("custom_text_angle", "Rotation:", value = 0, min = -180, max = 180, step = 1)
              )
            ),
            fluidRow(
              column(4,
                actionButton("add_custom_text", "Add Text", class = "btn-success", icon = icon("plus"))
              ),
              column(4,
                actionButton("clear_custom_texts", "Clear All Texts", class = "btn-warning", icon = icon("trash"))
              )
            ),
            hr(),
            h5("Added Texts:"),
            uiOutput("custom_texts_list")
          )
        ),

        # Custom Images Section
        fluidRow(
          box(
            title = tagList(icon("image"), " Image Overlay"),
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            p("Add images as an overlay at specific positions on the page."),

            fluidRow(
              column(4,
                fileInput("custom_image_file", "Select Image:",
                          accept = c("image/png", "image/jpeg", "image/gif", "image/svg+xml"))
              ),
              column(2,
                numericInput("custom_image_x", "X Position:", value = 0.5, min = 0, max = 1, step = 0.01)
              ),
              column(2,
                numericInput("custom_image_y", "Y Position:", value = 0.5, min = 0, max = 1, step = 0.01)
              ),
              column(2,
                numericInput("custom_image_width", "Width:", value = 0.2, min = 0.01, max = 1, step = 0.01)
              ),
              column(2,
                numericInput("custom_image_height", "Height (0=auto):", value = 0, min = 0, max = 1, step = 0.01)
              )
            ),
            fluidRow(
              column(4,
                actionButton("add_custom_image", "Add Image", class = "btn-success", icon = icon("plus"))
              ),
              column(4,
                actionButton("clear_custom_images", "Clear All Images", class = "btn-warning", icon = icon("trash"))
              )
            ),
            hr(),
            h5("Added Images:"),
            uiOutput("custom_images_list")
          )
        )
      ),

      # Download Tab
      tabItem(
        tabName = "download",
        fluidRow(
          box(
            title = "Output Options",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            textInput("individual_name", "Sample/Individual Name", value = "Sample1"),
            selectInput("output_format", "File Format",
                        choices = c("pdf", "png", "tiff", "svg", "jpeg"), selected = "pdf"),
            # v138: Page orientation option
            selectInput("page_orientation", "Page Orientation",
                        choices = c("Landscape" = "landscape", "Portrait" = "portrait"),
                        selected = "landscape"),
            # v179: Option to preserve plot proportions when changing orientation
            checkboxInput("keep_proportions", "Keep plot proportions (don't stretch)", value = TRUE),
            # v125: Default to A4 landscape size (297 x 210 mm)
            numericInput("output_width", "Width", value = 29.7),
            numericInput("output_height", "Height", value = 21),
            selectInput("output_units", "Units", choices = c("cm", "mm", "in"), selected = "cm"),
            textInput("output_path", "Output Directory", value = "./"),
            checkboxInput("replace_name", "Use Custom File Name", value = FALSE),
            conditionalPanel(
              condition = "input.replace_name == true",
              textInput("custom_name", "Custom File Name", value = "tree_plot")
            ),
            textInput("prefix_text", "Optional Text at Beginning", value = "Tree__"),
            textInput("suffix_text", "Optional Text at End", value = "")
          ),
          
          box(
            title = tagList(
              "Final Preview ",
              # v142: Status indicators for Download tab
              span(id = "download_status_waiting",
                style = "display: inline-block; padding: 3px 10px; border-radius: 12px; background-color: #f8f9fa; color: #6c757d; font-size: 12px;",
                icon("clock"), " Waiting for data"
              ),
              span(id = "download_status_processing",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #6c757d; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("spinner", class = "fa-spin"), " Processing..."
              ),
              span(id = "download_status_ready",
                style = "display: none; padding: 3px 10px; border-radius: 12px; background-color: #28a745; color: #ffffff; font-size: 12px; font-weight: bold;",
                icon("check-circle"), " Ready"
              )
            ),
            status = "primary",
            solidHeader = TRUE,
            width = 8,
            # v125: Changed to imageOutput to match other previews
            imageOutput("final_preview", height = "auto"),
            downloadButton("download_plot", "Download Plot", class = "btn-primary"),
            downloadButton("download_yaml", "Download YAML Configuration", class = "btn-success")
          )
        )
      ),
      
      # Configuration Tab
      tabItem(
        tabName = "config",
        fluidRow(
          box(
            title = "YAML Configuration",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            tags$p("Current YAML configuration for the tree visualization:"),
            verbatimTextOutput("yaml_output"),
            downloadButton("download_yaml_config", "Download Configuration", class = "btn-success")
          )
        )
      )
    )
  )
)




# Define server logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)  # Increase to 100MB

  # S2.292dev: Version endpoint - allows checking version via URL
  # Access at: http://localhost:PORT/session/SESSION_ID/dataobj/version
  # Or use the simpler output$version_info below
  output$version_info <- renderText({
    paste0('{"version": "', VERSION, '", "app": "EZLineagePlotter"}')
  })

  # Reactive values to store data
  values <- reactiveValues(
    tree = NULL,
    tree_data = NULL,
    csv_data = NULL,
    id_match = NULL,
    classifications = list(),
    highlights = list(),
    heatmaps = list(),
    heatmap_configs = list(),  # v55: New multi-heatmap configuration list
    rotation_settings = list(
      primary = list(),
      secondary = list()
    ),
    rotation1_config = list(),  # Store primary rotation configuration
    rotation2_config = list(),  # Store secondary rotation configuration
    manual_rotation_config = list(),  # Store manual rotation configuration (node list)
    current_plot = NULL,
    plot_counter = 0,  # Counter to force reactive updates
    progress_message = "",  # Current progress message
    progress_visible = FALSE,  # Whether to show progress bar
    plot_generating = FALSE,  # Whether plot is currently being generated
    plot_ready = FALSE,  # Whether plot is ready to display
    yaml_import_status = NULL,  # S1.62dev: Status message for YAML import
    # S1.62dev: RData CNV data storage
    rdata_cnv_env = NULL,       # Raw environment from loaded RData
    rdata_cnv_matrix = NULL,    # Processed CNV matrix (rows=positions, cols=samples)
    rdata_import_status = NULL, # Status message for RData import
    rdata_sample_names = NULL,  # Sample names from RData (column names)
    rdata_auto_match = FALSE,   # S2.0: Whether RData samples auto-match tree tips
    rdata_mapping_column = NULL, # S2.0: User-selected CSV column for RData sample mapping
    # v121: Legend settings
    # S1.5: Added all missing defaults for proper legend styling
    # S2.292dev: Added font_family setting
    legend_settings = list(
      position = "right",
      show_classification = TRUE,
      show_highlight = TRUE,
      show_bootstrap = TRUE,
      show_heatmap = TRUE,
      show_pvalue = TRUE,
      title_size = 12,
      text_size = 10,
      font_family = "sans",
      key_size = 1,
      key_width = 1,
      key_height = 1,
      spacing = 0.3,
      spacing_vertical = 1,
      title_key_spacing = 0.2,
      key_spacing = 0.2,
      reverse_order = FALSE,
      box_background = "transparent",
      margin = 0.2
    ),
    # v130: Extra tab - page title, custom texts, and images
    page_title = list(
      enabled = FALSE,
      text = "",
      x = 0.5,
      y = 0.95,
      size = 18,
      color = "#000000",
      bold = TRUE,
      underline = FALSE,
      hjust = 0.5
    ),
    custom_texts = list(),  # List of custom text annotations
    custom_images = list(),  # List of custom images
    # v141: Plot position offsets for Extra tab
    plot_offset_x = 0,  # Horizontal offset (negative = left, positive = right)
    plot_offset_y = 0,  # Vertical offset (negative = down, positive = up)
    # v146: Plot scale percentage (zoom without distorting proportions)
    plot_scale_percent = 100,
    # v179: Tree stretch controls (may distort proportions)
    tree_stretch_x = 1,  # Horizontal stretch (tree length)
    tree_stretch_y = 1,  # Vertical stretch (tree width)
    # v179: Background color control
    background_color = "#FFFFFF",
    # S2.0-PERF: Two-tier caching for p_list_of_pairs (Option 3A)
    # These cache the expensive Fisher test calculations - only recalculate when
    # tree structure, classification, FDR, or simulate.p.value change.
    # Visual-only changes (legend sizes, colors, fonts) skip recalculation.
    # S2.7-PERF: Multi-entry cache - stores multiple classifications so switching
    # between previously computed classifications is instant (no recomputation).
    # The cache is keyed by hash, with LRU eviction when exceeding max entries.
    p_list_cache = list(),            # Named list: hash -> p_list_of_pairs
    p_list_cache_order = character(), # LRU order tracking (most recent last)
    p_list_cache_max_size = 10,       # Maximum number of cached classifications
    # Legacy single-entry cache variables (kept for compatibility, still updated)
    cached_p_list_hash = NULL,        # Hash of inputs that affect p-values
    cached_p_list_of_pairs = NULL,    # The cached p-value list
    cached_classification_column = NULL,  # The classification column that was used
    # S2.9-PERF: Heatmap caching - avoid regenerating unchanged heatmaps
    # The cache stores pre-computed tile_df and fill_colors for each heatmap index.
    # On Apply, we compare the current config hash with the cached hash.
    # If unchanged, we reuse the cached data instead of expensive regeneration.
    heatmap_cache = list(),           # Named list: heatmap_index -> list(tile_df, hash, ...)
    heatmap_cache_max_age = 5         # Clear cache entries older than N regenerations
  )

  classification_loading <- reactiveVal(FALSE)

  # S1.4-PERF: Plot trigger mechanism to batch multiple rapid updates
  # Instead of calling generate_plot() directly, observers increment this trigger.
  # A debounced observer watches the trigger and calls generate_plot() once.
  # This prevents redundant recalculations when multiple inputs change rapidly.
  plot_trigger <- reactiveVal(0)

  # S2.0-PERF: Increased debounce from 100ms to 500ms to better batch rapid changes
  # and prevent cascading plot regenerations that can freeze the UI
  plot_trigger_debounced <- debounce(reactive(plot_trigger()), 500)  # 500ms debounce

  # S2.0-PERF: Cooldown mechanism - track when last plot finished to prevent rapid re-triggering
  last_plot_time <- reactiveVal(0)
  PLOT_COOLDOWN_MS <- 1000  # 1 second cooldown after plot generation

  # S1.4-PERF: Helper function to request plot regeneration (use instead of direct generate_plot())
  request_plot_update <- function() {
    plot_trigger(plot_trigger() + 1)
  }

  # v107: Trigger for heatmap UI regeneration - only fires when we need to rebuild the UI
  # (add/remove/reorder heatmaps, or when CSV data changes)
  # This prevents the UI from rebuilding every time a slider or input changes
  heatmap_ui_trigger <- reactiveVal(0)

  # S2.8-FIX: Flag to inhibit color saving during UI rebuilds
  # When TRUE, color picker observers will skip saving to prevent default colors
  # from overwriting user's custom colors during heatmap add/remove operations
  inhibit_color_save <- reactiveVal(FALSE)

  # ==========================================================================
  # S1-PERF: DEBOUNCED REACTIVE INPUTS
  # ==========================================================================
  # These debounced versions prevent plot regeneration on every tiny slider movement.
  # The delay (300ms) waits until the user stops adjusting before triggering updates.
  # Use these instead of direct input$ values for sliders that affect plot rendering.

  DEBOUNCE_MS <- 300  # Debounce delay in milliseconds

  # Legend Tab sliders
  legend_title_size_d <- debounce(reactive(input$legend_title_size), DEBOUNCE_MS)
  legend_text_size_d <- debounce(reactive(input$legend_text_size), DEBOUNCE_MS)
  legend_key_size_d <- debounce(reactive(input$legend_key_size), DEBOUNCE_MS)
  legend_key_width_d <- debounce(reactive(input$legend_key_width), DEBOUNCE_MS)
  legend_key_height_d <- debounce(reactive(input$legend_key_height), DEBOUNCE_MS)
  legend_spacing_d <- debounce(reactive(input$legend_spacing), DEBOUNCE_MS)
  legend_spacing_vertical_d <- debounce(reactive(input$legend_spacing_vertical), DEBOUNCE_MS)
  legend_title_key_spacing_d <- debounce(reactive(input$legend_title_key_spacing), DEBOUNCE_MS)
  legend_key_spacing_d <- debounce(reactive(input$legend_key_spacing), DEBOUNCE_MS)
  legend_margin_d <- debounce(reactive(input$legend_margin), DEBOUNCE_MS)

  # Extra Tab sliders
  plot_offset_x_d <- debounce(reactive(input$plot_offset_x), DEBOUNCE_MS)
  plot_offset_y_d <- debounce(reactive(input$plot_offset_y), DEBOUNCE_MS)
  plot_scale_percent_d <- debounce(reactive(input$plot_scale_percent), DEBOUNCE_MS)
  tree_stretch_x_d <- debounce(reactive(input$tree_stretch_x), DEBOUNCE_MS)
  tree_stretch_y_d <- debounce(reactive(input$tree_stretch_y), DEBOUNCE_MS)
  page_title_size_d <- debounce(reactive(input$page_title_size), DEBOUNCE_MS)
  page_title_x_d <- debounce(reactive(input$page_title_x), DEBOUNCE_MS)
  page_title_y_d <- debounce(reactive(input$page_title_y), DEBOUNCE_MS)

  # Highlight sliders
  highlight_offset_d <- debounce(reactive(input$highlight_offset), DEBOUNCE_MS)
  highlight_vertical_offset_d <- debounce(reactive(input$highlight_vertical_offset), DEBOUNCE_MS)
  highlight_adjust_height_d <- debounce(reactive(input$highlight_adjust_height), DEBOUNCE_MS)
  highlight_adjust_width_d <- debounce(reactive(input$highlight_adjust_width), DEBOUNCE_MS)

  # Bootstrap sliders
  man_boot_x_offset_d <- debounce(reactive(input$man_boot_x_offset), DEBOUNCE_MS)
  bootstrap_label_size_d <- debounce(reactive(input$bootstrap_label_size), DEBOUNCE_MS)

  # Heatmap sliders
  heatmap_global_gap_d <- debounce(reactive(input$heatmap_global_gap), DEBOUNCE_MS)

  # Tree display sliders
  tip_font_size_d <- debounce(reactive(input$tip_font_size), DEBOUNCE_MS)
  edge_width_d <- debounce(reactive(input$edge_width), DEBOUNCE_MS)
  node_number_font_size_d <- debounce(reactive(input$node_number_font_size), DEBOUNCE_MS)
  tip_length_d <- debounce(reactive(input$tip_length), DEBOUNCE_MS)

  # S1.62dev: Background color (debounced to prevent rapid picker updates from causing lag)
  background_color_d <- debounce(reactive(input$background_color), DEBOUNCE_MS)

  # ==========================================================================
  # END DEBOUNCED REACTIVE INPUTS
  # ==========================================================================

  # ==========================================================================
  # S1.4-PERF: DEBOUNCED PLOT REGENERATION OBSERVER
  # ==========================================================================
  # This single observer handles all debounced plot update requests.
  # Instead of calling generate_plot() directly from many observers,
  # use request_plot_update() which increments plot_trigger.
  # This observer batches rapid updates and calls generate_plot() once.
  observeEvent(plot_trigger_debounced(), {
    req(values$plot_ready)  # Only regenerate if a plot exists
    req(plot_trigger() > 0)  # Ignore initial value

    # S2.0-PERF: Cooldown check - skip if a plot was generated too recently
    time_since_last <- as.numeric(Sys.time()) * 1000 - last_plot_time()
    if (time_since_last < PLOT_COOLDOWN_MS) {
      cat(file=stderr(), sprintf("[PERF] Debounced trigger skipped - cooldown active (%.0fms since last plot)\n", time_since_last))
      return()
    }

    cat(file=stderr(), "[PERF] Debounced plot trigger fired - regenerating plot\n")
    generate_plot()
  }, ignoreInit = TRUE)
  # ==========================================================================

  # v59: Helper functions to toggle status indicator via shinyjs (immediate UI updates)
  # Updated to toggle status indicators on ALL preview tabs (Tree Display, Classification, Bootstrap, Highlighting, Heatmap)
  show_status_waiting <- function() {
    # Tree Display tab
    shinyjs::show("status_waiting")
    shinyjs::hide("status_processing")
    shinyjs::hide("status_ready")
    shinyjs::hide("status_click_to_generate")
    # Classification tab
    shinyjs::show("class_status_waiting")
    shinyjs::hide("class_status_processing")
    shinyjs::hide("class_status_ready")
    shinyjs::hide("class_status_click_to_generate")
    # Bootstrap tab
    shinyjs::show("boot_status_waiting")
    shinyjs::hide("boot_status_processing")
    shinyjs::hide("boot_status_ready")
    shinyjs::hide("boot_status_click_to_generate")
    # Highlighting tab
    shinyjs::show("high_status_waiting")
    shinyjs::hide("high_status_processing")
    shinyjs::hide("high_status_ready")
    shinyjs::hide("high_status_click_to_generate")
    # Heatmap tab
    shinyjs::show("heat_status_waiting")
    shinyjs::hide("heat_status_processing")
    shinyjs::hide("heat_status_ready")
    shinyjs::hide("heat_status_click_to_generate")
    # v127: Legend tab
    shinyjs::show("legend_status_waiting")
    shinyjs::hide("legend_status_processing")
    shinyjs::hide("legend_status_ready")
    # v142: Download tab
    shinyjs::show("download_status_waiting")
    shinyjs::hide("download_status_processing")
    shinyjs::hide("download_status_ready")
    # Extra tab
    shinyjs::show("extra_status_waiting")
    shinyjs::hide("extra_status_processing")
    shinyjs::hide("extra_status_ready")
  }

  show_status_processing <- function() {
    # Tree Display tab
    shinyjs::hide("status_waiting")
    shinyjs::show("status_processing")
    shinyjs::hide("status_ready")
    shinyjs::hide("status_click_to_generate")
    # Classification tab
    shinyjs::hide("class_status_waiting")
    shinyjs::show("class_status_processing")
    shinyjs::hide("class_status_ready")
    shinyjs::hide("class_status_click_to_generate")
    # Bootstrap tab
    shinyjs::hide("boot_status_waiting")
    shinyjs::show("boot_status_processing")
    shinyjs::hide("boot_status_ready")
    shinyjs::hide("boot_status_click_to_generate")
    # Highlighting tab
    shinyjs::hide("high_status_waiting")
    shinyjs::show("high_status_processing")
    shinyjs::hide("high_status_ready")
    shinyjs::hide("high_status_click_to_generate")
    # Heatmap tab
    shinyjs::hide("heat_status_waiting")
    shinyjs::show("heat_status_processing")
    shinyjs::hide("heat_status_ready")
    shinyjs::hide("heat_status_click_to_generate")
    # v127: Legend tab
    shinyjs::hide("legend_status_waiting")
    shinyjs::show("legend_status_processing")
    shinyjs::hide("legend_status_ready")
    # v142: Download tab
    shinyjs::hide("download_status_waiting")
    shinyjs::show("download_status_processing")
    shinyjs::hide("download_status_ready")
    # Extra tab
    shinyjs::hide("extra_status_waiting")
    shinyjs::show("extra_status_processing")
    shinyjs::hide("extra_status_ready")
  }

  show_status_ready <- function() {
    # Tree Display tab
    shinyjs::hide("status_waiting")
    shinyjs::hide("status_processing")
    shinyjs::show("status_ready")
    shinyjs::hide("status_click_to_generate")
    # Classification tab
    shinyjs::hide("class_status_waiting")
    shinyjs::hide("class_status_processing")
    shinyjs::show("class_status_ready")
    shinyjs::hide("class_status_click_to_generate")
    # Bootstrap tab
    shinyjs::hide("boot_status_waiting")
    shinyjs::hide("boot_status_processing")
    shinyjs::show("boot_status_ready")
    shinyjs::hide("boot_status_click_to_generate")
    # Highlighting tab
    shinyjs::hide("high_status_waiting")
    shinyjs::hide("high_status_processing")
    shinyjs::show("high_status_ready")
    shinyjs::hide("high_status_click_to_generate")
    # Heatmap tab
    shinyjs::hide("heat_status_waiting")
    shinyjs::hide("heat_status_processing")
    shinyjs::show("heat_status_ready")
    shinyjs::hide("heat_status_click_to_generate")
    # v127: Legend tab
    shinyjs::hide("legend_status_waiting")
    shinyjs::hide("legend_status_processing")
    shinyjs::show("legend_status_ready")
    # v142: Download tab
    shinyjs::hide("download_status_waiting")
    shinyjs::hide("download_status_processing")
    shinyjs::show("download_status_ready")
    # Extra tab
    shinyjs::hide("extra_status_waiting")
    shinyjs::hide("extra_status_processing")
    shinyjs::show("extra_status_ready")
  }

  show_status_click_to_generate <- function() {
    # Tree Display tab
    shinyjs::hide("status_waiting")
    shinyjs::hide("status_processing")
    shinyjs::hide("status_ready")
    shinyjs::show("status_click_to_generate")
    # Classification tab
    shinyjs::hide("class_status_waiting")
    shinyjs::hide("class_status_processing")
    shinyjs::hide("class_status_ready")
    shinyjs::show("class_status_click_to_generate")
    # Bootstrap tab
    shinyjs::hide("boot_status_waiting")
    shinyjs::hide("boot_status_processing")
    shinyjs::hide("boot_status_ready")
    shinyjs::show("boot_status_click_to_generate")
    # Highlighting tab
    shinyjs::hide("high_status_waiting")
    shinyjs::hide("high_status_processing")
    shinyjs::hide("high_status_ready")
    shinyjs::show("high_status_click_to_generate")
    # Heatmap tab
    shinyjs::hide("heat_status_waiting")
    shinyjs::hide("heat_status_processing")
    shinyjs::hide("heat_status_ready")
    shinyjs::show("heat_status_click_to_generate")
    # v127: Legend tab (no click_to_generate status, just show waiting)
    shinyjs::show("legend_status_waiting")
    shinyjs::hide("legend_status_processing")
    shinyjs::hide("legend_status_ready")
    # Extra tab (no click_to_generate status, just show waiting)
    shinyjs::show("extra_status_waiting")
    shinyjs::hide("extra_status_processing")
    shinyjs::hide("extra_status_ready")
  }

  # Initialize YAML data structure immediately
  values$yaml_data <- list(
    "Individual general definitions" = list(
      Individual = "Sample1",
      "tree path" = list(NA),
      "mapping csv file" = NA,
      "individual column" = "",
      "out_file" = list(
        "base_path" = "./",
        "file_type" = "pdf",
        "optional text at beggining" = "Tree__",
        "optional text at end" = "",
        "replace name" = list(
          flag = "no",
          name = "tree_plot"
        )
      )
    ),
    "Mapping exl renaming titles" = list(
      "ID column" = ""
    ),
    "visual definitions" = list(
      "classification" = list(),
      "Bootstrap" = list(
        display = "yes",
        format = "triangles",
        param = "1"
      ),
      "rotation1" = list(
        display = "no",
        according = list()
      ),
      "rotation2" = list(
        display = "no",
        according = list()
      ),
      "manual_rotation" = list(
        display = "no",
        nodes = list()
      ),
      "trim tips" = list(
        display = "yes",
        length = 0.05
      ),
      "edge_width_multiplier" = list(
        size = 1
      ),
      "font_size" = list(
        tips = 3,
        legend_title = 30,
        legend_text = 20,
        legend_box = 15,
        heat_map_title = 25,
        heat_map_legend = 3.8
      )
    )
  )
  
  
  
  
  # Progress bar UI
  output$progress_bar_ui <- renderUI({
    if (values$progress_visible && values$progress_message != "") {
      tagList(
        tags$div(
          style = "display: flex; align-items: center; gap: 15px;",
          tags$div(
            class = "progress-spinner",
            style = "width: 35px; height: 35px; 
                     border: 4px solid rgba(255,255,255,0.3); 
                     border-top: 4px solid white; 
                     border-radius: 50%; 
                     animation: spin 0.8s linear infinite;",
            tags$style("@keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }")
          ),
          tags$div(
            style = "flex: 1;",
            tags$strong(values$progress_message, style = "font-size: 18px; color: white; font-weight: 600;")
          )
        )
      )
    }
  })
  
  # Show/hide progress bar
  observe({
    if (values$progress_visible) {
      shinyjs::show("global_progress")
    } else {
      shinyjs::hide("global_progress")
    }
  })
  
  # Output to control button visibility
  output$files_loaded <- renderText({
    as.character(!is.null(values$tree) && !is.null(values$csv_data))
  })
  outputOptions(output, "files_loaded", suspendWhenHidden = FALSE)
  
  # Define match status outputs (at top level of server function)
  output$match_success <- renderText({
    as.character(isTRUE(values$match_success))
  })
  
  output$match_warning <- renderText({
    as.character(isTRUE(values$match_warning))
  })
  
  output$match_error <- renderText({
    as.character(isTRUE(values$match_error))
  })
  
  outputOptions(output, "match_success", suspendWhenHidden = FALSE)
  outputOptions(output, "match_warning", suspendWhenHidden = FALSE)
  outputOptions(output, "match_error", suspendWhenHidden = FALSE)

  # S1.62dev: YAML import status output
  output$yaml_import_status <- renderText({
    if (is.null(values$yaml_import_status)) {
      ""
    } else {
      values$yaml_import_status
    }
  })

  # S1.62dev: RData CNV import status output
  output$rdata_import_status <- renderText({
    if (is.null(values$rdata_import_status)) {
      ""
    } else {
      values$rdata_import_status
    }
  })

  # S1.62dev: RData CNV file upload observer
  observeEvent(input$rdata_file, {
    cat(file=stderr(), "[RDATA-IMPORT] Observer triggered\n")
    req(input$rdata_file)
    cat(file=stderr(), "[RDATA-IMPORT] File received:", input$rdata_file$datapath, "\n")

    # Show progress indicator
    values$progress_message <- "🧬 Loading CNV data from RData..."
    values$progress_visible <- TRUE
    values$rdata_import_status <- "Processing RData file..."

    # Require tree and CSV to be loaded first
    if (is.null(values$tree) || is.null(values$csv_data)) {
      cat(file=stderr(), "[RDATA-IMPORT] ERROR: tree or csv_data is NULL\n")
      values$rdata_import_status <- "Error: Please upload tree and CSV files first."
      values$progress_visible <- FALSE
      showNotification("Please upload tree and CSV files first", type = "error")
      return()
    }

    # Extract CNV data from RData file
    # Use user-specified downsample factor (default 10) from UI input
    import_downsample <- if (!is.null(input$rdata_import_downsample)) input$rdata_import_downsample else 10
    cat(file=stderr(), paste0("[RDATA-IMPORT] Using import downsample factor: ", import_downsample, "\n"))
    result <- func.extract.cnv.from.rdata(input$rdata_file$datapath, downsample_factor = import_downsample)

    if (!is.null(result$error)) {
      cat(file=stderr(), "[RDATA-IMPORT] ERROR:", result$error, "\n")
      values$rdata_import_status <- paste("❌ Error:", result$error)
      values$progress_visible <- FALSE
      showNotification(result$error, type = "error")
      return()
    }

    # Store the extracted data
    values$rdata_cnv_matrix <- result$matrix
    values$rdata_sample_names <- result$sample_names

    # S1.62dev: Store chromosome information if available
    values$rdata_chr_info <- result$chr_info
    values$rdata_start_info <- result$start_info
    values$rdata_end_info <- result$end_info

    cat(file=stderr(), "[RDATA-IMPORT] Successfully loaded CNV data\n")
    cat(file=stderr(), "[RDATA-IMPORT] Matrix dimensions:", nrow(result$matrix), "x", ncol(result$matrix), "\n")
    cat(file=stderr(), "[RDATA-IMPORT] Sample names:", paste(head(result$sample_names, 5), collapse=", "),
        if(length(result$sample_names) > 5) "..." else "", "\n")

    # S2.0: Check auto-match - do RData sample names match tree tip labels?
    tree_tips <- values$tree$tip.label
    rdata_samples <- result$sample_names

    # Smart match function - normalizes strings for comparison
    smart_normalize <- function(x) {
      x <- gsub("[._-]", "", x)  # Remove common separators
      x <- gsub("^X", "", x)      # Remove leading X (R adds this to numeric column names)
      tolower(x)
    }

    # Try different matching strategies
    auto_match_count <- 0
    tree_tips_norm <- smart_normalize(tree_tips)
    rdata_samples_norm <- smart_normalize(rdata_samples)

    # Strategy 1: Exact match after normalization
    auto_match_count <- sum(tree_tips_norm %in% rdata_samples_norm)

    # Strategy 2: Substring match (RData sample contains tree tip or vice versa)
    if (auto_match_count < length(tree_tips) * 0.5) {
      for (tip in tree_tips_norm) {
        for (rs in rdata_samples_norm) {
          if (grepl(tip, rs, fixed = TRUE) || grepl(rs, tip, fixed = TRUE)) {
            auto_match_count <- auto_match_count + 1
            break
          }
        }
      }
    }

    auto_match_success <- auto_match_count >= length(tree_tips) * 0.5  # At least 50% match
    values$rdata_auto_match <- auto_match_success
    values$rdata_mapping_column <- NULL  # Reset mapping column

    cat(file=stderr(), paste0("[RDATA-IMPORT] Auto-match check: ", auto_match_count, "/", length(tree_tips),
                              " tips match (", round(auto_match_count/length(tree_tips)*100, 1), "%)\n"))
    cat(file=stderr(), paste0("[RDATA-IMPORT] Auto-match success: ", auto_match_success, "\n"))

    # S2.0: Build enhanced status message with sample names
    status_msg <- paste0(
      "✅ CNV data loaded successfully!\n",
      "Samples: ", result$n_samples, "\n",
      "Genomic positions: ", result$n_positions
    )
    if (!is.null(result$chr_info)) {
      unique_chrs <- unique(result$chr_info)
      status_msg <- paste0(status_msg, "\nChromosomes: ", paste(unique_chrs, collapse=", "))
    } else {
      status_msg <- paste0(status_msg, "\nChromosome info: Not available")
    }

    # Add sample names preview
    sample_preview <- paste(head(result$sample_names, 5), collapse=", ")
    if (length(result$sample_names) > 5) sample_preview <- paste0(sample_preview, ", ...")
    status_msg <- paste0(status_msg, "\n\n📋 Sample names (first 5):\n", sample_preview)

    # Add auto-match status
    if (auto_match_success) {
      status_msg <- paste0(status_msg, "\n\n✅ AUTO-MATCH SUCCESS!\nRData samples match tree tips - ready to use!")
    } else {
      status_msg <- paste0(status_msg, "\n\n⚠️ MANUAL MAPPING NEEDED\nRData samples don't match tree tips.\nSelect a mapping column in the Heatmap settings.")
    }

    values$rdata_import_status <- status_msg

    values$progress_visible <- FALSE
    showNotification(paste("CNV data loaded:", result$n_samples, "samples,", result$n_positions, "positions"),
                     type = "message")
  })

  # Bootstrap checkbox observer
  # Bootstrap checkbox observer
  observeEvent(input$show_bootstrap, {
    # v53: debug_cat("\n===observeEvent show_bootstrap FIRED===\n")
    # v53: debug_cat("New value:", input$show_bootstrap, "\n")
    # v53: cat(file=stderr(), "plot_ready:", values$plot_ready, "\n")
    
    req(values$plot_ready)
    
    # Force update YAML before generating
    update_yaml()
    
    # v53: cat(file=stderr(), "Bootstrap display in YAML:", 
    #     values$yaml_data$`visual definitions`$Bootstrap$display, "\n")
    
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Bootstrap format observer
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$bootstrap_format, {
    req(values$plot_ready, input$show_bootstrap == TRUE)
    request_plot_update()
  }, ignoreInit = TRUE)

  # Bootstrap param observer
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$bootstrap_param, {
    req(values$plot_ready, input$show_bootstrap == TRUE)
    request_plot_update()
  }, ignoreInit = TRUE)

  # Bootstrap label size observer
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$bootstrap_label_size, {
    req(values$plot_ready, input$show_bootstrap == TRUE)
    request_plot_update()
  }, ignoreInit = TRUE)

  # v112: Bootstrap position offset observer - makes slider immediately responsive
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(man_boot_x_offset_d(), {
    req(values$plot_ready, input$show_bootstrap == TRUE)
    request_plot_update()
  }, ignoreInit = TRUE)

  # S1.62dev: YAML import - Settings Only mode
  # Applies visual settings from YAML to the currently loaded tree/CSV data
  # Does NOT try to load files from paths (those paths are session-specific temp files)
  observeEvent(input$yaml_config, {
    cat(file=stderr(), "[YAML-IMPORT] Observer triggered\n")
    req(input$yaml_config)
    cat(file=stderr(), "[YAML-IMPORT] File received:", input$yaml_config$datapath, "\n")

    # S1.62dev: Show progress indicator immediately
    values$progress_message <- "📋 Importing YAML settings..."
    values$progress_visible <- TRUE
    values$yaml_import_status <- "Processing YAML configuration..."

    # S1.62dev: Require tree and CSV to be loaded first
    if (is.null(values$tree) || is.null(values$csv_data)) {
      cat(file=stderr(), "[YAML-IMPORT] ERROR: tree or csv_data is NULL\n")
      values$yaml_import_status <- "Error: Please upload tree and CSV files first before importing settings."
      values$progress_visible <- FALSE
      showNotification("Please upload tree and CSV files first", type = "error")
      return()
    }
    cat(file=stderr(), "[YAML-IMPORT] Tree and CSV are loaded, proceeding with YAML parse\n")

    # S1.62dev: Use withProgress for visible feedback during import
    # Show progress bar at top of page with notification style
    withProgress(message = 'Importing YAML settings...', value = 0, style = 'notification', {
      setProgress(value = 0.05, message = "Importing YAML settings...", detail = "Reading file...")
      yaml_data <- parse_yaml_config(input$yaml_config$datapath)
      if ("error" %in% names(yaml_data)) {
        cat(file=stderr(), "[YAML-IMPORT] ERROR parsing YAML:", yaml_data$error, "\n")
        values$yaml_import_status <- paste("Error parsing YAML:", yaml_data$error)
        values$progress_visible <- FALSE
        showNotification(yaml_data$error, type = "error")
        return()
      }
      cat(file=stderr(), "[YAML-IMPORT] YAML parsed successfully\n")
      setProgress(value = 0.2, detail = "File parsed successfully...")
    cat(file=stderr(), "[YAML-IMPORT] YAML keys:", paste(names(yaml_data), collapse=", "), "\n")

    # S1.62dev: Track what settings were imported
    imported_settings <- c()

    # Load individual name if available
    if (!is.null(yaml_data$`Individual general definitions`$Individual)) {
      updateTextInput(session, "individual_name", value = yaml_data$`Individual general definitions`$Individual)
      imported_settings <- c(imported_settings, "Individual name")
    }

    # S1.62dev: Load individual column if available
    if (!is.null(yaml_data$`Individual general definitions`$`individual column`)) {
      updateSelectInput(session, "individual_column", selected = yaml_data$`Individual general definitions`$`individual column`)
      imported_settings <- c(imported_settings, "Individual column")
    }

    # S1.62dev: Load ID column if available
    if (!is.null(yaml_data$`Mapping exl renaming titles`$`ID column`)) {
      updateSelectInput(session, "id_column", selected = yaml_data$`Mapping exl renaming titles`$`ID column`)
      imported_settings <- c(imported_settings, "ID column")
    }

    # S1.62dev: Skip file loading - use existing tree and CSV data
    # The YAML file paths point to temporary R session directories that no longer exist

    # Load display settings
    if (!is.null(yaml_data$`visual definitions`)) {
      
      setProgress(value = 0.3, detail = "Loading display settings...")

      # Load Bootstrap settings
      if (!is.null(yaml_data$`visual definitions`$Bootstrap)) {
        bootstrap_display <- yaml_data$`visual definitions`$Bootstrap$display
        if (func.check.bin.val.from.conf(bootstrap_display)) {
          updateCheckboxInput(session, "show_bootstrap", value = TRUE)
          if (!is.null(yaml_data$`visual definitions`$Bootstrap$format)) {
            updateRadioButtons(session, "bootstrap_format", 
                               selected = yaml_data$`visual definitions`$Bootstrap$format)
          }
          if (!is.null(yaml_data$`visual definitions`$Bootstrap$param)) {
            updateSliderInput(session, "bootstrap_param",
                              value = as.numeric(yaml_data$`visual definitions`$Bootstrap$param))
          }
          # S1.62dev: Import bootstrap label size
          if (!is.null(yaml_data$`visual definitions`$Bootstrap$label_size)) {
            updateSliderInput(session, "bootstrap_label_size",
                              value = as.numeric(yaml_data$`visual definitions`$Bootstrap$label_size))
            cat(file=stderr(), "[YAML-IMPORT] Imported bootstrap_label_size:",
                yaml_data$`visual definitions`$Bootstrap$label_size, "\n")
          }
        } else {
          updateCheckboxInput(session, "show_bootstrap", value = FALSE)
        }
      }
      
      # Load Rotation settings
      # S1.62dev: Fixed - populate rotation1_config/rotation2_config instead of rotation_settings
      if (!is.null(yaml_data$`visual definitions`$rotation1)) {
        rotation1_display <- yaml_data$`visual definitions`$rotation1$display
        if (func.check.bin.val.from.conf(rotation1_display)) {
          updateCheckboxInput(session, "enable_rotation", value = TRUE)
          updateRadioButtons(session, "rotation_type", selected = "primary")

          # Load rotation classes into rotation1_config
          if (!is.null(yaml_data$`visual definitions`$rotation1$according) &&
              length(yaml_data$`visual definitions`$rotation1$according) > 0) {
            config <- list()
            for (i in seq_along(yaml_data$`visual definitions`$rotation1$according)) {
              r <- yaml_data$`visual definitions`$rotation1$according[[i]]
              if (!is.null(r$col) && !is.null(r$val)) {
                config[[i]] <- list(col = r$col, val = r$val)
              }
            }
            values$rotation1_config <- config

            # Update UI to show the rotation groups
            if (length(config) > 0) {
              updateNumericInput(session, "rotation1_num_groups", value = length(config))
            }
            cat(file=stderr(), "[YAML-IMPORT] Imported rotation1 with", length(config), "groups\n")
          }
        }
      }

      if (!is.null(yaml_data$`visual definitions`$rotation2)) {
        rotation2_display <- yaml_data$`visual definitions`$rotation2$display
        if (func.check.bin.val.from.conf(rotation2_display)) {
          updateRadioButtons(session, "rotation_type", selected = "secondary")

          # Load rotation classes into rotation2_config
          if (!is.null(yaml_data$`visual definitions`$rotation2$according) &&
              length(yaml_data$`visual definitions`$rotation2$according) > 0) {
            config <- list()
            for (i in seq_along(yaml_data$`visual definitions`$rotation2$according)) {
              r <- yaml_data$`visual definitions`$rotation2$according[[i]]
              if (!is.null(r$col) && !is.null(r$val)) {
                config[[i]] <- list(col = r$col, val = r$val)
              }
            }
            values$rotation2_config <- config

            # Update UI to show the rotation groups
            if (length(config) > 0) {
              updateNumericInput(session, "rotation2_num_groups", value = length(config))
            }
            cat(file=stderr(), "[YAML-IMPORT] Imported rotation2 with", length(config), "groups\n")
          }
        }
      }

      # Load Manual rotation settings
      if (!is.null(yaml_data$`visual definitions`$manual_rotation)) {
        manual_rotation_display <- yaml_data$`visual definitions`$manual_rotation$display
        if (func.check.bin.val.from.conf(manual_rotation_display)) {
          updateCheckboxInput(session, "enable_rotation", value = TRUE)
          updateRadioButtons(session, "rotation_type", selected = "manual")

          # Load node list into manual_rotation_config
          if (!is.null(yaml_data$`visual definitions`$manual_rotation$nodes) &&
              length(yaml_data$`visual definitions`$manual_rotation$nodes) > 0) {
            # Convert list to numeric vector
            nodes <- unlist(yaml_data$`visual definitions`$manual_rotation$nodes)
            values$manual_rotation_config <- as.numeric(nodes)

            # Update UI to show the selected nodes
            updateSelectizeInput(session, "nodes_to_rotate", selected = as.character(nodes))
            cat(file=stderr(), "[YAML-IMPORT] Imported manual_rotation with", length(nodes), "nodes:",
                paste(nodes, collapse = ", "), "\n")
          }
        }
      }

      # Load Trim tips settings
      if (!is.null(yaml_data$`visual definitions`$`trim tips`)) {
        trim_display <- yaml_data$`visual definitions`$`trim tips`$display
        if (func.check.bin.val.from.conf(trim_display)) {
          updateCheckboxInput(session, "trim_tips", value = TRUE)
          if (!is.null(yaml_data$`visual definitions`$`trim tips`$length)) {
            updateSliderInput(session, "tip_length", 
                              value = as.numeric(yaml_data$`visual definitions`$`trim tips`$length))
          }
        } else {
          updateCheckboxInput(session, "trim_tips", value = FALSE)
        }
      }
      
      # Load Edge width multiplier
      if (!is.null(yaml_data$`visual definitions`$edge_width_multiplier)) {
        edge_width <- yaml_data$`visual definitions`$edge_width_multiplier$size
        if (!is.na(edge_width)) {
          updateSliderInput(session, "edge_width", value = as.numeric(edge_width))
        }
      }
      
      # Load Font sizes
      if (!is.null(yaml_data$`visual definitions`$font_size)) {
        if (!is.null(yaml_data$`visual definitions`$font_size$tips)) {
          updateSliderInput(session, "tip_font_size",
                            value = as.numeric(yaml_data$`visual definitions`$font_size$tips))
        }

        # S1.62dev: Load heatmap font size
        if (!is.null(yaml_data$`visual definitions`$font_size$heat_map_legend)) {
          updateSliderInput(session, "heatmap_font_size",
                            value = as.numeric(yaml_data$`visual definitions`$font_size$heat_map_legend))
        }
      }

      # S1.62dev: Load node numbers settings
      if (!is.null(yaml_data$`visual definitions`$node_numbers)) {
        node_numbers_display <- yaml_data$`visual definitions`$node_numbers$display
        if (func.check.bin.val.from.conf(node_numbers_display)) {
          updateCheckboxInput(session, "display_node_numbers", value = TRUE)
          if (!is.null(yaml_data$`visual definitions`$node_numbers$font_size)) {
            updateSliderInput(session, "node_number_font_size",
                              value = as.numeric(yaml_data$`visual definitions`$node_numbers$font_size))
          }
          cat(file=stderr(), "[YAML-IMPORT] Imported node_numbers display=yes, font_size=",
              yaml_data$`visual definitions`$node_numbers$font_size, "\n")
        } else {
          updateCheckboxInput(session, "display_node_numbers", value = FALSE)
          cat(file=stderr(), "[YAML-IMPORT] Imported node_numbers display=no\n")
        }
      }

      setProgress(value = 0.5, detail = "Loading classifications...")

      # Load Classification settings
      if (!is.null(yaml_data$`visual definitions`$classification)) {
        values$classifications <- list()
        
        for (i in seq_along(yaml_data$`visual definitions`$classification)) {
          class_def <- yaml_data$`visual definitions`$classification[[i]]
          class_item <- list()
          
          if (!is.null(class_def[[as.character(i)]])) {
            def <- class_def[[as.character(i)]]
            
            # Get FDR
            if (!is.null(def$FDR_perc)) {
              class_item$fdr <- as.numeric(def$FDR_perc)
              updateSliderInput(session, "fdr_perc", value = as.numeric(def$FDR_perc))
            } else {
              class_item$fdr <- 0.1
            }
            
            # Get title
            if (!is.null(def$title)) {
              class_item$title <- def$title
              updateTextInput(session, "classification_title", value = def$title)
            } else {
              class_item$title <- "Cell type"
            }
            
            # Get no cluster settings
            if (!is.null(def$non_cluster_title)) {
              class_item$no_cluster_title <- def$non_cluster_title
            } else {
              class_item$no_cluster_title <- "No cluster"
            }
            
            if (!is.null(def$non_cluster_color)) {
              class_item$no_cluster_color <- def$non_cluster_color
              updateSelectInput(session, "no_cluster_color", selected = def$non_cluster_color)
            } else {
              class_item$no_cluster_color <- "gray"
            }
            
            # Process classifications
            if (!is.null(def$according)) {
              class_item$classes <- list()
              
              for (j in seq_along(def$according)) {
                acc <- def$according[[j]]
                j_ch <- as.character(j)
                
                if (!is.null(acc[[j_ch]])) {
                  acc_item <- acc[[j_ch]]
                  
                  class_entry <- list(
                    column = acc_item$title1,
                    value = acc_item$value1,
                    display_name = acc_item$display_name,
                    color = acc_item$color
                  )
                  
                  class_item$classes[[j]] <- class_entry
                }
              }
            }
            
            # Process highlighting
            if (!is.null(def$highlight)) {
              highlight_display <- def$highlight$display
              
              if (func.check.bin.val.from.conf(highlight_display)) {
                updateCheckboxInput(session, "enable_highlight", value = TRUE)
                
                class_item$highlight <- list(enabled = TRUE, items = list())
                
                if (!is.null(def$highlight$according)) {
                  for (j in seq_along(def$highlight$according)) {
                    hi_acc <- def$highlight$according[[j]]
                    j_ch <- as.character(j)
                    
                    if (!is.null(hi_acc[[j_ch]])) {
                      hi_item <- hi_acc[[j_ch]]
                      
                      highlight_entry <- list(
                        column = hi_item$title1,
                        value = hi_item$value1,
                        display_name = hi_item$display_name,
                        color = hi_item$color
                      )
                      
                      # Add display title if available
                      if (!is.null(hi_item$display_title)) {
                        updateTextInput(session, "highlight_title", value = hi_item$display_title)
                      }
                      
                      class_item$highlight$items[[j]] <- highlight_entry
                    }
                  }
                }
              } else {
                class_item$highlight <- list(enabled = FALSE)
              }
            }
            
            # Process heatmaps
            if (!is.null(def$heatmap_display)) {
              class_item$heatmaps <- list()
              
              for (j in seq_along(def$heatmap_display)) {
                heat_def <- def$heatmap_display[[j]]
                
                if (func.check.bin.val.from.conf(heat_def$display)) {
                  updateCheckboxInput(session, "enable_heatmap", value = TRUE)
                  
                  heat_entry <- list(
                    title = heat_def$title,
                    is_discrete = func.check.bin.val.from.conf(heat_def$is_discrete)
                  )
                  
                  # Update UI for discrete/continuous
                  updateCheckboxInput(session, "heatmap_is_discrete", value = heat_entry$is_discrete)
                  
                  # Handle custom colors
                  if (!is.null(heat_def$man_define_colors) && 
                      func.check.bin.val.from.conf(heat_def$man_define_colors)) {
                    heat_entry$use_custom_colors <- TRUE
                    updateCheckboxInput(session, "use_custom_colors", value = TRUE)
                    
                    if (!is.null(heat_def$color_scale_option)) {
                      heat_entry$colors <- heat_def$color_scale_option
                    }
                  } else {
                    heat_entry$use_custom_colors <- FALSE
                  }
                  
                  # Handle columns
                  if (!is.null(heat_def$according)) {
                    heat_entry$columns <- list()
                    for (k in seq_along(heat_def$according)) {
                      heat_entry$columns[[k]] <- heat_def$according[[k]]
                    }
                  }
                  
                  class_item$heatmaps[[j]] <- heat_entry
                }
              }
            }
            
            values$classifications[[i]] <- class_item
          }
        }
      }
    }
    
    # Update output settings
    if (!is.null(yaml_data$`Individual general definitions`$out_file)) {
      if (!is.null(yaml_data$`Individual general definitions`$out_file$base_path)) {
        updateTextInput(session, "output_path", 
                        value = yaml_data$`Individual general definitions`$out_file$base_path)
      }
      
      if (!is.null(yaml_data$`Individual general definitions`$out_file$file_type)) {
        updateSelectInput(session, "output_format", 
                          selected = yaml_data$`Individual general definitions`$out_file$file_type)
      }
      
      if (!is.null(yaml_data$`Individual general definitions`$out_file$`optional text at beggining`)) {
        updateTextInput(session, "prefix_text", 
                        value = yaml_data$`Individual general definitions`$out_file$`optional text at beggining`)
      }
      
      if (!is.null(yaml_data$`Individual general definitions`$out_file$`optional text at end`)) {
        updateTextInput(session, "suffix_text", 
                        value = yaml_data$`Individual general definitions`$out_file$`optional text at end`)
      }
      
      if (!is.null(yaml_data$`Individual general definitions`$out_file$`replace name`)) {
        replace_name <- yaml_data$`Individual general definitions`$out_file$`replace name`

        if (!is.null(replace_name$flag) && func.check.bin.val.from.conf(replace_name$flag)) {
          updateCheckboxInput(session, "replace_name", value = TRUE)
          if (!is.null(replace_name$name)) {
            updateTextInput(session, "custom_name", value = replace_name$name)
          }
        }
      }

      # S1.62dev: Import page orientation and dimensions
      if (!is.null(yaml_data$`Individual general definitions`$out_file$page_orientation)) {
        updateSelectInput(session, "page_orientation",
                          selected = yaml_data$`Individual general definitions`$out_file$page_orientation)
        cat(file=stderr(), "[YAML-IMPORT] page_orientation:",
            yaml_data$`Individual general definitions`$out_file$page_orientation, "\n")
      }
      if (!is.null(yaml_data$`Individual general definitions`$out_file$output_width)) {
        updateNumericInput(session, "output_width",
                           value = as.numeric(yaml_data$`Individual general definitions`$out_file$output_width))
        cat(file=stderr(), "[YAML-IMPORT] output_width:",
            yaml_data$`Individual general definitions`$out_file$output_width, "\n")
      }
      if (!is.null(yaml_data$`Individual general definitions`$out_file$output_height)) {
        updateNumericInput(session, "output_height",
                           value = as.numeric(yaml_data$`Individual general definitions`$out_file$output_height))
        cat(file=stderr(), "[YAML-IMPORT] output_height:",
            yaml_data$`Individual general definitions`$out_file$output_height, "\n")
      }
      if (!is.null(yaml_data$`Individual general definitions`$out_file$output_units)) {
        updateSelectInput(session, "output_units",
                          selected = yaml_data$`Individual general definitions`$out_file$output_units)
        cat(file=stderr(), "[YAML-IMPORT] output_units:",
            yaml_data$`Individual general definitions`$out_file$output_units, "\n")
      }
      if (!is.null(yaml_data$`Individual general definitions`$out_file$keep_proportions)) {
        val <- func.check.bin.val.from.conf(yaml_data$`Individual general definitions`$out_file$keep_proportions)
        updateCheckboxInput(session, "keep_proportions", value = val)
        cat(file=stderr(), "[YAML-IMPORT] keep_proportions:", val, "\n")
      }

      imported_settings <- c(imported_settings, "Output settings")
    }

    setProgress(value = 0.6, detail = "Loading heatmaps...")

    # S1.62dev: Import heatmaps from new format
    if (!is.null(yaml_data$`visual definitions`$heatmaps) &&
        length(yaml_data$`visual definitions`$heatmaps) > 0) {
      cat(file=stderr(), "[YAML-IMPORT] Found heatmaps section with", length(yaml_data$`visual definitions`$heatmaps), "heatmaps\n")

      # Enable heatmap display checkbox
      updateCheckboxInput(session, "enable_heatmap", value = TRUE)

      values$heatmap_configs <- list()
      for (i in seq_along(yaml_data$`visual definitions`$heatmaps)) {
        h <- yaml_data$`visual definitions`$heatmaps[[i]]

        new_config <- list(
          columns = if (!is.null(h$columns)) unlist(h$columns) else character(0),
          title = if (!is.null(h$title)) h$title else paste0("Heatmap ", i),
          # S1.62dev: Import auto_type flag, default to FALSE for backward compatibility
          auto_type = if (!is.null(h$auto_type)) func.check.bin.val.from.conf(h$auto_type) else FALSE,
          type = if (!is.null(h$is_discrete) && func.check.bin.val.from.conf(h$is_discrete)) "discrete" else "continuous",
          show_colnames = if (!is.null(h$show_colnames)) func.check.bin.val.from.conf(h$show_colnames) else TRUE,
          colnames_angle = if (!is.null(h$colnames_angle)) as.numeric(h$colnames_angle) else 45,
          distance = if (!is.null(h$distance)) as.numeric(h$distance) else 0.02,
          height = if (!is.null(h$height)) as.numeric(h$height) else 0.8,
          # S1.62dev: Import row label settings from YAML instead of hardcoding
          show_row_labels = if (!is.null(h$show_row_labels)) func.check.bin.val.from.conf(h$show_row_labels) else FALSE,
          row_label_source = if (!is.null(h$row_label_source)) h$row_label_source else "colnames",
          row_label_font_size = if (!is.null(h$row_label_font_size)) as.numeric(h$row_label_font_size) else 2.5,
          custom_row_labels = if (!is.null(h$custom_row_labels)) h$custom_row_labels else "",
          row_label_offset = if (!is.null(h$row_label_offset)) as.numeric(h$row_label_offset) else 1.0,
          row_label_align = if (!is.null(h$row_label_align)) h$row_label_align else "left",
          label_mapping = if (!is.null(h$label_mapping)) h$label_mapping else list(),
          discrete_palette = if (!is.null(h$discrete_palette)) h$discrete_palette else "Set1",
          # S1.62dev: Import custom discrete colors from YAML instead of hardcoding
          custom_discrete = if (!is.null(h$custom_discrete)) func.check.bin.val.from.conf(h$custom_discrete) else FALSE,
          custom_colors = if (!is.null(h$custom_colors) && length(h$custom_colors) > 0) h$custom_colors else list(),
          cont_palette = if (!is.null(h$cont_palette)) h$cont_palette else "Blues",
          low_color = if (!is.null(h$low_color)) h$low_color else "#FFFFCC",
          high_color = if (!is.null(h$high_color)) h$high_color else "#006837",
          use_midpoint = if (!is.null(h$use_midpoint)) func.check.bin.val.from.conf(h$use_midpoint) else FALSE,
          mid_color = if (!is.null(h$mid_color)) h$mid_color else "#FFFF99",
          midpoint = if (!is.null(h$midpoint)) as.numeric(h$midpoint) else 0,
          # S1.61dev: Guide line settings
          show_guide_lines = if (!is.null(h$show_guide_lines)) func.check.bin.val.from.conf(h$show_guide_lines) else FALSE,
          guide_color1 = if (!is.null(h$guide_color1)) h$guide_color1 else "#CCCCCC",
          guide_color2 = if (!is.null(h$guide_color2)) h$guide_color2 else "#EEEEEE",
          guide_alpha = if (!is.null(h$guide_alpha)) as.numeric(h$guide_alpha) else 0.3,
          guide_width = if (!is.null(h$guide_width)) as.numeric(h$guide_width) else 0.5,
          guide_linetype = if (!is.null(h$guide_linetype)) h$guide_linetype else "solid",
          # S1.62dev: Import grid settings from YAML
          show_grid = if (!is.null(h$show_grid)) func.check.bin.val.from.conf(h$show_grid) else FALSE,
          grid_color = if (!is.null(h$grid_color)) h$grid_color else "#000000",
          grid_size = if (!is.null(h$grid_size)) as.numeric(h$grid_size) else 0.5,
          # S2.8: Import row line settings (horizontal lines)
          show_row_lines = if (!is.null(h$show_row_lines)) func.check.bin.val.from.conf(h$show_row_lines) else FALSE,
          row_line_color = if (!is.null(h$row_line_color)) h$row_line_color else "#000000",
          row_line_size = if (!is.null(h$row_line_size)) as.numeric(h$row_line_size) else 0.5,
          # S2.8: Import column line settings (vertical lines)
          show_col_lines = if (!is.null(h$show_col_lines)) func.check.bin.val.from.conf(h$show_col_lines) else FALSE,
          col_line_color = if (!is.null(h$col_line_color)) h$col_line_color else "#000000",
          col_line_size = if (!is.null(h$col_line_size)) as.numeric(h$col_line_size) else 0.5,
          # S2.8: Import RData heatmap settings
          data_source = if (!is.null(h$data_source)) h$data_source else "csv",
          rdata_mapping_column = if (!is.null(h$rdata_mapping_column)) h$rdata_mapping_column else "",
          # S2.8: Import WGD normalization settings
          cnv_wgd_norm = if (!is.null(h$cnv_wgd_norm)) func.check.bin.val.from.conf(h$cnv_wgd_norm) else FALSE,
          cnv_wgd_per_cell = if (!is.null(h$cnv_wgd_per_cell)) func.check.bin.val.from.conf(h$cnv_wgd_per_cell) else FALSE,
          cnv_wgd_column = if (!is.null(h$cnv_wgd_column)) h$cnv_wgd_column else "",
          # S2.8: Import display mode (basic or detailed)
          cnv_display_mode = if (!is.null(h$cnv_display_mode)) h$cnv_display_mode else "basic",
          # Render downsample factor (second stage, 0 = keep all)
          cnv_render_downsample = if (!is.null(h$cnv_render_downsample)) as.numeric(h$cnv_render_downsample) else 10,
          # Height scale for detailed mode
          cnv_height_scale = if (!is.null(h$cnv_height_scale)) as.numeric(h$cnv_height_scale) else 1.0
        )
        values$heatmap_configs <- c(values$heatmap_configs, list(new_config))

        # S2.292dev: Debug output for heatmap import
        cat(file=stderr(), sprintf("[YAML-IMPORT] Heatmap %d config:\n", i))
        cat(file=stderr(), sprintf("  data_source: %s\n", new_config$data_source))
        cat(file=stderr(), sprintf("  low_color: %s\n", new_config$low_color))
        cat(file=stderr(), sprintf("  high_color: %s\n", new_config$high_color))
        cat(file=stderr(), sprintf("  mid_color: %s\n", new_config$mid_color))
        cat(file=stderr(), sprintf("  use_midpoint: %s\n", new_config$use_midpoint))
        cat(file=stderr(), sprintf("  midpoint: %s\n", new_config$midpoint))
        if (!is.null(h$low_color)) cat(file=stderr(), sprintf("  YAML h$low_color: %s\n", h$low_color))
        if (!is.null(h$high_color)) cat(file=stderr(), sprintf("  YAML h$high_color: %s\n", h$high_color))
        if (!is.null(h$mid_color)) cat(file=stderr(), sprintf("  YAML h$mid_color: %s\n", h$mid_color))

        # S2.8: Update global rdata_mapping_column if this heatmap has one
        if (!is.null(new_config$data_source) && new_config$data_source == "rdata" &&
            !is.null(new_config$rdata_mapping_column) && new_config$rdata_mapping_column != "") {
          values$rdata_mapping_column <- new_config$rdata_mapping_column
          cat(file=stderr(), sprintf("[YAML-IMPORT] Set global rdata_mapping_column to: %s\n", new_config$rdata_mapping_column))
        }
      }
      # Trigger heatmap UI regeneration
      heatmap_ui_trigger(heatmap_ui_trigger() + 1)

      # Update RData mapping column dropdowns after UI regeneration
      # This ensures the saved values are displayed correctly
      # IMPORTANT: Must delay to allow Shiny to render the UI first
      if (!is.null(values$csv_data)) {
        csv_cols <- names(values$csv_data)
        shinyjs::delay(500, {
          for (cfg_idx in seq_along(values$heatmap_configs)) {
            cfg <- values$heatmap_configs[[cfg_idx]]
            if (!is.null(cfg$data_source) && cfg$data_source == "rdata" &&
                !is.null(cfg$rdata_mapping_column) && cfg$rdata_mapping_column != "") {
              updateSelectInput(session, paste0("heatmap_rdata_mapping_col_", cfg_idx),
                                choices = c("-- Select a column --" = "", csv_cols),
                                selected = cfg$rdata_mapping_column)
              cat(file=stderr(), sprintf("[YAML-IMPORT] Updated dropdown %d with saved mapping column: %s\n",
                                         cfg_idx, cfg$rdata_mapping_column))
            }
          }
        })
      }

      # S1.62dev: Also populate values$heatmaps for immediate plot rendering
      # This converts heatmap_configs to the format expected by the plot function
      heatmaps_for_plot <- lapply(values$heatmap_configs, function(cfg) {
        # S2.8: Only filter out heatmaps with no columns if they're NOT RData heatmaps
        # RData heatmaps don't use columns - they use the RData matrix directly
        if ((is.null(cfg$columns) || length(cfg$columns) == 0) &&
            (is.null(cfg$data_source) || cfg$data_source != "rdata")) {
          return(NULL)
        }
        list(
          title = if (!is.null(cfg$title)) cfg$title else "Heatmap",
          columns = cfg$columns,
          is_discrete = if (!is.null(cfg$type) && cfg$type == "discrete") TRUE else FALSE,
          # S1.62dev: Include show_colnames for column name visibility
          show_colnames = if (!is.null(cfg$show_colnames)) cfg$show_colnames else TRUE,
          colnames_angle = if (!is.null(cfg$colnames_angle)) cfg$colnames_angle else 45,
          discrete_palette = if (!is.null(cfg$discrete_palette)) cfg$discrete_palette else "Set1",
          # S1.62dev: color_scheme is what the plot function actually reads for discrete heatmaps
          color_scheme = if (!is.null(cfg$discrete_palette)) cfg$discrete_palette else "Set1",
          cont_palette = if (!is.null(cfg$cont_palette)) cfg$cont_palette else "Blues",
          custom_discrete = if (!is.null(cfg$custom_discrete)) cfg$custom_discrete else FALSE,
          man_define_colors = if (!is.null(cfg$custom_discrete)) cfg$custom_discrete else FALSE,
          custom_colors = if (!is.null(cfg$custom_colors)) cfg$custom_colors else list(),
          low_color = if (!is.null(cfg$low_color)) cfg$low_color else "#FFFFCC",
          high_color = if (!is.null(cfg$high_color)) cfg$high_color else "#006837",
          use_midpoint = if (!is.null(cfg$use_midpoint)) cfg$use_midpoint else FALSE,
          mid_color = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFF99",
          midpoint = if (!is.null(cfg$midpoint)) cfg$midpoint else 0,
          distance = if (!is.null(cfg$distance)) cfg$distance else 0.02,
          height = if (!is.null(cfg$height)) cfg$height else 0.8,
          row_height = if (!is.null(cfg$row_height)) cfg$row_height else 1,
          show_grid = if (!is.null(cfg$show_grid)) cfg$show_grid else FALSE,
          grid_color = if (!is.null(cfg$grid_color)) cfg$grid_color else "white",
          grid_size = if (!is.null(cfg$grid_size)) cfg$grid_size else 0.5,
          # S1.62dev: Row line settings
          show_row_lines = if (!is.null(cfg$show_row_lines)) cfg$show_row_lines else FALSE,
          row_line_color = if (!is.null(cfg$row_line_color)) cfg$row_line_color else "#000000",
          row_line_size = if (!is.null(cfg$row_line_size)) cfg$row_line_size else 0.5,
          # S2.8: Column line settings
          show_col_lines = if (!is.null(cfg$show_col_lines)) cfg$show_col_lines else FALSE,
          col_line_color = if (!is.null(cfg$col_line_color)) cfg$col_line_color else "#000000",
          col_line_size = if (!is.null(cfg$col_line_size)) cfg$col_line_size else 0.5,
          # S2.8: RData heatmap settings
          data_source = if (!is.null(cfg$data_source)) cfg$data_source else "csv",
          rdata_mapping_column = if (!is.null(cfg$rdata_mapping_column)) cfg$rdata_mapping_column else "",
          show_guides = if (!is.null(cfg$show_guides)) cfg$show_guides else FALSE,
          guide_color1 = if (!is.null(cfg$guide_color1)) cfg$guide_color1 else "#CCCCCC",
          guide_color2 = if (!is.null(cfg$guide_color2)) cfg$guide_color2 else "#EEEEEE",
          guide_alpha = if (!is.null(cfg$guide_alpha)) cfg$guide_alpha else 0.3,
          guide_width = if (!is.null(cfg$guide_width)) cfg$guide_width else 0.5,
          guide_linetype = if (!is.null(cfg$guide_linetype)) cfg$guide_linetype else "solid",
          show_row_labels = if (!is.null(cfg$show_row_labels)) cfg$show_row_labels else FALSE,
          row_label_source = if (!is.null(cfg$row_label_source)) cfg$row_label_source else "colnames",
          row_label_font_size = if (!is.null(cfg$row_label_font_size)) cfg$row_label_font_size else 2.5,
          # S1.62dev: Added missing row_label_offset and row_label_align
          row_label_offset = if (!is.null(cfg$row_label_offset)) cfg$row_label_offset else 1.0,
          row_label_align = if (!is.null(cfg$row_label_align)) cfg$row_label_align else "left",
          custom_row_labels = if (!is.null(cfg$custom_row_labels)) cfg$custom_row_labels else "",
          # S2.8: WGD normalization settings
          cnv_wgd_norm = if (!is.null(cfg$cnv_wgd_norm)) cfg$cnv_wgd_norm else FALSE,
          cnv_wgd_per_cell = if (!is.null(cfg$cnv_wgd_per_cell)) cfg$cnv_wgd_per_cell else FALSE,
          cnv_wgd_column = if (!is.null(cfg$cnv_wgd_column)) cfg$cnv_wgd_column else "",
          # S2.8: Display mode (basic or detailed)
          cnv_display_mode = if (!is.null(cfg$cnv_display_mode)) cfg$cnv_display_mode else "basic",
          # Render downsample factor (second stage, 0 = keep all)
          cnv_render_downsample = if (!is.null(cfg$cnv_render_downsample)) cfg$cnv_render_downsample else 10,
          # Height scale for detailed mode
          cnv_height_scale = if (!is.null(cfg$cnv_height_scale)) as.numeric(cfg$cnv_height_scale) else 1.0
        )
      })
      # Remove NULL entries (configs without columns)
      heatmaps_for_plot <- heatmaps_for_plot[!sapply(heatmaps_for_plot, is.null)]
      if (length(heatmaps_for_plot) > 0) {
        values$heatmaps <- heatmaps_for_plot
        cat(file=stderr(), "[YAML-IMPORT] Also populated values$heatmaps with", length(heatmaps_for_plot), "heatmaps for plot\n")
      }

      imported_settings <- c(imported_settings, "Heatmaps")
      cat(file=stderr(), "[YAML-IMPORT] Imported", length(values$heatmap_configs), "heatmap configs\n")
    }

    setProgress(value = 0.75, detail = "Loading highlights...")

    # S1.62dev: Import highlights from new format
    if (!is.null(yaml_data$`visual definitions`$highlights) &&
        length(yaml_data$`visual definitions`$highlights) > 0) {
      cat(file=stderr(), "[YAML-IMPORT] Found highlights section with", length(yaml_data$`visual definitions`$highlights), "highlights\n")

      # Enable highlight checkbox
      updateCheckboxInput(session, "enable_highlight", value = TRUE)

      values$highlights <- list()
      for (i in seq_along(yaml_data$`visual definitions`$highlights)) {
        h <- yaml_data$`visual definitions`$highlights[[i]]

        # Build items list
        items_list <- list()
        if (!is.null(h$items) && length(h$items) > 0) {
          for (j in seq_along(h$items)) {
            item <- h$items[[j]]
            items_list[[j]] <- list(
              column = if (!is.null(item$column)) item$column else "",
              value = if (!is.null(item$value)) item$value else "",
              display_name = if (!is.null(item$display_name)) item$display_name else "",
              color = if (!is.null(item$color)) item$color else "#FF0000",
              transparency = if (!is.null(item$transparency)) as.numeric(item$transparency) else 0.5
            )
          }
        }

        new_highlight <- list(
          enabled = if (!is.null(h$enabled)) func.check.bin.val.from.conf(h$enabled) else TRUE,
          title = if (!is.null(h$title)) h$title else "Highlight",
          column = if (!is.null(h$column)) h$column else "",
          offset = if (!is.null(h$offset)) as.numeric(h$offset) else 0,
          vertical_offset = if (!is.null(h$vertical_offset)) as.numeric(h$vertical_offset) else 0,
          adjust_height = if (!is.null(h$adjust_height)) as.numeric(h$adjust_height) else 1,
          adjust_width = if (!is.null(h$adjust_width)) as.numeric(h$adjust_width) else 1,
          items = items_list
        )
        values$highlights <- c(values$highlights, list(new_highlight))
      }

      # Set active highlight index
      values$active_highlight_index <- length(values$highlights)

      imported_settings <- c(imported_settings, "Highlights")
      cat(file=stderr(), "[YAML-IMPORT] Imported", length(values$highlights), "highlights\n")
    }

    # S1.62dev: Import legend settings from new format
    if (!is.null(yaml_data$`visual definitions`$legend)) {
      cat(file=stderr(), "[YAML-IMPORT] Found legend section\n")
      leg <- yaml_data$`visual definitions`$legend

      if (!is.null(leg$position)) {
        values$legend_settings$position <- leg$position
        updateSelectInput(session, "legend_position", selected = leg$position)
      }
      # S1.62dev: Import legend visibility settings and update UI checkboxes
      if (!is.null(leg$show_classification)) {
        val <- func.check.bin.val.from.conf(leg$show_classification)
        values$legend_settings$show_classification <- val
        updateCheckboxInput(session, "legend_show_classification", value = val)
        cat(file=stderr(), "[YAML-IMPORT] legend_show_classification:", val, "\n")
      }
      if (!is.null(leg$show_highlight)) {
        val <- func.check.bin.val.from.conf(leg$show_highlight)
        values$legend_settings$show_highlight <- val
        updateCheckboxInput(session, "legend_show_highlight", value = val)
        cat(file=stderr(), "[YAML-IMPORT] legend_show_highlight:", val, "\n")
      }
      if (!is.null(leg$show_bootstrap)) {
        val <- func.check.bin.val.from.conf(leg$show_bootstrap)
        values$legend_settings$show_bootstrap <- val
        updateCheckboxInput(session, "legend_show_bootstrap", value = val)
        cat(file=stderr(), "[YAML-IMPORT] legend_show_bootstrap:", val, "\n")
      }
      if (!is.null(leg$show_pvalue)) {
        val <- func.check.bin.val.from.conf(leg$show_pvalue)
        values$legend_settings$show_pvalue <- val
        updateCheckboxInput(session, "legend_show_pvalue", value = val)
        cat(file=stderr(), "[YAML-IMPORT] legend_show_pvalue:", val, "\n")
      }
      if (!is.null(leg$show_heatmap)) {
        val <- func.check.bin.val.from.conf(leg$show_heatmap)
        values$legend_settings$show_heatmap <- val
        updateCheckboxInput(session, "legend_show_heatmap", value = val)
        cat(file=stderr(), "[YAML-IMPORT] legend_show_heatmap:", val, "\n")
      }
      if (!is.null(leg$title_size)) {
        values$legend_settings$title_size <- as.numeric(leg$title_size)
        updateSliderInput(session, "legend_title_size", value = as.numeric(leg$title_size))
      }
      if (!is.null(leg$text_size)) {
        values$legend_settings$text_size <- as.numeric(leg$text_size)
        updateSliderInput(session, "legend_text_size", value = as.numeric(leg$text_size))
      }
      # S2.292dev: Import font_family
      if (!is.null(leg$font_family)) {
        values$legend_settings$font_family <- leg$font_family
        updateSelectInput(session, "legend_font_family", selected = leg$font_family)
        cat(file=stderr(), "[YAML-IMPORT] legend_font_family:", leg$font_family, "\n")
      }
      if (!is.null(leg$box_background)) {
        values$legend_settings$box_background <- leg$box_background
      }
      if (!is.null(leg$margin)) {
        values$legend_settings$margin <- as.numeric(leg$margin)
        updateSliderInput(session, "legend_margin", value = as.numeric(leg$margin))
      }
      imported_settings <- c(imported_settings, "Legend")
    }

    # S1.62dev: Import extra settings (plot position, scale, stretch, background)
    if (!is.null(yaml_data$`visual definitions`$extra_settings)) {
      cat(file=stderr(), "[YAML-IMPORT] Found extra_settings section\n")
      extra <- yaml_data$`visual definitions`$extra_settings

      if (!is.null(extra$plot_offset_x)) {
        values$plot_offset_x <- as.numeric(extra$plot_offset_x)
        updateSliderInput(session, "plot_offset_x", value = as.numeric(extra$plot_offset_x))
        cat(file=stderr(), "[YAML-IMPORT] plot_offset_x:", extra$plot_offset_x, "\n")
      }
      if (!is.null(extra$plot_offset_y)) {
        values$plot_offset_y <- as.numeric(extra$plot_offset_y)
        updateSliderInput(session, "plot_offset_y", value = as.numeric(extra$plot_offset_y))
        cat(file=stderr(), "[YAML-IMPORT] plot_offset_y:", extra$plot_offset_y, "\n")
      }
      if (!is.null(extra$plot_scale)) {
        updateSliderInput(session, "plot_scale", value = as.numeric(extra$plot_scale))
        cat(file=stderr(), "[YAML-IMPORT] plot_scale:", extra$plot_scale, "\n")
      }
      if (!is.null(extra$tree_stretch_x)) {
        values$tree_stretch_x <- as.numeric(extra$tree_stretch_x)
        updateSliderInput(session, "tree_stretch_x", value = as.numeric(extra$tree_stretch_x))
        cat(file=stderr(), "[YAML-IMPORT] tree_stretch_x:", extra$tree_stretch_x, "\n")
      }
      if (!is.null(extra$tree_stretch_y)) {
        values$tree_stretch_y <- as.numeric(extra$tree_stretch_y)
        updateSliderInput(session, "tree_stretch_y", value = as.numeric(extra$tree_stretch_y))
        cat(file=stderr(), "[YAML-IMPORT] tree_stretch_y:", extra$tree_stretch_y, "\n")
      }
      if (!is.null(extra$background_color)) {
        colourpicker::updateColourInput(session, "background_color", value = extra$background_color)
        cat(file=stderr(), "[YAML-IMPORT] background_color:", extra$background_color, "\n")
      }
      imported_settings <- c(imported_settings, "Extra settings")
    }

    # S1.62dev: Track visual definitions imported
    if (!is.null(yaml_data$`visual definitions`)) {
      if (!is.null(yaml_data$`visual definitions`$Bootstrap)) {
        imported_settings <- c(imported_settings, "Bootstrap")
      }
      if (!is.null(yaml_data$`visual definitions`$rotation1) || !is.null(yaml_data$`visual definitions`$rotation2)) {
        imported_settings <- c(imported_settings, "Rotation")
      }
      if (!is.null(yaml_data$`visual definitions`$`trim tips`)) {
        imported_settings <- c(imported_settings, "Trim tips")
      }
      if (!is.null(yaml_data$`visual definitions`$edge_width_multiplier)) {
        imported_settings <- c(imported_settings, "Edge width")
      }
      if (!is.null(yaml_data$`visual definitions`$font_size)) {
        imported_settings <- c(imported_settings, "Font size")
      }
      if (!is.null(yaml_data$`visual definitions`$classification) &&
          length(yaml_data$`visual definitions`$classification) > 0) {
        imported_settings <- c(imported_settings, "Classification")
      }
    }

    # Force update of UI elements based on currently loaded data
    updateSelectInput(session, "classification_column", choices = names(values$csv_data), selected = character(0))
    updateSelectInput(session, "highlight_column", choices = names(values$csv_data), selected = character(0))

    cat(file=stderr(), "[YAML-IMPORT] imported_settings:", paste(imported_settings, collapse=", "), "\n")

    # Generate plot with imported settings
    request_plot_update()

    # S1.62dev: Set status message and finalize
    setProgress(value = 0.95, detail = "Finalizing...")

    if (length(imported_settings) > 0) {
      values$yaml_import_status <- paste0(
        "✅ Settings imported successfully!\n",
        "Applied: ", paste(imported_settings, collapse = ", ")
      )
      cat(file=stderr(), "[YAML-IMPORT] SUCCESS - Settings imported:", paste(imported_settings, collapse=", "), "\n")
      showNotification("YAML settings imported successfully", type = "message")
    } else {
      values$yaml_import_status <- "⚠️ YAML loaded but no applicable settings found."
      cat(file=stderr(), "[YAML-IMPORT] WARNING - No applicable settings found in YAML\n")
      showNotification("YAML loaded but no applicable settings found", type = "warning")
    }

    # S1.62dev: Hide progress indicator when done
    values$progress_visible <- FALSE
    values$progress_message <- ""

    }) # End withProgress
  })
  
  # When CSV file is uploaded
  observeEvent(input$csv_file, {
    cat(file=stderr(), "[DEBUG-CSV] CSV file observer triggered\n")
    req(input$csv_file)

    csv_file <- input$csv_file
    cat(file=stderr(), sprintf("[DEBUG-CSV] File path: %s\n", csv_file$datapath))

    # Show progress
    values$progress_message <- "[FILE]  Loading CSV file..."
    values$progress_visible <- TRUE
    
    # Read CSV file
    tryCatch({
      # S2.0-PERF: Use data.table::fread() for fast CSV reading
      # We call it with :: to avoid loading the library and masking dplyr functions
      start_time <- Sys.time()
      cat(file=stderr(), "[DEBUG-CSV] About to call fread()\n")

      # Read with fread - much faster than read.csv or readr
      # data.table = FALSE returns a data.frame instead of data.table
      # Note: We avoid check.names parameter for compatibility, then manually fix names
      csv_data_raw <- data.table::fread(csv_file$datapath, data.table = FALSE)
      cat(file=stderr(), "[DEBUG-CSV] fread() completed\n")

      # Manually make names syntactically valid (like read.csv does with check.names=TRUE)
      # This ensures column names work properly for classification coloring
      names(csv_data_raw) <- make.names(names(csv_data_raw), unique = TRUE)
      cat(file=stderr(), "[DEBUG-CSV] Column names validated\n")

      read_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat(file=stderr(), sprintf("[PERF] CSV read with fread(): %.3f sec (%d rows x %d cols)\n",
          read_time, nrow(csv_data_raw), ncol(csv_data_raw)))

      # S2.0-PERF: Filter out columns with empty/auto-generated names
      # Patterns: "V1", "V2" (fread), "...XXXX" (readr), "X", "X.1", "X.2" (base R), empty names
      col_names <- names(csv_data_raw)
      valid_cols <- !grepl("^V\\d+$", col_names) &           # fread pattern: V1, V2, V3
                    !grepl("^\\.\\.\\.", col_names) &        # readr pattern: ...15372
                    !grepl("^X(\\.\\d+)?$", col_names) &     # base R pattern: X, X.1, X.2
                    col_names != "" & !is.na(col_names)

      removed_count <- sum(!valid_cols)
      if (removed_count > 0) {
        cat(file=stderr(), sprintf("[PERF] Removed %d empty/auto-named columns (keeping %d)\n",
            removed_count, sum(valid_cols)))
        csv_data <- csv_data_raw[, valid_cols, drop = FALSE]
      } else {
        csv_data <- csv_data_raw
      }

      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat(file=stderr(), sprintf("[PERF] CSV load complete: %.3f sec total\n", total_time))

      values$csv_data <- csv_data
      # v107: Trigger heatmap UI regeneration when CSV data changes (new column choices)
      heatmap_ui_trigger(heatmap_ui_trigger() + 1)

      # Update ID column dropdown
      updateSelectInput(session, "id_column", choices = names(csv_data))
      
      # Update individual_column column dropdown
      updateSelectInput(session, "individual_column", choices = names(csv_data))
      
      # Update classification column dropdown
      updateSelectInput(session, "classification_column", 
                        choices = c("-- Select Column --" = "", names(csv_data)), 
                        selected = "")   
      
      # Update highlight column dropdown
      # Update highlight column dropdown  
      updateSelectInput(session, "highlight_column", 
                        choices = c("-- Select Column --" = "", names(csv_data)),
                        selected = "")
      
      # v55: heatmap_columns removed - now using individual column selects per heatmap
      # updateSelectizeInput(session, "heatmap_columns", choices = names(csv_data))
      
      # Update CSV summary
      output$csv_summary <- renderPrint({
        cat("CSV loaded successfully\n")
        cat("Number of rows:", nrow(csv_data), "\n")
        cat("Number of columns:", ncol(csv_data), "\n")
        cat("First columns:", paste(head(names(csv_data), 8), collapse=", "), "\n")
      })
      
      # Update YAML structure
      values$yaml_data$`Individual general definitions`$`mapping csv file` <- csv_file$datapath
      
      # v53: cat(file=stderr(), "\n=== CSV Loaded - Checking Columns ===\n")
      # v53: cat(file=stderr(), "Available columns:", paste(names(csv_data), collapse=", "), "\n")
      # v53: cat(file=stderr(), "Number of rows:", nrow(csv_data), "\n\n")
      
      # Hide progress
      values$progress_visible <- FALSE

      # v57: Update status indicator - CSV is loaded
      if (!is.null(values$tree)) {
        show_status_click_to_generate()
      }

      # Don't auto-match - user must click "Process Data" button
      # v53: cat(file=stderr(), "CSV loaded. Ready for user to select columns and click 'Process Data'.\n")
    }, error = function(e) {
      values$progress_visible <- FALSE
      showNotification(paste("Error loading CSV:", e$message), type = "error")
    })
  })

  # Update individual value dropdown when column is selected
  observeEvent(input$individual_column, {
    # v53: cat(file=stderr(), paste("individual_column changed to:", input$individual_column), "\n")
    
    # Only proceed if we have data and a valid column selection
    if (is.null(values$csv_data)) {
      # v53: cat(file=stderr(), "No CSV data available\n")
      return()
    }
    
    if (is.null(input$individual_column) || input$individual_column == "") {
      # v53: cat(file=stderr(), "No individual column selected\n")
      updateSelectizeInput(session, "individual_value", choices = NULL, selected = NULL, server = TRUE)
      return()
    }
    
    # Check if column exists in CSV
    if (!(input$individual_column %in% names(values$csv_data))) {
      # v53: cat(file=stderr(), paste("Column", input$individual_column, "not found in CSV\n"))
      showNotification(paste("Column not found:", input$individual_column), type = "warning")
      updateSelectizeInput(session, "individual_value", choices = NULL, selected = NULL, server = TRUE)
      return()
    }
    
    # Get unique values safely
    tryCatch({
      unique_values <- unique(values$csv_data[[input$individual_column]])
      unique_values <- unique_values[!is.na(unique_values)]
      
      # v53: cat(file=stderr(), paste("Found", length(unique_values), "unique values\n"))
      # v53: cat(file=stderr(), "Unique values:", paste(head(unique_values, 10), collapse=", "), "\n")
      
      # S2.0-PERF: Use server-side selectize for large option lists (eliminates browser lag)
      updateSelectizeInput(session, "individual_value",
                        choices = unique_values,
                        selected = if(length(unique_values) > 0) unique_values[1] else NULL,
                        server = TRUE)
      
      # Store the individual column name in the YAML structure
      if (!is.null(values$yaml_data)) {
        values$yaml_data$`Individual general definitions`$`individual column` <- input$individual_column
      }
    }, error = function(e) {
      # v53: cat(file=stderr(), paste("Error getting unique values:", e$message), "\n")
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Update individual name when value is selected
  observeEvent(input$individual_value, {
    req(input$individual_value)
    
    # Update individual name in YAML
    values$yaml_data$`Individual general definitions`$Individual <- input$individual_value
    
    # Don't auto-match - user must click "Process Data" button
    # v53: cat(file=stderr(), "Individual value selected:", input$individual_value, "\n")
  })
  
  
  # When tree file is uploaded
  # When tree file is uploaded
  observeEvent(input$tree_file, {
    cat(file=stderr(), "[DEBUG-TREE] Tree file observer triggered\n")
    req(input$tree_file)

    tree_file <- input$tree_file
    cat(file=stderr(), sprintf("[DEBUG-TREE] File path: %s\n", tree_file$datapath))

    # Show progress
    values$progress_message <- "📂 Loading tree file..."
    values$progress_visible <- TRUE

    # Read tree file
    tryCatch({
      cat(file=stderr(), "[DEBUG-TREE] About to call read.tree()\n")
      tree <- read.tree(tree_file$datapath)
      cat(file=stderr(), "[DEBUG-TREE] read.tree() completed\n")
      values$tree <- tree
      # v56b: Suppress harmless fortify warnings
      values$tree_data <- suppressWarnings(ggtree(tree))$data

      # Update tree summary
      output$tree_summary <- renderPrint({
        cat("Tree loaded successfully\n")
        cat("Number of tips:", length(tree$tip.label), "\n")
        cat("Number of internal nodes:", tree$Nnode, "\n")
        cat("First few tip labels:", paste(head(tree$tip.label, 5), collapse=", "), "\n")
      })
      
      # Update nodes for manual rotation
      if (!is.null(values$tree_data)) {
        updateSelectizeInput(session, "nodes_to_rotate", 
                             choices = values$tree_data$node[values$tree_data$isTip == FALSE])
      }
      
      # Update YAML structure
      values$yaml_data$`Individual general definitions`$`tree path` <- list(tree_file$datapath)

      # S2.7-PERF: Clear cache when new tree is loaded (tree structure changed)
      values$p_list_cache <- list()
      values$p_list_cache_order <- character()
      values$cached_p_list_of_pairs <- NULL
      values$cached_p_list_hash <- NULL
      cat(file=stderr(), "[PERF-CACHE] Cache cleared (new tree loaded)\n")

      # Hide progress
      values$progress_visible <- FALSE

      # v57: Update status indicator - tree is loaded, ready to generate
      if (!is.null(values$csv_data)) {
        show_status_click_to_generate()
      }

      # Don't auto-match - user must click "Process Data" button
      # v53: cat(file=stderr(), "Tree loaded. Ready for user to select columns and click 'Process Data'.\n")
    }, error = function(e) {
      values$progress_visible <- FALSE
      showNotification(paste("Error loading tree:", e$message), type = "error")
    })
  })
  
  # When ID column is selected
  # When ID column is selected - just update YAML, don't match yet
  observeEvent(input$id_column, {
    req(input$id_column)
    
    # Update YAML structure
    if (!is.null(input$id_column) && input$id_column != "") {
      values$yaml_data$`Mapping exl renaming titles`$`ID column` <- input$id_column
    }
  })
  
  # When Process Data button is clicked
  observeEvent(input$process_data, {
    # Validate all required inputs
    if (is.null(values$tree)) {
      showNotification("Please upload a tree file first", type = "error")
      return()
    }
    
    if (is.null(values$csv_data)) {
      showNotification("Please upload a CSV file first", type = "error")
      return()
    }
    
    if (is.null(input$id_column) || input$id_column == "") {
      showNotification("Please select an ID column", type = "error")
      return()
    }
    
    # Check individual selection if not using all data
    if (!input$use_all_data) {
      if (is.null(input$individual_column) || input$individual_column == "") {
        showNotification("Please select an Individual column or check 'Use all data'", type = "warning")
        return()
      }
      
      if (is.null(input$individual_value) || input$individual_value == "") {
        showNotification("Please select an Individual value or check 'Use all data'", type = "warning")
        return()
      }
    }
    
    # Show progress
    values$progress_message <- "Matching tree IDs with CSV data..."
    values$progress_visible <- TRUE
    
    # Now call the matching function
    match_tree_with_csv()
    
    # Hide progress (will be hidden in match_tree_with_csv success/error)
    values$progress_visible <- FALSE
  })
  
  # Update YAML structure based on all input changes
  # Update YAML structure based on all input changes
  # Update YAML structure based on all input changes
  # Update YAML structure based on all input changes
  update_yaml <- function() {
    # Update Individual settings
    if (!is.null(input$individual_value) && input$individual_value != "") {
      values$yaml_data$`Individual general definitions`$Individual <- input$individual_value
    } else {
      values$yaml_data$`Individual general definitions`$Individual <- input$individual_name
    }    
    # Update output file settings
    values$yaml_data$`Individual general definitions`$out_file$base_path <- input$output_path
    values$yaml_data$`Individual general definitions`$out_file$file_type <- input$output_format
    values$yaml_data$`Individual general definitions`$out_file$`optional text at beggining` <- input$prefix_text
    values$yaml_data$`Individual general definitions`$out_file$`optional text at end` <- input$suffix_text
    
    # Use exists() or null check for all inputs that might not exist yet
    values$yaml_data$`Individual general definitions`$out_file$`replace name`$flag <- 
      if (!is.null(input$replace_name) && input$replace_name) "yes" else "no"
    values$yaml_data$`Individual general definitions`$out_file$`replace name`$name <- 
      if (!is.null(input$custom_name)) input$custom_name else "tree_plot"
    
    # Update ID column
    if (!is.null(input$id_column)) {
      values$yaml_data$`Mapping exl renaming titles`$`ID column` <- input$id_column
    }
    
    # Update basic visualization settings
    if (!is.null(input$trim_tips)) {
      values$yaml_data$`visual definitions`$`trim tips`$display <- if (input$trim_tips) "yes" else "no"
    }
    if (!is.null(input$tip_length)) {
      values$yaml_data$`visual definitions`$`trim tips`$length <- input$tip_length
    }
    if (!is.null(input$edge_width)) {
      values$yaml_data$`visual definitions`$edge_width_multiplier$size <- input$edge_width
    }
    if (!is.null(input$tip_font_size)) {
      values$yaml_data$`visual definitions`$font_size$tips <- input$tip_font_size
    }
    
    # Update Bootstrap settings
    if (!is.null(input$show_bootstrap)) {
      values$yaml_data$`visual definitions`$Bootstrap$display <- 
        if (input$show_bootstrap) "yes" else "no"
    }
    if (!is.null(input$bootstrap_format)) {
      values$yaml_data$`visual definitions`$Bootstrap$format <- input$bootstrap_format
    }
    if (!is.null(input$bootstrap_param)) {
      values$yaml_data$`visual definitions`$Bootstrap$param <- as.character(input$bootstrap_param)
    }
    
    # Update Rotation settings - ALWAYS keep both rotation1 and rotation2 in YAML
    # rotation1 (Primary)
    if (!is.null(values$rotation1_config) && length(values$rotation1_config) > 0) {
      # Build according from saved config
      acc <- list()
      for (i in seq_along(values$rotation1_config)) {
        if (!is.null(values$rotation1_config[[i]])) {
          col <- values$rotation1_config[[i]]$col
          val <- values$rotation1_config[[i]]$val
          if (!is.null(col) && col != "" && !is.null(val) && val != "") {
            acc[[as.character(i)]] <- list(title = col, value = val)
          }
        }
      }
      values$yaml_data$`visual definitions`$rotation1$according <- acc
      
      values$yaml_data$`visual definitions`$bootstrap_label_size <- input$bootstrap_label_size
      
      # Set display based on current rotation type selection
      if (!is.null(input$enable_rotation) && input$enable_rotation && 
          !is.null(input$rotation_type) && (input$rotation_type == "primary" || input$rotation_type == "manual")) {
        values$yaml_data$`visual definitions`$rotation1$display <- "yes"
      } else {
        values$yaml_data$`visual definitions`$rotation1$display <- "no"
      }
    } else {
      # No config saved, set display to no but keep structure
      values$yaml_data$`visual definitions`$rotation1$display <- "no"
      values$yaml_data$`visual definitions`$rotation1$according <- list()
    }
    
    # rotation2 (Secondary)
    if (!is.null(values$rotation2_config) && length(values$rotation2_config) > 0) {
      # Build according from saved config
      acc <- list()
      for (i in seq_along(values$rotation2_config)) {
        if (!is.null(values$rotation2_config[[i]])) {
          col <- values$rotation2_config[[i]]$col
          val <- values$rotation2_config[[i]]$val
          if (!is.null(col) && col != "" && !is.null(val) && val != "") {
            acc[[as.character(i)]] <- list(title = col, value = val)
          }
        }
      }
      values$yaml_data$`visual definitions`$rotation2$according <- acc
      
      # Set display based on current rotation type selection
      if (!is.null(input$enable_rotation) && input$enable_rotation && 
          !is.null(input$rotation_type) && input$rotation_type == "secondary") {
        values$yaml_data$`visual definitions`$rotation2$display <- "yes"
      } else {
        values$yaml_data$`visual definitions`$rotation2$display <- "no"
      }
    } else {
      # No config saved, set display to no but keep structure
      values$yaml_data$`visual definitions`$rotation2$display <- "no"
      values$yaml_data$`visual definitions`$rotation2$according <- list()
    }

    # Manual rotation (node list)
    # Ensure manual_rotation structure exists
    if (is.null(values$yaml_data$`visual definitions`$manual_rotation)) {
      values$yaml_data$`visual definitions`$manual_rotation <- list(display = "no", nodes = list())
    }

    if (!is.null(values$manual_rotation_config) && length(values$manual_rotation_config) > 0 &&
        !all(is.na(values$manual_rotation_config))) {
      # Convert node list to YAML-friendly format
      values$yaml_data$`visual definitions`$manual_rotation$nodes <- as.list(values$manual_rotation_config)

      # Set display based on current rotation type selection
      if (!is.null(input$enable_rotation) && input$enable_rotation &&
          !is.null(input$rotation_type) && input$rotation_type == "manual") {
        values$yaml_data$`visual definitions`$manual_rotation$display <- "yes"
      } else {
        values$yaml_data$`visual definitions`$manual_rotation$display <- "no"
      }
      cat(file=stderr(), sprintf("[YAML-UPDATE] Saved manual_rotation: %d nodes (%s)\n",
                                 length(values$manual_rotation_config),
                                 paste(values$manual_rotation_config, collapse = ", ")))
    } else {
      # No manual rotation config, set display to no but keep structure
      values$yaml_data$`visual definitions`$manual_rotation$display <- "no"
      values$yaml_data$`visual definitions`$manual_rotation$nodes <- list()
    }

    # Update advanced settings only if they exist
    if (!is.null(input$simulate_p_value)) {
      values$yaml_data$`visual definitions`$simulate_p_value <- 
        if (input$simulate_p_value) "yes" else "no"
    }
    
    if (!is.null(input$ladderize)) {
      values$yaml_data$`visual definitions`$laderize_flag <- 
        if (input$ladderize) "yes" else "no"
    }
    
    if (!is.null(input$display_node_numbers)) {
      values$yaml_data$`visual definitions`$flag_display_nod_number_on_tree <- 
        if (input$display_node_numbers) "yes" else "no"
    }
    
    if (!is.null(input$tip_name_display_flag)) {
      values$yaml_data$`visual definitions`$tip_name_display_flag <- 
        if (input$tip_name_display_flag) "yes" else "no"
    }
    
    # Update ID trimming settings - BUT ONLY if we dont have inferred parameters!
    # If we have inferred parameters, those take precedence over UI inputs
    # v53: cat(file=stderr(), "\nâš™ï¸ === IN update_yaml() ID TRIMMING SECTION ===\n")
    # v53: cat(file=stderr(), "âš™ï¸ values$trimming_params is NULL:", is.null(values$trimming_params), "\n")
    # v53: cat(file=stderr(), "âš™ï¸ input$id_tip_trim_flag is NULL:", is.null(input$id_tip_trim_flag), "\n")
    if (!is.null(input$id_tip_trim_flag)) {
      # v53: cat(file=stderr(), "âš™ï¸ input$id_tip_trim_flag value:", input$id_tip_trim_flag, "\n")
    }
    
    if (is.null(values$trimming_params) && !is.null(input$id_tip_trim_flag)) {
      if (is.null(values$trimming_params) && !is.null(input$id_tip_trim_flag)) {
        # v53: cat(file=stderr(), "âš™ï¸ âœ“ UPDATING trimming params from UI inputs\n")
      } else {
        # v53: cat(file=stderr(), "âš™ï¸ ✗ SKIPPING trimming params update (inferred params exist)\n")
      }
      # v53: cat(file=stderr(), "âš™ï¸ ==========================================\n\n")
      
      values$yaml_data$`visual definitions`$id_tip_trim_flag <- 
        if (input$id_tip_trim_flag) "yes" else "no"
      
      if (!is.null(input$id_tip_trim_start)) {
        values$yaml_data$`visual definitions`$id_tip_trim_start <- input$id_tip_trim_start
      }
      
      if (!is.null(input$id_tip_trim_end)) {
        values$yaml_data$`visual definitions`$id_tip_trim_end <- input$id_tip_trim_end
      }
      
      if (!is.null(input$id_tip_prefix)) {
        values$yaml_data$`visual definitions`$id_tip_prefix <- input$id_tip_prefix
      }
    }
    
    # Update advanced layout settings only if they exist
    if (!is.null(input$man_adjust_elipse)) {
      values$yaml_data$`visual definitions`$man_adjust_elipse <- input$man_adjust_elipse
    }
    if (!is.null(input$man_multiply_elipse)) {
      values$yaml_data$`visual definitions`$man_multiply_elipse <- input$man_multiply_elipse
    }
    if (!is.null(input$man_adjust_elipse_a)) {
      values$yaml_data$`visual definitions`$man_adjust_elipse_a <- input$man_adjust_elipse_a
    }
    if (!is.null(input$man_adjust_elipse_b)) {
      values$yaml_data$`visual definitions`$man_adjust_elipse_b <- input$man_adjust_elipse_b
    }
    if (!is.null(input$man_adj_second_legend)) {
      values$yaml_data$`visual definitions`$man_adj_second_legend <- input$man_adj_second_legend
    }
    if (!is.null(input$man_space_second_legend)) {
      values$yaml_data$`visual definitions`$man_space_second_legend <- input$man_space_second_legend
    }
    if (!is.null(input$man_adjust_image_of_second_legend)) {
      values$yaml_data$`visual definitions`$man_adjust_image_of_second_legend <- input$man_adjust_image_of_second_legend
    }
    if (!is.null(input$man_multiply_second_legend)) {
      values$yaml_data$`visual definitions`$man_multiply_second_legend <- input$man_multiply_second_legend
    }
    if (!is.null(input$man_multiply_second_legend_text)) {
      values$yaml_data$`visual definitions`$man_multiply_second_legend_text <- input$man_multiply_second_legend_text
    }
    if (!is.null(input$man_multiply_first_legend_text)) {
      values$yaml_data$`visual definitions`$man_multiply_first_legend_text <- input$man_multiply_first_legend_text
    }
    if (!is.null(input$man_multiply_first_legend_title_size)) {
      values$yaml_data$`visual definitions`$man_multiply_first_legend_title_size <- input$man_multiply_first_legend_title_size
    }
    if (!is.null(input$man_space_second_legend_multiplier)) {
      values$yaml_data$`visual definitions`$man_space_second_legend_multiplier <- input$man_space_second_legend_multiplier
    }
    if (!is.null(input$man_offset_for_highlight_legend_x)) {
      values$yaml_data$`visual definitions`$man_offset_for_highlight_legend_x <- input$man_offset_for_highlight_legend_x
    }
    if (!is.null(input$man_offset_second_legend)) {
      values$yaml_data$`visual definitions`$man_offset_second_legend <- input$man_offset_second_legend
    }
    if (!is.null(input$man_boot_x_offset)) {
      values$yaml_data$`visual definitions`$man_boot_x_offset <- input$man_boot_x_offset
    }
    if (!is.null(input$man_adj_heat_loc)) {
      values$yaml_data$`visual definitions`$man_adj_heat_loc <- input$man_adj_heat_loc
    }
    if (!is.null(input$man_adj_heat_loc2)) {
      values$yaml_data$`visual definitions`$man_adj_heat_loc2 <- input$man_adj_heat_loc2
    }
    if (!is.null(input$man_adj_heat_loc3)) {
      values$yaml_data$`visual definitions`$man_adj_heat_loc3 <- input$man_adj_heat_loc3
    }
    if (!is.null(input$add_date_to_text_flag)) {
      values$yaml_data$`visual definitions`$add_date_to_text_flag <- 
        if (input$add_date_to_text_flag) "yes" else "no"
    }
    if (!is.null(input$a_4_output)) {
      values$yaml_data$`visual definitions`$a_4_output <- 
        if (input$a_4_output) "yes" else "no"
    }
    
    # Update debug settings
    if (!is.null(input$debug_mode)) {
      values$yaml_data$`visual definitions`$debug_mode <- 
        if (input$debug_mode) "yes" else "no"
    }
    if (!is.null(input$debug_print_data_tree)) {
      values$yaml_data$`visual definitions`$debug_print_data_tree <- 
        if (input$debug_print_data_tree) "yes" else "no"
    }
    if (!is.null(input$flag_calc_scores_for_tree)) {
      values$yaml_data$`visual definitions`$flag_calc_scores_for_tree <- 
        if (input$flag_calc_scores_for_tree) "yes" else "no"
    }
    
    # Update Newick export settings
    if (!is.null(input$flag_make_newick_file)) {
      values$yaml_data$`visual definitions`$flag_make_newick_file <- 
        if (input$flag_make_newick_file) "yes" else "no"
    }
    if (!is.null(input$path_out_newick)) {
      values$yaml_data$`visual definitions`$path_out_newick <- input$path_out_newick
    }
    
    # Update heat legend replace
    if (!is.null(input$heat_legend_replace) && nchar(input$heat_legend_replace) > 0) {
      pairs <- strsplit(input$heat_legend_replace, ",")[[1]]
      heat_legend_replace <- list()
      
      for (pair in pairs) {
        if (grepl(":", pair)) {
          kv <- strsplit(pair, ":")[[1]]
          if (length(kv) == 2) {
            key <- trimws(kv[1])
            value <- trimws(kv[2])
            heat_legend_replace[[key]] <- value
          }
        }
      }
      
      if (length(heat_legend_replace) > 0) {
        values$yaml_data$`visual definitions`$heat_legend_replace <- heat_legend_replace
      }
    } else {
      values$yaml_data$`visual definitions`$heat_legend_replace <- NULL
    }
    
    # Update classification section
    values$yaml_data$`visual definitions`$classification <- list()
    
    # === CRITICAL FIX: Determine which classification(s) to use ===
    # This fixes the preview and selection functionality
    classifications_to_use <- NULL
    
    if (!is.null(values$temp_classification_preview)) {
      # PREVIEW MODE: User clicked "Update Preview" - show temp classification
      # v53: cat(file=stderr(), "🎨 YAML GENERATION: Using temp_classification_preview (PREVIEW MODE)\n")
      classifications_to_use <- list(values$temp_classification_preview)
      
    } else if (!is.null(values$active_classification_index) && 
               !is.null(values$classifications) &&
               length(values$classifications) >= values$active_classification_index) {
      # SELECTION MODE: User selected a specific classification via radio button
      # v53: cat(file=stderr(), "Ã°Å¸â€Ëœ YAML GENERATION: Using classification #", values$active_classification_index, " (SELECTION MODE)\n")
      classifications_to_use <- list(values$classifications[[values$active_classification_index]])
      
    } else if (!is.null(values$classifications) && length(values$classifications) > 0) {
      # DEFAULT MODE: Use the last classification
      # v53: cat(file=stderr(), "📋 YAML GENERATION: Using last classification (DEFAULT MODE)\n")
      classifications_to_use <- list(values$classifications[[length(values$classifications)]])
    }
    
    # If user has defined classifications, use the selected one(s)
    if (!is.null(classifications_to_use) && length(classifications_to_use) > 0) {
      for (i in seq_along(classifications_to_use)) {
        class_item <- list()
        class_def <- classifications_to_use[[i]]
        
        class_item[[as.character(i)]] <- list(
          according = list(),
          FDR_perc = class_def$fdr,
          non_cluster_title = class_def$no_cluster_title,
          non_cluster_color = class_def$no_cluster_color,
          not_in_legend = list("No cluster", "0", 0, 0, "0"),
          title = class_def$title,
          na_name = "na"
        )
        
        # Add classes
        if (!is.null(class_def$classes)) {
          for (j in seq_along(class_def$classes)) {
            class_according <- list()
            class_entry <- class_def$classes[[j]]
            
            class_according[[as.character(j)]] <- list(
              title1 = class_entry$column,
              value1 = list(class_entry$value),
              display_name = class_entry$display_name,
              color = class_entry$color
            )
            
            class_item[[as.character(i)]]$according <- c(
              class_item[[as.character(i)]]$according,
              list(class_according)
            )
          }
        } else {
          # Add a default "according" entry if no classes are defined
          default_according <- list()
          default_according[["1"]] <- list(
            title1 = input$id_column,
            value1 = list("Default"),
            display_name = "Default",
            color = "gray"
          )
          
          class_item[[as.character(i)]]$according <- list(default_according)
        }
        
        # v56: Add heatmaps - check both class_def$heatmaps and values$heatmaps
        heatmaps_to_use <- NULL
        if (!is.null(class_def$heatmaps) && length(class_def$heatmaps) > 0) {
          heatmaps_to_use <- class_def$heatmaps
        } else if (!is.null(values$heatmaps) && length(values$heatmaps) > 0) {
          # v56: Use values$heatmaps from the new multi-heatmap interface
          heatmaps_to_use <- values$heatmaps
        }
        
        if (!is.null(heatmaps_to_use) && length(heatmaps_to_use) > 0) {
          # v56c: DEBUG
          debug_cat(paste0("\n=== v56c: Adding ", length(heatmaps_to_use), " heatmap(s) to classification ", i, " ===\n"))

          class_item[[as.character(i)]]$heatmap_display <- list()

          for (j in seq_along(heatmaps_to_use)) {
            heatmap_entry <- heatmaps_to_use[[j]]

            # v56c: DEBUG
            debug_cat(paste0("  Processing heatmap ", j, ": ", heatmap_entry$title, "\n"))
            debug_cat(paste0("    Columns: ", paste(heatmap_entry$columns, collapse=", "), "\n"))

            heatmap_item <- list()
            heatmap_item[[as.character(j)]] <- list(
              display = "yes",
              title = heatmap_entry$title,
              is_discrete = if (heatmap_entry$is_discrete) "yes" else "no",
              according = list()
            )
            
            if (heatmap_entry$is_discrete) {
              # v70: Fixed bug - use man_define_colors (not use_custom_colors) and custom_colors (not colors)
              if (!is.null(heatmap_entry$man_define_colors) && heatmap_entry$man_define_colors) {
                heatmap_item[[as.character(j)]]$man_define_colors <- "yes"
                # v104: Store custom colors as a list with explicit names to preserve value->color mapping
                custom_colors_as_list <- as.list(heatmap_entry$custom_colors)
                names(custom_colors_as_list) <- names(heatmap_entry$custom_colors)
                heatmap_item[[as.character(j)]]$color_scale_option <- custom_colors_as_list
              } else {
                heatmap_item[[as.character(j)]]$color_scale_option <- heatmap_entry$color_scheme
              }
              # v70: Add NA color (default white)
              heatmap_item[[as.character(j)]]$na_color <- if (!is.null(heatmap_entry$na_color)) heatmap_entry$na_color else "white"
            } else {
              heatmap_item[[as.character(j)]]$low <- heatmap_entry$low_color
              heatmap_item[[as.character(j)]]$mid <- if (!is.null(heatmap_entry$mid_color)) heatmap_entry$mid_color else heatmap_entry$low_color
              heatmap_item[[as.character(j)]]$high <- heatmap_entry$high_color
              heatmap_item[[as.character(j)]]$midpoint <- if (!is.null(heatmap_entry$midpoint)) heatmap_entry$midpoint else 0
              # v112: Add NA color for continuous heatmaps (default grey90)
              heatmap_item[[as.character(j)]]$na_color <- if (!is.null(heatmap_entry$na_color)) heatmap_entry$na_color else "grey90"
            }

            # v104: Add per-heatmap distance
            heatmap_item[[as.character(j)]]$distance <- if (!is.null(heatmap_entry$distance)) heatmap_entry$distance else 0.02

            # v105: Add per-heatmap height (now called Column Width in UI)
            heatmap_item[[as.character(j)]]$height <- if (!is.null(heatmap_entry$height)) heatmap_entry$height else 0.8

            # v110: Add per-heatmap row height
            heatmap_item[[as.character(j)]]$row_height <- if (!is.null(heatmap_entry$row_height)) heatmap_entry$row_height else 1.0

            # v111: Add grid settings
            heatmap_item[[as.character(j)]]$show_grid <- if (!is.null(heatmap_entry$show_grid) && heatmap_entry$show_grid) "yes" else "no"
            heatmap_item[[as.character(j)]]$grid_color <- if (!is.null(heatmap_entry$grid_color)) heatmap_entry$grid_color else "#000000"
            heatmap_item[[as.character(j)]]$grid_size <- if (!is.null(heatmap_entry$grid_size)) heatmap_entry$grid_size else 0.5

            # S1.62dev: Add row line settings (horizontal lines only)
            heatmap_item[[as.character(j)]]$show_row_lines <- if (!is.null(heatmap_entry$show_row_lines) && heatmap_entry$show_row_lines) "yes" else "no"
            heatmap_item[[as.character(j)]]$row_line_color <- if (!is.null(heatmap_entry$row_line_color)) heatmap_entry$row_line_color else "#000000"
            heatmap_item[[as.character(j)]]$row_line_size <- if (!is.null(heatmap_entry$row_line_size)) heatmap_entry$row_line_size else 0.5

            # S1.62dev: Add column line settings (vertical lines)
            heatmap_item[[as.character(j)]]$show_col_lines <- if (!is.null(heatmap_entry$show_col_lines) && heatmap_entry$show_col_lines) "yes" else "no"
            heatmap_item[[as.character(j)]]$col_line_color <- if (!is.null(heatmap_entry$col_line_color)) heatmap_entry$col_line_color else "#000000"
            heatmap_item[[as.character(j)]]$col_line_size <- if (!is.null(heatmap_entry$col_line_size)) heatmap_entry$col_line_size else 0.5

            # S1.62dev: Add vertical text labels settings
            heatmap_item[[as.character(j)]]$show_vertical_text <- if (!is.null(heatmap_entry$show_vertical_text) && heatmap_entry$show_vertical_text) "yes" else "no"
            heatmap_item[[as.character(j)]]$vertical_text_column <- if (!is.null(heatmap_entry$vertical_text_column)) heatmap_entry$vertical_text_column else ""
            heatmap_item[[as.character(j)]]$vertical_text_size <- if (!is.null(heatmap_entry$vertical_text_size)) heatmap_entry$vertical_text_size else 3
            heatmap_item[[as.character(j)]]$vertical_text_offset <- if (!is.null(heatmap_entry$vertical_text_offset)) heatmap_entry$vertical_text_offset else 0.5
            heatmap_item[[as.character(j)]]$vertical_text_color <- if (!is.null(heatmap_entry$vertical_text_color)) heatmap_entry$vertical_text_color else "#000000"

            # v118: Add guide line settings (was missing - this caused tip guide lines not to work!)
            heatmap_item[[as.character(j)]]$show_guides <- if (!is.null(heatmap_entry$show_guides) && heatmap_entry$show_guides) "yes" else "no"
            heatmap_item[[as.character(j)]]$guide_color1 <- if (!is.null(heatmap_entry$guide_color1)) heatmap_entry$guide_color1 else "#CCCCCC"
            heatmap_item[[as.character(j)]]$guide_color2 <- if (!is.null(heatmap_entry$guide_color2)) heatmap_entry$guide_color2 else "#EEEEEE"
            heatmap_item[[as.character(j)]]$guide_alpha <- if (!is.null(heatmap_entry$guide_alpha)) heatmap_entry$guide_alpha else 0.3
            heatmap_item[[as.character(j)]]$guide_width <- if (!is.null(heatmap_entry$guide_width)) heatmap_entry$guide_width else 0.5
            heatmap_item[[as.character(j)]]$guide_linetype <- if (!is.null(heatmap_entry$guide_linetype)) heatmap_entry$guide_linetype else "solid"

            # v109: Add colnames_angle
            heatmap_item[[as.character(j)]]$colnames_angle <- if (!is.null(heatmap_entry$colnames_angle)) heatmap_entry$colnames_angle else 45

            # S1.62dev: Add show_colnames flag for column name visibility
            heatmap_item[[as.character(j)]]$show_colnames <- if (!is.null(heatmap_entry$show_colnames) && heatmap_entry$show_colnames) "yes" else "no"

            # v105: Add row labels settings
            heatmap_item[[as.character(j)]]$show_row_labels <- if (!is.null(heatmap_entry$show_row_labels) && heatmap_entry$show_row_labels) "yes" else "no"
            heatmap_item[[as.character(j)]]$row_label_source <- if (!is.null(heatmap_entry$row_label_source)) heatmap_entry$row_label_source else "colnames"
            heatmap_item[[as.character(j)]]$row_label_font_size <- if (!is.null(heatmap_entry$row_label_font_size)) heatmap_entry$row_label_font_size else 2.5
            # v111: Add row label offset and alignment
            heatmap_item[[as.character(j)]]$row_label_offset <- if (!is.null(heatmap_entry$row_label_offset)) heatmap_entry$row_label_offset else 1.0
            heatmap_item[[as.character(j)]]$row_label_align <- if (!is.null(heatmap_entry$row_label_align)) heatmap_entry$row_label_align else "left"
            heatmap_item[[as.character(j)]]$custom_row_labels <- if (!is.null(heatmap_entry$custom_row_labels)) heatmap_entry$custom_row_labels else ""
            # v108: Add label mapping
            if (!is.null(heatmap_entry$label_mapping) && length(heatmap_entry$label_mapping) > 0) {
              heatmap_item[[as.character(j)]]$label_mapping <- heatmap_entry$label_mapping
            }

            # S2.0-RDATA: Add data source for RData heatmaps (was missing - caused RData heatmaps to fail in custom classification path!)
            # Note: cnv_matrix is NOT serialized to YAML - it's passed as a parameter to func.print.lineage.tree
            if (!is.null(heatmap_entry$data_source) && heatmap_entry$data_source == "rdata") {
              heatmap_item[[as.character(j)]]$data_source <- "rdata"
              heatmap_item[[as.character(j)]]$use_midpoint <- "yes"  # Always use midpoint for CNV
              # Store CNV settings (but NOT the matrix itself - that's passed separately)
              heatmap_item[[as.character(j)]]$cnv_downsample <- if (!is.null(heatmap_entry$cnv_downsample)) heatmap_entry$cnv_downsample else 10
              heatmap_item[[as.character(j)]]$cnv_render_downsample <- if (!is.null(heatmap_entry$cnv_render_downsample)) heatmap_entry$cnv_render_downsample else 10
              heatmap_item[[as.character(j)]]$cnv_wgd_norm <- if (!is.null(heatmap_entry$cnv_wgd_norm) && heatmap_entry$cnv_wgd_norm) "yes" else "no"
              # S2.12: Per-cell WGD normalization settings
              heatmap_item[[as.character(j)]]$cnv_wgd_per_cell <- if (!is.null(heatmap_entry$cnv_wgd_per_cell) && heatmap_entry$cnv_wgd_per_cell) "yes" else "no"
              heatmap_item[[as.character(j)]]$cnv_wgd_column <- if (!is.null(heatmap_entry$cnv_wgd_column)) heatmap_entry$cnv_wgd_column else ""
              # S2.8: Display mode (basic or detailed)
              heatmap_item[[as.character(j)]]$cnv_display_mode <- if (!is.null(heatmap_entry$cnv_display_mode)) heatmap_entry$cnv_display_mode else "basic"
              # Height scale for detailed mode
              heatmap_item[[as.character(j)]]$cnv_height_scale <- if (!is.null(heatmap_entry$cnv_height_scale)) heatmap_entry$cnv_height_scale else 1.0
              # S2.0: Store mapping column for sample name matching
              heatmap_item[[as.character(j)]]$rdata_mapping_column <- heatmap_entry$rdata_mapping_column
              debug_cat(paste0("    S2.0-RDATA: RData heatmap, mapping_column=", heatmap_entry$rdata_mapping_column, "\n"))
              debug_cat(paste0("    S2.12: Per-cell WGD: enabled=", heatmap_entry$cnv_wgd_per_cell, ", column=", heatmap_entry$cnv_wgd_column, "\n"))
              debug_cat(paste0("    Height scale: ", heatmap_entry$cnv_height_scale, "\n"))
              debug_cat(paste0("    S2.8: Display mode: ", heatmap_entry$cnv_display_mode, "\n"))
            } else {
              heatmap_item[[as.character(j)]]$data_source <- "csv"
              # Add columns - format must match expected YAML structure
              # Each column entry needs to be a named list like list("1" = "column_name")
              if (!is.null(heatmap_entry$columns)) {
                for (k in seq_along(heatmap_entry$columns)) {
                  column_entry <- list()
                  column_entry[[as.character(k)]] <- heatmap_entry$columns[k]
                  heatmap_item[[as.character(j)]]$according[[k]] <- column_entry
                }
              }
            }

            class_item[[as.character(i)]]$heatmap_display[[j]] <- heatmap_item
            # S2.3-DEBUG: Show according columns for verification (custom classification path)
            cat(file=stderr(), paste0("[YAML-DEBUG-CUSTOM] Heatmap ", j, " YAML according columns: "))
            if (!is.null(heatmap_item[[as.character(j)]]$according) && length(heatmap_item[[as.character(j)]]$according) > 0) {
              for (k_debug in seq_along(heatmap_item[[as.character(j)]]$according)) {
                acc_entry <- heatmap_item[[as.character(j)]]$according[[k_debug]]
                acc_key <- names(acc_entry)[1]
                cat(file=stderr(), paste0(acc_entry[[acc_key]], " "))
              }
              cat(file=stderr(), "\n")
            } else {
              cat(file=stderr(), "(none)\n")
            }
          }
        }

        # ========================================================
        # === ADD HIGHLIGHT TO THIS CLASSIFICATION ===
        # *** KEY FIX: MOVED HERE - BEFORE adding class_item to YAML ***
        # ========================================================
        highlight_to_apply <- NULL

        # v131: DEBUG - trace highlight decision
        debug_cat(paste0("\n=== v131: HIGHLIGHT DECISION (classification ", i, ") ===\n"))
        debug_cat(paste0("  temp_highlight_preview is NULL: ", is.null(values$temp_highlight_preview), "\n"))
        debug_cat(paste0("  active_highlight_index: ", values$active_highlight_index, "\n"))
        debug_cat(paste0("  highlights length: ", length(values$highlights), "\n"))

        # v53: cat(file=stderr(), "\n📝 === ADDING HIGHLIGHT TO CLASSIFICATION", i, "===\n")
        
        if (!is.null(values$temp_highlight_preview)) {
          # PREVIEW MODE: User clicked "Update Preview"
          # v53: cat(file=stderr(), "🎨 Using temp_highlight_preview (PREVIEW MODE)\n")
          # v53: cat(file=stderr(), "   Column:", values$temp_highlight_preview$column, "\n")
          # v53: cat(file=stderr(), "   Items:", length(values$temp_highlight_preview$items), "\n")
          highlight_to_apply <- values$temp_highlight_preview
          
        } else if (!is.null(values$active_highlight_index) && 
                   !is.null(values$highlights) &&
                   length(values$highlights) >= values$active_highlight_index) {
          # SAVED MODE: User selected a saved highlight
          # v53: cat(file=stderr(), "Ã°Å¸â€Ëœ Using saved highlight #", values$active_highlight_index, "\n")
          highlight_to_apply <- values$highlights[[values$active_highlight_index]]
        } else {
          debug_cat("  NO highlight source found (will set display='no')\n")
        }

        # Apply highlight to THIS classification
        if (!is.null(highlight_to_apply)) {
          debug_cat(paste0("  APPLYING highlight with ", length(highlight_to_apply$items), " items (display='yes')\n"))
          # v53: cat(file=stderr(), "Ã¢Å“â€œ Applying highlight to class_item\n")
          
          # Build highlight YAML structure
          highlight_yaml <- list(
            display = "yes",
            offset = highlight_to_apply$offset,
            vertical_offset = if(!is.null(highlight_to_apply$vertical_offset)) highlight_to_apply$vertical_offset else 0,
            adjust_height = highlight_to_apply$adjust_height,
            adjust_width = highlight_to_apply$adjust_width,
            according = list()
          )
          
          # Add each highlighted value
          for (j in seq_along(highlight_to_apply$items)) {
            item <- highlight_to_apply$items[[j]]
            acc_item <- list()
            acc_item[[as.character(j)]] <- list(
              title1 = highlight_to_apply$column,
              value1 = item$value,
              display_name = item$display_name,
              color = item$color,
              transparency = if (!is.null(item$transparency)) item$transparency else 0.5,  # v139: Add transparency
              display_title = highlight_to_apply$title
            )
            highlight_yaml$according <- c(highlight_yaml$according, list(acc_item))
          }
          
          # ADD TO CLASS_ITEM (the key fix!)
          class_item[[as.character(i)]]$highlight <- highlight_yaml
          
          # v53: cat(file=stderr(), "Ã¢Å“â€œ Highlight added to class_item[[", i, "]]\n")
          
        } else {
          # No highlight - disable it
          # v53: cat(file=stderr(), "Ã¢Å“â€” No highlight to apply\n")
          class_item[[as.character(i)]]$highlight <- list(display = "no")
        }
        # ========================================================
        # === END OF HIGHLIGHT CODE ===
        # ========================================================
        
        # NOW add class_item to YAML (with highlight already included)
        values$yaml_data$`visual definitions`$classification <- c(
          values$yaml_data$`visual definitions`$classification,
          list(class_item)
        )
      }
    } 
    # If no user-defined classifications, create a default classification
    else {
      # Create a default classification entry
      default_classification <- list()
      default_classification[["1"]] <- list(
        according = list(),
        FDR_perc = 0.1,
        non_cluster_title = "No cluster",
        non_cluster_color = "gray",
        not_in_legend = list("No cluster", "0", 0, 0, "0"),
        title = "Cell type",
        na_name = "na"
      )
      
      # Add a default "according" entry
      default_according <- list()
      default_according[["1"]] <- list(
        title1 = if (!is.null(input$id_column)) input$id_column else "id",
        value1 = list("Default"),
        display_name = "Default",
        color = "gray"
      )
      
      default_classification[["1"]]$according <- list(default_according)

      # v56: Add heatmaps to default classification if values$heatmaps is set
      if (!is.null(values$heatmaps) && length(values$heatmaps) > 0) {
        # v56c: DEBUG
        debug_cat(paste0("\n=== v56c: Adding ", length(values$heatmaps), " heatmap(s) to DEFAULT classification ===\n"))

        default_classification[["1"]]$heatmap_display <- list()

        for (j in seq_along(values$heatmaps)) {
          heatmap_entry <- values$heatmaps[[j]]

          # v56c: DEBUG
          debug_cat(paste0("  Processing heatmap ", j, ": ", heatmap_entry$title, "\n"))
          debug_cat(paste0("    Columns: ", paste(heatmap_entry$columns, collapse=", "), "\n"))
          # S1.62dev: Debug data_source
          debug_cat(paste0("    S1.62dev data_source: ", if (!is.null(heatmap_entry$data_source)) heatmap_entry$data_source else "NULL", "\n"))
          debug_cat(paste0("    S1.62dev all fields: ", paste(names(heatmap_entry), collapse=", "), "\n"))

          heatmap_item <- list()
          heatmap_item[[as.character(j)]] <- list(
            display = "yes",
            title = heatmap_entry$title,
            is_discrete = if (heatmap_entry$is_discrete) "yes" else "no",
            according = list()
          )

          if (heatmap_entry$is_discrete) {
            # v69: Check for custom colors (man_define_colors flag)
            if (!is.null(heatmap_entry$man_define_colors) && heatmap_entry$man_define_colors) {
              heatmap_item[[as.character(j)]]$man_define_colors <- "yes"
              # v104: Store custom colors as a list with explicit names to preserve value->color mapping
              # Named vectors can lose their names when passed through list structures
              custom_colors_as_list <- as.list(heatmap_entry$custom_colors)
              names(custom_colors_as_list) <- names(heatmap_entry$custom_colors)
              heatmap_item[[as.character(j)]]$color_scale_option <- custom_colors_as_list
              debug_cat(paste0("    v104: Storing ", length(custom_colors_as_list), " custom colors with names: ", paste(names(custom_colors_as_list), collapse=", "), "\n"))
            } else {
              heatmap_item[[as.character(j)]]$color_scale_option <- heatmap_entry$color_scheme
            }
            # v70: Add NA color (default white)
            heatmap_item[[as.character(j)]]$na_color <- if (!is.null(heatmap_entry$na_color)) heatmap_entry$na_color else "white"
          } else {
            heatmap_item[[as.character(j)]]$low <- heatmap_entry$low_color
            heatmap_item[[as.character(j)]]$mid <- if (!is.null(heatmap_entry$mid_color)) heatmap_entry$mid_color else heatmap_entry$low_color
            heatmap_item[[as.character(j)]]$high <- heatmap_entry$high_color
            heatmap_item[[as.character(j)]]$midpoint <- if (!is.null(heatmap_entry$midpoint)) heatmap_entry$midpoint else 0
            # v114: Add NA color for continuous heatmaps (was missing in default classification path)
            heatmap_item[[as.character(j)]]$na_color <- if (!is.null(heatmap_entry$na_color)) heatmap_entry$na_color else "grey90"
          }

          # v104: Add per-heatmap distance
          heatmap_item[[as.character(j)]]$distance <- if (!is.null(heatmap_entry$distance)) heatmap_entry$distance else 0.02

          # v105: Add per-heatmap height
          heatmap_item[[as.character(j)]]$height <- if (!is.null(heatmap_entry$height)) heatmap_entry$height else 0.8

          # v114: Add per-heatmap row height (was missing in default classification path)
          heatmap_item[[as.character(j)]]$row_height <- if (!is.null(heatmap_entry$row_height)) heatmap_entry$row_height else 1.0

          # v114: Add grid settings (were missing in default classification path)
          heatmap_item[[as.character(j)]]$show_grid <- if (!is.null(heatmap_entry$show_grid) && heatmap_entry$show_grid) "yes" else "no"
          heatmap_item[[as.character(j)]]$grid_color <- if (!is.null(heatmap_entry$grid_color)) heatmap_entry$grid_color else "#000000"
          heatmap_item[[as.character(j)]]$grid_size <- if (!is.null(heatmap_entry$grid_size)) heatmap_entry$grid_size else 0.5

          # S1.62dev: Add row line settings (horizontal lines only)
          heatmap_item[[as.character(j)]]$show_row_lines <- if (!is.null(heatmap_entry$show_row_lines) && heatmap_entry$show_row_lines) "yes" else "no"
          heatmap_item[[as.character(j)]]$row_line_color <- if (!is.null(heatmap_entry$row_line_color)) heatmap_entry$row_line_color else "#000000"
          heatmap_item[[as.character(j)]]$row_line_size <- if (!is.null(heatmap_entry$row_line_size)) heatmap_entry$row_line_size else 0.5

          # S1.62dev: Add column line settings (vertical lines)
          heatmap_item[[as.character(j)]]$show_col_lines <- if (!is.null(heatmap_entry$show_col_lines) && heatmap_entry$show_col_lines) "yes" else "no"
          heatmap_item[[as.character(j)]]$col_line_color <- if (!is.null(heatmap_entry$col_line_color)) heatmap_entry$col_line_color else "#000000"
          heatmap_item[[as.character(j)]]$col_line_size <- if (!is.null(heatmap_entry$col_line_size)) heatmap_entry$col_line_size else 0.5

          # S1.62dev: Add vertical text labels settings
          heatmap_item[[as.character(j)]]$show_vertical_text <- if (!is.null(heatmap_entry$show_vertical_text) && heatmap_entry$show_vertical_text) "yes" else "no"
          heatmap_item[[as.character(j)]]$vertical_text_column <- if (!is.null(heatmap_entry$vertical_text_column)) heatmap_entry$vertical_text_column else ""
          heatmap_item[[as.character(j)]]$vertical_text_size <- if (!is.null(heatmap_entry$vertical_text_size)) heatmap_entry$vertical_text_size else 3
          heatmap_item[[as.character(j)]]$vertical_text_offset <- if (!is.null(heatmap_entry$vertical_text_offset)) heatmap_entry$vertical_text_offset else 0.5
          heatmap_item[[as.character(j)]]$vertical_text_color <- if (!is.null(heatmap_entry$vertical_text_color)) heatmap_entry$vertical_text_color else "#000000"

          # v120: Add guide line settings (were missing in default classification path - caused tip guide lines not to work!)
          heatmap_item[[as.character(j)]]$show_guides <- if (!is.null(heatmap_entry$show_guides) && heatmap_entry$show_guides) "yes" else "no"
          heatmap_item[[as.character(j)]]$guide_color1 <- if (!is.null(heatmap_entry$guide_color1)) heatmap_entry$guide_color1 else "#CCCCCC"
          heatmap_item[[as.character(j)]]$guide_color2 <- if (!is.null(heatmap_entry$guide_color2)) heatmap_entry$guide_color2 else "#EEEEEE"
          heatmap_item[[as.character(j)]]$guide_alpha <- if (!is.null(heatmap_entry$guide_alpha)) heatmap_entry$guide_alpha else 0.3
          heatmap_item[[as.character(j)]]$guide_width <- if (!is.null(heatmap_entry$guide_width)) heatmap_entry$guide_width else 0.5
          heatmap_item[[as.character(j)]]$guide_linetype <- if (!is.null(heatmap_entry$guide_linetype)) heatmap_entry$guide_linetype else "solid"

          # v109: Add colnames_angle
          heatmap_item[[as.character(j)]]$colnames_angle <- if (!is.null(heatmap_entry$colnames_angle)) heatmap_entry$colnames_angle else 45

          # S1.62dev: Add show_colnames flag for column name visibility
          heatmap_item[[as.character(j)]]$show_colnames <- if (!is.null(heatmap_entry$show_colnames) && heatmap_entry$show_colnames) "yes" else "no"

          # v105: Add row labels settings
          heatmap_item[[as.character(j)]]$show_row_labels <- if (!is.null(heatmap_entry$show_row_labels) && heatmap_entry$show_row_labels) "yes" else "no"
          heatmap_item[[as.character(j)]]$row_label_source <- if (!is.null(heatmap_entry$row_label_source)) heatmap_entry$row_label_source else "colnames"
          heatmap_item[[as.character(j)]]$row_label_font_size <- if (!is.null(heatmap_entry$row_label_font_size)) heatmap_entry$row_label_font_size else 2.5
          # v114: Add row label offset and alignment (were missing in default classification path)
          heatmap_item[[as.character(j)]]$row_label_offset <- if (!is.null(heatmap_entry$row_label_offset)) heatmap_entry$row_label_offset else 1.0
          heatmap_item[[as.character(j)]]$row_label_align <- if (!is.null(heatmap_entry$row_label_align)) heatmap_entry$row_label_align else "left"
          heatmap_item[[as.character(j)]]$custom_row_labels <- if (!is.null(heatmap_entry$custom_row_labels)) heatmap_entry$custom_row_labels else ""
          # v108: Add label mapping
          if (!is.null(heatmap_entry$label_mapping) && length(heatmap_entry$label_mapping) > 0) {
            heatmap_item[[as.character(j)]]$label_mapping <- heatmap_entry$label_mapping
          }

          # S1.62dev: Add data source for RData heatmaps
          # Note: cnv_matrix is NOT serialized to YAML - it's passed as a parameter to func.print.lineage.tree
          if (!is.null(heatmap_entry$data_source) && heatmap_entry$data_source == "rdata") {
            heatmap_item[[as.character(j)]]$data_source <- "rdata"
            heatmap_item[[as.character(j)]]$use_midpoint <- "yes"  # Always use midpoint for CNV
            # Store CNV settings (but NOT the matrix itself - that's passed separately)
            heatmap_item[[as.character(j)]]$cnv_downsample <- if (!is.null(heatmap_entry$cnv_downsample)) heatmap_entry$cnv_downsample else 10
            heatmap_item[[as.character(j)]]$cnv_render_downsample <- if (!is.null(heatmap_entry$cnv_render_downsample)) heatmap_entry$cnv_render_downsample else 10
            heatmap_item[[as.character(j)]]$cnv_wgd_norm <- if (!is.null(heatmap_entry$cnv_wgd_norm) && heatmap_entry$cnv_wgd_norm) "yes" else "no"
            # S2.12: Per-cell WGD normalization settings
            heatmap_item[[as.character(j)]]$cnv_wgd_per_cell <- if (!is.null(heatmap_entry$cnv_wgd_per_cell) && heatmap_entry$cnv_wgd_per_cell) "yes" else "no"
            heatmap_item[[as.character(j)]]$cnv_wgd_column <- if (!is.null(heatmap_entry$cnv_wgd_column)) heatmap_entry$cnv_wgd_column else ""
            # S2.8: Display mode (basic or detailed)
            heatmap_item[[as.character(j)]]$cnv_display_mode <- if (!is.null(heatmap_entry$cnv_display_mode)) heatmap_entry$cnv_display_mode else "basic"
            # Height scale for detailed mode
            heatmap_item[[as.character(j)]]$cnv_height_scale <- if (!is.null(heatmap_entry$cnv_height_scale)) heatmap_entry$cnv_height_scale else 1.0
            # S2.0: Store mapping column for sample name matching
            heatmap_item[[as.character(j)]]$rdata_mapping_column <- heatmap_entry$rdata_mapping_column
            debug_cat(paste0("    S2.0: RData heatmap, mapping_column=", heatmap_entry$rdata_mapping_column, "\n"))
            debug_cat(paste0("    S2.12: Per-cell WGD: enabled=", heatmap_entry$cnv_wgd_per_cell, ", column=", heatmap_entry$cnv_wgd_column, "\n"))
            debug_cat(paste0("    S2.8: Display mode: ", heatmap_entry$cnv_display_mode, "\n"))
            debug_cat(paste0("    Height scale: ", heatmap_entry$cnv_height_scale, "\n"))
          } else {
            heatmap_item[[as.character(j)]]$data_source <- "csv"
            # Add columns - format must match expected YAML structure
            if (!is.null(heatmap_entry$columns)) {
              for (k in seq_along(heatmap_entry$columns)) {
                column_entry <- list()
                column_entry[[as.character(k)]] <- heatmap_entry$columns[k]
                heatmap_item[[as.character(j)]]$according[[k]] <- column_entry
              }
            }
          }

          default_classification[["1"]]$heatmap_display[[j]] <- heatmap_item
          # S1.62dev: Debug final heatmap_item structure
          debug_cat(paste0("    S1.62dev final heatmap_item[[", as.character(j), "]] fields: ",
                          paste(names(heatmap_item[[as.character(j)]]), collapse=", "), "\n"))
          debug_cat(paste0("    S1.62dev final data_source in heatmap_item: ",
                          if (!is.null(heatmap_item[[as.character(j)]]$data_source)) heatmap_item[[as.character(j)]]$data_source else "NULL", "\n"))
          # S2.3-DEBUG: Show according columns for verification
          cat(file=stderr(), paste0("[YAML-DEBUG] Heatmap ", j, " YAML according columns: "))
          if (!is.null(heatmap_item[[as.character(j)]]$according) && length(heatmap_item[[as.character(j)]]$according) > 0) {
            for (k_debug in seq_along(heatmap_item[[as.character(j)]]$according)) {
              acc_entry <- heatmap_item[[as.character(j)]]$according[[k_debug]]
              acc_key <- names(acc_entry)[1]
              cat(file=stderr(), paste0(acc_entry[[acc_key]], " "))
            }
            cat(file=stderr(), "\n")
          } else {
            cat(file=stderr(), "(none)\n")
          }
        }
      }

      # v129: Add highlighting support to default classification (was missing!)
      # v130: Fixed - removed values$preview_highlight_active check (was never set)
      # Determine which highlight to apply (same logic as custom classification path)
      highlight_to_apply <- NULL

      # v131: DEBUG - trace highlight decision for default classification
      debug_cat(paste0("\n=== v131: HIGHLIGHT DECISION (DEFAULT classification) ===\n"))
      debug_cat(paste0("  temp_highlight_preview is NULL: ", is.null(values$temp_highlight_preview), "\n"))
      debug_cat(paste0("  active_highlight_index: ", values$active_highlight_index, "\n"))
      debug_cat(paste0("  highlights length: ", length(values$highlights), "\n"))

      # Check if preview mode is active (match custom classification logic at line 8765)
      if (!is.null(values$temp_highlight_preview)) {
        # PREVIEW MODE: Apply temporary highlight for preview
        highlight_to_apply <- values$temp_highlight_preview
      } else if (!is.null(values$active_highlight_index) &&
                 !is.null(values$highlights) &&
                 length(values$highlights) >= values$active_highlight_index) {
        # SAVED MODE: User selected a saved highlight
        highlight_to_apply <- values$highlights[[values$active_highlight_index]]
      }

      # Apply highlight to default classification
      if (!is.null(highlight_to_apply)) {
        debug_cat(paste0("\n=== v129: Adding highlight to DEFAULT classification ===\n"))

        # Build highlight YAML structure
        highlight_yaml <- list(
          display = "yes",
          offset = highlight_to_apply$offset,
          vertical_offset = if(!is.null(highlight_to_apply$vertical_offset)) highlight_to_apply$vertical_offset else 0,
          adjust_height = highlight_to_apply$adjust_height,
          adjust_width = highlight_to_apply$adjust_width,
          according = list()
        )

        # Add each highlighted value
        for (j in seq_along(highlight_to_apply$items)) {
          item <- highlight_to_apply$items[[j]]
          acc_item <- list()
          acc_item[[as.character(j)]] <- list(
            title1 = highlight_to_apply$column,
            value1 = item$value,
            display_name = item$display_name,
            color = item$color,
            transparency = if (!is.null(item$transparency)) item$transparency else 0.5,  # v139: Add transparency
            display_title = highlight_to_apply$title
          )
          highlight_yaml$according <- c(highlight_yaml$according, list(acc_item))
        }

        # ADD TO DEFAULT CLASSIFICATION
        default_classification[["1"]]$highlight <- highlight_yaml
        debug_cat(paste0("  Highlight added with ", length(highlight_to_apply$items), " items\n"))
      } else {
        # No highlight - disable it
        debug_cat("  NO highlight to apply (setting display='no')\n")
        default_classification[["1"]]$highlight <- list(display = "no")
      }

      # S1.62dev: Debug final classification structure before storing
      if (!is.null(default_classification[["1"]]$heatmap_display)) {
        debug_cat(paste0("\n=== S1.62dev: Final heatmap_display in default_classification ===\n"))
        for (h_idx in seq_along(default_classification[["1"]]$heatmap_display)) {
          hm_item <- default_classification[["1"]]$heatmap_display[[h_idx]]
          hm_inner <- hm_item[[as.character(h_idx)]]
          debug_cat(paste0("  Heatmap ", h_idx, " in YAML structure:\n"))
          debug_cat(paste0("    fields: ", paste(names(hm_inner), collapse=", "), "\n"))
          debug_cat(paste0("    data_source: ", if (!is.null(hm_inner$data_source)) hm_inner$data_source else "NULL", "\n"))
        }
      }
      values$yaml_data$`visual definitions`$classification <- list(default_classification)
    }

    # Update font sizes
    # v128: Use Legend tab font size settings for legend_title and legend_text
    # This ensures bootstrap legend matches other legends' font sizes
    values$yaml_data$`visual definitions`$font_size <- list(
      tips = if(!is.null(input$tip_font_size)) input$tip_font_size else 3,
      legend_title = if(!is.null(input$legend_title_size)) input$legend_title_size else 12,
      legend_text = if(!is.null(input$legend_text_size)) input$legend_text_size else 10,
      legend_box = 15,
      heat_map_title = 25,
      heat_map_legend = if(!is.null(input$heatmap_font_size)) input$heatmap_font_size else 3.8
    )

    # v125: Add heatmap global settings (gap between heatmaps)
    values$yaml_data$`visual definitions`$heatmap_global <- list(
      gap = if(!is.null(values$heatmap_global_gap)) values$heatmap_global_gap else 0.05
    )

    # Add compare_two_trees flag (always FALSE for Shiny app)
    values$yaml_data$`visual definitions`$compare_two_trees <- "no"
    
    # Add diagnostics at the end of update_yaml
    # v54: cat("YAML data structure updated. Current configuration:\n")
    # v53: print(utils::str(values$yaml_data))
  }
  
  ###end update yaml
  
  # Match tree tip labels with CSV IDs
  # Match tree tip labels with CSV IDs
  # Match tree tip labels with CSV IDs
  match_tree_with_csv <- function() {
    # S1-PERF: Timing instrumentation to identify bottlenecks
    perf_start_total <- Sys.time()

    # Validate all required inputs before attempting match
    if (is.null(values$tree)) {
      # v53: cat(file=stderr(), "ERROR: Cannot match - No tree loaded\n")
      showNotification("Please upload a tree file first", type = "error")
      return()
    }
    
    if (is.null(values$csv_data)) {
      # v53: cat(file=stderr(), "ERROR: Cannot match - No CSV loaded\n")
      showNotification("Please upload a CSV file first", type = "error")
      return()
    }
    
    if (is.null(input$id_column) || input$id_column == "") {
      # v53: cat(file=stderr(), "ERROR: Cannot match - No ID column selected\n")
      showNotification("Please select an ID column", type = "error")
      return()
    }
    
    # Verify ID column exists in CSV
    if (!input$id_column %in% colnames(values$csv_data)) {
      # v53: cat(file=stderr(), "ERROR: ID column '", input$id_column, "' not found in CSV\n")
      showNotification(paste("ID column", input$id_column, "not found in CSV"), type = "error")
      return()
    }
    
    # Check individual selection if not using all data
    if (!input$use_all_data) {
      if (is.null(input$individual_column) || input$individual_column == "") {
        # v53: cat(file=stderr(), "ERROR: Cannot match - No individual column selected\n")
        showNotification("Please select an Individual column or check 'Use all data'", type = "error")
        return()
      }
      
      if (is.null(input$individual_value) || input$individual_value == "") {
        # v53: cat(file=stderr(), "ERROR: Cannot match - No individual value selected\n")
        showNotification("Please select an Individual value or check 'Use all data'", type = "error")
        return()
      }
      
      # Verify individual column exists in CSV
      if (!input$individual_column %in% colnames(values$csv_data)) {
        # v53: cat(file=stderr(), "ERROR: Individual column '", input$individual_column, "' not found in CSV\n")
        showNotification(paste("Individual column", input$individual_column, "not found in CSV"), type = "error")
        return()
      }
    }
    
    # v53: cat(file=stderr(), "\n=== Starting match_tree_with_csv ===\n")
    # v53: cat(file=stderr(), "All validations passed. Proceeding with matching...\n\n")
    
    req(values$tree, values$csv_data, input$id_column)
    
    # Debug output for function start
    # v53: print("Starting match_tree_with_csv function")
    # v53: print(paste("use_all_data =", input$use_all_data))
    # v53: print(paste("individual_column =", input$individual_column))
    # v53: print(paste("individual_value =", input$individual_value))
    
    # Filter CSV by individual if needed
    filtered_csv <- values$csv_data

    # S1-PERF: Time individual filtering
    perf_filter_start <- Sys.time()

    # Only filter if not using all data and individual column/value are selected
    if (!input$use_all_data && 
        !is.null(input$individual_column) && input$individual_column != "" &&
        !is.null(input$individual_value) && input$individual_value != "") {
      
      # Check if column exists
      if (!(input$individual_column %in% names(values$csv_data))) {
        showNotification(paste("Column not found:", input$individual_column), type = "error")
        return()
      }
      
      tryCatch({
        # Filter to only rows matching the selected individual
        filtered_csv <- values$csv_data[values$csv_data[[input$individual_column]] == input$individual_value, ]
        
        # Debug output
        # v53: print(paste("Filtered by individual:", input$individual_value))
        # v53: print(paste("Filtered CSV rows:", nrow(filtered_csv)))
        
        # Check if we have any rows left
        if (nrow(filtered_csv) == 0) {
          showNotification("No rows match the selected individual", type = "warning")
        }
      }, error = function(e) {
        # v53: print(paste("Error filtering by individual:", e$message))
        showNotification(paste("Error filtering by individual:", e$message), type = "error")
      })
    } else {
      # v53: print("Using all data (no individual filtering)")
    }

    # S1-PERF: End individual filtering timing
    perf_filter_end <- Sys.time()
    cat(file=stderr(), sprintf("[PERF] Individual filtering: %.3f sec\n",
        as.numeric(difftime(perf_filter_end, perf_filter_start, units="secs"))))

    # Store the filtered CSV for use in other functions
    values$filtered_by_individual <- filtered_csv

    # Get tree tips and CSV IDs
    tree_labels <- values$tree$tip.label
    csv_ids <- filtered_csv[[input$id_column]]

    # Debug output
    # v53: print("Tree tip labels (first few):")
    # v53: print(head(tree_labels))
    # v53: print("CSV IDs column:")
    # v53: print(head(csv_ids))

    # S1-PERF: Time ID matching
    perf_match_start <- Sys.time()

    # Attempt to match tree labels with CSV IDs
    match_result <- match_tree_ids_with_csv(tree_labels, csv_ids)
    values$id_match <- match_result

    # S1-PERF: End ID matching timing
    perf_match_end <- Sys.time()
    cat(file=stderr(), sprintf("[PERF] ID matching (match_tree_ids_with_csv): %.3f sec\n",
        as.numeric(difftime(perf_match_end, perf_match_start, units="secs"))))

    # STORE THE INFERRED TRIMMING PARAMETERS
    values$trimming_params <- match_result$trimming_params
    
    # Also update YAML with these parameters
    values$yaml_data$`visual definitions`$id_tip_trim_flag <- 
      if (match_result$trimming_params$id_tip_trim_flag) "yes" else "no"
    values$yaml_data$`visual definitions`$id_tip_trim_start <- 
      match_result$trimming_params$id_tip_trim_start
    values$yaml_data$`visual definitions`$id_tip_trim_end <- 
      match_result$trimming_params$id_tip_trim_end
    values$yaml_data$`visual definitions`$id_tip_prefix <- 
      match_result$trimming_params$id_tip_prefix
    
    # Debug output
    # v53: cat(file=stderr(), "\n=== Inferred Trimming Parameters ===\n")
    # v53: cat(file=stderr(), "id_tip_trim_flag:", match_result$trimming_params$id_tip_trim_flag, "\n")
    # v53: cat(file=stderr(), "id_tip_trim_start:", match_result$trimming_params$id_tip_trim_start, "\n")
    # v53: cat(file=stderr(), "id_tip_trim_end:", match_result$trimming_params$id_tip_trim_end, "\n")
    # v53: cat(file=stderr(), "id_tip_prefix:", match_result$trimming_params$id_tip_prefix, "\n\n")
    
    
    # Debug the match result
    # v53: print("Match result summary:")
    # v53: print(match_result$summary)

    # S1-PERF: Time filtered CSV creation
    perf_filtcsv_start <- Sys.time()

    # Create filtered CSV data with only the matched rows
    if (match_result$summary$exact_matches > 0 || 
        match_result$summary$numeric_matches > 0 || 
        match_result$summary$prefix_suffix_matches > 0) {
      
      # S1-PERF: Get all matched CSV IDs using unlist() instead of loop
      # Growing vectors in loops is O(n²) - unlist is O(n)
      matched_ids <- unique(unlist(match_result$mapping, use.names = FALSE))
      
      # Filter CSV to only include matched rows
      values$filtered_csv <- filtered_csv[filtered_csv[[input$id_column]] %in% matched_ids, ]
      
      # v53: print("Filtered CSV rows:")
      # v53: print(nrow(values$filtered_csv))
      # v53: print("Filtered CSV columns:")
      # v53: print(names(values$filtered_csv))
    } else {
      values$filtered_csv <- NULL
    }

    # S1-PERF: End filtered CSV creation timing
    perf_filtcsv_end <- Sys.time()
    cat(file=stderr(), sprintf("[PERF] Filtered CSV creation: %.3f sec\n",
        as.numeric(difftime(perf_filtcsv_end, perf_filtcsv_start, units="secs"))))

    # ADD THIS DEBUG OUTPUT HERE:
    # v53: cat(file=stderr(), "\n=== AFTER FILTERING SECTION ===\n")
    # v53: cat(file=stderr(), "Made it past the filtering block\n")
    # v53: cat(file=stderr(), "About to check match_success_val\n\n")
    
    # Define match status variables
    #match_success_val <- length(match_result$matched$unmatched) < length(tree_labels)
    #match_warning_val <- length(match_result$matched$unmatched) > 0 && 
    #  length(match_result$matched$unmatched) < length(tree_labels)
    #match_error_val <- length(match_result$matched$unmatched) == length(tree_labels)
    
    #cat(file=stderr(), "\n=== AFTER MATCH STATUS VARIABLES ===\n")
    #cat(file=stderr(), "match_success_val:", match_success_val, "\n")
    #cat(file=stderr(), "match_warning_val:", match_warning_val, "\n")
    #cat(file=stderr(), "match_error_val:", match_error_val, "\n\n")
    
    # Set reactive values directly
    #output$match_success <- renderText({
    #  as.character(match_success_val)
    #})
    
    
    #cat(file=stderr(), "=== AFTER renderText CALLS ===\n")
    #cat(file=stderr(), "Set all renderText outputs\n\n")
    
    #output$match_warning <- renderText({
    #  as.character(match_warning_val)
    #})
    
    #output$match_error <- renderText({
    #  as.character(match_error_val)
    #})
    
    #outputOptions(output, "match_success", suspendWhenHidden = FALSE)
    #outputOptions(output, "match_warning", suspendWhenHidden = FALSE)
    #outputOptions(output, "match_error", suspendWhenHidden = FALSE)
    
    
    # Define match status variables
    match_success_val <- length(match_result$matched$unmatched) < length(tree_labels)
    match_warning_val <- length(match_result$matched$unmatched) > 0 && 
      length(match_result$matched$unmatched) < length(tree_labels)
    match_error_val <- length(match_result$matched$unmatched) == length(tree_labels)
    
    # Store match status in reactive values
    values$match_success <- match_success_val
    values$match_warning <- match_warning_val
    values$match_error <- match_error_val
    
    # v53: cat(file=stderr(), "=== AFTER STORING MATCH STATUS ===\n")
    # v53: cat(file=stderr(), "Stored match status in values\n\n")
    
    # Display unmatched tips
    output$unmatched_tips <- renderPrint({
      cat("Unmatched tree tips:\n")
      # v53: print(match_result$matched$unmatched)
    })
    
    # Create matching table
    match_summary <- data.frame(
      Type = c("Exact Matches", "Numeric Matches", "Prefix/Suffix Matches", "Unmatched", "Total Tree Tips"),
      Count = c(
        match_result$summary$exact_matches,
        match_result$summary$numeric_matches,
        match_result$summary$prefix_suffix_matches,
        match_result$summary$unmatched,
        match_result$summary$total_tree_labels
      )
    )
    
    output$matching_table <- renderDT({
      datatable(match_summary, options = list(dom = 't', ordering = FALSE))
    })
    
    # Update YAML structure with the correct ID column
    values$yaml_data$`Mapping exl renaming titles`$`ID column` <- input$id_column
    
    # v53: cat(file=stderr(), "\n=== BEFORE TEMP CSV CREATION ===\n")
    # v53: cat(file=stderr(), "About to check if we should create temp CSV\n")
    # v53: cat(file=stderr(), "match_success_val:", match_success_val, "\n\n")
    
    
    # If matching is successful, enable other tabs and generate plot
    # If matching is successful, enable other tabs and generate plot
    if (match_success_val) {
      # v53: cat(file=stderr(), "=== INSIDE match_success_val IF BLOCK ===\n")
      
      # Update classification values UI based on matched data
      # v53: cat(file=stderr(), "About to call update_classification_values()\n")
      # Update classification values UI based on matched data
      # TEMPORARILY COMMENT OUT:
      #update_classification_values()
      # v53: cat(file=stderr(), "Finished update_classification_values()\n")
      
      # Create temp CSV file IMMEDIATELY after filtering, BEFORE generate_plot
      # v53: cat(file=stderr(), "About to check filtered_csv for temp file creation\n")
      # v53: cat(file=stderr(), "values$filtered_csv is NULL:", is.null(values$filtered_csv), "\n")
      
      # S1-PERF: Time temp CSV writing
      perf_tempcsv_start <- Sys.time()

      # Create temp CSV file IMMEDIATELY after filtering, BEFORE generate_plot
      if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
        # v53: cat(file=stderr(), "Creating temp CSV file...\n")
        values$temp_csv_path <- tempfile(fileext = ".csv")

        # S1-PERF: Filter out columns with empty/auto-generated names (like "...15355")
        # These are likely empty columns that slow down CSV read/write significantly
        csv_to_write <- values$filtered_csv
        col_names <- names(csv_to_write)
        # Keep only columns with real names (not starting with "..." or empty)
        valid_cols <- !grepl("^\\.\\.\\.", col_names) & col_names != "" & !is.na(col_names)
        if (sum(valid_cols) < length(col_names)) {
          cat(file=stderr(), sprintf("[PERF] Filtering out %d empty/unnamed columns from temp CSV\n",
              sum(!valid_cols)))
          csv_to_write <- csv_to_write[, valid_cols, drop = FALSE]
        }

        write.csv(csv_to_write, values$temp_csv_path, row.names = FALSE)
        # v53: cat(file=stderr(), "Created temp CSV at:", values$temp_csv_path, "\n")
      } else {
        # v53: cat(file=stderr(), "WARNING: No filtered CSV data available for temp file\n")
      }

      # S1-PERF: End temp CSV writing timing
      perf_tempcsv_end <- Sys.time()
      cat(file=stderr(), sprintf("[PERF] Temp CSV write: %.3f sec\n",
          as.numeric(difftime(perf_tempcsv_end, perf_tempcsv_start, units="secs"))))
      # v53: cat(file=stderr(), "About to call generate_plot()\n")
      
      # Debug: Check filtered CSV state
      # v53: cat(file=stderr(), "\nÃ°Å¸â€Â === PRE-GENERATE_PLOT DEBUG ===\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â values$filtered_csv is NULL:", is.null(values$filtered_csv), "\n")
      if (!is.null(values$filtered_csv)) {
        # v53: cat(file=stderr(), "Ã°Å¸â€Â values$filtered_csv rows:", nrow(values$filtered_csv), "\n")
      }
      # v53: cat(file=stderr(), "Ã°Å¸â€Â values$temp_csv_path:", values$temp_csv_path, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â temp CSV exists:", file.exists(values$temp_csv_path), "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â About to call generate_plot()...\n\n")
      
      # S1-PERF: Time generate_plot
      perf_genplot_start <- Sys.time()

      # v53: cat(file=stderr(), "About to call generate_plot()\n")
      # Generate initial plot
      generate_plot()
      # v53: cat(file=stderr(), "Finished generate_plot()\n")

      # S1-PERF: End generate_plot timing
      perf_genplot_end <- Sys.time()
      cat(file=stderr(), sprintf("[PERF] generate_plot(): %.3f sec\n",
          as.numeric(difftime(perf_genplot_end, perf_genplot_start, units="secs"))))

      # S1-PERF: Total time
      perf_end_total <- Sys.time()
      cat(file=stderr(), sprintf("[PERF] === TOTAL Process Data: %.3f sec ===\n",
          as.numeric(difftime(perf_end_total, perf_start_total, units="secs"))))

      # Show notification based on match quality
      if (match_warning_val) {
        showNotification(
          paste("Warning:", match_result$summary$unmatched, "tree tips couldn't be matched."), 
          type = "warning"
        )
      } else {
        showNotification("All tree tips successfully matched to CSV data", type = "message")
      }
    } else if (match_error_val) {
      showNotification("Error: No tree tips could be matched to CSV data", type = "error")
    }
  }
  
  ####end match_tree_with_csv func
  
  # Classification values UI
  update_classification_values <- function() {
    req(values$csv_data, input$classification_column)
    
    # Create color mapping UI
    output$classification_values_ui <- renderUI({
      # CRITICAL: Validate input first
      req(input$classification_column)
      
      # Set loading flag
      classification_loading(TRUE)
      on.exit(classification_loading(FALSE))
      
      # FIXED: Use filtered CSV if available, otherwise use full CSV
      csv_data <- isolate({
        if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
          values$filtered_csv
        } else {
          values$csv_data
        }
      })
      
      column <- isolate(input$classification_column)
      current_class <- isolate(values$current_classification)
      
      # VALIDATE: Make sure column exists in data
      if (is.null(csv_data) || is.null(column) || column == "" || !(column %in% names(csv_data))) {
        # v53: cat(file=stderr(), "Ã¢Å¡Â Ã¯Â¸Â Invalid column selection, skipping UI generation\n")
        return(NULL)
      }
      
      # Get unique values - ISOLATED
      unique_values <- isolate({
        unique(csv_data[[column]])
      })
      unique_values <- unique_values[!is.na(unique_values)]
      
      # v53: cat(file=stderr(), "Creating color UI for column '", column, "' with", length(unique_values), "unique values\n")
      # v53: cat(file=stderr(), "Using", if (!is.null(values$filtered_csv)) "filtered_csv" else "csv_data", "\n")
      
      # Safety check - prevent massive UI generation
      if (length(unique_values) > 100) {
        return(
          tagList(
            tags$div(
              style = "padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
              tags$h5("! Too many unique values!", style = "color: #856404; margin-top: 0;"),
              tags$p(paste("The selected column has", length(unique_values), "unique values, which is too many to display efficiently.")),
              tags$p("Please select a column with fewer unique values (recommended: less than 50).")
            )
          )
        )
      }
      
      # Add palette selection at the top
      palette_ui <- fluidRow(
        column(12,
               tags$div(style = "background-color: #f8f9fa; padding: 10px; margin-bottom: 15px; border-radius: 5px;",
                        tags$h5("Quick Palette Options:", style = "margin-top: 0;"),
                        fluidRow(
                          column(6,
                                 selectInput("color_palette_choice", "Apply Color Palette:",
                                             choices = c(
                                               "Rainbow" = "rainbow",
                                               "Heat Colors" = "heat.colors",
                                               "Terrain Colors" = "terrain.colors",
                                               "Topo Colors" = "topo.colors",
                                               "CM Colors" = "cm.colors",
                                               "Viridis" = "viridis",
                                               "Plasma" = "plasma",
                                               "Inferno" = "inferno",
                                               "Magma" = "magma",
                                               "Cividis" = "cividis"
                                             ),
                                             selected = "rainbow"
                                 )
                          ),
                          column(6, style = "padding-top: 25px;",
                                 actionButton("apply_palette", "Apply Palette to All", 
                                              icon = icon("palette"), class = "btn-primary btn-sm")
                          )
                        )
               )
        )
      )
      
      # R color names list - comprehensive collection
      r_colors <- c(
        "Custom" = "",
        # Basic colors
        "red", "blue", "green", "yellow", "orange", "purple", "pink", "brown",
        "gray", "black", "white", "cyan", "magenta",
        
        # Metallics
        "gold", "silver", "bronze" = "goldenrod4",
        
        # Dark variants
        "darkred", "darkblue", "darkgreen", "darkorange", "darkviolet", 
        "darkgray", "darkcyan", "darkmagenta", "darkgoldenrod",
        
        # Light variants
        "lightblue", "lightgreen", "lightyellow", "lightpink", "lightgray",
        "lightcyan", "lightcoral", "lightsalmon", "lightsteelblue",
        
        # Named colors
        "steelblue", "skyblue", "navy", "maroon", "olive", "teal", "coral",
        "tomato", "salmon", "khaki", "plum", "orchid", "tan", "wheat",
        "sienna", "peru", "chocolate", "saddlebrown", "sandybrown",
        
        # Greens
        "forestgreen", "limegreen", "seagreen", "springgreen", "mediumseagreen",
        "palegreen", "darkseagreen", "olivedrab", "yellowgreen",
        
        # Blues
        "royalblue", "dodgerblue", "deepskyblue", "cornflowerblue", 
        "cadetblue", "midnightblue", "slateblue", "mediumblue",
        
        # Reds/Pinks
        "crimson", "firebrick", "indianred", "lightcoral", "hotpink",
        "deeppink", "palevioletred", "mediumvioletred",
        
        # Purples
        "violet", "indigo", "mediumpurple", "blueviolet", "darkslateblue",
        
        # Oranges/Yellows
        "goldenrod", "darkgoldenrod", "rosybrown", "burlywood",
        "moccasin", "navajowhite", "peachpuff", "palegoldenrod",
        
        # Grays
        "dimgray", "slategray", "lightslategray", "darkslategray",
        "gainsboro", "whitesmoke", "snow", "ivory",
        
        # Aqua/Turquoise
        "aqua", "aquamarine", "turquoise", "darkturquoise", "mediumturquoise",
        "paleturquoise", "lightseagreen", "darkcyan"
      )
      
      # Create individual color pickers for each value
      class_values <- lapply(seq_along(unique_values), function(i) {
        value <- unique_values[i]
        
        fluidRow(
          column(4, tags$b(value)),
          column(4, 
                 colourInput(
                   inputId = paste0("class_color_", i),
                   label = NULL,
                   value = rainbow(length(unique_values))[i],
                   allowTransparent = FALSE
                 )
          ),
          column(4, 
                 selectInput(
                   inputId = paste0("class_color_name_", i),
                   label = NULL,
                   choices = r_colors,
                   selected = ""
                 )
          )
        )
      })
      
      tagList(
        palette_ui,
        tags$h5("Individual Colors:"),
        class_values
      )
    })    
  }
  
  # Observer to update color picker when R color name is selected
  observeEvent(input$classification_column, {
    req(values$csv_data, input$classification_column)
    
    # v52: Use the SAME data source as the UI (filtered_csv if available)
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
      values$filtered_csv
    } else {
      values$csv_data
    }
    
    unique_values <- unique(csv_to_use[[input$classification_column]])
    unique_values <- unique_values[!is.na(unique_values)]
    
    # Create observers for each color name dropdown
    lapply(seq_along(unique_values), function(i) {
      observeEvent(input[[paste0("class_color_name_", i)]], {
        color_name <- input[[paste0("class_color_name_", i)]]
        if (!is.null(color_name) && color_name != "") {
          # Update the color picker with the selected R color
          updateColourInput(session, paste0("class_color_", i), value = color_name)
        }
      }, ignoreInit = TRUE)
    })
    
    update_classification_values()
  }, ignoreInit = TRUE)  # S2.0-PERF: Prevent firing on tab switch
  #}
  
  # Observer for Apply Palette button
  observeEvent(input$apply_palette, {
    req(values$csv_data, input$classification_column, input$color_palette_choice)
    
    # v52: Use the SAME data source as the UI (filtered_csv if available)
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
      values$filtered_csv
    } else {
      values$csv_data
    }
    
    unique_values <- unique(csv_to_use[[input$classification_column]])
    unique_values <- unique_values[!is.na(unique_values)]
    n <- length(unique_values)
    
    # Generate colors based on selected palette with safe error handling
    colors <- tryCatch({
      switch(input$color_palette_choice,
             "rainbow" = rainbow(n),
             "heat.colors" = heat.colors(n),
             "terrain.colors" = terrain.colors(n),
             "topo.colors" = topo.colors(n),
             "cm.colors" = cm.colors(n),
             "viridis" = {
               if (requireNamespace("viridisLite", quietly = TRUE)) {
                 viridisLite::viridis(n)
               } else {
                 # Fallback to rainbow if viridisLite not available
                 rainbow(n)
               }
             },
             "plasma" = {
               if (requireNamespace("viridisLite", quietly = TRUE)) {
                 viridisLite::plasma(n)
               } else {
                 heat.colors(n)
               }
             },
             "inferno" = {
               if (requireNamespace("viridisLite", quietly = TRUE)) {
                 viridisLite::inferno(n)
               } else {
                 terrain.colors(n)
               }
             },
             "magma" = {
               if (requireNamespace("viridisLite", quietly = TRUE)) {
                 viridisLite::magma(n)
               } else {
                 topo.colors(n)
               }
             },
             "cividis" = {
               if (requireNamespace("viridisLite", quietly = TRUE)) {
                 viridisLite::cividis(n)
               } else {
                 cm.colors(n)
               }
             },
             rainbow(n)  # default
      )
    }, error = function(e) {
      # v53: cat(file=stderr(), "Error generating palette:", e$message, "\n")
      rainbow(n)  # Fallback to rainbow on any error
    })
    
    # Update all color pickers
    lapply(seq_along(unique_values), function(i) {
      updateColourInput(session, paste0("class_color_", i), value = colors[i])
    })
    
    showNotification(paste("Applied", input$color_palette_choice, "palette to all colors"), 
                     type = "message", duration = 2)
  })
  
  # Render list of added classifications with radio buttons
  output$classifications_list_ui <- renderUI({
    if (is.null(values$classifications) || length(values$classifications) == 0) {
      return(tags$p(style = "color: gray; font-style: italic;", "No classifications saved yet"))
    }
    
    # Create radio buttons for classification selection
    choices <- setNames(seq_along(values$classifications), 
                        sapply(seq_along(values$classifications), function(i) {
                          class_def <- values$classifications[[i]]
                          num_values <- if (!is.null(class_def$classes)) length(class_def$classes) else 0
                          paste0(class_def$title, " (", class_def$column, ", ", num_values, " values)")
                        }))
    
    # Set selected to current active classification, or last one if not set
    selected_index <- if (!is.null(values$active_classification_index)) {
      values$active_classification_index
    } else {
      length(values$classifications)
    }
    
    tagList(
      tags$div(
        style = "margin-bottom: 10px;",
        radioButtons(
          inputId = "selected_classification_index",
          label = "Active Classification:",
          choices = choices,
          selected = selected_index
        )
      )
    )
  })
  
  # ============================================================================
  # HIGHLIGHTING TAB OBSERVERS
  # ============================================================================
  
  # Update highlight values when column changes
  # ============================================================================
  # HIGHLIGHTING TAB OBSERVERS
  # ============================================================================
  
  # Update highlight values when column changes
  # ============================================================================
  # HIGHLIGHTING TAB OBSERVERS - COMPLETE REPLACEMENT
  # ============================================================================
  
  # Update highlight column dropdown when column changes
  observeEvent(input$highlight_column, {
    if (is.null(input$highlight_column) || input$highlight_column == "") {
      updateSelectizeInput(session, "highlight_values", choices = character(0), selected = character(0))
      output$highlight_values_settings_ui <- renderUI({
        tags$p(style = "color: gray; font-style: italic;", "Please select a column first")
      })
      return(NULL)
    }
    
    req(values$csv_data, input$highlight_column)
    
    if (!(input$highlight_column %in% names(values$csv_data))) {
      return(NULL)
    }
    
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
      values$filtered_csv
    } else {
      values$csv_data
    }
    
    unique_values <- unique(csv_to_use[[input$highlight_column]])
    unique_values <- unique_values[!is.na(unique_values)]
    
    if (length(unique_values) > 100) {
      output$highlight_values_settings_ui <- renderUI({
        tags$div(
          style = "padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
          tags$h5("! Too many unique values!", style = "color: #856404; margin-top: 0;"),
          tags$p(paste("This column has", length(unique_values), "unique values."))
        )
      })
      updateSelectizeInput(session, "highlight_values", choices = character(0))
      return(NULL)
    }
    
    updateSelectizeInput(session, "highlight_values", choices = unique_values, selected = character(0))
    
    output$highlight_values_settings_ui <- renderUI({
      tags$div(
        tags$p(tags$strong(paste("Available values:", length(unique_values)))),
        tags$p("Select values from the dropdown to highlight them."),
        tags$p(style = "color: #666;", "After selecting, you can customize color and settings for each value below.")
      )
    })
  }, ignoreInit = TRUE)  # S2.0-PERF: Prevent firing on tab switch

  # Update when transparency changes
  observeEvent(input$highlight_transparency, {
    req(input$highlight_values)
    if (length(input$highlight_values) > 0) {
      isolate({
        selected_values <- input$highlight_values
        highlight_color <- input$highlight_color
        transparency <- input$highlight_transparency
        
        output$highlight_values_ui <- renderUI({
          tags$div(
            tags$h5(tags$strong(paste("Selected:", length(selected_values), "values"))),
            tags$div(
              style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 15px; background-color: #f8f9fa; max-height: 400px; overflow-y: auto;",
              lapply(selected_values, function(val) {
                tags$div(
                  style = "padding: 8px; margin: 5px 0; border-radius: 3px; background-color: white; display: flex; justify-content: space-between; align-items: center;",
                  tags$span(style = "font-weight: bold; flex-grow: 1;", as.character(val)),
                  tags$div(
                    style = sprintf("width: 60px; height: 30px; border-radius: 3px; border: 2px solid #999; background-color: %s; opacity: %s;", 
                                    highlight_color, transparency)
                  )
                )
              })
            )
          )
        })
      })
    }
  }, ignoreInit = TRUE)
  
  # List of saved highlights
  # List of saved highlights
  output$highlights_list_ui <- renderUI({
    if (is.null(values$highlights) || length(values$highlights) == 0) {
      return(tags$p(style = "color: gray; font-style: italic;", "No highlights saved yet"))
    }
    
    choices <- setNames(seq_along(values$highlights), 
                        sapply(seq_along(values$highlights), function(i) {
                          hl <- values$highlights[[i]]
                          num_items <- length(hl$items)
                          paste0(hl$title, " (", hl$column, ", ", num_items, " values)")
                        }))
    
    selected_index <- if (!is.null(values$active_highlight_index)) {
      values$active_highlight_index
    } else {
      length(values$highlights)
    }
    
    tagList(
      radioButtons(
        inputId = "selected_highlight_index",
        label = "Active Highlight:",
        choices = choices,
        selected = selected_index
      )
    )
  })
  
  # Observer for highlight radio button selection
  observeEvent(input$selected_highlight_index, {
    req(input$selected_highlight_index)
    
    values$active_highlight_index <- as.numeric(input$selected_highlight_index)
    
    # Load the selected highlight settings into the UI
    hl <- values$highlights[[values$active_highlight_index]]
    
    updateSelectInput(session, "highlight_column", selected = hl$column)
    updateSelectizeInput(session, "highlight_values", 
                         selected = sapply(hl$items, function(x) x$value))
    updateTextInput(session, "highlight_title", value = hl$title)
    updateSliderInput(session, "highlight_offset", value = hl$offset)
    updateSliderInput(session, "highlight_vertical_offset", 
                      value = if(!is.null(hl$vertical_offset)) hl$vertical_offset else 0)
    updateSliderInput(session, "highlight_adjust_height", value = hl$adjust_height)
    updateSliderInput(session, "highlight_adjust_width", value = hl$adjust_width)
    
    # Update individual value settings
    lapply(seq_along(hl$items), function(i) {
      item <- hl$items[[i]]
      updateColourInput(session, paste0("highlight_color_", i), value = item$color)
      updateSliderInput(session, paste0("highlight_transparency_", i), 
                        value = if (!is.null(item$transparency)) item$transparency else 0.5)
    })
    
    generate_plot()
  })
  
  # Update Preview button
  # Update Preview button
  # Update Preview button - WITH VISIBLE DEBUGGING
  # Update Preview button - FIXED for Windows
  # Update Preview button - SIMPLIFIED DEBUGGING
  observeEvent(input$update_highlight_preview, {
    cat(file=stderr(), paste0("\n[DEBUG-2ND-HIGHLIGHT] *** UPDATE_HIGHLIGHT_PREVIEW CLICKED at ", format(Sys.time(), "%H:%M:%OS3"), " ***\n"))

    # NOTE: No cooldown check here - user button clicks should never be blocked
    # Cooldown only applies to automatic/cascading triggers

    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   enable_highlight=", input$enable_highlight, "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   highlight_column=", input$highlight_column, "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   highlight_values=", paste(input$highlight_values, collapse=", "), "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   num highlight_values=", length(input$highlight_values), "\n"))

    
    # v53: cat(file=stderr(), "\n\nÃ°Å¸â€Â´Ã°Å¸â€Â´Ã°Å¸â€Â´ UPDATE HIGHLIGHT PREVIEW CLICKED Ã°Å¸â€Â´Ã°Å¸â€Â´Ã°Å¸â€Â´\n")
    
    # VISIBLE NOTIFICATION
    # v53: cat(file=stderr(), "[DEBUG] Preview button clicked\n")
    
    # v53: cat(file=stderr(), "enable_highlight:", input$enable_highlight, "\n")
    # v53: cat(file=stderr(), "highlight_column:", input$highlight_column, "\n")
    # v53: cat(file=stderr(), "highlight_values:", paste(input$highlight_values, collapse=", "), "\n")
    
    # Check if req() conditions are met
    if (!input$enable_highlight) {
      showNotification("Highlighting is not enabled!", type = "error", duration = 5)
      # v53: cat(file=stderr(), "FAILED: enable_highlight is FALSE\n\n")
      return(NULL)
    }
    
    if (is.null(input$highlight_column) || input$highlight_column == "") {
      showNotification("No column selected!", type = "error", duration = 5)
      # v53: cat(file=stderr(), "FAILED: No highlight_column\n\n")
      return(NULL)
    }
    
    if (is.null(input$highlight_values) || length(input$highlight_values) == 0) {
      showNotification("No values selected!", type = "error", duration = 5)
      # v53: cat(file=stderr(), "FAILED: No highlight_values\n\n")
      return(NULL)
    }
    
    # v53: cat(file=stderr(), "Ã¢Å"â€œ All checks passed, collecting settings...\n")
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   All checks passed, collecting highlight_items...\n"))

    # Collect settings for each value
    highlight_items <- lapply(seq_along(input$highlight_values), function(i) {
      val <- input$highlight_values[i]
      color <- input[[paste0("highlight_color_", i)]]
      trans <- input[[paste0("highlight_transparency_", i)]]
      
      if (is.null(color)) color <- rainbow(length(input$highlight_values))[i]
      if (is.null(trans)) trans <- 0.5
      
      # v53: cat(file=stderr(), "  Item", i, "- Value:", val, "Color:", color, "Trans:", trans, "\n")
      
      list(
        column = input$highlight_column,
        value = val,
        display_name = as.character(val),
        color = color,
        transparency = trans
      )
    })
    
    # v53: cat(file=stderr(), "Ã¢Å"â€œ Collected", length(highlight_items), "items\n")
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   Collected ", length(highlight_items), " highlight_items\n"))

    # Store as temporary preview
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   Storing temp_highlight_preview...\n"))
    values$temp_highlight_preview <- list(
      enabled = TRUE,
      title = input$highlight_title,
      column = input$highlight_column,
      offset = input$highlight_offset,
      vertical_offset = input$highlight_vertical_offset,
      adjust_height = input$highlight_adjust_height,
      adjust_width = input$highlight_adjust_width,
      items = highlight_items
    )
    
    # v53: cat(file=stderr(), "Ã¢Å“â€œ Stored temp_highlight_preview\n")
    # v53: cat(file=stderr(), "Ã¢Å“â€œ Calling generate_plot()...\n\n")
    
    showNotification(
      paste0("Generating preview with ", length(highlight_items), " values..."),
      type = "message", 
      duration = 3
    )
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â´ === STORED temp_highlight_preview ===\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â´ Structure:\n")
    # v54: str(values$temp_highlight_preview)
    # v53: cat(file=stderr(), "Ã°Å¸â€Â´ Column:", values$temp_highlight_preview$column, "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â´ Number of items:", length(values$temp_highlight_preview$items), "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â´ =====================================\n\n")
    
    
    # v131: DEBUG - confirm temp_highlight_preview is set before generate_plot()
    debug_cat("\n=== v131: HIGHLIGHT BUTTON - BEFORE generate_plot() ===\n")
    debug_cat(paste0("  temp_highlight_preview is NULL: ", is.null(values$temp_highlight_preview), "\n"))
    if (!is.null(values$temp_highlight_preview)) {
      debug_cat(paste0("  temp_highlight_preview column: ", values$temp_highlight_preview$column, "\n"))
      debug_cat(paste0("  temp_highlight_preview items: ", length(values$temp_highlight_preview$items), "\n"))
    }
    debug_cat("  CALLING generate_plot() NOW...\n")
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   About to call generate_plot() at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
    values$debug_trace_id <- "HIGHLIGHT_BUTTON_PREVIEW"
    generate_plot()
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   generate_plot() returned at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
    values$debug_trace_id <- NULL
    
    # v53: cat(file=stderr(), "Ã¢Å“â€œ generate_plot() completed\n\n")
  })
  
  
  # Save Highlight button
  # Save Highlight button
  # Save Highlight button
  observeEvent(input$save_highlight, {
    req(input$enable_highlight, input$highlight_column, input$highlight_values,
        length(input$highlight_values) > 0)
    
    # Collect settings for each value
    highlight_items <- lapply(seq_along(input$highlight_values), function(i) {
      val <- input$highlight_values[i]
      color <- input[[paste0("highlight_color_", i)]]
      trans <- input[[paste0("highlight_transparency_", i)]]
      
      if (is.null(color)) color <- rainbow(length(input$highlight_values))[i]
      if (is.null(trans)) trans <- 0.5
      
      list(
        column = input$highlight_column,
        value = val,
        display_name = as.character(val),
        color = color,
        transparency = trans
      )
    })
    
    # Create new highlight
    new_highlight <- list(
      enabled = TRUE,
      title = input$highlight_title,
      column = input$highlight_column,
      offset = input$highlight_offset,
      vertical_offset = input$highlight_vertical_offset,
      adjust_height = input$highlight_adjust_height,
      adjust_width = input$highlight_adjust_width,
      items = highlight_items
    )
    
    # Add to highlights list
    if (is.null(values$highlights)) {
      values$highlights <- list()
    }
    values$highlights <- c(values$highlights, list(new_highlight))
    
    # Set as active
    values$active_highlight_index <- length(values$highlights)
    
    # Clear temp preview
    values$temp_highlight_preview <- NULL
    
    generate_plot()
    showNotification("Highlight saved!", type = "message")
  })
  
  # Remove Highlight button
  # Remove Highlight button
  observeEvent(input$remove_highlight, {
    req(values$highlights, length(values$highlights) > 0)
    req(values$active_highlight_index)
    
    index_to_remove <- values$active_highlight_index
    
    if (index_to_remove > 0 && index_to_remove <= length(values$highlights)) {
      values$highlights <- values$highlights[-index_to_remove]
      
      if (length(values$highlights) > 0) {
        values$active_highlight_index <- min(index_to_remove, length(values$highlights))
      } else {
        values$active_highlight_index <- NULL
      }
      
      generate_plot()
      showNotification("Selected highlight removed", type = "warning")
    }
  })
  
  # Update highlight values when column changes
  # Update highlight values when column changes
  # Update highlight values when column changes
  # Update highlight values when column changes
  observeEvent(input$highlight_column, {
    # Skip if empty/null
    if (is.null(input$highlight_column) || input$highlight_column == "") {
      updateSelectizeInput(session, "highlight_values", choices = character(0), selected = character(0))
      output$highlight_values_ui <- renderUI({
        tags$p(style = "color: gray; font-style: italic;", "Please select a column first")
      })
      return(NULL)
    }
    
    req(values$csv_data, input$highlight_column)
    
    # VALIDATE column exists
    if (!(input$highlight_column %in% names(values$csv_data))) {
      return(NULL)
    }
    
    # FIXED: Use filtered CSV if available
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
      values$filtered_csv
    } else {
      values$csv_data
    }
    
    # Get unique values in the selected highlight column
    unique_values <- unique(csv_to_use[[input$highlight_column]])
    unique_values <- unique_values[!is.na(unique_values)]
    
    # Safety check
    if (length(unique_values) > 100) {
      output$highlight_values_ui <- renderUI({
        tags$div(
          style = "padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
          tags$h5("! Too many unique values!", style = "color: #856404; margin-top: 0;"),
          tags$p(paste("This column has", length(unique_values), "unique values. Please select a column with fewer values."))
        )
      })
      updateSelectizeInput(session, "highlight_values", choices = character(0), selected = character(0))
      return(NULL)
    }
    
    # Update highlight values dropdown
    updateSelectizeInput(session, "highlight_values", choices = unique_values, selected = character(0))
    
    # Create initial display (will be updated when values are selected)
    output$highlight_values_ui <- renderUI({
      if (length(unique_values) == 0) {
        return(tags$p(style = "color: gray; font-style: italic;", "No values found in selected column"))
      }
      
      tags$div(
        tags$p(tags$strong(paste("Available values:", length(unique_values)))),
        tags$p("Select values from the dropdown above to highlight them on the tree."),
        tags$p(style = "color: #666; font-size: 0.9em;",
               "All selected values will use the color specified in the 'Highlight Color' picker on the left.")
      )
    })
  }, ignoreInit = TRUE)  # S2.0-PERF: Prevent firing on tab switch

  # Update highlight values settings when values are selected
  observeEvent(input$highlight_values, {
    req(input$highlight_column)
    
    selected_values <- input$highlight_values
    
    if (is.null(selected_values) || length(selected_values) == 0) {
      output$highlight_values_settings_ui <- renderUI({
        tags$div(
          tags$p(tags$strong("No values selected")),
          tags$p("Select values from the dropdown to customize their highlight settings.")
        )
      })
      return(NULL)
    }
    
    # Generate individual settings for each value
    output$highlight_values_settings_ui <- renderUI({
      value_settings <- lapply(seq_along(selected_values), function(i) {
        value <- selected_values[i]
        
        # Get existing settings if this is from a saved highlight
        default_color <- rainbow(length(selected_values))[i]
        default_transparency <- 0.5
        
        tagList(
          tags$hr(style = if (i == 1) "margin-top: 0;" else ""),
          fluidRow(
            column(12, tags$h5(tags$strong(as.character(value))))
          ),
          fluidRow(
            column(4, 
                   colourInput(
                     inputId = paste0("highlight_color_", i),
                     label = "Color",
                     value = default_color,
                     allowTransparent = FALSE
                   )
            ),
            column(4,
                   sliderInput(
                     inputId = paste0("highlight_transparency_", i),
                     label = "Transparency",
                     min = 0.1, max = 1, value = default_transparency, step = 0.1
                   )
            ),
            column(4,
                   # Color preview box
                   tags$div(
                     style = "margin-top: 25px;",
                     tags$div(
                       style = sprintf(
                         "width: 100%%; height: 40px; border-radius: 5px; border: 2px solid #999; background-color: %s; opacity: %s;",
                         default_color, default_transparency
                       ),
                       id = paste0("highlight_preview_box_", i)
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
  
  # v50: Observer to clear highlight when enable_highlight is toggled OFF
  observeEvent(input$enable_highlight, {
    if (!input$enable_highlight) {
      # v53: cat(file=stderr(), "[DEBUG] enable_highlight turned OFF - clearing highlight preview\n")
      values$temp_highlight_preview <- NULL
      generate_plot()
    }
  }, ignoreInit = TRUE)
  
  # v50: Observer to update plot when highlight_values are changed/removed
  observeEvent(input$highlight_values, {
    # If values are cleared while highlighting is enabled, clear the preview
    if (input$enable_highlight && (is.null(input$highlight_values) || length(input$highlight_values) == 0)) {
      # v53: cat(file=stderr(), "[DEBUG] highlight_values cleared - removing highlight preview\n")
      values$temp_highlight_preview <- NULL
      generate_plot()
    }
  }, ignoreInit = TRUE, ignoreNULL = FALSE)
  
  # Update preview boxes when color or transparency changes
  observe({
    req(input$highlight_values)
    
    lapply(seq_along(input$highlight_values), function(i) {
      color_input <- paste0("highlight_color_", i)
      trans_input <- paste0("highlight_transparency_", i)
      
      if (!is.null(input[[color_input]]) && !is.null(input[[trans_input]])) {
        shinyjs::runjs(sprintf(
          "$('#highlight_preview_box_%d').css({'background-color': '%s', 'opacity': '%s'});",
          i, input[[color_input]], input[[trans_input]]
        ))
      }
    })
  })
  
  
  # Update when color changes
  observeEvent(input$highlight_color, {
    req(input$highlight_values)
    if (length(input$highlight_values) > 0) {
      # Trigger re-render
      isolate({
        selected_values <- input$highlight_values
        highlight_color <- input$highlight_color
        transparency <- input$highlight_transparency
        
        output$highlight_values_ui <- renderUI({
          tags$div(
            tags$h5(tags$strong(paste("Selected:", length(selected_values), "values"))),
            tags$div(
              style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 15px; background-color: #f8f9fa; max-height: 400px; overflow-y: auto;",
              lapply(selected_values, function(val) {
                tags$div(
                  style = "padding: 8px; margin: 5px 0; border-radius: 3px; background-color: white; display: flex; justify-content: space-between; align-items: center;",
                  tags$span(style = "font-weight: bold; flex-grow: 1;", as.character(val)),
                  tags$div(
                    style = sprintf("width: 60px; height: 30px; border-radius: 3px; border: 2px solid #999; background-color: %s; opacity: %s;", 
                                    highlight_color, transparency)
                  )
                )
              })
            ),
            tags$hr(),
            tags$p(style = "color: #666; font-size: 0.9em;",
                   icon("info-circle"),
                   " Adjust settings on the left to change appearance.")
          )
        })
      })
    }
  }, ignoreInit = TRUE)
  
  # Update display when transparency changes
  # Update when transparency changes
  observeEvent(input$highlight_transparency, {
    req(input$highlight_values)
    if (length(input$highlight_values) > 0) {
      isolate({
        selected_values <- input$highlight_values
        highlight_color <- input$highlight_color
        transparency <- input$highlight_transparency
        
        output$highlight_values_ui <- renderUI({
          tags$div(
            tags$h5(tags$strong(paste("Selected:", length(selected_values), "values"))),
            tags$div(
              style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 15px; background-color: #f8f9fa; max-height: 400px; overflow-y: auto;",
              lapply(selected_values, function(val) {
                tags$div(
                  style = "padding: 8px; margin: 5px 0; border-radius: 3px; background-color: white; display: flex; justify-content: space-between; align-items: center;",
                  tags$span(style = "font-weight: bold; flex-grow: 1;", as.character(val)),
                  tags$div(
                    style = sprintf("width: 60px; height: 30px; border-radius: 3px; border: 2px solid #999; background-color: %s; opacity: %s;", 
                                    highlight_color, transparency)
                  )
                )
              })
            )
          )
        })
      })
    }
  }, ignoreInit = TRUE)
  
  # v54: Offset adjustment button observers
  observeEvent(input$offset_down_big, {
    current <- if (is.null(input$highlight_offset)) 0 else input$highlight_offset
    new_val <- max(-10, current - 1)
    updateNumericInput(session, "highlight_offset", value = new_val)
  })
  
  observeEvent(input$offset_down_small, {
    current <- if (is.null(input$highlight_offset)) 0 else input$highlight_offset
    new_val <- max(-10, round(current - 0.1, 2))
    updateNumericInput(session, "highlight_offset", value = new_val)
  })
  
  observeEvent(input$offset_down_tiny, {
    current <- if (is.null(input$highlight_offset)) 0 else input$highlight_offset
    new_val <- max(-10, round(current - 0.01, 3))
    updateNumericInput(session, "highlight_offset", value = new_val)
  })
  
  observeEvent(input$offset_up_big, {
    current <- if (is.null(input$highlight_offset)) 0 else input$highlight_offset
    new_val <- min(10, current + 1)
    updateNumericInput(session, "highlight_offset", value = new_val)
  })
  
  observeEvent(input$offset_up_small, {
    current <- if (is.null(input$highlight_offset)) 0 else input$highlight_offset
    new_val <- min(10, round(current + 0.1, 2))
    updateNumericInput(session, "highlight_offset", value = new_val)
  })
  
  observeEvent(input$offset_up_tiny, {
    current <- if (is.null(input$highlight_offset)) 0 else input$highlight_offset
    new_val <- min(10, round(current + 0.01, 3))
    updateNumericInput(session, "highlight_offset", value = new_val)
  })
  
  # v54: Auto-preview when settings change (if checkbox is enabled)
  # Use a reactive timer to debounce rapid changes
  auto_preview_timer <- reactiveVal(NULL)

  observe({
    # Watch these inputs
    input$highlight_offset
    input$highlight_vertical_offset
    input$highlight_adjust_height
    input$highlight_adjust_width

    # Only auto-update if checkbox is checked and we have highlight values
    if (isTRUE(input$auto_preview_highlight) &&
        !is.null(input$highlight_values) &&
        length(input$highlight_values) > 0 &&
        !is.null(values$temp_highlight_preview)) {

      # S2.0-PERF: Skip if a plot was generated too recently (prevents cascading regeneration)
      time_since_last <- as.numeric(Sys.time()) * 1000 - isolate(last_plot_time())
      if (time_since_last < PLOT_COOLDOWN_MS) {
        cat(file=stderr(), sprintf("[PERF] Auto-preview skipped - cooldown active (%.0fms since last plot)\n", time_since_last))
        return()
      }

      # Trigger preview update (same logic as update_highlight_preview button)
      isolate({
        values$temp_highlight_preview$offset <- input$highlight_offset
        values$temp_highlight_preview$vertical_offset <- input$highlight_vertical_offset
        values$temp_highlight_preview$adjust_height <- input$highlight_adjust_height
        values$temp_highlight_preview$adjust_width <- input$highlight_adjust_width

        generate_plot()
      })
    }
  })
  
  
  # Rotation Status Display Box
  output$rotation_status_box <- renderUI({
    req(values$csv_data, input$enable_rotation)
    
    # Build Primary Rotation status
    rotation1_status <- if (!is.null(values$rotation1_config) && length(values$rotation1_config) > 0) {
      col_name <- unique(sapply(values$rotation1_config, function(x) x$col))
      col_name <- col_name[col_name != ""][1]
      val_list <- sapply(values$rotation1_config, function(x) x$val)
      val_list <- val_list[val_list != ""]
      
      if (length(val_list) > 0) {
        tagList(
          tags$p(style = "margin: 5px 0;", 
                 icon("check-circle", style = "color: green;"), 
                 tags$strong(" Configured:"), 
                 sprintf(" %d groups", length(val_list))),
          tags$p(style = "margin: 5px 0; padding-left: 20px;", 
                 tags$strong("Column:"), col_name),
          tags$p(style = "margin: 5px 0; padding-left: 20px;", 
                 tags$strong("Values:"), paste(val_list, collapse = ", "))
        )
      } else {
        tags$p(style = "margin: 5px 0;", 
               icon("times-circle", style = "color: gray;"), 
               " Not configured")
      }
    } else {
      tags$p(style = "margin: 5px 0;", 
             icon("times-circle", style = "color: gray;"), 
             " Not configured")
    }
    
    # Build Secondary Rotation status
    rotation2_status <- if (!is.null(values$rotation2_config) && length(values$rotation2_config) > 0) {
      col_name <- unique(sapply(values$rotation2_config, function(x) x$col))
      col_name <- col_name[col_name != ""][1]
      val_list <- sapply(values$rotation2_config, function(x) x$val)
      val_list <- val_list[val_list != ""]
      
      if (length(val_list) > 0) {
        tagList(
          tags$p(style = "margin: 5px 0;", 
                 icon("check-circle", style = "color: green;"), 
                 tags$strong(" Configured:"), 
                 sprintf(" %d groups", length(val_list))),
          tags$p(style = "margin: 5px 0; padding-left: 20px;", 
                 tags$strong("Column:"), col_name),
          tags$p(style = "margin: 5px 0; padding-left: 20px;", 
                 tags$strong("Values:"), paste(val_list, collapse = ", "))
        )
      } else {
        tags$p(style = "margin: 5px 0;", 
               icon("times-circle", style = "color: gray;"), 
               " Not configured")
      }
    } else {
      tags$p(style = "margin: 5px 0;", 
             icon("times-circle", style = "color: gray;"), 
             " Not configured")
    }
    
    # Build Manual Rotation status
    manual_rotation_status <- if (!is.null(values$manual_rotation_config) && length(values$manual_rotation_config) > 0) {
      node_list <- values$manual_rotation_config
      tagList(
        tags$p(style = "margin: 5px 0;", 
               icon("check-circle", style = "color: green;"), 
               tags$strong(" Configured:"), 
               sprintf(" %d nodes", length(node_list))),
        tags$p(style = "margin: 5px 0; padding-left: 20px;", 
               tags$strong("Nodes:"), paste(node_list, collapse = ", "))
      )
    } else {
      tags$p(style = "margin: 5px 0;", 
             icon("times-circle", style = "color: gray;"), 
             " Not configured")
    }
    
    # Create collapsible box
    box(
      title = tagList(icon("info-circle"), " Rotation Status"),
      status = "info",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = FALSE,
      width = 12,
      tags$div(
        tags$h5(tags$strong("Primary Rotation (rotation1):")),
        rotation1_status,
        tags$hr(),
        tags$h5(tags$strong("Secondary Rotation (rotation2):")),
        rotation2_status,
        tags$hr(),
        tags$h5(tags$strong("Manual Node Rotation:")),
        manual_rotation_status
      )
    )
  })
  
  # Rotation classes UI
  output$rotation_classes_ui <- renderUI({
    req(values$csv_data, input$enable_rotation, input$rotation_type != "manual")
    if (input$rotation_type == "primary") {
      tagList(
        numericInput("rotation1_num_groups", "Number of rotation groups:", value = 2, min = 2, max = 10),
        uiOutput("rotation1_groups_ui"),
        br(),
        actionButton("apply_rotation1", "Apply Rotation", icon = icon("rotate"), class = "btn-primary"),
        actionButton("clear_rotation1", "Clear Primary Configuration", icon = icon("trash"), class = "btn-warning")
      )
    } else if (input$rotation_type == "secondary") {
      tagList(
        numericInput("rotation2_num_groups", "Number of rotation groups:", value = 2, min = 2, max = 10),
        uiOutput("rotation2_groups_ui"),
        br(),
        actionButton("apply_rotation2", "Apply Rotation", icon = icon("rotate"), class = "btn-primary"),
        actionButton("clear_rotation2", "Clear Secondary Configuration", icon = icon("trash"), class = "btn-warning")
      )
    }
  })
  
  output$rotation1_groups_ui <- renderUI({
    req(input$rotation1_num_groups, values$csv_data)
    col_names <- names(values$csv_data)
    group_uis <- lapply(1:input$rotation1_num_groups, function(i) {
      tagList(
        hr(),
        h5(paste("Rotation Group", i)),
        fluidRow(
          column(6, selectInput(inputId = paste0("rotation1_col_", i), label = "Column:", choices = c("-- Select --" = "", col_names))),
          column(6, uiOutput(paste0("rotation1_val_ui_", i)))
        )
      )
    })
    tagList(group_uis)
  })
  
  output$rotation2_groups_ui <- renderUI({
    req(input$rotation2_num_groups, values$csv_data)
    col_names <- names(values$csv_data)
    group_uis <- lapply(1:input$rotation2_num_groups, function(i) {
      tagList(
        hr(),
        h5(paste("Rotation Group", i)),
        fluidRow(
          column(6, selectInput(inputId = paste0("rotation2_col_", i), label = "Column:", choices = c("-- Select --" = "", col_names))),
          column(6, uiOutput(paste0("rotation2_val_ui_", i)))
        )
      )
    })
    tagList(group_uis)
  })
  
  observe({
    req(values$csv_data, input$rotation1_num_groups)
    for (i in 1:input$rotation1_num_groups) {
      local({
        my_i <- i
        output[[paste0("rotation1_val_ui_", my_i)]] <- renderUI({
          col_input <- input[[paste0("rotation1_col_", my_i)]]
          if (is.null(col_input) || col_input == "") {
            selectInput(inputId = paste0("rotation1_val_", my_i), label = "Value:", choices = c("Select column first" = ""))
          } else {
            # FIXED: Use filtered CSV if available
            csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
              values$filtered_csv
            } else {
              values$csv_data
            }
            unique_vals <- unique(csv_to_use[[col_input]])
            unique_vals <- unique_vals[!is.na(unique_vals)]
            selectInput(inputId = paste0("rotation1_val_", my_i), label = "Value:", choices = c("-- Select --" = "", unique_vals))
          }
        })
      })
    }
  })
  
  observe({
    req(values$csv_data, input$rotation2_num_groups)
    for (i in 1:input$rotation2_num_groups) {
      local({
        my_i <- i
        output[[paste0("rotation2_val_ui_", my_i)]] <- renderUI({
          col_input <- input[[paste0("rotation2_col_", my_i)]]
          if (is.null(col_input) || col_input == "") {
            selectInput(inputId = paste0("rotation2_val_", my_i), label = "Value:", choices = c("Select column first" = ""))
          } else {
            # FIXED: Use filtered CSV if available
            csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
              values$filtered_csv
            } else {
              values$csv_data
            }
            unique_vals <- unique(csv_to_use[[col_input]])
            unique_vals <- unique_vals[!is.na(unique_vals)]
            selectInput(inputId = paste0("rotation2_val_", my_i), label = "Value:", choices = c("-- Select --" = "", unique_vals))
          }
        })
      })
    }
  })
  
  observeEvent(input$rotation1_col_1, {
    req(input$rotation1_num_groups)
    first_col <- input$rotation1_col_1
    if (!is.null(first_col) && first_col != "" && input$rotation1_num_groups > 1) {
      for (i in 2:input$rotation1_num_groups) {
        current_val <- input[[paste0("rotation1_col_", i)]]
        if (is.null(current_val) || current_val == "") {
          updateSelectInput(session, paste0("rotation1_col_", i), selected = first_col)
        }
      }
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$rotation2_col_1, {
    req(input$rotation2_num_groups)
    first_col <- input$rotation2_col_1
    if (!is.null(first_col) && first_col != "" && input$rotation2_num_groups > 1) {
      for (i in 2:input$rotation2_num_groups) {
        current_val <- input[[paste0("rotation2_col_", i)]]
        if (is.null(current_val) || current_val == "") {
          updateSelectInput(session, paste0("rotation2_col_", i), selected = first_col)
        }
      }
    }
  }, ignoreInit = TRUE)
  
  # Restore rotation1 configuration when switching to primary
  observeEvent(input$rotation_type, {
    if (!is.null(input$rotation_type) && input$rotation_type == "primary") {
      if (!is.null(values$rotation1_config) && length(values$rotation1_config) > 0) {
        # Restore number of groups
        updateNumericInput(session, "rotation1_num_groups", value = length(values$rotation1_config))
        
        # Wait a bit for UI to render, then restore values
        Sys.sleep(0.1)
        for (i in seq_along(values$rotation1_config)) {
          if (!is.null(values$rotation1_config[[i]])) {
            col <- values$rotation1_config[[i]]$col
            val <- values$rotation1_config[[i]]$val
            if (!is.null(col) && col != "") {
              updateSelectInput(session, paste0("rotation1_col_", i), selected = col)
            }
            if (!is.null(val) && val != "") {
              updateSelectInput(session, paste0("rotation1_val_", i), selected = val)
            }
          }
        }
      }
    }
  }, ignoreInit = TRUE)
  
  # Restore rotation2 configuration when switching to secondary
  observeEvent(input$rotation_type, {
    if (!is.null(input$rotation_type) && input$rotation_type == "secondary") {
      if (!is.null(values$rotation2_config) && length(values$rotation2_config) > 0) {
        # Restore number of groups
        updateNumericInput(session, "rotation2_num_groups", value = length(values$rotation2_config))
        
        # Wait a bit for UI to render, then restore values
        Sys.sleep(0.1)
        for (i in seq_along(values$rotation2_config)) {
          if (!is.null(values$rotation2_config[[i]])) {
            col <- values$rotation2_config[[i]]$col
            val <- values$rotation2_config[[i]]$val
            if (!is.null(col) && col != "") {
              updateSelectInput(session, paste0("rotation2_col_", i), selected = col)
            }
            if (!is.null(val) && val != "") {
              updateSelectInput(session, paste0("rotation2_val_", i), selected = val)
            }
          }
        }
      }
    }
  }, ignoreInit = TRUE)
  
  # Update heatmap custom colors UI
  output$heatmap_custom_colors_ui <- renderUI({
    req(input$heatmap_columns)
    
    num_cols <- length(input$heatmap_columns)
    
    color_inputs <- lapply(1:min(num_cols, 10), function(i) {
      label_text <- ""
      
      if (i <= length(input$heatmap_columns)) {
        label_text <- input$heatmap_columns[i]
      } else {
        label_text <- paste("Color", i)
      }
      
      colourInput(
        inputId = paste0("heatmap_custom_color_", i),
        label = label_text,
        value = rainbow(min(num_cols, 10))[i]
      )
    })
    
    tagList(
      color_inputs
    )
  })
  
  # Heatmap table
  output$heatmap_table <- renderDT({
    req(input$enable_heatmap)
    
    # Get current heatmaps
    heatmaps <- values$heatmaps
    
    if (length(heatmaps) == 0) {
      return(data.frame(
        Title = character(0),
        `Column Count` = integer(0),
        Type = character(0)
      ))
    }
    
    heatmap_data <- lapply(heatmaps, function(hm) {
      type_text <- ""
      if (hm$is_discrete) {
        type_text <- "Discrete"
      } else {
        type_text <- "Continuous"
      }
      
      list(
        Title = hm$title,
        `Column Count` = length(hm$columns),
        Type = type_text
      )
    })
    
    heatmap_df <- do.call(rbind, lapply(heatmap_data, as.data.frame))
    rownames(heatmap_df) <- NULL
    
    datatable(heatmap_df, options = list(dom = 't', ordering = FALSE))
  })
  
  # ============================================
  # v55: NEW MULTI-HEATMAP SYSTEM
  # ============================================
  
  # Render heatmap cards UI
  output$heatmap_cards_ui <- renderUI({
    # v107: Use trigger to control when UI regenerates
    # This prevents rebuilding UI on every property change (slider, color, etc.)
    heatmap_ui_trigger()  # Take dependency on trigger only

    # Use isolate() to get values without creating reactive dependencies
    configs <- isolate(values$heatmap_configs)
    csv_data_local <- isolate(values$csv_data)

    if (length(configs) == 0) {
      return(tags$div(
        class = "alert alert-info",
        icon("info-circle"),
        " No heatmaps configured. Click 'Add New Heatmap' to create one."
      ))
    }

    # Build cards for each heatmap
    cards <- lapply(seq_along(configs), function(i) {
      cfg <- configs[[i]]
      card_id <- paste0("heatmap_card_", i)

      # Get column choices from CSV (using isolated value)
      col_choices <- if (!is.null(csv_data_local)) names(csv_data_local) else character(0)

      # S2.5-FIX: Reorder choices to put selected columns first (in their original order)
      # This ensures selectize preserves the user's column selection order when UI rebuilds
      # Without this, adding a second heatmap causes the first heatmap's columns to reorder
      # (selectize displays items in the order they appear in choices)
      if (!is.null(cfg$columns) && length(cfg$columns) > 0 && length(col_choices) > 0) {
        # Get columns that are both selected AND exist in current choices
        valid_selected <- cfg$columns[cfg$columns %in% col_choices]
        # Get remaining columns not in selection
        remaining_choices <- setdiff(col_choices, valid_selected)
        # Put selected first (preserving their order), then remaining alphabetically
        col_choices <- c(valid_selected, sort(remaining_choices))
      }
      
      # v56/v110/v119: Determine detected type based on first column (if multiple, they should be same type)
      # v119: Improved detection - also try to convert string columns to numeric
      detected_type <- "unknown"
      if (!is.null(cfg$columns) && length(cfg$columns) > 0 && !is.null(csv_data_local)) {
        # Check first column to determine type
        first_col <- cfg$columns[1]
        if (first_col %in% names(csv_data_local)) {
          col_data <- csv_data_local[[first_col]]
          if (!is.null(col_data)) {
            non_na_vals <- na.omit(col_data)
            unique_vals <- length(unique(non_na_vals))
            is_numeric <- is.numeric(col_data)

            # v119: If not already numeric, try to convert (handles "1", "2", "3.5" etc.)
            if (!is_numeric && length(non_na_vals) > 0) {
              # Try converting to numeric
              converted <- suppressWarnings(as.numeric(as.character(non_na_vals)))
              # Check if most values converted successfully (at least 80%)
              conversion_rate <- sum(!is.na(converted)) / length(non_na_vals)
              if (conversion_rate >= 0.8) {
                is_numeric <- TRUE
                non_na_vals <- converted[!is.na(converted)]
                unique_vals <- length(unique(non_na_vals))
              }
            }

            # v124: Harmonized detection with apply_heatmaps section
            # Logic: decimals = continuous, otherwise check unique values and range
            if (is_numeric && length(non_na_vals) > 0) {
              # Check if values have decimals (not all integers)
              has_decimals <- any(non_na_vals != floor(non_na_vals))
              # Calculate value range for additional check
              val_range <- diff(range(non_na_vals, na.rm = TRUE))

              if (has_decimals) {
                # Decimals ALWAYS mean continuous (measurements, percentages)
                detected_type <- "continuous"
              } else {
                # v124: Match apply section logic - boolean-like or small categorical = discrete
                is_boolean_like <- unique_vals <= 2 && val_range <= 1
                is_small_categorical <- unique_vals <= 3 && val_range <= 2
                if (is_boolean_like || is_small_categorical) {
                  detected_type <- "discrete"
                } else if (unique_vals >= 5 || val_range >= 5) {
                  # v124: >=5 unique values OR range>=5 = continuous
                  detected_type <- "continuous"
                } else {
                  detected_type <- "discrete"
                }
              }
            } else {
              detected_type <- "discrete"
            }
          }
        }
      }
      
      # Discrete palette options
      discrete_palettes <- c("Set1", "Set2", "Set3", "Paired", "Dark2", "Accent", "Pastel1", "Pastel2")
      
      # Continuous palette options
      continuous_palettes <- c("Blues", "Greens", "Reds", "Purples", "Oranges", 
                               "Viridis", "Plasma", "Inferno", "Magma",
                               "RdBu", "RdYlGn", "PiYG", "BrBG")
      
      tags$div(
        id = card_id,
        class = "well",
        style = "margin-bottom: 15px; padding: 15px; border-left: 4px solid #3c8dbc;",
        
        # Header with title and action buttons
        fluidRow(
          column(6,
                 tags$h4(
                   style = "margin-top: 0;",
                   icon("layer-group"),
                   paste(" Heatmap", i),
                   # v56: Show column count instead of single column name
                   if (!is.null(cfg$columns) && length(cfg$columns) > 0) {
                     tags$small(class = "text-muted", paste0(" (", length(cfg$columns), " columns)"))
                   }
                 )
          ),
          column(6, style = "text-align: right;",
                 # Move up/down buttons
                 if (i > 1) {
                   actionButton(paste0("heatmap_up_", i), "", icon = icon("arrow-up"), 
                                class = "btn-sm btn-default", style = "margin-right: 5px;")
                 },
                 if (i < length(configs)) {
                   actionButton(paste0("heatmap_down_", i), "", icon = icon("arrow-down"), 
                                class = "btn-sm btn-default", style = "margin-right: 5px;")
                 },
                 actionButton(paste0("heatmap_remove_", i), "", icon = icon("trash"), 
                              class = "btn-sm btn-danger")
          )
        ),
        
        hr(style = "margin: 10px 0;"),

        # S1.62dev: Data source selector (CSV or RData CNV)
        fluidRow(
          column(12,
                 radioButtons(paste0("heatmap_data_source_", i), "Data Source:",
                              choices = c("CSV Columns" = "csv", "RData CNV" = "rdata"),
                              selected = if (!is.null(cfg$data_source)) cfg$data_source else "csv",
                              inline = TRUE)
          )
        ),

        # S1.62dev: Conditional panel for CSV columns (default)
        conditionalPanel(
          condition = paste0("input.heatmap_data_source_", i, " == 'csv'"),
          # Main configuration
          fluidRow(
            column(8,
                   # v56: Changed to selectizeInput with multiple=TRUE for multiple columns
                   selectizeInput(paste0("heatmap_columns_", i), "Data Columns (select one or more)",
                                  choices = col_choices,
                                  selected = if (!is.null(cfg$columns)) cfg$columns else NULL,
                                  multiple = TRUE,
                                  options = list(placeholder = "Select columns..."))
            ),
            column(4,
                   textInput(paste0("heatmap_title_", i), "Legend Title",
                             value = if (!is.null(cfg$title)) cfg$title else paste0("Heatmap ", i))
            )
          )
        ),

        # S1.62dev: Conditional panel for RData CNV
        conditionalPanel(
          condition = paste0("input.heatmap_data_source_", i, " == 'rdata'"),
          fluidRow(
            column(4,
                   textInput(paste0("heatmap_title_rdata_", i), "Legend Title",
                             value = if (!is.null(cfg$title)) cfg$title else "CNV")
            ),
            column(4,
                   sliderInput(paste0("heatmap_cnv_downsample_", i), "Downsample Factor",
                               min = 1, max = 50, value = if (!is.null(cfg$cnv_downsample)) cfg$cnv_downsample else 10,
                               step = 1)
            ),
            column(4,
                   checkboxInput(paste0("heatmap_cnv_wgd_norm_", i), "WGD Normalization (All Cells)",
                                 value = if (!is.null(cfg$cnv_wgd_norm)) cfg$cnv_wgd_norm else FALSE)
            )
          ),
          # S2.12: Per-cell WGD normalization option
          fluidRow(
            column(4,
                   checkboxInput(paste0("heatmap_cnv_wgd_per_cell_", i), "Per-cell WGD (from CSV column)",
                                 value = if (!is.null(cfg$cnv_wgd_per_cell)) cfg$cnv_wgd_per_cell else FALSE)
            ),
            column(8,
                   conditionalPanel(
                     condition = paste0("input.heatmap_cnv_wgd_per_cell_", i),
                     selectizeInput(paste0("heatmap_cnv_wgd_column_", i), "WGD Column",
                                    choices = if (!is.null(values$csv_data)) {
                                      csv_cols <- names(values$csv_data)
                                      # Auto-detect WGD column
                                      wgd_cols <- csv_cols[grepl("^WGD$", csv_cols, ignore.case = TRUE)]
                                      if (length(wgd_cols) > 0) {
                                        c(csv_cols)  # All columns available
                                      } else {
                                        csv_cols
                                      }
                                    } else character(0),
                                    selected = if (!is.null(cfg$cnv_wgd_column)) {
                                      cfg$cnv_wgd_column
                                    } else if (!is.null(values$csv_data)) {
                                      # Auto-select WGD column if it exists
                                      csv_cols <- names(values$csv_data)
                                      wgd_match <- csv_cols[grepl("^WGD$", csv_cols, ignore.case = TRUE)]
                                      if (length(wgd_match) > 0) wgd_match[1] else NULL
                                    } else NULL,
                                    options = list(placeholder = "Select WGD column..."))
                   )
            )
          ),
          tags$div(
            style = "padding: 0 15px; margin-bottom: 10px;",
            tags$small(class = "text-muted",
                       icon("info-circle"),
                       " Per-cell WGD: normalize only cells where selected column value = 1. ",
                       "Auto-detects column named 'WGD'.")
          ),
          # S2.8: Display mode selector (basic vs detailed)
          fluidRow(
            column(6,
                   radioButtons(paste0("heatmap_cnv_display_mode_", i), "Display Mode:",
                                choices = c("Basic (faster)" = "basic", "Detailed (like pheatmap)" = "detailed"),
                                selected = if (!is.null(cfg$cnv_display_mode)) cfg$cnv_display_mode else "basic",
                                inline = TRUE)
            ),
            column(6,
                   tags$div(
                     style = "padding-top: 25px;",
                     tags$small(class = "text-muted",
                                icon("info-circle"),
                                " Detailed mode: smoother rendering with finer color gradations")
                   )
            )
          ),
          # Render downsample factor (second stage downsampling)
          fluidRow(
            column(6,
                   numericInput(paste0("heatmap_cnv_render_downsample_", i),
                               "Render Downsample Factor",
                               value = if (!is.null(cfg$cnv_render_downsample)) cfg$cnv_render_downsample else 10,
                               min = 0, max = 100, step = 1)
            ),
            column(6,
                   tags$div(
                     style = "padding-top: 25px;",
                     tags$small(class = "text-muted",
                                icon("info-circle"),
                                " 0 = keep all, N > 1 = keep every Nth position during render")
                   )
            )
          ),
          # Height scale for detailed mode (allows compressing the heatmap vertically)
          fluidRow(
            column(6,
                   sliderInput(paste0("heatmap_cnv_height_scale_", i),
                               "Height Scale (detailed mode)",
                               min = 0.1, max = 2.0,
                               value = if (!is.null(cfg$cnv_height_scale)) cfg$cnv_height_scale else 1.0,
                               step = 0.1)
            ),
            column(6,
                   tags$div(
                     style = "padding-top: 25px;",
                     tags$small(class = "text-muted",
                                icon("info-circle"),
                                " Controls total heatmap height. Lower values compress the heatmap.")
                   )
            )
          ),
          # S2.0: Sample mapping column selector (shown when auto-match fails)
          fluidRow(
            column(12,
                   uiOutput(paste0("heatmap_rdata_mapping_ui_", i))
            )
          ),
          fluidRow(
            column(12,
                   tags$div(
                     style = "background: #d1ecf1; padding: 10px; border-radius: 5px; margin-top: 5px;",
                     tags$small(
                       icon("info-circle"),
                       " CNV data will be displayed from the loaded RData file. ",
                       "Default color scale: Red (loss) - White (neutral) - Blue (gain)."
                     )
                   )
            )
          )
        ),

        # S1.62dev: CSV-only options (column range and detected type)
        conditionalPanel(
          condition = paste0("input.heatmap_data_source_", i, " == 'csv'"),
          # v121: Column range selector - allows quick selection of contiguous columns
          fluidRow(
            column(6,
                   tags$div(
                     style = "display: flex; align-items: flex-end; gap: 10px;",
                     tags$div(
                       style = "flex: 1;",
                       textInput(paste0("heatmap_col_range_", i), "Column Range (e.g., 2-10)",
                                 value = "",
                                 placeholder = "2-10 or 3-15")
                     ),
                     tags$div(
                       style = "padding-bottom: 15px;",
                       actionButton(paste0("heatmap_add_range_", i), "Add Range",
                                    class = "btn-sm btn-info",
                                    icon = icon("plus"))
                     )
                   )
            ),
            column(6,
                   tags$div(
                     style = "padding-top: 25px;",
                     tags$small(class = "text-muted",
                                icon("info-circle"),
                                " Enter column numbers (e.g., '2-10') to add columns by position")
                   )
            )
          ),

          # v127: Dynamic detected type display - updates when columns change
          fluidRow(
            column(4,
                   tags$label("Detected Type"),
                   uiOutput(paste0("heatmap_detected_type_display_", i))
            ),
            # v111: Removed "Data Columns Count" - not useful to the user
          column(4),
          column(4)
          )
        ),  # End CSV conditionalPanel

        # Type override and settings
        fluidRow(
          column(4,
                 checkboxInput(paste0("heatmap_auto_type_", i), "Auto-detect type", 
                               value = if (!is.null(cfg$auto_type)) cfg$auto_type else TRUE)
          ),
          column(4,
                 conditionalPanel(
                   condition = paste0("!input.heatmap_auto_type_", i),
                   radioButtons(paste0("heatmap_type_", i), "Force Type:",
                                choices = c("Discrete" = "discrete", "Continuous" = "continuous"),
                                selected = if (!is.null(cfg$type)) cfg$type else "discrete",
                                inline = TRUE)
                 )
          ),
          column(4,
                 sliderInput(paste0("heatmap_colnames_angle_", i), "Column name angle",
                             min = 0, max = 90,
                             value = if (!is.null(cfg$colnames_angle)) cfg$colnames_angle else 45,
                             step = 15)
          )
        ),

        # v105: Per-heatmap distance and height sliders
        fluidRow(
          column(4,
                 sliderInput(paste0("heatmap_distance_", i), "Distance from Tree",
                             min = 0, max = 1.0,
                             value = if (!is.null(cfg$distance)) cfg$distance else 0.02,
                             step = 0.01)
          ),
          column(4,
                 sliderInput(paste0("heatmap_height_", i), "Row Height",
                             min = 0.01, max = 3.0,
                             value = if (!is.null(cfg$height)) cfg$height else 0.8,
                             step = 0.01)
          ),
          column(4,
                 sliderInput(paste0("heatmap_row_height_", i), "Column Width",
                             min = 0.5, max = 3.0,
                             value = if (!is.null(cfg$row_height)) cfg$row_height else 1.0,
                             step = 0.1)
          )
        ),

        # v111: Grid options for heatmap squares
        fluidRow(
          column(4,
                 checkboxInput(paste0("heatmap_show_grid_", i), "Show grid around tiles",
                               value = if (!is.null(cfg$show_grid)) cfg$show_grid else FALSE)
          ),
          column(4,
                 colourInput(paste0("heatmap_grid_color_", i), "Grid color",
                             value = if (!is.null(cfg$grid_color)) cfg$grid_color else "#000000",
                             showColour = "background")
          ),
          column(4,
                 sliderInput(paste0("heatmap_grid_size_", i), "Grid line width",
                             min = 0.1, max = 2.0,
                             value = if (!is.null(cfg$grid_size)) cfg$grid_size else 0.5,
                             step = 0.1)
          )
        ),
        # S1.62dev: Horizontal row lines option
        fluidRow(
          column(4,
                 checkboxInput(paste0("heatmap_show_row_lines_", i), "Show horizontal row lines",
                               value = if (!is.null(cfg$show_row_lines)) cfg$show_row_lines else FALSE)
          ),
          column(4,
                 colourInput(paste0("heatmap_row_line_color_", i), "Row line color",
                             value = if (!is.null(cfg$row_line_color)) cfg$row_line_color else "#000000",
                             showColour = "background")
          ),
          column(4,
                 sliderInput(paste0("heatmap_row_line_size_", i), "Row line width",
                             min = 0.1, max = 2.0,
                             value = if (!is.null(cfg$row_line_size)) cfg$row_line_size else 0.5,
                             step = 0.1)
          )
        ),
        # S1.62dev: Vertical column lines option
        fluidRow(
          column(4,
                 checkboxInput(paste0("heatmap_show_col_lines_", i), "Show vertical column lines",
                               value = if (!is.null(cfg$show_col_lines)) cfg$show_col_lines else FALSE)
          ),
          column(4,
                 colourInput(paste0("heatmap_col_line_color_", i), "Column line color",
                             value = if (!is.null(cfg$col_line_color)) cfg$col_line_color else "#000000",
                             showColour = "background")
          ),
          column(4,
                 sliderInput(paste0("heatmap_col_line_size_", i), "Column line width",
                             min = 0.1, max = 2.0,
                             value = if (!is.null(cfg$col_line_size)) cfg$col_line_size else 0.5,
                             step = 0.1)
          )
        ),

        # S1.62dev: Vertical text labels below heatmap
        tags$div(
          style = "background-color: #f5f0e8; padding: 10px; border-radius: 5px; margin-top: 10px; margin-bottom: 10px;",
          tags$h5(icon("font"), " Vertical Text Labels (below heatmap)"),
          fluidRow(
            column(4,
                   checkboxInput(paste0("heatmap_show_vertical_text_", i), "Show vertical text labels",
                                 value = if (!is.null(cfg$show_vertical_text)) cfg$show_vertical_text else FALSE)
            ),
            column(4,
                   selectizeInput(paste0("heatmap_vertical_text_column_", i), "Text from CSV column:",
                                  choices = col_choices,
                                  selected = if (!is.null(cfg$vertical_text_column)) cfg$vertical_text_column else NULL,
                                  options = list(placeholder = "Select column"))
            ),
            column(4,
                   sliderInput(paste0("heatmap_vertical_text_size_", i), "Text size",
                               min = 1, max = 10,
                               value = if (!is.null(cfg$vertical_text_size)) cfg$vertical_text_size else 3,
                               step = 0.5)
            )
          ),
          fluidRow(
            column(4,
                   numericInput(paste0("heatmap_vertical_text_offset_", i), "Vertical offset",
                                value = if (!is.null(cfg$vertical_text_offset)) cfg$vertical_text_offset else 0.5,
                                min = -5, max = 10, step = 0.1)
            ),
            column(4,
                   colourInput(paste0("heatmap_vertical_text_color_", i), "Text color",
                               value = if (!is.null(cfg$vertical_text_color)) cfg$vertical_text_color else "#000000",
                               showColour = "background")
            ),
            column(4)
          )
        ),

        # v116: Tip Guide Lines - vertical lines from tips through heatmap
        tags$div(
          style = "background-color: #e8f4f8; padding: 10px; border-radius: 5px; margin-top: 10px; margin-bottom: 10px;",
          tags$h5(icon("grip-lines-vertical"), " Tip Guide Lines"),
          fluidRow(
            column(4,
                   checkboxInput(paste0("heatmap_show_guides_", i), "Show vertical guide lines",
                                 value = if (!is.null(cfg$show_guides)) cfg$show_guides else FALSE)
            ),
            column(4,
                   colourInput(paste0("heatmap_guide_color1_", i), "Guide color 1",
                               value = if (!is.null(cfg$guide_color1)) cfg$guide_color1 else "#CCCCCC",
                               showColour = "background")
            ),
            column(4,
                   colourInput(paste0("heatmap_guide_color2_", i), "Guide color 2",
                               value = if (!is.null(cfg$guide_color2)) cfg$guide_color2 else "#EEEEEE",
                               showColour = "background")
            )
          ),
          fluidRow(
            column(4,
                   sliderInput(paste0("heatmap_guide_alpha_", i), "Guide transparency",
                               min = 0.05, max = 1.0,
                               value = if (!is.null(cfg$guide_alpha)) cfg$guide_alpha else 0.3,
                               step = 0.05)
            ),
            column(4,
                   sliderInput(paste0("heatmap_guide_width_", i), "Guide line width",
                               min = 0.1, max = 10.0,
                               value = if (!is.null(cfg$guide_width)) cfg$guide_width else 0.5,
                               step = 0.1)
            ),
            column(4,
                   selectInput(paste0("heatmap_guide_linetype_", i), "Line type",
                               choices = c("solid" = "solid",
                                          "dashed" = "dashed",
                                          "dotted" = "dotted",
                                          "dotdash" = "dotdash",
                                          "longdash" = "longdash",
                                          "twodash" = "twodash"),
                               selected = if (!is.null(cfg$guide_linetype)) cfg$guide_linetype else "solid")
            )
          )
        ),

        # v105/v108: Row labels settings with per-column mapping
        tags$div(
          style = "background-color: #fff9e6; padding: 10px; border-radius: 5px; margin-top: 10px; margin-bottom: 10px;",
          tags$h5(icon("font"), " Row Labels (next to heatmap)"),
          fluidRow(
            column(4,
                   checkboxInput(paste0("heatmap_show_row_labels_", i), "Show row labels",
                                 value = if (!is.null(cfg$show_row_labels)) cfg$show_row_labels else FALSE)
            ),
            column(4,
                   selectInput(paste0("heatmap_row_label_source_", i), "Label source",
                               choices = c("Column names" = "colnames",
                                           "Custom mapping" = "mapping",
                                           "Comma-separated list" = "custom"),
                               selected = if (!is.null(cfg$row_label_source)) cfg$row_label_source else "colnames")
            ),
            column(4,
                   sliderInput(paste0("heatmap_row_label_font_size_", i), "Label font size",
                               min = 1, max = 12,  # v109: Increased max from 8 to 12
                               value = if (!is.null(cfg$row_label_font_size)) cfg$row_label_font_size else 2.5,
                               step = 0.5)
            )
          ),
          # v111: Row label offset and alignment options
          fluidRow(
            column(4,
                   sliderInput(paste0("heatmap_row_label_offset_", i), "Label offset from heatmap",
                               min = -8, max = 8,
                               value = if (!is.null(cfg$row_label_offset)) cfg$row_label_offset else 1.0,
                               step = 0.5)
            ),
            column(4,
                   selectInput(paste0("heatmap_row_label_align_", i), "Label alignment",
                               choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                               selected = if (!is.null(cfg$row_label_align)) cfg$row_label_align else "left")
            ),
            column(4)
          ),
          # v108: Dynamic UI for custom labels - either mapping or comma-separated
          uiOutput(paste0("heatmap_custom_labels_ui_", i))
        ),

        # v108: Replaced conditionalPanels with uiOutput for reactive type detection
        # This ensures the color settings show up immediately when columns are selected
        uiOutput(paste0("heatmap_type_settings_ui_", i))
      )
    })
    
    do.call(tagList, cards)
  })
  
  # Add new heatmap
  observeEvent(input$add_new_heatmap, {
    # v141: Increased max heatmaps from 6 to 10
    if (length(values$heatmap_configs) >= 10) {
      showNotification("Maximum 10 heatmaps allowed", type = "warning")
      return()
    }

    # S2.6-DEBUG: Log existing custom_colors BEFORE adding new heatmap
    cat(file=stderr(), paste0("\n[S2.6-DEBUG] ===== ADD NEW HEATMAP CLICKED =====\n"))
    cat(file=stderr(), paste0("[S2.6-DEBUG] Current number of heatmap_configs: ", length(values$heatmap_configs), "\n"))
    for (debug_idx in seq_along(values$heatmap_configs)) {
      cfg <- values$heatmap_configs[[debug_idx]]
      cat(file=stderr(), paste0("[S2.6-DEBUG] Heatmap ", debug_idx, " custom_colors:\n"))
      if (!is.null(cfg$custom_colors) && length(cfg$custom_colors) > 0) {
        cat(file=stderr(), paste0("[S2.6-DEBUG]   keys: ", paste(names(cfg$custom_colors), collapse=", "), "\n"))
        cat(file=stderr(), paste0("[S2.6-DEBUG]   values: ", paste(unlist(cfg$custom_colors), collapse=", "), "\n"))
      } else {
        cat(file=stderr(), paste0("[S2.6-DEBUG]   (no custom colors)\n"))
      }
    }
    cat(file=stderr(), paste0("[S2.6-DEBUG] ================================\n"))

    # v56: Add new empty config with columns (plural) for multiple column support
    # v107: Added distance and height initialization to prevent reactive loops
    new_config <- list(
      # S1.62dev: Data source selector (csv or rdata)
      data_source = "csv",
      columns = character(0),  # v56: Changed from column to columns
      title = paste0("Heatmap ", length(values$heatmap_configs) + 1),
      auto_type = TRUE,
      type = "discrete",
      show_colnames = TRUE,
      colnames_angle = 45,
      distance = 0.02,       # v107: Initialize distance
      height = 0.8,          # v107: Initialize height
      show_row_labels = FALSE,         # v107: Initialize row label settings
      row_label_source = "colnames",
      row_label_font_size = 2.5,
      custom_row_labels = "",
      discrete_palette = "Set1",
      custom_discrete = FALSE,
      custom_colors = list(),
      cont_palette = "Blues",
      low_color = "#FFFFCC",
      high_color = "#006837",
      use_midpoint = FALSE,
      mid_color = "#FFFF99",
      midpoint = 0,
      # S1.62dev: CNV-specific settings
      cnv_downsample = 10,  # Kept for YAML compatibility
      cnv_render_downsample = 10,  # Second stage downsampling (10 = default)
      cnv_height_scale = 1.0,  # Height scaling for detailed mode
      cnv_wgd_norm = FALSE,
      # S2.12: Per-cell WGD normalization settings
      cnv_wgd_per_cell = FALSE,
      cnv_wgd_column = NULL,
      # S2.8: Display mode (basic = geom_tile, detailed = geom_raster like pheatmap)
      cnv_display_mode = "basic"
    )
    
    values$heatmap_configs <- c(values$heatmap_configs, list(new_config))

    # S2.6-DEBUG: Log existing custom_colors AFTER adding new heatmap
    cat(file=stderr(), paste0("\n[S2.6-DEBUG] After adding heatmap, configs count: ", length(values$heatmap_configs), "\n"))
    for (debug_idx in seq_along(values$heatmap_configs)) {
      cfg <- values$heatmap_configs[[debug_idx]]
      cat(file=stderr(), paste0("[S2.6-DEBUG] Heatmap ", debug_idx, " custom_colors AFTER:\n"))
      if (!is.null(cfg$custom_colors) && length(cfg$custom_colors) > 0) {
        cat(file=stderr(), paste0("[S2.6-DEBUG]   keys: ", paste(names(cfg$custom_colors), collapse=", "), "\n"))
        cat(file=stderr(), paste0("[S2.6-DEBUG]   values: ", paste(unlist(cfg$custom_colors), collapse=", "), "\n"))
      } else {
        cat(file=stderr(), paste0("[S2.6-DEBUG]   (no custom colors)\n"))
      }
    }

    # S2.8-FIX: Set inhibit flag BEFORE triggering UI rebuild
    # This prevents color picker observers from overwriting custom colors with defaults
    inhibit_color_save(TRUE)
    cat(file=stderr(), "[S2.8-DEBUG] inhibit_color_save set to TRUE (add heatmap)\n")

    # v107: Trigger UI regeneration when heatmap is added
    heatmap_ui_trigger(heatmap_ui_trigger() + 1)

    # S2.9-FIX: Reset inhibit flag after UI has had time to rebuild
    # Increased delay from 500ms to 2500ms because Shiny's reactive system
    # takes time to process heatmap_ui_trigger, then renderUI runs, then
    # colourInput widgets are created, and finally change events fire.
    # The inhibit must stay TRUE until ALL of these steps complete.
    shinyjs::delay(2500, {
      inhibit_color_save(FALSE)
      cat(file=stderr(), "[S2.9-DEBUG] inhibit_color_save set to FALSE (after 2500ms delay)\n")
    })

    showNotification(paste("Heatmap", length(values$heatmap_configs), "added"), type = "message")
  })

  # Generic observer for heatmap removal buttons
  observe({
    lapply(1:6, function(i) {
      observeEvent(input[[paste0("heatmap_remove_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs <- values$heatmap_configs[-i]

          # S2.9-FIX: Set inhibit flag before UI rebuild
          inhibit_color_save(TRUE)

          # v107: Trigger UI regeneration when heatmap is removed
          heatmap_ui_trigger(heatmap_ui_trigger() + 1)

          # S2.9-FIX: Reset inhibit flag after UI rebuild (2500ms to cover full rebuild cycle)
          shinyjs::delay(2500, { inhibit_color_save(FALSE) })

          showNotification(paste("Heatmap", i, "removed"), type = "message")
        }
      }, ignoreInit = TRUE)
    })
  })

  # Generic observer for move up buttons
  observe({
    lapply(2:6, function(i) {
      observeEvent(input[[paste0("heatmap_up_", i)]], {
        if (i <= length(values$heatmap_configs) && i > 1) {
          configs <- values$heatmap_configs
          temp <- configs[[i]]
          configs[[i]] <- configs[[i-1]]
          configs[[i-1]] <- temp
          values$heatmap_configs <- configs
          # S2.9-FIX: Set inhibit flag before UI rebuild
          inhibit_color_save(TRUE)
          # v107: Trigger UI regeneration when heatmap is moved
          heatmap_ui_trigger(heatmap_ui_trigger() + 1)
          # S2.9-FIX: Reset inhibit flag after UI rebuild (2500ms to cover full rebuild cycle)
          shinyjs::delay(2500, { inhibit_color_save(FALSE) })
        }
      }, ignoreInit = TRUE)
    })
  })

  # Generic observer for move down buttons
  observe({
    lapply(1:5, function(i) {
      observeEvent(input[[paste0("heatmap_down_", i)]], {
        if (i < length(values$heatmap_configs)) {
          configs <- values$heatmap_configs
          temp <- configs[[i]]
          configs[[i]] <- configs[[i+1]]
          configs[[i+1]] <- temp
          values$heatmap_configs <- configs
          # S2.9-FIX: Set inhibit flag before UI rebuild
          inhibit_color_save(TRUE)
          # v107: Trigger UI regeneration when heatmap is moved
          heatmap_ui_trigger(heatmap_ui_trigger() + 1)
          # S2.9-FIX: Reset inhibit flag after UI rebuild (2500ms to cover full rebuild cycle)
          shinyjs::delay(2500, { inhibit_color_save(FALSE) })
        }
      }, ignoreInit = TRUE)
    })
  })
  
  # Observer to update heatmap config when inputs change
  observe({
    lapply(1:6, function(i) {
      # v56: Columns change (plural for multiple column support)
      observeEvent(input[[paste0("heatmap_columns_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$columns <- input[[paste0("heatmap_columns_", i)]]
        }
      }, ignoreInit = TRUE)
      
      # Title change
      observeEvent(input[[paste0("heatmap_title_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$title <- input[[paste0("heatmap_title_", i)]]
        }
      }, ignoreInit = TRUE)

      # S1.62dev: Data source change (csv or rdata)
      observeEvent(input[[paste0("heatmap_data_source_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$data_source <- input[[paste0("heatmap_data_source_", i)]]
          # If switching to RData, set type to continuous and appropriate colors
          if (input[[paste0("heatmap_data_source_", i)]] == "rdata") {
            cat(file=stderr(), paste0("\n[DEBUG-COLOR] Data source changed to 'rdata' for heatmap ", i, "\n"))
            values$heatmap_configs[[i]]$type <- "continuous"
            values$heatmap_configs[[i]]$auto_type <- FALSE
            # S1.62dev: Set red-white-blue color scheme for CNV (red=loss, blue=gain)
            values$heatmap_configs[[i]]$low_color <- "#FF0000"   # Red for deletion/loss
            values$heatmap_configs[[i]]$mid_color <- "#FFFFFF"   # White for neutral
            values$heatmap_configs[[i]]$high_color <- "#0000FF"  # Blue for amplification/gain
            values$heatmap_configs[[i]]$use_midpoint <- TRUE
            values$heatmap_configs[[i]]$midpoint <- 2  # Diploid baseline
            cat(file=stderr(), paste0("[DEBUG-COLOR] Set config colors: low=#FF0000, mid=#FFFFFF, high=#0000FF\n"))
            # S1.62dev: Update UI color pickers to reflect the new colors
            updateColourInput(session, paste0("heatmap_low_color_", i), value = "#FF0000")
            updateColourInput(session, paste0("heatmap_mid_color_", i), value = "#FFFFFF")
            updateColourInput(session, paste0("heatmap_high_color_", i), value = "#0000FF")
            updateCheckboxInput(session, paste0("heatmap_use_midpoint_", i), value = TRUE)
            updateNumericInput(session, paste0("heatmap_midpoint_", i), value = 2)
            cat(file=stderr(), "[DEBUG-COLOR] Called updateColourInput for low/mid/high colors\n")
          }
        }
      }, ignoreInit = TRUE)

      # S1.62dev: RData title change
      observeEvent(input[[paste0("heatmap_title_rdata_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$title <- input[[paste0("heatmap_title_rdata_", i)]]
        }
      }, ignoreInit = TRUE)

      # S1.62dev: CNV downsample factor change
      observeEvent(input[[paste0("heatmap_cnv_downsample_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$cnv_downsample <- input[[paste0("heatmap_cnv_downsample_", i)]]
        }
      }, ignoreInit = TRUE)

      # S1.62dev: CNV WGD normalization change
      # S2.13: Made mutually exclusive with per-cell WGD option
      observeEvent(input[[paste0("heatmap_cnv_wgd_norm_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_cnv_wgd_norm_", i)]]
          values$heatmap_configs[[i]]$cnv_wgd_norm <- new_val
          # If checked, uncheck per-cell WGD option (mutually exclusive)
          if (isTRUE(new_val)) {
            updateCheckboxInput(session, paste0("heatmap_cnv_wgd_per_cell_", i), value = FALSE)
            values$heatmap_configs[[i]]$cnv_wgd_per_cell <- FALSE
          }
        }
      }, ignoreInit = TRUE)

      # S2.12: Per-cell WGD checkbox change
      # S2.13: Made mutually exclusive with all-cells WGD option
      observeEvent(input[[paste0("heatmap_cnv_wgd_per_cell_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_cnv_wgd_per_cell_", i)]]
          values$heatmap_configs[[i]]$cnv_wgd_per_cell <- new_val
          # If checked, uncheck all-cells WGD option (mutually exclusive)
          if (isTRUE(new_val)) {
            updateCheckboxInput(session, paste0("heatmap_cnv_wgd_norm_", i), value = FALSE)
            values$heatmap_configs[[i]]$cnv_wgd_norm <- FALSE
          }
        }
      }, ignoreInit = TRUE)

      # S2.12: Per-cell WGD column selection change
      observeEvent(input[[paste0("heatmap_cnv_wgd_column_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$cnv_wgd_column <- input[[paste0("heatmap_cnv_wgd_column_", i)]]
        }
      }, ignoreInit = TRUE)

      # S2.8: Display mode change (basic vs detailed)
      observeEvent(input[[paste0("heatmap_cnv_display_mode_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$cnv_display_mode <- input[[paste0("heatmap_cnv_display_mode_", i)]]
          cat(file=stderr(), paste0("[HEATMAP-CONFIG] Heatmap ", i, " display mode changed to: ",
                                    input[[paste0("heatmap_cnv_display_mode_", i)]], "\n"))
        }
      }, ignoreInit = TRUE)

      # Render downsample factor change (second stage downsampling)
      observeEvent(input[[paste0("heatmap_cnv_render_downsample_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$cnv_render_downsample <- input[[paste0("heatmap_cnv_render_downsample_", i)]]
          cat(file=stderr(), paste0("[HEATMAP-CONFIG] Heatmap ", i, " render downsample changed to: ",
                                    input[[paste0("heatmap_cnv_render_downsample_", i)]], "\n"))
        }
      }, ignoreInit = TRUE)

      # Height scale change (for detailed mode)
      observeEvent(input[[paste0("heatmap_cnv_height_scale_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$cnv_height_scale <- input[[paste0("heatmap_cnv_height_scale_", i)]]
          cat(file=stderr(), paste0("[HEATMAP-CONFIG] Heatmap ", i, " height scale changed to: ",
                                    input[[paste0("heatmap_cnv_height_scale_", i)]], "\n"))
        }
      }, ignoreInit = TRUE)

      # Column names angle change
      observeEvent(input[[paste0("heatmap_colnames_angle_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$colnames_angle <- input[[paste0("heatmap_colnames_angle_", i)]]
        }
      }, ignoreInit = TRUE)

      # v107: Per-heatmap distance from tree (with guard to prevent reactive loop)
      # Changed ignoreInit to TRUE since values are now initialized in new_config
      observeEvent(input[[paste0("heatmap_distance_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_distance_", i)]]
          current_val <- values$heatmap_configs[[i]]$distance
          # Only update if value actually changed (prevents reactive loop)
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$distance <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # v107: Per-heatmap height (with guard to prevent reactive loop)
      # Changed ignoreInit to TRUE since values are now initialized in new_config
      observeEvent(input[[paste0("heatmap_height_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_height_", i)]]
          current_val <- values$heatmap_configs[[i]]$height
          # Only update if value actually changed (prevents reactive loop)
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$height <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # v110: Per-heatmap row height (controls visual height of each row)
      observeEvent(input[[paste0("heatmap_row_height_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_row_height_", i)]]
          current_val <- values$heatmap_configs[[i]]$row_height
          # Only update if value actually changed (prevents reactive loop)
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$row_height <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # v119: Grid settings observers (were missing)
      observeEvent(input[[paste0("heatmap_show_grid_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_show_grid_", i)]]
          current_val <- values$heatmap_configs[[i]]$show_grid
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$show_grid <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_grid_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_grid_color_", i)]]
          current_val <- values$heatmap_configs[[i]]$grid_color
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$grid_color <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_grid_size_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_grid_size_", i)]]
          current_val <- values$heatmap_configs[[i]]$grid_size
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$grid_size <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # S1.62dev: Row line settings observers
      observeEvent(input[[paste0("heatmap_show_row_lines_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_show_row_lines_", i)]]
          current_val <- values$heatmap_configs[[i]]$show_row_lines
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$show_row_lines <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_row_line_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_row_line_color_", i)]]
          current_val <- values$heatmap_configs[[i]]$row_line_color
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$row_line_color <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_row_line_size_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_row_line_size_", i)]]
          current_val <- values$heatmap_configs[[i]]$row_line_size
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$row_line_size <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # S1.62dev: Vertical column lines observers
      observeEvent(input[[paste0("heatmap_show_col_lines_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_show_col_lines_", i)]]
          current_val <- values$heatmap_configs[[i]]$show_col_lines
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$show_col_lines <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_col_line_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_col_line_color_", i)]]
          current_val <- values$heatmap_configs[[i]]$col_line_color
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$col_line_color <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_col_line_size_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_col_line_size_", i)]]
          current_val <- values$heatmap_configs[[i]]$col_line_size
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$col_line_size <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # S1.62dev: Vertical text labels observers
      observeEvent(input[[paste0("heatmap_show_vertical_text_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_show_vertical_text_", i)]]
          current_val <- values$heatmap_configs[[i]]$show_vertical_text
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$show_vertical_text <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_vertical_text_column_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_vertical_text_column_", i)]]
          current_val <- values$heatmap_configs[[i]]$vertical_text_column
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$vertical_text_column <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_vertical_text_size_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_vertical_text_size_", i)]]
          current_val <- values$heatmap_configs[[i]]$vertical_text_size
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$vertical_text_size <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_vertical_text_offset_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_vertical_text_offset_", i)]]
          current_val <- values$heatmap_configs[[i]]$vertical_text_offset
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$vertical_text_offset <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_vertical_text_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_vertical_text_color_", i)]]
          current_val <- values$heatmap_configs[[i]]$vertical_text_color
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$vertical_text_color <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # v119: Guide line settings observers (were missing - caused guide lines to not work)
      observeEvent(input[[paste0("heatmap_show_guides_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_show_guides_", i)]]
          current_val <- values$heatmap_configs[[i]]$show_guides
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$show_guides <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_guide_color1_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_guide_color1_", i)]]
          current_val <- values$heatmap_configs[[i]]$guide_color1
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$guide_color1 <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_guide_color2_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_guide_color2_", i)]]
          current_val <- values$heatmap_configs[[i]]$guide_color2
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$guide_color2 <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_guide_alpha_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_guide_alpha_", i)]]
          current_val <- values$heatmap_configs[[i]]$guide_alpha
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$guide_alpha <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_guide_width_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_guide_width_", i)]]
          current_val <- values$heatmap_configs[[i]]$guide_width
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$guide_width <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # S1.61: Guide line type observer
      observeEvent(input[[paste0("heatmap_guide_linetype_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_guide_linetype_", i)]]
          current_val <- values$heatmap_configs[[i]]$guide_linetype
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$guide_linetype <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # v107: Row labels settings (with guards to prevent reactive loop)
      # Changed ignoreInit to TRUE since values are now initialized in new_config
      observeEvent(input[[paste0("heatmap_show_row_labels_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_show_row_labels_", i)]]
          current_val <- values$heatmap_configs[[i]]$show_row_labels
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$show_row_labels <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_row_label_source_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_row_label_source_", i)]]
          current_val <- values$heatmap_configs[[i]]$row_label_source
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$row_label_source <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_row_label_font_size_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_row_label_font_size_", i)]]
          current_val <- values$heatmap_configs[[i]]$row_label_font_size
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$row_label_font_size <- new_val
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_custom_row_labels_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          new_val <- input[[paste0("heatmap_custom_row_labels_", i)]]
          current_val <- values$heatmap_configs[[i]]$custom_row_labels
          if (is.null(current_val) || !identical(new_val, current_val)) {
            values$heatmap_configs[[i]]$custom_row_labels <- new_val
          }
        }
      }, ignoreInit = TRUE)

      # Auto type change
      observeEvent(input[[paste0("heatmap_auto_type_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$auto_type <- input[[paste0("heatmap_auto_type_", i)]]
        }
      }, ignoreInit = TRUE)
      
      # Manual type change
      observeEvent(input[[paste0("heatmap_type_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$type <- input[[paste0("heatmap_type_", i)]]
        }
      }, ignoreInit = TRUE)
      
      # Discrete palette change
      # S1.62dev: Fixed - also update values$heatmaps and trigger plot regeneration
      observeEvent(input[[paste0("heatmap_discrete_palette_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$discrete_palette <- input[[paste0("heatmap_discrete_palette_", i)]]
          # S1.62dev: Also update values$heatmaps if it exists (heatmap was already applied)
          if (!is.null(values$heatmaps) && i <= length(values$heatmaps)) {
            values$heatmaps[[i]]$color_scheme <- input[[paste0("heatmap_discrete_palette_", i)]]
            request_plot_update()
          }
        }
      }, ignoreInit = TRUE)
      
      # Continuous palette change
      # S1.62dev: Fixed - also update values$heatmaps and trigger plot regeneration
      # S1.62dev: Also update Low/High/Mid color inputs based on selected palette
      observeEvent(input[[paste0("heatmap_cont_palette_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          palette_name <- input[[paste0("heatmap_cont_palette_", i)]]
          values$heatmap_configs[[i]]$cont_palette <- palette_name

          # S1.62dev: Don't override colors for RData heatmaps - they use red-white-blue by default
          current_data_source <- input[[paste0("heatmap_data_source_", i)]]
          if (!is.null(current_data_source) && current_data_source == "rdata") {
            cat(file=stderr(), paste0("[DEBUG-COLOR] Skipping palette color update for RData heatmap ", i, "\n"))
            return()  # Don't change colors for RData - user can manually change if needed
          }

          # S1.62dev: Define color mappings for each palette
          palette_colors <- list(
            # Sequential palettes (low=light, high=dark)
            Blues = list(low = "#DEEBF7", high = "#08519C", mid = "#6BAED6"),
            Greens = list(low = "#E5F5E0", high = "#006D2C", mid = "#74C476"),
            Reds = list(low = "#FEE0D2", high = "#A50F15", mid = "#FB6A4A"),
            Purples = list(low = "#EFEDF5", high = "#54278F", mid = "#9E9AC8"),
            Oranges = list(low = "#FEE6CE", high = "#A63603", mid = "#FD8D3C"),
            # Viridis family
            Viridis = list(low = "#FDE725", high = "#440154", mid = "#21918C"),
            Plasma = list(low = "#F0F921", high = "#0D0887", mid = "#CC4778"),
            Inferno = list(low = "#FCFFA4", high = "#000004", mid = "#BB3754"),
            Magma = list(low = "#FCFDBF", high = "#000004", mid = "#B63679"),
            # Diverging palettes (low=color1, mid=neutral, high=color2)
            RdBu = list(low = "#2166AC", high = "#B2182B", mid = "#F7F7F7"),
            RdYlGn = list(low = "#D73027", high = "#1A9850", mid = "#FFFFBF"),
            PiYG = list(low = "#C51B7D", high = "#4D9221", mid = "#F7F7F7"),
            BrBG = list(low = "#8C510A", high = "#01665E", mid = "#F5F5F5")
          )

          # Update color inputs if palette is recognized
          if (palette_name %in% names(palette_colors)) {
            colors <- palette_colors[[palette_name]]
            colourpicker::updateColourInput(session, paste0("heatmap_low_color_", i), value = colors$low)
            colourpicker::updateColourInput(session, paste0("heatmap_high_color_", i), value = colors$high)
            colourpicker::updateColourInput(session, paste0("heatmap_mid_color_", i), value = colors$mid)

            # Also update config
            values$heatmap_configs[[i]]$low_color <- colors$low
            values$heatmap_configs[[i]]$high_color <- colors$high
            values$heatmap_configs[[i]]$mid_color <- colors$mid
          }

          # S1.62dev: Also update values$heatmaps if it exists (heatmap was already applied)
          if (!is.null(values$heatmaps) && i <= length(values$heatmaps)) {
            values$heatmaps[[i]]$cont_palette <- palette_name
            if (palette_name %in% names(palette_colors)) {
              colors <- palette_colors[[palette_name]]
              values$heatmaps[[i]]$low_color <- colors$low
              values$heatmaps[[i]]$high_color <- colors$high
              values$heatmaps[[i]]$mid_color <- colors$mid
            }
            request_plot_update()
          }
        }
      }, ignoreInit = TRUE)
      
      # Color changes
      # S2.292dev: Added debug output to track when color observers fire
      # S2.292dev-FIX: Added inhibit_color_save check to prevent issues during UI rebuild
      observeEvent(input[[paste0("heatmap_low_color_", i)]], {
        # Check inhibit flag first
        if (isTRUE(isolate(inhibit_color_save()))) {
          cat(file=stderr(), sprintf("[COLOR-OBSERVER] heatmap_%d low_color BLOCKED by inhibit flag\n", i))
          return(NULL)
        }
        color_val <- input[[paste0("heatmap_low_color_", i)]]
        cat(file=stderr(), sprintf("[COLOR-OBSERVER] heatmap_%d low_color observer fired: '%s'\n", i, color_val))
        if (i <= length(values$heatmap_configs)) {
          old_val <- values$heatmap_configs[[i]]$low_color
          cat(file=stderr(), sprintf("[COLOR-OBSERVER]   Old value: '%s', New value: '%s'\n",
                                     ifelse(is.null(old_val), "NULL", old_val), color_val))
          values$heatmap_configs[[i]]$low_color <- color_val
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_high_color_", i)]], {
        # Check inhibit flag first
        if (isTRUE(isolate(inhibit_color_save()))) {
          cat(file=stderr(), sprintf("[COLOR-OBSERVER] heatmap_%d high_color BLOCKED by inhibit flag\n", i))
          return(NULL)
        }
        color_val <- input[[paste0("heatmap_high_color_", i)]]
        cat(file=stderr(), sprintf("[COLOR-OBSERVER] heatmap_%d high_color observer fired: '%s'\n", i, color_val))
        if (i <= length(values$heatmap_configs)) {
          old_val <- values$heatmap_configs[[i]]$high_color
          cat(file=stderr(), sprintf("[COLOR-OBSERVER]   Old value: '%s', New value: '%s'\n",
                                     ifelse(is.null(old_val), "NULL", old_val), color_val))
          values$heatmap_configs[[i]]$high_color <- color_val
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_mid_color_", i)]], {
        # Check inhibit flag first
        if (isTRUE(isolate(inhibit_color_save()))) {
          cat(file=stderr(), sprintf("[COLOR-OBSERVER] heatmap_%d mid_color BLOCKED by inhibit flag\n", i))
          return(NULL)
        }
        color_val <- input[[paste0("heatmap_mid_color_", i)]]
        cat(file=stderr(), sprintf("[COLOR-OBSERVER] heatmap_%d mid_color observer fired: '%s'\n", i, color_val))
        if (i <= length(values$heatmap_configs)) {
          old_val <- values$heatmap_configs[[i]]$mid_color
          cat(file=stderr(), sprintf("[COLOR-OBSERVER]   Old value: '%s', New value: '%s'\n",
                                     ifelse(is.null(old_val), "NULL", old_val), color_val))
          values$heatmap_configs[[i]]$mid_color <- color_val
        }
      }, ignoreInit = TRUE)

      # S2.292dev: Sync colourInput -> hex textInput for continuous colors
      observeEvent(input[[paste0("heatmap_low_color_", i)]], {
        color_val <- input[[paste0("heatmap_low_color_", i)]]
        if (!is.null(color_val)) {
          updateTextInput(session, paste0("heatmap_low_color_hex_", i), value = toupper(color_val))
        }
      }, ignoreInit = TRUE, priority = -1)

      observeEvent(input[[paste0("heatmap_high_color_", i)]], {
        color_val <- input[[paste0("heatmap_high_color_", i)]]
        if (!is.null(color_val)) {
          updateTextInput(session, paste0("heatmap_high_color_hex_", i), value = toupper(color_val))
        }
      }, ignoreInit = TRUE, priority = -1)

      observeEvent(input[[paste0("heatmap_mid_color_", i)]], {
        color_val <- input[[paste0("heatmap_mid_color_", i)]]
        if (!is.null(color_val)) {
          updateTextInput(session, paste0("heatmap_mid_color_hex_", i), value = toupper(color_val))
        }
      }, ignoreInit = TRUE, priority = -1)

      observeEvent(input[[paste0("heatmap_", i, "_cont_na_color")]], {
        color_val <- input[[paste0("heatmap_", i, "_cont_na_color")]]
        if (!is.null(color_val)) {
          updateTextInput(session, paste0("heatmap_", i, "_cont_na_color_hex"), value = toupper(color_val))
          # Also save to config - but check inhibit flag first
          if (!isTRUE(isolate(inhibit_color_save())) && i <= length(values$heatmap_configs)) {
            values$heatmap_configs[[i]]$na_color <- color_val
          }
        }
      }, ignoreInit = TRUE, priority = -1)

      # S2.292dev: Sync hex textInput -> colourInput for continuous colors
      observeEvent(input[[paste0("heatmap_low_color_hex_", i)]], {
        hex_val <- input[[paste0("heatmap_low_color_hex_", i)]]
        if (!is.null(hex_val) && nchar(trimws(hex_val)) > 0) {
          # Normalize hex value (add # if missing)
          hex_val <- trimws(hex_val)
          if (!grepl("^#", hex_val)) hex_val <- paste0("#", hex_val)
          # Validate hex format
          if (grepl("^#[0-9A-Fa-f]{6}$", hex_val)) {
            updateColourInput(session, paste0("heatmap_low_color_", i), value = hex_val)
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_high_color_hex_", i)]], {
        hex_val <- input[[paste0("heatmap_high_color_hex_", i)]]
        if (!is.null(hex_val) && nchar(trimws(hex_val)) > 0) {
          hex_val <- trimws(hex_val)
          if (!grepl("^#", hex_val)) hex_val <- paste0("#", hex_val)
          if (grepl("^#[0-9A-Fa-f]{6}$", hex_val)) {
            updateColourInput(session, paste0("heatmap_high_color_", i), value = hex_val)
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_mid_color_hex_", i)]], {
        hex_val <- input[[paste0("heatmap_mid_color_hex_", i)]]
        if (!is.null(hex_val) && nchar(trimws(hex_val)) > 0) {
          hex_val <- trimws(hex_val)
          if (!grepl("^#", hex_val)) hex_val <- paste0("#", hex_val)
          if (grepl("^#[0-9A-Fa-f]{6}$", hex_val)) {
            updateColourInput(session, paste0("heatmap_mid_color_", i), value = hex_val)
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_", i, "_cont_na_color_hex")]], {
        hex_val <- input[[paste0("heatmap_", i, "_cont_na_color_hex")]]
        if (!is.null(hex_val) && nchar(trimws(hex_val)) > 0) {
          hex_val <- trimws(hex_val)
          if (!grepl("^#", hex_val)) hex_val <- paste0("#", hex_val)
          if (grepl("^#[0-9A-Fa-f]{6}$", hex_val)) {
            updateColourInput(session, paste0("heatmap_", i, "_cont_na_color"), value = hex_val)
          }
        }
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("heatmap_use_midpoint_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$use_midpoint <- input[[paste0("heatmap_use_midpoint_", i)]]
        }
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("heatmap_midpoint_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$midpoint <- input[[paste0("heatmap_midpoint_", i)]]
        }
      }, ignoreInit = TRUE)

      # v121: Column range observer - allows adding columns by position range (e.g., "2-10")
      observeEvent(input[[paste0("heatmap_add_range_", i)]], {
        range_text <- input[[paste0("heatmap_col_range_", i)]]
        if (is.null(range_text) || nchar(trimws(range_text)) == 0) {
          showNotification("Please enter a column range (e.g., 2-10)", type = "warning")
          return()
        }

        # Parse the range
        range_text <- trimws(range_text)
        parts <- strsplit(range_text, "-")[[1]]
        if (length(parts) != 2) {
          showNotification("Invalid range format. Use 'start-end' (e.g., 2-10)", type = "error")
          return()
        }

        start_col <- suppressWarnings(as.integer(trimws(parts[1])))
        end_col <- suppressWarnings(as.integer(trimws(parts[2])))

        if (is.na(start_col) || is.na(end_col)) {
          showNotification("Column numbers must be integers (e.g., 2-10)", type = "error")
          return()
        }

        if (start_col < 1 || end_col < start_col) {
          showNotification("Invalid range: start must be >= 1 and end must be >= start", type = "error")
          return()
        }

        # Get current CSV column names
        csv_cols <- if (!is.null(values$csv_data)) names(values$csv_data) else character(0)
        if (length(csv_cols) == 0) {
          showNotification("No CSV data loaded", type = "error")
          return()
        }

        if (end_col > length(csv_cols)) {
          showNotification(paste0("End column (", end_col, ") exceeds number of columns (", length(csv_cols), ")"), type = "error")
          return()
        }

        # Get column names for the range
        range_cols <- csv_cols[start_col:end_col]

        # Get current selection and add range columns
        current_cols <- input[[paste0("heatmap_columns_", i)]]
        if (is.null(current_cols)) current_cols <- character(0)

        # Add new columns (avoiding duplicates)
        new_cols <- unique(c(current_cols, range_cols))

        # Update the selectize input
        updateSelectizeInput(session, paste0("heatmap_columns_", i), selected = new_cols)

        # Also update config
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$columns <- new_cols
        }

        # Clear the range input
        updateTextInput(session, paste0("heatmap_col_range_", i), value = "")

        showNotification(paste0("Added columns ", start_col, "-", end_col, " (", length(range_cols), " columns): ",
                                paste(range_cols, collapse = ", ")), type = "message")
      }, ignoreInit = TRUE)
    })
  })

  # v62: Render palette previews for discrete heatmaps
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_discrete_palette_preview_", i)]] <- renderUI({
        palette_name <- input[[paste0("heatmap_discrete_palette_", i)]]
        if (is.null(palette_name)) palette_name <- "Set1"

        # Get colors from the selected palette (8 colors for preview)
        n_colors <- 8
        colors <- tryCatch({
          if (palette_name %in% c("Set1", "Set2", "Set3", "Paired", "Dark2", "Accent", "Pastel1", "Pastel2")) {
            RColorBrewer::brewer.pal(min(n_colors, RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]), palette_name)
          } else {
            rainbow(n_colors)
          }
        }, error = function(e) rainbow(n_colors))

        # Create color swatches
        swatches <- lapply(colors, function(col) {
          tags$span(style = paste0(
            "display: inline-block; width: 20px; height: 20px; ",
            "background-color: ", col, "; margin-right: 2px; border: 1px solid #ccc; border-radius: 2px;"
          ))
        })

        tags$div(
          style = "margin-top: 5px;",
          do.call(tagList, swatches)
        )
      })
    })
  })

  # v62: Render palette previews for continuous heatmaps
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_cont_palette_preview_", i)]] <- renderUI({
        palette_name <- input[[paste0("heatmap_cont_palette_", i)]]
        if (is.null(palette_name)) palette_name <- "Blues"

        # Get colors from the selected palette (gradient preview with 10 colors)
        n_colors <- 10
        colors <- tryCatch({
          if (palette_name %in% c("Viridis", "Plasma", "Inferno", "Magma")) {
            viridis::viridis(n_colors, option = tolower(palette_name))
          } else if (palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
            colorRampPalette(RColorBrewer::brewer.pal(9, palette_name))(n_colors)
          } else {
            colorRampPalette(c("white", "blue"))(n_colors)
          }
        }, error = function(e) colorRampPalette(c("white", "blue"))(n_colors))

        # Create gradient preview as a bar
        gradient_css <- paste0(
          "background: linear-gradient(to right, ",
          paste(colors, collapse = ", "), ");"
        )

        tags$div(
          style = paste0(
            "margin-top: 5px; height: 20px; border-radius: 3px; border: 1px solid #ccc; ",
            gradient_css
          )
        )
      })
    })
  })

  # v64: Render palette previews for discrete heatmaps
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_discrete_palette_preview_", i)]] <- renderUI({
        palette_name <- input[[paste0("heatmap_discrete_palette_", i)]]
        if (is.null(palette_name)) palette_name <- "Set1"

        # Get colors from the selected palette
        n_colors <- 8  # Show 8 colors for discrete palette preview
        colors <- tryCatch({
          if (palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
            max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
            RColorBrewer::brewer.pal(min(n_colors, max_colors), palette_name)
          } else {
            rainbow(n_colors)
          }
        }, error = function(e) rainbow(n_colors))

        # Create discrete color swatches preview
        color_boxes <- lapply(colors, function(col) {
          tags$div(
            style = paste0(
              "display: inline-block; width: ", 100/length(colors), "%; ",
              "height: 20px; background-color: ", col, ";"
            )
          )
        })

        tags$div(
          style = "margin-top: 5px; border-radius: 3px; border: 1px solid #ccc; overflow: hidden;",
          do.call(tags$div, c(list(style = "display: flex;"), color_boxes))
        )
      })
    })
  })

  # v70: R color names list for dropdown menus (used in both classification and heatmap)
  heat_r_colors <- c(
    "Custom" = "",
    # Basic colors
    "red", "blue", "green", "yellow", "orange", "purple", "pink", "brown",
    "gray", "black", "white", "cyan", "magenta",
    # Dark variants
    "darkred", "darkblue", "darkgreen", "darkorange", "darkviolet",
    "darkgray", "darkcyan", "darkmagenta",
    # Light variants
    "lightblue", "lightgreen", "lightyellow", "lightpink", "lightgray",
    "lightcyan", "lightcoral", "lightsalmon",
    # Named colors
    "steelblue", "skyblue", "navy", "maroon", "olive", "teal", "coral",
    "tomato", "salmon", "khaki", "plum", "orchid", "tan",
    # Greens
    "forestgreen", "limegreen", "seagreen", "springgreen",
    # Blues
    "royalblue", "dodgerblue", "deepskyblue", "cornflowerblue",
    # Reds/Pinks
    "crimson", "firebrick", "indianred", "hotpink", "deeppink"
  )

  # v127: Dynamic detected type display - updates when column selection changes
  # This replaces the static label that was set at UI render time
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_detected_type_display_", i)]] <- renderUI({
        # Take dependency on column selection to reactively update
        cols_selected <- input[[paste0("heatmap_columns_", i)]]

        # If no columns selected, show placeholder
        if (is.null(cols_selected) || length(cols_selected) == 0 || is.null(values$csv_data)) {
          return(tags$div(
            class = "label label-default",
            style = "display: block; padding: 8px; text-align: center; margin-top: 0;",
            "Select column first"
          ))
        }

        # Compute detected type from the first column (same logic as heatmap_type_settings_ui_)
        first_col <- cols_selected[1]
        if (!(first_col %in% names(values$csv_data))) {
          return(tags$div(
            class = "label label-warning",
            style = "display: block; padding: 8px; text-align: center; margin-top: 0;",
            "Column not found"
          ))
        }

        col_data <- values$csv_data[[first_col]]
        unique_vals <- length(unique(na.omit(col_data)))
        is_numeric <- is.numeric(col_data)

        # v127: Check for decimals in string representation BEFORE conversion
        has_decimal_in_string <- FALSE
        if (is.character(col_data) || is.factor(col_data)) {
          char_data <- as.character(na.omit(col_data))
          char_data <- char_data[!toupper(trimws(char_data)) %in% c("NA", "N/A", "#N/A", "NULL", "")]
          if (length(char_data) > 0) {
            has_decimal_in_string <- any(grepl("\\.[0-9]+", char_data))
          }
        } else if (is.numeric(col_data)) {
          char_data <- as.character(na.omit(col_data))
          has_decimal_in_string <- any(grepl("\\.[0-9]+", char_data))
        }

        # Try to convert character columns to numeric
        originally_numeric <- is_numeric
        if (!is_numeric && is.character(col_data)) {
          clean_col_data <- col_data
          clean_col_data[toupper(trimws(clean_col_data)) %in% c("NA", "N/A", "#N/A", "NULL", "")] <- NA
          numeric_attempt <- suppressWarnings(as.numeric(clean_col_data))
          non_na_original <- sum(!is.na(clean_col_data))
          non_na_converted <- sum(!is.na(numeric_attempt))
          if (non_na_original > 0 && (non_na_converted / non_na_original) >= 0.5) {
            is_numeric <- TRUE
            col_data <- numeric_attempt
          }
        }

        # Determine type
        detected_type <- "discrete"
        if (is_numeric) {
          non_na_vals <- na.omit(col_data)
          unique_numeric_vals <- length(unique(non_na_vals))
          epsilon <- 1e-6
          has_decimals <- any(abs(non_na_vals - floor(non_na_vals)) > epsilon, na.rm = TRUE)

          if (!has_decimals && has_decimal_in_string) {
            has_decimals <- TRUE
          }

          val_range <- if (length(non_na_vals) > 0) diff(range(non_na_vals)) else 0

          if (has_decimals) {
            detected_type <- "continuous"
          } else if (originally_numeric) {
            is_boolean_like <- unique_numeric_vals <= 2 && val_range <= 1
            is_small_categorical <- unique_numeric_vals <= 3 && val_range <= 2
            detected_type <- if (is_boolean_like || is_small_categorical) "discrete" else "continuous"
          } else {
            detected_type <- if (unique_numeric_vals > 8 || val_range > 10) "continuous" else "discrete"
          }
        }

        # Debug output to console
        debug_cat(paste0("  v27 DETECTED TYPE DISPLAY: column=", first_col,
                                   ", detected=", detected_type, "\n"))

        # Return the styled label
        tags$div(
          class = if (detected_type == "continuous") "label label-info" else "label label-success",
          style = "display: block; padding: 8px; text-align: center; margin-top: 0;",
          if (detected_type == "continuous") "Continuous (numeric)" else "Discrete (categorical)"
        )
      })
    })
  })

  # v108: New reactive renderUI for type-specific settings (discrete vs continuous)
  # This replaces the old conditionalPanels which used a hardcoded detected_type
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_type_settings_ui_", i)]] <- renderUI({
        # Take dependency on column selection to reactively update when columns change
        cols_selected <- input[[paste0("heatmap_columns_", i)]]
        auto_type <- input[[paste0("heatmap_auto_type_", i)]]
        forced_type <- input[[paste0("heatmap_type_", i)]]

        # S1.62dev: Check if this is an RData heatmap - show continuous settings directly
        data_source <- input[[paste0("heatmap_data_source_", i)]]
        if (!is.null(data_source) && data_source == "rdata") {
          # RData CNV heatmaps are always continuous - show continuous color settings
          cfg <- isolate(values$heatmap_configs[[i]])
          # S2.292dev: Debug output for renderUI
          cat(file=stderr(), sprintf("[RENDER-UI] RData heatmap %d color settings from cfg:\n", i))
          cat(file=stderr(), sprintf("[RENDER-UI]   cfg is NULL: %s\n", is.null(cfg)))
          if (!is.null(cfg)) {
            cat(file=stderr(), sprintf("[RENDER-UI]   cfg$low_color: %s\n", ifelse(is.null(cfg$low_color), "NULL", cfg$low_color)))
            cat(file=stderr(), sprintf("[RENDER-UI]   cfg$mid_color: %s\n", ifelse(is.null(cfg$mid_color), "NULL", cfg$mid_color)))
            cat(file=stderr(), sprintf("[RENDER-UI]   cfg$high_color: %s\n", ifelse(is.null(cfg$high_color), "NULL", cfg$high_color)))
          }
          continuous_palettes <- c("Blues", "Greens", "Reds", "Purples", "Oranges",
                                   "Viridis", "Plasma", "Inferno", "Magma",
                                   "RdBu", "RdYlGn", "PiYG", "BrBG")
          # S2.292dev: Add manual hex input next to each color picker
          return(tags$div(
            style = "background-color: #e6f3ff; padding: 10px; border-radius: 5px; margin-top: 10px;",
            tags$h5(icon("dna"), " CNV Color Settings"),
            fluidRow(
              column(4,
                     selectInput(paste0("heatmap_cont_palette_", i), "Color Palette",
                                 choices = continuous_palettes,
                                 selected = if (!is.null(cfg$cont_palette)) cfg$cont_palette else "RdBu")
              ),
              column(4,
                     # S1.62dev: Default red for deletion/loss
                     tags$label("Low (Deletion)"),
                     tags$div(style = "display: flex; gap: 5px; align-items: center;",
                       colourInput(paste0("heatmap_low_color_", i), NULL,
                                   value = if (!is.null(cfg$low_color)) cfg$low_color else "#FF0000"),
                       textInput(paste0("heatmap_low_color_hex_", i), NULL,
                                 value = if (!is.null(cfg$low_color)) cfg$low_color else "#FF0000",
                                 width = "90px")
                     )
              ),
              column(4,
                     # S1.62dev: Default blue for amplification/gain
                     tags$label("High (Amplification)"),
                     tags$div(style = "display: flex; gap: 5px; align-items: center;",
                       colourInput(paste0("heatmap_high_color_", i), NULL,
                                   value = if (!is.null(cfg$high_color)) cfg$high_color else "#0000FF"),
                       textInput(paste0("heatmap_high_color_hex_", i), NULL,
                                 value = if (!is.null(cfg$high_color)) cfg$high_color else "#0000FF",
                                 width = "90px")
                     )
              )
            ),
            fluidRow(
              column(4,
                     checkboxInput(paste0("heatmap_use_midpoint_", i), "Use midpoint color",
                                   value = if (!is.null(cfg$use_midpoint)) cfg$use_midpoint else TRUE)
              ),
              column(4,
                     conditionalPanel(
                       condition = paste0("input.heatmap_use_midpoint_", i),
                       tags$label("Mid (Diploid)"),
                       tags$div(style = "display: flex; gap: 5px; align-items: center;",
                         colourInput(paste0("heatmap_mid_color_", i), NULL,
                                     value = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFFFF"),
                         textInput(paste0("heatmap_mid_color_hex_", i), NULL,
                                   value = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFFFF",
                                   width = "90px")
                       )
                     )
              ),
              column(4,
                     conditionalPanel(
                       condition = paste0("input.heatmap_use_midpoint_", i),
                       numericInput(paste0("heatmap_midpoint_", i), "Midpoint Value",
                                    value = if (!is.null(cfg$midpoint)) cfg$midpoint else 0,
                                    step = 0.5)
                     )
              )
            ),
            fluidRow(
              column(4,
                     tags$label("NA Color"),
                     tags$div(style = "display: flex; gap: 5px; align-items: center;",
                       colourInput(paste0("heatmap_", i, "_cont_na_color"), NULL,
                                   value = if (!is.null(cfg$na_color)) cfg$na_color else "#BEBEBE",
                                   showColour = "background"),
                       textInput(paste0("heatmap_", i, "_cont_na_color_hex"), NULL,
                                 value = if (!is.null(cfg$na_color)) cfg$na_color else "#BEBEBE",
                                 width = "90px")
                     )
              ),
              column(8,
                     tags$p(class = "text-muted", style = "margin-top: 25px;",
                            icon("info-circle"), " Blue-White-Red scale centered at 0")
              )
            )
          ))
        }

        # If no columns selected, show nothing
        if (is.null(cols_selected) || length(cols_selected) == 0 || is.null(values$csv_data)) {
          return(tags$p(class = "text-muted", style = "margin-top: 10px;",
                        icon("info-circle"), " Select column(s) to see color settings"))
        }

        # Compute detected type from the first column
        first_col <- cols_selected[1]
        if (!(first_col %in% names(values$csv_data))) {
          return(tags$p(class = "text-muted", "Column not found"))
        }

        col_data <- values$csv_data[[first_col]]
        unique_vals <- length(unique(na.omit(col_data)))
        is_numeric <- is.numeric(col_data)

        # v117: Check for decimals in string representation BEFORE conversion
        # This catches values like "15.7" that might lose precision
        # Also check originally numeric data by converting to string
        # IMPORTANT: Filter out "NA" strings which are not R's NA
        has_decimal_in_string <- FALSE
        if (is.character(col_data) || is.factor(col_data)) {
          char_data <- as.character(na.omit(col_data))
          # v117: Remove "NA" strings (case-insensitive) that are NOT R's NA
          # v118: Added #N/A to the list of NA-like strings
          char_data <- char_data[!toupper(trimws(char_data)) %in% c("NA", "N/A", "#N/A", "NULL", "")]
          if (length(char_data) > 0) {
            has_decimal_in_string <- any(grepl("\\.[0-9]+", char_data))
          }
        } else if (is.numeric(col_data)) {
          # v116: For originally numeric data, convert to string and check for decimals
          char_data <- as.character(na.omit(col_data))
          has_decimal_in_string <- any(grepl("\\.[0-9]+", char_data))
        }

        # v111: Better detection - also try to convert character columns to numeric
        originally_numeric <- is_numeric  # v112: Track if originally numeric
        if (!is_numeric && is.character(col_data)) {
          # v117: Filter out NA-like strings before attempting conversion
          # v118: Added #N/A to the list of NA-like strings
          clean_col_data <- col_data
          clean_col_data[toupper(trimws(clean_col_data)) %in% c("NA", "N/A", "#N/A", "NULL", "")] <- NA
          # Try converting to numeric - see what proportion succeeds
          numeric_attempt <- suppressWarnings(as.numeric(clean_col_data))
          non_na_original <- sum(!is.na(clean_col_data))  # v117: Count after cleaning NA strings
          non_na_converted <- sum(!is.na(numeric_attempt))
          # v116: Lower threshold to 50% for consistency with heatmap config
          if (non_na_original > 0 && (non_na_converted / non_na_original) >= 0.5) {
            is_numeric <- TRUE
            col_data <- numeric_attempt  # Use converted data for further analysis
          }
        }

        # v112: Better detection for numeric columns
        if (is_numeric) {
          # Recalculate unique values after numeric conversion (important!)
          non_na_vals <- na.omit(col_data)
          unique_numeric_vals <- length(unique(non_na_vals))

          # v116: More robust decimal detection with epsilon for floating-point precision
          epsilon <- 1e-6
          has_decimals <- any(abs(non_na_vals - floor(non_na_vals)) > epsilon, na.rm = TRUE)

          # v116: Also use string-based detection as fallback
          if (!has_decimals && has_decimal_in_string) {
            has_decimals <- TRUE
          }

          # Check the range of values
          val_range <- if (length(non_na_vals) > 0) diff(range(non_na_vals)) else 0

          # v116: Simplified logic - decimals ALWAYS mean continuous
          if (has_decimals) {
            detected_type <- "continuous"
          } else if (originally_numeric) {
            # Originally numeric without decimals
            is_boolean_like <- unique_numeric_vals <= 2 && val_range <= 1
            is_small_categorical <- unique_numeric_vals <= 3 && val_range <= 2
            detected_type <- if (is_boolean_like || is_small_categorical) "discrete" else "continuous"
          } else {
            # Converted from character without decimals - be more conservative
            detected_type <- if (unique_numeric_vals > 8 || val_range > 10) "continuous" else "discrete"
          }
        } else {
          detected_type <- "discrete"
        }

        # Determine actual type based on auto-detect checkbox
        actual_type <- if (isTRUE(auto_type) || is.null(auto_type)) detected_type else forced_type
        if (is.null(actual_type)) actual_type <- "discrete"

        # Get config for current values (isolated to prevent loops)
        cfg <- isolate(values$heatmap_configs[[i]])

        # Palette options
        discrete_palettes <- c("Set1", "Set2", "Set3", "Paired", "Dark2", "Accent", "Pastel1", "Pastel2")
        continuous_palettes <- c("Blues", "Greens", "Reds", "Purples", "Oranges",
                                 "Viridis", "Plasma", "Inferno", "Magma",
                                 "RdBu", "RdYlGn", "PiYG", "BrBG")

        if (actual_type == "discrete") {
          # Discrete settings UI
          tags$div(
            style = "background-color: #f9f9f9; padding: 10px; border-radius: 5px; margin-top: 10px;",
            tags$h5(icon("palette"), " Discrete Color Settings"),
            tags$div(
              style = "background-color: #f0f0f0; padding: 10px; margin-bottom: 10px; border-radius: 5px;",
              fluidRow(
                column(6,
                       selectInput(paste0("heatmap_discrete_palette_", i), "Color Palette",
                                   choices = discrete_palettes,
                                   selected = if (!is.null(cfg$discrete_palette)) cfg$discrete_palette else "Set1")
                ),
                column(6, style = "padding-top: 25px;",
                       actionButton(paste0("apply_heatmap_palette_", i), "Apply Palette to All",
                                    icon = icon("palette"), class = "btn-info btn-sm")
                )
              ),
              uiOutput(paste0("heatmap_discrete_palette_preview_", i))
            ),
            tags$h6("Value Colors:", style = "margin-top: 10px; margin-bottom: 5px;"),
            uiOutput(paste0("heatmap_discrete_colors_ui_", i))
          )
        } else {
          # Continuous settings UI
          # S2.292dev: Add manual hex input next to each color picker
          tags$div(
            style = "background-color: #f0f7ff; padding: 10px; border-radius: 5px; margin-top: 10px;",
            tags$h5(icon("sliders-h"), " Continuous Color Settings"),
            fluidRow(
              column(4,
                     selectInput(paste0("heatmap_cont_palette_", i), "Color Palette",
                                 choices = continuous_palettes,
                                 selected = if (!is.null(cfg$cont_palette)) cfg$cont_palette else "Blues"),
                     uiOutput(paste0("heatmap_cont_palette_preview_", i))
              ),
              column(4,
                     tags$label("Low Color"),
                     tags$div(style = "display: flex; gap: 5px; align-items: center;",
                       colourInput(paste0("heatmap_low_color_", i), NULL,
                                   value = if (!is.null(cfg$low_color)) cfg$low_color else "#FFFFCC"),
                       textInput(paste0("heatmap_low_color_hex_", i), NULL,
                                 value = if (!is.null(cfg$low_color)) cfg$low_color else "#FFFFCC",
                                 width = "90px")
                     )
              ),
              column(4,
                     tags$label("High Color"),
                     tags$div(style = "display: flex; gap: 5px; align-items: center;",
                       colourInput(paste0("heatmap_high_color_", i), NULL,
                                   value = if (!is.null(cfg$high_color)) cfg$high_color else "#006837"),
                       textInput(paste0("heatmap_high_color_hex_", i), NULL,
                                 value = if (!is.null(cfg$high_color)) cfg$high_color else "#006837",
                                 width = "90px")
                     )
              )
            ),
            fluidRow(
              column(4,
                     checkboxInput(paste0("heatmap_use_midpoint_", i), "Use midpoint color",
                                   value = if (!is.null(cfg$use_midpoint)) cfg$use_midpoint else FALSE)
              ),
              column(4,
                     conditionalPanel(
                       condition = paste0("input.heatmap_use_midpoint_", i),
                       tags$label("Mid Color"),
                       tags$div(style = "display: flex; gap: 5px; align-items: center;",
                         colourInput(paste0("heatmap_mid_color_", i), NULL,
                                     value = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFF99"),
                         textInput(paste0("heatmap_mid_color_hex_", i), NULL,
                                   value = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFF99",
                                   width = "90px")
                       )
                     )
              ),
              column(4,
                     conditionalPanel(
                       condition = paste0("input.heatmap_use_midpoint_", i),
                       numericInput(paste0("heatmap_midpoint_", i), "Midpoint Value",
                                    value = if (!is.null(cfg$midpoint)) cfg$midpoint else 0,
                                    step = 0.1)
                     )
              )
            ),
            # v112: NA color for continuous heatmaps
            fluidRow(
              column(4,
                     tags$label("NA Color"),
                     tags$div(style = "display: flex; gap: 5px; align-items: center;",
                       colourInput(paste0("heatmap_", i, "_cont_na_color"), NULL,
                                   value = if (!is.null(cfg$na_color)) cfg$na_color else "#BEBEBE",
                                   showColour = "background"),
                       textInput(paste0("heatmap_", i, "_cont_na_color_hex"), NULL,
                                 value = if (!is.null(cfg$na_color)) cfg$na_color else "#BEBEBE",
                                 width = "90px")
                     )
              ),
              column(8)
            )
          )
        }
      })
    })
  })

  # v108: Render custom labels UI (mapping or comma-separated)
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_custom_labels_ui_", i)]] <- renderUI({
        # Depend on label source selection and column selection
        label_source <- input[[paste0("heatmap_row_label_source_", i)]]
        cols_selected <- input[[paste0("heatmap_columns_", i)]]

        # Show nothing for "colnames" source
        if (is.null(label_source) || label_source == "colnames") {
          return(NULL)
        }

        # Show comma-separated text input for "custom" source
        if (label_source == "custom") {
          cfg <- isolate(values$heatmap_configs[[i]])
          return(fluidRow(
            column(12,
                   textInput(paste0("heatmap_custom_row_labels_", i), "Custom labels (comma-separated)",
                             value = if (!is.null(cfg$custom_row_labels)) cfg$custom_row_labels else "",
                             placeholder = "Label1, Label2, Label3...")
            )
          ))
        }

        # Show per-column mapping UI for "mapping" source
        if (label_source == "mapping") {
          if (is.null(cols_selected) || length(cols_selected) == 0) {
            return(tags$p(class = "text-muted",
                          icon("info-circle"), " Select columns first to map labels"))
          }

          cfg <- isolate(values$heatmap_configs[[i]])
          existing_mapping <- if (!is.null(cfg$label_mapping)) cfg$label_mapping else list()

          # Create a text input for each column
          mapping_rows <- lapply(seq_along(cols_selected), function(j) {
            col_name <- cols_selected[j]
            # Get existing mapping or default to column name
            existing_label <- if (!is.null(existing_mapping[[col_name]])) existing_mapping[[col_name]] else col_name

            fluidRow(
              style = "margin-bottom: 5px; padding: 3px 0;",
              column(5,
                     tags$div(
                       style = "padding-top: 7px; font-family: monospace; font-size: 11px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;",
                       title = col_name,
                       col_name
                     )
              ),
              column(1,
                     tags$div(style = "padding-top: 7px; text-align: center;", icon("arrow-right"))
              ),
              column(6,
                     textInput(paste0("heatmap_", i, "_label_", j), NULL,
                               value = existing_label,
                               placeholder = "Custom label...")
              )
            )
          })

          return(tags$div(
            style = "background-color: #fff; padding: 10px; border-radius: 5px; border: 1px solid #e0e0e0;",
            tags$h6("Column to Label Mapping:", style = "margin-top: 0; margin-bottom: 8px;"),
            tags$small(class = "text-muted", paste0(length(cols_selected), " column(s) - edit labels on the right")),
            tags$hr(style = "margin: 8px 0;"),
            do.call(tagList, mapping_rows)
          ))
        }

        return(NULL)
      })
    })
  })

  # S2.0: Render RData sample mapping column selector
  # Shows when auto-match fails and user needs to select which CSV column has the sample names
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_rdata_mapping_ui_", i)]] <- renderUI({
        # Only show for RData data source
        data_source <- input[[paste0("heatmap_data_source_", i)]]
        if (is.null(data_source) || data_source != "rdata") {
          return(NULL)
        }

        # Check if RData is loaded
        if (is.null(values$rdata_cnv_matrix)) {
          return(tags$div(
            style = "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            icon("exclamation-triangle"),
            " Please load an RData CNV file first (Data Import tab)"
          ))
        }

        # Check auto-match status
        if (isTRUE(values$rdata_auto_match)) {
          # Auto-match succeeded - show success message
          return(tags$div(
            style = "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            icon("check-circle"),
            " RData samples auto-matched to tree tips - ready to use!"
          ))
        }

        # Auto-match failed - show dropdown for manual mapping
        # Get CSV column names
        csv_cols <- if (!is.null(values$csv_data)) names(values$csv_data) else character(0)

        # S2.8: Get the saved mapping column from heatmap_configs if available
        saved_mapping_col <- NULL
        if (!is.null(values$heatmap_configs) && length(values$heatmap_configs) >= i) {
          saved_mapping_col <- values$heatmap_configs[[i]]$rdata_mapping_column
        }
        # Fall back to global value if no per-heatmap value
        if (is.null(saved_mapping_col) || saved_mapping_col == "") {
          saved_mapping_col <- values$rdata_mapping_column
        }

        # Show RData sample names preview and mapping dropdown
        sample_preview <- paste(head(values$rdata_sample_names, 3), collapse=", ")
        if (length(values$rdata_sample_names) > 3) {
          sample_preview <- paste0(sample_preview, ", ...")
        }

        tags$div(
          style = "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
          tags$div(
            style = "margin-bottom: 8px;",
            icon("exclamation-triangle"),
            tags$strong(" Manual mapping needed"),
            tags$br(),
            tags$small(
              style = "color: #856404;",
              "RData sample names (", sample_preview, ") don't match tree tip labels."
            )
          ),
          fluidRow(
            column(8,
                   selectInput(paste0("heatmap_rdata_mapping_col_", i),
                               "Select CSV column with sample names:",
                               choices = c("-- Select a column --" = "", csv_cols),
                               selected = if (!is.null(saved_mapping_col) && saved_mapping_col != "") saved_mapping_col else "")
            ),
            column(4,
                   tags$div(
                     style = "padding-top: 25px;",
                     tags$small(
                       class = "text-muted",
                       "Choose the column that maps tree IDs to RData sample names"
                     )
                   )
            )
          )
        )
      })
    })
  })

  # S2.0: Observer to update rdata_mapping_column when user selects a column
  observe({
    lapply(1:6, function(i) {
      observeEvent(input[[paste0("heatmap_rdata_mapping_col_", i)]], {
        selected_col <- input[[paste0("heatmap_rdata_mapping_col_", i)]]
        if (!is.null(selected_col) && selected_col != "") {
          values$rdata_mapping_column <- selected_col
          # Also update the per-heatmap config so it gets saved to YAML
          if (!is.null(values$heatmap_configs) && length(values$heatmap_configs) >= i) {
            values$heatmap_configs[[i]]$rdata_mapping_column <- selected_col
          }
          cat(file=stderr(), paste0("[RDATA-MAPPING] User selected mapping column for heatmap ", i, ": '", selected_col, "'\n"))
        }
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    })
  })

  # v70: Render discrete color pickers for each heatmap - with NA color and dropdown menus
  observe({
    lapply(1:6, function(i) {
      output[[paste0("heatmap_discrete_colors_ui_", i)]] <- renderUI({
        # v69: Use columns (plural) - get unique values from first column
        cols_selected <- input[[paste0("heatmap_columns_", i)]]
        if (is.null(cols_selected) || length(cols_selected) == 0 || is.null(values$csv_data)) {
          return(tags$p(class = "text-muted", "Select column(s) first to see value-color mappings"))
        }

        # Get unique values from first column (for discrete coloring)
        first_col <- cols_selected[1]
        if (!(first_col %in% names(values$csv_data))) {
          return(tags$p(class = "text-muted", "Column not found"))
        }

        # v125: Filter to only patient-specific data (rows matching tree tips)
        # This prevents "Too many unique values" error when CSV has more values than tree tips
        filtered_data <- values$csv_data
        if (!is.null(values$tree) && !is.null(input$id_column) && input$id_column %in% names(values$csv_data)) {
          tree_tips <- values$tree$tip.label
          id_col_data <- as.character(values$csv_data[[input$id_column]])
          matching_rows <- id_col_data %in% tree_tips
          if (any(matching_rows)) {
            filtered_data <- values$csv_data[matching_rows, , drop = FALSE]
            debug_cat(paste0("v125: Filtered unique values from ", nrow(values$csv_data), " to ", nrow(filtered_data), " rows (tree tips)\n"))
          }
        }

        unique_vals <- sort(unique(na.omit(filtered_data[[first_col]])))
        n_vals <- length(unique_vals)

        if (n_vals == 0) {
          return(tags$p(class = "text-muted", "No values found in selected column"))
        }

        if (n_vals > 30) {
          return(tags$p(class = "text-warning",
                        icon("exclamation-triangle"),
                        paste0(" Too many unique values (", n_vals, "). Use a palette - individual colors not shown.")))
        }

        # v69: Get current palette for default colors
        current_palette <- input[[paste0("heatmap_discrete_palette_", i)]]
        if (is.null(current_palette)) current_palette <- "Set1"

        # Generate default colors from palette
        default_colors <- tryCatch({
          max_colors <- RColorBrewer::brewer.pal.info[current_palette, "maxcolors"]
          if (n_vals <= max_colors) {
            RColorBrewer::brewer.pal(max(3, n_vals), current_palette)[1:n_vals]
          } else {
            # Use colorRampPalette for more colors
            colorRampPalette(RColorBrewer::brewer.pal(max_colors, current_palette))(n_vals)
          }
        }, error = function(e) {
          rainbow(n_vals)
        })

        # S2.4-FIX: Get stored colors from heatmap_configs (used as fallback when input doesn't exist yet)
        stored_colors <- isolate({
          if (i <= length(values$heatmap_configs) && !is.null(values$heatmap_configs[[i]]$custom_colors)) {
            values$heatmap_configs[[i]]$custom_colors
          } else {
            list()
          }
        })

        # S2.6-DEBUG: Log stored colors during UI rebuild
        cat(file=stderr(), paste0("\n[S2.6-DEBUG] renderUI for discrete colors heatmap ", i, "\n"))
        cat(file=stderr(), paste0("[S2.6-DEBUG] length(stored_colors) = ", length(stored_colors), "\n"))
        if (length(stored_colors) > 0) {
          cat(file=stderr(), paste0("[S2.6-DEBUG] stored_colors keys: ", paste(names(stored_colors), collapse=", "), "\n"))
          cat(file=stderr(), paste0("[S2.6-DEBUG] stored_colors values: ", paste(unlist(stored_colors), collapse=", "), "\n"))
        }

        # v70: Generate color pickers for each value WITH dropdown menu
        # S2.292dev: Add manual hex input for discrete colors
        color_pickers <- lapply(seq_along(unique_vals), function(j) {
          val <- as.character(unique_vals[j])

          # v69: Check if there's an existing custom color for this value
          # S2.4-FIX: Also check stored colors in heatmap_configs as fallback
          # S2.7-FIX: Prioritize stored_color over existing_color
          # stored_color is the source of truth (persisted in heatmap_configs)
          # This prevents UI rebuild from showing stale/default colors
          existing_color <- isolate(input[[paste0("heatmap_", i, "_color_", j)]])
          stored_color <- stored_colors[[val]]
          color_to_use <- if (!is.null(stored_color)) {
            stored_color  # S2.7-FIX: Always use stored color if available (source of truth)
          } else if (!is.null(existing_color)) {
            existing_color  # Fall back to existing input value
          } else {
            default_colors[j]
          }

          fluidRow(
            style = "margin-bottom: 3px;",
            column(3, tags$label(val, style = "padding-top: 5px; font-weight: normal; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;", title = val)),
            column(3,
                   tags$div(style = "display: flex; gap: 3px; align-items: center;",
                     colourInput(paste0("heatmap_", i, "_color_", j), NULL,
                                 value = color_to_use, showColour = "background"),
                     textInput(paste0("heatmap_", i, "_color_hex_", j), NULL,
                               value = color_to_use, width = "80px")
                   )
            ),
            column(3,
                   selectInput(paste0("heatmap_", i, "_color_name_", j), NULL,
                               choices = heat_r_colors, selected = "")
            ),
            column(3)
          )
        })

        # v70: NA color picker (always shown at the end)
        # S2.4-FIX: Also check stored NA color in heatmap_configs as fallback
        # S2.7-FIX: Prioritize stored_na_color over existing_na_color (same pattern as value colors)
        existing_na_color <- isolate(input[[paste0("heatmap_", i, "_na_color")]])
        stored_na_color <- isolate({
          if (i <= length(values$heatmap_configs) && !is.null(values$heatmap_configs[[i]]$na_color)) {
            values$heatmap_configs[[i]]$na_color
          } else {
            NULL
          }
        })
        na_color_to_use <- if (!is.null(stored_na_color)) {
          stored_na_color  # S2.7-FIX: Stored color takes priority
        } else if (!is.null(existing_na_color)) {
          existing_na_color
        } else {
          "white"
        }

        # S2.292dev: Add hex input to NA color row
        na_color_row <- fluidRow(
          style = "margin-bottom: 3px; background-color: #f8f8f8; padding: 5px; border-radius: 3px; margin-top: 10px;",
          column(3, tags$label("NA / Missing", style = "padding-top: 5px; font-weight: bold; font-style: italic;")),
          column(3,
                 tags$div(style = "display: flex; gap: 3px; align-items: center;",
                   colourInput(paste0("heatmap_", i, "_na_color"), NULL,
                               value = na_color_to_use, showColour = "background"),
                   textInput(paste0("heatmap_", i, "_na_color_hex"), NULL,
                             value = na_color_to_use, width = "80px")
                 )
          ),
          column(3,
                 selectInput(paste0("heatmap_", i, "_na_color_name"), NULL,
                             choices = heat_r_colors, selected = "")
          ),
          column(3)
        )

        # v104: Restructured with z-index fixes for dropdown menus
        # All dropdowns now float above other elements when opened
        tags$div(
          style = "padding: 10px; background: white; border-radius: 3px; border: 1px solid #ddd; position: relative;",
          # v104: CSS to make selectize dropdowns float above everything
          tags$style(HTML(paste0("
            #heatmap_discrete_colors_ui_", i, " .selectize-dropdown {
              z-index: 10000 !important;
              position: absolute !important;
            }
            #heatmap_discrete_colors_ui_", i, " .selectize-control {
              position: relative;
            }
          "))),
          tags$small(class = "text-muted", paste0(n_vals, " unique value(s) + NA color")),
          tags$hr(style = "margin: 5px 0;"),
          # v70: Header row
          fluidRow(
            style = "margin-bottom: 5px; font-weight: bold; font-size: 11px;",
            column(4, "Value"),
            column(4, "Color"),
            column(4, "R Color Name")
          ),
          # v104: Removed max-height constraint and overflow:auto to prevent dropdown clipping
          # The scrollable container was causing dropdown menus to be hidden
          tags$div(
            style = "max-height: 350px; overflow-y: visible; overflow-x: visible;",
            do.call(tagList, color_pickers)
          ),
          # v103: NA color row outside scrollable area - dropdown won't be clipped
          na_color_row
        )
      })
    })
  })

  # v70: Observer to update heatmap color pickers when R color name dropdown changes
  observe({
    lapply(1:6, function(i) {
      # Observe changes to color name dropdowns for each value (up to 30 values)
      lapply(1:30, function(j) {
        observeEvent(input[[paste0("heatmap_", i, "_color_name_", j)]], {
          color_name <- input[[paste0("heatmap_", i, "_color_name_", j)]]
          if (!is.null(color_name) && color_name != "") {
            updateColourInput(session, paste0("heatmap_", i, "_color_", j), value = color_name)
          }
        }, ignoreInit = TRUE)
      })

      # v70: Also observe NA color name dropdown
      observeEvent(input[[paste0("heatmap_", i, "_na_color_name")]], {
        color_name <- input[[paste0("heatmap_", i, "_na_color_name")]]
        if (!is.null(color_name) && color_name != "") {
          updateColourInput(session, paste0("heatmap_", i, "_na_color"), value = color_name)
        }
      }, ignoreInit = TRUE)
    })
  })

  # S2.4-FIX: Observer to save discrete color picker values to heatmap_configs
  # This fixes the bug where adding a second heatmap resets the first heatmap's colors
  # The issue was that colors were only stored in Shiny inputs, which get destroyed/recreated
  # when the UI rebuilds. Now we persist them in heatmap_configs.
  # S2.7-FIX: Added protection against UI rebuild overwriting custom colors with defaults.
  # When UI rebuilds, color pickers may fire change events with default palette colors.
  # We now detect when incoming color is a default palette color that would overwrite
  # an existing custom color, and skip saving in that case.
  # S2.8-FIX: Added inhibit_color_save flag to completely block saving during UI rebuilds.
  # This is more robust than color comparison because it prevents ALL saves during the
  # critical window when UI is rebuilding.
  observe({
    lapply(1:6, function(i) {
      lapply(1:30, function(j) {
        observeEvent(input[[paste0("heatmap_", i, "_color_", j)]], {
          # S2.8-FIX: Check inhibit flag first - skip ALL saves during UI rebuild
          inhibit_val <- isolate(inhibit_color_save())
          color_val <- input[[paste0("heatmap_", i, "_color_", j)]]
          cat(file=stderr(), paste0("[S2.11-DEBUG] Color observer fired: heatmap=", i,
                                   ", color=", j, ", inhibit=", inhibit_val,
                                   ", value=", ifelse(is.null(color_val), "NULL", color_val), "\n"))
          if (isTRUE(inhibit_val)) {
            cat(file=stderr(), paste0("[S2.8-DEBUG] BLOCKED color save for heatmap ", i,
                                     " color ", j, " (UI rebuild in progress)\n"))
            return(NULL)
          }

          if (i <= length(values$heatmap_configs)) {
            # Get current column to know the unique values
            cols_selected <- input[[paste0("heatmap_columns_", i)]]
            if (!is.null(cols_selected) && length(cols_selected) > 0 && !is.null(values$csv_data)) {
              first_col <- cols_selected[1]
              if (first_col %in% names(values$csv_data)) {
                # Get filtered unique values (same logic as renderUI)
                filtered_data <- values$csv_data
                if (!is.null(values$tree) && !is.null(input$id_column) && input$id_column %in% names(values$csv_data)) {
                  tree_tips <- values$tree$tip.label
                  id_col_data <- as.character(values$csv_data[[input$id_column]])
                  matching_rows <- id_col_data %in% tree_tips
                  if (any(matching_rows)) {
                    filtered_data <- values$csv_data[matching_rows, , drop = FALSE]
                  }
                }
                unique_vals <- sort(unique(na.omit(filtered_data[[first_col]])))

                if (j <= length(unique_vals)) {
                  # Initialize custom_colors if needed
                  if (is.null(values$heatmap_configs[[i]]$custom_colors)) {
                    values$heatmap_configs[[i]]$custom_colors <- list()
                  }
                  # Store the color with the value name as key
                  val_name <- as.character(unique_vals[j])
                  color_value <- input[[paste0("heatmap_", i, "_color_", j)]]

                  if (!is.null(color_value)) {
                    # S2.7-FIX: Check if this is a default palette color that would overwrite a custom color
                    # Get current palette to compute default colors
                    current_palette <- isolate(input[[paste0("heatmap_discrete_palette_", i)]])
                    if (is.null(current_palette)) current_palette <- "Set1"

                    # Generate default colors for comparison
                    n_vals <- length(unique_vals)
                    default_colors <- tryCatch({
                      max_colors <- RColorBrewer::brewer.pal.info[current_palette, "maxcolors"]
                      if (n_vals <= max_colors) {
                        RColorBrewer::brewer.pal(max(3, n_vals), current_palette)[1:n_vals]
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(max_colors, current_palette))(n_vals)
                      }
                    }, error = function(e) {
                      rainbow(n_vals)
                    })

                    default_color_for_j <- if (j <= length(default_colors)) default_colors[j] else "#808080"
                    existing_stored_color <- values$heatmap_configs[[i]]$custom_colors[[val_name]]

                    # S2.9-FIX: Skip saving if color is already stored (no-op during UI rebuild)
                    # This additional check helps when renderUI creates widgets with stored values
                    if (!is.null(existing_stored_color) && toupper(color_value) == toupper(existing_stored_color)) {
                      # Color is already saved with this exact value, skip redundant save
                      return(NULL)
                    }

                    # S2.7-FIX: Skip saving if:
                    # 1. The incoming color is the default palette color for this position
                    # 2. AND there's already a stored custom color that's different
                    # This prevents UI rebuild from resetting custom colors to defaults
                    is_default_color <- toupper(color_value) == toupper(default_color_for_j)
                    has_different_stored <- !is.null(existing_stored_color) &&
                                            toupper(existing_stored_color) != toupper(default_color_for_j)

                    if (is_default_color && has_different_stored) {
                      # S2.7-DEBUG: Log skipped saves
                      cat(file=stderr(), paste0("[S2.7-DEBUG] SKIPPED saving default color for heatmap ", i,
                                               " value '", val_name, "' (default=", default_color_for_j,
                                               ", keeping stored=", existing_stored_color, ")\n"))
                    } else {
                      # S2.11-DEBUG: Trace what's happening with the save
                      cat(file=stderr(), paste0("[S2.11-DEBUG] Color save for heatmap ", i, " color ", j,
                                               ": val_name='", val_name, "', color_value='", color_value,
                                               "', is_default=", is_default_color,
                                               ", existing_stored='", ifelse(is.null(existing_stored_color), "NULL", existing_stored_color),
                                               "', default_for_j='", default_color_for_j, "'\n"))
                      # Normal save: either it's a custom color, or there's no conflict
                      values$heatmap_configs[[i]]$custom_colors[[val_name]] <- color_value
                      values$heatmap_configs[[i]]$custom_discrete <- TRUE
                      # S2.6-DEBUG: Log when color is saved
                      cat(file=stderr(), paste0("[S2.6-DEBUG] SAVED color for heatmap ", i, " value '", val_name, "' = ", color_value, "\n"))
                    }
                  }
                }
              }
            }
          }
        }, ignoreInit = TRUE)
      })

      # S2.4-FIX: Also save NA color
      observeEvent(input[[paste0("heatmap_", i, "_na_color")]], {
        # S2.8-FIX: Check inhibit flag first
        if (isTRUE(isolate(inhibit_color_save()))) {
          return(NULL)
        }

        if (i <= length(values$heatmap_configs)) {
          na_color <- input[[paste0("heatmap_", i, "_na_color")]]
          if (!is.null(na_color)) {
            values$heatmap_configs[[i]]$na_color <- na_color
          }
        }
      }, ignoreInit = TRUE)

      # S2.292dev: Sync discrete colourInput -> hex textInput
      lapply(1:30, function(j) {
        observeEvent(input[[paste0("heatmap_", i, "_color_", j)]], {
          color_val <- input[[paste0("heatmap_", i, "_color_", j)]]
          if (!is.null(color_val)) {
            updateTextInput(session, paste0("heatmap_", i, "_color_hex_", j), value = toupper(color_val))
          }
        }, ignoreInit = TRUE, priority = -1)
      })

      # S2.292dev: Sync discrete NA colourInput -> hex textInput
      observeEvent(input[[paste0("heatmap_", i, "_na_color")]], {
        color_val <- input[[paste0("heatmap_", i, "_na_color")]]
        if (!is.null(color_val)) {
          updateTextInput(session, paste0("heatmap_", i, "_na_color_hex"), value = toupper(color_val))
        }
      }, ignoreInit = TRUE, priority = -1)

      # S2.292dev: Sync discrete hex textInput -> colourInput
      lapply(1:30, function(j) {
        observeEvent(input[[paste0("heatmap_", i, "_color_hex_", j)]], {
          hex_val <- input[[paste0("heatmap_", i, "_color_hex_", j)]]
          if (!is.null(hex_val) && nchar(trimws(hex_val)) > 0) {
            hex_val <- trimws(hex_val)
            if (!grepl("^#", hex_val)) hex_val <- paste0("#", hex_val)
            if (grepl("^#[0-9A-Fa-f]{6}$", hex_val)) {
              updateColourInput(session, paste0("heatmap_", i, "_color_", j), value = hex_val)
            }
          }
        }, ignoreInit = TRUE)
      })

      # S2.292dev: Sync discrete NA hex textInput -> colourInput
      observeEvent(input[[paste0("heatmap_", i, "_na_color_hex")]], {
        hex_val <- input[[paste0("heatmap_", i, "_na_color_hex")]]
        if (!is.null(hex_val) && nchar(trimws(hex_val)) > 0) {
          hex_val <- trimws(hex_val)
          if (!grepl("^#", hex_val)) hex_val <- paste0("#", hex_val)
          if (grepl("^#[0-9A-Fa-f]{6}$", hex_val)) {
            updateColourInput(session, paste0("heatmap_", i, "_na_color"), value = hex_val)
          }
        }
      }, ignoreInit = TRUE)
    })
  })

  # v69: Observer for "Apply Palette to All" buttons for heatmaps
  observe({
    lapply(1:6, function(i) {
      observeEvent(input[[paste0("apply_heatmap_palette_", i)]], {
        cols_selected <- input[[paste0("heatmap_columns_", i)]]
        if (is.null(cols_selected) || length(cols_selected) == 0 || is.null(values$csv_data)) {
          showNotification("Please select a column first", type = "warning")
          return()
        }

        # Get unique values
        first_col <- cols_selected[1]
        if (!(first_col %in% names(values$csv_data))) return()

        unique_vals <- sort(unique(na.omit(values$csv_data[[first_col]])))
        n_vals <- length(unique_vals)

        if (n_vals == 0 || n_vals > 30) return()

        # Get selected palette
        palette_name <- input[[paste0("heatmap_discrete_palette_", i)]]
        if (is.null(palette_name)) palette_name <- "Set1"

        # Generate colors from palette
        new_colors <- tryCatch({
          max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
          if (n_vals <= max_colors) {
            RColorBrewer::brewer.pal(max(3, n_vals), palette_name)[1:n_vals]
          } else {
            colorRampPalette(RColorBrewer::brewer.pal(max_colors, palette_name))(n_vals)
          }
        }, error = function(e) {
          rainbow(n_vals)
        })

        # Update all color pickers and hex inputs
        for (j in seq_along(unique_vals)) {
          updateColourInput(session, paste0("heatmap_", i, "_color_", j), value = new_colors[j])
          updateTextInput(session, paste0("heatmap_", i, "_color_hex_", j), value = toupper(new_colors[j]))
        }

        showNotification(paste("Applied", palette_name, "palette to", n_vals, "values"), type = "message")
      }, ignoreInit = TRUE)
    })
  })
  
  # Apply heatmaps button
  observeEvent(input$apply_heatmaps, {
    cat(file=stderr(), "\n[HEATMAP-APPLY] Apply Heatmaps button clicked!\n")

    # Convert heatmap_configs to the format expected by the plotting function
    if (length(values$heatmap_configs) == 0) {
      showNotification("No heatmaps configured", type = "warning")
      return()
    }

    cat(file=stderr(), paste0("[HEATMAP-APPLY] Processing ", length(values$heatmap_configs), " heatmap config(s)\n"))

    # S2.3-DEBUG: Show all heatmap column selections upfront
    cat(file=stderr(), "\n[HEATMAP-APPLY-DEBUG] === HEATMAP COLUMN MAPPING ===\n")
    for (debug_i in seq_along(values$heatmap_configs)) {
      debug_cols <- input[[paste0("heatmap_columns_", debug_i)]]
      cat(file=stderr(), paste0("[HEATMAP-APPLY-DEBUG] Heatmap ", debug_i,
                                 " columns from input: ",
                                 if (is.null(debug_cols)) "NULL" else paste(debug_cols, collapse=", "), "\n"))
    }
    cat(file=stderr(), "[HEATMAP-APPLY-DEBUG] ==============================\n\n")

    # v56a: Build heatmaps list from configs with multiple column support
    # Read directly from inputs to ensure we get current values (fixes ignoreInit issue)
    heatmaps_list <- lapply(seq_along(values$heatmap_configs), function(i) {
      cfg <- values$heatmap_configs[[i]]

      # S1.62dev: Check data source - CSV columns or RData CNV
      current_data_source <- input[[paste0("heatmap_data_source_", i)]]
      if (is.null(current_data_source)) current_data_source <- "csv"
      cat(file=stderr(), paste0("[HEATMAP-APPLY] Heatmap ", i, " data_source: '", current_data_source, "'\n"))

      # S1.62dev: Handle RData CNV source
      if (current_data_source == "rdata") {
        cat(file=stderr(), "[HEATMAP-APPLY] Entering RData CNV path\n")
        # Check if RData CNV matrix is available
        if (is.null(values$rdata_cnv_matrix)) {
          showNotification(paste("Heatmap", i, ": No RData CNV file loaded. Please upload an RData file first."), type = "error")
          return(NULL)
        }

        # Get CNV settings - these will be stored and applied later in func.print.lineage.tree
        # Use new field name, but keep old one for backwards compatibility
        cnv_render_downsample <- input[[paste0("heatmap_cnv_render_downsample_", i)]]
        if (is.null(cnv_render_downsample)) cnv_render_downsample <- 10
        cnv_downsample <- cnv_render_downsample  # For backwards compatibility
        cnv_wgd_norm <- input[[paste0("heatmap_cnv_wgd_norm_", i)]]
        if (is.null(cnv_wgd_norm)) cnv_wgd_norm <- FALSE
        # S2.12: Per-cell WGD settings
        cnv_wgd_per_cell <- input[[paste0("heatmap_cnv_wgd_per_cell_", i)]]
        if (is.null(cnv_wgd_per_cell)) cnv_wgd_per_cell <- FALSE
        cnv_wgd_column <- input[[paste0("heatmap_cnv_wgd_column_", i)]]
        if (is.null(cnv_wgd_column)) cnv_wgd_column <- NULL

        debug_cat(paste0("  S1.62dev: RData CNV heatmap ", i, ": render_downsample=", cnv_render_downsample, ", wgd_norm=", cnv_wgd_norm, "\n"))
        cat(file=stderr(), paste0("[RENDER-DOWNSAMPLE] Heatmap ", i, " render_downsample value: ", cnv_render_downsample, "\n"))
        debug_cat(paste0("  S2.12: Per-cell WGD: enabled=", cnv_wgd_per_cell, ", column=", ifelse(is.null(cnv_wgd_column), "NULL", cnv_wgd_column), "\n"))
        debug_cat(paste0("  S1.62dev: Raw CNV matrix: ", nrow(values$rdata_cnv_matrix), " samples x ", ncol(values$rdata_cnv_matrix), " positions\n"))

        # Build heatmap entry for RData CNV
        # Note: cnv_matrix is NOT stored here - it's passed as a parameter to func.print.lineage.tree
        # because large matrices don't serialize properly to YAML
        cat(file=stderr(), paste0("[DEBUG-COLOR] Building RData heatmap entry for heatmap ", i, "\n"))
        cat(file=stderr(), paste0("[DEBUG-COLOR] cfg$low_color = ", ifelse(is.null(cfg$low_color), "NULL", cfg$low_color), "\n"))
        cat(file=stderr(), paste0("[DEBUG-COLOR] cfg$mid_color = ", ifelse(is.null(cfg$mid_color), "NULL", cfg$mid_color), "\n"))
        cat(file=stderr(), paste0("[DEBUG-COLOR] cfg$high_color = ", ifelse(is.null(cfg$high_color), "NULL", cfg$high_color), "\n"))
        cat(file=stderr(), paste0("[DEBUG-COLLINES] show_col_lines input = ", ifelse(is.null(input[[paste0("heatmap_show_col_lines_", i)]]), "NULL", input[[paste0("heatmap_show_col_lines_", i)]]), "\n"))
        # S2.0: Get the mapping column (user-selected CSV column for sample name mapping)
        mapping_column <- input[[paste0("heatmap_rdata_mapping_col_", i)]]
        if (is.null(mapping_column) || mapping_column == "") {
          # Fall back to values$rdata_mapping_column if not set for this heatmap
          mapping_column <- values$rdata_mapping_column
        }
        cat(file=stderr(), paste0("[DEBUG-RDATA] Mapping column for heatmap ", i, ": '",
                                  ifelse(is.null(mapping_column), "NULL", mapping_column), "'\n"))

        heatmap_entry <- list(
          title = cfg$title,
          is_discrete = FALSE,  # CNV data is always continuous
          data_source = "rdata",
          # Store CNV settings (processing happens in func.print.lineage.tree)
          cnv_downsample = cnv_downsample,  # Backwards compatibility
          cnv_render_downsample = cnv_render_downsample,  # New field name (default 0)
          cnv_wgd_norm = cnv_wgd_norm,
          # S2.12: Per-cell WGD normalization settings
          cnv_wgd_per_cell = cnv_wgd_per_cell,
          cnv_wgd_column = cnv_wgd_column,
          # S2.8: Display mode (basic = geom_tile, detailed = geom_raster like pheatmap)
          cnv_display_mode = if (!is.null(input[[paste0("heatmap_cnv_display_mode_", i)]])) input[[paste0("heatmap_cnv_display_mode_", i)]] else "basic",
          # Height scale for detailed mode
          cnv_height_scale = if (!is.null(input[[paste0("heatmap_cnv_height_scale_", i)]])) input[[paste0("heatmap_cnv_height_scale_", i)]] else 1.0,
          # S2.0: Store mapping column for sample name matching
          rdata_mapping_column = mapping_column,
          columns = character(0),  # No columns for RData - data comes from parameter
          show_colnames = if (!is.null(cfg$show_colnames)) cfg$show_colnames else FALSE,  # Usually too many columns
          colnames_angle = if (!is.null(cfg$colnames_angle)) cfg$colnames_angle else 90,
          font_size = input$heatmap_global_font,
          # Distance and height settings
          distance = if (!is.null(input[[paste0("heatmap_distance_", i)]])) input[[paste0("heatmap_distance_", i)]] else 0.02,
          height = if (!is.null(input[[paste0("heatmap_height_", i)]])) input[[paste0("heatmap_height_", i)]] else 0.8,
          row_height = if (!is.null(input[[paste0("heatmap_row_height_", i)]])) input[[paste0("heatmap_row_height_", i)]] else 1.0,
          # Grid settings
          show_grid = if (!is.null(input[[paste0("heatmap_show_grid_", i)]])) input[[paste0("heatmap_show_grid_", i)]] else FALSE,
          grid_color = if (!is.null(input[[paste0("heatmap_grid_color_", i)]])) input[[paste0("heatmap_grid_color_", i)]] else "#000000",
          grid_size = if (!is.null(input[[paste0("heatmap_grid_size_", i)]])) input[[paste0("heatmap_grid_size_", i)]] else 0.5,
          # S1.62dev: Row line settings (horizontal lines only)
          show_row_lines = if (!is.null(input[[paste0("heatmap_show_row_lines_", i)]])) input[[paste0("heatmap_show_row_lines_", i)]] else FALSE,
          row_line_color = if (!is.null(input[[paste0("heatmap_row_line_color_", i)]])) input[[paste0("heatmap_row_line_color_", i)]] else "#000000",
          row_line_size = if (!is.null(input[[paste0("heatmap_row_line_size_", i)]])) input[[paste0("heatmap_row_line_size_", i)]] else 0.5,
          # S1.62dev: Column line settings (vertical lines)
          show_col_lines = if (!is.null(input[[paste0("heatmap_show_col_lines_", i)]])) input[[paste0("heatmap_show_col_lines_", i)]] else FALSE,
          col_line_color = if (!is.null(input[[paste0("heatmap_col_line_color_", i)]])) input[[paste0("heatmap_col_line_color_", i)]] else "#000000",
          col_line_size = if (!is.null(input[[paste0("heatmap_col_line_size_", i)]])) input[[paste0("heatmap_col_line_size_", i)]] else 0.5,
          # Guide lines (usually not useful for CNV with many columns)
          show_guides = FALSE,
          # Row labels
          show_row_labels = if (!is.null(input[[paste0("heatmap_show_row_labels_", i)]])) input[[paste0("heatmap_show_row_labels_", i)]] else FALSE,
          row_label_source = "colnames",
          row_label_font_size = if (!is.null(input[[paste0("heatmap_row_label_font_size_", i)]])) input[[paste0("heatmap_row_label_font_size_", i)]] else 2.5,
          row_label_offset = if (!is.null(input[[paste0("heatmap_row_label_offset_", i)]])) input[[paste0("heatmap_row_label_offset_", i)]] else 1.0,
          row_label_align = if (!is.null(input[[paste0("heatmap_row_label_align_", i)]])) input[[paste0("heatmap_row_label_align_", i)]] else "left",
          # S1.62dev: Color settings - red-white-blue for CNV (red=loss, blue=gain)
          # Note: Both 'low'/'mid'/'high' (for rendering) and 'low_color'/'mid_color'/'high_color' (for UI) are set
          low_color = if (!is.null(cfg$low_color)) cfg$low_color else "#FF0000",   # Red for deletion/loss
          mid_color = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFFFF",   # White for neutral
          high_color = if (!is.null(cfg$high_color)) cfg$high_color else "#0000FF", # Blue for amplification/gain
          low = if (!is.null(cfg$low_color)) cfg$low_color else "#FF0000",   # For rendering code
          mid = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFFFF",   # For rendering code
          high = if (!is.null(cfg$high_color)) cfg$high_color else "#0000FF", # For rendering code
          midpoint = if (!is.null(cfg$midpoint)) cfg$midpoint else 2,  # Center at diploid (2)
          use_midpoint = TRUE,  # Always use midpoint for CNV
          na_color = "grey90"
        )

        return(heatmap_entry)
      }

      # v56a: Read columns directly from input (fixes issue where ignoreInit=TRUE misses initial selection)
      current_columns <- input[[paste0("heatmap_columns_", i)]]
      if (is.null(current_columns) || length(current_columns) == 0) {
        return(NULL)
      }

      # Update config with current columns
      cfg$columns <- current_columns
      cfg$data_source <- "csv"  # Explicitly mark as CSV source

      # v122: Read auto_type from current input (not from stale cfg)
      current_auto_type <- input[[paste0("heatmap_auto_type_", i)]]
      if (is.null(current_auto_type)) current_auto_type <- TRUE

      # v122: Read forced type from current input
      current_forced_type <- input[[paste0("heatmap_type_", i)]]
      if (is.null(current_forced_type)) current_forced_type <- "discrete"

      # v116: Improved auto-detect logic for discrete vs continuous
      # Priority: decimals = continuous, non-numeric = discrete, then check unique values
      actual_type <- current_forced_type  # v122: Default to forced type
      first_col <- cfg$columns[1]
      debug_cat(paste0("  v22 AUTO-DETECT: auto_type=", current_auto_type,
                                 ", forced_type=", current_forced_type, "\n"))
      if (current_auto_type && !is.null(values$csv_data) && first_col %in% names(values$csv_data)) {
        col_data <- values$csv_data[[first_col]]
        col_data_clean <- na.omit(col_data)
        unique_vals <- length(unique(col_data_clean))
        is_numeric <- is.numeric(col_data)
        converted_to_numeric <- FALSE
        originally_numeric <- is_numeric

        # v117: Check for decimal points in string representation FIRST
        # This catches cases like "23.6" that might get converted to numeric
        # Also check originally numeric data by converting to string
        # IMPORTANT: Filter out "NA" strings which are not R's NA
        has_decimal_in_string <- FALSE
        if (is.character(col_data) || is.factor(col_data)) {
          char_data <- as.character(na.omit(col_data))
          # v117: Remove "NA" strings (case-insensitive) that are NOT R's NA
          # v118: Added #N/A to the list of NA-like strings
          char_data <- char_data[!toupper(trimws(char_data)) %in% c("NA", "N/A", "#N/A", "NULL", "")]
          if (length(char_data) > 0) {
            has_decimal_in_string <- any(grepl("\\.[0-9]+", char_data))
          }
        } else if (is.numeric(col_data)) {
          # v116: For originally numeric data, convert to string and check for decimals
          # This catches values like 15.7 that are already loaded as numeric
          char_data <- as.character(na.omit(col_data))
          has_decimal_in_string <- any(grepl("\\.[0-9]+", char_data))
        }
        debug_cat(paste0("  v17 AUTO-DETECT: has_decimal_in_string=", has_decimal_in_string, "\n"))

        # v117: More aggressive conversion - try to convert any non-numeric column to numeric
        # First clean up NA-like strings
        if (!is_numeric) {
          clean_col_data <- as.character(col_data)
          # v117: Convert NA-like strings to actual NA
          # v118: Added #N/A to the list of NA-like strings
          clean_col_data[toupper(trimws(clean_col_data)) %in% c("NA", "N/A", "#N/A", "NULL", "")] <- NA
          numeric_attempt <- suppressWarnings(as.numeric(clean_col_data))
          non_na_original <- sum(!is.na(clean_col_data))  # v117: Count after cleaning NA strings
          non_na_converted <- sum(!is.na(numeric_attempt))
          # v115: Lower threshold to 50% for numeric detection (was 80%)
          if (non_na_original > 0 && (non_na_converted / non_na_original) >= 0.5) {
            is_numeric <- TRUE
            converted_to_numeric <- TRUE
            col_data_clean <- na.omit(numeric_attempt)
            unique_vals <- length(unique(col_data_clean))
          }
        }

        # v117: Debug output for auto-detect troubleshooting
        debug_cat(paste0("  v17 AUTO-DETECT: column=", first_col,
                                   ", is_numeric=", is_numeric,
                                   ", originally_numeric=", originally_numeric,
                                   ", converted=", converted_to_numeric,
                                   ", unique_vals=", unique_vals, "\n"))

        # v117: Better heuristic - decimals ALWAYS mean continuous
        if (!is_numeric) {
          # Non-numeric data is always discrete
          actual_type <- "discrete"
          debug_cat(paste0("  v17 AUTO-DETECT: Result=discrete (non-numeric)\n"))
        } else {
          # v117: Check for decimal values with tolerance for floating-point precision
          epsilon <- 1e-6
          has_decimals <- any(abs(col_data_clean - floor(col_data_clean)) > epsilon, na.rm = TRUE)

          # v117: Also use the string-based detection result
          if (!has_decimals && has_decimal_in_string) {
            has_decimals <- TRUE
            debug_cat(paste0("  v17 AUTO-DETECT: decimal detected via string check\n"))
          }

          debug_cat(paste0("  v20 AUTO-DETECT: has_decimals=", has_decimals, "\n"))

          # v120: Harmonized logic with UI section for consistent detection
          # Get value range for range-based checks
          val_range <- if (length(col_data_clean) > 0) diff(range(col_data_clean, na.rm = TRUE)) else 0

          # v120: Decimals ALWAYS mean continuous (measurements, percentages, etc.)
          if (has_decimals) {
            actual_type <- "continuous"
            debug_cat(paste0("  v20 AUTO-DETECT: Result=continuous (has decimals)\n"))
          } else if (originally_numeric) {
            # v120: Originally numeric without decimals - match UI section logic
            is_boolean_like <- unique_vals <= 2 && val_range <= 1
            is_small_categorical <- unique_vals <= 3 && val_range <= 2
            if (is_boolean_like || is_small_categorical) {
              actual_type <- "discrete"
              debug_cat(paste0("  v20 AUTO-DETECT: Result=discrete (boolean-like or small categorical)\n"))
            } else {
              actual_type <- "continuous"
              debug_cat(paste0("  v20 AUTO-DETECT: Result=continuous (originally numeric, many values)\n"))
            }
          } else {
            # v123: Converted from character without decimals - more lenient for numeric-like data
            # If >=5 unique values OR range >=5 → treat as continuous
            # This catches integer sequences like 1,2,3,4,5 which are typically counts/scores
            if (unique_vals >= 5 || val_range >= 5) {
              actual_type <- "continuous"
              debug_cat(paste0("  v23 AUTO-DETECT: Result=continuous (>=5 unique or range>=5)\n"))
            } else {
              actual_type <- "discrete"
              debug_cat(paste0("  v23 AUTO-DETECT: Result=discrete (few unique, narrow range)\n"))
            }
          }
        }
      }

      # S1.62dev: Store the computed actual_type in config for export
      if (i <= length(values$heatmap_configs)) {
        values$heatmap_configs[[i]]$actual_type <- actual_type
      }

      # v105/v111: Read per-heatmap settings
      current_distance <- input[[paste0("heatmap_distance_", i)]]
      current_height <- input[[paste0("heatmap_height_", i)]]
      current_row_height <- input[[paste0("heatmap_row_height_", i)]]  # v111: Add row_height
      # v111: Grid settings
      show_grid <- input[[paste0("heatmap_show_grid_", i)]]
      grid_color <- input[[paste0("heatmap_grid_color_", i)]]
      grid_size <- input[[paste0("heatmap_grid_size_", i)]]
      # v116: Guide line settings
      show_guides <- input[[paste0("heatmap_show_guides_", i)]]
      guide_color1 <- input[[paste0("heatmap_guide_color1_", i)]]
      guide_color2 <- input[[paste0("heatmap_guide_color2_", i)]]
      guide_alpha <- input[[paste0("heatmap_guide_alpha_", i)]]
      guide_width <- input[[paste0("heatmap_guide_width_", i)]]
      guide_linetype <- input[[paste0("heatmap_guide_linetype_", i)]]
      show_row_labels <- input[[paste0("heatmap_show_row_labels_", i)]]
      row_label_source <- input[[paste0("heatmap_row_label_source_", i)]]
      row_label_font_size <- input[[paste0("heatmap_row_label_font_size_", i)]]
      row_label_offset <- input[[paste0("heatmap_row_label_offset_", i)]]  # v111: Label offset
      row_label_align <- input[[paste0("heatmap_row_label_align_", i)]]  # v111: Label alignment
      custom_row_labels <- input[[paste0("heatmap_custom_row_labels_", i)]]

      # v108: Collect label mapping if using "mapping" source
      label_mapping <- list()
      if (!is.null(row_label_source) && row_label_source == "mapping" && length(cfg$columns) > 0) {
        for (j in seq_along(cfg$columns)) {
          mapped_label <- input[[paste0("heatmap_", i, "_label_", j)]]
          if (!is.null(mapped_label) && nchar(trimws(mapped_label)) > 0) {
            label_mapping[[cfg$columns[j]]] <- mapped_label
          } else {
            label_mapping[[cfg$columns[j]]] <- cfg$columns[j]  # Default to column name
          }
        }
      }

      heatmap_entry <- list(
        title = cfg$title,
        is_discrete = (actual_type == "discrete"),
        columns = cfg$columns,  # v56: Now supports multiple columns
        show_colnames = cfg$show_colnames,
        colnames_angle = if (!is.null(cfg$colnames_angle)) cfg$colnames_angle else 45,
        font_size = input$heatmap_global_font,
        distance = if (!is.null(current_distance)) current_distance else 0.02,
        height = if (!is.null(current_height)) current_height else 0.8,  # v105: Per-heatmap height (column width)
        row_height = if (!is.null(current_row_height)) current_row_height else 1.0,  # v111: Per-heatmap row height
        show_grid = if (!is.null(show_grid)) show_grid else FALSE,  # v111: Grid around tiles
        grid_color = if (!is.null(grid_color)) grid_color else "#000000",  # v111
        grid_size = if (!is.null(grid_size)) grid_size else 0.5,  # v111
        # S1.62dev: Row line settings (horizontal lines only)
        show_row_lines = if (!is.null(input[[paste0("heatmap_show_row_lines_", i)]])) input[[paste0("heatmap_show_row_lines_", i)]] else FALSE,
        row_line_color = if (!is.null(input[[paste0("heatmap_row_line_color_", i)]])) input[[paste0("heatmap_row_line_color_", i)]] else "#000000",
        row_line_size = if (!is.null(input[[paste0("heatmap_row_line_size_", i)]])) input[[paste0("heatmap_row_line_size_", i)]] else 0.5,
        # S1.62dev: Column line settings (vertical lines)
        show_col_lines = if (!is.null(input[[paste0("heatmap_show_col_lines_", i)]])) input[[paste0("heatmap_show_col_lines_", i)]] else FALSE,
        col_line_color = if (!is.null(input[[paste0("heatmap_col_line_color_", i)]])) input[[paste0("heatmap_col_line_color_", i)]] else "#000000",
        col_line_size = if (!is.null(input[[paste0("heatmap_col_line_size_", i)]])) input[[paste0("heatmap_col_line_size_", i)]] else 0.5,
        # v116: Guide lines
        show_guides = if (!is.null(show_guides)) show_guides else FALSE,
        guide_color1 = if (!is.null(guide_color1)) guide_color1 else "#CCCCCC",
        guide_color2 = if (!is.null(guide_color2)) guide_color2 else "#EEEEEE",
        guide_alpha = if (!is.null(guide_alpha)) guide_alpha else 0.3,
        guide_width = if (!is.null(guide_width)) guide_width else 0.5,
        guide_linetype = if (!is.null(guide_linetype)) guide_linetype else "solid",
        show_row_labels = if (!is.null(show_row_labels)) show_row_labels else FALSE,
        row_label_source = if (!is.null(row_label_source)) row_label_source else "colnames",
        row_label_font_size = if (!is.null(row_label_font_size)) row_label_font_size else 2.5,
        row_label_offset = if (!is.null(row_label_offset)) row_label_offset else 1.0,  # v111
        row_label_align = if (!is.null(row_label_align)) row_label_align else "left",  # v111
        custom_row_labels = if (!is.null(custom_row_labels)) custom_row_labels else "",
        label_mapping = label_mapping  # v108: Per-column label mapping
      )
      
      # S2.11-DEBUG: Trace through discrete color handling
      cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " actual_type='", actual_type, "'\n"))

      if (actual_type == "discrete") {
        # v69: Get palette from current input (not stale cfg)
        # S1.62dev: Fall back to cfg$discrete_palette if input not available (e.g., after YAML import)
        current_palette <- input[[paste0("heatmap_discrete_palette_", i)]]
        heatmap_entry$color_scheme <- if (!is.null(current_palette)) {
          current_palette
        } else if (!is.null(cfg$discrete_palette)) {
          cfg$discrete_palette
        } else {
          "Set1"
        }

        # S2.11-DEBUG: Trace first_col and csv_data availability
        cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " first_col='", first_col, "'\n"))
        cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " csv_data is NULL: ", is.null(values$csv_data), "\n"))
        cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " first_col in names: ", first_col %in% names(values$csv_data), "\n"))

        # v69: Collect custom colors if they've been set
        # S2.1-FIX: Use filtered data (matching tree tips) to match renderUI color picker creation
        # Bug fix: Previously used unfiltered values$csv_data which caused wrong value-color mappings
        if (!is.null(values$csv_data) && first_col %in% names(values$csv_data)) {
          # S2.1-FIX: Filter to only patient-specific data (rows matching tree tips)
          # This MUST match the filtering in heatmap_discrete_colors_ui renderUI (lines 15342-15353)
          filtered_data <- values$csv_data
          if (!is.null(values$tree) && !is.null(input$id_column) && input$id_column %in% names(values$csv_data)) {
            tree_tips <- values$tree$tip.label
            id_col_data <- as.character(values$csv_data[[input$id_column]])
            matching_rows <- id_col_data %in% tree_tips
            if (any(matching_rows)) {
              filtered_data <- values$csv_data[matching_rows, , drop = FALSE]
              debug_cat(paste0("S2.1-FIX: Filtered data for color collection from ", nrow(values$csv_data), " to ", nrow(filtered_data), " rows (tree tips)\n"))
            }
          }

          unique_vals <- sort(unique(na.omit(filtered_data[[first_col]])))
          n_vals <- length(unique_vals)

          # S2.11-DEBUG: Trace unique values
          cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " unique_vals: ", paste(unique_vals, collapse=", "), "\n"))
          cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " n_vals=", n_vals, "\n"))

          if (n_vals > 0 && n_vals <= 30) {
            custom_colors <- c()
            has_custom_colors <- FALSE

            # S2.10-FIX: Prioritize stored colors from heatmap_configs over input values
            # This fixes the bug where UI rebuild during "Add Heatmap" could cause
            # input color pickers to have stale/default values, overwriting custom colors.
            # cfg$custom_colors is the source of truth, maintained by the color observer.
            stored_cfg_colors <- cfg$custom_colors

            # S2.11-DEBUG: Trace stored colors from config
            cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " cfg$custom_colors is NULL: ", is.null(stored_cfg_colors), "\n"))
            if (!is.null(stored_cfg_colors) && length(stored_cfg_colors) > 0) {
              cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " cfg$custom_colors keys: ", paste(names(stored_cfg_colors), collapse=", "), "\n"))
              cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " cfg$custom_colors values: ", paste(unlist(stored_cfg_colors), collapse=", "), "\n"))
            }

            for (j in seq_along(unique_vals)) {
              val_name <- as.character(unique_vals[j])

              # S2.10-FIX: First try stored color, then fall back to input
              stored_color <- if (!is.null(stored_cfg_colors)) stored_cfg_colors[[val_name]] else NULL
              color_input <- input[[paste0("heatmap_", i, "_color_", j)]]

              # S2.11-DEBUG: Trace individual color resolution
              cat(file=stderr(), paste0("[S2.11-DEBUG] Heatmap ", i, " val '", val_name, "': stored=",
                  ifelse(is.null(stored_color), "NULL", stored_color),
                  ", input=", ifelse(is.null(color_input), "NULL", color_input), "\n"))

              # Use stored color if available (source of truth), otherwise use input
              color_to_use <- if (!is.null(stored_color)) stored_color else color_input

              if (!is.null(color_to_use)) {
                custom_colors[val_name] <- color_to_use
                has_custom_colors <- TRUE
              }
            }

            # S2.10-DEBUG: Log color sources for debugging
            if (has_custom_colors) {
              cat(file=stderr(), paste0("[S2.10-DEBUG] Heatmap ", i, " custom_colors for Apply:\n"))
              for (val_name in names(custom_colors)) {
                cat(file=stderr(), paste0("[S2.10-DEBUG]   ", val_name, " = ", custom_colors[val_name], "\n"))
              }
            }

            if (has_custom_colors) {
              heatmap_entry$custom_colors <- custom_colors
              heatmap_entry$man_define_colors <- TRUE
            }
          }

          # v70: Get NA color from input (default to white if not set)
          # S2.10-FIX: Also check stored NA color in cfg
          na_color_input <- input[[paste0("heatmap_", i, "_na_color")]]
          stored_na_color <- cfg$na_color
          heatmap_entry$na_color <- if (!is.null(stored_na_color)) {
            stored_na_color
          } else if (!is.null(na_color_input)) {
            na_color_input
          } else {
            "white"
          }
        }
      } else {
        heatmap_entry$low_color <- cfg$low_color
        heatmap_entry$high_color <- cfg$high_color
        if (cfg$use_midpoint) {
          heatmap_entry$mid_color <- cfg$mid_color
          heatmap_entry$midpoint <- cfg$midpoint
        }
        # v112: Get NA color for continuous heatmaps from input (default to grey90)
        cont_na_color_input <- input[[paste0("heatmap_", i, "_cont_na_color")]]
        heatmap_entry$na_color <- if (!is.null(cont_na_color_input)) cont_na_color_input else "grey90"
      }

      heatmap_entry
    })
    
    # Remove NULLs
    heatmaps_list <- heatmaps_list[!sapply(heatmaps_list, is.null)]
    
    if (length(heatmaps_list) == 0) {
      showNotification("No valid heatmaps configured (select columns)", type = "warning")
      return()
    }
    
    values$heatmaps <- heatmaps_list

    # v125: Store global gap setting for use in plotting
    values$heatmap_global_gap <- if (!is.null(input$heatmap_global_gap)) input$heatmap_global_gap else 0.05

    # v56c: DEBUG - Print heatmap structure being applied
    debug_cat("\n=== v56c HEATMAP APPLY DEBUG ===\n")
    debug_cat("Number of heatmaps:", length(heatmaps_list), "\n")
    for (h_idx in seq_along(heatmaps_list)) {
      hm <- heatmaps_list[[h_idx]]
      debug_cat(paste0("  Heatmap ", h_idx, ":\n"))
      debug_cat(paste0("    Title: ", hm$title, "\n"))
      debug_cat(paste0("    Columns: ", paste(hm$columns, collapse=", "), "\n"))
      debug_cat(paste0("    Is discrete: ", hm$is_discrete, "\n"))
      # S1.62dev: Add data_source debug
      debug_cat(paste0("    S1.62dev data_source: ", if (!is.null(hm$data_source)) hm$data_source else "NULL", "\n"))
      debug_cat(paste0("    S1.62dev all fields: ", paste(names(hm), collapse=", "), "\n"))
    }
    debug_cat("================================\n\n")

    # Generate plot
    generate_plot()

    showNotification(paste(length(heatmaps_list), "heatmap(s) applied"), type = "message")
  })
  
  # ============================================
  # END v55: NEW MULTI-HEATMAP SYSTEM
  # ============================================

  # ============================================
  # v121: LEGEND SETTINGS SYSTEM
  # ============================================

  # Observer for Apply Legend Settings button
  # v180: Updated with new controls (key width/height, byrow, box background, margin)
  observeEvent(input$apply_legend_settings, {
    debug_cat("\n=== v180: APPLYING LEGEND SETTINGS ===\n")

    # Update legend settings in reactive values
    # S2.292dev: Added font_family setting
    values$legend_settings <- list(
      position = input$legend_position,
      # Visibility controls
      show_classification = input$legend_show_classification,
      show_highlight = input$legend_show_highlight,
      show_bootstrap = input$legend_show_bootstrap,
      show_pvalue = input$legend_show_pvalue,
      show_heatmap = input$legend_show_heatmap,  # v180: Simplified to single checkbox
      # Font settings
      title_size = input$legend_title_size,
      text_size = input$legend_text_size,
      font_family = if (!is.null(input$legend_font_family)) input$legend_font_family else "sans",
      key_size = input$legend_key_size,
      # v180: Key dimensions
      key_width = input$legend_key_width,
      key_height = input$legend_key_height,
      # Spacing controls
      spacing = input$legend_spacing,
      spacing_vertical = input$legend_spacing_vertical,
      title_key_spacing = input$legend_title_key_spacing,  # v180: Title to keys spacing
      key_spacing = input$legend_key_spacing,              # v180: Between keys spacing
      # Layout controls
      reverse_order = input$legend_reverse_order,
      # v180: Background controls
      box_background = input$legend_box_background,
      margin = input$legend_margin
    )

    debug_cat(paste0("  Position: ", input$legend_position, "\n"))
    debug_cat(paste0("  Show classification: ", input$legend_show_classification, "\n"))
    debug_cat(paste0("  Show highlight: ", input$legend_show_highlight, "\n"))
    debug_cat(paste0("  Show bootstrap: ", input$legend_show_bootstrap, "\n"))
    debug_cat(paste0("  Show P value: ", input$legend_show_pvalue, "\n"))
    debug_cat(paste0("  Show heatmap: ", input$legend_show_heatmap, "\n"))
    debug_cat(paste0("  S2.292dev: Font family: ", input$legend_font_family, "\n"))
    cat(file=stderr(), paste0("[LEGEND-SETTINGS] Font family set to: '", input$legend_font_family, "'\n"))
    cat(file=stderr(), paste0("[LEGEND-SETTINGS] Stored in values$legend_settings$font_family: '", values$legend_settings$font_family, "'\n"))
    debug_cat(paste0("  v80: Key width: ", input$legend_key_width, ", height: ", input$legend_key_height, "\n"))
    debug_cat(paste0("  v80: Title-key spacing: ", input$legend_title_key_spacing, "\n"))
    debug_cat(paste0("  v80: Between keys spacing: ", input$legend_key_spacing, "\n"))
    debug_cat(paste0("  v80: Reverse order: ", input$legend_reverse_order, "\n"))
    debug_cat(paste0("  v80: Box background: ", input$legend_box_background, ", Margin: ", input$legend_margin, "\n"))
    debug_cat("======================================\n\n")

    # Regenerate plot with new legend settings
    generate_plot()

    showNotification("Legend settings applied", type = "message")
  })

  # ============================================
  # END v121: LEGEND SETTINGS SYSTEM
  # ============================================

  # v139: Observer to swap width/height when page orientation changes
  # Changed to always swap dimensions when orientation changes (not just when they don't match)
  # v141: Also regenerate plot to show the orientation change visually
  # v180: Fixed "keep_proportions" - always swap PAGE dimensions, but plot stays proportional
  observeEvent(input$page_orientation, {
    current_width <- isolate(input$output_width)
    current_height <- isolate(input$output_height)
    keep_proportions <- isolate(input$keep_proportions)

    debug_cat(paste0("\n=== v180: PAGE ORIENTATION CHANGED ===\n"))
    debug_cat(paste0("  Orientation: ", input$page_orientation, "\n"))
    debug_cat(paste0("  Keep proportions: ", keep_proportions, "\n"))
    debug_cat(paste0("  Current width: ", current_width, ", height: ", current_height, "\n"))

    # v180: ALWAYS swap page dimensions when orientation changes
    # The keep_proportions option controls whether the PLOT stretches to fill the new page
    # (handled separately in the rendering code with cowplot positioning)
    if (!is.null(current_width) && !is.null(current_height)) {
      if (input$page_orientation == "landscape") {
        # Landscape: width should be > height
        if (current_width < current_height) {
          debug_cat(paste0("  Swapping to landscape: width=", current_height, ", height=", current_width, "\n"))
          updateNumericInput(session, "output_width", value = current_height)
          updateNumericInput(session, "output_height", value = current_width)
        } else {
          debug_cat("  Already in landscape orientation (width >= height)\n")
        }
      } else {
        # Portrait: height should be > width
        if (current_width > current_height) {
          debug_cat(paste0("  Swapping to portrait: width=", current_height, ", height=", current_width, "\n"))
          updateNumericInput(session, "output_width", value = current_height)
          updateNumericInput(session, "output_height", value = current_width)
        } else {
          debug_cat("  Already in portrait orientation (height >= width)\n")
        }
      }
    } else {
      debug_cat("  WARNING: width or height is NULL, cannot swap\n")
    }

    # v180: Log the effect of keep_proportions
    if (isTRUE(keep_proportions)) {
      debug_cat("  v180: Plot proportions preserved - plot won't stretch to fill new page\n")
    }
    debug_cat("=== END PAGE ORIENTATION ===\n")

    # v143: Regenerate plot after a small delay to ensure new dimensions are available
    # updateNumericInput doesn't update input$ values immediately - they need a flush cycle
    if (!is.null(values$tree_data)) {
      shinyjs::delay(100, {
        generate_plot()
        debug_cat("  v143: Plot regenerated for orientation preview (after delay)\n")
      })
    }
  }, ignoreInit = TRUE)

  # Update Preview (without saving)
  observeEvent(input$update_classification_preview, {
    classification_loading(TRUE)
    req(input$classification_column, input$classification_title)
    
    # v53: cat(file=stderr(), "\n🎨 UPDATE PREVIEW (without saving)\n")
    
    # v52: Use the SAME data source as the UI (filtered_csv if available)
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
      values$filtered_csv
    } else {
      values$csv_data
    }
    
    # Get unique values and their colors
    unique_values <- unique(csv_to_use[[input$classification_column]])
    unique_values <- unique_values[!is.na(unique_values)]
    
    # v51: Add debug output for classification color mapping
    # v53: cat(file=stderr(), "[DEBUG] Classification column:", input$classification_column, "\n")
    # v53: cat(file=stderr(), "[DEBUG] Using:", if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) "filtered_csv" else "csv_data", "\n")
    # v53: cat(file=stderr(), "[DEBUG] Unique values order:", paste(unique_values, collapse=" -> "), "\n")
    
    classes <- lapply(seq_along(unique_values), function(i) {
      value <- unique_values[i]
      color_input_id <- paste0("class_color_", i)
      color_value <- NULL
      
      if (!is.null(input[[color_input_id]])) {
        color_value <- input[[color_input_id]]
      } else {
        color_value <- rainbow(length(unique_values))[i]
      }
      
      # v51: Debug output for each class-color mapping
      # v53: cat(file=stderr(), "[DEBUG] Class", i, ":", value, "-> color:", color_value, "\n")
      
      list(
        column = input$classification_column,
        value = value,
        display_name = as.character(value),
        color = color_value
      )
    })
    
    # Create temporary classification for preview
    temp_classification <- list(
      title = input$classification_title,
      column = input$classification_column,
      classes = classes,
      fdr = input$fdr_perc,
      no_cluster_title = "No cluster",
      no_cluster_color = input$no_cluster_color,
      highlight = list(enabled = FALSE)
    )
    
    # Store as temporary (will be used by generate_plot if classifications list is empty)
    values$temp_classification_preview <- temp_classification
    
    classification_loading(FALSE)
    
    # Regenerate plot
    generate_plot()
    
    showNotification("Preview updated (not saved)", type = "message", duration = 3)
  })
  
  # Add classification
  observeEvent(input$add_classification, {
    req(input$classification_column, input$classification_title)
    
    # v52: Use the SAME data source as the UI (filtered_csv if available)
    csv_to_use <- if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
      values$filtered_csv
    } else {
      values$csv_data
    }
    
    # Get unique values and their colors
    unique_values <- unique(csv_to_use[[input$classification_column]])
    unique_values <- unique_values[!is.na(unique_values)]
    
    classes <- lapply(seq_along(unique_values), function(i) {
      value <- unique_values[i]
      color_input_id <- paste0("class_color_", i)
      color_value <- NULL
      
      if (!is.null(input[[color_input_id]])) {
        color_value <- input[[color_input_id]]
      } else {
        color_value <- rainbow(length(unique_values))[i]
      }
      
      list(
        column = input$classification_column,
        value = value,
        display_name = as.character(value),
        color = color_value
      )
    })
    
    # Create new classification
    new_classification <- list(
      title = input$classification_title,
      column = input$classification_column,
      classes = classes,
      fdr = input$fdr_perc,
      no_cluster_title = "No cluster",
      no_cluster_color = input$no_cluster_color,
      highlight = list(enabled = FALSE)
    )
    
    # Add to classifications list
    values$classifications <- c(values$classifications, list(new_classification))
    
    # Set as active classification
    values$active_classification_index <- length(values$classifications)
    
    # Clear temp preview
    values$temp_classification_preview <- NULL
    
    # Update plot
    generate_plot()
    
    showNotification("Classification saved", type = "message")
  })
  
  # Observer for classification radio button selection
  observeEvent(input$selected_classification_index, {
    req(input$selected_classification_index)
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Ëœ Classification selection changed to:", input$selected_classification_index, "\n")
    
    # Update active classification index
    values$active_classification_index <- as.numeric(input$selected_classification_index)
    
    # Regenerate plot with selected classification
    generate_plot()
  })
  
  # Remove classification
  observeEvent(input$remove_classification, {
    req(values$classifications, length(values$classifications) > 0)
    req(values$active_classification_index)
    
    # Remove the SELECTED classification (not just the last one)
    index_to_remove <- values$active_classification_index
    
    if (index_to_remove > 0 && index_to_remove <= length(values$classifications)) {
      values$classifications <- values$classifications[-index_to_remove]
      
      # Adjust active index
      if (length(values$classifications) > 0) {
        values$active_classification_index <- min(index_to_remove, length(values$classifications))
      } else {
        values$active_classification_index <- NULL
      }
      
      # Update plot
      generate_plot()
      
      showNotification("Selected classification removed", type = "warning")
    }
  })
  
  # Add highlight
  observeEvent(input$add_highlight, {
    cat(file=stderr(), paste0("\n[DEBUG-2ND-HIGHLIGHT] *** ADD_HIGHLIGHT OBSERVER TRIGGERED at ", format(Sys.time(), "%H:%M:%OS3"), " ***\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   input$enable_highlight=", input$enable_highlight, "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   input$highlight_column=", input$highlight_column, "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   input$highlight_values=", paste(input$highlight_values, collapse=", "), "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   num_classifications=", length(values$classifications), "\n"))

    req(input$enable_highlight, input$highlight_column, input$highlight_values,
        length(values$classifications) > 0)

    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   req() passed, proceeding...\n"))

    # Get selected values and their colors
    selected_values <- input$highlight_values
    
    highlight_items <- lapply(seq_along(selected_values), function(i) {
      value <- selected_values[i]
      value_index <- match(value, unique(values$csv_data[[input$highlight_column]]))
      
      color_input_id <- paste0("highlight_color_", value_index)
      color_value <- NULL
      
      if (!is.null(input[[color_input_id]])) {
        color_value <- input[[color_input_id]]
      } else {
        color_value <- rainbow(length(unique(values$csv_data[[input$highlight_column]])), alpha = 0.5)[value_index]
      }
      
      list(
        column = input$highlight_column,
        value = value,
        display_name = as.character(value),
        color = color_value
      )
    })
    
    # Update last classification with highlight info
    last_classification_index <- length(values$classifications)
    
    values$classifications[[last_classification_index]]$highlight <- list(
      enabled = TRUE,
      title = input$highlight_title,
      transparency = input$highlight_transparency,
      width = input$highlight_width,
      height = input$highlight_height,
      items = highlight_items
    )
    
    # Update plot
    generate_plot()
    
    showNotification("Highlight added to current classification", type = "message")
  })
  
  # Remove highlight
  observeEvent(input$remove_highlight, {
    req(values$classifications, length(values$classifications) > 0)
    
    # Disable highlighting for last classification
    last_classification_index <- length(values$classifications)
    
    if (length(values$classifications) > 0 && 
        !is.null(values$classifications[[last_classification_index]]$highlight)) {
      values$classifications[[last_classification_index]]$highlight$enabled <- FALSE
    }
    
    # Update plot
    generate_plot()
    
    showNotification("Highlight removed from current classification", type = "message")
  })
  
  # ============================================================================
  # OBSERVERS FOR TREE DISPLAY SETTINGS
  # ============================================================================
  # Manual "Apply Settings" button removed - settings auto-apply now
  # observeEvent(input$apply_tree_settings, {
  #   req(values$plot_ready)
  #   generate_plot()
  #   showNotification("Settings applied! Plot updated.", type = "message", duration = 3)
  # })
  
  # These trigger plot regeneration when user changes display settings

  # Tip font size
  # S1.4-PERF: Using debounced input + request_plot_update() for batched updates
  observeEvent(tip_font_size_d(), {
    req(values$plot_ready)  # Only if plot has been generated at least once
    request_plot_update()
  }, ignoreInit = TRUE)

  # Edge width
  # S1.4-PERF: Using debounced input + request_plot_update() for batched updates
  observeEvent(edge_width_d(), {
    req(values$plot_ready)
    request_plot_update()
  }, ignoreInit = TRUE)

  # Tip length
  # S1.4-PERF: Using debounced input + request_plot_update() for batched updates
  observeEvent(tip_length_d(), {
    req(values$plot_ready)
    request_plot_update()
  }, ignoreInit = TRUE)
  
  # Trim tips checkbox
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$trim_tips, {
    req(values$plot_ready)
    request_plot_update()
  }, ignoreInit = TRUE)

  # Display node numbers
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$display_node_numbers, {
    req(values$plot_ready)
    request_plot_update()
  }, ignoreInit = TRUE)

  # Ladderize
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$ladderize, {
    req(values$plot_ready)
    request_plot_update()
  }, ignoreInit = TRUE)
  
  validate_rotation <- function(num_groups, rotation_prefix) {
    errors <- c()
    for (i in 1:num_groups) {
      col <- input[[paste0(rotation_prefix, "_col_", i)]]
      val <- input[[paste0(rotation_prefix, "_val_", i)]]
      if (is.null(col) || col == "" || is.null(val) || val == "") {
        errors <- c(errors, paste0("Group ", i, " incomplete"))
      }
    }
    pairs <- list()
    for (i in 1:num_groups) {
      col <- input[[paste0(rotation_prefix, "_col_", i)]]
      val <- input[[paste0(rotation_prefix, "_val_", i)]]
      if (!is.null(col) && col != "" && !is.null(val) && val != "") {
        p <- paste0(col, "=", val)
        if (p %in% pairs) {
          errors <- c(errors, paste0("Duplicate: ", p))
        }
        pairs <- c(pairs, p)
      }
    }
    return(errors)
  }
  
  observeEvent(input$apply_rotation1, {
    req(values$plot_ready, input$rotation1_num_groups)
    errors <- validate_rotation(input$rotation1_num_groups, "rotation1")
    if (length(errors) > 0) {
      showNotification(paste(errors, collapse = "\n"), type = "error", duration = 10)
      return()
    }
    
    # Save rotation1 configuration
    config <- list()
    for (i in 1:input$rotation1_num_groups) {
      col <- input[[paste0("rotation1_col_", i)]]
      val <- input[[paste0("rotation1_val_", i)]]
      if (!is.null(col) && col != "" && !is.null(val) && val != "") {
        config[[i]] <- list(col = col, val = val)
      }
    }
    values$rotation1_config <- config
    
    generate_plot()
    showNotification("Primary rotation applied!", type = "message")
  })
  
  observeEvent(input$apply_rotation2, {
    req(values$plot_ready, input$rotation2_num_groups)
    errors <- validate_rotation(input$rotation2_num_groups, "rotation2")
    if (length(errors) > 0) {
      showNotification(paste(errors, collapse = "\n"), type = "error", duration = 10)
      return()
    }
    
    # Save rotation2 configuration
    config <- list()
    for (i in 1:input$rotation2_num_groups) {
      col <- input[[paste0("rotation2_col_", i)]]
      val <- input[[paste0("rotation2_val_", i)]]
      if (!is.null(col) && col != "" && !is.null(val) && val != "") {
        config[[i]] <- list(col = col, val = val)
      }
    }
    values$rotation2_config <- config
    
    generate_plot()
    showNotification("Secondary rotation applied!", type = "message")
  })
  
  # Clear rotation configuration buttons
  observeEvent(input$clear_rotation1, {
    values$rotation1_config <- list()
    updateNumericInput(session, "rotation1_num_groups", value = 2)
    showNotification("Primary rotation configuration cleared!", type = "warning")
  })
  
  observeEvent(input$clear_rotation2, {
    values$rotation2_config <- list()
    updateNumericInput(session, "rotation2_num_groups", value = 2)
    showNotification("Secondary rotation configuration cleared!", type = "warning")
  })
  
  # Manual rotation apply button
  observeEvent(input$apply_manual_rotation, {
    cat(file=stderr(), "\n[DEBUG-ROTATION] apply_manual_rotation button clicked\n")
    cat(file=stderr(), paste0("[DEBUG-ROTATION] values$plot_ready = ", values$plot_ready, "\n"))
    cat(file=stderr(), paste0("[DEBUG-ROTATION] input$nodes_to_rotate = ", paste(input$nodes_to_rotate, collapse=", "), "\n"))
    req(values$plot_ready, input$nodes_to_rotate)
    
    if (is.null(input$nodes_to_rotate) || length(input$nodes_to_rotate) == 0) {
      showNotification("Please select at least one node to rotate", type = "error", duration = 5)
      return()
    }
    
    # Save manual rotation configuration
    values$manual_rotation_config <- as.numeric(input$nodes_to_rotate)
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â§ MANUAL ROTATION APPLIED\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â§ Nodes to rotate:", paste(values$manual_rotation_config, collapse=", "), "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â§ Regenerating plot...\n")
    
    generate_plot()
    showNotification("Manual rotation applied!", type = "message")
  })
  
  # Clear manual rotation button
  observeEvent(input$clear_manual_rotation, {
    values$manual_rotation_config <- list()
    updateSelectizeInput(session, "nodes_to_rotate", selected = character(0))
    showNotification("Manual rotation configuration cleared!", type = "warning")
  })
  
  # === NEW: Reactive observer for node number font size ===
  # S1.4-PERF: Using debounced input + request_plot_update() for batched updates
  observeEvent(node_number_font_size_d(), {
    req(values$plot_ready)
    req(input$display_node_numbers == TRUE)  # Only regenerate if node numbers are displayed
    request_plot_update()
  }, ignoreInit = TRUE)
  
  # S1.4-PERF: Removed duplicate display_node_numbers observer (already handled earlier)
  
  # === Reactive observer for highlight checkbox ===
  # S1.4-PERF: Using request_plot_update() for batched updates
  observeEvent(input$highlight_selected_nodes, {
    req(values$plot_ready)
    req(input$nodes_to_rotate)  # Only if nodes are selected
    request_plot_update()
  }, ignoreInit = TRUE)
  
  # ============================================================================
  
  # Make sure font_size section has all required fields
  update_font_sizes <- function() {
    # Ensure these values exist and are properly set
    # v128: Use Legend tab font size settings for legend_title and legend_text
    # This ensures bootstrap legend matches other legends' font sizes
    values$yaml_data$`visual definitions`$font_size <- list(
      tips = if(!is.null(input$tip_font_size)) input$tip_font_size else 3,
      legend_title = if(!is.null(input$legend_title_size)) input$legend_title_size else 12,
      legend_text = if(!is.null(input$legend_text_size)) input$legend_text_size else 10,
      legend_box = 15,    # Fixed value
      heat_map_title = 25, # Fixed value
      heat_map_legend = if(!is.null(input$heatmap_font_size)) input$heatmap_font_size else 3.8
    )
  }

  # Generate plot based on current settings
  generate_plot <- function() {
    # S1.4-PERF: Simplified entry logging (reduced from verbose DEBUG-2ND-HIGHLIGHT)
    cat(file=stderr(), "[PERF] generate_plot() called\n")

    # === GUARD CHECKS (S1.4-PERF: consolidated - removed duplicate guards) ===

    # S2.0-PERF: Quick cooldown to prevent rapid consecutive calls (500ms)
    # This catches cascading calls from multiple observers
    time_since_last <- as.numeric(Sys.time()) * 1000 - last_plot_time()
    if (time_since_last < 500 && time_since_last > 0) {  # >0 to allow first call
      cat(file=stderr(), sprintf("[PERF] Skipping - rapid call (%.0fms since last plot)\n", time_since_last))
      return(NULL)
    }

    if (classification_loading()) {
      cat(file=stderr(), "[PERF] Skipping - classification UI loading\n")
      return(NULL)
    }

    # Dont generate if already generating (recursion guard)
    if (!is.null(values$plot_generating) && values$plot_generating) {
      cat(file=stderr(), "[PERF] Skipping - plot already generating\n")
      return(NULL)
    }
    # === END GUARD CHECKS ===

    # v57: Show processing status IMMEDIATELY via shinyjs (before R blocks)
    show_status_processing()

    # Set generating flag for status indicator
    values$plot_generating <- TRUE
    values$plot_ready <- FALSE  # Mark as not ready while generating
    #shiny::invalidateLater(0, session)

    # S1-PERF: Removed Sys.sleep(0.3) - unnecessary delay

    values$progress_message <- "🎨 Generating your beautiful plot..."
    values$progress_visible <- TRUE

    # S1-PERF: Removed Sys.sleep(0.2) - unnecessary delay

    # Check if we have the necessary data
    if (is.null(values$tree) || is.null(values$csv_data)) {
      # v53: cat(file=stderr(), "Missing required data (tree or csv_data)\n")
      values$progress_visible <- FALSE
      values$plot_generating <- FALSE
      showNotification("Missing required data for plot generation", type = "warning")
      return(NULL)
    }
    
    # Check if temp CSV exists (created during Process Data)
    if (is.null(values$temp_csv_path) || !file.exists(values$temp_csv_path)) {
      # v53: cat(file=stderr(), "Missing temp CSV file\n")
      values$progress_visible <- FALSE
      values$plot_generating <- FALSE
      showNotification("Please process data first", type = "warning")
      return(NULL)
    }
    
    # v53: cat(file=stderr(), "All checks passed, proceeding with plot generation\n")
    
    # Check if we have any matches
    if (values$id_match$summary$unmatched == values$id_match$summary$total_tree_labels) {
      values$progress_visible <- FALSE
      showNotification("No matching between tree tips and CSV data", type = "warning")
      return(NULL)
    }
    
    # v53: cat(file=stderr(), 'values$id_match$mapping=', "\n")
    if (length(values$id_match$mapping) > 0) {
      for (i in 1:min(5, length(values$id_match$mapping))) {
        # v53: cat(file=stderr(), names(values$id_match$mapping)[i], "->", 
        #     values$id_match$mapping[[i]], "\n")
      }
    } else {
      # v53: cat(file=stderr(), "No mappings found\n")
    }
    
    # S1-PERF: Filter CSV to only include matched rows - use unlist instead of loop
    matched_ids <- unique(unlist(values$id_match$mapping, use.names = FALSE))
    
    values$filtered_csv <- values$csv_data[values$csv_data[[input$id_column]] %in% matched_ids, ]
    # v53: print("Filtered CSV rows:")
    # v53: print(nrow(values$filtered_csv))
    
    # Check if filtered CSV has data
    if (nrow(values$filtered_csv) == 0) {
      showNotification("After filtering, no CSV rows matched tree tips", type = "warning")
      return(NULL)
    }
    
    # DEBUG: Check trimming_params before update_yaml
    # v53: cat(file=stderr(), "\nðŸ”§ === BEFORE update_yaml() ===\n")
    # v53: cat(file=stderr(), "ðŸ”§ values$trimming_params is NULL:", is.null(values$trimming_params), "\n")
    if (!is.null(values$trimming_params)) {
      # v53: cat(file=stderr(), "ðŸ”§ trimming_params$id_tip_trim_flag:", values$trimming_params$id_tip_trim_flag, "\n")
      # v53: cat(file=stderr(), "ðŸ”§ trimming_params$id_tip_trim_start:", values$trimming_params$id_tip_trim_start, "\n")
      # v53: cat(file=stderr(), "ðŸ”§ trimming_params$id_tip_trim_end:", values$trimming_params$id_tip_trim_end, "\n")
    }
    # v53: cat(file=stderr(), "ðŸ”§ ==============================\n\n")
    # DEBUG: Confirm trimming params were restored
    # v53: cat(file=stderr(), "\nðŸ” === AFTER RESTORING TRIMMING PARAMS ===\n")
    # v53: cat(file=stderr(), "ðŸ” values$yaml_data trimming values:\n")
    # v53: cat(file=stderr(), "ðŸ”   id_tip_trim_flag:", values$yaml_data$`visual definitions`$id_tip_trim_flag, "\n")
    # v53: cat(file=stderr(), "ðŸ”   id_tip_trim_start:", values$yaml_data$`visual definitions`$id_tip_trim_start, "\n")
    # v53: cat(file=stderr(), "ðŸ”   id_tip_trim_end:", values$yaml_data$`visual definitions`$id_tip_trim_end, "\n")
    # v53: cat(file=stderr(), "ðŸ” ========================================\n\n")
    
    
    # Update YAML with current UI settings
    update_yaml()
    
    # If we have inferred trimming parameters, make sure they're in the YAML
    if (!is.null(values$trimming_params)) {
      values$yaml_data$`visual definitions`$id_tip_trim_flag <- 
        if (values$trimming_params$id_tip_trim_flag) "yes" else "no"
      values$yaml_data$`visual definitions`$id_tip_trim_start <- 
        values$trimming_params$id_tip_trim_start
      values$yaml_data$`visual definitions`$id_tip_trim_end <- 
        values$trimming_params$id_tip_trim_end
      values$yaml_data$`visual definitions`$id_tip_prefix <- 
        values$trimming_params$id_tip_prefix
    }
    
    # Make sure font sizes are properly set
    update_font_sizes()
    
    # Ensure all required flags exist
    if (is.null(values$yaml_data$`visual definitions`$flag_make_newick_file)) {
      values$yaml_data$`visual definitions`$flag_make_newick_file <- "no"
    }
    if (is.null(values$yaml_data$`visual definitions`$path_out_newick)) {
      values$yaml_data$`visual definitions`$path_out_newick <- ""
    }
    if (is.null(values$yaml_data$`visual definitions`$tip_name_display_flag)) {
      values$yaml_data$`visual definitions`$tip_name_display_flag <- "yes"
    }
    if (is.null(values$yaml_data$`visual definitions`$simulate_p_value)) {
      values$yaml_data$`visual definitions`$simulate_p_value <- "yes"
    }
    if (is.null(values$yaml_data$`visual definitions`$flag_calc_scores_for_tree)) {
      values$yaml_data$`visual definitions`$flag_calc_scores_for_tree <- "no"
    }
    if (is.null(values$yaml_data$`visual definitions`$flag_display_nod_number_on_tree)) {
      values$yaml_data$`visual definitions`$flag_display_nod_number_on_tree <- "no"
    }
    
    # CRITICAL: Add compare_two_trees flag
    if (is.null(values$yaml_data$`visual definitions`$compare_two_trees)) {
      values$yaml_data$`visual definitions`$compare_two_trees <- "no"
    }
    
    # Add other potentially missing flags
    if (is.null(values$yaml_data$`visual definitions`$debug_mode)) {
      values$yaml_data$`visual definitions`$debug_mode <- "no"
    }
    if (is.null(values$yaml_data$`visual definitions`$debug_print_data_tree)) {
      values$yaml_data$`visual definitions`$debug_print_data_tree <- "no"
    }
    
    # Create temporary YAML file for func.print.lineage.tree
    temp_yaml_file <- tempfile(fileext = ".yaml")
    
    # Use the temp CSV path created during "Process Data"
    # Update YAML to point to this temp CSV
    yaml_data_modified <- values$yaml_data
    yaml_data_modified$`Individual general definitions`$`mapping csv file` <- values$temp_csv_path

    # S1.62dev: ALWAYS use PDF for internal file operations
    # The user's selected format (JPEG/PNG/etc.) is only used for final download.
    # JPEG/PNG are much slower to render than PDF for complex plots, causing crashes
    # when combined with portrait orientation and other intensive operations.
    yaml_data_modified$`Individual general definitions`$out_file$file_type <- "pdf"
    
    # DEBUG: Print classification structure
    if (!is.null(values$yaml_data$`visual definitions`$classification)) {
      # v53: cat(file=stderr(), "\n=== CLASSIFICATION YAML STRUCTURE ===\n")
      # v53: cat(file=stderr(), "Number of classifications:", length(values$yaml_data$`visual definitions`$classification), "\n")
      for (i in seq_along(values$yaml_data$`visual definitions`$classification)) {
        class_item <- values$yaml_data$`visual definitions`$classification[[i]]
        # v53: debug_cat(paste0("Classification ", i, ":\n"))
        # v53: cat(file=stderr(), "  Title:", class_item[[as.character(i)]]$title, "\n")
        if (!is.null(class_item[[as.character(i)]]$according)) {
          # v53: cat(file=stderr(), "  Number of classes:", length(class_item[[as.character(i)]]$according), "\n")
          # Print first class as example
          if (length(class_item[[as.character(i)]]$according) > 0) {
            first_class <- class_item[[as.character(i)]]$according[[1]]
            # v53: cat(file=stderr(), "  First class title1:", first_class[["1"]]$title1, "\n")
            # v53: cat(file=stderr(), "  First class value1:", paste(first_class[["1"]]$value1, collapse=", "), "\n")
            # v53: cat(file=stderr(), "  First class color:", first_class[["1"]]$color, "\n")
          }
        }
      }
      # v53: debug_cat("=====================================\n\n")
    }
    
    # DEBUG: Show trimming params before writing YAML
    # v53: cat(file=stderr(), "\nðŸ“ === BEFORE WRITING YAML FILE ===\n")
    # v53: cat(file=stderr(), "ðŸ“ yaml_data_modified trimming values:\n")
    # v53: cat(file=stderr(), "ðŸ“   id_tip_trim_flag:", yaml_data_modified$`visual definitions`$id_tip_trim_flag, "\n")
    # v53: cat(file=stderr(), "ðŸ“   id_tip_trim_start:", yaml_data_modified$`visual definitions`$id_tip_trim_start, "\n")
    # v53: cat(file=stderr(), "ðŸ“   id_tip_trim_end:", yaml_data_modified$`visual definitions`$id_tip_trim_end, "\n")
    # v53: cat(file=stderr(), "ðŸ“ ==================================\n\n")
    
    # Write modified YAML to temporary file
    writeLines(yaml::as.yaml(yaml_data_modified, indent.mapping.sequence = TRUE), temp_yaml_file)

    # v57: DEBUG - Show heatmap structure in YAML
    if (!is.null(yaml_data_modified$`visual definitions`$classification)) {
      debug_cat("\n=== v57 HEATMAP YAML DEBUG ===\n")
      for (ci in seq_along(yaml_data_modified$`visual definitions`$classification)) {
        class_entry <- yaml_data_modified$`visual definitions`$classification[[ci]]
        class_key <- names(class_entry)[1]
        debug_cat(paste0("Classification ", ci, " (key=", class_key, "):\n"))
        debug_cat(paste0("  Has heatmap_display: ", !is.null(class_entry[[class_key]]$heatmap_display), "\n"))
        if (!is.null(class_entry[[class_key]]$heatmap_display)) {
          debug_cat(paste0("  Number of heatmaps: ", length(class_entry[[class_key]]$heatmap_display), "\n"))
          for (hi in seq_along(class_entry[[class_key]]$heatmap_display)) {
            hm <- class_entry[[class_key]]$heatmap_display[[hi]]
            hm_key <- names(hm)[1]
            debug_cat(paste0("    Heatmap ", hi, " (key=", hm_key, "):\n"))
            debug_cat(paste0("      Title: ", hm[[hm_key]]$title, "\n"))
            debug_cat(paste0("      Display: ", hm[[hm_key]]$display, "\n"))
            # S1.62dev: Add data_source debug
            debug_cat(paste0("      S1.62dev data_source: ", if (!is.null(hm[[hm_key]]$data_source)) hm[[hm_key]]$data_source else "NULL", "\n"))
            debug_cat(paste0("      S1.62dev all fields: ", paste(names(hm[[hm_key]]), collapse=", "), "\n"))
            debug_cat(paste0("      According length: ", length(hm[[hm_key]]$according), "\n"))
            if (length(hm[[hm_key]]$according) > 0) {
              for (ai in seq_along(hm[[hm_key]]$according)) {
                acc <- hm[[hm_key]]$according[[ai]]
                acc_key <- names(acc)[1]
                debug_cat(paste0("        Column ", ai, ": ", acc[[acc_key]], "\n"))
              }
            }
          }
        }
      }
      debug_cat("==============================\n\n")
    }

    # DEBUG: Also print part of the YAML file
    # v53: cat(file=stderr(), "\n=== YAML FILE CONTENT (first 50 lines) ===\n")
    yaml_lines <- readLines(temp_yaml_file, n = 50)
    # v53: cat(file=stderr(), paste(yaml_lines, collapse = "\n"), "\n")
    # v53: debug_cat("==========================================\n\n")
    
    # Define viridis_option_list
    viridis_option_list <- c("A", "B", "C", "D", "E")
    
    # Call the tree visualization function with better error handling
    result <- tryCatch({
      # Set sensible defaults
      width_val <- if (!is.null(input$output_width) && !is.na(input$output_width)) input$output_width else 170
      height_val <- if (!is.null(input$output_height) && !is.na(input$output_height)) input$output_height else 60
      units_val <- if (!is.null(input$output_units) && !is.na(input$output_units)) input$output_units else "cm"
      
      # Add debug output
      # v53: print(paste("Calling func.print.lineage.tree with width:", width_val, 
      #             "height:", height_val, "units:", units_val))
      
      # === DEBUG CHECKPOINT 1: UI INPUT VALUES ===
      # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 1: READING UI INPUTS\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â input$node_number_font_size:", input$node_number_font_size, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â input$display_node_numbers:", input$display_node_numbers, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â input$highlight_selected_nodes:", input$highlight_selected_nodes, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â input$nodes_to_rotate:", paste(input$nodes_to_rotate, collapse=", "), "\n")
      
      # Prepare parameter values with debug output
      font_size_to_pass <- if (!is.null(input$node_number_font_size)) {
        input$node_number_font_size
      } else {
        3.5
      }
      
      display_nodes_to_pass <- if (!is.null(input$display_node_numbers)) {
        input$display_node_numbers
      } else {
        FALSE
      }
      
      highlight_flag_to_pass <- if (!is.null(input$highlight_selected_nodes) && 
                                    !is.null(input$nodes_to_rotate) && 
                                    length(input$nodes_to_rotate) > 0) {
        input$highlight_selected_nodes
      } else {
        FALSE
      }
      
      nodes_to_highlight_to_pass <- if (!is.null(input$nodes_to_rotate) && 
                                        length(input$nodes_to_rotate) > 0) {
        as.numeric(input$nodes_to_rotate)
      } else {
        NA
      }
      
      # v53: cat(file=stderr(), "Ã°Å¸â€Â font_size_to_pass:", font_size_to_pass, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â display_nodes_to_pass:", display_nodes_to_pass, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â highlight_flag_to_pass:", highlight_flag_to_pass, "\n")
      # v53: cat(file=stderr(), "Ã°Å¸â€Â nodes_to_highlight_to_pass:", paste(nodes_to_highlight_to_pass, collapse=", "), "\n")
      # v53: debug_cat("================================================\n\n")
      
      # Call func.print.lineage.tree with the temp YAML file
      # v54: Wrap in suppressWarnings to suppress -Inf and other harmless warnings
      cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CALLING func.print.lineage.tree at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
      cat(file=stderr(), paste0("[DEBUG-ROTATION] values$manual_rotation_config = ", paste(values$manual_rotation_config, collapse=", "), "\n"))
      cat(file=stderr(), paste0("[DEBUG-ROTATION] length(manual_rotation_config) = ", length(values$manual_rotation_config), "\n"))
      tree_result <- suppressWarnings(func.print.lineage.tree(
        conf_yaml_path = temp_yaml_file,
        width = width_val,
        height = height_val,
        units_out = units_val,
        debug_mode = FALSE,
        compare_two_trees = FALSE,
        list_nodes_to_rotate = if (!is.null(values$manual_rotation_config) && 
                                   length(values$manual_rotation_config) > 0) {
          values$manual_rotation_config
        } else {
          NA
        },
        flag_display_nod_number_on_tree = display_nodes_to_pass,
        node_number_font_size = font_size_to_pass,
        highlight_manual_nodes = highlight_flag_to_pass,
        manual_nodes_to_highlight = nodes_to_highlight_to_pass,
        #bootstrap_label_size=bootstrap_label_size,
        bootstrap_label_size = if (!is.null(input$bootstrap_label_size)) {
          input$bootstrap_label_size
        } else {
          3.5
        },
        # v111: Pass bootstrap position offset slider value
        man_boot_x_offset = if (!is.null(input$man_boot_x_offset)) {
          input$man_boot_x_offset
        } else {
          0
        },
        # v103: Pass heatmap tree distance slider value
        heatmap_tree_distance = if (!is.null(input$heatmap_tree_distance)) {
          input$heatmap_tree_distance
        } else {
          0.02
        },
        # v125: Pass global gap between heatmaps
        heatmap_global_gap = if (!is.null(values$heatmap_global_gap)) {
          values$heatmap_global_gap
        } else {
          0.05
        },
        # v135: Pass legend settings for highlight/bootstrap legends
        legend_settings = values$legend_settings,
        # S1.62dev: Pass RData CNV matrix for heatmaps with data_source="rdata"
        rdata_cnv_matrix = values$rdata_cnv_matrix,
        # S2.0-PERF: Two-tier caching - pass cached p-values (Option 3A)
        cached_p_list_of_pairs = values$cached_p_list_of_pairs,
        cached_p_list_hash = values$cached_p_list_hash,
        # S2.7-PERF: Multi-entry cache for instant switching between classifications
        p_list_cache = values$p_list_cache,
        # S2.9-PERF: Heatmap cache for reusing unchanged heatmaps
        heatmap_cache = values$heatmap_cache
      ))
      cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] RETURNED from func.print.lineage.tree at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))

      # S2.0-PERF: Extract plots and cache data from new return structure (Option 3A)
      # tree_result now has structure: list(plots = out_trees, cache_data = cache_data)
      plots_result <- NULL
      if (!is.null(tree_result) && is.list(tree_result)) {
        if (!is.null(tree_result$plots)) {
          plots_result <- tree_result$plots
          cat(file=stderr(), "[PERF-CACHE] Extracted plots from tree_result$plots\n")
        } else {
          # Fallback for backward compatibility if structure is different
          plots_result <- tree_result
          cat(file=stderr(), "[PERF-CACHE] Using tree_result directly (legacy format)\n")
        }

        # Store cache data for future use
        if (!is.null(tree_result$cache_data)) {
          # S2.0-PERF: Legacy single-entry cache (for backward compatibility)
          values$cached_p_list_of_pairs <- tree_result$cache_data$p_list_of_pairs
          values$cached_p_list_hash <- tree_result$cache_data$p_list_hash
          values$cached_classification_column <- input$classification_column

          # S2.7-PERF: Multi-entry cache - store in named list keyed by hash
          # This allows instant switching between previously computed classifications
          new_hash <- tree_result$cache_data$p_list_hash
          if (!is.null(new_hash) && !is.null(tree_result$cache_data$p_list_of_pairs)) {
            # Add new entry to cache (or update existing)
            values$p_list_cache[[new_hash]] <- tree_result$cache_data$p_list_of_pairs

            # Update LRU order: remove if exists, then add to end (most recent)
            values$p_list_cache_order <- c(
              setdiff(values$p_list_cache_order, new_hash),
              new_hash
            )

            # Evict oldest entries if cache exceeds max size
            max_size <- values$p_list_cache_max_size
            if (is.null(max_size)) max_size <- 10
            while (length(values$p_list_cache_order) > max_size) {
              oldest_hash <- values$p_list_cache_order[1]
              values$p_list_cache[[oldest_hash]] <- NULL
              values$p_list_cache_order <- values$p_list_cache_order[-1]
              cat(file=stderr(), sprintf("[PERF-CACHE] Evicted oldest cache entry (hash: %s)\n",
                                         substr(oldest_hash, 1, 8)))
            }

            cat(file=stderr(), sprintf("[PERF-CACHE] Stored in multi-cache (hash: %s, %d/%d entries, classification: %s)\n",
                                       substr(new_hash, 1, 8),
                                       length(values$p_list_cache_order), max_size,
                                       values$cached_classification_column))
          } else {
            cat(file=stderr(), sprintf("[PERF-CACHE] Stored cache (hash: %s, classification: %s)\n",
                                       substr(values$cached_p_list_hash, 1, 8),
                                       values$cached_classification_column))
          }

          # S2.9-PERF: Store updated heatmap cache
          if (!is.null(tree_result$cache_data$heatmap_cache)) {
            values$heatmap_cache <- tree_result$cache_data$heatmap_cache
            cat(file=stderr(), sprintf("[S2.9-CACHE] Updated heatmap cache (%d entries)\n",
                                       length(values$heatmap_cache)))
          }
        }
      }

      # Debug output
      # v53: cat(file=stderr(), "\n=== AFTER func.print.lineage.tree ===\n")
      # v53: cat(file=stderr(), "tree_result is NULL:", is.null(tree_result), "\n")
      if (!is.null(plots_result) && is.list(plots_result)) {
        # v53: cat(file=stderr(), "plots_result is a list with", length(plots_result), "element(s)\n")
        # v53: cat(file=stderr(), "List names:", paste(names(plots_result), collapse=", "), "\n")
        if (length(plots_result) > 0) {
          # v53: cat(file=stderr(), "First element class:", class(plots_result[[1]]), "\n")
          # v53: cat(file=stderr(), "First element inherits ggplot:", inherits(plots_result[[1]], "ggplot"), "\n")
          if (length(plots_result) > 1) {
            # v53: cat(file=stderr(), "NOTE: Multiple plots returned (", length(plots_result), "). Using the last one.\n")
          }
        }
      }
      # v53: debug_cat("====================================\n\n")

      # Extract the plot from out_trees list
      # The function returns out_trees which is a list indexed by numbers like "1", "2", etc.
      # When multiple classifications exist, use the LAST plot (most complete)
      tree_plot <- NULL
      if (!is.null(plots_result) && is.list(plots_result) && length(plots_result) > 0) {
        # Extract the LAST plot from the list (most recent classification)
        plot_index <- length(plots_result)
        tree_plot <- plots_result[[plot_index]]
        # v53: cat(file=stderr(), "Successfully extracted plot from plots_result[[", plot_index, "]]\n")
        # v53: cat(file=stderr(), "Plot class:", class(tree_plot), "\n")
        # v53: cat(file=stderr(), "Plot inherits ggplot:", inherits(tree_plot, "ggplot"), "\n")
      } else {
        # v53: cat(file=stderr(), "ERROR: Could not extract plot from plots_result\n")
        # v53: cat(file=stderr(), "plots_result structure:\n")
        # v54: str(plots_result)
      }
      
      # Check if we got a valid plot object
      if (is.null(tree_plot) || !inherits(tree_plot, "ggplot")) {
        # v53: cat(file=stderr(), "WARNING: Could not extract valid ggplot object from result\n")
        return(NULL)
      }
      
      # v53: cat(file=stderr(), "Returning valid plot object\n")
      
      # IMPORTANT: Build the plot to resolve all data references
      # v53: cat(file=stderr(), "Building plot to resolve data references...\n")
      built_plot <- tryCatch({
        ggplot2::ggplot_build(tree_plot)
        tree_plot  # Return the original plot if build succeeds
      }, error = function(e) {
        # v53: cat(file=stderr(), "Error building plot:", e$message, "\n")
        # Try to just return the plot anyway
        tree_plot
      })
      
      # Don't return yet - we need to save the plot first!
      # v53: cat(file=stderr(), "Plot object prepared\n")
      built_plot  # This will be the return value of tryCatch
      
    }, error = function(e) {
      # S1.62dev: Print full error details to console for debugging
      cat(file=stderr(), "\n=== ERROR in generate_plot tryCatch ===\n")
      cat(file=stderr(), paste0("Error message: ", e$message, "\n"))
      cat(file=stderr(), paste0("Error call: ", deparse(e$call), "\n"))

      # IMPORTANT: Hide progress bar when error occurs
      values$progress_visible <- FALSE
      values$progress_message <- ""
      values$plot_generating <- FALSE

      # Show user-friendly notification
      showNotification(
        paste("Plot generation failed:", e$message),
        type = "error",
        duration = 10
      )
      
      # Clean up temporary YAML file
      if (exists("temp_yaml_file")) {
        unlink(temp_yaml_file)
      }
      
      return(NULL)
    })
    
    # Clean up temporary YAML file
    unlink(temp_yaml_file)
    
    # DEBUG: Check what we got back
    # v53: cat(file=stderr(), "\n=== Plot Generation Result ===\n")
    # v53: cat(file=stderr(), "result is NULL:", is.null(result), "\n")
    if (!is.null(result)) {
      # v53: cat(file=stderr(), "result class:", class(result), "\n")
      # v53: cat(file=stderr(), "result inherits ggplot:", inherits(result, "ggplot"), "\n")
    }
    # v53: debug_cat("==============================\n\n")
    
    # If we got a valid result
    # If we got a valid result
    # If we got a valid result
    if (!is.null(result)) {
      # v53: cat(file=stderr(), "=== Attempting to save plot ===\n")

      # v121: Apply legend settings from the Legend tab
      # v180: Enhanced with key dimensions, byrow, box background, margin controls
      legend_settings <- values$legend_settings
      if (!is.null(legend_settings)) {
        debug_cat(paste0("\n=== v180: Applying legend settings to plot ===\n"))
        debug_cat(paste0("  Position: ", legend_settings$position, "\n"))

        # v125: Determine legend layout based on position
        # For top/bottom: legends arranged horizontally, but each legend has title above values
        # For left/right: legends arranged vertically, with title next to values
        is_horizontal_position <- legend_settings$position %in% c("top", "bottom")

        # v180: Get spacing values (use defaults if not set)
        h_spacing <- if (!is.null(legend_settings$spacing)) legend_settings$spacing else 0.3
        v_spacing <- if (!is.null(legend_settings$spacing_vertical)) legend_settings$spacing_vertical else 1
        title_key_spacing <- if (!is.null(legend_settings$title_key_spacing)) legend_settings$title_key_spacing else 0.2
        key_spacing <- if (!is.null(legend_settings$key_spacing)) legend_settings$key_spacing else 0.2

        # v180: Get key dimensions (use defaults if not set)
        key_width <- if (!is.null(legend_settings$key_width)) legend_settings$key_width else 1
        key_height <- if (!is.null(legend_settings$key_height)) legend_settings$key_height else 1

        # S2.292dev: Get font family setting
        font_family <- if (!is.null(legend_settings$font_family)) legend_settings$font_family else "sans"
        cat(file=stderr(), sprintf("[FONT-DEBUG] legend_settings$font_family = '%s'\n",
                                   ifelse(is.null(legend_settings$font_family), "NULL", legend_settings$font_family)))
        cat(file=stderr(), sprintf("[FONT-DEBUG] Using font_family = '%s'\n", font_family))

        # v180: Get background settings
        # S1.5: Fix RGBA colors from colourpicker - extract RGB portion if 8-char hex
        # colourpicker with allowTransparent=TRUE returns #RRGGBBAA format
        # where AA is alpha (00=transparent, FF=opaque). We convert to #RRGGBB for ggplot2.
        # S1.62dev: Fixed - check if alpha is 00 (fully transparent) and use "transparent" string
        box_bg_raw <- if (!is.null(legend_settings$box_background)) legend_settings$box_background else "transparent"
        if (!is.null(box_bg_raw) && nchar(box_bg_raw) == 9 && substr(box_bg_raw, 1, 1) == "#") {
          # 8-char hex with # prefix = #RRGGBBAA
          alpha_hex <- substr(box_bg_raw, 8, 9)
          if (alpha_hex == "00") {
            # S1.62dev: Fully transparent - use "transparent" string instead of stripping to black
            box_bg <- "transparent"
          } else {
            # Has some opacity - extract just #RRGGBB
            box_bg <- substr(box_bg_raw, 1, 7)
          }
        } else {
          box_bg <- box_bg_raw
        }
        legend_margin_val <- if (!is.null(legend_settings$margin)) legend_settings$margin else 0.2

        # v179: Use horizontal spacing for top/bottom, vertical spacing for left/right
        spacing_val <- if (is_horizontal_position) h_spacing else v_spacing

        # Build theme modifications for legend
        # v180: Added key width/height, title-key spacing, box background, margin
        # S1.5: Fixed - added legend.background and legend.key for proper background coloring
        # S2.292dev: Added font_family to legend text elements
        # S2.292dev-FIX: Also set text element as parent to ensure font propagates
        cat(file=stderr(), sprintf("[FONT-THEME] Setting legend.title and legend.text with family='%s'\n", font_family))
        legend_theme <- theme(
          # Set base text family for all text elements (parent)
          text = element_text(family = font_family),
          legend.position = legend_settings$position,
          legend.title = element_text(size = legend_settings$title_size, face = "bold", family = font_family),
          legend.text = element_text(size = legend_settings$text_size, family = font_family),
          legend.key.size = unit(legend_settings$key_size, "lines"),
          legend.key.width = unit(key_width, "lines"),    # v180: Custom key width
          legend.key.height = unit(key_height, "lines"),  # v180: Custom key height
          legend.key = element_rect(fill = box_bg, colour = NA),  # S1.5: Key backgrounds match legend bg
          legend.spacing = unit(spacing_val, "cm"),
          legend.spacing.x = unit(h_spacing, "cm"),
          legend.spacing.y = unit(v_spacing, "cm"),
          legend.key.spacing = unit(key_spacing, "cm"),           # v180: Between keys spacing
          legend.key.spacing.y = unit(key_spacing, "cm"),         # v180: Vertical between keys
          legend.title.position = "top",                           # v180: Title above keys
          legend.background = element_rect(fill = box_bg, colour = NA),      # S1.5: Individual legend backgrounds
          legend.box.background = element_rect(fill = box_bg, colour = NA),  # v180: Outer box background
          legend.margin = margin(legend_margin_val, legend_margin_val, legend_margin_val, legend_margin_val, "cm"),  # v180
          # v125: For top/bottom, arrange legends horizontally but stack items vertically
          legend.box = if (is_horizontal_position) "horizontal" else "vertical",
          legend.direction = if (is_horizontal_position) "vertical" else "vertical"
        )

        # Apply the legend theme
        result <- result + legend_theme

        debug_cat(paste0("  v80: Legend box=", if (is_horizontal_position) "horizontal" else "vertical",
                                   ", direction=vertical\n"))
        debug_cat(paste0("  v80: Key dims: ", key_width, "x", key_height, " lines\n"))
        debug_cat(paste0("  v80: Spacings - H:", h_spacing, ", V:", v_spacing,
                                   ", Title-key:", title_key_spacing, ", Key:", key_spacing, "\n"))
        debug_cat(paste0("  v80: Box bg:", box_bg, ", Margin:", legend_margin_val, "\n"))
        debug_cat(paste0("  S2.292dev: Font family: '", font_family, "'\n"))
        cat(file=stderr(), paste0("[LEGEND-FONT] Applying font_family='", font_family, "' to legend\n"))

        # v180: Apply visibility controls using guides()
        # Also apply reverse order if requested
        guides_list <- list()
        reverse_order <- isTRUE(legend_settings$reverse_order)

        # Hide specific legends based on visibility settings
        if (!isTRUE(legend_settings$show_classification)) {
          guides_list$colour <- "none"
        } else {
          guides_list$colour <- guide_legend(reverse = reverse_order)
        }

        # v180: P value uses size aesthetic
        if (!isTRUE(legend_settings$show_pvalue)) {
          guides_list$size <- "none"
        } else {
          guides_list$size <- guide_legend(reverse = reverse_order)
        }

        # v180: Heatmaps use fill aesthetic
        # NOTE: With ggnewscale, each heatmap has its own fill scale.
        # Using guides(fill = ...) only affects the LAST fill scale, which can break
        # earlier heatmaps. Only hide ALL fill legends if explicitly requested.
        # When show_heatmap is TRUE, we don't override - let each scale keep its own guide.
        if (!isTRUE(legend_settings$show_heatmap)) {
          # v180: To hide all fill legends with ggnewscale, we set fill to "none"
          # This will hide the last scale's legend; earlier scales need their guides
          # set to "none" when created. For now, this is a partial solution.
          guides_list$fill <- "none"
        }
        # v180: When show_heatmap is TRUE, we DON'T add fill to guides_list
        # This prevents overriding individual heatmap scale guides

        if (length(guides_list) > 0) {
          result <- result + do.call(guides, guides_list)
        }

        debug_cat(paste0("  v80: Reverse order=", reverse_order, "\n"))
        debug_cat(paste0("  Legend settings applied successfully\n"))
      }

      # v130: Apply Extra tab settings (page title, custom texts, images)
      # These are applied AFTER legend settings so they appear on top
      tryCatch({
        # v144: Plot position offsets are now applied AFTER all modifications
        # using cowplot's draw_plot for true translation (no squeezing)
        # Store offsets for later use during rendering
        plot_off_x <- if (!is.null(values$plot_offset_x)) values$plot_offset_x else 0
        plot_off_y <- if (!is.null(values$plot_offset_y)) values$plot_offset_y else 0
        # v146: Get scale percentage
        plot_scale <- if (!is.null(values$plot_scale_percent)) values$plot_scale_percent else 100
        # v179: Get tree stretch values
        tree_stretch_x <- if (!is.null(values$tree_stretch_x)) values$tree_stretch_x else 1
        tree_stretch_y <- if (!is.null(values$tree_stretch_y)) values$tree_stretch_y else 1
        # v179: Get background color
        bg_color <- if (!is.null(values$background_color)) values$background_color else "#FFFFFF"

        # v144: Store the offsets in result for later extraction
        # We'll apply them during the final rendering step
        attr(result, "plot_offset_x") <- plot_off_x
        attr(result, "plot_offset_y") <- plot_off_y
        attr(result, "plot_scale_percent") <- plot_scale
        # v179: Store tree stretch values
        attr(result, "tree_stretch_x") <- tree_stretch_x
        attr(result, "tree_stretch_y") <- tree_stretch_y
        attr(result, "background_color") <- bg_color
        # v180: Store keep_proportions setting for plot proportions preservation
        keep_proportions <- if (!is.null(input$keep_proportions)) input$keep_proportions else FALSE
        attr(result, "keep_proportions") <- keep_proportions

        if (plot_off_x != 0 || plot_off_y != 0 || plot_scale != 100) {
          debug_cat(paste0("\n=== v146: STORING PLOT POSITION & SCALE ===\n"))
          debug_cat(paste0("  X offset: ", plot_off_x, " (positive = right)\n"))
          debug_cat(paste0("  Y offset: ", plot_off_y, " (positive = up)\n"))
          debug_cat(paste0("  Scale: ", plot_scale, "%\n"))
          debug_cat(paste0("  v46: Offsets and scale will be applied via cowplot draw_plot during rendering\n"))
        }

        # v179: Apply background color
        if (bg_color != "#FFFFFF") {
          debug_cat(paste0("\n=== v179: APPLYING BACKGROUND COLOR ===\n"))
          debug_cat(paste0("  Background: ", bg_color, "\n"))
          result <- result + theme(
            plot.background = element_rect(fill = bg_color, colour = NA),
            panel.background = element_rect(fill = bg_color, colour = NA)
          )
        }

        # v179: Log tree stretch values (applied during rendering)
        if (tree_stretch_x != 1 || tree_stretch_y != 1) {
          debug_cat(paste0("\n=== v179: TREE STRETCH VALUES ===\n"))
          debug_cat(paste0("  X stretch (length): ", tree_stretch_x, "x\n"))
          debug_cat(paste0("  Y stretch (width): ", tree_stretch_y, "x\n"))
          debug_cat(paste0("  v79: Stretch will be applied via coord transformation\n"))
        }

        # Apply page title
        page_title_settings <- values$page_title

        # v131: DEBUG - trace page title settings
        debug_cat(paste0("\n=== v131: PAGE TITLE CHECK ===\n"))
        debug_cat(paste0("  page_title_settings is NULL: ", is.null(page_title_settings), "\n"))
        if (!is.null(page_title_settings)) {
          debug_cat(paste0("  enabled: ", page_title_settings$enabled, "\n"))
          debug_cat(paste0("  text: '", page_title_settings$text, "'\n"))
          debug_cat(paste0("  text length: ", nchar(page_title_settings$text), "\n"))
        }

        if (!is.null(page_title_settings) && isTRUE(page_title_settings$enabled) &&
            !is.null(page_title_settings$text) && nchar(page_title_settings$text) > 0) {
          debug_cat(paste0("\n=== v131: Applying page title ===\n"))
          debug_cat(paste0("  Title: ", page_title_settings$text, "\n"))
          debug_cat(paste0("  Size: ", page_title_settings$size, "\n"))
          debug_cat(paste0("  Color: ", page_title_settings$color, "\n"))

          fontface <- if (page_title_settings$bold) "bold" else "plain"
          hjust_val <- page_title_settings$hjust

          result <- result + labs(title = page_title_settings$text) +
            theme(
              plot.title = element_text(
                size = page_title_settings$size,
                colour = page_title_settings$color,
                face = fontface,
                hjust = hjust_val
              )
            )

          # Add underline if requested (using geom_segment as underline)
          if (isTRUE(page_title_settings$underline)) {
            debug_cat(paste0("  Adding underline\n"))
            # Note: Underline in ggplot title is complex, would need grid manipulation
            # For now, we document that underline is not fully supported
          }
          debug_cat(paste0("  Page title applied successfully\n"))
        } else {
          debug_cat(paste0("  Page title NOT applied (condition not met)\n"))
        }

        # v143: Apply custom text annotations as TRUE overlays using annotation_custom
        # annotation_custom with grid::textGrob is a true overlay that never affects plot limits
        # Unlike annotate() + coord_cartesian which can override ggtree's coordinate system
        custom_texts <- values$custom_texts
        if (!is.null(custom_texts) && length(custom_texts) > 0) {
          debug_cat(paste0("\n=== v143: Applying ", length(custom_texts), " custom text(s) as TRUE OVERLAY ===\n"))
          debug_cat(paste0("  Using annotation_custom with grid::textGrob (never affects plot limits)\n"))

          for (i in seq_along(custom_texts)) {
            txt <- custom_texts[[i]]

            debug_cat(paste0("  Text ", i, ": \"", substr(txt$content, 1, 20), "...\" at normalized (",
                                       round(txt$x, 2), ", ", round(txt$y, 2), ")\n"))

            # Create a text grob with proper formatting
            fontface_val <- switch(txt$fontface,
              "plain" = 1,
              "bold" = 2,
              "italic" = 3,
              "bold.italic" = 4,
              1  # default to plain
            )

            text_grob <- grid::textGrob(
              label = txt$content,
              x = grid::unit(txt$x, "npc"),
              y = grid::unit(txt$y, "npc"),
              gp = grid::gpar(
                fontsize = txt$size,
                col = txt$color,
                fontface = fontface_val
              ),
              hjust = txt$hjust,
              vjust = txt$vjust,
              rot = txt$angle
            )

            # annotation_custom with -Inf/Inf places in normalized page coordinates
            # The grob's own x/y units (npc) handle the actual positioning
            result <- result + annotation_custom(
              grob = text_grob,
              xmin = -Inf, xmax = Inf,
              ymin = -Inf, ymax = Inf
            )
          }

          debug_cat(paste0("  v43: Text overlays applied - plot coordinates unchanged\n"))
        }

        # v145: Store custom images for later application (after cowplot wrapping)
        # This ensures images are drawn on TOP of everything including cowplot canvas
        custom_images <- values$custom_images
        attr(result, "custom_images") <- custom_images
        if (!is.null(custom_images) && length(custom_images) > 0) {
          debug_cat(paste0("\n=== v145: ", length(custom_images), " custom image(s) queued for overlay ===\n"))
        }
      }, error = function(e) {
        debug_cat(paste0("  v30 Extra tab ERROR: ", e$message, "\n"))
      })

      # S1.3-PERF: This is now the ONLY call to func.move.tiplabels.to.front()
      # Layer reordering happens once here after all layers are added, not after each step
      # This reduces redundant calls from 3+ per render to just 1
      result <- func.move.tiplabels.to.front(result)

      # Store the plot with legend settings applied
      cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] generate_plot: storing plot in values$current_plot at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
      values$current_plot <- result

      cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT A: After storing plot, before legend coord extraction\n"))

      # v148: Extract and output all legend coordinates with coordinate system explanations
      tryCatch({
        debug_cat(paste0("\n=== v148: LEGEND COORDINATE SYSTEMS EXPLAINED ===\n"))
        debug_cat(paste0("\n"))
        debug_cat(paste0("NOTE: There are TWO different coordinate systems:\n"))
        debug_cat(paste0("\n"))
        debug_cat(paste0("1) PLOT COORDINATES (used by Highlight & Bootstrap legends):\n"))
        debug_cat(paste0("   - These are DATA coordinates within the plot area\n"))
        debug_cat(paste0("   - x = tree depth direction (after coord_flip: vertical position)\n"))
        debug_cat(paste0("   - y = tip position direction (after coord_flip: horizontal position)\n"))
        debug_cat(paste0("   - Negative x values place items below the tree baseline\n"))
        debug_cat(paste0("   - Controlled via: highlight_x_offset, highlight_y_offset,\n"))
        debug_cat(paste0("                     bootstrap_x_offset, bootstrap_y_offset in Legend tab\n"))
        debug_cat(paste0("\n"))
        debug_cat(paste0("2) GRID COORDINATES (used by Classification/Heatmap/P-value legends):\n"))
        debug_cat(paste0("   - These are LAYOUT CELL positions in the rendered gtable\n"))
        debug_cat(paste0("   - Not directly comparable to plot coordinates\n"))
        debug_cat(paste0("   - Positioned by ggplot's legend.position theme setting\n"))
        debug_cat(paste0("\n"))
        debug_cat(paste0("To ALIGN legends: Use the Legend tab offset controls to move\n"))
        debug_cat(paste0("Highlight/Bootstrap legends up/down/left/right until visually aligned.\n"))
        debug_cat(paste0("\n"))

        # Build the plot to extract grob information
        plot_build <- ggplot2::ggplot_build(result)
        plot_gtable <- ggplot2::ggplot_gtable(plot_build)

        # Get plot data range to help understand coordinate scale
        if (!is.null(plot_build$layout$panel_params) && length(plot_build$layout$panel_params) > 0) {
          pp <- plot_build$layout$panel_params[[1]]
          debug_cat(paste0("PLOT DATA RANGE (for reference):\n"))
          if (!is.null(pp$x.range)) {
            debug_cat(paste0("  x-axis range: ", round(pp$x.range[1], 2), " to ", round(pp$x.range[2], 2), "\n"))
          }
          if (!is.null(pp$y.range)) {
            debug_cat(paste0("  y-axis range: ", round(pp$y.range[1], 2), " to ", round(pp$y.range[2], 2), "\n"))
          }
          debug_cat(paste0("\n"))
        }

        debug_cat(paste0("GGPLOT LEGENDS (Grid Coordinates):\n"))
        # Find all legend grobs
        # S1.62dev: Added safer NA handling to prevent "missing value where TRUE/FALSE needed"
        layout_names <- plot_gtable$layout$name
        if (!is.null(layout_names) && length(layout_names) > 0) {
          # Use na.rm-safe grep that treats NA as FALSE
          legend_grobs <- which(sapply(layout_names, function(x) {
            !is.na(x) && grepl("guide-box", x)
          }))
        } else {
          legend_grobs <- integer(0)
        }
        if (length(legend_grobs) > 0) {
          for (leg_i in seq_along(legend_grobs)) {
            leg_idx <- legend_grobs[leg_i]
            leg_name <- plot_gtable$layout$name[leg_idx]
            leg_l <- plot_gtable$layout$l[leg_idx]
            leg_r <- plot_gtable$layout$r[leg_idx]
            leg_t <- plot_gtable$layout$t[leg_idx]
            leg_b <- plot_gtable$layout$b[leg_idx]
            debug_cat(paste0("  Legend Box ", leg_i, " ('", leg_name, "'):\n"))
            debug_cat(paste0("    Grid cell: column ", leg_l, "-", leg_r, ", row ", leg_t, "-", leg_b, "\n"))
          }
        } else {
          debug_cat(paste0("  No legend guide-boxes found in gtable\n"))
        }

        # Also extract legend titles from the built plot scales
        scales_info <- plot_build$plot$scales$scales
        if (!is.null(scales_info) && length(scales_info) > 0) {
          debug_cat(paste0("\n  Active Scales with Legends:\n"))
          for (scale_i in seq_along(scales_info)) {
            scale_obj <- scales_info[[scale_i]]
            # S1.62dev: Safer check for scale name - handle NA and NULL
            scale_name <- tryCatch(scale_obj$name, error = function(e) NULL)
            if (!is.null(scale_name) && !is.na(scale_name) && nchar(as.character(scale_name)) > 0) {
              aesthetics_str <- tryCatch(paste(scale_obj$aesthetics, collapse=", "), error = function(e) "unknown")
              debug_cat(paste0("    - '", scale_name, "' (", aesthetics_str, ")\n"))
            }
          }
        }

        debug_cat(paste0("\n=================================================\n"))
      }, error = function(e) {
        # S2.0-PERF: Suppress legend coord extraction errors - this is optional debug info
        # and the error "missing value where TRUE/FALSE needed" is common with some scale types
        # debug_cat(paste0("  v48: Error extracting legend coords: ", e$message, "\n"))
      })

      cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT B: After legend coord extraction, before temp file creation\n"))

      # Create a unique temp file with timestamp to force browser refresh
      # S2.0-PERF: Use SVG for faster preview rendering
      temp_plot_file <- file.path(tempdir(), paste0("shiny_plot_", Sys.getpid(), "_",
                                                    format(Sys.time(), "%Y%m%d_%H%M%S_%OS3"), ".svg"))
      # v53: cat(file=stderr(), "Temp file path:", temp_plot_file, "\n")
      
      # Clean up old plot files to avoid accumulation
      old_files <- list.files(tempdir(), pattern = paste0("^shiny_plot_", Sys.getpid()), full.names = TRUE)
      if (length(old_files) > 5) {  # Keep only last 5 files
        old_files_sorted <- old_files[order(file.info(old_files)$mtime)]
        files_to_delete <- old_files_sorted[1:(length(old_files_sorted) - 5)]
        sapply(files_to_delete, unlink)
      }

      cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT C: Before ggsave tryCatch, temp_plot_file=", temp_plot_file, "\n"))

      tryCatch({
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT D: Inside ggsave tryCatch\n"))
        # Save the plot as PNG
        # v145: Use the SAME proportions as the download tab settings
        # This ensures consistent appearance across all tabs (tree, bootstrap, heatmap, extra, download)
        # Get user's output dimensions and convert to inches for ggsave
        user_width <- if (!is.null(input$output_width) && !is.na(input$output_width)) input$output_width else 29.7
        user_height <- if (!is.null(input$output_height) && !is.na(input$output_height)) input$output_height else 21
        user_units <- if (!is.null(input$output_units)) input$output_units else "cm"

        # Convert to inches for ggsave
        if (user_units == "cm") {
          preview_width <- user_width / 2.54
          preview_height <- user_height / 2.54
        } else if (user_units == "mm") {
          preview_width <- user_width / 25.4
          preview_height <- user_height / 25.4
        } else {
          preview_width <- user_width
          preview_height <- user_height
        }

        # Cap preview size to reasonable limits while maintaining aspect ratio
        max_preview_dim <- 25  # Max 25 inches
        if (preview_width > max_preview_dim || preview_height > max_preview_dim) {
          scale_factor <- max_preview_dim / max(preview_width, preview_height)
          preview_width <- preview_width * scale_factor
          preview_height <- preview_height * scale_factor
        }

        debug_cat(paste0("  v45: Preview using DOWNLOAD TAB proportions\n"))
        debug_cat(paste0("  v45: User dimensions: ", user_width, " x ", user_height, " ", user_units, "\n"))
        debug_cat(paste0("  v45: Preview dimensions: ", round(preview_width, 2), " x ", round(preview_height, 2), " in\n"))

        # v146: Apply plot position offsets AND scale using cowplot for true transformation
        # v179: Also apply tree stretch (different x/y scaling)
        # This moves and scales the plot without squeezing/distorting proportions
        plot_to_save <- result
        offset_x <- attr(result, "plot_offset_x")
        offset_y <- attr(result, "plot_offset_y")
        scale_pct <- attr(result, "plot_scale_percent")
        if (is.null(scale_pct)) scale_pct <- 100
        # v179: Get tree stretch values
        tree_stretch_x <- attr(result, "tree_stretch_x")
        tree_stretch_y <- attr(result, "tree_stretch_y")
        if (is.null(tree_stretch_x)) tree_stretch_x <- 1
        if (is.null(tree_stretch_y)) tree_stretch_y <- 1

        # v180: Get keep_proportions setting
        keep_proportions <- attr(result, "keep_proportions")
        if (is.null(keep_proportions)) keep_proportions <- FALSE

        # v180: Calculate proportions adjustment when switching orientation
        # Default page is landscape (width > height, aspect ratio ~1.414)
        current_aspect <- preview_width / preview_height
        landscape_aspect <- 29.7 / 21  # A4 landscape ratio (~1.414)

        # v180: Check if we need proportions adjustment (portrait page with proportions preserved)
        proportion_adj_w <- 1
        proportion_adj_h <- 1
        if (isTRUE(keep_proportions) && current_aspect < 1) {
          # Portrait page - need to scale plot to fit landscape-shaped plot in portrait page
          # The plot should maintain landscape proportions (wider than tall)
          # Scale down to fit: width fits fully, height is proportionally smaller
          proportion_adj_w <- 1  # Plot spans full width
          proportion_adj_h <- current_aspect / landscape_aspect  # Scale height to maintain aspect ratio
          debug_cat(paste0("\n=== v180: PRESERVING PLOT PROPORTIONS ===\n"))
          debug_cat(paste0("  Current aspect: ", round(current_aspect, 3), " (portrait)\n"))
          debug_cat(paste0("  Landscape aspect: ", round(landscape_aspect, 3), "\n"))
          debug_cat(paste0("  Proportion adjustment: w=", round(proportion_adj_w, 3), ", h=", round(proportion_adj_h, 3), "\n"))
        }

        # Check if we need to apply any transformation
        needs_transform <- (!is.null(offset_x) && !is.null(offset_y) && (offset_x != 0 || offset_y != 0)) ||
                          scale_pct != 100 || tree_stretch_x != 1 || tree_stretch_y != 1 ||
                          proportion_adj_w != 1 || proportion_adj_h != 1

        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT E: needs_transform=", needs_transform, "\n"))

        if (needs_transform) {
          debug_cat(paste0("\n=== v180: APPLYING PLOT POSITION, SCALE, STRETCH AND PROPORTIONS WITH COWPLOT ===\n"))

          # Convert slider values to position offsets
          # X slider: -5 to 5 -> position offset of -0.25 to 0.25 (50% of canvas width total)
          # Y slider: -10 to 5 -> position offset of -0.5 to 0.25 (move down more than up)
          x_pos_offset <- if (!is.null(offset_x)) offset_x * 0.05 else 0  # Each unit moves 5% of canvas
          y_pos_offset <- if (!is.null(offset_y)) offset_y * 0.05 else 0  # Each unit moves 5% of canvas

          # v146: Convert scale percentage to draw_plot dimensions
          # 100% = width/height of 1 (full canvas)
          # 50% = width/height of 0.5 (half size, centered)
          # 200% = width/height of 2 (double size, centered)
          scale_factor <- scale_pct / 100

          # v179: Apply tree stretch to width and height separately
          # This allows stretching the tree longer (x) or wider (y) independently
          # v180: Also apply proportions adjustment
          final_width <- scale_factor * tree_stretch_x * proportion_adj_w
          final_height <- scale_factor * tree_stretch_y * proportion_adj_h

          # Calculate position to center the scaled plot
          # When scale_factor = 1, x = 0 + x_offset, y = 0 + y_offset (top-left)
          # When scale_factor = 0.5, x = 0.25 + x_offset (centered at 0.5)
          # Formula: x = (1 - width) / 2 + x_offset
          center_x <- (1 - final_width) / 2 + x_pos_offset
          center_y <- (1 - final_height) / 2 + y_pos_offset

          debug_cat(paste0("  Scale: ", scale_pct, "% (factor: ", round(scale_factor, 3), ")\n"))
          debug_cat(paste0("  v79: Tree stretch X: ", tree_stretch_x, "x, Y: ", tree_stretch_y, "x\n"))
          debug_cat(paste0("  v80: Proportion adj W: ", round(proportion_adj_w, 3), ", H: ", round(proportion_adj_h, 3), "\n"))
          debug_cat(paste0("  v80: Final dimensions: width=", round(final_width, 3), ", height=", round(final_height, 3), "\n"))
          debug_cat(paste0("  X position offset: ", round(x_pos_offset, 3), "\n"))
          debug_cat(paste0("  Y position offset: ", round(y_pos_offset, 3), "\n"))
          debug_cat(paste0("  Final position: (", round(center_x, 3), ", ", round(center_y, 3), ")\n"))

          # S1.62dev: Validate cowplot parameters to prevent potential crashes
          # Ensure dimensions are within reasonable bounds
          if (final_width < 0.01) {
            cat(file=stderr(), paste0("[WARN] final_width too small (", final_width, "), clamping to 0.01\n"))
            final_width <- 0.01
          }
          if (final_height < 0.01) {
            cat(file=stderr(), paste0("[WARN] final_height too small (", final_height, "), clamping to 0.01\n"))
            final_height <- 0.01
          }
          if (final_width > 3) {
            cat(file=stderr(), paste0("[WARN] final_width too large (", final_width, "), clamping to 3\n"))
            final_width <- 3
          }
          if (final_height > 3) {
            cat(file=stderr(), paste0("[WARN] final_height too large (", final_height, "), clamping to 3\n"))
            final_height <- 3
          }
          # Recalculate center after clamping
          center_x <- (1 - final_width) / 2 + x_pos_offset
          center_y <- (1 - final_height) / 2 + y_pos_offset

          # Use ggdraw to create a canvas and draw_plot to position and scale the plot
          # v180: Now includes proportions preservation for orientation changes
          cat(file=stderr(), paste0("[DEBUG] Creating cowplot canvas: x=", round(center_x,3), ", y=", round(center_y,3),
                                   ", w=", round(final_width,3), ", h=", round(final_height,3), "\n"))
          plot_to_save <- cowplot::ggdraw() +
            cowplot::draw_plot(result, x = center_x, y = center_y,
                              width = final_width, height = final_height)

          debug_cat(paste0("  v80: Plot wrapped in cowplot canvas with scale, stretch, proportions and offset positioning\n"))
        }

        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT F: After cowplot transformation\n"))

        # v145: Apply custom images as TRUE overlays using cowplot::draw_image
        # This is done AFTER cowplot wrapping to ensure images are on top
        custom_images <- attr(result, "custom_images")
        if (!is.null(custom_images) && length(custom_images) > 0) {
          debug_cat(paste0("\n=== v145: Applying ", length(custom_images), " custom image(s) as TRUE OVERLAY ===\n"))

          # Ensure we have a cowplot canvas
          if (!inherits(plot_to_save, "ggdraw")) {
            plot_to_save <- cowplot::ggdraw(plot_to_save)
          }

          for (i in seq_along(custom_images)) {
            img <- custom_images[[i]]
            if (file.exists(img$path)) {
              tryCatch({
                debug_cat(paste0("  Image ", i, ": ", img$name, "\n"))
                debug_cat(paste0("    Position: (", round(img$x, 2), ", ", round(img$y, 2), ")\n"))
                debug_cat(paste0("    Width: ", round(img$width, 3), "\n"))

                # Calculate height maintaining aspect ratio if not specified
                img_height <- img$height
                if (is.null(img_height) || img_height <= 0) {
                  # Read image to get aspect ratio
                  file_ext <- tolower(tools::file_ext(img$path))
                  if (file_ext %in% c("png")) {
                    img_data <- png::readPNG(img$path)
                  } else if (file_ext %in% c("jpg", "jpeg")) {
                    img_data <- jpeg::readJPEG(img$path)
                  } else {
                    img_data <- NULL
                  }
                  if (!is.null(img_data)) {
                    img_aspect <- dim(img_data)[1] / dim(img_data)[2]
                    img_height <- img$width * img_aspect
                    debug_cat(paste0("    Auto height: ", round(img_height, 3), " (aspect: ", round(img_aspect, 2), ")\n"))
                  } else {
                    img_height <- img$width  # Default to square if can't determine
                  }
                }

                # Use cowplot::draw_image for proper overlay
                plot_to_save <- plot_to_save +
                  cowplot::draw_image(img$path,
                                     x = img$x - img$width/2,  # Center horizontally
                                     y = img$y - img_height/2,  # Center vertically
                                     width = img$width,
                                     height = img_height)
                debug_cat(paste0("    Image added as overlay using cowplot::draw_image\n"))
              }, error = function(e) {
                debug_cat(paste0("  ERROR loading image ", i, ": ", e$message, "\n"))
              })
            } else {
              debug_cat(paste0("  Image ", i, ": File not found - ", img$path, "\n"))
            }
          }
          debug_cat(paste0("  v45: Custom images applied as overlays\n"))
        }

        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] CHECKPOINT G: After custom images, about to call ggsave\n"))

        # v53: cat(file=stderr(), "Calling ggsave...\n")
        # v54: Wrap in suppressWarnings to suppress scale warnings
        # S1.62dev: Added dimension logging for crash diagnosis
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] About to call ggsave at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   temp_plot_file=", temp_plot_file, "\n"))
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT]   preview dimensions: ", round(preview_width, 2), "x", round(preview_height, 2), " in\n"))

        # S1.62dev: Validate preview dimensions before ggsave
        if (is.na(preview_width) || is.na(preview_height) ||
            preview_width <= 0 || preview_height <= 0 ||
            preview_width > 100 || preview_height > 100) {
          cat(file=stderr(), paste0("[ERROR] Invalid preview dimensions! Using defaults.\n"))
          preview_width <- if (is.na(preview_width) || preview_width <= 0 || preview_width > 100) 11.69 else preview_width
          preview_height <- if (is.na(preview_height) || preview_height <= 0 || preview_height > 100) 8.27 else preview_height
        }

        suppressWarnings(ggsave(
          filename = temp_plot_file,
          plot = plot_to_save,
          width = preview_width,
          height = preview_height,
          units = "in",
          dpi = 150,
          limitsize = FALSE
        ))
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] ggsave completed at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
        
        # v53: cat(file=stderr(), "ggsave completed\n")
        
        # Wait a moment for file system to sync
        Sys.sleep(0.5)
        
        # v53: cat(file=stderr(), "File exists after ggsave:", file.exists(temp_plot_file), "\n")
        
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] Checking if file exists: ", temp_plot_file, "\n"))
        if (file.exists(temp_plot_file)) {
          cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] SUCCESS: File exists! Size=", file.info(temp_plot_file)$size, " bytes\n"))
          # v53: cat(file=stderr(), "File size:", file.info(temp_plot_file)$size, "bytes\n")

          # Store the file path - this will trigger renderImage
          values$temp_plot_file <- temp_plot_file
          values$plot_ready <- TRUE
          values$plot_counter <- values$plot_counter + 1  # Increment to force reactive update
          cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] plot_counter incremented to ", values$plot_counter, "\n"))

          # S2.0-PERF: Set cooldown timer to prevent rapid re-triggering
          last_plot_time(as.numeric(Sys.time()) * 1000)

          # Hide progress - plot is ready!
          values$progress_visible <- FALSE
          values$progress_message <- ""
          values$plot_generating <- FALSE  # Turn off generating indicator

          # v57: Show Ready status via shinyjs
          cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] Calling show_status_ready()\n"))
          show_status_ready()

          # v53: cat(file=stderr(), "Plot saved successfully and path stored\n")
          # v53: cat(file=stderr(), "Plot counter now:", values$plot_counter, "\n")
        } else {
          cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] ERROR: File does NOT exist after ggsave!\n"))
          # v53: cat(file=stderr(), "ERROR: File does not exist after ggsave\n")
          values$temp_plot_file <- NULL
          values$plot_ready <- FALSE
          values$progress_visible <- FALSE
          values$plot_generating <- FALSE
          # v57: Show click to generate on error
          show_status_click_to_generate()
        }

      }, error = function(e) {
        cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] ERROR in tryCatch: ", e$message, "\n"))
        # v53: cat(file=stderr(), "ERROR saving plot:", e$message, "\n")
        # v53: cat(file=stderr(), "Full error:\n")
        # v53: print(e)
        values$temp_plot_file <- NULL
        values$plot_ready <- FALSE
        values$progress_visible <- FALSE
        values$plot_generating <- FALSE
        # v57: Show click to generate on error
        show_status_click_to_generate()
      })
      
      # v53: debug_cat("==============================\n\n")
    }
    
    # Ensure plot_generating is always reset
    if (is.null(result)) {
      values$plot_generating <- FALSE
      values$progress_visible <- FALSE
    }
    
    # FAILSAFE: Always ensure plot_generating is reset
    values$plot_generating <- FALSE
    values$progress_visible <- FALSE

    # S2.0-PERF: Async conditional garbage collection (options 4A+4B)
    # - Only run gc() every 3 plots to reduce overhead
    # - Run asynchronously via later::later() to not block UI response
    # This saves ~0.1-0.3 sec per plot while still preventing memory accumulation
    # NOTE: Capture plot_counter in local var since later() callback runs outside reactive context
    # Use force() to ensure immediate evaluation before the closure is created
    current_plot_num <- values$plot_counter
    force(current_plot_num)  # Force immediate evaluation
    if (current_plot_num %% 3 == 0) {
      # Create callback with explicit local binding to avoid closure issues
      gc_callback <- local({
        plot_num <- current_plot_num  # Explicit local copy
        function() {
          gc_result <- gc(verbose = FALSE)
          mem_used_mb <- sum(gc_result[, 2])  # Used memory in MB
          cat(file=stderr(), paste0("[PERF-GC] Async gc() completed. Memory: ", round(mem_used_mb, 1), " MB (plot #", plot_num, ")\n"))
        }
      })
      later::later(gc_callback, delay = 0.1)  # 100ms delay to let UI update first
      cat(file=stderr(), "[PERF-GC] Scheduled async gc() (every 3rd plot)\n")
    } else {
      cat(file=stderr(), paste0("[PERF-GC] Skipped gc() (plot #", current_plot_num, ", runs every 3rd)\n"))
    }

    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] EXIT generate_plot() at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))
    cat(file=stderr(), paste0("[DEBUG-2ND-HIGHLIGHT] ========================================\n\n"))

    # S1.62dev: Flush stderr to ensure log messages are written before potential crash
    flush(stderr())

    # v53: cat(file=stderr(), "Finished generate_plot()\n")
  }  # End of generate_plot function
  
  
  
  # v57: tree_status_indicator replaced with static HTML + shinyjs for immediate updates
  # The status indicator is now controlled via show_status_waiting(), show_status_processing(),
  # show_status_ready(), and show_status_click_to_generate() helper functions
  
  # Output renderers (outside of generate_plot function)
  output$tree_preview <- renderImage({
    # S1.62dev: Added logging for crash diagnosis
    cat(file=stderr(), paste0("[RENDER] tree_preview renderImage called at ", format(Sys.time(), "%H:%M:%OS3"), "\n"))

    # Force reactive update by depending on plot_counter
    req(values$temp_plot_file, values$plot_counter)

    cat(file=stderr(), paste0("[RENDER] tree_preview: file=", values$temp_plot_file, ", counter=", values$plot_counter, "\n"))
    
    # v53: cat(file=stderr(), "Temp file:", values$temp_plot_file, "\n")
    # v53: cat(file=stderr(), "File exists:", file.exists(values$temp_plot_file), "\n")
    # v53: cat(file=stderr(), "Plot counter:", values$plot_counter, "\n")

    # S1.62dev: Added file existence check for crash diagnosis
    if (!file.exists(values$temp_plot_file)) {
      cat(file=stderr(), paste0("[RENDER] ERROR: temp_plot_file does not exist!\n"))
      return(NULL)
    }

    cat(file=stderr(), paste0("[RENDER] tree_preview returning image list\n"))
    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Tree plot preview"
    )
  }, deleteFile = FALSE)
  
  # v59: Removed renderUI status indicators for classification, bootstrap, highlight, and heatmap tabs
  # Status indicators are now static HTML elements toggled via shinyjs helper functions
  # (show_status_waiting, show_status_processing, show_status_ready, show_status_click_to_generate)

  output$classification_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    # v53: cat(file=stderr(), "\n=== renderImage called for classification_preview ===\n")
    
    req(values$temp_plot_file, values$plot_counter)
    
    # v53: cat(file=stderr(), "Temp file:", values$temp_plot_file, "\n")
    # v53: cat(file=stderr(), "File exists:", file.exists(values$temp_plot_file), "\n")
    # v53: cat(file=stderr(), "Plot counter:", values$plot_counter, "\n")
    
    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Classification preview"
    )
  }, deleteFile = FALSE)
  
  output$bootstrap_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    # v53: cat(file=stderr(), "\n=== renderImage called for bootstrap_preview ===\n")
    
    req(values$temp_plot_file, values$plot_counter)
    
    # v53: cat(file=stderr(), "Temp file:", values$temp_plot_file, "\n")
    # v53: cat(file=stderr(), "File exists:", file.exists(values$temp_plot_file), "\n")
    # v53: cat(file=stderr(), "Plot counter:", values$plot_counter, "\n")
    
    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Bootstrap preview"
    )
  }, deleteFile = FALSE)
  
  output$highlight_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    # v53: cat(file=stderr(), "\n=== renderImage called for highlight_preview ===\n")
    
    req(values$temp_plot_file, values$plot_counter)
    
    # v53: cat(file=stderr(), "Temp file:", values$temp_plot_file, "\n")
    # v53: cat(file=stderr(), "File exists:", file.exists(values$temp_plot_file), "\n")
    # v53: cat(file=stderr(), "Plot counter:", values$plot_counter, "\n")
    
    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Highlight preview"
    )
  }, deleteFile = FALSE)
  
  output$heatmap_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    # v53: cat(file=stderr(), "\n=== renderImage called for heatmap_preview ===\n")
    
    req(values$temp_plot_file, values$plot_counter)
    
    # v53: cat(file=stderr(), "Temp file:", values$temp_plot_file, "\n")
    # v53: cat(file=stderr(), "File exists:", file.exists(values$temp_plot_file), "\n")
    # v53: cat(file=stderr(), "Plot counter:", values$plot_counter, "\n")
    
    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Heatmap preview"
    )
  }, deleteFile = FALSE)
  
  # v125: Final preview - using renderImage like other tabs for consistent aspect ratio
  output$final_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    req(values$temp_plot_file, values$plot_counter)

    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Final preview"
    )
  }, deleteFile = FALSE)

  # v124: Legend preview - using renderImage like other tabs for consistent aspect ratio
  output$legend_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    req(values$temp_plot_file, values$plot_counter)

    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Legend preview"
    )
  }, deleteFile = FALSE)

  # v130: Extra tab preview
  output$extra_preview <- renderImage({
    # Force reactive update by depending on plot_counter
    req(values$temp_plot_file, values$plot_counter)

    list(
      src = values$temp_plot_file,
      contentType = "image/svg+xml",
      width = "100%",
      alt = "Extra preview"
    )
  }, deleteFile = FALSE)

  # v130: Extra tab status updates when data is loaded
  observe({
    if (!is.null(values$tree) && !is.null(values$csv_data)) {
      shinyjs::hide("extra_status_waiting")
      shinyjs::hide("extra_status_processing")
      shinyjs::show("extra_status_ready")
    } else {
      shinyjs::show("extra_status_waiting")
      shinyjs::hide("extra_status_processing")
      shinyjs::hide("extra_status_ready")
    }
  })

  # v130: Observer for page title settings
  observeEvent({
    input$enable_page_title
    input$page_title_text
    input$page_title_x
    input$page_title_y
    input$page_title_size
    input$page_title_color
    input$page_title_bold
    input$page_title_underline
    input$page_title_hjust
  }, {
    values$page_title <- list(
      enabled = isTRUE(input$enable_page_title),
      text = if(!is.null(input$page_title_text)) input$page_title_text else "",
      x = if(!is.null(input$page_title_x)) input$page_title_x else 0.5,
      y = if(!is.null(input$page_title_y)) input$page_title_y else 0.95,
      size = if(!is.null(input$page_title_size)) input$page_title_size else 18,
      color = if(!is.null(input$page_title_color)) input$page_title_color else "#000000",
      bold = isTRUE(input$page_title_bold),
      underline = isTRUE(input$page_title_underline),
      hjust = if(!is.null(input$page_title_hjust)) as.numeric(input$page_title_hjust) else 0.5
    )
  }, ignoreInit = TRUE)

  # v130: Add custom text annotation
  observeEvent(input$add_custom_text, {
    req(input$custom_text_content)
    if (nchar(trimws(input$custom_text_content)) > 0) {
      new_text <- list(
        id = paste0("text_", length(values$custom_texts) + 1, "_", Sys.time()),
        content = input$custom_text_content,
        x = if(!is.null(input$custom_text_x)) input$custom_text_x else 0.5,
        y = if(!is.null(input$custom_text_y)) input$custom_text_y else 0.5,
        size = if(!is.null(input$custom_text_size)) input$custom_text_size else 12,
        color = if(!is.null(input$custom_text_color)) input$custom_text_color else "#000000",
        fontface = if(!is.null(input$custom_text_fontface)) input$custom_text_fontface else "plain",
        hjust = if(!is.null(input$custom_text_hjust)) as.numeric(input$custom_text_hjust) else 0.5,
        vjust = if(!is.null(input$custom_text_vjust)) as.numeric(input$custom_text_vjust) else 0.5,
        angle = if(!is.null(input$custom_text_angle)) input$custom_text_angle else 0
      )
      values$custom_texts <- c(values$custom_texts, list(new_text))
      # Clear the text input
      updateTextInput(session, "custom_text_content", value = "")
    }
  })

  # v130: Clear all custom texts
  observeEvent(input$clear_custom_texts, {
    values$custom_texts <- list()
  })

  # v130: Render list of custom texts
  output$custom_texts_list <- renderUI({
    if (length(values$custom_texts) == 0) {
      return(tags$p(class = "text-muted", "No custom texts added yet."))
    }

    text_items <- lapply(seq_along(values$custom_texts), function(i) {
      txt <- values$custom_texts[[i]]
      tags$div(
        style = "padding: 5px; margin: 5px 0; border: 1px solid #ddd; border-radius: 4px; background-color: #f9f9f9;",
        tags$span(style = paste0("color: ", txt$color, "; font-weight: ",
                                  if(txt$fontface == "bold" || txt$fontface == "bold.italic") "bold" else "normal", ";"),
                  paste0(i, ". \"", substr(txt$content, 1, 30), if(nchar(txt$content) > 30) "..." else "", "\"")),
        tags$span(class = "text-muted",
                  paste0(" (x:", round(txt$x, 2), ", y:", round(txt$y, 2), ", size:", txt$size, ")")),
        actionButton(paste0("delete_text_", i), "", icon = icon("times"),
                     class = "btn-xs btn-danger", style = "float: right; padding: 2px 6px;")
      )
    })

    do.call(tagList, text_items)
  })

  # v130: Delete individual custom text (dynamic observers)
  observe({
    lapply(seq_along(values$custom_texts), function(i) {
      observeEvent(input[[paste0("delete_text_", i)]], {
        if (i <= length(values$custom_texts)) {
          values$custom_texts <- values$custom_texts[-i]
        }
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # v130: Add custom image
  observeEvent(input$add_custom_image, {
    req(input$custom_image_file)

    # Copy file to temp location to persist it
    temp_path <- file.path(tempdir(), paste0("custom_img_", length(values$custom_images) + 1, "_",
                                              basename(input$custom_image_file$name)))
    file.copy(input$custom_image_file$datapath, temp_path, overwrite = TRUE)

    new_image <- list(
      id = paste0("img_", length(values$custom_images) + 1, "_", Sys.time()),
      path = temp_path,
      name = input$custom_image_file$name,
      x = if(!is.null(input$custom_image_x)) input$custom_image_x else 0.5,
      y = if(!is.null(input$custom_image_y)) input$custom_image_y else 0.5,
      width = if(!is.null(input$custom_image_width)) input$custom_image_width else 0.2,
      height = if(!is.null(input$custom_image_height) && input$custom_image_height > 0) input$custom_image_height else NULL
    )
    values$custom_images <- c(values$custom_images, list(new_image))
  })

  # v130: Clear all custom images
  observeEvent(input$clear_custom_images, {
    values$custom_images <- list()
  })

  # v130: Render list of custom images
  output$custom_images_list <- renderUI({
    if (length(values$custom_images) == 0) {
      return(tags$p(class = "text-muted", "No custom images added yet."))
    }

    img_items <- lapply(seq_along(values$custom_images), function(i) {
      img <- values$custom_images[[i]]
      tags$div(
        style = "padding: 5px; margin: 5px 0; border: 1px solid #ddd; border-radius: 4px; background-color: #f9f9f9;",
        tags$span(icon("image"), paste0(i, ". ", img$name)),
        tags$span(class = "text-muted",
                  paste0(" (x:", round(img$x, 2), ", y:", round(img$y, 2),
                         ", w:", round(img$width, 2),
                         if(!is.null(img$height)) paste0(", h:", round(img$height, 2)) else "", ")")),
        actionButton(paste0("delete_image_", i), "", icon = icon("times"),
                     class = "btn-xs btn-danger", style = "float: right; padding: 2px 6px;")
      )
    })

    do.call(tagList, img_items)
  })

  # v130: Delete individual custom image (dynamic observers)
  observe({
    lapply(seq_along(values$custom_images), function(i) {
      observeEvent(input[[paste0("delete_image_", i)]], {
        if (i <= length(values$custom_images)) {
          values$custom_images <- values$custom_images[-i]
        }
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # v141: Observer for plot position X slider
  # S1-PERF: Using debounced version to prevent rapid updates
  observeEvent(plot_offset_x_d(), {
    values$plot_offset_x <- plot_offset_x_d()
  }, ignoreInit = TRUE)

  # v141: Observer for plot position Y slider
  # S1-PERF: Using debounced version to prevent rapid updates
  observeEvent(plot_offset_y_d(), {
    values$plot_offset_y <- plot_offset_y_d()
  }, ignoreInit = TRUE)

  # v141: Reset plot position button
  observeEvent(input$reset_plot_position, {
    updateSliderInput(session, "plot_offset_x", value = 0)
    updateSliderInput(session, "plot_offset_y", value = 0)
    values$plot_offset_x <- 0
    values$plot_offset_y <- 0
  })

  # v146: Observer for plot scale slider
  # S1-PERF: Using debounced version to prevent rapid updates
  observeEvent(plot_scale_percent_d(), {
    values$plot_scale_percent <- plot_scale_percent_d()
  }, ignoreInit = TRUE)

  # v146: Reset plot scale button
  observeEvent(input$reset_plot_scale, {
    updateSliderInput(session, "plot_scale_percent", value = 100)
    values$plot_scale_percent <- 100
  })

  # v179: Observer for tree stretch X slider (horizontal length)
  # S1-PERF: Using debounced version to prevent rapid updates
  observeEvent(tree_stretch_x_d(), {
    values$tree_stretch_x <- tree_stretch_x_d()
  }, ignoreInit = TRUE)

  # v179: Observer for tree stretch Y slider (vertical width)
  # S1-PERF: Using debounced version to prevent rapid updates
  observeEvent(tree_stretch_y_d(), {
    values$tree_stretch_y <- tree_stretch_y_d()
  }, ignoreInit = TRUE)

  # v179: Reset tree stretch button
  observeEvent(input$reset_tree_stretch, {
    updateSliderInput(session, "tree_stretch_x", value = 1)
    updateSliderInput(session, "tree_stretch_y", value = 1)
    values$tree_stretch_x <- 1
    values$tree_stretch_y <- 1
  })

  # v179: Observer for background color
  # S1.62dev: Use debounced version to prevent lag from rapid color picker changes
  # Also trigger automatic plot regeneration so users don't need to click Apply
  observeEvent(background_color_d(), {
    values$background_color <- background_color_d()
    # S1.62dev: Auto-regenerate plot when background color changes
    if (isTRUE(values$plot_ready)) {
      request_plot_update()
    }
  }, ignoreInit = TRUE)

  # v179: Reset background color button
  observeEvent(input$reset_background, {
    colourpicker::updateColourInput(session, "background_color", value = "#FFFFFF")
    values$background_color <- "#FFFFFF"
  })

  # v130: Apply Extra settings to plot
  # v139: Added processing indicator
  observeEvent(input$extra_apply, {
    req(values$plot_ready)

    # Show processing indicator
    shinyjs::hide("extra_status_waiting")
    shinyjs::show("extra_status_processing")
    shinyjs::hide("extra_status_ready")

    # Use a slight delay to ensure UI updates before plot generation
    shinyjs::delay(50, {
      generate_plot()

      # Show ready indicator after plot generation
      shinyjs::hide("extra_status_waiting")
      shinyjs::hide("extra_status_processing")
      shinyjs::show("extra_status_ready")
    })
  })

  ###################

  # S1.62dev: Define YAML content reactive - COMPLETE export of all visual settings
  yaml_content <- reactive({
    # Check if we have the necessary data
    req(values$tree, values$csv_data)

    # S1.62dev: Build classification list from values$classifications
    classification_list <- list()
    if (!is.null(values$classifications) && length(values$classifications) > 0) {
      for (i in seq_along(values$classifications)) {
        class_def <- values$classifications[[i]]

        # Build the "according" list for this classification
        according_list <- list()
        if (!is.null(class_def$classes) && length(class_def$classes) > 0) {
          for (j in seq_along(class_def$classes)) {
            cls <- class_def$classes[[j]]
            according_list[[j]] <- list()
            according_list[[j]][[as.character(j)]] <- list(
              title1 = if (!is.null(cls$column)) cls$column else "",
              value1 = if (!is.null(cls$value)) cls$value else "",
              display_name = if (!is.null(cls$display_name)) cls$display_name else "",
              color = if (!is.null(cls$color)) cls$color else "#000000"
            )
          }
        }

        # Build highlight settings for this classification
        highlight_list <- list(display = "no")
        if (!is.null(class_def$highlight) && isTRUE(class_def$highlight$enabled)) {
          highlight_list$display <- "yes"
          highlight_list$according <- list()
          if (!is.null(class_def$highlight$items) && length(class_def$highlight$items) > 0) {
            for (j in seq_along(class_def$highlight$items)) {
              hi <- class_def$highlight$items[[j]]
              highlight_list$according[[j]] <- list()
              highlight_list$according[[j]][[as.character(j)]] <- list(
                title1 = if (!is.null(hi$column)) hi$column else "",
                value1 = if (!is.null(hi$value)) hi$value else "",
                display_name = if (!is.null(hi$display_name)) hi$display_name else "",
                color = if (!is.null(hi$color)) hi$color else "#FF0000"
              )
            }
          }
        }

        classification_list[[i]] <- list()
        classification_list[[i]][[as.character(i)]] <- list(
          title = if (!is.null(class_def$title)) class_def$title else "Classification",
          column = if (!is.null(class_def$column)) class_def$column else "",
          FDR_perc = if (!is.null(class_def$fdr)) class_def$fdr else 0.1,
          non_cluster_title = if (!is.null(class_def$no_cluster_title)) class_def$no_cluster_title else "No cluster",
          non_cluster_color = if (!is.null(class_def$no_cluster_color)) class_def$no_cluster_color else "gray",
          according = according_list,
          highlight = highlight_list
        )
      }
    }

    # S1.62dev: Build heatmap list from values$heatmap_configs
    heatmap_list <- list()
    if (!is.null(values$heatmap_configs) && length(values$heatmap_configs) > 0) {
      for (i in seq_along(values$heatmap_configs)) {
        cfg <- values$heatmap_configs[[i]]

        heatmap_list[[i]] <- list(
          display = "yes",
          title = if (!is.null(cfg$title)) cfg$title else paste0("Heatmap ", i),
          # S1.62dev: Use actual_type (computed from auto-detect) if available
          # If auto_type is TRUE but actual_type not computed yet, default to "no" (continuous)
          # Only use cfg$type when auto_type is FALSE (user explicitly set the type)
          is_discrete = if (!is.null(cfg$actual_type)) {
            if (cfg$actual_type == "discrete") "yes" else "no"
          } else if (is.null(cfg$auto_type) || cfg$auto_type == TRUE) {
            # auto_type is enabled but actual_type not computed - can't reliably determine
            # Default to "no" (continuous) as safer fallback
            "no"
          } else if (!is.null(cfg$type) && cfg$type == "discrete") {
            "yes"
          } else {
            "no"
          },
          # S1.62dev: Export auto_type flag so import knows whether to re-detect
          auto_type = if (!is.null(cfg$auto_type) && cfg$auto_type) "yes" else "no",
          columns = if (!is.null(cfg$columns)) as.list(cfg$columns) else list(),
          distance = if (!is.null(cfg$distance)) cfg$distance else 0.02,
          height = if (!is.null(cfg$height)) cfg$height else 0.8,
          show_colnames = if (!is.null(cfg$show_colnames) && cfg$show_colnames) "yes" else "no",
          colnames_angle = if (!is.null(cfg$colnames_angle)) cfg$colnames_angle else 45,
          discrete_palette = if (!is.null(cfg$discrete_palette)) cfg$discrete_palette else "Set1",
          cont_palette = if (!is.null(cfg$cont_palette)) cfg$cont_palette else "Blues",
          low_color = if (!is.null(cfg$low_color)) cfg$low_color else "#FFFFCC",
          high_color = if (!is.null(cfg$high_color)) cfg$high_color else "#006837",
          use_midpoint = if (!is.null(cfg$use_midpoint) && cfg$use_midpoint) "yes" else "no",
          mid_color = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFF99",
          midpoint = if (!is.null(cfg$midpoint)) cfg$midpoint else 0,
          # S1.61dev: Guide line settings
          show_guide_lines = if (!is.null(cfg$show_guide_lines) && cfg$show_guide_lines) "yes" else "no",
          guide_color1 = if (!is.null(cfg$guide_color1)) cfg$guide_color1 else "#CCCCCC",
          guide_color2 = if (!is.null(cfg$guide_color2)) cfg$guide_color2 else "#EEEEEE",
          guide_alpha = if (!is.null(cfg$guide_alpha)) cfg$guide_alpha else 0.3,
          guide_width = if (!is.null(cfg$guide_width)) cfg$guide_width else 0.5,
          guide_linetype = if (!is.null(cfg$guide_linetype)) cfg$guide_linetype else "solid",
          # S1.62dev: Row label settings (were missing from export)
          show_row_labels = if (!is.null(cfg$show_row_labels) && cfg$show_row_labels) "yes" else "no",
          row_label_source = if (!is.null(cfg$row_label_source)) cfg$row_label_source else "colnames",
          row_label_font_size = if (!is.null(cfg$row_label_font_size)) cfg$row_label_font_size else 2.5,
          row_label_offset = if (!is.null(cfg$row_label_offset)) cfg$row_label_offset else 1.0,
          row_label_align = if (!is.null(cfg$row_label_align)) cfg$row_label_align else "left",
          custom_row_labels = if (!is.null(cfg$custom_row_labels)) cfg$custom_row_labels else "",
          label_mapping = if (!is.null(cfg$label_mapping)) cfg$label_mapping else list(),
          # S1.62dev: Grid settings (were missing from export)
          show_grid = if (!is.null(cfg$show_grid) && cfg$show_grid) "yes" else "no",
          grid_color = if (!is.null(cfg$grid_color)) cfg$grid_color else "#000000",
          grid_size = if (!is.null(cfg$grid_size)) cfg$grid_size else 0.5,
          # S1.62dev: Discrete custom color settings (were missing from export)
          custom_discrete = if (!is.null(cfg$custom_discrete) && cfg$custom_discrete) "yes" else "no",
          custom_colors = if (!is.null(cfg$custom_colors) && length(cfg$custom_colors) > 0) cfg$custom_colors else list(),
          # S2.8: Row line settings (horizontal lines) - were missing from export
          show_row_lines = if (!is.null(cfg$show_row_lines) && cfg$show_row_lines) "yes" else "no",
          row_line_color = if (!is.null(cfg$row_line_color)) cfg$row_line_color else "#000000",
          row_line_size = if (!is.null(cfg$row_line_size)) cfg$row_line_size else 0.5,
          # S2.8: Column line settings (vertical lines) - were missing from export
          show_col_lines = if (!is.null(cfg$show_col_lines) && cfg$show_col_lines) "yes" else "no",
          col_line_color = if (!is.null(cfg$col_line_color)) cfg$col_line_color else "#000000",
          col_line_size = if (!is.null(cfg$col_line_size)) cfg$col_line_size else 0.5,
          # S2.8: RData heatmap settings - were missing from export
          data_source = if (!is.null(cfg$data_source)) cfg$data_source else "csv",
          rdata_mapping_column = if (!is.null(cfg$rdata_mapping_column)) cfg$rdata_mapping_column else "",
          # S2.8: WGD normalization settings - were missing from export
          cnv_wgd_norm = if (!is.null(cfg$cnv_wgd_norm) && cfg$cnv_wgd_norm) "yes" else "no",
          cnv_wgd_per_cell = if (!is.null(cfg$cnv_wgd_per_cell) && cfg$cnv_wgd_per_cell) "yes" else "no",
          cnv_wgd_column = if (!is.null(cfg$cnv_wgd_column)) cfg$cnv_wgd_column else "",
          # S2.8: Display mode (basic or detailed) - for RData heatmaps
          cnv_display_mode = if (!is.null(cfg$cnv_display_mode)) cfg$cnv_display_mode else "basic",
          # Render downsample factor (second stage, 10 = default)
          cnv_render_downsample = if (!is.null(cfg$cnv_render_downsample)) cfg$cnv_render_downsample else 10,
          # Height scale for detailed mode
          cnv_height_scale = if (!is.null(cfg$cnv_height_scale)) cfg$cnv_height_scale else 1.0
        )
      }
    }

    # S1.62dev: Build highlights list from values$highlights
    # Also include temp_highlight_preview if it exists (unsaved preview)
    highlights_list <- list()

    # First add saved highlights
    if (!is.null(values$highlights) && length(values$highlights) > 0) {
      for (i in seq_along(values$highlights)) {
        hi <- values$highlights[[i]]

        # Build items list
        items_list <- list()
        if (!is.null(hi$items) && length(hi$items) > 0) {
          for (j in seq_along(hi$items)) {
            item <- hi$items[[j]]
            items_list[[j]] <- list(
              column = if (!is.null(item$column)) item$column else "",
              value = if (!is.null(item$value)) item$value else "",
              display_name = if (!is.null(item$display_name)) item$display_name else "",
              color = if (!is.null(item$color)) item$color else "#FF0000",
              transparency = if (!is.null(item$transparency)) item$transparency else 0.5
            )
          }
        }

        highlights_list[[i]] <- list(
          enabled = if (!is.null(hi$enabled) && hi$enabled) "yes" else "no",
          title = if (!is.null(hi$title)) hi$title else "Highlight",
          column = if (!is.null(hi$column)) hi$column else "",
          offset = if (!is.null(hi$offset)) hi$offset else 0,
          vertical_offset = if (!is.null(hi$vertical_offset)) hi$vertical_offset else 0,
          adjust_height = if (!is.null(hi$adjust_height)) hi$adjust_height else 1,
          adjust_width = if (!is.null(hi$adjust_width)) hi$adjust_width else 1,
          items = items_list
        )
      }
    }

    # S1.62dev: Also include temp_highlight_preview if it exists (unsaved preview work)
    if (!is.null(values$temp_highlight_preview) &&
        !is.null(values$temp_highlight_preview$items) &&
        length(values$temp_highlight_preview$items) > 0) {
      hi <- values$temp_highlight_preview

      # Build items list for preview
      items_list <- list()
      for (j in seq_along(hi$items)) {
        item <- hi$items[[j]]
        items_list[[j]] <- list(
          column = if (!is.null(item$column)) item$column else "",
          value = if (!is.null(item$value)) item$value else "",
          display_name = if (!is.null(item$display_name)) item$display_name else "",
          color = if (!is.null(item$color)) item$color else "#FF0000",
          transparency = if (!is.null(item$transparency)) item$transparency else 0.5
        )
      }

      # Add preview as a highlight (marked with title suffix if not already saved)
      preview_entry <- list(
        enabled = "yes",
        title = if (!is.null(hi$title)) hi$title else "Highlight (Unsaved)",
        column = if (!is.null(hi$column)) hi$column else "",
        offset = if (!is.null(hi$offset)) hi$offset else 0,
        vertical_offset = if (!is.null(hi$vertical_offset)) hi$vertical_offset else 0,
        adjust_height = if (!is.null(hi$adjust_height)) hi$adjust_height else 1,
        adjust_width = if (!is.null(hi$adjust_width)) hi$adjust_width else 1,
        items = items_list
      )
      highlights_list <- c(highlights_list, list(preview_entry))
    }

    # S1.62dev: Build legend settings from values$legend_settings
    # S2.292dev: Added font_family to export
    legend_settings_yaml <- list(
      position = if (!is.null(values$legend_settings$position)) values$legend_settings$position else "right",
      show_classification = if (!is.null(values$legend_settings$show_classification) && values$legend_settings$show_classification) "yes" else "no",
      show_highlight = if (!is.null(values$legend_settings$show_highlight) && values$legend_settings$show_highlight) "yes" else "no",
      show_bootstrap = if (!is.null(values$legend_settings$show_bootstrap) && values$legend_settings$show_bootstrap) "yes" else "no",
      # S1.62dev: Added show_pvalue to export
      show_pvalue = if (!is.null(values$legend_settings$show_pvalue) && values$legend_settings$show_pvalue) "yes" else "no",
      show_heatmap = if (!is.null(values$legend_settings$show_heatmap) && values$legend_settings$show_heatmap) "yes" else "no",
      title_size = if (!is.null(values$legend_settings$title_size)) values$legend_settings$title_size else 12,
      text_size = if (!is.null(values$legend_settings$text_size)) values$legend_settings$text_size else 10,
      font_family = if (!is.null(values$legend_settings$font_family)) values$legend_settings$font_family else "sans",
      box_background = if (!is.null(values$legend_settings$box_background)) values$legend_settings$box_background else "transparent",
      margin = if (!is.null(values$legend_settings$margin)) values$legend_settings$margin else 0.2
    )

    # Create the YAML structure based on current settings
    yaml_data <- list(
      "Individual general definitions" = list(
        Individual = input$individual_name,
        "individual column" = input$individual_column,
        "tree path" = if (!is.null(input$tree_file)) {
          list(input$tree_file$datapath)
        } else {
          list(NULL)
        },
        "mapping csv file" = if (!is.null(input$csv_file)) {
          input$csv_file$datapath
        } else {
          NULL
        },
        "out_file" = list(
          "base_path" = input$output_path,
          "file_type" = input$output_format,
          # S1.62dev: Export page orientation and dimensions
          "page_orientation" = if (!is.null(input$page_orientation)) input$page_orientation else "landscape",
          "output_width" = if (!is.null(input$output_width)) input$output_width else 29.7,
          "output_height" = if (!is.null(input$output_height)) input$output_height else 21,
          "output_units" = if (!is.null(input$output_units)) input$output_units else "cm",
          "keep_proportions" = if (!is.null(input$keep_proportions) && input$keep_proportions) "yes" else "no",
          "optional text at beggining" = input$prefix_text,
          "optional text at end" = input$suffix_text,
          "replace name" = list(
            flag = if (input$replace_name) "yes" else "no",
            name = input$custom_name
          )
        )
      ),
      "Mapping exl renaming titles" = list(
        "ID column" = input$id_column
      ),
      "visual definitions" = list(
        "classification" = classification_list,
        "heatmaps" = heatmap_list,
        "highlights" = highlights_list,
        "legend" = legend_settings_yaml,
        "Bootstrap" = list(
          display = if (input$show_bootstrap) "yes" else "no",
          format = input$bootstrap_format,
          param = as.character(input$bootstrap_param),
          label_size = if (!is.null(input$bootstrap_label_size)) input$bootstrap_label_size else 1.5
        ),
        "node_numbers" = list(
          display = if (!is.null(input$display_node_numbers) && input$display_node_numbers) "yes" else "no",
          font_size = if (!is.null(input$node_number_font_size)) input$node_number_font_size else 3.5
        ),
        "rotation1" = list(
          display = if (input$enable_rotation &&
                        (input$rotation_type == "primary" || input$rotation_type == "manual")) "yes" else "no",
          # S1.62dev: Fixed - read from rotation1_config instead of rotation_settings$primary
          according = if (!is.null(values$rotation1_config) && length(values$rotation1_config) > 0) {
            lapply(values$rotation1_config, function(r) list(col = r$col, val = r$val))
          } else list()
        ),
        "rotation2" = list(
          display = if (input$enable_rotation && input$rotation_type == "secondary") "yes" else "no",
          # S1.62dev: Fixed - read from rotation2_config instead of rotation_settings$secondary
          according = if (!is.null(values$rotation2_config) && length(values$rotation2_config) > 0) {
            lapply(values$rotation2_config, function(r) list(col = r$col, val = r$val))
          } else list()
        ),
        "manual_rotation" = list(
          display = if (!is.null(input$enable_rotation) && input$enable_rotation &&
                        !is.null(input$rotation_type) && input$rotation_type == "manual" &&
                        !is.null(values$manual_rotation_config) && length(values$manual_rotation_config) > 0) "yes" else "no",
          nodes = if (!is.null(values$manual_rotation_config) && length(values$manual_rotation_config) > 0 &&
                      !all(is.na(values$manual_rotation_config))) {
            as.list(values$manual_rotation_config)
          } else list()
        ),
        "trim tips" = list(
          display = if (input$trim_tips) "yes" else "no",
          length = input$tip_length
        ),
        "edge_width_multiplier" = list(
          size = input$edge_width
        ),
        "font_size" = list(
          tips = input$tip_font_size,
          legend_title = if (!is.null(values$legend_settings$title_size)) values$legend_settings$title_size else 13,
          legend_text = if (!is.null(values$legend_settings$text_size)) values$legend_settings$text_size else 10,
          heat_map_legend = input$heatmap_font_size
        ),
        # S1.62dev: Extra tab settings (plot position, scale, stretch)
        "extra_settings" = list(
          plot_offset_x = if (!is.null(values$plot_offset_x)) values$plot_offset_x else 0,
          plot_offset_y = if (!is.null(values$plot_offset_y)) values$plot_offset_y else 0,
          plot_scale = if (!is.null(input$plot_scale)) input$plot_scale else 100,
          tree_stretch_x = if (!is.null(input$tree_stretch_x)) input$tree_stretch_x else 1,
          tree_stretch_y = if (!is.null(input$tree_stretch_y)) input$tree_stretch_y else 1,
          background_color = if (!is.null(input$background_color)) input$background_color else "white"
        )
      )
    )

    # Convert to YAML text
    settings_to_yaml(yaml_data)
  })
  
  # Output YAML content
  output$yaml_output <- renderText({
    yaml_content()
  })
  
  # Download YAML configuration
  output$download_yaml <- downloadHandler(
    filename = function() {
      paste0(input$individual_name, "_config.yaml")
    },
    content = function(file) {
      # Use the reactive expression
      writeLines(yaml_content(), file)
    }
  )
  
  # Download YAML configuration from configuration tab
  output$download_yaml_config <- downloadHandler(
    filename = function() {
      paste0(input$individual_name, "_config.yaml")
    },
    content = function(file) {
      # Use the reactive expression
      writeLines(yaml_content(), file)
    }
  )
  
  # Download plot - v134: Use the actual generated plot (with heatmaps) instead of basic tree
  output$download_plot <- downloadHandler(
    filename = function() {
      if (input$replace_name) {
        paste0(input$custom_name, ".", input$output_format)
      } else {
        paste0(input$prefix_text, "_", input$individual_name, "_", input$suffix_text, ".", input$output_format)
      }
    },
    content = function(file) {
      # v134: Use the actual generated plot stored in values$current_plot
      # This includes all layers: tree, heatmaps, highlights, bootstrap, etc.

      # Check if we have a valid generated plot
      if (!is.null(values$current_plot)) {
        # Use ggsave for ggplot objects - handles all formats correctly

        # Determine width/height - convert to inches based on units
        width_val <- input$output_width
        height_val <- input$output_height

        # Convert to inches if needed (ggsave needs consistent units)
        if (input$output_units == "cm") {
          width_in <- width_val / 2.54
          height_in <- height_val / 2.54
        } else if (input$output_units == "mm") {
          width_in <- width_val / 25.4
          height_in <- height_val / 25.4
        } else {
          width_in <- width_val
          height_in <- height_val
        }

        # Set DPI based on format
        dpi_val <- if (input$output_format %in% c("pdf", "svg")) 300 else 300

        # v146: Apply cowplot positioning AND scale for downloads too
        # v180: Also apply tree stretch and proportions preservation
        plot_to_download <- values$current_plot
        offset_x <- attr(values$current_plot, "plot_offset_x")
        offset_y <- attr(values$current_plot, "plot_offset_y")
        scale_pct <- attr(values$current_plot, "plot_scale_percent")
        if (is.null(scale_pct)) scale_pct <- 100
        # v179: Get tree stretch values
        tree_stretch_x <- attr(values$current_plot, "tree_stretch_x")
        tree_stretch_y <- attr(values$current_plot, "tree_stretch_y")
        if (is.null(tree_stretch_x)) tree_stretch_x <- 1
        if (is.null(tree_stretch_y)) tree_stretch_y <- 1
        # v180: Get keep_proportions setting
        keep_proportions <- attr(values$current_plot, "keep_proportions")
        if (is.null(keep_proportions)) keep_proportions <- FALSE

        # v180: Calculate proportions adjustment for portrait pages
        current_aspect <- width_in / height_in
        landscape_aspect <- 29.7 / 21
        proportion_adj_w <- 1
        proportion_adj_h <- 1
        if (isTRUE(keep_proportions) && current_aspect < 1) {
          proportion_adj_w <- 1
          proportion_adj_h <- current_aspect / landscape_aspect
          debug_cat(paste0("v180: Download - preserving proportions, adj_h=", round(proportion_adj_h, 3), "\n"))
        }

        needs_transform <- (!is.null(offset_x) && !is.null(offset_y) && (offset_x != 0 || offset_y != 0)) ||
                          scale_pct != 100 || tree_stretch_x != 1 || tree_stretch_y != 1 ||
                          proportion_adj_w != 1 || proportion_adj_h != 1

        if (needs_transform) {
          x_pos_offset <- if (!is.null(offset_x)) offset_x * 0.05 else 0
          y_pos_offset <- if (!is.null(offset_y)) offset_y * 0.05 else 0
          scale_factor <- scale_pct / 100

          # v180: Apply tree stretch and proportions preservation
          final_width <- scale_factor * tree_stretch_x * proportion_adj_w
          final_height <- scale_factor * tree_stretch_y * proportion_adj_h

          # Calculate position to center the scaled plot
          center_x <- (1 - final_width) / 2 + x_pos_offset
          center_y <- (1 - final_height) / 2 + y_pos_offset

          plot_to_download <- cowplot::ggdraw() +
            cowplot::draw_plot(values$current_plot, x = center_x, y = center_y,
                              width = final_width, height = final_height)
          debug_cat(paste0("v180: Download plot with scale: ", scale_pct, "%, stretch: x=", tree_stretch_x, ", y=", tree_stretch_y, ", proportions adj: ", round(proportion_adj_h, 3), "\n"))
        }

        tryCatch({
          suppressWarnings(ggsave(
            filename = file,
            plot = plot_to_download,
            width = width_in,
            height = height_in,
            units = "in",
            dpi = dpi_val,
            device = input$output_format,
            limitsize = FALSE
          ))
          debug_cat(paste0("v134: Download saved successfully: ", file, "\n"))
        }, error = function(e) {
          debug_cat(paste0("v134: Error saving download: ", e$message, "\n"))
          # Fallback to basic plot if ggsave fails
          if (input$output_format %in% c("pdf", "svg")) {
            pdf(file, width = width_in, height = height_in)
          } else {
            png(file, width = width_val, height = height_val, units = input$output_units, res = 300)
          }
          plot(values$tree, main = "Tree Visualization")
          dev.off()
        })
      } else {
        # Fallback if no generated plot available
        debug_cat("v134: No current_plot available, using basic tree plot\n")
        if (input$output_format %in% c("pdf", "svg")) {
          pdf(file, width = input$output_width/2.54, height = input$output_height/2.54)
        } else {
          png(file, width = input$output_width, height = input$output_height, units = input$output_units, res = 300)
        }
        plot(values$tree, main = "Tree Visualization")
        dev.off()
      }
    }
  )
  
  # S1.62dev: Convert settings to YAML - simply converts the input structure to YAML text
  # The actual YAML structure is built in yaml_content() reactive
  settings_to_yaml <- function(settings) {
    yaml::as.yaml(settings, indent.mapping.sequence = TRUE)
  }

} # End of server function

# Run the application
#shinyApp(ui = ui, server = server)

#####part 6

####### part 6





#HEREEE

shiny::devmode(TRUE)
######part RUN

# Run the application
shinyApp(ui = ui, server = server)