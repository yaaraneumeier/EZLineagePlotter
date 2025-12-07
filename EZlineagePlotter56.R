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
  
  # Try exact matches first
  exact_matches <- intersect(tree_labels, csv_ids)
  matched$exact_matches <- exact_matches
  
  # Remove exact matches from consideration
  remaining_tree_labels <- setdiff(tree_labels, exact_matches)
  
  if (length(remaining_tree_labels) > 0) {
    # Try numeric conversion match
    numeric_tree <- suppressWarnings(as.numeric(remaining_tree_labels))
    numeric_csv <- suppressWarnings(as.numeric(csv_ids))
    
    valid_numeric_tree <- !is.na(numeric_tree)
    valid_numeric_csv <- !is.na(numeric_csv)
    
    if (any(valid_numeric_tree) && any(valid_numeric_csv)) {
      for (i in which(valid_numeric_tree)) {
        matches <- which(numeric_csv == numeric_tree[i])
        if (length(matches) > 0) {
          matched$numeric_matches[[remaining_tree_labels[i]]] <- csv_ids[matches]
        }
      }
    }
    
    # Try prefix/suffix matching
    for (tree_label in remaining_tree_labels) {
      if (tree_label %in% names(matched$numeric_matches)) next
      
      # Check if tree label is prefix or suffix of any csv id
      prefix_matches <- character(0)
      suffix_matches <- character(0)
      csv_as_prefix <- character(0)  # Initialize here
      csv_as_suffix <- character(0)  # Initialize here
      
      # Safely check for prefix/suffix matches
      for (csv_id in csv_ids) {
        # Skip NA values
        if (is.na(tree_label) || is.na(csv_id)) next
        
        # Check for prefix match
        if (grepl(paste0("^", tree_label), csv_id, fixed = FALSE)) {
          prefix_matches <- c(prefix_matches, csv_id)
        }
        
        # Check for suffix match
        if (grepl(paste0(tree_label, "$"), csv_id, fixed = FALSE)) {
          suffix_matches <- c(suffix_matches, csv_id)
        }
        
        # Check if csv id is prefix of tree label
        if (grepl(paste0("^", csv_id), tree_label, fixed = FALSE)) {
          csv_as_prefix <- c(csv_as_prefix, csv_id)
        }
        
        # Check if csv id is suffix of tree label
        if (grepl(paste0(csv_id, "$"), tree_label, fixed = FALSE)) {
          csv_as_suffix <- c(csv_as_suffix, csv_id)
        }
      }
      
      all_matches <- c(prefix_matches, suffix_matches, csv_as_prefix, csv_as_suffix)
      
      if (length(all_matches) > 0) {
        matched$prefix_suffix_matches[[tree_label]] <- unique(all_matches)
      } else {
        matched$unmatched <- c(matched$unmatched, tree_label)
      }
    }
  }
  
  # Generate a mapping between tree labels and csv ids
  final_mapping <- list()
  
  # Add exact matches
  for (label in matched$exact_matches) {
    final_mapping[[label]] <- label
  }
  
  # Add numeric matches (pick first if multiple)
  for (tree_label in names(matched$numeric_matches)) {
    final_mapping[[tree_label]] <- matched$numeric_matches[[tree_label]][1]
  }
  
  # Add prefix/suffix matches (pick first if multiple)
  for (tree_label in names(matched$prefix_suffix_matches)) {
    final_mapping[[tree_label]] <- matched$prefix_suffix_matches[[tree_label]][1]
  }
  
  # v53: print("final_mapping is")
  # v53: print(final_mapping)
  
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
      cat(file=stderr(), paste0("\n=== v82: Repairing corrupted plot mapping ===\n"))
      cat(file=stderr(), paste0("  Original mapping class: ", paste(class(p$mapping), collapse=", "), "\n"))
    }

    tryCatch({
      p$mapping <- aes()
      repaired <- TRUE
      if (verbose) {
        cat(file=stderr(), paste0("  Fixed mapping class: ", paste(class(p$mapping), collapse=", "), "\n"))
      }
    }, error = function(e) {
      cat(file=stderr(), paste0("  v82: Could not repair mapping: ", e$message, "\n"))
    })
  }

  # v82: REMOVED layer mapping repairs - these were breaking the heatmap fill mapping
  # The layer mappings are set correctly by gheatmap and should not be modified.
  # Only the top-level plot mapping sometimes gets corrupted to a data.frame.

  if (repaired && verbose) {
    cat(file=stderr(), paste0("================================\n"))
  }

  return(p)
}

# v71: Function to diagnose which layer is causing ggplot_build to fail
func.diagnose.layer.issues <- function(p, verbose = TRUE) {
  if (verbose) {
    cat(file=stderr(), paste0("\n=== v71: DIAGNOSING LAYER ISSUES ===\n"))
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
          cat(file=stderr(), paste0("  Layer ", i, " (", geom_class, "): FAILED - ", e$message, "\n"))
        }
        FALSE
      })

      if (!layer_ok) {
        problematic_layers <- c(problematic_layers, i)
      } else if (verbose) {
        cat(file=stderr(), paste0("  Layer ", i, " (", geom_class, "): OK\n"))
      }
    }
  }

  if (verbose) {
    if (length(problematic_layers) == 0) {
      cat(file=stderr(), paste0("  No obvious layer issues detected\n"))
    } else {
      cat(file=stderr(), paste0("  Problematic layers: ", paste(problematic_layers, collapse=", "), "\n"))
    }
    cat(file=stderr(), paste0("================================\n"))
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
        cat(file=stderr(), paste0("  v146: Using p$data for ellipse positioning (heat_flag=TRUE)\n"))
        cat(file=stderr(), paste0("  v146: high_nodes_table1 rows: ", nrow(high_nodes_table1), "\n"))
        if (nrow(high_nodes_table1) > 0) {
          cat(file=stderr(), paste0("  v146: y-range: ", round(min(high_nodes_table1$y, na.rm=TRUE), 2),
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
      cat(file=stderr(), paste0("\n=== v148: PLOT ELLIPSE (high1) ===\n"))
      cat(file=stderr(), paste0("  Alpha: ", alpha_val, "\n"))
      cat(file=stderr(), paste0("  Dimensions: a=", round(a, 4), ", b=", round(b, 4), "\n"))
      cat(file=stderr(), paste0("==================================\n"))

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

        cat(file=stderr(), paste0("\n=== v148: ELLIPSE X0 POSITIONING (heat mode) ===\n"))
        cat(file=stderr(), paste0("  tree_max_x (from p$data tips): ", round(tree_max_x_from_p, 4), "\n"))
        cat(file=stderr(), paste0("  tree_max_x (from pr440 - for reference): ", round(tree_max_x_original, 4), "\n"))
        cat(file=stderr(), paste0("  man_adjust_elipse: ", man_adjust_elipse, "\n"))
        cat(file=stderr(), paste0("  Sample node x values: ", paste(round(head(high_nodes_table1$x, 3), 4), collapse=", "), "\n"))
        cat(file=stderr(), paste0("  Sample calculated x0: ", paste(round((tree_max_x_from_p - head(high_nodes_table1$x, 3)) * (-1) - man_adjust_elipse, 4), collapse=", "), "\n"))
        cat(file=stderr(), paste0("=============================================\n"))

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

      p <- p +
        geom_ellipse(data = high_nodes_table2,
                     aes(x0 = ((max(p$data[,'x']) - x) * (x_range_min)),
                         y0 = y, a = a, b = b, angle = 0),
                     fill = high_color_list[[2]], alpha = alpha_val2, linetype = "blank", show.legend = FALSE)
    } else if (index_high == 3) {
      high_nodes_table3 <- p$data[tree_TRY$data$high3 == TRUE,]
      # v139: Use high_alpha_list for transparency
      alpha_val3 <- if (length(high_alpha_list) >= 3 && !is.null(high_alpha_list[[3]])) high_alpha_list[[3]] else 0.5

      p <- p +
        geom_ellipse(data = high_nodes_table3,
                     aes(x0 = ((max(pr440_short_tips_TRY$data[,'x']) - x) * (x_range_min)),
                         y0 = y, a = a, b = b, angle = 0),
                     fill = high_color_list[[3]], alpha = alpha_val3, linetype = "blank", show.legend = FALSE)
    }
  }
  
  # Reorder layers for proper display - BUT ONLY IF ENOUGH LAYERS EXIST
  # Bug #11 fix: Check layer count before reordering to prevent crashes
  num_layers <- length(p$layers)
  # v53: cat(file=stderr(), paste0("🔵 func_highlight: Plot has ", num_layers, " layers\n"))
  
  if (num_layers >= 8) {
    # Only do layer reordering if we have enough layers
    p1 <- p$layers[1]
    p2 <- p$layers[2]
    p7 <- p$layers[7]
    p$layers[2] <- p7 
    p$layers[7] <- p2
    
    p3 <- p$layers[3]
    p4 <- p$layers[4]
    p5 <- p$layers[5] 
    p6 <- p$layers[6]
    p7 <- p$layers[7]
    p8 <- p$layers[8]
    
    p$layers[3] <- p6
    p$layers[4] <- p7
    p$layers[5] <- p8
    p$layers[6] <- p3
    p$layers[7] <- p4  
    p$layers[8] <- p3
    
    # v53: cat(file=stderr(), "🔵 func_highlight: Layer reordering applied (8+ layers)\n")
  } else {
    # v53: cat(file=stderr(), paste0("🔵 func_highlight: Skipping layer reordering (only ", num_layers, " layers, need 8)\n"))
  }
  
  return(p)
}
# v158: Helper function to add custom legends to gtable's legend area
# This places highlight and bootstrap legends alongside ggplot legends (Classification, Heatmap)
# Uses absolute units and proper viewport to ensure visibility
func.add.custom.legends.to.gtable <- function(gt, legend_info) {
  if (is.null(legend_info)) {
    return(gt)
  }

  cat(file=stderr(), paste0("\n=== v158: ADDING CUSTOM LEGENDS TO GTABLE ===\n"))

  # v158: Find the guide-box that actually has content (based on legend position)
  # Try guide-boxes in order of priority: right, bottom, left, top
  guide_box_names <- c("guide-box-right", "guide-box-bottom", "guide-box-left", "guide-box-top", "guide-box")
  guide_box_idx <- NULL
  guide_box_name <- NULL

  for (gb_name in guide_box_names) {
    idx <- which(gt$layout$name == gb_name)
    if (length(idx) > 0) {
      # Check if this guide-box has meaningful dimensions
      gb_row <- gt$layout$t[idx]
      gb_col <- gt$layout$l[idx]

      # Get the width and height of this cell
      cell_width <- gt$widths[gb_col]
      cell_height <- gt$heights[gb_row]

      # Convert to numeric (cm) for comparison
      width_cm <- tryCatch(grid::convertWidth(cell_width, "cm", valueOnly = TRUE), error = function(e) 0)
      height_cm <- tryCatch(grid::convertHeight(cell_height, "cm", valueOnly = TRUE), error = function(e) 0)

      cat(file=stderr(), paste0("  Checking ", gb_name, " at row ", gb_row, ", col ", gb_col,
                                 " (w=", round(width_cm, 2), "cm, h=", round(height_cm, 2), "cm)\n"))

      # Use the first guide-box that has some size
      if (width_cm > 0.1 || height_cm > 0.1) {
        guide_box_idx <- idx
        guide_box_name <- gb_name
        break
      }
    }
  }

  if (is.null(guide_box_idx)) {
    cat(file=stderr(), paste0("  WARNING: No usable guide-box found, using guide-box-right as fallback\n"))
    guide_box_idx <- which(gt$layout$name == "guide-box-right")
    if (length(guide_box_idx) == 0) {
      cat(file=stderr(), paste0("  ERROR: No guide-box-right found either, skipping custom legends\n"))
      return(gt)
    }
    guide_box_name <- "guide-box-right"
  }

  guide_row <- gt$layout$t[guide_box_idx]
  guide_col <- gt$layout$l[guide_box_idx]
  cat(file=stderr(), paste0("  Using ", guide_box_name, " at row ", guide_row, ", col ", guide_col, "\n"))

  # v158: Build legend using simple grob list with fixed positioning
  # Calculate how many items we have
  n_highlight <- 0
  n_bootstrap <- 0
  if (!is.null(legend_info$highlight) && legend_info$show_highlight) {
    n_highlight <- length(legend_info$highlight$labels)
  }
  if (!is.null(legend_info$bootstrap) && legend_info$show_bootstrap) {
    n_bootstrap <- 3  # Always 3 bootstrap levels
  }

  if (n_highlight == 0 && n_bootstrap == 0) {
    cat(file=stderr(), paste0("  No custom legends to add\n"))
    return(gt)
  }

  # v158: Use fixed font sizes (ignore potentially corrupted passed values)
  title_fontsize <- 12
  text_fontsize <- 10

  # Build all grobs with absolute positioning from top
  legend_grobs <- grid::gList()
  current_y <- grid::unit(1, "npc")  # Start from top
  line_height <- grid::unit(18, "pt")

  # Highlight legend
  if (!is.null(legend_info$highlight) && legend_info$show_highlight) {
    hi <- legend_info$highlight
    cat(file=stderr(), paste0("  Highlight legend: ", n_highlight, " items\n"))

    # Title
    legend_grobs <- grid::gList(legend_grobs, grid::textGrob(
      label = hi$title,
      x = grid::unit(0.1, "cm"),
      y = current_y - grid::unit(5, "pt"),
      hjust = 0, vjust = 1,
      gp = grid::gpar(fontsize = title_fontsize, fontface = "bold")
    ))
    current_y <- current_y - line_height

    # Each highlight item
    for (i in seq_along(hi$labels)) {
      # Ellipse key
      theta <- seq(0, 2*pi, length.out = 30)
      key_size <- 6  # points

      legend_grobs <- grid::gList(legend_grobs, grid::polygonGrob(
        x = grid::unit(0.5, "cm") + grid::unit(key_size * cos(theta), "pt"),
        y = current_y - grid::unit(8, "pt") + grid::unit(key_size * 0.6 * sin(theta), "pt"),
        gp = grid::gpar(fill = hi$colors[[i]], alpha = hi$alphas[[i]], col = NA)
      ))

      # Label
      legend_grobs <- grid::gList(legend_grobs, grid::textGrob(
        label = hi$labels[[i]],
        x = grid::unit(1.2, "cm"),
        y = current_y - grid::unit(8, "pt"),
        hjust = 0, vjust = 0.5,
        gp = grid::gpar(fontsize = text_fontsize)
      ))
      current_y <- current_y - line_height
    }
    current_y <- current_y - grid::unit(10, "pt")  # Gap before bootstrap
  }

  # Bootstrap legend
  if (!is.null(legend_info$bootstrap) && legend_info$show_bootstrap) {
    cat(file=stderr(), paste0("  Bootstrap legend\n"))

    # Title
    legend_grobs <- grid::gList(legend_grobs, grid::textGrob(
      label = "Bootstrap",
      x = grid::unit(0.1, "cm"),
      y = current_y - grid::unit(5, "pt"),
      hjust = 0, vjust = 1,
      gp = grid::gpar(fontsize = title_fontsize, fontface = "bold")
    ))
    current_y <- current_y - line_height

    # Bootstrap items
    tri_labels <- c(">90%", ">80%", ">70%")
    tri_sizes <- c(5, 4, 3)  # points

    for (i in 1:3) {
      sz <- tri_sizes[i]

      # Triangle pointing up
      legend_grobs <- grid::gList(legend_grobs, grid::polygonGrob(
        x = grid::unit(0.5, "cm") + grid::unit(c(-sz, sz, 0), "pt"),
        y = current_y - grid::unit(8, "pt") + grid::unit(c(-sz*0.6, -sz*0.6, sz*0.8), "pt"),
        gp = grid::gpar(fill = "grey36", col = "grey20", alpha = 0.5)
      ))

      # Label
      legend_grobs <- grid::gList(legend_grobs, grid::textGrob(
        label = tri_labels[i],
        x = grid::unit(1.2, "cm"),
        y = current_y - grid::unit(8, "pt"),
        hjust = 0, vjust = 0.5,
        gp = grid::gpar(fontsize = text_fontsize)
      ))
      current_y <- current_y - line_height
    }
  }

  # Calculate total size (approximate)
  total_items <- n_highlight + n_bootstrap + (if(n_highlight > 0) 1 else 0) + (if(n_bootstrap > 0) 1 else 0)
  total_height_cm <- total_items * 0.5 + 0.5
  total_width_cm <- 3.0

  cat(file=stderr(), paste0("  Total items: ", total_items, ", size: ", total_height_cm, " x ", total_width_cm, " cm\n"))

  # Create a gTree with all the legend grobs
  custom_legend <- grid::gTree(children = legend_grobs, vp = grid::viewport())

  # v158: Add to appropriate position based on guide-box type
  is_horizontal <- grepl("bottom|top", guide_box_name)

  if (is_horizontal) {
    # For bottom/top legends, add a column to the LEFT of the guide-box
    cat(file=stderr(), paste0("  Adding column for horizontal legend layout\n"))
    gt <- gtable::gtable_add_cols(gt, widths = grid::unit(total_width_cm, "cm"), pos = guide_col - 1)

    # Add grob to the new column
    gt <- gtable::gtable_add_grob(
      gt,
      grobs = custom_legend,
      t = guide_row,
      l = guide_col,  # New column is now at this position
      name = "custom-highlight-bootstrap-legend"
    )
    cat(file=stderr(), paste0("  Custom legends added LEFT of guide-box at row ", guide_row, ", col ", guide_col, "\n"))
  } else {
    # For right/left legends, add a row ABOVE the guide-box
    cat(file=stderr(), paste0("  Adding row for vertical legend layout\n"))
    gt <- gtable::gtable_add_rows(gt, heights = grid::unit(total_height_cm, "cm"), pos = guide_row - 1)

    # Add grob to the new row
    gt <- gtable::gtable_add_grob(
      gt,
      grobs = custom_legend,
      t = guide_row,  # New row is now at this position
      l = guide_col,
      name = "custom-highlight-bootstrap-legend"
    )
    cat(file=stderr(), paste0("  Custom legends added ABOVE guide-box at row ", guide_row, ", col ", guide_col, "\n"))
  }

  cat(file=stderr(), paste0("===============================================\n"))

  return(gt)
}


# Function to create the second legend
# v145: Added high_alpha_list parameter for transparency control
# v158: Modified to store legend info as attribute instead of using annotation_custom
func.make.second.legend <- function(p, FLAG_BULK_DISPLAY, how_many_hi, heat_flag, how_many_boxes,
                                    how_mant_rows, boudariestt, y_off_base, high_title_list,
                                    size_font_legend_title, high_label_list, size_font_legend_text,
                                    high_color_list, a, b, x_range_min, show_boot_flag, size_90,
                                    size_80, size_70, man_adjust_image_of_second_legend,
                                    man_multiply_second_legend, man_multiply_second_legend_text,
                                    man_multiply_elipse, man_space_second_legend,
                                    man_space_second_legend_multiplier, man_offset_for_highlight_legend_x,
                                    debug_mode = FALSE, boot_values = NA, man_offset_second_legend = 0, width,
                                    bootstrap_label_size = 1.5,  # v129: Reduced from 3.5 for smaller default legend
                                    # v133: New highlight legend settings
                                    highlight_x_offset = 0, highlight_y_offset = 0,
                                    highlight_title_size = NULL, highlight_text_size = NULL,
                                    highlight_title_gap = 1, highlight_label_gap = 0.5,
                                    # v133: New bootstrap legend settings
                                    # v143: Added bootstrap_title_x_offset for moving title to the right
                                    bootstrap_x_offset = 0, bootstrap_y_offset = 0,
                                    bootstrap_title_x_offset = 2,  # v143: Extra offset for title (moves it right)
                                    bootstrap_title_size_mult = NULL, bootstrap_text_size_mult = NULL,
                                    bootstrap_title_gap = 2, bootstrap_label_gap = 2,
                                    # v138: Show/hide legend controls
                                    show_highlight_legend = TRUE, show_bootstrap_legend = TRUE,
                                    # v145: Transparency list for legend ellipses
                                    high_alpha_list = NULL) {

  # v158: GTABLE-BASED APPROACH
  # Store legend specifications as an attribute on the plot object.
  # The rendering code will add these to the gtable's legend area (alongside ggplot legends).

  cat(file=stderr(), paste0("\n=== v158: GTABLE-BASED LEGENDS ===\n"))
  cat(file=stderr(), paste0("  Storing legend info as attribute for gtable insertion\n"))
  cat(file=stderr(), paste0("  Legends will appear in same area as ggplot legends (Classification, Heatmap)\n"))

  # Initialize high_alpha_list if NULL
  if (is.null(high_alpha_list) || length(high_alpha_list) == 0) {
    high_alpha_list <- rep(0.5, how_many_hi)
  }

  # Title fontsize should match ggplot legend titles
  title_fontsize <- if (!is.null(highlight_title_size)) highlight_title_size else size_font_legend_title
  text_fontsize <- if (!is.null(highlight_text_size)) highlight_text_size else size_font_legend_text

  # Bootstrap title size
  boot_title_fontsize <- if (!is.null(bootstrap_title_size_mult)) bootstrap_title_size_mult else size_font_legend_title
  boot_text_fontsize <- if (!is.null(bootstrap_text_size_mult)) bootstrap_text_size_mult else size_font_legend_text

  cat(file=stderr(), paste0("  v158: Title fontsize: ", title_fontsize, " (matches ggplot legend.title)\n"))
  cat(file=stderr(), paste0("  v158: Text fontsize: ", text_fontsize, "\n"))

  # Build legend info structure
  legend_info <- list()

  # Highlight legend info
  if (FLAG_BULK_DISPLAY == TRUE && show_highlight_legend == TRUE && how_many_hi > 0) {
    # Get first highlight title (or use "Highlight" if not available)
    highlight_title <- if (!is.null(high_title_list) && length(high_title_list) > 0) {
      high_title_list[[1]]
    } else {
      "Highlight"
    }

    legend_info$highlight <- list(
      title = highlight_title,
      labels = high_label_list,
      colors = high_color_list,
      alphas = high_alpha_list,
      title_fontsize = title_fontsize,
      text_fontsize = text_fontsize
    )

    cat(file=stderr(), paste0("  v158: Highlight legend: title='", highlight_title,
                               "', ", length(high_label_list), " items\n"))
    for (i in seq_along(high_label_list)) {
      cat(file=stderr(), paste0("    Item ", i, ": label='", high_label_list[[i]],
                                 "', color=", high_color_list[[i]],
                                 ", alpha=", high_alpha_list[[i]], "\n"))
    }
  }

  # Bootstrap legend info
  if (show_boot_flag == TRUE && show_bootstrap_legend == TRUE) {
    # boot_values is a list, so check if it's a valid list with 'format' element
    if (!is.null(boot_values) && is.list(boot_values) && !is.null(boot_values$'format') && boot_values$'format' == 'triangles') {
      legend_info$bootstrap <- list(
        title_fontsize = boot_title_fontsize,
        text_fontsize = boot_text_fontsize
      )
      cat(file=stderr(), paste0("  v158: Bootstrap legend enabled\n"))
    }
  }

  legend_info$show_highlight <- show_highlight_legend
  legend_info$show_bootstrap <- show_bootstrap_legend

  # Store legend info as attribute on the plot
  attr(p, "custom_legend_info") <- legend_info

  cat(file=stderr(), paste0("  v158: Legend info stored as plot attribute\n"))
  cat(file=stderr(), paste0("  v158: Will be added to gtable during rendering\n"))
  cat(file=stderr(), paste0("=================================================\n"))

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
  # v53: cat(file=stderr(), "==========================================\n")
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
  # v53: cat(file=stderr(), "==========================================\n\n")
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
  cat(file=stderr(), paste0("\n=== v146: EXTRACTING HIGHLIGHT TRANSPARENCY ===\n"))
  high_alpha_list <- c()
  for (in_hi in indexes_hi) {
    alpha_val <- hi_def$according[[in_hi]][[as.character(in_hi)]]$transparency
    cat(file=stderr(), paste0("  Highlight ", in_hi, ":\n"))
    cat(file=stderr(), paste0("    hi_def$according[[", in_hi, "]][[\"", in_hi, "\"]]$transparency = ",
                              if(is.null(alpha_val)) "NULL" else alpha_val, "\n"))
    # Default to 0.5 if transparency not specified
    high_alpha_list[[in_hi]] <- if (!is.null(alpha_val)) alpha_val else 0.5
    cat(file=stderr(), paste0("    high_alpha_list[[", in_hi, "]] = ", high_alpha_list[[in_hi]], "\n"))
  }
  cat(file=stderr(), paste0("  Final high_alpha_list: ", paste(high_alpha_list, collapse=", "), "\n"))
  cat(file=stderr(), paste0("==============================================\n"))

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
        flag_need_to_fix <- TRUE
        ix <- 1
        
        children_weights <- func.make.children.weight.list(children, nod, list_weights_for_nodes_dx.rx)
        children_weights_SECOND <- func.make.children.weight.list(children, nod, list_weights_for_nodes_frac)
        
        children_weights_ordered <- sort(children_weights)
        
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
        
        for (i in 1:length(children)) {
          if (children[i] == children[dest[i]]) {
            # Don't flip
          } else {
            tree_return <- flip(tree_return, children[i], children[dest[i]])
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
  # Check if list_nodes_to_rotate is valid (not NA and has length > 0)
  if (!is.null(list_nodes_to_rotate) && length(list_nodes_to_rotate) > 0 && !all(is.na(list_nodes_to_rotate))) {
    # v53: print("rotate specific nodes")
    # v53: print(list_nodes_to_rotate)
    for (nod in list_nodes_to_rotate) {
      children <- which(tree_TRY1$data$parent == nod & tree_TRY1$data$node != nod)
      # v53: print(children)
      
      if (length(children) < 2) {
        tree_TRY1 <- flip(tree_TRY1, children[1], children[3]) # 321
        tree_TRY1 <- flip(tree_TRY1, children[2], children[3])
      } else {
        tree_TRY1 <- flip(tree_TRY1, children[1], children[2])
      }
    }
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
                                    legend_settings = NULL) {  # v135: Legend settings for highlight/bootstrap legends

  # === DEBUG CHECKPOINT 2: FUNCTION ENTRY ===
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 2: func.print.lineage.tree ENTRY\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â node_number_font_size received:", node_number_font_size, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â flag_display_nod_number_on_tree:", flag_display_nod_number_on_tree, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â highlight_manual_nodes received:", highlight_manual_nodes, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â manual_nodes_to_highlight received:", paste(manual_nodes_to_highlight, collapse=", "), "\n")
  # v53: cat(file=stderr(), "================================================\n\n")
  
  yaml_file<- func.read_yaml(conf_yaml_path)
  
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
  if (flag_csv_read_func=="fread"){
    fread_rownames(csv_path, row.var = rowname_param)
  } else {
    readfile <- read.csv(csv_path)
  }
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
    tree440 <- read.tree(tree_path)
    
    # Apply ladderize if flag is TRUE
    if (laderize_flag == TRUE) {
      # v53: print("Ladderizing tree...")
      tree440 <- ape::ladderize(tree440, right = TRUE)
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
      cat(file=stderr(), paste0("\n=== v57: Classification attributes (att1): ", paste(att1, collapse=", "), " ===\n"))

      if ('heatmap_display' %in% att1) {

        # v56c: DEBUG
        cat(file=stderr(), "\n=== v57: FOUND heatmap_display in classification ===\n")

        heat_definitions <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$heatmap_display
        heat_list_len <- length(heat_definitions)
        cat(file=stderr(), paste0("  Number of heatmaps found: ", heat_list_len, "\n"))
        cat(file=stderr(), paste0("  heat_definitions structure: ", class(heat_definitions), "\n"))
        if (heat_list_len > 0) {
          cat(file=stderr(), paste0("  First heatmap names: ", paste(names(heat_definitions[[1]]), collapse=", "), "\n"))
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

              # v109: Get colnames_angle (default 45)
              if ('colnames_angle' %in% names(heat_map_i_def)) {
                param[['colnames_angle']] <- as.numeric(heat_map_i_def[['colnames_angle']])
              } else {
                param[['colnames_angle']] <- 45
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
              cat(file=stderr(), paste0("\n=== v121: TIP GUIDE SETTINGS DEBUG (discrete) ===\n"))
              cat(file=stderr(), paste0("  heat_map_i_def names: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
              if ('show_guides' %in% names(heat_map_i_def)) {
                raw_val <- heat_map_i_def[['show_guides']]
                cat(file=stderr(), paste0("  show_guides raw value: ", raw_val, " (class: ", class(raw_val), ")\n"))
                param[['show_guides']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_guides']])
                cat(file=stderr(), paste0("  show_guides after conversion: ", param[['show_guides']], "\n"))
              } else {
                cat(file=stderr(), paste0("  show_guides NOT FOUND in heat_map_i_def\n"))
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

              # v109: Get colnames_angle (default 45) for continuous heatmaps too
              if ('colnames_angle' %in% names(heat_map_i_def)) {
                param[['colnames_angle']] <- as.numeric(heat_map_i_def[['colnames_angle']])
              } else {
                param[['colnames_angle']] <- 45
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
              cat(file=stderr(), paste0("\n=== v121: TIP GUIDE SETTINGS DEBUG (continuous) ===\n"))
              cat(file=stderr(), paste0("  heat_map_i_def names: ", paste(names(heat_map_i_def), collapse=", "), "\n"))
              if ('show_guides' %in% names(heat_map_i_def)) {
                raw_val <- heat_map_i_def[['show_guides']]
                cat(file=stderr(), paste0("  show_guides raw value: ", raw_val, " (class: ", class(raw_val), ")\n"))
                param[['show_guides']] <- func.check.bin.val.from.conf(heat_map_i_def[['show_guides']])
                cat(file=stderr(), paste0("  show_guides after conversion: ", param[['show_guides']], "\n"))
              } else {
                cat(file=stderr(), paste0("  show_guides NOT FOUND in heat_map_i_def\n"))
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
            }



            heat_display_params_list[[indx_for_sav]] <- param
            # print("B12")
            
            
            
            heat_display_vec <- c(heat_display_vec, TRUE)
            heat_flag <- TRUE
            acc_heat_list <- heat_map_i_def$according
            heat_map_title <- heat_map_i_def$title
            # v53: print("heat_map_title is")
            # v53: print(heat_map_title)
            heat_map_title_list <- c(heat_map_title_list,heat_map_title)
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
              cat(file=stderr(), paste0("\n=== v57: Extracting heatmap columns from 'according' ===\n"))
              cat(file=stderr(), paste0("  acc_heat_list length: ", length(acc_heat_list), "\n"))
              ind <-1
              for (j in acc_heat_list) {

                j1 <- names(j)
                cat(file=stderr(), paste0("  Column ", ind, ": j1=", j1, "\n"))

                ind<- ind+1
                j2<- j[[j1]]
                cat(file=stderr(), paste0("    j2 (column name)=", j2, "\n"))

                l_titles_for_heat <- c(l_titles_for_heat,j2)
              }
              cat(file=stderr(), paste0("  Final l_titles_for_heat: ", paste(l_titles_for_heat, collapse=", "), "\n"))
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
            cat(file=stderr(), paste0("\n=== v57: Validating heatmap columns ===\n"))
            cat(file=stderr(), paste0("  title.id: ", title.id, "\n"))
            cat(file=stderr(), paste0("  l_titles_for_heat: ", paste(l_titles_for_heat, collapse=", "), "\n"))
            cat(file=stderr(), paste0("  Available CSV columns: ", paste(head(names(readfile440), 10), collapse=", "), "...\n"))

            valid_columns <- c(title.id, l_titles_for_heat)
            valid_columns <- valid_columns[valid_columns %in% names(readfile440)]
            cat(file=stderr(), paste0("  Valid columns (after filtering): ", paste(valid_columns, collapse=", "), "\n"))
            cat(file=stderr(), paste0("  Number of valid columns: ", length(valid_columns), "\n"))

            # Select only the valid columns
            df_heat_temp <- readfile440[, valid_columns, drop = FALSE]
            cat(file=stderr(), paste0("  df_heat_temp dimensions: ", nrow(df_heat_temp), " rows x ", ncol(df_heat_temp), " cols\n"))
            
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
            cat(file=stderr(), paste0("\n=== v66: TIP LABEL EXTRACTION DEBUG ===\n"))
            cat(file=stderr(), paste0("  ggtree_labels NA count: ", sum(is.na(ggtree_labels)), " out of ", length(ggtree_labels), "\n"))
            cat(file=stderr(), paste0("  tree440$tip.label sample: ", paste(head(tree440$tip.label, 5), collapse=", "), "\n"))
            cat(file=stderr(), paste0("  g_check_tip$node sample: ", paste(head(g_check_tip$node, 5), collapse=", "), "\n"))

            # v66: Check if ANY labels are NA or empty - if so, use tree440$tip.label
            # This is more robust than only checking if ALL are NA
            if (any(is.na(ggtree_labels)) || any(ggtree_labels == "")) {
              cat(file=stderr(), paste0("  v66: Found NA/empty labels in ggtree, using tree440$tip.label\n"))

              # For tip nodes, node IDs 1 to Ntip correspond directly to tree$tip.label indices
              tip_node_ids <- g_check_tip$node
              ntips <- length(tree440$tip.label)

              # Validate that node IDs are in valid range
              if (all(tip_node_ids >= 1 & tip_node_ids <= ntips)) {
                ggtree_labels <- tree440$tip.label[tip_node_ids]
                cat(file=stderr(), paste0("  v66: Successfully extracted labels from tree440$tip.label\n"))
              } else {
                # Fallback: use tip labels in their original order from tree440
                cat(file=stderr(), paste0("  v66: WARNING - node IDs out of range, using tree tip order\n"))
                # Order g_check_tip by y coordinate (visual order) and assign labels
                tip_order <- order(g_check_tip$y)
                ggtree_labels <- tree440$tip.label[tip_order]
              }

              cat(file=stderr(), paste0("  v66: After fix, ggtree_labels sample: ", paste(head(ggtree_labels, 5), collapse=", "), "\n"))
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
            
            
            # v58: DEBUG - Show matching details
            cat(file=stderr(), paste0("\n=== v58: Matching heatmap data to tree tips ===\n"))
            cat(file=stderr(), paste0("  tip_list length: ", length(tip_list), "\n"))
            cat(file=stderr(), paste0("  tip_list sample: ", paste(head(tip_list, 5), collapse=", "), "\n"))
            cat(file=stderr(), paste0("  CSV ID column (", title.id, ") sample: ", paste(head(df_heat_temp[[title.id]], 5), collapse=", "), "\n"))

            # Check for matches before applying
            matches <- match(tip_list, df_heat_temp[[title.id]])
            num_matches <- sum(!is.na(matches))
            cat(file=stderr(), paste0("  Number of matches found: ", num_matches, " out of ", length(tip_list), " tips\n"))

            df_heat_temp <- df_heat_temp[matches,]
            cat(file=stderr(), paste0("  After match(): ", nrow(df_heat_temp), " rows\n"))

            ro= na.omit(df_heat_temp[[title.id]])
            df_heat_temp_filtered<- df_heat_temp[df_heat_temp[[title.id]] %in% (ro), ]

            df_heat_temp<- df_heat_temp_filtered
            cat(file=stderr(), paste0("  After filtering NAs: ", nrow(df_heat_temp), " rows\n"))
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
              cat(file=stderr(), paste0("  WARNING: No matching data for heatmap - skipping\n"))
              heat_display_vec <- c(heat_display_vec, FALSE)
              # v58: Reset heat_flag if no heatmaps have data
              if (indx_for_sav == 1) {
                heat_flag <- FALSE
                cat(file=stderr(), paste0("  Resetting heat_flag to FALSE\n"))
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
      # v53: cat(file=stderr(), "================================================\n\n")
      
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
        legend_settings = legend_settings  # v136: Pass legend settings for highlight/bootstrap legends
      )
      # }

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
                                         legend_settings = NULL) {  # v136: Legend settings for highlight/bootstrap legends

  # === DEBUG CHECKPOINT 4: INNER FUNCTION ENTRY ===
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 4: func.make.plot.tree.heat.NEW ENTRY\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â node_number_font_size:", node_number_font_size, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â flag_display_nod_number_on_tree:", flag_display_nod_number_on_tree, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â highlight_manual_nodes:", highlight_manual_nodes, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â manual_nodes_to_highlight:", paste(manual_nodes_to_highlight, collapse=", "), "\n")
  # v53: cat(file=stderr(), "================================================\n\n")
  
  if (debug_mode == TRUE) {
    # v53: print("In func.make.plot.tree.HEAT")
  }
  
  # Prepare data structures
  dx_rx_types1 <- dx_rx_types1_short
  # v53: print("A")
  # v53: print(dx_rx_types1)
  
  # v56b: Wrap in suppressWarnings to suppress harmless ggtree/ggplot2 fortify warnings
  pr440 <- suppressWarnings(ggtree(tree440))
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
  tree_with_group <- ggtree::groupOTU(tree440, list_node_by_class)
  
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
  tree_with_group_CPY <- ggtree::groupOTU(tree440, list_rename_by_class)
  # v56b: Wrap in suppressWarnings to suppress harmless fortify warnings
  levels_groups <- levels(suppressWarnings(ggtree(tree_with_group_CPY))$data$group)
  # v53: print("E")

  # Create tree with coloring
  # v56b: Wrap in suppressWarnings to suppress harmless fortify warnings
  tree_TRY <- suppressWarnings(ggtree(tree_with_group_CPY, aes(color = new_class, size = p_val_new),
                     ladderize = laderize_flag))
  
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
  
  # Calculate p-values
  p_list_of_pairs <- func.create.p_list_of_pairs(
    list_node_by_class, d440, dx_rx_types1_short,
    cc_nodss, tree_with_group, FDR_perc, tree, cc_tipss,
    tree_TRY, tree_size, no_name, simulate.p.value
  )
  
  p_PAIRS_pval_list <- func.create.p_val_list_FROM_LIST(FDR_perc, tree_TRY, p_list_of_pairs, op_list)
  
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
  pr440_short_tips_TRY_new <- pr440_short_tips_TRY + 
    scale_color_manual(values = colors_scale2) + 
    scale_size_manual(values = list_of_sizes)
  
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
  
  # Add node numbers if requested
  # v53: cat(file=stderr(), "\nÃ°Å¸â€Â DEBUG CHECKPOINT 5: NODE NUMBER RENDERING\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â flag_display_nod_number_on_tree is:", flag_display_nod_number_on_tree, "\n")
  # v53: cat(file=stderr(), "Ã°Å¸â€Â node_number_font_size is:", node_number_font_size, "\n")
  
  if (flag_display_nod_number_on_tree == TRUE) {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€œ ADDING NODE NUMBERS with size:", node_number_font_size, "\n")
    # v53: cat(file=stderr(), "================================================\n\n")
    
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
      geom_text(
        aes(label = node, angle = 90, colour = "black"), 
        hjust = -0.5, vjust = -0.4, size = node_number_font_size,
        show.legend = FALSE, colour = "black"
      )
  } else {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€” NODE NUMBERS NOT ENABLED\n")
    # v53: cat(file=stderr(), "================================================\n\n")
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
      # v53: cat(file=stderr(), "================================================\n\n")
      
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
      # v53: cat(file=stderr(), "================================================\n\n")
    }
  } else {
    # v53: cat(file=stderr(), "Ã°Å¸â€Â Ã¢Å“â€” HIGHLIGHTING NOT ENABLED\n")
    if (highlight_manual_nodes != TRUE) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â   Reason: highlight_manual_nodes is FALSE\n")
    }
    if (is.na(manual_nodes_to_highlight[1])) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Â   Reason: manual_nodes_to_highlight is NA\n")
    }
    # v53: cat(file=stderr(), "================================================\n\n")
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
    cat(file=stderr(), paste0("  v132: Error computing boudariestt: ", e$message, "\n"))
    # Return default values if computation fails
    list(xmin = min(p$data$x, na.rm = TRUE), xmax = max(p$data$x, na.rm = TRUE))
  })
  cat(file=stderr(), paste0("  v132: boudariestt initialized early: xmin=", boudariestt$xmin, ", xmax=", boudariestt$xmax, "\n"))

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
    p <- func_highlight(
      p, how_many_hi, heat_flag, high_color_list, a, b, man_adjust_elipse,
      pr440_short_tips_TRY, boudariestt, debug_mode, high_offset, high_vertical_offset,
      high_alpha_list
    )
  }

  if (length(b) == 0) {
    b <- 0.2
  }
  
  # Add heatmap if requested
  # v99: MANUAL HEATMAP - Replaced gheatmap() with manual geom_tile() approach
  # because gheatmap was corrupting the plot's @mapping property
  cat(file=stderr(), paste0("\n=== v99: HEATMAP RENDERING ===\n"))
  cat(file=stderr(), paste0("  heat_flag: ", heat_flag, "\n"))
  if (heat_flag == TRUE) {
    cat(file=stderr(), paste0("  heat_map_title_list length: ", length(heat_map_title_list), "\n"))
    cat(file=stderr(), paste0("  dxdf440_for_heat length: ", length(dxdf440_for_heat), "\n"))
  }

  # v99: Save a backup of p before heatmap for fallback recovery
  tt <- p

  # v132: Refresh boudariestt after highlighting changes (was already initialized at v132 block above)
  # This ensures func.make.second.legend has current plot boundaries
  boudariestt <- tryCatch({
    func.find.plot.boundaries(tt, debug_mode)
  }, error = function(e) {
    cat(file=stderr(), paste0("  v132: Error refreshing boudariestt: ", e$message, "\n"))
    # Return default values if computation fails
    list(xmin = min(p$data$x, na.rm = TRUE), xmax = max(p$data$x, na.rm = TRUE))
  })
  cat(file=stderr(), paste0("  v132: boudariestt refreshed: xmin=", boudariestt$xmin, ", xmax=", boudariestt$xmax, "\n"))

  # v122: MULTIPLE HEATMAPS IMPLEMENTATION
  # Refactored from v99 to support multiple heatmaps with spacing control
  # Each heatmap can have its own colors, type (discrete/continuous), and parameters
  if (heat_flag == TRUE && length(dxdf440_for_heat) > 0) {
    cat(file=stderr(), paste0("\n=== v122: ENTERING MULTIPLE HEATMAPS CODE ===\n"))
    cat(file=stderr(), paste0("  Number of heatmaps to render: ", length(dxdf440_for_heat), "\n"))

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
    cat(file=stderr(), paste0("  v125: heatmap_global_gap=", heatmap_global_gap, ", heatmap_spacing=", heatmap_spacing, "\n"))

    # v122: Loop through all heatmaps
    for (heat_idx in 1:length(dxdf440_for_heat)) {
      cat(file=stderr(), paste0("\n=== v122: PROCESSING HEATMAP ", heat_idx, " of ", length(dxdf440_for_heat), " ===\n"))

      # Get this heatmap's data
      heat_data <- dxdf440_for_heat[[heat_idx]]

      # v122: Validate heatmap data
      cat(file=stderr(), paste0("  Initial heat_data dimensions: ", nrow(heat_data), " x ", ncol(heat_data), "\n"))

      # v122: Check if heat_data is valid - skip this heatmap if invalid
      if (is.null(heat_data) || !is.data.frame(heat_data) || nrow(heat_data) == 0) {
        cat(file=stderr(), paste0("  ERROR: Invalid heatmap data - skipping heatmap ", heat_idx, "\n"))
        next  # v122: Continue to next heatmap instead of stopping all heatmaps
      }

      cat(file=stderr(), paste0("  Tree has ", length(tree_tips), " tips\n"))

      # v122: Check current row names
      current_rownames <- rownames(heat_data)
      cat(file=stderr(), paste0("  Current rownames: ", paste(head(current_rownames, 5), collapse=", "), "\n"))

      # v122: Verify row names are valid and match tree tips
      if (is.null(current_rownames) || all(current_rownames == as.character(1:nrow(heat_data)))) {
        cat(file=stderr(), paste0("  WARNING: Heat data has default numeric row names - skipping heatmap ", heat_idx, "\n"))
        next
      }

      # v122: Check how many row names match tree tips
      matching_tips <- sum(current_rownames %in% tree_tips)
      cat(file=stderr(), paste0("  Row names matching tree tips: ", matching_tips, " / ", nrow(heat_data), "\n"))

      if (matching_tips == 0) {
        cat(file=stderr(), paste0("  ERROR: No row names match tree tips - skipping heatmap ", heat_idx, "\n"))
        next
      } else if (matching_tips < nrow(heat_data)) {
        cat(file=stderr(), paste0("  WARNING: Only ", matching_tips, " row names match - filtering data\n"))
        heat_data <- heat_data[current_rownames %in% tree_tips, , drop = FALSE]
        cat(file=stderr(), paste0("  After filtering: ", nrow(heat_data), " rows\n"))
      }

      cat(file=stderr(), paste0("  heat_data dimensions: ", nrow(heat_data), " x ", ncol(heat_data), "\n"))
      cat(file=stderr(), paste0("  heat_data rownames sample: ", paste(head(rownames(heat_data), 5), collapse=", "), "\n"))
      cat(file=stderr(), paste0("  heat_data columns: ", paste(colnames(heat_data), collapse=", "), "\n"))

      # v122: Get heatmap parameters for THIS heatmap (use heat_idx, not 1)
      heat_param <- if (heat_idx <= length(heat_display_params_list)) heat_display_params_list[[heat_idx]] else list()
      is_discrete <- ifelse(!is.null(heat_param) && !is.na(heat_param['is_discrete']),
                            heat_param['is_discrete'] == TRUE, FALSE)
      cat(file=stderr(), paste0("  is_discrete: ", is_discrete, "\n"))

      # v122: Calculate heatmap positioning - use current_heatmap_x_start from previous heatmap
      # For first heatmap, this is tree_xmax; for subsequent, it's the end of the previous heatmap
      per_heatmap_distance <- if (!is.null(heat_param[['distance']])) heat_param[['distance']] else heatmap_tree_distance
      heatmap_offset <- tree_width * per_heatmap_distance
      cat(file=stderr(), paste0("  per_heatmap_distance: ", per_heatmap_distance, "\n"))
      cat(file=stderr(), paste0("  current_heatmap_x_start: ", current_heatmap_x_start, "\n"))
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
        cat(file=stderr(), paste0("  v116: tip spacing stats: min=", base_tip_spacing,
                                   ", median=", median_spacing,
                                   ", mean=", mean_spacing, "\n"))
      } else {
        base_tip_spacing <- 1.0  # Fallback for single tip
      }

      # tile_height: controls row height (data y units)
      # row_height_value=1.0 means tiles touch at minimum spacing (no overlap)
      # <1 means gaps, >1 means overlap
      tile_height <- base_tip_spacing * row_height_value

      cat(file=stderr(), paste0("  v116: tip_y range: [", min(tip_y_positions), ", ", max(tip_y_positions), "]\n"))
      cat(file=stderr(), paste0("  v116: n_tips: ", length(tip_y_positions), "\n"))
      cat(file=stderr(), paste0("  v116: base_tip_spacing (min): ", base_tip_spacing, "\n"))

      cat(file=stderr(), paste0("  column_width_value: ", column_width_value, "\n"))
      cat(file=stderr(), paste0("  row_height_value: ", row_height_value, "\n"))
      cat(file=stderr(), paste0("  tile_width (column spacing): ", tile_width, "\n"))
      cat(file=stderr(), paste0("  tile_height (row height): ", tile_height, "\n"))

      cat(file=stderr(), paste0("  tree_xmin: ", tree_xmin, "\n"))
      cat(file=stderr(), paste0("  tree_xmax: ", tree_xmax, "\n"))
      cat(file=stderr(), paste0("  tree_width: ", tree_width, "\n"))
      cat(file=stderr(), paste0("  heatmap_offset: ", heatmap_offset, "\n"))
      cat(file=stderr(), paste0("  tile_width: ", tile_width, "\n"))

      # v99: Build heatmap data frame for geom_tile
      # We need: x (column position), y (tip position), fill (value)
      cat(file=stderr(), paste0("\n=== v99: BUILDING HEATMAP TILE DATA ===\n"))

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
        cat(file=stderr(), paste0("  Created tile_df with ", nrow(tile_df), " tiles\n"))
        cat(file=stderr(), paste0("  x range: [", min(tile_df$x), ", ", max(tile_df$x), "]\n"))
        cat(file=stderr(), paste0("  y range: [", min(tile_df$y), ", ", max(tile_df$y), "]\n"))
        cat(file=stderr(), paste0("  Unique values: ", paste(unique(tile_df$value), collapse=", "), "\n"))

        # v99: Add heatmap tiles using geom_tile
        cat(file=stderr(), paste0("\n=== v99: ADDING GEOM_TILE LAYER ===\n"))

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

          cat(file=stderr(), paste0("  v113: show_grid=", show_grid, " (class: ", class(show_grid), ")\n"))
          cat(file=stderr(), paste0("  v113: grid_color=", grid_color, ", grid_size=", grid_size, "\n"))

          # v113: Ensure show_grid is properly evaluated as boolean
          show_grid_bool <- isTRUE(show_grid) || identical(show_grid, TRUE) || identical(show_grid, "yes") || identical(show_grid, "TRUE")
          cat(file=stderr(), paste0("  v113: show_grid_bool=", show_grid_bool, "\n"))

          # v123: For heatmaps after the first, add new_scale_fill() BEFORE adding geom_tile
          # This is critical - the scale must be reset before the new layer that uses fill
          if (heat_idx > 1) {
            cat(file=stderr(), paste0("  v123: Adding new_scale_fill() BEFORE geom_tile for heatmap ", heat_idx, "\n"))
            p <- p + ggnewscale::new_scale_fill()
          }

          if (show_grid_bool) {
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

            cat(file=stderr(), paste0("  v115: Added explicit grid lines: ", nrow(v_lines_df), " vertical, ", nrow(h_lines_df), " horizontal\n"))
          } else {
            p_with_tiles <- p + geom_tile(
              data = tile_df,
              aes(x = x, y = y, fill = value),
              width = tile_width,     # v112: Column spacing (fixed)
              height = tile_height,   # v112: Row height (slider-controlled)
              inherit.aes = FALSE
            )
          }

          cat(file=stderr(), paste0("  geom_tile added successfully\n"))

          # v100: Add color scale using user-selected colors
          na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "grey90"
          # v122: Use heat_idx to get the correct heatmap title
          heatmap_title <- if (heat_idx <= length(heat_map_title_list)) heat_map_title_list[[heat_idx]] else paste0("Heatmap ", heat_idx)

          # v123: new_scale_fill() is now added BEFORE geom_tile (see line ~4993)

          if (is_discrete) {
            cat(file=stderr(), paste0("  Adding discrete color scale\n"))

            # v100: Check for custom colors from user
            man_define_colors <- !is.null(heat_param['man_define_colors']) &&
                                 !is.na(heat_param['man_define_colors']) &&
                                 heat_param['man_define_colors'] == TRUE
            custom_colors_raw <- heat_param[['color_scale_option']]

            # v104: Convert from list to named vector if needed (list preserves names through YAML)
            if (is.list(custom_colors_raw) && !is.null(names(custom_colors_raw))) {
              # Convert named list to named character vector
              custom_colors <- unlist(custom_colors_raw)
              cat(file=stderr(), paste0("  v104: Converted custom_colors from named list to named vector\n"))
            } else {
              custom_colors <- custom_colors_raw
            }

            cat(file=stderr(), paste0("  man_define_colors: ", man_define_colors, "\n"))
            cat(file=stderr(), paste0("  custom_colors: ", paste(custom_colors, collapse=", "), "\n"))
            cat(file=stderr(), paste0("  na_color: ", na_color, "\n"))

            if (man_define_colors && !is.null(custom_colors) && length(custom_colors) > 0) {
              # v101: Use custom colors provided by user
              cat(file=stderr(), paste0("  Using ", length(custom_colors), " custom colors from user\n"))

              # v104: custom_colors should now be a named vector where names are the value labels
              # Debug: show what we received
              cat(file=stderr(), paste0("  custom_colors names: ", paste(names(custom_colors), collapse=", "), "\n"))
              cat(file=stderr(), paste0("  custom_colors values: ", paste(custom_colors, collapse=", "), "\n"))

              # Get unique values from data (excluding NA) to match with colors
              unique_vals <- unique(tile_df$value)
              unique_vals <- unique_vals[!is.na(unique_vals)]
              cat(file=stderr(), paste0("  Unique values in data: ", paste(unique_vals, collapse=", "), "\n"))

              # v101: Fix color mapping - custom_colors is already a named vector
              # Match colors by value name, not by position
              if (!is.null(names(custom_colors)) && length(names(custom_colors)) > 0) {
                # Use the named vector directly - lookup by value name
                color_vec <- custom_colors[unique_vals]
                # For any values not in custom_colors, use a fallback color
                missing_vals <- unique_vals[is.na(color_vec)]
                if (length(missing_vals) > 0) {
                  cat(file=stderr(), paste0("  WARNING: Missing colors for: ", paste(missing_vals, collapse=", "), "\n"))
                  # Use grey for missing values
                  color_vec[is.na(color_vec)] <- "grey50"
                }
              } else {
                # Fallback: custom_colors is unnamed, use positional assignment (sorted order)
                cat(file=stderr(), paste0("  custom_colors is unnamed, using positional assignment\n"))
                sorted_vals <- sort(unique_vals)
                if (length(custom_colors) >= length(sorted_vals)) {
                  color_vec <- setNames(custom_colors[1:length(sorted_vals)], sorted_vals)
                } else {
                  color_vec <- setNames(rep(custom_colors, length.out = length(sorted_vals)), sorted_vals)
                }
              }
              cat(file=stderr(), paste0("  Color mapping: ", paste(names(color_vec), "=", color_vec, collapse=", "), "\n"))

              p_with_tiles <- p_with_tiles + scale_fill_manual(
                values = color_vec,
                name = heatmap_title,
                na.value = na_color
              )
            } else if (!is.null(custom_colors) && length(custom_colors) == 1 &&
                       custom_colors %in% rownames(RColorBrewer::brewer.pal.info)) {
              # v100: Use RColorBrewer palette
              cat(file=stderr(), paste0("  Using RColorBrewer palette: ", custom_colors, "\n"))
              p_with_tiles <- p_with_tiles + scale_fill_brewer(
                palette = custom_colors,
                name = heatmap_title,
                na.value = na_color
              )
            } else {
              # v100: Fallback to viridis
              cat(file=stderr(), paste0("  Using default viridis palette\n"))
              p_with_tiles <- p_with_tiles + scale_fill_viridis_d(
                name = heatmap_title,
                na.value = na_color
              )
            }
          } else {
            cat(file=stderr(), paste0("  Adding continuous color scale\n"))

            # v100: Get continuous scale colors from parameters
            low_color <- if (!is.null(heat_param['low']) && !is.na(heat_param['low'])) heat_param['low'] else "beige"
            mid_color <- if (!is.null(heat_param['mid']) && !is.na(heat_param['mid'])) heat_param['mid'] else "seashell2"
            high_color <- if (!is.null(heat_param['high']) && !is.na(heat_param['high'])) heat_param['high'] else "firebrick4"
            midpoint <- if (!is.null(heat_param['midpoint']) && !is.na(heat_param['midpoint'])) as.numeric(heat_param['midpoint']) else 0.02
            limits <- heat_param[['limits']]

            # v113: Debug output for continuous scale colors including NA color
            cat(file=stderr(), paste0("  Colors: low=", low_color, ", mid=", mid_color, ", high=", high_color, "\n"))
            cat(file=stderr(), paste0("  Midpoint: ", midpoint, "\n"))
            cat(file=stderr(), paste0("  v113: na_color for continuous: ", na_color, "\n"))
            cat(file=stderr(), paste0("  v113: heat_param na_color value: ", ifelse(is.null(heat_param[['na_color']]), "NULL", heat_param[['na_color']]), "\n"))

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

          cat(file=stderr(), paste0("  Color scale added successfully\n"))

          # v105: Add row labels if enabled
          show_row_labels <- if (!is.null(heat_param[['show_row_labels']])) heat_param[['show_row_labels']] else FALSE
          if (show_row_labels) {
            cat(file=stderr(), paste0("  Adding row labels...\n"))

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
              cat(file=stderr(), paste0("  Using label mapping for ", length(labels_to_use), " columns\n"))
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

            cat(file=stderr(), paste0("  Row labels: ", paste(labels_to_use, collapse=", "), "\n"))

            # v113: Improved row labels positioning
            # Labels appear below the heatmap (at lower y values than the tips)
            # With scale_y_reverse + coord_flip, this places them visually to the right

            # v113: Get row label offset and alignment from heat_param
            row_label_offset <- if (!is.null(heat_param[['row_label_offset']])) as.numeric(heat_param[['row_label_offset']]) else 1.0
            row_label_align <- if (!is.null(heat_param[['row_label_align']])) heat_param[['row_label_align']] else "left"

            cat(file=stderr(), paste0("  v113: row_label_offset from heat_param: ", row_label_offset, "\n"))
            cat(file=stderr(), paste0("  v113: row_label_align from heat_param: ", row_label_align, "\n"))

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
            cat(file=stderr(), paste0("  Row labels added successfully\n"))
            cat(file=stderr(), paste0("  Label position: label_y_pos=", label_y_pos, ", angle=", colnames_angle, ", hjust=", hjust_val, ", vjust=", vjust_val, "\n"))
          }

          # v116/v119: Add tip guide lines (vertical lines from tips through heatmap)
          show_guides <- if (!is.null(heat_param[['show_guides']])) heat_param[['show_guides']] else FALSE
          show_guides_bool <- isTRUE(show_guides) || identical(show_guides, TRUE) || identical(show_guides, "yes") || identical(show_guides, "TRUE")

          # v119: Debug output to trace guide line settings
          cat(file=stderr(), paste0("\n=== v119: TIP GUIDE LINES CHECK ===\n"))
          cat(file=stderr(), paste0("  show_guides raw value: ", show_guides, " (class: ", class(show_guides), ")\n"))
          cat(file=stderr(), paste0("  show_guides_bool: ", show_guides_bool, "\n"))

          if (show_guides_bool) {
            cat(file=stderr(), paste0("\n=== v116: ADDING TIP GUIDE LINES ===\n"))

            # Get guide line settings
            guide_color1 <- if (!is.null(heat_param[['guide_color1']])) heat_param[['guide_color1']] else "#CCCCCC"
            guide_color2 <- if (!is.null(heat_param[['guide_color2']])) heat_param[['guide_color2']] else "#EEEEEE"
            guide_alpha <- if (!is.null(heat_param[['guide_alpha']])) as.numeric(heat_param[['guide_alpha']]) else 0.3
            guide_width <- if (!is.null(heat_param[['guide_width']])) as.numeric(heat_param[['guide_width']]) else 0.5

            cat(file=stderr(), paste0("  guide_color1: ", guide_color1, "\n"))
            cat(file=stderr(), paste0("  guide_color2: ", guide_color2, "\n"))
            cat(file=stderr(), paste0("  guide_alpha: ", guide_alpha, "\n"))
            cat(file=stderr(), paste0("  guide_width: ", guide_width, "\n"))

            # v125: Get tip positions from tip_data to start guide lines at actual tip locations
            # Each tip may have a different x position based on branch lengths
            n_tips <- nrow(tip_data)
            cat(file=stderr(), paste0("  Number of tips: ", n_tips, "\n"))

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

            cat(file=stderr(), paste0("  v125: Guide lines from individual tip x positions to x_max=", x_max, "\n"))
            cat(file=stderr(), paste0("  v125: Tip x range: [", min(tip_data$x), ", ", max(tip_data$x), "]\n"))

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
                inherit.aes = FALSE
              ) +
              geom_segment(
                data = guide_lines_df[guide_lines_df$color_idx == 1, ],
                aes(x = x, xend = xend, y = y, yend = yend),
                color = guide_color2_alpha,
                linewidth = guide_width,
                inherit.aes = FALSE
              )

            cat(file=stderr(), paste0("  v116: Added ", n_tips, " tip guide lines\n"))
          }

          cat(file=stderr(), paste0("  Final layers: ", length(p_with_tiles$layers), "\n"))
          p_with_tiles

        }, error = function(e) {
          cat(file=stderr(), paste0("  ERROR adding heatmap: ", e$message, "\n"))
          cat(file=stderr(), paste0("  Returning tree without heatmap\n"))
          p
        })

        cat(file=stderr(), paste0("=== v122: HEATMAP ", heat_idx, " COMPLETE ===\n"))
        cat(file=stderr(), paste0("  p layers after heatmap: ", length(p$layers), "\n"))
        cat(file=stderr(), paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))

        # v125: Update current_heatmap_x_start for next heatmap
        # Calculate the rightmost x position of this heatmap + global gap
        this_heatmap_x_end <- max(tile_df$x) + tile_width / 2
        current_heatmap_x_start <- this_heatmap_x_end + heatmap_spacing  # v125: Add gap between heatmaps
        cat(file=stderr(), paste0("  v125: Updated current_heatmap_x_start to ", current_heatmap_x_start, " (added gap=", heatmap_spacing, ")\n"))

      } else {
        cat(file=stderr(), paste0("  WARNING: No tile data created - skipping heatmap ", heat_idx, "\n"))
      }
    } # End of v122 for loop for this heatmap
  }
  # END v122 MULTIPLE HEATMAPS

  # v94: Track p right after heatmap block
  cat(file=stderr(), paste0("\n=== v94: Immediately after simplified heatmap block ===\n"))
  cat(file=stderr(), paste0("  heat_flag: ", heat_flag, "\n"))
  cat(file=stderr(), paste0("  p layers: ", length(p$layers), "\n"))
  cat(file=stderr(), paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))

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
    cat(file=stderr(), paste0("\n=== v63: Pre-scaling x range: [", x_range_before[1], ", ", x_range_before[2], "] ===\n"))

    for (i in cc_totss) {
      par <- tt$data$parent[i]
      parX <- tt$data$x[par]
      tt$data[tt$data$node[i], "x"] <- tt$data[tt$data$node[i], "x"] * 15
    }

    # v63: DEBUG - show x range AFTER scaling
    x_range_after <- range(tt$data$x, na.rm = TRUE)
    cat(file=stderr(), paste0("=== v63: Post-scaling x range: [", x_range_after[1], ", ", x_range_after[2], "] ===\n"))
    
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
      if (is.na(flag_colnames[1])) {
        for (a in 1:length(heat_map_title_list)) {
          flag_colnames[a] <- FALSE
        }
      }
      
      if (flag_colnames[j1] == FALSE) {
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
        cat(file=stderr(), paste0("\n=== v89: Converting discrete heatmap to factors ===\n"))
        for (col_idx in 1:ncol(dxdf440_for_heat[[j1]])) {
          col_name <- colnames(dxdf440_for_heat[[j1]])[col_idx]
          col_vals <- dxdf440_for_heat[[j1]][, col_idx]
          unique_vals <- sort(unique(na.omit(col_vals)))
          cat(file=stderr(), paste0("  Column '", col_name, "': ", length(unique_vals), " unique values\n"))
          cat(file=stderr(), paste0("  Unique values: ", paste(head(unique_vals, 10), collapse=", "),
                                    if(length(unique_vals) > 10) "..." else "", "\n"))
          # v89: RESTORED - Convert to factor with sorted levels (REQUIRED for discrete heatmaps)
          dxdf440_for_heat[[j1]][, col_idx] <- factor(col_vals, levels = unique_vals)
          cat(file=stderr(), paste0("  Converted to factor with ", length(unique_vals), " levels\n"))
        }
        cat(file=stderr(), paste0("================================\n"))
      }

      # Create the heatmap
      # v61: DEBUG - show data structure before gheatmap call
      cat(file=stderr(), paste0("\n=== v61: GHEATMAP DATA DEBUG ===\n"))
      cat(file=stderr(), paste0("  Heatmap index: j1=", j1, ", j=", j, "\n"))
      heat_data <- dxdf440_for_heat[[j1]]
      cat(file=stderr(), paste0("  heat_data dimensions: ", nrow(heat_data), " x ", ncol(heat_data), "\n"))
      cat(file=stderr(), paste0("  heat_data rownames sample: ", paste(head(rownames(heat_data), 5), collapse=", "), "\n"))
      cat(file=stderr(), paste0("  heat_data columns: ", paste(colnames(heat_data), collapse=", "), "\n"))
      tt_tips <- subset(tt$data, isTip == TRUE)
      cat(file=stderr(), paste0("  Tree tip labels sample: ", paste(head(tt_tips$label, 5), collapse=", "), "\n"))
      # Check if rownames match tree tip labels
      matching_tips <- sum(rownames(heat_data) %in% tt_tips$label)
      cat(file=stderr(), paste0("  Rownames matching tree tips: ", matching_tips, " / ", nrow(heat_data), "\n"))
      cat(file=stderr(), paste0("  offset (new_heat_x): ", new_heat_x, "\n"))
      cat(file=stderr(), paste0("  width (wi): ", wi, "\n"))
      cat(file=stderr(), paste0("================================\n"))

      # v83: RESTORED duplicate gheatmap pattern from v61 - THIS IS INTENTIONAL
      # The user confirmed this pattern was in the original lineage plotter code and is required.
      # DO NOT REMOVE THIS DUPLICATE CALL - it is necessary for gheatmap to work correctly.
      # First gheatmap call creates the initial structure
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
            cat(file=stderr(), paste0("  v85: Found tile layer ", layer_idx, " with values: ",
                                      paste(head(layer_vals, 5), collapse=", "), "\n"))
          }
        }
      }
      tile_values <- unique(tile_values)
      cat(file=stderr(), paste0("  v85: All tile values: ", paste(tile_values, collapse=", "), "\n"))

      # v82: DEBUG - verify gheatmap result and tile layer data
      cat(file=stderr(), paste0("\n=== v71: POST-GHEATMAP DEBUG ===\n"))
      cat(file=stderr(), paste0("  Number of layers in plot: ", length(pr440_short_tips_TRY_heat$layers), "\n"))
      gheatmap_xrange <- range(pr440_short_tips_TRY_heat$data$x, na.rm = TRUE)
      cat(file=stderr(), paste0("  Plot data x range: [", gheatmap_xrange[1], ", ", gheatmap_xrange[2], "]\n"))

      # Check if there's rect/tile data (heatmap)
      layer_types <- sapply(pr440_short_tips_TRY_heat$layers, function(l) class(l$geom)[1])
      cat(file=stderr(), paste0("  Layer geom types: ", paste(layer_types, collapse=", "), "\n"))

      # v72: Enhanced debugging to find and inspect the GeomTile layer
      for (layer_idx in seq_along(pr440_short_tips_TRY_heat$layers)) {
        layer <- pr440_short_tips_TRY_heat$layers[[layer_idx]]
        if (inherits(layer$geom, "GeomTile")) {
          cat(file=stderr(), paste0("  GeomTile found at layer ", layer_idx, "\n"))
          tryCatch({
            layer_data <- layer$data
            if (is.function(layer_data)) {
              cat(file=stderr(), paste0("    Layer data is a function, evaluating...\n"))
              layer_data <- layer_data(pr440_short_tips_TRY_heat$data)
            }
            if (!is.null(layer_data) && is.data.frame(layer_data)) {
              cat(file=stderr(), paste0("    Layer data rows: ", nrow(layer_data), "\n"))
              cat(file=stderr(), paste0("    Layer data columns: ", paste(names(layer_data), collapse=", "), "\n"))
              if ("x" %in% names(layer_data)) {
                x_vals <- layer_data$x
                cat(file=stderr(), paste0("    Tile x range: [", min(x_vals, na.rm=TRUE), ", ", max(x_vals, na.rm=TRUE), "]\n"))
              }
              if ("y" %in% names(layer_data)) {
                y_vals <- layer_data$y
                cat(file=stderr(), paste0("    Tile y range: [", min(y_vals, na.rm=TRUE), ", ", max(y_vals, na.rm=TRUE), "]\n"))
              }
              # v72: Check value column (this is what fill maps to in gheatmap)
              if ("value" %in% names(layer_data)) {
                cat(file=stderr(), paste0("    Value column (unique): ", paste(unique(layer_data$value), collapse=", "), "\n"))
                cat(file=stderr(), paste0("    Value column class: ", class(layer_data$value)[1], "\n"))
              }
              # v72: Check width column (tile width)
              if ("width" %in% names(layer_data)) {
                width_vals <- layer_data$width
                cat(file=stderr(), paste0("    Width values: [", min(width_vals, na.rm=TRUE), ", ", max(width_vals, na.rm=TRUE), "]\n"))
              }
              if ("fill" %in% names(layer_data)) {
                cat(file=stderr(), paste0("    Fill values (first 5): ", paste(head(layer_data$fill, 5), collapse=", "), "\n"))
              }
              # v72: Check the layer's aesthetic mapping
              if (!is.null(layer$mapping)) {
                cat(file=stderr(), paste0("    Layer mapping: ", paste(names(layer$mapping), collapse=", "), "\n"))
              }
            } else {
              cat(file=stderr(), paste0("    WARNING: Layer data is not a data.frame: ", class(layer_data)[1], "\n"))
            }
          }, error = function(e) {
            cat(file=stderr(), paste0("    ERROR accessing layer data: ", e$message, "\n"))
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
        cat(file=stderr(), paste0("\n=== v65: DISCRETE HEATMAP COLOR DEBUG ===\n"))
        cat(file=stderr(), paste0("  heat_param['is_discrete']: ", heat_param['is_discrete'], "\n"))
        cat(file=stderr(), paste0("  heat_param['man']: ", heat_param['man'], "\n"))
        cat(file=stderr(), paste0("  heat_param['man_define_colors']: ", heat_param['man_define_colors'], "\n"))
        cat(file=stderr(), paste0("  heat_param[['color_scale_option']]: ",
                                  if(is.null(heat_param[['color_scale_option']])) "NULL"
                                  else paste(class(heat_param[['color_scale_option']]), collapse=", "), "\n"))
        if (!is.null(heat_param[['color_scale_option']])) {
          cat(file=stderr(), paste0("  color_scale_option value: ",
                                    paste(heat_param[['color_scale_option']], collapse=", "), "\n"))
          if (is.list(heat_param[['color_scale_option']])) {
            cat(file=stderr(), paste0("  color_scale_option$color_scale_option: ",
                                      heat_param[['color_scale_option']]$color_scale_option, "\n"))
          }
        }
        cat(file=stderr(), paste0("========================================\n"))

        if (heat_param['man'] == FALSE) {
          if (heat_param['man_define_colors'] == FALSE) {
            # v67: SIMPLIFIED - color_scale_option is now directly a palette name string (e.g., "Set1")
            # No more nested list handling needed after the fix at line 2616
            palette_name <- heat_param[['color_scale_option']]

            cat(file=stderr(), paste0("\n=== v67: Discrete palette check ===\n"))
            cat(file=stderr(), paste0("  palette_name: ", if(is.null(palette_name)) "NULL" else palette_name, "\n"))
            if (!is.null(palette_name)) {
              cat(file=stderr(), paste0("  Is valid RColorBrewer palette: ",
                                        palette_name %in% rownames(RColorBrewer::brewer.pal.info), "\n"))
            }
            cat(file=stderr(), paste0("================================\n"))

            # v67: Use RColorBrewer palette if available, otherwise use default hue
            if (!is.null(palette_name) && is.character(palette_name) &&
                palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
              # Get unique values to determine number of colors needed
              heat_data_vals <- unique(na.omit(dxdf440_for_heat[[j1]][,1]))
              n_vals <- length(heat_data_vals)
              max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
              n_colors <- min(n_vals, max_colors)

              cat(file=stderr(), paste0("\n=== v67: Applying discrete palette ===\n"))
              cat(file=stderr(), paste0("  Palette: ", palette_name, "\n"))
              cat(file=stderr(), paste0("  Number of unique values: ", n_vals, "\n"))
              cat(file=stderr(), paste0("  Colors to use: ", n_colors, "\n"))
              cat(file=stderr(), paste0("================================\n"))

              # v70: Get NA color from heat_param (default white)
              na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "white"
              cat(file=stderr(), paste0("  NA color: ", na_color, "\n"))

              pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
                scale_fill_brewer(palette = palette_name, name = heat_map_title_list[[j1]], na.value = na_color)
              # v82: Removed aggressive mapping repair - will do single repair at end of heatmap loop
            } else {
              cat(file=stderr(), paste0("  v67: Using default hue scale (no valid palette specified)\n"))
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
            cat(file=stderr(), paste0("\n=== v70: Applying custom discrete colors ===\n"))
            cat(file=stderr(), paste0("  custom_colors class: ", paste(class(custom_colors), collapse=", "), "\n"))
            cat(file=stderr(), paste0("  custom_colors length: ", length(custom_colors), "\n"))
            cat(file=stderr(), paste0("  custom_colors names: ", paste(head(names(custom_colors), 10), collapse=", "), "\n"))
            cat(file=stderr(), paste0("  custom_colors values: ", paste(head(custom_colors, 10), collapse=", "), "\n"))

            # v85: Use tile_values collected from actual tile layer data (not dxdf440_for_heat)
            # This ensures we match exactly what gheatmap created
            heat_data_vals <- tile_values
            cat(file=stderr(), paste0("  v85: Using tile_values from layer: ", paste(head(heat_data_vals, 10), collapse=", "), "\n"))

            # v73: Fix - properly subset custom_colors to match tile values
            # This prevents NA names which cause scale_fill_manual to fail
            n_levels <- length(heat_data_vals)

            if (n_levels == 0) {
              # v85: If no tile values found, fall back to using dxdf440_for_heat
              cat(file=stderr(), paste0("  v85: WARNING - no tile values found, using source data\n"))
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
              cat(file=stderr(), paste0("  v73: WARNING - custom_colors has no names, subsetting and assigning by position\n"))
              cat(file=stderr(), paste0("  v73: Number of tile values: ", n_levels, "\n"))
              cat(file=stderr(), paste0("  v73: Number of custom colors: ", length(custom_colors), "\n"))

              # CRITICAL: Only use as many colors as there are values
              if (n_levels > 0 && length(custom_colors) >= n_levels) {
                # Subset colors to match values exactly
                colors_to_use <- custom_colors[1:n_levels]
                names(colors_to_use) <- as.character(heat_data_vals)
                custom_colors <- colors_to_use
                cat(file=stderr(), paste0("  v73: Subsetted to ", n_levels, " colors\n"))
                cat(file=stderr(), paste0("  v73: Final names: ", paste(names(custom_colors), collapse=", "), "\n"))
                cat(file=stderr(), paste0("  v73: Final values: ", paste(custom_colors, collapse=", "), "\n"))
              } else if (n_levels > 0) {
                # Not enough colors, recycle
                cat(file=stderr(), paste0("  v73: WARNING - not enough colors, recycling\n"))
                colors_to_use <- rep(custom_colors, length.out = n_levels)
                names(colors_to_use) <- as.character(heat_data_vals)
                custom_colors <- colors_to_use
              }
            } else {
              # custom_colors already has names - ensure they match values
              cat(file=stderr(), paste0("  v73: custom_colors already named: ", paste(names(custom_colors), collapse=", "), "\n"))
              # Only keep colors whose names are in heat_data_vals
              valid_names <- names(custom_colors) %in% as.character(heat_data_vals)
              if (any(valid_names)) {
                custom_colors <- custom_colors[valid_names]
                cat(file=stderr(), paste0("  v73: After filtering: ", paste(names(custom_colors), collapse=", "), "\n"))
              }
            }
            # v70: Get NA color from heat_param (default white)
            na_color <- if (!is.null(heat_param[['na_color']])) heat_param[['na_color']] else "white"
            cat(file=stderr(), paste0("  NA color: ", na_color, "\n"))
            cat(file=stderr(), paste0("================================\n"))

            # v88: CRITICAL FIX - Test ggplot_build BEFORE applying scale to diagnose root cause
            cat(file=stderr(), paste0("  v88: Testing ggplot_build BEFORE scale application...\n"))
            pre_scale_ok <- tryCatch({
              test_build_pre <- ggplot2::ggplot_build(pr440_short_tips_TRY_heat)
              cat(file=stderr(), paste0("  v88: Pre-scale ggplot_build: SUCCESS\n"))
              TRUE
            }, error = function(e) {
              cat(file=stderr(), paste0("  v88: Pre-scale ggplot_build: FAILED - ", e$message, "\n"))
              cat(file=stderr(), paste0("  v88: ERROR IS IN GHEATMAP OUTPUT, NOT SCALE\n"))
              # Get more info about the error
              cat(file=stderr(), paste0("  v88: Checking individual layers...\n"))
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
                    cat(file=stderr(), paste0("    Layer ", li, " (", class(layer$geom)[1], "): mapping=", paste(mapping_names, collapse=","), "\n"))
                  } else {
                    cat(file=stderr(), paste0("    Layer ", li, " (", class(layer$geom)[1], "): no mapping\n"))
                  }
                  TRUE
                }, error = function(e2) {
                  cat(file=stderr(), paste0("    Layer ", li, " (", class(pr440_short_tips_TRY_heat$layers[[li]]$geom)[1], "): FAILED - ", e2$message, "\n"))
                  FALSE
                })
              }
              FALSE
            })

            # v88: Apply scale regardless of pre-test result (the scale itself isn't the problem)
            cat(file=stderr(), paste0("  v88: Applying scale_fill_manual...\n"))
            tryCatch({
              pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat +
                scale_fill_manual(
                  values = custom_colors,
                  name = heat_map_title_list[[j1]],
                  na.value = na_color
                )
              cat(file=stderr(), paste0("  v88: scale_fill_manual applied\n"))

              # v88: Test after scale
              post_scale_ok <- tryCatch({
                test_build_post <- ggplot2::ggplot_build(pr440_short_tips_TRY_heat)
                cat(file=stderr(), paste0("  v88: Post-scale ggplot_build: SUCCESS\n"))
                TRUE
              }, error = function(e) {
                cat(file=stderr(), paste0("  v88: Post-scale ggplot_build: FAILED - ", e$message, "\n"))
                FALSE
              })

            }, error = function(e) {
              cat(file=stderr(), paste0("  v88: scale_fill_manual failed: ", e$message, "\n"))
              cat(file=stderr(), paste0("  v88: Trying scale_fill_discrete as fallback\n"))
              tryCatch({
                pr440_short_tips_TRY_heat <<- pr440_short_tips_TRY_heat +
                  scale_fill_discrete(name = heat_map_title_list[[j1]], na.value = na_color)
              }, error = function(e2) {
                cat(file=stderr(), paste0("  v88: scale_fill_discrete also failed: ", e2$message, "\n"))
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
      cat(file=stderr(), paste0("\n=== v68: After scale application ===\n"))
      cat(file=stderr(), paste0("  j=", j, ", j1=", j1, "\n"))
      cat(file=stderr(), paste0("================================\n"))

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
      cat(file=stderr(), paste0("\n=== v68: Before func.calc.min_col_x ===\n"))

      # Store position for next heatmap (wrapped in tryCatch to prevent errors from breaking the loop)
      min_col_x_of_frame_of_prev_heat <- tryCatch({
        func.calc.min_col_x_of_frame_of_prev_heat(pr440_short_tips_TRY_heat, colnames_sub_df_heat_j)
      }, error = function(e) {
        cat(file=stderr(), paste0("  v68: WARNING - func.calc.min_col_x_of_frame_of_prev_heat failed: ", e$message, "\n"))
        0  # Return default value
      })

      cat(file=stderr(), paste0("  v68: min_col_x_of_frame_of_prev_heat = ", min_col_x_of_frame_of_prev_heat, "\n"))
      cat(file=stderr(), paste0("================================\n"))
    }

    # v68: DEBUG - after for loop
    cat(file=stderr(), paste0("\n=== v68: After heatmap for loop ===\n"))

    # Update plot with heatmap
    p <- pr440_short_tips_TRY_heat

    # v71: Repair mapping before final state check
    p <- func.repair.ggtree.mapping(p, verbose = TRUE)

    # v71: Diagnose any layer issues
    problematic_layers <- func.diagnose.layer.issues(p, verbose = TRUE)

    # v71: Final heatmap state with layer analysis
    cat(file=stderr(), paste0("\n=== v71: FINAL HEATMAP STATE ===\n"))
    final_xrange <- range(p$data$x, na.rm = TRUE)
    cat(file=stderr(), paste0("  Tree data x range: [", final_xrange[1], ", ", final_xrange[2], "]\n"))
    cat(file=stderr(), paste0("  Number of layers: ", length(p$layers), "\n"))

    # v71: Check for heatmap layer (GeomTile) and get its x coordinates
    layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
    cat(file=stderr(), paste0("  Layer types: ", paste(layer_types, collapse=", "), "\n"))

    # v71: Find the GeomTile layer (heatmap) and get its data directly
    heatmap_xmax <- NULL
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      if (inherits(layer$geom, "GeomTile")) {
        cat(file=stderr(), paste0("  Found GeomTile at layer ", i, "\n"))
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
              cat(file=stderr(), paste0("    GeomTile x range: [", min(tile_x, na.rm = TRUE), ", ", tile_xmax, "]\n"))
              heatmap_xmax <- tile_xmax
            }
          } else {
            cat(file=stderr(), paste0("    GeomTile layer data does not have x column\n"))
          }
        }, error = function(e) {
          cat(file=stderr(), paste0("    Could not access GeomTile data: ", e$message, "\n"))
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
              cat(file=stderr(), paste0("    Built layer ", i, " x range: [",
                                        min(layer_x, na.rm = TRUE), ", ", layer_xmax, "]\n"))
              if (is.null(heatmap_xmax) || layer_xmax > heatmap_xmax) {
                heatmap_xmax <- layer_xmax
              }
            }
          }
        }
      }, error = function(e) {
        cat(file=stderr(), paste0("  v71: ggplot_build failed (will use fallback): ", e$message, "\n"))
      })
    }
    cat(file=stderr(), paste0("================================\n"))

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

    cat(file=stderr(), paste0("\n=== v71: EXPANDING X-AXIS FOR HEATMAP ===\n"))
    cat(file=stderr(), paste0("  Tree x range: [", final_xrange[1], ", ", final_xrange[2], "]\n"))
    cat(file=stderr(), paste0("  Tree width: ", tree_width, "\n"))
    cat(file=stderr(), paste0("  Heatmap offset (new_heat_x): ", new_heat_x, "\n"))
    cat(file=stderr(), paste0("  Heatmap width (wi): ", wi, "\n"))
    cat(file=stderr(), paste0("  Calculated x max: ", calculated_xmax, "\n"))
    cat(file=stderr(), paste0("  Proportional x max: ", proportional_xmax, "\n"))
    cat(file=stderr(), paste0("  Fixed margin x max: ", fixed_margin_xmax, "\n"))
    cat(file=stderr(), paste0("  Detected from build: ", ifelse(is.null(heatmap_xmax), "NULL", heatmap_xmax), "\n"))
    cat(file=stderr(), paste0("  Final expected max x: ", expected_xmax, "\n"))
    cat(file=stderr(), paste0("  Setting coord_flip xlim to: [", final_xrange[1], ", ", expected_xmax, "]\n"))
    cat(file=stderr(), paste0("========================================\n"))

    # v88: SIMPLIFIED EXPANSION - Skip expansion functions if plot can't be built
    # Test if plot is buildable before trying expansion
    plot_buildable <- tryCatch({
      ggplot2::ggplot_build(p)
      cat(file=stderr(), paste0("  v88: Plot is buildable, proceeding with expansion\n"))
      TRUE
    }, error = function(e) {
      cat(file=stderr(), paste0("  v88: Plot is NOT buildable: ", e$message, "\n"))
      cat(file=stderr(), paste0("  v88: Skipping expansion - will render plot as-is\n"))
      FALSE
    })

    if (plot_buildable) {
      # Calculate how much expansion is needed on the right side (positive x direction)
      expansion_ratio <- (expected_xmax - final_xrange[2]) / abs(final_xrange[1] - final_xrange[2])
      expansion_ratio <- max(0.3, expansion_ratio)

      cat(file=stderr(), paste0("  v88: Using hexpand() with ratio: ", expansion_ratio, "\n"))

      # v88: Try hexpand first
      expansion_success <- FALSE
      tryCatch({
        p <- p + ggtree::hexpand(ratio = expansion_ratio, direction = 1)
        expansion_success <- TRUE
        cat(file=stderr(), paste0("  v88: hexpand applied successfully\n"))
      }, error = function(e) {
        cat(file=stderr(), paste0("  v88: hexpand failed: ", e$message, "\n"))
      })

      # v88: Fallback to xlim_expand
      if (!expansion_success) {
        cat(file=stderr(), paste0("  v88: Trying xlim_expand fallback\n"))
        tryCatch({
          p <- p + ggtree::xlim_expand(c(0, expected_xmax), "right")
          expansion_success <- TRUE
          cat(file=stderr(), paste0("  v88: xlim_expand applied\n"))
        }, error = function(e2) {
          cat(file=stderr(), paste0("  v88: xlim_expand also failed: ", e2$message, "\n"))
        })
      }

      # v88: Final fallback - just add geom_blank with expanded limits
      if (!expansion_success) {
        tryCatch({
          # Use geom_blank to expand plot limits without modifying coordinate system
          p <- p + geom_blank(data = data.frame(x = c(final_xrange[1], expected_xmax), y = c(1, 1)),
                              aes(x = x, y = y))
          cat(file=stderr(), paste0("  v88: geom_blank expansion applied\n"))
        }, error = function(e3) {
          cat(file=stderr(), paste0("  v88: All expansion methods failed\n"))
        })
      }

      # Repair mapping after changes
      p <- func.repair.ggtree.mapping(p)
    } else {
      # v88: Plot is not buildable - don't add any expansion, just try to render
      cat(file=stderr(), paste0("  v88: WARNING - Plot cannot be built, skipping all expansion\n"))
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
      cat(file=stderr(), paste0("  v132: Error computing boudariestt (2nd block): ", e$message, "\n"))
      list(xmin = min(p$data$x, na.rm = TRUE), xmax = max(p$data$x, na.rm = TRUE))
    })
    cat(file=stderr(), paste0("  v132: boudariestt initialized (2nd block): xmin=", boudariestt$xmin, ", xmax=", boudariestt$xmax, "\n"))
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

    cat(file=stderr(), paste0("\n=== v146: ELLIPSE SCALING FOR HEATMAP ===\n"))
    cat(file=stderr(), paste0("  Tree width: ", round(tree_width, 2), "\n"))
    cat(file=stderr(), paste0("  Total plot width: ", round(total_width, 2), "\n"))
    cat(file=stderr(), paste0("  Tree ratio (tree/total): ", round(tree_ratio, 3), "\n"))

    # v146: Scale both height (a) and width (b) based on tree ratio
    # Use the same base calculation as non-heatmap case, then scale by tree_ratio
    base_a <- (tip_length * size_tip_text / 800 + man_adjust_elipse_a) * tree_ratio
    base_b <- (0.12 + man_adjust_elipse_b) * tree_ratio

    a <- base_a * adjust_height_ecliplse
    b <- base_b * adjust_width_eclipse

    cat(file=stderr(), paste0("  Base a (height): ", round(base_a, 4), "\n"))
    cat(file=stderr(), paste0("  Base b (width): ", round(base_b, 4), "\n"))
    cat(file=stderr(), paste0("  Final a (after user adjust): ", round(a, 4), "\n"))
    cat(file=stderr(), paste0("  Final b (after user adjust): ", round(b, 4), "\n"))
    cat(file=stderr(), paste0("=============================================\n"))

    # v139: Pass high_alpha_list for transparency
    p <- func_highlight(
      p, how_many_hi, heat_flag, high_color_list, a, b, man_adjust_elipse,
      pr440_short_tips_TRY, boudariestt, debug_mode, high_offset, high_vertical_offset,
      high_alpha_list
    )
  }
  
  if (length(b) == 0) {
    b <- 0.2
  }
  
  # v94: DEBUG - track layers after if(FALSE) block
  cat(file=stderr(), paste0("\n=== v94: AFTER if(FALSE) block ===\n"))
  cat(file=stderr(), paste0("  p layers: ", length(p$layers), "\n"))
  cat(file=stderr(), paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))
  cat(file=stderr(), paste0("================================\n"))

  # Add second legend if heatmap exists
  if (length(heat_map_title_list) > 0) {
    cat(file=stderr(), paste0("\n=== v95: Before func.make.second.legend ===\n"))
    cat(file=stderr(), paste0("  p layers: ", length(p$layers), "\n"))
    cat(file=stderr(), paste0("  heat_flag: ", heat_flag, "\n"))
    cat(file=stderr(), paste0("  FLAG_BULK_DISPLAY: ", FLAG_BULK_DISPLAY, "\n"))

    # v95: Debug boudariestt
    if (exists("boudariestt") && !is.null(boudariestt)) {
      cat(file=stderr(), paste0("  boudariestt$xmax: ", boudariestt$xmax, "\n"))
      cat(file=stderr(), paste0("  boudariestt$xmin: ", boudariestt$xmin, "\n"))
    } else {
      cat(file=stderr(), paste0("  WARNING: boudariestt is NULL or doesn't exist!\n"))
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

    cat(file=stderr(), paste0("  v133: highlight offsets - x:", highlight_x_off, ", y:", highlight_y_off, "\n"))
    cat(file=stderr(), paste0("  v133: bootstrap offsets - x:", bootstrap_x_off, ", y:", bootstrap_y_off, "\n"))
    cat(file=stderr(), paste0("  v138: show_highlight_legend:", show_highlight_leg, ", show_bootstrap_legend:", show_bootstrap_leg, "\n"))

    # v151: Calculate y_off_base to position legends on the RIGHT side
    # With coord_flip + scale_y_reverse: y values become horizontal positions
    # y = 1 is at visual LEFT, y = max_tips is at visual RIGHT
    # CRITICAL: Using y > max_tips EXPANDS the plot range and shrinks the tree!
    # Solution: Use y = max_tips (not max_tips + 5) to minimize expansion
    n_tips <- sum(p$data$isTip == TRUE, na.rm = TRUE)
    max_y_in_data <- max(p$data$y, na.rm = TRUE)
    y_off_base <- max(n_tips, max_y_in_data)  # Stay within existing data range
    cat(file=stderr(), paste0("  v151: y_off_base=", y_off_base, " (n_tips=", n_tips, ", max_y=", round(max_y_in_data, 2), ")\n"))

    # v95: Wrap in tryCatch to catch any errors
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
      cat(file=stderr(), paste0("  func.make.second.legend: SUCCESS\n"))
      result
    }, error = function(e) {
      cat(file=stderr(), paste0("  func.make.second.legend ERROR: ", e$message, "\n"))
      cat(file=stderr(), paste0("  Returning plot without second legend modifications\n"))
      p  # Return original plot on error
    })

    cat(file=stderr(), paste0("=== v95: After func.make.second.legend ===\n"))
    cat(file=stderr(), paste0("  p layers: ", length(p$layers), "\n"))
    cat(file=stderr(), paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))
  }

  # v94: DEBUG - before bootstrap triangles
  cat(file=stderr(), paste0("\n=== v94: Before bootstrap triangles ===\n"))
  cat(file=stderr(), paste0("  p layers: ", length(p$layers), "\n"))
  cat(file=stderr(), paste0("  Layer types: ", paste(sapply(p$layers, function(l) class(l$geom)[1]), collapse=", "), "\n"))

  # v96: FIX - Repair corrupted mapping BEFORE adding bootstrap triangles
  # gheatmap() corrupts the plot's @mapping attribute (changes it from ggplot2::mapping to data.frame)
  # This causes geom_nodepoint() calls below to crash silently
  # We must repair the mapping BEFORE adding any new layers
  cat(file=stderr(), paste0("\n=== v96: Repairing mapping before bootstrap triangles ===\n"))
  cat(file=stderr(), paste0("  Mapping class before repair: ", paste(class(p$mapping), collapse=", "), "\n"))
  p <- func.repair.ggtree.mapping(p, verbose = TRUE)
  cat(file=stderr(), paste0("  Mapping class after repair: ", paste(class(p$mapping), collapse=", "), "\n"))

  # v90: Add bootstrap triangles AFTER heatmap processing to avoid "missing x and y" error
  # When gheatmap transforms the plot data, any geom_nodepoint layers added beforehand
  # lose their x and y aesthetic mappings. By adding them here (after gheatmap),
  # the layers are created with the correct data structure.
  if (bootstrap_triangles_enabled && !is.null(bootstrap_triangles_params)) {
    cat(file=stderr(), paste0("\n=== v90: Adding bootstrap triangles after heatmap ===\n"))
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
      cat(file=stderr(), paste0("  Bootstrap triangles added successfully\n"))
    }, error = function(e) {
      cat(file=stderr(), paste0("  v96 ERROR adding bootstrap triangles: ", e$message, "\n"))
      cat(file=stderr(), paste0("  Continuing without bootstrap triangles\n"))
    })
    cat(file=stderr(), paste0("================================\n"))
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

  # v72: FINAL DEBUG - verify heatmap layer exists and inspect plot state
  cat(file=stderr(), paste0("\n=== v72: FINAL PLOT STATE BEFORE GGSAVE ===\n"))
  cat(file=stderr(), paste0("  Number of layers: ", length(p$layers), "\n"))
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  cat(file=stderr(), paste0("  Layer types: ", paste(layer_types, collapse=", "), "\n"))

  # Check for GeomTile (heatmap)
  geomtile_idx <- which(layer_types == "GeomTile")
  if (length(geomtile_idx) > 0) {
    cat(file=stderr(), paste0("  GeomTile found at layers: ", paste(geomtile_idx, collapse=", "), "\n"))
    for (idx in geomtile_idx) {
      tryCatch({
        tile_layer <- p$layers[[idx]]
        tile_data <- tile_layer$data
        if (is.function(tile_data)) {
          tile_data <- tile_data(p$data)
        }
        if (!is.null(tile_data) && is.data.frame(tile_data)) {
          cat(file=stderr(), paste0("    Layer ", idx, " data rows: ", nrow(tile_data), "\n"))
          if ("x" %in% names(tile_data)) {
            cat(file=stderr(), paste0("    Layer ", idx, " x range: [", min(tile_data$x, na.rm=TRUE), ", ", max(tile_data$x, na.rm=TRUE), "]\n"))
          }
          if ("value" %in% names(tile_data)) {
            cat(file=stderr(), paste0("    Layer ", idx, " values: ", paste(unique(tile_data$value), collapse=", "), "\n"))
          }
        }
      }, error = function(e) {
        cat(file=stderr(), paste0("    ERROR: ", e$message, "\n"))
      })
    }
  } else {
    cat(file=stderr(), paste0("  WARNING: No GeomTile layers found! Heatmap may not be displayed.\n"))
  }

  # Check coordinate system
  if (!is.null(p$coordinates)) {
    cat(file=stderr(), paste0("  Coordinate system: ", class(p$coordinates)[1], "\n"))
    if (inherits(p$coordinates, "CoordCartesian")) {
      if (!is.null(p$coordinates$limits$x)) {
        cat(file=stderr(), paste0("  X limits: [", p$coordinates$limits$x[1], ", ", p$coordinates$limits$x[2], "]\n"))
      }
    }
  }

  # Try a final ggplot_build to get computed values
  tryCatch({
    final_built <- ggplot2::ggplot_build(p)
    cat(file=stderr(), paste0("  ggplot_build successful\n"))

    # Find tile layer in built data
    for (i in seq_along(final_built$data)) {
      if ("fill" %in% names(final_built$data[[i]]) && "width" %in% names(final_built$data[[i]])) {
        bd <- final_built$data[[i]]
        cat(file=stderr(), paste0("  Built layer ", i, ": ", nrow(bd), " rows, x=[",
                                  min(bd$x, na.rm=TRUE), ", ", max(bd$x, na.rm=TRUE),
                                  "], fill=", paste(unique(bd$fill), collapse=","), "\n"))
      }
    }

    # Check panel ranges
    if (!is.null(final_built$layout$panel_params)) {
      for (panel_idx in seq_along(final_built$layout$panel_params)) {
        pp <- final_built$layout$panel_params[[panel_idx]]
        if (!is.null(pp$x.range)) {
          cat(file=stderr(), paste0("  Panel ", panel_idx, " x.range: [", pp$x.range[1], ", ", pp$x.range[2], "]\n"))
        }
      }
    }
  }, error = function(e) {
    cat(file=stderr(), paste0("  ggplot_build error: ", e$message, "\n"))
  })
  cat(file=stderr(), paste0("============================================\n"))

  # v88: Save the plot with comprehensive error handling and multiple fallbacks
  save_success <- FALSE

  # v158: Check if plot has custom legend info for gtable insertion
  custom_legend_info <- attr(p, "custom_legend_info")

  # v158: Primary attempt - with gtable legend insertion if needed
  tryCatch({
    if (!is.null(custom_legend_info) && (length(custom_legend_info$highlight) > 0 || length(custom_legend_info$bootstrap) > 0)) {
      # Build gtable from plot
      cat(file=stderr(), paste0("\n=== v158: Building gtable for custom legend insertion ===\n"))
      gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))

      # Add custom legends to gtable
      gt <- func.add.custom.legends.to.gtable(gt, custom_legend_info)

      # Save the modified gtable
      ggsave(out_file_path, plot = gt, width = width, height = height, units = units_out, limitsize = FALSE)
      save_success <- TRUE
      cat(file=stderr(), paste0("\n=== v158: Plot with custom legends saved successfully ===\n"))
    } else {
      # No custom legends, use standard ggsave
      ggsave(out_file_path, plot = p, width = width, height = height, units = units_out, limitsize = FALSE)
      save_success <- TRUE
      cat(file=stderr(), paste0("\n=== v88: Plot saved successfully ===\n"))
    }
  }, error = function(e) {
    cat(file=stderr(), paste0("\n=== v88: GGSAVE ERROR ===\n"))
    cat(file=stderr(), paste0("  Primary error: ", e$message, "\n"))
  })

  # v88: Fallback 1 - repair mapping and try again
  if (!save_success) {
    cat(file=stderr(), paste0("  v88: Trying fallback 1 - repair mapping\n"))
    tryCatch({
      p_repaired <- func.repair.ggtree.mapping(p, verbose = TRUE)
      ggsave(out_file_path, plot = p_repaired, width = width, height = height,
             units = units_out, limitsize = FALSE)
      save_success <- TRUE
      cat(file=stderr(), paste0("  v88: Fallback 1 succeeded\n"))
    }, error = function(e2) {
      cat(file=stderr(), paste0("  v88: Fallback 1 failed: ", e2$message, "\n"))
    })
  }

  # v88: Fallback 2 - remove heatmap scale and try with defaults
  if (!save_success && heat_flag) {
    cat(file=stderr(), paste0("  v88: Trying fallback 2 - reset fill scale to defaults\n"))
    tryCatch({
      # Create fresh plot with default scale
      p_default <- p + scale_fill_discrete(na.value = "white")
      p_default <- func.repair.ggtree.mapping(p_default)
      ggsave(out_file_path, plot = p_default, width = width, height = height,
             units = units_out, limitsize = FALSE)
      save_success <- TRUE
      cat(file=stderr(), paste0("  v88: Fallback 2 succeeded (using default colors)\n"))
    }, error = function(e3) {
      cat(file=stderr(), paste0("  v88: Fallback 2 failed: ", e3$message, "\n"))
    })
  }

  # v88: Fallback 3 - use print() to render instead of ggsave
  if (!save_success) {
    cat(file=stderr(), paste0("  v88: Trying fallback 3 - direct PNG rendering\n"))
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
      cat(file=stderr(), paste0("  v88: Fallback 3 succeeded\n"))
    }, error = function(e4) {
      cat(file=stderr(), paste0("  v88: Fallback 3 failed: ", e4$message, "\n"))
      tryCatch(dev.off(), error = function(x) {})  # Clean up device
    })
  }

  # v88: Final fallback - save tree without heatmap
  if (!save_success && heat_flag) {
    cat(file=stderr(), paste0("  v88: Trying fallback 4 - save tree without heatmap\n"))
    tryCatch({
      # Use the original tree plot (tt) without heatmap
      ggsave(out_file_path, plot = tt, width = width, height = height,
             units = units_out, limitsize = FALSE)
      save_success <- TRUE
      cat(file=stderr(), paste0("  v88: Fallback 4 succeeded (saved tree only, no heatmap)\n"))
      cat(file=stderr(), paste0("  v88: WARNING - Heatmap could not be rendered\n"))
    }, error = function(e5) {
      cat(file=stderr(), paste0("  v88: Fallback 4 failed: ", e5$message, "\n"))
      stop(e5)  # Re-throw if all fallbacks fail
    })
  }

  if (!save_success) {
    cat(file=stderr(), paste0("  v88: All save attempts failed\n"))
    stop("Could not save plot after multiple attempts")
  }

  return(p)
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
  dashboardHeader(title = "Lineage Tree Plotter v158"),
  
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
    ),

    # Configuration file input
    fileInput("yaml_config", "Import YAML Configuration", 
              accept = c(".yaml", ".yml"))
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
            title = "📌 Version Info",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tags$div(style = "background: #d4edda; padding: 15px; border-radius: 5px; border: 2px solid #28a745;",
                     tags$h4(style = "color: #155724; margin: 0;", "v158 Active!"),
                     tags$p(style = "margin: 10px 0 0 0; color: #155724;",
                            "New in v158:",
                            tags$ul(
                              tags$li("FIXED LEGEND VISIBILITY: Custom legends now properly visible using nested gtable approach"),
                              tags$li("Legends appear ABOVE the ggplot guide-box with proper sizing and positioning")
                            ),
                            "Previous fixes:",
                            tags$ul(
                              tags$li("v155: Gtable-based legend approach (alignment with ggplot legends)"),
                              tags$li("v154: Title fontsize matches ggplot legend titles"),
                              tags$li("v153: Plot distortion fixed")
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
              selectInput("individual_value", "Select Individual", 
                          choices = NULL, 
                          selected = NULL)
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
            checkboxInput("trim_tips", "Trim Tips", value = TRUE),
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

            # Legend Visibility box
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
              checkboxInput("legend_show_heatmap", "Heatmap Legends", value = TRUE)
            ),

            # Font Sizes box
            box(
              title = NULL,
              status = "warning",
              solidHeader = FALSE,
              width = 12,
              tags$h4(icon("text-height"), " Font Sizes", style = "margin-top: 0;"),
              sliderInput("legend_title_size", "Legend Title Size",
                          min = 4, max = 48, value = 12, step = 1),  # v122: Increased max from 24 to 48
              sliderInput("legend_text_size", "Legend Text Size",
                          min = 2, max = 36, value = 10, step = 1)   # v122: Increased max from 18 to 36
            ),

            # Symbol Settings box
            box(
              title = NULL,
              status = "success",
              solidHeader = FALSE,
              width = 12,
              tags$h4(icon("square"), " Symbol Settings", style = "margin-top: 0;"),
              sliderInput("legend_key_size", "Legend Key Size (symbols)",
                          min = 0.1, max = 5, value = 1, step = 0.1),  # v122: Increased max from 2 to 5
              sliderInput("legend_spacing", "Legend Spacing",
                          min = 0.05, max = 3, value = 0.3, step = 0.05)  # v122: Increased max from 1 to 3
            ),

            # v133: Highlight Legend Settings box
            box(
              title = NULL,
              status = "warning",
              solidHeader = FALSE,
              width = 12,
              collapsible = TRUE,
              collapsed = TRUE,
              tags$h4(icon("highlighter"), " Highlight Legend Settings", style = "margin-top: 0;"),
              tags$p(class = "text-muted", "Control highlight legend appearance (positioned below main legend):"),
              fluidRow(
                column(6,
                  # v138: Changed default from 0 to -2 to keep legend within image boundaries
                  sliderInput("highlight_legend_x_offset", "X Offset (Left/Right)",
                              min = -5, max = 5, value = -2, step = 0.1)
                ),
                column(6,
                  # v142: Changed default from 0 to -3 to position below main legend
                  sliderInput("highlight_legend_y_offset", "Y Offset (Up/Down)",
                              min = -10, max = 10, value = -3, step = 0.5)
                )
              ),
              fluidRow(
                column(6,
                  sliderInput("highlight_legend_title_size", "Title Size",
                              min = 0.2, max = 10, value = 0.3, step = 0.1)  # v144: smaller default (was 0.5)
                ),
                column(6,
                  sliderInput("highlight_legend_text_size", "Label Text Size",
                              min = 0.3, max = 8, value = 1, step = 0.1)  # v142: smaller default (was 1.5)
                )
              ),
              fluidRow(
                column(6,
                  sliderInput("highlight_legend_title_gap", "Gap: Title to Labels",
                              min = -5, max = 10, value = 1, step = 0.25)
                ),
                column(6,
                  sliderInput("highlight_legend_label_gap", "Gap: Between Labels",
                              min = 0.01, max = 3, value = 0.5, step = 0.05)
                )
              )
            ),

            # v133: Bootstrap Legend Settings box
            box(
              title = NULL,
              status = "info",
              solidHeader = FALSE,
              width = 12,
              collapsible = TRUE,
              collapsed = TRUE,
              tags$h4(icon("chart-line"), " Bootstrap Legend Settings", style = "margin-top: 0;"),
              tags$p(class = "text-muted", "Control bootstrap legend appearance (positioned below highlight legend):"),
              fluidRow(
                column(6,
                  # v138: Changed default from 0 to -2 to keep legend within image boundaries
                  sliderInput("bootstrap_legend_x_offset", "X Offset (Left/Right)",
                              min = -5, max = 5, value = -2, step = 0.1)
                ),
                column(6,
                  # v144: Changed default from -6 to -8 to position further below highlight legend
                  sliderInput("bootstrap_legend_y_offset", "Y Offset (Up/Down)",
                              min = -15, max = 10, value = -8, step = 0.5)
                )
              ),
              fluidRow(
                column(6,
                  sliderInput("bootstrap_legend_title_size", "Title Size",
                              min = 0.1, max = 5, value = 0.8, step = 0.05)
                ),
                column(6,
                  sliderInput("bootstrap_legend_text_size", "Label Text Size",
                              min = 0.1, max = 4, value = 0.6, step = 0.05)
                )
              ),
              fluidRow(
                column(6,
                  # v143: New control to move bootstrap title further right
                  sliderInput("bootstrap_legend_title_x_offset", "Title X Offset (Left/Right)",
                              min = -5, max = 10, value = 2, step = 0.5)
                ),
                column(6,
                  sliderInput("bootstrap_legend_title_gap", "Gap: Title to Triangles",
                              min = -5, max = 10, value = 2, step = 0.25)  # v133: smaller default gap
                )
              ),
              fluidRow(
                column(6,
                  sliderInput("bootstrap_legend_label_gap", "Gap: Labels to Triangles",
                              min = -5, max = 10, value = 2, step = 0.25)
                )
              )
            ),

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
                icon("spinner"), " Processing..."
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
    # v121: Legend settings
    legend_settings = list(
      position = "right",
      show_classification = TRUE,
      show_highlight = TRUE,
      show_bootstrap = TRUE,
      show_heatmap = TRUE,
      title_size = 12,
      text_size = 10,
      key_size = 1,
      spacing = 0.3
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
    plot_scale_percent = 100
  )
  
  classification_loading <- reactiveVal(FALSE)

  # v107: Trigger for heatmap UI regeneration - only fires when we need to rebuild the UI
  # (add/remove/reorder heatmaps, or when CSV data changes)
  # This prevents the UI from rebuilding every time a slider or input changes
  heatmap_ui_trigger <- reactiveVal(0)

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
  
  # Bootstrap checkbox observer
  # Bootstrap checkbox observer
  observeEvent(input$show_bootstrap, {
    # v53: cat(file=stderr(), "\n===observeEvent show_bootstrap FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$show_bootstrap, "\n")
    # v53: cat(file=stderr(), "plot_ready:", values$plot_ready, "\n")
    
    req(values$plot_ready)
    
    # Force update YAML before generating
    update_yaml()
    
    # v53: cat(file=stderr(), "Bootstrap display in YAML:", 
    #     values$yaml_data$`visual definitions`$Bootstrap$display, "\n")
    
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Bootstrap format observer
  observeEvent(input$bootstrap_format, {
    # v53: cat(file=stderr(), "\n===observeEvent bootstrap_format FIRED===\n")
    req(values$plot_ready, input$show_bootstrap == TRUE)
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Bootstrap param observer
  observeEvent(input$bootstrap_param, {
    # v53: cat(file=stderr(), "\n===observeEvent bootstrap_param FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$bootstrap_param, "\n")
    req(values$plot_ready, input$show_bootstrap == TRUE)
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Bootstrap label size observer
  observeEvent(input$bootstrap_label_size, {
    # v53: cat(file=stderr(), "\n===observeEvent bootstrap_label_size FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$bootstrap_label_size, "\n")
    req(values$plot_ready, input$show_bootstrap == TRUE)
    generate_plot()
  }, ignoreInit = TRUE)

  # v112: Bootstrap position offset observer - makes slider immediately responsive
  observeEvent(input$man_boot_x_offset, {
    cat(file=stderr(), "\n===observeEvent man_boot_x_offset FIRED===\n")
    cat(file=stderr(), "New value:", input$man_boot_x_offset, "\n")
    req(values$plot_ready, input$show_bootstrap == TRUE)
    generate_plot()
  }, ignoreInit = TRUE)

  # Process YAML configuration when loaded
  observeEvent(input$yaml_config, {
    req(input$yaml_config)
    
    yaml_data <- parse_yaml_config(input$yaml_config$datapath)
    if ("error" %in% names(yaml_data)) {
      showNotification(yaml_data$error, type = "error")
      return()
    }
    
    # Load the configuration into the app
    updateTextInput(session, "individual_name", value = yaml_data$`Individual general definitions`$Individual)
    
    # Load tree file if available
    if (!is.null(yaml_data$`Individual general definitions`$`tree path`) && 
        file.exists(yaml_data$`Individual general definitions`$`tree path`[[1]])) {
      # Load tree file
      tree_file <- yaml_data$`Individual general definitions`$`tree path`[[1]]
      tree <- read.tree(tree_file)
      values$tree <- tree
      # v56b: Suppress harmless fortify warnings
      values$tree_data <- suppressWarnings(ggtree(tree))$data

      # Update tree summary
      output$tree_summary <- renderPrint({
        cat("Tree loaded from YAML configuration\n")
        cat("Number of tips:", length(tree$tip.label), "\n")
        cat("Number of internal nodes:", tree$Nnode, "\n")
      })
    }
    
    # Load CSV file if available
    if (!is.null(yaml_data$`Individual general definitions`$`mapping csv file`) && 
        file.exists(yaml_data$`Individual general definitions`$`mapping csv file`)) {
      csv_file <- yaml_data$`Individual general definitions`$`mapping csv file`
      csv_data <- read.csv(csv_file)
      values$csv_data <- csv_data
      # v107: Trigger heatmap UI regeneration when CSV data changes (new column choices)
      heatmap_ui_trigger(heatmap_ui_trigger() + 1)

      # Update ID column
      id_col <- yaml_data$`Mapping exl renaming titles`$`ID column`
      if (id_col %in% names(csv_data)) {
        updateSelectInput(session, "id_column", choices = names(csv_data), selected = id_col)
      } else {
        updateSelectInput(session, "id_column", choices = names(csv_data))
      }
      
      # Update Individual column
      individual_col <- yaml_data$`Individual general definitions`$`individual column`
      if (!is.null(individual_col) && individual_col %in% names(csv_data)) {
        updateSelectInput(session, "individual_column", choices = names(csv_data), selected = individual_col)
      } else {
        updateSelectInput(session, "individual_column", choices = names(csv_data))
      }
      
      # Update CSV summary
      output$csv_summary <- renderPrint({
        cat("CSV loaded from YAML configuration\n")
        cat("Number of rows:", nrow(csv_data), "\n")
        cat("Number of columns:", ncol(csv_data), "\n")
        cat("ID column:", id_col, "\n")
      })
    }
    
    # Load display settings
    if (!is.null(yaml_data$`visual definitions`)) {
      
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
        } else {
          updateCheckboxInput(session, "show_bootstrap", value = FALSE)
        }
      }
      
      # Load Rotation settings
      if (!is.null(yaml_data$`visual definitions`$rotation1)) {
        rotation1_display <- yaml_data$`visual definitions`$rotation1$display
        if (func.check.bin.val.from.conf(rotation1_display)) {
          updateCheckboxInput(session, "enable_rotation", value = TRUE)
          updateRadioButtons(session, "rotation_type", selected = "primary")
          
          # Load rotation classes
          if (!is.null(yaml_data$`visual definitions`$rotation1$according)) {
            values$rotation_settings$primary <- yaml_data$`visual definitions`$rotation1$according
          }
        }
      }
      
      if (!is.null(yaml_data$`visual definitions`$rotation2)) {
        rotation2_display <- yaml_data$`visual definitions`$rotation2$display
        if (func.check.bin.val.from.conf(rotation2_display)) {
          # If both rotations are enabled, default to primary first
          values$rotation_settings$secondary <- yaml_data$`visual definitions`$rotation2$according
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
        
        # Other font sizes could be loaded here
      }
      
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
    }
    
    # Force update of UI elements based on loaded data
    updateSelectInput(session, "classification_column", choices = names(values$csv_data), selected = character(0))
    updateSelectInput(session, "highlight_column", choices = names(values$csv_data), selected = character(0))
    # v55: heatmap_columns removed - now using individual column selects per heatmap
    # updateSelectizeInput(session, "heatmap_columns", choices = names(values$csv_data))
    
    # Generate initial plot based on loaded configuration
    generate_plot()
    
    showNotification("YAML configuration loaded successfully", type = "message")
  })
  
  # When CSV file is uploaded
  observeEvent(input$csv_file, {
    req(input$csv_file)
    
    csv_file <- input$csv_file
    
    # Show progress
    values$progress_message <- "[FILE]  Loading CSV file..."
    values$progress_visible <- TRUE
    
    # Read CSV file
    tryCatch({
      csv_data <- read.csv(csv_file$datapath)
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
      updateSelectInput(session, "individual_value", choices = NULL, selected = NULL)
      return()
    }
    
    # Check if column exists in CSV
    if (!(input$individual_column %in% names(values$csv_data))) {
      # v53: cat(file=stderr(), paste("Column", input$individual_column, "not found in CSV\n"))
      showNotification(paste("Column not found:", input$individual_column), type = "warning")
      updateSelectInput(session, "individual_value", choices = NULL, selected = NULL)
      return()
    }
    
    # Get unique values safely
    tryCatch({
      unique_values <- unique(values$csv_data[[input$individual_column]])
      unique_values <- unique_values[!is.na(unique_values)]
      
      # v53: cat(file=stderr(), paste("Found", length(unique_values), "unique values\n"))
      # v53: cat(file=stderr(), "Unique values:", paste(head(unique_values, 10), collapse=", "), "\n")
      
      # Update individual value dropdown with these unique values
      updateSelectInput(session, "individual_value", 
                        choices = unique_values,
                        selected = if(length(unique_values) > 0) unique_values[1] else NULL)
      
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
    req(input$tree_file)
    
    tree_file <- input$tree_file
    
    # Show progress
    values$progress_message <- "📂 Loading tree file..."
    values$progress_visible <- TRUE
    
    # Read tree file
    tryCatch({
      tree <- read.tree(tree_file$datapath)
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
          cat(file=stderr(), paste0("\n=== v56c: Adding ", length(heatmaps_to_use), " heatmap(s) to classification ", i, " ===\n"))

          class_item[[as.character(i)]]$heatmap_display <- list()

          for (j in seq_along(heatmaps_to_use)) {
            heatmap_entry <- heatmaps_to_use[[j]]

            # v56c: DEBUG
            cat(file=stderr(), paste0("  Processing heatmap ", j, ": ", heatmap_entry$title, "\n"))
            cat(file=stderr(), paste0("    Columns: ", paste(heatmap_entry$columns, collapse=", "), "\n"))

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

            # v118: Add guide line settings (was missing - this caused tip guide lines not to work!)
            heatmap_item[[as.character(j)]]$show_guides <- if (!is.null(heatmap_entry$show_guides) && heatmap_entry$show_guides) "yes" else "no"
            heatmap_item[[as.character(j)]]$guide_color1 <- if (!is.null(heatmap_entry$guide_color1)) heatmap_entry$guide_color1 else "#CCCCCC"
            heatmap_item[[as.character(j)]]$guide_color2 <- if (!is.null(heatmap_entry$guide_color2)) heatmap_entry$guide_color2 else "#EEEEEE"
            heatmap_item[[as.character(j)]]$guide_alpha <- if (!is.null(heatmap_entry$guide_alpha)) heatmap_entry$guide_alpha else 0.3
            heatmap_item[[as.character(j)]]$guide_width <- if (!is.null(heatmap_entry$guide_width)) heatmap_entry$guide_width else 0.5

            # v109: Add colnames_angle
            heatmap_item[[as.character(j)]]$colnames_angle <- if (!is.null(heatmap_entry$colnames_angle)) heatmap_entry$colnames_angle else 45

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

            # Add columns - format must match expected YAML structure
            # Each column entry needs to be a named list like list("1" = "column_name")
            if (!is.null(heatmap_entry$columns)) {
              for (k in seq_along(heatmap_entry$columns)) {
                column_entry <- list()
                column_entry[[as.character(k)]] <- heatmap_entry$columns[k]
                heatmap_item[[as.character(j)]]$according[[k]] <- column_entry
              }
            }

            class_item[[as.character(i)]]$heatmap_display[[j]] <- heatmap_item
          }
        }
        
        # ========================================================
        # === ADD HIGHLIGHT TO THIS CLASSIFICATION ===
        # *** KEY FIX: MOVED HERE - BEFORE adding class_item to YAML ***
        # ========================================================
        highlight_to_apply <- NULL

        # v131: DEBUG - trace highlight decision
        cat(file=stderr(), paste0("\n=== v131: HIGHLIGHT DECISION (classification ", i, ") ===\n"))
        cat(file=stderr(), paste0("  temp_highlight_preview is NULL: ", is.null(values$temp_highlight_preview), "\n"))
        cat(file=stderr(), paste0("  active_highlight_index: ", values$active_highlight_index, "\n"))
        cat(file=stderr(), paste0("  highlights length: ", length(values$highlights), "\n"))

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
          cat(file=stderr(), "  NO highlight source found (will set display='no')\n")
        }

        # Apply highlight to THIS classification
        if (!is.null(highlight_to_apply)) {
          cat(file=stderr(), paste0("  APPLYING highlight with ", length(highlight_to_apply$items), " items (display='yes')\n"))
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
        cat(file=stderr(), paste0("\n=== v56c: Adding ", length(values$heatmaps), " heatmap(s) to DEFAULT classification ===\n"))

        default_classification[["1"]]$heatmap_display <- list()

        for (j in seq_along(values$heatmaps)) {
          heatmap_entry <- values$heatmaps[[j]]

          # v56c: DEBUG
          cat(file=stderr(), paste0("  Processing heatmap ", j, ": ", heatmap_entry$title, "\n"))
          cat(file=stderr(), paste0("    Columns: ", paste(heatmap_entry$columns, collapse=", "), "\n"))

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
              cat(file=stderr(), paste0("    v104: Storing ", length(custom_colors_as_list), " custom colors with names: ", paste(names(custom_colors_as_list), collapse=", "), "\n"))
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

          # v120: Add guide line settings (were missing in default classification path - caused tip guide lines not to work!)
          heatmap_item[[as.character(j)]]$show_guides <- if (!is.null(heatmap_entry$show_guides) && heatmap_entry$show_guides) "yes" else "no"
          heatmap_item[[as.character(j)]]$guide_color1 <- if (!is.null(heatmap_entry$guide_color1)) heatmap_entry$guide_color1 else "#CCCCCC"
          heatmap_item[[as.character(j)]]$guide_color2 <- if (!is.null(heatmap_entry$guide_color2)) heatmap_entry$guide_color2 else "#EEEEEE"
          heatmap_item[[as.character(j)]]$guide_alpha <- if (!is.null(heatmap_entry$guide_alpha)) heatmap_entry$guide_alpha else 0.3
          heatmap_item[[as.character(j)]]$guide_width <- if (!is.null(heatmap_entry$guide_width)) heatmap_entry$guide_width else 0.5

          # v109: Add colnames_angle
          heatmap_item[[as.character(j)]]$colnames_angle <- if (!is.null(heatmap_entry$colnames_angle)) heatmap_entry$colnames_angle else 45

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

          # Add columns - format must match expected YAML structure
          if (!is.null(heatmap_entry$columns)) {
            for (k in seq_along(heatmap_entry$columns)) {
              column_entry <- list()
              column_entry[[as.character(k)]] <- heatmap_entry$columns[k]
              heatmap_item[[as.character(j)]]$according[[k]] <- column_entry
            }
          }

          default_classification[["1"]]$heatmap_display[[j]] <- heatmap_item
        }
      }

      # v129: Add highlighting support to default classification (was missing!)
      # v130: Fixed - removed values$preview_highlight_active check (was never set)
      # Determine which highlight to apply (same logic as custom classification path)
      highlight_to_apply <- NULL

      # v131: DEBUG - trace highlight decision for default classification
      cat(file=stderr(), paste0("\n=== v131: HIGHLIGHT DECISION (DEFAULT classification) ===\n"))
      cat(file=stderr(), paste0("  temp_highlight_preview is NULL: ", is.null(values$temp_highlight_preview), "\n"))
      cat(file=stderr(), paste0("  active_highlight_index: ", values$active_highlight_index, "\n"))
      cat(file=stderr(), paste0("  highlights length: ", length(values$highlights), "\n"))

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
        cat(file=stderr(), paste0("\n=== v129: Adding highlight to DEFAULT classification ===\n"))

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
        cat(file=stderr(), paste0("  Highlight added with ", length(highlight_to_apply$items), " items\n"))
      } else {
        # No highlight - disable it
        cat(file=stderr(), "  NO highlight to apply (setting display='no')\n")
        default_classification[["1"]]$highlight <- list(display = "no")
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
    
    # Attempt to match tree labels with CSV IDs
    match_result <- match_tree_ids_with_csv(tree_labels, csv_ids)
    values$id_match <- match_result
    
    
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
    
    # Create filtered CSV data with only the matched rows
    if (match_result$summary$exact_matches > 0 || 
        match_result$summary$numeric_matches > 0 || 
        match_result$summary$prefix_suffix_matches > 0) {
      
      # Get all matched CSV IDs
      matched_ids <- c()
      for (mapping in match_result$mapping) {
        matched_ids <- c(matched_ids, mapping)
      }
      matched_ids <- unique(matched_ids)
      
      # Filter CSV to only include matched rows
      values$filtered_csv <- filtered_csv[filtered_csv[[input$id_column]] %in% matched_ids, ]
      
      # v53: print("Filtered CSV rows:")
      # v53: print(nrow(values$filtered_csv))
      # v53: print("Filtered CSV columns:")
      # v53: print(names(values$filtered_csv))
    } else {
      values$filtered_csv <- NULL
    }
    
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
      
      # Create temp CSV file IMMEDIATELY after filtering, BEFORE generate_plot
      if (!is.null(values$filtered_csv) && nrow(values$filtered_csv) > 0) {
        # v53: cat(file=stderr(), "Creating temp CSV file...\n")
        values$temp_csv_path <- tempfile(fileext = ".csv")
        write.csv(values$filtered_csv, values$temp_csv_path, row.names = FALSE)
        # v53: cat(file=stderr(), "Created temp CSV at:", values$temp_csv_path, "\n")
      } else {
        # v53: cat(file=stderr(), "WARNING: No filtered CSV data available for temp file\n")
      }
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
      
      # v53: cat(file=stderr(), "About to call generate_plot()\n")
      # Generate initial plot
      generate_plot()
      # v53: cat(file=stderr(), "Finished generate_plot()\n")
      
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
  })
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
  })
  
  
  
  
  
  
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
    
    # v53: cat(file=stderr(), "Ã¢Å“â€œ All checks passed, collecting settings...\n")
    
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
    
    # v53: cat(file=stderr(), "Ã¢Å“â€œ Collected", length(highlight_items), "items\n")
    
    # Store as temporary preview
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
    cat(file=stderr(), "\n=== v131: HIGHLIGHT BUTTON - BEFORE generate_plot() ===\n")
    cat(file=stderr(), paste0("  temp_highlight_preview is NULL: ", is.null(values$temp_highlight_preview), "\n"))
    if (!is.null(values$temp_highlight_preview)) {
      cat(file=stderr(), paste0("  temp_highlight_preview column: ", values$temp_highlight_preview$column, "\n"))
      cat(file=stderr(), paste0("  temp_highlight_preview items: ", length(values$temp_highlight_preview$items), "\n"))
    }
    cat(file=stderr(), "  CALLING generate_plot() NOW...\n")
    values$debug_trace_id <- "HIGHLIGHT_BUTTON_PREVIEW"
    generate_plot()
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
  })
  
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
        ),

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
        ),

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
                             min = 0, max = 90, value = 45, step = 15)
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
                             min = 0.1, max = 3.0,
                             value = if (!is.null(cfg$height)) cfg$height else 0.8,
                             step = 0.1)
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
            column(4)
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
    
    # v56: Add new empty config with columns (plural) for multiple column support
    # v107: Added distance and height initialization to prevent reactive loops
    new_config <- list(
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
      midpoint = 0
    )
    
    values$heatmap_configs <- c(values$heatmap_configs, list(new_config))
    # v107: Trigger UI regeneration when heatmap is added
    heatmap_ui_trigger(heatmap_ui_trigger() + 1)
    showNotification(paste("Heatmap", length(values$heatmap_configs), "added"), type = "message")
  })

  # Generic observer for heatmap removal buttons
  observe({
    lapply(1:6, function(i) {
      observeEvent(input[[paste0("heatmap_remove_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs <- values$heatmap_configs[-i]
          # v107: Trigger UI regeneration when heatmap is removed
          heatmap_ui_trigger(heatmap_ui_trigger() + 1)
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
          # v107: Trigger UI regeneration when heatmap is moved
          heatmap_ui_trigger(heatmap_ui_trigger() + 1)
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
          # v107: Trigger UI regeneration when heatmap is moved
          heatmap_ui_trigger(heatmap_ui_trigger() + 1)
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
      observeEvent(input[[paste0("heatmap_discrete_palette_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$discrete_palette <- input[[paste0("heatmap_discrete_palette_", i)]]
        }
      }, ignoreInit = TRUE)
      
      # Continuous palette change
      observeEvent(input[[paste0("heatmap_cont_palette_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$cont_palette <- input[[paste0("heatmap_cont_palette_", i)]]
        }
      }, ignoreInit = TRUE)
      
      # Color changes
      observeEvent(input[[paste0("heatmap_low_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$low_color <- input[[paste0("heatmap_low_color_", i)]]
        }
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("heatmap_high_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$high_color <- input[[paste0("heatmap_high_color_", i)]]
        }
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("heatmap_mid_color_", i)]], {
        if (i <= length(values$heatmap_configs)) {
          values$heatmap_configs[[i]]$mid_color <- input[[paste0("heatmap_mid_color_", i)]]
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
        cat(file=stderr(), paste0("  v127 DETECTED TYPE DISPLAY: column=", first_col,
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
                     colourInput(paste0("heatmap_low_color_", i), "Low Color",
                                 value = if (!is.null(cfg$low_color)) cfg$low_color else "#FFFFCC")
              ),
              column(4,
                     colourInput(paste0("heatmap_high_color_", i), "High Color",
                                 value = if (!is.null(cfg$high_color)) cfg$high_color else "#006837")
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
                       colourInput(paste0("heatmap_mid_color_", i), "Mid Color",
                                   value = if (!is.null(cfg$mid_color)) cfg$mid_color else "#FFFF99")
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
                     colourInput(paste0("heatmap_", i, "_cont_na_color"), "NA Color",
                                 value = if (!is.null(cfg$na_color)) cfg$na_color else "#BEBEBE",
                                 showColour = "background")
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
            cat(file=stderr(), paste0("v125: Filtered unique values from ", nrow(values$csv_data), " to ", nrow(filtered_data), " rows (tree tips)\n"))
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

        # v70: Generate color pickers for each value WITH dropdown menu
        color_pickers <- lapply(seq_along(unique_vals), function(j) {
          val <- as.character(unique_vals[j])

          # v69: Check if there's an existing custom color for this value
          existing_color <- isolate(input[[paste0("heatmap_", i, "_color_", j)]])
          color_to_use <- if (!is.null(existing_color)) existing_color else default_colors[j]

          fluidRow(
            style = "margin-bottom: 3px;",
            column(4, tags$label(val, style = "padding-top: 5px; font-weight: normal; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;", title = val)),
            column(4,
                   colourInput(paste0("heatmap_", i, "_color_", j), NULL,
                               value = color_to_use, showColour = "background")
            ),
            column(4,
                   selectInput(paste0("heatmap_", i, "_color_name_", j), NULL,
                               choices = heat_r_colors, selected = "")
            )
          )
        })

        # v70: NA color picker (always shown at the end)
        existing_na_color <- isolate(input[[paste0("heatmap_", i, "_na_color")]])
        na_color_to_use <- if (!is.null(existing_na_color)) existing_na_color else "white"

        na_color_row <- fluidRow(
          style = "margin-bottom: 3px; background-color: #f8f8f8; padding: 5px; border-radius: 3px; margin-top: 10px;",
          column(4, tags$label("NA / Missing", style = "padding-top: 5px; font-weight: bold; font-style: italic;")),
          column(4,
                 colourInput(paste0("heatmap_", i, "_na_color"), NULL,
                             value = na_color_to_use, showColour = "background")
          ),
          column(4,
                 selectInput(paste0("heatmap_", i, "_na_color_name"), NULL,
                             choices = heat_r_colors, selected = "")
          )
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

        # Update all color pickers
        for (j in seq_along(unique_vals)) {
          updateColourInput(session, paste0("heatmap_", i, "_color_", j), value = new_colors[j])
        }

        showNotification(paste("Applied", palette_name, "palette to", n_vals, "values"), type = "message")
      }, ignoreInit = TRUE)
    })
  })
  
  # Apply heatmaps button
  observeEvent(input$apply_heatmaps, {
    # Convert heatmap_configs to the format expected by the plotting function
    if (length(values$heatmap_configs) == 0) {
      showNotification("No heatmaps configured", type = "warning")
      return()
    }

    # v56a: Build heatmaps list from configs with multiple column support
    # Read directly from inputs to ensure we get current values (fixes ignoreInit issue)
    heatmaps_list <- lapply(seq_along(values$heatmap_configs), function(i) {
      cfg <- values$heatmap_configs[[i]]

      # v56a: Read columns directly from input (fixes issue where ignoreInit=TRUE misses initial selection)
      current_columns <- input[[paste0("heatmap_columns_", i)]]
      if (is.null(current_columns) || length(current_columns) == 0) {
        return(NULL)
      }

      # Update config with current columns
      cfg$columns <- current_columns

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
      cat(file=stderr(), paste0("  v122 AUTO-DETECT: auto_type=", current_auto_type,
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
        cat(file=stderr(), paste0("  v117 AUTO-DETECT: has_decimal_in_string=", has_decimal_in_string, "\n"))

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
        cat(file=stderr(), paste0("  v117 AUTO-DETECT: column=", first_col,
                                   ", is_numeric=", is_numeric,
                                   ", originally_numeric=", originally_numeric,
                                   ", converted=", converted_to_numeric,
                                   ", unique_vals=", unique_vals, "\n"))

        # v117: Better heuristic - decimals ALWAYS mean continuous
        if (!is_numeric) {
          # Non-numeric data is always discrete
          actual_type <- "discrete"
          cat(file=stderr(), paste0("  v117 AUTO-DETECT: Result=discrete (non-numeric)\n"))
        } else {
          # v117: Check for decimal values with tolerance for floating-point precision
          epsilon <- 1e-6
          has_decimals <- any(abs(col_data_clean - floor(col_data_clean)) > epsilon, na.rm = TRUE)

          # v117: Also use the string-based detection result
          if (!has_decimals && has_decimal_in_string) {
            has_decimals <- TRUE
            cat(file=stderr(), paste0("  v117 AUTO-DETECT: decimal detected via string check\n"))
          }

          cat(file=stderr(), paste0("  v120 AUTO-DETECT: has_decimals=", has_decimals, "\n"))

          # v120: Harmonized logic with UI section for consistent detection
          # Get value range for range-based checks
          val_range <- if (length(col_data_clean) > 0) diff(range(col_data_clean, na.rm = TRUE)) else 0

          # v120: Decimals ALWAYS mean continuous (measurements, percentages, etc.)
          if (has_decimals) {
            actual_type <- "continuous"
            cat(file=stderr(), paste0("  v120 AUTO-DETECT: Result=continuous (has decimals)\n"))
          } else if (originally_numeric) {
            # v120: Originally numeric without decimals - match UI section logic
            is_boolean_like <- unique_vals <= 2 && val_range <= 1
            is_small_categorical <- unique_vals <= 3 && val_range <= 2
            if (is_boolean_like || is_small_categorical) {
              actual_type <- "discrete"
              cat(file=stderr(), paste0("  v120 AUTO-DETECT: Result=discrete (boolean-like or small categorical)\n"))
            } else {
              actual_type <- "continuous"
              cat(file=stderr(), paste0("  v120 AUTO-DETECT: Result=continuous (originally numeric, many values)\n"))
            }
          } else {
            # v123: Converted from character without decimals - more lenient for numeric-like data
            # If >=5 unique values OR range >=5 → treat as continuous
            # This catches integer sequences like 1,2,3,4,5 which are typically counts/scores
            if (unique_vals >= 5 || val_range >= 5) {
              actual_type <- "continuous"
              cat(file=stderr(), paste0("  v123 AUTO-DETECT: Result=continuous (>=5 unique or range>=5)\n"))
            } else {
              actual_type <- "discrete"
              cat(file=stderr(), paste0("  v123 AUTO-DETECT: Result=discrete (few unique, narrow range)\n"))
            }
          }
        }
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
        # v116: Guide lines
        show_guides = if (!is.null(show_guides)) show_guides else FALSE,
        guide_color1 = if (!is.null(guide_color1)) guide_color1 else "#CCCCCC",
        guide_color2 = if (!is.null(guide_color2)) guide_color2 else "#EEEEEE",
        guide_alpha = if (!is.null(guide_alpha)) guide_alpha else 0.3,
        guide_width = if (!is.null(guide_width)) guide_width else 0.5,
        show_row_labels = if (!is.null(show_row_labels)) show_row_labels else FALSE,
        row_label_source = if (!is.null(row_label_source)) row_label_source else "colnames",
        row_label_font_size = if (!is.null(row_label_font_size)) row_label_font_size else 2.5,
        row_label_offset = if (!is.null(row_label_offset)) row_label_offset else 1.0,  # v111
        row_label_align = if (!is.null(row_label_align)) row_label_align else "left",  # v111
        custom_row_labels = if (!is.null(custom_row_labels)) custom_row_labels else "",
        label_mapping = label_mapping  # v108: Per-column label mapping
      )
      
      if (actual_type == "discrete") {
        # v69: Get palette from current input (not stale cfg)
        current_palette <- input[[paste0("heatmap_discrete_palette_", i)]]
        heatmap_entry$color_scheme <- if (!is.null(current_palette)) current_palette else "Set1"

        # v69: Collect custom colors if they've been set
        if (!is.null(values$csv_data) && first_col %in% names(values$csv_data)) {
          unique_vals <- sort(unique(na.omit(values$csv_data[[first_col]])))
          n_vals <- length(unique_vals)

          if (n_vals > 0 && n_vals <= 30) {
            custom_colors <- c()
            has_custom_colors <- FALSE
            for (j in seq_along(unique_vals)) {
              color_input <- input[[paste0("heatmap_", i, "_color_", j)]]
              if (!is.null(color_input)) {
                custom_colors[as.character(unique_vals[j])] <- color_input
                has_custom_colors <- TRUE
              }
            }
            if (has_custom_colors) {
              heatmap_entry$custom_colors <- custom_colors
              heatmap_entry$man_define_colors <- TRUE
            }
          }

          # v70: Get NA color from input (default to white if not set)
          na_color_input <- input[[paste0("heatmap_", i, "_na_color")]]
          heatmap_entry$na_color <- if (!is.null(na_color_input)) na_color_input else "white"
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
    cat(file=stderr(), "\n=== v56c HEATMAP APPLY DEBUG ===\n")
    cat(file=stderr(), "Number of heatmaps:", length(heatmaps_list), "\n")
    for (h_idx in seq_along(heatmaps_list)) {
      hm <- heatmaps_list[[h_idx]]
      cat(file=stderr(), paste0("  Heatmap ", h_idx, ":\n"))
      cat(file=stderr(), paste0("    Title: ", hm$title, "\n"))
      cat(file=stderr(), paste0("    Columns: ", paste(hm$columns, collapse=", "), "\n"))
      cat(file=stderr(), paste0("    Is discrete: ", hm$is_discrete, "\n"))
    }
    cat(file=stderr(), "================================\n\n")

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
  observeEvent(input$apply_legend_settings, {
    cat(file=stderr(), "\n=== v133: APPLYING LEGEND SETTINGS ===\n")

    # Update legend settings in reactive values
    values$legend_settings <- list(
      position = input$legend_position,
      show_classification = input$legend_show_classification,
      show_highlight = input$legend_show_highlight,
      show_bootstrap = input$legend_show_bootstrap,
      show_heatmap = input$legend_show_heatmap,
      title_size = input$legend_title_size,
      text_size = input$legend_text_size,
      key_size = input$legend_key_size,
      spacing = input$legend_spacing,
      # v133: Highlight legend settings
      highlight_x_offset = input$highlight_legend_x_offset,
      highlight_y_offset = input$highlight_legend_y_offset,
      highlight_title_size = input$highlight_legend_title_size,
      highlight_text_size = input$highlight_legend_text_size,
      highlight_title_gap = input$highlight_legend_title_gap,
      highlight_label_gap = input$highlight_legend_label_gap,
      # v133: Bootstrap legend settings
      # v143: Added bootstrap_title_x_offset for moving title to the right
      bootstrap_x_offset = input$bootstrap_legend_x_offset,
      bootstrap_y_offset = input$bootstrap_legend_y_offset,
      bootstrap_title_x_offset = input$bootstrap_legend_title_x_offset,  # v143
      bootstrap_title_size = input$bootstrap_legend_title_size,
      bootstrap_text_size = input$bootstrap_legend_text_size,
      bootstrap_title_gap = input$bootstrap_legend_title_gap,
      bootstrap_label_gap = input$bootstrap_legend_label_gap
    )

    cat(file=stderr(), paste0("  Position: ", input$legend_position, "\n"))
    cat(file=stderr(), paste0("  Show classification: ", input$legend_show_classification, "\n"))
    cat(file=stderr(), paste0("  Show highlight: ", input$legend_show_highlight, "\n"))
    cat(file=stderr(), paste0("  Show bootstrap: ", input$legend_show_bootstrap, "\n"))
    cat(file=stderr(), paste0("  Show heatmap: ", input$legend_show_heatmap, "\n"))
    cat(file=stderr(), paste0("  Title size: ", input$legend_title_size, "\n"))
    cat(file=stderr(), paste0("  Text size: ", input$legend_text_size, "\n"))
    cat(file=stderr(), paste0("  v133: Highlight legend - x_offset: ", input$highlight_legend_x_offset,
                               ", y_offset: ", input$highlight_legend_y_offset, "\n"))
    cat(file=stderr(), paste0("  v133: Bootstrap legend - x_offset: ", input$bootstrap_legend_x_offset,
                               ", y_offset: ", input$bootstrap_legend_y_offset, "\n"))
    cat(file=stderr(), "======================================\n\n")

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
  observeEvent(input$page_orientation, {
    current_width <- isolate(input$output_width)
    current_height <- isolate(input$output_height)

    cat(file=stderr(), paste0("\n=== v141: PAGE ORIENTATION CHANGED ===\n"))
    cat(file=stderr(), paste0("  Orientation: ", input$page_orientation, "\n"))
    cat(file=stderr(), paste0("  Current width: ", current_width, ", height: ", current_height, "\n"))

    if (!is.null(current_width) && !is.null(current_height)) {
      if (input$page_orientation == "landscape") {
        # Landscape: width should be > height
        if (current_width < current_height) {
          cat(file=stderr(), paste0("  Swapping to landscape: width=", current_height, ", height=", current_width, "\n"))
          updateNumericInput(session, "output_width", value = current_height)
          updateNumericInput(session, "output_height", value = current_width)
        } else {
          cat(file=stderr(), "  Already in landscape orientation (width >= height)\n")
        }
      } else {
        # Portrait: height should be > width
        if (current_width > current_height) {
          cat(file=stderr(), paste0("  Swapping to portrait: width=", current_height, ", height=", current_width, "\n"))
          updateNumericInput(session, "output_width", value = current_height)
          updateNumericInput(session, "output_height", value = current_width)
        } else {
          cat(file=stderr(), "  Already in portrait orientation (height >= width)\n")
        }
      }
    } else {
      cat(file=stderr(), "  WARNING: width or height is NULL, cannot swap\n")
    }
    cat(file=stderr(), "=== END PAGE ORIENTATION ===\n")

    # v143: Regenerate plot after a small delay to ensure new dimensions are available
    # updateNumericInput doesn't update input$ values immediately - they need a flush cycle
    if (!is.null(values$tree_data)) {
      shinyjs::delay(100, {
        generate_plot()
        cat(file=stderr(), "  v143: Plot regenerated for orientation preview (after delay)\n")
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
    req(input$enable_highlight, input$highlight_column, input$highlight_values, 
        length(values$classifications) > 0)
    
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
  observeEvent(input$tip_font_size, {
    # v53: cat(file=stderr(), "\n===observeEvent tip_font_size FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$tip_font_size, "\n")
    # v53: cat(file=stderr(), "plot_ready:", values$plot_ready, "\n")
    req(values$plot_ready)  # Only if plot has been generated at least once
    # v53: cat(file=stderr(), "Calling generate_plot()...\n")
    generate_plot()
    # v53: cat(file=stderr(), "===observer complete===\n\n")
  }, ignoreInit = TRUE)
  
  # Edge width
  observeEvent(input$edge_width, {
    # v53: cat(file=stderr(), "\n===observeEvent edge_width FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$edge_width, "\n")
    req(values$plot_ready)
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Tip length
  observeEvent(input$tip_length, {
    # v53: cat(file=stderr(), "\n===observeEvent tip_length FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$tip_length, "\n")
    req(values$plot_ready)
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Trim tips checkbox
  observeEvent(input$trim_tips, {
    # v53: cat(file=stderr(), "\n===observeEvent trim_tips FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$trim_tips, "\n")
    req(values$plot_ready)
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Display node numbers
  observeEvent(input$display_node_numbers, {
    # v53: cat(file=stderr(), "\n===observeEvent display_node_numbers FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$display_node_numbers, "\n")
    req(values$plot_ready)
    generate_plot()
  }, ignoreInit = TRUE)
  
  # Ladderize
  observeEvent(input$ladderize, {
    # v53: cat(file=stderr(), "\n===observeEvent ladderize FIRED===\n")
    # v53: cat(file=stderr(), "New value:", input$ladderize, "\n")
    req(values$plot_ready)
    generate_plot()
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
  observeEvent(input$node_number_font_size, {
    req(values$plot_ready)
    req(input$display_node_numbers == TRUE)  # Only regenerate if node numbers are displayed
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â§ Font size slider changed to:", input$node_number_font_size, "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â§ Regenerating plot...\n")
    
    generate_plot()
  }, ignoreInit = TRUE)  # Don't trigger on initialization
  
  # === NEW: Reactive observer for display node numbers checkbox ===
  observeEvent(input$display_node_numbers, {
    req(values$plot_ready)
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â§ Display node numbers changed to:", input$display_node_numbers, "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â§ Regenerating plot...\n")
    
    generate_plot()
  }, ignoreInit = TRUE)
  
  # === NEW: Reactive observer for highlight checkbox ===
  observeEvent(input$highlight_selected_nodes, {
    req(values$plot_ready)
    req(input$nodes_to_rotate)  # Only if nodes are selected
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Â§ Highlight checkbox changed to:", input$highlight_selected_nodes, "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Â§ Regenerating plot...\n")
    
    generate_plot()
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
  # Generate plot based on current settings
  # Generate plot based on current settings
  generate_plot <- function() {
    
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Âµ === generate_plot() ENTRY POINT ===\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ classification_loading():", classification_loading(), "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ values$plot_generating:", values$plot_generating, "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ values$tree is NULL:", is.null(values$tree), "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ values$csv_data is NULL:", is.null(values$csv_data), "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ values$temp_csv_path:", values$temp_csv_path, "\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ File exists:", file.exists(values$temp_csv_path), "\n\n")
    
    
    # v53: cat(file=stderr(), "ðŸ”µ TRACE ID:", if (!is.null(values$debug_trace_id)) values$debug_trace_id else "UNKNOWN", "\n")
    # v53: cat(file=stderr(), "\nÃ°Å¸â€Âµ === generate_plot() ENTRY ===\n")
    # v53: cat(file=stderr(), "Ã°Å¸â€Âµ temp_highlight_preview is NULL:", is.null(values$temp_highlight_preview), "\n")
    if (!is.null(values$temp_highlight_preview)) {
      # v53: cat(file=stderr(), "Ã°Å¸â€Âµ temp_highlight_preview has", length(values$temp_highlight_preview$items), "items\n")
    }
    # v53: cat(file=stderr(), "\n")
    
    if (classification_loading()) {
      # v53: cat(file=stderr(), "Ã¢ÂÂ¸Ã¯Â¸Â Skipping plot generation - classification UI loading\n")
      return(NULL)
    }
    
    # Don't generate if already generating (recursion guard)
    if (!is.null(values$plot_generating) && values$plot_generating) {
      # v53: cat(file=stderr(), "Ã¢Å¡Â Ã¯Â¸Â WARNING: Plot already generating. Preventing recursion.\n")
      return(NULL)
    }
    
    #if (classification_loading()) {
    #  cat(file=stderr(), "Ã¢ÂÂ¸Ã¯Â¸Â Skipping plot generation - classification UI loading\n")
    #  return(NULL)
    #}
    
    if (classification_loading()) {
      # v53: cat(file=stderr(), "Ã¢ÂÂ¸Ã¯Â¸Â Skipping plot generation - classification UI loading\n")
      return(NULL)
    }
    
    # Don't generate if already generating (recursion guard)
    if (!is.null(values$plot_generating) && values$plot_generating) {
      # v53: cat(file=stderr(), "Ã¢Å¡Â Ã¯Â¸Â WARNING: Plot already generating. Preventing recursion.\n")
      return(NULL)
    }
    

    # v53: cat(file=stderr(), "\n=== generate_plot() called ===\n")

    # v57: Show processing status IMMEDIATELY via shinyjs (before R blocks)
    show_status_processing()

    # Set generating flag for status indicator
    values$plot_generating <- TRUE
    values$plot_ready <- FALSE  # Mark as not ready while generating
    #shiny::invalidateLater(0, session)
    
    # Small delay to ensure progress shows
    # Give UI time to update
    Sys.sleep(0.3)
    
    values$progress_message <- "🎨 Generating your beautiful plot..."
    values$progress_visible <- TRUE
    
    # Another small delay to ensure progress bar shows
    Sys.sleep(0.2)
    
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
    
    # Filter CSV to only include matched rows
    matched_ids <- c()
    for (mapping in values$id_match$mapping) {
      matched_ids <- c(matched_ids, mapping)
    }
    matched_ids <- unique(matched_ids)
    
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
    
    # DEBUG: Print classification structure
    if (!is.null(values$yaml_data$`visual definitions`$classification)) {
      # v53: cat(file=stderr(), "\n=== CLASSIFICATION YAML STRUCTURE ===\n")
      # v53: cat(file=stderr(), "Number of classifications:", length(values$yaml_data$`visual definitions`$classification), "\n")
      for (i in seq_along(values$yaml_data$`visual definitions`$classification)) {
        class_item <- values$yaml_data$`visual definitions`$classification[[i]]
        # v53: cat(file=stderr(), paste0("Classification ", i, ":\n"))
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
      # v53: cat(file=stderr(), "=====================================\n\n")
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
      cat(file=stderr(), "\n=== v57 HEATMAP YAML DEBUG ===\n")
      for (ci in seq_along(yaml_data_modified$`visual definitions`$classification)) {
        class_entry <- yaml_data_modified$`visual definitions`$classification[[ci]]
        class_key <- names(class_entry)[1]
        cat(file=stderr(), paste0("Classification ", ci, " (key=", class_key, "):\n"))
        cat(file=stderr(), paste0("  Has heatmap_display: ", !is.null(class_entry[[class_key]]$heatmap_display), "\n"))
        if (!is.null(class_entry[[class_key]]$heatmap_display)) {
          cat(file=stderr(), paste0("  Number of heatmaps: ", length(class_entry[[class_key]]$heatmap_display), "\n"))
          for (hi in seq_along(class_entry[[class_key]]$heatmap_display)) {
            hm <- class_entry[[class_key]]$heatmap_display[[hi]]
            hm_key <- names(hm)[1]
            cat(file=stderr(), paste0("    Heatmap ", hi, " (key=", hm_key, "):\n"))
            cat(file=stderr(), paste0("      Title: ", hm[[hm_key]]$title, "\n"))
            cat(file=stderr(), paste0("      Display: ", hm[[hm_key]]$display, "\n"))
            cat(file=stderr(), paste0("      According length: ", length(hm[[hm_key]]$according), "\n"))
            if (length(hm[[hm_key]]$according) > 0) {
              for (ai in seq_along(hm[[hm_key]]$according)) {
                acc <- hm[[hm_key]]$according[[ai]]
                acc_key <- names(acc)[1]
                cat(file=stderr(), paste0("        Column ", ai, ": ", acc[[acc_key]], "\n"))
              }
            }
          }
        }
      }
      cat(file=stderr(), "==============================\n\n")
    }

    # DEBUG: Also print part of the YAML file
    # v53: cat(file=stderr(), "\n=== YAML FILE CONTENT (first 50 lines) ===\n")
    yaml_lines <- readLines(temp_yaml_file, n = 50)
    # v53: cat(file=stderr(), paste(yaml_lines, collapse = "\n"), "\n")
    # v53: cat(file=stderr(), "==========================================\n\n")
    
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
      # v53: cat(file=stderr(), "================================================\n\n")
      
      # Call func.print.lineage.tree with the temp YAML file
      # v54: Wrap in suppressWarnings to suppress -Inf and other harmless warnings
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
        legend_settings = values$legend_settings
      ))
      
      # Debug output
      # v53: cat(file=stderr(), "\n=== AFTER func.print.lineage.tree ===\n")
      # v53: cat(file=stderr(), "tree_result is NULL:", is.null(tree_result), "\n")
      if (!is.null(tree_result) && is.list(tree_result)) {
        # v53: cat(file=stderr(), "tree_result is a list with", length(tree_result), "element(s)\n")
        # v53: cat(file=stderr(), "List names:", paste(names(tree_result), collapse=", "), "\n")
        if (length(tree_result) > 0) {
          # v53: cat(file=stderr(), "First element class:", class(tree_result[[1]]), "\n")
          # v53: cat(file=stderr(), "First element inherits ggplot:", inherits(tree_result[[1]], "ggplot"), "\n")
          if (length(tree_result) > 1) {
            # v53: cat(file=stderr(), "NOTE: Multiple plots returned (", length(tree_result), "). Using the last one.\n")
          }
        }
      }
      # v53: cat(file=stderr(), "====================================\n\n")
      
      # Extract the plot from out_trees list
      # The function returns out_trees which is a list indexed by numbers like "1", "2", etc.
      # When multiple classifications exist, use the LAST plot (most complete)
      tree_plot <- NULL
      if (!is.null(tree_result) && is.list(tree_result) && length(tree_result) > 0) {
        # Extract the LAST plot from the list (most recent classification)
        plot_index <- length(tree_result)
        tree_plot <- tree_result[[plot_index]]
        # v53: cat(file=stderr(), "Successfully extracted plot from tree_result[[", plot_index, "]]\n")
        # v53: cat(file=stderr(), "Plot class:", class(tree_plot), "\n")
        # v53: cat(file=stderr(), "Plot inherits ggplot:", inherits(tree_plot, "ggplot"), "\n")
      } else {
        # v53: cat(file=stderr(), "ERROR: Could not extract plot from tree_result\n")
        # v53: cat(file=stderr(), "tree_result structure:\n")
        # v54: str(tree_result)
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
      # Print full error details to console for debugging
      # v53: cat(file=stderr(), "\n=== ERROR in generate_plot ===\n")
      # v53: cat(file=stderr(), "Error message:", e$message, "\n")
      # v53: cat(file=stderr(), "Full error:\n")
      # v53: print(e)
      # v53: cat(file=stderr(), "================================\n\n")
      
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
    # v53: cat(file=stderr(), "==============================\n\n")
    
    # If we got a valid result
    # If we got a valid result
    # If we got a valid result
    if (!is.null(result)) {
      # v53: cat(file=stderr(), "=== Attempting to save plot ===\n")

      # v121: Apply legend settings from the Legend tab
      legend_settings <- values$legend_settings
      if (!is.null(legend_settings)) {
        cat(file=stderr(), paste0("\n=== v121: Applying legend settings to plot ===\n"))
        cat(file=stderr(), paste0("  Position: ", legend_settings$position, "\n"))

        # v125: Determine legend layout based on position
        # For top/bottom: legends arranged horizontally, but each legend has title above values
        # For left/right: legends arranged vertically, with title next to values
        is_horizontal_position <- legend_settings$position %in% c("top", "bottom")

        # Build theme modifications for legend
        legend_theme <- theme(
          legend.position = legend_settings$position,
          legend.title = element_text(size = legend_settings$title_size, face = "bold"),
          legend.text = element_text(size = legend_settings$text_size),
          legend.key.size = unit(legend_settings$key_size, "lines"),
          legend.spacing = unit(legend_settings$spacing, "cm"),
          # v125: For top/bottom, arrange legends horizontally but stack items vertically
          legend.box = if (is_horizontal_position) "horizontal" else "vertical",
          legend.direction = if (is_horizontal_position) "vertical" else "vertical"
        )

        # Apply the legend theme
        result <- result + legend_theme

        cat(file=stderr(), paste0("  v125: Legend box=", if (is_horizontal_position) "horizontal" else "vertical",
                                   ", direction=vertical\n"))

        # Apply visibility controls using guides()
        guides_list <- list()

        # Hide specific legends based on visibility settings
        if (!isTRUE(legend_settings$show_classification)) {
          guides_list$colour <- "none"
        }
        if (!isTRUE(legend_settings$show_heatmap)) {
          guides_list$fill <- "none"
        }
        if (!isTRUE(legend_settings$show_bootstrap)) {
          guides_list$size <- "none"
        }

        if (length(guides_list) > 0) {
          result <- result + do.call(guides, guides_list)
        }

        cat(file=stderr(), paste0("  Legend settings applied successfully\n"))
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

        # v144: Store the offsets in result for later extraction
        # We'll apply them during the final rendering step
        attr(result, "plot_offset_x") <- plot_off_x
        attr(result, "plot_offset_y") <- plot_off_y
        attr(result, "plot_scale_percent") <- plot_scale

        if (plot_off_x != 0 || plot_off_y != 0 || plot_scale != 100) {
          cat(file=stderr(), paste0("\n=== v146: STORING PLOT POSITION & SCALE ===\n"))
          cat(file=stderr(), paste0("  X offset: ", plot_off_x, " (positive = right)\n"))
          cat(file=stderr(), paste0("  Y offset: ", plot_off_y, " (positive = up)\n"))
          cat(file=stderr(), paste0("  Scale: ", plot_scale, "%\n"))
          cat(file=stderr(), paste0("  v146: Offsets and scale will be applied via cowplot draw_plot during rendering\n"))
        }

        # Apply page title
        page_title_settings <- values$page_title

        # v131: DEBUG - trace page title settings
        cat(file=stderr(), paste0("\n=== v131: PAGE TITLE CHECK ===\n"))
        cat(file=stderr(), paste0("  page_title_settings is NULL: ", is.null(page_title_settings), "\n"))
        if (!is.null(page_title_settings)) {
          cat(file=stderr(), paste0("  enabled: ", page_title_settings$enabled, "\n"))
          cat(file=stderr(), paste0("  text: '", page_title_settings$text, "'\n"))
          cat(file=stderr(), paste0("  text length: ", nchar(page_title_settings$text), "\n"))
        }

        if (!is.null(page_title_settings) && isTRUE(page_title_settings$enabled) &&
            !is.null(page_title_settings$text) && nchar(page_title_settings$text) > 0) {
          cat(file=stderr(), paste0("\n=== v131: Applying page title ===\n"))
          cat(file=stderr(), paste0("  Title: ", page_title_settings$text, "\n"))
          cat(file=stderr(), paste0("  Size: ", page_title_settings$size, "\n"))
          cat(file=stderr(), paste0("  Color: ", page_title_settings$color, "\n"))

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
            cat(file=stderr(), paste0("  Adding underline\n"))
            # Note: Underline in ggplot title is complex, would need grid manipulation
            # For now, we document that underline is not fully supported
          }
          cat(file=stderr(), paste0("  Page title applied successfully\n"))
        } else {
          cat(file=stderr(), paste0("  Page title NOT applied (condition not met)\n"))
        }

        # v143: Apply custom text annotations as TRUE overlays using annotation_custom
        # annotation_custom with grid::textGrob is a true overlay that never affects plot limits
        # Unlike annotate() + coord_cartesian which can override ggtree's coordinate system
        custom_texts <- values$custom_texts
        if (!is.null(custom_texts) && length(custom_texts) > 0) {
          cat(file=stderr(), paste0("\n=== v143: Applying ", length(custom_texts), " custom text(s) as TRUE OVERLAY ===\n"))
          cat(file=stderr(), paste0("  Using annotation_custom with grid::textGrob (never affects plot limits)\n"))

          for (i in seq_along(custom_texts)) {
            txt <- custom_texts[[i]]

            cat(file=stderr(), paste0("  Text ", i, ": \"", substr(txt$content, 1, 20), "...\" at normalized (",
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

          cat(file=stderr(), paste0("  v143: Text overlays applied - plot coordinates unchanged\n"))
        }

        # v145: Store custom images for later application (after cowplot wrapping)
        # This ensures images are drawn on TOP of everything including cowplot canvas
        custom_images <- values$custom_images
        attr(result, "custom_images") <- custom_images
        if (!is.null(custom_images) && length(custom_images) > 0) {
          cat(file=stderr(), paste0("\n=== v145: ", length(custom_images), " custom image(s) queued for overlay ===\n"))
        }
      }, error = function(e) {
        cat(file=stderr(), paste0("  v130 Extra tab ERROR: ", e$message, "\n"))
      })

      # Store the plot with legend settings applied
      values$current_plot <- result

      # v148: Extract and output all legend coordinates with coordinate system explanations
      tryCatch({
        cat(file=stderr(), paste0("\n=== v148: LEGEND COORDINATE SYSTEMS EXPLAINED ===\n"))
        cat(file=stderr(), paste0("\n"))
        cat(file=stderr(), paste0("NOTE: There are TWO different coordinate systems:\n"))
        cat(file=stderr(), paste0("\n"))
        cat(file=stderr(), paste0("1) PLOT COORDINATES (used by Highlight & Bootstrap legends):\n"))
        cat(file=stderr(), paste0("   - These are DATA coordinates within the plot area\n"))
        cat(file=stderr(), paste0("   - x = tree depth direction (after coord_flip: vertical position)\n"))
        cat(file=stderr(), paste0("   - y = tip position direction (after coord_flip: horizontal position)\n"))
        cat(file=stderr(), paste0("   - Negative x values place items below the tree baseline\n"))
        cat(file=stderr(), paste0("   - Controlled via: highlight_x_offset, highlight_y_offset,\n"))
        cat(file=stderr(), paste0("                     bootstrap_x_offset, bootstrap_y_offset in Legend tab\n"))
        cat(file=stderr(), paste0("\n"))
        cat(file=stderr(), paste0("2) GRID COORDINATES (used by Classification/Heatmap/P-value legends):\n"))
        cat(file=stderr(), paste0("   - These are LAYOUT CELL positions in the rendered gtable\n"))
        cat(file=stderr(), paste0("   - Not directly comparable to plot coordinates\n"))
        cat(file=stderr(), paste0("   - Positioned by ggplot's legend.position theme setting\n"))
        cat(file=stderr(), paste0("\n"))
        cat(file=stderr(), paste0("To ALIGN legends: Use the Legend tab offset controls to move\n"))
        cat(file=stderr(), paste0("Highlight/Bootstrap legends up/down/left/right until visually aligned.\n"))
        cat(file=stderr(), paste0("\n"))

        # Build the plot to extract grob information
        plot_build <- ggplot2::ggplot_build(result)
        plot_gtable <- ggplot2::ggplot_gtable(plot_build)

        # Get plot data range to help understand coordinate scale
        if (!is.null(plot_build$layout$panel_params) && length(plot_build$layout$panel_params) > 0) {
          pp <- plot_build$layout$panel_params[[1]]
          cat(file=stderr(), paste0("PLOT DATA RANGE (for reference):\n"))
          if (!is.null(pp$x.range)) {
            cat(file=stderr(), paste0("  x-axis range: ", round(pp$x.range[1], 2), " to ", round(pp$x.range[2], 2), "\n"))
          }
          if (!is.null(pp$y.range)) {
            cat(file=stderr(), paste0("  y-axis range: ", round(pp$y.range[1], 2), " to ", round(pp$y.range[2], 2), "\n"))
          }
          cat(file=stderr(), paste0("\n"))
        }

        cat(file=stderr(), paste0("GGPLOT LEGENDS (Grid Coordinates):\n"))
        # Find all legend grobs
        legend_grobs <- which(grepl("guide-box", plot_gtable$layout$name))
        if (length(legend_grobs) > 0) {
          for (i in seq_along(legend_grobs)) {
            leg_idx <- legend_grobs[i]
            leg_name <- plot_gtable$layout$name[leg_idx]
            leg_l <- plot_gtable$layout$l[leg_idx]
            leg_r <- plot_gtable$layout$r[leg_idx]
            leg_t <- plot_gtable$layout$t[leg_idx]
            leg_b <- plot_gtable$layout$b[leg_idx]
            cat(file=stderr(), paste0("  Legend Box ", i, " ('", leg_name, "'):\n"))
            cat(file=stderr(), paste0("    Grid cell: column ", leg_l, "-", leg_r, ", row ", leg_t, "-", leg_b, "\n"))
          }
        } else {
          cat(file=stderr(), paste0("  No legend guide-boxes found in gtable\n"))
        }

        # Also extract legend titles from the built plot scales
        scales_info <- plot_build$plot$scales$scales
        if (!is.null(scales_info) && length(scales_info) > 0) {
          cat(file=stderr(), paste0("\n  Active Scales with Legends:\n"))
          for (i in seq_along(scales_info)) {
            scale_obj <- scales_info[[i]]
            if (!is.null(scale_obj$name) && nchar(as.character(scale_obj$name)) > 0) {
              cat(file=stderr(), paste0("    - '", scale_obj$name, "' (", paste(scale_obj$aesthetics, collapse=", "), ")\n"))
            }
          }
        }

        cat(file=stderr(), paste0("\n=================================================\n"))
      }, error = function(e) {
        cat(file=stderr(), paste0("  v148: Error extracting legend coords: ", e$message, "\n"))
      })

      # Create a unique temp file with timestamp to force browser refresh
      temp_plot_file <- file.path(tempdir(), paste0("shiny_plot_", Sys.getpid(), "_", 
                                                    format(Sys.time(), "%Y%m%d_%H%M%S_%OS3"), ".png"))
      # v53: cat(file=stderr(), "Temp file path:", temp_plot_file, "\n")
      
      # Clean up old plot files to avoid accumulation
      old_files <- list.files(tempdir(), pattern = paste0("^shiny_plot_", Sys.getpid()), full.names = TRUE)
      if (length(old_files) > 5) {  # Keep only last 5 files
        old_files_sorted <- old_files[order(file.info(old_files)$mtime)]
        files_to_delete <- old_files_sorted[1:(length(old_files_sorted) - 5)]
        sapply(files_to_delete, unlink)
      }
      
      tryCatch({
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

        cat(file=stderr(), paste0("  v145: Preview using DOWNLOAD TAB proportions\n"))
        cat(file=stderr(), paste0("  v145: User dimensions: ", user_width, " x ", user_height, " ", user_units, "\n"))
        cat(file=stderr(), paste0("  v145: Preview dimensions: ", round(preview_width, 2), " x ", round(preview_height, 2), " in\n"))

        # v146: Apply plot position offsets AND scale using cowplot for true transformation
        # This moves and scales the plot without squeezing/distorting proportions
        plot_to_save <- result
        offset_x <- attr(result, "plot_offset_x")
        offset_y <- attr(result, "plot_offset_y")
        scale_pct <- attr(result, "plot_scale_percent")
        if (is.null(scale_pct)) scale_pct <- 100

        # Check if we need to apply any transformation
        needs_transform <- (!is.null(offset_x) && !is.null(offset_y) && (offset_x != 0 || offset_y != 0)) || scale_pct != 100

        if (needs_transform) {
          cat(file=stderr(), paste0("\n=== v146: APPLYING PLOT POSITION AND SCALE WITH COWPLOT ===\n"))

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

          # Calculate position to center the scaled plot
          # When scale_factor = 1, x = 0 + x_offset, y = 0 + y_offset (top-left)
          # When scale_factor = 0.5, x = 0.25 + x_offset (centered at 0.5)
          # Formula: x = (1 - scale_factor) / 2 + x_offset
          center_x <- (1 - scale_factor) / 2 + x_pos_offset
          center_y <- (1 - scale_factor) / 2 + y_pos_offset

          cat(file=stderr(), paste0("  Scale: ", scale_pct, "% (factor: ", round(scale_factor, 3), ")\n"))
          cat(file=stderr(), paste0("  X position offset: ", round(x_pos_offset, 3), "\n"))
          cat(file=stderr(), paste0("  Y position offset: ", round(y_pos_offset, 3), "\n"))
          cat(file=stderr(), paste0("  Final position: (", round(center_x, 3), ", ", round(center_y, 3), ")\n"))

          # Use ggdraw to create a canvas and draw_plot to position and scale the plot
          # The plot is scaled uniformly (same width and height factor) to avoid distortion
          plot_to_save <- cowplot::ggdraw() +
            cowplot::draw_plot(result, x = center_x, y = center_y,
                              width = scale_factor, height = scale_factor)

          cat(file=stderr(), paste0("  v146: Plot wrapped in cowplot canvas with scale and offset positioning\n"))
        }

        # v145: Apply custom images as TRUE overlays using cowplot::draw_image
        # This is done AFTER cowplot wrapping to ensure images are on top
        custom_images <- attr(result, "custom_images")
        if (!is.null(custom_images) && length(custom_images) > 0) {
          cat(file=stderr(), paste0("\n=== v145: Applying ", length(custom_images), " custom image(s) as TRUE OVERLAY ===\n"))

          # Ensure we have a cowplot canvas
          if (!inherits(plot_to_save, "ggdraw")) {
            plot_to_save <- cowplot::ggdraw(plot_to_save)
          }

          for (i in seq_along(custom_images)) {
            img <- custom_images[[i]]
            if (file.exists(img$path)) {
              tryCatch({
                cat(file=stderr(), paste0("  Image ", i, ": ", img$name, "\n"))
                cat(file=stderr(), paste0("    Position: (", round(img$x, 2), ", ", round(img$y, 2), ")\n"))
                cat(file=stderr(), paste0("    Width: ", round(img$width, 3), "\n"))

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
                    cat(file=stderr(), paste0("    Auto height: ", round(img_height, 3), " (aspect: ", round(img_aspect, 2), ")\n"))
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
                cat(file=stderr(), paste0("    Image added as overlay using cowplot::draw_image\n"))
              }, error = function(e) {
                cat(file=stderr(), paste0("  ERROR loading image ", i, ": ", e$message, "\n"))
              })
            } else {
              cat(file=stderr(), paste0("  Image ", i, ": File not found - ", img$path, "\n"))
            }
          }
          cat(file=stderr(), paste0("  v145: Custom images applied as overlays\n"))
        }

        # v53: cat(file=stderr(), "Calling ggsave...\n")
        # v54: Wrap in suppressWarnings to suppress scale warnings
        suppressWarnings(ggsave(
          filename = temp_plot_file,
          plot = plot_to_save,
          width = preview_width,
          height = preview_height,
          units = "in",
          dpi = 150,
          limitsize = FALSE
        ))
        
        # v53: cat(file=stderr(), "ggsave completed\n")
        
        # Wait a moment for file system to sync
        Sys.sleep(0.5)
        
        # v53: cat(file=stderr(), "File exists after ggsave:", file.exists(temp_plot_file), "\n")
        
        if (file.exists(temp_plot_file)) {
          # v53: cat(file=stderr(), "File size:", file.info(temp_plot_file)$size, "bytes\n")

          # Store the file path - this will trigger renderImage
          values$temp_plot_file <- temp_plot_file
          values$plot_ready <- TRUE
          values$plot_counter <- values$plot_counter + 1  # Increment to force reactive update

          # Hide progress - plot is ready!
          values$progress_visible <- FALSE
          values$progress_message <- ""
          values$plot_generating <- FALSE  # Turn off generating indicator

          # v57: Show Ready status via shinyjs
          show_status_ready()

          # v53: cat(file=stderr(), "Plot saved successfully and path stored\n")
          # v53: cat(file=stderr(), "Plot counter now:", values$plot_counter, "\n")
        } else {
          # v53: cat(file=stderr(), "ERROR: File does not exist after ggsave\n")
          values$temp_plot_file <- NULL
          values$plot_ready <- FALSE
          values$progress_visible <- FALSE
          values$plot_generating <- FALSE
          # v57: Show click to generate on error
          show_status_click_to_generate()
        }
        
      }, error = function(e) {
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
      
      # v53: cat(file=stderr(), "==============================\n\n")
    }
    
    # Ensure plot_generating is always reset
    if (is.null(result)) {
      values$plot_generating <- FALSE
      values$progress_visible <- FALSE
    }
    
    # FAILSAFE: Always ensure plot_generating is reset
    values$plot_generating <- FALSE
    values$progress_visible <- FALSE
    
    # v53: cat(file=stderr(), "Finished generate_plot()\n")
  }  # End of generate_plot function
  
  
  
  # v57: tree_status_indicator replaced with static HTML + shinyjs for immediate updates
  # The status indicator is now controlled via show_status_waiting(), show_status_processing(),
  # show_status_ready(), and show_status_click_to_generate() helper functions
  
  # Output renderers (outside of generate_plot function)
  output$tree_preview <- renderImage({
    # v53: cat(file=stderr(), "\n=== renderImage called for tree_preview ===\n")
    
    # Force reactive update by depending on plot_counter
    req(values$temp_plot_file, values$plot_counter)
    
    # v53: cat(file=stderr(), "Temp file:", values$temp_plot_file, "\n")
    # v53: cat(file=stderr(), "File exists:", file.exists(values$temp_plot_file), "\n")
    # v53: cat(file=stderr(), "Plot counter:", values$plot_counter, "\n")
    
    list(
      src = values$temp_plot_file,
      contentType = "image/png",
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
      contentType = "image/png",
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
      contentType = "image/png",
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
      contentType = "image/png",
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
      contentType = "image/png",
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
      contentType = "image/png",
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
      contentType = "image/png",
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
      contentType = "image/png",
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
  observeEvent(input$plot_offset_x, {
    values$plot_offset_x <- input$plot_offset_x
  }, ignoreInit = TRUE)

  # v141: Observer for plot position Y slider
  observeEvent(input$plot_offset_y, {
    values$plot_offset_y <- input$plot_offset_y
  }, ignoreInit = TRUE)

  # v141: Reset plot position button
  observeEvent(input$reset_plot_position, {
    updateSliderInput(session, "plot_offset_x", value = 0)
    updateSliderInput(session, "plot_offset_y", value = 0)
    values$plot_offset_x <- 0
    values$plot_offset_y <- 0
  })

  # v146: Observer for plot scale slider
  observeEvent(input$plot_scale_percent, {
    values$plot_scale_percent <- input$plot_scale_percent
  }, ignoreInit = TRUE)

  # v146: Reset plot scale button
  observeEvent(input$reset_plot_scale, {
    updateSliderInput(session, "plot_scale_percent", value = 100)
    values$plot_scale_percent <- 100
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

  # Update YAML output
  output$yaml_output <- renderText({
    yaml::as.yaml(values$yaml_data, indent.mapping.sequence = TRUE)
  })

  ###################
  
  # Define YAML content reactive
  # Define YAML content reactive
  yaml_content <- reactive({
    # Check if we have the necessary data
    req(values$tree, values$csv_data)
    
    
    # Create the YAML structure based on current settings
    yaml_data <- list(
      "Individual general definitions" = list(
        Individual = input$individual_name,
        "tree path" = if (!is.null(input$tree_file)) {
          list(input$tree_file$datapath)
        } else {
          list(NULL)
        },  # Added missing comma here
        "mapping csv file" = if (!is.null(input$csv_file)) {
          input$csv_file$datapath
        } else {
          NULL
        },  # Added missing comma here
        "out_file" = list(
          "base_path" = input$output_path,
          "file_type" = input$output_format,
          "optional text at beggining" = input$prefix_text,
          "optional text at end" = input$suffix_text,
          "replace name" = list(
            flag = if (input$replace_name) {
              "yes"
            } else {
              "no"
            },
            name = input$custom_name
          )
        )
      ),
      "Mapping exl renaming titles" = list(
        "ID column" = input$id_column
      ),
      "visual definitions" = list(
        "classification" = list(),
        "Bootstrap" = list(
          display = if (input$show_bootstrap) {
            "yes"
          } else {
            "no"
          },
          format = input$bootstrap_format,
          param = as.character(input$bootstrap_param)
        ),
        "rotation1" = list(
          display = if (input$enable_rotation && 
                        (input$rotation_type == "primary" || input$rotation_type == "manual")) {
            "yes"
          } else {
            "no"
          },
          according = list()
        ),
        "rotation2" = list(
          display = if (input$enable_rotation && input$rotation_type == "secondary") {
            "yes"
          } else {
            "no"
          },
          according = list()
        ),
        "trim tips" = list(
          display = if (input$trim_tips) {
            "yes"
          } else {
            "no"
          },
          length = input$tip_length
        ),
        "edge_width_multiplier" = list(
          size = input$edge_width
        ),
        "font_size" = list(
          tips = input$tip_font_size,
          legend_title = 13,
          legend_text = 10,
          legend_box = 10,
          heat_map_title = 137,
          heat_map_legend = input$heatmap_font_size
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
        plot_to_download <- values$current_plot
        offset_x <- attr(values$current_plot, "plot_offset_x")
        offset_y <- attr(values$current_plot, "plot_offset_y")
        scale_pct <- attr(values$current_plot, "plot_scale_percent")
        if (is.null(scale_pct)) scale_pct <- 100

        needs_transform <- (!is.null(offset_x) && !is.null(offset_y) && (offset_x != 0 || offset_y != 0)) || scale_pct != 100

        if (needs_transform) {
          x_pos_offset <- if (!is.null(offset_x)) offset_x * 0.05 else 0
          y_pos_offset <- if (!is.null(offset_y)) offset_y * 0.05 else 0
          scale_factor <- scale_pct / 100

          # Calculate position to center the scaled plot
          center_x <- (1 - scale_factor) / 2 + x_pos_offset
          center_y <- (1 - scale_factor) / 2 + y_pos_offset

          plot_to_download <- cowplot::ggdraw() +
            cowplot::draw_plot(values$current_plot, x = center_x, y = center_y,
                              width = scale_factor, height = scale_factor)
          cat(file=stderr(), paste0("v146: Download plot positioned with cowplot (scale: ", scale_pct, "%)\n"))
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
          cat(file=stderr(), paste0("v134: Download saved successfully: ", file, "\n"))
        }, error = function(e) {
          cat(file=stderr(), paste0("v134: Error saving download: ", e$message, "\n"))
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
        cat(file=stderr(), "v134: No current_plot available, using basic tree plot\n")
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
  
  # Convert Shiny app settings to YAML format
  # Convert Shiny app settings to YAML format
  settings_to_yaml <- function(settings) {
    # Create basic structure with required fields
    yaml_structure <- list(
      "Individual general definitions" = list(
        Individual = input$individual_name,
        "tree path" = if (!is.null(input$tree_file)) {
          list(input$tree_file$datapath)
        } else {
          list(NA)
        },  # Added missing comma
        "mapping csv file" = if (!is.null(input$csv_file)) {
          input$csv_file$datapath
        } else {
          NA
        },  # Added missing comma
        "out_file" = list(
          "base_path" = input$output_path,
          "file_type" = input$output_format,
          "optional text at beggining" = input$prefix_text,
          "optional text at end" = input$suffix_text,
          "replace name" = list(
            flag = if (input$replace_name) {
              "yes"
            } else {
              "no"
            },
            name = input$custom_name
          )
        )
      ),
      "Mapping exl renaming titles" = list(
        "ID column" = input$id_column
      ),
      "visual definitions" = list(
        "font_size" = list(
          tips = 3,
          legend_title = 30,
          legend_text = 20,
          legend_box = 15,
          heat_map_title = 25,
          heat_map_legend = 3.8
        ),
        "compare_two_trees" = "no"
      )
    )
    
    # Add visual definitions
    yaml_structure$`visual definitions` <- list(
      "classification" = list(),
      "Bootstrap" = list(
        display = if (input$show_bootstrap) {
          "yes"
        } else {
          "no"
        },
        format = input$bootstrap_format,
        param = as.character(input$bootstrap_param)
      ),
      "rotation1" = list(
        display = if (input$enable_rotation && 
                      (input$rotation_type == "primary" || input$rotation_type == "manual")) {
          "yes"
        } else {
          "no"
        },
        according = list()
      ),
      "rotation2" = list(
        display = if (input$enable_rotation && input$rotation_type == "secondary") {
          "yes"
        } else {
          "no"
        },
        according = list()
      ),
      "trim tips" = list(
        display = if (input$trim_tips) {
          "yes"
        } else {
          "no"
        },
        length = input$tip_length
      ),
      "edge_width_multiplier" = list(
        size = input$edge_width
      ),
      "font_size" = list(
        tips = input$tip_font_size,
        legend_title = 30,
        legend_text = 20,
        legend_box = 15,
        heat_map_title = 25,
        heat_map_legend = input$heatmap_font_size
      )
    )
    
    # Create YAML text from the structure
    yaml_text <- yaml::as.yaml(yaml_structure, indent.mapping.sequence = TRUE)
    return(yaml_text)
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