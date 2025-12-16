# EZLineagePlotter Performance Optimization Guide

This document outlines performance issues identified in the EZLineagePlotter Shiny application and provides detailed recommendations for optimization.

## Current Performance Profile

- **File size**: 16,404 lines of R code
- **Observers**: 138 `observeEvent`/`observe` calls
- **Reactive expressions**: 31 `reactive()` calls
- **Render functions**: 61 `renderUI`/`renderImage`/`renderPlot` etc.
- **Typical plot generation time**: 3-8 seconds

---

## HIGH IMPACT Optimizations

### 1. Separate UI Updates from Plot Regeneration

**Problem**:
The app regenerates the entire plot (~3-8 seconds) whenever certain UI elements change. Users experience lag even when making small adjustments like changing colors, because they expect immediate visual feedback in the UI before applying to the plot.

**Current Behavior**:
- Some observers directly call `generate_plot()` or `request_plot_update()`
- UI responsiveness suffers during plot generation
- Users must wait for full regeneration to see each change

**Suggested Fix**:
1. Implement an "Apply Changes" button pattern for configuration sections
2. Allow users to make multiple changes, then apply all at once
3. Add a "Live Preview" toggle to let users choose between immediate updates and manual apply

**Implementation Details**:
```r
# Add apply button to heatmap configuration
actionButton("apply_heatmap_changes", "Apply Heatmap Settings", class = "btn-primary")

# Only regenerate plot when apply button is clicked
observeEvent(input$apply_heatmap_changes, {
  generate_plot()
})

# Remove generate_plot() calls from individual heatmap observers
```

**What Can Go Wrong**:
- Users may forget to click "Apply" and think their changes didn't save
- Need to provide visual indication that changes are pending
- Must ensure configuration is saved even without applying to plot

**Testing Checklist**:
- [ ] UI remains responsive while making changes
- [ ] Changes are correctly stored in `values$heatmap_configs`
- [ ] Apply button triggers correct plot regeneration
- [ ] Visual indicator shows when changes are pending
- [ ] Settings persist correctly on YAML export/import

---

### 2. Lazy Loading for Tabs

**Problem**:
All tab contents are rendered on app startup, even tabs the user may never visit. This increases initial load time and memory usage.

**Current Behavior**:
- All `renderUI` outputs are computed immediately
- Tab contents are fully built even if hidden
- Memory is allocated for all UI elements upfront

**Suggested Fix**:
1. Use conditional rendering based on which tab is active
2. Only compute expensive UI elements when the tab is first accessed
3. Cache computed UI after first access

**Implementation Details**:
```r
# Track which tabs have been visited
values$visited_tabs <- list()

# Lazy render for rotation tab
output$rotation_status_box <- renderUI({
  # Only render if tab has been visited
  req(input$main_tabs == "rotation_tab" || "rotation" %in% values$visited_tabs)
  values$visited_tabs <- c(values$visited_tabs, "rotation")

  # ... existing render code ...
})
```

**What Can Go Wrong**:
- First-time tab access may feel slow
- Need to handle cases where data isn't ready when tab is accessed
- State synchronization issues if tab content depends on other tabs

**Testing Checklist**:
- [ ] Initial app load is faster
- [ ] Each tab renders correctly on first access
- [ ] Tab content updates correctly when underlying data changes
- [ ] No errors when switching tabs rapidly
- [ ] Memory usage is reduced

---

### 3. Reduce Observer Count

**Problem**:
138 observers create significant overhead. Each observer maintains state and is checked on every reactive flush, even if its inputs haven't changed.

**Current Behavior**:
- Many small observers for related functionality (e.g., separate observer for each color input)
- Observers inside loops create dynamic observer proliferation
- Some observers may fire unnecessarily due to broad reactive dependencies

**Suggested Fix**:
1. Consolidate related observers using `eventReactive()` with multiple triggers
2. Use `bindEvent()` pattern for cleaner reactive chains
3. Move loop-based observer creation to a single parameterized observer

**Implementation Details**:
```r
# BEFORE: Multiple separate observers
observeEvent(input$heatmap_low_color_1, { ... })
observeEvent(input$heatmap_mid_color_1, { ... })
observeEvent(input$heatmap_high_color_1, { ... })

# AFTER: Single consolidated observer
observeEvent({
  list(
    input$heatmap_low_color_1,
    input$heatmap_mid_color_1,
    input$heatmap_high_color_1
  )
}, {
  if (1 <= length(values$heatmap_configs)) {
    values$heatmap_configs[[1]]$low_color <- input$heatmap_low_color_1
    values$heatmap_configs[[1]]$mid_color <- input$heatmap_mid_color_1
    values$heatmap_configs[[1]]$high_color <- input$heatmap_high_color_1
  }
}, ignoreInit = TRUE)
```

**What Can Go Wrong**:
- Consolidated observers may fire when only one input changes (unnecessary updates)
- Need to carefully handle `ignoreInit` and `ignoreNULL` parameters
- May break existing reactive chains if not done carefully

**Testing Checklist**:
- [ ] All functionality still works after consolidation
- [ ] No duplicate processing when single input changes
- [ ] `ignoreInit = TRUE` prevents startup cascade
- [ ] YAML export/import still works correctly
- [ ] Performance profiler shows reduced observer overhead

---

## MEDIUM IMPACT Optimizations

### 4. Add Debouncing to Color Pickers

**Problem**:
`colourpicker` inputs fire events continuously while the user drags through the color palette. This causes many rapid reactive updates.

**Current Behavior**:
- Each tiny color change triggers observer
- Observers update `values$heatmap_configs` many times per second
- Even though plot doesn't regenerate, the reactive system is stressed

**Suggested Fix**:
Add debouncing to color-related reactive values.

**Implementation Details**:
```r
# Define debounce delay
COLOR_DEBOUNCE_MS <- 300

# Create debounced reactive for each color input
heatmap_low_color_1_d <- debounce(
  reactive(input$heatmap_low_color_1),
  COLOR_DEBOUNCE_MS
)

# Use debounced version in observer
observeEvent(heatmap_low_color_1_d(), {
  if (1 <= length(values$heatmap_configs)) {
    values$heatmap_configs[[1]]$low_color <- heatmap_low_color_1_d()
  }
}, ignoreInit = TRUE)
```

**What Can Go Wrong**:
- 300ms delay may feel sluggish for some users
- Need to balance responsiveness vs. performance
- Debounced values may cause timing issues with dependent reactives

**Testing Checklist**:
- [ ] Color picker still feels responsive
- [ ] Final color value is correctly captured
- [ ] No "lost" color changes when user changes quickly then clicks away
- [ ] Works correctly with Apply button pattern (if implemented)

---

### 5. Cache Intermediate Results

**Problem**:
Expensive computations are repeated unnecessarily. Tree parsing, data transformations, and weight calculations happen on every plot regeneration.

**Current Behavior**:
- Tree is re-processed from YAML each time
- Heatmap data is transformed repeatedly
- Weight calculations for rotation run even when rotation settings haven't changed

**Suggested Fix**:
1. Use `memoise` package for pure function caching
2. Store computed results in `values$` and only recompute when inputs change
3. Add checksum-based cache invalidation

**Implementation Details**:
```r
library(memoise)

# Memoize expensive pure functions
func.create.weight_list_cached <- memoise(func.create.weight_list)

# Store computed results with cache keys
values$tree_cache <- list(
  tree = NULL,
  tree_path_hash = NULL
)

# Only recompute if tree file changed
current_hash <- digest::digest(tree_file$datapath)
if (is.null(values$tree_cache$tree_path_hash) ||
    values$tree_cache$tree_path_hash != current_hash) {
  values$tree_cache$tree <- read.tree(tree_file$datapath)
  values$tree_cache$tree_path_hash <- current_hash
}
```

**What Can Go Wrong**:
- Cache invalidation bugs (stale data shown)
- Memory usage increase from cached data
- Need to correctly identify cache keys (what triggers recomputation)

**Testing Checklist**:
- [ ] Cached results are correct
- [ ] Cache invalidates when source data changes
- [ ] Memory usage is acceptable
- [ ] No stale data displayed after file re-upload
- [ ] Performance improvement is measurable

---

### 6. Reduce renderUI Calls

**Problem**:
61+ `renderUI` calls rebuild UI elements dynamically. Each rebuild causes browser DOM manipulation and potential layout reflow.

**Current Behavior**:
- Classification, highlight, and rotation UIs are rebuilt on changes
- Select inputs are rebuilt instead of updated
- Some renderUI calls happen more frequently than necessary

**Suggested Fix**:
1. Use `updateSelectInput()` family functions instead of rebuilding
2. Pre-render UI templates and show/hide with `shinyjs`
3. Use `insertUI`/`removeUI` for truly dynamic content

**Implementation Details**:
```r
# BEFORE: Rebuilding select input
output$classification_values_ui <- renderUI({
  req(input$classification_column)
  values <- unique(values$csv_data[[input$classification_column]])
  selectInput("classification_value", "Value", choices = values)
})

# AFTER: Updating existing select input
observeEvent(input$classification_column, {
  req(input$classification_column)
  values <- unique(values$csv_data[[input$classification_column]])
  updateSelectInput(session, "classification_value", choices = values)
})
```

**What Can Go Wrong**:
- `updateSelectInput` may not handle all edge cases (e.g., selected value not in new choices)
- Show/hide approach requires pre-rendering more UI upfront
- May break existing CSS styling that depends on DOM structure

**Testing Checklist**:
- [ ] All select inputs update correctly
- [ ] Selected values persist appropriately
- [ ] No visual glitches during updates
- [ ] Accessibility is maintained
- [ ] Works correctly with dynamically created inputs (loop-based)

---

### 7. Plot Incremental Updates

**Problem**:
The entire ggplot is regenerated for every change, even when only a small visual element changed (e.g., legend position).

**Current Behavior**:
- `func.print.lineage.tree` rebuilds everything from scratch
- Plot generation takes 3-8 seconds regardless of change size
- No caching of unchanged plot layers

**Suggested Fix**:
1. Cache the base tree plot and layer additional elements
2. Consider using `plotly` with `plotlyProxy()` for incremental updates
3. Separate "structural" changes (requiring full rebuild) from "visual" changes (color, legend, etc.)

**Implementation Details**:
```r
# Cache base tree structure
if (is.null(values$base_tree_plot) || values$tree_structure_changed) {
  values$base_tree_plot <- create_base_tree_plot(tree_data)
  values$tree_structure_changed <- FALSE
}

# Add layers incrementally
final_plot <- values$base_tree_plot +
  add_heatmap_layer(heatmap_config) +
  add_legend_styling(legend_config)
```

**What Can Go Wrong**:
- ggplot objects can be large in memory
- Some changes may require full rebuild (e.g., tree rotation)
- Layer ordering dependencies may cause issues

**Testing Checklist**:
- [ ] Cached base plot produces identical results
- [ ] Incremental updates are visually correct
- [ ] Memory usage is acceptable
- [ ] Full rebuild still works when needed
- [ ] Performance improvement is significant (>50% for visual-only changes)

---

## LOWER IMPACT Optimizations (Quick Wins)

### 8. Reduce Debug Output

**Problem**:
Many `cat(file=stderr(), ...)` calls throughout the code add I/O overhead, especially in hot paths.

**Current Behavior**:
- Debug output on every plot generation
- Multiple debug statements in loops
- No conditional check before outputting

**Suggested Fix**:
Wrap debug output in conditional checks or use a debug flag.

**Implementation Details**:
```r
# Add debug flag at top of file
DEBUG_MODE <- FALSE

# Wrap debug output
if (DEBUG_MODE) {
  cat(file=stderr(), paste0("[DEBUG] Some message\n"))
}

# Or create helper function
debug_log <- function(msg) {
  if (getOption("ezlineage.debug", FALSE)) {
    cat(file=stderr(), paste0("[DEBUG] ", msg, "\n"))
  }
}
```

**What Can Go Wrong**:
- May lose useful debugging information in production
- Need to ensure flag is easily togglable for development
- Some debug output may be needed for user-facing error messages

**Testing Checklist**:
- [ ] Debug output is suppressed when flag is FALSE
- [ ] Important error messages still appear
- [ ] Debug mode can be easily enabled for troubleshooting
- [ ] Performance improvement is measurable (especially in loops)

---

### 9. Optimize Large Data Handling

**Problem**:
The debug output shows: `[PERF] Removed 16221 empty/auto-named columns on CSV load (keeping 99)`. This suggests heavy data cleaning on every CSV load, and potentially inefficient data structures.

**Current Behavior**:
- CSV loaded as data.frame
- Empty/auto-named columns filtered on load
- Data may be copied multiple times during processing

**Suggested Fix**:
1. Use `data.table` for large dataset operations (10-100x faster)
2. Do column filtering once on upload, not repeatedly
3. Pre-filter columns before storing in reactive values
4. Use column indices instead of names where possible

**Implementation Details**:
```r
library(data.table)

# On CSV upload - clean once and store
observeEvent(input$csv_file, {
  req(input$csv_file)

  # Read with data.table (faster)
  dt <- fread(input$csv_file$datapath)

  # Remove empty/unnamed columns ONCE
  # Identify columns to keep
  valid_cols <- names(dt)[!grepl("^V[0-9]+$|^\\.\\.\\.[0-9]+$", names(dt))]
  valid_cols <- valid_cols[sapply(dt[, ..valid_cols], function(x) !all(is.na(x) | x == ""))]

  # Store only valid columns
  values$csv_data <- as.data.frame(dt[, ..valid_cols])
  values$csv_columns <- valid_cols

  cat(file=stderr(), sprintf("[PERF] CSV loaded: %d rows, %d columns (removed %d invalid columns)\n",
      nrow(values$csv_data), length(valid_cols), ncol(dt) - length(valid_cols)))
})
```

**What Can Go Wrong**:
- `data.table` syntax differs from `data.frame` (may break existing code)
- Need to convert back to data.frame for compatibility with some functions
- Memory spike during initial load if dataset is very large

**Testing Checklist**:
- [ ] CSV loads correctly with various file formats
- [ ] Column filtering removes correct columns
- [ ] Data integrity is preserved
- [ ] Memory usage is reduced
- [ ] Downstream operations still work (heatmap, classification, etc.)
- [ ] Performance improvement is measurable

**Files to Modify**:
- CSV upload observer (around line 9200-9300)
- Any code that accesses `values$csv_data` with column filtering

---

### 10. Use `freezeReactiveValue()`

**Problem**:
When updating multiple related inputs programmatically, each update triggers reactive cascade. This causes unnecessary intermediate computations.

**Current Behavior**:
- `updateSelectInput()` and `updateColourInput()` calls trigger observers immediately
- Rapid sequential updates cause reactive "storm"
- Intermediate states may cause errors or warnings

**Suggested Fix**:
Use `freezeReactiveValue()` before batch updates.

**Implementation Details**:
```r
# When switching data source to RData, update multiple colors
observeEvent(input[[paste0("heatmap_data_source_", i)]], {
  if (input[[paste0("heatmap_data_source_", i)]] == "rdata") {
    # Freeze inputs before updating
    freezeReactiveValue(input, paste0("heatmap_low_color_", i))
    freezeReactiveValue(input, paste0("heatmap_mid_color_", i))
    freezeReactiveValue(input, paste0("heatmap_high_color_", i))

    # Now update all at once
    updateColourInput(session, paste0("heatmap_low_color_", i), value = "#FF0000")
    updateColourInput(session, paste0("heatmap_mid_color_", i), value = "#FFFFFF")
    updateColourInput(session, paste0("heatmap_high_color_", i), value = "#0000FF")
  }
}, ignoreInit = TRUE)
```

**What Can Go Wrong**:
- Frozen values won't trigger dependent observers until unfrozen
- Need to understand which observers depend on frozen inputs
- May cause unexpected behavior if freeze/unfreeze timing is wrong

**Testing Checklist**:
- [ ] Batch updates complete without intermediate observer firing
- [ ] Final values are correct after all updates
- [ ] Dependent UI elements update correctly
- [ ] No "stuck" frozen values

---

### 11. Profile with `profvis`

**Problem**:
Without profiling data, optimization efforts may target the wrong areas.

**Current Behavior**:
- Manual timing with `[PERF]` log statements
- No systematic profiling

**Suggested Fix**:
Add `profvis` profiling to identify actual bottlenecks.

**Implementation Details**:
```r
# Run app with profiling
library(profvis)

profvis({
  # Simulate typical user workflow
  runApp('EZlineagePlotter56.R')
})

# Or profile specific function
profvis({
  result <- func.print.lineage.tree(yaml_path, width=170, height=60)
})
```

**What Can Go Wrong**:
- Profiling adds overhead (results may not reflect real-world performance)
- Interactive Shiny apps are harder to profile than batch scripts
- May generate large profile files

**Testing Checklist**:
- [ ] Profile successfully captures app execution
- [ ] Bottlenecks are clearly identified
- [ ] Optimization targets are validated with profile data

---

## Additional Suggestions

### 12. Use `shiny.maxRequestSize` Wisely

**Problem**:
Large file uploads can crash the app or cause timeouts.

**Suggested Fix**:
```r
options(shiny.maxRequestSize = 50*1024^2)  # 50 MB limit
```

### 13. Implement Progress Indicators

**Problem**:
Users don't know if the app is working or frozen during long operations.

**Suggested Fix**:
Use `withProgress()` or `shiny::Progress` for long-running operations.

### 14. Consider Server-Side Processing for Large Tables

**Problem**:
Large DataTables rendered client-side cause browser lag.

**Suggested Fix**:
```r
output$data_table <- DT::renderDataTable({
  datatable(values$csv_data,
            options = list(
              processing = TRUE,
              serverSide = TRUE,  # Enable server-side processing
              pageLength = 25
            ))
})
```

### 15. Optimize Image Output

**Problem**:
High-resolution plot images are large and slow to transfer.

**Suggested Fix**:
- Use appropriate resolution for preview (72-96 DPI)
- Higher resolution only for export
- Consider SVG format for scalable graphics

---

## Implementation Priority

| Priority | Item | Effort | Impact |
|----------|------|--------|--------|
| 1 | #9 Optimize Large Data Handling | Medium | High |
| 2 | #4 Debounce Color Pickers | Low | Medium |
| 3 | #8 Reduce Debug Output | Low | Low-Medium |
| 4 | #10 Use freezeReactiveValue | Low | Medium |
| 5 | #1 Apply Button Pattern | Medium | High |
| 6 | #5 Cache Intermediate Results | High | High |
| 7 | #3 Reduce Observer Count | High | Medium |
| 8 | #6 Reduce renderUI Calls | Medium | Medium |
| 9 | #2 Lazy Loading for Tabs | Medium | Medium |
| 10 | #7 Plot Incremental Updates | Very High | Very High |

---

## Measuring Success

Before and after each optimization:

1. **Startup time**: Time from `runApp()` to interactive UI
2. **Plot generation time**: Time from clicking "Generate" to plot display
3. **UI responsiveness**: Time to update color picker, dropdown, etc.
4. **Memory usage**: Peak memory during typical workflow
5. **Observer count**: `length(reactlog::reactlog()$nodes)`

Use consistent test data and workflow for meaningful comparisons.
