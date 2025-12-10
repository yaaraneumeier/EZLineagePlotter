# Performance Optimization Plan for EZLineagePlotter

## Overview

This document outlines potential performance improvements for EZLineagePlotter, including risks and mitigation strategies for each approach.

---

## Option 1: Reactive Debouncing

### Description
Add `debounce()` to slider and numeric inputs to delay plot regeneration until the user stops adjusting values.

### Current Problem
Every tiny movement of a slider (e.g., legend offset, ellipse size) triggers a full plot regeneration, causing lag.

### Implementation
```r
# Instead of directly using input$legend_x_offset
legend_x_offset_debounced <- debounce(reactive(input$legend_x_offset), 300)

# Then use legend_x_offset_debounced() in observers
```

### Potential Problems

| Problem | Likelihood | Mitigation |
|---------|------------|------------|
| UI feels unresponsive if debounce too long | Medium | Use 200-400ms delay - fast enough to feel responsive |
| Some inputs need immediate response | Low | Only debounce sliders/numeric inputs, not buttons or checkboxes |
| Debounced value might be stale on first render | Low | Use `ignoreInit = TRUE` or provide sensible defaults |

### Risk Level: **Low**

### Estimated Impact: **High** - Significantly reduces unnecessary renders

---

## Option 2: Plot Caching with bindCache()

### Description
Cache expensive plot computations so identical inputs return cached results instantly.

### Current Problem
Even when inputs haven't changed, the plot may be rebuilt unnecessarily.

### Implementation
```r
output$tree_plot <- renderPlot({
  # plot generation code
}) %>% bindCache(
  input$tree_file,
  input$csv_file,
  values$classifications,
  values$heatmaps,
  # ... other relevant inputs
)
```

### Potential Problems

| Problem | Likelihood | Mitigation |
|---------|------------|------------|
| Cache key doesn't capture all dependencies | High | Carefully list ALL inputs that affect the plot; test thoroughly |
| Stale cache shows outdated plot | Medium | Include reactive values in cache key, not just inputs |
| Memory usage increases with cache size | Medium | Set `cache = cachem::cache_mem(max_size = 100 * 1024^2)` to limit |
| Complex nested reactive values hard to cache | High | May need to serialize/hash complex objects for cache key |
| Cache invalidation issues | Medium | Clear cache on major state changes (new file upload) |

### Risk Level: **Medium-High**

### Estimated Impact: **High** - But requires careful implementation

---

## Option 3: Data Processing Optimization

### Description
Optimize the "Process Data" step using faster data structures and vectorized operations.

### Current Problem
CSV matching and validation can be slow for large files due to row-by-row operations.

### Implementation
```r
# Use data.table for faster operations
library(data.table)
csv_dt <- as.data.table(csv_data)

# Vectorized matching instead of loops
matched_indices <- match(tree_tips, csv_dt[[id_column]])
```

### Potential Problems

| Problem | Likelihood | Mitigation |
|---------|------------|------------|
| data.table syntax differs from data.frame | High | Wrap in helper functions to isolate syntax differences |
| Existing code assumes data.frame behavior | High | Test all downstream code; some ggplot/ggtree functions may need data.frame |
| New dependency added | Low | data.table is stable and widely used |
| Memory copying during conversion | Low | Convert once at import, keep as data.table throughout |
| Factor handling differs | Medium | Explicitly convert factors where needed |

### Risk Level: **Medium**

### Estimated Impact: **Medium** - Most noticeable with large datasets (1000+ tips)

---

## Option 4: Reduce Redundant Recalculations

### Description
Use `isolate()`, `req()`, and `reactiveVal` to prevent unnecessary reactive cascade updates.

### Current Problem
Changing one input can trigger multiple observers in sequence, each regenerating the plot.

### Implementation
```r
# Store intermediate results
base_tree_plot <- reactiveVal(NULL)

# Only update when specific inputs change
observeEvent(input$tree_file, {
  # Build base tree
  base_tree_plot(new_plot)
}, ignoreInit = TRUE)

# Use req() to short-circuit
output$tree_plot <- renderPlot({
  req(values$tree)  # Don't run if tree is NULL
  req(values$csv_data)  # Don't run if CSV is NULL
  # ... rest of plot code
})

# Use isolate() to read without creating dependency
observeEvent(input$apply_button, {
  current_tree <- isolate(values$tree)  # Read but don't depend
  # ... process
})
```

### Potential Problems

| Problem | Likelihood | Mitigation |
|---------|------------|------------|
| Over-isolation breaks reactivity | High | Only isolate values that shouldn't trigger updates |
| req() blocks render permanently if condition never met | Medium | Ensure conditions can be satisfied; add error messages |
| reactiveVal not updated when expected | Medium | Track all places that should update it |
| Debugging becomes harder | Medium | Add debug_cat() statements to track reactive flow |
| Race conditions between observers | Low | Use `priority` parameter or consolidate observers |

### Risk Level: **Medium**

### Estimated Impact: **High** - Can dramatically reduce redundant work

---

## Option 5: Layer Management Optimization

### Description
Reduce calls to layer manipulation functions and optimize their implementation.

### Current Problem
`func.move.tiplabels.to.front()` is called multiple times per render cycle, each time iterating through all layers.

### Implementation
```r
# Call only once at the very end, right before rendering
# Remove intermediate calls

# Optimize the function itself
func.move.tiplabels.to.front <- function(p, verbose = FALSE) {
  if (is.null(p$layers) || length(p$layers) == 0) return(p)

  # Use which() instead of multiple conditionals
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  text_idx <- which(layer_classes %in% c("GeomText", "GeomLabel", "GeomTiplab"))

  if (length(text_idx) == 0) return(p)

  # Reorder in one operation
  other_idx <- setdiff(seq_along(p$layers), text_idx)
  p$layers <- p$layers[c(other_idx, text_idx)]

  return(p)
}
```

### Potential Problems

| Problem | Likelihood | Mitigation |
|---------|------------|------------|
| Removing intermediate calls breaks layer order | Medium | Test all code paths that add layers |
| Some layers added after "final" call | Medium | Ensure the final call is truly final |
| Layer class names may vary by ggplot2 version | Low | Test with multiple ggplot2 versions |
| vapply may fail if unexpected geom structure | Low | Add tryCatch wrapper |

### Risk Level: **Low-Medium**

### Estimated Impact: **Medium** - More noticeable with many layers (heatmaps + highlights)

---

## Option 6: Progress Indicators

### Description
Add progress bars for long operations to improve perceived performance.

### Current Problem
Users don't know if the app is working or frozen during long operations.

### Implementation
```r
# Use withProgress for long operations
withProgress(message = 'Generating plot...', value = 0, {

  incProgress(0.2, detail = "Building tree")
  # tree code

  incProgress(0.3, detail = "Adding heatmaps")
  # heatmap code

  incProgress(0.3, detail = "Adding legends")
  # legend code

  incProgress(0.2, detail = "Finalizing")
  # final code
})
```

### Potential Problems

| Problem | Likelihood | Mitigation |
|---------|------------|------------|
| Progress percentages inaccurate | Low | Use estimates; users care more about activity than accuracy |
| Nested withProgress calls conflict | Medium | Use single withProgress at top level |
| Progress bar doesn't appear for fast operations | Low | Only add for operations known to be slow |
| Slight overhead from progress updates | Very Low | Negligible compared to actual work |

### Risk Level: **Very Low**

### Estimated Impact: **Low** (actual) / **High** (perceived) - Doesn't speed up, but feels faster

---

## Recommended Implementation Order

Based on risk vs. impact:

1. **Option 6: Progress Indicators** - Zero risk, immediate UX improvement
2. **Option 1: Reactive Debouncing** - Low risk, high impact
3. **Option 5: Layer Management** - Low-medium risk, quick win
4. **Option 4: Reduce Redundant Recalculations** - Medium risk, high impact (do incrementally)
5. **Option 3: Data Processing** - Medium risk, good for large datasets
6. **Option 2: Plot Caching** - High risk, save for last (most complex)

---

## Testing Strategy

For each optimization:

1. **Before implementing**: Document current behavior and timing
2. **After implementing**:
   - Test all features still work
   - Measure timing improvement
   - Test edge cases (empty data, large data, rapid input changes)
3. **Regression tests**:
   - Upload tree + CSV
   - Add classification
   - Add 3 heatmaps
   - Add highlight
   - Add bootstrap
   - Change legend settings
   - Export PDF/PNG
