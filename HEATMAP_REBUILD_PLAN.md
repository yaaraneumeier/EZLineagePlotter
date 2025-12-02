# Heatmap Rebuild Plan

## Overview

This document contains the step-by-step plan for rebuilding the heatmap functionality in EZLineagePlotter.
The complex heatmap code (~800 lines) has been disabled in v91 and replaced with a simple proof-of-concept.

**Current Status**: v108 - Fixed heatmap UI issues: color settings now show immediately when columns selected, row height slider works properly, row labels appear on right side, added per-column custom label mapping

## Why We Simplified

The previous heatmap code had accumulated many fixes (v58-v90) trying to solve rendering issues:
- "Problem while setting up geom" errors
- Heatmap disappearing entirely
- Mapping corruption issues
- Complex duplicate gheatmap calls
- Multiple scale application conflicts

Rather than continuing to patch, we're rebuilding from scratch with a tested baseline.

---

## Phase 1: Validate Simple Heatmap Works (Current - v91)

**Goal**: Prove that a basic gheatmap call renders correctly

**Implementation** (done in v91):
- Single gheatmap call on tree object
- Use viridis_d for discrete data OR gradient2 for continuous data
- Simple hexpand for x-axis expansion
- Full diagnostic logging

**Test Cases**:
- [ ] Load tree + CSV with discrete heatmap column
- [ ] Load tree + CSV with continuous heatmap column
- [ ] Verify heatmap tiles appear next to tree
- [ ] Verify legend appears

**If v91 heatmap DOES NOT appear**: The issue is deeper (in data preparation, rowname matching, or ggtree version incompatibility)

**If v91 heatmap DOES appear**: Proceed to Phase 2

---

## Phase 2: Add User-Configurable Width and Offset

**Goal**: Allow users to adjust heatmap position

**Add**:
- Use `man_adj_heat_loc` for offset adjustment
- Use `width_heatmap` from parameters for width
- Read from `heat_display_params_list` for per-heatmap settings

**Code Location**: Lines 4567-4573 in v91

**Implementation**:
```r
# Replace hardcoded values:
new_heat_x <- 0.7  # -> use: off_base + man_adj_heat_loc
wi <- 1.4 * ncol(heat_data) / 23  # -> use: width_heatmap[1] if provided
```

---

## Phase 3: Add Custom Color Scales for Discrete Data

**Goal**: Support RColorBrewer palettes and custom color vectors

**Add**:
- Check `heat_param['man']` for manual palette mode
- Check `heat_param['man_define_colors']` for custom colors
- Check `heat_param[['color_scale_option']]` for palette name

**Implementation**:
```r
if (is_discrete) {
  if (heat_param['man_define_colors'] == TRUE) {
    # Use custom colors from heat_param[['color_scale_option']]
    p <- p + scale_fill_manual(values = custom_colors, ...)
  } else {
    palette_name <- heat_param[['color_scale_option']]
    if (palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      p <- p + scale_fill_brewer(palette = palette_name, ...)
    } else {
      p <- p + scale_fill_viridis_d(...)
    }
  }
}
```

---

## Phase 4: Add Custom Color Scales for Continuous Data

**Goal**: Support user-defined low/mid/high colors and limits

**Add**:
- Read `heat_param['low']`, `heat_param['mid']`, `heat_param['high']`
- Read `heat_param[['limits']]` for scale limits
- Support midpoint configuration

**Implementation** (already partially in v91, enhance):
```r
if (!is_discrete) {
  limits <- heat_param[['limits']]
  if (!is.na(limits[1])) {
    p <- p + scale_fill_gradient2(..., limits = limits)
  } else {
    p <- p + scale_fill_gradient2(...)
  }
}
```

---

## Phase 5: Support Multiple Heatmaps

**Goal**: Allow 2-3 heatmaps side by side

**Key Concept**: Each additional heatmap needs:
1. `new_scale_fill()` to reset the fill scale
2. Calculated offset based on previous heatmap position
3. Its own color scale

**Implementation**:
```r
for (j1 in 1:length(heat_map_title_list)) {
  if (heat_display_vec[j1] == FALSE) next

  if (j > 1) {
    p <- p + new_scale_fill()
    new_heat_x <- previous_heat_end + margin
  }

  p <- gheatmap(p, data = dxdf440_for_heat[[j1]], offset = new_heat_x, ...)
  p <- p + appropriate_color_scale(heat_display_params_list[[j1]])

  previous_heat_end <- calculate_heat_end(p)
}
```

**Critical**: Track heatmap counter `j` carefully (was source of many bugs in original code)

---

## Phase 6: Add Column Labels and Formatting

**Goal**: Support column name display with angles and custom labels

**Add**:
- `flag_colnames[j1]` - whether to show column names
- `heat_maps_titles_angles_vector[j1]` - rotation angle
- `custom_column_labels` - override column names
- `heat_legend_replace` - replace specific legend entries

**Implementation**:
```r
if (flag_colnames[j1] == TRUE) {
  custom_labels <- colnames(heat_data)
  # Apply heat_legend_replace if defined
} else {
  custom_labels <- rep("", ncol(heat_data))
}

p <- gheatmap(...,
  colnames = TRUE,
  colnames_angle = heat_maps_titles_angles_vector[j1],
  custom_column_labels = custom_labels
)
```

---

## Phase 7: X-Axis Expansion and Viewport

**Goal**: Ensure heatmap is visible within plot bounds

**Add**:
- Calculate expected x-range based on tree width + heatmap width
- Use `hexpand()` or `xlim_expand()` appropriately
- Handle edge cases where expansion fails

**Implementation**:
```r
# Calculate how much expansion needed
tree_xmax <- max(p$data$x, na.rm = TRUE)
expected_xmax <- tree_xmax + offset + width + margin

# Try expansion methods in order of preference
tryCatch({
  p <- p + ggtree::hexpand(ratio = expansion_ratio, direction = 1)
}, error = function(e) {
  tryCatch({
    p <- p + ggtree::xlim_expand(c(0, expected_xmax), "right")
  }, error = function(e2) {
    # Use geom_blank fallback
    p <- p + geom_blank(data = data.frame(x = expected_xmax, y = 1), aes(x=x, y=y))
  })
})
```

---

## Phase 8: Integration with Bootstrap Triangles

**Goal**: Ensure bootstrap display works with heatmap

**Key Issue**: Bootstrap triangle layers added BEFORE gheatmap lose x/y aesthetics

**Solution** (implemented in v90):
- Store bootstrap parameters but don't add layers until AFTER gheatmap
- Add `geom_nodepoint` layers at the end of the function

---

## Phase 9: Integration with Highlighting/Ellipses

**Goal**: Ensure highlighting ellipses work with heatmap

**Note**: This is handled separately after heatmap section in `func_highlight()`

---

## Testing Checklist

After each phase, test:

1. **Basic Tree Only** (no heatmap)
   - [ ] Tree renders correctly
   - [ ] Classifications colored correctly
   - [ ] Bootstrap values shown if enabled

2. **Tree + Single Discrete Heatmap**
   - [ ] Heatmap tiles visible
   - [ ] Colors map to categorical values
   - [ ] Legend shows categories

3. **Tree + Single Continuous Heatmap**
   - [ ] Heatmap tiles visible
   - [ ] Gradient colors correct
   - [ ] Legend shows scale

4. **Tree + Multiple Heatmaps**
   - [ ] Both heatmaps visible
   - [ ] Separate legends for each
   - [ ] No overlap issues

5. **All Features Combined**
   - [ ] Tree + heatmap + bootstrap + highlighting

---

## Reference: Original Working Code

From the original LineagePlotter (non-Shiny version):
https://github.com/yaaraneumeier/LineagePlotter/blob/main/LineagePlotter

Key patterns that worked:
```r
# Data preparation
df_heat_temp <- readfile440[, c(title.id, l_titles_for_heat)]
df_heat_temp <- df_heat_temp[match(tip_list, df_heat_temp[[title.id]]), ]
rownames(df_heat_temp) <- paste0(id_tip_prefix, df_heat_temp[[title.id]])

# gheatmap call
pr440_short_tips_TRY_heat <- gheatmap(tt, data = dxdf440_for_heat[[j1]],
                                      colnames_angle = colnames_angle,
                                      offset = new_heat_x,
                                      width = wi,
                                      font.size = size_font_heat_map_legend,
                                      colnames_offset_y = heat_names_offset,
                                      legend_title = heat_map_title_list[[j1]])
```

---

## Notes for Claude

When asked to work on heatmap:
1. Check which phase we're on (look at version number)
2. Read this plan to understand context
3. Make incremental changes, one phase at a time
4. Update version number after each change
5. Run tests to verify before moving to next phase

The complex old code is preserved in `if (FALSE) { ... }` block starting at line ~4631 in EZlineagePlotter56.R for reference.
