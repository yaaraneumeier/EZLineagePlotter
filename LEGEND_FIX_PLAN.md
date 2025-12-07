# EZLineagePlotter Legend Fix Plan (v151)

## Issue 1: Legend Ellipse Transparency NOT Working

### Root Cause Found:

**PLOT ellipse (line 1496, 1521):**
```r
fill = high_color_list[[1]], alpha = alpha_val, linetype = "blank"
```
Uses `alpha = alpha_val` as a **SEPARATE PARAMETER** to geom_ellipse.

**LEGEND ellipse (line 1768):**
```r
fill = scales::alpha(high_color_list[[index_high]], current_alpha), colour = NA, linetype = "blank"
```
Uses `scales::alpha()` to **EMBED alpha in fill color** - this does NOT work the same way!

### The Fix:
Change legend ellipse to use `alpha = current_alpha` as a separate parameter, EXACTLY like the plot ellipse.

**Before (line 1768):**
```r
fill = scales::alpha(high_color_list[[index_high]], current_alpha), colour = NA, linetype = "blank", show.legend = FALSE
```

**After:**
```r
fill = high_color_list[[index_high]], alpha = current_alpha, colour = NA, linetype = "blank", show.legend = FALSE
```

---

## Issue 2: Legend Alignment - Expanding Plot Range

### Root Cause Found:

From debug output:
```
Panel 1 x.range: [-95.525, 4.025]  <- HUGE range!
Bootstrap Legend: Title position: x=-1.44, y=91
Highlight Legend: Label position: x=1.28, y=80.8, Ellipse position: y=86
```

**The Problem:**
- Tree has 80 tips, so y-range should be ~[1, 80]
- Bootstrap legend at y=91 EXPANDS the plot's data range
- With `coord_flip`, y becomes visual x (horizontal)
- This makes the tree appear tiny because the panel range is expanded to accommodate legends

### Current Legend Positioning (lines 1637-1647, 1789-1795):

```r
# y_off_base calculation (line 6621-6627):
n_tips <- 80
y_off_base <- n_tips + 5  # = 85 (for RIGHT side)

# new_base_for_second_legend_normalized = 3 (for x position)

# Highlight legend:
label_x = 3 + extra  # = ~1.28
ellipse_y = y_off_base + 1 + adjustments  # = 86+

# Bootstrap legend:
boot_y_base = y_off_base + bootstrap_y_offset  # = 85+
boot_title_y = boot_y_base + 6  # = 91
boot_triangles_y = boot_title_y - 2  # = 89
```

### The Problem Explained:

With coord_flip + scale_y_reverse:
- Data `x` -> visual vertical (tree depth)
- Data `y` -> visual horizontal (tip positions)

When legend y-values exceed n_tips (80), they EXPAND the plot's y-range, which after coord_flip expands the HORIZONTAL range. This shrinks the tree.

### The Fix Strategy:

**Option A: Use annotation_custom() with grob objects**
- Place legends OUTSIDE the plot panel using gtable/grid coordinates
- Legends won't affect plot data range
- Most robust solution

**Option B: Constrain legend y-values to within tip range**
- Instead of y_off_base = n_tips + 5 = 85
- Use negative values or small offsets from tree
- Position legends at y = -5 to -20 (below tip 1)
- Then use coord_cartesian(clip = "off") to show them

**Option C: Use ggplot's native legend system**
- Convert highlight/bootstrap to proper ggplot scales
- Let ggplot position them in the legend box
- Most integrated solution but requires significant refactoring

### Recommended Fix (Option B - quickest):

1. Change y_off_base to position legends BELOW the tree (negative y):
   ```r
   y_off_base <- -5  # Below tree tips (tip 1 is at y=1)
   ```

2. Add coord_cartesian(clip = "off") to allow rendering outside panel

3. Adjust x positions to align vertically with ggplot legends

### Detailed Code Changes for v151:

#### Change 1: Fix legend ellipse transparency (EZlineagePlotter56.R line 1768)
```r
# BEFORE:
fill = scales::alpha(high_color_list[[index_high]], current_alpha), colour = NA, linetype = "blank"

# AFTER:
fill = high_color_list[[index_high]], alpha = current_alpha, colour = NA, linetype = "blank"
```

#### Change 2: Fix y_off_base calculation (line 6621-6627)
```r
# BEFORE:
y_off_base <- n_tips + 5  # Position beyond the rightmost tips

# AFTER:
# Position legends BELOW tree (at negative y, which is LEFT side visually after coord_flip)
# OR use positive small values that won't expand range significantly
y_off_base <- max(p$data$y, na.rm=TRUE)  # Stay within existing range
```

#### Change 3: Add coord_cartesian clip setting
Need to ensure legends can render outside panel if needed.

---

## Debug Output Additions Needed

Add explicit output for each legend element:
```r
cat(file=stderr(), "\n=== v151: LEGEND ELEMENT COORDINATES ===\n")
cat(file=stderr(), "HIGHLIGHT LEGEND:\n")
cat(file=stderr(), paste0("  Title: x=", title_x, ", y=", title_y, "\n"))
cat(file=stderr(), paste0("  Ellipse: x=", ellipse_x, ", y=", ellipse_y, ", a=", a, ", b=", b, "\n"))
cat(file=stderr(), paste0("  Label: x=", label_x, ", y=", label_y, "\n"))
cat(file=stderr(), "BOOTSTRAP LEGEND:\n")
cat(file=stderr(), paste0("  Title: x=", boot_title_x, ", y=", boot_title_y, "\n"))
cat(file=stderr(), paste0("  Triangles: x=", tri_x, ", y=", boot_triangles_y, "\n"))
cat(file=stderr(), paste0("  Labels: x=", lab_x, ", y=", boot_labels_y, "\n"))
cat(file=stderr(), "TREE DATA RANGE:\n")
cat(file=stderr(), paste0("  y-range: ", min_y, " to ", max_y, "\n"))
cat(file=stderr(), paste0("  x-range: ", min_x, " to ", max_x, "\n"))
cat(file=stderr(), "=========================================\n")
```

---

## Summary of Changes for v151

1. **Line 1768**: Change `fill = scales::alpha(...)` to `fill = ..., alpha = current_alpha`
2. **Lines 6621-6627**: Recalculate y_off_base to NOT exceed tip range
3. **Lines 1789-1795**: Adjust bootstrap y-positions to stay within range
4. **Add comprehensive debug output** for all legend elements
5. **Update version** to v151 in header and info box
