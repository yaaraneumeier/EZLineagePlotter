# Performance Optimization Analysis for EZLineagePlotter

**Date:** December 2024
**Based on:** Profiling data from EZlineagePlotter56.R

---

## Current Performance Baseline

| Operation | Time | % of Total |
|-----------|------|------------|
| CSV load | 0.36 sec | 11% |
| Tree plot generation (`func.make.plot.tree.heat.NEW`) | 1.4-1.7 sec | 50% |
| ggsave (SVG) | 0.2-0.3 sec | 8% |
| Post-plot operations (gc, status) | ~1 sec | 30% |
| **Total `generate_plot()`** | **~3.3 sec** | 100% |

---

## Identified Bottlenecks

### BOTTLENECK 1: `p_list_of_pairs` Calculation (0.4-0.5 sec per plot)

**Location:** `func.create.p_list_of_pairs()` at lines 1280-1344

**What it does:**
- Loops through every internal node (`cc_nodss`)
- For each node, calls `ggtree:::getSubtree()` to get all descendant nodes
- Builds a 2x2 contingency table
- Runs `fisher.test()` for each node

**Why it's slow:**
- For a tree with ~850 internal nodes, this means ~850 `getSubtree()` calls and ~850 `fisher.test()` calls
- Each `getSubtree()` traverses the tree structure
- Each `fisher.test()` has R function call overhead

**Optimization Options:**

| Option | What to do | Expected Savings | Risk | Effort |
|--------|-----------|------------------|------|--------|
| **1A: Caching** | Cache results when tree/classifications don't change. Only recalculate when `tree`, `list_node_by_class`, or `FDR_perc` change. Visual-only changes (legend size, colors, etc.) would skip this entirely. | **0.4-0.5 sec** on visual-only updates | Low | Medium |
| **1B: Precompute subtrees** | Compute all subtrees once when tree loads (in `values$tree_subtrees`), then look up during p-value calculation instead of calling `getSubtree()` repeatedly | **0.1-0.2 sec** | Low | Medium |
| **1C: Vectorize** | Replace the R loop with vectorized operations using `sapply`/`vapply` or `data.table` for the contingency table building | **0.1-0.15 sec** | Medium | High |

---

### BOTTLENECK 2: Multiple `ggtree()` Calls (0.25-0.3 sec total)

**Locations:**
- Line 5253: `ggtree(tree440)` → 0.08-0.1 sec
- Line 5324: `ggtree(tree_with_group_CPY)` just to get `$data$group` levels → 0.08-0.1 sec
- Line 5331: `ggtree(tree_with_group_CPY, aes(...))` → 0.085-0.1 sec

**Why it's wasteful:**
- The first `ggtree(tree440)` is called every plot, but the tree rarely changes
- The second call (line 5324) only extracts `levels(...)` - doesn't need full ggtree

**Optimization Options:**

| Option | What to do | Expected Savings | Risk | Effort |
|--------|-----------|------------------|------|--------|
| **2A: Cache tree data** | Store `ggtree(tree440)$data` when tree loads (already in `values$tree_data`), reuse it in `func.make.plot.tree.heat.NEW` instead of recalculating | **0.08-0.1 sec** | Low | Low |
| **2B: Optimize levels extraction** | Replace `levels(ggtree(tree_with_group_CPY)$data$group)` with `levels(groupOTU(...)@data$group)` or extract from `tree_with_group_CPY` directly without ggtree call | **0.08-0.1 sec** | Low | Low |

---

### BOTTLENECK 3: Redundant Full Recalculation

**Current behavior:**
Every plot update (even changing legend font size) triggers the full pipeline:
- Read CSV/YAML (0.01 sec)
- Build tree structures (0.3 sec)
- Calculate p-values (0.4-0.5 sec)
- Build ggplot layers (0.4 sec)
- Save to SVG (0.2 sec)

**The key insight:**
Changes fall into two categories:
1. **Data changes**: Tree file, CSV file, classifications, FDR threshold → require full recalculation
2. **Visual changes**: Legend size, colors, font sizes, scaling → only need re-rendering

**Optimization Option:**

| Option | What to do | Expected Savings | Risk | Effort |
|--------|-----------|------------------|------|--------|
| **3A: Two-tier caching** | Store intermediate results (`p_list_of_pairs`, `tree_TRY` with p-values assigned) in `values$`. Track what inputs changed. If only visual inputs changed, skip statistical calculations and reuse cached tree. | **1.0-1.5 sec** on visual-only updates | Medium | High |

**Detailed breakdown for 3A:**
- Store `values$cached_tree_TRY` after p-values are assigned (line 5385)
- Store `values$cached_inputs_hash` = hash of (tree, csv, FDR, classifications)
- On `generate_plot()`, check if inputs hash matches
- If match → skip to line ~5400 using cached `tree_TRY`
- If different → full recalculation

---

### BOTTLENECK 4: Post-Plot Operations (~1 sec)

**What happens after ggsave:**
- `gc()` garbage collection
- Status updates
- Memory logging

**Optimization Options:**

| Option | What to do | Expected Savings | Risk | Effort |
|--------|-----------|------------------|------|--------|
| **4A: Conditional gc()** | Only call `gc()` every 2-3 plots or when memory exceeds threshold | **0.1-0.3 sec** | Low | Low |
| **4B: Async gc()** | Move gc() to a `later::later()` callback so it doesn't block the UI response | **0.1-0.3 sec** perceived | Low | Medium |

---

## Summary Table: All Options Ranked by Impact

| Priority | Option | Savings | Conditions | Implementation Effort |
|----------|--------|---------|------------|----------------------|
| 1 | 3A: Two-tier caching | **1.0-1.5 sec** | Visual-only changes | High |
| 2 | 1A: Cache p_list_of_pairs | **0.4-0.5 sec** | Visual-only changes | Medium |
| 3 | 2A+2B: Reduce ggtree calls | **0.15-0.2 sec** | All updates | Low |
| 4 | 1B: Precompute subtrees | **0.1-0.2 sec** | All updates | Medium |
| 5 | 4A+4B: Conditional async gc() | **0.1-0.3 sec** | All updates | Low |

---

## Realistic Impact Assessment

**Current total time:** ~3.3 sec per plot

**If you implement options 1A + 2A + 2B + 4A:**
- Data changes: ~2.8-3.0 sec (modest improvement)
- Visual-only changes: **~1.5-2.0 sec** (significant improvement)

**If you implement the full 3A (two-tier caching):**
- Data changes: ~2.8-3.0 sec
- Visual-only changes: **~0.8-1.2 sec** (major improvement)

---

## Implementation Status

| Option | Status | Date | Notes |
|--------|--------|------|-------|
| 4A+4B: Async conditional gc() | **IMPLEMENTED** | Dec 2024 | gc() runs every 3 plots, async via `later::later()` |
| 3A: Two-tier caching | **IMPLEMENTED** | Dec 2024 | Caches `p_list_of_pairs` with hash-based invalidation |
| 2A: Cache tree data | Pending | | |
| 2B: Optimize levels extraction | Pending | | |
| 1A: Cache p_list_of_pairs | **SUPERSEDED** | Dec 2024 | Covered by 3A implementation |

---

## Implementation Details

### Option 3A: Two-tier Caching (IMPLEMENTED)

**What was implemented:**
1. **Cache storage** in `values$`:
   - `values$cached_p_list_of_pairs` - The cached Fisher test p-values
   - `values$cached_p_list_hash` - MD5 hash for cache validation
   - `values$cached_classification_column` - The classification column used (for debugging)

2. **Hash function** `func.create.p_list_cache_hash()`:
   - Computes MD5 hash from: tree structure (tips, edges), classification mapping, FDR threshold, simulate.p.value
   - **CRITICAL**: Classification changes invalidate the cache (as requested)

3. **Cache check** in `func.make.plot.tree.heat.NEW()`:
   - Computes current hash from inputs
   - If hash matches cached hash AND cached data exists → use cached p_list_of_pairs
   - If mismatch → recalculate and return new cache data

4. **Cache flow** through functions:
   - `generate_plot()` passes cached data to `func.print.lineage.tree()`
   - `func.print.lineage.tree()` passes to `func.make.plot.tree.heat.NEW()`
   - New cache data is returned and stored in `values$`

**Log messages to look for:**
- `[PERF-CACHE] Using cached p_list_of_pairs (hash: XXXXXXXX)` - Cache hit
- `[PERF-CACHE] Cache INVALIDATED - hash changed from XXXXXXXX to YYYYYYYY` - Cache miss (data changed)
- `[PERF-CACHE] No cache available - computing p_list_of_pairs` - First run
- `[PERF-CACHE] Stored cache (hash: XXXXXXXX, classification: column_name)` - Cache stored

**Cache invalidation triggers:**
- Tree file change
- Classification column change (column name or values)
- FDR threshold change
- simulate.p.value toggle

**Visual-only changes that use cache:**
- Legend font sizes
- Legend colors
- Output dimensions
- Bootstrap display settings
- Heatmap visual settings
- Node number display

---

## Code References

- `func.create.p_list_cache_hash()`: lines ~1283-1306
- `func.create.p_list_of_pairs()`: lines ~1310-1375
- `func.make.plot.tree.heat.NEW()`: lines ~5230-8200 (cache check at ~5394-5426)
- `func.print.lineage.tree()`: lines ~2877-5227
- `generate_plot()`: lines ~16382-17600
- Cache storage in values$: line ~9397-9403
- gc() call location: line ~17640
