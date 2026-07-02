# Comparison Mode — Implementation Plan

Adds a third app mode (**Comparison**) alongside Single Tree and Multiple Trees,
for comparing **two trees** as a **tanglegram** (facing trees + tip-connecting
lines) with an entanglement score and optional branch **untangling**. Ported and
cleaned up from the original notebook (`func.make.compare_fig`,
`func.make.X.score`, `func_tree_rotation`, …). Batch comparison of many pairs is
intentionally **out of scope for v1**.

## Guiding decisions (all confirmed with the user)

- **New mode, not a tab** — a third mode reachable via the corner switch /
  `?mode=compare`. Implemented as its own module (`EZlineagePlotter56_cmp.R`),
  mirroring the multi-tree module. The single-mode file is never modified.
- **Coloring**: trees start **gray**. A Classification tab (like the other modes)
  applies **one shared column + one shared palette to both trees**. The coloring
  itself uses the **notebook's method** (color by class; not every branch) — this
  is deliberate, it looked better.
- **Tip matching** uses the app's existing **ID mapping** (trim/prefix/ID column),
  which is more general than raw Newick labels.
- **Mismatched tips** are handled in the **Upload tab**: after uploading the two
  trees + CSV, a match-check step reports matched / mismatched tips and lets the
  user **prune to shared tips** or **keep all**.
- **No ladderize, ever** — rotation to match the trees is the whole point.
- **Untangle** is **off by default**, **button-triggered** (never auto), with
  options: which tree to rotate (A / B / Both), number of passes, and polytomy
  thoroughness (`PERM_BOUND`, only when multifurcations exist). It is slow, so it
  shows progress and a "this can take a while" message.
- **Versions**: the Original plus each untangle run are stored as **versions** the
  user can switch between instantly (no recompute). Each version is auto-labeled
  from its parameters and **renamable** and **deletable** ("Clear all" too).
- **Scores** are always reported **relative to the Original baseline**
  (`X_original → X_version`, normalized by tips², plus % reduction).
- **Untangle result vs coloring are separated**: a version stores only the
  **rotated ordering** (as Newick); classification coloring is applied at render.
  So recoloring is free and never invalidates versions. Uploading new trees or
  changing the prune choice clears versions (topology changed).

## Persistence & downloads

- **Save/Load** works like the other modes: Upload + Download + **Config** tabs.
- The config YAML has a `comparison` section that embeds, per input tree,
  **path + original name + Newick**, plus the shared mapping/classification, the
  prune choice, and the **versions** (each with its rotated Newick + scores +
  label). Loading restores versions instantly (Newick read back, no ladderize).
- **Size guard**: if the stored versions exceed a size threshold, show a message
  and block adding new versions until the user deletes some. (Large configs are
  not the standard case today; this just prevents runaway files.)
- **Download tab** exports, for the currently displayed version: the **plot
  (PDF/PNG)** and the **rotated Newick(s)** (reusing the app's existing rotated-
  tree Newick download).

## Draft `comparison` YAML shape

```
comparison:
  trees:
    A: { path, name, newick }
    B: { path, name, newick }
  csv:            # handled like the other modes
  mapping:        { id_column, id_tip_trim_flag, start, end, prefix }
  classification: { column, palette }
  prune:          shared | keep
  versions:
    - { id, label, params: { tree, passes, polytomy }, score_original, score_version, newick_A, newick_B }
  current_version: <id>
```

## Files

- **New `EZlineagePlotter56_cmp.R`** — module: `cmp_tabItem_*()` UI builders,
  `cmp_install_server()`, and the ported+cleaned comparison engine.
- **Edit `EZlineagePlotter56_combined.R`** — source the module, build
  `compare_ui`, make `ui()`/`server()` 3-way, extend the corner switch to cycle
  Single → Multiple → Comparison.
- **No Dockerfile change** — engine needs only `ape`, `ggtree`, `ggnewscale`,
  `scales`, `tidyverse`, all already in the image.

## Critical technical note — name collisions

The notebook engine and the app share some `func.*` names with different
signatures (e.g. `func.rotate.specific.nodes`,
`func.make.leaves_id_ordered_for_df440`). `combined.R` evals the single-mode file
into the global environment, so the comparison engine must be **namespaced**
(prefix `cmp.` and/or sourced `local = TRUE`) to avoid clobbering single-mode
functions. This is the main correctness risk.

## Data model

A **version** = rotated ordering, not a picture:
`{ id, label, params, score_original, score_version, newick_A, newick_B }`, held
in `reactiveValues`. Rendering = read ordering → `ggtree` (no ladderize) → apply
current classification coloring → tanglegram. Untangle (expensive, occasional)
and coloring (cheap, frequent) stay separate.

## Engine port (into the module, namespaced `cmp.`)

- `cmp.make.compare.fig` (from `func.make.compare_fig`) — restructured to build
  the tanglegram ggplot from two ordered + colored trees and **return the ggplot**
  (no `ggsave`, no debug prints).
- `cmp.make.X.score` + `cmp.check.lines.intersect` — the crossing score.
- `cmp.tree.rotation` + helper cluster — the untangler, parameterized by
  which-tree / passes / polytomy bound.
- Prune-to-shared helper, driven by the Upload choice + ID mapping.

## Tabs

1. **Upload** — Tree A, Tree B, one CSV + ID-mapping widgets → match-check panel
   → prune-or-keep. Builds "Original" (v0, gray).
2. **Classification** — shared column + palette; notebook-style coloring; gray by
   default.
3. **Compare / Untangle** — tanglegram preview; untangle controls + "Untangle
   now"; Versions panel (select / rename / delete / clear); score readout;
   progress messaging.
4. **Download** — current version's plot (PDF/PNG) + rotated Newick(s).
5. **Config** — Save/Load the `comparison` YAML.

## Milestones (each independently testable)

1. **Module skeleton + 3-way mode switch** — empty Comparison mode loads and is
   reachable via the corner button / `?mode=compare`. Touches only the new module
   and `combined.R`; changes no existing behavior. ← *in progress*
2. Upload + match-check + prune/keep; gray "Original" tanglegram preview.
3. Classification (shared column/palette, notebook coloring).
4. Untangle engine port → run → score → auto-switch to result.
5. Versions store (rename/delete/clear, baseline % scores, size guard).
6. Download (plot + rotated Newick).
7. YAML save/load of comparison + versions.
8. Docker test → version bump → PR → build → release.
