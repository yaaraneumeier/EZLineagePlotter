# headless — render EZLineagePlotter figures without the UI

Re-creates a figure from a saved `*_config.yaml` (the app's "Download
Configuration") plus the original input files, with no browser. It drives the
app's own `server()` through `shiny::testServer`, replaying what a user does in
Single Tree mode: upload tree → CSV → (CNV RData → chromosome-map RData) → YAML,
Process Data, add a classification, generate, download. No plotting code is
copied, so the output follows the app version in this repository.

## Usage

```bash
Rscript headless/render.R \
  --yaml  MyTree_config.yaml \
  --tree  tree.newick \
  --csv   mapping.csv \
  --out   Tree___MyTree.pdf \
  --rdata results_CNV_tool_final.RData \   # only if a heatmap has data_source: rdata
  --annot Data_Seg_Copy.RData \            # Annot table -> CNV chromosome lines/labels
  --patch patch.yaml \                     # optional, see below
  --class-column Cell.type.2               # tree colouring column (default Cell.type.2)
```

From R:

```r
source("headless/ez_headless.R")          # or Sys.setenv(EZ_APP_DIR = "<repo>") first
res <- ez_render(yaml, tree, csv, "out/Tree___x.pdf", rdata = NULL, annot = NULL,
                 patch = NULL, class_column = "Cell.type.2", class_title = "Cell type", seed = 1)
```

Needs the app's packages (see the top-level README), plus `xml2`.
`aricode` is not needed (it is loaded by the app but unused).

## What the saved config does not contain

| Missing from the YAML | How `ez_render` handles it |
|---|---|
| tree / CSV / RData paths (they point at the upload's temp files) | passed as arguments |
| the tree classification (`classification: []` unless one was *added*, not just previewed) | `class_column` coloured with `PAPER_PALETTE`; values without a palette colour get the app's rainbow default and are returned in `res$unmatched_classes` |
| colours left at their palette default in a "custom colours" discrete heatmap (only changed colours are saved) | filled from the palette by sorted value position, as the colour pickers do before Apply |
| heatmap row-label mapping, CNV chromosome line/label settings (older app versions) | supply them with `--patch` |
| chromosome boundaries for RData CNV heatmaps | `--annot` (an RData with an `Annot` data frame with a `Chr` column) |

Patch file format — 1-based heatmap index → fields, using the field names the
app's YAML import reads, plus optional whole heatmaps to insert
(`insert_heatmaps: [{at: 1, heatmap: {...}}]`, applied after the field edits);
merged into a temporary copy of the config. `ez_render()` also returns the
rendered ggplot (`res$plot`) and left-to-right tip order (`res$tip_order`).

```yaml
heatmaps:
  1: {label_mapping: {ZCA: Cluster}}
  9:
    cnv_chr_lines: 'yes'
    cnv_chr_labels: 'yes'
    cnv_chr_label_position: left   # heatmap is drawn flipped: 'left' ends up on the right
    cnv_chr_label_angle: 0
    cnv_chr_label_size: 2.5
```

## Mock-session workarounds (in `ez_headless.R`)

| Behaviour under `testServer` | Workaround |
|---|---|
| inputs start `NULL`, not at their UI defaults | defaults read from the rendered `ui` HTML and set first |
| `update*Input()` does not change `input` | `session$sendInputMessage` is replaced to feed values back via `setInputs()` |
| action-button observers (`process_data`, `add_classification`) do not fire | `match_tree_with_csv()` is called directly; the classification is stored in `values$classifications` as the handler does |
| the YAML import drops `cnv_chr_*` when building the render list | copied from `values$heatmap_configs` |
| "use all data" makes the individual selectize pick the first CSV individual and rename the plot | that update is dropped |
| the engine also `ggsave()`s a copy to `./` | runs in a temporary working directory |

## Notes

- The Fisher tests behind branch widths use `simulate.p.value = TRUE`;
  `seed` makes a run reproducible, not identical to an interactive session.
- The classification legend title is the classification title ("Cell type" by
  default); figures made before S3.15 show `new_class` instead.
- Runtime is dominated by the classification tests: roughly 30–90 s per figure
  for 50–450 tips.
