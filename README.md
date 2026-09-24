# ez_headless — EZLineagePlotter figures without the Shiny UI

Re-renders the tree figures (`data/Trees 04-2026`, `data/Trees 07-2026`) from their
saved `*_config.yaml` files with the **unmodified** app code, and builds the paper
figures (`paper/`). The repository is self-contained: app, inputs and fonts are
inside it; nothing is read from outside (verified 2026-09-24: a copy in another
directory rebuilt BRCA-795 / MM-412 / NSCLC-0401 pixel-identical).

## Deploy

```
ez_headless/
  vendor/EZLineagePlotter/   the app, branch headless-render @ 3a1b0c5 (see vendor/*.VENDORED.txt);
                             headless driver in its headless/ folder
  data/                      inputs (layout as in the klein/ project folder):
    masterlist_lineage_tree_v84_20260923.csv, patient history-5.pdf (banners),
    zca_leaf_cluster_map_20260909.csv, trees2025_ordering/*.newick,
    patients_samples_by_name/<patient>/RData/{results_CNV_tool_final*,Data_Seg_Copy}.RData,
    Trees 04-2026/ + Trees 07-2026/  (configs + reference renders for compare.py; the .pptx decks are not copied),
    mutation_labels.tsv  (mutation row labels, from mutation_labels.py)
                             NOT IN GIT (patient data): shipped separately as
                             ez_headless_data_<date>.tar.gz -> extract in the repo root (creates data/),
                             check with sha256sum -c ez_headless_data_<date>.tar.gz.sha256
  fonts/                     Nimbus Sans (URW base35, AGPL + font exception, LICENSE/COPYING) + fonts.conf
  env/                       environment.yml, install_cran_pins.R, check_env.R, activate.sh
  paper/                     paper figures: render_paper.R, paper_fig.R, make_all.sh
  individuals.R, patches.yaml, run_all.R, zca_*.R, compare.py, mutation_labels.py
```

Setup (linux-64; needs conda/mamba and internet once):

```bash
tar -xzf /path/to/ez_headless_data_20260924.tar.gz    # in the repo root -> data/ (patient data, not in git)
conda env create -f env/environment.yml          # env "ez_headless": R 4.3.3, ggtree 3.10, poppler, ...
conda activate ez_headless
Rscript env/install_cran_pins.R                  # ggplot2 4.0.2, S7, scales, gtable, shinyBS from CRAN
source env/activate.sh && Rscript env/check_env.R   # preflight: versions, pdftoppm, fonts, data
bash paper/make_all.sh                           # all 9 paper figures -> out/paper/final/ (~7 min, parallel)
```

`paper/make_all.sh` sources `env/activate.sh` itself (conda env + `FONTCONFIG_FILE=fonts/fonts.conf`)
and runs the preflight first. Overrides:

| Variable | Default | Use |
|---|---|---|
| `EZ_CONDA_ENV` | `ez_headless` | env name or prefix; `none` = use the R already on PATH |
| `CONDA_ROOT` | usual install locations | conda install to use |
| `EZ_ROOT` | current directory | repo root (scripts are run from it) |
| `EZ_APP_DIR` | `vendor/EZLineagePlotter` | another app checkout |
| `EZ_DATA` | `data` | another data folder with the same layout |

Run the other R scripts from the repo root after `source env/activate.sh`. The
bold font matters: without `fonts.conf` (or a system Nimbus Sans) the figures fall
back to a font without a bold face and every "bold" renders regular.
`data/` holds patient-derived data (masterlist, CNV, mutations, trees, clinical timelines) and is
`.gitignore`d; it travels only as the archive.

```bash
source env/activate.sh
Rscript run_all.R              # all 19 configs  (04-2026 ~35 s each, 07-2026 ~80 s each)
Rscript run_all.R MM-127_v82   # subset: configs whose file name contains the text
python3 compare.py             # out/compare/*.png  reference | re-render side by side
Rscript zca_trees.R            # 10 trees (latest config each) with the ZCA Cluster row -> out/ZCA_trees/
Rscript zca_check.R            # drawn Cluster row vs slide map -> out/ZCA_trees/zca_check_*.tsv
Rscript zca_topology.R         # ZCA vs tree topology -> out/ZCA_trees/zca_topology_*.tsv
```

`individuals.R` holds the shared input table (newick / RData per individual).

**Paper figures** (production layout, banners, one template for all patients): `paper/`
(`bash paper/make_all.sh` → `out/paper/final/`). Full documentation: `PAPER_FIGURES_HANDOFF.md` (project notes, kept outside git; `klein/PAPER_FIGURES_HANDOFF.md`).

## ZCA check (2026-09-23)

`zca_check.R` decodes the Cluster row from the rendered ggplot (tile colour →
cluster via the layer's fill scale, tile position → tip label), so it checks
what is drawn, not the masterlist. Result: all 1020 leaves in the 10 trees match
`zca_leaf_cluster_map_20260909.csv`; masterlist v84 `ZCA` equals the 2026-09-12
`ZCA_fixed` column. Leaf order equals the slides except BRCA-795 (the v10 config's
node-81 rotation) and NSCLC-2680 (the slides swap two sibling blocks — the known
exception in ZCA_CLUSTER_HANDOFF.md). The 07-2026 BRCA-795 configs give Cluster
1B and Cluster 2 the same blue; `zca_trees.R` recolours Cluster 2 orange.

Outputs: `out/<Trees folder>/<same file name as the reference figure>`, plus
`out/run_log*.tsv` (status, seconds, rdata used, classes without a paper colour).

## How it works

`headless/ez_headless.R` loads the app (everything except `shinyApp()`) and drives its
real `server()` through `shiny::testServer`, replaying what a user does:
upload tree → CSV → CNV RData → chromosome-map RData → YAML, process data,
add a classification, generate, download. No plotting code is copied, so the
output follows whatever the app version does.

Mock-session quirks it works around:

| Quirk | Workaround |
|---|---|
| inputs start `NULL`, not at their UI defaults | defaults read from the rendered `ui` HTML and set first |
| `update*Input()` does not change `input` | `session$sendInputMessage` is replaced to feed values back via `setInputs()` |
| action-button observers (`process_data`, `add_classification`) do not fire | their bodies are called / reproduced directly (`match_tree_with_csv()`, `values$classifications`) |
| the engine also `ggsave()`s a copy to `./` | runs in a temp working dir |
| with "use all data" the selectize picks the first individual and renames the plot | that one update is dropped |

## What the saved YAMLs do not contain (and where it comes from here)

| Missing from YAML | Source used |
|---|---|
| tree / CSV / RData paths (point at deleted `/tmp` uploads) | `run_all.R` → `INDIVIDUALS`: `trees2025_ordering/*_bootstrapped_rerooted.newick` (tip order and topology identical to the app's `*_rotated.newick` exports), `masterlist_lineage_tree_v84_20260923.csv`, `patients_samples_by_name/<patient>/RData/` |
| tree classification (`classification: []` in every config) | `Cell.type.2` coloured with the app's `PAPER_PALETTE` — the palette the Classification tab offers when a column's values match it |
| chromosome mapping for the CNV heatmap | `Data_Seg_Copy.RData` (`Annot`) from the same patient folder |
| mutation row labels ("ESR1 p.S463P") | `mutation_labels.py` → `mutation_labels.tsv` from `scIMPACT Mutations_Classifier-AllelicCounts_20250116ZCxlsx.xlsx` (SYMBOL + theVariant → `SYMBOL p.AAchange`), applied to every mutation heatmap by `mutation_label_patch()` in `individuals.R`. One correction: NRAS 1:114713908 T>C is Q61R (workbook says Q61P) |
| other heatmap row labels (`label_mapping`) and CNV chromosome lines/labels, 07-2026 configs | `patches.yaml`, read off the reference PDFs (their text is outlined, so by eye). Merged into a temp copy of the config; originals untouched |

Discrete heatmaps with custom colours: the YAML holds only the colours the user
changed (app S2.7 skips saving palette defaults); the driver fills the rest from
the palette by sorted value, as the browser's colour pickers do before Apply.

Also: the app's YAML import drops the `cnv_chr_*` fields when it builds the
render list (in the browser they come back when the heatmap tab is applied);
the driver copies them across.

## Known differences from the references

**04-2026 configs are incomplete and cannot be fully reproduced as-is.** Their
heatmaps 7–10 were saved with empty column lists. From the reference figures
they were: a large continuous (1–5, green) block of CSV columns without row
labels, plus WGD, Date, an NRAS survey row (MM-127) and extra discrete rows
(MM-412, MM-423). The continuous block's columns do not exist in masterlist v84
(no per-bin copy-number columns — only `Absolute CNV`). The re-renders show
heatmaps 1–6 only; tree, tip order, bootstrap and heatmaps 1–6 otherwise match.


- Classification legend title is "Cell type" (current app default) instead of
  `new_class` (the references predate the app's S3.15 legend-title fix), so the
  legend order at the bottom differs.
- 04-2026 references coloured the tree with the older rainbow default; the
  re-renders use the paper palette, as requested.
- ZCA (Cluster row) and other masterlist columns come from v84, so they can
  differ where the masterlist changed since v72/v82 (e.g. the 2026-09-12 ZCA fix).
- Fisher tests run with `simulate.p.value = TRUE`; `set.seed(1)` makes branch
  widths reproducible run to run, but not identical to the original session.
- MM15-127: `results_CNV_tool_final.RData` (not `_v2`) matches the reference CNV
  panel; `_v2` differs in two samples.

## ZCA vs tree topology (2026-09-23)

`zca_topology.R`, on the `trees2025_ordering` newicks. In all 10 trees the ZCA
partition is convex: exact Sankoff parsimony of ZCA = (clusters − 1), i.e. the
clusters are produced by cutting exactly k − 1 edges, each cluster one connected
piece. As drawn, every cluster is an exact clade except the "trunk" cluster
(Cluster 1, or 1A in 795-09 / MM16-412), which is what remains after the nested
subclades are cut off (paraphyletic as drawn, a clean split unrooted). Lettered
families: 1A+1B is a clade in 795-09 and MM16-412; in BC15-0267 2A–2D are a chain
of nested cuts with Cluster 3 at its tip, so "2A–2D" is a group only with Cluster 3.
Weakly supported cut edges (newick support < 0.5): 775-13 1|2 (0.30), 7|8 (0.35);
MM16-412 1A|1B (0.33), 1A|2 (0.22); BC14-2680 2|3 (0.33).
