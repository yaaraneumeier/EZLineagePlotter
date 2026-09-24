# EZLineagePlotter

A Shiny application for visualizing and analyzing phylogenetic trees with rich annotation capabilities.

## Version S3.16 (Stable Release)

### What's New in S3.16
- **Extra heatmap legends**: add one or more manual color bars for a continuous/CNV heatmap (same colors) with your own title and tick values/labels, from the Legend tab
- **Side-by-side legends**: a Legend-tab option to arrange legends in a horizontal row (so a heatmap's own legend and its extra legend sit next to each other) instead of stacked
- **Row mask**: gray out whole heatmap rows by a CSV column value, with solid / background / dots styles and a chosen mask color
- **Clade Clusters downloads**: export the tip names in current left-to-right tree order (one per line, comma-separated, or CSV), and a tip→clade mapping (TSV; tips with no clade written `NAN`)
- **Fixes**: heatmap column-name mapping now saves to the config; RData CNV heatmap legend value scale now matches the midpoint-centered colors; RData CNV default colors are now low = blue, high = red

### Earlier highlights
- **S3.15**: Classification 'Legend Title' fix, legend Title/Keys alignment (left/center/right), key-to-label spacing, dev mode off (no hard refresh needed)
- **S3.14**: Clade Clusters from a CSV column, Classification 'Paper palette', up to 15 heatmaps (fixed 7+ color/title), RData CNV custom row label fix, faster CSV loading (`data.table::fread`)
- **S3.13**: Clade Clusters overlay by tip range (bracket / shaded band / separator lines / colored tip strip), Mirror Tree (single + multi mode), multifurcating-node rotation with per-branch left/right ordering, multi-tree manual rotation fix
- **S2.0**: RData CNV heatmaps from QDNAseq/scIMPACT pipelines
- **S1.6**: Performance optimizations and legend/highlight bug fixes

## Features

- **Tree Visualization**: Display phylogenetic trees with customizable layouts
- **Tree Rotation**: Rotate and reorder tree branches - flip clades, rotate nodes, and reroot the tree to optimize visual presentation
- **Classification Coloring**: Color tree tips based on classification data from CSV files
- **Multiple Heatmaps**: Add multiple heatmaps alongside the tree with discrete or continuous color scales
- **Highlight Regions**: Draw ellipses around groups of related tips
- **Bootstrap Values**: Display bootstrap support values as triangles, percentages, or color-coded nodes
- **P-value Display**: Show statistical significance with size-scaled points
- **Flexible Legends**: Customizable legend positioning, styling, and visibility controls
- **Export Options**: Save plots as PDF or PNG with custom dimensions

## Requirements

The following R packages are required:

```r
shiny
shinydashboard
shinyjs
shinyWidgets
shinyBS
colourpicker
DT
ape
dplyr
ggplot2
ggtree
ggforce
yaml
stringr
scales
tidyverse
ggnewscale
gridExtra
cowplot
jpeg
combinat
infotheo
aricode
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yaaraneumeier/EZLineagePlotter.git
cd EZLineagePlotter
```

2. Install required R packages:
```r
install.packages(c("shiny", "shinydashboard", "shinyjs", "shinyWidgets",
                   "shinyBS", "colourpicker", "DT", "ape", "dplyr",
                   "ggplot2", "yaml", "stringr", "scales", "tidyverse",
                   "ggnewscale", "gridExtra", "cowplot", "jpeg",
                   "combinat", "infotheo", "aricode"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")

# ggforce from CRAN
install.packages("ggforce")
```

## Usage

Run the application:
```r
shiny::runApp("EZlineagePlotter56.R")
```

### Files / app structure

The app ships as several files that must live in the **same directory** (they
locate each other via the working directory):

| File | Purpose |
|------|---------|
| `EZlineagePlotter56.R` | Single-tree app (standalone) |
| `EZlineagePlotter56_mt.R` | Multiple-trees module (sourced by the launchers) |
| `EZlineagePlotter56_combined.R` | Both modes in one app (single + multi) |
| `app.R` | Entry point for Shiny Server / Connect / shinyapps.io — runs the **full** app |

### Running the full app (both modes)

From an R session, inside the project folder:
```r
shiny::runApp("EZlineagePlotter56_combined.R")
```
- Single Tree mode is the default.
- Multiple Trees mode: add `?mode=multi` to the URL, or use the in-app
  mode-switch button.

### Running locally with Docker (no R install needed)

For a single user who just wants to run the full app on their own computer
(Windows/macOS/Linux) without installing R or any packages, use the bundled
container image. See **[RUN_WITH_DOCKER.md](RUN_WITH_DOCKER.md)** for the
step-by-step guide. In short, once Docker Desktop is installed:

```bash
docker run --rm -p 3838:3838 ghcr.io/yaaraneumeier/ezlineageplotter:latest
```

then open http://localhost:3838. Double-click launchers
(`run-ezlineageplotter.bat` / `.command` / `.sh`) are provided so non-technical
users don't have to touch a terminal. The image is built and published
automatically by the `Build & publish Docker image` GitHub Action.

### Deploying to a server (Shiny Server / Posit Connect / shinyapps.io)

These platforms run an `app.R` in the application directory. A ready-made
`app.R` is included that launches the **full** app (both modes). To deploy:

1. Put all the `.R` files in one application directory, e.g.
   `/srv/shiny-server/ezlineageplotter/`:
   ```
   ezlineageplotter/
   ├── app.R
   ├── EZlineagePlotter56.R
   ├── EZlineagePlotter56_mt.R
   └── EZlineagePlotter56_combined.R
   ```
2. Ensure the required R packages are installed for the server's R (see
   Installation above; `ggtree`/`treeio` are Bioconductor).
3. Point the server at that directory. Shiny Server picks up `app.R`
   automatically; for shinyapps.io use
   `rsconnect::deployApp("ezlineageplotter")`.

Single Tree is the default URL; append `?mode=multi` for Multiple Trees mode.

### Input Files

1. **Tree File**: Newick format (.nwk, .tree, .newick)
2. **Classification CSV**: CSV file with tip labels and classification data

### Workflow

1. **Data Upload Tab**: Upload your tree and CSV files
2. **Classification Tab**: Configure how tips are colored based on CSV columns
3. **Heatmap Tab**: Add heatmaps to display additional data alongside the tree
4. **Highlight Tab**: Define highlight regions to emphasize groups of tips
5. **Bootstrap Tab**: Configure bootstrap value display
6. **Legend Tab**: Customize legend appearance and positioning
7. **Extra Tab**: Adjust tree scaling, background, and overlays
8. **Download Tab**: Export your plot as PDF or PNG

## Debug Mode

Debug output is disabled by default for performance. To enable verbose logging, set `DEBUG_VERBOSE <- TRUE` near the top of `EZlineagePlotter56.R`:

```r
DEBUG_VERBOSE <- TRUE
```

## License

This project is provided as-is for research and educational purposes.

## Author

Yaara Neumeier
