# EZLineagePlotter

A Shiny application for visualizing and analyzing phylogenetic trees with rich annotation capabilities.

## Version S1.6 (Stable Release)

### What's New in S1.6
- **Performance Optimizations**: Faster rendering with optimized layer management and batched updates
- **Bug Fixes**: Legend background color works correctly, fixed multiple highlights with heatmap issue

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

Debug output is disabled by default for performance. To enable verbose logging, edit line 87 in `EZlineagePlotter56.R`:

```r
DEBUG_VERBOSE <- TRUE
```

## License

This project is provided as-is for research and educational purposes.

## Author

Yaara Neumeier
