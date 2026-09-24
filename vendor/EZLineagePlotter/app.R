# ============================================================================
# app.R - Entry point for Shiny Server / Posit Connect / shinyapps.io.
#
# Runs the FULL EZLineagePlotter app (single-tree AND multiple-trees modes).
#   - Single Tree mode is the default view.
#   - Multiple Trees mode: add ?mode=multi to the URL, or use the in-app
#     mode-switch button in the corner.
#
# All the .R files (this file plus EZlineagePlotter56.R, EZlineagePlotter56_mt.R
# and EZlineagePlotter56_combined.R) must live in the SAME directory; the app
# locates them via the working directory.
#
# The line below evaluates the combined launcher (which ends in
# shinyApp(ui, server)) and returns that app object for the server to run.
# ============================================================================

source("EZlineagePlotter56_combined.R")$value
