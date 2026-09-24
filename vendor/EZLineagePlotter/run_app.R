# ============================================================================
# EZLineagePlotter Launcher
# Run this file to choose between Single Tree and Multiple Trees mode.
# Each mode is a completely separate app — no interference between them.
# ============================================================================

mode <- readline(prompt = "Enter mode (1 = Single Tree, 2 = Multiple Trees): ")

if (mode == "2") {
  cat("Launching Multiple Trees mode...\n")
  source(file.path(getwd(), "EZlineagePlotter56_multi.R"), local = TRUE)
} else {
  cat("Launching Single Tree mode...\n")
  source(file.path(getwd(), "EZlineagePlotter56.R"), local = TRUE)
}
