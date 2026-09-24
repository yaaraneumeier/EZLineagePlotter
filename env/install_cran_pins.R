# CRAN packages installed on top of the conda env, at the versions the figures
# were made with (conda-forge's R 4.3 builds stop at ggplot2 3.5.1; the layout
# code relies on ggplot2 4 / S7). Order matters: dependencies first.
# shinyBS: only on the Anaconda "r" channel, not conda-forge.
pins <- c(S7 = "0.2.1", scales = "1.4.0", gtable = "0.3.6", ggplot2 = "4.0.2", shinyBS = "0.61.1")
repos <- Sys.getenv("CRAN_REPO", "https://cloud.r-project.org")
for (p in names(pins)) {
  have <- tryCatch(as.character(packageVersion(p)), error = function(e) "")
  if (identical(have, pins[[p]])) { message(p, " ", have, " ok"); next }
  remotes::install_version(p, version = pins[[p]], repos = repos, upgrade = "never", dependencies = FALSE)
}
