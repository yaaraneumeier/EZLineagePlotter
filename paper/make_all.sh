#!/usr/bin/env bash
# make_all.sh - rebuild every paper tree figure from scratch.
#
#   bash paper/make_all.sh [IND ...]     (from anywhere; default: all 9)
#
# 1. render_paper.R <IND>   app render with the paper settings (parallel, ~3-4 min)
#      -> out/paper/cache/<IND>_plot.rds, <IND>_paper_config.yaml, <IND>_app.pdf
# 2. paper_fig.R <IND> final --pdf   production layout (~1 min)
#      -> out/paper/<IND>_final.{png,pdf}
# 3. collect into out/paper/final/: figures, the exact configs used, logs
# Environment: env/activate.sh (conda env EZ_CONDA_ENV, default ez_headless;
# EZ_CONDA_ENV=none to use the current R) + vendored fonts. See README "Deploy".
set -eo pipefail
cd "$(dirname "$0")/.."
source env/activate.sh
set -u          # after activate: conda's hooks use unset variables
Rscript env/check_env.R

ALL=(BRCA-795 BRCA-775 MM-127 MM-412 MM-423 NSCLC-0267 NSCLC-0401 NSCLC-0545 NSCLC-2680)
if [ $# -gt 0 ]; then INDS=("$@"); else INDS=("${ALL[@]}"); fi
LOG=out/paper/final/logs
mkdir -p out/paper/final/configs "$LOG"

for i in "${INDS[@]}"; do
  ( Rscript paper/render_paper.R "$i" > "$LOG/render_$i.log" 2>&1 \
    && Rscript paper/paper_fig.R "$i" final --pdf > "$LOG/fig_$i.log" 2>&1 \
    && echo "[make_all] $i ok" || echo "[make_all] $i FAILED (see $LOG)" ) &
done
wait

for i in "${INDS[@]}"; do
  cp "out/paper/${i}_final.png" "out/paper/${i}_final.pdf" out/paper/final/ 2>/dev/null || true
  cp "out/paper/cache/${i}_paper_config.yaml" out/paper/final/configs/ 2>/dev/null || true
  grep -hE '^\[render_paper\]|^\[paper\]' "$LOG/render_$i.log" "$LOG/fig_$i.log" > "$LOG/summary_$i.txt" || true
done
ls -la out/paper/final
