# source env/activate.sh   (from anywhere) - activates the conda env and the vendored fonts.
# EZ_CONDA_ENV: env name or prefix path (default ez_headless); EZ_CONDA_ENV=none uses the current R.
# CONDA_ROOT: conda/mamba install to use, if not one of the usual locations.
_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
_env="${EZ_CONDA_ENV:-ez_headless}"
if [ "$_env" != "none" ]; then
  _ok=""
  for c in "${CONDA_ROOT:-}" "$(conda info --base 2>/dev/null)" "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
    [ -n "$c" ] && [ -f "$c/etc/profile.d/conda.sh" ] || continue
    source "$c/etc/profile.d/conda.sh"
    if conda activate "$_env" 2>/dev/null; then _ok=1; break; fi
  done
  [ -n "$_ok" ] || { echo "env/activate.sh: cannot activate conda env '$_env' (set EZ_CONDA_ENV / CONDA_ROOT, or EZ_CONDA_ENV=none)" >&2; return 1 2>/dev/null || exit 1; }
fi
export FONTCONFIG_FILE="$_root/fonts/fonts.conf"
export EZ_ROOT="$_root"
unset _root _env _ok
