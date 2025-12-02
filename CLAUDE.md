# EZLineagePlotter - Project Memory

## Version Update Requirements

**IMPORTANT**: Every time changes are made to the app, you MUST:

1. Update the version box in `EZlineagePlotter56.R` (around line 4912)
2. Increment the version number (e.g., v57 -> v58)
3. Update the release notes in the version box to describe what changed
4. Tell the user the updated version name
5. Provide the git pull command to get the new version:
   ```bash
   git pull origin claude/analyze-r-code-014ywfqG5i9T5fPLFitF49fG
   ```
6. Provide the curl command to download the R file directly:
   ```bash
   curl -O https://raw.githubusercontent.com/yaaraneumeier/EZLineagePlotter/claude/analyze-r-code-014ywfqG5i9T5fPLFitF49fG/EZlineagePlotter56.R
   ```

## Current Version
**v115**

## Version Box Location
The version box is in the `tabItems` section under `tabName = "data_upload"`, approximately lines 6454-6475.
