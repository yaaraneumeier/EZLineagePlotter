#!/bin/bash
# macOS launcher. Double-click in Finder (you may need to run, once:
#   chmod +x run-ezlineageplotter.command
# and the first time, right-click -> Open to get past Gatekeeper).
echo "============================================================"
echo "  EZLineagePlotter"
echo "============================================================"
echo
echo "Make sure Docker Desktop is running (whale icon in the menu bar)."
echo
echo "Getting the latest version..."
echo "(The FIRST run downloads a few GB - this is normal and happens once.)"
echo
if ! docker pull ghcr.io/yaaraneumeier/ezlineageplotter:latest; then
  echo
  echo "Could not reach Docker. Is Docker Desktop running?"
  read -r -p "Press Enter to close."
  exit 1
fi
echo
echo "Starting the app. Your browser will open in a few seconds."
echo "If it does not, open this address yourself:  http://localhost:3838"
echo
echo "Keep this window open while you use the app. Press Ctrl+C to stop."
( sleep 25 && open http://localhost:3838 ) &
docker run --rm -p 3838:3838 ghcr.io/yaaraneumeier/ezlineageplotter:latest
