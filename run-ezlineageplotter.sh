#!/bin/bash
# Linux launcher. Make it executable once:  chmod +x run-ezlineageplotter.sh
# then run:  ./run-ezlineageplotter.sh
echo "============================================================"
echo "  EZLineagePlotter"
echo "============================================================"
echo
echo "Getting the latest version (first run downloads a few GB, once)..."
if ! docker pull ghcr.io/yaaraneumeier/ezlineageplotter:latest; then
  echo "Could not reach Docker. Is the Docker service running?"
  exit 1
fi
echo
echo "Starting the app. Open this address in your browser:  http://localhost:3838"
echo "Keep this terminal open while you use the app. Press Ctrl+C to stop."
( sleep 20 && (xdg-open http://localhost:3838 >/dev/null 2>&1 || true) ) &
docker run --rm -p 3838:3838 ghcr.io/yaaraneumeier/ezlineageplotter:latest
