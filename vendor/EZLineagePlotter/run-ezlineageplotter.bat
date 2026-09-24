@echo off
title EZLineagePlotter
echo ============================================================
echo   EZLineagePlotter
echo ============================================================
echo.
echo Make sure Docker Desktop is running (whale icon in the tray).
echo.
echo Getting the latest version...
echo (The FIRST run downloads a few GB - this is normal and happens once.)
echo.
docker pull ghcr.io/yaaraneumeier/ezlineageplotter:latest
if errorlevel 1 (
  echo.
  echo Could not reach Docker. Is Docker Desktop running?
  echo.
  pause
  exit /b 1
)
echo.
echo Starting the app. Your browser will open in a few seconds.
echo If it does not, open this address yourself:  http://localhost:3838
echo.
echo Keep THIS window open while you use the app.
echo To stop the app: close this window (or press Ctrl+C).
echo.
start "" /b cmd /c "ping -n 25 127.0.0.1 >nul && start http://localhost:3838"
docker run --rm -p 3838:3838 ghcr.io/yaaraneumeier/ezlineageplotter:latest
