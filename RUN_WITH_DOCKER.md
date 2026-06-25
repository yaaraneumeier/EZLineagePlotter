# Running EZLineagePlotter on your own computer (with Docker)

This lets you run the **full app** locally in your web browser — exactly like
the server version, but on your own machine. You do **not** need R, RStudio, or
to install any R packages. Everything is bundled in one container.

It works the same on **Windows, macOS, and Linux**.

---

## One-time setup

### 1. Install Docker Desktop
- Download from https://www.docker.com/products/docker-desktop/ and install it.
- Start Docker Desktop and wait until the whale icon says it is **running**.
- (Windows may ask to enable WSL2 / virtualization — accept. This needs
  administrator rights; if you can't install it, ask IT, or use the RStudio
  method in the main README instead.)

### 2. Get the launcher file for your system
From this repository, download the one that matches your computer:
- **Windows:** `run-ezlineageplotter.bat`
- **macOS:** `run-ezlineageplotter.command`
- **Linux:** `run-ezlineageplotter.sh`

Put it anywhere convenient (e.g. your Desktop).

---

## Every time you want to use the app

1. Make sure **Docker Desktop is running**.
2. Double-click the launcher:
   - **Windows:** double-click `run-ezlineageplotter.bat`
   - **macOS:** double-click `run-ezlineageplotter.command`
     (first time only: right-click → **Open** to get past the security prompt;
     if it won't run, open Terminal and run `chmod +x run-ezlineageplotter.command`)
   - **Linux:** `chmod +x run-ezlineageplotter.sh` once, then `./run-ezlineageplotter.sh`
3. A black window opens and the app starts. Your browser opens automatically to
   **http://localhost:3838** after a few seconds. If it doesn't, just open that
   address yourself.
4. Use the app normally — upload your tree/CSV and download your figures through
   the browser, exactly like the server version.

> The **first run downloads a few GB** (the app image). This is normal and
> happens only once. Later runs start in seconds.

### To stop the app
Close the black launcher window (or press **Ctrl+C** in it). Your data isn't
stored in the container — everything is uploaded/downloaded through the browser.

### To update to a newer version
Just run the launcher again — it automatically pulls the latest image.

---

## For a hands-on user (skip the launcher)
Once Docker is installed you can run it directly:

```bash
docker run --rm -p 3838:3838 ghcr.io/yaaraneumeier/ezlineageplotter:latest
```

Then open http://localhost:3838. Use `?mode=multi` in the URL for multiple-trees mode.

---

## Troubleshooting

- **"Could not reach Docker" / nothing happens** — Docker Desktop isn't running.
  Start it, wait for the whale icon to settle, then try again.
- **Browser shows "can't connect" for a few seconds** — the app is still
  starting up (it can take ~10–20 seconds). Refresh.
- **Port 3838 already in use** — something else is using that port. Edit the
  launcher and change both numbers in `-p 3838:3838` to e.g. `-p 8080:3838`,
  then open http://localhost:8080.
- **App feels low on memory with big trees/heatmaps** — give Docker more RAM:
  Docker Desktop → **Settings → Resources** → raise Memory to 4 GB (or more).
