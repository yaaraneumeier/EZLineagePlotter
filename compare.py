"""Side-by-side reference vs. re-rendered figure for every row of out/run_log.tsv.

    python3 compare.py            -> out/compare/<config stem>.png
"""
import csv, os, subprocess, tempfile
from PIL import Image, ImageDraw

KLEIN = os.environ.get("EZ_DATA", os.path.join(os.environ.get("EZ_ROOT", "."), "data"))  # run from the repo root
OUT = "out/compare"
DPI = 70


def raster(path, dpi=DPI):
    if path.lower().endswith(".png"):
        im = Image.open(path).convert("RGB")
        w = int(29.7 / 2.54 * dpi)  # references are A4 landscape
        return im.resize((w, int(im.height * w / im.width)))
    with tempfile.TemporaryDirectory() as d:
        subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-f", "1", "-l", "1", path, f"{d}/p"], check=True)
        return Image.open(os.path.join(d, sorted(os.listdir(d))[0])).convert("RGB")


os.makedirs(OUT, exist_ok=True)
import glob
rows = [r for f in sorted(glob.glob("out/run_log*.tsv")) for r in csv.DictReader(open(f), delimiter="\t")]
for row in rows:
    if row["ok"] != "TRUE":
        print("skip (failed):", row["config"])
        continue
    refs = [r for r in row["reference"].split(";") if r]
    ext = os.path.splitext(row["output"])[1]
    ref = next((r for r in refs if r.endswith(ext)), refs[0] if refs else None)
    if ref is None:
        print("no reference:", row["config"])
        continue
    a = raster(os.path.join(KLEIN, row["folder"], ref))
    b = raster(row["output"])
    c = Image.new("RGB", (a.width + b.width + 10, max(a.height, b.height) + 22), "white")
    d = ImageDraw.Draw(c)
    d.text((5, 4), "REFERENCE  " + ref, fill="red")
    d.text((a.width + 15, 4), "EZ_HEADLESS  " + os.path.basename(row["output"]), fill="blue")
    c.paste(a, (0, 22))
    c.paste(b, (a.width + 10, 22))
    name = row["config"].replace("_config", "").replace(".yaml", "")
    c.save(os.path.join(OUT, f"{row['folder'][-7:]}__{name}.png"))
    print("wrote", name)
