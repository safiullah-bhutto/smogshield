#!/usr/bin/env python3
"""
SmogShield - harmonize NASA FIRMS MODIS + VIIRS active-fire detections.

MODIS (1 km) and VIIRS (375 m) often detect the same fire, with different
overpass times and different confidence scales. This script:
  1. loads MODIS and VIIRS FIRMS CSV files
  2. puts confidence on one 0-1 scale
  3. links detections from DIFFERENT sensors that are close in space and time
  4. merges each linked group into one fire event
  5. writes merged_hotspots.csv, stats.json and map.html

Usage:
  python harmonize.py --modis modis.csv --viirs viirs.csv --out output
  (optional) --bbox 23.5 29.5 69.0 75.5   # south north west east (Punjab-ish)

Data: NASA FIRMS (https://firms.modaps.eosdis.nasa.gov/)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# VIIRS confidence is categorical (l/n/h). These numeric values are an
# ASSUMPTION used to put both sensors on one 0-1 scale.
VIIRS_CONF = {"l": 0.25, "low": 0.25, "n": 0.60, "nominal": 0.60,
              "h": 0.90, "high": 0.90}


def load(path):
    if str(path).lower().endswith(".dbf"):
        # Shapefile attribute table (FIRMS shapefile download): needs `pip install dbfread`
        from dbfread import DBF
        df = pd.DataFrame(list(DBF(str(path), load=True, encoding="utf-8")))
    else:
        df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    is_viirs = "bright_ti4" in df.columns
    out = pd.DataFrame({
        "lat": df["latitude"].astype(float),
        "lon": df["longitude"].astype(float),
        "date": df["acq_date"].astype(str),
        "frp": pd.to_numeric(df["frp"], errors="coerce") if "frp" in df else np.nan,
        "source": "VIIRS" if is_viirs else "MODIS",
    })
    hhmm = df["acq_time"].astype(int).astype(str).str.zfill(4)
    out["minutes"] = hhmm.str[:2].astype(int) * 60 + hhmm.str[2:].astype(int)
    if is_viirs:
        out["conf"] = df["confidence"].astype(str).str.strip().str.lower().map(VIIRS_CONF)
    else:
        out["conf"] = pd.to_numeric(df["confidence"], errors="coerce") / 100.0
    out["conf"] = out["conf"].fillna(0.5)
    return out


def harmonize(df, radius_km=1.0, window_min=180):
    df = df.reset_index(drop=True)
    parent = np.arange(len(df))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    lat0 = np.deg2rad(df["lat"].mean())
    x = df["lon"].values * 111.32 * np.cos(lat0)
    y = df["lat"].values * 110.57
    minutes = df["minutes"].values
    source = df["source"].values

    for _, idx in df.groupby("date").indices.items():
        if len(idx) < 2:
            continue
        tree = cKDTree(np.c_[x[idx], y[idx]])
        for a, b in tree.query_pairs(radius_km):
            ia, ib = idx[a], idx[b]
            # only link detections from different sensors (cross-sensor duplicates)
            if source[ia] != source[ib] and abs(minutes[ia] - minutes[ib]) <= window_min:
                ra, rb = find(ia), find(ib)
                if ra != rb:
                    parent[rb] = ra

    df["event_id"] = [find(i) for i in range(len(df))]
    g = df.groupby("event_id")
    ev = g.agg(lat=("lat", "mean"), lon=("lon", "mean"), date=("date", "first"),
               frp=("frp", "max"), n_detections=("lat", "size"),
               sources=("source", lambda s: "+".join(sorted(set(s))))).reset_index(drop=True)
    # agreement between sensors raises confidence: 1 - prod(1 - c_sensor)
    per_sensor = df.groupby(["event_id", "source"])["conf"].max()
    combined = per_sensor.groupby("event_id").apply(lambda c: 1 - np.prod(1 - c.values))
    ev["confidence"] = combined.reindex(g.size().index).values.round(3)
    ev["lat"] = ev["lat"].round(5)
    ev["lon"] = ev["lon"].round(5)
    return ev


MAP_HTML = """<!doctype html><html><head><meta charset="utf-8">
<title>SmogShield - Harmonized Fire Hotspots</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>html,body,#map{height:100%;margin:0}
.box{position:absolute;z-index:999;top:10px;right:10px;background:#fff;padding:10px 14px;
font:14px sans-serif;border-radius:8px;box-shadow:0 1px 6px #0004;max-width:260px}
.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px}</style></head>
<body><div id="map"></div>
<div class="box"><b>SmogShield</b><br>Harmonized MODIS + VIIRS fires<br><br>
<div id="stats"></div><br>
<span class="dot" style="background:#d62728"></span>Seen by both sensors<br>
<span class="dot" style="background:#ff7f0e"></span>VIIRS only<br>
<span class="dot" style="background:#1f77b4"></span>MODIS only<br>
<small>Data: NASA FIRMS</small></div>
<script>
const events = __EVENTS__; const stats = __STATS__;
const map = L.map('map').setView([__LAT__, __LON__], 7);
const esriStreet = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}',
  {attribution:'Tiles &copy; Esri',maxZoom:18}).addTo(map);
const esriSat = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  {attribution:'Tiles &copy; Esri',maxZoom:18});
L.control.layers({'Street map':esriStreet,'Satellite':esriSat}).addTo(map);
map.fitBounds(L.latLngBounds(events.map(e=>[e.lat,e.lon])), {padding:[30,30]});
const color = s => s.includes('+') ? '#d62728' : (s==='VIIRS' ? '#ff7f0e' : '#1f77b4');
events.forEach(e => L.circleMarker([e.lat,e.lon],
  {radius:4+Math.min(8,(e.frp||0)/15),color:color(e.sources),fillOpacity:.6,weight:1})
  .bindPopup(`<b>${e.sources}</b><br>Date: ${e.date}<br>Confidence: ${e.confidence}<br>FRP: ${e.frp}`).addTo(map));
document.getElementById('stats').innerHTML =
  `Raw detections: <b>${stats.raw_total}</b><br>Merged events: <b>${stats.events}</b><br>`+
  `Duplicates removed: <b>${stats.duplicates_removed}</b><br>Seen by both: <b>${stats.seen_by_both}</b>`;
</script></body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modis", nargs="+", required=True)
    ap.add_argument("--viirs", nargs="+", required=True)
    ap.add_argument("--out", default="output")
    ap.add_argument("--bbox", nargs=4, type=float, metavar=("S", "N", "W", "E"),
                    help="optional filter: south north west east")
    ap.add_argument("--radius-km", type=float, default=1.0)
    ap.add_argument("--window-min", type=int, default=180)
    a = ap.parse_args()

    raw = pd.concat([load(p) for p in a.modis + a.viirs], ignore_index=True)
    if a.bbox:
        s, n, w, e = a.bbox
        raw = raw[(raw.lat >= s) & (raw.lat <= n) & (raw.lon >= w) & (raw.lon <= e)]
    if raw.empty:
        raise SystemExit("No detections after filtering. Check your files / bbox.")

    ev = harmonize(raw, a.radius_km, a.window_min)
    stats = {
        "raw_modis": int((raw.source == "MODIS").sum()),
        "raw_viirs": int((raw.source == "VIIRS").sum()),
        "raw_total": int(len(raw)),
        "events": int(len(ev)),
        "duplicates_removed": int(len(raw) - len(ev)),
        "seen_by_both": int(ev["sources"].str.contains(r"\+").sum()),
        "radius_km": a.radius_km,
        "window_min": a.window_min,
    }

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ev.to_csv(out / "merged_hotspots.csv", index=False)
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    html = (MAP_HTML.replace("__EVENTS__", ev.to_json(orient="records"))
            .replace("__STATS__", json.dumps(stats))
            .replace("__LAT__", str(round(ev.lat.mean(), 3)))
            .replace("__LON__", str(round(ev.lon.mean(), 3))))
    (out / "map.html").write_text(html)
    print(json.dumps(stats, indent=2))
    print(f"Wrote {out}/merged_hotspots.csv, stats.json, map.html")


if __name__ == "__main__":
    main()