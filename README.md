# SmogShield

Harmonizing NASA MODIS and VIIRS active-fire detections to support crop-burning smog monitoring in Punjab.

**NASA Space Apps Challenge 2026 - Challenge: Harmonization of MODIS and VIIRS Hot Spots**

## Problem
MODIS (1 km) and VIIRS (375 m) often detect the same fire, but with different resolution, overpass time and confidence scales. Counting both double-counts fires and makes smoke-source maps unreliable.

## What it does
1. Loads NASA FIRMS MODIS and VIIRS files (CSV or shapefile `.dbf`).
2. Converts both confidence scales to one common 0-1 score (MODIS 0-100 / 100; VIIRS low/nominal/high mapped to 0.25/0.60/0.90 - an assumption).
3. Links detections from different sensors within 1 km and 3 hours (configurable).
4. Merges each group into one fire event; when both sensors agree, confidence is raised using 1 - prod(1 - c).
5. Outputs `merged_hotspots.csv`, `stats.json` (raw vs merged counts) and an interactive `map.html`.

## Example result (recent FIRMS sample, Punjab region)
- Raw detections: 1,965 (408 MODIS, 1,557 VIIRS)
- Merged fire events: 1,946
- Duplicates removed: 19
- Events seen by both sensors: 13

## Run
```
pip install -r requirements.txt
python harmonize.py --modis modis.dbf --viirs viirs.dbf --out output --bbox 28 33 69 76
```
Open `output/map.html` in a browser (internet is needed for the background map).

## Data (NASA FIRMS)
Download MODIS and VIIRS active-fire files for South Asia from https://firms.modaps.eosdis.nasa.gov/active_fire/ (a free Earthdata login is required).

## Next steps
Combine the harmonized hotspots with NASA POWER wind data to estimate smoke travel and produce city-level smog risk warnings (not implemented yet).

## Data sources / references
- NASA FIRMS (MODIS C6.1, VIIRS 375 m active fire)
- Python, pandas, NumPy, SciPy, dbfread, Leaflet
- Map tiles: Esri; map data: OpenStreetMap contributors

## AI disclosure
AI tools (Claude) were used to help write parts of the code and documentation. The challenge choice, data download, testing and review were done by the team.

## License
Apache 2.0
