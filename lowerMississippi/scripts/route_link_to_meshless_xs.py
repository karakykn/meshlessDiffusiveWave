# Create a small utility that reads an NWM Route_Link.nc, extracts compound trapezoid
# parameters for a selected chain of links (e.g., Baton Rouge -> Donaldsonville),
# and writes meshless cross-section files `xs*` (format: h, x_left, x_right).
#
# NOTE: This code assumes you will run it on your own machine where the Route_Link.nc
# exists. Internet access is not required.
#
# It also writes a CSV catalog so you can map from NWM feature_id to the xs file name.
#
# Usage examples (after saving this file as route_link_to_meshless_xs.py):
#
# 1) By gage IDs (recommended):
#    python route_link_to_meshless_xs.py \
#       --routelink /path/to/Route_Link.nc \
#       --start_gage 07374000 --end_gage 07374200 \
#       --out_dir ./xs --dh 0.25 --hmax 60
#
# 2) By a list of feature_ids (if you already know the links to use, upstream to downstream):
#    python route_link_to_meshless_xs.py \
#       --routelink /path/to/Route_Link.nc \
#       --feature_ids 12345678 12345679 12345680 \
#       --out_dir ./xs
#
# 3) By a lat/lon bounding box (rough cut; the script will attempt to order by the 'to' chain):
#    python route_link_to_meshless_xs.py \
#       --routelink /path/to/Route_Link.nc \
#       --bbox -91.25 30.05 -90.85 30.50 \
#       --out_dir ./xs
#
# The script will:
#  - Read BtmWdth (B), ChSlp (z), TopWdth (Tw), TopWdthCC (Twcc if available), Manning's n (n) and nCC (if available)
#  - Compute bankfull depth db from Tw = B + 2 z db  =>  db = (Tw - B)/(2 z)  (if z>0)
#  - Build cross-section files xs<k> with columns: h, x_left(h), x_right(h)
#    using trapezoid up to db and vertical walls to Twcc above db
#  - Write xs_catalog.csv listing: xs_id, feature_id, B, z, Tw, Twcc, db, n, nCC
#
# You can then set xsInfo in your model to reference these xs files along your segment.
#
import os
import argparse
import math
from typing import List, Dict, Any, Optional

try:
    import netCDF4 as nc
except Exception as e:
    nc = None

import csv
from pathlib import Path

SCRIPT_PATH = Path("/mnt/data/route_link_to_meshless_xs.py")
README_PATH = Path("/mnt/data/route_link_to_meshless_README.txt")

# script_code =
import os
import argparse
import math
from typing import List, Dict, Any, Optional
from pathlib import Path
import csv

try:
    import netCDF4 as nc
except Exception as e:
    nc = None

def _require(cond: bool, msg: str):
    if not cond:
        raise RuntimeError(msg)

def _decode_char_array(var) -> list:
    """
    Convert a netCDF4 char array (N x strlen) to a list of Python strings.
    If var is missing, return empty list.
    """
    if var is None:
        return []
    # netCDF4 returns numpy array of type 'S1' for char variables.
    # We join along last axis.
    import numpy as np
    arr = var[:]
    if arr.dtype.kind in ('S', 'U'):
        # Sometimes it's already bytes strings per-element
        # Ensure 2D (N, strlen)
        a = arr
        if a.ndim == 1:
            return [str(x, 'utf-8') if isinstance(x, (bytes, bytearray)) else str(x) for x in a.tolist()]
        elif a.ndim == 2:
            out = []
            for row in a:
                if row.dtype.kind in ('S', 'U'):
                    s = b''.join(row).decode('utf-8', errors='ignore') if row.dtype.kind == 'S' else ''.join(row.tolist())
                else:
                    s = ''.join([str(x) for x in row.tolist()])
                out.append(s.strip())
            return out
    # Fallback: assume last dim is strlen
    if arr.ndim >= 2:
        # join along last axis
        joined = [''.join([chr(c) if isinstance(c, (int,)) else (c.decode("utf-8") if isinstance(c, (bytes, bytearray)) else str(c))
                           for c in row]).strip()
                  for row in arr.reshape(arr.shape[0], -1)]
        return joined
    # Single string?
    val = arr.tolist()
    if isinstance(val, (bytes, bytearray)):
        return [val.decode('utf-8')]
    return [str(val)]

def read_routelink(path: str) -> Dict[str, Any]:
    _require(nc is not None, "netCDF4 is required. Please `pip install netCDF4`.")
    _require(os.path.exists(path), f"Route_Link file not found: {path}")
    ds = nc.Dataset(path, mode='r')

    def safe_var(name):
        return ds.variables[name] if name in ds.variables else None

    link = safe_var('link')[:] if safe_var('link') is not None else None
    to = safe_var('to')[:] if safe_var('to') is not None else None
    lat = safe_var('lat')[:] if safe_var('lat') is not None else None
    lon = safe_var('lon')[:] if safe_var('lon') is not None else None
    B = safe_var('BtmWdth')[:] if safe_var('BtmWdth') is not None else None
    z = safe_var('ChSlp')[:] if safe_var('ChSlp') is not None else None
    n_mann = safe_var('n')[:] if safe_var('n') is not None else None
    n_cc = safe_var('nCC')[:] if safe_var('nCC') is not None else None
    Tw = safe_var('TopWdth')[:] if safe_var('TopWdth') is not None else None
    Twcc = safe_var('TopWdthCC')[:] if safe_var('TopWdthCC') is not None else None
    length = safe_var('Length')[:] if safe_var('Length') is not None else None
    gages = _decode_char_array(safe_var('gages')) if safe_var('gages') is not None else []

    ds.close()

    return dict(link=link, to=to, lat=lat, lon=lon,
                B=B, z=z, n=n_mann, nCC=n_cc, Tw=Tw, TwCC=Twcc, Length=length, gages=gages)

def chain_links_by_gages(rl: Dict[str, Any], start_gage: str, end_gage: Optional[str]) -> List[int]:
    gages = rl['gages']
    _require(len(gages) == len(rl['link']), "gages length mismatch; this Route_Link may not include a 'gages' field.")
    # Normalize gage strings to just digits (strip spaces)
    gnorm = [g.strip() for g in gages]
    try:
        i_start = gnorm.index(start_gage.strip())
    except ValueError:
        raise RuntimeError(f"start_gage {start_gage} not found in Route_Link 'gages'.")
    i_end = None
    if end_gage:
        try:
            i_end = gnorm.index(end_gage.strip())
        except ValueError:
            i_end = None

    # Follow 'to' pointers from i_start until i_end or termination.
    to = rl['to']
    _require(to is not None, "Missing 'to' variable in Route_Link; cannot chain.")
    links = rl['link']
    id2idx = {int(links[i]): i for i in range(len(links))}

    chain = [i_start]
    cur = i_start
    visited = set([cur])
    while True:
        next_link_id = int(to[cur])
        if next_link_id == 0:
            break
        if next_link_id not in id2idx:
            break
        nxt = id2idx[next_link_id]
        if nxt in visited:
            break
        chain.append(nxt)
        visited.add(nxt)
        if i_end is not None and nxt == i_end:
            break
        cur = nxt
    return chain

def chain_links_by_bbox(rl: Dict[str, Any], bbox: List[float]) -> List[int]:
    """
    bbox = [minlon, minlat, maxlon, maxlat]; rough subset and then order along 'to' chain.
    We will select all links whose midpoints fall within bbox, then find a start (not pointed to by any selected) and follow 'to' within subset.
    """
    lon = rl['lon']; lat = rl['lat']
    _require(lon is not None and lat is not None, "Route_Link lacks 'lon'/'lat'.")
    idxs = [i for i in range(len(lon)) if (bbox[0] <= lon[i] <= bbox[2] and bbox[1] <= lat[i] <= bbox[3])]
    _require(len(idxs) > 0, "No links found in bbox.")
    selected = set(idxs)
    links = rl['link']; to = rl['to']
    id2idx = {int(links[i]): i for i in range(len(links))}
    # find potential starts: nodes in selected not pointed to by any other selected
    pointed = set()
    for i in idxs:
        tid = int(to[i]) if to is not None else 0
        if tid in id2idx:
            j = id2idx[tid]
            if j in selected:
                pointed.add(j)
    starts = [i for i in idxs if i not in pointed]
    if not starts:
        # fallback: pick the northernmost (max lat) as upstream
        starts = sorted(idxs, key=lambda i: (-lat[i], lon[i]))
    chain = []
    visited = set()
    cur = starts[0]
    while cur not in visited and cur in selected:
        chain.append(cur)
        visited.add(cur)
        nxt_id = int(to[cur]) if to is not None else 0
        if nxt_id == 0 or nxt_id not in id2idx:
            break
        nxt = id2idx[nxt_id]
        if nxt not in selected:
            break
        cur = nxt
    return chain

def chain_links_by_ids(rl: Dict[str, Any], feature_ids: List[int]) -> List[int]:
    links = rl['link']
    id2idx = {int(links[i]): i for i in range(len(links))}
    idxs = []
    for fid in feature_ids:
        if int(fid) not in id2idx:
            raise RuntimeError(f"feature_id {fid} not found in Route_Link.")
        idxs.append(id2idx[int(fid)])
    return idxs

def compute_db(Tw: float, B: float, z: float) -> float:
    # Bankfull depth from trapezoid geometry: Tw = B + 2*z*db
    if z is None or z <= 0:
        return max(0.01, (Tw - B)/2.0)  # fallback
    return max(0.01, (Tw - B) / (2.0 * z))

def write_xs_file(path: str, B: float, z: float, db: float, Twcc: float,
                  hmax: float = None, dh: float = 0.1):
    if hmax is None:
        hmax = max(db + 10.0, 2.0*db)
    nsteps = int(math.floor(hmax / dh)) + 1
    with open(path, 'w') as f:
        for k in range(nsteps):
            h = k * dh
            if h <= db:
                xL = -B/2.0 - z*h
                xR = +B/2.0 + z*h
            else:
                xL = -Twcc/2.0
                xR = +Twcc/2.0
            f.write(f"{h:.6f} {xL:.6f} {xR:.6f}\n")

def main():
    ap = argparse.ArgumentParser(description="Convert NWM Route_Link compound trapezoid parameters to meshless xs* files.")
    ap.add_argument("--routelink", required=True, help="Path to Route_Link.nc")
    ap.add_argument("--out_dir", required=True, help="Output directory for xs files")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--start_gage", help="USGS gage id to start chain (e.g., 07374000)")
    ap.add_argument("--end_gage", help="USGS gage id to stop chain (optional)")
    group.add_argument("--feature_ids", nargs="+", type=int, help="Explicit list of feature_ids to include (upstream to downstream)")
    group.add_argument("--bbox", nargs=4, type=float, help="minlon minlat maxlon maxlat")
    ap.add_argument("--dh", type=float, default=0.2, help="Vertical sampling [m]")
    ap.add_argument("--hmax", type=float, default=None, help="Max depth to tabulate [m]; default = max(db+10, 2*db) per link")
    ap.add_argument("--twcc_mult", type=float, default=3.0, help="If TopWdthCC missing, use multiplier * TopWdth")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rl = read_routelink(args.routelink)

    if args.feature_ids:
        chain = chain_links_by_ids(rl, args.feature_ids)
    elif args.start_gage:
        chain = chain_links_by_gages(rl, args.start_gage, args.end_gage)
    else:
        chain = chain_links_by_bbox(rl, args.bbox)

    links = rl['link']
    B = rl['B']; z = rl['z']; Tw = rl['Tw']; TwCC = rl['TwCC']
    n = rl['n']; nCC = rl['nCC']

    catalog_rows = [["xs_id", "feature_id", "BtmWdth", "ChSlp", "TopWdth", "TopWdthCC", "dbankfull", "n", "nCC"]]
    xs_id = 0
    for idx in chain:
        fid = int(links[idx])
        Bi = float(B[idx]) if B is not None else 0.0
        zi = max(1e-6, float(z[idx])) if z is not None else 1.0
        Twi = float(Tw[idx]) if Tw is not None else (Bi + 2*zi*5.0)
        if TwCC is not None:
            Twcci = float(TwCC[idx])
        else:
            Twcci = args.twcc_mult * Twi
        dbi = compute_db(Twi, Bi, zi)
        xs_path = os.path.join(args.out_dir, f"xs{xs_id}")
        write_xs_file(xs_path, Bi, zi, dbi, Twcci, hmax=args.hmax, dh=args.dh)
        catalog_rows.append([xs_id, fid, Bi, zi, Twi, Twcci, dbi, (float(n[idx]) if n is not None else ""), (float(nCC[idx]) if nCC is not None else "")])
        xs_id += 1

    with open(os.path.join(args.out_dir, "xs_catalog.csv"), "w", newline="") as cf:
        cw = csv.writer(cf)
        cw.writerows(catalog_rows)

    print(f"Wrote {xs_id} cross-section files to {args.out_dir}")
    print("Catalog:", os.path.join(args.out_dir, "xs_catalog.csv"))
    print("Tip: map your model nodes to xs_id via xsInfo using the Length of each link if you want per-link geometry.")

if __name__ == "__main__":
    main()


readme_text = r'''
route_link_to_meshless_xs.py
============================

What it does
------------
- Reads **NWM/WRF‑Hydro Route_Link.nc** (the reach parameter file).
- Pulls the per‑reach compound trapezoid parameters used by the National Water Model:
  * Bottom width (`BtmWdth`), side slope (`ChSlp`),
  * Bankfull top width (`TopWdth`), compound top width (`TopWdthCC`),
  * Roughness (`n`) and compound roughness (`nCC`, if present_baton).
- Computes **bankfull depth** as `db = (TopWdth - BtmWdth) / (2 * ChSlp)`.
- Writes one **meshless cross‑section file** (`xs*`) per selected reach with columns:
  `h  x_left(h)  x_right(h)`, using trapezoid up to `db` and vertical overbank walls to `TopWdthCC`.
- Also writes `xs_catalog.csv` so you can map NWM `feature_id` → your `xs*` files.

How to run
----------
Example A — by gage IDs (Baton Rouge → Donaldsonville):
'''