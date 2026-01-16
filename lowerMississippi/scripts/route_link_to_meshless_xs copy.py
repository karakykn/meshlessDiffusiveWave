# Modified route_link_to_meshless_xs.py that decodes the 'gages' char array more robustly,
# extracts USGS gage IDs (if present_baton), and writes gage URLs into xs_catalog.csv plus a separate gages.csv.
# It is saved as /mnt/data/route_link_to_meshless_xs_with_gages.py and is ready to run.
from pathlib import Path
code = r"""
#!/usr/bin/env python3
\"\"\"route_link_to_meshless_xs_with_gages.py

Reads NWM Route_Link.nc and writes meshless-model xs files (h, x_left, x_right).
Enhancement: robustly decodes the 'gages' field (often a char array of bytes)
and writes USGS gage URLs into the catalog output (xs_catalog.csv) and a separate gages.csv.

If the 'gages' entries are empty (e.g. b' '), the script leaves gage fields blank.
If a gage string contains non-digit characters, the script attempts to extract the station number
by taking the longest contiguous digit substring. This heuristic matches many NWM 'gages'
encodings where the COMID->USGS gage mapping is stored as a padded char field.

Output files:
  xs<N>                 (per-link cross section file: h x_left x_right)
  xs_catalog.csv        (catalog mapping xs_id -> feature_id -> geometry & gage info & urls)
  gages.csv             (compact mapping feature_id, xs_id -> gage_id, url)

Usage examples (same as original script):
  python route_link_to_meshless_xs_with_gages.py --routelink /path/Route_Link.nc \ 
        --start_gage 07374000 --end_gage 07374200 --out_dir ./xs --dh 0.25 --hmax 60

Dependencies:
  - netCDF4
\"\"\"

from __future__ import annotations
import argparse, os, math, re, csv
from pathlib import Path
from typing import List, Dict, Any, Optional

try:
    import netCDF4 as nc
except Exception as e:
    nc = None

def _require(cond: bool, msg: str):
    if not cond:
        raise RuntimeError(msg)

def _decode_char_array(var) -> List[str]:
    \"\"\"Convert various netCDF char/byte array formats to a list of python strings.
    Handles:
      - 2D byte arrays like shape (N, S) with type 'S1' (common)
      - 1D byte strings (already bytes objects)
      - missing/None -> returns empty list
    Returns one string per row; strips whitespace and nulls.
    \"\"\"
    if var is None:
        return []
    try:
        arr = var[:]
    except Exception:
        # fallback: treat as scalar
        try:
            v = var.getValue()
            return [str(v).strip()]
        except Exception:
            return []

    import numpy as np
    if isinstance(arr, np.ndarray):
        if arr.dtype.kind in ('S', 'U'):
            # If shape is (N,) of bytes/strings
            if arr.ndim == 1:
                out = []
                for v in arr.tolist():
                    if isinstance(v, (bytes, bytearray)):
                        s = v.decode('utf-8', errors='ignore').strip()
                    else:
                        s = str(v).strip()
                    out.append(s)
                return out
            # If shape is (N, strlen)
            if arr.ndim == 2:
                out = []
                for row in arr:
                    if row.dtype.kind == 'S':
                        # row is array of bytes (each element one char)
                        b = b''.join(row.tolist())
                        out.append(b.decode('utf-8', errors='ignore').strip().strip('\\x00').strip())
                    else:
                        # combine elements to string
                        try:
                            parts = [p.decode('utf-8', errors='ignore') if isinstance(p, (bytes, bytearray)) else str(p) for p in row.tolist()]
                            out.append(''.join(parts).strip())
                        except Exception:
                            out.append(str(row).strip())
                return out
        # Other dtypes: attempt converting each row to string
        if arr.ndim >= 2:
            out = []
            for row in arr.reshape(arr.shape[0], -1):
                try:
                    s = ''.join([chr(int(x)) if isinstance(x, (int,)) else str(x) for x in row.tolist()])
                except Exception:
                    s = ' '.join([str(x) for x in row.tolist()])
                out.append(s.strip())
            return out
        # 1D numeric array
        return [str(x).strip() for x in arr.tolist()]

    # scalar bytes or string
    if isinstance(arr, (bytes, bytearray)):
        return [arr.decode('utf-8', errors='ignore').strip()]
    return [str(arr).strip()]

def _extract_gage_id(s: str) -> Optional[str]:
    \"\"\"Heuristic: from a string, extract the most likely USGS station number (longest contiguous digit run).
    Returns None if no digits found.\"\"\"
    if not s:
        return None
    # remove obvious fillers
    s2 = s.strip().replace('\\x00','').replace('\\x20',' ').strip()
    # find all digit runs
    runs = re.findall(r\"(\\d{4,})\", s2)  # at least 4 digits (USGS ids are generally 8)
    if not runs:
        # fallback: any digit run
        runs = re.findall(r\"(\\d+)\", s2)
    if not runs:
        return None
    # pick the longest run (prefer 8+ digits if present_baton)
    runs_sorted = sorted(runs, key=lambda r: (-len(r), r))
    return runs_sorted[0]

def _usgs_monitoring_url(site_no: str) -> str:
    \"\"\"Return the USGS monitoring page URL for a site number (human-friendly).\"\"\"
    # canonical monitoring location page
    return f\"https://waterdata.usgs.gov/monitoring-location/{site_no}\"

def read_routelink(path: str) -> Dict[str, Any]:
    _require(nc is not None, \"netCDF4 is required. Please `pip install netCDF4`.\")
    _require(os.path.exists(path), f\"Route_Link file not found: {path}\")
    ds = nc.Dataset(path, mode='r')

    def safe_var(name):
        return ds.variables[name] if name in ds.variables else None

    # Typical variables; may vary by Route_Link version
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
    gages_var = safe_var('gages') if safe_var('gages') is not None else None
    gages = _decode_char_array(gages_var) if gages_var is not None else []

    ds.close()

    # normalize gages list length to number of links (if present_baton)
    if link is not None and len(gages) != len(link):
        # try to pad/truncate
        gages = (gages + [''] * len(link))[:len(link)]
    return dict(link=link, to=to, lat=lat, lon=lon,
                B=B, z=z, n=n_mann, nCC=n_cc, Tw=Tw, TwCC=Twcc, Length=length, gages=gages)

def chain_links_by_gages(rl: Dict[str, Any], start_gage: str, end_gage: Optional[str]) -> List[int]:
    gages = rl.get('gages', [])
    _require(len(gages) == len(rl['link']), \"gages length mismatch; this Route_Link may not include a 'gages' field.\")
    # normalize simple strings
    gnorm = [g.strip() if isinstance(g, str) else '' for g in gages]
    try:
        i_start = gnorm.index(start_gage.strip())
    except ValueError:
        raise RuntimeError(f\"start_gage {start_gage} not found in Route_Link 'gages'. Found: {','.join([g for g in gnorm if g])}\")
    i_end = None
    if end_gage:
        try:
            i_end = gnorm.index(end_gage.strip())
        except ValueError:
            i_end = None

    to = rl['to']
    _require(to is not None, \"Missing 'to' variable in Route_Link; cannot chain.\")
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
    lon = rl.get('lon', None); lat = rl.get('lat', None)
    _require(lon is not None and lat is not None, \"Route_Link lacks 'lon'/'lat'.\")
    idxs = [i for i in range(len(lon)) if (bbox[0] <= lon[i] <= bbox[2] and bbox[1] <= lat[i] <= bbox[3])]
    _require(len(idxs) > 0, \"No links found in bbox.\")
    selected = set(idxs)
    links = rl['link']; to = rl['to']
    id2idx = {int(links[i]): i for i in range(len(links))}
    pointed = set()
    for i in idxs:
        tid = int(to[i]) if to is not None else 0
        if tid in id2idx:
            j = id2idx[tid]
            if j in selected:
                pointed.add(j)
    starts = [i for i in idxs if i not in pointed]
    if not starts:
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
            raise RuntimeError(f\"feature_id {fid} not found in Route_Link.\")
        idxs.append(id2idx[int(fid)])
    return idxs

def compute_db_from_Tw_and_B(Tw: float, B: float, z: float) -> float:
    \"\"\"If TopWdth is given use db = (Tw - B) / (2*z) (handles z>0). Falls back to small depth if invalid.\"\"\"
    try:
        if z is None or float(z) <= 0:
            return max(0.01, (float(Tw) - float(B))/2.0)
        return max(0.01, (float(Tw) - float(B)) / (2.0 * float(z)))
    except Exception:
        return 0.01

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
            f.write(f\"{h:.6f} {xL:.6f} {xR:.6f}\\n\")


def main():
    ap = argparse.ArgumentParser(description=\"Convert NWM Route_Link compound trapezoid parameters to meshless xs* files (with gage URLs).\")
    ap.add_argument(\"--routelink\", required=True, help=\"Path to Route_Link.nc\")
    ap.add_argument(\"--out_dir\", required=True, help=\"Output directory for xs files\")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument(\"--start_gage\", help=\"USGS gage id to start chain (e.g., 07374000)\")
    ap.add_argument(\"--end_gage\", help=\"USGS gage id to stop chain (optional)\")
    group.add_argument(\"--feature_ids\", nargs=\"+\", type=int, help=\"Explicit list of feature_ids to include (upstream to downstream)\")
    group.add_argument(\"--bbox\", nargs=4, type=float, help=\"minlon minlat maxlon maxlat\")
    ap.add_argument(\"--dh\", type=float, default=0.2, help=\"Vertical sampling [m]\")
    ap.add_argument(\"--hmax\", type=float, default=None, help=\"Max depth to tabulate [m]; default = max(db+10, 2*db) per link\")
    ap.add_argument(\"--twcc_mult\", type=float, default=3.0, help=\"If TopWdthCC missing, use multiplier * TopWdth\")
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
    gages_list = rl.get('gages', [])

    catalog_rows = [[\"xs_id\", \"feature_id\", \"BtmWdth\", \"ChSlp\", \"TopWdth\", \"TopWdthCC\", \"dbankfull\", \"n\", \"nCC\", \"gage_raw\", \"gage_id\", \"gage_url\"]]
    gage_map_rows = [[\"xs_id\",\"feature_id\",\"gage_raw\",\"gage_id\",\"gage_url\"]]
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
        dbi = compute_db_from_Tw_and_B(Twi, Bi, zi)
        xs_path = os.path.join(args.out_dir, f\"xs{xs_id}\")
        write_xs_file(xs_path, Bi, zi, dbi, Twcci, hmax=args.hmax, dh=args.dh)

        # gage parsing
        raw_g = ''
        gid = None
        url = ''
        try:
            raw_g = gages_list[idx] if idx < len(gages_list) else ''
            if isinstance(raw_g, bytes):
                try:
                    raw_g = raw_g.decode('utf-8', errors='ignore')
                except Exception:
                    raw_g = str(raw_g)
            raw_g = str(raw_g).strip()
            gid = _extract_gage_id(raw_g)
            if gid:
                url = _usgs_monitoring_url(gid)
        except Exception:
            raw_g = ''
            gid = None
            url = ''

        catalog_rows.append([xs_id, fid, Bi, zi, Twi, Twcci, dbi, (float(n[idx]) if n is not None else ''), (float(nCC[idx]) if nCC is not None else ''), raw_g, (gid if gid else ''), url])
        gage_map_rows.append([xs_id, fid, raw_g, (gid if gid else ''), url])
        xs_id += 1

    catalog_path = os.path.join(args.out_dir, \"xs_catalog.csv\")
    with open(catalog_path, 'w', newline='') as cf:
        cw = csv.writer(cf)
        cw.writerows(catalog_rows)

    gages_path = os.path.join(args.out_dir, \"gages.csv\")
    with open(gages_path, 'w', newline='') as gf:
        gw = csv.writer(gf)
        gw.writerows(gage_map_rows)

    print(f\"Wrote {xs_id} cross-section files to {args.out_dir}\")
    print(\"Catalog:\", catalog_path)
    print(\"Gage map:\", gages_path)
    print(\"Tip: map your model nodes to xs_id via xsInfo using the Length of each link if you want per-link geometry.\")


if __name__ == '__main__':
    main()
"""

out_path = Path("/mnt/data/route_link_to_meshless_xs_with_gages.py")
out_path.write_text(code)
print("Wrote file:", out_path)
print("\nQuick usage example:")
print("python /mnt/data/route_link_to_meshless_xs_with_gages.py --routelink /path/to/Route_Link.nc --start_gage 07374000 --end_gage 07374200 --out_dir ./xs --dh 0.25 --hmax 60")
