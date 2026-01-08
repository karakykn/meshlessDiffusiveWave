import netCDF4 as nc
from typing import List, Dict, Any, Optional
import argparse
import os
from pathlib import Path
import csv
import numpy as np
import pandas as pd


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

def _require(cond: bool, msg: str):
    if not cond:
        raise RuntimeError(msg)

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

def upstream_ids_numpy(route, start_ids, depth=10):
    """
    route: dict with keys 'to' and 'link' (iterables/arrays of same length)
           meaning: link_arr[i] -> to_arr[i]  (i.e., link flows into to)
    start_ids: scalar or iterable of starting ID(s) to find upstream from
    depth: how many upstream steps to follow
    Returns: list of numpy arrays: result[0] = start_ids, result[1] = upstream at 1 step, ...
    """
    def upstream_ids_numpy(route, start_ids, depth):
        to_arr = np.array(route['to'])
        link_arr = np.array(route['link'])
        current = np.unique(np.atleast_1d(start_ids))
        results = [current.copy()]
        for level in range(depth):
            mask = np.isin(to_arr, current)
            found = link_arr[mask]
            if found.size == 0:
                break
            found = np.unique(found)
            results.append(found)
            current = found
        return results

    levels = upstream_ids_numpy(route, start_ids, depth)

    df = pd.DataFrame(route)

    # 3️⃣ Collect matching rows for each level
    rows = []
    for lvl, ids in enumerate(levels):
        subset = df[df['link'].isin(ids)].copy()
        subset['level'] = lvl
        rows.append(subset)

    filename = f'../data/batonUpstreams.csv'

    if rows:
        all_df = pd.concat(rows, ignore_index=True)
        all_df.to_csv(filename, index=False)
        print(f"Saved upstream info for {len(all_df)} links to '{filename}'")
        return all_df
    else:
        print("No upstream links found.")
        return pd.DataFrame()

SCRIPT_PATH = Path("/mnt/data/route_link_to_meshless_xs.py")
README_PATH = Path("/mnt/data/route_link_to_meshless_README.txt")
route_path = f'../../../Research/Diffusive Wave RBFCM/routelink/RouteLink_CONUS.nc'
route = read_routelink(route_path)

ids = 19088319 #baton

levels = upstream_ids_numpy(route, ids, depth=50)
# prints each level
for lvl, ids in enumerate(levels):
    print(f"level {lvl}: {ids}")