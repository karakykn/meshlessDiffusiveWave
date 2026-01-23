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
    alt = safe_var('alt')[:] if safe_var('alt') is not None else None
    S0 = safe_var('So')[:] if safe_var('So') is not None else None
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
                B=B, z=z, n=n_mann, nCC=n_cc, Tw=Tw, TwCC=Twcc, Length=length, alt = alt, So = S0, gages=gages)

import pandas as pd

def downstream_chain_to_csv(route, start_link, filename='downstream_chain.csv', save=True):
    """
    Follow downstream links from start_link using route['to'] mapping,
    collect all rows for each visited link (including the start), and save to CSV.

    Parameters
    ----------
    route : dict-like
        Must contain at least 'link' and 'to' (iterables of same length).
        Other columns in route will be preserved in the output csv.
    start_link : scalar
        The initial link id to start from.
    filename : str
        CSV filename to write (if save=True).
    save : bool
        If True, writes the CSV to filename. If False, only returns the dataframe.

    Returns
    -------
    df_chain : pandas.DataFrame
        DataFrame containing the collected rows in downstream order with extra columns:
         - step : integer step from start (0 -> start_link)
         - is_initial : True for the start link row(s)
         - cycle_detected : True if we stopped because of a cycle (set only for final row)
    """
    # make DataFrame from route dict
    df = pd.DataFrame(route)

    # quick checks
    if 'link' not in df.columns or 'to' not in df.columns:
        raise ValueError("route must contain at least 'link' and 'to' keys/columns")

    # build mapping link -> to (if multiple rows per link, keep first mapping but warn later)
    link_to_map = {}
    duplicates = df[df.duplicated(subset=['link'], keep=False)]
    if not duplicates.empty:
        # there are duplicate link entries; we will still produce rows for each,
        # but downstream stepping uses the first 'to' encountered per link id.
        dup_links = sorted(duplicates['link'].unique())
        print(f"Warning: duplicate rows found for link ids: {dup_links}. "
              "All rows will be included in output, but downstream following will use the first 'to' value for each link.")

    # use first occurrence as mapping
    for _, row in df.drop_duplicates(subset=['link']).iterrows():
        link_to_map[row['link']] = row['to']

    chain_rows = []
    visited = set()
    current = start_link
    step = 0
    cycle_detected = False

    while True:
        # If already visited -> cycle
        if current in visited:
            cycle_detected = True
            print(f"Cycle detected at link id {current}. Stopping traversal.")
            break

        visited.add(current)

        # collect all rows that belong to this link (could be multiple)
        rows_this_link = df[df['link'] == current].copy()
        if rows_this_link.empty:
            # no row for this link id in route: still record a minimal row
            # construct a minimal dict so CSV indicates missing attributes
            minimal = {'link': current, 'to': link_to_map.get(current, None)}
            rows_this_link = pd.DataFrame([minimal])
            print(f"Note: no row found in route for link id {current}; placing a minimal record.")

        # add bookkeeping columns
        rows_this_link['step'] = step
        rows_this_link['is_initial'] = (step == 0)
        rows_this_link['cycle_detected'] = False  # may be set later if final cause is cycle

        chain_rows.append(rows_this_link)
        step += 1

        # find downstream id to continue
        # prefer the mapping built earlier (first occurrence); fallback to the 'to' value from the collected row(s)
        downstream = link_to_map.get(current, None)
        if downstream is None:
            # try to get 'to' from the first collected row
            downstream = rows_this_link.iloc[0].get('to', None)

        # stop when downstream is missing or equals 0
        if pd.isna(downstream) or downstream == 0:
            break

        # continue to downstream link id
        current = downstream

    # mark cycle flag on last row if cycle detected
    if cycle_detected and chain_rows:
        chain_rows[-1]['cycle_detected'] = True

    # concat rows
    df_chain = pd.concat(chain_rows, ignore_index=True, sort=False)

    # reorder columns: put bookkeeping columns at front
    cols = list(df_chain.columns)
    for special in ['step', 'is_initial', 'cycle_detected']:
        if special in cols:
            cols.insert(0, cols.pop(cols.index(special)))
    df_chain = df_chain[cols]

    if save:
        df_chain.to_csv(filename, index=False)
        print(f"Saved downstream chain (length {len(df_chain['step'].unique())} steps, {len(df_chain)} rows) -> '{filename}'")

    return df_chain

route_path = f'../../../Research/Diffusive Wave RBFCM/routelink/RouteLink_CONUS.nc'
route = read_routelink(route_path)
id = 15024889
df_out = downstream_chain_to_csv(route, id, filename='../data/from_sacramento.csv')