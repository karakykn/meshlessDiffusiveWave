#!/usr/bin/env python3
"""
inspect_routelink_batonrouge.py

Usage:
    python inspect_routelink_batonrouge.py /path/to/RouteLink_CONUS.nc

What it does:
- Opens the given RouteLink (NetCDF) file
- Decodes the 'gages' variable robustly
- Finds entries where the gage == "07374000"
- Prints a readable table of key parameters for those entries
"""

import sys
from pathlib import Path
import numpy as np

try:
    from netCDF4 import Dataset
except Exception as e:
    raise SystemExit("This script requires the netCDF4 package. Install with: pip install netCDF4")

def decode_gages_var(g_raw):
    """Robustly convert a netCDF 'gages' variable slice to list[str]."""
    if g_raw is None:
        return []
    # If it's already a 1-D array of bytes/strings
    if getattr(g_raw, "ndim", 0) == 1:
        out = []
        for elem in g_raw:
            if isinstance(elem, (bytes, bytearray, np.bytes_)):
                out.append(elem.decode('utf-8', errors='ignore').strip())
            else:
                # sometimes netCDF returns numpy.str_ or numbers
                out.append(str(elem).strip())
        return out
    # If it's a 2-D char array (N x strlen)
    if getattr(g_raw, "ndim", 0) >= 2:
        out = []
        # join bytes across last axis per record
        for row in g_raw:
            # dtype 'S' (bytes)
            if row.dtype.kind == 'S':
                try:
                    joined = b"".join(row.tolist())
                    out.append(joined.decode('utf-8', errors='ignore').strip())
                except Exception:
                    out.append("".join([ (x.decode('utf-8', errors='ignore') if isinstance(x, (bytes, np.bytes_)) else str(x)) for x in row]).strip())
            elif row.dtype.kind == 'U':
                out.append("".join(row.tolist()).strip())
            else:
                # numeric codes or fallback
                try:
                    out.append("".join([chr(int(x)) for x in row.flatten().tolist()]).strip())
                except Exception:
                    out.append("".join([str(x) for x in row.flatten().tolist()]).strip())
        return out
    # final fallback
    try:
        return [ (e.decode('utf-8') if isinstance(e, (bytes, np.bytes_)) else str(e)).strip() for e in list(g_raw) ]
    except Exception:
        return [str(g_raw).strip()]

def get_var(ds, name):
    """Return variable array if exists, else None."""
    return ds.variables[name][:] if name in ds.variables else None

def maybe_get_value(arr, idx):
    """Return arr[idx] as Python scalar or '' if missing."""
    if arr is None:
        return ""
    try:
        val = arr[idx]
        # numpy types -> native
        if isinstance(val, (np.generic,)):
            return val.item()
        return val
    except Exception:
        return ""

def pretty_print_rows(rows, headers):
    # compute column widths
    widths = [max(len(str(h)), max(len(str(r[i])) for r in rows)) for i,h in enumerate(headers)]
    # header
    hdr = "  ".join(h.ljust(widths[i]) for i,h in enumerate(headers))
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print("  ".join(str(r[i]).ljust(widths[i]) for i in range(len(headers))))

def main(ncpath):
    ncpath = Path(ncpath)
    if not ncpath.exists():
        raise SystemExit(f"File not found: {ncpath}")

    ds = Dataset(str(ncpath), "r")
    try:
        # Load variables of interest (if present_baton)
        vars_of_interest = ["link", "to", "lat", "lon", "Length", "So", "Slope",
                            "BtmWdth", "ChSlp", "TopWdth", "TopWdthCC", "n", "nCC", "gages"]
        # Read raw arrays if available
        d = {}
        for v in vars_of_interest:
            d[v] = get_var(ds, v)

        # decode gages robustly
        g_raw = d.get("gages")
        gages = decode_gages_var(g_raw) if g_raw is not None else []

        # Find Baton Rouge gage entries
        target = "07374000"
        matches = [i for i,g in enumerate(gages) if g == target]

        if not matches:
            print(f"Gage {target} not found in 'gages' variable. Trying substring search...")
            matches = [i for i,g in enumerate(gages) if target in g]
        if not matches:
            print(f"No matches for {target}. You can inspect a slice of gages like gages[:50].")
            # Print a small sample so user can see format
            print("Example gages sample (first 50):")
            print(gages[:50])
            return

        # For each matching index, collect readable records
        rows = []
        headers = ["idx", "feature_id", "to", "lat", "lon", "Length_m", "So",
                   "BtmWdth_m", "ChSlp", "TopWdth_m", "TopWdthCC_m", "n", "nCC", "gage_str"]
        for idx in matches:
            link = maybe_get_value(d.get("link"), idx)
            tov  = maybe_get_value(d.get("to"), idx)
            lat  = maybe_get_value(d.get("lat"), idx)
            lon  = maybe_get_value(d.get("lon"), idx)
            length = maybe_get_value(d.get("Length"), idx)
            so = maybe_get_value(d.get("So"), idx) or maybe_get_value(d.get("Slope"), idx)
            Btm = maybe_get_value(d.get("BtmWdth"), idx)
            ChSlp = maybe_get_value(d.get("ChSlp"), idx)
            Tw = maybe_get_value(d.get("TopWdth"), idx)
            Twcc = maybe_get_value(d.get("TopWdthCC"), idx)
            n = maybe_get_value(d.get("n"), idx)
            ncc = maybe_get_value(d.get("nCC"), idx)
            gstr = gages[idx] if idx < len(gages) else ""
            rows.append([idx, link, tov, lat, lon, length, so, Btm, ChSlp, Tw, Twcc, n, ncc, gstr])

        pretty_print_rows(rows, headers)

        # Optionally: show downstream chain (one step)
        print("\nDownstream chain (feature_id -> to):")
        for r in rows:
            fid = r[1]
            tofid = r[2]
            print(f"{fid} -> {tofid}")

    finally:
        ds.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_routelink_batonrouge.py /path/to/RouteLink_CONUS.nc")
        sys.exit(1)
    main(sys.argv[1])
