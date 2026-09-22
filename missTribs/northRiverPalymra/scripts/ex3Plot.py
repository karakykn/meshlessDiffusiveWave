import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
from scipy.integrate import simpson
import matplotlib.dates as mdates


def collect_Q_values(base_dir="run", node=-1, verbose=False):
    """
    Collect values from 'Q' files in numeric subdirectories of base_dir.

    Returns sorted arrays (times, values).
    """
    # gather numeric subdirectories
    numeric_dirs = []
    try:
        for name in os.listdir(base_dir):
            full = os.path.join(base_dir, name)
            if not os.path.isdir(full):
                continue
            try:
                t = float(name)
                numeric_dirs.append((t, name))
            except ValueError:
                if verbose:
                    print(f"Skipping non-numeric directory: {name}")
                continue
    except FileNotFoundError:
        raise FileNotFoundError(f"Base directory not found: {base_dir}")

    numeric_dirs.sort(key=lambda x: x[0])

    times = []
    values = []
    skipped = []

    for t, name in numeric_dirs:
        qpath = os.path.join(base_dir, name, "Q")
        if not os.path.isfile(qpath):
            skipped.append((name, "no Q file"))
            if verbose:
                print(f"Skipping {name}: Q file not found")
            continue

        try:
            data = np.loadtxt(qpath)
        except Exception as e:
            skipped.append((name, f"load error: {e}"))
            if verbose:
                print(f"Skipping {name}: error loading Q: {e}")
            continue

        if data.size == 0:
            skipped.append((name, "empty Q"))
            if verbose:
                print(f"Skipping {name}: Q is empty")
            continue

        try:
            if data.ndim == 1:
                val = data[node]
            else:
                val = data[node, -1]
        except IndexError:
            skipped.append((name, "node index out of range"))
            if verbose:
                print(f"Skipping {name}: node {node} out of range (shape {data.shape})")
            continue
        except Exception as e:
            skipped.append((name, f"indexing error: {e}"))
            if verbose:
                print(f"Skipping {name}: indexing error: {e}")
            continue

        times.append(float(t))
        values.append(float(val))

    if verbose and skipped:
        print(f"Skipped {len(skipped)} entries. Example skips: {skipped[:5]}")

    return np.array(times), np.array(values)

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 16,
    "legend.fontsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "lines.linewidth": 1.5
})

# helper to build Meshless DataFrame for a node
def make_Meshless_df(run_dir, node_id, start_date):
    secs, dis = collect_Q_values(run_dir, node=node_id)
    # sort by seconds (times may already be sorted, but ensure safety)
    order = np.argsort(secs)
    secs = secs[order]
    dis = dis[order]

    df = pd.DataFrame({
        'seconds': secs,
        'discharge-cms': dis
    })
    # create Date column
    df['Date'] = [start_date + timedelta(seconds=float(s)) for s in df['seconds']]
    return df

# choose start date consistent with your pipeline
start_date = datetime(2008, 6, 1, 0, 0, 0)

# Meshless datasets
us = make_Meshless_df('../segment0/run', node_id=0, start_date=start_date)
ds = make_Meshless_df('../segment0/run', node_id=-1, start_date=start_date)

plt.plot(us['Date'], us['discharge-cms'])
plt.plot(ds['Date'], ds['discharge-cms'])
plt.ylabel('Discharge m3s')
plt.title('north river at palymra')
plt.show()

ds = ds.set_index("Date").sort_index()

# hourly resample
df_hourly = ds.resample("h").mean()

# interpolate discharge
df_hourly["discharge-cms"] = df_hourly["discharge-cms"].interpolate(method="time")

# save to CSV with Date as a normal column again
df_hourly_reset = df_hourly.reset_index()
df_hourly_reset.to_csv("northRiverPalymra.csv", index=False)

print(df_hourly_reset.head())

# plot
plt.plot(df_hourly_reset["Date"], df_hourly_reset["discharge-cms"])
plt.show()