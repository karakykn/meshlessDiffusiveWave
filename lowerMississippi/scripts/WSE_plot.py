import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from datetime import datetime, timedelta

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

def mpe(exact, approx):
    return np.sum(np.abs((exact - approx) / exact)) / len(exact) * 100

def compute_errors_wse(site_name, obs_df, present_df, cnx_df=None):
    obs_df = obs_df.dropna(subset=['stage-m']).reset_index(drop=True)

    # Mean Percentage Error (interpolate present and CNX onto observation time grid)
    present_interp = np.interp(obs_df['seconds'], present_df['seconds'], present_df['stage elev-m'])
    mpe_present = mpe(obs_df['stage-m'], present_interp)
    rmse_present = rmse(obs_df['stage-m'], present_interp)

    # mpe_cnx = None
    # if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
    #     cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    #     mpe_cnx = np.mean(np.abs((obs_df['Q-cms'] - cnx_interp) / obs_df['Q-cms'])) * 100

    print(f"---{site_name}---")
    print(f"Mean Percentage Error (Present): {mpe_present:.2f}%")
    print(f"Root Mean Square Error (Present): {rmse_present:.2f} m")
    print("")

def make_present_df(run_dir, node_id, start_date = pd.Timestamp(1500, 1, 1, 0, 0)):
    secs, h = collect_h_values(run_dir, node=node_id)
    secs, dis = collect_Q_values(run_dir, node=node_id)
    # sort by seconds (times may already be sorted, but ensure safety)
    order = np.argsort(secs)
    secs = secs[order]
    dis = dis[order]
    h = h[order]
    slopes = np.loadtxt(f'{run_dir}/../geo/slopes')
    locs = np.loadtxt(f'{run_dir}/../geo/nodes')
    slopes_av = np.zeros(len(slopes)-1)
    elev_thal = np.zeros(len(slopes))
    elev_thal[-1] = -22.0501 #thalweg elebation at the downstream
    slopes_av[:] = (slopes[1:]+ slopes[:-1]) / 2
    for i in range(len(slopes)-2, -1, -1):
        elev_thal[i] = (locs[i+1] - locs[i]) * slopes_av[i] + elev_thal[i+1]

    df = pd.DataFrame({
        'seconds': secs,
        'hours': secs / 3600,
        'flow depth-m': h,
        'stage elev-m': elev_thal[node_id] + h,
        'discharge-cms': dis
    })
    # create Date column
    if start_date != pd.Timestamp(1500, 1, 1, 0, 0):
        df['Date'] = [start_date + timedelta(seconds=float(s)) for s in df['seconds']]
    return df

def collect_h_values(base_dir="run", node=-1, verbose=False):
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
        qpath = os.path.join(base_dir, name, "h")
        if not os.path.isfile(qpath):
            skipped.append((name, "no h file"))
            if verbose:
                print(f"Skipping {name}: h file not found")
            continue

        try:
            data = np.loadtxt(qpath)
        except Exception as e:
            skipped.append((name, f"load error: {e}"))
            if verbose:
                print(f"Skipping {name}: error loading h: {e}")
            continue

        if data.size == 0:
            skipped.append((name, "empty h"))
            if verbose:
                print(f"Skipping {name}: h is empty")
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
    "font.family": "Helvetica",   # Global font
    "font.size": 12,              # Default font size
    "axes.labelsize": 12,         # Axis label size
    "axes.titlesize": 14,         # Title size
    "legend.fontsize": 11,        # Legend font size
    "xtick.labelsize": 11,        # X-tick label size
    "ytick.labelsize": 11,        # Y-tick label size
    "lines.linewidth": 1.2        # Default line width
})

# Case folder
caseName = '../'

fPath = f'../data/hec/'
observation_tarbert = pd.read_csv(f'{fPath}/../../segment0/geo/boundary_Q', sep='\t', header=None, names=['t', 'Q'])
observation_belle = pd.read_csv(f'{fPath}/../usgs_belle', sep='\t', skiprows=26)
observation_belle['stage-m'] = observation_belle['14n'] / 3.281 - 6.85 / 3.281
observation_belle['Q-cms'] = observation_belle['14n.1'] / 35.31466621266132
observation_belle['20d_1'] = pd.to_datetime(observation_belle['20d'])
observation_belle['seconds'] = (observation_belle['20d_1'] - observation_belle['20d_1'].iloc[0]).dt.total_seconds()
observation_baton = pd.read_csv(f'{fPath}/../usgs_baton', sep='\t', skiprows=27)
observation_baton['Q-cms'] = observation_baton['14n.3'] / 35.31466621266132
observation_baton['20d_1'] = pd.to_datetime(observation_baton['20d'])
observation_baton['stage-m'] = observation_baton['14n.2'] / 3.281
observation_baton['seconds'] = (observation_baton['20d_1'] - observation_baton['20d_1'].iloc[0]).dt.total_seconds()

start = pd.Timestamp(2011, 1, 1)
present_baton = make_present_df('../segment0/run', node_id=28, start_date = start)
present_belle = make_present_df('../segment0/run', node_id=72, start_date = start)

obs_headofpasses = pd.read_csv(f'../data/headofpasses_stage.csv', skiprows=25, sep=',', names=['d', 'stage-ft'])
obs_headofpasses['Date'] = pd.to_datetime(obs_headofpasses["d"], format='mixed')
obs_headofpasses['stage-ft'] = pd.to_numeric(obs_headofpasses['stage-ft'], errors='coerce')
obs_headofpasses = obs_headofpasses.dropna()
obs_headofpasses['stage-m'] = obs_headofpasses['stage-ft'] / 3.281

fig, axes = plt.subplots(2, 1, figsize=(7, 8/3*2), sharex=True)
axes[0].plot(observation_baton['20d_1'], observation_baton['stage-m'], label='USGS (Gage: 07374000)', marker='o', linestyle='None', markersize=2, color='b', zorder=9)
axes[0].plot(present_baton['Date'], present_baton['stage elev-m'], label='Present', linewidth=1, color='k', zorder=1)
axes[0].set_ylabel(r'Water surface elevation (m)', fontsize=12, fontname='Helvetica')
axes[0].legend()

axes[1].plot(observation_belle['20d_1'], observation_belle['stage-m'], label='USGS (Gage: 07374000)', marker='o', linestyle='None', markersize=2, color='b', zorder=9)
axes[1].plot(present_belle['Date'], present_belle['stage elev-m'], label='Present', linewidth=1, color='k', zorder=1)
axes[1].plot(obs_headofpasses['Date'], obs_headofpasses['stage-m'], label='Present', linewidth=1, color='r', zorder=1)
axes[1].set_ylabel(r'Water surface elevation (m)', fontsize=12, fontname='Helvetica')
axes[1].legend()

plt.xlabel('Date', fontsize=12, fontname='Helvetica')

compute_errors_wse("Baton Rouge", observation_baton, present_baton)
compute_errors_wse("Belle Chasse", observation_belle, present_belle)

plt.show()