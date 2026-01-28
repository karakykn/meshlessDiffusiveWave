import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from datetime import datetime, timedelta

def compute_errors(site_name, obs_df, Meshless_df, cnx_df=None):
    # drop NaN discharge values
    # lT = cnx_df['t'].iloc[-1]
    # lT = pd.to_timedelta(obs_df['Date'].iloc[-1] - obs_df['Date'].iloc[0], unit='s')
    # obs_df = obs_df[obs_df['seconds'] <= lT]
    # obs_df = obs_df.dropna(subset=['Q-cms']).reset_index(drop=True)
    # Meshless_df = Meshless_df[Meshless_df['seconds'] <= lT]
    # Meshless_df = Meshless_df.dropna(subset=['discharge-cms']).reset_index(drop=True)
    # if cnx_df is not None and 'Q' in cnx_df.columns:
    #     cnx_df = cnx_df.dropna(subset=['Q']).reset_index(drop=True)

    # # mass out from observed
    # mass_out_obs = simpson(obs_df['Q-cms'], obs_df['seconds'])
    # netflux_perc = ((mass_in - mass_out_obs) / mass_in) * 100
    #
    # # Meshless mass out
    # mass_out_Meshless = simpson(Meshless_df['discharge-cms'], Meshless_df['seconds'])
    # mass_balance_error_Meshless = ((mass_in - mass_out_Meshless) / mass_in) * 100
    #
    # # CNX mass out (if available)
    # mass_balance_error_cnx = None
    # if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
    #     mass_out_cnx = simpson(cnx_df['Q'], cnx_df['t'])
    #     mass_balance_error_cnx = ((mass_in - mass_out_cnx) / mass_in) * 100

    obs_df['seconds'] = (obs_df['Date'] - obs_df['Date'].iloc[0]).dt.total_seconds()

    # Mean Percentage Error (interpolate Meshless and CNX onto observation time grid)
    Meshless_interp = np.interp(obs_df['seconds'], Meshless_df['seconds'], Meshless_df['discharge-cms'])
    mpe_Meshless = mpe(obs_df['Q-cms'], Meshless_interp)
    rmse_Meshless = rmse(obs_df['Q-cms'], Meshless_interp)

    # mpe_cnx = None
    # if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
    #     cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    #     mpe_cnx = np.mean(np.abs((obs_df['Q-cms'] - cnx_interp) / obs_df['Q-cms'])) * 100

    # cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    # mpe_cnx = mpe(obs_df['Q-cms'], cnx_interp)
    # rmse_cnx = rmse(obs_df['Q-cms'], cnx_interp)

    print(f"---{site_name}---")
    print(f"Mean Percentage Error (Meshless): {mpe_Meshless:.2f}%")
    # print(f"Mean Percentage Error (CNS): {mpe_cnx:.2f}%")
    print(f"Root Mean Square Error (Meshless): {rmse_Meshless:.2f} cms")
    # print(f"Root Mean Square Error (CNS): {rmse_cnx:.2f} cms")
    print("")

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

def mpe(exact, approx):
    return np.sum(np.abs((exact - approx) / exact)) / len(exact) * 100

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


# --- Global Plot Settings ---
plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "legend.fontsize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "lines.linewidth": 1.2
})

# --- Case Folder ---
caseName = '../'
df_freeport = pd.read_csv(f'../data/hydros/freeport.csv', header=None, sep='\t', skiprows=25)
df_freeport['Date'] = pd.to_datetime(df_freeport.iloc[:,2], format='mixed')
df_freeport['Q-cms'] = df_freeport.iloc[:,3] / 35.31466621266132
df_freeport['Seconds'] = (df_freeport['Date'] - df_freeport['Date'].iloc[0]).dt.total_seconds()

meshless_freeport = make_Meshless_df(f'../segment2/run', node_id=8, start_date=df_freeport['Date'].iloc[0])
meshless_freeport['Seconds'] = (meshless_freeport['Date'] - meshless_freeport['Date'].iloc[0]).dt.total_seconds()
fig, axes = plt.subplots(1, 1, figsize=(7, 4), sharex='col')
axes.plot(df_freeport['Seconds'], df_freeport['Q-cms'], label='USGS (Gage: 11447650)', marker='o', linestyle='None', markersize=1.6, color='b', zorder=3)
axes.plot(meshless_freeport['Seconds'], meshless_freeport['discharge-cms'], label='Meshless', linewidth=1.5, color='k', zorder=2)

df_yolo = pd.read_csv(f'../data/hydros/yolobypass.csv', header=None, sep='\t', skiprows=25)
df_yolo['Date'] = pd.to_datetime(df_yolo.iloc[:,2], format='mixed')
df_yolo['Q-cms'] = df_yolo.iloc[:,3] / 35.31466621266132
df_yolo['Seconds'] = (df_yolo['Date'] - df_yolo['Date'].iloc[0]).dt.total_seconds()
axes.plot(df_yolo['Seconds'], df_yolo['Q-cms'], label='Yolo outflow', linewidth=1.5, color='r', zorder=2)

df_fairoaks = pd.read_csv(f'../data/hydros/fairoaks_usgs.csv', header=None, sep='\t', skiprows=25)
df_fairoaks['Date'] = pd.to_datetime(df_yolo.iloc[:,2], format='mixed')
df_fairoaks['Q-cms'] = df_fairoaks.iloc[:,3] / 35.31466621266132
df_fairoaks['Seconds'] = (df_fairoaks['Date'] - df_fairoaks['Date'].iloc[0]).dt.total_seconds()
axes.plot(df_fairoaks['Seconds'], df_fairoaks['Q-cms'], label='Yolo outflow', linewidth=1.5, color='r', zorder=2)

# compute_errors("Freeport", df_freeport, meshless_freeport)df_yolo

plt.show()