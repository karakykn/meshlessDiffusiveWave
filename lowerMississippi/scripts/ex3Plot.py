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

fPath = f'../data/hec/'
hec_n = pd.read_csv(f'{fPath}hec_ncc_belleChse.txt', sep='\t')
hec_ncc = pd.read_csv(f'{fPath}hec_ncc_belleChse.txt', sep='\t')
observation_tarbert = pd.read_csv(f'{fPath}/../../segment0/geo/boundary_Q', sep='\t', header=None, names=['t', 'Q'])
observation_belle = pd.read_csv(f'{fPath}/../usgs_belle', sep='\t', skiprows=26)
observation_belle['h-m'] = observation_belle['14n'] / 3.281
observation_belle['Q-cms'] = observation_belle['14n.1'] / 35.31466621266132
observation_belle['20d_1'] = pd.to_datetime(observation_belle['20d'])
observation_belle['seconds'] = (observation_belle['20d_1'] - observation_belle['20d_1'].iloc[0]).dt.total_seconds()
observation_baton = pd.read_csv(f'{fPath}/../usgs_baton', sep='\t', skiprows=27)
CNX_baton = pd.read_csv(f'{fPath}/../beg_data/baton_disc_beg', sep='\t')
CNX_belle = pd.read_csv(f'{fPath}/../beg_data/belle_disc_beg', sep='\t')
observation_baton['Q-cms'] = observation_baton['14n.3'] / 35.31466621266132
observation_baton['20d_1'] = pd.to_datetime(observation_baton['20d'])
observation_baton['seconds'] = (observation_baton['20d_1'] - observation_baton['20d_1'].iloc[0]).dt.total_seconds()
hec_n['Date'] = hec_n['Time and Date']
hec_ncc['Date'] = hec_ncc['Time and Date']

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
start_date = datetime(2011, 1, 1, 0, 0, 0)

# Meshless datasets
Meshless_baton = make_Meshless_df('../segment0/run', node_id=28, start_date=start_date)  # node 28 is baton
Meshless_belle = make_Meshless_df('../segment0/run', node_id=72, start_date=start_date)  # node 72 is belle chase

# ensure CNX Date columns exist like before (if t column exists)
if 't' in CNX_baton.columns:
    CNX_baton = CNX_baton.copy()
    CNX_baton['Date'] = [start_date + timedelta(seconds=float(s)) for s in CNX_baton['t']]

if 't' in CNX_belle.columns:
    CNX_belle = CNX_belle.copy()
    CNX_belle['Date'] = [start_date + timedelta(seconds=float(s)) for s in CNX_belle['t']]

# parse HEC dates (keep your original parsing safely)
for i in range(hec_n.shape[0]):
    hec_n.loc[i, 'Date'] = hec_n.loc[i, 'Time and Date'].split(' ')[0]
    hec_ncc.loc[i, 'Date'] = hec_ncc.loc[i, 'Time and Date'].split(' ')[0]

observation_belle['20d'] = pd.to_datetime(observation_belle['20d'])
observation_baton['20d'] = pd.to_datetime(observation_baton['20d'])
hec_n['DateO'] = pd.to_datetime(hec_n['Date'], format='%d%b%Y')
hec_ncc['DateO'] = pd.to_datetime(hec_ncc['Date'], format='%d%b%Y')

# ---- PLOT using subplots for Baton Rouge and Belle Chase ----
fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

# Baton Rouge (top)
ax = axes[0]
ax.plot(observation_baton['20d'], observation_baton['Q-cms'], label='USGS (Gage: 07374000)', marker='o', linestyle='None', markersize=2, color='b', zorder=9)
ax.plot(Meshless_baton['Date'], Meshless_baton['discharge-cms'], label='Meshless', linewidth=1.5, color='k', zorder=10)
if 'Q' in CNX_baton.columns:
    ax.plot(CNX_baton['Date'], CNX_baton['Q'], label='CNS - Beg $\mathit{et\,al.}$ 2023', color='r', linewidth=1.5, zorder=8)
ax.set_ylabel('Discharge (cms)')
ax.set_title('Mississippi River at Baton Rouge')
ax.legend()
# ax.grid(True)

# Belle Chase (bottom)
ax = axes[1]
ax.plot(observation_belle['20d'], observation_belle['Q-cms'], label='USGS (Gage: 07374525)', marker='o', linestyle='None', markersize=2, color='b', zorder=9)
ax.plot(Meshless_belle['Date'], Meshless_belle['discharge-cms'], label='Meshless', linewidth=1.5, color='k', zorder=10)
if 'Q' in CNX_belle.columns:
    ax.plot(CNX_belle['Date'], CNX_belle['Q'], label='CNS - Beg $\mathit{et\,al.}$ 2023', linewidth=1.5, color='r', zorder=8)
ax.set_ylabel('Discharge (cms)')
ax.set_title('Mississippi River at Belle Chasse')
ax.legend()
# ax.grid(True)

# x-limits from HEC if available
try:
    xmin = hec_n['DateO'].iloc[0]
    xmax = hec_n['DateO'].iloc[-1]
    axes[0].set_xlim(xmin, xmax)
except Exception:
    pass

plt.tight_layout()
ax.xaxis.set_major_locator(
    mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11])
)
ax.xaxis.set_minor_locator(mdates.MonthLocator())
# plt.savefig('ex3.pdf')
plt.show()

# ---- MASS BALANCE & ERROR CALCS for both sites ----
# mass_in computed from observation_tarbert as in your code
observation_tarbert = observation_tarbert[observation_tarbert['t'] < CNX_baton['t'].iloc[-1]]
mass_in = simpson(observation_tarbert['Q'], observation_tarbert['t'])

def compute_errors(site_name, obs_df, Meshless_df, cnx_df=None):
    # drop NaN discharge values
    lT = cnx_df['t'].iloc[-1]
    obs_df = obs_df[obs_df['seconds'] <= lT]
    obs_df = obs_df.dropna(subset=['Q-cms']).reset_index(drop=True)
    Meshless_df = Meshless_df[Meshless_df['seconds'] <= lT]
    Meshless_df = Meshless_df.dropna(subset=['discharge-cms']).reset_index(drop=True)
    if cnx_df is not None and 'Q' in cnx_df.columns:
        cnx_df = cnx_df.dropna(subset=['Q']).reset_index(drop=True)

    # mass out from observed
    mass_out_obs = simpson(obs_df['Q-cms'], obs_df['seconds'])
    netflux_perc = ((mass_in - mass_out_obs) / mass_in) * 100

    # Meshless mass out
    mass_out_Meshless = simpson(Meshless_df['discharge-cms'], Meshless_df['seconds'])
    mass_balance_error_Meshless = ((mass_in - mass_out_Meshless) / mass_in) * 100

    # CNX mass out (if available)
    mass_balance_error_cnx = None
    if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
        mass_out_cnx = simpson(cnx_df['Q'], cnx_df['t'])
        mass_balance_error_cnx = ((mass_in - mass_out_cnx) / mass_in) * 100

    # Mean Percentage Error (interpolate Meshless and CNX onto observation time grid)
    Meshless_interp = np.interp(obs_df['seconds'], Meshless_df['seconds'], Meshless_df['discharge-cms'])
    mpe_Meshless = mpe(obs_df['Q-cms'], Meshless_interp)
    rmse_Meshless = rmse(obs_df['Q-cms'], Meshless_interp)

    # mpe_cnx = None
    # if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
    #     cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    #     mpe_cnx = np.mean(np.abs((obs_df['Q-cms'] - cnx_interp) / obs_df['Q-cms'])) * 100

    cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    mpe_cnx = mpe(obs_df['Q-cms'], cnx_interp)
    rmse_cnx = rmse(obs_df['Q-cms'], cnx_interp)

    print(f"---{site_name}---")
    print(f"Mean Percentage Error (Meshless): {mpe_Meshless:.2f}%")
    print(f"Mean Percentage Error (CNS): {mpe_cnx:.2f}%")
    print(f"Root Mean Square Error (Meshless): {rmse_Meshless:.2f} cms")
    print(f"Root Mean Square Error (CNS): {rmse_cnx:.2f} cms")
    print("")

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

def mpe(exact, approx):
    return np.sum(np.abs((exact - approx) / exact)) / len(exact) * 100

# compute errors for Baton and Belle
compute_errors("Baton Rouge", observation_baton, Meshless_baton, CNX_baton)
compute_errors("Belle Chasse", observation_belle, Meshless_belle, CNX_belle)
