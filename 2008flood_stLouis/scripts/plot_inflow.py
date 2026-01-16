import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

def make_present_df(run_dir, node_id, start_date):
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

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

# --- Global Plot Settings ---
plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "legend.fontsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 11,
    "lines.linewidth": 1.2
})

# --- Case Folder ---
caseName = '../'
usgs_dir = os.path.join('..', 'data', 'hydros')
# cn_dir = os.path.join('..', 'data', 'cns')

# --- Discover All Segments ---
segments = sorted([d for d in os.listdir(caseName) if d.startswith("segment")])
n_segments = len(segments)

# fig, axes = plt.subplots(n_segments, 1, figsize=(7, 2 * n_segments), sharex='col')
#
# labels = ['Mississippi, Quincy to Grafton',
#           'Missouri, St. Charles to Madison',
#           'Illinois, Valley City to Grafton',
#           'Mississippi, Grafton to Madison',
#           'Mississippi, Madison to Chester']
start_date = datetime(2008, 1, 1, 0, 0, 0)
# for i in range(n_segments):
#     present_up = make_present_df(f'../segment{i}/run', node_id=0, start_date=start_date)
#     present_down = make_present_df(f'../segment{i}/run', node_id=-1, start_date=start_date)
#
#     axes[i].plot(present_up['Date'], present_up['discharge-cms'], color='b', label='upstream')
#     axes[i].plot(present_down['Date'], present_down['discharge-cms'], color='k', label='downstream')
#     axes[i].set_title(labels[i])
#     axes[i].legend()
#
# # plt.title(f'2008 Flood')
# plt.ylabel('Discharge (cms)')
# plt.xlabel('Date')


mis_up = make_present_df(f'../segment0/run', node_id=0, start_date=start_date)
ill_up = make_present_df(f'../segment1/run', node_id=0, start_date=start_date)
miso_up = make_present_df(f'../segment2/run', node_id=0, start_date=start_date)


usgs_stLouis = pd.read_csv(usgs_dir + f'/stLouis_mississippi', skiprows=24, sep='\t')
usgs_stLouis['Date'] = pd.to_datetime(usgs_stLouis['20d'])
usgs_stLouis['Q-cms'] = usgs_stLouis['14n'] / 35.31466621266132
# usgs_stLouis = usgs_stLouis[usgs_stLouis['Date'] <= present_down['Date'].iloc[-1]]
# dat = dat[dat['date'] <= pd.to_datetime('08-18-2016 00:00', format='%m-%d-%Y %H:%M')]
present_stLouis = make_present_df(f'../segment4/run', node_id=4, start_date=start_date)

usgs_grafton = pd.read_csv(usgs_dir + f'/grafton_mississippi', skiprows=26, sep='\t')
usgs_grafton['Date'] = pd.to_datetime(usgs_grafton['20d'])
usgs_grafton['Q-cms'] = usgs_grafton['14n'] / 35.31466621266132
present_grafton = make_present_df(f'../segment3/run', node_id=0, start_date=start_date)

usgs_stCharles = pd.read_csv(usgs_dir + f'/stCharles_mississippi', skiprows=26, sep='\t')
usgs_stCharles['Date'] = pd.to_datetime(usgs_stCharles['20d'])
usgs_stCharles['Q-cms'] = usgs_stCharles['14n'] / 35.31466621266132
present_stCharles = make_present_df(f'../segment1/run', node_id=-6, start_date=start_date)

fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex='col')

# axes[0].plot(usgs_grafton['Date'], usgs_grafton['Q-cms'], label='USGS (Gage: 05587450)', marker='o', linestyle='None', markersize=1.6, color='b', zorder=3)
# axes[0].plot(present_grafton['Date'], present_grafton['discharge-cms'], label='Present', linewidth=1.5, color='k', zorder=2)

axes[0].plot(mis_up['Date'], mis_up['discharge-cms'], label='Mississippi inflow', linewidth=1.6, color='k', alpha=1)
# axes[0].plot(ill_up['Date'], ill_up['discharge-cms'], label='Illinois inflow', linewidth=.75, color='b', alpha=.5)
# axes[0].plot(miso_up['Date'], miso_up['discharge-cms'], label='Mississippi inflow', linewidth=.75, color='r', alpha=.5)

axes[0].set_title('Mississippi River at Saverton')
axes[0].set_ylabel('Discharge (cms)')
# axes[0].legend()

# axes[1].plot(usgs_sutCharles['Date'], present_stCharles['discharge-cms'], label='Present', linewidth=1.5, color='k', zorder=2)

# axes[1].plot(mis_up['Date'], mis_up['discharge-cms'], label='Mississippi inflow', linewidth=.75, color='k', alpha=.5)
axes[1].plot(ill_up['Date'], ill_up['discharge-cms'], label='Illinois inflow', linewidth=1.6, color='k', alpha=1)
# axes[1].plot(miso_up['Date'], miso_up['discharge-cms'], label='Mississippi inflow', linewidth=.75, color='r', alpha=.5)

axes[1].set_title('Illinois River at Valley City')
axes[1].set_ylabel('Discharge (cms)')
# axes[1].legend()

# axes[2].plot(usgs_stLouis['Date'], usgs_stLouis['Q-cms'], label='USGS (Gage: 07010000)', marker='o', linestyle='None', markersize=1.6, color='b', zorder=3)
# axes[2].plot(present_stLouis['Date'], present_stLouis['discharge-cms'], label='Present', linewidth=1.5, color='k', zorder=2)

# axes[2].plot(mis_up['Date'], mis_up['discharge-cms'], label='Mississippi inflow', linewidth=.75, color='k', alpha=.5)
# axes[2].plot(ill_up['Date'], ill_up['discharge-cms'], label='Illinois inflow', linewidth=.75, color='b', alpha=.5)
axes[2].plot(miso_up['Date'], miso_up['discharge-cms'], label='Mississippi inflow', linewidth=1.6, color='k', alpha=1)

axes[2].set_title('Missouri River at Hermann')
# axes[2].legend()
axes[2].set_ylabel('Discharge (cms)')
plt.tight_layout()
plt.xlabel('Date')
try:
    xmin = usgs_stLouis['Date'].iloc[0]
    xmax = usgs_stLouis['Date'].iloc[-1]
    axes[0].set_xlim(xmin, xmax)
except Exception:
    pass

df = pd.read_csv(f'../data/hydros/saverton_stage_mississippi.csv', skiprows=18, names=['d', 'wse'])
# df['Date'] = pd.to_datetime(df['d'], format='mixed')
# dt = pd.to_datetime(1704801000, unit="s")

# df['WSE'] = pd.to_numeric(df['wse'], errors='coerce')


# axes[3].plot(mis_up['Date'], mis_up['discharge-cms'], label='Mississippi inflow', linewidth=1.6, color='k', alpha=1)
# axes[0].plot(ill_up['Date'], ill_up['discharge-cms'], label='Illinois inflow', linewidth=.75, color='b', alpha=.5)
# axes[0].plot(miso_up['Date'], miso_up['discharge-cms'], label='Mississippi inflow', linewidth=.75, color='r', alpha=.5)

# axes[3].set_title('Mississippi River at Saverton')
# axes[3].set_ylabel('Water Surface Elevation (m)')


def compute_errors(site_name, obs_df, present_df, cnx_df=None):
    # drop NaN discharge values
    # lT = cnx_df['t'].iloc[-1]
    # lT = pd.to_timedelta(obs_df['Date'].iloc[-1] - obs_df['Date'].iloc[0], unit='s')
    # obs_df = obs_df[obs_df['seconds'] <= lT]
    # obs_df = obs_df.dropna(subset=['Q-cms']).reset_index(drop=True)
    # present_df = present_df[present_df['seconds'] <= lT]
    # present_df = present_df.dropna(subset=['discharge-cms']).reset_index(drop=True)
    # if cnx_df is not None and 'Q' in cnx_df.columns:
    #     cnx_df = cnx_df.dropna(subset=['Q']).reset_index(drop=True)

    # # mass out from observed
    # mass_out_obs = simpson(obs_df['Q-cms'], obs_df['seconds'])
    # netflux_perc = ((mass_in - mass_out_obs) / mass_in) * 100
    #
    # # present mass out
    # mass_out_present = simpson(present_df['discharge-cms'], present_df['seconds'])
    # mass_balance_error_present = ((mass_in - mass_out_present) / mass_in) * 100
    #
    # # CNX mass out (if available)
    # mass_balance_error_cnx = None
    # if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
    #     mass_out_cnx = simpson(cnx_df['Q'], cnx_df['t'])
    #     mass_balance_error_cnx = ((mass_in - mass_out_cnx) / mass_in) * 100

    obs_df['seconds'] = (obs_df['Date'] - obs_df['Date'].iloc[0]).dt.total_seconds()

    # Mean Percentage Error (interpolate present and CNX onto observation time grid)
    present_interp = np.interp(obs_df['seconds'], present_df['seconds'], present_df['discharge-cms'])
    mpe_present = mpe(obs_df['Q-cms'], present_interp)
    rmse_present = rmse(obs_df['Q-cms'], present_interp)

    # mpe_cnx = None
    # if cnx_df is not None and 'Q' in cnx_df.columns and 't' in cnx_df.columns:
    #     cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    #     mpe_cnx = np.mean(np.abs((obs_df['Q-cms'] - cnx_interp) / obs_df['Q-cms'])) * 100

    # cnx_interp = np.interp(obs_df['seconds'], cnx_df['t'], cnx_df['Q'])
    # mpe_cnx = mpe(obs_df['Q-cms'], cnx_interp)
    # rmse_cnx = rmse(obs_df['Q-cms'], cnx_interp)

    print(f"---{site_name}---")
    print(f"Mean Percentage Error (Present): {mpe_present:.2f}%")
    # print(f"Mean Percentage Error (CNS): {mpe_cnx:.2f}%")
    print(f"Root Mean Square Error (Present): {rmse_present:.2f} cms")
    # print(f"Root Mean Square Error (CNS): {rmse_cnx:.2f} cms")
    print("")

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

def mpe(exact, approx):
    return np.sum(np.abs((exact - approx) / exact)) / len(exact) * 100

# compute errors for Baton and Belle
compute_errors("Grafton", usgs_grafton, present_grafton)
compute_errors("St. Charles", usgs_stCharles, present_stCharles)
compute_errors("St. Louis", usgs_stLouis, present_stLouis)

# df = pd.read_csv(f'../segment1/geo/boundary_Q', header=None, sep=' ', names=['sec', 'Q-cms'])
# df['Date'] = start_date + pd.to_timedelta(df['sec'], unit='s')
# print(df)

plt.savefig('ex4plot.pdf')
plt.show()

# plt.plot(df['Date'], df['Q-cms'])
# plt.show()