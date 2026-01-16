import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

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
hecras_dir = os.path.join('..', 'data', 'ex2_hecras.txt')
nwm_dir = os.path.join('..', 'data', 'ex2_nwm')

# --- Reference Elevations for HEC-RAS (m) ---
elevations = {
    0: {'us': 2.5, 'ds': 2.0},  # channel 1
    1: {'us': 3.0, 'ds': 2.0},  # channel 2
    2: {'us': 2.0, 'ds': 0.0}   # channel 3
}

# --- Discover All Segments ---
segments = sorted([d for d in os.listdir(caseName) if d.startswith("segment")])
n_segments = len(segments)

# --- Setup Subplots (2 columns: discharge & depth) ---
fig, axes = plt.subplots(n_segments, 1, figsize=(7, 3 * n_segments), sharex='col')
if n_segments == 1:
    axes = np.array([axes])  # make 2D

# Storage for error calculations
upstreamQs_all = {}
downstreamQs_all = {}
time_all = {}
hec_us_all = {}

# --- Loop Over Segments and Plot ---
for idx, seg in enumerate(segments):
    seg_num = idx + 1
    axQ = axes[idx]  # discharge
    # axH = axes[idx, 1]  # depth

    # --- HEC-RAS file names ---
    hec_file = os.path.join(hecras_dir)

    if os.path.exists(hec_file):
        hec = pd.read_csv(hec_file, sep=r'\s+', header=None, skiprows=13)
        hec = hec.iloc[2:,:7]
        hec['Q-cms'] = hec.iloc[:, 4].astype(float)
        hec['h-m'] = hec.iloc[:, 6].astype(float) - hec.iloc[:, 5].astype(float)
        hec['Date'] = pd.to_datetime(hec.iloc[:,2] + " " + hec.iloc[:,3])
        hec['Seconds'] = (hec['Date'] - hec['Date'][2]).dt.total_seconds()
        hec['Hours'] = hec['Seconds'] / 3600
        hec = hec[hec.iloc[:,0] == f'{idx+1}']

        usId = max(hec.iloc[:,1])

        us_Q = hec[hec.iloc[:,1] == usId]
        ds_Q = hec[hec.iloc[:,1] != usId]

        totalTime = 86400
        # --- Plot discharge (HEC-RAS) ---
        if idx == 2:
            axQ.plot(us_Q['Hours'], us_Q['Q-cms'], label='Upstream (Dynamic)', linewidth=.8, color='b', zorder=14)
        else:
            axQ.plot(us_Q['Hours'], us_Q['Q-cms'], 'k', label='Upstream', linewidth=.8, zorder=14)
        axQ.plot(ds_Q['Hours'], ds_Q['Q-cms'], label='Downstream (Dynamic)', linewidth=1.6, color='b', zorder=19)

        # --- Convert elevations → depths and plot (HEC-RAS) ---
        # h_us = us_h - elevations[idx]['us']
        # h_ds = ds_h - elevations[idx]['ds']

        # if idx == 2:
        #     axH.plot(us_Q['Hours'], us_Q['h-m'], 'k--', label='Upstream (Dynamic)', linewidth=.8)
        # else:
        #     axH.plot(us_Q['Hours'], us_Q['h-m'], 'k--', label='Upstream (Dynamic)', linewidth=.8)
        # axH.plot(ds_Q['Hours'], ds_Q['h-m'], 'k--', label='Downstream (Dynamic)', linewidth=1.6)

        # Save downstream HEC for MPE
        hec_us_all[idx] = (ds_Q['Seconds'].astype(float), ds_Q['Q-cms'].astype(float))

    # --- Present Model Results ---
    segmentPath = os.path.join(caseName, seg, 'run')
    if not os.path.exists(segmentPath):
        continue

    timeValues, downstreamQ, upstreamQ_present = [], [], []
    downstreamH, upstreamH_present = [], []

    timeDirs = [d for d in os.listdir(segmentPath) if os.path.isdir(os.path.join(segmentPath, d))]
    for tdir in timeDirs:
        if tdir.startswith('.'):
            continue
        try:
            timeVal = float(tdir)
        except ValueError:
            continue

        Qfile = os.path.join(segmentPath, tdir, 'Q')
        Hfile = os.path.join(segmentPath, tdir, 'h')
        if os.path.exists(Qfile):
            Qdata = np.loadtxt(Qfile)
            if Qdata.size > 0:
                downstreamQ.append(Qdata[-1])
                upstreamQ_present.append(Qdata[0])  # now included for all channels
                timeValues.append(timeVal)
        if os.path.exists(Hfile):
            Hdata = np.loadtxt(Hfile)
            if Hdata.size > 0:
                downstreamH.append(Hdata[-1])
                upstreamH_present.append(Hdata[0])  # now included for all channels

    if timeValues:
        timeValues = np.array(timeValues)
        sortIdx = np.argsort(timeValues)
        timeValues = timeValues[sortIdx]
        downstreamQ = np.array(downstreamQ)[sortIdx]
        downstreamH = np.array(downstreamH)[sortIdx]

        upstreamQ_present = np.array(upstreamQ_present)[sortIdx] if upstreamQ_present else []
        upstreamH_present = np.array(upstreamH_present)[sortIdx] if upstreamH_present else []

        # Store for error calculations
        downstreamQs_all[idx] = downstreamQ
        time_all[idx] = timeValues
        if len(upstreamQ_present):
            upstreamQs_all[idx] = upstreamQ_present

        # --- Plot Upstream (now for all channels) ---
        if idx == 2 and len(upstreamQ_present) == len(timeValues):
            axQ.plot(timeValues / 3600, upstreamQ_present, 'k', label='Upstream (Present)', linewidth=.8, zorder=15)
        # if len(upstreamH_present) == len(timeValues):
        #     axH.plot(timeValues / 3600, upstreamH_present, 'k', label='Upstream (Present)', linewidth=.8)

        # --- Plot Downstream ---
        axQ.plot(timeValues / 3600, downstreamQ, 'k', label='Downstream (Present)', linewidth=1.6, zorder = 20)
        # axH.plot(timeValues / 3600, downstreamH, 'k', label='Downstream (Present)', linewidth=1.6)

    # --- NWM Results (discharge only) ---
    x_file = os.path.join(nwm_dir, f't_dsQ_CNX')
    t_file = os.path.join(nwm_dir, f't_dsQ_CNT')

    if os.path.exists(x_file):
        begtime_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 0] / 3600
        nwmus_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 1 + idx]
        nwmds_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 4 + idx]
        begtime_t = pd.read_csv(t_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 0] / 3600
        nwmus_t = pd.read_csv(t_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 1 + idx]
        nwmds_t = pd.read_csv(t_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 4 + idx]

        if idx == 2:
            axQ.plot(begtime_x, nwmus_x, label=r'Upstream (CNS - Beg $\mathit{et\,al.}$ 2023)', linewidth=.8,color='r', zorder=18)
            axQ.plot(begtime_t, nwmus_t, label=r'Upstream (CNT - Beg $\mathit{et\,al.}$ 2023)',
                     linewidth=.8, color='orange', zorder=12)
        axQ.plot(begtime_x, nwmds_x, label=r'Downstream (CNS - Beg $\mathit{et\,al.}$ 2023)', linewidth=1.6,color='r',zorder=17)
        axQ.plot(begtime_t, nwmds_t, label=r'Downstream (CNT - Beg $\mathit{et\,al.}$ 2023)',
                 linewidth=1.6, color='orange', zorder=10)

    # --- Formatting ---
    axQ.set_xlim([0, 24])
    # axH.set_xlim([0, 24])
    axQ.set_ylabel(r'Discharge (cms)')
    # axH.set_ylabel(r'h (m)')
    axQ.set_title(f'Channel {seg_num}')
    # axH.set_title(f'Channel {seg_num}')
    axQ.legend(frameon=False, fontsize=9) if idx == 2 else axQ.legend(frameon=False, fontsize=12)
    # axH.legend(frameon=False)




# --- Global Axis Labels ---
axes[-1].set_xlabel('Time ($hr$)')
# axes[-1, 1].set_xlabel('Time ($hr$)')
plt.tight_layout()

# ======================================================
# === ERROR CALCULATIONS (Mass Balance + Composite MPE) ==
# ======================================================

# --- Mass balance: in = sum of upstreamQ of ch1 & ch2, out = downstream of ch3 ---
if 0 in downstreamQs_all and 1 in downstreamQs_all and 2 in downstreamQs_all:
    t3 = time_all[2]
    Qup1 = np.interp(t3, time_all[0], downstreamQs_all[0])
    Qup2 = np.interp(t3, time_all[1], downstreamQs_all[1])
    Qin_total = Qup1 + Qup2
    Qout = downstreamQs_all[2]

    mass_in = simpson(Qin_total, x=t3)
    mass_out = simpson(Qout, x=t3)
    mass_balance_error = ((mass_in - mass_out) / mass_in) * 100
else:
    mass_balance_error = np.nan

# --- Composite MPE ---
mpe_values = []
common_time_hr = np.linspace(0, 24, 300)
for idx in hec_us_all:
    if idx in upstreamQs_all:
        hec_t, hec_q = hec_us_all[idx]
        model_t = time_all[idx] / 3600
        model_q = upstreamQs_all[idx]

        hec_interp = np.interp(common_time_hr, hec_t, hec_q)
        model_interp = np.interp(common_time_hr, model_t, model_q)

        mpe = np.mean(np.abs((model_interp - hec_interp) / hec_interp)) * 100
        mpe_values.append(mpe)

composite_mpe = np.mean(mpe_values) if mpe_values else np.nan

# print(f"Mass Balance Error (Present): {mass_balance_error:.3f}%")
# print(f"Mean Percentage Error (MPE) at channel 3 upstream: {composite_mpe:.3f}%")
# plt.savefig('ex2.pdf')
plt.show()

for idx in range(3):
    print(f'---------- Channel {idx} ----------')

    mass_in = simpson(upstreamQs_all[idx], time_all[idx])
    mass_out = simpson(downstreamQs_all[idx], time_all[idx])
    mbe_p = (mass_in - mass_out) / mass_in * 100
    print(f'Mass balance error (Present): {mbe_p:.2f}%')

    x_file = os.path.join(nwm_dir, f't_dsQ_CNX')
    t_file = os.path.join(nwm_dir, f't_dsQ_CNT')

    if os.path.exists(x_file):
        begtime_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 0] / 3600
        nwmus_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 1 + idx]
        nwmds_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 4 + idx]

    mass_in = simpson(nwmus_x, begtime_x)
    mass_out = simpson(nwmds_x, begtime_x)
    mbe_p = (mass_in - mass_out) / mass_in * 100
    print(f'Mass balance error (CNS): {mbe_p:.2f}%')

    ds_file = os.path.join(hecras_dir, f'ds{idx+1}.txt')
    ds = pd.read_csv(ds_file, sep=r'\s+', header=None)
    ds_Q = pd.to_numeric(ds.iloc[1:, 4], errors='coerce').to_numpy()

    common_time_hr = np.linspace(0, 24, 200)
    ds_interp = np.interp(common_time_hr, timeHec, ds_Q)
    model_interp = np.interp(common_time_hr, time_all[idx] / 3600, downstreamQs_all[idx])
    cns_interp = np.interp(common_time_hr, begtime_x / 3600, nwmds_x)

    rmse_p = rmse(ds_interp, model_interp)
    rmse_cns = rmse(ds_interp, cns_interp)

    print(f'Root mean square error (Present): {rmse_p:.4f}cms')
    print(f'Root mean square error (CNS): {rmse_cns:.4f}cms')