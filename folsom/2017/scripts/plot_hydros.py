import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

def mpe(exact, approx):
    return np.sum(np.abs((exact - approx) / exact)) / len(exact) * 100

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
# hecras_dir = os.path.join('..', 'data', 'ex2')
# nwm_dir = os.path.join('..', 'data', 'ex2_nwm')

# --- Reference Elevations for HEC-RAS (m) ---
# elevations = {
#     0: {'us': 2.5, 'ds': 2.0},  # channel 1
#     1: {'us': 3.0, 'ds': 2.0},  # channel 2
#     2: {'us': 2.0, 'ds': 0.0}   # channel 3
# }

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
# hec_us_all = {}

titr = ['American River (Folsom to Sacramento)',
        'Sacramento River (Verona to Sacramento)',
        'Sacramento River (Sacramento to AB Delta)']
# --- Loop Over Segments and Plot ---
for idx, seg in enumerate(segments):
    seg_num = idx + 1
    axQ = axes[idx]  # discharge
    # axH = axes[idx, 1]  # depth

    if idx ==2:
        df = pd.read_csv(f'../data/hydros/abdelta_usgs.csv', skiprows=24, sep='\t')
        df['Date'] = pd.to_datetime(df.iloc[:, 2], format='mixed')
        df['Q-cms'] = df.iloc[:, 3] / 35.31466621266132
        df['Seconds'] = (df['Date'] - df['Date'].iloc[0]).dt.total_seconds()
        df['Hours'] = df['Seconds'] / 3600

        axQ.plot(df['Hours'], df['Q-cms'], label='USGS (Gage: 11447890)', marker='o', linestyle='None', markersize=1.6, color='b', zorder=3)

    # --- Meshless Model Results ---
    segmentPath = os.path.join(caseName, seg, 'run')
    if not os.path.exists(segmentPath):
        continue

    timeValues, downstreamQ, upstreamQ_Meshless = [], [], []
    downstreamH, upstreamH_Meshless = [], []

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
                upstreamQ_Meshless.append(Qdata[0])  # now included for all channels
                timeValues.append(timeVal)
        if os.path.exists(Hfile):
            Hdata = np.loadtxt(Hfile)
            if Hdata.size > 0:
                downstreamH.append(Hdata[-1])
                upstreamH_Meshless.append(Hdata[0])  # now included for all channels

    if timeValues:
        timeValues = np.array(timeValues)
        sortIdx = np.argsort(timeValues)
        timeValues = timeValues[sortIdx]
        downstreamQ = np.array(downstreamQ)[sortIdx]
        downstreamH = np.array(downstreamH)[sortIdx]

        upstreamQ_Meshless = np.array(upstreamQ_Meshless)[sortIdx] if upstreamQ_Meshless else []
        upstreamH_Meshless = np.array(upstreamH_Meshless)[sortIdx] if upstreamH_Meshless else []

        # Store for error calculations
        downstreamQs_all[idx] = downstreamQ
        time_all[idx] = timeValues
        if len(upstreamQ_Meshless):
            upstreamQs_all[idx] = upstreamQ_Meshless

        # --- Plot Upstream (now for all channels) ---
        if idx == 2 and len(upstreamQ_Meshless) == len(timeValues):
            axQ.plot(timeValues / 3600, upstreamQ_Meshless, 'k', label='Upstream (Meshless)', linewidth=.8, zorder=15)

        else:
            axQ.plot(timeValues / 3600, upstreamQ_Meshless, 'k', label='Upstream', linewidth=.8, zorder=14)
        # if len(upstreamH_Meshless) == len(timeValues):
        #     axH.plot(timeValues / 3600, upstreamH_Meshless, 'k', label='Upstream (Meshless)', linewidth=.8)


        # --- Plot Downstream ---
        axQ.plot(timeValues / 3600, downstreamQ, 'k', label='Downstream (Meshless)', linewidth=1.6, zorder = 20)
        # axH.plot(timeValues / 3600, downstreamH, 'k', label='Downstream (Meshless)', linewidth=

        if idx == 1:
            # yyy = np.ones(len(timeValues)) * 2200
            axes[2].plot(timeValues / 3600, downstreamQ, 'r', label='D', linewidth=1.6, zorder=20)

    # --- NWM Results (discharge only) ---
    # x_file = os.path.join(nwm_dir, f't_dsQ_CNX')
    # t_file = os.path.join(nwm_dir, f't_dsQ_CNT')

    # if os.path.exists(us_file) and os.path.exists(ds_file):
    #     begtime_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 0] / 3600
    #     nwmus_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 1 + idx]
    #     nwmds_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 4 + idx]
    #     begtime_t = pd.read_csv(t_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 0] / 3600
    #     nwmus_t = pd.read_csv(t_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 1 + idx]
    #     nwmds_t = pd.read_csv(t_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 4 + idx]

        # if idx == 2:
        #     axQ.plot(begtime_x, nwmus_x, label=r'Upstream (CNS - Beg $\mathit{et\,al.}$ 2023)', linewidth=.8,color='r', zorder=18)
        #     axQ.plot(begtime_t, nwmus_t, label=r'Upstream (CNT - Beg $\mathit{et\,al.}$ 2023)',
        #              linewidth=.8, color='orange', zorder=12)
        # axQ.plot(begtime_x, nwmds_x, label=r'Downstream (CNS - Beg $\mathit{et\,al.}$ 2023)', linewidth=1.6,color='r',zorder=17)
        # axQ.plot(begtime_t, nwmds_t, label=r'Downstream (CNT - Beg $\mathit{et\,al.}$ 2023)',
        #          linewidth=1.6, color='orange', zorder=10)

    # --- Formatting ---
    # axQ.set_xlim([0, 24])
    # axH.set_xlim([0, 24])
    axQ.set_ylabel(r'Discharge (cms)')
    # axH.set_ylabel(r'h (m)')
    axQ.set_title(f'{titr[idx]}')
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
# for idx in hec_us_all:
#     if idx in upstreamQs_all:
#         hec_t, hec_q = hec_us_all[idx]
#         model_t = time_all[idx] / 3600
#         model_q = upstreamQs_all[idx]
#
#         hec_interp = np.interp(common_time_hr, hec_t, hec_q)
#         model_interp = np.interp(common_time_hr, model_t, model_q)
#
#         mpe_ = np.mean(np.abs((model_interp - hec_interp) / hec_interp)) * 100
#         mpe_values.append(mpe_)

composite_mpe = np.mean(mpe_values) if mpe_values else np.nan

# print(f"Mass Balance Error (Meshless): {mass_balance_error:.3f}%")
# print(f"Mean Percentage Error (MPE) at channel 3 upstream: {composite_mpe:.3f}%")
# plt.savefig('ex2.pdf')
plt.show()

for idx in range(3):
    print(f'---------- Channel {idx} ----------')

    mass_in = simpson(upstreamQs_all[idx], time_all[idx])
    mass_out = simpson(downstreamQs_all[idx], time_all[idx])
    mbe_p = (mass_in - mass_out) / mass_in * 100
    print(f'Mass balance error (Meshless): {mbe_p:.2f}%')

    # x_file = os.path.join(nwm_dir, f't_dsQ_CNX')
    # t_file = os.path.join(nwm_dir, f't_dsQ_CNT')
    #
    # if os.path.exists(us_file) and os.path.exists(ds_file):
    #     begtime_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 0] / 3600
    #     nwmus_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 1 + idx]
    #     nwmds_x = pd.read_csv(x_file, header=None, skiprows=1, sep=' ').to_numpy()[:, 4 + idx]

    # mass_in = simpson(nwmus_x, begtime_x)
    # mass_out = simpson(nwmds_x, begtime_x)
    # mbe_p = (mass_in - mass_out) / mass_in * 100
    # print(f'Mass balance error (CNS): {mbe_p:.2f}%')

    # ds_file = os.path.join(hecras_dir, f'ds{idx+1}.txt')
    # ds = pd.read_csv(ds_file, sep=r'\s+', header=None)
    # ds_Q = pd.to_numeric(ds.iloc[1:, 4], errors='coerce').to_numpy()

    # common_time_hr = np.linspace(0, 24, 200)
    # ds_interp = np.interp(common_time_hr, timeHec, ds_Q)
    # model_interp = np.interp(common_time_hr, time_all[idx] / 3600, downstreamQs_all[idx])
    # cns_interp = np.interp(common_time_hr, begtime_x / 3600, nwmds_x)

    # rmse_p = rmse(ds_interp, model_interp)
    # rmse_cns = rmse(ds_interp, cns_interp)
    #
    # mpe_p = mpe(ds_interp, model_interp)
    # mpe_cns = mpe(ds_interp, cns_interp)
    #
    # print(f'Root mean square error (Meshless): {rmse_p:.4f}cms')
    # print(f'Root mean square error (CNS): {rmse_cns:.4f}cms')
    # print(f'MPE (Meshless): {mpe_p:.2f}%')
    # print(f'MPE (CNS): {mpe_cns:.2f}%')

# us_file = os.path.join(hecras_dir, f'us{3}.txt')
# us = pd.read_csv(us_file, sep=r'\s+', header=None)
# us_Q = pd.to_numeric(us.iloc[1:, 4], errors='coerce').to_numpy()
# us_interp = np.interp(common_time_hr, timeHec, us_Q)
# model_interp = np.interp(common_time_hr, time_all[2] / 3600, upstreamQs_all[2])
# cns_interp = np.interp(common_time_hr, begtime_x / 3600, nwmus_x)

# rmse_p = rmse(us_interp, model_interp)
# rmse_cns = rmse(us_interp, cns_interp)

# mpe_p = mpe(us_interp, model_interp)
# mpe_cns = mpe(us_interp, cns_interp)

# print(f'Root mean square error (Meshless): {rmse_p:.4f}cms')
# print(f'Root mean square error (CNS): {rmse_cns:.4f}cms')
# print(f'MPE (Meshless): {mpe_p:.2f}%')
# print(f'MPE (CNS): {mpe_cns:.2f}%')