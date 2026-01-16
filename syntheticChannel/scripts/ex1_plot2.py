import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

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

# Read us1 and ds1 (HEC-RAS outputs)
hec_file = os.path.join('..', 'data', 'ex1_hecras', 'ex1_hecras.txt')
hec = pd.read_csv(hec_file, skiprows=17, header = None, sep=',')
hec['Date'] = pd.to_datetime(hec.iloc[:, 2], format='mixed')
hec['Seconds'] = (hec['Date'] - hec['Date'][0]).dt.total_seconds()

hecU = hec[hec.iloc[:,1]==20000]
hecD = hec[hec.iloc[:,1]==0]
us_Q = hecU.iloc[:, 3]
ds_Q = hecD.iloc[:, 3]

# Initialize arrays
hecupQ = np.zeros(len(ds_Q))
hecdownQ = np.zeros(len(ds_Q))

hecupQ[:] = us_Q
hecdownQ[:] = ds_Q

# Time axis for HEC-RAS
totalTime = 86400
timeInHr = totalTime / 3600
nHec = len(ds_Q)
timeHec = np.linspace(0, timeInHr, nHec)

plt.figure()
plt.plot(timeHec, hecupQ, 'k', label='Upstream', linewidth=1, zorder=1)
plt.plot(timeHec, hecdownQ, 'k', label='Downstream (Dynamic)', linewidth=1.6, color='b', zorder=4)

# Read NWM file
nwm_file_x = os.path.join('..', 'data', 'ex1_nwm', 't_dsQ_CNX')
nwm_file_t = os.path.join('..', 'data', 'ex1_nwm', 't_dsQ_CNT')
nwm_x = pd.read_csv(nwm_file_x, header=None, skiprows=1, sep=' ').to_numpy()
nwm_t = pd.read_csv(nwm_file_t, header=None, skiprows=1, sep=' ').to_numpy()

# Loop over segments
segments = [d for d in os.listdir(caseName) if d.startswith("segment")]
for seg in segments:
    segmentPath = os.path.join(caseName, seg, 'run')
    if not os.path.exists(segmentPath):
        continue

    timeValues, upstreamQ, downstreamQ, upstreamH, downstreamH = [], [], [], [], []

    timeDirs = [d for d in os.listdir(segmentPath) if os.path.isdir(os.path.join(segmentPath, d))]
    for tdir in timeDirs:
        if tdir.startswith('.'):
            continue
        try:
            timeVal = float(tdir)
        except ValueError:
            continue

        timeValues.append(timeVal)

        # Q file
        Qfile = os.path.join(segmentPath, tdir, 'Q')
        if os.path.exists(Qfile):
            Qdata = np.loadtxt(Qfile)
            if Qdata.size > 0:
                upstreamQ.append(Qdata[0])
                downstreamQ.append(Qdata[-1])
            else:
                upstreamQ.append(np.nan)
                downstreamQ.append(np.nan)
        else:
            upstreamQ.append(np.nan)
            downstreamQ.append(np.nan)

        # h file
        Hfile = os.path.join(segmentPath, tdir, 'h')
        if os.path.exists(Hfile):
            Hdata = np.loadtxt(Hfile)
            if Hdata.size > 0:
                upstreamH.append(Hdata[0])
                downstreamH.append(Hdata[-1])
            else:
                upstreamH.append(np.nan)
                downstreamH.append(np.nan)
        else:
            upstreamH.append(np.nan)
            downstreamH.append(np.nan)

    # Sort by time
    timeValues = np.array(timeValues)
    sortIdx = np.argsort(timeValues)
    timeValues = timeValues[sortIdx]
    upstreamQ = np.array(upstreamQ)[sortIdx]
    downstreamQ = np.array(downstreamQ)[sortIdx]

    # Plot results
    plt.plot(timeValues/3600, downstreamQ, 'k', label='Downstream (Present)', linewidth=1.6, zorder=5)
    plt.plot(nwm_x[:, 0]/3600, nwm_x[:, 1], label=r'Downstream (CNS - Beg $\mathit{et\,al.}$ 2023)', linewidth=1.6, color='r', zorder=3)
    plt.plot(nwm_t[:, 0] / 3600, nwm_t[:, 1], label=r'Downstream (CNT - Beg $\mathit{et\,al.}$ 2023)', linewidth=1.6, color='orange', zorder=2)

# Axis formatting
mass_in = 7200 * 242.5
mass_out = simpson(downstreamQ, timeValues)
mass_balance_error = ((mass_in - mass_out) / mass_in) * 100

mass_in = 7200 * 242.5 / 3600
mass_out = simpson(ds_Q, timeHec)
mass_balance_error_d = ((mass_in - mass_out) / mass_in) * 100

mass_in = 7200 * 242.5
mass_out = simpson(nwm_x[:, 1], nwm_x[:, 0])
mass_balance_error_nwm = ((mass_in - mass_out) / mass_in) * 100

# --- Mean Percentage Error (vs HEC-RAS) ---
# Interpolate both datasets to common time grid (0–24 hr)
common_time_hr = np.linspace(0, 24, 200)
ds_interp = np.interp(common_time_hr, timeHec, ds_Q)
model_interp = np.interp(common_time_hr, timeValues/3600, downstreamQ)
cns_interp = np.interp(common_time_hr, nwm_x[:, 0]/3600, nwm_x[:, 1])

mpe = np.mean(np.abs((model_interp - ds_interp) / ds_interp)) * 100
rmse_p = rmse(ds_interp, model_interp)
rmse_cns = rmse(ds_interp, cns_interp)

# --- Print errors ---
print(f"Mass Balance Error (Dynamic): {mass_balance_error_d:.3f}%")
print(f"Mass Balance Error (Present): {mass_balance_error:.3f}%")
print(f"Mass Balance Error (CNS): {mass_balance_error_nwm:.3f}%")
print(f"Root Mean Square Error (Present): {rmse_p:.4f}cms")
print(f"Root Mean Square Error (CNS): {rmse_cns:.4f}cms")


plt.xlim([0, 24])
plt.ylim([19, 26])
plt.xlabel('Time (hr)', fontsize=12, fontname='Helvetica')
plt.ylabel(r'Discharge $(cms)$', fontsize=12, fontname='Helvetica')
plt.legend()
# plt.savefig('ex1.pdf')
plt.show()