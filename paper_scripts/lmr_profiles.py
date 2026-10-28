"""L.M.R. 2011 (lowerMississippi_paper) with the default source.py: profiles of Q, G and P = Q - G along the channel
at selected days, and the discharges at the nodes on both sides of each spillway, for the lateral-outflow figures.
The case is copied to a temporary folder, so the case folder is not changed.
Output (this folder): lmr_profiles.csv, lmr_series.csv
"""
import os, sys, shutil, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
sys.path.insert(0, CODE)
import source

CASE = 'lowerMississippi_paper'
DAYS = [120.0, 130.0, 134.0, 137.0]          # model time 0 = 1 January 2011 00:00
NODES = [10, 11, 58, 59, 60, 76, 77, 28, 72]  # 0-based: around Morganza (38.0 km), Bonnet Carre (283.4 km, node 59), Bohemia (390.0 km); gauges 28, 72

tmp = tempfile.mkdtemp(prefix='lmr_lat_')
shutil.copytree(os.path.join(CODE, CASE), os.path.join(tmp, CASE), ignore=shutil.ignore_patterns('scripts', 'data', '*.log'))
os.chdir(tmp)
net = source.Network(CASE)
net.initial_depths()
seg = net.segments[0]
x = seg.geo['nodes'].copy()
prof, rows = {}, []
targets = [d * 86400.0 for d in DAYS]
k = 0
while net.time < net.simEndTime:
    Qold, Gold, told = seg.Q.copy(), (seg.G_new.copy() if seg.G_new is not None else np.zeros_like(seg.Q)), net.time
    net.step(0.5)
    G = seg.G_new if seg.G_new is not None else np.zeros_like(seg.Q)
    rows.append([net.time] + list(seg.Q[NODES]) + list(G[NODES]))
    while k < len(targets) and net.time >= targets[k]:
        w = (targets[k] - told) / (net.time - told)          # linear in time between the two steps
        Qk = (1 - w) * Qold + w * seg.Q
        Gk = (1 - w) * Gold + w * G
        prof[DAYS[k]] = (Qk, Gk, seg.h.copy())
        k += 1
    if net.time > 200 * 86400:
        break
cols = ['x_m']
data = [x]
for d in DAYS:
    Qk, Gk, hk = prof[d]
    cols += [f'Q_d{d:g}', f'G_d{d:g}', f'P_d{d:g}']
    data += [Qk, Gk, Qk - Gk]
np.savetxt(os.path.join(HERE, 'lmr_profiles.csv'), np.column_stack(data), delimiter=',', header=','.join(cols), comments='', fmt='%.6f')
hdr = ['t_s'] + [f'Q_{n}' for n in NODES] + [f'G_{n}' for n in NODES]
np.savetxt(os.path.join(HERE, 'lmr_series.csv'), np.array(rows), delimiter=',', header=','.join(hdr), comments='', fmt='%.6f')
shutil.rmtree(tmp)
print('steps', len(rows), 'end day', net.time / 86400)
