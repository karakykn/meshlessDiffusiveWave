"""L.M.R. 2011, point lateral outflows (SM, Representation of Lateral Outflows): (a) G and (b) Q and Q - G
along the river on 18 May 2011; (c) decrease in discharge across Morganza and Bonnet Carre Spillways and the
spillway outflows. Data from lmr_profiles.py (default source.py, lowerMississippi_paper, theta 0.5).
    python3 fig_lateral_lmr.py [--out folder]
Output: lateral_lmr.pdf/.png (default: this folder)
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from paper_style import PAPER_STYLE, N_RAMP, REFERENCE_GREY, savefig_safe
import figure_tools as ft

out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else HERE
plt.rcdefaults()
plt.rcParams.update(PAPER_STYLE)
DAY = 137                                   # 18 May 2011 00:00 (model time 0 = 1 January 2011)
SPILL = [('Morganza', 38.0), ('Bonnet Carre', 283.426), ('Bohemia', 390.0)]

hdr = open(os.path.join(HERE, 'lmr_profiles.csv')).readline().strip().split(',')
pr = np.loadtxt(os.path.join(HERE, 'lmr_profiles.csv'), delimiter=',', skiprows=1)
x = pr[:, 0] / 1000
Q, G, P = (pr[:, hdr.index(f'{v}_d{DAY}')] for v in 'QGP')
sh = open(os.path.join(HERE, 'lmr_series.csv')).readline().strip().split(',')
se = np.loadtxt(os.path.join(HERE, 'lmr_series.csv'), delimiter=',', skiprows=1)
t = ft.to_dates(se[:, 0], 2011)
c = lambda n: se[:, sh.index(n)]
dM, qM = c('Q_10') - c('Q_11'), -(c('G_11') - c('G_10'))         # Morganza between nodes 10 and 11 (0-based)
dB, qB = c('Q_58') - c('Q_60'), -(c('G_60') - c('G_58'))         # Bonnet Carre at node 59

fig = plt.figure(figsize=(3.5, 6.5))
outer = fig.add_gridspec(2, 1, height_ratios=[2.05, 1.0], hspace=0.24, left=0.21, right=0.97, top=0.99, bottom=0.075)
top = outer[0].subgridspec(2, 1, hspace=0.07)
a = fig.add_subplot(top[0])
b = fig.add_subplot(top[1], sharex=a)
cc = fig.add_subplot(outer[1])
plt.setp(a.get_xticklabels(), visible=False)
for ax in (a, b):
    for _, xs in SPILL:
        ax.axvline(xs, color=REFERENCE_GREY, lw=0.8, zorder=0)
    ax.set_xlim(0, x[-1])
# (a) G
a.plot(x, G, color='black', lw=1.0, marker='o', ms=2.6, mfc='white', mec='black', mew=0.6, drawstyle='default', zorder=3)
a.text(45, -4700, 'Morganza', ha='left', va='bottom', fontsize=8)
a.text(278, -5800, 'Bonnet Carre', ha='right', va='top', fontsize=8)
a.text(385, -13450, 'Bohemia', ha='right', va='bottom', fontsize=8)
a.set_ylim(-17000, 3500)
a.set_ylabel(r'$G$ (m$^3$/s)')
# (b) Q and P
b.plot(x, Q, color=N_RAMP[3], lw=1.3, marker='o', ms=2.6, mfc=N_RAMP[3], mec=N_RAMP[3], label=r'$Q$', zorder=4)
b.plot(x, P, color=N_RAMP[1], lw=1.3, ls='-.', label=r'$Q-G$', zorder=3)
b.set_ylim(25000, 48000)
b.set_ylabel(r'Discharge (m$^3$/s)')
b.set_xlabel('Distance from Tarbert Landing (km)')
b.legend(loc='lower left', fontsize=8, handlelength=2.6)
# (c) decrease across the spillways
m = (se[:, 0] >= 120 * 86400) & (se[:, 0] <= 197 * 86400)
cc.plot(t[m], qM[m], color=REFERENCE_GREY, lw=3.4, solid_capstyle='butt', label=r'Outflow, $|Q_\ell|$', zorder=1)
cc.plot(t[m], qB[m], color=REFERENCE_GREY, lw=3.4, solid_capstyle='butt', zorder=1)
cc.plot(t[m], dM[m], color=N_RAMP[3], lw=1.2, label='Morganza', zorder=3)
cc.plot(t[m], dB[m], color=N_RAMP[1], lw=1.2, ls='-.', label='Bonnet Carre', zorder=3)
cc.set_xlim(ft.to_dates(120 * 86400, 2011), ft.to_dates(197 * 86400, 2011))
cc.set_ylim(-500, 11000)
cc.xaxis.set_major_locator(mdates.MonthLocator())
cc.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=[8, 15, 22]))
cc.xaxis.set_major_formatter(mdates.DateFormatter('%-d %b'))
cc.set_xlabel('2011')
cc.set_ylabel(r'Decrease in discharge (m$^3$/s)')
cc.legend(loc='upper right', fontsize=8, handlelength=2.6)
for ax, l in ((a, 'a'), (b, 'b'), (cc, 'c')):
    ax.text(0.03, 0.95, f'({l})', transform=ax.transAxes, fontweight='bold', va='top', ha='left')
os.makedirs(out, exist_ok=True)
for ext in ('pdf', 'png'):
    savefig_safe(fig, os.path.join(out, 'lateral_lmr.' + ext), bbox_inches='tight')
print('Saved', os.path.join(out, 'lateral_lmr.pdf'), 'x end', x[-1])
