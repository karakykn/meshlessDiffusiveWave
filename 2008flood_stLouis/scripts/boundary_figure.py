"""Midwest flood, 2008: boundary conditions of the model (figure 'inoutex4' of the paper), from the model input:
  ex4_inflows  inflows at Saverton (segment 0), Valley City (segment 2) and Hermann (segment 1), and flow depth at
               Chester (segment 4, outlet); two panels, one time axis

  python boundary_figure.py                plot (PDF and PNG next to this script)
  python boundary_figure.py --manuscript   also copy the PDF to the fig folder of the manuscript
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(CASE_DIR))
import figure_tools as ft

YEAR = 2008                                  # model time 0 = 1 January 2008
INFLOWS = {'Saverton': 0, 'Valley City': 2, 'Hermann': 1}      # name: segment with this upstream inflow
OUTLET = ('Chester', 4)


def geo(seg, f):
    return np.loadtxt(os.path.join(CASE_DIR, f'segment{seg}', 'geo', f))


def main():
    args = ft.arguments(__doc__)
    fig, (ax1, ax2) = ft.stacked_figure(2)
    qmax = 0
    for k, (name, seg) in enumerate(INFLOWS.items()):
        q = geo(seg, 'boundary_Q')
        ft.plot_input(ax1, ft.to_dates(q[:, 0], YEAR), q[:, 1], k, name)
        qmax = max(qmax, q[:, 1].max())
    ax1.set_ylim(0, 1.45 * qmax)
    ax1.set_ylabel('Discharge (m³/s)')
    name, seg = OUTLET
    h = geo(seg, 'boundary_h')
    ft.plot_input(ax2, ft.to_dates(h[:, 0], YEAR), h[:, 1], 0, name)
    lo, hi = h[:, 1].min(), h[:, 1].max()
    ax2.set_ylim(lo - 0.1 * (hi - lo), hi + 0.45 * (hi - lo))
    ax2.set_ylabel('Flow depth (m)')
    ft.date_axis(ax2, YEAR, every=2)
    ax2.set_xlabel(str(YEAR))
    for ax in (ax1, ax2):
        ax.legend(loc='upper right', fontsize=ft.LEGEND_SIZE, handlelength=2.6, borderaxespad=0.5)
    fig.align_ylabels([ax1, ax2])
    fig.tight_layout(h_pad=0.4)
    pdf = ft.save(fig, args.out or HERE, 'ex4_inflows')
    if args.manuscript:
        ft.copy_to_manuscript([pdf])


if __name__ == '__main__':
    main()
