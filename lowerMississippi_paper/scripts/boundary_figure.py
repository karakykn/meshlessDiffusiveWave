"""Lower Mississippi River, 2011: boundary conditions of the model (figure 'latex3' of the paper), from the model
input (segment0/geo/boundary_Q, boundary_h, laterals and lateralDatas):
  ex3_inflow   (a) inflow at Tarbert Landing and flow depth at Head of Passes (two panels, one time axis)
  ex3_lateral  (b) lateral outflows at the Morganza, Bonnet Carre and Bohemia spillways

  python boundary_figure.py                plot (PDF and PNG next to this script)
  python boundary_figure.py --manuscript   also copy the PDFs to the fig folder of the manuscript
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(CASE_DIR))
import figure_tools as ft

GEO = os.path.join(CASE_DIR, 'segment0', 'geo')
YEAR = 2011                                  # model time 0 = 1 January 2011
SPILLWAYS = {'morganza': 'Morganza Spillway', 'bonnet_carre': 'Bonnet Carre Spillway', 'bohemia': 'Bohemia Spillway'}


def laterals():
    """{label: (t (s), outflow (m3/s), positive)} from the laterals table of the case"""
    rows = [l.split(',') for l in open(os.path.join(GEO, 'laterals')).read().splitlines()[1:] if l.strip()]
    out = {}
    for name, _, qfile in rows:
        a = np.atleast_2d(np.loadtxt(os.path.join(GEO, 'lateralDatas', qfile.strip())))
        out[SPILLWAYS.get(name.strip(), name.strip())] = (a[:, 0], -a[:, 1])     # outflows are negative in the input
    return out


def main():
    args = ft.arguments(__doc__)
    q = np.loadtxt(os.path.join(GEO, 'boundary_Q'))
    h = np.loadtxt(os.path.join(GEO, 'boundary_h'))
    pdfs = []

    # (a) inflow and downstream flow depth
    fig, (ax1, ax2) = ft.stacked_figure(2)
    ft.plot_input(ax1, ft.to_dates(q[:, 0], YEAR), q[:, 1], 0, 'Tarbert Landing')
    ax1.set_ylim(0, 1.3 * q[:, 1].max())
    ax1.set_ylabel('Discharge (m³/s)')
    ft.plot_input(ax2, ft.to_dates(h[:, 0], YEAR), h[:, 1], 0, 'Head of Passes')
    lo, hi = h[:, 1].min(), h[:, 1].max()
    ax2.set_ylim(lo - 0.15 * (hi - lo), hi + 0.6 * (hi - lo))
    ax2.set_ylabel('Flow depth (m)')
    ft.date_axis(ax2, YEAR, every=2)
    ax2.set_xlabel(str(YEAR))
    for ax in (ax1, ax2):
        ax.legend(loc='upper right', fontsize=ft.LEGEND_SIZE, handlelength=2.6, borderaxespad=0.5)
    ax1.text(0.03, 0.95, '(a)', transform=ax1.transAxes, fontweight='bold', va='top', ha='left')
    fig.align_ylabels([ax1, ax2])
    fig.tight_layout(h_pad=0.4)
    pdfs.append(ft.save(fig, args.out or HERE, 'ex3_inflow'))

    # (b) lateral outflows
    fig, ax = ft.box_figure()
    lat = laterals()
    for k, (name, (t, qo)) in enumerate(lat.items()):
        ft.plot_input(ax, ft.to_dates(t, YEAR), qo, k, name)
    ax.set_ylim(0, 1.45 * max(v[1].max() for v in lat.values()))      # room for the legend above the peaks
    ax.set_ylabel('Lateral outflow (m³/s)')
    ft.date_axis(ax, YEAR, every=2)
    ax.set_xlabel(str(YEAR))
    ft.finish_box(fig, ax, 'b')
    pdfs.append(ft.save(fig, args.out or HERE, 'ex3_lateral'))

    if args.manuscript:
        ft.copy_to_manuscript(pdfs)


if __name__ == '__main__':
    main()
