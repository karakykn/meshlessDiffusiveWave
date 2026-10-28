"""Lower Mississippi River, 2011: hydrographs at Baton Rouge and Belle Chasse.
USGS daily discharge, CNS (Beg et al. 2023) and the meshless model. One file per gauge, close to a box:
ex3_1 (a, Baton Rouge) and ex3_2 (b, Belle Chasse).

  python hydrograph_figure.py                plot (the model runs only if its stored series is missing)
  python hydrograph_figure.py --run          rerun the meshless model first (about a minute)
  python hydrograph_figure.py --manuscript   also copy the PDFs to the fig folder of the manuscript
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(CASE_DIR))
import figure_tools as ft

CASE = os.path.basename(CASE_DIR)
DATA = os.path.join(CASE_DIR, 'data')
YEAR = 2011                                  # model time 0 = 1 January 2011
GAUGES = {'Baton Rouge': (0, 28), 'Belle Chasse': (0, 72)}            # name: (segment, node)
USGS = {'Baton Rouge': ('usgs_baton', '61176_00060_00003', '07374000'),  # file, discharge column, gage number
        'Belle Chasse': ('usgs_belle', '61182_00060_00003', '07374525')}
CNS = {'Baton Rouge': 'baton_disc_beg', 'Belle Chasse': 'belle_disc_beg'}  # t (s), Q (m3/s)
FILES = {'Baton Rouge': ('a', 'ex3_1'), 'Belle Chasse': ('b', 'ex3_2')}  # letter, file


def main():
    args = ft.arguments(__doc__)
    t, Q = ft.meshless_series(CASE, GAUGES, rerun=args.run)
    pdfs = []
    for k, (p, (letter, name)) in enumerate(FILES.items()):
        fig, ax = ft.box_figure()
        f, col, gage = USGS[p]
        t_obs, q_obs = ft.read_usgs(os.path.join(DATA, f), col, f'{YEAR}-01-01')
        cns = np.loadtxt(os.path.join(DATA, 'beg_data', CNS[p]), skiprows=1)
        ft.plot(ax, ft.to_dates(t_obs, YEAR), q_obs, 'obs')
        ft.plot(ax, ft.to_dates(cns[:, 0], YEAR), cns[:, 1], 'cns')
        ft.plot(ax, ft.to_dates(t, YEAR), Q[p], 'mesh')
        ax.set_ylim(0, ft.headroom(k) * max(q_obs.max(), Q[p].max(), cns[:, 1].max()))
        ax.set_ylabel('Discharge (m³/s)')
        ft.date_axis(ax, YEAR, every=2)
        ax.set_xlabel(str(YEAR))
        ft.finish_box(fig, ax, letter, legend=(k == 0))      # one legend, in (a)
        pdfs.append(ft.save(fig, args.out or HERE, name))
    if args.manuscript:
        ft.copy_to_manuscript(pdfs)


if __name__ == '__main__':
    main()
