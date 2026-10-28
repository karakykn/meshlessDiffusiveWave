"""Midwest flood, 2008: hydrographs at Grafton, St. Charles and St. Louis.
USGS daily discharge, CNS (run by the authors with the inputs of the meshless model) and the meshless model. One file per gauge, close to a box:
ex4_1 (a, Grafton), ex4_2 (b, St. Charles) and ex4_3 (c, St. Louis).
The CNS series is data/hydros/cns_gauges.csv (cns_series.py, from the CNS run of this case).

  python hydrograph_figure.py                plot (the model runs only if its stored series is missing)
  python hydrograph_figure.py --run          rerun the meshless model first (a few minutes)
  python hydrograph_figure.py --manuscript   also copy the PDFs to the fig folder of the manuscript
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(CASE_DIR))
import figure_tools as ft

CASE = os.path.basename(CASE_DIR)
DATA = os.path.join(CASE_DIR, 'data', 'hydros')
YEAR = 2008                                  # model time 0 = 1 January 2008
GAUGES = {'Grafton': (3, 0), 'St. Charles': (1, 22), 'St. Louis': (4, 4)}      # name: (segment, node)
USGS = {'Grafton': ('grafton_mississippi', '325268_00060_00003', '05587450'),  # file, discharge column, gage number
        'St. Charles': ('stCharles_mississippi', '75947_00060_00003', '06935965'),
        'St. Louis': ('stLouis_mississippi', '75971_00060_00003', '07010000')}
CNS = os.path.join(DATA, 'cns_gauges.csv')   # t_s, one column per gauge
FILES = {'Grafton': ('a', 'ex4_1'), 'St. Charles': ('b', 'ex4_2'), 'St. Louis': ('c', 'ex4_3')}   # letter, file


def read_cns():
    a = np.genfromtxt(CNS, delimiter=',', names=True)
    return {p: (a['t_s'], a[p.replace(' ', '_').replace('.', '')]) for p in GAUGES}


def main():
    args = ft.arguments(__doc__)
    t, Q = ft.meshless_series(CASE, GAUGES, rerun=args.run)
    cns = read_cns()
    pdfs = []
    for k, (p, (letter, name)) in enumerate(FILES.items()):
        fig, ax = ft.box_figure()
        f, col, gage = USGS[p]
        t_obs, q_obs = ft.read_usgs(os.path.join(DATA, f), col, f'{YEAR}-01-01')
        t_c, q_c = cns[p]
        ft.plot(ax, ft.to_dates(t_obs, YEAR), q_obs, 'obs')
        ft.plot(ax, ft.to_dates(t_c, YEAR), q_c, 'cns', label='CNS')      # CNS run by the authors for this case
        ft.plot(ax, ft.to_dates(t, YEAR), Q[p], 'mesh')
        ax.set_ylim(0, ft.headroom(k) * max(q_obs.max(), Q[p].max(), q_c.max()))
        ax.set_ylabel('Discharge (m³/s)')
        ft.date_axis(ax, YEAR, every=2)
        ax.set_xlabel(str(YEAR))
        ft.finish_box(fig, ax, letter, legend=(k == 0))      # one legend, in (a)
        pdfs.append(ft.save(fig, args.out or HERE, name))
    if args.manuscript:
        ft.copy_to_manuscript(pdfs)


if __name__ == '__main__':
    main()
