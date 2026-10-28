"""Single channel (verification case): hydrographs at the downstream end, in the style of all paper figures.

Upstream inflow, HEC-RAS dynamic wave solution, CNS (Beg et al. 2023) and the meshless model.
Run from any folder:
    python plot_single.py          plot; the meshless model runs only if hydrographs_meshless.csv is missing
    python plot_single.py --run    rerun the meshless model (source.py, under a second) before plotting
    python plot_single.py --cns beg   CNS results of Beg et al. (n = 0.04, 4 m outlet level) instead of the
                                      CNS rerun with n = 0.0125 and normal depth (data/cns_n0125)
    python plot_single.py --manuscript   also copy ex1_Q.pdf to the fig folder of the manuscript
Writes ex1_Q.pdf and ex1_Q.png (3.5 x 2.7 in, close to a box) next to this script.
Needs numpy and matplotlib; the style comes from paper_style.py in the meshlessDiffusive_Q folder."""
import os, sys, io, shutil, argparse, tempfile, contextlib, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))       # syntheticChannel/scripts
CASE_DIR = os.path.dirname(HERE)                         # syntheticChannel
ROOT = os.path.dirname(CASE_DIR)                         # meshlessDiffusive_Q: source.py, paper_style.py
sys.path.insert(0, ROOT)
from paper_style import PAPER_STYLE, HYDRO_SERIES, savefig_safe   # noqa: E402  (shared with the sensitivity figures)
MANUSCRIPT_FIG = os.path.expanduser('~/Documents/Doctorate/Research/Diffusive Wave RBFCM/manuscript revision/fig')

CASE = os.path.basename(CASE_DIR)
DATA = os.path.join(CASE_DIR, 'data')
SERIES_FILE = os.path.join(HERE, 'hydrographs_meshless.csv')   # t (s), Q at the downstream end
T_END = 86400.0                  # 24 h (the input file stops at 84,600 s)
# CNS rerun with the settings of the paper (n = 0.0125, rating curve at the outlet); output q.txt of CNS:
# t (min), channel, Q at the 21 nodes every 10 min
CNS_N0125 = os.path.join(DATA, 'cns_n0125', 'q.txt')
THETA = 0.5                      # Crank-Nicolson weight of the cases


# ---------------------------------------------------------------------------------------------
# Meshless model
# ---------------------------------------------------------------------------------------------
def run_meshless():
    """Runs the case with source.py and returns t (s) and Q at the downstream end at every time step.
    The run uses a temporary copy of the case, so nothing in the case folder is changed."""
    import source

    def model_input_only(folder, names):          # no stored results, data or scripts
        if os.path.basename(folder) == 'run':
            return [n for n in names if n != '0']
        if os.path.abspath(folder) == CASE_DIR:
            return [n for n in names if n in ('data', 'scripts')]
        return []

    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(CASE_DIR, os.path.join(tmp, CASE), ignore=model_input_only)
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)     # 0/0 in the conveyance table at zero depth
                net = source.Network(CASE)
                net.simEndTime = T_END
                net.initial_depths()
                seg = net.segments[0]
                t, q = [net.time], [seg.Q[-1]]
                while net.time < net.simEndTime:
                    net.step(THETA)
                    t.append(net.time)
                    q.append(seg.Q[-1])
        finally:
            os.chdir(cwd)
    return np.array(t), np.array(q)


def meshless(rerun):
    """Meshless series, stored in SERIES_FILE."""
    if rerun or not os.path.exists(SERIES_FILE):
        print('Running the meshless model ...')
        t, q = run_meshless()
        np.savetxt(SERIES_FILE, np.column_stack([t, q]), delimiter=',', header='t_s,Downstream end',
                   comments='', fmt='%.6f')
    a = np.loadtxt(SERIES_FILE, delimiter=',', skiprows=1)
    return a[:, 0], a[:, 1]


# ---------------------------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------------------------
def references(cns_source='n0125'):
    # HEC-RAS: discharge at river station 0 (downstream end), 145 profiles at 10-min intervals over 24 h
    hec = np.loadtxt(os.path.join(DATA, 'ex1_hecras', 'ex1_hec_n80.txt'), skiprows=12, usecols=(1, 2))
    q_hec = hec[hec[:, 0] == 0, 1]
    t_hec = np.linspace(0, 86400, len(q_hec))
    # CNS: t (s), Q at the downstream end
    if cns_source == 'n0125':
        a = np.loadtxt(CNS_N0125)
        cns = np.column_stack([a[:, 0] * 60, a[:, -1]])          # min -> s, last node
    else:
        cns = np.loadtxt(os.path.join(DATA, 'ex1_nwm', 't_dsQ_CNX'), skiprows=1)
    # prescribed inflow: t (s), Q
    inflow = np.loadtxt(os.path.join(CASE_DIR, 'segment0', 'geo', 'boundary_Q'))
    return (t_hec, q_hec), (cns[:, 0], cns[:, 1]), (inflow[:, 0], inflow[:, 1])


# ---------------------------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------------------------
MARK_PHASE = {'mesh': 0.0, 'hec': 0.5, 'cns': 0.25}    # markers of the series are staggered


def line(ax, x, y, kind, marks=12, label=None):
    """One series in the paper style (paper_style.HYDRO_SERIES); about 12 markers per line."""
    kw = dict(HYDRO_SERIES[kind])
    if label:
        kw['label'] = label
    x, y = np.asarray(x, float), np.asarray(y, float)
    if kw.get('marker') and kw.get('ls') != 'none':
        at = x.min() + (np.arange(marks) + MARK_PHASE.get(kind, 0.0)) * (x.max() - x.min()) / marks
        kw['markevery'] = sorted({int(np.abs(x - v).argmin()) for v in at})
    ax.plot(x, y, **kw)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run', action='store_true', help='rerun the meshless model before plotting')
    p.add_argument('--manuscript', action='store_true', help='also copy ex1_Q.pdf to the fig folder of the manuscript')
    p.add_argument('--cns', choices=('n0125', 'beg'), default='n0125',
                   help='CNS series: rerun with n = 0.0125 (default) or the results of Beg et al.')
    args = p.parse_args()

    t, q = meshless(args.run)
    hec, cns, inflow = references(args.cns)
    h = 1 / 3600                                            # s -> h

    plt.rcParams.update(PAPER_STYLE)
    fig, ax = plt.subplots(figsize=(3.5, 2.7))                     # one journal column, close to a box
    line(ax, inflow[0] * h, inflow[1], 'inflow')
    line(ax, hec[0] * h, hec[1], 'hec')
    line(ax, cns[0] * h, cns[1], 'cns', label='CNS' if args.cns == 'n0125' else None)   # n0125: CNS run by the authors
    line(ax, t * h, q, 'mesh')
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 4))
    ax.set_ylim(19.5, 25.5)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Discharge (m³/s)')

    ax.legend(loc='upper right', fontsize=8, handlelength=2.6, borderaxespad=0.5)   # boxed legend inside the axes
    fig.tight_layout()

    for ext in ('pdf', 'png'):
        savefig_safe(fig, os.path.join(HERE, f'ex1_Q.{ext}'), bbox_inches='tight')
    print('Saved', os.path.join(HERE, 'ex1_Q.pdf'))
    if args.manuscript:
        shutil.copyfile(os.path.join(HERE, 'ex1_Q.pdf'), os.path.join(MANUSCRIPT_FIG, 'ex1_Q.pdf'))
        print('Copied to', os.path.join(MANUSCRIPT_FIG, 'ex1_Q.pdf'))


if __name__ == '__main__':
    main()
