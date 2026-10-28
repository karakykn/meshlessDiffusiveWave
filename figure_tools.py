"""Tools shared by the hydrograph figures of the backwater network, Lower Mississippi and Midwest cases
(<case>/scripts/hydrograph_figure.py; the single channel has its own syntheticChannel/scripts/plot_single.py):
running the meshless model (source.py) at the gauges, reading USGS records, and the figure style.

The style (paper_style.py) is the one of the sensitivity figures: Helvetica, boxed axes, inward ticks,
no grid, a boxed legend inside the top panel, blues for the models and grey for the reference."""
import os, sys, io, shutil, argparse, tempfile, contextlib, warnings
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = os.path.dirname(os.path.abspath(__file__))
CFS = 35.31466621266132            # ft3/s in one m3/s


# The style is shared with the sensitivity figures (paper_style.py).
from paper_style import PAPER_STYLE, HYDRO_SERIES, N_RAMP, savefig_safe  # noqa: E402
MARK_PHASE = {'mesh': 0.0, 'hec': 0.5, 'cns': 0.25}    # markers of different series are staggered
BOX = (3.5, 2.7)                   # in: one panel per file, one journal column wide, close to a box
LEGEND_SIZE = 8                    # pt; with ticks 8 pt and labels 9 pt no text is below 8 pt at one column width
# measured input series (boundary conditions): blues from dark to light, each with its own line style
INPUT_SERIES = [dict(color=N_RAMP[3], ls='-'), dict(color=N_RAMP[2], ls='--'), dict(color=N_RAMP[1], ls='-.'),
                dict(color=N_RAMP[0], ls=':')]
MANUSCRIPT_FIG = os.path.expanduser('~/Documents/Doctorate/Research/Diffusive Wave RBFCM/manuscript revision/fig')


def arguments(doc, cns=False):
    p = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    if cns:
        p.add_argument('--cns', choices=('n0125', 'beg'), default='n0125',
                       help='CNS series: rerun with n = 0.0125 (default) or the results of Beg et al.')
    p.add_argument('--run', action='store_true', help='rerun the meshless model before plotting')
    p.add_argument('--out', default=None, help='folder for the figure (default: the case scripts folder)')
    p.add_argument('--manuscript', action='store_true', help='also copy the PDFs to the fig folder of the manuscript')
    return p.parse_args()


# =============================================================================================
# Meshless model at the gauges
# =============================================================================================
def run_meshless(case, gauges, t_end=None, theta=0.5):
    """Runs the case with source.py and returns t (s) and {gauge: Q (m3/s)} at every time step.
    gauges = {name: (segment, node)}. The model runs on a temporary copy of the case, so nothing in the
    case folder (e.g. the junction boundary files written by Network) is changed."""
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    import source
    src = os.path.join(ROOT, case)

    def model_input_only(folder, names):          # no stored results, data or scripts
        if os.path.basename(folder) == 'run':
            return [n for n in names if n != '0']
        if os.path.abspath(folder) == src:
            return [n for n in names if n in ('data', 'scripts')]
        return []

    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(src, os.path.join(tmp, case), ignore=model_input_only)
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)     # 0/0 in the conveyance table at zero depth
                net = source.Network(case)
                if t_end is not None:
                    net.simEndTime = t_end
                at = list(gauges.values())
                net.initial_depths()
                T, Q = [net.time], [[net.segments[s].Q[n] for s, n in at]]
                while net.time < net.simEndTime:
                    net.step(theta)
                    T.append(net.time)
                    Q.append([net.segments[s].Q[n] for s, n in at])
        finally:
            os.chdir(cwd)
    Q = np.array(Q)
    return np.array(T), {k: Q[:, i] for i, k in enumerate(gauges)}


def meshless_series(case, gauges, rerun=False, t_end=None):
    """Meshless discharge at the gauges, kept in <case>/scripts/hydrographs_meshless.csv. The model is run
    when the file is missing, holds other gauges, or rerun is True."""
    path = os.path.join(ROOT, case, 'scripts', 'hydrographs_meshless.csv')
    header = 't_s,' + ','.join(gauges)
    if rerun or not os.path.exists(path) or open(path).readline().strip() != header:
        print(f'Running the meshless model for {case} ...')
        t, Q = run_meshless(case, gauges, t_end)
        np.savetxt(path, np.column_stack([t] + [Q[k] for k in gauges]), delimiter=',', header=header,
                   comments='', fmt='%.6f')
    a = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    return a[:, 0], {k: a[:, i + 1] for i, k in enumerate(gauges)}


# =============================================================================================
# Reference data
# =============================================================================================
def read_usgs(path, column, start):
    """Daily discharge from a USGS rdb file: t (s after start, 'YYYY-MM-DD') and Q (m3/s)."""
    rows = [line.rstrip('\n').split('\t') for line in open(path) if not line.startswith('#')]
    i_t, i_q = rows[0].index('datetime'), rows[0].index(column)
    t0 = datetime.strptime(start, '%Y-%m-%d')
    t, q = [], []
    for r in rows[2:]:                             # rows[1] holds the column formats
        try:
            ti = (datetime.strptime(r[i_t][:10], '%Y-%m-%d') - t0).total_seconds()
            qi = float(r[i_q]) / CFS
        except (ValueError, IndexError):           # missing or flagged values
            continue
        t.append(ti)
        q.append(qi)
    return np.array(t), np.array(q)


def to_dates(t, year):
    """seconds after 1 January of year -> matplotlib date numbers"""
    return mdates.date2num(datetime(year, 1, 1)) + np.asarray(t, float) / 86400.0


def date_axis(ax, year, every=1):
    """one year; month ticks, a label every `every` months"""
    ax.set_xlim(to_dates(0, year), to_dates(0, year + 1))
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1, 13, every)))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))


# =============================================================================================
# Plotting
# =============================================================================================
def figure(nrows, width, height):
    """Figure with nrows panels sharing the time axis."""
    plt.rcdefaults()
    plt.rcParams.update(PAPER_STYLE)
    fig, axes = plt.subplots(nrows, 1, figsize=(width, height), sharex=True, squeeze=False)
    return fig, list(axes[:, 0])


def plot(ax, x, y, kind, label=None, marks=12):
    """One series. kind: mesh, hec, cns, inflow or obs."""
    kw = dict(HYDRO_SERIES[kind])
    if label:
        kw['label'] = label
    x, y = np.asarray(x, float), np.asarray(y, float)
    if kw.get('marker') and kw.get('ls') != 'none':      # a few markers, staggered between series
        x0, x1 = x.min(), x.max()
        at = x0 + (np.arange(marks) + MARK_PHASE.get(kind, 0.0)) * (x1 - x0) / marks
        kw['markevery'] = sorted({int(np.abs(x - v).argmin()) for v in at})
    return ax.plot(x, y, gid=kind, **kw)[0]


def panel_label(ax, text):
    ax.text(0.012, 0.95, text, transform=ax.transAxes, fontweight='bold', va='top', ha='left')


def headroom(panel=0):
    """y-axis factor above the highest value: room for the legend in the top panel, for the label in the others"""
    return 1.42 if panel == 0 else 1.15


def finish(fig, axes):
    """One boxed legend inside the top panel (upper right) with the series of all panels."""
    entries = {}
    for ax in axes:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            entries.setdefault(lab, h)
    axes[0].legend(list(entries.values()), list(entries.keys()), loc='upper right', handlelength=2.6)
    fig.tight_layout()


def box_figure(width=BOX[0], height=BOX[1]):
    """one panel, close to a box (the size of the hydrograph figures of the paper)"""
    plt.rcdefaults()
    plt.rcParams.update(PAPER_STYLE)
    fig, ax = plt.subplots(figsize=(width, height))
    return fig, ax


def finish_box(fig, ax, letter=None, legend=True, loc='upper right'):
    """panel letter (bold, upper left) and a boxed legend inside the panel"""
    if letter:
        ax.text(0.03, 0.95, f'({letter})', transform=ax.transAxes, fontweight='bold', va='top', ha='left')
    if legend:
        ax.legend(loc=loc, fontsize=LEGEND_SIZE, handlelength=2.6, borderaxespad=0.5)
    fig.tight_layout()


def plot_input(ax, x, y, k, label):
    """k-th measured input series (boundary condition): line only, style INPUT_SERIES[k]"""
    return ax.plot(x, y, lw=1.3, label=label, **INPUT_SERIES[k])[0]


def stacked_figure(nrows, width=BOX[0], height=3.6):
    """nrows panels sharing the time axis, one journal column wide"""
    return figure(nrows, width, height)


def save(fig, folder, name='hydrographs'):
    """<folder>/<name>.pdf and .png; returns the PDF path"""
    os.makedirs(folder, exist_ok=True)
    for ext in ('pdf', 'png'):
        savefig_safe(fig, os.path.join(folder, f'{name}.{ext}'), bbox_inches='tight')
    plt.close(fig)
    print('Saved', os.path.join(folder, f'{name}.pdf'))
    return os.path.join(folder, f'{name}.pdf')


def copy_to_manuscript(pdfs, fig_dir=MANUSCRIPT_FIG):
    """copies the PDFs to the fig folder of the manuscript (same file names)"""
    for f in pdfs:
        shutil.copyfile(f, os.path.join(fig_dir, os.path.basename(f)))
        print('Copied to', os.path.join(fig_dir, os.path.basename(f)))
