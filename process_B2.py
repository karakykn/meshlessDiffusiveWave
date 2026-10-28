"""
process_B2.py - automated B2 sensitivity sweep (reviewer Comment 2): Crank-Nicolson weight theta and
Courant number Cr.

Started as a copy of process_B1.py and works the same way. It runs the single channel (syntheticChannel)
as it is - its own 21 uniform nodes (dx = 1 km) and the c / r_av of its shape file (8) - with the final
model in source.py (Gamma from the depth, advective outflow row) for every pair of

      Crank-Nicolson weight  theta : 0, 0.15, 0.35, 0.5, 0.65, 0.85, 1   (0 explicit, 0.5 trapezoidal, 1 implicit)
    x Courant number         Cr    : 0.25, 0.5, 0.75, 1

and scores every run with the error measures of the paper (Error Calculations), from the downstream hydrograph
against HEC-RAS (MBE against the inflow), with the same functions as process_B1.py:
    MPE    mean percentage error on the HEC-RAS times
    MBE    (V_in - V_out) / V_in over 0-24 h
    dQp    difference in peak discharge (%), peak from a parabola through the highest sample and its neighbours
    dt_a   difference in arrival time (min); the definition is process_B1.ARRIVAL
The cases use theta = 0.5 and Cr = 1 (the submitted paper used theta = 0.85 and Cr = 0.99).
Cr = 0 is not run: the time step is Cr * min(dx / C), so it would be zero and the run would never advance.

Courant number: source.py sets dt = Cr * min_i(dx_i / C_i) at every time step (C the corrected celerity), so
Cr is the largest local Courant number. It is the 'CFL' line of the input file.

Crank-Nicolson weight: theta is not read from the input file; it is the argument of Network.solve(theta),
which passes it to SingleChannel.solve_discharge at every time step. The worker counts those calls and checks the
theta they receive; a run in which solve() did not pass the swept theta at every time step is recorded as an
error (this guards against a later change in source.py).

For every run the script
  1. builds a separate case folder under OUT_DIR/runs/, a copy of syntheticChannel (which is only read) with
     Cr in the input file,
  2. runs it with the solver in source.py, several runs at a time, each in its own process with a
     time limit. A run that blows up (theta < 0.5 can) is stopped at once and recorded as unstable,
  3. scores the downstream hydrograph with the four measures above,
  4. writes OUT_DIR/B2_results.csv (one row per run) and B2_figure.png and .pdf, one panel per measure in
     FIG_METRICS (the paper: MPE only) against theta, one line per Cr; failed runs are marked with a cross, and
     runs that finished far off (outside PLOT_RANGE) with a triangle and their value. The console shows the four
     as theta x Cr tables.

Finished runs are skipped when the script is started again; a run is repeated when its settings change.
The whole sweep takes well under a minute.

Usage (from this folder, or press Run in PyCharm):
    python process_B2.py                  # run what is missing, then score and plot
    python process_B2.py --jobs 6         # runs at the same time (default: CPU count - 1)
    python process_B2.py --match _Cr1     # only runs whose name contains _Cr1
    python process_B2.py --score-only     # rescore and replot finished runs, run nothing
    python process_B2.py --plot-only      # only redraw the figure from B2_results.csv
    python process_B2.py --dry-run        # list the runs and exit
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

# figure style and helpers shared with the B1 sweep
import hashlib
from process_B1 import N_RAMP, N_STYLES, PAPER_STYLE, errors, mbe_paper, peak, trapezoid
from paper_style import savefig_safe  # noqa: E402  (PDF fallback when TrueType embedding fails)

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------------
caseName = 'syntheticChannel'                       # base case: single segment, 21 uniform nodes; only read
OUT_DIR = os.path.join('sensitivity', 'B2_final_model')   # everything the sweep writes goes here

THETAS = [0, 0.15, 0.35, 0.5, 0.65, 0.85, 1]        # Crank-Nicolson weight
COURANTS = [0.25, 0.5, 0.75, 1.0]                   # Courant number (the 'CFL' line of the input file)


def case_shape_coef():
    """c / r_av the model uses for the base case: its shape file (select_shape_coefficients), else the input file"""
    f = os.path.join(HERE, caseName, 'shape')
    if os.path.exists(f):
        return float(np.atleast_1d(np.loadtxt(f))[0])
    with open(os.path.join(HERE, caseName, 'input')) as fh:
        return float(fh.read().splitlines()[7])


SHAPE_COEF = case_shape_coef()                      # c / r_av (8 for syntheticChannel)
CASE_THETA = 0.5                                    # theta of the cases: grey line in the figure

SIM_END = 86400.0            # s; as in B1: covers the whole 24 h inflow and HEC-RAS record
RBF_TYPE = 1                 # 1: MQ

RUN_TIMEOUT = 3600           # s; a run still going after this is stopped and recorded as 'timeout'
SERIES_ROWS = 5000           # at most this many time levels are saved in each run's series.csv
N_JOBS = max(1, (os.cpu_count() or 2) - 1)

HEC_FILE = os.path.join(caseName, 'data', 'ex1_hecras', 'ex1_hec_n80.txt')   # the file ex1_Q.py uses
# version of the model: runs made with another source.py are redone
SOURCE_MD5 = hashlib.md5(open(os.path.join(HERE, 'source.py'), 'rb').read()).hexdigest()
BLOWUP_FACTOR = 10           # |Q| above this multiple of the peak inflow, at any node, counts as unstable
MAX_STEPS = 100000           # a run needing more time steps than this is stopped as unstable (normal: ~130 / Cr)
DT_JUMP = 20                 # a last time step this many times the median one (NaN or zero celerity): unstable
# Figure: a run that finished but lies outside this range is not drawn on its line; it is marked at the top of the
# panel with a triangle and its value, so it does not squash the other runs. None: draw every run.
PLOT_RANGE = {'mpe_hec': (0, 1), 'mbe_pct': (-1, 1),                 # %
              'dQp_hec_pct': (-5, 5), 'dta_hec_min': (-60, 60)}       # %, min
MPE_SCALE = 'linear'         # 'linear' or 'log'
# figure: the measures drawn, one panel each, top to bottom. The paper shows the MPE only;
# ['mpe_hec', 'mbe_pct', 'dQp_hec_pct', 'dta_hec_min'] gives the full figure.
FIG_METRICS = ['mpe_hec']


# ------------------------------------------------------------------------------------------------
# Runs and their case folders
# ------------------------------------------------------------------------------------------------
def run_name(theta, cr):
    return f'th{theta:g}_Cr{cr:g}'


def sweep_runs():
    """all (theta, Cr) pairs; Cr <= 0 is left out, since the time step would be zero"""
    return [{'name': run_name(th, cr), 'theta': th, 'Cr': cr} for th in THETAS for cr in COURANTS if cr > 0]


def run_dir(run):
    return os.path.join(OUT_DIR, 'runs', run['name'])


def settings_of(run):
    """everything that defines a run: a finished run is reused only while this is unchanged"""
    return {'base_case': caseName, 'theta': run['theta'], 'courant': run['Cr'], 'c_rav': SHAPE_COEF,
            'sim_end': SIM_END, 'rbf_type': RBF_TYPE, 'model': SOURCE_MD5}


def build_case(run):
    """case folder for one run: the base case unchanged, except Cr, c / r_av and the run settings"""
    d = run_dir(run)
    seg_src, seg_dst = os.path.join(caseName, 'segment0'), os.path.join(d, 'segment0')
    if os.path.isdir(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join(seg_src, 'geo'), os.path.join(seg_dst, 'geo'),
                    ignore=shutil.ignore_patterns('.DS_Store'))
    os.makedirs(os.path.join(seg_dst, 'run', '0'))
    for f in ('Q', 'h'):                                # initial condition
        shutil.copy(os.path.join(seg_src, 'run', '0', f), os.path.join(seg_dst, 'run', '0', f))
    shutil.copy(os.path.join(caseName, 'network'), d)

    # input: the base file with the swept values in the lines Network.__init__ reads
    with open(os.path.join(caseName, 'input')) as fh:
        lines = fh.read().splitlines()
    lines[1], lines[4], lines[7] = f'{run["Cr"]:g}', f'{RBF_TYPE:d}', f'{SHAPE_COEF:g}'
    lines[10] = f'{SIM_END:g}'
    lines[13] = '1000000000'          # print step: solve() writes no time steps; the worker records them
    with open(os.path.join(d, 'input'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    x = np.loadtxt(os.path.join(seg_dst, 'geo', 'nodes'))
    dx = np.diff(x)
    params = dict(settings_of(run), N=len(x), dx_min=float(dx.min()), dx_max=float(dx.max()))
    with open(os.path.join(d, 'params.json'), 'w') as fh:
        json.dump(params, fh, indent=1)


def is_done(run):
    d = run_dir(run)
    try:
        with open(os.path.join(d, 'params.json')) as fh:
            saved = json.load(fh)
    except (OSError, ValueError):
        return False
    unchanged = all(saved.get(k) == v for k, v in settings_of(run).items())
    return unchanged and os.path.exists(os.path.join(d, 'solver_info.json'))


# ------------------------------------------------------------------------------------------------
# One run (in a child process)
# ------------------------------------------------------------------------------------------------
class Stop(Exception):
    """raised inside solve() to stop a run that blew up"""


def write_info(d, info):
    """solver_info.json is written last: it marks the run as finished"""
    with open(os.path.join(d, 'solver_info.json'), 'w') as fh:
        json.dump(info, fh, indent=1)


def save_series(d, rows):
    """time levels (t, Q and h at both ends) to series.csv, thinned to SERIES_ROWS rows; the peak of the
    downstream discharge and its neighbours are always kept"""
    s = np.asarray(rows, dtype=float).reshape(-1, 5)
    n = len(s)
    if n > SERIES_ROWS:
        keep = np.linspace(0, n - 1, SERIES_ROWS).astype(int)
        q = np.where(np.isfinite(s[:, 2]), s[:, 2], -np.inf)
        k = int(np.argmax(q))
        s = s[np.unique(np.r_[keep, max(k - 1, 0), k, min(k + 1, n - 1)])]
    pd.DataFrame(s, columns=['t', 'Q_up', 'Q_down', 'h_up', 'h_down']).to_csv(os.path.join(d, 'series.csv'), index=False)


def run_worker(d):
    """runs one case in this process; the sweep starts it in a child process (output: d/log.txt)"""
    from source import Network

    with open(os.path.join(d, 'params.json')) as fh:
        theta = float(json.load(fh)['theta'])
    info = {'status': 'ok', 'error': '', 'theta': theta}
    try:
        solver = Network(d)
    except Exception as e:        # e.g. the RBF matrices cannot be built
        traceback.print_exc()
        write_info(d, {'status': 'error', 'error': f'setup: {type(e).__name__}: {e}'})
        return
    seg = solver.segments[0]
    x = seg.geo['nodes']
    peak_in = seg.Q_boundaries[:, 1].max()

    # theta: Network.solve(theta) passes it to seg.solve_discharge(dt, theta) at every time step. The calls are
    # counted, and only those with the swept theta count.
    solve_cn = seg.solve_discharge
    cn_calls = [0]

    def solve_discharge_counted(dt, theta_passed=None):
        if theta_passed == theta:
            cn_calls[0] += 1
        return solve_cn(dt, theta_passed)

    seg.solve_discharge = solve_discharge_counted

    # Network.step() calls update_celerity_diffusivity once at the start of every time step, when the segment
    # holds the state of the previous time level. It is wrapped to record that state and to stop the run as
    # soon as it blows up.
    rows = []                                   # (t, Q_up, Q_down, h_up, h_down) at every time level
    chk = {'finite': True, 'max_abs_Q': 0.0}
    update_params = seg.update_celerity_diffusivity

    def record():
        Q, h = seg.Q, seg.h
        rows.append((solver.time, Q[0], Q[-1], h[0], h[-1]))
        if np.isfinite(Q).all() and np.isfinite(h).all():
            chk['max_abs_Q'] = max(chk['max_abs_Q'], float(np.abs(Q).max()))
        else:
            chk['finite'] = False

    def record_and_update_params():
        record()
        if not chk['finite'] or chk['max_abs_Q'] > BLOWUP_FACTOR * peak_in:
            raise Stop(f'blew up at t = {solver.time:.0f} s (max |Q| = {chk["max_abs_Q"]:.3g})')
        if len(rows) > MAX_STEPS:
            raise Stop(f'more than {MAX_STEPS} time steps by t = {solver.time:.0f} s')
        update_params()

    seg.update_celerity_diffusivity = record_and_update_params

    t0 = time.perf_counter()
    try:
        solver.initial_depths()       # depths of the initial discharges (solve() does it too: no change then)
        storage0 = trapezoid(seg.area, x)
        solver.solve(theta)
    except Exception as e:
        if not isinstance(e, Stop):
            traceback.print_exc()
        info.update(status='error', error=f'{type(e).__name__}: {e}')
    info['wall_time'] = time.perf_counter() - t0
    record()                                              # the last time level

    s = np.asarray(rows, dtype=float)
    t, q_down = s[:, 0], s[:, 2]
    t_cut = np.r_[t[t < SIM_END], SIM_END]
    v_out_T = trapezoid(np.interp(t_cut, t, q_down), t_cut) if t[-1] >= SIM_END else float('nan')
    dts = np.diff(t)
    n_steps = len(t) - 1
    info.update(t_end=float(t[-1]), n_steps=n_steps, cn_calls=cn_calls[0],
                dt_jump=float(dts[-1] / np.median(dts[:-1])) if n_steps >= 3 else 1.0,
                v_out=trapezoid(q_down, t),                # outflow volume from every time level
                v_out_T=v_out_T,                           # outflow volume over [0, SIM_END] (MBE of the paper)
                storage_change=trapezoid(seg.area, x) - storage0 if 'storage0' in locals() else float('nan'),
                all_finite=chk['finite'], max_abs_Q=chk['max_abs_Q'])
    if info['status'] == 'ok' and cn_calls[0] != n_steps:
        info.update(status='error', error=f'theta was passed in {cn_calls[0]} of {n_steps} time steps: '
                                          f'has Network.solve() in source.py stopped passing theta to solve_discharge?')
    elif blew_up(info, peak_in):
        info['status'] = 'unstable'
    save_series(d, rows)
    write_info(d, info)


def blew_up(info, peak_in):
    """the solution went non-finite or grew without bound (the worker stops such runs with Stop). A NaN or
    zero celerity makes the time step NaN, infinite or suddenly huge, and pinv then often fails with
    LinAlgError; all count as blow-ups. Setup errors and the theta check are not blow-ups."""
    error = info.get('error', '')
    if info.get('status') not in ('ok', 'error') or error.startswith(('setup', 'theta', 'worker', 'stopped')):
        return False
    return (error.startswith('Stop') or 'LinAlgError' in error or not info.get('all_finite', True)
            or info.get('max_abs_Q', 0) > BLOWUP_FACTOR * peak_in or not np.isfinite(info.get('t_end', 0))
            or info.get('dt_jump', 1) > DT_JUMP)


# ------------------------------------------------------------------------------------------------
# The sweep
# ------------------------------------------------------------------------------------------------
def execute_one(run, env):
    d = run_dir(run)
    build_case(run)
    info_path = os.path.join(d, 'solver_info.json')
    with open(os.path.join(d, 'log.txt'), 'w') as log:
        try:
            code = subprocess.run([sys.executable, os.path.abspath(__file__), '--worker', d], cwd=HERE, env=env,
                                  stdout=log, stderr=subprocess.STDOUT, timeout=RUN_TIMEOUT).returncode
        except subprocess.TimeoutExpired:
            code = None
    if not os.path.exists(info_path):
        if code is not None and code < 0:
            # stopped by a signal (Ctrl-C, killed): not recorded, so it runs again next time
            return {'status': 'interrupted'}
        if code is None:
            write_info(d, {'status': 'timeout', 'error': f'stopped after {RUN_TIMEOUT} s', 'wall_time': float(RUN_TIMEOUT)})
        else:
            write_info(d, {'status': 'error', 'error': f'worker exited with code {code}, see log.txt'})
    with open(info_path) as fh:
        return json.load(fh)


def execute(runs, jobs):
    todo = [r for r in runs if not is_done(r)]
    print(f'{len(runs) - len(todo)} of {len(runs)} runs already done; running {len(todo)}, {jobs} at a time')
    if not todo:
        return
    todo.sort(key=lambda r: r['Cr'])        # small Cr = most time steps first
    env = dict(os.environ, MPLBACKEND='Agg')
    if jobs > 1:   # one BLAS thread per run, or parallel runs compete for the cores
        env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')

    pool = ThreadPoolExecutor(max_workers=jobs)
    futures = {pool.submit(execute_one, r, env): r for r in todo}
    w = len(str(len(todo)))
    try:
        for k, fut in enumerate(as_completed(futures), 1):
            r = futures[fut]
            try:
                info = fut.result()
            except Exception as e:    # not recorded: the run is tried again next time
                info = {'status': 'error', 'error': f'sweep: {type(e).__name__}: {e}'}
            print(f'[{k:>{w}}/{len(todo)}] {r["name"]:<16} {info["status"]:<11} '
                  f'{info.get("wall_time", float("nan")):7.2f} s  {info.get("error", "")}', flush=True)
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        sys.exit('\nStopped. Finished runs are kept; the rest run next time.')
    pool.shutdown()


# ------------------------------------------------------------------------------------------------
# Scoring
# ------------------------------------------------------------------------------------------------
def load_hec():
    """HEC-RAS discharge at the downstream end (river station 0), read as in ex1_Q.py"""
    raw = pd.read_csv(HEC_FILE, sep=r'\s+', header=None, skiprows=12).to_numpy(dtype=float)
    q = raw[raw[:, 1] == 0, 2]
    return np.linspace(0, 24 * 3600, len(q)), q


def load_inflow():
    """upstream boundary discharge (t, Q) of the base case"""
    return np.atleast_2d(np.loadtxt(os.path.join(caseName, 'segment0', 'geo', 'boundary_Q')))


def score_all(runs):
    """one row per finished run: the error measures of the paper (process_B1.errors and mbe_paper)"""
    t_hec, q_hec = load_hec()
    bq = load_inflow()
    peak_in = bq[:, 1].max()
    rows = []
    for run in runs:
        if not is_done(run):
            continue
        d = run_dir(run)
        with open(os.path.join(d, 'solver_info.json')) as fh:
            info = json.load(fh)
        status = 'unstable' if blew_up(info, peak_in) else info['status']
        row = {'run': run['name'], 'theta': run['theta'], 'courant': run['Cr'], 'status': status,
               'error': info.get('error', ''), 'mpe_hec': np.nan, 'mbe_pct': np.nan,
               'dQp_hec_pct': np.nan, 'dta_hec_min': np.nan,
               'n_steps': info.get('n_steps', np.nan), 'wall_time_s': info.get('wall_time', np.nan)}
        if status == 'ok' and os.path.exists(os.path.join(d, 'series.csv')):
            s = pd.read_csv(os.path.join(d, 'series.csv'))
            t, qd = s['t'].to_numpy(), s['Q_down'].to_numpy()
            row.update(errors(t, qd, t_hec, q_hec, 'hec'))
            row['mbe_pct'] = mbe_paper(info['v_out_T'], bq) if np.isfinite(info.get('v_out_T', np.nan)) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------
# Figure and tables
# ------------------------------------------------------------------------------------------------
def cr_line_styles(crs):
    """plot keywords per Cr: the one-hue blue ramp of B1 (light: small Cr, dark: large Cr), plus a line style
    and a marker, so the lines read the same in colour and in grey scale"""
    from matplotlib.colors import LinearSegmentedColormap, to_hex
    ramp = LinearSegmentedColormap.from_list('ramp', N_RAMP)
    if len(crs) == len(N_RAMP):
        shades = N_RAMP
    else:
        shades = [to_hex(ramp(v)) for v in (np.linspace(0, 1, len(crs)) if len(crs) > 1 else [1.0])]
    out = {}
    for k, (cr, colour) in enumerate(zip(crs, shades)):
        st = N_STYLES[k % len(N_STYLES)]
        out[cr] = dict(color=colour, mec=colour, mfc=colour if st['filled'] else 'white', ls=st['ls'], marker=st['marker'])
    return out


def plot(table, path):
    """panels stacked vertically against theta (shared axis), in the style of the paper's figures, one per measure in
    FIG_METRICS: MPE against HEC-RAS, MBE, difference in peak discharge and difference in arrival time against HEC-RAS.
    A single panel is a near-square figure one journal column wide and has no panel letter.
    One line per Cr, told apart by colour lightness, line style and marker. Runs not on the lines are marked
    in a strip at the top of the panel, one row per Cr in the colour of its Cr: a cross for a failed run,
    a triangle and the value for a run that finished outside PLOT_RANGE. The grey vertical line is the
    theta of the cases (0.5), the grey horizontal line zero (no difference). Saved as PNG and PDF."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.transforms import blended_transform_factory

    if table.empty:
        return
    crs = sorted(table['courant'].unique())
    style = cr_line_styles(crs)
    all_panels = {'mpe_hec': ('MPE (%)', MPE_SCALE), 'mbe_pct': ('MBE (%)', 'linear'),
                  'dQp_hec_pct': (r'$\Delta Q_p$ (%)', 'linear'),
                  'dta_hec_min': (r'$\Delta t_a$ (min)', 'linear')}
    panels = [(col,) + all_panels[col] for col in FIG_METRICS]
    any_off = False

    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(len(panels), 1, figsize=(3.5, max(3.0, 0.6 + 1.7 * len(panels))), sharex=True,
                                 squeeze=False)
        axes = axes[:, 0]
        for j, (ax, (col, label, scale)) in enumerate(zip(axes, panels)):
            top = blended_transform_factory(ax.transData, ax.transAxes)
            lo_r, hi_r = (PLOT_RANGE or {}).get(col) or (-np.inf, np.inf)
            shown = []
            for k, cr in enumerate(crs):
                sub = table[table['courant'] == cr].sort_values('theta')
                fin = sub[(sub['status'] == 'ok') & np.isfinite(sub[col])]
                inside = fin[(fin[col] >= lo_r) & (fin[col] <= hi_r) & ((fin[col] > 0) | (scale != 'log'))]
                ax.plot(inside['theta'], inside[col], **style[cr])
                shown.append(inside[col])
                y_row = 0.94 - 0.055 * k                    # this Cr's row in the strip at the top
                bad = sub.loc[sub['status'] != 'ok', 'theta']
                if len(bad):
                    ax.plot(bad, np.full(len(bad), y_row), ls='', marker='x', ms=4, mew=1.0,
                            color=style[cr]['color'], transform=top, clip_on=False)
                for _, r in fin.drop(inside.index).iterrows():
                    any_off = True
                    ax.plot([r['theta']], [y_row], ls='', marker='^' if r[col] > 0 else 'v', ms=3.5,
                            color=style[cr]['color'], transform=top, clip_on=False)
                    ax.text(r['theta'] + 0.03, y_row, f'{r[col]:.3g}', transform=top, va='center', ha='left',
                            fontsize=7.5, color=style[cr]['color'])
            ax.set_yscale(scale)
            ax.set_ylabel(label, labelpad=2)
            ax.set_xlim(-0.05, 1.05)
            ax.set_xticks([0, 0.5, 1], ['0', '0.5', '1'])
            ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0.25, 0.75]))
            ax.axvline(CASE_THETA, color='0.7', lw=0.6, zorder=0)
            if scale == 'linear':
                ax.axhline(0, color='0.7', lw=0.6, zorder=0)
            # y range: the runs on the lines (and 0 on a linear axis), with room at the top for the strip
            v = pd.concat(shown).replace([np.inf, -np.inf], np.nan).dropna()
            v = v if scale == 'log' else pd.concat([v, pd.Series([0.0])])
            if len(v):
                lo, hi = v.min(), v.max()
                if scale == 'log':
                    span = max(np.log10(hi / lo), 0.3)
                    ax.set_ylim(lo / 10 ** (0.08 * span), hi * 10 ** (0.6 * span))
                else:
                    span = max(hi - lo, 1e-3)
                    ax.set_ylim(lo - 0.08 * span, hi + 0.6 * span)
            if len(panels) > 1:
                ax.text(0.97, 0.94, f'({chr(97 + j)})', transform=ax.transAxes, ha='right', va='top', fontweight='bold')
        handles = [plt.Line2D([], [], **style[cr]) for cr in crs]
        labels = [f'$C_n$ = {cr:g}' for cr in crs]               # C_n: Courant number, as in the paper
        if (table['status'] != 'ok').any():
            handles.append(plt.Line2D([], [], ls='', marker='x', ms=4, mew=1.0, color='0.35'))
            labels.append('failed')
        if any_off:
            handles.append(plt.Line2D([], [], ls='', marker='^', ms=3.5, color='0.35'))
            labels.append('off scale (value)')
        axes[-1].set_xlabel(r'$\theta$')
        fig.align_ylabels(axes)
        ncol = 3
        rows = -(-len(labels) // ncol)
        fig.legend(handles, labels, loc='upper center', ncol=ncol, bbox_to_anchor=(0.5, 1.0), fontsize=9,
                   handlelength=2.2, columnspacing=0.8, handletextpad=0.4)     # 9 pt: subscripts of C_n >= 6 pt
        fig.tight_layout(rect=(0, 0, 1, 1 - (0.1 + 0.2 * rows) / fig.get_figheight()), h_pad=0.6)
        savefig_safe(fig, path)
        savefig_safe(fig, os.path.splitext(path)[0] + '.pdf')
        plt.close(fig)


def print_tables(table):
    """theta x Cr tables of the four measures; failed runs show their status"""
    for col, title, fmt in (('mpe_hec', 'MPE vs HEC-RAS (%)', '{:.3f}'), ('mbe_pct', 'MBE (%)', '{:+.3f}'),
                            ('dQp_hec_pct', 'Difference in peak discharge vs HEC-RAS (%)', '{:+.2f}'),
                            ('dta_hec_min', 'Difference in arrival time vs HEC-RAS (min)', '{:+.1f}')):
        cells = table.apply(lambda r: fmt.format(r[col]) if r['status'] == 'ok' and np.isfinite(r[col]) else r['status'],
                            axis=1)
        tab = pd.DataFrame({'theta': table['theta'], 'Cr': table['courant'], 'v': cells}).pivot(
            index='theta', columns='Cr', values='v')
        tab.index = [f'theta = {v:g}' for v in tab.index]
        tab.columns = [f'Cr = {v:g}' for v in tab.columns]
        print(f'\n{title}:')
        print(tab.to_string())


# ------------------------------------------------------------------------------------------------
def main():
    os.chdir(HERE)
    p = argparse.ArgumentParser(description='B2 sensitivity sweep: Crank-Nicolson weight theta x Courant number.')
    p.add_argument('--jobs', type=int, default=N_JOBS, help=f'runs at the same time (default {N_JOBS})')
    p.add_argument('--match', default='', help='only run the runs whose name contains this text')
    p.add_argument('--score-only', action='store_true', help='rescore and replot finished runs, run nothing')
    p.add_argument('--dry-run', action='store_true', help='list the runs and exit')
    p.add_argument('--plot-only', action='store_true', help='redraw the figure from B2_results.csv, rescore nothing')
    p.add_argument('--worker', help=argparse.SUPPRESS)     # internal: run one case folder
    args = p.parse_args()
    if args.worker:
        run_worker(args.worker)
        return
    if args.plot_only:
        plot(pd.read_csv(os.path.join(OUT_DIR, 'B2_results.csv')), os.path.join(OUT_DIR, 'B2_figure.png'))
        print(f'Figure redrawn from {os.path.join(OUT_DIR, "B2_results.csv")}')
        return

    segments = [s for s in os.listdir(caseName) if s.startswith('segment') and os.path.isdir(os.path.join(caseName, s))]
    if len(segments) != 1:
        sys.exit(f'{caseName} has {len(segments)} segments; this sweep is written for a single channel.')
    x = np.loadtxt(os.path.join(caseName, 'segment0', 'geo', 'nodes'))
    print(f'{caseName}: {len(x)} nodes, dx = {np.diff(x).min():g} - {np.diff(x).max():g} m, c / r_av = {SHAPE_COEF:g}')
    if any(cr <= 0 for cr in COURANTS):
        print('Cr <= 0 left out: the time step would be zero and the run would never advance.')
    runs = sweep_runs()
    selected = [r for r in runs if args.match in r['name']]
    if args.dry_run:
        for r in selected:
            print(f'{r["name"]:<16} {"done" if is_done(r) else ""}')
        print(f'{len(selected)} runs, {sum(map(is_done, selected))} done')
        return

    os.makedirs(os.path.join(OUT_DIR, 'runs'), exist_ok=True)
    if not args.score_only:
        execute(selected, max(1, args.jobs))

    table = score_all(runs)       # every finished run of the grid, whatever --match was
    if table.empty:
        print('No finished runs to score.')
        return
    table.to_csv(os.path.join(OUT_DIR, 'B2_results.csv'), index=False)
    plot(table, os.path.join(OUT_DIR, 'B2_figure.png'))

    print(f'\n{len(table)} runs scored: ' + ', '.join(f'{k} {v}' for k, v in table['status'].value_counts().items()))
    print_tables(table)
    print(f'\nResults in {os.path.join(HERE, OUT_DIR)}')


if __name__ == '__main__':
    main()
