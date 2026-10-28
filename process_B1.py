"""
process_B1.py - automated B1 sensitivity sweep (reviewer Comment 2).

Started as a copy of process.py. Instead of one run, it runs the whole B1 grid on the single
channel (syntheticChannel):

      shape parameter coefficient  c / r_av : 0.5 ... 5000
    x number of nodes              N        : 11, 21, 41, 81
    x spacing ratio                R        : dx_max / dx_min = 1 (uniform) ... 1000
    x adjacent spacing ratio       rho      : largest ratio between two neighbouring gaps, 1 < rho <= R

Node layout for a given (N, R, rho), see gap_sizes(): every gap is different. In log scale the gaps
form a random walk between the shortest gap (1) and the longest (R), with neighbouring gaps at most
rho apart; one gap is exactly 1, one exactly R and one pair of neighbours exactly rho apart, so R and
rho are met exactly. rho = R lets neighbouring gaps be anything between 1 and R (log-uniform); a
small rho makes the gaps change gradually along the channel. The gaps are then scaled so they add
up to the channel length. The random numbers come from LAYOUT_SEED, N, R and rho, so a layout is
the same every time, and the same layout is used for every c / r_av.
A (N, R, rho) combination is skipped when N is too small to grow from 1 to R with ratio rho.

For every combination the script
  1. builds a separate case folder under OUT_DIR/runs/ from syntheticChannel
     (syntheticChannel itself is only read, never changed),
  2. runs it with the final model in source.py, Network(case).solve(THETA) (Gamma from the depth,
     advective outflow row), several runs at a time, each in its own process with a time limit; the
     state at the channel ends is recorded in memory at every time step (nothing is written per time step),
  3. scores the downstream hydrograph against HEC-RAS with the error measures of the paper (Error
     Calculations): MPE on the HEC-RAS times, MBE, difference in peak discharge dQp and in arrival time dt_a
     (definition set by ARRIVAL),
  4. writes OUT_DIR/B1_results.csv (one row per run), B1_summary.csv (one row per
     R, N, c / r_av with the rho values pooled) and a figure in the style of the paper,
     B1_quicklook.png and .pdf (rows = R, columns = FIG_METRICS; colour that also reads in black and white).
Every run also records c_selected, the c / r_av that select_shape_coefficients() (source.py) picks for
its node layout, i.e. the value the model would use for that channel.

Finished runs are skipped when the script is started again, so a sweep can be stopped and resumed.
A run is repeated when its settings change; delete its folder to force a rerun.

Usage (from this folder, or press Run in PyCharm):
    python process_B1.py                  # run what is missing, then score and plot
    python process_B1.py --jobs 6         # runs at the same time (default: CPU count - 1)
    python process_B1.py --match _R100_   # only runs whose name contains _R100_
    python process_B1.py --score-only     # rescore and replot finished runs, run nothing
    python process_B1.py --plot-only      # only redraw the figure from B1_results.csv
    python process_B1.py --dry-run        # list the runs (and the skipped node layouts) and exit

Cost: the time step follows the shortest gap. With N = 81 and R = 1000 the shortest gap is about
1.4 m and a run takes roughly 10^5 time steps; R = 1000 is about 80 % of the cost of the sweep.

Mass balance error (MBE): (V_in - V_out) / V_in, with the inflow volume V_in and the outflow volume V_out
over [0, SIM_END] (24 h), as in the paper. The flow starts and ends at the base flow.
"""
import argparse
import hashlib
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

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------------
caseName = 'syntheticChannel'                  # base case (single segment); only read
OUT_DIR = os.path.join('sensitivity', 'B1_final_model')   # everything the sweep writes goes here

SHAPE_COEFS = [0.5, 1, 2, 4, 8, 16, 24, 32, 64, 128, 512, 4096]   # c / r_av
NODE_COUNTS = [11, 21, 41, 81]
SPACING_RATIOS = [1, 10, 100, 1000, 5000]    # R = dx_max / dx_min; 1 = uniform
ADJACENT_RATIOS = [1.5, 2.5, 6]                  # rho values used below R; rho = R (abrupt) is always added

CFL = 1.0                    # Courant number of the cases
THETA = 0.5                  # Crank-Nicolson weight of the cases, passed to Network.solve()
BASELINE_C = 8               # c / r_av the model uses for the single channel (syntheticChannel/shape): grey line
SIM_END = 86400.0            # s; the case input has 84600, 86400 covers the whole 24 h inflow and HEC-RAS record
RBF_TYPE = 1                 # 1: MQ
LAYOUT_SEED = 2            # change to draw a different set of node layouts
NODE_LAYOUT = 'random-walk-v1'   # stored with every run: runs made with another layout are redone

RUN_TIMEOUT = 3600           # s; a run still going after this is stopped and recorded as 'timeout'
SERIES_ROWS = 5000           # at most this many time levels are saved in each run's series.csv
N_JOBS = max(1, (os.cpu_count() or 2) - 1)

# version of the model: runs made with another source.py are redone
SOURCE_MD5 = hashlib.md5(open(os.path.join(HERE, 'source.py'), 'rb').read()).hexdigest()
HEC_FILE = os.path.join(caseName, 'data', 'ex1_hecras', 'ex1_hec_n80.txt')   # the file ex1_Q.py uses
BLOWUP_FACTOR = 10           # |Q| above this multiple of the peak inflow, at any node, counts as unstable
# figure: fixed y range and tick spacing of a column (values outside run off the panel); MPE is log-scaled and fitted
Y_LIMITS = {'mbe_pct': ((-0.1, 0.1), 0.1),    # %
            'dQp_hec_pct': ((-2, 2), 2),      # %
            'dta_hec_min': ((-10, 10), 10)}   # min
# figure: the error measures drawn, one column each in this order, and whether a last column shows the failed runs.
# The paper shows the MPE only; ['mpe_hec', 'mbe_pct', 'dQp_hec_pct', 'dta_hec_min'] with FIG_FAILED_COLUMN = True
# gives the full figure.
FIG_METRICS = ['mpe_hec']
FIG_FAILED_COLUMN = False
FIG_COLUMNS = {'mpe_hec': ('MPE (%)', 'log'), 'mbe_pct': ('MBE (%)', 'linear'),     # y label, y scale
               'dQp_hec_pct': (r'$\Delta Q_p$ (%)', 'linear'), 'dta_hec_min': (r'$\Delta t_a$ (min)', 'linear')}
ARRIVAL = 'half-rise'        # arrival-time definition (see arrival_times()): 'half-rise' or 'rise99'


# ------------------------------------------------------------------------------------------------
# Node layout
# ------------------------------------------------------------------------------------------------
def grading_steps(R, rho):
    """number of gaps needed to grow from the shortest gap to R times its length, each gap at most
    rho times its neighbour"""
    return 0 if R == 1 else int(np.ceil(np.log(R) / np.log(rho) - 1e-9))


def feasible(N, R, rho):
    """True when the N - 1 gaps leave room to grow from the shortest gap to the longest with ratio rho"""
    return R == 1 or grading_steps(R, rho) <= N - 2


def gap_sizes(N, R, rho, seed):
    """N - 1 relative gap sizes with max / min = R and largest ratio between neighbours = rho, both exact.

    Works with y = ln(gap), which must stay in [0, ln R] and change by at most ln(rho) from one gap
    to the next:
      1. two gaps are drawn at random, at least grading_steps(R, rho) gaps apart; one is the
         shortest gap (y = 0), the other the longest (y = ln R),
      2. the neighbour of the shortest gap, on the side of the longest, is set to rho (y = ln rho),
         so at least one pair of neighbours has the ratio rho exactly,
      3. every other gap is drawn, from the first to the last, uniformly in y within the range
         that keeps it within ln(rho) of the previous gap and within reach of the fixed gaps.
    Large rho gives gaps that are close to independent and log-uniform between 1 and R; small rho
    gives a random walk that wanders gradually between the two."""
    M = N - 1
    if R == 1:
        return np.ones(M)
    Y, s, k = np.log(R), np.log(rho), grading_steps(R, rho)
    rng = np.random.default_rng([seed, N, int(round(R * 1000)), int(round(rho * 1000))])
    pairs = [(a, b) for a in range(M) for b in range(M) if abs(a - b) >= k]
    i_min, i_max = pairs[rng.integers(len(pairs))]
    fixed = {i_min: 0.0, i_max: Y, i_min + int(np.sign(i_max - i_min)): s}
    # range each gap may take, given the fixed gaps and the limit on neighbouring gaps
    i = np.arange(M)
    lo, hi = np.zeros(M), np.full(M, Y)
    for p, v in fixed.items():
        lo = np.maximum(lo, v - s * np.abs(i - p))
        hi = np.minimum(hi, v + s * np.abs(i - p))
    y = np.empty(M)
    for n in range(M):
        a, b = (lo[n], hi[n]) if n == 0 else (max(lo[n], y[n - 1] - s), min(hi[n], y[n - 1] + s))
        y[n] = fixed.get(n, rng.uniform(a, b))
    return np.exp(y)


def make_nodes(x_start, x_end, N, R, rho, seed=0):
    """N nodes on [x_start, x_end] with dx_max / dx_min = R and largest adjacent gap ratio = rho:
    the gaps from gap_sizes(), scaled to add up to x_end - x_start (which keeps both ratios)"""
    g = gap_sizes(N, R, rho, seed)
    x = x_start + (x_end - x_start) * np.concatenate(([0.0], np.cumsum(g) / g.sum()))
    x[-1] = x_end
    return x


def irregularity(x):
    """R = dx_max / dx_min (the RBF mesh ratio in 1D) and rho = the largest ratio of adjacent spacings"""
    dx = np.diff(x)
    adj = dx[1:] / dx[:-1]
    rho = float(np.max(np.maximum(adj, 1 / adj))) if len(adj) else 1.0
    return float(dx.max() / dx.min()), rho


# ------------------------------------------------------------------------------------------------
# Runs and their case folders
# ------------------------------------------------------------------------------------------------
def run_name(N, R, rho, c):
    return f'N{N:03d}_R{R:g}_rho{rho:g}_c{c:g}'


def rho_values(R):
    return [1.0] if R == 1 else sorted({float(r) for r in ADJACENT_RATIOS if 1 < r < R} | {float(R)})


def sweep_runs():
    """all runs of the grid, and the (N, R, rho) layouts left out because N is too small for them"""
    runs, skipped = [], []
    for R in SPACING_RATIOS:
        for N in NODE_COUNTS:
            for rho in rho_values(R):
                if not feasible(N, R, rho):
                    skipped.append((N, R, rho))
                    continue
                for c in SHAPE_COEFS:
                    runs.append({'name': run_name(N, R, rho, c), 'N': N, 'R': R, 'rho': rho, 'c': c})
    return runs, skipped


def run_dir(run):
    return os.path.join(OUT_DIR, 'runs', run['name'])


def settings_of(run):
    """everything that defines a run: a finished run is reused only while this is unchanged"""
    return {'base_case': caseName, 'N': run['N'], 'R': run['R'], 'rho': run['rho'], 'c_rav': run['c'],
            'CFL': CFL, 'theta': THETA, 'sim_end': SIM_END, 'rbf_type': RBF_TYPE, 'node_layout': NODE_LAYOUT,
            'layout_seed': LAYOUT_SEED, 'model': SOURCE_MD5}


def build_case(run):
    """case folder for one run: the base case with new nodes, c / r_av and run settings"""
    d = run_dir(run)
    seg_src, seg_dst = os.path.join(caseName, 'segment0'), os.path.join(d, 'segment0')
    if os.path.isdir(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join(seg_src, 'geo'), os.path.join(seg_dst, 'geo'),
                    ignore=shutil.ignore_patterns('.DS_Store'))
    os.makedirs(os.path.join(seg_dst, 'run', '0'))
    shutil.copy(os.path.join(caseName, 'network'), d)

    # input: the base file with the swept values in the lines Network.__init__ reads
    with open(os.path.join(caseName, 'input')) as fh:
        lines = fh.read().splitlines()
    lines[1], lines[4], lines[7] = f'{CFL:g}', f'{RBF_TYPE:d}', f'{run["c"]:g}'
    lines[10] = f'{SIM_END:g}'
    lines[13] = '1000000000'          # print step: solve() writes no time steps; the worker records them
    with open(os.path.join(d, 'input'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    # node-dependent files: the new nodes, everything else interpolated from the base case
    geo = os.path.join(seg_dst, 'geo')
    x_old = np.loadtxt(os.path.join(geo, 'nodes'))
    x = make_nodes(x_old[0], x_old[-1], run['N'], run['R'], run['rho'], LAYOUT_SEED)
    for f in ('slopes', 'mannings_n'):
        np.savetxt(os.path.join(geo, f), np.interp(x, x_old, np.loadtxt(os.path.join(geo, f))))
    xs_old = np.loadtxt(os.path.join(geo, 'xsInfo'), dtype=int)
    nearest = np.abs(x[:, None] - x_old[None, :]).argmin(axis=1)
    np.savetxt(os.path.join(geo, 'xsInfo'), xs_old[nearest], fmt='%d')
    np.savetxt(os.path.join(geo, 'nodes'), x)
    for f in ('Q', 'h'):                                # initial condition
        v = np.loadtxt(os.path.join(seg_src, 'run', '0', f))
        np.savetxt(os.path.join(seg_dst, 'run', '0', f), np.interp(x, x_old, v))

    R, rho = irregularity(x)
    assert np.isclose(R, run['R']) and np.isclose(rho, run['rho']), f'layout gives R = {R}, rho = {rho}'
    dx = np.diff(x)
    params = dict(settings_of(run), R_nodes=R, rho_nodes=rho, dx_min=float(dx.min()), dx_max=float(dx.max()),
                  c_m=run['c'] * (x[-1] - x[0]) / (run['N'] - 1))      # c = (c / r_av) * r_av, r_av = L / (N - 1) as in buildMQ
    with open(os.path.join(d, 'params.json'), 'w') as fh:
        json.dump(params, fh, indent=1)


def selected_c(d):
    """c / r_av that select_shape_coefficients (source.py) picks for this node layout (no shape file is written,
    so the run itself uses the swept c / r_av)"""
    import contextlib, io
    from source import select_shape_coefficients
    try:
        with contextlib.redirect_stdout(io.StringIO()), np.errstate(all='ignore'):
            return float(select_shape_coefficients(d, write=False)[0][0])
    except Exception:
        return float('nan')


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
def trapezoid(y, x):
    y, x = np.asarray(y, float), np.asarray(x, float)
    return float(np.sum((y[1:] + y[:-1]) * np.diff(x)) / 2) if len(x) > 1 else 0.0


def safe_cond(M):
    try:
        return float(np.linalg.cond(M))
    except (np.linalg.LinAlgError, ValueError):
        return float('nan')


def write_info(d, info):
    """solver_info.json is written last: it marks the run as finished"""
    with open(os.path.join(d, 'solver_info.json'), 'w') as fh:
        json.dump(info, fh, indent=1)


def save_series(d, rows):
    """time levels (t, Q and h at both ends) to series.csv, thinned to SERIES_ROWS rows; the peak of
    the downstream discharge and its neighbours are always kept"""
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

    info = {'status': 'ok', 'error': ''}
    # c / r_av that the model would select for this layout; computed here, in the run's own process, and not in
    # build_case: the sweep builds cases in parallel threads, and linear algebra there can hang (seen on macOS)
    info['c_selected'] = selected_c(d)
    try:
        solver = Network(d)
    except Exception as e:        # e.g. the RBF matrices cannot be built for this c / r_av
        traceback.print_exc()
        write_info(d, {'status': 'error', 'error': f'setup: {type(e).__name__}: {e}', 'c_selected': info['c_selected']})
        return
    seg = solver.segments[0]
    N, x = seg.nodeNo, seg.geo['nodes']

    # Network.step() calls update_celerity_diffusivity once at the start of every time step, when the
    # segment holds the state of the previous time level. It is wrapped to record that state, the
    # smallest Gamma = Q / Qbar(h) (< 0: reverse flow) and the residual of the depth solution
    # (Newton's method in depth_profile) of the previous backward sequence.
    rows = []                                   # (t, Q_up, Q_down, h_up, h_down) at every time level
    chk = {'finite': True, 'max_abs_Q': 0.0, 'min_Q': np.inf, 'min_gamma': np.inf, 'max_depth_residual': 0.0,
           'max_newton_iterations': 0}
    update_params = seg.update_celerity_diffusivity

    def record():
        Q, h = seg.Q, seg.h
        rows.append((solver.time, Q[0], Q[-1], h[0], h[-1]))
        if np.isfinite(Q).all() and np.isfinite(h).all():
            chk['max_abs_Q'] = max(chk['max_abs_Q'], float(np.abs(Q).max()))
            chk['min_Q'] = min(chk['min_Q'], float(Q.min()))
        else:
            chk['finite'] = False

    def record_and_update_params():
        record()
        chk['max_depth_residual'] = max(chk['max_depth_residual'], float(getattr(seg, 'depth_residual', 0.0)))
        chk['max_newton_iterations'] = max(chk['max_newton_iterations'], int(getattr(seg, 'depth_iterations', 0)))
        update_params()
        if np.isfinite(seg.Gamma[:-1]).all():
            chk['min_gamma'] = min(chk['min_gamma'], float(seg.Gamma[:-1].min()))

    seg.update_celerity_diffusivity = record_and_update_params

    # RBF interpolation matrix with the polynomial terms, [[Phi, P], [P^T, 0]]: depends on nodes and c only
    A = np.zeros((N + 3, N + 3))
    A[:N, :] = seg.f
    A[N:, :N] = seg.f[:, N:].T
    info['cond_interp'] = safe_cond(A)

    t0 = time.perf_counter()
    try:
        solver.initial_depths()       # depths of the initial discharges (solve() does it too: no change then)
        storage0 = trapezoid(seg.area, x)
        solver.solve(THETA)
    except Exception as e:
        traceback.print_exc()
        info.update(status='error', error=f'{type(e).__name__}: {e}')
    info['wall_time'] = time.perf_counter() - t0
    record()                                              # the last time level

    s = np.asarray(rows, dtype=float)
    t, q_down = s[:, 0], s[:, 2]
    t_cut = np.r_[t[t < SIM_END], SIM_END]
    v_out_T = trapezoid(np.interp(t_cut, t, q_down), t_cut) if t[-1] >= SIM_END else float('nan')
    info.update(t_end=float(t[-1]), n_steps=len(t) - 1,
                dt_min=float(np.diff(t).min()) if len(t) > 1 else float('nan'),
                v_out=trapezoid(q_down, t),                # outflow volume from every time level
                v_out_T=v_out_T,                           # outflow volume over [0, SIM_END] (MBE of the paper)
                min_Q_down=float(np.nanmin(q_down)) if np.isfinite(q_down).any() else float('nan'),
                cond_sys=safe_cond(seg.sys),               # Crank-Nicolson system of the last step
                storage_change=trapezoid(seg.area, x) - storage0 if 'storage0' in locals() else float('nan'),
                min_gamma=chk['min_gamma'], max_depth_residual=chk['max_depth_residual'],
                max_newton_iterations=chk['max_newton_iterations'],
                all_finite=chk['finite'], max_abs_Q=chk['max_abs_Q'], min_Q=chk['min_Q'])
    if blew_up(info, seg.Q_boundaries[:, 1].max()):
        info['status'] = 'unstable'
    save_series(d, rows)
    write_info(d, info)


def blew_up(info, peak_in):
    """the solution went non-finite or grew without bound. A NaN celerity makes the next time step
    huge (t_end far beyond SIM_END) and pinv then fails with LinAlgError; both count as blow-ups.
    Errors while setting the case up (before solve) are not blow-ups."""
    error = info.get('error', '')
    return (info.get('status') in ('ok', 'error') and not error.startswith('setup')
            and (not info.get('all_finite', True) or info.get('max_abs_Q', 0) > BLOWUP_FACTOR * peak_in
                 or info.get('t_end', 0) > 1.5 * SIM_END or 'LinAlgError' in error))


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


def expected_cost(run):
    """relative cost of a run: number of time steps (1 / dx_min) times the cost of a step (~N^2)"""
    x = make_nodes(0.0, 1.0, run['N'], run['R'], run['rho'], LAYOUT_SEED)
    return run['N'] ** 2 / np.diff(x).min()


def execute(runs, jobs):
    todo = [r for r in runs if not is_done(r)]
    print(f'{len(runs) - len(todo)} of {len(runs)} runs already done; running {len(todo)}, {jobs} at a time')
    if not todo:
        return
    # the most expensive runs first, so no long run is left alone at the end
    todo.sort(key=expected_cost, reverse=True)
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
            print(f'[{k:>{w}}/{len(todo)}] {r["name"]:<30} {info["status"]:<11} '
                  f'{info.get("wall_time", float("nan")):8.1f} s  {info.get("error", "")}', flush=True)
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


def peak(t, q):
    """time and value of the peak, refined by a parabola through the highest sample and its neighbours"""
    k = int(np.argmax(q))
    if 0 < k < len(q) - 1:
        a, b, c = np.polyfit(t[k - 1:k + 2] - t[k], q[k - 1:k + 2], 2)
        if a < 0:
            return t[k] - b / (2 * a), c - b * b / (4 * a)
    return t[k], q[k]


def crossing(t, q, thr, t_after):
    """first time after t_after at which q rises through thr (linear interpolation)"""
    m = t >= t_after
    tt, qq = t[m], q[m]
    j = np.where((qq[:-1] < thr) & (qq[1:] >= thr))[0]
    if len(j) == 0:
        return np.nan
    j = j[0]
    return tt[j] + (thr - qq[j]) * (tt[j + 1] - tt[j]) / (qq[j + 1] - qq[j])


def arrival(t, q, frac):
    """first time the rising limb reaches Q_b + frac (Q_p - Q_b), with the base flow Q_b (lowest value before
    the peak) and the peak Q_p of the same hydrograph (linear interpolation between time levels)"""
    k = int(np.argmax(q))
    i_b = int(np.argmin(q[:k + 1]))
    thr = q[i_b] + frac * (peak(t, q)[1] - q[i_b])
    j = i_b + np.where(q[i_b:k + 1] >= thr)[0]
    j = int(j[0]) if len(j) else k
    if j == 0:
        return t[0]
    return t[j - 1] + (thr - q[j - 1]) * (t[j] - t[j - 1]) / (q[j] - q[j - 1])


def arrival_times(t, q, t_ref, q_ref):
    """arrival times (s) of the model and of the reference, by ARRIVAL:
    'half-rise': the rising limb reaches Q_b + 0.5 (Q_p - Q_b) of the REFERENCE (one level for both);
    'rise99':    each hydrograph reaches Q_b + 0.99 (Q_p - Q_b) of ITS OWN base flow and peak"""
    if ARRIVAL == 'half-rise':
        k = int(np.argmax(q_ref))
        i_b = int(np.argmin(q_ref[:k + 1]))
        thr = q_ref[i_b] + 0.5 * (peak(t_ref, q_ref)[1] - q_ref[i_b])
        return crossing(t, q, thr, t_ref[i_b]), crossing(t_ref, q_ref, thr, t_ref[i_b])
    if ARRIVAL == 'rise99':
        return arrival(t, q, 0.99), arrival(t_ref, q_ref, 0.99)
    raise ValueError(f"ARRIVAL must be 'half-rise' or 'rise99', not {ARRIVAL!r}")


def errors(t, q, t_ref, q_ref, tag):
    """error measures of the paper (Error Calculations): MPE on the reference times, difference in peak
    discharge (parabola through the highest sample) and in arrival time (definition: ARRIVAL, see
    arrival_times()); RMSE and peak timing are kept as extras"""
    m = np.interp(t_ref, t, q)
    tp, qp = peak(t, q)
    tr, qr = peak(t_ref, q_ref)
    ta, ta_ref = arrival_times(t, q, t_ref, q_ref)
    return {f'mpe_{tag}': np.mean(np.abs((q_ref - m) / q_ref)) * 100,
            f'rmse_{tag}': np.sqrt(np.mean((m - q_ref) ** 2)),
            f'dQp_{tag}_m3s': qp - qr,
            f'dQp_{tag}_pct': (qp - qr) / qr * 100,
            f'dta_{tag}_min': (ta - ta_ref) / 60,
            f'peakT_err_{tag}_min': (tp - tr) / 60}


def mass_balance_error(t_end, v_out, bq, storage_change):
    """(V_in - V_out - dS) / V_in in %, over the simulated window [0, t_end]"""
    tb = np.union1d(bq[(bq[:, 0] > 0) & (bq[:, 0] < t_end), 0], [0.0, t_end])
    v_in = trapezoid(np.interp(tb, bq[:, 0], bq[:, 1]), tb)      # exact for the piecewise-linear inflow
    return (v_in - v_out - storage_change) / v_in * 100


def mbe_paper(v_out_T, bq):
    """MBE of the paper in %: (V_in - V_out) / V_in over [0, SIM_END]"""
    tb = np.union1d(bq[(bq[:, 0] > 0) & (bq[:, 0] < SIM_END), 0], [0.0, SIM_END])
    v_in = trapezoid(np.interp(tb, bq[:, 0], bq[:, 1]), tb)      # exact for the piecewise-linear inflow
    return (v_in - v_out_T) / v_in * 100


def score_all(runs):
    """one row per finished run"""
    t_hec, q_hec = load_hec()
    bq = np.atleast_2d(np.loadtxt(os.path.join(caseName, 'segment0', 'geo', 'boundary_Q')))
    base_flow, peak_in = bq[0, 1], bq[:, 1].max()
    rows = []
    for run in runs:
        if not is_done(run):
            continue
        d = run_dir(run)
        with open(os.path.join(d, 'params.json')) as fh:
            params = json.load(fh)
        with open(os.path.join(d, 'solver_info.json')) as fh:
            info = json.load(fh)

        status, error = info['status'], info.get('error', '')
        if blew_up(info, peak_in):
            status = 'unstable'
        with np.errstate(divide='ignore', invalid='ignore'):
            lc_interp, lc_sys = np.log10(info.get('cond_interp', np.nan)), np.log10(info.get('cond_sys', np.nan))
        row = {'run': run['name'], 'N': run['N'], 'R': run['R'], 'rho': run['rho'], 'c_rav': run['c'], 'c_m': params['c_m'], 'dx_min': params['dx_min'], 'dx_max': params['dx_max'],
               'status': status, 'error': error,
               'wall_time_s': info.get('wall_time', np.nan), 'n_steps': info.get('n_steps', np.nan),
               'log10_cond_interp': lc_interp, 'log10_cond_sys': lc_sys,
               'c_selected': info.get('c_selected', params.get('c_selected', np.nan)), 'min_gamma': info.get('min_gamma', np.nan),
               'max_depth_residual_m': info.get('max_depth_residual', np.nan),
               'max_newton_iterations': info.get('max_newton_iterations', np.nan),
               'min_Q_all_nodes': info.get('min_Q', np.nan)}
        series = os.path.join(d, 'series.csv')
        if status == 'ok' and os.path.exists(series):
            s = pd.read_csv(series)
            t, qd = s['t'].to_numpy(), s['Q_down'].to_numpy()
            tp, qp = peak(t, qd)
            row.update(dt_min_s=info['dt_min'], dt_mean_s=info['t_end'] / max(info['n_steps'], 1),
                       Q_peak=qp, t_peak_h=tp / 3600, undershoot_cms=max(0.0, base_flow - info['min_Q_down']),
                       mbe_pct=mbe_paper(info['v_out_T'], bq) if np.isfinite(info.get('v_out_T', np.nan))
                       else np.nan)
            row.update(errors(t, qd, t_hec, q_hec, 'hec'))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(table):
    """one row per (R, N, c / r_av), the rho values pooled"""
    keys = ['R', 'N', 'c_rav']
    stats = {'mpe_hec': ['median', 'max'], 'rmse_hec': ['median'],
             'mbe_pct': ['median', 'min', 'max'], 'peakQ_err_hec_pct': ['median'], 'peakT_err_hec_min': ['median'],
             'undershoot_cms': ['max'], 'max_depth_residual_m': ['max'], 'c_selected': ['median'],
             'log10_cond_interp': ['median'],
             'dx_min': ['min'], 'n_steps': ['median'], 'wall_time_s': ['median']}
    stats = {k: v for k, v in stats.items() if k in table}
    summary = table.groupby(keys).agg(stats)
    summary.columns = [f'{a}_{b}' for a, b in summary.columns]
    counts = table.groupby(keys).agg(runs=('status', 'size'), failed=('status', lambda s: int((s != 'ok').sum())),
                                     rho_values=('rho', lambda s: ' '.join(f'{v:g}' for v in sorted(s))))
    return counts.join(summary).reset_index()


# Lines for the N values. Each N differs in three ways, so the figure reads the same in colour, in
# grey scale and for colour-blind readers: a blue that gets darker with N, a line style and a marker.
# The colours, styles and PAPER_STYLE are shared by all figures of the paper (paper_style.py).
from paper_style import N_RAMP, N_STYLES, PAPER_STYLE  # noqa: E402


def n_line_styles(Ns):
    """plot keywords per N: colour from N_RAMP (resampled when there are not four N values),
    line style and marker from N_STYLES"""
    from matplotlib.colors import LinearSegmentedColormap, to_hex
    ramp = LinearSegmentedColormap.from_list('n_ramp', N_RAMP)
    shades = N_RAMP if len(Ns) == len(N_RAMP) else [to_hex(ramp(v)) for v in np.linspace(0, 1, len(Ns))]
    out = {}
    for k, (N, colour) in enumerate(zip(Ns, shades)):
        st = N_STYLES[k % len(N_STYLES)]
        out[N] = dict(color=colour, mec=colour, mfc=colour if st['filled'] else 'white', ls=st['ls'], marker=st['marker'])
    return out


def plot(table, path):
    """figure in the style of the paper's figures, in colour that also reads in black and white.

    One row per R (top: uniform nodes, down: larger R) and one column per error measure in FIG_METRICS (MPE, MBE,
    dQp, dt_a, all against HEC-RAS) against c / r_av, plus, with FIG_FAILED_COLUMN, a last column with the failed
    runs. The rows share their y axes, so a column reads from top to bottom as the effect of R. One line per N
    (median over the rho values of the runs that finished), told apart by colour lightness, line style and marker
    (see N_STYLES). A median outside the fixed range of its column (Y_LIMITS) is not drawn; a triangle on the upper
    or lower edge of the panel marks it (pointing up: above the range, down: below). In the failed-runs column a
    square marks an N and c / r_av where runs failed: filled when every rho value failed, open when only some did.
    One or two columns fit one journal column (3.5 in wide), more take the page width. Saved as PNG and PDF."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.transforms import blended_transform_factory

    if 'mpe_hec' not in table:
        return
    Ns, R_list = sorted(table['N'].unique()), sorted(table['R'].unique())
    style = n_line_styles(Ns)
    columns = [(col,) + FIG_COLUMNS[col] + (lambda v: v,) for col in FIG_METRICS]
    ok_all = table[table['status'] == 'ok']
    nc = len(columns) + int(FIG_FAILED_COLUMN)
    narrow = nc <= 2                                     # fits one journal column
    any_off = False

    with plt.rc_context(PAPER_STYLE):
        size = (3.5, 1.3 * len(R_list) + 0.7) if narrow else (7.2, 1.4 * len(R_list) + 0.9)
        ratios = [1] * len(columns) + ([0.75] if FIG_FAILED_COLUMN else [])
        fig, axes = plt.subplots(len(R_list), nc, figsize=size, sharex=True, sharey='col', squeeze=False,
                                 gridspec_kw={'width_ratios': ratios})
        for i, R in enumerate(R_list):
            for k, N in enumerate(Ns):
                sub = table[(table['R'] == R) & (table['N'] == N)]
                ok = sub[sub['status'] == 'ok']
                for j, (col, _, _, f) in enumerate(columns):
                    med = f(ok[col]).groupby(ok['c_rav']).median().dropna()
                    if not len(med):
                        continue
                    ax = axes[i, j]
                    if col in Y_LIMITS:
                        (lo, hi), _ = Y_LIMITS[col]
                        edge = blended_transform_factory(ax.transData, ax.transAxes)
                        for marker, y_edge, out in (('^', 1.0, med > hi), ('v', 0.0, med < lo)):
                            if out.any():
                                any_off = True
                                ax.plot(med.index[out], np.full(int(out.sum()), y_edge), ls='', marker=marker, ms=3.5,
                                        color=style[N]['color'], mec=style[N]['color'], transform=edge, clip_on=False,
                                        zorder=6)
                        med = med.where((med >= lo) & (med <= hi))        # off-scale values break the line
                    ax.plot(med.index, med.values, **style[N])
                if FIG_FAILED_COLUMN:
                    # failed runs: filled square when all rho values failed at this c / r_av, open when some did
                    runs = sub.groupby('c_rav')['status'].agg(total='size', failed=lambda s: int((s != 'ok').sum()))
                    for fill, pick in (('black', runs['failed'] == runs['total']),
                                       ('white', (runs['failed'] > 0) & (runs['failed'] < runs['total']))):
                        c_bad = runs.index[pick]
                        if len(c_bad):
                            axes[i, -1].plot(c_bad, np.full(len(c_bad), k), ls='', marker='s', ms=3.5,
                                             mfc=style[N]['color'] if fill == 'black' else 'white', mec=style[N]['color'])
            for j in range(nc):
                ax = axes[i, j]
                ax.set_xscale('log')
                ax.axvline(BASELINE_C, color='0.7', lw=0.6, zorder=0)   # c / r_av of the single channel
                ax.text(0.96, 0.94, f'({chr(97 + (i * nc + j) % 26)})', transform=ax.transAxes,
                        ha='right', va='top', fontweight='bold')
            for j, (col, label, scale, _) in enumerate(columns):
                axes[i, j].set_yscale(scale)
                axes[i, j].set_ylabel(label)
                if col != 'mpe_hec':
                    axes[i, j].axhline(0, color='0.7', lw=0.6, zorder=0)
            if FIG_FAILED_COLUMN:
                axes[i, -1].set_ylabel('N')
            axes[i, -1].annotate('uniform nodes' if R == 1 else f'R = {R:g}', xy=(1.04 if narrow else 1.08, 0.5),
                                 xycoords='axes fraction', rotation=270, ha='left', va='center', fontsize=9.5,
                                 fontweight='bold')
        # y limits: fixed ones (Y_LIMITS), else all data of the column with room at the top for the panel letter
        for j, (col, _, scale, f) in enumerate(columns):
            if col in Y_LIMITS:
                (lo, hi), step = Y_LIMITS[col]
                axes[0, j].set_ylim(lo, hi)
                axes[0, j].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(step))
                continue
            v = f(ok_all[col]).replace([np.inf, -np.inf], np.nan).dropna()
            v = v[v > 0] if scale == 'log' else v
            if len(v):
                lo, hi = v.min(), v.max()
                if scale == 'log':
                    span = np.log10(hi / lo)
                    axes[0, j].set_ylim(lo / 10 ** (0.08 * span), hi * 10 ** (0.3 * span))
                else:
                    axes[0, j].set_ylim(lo - 0.08 * (hi - lo), hi + 0.3 * (hi - lo))
        if FIG_FAILED_COLUMN:
            axes[0, -1].set_ylim(-0.7, len(Ns) - 0.3 + 0.6)
            axes[0, -1].set_yticks(range(len(Ns)), [str(N) for N in Ns])
            axes[0, -1].yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            axes[0, -1].set_title('failed runs', fontsize=9)
        # x axis: a tick at every decade, labelled at every other one when the panels are narrow
        x_lo, x_hi = table['c_rav'].min() / 1.8, table['c_rav'].max() * 1.8
        decades = range(int(np.ceil(np.log10(x_lo))), int(np.floor(np.log10(x_hi))) + 1)
        every = 1 if nc == 1 else 2
        for ax in axes[-1]:
            ax.set_xlim(x_lo, x_hi)
            ax.set_xticks([10.0 ** d for d in decades], [f'10$^{{{d}}}$' if d % every == 0 else '' for d in decades])
            ax.set_xlabel('c / r$_{av}$')
        fig.align_ylabels(axes[:, 0])
        handles = [plt.Line2D([], [], **style[N]) for N in Ns]
        labels = [f'N = {N}' for N in Ns]
        if any_off:
            handles.append(plt.Line2D([], [], ls='', marker='^', ms=3.5, color='0.35'))
            labels.append('off scale')
        if FIG_FAILED_COLUMN:
            handles += [plt.Line2D([], [], ls='', marker='s', ms=3.5, mfc='0.35', mec='0.35'),
                        plt.Line2D([], [], ls='', marker='s', ms=3.5, mfc='white', mec='0.35')]
            labels += [r'all $\rho$ failed', r'some $\rho$ failed']
        ncol = 4
        rows = -(-len(labels) // ncol)
        if narrow:
            fig.legend(handles, labels, loc='upper center', ncol=ncol, bbox_to_anchor=(0.5, 1.0),
                       handlelength=2.2, columnspacing=0.8, handletextpad=0.4)
        else:
            fig.legend(handles, labels, loc='upper center', ncol=ncol, bbox_to_anchor=(0.5, 1.0),
                       handlelength=3.0, columnspacing=1.2, handletextpad=0.5)
        fig.tight_layout(rect=(0, 0, 0.965 if not narrow else 0.95, 1 - (0.12 + 0.22 * rows) / fig.get_figheight()),
                         h_pad=0.4, w_pad=0.6)
        fig.savefig(path)
        fig.savefig(os.path.splitext(path)[0] + '.pdf')
        plt.close(fig)


# ------------------------------------------------------------------------------------------------
def main():
    os.chdir(HERE)
    p = argparse.ArgumentParser(description='B1 sensitivity sweep: c / r_av x N x spacing ratios R and rho.')
    p.add_argument('--jobs', type=int, default=N_JOBS, help=f'runs at the same time (default {N_JOBS})')
    p.add_argument('--match', default='', help='only run the runs whose name contains this text')
    p.add_argument('--score-only', action='store_true', help='rescore and replot finished runs, run nothing')
    p.add_argument('--dry-run', action='store_true', help='list the runs and exit')
    p.add_argument('--plot-only', action='store_true', help='redraw the figure from B1_results.csv, rescore nothing')
    p.add_argument('--worker', help=argparse.SUPPRESS)     # internal: run one case folder
    args = p.parse_args()
    if args.worker:
        run_worker(args.worker)
        return
    if args.plot_only:
        plot(pd.read_csv(os.path.join(OUT_DIR, 'B1_results.csv')), os.path.join(OUT_DIR, 'B1_quicklook.png'))
        print(f'Figure redrawn from {os.path.join(OUT_DIR, "B1_results.csv")}')
        return

    segments = [s for s in os.listdir(caseName) if s.startswith('segment') and os.path.isdir(os.path.join(caseName, s))]
    if len(segments) != 1:
        sys.exit(f'{caseName} has {len(segments)} segments; this sweep is written for a single channel.')
    runs, skipped = sweep_runs()
    selected = [r for r in runs if args.match in r['name']]
    if args.dry_run:
        for r in selected:
            print(f'{r["name"]:<30} {"done" if is_done(r) else ""}')
        print(f'{len(selected)} runs, {sum(map(is_done, selected))} done')
        for N, R, rho in skipped:
            print(f'skipped: N = {N}, R = {R:g}, rho = {rho:g} needs {grading_steps(R, rho) + 1} gaps '
                  f'to grow from 1 to R, N = {N} has {N - 1}')
        return

    os.makedirs(os.path.join(OUT_DIR, 'runs'), exist_ok=True)
    if not args.score_only:
        execute(selected, max(1, args.jobs))

    table = score_all(runs)       # every finished run of the grid, whatever --match was
    if table.empty:
        print('No finished runs to score.')
        return
    table.to_csv(os.path.join(OUT_DIR, 'B1_results.csv'), index=False)
    summarize(table).to_csv(os.path.join(OUT_DIR, 'B1_summary.csv'), index=False)
    plot(table, os.path.join(OUT_DIR, 'B1_quicklook.png'))

    print(f'\n{len(table)} runs scored: ' + ', '.join(f'{k} {v}' for k, v in table['status'].value_counts().items()))
    base = table[(table['N'] == 21) & (table['R'] == 1) & (table['c_rav'] == BASELINE_C) & (table['status'] == 'ok')]
    if len(base) and 'mpe_hec' in base:
        b = base.iloc[0]
        print(f'Baseline (N = 21, uniform, c/r_av = {BASELINE_C:g}): MPE {b["mpe_hec"]:.3f} %, MBE {b["mbe_pct"]:.3f} %, '
              f'dQp {b["dQp_hec_m3s"]:+.3f} m3/s ({b["dQp_hec_pct"]:+.2f} %), dt_a {b["dta_hec_min"]:+.1f} min')
    if skipped:
        print(f'{len(skipped)} (N, R, rho) layouts left out because N is too small (see --dry-run)')
    print(f'Results in {os.path.join(HERE, OUT_DIR)}')


if __name__ == '__main__':
    main()
