"""Error metrics of the error-summary table (Table 4) for all cases, and a LaTeX draft of the table.

MPE   on the reference times (model interpolated there)
MBE   prismatic cases: (V_in - V_out) / V_in over 0-24 h, V_in exact from the prescribed inflow;
      real cases: (V_ref - V_model) / V_ref over the MPE period (model sampled at the daily gage times)
dQp   model peak - reference peak (parabola through the highest sample); real cases: in the flood event
dt_a  arrival-time difference, model - reference, by ARRIVAL (real cases: within the flood event):
      'half-rise': the rising limb reaches Q_b + 0.5 (Q_p - Q_b) of the reference (one level for both)
      'rise99':    each hydrograph reaches Q_b + 0.99 (Q_p - Q_b) of its own base flow Q_b and peak Q_p
Meshless series: <case>/scripts/hydrographs_meshless.csv (written by the case plotting scripts). The backwater
network is scored by akanYenNetwork/scripts/hydrograph_figure.py (metrics_ex2.csv, read here). CPU times
(time_meshless.py and the CNS timing runs) are not computed here; the Runtime column is left empty.
Run from anywhere: python table_metrics.py   -> prints the metrics, writes table_error_summary.tex here."""
import csv, os, sys
from datetime import datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # meshlessDiffusive_Q
sys.path.insert(0, ROOT)
from figure_tools import read_usgs                               # noqa: E402
ARRIVAL = 'half-rise'                                            # or 'rise99'
FULL = (0.0, 86400.0)


def trap(y, x):
    return float(np.sum((y[1:] + y[:-1]) * (x[1:] - x[:-1])) / 2)


def peak(t, q):
    i = int(np.argmax(q))
    if 0 < i < len(q) - 1:
        a, b, c = np.polyfit(t[i - 1:i + 2] - t[i], q[i - 1:i + 2], 2)
        if a < 0 and abs(-b / (2 * a)) <= max(t[i] - t[i - 1], t[i + 1] - t[i]):
            return t[i] - b / (2 * a), c - b * b / (4 * a)
    return t[i], q[i]


def crossing(t, q, thr, t_after):
    m = t >= t_after
    tt, qq = t[m], q[m]
    j = np.where((qq[:-1] < thr) & (qq[1:] >= thr))[0]
    if len(j) == 0:
        return np.nan
    j = j[0]
    return tt[j] + (thr - qq[j]) * (tt[j + 1] - tt[j]) / (qq[j + 1] - qq[j])


def arrival(t, q, frac):
    k = int(np.argmax(q))
    i_b = int(np.argmin(q[:k + 1]))
    thr = q[i_b] + frac * (peak(t, q)[1] - q[i_b])
    j = i_b + np.where(q[i_b:k + 1] >= thr)[0]
    j = int(j[0]) if len(j) else k
    return t[0] if j == 0 else t[j - 1] + (thr - q[j - 1]) * (t[j] - t[j - 1]) / (q[j] - q[j - 1])


def arrival_difference(tr, qr, tm, qm):
    """model - reference arrival time (h); tr, qr and tm, qm restricted to the evaluation window"""
    if ARRIVAL == 'half-rise':
        k = int(np.argmax(qr))
        i_b = int(np.argmin(qr[:k + 1]))
        thr = qr[i_b] + 0.5 * (peak(tr, qr)[1] - qr[i_b])
        return (crossing(tm, qm, thr, tr[i_b]) - crossing(tr, qr, thr, tr[i_b])) / 3600
    if ARRIVAL == 'rise99':
        return (arrival(tm, qm, 0.99) - arrival(tr, qr, 0.99)) / 3600
    raise ValueError(ARRIVAL)


def vin(*files):
    tt = np.linspace(0, FULL[1], 864001)
    return sum(trap(np.interp(tt, *np.loadtxt(f).T[:2]), tt) for f in files)


def mbe_synthetic(v_in, t, q):
    m = (t > 0) & (t < FULL[1])
    tt = np.r_[0.0, t[m], FULL[1]]
    return (v_in - trap(np.interp(tt, t, q), tt)) / v_in * 100


def scores(tr, qr, tm, qm, win, mpe_win, mbe=None):
    m = (tr >= mpe_win[0]) & (tr <= mpe_win[1])
    qi = np.interp(tr[m], tm, qm)
    mpe = np.mean(np.abs(qr[m] - qi) / qr[m]) * 100
    if mbe is None:
        mbe = (trap(qr[m], tr[m]) - trap(qi, tr[m])) / trap(qr[m], tr[m]) * 100
    w, wm = (tr >= win[0]) & (tr <= win[1]), (tm >= win[0]) & (tm <= win[1])
    qp_r, qp_m = peak(tr[w], qr[w])[1], peak(tm[wm], qm[wm])[1]
    dta = arrival_difference(tr[w], qr[w], tm[wm], qm[wm])
    return dict(MPE=mpe, MBE=mbe, dQp=qp_m - qp_r, dQp_pct=(qp_m - qp_r) / qp_r * 100, dta=dta)


def meshless(case):
    path = os.path.join(ROOT, case, 'scripts', 'hydrographs_meshless.csv')
    names = open(path).readline().strip().split(',')[1:]
    a = np.loadtxt(path, delimiter=',', skiprows=1)
    return a[:, 0], {n: a[:, i + 1] for i, n in enumerate(names)}


def day(y, m, d):
    return (datetime(y, m, d) - datetime(y, 1, 1)).total_seconds()


def all_rows():
    rows = []
    # single channel
    d = os.path.join(ROOT, 'syntheticChannel')
    hec = np.loadtxt(f'{d}/data/ex1_hecras/ex1_hec_n80.txt', skiprows=12, usecols=(1, 2))
    qr = hec[hec[:, 0] == 0, 1]; tr = np.linspace(0, 86400, len(qr))
    a = np.loadtxt(f'{d}/data/cns_n0125/q.txt')                   # CNS, n = 0.0125, rating curve: t (min), channel, Q
    cns = np.column_stack([a[:, 0] * 60, a[:, -1]])
    v1 = vin(f'{d}/segment0/geo/boundary_Q')
    t, Q = meshless('syntheticChannel')
    rows.append(('Single channel', 'CNS', 21, 'Ds. end', scores(tr, qr, cns[:, 0], cns[:, 1], FULL, FULL, mbe_synthetic(v1, cns[:, 0], cns[:, 1]))))
    rows.append(('Single channel', 'Meshless', 21, 'Ds. end', scores(tr, qr, t, Q['Downstream end'], FULL, FULL, mbe_synthetic(v1, t, Q['Downstream end']))))
    # backwater network (Akan and Yen 1981): metrics of akanYenNetwork/scripts/hydrograph_figure.py
    f = os.path.join(ROOT, 'akanYenNetwork', 'scripts', 'metrics_ex2.csv')
    rd = list(csv.DictReader(open(f)))
    loc = {'Channel 1': 'Ch. 1 ds. end', 'Channels 2 and 3': 'Ch. 2 ds. end', 'Channel 6': 'Ch. 6 ds. end'}
    for model in ('CNS', 'Meshless'):
        for r in rd:
            if r['model'] == model and r['location (outflow)'] in loc:
                rows.append(('Backwater network', model, 66, loc[r['location (outflow)']],
                             dict(MPE=float(r['MPE %']), MBE=float(r['MBE %']), dQp=float(r['dQp m3/s']),
                                  dQp_pct=float(r['dQp %']), dta=float(r['dt_a min']) / 60)))
    # Lower Mississippi, 2011
    d = os.path.join(ROOT, 'lowerMississippi_paper', 'data')
    usgs = {'Baton Rouge': read_usgs(f'{d}/usgs_baton', '61176_00060_00003', '2011-01-01'),
            'Belle Chasse': read_usgs(f'{d}/usgs_belle', '61182_00060_00003', '2011-01-01')}
    cns = {p: np.loadtxt(f'{d}/beg_data/{f}', skiprows=1) for p, f in (('Baton Rouge', 'baton_disc_beg'), ('Belle Chasse', 'belle_disc_beg'))}
    mpe_win = (0.0, cns['Baton Rouge'][-1, 0])                     # common period: end of the CNS record
    win = (day(2011, 4, 11), day(2011, 6, 30))
    t, Q = meshless('lowerMississippi_paper')                    # one channel, spillways as point lateral outflows
    for p in usgs:
        rows.append(('LMR', 'CNS', 97, p, scores(*usgs[p], cns[p][:, 0], cns[p][:, 1], win, mpe_win)))
    for p in usgs:
        rows.append(('LMR', 'Meshless', 96, p, scores(*usgs[p], t, Q[p], win, mpe_win)))
    # Midwest, 2008
    d = os.path.join(ROOT, '2008flood_stLouis', 'data', 'hydros')
    usgs = {'Grafton': read_usgs(f'{d}/grafton_mississippi', '325268_00060_00003', '2008-01-01'),
            'St. Charles': read_usgs(f'{d}/stCharles_mississippi', '75947_00060_00003', '2008-01-01'),
            'St. Louis': read_usgs(f'{d}/stLouis_mississippi', '75971_00060_00003', '2008-01-01')}
    win = (day(2008, 5, 30), day(2008, 7, 19))
    t, Q = meshless('2008flood_stLouis')
    mpe_win = (0.0, t[-1])
    f = os.path.join(d, 'cns_gauges.csv')                          # CNS run of the same case (scripts/cns_series.py)
    if os.path.exists(f):
        c = np.loadtxt(f, delimiter=',', skiprows=1)
        cns = {p: (c[:, 0], c[:, i + 1]) for i, p in enumerate(usgs)}
        mpe_win = (0.0, min(t[-1], c[-1, 0]))                     # common period
        for p in usgs:
            rows.append(('Midwest 2008', 'CNS', 109, p, scores(*usgs[p], *cns[p], win, mpe_win)))
    for p in usgs:
        rows.append(('Midwest 2008', 'Meshless', 109, p, scores(*usgs[p], t, Q[p], win, mpe_win)))
    return rows


def fmt(v, dec, sign=True):
    s = f'{v:.{dec}f}'
    if float(s) == 0:
        return f'{abs(float(s)):.{dec}f}'
    return ('+' if sign and v > 0 else '') + s


def latex(rows):
    synthetic = {'Single channel', 'Backwater network'}
    best = {}
    for case, model, n, loc, r in rows:
        for k in ('MPE', 'MBE', 'dQp', 'dQp_pct', 'dta'):
            best.setdefault((case, loc, k), []).append((abs(r[k]), model))
    def cell(case, model, loc, k, text):
        c = best[(case, loc, k)]
        win = len({m for _, m in c}) > 1 and min(c)[1] == model
        return f'\\boldmath${text}$' if win else f'${text}$'
    out, last_case, last_model = [], None, None
    for case, model, n, loc, r in rows:
        syn = case in synthetic
        if case != last_case:
            if last_case is not None:
                out.append('\\midrule')
            first = case
        elif model != last_model:
            out.append('\\cmidrule(l){2-10}')
            first = ''
        else:
            first = ''
        mod = model if (case != last_case or model != last_model) else ''
        nn = str(n) if mod else ''
        cells = [cell(case, model, loc, 'MPE', f"{r['MPE']:.2f}"),
                 cell(case, model, loc, 'MBE', fmt(r['MBE'], 2)),
                 cell(case, model, loc, 'dQp', fmt(r['dQp'], 2) if syn else fmt(r['dQp'], 0)),
                 cell(case, model, loc, 'dQp_pct', fmt(r['dQp_pct'], 1)),
                 cell(case, model, loc, 'dta', fmt(r['dta'], 3 if case == 'Backwater network' else 2) if syn
                      else fmt(r['dta'], 0))]
        out.append(f'{first} & {mod} & {nn} &  & {loc} & ' + ' & '.join(cells) + ' \\\\')
        last_case, last_model = case, model
    head = r"""\begin{table}[htbp]
\caption{Error metrics at the control points (Section~\ref{sec:errors}). $N$: total number of nodes. Differences are model minus reference. The reference is HEC-RAS for the single channel, the dynamic wave solution of Akan and Yen (1981) for the backwater network, and the USGS gage record for the LMR and Midwest 2008. Bold: better value at the location.}
\label{tab:error-summary}
\centering
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llrrlrrrrr}
\toprule
Case & Model & $N$ & Runtime & Location & MPE & MBE & \multicolumn{2}{c}{$\Delta Q_p$} & $\Delta t_a$ \\
\cmidrule(lr){8-9}
 & & & (s) & & (\%) & (\%) & (m$^3$/s) & (\%) & (h) \\
\midrule
"""
    return head + '\n'.join(out) + '\n\\bottomrule\n\\end{tabular}\n\\end{table}\n'


if __name__ == '__main__':
    rows = all_rows()
    for case, model, n, loc, r in rows:
        print(f"{case:15s} {model:9s} {loc:13s} MPE {r['MPE']:7.3f}  MBE {r['MBE']:7.3f}  dQp {r['dQp']:9.3f} ({r['dQp_pct']:6.2f} %)  dt_a {r['dta']:7.2f} h")
    open(os.path.join(HERE, 'table_error_summary.tex'), 'w').write(latex(rows))
    print('Wrote', os.path.join(HERE, 'table_error_summary.tex'))
