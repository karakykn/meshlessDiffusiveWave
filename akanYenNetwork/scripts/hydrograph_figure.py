"""Akan and Yen (1981), Example I: outflow hydrographs of channels 1, 2 (= 3), 4 and 6.
Dynamic wave solution of Akan and Yen (1981, Fig. 8, digitized: data/akanYen_fig8_digitized.csv), CNS
(data/cns_q.txt, MESH_Code_irregular-master/CNS_timing/cases_akanYen/akan_yen_ex1) and the meshless model.
One file per panel, close to a box: ex2_1 (a, channel 1), ex2_2 (b, channels 2 and 3), ex2_3 (c, channel 6).
Channel 4 is not plotted (same channel and inflow as channel 1, nearly the same hydrograph; 4-panel version:
ex2_1..ex2_4 before 2026-10-08 evening). Also writes the error metrics against the dynamic wave solution for
channels 1, 2-3, 4 and 6 (metrics_ex2.csv).

  python hydrograph_figure.py                plot (the model runs only if its stored series is missing)
  python hydrograph_figure.py --run          rerun the meshless model first (a few seconds)
  python hydrograph_figure.py --manuscript   also copy the PDFs to the fig folder of the manuscript

Metrics (reference: Akan and Yen's dynamic wave on its digitized times, every 0.25 min over 0-40 min):
  MPE   Eq. 41, model interpolated to the reference times
  MBE   volume balance of the channel that ends at the location, (V_in - V_out - dS) / V_in over 0-40 min
        (dS: change of its storage). Not against the reference: Akan and Yen's printed curves pass 4.7 % less
        water than enters channel 1 over 0-40 min (their coarse 1981 scheme), so a volume comparison with them
        measures their mass loss. CNS storage from data/cns_area.txt.
  dQp   model peak - reference peak (parabola through the highest sample); channels 2-3: the maximum
  dt_a  arrival-time difference (min) where the rising limb crosses Q_base + 0.5 (Q_p - Q_base) of the reference;
        not defined for channels 2-3 (no rising limb from base flow)
"""
import os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(CASE_DIR))
import figure_tools as ft

CASE = os.path.basename(CASE_DIR)
DATA = os.path.join(CASE_DIR, 'data')
# name: (segment, node, panel of Fig. 8 of Akan and Yen, CNS channel)
GAUGES = {'Channel 1': (0, -1), 'Channels 2 and 3': (1, -1), 'Channel 4': (3, -1), 'Channel 6': (5, -1)}
FIG8 = {'Channel 1': 'a', 'Channels 2 and 3': 'b', 'Channel 4': 'c', 'Channel 6': 'd'}
CNS_CH = {'Channel 1': 1, 'Channels 2 and 3': 2, 'Channel 4': 4, 'Channel 6': 6}
FILES = {'Channel 1': ('a', 'ex2_1'), 'Channels 2 and 3': ('b', 'ex2_2'), 'Channel 6': ('c', 'ex2_3')}
YLIM = {'Channel 1': (2, 17), 'Channels 2 and 3': (0.8, 2.8), 'Channel 4': (2, 13), 'Channel 6': (8, 26)}


def inflow_triangle(t_min):
    return np.interp(t_min, [0, 10, 20, 1e3], [3, 12, 3, 3])


def akan_yen(panel, curve='dynamic'):
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'akanYen_fig8_digitized.csv')))
            if r['panel'] == panel and r['curve'] == curve]
    return np.array([float(r['t_min']) for r in rows]), np.array([float(r['Q_m3s']) for r in rows])


def cns(channel):
    a = np.loadtxt(os.path.join(DATA, 'cns_q.txt'))
    b = a[a[:, 1] == channel]
    return b[:, 0], b[:, -1]                        # time (min), discharge at the last node


def peak(t, q):
    i = int(np.argmax(q))
    if 0 < i < len(q) - 1:
        a, b, c = np.polyfit(t[i - 1:i + 2] - t[i], q[i - 1:i + 2], 2)
        if a < 0:
            x = -b / (2 * a)
            if abs(x) <= max(t[i] - t[i - 1], t[i + 1] - t[i]):
                return t[i] + x, c - b * b / (4 * a)
    return t[i], q[i]


def crossing(t, q, thr):
    j = np.where((q[:-1] < thr) & (q[1:] >= thr))[0]
    if len(j) == 0:
        return np.nan
    j = j[0]
    return t[j] + (thr - q[j]) * (t[j + 1] - t[j]) / (q[j + 1] - q[j])


def scores(tr, qr, tm, qm, arrival=True):
    qi = np.interp(tr, tm, qm)
    mpe = np.mean(np.abs(qr - qi) / qr) * 100
    vr = np.sum((qr[1:] + qr[:-1]) * np.diff(tr)) / 2
    vm = np.sum((qi[1:] + qi[:-1]) * np.diff(tr)) / 2
    tp_r, qp_r = peak(tr, qr)
    tp_m, qp_m = peak(tm, qm)
    dta = np.nan
    if arrival:
        thr = qr[0] + 0.5 * (qp_r - qr[0])
        dta = crossing(tm, qm, thr) - crossing(tr, qr, thr)
    return dict(MPE=mpe, MBE=(vr - vm) / vr * 100, dQp=qp_m - qp_r, dQp_pct=(qp_m - qp_r) / qp_r * 100, dta_min=dta,
                Qp_ref=qp_r, Qp_model=qp_m)


def meshless_balance():
    '''(V_in - V_out - dS) / V_in (%) of every channel over 0-40 min, from a run of the model'''
    import io, shutil, tempfile, contextlib, warnings
    import source
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(CASE_DIR, os.path.join(tmp, CASE),
                        ignore=lambda d, n: [x for x in n if x in ('data', 'scripts') or (os.path.basename(d) == 'run' and x != '0')])
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                net = source.Network(CASE)
                net.initial_depths()
                segs = net.segments
                x = [segs[i].geo['nodes'] for i in range(6)]
                stor = lambda i: float(np.sum((segs[i].area[1:] + segs[i].area[:-1]) * np.diff(x[i])) / 2)
                S0 = [stor(i) for i in range(6)]
                vin, vout = np.zeros(6), np.zeros(6)
                q0 = [(segs[i].Q[0], segs[i].Q[-1]) for i in range(6)]
                while net.time < net.simEndTime - 1e-9:
                    t_old = net.time
                    net.step(0.5)
                    dt = net.time - t_old
                    for i in range(6):
                        a, b = q0[i]
                        if net.time > 2400.0:            # last step: up to 40 min, linear in time
                            w = (2400.0 - t_old) / dt
                            qa = a + w * (segs[i].Q[0] - a)
                            qb = b + w * (segs[i].Q[-1] - b)
                            vin[i] += (a + qa) / 2 * (2400.0 - t_old)
                            vout[i] += (b + qb) / 2 * (2400.0 - t_old)
                        else:
                            vin[i] += (a + segs[i].Q[0]) / 2 * dt
                            vout[i] += (b + segs[i].Q[-1]) / 2 * dt
                        q0[i] = (segs[i].Q[0], segs[i].Q[-1])
                S1 = [stor(i) for i in range(6)]
        finally:
            os.chdir(cwd)
    return {i + 1: (vin[i] - vout[i] - (S1[i] - S0[i])) / vin[i] * 100 for i in range(6)}


def cns_balance():
    q = np.loadtxt(os.path.join(DATA, 'cns_q.txt'))
    A = np.loadtxt(os.path.join(DATA, 'cns_area.txt'))
    out = {}
    for j in range(1, 7):
        qj, Aj = q[q[:, 1] == j], A[A[:, 1] == j]
        t = qj[:, 0] * 60
        x = np.arange(Aj.shape[1] - 2) * 60.0
        trap = lambda tt, v: np.sum((v[1:] + v[:-1]) * np.diff(tt)) / 2
        vin, vout = trap(t, qj[:, 2]), trap(t, qj[:, -1])
        dS = trap(x, Aj[-1, 2:]) - trap(x, Aj[0, 2:])
        out[j] = (vin - vout - dS) / vin * 100
    return out


def main():
    args = ft.arguments(__doc__)
    t, Q = ft.meshless_series(CASE, GAUGES, rerun=args.run)
    t_min = t / 60.0
    bal_mesh, bal_cns = meshless_balance(), cns_balance()
    pdfs, rows = [], []
    for k, g in enumerate(GAUGES):
        tr, qr = akan_yen(FIG8[g])
        tc, qc = cns(CNS_CH[g])
        for model, (tm, qm) in (('CNS', (tc, qc)), ('Meshless', (t_min, Q[g]))):
            s = scores(tr, qr, tm, qm, arrival=(g != 'Channels 2 and 3'))
            s['MBE'] = (bal_cns if model == 'CNS' else bal_mesh)[CNS_CH[g]]
            rows.append([g, model] + [round(s[x], 3) for x in ('MPE', 'MBE', 'dQp', 'dQp_pct', 'dta_min',
                                                                'Qp_ref', 'Qp_model')])
        if g not in FILES:
            continue
        letter, name = FILES[g]
        fig, ax = ft.box_figure()
        ft.plot(ax, tr, qr, 'hec', label='Dynamic (Akan and Yen 1981)')
        ft.plot(ax, tc, qc, 'cns', label='CNS')      # CNS run by the authors for this case
        ft.plot(ax, t_min, Q[g], 'mesh')
        if g in ('Channel 1', 'Channel 4'):
            ft.plot(ax, tr, inflow_triangle(tr), 'inflow')
        elif g == 'Channels 2 and 3':
            ft.plot(ax, [0, 40], [2.0, 2.0], 'inflow')
        ax.set_xlim(0, 40)
        ax.set_xticks(range(0, 41, 10))
        ax.set_ylim(*YLIM[g])
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Discharge (m³/s)')
        ft.finish_box(fig, ax, letter, legend=(k == 0))      # one legend, in (a)
        pdfs.append(ft.save(fig, args.out or HERE, name))
    with open(os.path.join(args.out or HERE, 'metrics_ex2.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['location (outflow)', 'model', 'MPE %', 'MBE %', 'dQp m3/s', 'dQp %', 'dt_a min', 'Qp ref',
                    'Qp model'])
        w.writerows(rows)
    for r in rows:
        print(' | '.join(str(x) for x in r))
    if args.manuscript:
        ft.copy_to_manuscript(pdfs)


if __name__ == '__main__':
    main()
