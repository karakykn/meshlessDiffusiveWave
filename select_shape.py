"""Shape parameter of every channel of the application cases.

For each channel, select_shape_coefficients (source.py) evaluates the scheme's operators, dQ/dx, d2Q/dx2 and the
depth integration, on two smooth test profiles over the channel's own nodes and picks the c / r_av with the
smallest error. Here the candidates are restricted to the admissible range of the B1 sweep (ADMISSIBLE), so
the test function only chooses within the range where the model is stable and accurate for every N, R and rho.

Writes <case>/shape for every case (read by Network.load_segments) and, in sensitivity/shape_selection/,
table_shape.tex (selected values) and test_errors.csv (test-function error of every channel over a wide
range of c / r_av, for reference).

    python select_shape.py
"""
import contextlib
import io
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from source import select_shape_coefficients  # noqa: E402

ADMISSIBLE = (8, 64)        # c / r_av: B1 sweep, stable and accurate for every N, R and rho
GRID = (2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128)            # c / r_av values of the test function
CANDIDATES = tuple(c for c in GRID if ADMISSIBLE[0] <= c <= ADMISSIBLE[1])
OUT_DIR = os.path.join(HERE, 'sensitivity', 'shape_selection')

# case label, folder, channels in the order of the paper: (segment, name)
CASES = [('Single channel', 'syntheticChannel', [(0, '--')]),
         ('Backwater network', 'akanYenNetwork', [(0, 'Channels 1--6')]),       # six identical channels
         ('LMR', 'lowerMississippi_paper', [(0, 'Tarbert Landing--Head of Passes')]),
         ('Midwest 2008', '2008flood_stLouis', [(0, 'Mississippi, Saverton--Grafton'),
                                                (2, 'Illinois, Valley City--Grafton'),
                                                (3, 'Mississippi, Grafton--Madison'),
                                                (1, 'Missouri, Hermann--Madison'),
                                                (4, 'Mississippi, Madison--Chester')])]


def spacing(case, seg):
    """N, L (m), R = dx_max / dx_min and rho = largest ratio of adjacent spacings of a channel"""
    x = np.loadtxt(os.path.join(HERE, case, f'segment{seg}', 'geo', 'nodes'))
    dx = np.diff(x)
    rho = np.max(np.maximum(dx[1:] / dx[:-1], dx[:-1] / dx[1:])) if len(dx) > 1 else 1.0
    return len(x), x[-1] - x[0], dx.max() / dx.min(), rho


def select(case, candidates, write):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return select_shape_coefficients(os.path.join(HERE, case), candidates=candidates, write=write)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rows, curves = [], []
    for label, case, channels in CASES:
        coefs, _ = select(case, CANDIDATES, write=True)                 # writes <case>/shape
        _, errs = select(case, GRID, write=False)
        for seg, name in channels:
            N, L, R, rho = spacing(case, seg)
            rows.append((label, name, N, L, R, rho, coefs[seg], coefs[seg] * L / (N - 1)))
            curves.append([label, name, seg, N] + list(errs[seg]))
        print(f'{case}/shape: ' + ' '.join(f'{c:g}' for c in coefs))

    with open(os.path.join(OUT_DIR, 'test_errors.csv'), 'w') as fh:
        fh.write('case,channel,segment,N,' + ','.join(f'c{c:g}' for c in GRID) + '\n')
        for r in curves:
            fh.write(f'{r[0]},"{r[1]}",{r[2]},{r[3]},' + ','.join(f'{e:.3e}' for e in r[4:]) + '\n')

    tex = [r'\begin{table}[htbp]',
           r'\caption{Node distribution and shape parameter of the channels. $N$: number of nodes; $L$: channel '
           r'length; $R=\Delta x_{max}/\Delta x_{min}$; $\rho$: largest ratio of adjacent node spacings; '
           r'$c/r_{av}$: value selected by the test function within the admissible range '
           f'${ADMISSIBLE[0]}\\leq c/r_{{av}}\\leq {ADMISSIBLE[1]}$.' r'}',
           r'\label{tab:shape}', r'\centering', r'\small', r'\setlength{\tabcolsep}{4pt}',
           r'\begin{tabular}{llrrrrr}', r'\toprule',
           r'Case & Channel & $N$ & $L$ & $R$ & $\rho$ & $c/r_{av}$ \\',
           r' & & & (km) & & & \\', r'\midrule']
    last = None
    for label, name, N, L, R, rho, coef, c in rows:
        if last is not None and label != last:
            tex.append(r'\midrule')
        tex.append(f'{label if label != last else ""} & {name} & {N} & {L / 1e3:.1f} & {R:.1f} & {rho:.1f} & '
                   f'{coef:g} \\\\')
        last = label
    tex += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    path = os.path.join(OUT_DIR, 'table_shape.tex')
    open(path, 'w').write('\n'.join(tex) + '\n')
    print(f'candidates c / r_av: {CANDIDATES}\nwrote {os.path.relpath(path, HERE)} and test_errors.csv')


if __name__ == '__main__':
    main()
