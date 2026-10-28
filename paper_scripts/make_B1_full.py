"""Full shape-parameter figure of the Supplemental Materials: MPE, MBE, dQp, dt_a and the failed runs against
c / r_av, one row per R, redrawn from the finished sweep (sensitivity/B1_final_model/B1_results.csv) with plot()
of process_B1.py. Nothing in sensitivity/ is written.
    python3 make_B1_full.py
Output: B1_full.pdf/.png (this folder)
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                           # meshlessDiffusive_Q
sys.path.insert(0, ROOT)
import process_B1 as b1                                # noqa: E402

b1.FIG_METRICS = ['mpe_hec', 'mbe_pct', 'dQp_hec_pct', 'dta_hec_min']
b1.FIG_FAILED_COLUMN = True
table = pd.read_csv(os.path.join(ROOT, 'sensitivity', 'B1_final_model', 'B1_results.csv'))
b1.plot(table, os.path.join(HERE, 'B1_full.png'))
print('Saved', os.path.join(HERE, 'B1_full.pdf'))
