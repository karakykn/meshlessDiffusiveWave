"""
Runs meshless cases several times and reports the CPU time of the time loop (median over the runs).

    python3 time_meshless.py                                 # the four cases of the paper, 5 runs each
    python3 time_meshless.py 3 lowerMississippi_paper 2008flood_stLouis

Every run is a fresh python process with single-threaded numpy (BLAS/LAPACK), so that the CPU time is comparable
with the serial Fortran CNS code; its screen output goes to <case>/cpu_run_<k>.log. Like process.py, every run
rewrites the results in <case>/segment*/run. The table is also written to timing_meshless.csv.
"""
import csv
import json
import os
import statistics
import subprocess
import sys

THETA = 0.5
THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS',
               'NUMEXPR_NUM_THREADS')
KEYS = ('nodes', 'steps', 'cpu_loop', 'cpu_output', 'cpu_compute', 'wall_loop')
HERE = os.path.dirname(os.path.abspath(__file__))


def run_case(case, repeats):
    env = dict(os.environ, **{v: '1' for v in THREAD_VARS})
    code = (f"from source import *\nn = Network('{case}')\nn.solve({THETA})\n"
            "import json; print('CPU_JSON ' + json.dumps(n.cpu_time))")
    runs = []
    for k in range(1, repeats + 1):
        log = os.path.join(HERE, case, f'cpu_run_{k}.log')
        with open(log, 'w') as f:
            subprocess.run([sys.executable, '-c', code], cwd=HERE, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True)
        line = [l for l in open(log) if l.startswith('CPU_JSON ')][-1]
        r = json.loads(line[9:])
        runs.append(r)
        print(f"  {case} run {k}: loop {r['cpu_loop']:.4f} s, without output {r['cpu_compute']:.4f} s")
    med = {key: statistics.median(r[key] for r in runs) for key in KEYS}
    med.update(case=case, runs=repeats)
    return med


if __name__ == '__main__':
    repeats = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cases = sys.argv[2:] or ['syntheticChannel', 'akanYenNetwork', 'lowerMississippi_paper', '2008flood_stLouis']
    rows = [run_case(c, repeats) for c in cases]
    print(f"\nMeshless, median of {repeats} runs (CPU s)")
    print(f"{'case':<20}{'nodes':>6}{'steps':>8}{'loop':>10}{'output':>10}{'w/o out':>10}{'ms/step':>10}{'wall':>10}")
    for r in rows:
        print(f"{r['case']:<20}{r['nodes']:6.0f}{r['steps']:8.0f}{r['cpu_loop']:10.4f}{r['cpu_output']:10.4f}"
              f"{r['cpu_compute']:10.4f}{1e3 * r['cpu_compute'] / r['steps']:10.4f}{r['wall_loop']:10.4f}")
    with open(os.path.join(HERE, 'timing_meshless.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['case', 'runs'] + list(KEYS))
        w.writeheader()
        w.writerows(rows)
