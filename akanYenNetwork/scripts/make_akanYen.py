"""
Akan and Yen (1981), Example I: case folder for the meshless model (source.py).

Network (Fig. 5): channels 1, 2, 3 join channel 5; channels 5 and 4 join channel 6 (outlet).
Segment numbering: channel k -> segment k-1.
Table 1: L = 600 m, rectangular sections, slope, width and Manning n per channel.
Inflow (Fig. 6): channels 1 and 4 triangular, 3 -> 12 m3/s at 10 min -> 3 m3/s at 20 min;
channels 2 and 3 constant 2 m3/s. Initial condition: steady flow (3, 2, 2, 3, 7, 10 m3/s).
Downstream boundary: normal depth (uniform flow) at the exit of channel 6.
Akan and Yen: dx = 60 m, dt = 60 s. Here dx = 60 m (11 nodes per channel); dt from the Courant number.

usage: python make_akanYen.py <case folder> [dx] [end time s]
"""
import os, sys
import numpy as np

case = sys.argv[1] if len(sys.argv) > 1 else 'akanYen1981'
dx = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
t_end = float(sys.argv[3]) if len(sys.argv) > 3 else 2400.0

L = 600.0
S0 = [0.0005, 0.0005, 0.0005, 0.0005, 0.0010, 0.0010]
B = [5.0, 5.0, 5.0, 5.0, 8.0, 10.0]
n = [0.0138, 0.0207, 0.0207, 0.0138, 0.0141, 0.0125]
Q0 = [3.0, 2.0, 2.0, 3.0, 7.0, 10.0]
network = [(0, 4), (1, 4), (2, 4), (4, 5), (3, 5)]

tri = np.array([[0, 3.0], [600, 12.0], [1200, 3.0], [max(t_end, 2400) + 3600, 3.0]])
const2 = np.array([[0, 2.0], [max(t_end, 2400) + 3600, 2.0]])
inflow = {0: tri, 1: const2, 2: const2, 3: tri}


def normal_depth(Q, b, n_, S):
    h = 1.0
    for _ in range(100):
        A, P = b * h, b + 2 * h
        f = A * (A / P) ** (2 / 3) * np.sqrt(S) / n_ - Q
        dA, dP = b, 2.0
        df = np.sqrt(S) / n_ * (5 / 3 * (A / P) ** (2 / 3) * dA - 2 / 3 * (A / P) ** (5 / 3) * dP)
        h -= f / df
    return h


os.makedirs(case, exist_ok=True)
with open(os.path.join(case, 'input'), 'w') as f:
    f.write('CFL:\n1\n\nRBF type: (0 TPS, 1 MQ, 2 IMQ)\n1\n\n'
            'shape parameter coefficient MQ (c = beta * r_av), beta: (beta is beta of tps if tps is on)\n16\n\n'
            f'simulation end time (seconds):\n{int(t_end)}\n\nprint step:\n1\n\n'
            'downstream bc (0: normal flow, 1: flow depth specified):\n0\n')
np.savetxt(os.path.join(case, 'network'), np.array(network), fmt='%d')

x = np.arange(0, L + dx / 2, dx)
N = len(x)
for s in range(6):
    g = os.path.join(case, f'segment{s}', 'geo')
    os.makedirs(os.path.join(g, 'xSecs'), exist_ok=True)
    r0 = os.path.join(case, f'segment{s}', 'run', '0')
    os.makedirs(r0, exist_ok=True)
    np.savetxt(os.path.join(g, 'nodes'), x)
    np.savetxt(os.path.join(g, 'slopes'), np.full(N, S0[s]))
    np.savetxt(os.path.join(g, 'mannings_n'), np.full(N, n[s]))
    np.savetxt(os.path.join(g, 'xsInfo'), np.zeros(N, dtype=int), fmt='%d')
    b = B[s]
    # rectangular main channel up to 2.5 m (depths here stay below 1.5 m); nominal 1:100 overbank above,
    # same convention as the synthetic cases (never reached)
    xs = np.array([[0.0, 0.0, b], [2.5, 0.0, b], [2.6, -10.0, b + 10.0], [3.0, -10.0, b + 10.0]])
    np.savetxt(os.path.join(g, 'xSecs', 'xs0'), xs, fmt='%g')
    if s in inflow:
        np.savetxt(os.path.join(g, 'boundary_Q'), inflow[s], fmt='%g')
    else:
        np.savetxt(os.path.join(g, 'boundary_Q'), np.array([[0, Q0[s]]]), fmt='%g')
    hn = normal_depth(Q0[s], b, n[s], S0[s])
    np.savetxt(os.path.join(g, 'boundary_h'), np.array([[0, hn], [t_end + 3600, hn]]))
    np.savetxt(os.path.join(r0, 'Q'), np.full(N, Q0[s]))
    np.savetxt(os.path.join(r0, 'h'), np.full(N, hn))
    print(f'channel {s+1}: N {N}, normal depth at Q0 = {Q0[s]:g}: {hn:.4f} m, '
          f'at 12 m3/s: {normal_depth(12.0, b, n[s], S0[s]):.4f} m')
