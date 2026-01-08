import numpy as np

caseName = 'lowerMississippi'
segtmentId = 0

L = 85900
N = 201
S_0 = 1e-5
n = 2.8e-2
Q_0 = 8381.786712
h_0 = 16

np.savetxt(f'{caseName}/segment{segtmentId}/geo/nodes', np.linspace(0, L, N))
np.savetxt(f'{caseName}/segment{segtmentId}/geo/slopes', np.ones(N) * S_0)
np.savetxt(f'{caseName}/segment{segtmentId}/geo/mannings_n', np.ones(N) * n)
np.savetxt(f'{caseName}/segment{segtmentId}/run/Q', np.ones(N) * Q_0)
np.savetxt(f'{caseName}/segment{segtmentId}/run/h', np.ones(N) * h_0)
np.savetxt(f'{caseName}/segment{segtmentId}/geo/xsInfo', np.zeros(N, dtype=int), fmt='%d')