import numpy as np
import pandas as pd

mann = .0446
# mann = .0125
n = 20
L = 20000
baseF = 20
baseh = 4
slope = 1e-4

for i in range(1):
    outPath = f'../segment{i}/geo/'
    # locs = np.loadtxt(f'{outPath}nodes')
    # nn = locs.shape[0]
    nn = n + 1
    locs = np.linspace(0,L, nn)
    bQ, bh = np.ones(nn) * baseF, np.ones(nn) * baseh
    man = np.ones(nn) * mann
    slop = np.ones(nn) * slope
    xsinf = np.zeros(nn)
    for i in range(nn):
        xsinf[i] = 0
    slop[-1] = slop[-2]
    xsinf[-1] = xsinf[-2]

    np.savetxt(outPath + '../run/0/h', bh)
    np.savetxt(outPath + '../run/0/Q', bQ)
    np.savetxt(outPath + 'nodes', locs)
    np.savetxt(outPath + 'mannings_n', man)
    np.savetxt(outPath + 'slopes', slop)
    np.savetxt(outPath + 'xsInfo', xsinf, fmt='%d')
