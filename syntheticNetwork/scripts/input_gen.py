import numpy as np
import pandas as pd

mann = [.0446, .0446, .0446]
# mann = [.0125, .0125, .0125]
n = [20, 40, 40]
L = [5000, 10000, 20000]
baseF = [20, 20, 40]
baseh = [4, 4, 4]
slope = 1e-4

for i in range(3):
    outPath = f'../segment{i}/geo/'
    # locs = np.loadtxt(f'{outPath}nodes')
    # nn = locs.shape[0]
    nn = n[i] + 1
    locs = np.linspace(0,L[i], nn)
    bQ, bh = np.ones(nn) * baseF[i], np.ones(nn) * baseh[i]
    man = np.ones(nn) * mann[i]
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
