import numpy as np
import pandas as pd

mann = [.0446, .0446, .0446]
# mann = [.0125, .0125, .0125]

for i in range(3):
    outPath = f'../segment{i}/geo/'
    locs = np.loadtxt(f'{outPath}nodes')
    nn = locs.shape[0]
    man = np.ones(nn) * mann[i]
    slop = np.zeros(nn)
    xsinf = np.zeros(nn)
    for i in range(nn):
        xsinf[i] = 0
    slop[-1] = slop[-2]
    xsinf[-1] = xsinf[-2]

    np.savetxt(outPath + 'mannings_n', man)
    # np.savetxt(outPath + 'slopes', slop)
    # np.savetxt(outPath + 'xsInfo', xsinf, fmt='%d')
