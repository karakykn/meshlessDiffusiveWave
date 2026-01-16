import numpy as np
import pandas as pd

df = pd.read_csv(f'../data/fromTarbert_edit.csv').iloc[:95]
nn = df.shape[0]

nodes = np.zeros(nn+1)
man = np.zeros(nn+1)
slop = np.zeros(nn+1)
xsinf = np.zeros(nn+1)
for i in range(nn):
    nodes[i+1] = nodes[i] + df['Length'].iloc[i]
    slop[i] = df['So'].iloc[i]
    man[i] = df['n'].iloc[i]
slop[-1] = slop[-2]
man[-1] = man[-2]

np.savetxt('inputs/nodes', nodes)
np.savetxt('inputs/mannings_n', man)
np.savetxt('inputs/slopes', slop)
np.savetxt('inputs/xsInfo', xsinf)