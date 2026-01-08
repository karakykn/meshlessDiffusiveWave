import pandas as pd
import numpy as np

net = pd.read_csv(f'../network', header=None, sep=' ')
all = net.stack().unique()
downstreams = net.iloc[:,1].unique()

for i in all:
    if i not in downstreams:
        fPath = f'../segment{i}/geo/boundary_Q'
        nodeN = np.loadtxt(f'../segment{i}/geo/nodes').shape[0]
        df = pd.read_csv(fPath, names=['sec', 'Q'], sep=' ')
        Q = df.iloc[0,1]
        oPath = f'../segment{i}/run/0/Q'
        np.savetxt(oPath, np.ones(nodeN) * Q)
        np.savetxt(f'../segment{i}/run/0/h', np.ones(nodeN))
    else:
        ups = net[net.iloc[:,1] == i].iloc[:, 0]
        totQ = 0
        for j in ups:
            fPath = f'../segment{j}/geo/boundary_Q'
            df = pd.read_csv(fPath, names=['sec', 'Q'], sep=' ')
            totQ += df.iloc[0, 1]
        nodeN = np.loadtxt(f'../segment{i}/geo/nodes').shape[0]
        oPath = f'../segment{i}/run/0/Q'
        np.savetxt(oPath, np.ones(nodeN) * totQ)
        np.savetxt(f'../segment{i}/geo/boundary_Q', np.array([[0, totQ]]))
        np.savetxt(f'../segment{i}/run/0/h', np.ones(nodeN))


