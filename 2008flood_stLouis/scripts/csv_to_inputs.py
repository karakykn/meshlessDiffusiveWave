import numpy as np
import pandas as pd

filePs = [f'../data/from_saverton_mississippi_edit.csv',
          f'../data/from_hermann_missouri_edit.csv',
          f'../data/from_valleycity_illinois_edit.csv',
          f'../data/from_grafton_mississippi_edit.csv',
          f'../data/from_madison_mississippi_edit.csv']

for i in range(5):
    outPath = f'../segment{i}/geo/'
    df = pd.read_csv(filePs[i])
    df = df[df.iloc[:, 0].notna()]
    nn = df.shape[0]

    nodes = np.zeros(nn+1)
    man = np.zeros(nn+1)
    slop = np.zeros(nn+1)
    xsinf = np.zeros(nn+1)
    for i in range(nn):
        nodes[i+1] = nodes[i] + df['Length'].iloc[i]
        slop[i] = df['So'].iloc[i]
        man[i] = df['n'].iloc[i]
        xsinf[i] = i
    slop[-1] = slop[-2]
    man[-1] = man[-2]
    xsinf[-1] = xsinf[-2]

    np.savetxt(outPath + 'nodes', nodes)
    np.savetxt(outPath + 'mannings_n', man)
    np.savetxt(outPath + 'slopes', slop)
    np.savetxt(outPath + 'xsInfo', xsinf, fmt='%d')
    xs = np.zeros((4, 3))
    df['trapDepth'] = (df['Tw'] - df['B']) * df['z'] / 2
    # print(data.shape[0])
    for i in range(df.shape[0]):
        xs[0, 1] = -df['B'].loc[i] / 2
        xs[0, 2] = df['B'].loc[i] / 2
        xs[1, 0] = df['trapDepth'].loc[i]
        xs[2, 0] = df['trapDepth'].loc[i] + .001
        xs[1, 1] = -df['Tw'].loc[i] / 2
        xs[1, 2] = df['Tw'].loc[i] / 2
        xs[2, 1] = -df['TwCC'].loc[i] / 2
        xs[2, 2] = df['TwCC'].loc[i] / 2
        xs[3, 0] = 70
        xs[3, 1] = -df['TwCC'].loc[i] / 2
        xs[3, 2] = df['TwCC'].loc[i] / 2
        np.savetxt(f'{outPath}/xSecs/xs{i}', xs)