import numpy as np
import pandas as pd

data = pd.read_csv('../data/baton_belle.csv')
output_path = f'../data/xs/'
xs = np.zeros((4,3))
# print(data.shape[0])
for i in range(data.shape[0]):
# for i in range(1):
    xs[0, 1] = -data['BtmWdth'].loc[i] / 2
    xs[0, 2] = data['BtmWdth'].loc[i] / 2
    xs[1, 0] = data['trapDepth'].loc[i]
    xs[2, 0] = data['trapDepth'].loc[i] + .001
    xs[1, 1] = -data['TopWdth'].loc[i] / 2
    xs[1, 2] = data['TopWdth'].loc[i] / 2
    xs[2, 1] = -data['TopWdthCC'].loc[i] / 2
    xs[2, 2] = data['TopWdthCC'].loc[i] / 2
    xs[3, 0] = 70
    xs[3, 1] = -data['TopWdthCC'].loc[i] / 2
    xs[3, 2] = data['TopWdthCC'].loc[i] / 2
    np.savetxt(f'{output_path}xs{i}', xs)

# print(xs)
