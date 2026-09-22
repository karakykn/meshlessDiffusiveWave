import numpy as np
import pandas as pd

data = pd.read_csv(f'../../1dRoute/northRiverPalymra.csv')
output_path = f'../segment0/geo/xSecs/'
xs = np.zeros((4,3))
# print(data.shape[0])
data['trapDepth'] = (data['Tw'] - data['B']) * data['z'] / 2
for i in range(data.shape[0]):
    xs[0, 1] = -data['B'].loc[i] / 2
    xs[0, 2] = data['B'].loc[i] / 2
    xs[1, 0] = data['trapDepth'].loc[i]
    xs[2, 0] = data['trapDepth'].loc[i] + .001
    xs[1, 1] = -data['Tw'].loc[i] / 2
    xs[1, 2] = data['Tw'].loc[i] / 2
    xs[2, 1] = -data['TwCC'].loc[i] / 2
    xs[2, 2] = data['TwCC'].loc[i] / 2
    xs[3, 0] = 20
    xs[3, 1] = -data['TwCC'].loc[i] / 2
    xs[3, 2] = data['TwCC'].loc[i] / 2
    np.savetxt(f'{output_path}xs{i}', xs)

# print(xs)
