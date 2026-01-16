import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fPath = f'../data/fort_stage.csv'
df = pd.read_csv(fPath, skiprows=21).iloc[:-1,:]
df_a = pd.read_csv(f'../segment0/geo/lateralDatas/qlat1', sep='\t', names=['t', 'Q'])
df_fort = df_a.copy()
df = df.fillna(0)
# df.loc[df['synthetic D'] == '#VALUE!', 'synthetic D'] = 0
df.loc[df['Stage (Ft)'] == 'M', 'Stage (Ft)'] = 0
# print(df['Stage (Ft)'].max())
limS = 3
mask = df['Stage (Ft)'].astype(float) > limS
df['Q'] = df_a['Q']
peak = 100000 * 0.02831683199881
pwr = 3/2
df.loc[mask, 'Q'] = peak * (df['Stage (Ft)'].astype(float) - limS) ** pwr / (df['Stage (Ft)'].astype(float).max() - limS) ** pwr
df_fort['Q'] = -df['Q']
# print(df_fort)
df_fort.to_csv(f'../segment0/geo/lateralDatas/qlat7', sep='\t', header=None, columns=['t', 'Q'], index=False)

# plt.plot(df_fort['t'], df['Q'])
# plt.show()