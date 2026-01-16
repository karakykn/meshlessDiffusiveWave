import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fPath = f'../data/lahache_stage.csv'
df = pd.read_csv(fPath, skiprows=19).iloc[:-1, :]
df_a = pd.read_csv(f'../segment0/geo/lateralDatas/qlat1', sep='\t', names=['t', 'Q'])
df_bohemia = df_a.copy()
df = df.fillna(0)
# df.loc[df['synthetic D'] == '#VALUE!', 'synthetic D'] = 0
df.loc[df['Stage (Ft)'] == 'M', 'Stage (Ft)'] = 0
limS = 5.74
mask = df['Stage (Ft)'].astype(float) > limS
df['Q'] = df_a['Q']
peak = 40000 * 0.02831683199881
pwr = 3/2
df.loc[mask, 'Q'] = peak * (df['Stage (Ft)'].astype(float) - limS) ** pwr / (df['Stage (Ft)'].astype(float).max() - limS) ** pwr
df_bohemia['Q'] = -df['Q']
# print(df_bohemia)
df_bohemia.to_csv(f'../segment0/geo/lateralDatas/qlat6', sep='\t', header=None, columns=['t', 'Q'], index=False)

plt.plot(df_bohemia['t'], -df_bohemia['Q'])
plt.show()