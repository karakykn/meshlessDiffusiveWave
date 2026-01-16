import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fPath = f'../data/westbayou_stage.csv'
df = pd.read_csv(fPath, skiprows=20).iloc[:-2,:]
df_a = pd.read_csv(f'../segment0/geo/lateralDatas/qlat1', sep='\t', names=['t', 'Q'])
df_wb = df_a.copy()
df = df.fillna(0)
# df.loc[df['synthetic D'] == '#VALUE!', 'synthetic D'] = 0
df.loc[df['Stage (Ft)'] == 'M', 'Stage (Ft)'] = 0
peakS = 3.61
df['Q'] = df_a['Q']
peak = 51000 * 0.02831683199881
pwr = 3/2
df['Q'] = peak * (df['Stage (Ft)'].astype(float)) ** pwr / (df['Stage (Ft)'].astype(float).max()) ** pwr
df_wb['Q'] = -df['Q']
# print(df_wb)
df_wb.to_csv(f'../segment0/geo/lateralDatas/qlat7', sep='\t', header=None, columns=['t', 'Q'], index=False)

# plt.plot(df_wb['t'], df['Q'])
# plt.show()