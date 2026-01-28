import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_verona = pd.read_csv(f'../data/hydros/verona_usgs.csv', header=None, sep='\t')
df_yolo = pd.read_csv(f'../data/hydros/yolobypass.csv', header=None, sep='\t', skiprows=25)
df_fairoaks = pd.read_csv(f'../data/hydros/fairoaks_usgs.csv', header=None, sep='\t', skiprows=25)
df_freeport = pd.read_csv(f'../data/hydros/freeport.csv', header=None, sep='\t', skiprows=25)
df_deepDivert = pd.read_csv(f'../data/hydros/freeport.csv', header=None, sep='\t', skiprows=25)

df_verona['Date'] = pd.to_datetime(df_verona.iloc[:,2], format='mixed')
df_verona['Q-cms'] = df_verona.iloc[:,3] / 35.31466621266132
df_verona = df_verona[(df_verona['Date'] < pd.Timestamp(2017, 2, 26)) & (df_verona['Date'] > pd.Timestamp(2017, 1, 31))]
df_verona = df_verona.reset_index(drop=True)

df_yolo['Date'] = pd.to_datetime(df_yolo.iloc[:,2], format='mixed')
df_yolo['Q-cms'] = df_yolo.iloc[:,3] / 35.31466621266132
df_yolo['-Q-cms'] = df_yolo['Q-cms']
df_yolo['Seconds'] = (df_yolo['Date'] - df_yolo['Date'].iloc[0]).dt.total_seconds()
df_to_save = df_yolo[['Seconds', '-Q-cms']]
df_to_save.to_csv('../data/hydros/yolo', header=None, index=False, sep=' ')

df_fairoaks['Date'] = pd.to_datetime(df_fairoaks.iloc[:,2], format='mixed')
df_fairoaks['Q-cms'] = df_fairoaks.iloc[:,3] / 35.31466621266132
df_freeport['Date'] = pd.to_datetime(df_freeport.iloc[:,2], format='mixed')
df_freeport['Q-cms'] = df_freeport.iloc[:,3] / 35.31466621266132

df_deepDivert['Date'] = pd.to_datetime(df_deepDivert.iloc[:,2], format='mixed')
df_deepDivert['Q-cms'] = df_verona['Q-cms'] - df_yolo['Q-cms'] + df_fairoaks['Q-cms'] - df_freeport['Q-cms']
df_deepDivert['Seconds'] = (df_deepDivert['Date'] - df_deepDivert['Date'].iloc[0]).dt.total_seconds()

df_to_save = df_deepDivert[['Seconds', 'Q-cms']]
df_to_save.to_csv('../data/hydros/deepdivert', header=None, index=False, sep=' ')

df_sacramentoIn = pd.read_csv(f'../data/hydros/freeport.csv', header=None, sep='\t', skiprows=25)
df_sacramentoIn['Date'] = pd.to_datetime(df_sacramentoIn.iloc[:,2], format='mixed')
df_sacramentoIn['Q-cms'] = df_verona['Q-cms'] - df_yolo['Q-cms']
df_sacramentoIn['Seconds'] = (df_sacramentoIn['Date'] - df_sacramentoIn['Date'].iloc[0]).dt.total_seconds()

df_to_save = df_sacramentoIn[['Seconds', 'Q-cms']]
df_to_save.to_csv('../data/hydros/sacramentoIn', header=None, index=False, sep=' ')

plt.plot(df_verona['Date'], df_verona['Q-cms'])
plt.plot(df_yolo['Date'], df_yolo['Q-cms'])
plt.plot(df_yolo['Date'], df_verona['Q-cms']-df_yolo['Q-cms'])
plt.show()