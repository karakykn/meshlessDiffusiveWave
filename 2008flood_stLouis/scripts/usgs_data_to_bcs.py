import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(9, 5))

# segment 1 inflow hydro
df = pd.read_csv(f'../data/hydros/hermann_missouri', skiprows = 25, sep = '\t')
df['date'] = pd.to_datetime(df['20d'], format='%Y-%m-%d')
df['Q-cms'] = df['14n'] / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
df_to_save = df[['seconds', 'Q-cms']]
df_to_save.to_csv('../segment1/geo/boundary_Q', header=None, index=False, sep=' ')

ax.plot(df['seconds'], df['Q-cms'], 'b', zorder=1)

# segment 2 inflow hydro
df = pd.read_csv(f'../data/hydros/valleyCity_illinois', skiprows = 24, sep = '\t')
df['date'] = pd.to_datetime(df['20d'], format='%Y-%m-%d')
df['Q-cms'] = df['14n'] / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
# print(df)
df_to_save = df[['seconds', 'Q-cms']]
df_to_save.to_csv('../segment2/geo/boundary_Q', header=None, index=False, sep=' ')

ax.plot(df['seconds'], df['Q-cms'], 'r', zorder=2)

# # segment 0 inflow hydro
df = pd.read_csv(f'../data/hydros/saverton_mississippi', sep = '\t')
df['date'] = pd.to_datetime(df['Date / Time'], format='mixed')
df['Q-cms'] = df['Flow (CFS)'].str.replace(",", "").astype(float) / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
# print(df)
df_to_save = df[['seconds', 'Q-cms']]
df_to_save.to_csv('../segment0/geo/boundary_Q', header=None, index=False, sep=' ')

ax.plot(df['seconds'], df['Q-cms'], 'k', zorder=0)

plt.show()