import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates

fPath = f'../data/hydros/'
hydros = ['Saverton_mississippi', 'valleyCity_illinois', 'hermann_missouri']
stage = f'chester_stage_mississippi'

fig, ax1 = plt.subplots(1, figsize=(7,4))

ax2 = ax1.twinx()
ax2.set_ylabel('Stage (m)')

df = pd.read_csv(f'../data/hydros/saverton_mississippi', sep = '\t')
df['date'] = pd.to_datetime(df['Date / Time'], format='mixed')
df['Q-cms'] = df['Flow (CFS)'].str.replace(",", "").astype(float) / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
ax1.plot(df['date'], df['Q-cms'], 'b', zorder=1, label='Inflow at Saverton, Mississippi', linewidth=1.6, color='k', alpha=1)

df = pd.read_csv(f'../data/hydros/valleyCity_illinois', skiprows = 24, sep = '\t')
df['date'] = pd.to_datetime(df['20d'], format='%Y-%m-%d')
df['Q-cms'] = df['14n'] / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
# print(df)
df_to_save = df[['seconds', 'Q-cms']]
# df_to_save.to_csv('../segment2/geo/boundary_Q', header=None, index=False, sep=' ')
df_to_save.to_csv(f'../data/hydros/to_hec/valleyCity_flow', columns=['Q-cms'], header=None, index=False, sep=' ')
ax1.plot(df['date'], df['Q-cms'], 'r', zorder=2, label='Inflow at Valley City, Illinois', linewidth=1.6, color='b', alpha=1)

df = pd.read_csv(f'../data/hydros/hermann_missouri', skiprows = 25, sep = '\t')
df['date'] = pd.to_datetime(df['20d'], format='%Y-%m-%d')
df['Q-cms'] = df['14n'] / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
df_to_save = df[['seconds', 'Q-cms']]
ax1.plot(df['date'], df['Q-cms'], 'b', zorder=3, label='Inflow at Hermann, Missouri', linewidth=1.6, color='r', alpha=1)

df = pd.read_csv(f'../data/hydros/chester_stage_mississippi', skiprows=27, sep=',')
df['date'] = pd.to_datetime(df["20d"], format='mixed')
df['stage-ft'] = pd.to_numeric(df['14n'], errors='coerce')
df = df.dropna()
df['stage-m'] = df['stage-ft'] / 3.281
ax2.plot(df['date'], df['stage-m'], 'b', zorder=4, label='Stage at Chester, Mississippi', linewidth=1.6, color='orange', alpha=1)

xmin = df['date'].iloc[0]
xmax = df['date'].iloc[-1]
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(0, 18000)
ax2.set_ylim(0, 20)
ax1.xaxis.set_major_locator(
    mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11])
)
ax1.set_ylabel('Discharge (cms)')
plt.xlabel('Date')
plt.tight_layout()
ax1.legend(loc="upper left")
ax2.legend()
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
plt.savefig('ex4_inflows.pdf')
plt.show()
