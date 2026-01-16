import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "legend.fontsize": 9,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "lines.linewidth": 1.2
})

fig, ax = plt.subplots(figsize=(7, 4))

# segment 1 inflow hydro
df = pd.read_csv(f'../data/hydros/hermann_missouri', skiprows = 25, sep = '\t')
df['date'] = pd.to_datetime(df['20d'], format='%Y-%m-%d')
df['Q-cms'] = df['14n'] / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
df_to_save = df[['seconds', 'Q-cms']]
# df_to_save.to_csv('../segment1/geo/boundary_Q', header=None, index=False, sep=' ')
# df_to_save.to_csv(f'../data/hydros/to_hec/hermann_flow', columns=['Q-cms'], header=None, index=False, sep=' ')

ax.plot(df['date'], df['Q-cms'], 'b', zorder=1)

# segment 2 inflow hydro
df = pd.read_csv(f'../data/hydros/valleyCity_illinois', skiprows = 24, sep = '\t')
df['date'] = pd.to_datetime(df['20d'], format='%Y-%m-%d')
df['Q-cms'] = df['14n'] / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
# print(df)
df_to_save = df[['seconds', 'Q-cms']]
# df_to_save.to_csv('../segment2/geo/boundary_Q', header=None, index=False, sep=' ')
# df_to_save.to_csv(f'../data/hydros/to_hec/valleyCity_flow', columns=['Q-cms'], header=None, index=False, sep=' ')

ax.plot(df['date'], df['Q-cms'], 'r', zorder=2)

# # segment 0 inflow hydro
df = pd.read_csv(f'../data/hydros/saverton_mississippi', sep = '\t')
df['date'] = pd.to_datetime(df['Date / Time'], format='mixed')
df['Q-cms'] = df['Flow (CFS)'].str.replace(",", "").astype(float) / 35.31466621266132
df['seconds'] = (df['date'] - df['date'].iloc[0]).dt.total_seconds()
# print(df)
df_to_save = df[['seconds', 'Q-cms']]
# df_to_save.to_csv('../segment0/geo/boundary_Q', header=None, index=False, sep=' ')
# df_to_save.to_csv(f'../data/hydros/to_hec/saverton_flow', columns=['Q-cms'], header=None, index=False, sep=' ')

ax.plot(df['date'], df['Q-cms'], 'k', zorder=0)
xmin = df['date'].iloc[0]
xmax = df['date'].iloc[-1]
ax.set_xlim(xmin, xmax)
# plt.tight_layout()
plt.xlabel('Date')
ax.xaxis.set_major_locator(
    mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11])
)
ax.xaxis.set_minor_locator(mdates.MonthLocator())
plt.show()