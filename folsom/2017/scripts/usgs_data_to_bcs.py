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

fig, ax = plt.subplots(2, 1, figsize=(7, 8))

# folsom
df = pd.read_csv(f'../data/hydros/folsom_hourly_outflow.csv', sep = ',')
# df = df.iloc[::24]
df['Date'] = pd.to_datetime(df.iloc[:, 0], format='mixed')
df['Q-cms'] = df.iloc[:, 1] / 35.31466621266132
df_2017 = df[(df['Date'] < pd.Timestamp(2017, 2, 25)) & (df['Date'] > pd.Timestamp(2017, 2, 1))]
df_2017 = df[(df['Date'] < pd.Timestamp(2017, 2, 25)) & (df['Date'] > pd.Timestamp(2017, 2, 1))]
df_2017['Seconds'] = (df_2017['Date'] - df_2017['Date'].iloc[0]).dt.total_seconds()
df_2006 = df[(df['Date'] < pd.Timestamp(2006, 1, 15)) & (df['Date'] > pd.Timestamp(2005, 12, 15))]
df_2006 = df[(df['Date'] < pd.Timestamp(2006, 1, 15)) & (df['Date'] > pd.Timestamp(2005, 12, 15))]
df_2006['Seconds'] = (df_2006['Date'] - df_2006['Date'].iloc[0]).dt.total_seconds()

ax[0].plot(df_2017['Date'], df_2017['Q-cms'], 'k', zorder=1, label = 'Folsom, American River')
ax[1].plot(df_2006['Date'], df_2006['Q-cms'], 'k', zorder=1, label = 'Folsom, American River')

df_to_save = df_2017[['Seconds', 'Q-cms']]
df_to_save.to_csv('../segment0/geo/boundary_Q', header=None, index=False, sep=' ')

### verona
df = pd.read_csv(f'../data/hydros/verona_usgs.csv', sep = '\t', header = None)
df['Date'] = pd.to_datetime(df.iloc[:,2], format='mixed')
df['Q-cms'] = df.iloc[:,3] / 35.31466621266132
df_2017_v = df[(df['Date'] < pd.Timestamp(2017, 2, 25)) & (df['Date'] > pd.Timestamp(2017, 2, 1))]
df_2017_v = df[(df['Date'] < pd.Timestamp(2017, 2, 25)) & (df['Date'] > pd.Timestamp(2017, 2, 1))]
df_2017_v['Seconds'] = (df_2017_v['Date'] - df_2017_v['Date'].iloc[0]).dt.total_seconds()
df_2006 = df[(df['Date'] < pd.Timestamp(2006, 1, 15)) & (df['Date'] > pd.Timestamp(2005, 12, 15))]
df_2006 = df[(df['Date'] < pd.Timestamp(2006, 1, 15)) & (df['Date'] > pd.Timestamp(2005, 12, 15))]
df_2006['Seconds'] = (df_2006['Date'] - df_2006['Date'].iloc[0]).dt.total_seconds()

ax[0].plot(df_2017_v['Date'], df_2017_v['Q-cms'], 'b', zorder=1, label = 'Verona, Sacramento River')
ax[1].plot(df_2006['Date'], df_2006['Q-cms'], 'b', zorder=1, label = 'Verona, Sacramento River')

df_to_save = df_2017_v[['Seconds', 'Q-cms']]
df_to_save.to_csv('../segment1/geo/boundary_Q', header=None, index=False, sep=' ')

ax[0].set_title('2017 Event')
ax[1].set_title('2006 Event')
ax[0].set_ylabel('Discharge (cms)')
ax[1].set_ylabel('Discharge (cms)')
ax[1].set_xlabel('Date')
ax[0].legend()
ax[1].legend()
# plt.show()

df_2017["Year"] = df_2017["Date"].dt.year
df_2017["dayofyear"] = df_2017["Date"].dt.dayofyear

maxs = pd.DataFrame({
    'Year': df_2017['Year'].unique(),
    'Day': pd.NA,
    'max-Q': pd.NA
})
for year, group in df_2017.groupby("Year"):
    idx = maxs.index[maxs["Year"] == year]
    maxs['max-Q'].iloc[idx] = group['Q-cms'].max()
    maxs['Day'].iloc[idx] = group['dayofyear'].iloc[group['Q-cms'].argmax()]
    # ax.plot(group["dayofyear"], group["Q-cms"], label=year)
maxs_sorted = maxs.sort_values(by='max-Q', ascending=False)
print(maxs_sorted.head())