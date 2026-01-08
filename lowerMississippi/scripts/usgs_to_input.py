import numpy as np
import pandas as pd

file = f'../data/usgs'
dat = pd.read_csv(file, names = ['date', 'h-feet', 'Q-cfs', 'h-m', 'Q-cms'], header=None, sep='\t',)
dat['h-m'] = dat['h-feet'] / 3.281
dat['Q-m3s'] = dat['Q-cfs'] / 35.31466621266132
dat['date'] = pd.to_datetime(dat['date'], format='%d.%m.%Y')
dat['seconds'] = (dat['date'] - dat['date'].iloc[0]).dt.total_seconds()
df_to_save = dat[['seconds', 'Q-m3s']]
# Save to CSV
df_to_save.to_csv('boundary_Q_baton', header=None, index=False, sep=' ')
print(dat)