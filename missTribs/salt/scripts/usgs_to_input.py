import numpy as np
import pandas as pd

file = f'../../1dRoute/hidros/saltRivefnearnewLondon.csv'
dat = pd.read_csv(file, sep=',')
# dat['h-m'] = dat['h-feet'] / 3.281
dat['Q-m3s'] = dat['value'] / 35.31466621266132
dat['time'] = pd.to_datetime(dat['time'], format='%Y-%m-%d')
dat['seconds'] = (dat['time'] - dat['time'].iloc[-1]).dt.total_seconds()
df_to_save = dat[['seconds', 'Q-m3s']]
dat = dat.sort_values(by='seconds')

df_to_save = dat[['seconds', 'Q-m3s']]
df_to_save.to_csv('../segment0/geo/boundary_Q', header=None, index=False, sep=' ')
print(df_to_save)