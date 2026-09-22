import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = f'../../1dRoute/hidros/southFabiusnearTaylor.csv'
dat = pd.read_csv(file, sep=',')
# dat['h-m'] = dat['h-feet'] / 3.281
dat['Q-m3s'] = dat['value'] / 35.31466621266132
dat['time'] = pd.to_datetime(dat['time'], format='%d.%m.%Y')
dat['seconds'] = (dat['time'] - dat['time'].iloc[-1]).dt.total_seconds()
df_to_save = dat[['seconds', 'Q-m3s']]
dat = dat.sort_values(by='seconds')

df_to_save = dat[['seconds', 'Q-m3s']]
df_to_save.to_csv('../segment1/geo/boundary_Q', header=None, index=False, sep=' ')
plt.plot(df_to_save['seconds'], df_to_save['Q-m3s'])
plt.show()
print(df_to_save)