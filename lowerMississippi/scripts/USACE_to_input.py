import pandas as pd

file = '../data/tarbertQ.csv'
dat = pd.read_csv(file, sep=',')

df_to_save = dat[['seconds', 'cms']]
df_to_save.to_csv('../segment0/geo/boundaryQ', header=None, index=False, sep=' ')
