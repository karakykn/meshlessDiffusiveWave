import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, mu, sig):
    return (
        1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )

fpath = f'../segment0/geo/boundary_Q'
df = pd.read_csv(fpath, sep='\t', names=['t', 'Q'])

start = df.iloc[0, 1]
peak = df['Q'].max()
amp = peak - start
peak_T = df.iloc[df['Q'].argmax(), 0]

sigma = 10e-2
df['Q_art'] = gaussian(df['t'], peak_T, sigma * df.iloc[-1,0])
df['Q_art'] = (df['Q_art'] / df['Q_art'].max()) * amp + start
plt.plot(df['t'], df['Q_art'])
# plt.show()

df.to_csv(f'../segment0/geo/boundary_Q_normal', sep='\t', columns=['t', 'Q_art'], index=False, header = None)