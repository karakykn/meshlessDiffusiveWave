import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fPath = f'laterals/qlat1'
df = pd.read_csv(fPath, sep='\t', names=['t', 'Q'])
df['Q'] = df['Q']
df.to_csv(f'../segment0/geo/lateralDatas/qlat1', sep='\t', header=None, columns=['t', 'Q'], index=False)