import numpy as np
import pandas as pd

df = pd.read_csv(f'../data/headofpasses_stage.csv', skiprows=25, sep=',', names=['d', 'stage-ft'])
df['Date'] = pd.to_datetime(df["d"], format='mixed')
df['stage-ft'] = pd.to_numeric(df['stage-ft'], errors='coerce')
df = df.dropna()
df['stage-m'] = df['stage-ft'] / 3.281
df['seconds'] = (df['Date'] - df['Date'].iloc[0]).dt.total_seconds()
alpha0 = -22.0501
df['flowdepth-m'] = df['stage-m'] - alpha0

df.to_csv(f'../segment0/geo/boundary_h', columns=['seconds', 'flowdepth-m'], sep=' ', index = False, header = None)