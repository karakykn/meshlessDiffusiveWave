import numpy as np
import pandas as pd

df = pd.read_csv(f'../data/hydros/chester_stage_mississippi', skiprows=27, sep=',')
df['Date'] = pd.to_datetime(df["20d"], format='mixed')
df['stage-ft'] = pd.to_numeric(df['14n'], errors='coerce')
df = df.dropna()
df['stage-m'] = df['stage-ft'] / 3.281
df['seconds'] = (df['Date'] - df['Date'].iloc[0]).dt.total_seconds()
h0 = 1.446386185842965943e+01
alpha0 = h0 - df['stage-m'][0]
df['flowdepth-m'] = df['stage-m'] + alpha0

pass

df.to_csv(f'../segment4/geo/boundary_h', columns=['seconds', 'flowdepth-m'], sep=' ', index = False, header = None)
df.to_csv(f'../data/hydros/to_hec/chester_stage', columns=['stage-m'], index = False, header = None)