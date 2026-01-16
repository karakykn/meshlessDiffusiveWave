import numpy as np
import pandas as pd

df = pd.read_csv(f'../data/hydros/grafton_mississippi', skiprows=26, sep='\t')
df['Date'] = pd.to_datetime(df["20d"], format='mixed')
df['stage-ft'] = pd.to_numeric(df['14n.1'], errors='coerce')
df = df.dropna()
datum = 403.32
df['stage-ft'] = df['stage-ft'] + datum
df['stage-m'] = df['stage-ft'] / 3.281

df.to_csv(f'grafton_data.csv', columns=['Date', 'stage-m'], index=False)