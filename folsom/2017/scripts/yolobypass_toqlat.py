import pandas as pd
import pandas as pd

df = pd.read_csv(f'../data/hydros/abdeltaHeight', skiprows=27, sep='\t')
df['Date'] = pd.to_datetime(df.iloc[:, 2], format='mixed')
df['h-m'] = df.iloc[:, 4] / 3.281
df['Seconds'] = (df['Date'] - df['Date'].iloc[0]).dt.total_seconds()
df.to_csv(f'../../2017_single/segment0/geo/boundary_h', header=None, columns=['Seconds', 'h-m'], sep=' ', index=False)