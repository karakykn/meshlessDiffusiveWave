import pandas as pd
import pandas as pd

df = pd.read_csv(f'../data/hydros/yolobypass.csv', skiprows=24, sep='\t')
df['Date'] = pd.to_datetime(df.iloc[:, 2], format='mixed')
df['Q-cms'] = -df.iloc[:, 3] / 35.31466621266132
df['Seconds'] = (df['Date'] - df['Date'].iloc[0]).dt.total_seconds()
df.to_csv(f'../segment1/geo/lateralDatas/qlat0', header=None, columns=['Seconds', 'Q-cms'], sep=' ', index=False)