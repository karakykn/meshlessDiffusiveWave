import pandas as pd
import numpy as np
import datetime as dt

fpath = f'/Users/ismetkarakan/Documents/Doctorate/meshlessDiffusive/lowerMississippi/data/headofpasses_stage.csv'
df = pd.read_csv(fpath, skiprows=24)
df = df.iloc[:-2]

dum_df = pd.read_csv(f'../segment0/geo/lateralDatas/qlat0', names=['t','Q'], sep='\t')

# thalweg = 14
df['t-sec'] = dum_df['t']
df['Stage (Ft)'] = pd.to_numeric(df['Stage (Ft)'], errors='coerce')
# df.loc[df['Stage (Ft)'] == 'M', 'Stage (Ft)'] = (s.shift(1) + s.shift(-1)) / 2
h0 = 2.228387855117944483e+01
df['Stage (m)'] = df['Stage (Ft)'].astype(float) / 3.281
alpha0 = h0 - df['Stage (m)'][0]
df['flow depth-m'] = df['Stage (m)'] + alpha0
df.to_csv(f'../segment0/geo/boundar_h', index=False, header=None, sep='\t', columns=['t-sec', 'flow depth-m'])




# df = pd.read_csv(fpath, sep='\t', names=['time-min', 'waterL'])
#
# start_time = pd.Timestamp('2009-01-01 00:00:00')
#
# df['time-date'] = start_time + pd.to_timedelta(df['time-min'], unit='m')
# start_2011 = '2011-01-01'
# end_2011   = '2012-01-01'
# df_2011 = df[(df['time-date'] >= start_2011) & (df['time-date'] <= end_2011)]
#
# thalweg = 1.296922031000000075e+01
# df_2011['flow-depth'] = df_2011['waterL'] + thalweg
# df_2011['time-sec'] = (df_2011['time-date'] - df_2011['time-date'].iloc[0]).dt.total_seconds()
#
# df_2011[['time-sec', 'flow-depth']].to_csv(f'inputs/downstreamFlowDepth', index=False, header=None, sep='\t')