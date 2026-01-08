import pandas as pd
import numpy as np
import datetime as dt

fpath = f'/Users/ismetkarakan/Documents/Doctorate/Junk/MESH_Code_irregular-master/Mississippi_River_11_years_20200511/EJ_wl_2009_2019.txt'
df = pd.read_csv(fpath, sep='\t', names=['time-min', 'waterL'])

start_time = pd.Timestamp('2009-01-01 00:00:00')

df['time-date'] = start_time + pd.to_timedelta(df['time-min'], unit='m')
start_2011 = '2011-01-01'
end_2011   = '2012-01-01'
df_2011 = df[(df['time-date'] >= start_2011) & (df['time-date'] <= end_2011)]

thalweg = 1.296922031000000075e+01
df_2011['flow-depth'] = df_2011['waterL'] + thalweg
df_2011['time-sec'] = (df_2011['time-date'] - df_2011['time-date'].iloc[0]).dt.total_seconds()

df_2011[['time-sec', 'flow-depth']].to_csv(f'inputs/downstreamFlowDepth', index=False, header=None, sep='\t')