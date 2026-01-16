import pandas as pd
import numpy as np

lat_idx = [128,139,179,196,211,240,250,252,256,258,261,265,267,269,271,274,276]

for i in range(len(lat_idx)):
    filePath = f'../../../Junk/MESH_Code_irregular-master/Mississippi_River_11_years_20200511/lateral_in_min_diffusive/lateral_0{lat_idx[i]}.txt'
    df = pd.read_csv(filePath, sep='\t', names=['t-min', 'Q-cms'])
    df['t-sec'] = df['t-min'] * 60
    start = pd.Timestamp('2009-01-01 00:00:00')
    start_2011 = pd.Timestamp('2011-01-01 00:00:00')
    start_2012 = pd.Timestamp('2012-01-01 00:00:00')
    df['date'] = start + pd.to_timedelta(df['t-sec'], unit='s')
    df_2011 = df[(df['date'] >= start_2011) & (df['date'] <= start_2012)]
    df_2011['t_res-s'] = df_2011['t-sec'] - df_2011['t-sec'].iloc[0]
    df_2011.to_csv(f'laterals/qlat{i+1}', columns=['t_res-s', 'Q-cms'], index=False, sep='\t', header=None)

# print(df_2011.head())