import pandas as pd

fPath = '../data/otflow_real/morganza_floodway'

# Original data (stays untouched)
df_ydk = pd.read_csv(
    '../segment0/geo/lateralDatas/qlat1',
    sep='\t',
    names=['time', 'Q']
)

# Morganza CSV
df = pd.read_csv(fPath)
df['Q-cms'] = df['discharge'] * 0.02831683199881

# Make sure date is datetime
df['date'] = pd.to_datetime(df['date'])

# Start of 2011
start_2011 = pd.Timestamp('2011-01-01')

# Seconds from 2011-01-01
df['time'] = (df['date'] - start_2011).dt.total_seconds().astype(int)

# 1. Find index of first match in df_ydk
first_time = df['time'].iloc[0]
matches = df_ydk.index[df_ydk['time'].eq(first_time)]

if matches.empty:
    raise ValueError(f"No matching time={first_time} found in df_ydk['time']")

start_idx = matches[0]
print("First match index in df_ydk:", start_idx)

# 2. Create df_n as a copy of df_ydk and zero out Q
df_n = df_ydk.copy()
df_n['Q'] = 0.0   # or 0 if you prefer integers

# 3. Replace Q in df_n with Q-cms starting from start_idx
end_idx = start_idx + len(df) - 1

# Optional safety check
if end_idx >= len(df_n):
    raise ValueError(
        f"df_n is too short to insert all Q-cms values "
        f"(need up to index {end_idx}, but length is {len(df_n)})"
    )

df_n.loc[start_idx:end_idx, 'Q'] = -df['Q-cms'].to_numpy()
df_n.to_csv(f'../segment0/geo/lateralDatas/qlat0', sep='\t', header = None, index=False)