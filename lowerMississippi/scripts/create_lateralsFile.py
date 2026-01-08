import pandas as pd
import numpy as np

filePath = f'../../../Junk/MESH_Code_irregular-master/Mississippi_River_11_years_20200511/dx.txt'
df = pd.read_csv(filePath, sep='\t', names=['x', 'dx'])
presentL = 510897
begL = 397113.28
rat = presentL / begL
# rat = 1
df['x_rev'] = df['x'].iloc[-1] - df['x']
lat_idx = [128,139,179,196,211,240,250,252,256,258,261,265,267,269,271,274,276]
lat_p = [2,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1]
# lat_idx = [lat_idx[10]]
begLat = np.loadtxt(f'../../../Junk/MESH_Code_irregular-master/Mississippi_River_11_years_20200511/dx.txt')
n = len(lat_idx)
# print(df['x'].iloc[lat_idx] * rat)
names = np.array([])
qfilee = np.array([])
kernel = np.array([])
spread = begLat[lat_idx-np.ones_like(lat_idx),1] * 2
x = (df['x'].iloc[lat_idx] * rat).to_numpy()
for i in range(n):
    names = np.append(names, f'lateral_pos{i}')
    qfilee = np.append(qfilee, f'qlat{i}')
    kernel = np.append(kernel, f'tri')

laters = {'names': names,
          'x': x,
          'spread_m': spread,
          'kernel': kernel,
          'q_file': qfilee,
          'latPoints': lat_p
          }
later = pd.DataFrame(laters)
later.to_csv(f'../segment0/geo/laterals_n', index=False, float_format="%.0f")