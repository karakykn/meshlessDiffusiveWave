import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

draw = [0,1,6]
color = ['k','b','r','g']

inflow = pd.read_csv(f'../data/tarbertQ.csv')

# plt.plot(inflow['date'], inflow['cms'])

for i in draw:
    fPath = f'../segment0/geo/lateralDatas/qlat{i}'
    df = pd.read_csv(fPath, sep='\t', names=['t', 'Q'])
    start = pd.Timestamp('2011-01-01 00:00:00')
    start_2011 = pd.Timestamp('2011-01-01 00:00:00')
    start_2012 = pd.Timestamp('2012-01-01 00:00:00')
    df['date'] = start + pd.to_timedelta(df['t'], unit='s')
    df_2011 = df[(df['date'] >= start_2011) & (df['date'] <= start_2012)]
    plt.plot(df_2011['date'], -df_2011['Q'], label=f'{i}')

plt.plot(df['date'], inflow['cms']/10)
plt.legend()
plt.show()

# # fPath = f'../segment0/geo/lateralDatas/qlat1'
# # df = pd.read_csv(fPath, sep='\t', names=['t', 'Q'])
# start = pd.Timestamp('2011-01-01 00:00:00')
# df['date'] = start + pd.to_timedelta(df['t'], unit='s')
# # print(df[df['Q'] < -8890])
# # plt.plot(df['date'], -df['Q'])
# # plt.legend()
# # plt.show()
#
# df_bohemia = df.copy()
# df_bohemia['Q'] = 0
# peak_bohemia = 50000 * 0.02831683199881
# span_sec = 3 * 7 * 24 * 60 * 60
# peak_date = pd.Timestamp('2011-05-18 00:00:00')
# start_Date = pd.Timestamp('2011-01-01 00:00:00')
# peak_sec = (peak_date - start_Date).total_seconds()
#
# sigmaC = .5e-1
# # sigma = span_sec * sigmaC
# # Upstream side (x <= x0): decays quickly
# df_bohemia['Q'] = -np.exp(-0.5 * ((df_bohemia['t'] - peak_sec) / (df_bohemia['t'].iloc[-1] * sigmaC)) ** 2) * peak_bohemia
# df_bohemia.loc[(~df_bohemia['date'].between(pd.Timestamp('2011-05-01 00:00:00'), pd.Timestamp('2011-06-07 00:00:00'))), 'Q'] = 0
# df_bohemia.to_csv(f'../segment0/geo/lateralDatas/qlat6', sep='\t', header=None, columns=['t', 'Q'], index=False)
#
#
# plt.plot(df_bohemia['date'], df_bohemia['Q'])
# plt.show()

# fPath = f'../segment0/geo/lateralDatas/qlat1'
# df = pd.read_csv(fPath, sep='\t', names=['t', 'Q'])
# df['Q'] = df['Q'] / 2
# df.to_csv(f'../segment0/geo/lateralDatas/qlat1',header=None, index=False, sep='\t', )