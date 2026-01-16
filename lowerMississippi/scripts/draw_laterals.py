import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 12,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "legend.fontsize": 11,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "lines.linewidth": 1.5
})
draw = [0, 1, 6]
color = ['b','r','orange']
inflow = pd.read_csv('../segment0/geo/boundary_Q', names = ['seconds', 'cms'], sep='\t')

fig, ax1 = plt.subplots(figsize=(7, 4))

start = pd.Timestamp('2011-01-01 00:00:00')
start_2011 = pd.Timestamp('2011-01-01 00:00:00')
start_2012 = pd.Timestamp('2012-01-01 00:00:00')

inflow['Date'] = pd.to_datetime(start) + pd.to_timedelta(inflow.iloc[:,0].astype(float), unit="s")

# ---- Left axis: inflow ----
ax1.plot(inflow['Date'], inflow['cms'] ,
         color='k', label='Tarbert Landing', linewidth=1.5)
ax1.set_ylabel('Inflow (cms)')
ax1.set_xlabel('Date')

# ---- Right axis: Q from loop ----
ax2 = ax1.twinx()
ax2.set_ylabel('Lateral Outflow at Spillways (cms)')

lines = []
labels = ['Morganza Spillway', 'Bonnet Carre Spillway', 'Bohemia Spillway']

for j, i in enumerate(draw):
    fPath = f'../segment0/geo/lateralDatas/qlat{i}'
    df = pd.read_csv(fPath, sep='\t', names=['t', 'Q'])

    df['date'] = start + pd.to_timedelta(df['t'], unit='s')
    df_2011 = df[(df['date'] >= start_2011) & (df['date'] <= start_2012)]

    line, = ax2.plot(df_2011['date'], -df_2011['Q'], label=labels[j], color = color[j], linewidth=1.5)
    lines.append(line)
    labels.append(f'Qlat {i}')

# ---- Combine legends from both axes ----
lines1, labels1 = ax1.get_legend_handles_labels()
ax1.legend(lines1 + lines, labels1 + labels, loc='upper right')

xmin = df['date'].iloc[0]
xmax = df['date'].iloc[-1]
ax1.set_xlim(xmin, xmax)
ax1.xaxis.set_major_locator(
    mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11])
)
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
plt.tight_layout()
plt.savefig('ex3_lats.pdf')
plt.show()
