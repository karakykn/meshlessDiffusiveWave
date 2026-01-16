
#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

CASE = Path("/mnt/data/lowerMississippi")
GEO = CASE / "segment0" / "geo"
( CASE / "data" ).mkdir(exist_ok=True, parents=True)

USGS_Q_CSV = CASE / "data" / "baton_rouge_discharge.csv"
USACE_STAGE_CSV = CASE / "data" / "donaldsonville_stage.csv"

DT_SEC = 3600
H_OFFSET_M = 18.00

def to_uniform_time(df, col, dt_sec=3600):
    df = df.copy()
    df['datetime'] = pd.to_datetime(df['datetime'], utc=True).dt.tz_convert(None)
    df = df.set_index('datetime').sort_index()
    out = df.resample(f"{int(dt_sec/60)}min").interpolate('time')
    out['t_sec'] = (out.index - out.index[0]).total_seconds()
    out = out.reset_index(drop=False).rename(columns={'index':'datetime'})
    return out[['t_sec', col]]

def main():
    q = pd.read_csv(USGS_Q_CSV)
    q.rename(columns={q.columns[0]:'datetime', q.columns[1]:'discharge_cfs'}, inplace=True)
    qu = to_uniform_time(q, 'discharge_cfs', dt_sec=DT_SEC)
    qu['Q'] = qu['discharge_cfs'] * 0.028316846592
    qu[['t_sec','Q']].to_csv(GEO / "boundary_Q", sep=' ', header=False, index=False)

    s = pd.read_csv(USACE_STAGE_CSV)
    s.rename(columns={s.columns[0]:'datetime', s.columns[1]:'stage_ft'}, inplace=True)
    su = to_uniform_time(s, 'stage_ft', dt_sec=DT_SEC)
    su['h_m'] = su['stage_ft'] * 0.3048 + H_OFFSET_M
    su[['t_sec','h_m']].to_csv(GEO / "boundary_h", sep=' ', header=False, index=False)

    print("Wrote:", GEO / "boundary_Q")
    print("Wrote:", GEO / "boundary_h")

if __name__ == "__main__":
    main()
