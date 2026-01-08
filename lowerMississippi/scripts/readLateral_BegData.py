import pandas as pd

# Read the data file
lat1 = pd.read_csv(
    "../../../Junk/MESH_Code_irregular-master/Mississippi_River_11_years_20200511/lateral_in_min_diffusive/lateral_0179.txt",
    sep="\t",
    names=["time-m", "qlat"]
)


# Convert the "time-m" column (minutes since start) into datetime objects
start_date = pd.Timestamp("2009-01-01")
lat1["datetime"] = start_date + pd.to_timedelta(lat1["time-m"], unit="m")
lat1['qlat'] = lat1['qlat']

# Filter for data within the year 2011
lat1_2011 = lat1[lat1["datetime"].dt.year == 2011]
start_date = pd.Timestamp("2011-01-01")
lat1_2011["time-s"] = (lat1_2011["datetime"] - start_date).dt.total_seconds()

lat1_2011[["time-s", "qlat"]].to_csv("../segment0/geo/qlat2", sep="\t", index=False, header=False)

# Optionally, view or save
print(lat1_2011.head())
# lat1_2011.to_csv("lateral_2011.csv", index=False)
