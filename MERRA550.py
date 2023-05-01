import netCDF4 as nc
import pandas as pd
import numpy as np
import glob
import datetime
# import os

stations = pd.read_pickle("allstations.pkl")
# stations = pd.read_csv("allstations.csv", sep=',', parse_dates=[1])

# Some times the AERONET measurements are after 23:30, so they should be
# matched with MERRA2 files for the next day.
stations.loc[stations["Time"].dt.time > datetime.time(23, 30), "Time"] \
    += datetime.timedelta(minutes=30)

# MERRA2 files are one per each day. To minimize reading/writing overhead
# we will work day by day instead of station by station, so get all dates
dates = np.unique(stations["Time"].dt.date)

# Empty dataframe which will hold the data later
merra = pd.DataFrame(index=range(0, len(stations)),
                     columns=["Lat", "Lon", "Time", "AOD", "SAOD",
                              "DUAOD", "DUSAOD", "BCAOD", "BCSAOD",
                              "SSAOD", "SSSAOD", "SUAOD", "SUSAOD",
                              "OCAOD", "OCSAOD"])

# Just a random nc file, to read the latitudes, longitudes, hours
# ds = nc.Dataset("/home/christos/Data/MERRA2/" +
#                 "MERRA2_200.tavg1_2d_aer_Nx.19930617.nc4.nc4")
# lats = ds["lat"][:]  # Every 0.5 degrees
# lons = ds["lon"][:]  # Every 0.625 degrees
# times = ds["time"][:]  # 0, 1, 2, ..., 23 hours

for curr in dates:  # Current date
    print(curr)
    currstations = stations[stations["Time"].dt.date == curr]
    datestr = curr.strftime("%Y%m%d")

    # The filenames in MERRA2 are not uniform, so we have to search for the
    # correct day
    fn = glob.glob("/media/hdisk/MERRA2/MERRA2_*"+datestr+"*")
    if len(fn) != 1:
        raise Exception("More than one filename matched date")
    fn = fn[0]
    ds = nc.Dataset(fn)
    # The following 3 lines could be read out of the loop, but the nc files
    # have different grids
    lats = ds["lat"][:]  # Every 0.5 degrees
    lons = ds["lon"][:]  # Every 0.625 degrees
    times = ds["time"][:]  # 0, 1, 2, ..., 23 hours

    # For the current date, iterate among the stations that have measurements
    for idx, station in currstations.iterrows():
        time = station["Time"].hour*60. + station["Time"].minute + \
            station["Time"].second/60.  # In minute of the day

        latidx = (np.abs(lats - station["Lat"])).argmin()
        if (latidx == 0) or (latidx == len(lats)-1):
            print(station["Site"], station["Time"], "latitude close to edges")
        lonidx = (np.abs(lons - station["Lon"])).argmin()
        if (lonidx == 0) or (lonidx == len(lons)-1):
            print(station["Site"], station["Time"], "longitude close to edges")
        timeidx = (np.abs(times - time)).argmin()
#         if (timeidx == 0) or (timeidx == len(times)-1):
#             print(station["Site"], station["Time"], "time close to edges")

        merra.loc[idx]["Lat"] = lats[latidx]
        merra.loc[idx]["Lon"] = lons[lonidx]
        merra.loc[idx]["Time"] = times[timeidx]
        merra.loc[idx]["AOD"] = ds["TOTEXTTAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["SAOD"] = ds["TOTSCATAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["DUAOD"] = ds["DUEXTTAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["DUSAOD"] = ds["DUSCATAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["SSAOD"] = ds["SSEXTTAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["SSSAOD"] = ds["SSSCATAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["SUAOD"] = ds["SUEXTTAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["SUSAOD"] = ds["SUSCATAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["BCAOD"] = ds["BCEXTTAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["BCSAOD"] = ds["BCSCATAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["OCAOD"] = ds["OCEXTTAU"][timeidx, latidx, lonidx]
        merra.loc[idx]["OCSAOD"] = ds["OCSCATAU"][timeidx, latidx, lonidx]

merra.insert(loc=5, column="AAOD", value=merra["AOD"] - merra["SAOD"])
merra.insert(loc=6, column="SSA", value=merra["SAOD"]/merra["AOD"])
del lats, lons, lonidx, latidx, station, time, times, timeidx, idx, ds, fn
del dates, curr, currstations, datestr

merra.to_csv("MERRA2/MERRA2_550nm.csv", sep=",")
# merra.to_pickle("MERRA2/MERRA2_550nm.pkl")

# %% Code snippet to remove the unecessarily downloaded nc files
# dt_list = pd.date_range(start='1993-06-17', end='2009-12-31', freq='D')
# for date in dt_list:
#     if date not in dates:
#         datestr = date.strftime("%Y%m%d")
#         fn = glob.glob("/media/hdisk/MERRA2/MERRA2_*"+datestr+"*")
#         os.remove(fn[0])
