import os
import pandas as pd
import mylib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import haversine as hs
from matplotlib.ticker import StrMethodFormatter

# %% read the AERONET data and create the DateFrame 'stations'
file_dir = "/home/christos/Data/AERONET/INV/LEV20/ALL/ALL_POINTS"
for _, _, files in os.walk(file_dir, topdown=True):
    pass  # Just populate files

stations = pd.DataFrame()
for file in files:
    print(file)
    station = mylib.aeronet_read(file_dir, file)
    stations = pd.concat([stations, station])
del file, station, files, file_dir

# Index is a little messed up, let's reset it to sequential values
stations = stations.reset_index(drop=True)

# Use more user-friendly column names
stations = mylib.rename_columns(stations)

# Check if the calculated SSA is futher from the given SSA by an error margin
# mylib.check_ssa(stations, error_margin=2.4E-4)

# Check if the calculated AExp is futher from the given AE by an error margin
# mylib.check_ae(stations, error_margin=4.9E-2)

# Check if the calculated AAExp is futher from the given AAE by an error margin
# mylib.check_aae(stations, error_margin=3.9E-1)

# Uses the Angstrom Exponent to interpolate the AOD measurements to the
# standard wavelengths, 440, 675, 870 nm
mylib.fill_AOD(stations)

# Same thing for the AAOD
mylib.fill_AAOD(stations)

# Now fill SSA at the standard wavelengths
# mylib.fill_SSA(stations)

# %% Create the 550 nm data for AOD and AAOD

# Slow method for AOD and AAOD, using 2nd degree polynomial fit
stations.insert(loc=4, column="AOD550", value=np.nan)
stations.insert(loc=10, column="AAOD550", value=np.nan)
x = np.array(np.log([440, 675, 870, 1020])).T
for i in range(0, len(stations)):
    y = stations[["AOD440", "AOD675", "AOD870", "AOD1020"]].iloc[i]
    # polynomial fit with degree = 2
    model = np.poly1d(np.polyfit(x, y, 2))
    stations.at[stations.index[i], "AOD550"] = model(np.log(550))

    # Does the polynomial fit hold for AAOD?
    y = stations[["AAOD440", "AAOD675", "AAOD870", "AAOD1020"]].iloc[i]
    model = np.poly1d(np.polyfit(x, y, 2))
    stations.at[stations.index[i], "AAOD550"] = model(np.log(550))

# # Fast, vectorized method for AOD and AAOD using AE
# AOD550 = stations["AOD440"]*(440./550)**stations["AE440-870"]
# stations.insert(loc=2, column="AOD550_1", value=AOD550)
# AAOD550 = stations["AAOD440"]*(440./550)**stations["AAE440-870"]
# stations.insert(loc=12, column="AAOD550_1", value=AAOD550)

# Now the SSA at 550
SSA550 = (1. - stations["AAOD550"]/stations["AOD550"])
stations.insert(loc=15, column="SSA550", value=SSA550)

del i, x, y, model, SSA550

# %% Estimate SAE450-550 and AAE532-660 for Cappa et al. (2016)
# Using 440 - 675 for both cases, it doesn't make sense to interpolate

SAOD440 = stations["AOD440"] - stations["AAOD440"]
SAOD675 = stations["AOD675"] - stations["AAOD675"]
stations["SAE440-675"] = np.log(SAOD440 / SAOD675) / np.log(675./440.)
del SAOD440, SAOD675
stations["AAE440-675"] = np.log(stations["AAOD440"] / stations["AAOD675"]) \
    / np.log(675./440.)

# Only for checking the agreement between the 2nd degree polynomial and
# AE interpolations

# stations["AE440-675"] = np.log(stations["AOD440"] / stations["AOD675"]) \
#     / np.log(675./440.)
# AOD550 = stations["AOD440"]*(440./550)**stations["AE440-675"]
# stations.insert(loc=2, column="AOD550_2", value=AOD550)
# AAOD550 = stations["AAOD440"]*(440./550)**stations["AAE440-675"]
# stations.insert(loc=12, column="AAOD550_2", value=AAOD550)

stations["AerType"] = mylib.cappa(stations["SAE440-675"],
                                  stations["AAE440-675"])
plt.savefig("../cappa.pdf", bbox_inches="tight")


# %% Koppen Climate type Chen & Chen
# http://hanschen.org/koppen#permission-and-copyright

stations.insert(column="Climate", loc=27, value=np.nan)
climate = pd.read_csv("Koeppen/koppen_1901-2010.tsv", sep='\t')
stcoor = np.unique(stations[["Lat", "Lon"]], axis=0)
for i in range(0, len(stcoor)):
    distance = hs.haversine_vector(stcoor[i],
                                   climate[["latitude", "longitude"]],
                                   comb=True)
    j = distance.argmin()
    ktype = climate.iloc[j]["p1901_2010"]
    stations.loc[(stations["Lon"] == stcoor[i, 1]) &
                 (stations["Lat"] == stcoor[i, 0]), "Climate"] = ktype

del i, j, distance, climate, ktype, stcoor


# %% Cartopy map of the station locations and number of observations
coords = np.unique(stations[["Lat", "Lon"]], axis=0, return_counts=True)

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

ax.set_global()
ax.stock_img()
ax.coastlines()

# ax.plot(coords[0][:, 1], coords[0][:, 0], 'ro', markersize=coords[1],
#         transform=ccrs.PlateCarree())
sc = ax.scatter(coords[0][:, 1], coords[0][:, 0], s=coords[1]/40, color='r',
                transform=ccrs.PlateCarree())

kw = dict(prop="sizes", num=6, color='r', func=lambda s: s*40)
plt.legend(*sc.legend_elements(**kw), loc='lower left', labelcolor='r')
plt.savefig("../Locations_map.pdf", bbox_inches="tight")
plt.show()
del kw, sc, coords, ax, fig


# %% Cartopy map of the station locations, with AOD, AAOD, SSA
lat = stations.groupby('Site')["Lat"].mean()
lon = stations.groupby('Site')["Lon"].mean()
fig = plt.figure(figsize=(20, 5))

aod550 = stations.groupby('Site')["AOD550"].mean()
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson())
ax.set_global()
# ax.stock_img()
ax.coastlines()
# ax.plot(coords[0][:, 1], coords[0][:, 0], 'ro', markersize=coords[1],
#         transform=ccrs.PlateCarree())
sc = ax.scatter(lon, lat, c=aod550, s=10,
                norm=colors.LogNorm(aod550.min(), vmax=1.1),  # aod550.max()),
                cmap='viridis',
                transform=ccrs.PlateCarree())
plt.colorbar(sc, location="bottom", extend="max",
             ticks=np.linspace(0.3, 1, 8),
             format=StrMethodFormatter("{x:g}"))
plt.title("AERONET AOD")


aaod550 = stations.groupby('Site')["AAOD550"].mean()
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()
sc = ax.scatter(lon, lat, c=aaod550, s=10,
                norm=colors.LogNorm(vmin=aaod550.min(), vmax=aaod550.max()),
                cmap='viridis',
                transform=ccrs.PlateCarree())
plt.colorbar(sc, location="bottom",
             ticks=[0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1],
             format=StrMethodFormatter("{x:g}"))  # format="%.1f")
plt.title("AERONET AAOD")


ssa550 = stations.groupby('Site')["SSA550"].mean()
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()
sc = ax.scatter(lon, lat, c=ssa550, s=10,
                # norm=colors.LogNorm(vmin=ssa550.min(), vmax=ssa550.max()),
                cmap='viridis',
                transform=ccrs.PlateCarree())
plt.colorbar(sc, location="bottom", format="%.2f")
plt.title("AERONET SSA")
plt.savefig("../OD_SSA_Clim_map.pdf", bbox_inches="tight")
plt.show()

del aod550, aaod550, ssa550, sc, ax, fig, lat, lon


# %% Assign continent to each station
continents = pd.read_csv("sites_and_continents.csv", sep=',',
                         index_col=(0))
dictio = continents.set_index('Site')['continent'].to_dict()
# use the dictionary to add the continent on stations df
stations['Continent'] = stations['Site'].map(dictio)
del dictio, continents

# examine the stations that are categorized as nan
# stations[stations['Continent'].isna()]
# The following four stations are set to no continent, manual assignment

# Alboran :Gulf of gibraltar, between Europe and Africa.
stations.loc[stations.Site == "Alboran", "Continent"] = "Africa"
# Abu_Al_Bukhoosh : Persian Gulf : Asia
stations.loc[stations.Site == "Abu_Al_Bukhoosh", "Continent"] = "Asia"
# COVE: North America : US
stations.loc[stations.Site == "COVE", "Continent"] = "North America"
# Socheongcho : Asia : nan
stations.loc[stations.Site == "Socheongcho", "Continent"] = "Asia"
# COVE_SEAPRISM : North America : US
stations.loc[stations.Site == "COVE_SEAPRISM", "Continent"] = "North America"

# count how many measurements exist per continent
# num_of_measurements_per_continent = \
#     stations['Continent'].value_counts(dropna=False)
# print(num_of_measurements_per_continent)


# %% Cleanup and save allstations

stations.drop(columns=["AOD440", "AOD675", "AOD870", "AOD1020", "AE440-870",
                       # "SSA440", "SSA675", "SSA870", "SSA1020",
                       "AAOD440", "AAOD675", "AAOD870", "AAOD1020",
                       "AOD443", "AOD667", "AOD865",
                       # "SSA443", "SSA667", "SSA865",
                       "AAOD443", "AAOD667", "AAOD865"], inplace=True)
stations.to_pickle("allstations.pkl")
