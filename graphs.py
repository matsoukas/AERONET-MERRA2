import mylib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import numpy as np
import seaborn as sns
# import netCDF4 as nc
import cartopy.crs as ccrs
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import StrMethodFormatter


# stations = pd.read_csv("allstations.csv", sep=',', parse_dates=[1])
stations = pd.read_pickle("allstations.pkl")

merra = pd.read_csv("MERRA2/MERRA2_550nm.csv", index_col=0)

# # %% Correction with altitude, not used
# # station_names = np.unique(stations.index.get_level_values("Site"))
# station_names = np.unique(stations["Site"])
# ds = nc.Dataset("MERRA2/MERRA2_101.const_2d_asm_Nx.00000000.nc4")
# melev = np.squeeze(ds["PHIS"][:]/9.80665)
# lats = ds["lat"][:]
# lons = ds["lon"][:]
# # plt.imshow(np.flipud(melev))
# # plt.colorbar()
# # plt.show()

# for station_name in station_names:
#     # If stations is read with a Multiindex
#     # currstations = stations.xs(station_name, level="Site")
#     # else
#     currstations = stations[stations["Site"] == station_name]
#     elev = currstations["Elev"].mean()
#     lat = currstations["Lat"].mean()
#     lon = currstations["Lon"].mean()
#     latidx = (np.abs(lats - lat)).argmin()
#     lonidx = (np.abs(lons - lon)).argmin()

#     # If stations is read with a Multiindex
#     # stations.loc[station_name, "Diff"] = melev[latidx, lonidx] - elev
#     # stations.loc[station_name, "CorrTau"] = \
#     #    (currstations["AOD550"] * np.exp((currstations["Diff"])/2000.)).values
#     # else
#     # name, "AOD550"]
#     stations.loc[stations["Site"] == station_name, "Diff"] = \
#         melev[latidx, lonidx] - elev
#     factor = np.exp((melev[latidx, lonidx] - elev)/2100.)
#     stations.loc[stations["Site"] == station_name, "cAOD"] = \
#         currstations["AOD550"] * factor
#     stations.loc[stations["Site"] == station_name, "cAAOD"] = \
#         currstations["AAOD550"] * factor

# del ds, currstations, elev, lat, lon, lats, lons, latidx, lonidx, station_name
# del station_names, melev, factor


# %%  Histograms of MERRA2 and AERONET quantities
fig = plt.figure(figsize=(20, 11))
ax = fig.add_subplot(231)
sns.histplot(stations["AOD550"], log_scale=True, stat="density", color='blue',
             element='step')
sns.histplot(merra["AOD"], log_scale=True, stat="density", color='orange',
             element='step')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([0.03, 0.1, 0.3, 1, 3, 10])
ax.xaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
plt.xlabel("AOD")
plt.text(0.05, 0.95, 'a', transform=ax.transAxes)

ax = fig.add_subplot(232)
sns.histplot(stations["AAOD550"], log_scale=True, stat="density", color='blue',
             element='step')
sns.histplot(merra["AAOD"], log_scale=True, stat="density", color='orange',
             element='step')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([1E-3, 3E-3, 1E-2, 3E-2, 1E-1, 3E-1, 1])
ax.xaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
plt.xlabel("AAOD")
plt.text(0.05, 0.95, 'b', transform=ax.transAxes)

ax = fig.add_subplot(233)
sns.histplot(stations["SSA550"], stat="density", color='blue',
             binrange=(0.75, 1), element='step')
sns.histplot(merra["SSA"], stat="density", color='orange', binrange=(0.75, 1),
             element='step')
plt.xlabel("SSA")
plt.legend(["AERONET", "MERRA-2"])
plt.text(0.05, 0.95, 'c', transform=ax.transAxes)

ax = fig.add_subplot(234)
sns.histplot(merra["AOD"] - stations["AOD550"], binrange=(-1.1, 1.1),
             stat='density', element='step')
plt.axvline(0, color='red', linestyle="dashed")
plt.axvline((merra["AOD"] - stations["AOD550"]).mean(), color='black',
            linestyle="dotted")
plt.xlabel("AOD")
plt.text(0.05, 0.95, 'd', transform=ax.transAxes)

ax = fig.add_subplot(235)
sns.histplot(merra["AAOD"] - stations["AAOD550"], binrange=(-0.15, 0.15),
             stat='density', element='step')
plt.axvline(0, color='red', linestyle="dashed")
plt.axvline((merra["AAOD"] - stations["AAOD550"]).mean(), color='black',
            linestyle="dotted")
plt.xlabel("AAOD")
plt.text(0.05, 0.95, 'e', transform=ax.transAxes)

ax = fig.add_subplot(236)
sns.histplot(merra["SSA"] - stations["SSA550"], binrange=(-0.15, 0.15),
             stat='density', element='step')
plt.axvline(0, color='red', linestyle="dashed")
plt.axvline((merra["SSA"] - stations["SSA550"]).mean(), color='black',
            linestyle="dotted")
plt.xlabel("SSA")
plt.text(0.05, 0.95, 'f', transform=ax.transAxes)


plt.savefig("../histograms.pdf", bbox_inches="tight")
del fig, ax


# %% Cartopy map of the differences between AERONET and MERRA2
lat = stations.groupby('Site')["Lat"].mean()
lon = stations.groupby('Site')["Lon"].mean()
fig = plt.figure(figsize=(20, 5))

stations["AODDiff"] = merra["AOD"] - stations["AOD550"]
diff = stations.groupby('Site')["AODDiff"].mean()
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()
sc = ax.scatter(lon, lat, c=diff, s=10, cmap='seismic',
                transform=ccrs.PlateCarree(),
                norm=colors.SymLogNorm(linthresh=0.01, vmin=diff.min(),
                                       vmax=diff.max()))
# formatter = LogFormatter(10, labelOnlyBase=True, minor_thresholds=(.5,1))
# cb = plt.colorbar(sc, ticks=[.1,.5, 0.6], format=formatter)
plt.colorbar(sc, location="bottom")
plt.title("MERRA-2 - AERONET AOD")
# kw = dict(prop="sizes", num=6, color='r')
# plt.legend(*sc.legend_elements(**kw), loc='lower left', labelcolor='r')

stations["AAODDiff"] = merra["AAOD"] - stations["AAOD550"]
diff = stations.groupby('Site')["AAODDiff"].mean()
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()
sc = ax.scatter(lon, lat, c=diff, s=10, cmap='seismic',
                norm=colors.SymLogNorm(linthresh=0.01, vmin=diff.min(),
                                       vmax=diff.max()),
                transform=ccrs.PlateCarree())
plt.colorbar(sc, location="bottom")
plt.title("MERRA-2 - AERONET AAOD")

stations["SSADiff"] = merra["SSA"] - stations["SSA550"]
diff = stations.groupby('Site')["SSADiff"].mean()
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()
sc = ax.scatter(lon, lat, c=diff, s=10, cmap='seismic',
                norm=colors.SymLogNorm(linthresh=0.01, vmin=diff.min(),
                                       vmax=diff.max()),
                transform=ccrs.PlateCarree())
plt.colorbar(sc, location="bottom")
plt.title("MERRA-2 - AERONET SSA")

plt.savefig("../Diffs_OD_SSA_Clim_map.pdf", bbox_inches="tight")
plt.show()
del diff, ax, lat, lon, fig, sc

# %% Taylor diagrams AERONET-MERRA for AOD, AAOD, SSA
fig = plt.figure(figsize=(20, 6))
varis = ["AOD", "AAOD", "SSA"]
mylib.taylorWL(stations, merra, varis)
plt.savefig("../Taylor.pdf", bbox_inches="tight")

del fig, varis


# %% AOD Taylor diagrams for every AerType
fig = plt.figure(figsize=(30, 10))
ax = fig.add_subplot(131)
mylib.taylor_by_attribute(stations, "AOD550", merra, "AOD", "AerType")

ax = fig.add_subplot(132)
mylib.taylor_by_attribute(stations, "AAOD550", merra, "AAOD", "AerType")

ax = fig.add_subplot(133)
mylib.taylor_by_attribute(stations, "SSA550", merra, "SSA", "AerType")

# maybe larger font in the legend?
plt.savefig("../TaylorByType.pdf", bbox_inches="tight")
del fig, ax


# %% Scatterplots with MERRA2 on y and AERONET on x

# AOD
fig = plt.figure(figsize=(15, 5))
ax = fig.add_subplot(1, 3, 1)
mylib.scatter(stations, "AOD550", merra, "AOD", ax)

# AAOD
ax = fig.add_subplot(1, 3, 2)
mylib.scatter(stations, "AAOD550", merra, "AAOD", ax)

ax = fig.add_subplot(1, 3, 3)
mylib.scatter(stations, "SSA550", merra, "SSA", ax)

plt.savefig("../KDE_Scatter.pdf", bbox_inches="tight")
del fig, ax


# %% Koppen climate types

# If we use all 3 letters of the Climate we have 27 climate types
# If we use the first 2 letters of Climate we have 15 climate types
# If we use only the 1st letter of the Climate we have 5 climate types
# count_Climate = stations.Climate.unique().size

# Create a column only with the main climate type
stations['ClimateZone'] = stations['Climate'].str[:1]
# Separate Central Asian Deserts from type B to type B_As(ia)
mask = (stations["ClimateZone"] == "B") & (stations["Lat"] >= 35) & \
        (stations["Lon"] >= 47)
stations.loc[mask, "ClimateZone"] = "B_As"

# Create taylor diagrams for AOD 440-1020: each point represents a
# Koppen climate type
fig = plt.figure(figsize=(22, 7))
ax = fig.add_subplot(131)
mylib.taylor_by_attribute(stations, "AOD550", merra, "AOD", "ClimateZone")

ax = fig.add_subplot(132)
mylib.taylor_by_attribute(stations, "AAOD550", merra, "AAOD", "ClimateZone")

ax = fig.add_subplot(133)
mylib.taylor_by_attribute(stations, "SSA550", merra, "SSA", "ClimateZone")

plt.savefig("../TaylorByClimate.pdf", bbox_inches="tight")

del fig, ax


# %% Climate zones and aerosol types
# Creates a heatmap, showing climate zones on the vertical and percentages
# of aerosol types for each zone in the horizontal
aertypes = np.unique(stations["AerType"])
climes = np.unique(stations["ClimateZone"])
perc_map = pd.DataFrame(index=climes, columns=aertypes, dtype=float)
for clime in climes:
    stat = stations[stations["ClimateZone"] == clime]
    res = stat["AerType"].value_counts(normalize=True, sort=False)*100
    perc_map.loc[clime] = res

plt.figure()
sns.heatmap(perc_map, annot=True, fmt=".0f", cbar=False, cmap='viridis')
_ = plt.yticks(rotation=0)
plt.savefig("../ClimatesHeatmap1.pdf", bbox_inches="tight")

perc_map = pd.DataFrame(index=aertypes, columns=climes, dtype=float)
for aertype in aertypes:
    stat = stations[stations["AerType"] == aertype]
    res = stat["ClimateZone"].value_counts(normalize=True, sort=False)*100
    perc_map.loc[aertype] = res

plt.figure()
sns.heatmap(perc_map, annot=True, fmt=".0f", cbar=False, cmap='viridis')
_ = plt.yticks(rotation=0)
plt.savefig("../ClimatesHeatmap2.pdf", bbox_inches="tight")

del res, stat, perc_map, climes, aertypes, clime, aertype


# %% Assign MERRA2 dominant types
# Provides the dominant aerosol type for merra, by selecting the type with
# the maximum AOD among Dust, Sea salt, Sulfates, BC, OC
merratype = merra[["DUAOD", "SSAOD", "SUAOD", "BCAOD", "OCAOD"]].idxmax(axis=1)
merratype = merratype.str[:2]
merra.insert(loc=17, column="AerType", value=merratype)
aertypes = np.unique(stations["AerType"])
merratypes = np.unique(merra["AerType"])

# BC is dominant at very few places
pd.value_counts(merra["AerType"])

perc_map = pd.DataFrame(index=aertypes, columns=merratypes, dtype=float)
for aertype in aertypes:
    stat = merra[stations["AerType"] == aertype]
    res = stat["AerType"].value_counts(normalize=True, sort=False)*100
    perc_map.loc[aertype] = res

plt.figure()
sns.heatmap(perc_map, annot=True, fmt=".0f", cbar=False, cmap='viridis')
_ = plt.yticks(rotation=0)
plt.savefig("../TypesHeatmap1.pdf", bbox_inches="tight")


perc_map = pd.DataFrame(index=merratypes, columns=aertypes, dtype=float)
for aertype in merratypes:
    stat = stations[merra["AerType"] == aertype]
    res = stat["AerType"].value_counts(normalize=True, sort=False)*100
    perc_map.loc[aertype] = res

plt.figure()
sns.heatmap(perc_map, annot=True, fmt=".0f", cbar=False, cmap='viridis')
_ = plt.yticks(rotation=0)
plt.savefig("../TypesHeatmap2.pdf", bbox_inches="tight")

del res, stat, perc_map, aertypes, merratypes, aertype, merratype


# %% Create a taylor diagram for AOD550 for every continent
fig = plt.figure(figsize=(30, 10))
ax = fig.add_subplot(131)
mylib.taylor_by_attribute(stations, "AOD550", merra, "AOD", "Continent")

ax = fig.add_subplot(132)
mylib.taylor_by_attribute(stations, "AAOD550", merra, "AAOD", "Continent")

ax = fig.add_subplot(133)
mylib.taylor_by_attribute(stations, "SSA550", merra, "SSA", "Continent")
plt.savefig("../TaylorByContinent.pdf", bbox_inches="tight")


# %% Continent Heatmap
aertypes = np.unique(stations["AerType"])
continents = np.unique(stations["Continent"])
perc_map = pd.DataFrame(index=continents, columns=aertypes, dtype=float)
for continent in continents:
    stat = stations[stations["Continent"] == continent]
    res = stat["AerType"].value_counts(normalize=True, sort=False)*100
    perc_map.loc[continent] = res

plt.figure()
sns.heatmap(perc_map, annot=True, fmt=".0f", cbar=False, cmap='viridis')
_ = plt.yticks(rotation=0)
plt.savefig("../ContinentssHeatmap1.pdf", bbox_inches="tight")

# %% Creates AAOD and SSA histograms for every cappa type
# AAOD
fig = plt.figure(figsize=(20, 8))
aertypes = np.unique(stations["AerType"])
for aertype in aertypes:
    index = list(aertypes).index(aertype)+1
    ax = fig.add_subplot(2, 4, index)
    sns.histplot(stations[stations["AerType"] == aertype]["AAOD550"],
                 log_scale=True, stat="density", color='blue',
                 element='step')
    sns.histplot(merra[stations["AerType"] == aertype]["AAOD"],
                 log_scale=True, stat="density", color='orange',
                 element='step')
    plt.xlabel(aertype)
plt.legend(["AERONET", "MERRA-2"])

# SSA
fig = plt.figure(figsize=(20, 8))
aertypes = np.unique(stations["AerType"])
for aertype in aertypes:
    index = list(aertypes).index(aertype)+1
    ax = fig.add_subplot(2, 4, index)
    sns.histplot(stations[stations["AerType"] == aertype]["SSA550"],
                 stat="density", color='blue',
                 element='step')
    sns.histplot(merra[stations["AerType"] == aertype]["SSA"],
                 stat="density", color='orange',
                 element='step')
    plt.xlabel(aertype)
plt.legend(["AERONET", "MERRA-2"])
