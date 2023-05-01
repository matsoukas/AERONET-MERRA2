import pandas as pd
import numpy as np
from datetime import datetime
import os
from scipy.stats import linregress, gaussian_kde
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
import skill_metrics as sm


# %% Used in read_AERONET.py
def dateparse_aeronet(x):
    y = datetime.strptime(x, "%d:%m:%Y %H:%M:%S")
    return y


def aeronet_read(file_dir, file):
    # dateparse = lambda x: datetime.strptime(x, "%d:%m:%Y %H:%M:%S")
    aeronet = pd.read_csv(os.path.join(file_dir, file), skiprows=6,
                          usecols=[0, 1, 2, 3,
                                   10, 11, 12, 13,  # Extinction AOD
                                   22,              # Extinction Ang Exponent
                                   # 23, 24, 25, 26,  # SSA
                                   27, 28, 29, 30,  # AAOD
                                   31,              # Absorption Ang Exponent
                                   # 40, 41, 42, 43,  # Asymmetry parameter
                                   # 149,             # Instrument
                                   150, 151           # Lat, Lon
                                   # 152              # Elevation
                                   # 153, 154         # Qual. Level, Scan Type
                                   ],
                          na_values=[-999],
                          parse_dates={'Time': [1, 2]},
                          date_parser=dateparse_aeronet,
                          encoding='latin_1')
    # aeronet = aeronet.set_index(['Site', 'Time'])
    # Use the xs() method in DataFrames for selection using multiindices
    aeronet = aeronet.dropna(axis=0, how='any')
    return aeronet


def rename_columns(df):
    orig = list(df.columns)
    new = orig[:]
    new[2] = new[2].replace("Day_of_Year", "DOY")
    new = [w.replace("AOD_Extinction-Total", "AOD") for w in new]
    new = [w.replace("Absorption_AOD", "AAOD") for w in new]
    # new = [w.replace("Single_Scattering_Albedo", "SSA") for w in new]
    # new[7] = new[7].replace("Extinction_Angstrom_Exponent_440-870nm-Total",
    #                         "AE440-870")
    # new[16] = new[16].replace("Absorption_Angstrom_Exponent_440-870nm",
    #                           "AAE440-870")
    # new[17] = new[17].replace("Instrument_Number", "Instr")
    # new[18] = new[18].replace("Latitude(Degrees)", "Lat")
    # new[19] = new[19].replace("Longitude(Degrees)", "Lon")
    # new[20] = new[20].replace("Elevation(m)", "Elev")
    new[7] = new[7].replace("Extinction_Angstrom_Exponent_440-870nm-Total",
                            "AE440-870")
    new[12] = new[12].replace("Absorption_Angstrom_Exponent_440-870nm",
                              "AAE440-870")
    new[13] = new[13].replace("Latitude(Degrees)", "Lat")
    new[14] = new[14].replace("Longitude(Degrees)", "Lon")
    new = [w.replace("nm", "") for w in new]
    new = [w.replace("[", "") for w in new]
    new = [w.replace("]", "") for w in new]

    mapper = dict(zip(orig, new))

    df = df.rename(mapper, axis='columns')
    return df


def check_ssa(df, error_margin):
    wls = ["440", "675", "870", "1020", "443", "667", "865"]
    for wl in wls:
        mySSA = (df["AOD"+wl] - df["AAOD"+wl]) / df["AOD"+wl]
        diffs = abs(mySSA - df["SSA"+wl]) > error_margin
        if np.count_nonzero(diffs) > 0:
            raise Exception("SSA disagreement for " + wl + "nm")
    return None


def check_ae(df, error_margin):
    myAE = np.log(df["AOD870"]/df["AOD440"]) / np.log(440/870)
    diffs = abs(myAE - df["AE440-870"]) > error_margin
    if np.count_nonzero(diffs) > 0:
        raise Exception("AE disagreement")
    return None


def check_aae(df, error_margin):
    myAAE = np.log(df["AAOD870"]/df["AAOD440"]) / np.log(440/870)
    diffs = abs(myAAE - df["AAE440-870"]) > error_margin
    if np.count_nonzero(diffs) > 0:
        raise Exception("AAE disagreement")
    return None


def fill_AOD(df):
    wl = {"440": "443", "675": "667", "870": "865"}
    for wl1 in wl:
        mask = df["AOD"+wl1].isnull()
        df.loc[mask, "AOD"+wl1] = df.loc[mask, "AOD"+wl[wl1]] * \
            (float(wl1)/float(wl[wl1]))**-df.loc[mask, "AE440-870"]
    return None


def fill_AAOD(df):
    wl = {"440": "443", "675": "667", "870": "865"}
    for wl1 in wl:
        mask = df["AAOD"+wl1].isnull()
        df.loc[mask, "AAOD"+wl1] = df.loc[mask, "AAOD"+wl[wl1]] * \
            (float(wl1)/float(wl[wl1]))**-df.loc[mask, "AAE440-870"]
    return None


def fill_SSA(df):
    wl = ["440", "675", "870"]
    for wl1 in wl:
        mask = df["SSA"+wl1].isnull()
        df.loc[mask, "SSA"+wl1] = 1 - df.loc[mask, "AAOD"+wl1] \
            / df.loc[mask, "AOD"+wl1]
    return None


# %% Used in read_MERRA.py
# def dateparse_merra(x):
#     y = datetime.strptime(x, "%d %m %Y %H")
#     return y


# def AE_interp(merra_bef, merra_aft, var, start, end, middle):
#     """
#     A function that interpolates spectrally MERRA2 data and creates a new
#     column with the wavelength of middle
#     Parameters
#     ----------
#     merra_bef : DataFrame
#         The MERRA2 values before the AERONET measurement time
#     merra_aft : DataFrame
#         The MERRA2 values after the AERONET measurement time
#     var : string
#         The name of the quantity in merra_bef or _aft, e.g. "SAOD"
#     start : int
#         The smaller wavelength of MERRA2 closest to the AERONET wl
#     end : int
#         The larger wavelength of MERRA2 closest to the AERONET wl
#     middle : int
#         The wl, at which we want to spectrally interpolate the
#         MERRA2 data
#     Returns
#     -------
#     Nothing, it modifies merra_bef and merra_aft
#     """
#     AE_bef = np.log(merra_bef[var+str(start)] / merra_bef[var+str(end)]) / \
#         np.log(end/start)
#     AE_aft = np.log(merra_aft[var+str(start)] / merra_aft[var+str(end)]) / \
#         np.log(end/start)

#     merra_bef[var+str(middle)] = merra_bef[var+str(start)] * \
#         (start/middle)**AE_bef
#     merra_aft[var+str(middle)] = merra_aft[var+str(start)] * \
#         (start/middle)**AE_aft
#     return None


# def select3h(stations, merra_bef, merra_aft):
#     """
#     Function selecting the 3h value from MERRA closest to the AERONET time
#     for all relevant optical properties
#     Parameters
#     ----------
#     stations : DataFrame
#         All AERONET data
#     merra_bef : DataFrame
#         The MERRA2 values before the AERONET measurement time
#     merra_aft : DataFrame
#         The MERRA2 values after the AERONET measurement time
#     Returns
#     -------
#     merra : Similar to merra_bef and merra_aft, but for the correct 3hourly
#         value
#     """
#     qus = ['AOD440', 'AOD675', 'AOD870', 'AOD1020',
#            'SAOD440', 'SAOD675', 'SAOD870', 'SAOD1020',
#            'AAOD440', 'AAOD675', 'AAOD870', 'AAOD1020',
#            'SSA440', 'SSA675', 'SSA870', 'SSA1020',
#            'SAE450-650', 'AAE532-660']

#     merra = pd.DataFrame(index=range(0, len(stations)),
#                          columns=qus)

#     tdiff_bef = stations["Time"] - merra_bef["Time"]
#     tdiff_bef.rename("Before", inplace=True)
#     tdiff_aft = merra_aft["Time"] - stations["Time"]
#     tdiff_aft.rename("After", inplace=True)
#     # Create tdiffs a DataFrame with columns "Before" and "After"
#     # tdiffs holds the time differences between AERONET and MERRA
#     tdiffs = pd.concat([tdiff_bef, tdiff_aft], axis=1)
#     for qu in qus:
#         # If the MERRA2 "before" value is closer to the AERONET time keep the
#         # "before" value, else keep the "after" value
#         merra[qu] = np.where(tdiffs["Before"] < tdiffs["After"],
#                              merra_bef[qu], merra_aft[qu])
#     return merra


# %% Used in graph.py
def scatter(stations, s_var, merra, m_var, ax):
    """
    Function creating scatterplots between MERRA2 and AERONET
    Parameters
    ----------
    stations : DataFrame
        All AERONET data
    merra : DataFrame
        The MERRA2 3hourly values closest to the AERONET measurement time
    var : string
        AERONET and MERRA2 column name, has to be the same between the two
    ax : AxesSubplot object
        In which subplot the scatterplot should be placed
    Returns
    -------
    None.
    """
    # Capped number of points for development efficiency, REMOVE later
    x = stations[s_var]  # .iloc[:10000]
    y = merra[m_var]  # pd.Series(y, name=Mvar)
    # y = y  # .iloc[:10000]

    rmse = mean_squared_error(x, y, squared=False)  # Root Mean Square Error
    bias = y.mean()-x.mean()  # Systematic difference, i.e. bias
    rbias = bias/x.mean()*100  # Relative bias
    result = linregress(x, y)  # Object keeping the linear rergression results

    # The gaussian_kde take a looong time if many points are involved
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    sns.scatterplot(x=x, y=y, c=z, cmap="viridis", s=5, linewidths=0)

    # Alternatively, use this without kde
    # For some reason seaborne does not show all points ?! I switched to pyplot
    # plt.plot(x, y, '.')
    if "AOD" in s_var:  # If ploting AODs begin axes from 0
        plt.plot([0, x.max()], [0, x.max()], 'b:')
    else:  # begin axes from the minimum values
        plt.plot([x.min(), x.max()], [x.min(), x.max()], 'b:')
    ax.set_aspect('equal')
    plt.plot([x.min(), x.max()], result.intercept +
             result.slope*np.array([x.min(), x.max()]), 'r:')

    plt.text(0.05, 0.95, f'CC={result.rvalue:.3f}', horizontalalignment='left',
             size='medium', color='blue', transform=ax.transAxes)
    plt.text(0.05, 0.90, f'Slope={result.slope:.3f}$\pm${result.stderr:.3f}',
             horizontalalignment='left', size='medium', color='blue',
             transform=ax.transAxes)
    plt.text(0.05, 0.85, f'BIAS={bias:.3f}', horizontalalignment='left',
             size='medium', color='blue', transform=ax.transAxes)
    plt.text(0.05, 0.80, f'RBIAS={rbias:.1f}%', horizontalalignment='left',
             size='medium', color='blue', transform=ax.transAxes)
    plt.text(0.05, 0.75, f'RMSE={rmse:.3f}', horizontalalignment='left',
             size='medium', color='blue', transform=ax.transAxes)
    plt.xlabel("AERONET " + m_var)
    plt.ylabel("MERRA-2 " + m_var)
    return None


def cappa(sae, aae):
    """
    Function creating the scatterplot and performing the classification
    of Cappa et al. (2016)
    Parameters
    ----------
    sae : Series
        Scattering Angstrom Exponent
    aae : Series
        Absorption Angstrom Exponent
    Returns
    -------
    aertyp : Series
        Aerosol type by Cappa et al. (2016)
    """
    aertyp = np.empty(sae.shape[0])
    aertyp[:] = np.nan
    aertyp = pd.Series(aertyp)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')

    aertyp = np.where((sae <= 0.2) & (aae >= 2.0), "Dust", aertyp)
    aertyp = np.where(((sae <= 1.5) & (aae >= 1.5) & (aae < 2.0)) |
                      ((aae >= 2.0) & (sae > 0.2) & (sae <= 1.5)),
                      "Dust/BC/BrC", aertyp)
    aertyp = np.where((sae > 1.5) & (aae >= 2.0), "BrC", aertyp)
    aertyp = np.where((sae > 1.5) & (aae >= 1.5) & (aae < 2),
                      "BC/BrC", aertyp)
    aertyp = np.where(((sae <= 1.) & (aae >= 1.) & (aae < 1.5)) |
                      ((sae > 1) & (sae <= 1.5) & (aae >= sae) & (aae < 1.5)),
                      "Large BC mix", aertyp)
    aertyp = np.where(((sae > 1.5) & (aae >= 1.) & (aae < 1.5)) |
                      ((sae > 1.) & (sae <= 1.5) & (aae < sae) & (aae >= 1.)),
                      "BC dominated", aertyp)
    aertyp = np.where((sae <= 1.0) & (aae < 1.0),
                      "Large low abs./ Large black",
                      aertyp)
    aertyp = np.where((sae > 1.0) & (aae < 1.0),
                      "Small/Low abs. mix", aertyp)

    # The gaussian_kde take a looong time if many points are involved
    xy = np.vstack([sae[:], aae[:]])
    z = gaussian_kde(xy)(xy)
    sns.scatterplot(x=sae[:], y=aae[:], c=z, cmap="viridis", s=5,
                    linewidths=0)
    colo = (0.5, 0, 0)
    plt.text(1.6, -1, "Small/Low abs. mix", fontsize='small', color=colo)
    # plt.text(1.6, 1.17, "BC dominated", color=colo)
    plt.text(2.3, 1.17, "BC dominated", fontsize='small', color=colo)
    # plt.text(1.8, 1.67, "BC/BrC", color=colo)
    plt.text(2.2, 1.67, "BC/BrC", fontsize='small', color=colo)
    # plt.text(-0.3, 1.67, "Dust/BC/BrC", color=colo)
    plt.text(-1, 1.6, "Dust/\nBC/BrC", fontsize='small', color=colo)
    plt.text(-1, -1.1, "Large low abs./\nLarge black", fontsize='small',
             color=colo)
    plt.text(-1, 3, "Dust", fontsize='small', color=colo)
    plt.text(-1, 1.1, "Large BC\nmix", fontsize='small', color=colo)
    plt.text(1.6, 2.5, "BrC", fontsize='small', color=colo)
    # sns.scatterplot(x=sae, y=aae, s=5, hue=aertyp, legend=False)
    plt.vlines(x=[0.2, 1, 1.5], ymin=[2, aae.min(), 1.5],
               ymax=[aae.max(), 1, aae.max()],
               linestyle='dashed', color=colo)
    plt.hlines(y=[2, 2, 1.5, 1], xmin=[-1, 1.5, -1, -1],
               xmax=[0.2, 3.5, 3.5, 3.5], linestyle='dashed', color=colo)
    plt.plot([1, 1.5], [1, 1.5], linestyle='dashed', color=colo)

    return aertyp


def taylor_by_attribute(stations, s_var, merra, m_var, attribute_name):
    attrs = np.unique(stations[attribute_name])
    sdev = np.empty(len(attrs)+1)
    crmsd = np.empty(len(attrs)+1)
    ccoef = np.empty(len(attrs)+1)
    sdev[0] = 1.
    crmsd[0] = 0.
    ccoef[0] = 1.

    for attr in attrs:
        t_s = sm.taylor_statistics(merra.loc[stations[attribute_name] == attr,
                                             m_var],
                                   stations.loc[stations[attribute_name] ==
                                                attr, s_var])
        i = np.where(attrs == attr)[0] + 1
        sdev[i] = t_s['sdev'][1]/t_s['sdev'][0]
        crmsd[i] = t_s['crmsd'][1]/t_s['sdev'][0]
        ccoef[i] = t_s['ccoef'][1]

    attrs = np.insert(attrs, 0, 'AERONET')
    attrs = np.ndarray.tolist(attrs)

    if "SSA" in s_var:
        (title_STD, Legend, markerlabel) = ("off", "on", attrs)
    elif "AAOD" in s_var:
        (title_STD, Legend, markerlabel) = ("off", "on", [])
    else:
        (title_STD, Legend, markerlabel) = ("on", "on", [])
    sm.taylor_diagram(sdev, crmsd, ccoef, checkStats="on",
                      labelRMS="",  # titleRMSDangle=140,
                      titleSTD=title_STD,
                      titleOBS='AERONET', markerOBS="*", styleOBS='-',
                      markerSize=12, alpha=0.0, markerLegend=Legend,
                      markerLabel=markerlabel, markerColor="r")
    # plt.title(var + " nm")

    return None


def taylorWL(stations, merra, varis):
    """
    Produce the Taylor diagram
      Note that the first index corresponds to the reference series for
        the diagram. For example sdev[0] is the standard deviation of the
        reference series and sdev[1:4] are the standard deviations of the
        other 3 series. The value of sdev[0] is used to define the origin
        of the RMSD contours. The other values are used to plot the points
        (total of 3) that appear in the diagram.
    """
    sdev = np.empty(len(varis)+1)
    crmsd = np.empty(len(varis)+1)
    ccoef = np.empty(len(varis)+1)
    sdev[0] = 1.
    crmsd[0] = 0.
    ccoef[0] = 1.
    for var in varis:
        if "AOD" in var:
            taylor_stats = sm.taylor_statistics(merra[var],
                                                stations[var + "550"])
        elif "SSA" in var:
            taylor_stats = sm.taylor_statistics(merra[var],
                                                stations[var + "550"])
        i = varis.index(var)+1
        sdev[i] = taylor_stats['sdev'][1]/taylor_stats['sdev'][0]
        crmsd[i] = taylor_stats['crmsd'][1]/taylor_stats['sdev'][0]
        ccoef[i] = taylor_stats['ccoef'][1]

    varis.insert(0, "AERONET")

    # sm.taylor_diagram(sdev, crmsd, ccoef, checkStats="on",
    #                   labelRMS="NorUnRMS", titleRMSDangle=140,
    #                   titleOBS='AERONET', markerOBS="*", styleOBS='-',
    #                   markerSize=12, alpha=0.0, markerLegend='on',
    #                   markerLabel=varis, markerColor="r")

    sm.taylor_diagram(sdev, crmsd, ccoef, checkStats="on",
                      labelRMS="NorUnRMS", titleRMSDangle=140,
                      titleOBS='AERONET', markerOBS="*", styleOBS='-',
                      markersize=12, markerLabel=varis, markerLabelColor='r')
    return None
