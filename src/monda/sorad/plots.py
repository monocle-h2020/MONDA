#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""

Plot functions for monda So-Rad data download illustrating: (ir)radiance spectra, quality
control procedures applied to the processing chain, Rrs results, atmospheric conditions
(Ls/Ed(400) ratio) and coverage map.

------------------------------------------------------------------------------
Tom Jordan - tjor@pml.ac.uk - Feb 2022.
"""


import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

log = logging.getLogger('sorad-plotter')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


# sub routines
def wl_find(wl, target):
    "Finds wavelength index nearest to target"
    return np.argmin((np.array(wl)-target)**2)

def percentile_amplitude(values):
    "determine the 99th percentile of a range to limit y axis range limit"
    values =values[np.isfinite(values)]
    values = values[values>=0.0]
    return np.percentile(values, 99)

def outlier_limit(values):
    "determine outlier amplitude to limit y axis range limit"
    values = values[np.isfinite(values)]
    values = values[values>=0.0]
    outlier_limit = np.median(values) + 1.5 * np.std(values)
    ylimit = np.max(values[values<outlier_limit])
    return ylimit

def find_max_in_wl_range(wl, spectra, minwl, maxwl):
    "find peak amplitude wihtin a wavelength range for a set of spectra"
    wlmin = wl_find(wl, 400)
    wlmax = wl_find(wl, 800)
    maxima = np.amax(spectra[:, wlmin:wlmax], axis=0)
    return maxima

# plot functions
def plot_ed_ls_lt(ed, ls, lt, time, wl, file_id, target):
    """Spectral plots for ed ls lt"""

    plt.figure(figsize=(8,16))
    plt.rc('font', size=16)
    plt.suptitle(str(file_id))

    plt.subplot(3,1,1)
    plt.plot(wl,ed.T,linewidth=0.4,alpha=0.6)
    plt.xlim(350,900)
    plt.gca().set_ylim(bottom =-0.001)
    ymax = outlier_limit(find_max_in_wl_range(wl, ed, 400, 800))  # limit y axis based on maximum values < outlier range
    plt.gca().set_ylim(top = ymax)
    #plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$E_{d}$ [mW m$^{-2}$ nm$^{-1}$]')

    plt.subplot(3,1,2)
    plt.plot(wl,ls.T,linewidth=0.4,alpha=0.6)
    plt.xlim(350,900)
    plt.gca().set_ylim(bottom =-0.001)
    ymax = outlier_limit(find_max_in_wl_range(wl, ls, 400, 800))  # limit y axis based on maximum values < outlier range
    plt.gca().set_ylim(top = ymax)
    #plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$L_{s}$ [mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$]')

    plt.subplot(3,1,3)
    plt.plot(wl,lt.T,linewidth=0.4,alpha=0.6)
    plt.xlim(350,900)
    plt.gca().set_ylim(bottom =-0.001)
    ymax = outlier_limit(find_max_in_wl_range(wl, lt, 400, 800))  # limit y axis based on maximum values < outlier range
    plt.gca().set_ylim(top = ymax)
    #plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$L_{t}$ [mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$]')

    plt.tight_layout()
    plt.savefig(os.path.join(target, file_id + '_Ed-Ls-Lt-spectra.png'), format='png', dpi=150)

    return


def plot_rrs_qc_3c(rrs, time, wl, q_1, q_2, q_3, file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_3) > 0:
        plt.figure()
        plt.figure(figsize=(12,8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:,q_2==1]*1000))/1000;

        plt.subplot(2, 2, 1)
        plt.title(f"No QC (n = {int(len(q_1))})")
        plt.plot(wl, rrs.T, linewidth=0.4, alpha=0.6)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 2, 2)
        plt.title(f"Rad. QC (n = {int(np.sum(q_1))})")
        plt.plot(wl, rrs.T[:,q_1==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplot(2, 2, 3)
        plt.title(f"Rad. + 3C QC (n = {int(np.sum(q_2))})")
        plt.plot(wl,rrs.T[:,q_2==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$   [sr$^{-1}$]')

        plt.subplot(2, 2, 4)
        plt.title(f"Rad. + 3C + Rrs QC (n = {int(np.sum(q_3))})")
        plt.plot(wl,rrs.T[:,q_3==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplots_adjust(hspace=0.5)

        plt.savefig(os.path.join(target, file_id + '_Rrs_QC.png'), format='png', dpi=150)

    return


def plot_rrs_qc_fp(rrs, time, wl, q_1, q_2, file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_2) > 0:
        plt.figure()
        plt.figure(figsize=(12, 8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:,q_2==1]*1000))/1000;

        plt.subplot(2, 2, 1)
        plt.title(f"No QC (n = {int(len(q_1))})")
        plt.plot(wl, rrs.T, linewidth=0.4, alpha=0.6)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 2, 2)
        plt.title(f"Rad QC (n = {int(np.sum(q_1))})")
        plt.plot(wl, rrs.T[:,q_1==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplot(2, 2, 3)
        plt.title(f"Rad + Rrs QC (n = {int(np.sum(q_2))})")
        plt.plot(wl, rrs.T[:,q_2==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 950)
        plt.ylim(-0.002, ymax) # force axis limits
        #plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplots_adjust(hspace=0.5)

        plt.savefig(os.path.join(target, file_id + '_Rrs_QC.png'), format='png', dpi=150)

    return


def plot_coveragemap(lat, lon , q, file_id, target, map_resolution=11):
    """ coverage map showing quality control filtered data: color scheme matches
    `results' plot"""
    if np.sum(q) > 0:
        colors= cm.cool(np.linspace(0, 1, int(sum(q))))

        plt.figure(figsize=(15,10))
        extent = [np.floor(np.min(lon*100))/100, np.ceil(np.max(lon*100))/100, np.floor(np.min(lat*100))/100, np.ceil(np.max(lat*100))/100]
        request = cimgt.StamenTerrain()
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_image(request, map_resolution)
        ax.set_extent(extent, ccrs.PlateCarree())
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels_top = gl.right_labels = False
        gl.xformatter =  LONGITUDE_FORMATTER
        gl.yformatter =  LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12,  'rotation':45}
        gl.ylabel_style = {'size': 12,  'rotation': 0}

        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.tick_params(labelsize=10)
        plt.scatter(lon,lat,s=8,color='gray',transform=ccrs.PlateCarree(), label='Failed QC')
        plt.scatter(lon[q==1][0],lat[q==1][0],s=15,color=colors[0],transform=ccrs.PlateCarree(),label='Passed QC')
        for i in range(int(sum(q))):
            plt.scatter(lon[q==1][i],lat[q==1][i],s=15,color=colors[i],transform=ccrs.PlateCarree())

        plt.rc('font', size=14)
        plt.title(str(file_id))
        plt.legend()

        plt.savefig(os.path.join(target, file_id + '_coverage-map.png'), format='png', dpi=150)

    return


def plot_results(ed ,ls, wl_out, rrs, rrswl, time, q, file_id, target):
    """ Results plot showing: (i) sky measurement conditions using ls(400)/ed(400) ratio,
    (ii) qc-filtered rrs. Panel (i) is used to illustrate timestamps passing qc. Color scale matches
    between panels so time series and rrs can be visually referenced """

    if np.sum(q) > 0:
        lambda_400 = wl_find(wl_out ,400)  # atmopsheric conditions, clear (<.25) to overcast (>.8)
        ls_ed_400 = ls[:,lambda_400] / ed[:,lambda_400]

        colors = cm.cool(np.linspace(0, 1, int(sum(q))))  # colour map to match rrs with time series
        timestamp = np.array(time)  # convert to np array for masking

        plt.figure(figsize=(10, 14))
        plt.rc('font', size=19)

        plt.subplot(2, 1, 1)
        plt.suptitle(str(file_id))

        plt.title(f"Sky conditions: $\pi Ls(400)/Ed$(400) = {np.round(np.mean(np.pi*ls_ed_400),3)} +/- {np.round(np.std(np.pi*ls_ed_400), 3)}")
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
        plt.ylabel('Degree of cloudiness: \n $\pi L_{s}(400)/E_{d}$(400)')
        plt.plot_date(timestamp, np.pi*ls_ed_400, color='gray', ms=1, label = 'Failed QC')
        plt.plot_date(timestamp[q==1][0], np.pi*ls_ed_400[q==1][0], color=colors[0,:], ms=3 ,label = 'Passed QC')

        for i in range(int(sum(q))):
            plt.plot_date(timestamp[q==1][i], np.pi*ls_ed_400[q==1][i], color=colors[i,:], ms=3)

        plt.ylim(0, 1.6)
        plt.legend()
        #plt.grid()
        plt.xlabel('UTC time [hrs]')

        plt.subplot(2,1,2)
        plt.title('Remote-sensing reflectance: $R_{rs}$')

        for i in range(int(sum(q))):
            plt.plot(rrswl, rrs[q==1][i,:], color=colors[i,:], linewidth=0.6, alpha=0.6)

        plt.plot(rrswl, np.nanmean(rrs[q==1], axis=0), color='black', linewidth=2, alpha=1, label='Mean spectrum')
        plt.xlim(350, 900)
        plt.gca().set_ylim(bottom =-0.001)
        #plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$  [sr$^{-1}$]')
        plt.legend()

        plt.savefig(os.path.join(target, file_id + '_results.png'), format='png', dpi=150)

    return
