#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""

(i) Plot functions for monda HSP data download show irradiance spectra 
(ii) Method to extract AOT (Aerosol Optical Thicknes susing direct-beam method 
(based on wood et al. 2017).
(iii)  Plot functions for AOT spectra.

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - July 2022
John Wood -  john@peakdesign.co.uk - Jully 2022

"""


import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

from scipy.interpolate import interp1d # interpolation and filtering

import ephem # library for solar angle computations


log = logging.getLogger('sorad-plotter')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def wl_find(wl, target):
    "Finds wavelength index nearest to target"
    return np.argmin((np.array(wl)-target)**2)


def find_max_in_wl_range(wl, spectra, minwl, maxwl):
    "find peak amplitude wihtin a wavelength range for a set of spectra"
    wlmin = wl_find(wl, 400)
    wlmax = wl_find(wl, 800)
    maxima = np.amax(spectra[:, wlmin:wlmax], axis=0)
    return maxima



def plot_irradiance_hsp(wl, eds, ed, IDR, time, file_id, target):
    "plot function for eds, ed and IDR"
   
    eds = 1000*eds # convert to mW
    ed = 1000*ed
   
    plt.figure(figsize=(18,14))
    plt.rc('font', size=24)
    plt.suptitle(file_id)
    
    colors_1 = cm.autumn(np.linspace(0, 1, len(time))) 
    colors_2 = cm.cool(np.linspace(0, 1, len(time)))
    
    # 1 Spectral plot for Ed
    plt.subplot(2,2,1)
    plt.title('Total downwelling irradiance spectra')
    plt.ylabel('$E_{d}$ [mW m$^{-2}$ nm$^{-1}$]')  
    for i in range(len(ed)):
       plt.plot(wl, ed[i,:], linewidth=1.2, color=colors_1[i], alpha=0.5)
    plt.xlim(350,900)
    plt.xticks([400, 500, 600,700, 800, 900])
    plt.legend(loc=1)
    plt.grid()
    plt.xlabel('Wavelength [nm]')
    maxwl = max(find_max_in_wl_range(wl, ed, 400, 800))
    plt.gca().set_ylim(bottom = 0, top = maxwl)
                       
    # 2. Spectral plot for Eds
    plt.subplot(2,2,2)
    plt.title('Diffuse irradiance spectra')
    plt.ylabel('$E_{ds}$ [mW m$^{-2}$ nm$^{-1}$]')
    for i in range(len(eds)):
        plt.plot(wl,eds[i,:], linewidth=1.2, color = colors_2[i], alpha=0.5)
    plt.xlim(350,900)
    plt.xticks([400, 500, 600,700, 800, 900])
    plt.grid()
    plt.xlabel('Wavelength [nm]')
    maxwl = max(find_max_in_wl_range(wl, eds, 400, 800))
    plt.gca().set_ylim(bottom = 0, top = maxwl)
    
    # 3. Time series of spectral maxima
    plt.subplot(2,2,3)
    plt.title('Time series for spectral maxima')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.ylabel('Irradiance [mW m$^{-2}$ nm$^{-1}$]')
    for i in range(len(ed)):
        plt.plot_date(time[i],(np.nanmax(ed,axis=1))[i],color=colors_1[i],ms=5)
        plt.plot_date(time[i],(np.nanmax(eds,axis=1))[i],color=colors_2[i],ms=5)
    middle_index =  int(np.round(len(eds)/2))
    plt.plot_date(time[middle_index], (np.nanmax(ed,axis=1))[middle_index], color = colors_1[middle_index], ms = 5, label = '$E_{d}$')
    plt.plot_date(time[middle_index], (np.nanmax(eds,axis=1))[middle_index], color = colors_2[middle_index], ms = 5,label =' $E_{ds}$')
    plt.gca().set_ylim(bottom =0)
    plt.legend()
    plt.grid()
    plt.xlabel('UTC time [hrs]')
    
    # 4. Time series of IDR
    plt.subplot(2,2,4)
    plt.title('Time series for integrated diffuse ratio')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.ylabel('IDR [no units]')
    plt.plot_date(time,IDR,color='black',ms=5)
    plt.gca().set_ylim(bottom = 0 , top = 1)
    plt.grid()
    plt.xlabel('UTC time [hrs]')
    
    plt.tight_layout()
    
    plt.savefig(os.path.join(target, file_id + 'hsp_irradiance.png'), format='png', dpi=150)

    return


# sub routines to calculate AOT & spectral plots
def calc_solar_zenith(time, lat, lon):
    'solar zenith angle function from ephem library. Acts on timstamp, and lat-lon vectors'''
    solar_zenith = np.nan*np.ones(len(time))
    for i in range(len(time)):    
        obs = ephem.Observer()
        sun = ephem.Sun()
        obs.date = time[i]
        obs.lat, obs.lon = str(lat[i]), str(lon[i])
        sun.compute(obs)
        solar_zenith[i] = 90 - (sun.alt * 180. / np.pi)
   
    return solar_zenith

def calc_rayleigh(wl, pressure = 1013.25):
    'function for Rayleigh component of atmospheric optical thickness'
    a1 = 117.25942 
    a2 = -1.3215
    a3 = 0.00032073
    a4 = -0.000076842
    p0 = 1013.25 
    wl  = wl*1e-3 # convert wl from nano to micrometres (for Rayleigh formula)    
    
    tau_r = (pressure/p0)*(1/(a1*wl**4 + a2*wl**2 + a3 + a4/wl**2))
    
    return tau_r

def calc_air_mass(solar_zenith, pressure = 1013.25):
    'function for atmospheric air mass'
    a = 0.50572 
    b = 96.07995
    c = 1.6364
    C = np.cos(np.pi/180.0 * solar_zenith)
    mu = C + a*(b - solar_zenith)**(-c) # inverse of atm. air mass 
    
    m = pressure / 1013.25 / mu # used as defintion of air mass
    
    return m


def calc_aot_direct(ed, eds, e_solar, time, solar_zenith, wl):
    'function for total and Aerosol component of atmospheric optical thickness using'
    'direct beam method (Wood et al. 2017)'
   
    edd = ed - eds # direct component
    
    # initalize tau_r, tau_t, tau_a
    tau_r = calc_rayleigh(wl)   
    tau_t = np.nan*np.ones([len(edd),len(edd.T)]) # total OT 
    tau_a = np.nan*np.ones([len(edd),len(edd.T)]) # aerosol OT (AOT)
    tau_err = np.nan*np.ones([len(edd)]) # error on tau cpts - assumed same for tau_t and tau_a (no error on tau_r) and independent of wavelength
    
    # empirical correction from Wood et al. 2017
    wl_e = np.array([440, 500, 675, 870, 1020]) #  wl bins for empirical correction
    am_e = np.array([1, 2, 3, 6, 10])   # airmasses for empirical corrections (Wood 2017)
    offset_am = np.array([0.0097, 0.0177, -0.0033, -0.0067, -0.0117]) # regression coefficients
    offset_wl = np.array([0.0244, 0.0260, 0.0182, 0.0124, 0.0457])
    slope_wl = np.array([1.2701, 1.2893, 1.3549, 1.4522, 1.5237])
    
    offset_am = interp1d(am_e, offset_am, kind = 'linear', axis = 0, fill_value = "extrapolate") # piecewise linear interpolation
    offset_wl = interp1d(wl_e, offset_wl, kind = 'linear', axis = 0, fill_value = "extrapolate")(wl) 
    slope_wl = interp1d(wl_e, slope_wl, kind = 'linear', axis = 0, fill_value = "extrapolate")(wl) 
    
    # compute tau cpts using Beer's law
    epsilon = 0.01673 # eccentricity
    for i in range(len(edd)):
        day_of_year = time[i].timetuple().tm_yday
        r_au = 1 - epsilon*np.cos((2*np.pi/365.25)*(day_of_year  - 4)) 
        e_toa = e_solar/r_au**2 # top of atmosphere irradiance
        m = calc_air_mass(solar_zenith[i]) # solar air mass
        edd[i,:] = edd[i,:]/np.cos((np.pi/180)*solar_zenith[i]) # scaling sometimes referred to as `edd_ni' 
        tau_t[i,:] = -(1/m)*np.log(edd[i,:]/e_toa) # total optical thickness from beers' law
        tau_a[i,:] = tau_t[i,:] - tau_r  # tau_a that has not been emprically corrected
        tau_a[i,:] = (tau_a[i,:] - offset_am(m) - offset_wl) * slope_wl # tau_a that has been empirically corrected 
        tau_err[i] =  (1/m)*0.07  
        
    return tau_t, tau_a, tau_err


def plot_aot(wl, tau_a, IDR, time, file_id, target, IDR_threshold=0.5):
    "plot function for aot. IDR < 0.6 used as intial QC"
   
    plt.figure(figsize=(18,14))
    plt.rc('font', size=24)
    plt.suptitle(file_id)
    
    time = np.array(time)[IDR < IDR_threshold] # intial QC Step based on IDR
    tau_a = tau_a[IDR < IDR_threshold, :]
    
    colors_1 = cm.viridis(np.linspace(0, 1, len(time)))

    plt.subplot(2,2,1)
    plt.title('AOT spectra')
    for i in range(len(time)):
        if IDR[i] < IDR_threshold:
            plt.plot(wl, tau_a[i,:], linewidth=1.2, color=colors_1[i],alpha=0.5)
    plt.xlim(385,900)
    plt.xticks([400, 500, 600,700, 800, 900])
    plt.legend(loc=1)
    plt.ylim(0,1.5)
    plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('AOT spectra')
    plt.gca().set_ylim(bottom =0)
    
    plt.subplot(2,2,2)
    plt.title('Time series for AOT(500)')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.ylabel('AOT(500)')
    wl_500 =wl_find(wl,500)
    for i in range(len(time)):
        if IDR[i] < IDR_threshold:
            plt.plot_date(time[i], tau_a[i,wl_500], color=colors_1[i], ms=5) # assumes wl vector
    plt.ylim(0,1.5)
    plt.gca().set_ylim(bottom = 0)
    plt.grid()
    plt.xlabel('UTC time [hrs]')
    
    plt.tight_layout()
    plt.savefig(os.path.join(target, file_id + '_aot.png'), format='png', dpi=150)

    return