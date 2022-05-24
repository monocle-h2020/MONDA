#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""

Quality control functions for monda So-Rad data download

Quality control functions output a logical mask using convention:
 1 == valid timestamp, 0 == not valid timestamp.

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - Feb 2022
"""

import numpy as np
import datetime
import pandas as pd


# sub routines
def wl_find(wl, target):
    "Finds wavelength index nearest to target"
    return np.argmin((np.array(wl)-target)**2)


def combined_filter(q_1, q_2):
   """Combines 2 QC masks - tests for failures subject to logical OR condition"""

   q_combined = np.ones(len(q_1))
   for i in range(len(q_1)):
      if q_1[i] == 0 or q_2[i] == 0:
         q_combined[i] = 0

   return q_combined


def nearest(items, pivot):
     """Sub function used to locate nearest values indicies in list x relative to a pivot - used within rolling variability """

     nearest_value = min(items, key=lambda x: abs(x - pivot)) # finds nearest value
     nearest_index = items.index(min(items, key=lambda x: abs(x - pivot))) # finds nearest index

     return nearest_value, nearest_index


# Filters used in QC of l and e spectra (Step 1 in 3C and FP chain)
def qc_lt_ed_filter(ed, lt, time, wl, threshold = 0.020):
    """Funtion to filter by lt_ed ratio in NIR: basic implementation using absolute threshold defined in sr^-1 on [850, 950] nm"""

    lt_ed = np.divide(lt, ed)
    lambda_l = wl_find(800, wl)
    lambda_h = wl_find(950, wl)
    lw_ratio = np.nanmax(lt_ed[:,lambda_l:lambda_h], axis=1)

    q_ratio =  np.ones(len(ed))
    q_ratio[lw_ratio > threshold] = 0

    return q_ratio


def qc_ed_filter(ed, min_ed_threshold = 500):
    """Function to filter by the minimum (spectral max) ed value"""

    q_ed =  np.ones(len(ed)) # qc
    spec_max_ed = np.nanmax(ed.T, 0) # spectral max for Ed
    q_ed[spec_max_ed < min_ed_threshold] = 0

    return q_ed


def qc_ls_filter(ls, wl, threshold = 1):
    """Funtion to filter out ls spectra that are anomalously high in NIR"""

    lambda_1l = wl_find(550, wl) # central wavelength band
    lambda_1h = wl_find(580, wl)
    lambda_2l = wl_find(800, wl) # long wavelength band
    lambda_2h = wl_find(950, wl)

    ls_ratio = np.nanmax(ls[:,lambda_2l:lambda_2h], axis=1) / np.nanmax(ls[:,lambda_1l:lambda_1h], axis=1)

    q_ls =  np.ones(len(ls)) # qc mask: 1=valid, 0=not valid
    q_ls[ls_ratio > threshold] = 0

    return q_ls


# 3C-specific filtering (step 2 in 3C QC chain)
def qc_3cresidual(q_rad,resid,tol = 3):
    """3C residual fliter based on standard deviation of daily rmsd ditribution that has
       already passed radiometric quality control. tol is a standard deviation multiple
       with 3*sigma set to default"""

    q_rad_3c = q_rad
    if np.sum(q_rad) > 0:
        resid_max = np.ones(len(q_rad))*(tol*np.nanstd(resid[q_rad==1]) + np.nanmean(resid[q_rad==1]))
        q_rad_3c[resid > resid_max] = 0

    return q_rad_3c


def qc_3c_rho_filter(rho_ds, rho_dd, rho_s, upperbound = 0.1):
    """Function to remove 3C retreivals where rho factors terminate on upper optimzation bound
    (note: not strictly required for lower optimzation bound)"""

    q_rho3c =  np.ones(len(rho_ds))

    q_rho3c[np.abs(rho_ds  - upperbound) < 0.001] = 0
    q_rho3c[np.abs(rho_dd  - upperbound) < 0.001] = 0
    q_rho3c[np.abs(rho_s - upperbound) < 0.001] = 0

    return q_rho3c


# filters used in qc of Rrs spectra: Step 3 in 3C chain and 2 in FP chain)
def qc_SS_NIR_filter(wl, rrs, upperthreshold = 3, lowerthreshold = 0.5):
    """Filter based on NIR similarity spectrum (reflectance ratio of 779 and 865 nm) in Ruddick et al. 2006."""

    lambda779 = wl_find(wl, 779)
    lambda865 = wl_find(wl, 865)

    r779_865 = rrs[:,lambda779]/rrs[:,lambda865]

    q_ss =  np.ones(len(rrs)) # qc mask: 1=valid, 0=not valid
    q_ss[r779_865  > upperthreshold] = 0
    q_ss[r779_865  < lowerthreshold] = 0

    return q_ss


def qc_rrs_maxrange(rrs, upperthreshold = 0.1, lowerthreshold = 0.0):
    """ Filter based on range of rrs maximum """

    q_maxrange = np.ones(len(rrs))

    rrs_max = np.nanmax(rrs,axis=1)
    q_maxrange[rrs_max > upperthreshold] = 0
    q_maxrange[rrs_max < lowerthreshold] = 0

    return q_maxrange


def qc_rrs_min(rrs, wl):
    """ Replicates filtering step in FP processing that Rrs must be > 0 on [370,700] nm
    but applied to Rrs with offset included"""

    q_min = np.ones(len(rrs))

    lambda375 = wl_find(wl,375)
    lambda700= wl_find(wl,700)

    rrs_min = np.nanmin(rrs[:,lambda375:lambda700],axis=1)
    q_min[rrs_min <  0] = 0

    return q_min


# optional filters (not currently applied in FP and 3C chain)
def rolling_variability(spectrum, timestamp ,window_time, wl_range):
    """ Sub funtion to calculate rolling spectral variability metrics based on z-score (std normalized spectra) and percentage differences (mean normalized spectra).
    Based on function in Groetsch et. al 2017

    # Work flow:
    # (i) converts timestamp and spectra to data frame format (indexed by timestamp, rows wavelengths)
    # (ii) calculates rolling mean and standrad deviations using dataframe method
    # (iii) re-references rolling mean/std to window center
    # (iv) computes rolling z-score and percentage difference metrics
    # Input - window size in minutes and l or e spectrum
    # Notes - data gaps in rolling mean/std are handled via linear interpolation (default setting of rolling functions in pandas)
    """

    # (i) dataframe format
    window_time_secs = window_time*60 # converts to secs
    timeindex = [pd.Timestamp(timestamp[i]) for i in range(len(timestamp))] # converts to pandas format
    spectrum_DF = pd.DataFrame(spectrum,timeindex,columns = wl_range) # creates dataframe for spectum indexed by timestamp

    # (ii) Rolling means
    spectrum_Rmean = np.array(spectrum_DF.rolling(window=str(window_time)+'min').mean()) # rolling mean of spectrum
    spectrum_Rstd = np.array(spectrum_DF.rolling(window=str(window_time)+'min').std()) # rolling std of spectrum

    # (iii) Reference rolling mean and std to center of moving window (this is required as window center CANNOT be applied to ragged datetime index)
    timestamp_center = [timestamp[i]+datetime.timedelta(seconds=round(window_time_secs/2)) for i in range(len(timestamp))] # time-stamp shifted to window center
    spectrum_Rmean_centered = np.nan*np.ones([len(timestamp),len(wl_range)]) # centered versions of rolling mean and std
    spectrum_Rstd_centered = np.nan*np.ones([len(timestamp),len(wl_range)])
    for i in range(len(timestamp)):
        centertime = nearest(timestamp,timestamp_center[i]) # finds time and index of window center using nearest function
        centerindex = int(centertime[1]) # converts nearest index to int
        spectrum_Rmean_centered[i,:] = spectrum_Rmean[centerindex,:]
        spectrum_Rstd_centered[i,:] = spectrum_Rstd[centerindex,:]

    # (iv) Caculate moving z-score and percentage-differences
    pct_diff = np.nan*np.ones([len(timestamp),len(wl_range)]) # matrix for rolling percentage difference
    pct_diff_max = np.nan*np.ones(len(timestamp)) # spectral max of rolling pct difference
    pct_diff_mean = np.nan*np.ones(len(timestamp)) # spectral mean of rolling pct difference
    zscore = np.nan*np.ones([len(timestamp),len(wl_range)]) # matrix of (absolute) rolling z-score
    zscore_max =np.nan*np.ones(len(timestamp)) # spectral max of rolling z-score
    zscore_mean = np.nan*np.ones(len(timestamp)) # spectral mean of rolling z-score
    for i in range(len(timestamp)):
        pct_diff[i,:] = np.abs(100*(spectrum[i,:]-spectrum_Rmean_centered[i,:])/spectrum_Rmean_centered[i,:]) # rolling percentage difference
        pct_diff_max[i] = np.nanmax(pct_diff[i,:]) # spectral max
        pct_diff_mean[i] = np.nanmean(pct_diff[i,:]) # spectral mean
        zscore[i,:] = np.abs((spectrum[i,:]-spectrum_Rmean_centered[i,:])/spectrum_Rstd_centered[i,:]) # rollling absolute z-score
        zscore_max[i] = np.nanmax(zscore[i,:]) # spectral max
        zscore_mean[i] = np.nanmean(zscore[i,:]) # spectral mean

    return pct_diff_max, pct_diff_mean, zscore_max, zscore_mean


def qc_radiometric_variability(ed, lt, ls, time, wl, windowlength = 60, var_threshold = 1.1, var_metric = 'zscore_max'):
   """Filters out local variability anomalie of input l and e spectra based on rolling
      variability metric (spectral maxium of absolute z-score as default) """

   if var_metric ==  'pct_diff_max':
        index = 0
   elif var_metric ==  'pct_diff_mean':
        index = 1
   elif var_metric ==  'zscore_max':
       index =  2
   elif var_metric ==  'zscore_mean':
       index =  3

   ls_var = rolling_variability(ls, time, windowlength, wl)[index]  # caculates rolling z-score (spectral maximum)
   lt_var = rolling_variability(lt, time, windowlength, wl)[index]  # can also select pct_diff_max pct_diff_mean, zscore_max, zscore_mean using
   ed_var = rolling_variability(ed, time, windowlength, wl)[index]  # [0], [1], [3] indicies

   var_max = np.array([np.nanmax([ls_var[i], lt_var[i], ed_var[i]])  for i in range(len(ed))])  # takes `worst case variability metric' of 3 sensors

   q_var =  np.ones(len(ed))  # qc mask: 1=valid, 0=not valid
   q_var[var_max > var_threshold] = 0

   return q_var


def qc_coastalwater_rrsfilter(rrs, wl):
    """Rrs-shape based filter from Warren et al. 2019"""

    q_coastal =  np.ones(len(rrs)) # qc mask: 1=valid, 0=not valid

    lambda_350 = wl_find(350,wl)   # filter 1:  remove mean rrs on [350, 400] nm  < -0.0005 sr^-1
    lambda_400 = wl_find(400,wl)
    mean_rrs_350_400 = np.nanmean(rrs[:,lambda_350:lambda_400],axis=1)
    q_coastal[mean_rrs_350_400 < -0.0005] = 0

    lambda_800 = wl_find(800,wl)   # filter 2: remove mean rrs on [800, 900] nm < -0.0005 sr^-1
    lambda_900 = wl_find(900,wl)
    mean_rrs_800_900 = np.nanmean(rrs[:,lambda_800:lambda_900],axis=1)
    q_coastal[mean_rrs_800_900 < -0.0005] = 0

    max_rrs = np.max(rrs,axis=1)   # filter 3: remove max rrs > 0.015 sr^-1
    q_coastal[max_rrs > 0.015] = 0

    lambda_760 = wl_find(760,wl) # filter 4: oxygen band SNR filter - remove
    lambda_770 = wl_find(770,wl) # local peak intensities on [760, 770] nm that
    lambda_560 = wl_find(560,wl) # are > 10 %  of the green rrs peak
    lambda_600 = wl_find(600,wl)
    NIR_peak = np.max(rrs[:,lambda_760:lambda_770],axis=1) - np.min(rrs[:,lambda_760:lambda_770],axis=1)
    green_peak =  np.max(rrs[:,lambda_560:lambda_600],axis=1)
    q_coastal[NIR_peak/green_peak  > lambda_600] = 0

    max_rrs_index = np.argmax(rrs,axis=1) # filter 5: remove case when peak rrs occurs > 600 nm
    lambda_600 = wl_find(600,wl)
    q_coastal[max_rrs_index > lambda_600] = 0

    return q_coastal
