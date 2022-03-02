# -*- coding: utf-8 -*-
import sys
import argparse
print(sys.path)
"""

Monda test script to retrieve So-Rad Rrs and (ir)radiance spectra from PML Geoserver using WFS standard,
filter data using Quality Control principles, plot results, save Rrs spectra and core metadata (including QC masks)

Data access is provided from sorad.data_access.sorad_access, quality control is provided from sorad.qc.sorad_qc
and plots is provided from sorad.data_analysis.sorad_plots

---------------------------------------------------------------------------------

The Rrs spectra are available from two separate algorithms:
- Fingerprint (FP), Simis and Olsson 2013, RSE
- 3C, Groetsch et al. 2017, OE

--------------------------------------------------------------------------------

The FP and 3C outputs have different structures, quality control chains and plot functions: use the -a argument to specify which should be used.

Example usage via command line (see parse-args function in sorad.data_access.sorad_access)
    test_wfs_fpand3c.py  -p PML_SR004 -a 3c -i 2021-08-16 -f 2021-08-22 -c True -t test_folder

This will produce results from So-Rad platform with serial number PML_SR004 using 3C Rrs processing scheme
between 2021-08-16 and 2021-08-22 saved to test folder.

-------------------------------------------------------------------------------

Default quality control chain for 3C:
       (i)   Radiometric filtering (lt/ed glint filter, ed and ls anomaly filters)
       (ii)  Algorithmic filtering (3c residual and termination of rho_dd, rho_ds, or rho_ds at optimizion bounds)
       (iii) Rrs-based filtering (similarity spectrum, and filter based on range of Rrs)

Default quality control chain for FP:
       (i)   Radiometric filtering (lt/ed glint filter, ed and ls anomaly filters)
       (ii)  Rrs-based filtering (similarity spectrum, and filter based on max/min of Rrs)
       Note: `algorithmic filtering' is already applied to FP data on the geoserver (based on rho_s bounds)

Optional filters:
       (i) Rrs shape based filtering (currently includes `coastal water filter' from Warren et al. 2019). Users are encouraged to 
       apply filters for their local water type
       (ii) Variability filtering (based on z-score), following Groetsch et al. 2017.

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - Feb 2022
"""

import os
import numpy as np
import monda
from monda.sorad.qc import sorad_qc as qc
from monda.sorad.data_access import sorad_access as access
from monda.sorad.data_analysis import sorad_plots as plots
import datetime
import logging
import pandas as pd

log = logging.getLogger('sorad-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def run_example_fp(platform_id = 'PML_SR004',
                   final_day = datetime.date.today() - datetime.timedelta(days=1),
                   initial_day = datetime.date.today() - datetime.timedelta(days=1),
                   target='.', spectra = False):
    """ fp quality control, results plots, and data output (daily data binning) """

    days = (final_day - initial_day).days + 1
    for i in range(days):

        datetime_i = initial_day + datetime.timedelta(days = i)
        datetime_iplus1 = initial_day + datetime.timedelta(days = i+1)

        response = access.get_wfs(platform = platform_id, timewindow = (datetime_i, datetime_iplus1), layer='rsg:sorad_public_view_fp_full')
        file_id = platform_id + '_' + datetime_i.strftime('%Y-%m-%d') + '_FP'   # ID used for daily plot functions
        log.info(f"{response['length']} features received.")

        if response['length'] > 0:

            wl = np.arange(response['result'][0]['wl_grid'][0], response['result'][0]['wl_grid'][1]-1, response['result'][0]['wl_grid'][2])  # reconstruct wavelength grid for Rrs
            time = [response['result'][i]['time'] for i in range(len(response['result']))]
            lat = np.array([response['result'][i]['lat'] for i in range(len(response['result']))])
            lon = np.array([response['result'][i]['lon'] for i in range(len(response['result']))])
            rel_view_az = np.array([response['result'][i]['rel_view_az'] for i in range(len(response['result']))])

            ed = access.get_l1spectra(response, 'ed_', wl) # irradiance spectra in 2D matrix format: rows time index, columns wavelength
            ls = access.get_l1spectra(response, 'ls_', wl)
            lt = access.get_l1spectra(response, 'lt_', wl)

            rrs_fp = np.array([response['result'][i]['rrs'][:] for i in range(len(response['result']))]) # rrs spectra 2D matrix format: rows time index, columns wavelength
            offset_fp = np.array([response['result'][i]['nir_offset'] for i in range(len(response['result']))])
            rrs_fp  =  np.array([rrs_fp[i,:] - np.ones(len(wl))*offset_fp[i] for i in range(len(rrs_fp))]) # spectral offset (applied as default definition of FP rrs)

            plots.plot_ed_ls_lt(ed, ls, lt, time, wl, file_id, target)

            # Step (i) radiometric quality control filters (i.e. QC applied to l or e spectra)
            q_lt_ed = qc.qc_lt_ed_filter(ed, lt, time, wl, threshold = 0.020) # lt/Ed ratio (glint) filtering
            q_ed = qc.qc_ed_filter(ed, min_ed_threshold = 500) # filters on ed and ls anomalies
            q_ls = qc.qc_ls_filter(ls, wl, threshold = 1)
            q_rad = qc.combined_filter(qc.combined_filter(q_lt_ed, q_ed), q_ls) # combined `radiometric' qc mask

            # Step (ii)  addtional qc metrics that apply to rrs spectrum
            q_ss = qc.qc_SS_NIR_filter(wl, rrs_fp, upperthreshold = 3, lowerthreshold = 0.5) # similarity spectrum
            q_maxrange =  qc.qc_rrs_maxrange(rrs_fp, upperthreshold = 0.1, lowerthreshold = 0.00) # filters on max and min rrs
            q_min =  qc.qc_rrs_min(rrs_fp,wl)
            q_rad_rrs = qc.combined_filter(q_rad, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange)))) # recommended rrs qc mask for FP method (combines step (i) and (ii) QC)

            # Optional filters
            # q_coastal = qc_coastalwater_rrsfilter(rrs_fp, wl) #  filter based on expected shape of rrs - example from Warren 2019 used. users can input their own spectra here (will depend on water type)
            # q_var = qc_radiometric_variability(ed, lt, ls, time, wl, windowlength = 60, var_threshold =1.1, var_metric = 'zscore_max')

            plots.plot_rrs_qc_fp(rrs_fp, time, wl, q_rad, q_rad_rrs, file_id, target)
            plots.plot_results(ed, ls, rrs_fp, time, wl, q_rad_rrs, file_id, target)
            plots.plot_coveragemap(lat, lon, q_rad, file_id, target, map_resolution = 11)

            d = pd.DataFrame()   # store core metadata and qc flags in a data frame
            d['timestamp'] = time
            d['lat'] = lat
            d['lon'] = lon
            d['rel_view_az '] = rel_view_az
            d['q_1'] = q_rad     # Mask after step (i) QC
            d['q_2'] = q_rad_rrs # Mask after step (ii) QC. Currently recommended for FP rrs data analysis

            # optional to output all qc masks
            # q_keys = [i for i in locals() if i.startswith('q_')]
            # for i in range(len(q_keys)): # add all qc fields to the dataframe
            # d[str(q_keys[i])] = pd.Series(eval(q_keys[i]))

            if target != '.': # save daily output
                target_folder = os.path.join(os.getcwd() + '/' + target + '/')
                d_filename = target_folder + file_id + '_data.csv'
                d.to_csv(d_filename)
                r_filename = target_folder + file_id + '_rrs.csv'
                np.savetxt(r_filename, rrs_fp, delimiter=',',  header = str(list(wl))[1:-1], fmt='%.8f')
                if spectra == True:
                    ls_filename = target_folder + file_id + '_ls.csv'
                    np.savetxt(ls_filename, ls, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    lt_filename = target_folder + file_id + '_lt.csv'
                    np.savetxt(lt_filename, lt, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    ed_filename = target_folder + file_id + '_ed.csv'
                    np.savetxt(ed_filename, ed, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')

            elif target == '.':
                d_filename = file_id + '_data.csv'
                d.to_csv(d_filename)
                r_filename = file_id + '_rrs.csv'
                np.savetxt(r_filename, rrs_fp, delimiter=',',  header = str(list(wl))[1:-1], fmt='%.8f')
                if spectra == True:
                    ls_filename = file_id + '_ls.csv'
                    np.savetxt(ls_filename, ls, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    lt_filename = file_id + '_lt.csv'
                    np.savetxt(lt_filename, lt, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    ed_filename = file_id + '_ed.csv'
                    np.savetxt(ed_filename, ed, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')

    return response


def run_example_3c(platform_id = 'PML_SR004',
                   final_day = datetime.date.today() - datetime.timedelta(days=1),
                   initial_day = datetime.date.today() - datetime.timedelta(days=1),
                   target = '.', spectra = False):
    """ 3c quality control, results plots, and data output (daily data binning) """

    days = (final_day - initial_day).days + 1

    for i in range(days):
        datetime_i = initial_day + datetime.timedelta(days = i)
        datetime_iplus1 = initial_day + datetime.timedelta(days = i+1)

        response = access.get_wfs(platform=platform_id, timewindow = (datetime_i, datetime_iplus1), layer='rsg:sorad_public_view_3c_full')
        file_id = platform_id + '_' + datetime_i.strftime('%Y-%m-%d') + '_3C'
        log.info(f"{response['length']} features received.")

        if response['length'] > 0:

            wl = np.arange(response['result'][0]['c3_wl_grid'][0], response['result'][0]['c3_wl_grid'][1], response['result'][0]['c3_wl_grid'][2])  # core metadata
            time = [response['result'][i]['time'] for i in range(len(response['result']))]
            lat = np.array([response['result'][i]['lat'] for i in range(len(response['result']))])
            lon = np.array([response['result'][i]['lon'] for i in range(len(response['result']))])
            rel_view_az = np.array([response['result'][i]['rel_view_az'] for i in range(len(response['result']))])
            # gps_speed = [response['result'][i]['gps_speed'] for i in range(len(response['result']))]

            ed = access.get_l1spectra(response, 'ed_', wl) # # irradiance spectra in 2D matrix format: rows time index, columns wavelength
            ls = access.get_l1spectra(response, 'ls_', wl)
            lt = access.get_l1spectra(response, 'lt_', wl)
            rrs_3c = np.array([response['result'][i]['c3_rrs'][:] for i in range(len(response['result']))]) # 2D matrix format: rows time index, columns wavelength

            plots.plot_ed_ls_lt(ed, ls, lt, time, wl, file_id , target)

            # Step (i) radiometric quality control filters (i.e. QC applied to l or e spectra)
            q_lt_ed = qc.qc_lt_ed_filter(ed, lt, time, wl, threshold = 0.020) # lt/Ed ratio filtering
            q_ed = qc.qc_ed_filter(ed, min_ed_threshold = 500)  # ed anomalie
            q_ls = qc.qc_ls_filter(ls, wl, threshold = 1) #  ls anomalie
            q_rad = qc.combined_filter(qc.combined_filter(q_lt_ed, q_ed), q_ls) # combined `radiometric' qc mask

            # (ii) algorithmic qc filters specfic to 3C (rmsd or termination at rho bounds)
            rmsd_3c = np.array([response['result'][i]['c3_rmsd'] for i in range(len(response['result']))]) # rmsd values used in 3C residual filter
            rho_ds = np.array([response['result'][i]['c3_rho_ds'] for i in range(len(response['result']))]) # rho factors
            rho_dd = np.array([response['result'][i]['c3_rho_dd'] for i in range(len(response['result']))])
            rho_s = np.array([response['result'][i]['c3_rho_s'] for i in range(len(response['result']))])

            q_rho = qc.qc_3c_rho_filter(rho_ds,rho_dd,rho_s,upperbound = 0.1)
            q_rad_resid = qc.qc_3cresidual(q_rad, rmsd_3c, tol = 1.5)
            q_rad_3c = qc.combined_filter(q_rho, q_rad_resid)

            # Step (iii)  addtional qc metrics that apply to rrs spectrum
            q_ss = qc.qc_SS_NIR_filter(wl,rrs_3c, upperthreshold = 3, lowerthreshold = 0.5) # similarity spectrum
            q_maxrange =  qc.qc_rrs_maxrange(rrs_3c, upperthreshold = 0.1, lowerthreshold = 0.00) # max and min rrs filters
            q_min =  qc.qc_rrs_min(rrs_3c, wl)
            q_rad_3c_rrs = qc.combined_filter(q_rad_3c, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange)))) # recommended rrs qc mask for 3C method (combines step (i), (ii) and (iii) QC)

            #  (optional)
            # q_var = qc_radiometric_variability(ed, lt, ls, time, wl, windowlength = 60, var_threshold =1.1, var_metric = 'zscore_max')
            # q_coastal = qc_coastalwater_rrsfilter(rrs_3c,wl) #  filter based on expected shape of rrs - example from Warren 2019 used. users can input there own spectra here

            plots.plot_rrs_qc_3c(rrs_3c, time, wl, q_rad, q_rad_3c, q_rad_3c_rrs, file_id, target)
            plots.plot_results(ed, ls, rrs_3c, time, wl, q_rad_3c_rrs, file_id, target)
            plots.plot_coveragemap(lat, lon, q_rad_3c, file_id, target, map_resolution = 11)

            # metadata
            d = pd.DataFrame()
            d['timestamp'] = time
            d['lat'] = lat
            d['lon'] = lon
            d['rel_view_az'] = rel_view_az
            d['q_1'] = q_rad  # Mask after step (i) QC
            d['q_2'] =  q_rad_3c  # Mask after step (ii) QC
            d['q_3'] =  q_rad_3c_rrs # Mask after step (iii) QC: currently recommended for 3C  rrs data analysis
            # q_keys = [i for i in locals() if i.startswith('q_rad')]  #  code that can be used to output all qc fields
            # for i in range(len(q_keys)): # add all qc fields to the dataframe
            # d[str(q_keys[i])] = pd.Series(eval(q_keys[i]))

            if target != '.':
                target_folder = os.path.join(os.getcwd() + '/' + target + '/')
                d_filename = target_folder + file_id + '_data.csv'
                d.to_csv(d_filename)
                r_filename = target_folder + file_id + '_rrs.csv'
                np.savetxt(r_filename, rrs_3c, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                if spectra == True:
                    ls_filename = target_folder + file_id + '_ls.csv'
                    np.savetxt(ls_filename, ls, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    lt_filename = target_folder + file_id + '_lt.csv'
                    np.savetxt(lt_filename, lt, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    ed_filename = target_folder + file_id + '_ed.csv'
                    np.savetxt(ed_filename, ed, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')

            elif target == '.':
                d_filename = file_id + '_data.csv'
                d.to_csv(d_filename)
                r_filename = file_id + '_rrs.csv'
                np.savetxt(r_filename, rrs_3c, delimiter=',',  header = str(list(wl))[1:-1], fmt='%.8f')
                if spectra == True:
                    ls_filename = file_id + '_ls.csv'
                    np.savetxt(ls_filename, ls, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    lt_filename = file_id + '_lt.csv'
                    np.savetxt(lt_filename, lt, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')
                    ed_filename = file_id + '_ed.csv'
                    np.savetxt(ed_filename, ed, delimiter=',', header = str(list(wl))[1:-1], fmt='%.8f')

    return response


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--algorithm', required = False, type=str, default = 'fp', help = "rrs processing algorithm: fp or 3c")
    parser.add_argument('-p','--platform', required = False, type = str, default = 'PML_SR004', help = "Platform serial number, e.g. PML_SR004.")
    parser.add_argument('-i','--initial_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = datetime.date(2021,10,21) - datetime.timedelta(days=1),help = "Initial day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-f','--final_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = datetime.date(2021,10,22), help = "Final day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-t','--target', required = False, type = str, default='SoRad_testoutput', help = "Target folder for plots to be written to (defaults to current folder).")
    parser.add_argument('-s','--spectra', required = False, type = bool, default = False, help = "Output Ls, Lt, Ed spectra")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()

    if args.target != '.' :
        target_directory = os.path.join(os.getcwd() + '/' + args.target)
        if os.path.isdir(target_directory) == False:
            os.mkdir(os.path.join(os.getcwd() + '/' + args.target))

    if args.algorithm.lower() == 'fp':
        response = run_example_fp(args.platform, args.final_day, args.initial_day, args.target, args.spectra)
    elif args.algorithm.lower()  == '3c':
        response = run_example_3c(args.platform, args.final_day, args.initial_day, args.target, args.spectra)
