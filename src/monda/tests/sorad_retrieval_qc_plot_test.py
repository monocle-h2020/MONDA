# -*- coding: utf-8 -*-
"""

Monda example script to retrieve So-Rad Rrs and (ir)radiance spectra from PML Geoserver using WFS standard,
filter data using Quality Control principles, plot results, save Rrs spectra and some metadata (including QC masks).
The data and plots will be written into daily binned files.

---------------------------------------------------------------------------------

The Rrs spectra are available from two separate algorithms:
- Fingerprint (FP), Simis and Olsson 2013, RSE
- 3C, Groetsch et al. 2017, OE

--------------------------------------------------------------------------------

The FP and 3C outputs have different structures, quality control chains and plot functions: use the -a argument to specify which should be used.

For command line options please see
    python sorad_retrieval_qc_plot_test.py -h


For example,
    python sorad_retrieval_qc_plot_test.py -a 3c -p PML_SR004 -i "2021-08-16 00:00:00" -e "2021-08-22 23:59:59" -r -m -c -g --bbox 35 -10 40 -8
will plot and save all results from So-Rad platform PML_SR004 for the given time period, with Rrs derived from the 3C algorithm.
Files will be written to a default folder in the current directory (override this by providing an alternative path with the -t argument).
An optional bounding box is specified providing four space-separated coordinates (corner 1 lat lon corner 2 lat lon).

-------------------------------------------------------------------------------

Default quality control chain for 3C:
       (0)   Filters based on relative azimuth and tilt angle
       (i)   Radiometric filtering (Lt/Ed glint filter, Ed and Ls anomaly filters)
       (ii)  Algorithmic filtering (3c residual and termination of rho_dd, rho_ds, or rho_ds at optimizion bounds)
       (iii) Rrs-based filtering (similarity spectrum, and filter based on range of Rrs)

Default quality control chain for FP:
       (0)   Filters based on relative azimuth and tilt angle
       (i)   Radiometric filtering (Lt/Ed glint filter, Ed and Ls anomaly filters)
       (ii)  Rrs-based filtering (similarity spectrum, and filter based on max/min of Rrs)
       Note: `algorithmic filtering' is already applied to FP data on the geoserver (based on rho_s bounds)

Optional filters (commented out below):
        - Rrs shape based filtering (currently includes `coastal water filter' from Warren et al. 2019).
           Users are encouraged to apply filters for their local water type
        - Variability filtering (based on z-score), following Groetsch et al. 2017.

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - Feb 2022
                            - Spring 2024
"""

import sys
import os
import numpy as np

from monda.sorad import access, plots, qc 

import datetime
import logging
#import pandas as pd
import argparse
#from math import ceil

log = logging.getLogger('sorad-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def run_example(platform_id = 'PML_SR001',
                start_time = datetime.datetime(2023,10,10,0,0,0),
                end_time   = datetime.datetime(2023,10,10,23,59,59),
                bbox = None,
                target='.', rrsalgorithm='3c',
                output_radiance = True,
                output_metadata = True,
                output_rrs = True,
                output_plots = True):

    """
    Download data from a specific So-Rad platform processed to Rrs using the Fingerprint (fp) or 3C (3c) algorithm.
    Apply quality control filters and plot results, save output in daily bins between requested start and end times

    platform_id: the serial number of a So-Rad platform, or None
    start_time:  date/time to start collecting data from
    end_time:    date/time to start collecting data from
    bbox:        Bounding box corner coordinates (lat,lon,lat,lon)
    target:      destination folder for plots/data
    rrsalgorithm Fingerprint (fp) or 3c Rrs processing algorithm
    """

    wl_output = np.arange(300, 1001, 1)  # output range to interpolate (ir)radiance spectra. Note that Rrs is already interpolated to 1-nm over the available sensor range.

    if rrsalgorithm == 'fp':
        layer = 'rsg:sorad_public_view_fp_full'
    elif rrsalgorithm == '3c':
        layer = 'rsg:sorad_public_view_3c_full'
    else:
        log.error(f"{rrsalgorithm} is not a valid choice and must be one of [fp, 3c]")

    initial_day = start_time.date()
    final_day = end_time.date()
    assert final_day >= initial_day

    # split the request into discrete days
    days = (final_day - initial_day).days + 1
    log.info(f"Request spans {days} day(s)")

    for i in range(days):
        this_day      = initial_day + datetime.timedelta(days = i)
        datetime_i    = datetime.datetime(this_day.year, this_day.month, this_day.day, 0, 0, 0)
        datetime_e    = datetime.datetime(this_day.year, this_day.month, this_day.day, 23, 59, 59, 999999)

        if datetime_i < start_time:
            datetime_i = start_time
        if datetime_e > end_time:
            datetime_e = end_time

        log.info(f"Request timeframe {datetime_i.isoformat()} - {datetime_e.isoformat()}")

        response = access.get_wfs(platform = platform_id,
                                  timewindow = (datetime_i, datetime_e),
                                  layer=layer, bbox=bbox)

        log.info(f"{response['length']} features received.")
        if response['length'] == 0:
            continue
    
        file_id = f"{datetime_i.strftime('%Y-%m-%d')}_{rrsalgorithm.upper()}"   # labelling for output data files and plots

        rrswl, time, lat, lon, rel_view_az,\
               ed, ls, lt, rrs,\
               sample_uuids, platform_ids, platform_uuids,\
               gps_speeds, tilt_avgs, tilt_stds = access.unpack_response(response, rrsalgorithm, wl_output)

        if output_plots:
            log.info("Creating (ir)radiance plots")
            plots.plot_ed_ls_lt(ed, ls, lt, time, wl_output, file_id, target)

        # Step (0) QC filters based on relative aziumuth and tilt/tilt std
        q_az =  qc.rel_az_filter(rel_view_az, lower_azi_bound = 110, upper_azi_bound = 150)
        q_tilt, q_tilt_std =  qc.tilt_filter(tilt_avgs, tilt_stds, upper_tilt_bound=5, upper_tilt_std_bound = 2)
        q_0 =  qc.combined_filter(q_az , qc.combined_filter(q_tilt, q_tilt_std))
   
        # Step (i) radiometric quality control filters (i.e. QC applied to l or e spectra)
        q_lt_ed = qc.qc_lt_ed_filter(ed, lt, time, wl_output, threshold = 0.020) # lt/Ed ratio (glint) filtering
        q_ed =    qc.qc_ed_filter(ed, min_ed_threshold = 500) # filters on ed and ls anomalies
        q_ls =    qc.qc_ls_filter(ls, wl_output, threshold = 1)
        q_1 =     qc.combined_filter(qc.combined_filter(qc.combined_filter(q_lt_ed, q_ed), q_ls), q_0) # combined `radiometric' qc mask - inlcludes data passing step(0)

        # Step (ii): (3c) algorithmic qc filters specfic to 3C (rmsd or termination at rho bounds)
        if rrsalgorithm == '3c':
            rmsd_3c =     np.array([response['result'][i]['c3_rmsd'] for i in range(len(response['result']))]) # rmsd values used in 3C residual filter
            rho_ds =      np.array([response['result'][i]['c3_rho_ds'] for i in range(len(response['result']))]) # rho factors
            rho_dd =      np.array([response['result'][i]['c3_rho_dd'] for i in range(len(response['result']))])
            rho_s =       np.array([response['result'][i]['c3_rho_s'] for i in range(len(response['result']))])
            q_rho =       qc.qc_3c_rho_filter(rho_ds, rho_dd, rho_s, upperbound = 0.1) # removes data where rho terminates at optimization bounds
            q_1_resid = qc.qc_3cresidual(q_1, rmsd_3c, tol = 1.5) # removes data where residual parameter is above threshold standard-deivation multiple
            q_2 =    qc.combined_filter(q_rho, q_1_resid)
        elif rrsalgorithm == 'fp':
            q_2 = np.nan*np.ones(len(q_1)) # q2 stored as NaN for fp
           
        # Step (iii):  addtional qc metrics that apply to Rrs spectrum
        q_ss =        qc.qc_ss_nir_filter(rrswl, rrs, upperthreshold = 3, lowerthreshold = 0.5)  # similarity spectrum filter
        q_maxrange =  qc.qc_rrs_maxrange(rrs, upperthreshold = 0.1, lowerthreshold = 0.00)    # filters on max and min rrs
        q_min =       qc.qc_rrs_min(rrs, rrswl)
        
        # Optional filters
        # q_coastal = qc_coastalwater_rrsfilter(rrs, wl) #  filter based on expected shape of rrs - example from Warren 2019 used. users can input their own spectra here (will depend on water type)
        # q_var = qc_radiometric_variability(ed, lt, ls, time, wl, windowlength = 60, var_threshold =1.1, var_metric = 'zscore_max')

        if rrsalgorithm == '3c':
            q_3 = qc.combined_filter(q_2, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange)))) # recommended rrs qc mask for 3C method (combines step (i), (ii) and (iii) QC)
        elif rrsalgorithm == 'fp':
            q_3 = qc.combined_filter(q_1, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange)))) # recommended rrs qc mask for FP method (combines step (i) and (iii) QC)


        if output_plots and rrsalgorithm == 'fp':
            log.info("Creating Rrs plots")
            plots.plot_rrs_qc_fp(rrs, time, rrswl, q_1, q_3, file_id, target)
            plots.plot_coveragemap(lat, lon, q_3, file_id, target, map_resolution = 11)
            plots.plot_results(ed, ls, wl_output, rrs, rrswl, time, q_3, file_id, target)

        elif output_plots and rrsalgorithm == '3c':
            log.info("Creating Rrs plots")
            plots.plot_rrs_qc_3c(rrs, time, rrswl, q_0, q_1, q_2, q_3, file_id, target)
            plots.plot_coveragemap(lat, lon, q_3, file_id, target, map_resolution = 11)
            plots.plot_results(ed, ls, wl_output, rrs, rrswl, time, q_3, file_id, target)

        d = access.meta_dataframe(sample_uuids, platform_ids, time, lat, lon, gps_speeds, tilt_avgs, tilt_stds, rel_view_az, q_0, q_1, q_2, q_3)

        # Store outputs
        if output_metadata:
            d_filename = os.path.join(target, file_id + '_metadata.csv')
            if os.path.exists(d_filename):
                log.warning(f"File {d_filename} was overwritten")
            d.to_csv(d_filename)

        if output_rrs:
            r_filename = os.path.join(target, file_id + '_Rrs.csv')
            if os.path.exists(r_filename):
                log.warning(f"File {r_filename} was overwritten")
            np.savetxt(r_filename, rrs, delimiter=',',  header = ",".join([str(w) for w in rrswl]), fmt='%.8f')

        if output_radiance:
            header = ",".join([str(w) for w in wl_output])
            ls_filename = os.path.join(target, file_id + '_Ls.csv')
            lt_filename = os.path.join(target, file_id + '_Lt.csv')
            ed_filename = os.path.join(target, file_id + '_Ed.csv')
            if os.path.exists(ls_filename):
                log.warning(f"File {ls_filename} was overwritten")
            if os.path.exists(lt_filename):
                log.warning(f"File {lt_filename} was overwritten")
            if os.path.exists(ed_filename):
                log.warning(f"File {ed_filename} was overwritten")
            np.savetxt(ls_filename, ls, delimiter=',', header = header, fmt='%.8f')
            np.savetxt(lt_filename, lt, delimiter=',', header = header, fmt='%.8f')
            np.savetxt(ed_filename, ed, delimiter=',', header = header, fmt='%.8f')

    return response

def parse_args():
    """Interpret command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--algorithm',   required = False, type=str, default = 'fp', help = "rrs processing algorithm: fp or 3c")
    parser.add_argument('-p','--platform',    required = False, type = str, default = 'PML_SR004', help = "Platform serial number, e.g. PML_SR004.")
    parser.add_argument('-i','--start_time',  required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2021,10,21,0,0,0),
                                              help = "Initial UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-e','--end_time',    required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2021,10,22,23,59,59),
                                              help = "Final UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-b','--bbox',        required = False, type = float, nargs='+', default = None, help = "Restrict query to bounding box format [corner1lat corner1lon corner2lat corner2lon]")
    parser.add_argument('-t','--target',      required = False, type = str, default=None,
                                              help = "Path to target folder for plots (defaults to 'So-Rad_testoutput' in the current folder).")
    parser.add_argument('-r','--output_radiance',  required = False, action='store_true', help = "Output Ls, Lt, Ed spectra to csv file")
    parser.add_argument('-m','--output_metadata',  required = False, action='store_true', help = "Output metadata to csv file")
    parser.add_argument('-c','--output_rrs',  required = False, action='store_true', help = "Output Rrs spectra to csv file")
    parser.add_argument('-g','--output_plots',  required = False, action='store_true', help = "Output plots")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()

    if not any([args.output_radiance, args.output_metadata, args.output_rrs, args.output_plots]):
        log.warning("No plots or data outputs specified (see -h for help)")

    if args.target is None:
        args.target = os.path.join('.', 'So-Rad_test-output')

    if not os.path.isdir(args.target):
        os.mkdir(args.target)

    # response = run_example(args.platform, args.start_time, args.end_time, args.bbox, args.target, args.algorithm.lower(),  # use if running in terminal
    #                       args.output_radiance, args.output_metadata, args.output_rrs, args.output_plots)

    response = run_example() # use if