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
       (i)   Radiometric filtering (lt/ed glint filter, ed and ls anomaly filters)
       (ii)  Algorithmic filtering (3c residual and termination of rho_dd, rho_ds, or rho_ds at optimizion bounds)
       (iii) Rrs-based filtering (similarity spectrum, and filter based on range of Rrs)

Default quality control chain for FP:
       (i)   Radiometric filtering (lt/ed glint filter, ed and ls anomaly filters)
       (ii)  Rrs-based filtering (similarity spectrum, and filter based on max/min of Rrs)
       Note: `algorithmic filtering' is already applied to FP data on the geoserver (based on rho_s bounds)

Optional filters (commented out below):
       (i) Rrs shape based filtering (currently includes `coastal water filter' from Warren et al. 2019).
           Users are encouraged to apply filters for their local water type
       (ii) Variability filtering (based on z-score), following Groetsch et al. 2017.

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - Feb 2022
"""

import sys
import os
import numpy as np
from monda.sorad import access, plots, qc
import datetime
import logging
import pandas as pd
import argparse
from math import ceil

log = logging.getLogger('sorad-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def run_example(platform_id = 'PML_SR004',
                start_time = datetime.datetime(2021,10,21,0,0,0),
                end_time   = datetime.datetime(2021,10,22,23,59,59),
                bbox = None,
                target='.', rrsalgorithm='fp',
                output_radiance = False,
                output_metadata = False,
                output_rrs = False,
                output_plots = False):

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
               gps_speeds, tilt_avgs, tilt_stds = unpack_response(response, rrsalgorithm, wl_output)

        if output_plots:
            log.info("Creating (ir)radiance plots")
            plots.plot_ed_ls_lt(ed, ls, lt, time, wl_output, file_id, target)


        # Step (i) radiometric quality control filters (i.e. QC applied to l or e spectra)
        q_lt_ed = qc.qc_lt_ed_filter(ed, lt, time, wl_output, threshold = 0.020) # lt/Ed ratio (glint) filtering
        q_ed =    qc.qc_ed_filter(ed, min_ed_threshold = 500) # filters on ed and ls anomalies
        q_ls =    qc.qc_ls_filter(ls, wl_output, threshold = 1)
        q_rad =   qc.combined_filter(qc.combined_filter(q_lt_ed, q_ed), q_ls) # combined `radiometric' qc mask

        # Step (ii)  addtional qc metrics that apply to Rrs spectrum
        q_ss =        qc.qc_SS_NIR_filter(rrswl, rrs, upperthreshold = 3, lowerthreshold = 0.5)  # similarity spectrum
        q_maxrange =  qc.qc_rrs_maxrange(rrs, upperthreshold = 0.1, lowerthreshold = 0.00)    # filters on max and min rrs
        q_min =       qc.qc_rrs_min(rrs, rrswl)

        # (3c) algorithmic qc filters specfic to 3C (rmsd or termination at rho bounds)
        if rrsalgorithm == '3c':
            rmsd_3c =     np.array([response['result'][i]['c3_rmsd'] for i in range(len(response['result']))]) # rmsd values used in 3C residual filter
            rho_ds =      np.array([response['result'][i]['c3_rho_ds'] for i in range(len(response['result']))]) # rho factors
            rho_dd =      np.array([response['result'][i]['c3_rho_dd'] for i in range(len(response['result']))])
            rho_s =       np.array([response['result'][i]['c3_rho_s'] for i in range(len(response['result']))])
            q_rho =       qc.qc_3c_rho_filter(rho_ds, rho_dd, rho_s, upperbound = 0.1)
            q_rad_resid = qc.qc_3cresidual(q_rad, rmsd_3c, tol = 1.5)
            q_rad_3c =    qc.combined_filter(q_rho, q_rad_resid)

        if rrsalgorithm == '3c':
            q_rad_rrs = qc.combined_filter(q_rad_3c, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange)))) # recommended rrs qc mask for 3C method (combines step (i), (ii) and (iii) QC)
        elif rrsalgorithm == 'fp':
            q_rad_rrs = qc.combined_filter(q_rad,    qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange)))) # recommended rrs qc mask for FP method (combines step (i) and (ii) QC)

        # Optional filters
        # q_coastal = qc_coastalwater_rrsfilter(rrs, wl) #  filter based on expected shape of rrs - example from Warren 2019 used. users can input their own spectra here (will depend on water type)
        # q_var = qc_radiometric_variability(ed, lt, ls, time, wl, windowlength = 60, var_threshold =1.1, var_metric = 'zscore_max')

        if output_plots and rrsalgorithm == 'fp':
            log.info("Creating Rrs plots")
            plots.plot_rrs_qc_fp(rrs, time, rrswl, q_rad, q_rad_rrs,              file_id, target)
            plots.plot_coveragemap(lat, lon, q_rad, file_id, target, map_resolution = 11)
            plots.plot_results(ed, ls, wl_output, rrs, rrswl, time, q_rad_rrs, file_id, target)

        elif output_plots and rrsalgorithm == '3c':
            log.info("Creating Rrs plots")
            plots.plot_rrs_qc_3c(rrs, time, rrswl, q_rad, q_rad_3c, q_rad_rrs, file_id, target)
            plots.plot_coveragemap(lat, lon, q_rad_3c, file_id, target, map_resolution = 11)
            plots.plot_results(ed, ls, wl_output, rrs, rrswl, time, q_rad_rrs, file_id, target)


        d = pd.DataFrame()   # store core metadata and qc flags in a data frame for easy output formatting
        d['sample_uuid'] = sample_uuids
        d['platform_id'] = platform_ids
        d['platform_uuid'] = platform_uuids
        d['timestamp'] = time
        d['lat'] = lat
        d['lon'] = lon
        d['gps_speed'] = gps_speeds
        d['tilt_avg'] = tilt_avgs
        d['tilt_std'] = tilt_stds
        d['rel_view_az '] = rel_view_az
        d['q_1'] = q_rad     # Mask after step (i) QC

        if rrsalgorithm == 'fp':
            d['q_2'] = q_rad_rrs # Mask after step (ii) QC. Currently recommended for FP rrs data analysis

        elif rrsalgorithm == '3c':
            d['q_2'] = q_rad_3c   # Mask after step (ii) QC
            d['q_3'] = q_rad_rrs  # Mask after step (iii) QC: currently recommended for 3C  rrs data analysis

        # optional to output all qc masks
        # q_keys = [i for i in locals() if i.startswith('q_')]
        # for i in range(len(q_keys)): # add all qc fields to the dataframe
        #     d[str(q_keys[i])] = pd.Series(eval(q_keys[i]))

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


def unpack_response(response, rrsalgorithm, wl_out):
    """
    Unpack the WFS response
    """
    #log.info(response['result'][0].keys())   # uncomment to show all available fields

    time          = [response['result'][i]['time'] for i in range(len(response['result']))]
    lat           = np.array([response['result'][i]['lat'] for i in range(len(response['result']))])
    lon           = np.array([response['result'][i]['lon'] for i in range(len(response['result']))])
    rel_view_az   = np.array([response['result'][i]['rel_view_az'] for i in range(len(response['result']))])
    sample_uuid   = np.array([response['result'][i]['sample_uuid'] for i in range(len(response['result']))])
    platform_id   = np.array([response['result'][i]['platform_id'] for i in range(len(response['result']))])
    platform_uuid = np.array([response['result'][i]['platform_uuid'] for i in range(len(response['result']))])
    gps_speed     = np.array([response['result'][i]['id'] for i in range(len(response['result']))])
    tilt_avg      = np.array([response['result'][i]['tilt_avg'] for i in range(len(response['result']))])
    tilt_std      = np.array([response['result'][i]['tilt_std'] for i in range(len(response['result']))])

    ed = access.get_l1spectra(response, 'ed_', wl_out) # # irradiance spectra in 2D matrix format: rows time index, columns wavelength
    ls = access.get_l1spectra(response, 'ls_', wl_out)
    lt = access.get_l1spectra(response, 'lt_', wl_out)


    if rrsalgorithm == '3c':
        rrswl = np.arange(response['result'][0]['c3_wl_grid'][0], response['result'][0]['c3_wl_grid'][1], response['result'][0]['c3_wl_grid'][2])  # reconstruct wavelength grid for Rrs
        rrs = np.array([response['result'][i]['c3_rrs'][:] for i in range(len(response['result']))]) # 2D matrix format: rows time index, columns wavelength

    elif rrsalgorithm == 'fp':
        rrswl  = np.arange(response['result'][0]['wl_grid'][0], response['result'][0]['wl_grid'][1]-1, response['result'][0]['wl_grid'][2])  # reconstruct wavelength grid for Rrs
        rrs_    = np.array([response['result'][i]['rrs'][:] for i in range(len(response['result']))])  # rrs spectra 2D matrix format: rows time index, columns wavelength
        offset = np.array([response['result'][i]['nir_offset'] for i in range(len(response['result']))])
        rrs = np.array([rrs_[i,:] - np.ones(len(rrswl))*offset[i] for i in range(len(rrs_))]) # spectral offset (applied as default definition of FP rrs)

    return rrswl, time, lat, lon, rel_view_az, ed, ls, lt, rrs, sample_uuid, platform_id, platform_uuid, gps_speed, tilt_avg, tilt_std


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

    response = run_example(args.platform, args.start_time, args.end_time, args.bbox, args.target, args.algorithm.lower(),
                           args.output_radiance, args.output_metadata, args.output_rrs, args.output_plots)
