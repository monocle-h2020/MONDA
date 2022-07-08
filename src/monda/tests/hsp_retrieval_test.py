# -*- coding: utf-8 -*-
"""

Monda example script to retrieve HSP data from PML Geoserver using WFS standard,

---------------------------------------------------------------------------------

To see all command line options use:
    python hsp_retrieval_test.py -h

-------------------------------------------------------------------------------

Stefan Simis - stsi@pml.ac.uk - May 2022ll
Tom Jordan - tjor@pml.ac.uk - July 2022

"""

import sys
import os
import numpy as np
# from monda.hsp import access
path_hsp = '/users/rsg-new/tjor/monocle/monda_tests_hsp/MONDA/src/monda/hsp' 
sys.path.append(path_hsp) # temporary path fixes as monda.hsp did not work.
import access # temporary path fixes as monda.hsp did not work
import plots# temporary path fixes as monda.hsp did not work

import datetime
import logging
import argparse
import pandas as pd

log = logging.getLogger('sorad-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def run_example(sensor_id = None,
                start_time = datetime.datetime(2021,10,20,9,0,0),
                end_time   = datetime.datetime(2021,10,20,12,0,0),
                bbox = None,
                target='.', 
                output_irradiance = True,
                output_aot = True,
                output_metadata = True,
                output_irr_plots = True,
                output_aot_plots  = True,
                ):

    """
    Download data from a HSP instrument

    sensor_id: the serial number of a So-Rad platform, or None
    start_time:  date/time to start collecting data from
    end_time:    date/time to start collecting data from
    bbox:        Bounding box corner coordinates (lat,lon,lat,lon)
    """

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
       
       layer = 'rsg:hsp_public_view_full'

       response = access.get_wfs(sensor = sensor_id,
                                 timewindow = (datetime_i, datetime_e),
                                 layer=layer, bbox=bbox)

       log.info(f"{response['length']} features received.")

       if response['length'] == 0:
           continue
       
        # show first record
       for key, val in response['result'][0].items():
            log.info(f"{key}: {val}")
                 

       file_id = f"{datetime_i.strftime('%Y-%m-%d')}"   # labelling for output data files and plots

       time, lat, lon, sample_id, platform_id, ed, eds, IDR, wl  = unpack_response(response)   

       d = pd.DataFrame()   # store core metadata in a data frame for easy output formatting
       d['platform_id'] = platform_id
       d['sensor_id'] = platform_id
       d['timestamp'] = time
       d['lat'] = lat
       d['lon'] = lon
  
       if output_irradiance:    
           header = ",".join([str(w) for w in wl])
           ed_filename = os.path.join(target, file_id + '_Ed.csv')
           eds_filename = os.path.join(target, file_id + '_Eds.csv')      
           if os.path.exists(ed_filename):
              log.warning(f"File {ed_filename} was overwritten")
           if os.path.exists(eds_filename):
              log.warning(f"File {eds_filename} was overwritten")
           np.savetxt(ed_filename, ed, delimiter=',', header = header, fmt='%.8f')
           np.savetxt(eds_filename, eds, delimiter=',', header = header, fmt='%.8f')
           
       if output_aot:      
          solar_zenith = plots.calc_solar_zenith(time, lat, lon) # calculate olar zenith
          esolar = np.array(pd.read_table(path_hsp + '/SolarSpectrum.txt',skiprows=1)['Solar ET']) # solar spectrum
          tau_t, tau_a, tau_err = plots.calc_aot_direct(ed, eds, esolar, time, solar_zenith, wl, use_filter = True) # calculate optical thickness components

          header = ",".join([str(w) for w in wl])
          tau_t_filename = os.path.join(target, file_id + '_tau_t.csv') # total optical thickness
          tau_a_filename = os.path.join(target, file_id + '_tau_a.csv') # aerosol optical thickness
          if os.path.exists(tau_t_filename):
            log.warning(f"File {tau_t_filename} was overwritten")
          if os.path.exists(tau_a_filename):
            log.warning(f"File {tau_a_filename} was overwritten")
          np.savetxt(tau_t_filename, tau_t, delimiter=',', header = header, fmt='%.8f')
          np.savetxt(tau_a_filename, tau_a, delimiter=',', header = header, fmt='%.8f')
    
       if output_metadata: 
           if output_aot: 
              d['solar_zenith'] = solar_zenith # extra fields for AOT metadata
              d['tau_err'] = tau_err
           d_filename = os.path.join(target, file_id + '_metadata.csv')
           if os.path.exists(d_filename):
              log.warning(f"File {d_filename} was overwritten")
              d.to_csv(d_filename)

       if output_irr_plots:
           log.info("Creating direct and diffuse irradiance plots")
           plots.plot_irradiance_hsp(wl, eds, ed, IDR, time, file_id, target)

       if output_aot_plots:
           log.info("Creating aot plots")         
           plots.plot_aot(wl, tau_a, IDR, time, file_id, target)
    

    return response



def unpack_response(response):
    """
    Unpack the WFS response
    """
    #log.info(response['result'][0].keys())   # uncomment to show all available fields (core fields shown below)

    time          = [response['result'][i]['spectral_record_end_time'] for i in range(len(response['result']))] 
    lat           = np.array([response['result'][i]['lat'] for i in range(len(response['result']))])
    lon           = np.array([response['result'][i]['lon'] for i in range(len(response['result']))])
        
    sample_id   = np.array([response['result'][i]['sample_id'] for i in range(len(response['result']))])
    platform_id   = np.array([response['result'][i]['platform_id'] for i in range(len(response['result']))])
        
    ed = np.array([response['result'][i]['global_spectrum'][:] for i in range(len(response['result']))]) # 2D matrix format: rows time index, columns wavelength
    eds = np.zeros([len(ed),len(ed.T)]) # note: eds has 1 length wl bin than ed (it is missing a zero- padding digit so this is added manually). This fix will need checking for other deployments
    eds[:,1:] = np.array([response['result'][i]['diffuse_spectrum'][:] for i in range(len(response['result']))]) # 2D matrix format: rows time index, columns wavelength
      
    IDR =  np.array([response['result'][i]['diffuse_integral'] for i in range(len(response['result']))])/np.array([response['result'][i]['global_integral'] for i in range(len(response['result']))])
    wl = np.arange(response['result'][0]['wavelengths'][0], response['result'][0]['wavelengths'][1] + response['result'][0]['wavelengths'][2], response['result'][0]['wavelengths'][2]) 
     
    return time, lat, lon, sample_id, platform_id, ed, eds, IDR, wl 



def parse_args():
    """Interpret command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--sensor',      required = False, type = str, default = None, help = "Instrument serial number")
    parser.add_argument('-i','--start_time',  required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2021,10,20,9,0,0),
                                              help = "Initial UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-e','--end_time',    required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2021,10,20,12,0,0),
                                              help = "Final UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-b','--bbox',        required = False, type = float, nargs='+', default = None, help = "Restrict query to bounding box format [corner1lat corner1lon corner2lat corner2lon]")
    parser.add_argument('-t','--target',      required = False, type = str, default=None,
                                              help = "Path to target folder for plots (defaults to 'hsp_testoutput' in the current folder).")
    parser.add_argument('-r','--output_irradiance',  required = False, action='store_true', help = "Output Ed  and Eds spectra to csv file")
    parser.add_argument('-m','--output_metadata',  required = False, action='store_true', help = "Output metadata to csv file")
    parser.add_argument('-a','--output_aot',  required = False, action='store_true', help = "Output aerosol optical thickness to csv file")
    parser.add_argument('-p','--output_irr_plots',  required = False, action='store_true', help = "Output plots")
    parser.add_argument('-q','--output_aot_plots',  required = False, action='store_true', help = "Output plots")


    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()

    response = run_example(args.sensor, args.start_time, args.end_time, args.bbox)

    if args.target is None:
        args.target = os.path.join('.', 'HSP_test-output')
