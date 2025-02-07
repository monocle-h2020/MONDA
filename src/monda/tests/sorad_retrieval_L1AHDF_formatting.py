# -*- coding: utf-8 -*-
"""

Monda example script to retrieve So-Rad L0 and (ir)radiance spectra from PML Geoserver using WFS standard, and
translate to `L1A hdf format' for input to HyperCP.

To get the python modules, it is easiest to add monda (pypi) to the hypercp environment

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - Feb 2025
                            

"""

## imports from within monda
import sys
import os
import numpy as np
sys.path.append('/users/rsg/tjor/Monda_L0_expts/MONDA/src/monda/sorad')  # as with for 2024 changes, I have used a local import so than I can modify acces and plots
import qc 
import plots
import access
#from monda.sorad import access, plots, qc 
import datetime
import logging
import argparse

import pandas as pd
import json

log = logging.getLogger('sorad-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)

# imports from within HyperCP
sys.path.append('/users/rsg/tjor/HyperCP_Sorad/HyperCP')  #
from Source.HDFRoot import HDFRoot
from Source.HDFGroup import HDFGroup
from Source.Utilities import Utilities

from Source.ProcessL1aTriOS import ProcessL1aTriOS



def run_example(platform_id = 'PML_SR002',
                start_time = datetime.datetime(2023,10,10,0,0),
                end_time = datetime.datetime(2023,10,10,23,59,59),
                bbox = None,
                target='.',  
                rrsalgorithm='3c',
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
        layer = 'rsg:sorad_dev_l0_hypercp' # level 0 layer
        # layer = 'rsg:sorad_public_l0' # level 0 layer
    elif rrsalgorithm == '3c':
        # layer = 'rsg:sorad_public_l0' # level 0 layer
        layer = 'rsg:sorad_dev_l0_hypercp' # level 0 layer
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
    
        if response['length'] == 0:
            continue

        date_id = f"{datetime_i.strftime('%Y-%m-%d')}"   # labelling for output data files and plots

        ###################################
        #  L0-specific steps begin here
        ################################
        
        # Re-structure
        sensor_ids = ['SAM_874F', 'SAM_874E', 'SAM_874C'] 
        # THE SENSOR ORDER IS ASSUMED TO BE: ES (ed), LI (ls), LT (lt)
    
        # Retrieve L0 data from database
        time, lat, lon, rel_view_az,\
        ed, ls, lt, \
        ed_inttime, ls_inttime, lt_inttime, \
        sample_uuids, platform_ids, platform_uuids,\
        gps_speeds, tilt_avgs, tilt_stds = access.unpack_response_l0(response, rrsalgorithm, wl_output)
     
        
        ed = ed.astype(int)
        ls = ls.astype(int)
        lt = lt.astype(int)
        
   
        # Store core so-rad metadata and flags for easy formatting
        d = pd.DataFrame()   
        d['time'] = time
        d['lat'] = lat
        d['lon'] = lon
        d['sample_uuid'] = sample_uuids
        d['platform_id'] = platform_ids
        d['platform_uuid'] = platform_uuids
        d['gps_speed'] = gps_speeds
        d['tilt_avg'] = tilt_avgs
        d['tilt_std'] = tilt_stds
        d['rel_view_az'] = rel_view_az
        d['inttime_' + sensor_ids[0]] = ed_inttime
        d['inttime_' + sensor_ids[1]] = ls_inttime
        d['inttime_' + sensor_ids[2]] = lt_inttime
        

        # Hardcoded for now (wl coeffs are pulled from l1 wfs command)
        sensor_ids = ['SAM_874F', 'SAM_874E', 'SAM_874C'] # ES (ed), LI (ls), LT (lt)
  
        # paths to calfiles and HCP config
        # configPath = MainConfig.settings['cfgPath']
        # cal_path = configPath[0:configPath.rfind('.')] + '_Calibration/'
        cal_path = '/users/rsg/tjor/HyperCP_Sorad/HyperCP/Config/sample_TRIOS_sorad_Calibration/'
        config_path = '/users/rsg/tjor/HyperCP_Sorad/HyperCP/Config/sample_TRIOS_sorad.cfg' 
        # For each triplet, this creates an HDF - taken from TriOSL1A
       
        hours = [response['result'][i]['time'].hour for i in range(len(response['result']))]
        hours_of_sampling = np.unique(hours)

        # loop over hours of sampling in each day
        for h in range(len(hours_of_sampling)):        
    
            # sub-sample each hour of so-rad data
            d_h =  d[hours == hours_of_sampling[h]]
            ed_h = ed[hours == hours_of_sampling[h]]
            ls_h = ls[hours == hours_of_sampling[h]]
            lt_h = lt[hours == hours_of_sampling[h]]            
            
            # create dictionary for sensor data
            data = {sensor_ids[0]: ed_h,  sensor_ids[1]: ls_h,  sensor_ids[2]: lt_h}
    
            # create HDF root 
            root = HDFRoot()
            root.id = "/"
            
            # create HDF attributes (possible function?)
            root.attributes["WAVELENGTH_UNITS"] = "nm"
            root.attributes["LI_UNITS"] = "count"
            root.attributes["LT_UNITS"] = "count"
            root.attributes["ES_UNITS"] = "count"
            root.attributes["SATPYR_UNITS"] = "count"
            root.attributes["PROCESSING_LEVEL"] = "1a"
            root.attributes["RAW_FILE_NAME"] = "" # leave bmal
            # root.attributes["TIME-STAMP"] = a_name
            root.attributes["CAST"] = str(response['result'][1]['time'].date()) + '_' + str(hours_of_sampling[h]).zfill(2)
    
            # loop over sensor IDs to created HDF sensor-groups
            for sensor_id in sensor_ids:  # master-loop - add as function at later stage?
                
                with open(config_path, 'r', encoding="utf-8") as fc:
                    text = fc.read()
                    conf_json = json.loads(text)
                sensor = conf_json['CalibrationFiles'][sensor_id + '.ini']['frameType']
                print(sensor)
        
                gp =  HDFGroup()
                gp.id = sensor_id + '.ini'
                root.groups.append(gp)
                
                ProcessL1aTriOS.attr_ini(cal_path + sensor_id + '.ini', gp)
            
                # Reshape data - follows  fields in TriOSL1A
                n_samples = len(d_h) #
                rec_datetag  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) # defined in mlb: set to zero for now
                rec_datetag2  = ProcessL1aTriOS.reshape_data('NONE',n_samples, data=np.zeros(n_samples)) # defined in mlb: set to zero for now
                rec_inttime = ProcessL1aTriOS.reshape_data(sensor, n_samples, d_h['inttime_' +sensor_id])
                rec_check  = ProcessL1aTriOS.reshape_data('SUM', n_samples, data=np.zeros(n_samples))
                rec_darkave  = ProcessL1aTriOS.reshape_data(sensor, n_samples, data=np.zeros(n_samples))
                rec_darksamp  = ProcessL1aTriOS.reshape_data(sensor, n_samples, data=np.zeros(n_samples)) 
                rec_frame  = ProcessL1aTriOS.reshape_data('COUNTER', n_samples, data=np.zeros(n_samples)) 
                rec_posframe  = ProcessL1aTriOS.reshape_data('COUNT',n_samples, data=np.zeros(n_samples)) 
                rec_sample  = ProcessL1aTriOS.reshape_data('DELAY', n_samples, data=np.zeros(n_samples)) 
                rec_spectemp  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) 
                rec_thermalresp  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) 
                rec_time  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) 
                rec_timetag2  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples))  # defined in mlb: set to zero for now
        
                # 
                gp.attributes['CalFileName'] = 'SAM_'+  sensor_id  + '.ini'
                gp.addDataset('DATETAG')
                gp.datasets['DATETAG'].data=np.array(rec_datetag2, dtype=[('NONE', '<f8')])
                gp.addDataset('INTTIME')
                gp.datasets['INTTIME'].data=np.array(rec_inttime, dtype=[('NONE', '<f8')])
                gp.addDataset('CHECK')
                gp.datasets['CHECK'].data=np.array(rec_check, dtype=[('NONE', '<f8')])
                gp.addDataset('DARK_AVE')
                gp.datasets['DARK_AVE'].data=np.array(rec_darkave, dtype=[('NONE', '<f8')])
                gp.addDataset('DARK_SAMP')
                gp.datasets['DARK_SAMP'].data=np.array(rec_darksamp, dtype=[('NONE', '<f8')])
                gp.addDataset('FRAME')
                gp.datasets['FRAME'].data=np.array(rec_frame, dtype=[('NONE', '<f8')])
                gp.addDataset('POSFRAME')
                gp.datasets['POSFRAME'].data=np.array(rec_posframe, dtype=[('NONE', '<f8')])
                gp.addDataset('SAMPLE')
                gp.datasets['SAMPLE'].data=np.array(rec_sample, dtype=[('NONE', '<f8')])
                gp.addDataset('SPECTEMP')
                gp.datasets['SPECTEMP'].data=np.array(rec_spectemp, dtype=[('NONE', '<f8')])
                gp.addDataset('THERMAL_RESP')
                gp.datasets['THERMAL_RESP'].data=np.array(rec_thermalresp, dtype=[('NONE', '<f8')])
                gp.addDataset('TIMER')
                gp.datasets['TIMER'].data=np.array(rec_time, dtype=[('NONE', '<f8')])
                gp.addDataset('TIMETAG2')
                gp.datasets['TIMETAG2'].data=np.array(rec_timetag2, dtype=[('NONE', '<f8')])
        
        
                # Computing wavelengths (from ini file)
                c0 = float(gp.attributes['c0s'])
                c1 = float(gp.attributes['c1s'])
                c2 = float(gp.attributes['c2s'])
                c3 = float(gp.attributes['c3s'])
                wl = []
                for i in range(1,256):
                    wl.append(str(round((c0 + c1*(i+1) + c2*(i+1)**2 + c3*(i+1)**3), 2)))
               
                # Add level0 sensor data (irradiances/radiances) 
                ds_dt = np.dtype({'names': wl,'formats': [np.float64]*len(wl)})
                my_arr = np.array(data[sensor_id]).transpose()
                rec_arr = np.rec.fromarrays(my_arr, dtype=ds_dt)
                gp.addDataset(sensor)
                gp.datasets[sensor].data=np.array(rec_arr, dtype=ds_dt)
               
                

                
                # Add callibration data
                metacal, cal = ProcessL1aTriOS.read_cal(cal_path + 'Cal_' + sensor_id + '.dat')
                B1 = gp.addDataset('CAL_' + sensor)
                B1.columns["0"] = cal.values[:,1].astype(np.float64)
                B1.columnsToDataset()
                ProcessL1aTriOS.get_attr(metacal,B1)                
                breakpoint()
                metaback, back = ProcessL1aTriOS.read_cal(cal_path + 'Back_' + sensor_id + '.dat')
                C1 = gp.addDataset('BACK_'+ sensor)
                C1.columns["0"] = back.values[:,1]
                C1.columns["1"] = back.values[:,2]
                C1.columnsToDataset()
                
         
                breakpoint()
        

    return response

def parse_args():
    """Interpret command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--algorithm',   required = False, type=str, default = 'fp', help = "rrs processing algorithm: fp or 3c")
    parser.add_argument('-p','--platform',    required = False, type = str, default = 'PML_SR002', help = "Platform serial number, e.g. PML_SR004.")
    parser.add_argument('-i','--start_time',  required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2023,3,24,0,0,0),
                                              help = "Initial UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-e','--end_time',    required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2023,4,24,23,59,59),
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
        args.target = os.path.join('.', '.')

    if not os.path.isdir(args.target):
        os.mkdir(args.target)

    # response = run_example(args.platform, args.start_time, args.end_time, args.bbox, args.target, args.algorithm.lower(),
    # args.output_radiance, args.output_metadata, args.output_rrs, args.output_plots)

    response = run_example()
