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
import h5py

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


def get_date_time_tags(d):
    'gets date and time tags consistent with mlb file definitions,'
    'datetag: fractional days since 1900'
    'datetag2: YYYYJJJ,  (year-jday) float'
    'timetag HMMSSmSuSmS (hour-minute-millisec) float'
    
    start_time = datetime.datetime(1900,1,1,0,0,0)
    
    datetags   = [] 
    datetag2s  = []
    timetags   = []
    
    for i in range(len(d)):  
        
          timestamp_i = datetime.datetime.strptime(d['timestamp'][i], '%Y-%m-%d %H:%M:%S.%f')
          
          # fractional days since 1900
          delta_t_i = timestamp_i - start_time 
          
          #
          datetag_i = delta_t_i.days + delta_t_i.seconds/(60*60*24) + delta_t_i.microseconds/(10**6) 
          datetag2_i = float(str(timestamp_i.year) + str(timestamp_i.timetuple().tm_yday))  
          timetag_i = float(str(timestamp_i.hour) + str(timestamp_i.minute).zfill(2) + str(timestamp_i.second).zfill(2) + str(timestamp_i.microsecond)[0:3])
          
          datetags.append(datetag_i)
          datetag2s.append(datetag2_i)
          timetags.append(timetag_i)

    d['datetag'] = np.array(datetags)
    d['datetag2'] = np.array(datetag2s)
    d['timetag'] = np.array(datetags)

    return d
          

def init_attributes(root, response, hour):
    'Itialise HDF attributes'
    
    root.attributes["WAVELENGTH_UNITS"] = "nm"
    root.attributes["LI_UNITS"] = "count"
    root.attributes["LT_UNITS"] = "count"
    root.attributes["ES_UNITS"] = "count"
    root.attributes["SATPYR_UNITS"] = "count"
    root.attributes["PROCESSING_LEVEL"] = "1a"
    root.attributes["RAW_FILE_NAME"] = "" # leave bmal
    root.attributes["CAST"] = str(response['result'][1]['time'].date()) + '_' + str(hour).zfill(2)
    start_datetime = datetime.datetime.strptime(str(response['result'][0]['time'])[0:19],'%Y-%m-%d %H:%M:%S')
    root.attributes["TIME-STAMP"] = datetime.datetime.strftime(start_datetime,'%a %b %d %H:%M:%S %Y')

    print('root attributes added to HDF')
        

    return


def init_sorad_group(root, d_h):
    
    'Itialise HDF group for sorad ancillary data'
        
    # create sorad group
    gp =  HDFGroup()
    gp.id = 'sorad'
    root.groups.append(gp)
    
    # platform fields are stored as attributes
    gp.attributes['PLATFORM_ID'] = d_h['platform_id'][0]
    gp.attributes['PLATFORM_UUID'] = d_h['platform_uuid'][0]

    # all other fields are stored as datasets
    n_samples = len(d_h) #

    # datetag and timetag (duplicating formatting steps from TriOS L1A)
    rec_datetag  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=d_h['datetag']) # defined in TriOS L1a
    rec_datetag2  = ProcessL1aTriOS.reshape_data('NONE',n_samples, data=d_h['datetag2']) # defined in TriOS L1a
    rec_timetag2  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=d_h['timetag'])  # defined in TriOS L1a
   
    gp.addDataset('DATETAG')
    gp.datasets['DATETAG'].data=np.array(rec_datetag2, dtype=[('NONE', '<f8')])
    gp.addDataset('TIMETAG2')
    gp.datasets['TIMETAG2'].data=np.array(rec_timetag2, dtype=[('NONE', '<f8')])

    #  
    gp.addDataset('LATPOS')
    gp.datasets['LATPOS'].data=np.array(d_h['lat'].values, dtype=[('NONE', '<f8')])
    gp.addDataset('LONPOS')
    gp.datasets['LONPOS'].data=np.array(d_h['lon'].values, dtype=[('NONE', '<f8')])
    gp.addDataset('TILT')
    gp.datasets['TILT'].data=np.array(d_h['tilt_avg'].values, dtype=[('NONE', '<f8')])
    gp.addDataset('RELATIVE_AZIMUTH')
    gp.datasets['RELATIVE_AZIMUTH'].data=np.array(d_h['rel_view_az'].values, dtype=[('NONE', '<f8')])
    # gp.addDataset('SORAD_UUID') # STATION == UUID 
    # gp.datasets['SORAD_UUID'].data=np.array(d_h['sample_uuid'].values, dtype=[('NONE', '<U36')])

    gp.addDataset('TILT_STD') # POTENTIALLY NOT NEEDED?
    gp.datasets['TILT_STD'].data=np.array(d_h['tilt_std'].values, dtype=[('NONE', '<f8')])
    gp.addDataset('GPS_SPEED') # POTENTIALLY NOT NEEDED?
    gp.datasets['GPS_SPEED'].data=np.array(d_h['gps_speed'].values, dtype=[('NONE', '<f8')])           


    print('sorad group added to HDF')
    
    return

def init_sensor_group(root, l0_data, sensor_id, config_path, cal_path, d_h):
    
    'Itialise HDF group for each TriOS sensor - this largely follows format of' 
    'eProcessL1TriOS.py'
  
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
  
    d_h  =  get_date_time_tags(d_h)
 
    rec_datetag  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=d_h['datetag']) # defined in mlb: # we use different date convention
    rec_datetag2  = ProcessL1aTriOS.reshape_data('NONE',n_samples, data=d_h['datetag2']) # defined in mlb: # we use different time convention
    rec_inttime = ProcessL1aTriOS.reshape_data(sensor, n_samples, d_h['inttime_' + sensor_id])
    rec_check  = ProcessL1aTriOS.reshape_data('SUM', n_samples, data=np.zeros(n_samples))
    rec_darkave  = ProcessL1aTriOS.reshape_data(sensor, n_samples, data=np.zeros(n_samples))
    rec_darksamp  = ProcessL1aTriOS.reshape_data(sensor, n_samples, data=np.zeros(n_samples)) 
    rec_frame  = ProcessL1aTriOS.reshape_data('COUNTER', n_samples, data=np.zeros(n_samples)) 
    rec_posframe  = ProcessL1aTriOS.reshape_data('COUNT',n_samples, data=np.zeros(n_samples)) 
    rec_sample  = ProcessL1aTriOS.reshape_data('DELAY', n_samples, data=np.zeros(n_samples)) 
    rec_spectemp  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) 
    rec_thermalresp  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) 
    rec_time  = ProcessL1aTriOS.reshape_data('NONE', n_samples, data=np.zeros(n_samples)) 
    rec_timetag2  = ProcessL1aTriOS.reshape_data('NONE', n_samples, d_h['timetag'])  # defined in mlb: set to zero for now

    # note - the majority of these datasets are set to zero 
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
    my_arr = np.array(l0_data[sensor_id]).transpose()
    rec_arr = np.rec.fromarrays(my_arr, dtype=ds_dt)
    gp.addDataset(sensor)
    gp.datasets[sensor].data=np.array(rec_arr, dtype=ds_dt)
      
    # Add callibration data
    metacal, cal = ProcessL1aTriOS.read_cal(cal_path + 'Cal_' + sensor_id + '.dat')
    B1 = gp.addDataset('CAL_' + sensor)
    B1.columns["0"] = cal.values[:,1].astype(np.float64)
    B1.columnsToDataset()
    ProcessL1aTriOS.get_attr(metacal,B1)                
   
    metaback, back = ProcessL1aTriOS.read_cal(cal_path + 'Back_' + sensor_id + '.dat')
    C1 = gp.addDataset('BACK_'+ sensor)
    C1.columns["0"] = back.values[:,1]
    C1.columns["1"] = back.values[:,2]
    C1.columnsToDataset()
    
    print(str(sensor_id) +' ('  + str(sensor) + ')  added to HDF')
    
    return


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
        
        # For now, hardcoded   
        sensor_ids = ['SAM_874F', 'SAM_874E', 'SAM_874C'] # ES (ed), LI (LS), LT (lt)
        cal_path = '/users/rsg/tjor/HyperCP_Sorad/HyperCP/Config/sample_TRIOS_sorad_Calibration/'
        config_path = '/users/rsg/tjor/HyperCP_Sorad/HyperCP/Config/sample_TRIOS_sorad.cfg' 
        output_path = '/users/rsg/tjor/Monda_L0_expts/MONDA/src/monda/HDF_output/' # temporay fix 
 
        # Retrieve L0 data from database
        time, lat, lon, rel_view_az,\
        ed, ls, lt, \
        ed_inttime, ls_inttime, lt_inttime, \
        sample_uuid, platform_id, platform_uuid,\
        gps_speed, tilt_avg, tilt_std = access.unpack_response_l0(response, rrsalgorithm, wl_output)
     
        # Store core Sorad metadata, including sensor integration times in a dataframe 
        d = access.meta_l0_dataframe(sample_uuid, platform_id, platform_uuid, time, lat, lon, gps_speed, tilt_avg, tilt_std, rel_view_az, ed_inttime, ls_inttime, lt_inttime, sensor_ids)
      
        # hours of sampling within day
        hours = [response['result'][i]['time'].hour for i in range(len(response['result']))]
        hours_of_sampling = np.unique(hours)
        
        
        outFFP = [] # list for output file path
        for h in range(len(hours_of_sampling)): # loop over hours of sampling in each day to create hourly hdf files
     
            # sub-sample each hour of so-rad data (indexed by hour h)
            d_h =  d[hours == hours_of_sampling[h]].reset_index(drop=True)
            ed_h = ed[hours == hours_of_sampling[h]]
            ls_h = ls[hours == hours_of_sampling[h]]
            lt_h = lt[hours == hours_of_sampling[h]]      
   
            # get date and time flags following conventions in ProcessTriOSL1A
            d_h = get_date_time_tags(d_h) 
   
            # create HDF root structure
            root = HDFRoot()
            root.id = "/"
            
            # create HDF attributes
            init_attributes(root, response, hours_of_sampling[h])

            # create sorad group
            init_sorad_group(root, d_h)
        
            # create sensor l0 groups
            l0_data = {sensor_ids[0]: ed_h,  sensor_ids[1]: ls_h, sensor_ids[2]: lt_h} # l0 data in dict
            for s in range(len(sensor_ids)):  # loop over sensor types
                init_sensor_group(root, l0_data, sensor_ids[s], config_path, cal_path, d_h)
         
            # write to file
            hdf_filename =  str(platform_id[0]) + '_' + f"{datetime_i.strftime('%Y-%m-%d')}"  + '_' + str(hours_of_sampling[h]).zfill(2) + '_L1A.hdf'
    
         
            outFFP.append(os.path.join(output_path, hdf_filename))
            root.attributes["L1A_FILENAME"] = outFFP[-1]
        
            try:
                # root.writeHDF5(new_name)
                root.writeHDF5(outFFP[-1])

            except Exception:
                msg = 'Unable to write L1A file. It may be open in another program.'
                Utilities.errorWindow("File Error", msg)
                print(msg)
                Utilities.writeLogFile(msg)
                return None, None
                
            d_h = []
            ed_h = []
            ls_h = []
            lt_h = []
   

    return 

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
