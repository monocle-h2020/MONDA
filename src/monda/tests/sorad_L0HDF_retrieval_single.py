#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is an example downloader for So-Rad data that has been structured
in HyperCP compatible `L0_HDF' files.

Sections of the code are harmonised with the `direct' L0_HDF So-Rad download.
https://github.com/monocle-h2020/so-rad/blob/master/bin/functions/download_functions.py

For now, this is a self

It is recommended that L0_HDF files in HyperCP are hourly.

Tom Jordan, Oct 2026, tjor@pml.ac.uk


"""

import os
import numpy as np
# from monda.sorad import access, plots, qc 
import datetime
import logging
from monda.sorad import access
#import argparse
import h5py # this is in the

log = logging.getLogger('download')

#
def filename_from_dates(platform_id, start_time, end_time, format='hdf'):
    
    """Generate filename from dates"""
    
    start_str = datetime.datetime.strftime(start_time, "%Y%m%dT%H%M%S")
    end_str =   datetime.datetime.strftime(end_time,   "%Y%m%dT%H%M%S")
    out_filepath = f"{platform_id}_{start_str}-{end_str}_L0.{format}"
    return out_filepath

# This could be moved to monda access
def get_l0spectra(response, spec_id):
    
    """
    Retrieves (Ir)radiance spectrum in digital counts at level 0 from WFS response
    
    Pads TriOS Gen 2 sensors (250 pixels on geoserver)
    to 256 pixels (as required for HyperCP input)
     
    """
    n_records = len(response['result'])
    n_pixels = len(response['result'][0]['l0_' + spec_id + '_spectrum'])
    
    if n_pixels == 256: # Gen 1 have 256 pixels
        spec_matrix = np.nan*np.ones([n_records, n_pixels]) 
        i = 0
        for i, res in enumerate(response['result']):
            spec = res['l0_' + spec_id + '_spectrum']
            spec_matrix[i,:] = spec
   
    elif n_pixels == 250: # Gen 2 have 250 pixels
        spec_matrix = np.nan*np.ones([n_records, n_pixels + 6]) # hardcoded (250 + 6 = 256)
        i = 0
        for i, res in enumerate(response['result']):
            spec = res['l0_' + spec_id + '_spectrum']
            spec_matrix[i,:-6] = spec # hardcoded for now
            spec_matrix[i,-6:] = 0 # we can't use NaN for no data; hence this is zero
    else:
        print('number of pixels is not 256 (G1) or 250 (G2)')
      
    spec_matrix = spec_matrix.astype(int) # convert floats to ints
        
    return spec_matrix


# This could be moved to monda access
def unpack_response_meta_L0(response):
    """
    Unpacks the WFS response metadata ready for L0 HDF formatting
    
    """

    time          = [response['result'][i]['time'] for i in range(len(response['result']))]
    lat           = np.array([response['result'][i]['lat'] for i in range(len(response['result']))])
    lon           = np.array([response['result'][i]['lon'] for i in range(len(response['result']))])
    rel_view_az   = np.array([response['result'][i]['rel_view_az'] for i in range(len(response['result']))])
    sample_uuid   = [response['result'][i]['sample_uuid'] for i in range(len(response['result']))]
    # platform_id   = np.array([response['result'][i]['platform_id'] for i in range(len(response['result']))]) 
    # platform_uuid = np.array([response['result'][i]['platform_uuid'] for i in range(len(response['result']))]) # not needed
    gps_speed     = np.array([response['result'][i]['gps_speed'] for i in range(len(response['result']))])
    tilt_avg      = np.array([response['result'][i]['tilt_avg'] for i in range(len(response['result']))])
    tilt_std      = np.array([response['result'][i]['tilt_std'] for i in range(len(response['result']))])

    return time, lat, lon, rel_view_az, sample_uuid, gps_speed, tilt_avg, tilt_std


def unpack_response_l0_sensor(response, spec_id):
    """
    Unpacks the WFS response level 0 sensor data. Includes group and frame fields 
    using TriOS SAM IDs
    
    Wavelength-pixel assignment is done within HyperCP
    """
    
    # names of sensor groups in L0 HDF
    sensor_group = response['result'][0]['l0_' + spec_id + '_SAM_id'] + '.ini'

    # SAM frame codes
    sensor_frame = response['result'][0]['l0_' + spec_id + '_SAM_id'].split('_')[1]
    
    # sensor inttime
    sensor_inttime = np.array([response['result'][i]['l0_' + spec_id + '_inttime'] for i in range(len(response['result']))])

    # sensor level 0 data
    sensor_l0 = get_l0spectra(response, spec_id) 
    
    return sensor_group, sensor_frame, sensor_inttime, sensor_l0


def save_to_hdf_from_GS(response, platform_id, destination_file):
    
    """
    Save records to a L0_hdf format (for ingestion by HyperCP). 
    
    Calls `unpack_response' functions which have already done some data 
    restructuring from the geoserver response
    
    """
    # create HDF root structure 
    f = h5py.File(destination_file, "w")

    # extract sorad metadata fields from geoserver response 
    record_time, latitude, longitude, rel_view_az, \
    sample_uuid, gps_speed, tilt_avg, tilt_std = unpack_response_meta_L0(response)
  
    # compute HyperCP datetag & timetag fields
    datetag2 = [float(datetime.datetime.strftime(record_time[i], '%Y%j')) for i in range(len(record_time))]
    timetag2 = [float(datetime.datetime.strftime(record_time[i], "%H%M%S.%f")[:-3]) for i in range(len(record_time))]
  
    # root attributes 
    # f.attrs["WAVELENGTH_UNITS"] = "nm" 
    f.attrs["LI_UNITS"] = "count"
    f.attrs["LT_UNITS"] = "count"
    f.attrs["ES_UNITS"] = "count"
    # f.attrs["SATPYR_UNITS"] = "count"  # include if needed, but no relation to So-Rad
    f.attrs["RAW_FILE_NAME"] = ""       # there is no upstream file
    f.attrs["PROCESSING_LEVEL"] = "0"
    start_time = record_time[0]
    f.attrs["CAST"] = datetime.datetime.strftime(start_time, "%Y%m%d_%H")
    f.attrs["TIME-STAMP"] = datetime.datetime.strftime(start_time, "%a %b %d %H:%M:%S %Y")

    # metadata group
    meta = f.create_group("sorad")
    meta.attrs['PLATFORM_ID'] = platform_id
    # meta.attrs['CalFileName'] = "n/a"
    meta.attrs['FrameType'] = 'Not Required'

    # add datasets
    meta.create_dataset('DATETAG', data=datetag2, dtype='f')
    meta.create_dataset('TIMETAG2', data=timetag2, dtype='f')
    meta.create_dataset('LATITUDE', data=latitude, dtype='f')
    meta.attrs['LATITUDE_UNITS'] = 'degrees'
    meta.create_dataset('LONGITUDE', data=longitude, dtype='f')
    meta.attrs['LONGITUDE_UNITS'] = 'degrees'
    meta.create_dataset('REL_AZ', data=rel_view_az, dtype='f')
    meta.attrs['REL_AZ_UNITS'] = 'degrees'
    meta.create_dataset('TILT', data=tilt_avg, dtype='f')
    meta.attrs['TILT_UNITS'] = 'degrees'
    meta.create_dataset('TILT_STD', data=tilt_std, dtype='f')
    meta.attrs['TILT_STD_UNITS'] = 'degrees'
    meta.create_dataset('GPS_SPEED', data=gps_speed, dtype='f')
    meta.attrs['GPS_SPEED_UNITS'] = 'm/s'
    meta.create_dataset('SAMPLE_UUID', data=sample_uuid, dtype = h5py.string_dtype())
    
    # extract sorad sensor L0 fields from geoserver response 
    ls_group, ls_frame, ls_intime, ls_l0 = unpack_response_l0_sensor(response, 'ls')
    ed_group, ed_frame, ed_intime, ed_l0 = unpack_response_l0_sensor(response, 'ed')
    lt_group, lt_frame, lt_intime, lt_l0 = unpack_response_l0_sensor(response, 'lt')
    
    # Sensor groups
    # naming convention: ES (ed), LI (ls), LT (lt)
    LI  = f.create_group(ls_group)
    LI.create_dataset('DATETAG',  data=datetag2, dtype='f')
    LI.create_dataset('TIMETAG2', data=timetag2, dtype='f')
    LI.attrs['FrameType'] = ls_frame
    LI.attrs['RadianceTerm1'] = 'LI'   # Satlantic naming legacy
    LI.attrs['RadianceTerm2'] = 'Ls'   # Gordon/Mobley naming legacy
    LI.create_dataset('L0', data = ls_l0, dtype='i8')
    LI.attrs['L0_units'] = 'count'
    LI.create_dataset('INTTIME', data=ls_intime, dtype='i8')
    LI.attrs['INTTIME_UNITS'] = 'ms'

    ES  = f.create_group(ed_group)
    ES.create_dataset('DATETAG', data=datetag2, dtype='f')
    ES.create_dataset('TIMETAG2', data=timetag2, dtype='f')
    ES.attrs['FrameType'] = ed_frame
    ES.attrs['RadianceTerm1'] = 'ES'   # Satlantic naming legacy
    ES.attrs['RadianceTerm2'] = 'Ed'   # Gordon/Mobley naming legacy
    ES.create_dataset('L0', data = ed_l0, dtype='i8')
    ES.attrs['L0_units'] = 'count'
    ES.create_dataset('INTTIME', data = ed_intime, dtype='i8')
    ES.attrs['INTTIME_UNITS'] = 'ms'

    LT  = f.create_group(lt_group)
    LT.create_dataset('DATETAG', data=datetag2, dtype='f')
    LT.create_dataset('TIMETAG2', data=timetag2, dtype='f')
    LT.attrs['FrameType'] = lt_frame
    LT.attrs['RadianceTerm1'] = 'LT'   # Satlantic naming legacy
    LT.attrs['RadianceTerm2'] = 'Lt'   # Gordon/Mobley naming legacy
    LT.create_dataset('L0', data = lt_l0 , dtype='i8')
    LT.attrs['L0_units'] = 'count'
    LT.create_dataset('INTTIME', data = lt_intime, dtype='i8')
    LT.attrs['INTTIME_UNITS'] = 'ms'

    # write to file
    print('Saving Level0 HDF file: ' + str(destination_file))
    f.attrs["L0_FILENAME"] = os.path.basename(destination_file)
    f.close()

if __name__ == '__main__':

    platform_id = 'PML_SR002'
    start_time = datetime.datetime(2024,8,15,11,0,0)
    end_time   = datetime.datetime(2024,8,15,11,59,59)
    layer_L0 = 'rsg:sorad_dev_l0_hypercp'    

    # WFS response
    response = access.get_wfs(platform = platform_id,
                              timewindow = (start_time, end_time),
                              layer = layer_L0,
                              bbox = None)
    

    # this is a temporary hard coding of SAM id field.
    # if these are added to the 'rsg:sorad_dev_l0_hypercp' layer, then the rest
    # of the functions can be used as they are
    n_samples = len(response['result'])
    for i in range(n_samples):
        response['result'][i]['l0_ed_SAM_id'] = 'SAM_874F'    
        response['result'][i]['l0_ls_SAM_id'] = 'SAM_874C'
        response['result'][i]['l0_lt_SAM_id'] = 'SAM_874E'

    # show first record
    for key, val in response['result'][0].items():
        print(f"{key}: {val}")

    first_time = response['result'][0]['time']
    last_time = response['result'][-1]['time']
    destination_file = filename_from_dates(platform_id, first_time, last_time, format='hdf')
    save_to_hdf_from_GS(response, platform_id, destination_file)
    
  