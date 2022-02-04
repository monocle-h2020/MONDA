#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""

Monda access functions to retrieve So-Rad data from PML Geoserver using WFS standard.

Includes functions to convert level 0 spectra to level 1 spectra

------------------------------------------------------------------------------

Tom Jordan - tjor@pml.ac.uk - Feb 2022.

"""

import sys
sys.path.append("..")
import logging
import numpy as np
import datetime
import argparse
import json
import urllib.request
import urllib.parse


log = logging.getLogger('test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)  
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def get_wfs(count=100, platform=None, timewindow=None, layer='rsg:sorad_public_view_fp_rrs'): # note: previously in geoserver_functions.py
    """Get features in json format and convert data into python friendly format"""
    time_start = datetime.datetime.strftime(timewindow[0], '%Y-%m-%dT%H:%M:%SZ')
    time_end   = datetime.datetime.strftime(timewindow[1], '%Y-%m-%dT%H:%M:%SZ')

    # build CQL filter
    cql = ''
    if timewindow is not None:
        cql += f"time between {time_start} AND {time_end}"
    if platform is not None:
        if len(cql)>0:
            cql += " AND "
        cql += f"""platform_id='{platform}'"""

    # sanity-check the request count
    layer_limit = 10000  # sadly there isn't a clear way to request this from the server
    if count > layer_limit:
        log.warning(f"Feature request count adjusted to {layer_limit} as supported by this server")
        count = layer_limit

    # build base URL
    url = 'https://rsg.pml.ac.uk/geoserver/rsg/wfs'
    req = 'GetFeature'
    version = '2.0.0'
    service = 'WFS'
    typeName = layer
    srs = 'EPSG:4326'
    fmt = 'json'
    wfs_url = f"{url}?request={req}&version={version}&service={service}&typeName={typeName}"
    wfs_url += f"&count={count}&srsname={srs}&outputFormat={fmt}&sortBy=time+A"
    if cql:
        wfs_url += f"&CQL_FILTER={urllib.parse.quote(cql)}"

    # prepare to page results by first getting number of hits on query, without returning features
    try:
        resp = urllib.request.urlopen(wfs_url + "&resultType=hits").read()
    except urllib.request.HTTPError as err:
        log.error(f"""WFS Server <a href="{url}">{url}</a> error""")
        log.exception(err)
        return None

    hitsdata = str(resp)
    if 'numberMatched' in hitsdata:
        hits = int(hitsdata.split('numberMatched="')[1].split('"')[0])
        log.info(f"{hits} features matched")
    else:
        log.error(f"No features match the request")
        return None

    if count < hits:
        npages = int(np.ceil(hits/count))
        log.info(f"Need to page the request: {npages} pages")
    else:
        npages = 1

    paged_result = []
    for page in range(npages):
        startIndex = page * count
        paged_wfs_url = wfs_url + f"&startIndex={startIndex}"
        log.info(paged_wfs_url)  # change to debug later

        try:
            resp = urllib.request.urlopen(paged_wfs_url).read()

        except urllib.request.HTTPError as err:
            log.error(f"""WFS Server <a href="{url}">{url}</a> error""")
            log.exception(err)
            return None

        try:
            data = json.loads(resp)
        except json.JSONDecodeError as err:
            log.error(f"Could not decode JSON response")
            log.exception(err)
            return None

        # flatten the response and deal with special cases
        features = [None]*len(data['features'])
        for i, datum in enumerate(data['features']):
            feature = datum['properties']
            feature['lon'] = datum['geometry']['coordinates'][0]
            feature['lat'] = datum['geometry']['coordinates'][1]
            for key, val in feature.items():
                if isinstance(val, str):
                    if val.count(',') > 0:  # convert comma-separated strings to numpy arrays
                        feature[key] = np.array([float(w) for w in val.split(',')])
                if ('date' in key) or ('time' in key):  # try to convert into timestamp
                    try:
                        feature[key] = datetime.datetime.strptime(val, '%Y-%m-%dT%H:%M:%S.%fZ')
                    except ValueError:
                        try:
                            feature[key] = datetime.datetime.strptime(val, '%Y-%m-%dT%H:%M:%SZ')
                        except ValueError:
                            log.error(f"Error parsing date or time field: {val}")
                    except Exception as err:
                            log.exception(err)
            features[i] = feature

        log.info(f"Page {page}, starting at count {startIndex}: {len(features)} features")
        paged_result = paged_result + features

    return {'result':paged_result, 'length':len(paged_result)}


def parse_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--algorithm', required = False, type=str, default = 'fp', help = "rrs processing algorithm: fp or 3c")
    parser.add_argument('-p','--platform', required = False, type = str, default = 'PML_SR004', help = "Platform serial number, e.g. PML_SR004.")
    parser.add_argument('-i','--intial_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = datetime.date(2021,10,21) - datetime.timedelta(days=1),help = "Intial day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-f','--final_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = datetime.date(2021,10,22), help = "Intial day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-t','--target', required = False, type = str, default='SoRad_testoutput', help = "Target folder for plots to be written to (defaults to current folder).")
    parser.add_argument('-s','--spectra', required = False, type = bool, default = False, help = "Option to output ls, lt, ed spectra")
   
    args = parser.parse_args()
    
    return args


def get_wl(res, spec_id):
    """Reconstruct level 1 wavelength grid for level 0 l and e spectra based on calibration records"""

    nbands = len(res[spec_id + 'spectrum'])
    wave = [0.0] * nbands
    c0s = res[spec_id + 'wl'][0]
    c1s = res[spec_id + 'wl'][1]
    c2s = res[spec_id + 'wl'][2]
    c3s = res[spec_id + 'wl'][3]
    for i in range(1, len(wave)+1, 1):
        wave[i-1] = (c0s) + (c1s*(i+1)) +\
            (c2s*(i+1)**2) + (c3s*(i+1)**3)

    return wave


def get_l1spectra(response, spec_id, wl):
    """
    (Ir)radiance spectrum as a 2D matrix sampled on common wavelength grid

    response: response from WFS
    spec_id: 
    wl: wavelength grid

    output:
    ndarray with rows timestamp index, columns wavelength
    """

    spec_matrix = np.nan*np.ones([len(response['result']), len(wl)])
    i = 0
    for res in response['result']:
        spec = res[spec_id + 'spectrum']
        spec_wl = np.array(get_wl(res, spec_id))
        spec_matrix[i,:] = np.interp(wl, spec_wl, spec)
        i = i + 1

    return spec_matrix
