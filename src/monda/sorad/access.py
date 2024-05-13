#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Monda access functions to retrieve So-Rad data from PML Geoserver using WFS standard.
------------------------------------------------------------------------------
Tom Jordan - tjor@pml.ac.uk - Feb 2022
Stefan Simis - stsi@pml.ac.uk - Feb 2022

"""

import logging
import sys
import numpy as np
import datetime
import json
import urllib.request
import urllib.parse

log = logging.getLogger('sorad-downloader')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def get_wfs(count=1000, platform=None, timewindow=None, layer='rsg:sorad_public_view_fp_rrs', bbox=None):
    """
    Get features from geoserver layer, parse json response into python friendly format

    : int count:        Number of features to retrieve per request. There is typically a request feature limit for each geoserver layer. This enables paging the response. (default 1000)
    : str platform:     Optionally specify a sensor platform (default None)
    : tuple timewindow: Specify a time window as a tuple of datetime or iso formatted strings. The default None will request data from the last 24h
    : str layer:        Specify the geoserver layer. For a list of available layers see https://rsg.pml.ac.uk/geoserver/. (default rsg:sorad_public_view_fp_rrs)

    """
    if isinstance(timewindow, tuple):
        if (isinstance(timewindow[0], str)) and (isinstance(timewindow[1], str)):
            try:
                time_start = datetime.datetime.strftime(timewindow[0], '%Y-%m-%dT%H:%M:%SZ')
                time_end   = datetime.datetime.strftime(timewindow[1], '%Y-%m-%dT%H:%M:%SZ')
            except:
                log.error("timewindow should be a tuple of type datetime.datetime or iso-formatted string i.e. 'YYYY-mm-ddTHH:MM:SSZ'")
                return None
        elif (isinstance(timewindow[0], datetime.datetime)) and (isinstance(timewindow[1], datetime.datetime)):
            time_start = timewindow[0]
            time_end = timewindow[1]
        else:
            log.error("timewindow should be a tuple of type datetime.datetime or iso-formatted string i.e. 'YYYY-mm-ddTHH:MM:SSZ'")
            return None
    elif timewindow is None:
        # default to the last 24H
        time_start = datetime.datetime.now() - datetime.timedelta(days=1)
        time_end = datetime.datetime.now()

    # build a CQL filter to get a time slice
    cql = ''

    cql += f"time between {time_start.isoformat()} AND {time_end.isoformat()}"

    if platform is not None:
        if len(cql)>0:
            cql += f" AND platform_id='{platform}'"

    if bbox is not None:
        if not len(bbox) == 4:
           log.error("Bounding box expects a tuple of four values (two corner coordinate pairs)")
        if len(cql)>0:
            cql += " AND "
        cql += f"""BBOX(location, {",".join([str(b) for b in bbox])})"""

    # sanity-check the request count
    layer_limit = 10000  # geoserver layer is limited to 10,000 items per request
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
        hits_request = wfs_url + "&resultType=hits"
        log.info(f"Pre-paging request: {hits_request}")
        resp = urllib.request.urlopen(hits_request).read()
    except urllib.request.HTTPError as err:
        log.error(f"""WFS Server <a href="{url}">{url}</a> error""")
        log.exception(err)
        return None

    hitsdata = str(resp)
    if 'numberMatched' in hitsdata:
        hits = int(hitsdata.split('numberMatched="')[1].split('"')[0])
        log.info(f"{hits} features matched")
    else:
        log.error(f"No features match this request. The server response was {resp}")
        return None

    if count < hits:
        # the number of features matching the request exceeds the volume we have set to retreive in one go
        npages = int(np.ceil(hits/count))
        log.info(f"Paging the request: {npages} pages")
    else:
        npages = 1

    paged_result = []
    for page in range(npages):
        startIndex = page * count
        paged_wfs_url = wfs_url + f"&startIndex={startIndex}"
        log.debug(paged_wfs_url)

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


def get_wl(res, spec_id):
    """Reconstruct wavelength grid from information from most recent calibration records (only for L1 (ir)radiance spectra)"""

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


def get_l1spectra(response, spec_id, wl_out=np.arange(350, 951, 1)):
    """
    (Ir)radiance spectrum as a 2D matrix sampled on common wavelength grid

    response: response from WFS
    spec_id:  identifier
    wl:       output wavelength grid

    output:
    ndarray with rows timestamp index, columns wavelength
    """

    spec_matrix = np.nan*np.ones([len(response['result']), len(wl_out)])
    i = 0
    for i, res in enumerate(response['result']):
        spec = res[spec_id + 'spectrum']
        spec_wl = np.array(get_wl(res, spec_id))  # get the wavelength grid for this particular spectrum
        spec_matrix[i,:] = np.interp(wl_out, spec_wl, spec)  # interpolate to common wavelength grid

    return spec_matrix


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
    gps_speed     = np.array([response['result'][i]['gps_speed'] for i in range(len(response['result']))])
    tilt_avg      = np.array([response['result'][i]['tilt_avg'] for i in range(len(response['result']))])
    tilt_std      = np.array([response['result'][i]['tilt_std'] for i in range(len(response['result']))])

    ed = get_l1spectra(response, 'ed_', wl_out) # # irradiance spectra in 2D matrix format: rows time index, columns wavelength
    ls = get_l1spectra(response, 'ls_', wl_out)
    lt = get_l1spectra(response, 'lt_', wl_out)

    if rrsalgorithm == '3c':
        rrswl = np.arange(response['result'][0]['c3_wl_grid'][0], response['result'][0]['c3_wl_grid'][1], response['result'][0]['c3_wl_grid'][2])  # reconstruct wavelength grid for Rrs
        rrs = np.array([response['result'][i]['c3_rrs'][:] for i in range(len(response['result']))]) # 2D matrix format: rows time index, columns wavelength

    elif rrsalgorithm == 'fp':
        rrswl  = np.arange(response['result'][0]['wl_grid'][0], response['result'][0]['wl_grid'][1]-1, response['result'][0]['wl_grid'][2])  # reconstruct wavelength grid for Rrs
        rrs_    = np.array([response['result'][i]['rrs'][:] for i in range(len(response['result']))])  # rrs spectra 2D matrix format: rows time index, columns wavelength
        offset = np.array([response['result'][i]['nir_offset'] for i in range(len(response['result']))])
        rrs = np.array([rrs_[i,:] - np.ones(len(rrswl))*offset[i] for i in range(len(rrs_))]) # spectral offset (applied as default definition of FP rrs)

    return rrswl, time, lat, lon, rel_view_az, ed, ls, lt, rrs, sample_uuid, platform_id, platform_uuid, gps_speed, tilt_avg, tilt_std
