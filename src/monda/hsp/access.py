#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Monda access functions to retrieve HSP data from PML Geoserver using WFS standard.
------------------------------------------------------------------------------
Tom Jordan - tjor@pml.ac.uk - May 2022
Stefan Simis - stsi@pml.ac.uk - May 2022

"""

import logging
import sys
import numpy as np
import datetime
import json
import urllib.request
import urllib.parse

log = logging.getLogger('hsp-downloader')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def get_wfs(count=1000, sensor=None, timewindow=None, layer='rsg:hsp_public_view_full', bbox=None):
    """
    Get features from geoserver layer, parse json response into python friendly format

    : int count:        Number of features to retrieve per request. There is typically a request feature limit for each geoserver layer. This enables paging the response. (default 1000)
    : str sensor:       Optionally specify a sensor id (default None)
    : tuple timewindow: Specify a time window as a tuple of datetime or iso formatted strings. The default None will request data from the last 24h
    : str layer:        Specify the geoserver layer. For a list of available layers see https://rsg.pml.ac.uk/geoserver/. (default rsg:hsp_public_view_full)

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

    cql += f"spectral_record_end_time between {time_start.isoformat()} AND {time_end.isoformat()}"

    if sensor is not None:
        if len(cql)>0:
            cql += f" AND sensor_id='{sensor}'"

    if bbox is not None:
        if not len(bbox) == 4:
           log.error("Bounding box expects a tuple of four values (two corner coordinate pairs)")
        if len(cql)>0:
            cql += " AND "
        cql += f"""BBOX(spectral_record_location, {",".join([str(b) for b in bbox])})"""

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
    wfs_url += f"&count={count}&srsname={srs}&outputFormat={fmt}&sortBy=spectral_record_end_time+A"
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
                    if key in ['timezone']:
                        continue
                    if feature[key] is None:
                        continue
                    try:
                        feature[key] = datetime.datetime.strptime(val, '%Y-%m-%dT%H:%M:%S.%fZ')
                    except ValueError:
                        try:
                            feature[key] = datetime.datetime.strptime(val, '%Y-%m-%dT%H:%M:%SZ')
                        except ValueError:
                            log.error(f"Error parsing date or time field: {val}")
                    except TypeError:
                        pass
                    except Exception as err:
                            log.exception(err)
            features[i] = feature

        log.info(f"Page {page}, starting at count {startIndex}: {len(features)} features")
        paged_result = paged_result + features

    return {'result':paged_result, 'length':len(paged_result)}
