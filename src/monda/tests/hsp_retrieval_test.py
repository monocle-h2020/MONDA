# -*- coding: utf-8 -*-
"""

Monda example script to retrieve HSP data from PML Geoserver using WFS standard,

---------------------------------------------------------------------------------

To see all command line options use:
    python hsp_retrieval_test.py -h

-------------------------------------------------------------------------------

Stefan Simis - stsi@pml.ac.uk - May 2022
"""

import sys
import os
import numpy as np
from monda.hsp import access
import datetime
import logging
import argparse

log = logging.getLogger('sorad-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)

def run_example(sensor_id = None,
                start_time = datetime.datetime(2021,10,21,12,0,0),
                end_time   = datetime.datetime(2021,10,21,14,0,0),
                bbox = None):

    """
    Download data from a HSP instrument

    sensor_id: the serial number of a So-Rad platform, or None
    start_time:  date/time to start collecting data from
    end_time:    date/time to start collecting data from
    bbox:        Bounding box corner coordinates (lat,lon,lat,lon)
    """

    log.info(f"Request timeframe {start_time.isoformat()} - {end_time.isoformat()}")

    layer = 'rsg:hsp_public_view_full'

    response = access.get_wfs(sensor = sensor_id,
                              timewindow = (start_time, end_time),
                              layer=layer, bbox=bbox)


    log.info(f"{response['length']} features received.")

    if response['length'] == 0:
        return None



    # show first record
    for key, val in response['result'][0].items():
        log.info(f"{key}: {val}")


    return response


def parse_args():
    """Interpret command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--sensor',      required = False, type = str, default = None, help = "Instrument serial number")
    parser.add_argument('-i','--start_time',  required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2021,10,21,12,0,0),
                                              help = "Initial UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-e','--end_time',    required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S'),
                                              default = datetime.datetime(2021,10,21,14,0,0),
                                              help = "Final UTC date/time in format 'YYYY-mm-dd HH:MM:SS'")
    parser.add_argument('-b','--bbox',        required = False, type = float, nargs='+', default = None, help = "Restrict query to bounding box format [corner1lat corner1lon corner2lat corner2lon]")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()

    response = run_example(args.sensor, args.start_time, args.end_time, args.bbox)
