#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Monda access functions to retrieve Mini-Secchi disk from PML Geoserver using WFS standard.

To see all available command line options use
 python minisecchi_retrieval_test.py -h

Example use
 python minisecchi_retrieval_test.py -i "2020-01-01" -f "2021-12-31" -t ./secchitest.csv --bbox 50 20 40 30

This will retrieve all records from 2020-2021 within a bounding box spanning 50N20E - 40N30E and save results to secchitest.csv in the current folder.


------------------------------------------------------------------------------
Stefan Simis - stsi@pml.ac.uk - Feb 2022
"""

import os
import sys
import logging
import datetime
import argparse
from monda.minisecchi import access

log = logging.getLogger('secchi-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def parse_args():
    default_start = datetime.datetime.strptime('2022-01-01', '%Y-%m-%d')
    default_end = datetime.datetime.strptime('2022-07-20', '%Y-%m-%d')

    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--initial_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = default_start, help = "Intial day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-f','--final_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = default_end, help = "Intial day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-t','--target_file', required = False, type = str, default = None, help = "Target file path to store results in csv format")
    parser.add_argument('-b','--bbox', required = False, type = float, nargs='+', default = None, help = "Restrict query to bounding box format [corner1lat corner1lon corner2lat corner2lon]")
    args = parser.parse_args()

    return args


def run_example(start_date, end_date, geoserver_layer='rsg:miniscecchi_public_view', target=None, bbox=None):

    response = access.get_wfs(count=1000, timewindow=(start_date, end_date), layer=geoserver_layer, bbox=bbox)
    data = response['result']

    try:
        datafields = data[0].keys()
    except IndexError:
        print("no data found for time and region specified")
    else:
        header = ",".join([k for k in datafields])

    if target is not None:
        # attempt to write data to target file
        if os.path.isfile(target):
            log.warning(f"Existing data file will be overwritten")

        with open(target, 'w') as outfile:
            outfile.write(header + '\n')
            for i, record in enumerate(data):
                record_str = ",".join([str(v) for v in record.values()])
                if i < len(data):
                    outfile.write(record_str + '\n')
                else:
                    outfile.write(record_str)


if __name__ == '__main__':
    args = parse_args()
    run_example(args.initial_day, args.final_day, geoserver_layer='rsg:minisecchi_public_view', target=args.target_file, bbox=args.bbox)
