#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Monda access functions to retrieve Mini-Secchi disk from PML Geoserver using WFS standard.
------------------------------------------------------------------------------
Stefan Simis - stsi@pml.ac.uk - Feb 2022
"""

import logging
import datetime
import argparse
from monda.minisecchi import secchi_access

log = logging.getLogger('secchi-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def parse_args():
    default_start = datetime.datetime.now().date() - datetime.timedelta(days=7)
    default_end = datetime.datetime.now().date()

    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--initial_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = default_start, help = "Intial day in format ['yyyy-mm-dd]'] ")
    parser.add_argument('-f','--final_day', required = False, type = lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default = default_end, help = "Intial day in format ['yyyy-mm-dd]'] ")

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()
    response = secchi_access.get_wfs(count=100, timewindow=(args.initial_day, args.final_day), layer='rsg:minisecchi_public_view')
    data = response['result']
    log.info(data)
