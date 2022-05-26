# -*- coding: utf-8 -*-
"""
Example/test downloading WISPstation data and optionally plotting Rrs and (ir)radiance curves and/or data download reports from WISPcloud

For full command line argument options see:
    python WISP_retrieval_plot_test.py -h

Example:
    python WISP_retrieval_plot_test.py -c -o -r -s '08:00:00' -e '16:00:00' -d '2019-06-21' -i 'WISPstation001' -t './test'


"""


from monda.WISP import plots
import os
import sys
import datetime
import logging
import argparse

log = logging.getLogger('WISP-test')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def run_example(instrument, date, start_time, end_time, target, output_radiance, output_rrs, output_data):
    log.info('running WISP test')

    if output_data:
        output_file_rrs = os.path.join(target, 'WISP_test_data_Rrs.txt')
        output_file_rad = os.path.join(target, 'WISP_test_data_Ed-Lt-Ls.txt')
    else:
        output_file_rrs = None
        output_file_rad = None

    plot1 = plots.plot_reflectance_station(instrument, date, start_time, end_time, output_file_rrs)
    plot2 = plots.plot_radiances_station(instrument, date, start_time, end_time, output_file_rad)

    if output_data:
        log.info(f"Rrs data saved to {output_file_rrs}")
        log.info(f"(ir)radiance data saved to {output_file_rad}")

    if plot1 is None:
       log.info("No data retrieved")
    elif output_rrs:
        plot_file_rrs = os.path.join(target, 'WISP_test_plot_Rrs.png')
        log.info(f"Saving Rrs plot to {plot_file_rrs}")
        plot1.savefig(plot_file_rrs, format='png', dpi=150)


    if plot1 is None:
       log.info("No data retrieved")
    elif output_radiance:
        plot_file_rad = os.path.join(target, 'WISP_test_plot_Ed-Lt-Ls.png')
        log.info(f"Saving (ir)radiances plot to {plot_file_rad}")
        plot2.savefig(plot_file_rad, format='png', dpi=150)


def parse_args():
    """Interpret command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--instrument',  required = False, type = str, default = 'WISPstation001', help = "WISP Instrument to retrieve data from (default WISPstation001)")
    parser.add_argument('-d','--date',        required = False, type = str, default = "2019-06-21", help = "Date to retrieve observation data from in format 'YYYY-mm-dd'")
    parser.add_argument('-s','--start_time',  required = False, type = str, default = "09:00:00", help = "First UTC Time to include in data retrieval in format 'HH:MM:SS' (default 09:00:00)")
    parser.add_argument('-e','--end_time',    required = False, type = str, default = "15:00:00", help = "Last UTC time to include in data retrieval in format 'HH:MM:SS' (default 15:00:00)")

    parser.add_argument('-t','--target',      required = False, type = str, default=None, help = "Path to target folder for plots (defaults to 'WISP_testoutput' in the current folder).")

    parser.add_argument('-r','--output_radiance',required = False, action='store_true', help = "Save plot of Ls, Lt, Ed spectra")
    parser.add_argument('-c','--output_rrs',     required = False, action='store_true', help = "Save plot of Rrs spectra")
    parser.add_argument('-o','--output_data',    required = False, action='store_true', help = "Save data reports")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()

    if not any([args.output_radiance, args.output_rrs]):
        log.warning("No plot outputs specified (see -h for help)")

    if args.target is None:
        args.target = os.path.join('.', 'WISP_test-output')

    if not os.path.isdir(args.target):
        os.mkdir(args.target)

    run_example(args.instrument, args.date, args.start_time, args.end_time,
                args.target, args.output_radiance, args.output_rrs, args.output_data)
