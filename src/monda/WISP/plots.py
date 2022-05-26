import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np

from monda.WISP import access

import sys
import logging

log = logging.getLogger('WISP-access')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)



def list_to_array(lstring):
   try:
        arr = np.array(lstring.lstrip('[').rstrip(']').split(',')).astype(np.float64)
        return arr
   except:
        return None

def plot_reflectance_station(instrument, day, start, stop, output_rrs):
    ''' Plot reflectance spectra for one WISPstation for one day
        Retrieve data from WISPcloud API
    '''
    style.use('seaborn-whitegrid')
    REQUEST = 'REQUEST=GetData&INSTRUMENT={}&include=measurement.date,measurement.id,level2.reflectance,site.name,level2.quality&TIME={}T{},{}T{}'.format(
        instrument, day, start, day, stop)

    l2r = access.WISP_data_API_call(REQUEST, output_rrs)
    rows = len(l2r)
    log.info(f"{rows} records retrieved")

    if rows == 0:
        return None

    color = iter(plt.cm.viridis(np.linspace(0, 1, rows)))
    wl = np.array(list(range(350, 901, 1)))
    fig = plt.figure(figsize=(8, 6))
    for meas in l2r:
        #print(meas.keys())
        Rrs = list_to_array(meas['level2.reflectance'])
        if Rrs is not None:
            c = next(color)
            plt.plot(wl, Rrs, c=c, label=meas['measurement.date'][10:16])
            station = meas['site.name']
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Rrs [sr-1]')
    plt.xlim(wl[0], wl[-1])
    plt.legend(loc='upper left', prop={'size': 8}, bbox_to_anchor=(1, 1))
    plt.title('{} at {} on {}'.format(instrument, station, day))
    plt.tight_layout(pad=7)

    return fig


def plot_radiances_station(instrument, day, start, stop, output_rad):
    ''' Plot (ir)radiance spectra for one WISPstation for one day
        Retrieve data from WISPcloud API
    '''
    style.use('seaborn-whitegrid')
    REQUEST = 'REQUEST=GetData&INSTRUMENT={}&include=measurement.date,measurement.id,site.name,ed.irradiance,ld.radiance,lu.radiance,level2.quality&TIME={}T{},{}T{}'.format(
        instrument, day, start, day, stop)

    l1b = access.WISP_data_API_call(REQUEST, output_rad)
    rows = len(l1b)
    log.info(f"{rows} records retrieved")

    if rows == 0:
        return None

    color = iter(plt.cm.Set2(np.linspace(0, 1, rows)))
    wl = np.array(list(range(350, 901, 1)))
    fig = plt.figure(figsize=(10, 8))
    for meas in l1b:
        Ed = list_to_array(meas['ed.irradiance'])
        Lu = list_to_array(meas['lu.radiance'])
        Ld = list_to_array(meas['ld.radiance'])
        if Ed is not None and Lu is not None and Ld is not None:
            #plot by quality ? Currently only read l2 quality (but this is a nice top level check)
            c = next(color)
            plt.subplot(311)
            plt.plot(wl, Ed, c=c, label=meas['measurement.date'][10:16])
            plt.subplot(312)
            plt.plot(wl, Ld, c=c)
            plt.subplot(313)
            plt.plot(wl, Lu, c=c)
            station = meas['site.name']
    plt.subplot(311)
    plt.ylabel('Ed [W/(m2*nm)]')
    plt.xlim(wl[0], wl[-1])
    plt.title('{} at {} on {}'.format(instrument, station, day))
    plt.legend(loc='upper left', prop={'size': 8}, bbox_to_anchor=(1.01, 1))
    plt.subplot(312)
    plt.ylabel('Ld [W/(m2*nm*sr)]')
    plt.xlim(wl[0], wl[-1])
    plt.subplot(313)
    plt.ylabel('Lu [W/(m2*nm*sr)]')
    plt.xlim(wl[0], wl[-1])
    plt.xlabel('Wavelength [nm]')
    return fig
