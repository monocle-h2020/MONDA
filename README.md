# MONDA (MONocle Data Analysis)

![monocle_logo](https://avatars1.githubusercontent.com/u/36449994?s=200&v=4)

## Package Description
This package contains a suite of tools for retrieving, apply quality checks to, analysing and plotting data from the 
sensors and platforms included in the MONOCLE observation network. This project has received funding from the 
European Union’s Horizon 2020 research and innovation programme under grant agreement No 776480

The MONOCLE project created a framework for building water quality sensor and platforms, networked to enhance the 
utility and accessibility of data from multiple sources, giving a more complete data landscape to support satellite
observation of water quality in optically complex coastal waters, lakes and estuaries.

For more information on the MONOCLE project see the [project website](https://www.monocle-h2020.eu/Home)

## Dependencies
This code requires:
- Python (>= 3.8)
- NumPy (>= 1.13.3)
- scikit-learn(>=0.23.2)
- Matplotlib (>=3.3.3)
- requests (>=2.27.1)
- cartopy (>=0.20.2)
## Installation
**NOTE:** Some users have encountered issues installing into a fresh conda environment with pip due to GEOS 
versions.  This can be solved by installing cartopy with conda (```conda install cartopy```) before installing monda using pip.

```pip install monda```

Example creating MONOCLE conda environment and then installing monda package:
```
conda create --name monocle_test python=3.8 cartopy 
conda activate monocle_test
pip install monda 
```
# Source code
To get the most up to date version of the source code please see the repository at:
```
https://github.com/monocle-h2020/MONDA.git
```

## Citation
If you use MONDA in a scientific publication, we would appreciate citations.

To cite the package as a whole you can use:
Simis, S., Jackson, T., Jordan, T., Peters, S., and Ghebrehiwot, S. (2022) Monda: Monocle Data Analysis python package, https://github.com/monocle-h2020/MONDA 

For single submodules (such as WISP or sorad) please use:
[submodule] In Simis, S., Jackson, T., Jordan, T., Peters, S., and Ghebrehiwot, S. (2022) Monda: Monocle Data Analysis python package, https://github.com/monocle-h2020/MONDA

## Contributors
This code was developed with input from Plymouth Marine laboratory (thja-pml@github, tjor@github, StefanSimis@github) and 
Water Insight (Semhar-Ghe@github, waterthing@github). 

## Submodule Information
The package contains access, quality control and visualisation tools for a number of sensor systems, for which details are provided below.

## WISP (station)
The WISPstation is a fixed position optical instrument used for measuring water-leaving reflectance.
It records radiance and irradiance with an extended wavelength range of 350nm to 1100nm in two viewing directions,
which enables continuous and autonomous high-quality measurements for water quality monitoring and satellite validation. 
The reflectance observations are used to validate satellite measurements of water-leaving reflectance. 
Concentrations of the most important bio-physical water quality parameters such as chlorophyll-a, cyanobacterial pigment, turbidity and suspended matter, are derived from the reflectance measurement. The WISPstation sends the measurements automatically over 3G/4G/5G to the “WISPcloud” cloud database which makes the results available via an API. Measurement frequency is by default a 15 min interval but be adjusted to suit user requirements. 

#### About WISPcloud
WISPcloud is a scalable Postgres database that autonomously receives, stores, performs quality control and 
applies water quality algorithms to all WISPstation measurements. It has an advanced API to serve data requests directly to customers. A separate online documentation can be found here. 

#### Acknowledgement 
The WISPstation public data were collected by users participating on H2020 funded projects such as EOMORES(http://eomores-h2020.eu), TAPAS(http://tapas-h2020.eu/) and MONOCLE(https://monocle-h2020.eu/). 

#### Example data availability
Please use the instrument identification serial number and date when searching for data using the WISPcloud API
 
| Instrument ID  | Country   | Station         | Longitude | Latitude | Start Date | End Date   |
|----------------|-----------|-----------------|-----------|----------|------------|------------|
| WISPstation001 | Italy     | Lake Trasimeno  | 12.344    | 43.1223  | 2018-04-30 | 2018-10-14 |
| WISPstation001 | Italy     | Lake Trasimeno  | 12.344    | 43.1223  | 2019-06-20 | 2021-05-04 |
| WISPstation004 | Greece    | Souda           | 24.1112   | 35.4800  | 2018-07-17 | 2019-08-09 |
| WISPstation005 | Estonia   | Lake Vortsjarv  | 26.1074   | 58.2109  | 2018-05-28 | 2018-10-26 |
| WISPstation005 | Estonia   | Lake Vortsjarv  | 26.1074   | 58.2109  | 2019-05-31 | 2019-11-01 |
| WISPstation006 | Lithuania | Curonian Lagoon | 21.1002   | 55.4126  | 2018-08-09 | 2019-10-14 |
| WISPstation007 | Lithuania | Klaipeda Harbor | 21.1016   | 55.7195  | 2018-08-13 | 2019-09-11 |
| WISPstation009 | Hungary   | Lake Balaton    | 17.8936   | 46.9143  | 2019-06-17 | 2019-07-12 |
| WISPstation009 | Hungary   | Halasto         | 17.6167   | 46.6342  | 2019-07-23 | 2019-10-07 |
 
#### Functionality of the submodule
An example script is provided to connect with the WISPcloud API and subsequently plot Rrs and (ir)radiance measurements using date and instrument serial number as input arguments. 

#### Minimum code example 
```
from monda.WISP import access
import numpy as np

instrument = "WISPstation001"
day = "2019-08-16"
start = "10:00:00"
stop = "15:00:00"

def list_to_array(lstring):
    try:
        arr = np.array(lstring.lstrip('[').rstrip(']').split(',')).astype(np.float64)
        return arr
    except:
        return None

REQUEST = 'REQUEST=GetData&INSTRUMENT={}&include=measurement.date,measurement.id,level2.reflectance,site.name,level2.quality&TIME={}T{},{}T{}'\
          .format(instrument, day, start, day, stop)

l2r = access.WISP_data_API_call(REQUEST)

print(l2r[0]) # data header

rrs = [list_to_array(meas['level2.reflectance']) for meas in l2r[1:]]
wl = np.array(list(range(350, 901, 1)))
stations = [meas['site.name'] for meas in l2r[1:]]
times = [meas['measurement.date'][10:16] for meas in l2r[1:]]
```


## Solar-Tracking Radiometry platform (So-Rad)
The So-Rad is a low-power, low cost autonomous platform to obtain high-frequency water-leaving reflectance from 
non-stationary platforms such as ships and buoys. So-Rad software is highly configurable and open-source. 
So-Rad optimizes the measurement geometry of commercially available sensors which increases the number of successful 
observations of water colour obtained from moving platforms (concept as in Simis and Olsson 2013). 

Hyperspectral water-leaving reflectance is used to determine diagnostic features in water colour that can be 
associated with phytoplankton biomass, suspended solids and dissolved organic matter concentration.

Observing in situ reflectance with sensors on the So-Rad is used to validate satellite observations, particularly the 
performance of algorithms that separate atmospheric and water-leaving radiance, which have high uncertainty in 
optically complex waters such as coastal seas and inland waters.  High-quality reference measurements are required, 
collected under optimal observation conditions (solar and viewing azimuth, sun elevation).

#### Added Value of So-Rad
* Off-shore satellite validation is currently limited to research vessels and fixed moorings that are costly to 
maintain. The So-Rad can be installed on non-stationary platforms and is ideally suited to be included on merchant 
vessels. Ferry routes are recommended because of predictable routes and schedules. Periodic sensor maintenance can be 
easily carried out by non-expert crew.
* A high degree of automation and low-power components means the platform can be installed in remote locations for 
autonomous operation.

#### Functionality of the submodule
The test script provided demonstrates how to download paged data from the So-Rad Geoserver layers hosted at PML. 
These layers offer unfiltered, calibrated (ir)radiance and reflectance spectra. The reflectance data are processed either with the Fingerprint or the 3C method. 
Subsequently, quality control filters can be applied and data visualized. The scripts allow downloads per time window and per instrument. 

#### Minimum code example

```
from monda.sorad import access, plots, qc
import datetime
import numpy as np

# collect data from any So-Rad platform obtained in the past 24h (adjust time range if this returns no data)
start_time =   datetime.datetime.now() - datetime.timedelta(days=1)
end_time   =   datetime.datetime.now()

platform_id = None  # Choose None for any platform

# choose 3c or fp as algorithm source in the layer name below
response = access.get_wfs(platform = platform_id,
                          timewindow = (start_time, end_time),
                          layer='rsg:sorad_public_view_3c_full',
                          bbox=None)

print(response['result'][0].keys())   # show all available fields

# extract (meta)data from the response
time          = [response['result'][i]['time'] for i in range(len(response['result']))]
lat           = np.array([response['result'][i]['lat'] for i in range(len(response['result']))])
lon           = np.array([response['result'][i]['lon'] for i in range(len(response['result']))])
rel_view_az   = np.array([response['result'][i]['rel_view_az'] for i in range(len(response['result']))])
sample_uuid   = np.array([response['result'][i]['sample_uuid'] for i in range(len(response['result']))])
platform_id   = np.array([response['result'][i]['platform_id'] for i in range(len(response['result']))])
platform_uuid = np.array([response['result'][i]['platform_uuid'] for i in range(len(response['result']))])
gps_speed     = np.array([response['result'][i]['id'] for i in range(len(response['result']))])
tilt_avg      = np.array([response['result'][i]['tilt_avg'] for i in range(len(response['result']))])
tilt_std      = np.array([response['result'][i]['tilt_std'] for i in range(len(response['result']))])

# define common wavelength grid for (ir)radiance data (these are provided at their native, sensor-specific resolution
wl_out = np.arange(300, 1001, 1)
ed = access.get_l1spectra(response, 'ed_', wl_out) # # irradiance spectra in 2D matrix format: rows time index, columns wavelength
ls = access.get_l1spectra(response, 'ls_', wl_out)
lt = access.get_l1spectra(response, 'lt_', wl_out)

rrswl = np.arange(response['result'][0]['c3_wl_grid'][0], response['result'][0]['c3_wl_grid'][1], response['result'][0]['c3_wl_grid'][2])  # reconstruct wavelength grid for Rrs
rrs = np.array([response['result'][i]['c3_rrs'][:] for i in range(len(response['result']))]) # 2D matrix format: rows time index, columns wavelength
```

## Hyperspectral pyranometer (HSP)

The HSP1 measures the spectrum of downwelling solar radiation, and how this is partitioned between Direct, Diffuse and Global Irradiance. 
This sensor provides a reference for the colour or spectral distribution of sunlight near the water's surface.

HSP1 can be used to improve the accuracy of other sensors making direct measurements of the reflected light from the water body.

#### Added Value of the HSP

* Measurements can be used to improve satellite water quality products by removing or correcting for atmospheric effects, especially in coastal or inland waters.
* Improving data products from other surface-based instruments (such as other MONOCLE instruments) where a reference is normally too expensive
* In situations where existing equipment is difficult or impossible to use such as moving platforms, boats, buoys or aircraft.

#### Functionality of the submodule

The access function downloads paged data from the HSP Geoserver layer hosted at PML. Queries can be refined based on time, sensor_id or bounding box.
Extend the access function to include any additional CQL filters needed. 

The HSP geoserver layer offers unfiltered, calibrated total and diffused irradiance spectra. 

#### Minimum code example

```
from monda.hsp import access
import datetime

sensor_id = None
start_time = datetime.datetime(2021,10,21,12,0,0)
end_time   = datetime.datetime(2021,10,21,14,0,0)
bbox = None
layer = 'rsg:hsp_public_view_full'

response = access.get_wfs(sensor = sensor_id,
                          timewindow = (start_time, end_time),
                          layer=layer, bbox=bbox)

# show first record
for key, val in response['result'][0].items():
    print(f"{key}: {val}")
```


## Mini-Secchi disk
Secchi disk measurements have been used to record water transparency over hundreds of years, and are used as a proxy for eutrophication 
in international water quality monitoring programmes. Recently, a simple hand-held device was designed to measure the Secchi depth 
and water colour (Forel-Ule scale) of lake, estuarine and nearshore regions. The device additionally comes shipped with a clip for pH paper. 

The mini-secchi App, available from Google play and Apple stores, aids with data collection through smartphones. The app will:
- provide instructions for safe and correct use
- gather all measurement data, including (optional) photos and quality control questions
- upload measurement data to the PML server
- in future, the app will show the user the results from nearby and recent observations. 

Once data are received at PML they are immediately processed and made available in WMS/WFS format. The code example below illustrates how to query the WFS. 

The Mini-Secchi disk was developed as a school project in collaboration with PML. 
Due to increased demand for larger quantities of Mini-Secchi disks in citizen science projects, the developers founded Brewtek to produce the devices at scale. 


#### Added value
The device is manufactured with marine resistant materials (mostly biodegradable) using a 3D printer and basic workshop tools. 
It is inexpensive to manufacture, lightweight, easy to use, and accessible to a wide range of users. 

It builds on a long tradition in optical limnology and oceanography, but is modified for ease of operation in smaller water bodies, and from small watercraft and platforms.

The 3D-printable design is particularly useful in community, educational and hobbyist settings.
For more information see the Mini-Secchi disk pages on the [project website](https://monocle-h2020.eu/Sensors_and_services/Mini-secchi_disk)


#### Functionality of the Mini-Secchi disk and app

Note: Mini-Secchi disk data on the geoserver are still a work in progress. Data formatting is likely to change, and test records have not yet been removed. 

#### Minimum code example

```
from monda.minisecchi import access
import datetime

start_date = datetime.datetime(2022,1,1,0,0,0)
end_date = datetime.datetime(2022,8,1,0,0,0)

geoserver_layer = 'rsg:minisecchi_public_view'
bbox = None  # alternatively supply a tuple of (corner1 lat, corner1 lon, corner2 lat, corner2 lon)

response = access.get_wfs(count=1000, timewindow=(start_date, end_date), layer=geoserver_layer, bbox=bbox)

data = response['result']
print(data[0].keys())
```
