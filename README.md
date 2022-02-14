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
- Python (>= 3.6)
- NumPy (>= 1.13.3)
- scikit-learn(>=0.23.2)
- Matplotlib (>=3.3.3)
- requests (>=2.27.1)

# Source code
To get the most up to date version of the source code please see the repository at:
```
https://github.com/monocle-h2020/MONDA.git
```

## Installation
You can either clone the repo here and then run pip install from that directory or you can install directly from pypi 
using:

```python3 -m pip install monda```

For now while we are in test mode it would be just pip installed from the local code directory:
```pip install .```

## Citation
If you use MONDA in a scientific publication, we would appreciate citations: 
*add in how we want this citing*

##Contributors
This code was developed with input from Plymouth Marine laboratory (thja-pml@github, tjor@github, StefaSimis@github) and 
Water Insight (Semhar-Ghe@github, waterthing@github). 

## Submodule Information
The package contains access, quality control and visualisation tools for a number of sensor systems, for which details are provided below.

### WISP (station)
The WISPstation is a fixed position optical instrument used for measuring water-leaving reflectance.
It records radiance and irradiance with an extended wavelength range of 350nm to 1100nm in two viewing directions,
which enables continuous and autonomous high-quality measurements for water quality monitoring and satellite validation. 
The reflectance observations are used to validate satellite measurements of water-leaving reflectance. 
Concentrations of the most important bio-physical water quality parameters such as chlorophyll-a, cyanobacterial pigment, turbidity and suspended matter, are derived from the reflectance measurement. The WISPstation sends the measurements automatically over 3G/4G/5G to the “WISPcloud” cloud database which makes the results available via an API. Measurement frequency is by default a 15 min interval but be adjusted to suit user requirements. 

### About WISPcloud
WISPcloud is a scalable Postgres database that autonomously receives, stores, performs quality control and 
applies water quality algorithms to all WISPstation measurements. It has an advanced API to serve data requests directly to customers. A separate online documentation can be found here. 

### Acknowledgement 
The WISPstation public data were collected by users participating on H2020 funded projects such as EOMORES(http://eomores-h2020.eu), TAPAS(http://tapas-h2020.eu/) and MONOCLE(https://monocle-h2020.eu/). 

### Example data availability
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
 
### About the api_public access example script
An example script is provided to connect with the WISPcloud API and subsequently plot Rrs and (ir)radiance measurements using date and instrument serial number as input arguments. 


### So-Rad
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

#### Added Value of So-Rad ####
* Off-shore satellite validation is currently limited to research vessels and fixed moorings that are costly to 
maintain. The So-Rad can be installed on non-stationary platforms and is ideally suited to be included on merchant 
vessels. Ferry routes are recommended because of predictable routes and schedules. Periodic sensor maintenance can be 
easily carried out by non-expert crew.
* A high degree of automation and low-power components means the platform can be installed in remote locations for 
autonomous operation.
