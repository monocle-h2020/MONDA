# MONDA (MONocle DAta package)

## Package Description
This package contains a suite of tools for retrieving, apply quality checks, analysing and plotting data from the 
sensors that are included in the MONOCLE observation network.  The tools have been developed within the work of the 
MONOCLE project.

The MONOCLE has project created a framework for building integrated water quality sensor networks to enhance the 
utility and accessiblity of data from multiple sources, giving a more complete data landscape.

For more information on the MONOCLE project see:
[MONOCLE](https://www.monocle-h2020.eu/Home)

## Dependencies
This code requires:
- Python (>= 3.6)
- NumPy (>= 1.13.3)
- scikit-learn(>=0.23.2)
- Matplotlib (>=3.3.3)
- requests (>=2.27.1)

## Installation
pip install monda

# Source code
To get the most up to date version of the code then checkout the repo with:
```
https://github.com/monocle-h2020/MONDA.git
```

## Citation
If you use MONDA in a scientific publication, we would appreciate citations: 
*add in how we want this citing*

## Submodule Information
The package contains access, qc and analysis tools for different types of sensors that are all used for water quality 
measurement.  Details of these are provided below.


### WISP (station)
The WISPstation is a fixed position optical instrument used for measuring water-leaving reflectance.
It records radiance and irradiance with an extended wavelength range of 350nm to 1100nm in two carefully chosen viewing directions,
which enables continuous and autonomous high-quality measurements for water quality monitoring and satellite validation. 
The reflectance observations are used to validate satellite measurements of water-leaving reflectance. 
In addition to the reflectance it determines the concentrations of the most important bio-physical water quality parameters such as chlorophyll, cyanobacterial pigment, turbidity and suspended matter.
The WISPstation sends the measurements automatically over 3G/4G/5G to our cloud database “WISPcloud”, which makes the results available via an API. Measurement frequency is by default is 15min interval ,but it can easily be adjusted to suit the user’s requirements. 
### About WISPcloud
WISPcloud is a scalable Postgres database that autonomously receives, stores, performs quality control and 
applies water quality algorithms to all WISPstation measurements. It has an advanced API to serve data requests directly to customers. A separate online documentation can be found here. 

### Acknowledgement 
The WISPstation public data were collected by users participating on H2020 funded projects such as EOMORES(http://eomores-h2020.eu), TAPAS(http://tapas-h2020.eu/) and MONOCLE(https://monocle-h2020.eu/). 

### Table that shows the overview of available dataset
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
We made an example to demonstrate how to connect with the WISPcloud api and plot Rrs and ir/radiance measurements using date and instrument serial as an input. 

#### username & password to access WISPcloud public API 
Username : public_access
Password:  WISPstation


### So-Rad

The So-Rad is a low-power, low cost autonomous system to obtain high-frequency water-leaving reflectance from 
non-stationary platforms such as ships and buoys. So-Rad software is highly configurable and open-source. 
So-Rad optimizes the measurement geometry of commercially available sensors which increases the number of successful 
observations of water colour obtained from moving platforms (Simis and Olsson 2013). 

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
