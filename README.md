# blender
Python scripts for blending stability classes used in Gaussian Dispersion Models.

Created by Shoma Yamanouchi (shoma.yamanouchi@ec.gc.ca)

Requirements:
Python 3+ (3.9+ recommended), numpy, datetime, metpy, herbie, grib2io

IMPORTANT NOTICE: This package will ONLY work in North America


Package Overview:
This Python package uses HRRRv4 (https://rapidrefresh.noaa.gov/hrrr/) meteorological data to obtain the Obukhov length and surface roughness to estimate the stability class of the atmosphere, to be used in Gaussian dispersion models. 

Other meteorological parameters (e.g., wind speed and direction) can also be obtained.

The main function of this package is: get_stability(lat,lon,xtime,HRRRpath).

The inputs of this function are: latitude (float or string of numbers), longitude (float or string of numbers), timestamp (string or datetime.datetime), path_to_HRRRpath (string), and it will return a list of length 4, of the form [class1 (str), weight1 (float), class2 (str), weight2 (float)] (weights will be 0.<=x<=1.). The timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object. If no blending, then weight1 = 1.0 and class2 = 'FALSE'.

Run the main.py to test.
