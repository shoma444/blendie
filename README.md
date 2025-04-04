# BLENDIE
Python script for blending stability classes used in Gaussian Dispersion Models.

Created by Shoma Yamanouchi (shoma.yamanouchi@ec.gc.ca)

Requirements:
Python 3.9+,  numpy, datetime, metpy, herbie, grib2io

Also requires a Linux environment to run

IMPORTANT NOTICE: This package will ONLY work in North America

#
#

Package Overview:
This Python package uses HRRRv4 (https://rapidrefresh.noaa.gov/hrrr/) meteorological data to obtain the Obukhov length and surface roughness to estimate the stability class of the atmosphere, to be used in Gaussian dispersion models. 

Other meteorological parameters (e.g., wind speed and direction) can also be obtained.

#
#

INSTALLATION:

NOTE: The grib2io package is a required dependency, and it requires an external NCEPLIBS-g2c library to work (https://pypi.org/project/grib2io/). For easiest access we suggest installing these via the Anaconda ecosystem:

conda install -c conda-forge nceplibs-g2c

Installing blendie can simply be done via the pip install command:

pip install ./dist/blendie-0.0.2.tar.gz

while in the same directory as the /dist folder. If errors regarding missing libraries (e.g., libg2c.so) are raised, try reinstallilng nceplibs-g2c, as that is likely the culprit.

#
#

USAGE:

The main function of this package is: get_stability(latitude, longitude, timestamp, path_to_save_HRRR_data).

The inputs of this function are: 

latitude (float or string of numbers), 

longitude (float or string of numbers), 

timestamp (string or datetime.datetime), and 

path_to_save_HRRR_data (string).



The function will return a list of length 4, of the form: 

[class1 (str), weight1 (float), class2 (str), weight2 (float)] (weights will be 0.0 <= x <= 1.0). 

The timestamp should be of the form yyyy-mm-ddthh (note that the letter "t" separates the day and hour, e.g., 2021-12-25t19) OR a datetime.datetime object. If the result is that no blending should be done, then weight1 = 1.0 and class2 = 'FALSE'.

#

Run the /src/blendie/blendie.py file to test the package. The results should look something like:

Testing...

  Test conducted for: 43.781092539863025 (lat), -79.4682135538496 (lon), on: 2024-11-11t10
  
HRRR data not found, attempting to download...
✅ Found ┊ model=hrrr ┊ product=nat ┊ 2024-Nov-11 10:00 UTC F00 ┊ GRIB2 @ google ┊ IDX @ google
👨🏻‍🏭 Created directory: [hrrr/20241111]
                
   PBL, frictional velocity, convective velocity, LMO, coriolis, z0:  555.89825 0.42 1.0329707536125516 -37.366211195124826 0.0001 0.5
   
   Wind speed, direction:  5.218844608196581 236.46394794175023

   ['B', 0.7924253315186035, 'C', 0.2075746684813964]

