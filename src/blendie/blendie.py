# Created by Shoma Yamanouchi (shoma.yamanouchi@ec.gc.ca)

######## Requirements ########
# Python 3+ (3.9+ recommended) 
# numpy
# datetime
# metpy
# herbie
# grib2io
##############################

# IMPORTANT NOTICE: This package will ONLY work in North America #



import sys, os, time
import numpy as np
import datetime as dt
import grib2io
from herbie import Herbie
from metpy.calc import wind_direction
from metpy.units import units


# =====================================================================
# MODEL DOWNLOADERS
# =====================================================================

def HRRRdownloader(xtime, HRRRpath):
    # Input: timestamp (string or datetime.datetime), path_to_HRRRpath (string)
    #
    # timestamp should be of the form yyyy-mm-ddthh
    # OR datetime.datetime object
    #
    # Function: Downloads and saves HRRRv4 meteorological data of
    # given timestamp to desired path
    #
    # Returns: True

    if type(xtime) == dt.datetime:

        xtime = xtime.strftime('%Y-%m-%dt%H')

    elif type(xtime) == str:

        if not xtime[4] == '-':
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

        if not xtime[7] == '-':
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

        try:
            int(xtime[:4])
            int(xtime[5:7])
            int(xtime[8:10])
            int(xtime[-2:])
        except:
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

    else:

        raise Exception(
            'Timestamp format error: xtime should be a '
            'datetime.datetime object or string of the form yyyy-mm-ddthh.'
        )

    if not os.path.isdir(HRRRpath):
        raise Exception(
            'Directory to save meteorological data not valid.'
        )

    H = Herbie(
        xtime[:10] + ' ' + xtime[-2:] + ':00:00',
        model='hrrr',
        product='nat',
        save_dir=HRRRpath,
        verbose=True,
        priority='google'
    )

    path = str(H.download())

    filename = (
        'hrrr.' +
        xtime[:10].replace('-', '') +
        '.t' +
        xtime[-2:] +
        'z.wrfnatf00.grib2'
    )

    # Preserve the original BLENDIE filename convention
    target = os.path.join(HRRRpath, filename)

    if os.path.abspath(path) != os.path.abspath(target):

        try:
            os.replace(path, target)
        except OSError:
            # Fall back to the original behavior if necessary
            os.system('mv "' + path + '" "' + target + '"')

    return True


def RAPdownloader(xtime, RAPpath):
    # Input: timestamp (string or datetime.datetime), path_to_RAPpath (string)
    #
    # timestamp should be of the form yyyy-mm-ddthh
    # OR datetime.datetime object
    #
    # Function: Downloads RAP 13-km native-grid meteorological data
    #()
    # Returns: True

    if type(xtime) == dt.datetime:

        xtime = xtime.strftime('%Y-%m-%dt%H')

    elif type(xtime) == str:

        if not xtime[4] == '-':
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

        if not xtime[7] == '-':
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

        try:
            int(xtime[:4])
            int(xtime[5:7])
            int(xtime[8:10])
            int(xtime[-2:])
        except:
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

    else:

        raise Exception(
            'Timestamp format error: xtime should be a '
            'datetime.datetime object or string of the form yyyy-mm-ddthh.'
        )

    if not os.path.isdir(RAPpath):
        raise Exception(
            'Directory to save meteorological data not valid.'
        )

    H = Herbie(
        xtime[:10] + ' ' + xtime[-2:] + ':00:00',
        model='rap',
        product='wrfnat',
        fxx=0,
        save_dir=RAPpath,
        verbose=True
    )

    path = str(H.download())

    filename = (
        'rap.' +
        xtime[:10].replace('-', '') +
        '.t' +
        xtime[-2:] +
        'z.wrfnatf00.grib2'
    )

    target = os.path.join(RAPpath, filename)

    # Herbie normally returns the correctly named RAP file.
    # Rename only if necessary.
    if os.path.abspath(path) != os.path.abspath(target):

        try:
            os.replace(path, target)
        except OSError:
            os.system('mv "' + path + '" "' + target + '"')

    return True


def haversine_distance(coord1, coord2):
    """
    Calculates the great-circle distance between two points 
    on the Earth's surface using NumPy.
    
    Parameters:
    coord1, coord2 : tuples or arrays of (latitude, longitude) in degrees.
    
    Returns:
    distance : float (Distance in kilometers)
    """
    # Earth's mean radius in kilometers
    EARTH_RADIUS = 6371.0
    
    # Extract latitudes and longitudes
    lat1, lon1 = coord1
    lat2, lon2 = coord2
    
    # Convert degrees to radians
    lat1, lon1, lat2, lon2 = np.radians([lat1, lon1, lat2, lon2])
    
    # Haversine formula core steps
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    # Calculate final distance
    distance_km = EARTH_RADIUS * c
    return distance_km


# =====================================================================
# METEOROLOGICAL CALCULATION
# =====================================================================

def get_meteo_all(
        lat,
        lon,
        xtime,
        HRRRpath,
        model='hrrr'):
    """
    Get meteorological parameters and calculate Obukhov length.

    Parameters
    ----------
    lat : float or str
        Latitude.

    lon : float or str
        Longitude.

    xtime : str or datetime.datetime
        Timestamp in the form yyyy-mm-ddthh.

    HRRRpath : str
        Directory containing meteorological GRIB2 files.

        For model='hrrr', this is the HRRR directory.
        For model='rap', this is the RAP directory.

        The parameter name HRRRpath is deliberately retained for
        backwards compatibility with BLENDIE v1.0.0.

    model : str, optional
        Meteorological model to use.

        'hrrr' (default)
            NOAA HRRR 3-km data.

        'rap'
            NOAA RAP 13-km data.

    Returns
    -------
    PBLH, friction, convective, LMO, coriolis, ws, wd, z0

    All returned quantities have the same meaning and ordering as
    BLENDIE v1.0.0.
    """

    # -------------------------------------------------------------
    # Validate model
    # -------------------------------------------------------------

    model = str(model).lower()

    if model not in ['hrrr', 'rap']:
        raise Exception(
            "Invalid model. model must be either 'hrrr' or 'rap'."
        )

    # -------------------------------------------------------------
    # Latitude / longitude
    # -------------------------------------------------------------

    try:
        lat, lon = float(lat), float(lon)

    except:

        raise Exception(
            'Invalid latitude and longitude form. '
            'They must be float or string of numbers.'
        )

    # -------------------------------------------------------------
    # Timestamp
    # -------------------------------------------------------------

    if type(xtime) == dt.datetime:

        xtime = xtime.strftime('%Y-%m-%dt%H')

    elif type(xtime) == str:

        if not xtime[4] == '-':
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

        if not xtime[7] == '-':
            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

        try:

            int(xtime[:4])
            int(xtime[5:7])
            int(xtime[8:10])
            int(xtime[-2:])

        except:

            raise Exception(
                'Timestamp format error: xtime should be a '
                'datetime.datetime object or string of the form yyyy-mm-ddthh.'
            )

    else:

        raise Exception(
            'Timestamp format error: xtime should be a '
            'datetime.datetime object or string of the form yyyy-mm-ddthh.'
        )

    # -------------------------------------------------------------
    # Validate data directory
    # -------------------------------------------------------------

    if not os.path.isdir(HRRRpath):

        raise Exception(
            'Directory to save meteorological data not valid.'
        )

    # -------------------------------------------------------------
    # Model-specific filename and GRIB record numbers
    # -------------------------------------------------------------

    if model == 'hrrr':

        filename = (
            'hrrr.' +
            xtime[:10].replace('-', '') +
            '.t' +
            xtime[-2:] +
            'z.wrfnatf00.grib2'
        )

        # HRRR v4 wrfnat analysis
        #
        # These are the same records used by BLENDIE v1.0.0,
        # except that sensible heat flux is corrected from
        # record 1067 to record 1065.

        GRIB = {
            'pressure': 1031,
            'humidity': 1042,
            'temperature': 1033,
            'pblh': 1112,
            'friction': 1064,
            'sensible_heat_flux': 1065,   #  Improved from v1.0.0, which used ground heat flux instead of sensible net heat flux.
            'ground_heat_flux': 1067, # Use the above sensible heat flux (1065) instead.
            'roughness': 1063,
            'u': 1046,
            'v': 1047
        }

    else:

        print('Using RAP model as specified')

        filename = (
            'rap.' +
            xtime[:10].replace('-', '') +
            '.t' +
            xtime[-2:] +
            'z.wrfnatf00.grib2'
        )

        # RAP 13-km wrfnat analysis
        #
        # NOAA RAP inventory:
        #
        # 1005 surface pressure
        # 1007 surface temperature
        # 1031 2-m specific humidity
        # 1034 10-m U wind
        # 1035 10-m V wind
        # 1049 surface roughness
        # 1050 friction velocity
        # 1051 sensible heat flux
        # 1053 ground heat flux
        # 1085 PBL height

        GRIB = {
            'pressure': 1004,
            'humidity': 1030,
            'temperature': 1006,
            'pblh': 1084,
            'friction': 1049,
            'sensible_heat_flux': 1050,
            'ground_heat_flux': 1052,
            'roughness': 1048,
            'u': 1033,
            'v': 1034
        }

    filepath = os.path.join(
        HRRRpath,
        filename
    )

    # -------------------------------------------------------------
    # Read file, downloading if necessary
    # -------------------------------------------------------------

    if os.path.isfile(filepath):

        g = grib2io.open(
            filepath,
            'r'
        )

    else:

        if model == 'hrrr':

            print(
                'HRRR data not found, attempting to download...'
            )

            HRRRdownloader(
                xtime,
                HRRRpath
            )

        elif model == 'rap':

            print(
                'RAP data not found, attempting to download...'
            )

            RAPdownloader(
                xtime,
                HRRRpath
            )

        if not os.path.isfile(filepath):

            raise Exception(
                model.upper() +
                ' data download completed, but expected file '
                'was not found: ' +
                filepath
            )

        g = grib2io.open(
            filepath,
            'r'
        )

    # -------------------------------------------------------------
    # Latitude / longitude
    # -------------------------------------------------------------

    lats, lons = g[0].latlons()

    # -------------------------------------------------------------
    # Distance function
    # -------------------------------------------------------------

    def dist(A, B):

        return (
            (A[0] - B[0]) ** 2 +
            (A[1] - B[1]) ** 2
        ) ** 0.5

    # -------------------------------------------------------------
    # Find nearest grid cells
    #
    # Keep the same nearest-neighbor methodology as v1.0.0.
    # -------------------------------------------------------------

    dist_xy = []

    for i in range(len(lats)):

        for j in range(len(lats[0])):

            dist_xy.append(
                [
                    dist(
                        (lats[i][j], lons[i][j]),
                        (lat, lon)
                    ),
                    i,
                    j
                ]
            )

    dist_xy = sorted(
        dist_xy,
        key=lambda x: x[0]
    )

    # -------------------------------------------------------------
    # Constants
    # -------------------------------------------------------------

    def cube(x):

        if 0 <= x:

            return x ** (1. / 3.)

        else:

            return -(-x) ** (1. / 3.)

    grav = 9.81
    R = 287.0
    kappa = 0.40
    Cp = 1003.5

    # -------------------------------------------------------------
    # Find nearest valid meteorological pixel
    # -------------------------------------------------------------

    nearest_valid_pixel = 0
    convective = 0

    while convective == 0:

        if nearest_valid_pixel >= len(dist_xy):

            raise Exception(
                model.upper() +
                ' data invalid... Please check the ' +
                model.upper() +
                ' data'
            )

        finX = dist_xy[nearest_valid_pixel][1]
        finY = dist_xy[nearest_valid_pixel][2]

        # ---------------------------------------------------------
        # Read meteorological quantities
        # ---------------------------------------------------------

        pressure = g[
            GRIB['pressure']
        ].data[finX, finY]

        humidity_spec = g[
            GRIB['humidity']
        ].data[finX, finY]

        Temp = g[
            GRIB['temperature']
        ].data[finX, finY]

        PBLH = g[
            GRIB['pblh']
        ].data[finX, finY]

        friction = g[
            GRIB['friction']
        ].data[finX, finY]

        # ---------------------------------------------------------
        # IMPORTANT FIX:
        #
        # Use SENSIBLE heat flux rather than ground heat flux.
        #
        # HRRR:
        #   1066 = SHTFL
        #   1068 = GFLUX
        #
        # RAP:
        #   1051 = SHTFL
        #   1053 = GFLUX
        # ---------------------------------------------------------

        S = g[
            GRIB['sensible_heat_flux']
        ].data[finX, finY]

        z0 = g[
            GRIB['roughness']
        ].data[finX, finY]

        coriolis = 1e-4

        # ---------------------------------------------------------
        # Air density
        # ---------------------------------------------------------

        rho = (
            pressure * 0.0289652
        ) / (
            8.31446261815324 * Temp
        )

        # ---------------------------------------------------------
        # Convert sensible heat flux to K*m/s
        #
        # S has units W/m^2.
        #
        # q = H / (Cp * rho)
        #
        # This has units K*m/s.
        # ---------------------------------------------------------

        q = S / (Cp * rho)

        # ---------------------------------------------------------
        # Virtual temperature
        # ---------------------------------------------------------

        T_virt = Temp * (
            0.608 * humidity_spec + 1
        )

        # ---------------------------------------------------------
        # Virtual potential temperature
        # ---------------------------------------------------------

        theta_virt = T_virt * (
            pressure / 100000.0
        ) ** (-2. / 7.)

        # ---------------------------------------------------------
        # Convective velocity scale
        # ---------------------------------------------------------

        Zm = PBLH

        convective = cube(
            (
                q *
                grav *
                Zm
            ) /
            theta_virt
        )

        nearest_valid_pixel += 1

    print(
        'Using the number ' +
        str(nearest_valid_pixel) +
        ' closest ' +
        model.upper() +
        ' pixel'
    )
    distance_delta = haversine_distance((lat,lon),(lats[dist_xy[nearest_valid_pixel][1]][dist_xy[nearest_valid_pixel][2]], lons[dist_xy[nearest_valid_pixel][1]][dist_xy[nearest_valid_pixel][2]]))
    if model == 'hrrr':
        if distance_delta > 5.0:
            print('WARNING!!! Closest model grid is more than 5km away from specified location. Use with caution! (Blendie only works in North America. Try switching the model to RAP for wider coverage over North America.)')
    else:
        if distance_delta > 20.0:
            print('WARNING!!! Closest model grid is more than 5km away from specified location. Use with caution! (Blendie only works in North America.)')
    print(
        'Distance between observation and pixel used: ',distance_delta,'km'
    )
    # -------------------------------------------------------------
    # Obukhov length
    # -------------------------------------------------------------

    LMO = -Zm * (
        (
            friction /
            convective
        ) ** 3
    )

    # -------------------------------------------------------------
    # Wind
    # -------------------------------------------------------------

    u = g[
        GRIB['u']
    ].data[finX, finY]

    v = g[
        GRIB['v']
    ].data[finX, finY]

    ws = float(
        (
            u ** 2 +
            v ** 2
        ) ** 0.5
    )

    wd = float(
        wind_direction(
            u * units('m/s'),
            v * units('m/s')
        ).magnitude
    )

    # -------------------------------------------------------------
    # Return exactly the same structure as BLENDIE v1.0.0
    # -------------------------------------------------------------

    return (
        PBLH,
        friction,
        convective,
        LMO,
        coriolis,
        ws,
        wd,
        z0
    )


# =====================================================================
# WRAPPER FUNCTIONS
# =====================================================================

def get_meteo_wind(
        lat,
        lon,
        xtime,
        HRRRpath,
        model='hrrr'):

    # Input: latitude, longitude, timestamp, data path
    #
    # model:
    #   'hrrr' = default, preserves original behavior
    #   'rap'  = use RAP 13-km data
    #
    # Returns: ws, wd

    PBLH, friction, convective, LMO, coriolis, ws, wd, z0 = \
        get_meteo_all(
            lat,
            lon,
            xtime,
            HRRRpath,
            model=model
        )

    return ws, wd


def get_meteo_obukhov(
        lat,
        lon,
        xtime,
        HRRRpath,
        model='hrrr'):

    # Input: latitude, longitude, timestamp, data path
    #
    # model:
    #   'hrrr' = default, preserves original behavior
    #   'rap'  = use RAP 13-km data
    #
    # Returns: LMO (Obukhov length)

    PBLH, friction, convective, LMO, coriolis, ws, wd, z0 = \
        get_meteo_all(
            lat,
            lon,
            xtime,
            HRRRpath,
            model=model
        )

    return LMO


def get_meteo_roughness(
        lat,
        lon,
        xtime,
        HRRRpath,
        model='hrrr'):

    # Input: latitude, longitude, timestamp, data path
    #
    # model:
    #   'hrrr' = default, preserves original behavior
    #   'rap'  = use RAP 13-km data
    #
    # Returns: z0 (surface roughness length)

    PBLH, friction, convective, LMO, coriolis, ws, wd, z0 = \
        get_meteo_all(
            lat,
            lon,
            xtime,
            HRRRpath,
            model=model
        )

    return z0


def calculate_weights(z0,L):
    # Input: z0 (float), LMO (float)
    # Input must be float, surface roughness and Obukhov length
    # Returns: list of length 4 
    # Returned list is of the form [class1 (str), weight1 (float), class2 (str), weight2 (float)] (weights will be 0.<=x<=1.)
    # If no blending, then weight1 = 1.0 and class2 = 'FALSE'

    x = 1./L # 1/L
    # center lines of the stability classes following Sucevic & Djurisic (2012)
    A = [-11.4, 0.1]#(-11.4*x)**(-1./0.1)
    B = [-26., 0.17]#(-26.*x)**(-1./0.17)
    C = [-123., 0.3]#(-123.*x)**(-1./0.3)
    D = 0
    E = [123., 0.3]#(123.*x)**(-1./0.3)
    F = [26., 0.17]#(26.*x)**(-1./0.17)

    AB = [-15.59564831,   0.17468621]  # <- tweaked value. og value = [-15.25470547,   0.1664827]
    BC = [-75.,   0.35]                # <- tweaked value. og value = [-77.60063996,   0.42104886]
    CD = [-419.61923688,   0.46813927] # <- tweaked value. og value = [-1149.41858,  0.892895695]
    DE = [419.61923688,   0.46813927]  # <- tweaked value. og value = [419.61923688,   0.46813927]
    EF = [75.,  0.35]                  # <- tweaked value. og value = [80.57606239,  0.44978964]

    # Find 1/L given roughness and a, b values:
    def z0tox(a,b,z):
        return 1./(a*(z**b))

    zA  = z0tox(A[0],A[1],z0)
    zAB = z0tox(AB[0],AB[1],z0)
    zB  = z0tox(B[0],B[1],z0)
    zBC = z0tox(BC[0],BC[1],z0)
    zC  = z0tox(C[0],C[1],z0)
    zCD = z0tox(CD[0],CD[1],z0)
    zD  = 0
    zDE = z0tox(DE[0],DE[1],z0)
    zE  = z0tox(E[0],E[1],z0)
    zEF = z0tox(EF[0],EF[1],z0)
    zF  = z0tox(F[0],F[1],z0)

    #print(x,z0)
    # No blending if z0 < 0.05 OR z0 > 0.5
    # Blend if between those values
    if z0 < 0.002 or z0 > 0.5:
        if x <= zA:
            return ['A','1.0','FALSE','0.0']
        if zA < x < zAB:
            return ['A','1.0','FALSE','0.0']
        elif zAB <= x < zB:
            return ['B','1.0','FALSE','0.0']
        elif zB <= x < zBC:
            return ['B','1.0','FALSE','0.0']
        elif zBC <= x < zC:
            return ['C','1.0','FALSE','0.0']
        elif zC <= x < zCD:
            return ['C','1.0','FALSE','0.0']
        elif zCD <= x < zD:
            return ['D','1.0','FALSE','0.0']
        elif zD <= x < zDE:
            return ['D','1.0','FALSE','0.0']
        elif zDE <= x < zE:
            return ['E','1.0','FALSE','0.0']
        elif zE <= x < zEF:
            return ['E','1.0','FALSE','0.0']
        elif zEF <= x < zF:
            return ['F','1.0','FALSE','0.0']
    else:
        # no blending if x lies left of A or right of F
        if x <= zA:
            return ['A','1.0','FALSE','0.0']
        elif x >= zF:
            return ['F','1.0','FALSE','0.0']
        else:
            # Blend

            def weighingfunction(d1,std1,d2,std2):
                # returns averaging weights given distances and std deviations
                # assumes d1, std2 is for the point within the class, d2 and std2 are the closest class outside
                weight1 = 1./(d1/std1)
                weight2 = 1./(d2/std2)
                normalized_w1 = weight1/(weight1+weight2)
                normalized_w2 = weight2/(weight1+weight2)
                return normalized_w1, normalized_w2

            if zA < x < zAB:
                d1 = np.abs(x-zA)
                d2 = np.abs(x-zB)
                std1 = np.abs(zA-zAB)
                std2 = np.abs(zAB-zB)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['A',w1,'B',w2]
            
            elif zAB <= x < zB:
                d1 = np.abs(x-zB)
                d2 = np.abs(x-zA)
                std1 = np.abs(zB-zAB)
                std2 = np.abs(zAB-zA)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['B',w1,'A',w2]
            
            elif zB <= x < zBC:
                d1 = np.abs(x-zB)
                d2 = np.abs(x-zC)
                std1 = np.abs(zB-zBC)
                std2 = np.abs(zBC-zC)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['B',w1,'C',w2]  

            elif zBC <= x < zC:
                d1 = np.abs(x-zC)
                d2 = np.abs(x-zB)
                std1 = np.abs(zC-zBC)
                std2 = np.abs(zBC-zB)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['C',w1,'B',w2]  
            
            elif zC <= x < zCD:
                d1 = np.abs(x-zC)
                d2 = np.abs(x-zD)
                std1 = np.abs(zC-zCD)
                std2 = np.abs(zCD-zD)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['C',w1,'D',w2]  
                
            elif zCD <= x < zD:
                d1 = np.abs(x-zD)
                d2 = np.abs(x-zC)
                std1 = np.abs(zD-zCD)
                std2 = np.abs(zCD-zC)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['D',w1,'C',w2]  
            
            elif zD <= x < zDE:
                d1 = np.abs(x-zD)
                d2 = np.abs(x-zE)
                std1 = np.abs(zD-zDE)
                std2 = np.abs(zDE-zE)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['D',w1,'E',w2] 
            
            elif zDE <= x < zE:
                d1 = np.abs(x-zE)
                d2 = np.abs(x-zD)
                std1 = np.abs(zE-zDE)
                std2 = np.abs(zDE-zD)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['E',w1,'D',w2] 
            
            elif zE <= x < zEF:
                d1 = np.abs(x-zE)
                d2 = np.abs(x-zF)
                std1 = np.abs(zE-zEF)
                std2 = np.abs(zEF-zF)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['E',w1,'F',w2] 
            
            elif zEF <= x < zF:
                d1 = np.abs(x-zF)
                d2 = np.abs(x-zE)
                std1 = np.abs(zF-zEF)
                std2 = np.abs(zEF-zE)
                w1,w2 = weighingfunction(d1,std1,d2,std2)
                return ['F',w1,'E',w2] 


def get_stability(
        lat,
        lon,
        xtime,
        HRRRpath,
        model='hrrr'):

    # Input: latitude, longitude, timestamp, data path
    #
    # model:
    #   'hrrr' = default
    #   'rap'  = use RAP 13-km data
    #
    # Returns: list of length 4

    pbl, friction, convective, LMO, coriolis, ws, wd, z0 = \
        get_meteo_all(
            lat,
            lon,
            xtime,
            HRRRpath,
            model=model
        )

    results = calculate_weights(z0, LMO)

    return results


if __name__ == "__main__":
    ####################################
    # Run this file as main.py to test #
    ####################################

    print('Testing...')
    lati = '43.781092539863025'
    longi = '-79.4682135538496'
    testtime = '2024-11-11t10'
    print('\tTest conducted for: '+str(lati)+' (lat), '+str(longi)+' (lon), on: '+testtime)

    pbl, friction, convective, LMO, coriolis, ws, wd, z0 = get_meteo_all(lati,longi,testtime,'./',model='hrrr')

    print('\t\tPBL, frictional velocity, convective velocity, LMO, coriolis, z0: ', pbl, friction, convective, LMO, coriolis, z0)
    print('\t\tWind speed, direction: ', ws, wd)

    results = calculate_weights(z0,LMO)

    print('\t\t',results)
