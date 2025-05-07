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

def HRRRdownloader(xtime,HRRRpath):
    # Input: timestamp (string or datime.datetime), path_to_HRRRpath (string)
    # timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object
    # Function: Downloads and saves HRRRv4 meteorological data of given timestamp to desired path
    # Returns: None

    if type(xtime) == dt.datetime:
        xtime = xtime.strftime('%Y-%m-%dt%H')
    elif type(xtime) == str: 
        if not xtime[4] == '-':
            raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
        if not xtime[7] == '-':
            raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
        try:
            int(xtime[:4])
            int(xtime[5:7])
            int(xtime[8:10])
            int(xtime[-2:])
        except:
            raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
    else:
        raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
    
    if not os.path.isdir(HRRRpath):
        raise Exception('Directory to save meteorological data not valid.')

    H = Herbie(xtime[:10] + ' '+xtime[-2:]+':00:00', model='hrrr', product='nat', save_dir=HRRRpath,verbose=True, priority='google')
    path = str(H.download())
    filename = 'hrrr.' + xtime[:10].replace('-','') + '.t'+xtime[-2:]+'z.wrfnatf00.grib2'
    os.system('mv '+path+' '+path[:-25]+filename)
    os.system('mv '+path[:-25]+filename+' '+HRRRpath)
    return True



def get_meteo_all(lat,lon,xtime,HRRRpath):
    # Input: latitude (float or string of numbers), longitude (float or string of numbers), timestamp (string or datetime.datetime), path_to_HRRRpath (string)
    # timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object
    # Returns: PBLH, friction, convective, LMO, coriolis, ws, wd, z0 (all float)

    try:
        lat,lon = float(lat),float(lon)
    except:
        raise Exception('Invalid latitude and longitude form. They must be float or string of numbers.')

    if type(xtime) == dt.datetime:
        xtime = xtime.strftime('%Y-%m-%dt%H')
    elif type(xtime) == str: 
        if not xtime[4] == '-':
            raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
        if not xtime[7] == '-':
            raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
        try:
            int(xtime[:4])
            int(xtime[5:7])
            int(xtime[8:10])
            int(xtime[-2:])
        except:
            raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')
    else:
        raise Exception('Timestamp format error: xtime should be a datetime.datetime object or string of the form yyyy-mm-ddthh.')

    if not os.path.isdir(HRRRpath):
        raise Exception('Directory to save meteorological data not valid.')


    # hrrr file name example: hrrr.20230214_conus_hrrr.t00z.wrfnatf00.grib2
    filename = 'hrrr.' + xtime[:10].replace('-','') + '.t'+xtime[-2:]+'z.wrfnatf00.grib2'
    # read file: if it does not exist, download
    listdir = os.listdir(HRRRpath)
    if filename in listdir:
        g = grib2io.open(HRRRpath+filename,'r')
    else:
        timestamp = xtime
        print('HRRR data not found, attempting to download...')
        HRRRdownloader(timestamp,HRRRpath)
        g = grib2io.open(HRRRpath+filename,'r')

    lats, lons = g[0].latlons() 

    def dist(A,B):
        return ((A[0]-B[0])**2 + (A[1]-B[1])**2) ** 0.5
    
    x,y,finX,finY = 0,0,0,0
    L = 1e10
    while x < len(lats):
        while y < len(lats[0]):
            d = dist( (lats[x][y],lons[x][y]), (lat,lon) )
            if d<L:
                L = d
                finX,finY = x,y
            else:
                pass
            y+=1
        y = 0
        x+=1

    
    grav = 9.81 # m/s/s, gravitational acceleration
    R = 287.0 # specific gas constant
    kappa = 0.40 # Von Karman constant
    Cp = 1003.5 # specific heat of air at constant pressure

    pressure = g[1031].data[finX,finY] # read pressure
    humidity_spec = g[1042].data[finX,finY] # read specific humidity
    Temp = g[1033].data[finX,finY] # read temp
    PBLH = g[1112].data[finX,finY] # read PHL height
    friction = g[1064].data[finX,finY] # read friction velocity
    S = g[1067].data[finX,finY] # read ground heat flux
    z0 = g[1063].data[finX,finY]
    coriolis = 1e-4
    #rho = 1 # read density

    def cube(x):
        if 0<=x: 
            return x**(1./3.)
        else:
            return -(-x)**(1./3.)
        
    bad_output = False    
    
    rho = (pressure * 0.0289652) / (8.31446261815324 * Temp)
    q = S / (Cp * rho) # derive sensible heat flux (has units of K*m/s)
    T_virt = Temp * (0.608*humidity_spec + 1) # virtual temperature
    theta_virt = T_virt * (pressure/100000.0)**( -2./7.)# virtual potential temperature
    Zm = PBLH
    #print(Temp,S,q,humidity_spec,pressure,T_virt,theta_virt)
    convective = cube((q*grav*Zm) / (theta_virt))
    LMO = -Zm * ( ( (friction) / (convective) ) ** 3 )

    if(convective != 0):
        LMO = -Zm * ( ( (friction) / (convective) ) ** 3 )
    else:
        LMO = 0
        print('Warning: Ground heat flux from HRRR is zero. LMO will be set to 0.')
        bad_output = True
    u = g[1046].data[finX,finY]
    v = g[1047].data[finX,finY]

    ws = float((u**2 + v**2) ** 0.5)

    wd = float( wind_direction(u*(units('m/s')),v*(units('m/s'))).magnitude )


    return PBLH, friction, convective, LMO, coriolis, ws, wd, z0, bad_output


def get_meteo_wind(lat,lon,xtime,HRRRpath):
    # Input: latitude (float or string of numbers), longitude (float or string of numbers), timestamp (string or datetime.datetime), path_to_HRRRpath (string)
    # timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object
    # Returns: ws, wd (both float)

    PBLH, friction, convective, LMO, coriolis, ws, wd, z0, bad_output = get_meteo_all(lat,lon,xtime,HRRRpath)
    return ws, wd

def get_meteo_obukhov(lat,lon,xtime,HRRRpath):
    # Input: latitude (float or string of numbers), longitude (float or string of numbers), timestamp (string or datetime.datetime), path_to_HRRRpath (string)
    # timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object
    # Returns: LMO (Obukhov length) (float)

    PBLH, friction, convective, LMO, coriolis, ws, wd, z0, bad_output = get_meteo_all(lat,lon,xtime,HRRRpath)
    return LMO

def get_meteo_roughness(lat,lon,xtime,HRRRpath):
    # Input: latitude (float or string of numbers), longitude (float or string of numbers), timestamp (string or datetime.datetime), path_to_HRRRpath (string)
    # timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object
    # Returns: z0 (surface roughness length) (float)

    PBLH, friction, convective, LMO, coriolis, ws, wd, z0, bad_output = get_meteo_all(lat,lon,xtime,HRRRpath)
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

def get_stability(lat,lon,xtime,HRRRpath):
    # Input: latitude (float or string of numbers), longitude (float or string of numbers), timestamp (string or datetime.datetime), path_to_HRRRpath (string)
    # timestamp should be of the form yyyy-mm-ddthh OR datetime.datetime object
    # Returns: list of length 4 
    # Returned list is of the form [class1 (str), weight1 (float), class2 (str), weight2 (float)] (weights will be 0.<=x<=1.)
    # If no blending, then weight1 = 1.0 and class2 = 'FALSE

    pbl, friction, convective, LMO, coriolis, ws, wd, z0, bad_output = get_meteo_all(lat,lon,xtime,HRRRpath)
    if bad_output:
        results = ['D','1.0','FALSE','0.0']
        print('Warning: Ground heat flux from HRRR is zero. Returning the default stability class D.')
    else:
        results = calculate_weights(z0,LMO)
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

    pbl, friction, convective, LMO, coriolis, ws, wd, z0, bad_output = get_meteo_all(lati,longi,testtime,'./')

    print('\t\tPBL, frictional velocity, convective velocity, LMO, coriolis, z0: ', pbl, friction, convective, LMO, coriolis, z0)
    print('\t\tWind speed, direction: ', ws, wd)

    results = calculate_weights(z0,LMO)

    print('\t\t',results)
