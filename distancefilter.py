# point_distance_filter

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:05:17 2017

@author: lpalao
removes points location based on specified distance
This is the recent version of filtering points location based on a specified distance. The computation is based on haversine and vincentry
What it does: an algorithm that removes points that are within a specified distance to each other
How to use:
    Run it in command prompt or in any IDLE. Preferrably install python using Anaconda
    Make sure to have the following packages installed (easygui, shapefile, pandas)
    In line 183 and 184, choose the preferred computation for distance, haversine or vicentry. Put a comment (#) on
    the algorithm that WILL NOT be used
    Specify the shapefile in the file browser (*.shp) in GCS:WGS84 coordinate reference system
    Specify the Distance filter value in kilometers
    Input the name of the crop

In the same directory, a csv file will be created that is ready to be used in Maxent
"""

import math, easygui, shapefile, itertools, os
import pandas as pd
import numpy as np
from time import time

filepath = easygui.fileopenbox()

input_dist = int(raw_input("Distance Filter Value (in km)?: "))
input_crop = raw_input("what crop?: ")

directory = os.path.split(filepath)[0]

os.chdir(directory)

def dist_haversine(filepath,input_dist,input_crop):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """

    t0 = time()

    r = shapefile.Reader(filepath)
    idx = np.arange(len(r.records()))
    coordinates = []
    for i in idx:
        geom = r.shape(i)
        coordinates.append(geom.points[0])       
    
    acoords = np.array(coordinates)
    
    index = []        
    for r,n,l in itertools.izip(acoords[:,0],acoords[:,1],idx):
        if l in index:
            continue
        else:
            for i,j,k in itertools.izip(acoords[:,0],acoords[:,1], idx):
                if k in index:
                    continue
                
                else:
                
                    lon1=r
                    lat1=n
                    lon2=i
                    lat2=j
                    
                    coord_check = ((lon1 == lon2) & (lat1 == lat2))*1
                    
                    if coord_check == 1:
                        continue
                    
                    else:
                        lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
                        
                        # haversine formula
                        dlon = lon2 - lon1 
                        dlat = lat2 - lat1 
                        a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
                        c = 2 * math.asin(math.sqrt(a)) 
                        km = c*6371 #/1000.0
                    
                    if km < input_dist:
                        if k in index:
                            continue
                        else:
                            index.append(k)

    filterList = [i for j, i in enumerate(coordinates) if j not in index]
            
    df_coords = pd.DataFrame(filterList, columns=['Lon','Lat'])
    
    df_coords.insert(0, 'Crop', input_crop)
    
    print "\ntime elapsed: %.2fs" % (time() - t0)
    
    return df_coords.to_csv(directory + "\\" + "%s_distFilter_%skm.csv" % (input_crop, input_dist), sep=",", index=None)

def dist_vincentry(filepath,input_dist,input_crop):
    
    t0 = time()

    r = shapefile.Reader(filepath)
    idx = np.arange(len(r.records()))
    coordinates = []
    for i in idx:
        geom = r.shape(i)
        coordinates.append(geom.points[0])
    
    acoords=np.array(coordinates)
    index = []
    
    for r,n,l in itertools.izip(acoords[:,0],acoords[:,1],idx):
        if l in index:
            continue
        else:
            for i,j,k in itertools.izip(acoords[:,0],acoords[:,1],idx):
                if k in index:
                    continue
                else:
                    lon1=r
                    lat1=n
                    lon2=i
                    lat2=j
                    
                    coord_check = ((lon1 == lon2) & (lat1 == lat2))*1
                    
                    if coord_check == 1:
                        continue
                    
                    else:
            
                        # Ellipsoid Parameters
                        # need to change the semi-major axis value and flattening ratio to fit definition of any ellipsoid
                        a = 6378137  # semi-major axis
                        f = 1/298.257222101  # inverse flattening
                        b = abs((f*a)-a)  # semi-minor axis
                        L = math.radians(lon2-lon1)
                        U1 = math.atan((1-f) * math.tan(math.radians(lat1)))
                        U2 = math.atan((1-f) * math.tan(math.radians(lat2)))
                        sinU1 = math.sin(U1)
                        cosU1 = math.cos(U1)
                        sinU2 = math.sin(U2)
                        cosU2 = math.cos(U2)
                        lam = L
                    
                        for i in range(100):
                            sinLam = math.sin(lam)
                            cosLam = math.cos(lam)
                            sinSigma = math.sqrt((cosU2*sinLam)**2 +
                                                 (cosU1*sinU2-sinU1*cosU2*cosLam)**2)
                            if (sinSigma == 0):
                                distance = 0  # coincident points
                                break
                            cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLam
                            sigma = math.atan2(sinSigma, cosSigma)
                            sinAlpha = cosU1 * cosU2 * sinLam / sinSigma
                            cosSqAlpha = 1 - sinAlpha**2
                            cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha
                            if math.isnan(cos2SigmaM):
                                cos2SigmaM = 0  # equatorial line
                            C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
                            LP = lam
                            lam = L + (1-C) * f * sinAlpha * \
                                (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma *
                                                     (-1+2*cos2SigmaM*cos2SigmaM)))
                            if not abs(lam-LP) > 1e-12:
                                break
                        uSq = cosSqAlpha * (a**2 - b**2) / b**2
                        A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
                        B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
                        deltaSigma = B*sinSigma*(cos2SigmaM+B/4 *
                                                 (cosSigma*(-1+2*cos2SigmaM*cos2SigmaM) -
                                                  B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma) *
                                                  (-3+4*cos2SigmaM*cos2SigmaM)))
                        km = b*A*(sigma-deltaSigma)

                    if km < input_dist:
                        if k in index:
                            continue
                        else:
                            index.append(k)
                        
    filterList = [i for j, i in enumerate(coordinates) if j not in index]
            
    df_coords = pd.DataFrame(filterList, columns=['Lon','Lat'])
    
    df_coords.insert(0, 'Crop', input_crop)  
    
    print "\ntime elapsed: %.2fs" % (time() - t0)    
    
    return df_coords.to_csv(directory + "\\" + "%s_distFilter_%skm.csv" % (input_crop, input_dist), sep=",", index=None)

dist_haversine(filepath,input_dist,input_crop)
#dist_vincentry(filepath,input_dist,input_crop)
