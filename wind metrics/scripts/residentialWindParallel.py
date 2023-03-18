############### residentialWindParallel.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: create rose-wind shapefiles (1 degree radial resolution) for each maternal address
# Steps to create rose wind shapefiles include :
# 1) load 2 years of daily wind directions from storage
# 2) load maternal residence metadata including coordinates, birth date, and conception date
# 3) restrict daily wind directions to [conception date, birth date]
# 4) sum daily wind directions for pregnancy period 
# 5) create road wind shapefiles in ArcGIS


############### Setup: import libraries and define constants ##############
sucessfulImport = False
import time

# needed when using the ArcGIS license for a large # of parallel threads
while(sucessfulImport == False):
    try:
        import arcpy
        from arcpy.ia import * 
        from arcpy import env
        from arcpy.sa import *
        arcpy.env.overwriteOutput = True
        sucessfulImport = True
    except Exception as e:
        print("couldn't import arcgis: %s" %(str(e)))
        time.sleep(1)

import numpy as np
import math
import pandas as ps
import os
import datetime
import gConst as const
from multiprocessing import Pool
arcpy.env.overwriteOutput = True

YEAR = 2016
PARENT_FOLDER = const.WIND_FOLDER
masterIdFolder = PARENT_FOLDER + "/masterIds/" + str(YEAR) + '/'
windPartitions = PARENT_FOLDER + "/windPartitions/" + str(YEAR) + "/"
BIRTH_DATA = ps.read_csv(PARENT_FOLDER + "Birth_Addresses_Wind/births_by_year/csvs/births_" + str(YEAR) + ".csv")
OUTPUT_SHAPEOLDER = PARENT_FOLDER + "/outputShapefiles/" + str(YEAR)
if not (os.path.exists(OUTPUT_SHAPEOLDER)):
    os.mkdir(OUTPUT_SHAPEOLDER)

sr = arcpy.SpatialReference()   # WGS84 spatial reference
sr.factoryCode = 4326
sr.create()


# the following components are used for the Haversine formula.  For more details go to
# http://www.movable-type.co.uk/scripts/latlong.html
# https://www.eol.ucar.edu/content/wind-direction-quick-reference
EARTH_RADIUS = 6378140.0                # used to calculate wind direction
distance = 55000.0/EARTH_RADIUS         # used to calculate bearing and wind direction
N_ANGLES = 360



######################### HELPER FUNCTIONS #######################


# calculate latitude and longitude coordinates on the outside of the buffer region
# INPUTS:
#   origLat (float) - latitude coordinate
#   origLong (float) - longitude coordinate
#   bearing (float) - angular direction relative to true north
#   distance (float) - distance that new coordinates should be from the center
# OUTPUTS:
#   newLatit (float) - latitude coordinate of outer buffer edge
#   newLongit (float) - longitude coordinate of outer buffer edge
def calcCoords(origLat,origLong,bearing,distance):    
    latRadians = origLat*0.0174533          # convert latitude from degrees to radians
    longRadians = origLong*0.0174533        # convert longitude from degrees to radians
    
    # calculate the new latitude coordinate for a point that is 'distance' units away using the haversine formula (see weblink above)
    newLatit = math.asin(math.sin(latRadians)*math.cos(distance) + math.cos(latRadians)*math.sin(distance)*math.cos(bearing))/0.0174533
    
    # calculate the new longitude coordinate for a point that is 'distance' units away using the haversine formula
    longitP1 = math.sin(bearing)*math.sin(distance)*math.cos(latRadians)
    longitP2 = math.cos(distance) - math.sin(latRadians)*math.sin(newLatit*0.0174533)    
    newLongit = longRadians + math.atan2(longitP1,longitP2)
    newLongit = newLongit/0.0174533
    return([newLatit,newLongit])

# create a triangle polygon in arcgis using 3 coordinate points.  The triangle will represent one angular degree of coverage
# INPUTS:
#   wndDrt (float) - angular wind direction, with 0 degrees corresponding to true north
#   strtLat (float) - latitude of the centroid
#   strtLong (float) - longitude of the centroid
#   inputFilename (string) - filepath to the shapefile to which the triangle will be added
#   timestamp (int) - time that corresponds to the wind direction, in epoch format
def calcTriangle(wndDrt,strtLat, strtLong,distance,inputFilename,colNames,colVals):
    bearing1 = ((wndDrt + 0.5)%360.0)*0.0174533       # center of triangle + 0.5 degrees
    bearing2 = ((wndDrt - 0.5)%360.0)*0.0174533       # center of triangle - 0.5 degrees
    airMonitorLoc = arcpy.Point(strtLong,strtLat)       # defnie the air monitor location as one point of the triangle
    coordP1 = calcCoords(strtLat,strtLong,bearing1,distance)  # create another coordinate point
    upPoint1 = arcpy.Point(coordP1[1],coordP1[0])
    coordP2 = calcCoords(strtLat,strtLong,bearing2,distance)  # create another coordinate point
    upPoint2 = arcpy.Point(coordP2[1],coordP2[0])
    # Create a polygon geometry
    array = arcpy.Array([airMonitorLoc,upPoint1,upPoint2])
    polygon = arcpy.Polygon(array)
    # Open an InsertCursor and insert the new geometry
    cursor = arcpy.da.InsertCursor(inputFilename, ['SHAPE@'] + colNames)
    rowVals = tuple([polygon] + list([colVals]))
    cursor.insertRow(rowVals)
    # Delete cursor object
    del cursor

# add syntax necessary to create multiple attributes in an ArcGIS attribute table
# INPUTS:
#    fieldName (str) - name of the attribute field to add
# OUTPUTS:
#    array with [fieldName, datatype]
def mapCreateFieldNames(fieldName):
    return( [fieldName,"Short"])

# given a specific date and year, get indeces the daily wind direction is stored in the 2 year wind direciton array
# INPUTS:
#    inDate (str) - date of interest in Y-m-d format
#    bYear (int) - birth year
# OUTPUTS: 
#    array indeces where the daily wind direction is stored
def getBirthIndexes(inDate,bYear):
    datetime_object = datetime.datetime.strptime(inDate, '%Y-%m-%d')
    monthIndex = (datetime_object.month + 11) if datetime_object.year == bYear else (datetime_object.month -1)
    dayIndex = datetime_object.day - 1
    return([monthIndex,dayIndex])

# given data for a single birth, restrict wind values to the pregnancy period and create a wind rose matrix
# INPUTS:
#    windData (numpy array) - 2 years of daily wind direction for a single maternal residence
#    bdate (str) - birth date in Y-m-d format
#    cdate (str) - conception date in Y-m-d format
#    year (int) - birth year
# OUTPUTS:
#    wind rose numpy array - # of hours maternal residence is downwind of each 1 radial degree (360 total)
def processAnnualVals(windData,bdate,cdate,byear):
    bmonth,bday = getBirthIndexes(bdate,byear)
    cmonth,cday = getBirthIndexes(cdate,byear)
    startVals = np.sum(np.maximum(windData[cmonth][cday:],0),axis=0)
    startVals += np.sum(np.maximum(windData[bmonth][:bday],0),axis=0)
    for monthData in windData[cmonth+1:bmonth]:
        startVals += np.sum(np.maximum(monthData,0),axis=0)
    return(startVals)

# given maternal residence information and a wind rose matrix, create a rose wind shapefile in ArcGIS
# INPUTS:
#    lat (float) - maternal residence latitude
#    lon (float) - maternal residence longitude
#    angleVals (numpy array) - wind rose matrix for pregnancy period
#    uniqueId (str) - unique identifier for the birth record
def calcAnnualShapefile(lat,lon,year,angleVals,uniqueId):
    colNames = ['sum_' + str(year)]
    fieldNames = list(map(mapCreateFieldNames,colNames))
    yearShapefile = uniqueId + ".shp"
    arcpy.CreateFeatureclass_management(OUTPUT_SHAPEOLDER,yearShapefile,"POLYGON",'#','#','#',sr)
    arcpy.management.AddFields(OUTPUT_SHAPEOLDER + "/" + yearShapefile,fieldNames)
    for angle in range(360):
        calcTriangle(angle+0.5,lat,lon,distance,OUTPUT_SHAPEOLDER + "/" + yearShapefile,colNames,angleVals[angle])

# given 2 years of daily wind direction and all metadata about the maternal residence and birth date, 
# perform all steps necessary to create a wind rose shapefile
# INPUTS:
#    dataTuple (tuple) - contains wind rose matrix, maternal resdience info, birth date, 
#                        conception date, and unique identifier
def processSingleResidence(dataTuple):
    outputFile = OUTPUT_SHAPEOLDER + "/" + dataTuple[5] + ".shp"
    if not(os.path.exists(outputFile)):
        sumWindVals = processAnnualVals(dataTuple[0],dataTuple[1],dataTuple[2],YEAR)
        calcAnnualShapefile(dataTuple[3],dataTuple[4],YEAR,sumWindVals,dataTuple[5])

# given multiple sources of raw data, repackage the data to facilitate parallel processing 
# INPUTS:
#    windData (numpy array) = 2 years of daily wind direction for multiple maternal residences
#    idData (numpy array) - unique identifiers, sorted in the same order as the windData
#    personalData (pandas dataframe) - maternal residence and birth info
# OUTPUTS:
#    list of tuples - each tuple contains all info needed to create a wind rose for 1 maternal residence
def prepParallel(windData,idData,personalData):
    parallelArray = []
    for i in range(idData.count()[0]):
        curId = list(idData.iloc[i])
        curPersonal = personalData[personalData['uniqueid'].isin(curId)]
        if(int(curPersonal.count()[0])>0):
            curTuple = (windData[i],list(curPersonal['bdate'])[0],list(curPersonal['cdate'])[0],
            list(curPersonal['b_lat'])[0],list(curPersonal['b_long'])[0],curId[0])
            parallelArray.append(curTuple)
    return(parallelArray)



################# Main function ######################
if __name__ == '__main__':
    masterIds = os.listdir(masterIdFolder)
    for index in range(len(masterIds)):
        windData = np.load(windPartitions + "w_ " + str(index) + ".npy")
        idData = ps.read_csv(masterIdFolder + "id_" + str(index) + ".csv")
        parallelData = prepParallel(windData,idData,BIRTH_DATA)
        print("finished prepping data")
        with Pool(processes=64) as pool:
            pool.map(processSingleResidence,parallelData)