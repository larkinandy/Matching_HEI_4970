############### downloadWindVectors.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: 
# 1) download daily 10m u- and v- wind vector components from the ERA5 reanalysis product
# 2) calculate daily wind speed from u and v wind vectors
# 3) reference: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation


############### Setup: import libraries and define constants ##############
#import cdsapi
import os
import math
import numpy as np
from multiprocessing import Pool 
import arcpy
import gConst as const
arcpy.env.overwriteOutput = True
N_PROCESSES = 96 # number of parallel processes to run.  Empirically 96 is ideal for a 64-core multithreaded processor

# required API to automate ERA5 data download.   See https://cds.climate.copernicus.eu/api-how-to
#c = cdsapi.Client()

PARENT_FOLDER = const.WIND_FOLDER
WIND_VECTOR_FOLDER = PARENT_FOLDER + "wind_surfaces/wind_components/"
WIND_DIRECTION_FOLDER = PARENT_FOLDER + "wind_surfaces/wind_direction/"
DAYS_IN_MONTH = [31,28,31,30,31,30,31,31,30,31,30,31]
LEAP_YEARS = [2008,2012,2016]
COMPONENT_DICT = {
    'u_comp': '10U@SFC',
    'v_comp': '10V@SFC'
}
YEARS = list(range(2006,2007))
PROJECTION = "GEOGCS['GCS_Coordinate_System_imported_from_GRIB_file',DATUM['D_unknown',SPHEROID['Sphere',6367470.0,0.0]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"


################ Helper functions called by the main function ##############



# download hourly wind components for one 24 hour period
# INPUTS:
#    year (int) - year of interest
#    month (int) - month of interest
#    day (int) - day of interest
#    filename (string) - full filepath where the downloaded data should be stored
def downloadHourlyDataOneData(year,month,day,filename):
    d = c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind',
        ],
        'year': str(year),
        'month': str(month).zfill(2),
        'day': str(day).zfill(2),
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'format': 'grib',
        'area': [
            40, -120, 20, 
            -90,
        ],
    },
    filename)

# download one year of wind components (hourly temporal resolution)
# INPUTS:
#    year (int) - year of interest
def downloadSingleYearData(year):
    yearFolder = WIND_VECTOR_FOLDER + str(year) + "/"
    if not(os.path.exists(yearFolder)): 
        os.mkdir(yearFolder)

    # for each month in the year
    for month in list(range(1,13)):
        curDays = DAYS_IN_MONTH[month-1]
        if(year % 4 == 0 and month == 2): # adjust days in Feb for leap year if necessary
            curDays +=1
        
        # for each day in the month
        for day in list(range(1,curDays+1)):
            outputFile = yearFolder + str(year) + "_" + str(month).zfill(2) + "_" + str(day).zfill(2) + ".grib"
            if not os.path.exists(outputFile):
                downloadHourlyDataOneData(year,month,day,outputFile)



# create string to extract time value from multivariate timeseries raster
# INPUTS:
#    year (int) - year of interest
#    month (int) - month of interest
#    day (int) - day of interest
#    hour (int) - hour of interest
# OUTPUTS:
#    queryString (str) - string in format needed to extract hourly value from timeseries raster
def createTimeRangeString(year,month,day,hour):
    prefix = "StdTime "
    dateString = str(year) + "-" + str(month).zfill(2) + "-" + str(day).zfill(2)
    hourString = 'T' + str(hour-1).zfill(2) + ":00:00"
    return(prefix + dateString + hourString)

# create wind direction raster from u- and v- wind component rasters
# INPUTS:
#    uFile (string) - full filepath to the raster containing u-components
#    vFile (string) - full filepath to the raster containing v-components
# OUTPUTS:
#    windRaster (obj) - arcpy object of wind direction raster estimates
def calcWindDirection(uFile,vFile):
    u_comp = arcpy.Raster(uFile)
    v_comp = arcpy.Raster(vFile)
    uArr = arcpy.RasterToNumPyArray(u_comp)
    vArr = arcpy.RasterToNumPyArray(v_comp)
    windArr = 180 + (180/math.pi)*np.arctan2(uArr,vArr)
    windRaster = newRaster = arcpy.NumPyArrayToRaster(
        windArr,
        arcpy.Point(u_comp.extent.XMin,u_comp.extent.YMin),
        u_comp.meanCellWidth
    )
    return(windRaster)

# extract u- and v- wind components for a single hour from a 
# a multivariate timeseries raster and create a wind direction
# INPUTS:
#    inputFile (str) - full filepath to the multivariate raster containing u- and v-components
#    year (int) - year of interest
#    month (int) - month of interest
#    day (int) - day of interest
#    hour (int) - hour of interest
#    outputFile (str) - full filepath to the location where the wind direction raster will be saved
def convertComponentsToDirectionOneHour(inputFile,year,month,day,hour,outputFile):
    tempUComponent = PARENT_FOLDER + "temp_u_" + str(year) + ".crf"
    tempVComponent = PARENT_FOLDER + "temp_v_" + str(year) + ".crf"

    # extract u- and v-componenets from multivariate raster
    queryString = createTimeRangeString(year,month,day,hour)
    arcpy.md.SubsetMultidimensionalRaster(
        inputFile,tempUComponent,COMPONENT_DICT['u_comp'],
        "BY_VALUE",None,queryString,
        '', '', '', None, ''
    )
    arcpy.md.SubsetMultidimensionalRaster(
        inputFile,tempVComponent,COMPONENT_DICT['v_comp'],
        "BY_VALUE",None,queryString,
        '', '', '', None, ''
    )

    # create wind direction raster and save
    windDirection = calcWindDirection(tempUComponent,tempVComponent)
    windDirection.save(outputFile)
    arcpy.management.DefineProjection(outputFile, PROJECTION)

    # cleanup temporary files
    arcpy.Delete_management(tempUComponent)
    arcpy.Delete_management(tempVComponent)
    

# convert hourly u- and v- wind components for an entire year to wind direction estimates
# save wind direction estimates as .tif files, one file for each hour in the year
# INPUTS:
#    year (int) - year of interest
def convertComponentsToDirectionOneYear(year):
    
    # setup
    inputFolder = WIND_VECTOR_FOLDER + str(year) + "/"
    outputFolder = WIND_DIRECTION_FOLDER + str(year) + "/"
    if not(os.path.exists(outputFolder)):
        os.mkdir(outputFolder)
    
    # for each month in the year
    for month in list(range(1,13)):
        curDays = DAYS_IN_MONTH[month-1]
        if(year % 4 == 0 and month == 2): # adjust days in Feb for leap year if necessary
            curDays +=1
        
        # for each day in the month
        for day in list(range(1,curDays+1)):
            print("processing day %s" %(day))

            # for each hour of the day
            for hour in list(range(1,25)):
                outputFile = outputFolder + str(year) + "_" + str(month).zfill(2) + "_" + str(day).zfill(2) + "_" + str(hour).zfill(2) + ".tif"

                if not os.path.exists(outputFile):
                    inputFile = inputFolder + str(year) + "_" + str(month).zfill(2) + "_" + str(day).zfill(2) + ".grib"
                    convertComponentsToDirectionOneHour(inputFile,year,month,day,hour,outputFile)



################# Main function ######################

if __name__ == '__main__':
    for year in YEARS:
        downloadSingleYearData(year)
    with Pool(processes=N_PROCESSES) as pool:
        res = pool.map(convertComponentsToDirectionOneYear,YEARS)