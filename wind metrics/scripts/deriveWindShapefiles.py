############### deriveWindShapefiles.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: given hourly wind rasters and annual birth residences, extract hourly
#          values for two years and combine them to into a single csv


# Steps include:
# 1) extracting values to points
# 2) extracting date and hour from filename
# 3) combining results to create a single comprehensive csv

############### Setup: import libraries and define constants ##############

import arcpy
import pandas as ps
import numpy as np
import os
import gConst as const
arcpy.env.overwriteOutput = True
from multiprocessing import Pool

YEARS = list(range(2007,2017))
PARENT_FOLDER = const.WIND_FOLDER
BIRTH_FOLDER = PARENT_FOLDER + "Birth_Addresses_Wind/births_by_year/"
WIND_DIRECTION_FOLDER = PARENT_FOLDER + "wind_surfaces/wind_direction/"


################ helper functions ################

# check if 24 rasters are available for one day.  If so, return
# array with filepaths to the rasters
# INPUTS:
#    windDailyPrefix (str) - first half of filepath to the 24 .tifs
# OUTPUTS:
#    array containing filepaths to 24 hourly rasters in a single day
def createDailyFileArray(windDailyPrefix):
    dailyArray = []
    for hour in range(1,25):
        dailyArray.append([windDailyPrefix + "_" + str(hour).zfill(2) + ".tif", str(hour).zfill(2)])
    nCorrect = 0
    for arr in dailyArray:
        if(os.path.exists(arr[0])):
            nCorrect +=1
    if nCorrect == 24:
        return(dailyArray)
    return([])
        
# convert arcgis attribute table to pandas dataframe
# INPUTS:
#    in_fc (str) - filepath to the attribute table/shapefile
# OUTPUTS:
#    fc_dataframe (pandas dataframe) - attribute table in pandas format
def arcgis_table_to_df(in_fc):
    print("converting to dataframe")
    OIDFieldName = arcpy.Describe(in_fc).OIDFieldName
    final_fields = [field.name for field in arcpy.ListFields(in_fc)]
    data = [row for row in arcpy.da.SearchCursor(in_fc,final_fields)]
    fc_dataframe = ps.DataFrame(data,columns=final_fields)
    fc_dataframe = fc_dataframe.set_index(OIDFieldName,drop=True)
    valsList = []

    # for each hour in the day, prefix the variable name with a letter to avoid 
    # a bug issue later with variables that start with numbers
    for i in range(1,25):
        fc_dataframe[str(i).zfill(2)] = fc_dataframe["F" + str(i).zfill(2)].astype(int)
        valsList.append(str(i).zfill(2))
    valsList.append("uniqueid")
    fc_dataframe = fc_dataframe[valsList]
    return fc_dataframe

# for all births within a single year, extract 2 years of hourly wind directions
# INPUTS:
#    year (int) - year of births
def extractPointVals(year):
    birthShapefile = BIRTH_FOLDER + "births_" + str(year) + ".shp"
    masterBirthFile = "memory/births_" + str(year)

    # create a temp copy of the birth shapefile to prevent accidentally altering the original
    arcpy.CopyFeatures_management(birthShapefile, masterBirthFile)
    outputFolder = const.WIND_FOLDER + str(year)
    if not(os.path.exists(outputFolder)):
        os.mkdir(outputFolder)

    # for 2 years of wind direction, starting with the year before birth,
    for curYear in range(year-1,year+1):
        windFolder = WIND_DIRECTION_FOLDER + str(curYear) + "/"

        # for each month in the year and each day in the monthh, extract raster values to points and 
        # save the results in csv format
        for month in range(1,13):
            for day in range(1,32):
                outputFile = outputFolder + "/" + str(curYear) + "_" + str(month) + "_" + str(day) + ".csv"
                if not(os.path.exists(outputFile)):
                    birthCopy = "memory/t_births_" + str(year) + str(month) + str(day) 
                    arcpy.CopyFeatures_management(masterBirthFile, birthCopy)
                    windDayPrefix = windFolder + str(curYear) + "_" + str(month).zfill(2) + "_" + str(day).zfill(2) 
                    dailyArray = createDailyFileArray(windDayPrefix)
                    if(len(dailyArray)==24):
                        try:
                            arcpy.sa.ExtractMultiValuesToPoints(birthCopy,dailyArray)
                            pandasDf = arcgis_table_to_df(birthCopy)
                            pandasDf.to_csv(outputFile,index=False)
                            arcpy.Delete_management(birthCopy)
                        except Exception as e:
                            print("couldn't extract value to points for file %s: %s" %(windDayPrefix,str(e)))


# given a filename with the date in it, extract the month and day
# INPUTS:
#    filename (str) - can be absolute or relative filepath
# OUTPUTS:
#    month (int), day (int)
def getMonthDayFromFilename(filename):
    monthStart = filename.find("_")
    monthEnd = filename[monthStart+1:].find("_")
    month = filename[monthStart+1:][:monthEnd]
    dayStart = filename.rfind("_")
    dayEnd = filename[dayStart+1:].find(".csv")
    day = filename[dayStart+1:][:dayEnd]
    return([month,day])         

# given a csv of measurements, denormalize the variables and store in dataframe format
# INPUTS:
#    folder (str) - absolute filepath to the folder containing the csv
#    csv (str) - relative filepath to the csv file
#    curYear (int) - year of birth
# OUTPUTS:
#    pandas dataframe containing denormalized records
def transformCSVtoDF(folder,csv,curYear):
    df = ps.read_csv(folder + "/" + csv)
    curMonth,curDay =getMonthDayFromFilename(csv)
    hour,day,month,year,wnd_dir,master_id = [[] for x in range(6)]

    # denormalize dataset
    for curHour in range(1,25):
        tempVals = list(df[str(curHour).zfill(2)])
        nVals = len(tempVals)
        wnd_dir += tempVals
        hour += [curHour-1 for x in range(nVals)]
        month += [curMonth for x in range(nVals)]
        day += [curDay for x in range(nVals)]
        year += [curYear for x in range(nVals)]
        master_id += list(df['uniqueid'])
    df = ps.DataFrame({
        'year':year,
        'month':month,
        'day':day,
        'hour':hour,
        'wnd_dir':wnd_dir,
        'master_id':master_id    
    })
    return(df)

# given a single csv for each day, combine all records for a year of births into a single 
# comprehensive csv
# INPUTS:
#    year (int) - year of birth
def combineAnnualCSVs(year):
    outputFile = const.WIND_FOLDER + "/comb_" + str(year) + ".csv"
    if(os.path.exists(outputFile)):
        return
    inputFolder = const.WIND_FOLDER + str(year)
    csvsToCombine = os.listdir(inputFolder)
    pandasData = transformCSVtoDF(inputFolder,csvsToCombine[0],year)
    index = 1
    allPandas = []

    # for each daily csv file of wind direction estimates, add to the dataframe
    for file in csvsToCombine[1:]:
        if(index%30==0):
            allPandas.append(pandasData)
            pandasData = transformCSVtoDF(inputFolder,csvsToCombine[0],year)
        else:
            pandasData = ps.concat([pandasData,transformCSVtoDF(inputFolder,file,year)])
        print(index)
        index+=1
    for tempData in allPandas:
        pandasData = ps.concat([pandasData,tempData])
    pandasData.to_csv(outputFile,index=False)
            




################## Main Fucntion #################


if __name__ == '__main__':
    with Pool(processes=len(YEARS)) as pool:
        res = pool.map(extractPointVals,YEARS)
    
    for year in list(range(2007,2018)):
        combineAnnualCSVs(year)
