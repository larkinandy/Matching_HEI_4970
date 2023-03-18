############## deriveMatchParallel.py ###############
# Author: Andrew Larkin
# Date Created: January 24th, 2022
# Summary: Given shapefiles of residents in the bottom and top quartile of time spent downwind during prengnacy,
# Match the residents to the max downwind road (the road that's creating the biggest wind exposure for the exposed residence), 
# and calculate match score criteria metrics, including difference in distance to the max downwind road between exposed and control,
# difference in birth year between exposed and control, and difference in distance to the nearest road


############### Setup: import libraries and define constants ##############
sucessfulImport = False
import time
import os
from multiprocessing import Pool
import numpy as np
import pandas as ps
import random
import gConst as const

# needed when using the ArcGIS license for a large # of parallel threads
while(sucessfulImport == False):
    try:
        #print("try again")
        import arcpy
        from arcpy.ia import * 
        from arcpy import env
        from arcpy.sa import *
        arcpy.env.overwriteOutput = True
        arcpy.CheckOutExtension("Spatial")
        sucessfulImport = True
    except Exception as e:
        print("couldn't import arcgis: %s" %(str(e)))
        time.sleep(1)


# Define global constants used by all parallel processing threads
BUFFER_DISTANCE = 500
PARENT_FOLDER = const.WIND_ESTIMATES_FOLDER
BIRTHS_FOLDER = PARENT_FOLDER + "/births_by_year/"
MATCHING_FOLDER = PARENT_FOLDER + "matching/"
OLDER_TRAFFIC = const.TRAFFIC_FOLDER + "rhino_vmt_1995_2009.shp"
NEWER_TRAFFIC = const.TRAFFIC_FOLDER + "rhino_vmt_2010_2016.shp"



# reformat exposed cohort data from a table to tuples to facilitate high throughput parallel processing.  
# each tuple can be independently sent to a thread for parallel processing
# INPUTS:
#    exposed (pandas dataframe) - contains a unique id, birth year, and estimated gestational age
# OUTPUTS:
#    the input data transformed into an array of tuples
def prepParallel(exposed):
    dataTuples = []
    nExposed = exposed.count()[0]
    print("n exposure records: %i" %(nExposed))

    # create one tuple of data, with one tuple sent as inputs to each thread in a parallel processing acrhitecture
    for recordIndex in range(nExposed):
        curRecord = exposed.iloc[recordIndex]
        curId = curRecord['uniqueid']
        curYear = curRecord['byear']
        curCutoff = curRecord['b_es_ges']*18.9 # to save on computational cost, only consider road segments which empirically are in the top 50% of time upwind from residence
        dataTuples.append((curId,curYear,curCutoff))
    return(dataTuples)


# for a single maternal residence in the exposed group (i.e. highest quartile of wind epoxusre),
# identify all control residences that are close enough to be a match and calculate distance form control points to road segements
# INPUTS:
#    uniqueId (str) - unique id for each residence
#    year (str) - birth year of the exposed maternal residence
#    wndCutofff (float) - only consider road segments with are in the top 50% of time upwind from exposed residence
#    buffer (int) - maximum distance roads can be from exposed
# OUTPUTS:
#    distance from control points to nearby roads are stored in the file path defined by the variable controlsNearRd
#    distance fomr exposed residence to nearest road is directly returned by the function
def processOneResidence(uniqueId,year,wndCutoff,buffer):

    # name of the file that will contain the results
    controlsNearRd = MATCHING_FOLDER + "ctrlsNrRd_" + uniqueId + ".csv"

    # roads near the exposed residence, with hours downwind attached to each road segment.  There's a unique shapefile for each exposed residence
    roadWndShapefile = PARENT_FOLDER + "shpFile/" + str(year) + "/" + uniqueId + ".shp"

    # all control points, not just the ones near the exposed residence
    controlPoints = MATCHING_FOLDER + "/bottom_" + str(buffer) + ".shp"

    # roads that are selected because wind exposure to the maternal residence is above the wndCutoff
    selectedRds = MATCHING_FOLDER + "tmpRds" + uniqueId + ".shp"

    # road broken into 10m segments
    selectedRdsDissolved =  MATCHING_FOLDER + "tmpRds2" + uniqueId  + ".shp"

    # query used to select roads with wind exposure above windCutoff
    expression = '"' + uniqueId + "_su" + '"' + '>=' + str(wndCutoff) + ' AND ' + '"a_bufferDi" <=' + str(buffer)

    # select roads above wind cutoff and convert to numpy array
    try:
        arcpy.Select_analysis(roadWndShapefile,selectedRds,expression)
    except Exception as e:
        print(str(e))
    arr = arcpy.da.FeatureClassToNumPyArray(selectedRds, 'a_bufferDi')
    arr = [item for sublist in arr for item in sublist]
    distToRoad = min(arr) # get distance to nearest road

    # difference between exposed and control in distance to the matched road cannot be more than 0.5*buffer
    maxSearchDistance = distToRoad+0.5*buffer

    # break road network into 10m segements
    arcpy.management.Dissolve(selectedRds,selectedRdsDissolved)

    # get distance from 10m road segemetns to all nearby control points.  Results are stored in table
    arcpy.GenerateNearTable_analysis(
            selectedRdsDissolved,controlPoints,controlsNearRd,str(maxSearchDistance) + " Meters", "NO_LOCATION",
            "NO_ANGLE","ALL","#","GEODESIC"
        )
    arcpy.Delete_management(selectedRdsDissolved)
    return(distToRoad)



# calculate match scores for candidate control matches near a single exposed residence
# sort the results so the best match is near the top
# INPUTS:
#    controlData (dataframe) - contains metrics of distance from exposed and controls to roads, and birth years
#    maxDiff (int) - the threshold for the largest difference between exposed and control that will be considered 
#                    a viable match
# OUTPUTS:
#    controlData (dataframe) - scores appended to the input data, and sorted with the best match at the top
def screenAndRankControls(controlData,maxDiff):
    controlData['dist_diff'] = controlData['exp_dist'] - controlData['ctrl_dist']
    controlData = controlData[controlData['dist_diff'] < maxDiff]
    controlData = controlData[controlData['dist_diff'] > -1*maxDiff]
    controlData['year_diff'] = controlData['ctrl_year'] - controlData['exp_year']
    controlData['abs_dist'] = controlData['dist_diff'].abs()
    controlData['abs_year'] = controlData['year_diff'].abs()
    controlData.sort_values(by=['abs_dist','abs_year'],inplace=True)
    return(controlData)



# calculate the difference in distance to match road between exposure residence and nearby controls
# INPUTS:
#    bufferSize (int) - max distance from residence to road
#    controlData (pandas dataframe) - contains records of all controls across texas
#    expId (str) - unique identifier for the exposed residence
#    expDist (str) - distance from exposed residence to nearest road
#    expYear (int) - birth year at exposed residence
# OUTPUTS:
#    results are stored in a csv file at the aboslute filepath defined by the variable 'outputCSVFilepath'
def processNearCandidates(bufferSize,controlData,expId,expDist,expYear):

    # where results will be stored
    outputCSVFilepath = MATCHING_FOLDER + str(BUFFER_DISTANCE) + "/" +  expId + ".csv"
    if(os.path.exists(outputCSVFilepath)):
        return
    
    # load file containing a list of control residence near the exposed residence road network
    # and join with the more comprehensive data about control residence
    controlsNearRd = MATCHING_FOLDER + "ctrlsNrRd_" + expId + ".csv"
    nearCandDist = ps.read_csv(controlsNearRd)
    nearCandDist = nearCandDist.merge(controlData, left_on='NEAR_FID', right_on='FID')

    # select only the columns of interest and rename them
    nearCandDist = nearCandDist[nearCandDist.columns[nearCandDist.columns.isin(['NEAR_DIST', 'byear','uniqueid'])]]
    nearCandDist.columns = ['ctrl_dist','ctrl_id','ctrl_year']

    # if there are no candidate control matches nearby then stop processing the exposed residence
    nCompares = nearCandDist.count()[0]
    if(nCompares==0):
        return

    # crate a matrix structure that will allow us to subtract exposure values from control values by substracting two arrays
    expId = [expId for x in range(nCompares)]
    expDist = [expDist for x in range(nCompares)]
    expYear = [expYear for x in range(nCompares)]
    nearCandDist['exp_id'] = expId
    nearCandDist['exp_dist'] = expDist
    nearCandDist['exp_year'] = expYear

    # calculate match quality scores and sort candidate controls by score
    controlCandidates = screenAndRankControls(nearCandDist,bufferSize*0.5)
    controlCandidates.to_csv(outputCSVFilepath,index=False)

    # clean up
    selectedRds = MATCHING_FOLDER + "tmpRds" + expId[0] + ".shp"
    arcpy.Delete_management(selectedRds)
    arcpy.Delete_management(controlsNearRd)


# find the controls near an exposed residence and calculate match criteria metrics, including distance from control to the match road and 
# distance to the nearest road
# INPUTS:
#    dataTuple:
#       1) unique id of the exposed residence to match to (str)
#       2) birth year of the exposed residence to match to (int)
#       3) minimum hours upwind for road segements that should be considered (float)
# OUTPUTS:
#     results are stored in csv file defined by the variable 'outputCSV'
def matchOneResidence(dataTuple):
    outputCSV = MATCHING_FOLDER + str(BUFFER_DISTANCE) + "/" + dataTuple[0] + ".csv"
    if(os.path.exists(outputCSV)):
        print("id %s already processed" %(dataTuple[0]))
        return
    
    # identify all control points close enough to be a candidate match and calculate distances 
    # from control points to road segements near the exposed 
    try:
        control = ps.read_csv(MATCHING_FOLDER + "bottom_" + str(BUFFER_DISTANCE) + ".csv")
        distToRoad = processOneResidence(dataTuple[0],dataTuple[1],dataTuple[2],BUFFER_DISTANCE)
    except Exception as e:
        print("couldn't calc dist to road: " + str(e))
        return
    
    # calculate the difference in distance to match road between exposure residence and nearby controls
    try:
        processNearCandidates(BUFFER_DISTANCE,control,dataTuple[0],distToRoad,dataTuple[1])
    except Exception as e:
        print("couldn't compare cats: " + str(e))

    print("completed id %s" %(dataTuple[0]))


if __name__ == '__main__':
    
    # prepare data for high through parallel analyses    
    exposed = ps.read_csv(MATCHING_FOLDER + "top_" + str(BUFFER_DISTANCE) + ".csv")
    print(exposed.head())
    dataTuples = prepParallel(exposed)
    random.shuffle(dataTuples)
    print("completed prepping data for paralell processing")

    # run matches, 96 maternal residences at a time on 96 threads
    pool = Pool(processes=96) # using 96/128 threads seems to work best on a 64-core multithreaded workstation
    res = pool.map_async(matchOneResidence,dataTuples)
    res.get()

# end of deriveMatchParallel.py