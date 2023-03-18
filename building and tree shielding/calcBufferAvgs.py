############### calculateBufferAvgs.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: given a shapefile containing wind, road, and shielding estimates,
#          calculate buffer averages for each exposure metric

# Steps include:
# 1) exporting shapefile into csv format
# 2) renaming variables to match the final data dictionary
# 3) calculating buffer averages 
# 4) exporting buffer averages to csv


############### Setup: import libraries and define constants ##############
import os
from multiprocessing import Pool
sucessfulImport = False
import time
import numpy as np
import pandas as ps
import sys
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


YEAR = sys.argv[1]



# given one shapefile containing exposure metrics for a single materal residence,
# transform results to a dataframe and then export to a csv
# also, reanme variables to match the final data dictionary
# INPUTS:
#    shpFilepath (str) - absolute filepath to the shapefile to process
def processShp(shpFilepath,birthData):
    outputFile = const.WIND_FOLDER + "bufferAvgs/" + str(YEAR) + "/" + shpFilepath[:-4] + ".csv"
    print(outputFile)
    if os.path.exists(outputFile): return
    arcpy.management.CalculateGeometryAttributes(const.WIND_FOLDER + "shpFile/" + str(YEAR) + "/" + shpFilepath, "rdsSplit_4 LENGTH_GEODESIC", "METERS", '', None, "SAME_AS_INPUT")
    # todo: Update filepath to work with maternal resiences rather than air monitors
    curRecord = (birthData[birthData['uniqueid'] == shpFilepath[:-4]]).iloc[0]
    n_hours_preg = int(curRecord['b_es_ges'])*24*7
    if(n_hours_preg==0):return
    # convert attribute table of the shapefile into a numpy array
    npArr = arcpy.da.TableToNumPyArray(
        const.WIND_FOLDER + "shpFile/" + str(YEAR) + "/" + shpFilepath,
        ['rdsSplit_4','rdsSplit_1','a_trSh','a_bufferDi','a_bldgSh',shpFilepath[:-4] + '_su']
    )
    df = ps.DataFrame(npArr)

    # rename columns
    df.columns = ['rd_len','rd_type','tr_sh','b_dist','bld_sh','an_wnd']
    df['rd_tr'] = df['rd_len']*((100-df['tr_sh'])/100.0)
    df['rd_bl'] = df['rd_len']*((100-df['bld_sh'])/100.0)
    df['rd_wnd'] = df['rd_len']*df['an_wnd']

    # get unique buffers.  the bounds are [10,500], but the subset can differ for each maternal residence!
    
    #uniqueBuffers = list(set(df['b_dist']))
    #uniqueBuffers.sort()
    uniqueBuffers = [50,100,200,300,500]
    uniqueRds = list(set(df['rd_type']))
    uniqueRds.sort()
    outputs = [[] for x in range(15)]
    
    # for each buffer, extract values and store in arrays to later create a dataframe with
    for buffer in uniqueBuffers:
        tempSet = df[df['b_dist']<=buffer]
        for road in uniqueRds:
            tempSet2 = tempSet[tempSet['rd_type']==road]
            if(tempSet2.count()[0]>0):
                outputs[0].append(buffer)
                outputs[1].append(road)
                totalRd = np.sum(tempSet2['rd_len'])
                totalWnd = np.sum(tempSet2['rd_wnd'])
                pWind = 100*tempSet2['an_wnd']/n_hours_preg
                percentWnd = (1.0/(totalRd))*np.sum(pWind*tempSet2['rd_len'])
                maxVal = np.max(tempSet2['an_wnd'])
                if(maxVal==0):
                    maxWind = 0
                else:
                    maxWind = (100.0/n_hours_preg)*maxVal
                cutVal = maxWind*0.9
                cutSubset = tempSet2[pWind >=cutVal]
                cutLen = np.sum(cutSubset['rd_len'])
                cutBuilding = np.mean(cutSubset['bld_sh'])
                cutTree = np.mean(cutSubset['tr_sh'])
                outputs[2].append(totalRd)
                outputs[3].append(np.sum(tempSet2['rd_tr']))
                outputs[4].append(np.sum(tempSet2['rd_bl']))
                outputs[5].append(totalWnd)
                outputs[6].append(np.mean(tempSet2['tr_sh']))
                outputs[7].append(np.mean(tempSet2['bld_sh']))
                outputs[8].append(totalWnd/totalRd)
                outputs[9].append(percentWnd)
                outputs[10].append(maxWind)
                outputs[11].append(cutVal)
                outputs[12].append(cutLen)
                outputs[13].append(cutBuilding)
                outputs[14].append(cutTree)
        tempSet2 = tempSet[tempSet['rd_type'].isin([0,1])]
        if(tempSet2.count()[0]>0):
            outputs[0].append(buffer)
            outputs[1].append(4)
            totalRd = np.sum(tempSet2['rd_len'])
            totalWnd = np.sum(tempSet2['rd_wnd'])
            pWind = 100*tempSet2['an_wnd']/n_hours_preg
            percentWnd = (1.0/(totalRd))*np.sum(pWind*tempSet2['rd_len'])
            maxVal = np.max(tempSet2['an_wnd'])
            if(maxVal==0):
                maxWind = 0
            else:
                maxWind = (100.0/n_hours_preg)*maxVal
            cutVal = maxWind*0.9
            cutSubset = tempSet2[pWind >=cutVal]
            cutLen = np.sum(cutSubset['rd_len'])
            cutBuilding = np.mean(cutSubset['bld_sh'])
            cutTree = np.mean(cutSubset['tr_sh'])
            outputs[2].append(totalRd)
            outputs[3].append(np.sum(tempSet2['rd_tr']))
            outputs[4].append(np.sum(tempSet2['rd_bl']))
            outputs[5].append(totalWnd)
            outputs[6].append(np.mean(tempSet2['tr_sh']))
            outputs[7].append(np.mean(tempSet2['bld_sh']))
            outputs[8].append(totalWnd/totalRd)
            outputs[9].append(percentWnd)
            outputs[10].append(maxWind)
            outputs[11].append(cutVal)
            outputs[12].append(cutLen)
            outputs[13].append(cutBuilding)
            outputs[14].append(cutTree)
        

        ####################################################################################################
        # HERE IS THE SPOT TO RENAME VARIABLES ONCE THE FINAL NAMES FOR THE DATA DICTONARY HAVE BEEN DECIDED
        tempDict = {
        'buffer':outputs[0],        # buffer distance
        'roadType':outputs[1],      # road category
        'roadLen':outputs[2],       # road length
        'tree':outputs[3],          # tree coverage between roads and maternal residence
        'bldg':outputs[4],          # building coverage between roads and maternal residence
        'wind':outputs[5],          # number of hours maternal residence is downwind of roads
        'treeShield':outputs[6],    # percent roads shielded by trees
        'bldgShield':outputs[7],    # percent roads shielded by buildings
        'annWind':outputs[8],       # annual wind direction
        'percWind':outputs[9],
        'maxWind':outputs[10],
        'cutoff':outputs[11],
        'cutoffLen':outputs[12],
        'cutoffBuilding':outputs[13],
        'cutoffTree':outputs[14]
    }
    print(outputFile)
    ps.DataFrame(tempDict).to_csv(outputFile,index=False)


# test if file has a csv extension
# INPUTS:
# filename (str) - can be absolute or relative filename
# OUTPUTS:
#    returns True if extension is csv, False otherwise
def isCSVFile(filename):
    return(filename[-3:]=="csv")

# test if file has an shp extension
# INPUTS:
#    filename (str) - can be absolute or relative filename
# OUTPUTS:
#    returns True if extension is shp, False otherwise
def isShapefile(filename):
    return(filename[-3:]=="shp")

# given a folder, return files within the folder that have a csv extension
# INPUTS:
#    folder (str) - absolute filepath to a folder
# OUTPUTS:
#    filteredList (str array) - relative filepaths to csv files within the folder 
def getCSVFilesToAvg(folder):
    candFiles = os.listdir(folder)
    csvfileBool = np.array(list(map(isCSVFile,candFiles)))
    filteredList = list(np.array(candFiles)[csvfileBool])
    return(filteredList)


# given a folder, return files within the folder that have a shp extension
# INPUTS:
#    folder (str) - absoluste filepath to a folder
# OUTPUTS:
#    filteredList (str array) - relative filepaths to shp files wihtin the folder
def getShpFilesToConvert(folder):
    candFiles = os.listdir(folder)
    shapefileBool = np.array(list(map(isShapefile,candFiles)))
    filteredList = list(np.array(candFiles)[shapefileBool])
    return(filteredList)

# convert all shapefiles within a folder to csv format
# INPUTS:
#     folder (str) - aboslute filepath to a folder
def convertShpsToCSV(folder,birthData):
    shpList = getShpFilesToConvert(folder)
    for file in shpList:
        processShp(file,birthData)

# if a maternal residence dataframe has no roads (and other exposure metrics)
# for a given buffer distance, add a row with zeros to make it easier to c
# combine with other maternal residence dataframes
# INPUTS:
#    bufferDist (int) - buffer distance with not metrics
#    roadType (int) - category of road
def createZero(bufferDist,roadType):
    zeroData = {
        'buffer':[bufferDist],
        'roadType':[roadType],
        'roadLen':[0],
        'tree':[0],
        'bldg':[0],
        'wind':[0],
        'treeShield':[0],
        'bldgShield':[0],
        'annWind':[0],
        'percWind':[0],
        'maxWind':[0],
        'cutoff':[0],
        'cutoffLen':[0],
        'cutoffBuilding':[0],
        'cutoffTree':[0]
    }
    df = ps.DataFrame(zeroData)
    return(df)

# given a maternal residence dataframe and a buffer distance, denormalize the dataframe by combining
# the roadtype with the exposure type (e.g. trees, building)
# INPUTS:
#    dataset (pandas dataframe) - contains all exposures for one maternal residence
#    buffferDist (int) - buffer distance to extract values for
# OUTPUTS:
#    exactBuffer (pandas dataframe) - denormalized and renamed buffer exposures
def getBufferValues(dataset,bufferDist):

    exactBuffer = dataset[dataset['buffer']==bufferDist]
    exactBuffer = exactBuffer[exactBuffer['roadType']!=2]

    coveredRoads = list(exactBuffer['roadType'])

    for roadType in [0,1,2,3]:
        if roadType not in coveredRoads:
            otherRoads = dataset[dataset['roadType']==roadType]
            otherRoadsSmall = otherRoads[otherRoads['buffer']<bufferDist]
            if otherRoadsSmall.count()[0] == 0:
                zeroDF = createZero(bufferDist,roadType)
                if(exactBuffer.count()[0]==0):
                    exactBuffer = zeroDF
                else:
                    exactBuffer = exactBuffer.append(zeroDF)
            else:
                closestBuffer = max(list(otherRoadsSmall['buffer']))
                closestDF = otherRoadsSmall[otherRoadsSmall['buffer']==closestBuffer]
                closestDF['buffer'] = bufferDist
                exactBuffer = exactBuffer.append(closestDF)
    exactBuffer = injectRoadIntoNames(exactBuffer)
    return(exactBuffer)


def createRdVarNames(curName):
    tempNames = [
        'buffer',
        'roadType',
        curName + 'len',
        curName + 'tr',
        curName + 'bl',
        curName + 'wn',
        curName + 'trsh',
        curName + 'blsh',
        curName + 'anwn',
        curName + 'percwn',
        curName + "maxwn",
        curName + 'cutoff',
        curName + 'cutoffLen',
        curName + 'cutBuilding',
        curName + 'cutTree'
    ]
    return(tempNames)

# denomarlize variables by adding road names.  Need to rewrite for the aadt categories
# INPUTS:
#    dataset (pandas dataframe) - exposures for a single buffer distance with normalized names
# OUTPUTS:
#    mjRds (pandas dataframe) - same exposures as input, but with denormalized variable names
def injectRoadIntoNames(dataset):
    rdNames = ['a0','a1','t0','t1','all']
    rdSets = []
    for index in range(len(rdNames)):
        curName = rdNames[index]
        tempNames = createRdVarNames(curName)
        tempRds = dataset[dataset['roadType']==index]
        try:
            tempRds.columns = tempNames
        except Exception as e:
            print("problem renaming: " + str(e))
            print(tempRds.head())
        tempRds = tempRds.drop(columns=['roadType'])
        if(tempRds.count()[0]>0):
            rdSets.append(tempRds)
    allRds = rdSets[0]
    for index in range(1,len(rdSets)):
        allRds = allRds.merge(rdSets[index])
    return(allRds)
   

# given a single csv file for a single maternal residence, load values into memory and prep
# for combining with other maternal residence csvs
# INPUTS:
#    filename (str) - relative filepath to csv
#    folder (str) - absolute filepath to folder containing csv
# OUTPUTS:
#    startBuffer (pandas dataframe) - csv values in dataframe format
def processSingleCSV(filename,folder):
    bufferDists = [50,100,200,300,500]
    curData = ps.read_csv(folder + filename)
    startBuffer = getBufferValues(curData,bufferDists[0])
    for buffer in bufferDists[1:]:
        try:
            tempVals = getBufferValues(curData,buffer)
            startBuffer = startBuffer.append(tempVals)
        except Exception as e:
            print(str(e))
    startBuffer['uniqueId_id'] = [filename[:-4] for x in range(startBuffer.count()[0])]
    return(startBuffer)

# for all csv files in a set (e.g. all maternal residences), load values into 
# memory and combine into a single csv
def calcAvgs(folder):
    csvFiles = getCSVFilesToAvg(folder)
    startBuffer = processSingleCSV(csvFiles[0],folder)
    index = 0
    for file in csvFiles[1:]:
        startBuffer = startBuffer.append(processSingleCSV(file,folder))
        if index % 250 == 0:
            print(index)
        index+=1
    startBuffer.to_csv(const.WIND_FOLDER + "bufferAvgs/bu_" + str(YEAR) + ".csv",index=False)


####################### MAIN FUNCTION ##################

if __name__ == '__main__':
    birthData = ps.read_csv(const.WIND_FOLDER + "Birth_Addresses_Wind/births_by_year/csvs/births_" + str(YEAR) + ".csv")
    # update filepaths for maternal residences rather than air monitors
    convertShpsToCSV(const.WIND_FOLDER + "shpFile/" + str(YEAR) + "/",birthData)
    calcAvgs(const.WIND_FOLDER + "bufferAvgs/" + str(YEAR) + "/")
