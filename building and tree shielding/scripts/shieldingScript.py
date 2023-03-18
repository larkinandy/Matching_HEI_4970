############### shieldingScriptMaternalResidence.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: calculate road, wind, and shielding exposures.
# THIS IS A WORK IN PROGRESS, NOT READY FOR DEPLOYMENT

# Steps include:
# 1) loading hourly wind estimates from storage
# 2) reformatting data for matrix algebra operations on GPUs
# 3) push data to GPU memory
# 4) calculate number of days downwind for each radial degree using GPUS
# 5) push wind matrix back to RAM 
# 6) save wind matrix to storage



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


YEAR = 2009

# needed to create geometry points for each maternal residence lat/lon
GCS = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]" 

# read c_date, b_date from birth records file
BIRTH_META = const.WIND_FOLDER + "Birth_Addresses_Wind/births_by_year/csvs/births_" + str(YEAR) + ".csv"

SMALLEST_BUFFER, BIGGEST_BUFFER = 10, 500 # analysis is restricted to 0.5km
PARENT_FOLDER = const.WIND_FOLDER

# given point of maternal residence, create buffers from 10 to 500 m at 10m increments
# INPUTS:
#    geometryPoint (arcpy geometry object) - point cenetered at the maternal residence
#    paths (dictionary) - set of input and output filepaths
def createBuffers(geometryPoint,paths):
    if not os.path.exists(paths['bufferFolder']):
        os.mkdir(paths['bufferFolder'])
    curBuffer = SMALLEST_BUFFER
    while (curBuffer <= BIGGEST_BUFFER):
        outputShapefilepath = paths['bufferFolder'] + "/b" + str(curBuffer) + paths['uniqueid'] + ".shp"
        if not(os.path.exists(outputShapefilepath)):
            arcpy.Buffer_analysis(
                geometryPoint, outputShapefilepath, str(curBuffer) + " Meters", "FULL", "ROUND", method="GEODESIC"
            )
        curBuffer +=10


# given point of maternal residence, create buffers from 10 to 500 m at 10m increments
# INPUTS:
#    geometryPoint (arcpy geometry object) - point cenetered at the maternal residence
#    paths (dictionary) - set of input and output filepaths
def deleteBuffers(paths):
    curBuffer = SMALLEST_BUFFER
    while (curBuffer <= BIGGEST_BUFFER):
        outputShapefilepath = paths['bufferFolder'] + "/b" + str(curBuffer) + paths['uniqueid'] + ".shp"
        if (os.path.exists(outputShapefilepath)):
            arcpy.Delete_management(outputShapefilepath)
        curBuffer +=10

# given list of input and output filepaths, create clips of the wind shapefile at increments of 10m
# INPUTS:
#     paths (dictionary) - set of input and output filepaths
def createWindBuffers(paths):
    if not os.path.exists(paths['windBufferFolder']):
        os.mkdir(paths['windBufferFolder'])
    curBuffer = SMALLEST_BUFFER
    while (curBuffer <= BIGGEST_BUFFER):
        outputShapefilepath = paths['windBufferFolder'] + "/bw" + str(curBuffer) + paths['uniqueid'] + ".shp"
        bufferShapefile = paths['bufferFolder'] + "/b" + str(curBuffer) + paths['uniqueid'] + ".shp"
        arcpy.Clip_analysis(paths['windShp'], bufferShapefile, outputShapefilepath)
        curBuffer +=10


# given list of input and output filepaths, create clips of the wind shapefile at increments of 10m
# INPUTS:
#     paths (dictionary) - set of input and output filepaths
def deleteWindBuffers(paths):
    curBuffer = SMALLEST_BUFFER
    while (curBuffer <= BIGGEST_BUFFER):
        outputShapefilepath = paths['windBufferFolder'] + "/bw" + str(curBuffer) + paths['uniqueid'] + ".shp"
        if(os.path.exists(outputShapefilepath)):
            arcpy.Delete_management(outputShapefilepath)
        
        curBuffer +=10


# in case of missing buffer values, find the nearest buffer with a valid value
# INPUTS:
#    num1 (int) - buffer distance where the value is missing
#    numList (int list) - list of distances with valid values
# OUTPUT:
#    closest distance with a valid value
def findNearestIndex(num1,numList):
    diffs = np.abs(np.subtract(numList,num1))
    return(np.argmin(diffs))

# fill in missing values for analyses where the tree rasters have missing values
# INPUTS:
#    intable (filepath) - absolute filepath of shapefile with an attribute table
#    bufferDist (int) - buffer distance that is missing a value
# OUTPUTS: 
#    pandas dataframe representation of the intable, updated with missing vals filled
def fillMisingVals(intable,bufferDist):
    rawData = table_to_data_frame(intable)
    coveredIndexes = list(rawData['FID_'])
    vals = list(rawData['MEAN'])

    # QA check
    if(len(vals)<360):

        # for each radial degree, if the value is missing then fill with the nearest bufffer distance of the 
        # same radial degree with a valid value
        for index in range(360):
            if index not in coveredIndexes:
                closestIndex = findNearestIndex(index,coveredIndexes)
                coveredIndexes.append(index)
                vals.append(vals[closestIndex])

    # combine updated lists into a single dataframe
    updatedDF = ps.DataFrame({
        'FID':coveredIndexes,
        'trSh' :vals,
        'bufferDist':[bufferDist for a in range(len(vals))]
    })
    return(updatedDF)

# convert arcgis table of building shielding estimates into pandas dataframe
# INPUTS:
#    inTable (str) - full filepath to the arcgis table
#    bufferDist (int) - buffer distance estimates in the inTable
# OUTPUTS:
#    updatedDF (pandas dataframe) - in Table in pandas format
def buildingToDF(inTable,bufferDist):
    rawData =  table_to_data_frame(inTable)
    coveredIndexes = list(rawData['FID_'])
    vals = list(rawData['PERCENTAGE'])
    updatedDF = ps.DataFrame({
        'FID':coveredIndexes,
        'bldgSh' :vals,
        'bufferDist':[bufferDist for a in range(len(vals))]
    })
    return(updatedDF)

# convert arcgis table into a pandas df with an object ID index
# INPUTFS:
#    in_table (str) - full filepath to the table to convert
# OUTPUTS:
#    pandas dataframe version of in_table
def table_to_data_frame(in_table):
    OIDFieldName = arcpy.Describe(in_table).OIDFieldName
    final_fields = [field.name for field in arcpy.ListFields(in_table)]
    data = [row for row in arcpy.da.SearchCursor(in_table, final_fields)]
    fc_dataframe = ps.DataFrame(data, columns=final_fields)
    fc_dataframe = fc_dataframe.set_index(OIDFieldName, drop=True)
    return fc_dataframe

# calculate percent tree between roads and maternal residence.  Repeat every 10m
# to apply to road networks at 10m resolution
# INPUTS:
#    bufferDist (int) - buffer distance to calculate percent trees
#    paths (dictionary) - list of input and output filepaths
# OUTPUTS:
#    dataframe containing estimated percent tree cover for each radial degree for the
#    given buffer distance
def calcZonalIntersectTree(bufferDist,paths):
    try:
        windShapefile = paths['windBufferFolder'] + "/bw" + str(bufferDist) + paths['uniqueid'] + ".shp"
        zonalFile =   paths['tmpTreeFolder'] + "/z" + str(bufferDist) + "_" + paths['uniqueid'] + ".dbf" # temporary file.  Deleted later

        # calcualte zonal statistics, fill in mising vals, and return as dataframe
        if not(os.path.exists(zonalFile)):
            arcpy.sa.ZonalStatisticsAsTable(
                windShapefile, "FID", paths['resampled_tree'], zonalFile, "DATA", "MEAN"
            )
        df = fillMisingVals(zonalFile,bufferDist)
        arcpy.Delete_management(zonalFile)
        return(df)
    except Exception as e:
        print("failed to calculate zonal statistics for trees")
        print(str(e))
        print(str(bufferDist))

# calculate percent building between roads and maternal residence.  Repeat every 10m
# to apply to road networks at 10m resolution
# INPUTS:
#    bufferDist (int) - buffer distance to calculate percent building
#    paths (dictionary) - list of input and output filepaths
# OUTPUTS:
#    dataframe containing estimates percent building cover for each radial degree for the given 
#    buffer distances
def calcZonalIntersectBuilding(bufferDist,paths):
    try:
        windShapefile = paths['windBufferFolder'] + "/bw" + str(bufferDist) + paths['uniqueid'] + ".shp"
        outputZonalFile = paths['tmpBuildingFolder'] + "/i" + str(bufferDist) + "_" + paths['uniqueid'] + ".dbf" # temporary file.  Deleted later

        # calculate zonal statistics, return as dataframe
        if not(os.path.exists(outputZonalFile)):
            arcpy.TabulateIntersection_analysis(
                windShapefile, "FID", paths['clipped_building'], outputZonalFile
            )
        df = buildingToDF(outputZonalFile,bufferDist)
        arcpy.Delete_management(outputZonalFile)
        return(df)
    except Exception as e:
        print("failed to calcualte zonal statistics for buildings")
        print(str(e))
        print(str(bufferDist))

# given one dataframe for each buffer distance, combine all dataframes into a single dataframe
# INPUTS:
#   dfArr (list of dataframes) - set of multiple dataframes, one for each buffer distance
# OUTPUTS:
#   mergedSet (pandas dataframe) - single dataframe containing all data from dfArr
def mergeSet(dfArr):
    mergedSet = dfArr[0]
    for curDF in dfArr[1:]:
        mergedSet = mergedSet.append(curDF)
    return(mergedSet)

# given shielding estimates for building and trees, combine them into a single dataset
# INPUTS:
#    shield1 (pandas dataframe) - pandas dataframe containing tree shielding estimates
#    shield2 (pandas dataframe) - contains building shielding estimates
# OUTPUTS:
#    mergedSet (pandas dataframe) - combined data from shield1 and shield2
def mergeShields(shield1,shield2,paths):
    mergedSet = ps.merge(
        shield1,shield2,how='outer',on=["FID","bufferDist"]
    )
    mergedSet = mergedSet.fillna(0)
    mergedSet['joinStr'] = mergedSet['FID'].astype(str).str.zfill(3) + mergedSet['bufferDist'].astype(str).str.zfill(3)
    outputFilename = paths['bufferFolder'] + "/shieldMeas" + paths['uniqueid'] + ".csv"
    mergedSet.to_csv(outputFilename,index=False)
    return(mergedSet)

# restrict road networks to within 500m of the maternal residence.  Greatly speeds up downstream computations
# INPUTS:
#    paths (dictionary) - list of filepaths
def clipRds(paths):
    if not(os.path.exists(paths['tmpRdFolder'])):
        os.mkdir(paths['tmpRdFolder'])
    tmpFiles = []

    # for each road classification, restrict road networks to wihtin 500m of residence and create temporary files of restricted data
    for index in range(len(paths['rdArray'])):
        tempFile = paths['tmpRdFolder'] + "/a" + str(index) + paths['uniqueid'] + ".shp"
        curRdInFile = paths['rdArray'][index]
        curRdBufferFile = paths['bufferFolder'] + "/b" + str(BIGGEST_BUFFER) + paths['uniqueid'] + ".shp"
        arcpy.Clip_analysis(curRdInFile, curRdBufferFile, tempFile)
        arcpy.AddField_management(tempFile, "rdType", "SHORT")
        arcpy.CalculateField_management(tempFile, "rdType", str(index),"PYTHON")
        tmpFiles.append(tempFile)
    
    # combine temporary road files into a single file with an additional attribute field to desginate road type
    arcpy.Merge_management(tmpFiles, paths['tmpRdFolder'] + "/rdsCombined" + paths['uniqueid'] + ".shp")
    arcpy.AddField_management(paths['tmpRdFolder'] + "/rdsCombined" + paths['uniqueid'] + ".shp", "wndFID", "SHORT")
    
    # cleanup temporary files created during this function
    for filepath in tmpFiles:
        arcpy.Delete_management(filepath)

# partition the road network around a maternal residence based on each of the 360 radial degrees
# INPUTS:
#    paths (dictionary) - list of filepaths
def createRoadIntersects(paths):

    # make a new feature called stateslyr
    arcpy.MakeFeatureLayer_management (paths['windBufferFolder'] + "/bw500" + paths['uniqueid'] + ".shp", "stateslyr")
    outputArray = []

    # for each radial degree, create an in memory copy of the road network clipped to the radial coverage.  
    # at the end, merge in memory files and write to disk
    for i in range(360):
        outputRds = "in_memory" + "/rds" + str(i) 
        outputArray.append(outputRds)
        FID_String = ' "FID" = ' + str(i) + ' '
        arcpy.SelectLayerByAttribute_management("stateslyr", "NEW_SELECTION", FID_String)
        arcpy.Clip_analysis(paths['tmpRdFolder'] + "/rdsCombined" + paths['uniqueid'] + ".shp", "stateslyr", outputRds)
        arcpy.CalculateField_management(outputRds, "wndFID", str(i),"PYTHON")
    arcpy.Merge_management(outputArray,paths['tmpRdFolder'] + "/rdsSplit" + paths['uniqueid'] + ".shp")
    
    # cleanup, delete in memory files
    for curFile in outputArray:
        arcpy.Delete_management(curFile)
    
# get the distances of road segments within 500m of residences.  Only distances with roads need to be processsed, all other 
# buffer sizes can be ignored to improve computational speed
# INPUTS:
#    paths (dictionary) - list of filepaths
# OUTPUTS:
#    validDists (list of ints) - buffer distances that need to be processed
def getBufferDistsToProcess(paths):
    validDists = []
    rdsShp = paths['tmpRdFolder'] + "/rdsSplit" + paths['uniqueid'] + ".shp"
    cursor = arcpy.SearchCursor(rdsShp)
    row = cursor.next()

    # search through all road segments.  Round distance from road to residence to nearest 10m, and return set of 
    # road segments distances
    while row:
        newDist = int(int(row.getValue("NEAR_DIST"))/10)*10
        if(newDist <= BIGGEST_BUFFER and newDist not in validDists):
            validDists.append(int(int(row.getValue("NEAR_DIST"))/10)*10)
        row = cursor.next()
    validDists = list(set(validDists))
    del cursor, row
    return(validDists)

# create a new attribute field in two shapefiles that will be used to link the two shapefiles together
# INPUTS:
#    shapefile (str) - filepath to a shapefile representing the maternal residence
#    paths (dictionary) - list of other filepaths used in the analysis
def linkFiles(shapefile,paths):
    arcpy.Near_analysis(paths['tmpRdFolder'] + "/rdsSplit" + paths['uniqueid'] + ".shp", shapefile, "", "NO_LOCATION", "NO_ANGLE", "GEODESIC")
    arcpy.AddField_management(paths['tmpRdFolder'] + "/rdsSplit" + paths['uniqueid'] + ".shp", "rndDist", "TEXT")
    arcpy.CalculateField_management(paths['tmpRdFolder'] + "/rdsSplit" + paths['uniqueid'] + ".shp", "rndDist", "str(!wndFID!).zfill(3) + str(int(!NEAR_DIST!/10)*10).zfill(3)", "PYTHON")

# calculate tree and building shielding for all buffer distances
# INPUTS:
#    paths (dictionary) - list of filepaths
# outputs:
#    pandas dataframe containing all tree and shielding estimates 
def calcShields(paths):
    
    # get buffer distances with a road segment
    bufferDists = getBufferDistsToProcess(paths)
    arr = []
    bufferDists.sort()
    
    # for all buffer dsitances with a road sgement, calculate tree shielding
    for buf in bufferDists:
        arr.append(calcZonalIntersectTree(buf,paths))
    combinedTree = mergeSet(arr)
    arr = []

    # for all buffer distances with a road segment, calcualte building shielding
    for buf in bufferDists:
        arr.append(calcZonalIntersectBuilding(buf,paths))

    # combine building and tree shielding estimates and return as pandas df
    combinedBuilding = mergeSet(arr)
    combinedAll = mergeShields(combinedTree,combinedBuilding,paths)
    return(combinedAll)

# join shield, wind, and road len metrics into a single shapefile
# INPUTS:
#     paths (pandas dataframe) - list of filepaths
def joinFiles(paths):

    # load shielding metrics into memory
    tempname = 'a' #+ paths['uniqueid']
    arcpy.conversion.TableToTable(paths['metricsCSV'], "in_memory", tempname)
    # add a new variable 'join_str2' used to combine shielding metrics with other metrics
    arcpy.AddField_management("in_memory/" + tempname, "join_str2", "TEXT",'#','#',6)
    arcpy.management.CalculateField(
        "in_memory/" + tempname, 
        "join_str2", 
        "str(!FID!).zfill(3) + str(!bufferDist!).zfill(3)", 
        "PYTHON3", 
        '', 
        "TEXT"
    )
    # join road and shielding metrics into new shapefile
    joinTable = arcpy.AddJoin_management(paths['rdsShp'], "rndDist", "in_memory/" + tempname, "join_str2", "KEEP_ALL")
    #arcpy.conversion.TableToTable(joinTable, paths['bufferFolder'],paths['tempMetrics'])
    # joint road,shielding,and wind metrics into a new shapefile
    joinTable2 = arcpy.AddJoin_management(joinTable,"wndFID",paths['windShp'],"FID","KEEP_ALL")
    arcpy.CopyFeatures_management(joinTable2, paths['shieldShp'])

# define input, temporary, and output filepaths used throughout this entire script
# INPUTS:
#    birth record (pandas dataframe) - contains only 1 row, one unique birth record
# OUTPUTS:
#    dictionary of filepaths
def defineFilepaths(birthRecord):
    paths = {}
    TREE_FOLDER = PARENT_FOLDER + "tr/"
    BUILDING_GEODATABASE = PARENT_FOLDER + "AnnualFootprints/MyProject6.gdb"

    
    ROADS_FOLDER = PARENT_FOLDER + "roads/" + str(YEAR) + "/"
    
    paths['uniqueid'] = birthRecord['uniqueid']
    paths['bufferFolder'] = PARENT_FOLDER + "buffers"
    paths['windBufferFolder'] = paths['bufferFolder'] + "/wind"
    paths['tmpTreeFolder'] = paths['bufferFolder'] + "/tmpTree"
    paths['tmpBuildingFolder'] = paths['bufferFolder'] + "/tmpBuilding"
    paths['tmpRdFolder'] = paths['bufferFolder'] + "/tmpRds"
    paths['tree_raster'] = TREE_FOLDER + "texas_tr2015.tif" if (YEAR> 2012) else TREE_FOLDER + "texas_tr2010.tif"
    paths['clipped_tree'] = TREE_FOLDER + paths['uniqueid'] + "_trclip.tif"
    paths['building_shapefile'] = BUILDING_GEODATABASE + "/footprints_" + str(YEAR)
    paths['resampled_tree'] = paths['tmpTreeFolder'] + "/t_" + paths['uniqueid'].lower() + "r.tif" 
    paths['clipped_building'] = paths['tmpBuildingFolder'] + "/bu" + paths['uniqueid'] + ".shp"
    paths['bufferExtent'] = paths['bufferFolder'] + "/b" + str(BIGGEST_BUFFER) + paths['uniqueid'] + ".shp"
    paths['rdsShp'] = paths['tmpRdFolder'] + "/rdsSplit" + birthRecord['uniqueid'] + ".shp"
    paths['metricsCSV'] = paths['bufferFolder'] + "/shieldMeas" + birthRecord['uniqueid'] + ".csv"
    paths['tempMetrics'] = birthRecord['uniqueid'] + "a.csv"
    paths['shieldShp'] = PARENT_FOLDER + "shpFile/" + str(YEAR) + "/" + birthRecord['uniqueid'] + ".shp"
    paths['rdArray'] = [
        ROADS_FOLDER + "/aadt0.shp",
        ROADS_FOLDER + "/aadt1.shp",
        ROADS_FOLDER + "/taadt0.shp",
        ROADS_FOLDER + "/taadt1.shp"    
    ]
    paths['windShp'] = birthRecord['windShp']

    return(paths)

# given a birth record with lat/lon coordinates, create an arcgis geometry object
# IMPORTANT: MATERNAL RESIDENCE COORDINATES MUST BE IN GCS WGS84 !!!
# INPUTS:
#    birthRecord (pandas dataframe) - contains one row, one unique birth record
# OUTPUTS:
#    geometry object centered on the maternal residence coordinates
def createGeometryPoints(birthRecord):
    try:
        tempPoint = arcpy.Point()
        tempPoint.X = float(birthRecord['b_long'])
        tempPoint.Y = float(birthRecord['b_lat'])
        tempMonitorShp = arcpy.PointGeometry(tempPoint,GCS)
        return(tempMonitorShp)
    except Exception as e:
        print("couldn't create point: " + str(e))
    return(None)

# create buffers from 10m to 500m centered on maternal residence
# INPUTS:
#    geometryPoint (arcpy geometry object) - centered on maternal residence
#    paths (dictionary) - list of filepaths
def createAllBuffers(geometryPoint,paths):
    createBuffers(geometryPoint,paths)
    createWindBuffers(paths)

# create buffers from 10m to 500m centered on maternal residence
# INPUTS:
#    geometryPoint (arcpy geometry object) - centered on maternal residence
#    paths (dictionary) - list of filepaths
def deleteAllBuffers(paths):
    deleteBuffers(paths)
    deleteWindBuffers(paths)

# prepare building and tree dtasets for shielding analysis.  Operations include clipping datasets
# to 500m and resamling the tree raster 
# INPUTS:
#    paths (dict) - list of filepaths
def preprocessShields(paths):
    outExtractByMask = arcpy.sa.ExtractByMask(paths['tree_raster'], paths['bufferExtent'])
    outExtractByMask.save(paths['clipped_tree'])
    arcpy.Resample_management(in_raster=paths['clipped_tree'], out_raster=paths['resampled_tree'], cell_size="0.000026949456 0.000026949456", resampling_type="BILINEAR")
    arcpy.Clip_analysis(paths['building_shapefile'], paths['bufferExtent'], paths['clipped_building'])

# peform all steps to calculate wind, shielding, and road metrics for one single maternal residence
# INPUTS:
#    birthRecord (pandas dataframe) - single row, unique birth record
def processSingleResidence(birthRecord):
    #print(birthRecord)
    # create temporary and output filepaths that incorporate the birth reqcord unique identifier
    paths = defineFilepaths(birthRecord)
    if(os.path.exists(paths['shieldShp'])):
        print("estimates already derived for record " + paths['uniqueid'])
        #try:
        #    deleteAllBuffers(paths)
        #except Exception as e:
        #    print("couldn't delete buffers: " + str(e))
        return
    if not(os.path.exists(paths['windShp'])):
        print("no wind estimates for record " + paths['uniqueid'])
        return
    # create an arcgis point geometry object centered on the maternal residence
    geometryPoint = createGeometryPoints(birthRecord)
    if(geometryPoint == None): return
    print("created geometry")

    # create buffers from 10m to 500m
    try:
        createAllBuffers(geometryPoint,paths)    
    except Exception as e:
        print("couldn't create buffers: " + str(e))
        return

    # prepare tree and building datasets for shielding metrics
    try:
        preprocessShields(paths)
        print("finished prep")
    except Exception as e:
        print("couldn't preprocess shields: " + str(e))
        return
    try:
        clipRds(paths)
        createRoadIntersects(paths)
        linkFiles(geometryPoint,paths)
        calcShields(paths)
        joinFiles(paths)
    except Exception as e:
        print("culdn't derive road weighted shield and wind: " + str(e))
        return
    try:
        deleteAllBuffers(paths)
    except Exception as e:
        print("couldn't delete buffers: " + str(e))
    print("completed birth record %s" %(birthRecord['uniqueid']))

# load birth records and use unique id to get the filepath for the maternal
# residence rose wind shapefile
# OUTPUTS:
#    vbariables from birth records in pandas dtaframe format
def getBirthMeta():
    WIND_FOLDER = PARENT_FOLDER + "outputShapefiles/" + str(YEAR) + "/"
    read = True
    while(read):
        try:
            birthMeta = ps.read_csv(BIRTH_META)
            read = False
        except Exception as e:
            time.sleep(1)
    nRows = birthMeta.count()[0]
    windShp = ['a' for x in range(nRows)]
    for rowNum in range(nRows):
        curRow = birthMeta.iloc[rowNum]
        curId = curRow['uniqueid']
        # get the filepath to the maternal residence rosewind shapefile
        # using the unique identifiers
        shpPath = WIND_FOLDER + curId + ".shp"
        windShp[rowNum] = shpPath
    birthMeta['windShp'] = windShp
    print(birthMeta.head())
    return(birthMeta) 

# transform data from table to tuple format to facilitate parallel processing
# INPUTS:
#    birthRecords (pandas dataframe) - data to transform
# OUTPUTS:
#    list of tuples, where threads can independently process data in each tuple
def prepForParallel(birthRecords):
    nRows = birthRecords.count()[0]
    parallelList = []
    for rowNum in range(nRows):
        curRow = birthRecords.iloc[rowNum]
        parallelList.append(curRow)
    return(parallelList)


####################### MAIN FUNCTION ##################
if __name__ == '__main__':
    
    # load birth records and convert to tuple for parallel processing
    birthMeta = getBirthMeta()
    result = prepForParallel(birthMeta)
    random.shuffle(result)

    pool = Pool(processes=64)
    res = pool.map_async(processSingleResidence,result)
    res.get()
   
    
    