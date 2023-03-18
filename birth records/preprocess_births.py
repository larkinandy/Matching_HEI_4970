############### preprocess_births.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: preprocess birth records to prepare them for wind analysis
# Prepreocessing actions include :
# 1) creating a birthdate timestamp, 
# 2) creating a conception date timestamp,
# 3) partitioning birth records by year of birth,
# 4) exporting select variables into csv format


############### Setup: import libraries and define constants ##############
import arcpy
from datetime import datetime
import os
import pandas as ps
import gConst as const

# filepaths
PARENT_FOLDER = const.WIND_FOLDER
BIRTH_FOLDER = PARENT_FOLDER + "Birth_Addresses_Wind/"
BIRTH_INPUT = BIRTH_FOLDER + "original_shapefile/Births0716_Wind.shp"
BY_YEAR_FOLDER = BIRTH_FOLDER + "/births_by_year/"
ROAD_FOLER = const.ROAD_FOLDER 


# attribute fields for shapefile 
ID_COLUMN = "uniqueid"
AGE_COLUMN = "b_es_ges"
B_DAY,B_MONTH,B_YEAR = "bday","bmonth","byear"
B_DATE, C_DATE = "bdate", "cdate" # conception date
FIELDS_TO_DELETE = [
    'mothersres',
    'mothersr_1',
    'mothersr_2',
    'mothersr_3',
    'b_geocod',
    'RDS123_500',
    'RDS50k_500',
    'Const_Proj'
]

################ Helper functions called by the main function ##############

# create date type attribute fields to a shapefile
# INPUTS:
#    inTable (string) - full filepath of the shapefile
#    fieldNames (list of strings) - names of the attribute fields to create
def createAttributeFields(inTable,fieldNames):
    for fieldName in fieldNames:
        arcpy.AddField_management(inTable, fieldName, "DATE")

# calculate birth date timestamp from year, month, and day stored in seperate attribute fields
# INPUTS:
#    inTable (string) - full filepath of the shapefile
#    outputColumn (string) - name of the attribute field to store the birth date in
def calcBirthDate(inTable,outputColumn):
    # write an inline function to be called by arcpy.  Similar to a code block in the ArcPro GUI
    exp = '''def add_date(year,month,day):
        from datetime import datetime
        dateString = str(year) + "-" + str(month) + "-" + str(day)
        date = datetime.strptime(dateString, "%Y-%m-%d")
        return date
    '''
    call = "add_date(!" + B_YEAR + "!,!" + B_MONTH + "!,!" + B_DAY + "!)"
    arcpy.CalculateField_management(inTable,outputColumn,call,"PYTHON",exp)

# calculate conception date by subtracting estimated gestational age from the birth date
# INPUTS:
#    inTable (string) - full filepath of the shapefile
#    outputColumn (string) - name of the attribute field to store the conception date in
def calcConceptionDate(inTable,outputColumn):
    # write an inline function to be called by arcpy.  Similar to a code block in the ArcPro GUI
    exp = '''def add_date(bdate,gest_age):
        from datetime import datetime
        from datetime import timedelta
        b_timestamp = datetime.strptime(bdate, "%m/%d/%Y")
        diff = timedelta(weeks=gest_age)
        conception_date = b_timestamp - diff
        return conception_date
    '''
    call = "add_date(!" + B_DATE + "!,!" + AGE_COLUMN + "!)"
    arcpy.CalculateField_management(inTable,outputColumn,call,"PYTHON",exp)


# convert the births attrribute tables within a shapefile into a pandas dataframe
# INPUTS:
#    in_fc (str) - full filepath of the shapefile
# OUTPUTS:
#    fc_dataframe (pandas dataframe) - selected attributes from the 
#                                      shapefile attribute table
def arcgis_table_to_df(in_fc):
    print("converting to dataframe")
    
    # get names of all attributes
    OIDFieldName = arcpy.Describe(in_fc).OIDFieldName
    final_fields = [field.name for field in arcpy.ListFields(in_fc)]

    # extract data from attribute table into dataframe
    data = [row for row in arcpy.da.SearchCursor(in_fc,final_fields)]
    fc_dataframe = ps.DataFrame(data,columns=final_fields)
    fc_dataframe = fc_dataframe.set_index(OIDFieldName,drop=True)

    # reduce data to vars of interest
    valsList = ["uniqueid","b_long","b_lat","bdate","cdate"]
    fc_dataframe = fc_dataframe[valsList]
    return fc_dataframe


# export a birth shapefile into a csv file
# INPUTS:
#     outputFolder (str) - full filepath for the folder where the birth records are stored
#     year (str) - year of the birth records to convert to csv
def exportIdsToCSV(outputFolder,year):
    inShapefile = outputFolder + "births_r_" + str(year) + ".shp"
    outCSV = outputFolder + "/csvs/births_" +str(year) + ".csv"
    if not (os.path.exists(outCSV)):
        df = arcgis_table_to_df(inShapefile)
        df.to_csv(outCSV,index=False)


# partition the shapefile by birth year
# INPUTS:
#    inTable (string) - full filepath of the shapefile
#    years (int array) - years to partition by
#    outputFolder (string) - full filepath to the folder for the partitioned shapefiles
def partitionByYear(inTable,years,outputFolder):
    for year in years:
        # write an inline function to be called by arcpy.  Similar to a code block in the ArcPro GUI
        exp = '''"''' + B_YEAR + '''"=''' + str(year)
        tempFilename = outputFolder + "births_" + str(year) + ".shp"
        outputFilename = outputFolder + "births_r_" + str(year) + ".shp"
        birthsInYear = arcpy.SelectLayerByAttribute_management(
            inTable, 'NEW_SELECTION', exp
        )
        roadFile = ROAD_FOLER + str(year) + "/" + "merged" + str(year) + ".shp"
        arcpy.CopyFeatures_management(birthsInYear,tempFilename)     
        birthsInYear2 = arcpy.SelectLayerByLocation_management(
            tempFilename, 'WITHIN_A_DISTANCE', roadFile, '500 Meters'
        )
        arcpy.CopyFeatures_management(birthsInYear2,outputFilename)
        exportIdsToCSV(outputFolder,year)
        



################# Main function ######################

if __name__ == '__main__':
    createAttributeFields(BIRTH_INPUT,[B_DATE,C_DATE])
    calcBirthDate(BIRTH_INPUT,B_DATE)
    calcConceptionDate(BIRTH_INPUT,C_DATE)
    arcpy.DeleteField_management(BIRTH_INPUT,FIELDS_TO_DELETE)
    partitionByYear(BIRTH_INPUT,list(range(2007,2016)),BY_YEAR_FOLDER)
