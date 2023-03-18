############### createWindMatrix.py #############
# Author: Andrew Larkin
# Developed for HEI Transit Study
# Summary: using hourly wind data, create daily estimates of hours upwind for each of 360 radial degrees
#          around each maternal residence in a yearly subset of the Texas cohort
# Steps to create wind estimates include :
# 1) loading hourly wind estimates from storage
# 2) reformatting data for matrix algebra operations on GPUs
# 3) push data to GPU memory
# 4) calculate number of days downwind for each radial degree using GPUS
# 5) push wind matrix back to RAM 
# 6) save wind matrix to storage


# important information:
# tested on RTX 3090 and RTX Titan GPUs
# 24GB GPU ram is necesssary for a batch size of 1000.  Decrease batch size for GPUs with less RAM
# CUDA version 11.3.  Newer version of CUDA are not comptabile with numba at this time (Fall 2021).


############### Setup: import libraries and define constants ##############
from __future__ import division
from numba import cuda, float32
import numpy as np
import math
import numba
import pandas as ps
import time
import os
import gConst as const

BIRTH_FOLDER = const.WIND_FOLDER + "Birth_Addresses_Wind/births_by_year/"
YEAR = 2009


# constants necessary for GPU oeprations.  MUST BE SET/HARD-CODED BEFORE RUNTIME!!!!
N_ANGLES = 360
N_BIRTHS = 624 # input data is preprocessed into batch sizes of 1000
thread_x_dim = math.ceil(360/128)*128 # number of threads per streaming multiprocessor
thread_dim = (thread_x_dim,1) # cuda requires 2d tuple
block_x_dim = 24 # months in 2 years, one of the dimensions of the output wind matrix
block_y_dim = 31 # max days in month, one ofthe dimensions of the output wind matrix
block_z_dim = N_BIRTHS # number of 
block_dim = (block_x_dim,block_y_dim,block_z_dim)
TOTAL_THREADS = thread_x_dim*block_x_dim*block_y_dim*block_z_dim



################ Helper functions called by the main function ##############

# create a uniqueId->index map.  Called in combination with the python map function
# INPUTS:
#    idToMap (str) - unique id too map
# OUTPUTS:
#    index of the unique id (int)
def mapIds(idToMap):
    return masterList.index(idToMap)

# read hourly wind estimates extracted for all residences in the birth year
# INPUTS:
#    year (int) - birth year
# OUTPUTS:
#    windTuples (array) - arrray of tuples.  Each tuple contains sorted 
#                         hourly wind estimates and unique ids
def readRawData(year):
    windTuples = [] # each tuple contains wind estimates and unique ids for 1000 births
    windData= ps.read_csv(const.WIND_FOLDER + "/comb_" + str(year) + ".csv")
    print("finished loading wind data")
    uniqueIds = list(set(windData['master_id'])) # all unique ids in the curr batch
    numIds = len(uniqueIds) 
    curIndex = 0

    # for each 1000 unique ids, get wind records for those ids, sort the ids and wind records
    # in the same order, and create a tuple.  Adjust the month variable to be relative to the first 
    # month of coverage (january year before birth)
    while(curIndex < numIds):
        subsetids = uniqueIds[curIndex: min(curIndex+1000,numIds)] # get next 1000 unique ids
        subsetData = windData[windData['master_id'].isin(subsetids)] # get records for this subset of ids
        subsetData.sort_values(by=['master_id'],inplace=True) # sort wind records by unique ids
        subsetids.sort() # sort unique ids
        
        # for records in the year of the birth, add 12 to the months value.  This makes the month variable
        # relative to the first month of coverage (january year before birth)
        tempVals = (subsetData['year'] - year +1)*12 # 
        subsetData['month'] = subsetData['month'] + tempVals
        
        windTuples.append((subsetData,subsetids))
        curIndex+=1000 # go to next batch of 1000 unique ids
        print(curIndex)
    print("completed partitioning wind data")
    return(windTuples)



################### precombiled GPU operations #####################
# cuda function.  Must be compiled before runtime (@cuda.jit decorator)



# rearrange data so hourly records in the matrix are organized by:
#      matrix[hourIndex][dayIndex][monthIndex][uniqueIdIndex]
# INPUTS:
#    idList (numpy array) - indeces for unique ids
#    hour (numpy arrray) - hour of the wind direction measures
#    day (numpy array) - day of the wind direction measures
#    month (numpy array) - month of the wind direction measures
#    wind_dir (numpy array) - wind direction measures
#    outData (numpy array) - stores reorganized wind direction measures
#    maxThread (int) - number of wind direction measures
@cuda.jit
def rearrangeData(idList,hour,day,month,wnd_dir,outData,maxThread):
    tid = cuda.blockIdx.x*cuda.blockDim.x + cuda.threadIdx.x
    if(tid < maxThread):
        outData[hour[tid],day[tid]-1,month[tid]-1,idList[tid]] = wnd_dir[tid]
    

# calculate number of hours each radial degree is upwind of the maternal residence
# results are stored in place in the output dataset
# INPUTS:
#    inputData (numpy array) - hourly wnd direction matrix
#    outputData (numpy array) - where results are stored in place
@cuda.jit
def calcDailyDownwind(inputData,outputData):
    
    
    # inline function.  Must be within compiled function so it can be sent to the 
    # GPU driver and implemented by each streaming multiprocessor at runtime
    # test if a radial degree is within the +/- 15 degree angular window for being downwind
    # INPUTS:
    #    val1 (int) - radial degree
    #    val2 (int) - wind direction
    # OUTPUTS:
    #     True if radial degree is within 15 degrees of wind direction. False otherwise.
    def isDownwind(val1,val2):
        diff1 = val1 - val2 if val1-val2 > 0 else (val1-val2 + 360)
        diff2 = val2 - val1 if val2-val1 >0 else (val2 - val1 + 360)
        minDiff = min(diff1,diff2)
        if(minDiff < 15):
            return 1
        return 0
    
    # each of the 360 degrees is assigned a unique thread id along the x axis of the cuda block
    angle = cuda.threadIdx.x
    
    # cuda has thousands of threads per sm along the x axis.  No need to go beyond 360
    if(angle >= 360):
        return

    # important that these dimensions/indeces match those assigned in the funtion 'rearrange data'
    day = cuda.blockIdx.y
    month = cuda.blockIdx.x
    station = cuda.blockIdx.z

    # for each hour in 24 hours of a day, test if the current radial degree is downwind of the wind direction
    # add to the hours downwind if yes, do not change hoursDownwind otherwise
    hoursDownwind = 0
    for hourIndex in range(24):
        wndAngle = inputData[hourIndex,day,month,station]
        if(wndAngle>=0):
            hoursDownwind += isDownwind(angle,wndAngle)
        else:
            hoursDownwind = -999
    
    # WARNING: IT IS ESSENTIAL THAT EACH THREAD HAS GUARANTEED EXCLUSIVE ACCESS TO IT's INDEX IN THE 
    # outputData matrix.  One thread, one index exactly.

    # update the output dataset in place
    outputData[station,month,day,angle] = hoursDownwind



################# Main function ######################
windTuples = readRawData(YEAR)
index = 0


# for each batch of 1000 maternal residences
# rearrange data into matrix organized by day, month, year, and unique id
# calculate daily number of hours downwind for all 360 degrees 
# save data as npy file
for tupleVal in windTuples:
    outputFile = const.WIND_FOLDER + 'windPartitions/' + str(YEAR) + '/w_ ' + str(index) + '.npy'
    
    # no need to process data if output file was already created
    if not(os.path.exists(outputFile)):
        print("starting index %i" %(index))

        windData, masterList = tupleVal[0], tupleVal[1]

        # get index for each unique id
        inputIds = list(map(mapIds,windData['master_id']))

        # transfer data from RAM to GPU memory.  Part of this process
        # includes creating variables on the GPU
        inputList = cuda.to_device(np.array(inputIds))
        cuda.synchronize()
        inputHour = cuda.to_device(np.array(windData['hour']))
        cuda.synchronize()
        inputDay = cuda.to_device(np.array(windData['day']))
        cuda.synchronize()
        inputMonth = cuda.to_device(np.array(windData['month']))
        cuda.synchronize()
        inputDir = cuda.to_device(np.array(windData['wnd_dir']))
        cuda.synchronize()        
        testInput = np.full((24,31,24,N_BIRTHS),-999,np.int16)
        outData = cuda.to_device(testInput)
        cuda.synchronize()
        nThreads = len(windData['hour'])
        cuda.synchronize()
        nBlocks = math.ceil(nThreads/thread_x_dim)
        cuda.synchronize()

        # rearrange wind direction so its organized in a matrix by day,month,year,unique id
        rearrangeData[(nBlocks,1),(thread_x_dim,1)](inputList,inputHour,inputDay,inputMonth,inputDir,outData,nThreads)
        cuda.synchronize()
        print("finished rearranging data")

        blankFile = np.full((N_BIRTHS,24,31,360),-999,np.float) # blank output matrix.  CUDA writes in place
        cuda.synchronize()

        station_mem = cuda.to_device(blankFile) # version of blank file to store on the GPU
        cuda.synchronize()

        # calcualte number or hours downwind for each of 360 degrees, for each day, each maternal residence
        calcDailyDownwind[block_dim,(thread_dim)](outData,station_mem)
        print("calcualted daily averages")
        cuda.synchronize()
     
        # transfer results from GPU to RAM
        resStation = station_mem.copy_to_host()

        # convert results to pandas dataframe and save 
        with open(outputFile, 'wb') as f:
            np.save(f, resStation)
        masterDF = ps.DataFrame({
            'masterIds':masterList
        })
        masterDF.to_csv(const.WIND_FOLDER + "masterIds/" + str(YEAR) + "/id_" + str(index) + ".csv",index=False)
    index+=1    
    
