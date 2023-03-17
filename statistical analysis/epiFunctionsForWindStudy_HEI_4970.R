### epiFunctionsForWindStudy_HEI_4970.R ###
# Author: Andrew Larkin
# Date created: Feb 1, 2023
# Summary: This R file contains functions for performing the wind based epidemiological analyses for HEI contract 4970, RFA 18-1/20-4-2.  
# Operations performed by these functions include running linear and logistic regression models, calculating tertiles and indicators 
# from continous variables, and restricting dataset by a given property such as distance to road for sensitivity models.  

# for details about cohort records and vital statistics, see "Data Dictionary - HEI Birth Data - 1996 to 2016 V1.docx"
# for details about deriving the wind-based exposure metrics, see "Deriving Wind Exposure Metrics_HEI_4970.docx"
# for details about matching upwinwd/downwind neighbords, see "Matching by Wind_HEI_4970.dox"
# the scripts which implements these functions for the epidemiologial analysis are titled:
#         - "epiAnalysisWindBirthWeight_HEI_4970.R"
#         - "epiAnalysisWindLowBirthWeight_HEI_4970.R"
#         - "epiAnalysisWindPretermBirth_HEI_4970.R"
#         - "epiAnalysisWindVeryPretermBirth_HEI_4970.R"


######### LOAD LIBRARIES AND SECRETS #############
library(rstudioapi)
library(survival)
setwd("C:/users/larki_9x8fs5f/Desktop/TexasHEIManuscript")
source("epiWindSecrets_HEI_4970.R") # contains filepaths and variables that should be kept secret


######## LOADING AND PREPROCESSING RECORDS ########


# loadSingleMatchSet - given match criteria conditions, load the upwind/downwind cohort match records and preprocess data for epidemiological analyses
# INPUTS
#    matchCriteria (int) - distance to road threshold used to derive upwind/downwind matches
#    multipleMatches (boolean) - whether matches are 1:1 (FALSE) or up to 1:4 (TRUE)
#    restricted (boolean) - whether the cohort was restricted to 37 to 42 weeks when deriving matches
# OUTPUTS:
#   R dataframe containing upwind/downwind matches for the specific input match criteria
loadMatchData <- function(matchCriteria, multipleMatches = FALSE,restricted=FALSE) {
  
  
  # filenames contain a 4 if multiple matches were allowed, 1 otherwise
  multipleMatchIndicator <- ifelse(multipleMatches==TRUE,4,1)
  restrictedIndicator <- ifelse(restricted==TRUE,"37to42","all")
      
  # define the full filepath to where matches created under the input constraints are stored and load the data from file
  inFile <- paste(epiFolder,"datasets/wind_epi_", multipleMatchIndicator, "_", matchCriteria,"_",restrictedIndicator,"_","Jan27.csv",sep="")
  inDataset <- read.csv(inFile)
  
  # calculate tertiles of neighborhood income
  inDataset$neigh_inc_tertile <- with(inDataset,ave(median_income_imputeavg5,byear, FUN=function(x).bincode(x,quantile(x,c(0:3/3),na.rm=TRUE), T,T)))
  inDataset[!is.na(inDataset$neigh_inc_tertile),]
  # create indicator variables for low term birth weight, preterm birth, and very perterm birth
  inDataset$ltbw <- ifelse(inDataset$b_wt_cgr < 2500, 1,0)
  inDataset$ptb <- ifelse(inDataset$b_es_ges < 37, 1,0)
  inDataset$vptb <- ifelse(inDataset$b_es_ges < 32, 1,0)
  
  # reset the upwind/downwind classification variable so that the upwind category (value = 1) is the reference category
  inDataset$windCat <- factor(inDataset$windCat)
  inDataset$windCat <- relevel(inDataset$windCat,ref="1")
  return(inDataset)
}


####### RESTRICTING DATASETS #########

# given a set of unique exposures, reduce the number of controls to equal the number of exposed with the same match id
# INPUTS:
#    ctrl (dataframe) - cohort members who are in the control group
#    exposed (dataframe) - cohort members in the exposed group
#    matchId (string) - current match group to get a balanced set of matches for
# OUTPUTS: 
#    a subset of the ctrl group, with n records equaling n exposed in the exposed dataset with the same matchid
getBalancedMatch <- function(ctrl,exposed,matchId) {
  tempExposed <- subset(exposed,matchid ==matchId) # restrict exposed dataset to a specific match id
  tempSubset <- subset(ctrl,matchid == matchId) # restrict control dataset to a specific match id
  # further restrict control dataset to at most n records, where n = number of records in tempSubset
  return(tempSubset[1:length(tempExposed[,1]),])
}

# restrict dataset to matches where both the exposed and control in the match have the same race/ethnicity status
# INPUTS:
#    inDataset (dataframe) - cohort records that contain both control and exposed
#    value (int) - categorical number that indicates the race category to restrict records to 
# OUTPUTS:
#    matchedSet (dataset) - restricted matches for the indicated race/ethnicity strata level 
getRaceSubset <- function(inDataset,value) {
  # restrict records to the race category of interest
  raceSubset <- subset(get(inDataset),b_m_race_eth ==value)
  
  # partition dataset to exposed and control subgroups
  exposed <- subset(raceSubset,windCat==0)
  ctrl <- subset(raceSubset,windCat==1)
  
  # find match ids where both the control and exposed records belong to the indicated race/ethnicity category
  ctrl <- subset(ctrl,matchid %in% unique(exposed$matchid))
  uniqueMatchIds <- unique(ctrl$matchid)
  
  # restrict exposed and control datasets to those with match ids that belong to the indicated race/ethnicity category
  exposed <- subset(exposed,matchid %in% uniqueMatchIds)
  ctrlMatch <- getBalancedMatch(ctrl,exposed,uniqueMatchIds[1])
  for (i in 2:length(uniqueMatchIds)) {
    ctrlMatch <- rbind(ctrlMatch,getBalancedMatch(ctrl,exposed,uniqueMatchIds[i]))
  }
  
  # combine exposed and control into a single dataset
  matchedSet <- rbind(ctrlMatch,exposed)
  return(matchedSet)
}


# restrict dataset to matches where both exposed and control in the match have the same ethnicity
# INPUTS:
#    inDataset (dataframe) - cohort records that contain both control and exposed
#    value (int) - categorical number that indicates the ethnicity category to restrict records to 
# OUTPUTS:
#    matchedSet (dataset) - restricted matches for the indicated ethnicity strata level 
getEthnicitySubset <- function(inDataset,value) {
  # restrict records to the ethnicity category of interest
  ethnicitySubset <- subset(get(inDataset),b_m_hispanic==value)
  
  # partition dataset to exposed and control subgroups
  exposed <- subset(ethnicitySubset,windCat==0)
  ctrl <- subset(ethnicitySubset,windCat==1)
  
  # find match ids where both the control and exposed records belong to the indicated ethnicity category
  ctrl <- subset(ctrl,matchid %in% unique(exposed$matchid))
  uniqueMatchIds <- unique(ctrl$matchid)
  
  # restrict exposed and control datasets to those with match ids that belong to the indicated ethnicity category
  exposed <- subset(exposed,matchid %in% uniqueMatchIds)
  ctrlMatch <- getBalancedMatch(ctrl,exposed,uniqueMatchIds[1])
  for (i in 2:length(uniqueMatchIds)) {
    ctrlMatch <- rbind(ctrlMatch,getBalancedMatch(ctrl,exposed,uniqueMatchIds[i]))
  }
  
  # combine exposed and control into a single dataset
  matchedSet <- rbind(ctrlMatch,exposed)
  return(matchedSet)
}


# restrict dataset to matches where both exposed and control in the match have high school or less education
# INPUTS:
#    inDataset (dataframe) - cohort records that contain both control and exposed
# OUTPUTS:
#    matchedSet (dataset) - restricted matches for the high school or less strata level 
getLowerEducSubset <- function(inDataset) {
  # restrict records to high school or less education
  lowerIncSubset <- subset(get(inDataset),b_m_educ<=2)
  
  # partition dataset to exposed and control subgroups
  exposed <- subset(lowerIncSubset,windCat==0)
  ctrl <- subset(lowerIncSubset,windCat==1)
  
  # find match ids where both the control and exposed records belong to high school or less education
  ctrl <- subset(ctrl,matchid %in% unique(exposed$matchid))
  uniqueMatchIds <- unique(ctrl$matchid)
  
  # restrict exposed and control datasets to those with match ids that belong to high school or less education
  exposed <- subset(exposed,matchid %in% uniqueMatchIds)
  ctrlMatch <- getBalancedMatch(ctrl,exposed,uniqueMatchIds[1])
  for (i in 2:length(uniqueMatchIds)) {
    ctrlMatch <- rbind(ctrlMatch,getBalancedMatch(ctrl,exposed,uniqueMatchIds[i]))
  }
  
  # combine exposed and control into a single dataset
  matchedSet <- rbind(ctrlMatch,exposed)
  return(matchedSet)
}


# restrict dataset to matches where both exposed and control in the match have more than high school or less education
# INPUTS:
#    inDataset (dataframe) - cohort records that contain both control and exposed
# OUTPUTS:
#    matchedSet (dataset) - restricted matches for the more than high school strata level 
getHigherEducSubset <- function(inDataset) {
  # restrict records to more than high school education
  higherIncSubset <- subset(get(inDataset),b_m_educ>2)
  
  # partition datasest to exposed and control subgroups
  exposed <- subset(higherIncSubset,windCat==0)
  ctrl <- subset(higherIncSubset,windCat==1)
  
  # find match ids where both the control and exposed records belong to more than high school education
  ctrl <- subset(ctrl,matchid %in% unique(exposed$matchid))
  uniqueMatchIds <- unique(ctrl$matchid)
  
  # restrict exposed and ccontrol datasets to those with match ids that belong to more than high school education
  exposed <- subset(exposed,matchid %in% uniqueMatchIds)
  ctrlMatch <- getBalancedMatch(ctrl,exposed,uniqueMatchIds[1])
  for (i in 2:length(uniqueMatchIds)) {
    ctrlMatch <- rbind(ctrlMatch,getBalancedMatch(ctrl,exposed,uniqueMatchIds[i]))
  }
  
  # combine exposed and control into a single dataset
  matchedSet <- rbind(ctrlMatch,exposed)
  return(matchedSet)
}


# restrict dataset to matches where both exposed and control in the match belong to the indicated foreign category
# INPUTS:
#    inDataset (dataframe) - cohort records that contain both control and exposed
#    value (int) - categorical number that indicates the category of of interest
# OUTPUTS:
#    matchedSet (dataset) - restricted matches for the indicated foreign category strata level 
getForeignCat <- function(inDataset,value) {
  # restrict records to the indicated foreign category of interest
  foreignSubset <- subset(get(inDataset),b_m_foreign==value)
  
  # partition dataset to exposed and control groups
  exposed <- subset(foreignSubset,windCat==0)
  ctrl <- subset(foreignSubset,windCat==1)
  
  # find match ids where both the control and exposed records belong to the indicated foreign categorical level
  ctrl <- subset(ctrl,matchid %in% unique(exposed$matchid))
  uniqueMatchIds <- unique(ctrl$matchid)
  
  # restrict exposed and control datasets to those with match ids that belong  to the indicated foreign category
  exposed <- subset(exposed,matchid %in% uniqueMatchIds)
  ctrlMatch <- getBalancedMatch(ctrl,exposed,uniqueMatchIds[1])
  for (i in 2:length(uniqueMatchIds)) {
    ctrlMatch <- rbind(ctrlMatch,getBalancedMatch(ctrl,exposed,uniqueMatchIds[i]))
  }
  
  # combine exposed and control into a single dataset
  matchedSet <- rbind(ctrlMatch,exposed)
  return(matchedSet)
}


# restrict dataset to matches where both exposed and control in the match belong to the neighborhood income category
# INPUTS:
#    inDataset (dataframe) - cohort records that contain both control and exposed
#    value (int) - categorical number that indicates the category of of interest
# OUTPUTS:
#    matchedSet (dataset) - restricted matches for the indicated neighborhood income strata level 
getIncomeCat <- function(inDataset,value) {
  # restrict records to the indicated neighborhood income level
  incomeSubset <- subset(get(inDataset),neigh_inc_tertile==value)
  
  # partition dataset to control and exposed groups
  exposed <- subset(incomeSubset,windCat==0)
  ctrl <- subset(incomeSubset,windCat==1)
  
  # find match ids where both the control and exposed records belong to the indicated neighborhood income level
  ctrl <- subset(ctrl,matchid %in% unique(exposed$matchid))
  uniqueMatchIds <- unique(ctrl$matchid)
  
  # restrict exposed and control datasets to those with match ids that belong to the indicated neighborhood income level
  exposed <- subset(exposed,matchid %in% uniqueMatchIds)
  ctrlMatch <- getBalancedMatch(ctrl,exposed,uniqueMatchIds[1])
  for (i in 2:length(uniqueMatchIds)) {
    ctrlMatch <- rbind(ctrlMatch,getBalancedMatch(ctrl,exposed,uniqueMatchIds[i]))
  }
  
  # combine exposed and control into a single dataset
  matchedSet <- rbind(ctrlMatch,exposed)
  return(matchedSet)
}


# get the number of cohort pairs within a [lowerDist, upperDist] distance range 
# INPUTS:
#    inDataset (dataframe) - cohort records from all distances
#    lowerDist (float) - lower threshold to restrict records to (inclusive)
#    upperDist (float) - upper threshold to restrict records to (exclusive)
# OUTPUTS:
#    number of matched pairs within the restricted distance interval
getSampleSizeRestrictedDistance <- function(inDataset,lowerDist,upperDist) {
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDist) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDist) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  return(length(restrictedData[,1])/2)
}


# get the number of events within a [lowerDist, upperDist] distance interval
# INPUTS:
#    inDataset (dataframe) - cohort records from all distances
#    lowerDist (float) - lower threshold to restrict records to (inclusive)
#    upperDist (float) - upper threshold to restrict records to (exclusive)
# OUTPUTS:
#    number of events within the restricted distance interval
getEventsRestrictedDistance <- function(inDataset,lowerDist,upperDist,outcome) {
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDist) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDist) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  return(sum(restrictedData[,outcome]))
}


# get the number of cohort pairs to births within a [firstYear, lastYear] time interval
# INPUTS:
#    inDataset (dataframe) - cohort records from all years
#    firstYear (int) - first birth year to restrict records to (inclusive)
#    lastYear (int) - last birth year to restrict records to (inclusive)
# OUTPUTS:
#    number of matched pairs within the restricted distance interval
getSampleSizeRestrictedTime <- function(inDataset,firstYear,lastYear) {
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,byear>=firstYear) # restrict records to those with birth years at least greater than or equal to firstYear
  tempData <- subset(tempData,byear<=lastYear) # further restrict records to those with birth years at least less than or equal to lastYear
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  return(length(restrictedData[,1])/2)
}


# get the number of events within a [firstYear, lastYear] time interval
# INPUTS:
#    inDataset (dataframe) - cohort records from all year
#    firstYear (int) - first birth year to restrict records to (inclusive)
#    lastYear (int) - last birth year to restrict records to (inclusive)
# OUTPUTS:
#    number of events within the restricted distance interval
getEventsestrictedTime <- function(inDataset,firstYear,lastYear,outcome) {
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,byear>=firstYear) # restrict records to those with birth years at least greater than or equal to firstYear
  tempData <- subset(tempData,byear<=lastYear) # further restrict records to those with birth years at least less than or equal to lastYear
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  return(sum(restrictedData[,outcome]))
}

########### FUNCTIONS FOR LINEAR REGRESSION MODELS ############

# get the coefficient and 95% CI for the wind treatment group variable in the linear regression model 
# INPUTS:
#    matchCriteria (int) - criteria used to match upwind/downwind participants. Can be 15,25,50, or 100
#    outcome (string) - variable name of the outcome of interest
#    resample (boolean) - whether controls were matched to multiple exposed
#    fullModel (boolean) - whether to run the base model with minimal covariates or fully adjusted model
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - indicator variable of whether the base (0) or fully adjusted model (1) was run
#       - indicator variable of whether controls were matched to multiple exposed (4) or not (0)
#       - wind coefficient and 95% CI
getWindCoefficientsLinear <- function(matchCriteria,outcome,resample=FALSE,fullModel=FALSE) {
  sampleCat <- ifelse(resample==TRUE,4,1) # 4 if controls were matched to multiple exposed, 0 ow
  fullCat <- ifelse(fullModel==TRUE,1,0) # 1 if running the fully adjusted model, 0 ow
  datasetName <- paste("epi_37_42_",sampleCat,"_",matchCriteria,sep="") # which dataset to pull data from.  Each combination 
                                                                    # of match criteria and resampling has a unique dataset name
  # run the fully adjusted model
  if(fullModel)
    { 
      model <- runFullLinearModel(datasetName,outcome)
  }
  
  
  # run the base model
  else {
    model <- runBaseLinearModel(datasetName,outcome)
  }
    return(c(
      length(get(datasetName)[,1])/2,
      fullCat,
      sampleCat,
      matchCriteria,
      summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
      ))
}

# get wind coefficient and 95% CI for the wind treatment group variable in the linear regression model, for a subset restricted to a distance from nearest road 
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name containing the birth outcome of interest
#    lowerDist (int) - lower distance threshold to restrict cohort records by
#    upperDist (int) - upper distance threshold to restrict cohort records by
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - lowerDist
#       - upperDist
#       - wind coefficient and 95% CI
getWindCoefficientsDistanceLinear <- function(inDataset,outcome,lowerDist,upperDist) {
  model <- runDistanceModelLinear(inDataset,outcome,lowerDist,upperDist)
  sampleSize <- getSampleSizeRestrictedDistance(inDataset,lowerDist,upperDist)
  return(c(
    sampleSize,
    lowerDist,
    upperDist,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the wind treatment group variable in the linear regression model with shielding covariates,
# for a subset restricted to a distance from nearest road
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name containing the birth outcome of interest
#    lowerDist (int) - lower distance threshold to restrict cohort records by
#    upperDist (int) - upper distance threshold to restrict cohort records by
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - lowerDist
#       - upperDist
#       - wind coefficient and 95% CI
getWindCoefficientsShieldingDistanceLinear <- function(inDataset,outcome,lowerDist,upperDist) {
  model <- runDistanceShieldingModelLinear(inDataset,outcome,lowerDist,upperDist)
  sampleSize <- getSampleSizeRestrictedDistance(inDataset,lowerDist,upperDist)
  return(c(
    sampleSize,
    lowerDist,
    upperDist,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the continuous % hours downwind variable in the linear regression model
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name containing the birth outcome of interest
#    lowerDist (int) - lower distance threshold to restrict cohort records by
#    upperDist (int) - upper distance threshold to restrict cohort records by
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - lowerDist
#       - upperDist
#       - wind coefficient and 95% CI
getPercentWindCoefficientsDistanceLinear <- function(inDataset,outcome,lowerDist,upperDist) {
  model <- runPercWindDistanceModelLinear(inDataset,outcome,lowerDist,upperDist)
  sampleSize <- getSampleSizeRestrictedDistance(inDataset,lowerDist,upperDist)
  return(c(
    sampleSize,
    lowerDist,
    upperDist,
    summary(model)$coefficients['allmwn',1]*10 + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['allmwn',2]*10
    ))
}

# get wind coefficient and 95% CI for the linear regression models for a subset of records restricted by time
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name containing the birth outcome of interest
#    firstYear (int) - first year (inclusive) of the restricted subset
#    lastYear (int) - last year (inclusive) of the restricted subset
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - firstYear
#       - lastYear
#       - wind coefficient and 95% CI
getWindCoefficientsTimeLinear <- function(inDataset,outcome,firstYear,lastYear) {
  model <- runTimeSensitivityModelLinear(inDataset,outcome,firstYear,lastYear)
  sampleSize <- getSampleSizeRestrictedTime(inDataset,firstYear,lastYear)
  return(c(
    sampleSize,
    firstYear,
    lastYear,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the linear regression model without the foreign born covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing restricted cohort records
#    outcome (string) - variable name containing the birth outcome of interest
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - wind coefficient and 95% CI
runForeignSensitivityModelLinear <- function(inDataset,outcome) {
  model <- (lm(get(outcome) ~ 
                 windCat + 
                 factor(b_es_ges) + 
                 b_m_age + 
                 factor(b_m_race_eth) + 
                 factor(b_m_hispanic) + 
                 factor(b_m_educ2) +
                 #factor(b_m_foreign) + 
                 factor(pay) + 
                 factor(b_m_wic) + 
                 factor(b_m_cig) + 
                 factor(b_mo_pcb) + 
                 b_m_wtg + 
                 factor(neigh_inc_tertile) +
                 factor(bmonth) + 
                 factor(byear) + 
                 factor(county_code) + 
                 near_dist,
            data=get(inDataset)))
  sampleSize <- length(get(inDataset)[,1])/2
  return(c(
    sampleSize,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the linear regression model without the neighborhood income covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing restricted cohort records
#    outcome (string) - variable name containing the birth outcome of interest
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - wind coefficient and 95% CI
runIncSensitivityModelLinear <- function(inDataset,outcome) {
  model <- (lm(get(outcome) ~ 
                 windCat + 
                 factor(b_es_ges) + 
                 b_m_age + 
                 factor(b_m_race_eth) + 
                 factor(b_m_hispanic) + 
                 factor(b_m_educ2) +
                 factor(b_m_foreign) + 
                 factor(pay) + 
                 factor(b_m_wic) + 
                 factor(b_m_cig) + 
                 factor(b_mo_pcb) + 
                 b_m_wtg + 
                 #factor(neigh_inc_tertile) +
                 factor(bmonth) + 
                 factor(byear) + 
                 factor(county_code) + 
                 near_dist,
               data=get(inDataset)))
  sampleSize <- length(get(inDataset)[,1])/2
  return(c(
    sampleSize,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the linear regression model without the education covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing restricted cohort records
#    outcome (string) - variable name containing the birth outcome of interest
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - wind coefficient and 95% CI
runEducSensitivityModelLinear <- function(inDataset,outcome) {
  model <- (lm(get(outcome) ~ 
                 windCat + 
                 factor(b_es_ges) + 
                 b_m_age + 
                 factor(b_m_race_eth) + 
                 factor(b_m_hispanic) + 
                 #factor(b_m_educ2) +
                 factor(b_m_foreign) + 
                 factor(pay) + 
                 factor(b_m_wic) + 
                 factor(b_m_cig) + 
                 factor(b_mo_pcb) + 
                 b_m_wtg + 
                 factor(neigh_inc_tertile) +
                 factor(bmonth) + 
                 factor(byear) + 
                 factor(county_code) + 
                 near_dist,
               data=get(inDataset)))
  sampleSize <- length(get(inDataset)[,1])/2
  return(c(
    sampleSize,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the linear regression model without the race or ethnicity covariates
# INPUTS:
#    inDataset (string) - name of the dataset containing restricted cohort records
#    outcome (string) - variable name containing the birth outcome of interest
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - wind coefficient and 95% CI
runRaceSensitivityModelLinear <- function(inDataset,outcome) {
  model <- (lm(get(outcome) ~ 
                 windCat + 
                 factor(b_es_ges) + 
                 b_m_age + 
                 #factor(b_m_race_eth) + 
                 #factor(b_m_hispanic) + 
                 factor(b_m_educ2) +
                 factor(b_m_foreign) + 
                 factor(pay) + 
                 factor(b_m_wic) + 
                 factor(b_m_cig) + 
                 factor(b_mo_pcb) + 
                 b_m_wtg + 
                 factor(neigh_inc_tertile) +
                 factor(bmonth) + 
                 factor(byear) + 
                 factor(county_code) + 
                 near_dist,
               data=get(inDataset)))
  sampleSize <- length(get(inDataset)[,1])/2
  return(c(
    sampleSize,
    summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# get wind coefficient and 95% CI for the linear regression model without the ethnicity covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing restricted cohort records
#    outcome (string) - variable name containing the birth outcome of interest
# OUTPUTS:
#    vector with the following information:
#       - number of pairs/matches
#       - wind coefficient and 95% CI
runEthnicitySensitivityModelLinear <- function(inDataset,outcome) {
  model <- (lm(get(outcome) ~ 
                 windCat + 
                 factor(b_es_ges) + 
                 b_m_age + 
                 factor(b_m_race_eth) + 
                 #factor(b_m_hispanic) + 
                 factor(b_m_educ2) +
                 factor(b_m_foreign) + 
                 factor(pay) + 
                 factor(b_m_wic) + 
                 factor(b_m_cig) + 
                 factor(b_mo_pcb) + 
                 b_m_wtg + 
                 factor(neigh_inc_tertile) +
                 factor(bmonth) + 
                 factor(byear) + 
                 factor(county_code) + 
                 near_dist,
               data=get(inDataset)))
  sampleSize <- length(get(inDataset)[,1])/2
  return(c(
    sampleSize,summary(model)$coefficients['windCat0',1] + qnorm(c(0.025,0.5,0.975))*summary(model)$coefficients['windCat0',2]
    ))
}

# run the linear regression model which includes the wind shielding covariates
# INPUTS:
#    inDataset (string) - name of the dataset containing restricted cohort records
#    outcome (string) - variable name containing the birth outcome of interest
# OUTPUTS:
#    lm object
runShieldingModelLinear <- function(inDataset,outcome) {
  return(lm(get(outcome) ~ 
              windCat + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) +
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) +
              factor(bmonth) + 
              factor(byear) + 
              factor(county_code) + 
              near_dist + 
              allcutTr +
              allcutBd,
            data=get(inDataset)))
}


# run a linear regression model with the "base" set of covariates
# INPUTS:
#    inDataset (string) - name of the dataset to use for the linear regression
#    outcome (string) - variable name for the outcome of interest
# OUTPUTS:
#    lm object
runBaseLinearModel <- function(inDataset,outcome) {
  return (lm(get(outcome) ~ 
               windCat + 
               factor(bs_sex) + 
               b_m_age + 
               factor(bmonth) +
               factor(byear) + 
               factor(county_code) + 
               factor(b_es_ges) +
               near_dist, 
             data=get(inDataset)))
}


# run a linear regression model with the "full" set of covariates
# INPUTS:
#    inDataset (string) - name of the dataset to use for the linear regression
#    outcome (string) - variable name for the outcome of interest
# OUTPUTS:
#    lm object
runFullLinearModel <- function(inDataset,outcome) {
  return(lm(get(outcome) ~ 
              windCat + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) +
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) +
              factor(bmonth) + 
              factor(byear) + 
              factor(county_code) + 
              near_dist,
            data=get(inDataset)))
} 


# model for generating linear model predictions. The covariates are the same as the 'full', with the exception of the variable 'county code'
# which had to be dropped to create predictions
# INPUTS:
#    inDataset (string) - name of the dataset to use for the linear regression
#    outcome (string) - variable name for the outcome of interest
# OUTPUTS:
#    lm object
runPredictLinearModel <- function(inDataset,outcome) {
  return(lm(get(outcome) ~ 
              windCat + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) +
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) +
              factor(bmonth) + 
              factor(byear) + 
              near_dist,
            data=get(inDataset)))
} 


# linear regression model with the full set of covariates, but restricted to a distance to the nearest road interval
# INPUTS:
#    inDataset (string) - name of the dataframe to run the linear regression model on
#    outcome (string) - name of the outcome variable
#    lowerDistance (int) - lower bound threshold for the interval of interest
#    upperDistance (int) - upper bound threshold for the interval of interest
# OUTPUTS:
#    lm model
runDistanceModelLinear <- function(inDataset,outcome,lowerDistance,upperDistance) {
  
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDistance) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDistance) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  

  # run linear regression model on restricted data
  return(lm(get(outcome) ~ 
              windCat + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) +
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) +
              factor(bmonth) + 
              factor(byear) + 
              factor(county_code) + 
              near_dist,
            data=restrictedData))
}

# linear regression model with the full set of covariates, but restricted to a distance to the nearest road interval
# INPUTS:
#    inDataset (string) - name of the dataframe to run the linear regression model on
#    outcome (string) - name of the outcome variable
#    lowerDistance (int) - lower bound threshold for the interval of interest
#    upperDistance (int) - upper bound threshold for the interval of interest
# OUTPUTS:
#    lm model
runPercWindDistanceModelLinear <- function(inDataset,outcome,lowerDistance,upperDistance) {
  
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDistance) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDistance) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  
  
  # run linear regression model on restricted data
  return(lm(get(outcome) ~ 
              allmwn + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) +
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) +
              factor(bmonth) + 
              factor(byear) + 
              factor(county_code) + 
              near_dist,
            data=restrictedData))
}



# linear regression model with the full set of covariates, but restricted to a distance to the nearest road interval
# INPUTS:
#    inDataset (string) - name of the dataframe to run the linear regression model on
#    outcome (string) - name of the outcome variable
#    lowerDistance (int) - lower bound threshold for the interval of interest
#    upperDistance (int) - upper bound threshold for the interval of interest
# OUTPUTS:
#    lm model
runDistanceShieldingModelLinear <- function(inDataset,outcome,lowerDistance,upperDistance) {
  
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDistance) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDistance) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  

  # run linear regression model on restricted data
  return(lm(get(outcome) ~ 
              windCat  + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) +
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) +
              factor(bmonth) + 
              factor(byear) + 
              factor(county_code) + 
              allcutTr + allcutBd +
              near_dist,
          data=restrictedData))
  }


# restrict records to a time subset of interest and run a linear regression model on the subset
# INPUTS:
#    inDataset (string) - name of the dataset containing the cohort records
#    outcome (string) - variable name for the outcome of interest
#    firstYear (int) - first year (inclusive) for the restricted time range of interest
#    lastYear (int) - last year (inclusive) for the restricted time range of interest
# OUTPUTS:
#    lm model
runTimeSensitivityModelLinear <- function(inDataset,outcome,firstYear,lastYear) {
  inData <- get(inDataset)
  tempData <- subset(inData,byear>=firstYear)
  tempData <- subset(tempData,byear<=lastYear)
  matchesInYear <- unique(tempData$matchid)
  restrictedData <- subset(inData, matchid %in% matchesInYear)
  cat(length(restrictedData[,1]))
  return(lm(get(outcome) ~ 
              windCat + 
              factor(b_es_ges) + 
              b_m_age + 
              factor(b_m_race_eth) + 
              factor(b_m_hispanic) + 
              factor(b_m_educ2) + 
              factor(b_m_foreign) + 
              factor(pay) + 
              factor(b_m_wic) + 
              factor(b_m_cig) + 
              factor(b_mo_pcb) + 
              b_m_wtg + 
              factor(neigh_inc_tertile) + 
              factor(bmonth) + 
              factor(byear) + 
              factor(county_code),
            data=restrictedData))
}


############# LOGISTIC REGRESSION MODELS #############


# get the coefficient and 95% CI for the wind treatment variable in a logistic regression models, restricted to a 
# subset of records within a distance interval from nearest road
# INPUTS:
#    inDataset (string) - name of the dataset containing the cohort records
#    outcome (string) - variable name of the outcome of interest
#    lowerDist (int) - lower distance (inclusive) threshold for the subset group
#    upperDist (int) - upper distance (exclusive) threshold for the subset group
# OUTPUTS:
#    a vector containing the following information:
#     - number of matches in the restricted group
#     - number of outcome events in the restricted group
#     - lowerDist
#     - upperDist
#     - OR and 95CI for the wind treatment variable
getWindCoefficientsDistanceLogistic <- function(inDataset,outcome,lowerDist,upperDist) {
  model <- runDistanceModelLogistic(inDataset,outcome,lowerDist,upperDist)
  sampleSize <- getSampleSizeRestrictedDistance(inDataset,lowerDist,upperDist)
  n_events <- getEventsRestrictedDistance(inDataset,lowerDist,upperDist,outcome)
  return(c(sampleSize,n_events,lowerDist,upperDist,summary(model)$conf.int['windCat0',]))
}


# get the OR and 95% CI for the wind treatment variable in a logistic regression model
# INPUTS:
#    matchCriteria (int) - criteria used to create upwind/downwind matches.  Can be 15,25,50, or 100
#    outcome (string) - variable name of the outcome of interest
#    resample (boolean) - whether to match multiple exposed to a single control
#    fullModel (boolean) - whether to run a model with minimal or a full set of covariates
#    restricted (boolean) - whether cohort records should be restricted to 37 to 42 weeks of gestational age
# OUTPUTS:
#    vector containing the following information:
#     - number of matches
#     - number of events 
#     - whether model was the basic set of covariates or the full set
#     - whether controls were matched to multiple exposed 
#     - matchCriteria
#     - estimated OR with 95% CI
getWindCoefficientsLogistic <- function(matchCriteria,outcome,resample=FALSE,fullModel=FALSE,restricted=FALSE) {
  sampleCat <- ifelse(resample==TRUE,4,1)
  fullCat <- ifelse(fullModel==TRUE,1,0)
  restricted <- ifelse(restricted==TRUE,"epi_37_42_","epi_all_")
  datasetName <- paste(restricted,sampleCat,"_",matchCriteria,sep="")
  n_events <- sum(get(datasetName)[,outcome])
  if(fullModel)
  { 
    model <- runFullLogisticModel(datasetName,outcome)
  }
  else {
    model <- runBaseLogisticModel(datasetName,outcome)
  }
  return(c(length(get(datasetName)[,1])/2,n_events,fullCat,sampleCat,matchCriteria,summary(model)$conf.int['windCat0',]))
}

# restrict a cohort to a time subset and get the OR and 95% CI for the wind treatment variable in the restricted subset
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
#    firstYear (int) - first year (inclusive) to restrict cohort records to
#    lastYear (int) - last year (inclusive) to restrict cohort records to
# OUTPUTS:
#    vector containing the following information:
#     - number of matches in the restricted subset
#     - number of events in the restricted subset
#     - firstYear
#     - lastYear
#     - wind treatment OR and 95% CI
getWindCoefficientsTimeLogistic <- function(inDataset,outcome,firstYear,lastYear) {
  model <- runTimeSensitivityModelLogistic(inDataset,outcome,firstYear,lastYear)
  sampleSize <- getSampleSizeRestrictedTime(inDataset,firstYear,lastYear)
  nEvents <- getEventsestrictedTime(inDataset,firstYear,lastYear,outcome)
  return(c(
    sampleSize,
    nEvents,
    firstYear,
    lastYear,
    summary(model)$conf.int['windCat0',]
    ))
}


# get wind treatment group OR and 95% CI for the logistic regression model that includes shielding covariates
# INPUTS:
#    inDataset (string) - name of the dataset containing the cohort records
#    outcome (string) - variable name of the outcome of interest
#    lowerDist (int) - lower distance (inclusive) threshold for the subset group
#    upperDist (int) - upper distance (exclusive) threshold for the subset group
# OUTPUTS:
#    a vector containing the following information:
#     - number of matches in the restricted group
#     - number of outcome events in the restricted group
#     - lowerDist
#     - upperDist
#     - OR and 95CI for the wind treatment variable
getWindCoefficientsShieldingDistanceLogistic <- function(inDataset,outcome,lowerDist,upperDist) {
  model <- runDistanceShieldingModelLogistic(inDataset,outcome,lowerDist,upperDist)
  sampleSize <- getSampleSizeRestrictedDistance(inDataset,lowerDist,upperDist)
  n_events <- getEventsRestrictedDistance(inDataset,lowerDist,upperDist,outcome)
  return(c(
    sampleSize,
    n_events,
    lowerDist,
    upperDist,
    summary(model)$conf.int['windCat0',]
    ))
}

# get the OR and 95% CI for a logistic regression models that replaces the wind treatment group with the % of pregnancy downwind variable,
# restricted to a distance interval of interest
# INPUTS:
#    inDataset (string) - name of the dataset containing the cohort records
#    outcome (string) - variable name of the outcome of interest
#    lowerDist (int) - lower distance (inclusive) threshold for the subset group
#    upperDist (int) - upper distance (exclusive) threshold for the subset group
# OUTPUTS:
#    a vector containing the following information:
#     - number of matches in the restricted group
#     - number of outcome events in the restricted group
#     - lowerDist
#     - upperDist
#     - OR and 95CI for the % hours downwind variable
getPercentWindCoefficientsDistanceLogistic <- function(inDataset,outcome,lowerDist,upperDist) {
  model <- runPercWindDistanceModelogistic(inDataset,outcome,lowerDist,upperDist)
  sampleSize <- getSampleSizeRestrictedDistance(inDataset,lowerDist,upperDist)
  n_events <- getEventsRestrictedDistance(inDataset,lowerDist,upperDist,outcome)
  return(c(
    sampleSize,
    n_events,
    lowerDist,
    upperDist,
    summary(model)$conf.int['allmwn',]
    ))
}

# logistic regression model with the full set of covariates, but restricted to a distance to the nearest road interval
# INPUTS:
#    inDataset (string) - name of the dataframe to run the linear regression model on
#    outcome (string) - name of the outcome variable
#    lowerDistance (int) - lower bound threshold for the interval of interest
#    upperDistance (int) - upper bound threshold for the interval of interest
# OUTPUTS:
#    lm model
runDistanceShieldingModelLogistic <- function(inDataset,outcome,lowerDistance,upperDistance) {
  
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDistance) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDistance) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  
  
  # run linear regression model on restricted data
  return(
    clogit(get(outcome) ~ 
             windCat + 
             b_m_age + 
             factor(b_m_race_eth) + 
             factor(b_m_hispanic) + 
             factor(b_m_educ2) + 
             factor(b_m_foreign) + 
             factor(pay) + 
             factor(b_m_wic) + 
             factor(b_m_cig) + 
             factor(b_mo_pcb) + 
             b_m_wtg + 
             factor(neigh_inc_tertile) + 
             factor(bmonth) + 
             factor(byear) + 
             near_dist + 
             strata(matchid) + 
             factor(county_code) + 
             allcutTr + 
             allcutBd,
           data=restrictedData,method="approximate"))
}


# logistic regression model with the full set of covariates, but restricted to a distance to the nearest road interval
# INPUTS:
#    inDataset (string) - name of the dataframe to run the linear regression model on
#    outcome (string) - name of the outcome variable
#    lowerDistance (int) - lower bound threshold for the interval of interest
#    upperDistance (int) - upper bound threshold for the interval of interest
# OUTPUTS:
#    lm model
runPercWindDistanceModelogistic <- function(inDataset,outcome,lowerDistance,upperDistance) {
  
  # restrict the dataset to records within the interval of interest
  inData <- get(inDataset)
  tempData <- subset(inData,near_dist>=lowerDistance) # get records with distances above the lower threshold
  tempData <- subset(tempData,near_dist<upperDistance) # further restrict records to those also below the upper threshold
  matchesInDistance <- unique(tempData$matchid) # get match ids that meet restriction criteria
  restrictedData <- subset(inData, matchid %in% matchesInDistance) # get records with restricted match ids
  
  
  # run linear regression model on restricted data
  return(
    clogit(get(outcome) ~ 
             allmwn + 
             b_m_age + 
             factor(b_m_race_eth) + 
             factor(b_m_hispanic) + 
             factor(b_m_educ2) + 
             factor(b_m_foreign) + 
             factor(pay) + 
             factor(b_m_wic) + 
             factor(b_m_cig) + 
             factor(b_mo_pcb) + 
             b_m_wtg + 
             factor(neigh_inc_tertile) + 
             factor(bmonth) + 
             factor(byear) + 
             near_dist + 
             strata(matchid) + 
             factor(county_code), 
           data=restrictedData,method="approximate"))
}

# run the logistic regression model with a minimal or "base" set of covariates
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name for the outcome of interest
# OUTPUTS:
#    conditional logistic regression object
runBaseLogisticModel <- function(inDataset,outcome) {
  return(clogit(get(outcome) ~ 
                  windCat + 
                  factor(bs_sex) + 
                  b_m_age + 
                  factor(bmonth) + 
                  factor(byear) + 
                  near_dist + 
                  strata(matchid) + 
                  factor(county_code), 
                data = get(inDataset),method="approximate"))
} 

# run the logistic regression model with a full set of covariates
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name for the outcome of interest
# OUTPUTS:
#    conditional logistic regression object
runFullLogisticModel <- function(inDataset,outcome) {
  return(clogit(get(outcome) ~ 
                  windCat + 
                  b_m_age + 
                  factor(b_m_race_eth) + 
                  factor(b_m_hispanic) + 
                  factor(b_m_educ2) + 
                  factor(b_m_foreign) + 
                  factor(pay) + 
                  factor(b_m_wic) + 
                  factor(b_m_cig) + 
                  factor(b_mo_pcb) + 
                  b_m_wtg + 
                  factor(neigh_inc_tertile) + 
                  factor(bmonth) + 
                  factor(byear) + 
                  near_dist + 
                  strata(matchid) + 
                  factor(county_code),
                data=get(inDataset),method="approximate"))
} 


# restrict cohort records to a distance to nearest road range and run a logistic regression model on the restricted dataset
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name for the outcome of interest
#    lowerDistance (int) - lower distance (inclusive) for the restricted dataset
#    upperDistance (int) - upper distance (inclusive) for the restricted dataset
# OUTPUTS:
#    conditional logistic regression model object
runDistanceModelLogistic <- function(inDataset,outcome,lowerDistance,upperDistance) {
  inData <- get(inDataset)
  
  tempData <- subset(inData,near_dist>=lowerDistance) # restrict to records with distance to nearest road above or equal to lower threshold
  tempData <- subset(tempData,near_dist<upperDistance) # further restrict records to distance to nearest road below upper threshold
  
  matchesInDistance <- unique(tempData$matchid) # get match ids where at least one record in the match is within the distance interval
  restrictedData <- subset(inData, matchid %in% matchesInDistance)
  
  return(clogit(get(outcome) ~ 
                  windCat + 
                  b_m_age + 
                  factor(b_m_race_eth) + 
                  factor(b_m_hispanic) + 
                  factor(b_m_educ2) + 
                  factor(b_m_foreign) + 
                  factor(pay) + 
                  factor(b_m_wic) + 
                  factor(b_m_cig) + 
                  factor(b_mo_pcb) + 
                  b_m_wtg + 
                  factor(neigh_inc_tertile) + 
                  factor(bmonth) + 
                  factor(byear) + 
                  near_dist + 
                  strata(matchid) + 
                  factor(county_code),
                data=restrictedData,method="approximate"))
}


# run a logistic regression model without the race or ethnicity covariates
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
# OUTPUTS:
#    a vector containing the following information:
#    - number of pairs 
#    - number of events
#    - OR and 95% CI for the wind treatment group variable
runRaceSensitivityModelLogistic <- function(inDataset,outcome) {
  model <- clogit(get(outcome) ~ 
                    windCat + 
                    b_m_age + 
                    factor(b_m_educ2) + 
                    factor(b_m_foreign) + 
                    factor(pay) + 
                    factor(b_m_wic) + 
                    factor(b_m_cig) + 
                    factor(b_mo_pcb) + 
                    b_m_wtg + 
                    factor(neigh_inc_tertile) + 
                    factor(bmonth) + 
                    factor(byear) + 
                    near_dist +
                    strata(matchid) + 
                    factor(county_code),
                  data=get(inDataset),method="approximate")
  
  sampleSize <- length(get(inDataset)[,1])/2
  n_events <- getEventsRestrictedDistance(inDataset,0,600,outcome)
  return(c(
    sampleSize,
    n_events,
    summary(model)$conf.int['windCat0',]
    ))
}

# run a logistic regression model without the ethnicity covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
# OUTPUTS:
#    a vector containing the following information:
#    - number of pairs 
#    - number of events
#    - OR and 95% CI for the wind treatment group variable
runEthnicitySensitivityModelLogistic <- function(inDataset,outcome) {
  model <- clogit(get(outcome) ~ 
                    windCat + 
                    b_m_age + 
                    factor(b_m_race_eth) + 
                    factor(b_m_educ2) + 
                    factor(b_m_foreign) + 
                    factor(pay) + 
                    factor(b_m_wic) + 
                    factor(b_m_cig) + 
                    factor(b_mo_pcb) + 
                    b_m_wtg + 
                    factor(neigh_inc_tertile) + 
                    factor(bmonth) + 
                    factor(byear) + 
                    near_dist +
                    strata(matchid) + 
                    factor(county_code),
                  data=get(inDataset),method="approximate")
  
  sampleSize <- length(get(inDataset)[,1])/2
  n_events <- getEventsRestrictedDistance(inDataset,0,600,outcome)
  return(c(
    sampleSize,
    n_events,
    summary(model)$conf.int['windCat0',]
    ))
}


# run a logistic regression model without the education covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
# OUTPUTS:
#    a vector containing the following information:
#    - number of pairs 
#    - number of events
#    - OR and 95% CI for the wind treatment group variable
runEducSensitivityModelLogistic <- function(inDataset,outcome) {
  model <- clogit(get(outcome) ~ 
                    windCat + 
                    b_m_age + 
                    factor(neigh_inc_tertile) + 
                    factor(b_m_race_eth) + 
                    factor(b_m_hispanic) + 
                    factor(b_m_foreign) + 
                    factor(pay) + 
                    factor(b_m_wic) + 
                    factor(b_m_cig) + 
                    factor(b_mo_pcb) + 
                    b_m_wtg + 
                    factor(bmonth) + 
                    factor(byear) + 
                    near_dist +
                    strata(matchid) + 
                    factor(county_code),
                  data=get(inDataset),method="approximate")
  
  sampleSize <- length(get(inDataset)[,1])/2
  n_events <- getEventsRestrictedDistance(inDataset,0,600,outcome)
  return(c(
    sampleSize,
    n_events,
    summary(model)$conf.int['windCat0',]
    ))
}


# run a logistic regression model without the neighborhood income covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
# OUTPUTS:
#    a vector containing the following information:
#    - number of pairs 
#    - number of events
#    - OR and 95% CI for the wind treatment group variable
runIncomeSensitivityModelLogistic <- function(inDataset,outcome) {
  model <- clogit(get(outcome) ~ 
                    windCat + 
                    b_m_age + 
                    factor(b_m_race_eth) + 
                    factor(b_m_hispanic) + 
                    factor(b_m_educ2) + 
                    factor(b_m_foreign) + 
                    factor(pay) + 
                    factor(b_m_wic) + 
                    factor(b_m_cig) + 
                    factor(b_mo_pcb) + 
                    b_m_wtg + 
                    factor(bmonth) + 
                    factor(byear) + 
                    near_dist +
                    strata(matchid) + 
                    factor(county_code),
                  data=get(inDataset),method="approximate")
  
  sampleSize <- length(get(inDataset)[,1])/2
  n_events <- getEventsRestrictedDistance(inDataset,0,600,outcome)
  return(c(
    sampleSize,
    n_events,
    summary(model)$conf.int['windCat0',]
    ))
}


# run a logistic regression model without the mother born in the US covariate
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
# OUTPUTS:
#    a vector containing the following information:
#    - number of pairs 
#    - number of events
#    - OR and 95% CI for the wind treatment group variable
runForeignSensitivityModelLogistic <- function(inDataset,outcome) {
  model <- clogit(get(outcome) ~ 
                    windCat + 
                    b_m_age + 
                    factor(neigh_inc_tertile) + 
                    factor(b_m_race_eth) + 
                    factor(b_m_hispanic) + 
                    factor(b_m_educ2) + 
                    factor(pay) + 
                    factor(b_m_wic) + 
                    factor(b_m_cig) + 
                    factor(b_mo_pcb) + 
                    b_m_wtg + 
                    factor(bmonth) + 
                    factor(byear) + 
                    near_dist +
                    strata(matchid) + 
                    factor(county_code),
                  data=get(inDataset),method="approximate")
  sampleSize <- length(get(inDataset)[,1])/2
  
  n_events <- getEventsRestrictedDistance(inDataset,0,600,outcome)
  return(c(
    sampleSize,
    n_events,
    summary(model)$conf.int['windCat0',]
    ))
}

# restrict cohort to a time interval of interest and run a logistic regression model on the restricted subset
# INPUTS:
#    inDataset (string) - name of the dataset containing cohort records
#    outcome (string) - variable name of the outcome of interest
#    firstYear (int) - first year (inclusive) to restrict records to
#    lastYear (int) - last year (inclusive) to restrict records to
# OUTPUTS:
#    considitional logistic regression object
runTimeSensitivityModelLogistic <- function(inDataset,outcome,firstYear,lastYear) {
  inData <- get(inDataset)
  
  tempData <- subset(inData,byear>=firstYear) # restrict cohort records to those where the birth year is >= firstYear
  tempData <- subset(tempData,byear<=lastYear) # further restrict the cohort to births <= lastYear
  
  
  matchesInYear <- unique(tempData$matchid) # get match ids where at least one of the records in the match is within the time interval
  restrictedData <- subset(inData, matchid %in% matchesInYear)
  
  return(clogit(get(outcome) ~ 
                  windCat + 
                  b_m_age + 
                  factor(b_m_race_eth) + 
                  factor(b_m_hispanic) + 
                  factor(b_m_educ2) +
                  factor(b_m_foreign) + 
                  factor(pay) + 
                  factor(b_m_wic) + 
                  factor(b_m_cig) + 
                  factor(b_mo_pcb) + 
                  b_m_wtg + 
                  factor(neigh_inc_tertile) +
                  factor(bmonth) + 
                  factor(byear) + 
                  factor(county_code) + 
                  strata(matchid),
                data=restrictedData,method="approximate"))
}



######## END OF epiFunctionsForWindStudy_HEI_4970.R #########