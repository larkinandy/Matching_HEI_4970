### epiAnalysisWindLowBirthWeight_HEI_4970.R ###
# Author: Andrew Larkin
# Date created: Feb 3, 2023
# Summary: This R file contains the epidemiological analysis comparing low term birth weights between neighbors near
# highways who are matched based on whether they are upwind/downwind of the highway. The main analyses are generalized
# linear regression models.  Sensitivity analyses are stratified by demographics (e.g. race) or exposure 
# (e.g. distance to nearest road)

# for details about cohort records and vital statistics, see "Data Dictionary - HEI Birth Data - 1996 to 2016 V1.docx"
# for details about deriving the wind-based exposure metrics, see "Deriving Wind Exposure Metrics_HEI_4970.docx"
# for details about matching upwinwd/downwind neighbors, see "Matching by Wind_HEI_4970.docx"


############# STEP 1: SETUP ############# 

setwd("C:/users/larki_9x8fs5f/Desktop/TexasHEIManuscript") # folder where are scripts are stored
source("epiFunctionsForWindStudy.R") # contains functions for running linear and logistic models
resultsFolder <- paste(epiFolder, "results/lowBirthWeight/", sep="") # folderpath where epi results will be stored

### load datasets ###
### restricted to 37 to 42 weeks (for term birth weight and low birth weight outcomes) ###

# datasets with single matches
epi_37_42_1_15 <- loadMatchData(15,multipleMatches = FALSE, restricted = TRUE) #15m, restricted, single matches
epi_37_42_1_25 <- loadMatchData(25, FALSE, TRUE) #25m, restricted, single matches
epi_37_42_1_50 <- loadMatchData(50, FALSE, TRUE) #50m, restricted, single matches
epi_37_42_1_100 <- loadMatchData(100, FALSE, TRUE) #100m, restricted, single matches

# datasets with multiple matches
epi_37_42_4_15 <- loadMatchData(15, multipleMatches = TRUE, restricted= TRUE) #15m, restricted, multiple matches
epi_37_42_4_25 <- loadMatchData(25, TRUE, TRUE) #25m, restricted, multiple matches
epi_37_42_4_50 <- loadMatchData(50, TRUE, TRUE) #50m, restricted, multiple matches
epi_37_42_4_100 <- loadMatchData(100, TRUE, TRUE) #100m, restricted, multiple matches


############# STEP 2: MATCH CRITERIA ANALYSIS ############# 

# single match, base models
lowbirth_0_1_15 <- getWindCoefficientsLogistic(15,"ltbw",resample=FALSE,fullModel=FALSE,restricted=TRUE)
lowbirth_0_1_25 <- getWindCoefficientsLogistic(25,"ltbw",resample=FALSE,fullModel=FALSE,restricted=TRUE)
lowbirth_0_1_50 <- getWindCoefficientsLogistic(50,"ltbw",resample=FALSE,fullModel=FALSE,restricted=TRUE)
lowbirth_0_1_100 <- getWindCoefficientsLogistic(100,"ltbw",resample=FALSE,fullModel=FALSE,restricted=TRUE)

# multiple match, base models
lowbirth_0_4_15 <- getWindCoefficientsLogistic(15,"ltbw",resample=TRUE,fullModel=FALSE,restricted=TRUE)
lowbirth_0_4_25 <- getWindCoefficientsLogistic(25,"ltbw",resample=TRUE,fullModel=FALSE,restricted=TRUE)
lowbirth_0_4_50 <- getWindCoefficientsLogistic(50,"ltbw",resample=TRUE,fullModel=FALSE,restricted=TRUE)
lowbirth_0_4_100 <- getWindCoefficientsLogistic(100,"ltbw",resample=TRUE,fullModel=FALSE,restricted=TRUE)

# single match, full models
lowbirth_1_1_15 <- getWindCoefficientsLogistic(15,"ltbw",resample=FALSE,fullModel=TRUE,restricted=TRUE)
lowbirth_1_1_25 <- getWindCoefficientsLogistic(25,"ltbw",resample=FALSE,fullModel=TRUE,restricted=TRUE)
lowbirth_1_1_50 <- getWindCoefficientsLogistic(50,"ltbw",resample=FALSE,fullModel=TRUE,restricted=TRUE)
lowbirth_1_1_100 <- getWindCoefficientsLogistic(100,"ltbw",resample=FALSE,fullModel=TRUE,restricted=TRUE)

# multiple match, full models
lowbirth_1_4_15 <- getWindCoefficientsLogistic(15,"ltbw",resample=TRUE,fullModel=TRUE,restricted=TRUE)
lowbirth_1_4_25 <- getWindCoefficientsLogistic(25,"ltbw",resample=TRUE,fullModel=TRUE,restricted=TRUE)
lowbirth_1_4_50 <- getWindCoefficientsLogistic(50,"ltbw",resample=TRUE,fullModel=TRUE,restricted=TRUE)
lowbirth_1_4_100 <- getWindCoefficientsLogistic(100,"ltbw",resample=TRUE,fullModel=TRUE,restricted=TRUE)


# combine match criteria results into a single dataset
lowbirthWeightMatchAnalysis <- rbind(
  
  # single match, base models
  lowbirth_0_1_15,
  lowbirth_0_1_25,
  lowbirth_0_1_50,
  lowbirth_0_1_100,
  
  # multiple match, base models
  lowbirth_0_4_15,
  lowbirth_0_4_25,
  lowbirth_0_4_50,
  lowbirth_0_4_100,
  
  # single match, full models
  lowbirth_1_1_15,
  lowbirth_1_1_25,
  lowbirth_1_1_50,
  lowbirth_1_1_100,
  
  # multiple match, full models
  lowbirth_1_4_15,
  lowbirth_1_4_25,
  lowbirth_1_4_50,
  lowbirth_1_4_100
)

colnames(lowbirthWeightMatchAnalysis) <- c("n_pairs","n_events","modelType","multipleMatch","matchCriteria","exp(coef)","exp(-coef)","lowerCI","upperCI")

# write match criteria results to csv
write.csv(lowbirthWeightMatchAnalysis,paste(resultsFolder,"lowbirthWeightMatchAnalysis.csv",sep=""),row.names = FALSE)



############# STEP 3: DISTANCE TO NEAREST ROAD ANALYSIS ############# 

# run models stratified by distance to nearest road
lowbirth_0_50 <- getWindCoefficientsDistanceLogistic("epi_37_42_4_100","ltbw",0,50)
lowbirth_50_100 <- getWindCoefficientsDistanceLogistic("epi_37_42_4_100","ltbw",50,100)
lowbirth_100_300 <- getWindCoefficientsDistanceLogistic("epi_37_42_4_100","ltbw",100,300)
lowbirth_300_400 <- getWindCoefficientsDistanceLogistic("epi_37_42_4_100","ltbw",300,400)
lowbirth_400_500 <- getWindCoefficientsDistanceLogistic("epi_37_42_4_100","ltbw",400,500)

# combine stratified distance analyses into a single dataset
lowbirthwt_distance <- rbind(
  lowbirth_0_50,
  lowbirth_50_100,
  lowbirth_100_300,
  lowbirth_300_400,
  lowbirth_400_500
)

colnames(lowbirthwt_distance) <- c("n_pairs","n_events","lower threshold","upper threshold", "exp(coef)","exp(-coef)","lowerCI","upperCI")

# write distance to road results to csv
write.csv(lowbirthwt_distance,paste(resultsFolder,"lowbirthWeightDistanceAnalysis.csv",sep=""),row.names = FALSE)




############# STEP 4: THREE YEAR ROLLNG AVERAGE ANALYSIS ############# 

# run three year rolling average models from 2007 to 2016
lowbirth_07_09 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2007,2009)
lowbirth_08_10 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2008,2010)
lowbirth_09_11 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2009,2011)
lowbirth_10_12 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2010,2012)
lowbirth_11_13 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2011,2013)
lowbirth_12_14 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2012,2014)
lowbirth_13_15 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2013,2015)
lowbirth_14_16 <- getWindCoefficientsTimeLogistic("epi_37_42_4_100","ltbw",2014,2016)

# combine three year rolling average results into a single dataset
lowbirthwt_time <- rbind(
  lowbirth_07_09,
  lowbirth_08_10,
  lowbirth_09_11,
  lowbirth_10_12,
  lowbirth_11_13,
  lowbirth_12_14,
  lowbirth_13_15,
  lowbirth_14_16
)

colnames(lowbirthwt_time) <- c("n_pairs","n_events","first year","last year", "exp(coef)","exp(-coef)","lowerCI","upperCI")

# write three year rolling average results to csv
write.csv(lowbirthwt_time,paste(resultsFolder,"lowbirthWeightTimeAnalysis.csv",sep=""),row.names = FALSE)



############# STEP 5: BUILDING AND TREE SHIELDING ANALYSIS ############# 


lowbirth_shield_0_50 <- getWindCoefficientsShieldingDistanceLogistic("epi_37_42_4_100",'ltbw',0,50)
lowbirth_shield_50_100 <- getWindCoefficientsShieldingDistanceLogistic("epi_37_42_4_100",'ltbw',50,100)
lowbirth_shield_100_300 <- getWindCoefficientsShieldingDistanceLogistic("epi_37_42_4_100",'ltbw',100,300)
lowbirth_shield_300_400 <- getWindCoefficientsShieldingDistanceLogistic("epi_37_42_4_100",'ltbw',300,400)
lowbirth_shield_400_500 <- getWindCoefficientsShieldingDistanceLogistic("epi_37_42_4_100",'ltbw',400,500)
lowbirth_shield_0_500 <- getWindCoefficientsShieldingDistanceLogistic("epi_37_42_4_100",'ltbw',0,600)

lowbirthwt_shield <- rbind(
  lowbirth_shield_0_50,
  lowbirth_shield_50_100,
  lowbirth_shield_100_300,
  lowbirth_shield_300_400,
  lowbirth_shield_400_500,
  lowbirth_shield_0_500
)

colnames(lowbirthwt_shield) <- c("n_pairs","n_events","lower dist","upper dist","exp(coef)","exp(-coef)","lowerCI","upperCI")
write.csv(lowbirthwt_shield,paste(resultsFolder,"lowbirthWeightShieldingAnalysis.csv",sep=""),row.names = FALSE)


############# STEP 6: CONTINUOUS WIND ANALYSIS ############# 

lowbirth_percWind_0_50 <- getPercentWindCoefficientsDistanceLogistic("epi_37_42_4_100",'ltbw',0,50)
lowbirth_percWind_50_100 <- getPercentWindCoefficientsDistanceLogistic("epi_37_42_4_100",'ltbw',50,100)
lowbirth_percWind_100_300 <- getPercentWindCoefficientsDistanceLogistic("epi_37_42_4_100",'ltbw',100,300)
lowbirth_percWind_300_400 <- getPercentWindCoefficientsDistanceLogistic("epi_37_42_4_100",'ltbw',300,400)
lowbirth_percWind_400_500 <- getPercentWindCoefficientsDistanceLogistic("epi_37_42_4_100",'ltbw',400,500)
lowbirth_percWind_0_500 <- getPercentWindCoefficientsDistanceLogistic("epi_37_42_4_100",'ltbw',0,500)

lowbirthwt_percWind <- rbind(
  lowbirth_percWind_0_50,
  lowbirth_percWind_50_100,
  lowbirth_percWind_100_300,
  lowbirth_percWind_300_400,
  lowbirth_percWind_400_500,
  lowbirth_percWind_0_500
)

colnames(lowbirthwt_percWind) <-  c("n_pairs","n_events","lower dist","upper dist","exp(coef)","exp(-coef)","lowerCI","upperCI")
write.csv(lowbirthwt_percWind,paste(resultsFolder,"lowbirthWeightPercentWindAnalysis.csv",sep=""),row.names = FALSE)


############# STEP 7: DEMOGRAPHIC SENSITIVITY MODELS ############# 

# restrict to black non-hispanic
blackSubset <- getRaceSubset("epi_37_42_4_100",2)
blackNonHisp <- getEthnicitySubset("blackSubset",0)
blackNonHisp_4_100 <- runRaceSensitivityModelLogistic("blackNonHisp","ltbw")

# restrict to white non-hispanic
whiteSubset <- getRaceSubset("epi_37_42_4_100",1)
whiteNonHisp <- getEthnicitySubset("whiteSubset",0)
whiteNonHisp_4_100 <- runRaceSensitivityModelLogistic("whiteNonHisp","ltbw")

# restrict to hispanic or latina
hispanic <- getEthnicitySubset("epi_37_42_4_100",0)
hisp_4_100 <- runEthnicitySensitivityModelLogistic("hispanic","ltbw")

# restrict to high school education or less
lowerEduc <- getLowerEducSubset("epi_37_42_4_100")
lowerEduc_4_100 <- runEducSensitivityModelLogistic("lowerEduc","ltbw")

# restrict to college diploma or greater
higherEduc <- getHigherEducSubset("epi_37_42_4_100")
higherEduc_4_100 <- runEducSensitivityModelLogistic("higherEduc","ltbw")

# restrict to US born mother
usBorn <- getForeignCat("epi_37_42_4_100",0)
usBorn_4_100 <- runForeignSensitivityModelLogistic("usBorn","ltbw")

# restrict to foreign born mother
foreignBorn <- getForeignCat("epi_37_42_4_100",1)
foreignBorn_4_100 <- runForeignSensitivityModelLogistic("foreignBorn","ltbw")

# restrict to lowest income neighborhood tertile
lowerIncome <- getIncomeCat("epi_37_42_4_100",1)
lowerIncome_4_100 <- runIncomeSensitivityModelLogistic("lowerIncome","ltbw")

# restrict to highest income neighborhood tertile
higherIncome <- getIncomeCat("epi_37_42_4_100",3)
higherIncome_4_100 <- runIncomeSensitivityModelLogistic("higherIncome","ltbw")


# combine stratifications into single dataset
demographicStrata <- rbind(
  blackNonHisp_4_100,
  whiteNonHisp_4_100,
  hisp_4_100,
  lowerEduc_4_100,
  higherEduc_4_100,
  usBorn_4_100,
  foreignBorn_4_100,
  lowerIncome_4_100,
  higherIncome_4_100
)

# create a column of labels for each stratification label
demographicCategories <- c(
  "black non-hispanic",
  "white non-hispanic",
  "hispanic or latina",
  "high school graduate or less",
  "greater than high school graduate",
  "US born",
  "foreign born",
  "low income neighborhood",
  "high income neighborhood"
)
demographicStrata <- cbind(demographicCategories,demographicStrata)

colnames(demographicStrata) <- c("category","n_pairs","n_events","exp(coef)","exp(-coef)","lowerCI","upperCI")

# write stratification results to csv
write.csv(demographicStrata,paste(resultsFolder,"lowbirthWeightDemographicSensitivityAnalysis.csv",sep=""),row.names = FALSE)





#############  end of epiAnalysisLowWindBirthWeight_HEI_4970.R ############# 