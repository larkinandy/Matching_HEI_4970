### epiAnalysisWindVeryPrematureBirth_HEI_4970.R ###
# Author: Andrew Larkin
# Date created: Feb 3, 2023
# Summary: This R file contains the epidemiological analysis comparing very premature birth events birth weights between neighbors near
# highways who are matched based on whether they are upwind/downwind of the highway. The main analyses are generalized
# linear regression models.  Sensitivity analyses are stratified by demographics (e.g. race) or exposure 
# (e.g. distance to nearest road)

# for details about cohort records and vital statistics, see "Data Dictionary - HEI Birth Data - 1996 to 2016 V1.docx"
# for details about deriving the wind-based exposure metrics, see "Deriving Wind Exposure Metrics_HEI_4970.docx"
# for details about matching upwinwd/downwind neighbors, see "Matching by Wind_HEI_4970.docx"


############# STEP 1: SETUP ############# 

setwd("C:/users/larki_9x8fs5f/Desktop/TexasHEIManuscript") # folder where are scripts are stored
source("epiFunctionsForWindStudy.R") # contains functions for running linear and logistic models
resultsFolder <- paste(epiFolder, "results/veryPrematureBirth/", sep="") # folderpath where epi results will be stored

### load datasets ###
### restricted to 37 to 42 weeks (for term birth weight and low birth weight outcomes) ###

# datasets with single matches
epi_all_1_15 <- loadMatchData(15,multipleMatches = FALSE, restricted = FALSE) #15m, restricted, single matches
epi_all_1_25 <- loadMatchData(25, FALSE, FALSE) #25m, restricted, single matches
epi_all_1_50 <- loadMatchData(50, FALSE, FALSE) #50m, restricted, single matches
epi_all_1_100 <- loadMatchData(100, FALSE, FALSE) #100m, restricted, single matches

# datasets with multiple matches
epi_all_4_15 <- loadMatchData(15, multipleMatches = TRUE, restricted= FALSE) #15m, restricted, multiple matches
epi_all_4_25 <- loadMatchData(25, TRUE, FALSE) #25m, restricted, multiple matches
epi_all_4_50 <- loadMatchData(50, TRUE, FALSE) #50m, restricted, multiple matches
epi_all_4_100 <- loadMatchData(100, TRUE, FALSE) #100m, restricted, multiple matches


############# STEP 2: MATCH CRITERIA ANALYSIS ############# 

# single match, base models
vptb_0_1_15 <- getWindCoefficientsLogistic(15,"vptb",resample=FALSE,fullModel=FALSE,restricted=FALSE)
vptb_0_1_25 <- getWindCoefficientsLogistic(25,"vptb",resample=FALSE,fullModel=FALSE,restricted=FALSE)
vptb_0_1_50 <- getWindCoefficientsLogistic(50,"vptb",resample=FALSE,fullModel=FALSE,restricted=FALSE)
vptb_0_1_100 <- getWindCoefficientsLogistic(100,"vptb",resample=FALSE,fullModel=FALSE,restricted=FALSE)

# multiple match, base models
vptb_0_4_15 <- getWindCoefficientsLogistic(15,"vptb",resample=TRUE,fullModel=FALSE,restricted=FALSE)
vptb_0_4_25 <- getWindCoefficientsLogistic(25,"vptb",resample=TRUE,fullModel=FALSE,restricted=FALSE)
vptb_0_4_50 <- getWindCoefficientsLogistic(50,"vptb",resample=TRUE,fullModel=FALSE,restricted=FALSE)
vptb_0_4_100 <- getWindCoefficientsLogistic(100,"vptb",resample=TRUE,fullModel=FALSE,restricted=FALSE)

# single match, full models
#vptb_1_1_15 <- getWindCoefficientsLogistic(15,"vptb",resample=FALSE,fullModel=TRUE,restricted=FALSE) # don't run - this model doesn't converge
vptb_1_1_25 <- getWindCoefficientsLogistic(25,"vptb",resample=FALSE,fullModel=TRUE,restricted=FALSE)
vptb_1_1_50 <- getWindCoefficientsLogistic(50,"vptb",resample=FALSE,fullModel=TRUE,restricted=FALSE)
vptb_1_1_100 <- getWindCoefficientsLogistic(100,"vptb",resample=FALSE,fullModel=TRUE,restricted=FALSE)

# multiple match, full models
vptb_1_4_15 <- getWindCoefficientsLogistic(15,"vptb",resample=TRUE,fullModel=TRUE,restricted=FALSE)
vptb_1_4_25 <- getWindCoefficientsLogistic(25,"vptb",resample=TRUE,fullModel=TRUE,restricted=FALSE)
vptb_1_4_50 <- getWindCoefficientsLogistic(50,"vptb",resample=TRUE,fullModel=TRUE,restricted=FALSE)
vptb_1_4_100 <- getWindCoefficientsLogistic(100,"vptb",resample=TRUE,fullModel=TRUE,restricted=FALSE)


# combine match criteria results into a single dataset
vptMatchAnalysis <- rbind(
  
  # single match, base models
  vptb_0_1_15,
  vptb_0_1_25,
  vptb_0_1_50,
  vptb_0_1_100,
  
  # multiple match, base models
  vptb_0_4_15,
  vptb_0_4_25,
  vptb_0_4_50,
  vptb_0_4_100,
  
  # single match, full models
  #vptb_1_1_15, don't save - this model doesn't converge
  vptb_1_1_25,
  vptb_1_1_50,
  vptb_1_1_100,
  
  # multiple match, full models
  vptb_1_4_15,
  vptb_1_4_25,
  vptb_1_4_50,
  vptb_1_4_100
)

colnames(vptMatchAnalysis) <- c("n_pairs","n_events","modelType","multipleMatch","matchCriteria","exp(coef)","exp(-coef)","lowerCI","upperCI")

# write match criteria results to csv
write.csv(vptMatchAnalysis,paste(resultsFolder,"vptbMatchAnalysis.csv",sep=""),row.names = FALSE)



############# STEP 3: DISTANCE TO NEAREST ROAD ANALYSIS ############# 

# run models stratified by distance to nearest road
#vptb_0_50 <- getWindCoefficientsDistanceLogistic("epi_all_4_100","vptb",0,50) # dont' run - this model doesn't converge
vptb_50_100 <- getWindCoefficientsDistanceLogistic("epi_all_4_100","vptb",50,100)
vptb_100_300 <- getWindCoefficientsDistanceLogistic("epi_all_4_100","vptb",100,300)
vptb_300_400 <- getWindCoefficientsDistanceLogistic("epi_all_4_100","vptb",300,400)
vptb_400_500 <- getWindCoefficientsDistanceLogistic("epi_all_4_100","vptb",400,500)

# combine stratified distance analyses into a single dataset
vptb_distance <- rbind(
  #ptb_0_50, # don't save - this model doesn't converge
  vptb_50_100,
  vptb_100_300,
  vptb_300_400,
  vptb_400_500
)

colnames(vptb_distance) <- c("n_pairs","n_events","lower threshold","upper threshold", "exp(coef)","exp(-coef)","lowerCI","upperCI")

# write distance to road results to csv
write.csv(vptb_distance,paste(resultsFolder,"vptbDistanceAnalysis.csv",sep=""),row.names = FALSE)




############# STEP 4: THREE YEAR ROLLNG AVERAGE ANALYSIS ############# 

# run three year rolling average models from 2007 to 2016
vptb_07_09 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2007,2009)
vptb_08_10 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2008,2010)
vptb_09_11 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2009,2011)
vptb_10_12 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2010,2012)
vptb_11_13 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2011,2013)
vptb_12_14 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2012,2014)
vptb_13_15 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2013,2015)
vptb_14_16 <- getWindCoefficientsTimeLogistic("epi_all_4_100","vptb",2014,2016)

# combine three year rolling average results into a single dataset
vptb_time <- rbind(
  vptb_07_09,
  vptb_08_10,
  vptb_09_11,
  vptb_10_12,
  vptb_11_13,
  vptb_12_14,
  vptb_13_15,
  vptb_14_16
)

colnames(vptb_time) <- c("n_pairs","n_events","first year","last year", "exp(coef)","exp(-coef)","lowerCI","upperCI")

# write three year rolling average results to csv
write.csv(vptb_time,paste(resultsFolder,"vptbTimeAnalysis.csv",sep=""),row.names = FALSE)



############# STEP 5: BUILDING AND TREE SHIELDING ANALYSIS ############# 


# vptb_shield_0_50 <- getWindCoefficientsShieldingDistanceLogistic("epi_all_4_100",'vptb',0,50) don't run - this model doesn't converge
vptb_shield_50_100 <- getWindCoefficientsShieldingDistanceLogistic("epi_all_4_100",'vptb',50,100)
vptb_shield_100_300 <- getWindCoefficientsShieldingDistanceLogistic("epi_all_4_100",'vptb',100,300)
vptb_shield_300_400 <- getWindCoefficientsShieldingDistanceLogistic("epi_all_4_100",'vptb',300,400)
vptb_shield_400_500 <- getWindCoefficientsShieldingDistanceLogistic("epi_all_4_100",'vptb',400,500)
vptb_shield_0_500 <- getWindCoefficientsShieldingDistanceLogistic("epi_all_4_100",'vptb',0,600)

vptb_shield <- rbind(
  #vptb_shield_0_50, # don't save - this model doesn't converge
  vptb_shield_50_100,
  vptb_shield_100_300,
  vptb_shield_300_400,
  vptb_shield_400_500,
  vptb_shield_0_500
)

colnames(vptb_shield) <- c("n_pairs","n_events","lower dist","upper dist","exp(coef)","exp(-coef)","lowerCI","upperCI")
write.csv(vptb_shield,paste(resultsFolder,"vptbShieldingAnalysis.csv",sep=""),row.names = FALSE)


############# STEP 6: CONTINUOUS WIND ANALYSIS ############# 

# vptb_percWind_0_50 <- getPercentWindCoefficientsDistanceLogistic("epi_all_4_100",'vptb',0,50) don't run - this model doesn't converge
vptb_percWind_50_100 <- getPercentWindCoefficientsDistanceLogistic("epi_all_4_100",'vptb',50,100)
vptb_percWind_100_300 <- getPercentWindCoefficientsDistanceLogistic("epi_all_4_100",'vptb',100,300)
vptb_percWind_300_400 <- getPercentWindCoefficientsDistanceLogistic("epi_all_4_100",'vptb',300,400)
vptb_percWind_400_500 <- getPercentWindCoefficientsDistanceLogistic("epi_all_4_100",'vptb',400,500)
vptb_percWind_0_500 <- getPercentWindCoefficientsDistanceLogistic("epi_all_4_100",'vptb',0,500)

vptb_percWind <- rbind( 
  #vptb_percWind_0_50, # don't save - this model doesn't converge
  vptb_percWind_50_100,
  vptb_percWind_100_300,
  vptb_percWind_300_400,
  vptb_percWind_400_500,
  vptb_percWind_0_500
)

colnames(vptb_percWind) <-  c("n_pairs","n_events","lower dist","upper dist","exp(coef)","exp(-coef)","lowerCI","upperCI")
write.csv(vptb_percWind,paste(resultsFolder,"vptbPercentWindAnalysis.csv",sep=""),row.names = FALSE)


############# STEP 7: DEMOGRAPHIC SENSITIVITY MODELS ############# 

# restrict to black non-hispanic
blackSubset <- getRaceSubset("epi_all_4_100",2)
blackNonHisp <- getEthnicitySubset("blackSubset",0)
blackNonHisp_4_100 <- runRaceSensitivityModelLogistic("blackNonHisp","vptb")

# restrict to white non-hispanic
whiteSubset <- getRaceSubset("epi_all_4_100",1)
whiteNonHisp <- getEthnicitySubset("whiteSubset",0)
whiteNonHisp_4_100 <- runRaceSensitivityModelLogistic("whiteNonHisp","vptb")

# restrict to hispanic or latina
hispanic <- getEthnicitySubset("epi_all_4_100",0)
hisp_4_100 <- runEthnicitySensitivityModelLogistic("hispanic","vptb")

# restrict to high school education or less
lowerEduc <- getLowerEducSubset("epi_all_4_100")
lowerEduc_4_100 <- runEducSensitivityModelLogistic("lowerEduc","vptb")

# restrict to college diploma or greater
higherEduc <- getHigherEducSubset("epi_all_4_100")
higherEduc_4_100 <- runEducSensitivityModelLogistic("higherEduc","vptb")

# restrict to US born mother
usBorn <- getForeignCat("epi_all_4_100",0)
usBorn_4_100 <- runForeignSensitivityModelLogistic("usBorn","vptb")

# restrict to foreign born mother
foreignBorn <- getForeignCat("epi_all_4_100",1)
foreignBorn_4_100 <- runForeignSensitivityModelLogistic("foreignBorn","vptb")

# restrict to lowest income neighborhood tertile
lowerIncome <- getIncomeCat("epi_all_4_100",1)
lowerIncome_4_100 <- runIncomeSensitivityModelLogistic("lowerIncome","vptb")

# restrict to highest income neighborhood tertile
higherIncome <- getIncomeCat("epi_all_4_100",3)
higherIncome_4_100 <- runIncomeSensitivityModelLogistic("higherIncome","vptb")


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
write.csv(demographicStrata,paste(resultsFolder,"vptbDemographicSensitivityAnalysis.csv",sep=""),row.names = FALSE)





#############  end of epiAnalysisVeryPrematureBirthWeight_HEI_4970.R ############# 