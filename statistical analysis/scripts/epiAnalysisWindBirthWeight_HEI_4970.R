### epiAnalysisWindBirthWeight_HEI_4970.R ###
# Author: Andrew Larkin
# Date created: Feb 3, 2023
# Summary: This R file contains the epidemiological analysis comparing term birth weights between neighbors near
# highways who are matched based on whether they are upwind/downwind of the highway. The main analyses are generalized
# linear regression models.  Sensitivity analyses are stratified by demographics (e.g. race) or exposure 
# (e.g. distance to nearest road)

# for details about cohort records and vital statistics, see "Data Dictionary - HEI Birth Data - 1996 to 2016 V1.docx"
# for details about deriving the wind-based exposure metrics, see "Deriving Wind Exposure Metrics_HEI_4970.docx"
# for details about matching upwinwd/downwind neighbors, see "Matching by Wind_HEI_4970.docx"


############# STEP 1: SETUP ############# 

setwd("C:/users/larki_9x8fs5f/Desktop/TexasHEIManuscript") # folder where are scripts are stored
source("epiFunctionsForWindStudy.R") # contains functions for running linear and logistic models
resultsFolder <- paste(epiFolder, "results/birthWeight/", sep="") # folderpath where epi results will be stored

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
birth_0_1_15 <- getWindCoefficientsLinear(15,"b_wt_cgr",resample=FALSE,fullModel=FALSE)
birth_0_1_25 <- getWindCoefficientsLinear(25,"b_wt_cgr",resample=FALSE,fullModel=FALSE)
birth_0_1_50 <- getWindCoefficientsLinear(50,"b_wt_cgr",resample=FALSE,fullModel=FALSE)
birth_0_1_100 <- getWindCoefficientsLinear(100,"b_wt_cgr",resample=FALSE,fullModel=FALSE)

# multiple match, base models
birth_0_4_15 <- getWindCoefficientsLinear(15,"b_wt_cgr",resample=TRUE,fullModel=FALSE)
birth_0_4_25 <- getWindCoefficientsLinear(25,"b_wt_cgr",resample=TRUE,fullModel=FALSE)
birth_0_4_50 <- getWindCoefficientsLinear(50,"b_wt_cgr",resample=TRUE,fullModel=FALSE)
birth_0_4_100 <- getWindCoefficientsLinear(100,"b_wt_cgr",resample=TRUE,fullModel=FALSE)

# single match, full models
birth_1_1_15 <- getWindCoefficientsLinear(15,"b_wt_cgr",resample=FALSE,fullModel=TRUE)
birth_1_1_25 <- getWindCoefficientsLinear(25,"b_wt_cgr",resample=FALSE,fullModel=TRUE)
birth_1_1_50 <- getWindCoefficientsLinear(50,"b_wt_cgr",resample=FALSE,fullModel=TRUE)
birth_1_1_100 <- getWindCoefficientsLinear(100,"b_wt_cgr",resample=FALSE,fullModel=TRUE)

# multiple match, full models
birth_1_4_15 <- getWindCoefficientsLinear(15,"b_wt_cgr",resample=TRUE,fullModel=TRUE)
birth_1_4_25 <- getWindCoefficientsLinear(25,"b_wt_cgr",resample=TRUE,fullModel=TRUE)
birth_1_4_50 <- getWindCoefficientsLinear(50,"b_wt_cgr",resample=TRUE,fullModel=TRUE)
birth_1_4_100 <- getWindCoefficientsLinear(100,"b_wt_cgr",resample=TRUE,fullModel=TRUE)


# combine match criteria results into a single dataset
birthWeightMatchAnalysis <- rbind(
  
  # single match, base models
  birth_0_1_15,
  birth_0_1_25,
  birth_0_1_50,
  birth_0_1_100,
  
  # multiple match, base models
  birth_0_4_15,
  birth_0_4_25,
  birth_0_4_50,
  birth_0_4_100,
  
  # single match, full models
  birth_1_1_15,
  birth_1_1_25,
  birth_1_1_50,
  birth_1_1_100,
  
  # multiple match, full models
  birth_1_4_15,
  birth_1_4_25,
  birth_1_4_50,
  birth_1_4_100
)

colnames(birthWeightMatchAnalysis) <- c("n_pairs","modelType","multipleMatch","matchCriteria","lowerCI","coefficient","upperCI")

# write match criteria results to csv
write.csv(birthWeightMatchAnalysis,paste(resultsFolder,"birthWeightMatchAnalysis.csv",sep=""),row.names = FALSE)



############# STEP 3: DISTANCE TO NEAREST ROAD ANALYSIS ############# 

# run models stratified by distance to nearest road
birth_0_50 <- getWindCoefficientsDistanceLinear("epi_37_42_4_100","b_wt_cgr",0,50)
birth_50_100 <- getWindCoefficientsDistanceLinear("epi_37_42_4_100","b_wt_cgr",50,100)
birth_100_300 <- getWindCoefficientsDistanceLinear("epi_37_42_4_100","b_wt_cgr",100,300)
birth_300_400 <- getWindCoefficientsDistanceLinear("epi_37_42_4_100","b_wt_cgr",300,400)
birth_400_500 <- getWindCoefficientsDistanceLinear("epi_37_42_4_100","b_wt_cgr",400,500)

# combine stratified distance analyses into a single dataset
birthwt_distance <- rbind(
  birth_0_50,
  birth_50_100,
  birth_100_300,
  birth_300_400,
  birth_400_500
)

colnames(birthwt_distance) <- c("n_pairs","lower threshold","upper threshold", "lowerCI","coefficient","upperCI")

# write distance to road results to csv
write.csv(birthwt_distance,paste(resultsFolder,"birthWeightDistanceAnalysis.csv",sep=""),row.names = FALSE)




############# STEP 4: THREE YEAR ROLLNG AVERAGE ANALYSIS ############# 

# run three year rolling average models from 2007 to 2016
birth_07_09 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2007,2009)
birth_08_10 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2008,2010)
birth_09_11 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2009,2011)
birth_10_12 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2010,2012)
birth_11_13 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2011,2013)
birth_12_14 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2012,2014)
birth_13_15 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2013,2015)
birth_14_16 <- getWindCoefficientsTimeLinear("epi_37_42_4_100","b_wt_cgr",2014,2016)

# combine three year rolling average results into a single dataset
birthwt_time <- rbind(
  birth_07_09,
  birth_08_10,
  birth_09_11,
  birth_10_12,
  birth_11_13,
  birth_12_14,
  birth_13_15,
  birth_14_16
)

colnames(birthwt_time) <- c("n_pairs","first year","last year", "lowerCI","coefficient","upperCI")

# write three year rolling average results to csv
write.csv(birthwt_time,paste(resultsFolder,"birthWeightTimeAnalysis.csv",sep=""),row.names = FALSE)



############# STEP 5: BUILDING AND TREE SHIELDING ANALYSIS ############# 


birth_shield_0_50 <- getWindCoefficientsShieldingDistanceLinear("epi_37_42_4_100",'b_wt_cgr',0,50)
birth_shield_50_100 <- getWindCoefficientsShieldingDistanceLinear("epi_37_42_4_100",'b_wt_cgr',50,100)
birth_shield_100_300 <- getWindCoefficientsShieldingDistanceLinear("epi_37_42_4_100",'b_wt_cgr',100,300)
birth_shield_300_400 <- getWindCoefficientsShieldingDistanceLinear("epi_37_42_4_100",'b_wt_cgr',300,400)
birth_shield_400_500 <- getWindCoefficientsShieldingDistanceLinear("epi_37_42_4_100",'b_wt_cgr',400,500)
birth_shield_0_500 <- getWindCoefficientsShieldingDistanceLinear("epi_37_42_4_100",'b_wt_cgr',0,500)

birthwt_shield <- rbind(
  birth_shield_0_50,
  birth_shield_50_100,
  birth_shield_100_300,
  birth_shield_300_400,
  birth_shield_400_500,
  birth_shield_0_500
)

colnames(birthwt_shield) <- c("n_pairs","lower dist","upper dist", "lowerCI","coefficient","upperCI")
write.csv(birthwt_shield,paste(resultsFolder,"birthWeightShieldingAnalysis.csv",sep=""),row.names = FALSE)


############# STEP 6: CONTINUOUS WIND ANALYSIS ############# 

birth_percWind_0_50 <- getPercentWindCoefficientsDistanceLinear("epi_37_42_4_100",'b_wt_cgr',0,50)
birth_percWind_50_100 <- getPercentWindCoefficientsDistanceLinear("epi_37_42_4_100",'b_wt_cgr',50,100)
birth_percWind_100_300 <- getPercentWindCoefficientsDistanceLinear("epi_37_42_4_100",'b_wt_cgr',100,300)
birth_percWind_300_400 <- getPercentWindCoefficientsDistanceLinear("epi_37_42_4_100",'b_wt_cgr',300,400)
birth_percWind_400_500 <- getPercentWindCoefficientsDistanceLinear("epi_37_42_4_100",'b_wt_cgr',400,500)
birth_percWind_0_500 <- getPercentWindCoefficientsDistanceLinear("epi_37_42_4_100",'b_wt_cgr',0,500)

birthwt_percWind <- rbind(
  birth_percWind_0_50,
  birth_percWind_50_100,
  birth_percWind_100_300,
  birth_percWind_300_400,
  birth_percWind_400_500,
  birth_percWind_0_500
)

colnames(birthwt_percWind) <- c("n_pairs","lower dist","upper dist", "lowerCI","coefficient","upperCI")
write.csv(birthwt_percWind,paste(resultsFolder,"birthWeightPercentWindAnalysis.csv",sep=""),row.names = FALSE)


############# STEP 7: DEMOGRAPHIC SENSITIVITY MODELS ############# 

# restrict to black non-hispanic
blackSubset <- getRaceSubset("epi_37_42_4_100",2)
blackNonHisp <- getEthnicitySubset("blackSubset",0)
blackNonHisp_4_100 <- runRaceSensitivityModelLinear("blackNonHisp","b_wt_cgr")

# restrict to white non-hispanic
whiteSubset <- getRaceSubset("epi_37_42_4_100",1)
whiteNonHisp <- getEthnicitySubset("whiteSubset",0)
whiteNonHisp_4_100 <- runRaceSensitivityModelLinear("whiteNonHisp","b_wt_cgr")

# restrict to hispanic or latina
hispanic <- getEthnicitySubset("epi_37_42_4_100",0)
hisp_4_100 <- runEthnicitySensitivityModelLinear("hispanic","b_wt_cgr")

# restrict to high school education or less
lowerEduc <- getLowerEducSubset("epi_37_42_4_100")
lowerEduc_4_100 <- runEducSensitivityModelLinear("lowerEduc","b_wt_cgr")

# restrict to college diploma or greater
higherEduc <- getHigherEducSubset("epi_37_42_4_100")
higherEduc_4_100 <- runEducSensitivityModelLinear("higherEduc","b_wt_cgr")

# restrict to US born mother
usBorn <- getForeignCat("epi_37_42_4_100",0)
usBorn_4_100 <- runForeignSensitivityModelLinear("usBorn","b_wt_cgr")

# restrict to foreign born mother
foreignBorn <- getForeignCat("epi_37_42_4_100",1)
foreignBorn_4_100 <- runForeignSensitivityModelLinear("foreignBorn","b_wt_cgr")

# restrict to lowest income neighborhood tertile
lowerIncome <- getIncomeCat("epi_37_42_4_100",1)
lowerIncome_4_100 <- runIncSensitivityModelLinear("lowerIncome","b_wt_cgr")

# restrict to highest income neighborhood tertile
higherIncome <- getIncomeCat("epi_37_42_4_100",3)
higherIncome_4_100 <- runIncSensitivityModelLinear("higherIncome","b_wt_cgr")


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

colnames(demographicStrata) <- c("category","n_pairs","lowerCI","coefficient","upperCI")

# write stratification results to csv
write.csv(demographicStrata,paste(resultsFolder,"birthWeightDemographicSensitivityAnalysis.csv",sep=""),row.names = FALSE)





#############  end of epiAnalysisWindBirthWeight_HEI_4970.R ############# 