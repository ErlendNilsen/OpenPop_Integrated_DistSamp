library(tidyverse)

if(!require("LivingNorwayR", character.only = TRUE)){
  devtools::install_github("LivingNorway/LivingNorwayR")
  library(LivingNorwayR,, character.only = TRUE)
}

# SETUP #
#-------#

## Define seed for initial value simulation and MCMC
mySeed <- 32

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Set switches 

# (Re-)downloading data
# downloadData <- FALSE
downloadData <- TRUE

# Aggregation to area level
areaAggregation <- FALSE

# Recruitment per adult or per adult female
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE

# Aggregation level for reproduction data
# NOTE: if this is not defined, will default to group level
sumR.Level <- "line" # Summing at the line level

# Time variation in survival
survVarT <- FALSE

# Rodent covariate on reproduction
fitRodentCov <- TRUE

## Switch for telemetry data (?)

telemetryData <- TRUE

# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData){
  #Rype_arkiv <- downloadLN(datasets = "Fjellstyrene", versions = 1.6, save = TRUE)
  Rype_arkiv <- downloadLN(datasets = c("Fjellstyrene"), versions = 1.7, save = TRUE)
}else{
  stop("downloadData = FALSE not supported yet. There is an issue with encoding when using LivingNorwayR::initializeDwCArchive() that needs to be resolved first.")
  #Rype_arkiv <- initializeDwCArchive("data/Rype_arkiv.zip")
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities/areas and time period of interest
localities <- listLocations()
areas <- listAreas()
minYear <- 2007
maxYear <- 2021

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
LT_data <- wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                                 duplTransects = duplTransects,
                                 localities = localities,
                                 #areas = areas,
                                 areaAggregation = areaAggregation,
                                 minYear = minYear, maxYear = maxYear)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

## Read in and reformat CMR data
d_cmr <- wrangleData_CMR(minYear = minYear)



# WRANGLE RODENT DATA #
#---------------------#

## Load and reformat rodent data
d_rodent <- wrangleData_Rodent(duplTransects = duplTransects,
                               localities = localities,
                               #areas = areas,
                               areaAggregation = areaAggregation,
                               minYear = minYear, maxYear = maxYear)


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = d_cmr,
                               d_rodent = d_rodent,
                               localities = localities, 
                               #areas = areas,
                               areaAggregation = areaAggregation,
                               excl_neverObs = TRUE,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               sumR.Level = "line",
                               dataVSconstants = TRUE,
                               save = TRUE)


# MODEL SETUP #
#-------------#

## Expand seeds for simulating initial values
nchains <- 3
MCMC.seeds <- expandSeed_MCMC(seed = mySeed, 
                              nchains = nchains)

## Write model code
modelCode <- writeModelCode(survVarT = survVarT)

## Setup for model using nimbleDistance::dHN
model_setup <- setupModel(modelCode = modelCode,
                          R_perF = R_perF,
                          survVarT = survVarT, 
                          fitRodentCov = fitRodentCov,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          testRun = TRUE, nchains = 3,
                          niter = 20000, nthin = 30, 
                          nburn = 11100,
                          initVals.seed = MCMC.seeds)

# MODEL (TEST) RUN #
#------------------#
t.start <- Sys.time()
IDSM.out <- nimbleMCMC(code = model_setup$modelCode,
                       data = input_data$nim.data, 
                       constants = input_data$nim.constants,
                       inits = model_setup$initVals, 
                       monitors = model_setup$modelParams,
                       nchains = model_setup$mcmcParams$nchains, 
                       niter = model_setup$mcmcParams$niter, 
                       nburnin = model_setup$mcmcParams$nburn, 
                       thin = model_setup$mcmcParams$nthin, 
                       samplesAsCodaMCMC = TRUE,
                       setSeed = mySeed)
Sys.time() - t.start

saveRDS(IDSM.out, file = 'rypeIDSM_dHN_multiArea_realData_Lierne.rds')


# TIDY UP POSTERIOR SAMPLES #
#---------------------------#

IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out, save = FALSE)
saveRDS(IDSM.out.tidy, file = 'rypeIDSM_dHN_multiArea_realData_Lierne_tidy.rds')


# OPTIONAL: MCMC TRACE PLOTS #
#----------------------------#

plotMCMCTraces(mcmc.out = IDSM.out.tidy,
               fitRodentCov = fitRodentCov,
               survVarT = survVarT)


# OPTIONAL: TIME SERIES PLOTS #
#-----------------------------#

plotTimeSeries(mcmc.out = IDSM.out.tidy, 
               N_areas = input_data$nim.constant$N_areas, 
               area_names = input_data$nim.constant$area_names, 
               N_sites = input_data$nim.constant$N_sites, 
               min_years = input_data$nim.constant$min_years, 
               max_years = input_data$nim.constant$max_years, 
               minYear = minYear, maxYear = maxYear,
               VitalRates = TRUE, DetectParams = TRUE, Densities = TRUE,
               showDataWindow = FALSE)


# OPTIONAL: PLOT VITAL RATE POSTERIORS #
#--------------------------------------#

plotPosteriorDens_VR(mcmc.out = IDSM.out.tidy,
                     N_areas = input_data$nim.constant$N_areas, 
                     area_names = input_data$nim.constant$area_names, 
                     N_years = input_data$nim.constant$N_years,
                     minYear = minYear,
                     survAreaIdx = input_data$nim.constants$SurvAreaIdx,
                     survVarT = survVarT,
                     fitRodentCov = fitRodentCov) 


# OPTIONAL: PLOT COVARIATE PREDICTIONS #
#--------------------------------------#

if(fitRodentCov){
  plotCovPrediction(mcmc.out = IDSM.out.tidy,
                    effectParam = "betaR.R",
                    covName = "Rodent occupancy",
                    minCov = 0, 
                    maxCov = 1,
                    meanCov = d_rodent$meanCov,
                    sdCov = d_rodent$sdCov,
                    N_areas = input_data$nim.constant$N_areas, 
                    area_names = input_data$nim.constant$area_names,
                    fitRodentCov = fitRodentCov)
}


# OPTIONAL: PLOT DETECTION FUNCTIONS #
#------------------------------------#

plotDetectFunction(mcmc.out = IDSM.out.tidy,
                   maxDist = input_data$nim.constants$W,
                   N_areas = input_data$nim.constant$N_areas, 
                   area_names = input_data$nim.constant$area_names)
