##################################################
#                                                #
#     Dynamic Community Abundance Model          #
#                                                #
#     Last updated: 29 APR 2024                  #
#                                                #
##################################################
# This model is a community abundance model with fully modeled 'recruitment' (births + 
# immigration) and 'survival' (1-[deaths + emigration]).

#### MODEL SET-UP ####
# load functions used to clean data & other processes
functions_to_load <- list.files("./functions/", full.names = TRUE)
for(fn in functions_to_load){
  source(fn)
}
# load some libraries to tidy-up the data
library(wrapr)
library(tidyr)
library(readr)

#### Read & Process Count Data ####
# Abundance data needs to be a 4-D array [species, sites, nights, seasons]
# read in the data
data_sets <- list.files("./data/abunTables/", full.names = TRUE)
df <- readr::read_csv(data_sets, id = "Species")
dets <- df[ ,3:20]

# stack the data, long-style
# note: wide_to_stacked is a custom function that does require you to have
# an object 'df' loaded with columns named 'Site' & 'Species'
stack <- wide_to_stacked(dets, 3)

# stack the replicates (nights) on top of each other and order the datasheet. 
# Data should be in order by species, then season, then site
n_long <- gather(stack, key = "night", value = "N", night1:night3, factor_key = TRUE)
n_long <- n_long[orderv(c(
  n_long[,"Species"], 
  n_long[,"Season"],
  n_long[,"Site"])), ]
n_long$night <- as.numeric(n_long$night)
n_long$Site <- as.numeric(as.factor(n_long$Site))
n_long$Species <- as.numeric(as.factor(n_long$Species))

# create a blank array that has the correct dimensions
my_array <- array(
  NA,
  dim = c(
    max(n_long$Species),
    max(n_long$Site),
    max(n_long$night),
    max(n_long$Season)
  )
)
# fill in the array with the information from the long-form data.frame
for(i in 1:nrow(n_long)){
  my_array[
    n_long$Species[i],
    n_long$Site[i],
    n_long$night[i],
    n_long$Season[i]
  ] <- n_long$N[i]
}
rm(dets, stack)

#### Read & Process Site Covariates ####
# Site covariates are a 2-D array [sites, covariates]
site_covs <- read.csv("./data/siteCovs.csv", stringsAsFactors = FALSE)
siteCovs <- scale(site_covs[ ,2:5])
# a PCA collapses the co-linear variables fairly neatly
pca <- princomp(siteCovs)
pca$loadings
# PC1: more + = canopy/herb/shrub; more - = humanMod
# let's make our covariate data.frame
covs_for_model <- cbind.data.frame("PC1" = pca$scores[ ,1],
                                   "contag" = scale(site_covs$contag),
                                   "site" = site_covs[ ,1])
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
# saving the results of the PCA for later use
write_csv(covs_for_model, "./data/covs_for_model.csv")
# turn your covariate data.frame into a matrix
covMatrix <- as.matrix(covs_for_model)
rm(siteCovs, covs_for_model)

#### Read & Process Observation Covariates ####
# Observation covariates are 3-D array [sites, nights, seasons]
# We have 3 separate arrays, one for each covariate
moon <- read.csv("./data/obsVars/Moon.csv", stringsAsFactors = FALSE)
moon1 <- as.data.frame(moon[,2:19])
# 'wideObs_to_stacked' is a custom function for stacking the data long-style
# it does require 'moon' to be loaded with a column named 'Site'
moonStack <- wideObs_to_stacked(moon1, 3) 
moon_long <- gather(moonStack, key = "night", value = "moon", night1:night3, factor_key = TRUE)
moon_long$moon <- as.numeric(scale(moon_long$moon))
moon_long <- moon_long[order(
  moon_long[,"Season"],
  moon_long[,"Site"]),]
moon_long$night <- as.numeric(moon_long$night)
moon_long$Site <- as.numeric(as.factor(moon_long$Site))
moonArray <- array(
  NA,
  dim = c(max(moon_long$Site),
          max(moon_long$night),
          max(moon_long$Season)
          )
  )
for(i in 1:nrow(moon_long)){
  moonArray[
    moon_long$Site[i],
    moon_long$night[i],
    moon_long$Season[i]
  ] <- moon_long$moon[i]
}
rm(moon1, moonStack, moon_long)

jDate <- read.csv("./data/obsVars/jDate.csv", stringsAsFactors = FALSE)
jDate1 <- as.data.frame(jDate[,2:19]) # to 19 for full dataset, to 10 for test
jDateStack <- wideObs_to_stacked(jDate1, 3) 
jDate_long <- gather(jDateStack, key = "night", value = "jDate", night1:night3, factor_key = TRUE)
jDate_long$jDate <- as.numeric(scale(jDate_long$jDate))
jDate_long <- jDate_long[order(
  jDate_long[,"Season"],
  jDate_long[,"Site"]),]
jDate_long$night <- as.numeric(jDate_long$night)
jDate_long$Site <- as.numeric(as.factor(jDate_long$Site))
jDateArray <- array(
  NA,
  dim = c(max(jDate_long$Site),
          max(jDate_long$night),
          max(jDate_long$Season)
  )
)
for(i in 1:nrow(jDate_long)){
  jDateArray[
    jDate_long$Site[i],
    jDate_long$night[i],
    jDate_long$Season[i]
  ] <- jDate_long$jDate[i]
}
rm(jDate1, jDateStack, jDate_long, jDate)

effort <- read.csv("./data/obsVars/Effort.csv", stringsAsFactors = FALSE)
effort1 <- as.data.frame(effort[,2:19]) # to 19 for full dataset, to 10 for test
effortStack <- wideObs_to_stacked(effort1, 3) 
effort_long <- gather(effortStack, key = "night", value = "effort", night1:night3, factor_key = TRUE)
effort_long$effort <- as.numeric(scale(effort_long$effort))
effort_long <- effort_long[order(
  effort_long[,"Season"],
  effort_long[,"Site"]),]
effort_long$night <- as.numeric(effort_long$night)
effort_long$Site <- as.numeric(as.factor(effort_long$Site))
effortArray <- array(
  NA,
  dim = c(max(effort_long$Site),
          max(effort_long$night),
          max(effort_long$Season)
  )
)
for(i in 1:nrow(effort_long)){
  effortArray[
    effort_long$Site[i],
    effort_long$night[i],
    effort_long$Season[i]
  ] <- effort_long$effort[i]
}
rm(effort1, effortStack, effort_long, effort, moon)

#### Generate "Season" Data ####
# will assign a correct description for each season of data
seasonData <- as.factor(c("spring","summer","fall","spring","summer","fall"))

#### RUN MODEL ####
# Data list for model
data_list <- list(
  nspec = max(n_long$Species),
  nsite = max(n_long$Site),
  PC1 = covMatrix[,1],
  contag = covMatrix[,2],
  nNight = max(n_long$night),
  moon = moonArray,
  jDate = jDateArray,
  effort = effortArray,
  y = my_array, 
  nSP = max(n_long$Season),
  season = seasonData,
  nseason = 3
)
# specify initial values for N (abundance) and R (recruits)
# initial values for N should be maximum possible counts, in our case 10
N_init = array(10, dim = c(data_list$nspec, data_list$nsite, data_list$nseason))
# set to NA for all seasons t>1 because we cannot know this (needs to be specified by R + S)
N_init[,,2:6] <- NA
# initial values for recruits can also be set to maximum possible counts
R_init = array(10, dim = c(data_list$nspec, data_list$nsite, data_list$nseason))
# providing null values for t = 1 (i.e., no recruits yet)
R_init[,,1] <- NA

my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      mu.beta0 = rnorm(1),
      tau.beta0 = rgamma(1,1,1),
      mu.beta1 = rnorm(1),
      tau.beta1 = rgamma(1,1,1),
      mu.alpha0 = rnorm(1),
      tau.alpha0 = rgamma(1,1,1),
      alpha1 = rnorm(1),
      alpha2 = rnorm(1),
      alpha3 = rnorm(1),
      mu.phi0 = rnorm(1),
      tau.phi0 = rgamma(1,1,1),
      mu.phi1 = rnorm(1),
      tau.phi1 = rgamma(1,1,1),
      mu.phi2 = rnorm(1),
      tau.phi2 = rgamma(1,1,1),
      mu.gamma0 = rnorm(1),
      tau.gamma0 = rgamma(1,1,1),
      mu.gamma1 = rnorm(1),
      tau.gamma1 = rgamma(1,1,1),
      tmp_gamma = rnorm(2),
      beta0 = rnorm(data_list$nspec),
      beta1 = rnorm(data_list$nspec),
      alpha0 = rnorm(data_list$nspec),
      phi0 = rnorm(data_list$nspec),
      phi1 = rnorm(data_list$nspec),
      phi2 = rnorm(data_list$nspec),
      gamma0 = rnorm(data_list$nspec),
      gamma1 = rnorm(data_list$nspec),
      R = R_init,
      # no init for S needed because N = R + S so S is fixed
      N = N_init,
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(
    switch(chain,
           "1" = gen_list(chain),
           "2" = gen_list(chain),
           "3" = gen_list(chain),
           "4" = gen_list(chain),
           "5" = gen_list(chain),
           "6" = gen_list(chain),
           "7" = gen_list(chain),
           "8" = gen_list(chain)
    )
  )
}

# Fit the model
library(runjags)
my_start <- Sys.time()
my_mod <- runjags::run.jags(
  model = "./JAGS/dynamicCommunityModel.R",
  monitor = c(# hyperprior (& non-species-specific) parameters
              "mu.beta0", "sig.beta0", "mu.beta1", "sig.beta1", 
              "mu.alpha0", "sig.alpha0", "alpha1", "alpha2", "alpha3", 
              "mu.phi0", "sig.phi0", "mu.phi1", "sig.phi1", "mu.phi2", "sig.phi2",
              "mu.gamma0", "sig.gamma0", "mu.gamma1", "sig.gamma1", "gamma2",
              # species-specific coefficients
              "beta0", "beta1", "alpha0", "phi0", "phi1", "phi2", "gamma0", "gamma1",
              # derived quantities
              "N"
              ),
  data = data_list,
  n.chains = 3,
  inits = my_inits,
  burnin = 150000,
  sample = 20000,
  adapt = 3500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()
#create a results folder and save model output there
dir.create("results")
saveRDS(my_mod, "./results/my_mod.RDS")

#### SUMMARIZE MODEL ####
summary(my_mod)
plot(my_mod)

library(coda)
library(dplyr)
library(MCMCvis)
mc <- as.mcmc(my_mod)
model_sum <- MCMCvis::MCMCsummary(my_mod$mcmc, 
                                  params = c("mu.alpha0", "sig.alpha0",
                                             "alpha0", "alpha1", "alpha2", "alpha3",
                                             "mu.beta0", "sig.beta0", "mu.beta1", "sig.beta1",
                                             "beta0", "beta1",
                                             "mu.phi0", "sig.phi0", "mu.phi1", "sig.phi1",
                                             "phi0", "phi1", "phi2",
                                             "mu.gamma0", "sig.gamma0", "mu.gamma1", "sig.gamma1",
                                             "gamma0", "gamma1", "gamma2"),
                                  probs = c(0.025, 0.5, 0.975), 
                                  round = 2)
model_sum

