##################################################
#                                                #
#     Autologistic Community Abundance Model     #
#                                                #
#     Last updated: 26 Jul 2023                  #
#                                                #
##################################################
# This model is a community abundance model, made temporally dependent through the use of an
# autologistic term

#### MODEL SET-UP ####
# load functions used to clean data
functions_to_load <- list.files("./functions/", full.names = TRUE)
for(fn in functions_to_load){
  source(fn)
}
library(tidyr)
library(readr)

# Abundance data needs to be a 4-D array [species, sites, nights, seasons]
# read in the data
data_sets <- list.files("./data/abunTables/", full.names = TRUE)
df <- readr::read_csv(data_sets, id = "Species")
dets <- df[ ,3:20] # 3:20 for full dataset, 3:11 for 2-season test

# stack the data, long-style
stack <- wide_to_stacked(dets, 3)

# stack the replicates (nights) on top of each other and order the datasheet. 
# Data should be in order by species, then season, then site
n_long <- gather(stack, key = "night", value = "N", night1:night3, factor_key = TRUE)
n_long <- n_long[order(
  n_long[,"Species"], 
  n_long[,"Season"],
  n_long[,"Site"]), ]
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

# Site covariates are a 2-D array [sites, covariates]
site_covs <- read.csv("./data/siteCovs.csv", stringsAsFactors = FALSE)
siteCovs <- scale(site_covs[ ,2:5])
# a PCA collapses the co-linear variables fairly neatly
pca <- princomp(siteCovs)
# PC1: more + = canopy/herb/shrub; more - = humanMod
# PC2: more - = invasive shrub dominated sites; more + = other sites
covs_for_model <- cbind.data.frame("PC1" = pca$scores[ ,1],
                                   "PC2" = pca$scores[ ,2],
                                   "contag" = scale(site_covs$contag),
                                   "site" = site_covs[ ,1])
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
write_csv(covs_for_model, "./data/covs_for_model.csv")
covMatrix <- as.matrix(covs_for_model)
rm(siteCovs, covs_for_model)

# Observation covariates are 3-D array [sites, nights, seasons]
# We have 3 separate arrays, one for each covariate
moon <- read.csv("./data/obsVars/Moon.csv", stringsAsFactors = FALSE)
moon1 <- as.data.frame(moon[,2:19]) # to 19 for full dataset, to 10 for test
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

seasonData <- as.factor(c("spring","summer","fall","spring","summer","fall"))

#### RUN MODEL ####
# Data list for model
data_list <- list(
  nspec = max(n_long$Species),
  nsite = max(n_long$Site),
  PC1 = covMatrix[,1],
  PC2 = covMatrix[,2],
  contag = covMatrix[,3],
  nNight = max(n_long$night),
  moon = moonArray,
  jDate = jDateArray,
  effort = effortArray,
  y = my_array, 
  nseason = max(n_long$Season),
  season = seasonData
)

# Fit the model
library(runjags)
my_start <- Sys.time()
my_mod <- runjags::run.jags(
  model = "./JAGS/dynamicCommunityModel.R",
  monitor = c(# hyperprior parameters
              "mu.beta0", "sig.beta0", "mu.beta1", "sig.beta1", "mu.beta2", "sig.beta2", 
              "mu.beta3", "sig.beta3", "gamma",
              "mu.alpha0", "sig.alpha0", "alpha1", "alpha2", "alpha3",
              "mu.phi", "sig.phi",
              # species-specific coefficients
              "beta0", "beta1", "beta2", "beta3", "alpha0", "phi"
              ),
  data = data_list,
  n.chains = 3,
  inits = my_inits,
  burnin = 70000,
  sample = 20000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()
saveRDS(my_mod, "./results/autolog_mod.RDS")

#### SUMMARIZE MODEL ####
summary(my_mod)
plot(my_mod)

library(coda)
library(MCMCvis)
mc <- as.mcmc(my_mod)
model_sum <- MCMCvis::MCMCsummary(mc, 
                                  params = c("beta0", "beta1", "beta2", "beta3", 
                                             "alpha0", "alpha1", "alpha2", "alpha3",
                                             "gamma", "phi"),
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                                  round = 2)
model_sum

