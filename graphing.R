library(runjags)
library(coda)
library(dplyr)
library(ggplot2)
# load functions used to clean data
setwd("~/Documents/For School/Iowa/Dissertation Research/Small Mammals/SMamm Community")
functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}

# read in the results
my_mod <- readRDS("./results/autolog_mod.RDS")
site_covs <- read.csv("./data/siteCovs.csv", stringsAsFactors = FALSE)
covs_for_model <- read.csv("./data/covs_for_model.csv", stringsAsFactors = FALSE)

# convert the model file to an mcmc object; use MCMCvis to get summary statistics
mc <- as.mcmc(my_mod)
library(MCMCvis)
model_sum <- MCMCvis::MCMCsummary(mc, 
                                  params = c("alpha1","alpha2","alpha3"),
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                                  round = 2)
model_sum

# convert the mcmc object to a matrix
mc <- as.matrix(mc)
# sub-sample the mcmc matrix a bit as we don't really need to make predictions with 
# all 60K samples
set.seed(554)
mc_sub <- mc[sample(1:nrow(mc), 10000), ]
# and use split_mcmc
mc <- split_mcmc(mc_sub)
rm(mc_sub)
# check out the dimensions of one of the list elements. for this example (a slope term for my 
# occupancy covariate), it should have a number of rows equal to the number of MCMC samples 
# (10000) & a number of columns equal to the number of species in the model.
dim(mc$beta1)
# [1] 10000  6
# the among-species terms (e.g., mu_beta0) will have a 1 as the 2nd term

##### PC1 ( - humanMod -- + Other) #####
# generate a sequence of contag values
range(covs_for_model$PC1)
# [1] -3.768367  1.913313
# I need to choose some 'pretty' numbers based on that
predV <- seq(-3.7, 2.0, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  predV,
  mean(covs_for_model$PC2),
  mean(scale(site_covs$contag))
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  predV,
  mean(covs_for_model$PC2),
  mean(scale(site_covs$contag)),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_lambda <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_lambda[,,1] <- exp(cbind(mc$beta0[,1],mc$beta1[,1],mc$beta2[,1],mc$beta3[,1]) %*% pred_mat)
pred_lambda[,,2] <- exp(cbind(mc$beta0[,2],mc$beta1[,2],mc$beta2[,2],mc$beta3[,2]) %*% pred_mat)
pred_lambda[,,3] <- exp(cbind(mc$beta0[,3],mc$beta1[,3],mc$beta2[,3],mc$beta3[,3]) %*% pred_mat)
pred_lambda[,,4] <- exp(cbind(mc$beta0[,4],mc$beta1[,4],mc$beta2[,4],mc$beta3[,4]) %*% pred_mat)
pred_lambda[,,5] <- exp(cbind(mc$beta0[,5],mc$beta1[,5],mc$beta2[,5],mc$beta3[,5]) %*% pred_mat)
pred_lambda[,,6] <- exp(cbind(mc$beta0[,6],mc$beta1[,6],mc$beta2[,6],mc$beta3[,6]) %*% pred_mat)

# create the matrix of the coefficients from the model
pred_lambda_auto <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat_auto),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_lambda_auto[,,1] <- exp(cbind(mc$beta0[,1],mc$beta1[,1],mc$beta2[,1],mc$beta3[,1],mc$phi[,1]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,2] <- exp(cbind(mc$beta0[,2],mc$beta1[,2],mc$beta2[,2],mc$beta3[,2],mc$phi[,2]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,3] <- exp(cbind(mc$beta0[,3],mc$beta1[,3],mc$beta2[,3],mc$beta3[,3],mc$phi[,3]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,4] <- exp(cbind(mc$beta0[,4],mc$beta1[,4],mc$beta2[,4],mc$beta3[,4],mc$phi[,4]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,5] <- exp(cbind(mc$beta0[,5],mc$beta1[,5],mc$beta2[,5],mc$beta3[,5],mc$phi[,5]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,6] <- exp(cbind(mc$beta0[,6],mc$beta1[,6],mc$beta2[,6],mc$beta3[,6],mc$phi[,6]) 
                             %*% pred_mat_auto)

# calculating average lambda
trueCount <- pred_lambda / (
  pred_lambda + (1 - pred_lambda_auto)
)
trueCount <- apply(
  trueCount,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
PC1 <- rbind(data.frame(group = "Northern short-tailed shrew",
                           pc1 = seq(-3.7, 2.0, length.out = 200),
                           lambda = trueCount[3,,1],
                           upper = trueCount[5,,1],
                           lower = trueCount[1,,1]),
                data.frame(group = "Prairie vole",
                           pc1 = seq(-3.7, 2.0, length.out = 200),
                           lambda = trueCount[3,,2],
                           upper = trueCount[5,,2],
                           lower = trueCount[1,,2]),
                data.frame(group = "Meadow vole",
                           pc1 = seq(-3.7, 2.0, length.out = 200),
                           lambda = trueCount[3,,3],
                           upper = trueCount[5,,3],
                           lower = trueCount[1,,3]),
                data.frame(group = "Deer mouse spp.",
                           pc1 = seq(-3.7, 2.0, length.out = 200),
                           lambda = trueCount[3,,4],
                           upper = trueCount[5,,4],
                           lower = trueCount[1,,4]),
                data.frame(group = "Harvest mouse",
                           pc1 = seq(-3.7, 2.0, length.out = 200),
                           lambda = trueCount[3,,5],
                           upper = trueCount[5,,5],
                           lower = trueCount[1,,5]),
                data.frame(group = "Meadow jumping mouse",
                           pc1 = seq(-3.7, 2.0, length.out = 200),
                           lambda = trueCount[3,,6],
                           upper = trueCount[5,,6],
                           lower = trueCount[1,,6]))
jpeg("./results/PC1.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(PC1, aes(x=pc1, y=lambda)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.2, show.legend=FALSE) +
  scale_fill_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  geom_line(aes(color=group), linewidth = 1.05) +
  scale_color_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  ylim(0,6) +
  labs(x="PC1", y="Local Abundance", color="Species") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title=element_blank(), legend.position=c(0.15,0.89), legend.background=element_blank(),
        legend.text=element_text(size=6)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
dev.off()

##### PC2 ( - shrubby -- + Other) #####
# generate a sequence of contag values
range(covs_for_model$PC2)
# [1] -3.184049  1.612469
# I need to choose some 'pretty' numbers based on that
predV <- seq(-3.2, 1.7, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  mean(covs_for_model$PC1),
  predV,
  mean(scale(site_covs$contag))
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  mean(covs_for_model$PC1),
  predV,
  mean(scale(site_covs$contag)),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_lambda <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_lambda[,,1] <- exp(cbind(mc$beta0[,1],mc$beta1[,1],mc$beta2[,1],mc$beta3[,1]) %*% pred_mat)
pred_lambda[,,2] <- exp(cbind(mc$beta0[,2],mc$beta1[,2],mc$beta2[,2],mc$beta3[,2]) %*% pred_mat)
pred_lambda[,,3] <- exp(cbind(mc$beta0[,3],mc$beta1[,3],mc$beta2[,3],mc$beta3[,3]) %*% pred_mat)
pred_lambda[,,4] <- exp(cbind(mc$beta0[,4],mc$beta1[,4],mc$beta2[,4],mc$beta3[,4]) %*% pred_mat)
pred_lambda[,,5] <- exp(cbind(mc$beta0[,5],mc$beta1[,5],mc$beta2[,5],mc$beta3[,5]) %*% pred_mat)
pred_lambda[,,6] <- exp(cbind(mc$beta0[,6],mc$beta1[,6],mc$beta2[,6],mc$beta3[,6]) %*% pred_mat)

# create the matrix of the coefficients from the model
pred_lambda_auto <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat_auto),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_lambda_auto[,,1] <- exp(cbind(mc$beta0[,1],mc$beta1[,1],mc$beta2[,1],mc$beta3[,1],mc$phi[,1]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,2] <- exp(cbind(mc$beta0[,2],mc$beta1[,2],mc$beta2[,2],mc$beta3[,2],mc$phi[,2]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,3] <- exp(cbind(mc$beta0[,3],mc$beta1[,3],mc$beta2[,3],mc$beta3[,3],mc$phi[,3]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,4] <- exp(cbind(mc$beta0[,4],mc$beta1[,4],mc$beta2[,4],mc$beta3[,4],mc$phi[,4]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,5] <- exp(cbind(mc$beta0[,5],mc$beta1[,5],mc$beta2[,5],mc$beta3[,5],mc$phi[,5]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,6] <- exp(cbind(mc$beta0[,6],mc$beta1[,6],mc$beta2[,6],mc$beta3[,6],mc$phi[,6]) 
                             %*% pred_mat_auto)

# calculating average lambda
trueCount <- pred_lambda / (
  pred_lambda + (1 - pred_lambda_auto)
)
trueCount <- apply(
  trueCount,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
PC2 <- rbind(data.frame(group = "Northern short-tailed shrew",
                        pc2 = seq(-3.2, 1.7, length.out = 200),
                        lambda = trueCount[3,,1],
                        upper = trueCount[5,,1],
                        lower = trueCount[1,,1]),
             data.frame(group = "Prairie vole",
                        pc2 = seq(-3.2, 1.7, length.out = 200),
                        lambda = trueCount[3,,2],
                        upper = trueCount[5,,2],
                        lower = trueCount[1,,2]),
             data.frame(group = "Meadow vole",
                        pc2 = seq(-3.2, 1.7, length.out = 200),
                        lambda = trueCount[3,,3],
                        upper = trueCount[5,,3],
                        lower = trueCount[1,,3]),
             data.frame(group = "Deer mouse spp.",
                        pc2 = seq(-3.2, 1.7, length.out = 200),
                        lambda = trueCount[3,,4],
                        upper = trueCount[5,,4],
                        lower = trueCount[1,,4]),
             data.frame(group = "Harvest mouse",
                        pc2 = seq(-3.2, 1.7, length.out = 200),
                        lambda = trueCount[3,,5],
                        upper = trueCount[5,,5],
                        lower = trueCount[1,,5]),
             data.frame(group = "Meadow jumping mouse",
                        pc2 = seq(-3.2, 1.7, length.out = 200),
                        lambda = trueCount[3,,6],
                        upper = trueCount[5,,6],
                        lower = trueCount[1,,6]))
jpeg("./results/PC2.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(PC2, aes(x=pc2, y=lambda)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.2, show.legend=FALSE) +
  scale_fill_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  geom_line(aes(color=group), linewidth = 1.05) +
  scale_color_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  ylim(0,6) +
  labs(x="PC2", y="Local Abundance", color="Species") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title=element_blank(), legend.position=c(0.15,0.89), legend.background=element_blank(),
        legend.text=element_text(size=6)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
dev.off()
##### CONTAG #####
# generate a sequence of contag values
range(site_covs$contag)
# [1] 13.70055 100.00000
# I need to choose some 'pretty' numbers based on that
predV <- seq(13, 100, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  mean(covs_for_model$PC1),
  mean(covs_for_model$PC2),
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  mean(covs_for_model$PC1),
  mean(covs_for_model$PC2),
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_lambda <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_lambda[,,1] <- exp(cbind(mc$beta0[,1],mc$beta1[,1],mc$beta2[,1],mc$beta3[,1]) %*% pred_mat)
pred_lambda[,,2] <- exp(cbind(mc$beta0[,2],mc$beta1[,2],mc$beta2[,2],mc$beta3[,2]) %*% pred_mat)
pred_lambda[,,3] <- exp(cbind(mc$beta0[,3],mc$beta1[,3],mc$beta2[,3],mc$beta3[,3]) %*% pred_mat)
pred_lambda[,,4] <- exp(cbind(mc$beta0[,4],mc$beta1[,4],mc$beta2[,4],mc$beta3[,4]) %*% pred_mat)
pred_lambda[,,5] <- exp(cbind(mc$beta0[,5],mc$beta1[,5],mc$beta2[,5],mc$beta3[,5]) %*% pred_mat)
pred_lambda[,,6] <- exp(cbind(mc$beta0[,6],mc$beta1[,6],mc$beta2[,6],mc$beta3[,6]) %*% pred_mat)

# create the matrix of the coefficients from the model
pred_lambda_auto <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat_auto),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_lambda_auto[,,1] <- exp(cbind(mc$beta0[,1],mc$beta1[,1],mc$beta2[,1],mc$beta3[,1],mc$phi[,1]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,2] <- exp(cbind(mc$beta0[,2],mc$beta1[,2],mc$beta2[,2],mc$beta3[,2],mc$phi[,2]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,3] <- exp(cbind(mc$beta0[,3],mc$beta1[,3],mc$beta2[,3],mc$beta3[,3],mc$phi[,3]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,4] <- exp(cbind(mc$beta0[,4],mc$beta1[,4],mc$beta2[,4],mc$beta3[,4],mc$phi[,4]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,5] <- exp(cbind(mc$beta0[,5],mc$beta1[,5],mc$beta2[,5],mc$beta3[,5],mc$phi[,5]) 
                             %*% pred_mat_auto)
pred_lambda_auto[,,6] <- exp(cbind(mc$beta0[,6],mc$beta1[,6],mc$beta2[,6],mc$beta3[,6],mc$phi[,6]) 
                             %*% pred_mat_auto)

# calculating average lambda
trueCount <- pred_lambda / (
  pred_lambda + (1 - pred_lambda_auto)
)
trueCount <- apply(
  trueCount,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
Contag <- rbind(data.frame(group = "Northern short-tailed shrew",
                           contag = seq(13, 100, length.out = 200),
                           lambda = trueCount[3,,1],
                           upper = trueCount[5,,1],
                           lower = trueCount[1,,1]),
                data.frame(group = "Prairie vole",
                           contag = seq(13, 100, length.out = 200),
                           lambda = trueCount[3,,2],
                           upper = trueCount[5,,2],
                           lower = trueCount[1,,2]),
                data.frame(group = "Meadow vole",
                           contag = seq(13, 100, length.out = 200),
                           lambda = trueCount[3,,3],
                           upper = trueCount[5,,3],
                           lower = trueCount[1,,3]),
                data.frame(group = "Deer mouse spp.",
                           contag = seq(13, 100, length.out = 200),
                           lambda = trueCount[3,,4],
                           upper = trueCount[5,,4],
                           lower = trueCount[1,,4]),
                data.frame(group = "Harvest mouse",
                           contag = seq(13, 100, length.out = 200),
                           lambda = trueCount[3,,5],
                           upper = trueCount[5,,5],
                           lower = trueCount[1,,5]),
                data.frame(group = "Meadow jumping mouse",
                           contag = seq(13, 100, length.out = 200),
                           lambda = trueCount[3,,6],
                           upper = trueCount[5,,6],
                           lower = trueCount[1,,6]))
jpeg("./results/contag.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(Contag, aes(x=contag, y=lambda)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.2, show.legend=FALSE) +
  scale_fill_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  geom_line(aes(color=group), linewidth = 1.05) +
  scale_color_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  ylim(0,6) +
  labs(x="Contagion Index", y="Local Abundance", color="Species") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title=element_blank(), legend.position=c(0.15,0.89), legend.background=element_blank(),
        legend.text=element_text(size=6)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
dev.off()
##### Moon Illumination #####
# generate a sequence of moon illumination values
range(moon_long$moon)
# [1] -1.343997  1.474745
predV <- seq(-1.343997, 1.474745, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  predV,
  0,  # mean of scale(temp)
  0   # mean of scale(effort)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# create the matrix of the coefficients from the model
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_p[,,1] <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,2] <- cbind(mc$alpha0[,2],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,3] <- cbind(mc$alpha0[,3],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,4] <- cbind(mc$alpha0[,4],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,5] <- cbind(mc$alpha0[,5],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,6] <- cbind(mc$alpha0[,6],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

# convert to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
prob_p <- apply(pred_p,
                c(2,3),
                logit2prob)
prob_p <- apply(
  prob_p,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
moon <- rbind(data.frame(group = "Northern short-tailed shrew",
                         moon = seq(0, 1, length.out = 200),
                         p = prob_p[3,,1],
                         upper = prob_p[5,,1],
                         lower = prob_p[1,,1]),
                data.frame(group = "Prairie vole",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,2],
                           upper = prob_p[5,,2],
                           lower = prob_p[1,,2]),
                data.frame(group = "Meadow vole",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,3],
                           upper = prob_p[5,,3],
                           lower = prob_p[1,,3]),
                data.frame(group = "Deer mouse spp.",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,4],
                           upper = prob_p[5,,4],
                           lower = prob_p[1,,4]),
                data.frame(group = "Harvest mouse",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,5],
                           upper = prob_p[5,,5],
                           lower = prob_p[1,,5]),
                data.frame(group = "Meadow jumping mouse",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,6],
                           upper = prob_p[5,,6],
                           lower = prob_p[1,,6]))
jpeg("./results/moon.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(moon, aes(x=moon, y=p)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.2, show.legend=FALSE) +
  scale_fill_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  geom_line(aes(color=group), linewidth = 1.05) +
  scale_color_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  ylim(0,1) +
  labs(x="Moon Illumination (proportion full)", y="Capture Probability", color="Species") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title=element_blank(), legend.position=c(0.15,0.89), legend.background=element_blank(),
        legend.text=element_text(size=6)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
dev.off()

##### Effort #####
# generate a sequence of moon illumination values
effort <- read.csv("./data/obsVars/Effort.csv", stringsAsFactors = FALSE)
effortScale <- as.data.frame(scale(effort[,2:19]))
range(effortScale)
# [1] -2.544777  1.402502
predV <- c(-2.5447767, -2.1500488, -1.7553209, -1.3605930, -0.9658651, 
           -0.5711372, -0.1764093, 0.2183186, 0.6130466, 1.0077745, 1.4025024)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  0,     # mean of scale(moon)
  0,     # mean of scale(temp)
  predV
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# create the matrix of the coefficients from the model
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat),
    6
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_p[,,1] <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,2] <- cbind(mc$alpha0[,2],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,3] <- cbind(mc$alpha0[,3],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,4] <- cbind(mc$alpha0[,4],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,5] <- cbind(mc$alpha0[,5],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,6] <- cbind(mc$alpha0[,6],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

# convert to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
prob_p <- apply(pred_p,
                c(2,3),
                logit2prob)
prob_p <- apply(
  prob_p,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
effort <- rbind(data.frame(group = "Northern short-tailed shrew",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,1],
                         upper = prob_p[5,,1],
                         lower = prob_p[1,,1]),
              data.frame(group = "Prairie vole",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,2],
                         upper = prob_p[5,,2],
                         lower = prob_p[1,,2]),
              data.frame(group = "Meadow vole",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,3],
                         upper = prob_p[5,,3],
                         lower = prob_p[1,,3]),
              data.frame(group = "Deer mouse spp.",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,4],
                         upper = prob_p[5,,4],
                         lower = prob_p[1,,4]),
              data.frame(group = "Harvest mouse",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,5],
                         upper = prob_p[5,,5],
                         lower = prob_p[1,,5]),
              data.frame(group = "Meadow jumping mouse",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,6],
                         upper = prob_p[5,,6],
                         lower = prob_p[1,,6]))
jpeg("./results/effort.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(effort, aes(x=effort, y=p, color=group)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), show.legend=FALSE) +
  scale_fill_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  geom_point() +
  scale_color_manual(values = c("#cc0c00","#7C878E","#84BD00","#786197","#5C88DA","#FFCD00")) +
  ylim(0,1) +
  labs(x="Trap Effort", y="Capture Probability", color="Species") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title=element_blank(), legend.position=c(0.15,0.89), legend.background=element_blank(),
        legend.text=element_text(size=6)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
dev.off()
