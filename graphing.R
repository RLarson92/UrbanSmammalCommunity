####################################################################################
#                                                                                  #
#                     Code used to create Figures 2 - 4                            #
#                                                                                  #
#                        Last Updated: 25 Jan 2024                                 #
#                                                                                  #
####################################################################################

library(runjags)
library(coda)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# load functions used to clean data
functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}

# read in the results & covariates for graphing
my_mod <- readRDS("./results/my_mod.RDS")
covs_for_model <- read.csv("./data/covs_for_model.csv", stringsAsFactors = FALSE)

# convert the model file to an mcmc object
# note: both 'runjags' and 'coda' need to be loaded for this to work
mc <- as.mcmc(my_mod)
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
# abundance covariate), it should have a number of rows equal to the number of MCMC samples 
# (10000) & a number of columns equal to the number of species in the model.
dim(mc$beta1)
# [1] 10000  6
# the among-species terms (e.g., mu_beta0) will have a 1 as the 2nd term

#####                            Persistence (phi[i,j,t])                                   #####
##### PC1 ( - humanMod -- + Other) #####
# the range of PC1 is -3.768367  1.913313
# I need to choose some 'pretty' numbers based on that
predV <- seq(-3.7, 2.0, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# PC1 values, and adding a row of 1's for the intercept & the mean (0) for other predictor values
pred_mat <- cbind(
  1,      # intercept
  predV,  # PC1 values
  0       # contag, held at its mean
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence from above, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                               covariates except your variable of interest)

# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat),
    7               # Number of species + 1 (for community average value)
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_phi[,,1] <- exp(cbind(mc$mu.phi0[,1],mc$mu.phi1[,1],mc$mu.phi2[,1]) %*% pred_mat)
pred_phi[,,2] <- exp(cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1]) %*% pred_mat)
pred_phi[,,3] <- exp(cbind(mc$phi0[,2],mc$phi1[,2],mc$phi2[,2]) %*% pred_mat)
pred_phi[,,4] <- exp(cbind(mc$phi0[,3],mc$phi1[,3],mc$phi2[,3]) %*% pred_mat)
pred_phi[,,5] <- exp(cbind(mc$phi0[,4],mc$phi1[,4],mc$phi2[,4]) %*% pred_mat)
pred_phi[,,6] <- exp(cbind(mc$phi0[,5],mc$phi1[,5],mc$phi2[,5]) %*% pred_mat)
pred_phi[,,7] <- exp(cbind(mc$phi0[,6],mc$phi1[,6],mc$phi2[,6]) %*% pred_mat)

# convert to probability
prob_phi <- apply(pred_phi,
                c(2,3),
                logit2prob)

truePhi <- apply(
  prob_phi,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

## Graphing
PC1 <- rbind(data.frame(group = "Community Average",
                        species = "Community Average",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,1],
                        upper = truePhi[5,,1],
                        lower = truePhi[1,,1]),
             data.frame(group = "Species-Specific",
                        species = "Short-tailed shrew",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,2],
                        upper = 0.5, #truePhi[4,,2],
                        lower = 0.5), #truePhi[2,,2]),
             data.frame(group = "Species-Specific",
                        species = "Prairie vole",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,3],
                        upper = 0.5, #truePhi[4,,3],
                        lower = 0.5), #truePhi[2,,3]),
             data.frame(group = "Species-Specific",
                        species = "Meadow vole",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,4],
                        upper = 0.5, #truePhi[4,,4],
                        lower = 0.5), #truePhi[2,,4]),
             data.frame(group = "Species-Specific",
                        species = "Deer mouse spp.",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,5],
                        upper = 0.5, #truePhi[4,,5],
                        lower = 0.5), #truePhi[2,,5]),
             data.frame(group = "Species-Specific",
                        species = "Harvest mouse",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,6],
                        upper = 0.5, #truePhi[4,,6],
                        lower = 0.5), #truePhi[2,,6]),
             data.frame(group = "Species-Specific",
                        species = "Jumping mouse",
                        pc1 = seq(-3.7, 2.0, length.out = 200),
                        phi = truePhi[3,,7],
                        upper = 0.5, #truePhi[4,,7],
                        lower = 0.5)) #truePhi[2,,7]))
p1<-ggplot(PC1, aes(x=pc1, y=phi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=species), alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=species), show.legend = FALSE) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5,0.5,0.5,0.5,0.5,0.5), guide="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0.5,1), expand = c(0,0)) +
  labs(x="PC1", y="Persistence Probability", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9))

##### CONTAG #####
# instructions are more or less the same as for PC1, just subbing in prediction values for
# contagion and setting PC1 to its mean value
predV <- seq(-1.33, 2.164393, length.out = 200)
pred_mat <- cbind(
  1,
  0,
  predV
)
pred_mat <- t(pred_mat)
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$beta0),
    ncol(pred_mat),
    7
  )
)
pred_phi[,,1] <- exp(cbind(mc$mu.phi0[,1],mc$mu.phi1[,1],mc$mu.phi2[,1]) %*% pred_mat)
pred_phi[,,2] <- exp(cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1]) %*% pred_mat)
pred_phi[,,3] <- exp(cbind(mc$phi0[,2],mc$phi1[,2],mc$phi2[,2]) %*% pred_mat)
pred_phi[,,4] <- exp(cbind(mc$phi0[,3],mc$phi1[,3],mc$phi2[,3]) %*% pred_mat)
pred_phi[,,5] <- exp(cbind(mc$phi0[,4],mc$phi1[,4],mc$phi2[,4]) %*% pred_mat)
pred_phi[,,6] <- exp(cbind(mc$phi0[,5],mc$phi1[,5],mc$phi2[,5]) %*% pred_mat)
pred_phi[,,7] <- exp(cbind(mc$phi0[,6],mc$phi1[,6],mc$phi2[,6]) %*% pred_mat)
prob_phi <- apply(pred_phi,
                  c(2,3),
                  logit2prob)
truePhi <- apply(
  prob_phi,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

## Graphing
# I converted the centered & scaled values for contagion back onto a 0-100 scale
# for easier interpretation
CONTAG <- rbind(data.frame(group = "Community Average",
                        species = "Community Average",
                        contag = seq(13, 100, length.out = 200), 
                        phi = truePhi[3,,1],
                        upper = truePhi[5,,1],
                        lower = truePhi[1,,1]),
             data.frame(group = "Species-Specific",
                        species = "Short-tailed shrew",
                        contag = seq(13, 100, length.out = 200),
                        phi = truePhi[3,,2],
                        upper = 0.5, #truePhi[4,,2],
                        lower = 0.5), #truePhi[2,,2]),
             data.frame(group = "Species-Specific",
                        species = "Prairie vole",
                        contag = seq(13, 100, length.out = 200),
                        phi = truePhi[3,,3],
                        upper = 0.5, #truePhi[4,,3],
                        lower = 0.5), #truePhi[2,,3]),
             data.frame(group = "Species-Specific",
                        species = "Meadow vole",
                        contag = seq(13, 100, length.out = 200),
                        phi = truePhi[3,,4],
                        upper = 0.5, #truePhi[4,,4],
                        lower = 0.5), #truePhi[2,,4]),
             data.frame(group = "Species-Specific",
                        species = "Deer mouse spp.",
                        contag = seq(13, 100, length.out = 200),
                        phi = truePhi[3,,5],
                        upper = 0.5, #truePhi[4,,5],
                        lower = 0.5), #truePhi[2,,5]),
             data.frame(group = "Species-Specific",
                        species = "Harvest mouse",
                        contag = seq(13, 100, length.out = 200),
                        phi = truePhi[3,,6],
                        upper = 0.5, #truePhi[4,,6],
                        lower = 0.5), #truePhi[2,,6]),
             data.frame(group = "Species-Specific",
                        species = "Jumping mouse",
                        contag = seq(13, 100, length.out = 200),
                        phi = truePhi[3,,7],
                        upper = 0.5, #truePhi[4,,7],
                        lower = 0.5)) #truePhi[2,,7]))
p2<- ggplot(CONTAG, aes(x=contag, y=phi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=species), alpha=.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=species)) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5,0.5,0.5,0.5,0.5,0.5), guide="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0.5,1), expand = c(0,0)) +
  labs(x="Contagion Index", y="", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank(), axis.text.x=element_text(size=7), axis.title=element_text(size=9),
        legend.title=element_blank(), legend.position=c(0.75,0.875), legend.background=element_blank(),
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))

jpeg("./results/phi.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2, ncol=2)
dev.off()

#####                             Recruitment (gamma[i,j,t])                               #####
##### CONTAG #####
predV <- seq(-1.33, 2.164393, length.out = 200)
# the predition variable matrix is a little different because of the 'season' covariate for
# recruitment. instead of holding that value at its mean (not sure what that would even be),
# we will compare recruitment at differing contagion values in different seasons
pred_mat <- cbind(
  1,
  predV,
  3      # transition from spring to summer 
)
pred_mat <- t(pred_mat)
pred_gamma <- array(
  NA,
  dim = c(
    nrow(mc$gamma0),
    ncol(pred_mat),
    7
  )
)
pred_gamma[,,1] <- exp(cbind(mc$mu.gamma0[,1],mc$mu.gamma1[,1],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,2] <- exp(cbind(mc$gamma0[,1],mc$gamma1[,1],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,3] <- exp(cbind(mc$gamma0[,2],mc$gamma1[,2],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,4] <- exp(cbind(mc$gamma0[,3],mc$gamma1[,3],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,5] <- exp(cbind(mc$gamma0[,4],mc$gamma1[,4],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,6] <- exp(cbind(mc$gamma0[,5],mc$gamma1[,5],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,7] <- exp(cbind(mc$gamma0[,6],mc$gamma1[,6],mc$gamma2[,1]) %*% pred_mat)

trueGamma <- apply(
  pred_gamma,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
ContagSP <- rbind(data.frame(group = "Community Average",
                           species = "Community Average",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,1],
                           upper = trueGamma[5,,1],
                           lower = trueGamma[1,,1]),
                data.frame(group = "Species-Specific",
                           species = "Short-tailed shrew",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,2],
                           upper = 0, #trueGamma[5,,2],
                           lower = 0), #trueGamma[1,,2]),
                data.frame(group = "Species-Specific",
                           species = "Prairie vole",
                            contag = seq(13, 100, length.out = 200),
                            gamma = trueGamma[3,,3],
                            upper = 0, #trueGamma[4,,3],
                            lower = 0), #trueGamma[2,,3]),
                data.frame(group = "Species-Specific",
                           species = "Meadow vole",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,4],
                           upper = 0, #trueGamma[5,,4],
                           lower = 0), #trueGamma[1,,4]),
                 data.frame(group = "Species-Specific",
                            species = "Deer mouse spp.",
                            contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,5],
                            upper = 0, #trueGamma[4,,5],
                            lower = 0), #trueGamma[2,,5]),
                data.frame(group = "Species-Specific",
                           species = "Harvest mouse",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,6],
                           upper = 0, #trueGamma[4,,6],
                           lower = 0), #trueGamma[2,,6]),
                data.frame(group = "Species-Specific",
                           species = "Jumping mouse",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,7],
                           upper = 0, #trueGamma[4,,7],
                           lower = 0)) #trueGamma[2,,7]))

p3 <- ggplot(ContagSP, aes(x=contag, y=gamma)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=group), show.legend=FALSE) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5), guide="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  labs(x="", y="", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=7), axis.text.y=element_blank())

pred_mat <- cbind(
  1,
  predV,
  1      # transition from summer to fall
)
pred_mat <- t(pred_mat)
pred_gamma <- array(
  NA,
  dim = c(
    nrow(mc$gamma0),
    ncol(pred_mat),
    7
  )
)
pred_gamma[,,1] <- exp(cbind(mc$mu.gamma0[,1],mc$mu.gamma1[,1],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,2] <- exp(cbind(mc$gamma0[,1],mc$gamma1[,1],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,3] <- exp(cbind(mc$gamma0[,2],mc$gamma1[,2],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,4] <- exp(cbind(mc$gamma0[,3],mc$gamma1[,3],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,5] <- exp(cbind(mc$gamma0[,4],mc$gamma1[,4],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,6] <- exp(cbind(mc$gamma0[,5],mc$gamma1[,5],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,7] <- exp(cbind(mc$gamma0[,6],mc$gamma1[,6],mc$gamma2[,1]) %*% pred_mat)
trueGamma <- apply(
  pred_gamma,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
ContagSU <- rbind(data.frame(group = "Community Average",
                           species = "Community Average",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,1],
                           upper = trueGamma[5,,1],
                           lower = trueGamma[1,,1]),
                data.frame(group = "Species-Specific",
                           species = "Short-tailed shrew",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,2],
                           upper = 0, #trueGamma[5,,2],
                           lower = 0), #trueGamma[1,,2]),
                data.frame(group = "Species-Specific",
                           species = "Prairie vole",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,3],
                           upper = 0, #trueGamma[4,,3],
                           lower = 0), #trueGamma[2,,3]),
                data.frame(group = "Species-Specific",
                           species = "Meadow vole",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,4],
                           upper = 0, #trueGamma[5,,4],
                           lower = 0), #trueGamma[1,,4]),
                data.frame(group = "Species-Specific",
                           species = "Deer mouse spp.",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,5],
                           upper = 0, #trueGamma[4,,5],
                           lower = 0), #trueGamma[2,,5]),
                data.frame(group = "Species-Specific",
                           species = "Harvest mouse",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,6],
                           upper = 0, #trueGamma[4,,6],
                           lower = 0), #trueGamma[2,,6]),
                data.frame(group = "Species-Specific",
                           species = "Jumping mouse",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,7],
                           upper = 0, #trueGamma[4,,7],
                           lower = 0)) #trueGamma[2,,7]))

p4 <- ggplot(ContagSU, aes(x=contag, y=gamma)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=group)) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5), guide="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  labs(x="", y="", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank(), axis.text.x = element_text(size=7),
        axis.title=element_text(size=9),
        legend.title=element_blank(), legend.position=c(0.6,0.875), legend.background=element_blank(),
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))

pred_mat <- cbind(
  1,
  predV,
  2      # transition from fall to next spring
)
pred_mat <- t(pred_mat)
pred_gamma <- array(
  NA,
  dim = c(
    nrow(mc$gamma0),
    ncol(pred_mat),
    7
  )
)
pred_gamma[,,1] <- exp(cbind(mc$mu.gamma0[,1],mc$mu.gamma1[,1],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,2] <- exp(cbind(mc$gamma0[,1],mc$gamma1[,1],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,3] <- exp(cbind(mc$gamma0[,2],mc$gamma1[,2],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,4] <- exp(cbind(mc$gamma0[,3],mc$gamma1[,3],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,5] <- exp(cbind(mc$gamma0[,4],mc$gamma1[,4],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,6] <- exp(cbind(mc$gamma0[,5],mc$gamma1[,5],mc$gamma2[,1]) %*% pred_mat)
pred_gamma[,,7] <- exp(cbind(mc$gamma0[,6],mc$gamma1[,6],mc$gamma2[,1]) %*% pred_mat)

trueGamma <- apply(
  pred_gamma,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

## Graphing
ContagFA <- rbind(data.frame(group = "Community Average",
                           species = "Community Average",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,1],
                           upper = trueGamma[5,,1],
                           lower = trueGamma[1,,1]),
                data.frame(group = "Species-Specific",
                           species = "Short-tailed shrew",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,2],
                           upper = 0, #trueGamma[5,,2],
                           lower = 0), #trueGamma[1,,2]),
                data.frame(group = "Species-Specific",
                           species = "Prairie vole",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,3],
                           upper = 0, #trueGamma[4,,3],
                           lower = 0), #trueGamma[2,,3]),
                data.frame(group = "Species-Specific",
                           species = "Meadow vole",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,4],
                           upper = 0, #trueGamma[5,,4],
                           lower = 0), #trueGamma[1,,4]),
                data.frame(group = "Species-Specific",
                           species = "Deer mouse spp.",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,5],
                           upper = 0, #trueGamma[4,,5],
                           lower = 0), #trueGamma[2,,5]),
                data.frame(group = "Species-Specific",
                           species = "Harvest mouse",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,6],
                           upper = 0, #trueGamma[4,,6],
                           lower = 0), #trueGamma[2,,6]),
                data.frame(group = "Species-Specific",
                           species = "Jumping mouse",
                           contag = seq(13, 100, length.out = 200),
                           gamma = trueGamma[3,,7],
                           upper = 0, #trueGamma[4,,7],
                           lower = 0)) #trueGamma[2,,7]))

p5 <- ggplot(ContagFA, aes(x=contag, y=gamma)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=group), show.legend = FALSE) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5), guide="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  labs(x="", y="Recruitment (individuals)", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))

jpeg("./results/gamma.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p5,p3,p4, ncol=3,
                        bottom=textGrob("Contagion Index", gp=gpar(fontsize=9)))
dev.off()

#####                                Supplemental Figures                                  #####
#####                               Detection (p[i,j,k,t])                                 #####
##### Moon Illumination #####
# generate a sequence of moon illumination values
# these roughly equate to 0% through 100% illumated
predV <- seq(-1.343997, 1.474745, length.out = 200)

# Generating a matrix that will contain values for predicting
pred_mat <- cbind(
  1,
  predV,
  0,     # mean of jDate
  0      # mean of effort
)
pred_mat <- t(pred_mat)
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat),
    7
  )
)
# fill it in
pred_p[,,1] <- cbind(mc$mu.alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,2] <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,3] <- cbind(mc$alpha0[,2],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,4] <- cbind(mc$alpha0[,3],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,5] <- cbind(mc$alpha0[,4],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,6] <- cbind(mc$alpha0[,5],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,7] <- cbind(mc$alpha0[,6],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

# convert to probability
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
moon <- rbind(data.frame(group = "Community Average",
                         species = "Community Average",
                         moon = seq(0, 1, length.out = 200),
                         p = prob_p[3,,1],
                         upper = prob_p[5,,1],
                         lower = prob_p[1,,1]),
                data.frame(group = "Species-Specific",
                           species = "Short-tailed shrew",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,2],
                           upper = 0, #prob_p[5,,2],
                           lower = 0), #prob_p[1,,2]),
                data.frame(group = "Species-Specific",
                           species = "Prairie vole",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,3],
                           upper = 0, #prob_p[5,,3],
                           lower = 0), #prob_p[1,,3]),
                data.frame(group = "Species-Specific",
                           species = "Meadow vole",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,4],
                           upper = 0, #prob_p[5,,4],
                           lower = 0), #prob_p[1,,4]),
                data.frame(group = "Species-Specific",
                           species = "Deer mouse spp.",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,5],
                           upper = 0, #prob_p[5,,5],
                           lower = 0), #prob_p[1,,5]),
                data.frame(group = "Species-Specific",
                           species = "Harvest mouse",
                           moon = seq(0, 1, length.out = 200),
                           p = prob_p[3,,6],
                           upper = 0, #prob_p[5,,6],
                           lower = 0), #prob_p[1,,6]),
              data.frame(group = "Species-Specific",
                         species = "Jumping mouse",
                         moon = seq(0, 1, length.out = 200),
                         p = prob_p[3,,7],
                         upper = 0, #prob_p[5,,7],
                         lower = 0)) #prob_p[1,,7]))

library(scales)
p6 <- ggplot(moon, aes(x=moon, y=p)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=group), show.legend=FALSE) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5)) +
  scale_x_continuous(labels = label_number(accuracy = 0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="Moon Illumination\n(proportion full)", y="Capture Probability", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9))

##### Julian Date #####
# generate a sequence of Julian date values that overlap with trapping timeframe
predV <- seq(-1.6, 1.6, length.out = 200)

# Generating a matrix that will contain values for predicting
pred_mat <- cbind(
  1,
  0,     # mean of scale(moon)
  predV,
  0      # mean of scale(effort)
)
pred_mat <- t(pred_mat)
# create the matrix
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat),
    7
  )
)
# fill it in
pred_p[,,1] <- cbind(mc$mu.alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,2] <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,3] <- cbind(mc$alpha0[,2],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,4] <- cbind(mc$alpha0[,3],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,5] <- cbind(mc$alpha0[,4],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,6] <- cbind(mc$alpha0[,5],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,7] <- cbind(mc$alpha0[,6],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

# convert to probability
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
# I converted jDate back into whole numbers for ease of interpretation
jDate <- rbind(data.frame(group = "Community Average",
                         species = "Community Average",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,1],
                         upper = prob_p[5,,1],
                         lower = prob_p[1,,1]),
              data.frame(group = "Species-Specific",
                         species = "Short-tailed shrew",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,2],
                         upper = 0, #prob_p[5,,2],
                         lower = 0), #prob_p[1,,2]),
              data.frame(group = "Species-Specific",
                         species = "Prairie vole",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,3],
                         upper = 0, #prob_p[5,,3],
                         lower = 0), #prob_p[1,,3]),
              data.frame(group = "Species-Specific",
                         species = "Meadow vole",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,4],
                         upper = 0, #prob_p[5,,4],
                         lower = 0), #prob_p[1,,4]),
              data.frame(group = "Species-Specific",
                         species = "Deer mouse spp.",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,5],
                         upper = 0, #prob_p[5,,5],
                         lower = 0), #prob_p[1,,5]),
              data.frame(group = "Species-Specific",
                         species = "Harvest mouse",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,6],
                         upper = 0, #prob_p[5,,6],
                         lower = 0), #prob_p[1,,6]),
              data.frame(group = "Species-Specific",
                         species = "Jumping mouse",
                         jDate = seq(101, 310, length.out = 200),
                         p = prob_p[3,,7],
                         upper = 0, #prob_p[5,,7],
                         lower = 0)) #prob_p[1,,7]))
p7 <- ggplot(jDate, aes(x=jDate, y=p)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values = c("#dddddd","#FFFFFF")) +
  geom_line(aes(color=species, linewidth=group), show.legend=FALSE) +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_linewidth_manual(values = c(1,0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="Julian Date\n ", y="", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9))

##### Effort #####
# the range of the 'effort' covariate is
# [1] -2.544777  1.402502
# however, rather than use a continuous scale, we need discrete counts of traps.
# I centered and scaled values from 1 - 10 "by hand" to use as predictors based on the mean
# and standard deviation from my dataset
predV <- c(-2.5447767, -2.1500488, -1.7553209, -1.3605930, -0.9658651, 
           -0.5711372, -0.1764093, 0.2183186, 0.6130466, 1.0077745, 1.4025024)
# these are basically the centered & scaled versions of 0,1,2,3,4,5,6,7,8,9,10

# Generating a matrix that will contain values for predicting
pred_mat <- cbind(
  1,
  0,     # mean of scale(moon)
  0,     # mean of scale(temp)
  predV
)
pred_mat <- t(pred_mat)
# create the matrix of the coefficients from the model
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat),
    7
  )
)
# fill it in
pred_p[,,1] <- cbind(mc$mu.alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,2] <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,3] <- cbind(mc$alpha0[,2],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,4] <- cbind(mc$alpha0[,3],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,5] <- cbind(mc$alpha0[,4],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,6] <- cbind(mc$alpha0[,5],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat
pred_p[,,7] <- cbind(mc$alpha0[,6],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

# convert to probability
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
effort <- rbind(data.frame(group = "Community Average",
                           species = "Community Average",
                           effort = c(0,1,2,3,4,5,6,7,8,9,10),
                           p = prob_p[3,,1],
                           upper = prob_p[5,,1],
                           lower = prob_p[1,,1]),
  data.frame(group = "Individual Species",
                          species = "Short-tailed shrew",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,2],
                         upper = prob_p[3,,2],
                         lower = prob_p[3,,2]),
              data.frame(group = "Individual Species",
                species = "Prairie vole",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,3],
                         upper = prob_p[3,,3],
                         lower = prob_p[3,,3]),
              data.frame(group = "Individual Species",
                species = "Meadow vole",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,4],
                         upper = prob_p[3,,4],
                         lower = prob_p[3,,4]),
              data.frame(group = "Individual Species",
                species = "Deer mouse spp.",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,5],
                         upper = prob_p[3,,5],
                         lower = prob_p[3,,5]),
              data.frame(group = "Individual Species",
                species = "Harvest mouse",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,6],
                         upper = prob_p[3,,6],
                         lower = prob_p[3,,6]),
              data.frame(group = "Individual Species",
                species = "Jumping mouse",
                         effort = c(0,1,2,3,4,5,6,7,8,9,10),
                         p = prob_p[3,,7],
                         upper = prob_p[3,,7],
                         lower = prob_p[3,,7]))
effort<-arrange(effort, desc(species))

p8<-ggplot(effort, aes(x=effort, y=p, color=species, size=group)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), show.legend=FALSE) +
  geom_point() +
  scale_color_manual(values = c("#000000","#77aadd","#ee8866","#eedd88","#ffaabb","#99ddff","#44bb99")) +
  scale_size_manual(values = c(0.6,0.3), guide="none") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  labs(x="Trap Effort\n ", y="", size="", color="Group") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.title=element_blank(), legend.position=c(0.6,0.875), legend.background=element_blank(),
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))

jpeg("./results/p.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p6,p7,p8, ncol=3)
dev.off()
