library(ggplot2)
library(coda)
library(vegan)

functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}
my_mod <- readRDS("./results/my_mod.RDS")
covs_for_model <- read.csv("./data/covs_for_model.csv", stringsAsFactors = FALSE)

mc <- coda::as.mcmc(my_mod)
# convert the mcmc object to a matrix
mc <- as.matrix(mc, iters = FALSE)
# sub-sample the mcmc matrix a bit as we don't really need to make predictions with 
# all 60K samples
set.seed(554)
mc_sub <- mc[sample(1:nrow(mc), 10000), ]
# and use split_mcmc
mc <- split_mcmc(mc_sub)
rm(mc_sub)

#### Species Richness ####
# collapse abundance estimates into occurrence data
tmp <- mc$N
tmp[tmp>0] <- 1
# determine average occurrence probability across time
tmp_avg <- apply(
  tmp,
  c(2,3),
  mean
)
# sum occurrence probabilities across species
Rich <- apply(
  tmp_avg,
  2,
  sum
)
range(Rich)
mean(Rich)
# using the hist function to check for normality
hist(Rich)

model <- lm(Rich ~ covs_for_model$PC1 + covs_for_model$contag)
summary(model)

#### Shannon's Diversity Index ####
tmp <- mc$N
# determine average abundance of each species across time
tmp_avg <- apply(
  tmp,
  c(2,3),
  mean
)
comms <- t(tmp_avg)
# calculate indexes;
H <- vegan::diversity(comms, index = "shannon")
# model
model <- lm(H ~ covs_for_model$PC1 + covs_for_model$contag)
summary(model)
