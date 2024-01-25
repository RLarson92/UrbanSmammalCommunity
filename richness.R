library(ggplot2)
library(coda)
library(vegan)
setwd("~/Documents/For School/Iowa/Dissertation Research/Small Mammals/SMamm Community")
functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}
my_mod <- readRDS("./results/colExt_mod.RDS")
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

data<-data.frame(cbind(covs_for_model$PC1,Rich))
jpeg("./results/richness.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(data = data, aes(x=V1,y=Rich)) +
  geom_point(color="#5f5f5f")+
  geom_smooth(method = 'lm', color='#000000') +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(limits=c(1,4), expand=c(0,0)) +
  labs(x="PC1", y="Species Richness") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7), axis.text.x=element_text(size=7), 
        axis.title=element_text(size=9))
dev.off()

#### Shannon's Diversity Index? ####
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
D <- vegan::diversity(comms, index = "simpson")
# models
model <- lm(H ~ covs_for_model$PC1 + covs_for_model$contag)
summary(model)
model <- lm(D ~ covs_for_model$PC1 + covs_for_model$contag)
summary(model)

data<-rbind(data.frame(x = covs_for_model$PC1,
                       y = H,
                       g = "Shannon's"),
            data.frame(x = covs_for_model$PC1,
                       y = D,
                       g = "Simpson's"))
# jpeg("./results/richness.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(data = data, aes(x=x,y=y,fill=g,color=g)) +
  geom_point()+
  geom_smooth(method = 'lm') +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="PC1", y="Index") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7), axis.text.x=element_text(size=7), 
        axis.title=element_text(size=9))

#### Average Abundance? ####
tmp <- mc$N
# determine average abundance of each species across time
tmp_avg <- apply(
  tmp,
  c(2,3),
  mean
)
comms <- t(tmp_avg)
pc1 <- rep(covs_for_model$PC1, 6)
species <- c(rep("BLBR",))

ggplot(data = comms, x=PC1) +
  geom_bar(position="stack") +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0,0)) +
  labs(x="PC1", y="Average Abundance (# of Individuals)") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7), axis.text.x=element_text(size=7), 
        axis.title=element_text(size=9))
