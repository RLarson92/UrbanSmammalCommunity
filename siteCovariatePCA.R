siteCovs <- read.csv("~/Documents/For School/Iowa/Dissertation Research/Small Mammals/SMamm Community/data/siteCovs.csv")
siteCovs <- scale(siteCovs[,2:5])
cor(siteCovs)
# canopy X shrub; not correlated p = 0.93
cor.test(siteCovs[,1], siteCovs[,2],
         alternative = "two.sided",
         method = "spearman")

######## canopy X humanMod; correlated negatively p < 0.01
cor.test(siteCovs[,1], siteCovs[,4],
         alternative = "two.sided",
         method = "spearman")

######## shrub X humanMod; correlated negatively p = 0.05
cor.test(siteCovs[,2], siteCovs[,4],
         alternative = "two.sided",
         method = "spearman")


pca <- princomp(siteCovs)
summary(pca)
pca$loadings[, 1:3]
pca$scores

library(ggfortify)
autoplot(pca, data = siteCovs,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
