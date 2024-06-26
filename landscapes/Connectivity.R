###########################################################################################################
#                                                                                                         #
#                                 Code for calculating contagion index values                             #
#                                          using `landscapemetrics`                                       #
#                                                                                                         #
#                                         Last edited: 6 Jun 2024                                        #
#                                                                                                         #
###########################################################################################################
# You'll need to load the 'raster' and 'landscapemetrics' packages to perform this analysis
library(raster)
library(landscapemetrics)

# Code to list all the files in the '.tifs' subfolder and import them in the R workspace
rasterList <- list.files("./landscapes/tifs/", full.names = TRUE)
rasterList <- rasterList[-46]  # this is some Windows file not part of the dataset

# telling R that all the files in rasterList are, in fact, rasters
allRasters <- lapply(rasterList, raster)

# creating an empty list item to store the results of the analysis in
resultsList_contag <- list()

# automating the calculation the contagion values for all 45 sites
for(i in 1:length(allRasters)){
  resultsList_contag[[i]] <- lsm_l_contag(allRasters[[i]])
}
# a couple of the trap sites contain only 1 land cover class, and will cause a warning
# from landscapemetrics. When inspected visually in ArcMap or the like, you'll see their
# contagion value should be 0 (i.e., low interspersion and high "clumping"). You'll have
# to tell R to replace the resulting 'NA's with '0'

# writing the results into that blank list and adding the site name
results_contag <- rlist::list.rbind(resultsList_contag)
results_contag$site <- c(1,11,111,117,12,120,122,123,124,135,136,140,142,143,147,148,161,162,17,
                         205,214,216,221,222,223,225,227,229,23,245,246,250,255,256,260,264,269,274,288,289,295,
                         8,83,88,90)
# subsetting just the contagion values
results_contag <- results_contag[,6:7]
results_contag$value[is.na(results_contag$value)] <- 0
colnames(results_contag) <- c("contag","site")

# writing the results to a .csv
write.csv(results_contag, "./data/contag.csv", row.names = FALSE)
