library(raster)
library(landscapemetrics)

rasterList <- list.files("Z://Rachel_SmallMammals/ArcMap/Connectivity_Landscapes/TIFFS/", full.names = TRUE)
rasterList <- rasterList[-46]

allRAsters <- lapply(rasterList, raster)

resultsList_contag <- list()

for(i in 1:length(allRAsters)){
  resultsList_contag[[i]] <- lsm_l_contag(allRAsters[[i]])
}

results_contag <- rlist::list.rbind(resultsList_contag)
results_contag$site <- rasterList

write.csv(results_contag, "Z://Rachel_SmallMammals/R/PlotConnectivity_Contag.csv", row.names = FALSE)

resultsList_aggri <- list()

for(i in 1:length(allRAsters)){
  resultsList_aggri[[i]] <- lsm_c_ai(allRAsters[[i]])
}

results_aggri <- rlist::list.rbind(resultsList_aggri)
results_aggri <- subset(results_aggri, class!=0)
results_aggri$site <- rasterList

write.csv(results_aggri, "Z://Rachel_SmallMammals/R/PlotConnectivity_aggri.csv", row.names = FALSE)
