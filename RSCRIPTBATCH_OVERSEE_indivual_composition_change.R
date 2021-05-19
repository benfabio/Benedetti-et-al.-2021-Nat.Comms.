
# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("RColorBrewer")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("geosphere")
library("rgdal")
library("raster")
library("betapart")

world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------

### Set the working directories, vectors etc.
WD <- getwd()
rcp <- "rcp85"
# Vectors
SDMs <- c('GAM','GLM','RF','ANN')
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
# Vector of earth system models
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")

# --------------------------------------------------------------------------------------------------------------------------------

### Need to compute some centroids (baseline and future) to derive a species centroids shift. First develop a test code based on a susbet of the simulated communities, for zooplankton (phytoplankton communities still being extracted)

setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# dir()[grep(".Rdata",dir())]
esm <- "CESM-BEC"
sdm <- "GAM"
p <- "p2"

for(sdm in SDMs) {
    
        message(paste("Extracting annual compositions for ",sdm, sep = ""))
        
        for(p in pools) {
            
            message(paste("based on predictors of ",p, sep = ""))
            base.phyto <- get(load(paste("table_ann_compo_phyto_baseline_",sdm,"_",p,".Rdata", sep = "")))
            base.zoo <- get(load(paste("table_ann_compo_zoo_baseline_",sdm,"_",p,".Rdata", sep = "")))
            commons.base.tot <- intersect(unique(base.phyto$cell_id), unique(base.zoo$cell_id)) # length(commons.base)
            base <- cbind(base.phyto[which(base.phyto$cell_id %in% commons.base.tot),], base.zoo[which(base.zoo$cell_id %in% commons.base.tot),c(4:length(base.zoo))])
            
            for(esm in ESMs) {
                
                message(paste("and getting future projections for ",esm, sep = ""))
                message(paste("",sep = ""))
                fut.phyto <- get(load(paste("table_ann_compo_phyto_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = "")))
                fut.zoo <- get(load(paste("table_ann_compo_zoo_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = "")))
                # Identify commons cells
                commons.fut.tot <- intersect(unique(fut.phyto$cell_id), unique(fut.zoo$cell_id)) # length(commons.fut)
                commons.phyto <- intersect(unique(base.phyto$cell_id), unique(fut.phyto$cell_id)) # length(commons.phyto)
                commons.zoo <- intersect(unique(base.zoo$cell_id), unique(fut.zoo$cell_id)) # length(commons.zoo)
                # Bind future comp
                fut <- cbind(fut.phyto[which(fut.phyto$cell_id %in% commons.fut.tot),], fut.zoo[which(fut.zoo$cell_id %in% commons.fut.tot),c(4:length(fut.zoo))])
                commons <- intersect(unique(base$cell_id), unique(fut$cell_id)) # length(commons)
                
                ### Compute changes in annual composition
                if(sdm %in% c("GLM","GAM","ANN")) {
                    thresholds <- seq(from = 0.25, to = 0.39, by = 0.01) 
                } else {
                    thresholds <- seq(from = 0.10, to = 0.24, by = 0.01) 
                } # eo if else loop for thresholds
                
                # t <- 0.35
                thresh <- lapply(thresholds, function(t) {

                              message(paste("Using t = ",t, sep = ""))
                              # Convert to P/1 using rbinom and initial probability
                              for(sp in colnames(base2)[c(4:length(base2))] ) {
                                      base2[,c(sp)][base2[,c(sp)] > t] <- 1
                                      base2[,c(sp)][base2[,c(sp)] <= t] <- 0
                                      # Future 2 now
                                      fut2[,c(sp)][fut2[,c(sp)] > t] <- 1
                                      fut2[,c(sp)][fut2[,c(sp)] <= t] <- 0
                              } # eo for loop
                             
                              # Return a list
                              list.div <- list(div, div.phyto, div.zoo)
                              rm(div, div.zoo, div.phyto); gc()
                              return(list.div)

                      } # eo 2nd FUN

                )  # eo 2nd lapply - annual
                
            } # eo 3rd for loop - esm in ESMs
            
        } # eo 2nd for loop - p in pools
    
} # eo 1st for loop - sdm in SDMs


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

    