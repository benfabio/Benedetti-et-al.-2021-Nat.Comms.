
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

world2 <- map_data("world2")
world <- map_data("world")

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

# setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# dir()[grep(".Rdata",dir())]
# esm <- "IPSL-PISCES"
# sdm <- "GAM"
# p <- "p1"
#
# base.phyto <- get(load(paste("table_ann_compo_phyto_baseline_",sdm,"_",p,".Rdata", sep = "")))
# base.zoo <- get(load(paste("table_ann_compo_zoo_baseline_",sdm,"_",p,".Rdata", sep = "")))
# fut.phyto <- get(load(paste("table_ann_compo_phyto_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = "")))
# fut.zoo <- get(load(paste("table_ann_compo_zoo_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = "")))
# dim(base.phyto); dim(base.zoo); dim(fut.phyto); dim(fut.zoo)
# colnames(fut.phyto)
#
# # Cbind according to their common cells
# commons.base.tot <- intersect(unique(base.phyto$cell_id), unique(base.zoo$cell_id)) # length(commons.base)
# commons.fut.tot <- intersect(unique(fut.phyto$cell_id), unique(fut.zoo$cell_id)) # length(commons.fut)
# # And per group
# commons.phyto <- intersect(unique(base.phyto$cell_id), unique(fut.phyto$cell_id)) # length(commons.phyto)
# commons.zoo <- intersect(unique(base.zoo$cell_id), unique(fut.zoo$cell_id)) # length(commons.zoo)
#
# base <- cbind(base.phyto[which(base.phyto$cell_id %in% commons.base.tot),], base.zoo[which(base.zoo$cell_id %in% commons.base.tot),c(4:length(base.zoo))])
# fut <- cbind(fut.phyto[which(fut.phyto$cell_id %in% commons.fut.tot),], fut.zoo[which(fut.zoo$cell_id %in% commons.fut.tot),c(4:length(fut.zoo))])
# # dim(base); dim(fut)
# # And commons between base and fut
# commons <- intersect(unique(base$cell_id), unique(fut$cell_id)) ; length(commons)
#
# # Computing species range centroid at T0: weighted average lon and lat
# m.base.zoo <- melt(base.zoo[base.zoo$cell_id %in% commons.zoo,], id.vars = c("cell_id","x","y") )
# m.fut.zoo <- melt(fut.zoo[fut.zoo$cell_id %in% commons.zoo,], id.vars = c("cell_id","x","y") )
# # dim(m.base.zoo); dim(m.fut.zoo)
# colnames(m.base.zoo)[c(4,5)] <- c("species","HSI")
# colnames(m.fut.zoo)[c(4,5)] <- c("species","HSI")
#
# # Move longitudes from -180°/180° # unique(m.base.zoo$x)
# #m.base.phyto$x2 <- m.base.phyto$x
# m.base.zoo$x2 <- m.base.zoo$x
# #m.fut.phyto$x2 <- m.fut.phyto$x
# m.fut.zoo$x2 <- m.fut.zoo$x
# #m.base.phyto[which(m.base.phyto$x > 180),c("x2")] <- m.base.phyto[which(m.base.phyto$x > 180),"x"] - 360
# m.base.zoo[which(m.base.zoo$x > 180),c("x2")] <- m.base.zoo[which(m.base.zoo$x > 180),"x"] - 360
# #m.fut.phyto[which(m.fut.phyto$x > 180),c("x2")] <- m.fut.phyto[which(m.fut.phyto$x > 180),"x"] - 360
# m.fut.zoo[which(m.fut.zoo$x > 180),c("x2")] <- m.fut.zoo[which(m.fut.zoo$x > 180),"x"] - 360
# # summary(m.base.zoo$x2); summary(m.fut.phyto$x2)
#
#
# # Use dplyr
# require("dplyr")
# centroids.base <- data.frame(m.base.zoo %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = abs(y), w = HSI) ) )
# centroids.fut <- data.frame(m.fut.zoo %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = abs(y), w = HSI) ) )
# # N Hemis
# centroids.base.NH <- data.frame(m.base.zoo[m.base.zoo$y >= 0,] %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
# centroids.fut.NH <- data.frame(m.fut.zoo[m.base.zoo$y >= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
# # S Hemis
# centroids.base.SH <- data.frame(m.base.zoo[m.base.zoo$y <= 0,] %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
# centroids.fut.SH <- data.frame(m.fut.zoo[m.base.zoo$y <= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
# # Check
# summary(centroids.base); dim(centroids.base)
# summary(centroids.fut); dim(centroids.fut)
# summary(centroids.base.NH); summary(centroids.base.SH)
# # Bind
# centroids <- data.frame(species = centroids.fut$species,
#          x_base = centroids.base$lon, y_base = centroids.base$lat,
#          x_fut = centroids.fut$lon, y_fut = centroids.fut$lat,
#          x_base_NH = centroids.base.NH$lon, y_base_NH = centroids.base.NH$lat,
#          x_fut_NH = centroids.fut.NH$lon, y_fut_NH = centroids.fut.NH$lat,
#          x_base_SH = centroids.base.SH$lon, y_base_SH = centroids.base.SH$lat,
#          x_fut_SH = centroids.fut.SH$lon, y_fut_SH = centroids.fut.SH$lat
# )
# # summary(centroids); dim(centroids)
# # Compute distance in km
# centroids$distm <- NA
# centroids$distm_NH <- NA
# centroids$distm_SH <- NA
# for(i in 1:nrow(centroids)) {
#      centroids[i,"distm"] <- distm(centroids[i,c("x_base","y_base")], centroids[i,c("x_fut","y_fut")], fun = distHaversine)
#      centroids[i,"distm_NH"] <- distm(centroids[i,c("x_base_NH","y_base_NH")], centroids[i,c("x_fut_NH","y_fut_NH")], fun = distHaversine)
#      centroids[i,"distm_SH"] <- distm(centroids[i,c("x_base_SH","y_base_SH")], centroids[i,c("x_fut_SH","y_fut_SH")], fun = distHaversine)
# } # eo dist for loop
# # Convert to kilometers
# centroids$distm <- (centroids$distm)/1000
# centroids$distm_NH <- (centroids$distm_NH)/1000
# centroids$distm_SH <- (centroids$distm_SH)/1000
#
# centroids[order(centroids$distm, decreasing = T),c("species","distm")]
#
# # Try map
# ggplot() + geom_point(aes(x = x_base, y = y_base), data = centroids, colour = "#3288bd") +
#      geom_point(aes(x = x_fut, y = y_fut), data = centroids, colour = "#d53e4f") + coord_quickmap() +
#      theme_bw()
#
# #  For a species, map base or fut pattern in HSI with centroid on top (to make sure these centroids make sense!)
# # Need to shift coords again in.base or fut.zoo
# base.zoo[which(base.zoo$x > 180),c("x2")] <- base.zoo[which(base.zoo$x > 180),"x"] - 360
#
# sp <- "Oithona_simplex"
# #quartz()
# ggplot() + geom_raster(aes(x = x2, y = y, fill = get(sp)), data = base.zoo) +
#       geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#       scale_fill_viridis(name = "Annual HSI", limits = c(0,1)) +
#       geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x2, y = y, z = get(sp)), data = base.zoo) +
#       geom_point(aes(x = x_base, y = y_base), data = centroids[centroids$species == sp,], colour = "black", pch = 23, fill = "#d53e4f") +
#       geom_point(aes(x = x_base_NH, y = y_base_NH), data = centroids[centroids$species == sp,], colour = "black", pch = 21, fill = "#d53e4f") +
#       geom_point(aes(x = x_base_SH, y = y_base_SH), data = centroids[centroids$species == sp,], colour = "black", pch = 22, fill = "#d53e4f") +
#          coord_quickmap() + theme_bw() + ylab("Latitude") + xlab("Longitude")
#
# #quartz()
# ggplot() + geom_raster(aes(x = x2, y = y, fill = get(sp)), data = base.zoo) +
#       geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x2, y = y, z = get(sp)), data = base.zoo) +
#       scale_fill_viridis(name = "Annual HSI", limits = c(0,1)) +
#       geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#       geom_point(aes(x = x_fut, y = y_fut), data = centroids[centroids$species == sp,], colour = "black", pch = 23, fill = "#d53e4f", size = 2) +
#       geom_point(aes(x = x_fut_NH, y = y_fut_NH), data = centroids[centroids$species == sp,], colour = "black", pch = 21, fill = "#d53e4f", size = 2) +
#       geom_point(aes(x = x_fut_SH, y = y_fut_SH), data = centroids[centroids$species == sp,], colour = "black", pch = 22, fill = "#d53e4f", size = 2) +
#       geom_point(aes(x = x_base, y = y_base), data = centroids[centroids$species == sp,], colour = "black", pch = 23, fill = "#3288bd", size = 2) +
#       geom_point(aes(x = x_base_NH, y = y_base_NH), data = centroids[centroids$species == sp,], colour = "black", pch = 21, fill = "#3288bd", size = 2) +
#       geom_point(aes(x = x_base_SH, y = y_base_SH), data = centroids[centroids$species == sp,], colour = "black", pch = 22, fill = "#3288bd", size = 2) +
#          coord_quickmap() + theme_bw() + ylab("Latitude") + xlab("Longitude") + ggtitle(sp)

# -------------------------------------------------------------------

### Incorporate those scripts in for loops or mclapply...maybe with a big mclapply based on file names
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# dir()[grep("table_ann_compo_",dir())]
#esm <- "IPSL-PISCES"
#sdm <- "RF"
#p <- "p2"

# base[base$x < 0 ,"x2"] <- (base[base$x < 0 ,"x"]) + 360
# --> x which are 180-360 --> -360

for(sdm in SDMs) {
    
        message(paste("Extracting annual compositions for ",sdm, sep = ""))
        
        for(p in pools) {
            
            setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
            message(paste("based on predictors of ",p, sep = ""))
            base.phyto <- get(load(paste("table_ann_compo_phyto_baseline_",sdm,"_",p,".Rdata", sep = "")))
            base.zoo <- get(load(paste("table_ann_compo_zoo_baseline_",sdm,"_",p,".Rdata", sep = "")))     
            
            for(esm in ESMs) {
                
                message(paste("and getting future projections for ",esm, sep = ""))
                fut.phyto <- get(load(paste("table_ann_compo_phyto_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = "")))
                fut.zoo <- get(load(paste("table_ann_compo_zoo_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = "")))
              
                # Cbind according to their common cells
                commons.base.tot <- intersect(unique(base.phyto$cell_id), unique(base.zoo$cell_id)) # length(commons.base)
                commons.fut.tot <- intersect(unique(fut.phyto$cell_id), unique(fut.zoo$cell_id)) # length(commons.fut)
                commons.phyto <- intersect(unique(base.phyto$cell_id), unique(fut.phyto$cell_id)) # length(commons.phyto)
                commons.zoo <- intersect(unique(base.zoo$cell_id), unique(fut.zoo$cell_id)) # length(commons.zoo)

                # Computing species range centroid at T0: weighted average lon and lat 
                m.base.zoo <- melt(base.zoo[base.zoo$cell_id %in% commons.zoo,], id.vars = c("cell_id","x","y") )
                m.fut.zoo <- melt(fut.zoo[fut.zoo$cell_id %in% commons.zoo,], id.vars = c("cell_id","x","y") )
                colnames(m.base.zoo)[c(4,5)] <- c("species","HSI")
                colnames(m.fut.zoo)[c(4,5)] <- c("species","HSI")
                m.base.phyto <- melt(base.phyto[base.phyto$cell_id %in% commons.phyto,], id.vars = c("cell_id","x","y") )
                m.fut.phyto <- melt(fut.phyto[fut.phyto$cell_id %in% commons.phyto,], id.vars = c("cell_id","x","y") )
                colnames(m.base.phyto)[c(4,5)] <- c("species","HSI")
                colnames(m.fut.phyto)[c(4,5)] <- c("species","HSI")
                
                # Move longitudes from -180°/180° # unique(m.base.zoo$x)
                m.base.phyto$x2 <- m.base.phyto$x
                m.base.zoo$x2 <- m.base.zoo$x
                m.fut.phyto$x2 <- m.fut.phyto$x
                m.fut.zoo$x2 <- m.fut.zoo$x
                m.base.phyto[which(m.base.phyto$x > 180),c("x2")] <- m.base.phyto[which(m.base.phyto$x > 180),"x"] - 360
                m.base.zoo[which(m.base.zoo$x > 180),c("x2")] <- m.base.zoo[which(m.base.zoo$x > 180),"x"] - 360
                m.fut.phyto[which(m.fut.phyto$x > 180),c("x2")] <- m.fut.phyto[which(m.fut.phyto$x > 180),"x"] - 360
                m.fut.zoo[which(m.fut.zoo$x > 180),c("x2")] <- m.fut.zoo[which(m.fut.zoo$x > 180),"x"] - 360
                # summary(m.base.zoo$x2); summary(m.fut.phyto$x2)
                  
                # Use dplyr to compute average coordinates
                require("dplyr")
                ### For zoo
                centroids.base.zoo <- data.frame(m.base.zoo %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = abs(y), w = HSI) ) )
                centroids.fut.zoo <- data.frame(m.fut.zoo %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = abs(y), w = HSI) ) )
                # N Hemis
                centroids.base.NH.zoo <- data.frame(m.base.zoo[m.base.zoo$y >= 0,] %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                centroids.fut.NH.zoo <- data.frame(m.fut.zoo[m.base.zoo$y >= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                # S Hemis
                centroids.base.SH.zoo <- data.frame(m.base.zoo[m.base.zoo$y <= 0,] %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                centroids.fut.SH.zoo <- data.frame(m.fut.zoo[m.base.zoo$y <= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
               
                ### For phyto
                centroids.base.phyto <- data.frame(m.base.phyto %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = abs(y), w = HSI) ) )
                centroids.fut.phyto <- data.frame(m.fut.phyto %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = abs(y), w = HSI) ) )
                # N Hemis
                centroids.base.NH.phyto <- data.frame(m.base.phyto[m.base.phyto$y >= 0,] %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                centroids.fut.NH.phyto <- data.frame(m.fut.phyto[m.base.phyto$y >= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                # S Hemis
                centroids.base.SH.phyto <- data.frame(m.base.phyto[m.base.phyto$y <= 0,] %>% group_by(species) %>%summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                centroids.fut.SH.phyto <- data.frame(m.fut.phyto[m.base.phyto$y <= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x2, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
                
                # Bind
                centroids.zoo <- data.frame(species = centroids.fut.zoo$species, 
                        x_base = centroids.base.zoo$lon, y_base = centroids.base.zoo$lat,
                        x_fut = centroids.fut.zoo$lon, y_fut = centroids.fut.zoo$lat,
                        x_base_NH = centroids.base.NH.zoo$lon, y_base_NH = centroids.base.NH.zoo$lat,
                        x_fut_NH = centroids.fut.NH.zoo$lon, y_fut_NH = centroids.fut.NH.zoo$lat, 
                        x_base_SH = centroids.base.SH.zoo$lon, y_base_SH = centroids.base.SH.zoo$lat,
                        x_fut_SH = centroids.fut.SH.zoo$lon, y_fut_SH = centroids.fut.SH.zoo$lat
                )
                centroids.phyto <- data.frame(species = centroids.fut.phyto$species, 
                        x_base = centroids.base.phyto$lon, y_base = centroids.base.phyto$lat,
                        x_fut = centroids.fut.phyto$lon, y_fut = centroids.fut.phyto$lat,
                        x_base_NH = centroids.base.NH.phyto$lon, y_base_NH = centroids.base.NH.phyto$lat,
                        x_fut_NH = centroids.fut.NH.phyto$lon, y_fut_NH = centroids.fut.NH.phyto$lat, 
                        x_base_SH = centroids.base.SH.phyto$lon, y_base_SH = centroids.base.SH.phyto$lat,
                        x_fut_SH = centroids.fut.SH.phyto$lon, y_fut_SH = centroids.fut.SH.phyto$lat
                )
                
                # Rbind
                centroids.phyto$group <- "Phytoplankton"
                centroids.zoo$group <- "Zooplankton"
                centroids <- rbind(centroids.phyto, centroids.zoo)
                # dim(centroids); summary(centroids)
                
                # Compute distance in km
                message(paste("Computing centroids shifts",sep = ""))
                centroids$distm <- NA
                centroids$distm_NH <- NA
                centroids$distm_SH <- NA
                centroids$distm_lat <- NA
                centroids$distm_lon <- NA
                centroids$distm_NH_lat <- NA
                centroids$distm_NH_lon <- NA
                centroids$distm_SH_lat <- NA
                centroids$distm_SH_lon <- NA
                
                # 26/03/2020: Also compute longitudonal ad latitudinal components of the distances:
                # - lat distance: keep x1 for future centroid (x1,y2)
                # - long distance: keep y1 for future centroid (x2,y1)
                for(i in 1:nrow(centroids)) {
                    
                    centroids[i,"distm"] <- distm(centroids[i,c("x_base","y_base")], centroids[i,c("x_fut","y_fut")], fun = distHaversine)
                    centroids[i,"distm_NH"] <- distm(centroids[i,c("x_base_NH","y_base_NH")], centroids[i,c("x_fut_NH","y_fut_NH")], fun = distHaversine)
                    centroids[i,"distm_SH"] <- distm(centroids[i,c("x_base_SH","y_base_SH")], centroids[i,c("x_fut_SH","y_fut_SH")], fun = distHaversine)
                    
                    centroids[i,"distm_lat"] <- distm(centroids[i,c("x_base","y_base")], centroids[i,c("x_base","y_fut")], fun = distHaversine)
                    centroids[i,"distm_NH_lat"] <- distm(centroids[i,c("x_base_NH","y_base_NH")], centroids[i,c("x_base_NH","y_fut_NH")], fun = distHaversine)
                    centroids[i,"distm_SH_lat"] <- distm(centroids[i,c("x_base_SH","y_base_SH")], centroids[i,c("x_base_SH","y_fut_SH")], fun = distHaversine)
                    
                    centroids[i,"distm_lon"] <- distm(centroids[i,c("x_base","y_base")], centroids[i,c("x_fut","y_base")], fun = distHaversine)
                    centroids[i,"distm_NH_lon"] <- distm(centroids[i,c("x_base_NH","y_base_NH")], centroids[i,c("x_fut_NH","y_base_NH")], fun = distHaversine)
                    centroids[i,"distm_SH_lon"] <- distm(centroids[i,c("x_base_SH","y_base_SH")], centroids[i,c("x_fut_NH","y_base_SH")], fun = distHaversine)
                    
                } # eo dist for loop
                # Convert to kilometers
                centroids$distm <- (centroids$distm)/1000
                centroids$distm_NH <- (centroids$distm_NH)/1000
                centroids$distm_SH <- (centroids$distm_SH)/1000
                centroids$distm_lat <- (centroids$distm_lat)/1000
                centroids$distm_NH_lat <- (centroids$distm_NH_lat)/1000
                centroids$distm_SH_lat <- (centroids$distm_SH_lat)/1000
                centroids$distm_lon <- (centroids$distm_lon)/1000
                centroids$distm_NH_lon <- (centroids$distm_NH_lon)/1000
                centroids$distm_SH_lon <- (centroids$distm_SH_lon)/1000
                # centroids[order(centroids$distm, decreasing = T),c("species","distm","distm_lat","distm_NH_lat")]
                # centroids[order(centroids$distm_lat, decreasing = T),c("species","distm","distm_lat","distm_NH_lat")]
                # summary(centroids)
                # summary( (centroids$distm_lat)/7 ) # 7 decades
                # summary( (centroids$distm)/7 ) # 7 decades
                # summary( (centroids$distm_NH_lat)/7 ) # 7 decades

                ### Save
                setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/species_shifts")   
                save(centroids, file = paste("table_species_shifts_",esm,"_",sdm,"_",p,".Rdata",sep = "") )   
                
                # Clean stuff and move to next level
                rm(centroids,centroids.zoo,centroids.phyto,centroids.fut.SH.phyto,centroids.base.SH.phyto,centroids.fut.NH.phyto,centroids.base.NH.phyto,
                    centroids.base.phyto,centroids.fut.phyto,centroids.fut.SH.zoo,centroids.base.SH.zoo,centroids.fut.NH.zoo,centroids.base.NH.zoo,
                    centroids.fut.zoo,centroids.base.zoo,m.fut.phyto,m.base.phyto,m.fut.zoo,m.base.zoo)
                gc()          
                
                setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
                message(paste("  ", sep = ""))

            } # eo 3rd for loop - esm in ESMs
            
        } # eo 2nd for loop - p in pools
    
} # eo 1st for loop - sdm in SDMs


# --------------------------------------------------------------------------------------------------------------------------------

### 26/03/2020: Summarize results from above (distrbution plot and PCA etc.)
library("vegan")
library("FactoMineR")
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/species_shifts") 
# dir()
### Rbind after a lapply, don't forget to specify the pool/ESM/SDM along the way
# f <- dir()[5]
res <- lapply(dir()[grep("table",dir())], function(f) {
            d <- get(load(f)) #  head(d)
            # Extract terms from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_")) # terms
            d$ESM <- terms[4,1]
            d$SDM <- terms[5,1]
            return(d)
    } # eo fun
) # eo lapply
# Rbind
ddf <- bind_rows(res)
dim(ddf); head(ddf)
rm(res); gc()
summary(ddf)
### Recall the variables' meaning:
# - distm = distance (km) between the 2 centroids based on absolute latitude
# - distm_lat = distance in latitudinal dimension (poleward of equatorward) between the 2 centroids based on absolute latitude
# - distm_lon = distance in longitudinal dimension (poleward of equatorward) between the 2 centroids based on absolute latitude

### Compute ensembles per SDM and ESM
ens <- data.frame(ddf %>% group_by(species) %>% summarize(group = unique(group),
            x_base = mean(x_base), y_base = mean(y_base), x_fut = mean(x_fut), y_fut = mean(y_fut),
            distm = mean(distm), distm_lat = mean(distm_lat), distm_lon = mean(distm_lon) )
) # eo ddf
summary(ens)
summary(ens[ens$group == "Phytoplankton",c("distm","distm_lat","distm_lon")])
summary(ens[ens$group == "Zooplankton",c("distm","distm_lat","distm_lon")])
# Convert to the shift speed (km/decade) by dividing by 7.
ens$shift <- (ens$distm)/7
ens$shift_lat <- (ens$distm_lat)/7
ens$shift_lon <- (ens$distm_lon)/7

ens[order(ens$shift_lat, decreasing = T),c("species","group","shift_lat")]

# Plot distrbution of per group
p1 <- ggplot(data = ens, aes(x = factor(group), y = shift, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Shift speed (km/dec)\nbased on the absolute centroid") + theme_classic() 
#
p2 <- ggplot(data = ens, aes(x = factor(group), y = shift_lat, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Latitudinal shift speed (km/dec)\nbased on the absolute centroid") + theme_classic() 
#
p3 <- ggplot(data = ens, aes(x = factor(group), y = shift_lon, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Longitudinal shift speed (km/dec)\nbased on the absolute centroid") + theme_classic() 

# Same with distance between centroids
ggplot(data = ens, aes(x = factor(group), y = distm, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Distance (km)\nbetween baseline and future centroid") + theme_classic() 
#
ggplot(data = ens, aes(x = factor(group), y = distm_lat, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Latitudinal distance (km)\nbetween baseline and future centroid") + theme_classic() 
#
ggplot(data = ens, aes(x = factor(group), y = distm_lon, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Longitudinal distance (km)\nbetween baseline and future centroid") + theme_classic() 
    
setwd(WD)
ggsave(plot = p1, filename = "plot_distrob_shift_absy_ensemble.jpg", dpi = 300, width = 4, height = 3)
ggsave(plot = p2, filename = "plot_distrob_shift_lat_absy_ensemble.jpg", dpi = 300, width = 4, height = 3)
ggsave(plot = p3, filename = "plot_distrob_shift_lon_absy_ensemble.jpg", dpi = 300, width = 4, height = 3)

# Add 2 factors that specify whetehr the shift is poleward/ equatorward and westward/ eastward
# REMEMBER: here coordibates are the common -180/+180 CRS 
ens$lat_dir <- NA
ens$lon_dir <- NA
for(i in c(1:nrow(ens)) ) {
    
    xbase <- ens[i,"x_base"]
    xfut <- ens[i,"x_fut"]
    ybase <- ens[i,"y_base"]
    yfut <- ens[i,"y_fut"]
    
    if(xbase >= xfut) {
        ens[i,"lon_dir"] <- "Eastward"
    } else {
        ens[i,"lon_dir"] <- "Westward"
    } # eo 1st if loop
    
    if(ybase >= yfut) {
        ens[i,"lat_dir"] <- "Equatorward"
    } else {
        ens[i,"lat_dir"] <- "Poleward"
    } # eo 1st if loop
    
} # for loop

### Check number per categories
nrow(ens[which(ens$lat_dir == "Poleward" & ens$lon_dir == "Eastward"),]) / nrow(ens)
nrow(ens[which(ens$lat_dir == "Poleward" & ens$lon_dir == "Westward"),]) / nrow(ens)
nrow(ens[which(ens$lat_dir == "Equatorward" & ens$lon_dir == "Eastward"),]) / nrow(ens)
nrow(ens[which(ens$lat_dir == "Equatorward" & ens$lon_dir == "Westward"),]) / nrow(ens)

median(ens[ens$lat_dir == "Poleward",c("shift_lat")]); IQR(ens[ens$lat_dir == "Poleward",c("shift_lat")])
median(ens[ens$lat_dir == "Equatorward",c("shift_lat")]); IQR(ens[ens$lat_dir == "Equatorward",c("shift_lat")])

median(ens[ens$lat_dir == "Poleward",c("shift")]); IQR(ens[ens$lat_dir == "Poleward",c("shift")])
median(ens[ens$lat_dir == "Equatorward",c("shift")]); IQR(ens[ens$lat_dir == "Equatorward",c("shift")])

### Separate phyto from zoo
median(ens[ens$lat_dir == "Poleward" & ens$group == "Phytoplankton",c("shift")]); IQR(ens[ens$lat_dir == "Poleward",c("shift")])
median(ens[ens$lat_dir == "Equatorward" & ens$group == "Phytoplankton",c("shift")]); IQR(ens[ens$lat_dir == "Equatorward",c("shift")])

median(ens[ens$lat_dir == "Poleward" & ens$group == "Phytoplankton",c("shift")]); IQR(ens[ens$lat_dir == "Poleward",c("shift")])
median(ens[ens$lat_dir == "Equatorward" & ens$group == "Phytoplankton",c("shift")]); IQR(ens[ens$lat_dir == "Equatorward",c("shift")])


# 34% of taxa shift Poleward and Eastward
# 45% of taxa shift Poleward and Westward
# 9% of taxa shift Equatorward and Eastward
# 12% of taxa shift Equatorward and Westward
# ---------------------------------------------
# 79% of taxa shift poleward/ 21% shift equatorward
# 43% of taxa shift eastward/ 57% shift westward

### Phyto vs. zoo
nrow(ens[which(ens$lat_dir == "Poleward" & ens$lon_dir == "Eastward" & ens$group == "Phytoplankton"),]) / nrow(ens[ens$group == "Phytoplankton",])
nrow(ens[which(ens$lat_dir == "Poleward" & ens$lon_dir == "Westward" & ens$group == "Phytoplankton"),]) / nrow(ens[ens$group == "Phytoplankton",])
nrow(ens[which(ens$lat_dir == "Equatorward" & ens$lon_dir == "Eastward" & ens$group == "Phytoplankton"),]) / nrow(ens[ens$group == "Phytoplankton",])
nrow(ens[which(ens$lat_dir == "Equatorward" & ens$lon_dir == "Westward" & ens$group == "Phytoplankton"),]) / nrow(ens[ens$group == "Phytoplankton",])
# 26% of taxa shift Poleward and Eastward
# 41% of taxa shift Poleward and Westward
### --> 67% phyto species migrate polewards !
# 16% of taxa shift Equatorward and Eastward
# 16% of taxa shift Equatorward and Westward

nrow(ens[which(ens$lat_dir == "Poleward" & ens$lon_dir == "Eastward" & ens$group == "Zooplankton"),]) / nrow(ens[ens$group == "Zooplankton",])
nrow(ens[which(ens$lat_dir == "Poleward" & ens$lon_dir == "Westward" & ens$group == "Zooplankton"),]) / nrow(ens[ens$group == "Zooplankton",])
nrow(ens[which(ens$lat_dir == "Equatorward" & ens$lon_dir == "Eastward" & ens$group == "Zooplankton"),]) / nrow(ens[ens$group == "Zooplankton",])
nrow(ens[which(ens$lat_dir == "Equatorward" & ens$lon_dir == "Westward" & ens$group == "Zooplankton"),]) / nrow(ens[ens$group == "Zooplankton",])
# 40% of taxa shift Poleward and Eastward
# 47% of taxa shift Poleward and Westward
### --> 87% zoo species migrate polewards !
# 5% of taxa shift Equatorward and Eastward
# 8% of taxa shift Equatorward and Westward


### Phyto lat shfts vs zoo lat shifts
summary(ens[ens$group == "Phytoplankton" & ens$lat_dir == "Poleward",c("shift","shift_lat","shift_lon")])
summary(ens[ens$group == "Phytoplankton" & ens$lat_dir == "Equatorward",c("shift","shift_lat","shift_lon")])
summary(ens[ens$group == "Zooplankton" & ens$lat_dir == "Poleward",c("shift","shift_lat","shift_lon")])
summary(ens[ens$group == "Zooplankton" & ens$lat_dir == "Equatorward",c("shift","shift_lat","shift_lon")])

median(ens[ens$group == "Phytoplankton" & ens$lat_dir == "Poleward",c("shift_lat")]); IQR(ens[ens$group == "Phytoplankton" & ens$lat_dir == "Poleward",c("shift_lat")])
median(ens[ens$group == "Zooplankton" & ens$lat_dir == "Poleward",c("shift_lat")]); IQR(ens[ens$group == "Zooplankton" & ens$lat_dir == "Poleward",c("shift_lat")])

median(ens[ens$group == "Phytoplankton" & ens$lat_dir == "Equatorward",c("shift_lat")]); IQR(ens[ens$group == "Phytoplankton" & ens$lat_dir == "Equatorward",c("shift_lat")])
median(ens[ens$group == "Zooplankton" & ens$lat_dir == "Equatorward",c("shift_lat")]); IQR(ens[ens$group == "Zooplankton" & ens$lat_dir == "Equatorward",c("shift_lat")])

median(ens[ens$group == "Phytoplankton" & ens$lon_dir == "Westward",c("shift_lat")]); IQR(ens[ens$group == "Phytoplankton" & ens$lon_dir == "Westward",c("shift_lon")])
median(ens[ens$group == "Zooplankton" & ens$lon_dir == "Westward",c("shift_lat")]); IQR(ens[ens$group == "Zooplankton" & ens$lon_dir == "Westward",c("shift_lon")])

median(ens[ens$group == "Phytoplankton" & ens$lon_dir == "Eastward",c("shift_lat")]); IQR(ens[ens$group == "Phytoplankton" & ens$lon_dir == "Eastward",c("shift_lat")])
median(ens[ens$group == "Zooplankton" & ens$lon_dir == "Eastward",c("shift_lat")]); IQR(ens[ens$group == "Zooplankton" & ens$lon_dir == "Eastward",c("shift_lat")])


# Ordinate in a PCA?
colnames(ens)
pca <- PCA(X = ens[,c(7:12)], scale.unit = TRUE, ncp = 3)
summary(pca)
ens$PC1 <- pca$ind$coord[,1]
ens$PC2 <- pca$ind$coord[,2]
ens[order(ens$PC1,decreasing = T),c("species","group","PC1","PC2")]


### Ensembles per SDM
ens.SDM <- data.frame(ddf %>% group_by(species,SDM) %>% summarize(group = unique(group),
            x_base = mean(x_base), y_base = mean(y_base), x_fut = mean(x_fut), y_fut = mean(y_fut),
            distm = mean(distm), distm_lat = mean(distm_lat), distm_lon = mean(distm_lon) )
) # eo ddf
ens.SDM$shift <- (ens.SDM$distm)/7
ens.SDM$shift_lat <- (ens.SDM$distm_lat)/7
ens.SDM$shift_lon <- (ens.SDM$distm_lon)/7

# Plot distribution by facetting per SDM
ggplot(data = ens.SDM, aes(x = factor(group), y = shift, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Shift speed (km/dec)\nbased on the absolute centroid") + theme_classic() +
    facet_wrap(~factor(ens.SDM$SDM), ncol = 2, scales = "fixed")

### Same, but per ESM
ens.ESM <- data.frame(ddf %>% group_by(species,ESM) %>% summarize(group = unique(group),
            x_base = mean(x_base), y_base = mean(y_base), x_fut = mean(x_fut), y_fut = mean(y_fut),
            distm = mean(distm), distm_lat = mean(distm_lat), distm_lon = mean(distm_lon) )
) # eo ddf
ens.ESM$shift <- (ens.ESM$distm)/7
ens.ESM$shift_lat <- (ens.ESM$distm_lat)/7
ens.ESM$shift_lon <- (ens.ESM$distm_lon)/7
# Plot distribution by facetting per ESM
ggplot(data = ens.ESM, aes(x = factor(group), y = shift, fill = factor(group))) + 
    scale_fill_manual(name = "", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) +
    xlab("") + ylab("Shift speed (km/dec)\nbased on the absolute centroid") + theme_classic() +
    facet_wrap(~factor(ens.ESM$ESM), ncol = 2, scales = "fixed")




# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
