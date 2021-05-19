
### ----------------------------------------------------------------------------------------------------------------------------

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("marmap")
library("dplyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("geosphere")
library("parallel")

### ----------------------------------------------------------------------------------------------------------------------------


##### 3°) BONUS: compute shortest distance to coast for each marine cell ! 
# https://stackoverflow.com/questions/21295302/calculating-minimum-distance-between-a-point-and-the-coast
#epsg.2062 <- "+proj=lcc +lat_1=40 +lat_0=40 +lon_0=0 +k_0=0.9988085293 +x_0=600000 +y_0=600000 +a=6378298.3 +b=6356657.142669561 +pm=madrid +units=m +no_defs"
# epsg.2062 projection system is in meters
#wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# wgs.84 projection are in degrees
### Get 10m coastline shapefile
# setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
# ?readOGR
# coast <- raster::shapefile(x = "ne_10m_coastline.shp")
# coast <- rgdal::readOGR("ne_10m_coastline", crs(wgs.84))

### Get bathy data
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 15)
# resolution = 15 mins --> 1/4° resolution
ras <- as.xyz(bathy)
colnames(ras) <- c("x","y","z")
#head(ras)
# summary(ras$z)
# Add a bolean desribing who's coastline and who's not
ras$coast <- NA
ras[ras$z <= -5,"coast"] <- FALSE # marine cells basically
ras[ras$z > -5,"coast"] <- TRUE # cells that are coastline or positive altitude

# For each cell of ras that has ras$coast == FALSE, compute Haversine distance to all other cells that have ras$coast == TRUE, and keep shortest distance
# Might need id_cell for this
ras$id_cell <- paste(ras$x, ras$y, sep = "_")
# length(unique(ras$id_cell)) # 1'036'800
# 64800 for 60min resolution
# 259200 for 30min resolution
# dim(ras) # 1'036'800

marine_ids <- unique(ras[ras$coast == FALSE, "id_cell"])
land_ids <- unique(ras[ras$coast == TRUE, "id_cell"])
# length(marine_ids) # 680'880 for 15min res
# 170139 for 30min 
# length(land_ids) # 355'920 ; 89'061
rownames(ras) <- ras$id_cell

### 02/05/2018: identify the marine_ids cells (coast == FALSE) that are comprised between 49.375°N & 76.875°N
todo <- unique(ras[which(ras$y <= 76.875 & ras$y >= 49.375 & ras$coast == FALSE),"id_cell"])
# length(unique(ras[which(ras$y <= 76.875& ras$y >= 49.375 & ras$coast == FALSE),"id_cell"])) # 68766

### To divide the work when calculating the 1/4° dist2coast raster
# 1:68088 #p1 - done
# 68089:136176 #p2 - done
# 136177:170221 #p3 - done --> redundant with p2... unique(ras[which(ras$y <= 76.875 & ras$y >= 49.375 & ras$coast == FALSE),"id_cell"])
# 170222:204264 #p4 - done
# 204265:272352 #p5 - done
# 272353:306397 #p6 - done
# 306398:340440 #p7 - done
# 340441:374485 #p8 - done
# 374485:408528 #p9 - done
# 408529:476616 #p10 - done
# 476617:510661 #p11 - done
# 510662:544705 #p12 - done
# 544706:578750 #p13 - done
# 578751:612795 #p14 - done
# 612796:646840 #p15 - done
# 646841:680880 #p16 - done

# distm(centroids[i,c("mean_lon_T1", "mean_lat_T1")], centroids[i,c("mean_lon_T2", "mean_lat_T2")], fun = distHaversine)
# ?distm
# i <- marine_ids[1] # For testing

distances <- mclapply(todo, function(i) {
					# Get id coords
					message(paste(i, spe = ""))
					xy <- ras[i,c("x","y")]
					distances <- t(distm(xy, ras[land_ids, c("x","y")], fun = distHaversine)/1000) # /1000 to get dist in km
					# summary(t(distances))
					# min(distances) # to ge the minimal distance to ANY land point
					return(data.frame(id = i, distkm = min(distances) ) )
}, mc.cores = 25 ) # eo mclapply


dist2coast <- do.call(rbind, distances)
# dim(dist2coast)

save(dist2coast, file = "dist2coast_15min_pmissing.Rdata")

#marine_ids[marine_ids == "-69.875_75.625"]
#match(c("-62.625_75.375","-4.375_75.625","128.625_75.875","-179.875_75.375","11.875_75.375"), marine_ids)







