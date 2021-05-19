
##### 24/04/2017: R Script to remove the zooplankton occurrences that are located near the coast © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### For each large zooplankton group, aims to:
#	- Attribute a 'z' value to each and every occurrence using a 1/4°, 1/2° or 1° bathymetric chart from NOAA ('marmap' package)
#	- Remove all occurrences that fall within the first 200m depth (too much coastal influence)
#	- Assess the importance of the resolution of the bathymetric chart to data removal (based on copepods for instance)
#	- BONUS: compute distance to coast for each cell of the ocean (i.e. cells for which z is < -200m)

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 24/04/2017

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("marmap")
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("geosphere")
library("parallel")

### ----------------------------------------------------------------------------------------------------------------------------


##### 1°) Get bathymetric chart from NOAA using the 'marmap' package
?getNOAA.bathy
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 15)
# resolution = 15 mins --> 1/4°
# resolution = 30 mins --> 1/2°
# resolution = 60 mins --> 1°
str(bathy)
# plot(bathy) # prints a plot on kryo-old...

# # For nicer plot:
# newcol <- palette.bathy(mat = bathy, layers = list(
# 	c(min(bathy),0,"purple","blue","lightblue"),
# 	c(0,max(bathy), "gray90", "gray10")), land = TRUE)
# plot(bathy, land = TRUE, n = 10, lwd = 0.5, image = TRUE, bpal = newcol)

### Convert to ddf
ras <- as.xyz(bathy)
colnames(ras) <- c("x","y","z")
# summary(ras)

### For a ggplot2 style map:
#quartz()
#ggplot() + theme_linedraw() + geom_raster(aes(x = x , y = y, fill = z), data = ras) + 
#	geom_contour(aes(x = x, y = y, z = z), data = ras, colour = "black", breaks = c(-2000,-200,-50,0)) +
#	scale_fill_distiller(name = "Bathymetry and altitude (m)", palette = "BrBG") + coord_quickmap() 

### Convert to raster
coordinates(ras) <- ~ x + y
gridded(ras) <- TRUE
# coerce to raster
r <- raster(ras)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#proj4string(r) = CRS("+init=EPSG:4326") 
r # Looks alright by me
# plot(r[r < -200])


##### 2°) Get copepod occurrences and use extract() function to attribute cooresponding z layer from 'r'
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/Copepoda_OBIS_23_04_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude")])
}) # eo lapply
# Check
d <- do.call(rbind, res)

#d <- get(load("Hydrozoa_OBIS_23_04_18.Rdata"))
dim(d) # 1'779'596 occurrences
str(d)
#unique(numeric(d$decimallongitude))
summary(d[,c("decimalLongitude","decimalLatitude")])
rm(res) ; gc()

sp <- SpatialPoints(d[,c("decimalLongitude","decimalLatitude")])
crs(sp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sp # gut gut
# ?extract
#      Extract values from a Raster* object at the locations of other
#      spatial data. You can use coordinates (points), lines, polygons or
#      an Extent (rectangle) object. You can also use cell numbers to
#      extract values.
#      If ‘y’ represents points, ‘extract’ returns the values of a
#      Raster* object for the cells in which a set of points fall. If ‘y’
#      represents lines, the ‘extract’ method returns the values of the
#      cells of a Raster* object that are touched by a line. If ‘y’
#      represents polygons, the ‘extract’ method returns the values of
#      the cells of a Raster* object that are covered by a polygon. A
#      cell is covered if its center is inside the polygon (but see the
#      ‘weights’ option for considering partly covered cells; and
#      argument ‘small’ for getting values for small polygons anyway).

# head( raster::extract(x = r, y = sp, method = 'bilinear') )
dd <- raster::extract(x = r, y = sp, method = 'bilinear')
length(dd) # 1779596 values ok.
summary(dd)
# Provide to 'd'
d$z <- dd

### How many point will you loose then...?
nrow(d[d$z > -200,]) / nrow(d) # 0.475 for res15 ; 0.464 for res30 ; 0.439 for res 60 ; 0.47 for res5
nrow(d[d$z > -100,]) / nrow(d) # 0.344 ; 0.3412 ; 0.318 ; 0.347
nrow(d[d$z > -50,]) / nrow(d) # 0.20 ; 0.226 ; 0.2203 ; 0.203
nrow(d[d$z > -20,]) / nrow(d) # 0.123 ; 0.123 ; 0.144 ; 
nrow(d[d$z > 0,]) / nrow(d) # 0.085 --> points on land ; 0.092 ; 0.0949 ; 0.062

### Wouldn't you rather just use a salinity mask :-/ ? Meike says no but feel free to test it


##### 03/05/18: Apply the 200m depth threshold anyways...
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/Copepoda_OBIS_23_04_18/")
files <- dir()
# f <- files[10] # For testing
### Apply in a for loop for 
for(f in files) {		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/Copepoda_OBIS_23_04_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# summary(data[,c("decimalLongitude","decimalLatitude")])
		# summary(raster::extract(x = r, y = data[,c("decimalLongitude","decimalLatitude")], method = 'bilinear'))
		# summary(raster::extract(x = r, y = SpatialPoints(data[,c("decimalLongitude","decimalLatitude")]), method = 'bilinear'))
		# colnames(data)
		data$bathymetry <- raster::extract(x = r, y = data[,c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data2 <- data[-which(data$bathymetry >= -200),]
		# Go to v3.1 to save if you have point left
		if( nrow(data2) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/Copepoda_OBIS_03_05_18/")
				save(data2, file = paste(spname, "_OBIS_03_05_18.Rdata", sep = ""))
				# Clean
				rm(data, data2, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, spname)
				gc()
		}
} # eo for loop

### Count remaining obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/Copepoda_OBIS_03_05_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d)


### To the same for all 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
dir()
data <- get(load("Cubozoa_GBIF_23_04_18.Rdata"))
dim(data)
#summary( as.numeric(levels(data$decimallongitude))[data$decimallongitude] )
#data$decimallongitude <- as.numeric(levels(data$decimallongitude))[data$decimallongitude] 
#data$decimallatitude <- as.numeric(levels(data$decimallatitude))[data$decimallatitude] 
#summary(data[,c("decimallongitude","decimallatitude")])
#data <- data %>% drop_na(decimallongitude, decimallatitude)

data$bathymetry <- raster::extract(x = r, y = data[,c("decimallongitude","decimallatitude")], method = 'bilinear')
data2 <- data[-which(data$bathymetry >= -200),]
dim(data2)
summary(data2) # check NAs
data3 <- data2 %>% drop_na(bathymetry)
dim(data3)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
save(data3, file = "Cubozoa_GBIF_03_05_18.Rdata")





##### 3°) BONUS: compute shortest distance to coast for each marine cell ! 
# https://stackoverflow.com/questions/21295302/calculating-minimum-distance-between-a-point-and-the-coast
epsg.2062 <- "+proj=lcc +lat_1=40 +lat_0=40 +lon_0=0 +k_0=0.9988085293 +x_0=600000 +y_0=600000 +a=6378298.3 +b=6356657.142669561 +pm=madrid +units=m +no_defs"
# epsg.2062 projection system is in meters
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
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
head(ras)
# summary(ras$z)
# Add a bolean desribing who's coastline and who's not
ras$coast <- NA
ras[ras$z <= -5,"coast"] <- FALSE # marince cells basically
ras[ras$z > -5,"coast"] <- TRUE # cells that are coastline or positive altitude

# For each cell of ras that has ras$coast == FALSE, compute Haversine distance to all other cells that have ras$coast == TRUE, and keep shortest distance
# Might need id_cell for this
ras$id_cell <- paste(ras$x, ras$y, sep = "_")
length(unique(ras$id_cell)) # 1'036'800
dim(ras) # 1'036'800

marine_ids <- unique(ras[ras$coast == FALSE,"id_cell"])
land_ids <- unique(ras[ras$coast == TRUE,"id_cell"])
length(marine_ids) # 680'880
length(land_ids) # 355'920
rownames(ras) <- ras$id_cell

# distm(centroids[i,c("mean_lon_T1", "mean_lat_T1")], centroids[i,c("mean_lon_T2", "mean_lat_T2")], fun = distHaversine)
# ?distm
# i <- marine_ids[1] # For testing
system.time(
distances <- mclapply(marine_ids[1:1000], function(i) {
					# Get id coords
					message(paste(i, spe = ""))
					xy <- ras[i,c("x","y")]
					distances <- t(distm(xy, ras[land_ids, c("x","y")], fun = distHaversine)/1000) # /1000 to get dist in km
					# summary(t(distances))
					# min(distances) # to ge the minimal distance to ANY land point
					return(data.frame(id = i, distkm = min(distances) ) )
}, mc.cores = 22
) # eo mclapply on a 20 cores machine
)

dist2coast <- do.call(rbind, distances)
dim(dist2coast)

# ### Convert to raster
# coordinates(ras) <- ~ x + y
# gridded(ras) <- TRUE
# # coerce to raster
# r <- raster(ras)
# crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# #proj4string(r) = CRS("+init=EPSG:4326")
# r





### ----------------------------------------------------------------------------------------------------------------------------

##### 26/04/18: Make a raster out of the 'dist2coast' object
# Get world coast
# setwd("/Users/fabiobenedetti/Desktop/PostDocs/TARA/")
setwd("/Users/fabiobenedetti/Desktop/")
dist <- get(load("dist2coast_60min.Rdata"))
dim(dist)
head(dist)
summary(dist)

# ### 02/05/18: For the 1/4° distance to coast raster
# distances <- dir()[grep("dist2coast_15min", dir())]
# res <- lapply(distances, function(f) {
# 			message(paste(f, sep = ""))
# 			data <- get(load(f))
# 			data$pice <- f
# 			return(data)
# }) # eo lapply
# # Check
# dist <- do.call(rbind, res)
# dim(dist)
# head(dist)
# tail(dist)
# summary(dist)
# count(dist$pice)
# rm(res) ; gc()
# Provide coordinates from ids (x_y)
dim(do.call(rbind, strsplit(as.character(dist$id), "_")))
head( do.call(rbind, strsplit(as.character(dist$id), "_")) )
dist$x <- as.numeric(do.call(rbind, strsplit(as.character(dist$id), "_"))[,1])
dist$y <- as.numeric(do.call(rbind, strsplit(as.character(dist$id), "_"))[,2])
summary(dist)

# library("dplyr")
# ddf <- data.frame(dist %>% group_by(pice) %>% summarise(count = n() ))
# sum(ddf$count)
# for(p in unique(dist$pice)) {
# 		message(paste(min(dist[dist$pice == p,"y"]), max(dist[dist$pice == p,"y"]), sep = "   "))
# 		message(paste(min(dist[dist$pice == p,"y"]), max(dist[dist$pice == p,"y"]), sep = "   "))
# }

# Map
quartz()
ggplot() + theme_bw() + geom_raster(aes(x = x, y = y, fill = log(distkm)), data = dist) + 
	scale_fill_distiller(name = "Distance to the\nnearest coast log(km)", palette = "Blues", direction = 1) + 
	coast + coord_quickmap()
# Nice

# Plot the distrbution of the distance data (log and not log)
quartz()
ggplot(dist, aes(distkm)) + geom_density(fill = "gray") + xlab("Distance to coast (km)") + ylab("Density") + theme_linedraw()
quartz()
ggplot(dist, aes(log(distkm))) + geom_density(fill = "gray") + xlab("Distance to coast log(km)") + ylab("Density") + theme_linedraw()

### Convert 'dist' to raster
library("raster")
colnames(dist)[2] <- c("z")
# Put x, y and z in the right order
ddf <- dist[,c("x","y","z")]
ras <- rasterFromXYZ(ddf)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ras # looks good !
# plot(ras)

save(ras, file = "dist2coast_raster_60min.Rdata")





### ----------------------------------------------------------------------------------------------------------------------------


##### 03/05/18: Apply v3.2 on v2 -> Distance to coast (10kms/20kms)
### Get distance to coast raster, 15min (1/4°)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/dist2coast")
r <- get(load("dist2coast_raster_15min.Rdata"))
r #

### Quickly examine loss of information with Copepoda OBIS files
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/Copepoda_OBIS_23_04_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
#unique(numeric(d$decimallongitude))
summary(d[,c("decimalLongitude","decimalLatitude")])
rm(res) ; gc()
sp <- SpatialPoints(d[,c("decimalLongitude","decimalLatitude")])
crs(sp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# head( raster::extract(x = r, y = sp, method = 'bilinear') )
dd <- raster::extract(x = r, y = sp, method = 'bilinear')
length(dd) # 1779596 values ok.
summary(dd)
# Provide to 'd'
d$z <- dd

### How many point will you loose then...?
### IMPORTANT: 2 1/4° cells are separated by nearly 28km already

nrow(d[d$z < 5,]) / nrow(d) 	# 0.0150
nrow(d[d$z < 10,]) / nrow(d) 	# 0.0157
nrow(d[d$z < 20,]) / nrow(d) 	# 0.0835
nrow(d[d$z < 25,]) / nrow(d) 	# 0.1439
nrow(d[d$z < 50,]) / nrow(d) 	# 0.3457
nrow(d[is.na(d$z),]) / nrow(d) 	# 0.0149

# Let's go with 25kms (and don't forget to remove NAs as well)


##### 03/05/18: Apply the 25kms
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/Copepoda_OBIS_23_04_18/")
files <- dir()
# f <- files[10] # For testing
### Apply in a for loop for 
for(f in files) {		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/Copepoda_OBIS_23_04_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# summary(data[,c("decimalLongitude","decimalLatitude")])
		# summary(raster::extract(x = r, y = data[,c("decimalLongitude","decimalLatitude")], method = 'bilinear'))
		# summary(raster::extract(x = r, y = SpatialPoints(data[,c("decimalLongitude","decimalLatitude")]), method = 'bilinear'))
		# colnames(data)
		data$dist2coast <- raster::extract(x = r, y = data[,c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data2 <- data[-which(data$dist2coast < 25),]
		summary(data2$dist2coast)
		# Remove NAs --> land cells
		data3 <- data2 %>% drop_na(dist2coast)
		dim(data3)
		
		# Go to v3.1 to save if you have point left
		if( nrow(data3) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/Copepoda_OBIS_03_05_18/")
				save(data3, file = paste(spname, "_OBIS_03_05_18.Rdata", sep = ""))
				# Clean
				rm(data, data2, data3, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, data3, spname)
				gc()
		}
} # eo for loop

### Count remaining obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/Copepoda_OBIS_03_05_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d)


### To the same for all 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
dir()
data <- get(load("Thaliacea_OBIS_23_04_18.Rdata"))
dim(data)

#summary( as.numeric(levels(data$decimallongitude))[data$decimallongitude] )
#data$decimallongitude <- as.numeric(levels(data$decimallongitude))[data$decimallongitude] 
#data$decimallatitude <- as.numeric(levels(data$decimallatitude))[data$decimallatitude] 
#summary(data[,c("decimallongitude","decimallatitude")])
#data <- data %>% drop_na(decimallongitude, decimallatitude)

data$dist2coast <- raster::extract(x = r, y = data[,c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data2 <- data[-which(data$dist2coast < 25),]
summary(data2$dist2coast)
# Remove NAs --> land cells
data3 <- data2 %>% drop_na(dist2coast)
dim(data3)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
save(data3, file = "Thaliacea_OBIS_03_05_18.Rdata")


























