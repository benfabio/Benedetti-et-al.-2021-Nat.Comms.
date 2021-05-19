
##### 04/05/2018: R Script to remove the zooplankton occurrences coming from net tows going deeper than 500m and/or those who exceed MLD © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### For each large zooplankton group, aims to:
#	- remove any occurrence coming from a net tow that went deeper than 500m
#	- load MLD monthly climatologies (De Boyer Montégut - 2°x2°) and apply MLD criterion
#	- report number of lost occurrences on excel sheet at each step

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 15/05/2018

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")

### ----------------------------------------------------------------------------------------------------------------------------

##### 1°) 500m depth criterion (v5.1v3.1 & v5.1v3.2)

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/Copepoda_OBIS_03_05_18/")
files <- dir()
# f <- files[217] # For testing

### Apply in a for loop for 
for(f in files) {		
		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/Copepoda_OBIS_03_05_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# summary(data$depth)
		# Sometimes you still have negative depths for some reason... -> absolute
		data$depth <- abs(data$depth)
		# Remove obs with depth > 500
		data2 <- data[!(data$depth > 500),]
		# summary(data2$depth)
		# dim(data2)
		
		# Go to v3.1 to save if you have point left
		if( nrow(data2) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/Copepoda_OBIS_04_05_18/")
				save(data2, file = paste(spname, "_OBIS_04_05_18.Rdata", sep = ""))
				# Clean
				rm(data, data2, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, spname)
				gc()
		}
		
} # eo for loop


### Check the resulting number of obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/Copepoda_OBIS_04_05_18/")
files <- dir() ; files
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude","depth")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 863'292
summary(d)
 
rm(res,d)
gc()


### Do the same for the other datasets
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
dir()
data <- get(load("Copepoda_GBIF_03_05_18.Rdata"))
dim(data)
colnames(data)
summary(data$depth)

### Since depth has the same colname in both OBIS and GBIF files, apply depth criterion in a for loop
files <- dir()[c(1:4,6:34)]

for(f in files) {
	
		# Get/ load/ message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
		message(paste("Doing ",f, sep = ""))
		data <- get(load(f))
		
		# summary(data$depth)
		# Sometimes you still have negative depths for some reason... -> absolute
		data$depth <- abs(data$depth)
		# Remove obs with depth > 500
		data2 <- data[!(data$depth > 500),]
		# summary(data2$depth)
		# dim(data2)
		message(paste("DIM = ", nrow(data2), sep = ""))
		
		# Go to v5.1v3.1 to save if you have point left
		if( nrow(data2) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
				require("stringr")
				save(data2, file = paste(str_replace_all(string = f, pattern = "_03_05", replacement = "_04_05"), sep = ""))
				# Clean
				rm(data, data2, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, spname)
				gc()
		}
	
} # eo for loop

### OKay, check some results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
dir()
files <- dir()[c(1:4,6:34)]
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data.frame(data[,c("depth")]))
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 
summary(d) # Ok gut



### ----------------------------------------------------------------------------------------------------------------------------

##### Do the same as above but for the v3.2 (dist2coast) data

### First, Copepoda OBIS
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/Copepoda_OBIS_04_05_18/")
files <- dir()
# f <- files[217] # For testing

### Apply in a for loop for 
for(f in files) {		
		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/Copepoda_OBIS_04_05_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# summary(data$depth)
		# Sometimes you still have negative depths for some reason... -> absolute
		data$depth <- abs(data$depth)
		# Remove obs with depth > 500
		data2 <- data[!(data$depth > 500),]
		# summary(data2$depth)
		# dim(data2)
		
		# Go to v3.1 to save if you have point left
		if( nrow(data2) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/Copepoda_OBIS_04_05_18/")
				save(data2, file = paste(spname, "_OBIS_04_05_18.Rdata", sep = ""))
				# Clean
				rm(data, data2, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, spname)
				gc()
		}
		
} # eo for loop


### Check the resulting number of obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/Copepoda_OBIS_04_05_18/")
files <- dir() ; files
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude","depth")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 1344909
summary(d) # ok gut gut
 
rm(res,d)
gc()



### Do the same for the other datasets
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
dir()

### Since depth has the same colname in both OBIS and GBIF files, apply depth criterion in a for loop
files <- dir()[c(1:4,6:34)]

for(f in files) {
	
		# Get/ load/ message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
		message(paste("Doing ",f, sep = ""))
		data <- get(load(f))
		
		# summary(data$depth)
		# Sometimes you still have negative depths for some reason... -> absolute
		data$depth <- abs(data$depth)
		# Remove obs with depth > 500
		data2 <- data[!(data$depth > 500),]
		# summary(data2$depth)
		# dim(data2)
		message(paste("DIM = ", nrow(data2), sep = ""))
		
		# Go to v5.1v3.1 to save if you have point left
		if( nrow(data2) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
				require("stringr")
				save(data2, file = paste(str_replace_all(string = f, pattern = "_03_05", replacement = "_04_05"), sep = ""))
				# Clean
				rm(data, data2, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, spname)
				gc()
		}
	
} # eo for loop

### OKay, check some results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
dir()
files <- dir()[c(1:4,6:34)]
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data.frame(data[,c("depth")]))
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 
summary(d) # Ok gut




### ----------------------------------------------------------------------------------------------------------------------------


##### 2°) MLD criterion (v5.2v3.1 & v3.2)

### Read netcdf file containing the MLD monthly climatologies and turn them into rasters
nc <- nc_open("mld_DT02_c1m_reg2.0.nc")
nc
names(nc$var)
# mld_da_mean = average mld based on density algorithm
# mld_dt_mean = average mld based on temperature algorithm
# print(nc)

# Extract data
lons <- ncvar_get(nc, "lon")
dim(lons) # 180 
lons # Need to - 180
lons2 <- lons - 180
lons2
lats <- ncvar_get(nc, "lat")
dim(lats) # 90
head(lats)
lats # This is OK

mld <- ncvar_get(nc, "mld")
dim(mld) # 180  90  12 -> lon/ lat / months
class(mld) # 3D array
length(mld) # 12*180*90
dimnames(mld)
dimnames(mld)[[1]] <- as.character(lons2)
dimnames(mld)[[2]] <- as.character(lats)
dimnames(mld)[[3]] <- as.character(c(1:12))
# Melt it and rename columns
mlds <- melt(mld)
dim(mlds)
head(mlds)
colnames(mlds)[1:4] <- c("x","y","Month","MLD")
summary(mlds) # 
# missing_value: -9999
# min_value: 10.2676773071289
# max_value: 772.370483398438
# mask_value: 1e+09

# Change mask value to NA
mlds$MLD[mlds$MLD == 1e+09] <- NA
summary(mlds) # Gut gut
str(mlds)
mlds$Month <- as.factor(mlds$Month)

### Convert a monthly climatology to raster for a try
ras6 <- mlds[which(mlds$Month == 12),c("x","y","MLD")]
colnames(ras6)[3] <- c("z")
coordinates(ras6) <- ~ x + y
# create an empty raster object to the extent of the points
rast <- raster(ext = extent(ras6), resolution = 2)
rasOut <- rasterize(ras6, rast, ras6$z, fun = mean)
rasOut
crs(rasOut) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(rasOut) # Seemed to have worked out

### Now, in a for loop, perform this and save 12 monthly reasters of MLD
### NOTE: to be ran on personal computer of you want to see the plots appearing...
for(i in unique(mlds$Month)) {
	
		# Message
		message(paste("Doing MLD raster for month ", i, sep = ""))
		ras <- mlds[which(mlds$Month == i),c("x","y","MLD")]
		colnames(ras)[3] <- c("z")
		coordinates(ras) <- ~ x + y
		# create an empty raster object to the extent of the points
		rast <- raster(ext = extent(ras), resolution = 2)
		rasOut <- rasterize(ras, rast, ras$z, fun = mean)
		#rasOut
		crs(rasOut) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
		quartz()
		plot(rasOut) # Seemed to have worked out
		# Save
		save(rasOut, file = paste("raster_mld_2d_",i,"_05_05_18.Rdata", sep = ""))
		
} # eo for loop

### OK


### ----------------------------------------------------------------------------------------------------------------------------


##### 07/05/18: Apply MLD crieria to v4v3.1 and v4v3.2 data --> v5.2v3.1 & v5.2v3.2

### First, get load the MLD 1° monthly climatologies from DR
stack <- get(load("deBoyerMontegut_MLD_RasterStack_V1.RData"))
class(stack)
stack
ras1 <- subset(stack, 1)
ras2 <- subset(stack, 2)
ras3 <- subset(stack, 3)
ras4 <- subset(stack, 4)
ras5 <- subset(stack, 5)
ras6 <- subset(stack, 6)
ras7 <- subset(stack, 7)
ras8 <- subset(stack, 8)
ras9 <- subset(stack, 9)
ras10 <- subset(stack, 10)
ras11 <- subset(stack, 11)
ras12 <- subset(stack, 12)


setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/Copepoda_OBIS_04_05_18/")
files <- dir()
# f <- files[558] # For testing
### Apply in a for loop for 
for(f in files) {		
		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/Copepoda_OBIS_04_05_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# Get month from (and year while you're at it) eventDate and 'lubridate'
		require("lubridate")
		data$month <- lubridate::month(data$eventDate)
		data$year <- lubridate::year(data$eventDate)
		data <- data %>% drop_na(depth) # Remove depth NAs (if some are left for SOME STRANGE reason)
		data$depth <- abs(data$depth)
		# Create salinity vector
		data$MLD <- NA
		# Fill it according to month
		### NOTE: IMPORTANT to specify the extract fun from the 'raster' package
		data[which(data$month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		# Check
		# summary(data)
		# Remove obs based on mld
		data2 <- data[-which(data$depth > (data$MLD+50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
		# summary(data2$MLD)
		data3 <- data2[!(data2$depth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
		#data3 <- data2 %>% drop_na(MLD) # Remove NAs
		# dim(data3)
		# summary(data3) ### Strange, sometimes NAs are created...use drop_na() to remove those
		data3 <- data3 %>% drop_na(year)
		
		# Go to v3.1 to save if you have point left
		if( nrow(data3) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/Copepoda_OBIS_15_05_18/")
				save(data3, file = paste(spname, "_OBIS_15_05_18.Rdata", sep = ""))
				# Clean
				rm(data, data2, data3, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, data3, spname)
				gc()
		}
		
} # eo for loop



### Check the resulting number of obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/Copepoda_OBIS_15_05_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude","species","depth","MLD","salinityWOA13")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 622206 ; 1072621
summary(d)

rm(res, d)
gc()





##### And do the same for the other datasets
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
dir()
### Since depth has the same colname in both OBIS and GBIF files, apply depth criterion in a for loop
files <- dir()[c(1:4,6:34)] # Remove Copepoda_OBIS which you just did
# f <- "Sagittoidea_OBIS_03_05_18.Rdata"

for(f in files) {
	
		# Get/ load/ message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
		message(paste("Doing ",f, sep = ""))
		data <- get(load(f))
		# Sometimes you still have negative depths for some reason... -> absolute
		data$depth <- abs(data$depth)
		data <- data %>% drop_na(depth) # Remove depth NAs
		
		# If OBIS type of file -> re-provide year and month from eventDate and use decimalLongitude and decimalLatitude
		if( grepl("OBIS", f) ) {
			
				require("lubridate")
				data$month <- lubridate::month(data$eventDate)
				data$year <- lubridate::year(data$eventDate)
		
				### Fill it according to month
				data$MLD <- NA 
				# NOTE: IMPORTANT to specify the extract fun from the 'raster' package
				# require("raster")
				data[which(data$month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
				data[which(data$month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		
		} else {
			
				### Fill it according to month
				data$MLD <- NA 
				# require("raster")
				data[which(data$month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimallongitude","decimallatitude")], method = 'bilinear')
				data[which(data$month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimallongitude","decimallatitude")], method = 'bilinear')
			
		}
		
		# Check
		# summary(data$MLD)
		# summary(data$depth)
		# Remove obs with SSS <= 20
		data2 <- data[-which(data$depth > (data$MLD+50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
		# summary(data2$MLD)
		data3 <- data2[!(data$depth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
		# dim(data3)
		data3 <- data3 %>% drop_na(year)
		message(paste("DIM = ", nrow(data3), sep = ""))
		
		# Go to v5.1v3.1 to save if you have point left
		if( nrow(data3) > 0 ) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
				require("stringr")
				save(data3, file = paste(str_replace_all(string = f, pattern = "_04_05", replacement = "_15_05"), sep = ""))
				# Clean
				rm(data, data2, data3, spname)
				gc()
		} else {
				# Clean
				rm(data, data2, data3, spname)
				gc()
		}
	
} # eo for loop


### Check results
