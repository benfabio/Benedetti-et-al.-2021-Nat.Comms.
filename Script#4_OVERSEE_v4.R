
##### 03/05/2018: R Script to remove the zooplankton occurrences located in cells with salinity < 20 © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### For each large zooplankton group, aims to:
#	- associate each zooplankton occurrence to a SSS salinity based on the monthly outputs of WOA2013v2 (1/4° res)
#	- remove any occurrence that is associated to a SSS value below 20 (not oceanic)
#	- report number of lost occurrences on excel sheet

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 04/05/2018

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


##### 1°) Check how to load the SSS rasters from the nc files
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/salinity_masks_WOA13v2")
# Try with :
#ras <- raster("woa13_decav_s01_04v2.nc")
#ras
#plot(ras)
### Okay, works

### Load all monthly SSS layers
ras1 <- raster("woa13_decav_s01_04v2.nc")
ras2 <- raster("woa13_decav_s02_04v2.nc")
ras3 <- raster("woa13_decav_s03_04v2.nc")
ras4 <- raster("woa13_decav_s04_04v2.nc")
ras5 <- raster("woa13_decav_s05_04v2.nc")
ras6 <- raster("woa13_decav_s06_04v2.nc")
ras7 <- raster("woa13_decav_s07_04v2.nc")
ras8 <- raster("woa13_decav_s08_04v2.nc")
ras9 <- raster("woa13_decav_s09_04v2.nc")
ras10 <- raster("woa13_decav_s10_04v2.nc")
ras11 <- raster("woa13_decav_s11_04v2.nc")
ras12 <- raster("woa13_decav_s12_04v2.nc")


##### 2°) Apply SSS mask to Copepoda OBIS and then all other datasets
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/Copepoda_OBIS_03_05_18/")
files <- dir()
#f <- files[217] # For testing
### Apply in a for loop for 
for(f in files) {		
		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/Copepoda_OBIS_03_05_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# Cretae salinity vector
		data$salinityWOA13 <- NA
		# Get month from (and year while you're at it) eventDate and 'lubridate'
		require("lubridate")
		#head(data$eventDate)
		# str(lubridate::month(data$eventDate))
		# head(lubridate::year(data$eventDate))
		data$month <- lubridate::month(data$eventDate)
		data$year <- lubridate::year(data$eventDate)
		
		# Fill it according to month
		data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		# Check
		# summary(data$salinityWOA13)
		# str(data$salinityWOA13)
		# Remove obs with SSS <= 20
		# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
		data2 <- data[!(data$salinityWOA13 < 20),]
		# summary(data2$salinityWOA13)
		# Remove NAs 
		data3 <- data2 %>% drop_na(salinityWOA13)
		# dim(data3)
		
		# Go to v3.1 to save if you have point left
		if( nrow(data3) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/Copepoda_OBIS_03_05_18/")
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


### Check the resulting number of obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/Copepoda_OBIS_03_05_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 893191
rm(res)
gc()




##### To the same for all groups
### A) For OBIS data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
dir()
data <- get(load("Podonidae_OBIS_03_05_18.Rdata"))
dim(data)

# Get month from (and year while you're at it) eventDate and 'lubridate' FOR OBIS data only !
require("lubridate")
data$month <- lubridate::month(data$eventDate)
data$year <- lubridate::year(data$eventDate)
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# Remove obs with SSS <= 20
# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2)
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3)
summary(data3$salinityWOA13)

### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
save(data3, file = "Podonidae_OBIS_03_05_18.Rdata")



### B) For GBIF data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
dir()
data <- get(load("Hyperiidea_GBIF_03_05_18.Rdata"))
dim(data)
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimallongitude","decimallatitude")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# str(data$salinityWOA13)
# Remove obs with SSS <= 20
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2)
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3)
summary(data3$salinityWOA13)

### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
save(data3, file = "Hyperiidea_GBIF_03_05_18.Rdata")



### ----------------------------------------------------------------------------------------------------------------------------


##### 04/05/18: Apply salinity mask to v3.2 data (same as above)

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/Copepoda_OBIS_03_05_18/")
files <- dir()
#f <- files[202] # For testing
### Apply in a for loop for 
for(f in files) {		
		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/Copepoda_OBIS_03_05_18/")
		message(paste(f, sep = ""))
		data <- get(load(f))
		spname <- unique(data$species)
		# Cretae salinity vector
		data$salinityWOA13 <- NA
		# Get month from (and year while you're at it) eventDate and 'lubridate'
		require("lubridate")
		#head(data$eventDate)
		# str(lubridate::month(data$eventDate))
		# head(lubridate::year(data$eventDate))
		data$month <- lubridate::month(data$eventDate)
		data$year <- lubridate::year(data$eventDate)
		
		# Fill it according to month
		data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
		# Check
		# summary(data$salinityWOA13)
		# str(data$salinityWOA13)
		# Remove obs with SSS <= 20
		# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
		data2 <- data[!(data$salinityWOA13 < 20),]
		# summary(data2$salinityWOA13)
		# Remove NAs 
		data3 <- data2 %>% drop_na(salinityWOA13)
		# dim(data3)
		
		# Go to v3.1 to save if you have point left
		if( nrow(data3) > 0) {
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/Copepoda_OBIS_04_05_18/")
				save(data3, file = paste(spname, "_OBIS_04_05_18.Rdata", sep = ""))
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
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/Copepoda_OBIS_04_05_18/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,c("decimalLongitude","decimalLatitude")])
}) # eo lapply
# Check
d <- do.call(rbind, res)
dim(d) # 
rm(res, d)
gc()



##### To the same for all groups
### A) For OBIS data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
dir()

data <- get(load("Podonidae_OBIS_03_05_18.Rdata"))
dim(data)

# Get month from (and year while you're at it) eventDate and 'lubridate' FOR OBIS data only !
require("lubridate")
data$month <- lubridate::month(data$eventDate)
data$year <- lubridate::year(data$eventDate)
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimalLongitude","decimalLatitude")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# Remove obs with SSS <= 20
# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2)
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3)
summary(data3$salinityWOA13)

### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
save(data3, file = "Podonidae_OBIS_04_05_18.Rdata")



### B) For GBIF data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
dir()
data <- get(load("Hyperiidea_GBIF_03_05_18.Rdata"))
dim(data)
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("decimallongitude","decimallatitude")], method = 'bilinear')
data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("decimallongitude","decimallatitude")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# str(data$salinityWOA13)
# Remove obs with SSS <= 20
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2)
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3)
summary(data3$salinityWOA13)

### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
save(data3, file = "Hyperiidea_GBIF_04_05_18.Rdata")


























