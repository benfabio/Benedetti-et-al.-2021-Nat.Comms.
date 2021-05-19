
##### 01/06/2018: R Script to 

### Aims to:
#	- Load the 2 datasets sent by Astrid Cornils (Southern Ocena copepods and others)
#	- Examine how datasets are structured
#	- Upload them on kryo
#	- Clean it by applying the same criteria you applied for the GBIF & the OBIS datasets :
#		• missing coordinates
#		• missing date
#		• missing depth
#		• date < 1800 (shouldn't be the case though)
#		• not species level
#		• land mask and dist2coast mask
#		• salinity mask
#		• sampling depth or MLD mask
#		• correct species names

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 04/06/2018

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

##### 1°) Load a .tab file and examine its structure, then check if they all the same
data <- read.table("antarctic_out_AC.txt", sep = "\t", h = TRUE, fill = NA)
dim(data)
str(data)
colnames(data)

# Before melting the dataset, replace double dots by "_"
colnames(data)[11:490] <- gsub("...", "_", as.character(colnames(data)[11:490]), fixed = TRUE)
colnames(data)[11:490] <- gsub("..", "_", as.character(colnames(data)[11:490]), fixed = TRUE)
colnames(data)[11:490] <- gsub(".", "_", as.character(colnames(data)[11:490]), fixed = TRUE)
# Check

# OK, they don't have the same labels when it comes to species names (logic), so you could treat them separately 


### A) Merge both datasets after melting them and giving them the same colnames

## 1) Southern Ocean data
data <- read.table("antarctic_out_AC.txt", sep = "\t", h = TRUE, fill = NA)
dim(data)
str(data)
colnames(data)

# Before melting the dataset, replace double dots by "_"
colnames(data)[11:490] <- gsub("...", "_", as.character(colnames(data)[11:490]), fixed = TRUE)
colnames(data)[11:490] <- gsub("..", "_", as.character(colnames(data)[11:490]), fixed = TRUE)
colnames(data)[11:490] <- gsub(".", "_", as.character(colnames(data)[11:490]), fixed = TRUE)

colnames(data)
### Remove any colnames that has sp or spp in it
# colnames(data)[grep("_sp_",colnames(data))]
# colnames(data)[grep("_spp",colnames(data))]
data <- data[,- which(names(data) %in% c(colnames(data)[grep("_spp",colnames(data))],colnames(data)[grep("_sp_",colnames(data))]))]
dim(data)

# Remove some other columns
data <- data[,- which(names(data) %in% 
		c("Copepoda_nauplii","Calanidae_copepodites_not_identified_","Harpacticoida_not_identified_","Monstrilloida","Mormonilloida",
		"Date_Time_of_event","Lubbockiidae","Harpacticoida","Mormonillidae","Oncaeidae","Oithonidae",""))]

# Now, melt to put species as rows and the abund + metadata as variable
d1 <- melt(data, id.vars = colnames(data)[c(1:10,397)])
d1[1:100,]
colnames(d1)
colnames(d1)[13] <- "abundance"
colnames(d1)[12] <- "species"
### Looks nice, extract day, month and year from date
library("lubridate")
d1$day <- lubridate::day(d1$date)
d1$month <- lubridate::month(d1$date)
d1$year <- lubridate::year(d1$date)
summary(d1)
str(d1)
# Remove NA abund
d1 <- d1 %>% drop_na(abundance)
summary(d1)
# And remove abund == zero
d1 <- d1[which(d1$abundance > 0),]
dim(d1) # 50'378

colnames(d1) # Use these colnames as reference for the other datasets, before rbinding



## 2) Non Southern Ocean
data <- read.table("others_out_AC.txt", sep = "\t", h = TRUE, fill = NA)
dim(data) #  464 705
str(data)
colnames(data)

# Before melting the dataset, replace double dots by "_"
colnames(data)[11:704] <- gsub("...", "_", as.character(colnames(data)[11:704]), fixed = TRUE)
colnames(data)[11:704] <- gsub("..", "_", as.character(colnames(data)[11:704]), fixed = TRUE)
colnames(data)[11:704] <- gsub(".", "_", as.character(colnames(data)[11:704]), fixed = TRUE)
colnames(data)[1:10] <- colnames(d1)[1:10]
colnames(data)
# Remove any colnames that has sp or spp in it
data <- data[,- which(names(data) %in% c(colnames(data)[grep("_spp",colnames(data))],colnames(data)[grep("_sp_",colnames(data))]))]
dim(data)

# Remove some other columns
data <- data[,- which(names(data) %in% 
		c("Copepoda_nauplii","Calanidae_copepodites_not_identified_","Harpacticoida","Monstrilloida","Mormonillidae",
		"Date_Time","Lubbockiidae","Harpacticoida","Mormonillidae","Oncaeidae","Oithonidae","Calanidae_c1_c3_not_identified_","Scolecitrichidae",
		"Calanidae_not_identified_","Optional_event_label","Calanoida_indeterminata_male","","",""))]

# OK. Now, melt to put species as rows and the abund + metadata as variable
d2 <- melt(data, id.vars = colnames(data)[c(1:10,552)])
d2[1:100,] # Ok looks gut
colnames(d2)
colnames(d2)[13] <- "abundance"
colnames(d2)[12] <- "species"
### Looks nice, extract day, month and year from date
library("lubridate")
d2$day <- lubridate::day(d2$date)
d2$month <- lubridate::month(d2$date)
d2$year <- lubridate::year(d2$date)
summary(d2)
str(d2)
# Remove NA abund
d2 <- d2 %>% drop_na(abundance)
summary(d2)
# And remove abund == zero
d2 <- d2[which(d2$abundance > 0),]
dim(d2) # 18'961
# Check if colnames match before rbinding
colnames(d2) ; colnames(d1)

ddf <- rbind(d1, d2)
dim(ddf) # 69'339
unique(ddf$species)
### OK, export species names and correct them using the 

save(ddf, file = "Copepoda_AC_01_06_18.Rdata")
write.csv(unique(ddf$species), "Coeppoda_species_AC_01_06_18.csv")



### B) Load corrected species names table and use to correct the Southern Ocean copepod dataset
dir()
names <- read.csv("Copepoda_species_AC_04_06_18.csv", h = TRUE, sep = ";")
dim(names)
str(names)
head(names)

data <- get(load("Copepoda_AC_01_06_18_toclean.Rdata"))
dim(data)
head(data)

# 1) Remove the species that are marked as 'remove'
toremove <- unique(names[which(names$action == "remove"),"correct_name"])
# Only 5 levels

data2 <- data[-which(data$species %in% toremove),]
dim(data) ; dim(data2) # Ok, not that many obs anyways

# 2) Correct species names...what would be the simplest strategy..for loop?
tocorrect <- unique(names[which(names$action == "correct"),"correct_name"])
### About issues when changing factor levels...
# https://stackoverflow.com/questions/11810605/replace-contents-of-factor-column-in-r-dataframe
library("dplyr")
library("forcats")


# For testing: 
# sp <- tocorrect[3]
for(sp in tocorrect) {
		# Useless message
		message(paste(sp, sep = ""))
		# Find the wrong names that correspond to 'sp', the real name
		wrongnames <- names[names$correct_name == sp,"species"]
		# wrongnames
		# Correct
		data2 <- as.matrix(data2) # dim(data2)
		data2[which(data2[,"species"] %in% wrongnames),"species"] <- as.character(sp)
		data2 <- as.data.frame(data2)
		# str(data2)
		# unique(data2$species)
} # eo for loop

# Check results
unique(data2$species)
### Great, seems to have worked ! Now you can focusq on cleaning the data, save & upload on kryo !
save(data2, file = "Copepoda_AC_04_06_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------


##### 2°) Load the "Copepoda_AC_04_06_18.Rdata" and apply all cleaning filters !

### 1) v2
#		• missing coordinates
#		• missing date
#		• missing depth
#		• date < 1800 (shouldn't be the case though)

d <- get(load("Copepoda_AC_04_06_18.Rdata"))
dim(d) # 69328   16
str(d)
d$abundance <- as.numeric(as.character(d$abundance)) # convert factor to numeric through characters
require("lubridate")
d$date <- as.Date(d$date) 
d$day <- lubridate::day(d$date)
d$month <- lubridate::month(d$date)
d$year <- lubridate::year(d$date)

# Check dates
summary(d[,c("day","month","year")]) # ok

# Check coordinates, depth etc
d$y <- as.numeric(as.character(d$y))
d$x <- as.numeric(as.character(d$x))
d$elevation <- as.numeric(as.character(d$elevation))
d$depth <- as.numeric(as.character(d$depth))
d$mindepth <- as.numeric(as.character(d$mindepth))
d$maxdepth <- as.numeric(as.character(d$maxdepth))
d$volume <- as.numeric(as.character(d$volume))
str(d)
summary(d) # ok

### OK, save to v2
dim(d) 
length(unique(d$species)) # 288 species
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ---------------------------------------------------------

### 2) v3.1: apply land mask
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
d <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
str(d)

### Get the bathymetry raster
require("marmap")
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 15)
str(bathy)
# Convert to ddf
ras <- as.xyz(bathy)
colnames(ras) <- c("x","y","z")
# summary(ras)
# Convert to raster
coordinates(ras) <- ~ x + y # takes a few seconds
gridded(ras) <- TRUE
# coerce to raster
r <- raster(ras)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r # Looks alright by me

### extract() the bathymetry data and remove what is shallower than 200m
d$bathymetry <- raster::extract(x = r, y = d[,c("x","y")], method = 'bilinear')
d2 <- d[-which(d$bathymetry >= -200),]
dim(d2) # 62893
# dim(d)
summary(d2) # check NAs
d3 <- d2 %>% drop_na(bathymetry)
dim(d3) # 62893
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
save(d3, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### 3) v3.2: apply dist2coast
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/dist2coast")
r <- get(load("dist2coast_raster_15min.Rdata"))
r #
### extract() the dist2coast data and remove what is within the first 25 kms
d$dist2coast <- raster::extract(x = r, y = d[,c("x","y")], method = 'bilinear')
d2 <- d[-which(d$dist2coast < 25),]
summary(d2)
# Remove NAs --> land cells
d3 <- d2 %>% drop_na(dist2coast)
dim(d3) # 55151
summary(d3)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
save(d3, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ---------------------------------------------------------

### 3) v4v3.1 & v4v3.2: Apply salinity mask
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/salinity_masks_WOA13v2")
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


### A) Get v3.1 data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("x","y")], method = 'bilinear')
data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("x","y")], method = 'bilinear')
data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("x","y")], method = 'bilinear')
data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("x","y")], method = 'bilinear')
data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("x","y")], method = 'bilinear')
data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("x","y")], method = 'bilinear')
data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("x","y")], method = 'bilinear')
data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("x","y")], method = 'bilinear')
data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("x","y")], method = 'bilinear')
data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("x","y")], method = 'bilinear')
data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("x","y")], method = 'bilinear')
data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# Remove obs with SSS <= 20
# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2) # 62893
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3) # 61442
summary(data3$salinityWOA13)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
save(data3, file = "Copepoda_PANGAEA_04_06_18.Rdata")



### B) Get v3.2 data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("x","y")], method = 'bilinear')
data[which(data$month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("x","y")], method = 'bilinear')
data[which(data$month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("x","y")], method = 'bilinear')
data[which(data$month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("x","y")], method = 'bilinear')
data[which(data$month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("x","y")], method = 'bilinear')
data[which(data$month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("x","y")], method = 'bilinear')
data[which(data$month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("x","y")], method = 'bilinear')
data[which(data$month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("x","y")], method = 'bilinear')
data[which(data$month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("x","y")], method = 'bilinear')
data[which(data$month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("x","y")], method = 'bilinear')
data[which(data$month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("x","y")], method = 'bilinear')
data[which(data$month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# Remove obs with SSS <= 20
# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2) # 55151
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3) # 55151
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
save(data3, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ---------------------------------------------------------

### 4) v5.1v3.1 & v5.1v3.2: Apply depth mask !

### A) Remove any obs that has a sampling depth below 500m (v5.1)
# 1) v4v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
summary(data$maxdepth)
# Sometimes you still have negative depths for some reason... -> absolute
# Remove obs with depth > 500
data2 <- data[!(data$maxdepth > 500),]
summary(data2$maxdepth)
dim(data2) # 46035
# Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
save(data2, file = "Copepoda_PANGAEA_04_06_18.Rdata")

# 2) v4v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
summary(data$maxdepth)
# Remove obs with depth > 500
data2 <- data[!(data$maxdepth > 500),]
summary(data2$maxdepth) 
dim(data2) # 41019
# Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
save(data2, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ---------------------------------------------------------


### 5) v5.2v3.1 & v5.2v3.2: Apply MLD mask !
# First, load the MLD climatologies
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/")
stack <- get(load("deBoyerMontegut_MLD_RasterStack_V1.RData"))
# class(stack)
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


# 1) v4v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
# Create salinity vector
data$MLD <- NA
# Fill it according to month
data[which(data$month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("x","y")], method = 'bilinear')
data[which(data$month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("x","y")], method = 'bilinear')
data[which(data$month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("x","y")], method = 'bilinear')
data[which(data$month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("x","y")], method = 'bilinear')
data[which(data$month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("x","y")], method = 'bilinear')
data[which(data$month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("x","y")], method = 'bilinear')
data[which(data$month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("x","y")], method = 'bilinear')
data[which(data$month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("x","y")], method = 'bilinear')
data[which(data$month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("x","y")], method = 'bilinear')
data[which(data$month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("x","y")], method = 'bilinear')
data[which(data$month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("x","y")], method = 'bilinear')
data[which(data$month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data)
# Remove obs based on MLD
data2 <- data[-which(data$maxdepth > (data$MLD + 50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
summary(data2$MLD)
dim(data2)

data3 <- data2[!(data2$maxdepth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
#data3 <- data2 %>% drop_na(MLD) # Remove NAs
dim(data3) # 14047
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
save(data3, file = "Copepoda_PANGAEA_04_06_18.Rdata")


# 2) v4v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
# Create salinity vector
data$MLD <- NA
# Fill it according to month
data[which(data$month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$month == 1),c("x","y")], method = 'bilinear')
data[which(data$month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$month == 2),c("x","y")], method = 'bilinear')
data[which(data$month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$month == 3),c("x","y")], method = 'bilinear')
data[which(data$month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$month == 4),c("x","y")], method = 'bilinear')
data[which(data$month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$month == 5),c("x","y")], method = 'bilinear')
data[which(data$month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$month == 6),c("x","y")], method = 'bilinear')
data[which(data$month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$month == 7),c("x","y")], method = 'bilinear')
data[which(data$month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$month == 8),c("x","y")], method = 'bilinear')
data[which(data$month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$month == 9),c("x","y")], method = 'bilinear')
data[which(data$month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$month == 10),c("x","y")], method = 'bilinear')
data[which(data$month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$month == 11),c("x","y")], method = 'bilinear')
data[which(data$month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data)

# Remove obs based on MLD
data2 <- data[-which(data$maxdepth > (data$MLD + 50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
summary(data2$MLD)
data3 <- data2[!(data2$maxdepth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
#data3 <- data2 %>% drop_na(MLD) # Remove NAs
dim(data3) # 11899

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
save(data3, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------

##### 6) Prepare the data to be merged with OBIS & GBIF at v6

### v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
dir()
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)

# Re-name comment to basisofrecord
colnames(data)[10] <- c("basisofrecord")
### Remove event, abund, volume, depth, mindepth, elevation
colnames(data) # so columns 1,5,6,7,9,13
data <- data[,c(2,3,4,8,10,11,12,14:length(data))]
# Re-name maxdepth -> depth
colnames(data)[4] <- c("depth")

# Re-name DOI -> rightsholder
colnames(data)[6] <- c("rightsholder")
data$source <- "Cornisl&al._2018"

# Remove some species that were left: 
# Mixtocalanus_alter_female_Scolecithricella_altera_
# Mixtocalanus_alter_copepodites_Scolecithricella_altera_
# Idomenella_antarctica_Idomene_antarctica_
# Mospicalanus_schielae_copepodites
# Calanoida_indeterminata_copepodites                    
# Calanoida_indeterminata_female                         
# Calanoida_indeterminata_male 
# Augaptilidae_copepodites
dim(data)
data <- data[-which(data$species %in% c("Mixtocalanus_alter_female_Scolecithricella_altera_","Mixtocalanus_alter_copepodites_Scolecithricella_altera_",
			"Idomenella_antarctica_Idomene_antarctica_","Mospicalanus_schielae_copepodites","Calanoida_indeterminata_copepodites","Lophotrix_frontalis",
			"Calanoida_indeterminata_female","Calanoida_indeterminata_male","Augaptilidae_copepodites","Aetideidae_copepodites")),]
#
dim(data) # 45770

# Correct 'Pleuromamman_xiphias'...
levels(data$species) <- c(levels(data$species), "Pleuromamma_xiphias")
data$species[data$species == 'Pleuromamman_xiphias'] <- 'Pleuromamma_xiphias'

### Add genus name using strsplit on species name
data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)
head(data$genus)

# Add class and phylum
data$class <- "Hexanauplia"
data$phylum <- "Arthropoda"


### Add higher classification: family and class
# unique(data$genus)
### Try using the 'taxize' package
library("taxize")
#?classification
#taxize::classification("Calanus_propinquus", db = "worms")
#str(taxize::classification("Calanus_propinquus", db = "worms"))
#class(taxize::classification("Calanus_propinquus", db = "worms"))
#classif <- data.frame(taxize::classification("Calanus_propinquus", db = "worms")[[1]][1:2])
#classif # ok, that's how you can extract information
#rm(classif); gc()
# sp <- "Lophotrix_frontalis"
for(sp in unique(data$species)) {
		
		# Retrieve the full classif
		classif <- data.frame(taxize::classification(sp, db = "worms")[[1]][1:2])
		
		# Get family and order
		fam <- classif[classif$rank == "Family","name"]
		order <- classif[classif$rank == "Order","name"]
		
		# Provide to data
		data[data$species == sp,"family"] <- fam
		data[data$species == sp,"order"] <- order
		rm(order, fam, classif)
		gc()
		
} # eo for loop
1

### Check results
head(data)
unique(data$family)
colnames(data)

### Check colnames
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
ref <- colnames(get(load("Copepoda_15_05_18.Rdata")))
ref
colnames(data)
# OK save in proper dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
save(data, file = "Copepoda_PANGAEA_04_06_18.Rdata")



### ---------------------------------------------------------

### v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
dir()
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)
# Re-name comment to basisofrecord
colnames(data)[10] <- c("basisofrecord")
### Remove event, abund, volume, depth, mindepth, elevation
colnames(data) # so columns 1,5,6,7,9,13
data <- data[,c(2,3,4,8,10,11,12,14:length(data))]
# Re-name maxdepth -> depth
colnames(data)[4] <- c("depth")
# Re-name DOI -> rightsholder
colnames(data)[6] <- c("rightsholder")
data$source <- "Cornisl&al._2018"
dim(data) # 41019
data <- data[-which(data$species %in% c("Mixtocalanus_alter_female_Scolecithricella_altera_","Mixtocalanus_alter_copepodites_Scolecithricella_altera_",
			"Idomenella_antarctica_Idomene_antarctica_","Mospicalanus_schielae_copepodites","Calanoida_indeterminata_copepodites","Lophotrix_frontalis",
			"Calanoida_indeterminata_female","Calanoida_indeterminata_male","Augaptilidae_copepodites","Aetideidae_copepodites")),]
#
dim(data) # 40779

# Correct 'Pleuromamman_xiphias'...
levels(data$species) <- c(levels(data$species), "Pleuromamma_xiphias")
data$species[data$species == 'Pleuromamman_xiphias'] <- 'Pleuromamma_xiphias'

### Add genus name using strsplit on species name
data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)
head(data$genus)

# Add class and phylum
data$class <- "Hexanauplia"
data$phylum <- "Arthropoda"


### Add higher classification: family and class
library("taxize")

for(sp in unique(data$species)) {
		# Retrieve the full classif
		classif <- data.frame(taxize::classification(sp, db = "worms")[[1]][1:2])
		# Get family and order
		fam <- classif[classif$rank == "Family","name"]
		order <- classif[classif$rank == "Order","name"]
		# Provide to data
		data[data$species == sp,"family"] <- fam
		data[data$species == sp,"order"] <- order
		rm(order, fam, classif)
		gc()
} # eo for loop

### Check results
head(data)
unique(data$family)
unique(data$order)
### Check colnames
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
ref <- colnames(get(load("Copepoda_15_05_18.Rdata")))
ref
colnames(data)

# Change institutioncode to institution & remove bathymetry
colnames(data)[22] <- "institution"
data <- data[,c(1:10,12:length(data))]

### Save in proper directory
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
save(data, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ---------------------------------------------------------

### v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
dir()
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)

# Re-name comment to basisofrecord
colnames(data)[10] <- c("basisofrecord")
### Remove event, abund, volume, depth, mindepth, elevation
colnames(data) # so columns 1,5,6,7,9,13
data <- data[,c(2,3,4,8,10,11,12,14:length(data))]
# Re-name maxdepth -> depth
colnames(data)[4] <- c("depth")
# Re-name DOI -> rightsholder
colnames(data)[6] <- c("rightsholder")
data$source <- "Cornisl&al._2018"
dim(data) # 14047
data <- data[-which(data$species %in% c("Mixtocalanus_alter_female_Scolecithricella_altera_","Mixtocalanus_alter_copepodites_Scolecithricella_altera_",
			"Idomenella_antarctica_Idomene_antarctica_","Mospicalanus_schielae_copepodites","Calanoida_indeterminata_copepodites","Lophotrix_frontalis",
			"Calanoida_indeterminata_female","Calanoida_indeterminata_male","Augaptilidae_copepodites","Aetideidae_copepodites")),]
#
dim(data) # 14014

# Correct 'Pleuromamman_xiphias'...
levels(data$species) <- c(levels(data$species), "Pleuromamma_xiphias")
data$species[data$species == 'Pleuromamman_xiphias'] <- 'Pleuromamma_xiphias'

### Add genus name using strsplit on species name
data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)
head(data$genus)
# Kep the full NAs rows because of the MLD filter
data <- data[!is.na(data$species),]
dim(data) # 11878

# Add class and phylum
data$class <- "Hexanauplia"
data$phylum <- "Arthropoda"

### Add higher classification: family and class
library("taxize")
for(sp in unique(data$species)) {
		# Retrieve the full classif
		classif <- data.frame(taxize::classification(sp, db = "worms")[[1]][1:2])
		# Get family and order
		fam <- classif[classif$rank == "Family","name"]
		order <- classif[classif$rank == "Order","name"]
		# Provide to data
		data[data$species == sp,"family"] <- fam
		data[data$species == sp,"order"] <- order
		rm(order, fam, classif)
		gc()
} # eo for loop

### Check results
head(data)
unique(data$family)
unique(data$order)
### Check colnames
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
ref <- colnames(get(load("Copepoda_15_05_18.Rdata")))
ref
colnames(data)

# Change institutioncode to institution & remove MLD
colnames(data)[22] <- "institution"
data <- data[,c(1:12,14:length(data))]

# Save in dir()
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
save(data, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ---------------------------------------------------------

### v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
dir()
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)

# Re-name comment to basisofrecord
colnames(data)[10] <- c("basisofrecord")
### Remove event, abund, volume, depth, mindepth, elevation
colnames(data) # so columns 1,5,6,7,9,13
data <- data[,c(2,3,4,8,10,11,12,14:length(data))]
# Re-name maxdepth -> depth
colnames(data)[4] <- c("depth")
# Re-name DOI -> rightsholder
colnames(data)[6] <- c("rightsholder")
data$source <- "Cornisl&al._2018"
dim(data) # 11899
data <- data[-which(data$species %in% c("Mixtocalanus_alter_female_Scolecithricella_altera_","Mixtocalanus_alter_copepodites_Scolecithricella_altera_",
			"Idomenella_antarctica_Idomene_antarctica_","Mospicalanus_schielae_copepodites","Calanoida_indeterminata_copepodites","Lophotrix_frontalis",
			"Calanoida_indeterminata_female","Calanoida_indeterminata_male","Augaptilidae_copepodites","Aetideidae_copepodites")),]
#
dim(data) # 11870

# Correct 'Pleuromamman_xiphias'...
levels(data$species) <- c(levels(data$species), "Pleuromamma_xiphias")
data$species[data$species == 'Pleuromamman_xiphias'] <- 'Pleuromamma_xiphias'

### Add genus name using strsplit on species name
data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)
head(data$genus)
# Kep the full NAs rows because of the MLD filter
data <- data[!is.na(data$species),]
dim(data) # 10591

# Add class and phylum
data$class <- "Hexanauplia"
data$phylum <- "Arthropoda"

### Add higher classification: family and class
library("taxize")
for(sp in unique(data$species)) {
		# Retrieve the full classif
		classif <- data.frame(taxize::classification(sp, db = "worms")[[1]][1:2])
		# Get family and order
		fam <- classif[classif$rank == "Family","name"]
		order <- classif[classif$rank == "Order","name"]
		# Provide to data
		data[data$species == sp,"family"] <- fam
		data[data$species == sp,"order"] <- order
		rm(order, fam, classif)
		gc()
} # eo for loop

### Check results
head(data)
unique(data$family)
unique(data$order)
### Check colnames
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
ref <- colnames(get(load("Copepoda_15_05_18.Rdata")))
ref
colnames(data)

# Change institutioncode to institution & remove MLD + bathymetry
colnames(data)[23] <- "institution"
data <- data[,c(1:10,12,13,15:length(data))]

# Save in dir()
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
save(data, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------

##### 7) Last step, remove ducplicates for each v6 dataset

### v6-v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
# dir()
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)
# Change institutioncode
colnames(data)[21] <- "institution"

### Define an id and remove duplicates
data$id <- paste(data$x, data$y, data$depth, data$day, data$month, data$year, data$species, sep = "_")
length(unique(data$id)) # 18111 --> not 45770 so there should be 27659 duplicates
dim(data[duplicated(data$id),]) # 27659, yep
data[1:150,]
data2 <- data[!duplicated(data$id),]
dim(data2) # 18111
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
save(data2, file = "Copepoda_PANGAEA_04_06_18.Rdata")




### v6-v5.1v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)

### Define an id and remove duplicates
data$id <- paste(data$x, data$y, data$depth, data$day, data$month, data$year, data$species, sep = "_")
length(unique(data$id)) # 16788 --> not 40779 so there should be 23991 duplicates
dim(data[duplicated(data$id),]) # 23991, yep
data[1:150,]
data2 <- data[!duplicated(data$id),]
dim(data2) # 16788

save(data2, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### v6-v5.2v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)

### Define an id and remove duplicates
data$id <- paste(data$x, data$y, data$depth, data$day, data$month, data$year, data$species, sep = "_")
length(unique(data$id)) # 4396 --> not 11878 so there should be 27659 duplicates
dim(data[duplicated(data$id),]) # 7482, yep
data[1:150,]
data2 <- data[!duplicated(data$id),]
dim(data2) # 4396
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
save(data2, file = "Copepoda_PANGAEA_04_06_18.Rdata")


### v6-v5.2v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
data <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
colnames(data)

### Define an id and remove duplicates
data$id <- paste(data$x, data$y, data$depth, data$day, data$month, data$year, data$species, sep = "_")
length(unique(data$id)) # 4055 --> not 10591 so there should be 6536 duplicates
dim(data[duplicated(data$id),]) # 6536, yep
data[1:150,]
data2 <- data[!duplicated(data$id),]
dim(data2) # 4055

save(data2, file = "Copepoda_PANGAEA_04_06_18.Rdata")


















