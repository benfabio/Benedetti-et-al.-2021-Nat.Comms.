
##### 06/06/2018: R Script to 

### Aims to:
#	- Load the Foraminifera data from GBIF & OBIS 
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

### Latest update: 12/06/2018

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

### ----------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Load all datasets and cbind them

### For OBIS: 
group <- dir()
# g <- "Globigerinoidea_OBIS"
res <- lapply(group, function(g) {
			setwd(paste("/Users/fabiobenedetti/Desktop/Foraminifera/OBIS/",g,"/", sep = ""))
			d <- read.csv(paste(g,".csv", sep = ""), h = T, sep = ",")
			message(paste(dim(d), sep = " "))
			return(d)
}
) # eo lapply
# rbind ?
table <- do.call(rbind, res)
dim(table)
table[1:1000,]
rm(res) ; gc()

setwd("/Users/fabiobenedetti/Desktop/Foraminifera/OBIS/")
save(table, file = "Foraminifera_OBIS_07_06_18.Rdata")

### For GBIF
setwd("/Users/fabiobenedetti/Desktop/Foraminifera/GBIF/")
files <- dir()
# g <- "Globigerinoidea_OBIS"
res <- lapply(files, function(f) {
			d <- read.csv(f, h = T, sep = "\t")
			message(paste(dim(d), sep = " "))
			return(d)
}
) # eo lapply
# rbind ?
table <- do.call(rbind, res)
dim(table)
str(table)
table[1:100,]
rm(res) ; gc()

setwd("/Users/fabiobenedetti/Desktop/Foraminifera/GBIF/")
save(table, file = "Foraminifera_GBIF_07_06_18.Rdata")


### 2°) Load the OBIS and GBIF data and apply primary cleaning steps:

# OBIS
setwd("/Users/fabiobenedetti/Desktop/Foraminifera/OBIS/")
# dir()
d <- get(load("Foraminifera_OBIS_07_06_18.Rdata"))
dim(d) # 678'921
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 678'921
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 624'571
d3 <- subset(d2, year > 1800)
dim(d3) # 624'567
summary(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 602'672
summary(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 592'933
unique(d5$species) 

summary(d5)
str(d5)

# What are the metadata you can use to identify the sediment records:
# bibliographicCitation
### unique(d5$bibliographicCitation)
# But how many NAs? 
nrow(d5[d5$bibliographicCitation == "",]) # 5026 ; (5026/ nrow(d5))*100 less than 1% pf the data have no bibliographicCitation --> remove
str(d5$bibliographicCitation)
### !!! Use these levels to identify what comes from CTD or net samples ! 
# Strategy: View the bibliographicCitation levels that are associated to sediment corezs, and use some character codes to identify
# levels(d5$bibliographicCitation)
codes <- c("sediment","paleo","Hole","core","CLIMAP","ODP","temperature reconstruction","DSDP","Oligocene","Neogene","Miocene",
		"Pliocene","Holocene","Paleocene","Pleistocene","abundance of Hole","SST","sedimentological","GeoB","sediments",
		"stratigraphically","time slice","Sediment","reconstruction","last glacial maximum","Site GIK","fossils","Shipboard Scientific Party",
		"DSDP Site","lithic","meiofauna","Sedimentological","sediment corev","Stable oxygen isotope","Paleocene-Eocene",
		"Cretaceous","Maastrichtian","K/T","Benth","benthic","AMK21","Discovery Reports","NU2_trap","trap L2")

tokeep <- levels(d5$bibliographicCitation)[!grepl(paste(codes, collapse = "|"), levels(d5$bibliographicCitation))] # levels to keep
levels(d5$bibliographicCitation)[grepl(paste(codes, collapse = "|"), levels(d5$bibliographicCitation))]  # levels NOT to keep

### OK check these levels: 
# [3] "Darling,K.F., Wade,C.M., Stewart,I.A., Kroo Nature 405 (6782), 43-47 (2000)"
# [4] "Darling,K.F., Wade,C.M., Stewart,I.A., Kroon,D., D Nature 405 (6782), 43-47 (2000)"
### --> Keep
nrow(d5[d5$bibliographicCitation == "Darling,K.F., Wade,C.M., Stewart,I.A., Kroon,D., D Nature 405 (6782), 43-47 (2000)",])
# But have been removed by preceding steps

# [5] "Earland, A. 1934. Foraminifera. Part III. The Falklands sector of the Antarctic (excluding South Georgia). Discovery Reports Vol. X, pp. 1-208."
nrow(d5[d5$bibliographicCitation == "Earland, A. 1934. Foraminifera. Part III. The Falklands sector of the Antarctic (excluding South Georgia). Discovery Reports Vol. X, pp. 1-208.",])
# 75 obs
head(d5[d5$bibliographicCitation == "Earland, A. 1934. Foraminifera. Part III. The Falklands sector of the Antarctic (excluding South Georgia). Discovery Reports Vol. X, pp. 1-208.",])
d5[d5$bibliographicCitation == "Earland, A. 1934. Foraminifera. Part III. The Falklands sector of the Antarctic (excluding South Georgia). Discovery Reports Vol. X, pp. 1-208.","occurrenceRemarks"]
### To be removed

# [6] "Guillem Mateu - gmateu@marum.de"
dim(d5[d5$bibliographicCitation == "Guillem Mateu - gmateu@marum.de",])
# But have been removed by preceding steps

# [64] "Laboratory of bottom fauna, P.P.Shirshov Institute of Oceanology of Russian Academy of Science"     
dim(d5[d5$bibliographicCitation == "Laboratory of bottom fauna, P.P.Shirshov Institute of Oceanology of Russian Academy of Science",])
# But have been removed by preceding steps
                                                                                                                                                                               
# [65] "Romero, Oscar E; Boeckel, Babette; Donner, Barbara; Lavik, Gaute; Fischer, Gerhard; Wefer, Gerold (2002): Flux data of NU2_trap (Table A1), doi:10.1594/PANGAEA.115823"   
# sediment trap data --> remove

### Remove empty levels
d6 <- d5[-which(d5$bibliographicCitation == ""),]
dim(d5) ; dim(d6) # removed 5026

### And filter the levels 
d6 <- d6[which(d6$bibliographicCitation %in% tokeep),]
dim(d6) # 69'584 observations !

# Save in v2 dir
save(d6, file = "Foraminifera_OBIS_08_06_18.Rdata" )


### --------------------------------------------------------------------------------------

# GBIF
setwd("/Users/fabiobenedetti/Desktop/Foraminifera/GBIF/")
dir()

d <- get(load("Foraminifera_GBIF_07_06_18.Rdata"))
dim(d) # 1'806'792
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 1'796'304
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 1'199'260
summary(d2$year)
d3 <- subset(d2, year > 1800)
dim(d3) # 62410
d4 <- d3 %>% drop_na(depth)
dim(d4) # 1'024'099
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 977'415
unique(d5$species) 

# What are the metadata you can use to identify the sediment records:
colnames(d5)
# datasetkey ?
head(d5$datasetkey) # nay
unique(d5$datasetkey) # nay
# publishingorgkey
head(d5$publishingorgkey) # same
# basisofrecord
unique(d5$basisofrecord) # "FOSSIL_SPECIMEN"
sum(is.na(d5$basisofrecord)) # no NA?
summary(d5$basisofrecord) # OK, 159 fossils
# issue
unique(d5$issue)
# collectioncode
str(d5$collectioncode)
levels(d5$collectioncode)

##### Now, how to retrieve the datasets' name ?

### Check out the datasets function of the rgbif package
library("rgbif")
dataset <- rgbif::datasets(type = 'occurrence', data = 'all', uuid = "880d5c76-f762-11e1-a439-00145eb45e9a")
str(dataset)
dataset$data$title

### Okay, get the title of each datasetkey (uuid) and provide to d5 !
keys <- data.frame(unique(d5$datasetkey))
colnames(keys) <- "key"
head(keys)
keys$name <- NA
### With lapply, get all the names correspondiong to each key 
names <- lapply(keys$key, function(k) {
			# Get dataset info 
			dataset <- rgbif::datasets(type = 'occurrence', data = 'all', uuid = k)
			# get title & return
			title <- dataset$data$title
			return(title)
}
) #eo lapply, takes a few minutes
str(names)
# rbind
namestable <- do.call(rbind, names)
length(namestable) #  == dim(keys)

keys$name <- namestable
head(keys)

### Now, identify the keys that correspond to datasets you want to keep. Like above, use the codes
codes <- c("sediment","paleo","Hole","core","CLIMAP","ODP","temperature reconstruction","DSDP","Oligocene","Neogene","Miocene","Mg/Ca",
		"Pliocene","Holocene","Paleocene","Pleistocene","abundance of Hole","SST","sedimentological","GeoB","sediments","sub-ice fauna",
		"stratigraphically","time slice","Sediment","reconstruction","last glacial maximum","Site GIK","fossils","Shipboard Scientific Party",
		"DSDP Site","lithic","meiofauna","Sedimentological","sediment corev","Stable oxygen isotope","Paleocene-Eocene","Dissolution index",
		"Cretaceous","Maastrichtian","K/T","Benth","benthic","AMK21","Discovery Reports","NU2_trap","petroleum","Paleobiology","LGM")
		
names2keep <- unique(keys$name)[!grepl(paste(codes, collapse = "|"), unique(keys$name))] # levels to keep
length(names2keep)
names2keep

### Check the validity of these datasets:
# "Biota occurrence data from plankton surveys around New Zealand"
### --> already present in OBIS, will be removed ultimately, keep

# "(Table 2) Dissolution index (XDX) values and average test mass for four species of planktonic foraminifera"  
### --> from a German PhD thesis studying sediment forams, remove

# "NMNH Extant Specimen Records"
### Not sure, but could be ok

# "Electron Micrograph Database - Marine Specimens"
# https://gcmd.nasa.gov/KeywordSearch/Metadata.do?Portal=amd_au&MetadataView=Full&MetadataType=0&KeywordPath=&OrigMetadataNode=AADC&EntryId=em_database
# https://data.aad.gov.au/metadata/records/em_database
### Hard to tell, keep for now

# "Arctic Ocean Diversity" 
# http://www.arcodiv.org/watercolumn/protist/Protist.html
# http://www.arcodiv.org/Database/Plankton_datasets.html
### -> Forams from the water column, keep

# "Auckland Museum NZ  Marine Collection"
# Hard to tell, keep for now

# "Archives of the Arctic Seas Zooplankton"
# https://www.gbif.org/dataset/451eb991-c1f4-479f-b1f8-7c1b4e8f9114
# Data collected from scientific cruises from 1900-1973 in the Eurasian Arctic Seas, Polar Basin and the North-West Pacific.
### --> keep

# "Seasonal dynamics of sub-ice fauna below pack ice in the Arctic (Fram Strait)"
# https://www.sciencedirect.com/science/article/pii/S0967063705002724
### sub-ice fauna --> remove

# "(Table 2) Dissolution index (XDX) values and average test mass for four species of planktonic foraminifera" 
### --> from the same German PhD thesis studying sediment forams, remove

# "(Table 1) CTD data, mean spring and annual temperatures and salinities from World Ocean Atlas (WOA, 2001), and Mg/Ca data, sampled with a multinet at station 74 and 81 from RV Poseidon cruise 334"
# https://www.pangaea.de/expeditions/cr.php/Poseidon
# Cruise report: http://epic.awi.de/29147/1/Sch2006bx.pdf
# sediment sampling --> remove !

### Use names2keep to get keys2keep, which you will use to clean d5
keys2keep <- keys[which(keys$name %in% names2keep),"key"]
length(keys2keep) # 628 datasets

d6 <- d5[which(d5$datasetkey %in% keys2keep),]
dim(d6) # 74'345
d6$datasetName <- NA

### Provide dataset name !
for(n in names2keep) {
		# useless message:
		message(paste(n, sep = ""))
		# Get the corresponding key
		k <- keys[which(keys$name == n),"key"]
		# Supply to the observations in d6
		d6[which(d6$datasetkey == k),"datasetName"] <- as.character(n)
}

unique(d6$datasetName)


### Save !
save(d6, file = "Foraminifera_GBIF_08_06_18.Rdata" )


### ----------------------------------------------------------------------------------------------------------------------------------------------------------

##### 08/06/2018: Plot sampling effort

### Get coastline			
setwd("/Users/fabiobenedetti/Desktop/PostDocs/TARA/")					
# World coastline:
cl <- read.csv("world_coast.csv", h = TRUE)

coast <- list(
 	   # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	 	geom_polygon(aes(x = lon, y = lat), data = cl, fill= "grey60"),
  	  	geom_path(aes(x = lon, y = lat), data = cl, colour = "black", linetype = 1),
  	  	# appropriate projection
  	  	coord_quickmap(),
  	  	# remove extra space around the coast
  	  	scale_x_continuous(name = "Longitude", 
                     breaks = c(-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180), 
                     labels = c("-180°E","-150°E","-120°E","-90°E","-60°E","-30°E","0°E","30°E","60°E","90°E","120°E","150°E","180°E"), 
                     expand = c(0,0)), 
  
  		scale_y_continuous(name = "Latitude", 
                     breaks = c(-90,-60,-30,0,30,60,90), 
                     labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"),
                     expand = c(0,0)),
  		# dark gray background for the panel and legend
  	  	theme(
    		panel.background = element_rect(fill = "white"),  # background
    		legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70")
  		 )
)

setwd("/Users/fabiobenedetti/Desktop/Foraminifera/OBIS/")
data <- get(load("Foraminifera_OBIS_08_06_18.Rdata"))
unique(data$species)
#quartz()
map <- ggplot() + geom_point(aes(x = decimalLongitude, y = decimalLatitude), data = data, alpha = 0.6, pch = 21, colour = "black", fill = "#3288bd") + 
		coast + coord_quickmap() + theme_bw()

ggsave(plot = map, filename = "Sampling_OBIS_08_08_18.pdf", dpi = 300, height = 10, width = 8)



setwd("/Users/fabiobenedetti/Desktop/Foraminifera/GBIF/")
data <- get(load("Foraminifera_GBIF_08_06_18.Rdata"))
unique(data$species)
#quartz()
map <- ggplot() + geom_point(aes(x = decimallongitude, y = decimallatitude), data = data, alpha = 0.6, pch = 21, colour = "black", fill = "#d53e4f") + 
		coast + coord_quickmap() + theme_bw()

ggsave(plot = map, filename = "Sampling_GBIF_08_08_18.pdf", dpi = 300, height = 10, width = 8)


### Print Ralph Schiebel dataset's list from OBIS
setwd("/Users/fabiobenedetti/Desktop/Foraminifera/OBIS/")
data <- get(load("Foraminifera_OBIS_08_06_18.Rdata"))
unique(data$bibliographicCitation)

write.table(unique(data$bibliographicCitation)[grepl(paste("Schiebel", collapse = "|"), unique(data$bibliographicCitation))], "Datasets_RSchiebel_OBIS.txt", sep = "\t ")


### ----------------------------------------------------------------------------------------------------------------------------------------------------------


##### 12/06/2018: Apply data cleaning steps from v3 to v8

### 1°) Apply both v3 criteria: v3.1 (land mask) and v3.2 (distance to coast)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
data <- get(load("Foraminifera_GBIF_08_06_18.Rdata"))
dim(data)

### A) get the bathymetry raster
require("marmap")
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 15)
str(bathy)
# Convert to ddf
ras <- as.xyz(bathy)
colnames(ras) <- c("x","y","z")
coordinates(ras) <- ~ x + y
gridded(ras) <- TRUE
r <- raster(ras)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r # Looks alright

### extract() the bathymetry data and remove what is shallower than 200m
colnames(data)[c(18,17)] <- c("x","y")
data$bathymetry <- raster::extract(x = r, y = data[,c("x","y")], method = 'bilinear')
data2 <- data[-which(data$bathymetry >= -200),]
dim(data2) # 
summary(data2) # check NAs
data3 <- data2 %>% drop_na(bathymetry)
dim(data3) # 
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
save(data3, file = "Foraminifera_GBIF_12_06_18.Rdata")


### B) Get the dist2coast raster
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/dist2coast")
r <- get(load("dist2coast_raster_15min.Rdata"))
r #
### extract() the dist2coast data and remove what is within the first 25 kms
colnames(data)[c(18,17)] <- c("x","y")
data$dist2coast <- raster::extract(x = r, y = data[,c("x","y")], method = 'bilinear')
data2 <- data[-which(data$dist2coast < 25),]
summary(data2$dist2coast)
# Remove NAs --> land cells
data3 <- data2 %>% drop_na(dist2coast)
dim(data3)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
save(data3, file = "Foraminifera_GBIF_12_06_18.Rdata")


### --------------------------------------------------------------------------

### 2°) Apply v4 to the two v3 datasets

### First, get salinity masks
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
data <- get(load("Foraminifera_GBIF_12_06_18.Rdata"))
# Create salinity vector
data$salinityWOA13 <- NA
#require("lubridate")
#data$month <- lubridate::month(data$eventDate)
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
dim(data2)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3) # 
summary(data3$salinityWOA13)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
save(data3, file = "Foraminifera_GBIF_12_06_18.Rdata")


### B) Get v3.2 data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
data <- get(load("Foraminifera_GBIF_12_06_18.Rdata"))
# Create salinity vector
data$salinityWOA13 <- NA
#require("lubridate")
#data$month <- lubridate::month(data$eventDate)
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
dim(data2) # 
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3) # 
summary(data3$salinityWOA13)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
save(data3, file = "Foraminifera_GBIF_12_06_18.Rdata")


### --------------------------------------------------------------------------

### 3°) Apply v5.1 and v5.2 to the two v4 datasets

### A) Remove any obs that has a sampling depth below 500m (v5.1)
# 1) v4v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
data <- get(load("Foraminifera_GBIF_12_06_18.Rdata"))
summary(data$depth)
# Sometimes you still have negative depths for some reason... -> absolute
data$depth <- abs(data$depth)
# Remove obs with depth > 500
data2 <- data[!(data$depth > 500),]
summary(data2$depth)
dim(data2)
# Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
save(data2, file = "Foraminifera_GBIF_12_06_18.Rdata")


# 2) v4v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
data <- get(load("Foraminifera_GBIF_12_06_18.Rdata"))
summary(data$depth)
# Sometimes you still have negative depths for some reason... -> absolute
data$depth <- abs(data$depth)
# Remove obs with depth > 500
data2 <- data[!(data$depth > 500),]
summary(data2$depth)
dim(data2) # 
# Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
save(data2, file = "Foraminifera_GBIF_12_06_18.Rdata")


### B) Remove any obs that is belwo the monthly MLD (v5.2)

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
data <- get(load("Foraminifera_GBIF_12_06_18.Rdata"))
# Create MLD vector
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
data2 <- data[-which(data$depth > (data$MLD + 50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
summary(data2$MLD)
data3 <- data2[!(data2$depth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
#data3 <- data2 %>% drop_na(MLD) # Remove NAs
dim(data3) # 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
save(data3, file = "Foraminifera_GBIF_12_06_18.Rdata")


# 2) v4v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
data <- get(load("Foraminifera_GBIF_12_06_18.Rdata"))
# Create MLD vector
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
data2 <- data[-which(data$depth > (data$MLD + 50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
summary(data2$MLD)
data3 <- data2[!(data2$depth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
#data3 <- data2 %>% drop_na(MLD) # Remove NAs
dim(data3) # 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
save(data3, file = "Foraminifera_GBIF_12_06_18.Rdata")


### --------------------------------------------------------------------------

### 4°) Merge OBIS & GBIF data (v6)

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
dir()
g <- "Foraminifera"

colnamesOBIS <- c("species","genus","family","order","class","phylum",
			"x","y","eventDate","day","month","year","depth",
			"identifiedBy","recordedBy","rightsHolder","institutionID","basisOfRecord",
			"dist2coast","salinityWOA13")
			

colnamesGBIF <- c("species","genus","family","order","class","phylum",
				"x","y","eventdate","day","month","year","depth",
				"identifiedby","recordedby","rightsholder","institutioncode","basisofrecord",
				"dist2coast","salinityWOA13")			

message(paste("Merging the ", g, sep = ""))
# Load OBIS & GBIF files
obis <- get(load(paste(g,"_OBIS_12_06_18.Rdata", sep = "")))
gbif <- get(load(paste(g,"_GBIF_12_06_18.Rdata", sep = "")))
summary(obis)
summary(gbif)
# Add day
require("lubridate")
obis$day <- lubridate::day(obis$eventDate)
# Restrict to the selected colnames 
# setdiff(colnamesOBIS, colnames(obis))
# setdiff(colnamesGBIF, colnames(gbif))
obis <- obis[,colnamesOBIS]
gbif <- gbif[,colnamesGBIF]
# Rename the columns in both cases
# colnames(obis) ; colnames(gbif)
colnames(obis)[9] <- c("date")
colnames(obis)[14:18] <- c("identifiedby","recordedby","rightsholder","institution","basis")
#
colnames(gbif)[7:9] <- c("x","y","date")
colnames(gbif)[17:18] <- c("institution","basis")
# Add original database
gbif$source <- "GBIF"
obis$source <- "OBIS"
# Merge
table <- rbind(obis, gbif)
# dim(table)
# Go to proper v6 directory and save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
save(table, file = paste(g,"_12_06_18.Rdata", sep = ""))
		
# Clean
rm(table, obis, gbif)
gc()



### Check results (n obs and n species)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
table <- get(load("Foraminifera_12_06_18.Rdata"))
dim(table) # 
length(unique(table$species)) # 



### Also, extract the full species name list for the other 3 v6 datasets, and use the setdiff() function to extract the species that are not common 
# a) full species names list from v6-v5.1v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
table <- get(load("Foraminifera_12_06_18.Rdata"))
full_sp_list <- unique(table$species)
full_sp_list <- str_replace_all(as.character(full_sp_list), " ", "_")
full_sp_list


# b) species list from v6-v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
table2 <- get(load("Foraminifera_12_06_18.Rdata"))
sp_list2 <- unique(table2$species)
sp_list2 <- str_replace_all(as.character(sp_list2), " ", "_")
# setdiff
setdiff(sp_list2, full_sp_list) # 0 species  -> no need


# c) species list from v6-v5.2v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
table3 <- get(load("Foraminifera_12_06_18.Rdata"))
sp_list3 <- unique(table3$species)
sp_list3 <- str_replace_all(as.character(sp_list3), " ", "_")
# setdiff
setdiff(sp_list3, full_sp_list) # 0 species  -> no need

# d) species list from v6-v5.2v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
table4 <- get(load("Foraminifera_12_06_18.Rdata"))
sp_list4 <- unique(table4$species)
sp_list4 <- str_replace_all(as.character(sp_list4), " ", "_")
# setdiff
setdiff(sp_list4, full_sp_list) # 0 species  -> no need

### Just write full_sp_list
species <- data.frame(table %>% 
		   	group_by(species) %>%
		   	summarise(genus = unique(genus), family = unique(family)[1], order = unique(order)[1], class = unique(class)[1] ) 
)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
write.table(x = species, file = "foram_species_12_06_18.txt", sep = "\t") 



### --------------------------------------------------------------------------


### 5°) v7: nomenclature cleaner
### First, load corrected species names
setwd("/UP_home/fabioben/Desktop/OVERSEE/data")
names <- read.csv("foram_species_12_06_18.csv", h = TRUE, sep = ";")
dim(names)
str(names)
head(names)


### Second, load the observation data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
# Load the data
data <- get(load("Foraminifera_12_06_18.Rdata"))
# Add underscore to species names if you have not done so already
data$species <- str_replace_all(as.character(data$species), " ", "_")
# unique(data$species)
# class(data$species)

# Remove the species that are marked as 'remove'
toremove <- as.character(unique(names[which(names$action == "remove"),"species"])) # 23 spp to remove
# class(toremove)
data2 <- data[!(data$species %in% toremove),]
dim(data) ; dim(data2) 

# And correct species labels when necessary
tocorrect <- unique(names[which(names$action == "correct"),"correct_name"]) # 4 levels
### BEWARE: 'tocorrect' contains correct labels but only for the species to be corrected, 
###			 the current wrong species names will be in the 'wrongnames' string

### For each label to be corrected, find the wrong labels in 'data2' and replace them
data3 <- as.matrix(data2) # needed to replace factor levels...
for(sp in tocorrect) {
		# Useless message, again
		message(paste(sp, sep = ""))
		# Find the wrong names that correspond to 'sp', the real name
		wrongnames <- names[names$correct_name == sp,"species"]
		# Correct
		data3[which(data3[,"species"] %in% wrongnames),"species"] <- as.character(sp)
} # eo for loop

# Check data3 if necessary
data3 <- as.data.frame(data3)
dim(data3)
str(data3)
head(data3)
unique(data3$basis)
unique(data3$species) # 43
length(unique(data3$species) )

### Save in proper v7 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
save(data3, file = "Foraminifera_12_06_18.Rdata")



### --------------------------------------------------------------------------


### 6°) v8: duplicates cleaner

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
dir()

# Load
data <- get(load("Foraminifera_12_06_18.Rdata"))
head(data)
str(data)

data$x <- as.numeric(as.character(data$x))
data$y <- as.numeric(as.character(data$y))
data$depth <- as.numeric(as.character(data$depth))
#data$dist2coast <- as.numeric(as.character(data$dist2coast))
#data$salinityWOA13 <- as.numeric(as.character(data$salinityWOA13))
data$date <- as.Date(data$date)
data$day <- as.numeric(as.character(data$day))
data$month <- as.numeric(as.character(data$month))
data$year <- as.numeric(as.character(data$year))

summary(data)

### Round the coordinates to the 1/4° degree
data$xbin <- as.factor(round(data$x, 3))
data$ybin <- as.factor(round(data$y, 3))
#head(data)
#tail(data)

### Separate OBIS from GBIF again
obis <- data[data$source == "OBIS",]
gbif <- data[data$source == "GBIF",]
dim(obis) ; dim(gbif)

### 125 observations have NA value for day column in GBIF...check them and remove them
gbif[is.na(gbif$day),] # 
#unique(gbif[is.na(gbif$day),"species"]) # common copepods...
gbif$day[is.na(gbif$day)] <- 99
obis$day[is.na(obis$day)] <- 99

### Create an ID based on: xbin, ybin, day, month, year, species name
obis$id <- paste(obis$xbin, obis$ybin, obis$depth, obis$day, obis$month, obis$year, obis$species, sep = "_") 
gbif$id <- paste(gbif$xbin, gbif$ybin, gbif$depth, gbif$day, gbif$month, gbif$year, gbif$species, sep = "_") 

length(unique(obis$id)) ; length(obis$id) 

length(unique(gbif$id)) ; length(gbif$id)  #  

### Remove duplicates from both obis and gbif
gbif2 <- gbif[!duplicated(gbif$id),]
dim(gbif) ; dim(gbif2)

obis2 <- obis[!duplicated(obis$id),]
dim(obis) ; dim(obis2)

### Now, compare the ids between OBIS & GBIF and remove those in GBIF
# Basically, setdiff(bigFrame, smallFrame) gets you the extra records in the first table.

diff <- setdiff(gbif2$id, obis2$id)
# class(diff) # character
length(diff)
(length(diff) / nrow(gbif2))*100 # 17.68% are not duplicates

### Use 'diff' to identify the obs of gbif2 that are to be added to obis2
tokeep <- gbif2[gbif2$id %in% diff,]

data2 <- rbind(obis2, tokeep)
dim(data2) # 
# And check if there duplicates in data2$id
length(unique(data2$id)) # ok

### Great, use the code above to identify and remove duplicates across datasets that combine OBIS & GBIF
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.2/")
save(data2, file = "Foraminifera_12_06_18.Rdata" )
# Clean 
rm(data2, data, tokeep, obis2, gbif2, obis, gbif) ; gc()

# Check results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.1/")
data <- get(load("Foraminifera_12_06_18.Rdata"))
nrow(data)
length(unique(data$species))













