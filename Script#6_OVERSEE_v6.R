
##### 15/05/2018: R Script to merge OBIS and GBIF datasets

### Aims to:
#	- identify the columns that are common to both OBIS & GBIF datasets
#	- re-name some columns so both types of datasets have common names (species, x, y, month, day, year, depth etc.)
#	- merge OBIS and GBIF datasets, for each major zooplankton group
#	- do it for 4 types of datasets: v5.1v3.1, v5.1v3.2, v5.2v3.1, & v5.2v3.2

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 17/05/2018

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


##### 1°) Load a GBIF and a OBIS dataset, examine column names and identify the ones you will keep common
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
dir()

# Let's use the Thecosomata for a change
gbif <- get(load("Thecosomata_GBIF_04_05_18.Rdata"))
obis <- get(load("Thecosomata_OBIS_04_05_18.Rdata"))
dim(gbif) # 31'099    46
dim(obis) # 32'131   144

# Since you have way less columns in GBIF than OBIS, GBIF datasets will be the ones setting the columns
colnames(gbif)
# summary(gbif)
#  [1] "gbifid"                        "datasetkey"
#  [3] "occurrenceid"                  "kingdom"
#  [5] "phylum"                        "class"
#  [7] "order"                         "family"
#  [9] "genus"                         "species"
# [11] "infraspecificepithet"          "taxonrank"
# [13] "scientificname"                "countrycode"
# [15] "locality"                      "publishingorgkey"
# [17] "decimallatitude"               "decimallongitude"
# [19] "coordinateuncertaintyinmeters" "coordinateprecision"
# [21] "elevation"                     "elevationaccuracy"
# [23] "depth"                         "depthaccuracy"
# [25] "eventdate"                     "day"
# [27] "month"                         "year"
# [29] "taxonkey"                      "specieskey"
# [31] "basisofrecord"                 "institutioncode"
# [33] "collectioncode"                "catalognumber"
# [35] "recordnumber"                  "identifiedby"
# [37] "license"                       "rightsholder"
# [39] "recordedby"                    "typestatus"
# [41] "establishmentmeans"            "lastinterpreted"
# [43] "mediatype"                     "issue"
# [45] "bathymetry"                    "salinityWOA13"

unique(gbif$mediatype) # not to keep
summary(gbif$issue) # keep
summary(gbif$lastinterpreted) # not to keep
summary(gbif$establishmentmeans) # not to keep
summary(gbif$typestatus) # not to keep
summary(gbif$recordedby) # keep just in case
summary(gbif$rightsholder) # keep just in case
summary(gbif$license) # nope
summary(gbif$identifiedby) # keep just in case
summary(gbif$recordnumber) # nein
summary(gbif$catalognumber) # nein
summary(gbif$collectioncode) # not really helpful most of the time
summary(gbif$institutioncode) # keep !
summary(gbif$basisofrecord) # keep !
summary(gbif$taxonkey) # nope
summary(gbif$specieskey) # nope
summary(gbif$publishingorgkey) # nope

### OKAY, so final choices are: 
# species and upper classif (genus, family, order, class, phylum)
# decimallongitude & decimallatitude --> x, y
# eventdate, day, (even NAs), month, year
# depth
# identifiedby, recordedby, rightsholder, institutioncode, basisofrecord
# bathymetry and salinityWOA13 which will be replaced by MLD and dist2coast when using other criteria

### So 20 columns. Now, identify the equivalenbt in OBIS datasets
colnames(obis)
summary(obis$identifiedBy) # ok
summary(obis$recordedBy) 
summary(obis[,c("rights", "rightsHolder")])
summary(obis$institutionID)
summary(obis$basisOfRecord)
summary(obis$associatedReferences) # not really what I expected

summary(obis$references) # 
summary(obis$associatedReferences) 
summary(obis$bibliographicCitation)

summary(obis$eventRemarks) # nope
unique(obis$occurrenceRemarks)
summary(obis$source) # or 'resource_id' -> nope
summary(obis[,c("institutionCode","institutionID")])

# species and upper classif (genus, family, order, class, phylum) --> same colnames
# decimalLongitude & decimalLatitude --> x, y
# eventdate, day, (even NAs), month, year (need to add day with lubridate though)
# depth
# identifiedBy --> identifiedby
# recordedBy --> recordedby
# rightsHolder --> rightsholder
# institutionID --> institutioncode
# basisOfRecord --> basisofrecord
# bathymetry and salinityWOA13


##### 2°) Now that you have identified the columns to keep, restrict the current datasets to these columns and merge OBIS & GBIF datasets

### A) Merge the Copepoda OBIS files from v5.1v3.1, v5.1v3.2, v5.2v3.1, & v5.2v3.2

colnames <- c("species","genus","family","order","class","phylum",
			"decimalLongitude","decimalLatitude","eventDate","day","month","year","depth",
			"identifiedBy","recordedBy","rightsHolder","institutionID","basisOfRecord",
			"dist2coast","salinityWOA13")
					
### First, you need to combine all the Copepoda OBIS species datasets (not all have the same columns)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/Copepoda_OBIS_15_05_18/")
files <- dir()
# Get-load each file and return column of interest
# f <- files[1]
res <- lapply(files, function(f) {
			# Useless message
			message(paste("Doing ", f, sep = ""))
			d <- get(load(f))
			# colnames(d)
			# summary(d)
			# Restrict to columns of interest
			require("lubridate")
			d$day <- lubridate::day(d$eventDate) 
			# Sometimes, some species files do not have some of the selected colnames (such as rightsHolder or institutionID)
			# Need to identify them with setdiff() and then if length(missingcols) > 0, add those missing cols with NAs 
			missingcols <- setdiff(colnames, colnames(d))
			if( length(missingcols) > 0 ) {
				d[,missingcols] <- NA
			}
			
			# Then restrict to selected column names
			dd <- d[,colnames]
			# Return
			return(dd)
}
) # eo lapply

table <- do.call(rbind, res)
dim(table)
summary(table) # ok, save it
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
save(table, file = "Copepoda_OBIS_15_05_18.Rdata")
### Repeat for the 3 other v5



### B) Merge OBIS & GBIF files for each zooplankton groups

# a) v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
dir()
groups <- c("Copelata","Copepoda","Ctenophora","Cubozoa","Euphausiidae","Gymnosomata","Hydrozoa","Hyperiidea","Myodocopina",
			"Mysidae","Sagittoidea","Scyphozoa","Thaliacea","Thecosomata")

colnamesOBIS <- c("species","genus","family","order","class","phylum",
			"decimalLongitude","decimalLatitude","eventDate","day","month","year","depth",
			"identifiedBy","recordedBy","rightsHolder","institutionID","basisOfRecord",
			"bathymetry","salinityWOA13")
			

colnamesGBIF <- c("species","genus","family","order","class","phylum",
				"decimallongitude","decimallatitude","eventdate","day","month","year","depth",
				"identifiedby","recordedby","rightsholder","institutioncode","basisofrecord",
				"bathymetry","salinityWOA13")			

### For each zooplankton group
# g <- "Hyperiidea"
for(g in groups) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_04_05_18.Rdata", sep = "")))
		gbif <- get(load(paste(g,"_GBIF_04_05_18.Rdata", sep = "")))
		#summary(obis)
		#summary(gbif)
		# Add day
		require("lubridate")
		obis$day <- lubridate::day(obis$eventDate)
		# Restrict to the selected colnames 
		# setdiff(colnamesOBIS, colnames(obis))
		obis <- obis[,colnamesOBIS]
		gbif <- gbif[,colnamesGBIF]
		
		# Rename the columns in both cases
		# colnames(obis) ; colnames(gbif)
		colnames(obis)[7:9] <- c("x","y","date")
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
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
		save(table, file = paste(g,"_15_05_18.Rdata", sep = ""))
		
		# Clean
		rm(table, obis, gbif)
		gc()
	
} # eo for loop

### And then just restrict the OBIS files of the Alciopidae, Lopadorrhynchidae, P.avirostris, Podonidae, Tomopteridae, Typhloscolecidae
groups2 <- c("Alciopidae","Lopadorrhynchidae","P.avirostris","Podonidae","Tomopteridae","Typhloscolecidae")
# g <- "Alciopidae"
for(g in groups2) {
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_04_05_18.Rdata", sep = "")))
		# Add day
		require("lubridate")
		obis$day <- lubridate::day(obis$eventDate)
		# Restrict to the selected colnames 
		# setdiff(colnamesOBIS, colnames(obis))
		obis <- obis[,colnamesOBIS]
		# Rename the columns in both cases
		# colnames(obis)
		colnames(obis)[7:9] <- c("x","y","date")
		colnames(obis)[14:18] <- c("identifiedby","recordedby","rightsholder","institution","basis")
		# Add original database
		obis$source <- "OBIS"
		# Go to proper v6 directory and save
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
		save(obis, file = paste(g,"_15_05_18.Rdata", sep = ""))
		# Check it can be loaded again...
		# ddd <- get(load(paste(g,"_15_05_18.Rdata", sep = ""))) 
		# str(ddd)
		# Clean
		rm(obis) ; gc()
	
} # eo for loop


### Gut, repeat for other v6 datasets

### ------------------------------------------------------------------------------------------------------------------------

# b) v5.1v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
dir()
groups <- c("Copelata","Copepoda","Ctenophora","Cubozoa","Euphausiidae","Gymnosomata","Hydrozoa","Hyperiidea","Myodocopina",
			"Mysidae","Sagittoidea","Scyphozoa","Thaliacea","Thecosomata")

# Replace bathymetry by dist2coast because of v3.2
colnamesOBIS <- c("species","genus","family","order","class","phylum",
			"decimalLongitude","decimalLatitude","eventDate","day","month","year","depth",
			"identifiedBy","recordedBy","rightsHolder","institutionID","basisOfRecord",
			"dist2coast","salinityWOA13")
			

colnamesGBIF <- c("species","genus","family","order","class","phylum",
				"decimallongitude","decimallatitude","eventdate","day","month","year","depth",
				"identifiedby","recordedby","rightsholder","institutioncode","basisofrecord",
				"dist2coast","salinityWOA13")			

### For each zooplankton group
# g <- "Sagittoidea"
for(g in groups) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_04_05_18.Rdata", sep = "")))
		gbif <- get(load(paste(g,"_GBIF_04_05_18.Rdata", sep = "")))
		# summary(obis)
		# summary(gbif)
		# Add day
		require("lubridate")
		obis$day <- lubridate::day(obis$eventDate)
		# Restrict to the selected colnames 
		# setdiff(colnamesOBIS, colnames(obis))
		obis <- obis[,colnamesOBIS]
		gbif <- gbif[,colnamesGBIF]
		
		# Rename the columns in both cases
		# colnames(obis) ; colnames(gbif)
		colnames(obis)[7:9] <- c("x","y","date")
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
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
		save(table, file = paste(g,"_15_05_18.Rdata", sep = ""))
		
		# Clean
		rm(table, obis, gbif)
		gc()
	
} # eo for loop


### And then just restrict the OBIS files of the Alciopidae, Lopadorrhynchidae, P.avirostris, Podonidae, Tomopteridae, Typhloscolecidae
groups2 <- c("Alciopidae","Lopadorrhynchidae","P.avirostris","Podonidae","Tomopteridae","Typhloscolecidae")

for(g in groups2) {
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_04_05_18.Rdata", sep = "")))
		# Add day
		require("lubridate")
		obis$day <- lubridate::day(obis$eventDate)
		# Restrict to the selected colnames 
		# setdiff(colnamesOBIS, colnames(obis))
		obis <- obis[,colnamesOBIS]
		# Rename the columns in both cases
		# colnames(obis) ; colnames(gbif)
		colnames(obis)[7:9] <- c("x","y","date")
		colnames(obis)[14:18] <- c("identifiedby","recordedby","rightsholder","institution","basis")
		# Add original database
		obis$source <- "OBIS"
		# Go to proper v6 directory and save
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
		save(obis, file = paste(g,"_15_05_18.Rdata", sep = ""))
		# Clean
		rm(obis) ; gc()
	
} # eo for loop



### ------------------------------------------------------------------------------------------------------------------------

# c) v5.2v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
dir()
groups <- c("Copelata","Copepoda","Ctenophora","Cubozoa","Euphausiidae","Gymnosomata","Hydrozoa","Hyperiidea","Myodocopina",
			"Mysidae","Sagittoidea","Scyphozoa","Thaliacea","Thecosomata")

# Replace bathymetry by dist2coast because of v3.2
colnamesOBIS <- c("species","genus","family","order","class","phylum",
			"decimalLongitude","decimalLatitude","eventDate","day","month","year","depth",
			"identifiedBy","recordedBy","rightsHolder","institutionID","basisOfRecord",
			"bathymetry","salinityWOA13")
			

colnamesGBIF <- c("species","genus","family","order","class","phylum",
				"decimallongitude","decimallatitude","eventdate","day","month","year","depth",
				"identifiedby","recordedby","rightsholder","institutioncode","basisofrecord",
				"bathymetry","salinityWOA13")			

### For each zooplankton group
# g <- "Copepoda"
for(g in groups) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_15_05_18.Rdata", sep = "")))
		gbif <- get(load(paste(g,"_GBIF_15_05_18.Rdata", sep = "")))
		# summary(gbif)
		# summary(obis)
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
		colnames(obis)[7:9] <- c("x","y","date")
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
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
		save(table, file = paste(g,"_15_05_18.Rdata", sep = ""))
		
		# Clean
		rm(table, obis, gbif)
		gc()
	
} # eo for loop


### And then just restrict the OBIS files of the Alciopidae, Lopadorrhynchidae, P.avirostris, Podonidae, Tomopteridae, Typhloscolecidae
groups2 <- c("Alciopidae","Lopadorrhynchidae","P.avirostris","Podonidae","Tomopteridae","Typhloscolecidae")

for(g in groups2) {
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_15_05_18.Rdata", sep = "")))
		# Add day
		require("lubridate")
		obis$day <- lubridate::day(obis$eventDate)
		# Restrict to the selected colnames 
		# setdiff(colnamesOBIS, colnames(obis))
		obis <- obis[,colnamesOBIS]
		# Rename the columns in both cases
		# colnames(obis) ; colnames(gbif)
		colnames(obis)[7:9] <- c("x","y","date")
		colnames(obis)[14:18] <- c("identifiedby","recordedby","rightsholder","institution","basis")
		# Add original database
		obis$source <- "OBIS"
		# Go to proper v6 directory and save
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
		save(obis, file = paste(g,"_15_05_18.Rdata", sep = ""))
		# Clean
		rm(obis) ; gc()
	
} # eo for loop



### ------------------------------------------------------------------------------------------------------------------------

# d) v5.2v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
dir()
groups <- c("Copelata","Copepoda","Ctenophora","Cubozoa","Euphausiidae","Gymnosomata","Hydrozoa","Hyperiidea","Myodocopina",
			"Mysidae","Sagittoidea","Scyphozoa","Thaliacea","Thecosomata")

# Replace bathymetry by dist2coast because of v3.2
colnamesOBIS <- c("species","genus","family","order","class","phylum",
			"decimalLongitude","decimalLatitude","eventDate","day","month","year","depth",
			"identifiedBy","recordedBy","rightsHolder","institutionID","basisOfRecord",
			"dist2coast","salinityWOA13")
			

colnamesGBIF <- c("species","genus","family","order","class","phylum",
				"decimallongitude","decimallatitude","eventdate","day","month","year","depth",
				"identifiedby","recordedby","rightsholder","institutioncode","basisofrecord",
				"dist2coast","salinityWOA13")			

### For each zooplankton group
# g <- "Copepoda"
for(g in groups) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_15_05_18.Rdata", sep = "")))
		gbif <- get(load(paste(g,"_GBIF_15_05_18.Rdata", sep = "")))
		# summary(gbif)
		# summary(obis)
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
		colnames(obis)[7:9] <- c("x","y","date")
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
		save(table, file = paste(g,"_15_05_18.Rdata", sep = ""))
		
		# Clean
		rm(table, obis, gbif)
		gc()
	
} # eo for loop


### And then just restrict the OBIS files of the Alciopidae, Lopadorrhynchidae, P.avirostris, Podonidae, Tomopteridae, Typhloscolecidae
groups2 <- c("Alciopidae","Lopadorrhynchidae","P.avirostris","Podonidae","Tomopteridae","Typhloscolecidae")

for(g in groups2) {
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
		message(paste("Merging the ", g, sep = ""))
		# Load OBIS & GBIF files
		obis <- get(load(paste(g,"_OBIS_15_05_18.Rdata", sep = "")))
		# Add day
		require("lubridate")
		obis$day <- lubridate::day(obis$eventDate)
		# Restrict to the selected colnames 
		# setdiff(colnamesOBIS, colnames(obis))
		obis <- obis[,colnamesOBIS]
		# Rename the columns in both cases
		# colnames(obis) ; colnames(gbif)
		colnames(obis)[7:9] <- c("x","y","date")
		colnames(obis)[14:18] <- c("identifiedby","recordedby","rightsholder","institution","basis")
		# Add original database
		obis$source <- "OBIS"
		# Go to proper v6 directory and save
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
		save(obis, file = paste(g,"_15_05_18.Rdata", sep = ""))
		# Clean
		rm(obis) ; gc()
	
} # eo for loop




### ------------------------------------------------------------------------------------------------------------------------

### Summarise resulting number of occurrences and species per grouo and v6 dataset
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
files <- dir()
# f <- files[10]
for(f in files) {
		# Useless message
		message(paste(f, sep = ""))
		data <- get(load(f))
		# str(d)
		message(paste("n = ",nrow(data), sep = ""))
		message(paste("nsp = ", length(unique(data$species)), sep = ""))
} # eo for 

### And for sum of species: 
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table)
length(unique(table$species))


### ----------------------------------------------------------------------------------------------------------------------------

##### 16/05/2018: Map the sampling effort (and species richness) in 1/4° degree cells for each v6 dataset

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
library("ggplot2")
library("RColorBrewer")

### Load coastline for mapping afterwards
cl <- read.csv("world_coast.csv", h = T)

coast <- list(
  # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  geom_polygon(aes(x = lon, y = lat), data = cl, fill = "grey50"),
  geom_path(aes(x = lon, y = lat), data = cl, colour = "black", linetype = 1, size = 0.5),
  # appropriate projection
  coord_quickmap(),
  # remove extra space around the coast
  scale_x_continuous(name = "Longitude", 
                     breaks = c(-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180), 
                     labels = c("-180°E","-150°E","-120°E","-90°E","-60°E","-30°E","0°","30°E","60°E","90°E","120°E","150°E","180°E"), 
                     expand = c(0,0)), 
  
  scale_y_continuous(name = "Latitude", 
                     breaks = c(-90,-60,-30,0,30,60,90), 
                     labels = c("-90°N","-60°N","-30°N","0°","30°N","60°N","90°N"),
                     expand = c(0,0)),
  # dark gray background for the panel and legend
  theme(
    panel.background = element_rect(fill = "white"),  # background
    legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "white")
  )
)


### Load and plot
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table)
length(unique(table$species))


### Bin x and y in 1/4°, use dplyr() to compute sampling effort and then map with 'coast'
colnames(table)
table$xbin <- round(table$x, digits = 0.1)
table$ybin <- round(table$y, digits = 0.1)
summary(table)
table[1:1000,c("x","y","xbin","ybin")]
table$id <- paste(table$xbin, table$ybin, sep = "_")

### dplyr
require("dplyr")
require("viridis")
table2 <- data.frame(table %>% 
			group_by(id) %>% 
			summarise(x = unique(xbin), y = unique(ybin), n = n(), richness = length(unique(species)))
) # eo ddf
dim(table2)
head(table2)

map_sampling <- ggplot() + geom_tile(aes(x = x, y = y, fill = log(n)), data = table2) + 
		scale_fill_distiller(name = "Sampling effort\nlog(n occurrences)", palette = "YlGnBu") + coast + coord_quickmap()

map_rich <- ggplot() + geom_tile(aes(x = x, y = y, fill = richness), data = table2) + 
		scale_fill_distiller(name = "Species richness", palette = "YlGnBu") + coast + coord_quickmap()

# Save
ggsave(plot = map_sampling, filename = "map_effort_v6v5.2v3.2.pdf", dpi = 300, width = 12, height = 9)
ggsave(plot = map_rich, filename = "map_richness_v6v5.2v3.2.pdf", dpi = 300, width = 12, height = 9)


# OKAY


### ----------------------------------------------------------------------------------------------------------------------------

##### 17/05/2018: For each group, load the species names list for the v6-v5.1v3.2 dataset (the one with highest diversity, is likely to contain all the species found in the other v6 datasets) and save them as .txt files to manually clean them.

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
# dir()
files <- dir()[grep(".Rdata", dir())] # don't take the .pdf...
# f <- files[3]
for(f in files) {
	
	# Lock & load
	setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
	data <- get(load(f))
	require("dplyr")
	# Retrive classif for each species name
	species <- data.frame(data %>% 
			   	group_by(species) %>%
			   	summarise(genus = unique(genus), family = unique(family)[1], order = unique(order)[1], class = unique(class)[1] ) 
				)
				
	# Save as .txt file
	setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
	write.table(x = species, file = paste(str_replace(f, "15_05_18.Rdata", "17_05_18.txt"), sep = "" ) )
	
} # eo for loop

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
dir()

### Also, extract the full species name list for the other 3 v6 datasets, and use the setdiff() function to extract the species that are not common 

# a) full species names list from v6-v5.1v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
files <- dir()[grep(".Rdata", dir())] # don't take the .pdf...
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table) # 2'410'080
length(unique(table$species)) # 2590

full_sp_list <- unique(table$species)
full_sp_list <- str_replace_all(as.character(full_sp_list), " ", "_")





# b) species list from v6-v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
files <- dir()[grep(".Rdata", dir())] # don't take the .pdf...
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data)
}
) # eo lapply
table2 <- do.call(rbind, res)
dim(table2) # 1'539'472
length(unique(table2$species)) # 2192
sp_list2 <- unique(table2$species)
sp_list2 <- str_replace_all(as.character(sp_list2), " ", "_")
table2$species <- str_replace_all(as.character(table2$species), " ", "_")

# setdiff
setdiff(sp_list2, full_sp_list) # 34 species !
diff <- setdiff(sp_list2, full_sp_list)
species2 <- data.frame(table2[which(table2$species %in% diff),] %>% 
		   	group_by(species) %>%
		   	summarise(genus = unique(genus), family = unique(family)[1], order = unique(order)[1], class = unique(class)[1] ) 
)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
write.table(x = species2, file = "species_v6-v5.1v3.1.txt", sep = "\t") 


# c) species list from v6-v5.2v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
files <- dir()[grep(".Rdata", dir())] # don't take the .pdf...
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data)
}
) # eo lapply
table3 <- do.call(rbind, res)
dim(table3) # 1077560
length(unique(table3$species)) # 1427
sp_list3 <- unique(table3$species)
sp_list3 <- str_replace_all(as.character(sp_list3), " ", "_")
table3$species <- str_replace_all(as.character(table3$species), " ", "_")

# setdiff
setdiff(sp_list3, full_sp_list) # 17 species !
diff <- setdiff(sp_list3, full_sp_list)
species3 <- data.frame(table3[which(table3$species %in% diff),] %>% 
		   	group_by(species) %>%
		   	summarise(genus = unique(genus), family = unique(family)[1], order = unique(order)[1], class = unique(class)[1] ) 
)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
write.table(x = species3, file = "species_v6-v5.2v3.1.txt", sep = "\t") 




# d) species list from v6-v5.2v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
files <- dir()[grep(".Rdata", dir())] # don't take the .pdf...
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data)
}
) # eo lapply
table4 <- do.call(rbind, res)
dim(table4) # 1877047
length(unique(table4$species)) # 1967
sp_list4 <- unique(table4$species)
sp_list4 <- str_replace_all(as.character(sp_list4), " ", "_")

# setdiff
setdiff(sp_list4, full_sp_list)
table4$species <- str_replace_all(as.character(table4$species), " ", "_")

# setdiff
setdiff(sp_list4, full_sp_list) # 4 species !
diff <- setdiff(sp_list4, full_sp_list)
species4 <- data.frame(table4[which(table4$species %in% diff),] %>% 
		   	group_by(species) %>%
		   	summarise(genus = unique(genus), family = unique(family)[1], order = unique(order)[1], class = unique(class)[1] ) 
)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
write.table(x = species4, file = "species_v6-v5.2v3.2.txt", sep = "\t") 













