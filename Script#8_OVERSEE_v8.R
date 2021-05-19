
##### 06/06/2018: R Script to identify and remove the duplicate observations in the v7 datasets 

### Aims to
#	- load the v7 dataset, identify whether you have both OBIS & GBIF (unique(source) > 1), and separate dataset if you need to
#	- identify duplicate observations uqing and is based on coords, depth, date and species name
# 	- remove the duplicates within each dataset
#	- remove the duplicates between the 2 data sources (OBIS & GBIF) when needed
#	- rbind and save in v8 dir

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 25/07/2018

library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")

### ----------------------------------------------------------------------------------------------------------------------------

##### Script to identify and remive duplicates observations between GIBF and OBIS

# https://stackoverflow.com/questions/3171426/compare-two-data-frames-to-find-the-rows-in-data-frame-1-that-are-not-present-in

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.2/")
dir()

# Load
data <- get(load("Copepoda_05_06_18.Rdata"))
head(data)
str(data)

data$x <- as.numeric(as.character(data$x))
data$y <- as.numeric(as.character(data$y))
data$depth <- as.numeric(as.character(data$depth))
data$dist2coast <- as.numeric(as.character(data$dist2coast))
data$salinityWOA13 <- as.numeric(as.character(data$salinityWOA13))
data$date <- as.Date(data$date)
data$day <- as.numeric(as.character(data$day))
data$month <- as.numeric(as.character(data$month))
data$year <- as.numeric(as.character(data$year))

summary(data)

### Round the coordinates to the 1/4° degree
data$xbin <- as.factor(round(data$x, 3))
data$ybin <- as.factor(round(data$y, 3))
head(data)
tail(data)

### Separate OBIS from GBIF again
obis <- data[data$source == "OBIS",]
gbif <- data[data$source == "GBIF",]

dim(obis) ; dim(gbif)
summary(obis)
summary(gbif)

### 125 observations have NA value for day column in GBIF...check them and remove them
gbif[is.na(gbif$day),] # 
unique(gbif[is.na(gbif$day),"species"]) # common copepods...
gbif$day[is.na(gbif$day)] <- 99
obis$day[is.na(obis$day)] <- 99

summary(obis)
summary(gbif)

### Create an ID based on: xbin, ybin, day, month, year, species name
obis$id <- paste(obis$xbin, obis$ybin, obis$depth, obis$day, obis$month, obis$year, obis$species, sep = "_") 
gbif$id <- paste(gbif$xbin, gbif$ybin, gbif$depth, gbif$day, gbif$month, gbif$year, gbif$species, sep = "_") 

length(obis$id) ; length(unique(obis$id)) # 793409 / 1342103 --> 40% are ducplicates
obis$id[1:1000]

length(gbif$id) ; length(unique(gbif$id)) #  521829/ 598290 , 13% are duplicates
gbif$id[1:1000]

### Remove ducplicates from both obis and gbif
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
dim(data2)
# And check if there duplicates in data2$id
length(unique(data2$id)) # ok

### Great, use the code above to identify and remove duplicates across datasets that combine OBIS & GBIF


### ----------------------------------------------------------------------------------------------------------------------------

##### 1°) v7-5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.1/")
dir()

# For testing:
# f <- dir()[2]
for(f in dir()) {
	
		# Useless message()
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.1/")
		message(paste("Doing ", f, sep = ""))
		data <- get(load(f))
		
		# Converting some stuff if needed
		data$x <- as.numeric(as.character(data$x))
		data$y <- as.numeric(as.character(data$y))
		data$depth <- as.numeric(as.character(data$depth))
		if( f == "Thecosomata_MAREDAT_05_06_18.Rdata" ) {
				colnames(data)[7] <- "date"
		} else {
				data$date <- as.Date(data$date)
		} # eo if else loop
		
		data$day <- as.numeric(as.character(data$day))
		data$month <- as.numeric(as.character(data$month))
		data$year <- as.numeric(as.character(data$year))
		
		# Round the coordinates to the 1/4° degree
		data$xbin <- as.factor(round(data$x, 3))
		data$ybin <- as.factor(round(data$y, 3))
		# head(data)
		# tail(data)
		
		# If length(unique(data$source)) > 1 -> you have OBIS, GBIF
		if( length(unique(data$source)) > 1 ) {
			
				# Separate OBIS from GBIF again
				obis <- data[data$source == "OBIS",]
				gbif <- data[data$source == "GBIF",]
				# Convert NA days to 999
				gbif$day[is.na(gbif$day)] <- 999
				obis$day[is.na(obis$day)] <- 999
				# Create an ID based on: xbin, ybin, depth, day, month, year, species name
				obis$id <- paste(obis$xbin, obis$ybin, obis$depth, obis$day, obis$month, obis$year, obis$species, sep = "_") 
				gbif$id <- paste(gbif$xbin, gbif$ybin, gbif$depth, gbif$day, gbif$month, gbif$year, gbif$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				gbif2 <- gbif[!duplicated(gbif$id),]
				obis2 <- obis[!duplicated(obis$id),]
				# nrow(gbif2) ; nrow(gbif)
				# nrow(obis2) ; nrow(obis)
				# Now, compare the ids between OBIS & GBIF and remove those in GBIF
				diff <- setdiff(gbif2$id, obis2$id)
				# length(diff)
				# (length(diff) / nrow(gbif2))*100 # % of non duplicates
				# Use 'diff' to identify the obs of gbif2 that are to be added to obis2
				tokeep <- gbif2[gbif2$id %in% diff,]
				# rbind
				data2 <- rbind(obis2, tokeep)
				# nrow(data2) / nrow(data)
				# length(unique(data2$id)) # ok
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.1/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18") )
				# Clean 
				rm(data2, data, tokeep, obis2, gbif2, obis, gbif) ; gc()
			
		} else {
			
				# Simply apply the same steps as above but on 'data'
				data$day[is.na(data$day)] <- 999
				# Create an ID based on: xbin, ybin, day, month, year, species name
				data$id <- paste(data$xbin, data$ybin, data$depth, data$day, data$month, data$year, data$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				data2 <- data[!duplicated(data$id),]
				# dim(data2) ; dim(data)
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.1/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18"))
				# Clean 
				rm(data2, data) ; gc()
			
		} # eo if else loop
		
	
} # eo for loop


### Check results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.1/")
# dir()
for(f in dir()) {
		# Useless message
		message(paste(str_replace(f, "_06_06_18.Rdata", ""), sep = ""))
		d <- get(load(f))
		message(paste(nrow(d), sep = ""))
		message(paste(length(unique(d$species)), sep = ""))
		message(paste("", sep = ""))
		rm(d)
} # eo for loop



### ----------------------------------------------------------------------------------------------------------------------------

##### 2°) v7-5.1v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.2/")
dir()

# For testing:
# f <- dir()[2]
for(f in dir()) {
	
		# Useless message()
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.2/")
		message(paste("Doing ", f, sep = ""))
		data <- get(load(f))
		
		# Converting some stuff if needed
		data$x <- as.numeric(as.character(data$x))
		data$y <- as.numeric(as.character(data$y))
		data$depth <- as.numeric(as.character(data$depth))
		if( f == "Thecosomata_MAREDAT_05_06_18.Rdata" ) {
				colnames(data)[7] <- "date"
		} else {
				data$date <- as.Date(data$date)
		} # eo if else loop
		
		data$day <- as.numeric(as.character(data$day))
		data$month <- as.numeric(as.character(data$month))
		data$year <- as.numeric(as.character(data$year))
		
		# Round the coordinates to the 1/4° degree
		data$xbin <- as.factor(round(data$x, 3))
		data$ybin <- as.factor(round(data$y, 3))
		# head(data)
		# tail(data)
		
		# If length(unique(data$source)) > 1 -> you have OBIS, GBIF
		if( length(unique(data$source)) > 1 ) {
			
				# Separate OBIS from GBIF again
				obis <- data[data$source == "OBIS",]
				gbif <- data[data$source == "GBIF",]
				# Convert NA days to 999
				gbif$day[is.na(gbif$day)] <- 999
				obis$day[is.na(obis$day)] <- 999
				# Create an ID based on: xbin, ybin, depth, day, month, year, species name
				obis$id <- paste(obis$xbin, obis$ybin, obis$depth, obis$day, obis$month, obis$year, obis$species, sep = "_") 
				gbif$id <- paste(gbif$xbin, gbif$ybin, gbif$depth, gbif$day, gbif$month, gbif$year, gbif$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				gbif2 <- gbif[!duplicated(gbif$id),]
				obis2 <- obis[!duplicated(obis$id),]
				# nrow(gbif2) ; nrow(gbif)
				# nrow(obis2) ; nrow(obis)
				# Now, compare the ids between OBIS & GBIF and remove those in GBIF
				diff <- setdiff(gbif2$id, obis2$id)
				# length(diff)
				# (length(diff) / nrow(gbif2))*100 # % of non duplicates
				# Use 'diff' to identify the obs of gbif2 that are to be added to obis2
				tokeep <- gbif2[gbif2$id %in% diff,]
				# rbind
				data2 <- rbind(obis2, tokeep)
				# nrow(data2) / nrow(data)
				# length(unique(data2$id)) # ok
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18") )
				# Clean 
				rm(data2, data, tokeep, obis2, gbif2, obis, gbif) ; gc()
			
		} else {
			
				# Simply apply the same steps as above but on 'data'
				data$day[is.na(data$day)] <- 999
				# Create an ID based on: xbin, ybin, day, month, year, species name
				data$id <- paste(data$xbin, data$ybin, data$depth, data$day, data$month, data$year, data$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				data2 <- data[!duplicated(data$id),]
				# dim(data2) ; dim(data)
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18"))
				# Clean 
				rm(data2, data) ; gc()
			
		} # eo if else loop
		
	
} # eo for loop


### Check results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
# dir()
for(f in dir()) {
		# Useless message
		message(paste(str_replace(f, "_06_06_18.Rdata", ""), sep = ""))
		d <- get(load(f))
		message(paste(nrow(d), sep = ""))
		message(paste(length(unique(d$species)), sep = ""))
		message(paste("", sep = ""))
		rm(d)
} # eo for loop





### ----------------------------------------------------------------------------------------------------------------------------

##### 3°) v7-5.2v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.1/")
dir()

# For testing:
# f <- dir()[2]
for(f in dir()) {
	
		# Useless message()
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.1/")
		message(paste("Doing ", f, sep = ""))
		data <- get(load(f))
		
		# Converting some stuff if needed
		data$x <- as.numeric(as.character(data$x))
		data$y <- as.numeric(as.character(data$y))
		data$depth <- as.numeric(as.character(data$depth))
		if( f == "Thecosomata_MAREDAT_05_06_18.Rdata" ) {
				colnames(data)[7] <- "date"
		} else {
				data$date <- as.Date(data$date)
		} # eo if else loop
		
		data$day <- as.numeric(as.character(data$day))
		data$month <- as.numeric(as.character(data$month))
		data$year <- as.numeric(as.character(data$year))
		
		# Round the coordinates to the 1/4° degree
		data$xbin <- as.factor(round(data$x, 3))
		data$ybin <- as.factor(round(data$y, 3))
		# head(data)
		# tail(data)
		
		# If length(unique(data$source)) > 1 -> you have OBIS, GBIF
		if( length(unique(data$source)) > 1 ) {
			
				# Separate OBIS from GBIF again
				obis <- data[data$source == "OBIS",]
				gbif <- data[data$source == "GBIF",]
				# Convert NA days to 999
				gbif$day[is.na(gbif$day)] <- 999
				obis$day[is.na(obis$day)] <- 999
				# Create an ID based on: xbin, ybin, depth, day, month, year, species name
				obis$id <- paste(obis$xbin, obis$ybin, obis$depth, obis$day, obis$month, obis$year, obis$species, sep = "_") 
				gbif$id <- paste(gbif$xbin, gbif$ybin, gbif$depth, gbif$day, gbif$month, gbif$year, gbif$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				gbif2 <- gbif[!duplicated(gbif$id),]
				obis2 <- obis[!duplicated(obis$id),]
				# nrow(gbif2) ; nrow(gbif)
				# nrow(obis2) ; nrow(obis)
				# Now, compare the ids between OBIS & GBIF and remove those in GBIF
				diff <- setdiff(gbif2$id, obis2$id)
				# length(diff)
				# (length(diff) / nrow(gbif2))*100 # % of non duplicates
				# Use 'diff' to identify the obs of gbif2 that are to be added to obis2
				tokeep <- gbif2[gbif2$id %in% diff,]
				# rbind
				data2 <- rbind(obis2, tokeep)
				# nrow(data2) / nrow(data)
				# length(unique(data2$id)) # ok
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.1/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18") )
				# Clean 
				rm(data2, data, tokeep, obis2, gbif2, obis, gbif) ; gc()
			
		} else {
			
				# Simply apply the same steps as above but on 'data'
				data$day[is.na(data$day)] <- 999
				# Create an ID based on: xbin, ybin, day, month, year, species name
				data$id <- paste(data$xbin, data$ybin, data$depth, data$day, data$month, data$year, data$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				data2 <- data[!duplicated(data$id),]
				# dim(data2) ; dim(data)
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.1/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18"))
				# Clean 
				rm(data2, data) ; gc()
			
		} # eo if else loop
		
	
} # eo for loop



### Check results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.1/")
# dir()
for(f in dir()) {
		# Useless message
		message(paste(str_replace(f, "_06_06_18.Rdata", ""), sep = ""))
		d <- get(load(f))
		message(paste(nrow(d), sep = ""))
		message(paste(length(unique(d$species)), sep = ""))
		message(paste("", sep = ""))
		rm(d)
} # eo for loop




### ----------------------------------------------------------------------------------------------------------------------------

##### 4°) v7-5.2v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
dir()

# For testing:
# f <- dir()[2]
for(f in dir()) {
	
		# Useless message()
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
		message(paste("Doing ", f, sep = ""))
		data <- get(load(f))
		
		# Converting some stuff if needed
		data$x <- as.numeric(as.character(data$x))
		data$y <- as.numeric(as.character(data$y))
		data$depth <- as.numeric(as.character(data$depth))
		if( f == "Thecosomata_MAREDAT_05_06_18.Rdata" ) {
				colnames(data)[7] <- "date"
		} else {
				data$date <- as.Date(data$date)
		} # eo if else loop
		
		data$day <- as.numeric(as.character(data$day))
		data$month <- as.numeric(as.character(data$month))
		data$year <- as.numeric(as.character(data$year))
		
		# Round the coordinates to the 1/4° degree
		data$xbin <- as.factor(round(data$x, 3))
		data$ybin <- as.factor(round(data$y, 3))
		# head(data)
		# tail(data)
		
		# If length(unique(data$source)) > 1 -> you have OBIS, GBIF
		if( length(unique(data$source)) > 1 ) {
			
				# Separate OBIS from GBIF again
				obis <- data[data$source == "OBIS",]
				gbif <- data[data$source == "GBIF",]
				# Convert NA days to 999
				gbif$day[is.na(gbif$day)] <- 999
				obis$day[is.na(obis$day)] <- 999
				# Create an ID based on: xbin, ybin, depth, day, month, year, species name
				obis$id <- paste(obis$xbin, obis$ybin, obis$depth, obis$day, obis$month, obis$year, obis$species, sep = "_") 
				gbif$id <- paste(gbif$xbin, gbif$ybin, gbif$depth, gbif$day, gbif$month, gbif$year, gbif$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				gbif2 <- gbif[!duplicated(gbif$id),]
				obis2 <- obis[!duplicated(obis$id),]
				# nrow(gbif2) ; nrow(gbif)
				# nrow(obis2) ; nrow(obis)
				# Now, compare the ids between OBIS & GBIF and remove those in GBIF
				diff <- setdiff(gbif2$id, obis2$id)
				# length(diff)
				# (length(diff) / nrow(gbif2))*100 # % of non duplicates
				# Use 'diff' to identify the obs of gbif2 that are to be added to obis2
				tokeep <- gbif2[gbif2$id %in% diff,]
				# rbind
				data2 <- rbind(obis2, tokeep)
				# nrow(data2) / nrow(data)
				# length(unique(data2$id)) # ok
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.2/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18") )
				# Clean 
				rm(data2, data, tokeep, obis2, gbif2, obis, gbif) ; gc()
			
		} else {
			
				# Simply apply the same steps as above but on 'data'
				data$day[is.na(data$day)] <- 999
				# Create an ID based on: xbin, ybin, day, month, year, species name
				data$id <- paste(data$xbin, data$ybin, data$depth, data$day, data$month, data$year, data$species, sep = "_") 
				# Remove ducplicates from both obis and gbif
				data2 <- data[!duplicated(data$id),]
				# dim(data2) ; dim(data)
				# Got to proper dir and save
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.2/")
				save(data2, file = str_replace(f, "05_06_18", "06_06_18"))
				# Clean 
				rm(data2, data) ; gc()
			
		} # eo if else loop
		
	
} # eo for loop


### Check results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.2v3.2/")
# dir()
for(f in dir()) {
		# Useless message
		message(paste(str_replace(f, "_06_06_18.Rdata", ""), sep = ""))
		d <- get(load(f))
		message(paste(nrow(d), sep = ""))
		message(paste(length(unique(d$species)), sep = ""))
		message(paste("", sep = ""))
		rm(d)
} # eo for loop



### ----------------------------------------------------------------------------------------------------------------------------

##### 06/06/18: Combine all datasets and map sampling effort
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
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data[,c("id","x","y","species","family","order")])
}
) # eo lapply
table <- do.call(rbind, res)
dim(table)
str(table)
length(unique(table$species)) 

### Counts per sp/ family/ orders
counts <- data.frame(table %>% count(species) )
head(counts)
counts[order(counts$n, decreasing = T),]
# Assess richness of sampling effort thresholds
dim(counts[counts$n >= 1000,]) # 135	; 168 	; 99 	; 133
dim(counts[counts$n >= 100,])  # 517 	; 608 	; 314 	; 408
dim(counts[counts$n >= 50,])   # 689 	; 786 	; 447 	; 529
dim(counts[counts$n >= 30,])   # 926

# counts <- data.frame(table %>% count(order) )
# counts[order(counts$n, decreasing = T),]

### Bin x and y in 1/4°, use dplyr() to compute sampling effort and then map with 'coast'
colnames(table)
table$xbin <- round(table$x, digits = 0.1)
table$ybin <- round(table$y, digits = 0.1)
summary(table)
table[1:1000,c("x","y","xbin","ybin")]
table$id2 <- paste(table$xbin, table$ybin, sep = "_")

### dplyr
require("dplyr")
require("viridis")
table2 <- data.frame(table %>% 
			group_by(id2) %>% 
			summarise(x = unique(xbin), y = unique(ybin), n = n(), richness = length(unique(species)))
) # eo ddf
dim(table2)
head(table2)

map_sampling <- ggplot() + geom_tile(aes(x = x, y = y, fill = log(n)), data = table2) + 
		scale_fill_distiller(name = "Sampling effort\nlog(n occurrences)", palette = "YlGnBu") + coast + coord_quickmap()

map_rich <- ggplot() + geom_tile(aes(x = x, y = y, fill = richness), data = table2) + 
		scale_fill_distiller(name = "Species richness", palette = "YlGnBu") + coast + coord_quickmap()

# Save
ggsave(plot = map_sampling, filename = "map_effort_v8v5.2v3.2.pdf", dpi = 300, width = 12, height = 9)
ggsave(plot = map_rich, filename = "map_richness_v8v5.2v3.2.pdf", dpi = 300, width = 12, height = 9)


# OKAY



### ----------------------------------------------------------------------------------------------------------------------------

##### 25/07/18: Map sampling effort of each major zooplankton group, to help choose main background selection groups

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
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			return(data[,c("id","x","y","species","family","order","class","phylum")])
}
) # eo lapply
table <- do.call(rbind, res)
dim(table) # 766'033      6
str(table)
colnames(table)
unique(table$order)
unique(table$class)
unique(table$phylum)

### Map sampling effort for unique phylum !
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/")
for(p in unique(table$phylum)) {
	
			message(paste("Mapping ", p, sep = ""))
			data <- table[table$phylum == p,]
			
			# provide 1° bins & ID
			data$xbin <- round(data$x, digits = 0.1)
			data$ybin <- round(data$y, digits = 0.1)
			data$id2 <- paste(data$xbin, data$ybin, sep = "_")
			
			# compute sampling effort 
			table2 <- data.frame(data %>% 
						group_by(id2) %>% 
						summarise(x = unique(xbin), y = unique(ybin), n = n(), richness = length(unique(species)))
			) # eo ddf
			
			# Make map
			map_sampling <- ggplot() + geom_tile(aes(x = x, y = y, fill = log(n)), data = table2) + 
								scale_fill_distiller(name = "Sampling effort\nlog(n occurrences)", palette = "YlGnBu") + 
								xlab("Longitude") + ylab("Latitude") + coast + coord_quickmap() + ggtitle(p)
								
			# Save
			ggsave(plot = map_sampling, filename = paste("map_sampling_effort_", p, "_v8-5.1v3.2.pdf", sep = ""), dpi = 300, width = 10, height = 8)
			
			# Make space
			rm(data, table2, map_sampling)
			gc()
	
} # eo for loop mapping















