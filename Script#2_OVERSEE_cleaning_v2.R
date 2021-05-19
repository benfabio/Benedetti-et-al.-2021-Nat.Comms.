
##### 23/04/2017: R Script to start cleaning the global zooplankton occurrences © Fabio Benedetti, ETH Zürich, IBP, UP Group

### For each large zooplankton group (20x2 datasets), aims to:
#	- Remove occurrences without coords
#	- Remove occurrences without date
#	- Remove occurrences with date < 1800
#	- Remove occurrences without depth data
#	- Remove occurrences without species level identification
#	- Note number of remaining occurrences each time ; report to excel sheet

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 23/04/2017

library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")

### ----------------------------------------------------------------------------------------------------------------------------

### First wd
wd <- getwd()
#  [1] "Alciopidae_OBIS_12_04_18"          "Copelata_GBIF_13_04_18.csv"
#  [3] "Copelata_OBIS_12_04_18"            "Copepoda_families_GBIF_17_04_18"
#  [5] "Copepoda_OBIS_16_04_18"            "Ctenophora_GBIF_13_04_18.csv"
#  [7] "Ctenophora_OBIS_12_04_18"          "Cubozoa_GBIF_13_04_18.csv"
#  [9] "Cubozoa_OBIS_12_04_18"             "Euphausiidae_GBIF_13_04_18.csv"
# [11] "Euphausiidae_OBIS_12_04_2018"      "Gammaridae_OBIS_12_04_18"
# [13] "Gymnosomata_GBIF_13_04_18.csv"     "Gymnosomata_OBIS_12_04_18"
# [15] "Hydrozoa_GBIF_20_04_18.Rdata"      "Hydrozoa_OBIS_12_04_18"
# [17] "Hyperiid_families_GBIF_13_04_18"   "Hyperiidae_OBIS_12_04_18"
# [19] "Lopadorrhynchidae_OBIS_12_04_18"   "Myodocopida_OBIS_12_04_18"
# [21] "Myodocopina_GBIF_13_04_18.csv"     "Mysidae_GBIF_13_04_18.csv"
# [23] "Mysidae_OBIS_12_04_18"             "Ostracoda_OBIS_12_04_18"
# [25] "Penilia_avirostris_OBIS_18_04_18"  "planktonic_copepoda_species.RData"
# [27] "Podonidae_OBIS_18_04_18"           "Sagittoidea_GBIF_13_04_18.csv"
# [29] "Sagittoidea_OBIS_12_04_18"         "Scyphozoa_GBIF_13_04_18.csv"
# [31] "Scyphozoa_OBIS_12_04_18"           "Thaliacea_GBIF_13_04_18.csv"
# [33] "Thaliacea_OBIS_12_04_18"           "Thecosomata_GBIF_13_04_18.csv"
# [35] "Thecosomata_OBIS_12_04_18"         "Tomopteridae_OBIS_12_04_18"
# [37] "Typhloscolecidae_OBIS_12_04_18"


### Copepoda
# OBIS
setwd(paste(wd,"/","Copepoda_OBIS_16_04_18","/", sep = ""))
# dir()
files <- dir()
# For each file, load, apply data cleaning and print in v2 directory
# For tetsing:
#f <- "Calanus agulhensis_OBIS_16_04_18.RData"
f <- "Acartia (Acartiura) hongi_OBIS_16_04_18.RData"
for(f in files) {
		# Message
		message(paste("Cleaning ", f, sep = ""))
		# Load
		d <- get(load(f))
		# dim(d)
		# colnames(d)
		# summary(d)
		
		# Clean
		d1 <- d %>% drop_na(decimalLongitude, decimalLatitude)
		# dim(d1)
		d2 <- d1 %>% drop_na(eventDate)
		# dim(d2)
		d3 <- subset(d2, yearcollected > 1800)
		# dim(d3)
		# Some times there are no depth col at all
		if( "depth" %in% colnames(d3) ) {
				
				d4 <- d3 %>% drop_na(depth)
				# dim(d4)
				d5 <- d4 %>% drop_na(species)
				# Save in v2 dir
				setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
				save(d5, file = paste(unique(d5$species),"OBIS_23_04_18.Rdata", sep = "_"))
				setwd(paste(wd,"/","Copepoda_OBIS_16_04_18","/", sep = "")) # got back to initial dir
				rm(d, d1, d2, d3, d4, d5)
				
		} else {
			
			setwd(paste(wd,"/","Copepoda_OBIS_16_04_18","/", sep = "")) # got back to initial dir
			rm(d, d1, d2, d3)
			
		}

} # eo for loop

# Go to v2 dir and see how many occ you've cleaned
files <- dir()
dims <- c("id","decimalLongitude","decimalLatitude","eventDate","datasetName","phylum","order","family",
		  "genus","scientificName","species","originalScientificName","yearcollected" )
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- get(load(f))
			#data <- read.csv(f, sep = "\t")
			return(data[,dims])
}) # eo lapply
# Check
data <- do.call(rbind, res)
rm(res) ; gc()
dim(data)
unique(data$species)


# GBIF
setwd(paste(wd,"/","Copepoda_families_GBIF_17_04_18","/", sep = ""))
files <- dir()
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			#data <- get(load(f))
			data <- read.csv(f, sep = "\t")
			return(data) #[,dims])
}) # eo lapply
# Check
d <- do.call(rbind, res)
rm(res) ; gc()
dim(d)
colnames(d)
unique(d$species)

# Clean
d1 <- d %>% drop_na(decimallongitude, decimallatitude)
dim(d1) # 1269366
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 1250033
d3 <- subset(d2, year > 1800)
dim(d3) # 1250033
d4 <- d3 %>% drop_na(depth)
dim(d4) # 1180194
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 671913
unique(d5$species) # Lots of cleaning ahead for this dataset ^^'
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Copepoda_GBIF_23_04_18.Rdata" )


### Thecosomata
# OBIS
setwd(paste(wd,"/","Thecosomata_OBIS_12_04_18","/", sep = ""))
# dir()
d <- read.csv("Thecosomata_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 100215
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 100215
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 94611
d3 <- subset(d2, year > 1800)
dim(d3) # 94607
d4 <- d3 %>% drop_na(depth)
dim(d4) # 90779
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 59879
unique(d5$species) 
### Check the species names associated to negative depths
summary(d5[which(d5$depth < 0),])
head(d5[which(d5$depth < 0),])
unique(d5[which(d5$depth < 0),"species"]) # Limacina helicina, so gut
d5$depth <- abs(d5$depth)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Thecosomata_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Thecosomata_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 93805
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 69824
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 62410
d3 <- subset(d2, year > 1800)
dim(d3) # 62410
d4 <- d3 %>% drop_na(depth)
dim(d4) # 59025
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 54419
unique(d5$species) 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Thecosomata_GBIF_23_04_18.Rdata" )


### Gymnosomata
# OBIS
setwd(paste(wd,"/","Gymnosomata_OBIS_12_04_18","/", sep = ""))
# dir()
d <- read.csv("Gymnosomata_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 15323
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 15323
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 13730
d3 <- subset(d2, year > 1800)
dim(d3) # 13730
d4 <- d3 %>% drop_na(depth)
dim(d4) # 13521
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 9722
unique(d5$species) 
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Gymnosomata_OBIS_23_04_18.Rdata" )

# GBIF
d <- read.csv("Gymnosomata_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 11671
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 10080
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 9168
d3 <- subset(d2, year > 1800)
dim(d3) # 9168
d4 <- d3 %>% drop_na(depth)
dim(d4) # 8668
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 7782
unique(d5$species) 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Gymnosomata_GBIF_23_04_18.Rdata" )


### Euphausiidae
# OBIS
setwd(paste(wd,"/","Euphausiidae_OBIS_12_04_2018","/", sep = ""))
# dir()
d <- read.csv("Euphausiidae_OBIS_12_04_2018.csv", h = T, sep = ",")
dim(d) # 179844
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 179844
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 176997
d3 <- subset(d2, year > 1800)
dim(d3) # 176997
d4 <- d3 %>% drop_na(depth)
dim(d4) # 94190
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 58752
unique(d5$species) 
summary(d5$depth)
d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Euphausiidae_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Euphausiidae_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 121928
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 119500
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 50810
d3 <- subset(d2, year > 1800)
dim(d3) # 50810
d4 <- d3 %>% drop_na(depth)
dim(d4) # 34671
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 31746
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Euphausiidae_GBIF_23_04_18.Rdata" )


### Ctenophora
# OBIS
setwd(paste(wd,"/","Ctenophora_OBIS_12_04_18","/", sep = ""))
# dir()
d <- read.csv("Ctenophora_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 31118
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 31118
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 29305
d3 <- subset(d2, year > 1800)
dim(d3) # 29305
d4 <- d3 %>% drop_na(depth)
dim(d4) # 26570
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 9481
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Ctenophora_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Ctenophora_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 40841
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 40353
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 7525
d3 <- subset(d2, year > 1800)
dim(d3) # 7525
d4 <- d3 %>% drop_na(depth)
dim(d4) # 3975
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 2279
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Ctenophora_GBIF_23_04_18.Rdata" )


### Scyphozoa
# OBIS
setwd(paste(wd,"/","Scyphozoa_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Scyphozoa_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 41173
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 41173
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 38635
d3 <- subset(d2, year > 1800)
dim(d3) # 38631
d4 <- d3 %>% drop_na(depth)
dim(d4) # 25878
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 12368
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Scyphozoa_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Scyphozoa_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 64991
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 57882
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 53865
d3 <- subset(d2, year > 1800)
dim(d3) # 53865
d4 <- d3 %>% drop_na(depth)
dim(d4) # 18633
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 14427
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Scyphozoa_GBIF_23_04_18.Rdata" )


### Hydrozoa
# OBIS
setwd(paste(wd,"/","Hydrozoa_OBIS_12_04_18","/", sep = ""))
dir()
d <- get(load("Hydrozoa_OBIS_20_04_18.Rdata"))
dim(d) # 67789
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 67789
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 63671
d3 <- subset(d2, year > 1800)
dim(d3) # 63659
d4 <- d3 %>% drop_na(depth)
dim(d4) # 60116
unique(d4$species)
#d5 <- d4[-which(d4$species == ""),]
dim(d4) # 60116
unique(d4$species)
summary(d4$depth)
d4$depth <- abs(d4$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d4, file = "Hydrozoa_OBIS_23_04_18.Rdata" )


# GBIF
d <- get(load("Hydrozoa_GBIF_20_04_18.Rdata"))
dim(d) # 24101
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 24101
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 17727
d3 <- subset(d2, year > 1800)
dim(d3) # 17727
d4 <- d3 %>% drop_na(depth)
dim(d4) # 13232
unique(d4$species)
#d5 <- d4[-which(d4$species == ""),]
dim(d4) # 13232
summary(d4$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d4, file = "Hydrozoa_GBIF_23_04_18.Rdata" )


### Cubozoa
# OBIS
setwd(paste(wd,"/","Cubozoa_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Cubozoa_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 315
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 315
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 246
d3 <- subset(d2, year > 1800)
dim(d3) # 246
d4 <- d3 %>% drop_na(depth)
dim(d4) # 74
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 66
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Cubozoa_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Cubozoa_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 1278
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 829
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 685
d3 <- subset(d2, year > 1800)
dim(d3) # 685
d4 <- d3 %>% drop_na(depth)
dim(d4) # 356
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 340
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Cubozoa_GBIF_23_04_18.Rdata" )


### Copelata
# OBIS
setwd(paste(wd,"/","Copelata_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Copelata_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 124549
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 124549
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 122781
d3 <- subset(d2, year > 1800)
dim(d3) # 122781
d4 <- d3 %>% drop_na(depth)
dim(d4) # 117646
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 20396
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Copelata_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Copelata_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 34729
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 34332
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 33411
d3 <- subset(d2, year > 1800)
dim(d3) # 33411
d4 <- d3 %>% drop_na(depth)
dim(d4) # 31346
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 3267
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Copelata_GBIF_23_04_18.Rdata" )



### Sagittoidea
# OBIS
setwd(paste(wd,"/","Sagittoidea_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Sagittoidea_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 210438
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 210438
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 188155
d3 <- subset(d2, year > 1800)
dim(d3) # 188155
d4 <- d3 %>% drop_na(depth)
dim(d4) # 181397
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 115677
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Sagittoidea_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Sagittoidea_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 26839
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 25914
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 21793
d3 <- subset(d2, year > 1800)
dim(d3) # 21793
d4 <- d3 %>% drop_na(depth)
dim(d4) # 18974
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 17529
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Sagittoidea_GBIF_23_04_18.Rdata" )



### Thaliacea
# OBIS
setwd(paste(wd,"/","Thaliacea_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Thaliacea_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 55282
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 55282
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 51863
d3 <- subset(d2, year > 1800)
dim(d3) # 51847
d4 <- d3 %>% drop_na(depth)
dim(d4) # 49538
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 6306
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Thaliacea_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Thaliacea_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 38755
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 37144
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 27937
d3 <- subset(d2, year > 1800)
dim(d3) # 27937
d4 <- d3 %>% drop_na(depth)
dim(d4) # 26482
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 9772
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Thaliacea_GBIF_23_04_18.Rdata" )


### Myodocopina
# OBIS
setwd(paste(wd,"/","Myodocopida_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Myodocopida_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 22850
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 22850
d2 <- d1[-which(d1$eventDate == ""),]
dim(d2) # 17067
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 17067
summary(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 7523
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 6502
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Myodocopina_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Myodocopina_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 16538
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 13427
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 11085
d3 <- subset(d2, year > 1800)
dim(d3) # 11085
d4 <- d3 %>% drop_na(depth)
dim(d4) # 7271
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 4892
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Myodocopina_GBIF_23_04_18.Rdata" )



### Mysidae
# OBIS
setwd(paste(wd,"/","Mysidae_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Mysidae_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 68762
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 68762
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 60901
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 60901
summary(d3$depth) # -9999.00     9.80    21.55   247.66   185.58 12500.00    30128 
summary(d3[d3$depth < 0,"depth"]) # -9999.00   -21.84   -15.43 -1260.03   -12.72    -7.94    30128 
### CONVERT -9999.00 to NAs
d3[,"depth"][d3[,"depth"] == -9999] <- NA
d3$depth <- abs(d3$depth)
summary(d3$depth) #  0.00    10.00    21.84   260.61   188.00 12500.00    30166 
d4 <- d3 %>% drop_na(depth)
dim(d4) # 30735
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 21794
unique(d5$species)
summary(d5$depth)
#d5$depth <- abs(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Mysidae_OBIS_23_04_18.Rdata" )


# GBIF
d <- read.csv("Mysidae_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 34015
colnames(d)
d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 25200
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 23164
d3 <- subset(d2, year > 1800)
dim(d3) # 23164
#summary(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 7270
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 6170
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Mysidae_GBIF_23_04_18.Rdata" )



### Hyperiidea
# OBIS
setwd(paste(wd,"/","Hyperiidae_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Hyperiidae_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 104557
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 104557
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 102459
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 102459
d3$depth <- abs(d3$depth)
summary(d3$depth) #  0.00    7.50    7.50   74.82   27.50 8000.00    4129 
d4 <- d3 %>% drop_na(depth)
dim(d4) # 98330
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 24031
unique(d5$species)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Hyperiidea_OBIS_23_04_18.Rdata" )


# GBIF
setwd(paste(wd,"/","Hyperiid_families_GBIF_13_04_18","/", sep = ""))
dir()
files <- dir()[grep(".csv",dir())]
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			data <- read.csv(f, sep = "\t")
			return(data) 
}) # eo lapply
# Check
d <- do.call(rbind, res)
rm(res) ; gc()
dim(d) # 31393

d1 <- d %>% drop_na(decimallatitude, decimallongitude)
dim(d1) # 28552
d2 <- d1 %>% drop_na(month, year)
dim(d2) # 21592
d3 <- subset(d2, year > 1800)
dim(d3) # 21592 
# summary(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 18435
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 14484
unique(d5$species) 
summary(d5$depth)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Hyperiidea_GBIF_23_04_18.Rdata" )


### Annelida 
# Tomopteridae
setwd(paste(wd,"/","Tomopteridae_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Tomopteridae_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 18160
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 18160
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 16541
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 16541
summary(d3$depth) # 
#d3$depth <- abs(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 15788
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 2463
unique(d5$species)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Tomopteridae_OBIS_23_04_18.Rdata" )


# Alciopidae
setwd(paste(wd,"/","Alciopidae_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Alciopidae_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 2285
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 2285
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 2218
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 2218
summary(d3$depth) # 
#d3$depth <- abs(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 1321
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 1125
unique(d5$species)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Alciopidae_OBIS_23_04_18.Rdata" )


# Lopadorrhynchidae
setwd(paste(wd,"/","Lopadorrhynchidae_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Lopadorrhynchidae_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 
summary(d3$depth) 
#d3$depth <- abs(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 
unique(d5$species)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Lopadorrhynchidae_OBIS_23_04_18.Rdata" )


# Typhloscolecidae
setwd(paste(wd,"/","Typhloscolecidae_OBIS_12_04_18","/", sep = ""))
dir()
d <- read.csv("Typhloscolecidae_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 1328
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 1328
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 1255
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 1255
summary(d3$depth) 
#d3$depth <- abs(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 1106
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 949
unique(d5$species)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Typhloscolecidae_OBIS_23_04_18.Rdata" )


### Cladocera
# Penilia avirostris
setwd(paste(wd,"/","Penilia_avirostris_OBIS_18_04_18","/", sep = ""))
dir()
d <- read.csv("Penilia_avirostris_OBIS_18_04_18.csv", h = T, sep = ",")
dim(d) # 5334
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 5334
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 5316
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 5316
summary(d3$depth) 
#d3$depth <- abs(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 4197
str(d4$species)
#unique(d4$species)
summary(d4$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d4, file = "P.avirostris_OBIS_23_04_18.Rdata" )


# Podonidae
setwd(paste(wd,"/","Podonidae_OBIS_18_04_18","/", sep = ""))
dir()
d <- read.csv("Podonidae_OBIS_18_04_18.csv", h = T, sep = ",")
dim(d) # 71205
colnames(d)
head(d)
# Clean
d1 <- d %>% drop_na(decimalLatitude, decimalLongitude)
dim(d1) # 71205
d2 <- d1[- which(d1$eventDate == ""),]
dim(d2) # 71134
summary(d2)
d3 <- subset(d2, year > 1800)
dim(d3) # 71134
summary(d3$depth) 
#d3$depth <- abs(d3$depth)
d4 <- d3 %>% drop_na(depth)
dim(d4) # 67303
str(d4$species)
d5 <- d4[-which(d4$species == ""),]
dim(d5) # 24522
unique(d5$species)
summary(d5$depth)
# Save in v2 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
save(d5, file = "Podonidae_OBIS_23_04_18.Rdata" )


### For now, v2 data comprise 3'071'283 observations, which roughly corresponds to 63.7% of the v1 data













