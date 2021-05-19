

setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v1/")
dir()

### For OBIS: 
setwd(paste(getwd(),"/","Copepoda_families_GBIF_17_04_18","/", sep = ""))
dir()

# read data
d <- read.csv("Podonidae_OBIS_18_04_18.csv", h = T, sep = ",")
dim(d)
head(d) 

unique(d$originalscientificname)
unique(d$species)


### For GBIF: 
d <- read.csv("Mysidae_GBIF_13_04_18.csv", sep = "\t")
dim(d)
head(d)
unique(d$species)
#unique(d$scientificname)
unique(d$basisofrecord)
nrow(d[which(d$basisofrecord == "FOSSIL_SPECIMEN"),])


# For GBIF's Hyperiids or Copepods
files <- dir()[grep(".csv", dir())]
#data <- get(load("Calocalanus indicus_OBIS_16_04_18.RData"))
#colnames(data)

### Dimensions for OBIS copepod data
#dims <- c("id","decimalLongitude","decimalLatitude","eventDate","datasetName","phylum","order","family",
		  #"genus","scientificName","species","originalScientificName","yearcollected" )
		  
res <- lapply(files, function(f) {
			message(paste(f, sep = ""))
			#data <- get(load(f))
			data <- read.csv(f, sep = "\t")
			return(data)
}) # eo lapply
# Check
data <- do.call(rbind, res)
rm(res) ; gc()
dim(data)
unique(data$species)


# Thecosomata OBIS : 100215 ; 55 species
# Thecosomata GBIF : 93805 ; 269 species (with fossils, but already have the nomenclature prepared to remove those)
# Gymnosomata OBIS : 15323 ; 16
# Gymnosomata GBIF : 11671 ; 37
# Euphausiidae OBIS : 179844 ; 82
# Euphausiidae GBIF : 121928 ; 85
# Ctenophora OBIS : 31118 ; 54
# Ctenophora GBIF : 40841 ; 65
# Scyphozoa OBIS : 41173 ; 109
# Scyphozoa GBIF : 64991 ; 223
# Hydrozoa OBIS : 67789 ; 249
# Hydrozoa GBIF : 24101 ; 292 
# Cubozoa OBIS : 315 ; 22
# Cubozoa GBIF : 1278 ; 40
# Copelata OBIS : 124549 ; 45
# Copelata GBIF : 34729 ; 52
# Sagittoidea OBIS : 210438 ; 70
# Sagittoidea GBIF : 26839 ; 87
# Thaliacea OBIS : 55282 ; 67
# Thaliacea GBIF : 38755 ; 73
# Myodocopia OBIS : 22850 ; 494
# Myodocopia GBIF : 16538 ; 753
# Mysidae OBIS : 68762 ; 590
# Mysidae GBIF : 34015 ; 900
# Hyperiidae OBIS : 104557 ; 244
# Hyperiidae GBIF : 31393 ; 191
# Tomopteridae : 18160 ; 19
# Alciopidae : 2285 ; 32
# Lopadorrhynchidae : 1456 ; 13
# Typhloscolecidae : 1328 ; 7
# Penilia avisrostris : 5334 ; 1
# Podonidae : 71205 ; 10


# -------------------------------------------------------------------------------------------------------------------


# Copepoda OBIS : 1865514 ; 669 species with at least 20 occurrences
# Copepoda GBIF : 1292633 ; 2520 because of many species with a few occ and synonyms...


1865514+1292633+100215+93805+15323+11671+179844+121928+31118+40841+41173+64991+67789+24101+315+1278+124549+34729+210438+26839+
55282+38755+22850+16538+68762+34015+104557+31393+18160+2285+1456+1328+5334+71205


# -------------------------------------------------------------------------------------------------------------------


##### 20/04/2018: Remove benthic and meroplanktonic species within the Hydrozoa, using the data from Gibbons et al. (2010) - J Biogeo.

### 0) Remove data with no species names

### 1) Remove all Anthoathecata but 2 families: Margelospidae & Porpitidae (Vellela)

### 2) Remove all Leptothecata and Limnomedusae

### 3) Keep all Siphonophora but 1 family: the Rhodaliidae 

### 4) Keep all Narcomedusae and Trachymedusae but 1 family: the Ptychogastriidae

### 5) Remove Actinulidae


### A) OBIS data
d <- read.csv("Hydrozoa_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 384 676 observations
# head(d) 

dd <- get(load("Hydrozoa_OBIS_20_04_18.Rdata"))
summary(dd$family)
summary(dd$order)

colnames(d)
unique(d$originalscientificname)
unique(d$species)

summary(d$family)
summary(d$order)

library("stringr")
d$species <- str_replace_all(as.character(d$species), " ", "_")
nrow(d[d$species == "",]) # 206 730
summary(d[d$species == "",])

# 0)
library("tidyr")
d <- d[- which(d$species == ""),]
dim(d) # 177 946

# 1)
unique(d[d$order == "Anthoathecata","family"]) # Remove these but Margelopsidae & Porpitidae
fam_2_rm <- unique(d[d$order == "Anthoathecata","family"])
fam_2_rm <- fam_2_rm[-c(16,21)]
d1 <- d[-which(d$order == "Anthoathecata" & d$family %in% fam_2_rm),]
summary(d1[d1$order == "Anthoathecata","family"]) # ok
dim(d1) # 146 601


# 2) 
d2 <- d1[-which(d1$order == "Leptothecata"),]
dim(d2) # 68 297
unique(d2$family) # No Limnomedusae in the data

# 3) 
d3 <- d2[-which(d2$family == "Rhodaliidae"),]
dim(d3) # 68 247 ok

# 4) 
unique(d3$family)
d4 <- d3[-which(d3$family == "Ptychogastriidae"),]
dim(d4) # 68 203
summary(d4$order)

# 5) 
d5 <- d4[-which(d4$order == "Actinulida"),]
dim(d5) # 68 179
summary(d5$order)
summary(d5$family)

nrow(d5) # 68179
unique(d5$species)

save(d5, file = "Hydrozoa_OBIS_20_04_18.Rdata")

dd <- dd[-which(dd$order == "Limnomedusae"),]
dim(dd)
save(dd, file = "Hydrozoa_OBIS_20_04_18.Rdata")
unique(dd$species)

### B) GBIF data
d <- read.csv("Hydrozoa_GBIF_13_04_18.csv", sep = "\t")
dim(d) # 190 758

colnames(d)
unique(d$species)
d$species <- str_replace_all(as.character(d$species), " ", "_")
d <- d[- which(d$species == ""),]
dim(d) # 134 951

unique(d$order)
unique(d$family)

# 1) 
unique(d[d$order == "Anthoathecata","family"]) # Remove these but Margelopsidae & Porpitidae
fam_2_rm <- unique(d[d$order == "Anthoathecata","family"])
fam_2_rm <- fam_2_rm[-c(14,32)]
d1 <- d[-which(d$order == "Anthoathecata" & d$family %in% fam_2_rm),]
summary(d1[d1$order == "Anthoathecata","family"]) # ok
dim(d1) # 100 923

# 2) 
d2 <- d1[-which(d1$order == "Leptothecata"),]
dim(d2) # 25 250
unique(d2$order) 
d2 <- d2[-which(d2$order == "Limnomedusae"),]
dim(d2) # 24 274

# 3) 
unique(d2$family)
d3 <- d2[-which(d2$family == "Rhodaliidae"),]
dim(d3) # 24 142

# 4) 
unique(d3$family)
d4 <- d3[-which(d3$family == "Ptychogastriidae"),]
dim(d4) # 24 122
# summary(d4$order)

# 5) 
d5 <- d4[-which(d4$order == "Actinulida"),]
dim(d5) # 24 101
summary(d5$order)
summary(d5$family)
unique(d5$species)

nrow(d5) # 24101
unique(d5$species) # 292 spp

save(d5, file = "Hydrozoa_GBIF_20_04_18.Rdata")



# -------------------------------------------------------------------------------------------------------------------


##### 20/04/2018: Check why there are so many species in the Thecosomata data from GBIF

d <- read.csv("Thecosomata_GBIF_13_04_18.csv", sep = "\t")
dim(d)
head(d)
colnames(d)

d <- d[-which(d$species == ""),]

unique(d$species)
unique(d$family)
unique(d$genus)

library("dplyr")

dd <- data.frame(d %>%
	group_by(species) %>%
	summarise(count = n() )
)

dd

write.table(dd, "Thecosomata_species_names_GBIF.txt")


##### 23/04/2018: Compare with the species names in OBIS - Correct species names

d <- read.csv("Thecosomata_OBIS_12_04_18.csv", h = T, sep = ",")
dim(d) # 100 215 observations
colnames(d)

length(unique(d$species)) # 56
unique(d$species)
library("dplyr")
dd <- data.frame(d %>%
	group_by(species) %>%
	summarise(count = n() )
)

dd

write.table(dd, "Thecosomata_species_names_OBIS.txt")
























