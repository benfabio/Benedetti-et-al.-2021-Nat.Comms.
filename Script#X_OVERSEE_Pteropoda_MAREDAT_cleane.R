
##### 30/05/2018: R Script to 

### Aims to:
#	- Load the Pteropod data from MAREDAT (taxon-resolved abundances and biomass)
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

### Latest update: 31/05/2018

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


##### 1°) Load the .txt file in v1 directory and apply v2 criteria (basic stuff)

d <- read.table("Pteropoda_MAREDAT_30_05_2018.txt", h = TRUE, sep = "\t")
dim(d) # 8809   13
str(d)
summary(d)
unique(d$Taxon) # dirty dirty...
# Strange and very high latitudes

### Check coordinates
summary(d[,c("y","x")])
d[d$y > 90,] # Need to reverse y and x here
d[d$y > 90,"y"] <- 54.08
d[d$y > 90,"x"] <- 138.12
# OK.

### Check dates
summary(d[,c("Day","Month","Year")])
# No NAs, gut. No date before 1800 too

### Check depth
summary(d$Depth) # No NAs, gut

### Remove anything that is not species-resolved
labels <- unique(d$Taxon)
labels
tokeep <- labels[c(2:8,10:15,17,22:26,28:43,45:57,59,60,62:82,84,85,87:95,97:105)]
unique(tokeep) # 91 

dd <- d[d$Taxon %in% tokeep,]
dim(dd) # 4765 obs @ species level

### Need to clean those labels before though
dd$Taxon <- str_replace_all(as.character(dd$Taxon), " ", "_")
unique(dd$Taxon) 
#  [1] "Clio_pyramidata"								--> Clio_pyramidata, keep

#  [2] "L._helicina_antarctica"							--> Limacina_rangii
levels(dd$Taxon) <- c(levels(dd$Taxon),"Limacina_rangii")
dd[dd$Taxon == "L._helicina_antarctica","Taxon"] <- "Limacina_rangii"
#  [3] "Limacina_retroversa"							--> Limacina_retroversa, keep

#  [4] "Clione_limacina_veliger"						--> Clione_limacina
levels(dd$Taxon) <- c(levels(dd$Taxon),"Clione_limacina")
dd[dd$Taxon == "Clione_limacina_veliger","Taxon"] <- "Clione_limacina"
#  [5] "L._helicina_veliger"							--> Limacina_helicina
levels(dd$Taxon) <- c(levels(dd$Taxon),"Limacina_helicina")
dd[dd$Taxon == "L._helicina_veliger","Taxon"] <- "Limacina_helicina"
#  [6] "Clione_limacina"								--> Clione_limacina, keep

#  [7] "Limacina_helicina"								--> Limacina_helicina, keep

#  [8] "Creseis_acicula"								--> Creseis_clava
levels(dd$Taxon) <- c(levels(dd$Taxon),"Creseis_clava")
dd[dd$Taxon == "Creseis_acicula","Taxon"] <- "Creseis_clava"
#  [9] "Limacina_helicina_"								--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina_","Taxon"] <- "Limacina_helicina"
# [10] "Limacina_helicina_______"						--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina_______","Taxon"] <- "Limacina_helicina"
# [11] "Clio_piatkowskii"								--> Clio_piatkowskii, keep

# [12] "Spongiobranchaea_australis"						--> Spongiobranchaea_australis, keep	

# [13] "Clione_antarctica"								--> Clione_limacina
dd[dd$Taxon == "Clione_antarctica","Taxon"] <- "Clione_limacina"
# [14] "L._trochiformis"								--> Limacina_trochiformis
levels(dd$Taxon) <- c(levels(dd$Taxon),"Limacina_trochiformis")
dd[dd$Taxon == "L._trochiformis","Taxon"] <- "Limacina_trochiformis"
# [15] "L._inflata"										--> Heliconoides_inflatus
levels(dd$Taxon) <- c(levels(dd$Taxon),"Heliconoides_inflatus")
dd[dd$Taxon == "L._inflata","Taxon"] <- "Heliconoides_inflatus"
# [16] "S._subula"										--> Styliola_subula
levels(dd$Taxon) <- c(levels(dd$Taxon),"Styliola_subula")
dd[dd$Taxon == "S._subula","Taxon"] <- "Styliola_subula"
# [17] "L.inflata"										--> Heliconoides_inflatus
dd[dd$Taxon == "L.inflata","Taxon"] <- "Heliconoides_inflatus"
# [18] "L.trochiformis"									--> Limacina_trochiformis
dd[dd$Taxon == "L.trochiformis","Taxon"] <- "Limacina_trochiformis"
# [19] "S.subula"										--> Styliola_subula
dd[dd$Taxon == "S.subula","Taxon"] <- "Styliola_subula"
# [20] "Clione_limacina_(polytrochous_larvae)"			--> Clione_limacina
dd[dd$Taxon == "Clione_limacina_(polytrochous_larvae)","Taxon"] <- "Clione_limacina"
# [21] "Limacina_helicina_(adults)"						--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina_(adults)","Taxon"] <- "Limacina_helicina"
# [22] "Limacina_helicina_(juveniles_and_larvae)"		--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina_(juveniles_and_larvae)","Taxon"] <- "Limacina_helicina"
# [23] "Clione_limacina__(veliger)"						--> Clione_limacina
dd[dd$Taxon == "Clione_limacina__(veliger)","Taxon"] <- "Clione_limacina"
# [24] "Limacina_helicina_(larvae)"						--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina_(larvae)","Taxon"] <- "Limacina_helicina"
# [25] "Clione_limacina_(veliger)"						--> Clione_limacina
dd[dd$Taxon == "Clione_limacina_(veliger)","Taxon"] <- "Clione_limacina"
# [26] "Pteropoda_larvae"								--> remove !!
dim(dd) # 4765
dd <- dd[-which(dd$Taxon == "Pteropoda_larvae"),] # 4739 obs after

# [27] "Clio_cuspidata_juveniles"						--> Clio_cuspidata
levels(dd$Taxon) <- c(levels(dd$Taxon),"Clio_cuspidata")
dd[dd$Taxon == "Clio_cuspidata_juveniles","Taxon"] <- "Clio_cuspidata"

# [28] "Clio_pyramidata_f._lanceolata_juveniles"		--> Clio_pyramidata
levels(dd$Taxon) <- c(levels(dd$Taxon),"Clio_pyramidata")
dd[dd$Taxon == "Clio_pyramidata_f._lanceolata_juveniles","Taxon"] <- "Clio_pyramidata"

# [29] "Clio_cuspidata"									--> Clio_cuspidata, keep

# [30] "Cuvierina_columnella"							--> Cuvierina_columnella, keep

# [31] "Hyalocylis_striata"								--> Hyalocylis_striata, keep

# [32] "Cavolinia_globulosa"							--> Cavolinia_globulosa, keep

# [33] "Cavolinia_longirostris_f.._Longirostris"		--> Diacavolinia_longirostris
levels(dd$Taxon) <- c(levels(dd$Taxon),"Diacavolinia_longirostris")
dd[dd$Taxon == "Cavolinia_longirostris_f.._Longirostris","Taxon"] <- "Diacavolinia_longirostris"

# [34] "Clio_convexa"									--> Clio_convexa, keep

# [35] "Diacria_rampali"								--> Diacria_trispinosa
levels(dd$Taxon) <- c(levels(dd$Taxon),"Diacria_trispinosa")
dd[dd$Taxon == "Diacria_rampali","Taxon"] <- "Diacria_trispinosa"

# [36] "Clio_convexa_juveniles"							--> Clio_convexa
dd[dd$Taxon == "Clio_convexa_juveniles","Taxon"] <- "Clio_convexa"

# [37] "Diacria_quadridentata_group_juveniles"			--> Diacria_quadridentata	
levels(dd$Taxon) <- c(levels(dd$Taxon),"Diacria_quadridentata")
dd[dd$Taxon == "Diacria_quadridentata_group_juveniles","Taxon"] <- "Diacria_quadridentata"

# [38] "Cavolinia_longirostris_f._strangulata"			--> Diacavolinia_longirostris
dd[dd$Taxon == "Cavolinia_longirostris_f._strangulata","Taxon"] <- "Diacavolinia_longirostris"

# [39] "Cavolinia_unicata_unicata_f._pusilla"			--> Cavolinia_uncinata
levels(dd$Taxon) <- c(levels(dd$Taxon),"Cavolinia_uncinata")
dd[dd$Taxon == "Cavolinia_unicata_unicata_f._pusilla","Taxon"] <- "Cavolinia_uncinata"

# [40] "Diacria_quadridentata"							--> Diacria_quadridentata, keep

# [41] "Clio_pyramidata_f._lanceolata"					--> Clio_pyramidata
dd[dd$Taxon == "Clio_pyramidata_f._lanceolata","Taxon"] <- "Clio_pyramidata"

# [42] "Diacria_danae"									--> Diacria_danae, keep

# [43] "Creseis_virgula"								--> Creseis_virgula, keep

# [44] "Cavolinia_longirostris_f.._angulosa"			--> Diacavolinia_longirostris
dd[dd$Taxon == "Cavolinia_longirostris_f.._angulosa","Taxon"] <- "Diacavolinia_longirostris"

# [45] "Cavolinia_longirostris_f.._longirostris"		--> Diacavolinia_longirostris
dd[dd$Taxon == "Cavolinia_longirostris_f.._longirostris","Taxon"] <- "Diacavolinia_longirostris"

# [46] "Cavolinia_longirostris_f._angulosa"				--> Diacavolinia_longirostris
dd[dd$Taxon == "Cavolinia_longirostris_f._angulosa","Taxon"] <- "Diacavolinia_longirostris"

# [47] "Diacria_costata"								--> Diacria_costata, keep

# [48] "Limacina_inflata"								--> Heliconoides_inflatus
dd[dd$Taxon == "Limacina_inflata","Taxon"] <- "Heliconoides_inflatus"

# [49] "C._uncinata"									--> Cavolinia_uncinata
dd[dd$Taxon == "C._uncinata","Taxon"] <- "Cavolinia_uncinata"

# [50] "Cavolinia_longirostris"							--> Diacavolinia_longirostris
dd[dd$Taxon == "Cavolinia_longirostris","Taxon"] <- "Diacavolinia_longirostris"

# [51] "C._virgula_virgula"								--> Creseis_virgula
dd[dd$Taxon == "C._virgula_virgula","Taxon"] <- "Creseis_virgula"

# [52] "Cavolinia_l._juveniles"							--> Diacavolinia_longirostris
dd[dd$Taxon == "Cavolinia_l._juveniles","Taxon"] <- "Diacavolinia_longirostris"

# [53] "C._virgula_conica"								--> Creseis_virgula
dd[dd$Taxon == "C._virgula_conica","Taxon"] <- "Creseis_virgula"

# [54] "Diacria_q._juveniles"							--> Diacria_quadridentata
dd[dd$Taxon == "Diacria_q._juveniles","Taxon"] <- "Diacria_quadridentata"

# [55] "C._virgula_constricta"							--> Creseis_virgula
dd[dd$Taxon == "C._virgula_constricta","Taxon"] <- "Creseis_virgula"

# [56] "Hyalocylix_striata"								--> Hyalocylis_striata
levels(dd$Taxon) <- c(levels(dd$Taxon),"Hyalocylis_striata")
dd[dd$Taxon == "Hyalocylix_striata","Taxon"] <- "Hyalocylis_striata"

# [57] "Clio_pyramidata_juveniles"						--> Clio_pyramidata
dd[dd$Taxon == "Clio_pyramidata_juveniles","Taxon"] <- "Clio_pyramidata"

# [58] "C._unicata_total"								--> Cavolinia_uncinata
dd[dd$Taxon == "C._unicata_total","Taxon"] <- "Cavolinia_uncinata"

# [59] "Limacina_bulimoides"							--> Limacina_bulimoides, keep

# [60] "C._virgula_virgula_(total)"						--> Creseis_virgula
dd[dd$Taxon == "C._virgula_virgula_(total)","Taxon"] <- "Creseis_virgula"

# [61] "Cuvierina_columnella_(total)"					--> Cuvierina_columnella
dd[dd$Taxon == "Cuvierina_columnella_(total)","Taxon"] <- "Cuvierina_columnella"

# [62] "C._virgula_conica_(total)"						--> Creseis_virgula
dd[dd$Taxon == "C._virgula_conica_(total)","Taxon"] <- "Creseis_virgula"

# [63] "C._pyramidata_juveniles_(total)"				--> Clio_pyramidata
dd[dd$Taxon == "C._pyramidata_juveniles_(total)","Taxon"] <- "Clio_pyramidata"

# [64] "Hyalocylix_striata_(total)"						--> Hyalocylis_striata
dd[dd$Taxon == "Hyalocylix_striata_(total)","Taxon"] <- "Hyalocylis_striata"

# [65] "L._inflata_(toal)"								--> Heliconoides_inflatus
dd[dd$Taxon == "L._inflata_(toal)","Taxon"] <- "Heliconoides_inflatus"

# [66] "C._virgula_constricta_(total)"					--> Creseis_virgula
dd[dd$Taxon == "C._virgula_constricta_(total)","Taxon"] <- "Creseis_virgula"

# [67] "Diacria_quadridentata_8total)"					--> Diacria_quadridentata
dd[dd$Taxon == "Diacria_quadridentata_8total)","Taxon"] <- "Diacria_quadridentata"

# [68] "Cavolinia_(total)"								--> REMOVE !!!
dd <- dd[-which(dd$Taxon == "Cavolinia_(total)"),]
# dim(dd) # 4737

# [69] "Limacina_bulimoides_(total)"					--> Limacina_bulimoides
dd[dd$Taxon == "Limacina_bulimoides_(total)","Taxon"] <- "Limacina_bulimoides"

# [70] "L._trochiformis_(total)"						--> Limacina_trochiformis
dd[dd$Taxon == "L._trochiformis_(total)","Taxon"] <- "Limacina_trochiformis"

# [71] "Creseis_acicula_(total)"						--> Creseis_clava
dd[dd$Taxon == "Creseis_acicula_(total)","Taxon"] <- "Creseis_clava"

# [72] "Clione_Limacina"								--> Clione_limacina, keep

# [73] "Diacria_quadrientata"							--> Diacria_quadridentata
dd[dd$Taxon == "Diacria_quadrientata","Taxon"] <- "Diacria_quadridentata"

# [74] "Styliola_subula"								--> Styliola_subula, keep

# [75] "Spiratella_inflata"								--> Heliconoides_inflatus
dd[dd$Taxon == "Spiratella_inflata","Taxon"] <- "Heliconoides_inflatus"

# [76] "Spiratella_lesueuri"							--> Limacina_lesueurii
dd[dd$Taxon == "Spiratella_lesueuri","Taxon"] <- "Limacina_lesueurii"

# [77] "Spiratella_bulimoides"							--> Limacina_bulimoides	
dd[dd$Taxon == "Spiratella_bulimoides","Taxon"] <- "Limacina_bulimoides"

# [78] "Spiratella_trochiformis"						--> Limacina_trochiformis
dd[dd$Taxon == "Spiratella_trochiformis","Taxon"] <- "Limacina_trochiformis"

# [79] "Diacria_trispinosa"								--> Diacria_trispinosa, keep

# [80] "Cavolinia_inflexa"								--> Cavolinia_inflexa, keep	

# [81] "Creseis_virgula_virgula"						--> Creseis_virgula
dd[dd$Taxon == "Creseis_virgula_virgula","Taxon"] <- "Creseis_virgula"

# [82] "Creseis_virgula_conica"							--> Creseis_virgula
dd[dd$Taxon == "Creseis_virgula_conica","Taxon"] <- "Creseis_virgula"

# [83] "Limacina__helicina_______"						--> Limacina_helicina
dd[dd$Taxon == "Limacina__helicina_______","Taxon"] <- "Limacina_helicina"

# [84] "Limacina_helicina_______________"				--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina_______________","Taxon"] <- "Limacina_helicina"

# [85] "Limacina_helicina___"							--> Limacina_helicina
dd[dd$Taxon == "Limacina_helicina___","Taxon"] <- "Limacina_helicina"

# [86] "Cavolinia_immatures"							--> remove, fossil
dd <- dd[-which(dd$Taxon == "Cavolinia_immatures"),]
# dim(dd) # 4688

# [87] "L._lesueuri"									--> Limacina_lesueurii
dd[dd$Taxon == "L._lesueuri","Taxon"] <- "Limacina_lesueurii"

# [88] "L._retroversa"									--> Limacina_retroversa
dd[dd$Taxon == "L._retroversa","Taxon"] <- "Limacina_retroversa"

# [89] "L._bulimoides"									--> Limacina_bulimoides
dd[dd$Taxon == "L._bulimoides","Taxon"] <- "Limacina_bulimoides"

# [90] "Cresesis_acicula"								--> Creseis_clava
dd[dd$Taxon == "Cresesis_acicula","Taxon"] <- "Creseis_clava"

# [91] "Creseis_vir._Virgula"							--> Creseis_virgula
dd[dd$Taxon == "Creseis_vir._Virgula","Taxon"] <- "Creseis_virgula"

### OKAY, check results : 
unique(dd$Taxon) # 27 different species in total...oh wait, one correction left:
dd[dd$Taxon == "Clione_Limacina","Taxon"] <- "Clione_limacina"
unique(dd$Taxon) # 26 then
library("dplyr")
counts <- data.frame(dd %>%
			group_by(Taxon) %>%
			summarise(n = n() )
)
counts

### Check
dd[1:50,]
### 2 things to clean before you go any further: 
#	- are zero abund/ biomass to be removed because they can be considered as absences ? 
#	- there seems to be a lot of ducplicates, remove them by using an id that combines: species name, day, month, mon, lat, depth.

# Check if zero biomass necessarily corresponds to zero abund
summary(dd[dd$Abund == 0,"Biomass"])
summary(dd[dd$Biomass == 0,"Abund"])
# Ok yes.

ddd <- dd[which(dd$Abund > 0),]
nrow(ddd) # 3583 obs

### Create and id
ddd$id <- paste(ddd$x, ddd$y, ddd$depth, ddd$Day, ddd$Month, ddd$Year, ddd$Taxon, sep = "_")
unique(ddd$id) # 2882 --> not 3583 so there should be 701 duplicates
ddd[1:100,] # Yep, easy to spot. Remove

ddd <- ddd[!duplicated(ddd$id),]
nrow(ddd)
# 2882 ! Nice
head(ddd)
ddd[1:100,]

### Got to v2 directory and save ! 
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
dir()
save(ddd, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------


##### 2°) Apply both v3 criteria: v3.1 (land mask) and v3.2 (distance to coast)
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v2/")
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
dim(data)

### A) get the bathymetry raster
require("marmap")
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 15)
str(bathy)
# plot(bathy) # prints a plot on kryo-old...
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
data$bathymetry <- raster::extract(x = r, y = data[,c("x","y")], method = 'bilinear')
data2 <- data[-which(data$bathymetry >= -200),]
dim(data2) # 
summary(data2) # check NAs
data3 <- data2 %>% drop_na(bathymetry)
dim(data3) # 2198 obs
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.1/")
save(data3, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


### B) Get the dist2coast raster
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/dist2coast")
r <- get(load("dist2coast_raster_15min.Rdata"))
r #
### extract() the dist2coast data and remove what is within the first 25 kms
data$dist2coast <- raster::extract(x = r, y = data[,c("x","y")], method = 'bilinear')
data2 <- data[-which(data$dist2coast < 25),]
summary(data2$dist2coast)
# Remove NAs --> land cells
data3 <- data2 %>% drop_na(dist2coast)
dim(data3)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
save(data3, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------

##### 3°) Apply v4 to the two v3 datasets

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
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$Month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$Month == 1),c("x","y")], method = 'bilinear')
data[which(data$Month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$Month == 2),c("x","y")], method = 'bilinear')
data[which(data$Month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$Month == 3),c("x","y")], method = 'bilinear')
data[which(data$Month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$Month == 4),c("x","y")], method = 'bilinear')
data[which(data$Month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$Month == 5),c("x","y")], method = 'bilinear')
data[which(data$Month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$Month == 6),c("x","y")], method = 'bilinear')
data[which(data$Month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$Month == 7),c("x","y")], method = 'bilinear')
data[which(data$Month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$Month == 8),c("x","y")], method = 'bilinear')
data[which(data$Month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$Month == 9),c("x","y")], method = 'bilinear')
data[which(data$Month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$Month == 10),c("x","y")], method = 'bilinear')
data[which(data$Month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$Month == 11),c("x","y")], method = 'bilinear')
data[which(data$Month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$Month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# Remove obs with SSS <= 20
# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2)
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3) # 2190
summary(data3$salinityWOA13)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
save(data3, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


### B) Get v3.2 data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v3.2/")
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
# Create salinity vector
data$salinityWOA13 <- NA
# Fill it according to month
data[which(data$Month == 1),"salinityWOA13"] <- raster::extract(x = ras1, y = data[which(data$Month == 1),c("x","y")], method = 'bilinear')
data[which(data$Month == 2),"salinityWOA13"] <- raster::extract(x = ras2, y = data[which(data$Month == 2),c("x","y")], method = 'bilinear')
data[which(data$Month == 3),"salinityWOA13"] <- raster::extract(x = ras3, y = data[which(data$Month == 3),c("x","y")], method = 'bilinear')
data[which(data$Month == 4),"salinityWOA13"] <- raster::extract(x = ras4, y = data[which(data$Month == 4),c("x","y")], method = 'bilinear')
data[which(data$Month == 5),"salinityWOA13"] <- raster::extract(x = ras5, y = data[which(data$Month == 5),c("x","y")], method = 'bilinear')
data[which(data$Month == 6),"salinityWOA13"] <- raster::extract(x = ras6, y = data[which(data$Month == 6),c("x","y")], method = 'bilinear')
data[which(data$Month == 7),"salinityWOA13"] <- raster::extract(x = ras7, y = data[which(data$Month == 7),c("x","y")], method = 'bilinear')
data[which(data$Month == 8),"salinityWOA13"] <- raster::extract(x = ras8, y = data[which(data$Month == 8),c("x","y")], method = 'bilinear')
data[which(data$Month == 9),"salinityWOA13"] <- raster::extract(x = ras9, y = data[which(data$Month == 9),c("x","y")], method = 'bilinear')
data[which(data$Month == 10),"salinityWOA13"] <- raster::extract(x = ras10, y = data[which(data$Month == 10),c("x","y")], method = 'bilinear')
data[which(data$Month == 11),"salinityWOA13"] <- raster::extract(x = ras11, y = data[which(data$Month == 11),c("x","y")], method = 'bilinear')
data[which(data$Month == 12),"salinityWOA13"] <- raster::extract(x = ras12, y = data[which(data$Month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data$salinityWOA13)
# Remove obs with SSS <= 20
# data2 <- data[-which(data$salinityWOA13 < 20),] # No idea why this particular line did not work while the following did....
data2 <- data[!(data$salinityWOA13 < 20),]
dim(data2) # 2805
# summary(data2$salinityWOA13)
# Remove NAs 
data3 <- data2 %>% drop_na(salinityWOA13)
dim(data3) # 2765
summary(data3$salinityWOA13)
### Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
save(data3, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------

##### 4°) Apply v5.1 and v5.2 to the two v4 datasets

### A) Remove any obs that has a sampling depth below 500m (v5.1)
# 1) v4v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.1/")
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
summary(data$Depth)
# Sometimes you still have negative depths for some reason... -> absolute
data$Depth <- abs(data$Depth)
# Remove obs with depth > 500
data2 <- data[!(data$Depth > 500),]
summary(data2$Depth)
dim(data2)
# Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
save(data2, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


# 2) v4v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
summary(data$Depth)
# Sometimes you still have negative depths for some reason... -> absolute
data$Depth <- abs(data$Depth)
# Remove obs with depth > 500
data2 <- data[!(data$Depth > 500),]
summary(data2$Depth)
dim(data2) # 2746
# Save
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
save(data2, file = "Pteropoda_MAREDAT_31_05_18.Rdata")



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
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
# Create salinity vector
data$MLD <- NA
# Fill it according to month
data[which(data$Month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$Month == 1),c("x","y")], method = 'bilinear')
data[which(data$Month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$Month == 2),c("x","y")], method = 'bilinear')
data[which(data$Month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$Month == 3),c("x","y")], method = 'bilinear')
data[which(data$Month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$Month == 4),c("x","y")], method = 'bilinear')
data[which(data$Month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$Month == 5),c("x","y")], method = 'bilinear')
data[which(data$Month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$Month == 6),c("x","y")], method = 'bilinear')
data[which(data$Month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$Month == 7),c("x","y")], method = 'bilinear')
data[which(data$Month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$Month == 8),c("x","y")], method = 'bilinear')
data[which(data$Month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$Month == 9),c("x","y")], method = 'bilinear')
data[which(data$Month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$Month == 10),c("x","y")], method = 'bilinear')
data[which(data$Month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$Month == 11),c("x","y")], method = 'bilinear')
data[which(data$Month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$Month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data)
# Remove obs based on MLD
data2 <- data[-which(data$Depth > (data$MLD + 50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
summary(data2$MLD)
data3 <- data2[!(data2$Depth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
#data3 <- data2 %>% drop_na(MLD) # Remove NAs
dim(data3) # 844
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
save(data3, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


# 2) v4v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v4v3.2/")
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
# Create salinity vector
data$MLD <- NA
# Fill it according to month
data[which(data$Month == 1),"MLD"] <- raster::extract(x = ras1, y = data[which(data$Month == 1),c("x","y")], method = 'bilinear')
data[which(data$Month == 2),"MLD"] <- raster::extract(x = ras2, y = data[which(data$Month == 2),c("x","y")], method = 'bilinear')
data[which(data$Month == 3),"MLD"] <- raster::extract(x = ras3, y = data[which(data$Month == 3),c("x","y")], method = 'bilinear')
data[which(data$Month == 4),"MLD"] <- raster::extract(x = ras4, y = data[which(data$Month == 4),c("x","y")], method = 'bilinear')
data[which(data$Month == 5),"MLD"] <- raster::extract(x = ras5, y = data[which(data$Month == 5),c("x","y")], method = 'bilinear')
data[which(data$Month == 6),"MLD"] <- raster::extract(x = ras6, y = data[which(data$Month == 6),c("x","y")], method = 'bilinear')
data[which(data$Month == 7),"MLD"] <- raster::extract(x = ras7, y = data[which(data$Month == 7),c("x","y")], method = 'bilinear')
data[which(data$Month == 8),"MLD"] <- raster::extract(x = ras8, y = data[which(data$Month == 8),c("x","y")], method = 'bilinear')
data[which(data$Month == 9),"MLD"] <- raster::extract(x = ras9, y = data[which(data$Month == 9),c("x","y")], method = 'bilinear')
data[which(data$Month == 10),"MLD"] <- raster::extract(x = ras10, y = data[which(data$Month == 10),c("x","y")], method = 'bilinear')
data[which(data$Month == 11),"MLD"] <- raster::extract(x = ras11, y = data[which(data$Month == 11),c("x","y")], method = 'bilinear')
data[which(data$Month == 12),"MLD"] <- raster::extract(x = ras12, y = data[which(data$Month == 12),c("x","y")], method = 'bilinear')
# Check
summary(data)
# Remove obs based on MLD
data2 <- data[-which(data$Depth > (data$MLD + 50)),] # Remove all obs whose net tows went deeper than the average monthly MLD + 50m
summary(data2$MLD)
data3 <- data2[!(data2$Depth > 150 & data2$MLD == NA),] # If some MLD have NAs, then just stick to a rather strick depth threshold of 150m depth
#data3 <- data2 %>% drop_na(MLD) # Remove NAs
dim(data3) # 1018
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
save(data3, file = "Pteropoda_MAREDAT_31_05_18.Rdata")


### ----------------------------------------------------------------------------------------------------------------------------

##### 5°) Prepare the data to be merged with OBIS & GBIF at v6

### v5.1v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.1/")
dir()
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
colnames(data)

# Re-name depth, day, month, year, date and taxon/ species
colnames(data)[3:8] <- c("depth","day","month","year","eventdate","species")
# Re-name Reference to rightsholder
colnames(data)[12] <- c("rightsholder")
### Remove abund & biomass: columns 9 & 10
colnames(data) # columns 9 & 10
data <- data[,c(1:8,11:16)]
# Remove id: col 12
colnames(data) # col 12
data <- data[,c(1:11,13:14)]
# Re-name Additional information -> basisofrecord
colnames(data)[9] <- c("basisofrecord")
# Re-name Source -> source
colnames(data)[11] <- c("source")
data$source <- "MAREDAT"

### Add genus name using strsplit on species name
data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)

data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames(data)
head(data)

# Add class and phylum
data$class <- "Gastropoda"
data$phylum <- "Mollusca"

### Add higher classification: genus, family, order, class (Gastropoda)
# unique(data$genus)
data[data$genus == "Clio","family"] <- "Cliidae"
data[data$genus == "Clio","order"] <- "Thecosomata"
data[data$genus == "Limacina","family"] <- "Limacinidae"
data[data$genus == "Limacina","order"] <- "Thecosomata"
data[data$genus == "Clione","family"] <- "Clionidae"
data[data$genus == "Clione","order"] <- "Gymnosomata"
data[data$genus == "Spongiobranchaea","family"] <- "Pneumodermatidae"
data[data$genus == "Spongiobranchaea","order"] <- "Gymnosomata"
data[data$genus == "Heliconoides","family"] <- "Limacinidae"
data[data$genus == "Heliconoides","order"] <- "Thecosomata"
data[data$genus == "Styliola","family"] <- "Creseidae"
data[data$genus == "Styliola","order"] <- "Thecosomata"
data[data$genus == "Diacria","family"] <- "Cavoliniidae"
data[data$genus == "Diacria","order"] <- "Thecosomata"
data[data$genus == "Cavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Cavolinia","order"] <- "Thecosomata"
data[data$genus == "Diacavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Diacavolinia","order"] <- "Thecosomata"
data[data$genus == "Hyalocylis","family"] <- "Creseidae"
data[data$genus == "Hyalocylis","order"] <- "Thecosomata"
data[data$genus == "Creseis","family"] <- "Creseidae"
data[data$genus == "Creseis","order"] <- "Thecosomata"

head(data)
colnames(data)

### Finish by saving in proper directory
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
save(data, file = "Thecosomata_MAREDAT_31_05_18.Rdata")





### v5.1v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.1v3.2/")
dir()
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
colnames(data) 
# Re-name depth, day, month, year, date and taxon/ species
colnames(data)[3:8] <- c("depth","day","month","year","eventdate","species")
# Re-name Reference to rightsholder
colnames(data)[12] <- c("rightsholder")
### Remove abund & biomass: columns 9 & 10
colnames(data) # columns 9 & 10
data <- data[,c(1:8,11:length(data))]
# Remove id: col 12
colnames(data) # col 12
data <- data[,c(1:11,13:length(data))]
# Re-name Additional information -> basisofrecord
colnames(data)[9] <- c("basisofrecord")
# Re-name Source -> source
colnames(data)[11] <- c("source")
data$source <- "MAREDAT"

### Add genus name using strsplit on species name
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)

data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames(data)
head(data)

# Add class and phylum
data$class <- "Gastropoda"
data$phylum <- "Mollusca"

### Add higher classification: genus, family, order, class (Gastropoda)
unique(data$genus) # Some NAS ?
data <- data %>% drop_na(species)
data[data$genus == "Clio","family"] <- "Cliidae"
data[data$genus == "Clio","order"] <- "Thecosomata"
data[data$genus == "Limacina","family"] <- "Limacinidae"
data[data$genus == "Limacina","order"] <- "Thecosomata"
data[data$genus == "Clione","family"] <- "Clionidae"
data[data$genus == "Clione","order"] <- "Gymnosomata"
data[data$genus == "Spongiobranchaea","family"] <- "Pneumodermatidae"
data[data$genus == "Spongiobranchaea","order"] <- "Gymnosomata"
data[data$genus == "Heliconoides","family"] <- "Limacinidae"
data[data$genus == "Heliconoides","order"] <- "Thecosomata"
data[data$genus == "Styliola","family"] <- "Creseidae"
data[data$genus == "Styliola","order"] <- "Thecosomata"
data[data$genus == "Diacria","family"] <- "Cavoliniidae"
data[data$genus == "Diacria","order"] <- "Thecosomata"
data[data$genus == "Cavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Cavolinia","order"] <- "Thecosomata"
data[data$genus == "Diacavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Diacavolinia","order"] <- "Thecosomata"
data[data$genus == "Hyalocylis","family"] <- "Creseidae"
data[data$genus == "Hyalocylis","order"] <- "Thecosomata"
data[data$genus == "Creseis","family"] <- "Creseidae"
data[data$genus == "Creseis","order"] <- "Thecosomata"

head(data)
colnames(data)

### Finish by saving in proper directory
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
save(data, file = "Thecosomata_MAREDAT_31_05_18.Rdata")






### v5.2v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.1/")
dir()
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
colnames(data)
# Re-name depth, day, month, year, date and taxon/ species
colnames(data)[3:8] <- c("depth","day","month","year","eventdate","species")
# Re-name Reference to rightsholder
colnames(data)[12] <- c("rightsholder")
### Remove abund & biomass: columns 9 & 10
colnames(data) # columns 9 & 10
data <- data[,c(1:8,11:length(data))]
# Remove id: col 12
colnames(data) # col 12
data <- data[,c(1:11,13:length(data))]
# Re-name Additional information -> basisofrecord
colnames(data)[9] <- c("basisofrecord")
# Re-name Source -> source
colnames(data)[11] <- c("source")
data$source <- "MAREDAT"

### Add genus name using strsplit on species name
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)

data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames(data)
head(data)

# Add class and phylum
data$class <- "Gastropoda"
data$phylum <- "Mollusca"

### Add higher classification: genus, family, order, class (Gastropoda)
unique(data$genus) # Some NAS ?
data <- data %>% drop_na(species)
data[data$genus == "Clio","family"] <- "Cliidae"
data[data$genus == "Clio","order"] <- "Thecosomata"
data[data$genus == "Limacina","family"] <- "Limacinidae"
data[data$genus == "Limacina","order"] <- "Thecosomata"
data[data$genus == "Clione","family"] <- "Clionidae"
data[data$genus == "Clione","order"] <- "Gymnosomata"
data[data$genus == "Spongiobranchaea","family"] <- "Pneumodermatidae"
data[data$genus == "Spongiobranchaea","order"] <- "Gymnosomata"
data[data$genus == "Heliconoides","family"] <- "Limacinidae"
data[data$genus == "Heliconoides","order"] <- "Thecosomata"
data[data$genus == "Styliola","family"] <- "Creseidae"
data[data$genus == "Styliola","order"] <- "Thecosomata"
data[data$genus == "Diacria","family"] <- "Cavoliniidae"
data[data$genus == "Diacria","order"] <- "Thecosomata"
data[data$genus == "Cavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Cavolinia","order"] <- "Thecosomata"
data[data$genus == "Diacavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Diacavolinia","order"] <- "Thecosomata"
data[data$genus == "Hyalocylis","family"] <- "Creseidae"
data[data$genus == "Hyalocylis","order"] <- "Thecosomata"
data[data$genus == "Creseis","family"] <- "Creseidae"
data[data$genus == "Creseis","order"] <- "Thecosomata"

head(data)
colnames(data)

### Finish by saving in proper directory
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
save(data, file = "Thecosomata_MAREDAT_31_05_18.Rdata")






### v5.2v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v5.2v3.2/")
dir()
data <- get(load("Pteropoda_MAREDAT_31_05_18.Rdata"))
colnames(data)
# Re-name depth, day, month, year, date and taxon/ species
colnames(data)[3:8] <- c("depth","day","month","year","eventdate","species")
# Re-name Reference to rightsholder
colnames(data)[12] <- c("rightsholder")
### Remove abund & biomass: columns 9 & 10
colnames(data) # columns 9 & 10
data <- data[,c(1:8,11:length(data))]
# Remove id: col 12
colnames(data) # col 12
data <- data[,c(1:11,13:length(data))]
# Re-name Additional information -> basisofrecord
colnames(data)[9] <- c("basisofrecord")
# Re-name Source -> source
colnames(data)[11] <- c("source")
data$source <- "MAREDAT"

### Add genus name using strsplit on species name
colnames2add <- c("family","order","class","phylum","identifiedby","recordedby","institutioncode")
data[setdiff(colnames2add, colnames(data))] <- NA
# Check
colnames(data)

data$genus <- do.call(rbind, str_split(as.character(data$species), "_"))[,1]
colnames(data)
head(data)

# Add class and phylum
data$class <- "Gastropoda"
data$phylum <- "Mollusca"

### Add higher classification: genus, family, order, class (Gastropoda)
unique(data$genus) # Some NAS ?
data <- data %>% drop_na(species)
data[data$genus == "Clio","family"] <- "Cliidae"
data[data$genus == "Clio","order"] <- "Thecosomata"
data[data$genus == "Limacina","family"] <- "Limacinidae"
data[data$genus == "Limacina","order"] <- "Thecosomata"
data[data$genus == "Clione","family"] <- "Clionidae"
data[data$genus == "Clione","order"] <- "Gymnosomata"
data[data$genus == "Spongiobranchaea","family"] <- "Pneumodermatidae"
data[data$genus == "Spongiobranchaea","order"] <- "Gymnosomata"
data[data$genus == "Heliconoides","family"] <- "Limacinidae"
data[data$genus == "Heliconoides","order"] <- "Thecosomata"
data[data$genus == "Styliola","family"] <- "Creseidae"
data[data$genus == "Styliola","order"] <- "Thecosomata"
data[data$genus == "Diacria","family"] <- "Cavoliniidae"
data[data$genus == "Diacria","order"] <- "Thecosomata"
data[data$genus == "Cavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Cavolinia","order"] <- "Thecosomata"
data[data$genus == "Diacavolinia","family"] <- "Cavoliniidae"
data[data$genus == "Diacavolinia","order"] <- "Thecosomata"
data[data$genus == "Hyalocylis","family"] <- "Creseidae"
data[data$genus == "Hyalocylis","order"] <- "Thecosomata"
data[data$genus == "Creseis","family"] <- "Creseidae"
data[data$genus == "Creseis","order"] <- "Thecosomata"

head(data)
colnames(data)

### Finish by saving in proper directory
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
save(data, file = "Thecosomata_MAREDAT_31_05_18.Rdata")

# Check number of observations per v6 dataset
datasets <- c("v6-v5.1v3.1","v6-v5.1v3.2","v6-v5.2v3.1","v6-v5.2v3.2")

for(d in datasets) {
		setwd(paste("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/",d,"/", sep = ""))
		data <- get(load("Thecosomata_MAREDAT_31_05_18.Rdata"))
		message(paste(nrow(data), sep = ""))
		message("nsp = ",paste(length(unique(data$species)), sep = ""))
		rm(data)
		gc()
} # eo for loop

















