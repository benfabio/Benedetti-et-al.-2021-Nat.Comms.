
### ==============================================================================================================================

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")

### ==============================================================================================================================

### 1°) Get all zooplankton data fitted with env layers from v9 datasets

# Go to proper dir
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
WD <- getwd()
files <- dir()[grep("30_04_19",dir())]

# Get colnames to retrive and rbind per dataset
d <- get(load(dir()[6]))
names <- colnames(d)[c(1:8,10:13,21:length(d))]
rm(d)

# Get all obs data
require("parallel")
# f <- "Alciopidae_matched_30_04_19.Rdata"
matched <- mclapply(files, function(f) {
				data <- get(load(f))
				message(paste("Reading ",f, sep = ""))
				return( subset(data, select = names) )	
			}, mc.cores = 15
) # eo lapply
match <- do.call(rbind, matched)
rm(matched, names) ; gc()
dim(match) # 1 184 690 for v3.2 ; 764'159 for 3.1


### 13/09/2018: Correct or remove some species' occurrences
# Remove the following because their distribution are way too biased
# Species with n occs >= 200
match <- match[which(match$species != "Drepanopus_pectinatus"),]
match <- match[which(match$species != "Tethys_vagina"),]
match <- match[which(match$species != "Agetus_limbatus"),]
match <- match[which(match$species != "Beroe_ovata"),]
match <- match[which(match$species != "Beroe_forskalii"),]
match <- match[which(match$species != "Bolinopsis_vitrea"),]
match <- match[which(match$species != "Doliolina_(Doliolina)_muelleri"),]
match <- match[which(match$species != "Bradydius_armatus"),]
match <- match[which(match$species != "Doliolum_nationalis"),]
# Species with n occs < 200
match <- match[which(match$species != "Acanthomysis_longicornis"),]
match <- match[which(match$species != "Acartia_(Acanthacartia)_bifilosa"),]
match <- match[which(match$species != "Acartia_(Acanthacartia)_fossae"),]
match <- match[which(match$species != "Acartia_(Acanthacartia)_tonsa"),]
match <- match[which(match$species != "Acartia_(Odontacartia)_amboinensis"),]
match <- match[which(match$species != "Alternochelata_sikorai"),]
match <- match[which(match$species != "Amallothrix_propinqua"),]
match <- match[which(match$species != "Bradyidius_armatus"),]
match <- match[which(match$species != "Dolioletta_gegenbauri"),]
match <- match[which(match$species != "Fritillaria_formica"),]
match <- match[which(match$species != "Globoquadrina_conglomerata"),]
match <- match[which(match$species != "Globorotalia_tumida"),]
match <- match[which(match$species != "Labidocera_aestiva"),]
match <- match[which(match$species != "Mesopodopsis_slabberi"),]
match <- match[which(match$species != "Oikopleura_(Coecaria)_fusiformis"),]
match <- match[which(match$species != "Oncaea_waldemari"),]
match <- match[which(match$species != "Parasagitta_euneritica"),]
match <- match[which(match$species != "Parasagitta_friderici"),]
match <- match[which(match$species != "Parvocalanus_elegans"),]
match <- match[which(match$species != "Schistomysis_kervillei"),]
match <- match[which(match$species != "Schistomysis_spiritus"),]
match <- match[which(match$species != "Scolecithricella_orientalis"),]
match <- match[which(match$species != "Tetrathyrus_forcipatus"),]

# Correct C.finn occurrences:
# - remove occ below 0° long
match2 <- match[!(match$species == "Calanus_finmarchicus" & match$y < 0),]
# dim(match2)
# - remove the 3 Mediterranean occurrences
match2 <- match2[!(match2$species == "Calanus_finmarchicus" & match2$x < 50 & match2$x > 0 & match2$y < 50 & match2$y > 25),]
# dim(match2)
# - remove the N. Pacific occurrences (in two times because positive and negative longitudes)
match2 <- match2[!(match2$species == "Calanus_finmarchicus" & match2$x < -100 & match2$x > -180 & match2$y < 70 & match2$y > 30),]
# dim(match2)
match2 <- match2[!(match2$species == "Calanus_finmarchicus" & match2$x < 180 & match2$x > 120 & match2$y < 70 & match2$y > 30),]
# dim(match2)
# Correct C. helgolandicus
# - remove occ below 25° long
match2 <- match2[!(match2$species == "Calanus_helgolandicus" & match2$y < 25),]
# dim(match2)
# - remove the Pacific occurrences (in two times because positive and negative longitudes)
match2 <- match2[!(match2$species == "Calanus_helgolandicus" & match2$x < -100),]
match2 <- match2[!(match2$species == "Calanus_helgolandicus" & match2$x > 60),]
# dim(match2)
# Correct C. armata, apply same criyteria as C. helgolandicus 
match2 <- match2[!(match2$species == "Candacia_armata" & match2$y < 25),]
# dim(match2)
# - remove the Pacific occurrences (in two times because positive and negative longitudes)
match2 <- match2[!(match2$species == "Candacia_armata" & match2$x < -100),]
match2 <- match2[!(match2$species == "Candacia_armata" & match2$x > 60),]
# dim(match2)
# Correct C. typicus : remove latitudes < 0 (n = 2)
match2 <- match2[!(match2$species == "Centropages_typicus" & match2$y < 0),]
# Correct O. frigida : remove latitudes > 0 (n = 1)
match2 <- match2[!(match2$species == "Oithona_frigida" & match2$y > 0),]
# Correct Metridia_curticauda : remove positive latitudes
match2 <- match2[!(match2$species == "Metridia_curticauda" & match2$y > 0),]
# Correct Vibilia_antarctica : remove positive latitudes
match2 <- match2[!(match2$species == "Vibilia_antarctica" & match2$y > 0),]
# Correct Aetideopsis armatus
match2 <- match2[!(match2$species == "Aetideopsis_armatus" & match2$x < -100),]
match2 <- match2[!(match2$species == "Aetideopsis_armatus" & match2$x > 100),]
# Correct Augaptilus glacialis armatus
match2 <- match2[!(match2$species == "Augaptilus_glacialis" & match2$y > 0 & match2$y < 50),]
# Correct Clausocalanus_brevipes armatus
match2 <- match2[!(match2$species == "Clausocalanus_brevipes" & match2$y > 0),]
# Correct Ditrichocorycaeus_anglicus:
match2 <- match2[!(match2$species == "Ditrichocorycaeus_anglicus" & match2$y < 0),]
match2 <- match2[!(match2$species == "Ditrichocorycaeus_anglicus" & match2$x < -100),]
# Correct Euphausia_frigida:
match2 <- match2[!(match2$species == "Euphausia_frigida" & match2$y > 0),]
# Correct Pseudocalanus_elongatus:
match2 <- match2[!(match2$species == "Pseudocalanus_elongatus" & match2$y < 20),]
match2 <- match2[!(match2$species == "Pseudocalanus_elongatus" & match2$x < -100),]
# Correct Rhincalanus_gigas:
match2 <- match2[!(match2$species == "Rhincalanus_gigas" & match2$y > -20),]
# Correct Temorites_brevis:
match2 <- match2[!(match2$species == "Temorites_brevis" & match2$y > 0 & match2$y < 50),]
# Correct Tomopteris_(Johnstonella)_helgolandica:
match2 <- match2[!(match2$species == "Tomopteris_(Johnstonella)_helgolandica" & match2$y < 0),]
# Correct Triconia borealis:
match2 <- match2[!(match2$species == "Triconia_borealis" & match2$y < 37),]

### 03/12/18
# Correct Aetideopsis_rostrata:
match2 <- match2[!(match2$species == "Aetideopsis_rostrata" & match2$y < 0),]
# Correct Candacia_columbiae:
match2 <- match2[!(match2$species == "Candacia_columbiae" & match2$y < 0),]
# Correct Candacia_discaudata:
match2 <- match2[!(match2$species == "Candacia_discaudata" & match2$x < 0),]
# Correct Epilabidocera_amphitrites:
match2 <- match2[!(match2$species == "Epilabidocera_amphitrites" & match2$x > 0),]
# Correct Hyperia_galba:
match2 <- match2[!(match2$species == "Hyperia_galba" & match2$y < 0),]
# Correct Neocalanus_tonsus:
match2 <- match2[!(match2$species == "Neocalanus_tonsus" & match2$y > 0),]
# Correct Oikopleura_(Vexillaria)_labradoriensis:
match2 <- match2[!(match2$species == "Oikopleura_(Vexillaria)_labradoriensis" & match2$y < 0),]
# Correct Oithona_attenuata:
match2 <- match2[!(match2$species == "Oithona_attenuata" & match2$x < 0),]
# Correct Paraheterorhabdus_compactus:
match2 <- match2[!(match2$species == "Paraheterorhabdus_compactus" & match2$y < 0),]
# Correct Pleurobrachia_pileus:
match2 <- match2[!(match2$species == "Pleurobrachia_pileus" & match2$x < 50 & match2$x > 30 & match2$y < 50 & match2$y > 25),]
# Correct Pleuromamma_scutullata:
match2 <- match2[!(match2$species == "Pleuromamma_scutullata" & match2$x < 0 & match2$x > -75),]
# Correct Pneumodermopsis_paucidens:
match2 <- match2[!(match2$species == "Pneumodermopsis_paucidens" & match2$y < 0),]
# Correct Travisiopsis_levinseni:
match2 <- match2[!(match2$species == "Travisiopsis_levinseni" & match2$y > 0),]
# Correct Tryphana_malmii:
match2 <- match2[!(match2$species == "Tryphana_malmii" & match2$y < 0),]
# Correct Vettoria_granulosa:
match2 <- match2[!(match2$species == "Vettoria_granulosa" & match2$y < 0),]

# Get vector of species names
all_spp_names <- unique(match2$species)

### Compute n obs per species
counts <- data.frame(match2 %>% 
		group_by(species) %>%
		summarise(phylum = unique(phylum), n = n()) 
) # eo ddf

### Define the pool of species with enough observations to draw the psAbs, let's start with n >= 300
# length(counts[counts$n >= 100,"species"])
# length(counts[counts$n >= 75,"species"])
# length(counts[counts$n >= 50,"species"])

species <- counts[counts$n >= 85,"species"]
# species <- counts[counts$n < 300 & counts$n >= 200,"species"]
species

### Define the background sampling strategy (either "total" or "target_group") and draw psAbs 
#strategy <- "group" # corresponds to d = 4 in Damiano's code
strategy <- "total"  # corresponds to d = 10 in Damiano's code

# Specify parameters to stratify the sampled environment: ** ADJUST IN CASE**
vec.strat <- c("SST","MLD1") # same as Damiano

### Define the function that you will use in a mclapply to derive the psAbs of each spp 
# sp <- "Metridia_lucens"

### 19/09/18: Vector os species to re-run
species2redo <- counts[counts$n >= 85 & counts$phylum == "Foraminifera","species"]
				
for(sp in species2redo ) {
	
			message(paste("Drawing psAbs for ",sp, " ============================================ ", sep = ""))
			group <- unique(match2[match2$species == sp,"phylum"])
								
			# Get species data
			all_id <- match2[match2$species == sp ,]
			n <- nrow(all_id)

			# Remove NAs from background data regarding these two variables chosen
			if( length(which(is.na(all_id[,which(names(all_id) == vec.strat[1])] ))) != 0 ) {
					all_id <- all_id[-which(is.na(all_id[,which(names(all_id) == vec.strat[1])])),]
			} # remove NAs regarding Var 1

			if( length(which(is.na(all_id[,which(names(all_id) == vec.strat[2])]))) != 0 ) { 
					all_id <- all_id[-which(is.na(all_id[,which(names(all_id) == vec.strat[2])])),]
			} # remove NAs regarding Var 2

			if ( length(which(is.na(all_id[,which(names(all_id) == vec.strat[1])]))) == 0 & length(which(is.na(all_id[,which(names(all_id)==vec.strat[2])]))) == 0 ) {
					all_id <- all_id	
			} #

			if( nrow(all_id) >= 85 ) {
			
				# Choice of background data depending on the "strategy"
				if( strategy == "total" ) {					
						bckgrnd <- match2[which(match2$species != sp & match2$SSS > 20 & match2$Bathy < -175),] 					
				} else {
					
					# Then use the groups' data as bckgrnd
					if(group %in% c("Arthropoda","Chaetognatha","Cnidaria","Mollusca"))
							bckgrnd <- match2[match2$phylum == group,]	
					else if (group == "Ctenophora") {
							bckgrnd <- match2[match2$phylum == "Cnidaria",]	
					} else {
							bckgrnd <- match2
					}
							
				} # eo if else loop
								
				### Specify range of variables and strata into which the values fall to drive sampling of absences proportionally to the overall presences points in the strata
				x_envir <- bckgrnd[,c(vec.strat[1])]
				y_envir <- bckgrnd[,c(vec.strat[2])]
				# Split ranges into environmental strata
				breaks <- 9

				### Create a matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
				x_breaks <- classIntervals(na.omit(x_envir), breaks, style = "equal")
				x_matrix <- cbind(x_breaks$brks[1:breaks], x_breaks$brks[2:(breaks + 1)], ID = 1:breaks )
				colnames(x_matrix) <- c("low","up","ID")
				y_breaks <- classIntervals(na.omit(y_envir), breaks, style = "equal")
				y_matrix <- cbind(y_breaks$brks[1:breaks], y_breaks$brks[2:(breaks + 1)], ID = 1:breaks )
				colnames(y_matrix) <- c("low","up","ID")
				# Define vector of length of total points of environmental variable
				x_reclass <- c(1:length(x_envir))
				y_reclass <- c(1:length(y_envir))
								
				# Allocate points from full data to one of the nine environmental strata per variable
				for(i in 1:breaks) {	
						x_reclass[which(x_envir >= x_matrix[i,"low"] & x_envir <= x_matrix[i,"up"] )] <- x_matrix[i,"ID"]
						y_reclass[which(y_envir >= y_matrix[i,"low"] & y_envir <= y_matrix[i,"up"] )] <- y_matrix[i,"ID"] 		
				} # eo for loop

				### Create an ID indicating the stratum (unique combination of variables) into which each point falls in full data-frame
				bckgrnd$x_rcls <- x_reclass
				bckgrnd$y_rcls <- y_reclass
				bckgrnd$xy_rcls <- x_reclass+10*y_reclass
								
				print( paste0(sp,", ", group, " | n = ", n, " | drawing psAbs")) # eo print
								
				### Extract frequencies by which points/sites of the target group fall into environmental strata. 
				# Then, derive the number of desired absences for the focal model species per stratum. 
				xy_rcls_freq <- data.frame( table(bckgrnd$xy_rcls) / length(bckgrnd$x_rcls) ) 		
				# Give name to column
				colnames(xy_rcls_freq)[1] <- "xy_rcls" 
				# Convert to numeric
				xy_rcls_freq$xy_rcls <- as.numeric(as.character(xy_rcls_freq$xy_rcls)) 

				### Add desired background points to be produced per stratum: generally 10 x more absences than presences
				xy_rcls_freq$prop_abs <- (nrow(all_id)*10)*xy_rcls_freq$Freq

				# To round desired absences to integer: adds column difference between smaller closest integer and desired number
				xy_rcls_freq$prop_abs_0 <- ( xy_rcls_freq$prop_abs - floor(xy_rcls_freq$prop_abs) )

				# To add column with random number between 0 and 1 (with steps of 0.01)
				xy_rcls_freq$prob <- sample(seq(0, 1, 0.01), nrow(xy_rcls_freq), replace = TRUE)

				# To add column with "1"
				xy_rcls_freq$absences <- 1

				# Round up absences for random subset
				xy_rcls_freq$absences[which(xy_rcls_freq$prop_abs_0 > xy_rcls_freq$prob)] <- ceiling(xy_rcls_freq$prop_abs[which(xy_rcls_freq$prop_abs_0 > xy_rcls_freq$prob)])
			
				# Round absences down for random subset
				xy_rcls_freq$absences[which(xy_rcls_freq$prop_abs_0 < xy_rcls_freq$prob)] <- floor(xy_rcls_freq$prop_abs[which(xy_rcls_freq$prop_abs_0 < xy_rcls_freq$prob)]) 

				# Skip strata without presences
				absence_groups <- xy_rcls_freq[xy_rcls_freq$absences > 0,] 

				### Select backround data, here including the points/sites of the focal species ('overlapping background')
				absence_table <- bckgrnd
				gc()
			
				### Randomly Select background pts for target species in each stratum proportionally to the density of samples/points in the background
				nnn <- nrow(absence_groups) 
				require("parallel")
				### Need to parallel it
				psAbs <- mclapply(X = c(1:nnn), mc.cores = 30, FUN = function(i) {
	
						# 1. Select available absences within stratum in question
						message(paste(i, sep = ""))
						grp_abs_table <- absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]

						# 2. Define the max nb of absences that can be drawn 
						absence_num <- ifelse(
								# Test if the number of desired background pts is bigger than the available background points
								absence_groups[i,"absences"] > nrow( absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]),		
								# if TRUE the potential points are insufficient - however, save the number of available points as absence_num
								nrow(absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]),		
								# ELSE: save the number of desired background points as absence_num   		
								absence_groups[i,"absences"]
						) # eo if else loop
						
						# 3. Randomly sample the background points from the table containing all possible absences for the stratum in question
						sampled_grp_abs_table <- grp_abs_table[sample(1:nrow(grp_abs_table), size = absence_num),]
				
						return(sampled_grp_abs_table)
				
					} # eo fun
			
				) # eo lapply
			
				### Merge presences (obs = 1) with absences (obs = 0)
				pseudoabs <- data.frame(do.call("rbind", psAbs), obs = 0)
				occ_table <- rbind(data.frame(all_id, obs = 1),  pseudoabs[,c(colnames(all_id),"obs")])
				rm(psAbs, nnn)
				gc()
				
				# Remove some psAbs that were drawn in SSS < 20
				occ_table <- occ_table[occ_table$SSS > 20,]

				### Create column with weights = 1; weights are associated with presences and absences for modelling
				occ_table$weights <- 1
				# Compute ratio of presences to absences
				abs_ratio <- nrow(occ_table[occ_table$obs == 1,]) / nrow(occ_table[occ_table$obs == 0,])
			
				# Add the ratio as weight for the psAbs
				occ_table$weights[occ_table$obs == 0] <- abs_ratio # For observation that are absences we add the ratio
				row.names(occ_table) <- c(1:nrow(occ_table)) # # Add row ID (order is important for later cross-validation procedure/TSS calculation)
			
				### Plot spatial distrib and save it (for info)
				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.1/total_background/maps/")
			
				require("ggplot2")
				map <- ggplot() +
 					geom_point(aes(x = x, y = y), data = occ_table[which(occ_table$obs == 0),], fill = "#d73027", pch = 21, colour = "black", alpha = 0.5) +
 					geom_point(aes(x = x, y = y), data = occ_table[which(occ_table$obs == 1),], fill = "#4575b4", pch = 21, colour = "black") +
 					scale_x_continuous(limits = c(-180,180)) + scale_y_continuous(limits = c(-90,90)) +
 					theme_bw() + xlab("Longitude") + ylab("Latitude") + coord_quickmap()
		
				ggsave(plot = map, filename = paste("map_psAbs_",sp,"_",strategy,"_","v3.1.jpg", sep = ""), dpi = 300, width = 12, height = 9)

				### Save the data to train some ENMs later
				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.1/total_background/")
				message(paste("Saving species dataset for ",sp, " ============================================ ", sep = ""))
				write.table(occ_table, paste("data_",strategy,"_",sp,".txt", sep = ""), sep = ";")
			
				### Clean some stuff 
				rm(occ_table, map, abs_ratio, pseudoabs, absence_table, xy_rcls_freq, y_reclass, x_reclass, y_matrix, x_matrix, y_envir, x_envir)
				gc()
				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
				
			} else {
				
				### Clean some stuff 
				rm(group, all_id)
				gc()
				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
				
			} # eo else if loop based on nrow(all_id)
			
} # eo for loop



### ==============================================================================================================================
### ==============================================================================================================================

# library("rgeos")
# library("raster")
# library("maptools")
# library("rgdal")
# library("tidyverse")
# library("stringr")
# library("reshape2")
# library("geosphere")
# library("ncdf4")
# library("classInt")
#
#
# ### 17/04/2019: Draw psAbs from rarefied dataset of occurrences
# ###	lower the threshold for drawing psAbs
#
# setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/")
# WD <- getwd()
# # dir()
#
# match <- read.table("rarefied_data_zooplankton_17_04_19.txt", sep = ";", h = T)
# dim(match)
#
# match <- match[which(match$species != "Drepanopus_pectinatus"),]
# match <- match[which(match$species != "Tethys_vagina"),]
# match <- match[which(match$species != "Agetus_limbatus"),]
# match <- match[which(match$species != "Beroe_ovata"),]
# match <- match[which(match$species != "Beroe_forskalii"),]
# match <- match[which(match$species != "Bolinopsis_vitrea"),]
# match <- match[which(match$species != "Doliolina_(Doliolina)_muelleri"),]
# match <- match[which(match$species != "Bradydius_armatus"),]
# match <- match[which(match$species != "Doliolum_nationalis"),]
# # Species with n occs < 200
# match <- match[which(match$species != "Acanthomysis_longicornis"),]
# match <- match[which(match$species != "Acartia_(Acanthacartia)_bifilosa"),]
# match <- match[which(match$species != "Acartia_(Acanthacartia)_fossae"),]
# match <- match[which(match$species != "Acartia_(Acanthacartia)_tonsa"),]
# match <- match[which(match$species != "Acartia_(Odontacartia)_amboinensis"),]
# match <- match[which(match$species != "Alternochelata_sikorai"),]
# match <- match[which(match$species != "Amallothrix_propinqua"),]
# match <- match[which(match$species != "Bradyidius_armatus"),]
# match <- match[which(match$species != "Dolioletta_gegenbauri"),]
# match <- match[which(match$species != "Fritillaria_formica"),]
# match <- match[which(match$species != "Globoquadrina_conglomerata"),]
# match <- match[which(match$species != "Globorotalia_tumida"),]
# match <- match[which(match$species != "Labidocera_aestiva"),]
# match <- match[which(match$species != "Mesopodopsis_slabberi"),]
# match <- match[which(match$species != "Oikopleura_(Coecaria)_fusiformis"),]
# match <- match[which(match$species != "Oncaea_waldemari"),]
# match <- match[which(match$species != "Parasagitta_euneritica"),]
# match <- match[which(match$species != "Parasagitta_friderici"),]
# match <- match[which(match$species != "Parvocalanus_elegans"),]
# match <- match[which(match$species != "Schistomysis_kervillei"),]
# match <- match[which(match$species != "Schistomysis_spiritus"),]
# match <- match[which(match$species != "Scolecithricella_orientalis"),]
# match <- match[which(match$species != "Tetrathyrus_forcipatus"),]
#
# # Correct C.finn occurrences:
# # - remove occ below 0° long
# match2 <- match[!(match$species == "Calanus_finmarchicus" & match$y < 0),]
# # dim(match2)
# # - remove the 3 Mediterranean occurrences
# match2 <- match2[!(match2$species == "Calanus_finmarchicus" & match2$x < 50 & match2$x > 0 & match2$y < 50 & match2$y > 25),]
# # dim(match2)
# # - remove the N. Pacific occurrences (in two times because positive and negative longitudes)
# match2 <- match2[!(match2$species == "Calanus_finmarchicus" & match2$x < -100 & match2$x > -180 & match2$y < 70 & match2$y > 30),]
# # dim(match2)
# match2 <- match2[!(match2$species == "Calanus_finmarchicus" & match2$x < 180 & match2$x > 120 & match2$y < 70 & match2$y > 30),]
# # dim(match2)
# # Correct C. helgolandicus
# # - remove occ below 25° long
# match2 <- match2[!(match2$species == "Calanus_helgolandicus" & match2$y < 25),]
# # dim(match2)
# # - remove the Pacific occurrences (in two times because positive and negative longitudes)
# match2 <- match2[!(match2$species == "Calanus_helgolandicus" & match2$x < -100),]
# match2 <- match2[!(match2$species == "Calanus_helgolandicus" & match2$x > 60),]
# # dim(match2)
# # Correct C. armata, apply same criyteria as C. helgolandicus
# match2 <- match2[!(match2$species == "Candacia_armata" & match2$y < 25),]
# # dim(match2)
# # - remove the Pacific occurrences (in two times because positive and negative longitudes)
# match2 <- match2[!(match2$species == "Candacia_armata" & match2$x < -100),]
# match2 <- match2[!(match2$species == "Candacia_armata" & match2$x > 60),]
# # dim(match2)
# # Correct C. typicus : remove latitudes < 0 (n = 2)
# match2 <- match2[!(match2$species == "Centropages_typicus" & match2$y < 0),]
# # Correct O. frigida : remove latitudes > 0 (n = 1)
# match2 <- match2[!(match2$species == "Oithona_frigida" & match2$y > 0),]
# # Correct Metridia_curticauda : remove positive latitudes
# match2 <- match2[!(match2$species == "Metridia_curticauda" & match2$y > 0),]
# # Correct Vibilia_antarctica : remove positive latitudes
# match2 <- match2[!(match2$species == "Vibilia_antarctica" & match2$y > 0),]
# # Correct Aetideopsis armatus
# match2 <- match2[!(match2$species == "Aetideopsis_armatus" & match2$x < -100),]
# match2 <- match2[!(match2$species == "Aetideopsis_armatus" & match2$x > 100),]
# # Correct Augaptilus glacialis armatus
# match2 <- match2[!(match2$species == "Augaptilus_glacialis" & match2$y > 0 & match2$y < 50),]
# # Correct Clausocalanus_brevipes armatus
# match2 <- match2[!(match2$species == "Clausocalanus_brevipes" & match2$y > 0),]
# # Correct Ditrichocorycaeus_anglicus:
# match2 <- match2[!(match2$species == "Ditrichocorycaeus_anglicus" & match2$y < 0),]
# match2 <- match2[!(match2$species == "Ditrichocorycaeus_anglicus" & match2$x < -100),]
# # Correct Euphausia_frigida:
# match2 <- match2[!(match2$species == "Euphausia_frigida" & match2$y > 0),]
# # Correct Pseudocalanus_elongatus:
# match2 <- match2[!(match2$species == "Pseudocalanus_elongatus" & match2$y < 20),]
# match2 <- match2[!(match2$species == "Pseudocalanus_elongatus" & match2$x < -100),]
# # Correct Rhincalanus_gigas:
# match2 <- match2[!(match2$species == "Rhincalanus_gigas" & match2$y > -20),]
# # Correct Temorites_brevis:
# match2 <- match2[!(match2$species == "Temorites_brevis" & match2$y > 0 & match2$y < 50),]
# # Correct Tomopteris_(Johnstonella)_helgolandica:
# match2 <- match2[!(match2$species == "Tomopteris_(Johnstonella)_helgolandica" & match2$y < 0),]
# # Correct Triconia borealis:
# match2 <- match2[!(match2$species == "Triconia_borealis" & match2$y < 37),]
#
# ### 03/12/18
# # Correct Aetideopsis_rostrata:
# match2 <- match2[!(match2$species == "Aetideopsis_rostrata" & match2$y < 0),]
# # Correct Candacia_columbiae:
# match2 <- match2[!(match2$species == "Candacia_columbiae" & match2$y < 0),]
# # Correct Candacia_discaudata:
# match2 <- match2[!(match2$species == "Candacia_discaudata" & match2$x < 0),]
# # Correct Epilabidocera_amphitrites:
# match2 <- match2[!(match2$species == "Epilabidocera_amphitrites" & match2$x > 0),]
# # Correct Hyperia_galba:
# match2 <- match2[!(match2$species == "Hyperia_galba" & match2$y < 0),]
# # Correct Neocalanus_tonsus:
# match2 <- match2[!(match2$species == "Neocalanus_tonsus" & match2$y > 0),]
# # Correct Oikopleura_(Vexillaria)_labradoriensis:
# match2 <- match2[!(match2$species == "Oikopleura_(Vexillaria)_labradoriensis" & match2$y < 0),]
# # Correct Oithona_attenuata:
# match2 <- match2[!(match2$species == "Oithona_attenuata" & match2$x < 0),]
# # Correct Paraheterorhabdus_compactus:
# match2 <- match2[!(match2$species == "Paraheterorhabdus_compactus" & match2$y < 0),]
# # Correct Pleurobrachia_pileus:
# match2 <- match2[!(match2$species == "Pleurobrachia_pileus" & match2$x < 50 & match2$x > 30 & match2$y < 50 & match2$y > 25),]
# # Correct Pleuromamma_scutullata:
# match2 <- match2[!(match2$species == "Pleuromamma_scutullata" & match2$x < 0 & match2$x > -75),]
# # Correct Pneumodermopsis_paucidens:
# match2 <- match2[!(match2$species == "Pneumodermopsis_paucidens" & match2$y < 0),]
# # Correct Travisiopsis_levinseni:
# match2 <- match2[!(match2$species == "Travisiopsis_levinseni" & match2$y > 0),]
# # Correct Tryphana_malmii:
# match2 <- match2[!(match2$species == "Tryphana_malmii" & match2$y < 0),]
# # Correct Vettoria_granulosa:
# match2 <- match2[!(match2$species == "Vettoria_granulosa" & match2$y < 0),]
#
# # Get vector of species names
# all_spp_names <- unique(match2$species)
#
# ### Compute n obs per species
# counts <- data.frame(match2 %>%
# 		group_by(species) %>%
# 		summarise(n = n())
# ) # eo ddf
#
# ### Define the pool of species with enough observations to draw the psAbs, let's start with n >= 300
# length(counts[counts$n >= 100,"species"])
# length(counts[counts$n >= 75,"species"])
# length(counts[counts$n >= 50,"species"])
#
# species <- counts[counts$n >= 50,"species"]
# # species <- counts[counts$n < 300 & counts$n >= 200,"species"]
# species
#
# strategy <- "total"  # corresponds to d = 10 in Damiano's code
# vec.strat <- c("SST","MLD1") # same as Damiano
#
# # sp <- "Acartia_(Acartiura)_clausi"
#
# # Remove "obs" column
# match2 <- match2[,c(1:34,36:37)]
#
# ### Generate psAbs
# for(sp in species[c(21:length(species))] ) {
#
# 			message(paste("Drawing psAbs for ",sp, " ============================================ ", sep = ""))
# 			group <- unique(match2[match2$species == sp,"phylum"])
#
# 			# Get species data
# 			all_id <- match2[match2$species == sp,]
# 			n <- nrow(all_id)
#
# 			# Remove NAs from background data regarding these two variables chosen
# 			if( length(which(is.na(all_id[,which(names(all_id) == vec.strat[1])] ))) != 0 ) {
# 					all_id <- all_id[-which(is.na(all_id[,which(names(all_id) == vec.strat[1])])),]
# 			} # remove NAs regarding Var 1
#
# 			if( length(which(is.na(all_id[,which(names(all_id) == vec.strat[2])]))) != 0 ) {
# 					all_id <- all_id[-which(is.na(all_id[,which(names(all_id) == vec.strat[2])])),]
# 			} # remove NAs regarding Var 2
#
# 			if ( length(which(is.na(all_id[,which(names(all_id) == vec.strat[1])]))) == 0 & length(which(is.na(all_id[,which(names(all_id)==vec.strat[2])]))) == 0 ) {
# 					all_id <- all_id
# 			} #
#
# 			if( nrow(all_id) >= 50 ) {
#
# 				# Choice of background data depending on the "strategy"
# 				if( strategy == "total" ) {
# 						bckgrnd <- match2[match2$species != sp,] # then get all zoo data
# 				} else {
#
# 					# Then use the groups' data as bckgrnd
# 					if(group %in% c("Arthropoda","Chaetognatha","Cnidara","Mollusca"))
# 							bckgrnd <- match2[match2$phylum == group,]
# 					else if (group == "Ctenophora") {
# 							bckgrnd <- match2[match2$phylum == "Cnidaria",]
# 					} else {
# 							bckgrnd <- match2
# 					}
#
# 				} # eo if else loop
#
# 				### Specify range of variables and strata into which the values fall to drive sampling of absences proportionally to the overall presences points in the strata
# 				x_envir <- bckgrnd[,c(vec.strat[1])]
# 				y_envir <- bckgrnd[,c(vec.strat[2])]
# 				# Split ranges into environmental strata
# 				breaks <- 9
#
# 				### Create a matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
# 				x_breaks <- classIntervals(na.omit(x_envir), breaks, style = "equal")
# 				x_matrix <- cbind(x_breaks$brks[1:breaks], x_breaks$brks[2:(breaks + 1)], ID = 1:breaks )
# 				colnames(x_matrix) <- c("low","up","ID")
# 				y_breaks <- classIntervals(na.omit(y_envir), breaks, style = "equal")
# 				y_matrix <- cbind(y_breaks$brks[1:breaks], y_breaks$brks[2:(breaks + 1)], ID = 1:breaks )
# 				colnames(y_matrix) <- c("low","up","ID")
# 				# Define vector of length of total points of environmental variable
# 				x_reclass <- c(1:length(x_envir))
# 				y_reclass <- c(1:length(y_envir))
#
# 				# Allocate points from full data to one of the nine environmental strata per variable
# 				for(i in 1:breaks) {
# 						x_reclass[which(x_envir >= x_matrix[i,"low"] & x_envir <= x_matrix[i,"up"] )] <- x_matrix[i,"ID"]
# 						y_reclass[which(y_envir >= y_matrix[i,"low"] & y_envir <= y_matrix[i,"up"] )] <- y_matrix[i,"ID"]
# 				} # eo for loop
#
# 				### Create an ID indicating the stratum (unique combination of variables) into which each point falls in full data-frame
# 				bckgrnd$x_rcls <- x_reclass
# 				bckgrnd$y_rcls <- y_reclass
# 				bckgrnd$xy_rcls <- x_reclass+10*y_reclass
#
# 				print( paste0(sp,", ", group, " | n = ", n, " | drawing psAbs")) # eo print
#
# 				### Extract frequencies by which points/sites of the target group fall into environmental strata.
# 				# Then, derive the number of desired absences for the focal model species per stratum.
# 				xy_rcls_freq <- data.frame( table(bckgrnd$xy_rcls) / length(bckgrnd$x_rcls) )
# 				# Give name to column
# 				colnames(xy_rcls_freq)[1] <- "xy_rcls"
# 				# Convert to numeric
# 				xy_rcls_freq$xy_rcls <- as.numeric(as.character(xy_rcls_freq$xy_rcls))
#
# 				### Add desired background points to be produced per stratum: generally 10 x more absences than presences
# 				xy_rcls_freq$prop_abs <- (nrow(all_id)*10)*xy_rcls_freq$Freq
#
# 				# To round desired absences to integer: adds column difference between smaller closest integer and desired number
# 				xy_rcls_freq$prop_abs_0 <- ( xy_rcls_freq$prop_abs - floor(xy_rcls_freq$prop_abs) )
#
# 				# To add column with random number between 0 and 1 (with steps of 0.01)
# 				xy_rcls_freq$prob <- sample(seq(0, 1, 0.01), nrow(xy_rcls_freq), replace = TRUE)
#
# 				# To add column with "1"
# 				xy_rcls_freq$absences <- 1
#
# 				# Round up absences for random subset
# 				xy_rcls_freq$absences[which(xy_rcls_freq$prop_abs_0 > xy_rcls_freq$prob)] <- ceiling(xy_rcls_freq$prop_abs[which(xy_rcls_freq$prop_abs_0 > xy_rcls_freq$prob)])
#
# 				# Round absences down for random subset
# 				xy_rcls_freq$absences[which(xy_rcls_freq$prop_abs_0 < xy_rcls_freq$prob)] <- floor(xy_rcls_freq$prop_abs[which(xy_rcls_freq$prop_abs_0 < xy_rcls_freq$prob)])
#
# 				# Skip strata without presences
# 				absence_groups <- xy_rcls_freq[xy_rcls_freq$absences > 0,]
#
# 				### Select backround data, here including the points/sites of the focal species ('overlapping background')
# 				absence_table <- bckgrnd
# 				gc()
#
# 				### Randomly Select background pts for target species in each stratum proportionally to the density of samples/points in the background
# 				nnn <- nrow(absence_groups)
# 				require("parallel")
# 				### Need to parallel it
# 				psAbs <- mclapply(X = c(1:nnn), mc.cores = 25, FUN = function(i) {
#
# 						# 1. Select available absences within stratum in question
# 						message(paste(i, sep = ""))
# 						grp_abs_table <- absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]
#
# 						# 2. Define the max nb of absences that can be drawn
# 						absence_num <- ifelse(
# 								# Test if the number of desired background pts is bigger than the available background points
# 								absence_groups[i,"absences"] > nrow( absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]),
# 								# if TRUE the potential points are insufficient - however, save the number of available points as absence_num
# 								nrow(absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]),
# 								# ELSE: save the number of desired background points as absence_num
# 								absence_groups[i,"absences"]
# 						) # eo if else loop
#
# 						# 3. Randomly sample the background points from the table containing all possible absences for the stratum in question
# 						sampled_grp_abs_table <- grp_abs_table[sample(1:nrow(grp_abs_table), size = absence_num),]
#
# 						return(sampled_grp_abs_table)
#
# 					} # eo fun
#
# 				) # eo lapply
#
# 				### Merge presences (obs = 1) with absences (obs = 0)
# 				pseudoabs <- data.frame(do.call("rbind", psAbs), obs = 0)
# 				occ_table <- rbind(data.frame(all_id, obs = 1),  pseudoabs[,c(colnames(all_id),"obs")])
# 				rm(psAbs, nnn)
# 				gc()
#
# 				# Remove some psAbs that were drawn in SSS < 20
# 				occ_table <- occ_table[occ_table$SSS > 20,]
#
# 				### Create column with weights = 1; weights are associated with presences and absences for modelling
# 				occ_table$weights <- 1
# 				# Compute ratio of presences to absences
# 				abs_ratio <- nrow(occ_table[occ_table$obs == 1,]) / nrow(occ_table[occ_table$obs == 0,])
#
# 				# Add the ratio as weight for the psAbs
# 				occ_table$weights[occ_table$obs == 0] <- abs_ratio # For observation that are absences we add the ratio
# 				row.names(occ_table) <- c(1:nrow(occ_table)) # # Add row ID (order is important for later cross-validation procedure/TSS calculation)
#
# 				### Plot spatial distrib and save it (for info)
# 				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1_rarefied/maps")
#
# 				require("ggplot2")
# 				require("scales")
# 				require("RColorBrewer")
# 				map <- ggplot() +
#  					geom_point(aes(x = x, y = y), data = occ_table[which(occ_table$obs == 0),], fill = "#d73027", pch = 21, colour = "black", alpha = 0.5) +
#  					geom_point(aes(x = x, y = y), data = occ_table[which(occ_table$obs == 1),], fill = "#4575b4", pch = 21, colour = "black") +
#  					scale_x_continuous(limits = c(-180,180)) + scale_y_continuous(limits = c(-90,90)) +
#  					theme_bw() + xlab("Longitude") + ylab("Latitude") + coord_quickmap()
#
# 				ggsave(plot = map, filename = paste("map_psAbs_",sp,"_",strategy,"_","v3.1_rarefied.jpg", sep = ""), dpi = 300, width = 13, height = 10)
#
# 				### Save the data to train some ENMs later
# 				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1_rarefied/")
# 				message(paste("Saving species dataset for ",sp, " ============================================ ", sep = ""))
# 				write.table(occ_table, paste("data_",strategy,"_",sp,".txt", sep = ""), sep = ";")
#
# 				### Clean some stuff
# 				rm(occ_table, map, abs_ratio, pseudoabs, absence_table, xy_rcls_freq, y_reclass, x_reclass, y_matrix, x_matrix, y_envir, x_envir)
# 				gc()
# 				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/")
#
# 			} else {
#
# 				### Clean some stuff
# 				rm(group, all_id)
# 				gc()
# 				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/")
#
# 			} # eo else if loop based on nrow(all_id)
#
# } # eo for loop


