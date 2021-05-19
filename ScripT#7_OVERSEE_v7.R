
##### 04/06/2018: R Script to homogenize the taxonomic nomenclature (species names) across all v6 zooplankton datasets

### Aims to
#		- load the classification keys (excel sheets of corrected species names) and the v6 datasets
#		- use them to correct the species names in the v6 datasets
#		- check resulting labels/ species names
#		- fill in the potential gaps at the higher taxonomic levels (order? family ?)

module load R/3.4.3/ # To load latest R version on kryo

### Latest update: 05/06/2018

library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")


### ----------------------------------------------------------------------------------------------------------------------------


##### 1) v6-v5.1v3.1  ----------------------------------------------------------------------------------------------------------

### First, load the classif file containing the corrected species labels
setwd("/UP_home/fabioben/Desktop/OVERSEE/data")
names_main <- read.csv("species_v6-v5.1v3.2.csv", h = TRUE, sep = ";")
dim(names_main)
str(names_main)
head(names_main)
# Plus the v6-v5.1v3.1 species
names <- read.csv("species_v6-v5.1v3.1.csv", h = TRUE, sep = ";")
dim(names)
str(names)
head(names)

# rbind both
names <- rbind(names_main, names)
rm(names_main)

# unique(names$correct_name)

### Second, load the observation data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
# dir()
# Identify the group files to clean, apply correction and save in v7 directory per file
files <- dir()[c(1:3,5,7:19,21:22)]
# files
### For each file:
#	- remove the obs that correspond to species names that are labelled as "to remove"
#	- correct the labels that are labelled as 'to correct
#	- save in v7 dir, check results

# f <- files[3]
for(f in files) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.1/")
		message(paste(f, sep = ""))
		
		# Load the data
		data <- get(load(f))
		# Add underscore to species names if you have not done so already
		data$species <- str_replace_all(as.character(data$species), " ", "_")
		# unique(data$species)
		# class(data$species)
		
		# Remove the species that are marked as 'remove'
		toremove <- as.character(unique(names[which(names$action == "remove"),"species"]))
		# class(toremove)
		data2 <- data[!(data$species %in% toremove),]
		# dim(data) ; dim(data2) 
		
		# And correct species labels when necessary
		tocorrect <- unique(names[which(names$action == "correct"),"correct_name"])
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
		# dim(data3)
		# str(data3)
		# head(data3)
		# unique(data3$species)
		
		### Save in proper v7 dir
		message(paste("------------------------------------------------------------------------------------------", sep = ""))
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.1/")
		save(data3, file = str_replace(f, "15_05_18", "05_06_18") )
	
} # eo for loop


### Check v7 results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.1/")
files <- dir()
res <- lapply(files, function(f) {
			dd <- get(load(f))
			return(dd)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table) # 1'532'615
str(table)
length(unique(table$species)) # 1973 species
rm(res,data3,data2,data,tocorrect,toremove,wrongnames) ; gc()

### Check if they match with the correct_names from names
unique(table$species)
dplyr::setdiff(unique(table$species), names$correct_name)
table[table$species == "Caesaromysis_hispida",]
names[names$correct_name == "Caesaromysis_hispida",]
# OK

unique(table[table$class == "Scyphozoa","species"])
# OK only the 3 holoplanktonic species :)
unique(table[table$class == "Hexanauplia","species"])

# count species and order per n
data.frame(table[table$class == "Hexanauplia",] %>% count(species))
# OK gut !!
rm(table)
gc()

### You can keep these v7 datasets, just copy/paste the two PANGAEA datasets





##### 2) v6-v5.1v3.2  ----------------------------------------------------------------------------------------------------------
### First, load the classif file containing the corrected species labels
setwd("/UP_home/fabioben/Desktop/OVERSEE/data")
names <- read.csv("species_v6-v5.1v3.2.csv", h = TRUE, sep = ";")
dim(names)
str(names)
head(names)


### Second, load the observation data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
# dir()
# Identify the group files to clean, apply correction and save in v7 directory per file
files <- dir()[c(1:3,5,7:19,21:22)]
# files
### For each file:
#	- remove the obs that correspond to species names that are labelled as "to remove"
#	- correct the labels that are labelled as 'to correct
#	- save in v7 dir, check results

# f <- files[2]
for(f in files) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.1v3.2/")
		message(paste(f, sep = ""))
		
		# Load the data
		data <- get(load(f))
		# Add underscore to species names if you have not done so already
		data$species <- str_replace_all(as.character(data$species), " ", "_")
		
		# Remove the species that are marked as 'remove'
		toremove <- as.character(unique(names[which(names$action == "remove"),"species"]))
		# class(toremove)
		data2 <- data[!(data$species %in% toremove),]
		# dim(data) ; dim(data2) 
		
		# And correct species labels when necessary
		tocorrect <- unique(names[which(names$action == "correct"),"correct_name"])
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
		# dim(data3)
		# str(data3)
		# head(data3)
		# unique(data3$species)
		
		### Save in proper v7 dir
		message(paste("------------------------------------------------------------------------------------------", sep = ""))
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.2/")
		save(data3, file = str_replace(f, "15_05_18", "05_06_18") )
	
} # eo for loop


### Check v7 results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.1v3.2/")
files <- dir()
res <- lapply(files, function(f) {
			dd <- get(load(f))
			return(dd)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table) # 2'393'495
str(table)
length(unique(table$species)) # 2211 species
# Clean 
rm(res, data3, data2, data, tocorrect, toremove, wrongnames) ; gc()

### Check if they match with the correct_names from names
unique(table$species)
dplyr::setdiff(unique(table$species), names$correct_name)
table[table$species == "Caesaromysis_hispida",]
names[names$correct_name == "Caesaromysis_hispida",]
# OK

unique(table[table$class == "Scyphozoa","species"])
# OK only the 3 holoplanktonic species :)
unique(table[table$class == "Hexanauplia","species"])

# count species and order per n
data.frame(table[table$class == "Hexanauplia",] %>% count(species))
# OK gut !!

rm(table)
gc()




##### 3) v6-v5.2v3.1  ----------------------------------------------------------------------------------------------------------

### First, load the classif file containing the corrected species labels
setwd("/UP_home/fabioben/Desktop/OVERSEE/data")
names_main <- read.csv("species_v6-v5.1v3.2.csv", h = TRUE, sep = ";")
dim(names_main)
str(names_main)
head(names_main)
# Plus the v6-v5.1v3.1 species
names <- read.csv("species_v6-v5.2v3.1.csv", h = TRUE, sep = ";")
dim(names)
str(names)
head(names)

# rbind both
names <- rbind(names_main, names)
rm(names_main)

### Second, load the observation data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
# dir()
# Identify the group files to clean, apply correction and save in v7 directory per file
files <- dir()[c(1:3,5,7:19,21:22)]
# files
### For each file:
#	- remove the obs that correspond to species names that are labelled as "to remove"
#	- correct the labels that are labelled as 'to correct
#	- save in v7 dir, check results

# f <- files[2]
for(f in files) {
		
		# Load the data
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.1/")
		data <- get(load(f))
		# Add underscore to species names if you have not done so already
		data$species <- str_replace_all(as.character(data$species), " ", "_")
	
		# Remove the species that are marked as 'remove'
		toremove <- as.character(unique(names[which(names$action == "remove"),"species"]))
		# class(toremove)
		data2 <- data[!(data$species %in% toremove),]
		# dim(data) ; dim(data2) 
		
		# And correct species labels when necessary
		tocorrect <- unique(names[which(names$action == "correct"),"correct_name"])
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
		# dim(data3)
		# str(data3)
		# head(data3)
		# unique(data3$species)
		
		### Save in proper v7 dir
		message(paste("------------------------------------------------------------------------------------------", sep = ""))
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.1/")
		save(data3, file = str_replace(f, "15_05_18", "05_06_18") )
	
} # eo for loop


### Check v7 results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.1/")
files <- dir()
res <- lapply(files, function(f) {
			dd <- get(load(f))
			return(dd)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table) # 1'072'102
str(table)
length(unique(table$species)) # 1313 species
rm(res, data3, data2, data, tocorrect, toremove, wrongnames) ; gc()

### Check if they match with the correct_names from names
unique(table$species)
dplyr::setdiff(unique(table$species), names$correct_name)
table[table$species == "Caesaromysis_hispida",]
names[names$correct_name == "Caesaromysis_hispida",]
# OK

unique(table[table$class == "Scyphozoa","species"])
# OK only the 3 holoplanktonic species :)
unique(table[table$class == "Hexanauplia","species"])

# count species and order per n
data.frame(table[table$class == "Hexanauplia",] %>% count(species))
# OK gut !!

rm(table)
gc()



##### 4) v6-v5.2v3.2  ----------------------------------------------------------------------------------------------------------

### First, load the classif file containing the corrected species labels
setwd("/UP_home/fabioben/Desktop/OVERSEE/data")
names_main <- read.csv("species_v6-v5.1v3.2.csv", h = TRUE, sep = ";")
dim(names_main)
str(names_main)
head(names_main)
# Plus the v6-v5.1v3.1 species
names <- read.csv("species_v6-v5.2v3.2.csv", h = TRUE, sep = ";")
dim(names)
str(names)
head(names)

# rbind both
names <- rbind(names_main, names)
rm(names_main)

### Second, load the observation data
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
# dir()
# Identify the group files to clean, apply correction and save in v7 directory per file
files <- dir()[c(1:3,5,7:19,21:22)]
# files
### For each file:
#	- remove the obs that correspond to species names that are labelled as "to remove"
#	- correct the labels that are labelled as 'to correct
#	- save in v7 dir, check results

# f <- files[2]
for(f in files) {
		
		# Useless message
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
		message(paste(f, sep = ""))
		
		# Load the data
		data <- get(load(f))
		# Add underscore to species names if you have not done so already
		data$species <- str_replace_all(as.character(data$species), " ", "_")
		
		# Remove the species that are marked as 'remove'
		toremove <- as.character(unique(names[which(names$action == "remove"),"species"]))
		# class(toremove)
		data2 <- data[!(data$species %in% toremove),]
		# dim(data) ; dim(data2) 
		
		# And correct species labels when necessary
		tocorrect <- unique(names[which(names$action == "correct"),"correct_name"])
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
		# dim(data3)
		# str(data3)
		# head(data3)
		# unique(data3$species)
		
		### Save in proper v7 dir
		message(paste("------------------------------------------------------------------------------------------", sep = ""))
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
		save(data3, file = str_replace(f, "15_05_18", "05_06_18") )
	
} # eo for loop


### Check v7 results
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
files <- dir()
res <- lapply(files, function(f) {
			dd <- get(load(f))
			return(dd)
}
) # eo lapply
table <- do.call(rbind, res)
dim(table) # 1'864'450 obs
str(table)
length(unique(table$species)) # 1663 species
rm(res, data3, data2, data, tocorrect, toremove, wrongnames) ; gc()

### Check if they match with the correct_names from names
unique(table$species)
dplyr::setdiff(unique(table$species), names$correct_name)
table[table$species == "Heteromysis_panamaensis",]
names[names$correct_name == "Heteromysis_panamaensis",]
# OK

unique(table[table$class == "Scyphozoa","species"])
# OK only the 3 holoplanktonic species :)
unique(table[table$class == "Hexanauplia","species"])

# count species and order per n
data.frame(table[table$class == "Hexanauplia",] %>% count(species))
# OK gut !!

rm(table)
gc()


##### ----------------------------------------------------------------------------------------------------------

### For each v7 dataset, report n obs and n species onto the workflow excel sheet
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
files <- dir()
for(f in files) {
		message(paste(f, sep = ""))
		d <- get(load(f))
		message(paste("n obs = ", nrow(d), sep = ""))
		message(paste("n species  = ", length(unique(d$species)), sep = ""))
		message(paste("", sep = ""))
		message(paste("", sep = ""))
		rm(d)
		gc()
} # eo for loop


##### ----------------------------------------------------------------------------------------------------------

### Finally, just c/p the 2 PANGAEA datasets in the v7 directories
# Load
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
d1 <- get(load("Copepoda_PANGAEA_04_06_18.Rdata"))
d2 <- get(load("Thecosomata_MAREDAT_31_05_18.Rdata"))
# Save in the v7 dir
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
save(d1, file = "Copepoda_PANGAEA_05_06_18.Rdata")
save(d2, file = "Thecosomata_MAREDAT_05_06_18.Rdata")


setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v7-5.2v3.2/")
d1 <- get(load("Copepoda_PANGAEA_05_06_18.Rdata"))
nrow(d1)
unique(d1$species)

d2 <- get(load("Thecosomata_MAREDAT_05_06_18.Rdata"))
nrow(d2)
unique(d2$species)

### And gather all datasets to report total nb of occurrences and species number
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v6-v5.2v3.2/")
dir()
res <- lapply(dir(), function(f) {
			dd <- get(load(f))
			return( dd[,c("x","y","species")] )
}
) # eo lapply
table <- do.call(rbind, res)
dim(table)
str(table)
head(table)
length(unique(table$species))











