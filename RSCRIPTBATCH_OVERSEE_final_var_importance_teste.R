
##### 20/03/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, D-USYS, ETH Zürich
##### Script for : 
#	- Defining a list containing the names and variables of all the variables pools defined
#	- Performing the two tests of variable importances (GLMs and RF) and extact normalized ranking of vars in them
#	- Extracting adjusted R2 of the models too as a complementary test

### Last update : 25/03/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("tidyverse")
library("reshape2")
library("caret")
library("lmerTest")

WD <- getwd()

### Define function for ranking variables 
vaRankerPHYTO <- function(f = files) {
	
				# Load the species dataset
				setwd(second.wd)
				data <- read.table(f, h = T, sep = ";")
				sp <- unique(data[1,"species"])
				grp <- unique(data[1,"group"])
				d <- data[,c(pool,"obs")]
				message(paste("Ranking variables importance for ",sp, " =================", sep = ""))
				message(paste("", sep = ""))
		
				### ================================================
				# Strat 1: linear model with forward selection of variables
				#colnames(data)[21] <- "dSST"
				lm <- lm(obs ~ ., data = d, weights = data$weights)
				rsq_lm <- summary(lm)$adj.r.squared
				ranklm <- varImp(lm)
				ranklm$norm <- ranklm$Overall/ max(ranklm$Overall)
				ranklm$species <- str_replace_all(sp, " ", "_")
				ranklm$grouping <- grp
				ranklm$var <- rownames(ranklm)
				colnames(ranklm)[1] <- "rank"
				rm(lm)
				ranklm <- ranklm[order(ranklm$rank, decreasing = T),]
				ranklm$R2 <- rsq_lm
				ranklm$test <- "test1"
				ranklm$pool <- n
	
				### ================================================
				# Strat 2: fast random forest
				require("ranger")
				rf <- ranger(obs ~ ., data = na.omit(d), num.trees = 700, min.node.size = 10, importance = "impurity")
				#rf$variable.importance
				rsq_rf <- rf$r.squared
				rankRF <- data.frame(var = names(rf$variable.importance), rank = rf$variable.importance)
				rankRF$norm <- rankRF$rank/ max(rankRF$rank)
				rankRF$species <- str_replace_all(sp, " ", "_")
				rankRF$grouping <- grp
				rm(rf)
				rankRF <- rankRF[order(rankRF$rank, decreasing = T),]
				rankRF$R2 <- rsq_rf
				rankRF$test <- "test2"
				rankRF$pool <- n
				
				### ================================================
				### Go save these tables in a dir
				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/varimport_tests_25_03_20/phyto")
				save(ranklm, file = paste("var_import_phyto_test1_pool",n,"_",str_replace_all(sp," ","_"),".Rdata", sep = "") )
				save(rankRF, file = paste("var_import_phyto_test2_pool",n,"_",str_replace_all(sp," ","_"),".Rdata", sep = "") )
                #return(list(ranklm,rankRF))
				rm(sp, d, rsq_rf, rsq_lm)
	
} # eo vaRankerPHYTO

vaRankerZOO <- function(f = files) {
	
				# Load the species dataset
                # f <- files[11]
				data <- read.table(f, h = T, sep = ";")
                
				### Provide the taxonomic grouping values you use for mapping richness patterns
				#data$grouping <- NA
				# Condition loop to fill the new grouping vector
				if( unique(data[1,"class"]) %in% c("Hexanauplia","Maxillopoda") ) {
						grp <- "Copepoda"
				} else if (unique(data[1,"phylum"]) == "Chaetognatha") {
						grp <- "Chaetognatha"
				} else if (unique(data[1,"phylum"]) == "Mollusca") {
						grp <- "Pteropoda"
				} else if (unique(data[1,"phylum"]) %in% c("Cnidaria","Ctenophora") ) {
						grp <- "Jellyfish"
				} else if (unique(data[1,"phylum"]) == "Chordata") {
						grp <- "Chordata"
				} else if (unique(data[1,"phylum"]) == "Annelida") {
						grp <- "Annelida"
				} else if (unique(data[1,"class"]) %in% c("Ostracoda","Branchiopoda") ) {
						grp <- "Other_arthropoda"
				} else if (unique(data[1,"class"]) == "Malacostraca") {
						grp <- "Malacostraca"
				} else if (unique(data[1,"phylum"]) == "Foraminifera") {
						grp <- "Foraminifera"
				} # eo else if loop
				
				sp <- unique(data[1,"species"])
				phyl <- unique(data[1,"phylum"])
				#grp <- unique(data$grouping)
				message(paste("Ranking variables importance for ",sp, " =================", sep = ""))
				message(paste("", sep = ""))
		
				d <- data[,c(pool,"obs")]
		
				### ================================================
				# Strat 1: linear model with forward selection of variables
				# colnames(data)[21] <- "dSST"
				lm <- lm(obs ~ ., data = d, weights = data$weights)
				rsq_lm <- summary(lm)$adj.r.squared
				ranklm <- varImp(lm)
				ranklm$norm <- ranklm$Overall/ max(ranklm$Overall)
				ranklm$species <- sp
				ranklm$phyl <- phyl
				ranklm$grouping <- grp
				ranklm$var <- rownames(ranklm)
				colnames(ranklm)[1] <- "rank"
				rm(lm)
				ranklm <- ranklm[order(ranklm$rank, decreasing = T),]
				ranklm$R2 <- rsq_lm
				ranklm$test <- "test1"
				ranklm$pool <- n
	
				### ================================================
				# Strat 2: fast random forest
				require("ranger")
				rf <- ranger(obs ~ ., data = na.omit(d), num.trees = 700, min.node.size = 10, importance = "impurity")
				rsq_rf <- rf$r.squared
				rankRF <- data.frame(var = names(rf$variable.importance), rank = rf$variable.importance)
				rankRF$norm <- rankRF$rank/ max(rankRF$rank)
				rankRF$species <- sp
				rankRF$phyl <- phyl
				rankRF$grouping <- grp
				rm(rf)
				rankRF <- rankRF[order(rankRF$rank, decreasing = T),]
				rankRF$R2 <- rsq_rf
				rankRF$test <- "test2"
				rankRF$pool <- n
		
				### Go save these tables in a dir
				setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/varimport_tests_25_03_20/zoo")
				save(ranklm, file = paste("var_import_zoo_test1_pool",n,"_",str_replace_all(sp," ","_"),".Rdata", sep = "") )
				save(rankRF, file = paste("var_import_zoo_test2_pool",n,"_",str_replace_all(sp," ","_"),".Rdata", sep = "") )
                #test <- rbind(ranklm,rankRF) 
                #return(test)
				#rm(phyl, sp, d, rankRF, ranklm)
				setwd(second.wd)
	
} # eo vaRankerZOO


# --------------------------------------------------------------------------------------------------------------------------------

### A°) Phytoplankton
setwd(paste(WD,"/phytoplankton_15_01_19/total_background/species_data/",sep=""))
second.wd <- getwd()
setwd(paste(second.wd,"/var_import_tests_03_19/",sep=""))
third.wd <- getwd()
setwd(second.wd)

# Define list of variable pools
pool1 <- c("SST","dSST","logNO3","logChl","logSiO2","logEKE","Nstar","Sistar","Wind")
pool2 <- c("SST","dSST","logNO3","logChl","logSiO2","logEKE","Nstar","Sistar","MLPAR1")
pool3 <- c("SST","dSST","logNO3","logChl","logSiO2","logEKE","Nstar","Sistar","MLD1","PAR")
### p1-3 : tests the effect of WindxMLPAR1x"MLD1","PAR"
pool4 <- c("SST","dSST","logNO3","logChl","logSiO2","logEKE","Nstar","PAR")
### p3-4: tests the effect of Sistar vs no Sistar
pool5 <- c("SST","dSST","logNO3","logChl","logEKE","Nstar","PAR")
pool6 <- c("SST","dSST","logChl","logSiO2","logEKE","Nstar","PAR")
### p5-6 -> tests the effect of logNO3xlogSiO2
### p4 vs p5&p6 --> is logNO3+logSiO2 better than logSiO2 OR logNO3 ? 
pool7 <- c("SST","dSST","logChl","logSiO2","Nstar","PAR")
### p6 vs p7 -> tests the effect of logEKE
pool8 <- c("SST","dSST","logChl","logSiO2","logEKE","PAR")
pool9 <- c("SST","dSST","logChl","logSiO2","PAR")
### p8 vs p9 -> tests the effetc of logEKE
### p7 vs p9 -> tests the effect of Nstar
### p6 vs p8 -> tests the effect of Nstar
### p7 vs p8 -> tests Nstar vs logEKE

### And then probably test the impact of each 'essential' variables individually (just to make sure they do decrease the r2)
# Start from the test 9/ pool7 pool: SST, dSST, Nstar, logSiO2, logChl, PAR  and compare it to:
# p10: remove SST
pool10 <- c("dSST","logChl","logSiO2","Nstar","PAR")
# p11: remove dSST
pool11 <- c("SST","logChl","logSiO2","Nstar","PAR")
# p12: remove logSiO2
pool12 <- c("SST","dSST","logChl","Nstar","PAR")
# p13: remove logChl
pool13 <- c("SST","dSST","logSiO2","Nstar","PAR")
# p14: remove PAR
pool14 <- c("SST","dSST","logChl","logSiO2","Nstar")


# ...
list_pools_phyto <- list(pool1,pool2,pool3,pool4,pool5,pool6,pool7,pool8,pool9,pool10,pool11,pool12,pool13,pool14)
names(list_pools_phyto) <- c("pool1","pool2","pool3","pool4","pool5","pool6","pool7","pool8","pool9","pool10","pool11","pool12","pool13","pool14")

### For each pool of vars in "list_pools_phyto", perform tests and save
#n <- 7

for(n in c(10:14) ) {
		
		message(paste("Performing tests for phytoplankton - pool ||  ", n, sep = ""))
		
		# Get corresponding pool and its label
		pool <- list_pools_phyto[[n]]
		label <- names(list_pools_phyto)[n]
		
		# Get files names in phyto.wd
		files <- dir()[grep("data_total", dir())]
		
		require("parallel")
		mclapply(X = files, FUN = vaRankerPHYTO, mc.cores = 25)
	
} # eo first for loop

### Go to third.wd and rbind all .Rdata files in a table 
setwd(third.wd)
files <- dir()
res <- mclapply(files, function(f) {
				data <- get(load(f))
				return(data)
			}, mc.cores = 25
) # eo lapply
# Rbind results
data <- do.call(rbind, res)
data$var <- factor(data$var)
str(data)
rm(res)

summary(data)
unique(data$pool)

# Save as .txt file and make plots from your comp
setwd(WD)
write.table(data, file = "poolz_skillz_table_phyto_tot_26_03_19.txt", sep = ";")

# --------------------------------------------------------------


### B°) Zooplankton
setwd(paste(WD,"/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/",sep=""))
second.wd <- getwd()
setwd(paste(second.wd,"/var_import_tests_03_19/",sep=""))
third.wd <- getwd()
setwd(second.wd)

# Define list of variable pools
pool1 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE","Wind","Nstar","Sistar")
pool2 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE","MLPAR","Nstar","Sistar")
pool3 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE","MLD1","PAR","Nstar","Sistar")
### p1-3: tests the effect of WindxMLPARxMLD+PAR

pool4 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE","Nstar","PAR")
pool5 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE","Nstar","MLD1")
### p4-5: tests the effect of PARvsMLD

### p3 vs p4 and p5 --> tets the effect of Sistar

pool6 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE","Nstar")
pool7 <- c("SST","dSST","dO2","logChl","logNO3","logEKE","Nstar")
### p6 vs p4 --> test the impact of PAR
### p5 vs p6 --> test the impact of MLD
### p6-7 test the effect of logSiO2xlogNO3

### And then probably test the impact fo removing logEKE and Nstar ? 
pool8 <- c("SST","dSST","dO2","logChl","logSiO2","Nstar")
pool9 <- c("SST","dSST","dO2","logChl","logSiO2","logEKE")
pool10 <- c("SST","dSST","dO2","logChl","logSiO2")
### p6 vs p8 --> tests the effect of logEKE
### p6 vs p9 --> tests the effect of Nstar
### p8 vs p10 --> tests the effect of Nstar too
### p9 vs p10 --> tests the effect of logEKE too

### And then probably test the impact of each 'essential' variables individually (just to make sure they do decrease the r2)
# Start from the test 14/ pool9 : SST, dSST, dO2, logSiO2, logChl, logEKE and compare it to:
# p11: remove SST
pool11 <- c("dSST","dO2","logChl","logSiO2","logEKE")
# p12: remove dSST
pool12 <- c("SST","dO2","logChl","logSiO2","logEKE")
# p13: remove dO2
pool13 <- c("SST","dSST","logChl","logSiO2","logEKE")
# p14: remove logChl
pool14 <- c("SST","dSST","dO2","logSiO2","logEKE")
# p15: remove logSiO2
pool15 <- c("SST","dSST","dO2","logChl","logEKE")

# ...
list_pools_zoo <- list(pool1,pool2,pool3,pool4,pool5,pool6,pool7,pool8,pool9,pool10,pool11,pool12,pool13,pool14,pool15)
names(list_pools_zoo) <- c("pool1","pool2","pool3","pool4","pool5","pool6","pool7","pool8","pool9","pool10","pool11","pool12","pool13","pool14","pool15")

### For each pool of vars in "list_pools_phyto", perform tests and save
# n <- 1

for(n in c(11:15) ) {
		
		message(paste("Performing tests for zooplankton - pool ||  ", n, sep = ""))
		
		# Get corresponding pool and its label
		pool <- list_pools_zoo[[n]]
		label <- names(list_pools_zoo)[n]
		
		# Get files names in phyto.wd
		files <- dir()[grep("data_total", dir())]
		
		require("parallel")
		mclapply(X = files, FUN = vaRankerZOO, mc.cores = 25)
	
} # eo first for loop

### Go to third.wd and rbind all .Rdata files in a table 
setwd(third.wd)
files <- dir()
require("parallel")
res <- mclapply(files, function(f) {
				data <- get(load(f))
				return(data)
			}, mc.cores = 25
) # eo lapply
# Rbind results
data <- do.call(rbind, res)
dim(data)
rm(res)

data$var <- factor(data$var)
data$grouping <- factor(data$grouping)
data$test <- factor(data$test)
# head(data)
summary(data)
unique(data$pool)

# Save as .txt file and make plots from your comp
setwd(WD)
write.table(data, file = "poolz_skillz_table_zoo_tot_26_03_19.txt", sep = ";")


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

### 25/03/2020: For phyto- and zooplankton, check the importance of predictors for the 4 final pools 

### A°) Phytoplankton
setwd(paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/",sep=""))
second.wd <- getwd()

# Define list of variable pools
pool1 <- c("SST","dSST","logNO3","logChl","Nstar","PAR")
pool2 <- c("SST","dSST","logSiO2","logChl","Nstar","PAR")
pool3 <- c("SST","dSST","logNO3","logChl","Sistar","Nstar","PAR")
pool4 <- c("SST","dSST","logNO3","logChl","PAR")
list_pools_phyto <- list(pool1,pool2,pool3,pool4)
names(list_pools_phyto) <- c("pool1","pool2","pool3","pool4")

### For each pool of vars in "list_pools_phyto", perform tests and save
# n <- 1
for(n in c(1:4)) {
    
        message(paste("Performing tests for phytoplankton - pool || ",n, sep = ""))
        message(paste("  ", sep = ""))
	    # Get corresponding pool and its label
	    pool <- list_pools_phyto[[n]]
	    label <- names(list_pools_phyto)[n]
	    
        # Get files names in phyto.wd
        setwd(second.wd)
	    files <- dir()[grep("data_total", dir())]
	    
        require("parallel")
	    mclapply(X = files, FUN = vaRankerPHYTO, mc.cores = 30)
    
} # eo for loop - n 

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/varimport_tests_25_03_20/phyto")
### retrieve results
files <- dir()
require("parallel")
res <- mclapply(files, function(f) {
				data <- get(load(f))
				return(data)
			}, mc.cores = 25
) # eo lapply
# Rbind results
table.phyto <- bind_rows(res)
dim(table.phyto); head(table.phyto)
rm(res)

# Rename vars 
unique(table.phyto$var)
table.phyto$var <- factor(table.phyto$var)
levels(table.phyto$var)[levels(table.phyto$var) == "Nstar"] <- "N*"
levels(table.phyto$var)[levels(table.phyto$var) == "Sistar"] <- "Si*"
# And groupings?
table.phyto$grouping <- factor(table.phyto$grouping)
levels(table.phyto$grouping)
levels(table.phyto$grouping)[levels(table.phyto$grouping) == "bacillariophyceae"] <- "Diatoms"
levels(table.phyto$grouping)[levels(table.phyto$grouping) == "dinoflagellata"] <- "Dinoflagellates"
levels(table.phyto$grouping)[levels(table.phyto$grouping) == "haptophyta"] <- "Haptophytes"

table.phyto <- table.phyto[table.phyto$grouping %in% c("Diatoms","Dinoflagellates","Haptophytes"),]

### Plot (violins + boxplots) distrbution of noramlized rank for LM tests. Also do it for groups by facetting
p1 <- ggplot(data = table.phyto[table.phyto$test == "test1" & table.phyto$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in multiple LM") + 
     theme_classic() 
#
p2 <- ggplot(data = table.phyto[table.phyto$test == "test1" & table.phyto$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in multiple LM") + 
     theme_classic() + facet_wrap(~factor(table.phyto[table.phyto$test == "test1" & table.phyto$R2 > 0.1,"grouping"]),
             ncol = 2, scales = "fixed")
 
### Same with RF tests
p3 <- ggplot(data = table.phyto[table.phyto$test == "test2" & table.phyto$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in RF") + 
     theme_classic() 
#
p4 <- ggplot(data = table.phyto[table.phyto$test == "test2" & table.phyto$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in RFM") + 
     theme_classic() + facet_wrap(~factor(table.phyto[table.phyto$test == "test2" & table.phyto$R2 > 0.1,"grouping"]),
             ncol = 2, scales = "fixed")
#
setwd(WD)
ggsave(plot = p1, filename = "plot_violins_ranks_test1_phyto.jpg", dpi = 300, width = 7, height = 3)             
ggsave(plot = p2, filename = "plot_violins_ranks_test1_phyto_groups.jpg", dpi = 300, width = 9.5, height = 5.5)                
ggsave(plot = p3, filename = "plot_violins_ranks_test2_phyto.jpg", dpi = 300, width = 7, height = 3)   
ggsave(plot = p4, filename = "plot_violins_ranks_test2_phyto_groups.jpg", dpi = 300, width = 9.5, height = 5.5)   


# ----------------------------------------------------------------

### B°) Zooplankton
setwd(paste(WD,"/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/",sep=""))
second.wd <- getwd()

# Define list of variable pools
pool1 <- c("SST","dSST","logNO3","logChl","dO2")
pool2 <- c("SST","dSST","logSiO2","logChl","dO2")
pool3 <- c("SST","dSST","logSiO2","logChl","dO2","Nstar")
pool4 <- c("SST","dSST","logNO3","logChl","dO2","Sistar")
list_pools_zoo <- list(pool1,pool2,pool3,pool4)
names(list_pools_zoo) <- c("pool1","pool2","pool3","pool4")

### For each pool of vars in "list_pools_phyto", perform tests and save
# n <- 1
for(n in c(1:4)) {
    
        message(paste("Performing tests for zooplankton - pool ||  ", n, sep = ""))
        message(paste("  ", sep = ""))
	    # Get corresponding pool and its label
	    pool <- list_pools_zoo[[n]]
	    label <- names(list_pools_zoo)[n]
	    
        # Get files names in zoo.wd
        setwd(second.wd)
	    files <- dir()[grep("data_total", dir())]
	    
        require("parallel")
	    mclapply(X = files, FUN = vaRankerZOO, mc.cores = 30)
    
} # eo for loop - n 

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/varimport_tests_25_03_20/zoo")
### retrieve results
files <- dir()
require("parallel")
res <- mclapply(files, function(f) {
				data <- get(load(f))
				return(data)
			}, mc.cores = 30
) # eo lapply
# Rbind results
table.zoo <- bind_rows(res)
dim(table.zoo); summary(table.zoo)
rm(res)

# Rename vars 
table.zoo$var <- factor(table.zoo$var)
levels(table.zoo$var)
levels(table.zoo$var)[levels(table.zoo$var) == "Nstar"] <- "N*"
levels(table.zoo$var)[levels(table.zoo$var) == "Sistar"] <- "Si*"
# And groupings?
table.zoo$grouping <- factor(table.zoo$grouping)
levels(table.zoo$grouping)

table.zoo <- table.zoo[!(table.zoo$grouping %in% c("Other_arthropoda","Annelida")),]

setwd(WD)
### Plot (violins + boxplots) distrbution of noramlized rank for LM tests. Also do it for groups by facetting
p1 <- ggplot(data = table.zoo[table.zoo$test == "test1" & table.zoo$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in multiple LM") + 
     theme_classic() 
#
p2 <- ggplot(data = table.zoo[table.zoo$test == "test1" & table.zoo$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in multiple LM") + 
     theme_classic() + facet_wrap(~factor(table.zoo[table.zoo$test == "test1" & table.zoo$R2 > 0.1,"grouping"]),
             ncol = 3, scales = "fixed")
 
### Same with RF tests
p3 <- ggplot(data = table.zoo[table.zoo$test == "test2" & table.zoo$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in RF") + 
     theme_classic() 
#
p4 <- ggplot(data = table.zoo[table.zoo$test == "test2" & table.zoo$R2 > 0.1,],
            aes(x = factor(var), y = norm, fill = factor(var))) + 
     geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
     scale_fill_brewer(name = "", palette = "Paired") + xlab("") + ylab("Normalized rank in RFM") + 
     theme_classic() + facet_wrap(~factor(table.zoo[table.zoo$test == "test2" & table.zoo$R2 > 0.1,"grouping"]),
             ncol = 3, scales = "fixed")
#
setwd(WD)
ggsave(plot = p1, filename = "plot_violins_ranks_test1_zoo.jpg", dpi = 300, width = 7, height = 3)             
ggsave(plot = p2, filename = "plot_violins_ranks_test1_zoo_groups.jpg", dpi = 300, width = 12, height = 8)                
ggsave(plot = p3, filename = "plot_violins_ranks_test2_zoo.jpg", dpi = 300, width = 7, height = 3)   
ggsave(plot = p4, filename = "plot_violins_ranks_test2_zoo_groups.jpg", dpi = 300, width = 12, height = 8)   



