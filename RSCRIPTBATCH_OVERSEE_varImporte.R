
### ==============================================================================================================================

library("raster")
library("sp")
library("dplyr")
library("stringr")
library("reshape2")
library("tidyverse")
library("caret")
library("MuMIn")
#library("randomUniformForest")
#library("randomForest")
library("lmerTest")

### ==============================================================================================================================

WD <- getwd()

setwd(paste(WD,"/","species_data_v9v3.1","/","total_bckgrnd","/", sep = ""))
files <- dir()[grep("data_total_", dir())]
 
# 02/10/2018: To identify the species for which variable importance tabels have already been printed
 species <- str_replace_all(files, ".txt", "")
 species <- str_replace_all(species, "data_total_", "")
 setwd( paste(WD,"/","species_data_v9v3.1","/","total_bckgrnd","/","var_import_tests3","/", sep = "") )
 donespp <- dir()[grep("_test1_", dir())]
 donespp <- str_replace_all(donespp, ".Rdata", "")
 donespp <- str_replace_all(donespp, "var_import_test1_", "")
 species2todo <- species[!(species %in% donespp)]
 files <- paste("data_total_",species2todo,".txt", sep = "")

setwd(paste(WD,"/","species_data_v9v3.1","/","total_bckgrnd","/", sep = ""))

#f <- files[2]

vaRanker <- function(f = files) {		
				
				# Load the species dataset
				data <- read.table(f, h = T, sep = ";")
				sp <- unique(data[data$obs == 1,"species"])
				phyl <- unique(data[data$obs == 1,"phylum"])
				message(paste("Ranking variables importance for ",sp, " =================", sep = ""))
				message(paste("", sep = ""))
				
				# Strat 1: fit GLM and extract ranking based on t stats
				fit1 <- glm(formula = obs ~ SST+SSS+MLD1+logChl+deltaT_1d+logSiO2+PAR, weights = weights, family = "binomial", data = data )
				require("caret")
				rank1 <- caret::varImp(fit1)
				rank1 <- rank1/ max(rank1)
				# And make a table out of it
				table1 <- data.frame(var = rownames(rank1), rank = rank1)
				table1$species <- sp
				table1$phyl <- phyl
				colnames(table1)[2] <- "rank"
				rm(rank1, fit1)
				table1 <- table1[order(table1$rank, decreasing = T),]
				
				### ================================================
				# Strat 2: Sum of squares from lmer, same formula as above but out y as random effect to account for repetitivity
				# require("lme4")
# 				rand_fit <- lme4::glmer(formula = as.numeric(obs) ~ SST+SSS+MLD1+logChl+deltaT_1d+logSiO2+PAR + (1|y),
# 										weights = weights, family = "binomial", data = data)  # takes a bit longer than the other
#
# 				table2 <- data.frame(anova(rand_fit))
# 				colnames(table2) <- c("Df","SSQ","MeanSQ","F")
#
# 				# And compute ranking (has to span the 0-1 range)
# 				max <- max(table2$SSQ)
# 				table2$rank <- table2$SSQ/ max
# 				table2$species <- sp
# 				table2$phyl <- phyl
# 				table2$var <- rownames(table2)
# 				rm(rand_fit, max)
# 				table2 <- table2[order(table2$rank, decreasing = T),]

				### ================================================
				# Strat 3: Extract variable importance from Random Forests
				data4rf <- na.omit(data[,c(8,40,41,36:38,26,27,30,31,19,21,22)]) 
				# Fit RF with 1000 trees
				require("randomForest")
				# ?randomForest
				RF <- randomForest(formula = obs ~ SST+SSS+MLD1+logChl+deltaT_1d+logSiO2+PAR, ntree = 999, data = data4rf)
				#RF <- randomUniformForest(Y = data4rf$obs, X = data4rf[,c(4:13)], ntree = 999, BreimanBounds = F )
				# plot(RF)
				# importance(RF)
				#table3 <- data.frame(RF$forest$variableImportance)
				table3 <- data.frame(var = rownames(importance(RF)), IncNodePurity = importance(RF))
				#colnames(table3) <- c("var","score","per","import")
				# And add a 'rank' ranging between 1 and 0
				table3$rank <- table3$IncNodePurity / max(table3$IncNodePurity) # or simply 100...
				table3$species <- sp
				table3$phyl <- phyl
				rownames(table3) <- table3$var
				table3 <- table3[order(table3$rank, decreasing = T),]
				rm(RF)
				
				### Go save these tables in a dir
				setwd( paste(WD,"/","species_data_v9v3.1","/","total_bckgrnd","/","var_import_tests3","/", sep = "") )
				save(table1, file = paste("var_import_test1_", sp, ".Rdata", sep = "") )
				#save(table2, file = paste("var_import_test2_", sp, ".Rdata", sep = "") )
				save(table3, file = paste("var_import_test3_", sp, ".Rdata", sep = "") )
				
				setwd(paste(WD,"/","species_data_v9v3.1","/","total_bckgrnd","/", sep = ""))
				rm(table1, table3) ; gc()
				
} # eo FUN


library("parallel")
mclapply(X = files, FUN = vaRanker, mc.cores = 20)



