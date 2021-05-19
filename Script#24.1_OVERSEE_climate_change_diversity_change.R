
##### 21/05/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#		- loading the results of the niche models for both zoo and phyto (similar predictors set)
#		- exclude taxa with aveg TSS < 0.3
#		- combine probabilities to compute sum of mean HSI
#		- from each pool annual div, plot ensemble diversity, for the present and the future !
#		- Combine and compute change sin diversity between 2100-2071 and 2000-1971
#		- Also examine changes in composition (Bray-Curtis)
 
### Last update: 22/08/2019

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("dplyr")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")
library("viridis")
library("scales")
library("maps")

world2 <- map_data(map = "world2")

# --------------------------------------------------------------------------------------------------------------------------------


### 1°) Set the working directories, vectors etc.
WD <- getwd()

setwd( paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future_rcp85/", sep = "") )
zoo.wd <- getwd()
setwd( paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future_rcp85/", sep = "") )
phyto.wd <- getwd()
setwd(WD)

# Vector of SDMs
SDMs <- c('GLM','GAM','RF','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
# p <- "p4"

### 2°) Plot distrbution of models'TSS values for zoo and phyto separately, using boxplots per SDM and facet per pool
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(zoo.wd,"/eval_scores_",p,"/", sep = "") )
					zoo_spp <- str_replace_all(dir(), "eval_scores_", "")
					zoo_spp <- str_replace_all(zoo_spp, ".Rdata", "")
					scores <- lapply(zoo_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- dplyr::bind_rows(scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}
) # eo lapply
# Rbind
table_scores_zoo <- do.call(rbind,res)
table_scores_zoo$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_zoo)), pattern = "_", n = 2)))[,1])
table_scores_zoo$kingdom <- "Zooplankton"
dim(table_scores_zoo); head(table_scores_zoo)
summary(table_scores_zoo)
rm(res)
table_scores_zoo[is.na(table_scores_zoo$TSS),]

# Plot distribution of TSS
#bplot <- ggplot(na.omit(table_scores_zoo), aes(x = factor(SDM), y = TSS, fill = factor((SDM)))) + 
			#geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "YlGnBu") + theme_bw() + xlab("") + 
			#scale_y_continuous(limits = c(0,1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			#facet_grid(~ factor(pool) )
setwd(WD)
#ggsave(plot = bplot, filename = "boxplot_scores_zoo_21_05_19.jpg", dpi = 300, height = 6, width = 10)


### And for phyto
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(phyto.wd,"/eval_scores_",p,"/", sep = "") )
					phyto_spp <- str_replace_all(dir(), "eval_scores_", "")
					phyto_spp <- str_replace_all(phyto_spp, ".Rdata", "")
					scores <- lapply(phyto_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind, scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}
) # eo lapply
# Rbind
table_scores_phyto <- do.call(rbind,res)
table_scores_phyto$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_phyto)), pattern = "_", n = 2)))[,1])
table_scores_phyto$kingdom <- "Phytoplankton"
dim(table_scores_phyto); head(table_scores_phyto)
summary(table_scores_phyto)
rm(res)
table_scores_phyto[is.na(table_scores_phyto$TSS),]

# Plot distrbution of TSS
#bplot <- ggplot(na.omit(table_scores_phyto), aes(x = factor(SDM), y = TSS, fill = factor((SDM)))) + 
			#geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "YlGnBu") + theme_bw() + xlab("") + 
			#scale_y_continuous(limits = c(0,1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			#facet_grid(~ factor(pool) )
setwd(WD)
#ggsave(plot = bplot, filename = "boxplot_scores_phyto_21_05_19.jpg", dpi = 300, height = 6, width = 10)


### Identify the species to use for ensemble diversity
require("dplyr")
scZ <- data.frame(na.omit(table_scores_zoo) %>%
  		group_by(species) %>%
  		summarise(avg_TSS = mean(TSS), sd_TSS = sd(TSS))
) # eo ddf
scZ[scZ$avg_TSS <= 0.3,] # none, only probel with ONE RUNE of S.magnus
zoo_spp <- unique(scZ$species)
zoo_spp <- zoo_spp[!(zoo_spp %in% c("Spinocalanus_magnus"))]
# For phytoplankton
scP <- data.frame(table_scores_phyto %>%
  		group_by(species) %>%
  		summarise(avg_TSS = mean(TSS), sd_TSS = sd(TSS))
) # eo ddf
phyto_spp <- unique(scP$species)
sp2rm <- scP[scP$avg_TSS <= 0.3,"species"] # generally taxa with presences all over the place
sp2rm <- c(sp2rm,"Actiniscus_pentasterias")
phyto_spp <- phyto_spp[!(phyto_spp %in% sp2rm)]


### 3°) For each kingdom, separately: 
# - load each pool's species HSI, compute sum of mean HSI to obtain monthly diversity estimates
# - compute annual diversity estimate for each pool separately
# - exclude grid cells where extrapolation is too frequent
# - compute ensemble projection of diversity (species richness and Shannon H')
### NOTE: will need 2 lapply (one for pool, one for months) and 1 mclapply (for extracting the species' HSI)
zoo_divs <- lapply(pools, function(p) {

 				message(paste(" ", sep = ""))
 				message(paste("Computing annual diversity for pool || ", p, sep = ""))
 				message(paste(" ", sep = ""))

 				# Extract monthly probabilities
 				res <- lapply(months, function(m)
 				{
 								message(paste(" ", sep = ""))
 								message(paste("Retrieving probabilities for ", m, sep = ""))
 								message(paste(" ", sep = ""))
 								# Load env variables
 								setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
 								env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
 								env <- env[-which(env$SSS < 20),]
 								env <- env[-which(env$Bathy > -175),]
 								# Load HSI
 								message(paste("Loading zoo projections ================================  ", sep = ""))
 								setwd(paste(zoo.wd,"/",p,"/", sep = ""))
 								require("parallel")
 								zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
 								# sp <- "Acartia.Acartia.danae"
 								probas_zoo <- lapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {

 									# Got to species dir
 									setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
 									message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

 									# Need to modify sp when there are 2 names and add brackets
 									if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
 										# Then add brackets around the second piece
 										sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
 												strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
 												strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
 									} # eo if loop

 									# If the 12 monthly projections are done for present & future
 									if( sum(grepl("proj_projection_", dir())) == 24 ) {

 										# Load projections for each SDM
 										setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
 										d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
 										# GAM
 										resModelGam <- d[,"GAM",,]
 										resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
 										resModelGam <- (resModelGam/1000)
 										# GLM
 										resModelGlm <- d[,"GLM",,]
 										resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
 										resModelGlm <- (resModelGlm/1000)
 										# RF
 										resModelRF <- d[,"RF",,]
 										resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
 										resModelRF <- (resModelRF/1000)
 										# ANN
 										resModelANN <- d[,"ANN",,]
 										resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
 										resModelANN <- (resModelANN/1000)

 										# Return
 										return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
 												GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

 									} else {

 										message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

 									} # eo if else loop

 								} ) # eo lapply
 								# cbind SDMs' mean HSI
 								tbl_zoo <- dplyr::bind_rows(probas_zoo)
 								rm(probas_zoo); gc()
 								# Compute average HSI (average across SDMs)
 								tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:length(tbl_zoo))]) )
 								# Compute richness from sum of HSI
 								rich_zoo <- data.frame(tbl_zoo %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
 								colnames(rich_zoo)[4] <- m
 								rm(tbl_zoo)
 								return(rich_zoo)
 							}
 				) # eo lapply
 				for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop

 				# rbind into one table ? (then average per cell_id for average )
 				table.zoo <- dplyr::bind_rows(res)
 				# rm(res)
 				# Compute annual div
 				zoo.div <- data.frame(table.zoo %>% group_by(cell_id) %>%
 						  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
 						  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
 				) # eo ddf
 				rm(table.zoo)
 				# summary(zoo.div)

 				# Return zoo.div
 				zoo.div$pool <- p
 				return(zoo.div)

 		} # eo fun

 ) # eo 1st lapply - pools
# Rbind
zoo_div <- dplyr::bind_rows(zoo_divs) 
dim(zoo_div); head(zoo_div)
summary(zoo_div)
unique(zoo_div$pool)

# Considering how long (>1h) it takes to compute every the diversity estimates, you should save the 'zoo_div' object 
setwd(zoo.wd)
write.table(zoo_div, file = "table_zoo_annual_div_2000-1971.txt", sep = "\t")
rm(zoo_divs)

## For phyto now
phyto_divs <- lapply(pools, function(p) {

 				message(paste(" ", sep = ""))
 				message(paste("Computing annual diversity for pool || ", p, sep = ""))
 				message(paste(" ", sep = ""))

 				# Extract monthly probabilities
 				res <- lapply(months, function(m)
 				{
 								message(paste(" ", sep = ""))
 								message(paste("Retrieving probabilities for ", m, sep = ""))
 								message(paste(" ", sep = ""))
 								# Load env variables
 								setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
 								env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
 								env <- env[-which(env$SSS < 20),]
 								env <- env[-which(env$Bathy > -175),]
 								# Load phyto probas
 								message(paste("Loading zoo projections ================================  ", sep = ""))
 								setwd(paste(phyto.wd,"/",p,"/", sep = ""))
 								require("parallel")
								# sp <- "Tripos.vultur"
 								probas_phyto <- lapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {

 									# Got to species dir
 									setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
 									message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

 									# Need to modify sp when there are 2 names and add brackets
 									if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
 										 # Then add brackets around the second piece
 										sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
 												strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
 												strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
 									} # eo if loop

 									# If the 12 monthly projections are done, for present and future
 									if( sum(grepl("proj_projection_", dir())) == 24 ) {

 										# Load projections for each SDM
 										setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
 										d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
 										# GAM
 										resModelGam <- d[,"GAM",,]
 										resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
 										resModelGam <- (resModelGam/1000)
 										# GLM
 										resModelGlm <- d[,"GLM",,]
 										resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
 										resModelGlm <- (resModelGlm/1000)
 										# RF
 										resModelRF <- d[,"RF",,]
 										resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
 										resModelRF <- (resModelRF/1000)
 										# ANN
 										resModelANN <- d[,"ANN",,]
 										resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
 										resModelANN <- (resModelANN/1000)
 										# Return
 										return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
 												GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

 									} else {

 										message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

 									} # eo if else loop

 								} ) # eo lapply
 								# cbind SDMs' mean HSI
 								tbl_phyto <- dplyr::bind_rows(probas_phyto)
 								rm(probas_phyto); gc()
 								# Compute average HSI (average across SDMs)
 								tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5:length(tbl_phyto))]) )
 								# Compute richness from sum of HSI
 								rich_phyto <- data.frame(tbl_phyto %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
 								colnames(rich_phyto)[4] <- m
 								rm(tbl_phyto)
 								return(rich_phyto)
 							}
 				) # eo lapply
 				for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop
 				# rbind into one table ? (then average per cell_id for average )
 				table.phyto <- dplyr::bind_rows(res)
 				rm(res)
 				# Compute annual div
 				phyto.div <- data.frame(table.phyto %>% group_by(cell_id) %>%
 						  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
 						  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
 				)#  eo ddf
 				rm(table.phyto)
 				# summary(phyto.div)
				
 				phyto.div$pool <- p
 				return(phyto.div)

 		} # eo fun

) # eo 1st lapply - pools
# Rbind
phyto_div <- dplyr::bind_rows(phyto_divs) 
dim(phyto_div); head(phyto_div)
summary(phyto_div)
unique(phyto_div$pool)


setwd(phyto.wd)
write.table(phyto_div, file = "table_phyto_annual_div_2000-1971.txt", sep = "\t")
rm(phyto_divs)
gc()


setwd(WD)

### And now, the same but for future projections based on GFDL-ESM2M
rcp <- "rcp85"

zoo_divs <- lapply(pools, function(p) {

 				message(paste(" ", sep = ""))
 				message(paste("Computing annual diversity for pool || ", p, sep = ""))
 				message(paste(" ", sep = ""))

 				# Extract monthly probabilities
 				res <- lapply(months, function(m) {
 								message(paste(" ", sep = ""))
 								message(paste("Retrieving probabilities for ", m, sep = ""))
 								message(paste(" ", sep = ""))
								if(m == "jan") {
									mm <- "Jan"
								} else if(m == "feb") {
									mm <- "Feb"
								} else if(m == "mar") {
									mm <- "Mar"
								} else if(m == "apr") {
									mm <- "Apr"
								} else if(m == "may") {
									mm <- "May"
								} else if(m == "jun") {
									mm <- "Jun"
								} else if(m == "jul") {
									mm <- "Jul"
								} else if(m == "aug") {
									mm <- "Aug"
								} else if(m == "sep") {
									mm <- "Sep"
								} else if(m == "oct") {
									mm <- "Oct"
								} else if(m == "nov") {
									mm <- "Nov"
								} else if(m == "dec") {
									mm <- "Dec"
								} #
								
 								# Load env variables
 								setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
 								env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
 
 								# Load HSI
 								message(paste("Loading zoo projections ================================  ", sep = ""))
 								setwd(paste(zoo.wd,"/",p,"/", sep = ""))
 								require("parallel")
 								zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
 								# sp <- "Acartia.Acartia.danae"
 								probas_zoo <- lapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {

 									# Got to species dir
 									setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
 									message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

 									# Need to modify sp when there are 2 names and add brackets
 									if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
 										# Then add brackets around the second piece
 										sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
 												strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
 												strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
 									} # eo if loop

 									# If the 12 monthly projections are done for present & future
 									if( sum(grepl("proj_projection_", dir())) == 24 ) {

 										# Load projections for each SDM
 										setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M", sep = ""), sep = "") )
 										d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M_", gsub("\\(|)","",sp), ".RData", sep = "") ))
 										# GAM
 										resModelGam <- d[,"GAM",,]
 										resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
 										resModelGam <- (resModelGam/1000)
 										# GLM
 										resModelGlm <- d[,"GLM",,]
 										resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
 										resModelGlm <- (resModelGlm/1000)
 										# RF
 										resModelRF <- d[,"RF",,]
 										resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
 										resModelRF <- (resModelRF/1000)
 										# ANN
 										resModelANN <- d[,"ANN",,]
 										resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
 										resModelANN <- (resModelANN/1000)

 										# Return
 										return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
 												GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

 									} else {

 										message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

 									} # eo if else loop

 								} ) # eo lapply
 								# cbind SDMs' mean HSI
 								tbl_zoo <- dplyr::bind_rows(probas_zoo)
 								rm(probas_zoo); gc()
 								# Compute average HSI (average across SDMs)
 								tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:length(tbl_zoo))]) )
								# unique(tbl_zoo$species)
 								# Compute richness from sum of HSI
 								rich_zoo <- data.frame(tbl_zoo %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
 								colnames(rich_zoo)[4] <- m
 								rm(tbl_zoo)
 								return(rich_zoo)
 							}
 				) # eo lapply
 				for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop

 				# rbind into one table ? (then average per cell_id for average )
 				table.zoo <- dplyr::bind_rows(res)
 				# rm(res)
 				# Compute annual div
 				zoo.div <- data.frame(table.zoo %>% group_by(cell_id) %>%
 						  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
 						  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
 				) # eo ddf
 				rm(table.zoo)
 				# summary(zoo.div)

 				# Return zoo.div
 				zoo.div$pool <- p
 				return(zoo.div)

 		} # eo fun

) # eo 1st lapply - pools
# Rbind
zoo_div <- dplyr::bind_rows(zoo_divs) 
dim(zoo_div); head(zoo_div)
summary(zoo_div)
unique(zoo_div$pool)

# Considering how long (>1h) it takes to compute every the diversity estimates, you should save the 'zoo_div' object 
setwd(zoo.wd)
write.table(zoo_div, file = "table_zoo_annual_div_2100-2071.txt", sep = "\t")
rm(zoo_divs)
gc()


## For phyto now
phyto_divs <- lapply(pools, function(p) {

 				message(paste(" ", sep = ""))
 				message(paste("Computing annual diversity for pool || ", p, sep = ""))
 				message(paste(" ", sep = ""))

 				# Extract monthly probabilities
 				res <- lapply(months, function(m) {
 								message(paste(" ", sep = ""))
 								message(paste("Retrieving probabilities for ", m, sep = ""))
 								message(paste(" ", sep = ""))
								if(m == "jan") {
									mm <- "Jan"
								} else if(m == "feb") {
									mm <- "Feb"
								} else if(m == "mar") {
									mm <- "Mar"
								} else if(m == "apr") {
									mm <- "Apr"
								} else if(m == "may") {
									mm <- "May"
								} else if(m == "jun") {
									mm <- "Jun"
								} else if(m == "jul") {
									mm <- "Jul"
								} else if(m == "aug") {
									mm <- "Aug"
								} else if(m == "sep") {
									mm <- "Sep"
								} else if(m == "oct") {
									mm <- "Oct"
								} else if(m == "nov") {
									mm <- "Nov"
								} else if(m == "dec") {
									mm <- "Dec"
								} #
								
 								# Load env variables
 								setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
 								env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
 
 								# Load phyto probas
 								message(paste("Loading zoo projections ================================  ", sep = ""))
 								setwd(paste(phyto.wd,"/",p,"/", sep = ""))
 								require("parallel")
								# sp <- "Tripos.vultur"
 								probas_phyto <- lapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {

 									# Got to species dir
 									setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
 									message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

 									# Need to modify sp when there are 2 names and add brackets
 									if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
 										 # Then add brackets around the second piece
 										sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
 												strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
 												strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
 									} # eo if loop

 									# If the 12 monthly projections are done, for present and future
 									if( sum(grepl("proj_projection_", dir())) == 24 ) {

 										# Load projections for each SDM
 										setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M", sep = ""), sep = "") )
 										d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M_", gsub("\\(|)","",sp), ".RData", sep = "") ))
 										# GAM
 										resModelGam <- d[,"GAM",,]
 										resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
 										resModelGam <- (resModelGam/1000)
 										# GLM
 										resModelGlm <- d[,"GLM",,]
 										resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
 										resModelGlm <- (resModelGlm/1000)
 										# RF
 										resModelRF <- d[,"RF",,]
 										resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
 										resModelRF <- (resModelRF/1000)
 										# ANN
 										resModelANN <- d[,"ANN",,]
 										resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
 										resModelANN <- (resModelANN/1000)
										
 										# Return
 										return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
 												GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

 									} else {

 										message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

 									} # eo if else loop

 								} ) # eo lapply
 								# cbind SDMs' mean HSI
 								tbl_phyto <- dplyr::bind_rows(probas_phyto)
 								rm(probas_phyto); gc()
 								# Compute average HSI (average across SDMs)
 								tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5:length(tbl_phyto))]) )
								# unique(tbl_phyto$species)
 								# Compute richness from sum of HSI
 								rich_phyto <- data.frame(tbl_phyto %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
 								colnames(rich_phyto)[4] <- m
 								rm(tbl_phyto)
 								return(rich_phyto)
 							}
 				) # eo lapply
 				for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop

 				# rbind into one table ? (then average per cell_id for average )
 				table.phyto <- dplyr::bind_rows(res)
 				rm(res)
 				# Compute annual div
 				phyto.div <- data.frame(table.phyto %>% group_by(cell_id) %>%
 						  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
 						  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
 				) #  eo ddf
 				rm(table.phyto)
 				# summary(phyto.div)
				
 				phyto.div$pool <- p
 				return(phyto.div)

 		} # eo fun

 ) # eo 1st lapply - pools
# Rbind
phyto_div <- dplyr::bind_rows(phyto_divs) 
dim(phyto_div); head(phyto_div)
summary(phyto_div)
unique(phyto_div$pool)


setwd(phyto.wd)
write.table(phyto_div, file = "table_phyto_annual_div_2100-2071.txt", sep = "\t")
rm(phyto_divs)
gc()


# -------------------------------------------------------

### 12/09/2019: Meanwhile, examine preliminary changes in annual plankton diversity

require("maps")
world2 <- map_data(map = "world2")

p <- "p1"
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

# Extract zoo annual HSI for the baseline period
res <- lapply(months, function(m) {
				message(paste(" ", sep = ""))
				message(paste("Retrieving probabilities for ", m, sep = ""))
				message(paste(" ", sep = ""))
				# Load env variables
				setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
				env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
				env <- env[-which(env$SSS < 20),]
				env <- env[-which(env$Bathy > -175),]
				# Load HSI
				message(paste("Loading zoo projections ================================  ", sep = ""))
				setwd(paste(zoo.wd,"/",p,"/", sep = ""))
				require("parallel")
				zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
				# sp <- "Acartia.Acartia.danae"
				probas_zoo <- mclapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {

					# Got to species dir
					setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
					message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

					# Need to modify sp when there are 2 names and add brackets
					if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
						# Then add brackets around the second piece
						sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
								strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
								strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
					} # eo if loop

					# If the 12 monthly projections are done for present & future
					if( sum(grepl("proj_projection_", dir())) == 24 ) {

						# Load projections for each SDM
						setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
						d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
						# GAM
						resModelGam <- d[,"GAM",,]
						resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
						resModelGam <- (resModelGam/1000)
						# GLM
						resModelGlm <- d[,"GLM",,]
						resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
						resModelGlm <- (resModelGlm/1000)
						# RF
						resModelRF <- d[,"RF",,]
						resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
						resModelRF <- (resModelRF/1000)
						# ANN
						resModelANN <- d[,"ANN",,]
						resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
						resModelANN <- (resModelANN/1000)

						# Return
						return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
								GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

					} else {

						message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

					} # eo if else loop

				}, mc.cores = 25 ) # eo lapply
				# cbind SDMs' mean HSI
				tbl_zoo <- dplyr::bind_rows(probas_zoo)
				rm(probas_zoo); gc()
				# Compute average HSI (average across SDMs)
				tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:length(tbl_zoo))]) )
				# Compute richness from sum of HSI
				rich_zoo <- data.frame(tbl_zoo %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
				colnames(rich_zoo)[4] <- m
				rm(tbl_zoo)
				return(rich_zoo)
			}
) # eo lapply
for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop
# rbind into one table ? (then average per cell_id for average )
table.zoo <- dplyr::bind_rows(res)
#dim(table.zoo)
#head(table.zoo)
# Compute annual div
zoo.div.baseline <- data.frame(table.zoo %>% group_by(cell_id) %>%
		  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
		  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
) # eo ddf
summary(zoo.div.baseline)
rm(table.zoo, res)

# Extract zoo annual HSI for the future period period
rcp <- "rcp85"
res <- lapply(months, function(m) {
				message(paste(" ", sep = ""))
				message(paste("Retrieving probabilities for ", m, sep = ""))
				message(paste(" ", sep = ""))
				if(m == "jan") {
					mm <- "Jan"
				} else if(m == "feb") {
					mm <- "Feb"
				} else if(m == "mar") {
					mm <- "Mar"
				} else if(m == "apr") {
					mm <- "Apr"
				} else if(m == "may") {
					mm <- "May"
				} else if(m == "jun") {
					mm <- "Jun"
				} else if(m == "jul") {
					mm <- "Jul"
				} else if(m == "aug") {
					mm <- "Aug"
				} else if(m == "sep") {
					mm <- "Sep"
				} else if(m == "oct") {
					mm <- "Oct"
				} else if(m == "nov") {
					mm <- "Nov"
				} else if(m == "dec") {
					mm <- "Dec"
				} #
				
				# Load env variables
				setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
				env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")

				# Load HSI
				message(paste("Loading zoo projections ================================  ", sep = ""))
				setwd(paste(zoo.wd,"/",p,"/", sep = ""))
				require("parallel")
				zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
				# sp <- "Acartia.Acartia.danae"
				probas_zoo <- mclapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {

					# Got to species dir
					setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
					message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

					# Need to modify sp when there are 2 names and add brackets
					if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
						# Then add brackets around the second piece
						sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
								strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
								strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
					} # eo if loop

					# If the 12 monthly projections are done for present & future
					if( sum(grepl("proj_projection_", dir())) == 24 ) {

						# Load projections for each SDM
						setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M", sep = ""), sep = "") )
						d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M_", gsub("\\(|)","",sp), ".RData", sep = "") ))
						# GAM
						resModelGam <- d[,"GAM",,]
						resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
						resModelGam <- (resModelGam/1000)
						# GLM
						resModelGlm <- d[,"GLM",,]
						resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
						resModelGlm <- (resModelGlm/1000)
						# RF
						resModelRF <- d[,"RF",,]
						resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
						resModelRF <- (resModelRF/1000)
						# ANN
						resModelANN <- d[,"ANN",,]
						resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
						resModelANN <- (resModelANN/1000)

						# Return
						return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
								GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

					} else {

						message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

					} # eo if else loop

				}, mc.cores = 25 ) # eo lapply
				# cbind SDMs' mean HSI
				tbl_zoo <- dplyr::bind_rows(probas_zoo)
				rm(probas_zoo); gc()
				# Compute average HSI (average across SDMs)
				tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:length(tbl_zoo))]) )
				# unique(tbl_zoo$species)
				# Compute richness from sum of HSI
				rich_zoo <- data.frame(tbl_zoo %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
				colnames(rich_zoo)[4] <- m
				rm(tbl_zoo)
				return(rich_zoo)
				
			} # eo lapply
			
) # eo lapply
for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop
# rbind into one table ? (then average per cell_id for average )
table.zoo <- dplyr::bind_rows(res)
# dim(table.zoo)
# head(table.zoo)
# Compute annual div
zoo.div.fut <- data.frame(table.zoo %>% group_by(cell_id) %>%
		  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
		  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
) # eo ddf
summary(zoo.div.fut)
rm(table.zoo,res)


### Make maps and compute difference
dim(zoo.div.baseline); dim(zoo.div.fut)
# And flip x coordinates for mapping
zoo.div.baseline$x2 <- zoo.div.baseline$x 
zoo.div.baseline[zoo.div.baseline$x < 0 ,"x2"] <- (zoo.div.baseline[zoo.div.baseline$x < 0 ,"x"]) + 360
#
zoo.div.fut$x2 <- zoo.div.fut$x 
zoo.div.fut[zoo.div.fut$x < 0 ,"x2"] <- (zoo.div.fut[zoo.div.fut$x < 0 ,"x"]) + 360

map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = na.omit(zoo.div.baseline)) +
 	scale_fill_viridis(name = "Zooplankton richness\n(baseline)", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = na.omit(zoo.div.fut)) +
 	scale_fill_viridis(name = "Zooplankton richness\n(2071-2100)", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
		
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_zoo_rich_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_zoo_rich_future_rcp85.jpg", dpi = 300, width = 7, height = 5)

### Gut, now calculate differences, but first make sure both ddf are following the same coordinates
summary(zoo.div.baseline)
# x = -179.5 - 179.5
summary(zoo.div.fut)
# x = 0.50 - 359.50
### OK, need to convert the longitudes of obs.clims to re-center them on the Pacific Ocean (== modelled)
zoo.div.baseline$x2 <- zoo.div.baseline$x 
zoo.div.baseline[zoo.div.baseline$x < 0 ,"x2"] <- (zoo.div.baseline[zoo.div.baseline$x < 0 ,"x"]) + 360
# zoo.div.baseline$x2 = zoo.div.fut$x ?

# Change ids
zoo.div.baseline$cell_id <- factor(paste(zoo.div.baseline$x2, zoo.div.baseline$y, sep = "_"))
zoo.div.fut$cell_id <- factor(paste(zoo.div.fut$x, zoo.div.fut$y, sep = "_"))
head(zoo.div.baseline$cell_id); head(zoo.div.fut$cell_id)
# Make them
zoo.div.fut <- zoo.div.fut[order(zoo.div.fut$cell_id),]
zoo.div.baseline <- zoo.div.baseline[order(zoo.div.baseline$cell_id),]
head(zoo.div.baseline); head(zoo.div.fut)


# Restrict to the same cells
# The elements of setdiff(x,y) are those elements in x but not in y.
length(setdiff(zoo.div.fut$cell_id, zoo.div.baseline$cell_id)) # 23247 cells that are not in baseline
length(setdiff(zoo.div.baseline$cell_id, zoo.div.fut$cell_id)) # 0 cells not future 
zoo.div.fut2 <- zoo.div.fut
zoo.div.baseline2 <- zoo.div.baseline
#
zoo.div.fut2 <- zoo.div.fut2[zoo.div.fut2$cell_id %in% unique(zoo.div.baseline2$cell_id),]
zoo.div.baseline2 <- zoo.div.baseline2[zoo.div.baseline2$cell_id %in% unique(zoo.div.fut2$cell_id),]
zoo.div.fut2 <- zoo.div.fut2[order(zoo.div.fut2$cell_id),]
zoo.div.baseline2 <- zoo.div.baseline2[order(zoo.div.baseline2$cell_id),]
# dim(zoo.div.baseline2); dim(zoo.div.fut2)
# summary(zoo.div.baseline2); summary(zoo.div.fut2)
# zoo.div.baseline2[1:100,]
# zoo.div.fut2[1:100,]

diff <- data.frame(x = zoo.div.baseline2$x2, y = zoo.div.baseline2$y, diff = (zoo.div.fut2$mean)-(zoo.div.baseline2$mean) ) # eo ddf
diff$perc <- (diff$diff)/(zoo.div.baseline2$mean)
diff$perc <- (diff$perc)*100
summary(diff)
# Flip x coordinates
diff$x2 <- diff$x 
diff[diff$x < 0 ,"x2"] <- (diff[diff$x < 0 ,"x"]) + 360

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = diff) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = diff) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

setwd(WD)
ggsave(plot = map3, filename = "map_zoo_rich_diff_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map4, filename = "map_zoo_rich_diff_perc_rcp85.jpg", dpi = 300, width = 7, height = 5)



### Ok, gut, now the same for phytoplankton ?
# Extract phyto annual HSI for the baseline period
res <- lapply(months, function(m) {
				message(paste(" ", sep = ""))
				message(paste("Retrieving probabilities for ", m, sep = ""))
				message(paste(" ", sep = ""))
				# Load env variables
				setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
				env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
				env <- env[-which(env$SSS < 20),]
				env <- env[-which(env$Bathy > -175),]
				# Load HSI
				message(paste("Loading zoo projections ================================  ", sep = ""))
				setwd(paste(phyto.wd,"/",p,"/", sep = ""))
				require("parallel")
				probas_phyto <- mclapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {

					# Got to species dir
					setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
					message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

					# Need to modify sp when there are 2 names and add brackets
					if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
						# Then add brackets around the second piece
						sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
								strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
								strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
					} # eo if loop

					# If the 12 monthly projections are done for present & future
					if( sum(grepl("proj_projection_", dir())) == 24 ) {

						# Load projections for each SDM
						setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
						d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
						# GAM
						resModelGam <- d[,"GAM",,]
						resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
						resModelGam <- (resModelGam/1000)
						# GLM
						resModelGlm <- d[,"GLM",,]
						resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
						resModelGlm <- (resModelGlm/1000)
						# RF
						resModelRF <- d[,"RF",,]
						resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
						resModelRF <- (resModelRF/1000)
						# ANN
						resModelANN <- d[,"ANN",,]
						resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
						resModelANN <- (resModelANN/1000)

						# Return
						return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
								GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

					} else {

						message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

					} # eo if else loop

				}, mc.cores = 25 ) # eo lapply
				# cbind SDMs' mean HSI
				tbl_phyto <- dplyr::bind_rows(probas_phyto)
				rm(probas_phyto); gc()
				# Compute average HSI (average across SDMs)
				tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5:length(tbl_phyto))]) )
				# Compute richness from sum of HSI
				rich_phyto <- data.frame(tbl_phyto %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
				colnames(rich_phyto)[4] <- m
				rm(tbl_phyto)
				return(rich_phyto)
			}
) # eo lapply
for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop
# rbind into one table ? (then average per cell_id for average )
table.phyto <- dplyr::bind_rows(res)
# Compute annual div
phyto.div.baseline <- data.frame(table.phyto %>% group_by(cell_id) %>%
		  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
		  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
) # eo ddf
summary(phyto.div.baseline)
rm(table.phyto, res)

# Extract zoo annual HSI for the future period period
rcp <- "rcp85"
res <- lapply(months, function(m) {
				message(paste(" ", sep = ""))
				message(paste("Retrieving probabilities for ", m, sep = ""))
				message(paste(" ", sep = ""))
				if(m == "jan") {
					mm <- "Jan"
				} else if(m == "feb") {
					mm <- "Feb"
				} else if(m == "mar") {
					mm <- "Mar"
				} else if(m == "apr") {
					mm <- "Apr"
				} else if(m == "may") {
					mm <- "May"
				} else if(m == "jun") {
					mm <- "Jun"
				} else if(m == "jul") {
					mm <- "Jul"
				} else if(m == "aug") {
					mm <- "Aug"
				} else if(m == "sep") {
					mm <- "Sep"
				} else if(m == "oct") {
					mm <- "Oct"
				} else if(m == "nov") {
					mm <- "Nov"
				} else if(m == "dec") {
					mm <- "Dec"
				} #
				# Load env variables
				setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
				env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
				# Load HSI
				message(paste("Loading zoo projections ================================  ", sep = ""))
				setwd(paste(phyto.wd,"/",p,"/", sep = ""))
				require("parallel")
				probas_phyto <- mclapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {

					# Got to species dir
					setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
					message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

					# Need to modify sp when there are 2 names and add brackets
					if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
						# Then add brackets around the second piece
						sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
								strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
								strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
					} # eo if loop

					# If the 12 monthly projections are done for present & future
					if( sum(grepl("proj_projection_", dir())) == 24 ) {

						# Load projections for each SDM
						setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M", sep = ""), sep = "") )
						d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M_", gsub("\\(|)","",sp), ".RData", sep = "") ))
						# GAM
						resModelGam <- d[,"GAM",,]
						resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
						resModelGam <- (resModelGam/1000)
						# GLM
						resModelGlm <- d[,"GLM",,]
						resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
						resModelGlm <- (resModelGlm/1000)
						# RF
						resModelRF <- d[,"RF",,]
						resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
						resModelRF <- (resModelRF/1000)
						# ANN
						resModelANN <- d[,"ANN",,]
						resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
						resModelANN <- (resModelANN/1000)

						# Return
						return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
								GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )

					} else {

						message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

					} # eo if else loop

				}, mc.cores = 25 ) # eo lapply
				# cbind SDMs' mean HSI
				tbl_phyto <- dplyr::bind_rows(probas_phyto)
				rm(probas_phyto); gc()
				# Compute average HSI (average across SDMs)
				tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5:length(tbl_phyto))]) )
				# Compute richness from sum of HSI
				rich_phyto <- data.frame(tbl_phyto %>% group_by(cell_id) %>% summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI)) )
				colnames(rich_phyto)[4] <- m
				rm(tbl_phyto)
				return(rich_phyto)
				
			} # eo lapply
			
) # eo lapply
for(i in 1:12) { names(res[[i]])[4] <- "div" } # eo for loop
# rbind into one table ? (then average per cell_id for average )
table.phyto <- dplyr::bind_rows(res)
# Compute annual div
phyto.div.fut <- data.frame(table.phyto %>% group_by(cell_id) %>%
		  summarise(x = unique(x), y = unique(y), cumul = sum(div, na.rm = T),
		  mean = mean(div, na.rm = T), NAs = sum(is.na(div)) )
) # eo ddf
summary(phyto.div.fut)
rm(table.phyto,res)

### Make maps and compute difference
dim(phyto.div.baseline); dim(phyto.div.fut)
# And flip x coordinates for mapping
phyto.div.baseline$x2 <- phyto.div.baseline$x 
phyto.div.baseline[phyto.div.baseline$x < 0 ,"x2"] <- (phyto.div.baseline[phyto.div.baseline$x < 0 ,"x"]) + 360
#
phyto.div.fut$x2 <- phyto.div.fut$x 
phyto.div.fut[phyto.div.fut$x < 0 ,"x2"] <- (phyto.div.fut[phyto.div.fut$x < 0 ,"x"]) + 360

map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = na.omit(phyto.div.baseline)) +
 	scale_fill_viridis(name = "Phytoplankton richness\n(baseline)", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = na.omit(phyto.div.fut)) +
 	scale_fill_viridis(name = "Phytoplankton richness\n(2071-2100)", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
		
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_phyto_rich_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_phyto_rich_future_rcp85.jpg", dpi = 300, width = 7, height = 5)

### Gut, now calculate differences, but first make sure both ddf are following the same coordinates
summary(phyto.div.baseline)
# x = -179.5 - 179.5
summary(phyto.div.fut)
# x = 0.50 - 359.50
### OK, need to convert the longitudes of obs.clims to re-center them on the Pacific Ocean (== modelled)
phyto.div.baseline$x2 <- phyto.div.baseline$x 
phyto.div.baseline[phyto.div.baseline$x < 0 ,"x2"] <- (phyto.div.baseline[phyto.div.baseline$x < 0 ,"x"]) + 360
# zoo.div.baseline$x2 = zoo.div.fut$x ?

# Change ids
phyto.div.baseline$cell_id <- factor(paste(phyto.div.baseline$x2, phyto.div.baseline$y, sep = "_"))
phyto.div.fut$cell_id <- factor(paste(phyto.div.fut$x, phyto.div.fut$y, sep = "_"))
#head(phyto.div.baseline$cell_id); head(phyto.div.fut$cell_id)
# Make them follow the same order
phyto.div.fut <- phyto.div.fut[order(phyto.div.fut$cell_id),]
phyto.div.baseline <- phyto.div.baseline[order(phyto.div.baseline$cell_id),]
head(phyto.div.baseline); head(phyto.div.fut)

# Restrict to the same cells
# The elements of setdiff(x,y) are those elements in x but not in y.
length(setdiff(phyto.div.fut$cell_id, phyto.div.baseline$cell_id)) # 23247 cells that are not in baseline
length(setdiff(phyto.div.baseline$cell_id, phyto.div.fut$cell_id)) # 0 cells not future 
phyto.div.fut2 <- phyto.div.fut
phyto.div.baseline2 <- phyto.div.baseline
#
phyto.div.fut2 <- phyto.div.fut2[phyto.div.fut2$cell_id %in% unique(phyto.div.baseline2$cell_id),]
phyto.div.baseline2 <- phyto.div.baseline2[phyto.div.baseline2$cell_id %in% unique(phyto.div.fut2$cell_id),]
phyto.div.fut2 <- phyto.div.fut2[order(phyto.div.fut2$cell_id),]
phyto.div.baseline2 <- phyto.div.baseline2[order(phyto.div.baseline2$cell_id),]
# dim(phyto.div.baseline2); dim(phyto.div.fut2)
# summary(phyto.div.baseline2); summary(phyto.div.fut2)

diff <- data.frame(x = phyto.div.baseline2$x2, y = phyto.div.baseline2$y, diff = (phyto.div.fut2$mean)-(phyto.div.baseline2$mean) ) # eo ddf
diff$perc <- (diff$diff)/(phyto.div.baseline2$mean)
diff$perc <- (diff$perc)*100
summary(diff)
# Flip x coordinates
diff$x2 <- diff$x 
diff[diff$x < 0 ,"x2"] <- (diff[diff$x < 0 ,"x"]) + 360

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = diff) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = diff) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

setwd(WD)
ggsave(plot = map3, filename = "map_phyto_rich_diff_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map4, filename = "map_phyto_rich_diff_perc_rcp85.jpg", dpi = 300, width = 7, height = 5)



### 13/08/19: Compute and plot total plankton diversity
head(zoo.div.baseline2); head(phyto.div.baseline2)
head(zoo.div.fut2); head(phyto.div.fut2)
# dim(zoo.div.baseline2); dim(phyto.div.baseline2)
# dim(zoo.div.fut2); dim(phyto.div.fut2)

total <- data.frame(x = zoo.div.baseline2$x2, y = zoo.div.baseline2$y, 
	base = (phyto.div.baseline2$mean)+(zoo.div.baseline2$mean),
	fut =  (phyto.div.fut2$mean)+(zoo.div.fut2$mean),
	zoo_base = zoo.div.baseline2$mean, zoo_fut = zoo.div.fut2$mean,
	phyto_base = phyto.div.baseline2$mean, phyto_fut = phyto.div.fut2$mean
) # eo ddf
total$diff <- (total$fut) - (total$base)
total$perc <- (total$diff)/(total$base)
total$perc <- (total$perc)*100
summary(total)

map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = base), data = na.omit(total)) +
 	scale_fill_viridis(name = "Plankton richness\n(baseline)", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = fut), data = na.omit(total)) +
 	scale_fill_viridis(name = "Plankton richness\n(2071-2100)", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = na.omit(total)) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = na.omit(total)) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

# Save
setwd(WD)
ggsave(plot = map1, filename = "map_plankton_rich_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_plankton_rich_fut.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_plankton_rich_diff_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map4, filename = "map_plankton_rich_diff_perc_rcp85.jpg", dpi = 300, width = 7, height = 5)

# And make maps of anomalies between zoo and phyto?
lm <- lm(zoo_base ~ phyto_base, data = total, na.action = na.exclude)
summary(lm) # Adjusted R-squared: R-squared:  0.5107 
#total$fit <- predict(lm)
total$anom_base <- residuals(lm)

lm <- lm(zoo_fut ~ phyto_fut, data = total, na.action = na.exclude)
summary(lm) # Adjusted R-squared: 0.5144 
#ensembles$fit <- predict(lm)
total$anom_fut <- residuals(lm)

summary(total)

map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = anom_base), data = na.omit(total)) +
 	scale_fill_gradient2(name = "Richness anomaly\n(baseline)", low = "#3288bd", high = "#d53e4f", mid = "white",
	limits = c(-200,100)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = anom_fut), data = na.omit(total)) +
 	scale_fill_gradient2(name = "Richness anomaly\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white",
	limits = c(-200,100)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

# Save
setwd(WD)
ggsave(plot = map5, filename = "map_zoo_anom_rich_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map6, filename = "map_zoo_anom_rich_future_rcp85.jpg", dpi = 300, width = 7, height = 5)


# And finally some zonal plots
colnames(total)
p <- ggplot() + geom_point(aes(x = y, y = base), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = base), data = total, colour = "#d53e4f" ) + 
		ylab("Plankton richness\n(baseline)") + xlab("Latitude (°)") + theme_bw() +
		scale_y_continuous(limits = c(0,340)) + coord_flip()
ggsave(plot = p, filename = "p1.jpg", dpi = 300, width = 3, height = 7)

p2 <- ggplot() + geom_point(aes(x = y, y = fut), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = fut), data = total, colour = "#d53e4f" ) + 
		ylab("Plankton richness\n(2100)") + xlab("Latitude (°)") + theme_bw() +
		scale_y_continuous(limits = c(0,340)) + coord_flip()
ggsave(plot = p2, filename = "p2.jpg", dpi = 300, width = 3, height = 7)

p3 <- ggplot() + geom_point(aes(x = y, y = zoo_base), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = zoo_base), data = total, colour = "#d53e4f" ) + 
		ylab("Zooplankton richness\n(baseline)") + xlab("Latitude (°)") + 
		scale_y_continuous(limits = c(0,220)) + theme_bw() + coord_flip()
ggsave(plot = p3, filename = "p3.jpg", dpi = 300, width = 3, height = 7)

p4 <- ggplot() + geom_point(aes(x = y, y = zoo_fut), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = zoo_fut), data = total, colour = "#d53e4f" ) + 
		ylab("Zooplankton richness\n(2100)") + xlab("Latitude (°)") + 
		scale_y_continuous(limits = c(0,220)) + theme_bw() + coord_flip()
ggsave(plot = p4, filename = "p4.jpg", dpi = 300, width = 3, height = 7)

p5 <- ggplot() + geom_point(aes(x = y, y = phyto_base), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = phyto_base), data = total, colour = "#d53e4f" ) + 
		ylab("Phytoplankton richness\n(baseline)") + xlab("Latitude (°)") + 
		scale_y_continuous(limits = c(0,220)) + theme_bw() + coord_flip()
ggsave(plot = p5, filename = "p5.jpg", dpi = 300, width = 3, height = 7)

p6 <- ggplot() + geom_point(aes(x = y, y = phyto_fut), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = phyto_fut), data = total, colour = "#d53e4f" ) + 
		ylab("Phytoplankton richness\n(2100)") + xlab("Latitude (°)") + 
		scale_y_continuous(limits = c(0,220)) + theme_bw() + coord_flip()
ggsave(plot = p6, filename = "p6.jpg", dpi = 300, width = 3, height = 7)

p7 <- ggplot() + geom_point(aes(x = y, y = diff), alpha = 0.3, data = total, colour = "grey75") + 
		geom_smooth(aes(x = y, y = diff), data = total, colour = "#d53e4f" ) + 
		ylab("Richness difference\n(2100-2000)") + xlab("Latitude (°)") + 
		theme_bw() + coord_flip()
ggsave(plot = p7, filename = "p7.jpg", dpi = 300, width = 3, height = 7)



# --------------------------------------------------------------------------------------------------------------------------------

### 16/08/19: Examine changes in species diversity and composition together 
library("vegan")
zoo.base <- read.table("table_zoo_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
zoo.fut <- read.table("table_zoo_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
# dim(zoo.base); dim(zoo.fut)
#colnames(zoo.base)
#colnames(zoo.fut)
# summary(zoo.fut)

# Restrict to same cells before computing turn-over
zoo.base <- na.omit(zoo.base)
zoo.fut <- na.omit(zoo.fut)
summary(zoo.base$x)
# x = -179.5 - 179.5
summary(zoo.fut$x)
# x = 0-360
zoo.base$x2 <- zoo.base$x 
zoo.base[zoo.base$x < 0 ,"x2"] <- (zoo.base[zoo.base$x < 0 ,"x"]) + 360

# Change ids
zoo.base$cell_id <- factor(paste(zoo.base$x2, zoo.base$y, sep = "_"))
zoo.fut$cell_id <- factor(paste(zoo.fut$x, zoo.fut$y, sep = "_"))
# Make them follow the same order
zoo.base <- zoo.base[order(zoo.base$cell_id),]
zoo.fut <- zoo.fut[order(zoo.fut$cell_id),]
head(zoo.base$cell_id); head(zoo.fut$cell_id)

# Restrict to the same cells
# The elements of setdiff(x,y) are those elements in x but not in y.
length(setdiff(zoo.base$cell_id, zoo.fut$cell_id)) # 118 cells that are not in baseline
length(setdiff(zoo.fut$cell_id, zoo.base$cell_id)) # 223 cells not future 
zoo.base <- zoo.base[zoo.base$cell_id %in% unique(zoo.fut$cell_id),]
zoo.fut <- zoo.fut[zoo.fut$cell_id %in% unique(zoo.base$cell_id),]
length(setdiff(zoo.base$cell_id, zoo.fut$cell_id)) # 0
length(setdiff(zoo.fut$cell_id, zoo.base$cell_id)) # 0
# Ok, rbind
zoo.fut$x2 <- zoo.fut$x

### Compute Bray-Curtis dissimilarity index based on the species' average probabilities
ddf <- rbind(zoo.base, zoo.fut)
dim(ddf)
dissim <- mclapply(X = unique(ddf$cell_id), FUN = function(i) {
				message(paste(i, sep = ""))
				diss <- as.numeric(vegdist(x = ddf[ddf$cell_id == i,c(4:length(ddf))],"bray"))
				return(diss)
	}, mc.cores = 25
) # eo dissim
d <- data.frame(do.call(rbind, dissim))
d$x <- zoo.base$x2
d$y <- zoo.base$y
colnames(d)[1] <- "bray"
head(d); summary(d)
require("maps")
world2 <- map_data(map = "world2")

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = bray), data = d) +
		scale_fill_viridis(name = paste("Bray-Curtis index\n","(2100-2000)", sep = ""), limits = c(0,0.25) ) + 
		geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
		 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
		        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
		 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		 		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
		 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )				
#setwd(WD)
ggsave(plot = map, filename = paste("map_bray_zoo_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

### Just curious, does it work with that: 
#require("betapart")
#beta.div <- beta.temp(zoo.base[,c(4:527)], zoo.fut[,c(4:527)], "jaccard")
# Nope, requires 1/0

### And compute diff in richness too
require("matrixStats")
zoo.base$rich <- rowSums(as.matrix(zoo.base[,c(4:527)]))
zoo.fut$rich <- rowSums(as.matrix(zoo.fut[,c(4:527)]))
# summary(zoo.base$rich); summary(zoo.fut$rich)
zoo.fut$diff <- (zoo.fut$rich) - (zoo.base$rich)
summary(zoo.fut$diff)
# %
zoo.fut$perc <- zoo.fut$diff / (zoo.base$rich)
summary(zoo.fut$perc)

map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = diff), data = zoo.fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map4 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc), data = zoo.fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map3, filename = paste("map_richdiff_zoo_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map4, filename = paste("map_richperc_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	


### Combine with bray curtis and plot regression? 
dim(zoo.fut); dim(d)
zoo.fut$bray <- d$bray

# plot <- ggplot() + geom_point(aes(x = diff, y = bray), data = zoo.fut) +
#     ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") +
#     theme_light()
# ggsave(plot = plot, filename = paste("plot_richxbray_zoo_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 4)
#
# plot <- ggplot() + geom_point(aes(x = abs(diff), y = bray), data = zoo.fut) +
#     ylab("Change in composition\n(Bray-Curtis index)") + xlab("Absolute difference in species richness") +
#     theme_light()
# ggsave(plot = plot, filename = paste("plot_richxbray_zoo_2100-2000_rcp85_abs.jpg", sep = ""), dpi = 300, width = 7, height = 4)

### Super interesting, because you can clearly see places where relatively high changes in composition are achieved by relatively small changes in richness --> high latitudes ?

plot <- ggplot() + geom_point(aes(x = diff, y = bray, fill = abs(y)), data = zoo.fut, pch = 21, colour = "black") + 
    scale_fill_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light()    
    
ggsave(plot = plot, filename = paste("plot_richxbray_zoo_2100-2000_rcp85_lat.jpg", sep = ""), dpi = 300, width = 7, height = 4)	   
#
plot <- ggplot() + geom_point(aes(x = abs(diff), y = bray, fill = abs(y)), data = zoo.fut, pch = 21, colour = "black") + 
   scale_fill_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light()    
ggsave(plot = plot, filename = paste("plot_richxbray_zoo_2100-2000_rcp85_abs_lat.jpg", sep = ""), dpi = 300, width = 7, height = 4)	   

### Another option for plotting this, like you used to, add a domain factor (Tropics/Extratropics and High lats) and facet_grid
zoo.fut$domain <- NA
zoo.fut[which(abs(zoo.fut$y) < 30),"domain"] <- "Tropical (<30°)"
zoo.fut[which(abs(zoo.fut$y) >= 30 & abs(zoo.fut$y) < 60),"domain"] <- "Temperate (30°-60°)"
zoo.fut[which(abs(zoo.fut$y) >= 60),"domain"] <- "Polar (>60°)"
# levels(factor(zoo.fut$domain))
plot <- ggplot(data=zoo.fut) + geom_point(aes(x = diff, y = bray, fill = factor(domain)), pch = 21, colour = "black") + 
    scale_fill_manual(name = "Domain", values = c("#4575b4","#abdda4","#d73027")) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light() + facet_grid(. ~ factor(domain), scales = "fixed")
    
ggsave(plot = plot, filename = paste("plot_diffxbray_zoo_2100-2000_rcp85_facet.jpg", sep = ""), dpi = 300, width = 8, height = 3)	   

### Plot distrbution of % change and Bray-Curtis index per domain
plot2 <- ggplot(zoo.fut, aes(x=bray)) + geom_density(aes(fill = factor(domain)), alpha = 0.5) + 
        scale_fill_manual(name = "Domain", values = c("#4575b4","#abdda4","#d73027")) +
        ylab("Density") + xlab("Change in composition (Bray-Curtis index)") + 
        theme_light() + facet_grid(. ~ factor(domain), scales = "free")
#
ggsave(plot = plot2, filename = paste("plot_dens_bray_zoo_2100-2000_rcp85_facet.jpg", sep = ""), dpi = 300, width = 8, height = 3)	  


# --------------------------------------------------------------------------------------------------------------------------------

### 19/08/19/ Same as above but with phytoplankton
library("vegan")
phyto.base <- read.table("table_phyto_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
phyto.fut <- read.table("table_phyto_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
# dim(phyto.base); dim(phyto.fut)
# summary(phyto.fut)

# Restrict to same cells before computing turn-over
phyto.base <- na.omit(phyto.base)
phyto.fut <- na.omit(phyto.fut)
summary(phyto.base$x)
# x = -179.5 - 179.5
summary(phyto.fut$x)
# x = 0-360
phyto.base$x2 <- phyto.base$x 
phyto.base[phyto.base$x < 0 ,"x2"] <- (phyto.base[phyto.base$x < 0 ,"x"]) + 360

# Change ids
phyto.base$cell_id <- factor(paste(phyto.base$x2, phyto.base$y, sep = "_"))
phyto.fut$cell_id <- factor(paste(phyto.fut$x, phyto.fut$y, sep = "_"))
# Make them follow the same order
phyto.base <- phyto.base[order(phyto.base$cell_id),]
phyto.fut <- phyto.fut[order(phyto.fut$cell_id),]
head(phyto.base$cell_id); head(phyto.fut$cell_id)

# Restrict to the same cells
# The elements of setdiff(x,y) are those elements in x but not in y.
length(setdiff(phyto.base$cell_id, phyto.fut$cell_id)) # 32 cells that are not in baseline
length(setdiff(phyto.fut$cell_id, phyto.base$cell_id)) # 223 cells not future 
phyto.base <- phyto.base[phyto.base$cell_id %in% unique(phyto.fut$cell_id),]
phyto.fut <- phyto.fut[phyto.fut$cell_id %in% unique(phyto.base$cell_id),]
length(setdiff(phyto.base$cell_id, phyto.fut$cell_id)) # 0
length(setdiff(phyto.fut$cell_id, phyto.base$cell_id)) # 0
phyto.fut$x2 <- phyto.fut$x

### Compute Bray-Curtis dissimilarity index based on the species' average probabilities
ddf <- rbind(phyto.base, phyto.fut)
dim(ddf)
dissim <- mclapply(X = unique(ddf$cell_id), FUN = function(i) {
				message(paste(i, sep = ""))
				diss <- as.numeric(vegdist(x = ddf[ddf$cell_id == i,c(4:length(ddf))],"bray"))
				return(diss)
	}, mc.cores = 27
) # eo dissim
d <- data.frame(do.call(rbind, dissim))
d$x <- phyto.base$x2
d$y <- phyto.base$y
colnames(d)[1] <- "bray"
head(d); summary(d)

# mapping
require("maps")
world2 <- map_data(map = "world2")

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = bray), data = d) +
		scale_fill_viridis(name = paste("Bray-Curtis index\n","(2100-2000)", sep = ""), limits = c(0,0.25) ) + 
		geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
		 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
		        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
		 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		 		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
		 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )				
#setwd(WD)
ggsave(plot = map, filename = paste("map_bray_phyto_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

### And compute diff in richness too
phyto.base[5000:5100,c("x2","y")]
phyto.fut[5000:5100,c("x2","y")]

require("matrixStats")
phyto.base$rich <- rowSums(as.matrix(phyto.base[,c(4:342)]))
phyto.fut$rich <- rowSums(as.matrix(phyto.fut[,c(4:342)]))
# summary(zoo.base$rich); summary(zoo.fut$rich)
phyto.fut$diff <- (phyto.fut$rich) - (phyto.base$rich)
summary(phyto.fut$diff)
# %
phyto.fut$perc <- phyto.fut$diff / (phyto.base$rich)
summary(phyto.fut$perc)

map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = diff), data = phyto.fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map4 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc), data = phyto.fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map3, filename = paste("map_richdiff_phyto_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map4, filename = paste("map_richperc_phyto_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

### Combine with bray curtis and plot regression? 
dim(phyto.fut); dim(d)
phyto.fut$bray <- d$bray
 
plot <- ggplot() + geom_point(aes(x = diff, y = bray, fill = abs(y)), data = phyto.fut, pch = 21, colour = "black") + 
    scale_fill_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light()    
    
ggsave(plot = plot, filename = paste("plot_richxbray_phyto_2100-2000_rcp85_lat.jpg", sep = ""), dpi = 300, width = 7, height = 4)	   
#
plot <- ggplot() + geom_point(aes(x = abs(diff), y = bray, fill = abs(y)), data = phyto.fut, pch = 21, colour = "black") + 
   scale_fill_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light()    
ggsave(plot = plot, filename = paste("plot_richxbray_phyto_2100-2000_rcp85_abs_lat.jpg", sep = ""), dpi = 300, width = 7, height = 4)	   

### Another option for plotting this, like you used to, add a domain factor (Tropics/Extratropics and High lats) and facet_grid
phyto.fut$domain <- NA
phyto.fut[which(abs(phyto.fut$y) < 30),"domain"] <- "Tropical (<30°)"
phyto.fut[which(abs(phyto.fut$y) >= 30 & abs(phyto.fut$y) < 60),"domain"] <- "Temperate (30°-60°)"
phyto.fut[which(abs(phyto.fut$y) >= 60),"domain"] <- "Polar (>60°)"
# levels(factor(zoo.fut$domain))
plot <- ggplot(data = phyto.fut) + geom_point(aes(x = diff, y = bray, fill = factor(domain)), pch = 21, colour = "black") + 
    scale_fill_manual(name = "Domain", values = c("#4575b4","#abdda4","#d73027")) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light() + facet_grid(. ~ factor(domain), scales = "fixed")
    
ggsave(plot = plot, filename = paste("plot_diffxbray_phyto_2100-2000_rcp85_facet.jpg", sep = ""), dpi = 300, width = 8, height = 3)	   

### Plot distrbution of % change and Bray-Curtis index per domain
plot2 <- ggplot(phyto.fut, aes(x=bray)) + geom_density(aes(fill = factor(domain)), alpha = 0.5) + 
        scale_fill_manual(name = "Domain", values = c("#4575b4","#abdda4","#d73027")) +
        ylab("Density") + xlab("Change in composition (Bray-Curtis index)") + 
        theme_light() + facet_grid(. ~ factor(domain), scales = "free")
#
ggsave(plot = plot2, filename = paste("plot_dens_bray_phyto_2100-2000_rcp85_facet.jpg", sep = ""), dpi = 300, width = 8, height = 3)	  

 
### And combine both phyto & zoo
phyto.base <- read.table("table_phyto_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
phyto.fut <- read.table("table_phyto_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
zoo.base <- read.table("table_zoo_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
zoo.fut <- read.table("table_zoo_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
# dim(phyto.base); dim(zoo.base); dim(phyto.fut); dim(zoo.fut)
# Ok, they have dimensions so you can just cbind, avoid repeating the coordinates
base <- cbind(phyto.base,zoo.base[,c(4:length(zoo.base))])
fut <- cbind(phyto.fut,zoo.fut[,c(4:length(zoo.fut))])
dim(base); dim(fut)
# colnames(base)
# colnames(fut)

# Restrict to same cells before computing turn-over
base <- na.omit(base)
fut <- na.omit(fut)
summary(base$x)
# x = -179.5 - 179.5
summary(fut$x)
# x = 0-360
base$x2 <- base$x 
base[base$x < 0 ,"x2"] <- (base[base$x < 0 ,"x"]) + 360

# Change ids
base$cell_id <- factor(paste(base$x2, base$y, sep = "_"))
fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
# Make them follow the same order
base <- base[order(base$cell_id),]
fut <- fut[order(fut$cell_id),]
head(base$cell_id); head(fut$cell_id)

# Restrict to the same cells
# The elements of setdiff(x,y) are those elements in x but not in y.
length(setdiff(base$cell_id, fut$cell_id)) # 115 cells that are not in baseline
length(setdiff(fut$cell_id, base$cell_id)) # 218 cells not future 
base <- base[base$cell_id %in% unique(fut$cell_id),]
fut <- fut[fut$cell_id %in% unique(base$cell_id),]
length(setdiff(base$cell_id, fut$cell_id)) # 0
length(setdiff(fut$cell_id, base$cell_id)) # 0
fut$x2 <- fut$x

### Compute Bray-Curtis dissimilarity index based on the species' average probabilities
ddf <- rbind(base, fut)
dim(ddf)
# i <- "0.5_-10.5"
dissim <- mclapply(X = unique(ddf$cell_id), FUN = function(i) {
				message(paste(i, sep = ""))
				bray <- as.numeric(vegdist(x = ddf[ddf$cell_id == i,c(4:length(ddf))],"bray"))
                jaccard <- as.numeric(vegdist(x = ddf[ddf$cell_id == i,c(4:length(ddf))],"jaccard"))
				return(data.frame(bray = bray, jac = jaccard))
	}, mc.cores = 27
) # eo dissim
d <- data.frame(do.call(rbind, dissim))
d$x <- base$x2
d$y <- base$y
colnames(d)
head(d); summary(d)
 
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = bray), data = d) +
		scale_fill_viridis(name = paste("Bray-Curtis index\n","(2100-2000)", sep = ""), limits = c(0,0.25) ) + 
		geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
		 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
		        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
		 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		 		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
		 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )				
#setwd(WD)
ggsave(plot = map, filename = paste("map_bray_plankton_2100-2000_rcp85_HSI.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = d) +
		scale_fill_viridis(name = paste("Jaccard index\n","(2100-2000)", sep = ""), limits = c(0,0.40) ) + 
		geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
		 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
		        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
		 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		 		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
		 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )				
#setwd(WD)
ggsave(plot = map, filename = paste("map_jac_plankton_2100-2000_rcp85_HSI.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
 
 

### And compute diff in richness too
# require("matrixStats")
base$rich <- rowSums(as.matrix(base[,c(4:865)]))
fut$rich <- rowSums(as.matrix(fut[,c(4:865)]))
# summary(base$rich); summary(fut$rich)
fut$diff <- (fut$rich) - (base$rich)
summary(fut$diff)
# %
fut$perc <- fut$diff / (base$rich)
summary(fut$perc)

map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = diff), data = fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save
map4 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc), data = fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map3, filename = paste("map_richdiff_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map4, filename = paste("map_richperc_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
 
### Combine with bray curtis and plot regression? 
dim(fut); dim(d)
fut$bray <- d$bray
 
plot <- ggplot() + geom_point(aes(x = diff, y = bray, fill = abs(y)), data = fut, pch = 21, colour = "black") + 
    scale_fill_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light()    
ggsave(plot = plot, filename = paste("plot_richxbray_plankton_2100-2000_rcp85_lat.jpg", sep = ""), dpi = 300, width = 7, height = 4)	   

#
plot <- ggplot() + geom_point(aes(x = abs(diff), y = bray, fill = abs(y)), data = fut, pch = 21, colour = "black") + 
   scale_fill_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light()    
ggsave(plot = plot, filename = paste("plot_richxbray_plankton_2100-2000_rcp85_abs_lat.jpg", sep = ""), dpi = 300, width = 7, height = 4)	   

### Another option for plotting this, like you used to, add a domain factor (Tropics/Extratropics and High lats) and facet_grid
fut$domain <- NA
fut[which(abs(fut$y) < 30),"domain"] <- "Tropical (<30°)"
fut[which(abs(fut$y) >= 30 & abs(fut$y) < 60),"domain"] <- "Temperate (30°-60°)"
fut[which(abs(fut$y) >= 60),"domain"] <- "Polar (>60°)"
# levels(factor(zoo.fut$domain))
plot <- ggplot(data = fut) + geom_point(aes(x = diff, y = bray, fill = factor(domain)), pch = 21, colour = "black") + 
    scale_fill_manual(name = "Domain", values = c("#4575b4","#abdda4","#d73027")) +
    ylab("Change in composition\n(Bray-Curtis index)") + xlab("Difference in species richness") + 
    theme_light() + facet_grid(. ~ factor(domain), scales = "fixed")
    
ggsave(plot = plot, filename = paste("plot_diffxbray_plankton_2100-2000_rcp85_facet.jpg", sep = ""), dpi = 300, width = 8, height = 3)	   

### Plot distrbution of % change and Bray-Curtis index per domain
plot2 <- ggplot(fut, aes(x=bray)) + geom_density(aes(fill = factor(domain)), alpha = 0.5) + 
        scale_fill_manual(name = "Domain", values = c("#4575b4","#abdda4","#d73027")) +
        ylab("Density") + xlab("Change in composition (Bray-Curtis index)") + 
        theme_light() + facet_grid(. ~ factor(domain), scales = "free")
#
ggsave(plot = plot2, filename = paste("plot_dens_bray_plankton_2100-2000_rcp85_facet.jpg", sep = ""), dpi = 300, width = 8, height = 3)	  



### 19/08/19: Use the cooccur library (Veech et al., 2013) to assess baseline and future predicted co-occurrences
# https://cran.r-project.org/web/packages/cooccur/cooccur.pdf 
# This R package applies the probabilistic model of species co-occurrence (Veech 2013) to a set of species distributed among a set of survey or
# sampling sites. The algorithm calculates the observed and expected frequencies of co-occurrence between each pair of species.
# The expected frequency is based on the distribution of each species being random and independent of the other species.
# The analysis returns the probabilities that a more extreme (either low or high) value of co-occurrence could have been obtained by chance.
# The package also includes functions for visualizing species co-occurrence results and preparing data for downstream analyses.
# This function takes a community dataset (data frame or matrix) of species by site presence-absence data and classifies species pairs
# as having positive, negative, and random associations based on the probabilistic model of specie co-occurrence from Veech (2013).

# require("cooccur")
# rownames(base) <- base$cell_id
# cooccur.base <- cooccur(mat = base[,c(4:length(base))], type = "site_spp", thresh = T, spp_names = T, prob = "hyper")
# summary(cooccur.base)
#


# --------------------------------------------------------------------------------------------------------------------------------

### 21/08/19: Compute nb of species that are lost/ gained between the assemblages, and compute the ratio of both per cell
# First, load and combine both phyto & zoo
phyto.base <- read.table("table_phyto_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
phyto.fut <- read.table("table_phyto_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
zoo.base <- read.table("table_zoo_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
zoo.fut <- read.table("table_zoo_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
# dim(phyto.base); dim(zoo.base); dim(phyto.fut); dim(zoo.fut)
base <- cbind(phyto.base,zoo.base[,c(4:length(zoo.base))])
fut <- cbind(phyto.fut,zoo.fut[,c(4:length(zoo.fut))])
dim(base); dim(fut)
base <- na.omit(base)
fut <- na.omit(fut)
base$x2 <- base$x 
base[base$x < 0 ,"x2"] <- (base[base$x < 0 ,"x"]) + 360
# Change ids
base$cell_id <- factor(paste(base$x2, base$y, sep = "_"))
fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
# Make them follow the same order
base <- base[order(base$cell_id),]
fut <- fut[order(fut$cell_id),]
head(base$cell_id); head(fut$cell_id)
# Restrict to the same cells
base <- base[base$cell_id %in% unique(fut$cell_id),]
fut <- fut[fut$cell_id %in% unique(base$cell_id),]
fut$x2 <- fut$x
### Add a factor specifying the time period
base$period <- factor("baseline")
fut$period <- factor("future")
# Rbind
ddf <- rbind(base, fut)
# dim(ddf); colnames(ddf)

# Ok, melt and id per coords 
mddf <- melt(ddf, id.vars = colnames(ddf)[c(1,2,3,866,867)])
head(mddf)
colnames(mddf)[c(6,7)] <- c("species","HSI")
# Dcast to put periods as columns
cast <- dcast(mddf[,c(1,3:7)], cell_id + x2 + y + species ~ period, value.var = "HSI")
head(cast) # nice

# For every species and per cell_id, compute difference in mean annual HSI
library("dplyr")
diff <- data.frame(cast %>% group_by(cell_id,species) %>% summarise(x = unique(x2), y = unique(y), diff = future - baseline) ) # eo ddf
summary(diff)
# Add a factor level showing the sign of the change (positive or negative)
diff$sign <- NA
diff[diff$diff > 0,"sign"] <- "positive"
diff[diff$diff < 0,"sign"] <- "negative"
diff[diff$diff == 0,"sign"] <- "null"
# summary(factor(diff$sign))
# head(diff)

# proportions
nrow(diff[diff$sign == "positive",]) # 18090381
nrow(diff[diff$sign == "negative",]) # 12361804
# 18090381/ 12361804

# Use tally() from dplyr to count pos/neg per cell_id
# https://dplyr.tidyverse.org/reference/tally.html
diff_cell <- data.frame(diff %>% group_by(cell_id,sign) %>% tally())
#dim(diff_cell); head(diff_cell); str(diff_cell)
#diff_cell[1:100,]
diff_cell <- dcast(diff_cell, cell_id ~ factor(sign), value.var = "n")
#dim(diff_cell); head(diff_cell) # nice
# summary(diff_cell)
# And compute ratio of positive/negative per cell
ratio <- data.frame(diff_cell %>% group_by(cell_id) %>% 
        summarise(ratio = positive/negative, ratio_posi = positive/862, ratio_nega = negative/862) 
) # eo ddf
summary(ratio)
### Provide x and y from cell_id and map
#summary( as.numeric(do.call(rbind,strsplit(as.character(ratio$cell_id), "_"))[,1] ))# lon
#summary( as.numeric(do.call(rbind,strsplit(as.character(ratio$cell_id), "_"))[,2] ))# lat
ratio$x <- as.numeric(do.call(rbind,strsplit(as.character(ratio$cell_id), "_"))[,1]) 
ratio$y <- as.numeric(do.call(rbind,strsplit(as.character(ratio$cell_id), "_"))[,2])

map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = ratio_nega), data = ratio) +
 	scale_fill_viridis(name = "Ratio of losses", option = "A", limits = c(0,1)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = ratio_posi), data = ratio) +
 	scale_fill_viridis(name = "Ratio of gains", option = "A", limits = c(0,1)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map1, filename = paste("map_ratio_loss_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 4)	
ggsave(plot = map2, filename = paste("map_ratio_gain_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 4)	

# Ratio of both 
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = ratio), data = ratio) +
 	scale_fill_gradient2(name = "Gains/Losses", midpoint = 1, low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map3, filename = paste("map_ratio_gain:loss_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 4)	

### And make bivariate plots after combining with estimate of difference in richness
base$rich <- rowSums(as.matrix(base[,c(4:865)]))
fut$rich <- rowSums(as.matrix(fut[,c(4:865)]))
fut$diff <-  (fut$rich) - (base$rich)
summary(fut$diff)

map4 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = diff), data = fut) +
 	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", midpoint = 0, low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map4, filename = paste("map_richdiff_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 4)	

### Make sure ratio and fut$diff follow the same grid
length(setdiff(ratio$cell_id, fut$cell_id)) # 0
length(setdiff(fut$cell_id, ratio$cell_id)) # 0 

ratio$diff <- fut$diff

plot <- ggplot(ratio) + geom_point(aes(x = diff, y = ratio), colour = "grey70") + 
        geom_hline(yintercept=1,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") + 
        xlab("Difference in species richness") + ylab("Ratio (gains/losses)") + 
        theme_classic()
ggsave(plot = plot, filename = paste("plot_richdiffxratio_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 4)	

### Classify into the 4 type sof changes
# - dSR > 0 and ratio > 1 --> increased potential diversity with more species having an increase in HSI (invasion)
# - dSR > 0 and ratio < 1 --> increased potential diversity with more species having a decrease in HSI (more spp loose HSI but overall HSI increases)
# - dSR < 0 and ratio < 1 --> decreased potential diversity with more species having a decrease in HSI (extirpation)
# - dSR < 0 and ratio > 1 --> decreased potential diversity with more species having an increase in HSI (less spp loose HSI but overall HSI decreases)
ratio$type <- NA
ratio[which(ratio$diff > 0 & ratio$ratio > 1),"type"] <- "Full gains"
ratio[which(ratio$diff > 0 & ratio$ratio <= 1),"type"] <- "Buffered gains"
ratio[which(ratio$diff < 0 & ratio$ratio > 1),"type"] <- "Buffered losses"
ratio[which(ratio$diff < 0 & ratio$ratio <= 1),"type"] <- "Full losses"
# summary(factor(ratio$type))
# levels(factor(ratio$type))

map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(type)), data = ratio) +
	scale_fill_manual(name = "Biodiversity change", values = c("#f1b6da","#e6f5d0","#de77ae","#7fbc41") ) + 
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map5, filename = paste("map_biodiv_changes_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	


### Compute mean/median diff per species
diff_spp <- data.frame(diff %>% group_by(species) %>% summarise(med_diff = median(diff), avg_diff = mean(diff)) ) # eo ddf
summary(diff_spp)


### Another strategy: subjectively convert annual probabilities higher than 0.333 to 1 and the others to 0 and compute Jaccard
### (like Beaugrand et al., 2015)
base2 <- base
fut2 <- fut
# Convert to 1.0 with 0.3 threshold
base2[,c(4:865)][base2[,c(4:865)] > 0.3] <- 1
base2[,c(4:865)][base2[,c(4:865)] <= 0.3] <- 0
# Future 2 now
fut2[,c(4:865)][fut2[,c(4:865)] > 0.3] <- 1
fut2[,c(4:865)][fut2[,c(4:865)] <= 0.3] <- 0
# Compute richness
base2$rich <- rowSums(as.matrix(base2[,c(4:865)]))
fut2$rich <- rowSums(as.matrix(fut2[,c(4:865)]))
fut2$diff <-  (fut2$rich) - (base2$rich)
summary(fut2$diff)

map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = fut2) +
	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map6, filename = paste("map_richdiff_binom_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

fut2$perc <- fut2$diff / (base2$rich)
summary(fut2$perc)
mapX <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = fut2) +
	scale_fill_gradient2(name = "Richness difference\n(2100-2000)", low = "#3288bd", high = "#d53e4f", mid = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = mapX, filename = paste("map_richperc_binom_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	


map7 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = rich), data = base2) +
	scale_fill_viridis(name = "Species richness\n(baseline)") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map7, filename = paste("map_rich_binom_plankton_baseline.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

map8 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = rich), data = fut2) +
	scale_fill_viridis(name = "Species richness\n(future)" ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map8, filename = paste("map_rich_binom_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

### Ok, looks nice, compute jaccard and its components
### PARALLELED VERSION
m = nrow(base2[,c(4:865)])
n = ncol(base2[,c(4:865)])
### With paralelling : 
require("doParallel")
require("plyr")
registerDoParallel(cores=20)
r <- data.frame()
d <- cbind(base2[,c(4:865)], fut2[,c(4:865)]) 
d$bit <- cut(1:m, 50, labels = F)
	system.time(r <- ddply(d, ~bit, function(x) {
	 	beta.temp(x[,1:n], x[,(n+1):(2*n)],"jaccard")
},.parallel = T))
gc()
dim(r)
summary(r)

r$x <- base2$x2
r$y <- base2$y
r$diff <- fut2$diff
r$perc <- fut2$perc

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = beta.jac), data = r) +
	scale_fill_viridis(name = "Jaccard index", limits = c(0,1)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_jac_t0.3_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = beta.jtu), data = r) +
	scale_fill_viridis(name = "Turn-over", option = "A", limits = c(0,1)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_jtu_t0.3_plankton_baseline.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = beta.jne), data = r) +
	scale_fill_viridis(name = "Nestedness", option = "A", limits = c(0,1)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_jne_t0.3_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	


### ANd Jratio 
r$jratio <- (r$beta.jne)/(r$beta.jac)
summary(r)

plot <- ggplot() + geom_point(aes(x= perc, y = beta.jac, fill = beta.jtu), data = r, pch = 21, colour = "black") +
            scale_fill_viridis(name = "Turn-over", guide = "colorbar", limits = c(0,0.9) ) + 
			geom_vline(xintercept = 0, linetype = "dotted") + 
			xlab("Difference in species richness (%)") + ylab("Jaccard dissimilarity index") + 
            theme_bw()
        
ggsave(plot = plot, filename = paste("plot_diffxjacxjtu_t0.3_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 3)	


plot <- ggplot() + geom_point(aes(x= perc, y = beta.jac, fill = beta.jne), data = r, pch = 21, colour = "black") +
            scale_fill_viridis(name = "Nestedness", guide = "colorbar", limits = c(0,0.6), option = "A") + 
			geom_vline(xintercept = 0, linetype = "dotted") + 
			xlab("Difference in species richness (%)") + ylab("Jaccard dissimilarity index") + 
            theme_bw()
        
ggsave(plot = plot, filename = paste("plot_diffxjacxjne_t0.3_plankton_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 6, height = 3)	



# --------------------------------------------------------------------------------------------------------------------------------

### 22/08/19: Compute species' turn-over through time using the 'codyn' R package
require("codyn")
?turnover
# Computes species turnover between time periods as the proportion of species either gained or lost relative to the total number of
# species observed across both time periods. Includes an option to compute turnover as just the proportion of species gained (i.e.,
# "appearances") or lost (i.e., "disappearances").

# First, load and combine both phyto & zoo
phyto.base <- read.table("table_phyto_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
phyto.fut <- read.table("table_phyto_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
zoo.base <- read.table("table_zoo_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
zoo.fut <- read.table("table_zoo_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
# dim(phyto.base); dim(zoo.base); dim(phyto.fut); dim(zoo.fut)
base <- cbind(phyto.base,zoo.base[,c(4:length(zoo.base))])
fut <- cbind(phyto.fut,zoo.fut[,c(4:length(zoo.fut))])
dim(base); dim(fut)
base <- na.omit(base)
fut <- na.omit(fut)
base$x2 <- base$x 
base[base$x < 0 ,"x2"] <- (base[base$x < 0 ,"x"]) + 360
# Change ids
base$cell_id <- factor(paste(base$x2, base$y, sep = "_"))
fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
# Make them follow the same order
base <- base[order(base$cell_id),]
fut <- fut[order(fut$cell_id),]
head(base$cell_id); head(fut$cell_id)
# Restrict to the same cells
base <- base[base$cell_id %in% unique(fut$cell_id),]
fut <- fut[fut$cell_id %in% unique(base$cell_id),]
fut$x2 <- fut$x
### Add a factor specifying the time period
base$time <- factor("2000")
fut$time <- factor("2100")
# Rbind
ddf <- rbind(base, fut)
# dim(ddf); colnames(ddf)

# Remove the points in the colnames
colnames(ddf)[c(4:865)] <- gsub("[.]","",as.character(colnames(ddf)[c(4:865)])) 

# Ok, melt and id per coords 
mddf <- melt(ddf, id.vars = colnames(ddf)[c(1,2,3,866,867)])
colnames(mddf)[c(6,7)] <- c("species","HSI")
head(mddf)
str(mddf)
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# The metric was introduced by MacArthur and Wilson (1963) and modified by Diamond (1969) to be expressed as
# a proportion in order to compare turnover between sites that differ in species richness.
# turnover() requires a data frame with columns for species, time and abundance, and includes an optional argument 
# to specify a column for spatial replicates.

# Time var needs to be numeric apparently
unique( as.numeric(as.character(mddf$time)) )
mddf$time <- as.numeric(as.character(mddf$time))
# Compute turnover
# ?turnover
turnover <- turnover(df = mddf, 
                    time.var = "time", 
                    species.var = "species", 
                    abundance.var = "HSI", 
                    replicate.var = "cell_id",
                    metric = "appearance"
) # might take some time
str(turnover) # is a ddf
summary(turnover)
# dim(turnover)
head(turnover)
# provide coordinates
turnover$x <- as.numeric(do.call(rbind,strsplit(as.character(turnover$cell_id), "_"))[,1]) 
turnover$y <- as.numeric(do.call(rbind,strsplit(as.character(turnover$cell_id), "_"))[,2])
# Map

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = appearance), data = turnover) +
	scale_fill_viridis(name = "Appearance rate") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_turnover_plankton_HSI_2100-2000_rcp85.jpg", sep = ""), dpi = 300, width = 7, height = 5)	


# --------------------------------------------------------------------------------------------------------------------------------
