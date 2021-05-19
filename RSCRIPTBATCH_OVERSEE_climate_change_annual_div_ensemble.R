
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

# --------------------------------------------------------------------------------------------------------------------------------

firstup <- function(x) {
  substr(x,1,1) <- toupper(substr(x,1,1))
  x
}

### 1°) Set the working directories, vectors etc.
WD <- getwd()
setwd( paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future_bis_rcp85/", sep = "") )
zoo.wd <- getwd()
setwd( paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future_bis_rcp85/", sep = "") )
phyto.wd <- getwd()
setwd(WD)

# Vector of SDMs
SDMs <- c('GAM','RF')#,'GLM','RF','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5") #,"RUN6","RUN7","RUN8","RUN9","RUN10") 
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3")#,"p4")
rcp <- "rcp85"

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
# dim(table_scores_zoo); head(table_scores_zoo)
# summary(table_scores_zoo)
rm(res)
# table_scores_zoo[is.na(table_scores_zoo$TSS),]

setwd(WD)

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
#dim(table_scores_phyto); head(table_scores_phyto)
# summary(table_scores_phyto)
rm(res)
#table_scores_phyto[is.na(table_scores_phyto$TSS),]

setwd(WD)

### Identify the species to use for ensemble diversity
require("dplyr")
scZ <- data.frame(na.omit(table_scores_zoo) %>%
  		group_by(species) %>%
  		summarise(avg_TSS = mean(TSS), sd_TSS = sd(TSS))
) # eo ddf
# scZ
scZ[scZ$avg_TSS <= 0.3,] # none, only probel with ONE RUNE of S.magnus
zoo_spp <- unique(scZ$species)
zoo_spp <- zoo_spp[!(zoo_spp %in% c("Spinocalanus_magnus"))]

# For phytoplankton
scP <- data.frame(table_scores_phyto %>%
  		group_by(species) %>%
  		summarise(avg_TSS = mean(TSS), sd_TSS = sd(TSS))
) # eo ddf
# scP
phyto_spp <- unique(scP$species)
sp2rm <- scP[scP$avg_TSS <= 0.3,"species"] # generally taxa with presences all over the place
sp2rm <- c(sp2rm,"Actiniscus_pentasterias")
phyto_spp <- phyto_spp[!(phyto_spp %in% sp2rm)]


### 03/04/20: Get the modelled species class and rbind phyto and zooplato in one table (for Appendix S8)
# setwd(paste(WD,"/","biology", sep=""))
# #phyto_spp
# #zoo_spp
# phyto <- data.frame(species = phyto_spp); zoo <- data.frame(species = zoo_spp)
# classifP <- get(load("classif_phytoplankton.Rdata"))
# classifZ <- get(load("classif_zooplankton.Rdata"))
# colnames(classifP) ; colnames(classifZ)
# unique(classifP$species) ; unique(phyto) # OK phyto is clear
# unique(classifZ$species) ; unique(zoo) # OK zoo is clear
# # Add classif columns for phyto
# phyto$genus <- NA
# phyto$family <- NA
# phyto$order <- NA
# phyto$class <- NA
# phyto$phylum <- NA
# phyto$group <- NA
# # Add classif columns for zoo
# zoo$genus <- NA
# zoo$family <- NA
# zoo$order <- NA
# zoo$class <- NA
# zoo$phylum <- NA
# zoo$group <- NA
# #sp <- unique(phyto$species)[3]
# for(sp in unique(phyto$species)) {
#     message(paste(sp, sep = ""))
#     phyto[phyto$species == sp,"genus"] <- as.character(classifP[classifP$species == sp,"genus"])
#     phyto[phyto$species == sp,"family"] <- as.character(classifP[classifP$species == sp,"family"])
#     phyto[phyto$species == sp,"order"] <- as.character(classifP[classifP$species == sp,"order"])
#     phyto[phyto$species == sp,"class"] <- as.character(classifP[classifP$species == sp,"class"])
#     phyto[phyto$species == sp,"phylum"] <- as.character(classifP[classifP$species == sp,"phylum"])
#     phyto[phyto$species == sp,"group"] <- as.character(classifP[classifP$species == sp,"group"])
# }
#
# for(sp in unique(zoo$species)) {
#     message(paste(sp, sep = ""))
#     zoo[zoo$species == sp,"genus"] <- as.character(classifZ[classifZ$species == sp,"genus"])
#     zoo[zoo$species == sp,"family"] <- as.character(classifZ[classifZ$species == sp,"family"])
#     zoo[zoo$species == sp,"order"] <- as.character(classifZ[classifZ$species == sp,"order"])
#     zoo[zoo$species == sp,"class"] <- as.character(classifZ[classifZ$species == sp,"class"])
#     zoo[zoo$species == sp,"phylum"] <- as.character(classifZ[classifZ$species == sp,"phylum"])
#     zoo[zoo$species == sp,"group"] <- as.character(classifZ[classifZ$species == sp,"group"])
# }
#
# ### Add trophic level
# phyto$level <- "Phytoplankton"
# zoo$level <- "Zooplankton"
#
# # Rbind
# table <- rbind(phyto, zoo)
# # Replace "_" in the species names
# table$species <- str_replace_all(as.character(table$species),"_"," ")
# write.table(table, file = "table_appendixS8.txt", sep = ";")

### 3°) For each kingdom and for each time periods, extratc species monthly HSI, average to obtain mean annual composition and save
# zoo_divs_baseline <- lapply(pools, function(p) {
#
#                  message(paste(" ", sep = ""))
#                  message(paste("Computing annual Shannon diversity index for pool || ", p, sep = ""))
#                  message(paste(" ", sep = ""))
#                  # Extract monthly probabilities
#                  res <- lapply(months, function(m)
#                  {
#                                  message(paste(" ", sep = ""))
#                                  message(paste("Retrieving probabilities for ", m, sep = ""))
#                                  message(paste(" ", sep = ""))
#                                  # Load env variables
#                                  setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
#                                  env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
#                                  env <- env[-which(env$SSS < 20),]
#                                  env <- env[-which(env$Bathy > -175),]
#                                  # Load phyto probas
#                                  message(paste("Loading zoo projections ================================  ", sep = ""))
#                                  setwd(paste(zoo.wd,"/",p,"/", sep = ""))
#                                  require("parallel")
#                                  zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
#                                  probas_zoo <- lapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {
#
#                                      # Got to species dir
#                                      setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
#                                      message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                      # Need to modify sp when there are 2 names and add brackets
#                                      if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
#                                          # Then add brackets around the second piece
#                                          sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
#                                      } # eo if loop
#
#                                      # If the 4 seasonal projections are done
#                                      if( sum(grepl("proj_projection_", dir())) == 24 )
#                                     {
#
#                                          ### Load projections for each SDM
#                                          setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
#                                          d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                          # GAM
#                                          resModelGam <- d[,"GAM",,]
#                                          resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
#                                          resModelGam <- (resModelGam/1000)
#                                          # GLM
#                                          resModelGlm <- d[,"GLM",,]
#                                          resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
#                                          resModelGlm <- (resModelGlm/1000)
#                                          # RF
#                                          resModelRF <- d[,"RF",,]
#                                          resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
#                                          resModelRF <- (resModelRF/1000)
#                                          # ANN
#                                          resModelANN <- d[,"ANN",,]
#                                          resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
#                                          resModelANN <- (resModelANN/1000)
#
#                                          # Return
#                                          return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
#                                                  GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
#
#                                      } else {
#
#                                          message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#
#                                      } # eo if else loop
#
#                                  } ) # eo lapply
#                                  # cbind SDMs' mean HSI
#                                  tbl_zoo <- dplyr::bind_rows(probas_zoo)
#                                  rm(probas_zoo); gc()
#                                  # Compute average HSI (average across SDMs)
#                                  tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:length(tbl_zoo))]) )
#                                  return( tbl_zoo[,c(1:4,9)] )
#                              }
#
#                  ) # eo lapply
#                  # rbind into one table ? (then average per cell_id for average )
#                  table.zoo <- dplyr::bind_rows(res)
#                  rm(res)
#                  # Compute annual HSI for each species
#                 comm <- dcast(table.zoo, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#                 # dim(comm); summary(comm)
#                  rm(table.zoo)
#                  # Return zoo.div
#                  comm$pool <- p
#                  return(comm)
#
#          }  # eo fun
#
# ) # eo 1st lapply - pools
# # Rbind
# zoo_div_base <- dplyr::bind_rows(zoo_divs_baseline)
#
# # Melt, average, and dcast?
# m_zoo_div_base <- melt(zoo_div_base, id = c("cell_id","x","y","pool"))
# colnames(m_zoo_div_base)[c(5:6)] <- c("species","HSI")
#
# # Re-dcast
# d_zoo_div_base <- dcast(m_zoo_div_base, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#
# # Save
# setwd(WD)
# write.table(d_zoo_div_base, file = "table_zoo_annual_composition_baseline.txt", sep = "\t")
# rm(m_zoo_div_base,zoo_div_base,zoo_divs_baseline)
# gc()
#
#
# # -------------------------------------------------------
#
# rcp <- "rcp85"
#
# ### Next, same but for future zoopl composition
# zoo_divs_fut <- lapply(pools, function(p) {
#
#                  message(paste(" ", sep = ""))
#                  message(paste("Computing annual Shannon diversity index for pool || ", p, sep = ""))
#                  message(paste(" ", sep = ""))
#                  # Extract monthly probabilities
#                  res <- lapply(months, function(m)
#                  {
#                                  message(paste(" ", sep = ""))
#                                  message(paste("Retrieving probabilities for ", m, sep = ""))
#                                  message(paste(" ", sep = ""))
#                                 if(m == "jan") {
#                                     mm <- "Jan"
#                                 } else if(m == "feb") {
#                                     mm <- "Feb"
#                                 } else if(m == "mar") {
#                                     mm <- "Mar"
#                                 } else if(m == "apr") {
#                                     mm <- "Apr"
#                                 } else if(m == "may") {
#                                     mm <- "May"
#                                 } else if(m == "jun") {
#                                     mm <- "Jun"
#                                 } else if(m == "jul") {
#                                     mm <- "Jul"
#                                 } else if(m == "aug") {
#                                     mm <- "Aug"
#                                 } else if(m == "sep") {
#                                     mm <- "Sep"
#                                 } else if(m == "oct") {
#                                     mm <- "Oct"
#                                 } else if(m == "nov") {
#                                     mm <- "Nov"
#                                 } else if(m == "dec") {
#                                     mm <- "Dec"
#                                 } #
#
#                                  # Load env variables
#                                  setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
#                                  env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
#
#                                  # Load zoo probas
#                                  message(paste("Loading zoo projections ================================  ", sep = ""))
#                                  setwd(paste(zoo.wd,"/",p,"/", sep = ""))
#                                  require("parallel")
#                                  zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
#                                  probas_zoo <- lapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {
#
#                                      # Got to species dir
#                                      setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
#                                      message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                      # Need to modify sp when there are 2 names and add brackets
#                                      if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
#                                          # Then add brackets around the second piece
#                                          sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
#                                      } # eo if loop
#
#                                      # If the 12 monthly projections are done for present & future
#                                      if( sum(grepl("proj_projection_", dir())) == 24 ) {
#
#                                          # Load projections for each SDM
#                                          setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M", sep = ""), sep = "") )
#                                          d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                          # GAM
#                                          resModelGam <- d[,"GAM",,]
#                                          resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
#                                          resModelGam <- (resModelGam/1000)
#                                          # GLM
#                                          resModelGlm <- d[,"GLM",,]
#                                          resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
#                                          resModelGlm <- (resModelGlm/1000)
#                                          # RF
#                                          resModelRF <- d[,"RF",,]
#                                          resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
#                                          resModelRF <- (resModelRF/1000)
#                                          # ANN
#                                          resModelANN <- d[,"ANN",,]
#                                          resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
#                                          resModelANN <- (resModelANN/1000)
#
#                                          # Return
#                                          return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
#                                                  GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
#
#                                      } else {
#
#                                          message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#
#                                      } # eo if else loop
#
#                                  } ) # eo lapply
#                                  # cbind SDMs' mean HSI
#                                  tbl_zoo <- dplyr::bind_rows(probas_zoo)
#                                  rm(probas_zoo); gc()
#                                  # Compute average HSI (average across SDMs)
#                                  tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:length(tbl_zoo))]) )
#                                  return( tbl_zoo[,c(1:4,9)] )
#                              }
#
#                  ) # eo lapply
#                  # rbind into one table ? (then average per cell_id for average )
#                  table.zoo <- dplyr::bind_rows(res)
#                  rm(res)
#                  # Compute annual HSI for each species
#                 comm <- dcast(table.zoo, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#                 # dim(comm); summary(comm)
#                  rm(table.zoo)
#                  # Return zoo.div
#                  comm$pool <- p
#                  return(comm)
#
#          }  # eo fun
#
# ) # eo 1st lapply - pools
# # Rbind
# zoo_div_fut <- dplyr::bind_rows(zoo_divs_fut)
#
# # Melt, average, and dcast?
# m_zoo_div_fut <- melt(zoo_div_fut, id = c("cell_id","x","y","pool"))
# colnames(m_zoo_div_fut)[c(5:6)] <- c("species","HSI")
#
# # Re-dcast
# d_zoo_div_fut <- dcast(m_zoo_div_fut, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#
# # Save
# setwd(WD)
# write.table(d_zoo_div_fut, file = "table_zoo_annual_composition_2100-2071_rcp85.txt", sep = "\t")
# rm(m_zoo_div_fut,zoo_div_fut,zoo_divs_fut)
# gc()
#
#
#
# # -------------------------------------------------------
#
# ### And now, phytoplankton baseline composition
# phyto_divs_baseline <- lapply(pools, function(p) {
#
#                  message(paste(" ", sep = ""))
#                  message(paste("Computing annual Shannon diversity index for pool || ", p, sep = ""))
#                  message(paste(" ", sep = ""))
#                  # Extract monthly probabilities
#                  res <- lapply(months, function(m)
#                  {
#                                  message(paste(" ", sep = ""))
#                                  message(paste("Retrieving probabilities for ", m, sep = ""))
#                                  message(paste(" ", sep = ""))
#                                  # Load env variables
#                                  setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
#                                  env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
#                                  env <- env[-which(env$SSS < 20),]
#                                  env <- env[-which(env$Bathy > -175),]
#
#                                  # Load phyto probas
#                                  message(paste("Loading zoo projections ================================  ", sep = ""))
#                                  setwd(paste(phyto.wd,"/",p,"/", sep = ""))
#                                  require("parallel")
#                                  probas <- lapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {
#
#                                      # Got to species dir
#                                      setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
#                                      message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                      # Need to modify sp when there are 2 names and add brackets
#                                      if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
#                                          # Then add brackets around the second piece
#                                          sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
#                                      } # eo if loop
#
#                                      # If the 4 seasonal projections are done
#                                      if( sum(grepl("proj_projection_", dir())) == 24 ) {
#
#                                          ### Load projections for each SDM
#                                          setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
#                                          d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                          # GAM
#                                          resModelGam <- d[,"GAM",,]
#                                          resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
#                                          resModelGam <- (resModelGam/1000)
#                                          # GLM
#                                          resModelGlm <- d[,"GLM",,]
#                                          resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
#                                          resModelGlm <- (resModelGlm/1000)
#                                          # RF
#                                          resModelRF <- d[,"RF",,]
#                                          resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
#                                          resModelRF <- (resModelRF/1000)
#                                          # ANN
#                                          resModelANN <- d[,"ANN",,]
#                                          resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
#                                          resModelANN <- (resModelANN/1000)
#
#                                          # Return
#                                          return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
#                                                  GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
#
#                                      } else {
#
#                                          message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#
#                                      } # eo if else loop
#
#                                  } ) # eo lapply
#                                  # cbind SDMs' mean HSI
#                                  tbl_phyto <- dplyr::bind_rows(probas)
#                                  rm(probas); gc()
#                                  # Compute average HSI (average across SDMs)
#                                  tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5:length(tbl_phyto))]) )
#                                  return( tbl_phyto[,c(1:4,9)] )
#                              }
#                  ) # eo lapply
#                  # rbind into one table ? (then average per cell_id for average )
#                  table.phyto <- dplyr::bind_rows(res)
#                  rm(res)
#                  # Compute annual HSI for each species
#                 comm <- dcast(table.phyto, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#                  rm(table.phyto)
#                  # Return zoo.div
#                  comm$pool <- p
#                  return(comm)
#
#          }  # eo fun
#
# ) # eo 1st lapply - pools
# # Rbind
# phyto_div_base <- dplyr::bind_rows(phyto_divs_baseline)
#
# # Melt, average, and dcast?
# m_phyto_div_base <- melt(phyto_div_base, id = c("cell_id","x","y","pool"))
# colnames(m_phyto_div_base)[c(5:6)] <- c("species","HSI")
# # Re-dcast
# d_phyto_div_base <- dcast(m_phyto_div_base, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#
# # Save
# setwd(WD)
# write.table(d_phyto_div_base, file = "table_phyto_annual_composition_baseline.txt", sep = "\t")
# rm(phyto_divs_baseline,m_phyto_div_base,phyto_div_base); gc()
#
#
# # -------------------------------------------------------
#
# phyto_divs_fut <- lapply(pools, function(p) {
#
#                  message(paste(" ", sep = ""))
#                  message(paste("Computing annual Shannon diversity index for pool || ", p, sep = ""))
#                  message(paste(" ", sep = ""))
#                  # Extract monthly probabilities
#                  res <- lapply(months, function(m)
#                  {
#                                  message(paste(" ", sep = ""))
#                                  message(paste("Retrieving probabilities for ", m, sep = ""))
#                                  message(paste(" ", sep = ""))
#                                 if(m == "jan") {
#                                     mm <- "Jan"
#                                 } else if(m == "feb") {
#                                     mm <- "Feb"
#                                 } else if(m == "mar") {
#                                     mm <- "Mar"
#                                 } else if(m == "apr") {
#                                     mm <- "Apr"
#                                 } else if(m == "may") {
#                                     mm <- "May"
#                                 } else if(m == "jun") {
#                                     mm <- "Jun"
#                                 } else if(m == "jul") {
#                                     mm <- "Jul"
#                                 } else if(m == "aug") {
#                                     mm <- "Aug"
#                                 } else if(m == "sep") {
#                                     mm <- "Sep"
#                                 } else if(m == "oct") {
#                                     mm <- "Oct"
#                                 } else if(m == "nov") {
#                                     mm <- "Nov"
#                                 } else if(m == "dec") {
#                                     mm <- "Dec"
#                                 } #
#
#                                  # Load env variables
#                                  setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
#                                  env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
#
#                                  # Load phyto probas
#                                  message(paste("Loading zoo projections ================================  ", sep = ""))
#                                  setwd(paste(phyto.wd,"/",p,"/", sep = ""))
#                                  require("parallel")
#                                  probas <- lapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {
#
#                                      # Got to species dir
#                                      setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
#                                      message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                      # Need to modify sp when there are 2 names and add brackets
#                                      if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
#                                          # Then add brackets around the second piece
#                                          sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
#                                                  strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
#                                      } # eo if loop
#
#                                      # If the 12 monthly projections are done for present & future
#                                      if( sum(grepl("proj_projection_", dir())) == 24 ) {
#
#                                          # Load projections for each SDM
#                                          setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M", sep = ""), sep = "") )
#                                          d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_GFDLESM2M_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                          # GAM
#                                          resModelGam <- d[,"GAM",,]
#                                          resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
#                                          resModelGam <- (resModelGam/1000)
#                                          # GLM
#                                          resModelGlm <- d[,"GLM",,]
#                                          resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
#                                          resModelGlm <- (resModelGlm/1000)
#                                          # RF
#                                          resModelRF <- d[,"RF",,]
#                                          resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
#                                          resModelRF <- (resModelRF/1000)
#                                          # ANN
#                                          resModelANN <- d[,"ANN",,]
#                                          resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
#                                          resModelANN <- (resModelANN/1000)
#
#                                          # Return
#                                          return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
#                                                  GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
#
#                                      } else {
#
#                                          message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#
#                                      } # eo if else loop
#
#                                  } ) # eo lapply
#                                  # cbind SDMs' mean HSI
#                                  tbl_phyto <- dplyr::bind_rows(probas)
#                                  rm(probas); gc()
#                                  # Compute average HSI (average across SDMs)
#                                  tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5:length(tbl_phyto))]) )
#                                  return( tbl_phyto[,c(1:4,9)] )
#                              }
#
#                  ) # eo lapply
#                  # rbind into one table ? (then average per cell_id for average )
#                  table.phyto <- dplyr::bind_rows(res)
#                  rm(res)
#                  # Compute annual HSI for each species
#                 comm <- dcast(table.phyto, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#                  rm(table.phyto)
#                  # Return zoo.div
#                  comm$pool <- p
#                  return(comm)
#
#          }  # eo fun
#
# ) # eo 1st lapply - pools
# # Rbind
# phyto_div_fut <- dplyr::bind_rows(phyto_divs_fut)
#
# # Melt, average, and dcast?
# m_phyto_div_fut <- melt(phyto_div_fut, id = c("cell_id","x","y","pool"))
# colnames(m_phyto_div_fut)[c(5:6)] <- c("species","HSI")
#
# # Re-dcast
# d_phyto_div_fut <- dcast(m_phyto_div_fut, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T)
#
# # Save
# setwd(WD)
# write.table(d_phyto_div_fut, file = "table_phyto_annual_composition_2100-2071_rcp85.txt", sep = "\t")
# rm(m_phyto_div_fut,phyto_div_fut,phyto_divs_fut)
# gc()


# -------------------------------------------------------


### 05/09/19: Same as above, but extract MONTHLY baseline and future community compositons 
# for testing
# m <- "apr"
# p <- "p1"

### 25/09/19: Extract baseline and future species composition (phyto & zoo) for bis projections ! (keep pools information)
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
rcp <- "rcp85"

for(m in months) {

         message(paste(" ", sep = ""))
         message(paste("Retrieving probabilities for ", m, sep = ""))
         message(paste(" ", sep = ""))
         # Load env variables
         setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
         env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
         env <- env[-which(env$SSS < 20),]
         env <- env[-which(env$Bathy > -175),]
         env$x2 <- env$x 
         env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
         env$id <- paste(env$x2, env$y, sep = "_")
         env <- env[order(env$id),]
         
         zoo_divs_baseline_monthly <- mclapply(pools, function(p) {

                      message(paste(" ", sep = ""))
                      message(paste("Retrieving monthly community composition for pool || ", p, sep = ""))
                      message(paste(" ", sep = ""))

                      # Load species probas
                      message(paste("Loading zoo projections ================================  ", sep = ""))
                      setwd(paste(zoo.wd,"/",p,"/", sep = ""))
                      require("parallel")
                      zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
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

                                          # If the 4 seasonal projections are done
                                          if( sum(grepl("proj_projection_", dir())) == 24 )
                                          {
                                              # Load projections for each SDM
                                              setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
                                              d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
                                              # GAM
                                              resModelGam <- d[,"GAM",,]
                                              resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
                                              resModelGam <- (resModelGam/1000)
                                              # GLM
                                              #resModelGlm <- d[,"GLM",,]
                                              #resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
                                              #resModelGlm <- (resModelGlm/1000)
                                              # RF
                                              resModelRF <- d[,"RF",,]
                                              resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
                                              resModelRF <- (resModelRF/1000)
                                              #ANN
                                              #resModelANN <- d[,"ANN",,]
                                              #resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
                                              #resModelANN <- (resModelANN/1000)
                                              # Return
                                              return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
                                                      GAM = resModelGam, RF = resModelRF) ) #, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )
                                          } else {

                                              message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

                                          } # eo if else loop

                            } # eo FUN

                     ) # eo lapply
                     # cbind SDMs' mean HSI
                     tbl_zoo <- dplyr::bind_rows(probas_zoo)
                     rm(probas_zoo); gc()
                     # Compute monthly average HSI (average across SDMs)
                     tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:6)]) )
                     tbl_zoo$pool <- p
                     # Return
                     return( tbl_zoo[,c(1:4,7:8)] )

                  }, mc.cores = 5

          ) # eo mclapply
          # rbind into one table ? (then average per cell_id for average )
          table.zoo <- dplyr::bind_rows(zoo_divs_baseline_monthly)
          # Compute species' average HSI across pool and dcast to have species as columns
          d_zoo_div_base <- dcast(table.zoo, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
          # head(d_zoo_div_base); dim(d_zoo_div_base); summary(d_zoo_div_base)
          # And print monthly baseline composition as .txt file
          setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition")
          write.table(d_zoo_div_base, file = paste("table_zoo_monthly_composition_baseline_",m,"_bis.txt", sep = ""), sep = "\t")
          rm(d_zoo_div_base, table.zoo, zoo_divs_baseline_monthly, env)
          gc()
          setwd(WD)
     
} # eo for loop - m in months


# Same as above, but for the future conditions
rcp <- "rcp85"

for(m in months) {

         message(paste(" ", sep = ""))
         message(paste("Retrieving probabilities for ",m, sep = ""))
         message(paste(" ", sep = ""))
         # Load env variables
         #mm <- firstup(m)
         #setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
         #env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
         
         ### 25/09/19: for bis projections, you must follow the same coordinates as in the baseline monthly clim
         setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
         env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
         env <- env[-which(env$SSS < 20),]
         env <- env[-which(env$Bathy > -175),]
         #
         env$x2 <- env$x 
         env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
         env$id <- paste(env$x2, env$y, sep = "_")
         # 
         env <- env[order(env$id),]
       
         
         # Load probability tables for various pools
         zoo_divs_fut_monthly <- mclapply(pools, function(p) {
                     message(paste(" ", sep = ""))
                     message(paste("Retrieving monthly future composition for pool || ",p,sep = ""))
                     message(paste(" ", sep = ""))
                     setwd(paste(zoo.wd,"/",p,"/", sep = ""))
                     require("parallel")
                     zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
                     probas_zoo <- lapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {
                                 # Go to species dir
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
                                         #resModelGlm <- d[,"GLM",,]
                                         #resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
                                         #resModelGlm <- (resModelGlm/1000)
                                         # RF
                                          resModelRF <- d[,"RF",,]
                                          resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
                                          resModelRF <- (resModelRF/1000)
                                         # ANN
                                         #resModelANN <- d[,"ANN",,]
                                         #resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
                                         #resModelANN <- (resModelANN/1000)

                                         # Return
                                         return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
                                          GAM = resModelGam, RF = resModelRF) ) #, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )

                                 } else {

                                         message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

                                 } # eo if else loop

                             } # eo FUN

                     ) # eo lapply
                     # rbind results from sdms
                     tbl_zoo <- dplyr::bind_rows(probas_zoo)
                     rm(probas_zoo); gc()
                     # Compute monthly average HSI (average across SDMs)
                     tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:6)]) )
                     tbl_zoo$pool <- p
                     return( tbl_zoo[,c(1:4,7,8)] )
                     
                  }, mc.cores = 5
         ) # eo mclapply
         # rbind into one table ? (then average per cell_id for average )
         table.zoo <- dplyr::bind_rows(zoo_divs_fut_monthly)
         # Compute species' average HSI across pool and dcast to have species as columns
         d_zoo_div_fut <- dcast(table.zoo, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
         # dim(d_zoo_div_fut); head(d_zoo_div_fut); summary(d_zoo_div_fut)
         # And print monthly baseline composition as .txt file
         setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition")
         write.table(d_zoo_div_fut, file = paste("table_zoo_monthly_composition_2100-2000","_GFDL-ESM2M_",rcp,"_",m,"_bis.txt", sep = ""), sep = "\t")
         rm(d_zoo_div_fut, zoo_divs_fut_monthly, table.zoo, env)
         gc()
         setwd(WD)
     
} # eo for loop - m in months


#  -------------------------------------------------------
  
# m <- "apr"
#  And for phytoplankton
for(m in months) {

         message(paste(" ", sep = ""))
         message(paste("Retrieving probabilities for ", m, sep = ""))
         message(paste(" ", sep = ""))
         # Load env variables
         setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
         env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
         env <- env[-which(env$SSS < 20),]
         env <- env[-which(env$Bathy > -175),]
         env$x2 <- env$x 
         env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
         env$id <- paste(env$x2, env$y, sep = "_")
         env <- env[order(env$id),]
         
         phyto_divs_baseline_monthly <- mclapply(pools, function(p) {
             
                      message(paste(" ", sep = ""))
                      message(paste("Retrieving monthly community composition for pool || ", p, sep = ""))
                      message(paste(" ", sep = ""))
                      # Load species probas
                      setwd(paste(phyto.wd,"/",p,"/", sep = ""))
                      require("parallel")
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

                                         # If the 4 seasonal projections are done
                                          if( sum(grepl("proj_projection_", dir())) == 24 ) {
                                              # Load projections for each SDM
                                              setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
                                              d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
                                              # GAM
                                              resModelGam <- d[,"GAM",,]
                                              resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
                                              resModelGam <- (resModelGam/1000)
                                              # GLM
                                              # resModelGlm <- d[,"GLM",,]
                                              #resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
                                              #resModelGlm <- (resModelGlm/1000)
                                              # RF
                                              resModelRF <- d[,"RF",,]
                                              resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
                                              resModelRF <- (resModelRF/1000)
                                              # ANN
                                              #resModelANN <- d[,"ANN",,]
                                              #resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
                                              #resModelANN <- (resModelANN/1000)
                                              # Return
                                              return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
                                                      GAM = resModelGam, RF = resModelRF) ) #, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )

                                          } else {

                                              message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = "") )

                                          } # eo if else loop

                            } # eo FUN

                     ) # eo lapply
                     # cbind SDMs' mean HSI
                     tbl_phyto <- dplyr::bind_rows(probas_phyto)
                     rm(probas_phyto); gc()
                     # Compute monthly average HSI (average across SDMs)
                     tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5,6)]) )
                     tbl_phyto$pool <- p
                     # Return
                     return( tbl_phyto[,c(1:4,7,8)] )
                  }, mc.cores = 5

          ) #  eo mclapply
          # rbind into one table ? (then average per cell_id for average )
          table.phyto <- dplyr::bind_rows(phyto_divs_baseline_monthly)
          # Compute species' average HSI across pool and dcast to have species as columns
          d_phyto_div_base <- dcast(table.phyto, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
          # dim(d_phyto_div_base); head(d_phyto_div_base); summary(d_phyto_div_base)
          # And print
          setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition")
          write.table(d_phyto_div_base, file = paste("table_phyto_monthly_composition_baseline_",m,"_bis.txt", sep = ""), sep = "\t")
          rm(d_phyto_div_base, table.phyto, phyto_divs_baseline_monthly, env)
          gc()
          setwd(WD)

 } # eo for loop - m in months


### And for future phyto monthly composition

for(m in months) {
    
        message(paste(" ", sep = ""))
        message(paste("Retrieving probabilities for ",m, sep = ""))
        message(paste(" ", sep = ""))
        # Load env variables
        #mm <- firstup(m) 
        #setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
        #env <- read.table(paste("clim_2100-2071_rcp85_diff_",mm,"_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
        
        ### 25/09/19: for bis projections, you must follow the same coordinates as in the baseline monthly clim
        setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
        env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
        env <- env[-which(env$SSS < 20),]
        env <- env[-which(env$Bathy > -175),]
        env$x2 <- env$x 
        env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
        env$id <- paste(env$x2, env$y, sep = "_")
        env <- env[order(env$id),]
        
        ### Load probability tables for various pools
        phyto_divs_fut_monthly <- mclapply(pools, function(p) {
     				message(paste(" ", sep = ""))
     				message(paste("Retrieving monthly future composition for pool || ",p,sep = ""))
     				message(paste(" ", sep = ""))
                    setwd(paste(phyto.wd,"/",p,"/", sep = ""))
                    require("parallel")
                    probas_phyto <- lapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {
                                
                                # Go to species dir
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
                                        #resModelGlm <- d[,"GLM",,]
                                        #resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
                                        #resModelGlm <- (resModelGlm/1000)
                                        # RF
                                        resModelRF <- d[,"RF",,]
                                        resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
                                        resModelRF <- (resModelRF/1000)
                                        # ANN
                                        #resModelANN <- d[,"ANN",,]
                                        #resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
                                        #resModelANN <- (resModelANN/1000)
                    
                                        # Return
                                        return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
                                         GAM = resModelGam, RF = resModelRF)) #, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )
                    
                                } else {
                    
                                        message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
                    
                                } # eo if else loop
                    
                            } # eo FUN
                             
                    ) # eo lapply
                    # rbind results from sdms
     				tbl_phyto <- dplyr::bind_rows(probas_phyto)
     				rm(probas_phyto); gc()
     				# Compute monthly average HSI (average across SDMs)
     				tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5,6)]) )
                    tbl_phyto$pool <- p
     				return( tbl_phyto[,c(1:4,7,8)] )
     			}, mc.cores = 5
     	)  # eo mclapply
     	# rbind into one table ? (then average per cell_id for average )
     	table.phyto <- dplyr::bind_rows(phyto_divs_fut_monthly)
        # Compute species' average HSI across pool and dcast to have species as columns
        d_phyto_div_fut <- dcast(table.phyto, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
        # dim(d_phyto_div_fut); head(d_phyto_div_fut); summary(d_phyto_div_fut)
        # And print monthly baseline composition as .txt file
        setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition")
        write.table(d_phyto_div_fut, file = paste("table_phyto_monthly_composition_2100-2000","_GFDL-ESM2M_",rcp,"_",m,"_bis.txt", sep = ""), sep = "\t")
        rm(d_phyto_div_fut, phyto_divs_fut_monthly, table.phyto, env)
        gc()
        setwd(WD)
    
} # eo for loop - m in months

### 02/09/19: Test maps to check coords
# ddd <- dcast(tbl_phyto, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
# dim(ddd); summary(ddd)
# map <- ggplot() + geom_raster(aes(x = x, y = y, fill = Amphisolenia_bidentata), data = ddd) +
#           scale_fill_viridis(name = "HSI") + coord_quickmap() + theme_bw()
# setwd("/net/kryo/work/fabioben/OVERSEE/data")
# ggsave(plot = map, filename = paste("test_map3.pdf", sep = ""), dpi = 300, width = 7, height = 5)




# ---------------------------------------------------------------------------------------------------------------------------------------------------


