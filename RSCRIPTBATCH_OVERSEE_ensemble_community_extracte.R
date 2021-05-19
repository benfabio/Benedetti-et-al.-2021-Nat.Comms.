
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
} # eo firstup fun

### 1°) Set the working directories, vectors etc.
WD <- getwd()
setwd( paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future_rcp85_ensemble/", sep = "") )
zoo.wd <- getwd()
setwd( paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future_rcp85_ensemble/", sep = "") )
phyto.wd <- getwd()
setwd(WD)

# Vector of SDMs
SDMs <- c('GAM','GLM','RF','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
rcp <- "rcp85"
# Vector of earth system models
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")

### Plot distrbution of models'TSS values for zoo and phyto separately, using boxplots per SDM and facet per pool
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(zoo.wd,"/eval_score_",p,"/", sep = "") )
					zoo_spp <- str_replace_all(dir(), "eval_scores_", "")
					zoo_spp <- str_replace_all(zoo_spp, ".Rdata", "")
					scores <- lapply(zoo_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind,scores)
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
# table_scores_zoo[is.na(table_scores_zoo$TSS),] # Spinocalanus_magnus

setwd(WD)

### And for phyto
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(phyto.wd,"/eval_score_",p,"/", sep = "") )
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
# dim(table_scores_phyto); head(table_scores_phyto)
# summary(table_scores_phyto)
rm(res)
# table_scores_phyto[is.na(table_scores_phyto$TSS),] # Actiniscus_pentasterias

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
sp2rm <- scP[scP$avg_TSS <= 0.3,"species"] 
sp2rm <- c(sp2rm,"Actiniscus_pentasterias")
phyto_spp <- phyto_spp[!(phyto_spp %in% sp2rm)]


### 2°) Extract the monthly BASELINE phyto- and zooplankton communities
# # m <- "apr"
for(m in months) {

          message(paste(" ", sep = ""))
          message(paste("Retrieving baseline probabilities for ", m, sep = ""))
          message(paste(" ", sep = ""))
          # Load env variables
          setwd("/net/hydro/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
          env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
          env <- env[-which(env$SSS < 20),]
          env <- env[-which(env$Bathy > -175),]
          env$x2 <- env$x
          env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
          env$id <- paste(env$x2, env$y, sep = "_")

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
                                           if( sum(grepl("proj_projection_", dir())) == 72 )
                                           {
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
                                               #ANN
                                               resModelANN <- d[,"ANN",,]
                                               resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
                                               resModelANN <- (resModelANN/1000)
                                               # Return
                                               return( data.frame(cell_id = paste(env$x2, env$y, sep = "_"), x = env$x2, y = env$y, species = gsub("\\.","_",sp),
                                                       GAM = resModelGam, RF = resModelRF, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )
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
                      return( tbl_zoo[,c(1:4,10,11)] )

                   }, mc.cores = 10

           ) # eo mclapply
           # rbind into one table ? (then average per cell_id for average )
           table.zoo <- dplyr::bind_rows(zoo_divs_baseline_monthly)
           # Compute species' average HSI across pool and dcast to have species as columns
           d_zoo_div_base <- dcast(table.zoo, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
           # head(d_zoo_div_base); dim(d_zoo_div_base); summary(d_zoo_div_base)
           # And print monthly baseline composition as .txt file
           setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
           write.table(d_zoo_div_base, file = paste("table_zoo_mon_composition_baseline_",m,".txt", sep = ""), sep = "\t")
           rm(d_zoo_div_base, table.zoo, zoo_divs_baseline_monthly, env)
           gc()
           setwd(WD)

} # eo for loop - m in months

for(m in months) {

          message(paste(" ", sep = ""))
          message(paste("Retrieving baseline phytoplankton probabilities for ", m, sep = ""))
          message(paste(" ", sep = ""))
          # Load env variables
          setwd("/net/hydro/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
          env <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
          env <- env[-which(env$SSS < 20),]
          env <- env[-which(env$Bathy > -175),]
          env$x2 <- env$x
          env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
          env$id <- paste(env$x2, env$y, sep = "_")
   
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
                                           if( sum(grepl("proj_projection_", dir())) == 72 ) {
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
                                               return( data.frame(cell_id = paste(env$x2, env$y, sep = "_"), x = env$x2, y = env$y, species = gsub("\\.","_",sp),
                                                       GAM = resModelGam, RF = resModelRF, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )

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
                      return( tbl_phyto[,c(1:4,10,11)] )
                   }, mc.cores = 10
           ) #  eo mclapply
           # rbind into one table ? (then average per cell_id for average )
           table.phyto <- dplyr::bind_rows(phyto_divs_baseline_monthly)
           # Compute species' average HSI across pool and dcast to have species as columns
           d_phyto_div_base <- dcast(table.phyto, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
           # dim(d_phyto_div_base); head(d_phyto_div_base); summary(d_phyto_div_base)
           # And print
           setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
           write.table(d_phyto_div_base, file = paste("table_phyto_mon_composition_baseline_",m,".txt", sep = ""), sep = "\t")
           rm(d_phyto_div_base, table.phyto, phyto_divs_baseline_monthly, env)
           gc()
           setwd(WD)

} # eo for loop - m in months



### 3°) Extract the monthly FUTURE phyto- and zooplankton communities for each ESMs. Use mclapply per ESM and for loops
# esm <- "MRI-NEMURO"
# m <- "apr"
# p <- "p1"
# extract.community <- function(esm = ESMs) {
#
#              message(paste(" ", sep = ""))
#              message(paste("Extracting future Phyto- and Zooplankton monthly communities for ",esm, sep = ""))
#              message(paste(" ", sep = ""))
#
#              for(m in months) {
#
#                        message(paste(" ", sep = ""))
#                        message(paste("Retrieving probabilities for ",m, " || ", esm, sep = ""))
#                        message(paste(" ", sep = ""))
#
#                        ### 25/09/19: for bis projections, you must follow the same coordinates as in the baseline monthly clim
#                        setwd(paste("/net/hydro/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims/",esm, sep = ""))
#                        mm <- firstup(m)
#                        env <- read.table(paste("clims_mon_",mm,"_",esm,"_rcp85_base+2100-2031.txt", sep = ""), h = T, sep = "\t")
#
#                        # Load probability tables for various pools
#                        zoo_divs_fut_monthly <- mclapply(pools, function(p) {
#
#                                    message(paste(" ", sep = ""))
#                                    message(paste("Retrieving monthly future composition for pool || ",p,sep = ""))
#                                    message(paste(" ", sep = ""))
#                                    setwd(paste(zoo.wd,"/",p,"/", sep = ""))
#
#                                    require("parallel")
#                                    zoo_spp <- gsub("\\(|\\)", "", zoo_spp)
#
#                                    probas_zoo <- lapply(X = str_replace_all(zoo_spp, "_", "."), FUN = function(sp) {
#
#                                                # Go to species dir
#                                                setwd(paste(zoo.wd,"/",p,"/",sp,"/", sep = ""))
#                                                message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                                # Need to modify sp when there are 2 names and add brackets
#                                                if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
#                                                    # Then add brackets around the second piece
#                                                    sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
#                                                    strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
#                                                    strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
#                                                } # eo if loop
#
#                                                # If the 12 monthly projections are done for present & future
#                                                if( sum(grepl("proj_projection_", dir())) == 72 ) {
#
#                                                        # Load projections for each SDM
#                                                        setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100__",esm, sep = ""), sep = "") )
#                                                        d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100__",esm,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                                        # GAM
#                                                        resModelGam <- d[,"GAM",,]
#                                                        resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
#                                                        resModelGam <- (resModelGam/1000)
#                                                        # GLM
#                                                        resModelGlm <- d[,"GLM",,]
#                                                        resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
#                                                        resModelGlm <- (resModelGlm/1000)
#                                                        # RF
#                                                        resModelRF <- d[,"RF",,]
#                                                        resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
#                                                        resModelRF <- (resModelRF/1000)
#                                                        # ANN
#                                                        resModelANN <- d[,"ANN",,]
#                                                        resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
#                                                        resModelANN <- (resModelANN/1000)
#
#                                                        # Return
#                                                        return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
#                                                         GAM = resModelGam, RF = resModelRF, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )
#
#                                                } else {
#
#                                                        message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#
#                                                } # eo if else loop
#
#                                            } # eo FUN
#
#                                    ) # eo lapply
#
#                                    # rbind results from sdms
#                                    tbl_zoo <- dplyr::bind_rows(probas_zoo)
#                                    rm(probas_zoo); gc()
#                                    # Compute monthly average HSI (average across SDMs)
#                                    tbl_zoo$mean_HSI <- rowMeans( as.matrix(tbl_zoo[,c(5:6)]) )
#                                    tbl_zoo$pool <- p
#                                    return( tbl_zoo[,c(1:4,10:11)] )
#
#                                 }, mc.cores = 5
#
#                        ) # eo mclapply
#                        # rbind into one table ? (then average per cell_id for average )
#                        table.zoo <- dplyr::bind_rows(zoo_divs_fut_monthly)
#                        # Compute species' average HSI across pool and dcast to have species as columns
#                        d_zoo_div_fut <- dcast(table.zoo, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
#                        # dim(d_zoo_div_fut); head(d_zoo_div_fut); summary(d_zoo_div_fut)
#                        # And print monthly baseline composition as .txt file
#                        setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
#                        message(paste("Printing zooplankton community table for ",m, " under  ",esm, sep = ""))
#                        write.table(d_zoo_div_fut, file = paste("table_zoo_mon_composition_2100-2000","_",esm,"_",rcp,"_",m,".txt", sep = ""), sep = "\t")
#                        rm(d_zoo_div_fut, zoo_divs_fut_monthly, table.zoo, env)
#                        gc()
#                        setwd(WD)
#
#              } # eo for loop - m in months
#
#
#              ### And for future phyto monthly composition
#              for(m in months) {
#
#                       message(paste(" ", sep = ""))
#                       message(paste("Retrieving probabilities for ",m, " || ", esm, sep = ""))
#                       message(paste(" ", sep = ""))
#
#                       ### 25/09/19: for bis projections, you must follow the same coordinates as in the baseline monthly clim
#                       setwd(paste("/net/hydro/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims/",esm, sep = ""))
#                       mm <- firstup(m)
#                       env <- read.table(paste("clims_mon_",mm,"_",esm,"_rcp85_base+2100-2031.txt", sep = ""), h = T, sep = "\t")
#
#                       ### Load probability tables for various pools
#                       phyto_divs_fut_monthly <- mclapply(pools, function(p) {
#
#                                   message(paste(" ", sep = ""))
#                                   message(paste("Retrieving monthly future composition for pool || ",p,sep = ""))
#                                   message(paste(" ", sep = ""))
#
#                                   setwd(paste(phyto.wd,"/",p,"/", sep = ""))
#                                   require("parallel")
#                                   probas_phyto <- lapply(X = str_replace_all(phyto_spp, "_", "."), FUN = function(sp) {
#
#                                               # Go to species dir
#                                               setwd(paste(phyto.wd,"/",p,"/",sp,"/", sep = ""))
#                                               message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                               # Need to modify sp when there are 2 names and add brackets
#                                               if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3) {
#                                                   # Then add brackets around the second piece
#                                                   sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
#                                                   strsplit(sp, ".", fixed = TRUE)[[1]][2],").",
#                                                   strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )
#                                               } # eo if loop
#
#                                               # If the 12 monthly projections are done for present & future
#                                               if( sum(grepl("proj_projection_", dir())) == 72 ) {
#
#                                                       # Load projections for each SDM
#                                                       setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100__",esm, sep = ""), sep = "") )
#                                                       d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100__",esm,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                                       # GAM
#                                                       resModelGam <- d[,"GAM",,]
#                                                       resModelGam <- apply(resModelGam, 1, mean, na.rm = F)
#                                                       resModelGam <- (resModelGam/1000)
#                                                       # GLM
#                                                       resModelGlm <- d[,"GLM",,]
#                                                       resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F)
#                                                       resModelGlm <- (resModelGlm/1000)
#                                                       # RF
#                                                       resModelRF <- d[,"RF",,]
#                                                       resModelRF <- apply(resModelRF, 1, mean, na.rm = F)
#                                                       resModelRF <- (resModelRF/1000)
#                                                       # ANN
#                                                       resModelANN <- d[,"ANN",,]
#                                                       resModelANN <- apply(resModelANN, 1, mean, na.rm = F)
#                                                       resModelANN <- (resModelANN/1000)
#
#                                                       # Return
#                                                       return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp),
#                                                        GAM = resModelGam, RF = resModelRF, GLM = resModelGlm, RF = resModelRF, ANN = resModelANN ) )
#
#                                               } else {
#
#                                                       message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#
#                                               } # eo if else loop
#
#                                           } # eo FUN
#
#                                   ) # eo lapply
#
#                                   # rbind results from sdms
#                                   tbl_phyto <- dplyr::bind_rows(probas_phyto)
#                                   rm(probas_phyto); gc()
#                                   # Compute monthly average HSI (average across SDMs)
#                                   tbl_phyto$mean_HSI <- rowMeans( as.matrix(tbl_phyto[,c(5,6)]) )
#                                   tbl_phyto$pool <- p
#                                   return( tbl_phyto[,c(1:4,10:11)] )
#
#                                }, mc.cores = 5
#
#                       )  # eo mclapply
#
#                       # rbind into one table ? (then average per cell_id for average )
#                       table.phyto <- dplyr::bind_rows(phyto_divs_fut_monthly)
#                       # Compute species' average HSI across pool and dcast to have species as columns
#                       d_phyto_div_fut <- dcast(table.phyto, cell_id + x + y + pool ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean_HSI")
#                       # dim(d_phyto_div_fut); head(d_phyto_div_fut); summary(d_phyto_div_fut)
#                       # And print monthly baseline composition as .txt file
#                       setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
#                       message(paste("Printing phytoplankton community table for ",m, " under  ",esm, sep = ""))
#                       write.table(d_phyto_div_fut, file = paste("table_phyto_mon_composition_2100-2000","_",esm,"_",rcp,"_",m,".txt", sep = ""), sep = "\t")
#                       rm(d_phyto_div_fut, phyto_divs_fut_monthly, table.phyto, env)
#                       gc()
#                       setwd(WD)
#
#              } # eo for loop - m in months
#
# } # eo FUN - extract.community
#
# ### Apply in parallel
# require("parallel")
# mclapply(X = ESMs, FUN = extract.community, mc.cores = 5)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
