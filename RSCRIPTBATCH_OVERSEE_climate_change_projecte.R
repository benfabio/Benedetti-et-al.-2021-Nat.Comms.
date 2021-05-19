
### ================================================================================================================

library("raster")
library("sp")
library("rJava") 
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")

### ================================================================================================================


### Set the dataset and the background selection strategy
set <- "v9v3.1"
strat <- "total"

### Set the working directories
# For zooplankton species data : 
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.1/total_background/")
WD <- getwd()
# Niche modelling dir

# Vector of SDMs
SDMs <- c('GLM','GAM','RF','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 

# Vector of spp names
files <- dir()[grep(".txt",dir())]
species <- stringr::str_replace_all(string = files, pattern = paste("data_",strat,"_", sep =""), replacement = "")
species <- stringr::str_replace_all(string = species, pattern = ".txt", replacement = "")
# species

### Load the env stacks (winter  and summer, like January & August) for projection
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
apr <- read.table("glob_stack_month_apr_21_02_19.txt", h = T, sep = ";")
jul <- read.table("glob_stack_month_jul_21_02_19.txt", h = T, sep = ";")
oct <- read.table("glob_stack_month_oct_21_02_19.txt", h = T, sep = ";")
jan <- read.table("glob_stack_month_jan_21_02_19.txt", h = T, sep = ";")
colnames(apr)[7] <- "dSST"
colnames(jul)[7] <- "dSST"
colnames(oct)[7] <- "dSST"
colnames(jan)[7] <- "dSST"
apr <- apr[-which(apr$SSS < 20),] 
apr <- apr[-which(apr$Bathy > -175),] 
jul <- jul[-which(jul$SSS < 20),] 
jul <- jul[-which(jul$Bathy > -175),] 
oct <- oct[-which(oct$SSS < 20),] 
oct <- oct[-which(oct$Bathy > -175),] 
jan <- jan[-which(jan$SSS < 20),] 
jan <- jan[-which(jan$Bathy > -175),] 

### Load the remaining 8 climatologies for annual proj of diversity
annual <- TRUE
# annual <- FALSE
if(annual == TRUE) {
		
		feb <- read.table("glob_stack_month_feb_21_02_19.txt", h = T, sep = ";")
		mar <- read.table("glob_stack_month_mar_21_02_19.txt", h = T, sep = ";")
		may <- read.table("glob_stack_month_may_21_02_19.txt", h = T, sep = ";")
		jun <- read.table("glob_stack_month_jun_21_02_19.txt", h = T, sep = ";")
		aug <- read.table("glob_stack_month_aug_21_02_19.txt", h = T, sep = ";")
		sep <- read.table("glob_stack_month_sep_21_02_19.txt", h = T, sep = ";")
		nov <- read.table("glob_stack_month_nov_21_02_19.txt", h = T, sep = ";")
		dec <- read.table("glob_stack_month_dec_21_02_19.txt", h = T, sep = ";")

		colnames(feb)[7] <- "dSST"
		colnames(mar)[7] <- "dSST"
		colnames(may)[7] <- "dSST"
		colnames(jun)[7] <- "dSST"
		colnames(aug)[7] <- "dSST"
		colnames(sep)[7] <- "dSST"
		colnames(nov)[7] <- "dSST"
		colnames(dec)[7] <- "dSST"
	
		feb <- feb[-which(feb$SSS < 20),] 
		feb <- feb[-which(feb$Bathy > -175),] 
		mar <- mar[-which(mar$SSS < 20),] 
		mar <- mar[-which(mar$Bathy > -175),] 
		may <- may[-which(may$SSS < 20),] 
		may <- may[-which(may$Bathy > -175),] 
		jun <- jun[-which(jun$SSS < 20),] 
		jun <- jun[-which(jun$Bathy > -175),] 
		aug <- aug[-which(aug$SSS < 20),] 
		aug <- aug[-which(aug$Bathy > -175),] 
		sep <- sep[-which(sep$SSS < 20),] 
		sep <- sep[-which(sep$Bathy > -175),] 
		nov <- nov[-which(nov$SSS < 20),] 
		nov <- nov[-which(nov$Bathy > -175),] 
		dec <- dec[-which(dec$SSS < 20),] 
		dec <- dec[-which(dec$Bathy > -175),]

} # eo if loop : annual switch

future <- TRUE
rcp <- "rcp26"
if(future == TRUE) {
		
		setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/future_monthly_clims/diff/",rcp, sep = ""))
		future.jan <- read.table(paste("clim_2100-2071_",rcp,"_diff_Jan_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.feb <- read.table(paste("clim_2100-2071_",rcp,"_diff_Feb_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.mar <- read.table(paste("clim_2100-2071_",rcp,"_diff_Mar_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.apr <- read.table(paste("clim_2100-2071_",rcp,"_diff_Apr_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.may <- read.table(paste("clim_2100-2071_",rcp,"_diff_May_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.jun <- read.table(paste("clim_2100-2071_",rcp,"_diff_Jun_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.jul <- read.table(paste("clim_2100-2071_",rcp,"_diff_Jul_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.aug <- read.table(paste("clim_2100-2071_",rcp,"_diff_Aug_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.sep <- read.table(paste("clim_2100-2071_",rcp,"_diff_Sep_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.oct <- read.table(paste("clim_2100-2071_",rcp,"_diff_Oct_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.nov <- read.table(paste("clim_2100-2071_",rcp,"_diff_Nov_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		future.dec <- read.table(paste("clim_2100-2071_",rcp,"_diff_Dec_GFDL-ESM2M_24_07_19_v3.txt", sep = ""), h = T, sep = "\t")
		setwd(WD)
		
} # eo if loop : future vs baseline switch


### 25/09/18: Sometimes the server shuts off and interrupts the niche modelling...need to identify the species that are already modelled and projected and remove them from the 'species' vector
#  setwd(second.wd)
#  dir() # just need to replace the dot by an underscore and there you go
#  donespp <- gsub("\\.", "_", dir()) # But also need to add brackets for some annoying species names
#  for(i in 1:length(donespp)) {
#     	sp <- donespp[i]
#     	if( length( strsplit(sp, "_", fixed = TRUE)[[1]] ) == 3) {
#     		  # Then add brackets around the second piece
#     		sp <- paste(strsplit(sp, "_", fixed = TRUE)[[1]][1],"_(",
#     				strsplit(sp, "_", fixed = TRUE)[[1]][2],")_",
#     				strsplit(sp, "_", fixed = TRUE)[[1]][3], sep = "")
#
#     		donespp[i] <- sp
#    		}
#  }  #eo for loop
#  # OK Now remove the 'donespp' from 'species'
#  species2 <- species[!(species %in% donespp)]
# #
# # # 25/09/18: And SOMETIMES, the models wouldn't get projected even though they're done !
#  setwd(second.wd)
#  modelled_spp <- dir()
#  list <- lapply(modelled_spp, function(s) {
#     				 # Go to species dir and list the number of elements
#     				 # First if loop to add brackets when necessary...
#     				if( length( strsplit(s, "_", fixed = TRUE)[[1]] ) == 3) {
#     					 # Then add brackets around the second piece
#     					s <- paste(strsplit(s, "_", fixed = TRUE)[[1]][1],"_(",
#     							strsplit(s, "_", fixed = TRUE)[[1]][2],")_",
#     							strsplit(s, "_", fixed = TRUE)[[1]][3], sep = "")
#     				}
#
#     				setwd( paste(second.wd, "/", str_replace_all(s, "_", "."), "/", sep = "") )
#     				n <- length(dir())
#     				 # If there are more than 5 elements (n >= 5), then it means the the projections have already been carried out sucessfully
#     				if(n >= 5) {
#     					rm(n)
#     				} else {
#     					return(s)
#     				}  # eo seond if else loop
# } ) # eo lapply
#
#  unproj_spp <- do.call(rbind, list)[,1]
#  #  So species2model is a combination of species2 & unproj_spp
#  species2model <- gsub("\\.", "_", c(species2, unproj_spp))
# # species2model <- gsub("\\.", "_", unproj_spp[13:24])
# # species2model <- gsub("\\.", "_", unproj_spp)
#  rm(modelled_spp, donespp)

setwd(WD)

### Set modelling options
myBiomodOption <- BIOMOD_ModelingOptions(
	
						GLM = list( type = 'quadratic',
    								interaction.level = 0,
   			 						myFormula = NULL,
    								test = 'AIC',
    								family = binomial("logit"),
    								mustart = 0.5,
    								control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),

						GAM = list(algo = 'GAM_mgcv',
								type = 's_smoother',
								k = 5,
								interaction.level = 0,
								myFormula = NULL,
								family = binomial("logit"),
								method = 'GCV.Cp',
								optimizer = c('outer','newton'),
								select = FALSE,
								knots = NULL,
								paraPen = NULL,
								control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
								, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
								, rank.tol = 1.49011611938477e-08
								, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
								, optim = list(factr=1e+07)
								, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) ), 
						
						RF = list(do.classif = T, 
							ntree = 750, 
							nodesize = 10 )		
								
) # eo modelling options

SpeciesNicheModelling <- function(sp = species) {
	
							# Get the data
							setwd(WD)
							message(paste("Modelling ",sp, " =========================================================", sep = ""))
							data <- read.table(paste("data_",strat,"_",sp,".txt", sep = ""), h = T, sep = ";")
							
							data2 <- na.omit( data[,c(vars,"obs")] )
							n <- nrow( data2[data2$obs == 1,] )
							rm(data2)
							
							# If n >= 85, continue, otherwise got to next species
							if(n >= 75) {
								
	  			  				### Initialisation: data formatting
								myRespName <- str_replace_all(sp, "_", ".")
								myRespName <- gsub("\\(|\\)", "", myRespName)
							
								# the presence/absences data for our species 
								myResp <- as.numeric(data$obs)

								# the XY coordinates of species data
								myRespXY <- data[,c("x","y")]

								# Environmental variables
								# colnames(data)[c(21)] <- c("deltaT")
								myExpl <- data[,vars]

								# weights vector
								weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
								# formatage des données
	  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
	                                      	 			 	expl.var = myExpl,
	                                       			  	 	resp.xy = myRespXY,
	                                       			  	 	resp.name = myRespName 
								)
							
	 			   				### Modelling
								setwd(second.wd)
								myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
								                                     models = SDMs, 
								                                     models.options = myBiomodOption, 
								                                     NbRunEval = length(eval_runs), 
								                                     DataSplit = 80, 
								                                     Yweights = weights,
								                                     VarImport = 0, 
								                                     models.eval.meth = c('ROC','TSS','KAPPA'),
								                                     SaveObj = FALSE,
								                                     do.full.models = FALSE
								) # eo modelling					
  
								### Get evaluation scores
								scores <- data.frame(get_evaluations(myBiomodModelOut))
								# Use lapply to summarize all scores (10 scores per 6 SDMs) AND PROBABILITY THRESH
								scores.ls <- lapply(SDMs, function(sdm) {
												# Retrieve all scores for the sdm 's' 
												tss <- scores["TSS",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
												tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
												colnames(tss) <- paste(sdm, eval_runs, sep = "_")
												colnames(tss) <- paste(sdm, eval_runs, sep = "_")
												tss <- t(tss)
												tss_cutoffs <- t(tss_cutoffs)
												colnames(tss_cutoffs) <- "Cutoff_TSS"
											
												# Same with AUC
												auc <- scores["ROC",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
												auc_cutoffs <- scores["ROC",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
												colnames(auc) <- paste(sdm, eval_runs, sep = "_")
												colnames(auc_cutoffs) <- paste(sdm, eval_runs, sep = "_")
												auc <- t(auc)
												auc_cutoffs <- t(auc_cutoffs)
												colnames(auc_cutoffs) <- "Cutoff_AUC"
											
												# Same with Cohen's KAPPA
												kappa <- scores["KAPPA",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
												kappa_cutoffs <- scores["KAPPA",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
												colnames(kappa) <- paste(sdm, eval_runs, sep = "_")
												colnames(kappa_cutoffs) <- paste(sdm, eval_runs, sep = "_")
												kappa <- t(kappa)
												kappa_cutoffs <- t(kappa_cutoffs)
												colnames(kappa_cutoffs) <- "Cutoff_Kappa"
											
												return(cbind(tss, tss_cutoffs, auc, auc_cutoffs, kappa, kappa_cutoffs))
								} ) # eo lapply
								scores.tbl <- data.frame(do.call(rbind, scores.ls))
								# scores.tbl
								rm(scores, scores.ls)
							
	  			  				### Make projections (they will be printed in the species folder)
								message(paste("Projecting ",sp, " =========================================================", sep = ""))
								# Project niches in April conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = apr[,vars],
													proj.name = paste("projection",sp,"apr", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection
								# Project niches in July conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = jul[,vars],
													proj.name = paste("projection",sp,"jul", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection			
								# Project niches in Octpber conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = oct[,vars],
													proj.name = paste("projection",sp,"oct", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection			
								# Project niches in January conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = jan[,vars],
													proj.name = paste("projection",sp,"jan", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection
							
								### Create a switch for annual projections (annual = T/F)
								if( annual == TRUE ) {
									# Feb
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = feb[,vars],
														proj.name = paste("projection",sp,"feb", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Mar
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = mar[,vars],
														proj.name = paste("projection",sp,"mar", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# May
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = may[,vars],
														proj.name = paste("projection",sp,"may", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Jun
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = jun[,vars],
														proj.name = paste("projection",sp,"jun", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Aug
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = aug[,vars],
														proj.name = paste("projection",sp,"aug", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection	
									# Sep
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = sep[,vars],
														proj.name = paste("projection",sp,"sep", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Nov
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = nov[,vars],
														proj.name = paste("projection",sp,"nov", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Dec
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = dec[,vars],
														proj.name = paste("projection",sp,"dec", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									
								} else { 
									message(paste("============== NOT PERFORMING ANNUAL PROJECTIONS ==============", sep = ""))
								} # eo if else loop
								
								if(future == TRUE) {
									# Project niches in future April conditions
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.apr[,vars],
														proj.name = paste("projection",sp,"apr_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Project niches in future July conditions
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.jul[,vars],
														proj.name = paste("projection",sp,"jul_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection			
									# Project niches in future October conditions
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.oct[,vars],
														proj.name = paste("projection",sp,"oct_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection			
									# Project niches in future January conditions
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.jan[,vars],
														proj.name = paste("projection",sp,"jan_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Feb
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.feb[,vars],
														proj.name = paste("projection",sp,"feb_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Mar
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.mar[,vars],
														proj.name = paste("projection",sp,"mar_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# May
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.may[,vars],
														proj.name = paste("projection",sp,"may_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Jun
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.jun[,vars],
														proj.name = paste("projection",sp,"jun_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Aug
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.aug[,vars],
														proj.name = paste("projection",sp,"aug_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection	
									# Sep
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.sep[,vars],
														proj.name = paste("projection",sp,"sep_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Nov
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.nov[,vars],
														proj.name = paste("projection",sp,"nov_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Dec
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = future.dec[,vars],
														proj.name = paste("projection",sp,"dec_2100_GFDLESM2M", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
								} else {
									message(paste("============== NOT PERFORMING FUTURE PROJECTIONS ==============", sep = ""))
								}

								### Save evaluation scores
								setwd(paste(WD,"/","niche.modelling_future_",rcp,"/eval_scores_",p,"/", sep = ""))
								save(scores.tbl, file = paste("eval_scores_",sp,".Rdata", sep = ""))
								setwd(WD)  			
								
							} else {
							
								message(paste("NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling

pools <- paste("p",c(1:4),sep = "")

# For testing:
#p <- "p1"

for(p in pools) {
	
	# Set pool of variables
	if(p == "p1") {
		vars <- c("SST","dSST","dO2","logNO3","logChl")
	} else if (p == "p2") {
		vars <- c("SST","dSST","dO2","logSiO2","logChl")
	} else if (p == "p3") {
		vars <- c("SST","dSST","dO2","logSiO2","logChl","Nstar")
	} else if (p == "p4") {
		vars <- c("SST","dSST","dO2","logNO3","logChl","Sistar")
	} 

	require("parallel")
	message(paste(" ", sep = ""))	
	message(paste("RUNNING ZOOPLANKTON PROJECTIONS FOR VARIABLE POOL ", p, sep = ""))
	message(paste(" ", sep = ""))
	
	setwd(paste(WD,"/niche.modelling_future_",rcp,"/",p,"/", sep = ""))
	second.wd <- getwd()
	setwd(WD)
	
	mclapply(X = species, SpeciesNicheModelling, mc.cores = 26)
	
} # eo for p in pools



### ==============================================================

### And now for phytoplankton 
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
WD <- getwd()
# Niche modelling dir
setwd(paste(WD,"/niche.modelling_future_",rcp,"/",p,"/", sep = ""))
second.wd <- getwd()
setwd(WD)
# Vector of SDMs
SDMs <- c('GLM','GAM','RF','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
### Vector of species names
files <- dir()[grep(".txt",dir())]
files <- files[-grep("skillz",files)]
# Extract species name sfrom file names
species <- stringr::str_replace_all(string = files, pattern = paste("data_",strat,"_", sep = ""), replacement = "")
species <- stringr::str_replace_all(string = species, pattern = ".txt", replacement = "")

pools <- paste("p",c(1:4),sep = "")

setwd(WD)

for(p in pools) {
	
	if(p == "p1") {
		vars <- c("SST","dSST","PAR","logNO3","logChl","Nstar")
	} else if (p == "p2") {
		vars <- c("SST","dSST","PAR","logSiO2","logChl","Nstar")
	} else if (p == "p3") {
		vars <- c("SST","dSST","PAR","logNO3","logChl","Nstar","Sistar")
	} else if (p == "p4") {
		vars <- c("SST","dSST","logChl","PAR","logNO3")
	} # eo else if loop
	 
	require("parallel")
	message(paste(" ", sep = ""))	
	message(paste("RUNNING PHYTOPLANKTON PROJECTIONS FOR VARIABLE POOL ", p, sep = ""))
	message(paste(" ", sep = ""))
	
	setwd(paste(WD,"/niche.modelling_future_",rcp,"/",p,"/", sep = ""))
	second.wd <- getwd()
	setwd(WD)
	
	mclapply(X = species, SpeciesNicheModelling, mc.cores = 26)
	
} # eo for pool





### ==============================================================================================================================
### ==============================================================================================================================
### ==============================================================================================================================

