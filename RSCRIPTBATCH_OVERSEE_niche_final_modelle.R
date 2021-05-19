
### ==============================================================================================================================

library("raster")
library("rgdal")
library("sp")
library("rJava") 
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")

### ==============================================================================================================================


### Set the dataset and the background selection strategy
#set <- "v9v3.2"
set <- "v9v3.1"
strat <- "total"
# Or: 
#strat <- "total2"

### Set the working directories
if (strat == "total") {
	
	setwd( paste(getwd(),"/","species_data_",set,"/","total_background/", sep = "") )
	WD <- getwd()
	# Niche modelling dir
	setwd(paste(WD,"/","niche.modelling","/", sep = ""))
	second.wd <- getwd()
	
} else if ( strat == "total2" ) {
	
	setwd( paste(getwd(),"/","species_data_",set,"/","total_background_rar/", sep = "") )
	WD <- getwd()
	setwd(paste(WD,"/","niche.modelling","/", sep = ""))
	second.wd <- getwd()
	
} # eo if else loop

setwd(WD)

### List of SDMs
SDMs <- c('GLM','GAM','RF','ANN')
#sdm <- "MAXENT.Phillips"
### List of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 

### List of tets species names
#files <- dir()[grep(".txt",dir())]
# Extract species name sfrom file names
#species <- stringr::str_replace_all(string = files, pattern = paste("data_",strat,"_", sep =""), replacement = "")
#species <- stringr::str_replace_all(string = species, pattern = ".txt", replacement = "")
# species

### To identify the species that have not been modelled and those that have been but for which eval scores are missing
#setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_",set,"/total_background_rar/eval_scores", sep = "") )
#donespp <- dir()
# Remove .Rdata and eval_scores from ch string
#donespp <- str_replace_all(donespp, ".Rdata", "")
#donespp <- str_replace_all(donespp, "eval_scores_", "")

#for(i in 1:length(donespp)) {
    	#sp <- donespp[i]
    	#if( length( strsplit(sp, "_", fixed = TRUE)[[1]] ) == 3) {
    		  # Then add brackets around the second piece
    		#sp <- paste(strsplit(sp, "_", fixed = TRUE)[[1]][1],"_(",
    				#strsplit(sp, "_", fixed = TRUE)[[1]][2],")_",
    				#strsplit(sp, "_", fixed = TRUE)[[1]][3], sep = "")

    		#donespp[i] <- sp
			#}
#}  # eo for loop
#  OK Now remove the 'donespp' from 'species'
#species2 <- species[!(species %in% donespp)]



### Load the env stacks (winter and summer, like January & August) for projection
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
jan <- read.table("glob_stack_month_jan_23_10_18.txt", h = T, sep = ";")
feb <- read.table("glob_stack_month_feb_23_10_18.txt", h = T, sep = ";")
mar <- read.table("glob_stack_month_mar_23_10_18.txt", h = T, sep = ";")
apr <- read.table("glob_stack_month_apr_23_10_18.txt", h = T, sep = ";")
may <- read.table("glob_stack_month_may_23_10_18.txt", h = T, sep = ";")
jun <- read.table("glob_stack_month_jun_23_10_18.txt", h = T, sep = ";")
jul <- read.table("glob_stack_month_jul_23_10_18.txt", h = T, sep = ";")
aug <- read.table("glob_stack_month_aug_23_10_18.txt", h = T, sep = ";")
sep <- read.table("glob_stack_month_sep_23_10_18.txt", h = T, sep = ";")
oct <- read.table("glob_stack_month_oct_23_10_18.txt", h = T, sep = ";")
nov <- read.table("glob_stack_month_nov_23_10_18.txt", h = T, sep = ";")
dec <- read.table("glob_stack_month_dec_23_10_18.txt", h = T, sep = ";")
# Vector of variables name
vars <- c("SST","deltaT","SSS","PAR","MLD1","logChl","logSiO2")

### Set modelling options
myBiomodOption <- BIOMOD_ModelingOptions(
	
						GLM = list( type = 'quadratic',
    								interaction.level = 0,
   			 						myFormula = NULL,
    								test = 'AIC',
    								family = binomial("logit"),
    								mustart = 0.5,
    								control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),
				
						 MAXENT.Phillips = list(path_to_maxent.jar = second.wd,
									maximumiterations = 500,
									visible = FALSE,
									linear = TRUE,
									quadratic = TRUE,
									product = TRUE,
									threshold = FALSE,
									hinge = TRUE,
									lq2lqptthreshold = 80,
									l2lqthreshold = 10,
									hingethreshold = 15,
									beta_threshold = -1,
									beta_categorical = -1,
									beta_lqp = -1,
									beta_hinge = -1,
									defaultprevalence = 0.5),

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
									ntree = 5000, 
									nodesize = 10
						)		
								
) # eo modelling options

setwd(WD)
### Define the function you will apply in parallel

# sp <- "Calanus_finmarchicus"

SpeciesNicheModelling <- function(sp = sp) {				  
							
							# Get the data
							setwd(WD)
							message(paste("Modelling ",sp, " =========================================================", sep = ""))
							data <- read.table(paste("data_",strat,"_",sp,".txt", sep = ""), h = T, sep = ";")
							# summary(data[data$obs == 1,])
							
  			  				### Initialisation: data formatting ==================================================================================
							myRespName <- str_replace_all(sp, "_", ".")
							myRespName <- gsub("\\(|\\)", "", myRespName)
							
							# the presence/absences data for our species 
							myResp <- as.numeric(data$obs)

							# the XY coordinates of species data
							myRespXY <- data[,c("x","y")]

							# Environmental variables
							colnames(data)[c(20)] <- c("deltaT")
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
  
  							### Save evaluation scores ==================================================================================
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
							} 
							) # eo lapply
							scores.tbl <- data.frame(do.call(rbind, scores.ls))
							# scores.tbl
							rm(scores, scores.ls)
							
  			  				### Make projections (they will be printed in the species folder) ============================================
							message(paste("Projecting ",sp, " =========================================================", sep = ""))
							
							### Project niches in jan conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = jan[,vars],
												proj.name = paste("projection",sp,"jan", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection

							### Project niches in feb conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = feb[,vars],
												proj.name = paste("projection",sp,"feb", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in mar conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = mar[,vars],
												proj.name = paste("projection",sp,"mar", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in apr conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = apr[,vars],
												proj.name = paste("projection",sp,"apr", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in may conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = may[,vars],
												proj.name = paste("projection",sp,"may", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in jun conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = jun[,vars],
												proj.name = paste("projection",sp,"jun", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in jul conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = jul[,vars],
												proj.name = paste("projection",sp,"jul", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in aug conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = aug[,vars],
												proj.name = paste("projection",sp,"aug", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in sep conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = sep[,vars],
												proj.name = paste("projection",sp,"sep", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in oct conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = oct[,vars],
												proj.name = paste("projection",sp,"oct", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in nov conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = nov[,vars],
												proj.name = paste("projection",sp,"nov", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
							### Project niches in dec conditions
							myBiomodProj <- BIOMOD_Projection(
												modeling.output = myBiomodModelOut,
												new.env = dec[,vars],
												proj.name = paste("projection",sp,"dec", sep = "_"),
												selected.models = 'all',
												binary.meth = 'TSS',
												compress = 'xz',
												clamping.mask = FALSE 
							) # eo projection
							
  			  				
							### Save evaluation scores ===================================================================================
							setwd(paste(WD,"/","eval_scores","/", sep = ""))
							save(scores.tbl, file = paste("eval_scores_",sp,".Rdata", sep = ""))
							setwd(WD)  		
  
} # eo NicheModelling


library("parallel")
mclapply(X = sp, SpeciesNicheModelling, mc.cores = 25)


### ==============================================================================================================================

