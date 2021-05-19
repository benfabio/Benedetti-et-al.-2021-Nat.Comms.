
##### 25/08/2018: R Script to match the species observations to values of the global monthly climatologies of the 21 predictors chosen

### Aims to
#	- load the monthly 1°x1° climatologies (txt files)
#	- convert env data frames to raster stacks to 
#	- load the v8 datasets (v8v5.1v3.1 and v8v5.1v3.2)
#	- use the extract() function of 'raster' to match the species obs to env predictors
#	- also get the 1°x1° coordinates
#	- examine distribution of the matched env data (boxplots/ frequency histograms) to assess how much is covered by the obs
#	- prepare data for psAbs generation using the target-group technique


### Latest update: 20/03/2019

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
library("ggplot2")
library("RColorBrewer")
library("viridis")


# ==============================================================================================================================

WD <- getwd()

##### 1°) Load the env data and convert to raster stacks
setwd(paste(WD,"/","env_predictors","/","global_monthly_clims_1d","/", sep = ""))
cl <- dir()[grep("23_10_18", dir())] ; cl
# Load each clim, turn into raster stack and store them in a list
clims <- lapply(cl, function(c) {
			# read table
			message(paste("Reading ",c, sep = ""))
			dat <- read.table(c, sep = ";", h = TRUE)
			# convert to stack
			coordinates(dat) <- ~ x + y
			gridded(dat) <- TRUE
			# coerce to raster
			stak <- stack(dat)
			crs(stak) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
			# return
			return(stak)
		}
) # eo lapply

# Examine resulting list of stacks
# str(clims) # class(clims)
# To extract one clim from the list 'clims'
# Change thr name of the bathymetry layer
#names(clims[[1]]) <- "bathymetry"
names(clims[[2]])
names(clims[[10]])
length(clims) # 12 elements, one per month, containing 24 fields 
for(i in c(1:24)) {
		message(paste(names(clims[[i]]), sep = ""))
		message(paste(".  ", sep = ""))
} # eo for loop

# 1 = Bathy, 2 = Chl, 3 = deltaT, 4 = logChla, 5 = logMLD, 6 = logNO3, 7 = logPO4, 8 = logSiO2
# 9 = MLD1, 10 = MLPAR, 11 = NO3, 12 = Nstar, 13 = PAR, 14 = PO4, 15 = Pstar, 16 = SiO2
# 17 = Sistar, 18 = SLA, 19 = SSS, 20 = SST, 21 = Wind speed




##### 2°) For each v8 dataset, load and use extract() to mact each species occurrence to env layers

### A°) v8-5.1v3.1 datasets
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.1/")
dir()
### Use double for loop to load data and then match each 21 dimensions ?
groups <- dir() ; groups
layers <- c(1:12) # 1 per month

# For testing:
# g <- "Copepoda_06_06_18.Rdata"
# l <- 1 # January 

for(g in groups) {
		
		# Useless message
		message(paste("Matching ", g, sep = ""))
		message(paste("                 ", sep = ""))
		
		# Load species data
		data <- get(load(g))
		# colnames(data)
		
		# Just to match annual data (bathy & dSST) + 1d coordinates
		envdata <- clims[[1]]
		# To extract a sub layer ?
		# subset(x = envdata, subset = 'SST')
		# To extract the 1d coordinates corresponding to EACH occurrence
		cells <- cellFromXY(subset(x = envdata, subset = 1), as.matrix(data[,c("x","y")]))
		coords <- data.frame(xyFromCell(subset(x = envdata, subset = 1), cells), 
				   value = raster::extract(subset(x = envdata, subset = 1), cells) ) 
		# 
		if( nrow(coords) == nrow(data) ) {
			data$xbin_1d <- coords$x
			data$ybin_1d <- coords$y
		} # eo if else loop
		
		# Treat the 2 non monthly clims (bathy & deltaT) differently from the others
		# Add bathymetry and dSST (annual)
		data$Bathy <- raster::extract(x = subset(x = envdata, subset = 'Bathy'), y = data[,c("x","y")], method = 'bilinear')
		data$dSST <- raster::extract(x = subset(x = envdata, subset = 'deltaT'), y = data[,c("x","y")], method = 'bilinear')					
		# summary(data)
		# nrow(data[data$Bathy > 0,])
		# Filter the points that were given positive bathymetry values because it is a 1° product
		data <- data[data$Bathy < 0,]
		
		# Create empty vectors for each predictor (20 in total since you already provided Bathy & dSST)
		vars <- names(clims[[2]])[c(2,6:24)]
		data[vars] <- NA

		### And fill those (according to month) in a for loop
		for(l in layers) {
				
				# Get the monthly clims 
				message(paste("Matching month ", l, sep = ""))
				envdata <- clims[[l]]
				
				# Fill it according to month
				data[which(data$month == l),"Chl"] <- raster::extract(x = subset(x = envdata, subset = 'Chl'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logChl"] <- raster::extract(x = subset(x = envdata, subset = 'logChl'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logMLD"] <- raster::extract(x = subset(x = envdata, subset = 'logMLD'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logNO3"] <- raster::extract(x = subset(x = envdata, subset = 'logNO3'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logPO4"] <- raster::extract(x = subset(x = envdata, subset = 'logPO4'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logSiO2"] <- raster::extract(x = subset(x = envdata, subset = 'logSiO2'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"MLD1"] <- raster::extract(x = subset(x = envdata, subset = 'MLD1'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"MLPAR"] <- raster::extract(x = subset(x = envdata, subset = 'MLPAR'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Nstar"] <- raster::extract(x = subset(x = envdata, subset = 'Nstar'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"NO3"] <- raster::extract(x = subset(x = envdata, subset = 'NO3'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"PAR"] <- raster::extract(x = subset(x = envdata, subset = 'PAR'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"PO4"] <- raster::extract(x = subset(x = envdata, subset = 'PO4'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Pstar"] <- raster::extract(x = subset(x = envdata, subset = 'Pstar'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SiO2"] <- raster::extract(x = subset(x = envdata, subset = 'SiO2'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Sistar"] <- raster::extract(x = subset(x = envdata, subset = 'Sistar'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SLA"] <- raster::extract(x = subset(x = envdata, subset = 'SLA'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SSS"] <- raster::extract(x = subset(x = envdata, subset = 'SSS'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SST"] <- raster::extract(x = subset(x = envdata, subset = 'SST'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Wind"] <- raster::extract(x = subset(x = envdata, subset = 'Wind'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"pCO2"] <- raster::extract(x = subset(x = envdata, subset = 'pCO2'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')

		} # eo second for loop
		
		# Check...
		# colnames(data)
		# summary(data)
		# str(data)
		
		# Save in v9 dir()
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
		save(data, file = paste(str_replace_all(g, "_06_06_18.Rdata", "_matched_28_11_18.Rdata"), sep = "") )
		
		# Clean and go back to groups' dir
		rm(data)
		gc()
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.1/")
	
} # eo first for loop
 
 
# -------------------------------------------------------------------------------------------------------------------------


### B°) v8-5.1v3.2 datasets
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
dir()

### Use double for loop to load data and then match each 21 dimensions ?
groups <- dir() ; groups
layers <- c(1:12)

# For testing:
# g <- "Copepoda_06_06_18.Rdata"
# l <- 20 # SST 

for(g in groups) {
		
		# Useless message
		message(paste("Matching ",g, sep = ""))
		message(paste("                 ", sep = ""))
		
		# Load species data
		data <- get(load(g))
		# colnames(data)
		
		envdata <- clims[[l]]
		# To extract a sub layer ?
		# subset(x = envdata, subset = 'SST')
		# To extract the 1d coordinates corresponding to EACH occurrence
		cells <- cellFromXY(subset(x = envdata, subset = 1), as.matrix(data[,c("x","y")]))
		coords <- data.frame(xyFromCell(subset(x = envdata, subset = 1), cells), 
				   value = raster::extract(subset(x = envdata, subset = 1), cells) ) 
		# 
		if( nrow(coords) == nrow(data) ) {
			data$xbin_1d <- coords$x
			data$ybin_1d <- coords$y
		} # eo if else loop
		
		# Treat the 2 non monthly clims (bathy & deltaT) differently from the others
		# Add bathymetry and dSST (annual)
		data$Bathy <- raster::extract(x = subset(x = envdata, subset = 'Bathy'), y = data[,c("x","y")], method = 'bilinear')
		data$dSST <- raster::extract(x = subset(x = envdata, subset = 'deltaT'), y = data[,c("x","y")], method = 'bilinear')					
		# summary(data)
		# nrow(data[data$Bathy > 0,])
		# Filter the points that were given positive bathymetry values because it is a 1° product
		data <- data[data$Bathy < 0,]
		
		# Create empty vectors for each predictor (20 in total since you already provided Bathy & dSST)
		vars <- names(clims[[2]])[c(2,6:24)]
		data[vars] <- NA

		### And fill those (according to month) in a for loop
		for(l in layers) {
				
				# Get the monthly clims 
				message(paste("Matching month ", l, sep = ""))
				envdata <- clims[[l]]
				
				# Fill it according to month
				data[which(data$month == l),"Chl"] <- raster::extract(x = subset(x = envdata, subset = 'Chl'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logChl"] <- raster::extract(x = subset(x = envdata, subset = 'logChl'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logMLD"] <- raster::extract(x = subset(x = envdata, subset = 'logMLD'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logNO3"] <- raster::extract(x = subset(x = envdata, subset = 'logNO3'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logPO4"] <- raster::extract(x = subset(x = envdata, subset = 'logPO4'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"logSiO2"] <- raster::extract(x = subset(x = envdata, subset = 'logSiO2'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"MLD1"] <- raster::extract(x = subset(x = envdata, subset = 'MLD1'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"MLPAR"] <- raster::extract(x = subset(x = envdata, subset = 'MLPAR'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Nstar"] <- raster::extract(x = subset(x = envdata, subset = 'Nstar'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"NO3"] <- raster::extract(x = subset(x = envdata, subset = 'NO3'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"PAR"] <- raster::extract(x = subset(x = envdata, subset = 'PAR'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"PO4"] <- raster::extract(x = subset(x = envdata, subset = 'PO4'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Pstar"] <- raster::extract(x = subset(x = envdata, subset = 'Pstar'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SiO2"] <- raster::extract(x = subset(x = envdata, subset = 'SiO2'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Sistar"] <- raster::extract(x = subset(x = envdata, subset = 'Sistar'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SLA"] <- raster::extract(x = subset(x = envdata, subset = 'SLA'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SSS"] <- raster::extract(x = subset(x = envdata, subset = 'SSS'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"SST"] <- raster::extract(x = subset(x = envdata, subset = 'SST'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"Wind"] <- raster::extract(x = subset(x = envdata, subset = 'Wind'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')
				data[which(data$month == l),"pCO2"] <- raster::extract(x = subset(x = envdata, subset = 'pCO2'), 
				y = data[which(data$month == l),c("x","y")], method = 'bilinear')

		} # eo second for loop
		
		# Check...
		# colnames(data)
		# summary(data)
		# unique(data$logSiO2)
		# str(data)
		
		# Save in v9 dir()
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.2/")
		save(data, file = paste(str_replace_all(g, "_06_06_18.Rdata", "_matched_28_11_18.Rdata"), sep = "") )
		
		# Clean and go back to groups' dir
		rm(data)
		gc()
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v8-5.1v3.2/")
	
} # eo first for loop
 



# ==============================================================================================================================


##### 27/08/2018: Examine how the env space is covered by the biological observations.
### Plot the two distributions (total and covered) on the same plot to visually assess how much is covered

# First, rbind all the env data matched with the species obs
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.2/")
d <- get(load(dir()[1]))
names <- colnames(d)[c(1:6,10:16,25:47)]
names
# Get all obs data
matched <- lapply(dir(), function(d) {
				data <- get(load(d))
				message(paste(d, sep = ""))
				return(subset(data, select = names ))	
}
) # eo lapply
match <- do.call(rbind, matched)
dim(match) # 766 033 obs for 3.1 ; 1 212 771 for 3.2
rm(matched) ; gc()
# Change 2 colnames to match vector of predictors below
colnames(match)[c(16,18)] <- c("Bathy","deltaT")

### For each predictor, load the data, melt it so you have all values in a vector, add a condition, then rbind with matched values that are given another condition, and plot the distribution of the predictor according to the two conditions with dodged histograms  
predictors <- c("SST","dSST","SSS","Bathy","Wind","SLA","MLD1","logMLD","PAR","MLPAR","pCO2",
				"Chl","logChl","NO3","logNO3","PO4","logPO4","SiO2","logSiO2","Nstar","Pstar","Sistar")

#
for(p in predictors) {
		
		# Go to predictors directory
		message(paste(p, sep = ""))
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
		# Read p
		pred <- read.table(paste(p,"_stack_1d.txt", sep = ""), h = T, sep = ";")
		# Melt
		m_pred <- melt(pred, id.vars = c("x","y"))
		# summary(m_pred)
		m_pred$cond <- "total" # provide a condition for plotting later
		colnames(m_pred)[4] <- p
		m_pred <- m_pred[,c(1,2,4,5)]
		
		# Get the matched data for this same pred
		matched <- match[,c("xbin_1d","ybin_1d",p)]
		matched$cond <- "matched"
		colnames(matched)[c(1:2)] <- c("x","y")
		
		# rbind for plot
		tot <- rbind(matched, m_pred)
		
		if( p == "Bathy" ) {
				tot <- tot[tot$Bathy < 0,]
		} # eo if loop
		
		### Compute % overlap between the areas of the 2 density kernels
		require("overlapping")
		over <- list(X1 = na.omit(tot[which(tot$cond == "total"),p]), X2 = na.omit(tot[which(tot$cond == "matched"),p]))
		areaoverlap <- round(overlap(over, plot = F)$OV, 3)
		# area of overlap is in @ overlap(over, plot = F)$OV
		
		density <- ggplot(data = tot, aes(x = tot[,3], fill = factor(cond))) + 
					geom_density(alpha = .4) + theme_bw() + xlab(paste(p," (overlap = ",areaoverlap,")", sep = ""))
					
		violin <- ggplot(data = tot, aes(x = factor(cond), y = tot[,3], fill = factor(cond))) + 
					geom_violin() + theme_bw() + xlab(paste(p," (overlap = ",areaoverlap,")", sep = "")) + ylab(p)
					
		# Go save plots			
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/")
		ggsave(plot = density, filename = paste("density_matched_",p,".pdf", sep = ""), dpi = 300, height = 5, width = 7)
		ggsave(plot = violin, filename = paste("violin_matched_",p,".pdf", sep = ""), dpi = 300, height = 5, width = 7)
	
} # eo for loop
	
### Using the for loop above, compute average density overlap
overlaps <- lapply(predictors, function(p) {
					
					message(paste(p, sep = ""))
					setwd("/UP_home/fabioben/Desktop/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
					# Read p
					pred <- read.table(paste(p,"_stack_1d.txt", sep = ""), h = T, sep = ";")
					# Melt
					m_pred <- melt(pred, id.vars = c("x","y"))
					# summary(m_pred)
					m_pred$cond <- "total" # provide a condition for plotting later
					colnames(m_pred)[4] <- p
					m_pred <- m_pred[,c(1,2,4,5)]
	
					# Get the matched data for this same pred
					matched <- match[,c("xbin_1d","ybin_1d",p)]
					matched$cond <- "matched"
					colnames(matched)[c(1:2)] <- c("x","y")
	
					# rbind for plot
					tot <- rbind(matched, m_pred)
	
					if( p == "Bathy" ) {
							tot <- tot[tot$Bathy < 0,]
					} # eo if loop
	
					### Compute % overlap between the areas of the 2 density kernels
					require("overlapping")
					over <- list(X1 = na.omit(tot[which(tot$cond == "total"),p]), X2 = na.omit(tot[which(tot$cond == "matched"),p]))
					areaoverlap <- round(overlap(over, plot = F)$OV, 3)
					# area of overlap is in @ overlap(over, plot = F)$OV
					# And compute number of NAs in matched
					nNA <- nrow(matched[is.na(matched[,p]),])
					# and % of occurrences lost because of this
					per <- nNA / nrow(matched)
					# return
					return(data.frame(overlap = areaoverlap, nNA = nNA, percentage = round(per,3)*100 ) )
	
}) # eo lapply

# rbind
kde <- data.frame(do.call(rbind, overlaps))
rownames(kde) <- predictors
kde
# Compute mean + sd
mean(kde$overlap); sd(kde$overlap)
# rm(overlaps) ; gc()

# For 3.1 -> 0.8398571 (0.84) ± 0.065
# For 3.2 -> 0.7961905 (0.80) ± 0.097


### From 'match' object, plot correlogram between the 21 predictors 

# v9v8v5v3.1
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd")
d <- read.table(dir()[1], sep = ";", h = T)
# colnames(d)
#d <- read.table("data_total_Calanus_finmarchicus.txt", sep = ";", h = T)
names <- colnames(d)[c(1:6,19:39,43:45)]
names
# Get all obs data
require("parallel")
matched <- mclapply( dir()[grep("data_total",dir())], function(d) {
				data <- read.table(d, sep = ";", h = T)
				message(paste(d, sep = ""))
				return(data[data$obs == 1,names])	
			}, mc.cores = 23
) # eo lapply
match <- do.call(rbind, matched)
dim(match) # 764 159 obs for 3.1 
rm(matched) ; gc()

### Plot correlogram
library("corrplot")
colnames(match)
p.mat <- cor.mtest(na.omit(match[,c(7:30)]))
p.mat 
# Correlation matrix
M <- cor(na.omit(match[,c(7:30)]), method = "spearman")
M

### Plot correlogram
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/")
pdf("corrplot_zoo_v9v5.1.3.1_20_03_19_Ponly.pdf", height = 15, width = 15)
corrplot(M, method = "color", col = col(200), type = "upper", order = "hclust", 
         	   	addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE )
dev.off()


# ----------------------------------------------

# v9v8v5v3.2
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.2/")

# Get all obs data
matched <- lapply(dir(), function(d) {
				data <- get(load(d))
				message(paste(d, sep = ""))
				return(subset(data, select = names ))	
}
) # eo lapply
match <- do.call(rbind, matched)
dim(match) # 1 184 690 for 3.2
rm(matched) ; gc()

### Plot correlogram
library("corrplot")
p.mat <- cor.mtest(na.omit(match[,c(14:35)]))
#p.mat 
# Correlation matrix
M <- cor(na.omit(match[,c(14:35)]), method = "spearman")
#M

### Plot correlogram
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot
setwd(WD)
pdf("corrplot_v9v5.1.3.2_28_11_18.pdf", height = 15, width = 15)

corrplot(M, method = "color", col = col(200),  
         		type = "upper", order = "hclust", 
         	   	addCoef.col = "black", # Ajout du coefficient de corrélation
         	  	tl.col = "black", tl.srt = 45, #Rotation des etiquettes de textes
         	 	# Combiner avec le niveau de significativité
         		#p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         		# Cacher les coefficients de corrélation sur la diagonale
         		diag = FALSE )
dev.off()






### 28/08/18: Also plot the distrbution of the matched data to help choose between 2 colinear variables (opt for the one that is most likely to b normally distributed)

### A) v9v3.1
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
d <- get(load(dir()[1]))
names <- colnames(d)[c(1:6,10:16,25:47)]
names
# Get all obs data
matched <- lapply(dir(), function(d) {
				data <- get(load(d))
				message(paste(d, sep = ""))
				return(subset(data, select = names ))	
}
) # eo lapply
match <- do.call(rbind, matched)
dim(match) # 766 033 obs for 3.1 ; 1 212 771 for 3.2
rm(matched) ; gc()
# Change 2 colnames to match vector of predictors below
colnames(match)[c(16,18)] <- c("Bathy","deltaT")
colnames(match)

# Vector of env predictors
predictors <- colnames(match)[c(14:37)] ; predictors

for(p in predictors) {
	
		message(paste("Plotting ", p, sep = ""))
		
		### If else loop to determine appropriate bin for dist plot
		if(p == "Bathy") {
				bin <- 100
		} else if (p == "deltaT") {
				bin <- 0.5
		} else if (p == "Chl") {
				bin <- 0.1
		} else if (p == "logChl") {
				bin <- 0.1
		} else if (p == "MLD1") {
				bin <- 10
		} else if (p == "logMLD") {
				bin <- 0.1
		} else if (p == "SLA") {
				bin <- 0.01
		} else if (p == "Wind") {
				bin <- 1
		} else if (p == "SST") {
				bin <- 1
		} else if (p == "SSS") {
				bin <- 1
		} else if (p == "PAR") {
				bin <- 1
		} else if (p == "MLPAR") {
				bin <- 1	
		} else if (p == "SiO2") {
				bin <- 1
		} else if (p == "Sistar") {
				bin <- 1
		} else if (p == "Nstar") {
				bin <- 1
		} else if (p == "NO3") {
				bin <- 1
		} else {
				bin <- 0.1
		} # eo if else loop
		
		# plot
		plot <- ggplot(match, aes(x = match[,p])) +
		   	 	geom_histogram(binwidth = bin, colour = "black", fill = "grey60") +
		    	geom_vline(aes(xintercept = median(match[,p], na.rm = T)), color = "red", linetype = "dashed", size = 1) +
				theme_bw() + xlab(paste(p, sep = "")) 
		
		### Save plot
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v9/")
		ggsave(plot = plot, filename = paste("plot_distrib_matched_",p,"_v9v3.1.pdf", sep = ""), dpi = 300, width = 6, height = 5) 
		
	
} # eo for loop




### B) v9v3.2
setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.2/")
d <- get(load(dir()[1]))
names <- colnames(d)[c(1:6,10:16,25:47)]
names
# Get all obs data
matched <- lapply(dir(), function(d) {
				data <- get(load(d))
				message(paste(d, sep = ""))
				return(subset(data, select = names ))	
}
) # eo lapply
match <- do.call(rbind, matched)
dim(match) # 766 033 obs for 3.1 ; 1 212 771 for 3.2
rm(matched) ; gc()
# Change 2 colnames to match vector of predictors below
colnames(match)[c(16,18)] <- c("Bathy","deltaT")
colnames(match)

# Vector of env predictors
predictors <- colnames(match)[c(16:36)] ; predictors

for(p in predictors) {
	
		message(paste("Plotting ", p, sep = ""))
		
		### If else loop to determine appropriate bin for dist plot
		if(p == "Bathy") {
				bin <- 100
		} else if (p == "deltaT") {
				bin <- 0.5
		} else if (p == "Chl") {
				bin <- 0.1
		} else if (p == "logChl") {
				bin <- 0.1
		} else if (p == "MLD1") {
				bin <- 10
		} else if (p == "logMLD") {
				bin <- 0.1
		} else if (p == "SLA") {
				bin <- 0.01
		} else if (p == "Wind") {
				bin <- 1
		} else if (p == "SST") {
				bin <- 1
		} else if (p == "SSS") {
				bin <- 1
		} else if (p == "PAR") {
				bin <- 1
		} else if (p == "MLPAR") {
				bin <- 1	
		} else if (p == "SiO2") {
				bin <- 1
		} else if (p == "Sistar") {
				bin <- 1
		} else if (p == "Nstar") {
				bin <- 1
		} else if (p == "NO3") {
				bin <- 1
		} else {
				bin <- 0.1
		} # eo if else loop
		
		# plot
		plot <- ggplot(match, aes(x = match[,p])) +
		   	 	geom_histogram(binwidth = bin, colour = "black", fill = "grey60") +
		    	geom_vline(aes(xintercept = median(match[,p], na.rm = T)), color = "red", linetype = "dashed", size = 1) +
				theme_bw() + xlab(paste(p, sep = "")) 
		
		### Save plot
		setwd("/UP_home/fabioben/Desktop/OVERSEE/data/biology/occurence_data_groups/v9/")
		ggsave(plot = plot, filename = paste("plot_distrib_matched_",p,"_v9v3.2.pdf", sep = ""), dpi = 300, width = 6, height = 5)
	
} # eo for loop


# ----------------------------------------------

### 11/03/2019: Match each P/A data with the logEKE climatologies (for phyto- and zooplankton), plot distribution and 
### heatmap of Spearman's correlations coefficients

### 1°) Get the monththly clims of logEKE
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
stack <- get(load("EKE_stack_1d.Rdata"))
class(stack)
stack

### 2°) get phyto data
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
files <- dir()[grep("data_total",dir())]
require("parallel")
res <- mclapply(files, function(f) {
			message(paste("Reading ", f, sep = ""))
			d <- read.table(f, sep = ";", h = T)
			return(d)
	}, mc.cores = 20
) # eo mclapply
phyto.data <- do.call(rbind, res)
rm(res)
dim(phyto.data)
head(phyto.data)

phyto.data$EKE <- NA
phyto.data[phyto.data$month == 1,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jan'), 
												y = phyto.data[which(phyto.data$month == 1),c("x","y")], method = 'bilinear')
												
phyto.data[phyto.data$month == 2,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Feb'), 
												y = phyto.data[which(phyto.data$month == 2),c("x","y")], method = 'bilinear')

phyto.data[phyto.data$month == 3,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Mar'), 
												y = phyto.data[which(phyto.data$month == 3),c("x","y")], method = 'bilinear')

phyto.data[phyto.data$month == 4,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Apr'), 
												y = phyto.data[which(phyto.data$month == 4),c("x","y")], method = 'bilinear')
												
phyto.data[phyto.data$month == 5,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_May'), 
												y = phyto.data[which(phyto.data$month == 5),c("x","y")], method = 'bilinear')

phyto.data[phyto.data$month == 6,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jun'), 
												y = phyto.data[which(phyto.data$month == 6),c("x","y")], method = 'bilinear')

phyto.data[phyto.data$month == 7,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jul'), 
												y = phyto.data[which(phyto.data$month == 7),c("x","y")], method = 'bilinear')
												
phyto.data[phyto.data$month == 8,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Aug'), 
												y = phyto.data[which(phyto.data$month == 8),c("x","y")], method = 'bilinear')

phyto.data[phyto.data$month == 9,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Sep'), 
												y = phyto.data[which(phyto.data$month == 9),c("x","y")], method = 'bilinear')
												
phyto.data[phyto.data$month == 10,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Oct'), 
												y = phyto.data[which(phyto.data$month == 10),c("x","y")], method = 'bilinear')
												
phyto.data[phyto.data$month == 11,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Nov'), 
												y = phyto.data[which(phyto.data$month == 11),c("x","y")], method = 'bilinear')

phyto.data[phyto.data$month == 12,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Dec'), 
												y = phyto.data[which(phyto.data$month == 12),c("x","y")], method = 'bilinear')

# And log transform
summary(phyto.data$EKE)
summary(log(phyto.data$EKE))
summary(log1p(phyto.data$EKE))

phyto.data$log1pEKE <- log1p(phyto.data$EKE)
phyto.data$logEKE <- log(phyto.data$EKE)
setwd(WD)

# Check
plot <- ggplot(phyto.data, aes(x = logEKE)) +
   	 	geom_histogram(binwidth = .1, colour = "black", fill = "grey60") +
    	geom_vline(aes(xintercept = median(phyto.data[,"logEKE"], na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("logEKE", sep = "")) 
ggsave(plot = plot, filename = paste("plot_distrib_matched_phyto_","logEKE","_v9v3.1.pdf", sep = ""), dpi = 300, width = 6, height = 5)

plot <- ggplot(phyto.data, aes(x = log1pEKE)) +
   	 	geom_histogram(binwidth = .01, colour = "black", fill = "grey60") +
    	geom_vline(aes(xintercept = median(phyto.data[,"log1pEKE"], na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("log1pEKE", sep = "")) 
ggsave(plot = plot, filename = paste("plot_distrib_matched_phyto_","log1PEKE","_v9v3.1.pdf", sep = ""), dpi = 300, width = 6, height = 5)


### Keep logEKE & plot correlogram
library("corrplot")
colnames(phyto.data)
p.mat <- cor.mtest(na.omit(phyto.data[phyto.data$obs == 1,c(6,7,9,10:12,15:22,35:37,41:43,50,79,81)]))
p.mat 
# Correlation matrix
M <- round(cor(na.omit(phyto.data[phyto.data$obs == 1,c(6,7,9,10:12,15:22,35:37,41:43,50,79,81)]), method = "spearman"),3)
M

### Plot correlogram
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot
setwd(WD)
pdf("corrplot_phyto_presences_only_11_03_19.pdf", height = 17, width = 17)

corrplot(M, method = "color", col = col(200), type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = F)
dev.off()

### Match each species dataset with logEKE (like above) and save as .txt
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/group_background/species_data/")
files <- dir()[grep("data_",dir())]

for(f in files) {
	
		message(paste("Reading ", f, sep = ""))
		d <- read.table(f, sep = ";", h = T)
		
		# Match logEKE data
		message(paste("Matching EKE values for ", unique(d[1,"species"]), sep = ""))
		d$EKE <- NA
		d[d$month == 1,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jan'), 
								y = d[which(d$month == 1),c("x","y")], method = 'bilinear')												
		d[d$month == 2,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Feb'), 
								y = d[which(d$month == 2),c("x","y")], method = 'bilinear')
		d[d$month == 3,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Mar'), 
								y = d[which(d$month == 3),c("x","y")], method = 'bilinear')
		d[d$month == 4,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Apr'), 
								y = d[which(d$month == 4),c("x","y")], method = 'bilinear')												
		d[d$month == 5,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_May'), 
								y = d[which(d$month == 5),c("x","y")], method = 'bilinear')
		d[d$month == 6,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jun'), 
								y = d[which(d$month == 6),c("x","y")], method = 'bilinear')
		d[d$month == 7,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jul'), 
								y = d[which(d$month == 7),c("x","y")], method = 'bilinear')												
		d[d$month == 8,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Aug'), 
								y = d[which(d$month == 8),c("x","y")], method = 'bilinear')
		d[d$month == 9,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Sep'), 
								y = d[which(d$month == 9),c("x","y")], method = 'bilinear')											
		d[d$month == 10,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Oct'), 
								y = d[which(d$month == 10),c("x","y")], method = 'bilinear')												
		d[d$month == 11,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Nov'), 
								y = d[which(d$month == 11),c("x","y")], method = 'bilinear')
		d[d$month == 12,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Dec'), 
								y = d[which(d$month == 12),c("x","y")], method = 'bilinear')

		# And log transform
		d$logEKE <- log(d$EKE)
		
		# Save
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/group_background/")
		write.table(d, file = f, sep = ";")
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/group_background/species_data/")
	
} # eo for loop

# Check some files
d <- read.table("data_total_Tripos_arietinus.txt", sep = ";", h = T)


### 3°) Same but for zoo data
# v9v8v5v3.1
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
d <- get(load(dir()[1]))
names <- colnames(d)[c(1:8,10:16,27:48)]
names
# Get all obs data
matched <- lapply(dir(), function(d) {
				data <- get(load(d))
				message(paste(d, sep = ""))
				return(subset(data, select = names ))	
}
) # eo lapply
match <- do.call(rbind, matched)
dim(match) # 764 159 obs for 3.1 
rm(matched) ; gc()

### Match with logEKE data
match$EKE <- NA
match[match$month == 1,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jan'), 
												y = match[which(match$month == 1),c("x","y")], method = 'bilinear')												
match[match$month == 2,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Feb'), 
												y = match[which(match$month == 2),c("x","y")], method = 'bilinear')
match[match$month == 3,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Mar'), 
												y = match[which(match$month == 3),c("x","y")], method = 'bilinear')
match[match$month == 4,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Apr'), 
												y = match[which(match$month == 4),c("x","y")], method = 'bilinear')												
match[match$month == 5,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_May'), 
												y = match[which(match$month == 5),c("x","y")], method = 'bilinear')
match[match$month == 6,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jun'), 
												y = match[which(match$month == 6),c("x","y")], method = 'bilinear')
match[match$month == 7,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jul'), 
												y = match[which(match$month == 7),c("x","y")], method = 'bilinear')												
match[match$month == 8,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Aug'), 
												y = match[which(match$month == 8),c("x","y")], method = 'bilinear')
match[match$month == 9,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Sep'), 
												y = match[which(match$month == 9),c("x","y")], method = 'bilinear')												
match[match$month == 10,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Oct'), 
												y = match[which(match$month == 10),c("x","y")], method = 'bilinear')												
match[match$month == 11,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Nov'), 
												y = match[which(match$month == 11),c("x","y")], method = 'bilinear')
match[match$month == 12,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Dec'), 
												y = match[which(match$month == 12),c("x","y")], method = 'bilinear')
# And log transform
summary(match$EKE)
summary(log(match$EKE))

#match$log1pEKE <- log1p(match$EKE)
match$logEKE <- log(match$EKE)
setwd(WD)

# Check
plot <- ggplot(match, aes(x = logEKE)) +
   	 	geom_histogram(binwidth = .1, colour = "black", fill = "grey60") +
    	geom_vline(aes(xintercept = median(match[,"logEKE"], na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("logEKE", sep = "")) 
ggsave(plot = plot, filename = paste("plot_distrib_matched_zoo_","logEKE","_v9v3.1.pdf", sep = ""), dpi = 300, width = 6, height = 5)


### Plot correlogram
library("corrplot")
colnames(match)
p.mat <- cor.mtest(na.omit(match[,c(16:39)]))
p.mat 
# Correlation matrix
M <- round(cor(na.omit(match[,c(16:39)]), method = "spearman"),3)
M
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot
setwd(WD)
pdf("corrplot_zoo_presences_v9v5.1.3.1_11_03_19.pdf", height = 15, width = 15)
corrplot(M, method = "color", col = col(200),type = "upper", order = "hclust", 
         	   	addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE )
dev.off()

### Ok, so match each species dataset with logEKE (like above) and save as .txt
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/group_bckgrnd")
files <- dir()[grep("data",dir())]

for(f in files) {
	
		message(paste("Reading ", f, sep = ""))
		d <- read.table(f, sep = ";", h = T)
		
		# Match logEKE data
		message(paste("Matching EKE values for ", unique(d[1,"species"]), sep = ""))
		d$EKE <- NA
		d[d$month == 1,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jan'),
                                y = d[which(d$month == 1),c("x","y")], method = 'bilinear')												
		d[d$month == 2,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Feb'), 
								y = d[which(d$month == 2),c("x","y")], method = 'bilinear')
		d[d$month == 3,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Mar'), 
								y = d[which(d$month == 3),c("x","y")], method = 'bilinear')
		d[d$month == 4,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Apr'), 
								y = d[which(d$month == 4),c("x","y")], method = 'bilinear')												
		d[d$month == 5,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_May'), 
								y = d[which(d$month == 5),c("x","y")], method = 'bilinear')
		d[d$month == 6,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jun'), 
								y = d[which(d$month == 6),c("x","y")], method = 'bilinear')
		d[d$month == 7,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Jul'), 
								y = d[which(d$month == 7),c("x","y")], method = 'bilinear')												
		d[d$month == 8,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Aug'), 
								y = d[which(d$month == 8),c("x","y")], method = 'bilinear')
		d[d$month == 9,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Sep'), 
								y = d[which(d$month == 9),c("x","y")], method = 'bilinear')											
		d[d$month == 10,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Oct'), 
								y = d[which(d$month == 10),c("x","y")], method = 'bilinear')												
		d[d$month == 11,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Nov'), 
								y = d[which(d$month == 11),c("x","y")], method = 'bilinear')
		d[d$month == 12,"EKE"] <- raster::extract(x = subset(x = stack, subset = 'EKE_Dec'), 
								y = d[which(d$month == 12),c("x","y")], method = 'bilinear')

		# And log transform
		d$logEKE <- log(d$EKE)
		
		# Save
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/group_bckgrnd")
		write.table(d, file = f, sep = ";")
	
} # eo for loop














