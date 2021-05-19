
##### 18/10/2018: R Script to prepare monthly climatologies of pCO2, match them with a P/A dataset (v9v3.2 for instance) and 
#####			  assess collinearity with other env predictors and perform a variable importance test (forward LM and RF)

### Aims to
#	- load the nc file from SOCAT (https://www.socat.info/index.php/data-access/)
#	- or the one from Dr. Peter Landschützer: http://cdiac.ess-dive.lbl.gov/ftp/oceans/spco2_1998_2011_ETH_SOM-FFN/
#	- convert to monthly rasters with the same grid as WOA13 - maps
#	- load the v9v3.2 datasets (total background)
#	- use the extract() function of 'raster' to match the species obs to pCO2 values
#	- examine distribution of the matched pCO2 data
#	- assess colinearity with other predictors
#	- perform variable importance using forward LM and RF, also add MLPAR
#	- add pCO2 data to raster stacks

### Latest update: 23/10/2018

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("ggplot2")
library("RColorBrewer")
library("viridis")

# ==============================================================================================================================

# You may need this
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/")
### World coastline:
cl <- read.csv("world_coast.csv", h = TRUE)
coast <- list(
 	   # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	 	geom_polygon(aes(x = lon, y = lat), data = cl, fill= "grey50"),
  	  	geom_path(aes(x = lon, y = lat), data = cl, colour = "black", linetype = 1),
  	  	# appropriate projection
  	  	coord_quickmap(),
  	  	# remove extra space around the coast
  	  	scale_x_continuous(name = "Longitude", 
                     breaks = c(-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180), 
                     labels = c("-180°E","-150°E","-120°E","-90°E","-60°E","-30°E","0°E","30°E","60°E","90°E","120°E","150°E","180°E"), 
                     expand = c(0,0)), 
  
  		scale_y_continuous(name = "Latitude", 
                     breaks = c(-90,-60,-30,0,30,60,90), 
                     labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"),
                     expand = c(0,0)),
  		# dark gray background for the panel and legend
  	  	theme(
    		panel.background = element_rect(fill = "white"),  # background
    		legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70")
  		 )
) 


##### 1°) Load the netCDF from SOCAT, extract climatologies etc.
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/")
dir()

### A) Latest SOCAT downloaded at https://www.socat.info/index.php/data-access/
nc <- nc_open("SOCATv6_tracks_gridded_monthly.nc")
nc
names(nc$var)
#  [1] "tmnth_bnds"            "count_ncruise"         "fco2_count_nobs"
#  [4] "fco2_ave_weighted"     "fco2_ave_unwtd"        "fco2_min_unwtd"
#  [7] "fco2_max_unwtd"        "fco2_std_weighted"     "fco2_std_unwtd"
# [10] "sst_count_nobs"        "sst_ave_weighted"      "sst_ave_unwtd"
# [13] "sst_min_unwtd"         "sst_max_unwtd"         "sst_std_weighted"
# [16] "sst_std_unwtd"         "salinity_count_nobs"   "salinity_ave_weighted"
# [19] "salinity_ave_unwtd"    "salinity_min_unwtd"    "salinity_max_unwtd"
# [22] "salinity_std_weighted" "salinity_std_unwtd"    "lat_offset_unwtd"
# [25] "lon_offset_unwtd"
#nc_close(nc)

#         xlon  Size:360
#             units: degrees_east
#             point_spacing: even
#             axis: X
#             modulo: 360
#             standard_name: longitude
#         ylat  Size:180
#             units: degrees_north
#             point_spacing: even
#             axis: Y
#             standard_name: latitude
#         tmnth  Size:576   *** is unlimited ***
#             units: days since 1970-01-01 00:00:00
#             axis: T
#             bounds: tmnth_bnds
#             time_origin: 1-JAN-1970
#             standard_name: time

lon <-  ncvar_get(nc, "lon_offset_unwtd") 
lat <-  ncvar_get(nc, "lat_offset_unwtd") 
#class(lon) # arrays
#str(lon)
lon <- melt(lon)
lat <- melt(lat)
summary(lon)
summary(lat)

time <- ncvar_get(nc, "tmnth_bnds") 
time <- melt(time)
summary(time)

p <- ncvar_get(nc, "fco2_ave_weighted") 
# Mean of fco2 recomputed computed by calculating the arithmetic mean value 
# for each cruise passing through the cell and then averaging these datasets.
p <- melt(p)
summary(p)
colnames(p)[1:4] <- c("x","y","time","pCO2")
# Transform lon and lat
p$x <- p$x - 180
p$y <- p$y - 90
# Convert time (days since 1970-01-01 00:00:00)
unique(p$time)
#as.Date(166:407, origin = '1970-01-01')
#summary (as.Date(p$time, origin = '1970-01-01') )
#summary( as.Date(as.POSIXct(p$time, origin = "1970-01-01")) )
### Weeeeiiirrddd...

#p$time <- as.Date(p$time, origin = '1970-01-01')
# # Restrict to one time value and map to see what's behind
# pp <- p[p$time == "1970-03-19",]
# dim(pp)
# summary(pp)
# map <- ggplot() + geom_raster(aes(x = x, y = y, fill = pCO2), data = pp) + coast +
# 			scale_fill_viridis(name = "pCO2 (µatm)") + coord_quickmap() + theme_bw()
# ggsave(plot = map, filename = "map_pCO2_1970-03-19.pdf", dpi = 300, width = 10, height = 7)

### Input from Cara: Each time slice corresponds to all the data that have been collected in a month
unique(p$time)
# For instance: 1970-01-02 -> January 1970 ; 1970-01-03 --> Feb 1970 ; 1970-01-04 --> Mar 1970
length(unique(p$time)) # 576
length(unique(p$time)) / 12 # 48
# 1970+58 --> 2018 !! 
rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), times = 48)
length(rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), times = 48)) # 576

time_values <- data.frame(ref = unique(p$time), months = rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), times = 48))
head(time_values)

# Provide months to p
p$month <- NA
for(i in unique(p$time)) {
		m <- unique(time_values[time_values$ref == i,"months"])
		p[p$time == i,"month"] <- m
} # eo for loop
head(p)
str(p)

### Use dplyr to get monthly data (12 time slices then)
p$cell_id <- paste(p$x, p$y, sep = "_")
table <- p %>% group_by(cell_id, month) %>% summarise(pCO2 = mean(pCO2, na.rm = T), x = unique(x), y = unique(y) )
table <- data.frame(table)
summary(table)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = pCO2), data = table[table$month == 10,]) + coast +
 			scale_fill_viridis(name = "pCO2 (µatm) - October") + coord_quickmap() + theme_bw()
ggsave(plot = map, filename = "map_pCO2_SOCAT_Oct.pdf", dpi = 300, width = 10, height = 7)





# ==========================================================

### B) pCO2 Landschützer et al. Data - data file downloaded from net/kryo/work/up/updata/ or  http://cdiac.ornl.gov/ftp/oceans/spco2_1998_2011_ETH_SOM-FFN
### Same file used by Damiano Roghetti
ncid <- nc_open("spco2_clim_ETH_SOM-FFN_CDIAC_G05.nc")
class(ncid)
names(ncid$var)
lon <-  ncvar_get(ncid,"lon") # not existing
lat <-  ncvar_get(ncid,"lat") # not existing
# month <- ncvar_get(ncid,"time") #  not existing
# para <- ncvar_get(ncid, "spco2_clim") # ARRAY, 360,180, 12
nc_close(ncid)

# Load ncdf file as RasterStack (PROBLEM FOR LATE USE: Source code is maintained)
uu <- stack("spco2_clim_ETH_SOM-FFN_CDIAC_G05.nc") 
uu
table <- as.data.frame(uu, xy = T)
colnames(table)[3:14] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
summary(table)
head(table)

# map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = pCO2_Nov), data = table) + coast +
# 			scale_fill_viridis(name = "pCO2 (µatm)") + coord_quickmap() + theme_bw()
# #
# ggsave(plot = map1, filename = "map_pCO2_Nov.pdf", dpi = 300, width = 10, height = 7)

mtbl <- melt(table, id.vars = c("x","y"))
colnames(mtbl) <- c("x","y","month","value")
head(mtbl)
plot <- ggplot() + geom_raster(aes(x = x, y = y, fill = value), data = mtbl) + coast + scale_fill_viridis() + facet_wrap(~ month) 	
ggsave(plot = plot, filename = paste("map_","pCO2",".pdf", sep = ""), dpi = 300, width = 25, height = 15)	
	
### Use uu for the next steps :-)			
names(uu) <- c("pCO2_Jan","pCO2_Feb","pCO2_Mar","pCO2_Apr","pCO2_May","pCO2_Jun","pCO2_Jul","pCO2_Aug","pCO2_Sep","pCO2_Oct","pCO2_Nov","pCO2_Dec")

			
# ==========================================================			
			
### C) Load the v9v3.2 datasets (total background)	
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.2/total_bckgrnd/")			
files <- dir()[grep("data_total_",dir())]
library("parallel")
data <- mclapply(X = files, FUN = function(f) {
				message(paste("Loading ", f, sep = ""))
				d <- read.table(f, h = T, sep = ";")
				return(d)
		}, mc.cores = 20
) # eo mclapply
# rbind
bio <- do.call(rbind, data)	
dim(bio) # 12227200       41
str(bio)			
			
			
### D) Match to uu data

# Create logChl_1d vector
bio$pCO2 <- NA
# Fill it according to month
bio[which(bio$month == 1),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Jan'), 
		y = bio[which(bio$month == 1),c("x","y")], method = 'bilinear')
bio[which(bio$month == 2),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Feb'), 
		y = bio[which(bio$month == 2),c("x","y")], method = 'bilinear')
bio[which(bio$month == 3),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Mar'), 
		y = bio[which(bio$month == 3),c("x","y")], method = 'bilinear')
bio[which(bio$month == 4),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Apr'), 
		y = bio[which(bio$month == 4),c("x","y")], method = 'bilinear')
bio[which(bio$month == 5),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_May'), 
		y = bio[which(bio$month == 5),c("x","y")], method = 'bilinear')
bio[which(bio$month == 6),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Jun'), 
		y = bio[which(bio$month == 6),c("x","y")], method = 'bilinear')
bio[which(bio$month == 7),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Jul'), 
		y = bio[which(bio$month == 7),c("x","y")], method = 'bilinear')
bio[which(bio$month == 8),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Aug'), 
		y = bio[which(bio$month == 8),c("x","y")], method = 'bilinear')
bio[which(bio$month == 9),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Sep'), 
		y = bio[which(bio$month == 9),c("x","y")], method = 'bilinear')
bio[which(bio$month == 10),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Oct'), 
		y = bio[which(bio$month == 10),c("x","y")], method = 'bilinear')
bio[which(bio$month == 11),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Nov'), 
		y = bio[which(bio$month == 11),c("x","y")], method = 'bilinear')
bio[which(bio$month == 12),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Dec'), 
		y = bio[which(bio$month == 12),c("x","y")], method = 'bilinear')

#
summary(bio$pCO2)
colnames(bio)
colnames(bio)[c(19,21)] <- c("Bathy","dSST")


### !! Input from Meike: tets the LOG of bathymetry to suppress the influence of higher values in the midelle of the ocean
bio$logBathy <- log(abs(bio$Bathy))
summary(bio$logBathy)
### For a better map:
require("marmap")
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, res = 15)
bathy <- as.xyz(bathy)
colnames(bathy) <- c("x","y","z")
summary(bathy)
bathy <- bathy[bathy$z < 0,]
bathy$z <- abs(bathy$z)
bathy$logged <- log(bathy$z)
plot <- ggplot() + geom_point(aes(x = x, y = y, colour = logged), data = bathy) + coast + scale_colour_viridis() 
ggsave(plot = plot, filename = "map_bathy_logged.pdf", dpi = 300, width = 10, height = 7)	




### E) Assess spearman's rank correlations with other env data (cols 19 to 39)
cormat <- round(cor(na.omit(bio[,c(19:39,42:43)]), method = "spearman"), 2)
head(cormat)

# Necessary FUNs
get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
	  
reorder_cormat <- function(cormat) {
		# Utiliser la corrélation entre les variables
		# comme mésure de distance
		dd <- as.dist((1-cormat) / 2)
		hc <- hclust(dd)
		cormat <-cormat[hc$order, hc$order]
}	
# Reordonner la matrice de corrélation
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
			theme_minimal() + # minimal theme
			theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

### Okray, let's add some text
heatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
			theme( axis.title.x = element_blank(),
  		  			axis.title.y = element_blank(),
  				  	panel.grid.major = element_blank(),
  				  	panel.border = element_blank(),
 				   	panel.background = element_blank(),
  				  	axis.ticks = element_blank(),
  				  	legend.justification = c(1, 0),
  				  	legend.position = c(0.6, 0.7),
  				  	legend.direction = "horizontal") +
  			guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
#quartz()
#heatmap		
setwd("/net/kryo/work/fabioben/OVERSEE/data/")
ggsave("heatmap.pdf", plot = heatmap, dpi = 300, height = 15, width = 15)

### pCO2 seems to be significantly correlated to log(Chl) but only at the -0.6 level...


### F) OK, now perform the variables' importance test like in script#14.3 !

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.2/total_bckgrnd/")			
files <- dir()[grep("data_total", dir())]
# f <- files[45]

for(f in files) {
	
		# Load the species dataset
		data <- read.table(f, h = T, sep = ";")
		sp <- unique(data[data$obs == 1,"species"])
		phyl <- unique(data[data$obs == 1,"phylum"])
		message(paste("Ranking variables importance for ",sp, " =================", sep = ""))
		message(paste("", sep = ""))
		
		# Create pCO2 vector
		data$logBathy <- log(abs(data$bathymetry_1d))
		
		# Create pCO2 vector
		data$pCO2 <- NA
		# Fill it according to month
		data[which(data$month == 1),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Jan'), 
				y = data[which(data$month == 1),c("x","y")], method = 'bilinear')
		data[which(data$month == 2),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Feb'), 
				y = data[which(data$month == 2),c("x","y")], method = 'bilinear')
		data[which(data$month == 3),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Mar'), 
				y = data[which(data$month == 3),c("x","y")], method = 'bilinear')
		data[which(data$month == 4),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Apr'), 
				y = data[which(data$month == 4),c("x","y")], method = 'bilinear')
		data[which(data$month == 5),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_May'), 
				y = data[which(data$month == 5),c("x","y")], method = 'bilinear')
		data[which(data$month == 6),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Jun'), 
				y = data[which(data$month == 6),c("x","y")], method = 'bilinear')
		data[which(data$month == 7),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Jul'), 
				y = data[which(data$month == 7),c("x","y")], method = 'bilinear')
		data[which(data$month == 8),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Aug'), 
				y = data[which(data$month == 8),c("x","y")], method = 'bilinear')
		data[which(data$month == 9),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Sep'), 
				y = data[which(data$month == 9),c("x","y")], method = 'bilinear')
		data[which(data$month == 10),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Oct'), 
				y = data[which(data$month == 10),c("x","y")], method = 'bilinear')
		data[which(data$month == 11),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Nov'), 
				y = data[which(data$month == 11),c("x","y")], method = 'bilinear')
		data[which(data$month == 12),"pCO2"] <- raster::extract(x = subset(x = uu, subset = 'pCO2_Dec'), 
				y = data[which(data$month == 12),c("x","y")], method = 'bilinear')
		
		### ================================================
		# Strat 1: linear model with forward selection of variables
		d <- data[,c("obs","SST","SSS","logChl","logSiO2","pCO2","Nstar","SLA","MLPAR","logBathy")]
		
		lm <- lm(obs ~ ., data = d, weights = data$weights)
		require("caret")
		ranklm <- varImp(lm)
		ranklm$norm <- ranklm$Overall/ max(ranklm$Overall)
		ranklm$species <- sp
		ranklm$phyl <- phyl
		ranklm$var <- rownames(ranklm)
		colnames(ranklm)[1] <- "rank"
		rm(lm)
		ranklm <- ranklm[order(ranklm$rank, decreasing = T),]
	
		### ================================================
		# Strat 2: fast random forest
		require("ranger")
		rf <- ranger(obs ~ ., data = na.omit(d), num.trees = 1000, min.node.size = 100, importance = "impurity")
		#rf$variable.importance
		rankRF <- data.frame(var = names(rf$variable.importance), rank = rf$variable.importance)
		rankRF$norm <- rankRF$rank/ max(rankRF$rank)
		rankRF$species <- sp
		rankRF$phyl <- phyl
		rm(rf)
		rankRF <- rankRF[order(rankRF$rank, decreasing = T),]
		
		### Go save these tables in a dir
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.2/total_bckgrnd/var_import_tests2/")
		save(ranklm, file = paste("var_import_test1_nodSST_", sp, ".Rdata", sep = "") )
		save(rankRF, file = paste("var_import_test2_nodSST_", sp, ".Rdata", sep = "") )
		rm(rankRF, ranklm, phyl, sp, d)
		setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.2/total_bckgrnd/")
	
} # eo for loop




# ==========================================================

##### 23/10/18: Load Landschüzter's monthly clims
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/")
# ncid <- nc_open("spco2_clim_ETH_SOM-FFN_CDIAC_G05.nc")
# class(ncid)
# names(ncid$var)
# lon <-  ncvar_get(ncid,"lon") # not existing
# lat <-  ncvar_get(ncid,"lat") # not existing
# # month <- ncvar_get(ncid,"time") #  not existing
# # para <- ncvar_get(ncid, "spco2_clim") # ARRAY, 360,180, 12
# nc_close(ncid)

# Load ncdf file as RasterStack (PROBLEM FOR LATE USE: Source code is maintained)
uu <- stack("spco2_clim_ETH_SOM-FFN_CDIAC_G05.nc") 
uu
names(uu) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# Turn into ddf
table <- as.data.frame(uu, xy = T)
colnames(table)
summary(table)
head(table)

# Combine with prexisting glob_stack_month files
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
# Example: april
apr <- read.table("glob_stack_month_Apr.txt", h = T, sep = ";")
colnames(apr)
head(apr)
apr$pCO2 <- table$Apr
# OK, do the same but 12 times
dir()

for(m in names(uu)) {
	
		# Load glob_stack_month
		message(paste("Loading clim for ", m, sep = ""))
		stck <- read.table(paste("glob_stack_month_",m,".txt", sep = ""), h = T, sep = ";")
		stck$pCO2 <- table[,c(m)]
		# Save 
		write.table(stck, file = paste("glob_stack_month_",m,"_23_10_18.txt", sep = ""), sep = ";")
	
} # eo for loop

### Check
d <- read.table("glob_stack_month_Apr_23_10_18.txt", h = T, sep = ";")
colnames(d)
summary(d)





			
			

