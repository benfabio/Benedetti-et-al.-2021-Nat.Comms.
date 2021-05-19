
##### 14/02/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#		- loading the netCDF files of altimetry data downloaded from AVISO/Copernicus
#		- examine content, structure and variables (absolute geostrophic velocities, anomaie,s sea level hight etc.)
#		- try to compute (for one .nc file) the mean eddy kinetic energy from the zonal and meridional components (U & V)
#		- do so for each nc file and compute the 12 monthly climatologies, 1°x1° resolution, of EKE or MEKE. 
 
### Last update : 21/02/2019

# ----------------------------------------------------------------------------------------

library("ncdf4")
library("raster")
library("rgdal")
library("sp")
library("tidyverse")
library("stringr")
library("reshape2")
library("viridis")

# ----------------------------------------------------------------------------------------

### Coastline fir later plotting
cl <- read.csv("world_coast.csv", h = TRUE)
coast <- list(
 	   # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	 	geom_polygon(aes(x = lon, y = lat), data = cl, fill = "grey55"),
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
  	  	theme(panel.background = element_rect(fill = "white"),  # background
    		legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70")
  		 )
)  # eo list

# Go to dir with the daily altimetry data
setwd(paste(getwd(),"/","AVISO_altimetry_data_14_02_19","/", sep = ""))
files <- dir()[grep(".nc", dir())]
# Get first nc file for testing
f <- files[1]

# Open it and check it out
nc <- nc_open(f)
nc
names(nc$var)
#  [1] "crs"      "lat_bnds" "lon_bnds" "ugos"     "vgos"     "err"     
#  [7] "vgosa"    "adt"      "sla"      "ugosa" 

# ugos = surface_geostrophic_eastward_sea_water_velocity
# vgos = surface_geostrophic_northward_sea_water_velocity
# vgosa = surface_geostrophic_northward_sea_water_velocity_assuming_sea_level_for_geoid ; 
# ugosa = surface_geostrophic_eastward_sea_water_velocity_assuming_sea_level_for_geoid
# The geostrophic velocity anomalies are referenced to the [1993, 2012] period

# sla = The sea level anomaly is the sea surface height above mean sea surface; it is referenced to the [1993, 2012] period

# Check one variable out, sla for instance
sla <- ncvar_get(nc, "sla") 
class(sla)
dim(sla)
head(sla)
# Can you rasterize it directly? 
ras <- raster::stack(f, varname = "ugos")
ras # ok...
ras2 <- raster::stack(f, varname = "vgos")
ras2 # ok.
sla <- raster::stack(f, varname = "sla")
sla <- as.data.frame(sla, xy = T)
dim(sla)
summary(sla)

rm(sla,ras,ras2)

### According to: Qiu & Chen (2004)
# https://journals.ametsoc.org/doi/pdf/10.1175/1520-0485(2004)034%3C1515%3ASMITEF%3E2.0.CO%3B2
# Eddy Kinetic Energy (EKE) has to be computed from the 2 components of sea surface height (SSH) anomalies as follows: 
# EKE = (1/2)*[(ugosa)^2 + (vgosa)^2]

# So, cbind ugosa and vgosa and compute EKE for file f
# Might to rotate too

ugosa <- as.data.frame(rotate(raster::stack(f, varname = "ugosa")), xy = T)
colnames(ugosa)[3] <- "ugosa"
vgosa <- as.data.frame(rotate(raster::stack(f, varname = "vgosa")), xy = T)
colnames(vgosa)[3] <- "vgosa"
summary(ugosa)
summary(vgosa)
# cbind, 
table <- cbind(ugosa, vgosa[,"vgosa"])
colnames(table)[4] <- "vgosa"

# Compute EKE:
table$eke <- 0.5*( (table$ugosa)^2 + (table$vgosa)^2 )
summary(table)

# Map it
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(eke)), data = na.omit(table)) + coast + 
		scale_fill_distiller(name = "EKE (m2/s2)", palette = "RdYlBu", limits = c(-10,1)) + 
		coord_quickmap() + theme_bw()
ggsave(plot = map, filename = "test_map.pdf", dpi = 300, height = 6, width = 8)

# Map it
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(eke)), data = na.omit(table)) + coast + 
		scale_fill_viridis(name = "EKE (m2/s2)", limits = c(-10,1)) + coord_quickmap() + theme_bw()
ggsave(plot = map2, filename = "test_map2.pdf", dpi = 300, height = 6, width = 8)


### Amazing ! 

# Now, the same but by degrading the resolution to 1°x1°
head(round(table$x)) ; head(round(table$y))
table$xbin <- round(table$x)
table$ybin <- round(table$y)
table$id <- paste(table$xbin, table$ybin, sep = "_")
# Compute evarge EKE using dplyr
require("dplyr")
clim <- data.frame(table %>%
			group_by(id) %>%
			summarise(x = unique(xbin), y = unique(ybin), mean = mean(eke) )
) # eo ddf
summary(clim)
str(clim)

# Map it
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(mean)), data = na.omit(clim)) + coast + 
		scale_fill_distiller(name = "Mean EKE", palette = "RdYlBu", limits = c(-10,.5))  + 
		coord_quickmap() + theme_bw()
ggsave(plot = map, filename = "test_map_1d.pdf", dpi = 300, height = 6, width = 8)

# Map it
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(mean)), data = na.omit(clim)) + coast + 
		scale_fill_viridis(name = "Mean EKE", limits = c(-10,.5)) + coord_quickmap() + theme_bw()
ggsave(plot = map2, filename = "test_map2_1d.pdf", dpi = 300, height = 6, width = 8)


### Do the same, but with ugos and vgos

ugos <- as.data.frame(rotate(raster::stack(f, varname = "ugos")), xy = T)
colnames(ugos)[3] <- "ugos"
vgos <- as.data.frame(rotate(raster::stack(f, varname = "vgos")), xy = T)
colnames(vgos)[3] <- "vgos"
summary(ugos)
summary(vgos)
# cbind, 
table <- cbind(ugos, vgos[,"vgos"])
colnames(table)[4] <- "vgos"

# Compute EKE:
table$eke <- 0.5*( (table$ugos)^2 + (table$vgos)^2 )
summary(table)

# Map it
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(eke)), data = na.omit(table)) + coast + 
		scale_fill_distiller(name = "MKE (m2/s2)", palette = "RdYlBu", limits = c(-11,1)) + 
		coord_quickmap() + theme_bw()
ggsave(plot = map, filename = "test_map_mke.pdf", dpi = 300, height = 6, width = 8)

# Map it
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(eke)), data = na.omit(table)) + coast + 
		scale_fill_viridis(name = "MKE (m2/s2)", limits = c(-11,1)) + coord_quickmap() + theme_bw()
ggsave(plot = map2, filename = "test_map2_mke.pdf", dpi = 300, height = 6, width = 8)

### Patterns are clearer with the EKE, but let's compute the climatologies for both anyways.


# ----------------------------------------------------------------------------------------

### 15/02/2019: Compute the 12 monthly climatologies of MKE and EKE.
### Use the file names in 'files' to identify the data associated to each months, and perform the same code as above in a for loop 
### to compute the 12 monthly climatologies of MKE and EKE.

f
unlist(strsplit(f, "_"))[6] # this is how you get the date, now extract the month

str_match_all(unlist(strsplit(f, "_"))[6], ".{2}") # noice! month is the third element
m <- unlist(str_match_all(unlist(strsplit(f, "_"))[6], ".{2}"))[3]
m

# Using lapply, associate a month to each file name 
res <- lapply(files, function(f) {
			m <- unlist(str_match_all(unlist(strsplit(f, "_"))[6], ".{2}"))[3]
			return(data.frame(file = f, month = m))
		} # eo fun
) # eo lapply
table <- do.call(rbind, res)
rm(res)
# head(table) ; dim(table) # Great.

# For m in unique(table$month), get corresponding files and compute climatology from them
# m <- "08"
for(m in unique(table$month) ) {
	
		message(paste("Doing month ", m, sep = ""))
		# Get corresponding files
		filenames <- as.character(table[table$month == m,"file"])
		
		# Extract MKE and EKE data from each of these files
		require("parallel")
		res <- lapply(filenames, function(f) {
					
					# Extract ugos, vgos, ugosa and vgosa
					require("raster")
					ugos <- as.data.frame(rotate(raster::stack(x = f, varname = "ugos")), xy = T)
					colnames(ugos)[3] <- "ugos"
					vgos <- as.data.frame(rotate(raster::stack(x = f, varname = "vgos")), xy = T)
					colnames(vgos)[3] <- "vgos"
					ugosa <- as.data.frame(rotate(raster::stack(x = f, varname = "ugosa")), xy = T)
					colnames(ugosa)[3] <- "ugosa"
					vgosa <- as.data.frame(rotate(raster::stack(x = f, varname = "vgosa")), xy = T)
					colnames(vgosa)[3] <- "vgosa"
					# cbind
					d <- cbind(ugos, vgos[,"vgos"], ugosa[,"ugosa"], vgosa[,"vgosa"])
					colnames(d)[c(4:6)] <- c("vgos","ugosa","vgosa")
					rm(ugos, vgos, vgosa, ugosa)
					# summary(d)
					# Compute MKE + EKE
					d$mke <- 0.5*( (d$ugos)^2 + (d$vgos)^2 )
					d$eke <- 0.5*( (d$ugosa)^2 + (d$vgosa)^2 )
					# Return the mke & eke (with coords)
					return( d[,c(1,2,7,8)] )
				}
				
		) # eo lapply
		# rbind
		ddf <- do.call(rbind, res)
		dim(ddf); head(ddf)
		rm(res)
		
		# Add cell id and use dplyr to compute climatology
		ddf$id <- as.factor(paste(ddf$x, ddf$y, sep = "_"))
		# length(unique(ddf$id))
		require("dplyr")
		clim <- data.frame(ddf %>%
					group_by(id) %>%
					summarise(x = unique(x), y = unique(y), mean_mke = mean(mke, na.rm = T), mean_eke = mean(eke, na.rm = T) )
		) # eo ddf
		# str(clim)
		# summary(clim)
		
		# Map eke and save climatology
		#map <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(mean_eke)), data = na.omit(clim)) + coast + 
				#scale_fill_viridis(name = "Mean EKE\nlog(m2/s2)", limits = c(-10,.5)) + coord_quickmap() + theme_bw()
		#ggsave(plot = map, filename = paste("map_clim_eke_",m,".pdf", sep = ""), dpi = 300, height = 6, width = 8)
		
		# Save as .txt
		message(paste("Saving climatology for month ", m, sep = ""))
		write.table(clim, file = paste("clim_eke_",m,".txt", sep = ""), sep = ";")
		
		# And now as raster stack
		clim <- clim[,c(2:5)]
		coordinates(clim) <- ~ x + y
		gridded(clim) <- TRUE
		# coerce to raster
		raster <- raster(clim)
		save(raster, file = paste("raster_clim_eke_",m,".Rdata", sep = "") )
		
		# Clean stuff
		rm(raster, clim, map, ddf)
		gc()
		
}


### Make a raster stack out of the 12 rasters
bands <- c("01","02","03","04","05","06","07","08","09","10","11","12")
# b <- "04"
res <- lapply(bands, function(b) {
			ras <- get(load(paste("raster_clim_eke_",b,".Rdata", sep ="")))
			ddf <- as.data.frame(ras, xy = T)
			colnames(ddf) <- c("x","y","EKE")
			if(b == "01") {
				return(ddf)
			} else {
				return(data.frame(EKE = ddf[,"EKE"]))
			} # eo if else loop
	
		} # eo fun
) # eo lapply
table <- do.call(cbind, res)
dim(table) ; colnames(table)
colnames(table)[c(3:14)] <- c("EKE_Jan","EKE_Feb","EKE_Mar","EKE_Apr","EKE_May","EKE_Jun",
							"EKE_Jul","EKE_Aug","EKE_Sep","EKE_Oct","EKE_Nov","EKE_Dec")

# Save as .txt
write.table(table, file = "EKE_stack_full_res.txt", sep = ";")

# And now as raster stack
coordinates(table) <- ~ x + y
gridded(table) <- TRUE
# coerce to raster
stack <- stack(table)
# Save
save(stack, file = "EKE_stack_full_res.Rdata")



### Same, but by degrading the res to 1d and compute average in the 1° cells with dplyr
res <- lapply(bands, function(b) {
			ras <- get(load(paste("raster_clim_eke_",b,".Rdata", sep ="")))
			ddf <- as.data.frame(ras, xy = T)
			colnames(ddf) <- c("x","y","EKE")
			# Compute mean EKE in 1°x1° cells
			ddf$xbin <- round(ddf$x)
			ddf$ybin <- round(ddf$y)
			ddf$id <- paste(ddf$xbin, ddf$ybin, sep = "_")
			# Compute evarge EKE using dplyr
			require("dplyr")
			clim <- data.frame(ddf %>%
						group_by(id) %>%
						summarise(x = unique(xbin), y = unique(ybin), mean = mean(EKE, na.rm = T) )
			) # eo ddf
			# summary(clim)
			if(b == "01") {
				return(clim)
			} else {
				return(data.frame(EKE = clim[,"mean"]))
			} # eo if else loop
	
		} # eo fun
) # eo lapply
table <- do.call(cbind, res)
dim(table) ; colnames(table)
colnames(table)[c(4:15)] <- c("EKE_Jan","EKE_Feb","EKE_Mar","EKE_Apr","EKE_May","EKE_Jun",
							"EKE_Jul","EKE_Aug","EKE_Sep","EKE_Oct","EKE_Nov","EKE_Dec")

# Save as .txt
write.table(table, file = "EKE_stack_1d.txt", sep = ";")
summary(table)

# And now as raster stack
table <- table[,c(2:15)]
coordinates(table) <- ~ x + y
gridded(table) <- TRUE
# coerce to raster
stack <- stack(table)
crs(stack) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
# Save
save(stack, file = "EKE_stack_1d.Rdata")
# summary(as.data.frame(stack, xy = T))


### 20/02/2019: Ok, but still need to re-interpolate on the WOA13 s grid to be consistent 
# Load a raster with the WOA grid
setwd(paste(getwd(),"/","global_monthly_clims_1d","/", sep = ""))
ras <- get(load("dO2_stack_1d.Rdata"))
# Check coords
ras
head(as.data.frame(ras, xy = T)) # ok

xy <- as.data.frame(ras, xy = T)[,c("x","y")]

# Load the newly defined stack of mean EKE
setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/env_predictors","/","AVISO_altimetry_data_14_02_19","/", sep = ""))
eke <- get(load("EKE_stack_1d.Rdata"))
eke
names(eke)

# Use bilinear extraction to get the EKE values on the WOA13's cell grid
raster::extract(x = subset(x = eke, subset = 'EKE_Jan'), y = xy, method = 'bilinear')

# Do it for each 12 months
ddf <- data.frame(
			x = xy$x, 
			y = xy$y, 
			EKE_Jan = raster::extract(x = subset(x = eke, subset = 'EKE_Jan'), y = xy, method = 'bilinear'),
			EKE_Feb = raster::extract(x = subset(x = eke, subset = 'EKE_Feb'), y = xy, method = 'bilinear'),
			EKE_Mar = raster::extract(x = subset(x = eke, subset = 'EKE_Mar'), y = xy, method = 'bilinear'),
			EKE_Apr = raster::extract(x = subset(x = eke, subset = 'EKE_Apr'), y = xy, method = 'bilinear'),
			EKE_May = raster::extract(x = subset(x = eke, subset = 'EKE_May'), y = xy, method = 'bilinear'),
			EKE_Jun = raster::extract(x = subset(x = eke, subset = 'EKE_Jun'), y = xy, method = 'bilinear'),
			EKE_Jul = raster::extract(x = subset(x = eke, subset = 'EKE_Jul'), y = xy, method = 'bilinear'),
			EKE_Aug = raster::extract(x = subset(x = eke, subset = 'EKE_Aug'), y = xy, method = 'bilinear'),
			EKE_Sep = raster::extract(x = subset(x = eke, subset = 'EKE_Sep'), y = xy, method = 'bilinear'),
			EKE_Oct = raster::extract(x = subset(x = eke, subset = 'EKE_Oct'), y = xy, method = 'bilinear'),
			EKE_Nov = raster::extract(x = subset(x = eke, subset = 'EKE_Nov'), y = xy, method = 'bilinear'),
			EKE_Dec = raster::extract(x = subset(x = eke, subset = 'EKE_Dec'), y = xy, method = 'bilinear')
	
) # eo ddf

head(ddf)
summary(ddf)

# Save as table and then as stack
write.table(ddf, file = "EKE_stack_1d.txt", sep = ";")
coordinates(ddf) <- ~ x + y
gridded(ddf) <- TRUE
# coerce to raster
stack <- stack(ddf)
crs(stack) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
stack
# head(as.data.frame(stack, xy = T))
# Save
save(stack, file = "EKE_stack_1d.Rdata")


### 21/02/2019: cbind 1d EKE monthly clims to each glob_stack_month -> glob_stack_month_21_02_19.txt
eke <- get(load("EKE_stack_1d.Rdata"))
eke
eke <- as.data.frame(eke, xy = TRUE)
dim(eke)
head(eke)
summary(log1p(eke$EKE_Oct))

# Plots of different logs to assess which has best distribution
# Non logged
ggplot(eke, aes(x = eke[,3])) +
		geom_histogram(binwidth = .01, colour = "black", fill = "grey60") +
		geom_vline(aes(xintercept = median(eke[,3], na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("Mean EKE", sep = "")) 

# log
quartz()
ggplot(eke, aes(x = log(eke[,3]) )) +
		geom_histogram(binwidth = .1, colour = "black", fill = "grey60") +
		geom_vline(aes(xintercept = median(log(eke[,3]), na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("Mean EKE - natural log", sep = "")) 


# log1p
quartz()
ggplot(eke, aes(x = log1p(eke[,3]) )) +
		geom_histogram(binwidth = .001, colour = "black", fill = "grey60") +
		geom_vline(aes(xintercept = median(log1p(eke[,3]), na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("Mean EKE - log1p", sep = "")) 

# log10
quartz()
ggplot(eke, aes(x = log10(eke[,3]))) +
		geom_histogram(binwidth = .1, colour = "black", fill = "grey60") +
		geom_vline(aes(xintercept = median(log10(eke[,3]), na.rm = T)), color = "red", linetype = "dashed", size = 1) +
		theme_bw() + xlab(paste("Mean EKE - log10", sep = "")) 

### Pick either the natural log or the power log


### Load each of the 12 monthly climatologies 
clim <- read.table("glob_stack_month_jan_31_01_19.txt", h = TRUE, sep = ";")
#dim(clim)
#head(clim)
clim$EKE <- eke$EKE_Jan
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_jan_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_feb_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Feb
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_feb_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_mar_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Mar
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_mar_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_apr_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Apr
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_apr_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_may_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_May
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_may_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_jun_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Jun
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_jun_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_jul_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Jul
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_jul_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_aug_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Aug
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_aug_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_sep_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Sep
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_sep_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_oct_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Oct
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_oct_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_nov_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Nov
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_nov_21_02_19.txt", sep = ";")

clim <- read.table("glob_stack_month_dec_31_01_19.txt", h = TRUE, sep = ";")
clim$EKE <- eke$EKE_Dec
clim$logEKE <- log(clim$EKE)
write.table(clim, file = "glob_stack_month_dec_21_02_19.txt", sep = ";")

### Pick a few climatologies and examine global correlations between variables to check if mean EKE (logged and not) covary strongly with the other variables
library("corrplot")
# Compute p-val matrix
colnames(clim)
p.mat <- cor.mtest( na.omit(clim[,c(3:4,7:29)]) )
M <- cor(na.omit(clim[,c(3:4,7:29)]), method = "spearman")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Go to plots directory
pdf(paste("heatmap_rich_","Dec","_12_02_19.pdf", sep = ""), height = 17, width = 17)
corrplot(M, method = "color", col = col(200), type = "upper", order = "hclust", 
         	   	addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE )
dev.off()





