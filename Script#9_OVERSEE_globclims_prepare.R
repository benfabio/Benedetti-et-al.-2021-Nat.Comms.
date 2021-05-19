
##### 23/08/18: Script to examine the environmental data (global monthly climatologies, 1°x1°) provided by Damiano and the ones you used during your PhD (also provided by Damiano), and to prepare the rasters before fitting them to the species observations.

library("ncdf4")
library("maptools")
library("sp")
library("raster")
library("rgeos")
library("fields")
library("reshape2")
library("dplyr")

WD <- getwd()

# ==============================================================================================================================

# dir()

##### 1°) Bathy
setwd(paste(WD,"/","Bathymetry NOAA ETOPO1","/", sep = ""))
bathy <- get(load("ETOPO_1by1.RData"))
class(bathy) # OK
plot(bathy) # OK
head(as.data.frame(bathy, xy = T))

##### 2°) SST
setwd(paste(WD,"/","WOA13 Temperature","/","Decadal average","/", sep = ""))
dir()
ras1 <- raster("woa13_decav_t01_01.nc")
class(ras1)
ras1
#plot(ras1) # OK...need to make a stack out of the 12 monthly climatologies
ras2 <- raster("woa13_decav_t02_01.nc")
ras3 <- raster("woa13_decav_t03_01.nc")
ras4 <- raster("woa13_decav_t04_01.nc")
ras5 <- raster("woa13_decav_t05_01.nc")
ras6 <- raster("woa13_decav_t06_01.nc")
ras7 <- raster("woa13_decav_t07_01.nc")
ras8 <- raster("woa13_decav_t08_01.nc")
ras9 <- raster("woa13_decav_t09_01.nc")
ras10 <- raster("woa13_decav_t10_01.nc")
ras11 <- raster("woa13_decav_t11_01.nc")
ras12 <- raster("woa13_decav_t12_01.nc")
#plot(ras12)
# Create stack and rename properly
sst_stak <- stack(ras1, ras2, ras3, ras4, ras5, ras6, ras7, ras8, ras9, ras10, ras11, ras12)
names(sst_stak) <- c("SST_Jan","SST_Feb","SST_Mar","SST_Apr","SST_May","SST_Jun","SST_Jul","SST_Aug","SST_Sep","SST_Oct","SST_Nov","SST_Dec")
plot(sst_stak)
head(as.data.frame(sst_stak, xy = T))
# Save
save(sst_stak, file = "SST_stack_1d.Rdata")
write.table(as.data.frame(sst_stak, xy = T), file = "SST_stack_1d.txt", sep = ";")

##### 2°) Do the same for SSS
setwd(paste(WD,"/","WOA13 Salinity","/","Decadal average", sep = ""))
dir()
ras1 <- raster("woa13_decav_s01_01.nc")
ras2 <- raster("woa13_decav_s02_01.nc")
ras3 <- raster("woa13_decav_s03_01.nc")
ras4 <- raster("woa13_decav_s04_01.nc")
ras5 <- raster("woa13_decav_s05_01.nc")
ras6 <- raster("woa13_decav_s06_01.nc")
ras7 <- raster("woa13_decav_s07_01.nc")
ras8 <- raster("woa13_decav_s08_01.nc")
ras9 <- raster("woa13_decav_s09_01.nc")
ras10 <- raster("woa13_decav_s10_01.nc")
ras11 <- raster("woa13_decav_s11_01.nc")
ras12 <- raster("woa13_decav_s12_01.nc")
# Create stack and rename properly
sss_stak <- stack(ras1, ras2, ras3, ras4, ras5, ras6, ras7, ras8, ras9, ras10, ras11, ras12)
names(sss_stak) <- c("SSS_Jan","SSS_Feb","SSS_Mar","SSS_Apr","SSS_May","SSS_Jun","SSS_Jul","SSS_Aug","SSS_Sep","SSS_Oct","SSS_Nov","SSS_Dec")
plot(sss_stak)
# Save
save(sss_stak, file = "SSS_stack_1d.Rdata")
write.table(as.data.frame(sss_stak, xy = T), file = "SSS_stack_1d.txt", sep = ";")



##### 3°) Do the same for Silicates, Nitrates and Phosphates 
setwd(paste(WD,"/","WOA13 Phosphate", sep = ""))
dir()
ras1 <- raster("woa13_all_p01_01.nc")
ras2 <- raster("woa13_all_p02_01.nc")
ras3 <- raster("woa13_all_p03_01.nc")
ras4 <- raster("woa13_all_p04_01.nc")
ras5 <- raster("woa13_all_p05_01.nc")
ras6 <- raster("woa13_all_p06_01.nc")
ras7 <- raster("woa13_all_p07_01.nc")
ras8 <- raster("woa13_all_p08_01.nc")
ras9 <- raster("woa13_all_p09_01.nc")
ras10 <- raster("woa13_all_p10_01.nc")
ras11 <- raster("woa13_all_p11_01.nc")
ras12 <- raster("woa13_all_p12_01.nc")
# Create stack and rename properly
PO4_stak <- stack(ras1, ras2, ras3, ras4, ras5, ras6, ras7, ras8, ras9, ras10, ras11, ras12)
names(PO4_stak) <- c("PO4_Jan","PO4_Feb","PO4_Mar","PO4_Apr","PO4_May","PO4_Jun","PO4_Jul","PO4_Aug","PO4_Sep","PO4_Oct","PO4_Nov","PO4_Dec")
plot(PO4_stak)
# Save
PO4_stak
save(PO4_stak, file = "PO4_stack_1d.Rdata")
write.table(as.data.frame(PO4_stak, xy = T), file = "PO4_stack_1d.txt", sep = ";")



##### 4°) For SeaWifs' PAR
setwd(paste(WD,"/","SeaWiFS PAR and CHL","/", sep = ""))
ras <- stack("SeaWiFS_PAR_MO_climatology_1deg.nc")
ras
names(ras) <- c("PAR_Jan","PAR_Feb","PAR_Mar","PAR_Apr","PAR_May","PAR_Jun","PAR_Jul","PAR_Aug","PAR_Sep","PAR_Oct","PAR_Nov","PAR_Dec")
plot(ras)
# Save
save(ras, file = "PAR_stack_1d.Rdata")
write.table(as.data.frame(ras, xy = T), file = "PAR_stack_1d.txt", sep = ";")


##### 5°) For SeaWifs' chlorophyll a
setwd(paste(WD,"/","SeaWiFS PAR and CHL","/", sep = ""))
ras <- stack("SeaWiFS_CHL_MO_climatology_1deg.nc")
ras
names(ras) <- c("Chl_Jan","Chl_Feb","Chl_Mar","Chl_Apr","Chl_May","Chl_Jun","Chl_Jul","Chl_Aug","Chl_Sep","Chl_Oct","Chl_Nov","Chl_Dec")
plot(ras)
# Save
save(ras, file = "Chl_stack_1d.Rdata")
write.table(as.data.frame(ras, xy = T), file = "Chl_stack_1d.txt", sep = ";")


##### 6°) For DBM's MLD
setwd(paste(WD,"/","MLD De Boyer Montegut","/","Output", sep = ""))
ras <- get(load("deBoyerMontegut_mld_RasterStack_V1.RData"))
ras
names(ras) <- c("MLD_Jan","MLD_Feb","MLD_Mar","MLD_Apr","MLD_May","MLD_Jun","MLD_Jul","MLD_Aug","MLD_Sep","MLD_Oct","MLD_Nov","MLD_Dec")
plot(ras)
# Save
save(ras, file = "MLD2_stack_1d.Rdata")
write.table(as.data.frame(ras, xy = T), file = "MLD1_stack_1d.txt", sep = ";")


##### 7°) For MLPAR & ∆SST: go to subfolder
setwd(paste(WD,"/","For MLPAR and dT","/", sep = ""))
# Get multistacks
d <- get(load("multistack_1.Rdata"))
class(d) # ddf
colnames(d) # 
head(d)

# Compare spatial resolution to the one of the other clims
setwd("/Users/fabiobenedetti/Desktop/PostDocs/ETHZ/OVERSEE/data/")
#sst_stak <- get(load("SST_stack_1d.Rdata"))
#head(as.data.frame(sst_stak, xy = TRUE)) # OK, it's the same. Let's extract monthly MLPAR and ∆SST

months <- c(1:12)
setwd(paste(WD,"/","For MLPAR and dT","/", sep = ""))
stacks <- lapply(months, function(m) {
		s <- get(load(paste("multistack_",m,".RData", sep = "")))
		return(s[,c(3:32)])
} # eo FUN
) # eo lapply

# bind
ddf <- do.call(cbind, stacks)
dim(ddf)
colnames(ddf)
# From ddf, retrieve the coordinates xy and the monthly clims of MLPAR and deltaT, then convert them to rasterStacks

# First, MLPAR:
mlpar_stack <- ddf[,grep("MLPAR1",colnames(ddf))]
mlpar_stack$x <- ddf$x
mlpar_stack$y <- ddf$y
summary(mlpar_stack)

colnames(mlpar_stack)[1:12] <- c("MLPAR_Jan","MLPAR_Feb","MLPAR_Mar","MLPAR_Apr","MLPAR_May","MLPAR_Jun","MLPAR_Jul","MLPAR_Aug",
								"MLPAR_Sep","MLPAR_Oct","MLPAR_Nov","MLPAR_Dec")

# Convert to stack
coordinates(mlpar_stack) <- ~ x + y
gridded(mlpar_stack) <- TRUE
# coerce to raster
raster <- stack(mlpar_stack)
raster
plot(raster)
# Save
save(raster, file = "MLPAR_stack_1d.Rdata")

# Second, deltaT
# summary(ddf[,grep("deltaT",colnames(ddf))]) # !!! NOT monthly
dT_stack <- data.frame(ddf[,c(22)])
dT_stack$x <- ddf$x
dT_stack$y <- ddf$y
summary(dT_stack)

colnames(dT_stack)[1] <- c("deltaT")

# Convert to stack
coordinates(dT_stack) <- ~ x + y
gridded(dT_stack) <- TRUE
# coerce to raster
raster <- stack(dT_stack)
raster
plot(raster)
# Save
save(raster, file = "deltaT_stack_1d.Rdata")


# ==============================================================================================================================

### 24/08/2018: Now, use the clims you created above to derive other climatologies of potentially important predictors: 
# v logChla
# v logMLD
# v logNO3
# v logPO4
# v logSiO2
# v Nstar
# v Pstar
# v Sistar

dir()[grep(".Rdata",dir())]

##### logChla
#chla <- get(load("Chl_stack_1d.Rdata"))
chla <- read.table("Chl_stack_1d.txt", h = T, sep = ";")
logchla <- cbind(chla[,c(1,2)], log10(chla[,c(3:14)]))
summary(logchla)
colnames(logchla)[3:14] <- c("logChl_Jan","logChl_Feb","logChl_Mar","logChl_Apr","logChl_May","logChl_Jun","logChl_Jul",
				"logChl_Aug","logChl_Sep","logChl_Oct","logChl_Nov","logChl_Dec")
				
write.table(logchla, "logChl_stack_1d.txt", sep = ";")


# Save
save(logchla, file = "logChl_stack_1d.Rdata")
# logchla <- get(load("logChl_stack_1d.Rdata"))
# logchla


##### logMLD
mld <- get(load("MLD1_stack_1d.Rdata"))
logmld <- log10(mld)
logmld ; plot(logmld)
names(logmld) <- c("logMLD_Jan","logMLD_Feb","logMLD_Mar","logMLD_Apr","logMLD_May","logMLD_Jun","logMLD_Jul",
				"logMLD_Aug","logMLD_Sep","logMLD_Oct","logMLD_Nov","logMLD_Dec")

# Save
save(logmld, file = "logMLD_stack_1d.Rdata")


##### 27/08/18 : need to use log1p instead of log10 for nutrients
##### Nutrients (NO3/PO4/SiO2)
#SiO2 <- get(load("SiO2_stack_1d.Rdata"))
PO4 <- read.table("PO4_stack_1d.txt", h = T, sep = ";")
colnames(PO4)
summary(PO4)
logPO4 <- cbind(PO4[,c(1,2)], log1p(PO4[,c(3:14)]))
summary(logPO4)
colnames(logPO4)[c(3:14)] <- c("logPO4_Jan","logPO4_Feb","logPO4_Mar","logPO4_Apr","logPO4_May","logPO4_Jun","logPO4_Jul",
				"logPO4_Aug","logPO4_Sep","logPO4_Oct","logPO4_Nov","logPO4_Dec")

write.table(logPO4, "logPO4_stack_1d.txt", sep = ";")

#logSiO2 <- log10(SiO2)
names(logSiO2) <- c("logSiO2_Jan","logSiO2_Feb","logSiO2_Mar","logSiO2_Apr","logSiO2_May","logSiO2_Jun","logSiO2_Jul",
				"logSiO2_Aug","logSiO2_Sep","logSiO2_Oct","logSiO2_Nov","logSiO2_Dec")
				
logSiO2 ; plot(logSiO2)

# Save
save(logSiO2, file = "logSiO2_stack_1d.Rdata")


#### Nstar
no3 <- get(load("NO3_stack_1d.Rdata"))
po4 <- get(load("PO4_stack_1d.Rdata"))
nstar <- no3 - (16*(po4))
names(nstar) <- c("Nstar_Jan","Nstar_Feb","Nstar_Mar","Nstar_Apr","Nstar_May","Nstar_Jun","Nstar_Jul",
				"Nstar_Aug","Nstar_Sep","Nstar_Oct","Nstar_Nov","Nstar_Dec")
				
nstar ; plot(nstar)
# Save
save(nstar, file = "Nstar_stack_1d.Rdata")


#### Pstar
pstar <- po4 - ((1/16)*no3)
names(pstar) <- c("Pstar_Jan","Pstar_Feb","Pstar_Mar","Pstar_Apr","Pstar_May","Pstar_Jun","Pstar_Jul",
				"Pstar_Aug","Pstar_Sep","Pstar_Oct","Pstar_Nov","Pstar_Dec")
				
pstar ; plot(pstar)
# Save
save(pstar, file = "Pstar_stack_1d.Rdata")

#### Sistar
sio2 <- get(load("SiO2_stack_1d.Rdata"))
sistar <- sio2 - no3
names(sistar) <- c("Sistar_Jan","Sistar_Feb","Sistar_Mar","Sistar_Apr","Sistar_May","Sistar_Jun","Sistar_Jul",
				"Sistar_Aug","Sistar_Sep","Sistar_Oct","Sistar_Nov","Sistar_Dec")
				
sistar ; plot(sistar)
# Save
save(sistar, file = "Sistar_stack_1d.Rdata")


# ==============================================================================================================================

### Check out the other drivers of zooplankton distribution you may find interesting...

setwd(paste(WD,"/","For MLPAR and dT","/", sep = ""))
months <- c(1:12)
stacks <- lapply(months, function(m) {
		s <- get(load(paste("multistack_",m,".RData", sep = "")))
		return(s[,c(3:32)])
} # eo FUN
) # eo lapply

# bind
ddf <- do.call(cbind, stacks)
dim(ddf)
colnames(ddf)

### In addition to the predictors chosen by Damiano, you may try to use : 
# - Wind speed
# - SLA


sla_stack <- ddf[,grep("SLA",colnames(ddf))]
sla_stack$x <- ddf$x
sla_stack$y <- ddf$y
summary(sla_stack)

colnames(sla_stack)[1:12] <- c("SLA_Jan","SLA_Feb","SLA_Mar","SLA_Apr","SLA_May","SLA_Jun","SLA_Jul","SLA_Aug",
								"SLA_Sep","SLA_Oct","SLA_Nov","SLA_Dec")

# Convert to stack
coordinates(sla_stack) <- ~ x + y
gridded(sla_stack) <- TRUE
# coerce to raster
raster <- stack(sla_stack)
raster
plot(raster)
# Save
save(raster, file = "SLA_stack_1d.Rdata")


# ==============================================================================================================================

### Now, map the each monthly climatology (one panel figure per variable, 12*22 = 264 plots)
### In the process, make sure all climatologies follow the same spatial cell grid of 1°x1°

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

### Get world coastline
cl <- read.csv("world_coast.csv", h = T)
coast <- list(
  # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  geom_polygon(aes(x = lon, y = lat), data = cl, fill = "grey50"),
  geom_path(aes(x = lon, y = lat), data = cl, colour = "black", linetype = 1, size = 0.5),
  # appropriate projection
  coord_quickmap(),
  # remove extra space around the coast
  scale_x_continuous(name = "Longitude", 
                     breaks = c(-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180), 
                     labels = c("-180°E","-150°E","-120°E","-90°E","-60°E","-30°E","0°","30°E","60°E","90°E","120°E","150°E","180°E"), 
                     expand = c(0,0)), 
  
  scale_y_continuous(name = "Latitude", 
                     breaks = c(-90,-60,-30,0,30,60,90), 
                     labels = c("-90°N","-60°N","-30°N","0°","30°N","60°N","90°N"),
                     expand = c(0,0)),
  # dark gray background for the panel and legend
  theme(
    panel.background = element_rect(fill = "white"),  # background
    legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "white")
  )
)

# ggplot() + coast + coord_quickmap()

### For each rasterStack :
# - load as stack
# - convert to data.frame
# - melt so you have months as factors
# - map with facet_grid using months as factors

WD <- getwd()
setwd(paste(WD,"/","global_monthly_clims_1d","/", sep = ""))
clims <- dir()[grep("23_10_18",dir())]
clims # ok
# c <- "glob_stack_month_Jul_23_10_18.txt"
for(c in clims) {
	
		message(paste("Mapping ", c, sep = ""))
		setwd(paste(WD,"/","global_monthly_clims_1d","/", sep = ""))
		### Load clim
		if( grepl(".Rdata",c) ) {
			cl <- get(load(c))
		} else {
			cl <- read.table(c, sep = ";", h = TRUE)
		} # eo first if else loop
		
		### Convert to ddf if cl is of class raster
		if( class(cl) != "data.frame" ) {
			cl_ddf <- as.data.frame(cl, xy = TRUE)
			# Melt
			cl_ddf_m <- melt(cl_ddf, id.vars = c("x","y"))
			colnames(cl_ddf_m) <- c("x","y","var","value")
		} else {
			cl_ddf_m <- melt(cl, id.vars = c("x","y"))
			colnames(cl_ddf_m) <- c("x","y","var","value")
		} # eo second if else loop
		
		### Map in panel
		if( grepl("SSS",c) ) {
				plot <- ggplot() + geom_raster(aes(x = x, y = y, fill = value), data = cl_ddf_m[cl_ddf_m$value > 30,]) + coast + 
					scale_fill_viridis() + facet_wrap(~ var) 
				
		} else {
				plot <- ggplot(data = cl_ddf_m) + geom_raster(aes(x = x, y = y, fill = value) ) + coast + 
						scale_fill_viridis() + facet_grid(~ var, scales = "free") 
		} # eo third if else loop		 
		
		### Save plot
		setwd(WD)
		ggsave(plot = plot, filename = paste("map_",c,".pdf", sep = ""), dpi = 300, width = 25, height = 15)
	
} # eo for loop

### 27/11/18: After correcting July and August log(Chl) climatologies: for individual re-mapping 
setwd(paste(WD,"/","global_monthly_clims_1d","/", sep = ""))
clim_jul <- read.table("glob_stack_month_Jul_23_10_18.txt", h = T, sep = ";")
clim_aug <- read.table("glob_stack_month_Aug_23_10_18.txt", h = T, sep = ";")
head(clim_jul)
head(clim_aug)

# Map logChl
plot1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = logChl), data = clim_jul) + coast + 
			scale_fill_viridis() + theme_bw() + coord_quickmap()
#
plot2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = logChl), data = clim_aug) + coast + 
			scale_fill_viridis() + theme_bw() + coord_quickmap()
#
ggsave(plot = plot1, filename = paste("map_logChl_Jul.pdf", sep = ""), dpi = 300, width = 10, height = 8)
ggsave(plot = plot2, filename = paste("map_logChl_Aug.pdf", sep = ""), dpi = 300, width = 10, height = 8)



### 24/08/18: Same as above, but plot distribution instead
for(c in clims) {
	
		message(paste("Plotting ", c, sep = ""))
		setwd(paste(WD,"/","global_monthly_clims","/", sep = ""))
		### Load clim
		if( grepl(".Rdata",c) ) {
			cl <- get(load(c))
		} else {
			cl <- read.table(c, sep = ";", h = TRUE)
		} # eo first if else loop
		
		### Convert to ddf if cl is of class raster
		if( class(cl) != "data.frame" ) {
			cl_ddf <- as.data.frame(cl, xy = TRUE)
			# Melt
			cl_ddf_m <- melt(cl_ddf, id.vars = c("x","y"))
			colnames(cl_ddf_m) <- c("x","y","month","value")
		} else {
			cl_ddf_m <- melt(cl, id.vars = c("x","y"))
			colnames(cl_ddf_m) <- c("x","y","month","value")
		} # eo second if else loop
		
		plot <- ggplot(cl_ddf_m, aes(x = value)) +
		   	 	geom_histogram(binwidth = .1 , colour = "black", fill = "white") +
		    	geom_vline(aes(xintercept = median(value, na.rm = T)), color = "red", linetype = "dashed", size = 1) +
				facet_wrap(~ month) 
		
		### Save plot
		setwd(WD)
		ggsave(plot = plot, filename = paste("plot_distrib_",c,".pdf", sep = ""), dpi = 300, width = 15, height = 10)
	
} # eo for loop


### And check if they all follow the same grid
for(c in clims) {
	
		message(paste("Checking ", c, sep = ""))
		setwd(paste(WD,"/","global_monthly_clims","/", sep = ""))
		### Load clim
		if( grepl(".Rdata",c) ) {
			cl <- get(load(c))
		} else {
			cl <- read.table(c, sep = ";", h = TRUE)
		} # eo first if else loop
		
		### Convert to ddf if cl is of class raster
		if( class(cl) != "data.frame" ) {
			cl_ddf <- as.data.frame(cl, xy = TRUE)
			# Melt
			cl_ddf_m <- melt(cl_ddf, id.vars = c("x","y"))
			colnames(cl_ddf_m) <- c("x","y","month","value")
		} else {
			cl_ddf_m <- melt(cl, id.vars = c("x","y"))
			colnames(cl_ddf_m) <- c("x","y","month","value")
		} # eo second if else loop
		
		message(head(cl_ddf_m[,c("x","y")]))
		message(paste(" ", sep = ""))
		message(paste(" ", sep = ""))
		
} # eo for loop

### OK, need to slightly re-interpolate the following climatologies: 
# SLA
# Wind
# deltaT
# MLPAR

#?interpolate

setwd(paste(WD,"/","global_monthly_clims","/", sep = ""))

# Get raster to interpolate: 
ras2interp <- get(load("Wind_stack_1d.Rdata"))
crs(ras2interp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
head(as.data.frame(ras2interp, xy = TRUE))
summary(as.data.frame(ras2interp, xy = TRUE))

data <- as.data.frame(ras2interp, xy = TRUE)
head(data)
data$y


# ==============================================================================================================================

##### 25/08/18: For SOME reason, KRYO won't accept .Rdata from my computer but only .txt files...convert the 14 climatologies that 
##### are not .Rdata to .txt files and upload on kryo...

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

clims <- dir()[grep(".Rdata",dir())]  # select the 14 Rdata files
# Get rid of MLD2...
clims <- clims[c(1:8,10:15)]

# For testing
c <- "Bathy_ETOPO_1d.Rdata"
c <- "MLD1_stack_1d.Rdata"

# For each clim, load and turn into ddf
for(c in clims) { 
		cl <- get(load(c))
		cl_ddf <- as.data.frame(cl, xy = TRUE)
		# colnames(cl_ddf) 
		write.table(cl_ddf, file = paste(str_replace_all(c,".Rdata",".txt"), sep = ""), sep = ";")
} # eo for loop


# ==============================================================================================================================

##### 29/08/18: For projecting the ENMs, create monthly stacks of env predictors (1 stack per month of the year; containing all 21 predictors) 
### Beware that deltaT & bathymetry are not monthly resolved
files <- dir()[grep("stack_1d.txt", dir())] ; files
months <- c(3:14) # keeping the indices of the columns in the files (1:2 correspond to x & y) so it is easier to retrive the columns
f <- files[20]
m <- 4

for(m in months) {
	
		# Useless message
		message(paste("Making stack for month = ", m, sep = ""))
	
		# Get all the predictors' values for month 'm'
		preds <- lapply(files, function(f) {
					# Load the file corresponding to the pred
					message(paste("Loading ", f, sep = ""))
					p <- read.table(f, h = T, sep = ";")
					# Get the name of the pred
					name <- str_replace_all(f, "_stack_1d.txt", "")
					# Choose the columns to return (month m & coordinates if p = bathymetry)
					if( f == "Bathy_stack_1d.txt" ) {
							ddf <- p[,c(1,2,3)]
							colnames(ddf)[1:3] <- c("x","y",name)
							return(ddf)
					} else if ( f == "deltaT_stack_1d.txt" ) {
							ddf <- p[,c(1,2,3)]
							colnames(ddf)[1:3] <- c("x","y",name)
							return(ddf)
					} else {	
							ddf <- data.frame(p[,c(m)])
							colnames(ddf)[1] <- name
							return(ddf)
					} # eo if else loop
		}
		) # eo lapply
		# cbind
		stck <- do.call(cbind, preds)
		# dim(stck) ; str(stck) ; summary(stck)
		rm(preds)
		
		# Save monthly stck
		write.table(stck, paste("glob_stack_month_",m,".txt", sep = ""), sep = ";")
		rm(stck) ; gc()
		message(paste("                 ", sep = ""))
	
} # eo for loop 


### Check the dimensions of each monthly stack
files <- dir()[grep("glob_stack_month", dir())] ; files
for(f in files) {	
		stck <- read.table(f, h = T, sep = ";")
		message(paste(dim(stck), sep = ""))
} # eo for loop


### 31/01/2019: Add dO2 (@175m depth) to the 12 monthly climatology txt files ("glob_monthly_clims")

library("sp")
library("raster")
library("fields")
library("reshape2")
library("tidyverse")

# Get file ame sof each of the 12 clim
clims <- dir()[grep("glob_stack_month_", dir())]
clims

# Get dO2 stack
dO2 <- read.table("dO2_stack_1d.txt", h = T, sep = ";")
dim(dO2)
head(dO2)
summary(dO2)

# c <- "glob_stack_month_apr_23_10_18.txt"
for(c in clims) {
	
	clim <- read.table(c, sep = ";", h = T)
	
	if(c == "glob_stack_month_apr_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Apr"]
	} else if (c == "glob_stack_month_aug_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Aug"]
	} else if (c == "glob_stack_month_dec_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Dec"]
	} else if (c == "glob_stack_month_feb_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Feb"]
	} else if (c == "glob_stack_month_jan_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Jan"]
	} else if (c == "glob_stack_month_jul_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Jul"]
	} else if (c == "glob_stack_month_jun_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Jun"]
	} else if (c == "glob_stack_month_mar_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Mar"]
	} else if (c == "glob_stack_month_may_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_May"]
	} else if (c == "glob_stack_month_nov_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Nov"]
	} else if (c == "glob_stack_month_oct_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Oct"]
	} else if (c == "glob_stack_month_sep_23_10_18.txt") {
			clim$dO2 <- dO2[,"dO2_Sep"]
	} # else if loop
	
	# Save by changing name (date update)
	write.table(x = clim, file = str_replace_all(c, "23_10_18", "31_01_19"), sep = ";")
	
}








