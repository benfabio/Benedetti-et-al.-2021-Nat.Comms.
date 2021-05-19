
##### 09/11/2018 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#		- loading and examining the future (2100) projections from the MAREMIP data (all models and the GFDL-TOPAZ model)
#		- extract and format (same grid etc.) the monthly data to perform the quick climate change projections 
#		- for GFDL-TOPAZ: use the MLD from "ALL MODELS" and for the latter use the surface salinity from GFDL-TOPAZ
#		- for PAR (irradiance): just use the same values from the contemporary monthly cliamtologies; shouldn't change in 80yrs
 
### Last update : 12/11/2018

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("stringr")
library("reshape2")
library("tidyverse")
library("viridis")
library("ncdf4")
library("RColorBrewer")
library("oce")
library("viridis")

# --------------------------------------------------------------------------------------------------------------------------------

WD <- getwd()

# Get world coastline
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/")
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
setwd(WD)


# --------------------------------------------------------------------------------------------------------------------------------


##### 1°) Load the data and examine them

### A) Start with ALL MODELS
setwd(paste(WD,"/","ALL_MODELS_Charlotte","/", sep = ""))
nc <- nc_open("SST_all_models_MAREMIP.nc")
print(nc)
names(nc$var) # contains the data from 9 different models 

# Let's start with PISCES (because French)
pisces <- ncvar_get(nc, "DATA_pisces")
class(pisces) # array
dim(pisces) # lon - lat - time
head(pisces)
# melt 
mpisces <- melt(pisces)
dim(mpisces)
head(mpisces) # lon/lat/var
colnames(mpisces) <- c("x","y","time","SST")
summary(mpisces)
# Convert x and y and add id
#mpisces$x <- mpisces$x - 180
mpisces$y <- mpisces$y - 90
mpisces$id <- factor(paste(mpisces$x, mpisces$y, sep = "_"))

# What does the time represent ?
# 2100-2011 --> 89 years. So value = 89 should be SST for 2100
# Let's commpare the sst distribution for 2011 and then for 2100 (should have increased by several degrees) 
#summary(mpisces[mpisces$time == 89,"SST"])
#summary(mpisces[mpisces$time == 1,"SST"]) # average has increased by 3°C ! and max by 4°!

# Compute difference
SST <- dcast(mpisces, id + x + y ~ time, value.var = "SST")
#dim(SST)
#str(SST)
# Change colnames
colnames(SST)[c(4:length(SST))] <- paste("year_",factor(c(2012:2100)), sep = "")
SST$diff <- SST$year_2100 - SST$year_2012
summary(SST$diff)
# map this 

# But first, need to convert it so the 0° lon is centered on the Atl
coordinates(SST) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(SST) <- TRUE
# coerce to raster
stck <- stack(SST)
crs(stck) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
stck <- rotate(stck)

### And re-convert to DDF
SST <- as.data.frame(stck, xy = T)
colnames(SST)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = SST) + coast + 
		scale_fill_gradient2(name = "SST change\n(°C)", low = "blue", high = "red") + 
		coord_quickmap() + xlab("Longitude") + ylab("Latitude") + theme_bw()

ggsave(plot = map, filename = "SST_diff_2100.pdf", dpi = 300, width = 8, height = 6)
# Nice.


### Try another variable (CHL or MLD)
nc <- nc_open("CHL_all_models_MAREMIP.nc")
print(nc)
names(nc$var) # contains the data from 9 different models 
# Take PISCES again
pisces <- ncvar_get(nc, "DATA_pisces")
class(pisces) # array
dim(pisces) # lon - lat - time
head(pisces)
# melt 
mpisces <- melt(pisces)
dim(mpisces)
head(mpisces) # lon/lat/var
colnames(mpisces) <- c("x","y","time","CHL")
summary(mpisces) 
### VERY LOW CHL values because in kg.m^-3
# Convert x and y and add id
#mpisces$x <- mpisces$x - 180
mpisces$y <- mpisces$y - 90
mpisces$id <- factor(paste(mpisces$x, mpisces$y, sep = "_"))
# And convert to mg/m3 -> multiply by 1000000
mpisces$CHL <- (mpisces$CHL)*10^6
summary(mpisces$CHL) 
# ANd log it !
mpisces$CHL <- log10(mpisces$CHL)
summary(mpisces$CHL) 


# DCAST
chl <- dcast(mpisces, id + x + y ~ time, value.var = "CHL")
#dim(SST)
#str(SST)

# Change colnames
colnames(chl)[c(4:length(chl))] <- paste("year_",factor(c(2012:2100)), sep = "")
chl$diff <- chl$year_2100 - chl$year_2012
summary(chl$diff)
# map this 

# But first, need to convert it so the 0° lon is centered on the Atl
coordinates(chl) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(chl) <- TRUE
# coerce to raster
stck <- stack(chl)
crs(stck) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
stck <- rotate(stck)

### And re-convert to DDF
chl <- as.data.frame(stck, xy = T)
colnames(chl)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = chl) + coast + 
		scale_fill_gradient2(name = "[Chl-a] change\n(mg/m3)", low = "orange", high = "green", mid = "white") + 
		coord_quickmap() + xlab("Longitude") + ylab("Latitude") + theme_bw()

ggsave(plot = map, filename = "chl_diff_2100.pdf", dpi = 300, width = 8, height = 6)
# Nice.
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = year_2012), data = chl) + coast + 
		scale_fill_viridis(name = "[Chl-a]\n(mg/m3)" ) + 
		coord_quickmap() + xlab("Longitude") + ylab("Latitude") + theme_bw()

ggsave(plot = map2, filename = "chl_2012.pdf", dpi = 300, width = 8, height = 6)
# Nice.





### And what about MLD ?
nc <- nc_open("MLD_all_models_MAREMIP.nc")
print(nc)
names(nc$var) # contains the data from 9 different models 
# Take PISCES again
pisces <- ncvar_get(nc, "DATA_pisces")
class(pisces) # array
dim(pisces) # lon - lat - time
head(pisces)
# melt 
mpisces <- melt(pisces)
dim(mpisces)
head(mpisces) # lon/lat/var
colnames(mpisces) <- c("x","y","time","MLD")
summary(mpisces) 
### VERY LOW MLD values because in kg.m^-3
# Convert x and y and add id
#mpisces$x <- mpisces$x - 180
mpisces$y <- mpisces$y - 90
mpisces$id <- factor(paste(mpisces$x, mpisces$y, sep = "_"))

# DCAST
mld <- dcast(mpisces, id + x + y ~ time, value.var = "MLD")
# Change colnames
colnames(mld)[c(4:length(mld))] <- paste("year_",factor(c(2012:2100)), sep = "")
mld$diff <- mld$year_2100 - mld$year_2012
summary(mld$diff)
# map this 

# But first, need to convert it so the 0° lon is centered on the Atl
coordinates(mld) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(mld) <- TRUE
# coerce to raster
stck <- stack(mld)
crs(stck) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
stck <- rotate(stck)

### And re-convert to DDF
mld <- as.data.frame(stck, xy = T)
colnames(mld)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = mld) + coast + 
		scale_fill_gradient2(name = "MLD change\n(m)", low = "blue", high = "red", mid = "white") + 
		coord_quickmap() + xlab("Longitude") + ylab("Latitude") + theme_bw()

ggsave(plot = map, filename = "mld_diff_2100.pdf", dpi = 300, width = 8, height = 6)
# Nice.


### And finally, SiO3
nc <- nc_open("SiO3_all_models_MAREMIP.nc")
print(nc)
names(nc$var) # contains the data from 9 different models 
# Take PISCES again
pisces <- ncvar_get(nc, "DATA_pisces")
class(pisces) # array
dim(pisces) # lon - lat - time
head(pisces)
# melt 
mpisces <- melt(pisces)
dim(mpisces)
head(mpisces) # lon/lat/var
colnames(mpisces) <- c("x","y","time","SiO3")
summary(mpisces) 
### VERY LOW SiO3 values because in kg.m^-3
# Convert x and y and add id
#mpisces$x <- mpisces$x - 180
mpisces$y <- mpisces$y - 90
mpisces$id <- factor(paste(mpisces$x, mpisces$y, sep = "_"))

# Silicate are expressed as mol/m3, you need to convert them to µmol/L
mpisces$SiO3 <- (mpisces$SiO3)*10^3 
# And log it
mpisces$SiO3 <- log1p(mpisces$SiO3)
summary(mpisces$SiO3)

# DCAST
sio3 <- dcast(mpisces, id + x + y ~ time, value.var = "SiO3")
# Change colnames
colnames(sio3)[c(4:length(sio3))] <- paste("year_",factor(c(2012:2100)), sep = "")
sio3$diff <- sio3$year_2100 - sio3$year_2012
summary(sio3$diff)
# map this 

# But first, need to convert it so the 0° lon is centered on the Atl
coordinates(sio3) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(sio3) <- TRUE
# coerce to raster
stck <- stack(sio3)
crs(stck) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
stck <- rotate(stck)

### And re-convert to DDF
sio3 <- as.data.frame(stck, xy = T)
colnames(sio3)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff), data = sio3) + coast + 
		scale_fill_gradient2(name = "SiO3 change\n(mumole/L)", low = "blue", high = "red", mid = "white") + 
		coord_quickmap() + xlab("Longitude") + ylab("Latitude") + theme_bw()

ggsave(plot = map, filename = "sio3_diff_2100.pdf", dpi = 300, width = 8, height = 6)
# Nice.

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = year_2012), data = sio3) + coast + 
		scale_fill_viridis(name = "SiO3\nlog(mumole/L)" ) + 
		coord_quickmap() + xlab("Longitude") + ylab("Latitude") + theme_bw()

ggsave(plot = map2, filename = "sio3_2012.pdf", dpi = 300, width = 8, height = 6)
# Nice.



# ---------------------------------------------------------------


### B) And what about GFDL-TOPAZ?
setwd(paste(WD,"/","GFDL_TOPAZ_monthly","/", sep = ""))
dir()


### SST
nc <- nc_open("GFDL_ESM2M_monthly_SST_2100.nc")
print(nc)
names(nc$var) # contains the data from 9 different models 
sst <- ncvar_get(nc, "SST")
class(sst) # array
dim(sst) # lon - lat - months
head(sst)
# melt 
msst <- melt(sst)
dim(msst)
head(msst) # lon/lat/var
colnames(msst) <- c("x","y","month","SST")
summary(msst) 

# Convert lats and provide cell id
msst$y <- msst$y - 90
msst$id <- factor(paste(msst$x, msst$y, sep = "_"))

# DCAST
sst2 <- dcast(msst, id + x + y ~ month, value.var = "SST")
# Change colnames
colnames(sst2)[c(4:length(sst2))] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# But first, need to convert it so the 0° lon is centered on the Atl
coordinates(sst2) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(sst2) <- TRUE
# coerce to raster
stck <- stack(sst2)
crs(stck) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
stck <- rotate(stck)

### And re-convert to DDF
sst2 <- as.data.frame(stck, xy = T)
colnames(sst2)
dim(sst2)

### Need to prepare 2 types of raster stacks: annual and monthly (Jan/Apr/Jul/Oct) from the GFDL-TOPAZ data.
### Strategy: 
#	- for monthly: extract the 4 months' values of SST, SST, CHL, Silicates, compute annual ∆SST in the process
#	- do it for 2 years: 2012 and 2100
#	- cbind the climatological monthly PAR (Jan/Apr/Jul/Oct) from contemporary period (should not change by 2100 since it depends on the Sun)

#	- for annual: use the monthly climatologies from above to compute annual climatologies for 2012 and 2100. 
#	- simply cbind the same annual SST range (max - min) you will have computed for the monthly stacks above
#	- again, compute an annual estimate of PAR and cbind to both present and future stacks
 
### You should obtain: 
#	- 2 times (2012 and 2100) 4 monthly stacks (Jan/Apr/Jul/Oct) of: SST, dSST, SST, PAR, log10(CHL), log1p(SiO3)
#	- 2 times (2012 and 2100) 1 annual stack for the same variables

### Don't forget to plot the values (present, future and difference) along the way 


# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("stringr")
library("reshape2")
library("tidyverse")
library("viridis")
library("ncdf4")
library("RColorBrewer")
library("oce")
library("viridis")

# --------------------------------------------------------

WD <- getwd()

# Get world coastline
setwd("/Users/fabiobenedetti/Desktop/PostDocs/ETHZ/OVERSEE/data/")
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
ggplot() + coast + theme_linedraw()


##### 1°) First, take care of the monthly climatologiesd for 2012 & 2100

### Go back to netcdf files directory
setwd(WD)
vars <- c("SST","SSS","STOTCHL","SSiO3")
years <- c("2012","2100")
months <- c("Jan","Apr","Jul","Oct")

# For testing the code below:
#v <- "SSiO3"
#y <- "2012"

for(y in years) {

		data <- lapply(vars, function(v) {
	
				# Useless message
				message(paste("Preparing cliamtologies for ", v, " for year ", y, sep = ""))
				nc <- nc_open(paste("GFDL_ESM2M_monthly_",v,"_",y,".nc", sep = ""))
				# print(nc)
				# names(nc$var) # contains the data from 9 different models 
				var <- ncvar_get(nc, v)
				# class(var) # array
				# dim(var) # lon - lat - months

				# Melt 
				mvar <- melt(var)
				#dim(msst)
				#head(msst) # lon/lat/var
				colnames(mvar) <- c("x","y","month",v)
				# Convert lats and provide cell id
				mvar$y <- mvar$y - 90
				mvar$id <- factor(paste(mvar$x, mvar$y, sep = "_"))
				# summary(mvar)
				
				# if v == Chl or Silicates, need to convert them to proper units and log-transform them 
				if(v == "STOTCHL") {
						### NOTE: STOTCHL in GFDL_ESM2M is expressed in mgChla.m-3 while siloicates climatologies are expressed in µmol/L
						### SO need to multiply by 10^3 and log1p
						mvar$STOTCHL <- log10(mvar$STOTCHL)
						# summary(mvar)
				} # eo 1st if loop
				
				if(v == "SSiO3") {
						### NOTE: SSiO3 in GFDL_ESM2M is expressed in mol.m-3 which is the unit you are using to train SDMs
						mvar$SSiO3 <- (mvar$SSiO3)*10^3 
						# And log it
						mvar$SSiO3 <- log1p(mvar$SSiO3)
						# summary(mvar$SSiO3)
				} # eo 2nd if loop

				# DCAST
				var2 <- dcast(mvar, id + x + y ~ month, value.var = v)
				# Change colnames
				colnames(var2)[c(4:length(var2))] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

				# But first, need to convert it so the 0° lon is centered on the Atl
				coordinates(var2) <- ~ x + y
				gridded(var2) <- TRUE
				stck <- stack(var2)
				crs(stck) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
				stck <- rotate(stck) # to center the projetcion on the Atlantic Ocean

				# And re-convert to DDF
				var2 <- as.data.frame(stck, xy = T)
				# colnames(var2)
				# dim(var2)
				# summary(var2)
				
				# Return all months to start with
				return(var2)
	
			}
	
		) # eo lapply
		
		# Cbind all?
		table <- do.call(cbind, data)
		rm(data)
		gc()
		# dim(table)
		# head(table)
		# colnames(table)
		# Change colnames: 4 to 15 -> SST
		# 19 - 30 -> SSS
		# 34 - 45 -> STOTCHL
		# 49 - 60 -> SSiO3 
		colnames(table)[c(4:15)] <- c("SST_Jan","SST_Feb","SST_Mar","SST_Apr","SST_May","SST_Jun","SST_Jul","SST_Aug","SST_Sept","SST_Oct","SST_Nov","SST_Dec")
		colnames(table)[c(19:30)] <- c("SSS_Jan","SSS_Feb","SSS_Mar","SSS_Apr","SSS_May","SSS_Jun","SSS_Jul","SSS_Aug","SSS_Sept","SSS_Oct","SSS_Nov","SSS_Dec")
		colnames(table)[c(34:45)] <- c("logChl_Jan","logChl_Feb","logChl_Mar","logChl_Apr","logChl_May","logChl_Jun","logChl_Jul","logChl_Aug","logChl_Sept","logChl_Oct","logChl_Nov","logChl_Dec")
		colnames(table)[c(49:60)] <- c("logSiO2_Jan","logSiO2_Feb","logSiO2_Mar","logSiO2_Apr","logSiO2_May","logSiO2_Jun","logSiO2_Jul","logSiO2_Aug","logSiO2_Sept","logSiO2_Oct","logSiO2_Nov","logSiO2_Dec")
		
		### Split into 4 monthly dataframes
		clim_jan <- table[,c("id","x","y",colnames(table)[grep("Jan",colnames(table))])]
		clim_apr <- table[,c("id","x","y",colnames(table)[grep("Apr",colnames(table))])]
		clim_jul <- table[,c("id","x","y",colnames(table)[grep("Jul",colnames(table))])]
		clim_oct <- table[,c("id","x","y",colnames(table)[grep("Oct",colnames(table))])]
		
		### Compute annual climatologies for all 4 variables and compute difference between max and min values 
		require("matrixStats")
		table$annual_SST <- rowMeans(table[,c(4:15)])
		table$annual_SSS <- rowMeans(table[,c(19:30)])
		table$annual_logChl <- rowMeans(table[,c(34:45)])
		table$annual_logSiO2 <- rowMeans(table[,c(49:60)])
		# Find max and min values of SST to compute diff	
		table$maxSST <- apply(table[,4:15], 1, max) 
		table$minSST <- apply(table[,4:15], 1, min)
		table$dSST <- table$maxSST - table$minSST
		# Provide it to the 4 monthly climatologies too
		clim_jan$dSST <- table$dSST 
		clim_apr$dSST <- table$dSST 
		clim_jul$dSST <- table$dSST 
		clim_oct$dSST <- table$dSST 
		# Gather annual climatologies in one data.frame
		clim_annual <- table[,c("id","x","y","annual_SST","annual_SSS","annual_logChl","annual_logSiO2","dSST")]
		
		### OK, now that you have the annual clim and the 4 monthly clims, simply map the values
		# First, the annual climatologies
		map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = annual_SST), data = clim_annual) + coast + 
				scale_fill_viridis(name = "Annual SST\n(°C)") + coord_quickmap() + theme_bw()
				
		map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = annual_SSS), data = clim_annual) + coast + 
				scale_fill_viridis(name = "Annual SSS") + coord_quickmap() + theme_bw()
				
		map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = annual_logChl), data = clim_annual) + coast + 
				scale_fill_viridis(name = "Annual Chl-a\nlog(mg.m3)", limits = c(-2,2)) + coord_quickmap() + theme_bw()
	
		map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = annual_logSiO2), data = clim_annual) + coast + 
				scale_fill_viridis(name = "Annual SiO2\nlog(mumoles/L)") + coord_quickmap() + theme_bw()
				
		map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = dSST), data = clim_annual) + coast + 
				scale_fill_viridis(name = "SST range\n(°C)") + coord_quickmap() + theme_bw()
				
		# Save them in proper dir
		setwd(WD)
		ggsave(plot = map1, filename = paste("map_annual_SST_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
		ggsave(plot = map2, filename = paste("map_annual_SSS_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
		ggsave(plot = map3, filename = paste("map_annual_logChl_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
		ggsave(plot = map4, filename = paste("map_annual_logSiO2_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
		ggsave(plot = map5, filename = paste("map_annual_dSST_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
		rm(map1, map2, map3, map4, map5)
		
		# Second, map the monthly climatologies (using a for loop)
		for(m in months) {
			
				if(m == "Jan") {
					data <- clim_jan
				} else if (m == "Apr") {
					data <- clim_apr
				} else if (m == "Jul") {
					data <- clim_jul
				} else {
					data <- clim_oct
				} # eo if else loop
			
				# Change colnames to remove month 
				colnames(data)[c(4:7)] <- c("SST","SSS","logChl","logSiO2")
			
				# Map like above
				map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = SST), data = data) + coast + 
						scale_fill_viridis(name = "SST\n(°C)") + coord_quickmap() + theme_bw()
				
				map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = SSS), data = data) + coast + 
						scale_fill_viridis(name = "SSS") + coord_quickmap() + theme_bw()
				
				map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = logChl), data = data) + coast + 
						scale_fill_viridis(name = "Chl-a\nlog(mg.m3)", limits = c(-2,2)) + coord_quickmap() + theme_bw()
					
				map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = logSiO2), data = data) + coast + 
						scale_fill_viridis(name = "SiO2\nlog(mumoles/L)") + coord_quickmap() + theme_bw()
						
				# Save maps
				ggsave(plot = map1, filename = paste("map_monthly_SST_",m,"_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
				ggsave(plot = map2, filename = paste("map_monthly_SSS_",m,"_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
				ggsave(plot = map3, filename = paste("map_monthly_logChl_",m,"_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
				ggsave(plot = map4, filename = paste("map_monthly_logSiO2_",m,"_",y,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
				
				rm(data, map1, map2, map3, map4)
			
		} # eo for loop - m in months
		
		# Save annual and monthly climatologies as .txt files
		message(paste("Saving cliamtologies for year ", y, sep = ""))
		setwd(WD)
		write.table(clim_annual, paste("clim_GFDL_ESM2M_annual_",y,".txt", sep = ""), sep = ";")
		write.table(clim_jan, paste("clim_GFDL_ESM2M_monthly_Jan_",y,".txt", sep = ""), sep = ";")
		write.table(clim_apr, paste("clim_GFDL_ESM2M_monthly_Apr_",y,".txt", sep = ""), sep = ";")
		write.table(clim_jul, paste("clim_GFDL_ESM2M_monthly_Jul_",y,".txt", sep = ""), sep = ";")
		write.table(clim_oct, paste("clim_GFDL_ESM2M_monthly_Oct_",y,".txt", sep = ""), sep = ";")
		
		# Make room
		rm(table, clim_jan, clim_apr, clim_jul, clim_oct, clim_annual)
		gc()
		
} # eo for loop - y in years



# --------------------------------------------------------

##### 12/11/10: Using the climatologies obtained from the code above, compute the difference between 2100 and 2012 for each variable
setwd(WD)
# There are 5 time periods (annual, Jan, Apr, Jul, Oct) and 5 variables (SST, dSST, SSS, logChl, logSiO2)
periods <- c("annual","Jan","Apr","Jul","Oct")
# p <- "annual"
# p <- "Jan"

for (p in periods) {
		
		# Message
		message(paste("Plotting difference in parameters for ", p, sep = ""))
		if(p == "annual") {
			
			present <- read.table(paste("clim_GFDL_ESM2M_annual_2012.txt", sep = ""), h = T, sep = ";")
			future <- read.table(paste("clim_GFDL_ESM2M_annual_2100.txt", sep = ""), h = T, sep = ";")
			
		} else if (p == "Jan") {
			
			present <- read.table(paste("clim_GFDL_ESM2M_monthly_Jan_2012.txt", sep = ""), h = T, sep = ";")
			future <- read.table(paste("clim_GFDL_ESM2M_monthly_Jan_2100.txt", sep = ""), h = T, sep = ";")
			
		} else if (p == "Apr") {
			
			present <- read.table(paste("clim_GFDL_ESM2M_monthly_Apr_2012.txt", sep = ""), h = T, sep = ";")
			future <- read.table(paste("clim_GFDL_ESM2M_monthly_Apr_2100.txt", sep = ""), h = T, sep = ";")
			
		} else if (p == "Jul") {
			
			present <- read.table(paste("clim_GFDL_ESM2M_monthly_Jul_2012.txt", sep = ""), h = T, sep = ";")
			future <- read.table(paste("clim_GFDL_ESM2M_monthly_Jul_2100.txt", sep = ""), h = T, sep = ";")
			
		} else {
			
			present <- read.table(paste("clim_GFDL_ESM2M_monthly_Oct_2012.txt", sep = ""), h = T, sep = ";")
			future <- read.table(paste("clim_GFDL_ESM2M_monthly_Oct_2100.txt", sep = ""), h = T, sep = ";")
			
		} # eo else if loop
		
		# Compute differences and map (use 'present' ddf as main table)
		if(p == "annual") {
			
			# Compute differences
			present$delta_SST <- (future$annual_SST) - (present$annual_SST)
			present$delta_SSS <- (future$annual_SSS) - (present$annual_SSS)
			present$delta_dSST <- (future$dSST) - (present$dSST)
			present$delta_logChl <- (future$annual_logChl) - (present$annual_logChl)
			present$delta_logSiO2 <- (future$annual_logSiO2) - (present$annual_logSiO2)
			# Map
			map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_SST), data = present) + coast + 
					scale_fill_gradient2(name = "Delta SST\n(°C)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
			
			map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_SSS), data = present) + coast + 
					scale_fill_gradient2(name = "Delta SSS", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
			
			map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_logChl), data = present) + coast + 
					scale_fill_gradient2(name = "Delta Chl-a\nlog(mg.m3)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
				
			map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_logSiO2), data = present) + coast + 
					scale_fill_gradient2(name = "Delta SiO2\nlog(mumoles/L)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()

			map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_dSST), data = present) + coast + 
					scale_fill_gradient2(name = "Delta dSST\n(°C)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
					
			# Save maps
			ggsave(plot = map1, filename = paste("map_diff_annual_SST.pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map2, filename = paste("map_diff_annual_SSS.pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map3, filename = paste("map_diff_annual_logChl.pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map4, filename = paste("map_diff_annual_logSiO2.pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map5, filename = paste("map_diff_annual_dSST.pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			
			rm(map1, map2, map3, map4, map5)
	
			
		} else {
			
			# Change colnames
			colnames(present)[c(4:7)] <- c("SST","SSS","logChl","logSiO2")
			colnames(future)[c(4:7)] <- c("SST","SSS","logChl","logSiO2")
			# Compute difference (except for dSST because it is annual, not monthly)
			present$delta_SST <- (future$SST) - (present$SST)
			present$delta_SSS <- (future$SSS) - (present$SSS)
			present$delta_logChl <- (future$logChl) - (present$logChl)
			present$delta_logSiO2 <- (future$logSiO2) - (present$logSiO2)
			
			# Map
			map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_SST), data = present) + coast + 
					scale_fill_gradient2(name = "Delta SST\n(°C)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
			
			map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_SSS), data = present) + coast + 
					scale_fill_gradient2(name = "Delta SSS", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
			
			map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_logChl), data = present) + coast + 
					scale_fill_gradient2(name = "Delta Chl-a\nlog(mg.m3)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
				
			map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta_logSiO2), data = present) + coast + 
					scale_fill_gradient2(name = "Delta SiO2\nlog(mumoles/L)", low = "#2166ac", high = "#b2182b") + coord_quickmap() + theme_bw()
					
			# Save maps
			ggsave(plot = map1, filename = paste("map_diff_monthly_SST_",p,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map2, filename = paste("map_diff_monthly_SSS_",p,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map3, filename = paste("map_diff_monthly_logChl_",p,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			ggsave(plot = map4, filename = paste("map_diff_monthly_logSiO2_",p,".pdf", sep = ""), dpi = 300, height = 6, width = 8)	
			
			rm(map1, map2, map3, map4)
			
		} # eo if loop - annual versus monthly 
		
	
} # eo for loop p in periods

### Check out the maps printed in your dir.
### All gut.

# --------------------------------------------------------

##### 12/11/18: Add contemporary annual (after calculating it) and monthly PAR to your 10 climatologies.

WD <- getwd()

# Get world coastline
setwd("/Users/fabiobenedetti/Desktop/PostDocs/ETHZ/OVERSEE/data/")
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

### Go to the directory with the global 1° climatologies
setwd("/Users/fabiobenedetti/Desktop/PostDocs/ETHZ/OVERSEE/data/global_monthly_clims/")
par <- read.table("PAR_stack_1d.txt", h = T, sep = ";")
dim(par)
str(par)

# First, compute annual PAR climatology and map it 
require("matrixStats")
par$annual <- rowMeans(as.matrix(par[,c(3:14)]), na.rm = T)
summary(par)
quartz()
ggplot() + geom_raster(aes(x = x, y = y, fill = annual), data = par) + 
	coast + scale_fill_viridis(name = "Annual PAR") + coord_quickmap() + theme_bw()
### Some weird values at high latitudes...identify the rows which have less than 4 or 5 monthly PAR values and remove the from the annual climatology (not robust enough)
par$na_count <- apply(par[,c(3:14)], 1, function(x) sum(is.na(x)))

# Rows with par$na_count >= 5 --> NA
par[par$na_count >= 5,"annual"] <- NA
### Re-map
quartz()
ggplot() + geom_raster(aes(x = x, y = y, fill = annual), data = par) + 
	coast + scale_fill_viridis(name = "Annual PAR") + coord_quickmap() + theme_bw()


### Now, load climatologies from GFDL_TOPAZ and cbind the PAR annual and monthly clims
setwd(WD)

clims <- dir()[grep("clim_GFDL_ESM2M_", dir())]
# c <- "clim_GFDL_ESM2M_monthly_Jan_2012.txt"
# c <- "clim_GFDL_ESM2M_annual_2100.txt"

for(c in clims) {
	
		# Load climatology
		clim <- read.table(c, h = T, sep = ";")
		
		# As a function of the clim name, load monthly or annual PAR
		if(c %in% c("clim_GFDL_ESM2M_annual_2100.txt","clim_GFDL_ESM2M_annual_2012.txt")) {
			
			clim$PAR <- par$annual
			
		} else if(c %in% c("clim_GFDL_ESM2M_monthly_Jan_2012.txt","clim_GFDL_ESM2M_monthly_Jan_2100.txt")) {
			
			clim$PAR <- par$PAR_Jan
			
		} else if (c %in% c("clim_GFDL_ESM2M_monthly_Apr_2012.txt","clim_GFDL_ESM2M_monthly_Apr_2100.txt")) {
			
			clim$PAR <- par$PAR_Apr
			
		} else if(c %in% c("clim_GFDL_ESM2M_monthly_Jul_2012.txt","clim_GFDL_ESM2M_monthly_Jul_2100.txt")) {
			
			clim$PAR <- par$PAR_Jul
			
		} else if(c %in% c("clim_GFDL_ESM2M_monthly_Oct_2012.txt","clim_GFDL_ESM2M_monthly_Oct_2100.txt")) {
			
			clim$PAR <- par$PAR_Oct
			
		} # eo if else loop
		
		# Change the colnames and save as table
		colnames(clim)[c(4:7)] <- c("SST","SSS","logChl","logSiO2")
		write.table(x = clim, file = str_replace(c, "clim_", "clim2_"), sep = ";")
	
} # eo for loop


### Ok, you seem to be set for some quick climate change predictions




