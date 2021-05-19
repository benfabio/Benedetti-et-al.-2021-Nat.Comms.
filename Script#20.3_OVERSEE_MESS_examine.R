
##### 20/05/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### RSCRIPTBATCH CMD for : 
#	- loading the MESS maps for each species, for each pool, and for the 2 kingdoms separately
#	- compute the frequency of extrapolation (average TOTAL) and identify the main env variables causing the extrapolation
#	- map avergae TOTAL value, for each month, of each env variable pool (12*5*2 = 120 maps)
 
### Last update : 20/05/2019

# --------------------------------------------------------------------------------------------------------------------------------

library("reshape2")
library("tidyverse")
library("viridis")
library("RColorBrewer")
library("modEvA")
library("raster")
library("maptools")
library("marmap")

# --------------------------------------------------------------------------------------------------------------------------------

# Master directory 
WD <- getwd()
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4","p5")

# For testing
# p <- "p1"
# m <- "apr"

### For zooplankton first
setwd(paste(WD,"/species_data_v9v3.1/total_background/niche.modelling/", sep = ""))
zoo.wd <- getwd()
for(p in pools) {
	
	# Go to 
	setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
	for(m in months) {
		
			message(paste("Plotting MESS for ", m, " of pool ", p, sep = ""))
			# Vector of species files corresponding to month m
			species <- dir()[grep(m,dir())]
			# Load species files 
			require("parallel")
			# sp <- "Vibilia_armata_apr.Rdata"
			res <- mclapply(species, function(sp) {
						d <- get(load(sp))
						# Add sp name (just in case)
						spp <- str_replace_all(sp, ".Rdata", "")
						spp <- str_replace_all(spp, paste("_",m,sep=""), "")
						d$species <- spp
						# Add cell id
						d$id <- paste(d$x, d$y, sep = "_")
						return(d)
					}, mc.cores = 25
			) # eo mclapply
			# Rbind
			table <- do.call(rbind, res)
			# dim(table); head(table)		
			rm(res)
			
			# Compute: mean TOTAL and frequency of TOTAL < 0, for each cell id
			require("dplyr")
			ddf <- data.frame( table %>% group_by(id) %>% 
					summarise(x = unique(x), y = unique(y), 
						mean = mean(TOTAL), f = (length(TOTAL[TOTAL < 0]))/length(unique(species)) ) 
			) # eo ddf
			# head(ddf); summary(ddf) # Ok, ready for mapping
			
			# Load world coastline etc.
			library("maps")
			pac <- map_data(map ="world2")
			ddf$x2 <- ddf$x # not compulsory, but just in case you want to keep the original vector of longitudes
			ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360
			
			map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = ddf) + 
					scale_fill_gradient2(name = "Mean MESS", high = "#313695", mid = "white", low = "#a50026") + 
					geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
					coord_quickmap() + 
  					scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
             	   		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
					scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
		     	   		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
  					theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
						panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

			
			map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = f), data = ddf) + 
					scale_fill_viridis(name = "Frequency of\nMESS < 0", discrete = F) + 
					geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
					coord_quickmap() + 
  					scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
              	  		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
					scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
		      	  		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
  					theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
						panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

			
			# Compute, for each variable, hwo many times they appear as driving predictor where TOTAl < 0
			vars <- table[table$TOTAL < 0,]
			# Plot density of obs per var
			plot <- ggplot(vars, aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
						scale_fill_brewer(name = "Basin", palette = "Spectral") + 
						geom_density(position = "stack", alpha = 0.5) + ylab("Density") + xlab("") + 
						theme_linedraw()
			#
			setwd(zoo.wd)
			ggsave(plot = plot, filename = paste("plot_MoD_",m,"_",p,".jpg", sep = ""), dpi = 300, height = 4, width = 6)
			ggsave(plot = map1, filename = paste("map_avg_MESS_",m,"_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
			ggsave(plot = map2, filename = paste("map_freq_",m,"_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
			
			# Clean 
			rm(plot, map1, map2, vars, ddf, table); gc()
			setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
		
	} # eo second for loop
	
} # eo first loop


### For phytoplankton 
setwd(paste(WD,"/phytoplankton_15_01_19/total_background/species_data/niche.modelling/", sep = ""))
phyto.wd <- getwd()

for(p in pools) {
	
	# Go to 
	setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
	for(m in months) {
			
			# Vector of species files corresponding to month m
			message(paste("Plotting MESS for ", m, " of pool ", p, sep = ""))
			species <- dir()[grep(m,dir())]
			# Load species files 
			require("parallel")
			# sp <- "Vibilia_armata_apr.Rdata"
			res <- mclapply(species, function(sp) {
						d <- get(load(sp))
						# Add sp name (just in case)
						spp <- str_replace_all(sp, ".Rdata", "")
						spp <- str_replace_all(spp, paste("_",m,sep=""), "")
						d$species <- spp
						# Add cell id
						d$id <- paste(d$x, d$y, sep = "_")
						return(d)
					}, mc.cores = 25
			) # eo mclapply
			# Rbind
			table <- do.call(rbind, res)
			# dim(table); head(table)		
			rm(res)
			
			# Compute: mean TOTAL and frequency of TOTAL < 0, for each cell id
			require("dplyr")
			ddf <- data.frame( table %>% group_by(id) %>% 
					summarise(x = unique(x), y = unique(y), 
						mean = mean(TOTAL), f = (length(TOTAL[TOTAL < 0]))/length(unique(species)) ) 
			) # eo ddf
			# head(ddf); summary(ddf) # Ok, ready for mapping
			
			# Load world coastline etc.
			library("maps")
			pac <- map_data(map ="world2")
			ddf$x2 <- ddf$x # not compulsory, but just in case you want to keep the original vector of longitudes
			ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360
			
			map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = ddf) + 
					scale_fill_gradient2(name = "Mean MESS", high = "#313695", mid = "white", low = "#a50026") + 
					geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
					coord_quickmap() + 
  					scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
             	   		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
					scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
		     	   		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
  					theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
						panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

			
			map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = f), data = ddf) + 
					scale_fill_viridis(name = "Frequency of\nMESS < 0", discrete = F) + 
					geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
					coord_quickmap() + 
  					scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
              	  		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
					scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
		      	  		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
  					theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
						panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

			
			# Compute, for each variable, hwo many times they appear as driving predictor where TOTAl < 0
			vars <- table[table$TOTAL < 0,]
			# Plot density of obs per var
			plot <- ggplot(vars, aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
						scale_fill_brewer(name = "Basin", palette = "Spectral") + 
						geom_density(position = "stack", alpha = 0.5) + ylab("Density") + xlab("") + 
						theme_linedraw()
			#
			setwd(phyto.wd)
			ggsave(plot = plot, filename = paste("plot_MoD_",m,"_",p,".jpg", sep = ""), dpi = 300, height = 4, width = 6)
			ggsave(plot = map1, filename = paste("map_avg_MESS_",m,"_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
			ggsave(plot = map2, filename = paste("map_freq_",m,"_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
			
			# Clean 
			rm(plot, map1, map2, vars, ddf, table); gc()
			setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
		
	} # eo second for loop
	
	
} # eo first loop


### Compute annual mean of TOTAL and frequency of < 0
# Zoo first again
setwd(paste(WD,"/species_data_v9v3.1/total_background/niche.modelling/", sep = ""))
zoo.wd <- getwd()
for(p in pools) {
	
		# Go to pool's dir
		message(paste("Plotting annual MESS for pool ", p, sep = ""))
		setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
		species <- dir()
		
		require("parallel")
		res <- mclapply(species, function(sp) {
					d <- get(load(sp))
					d$id <- paste(d$x, d$y, sep = "_")
					return(d)
				}, mc.cores = 27
		) # eo mclapply
		# Rbind
		# table <- do.call(rbind, res)
		# ??rbind.fill; ??bind_rows
		require("dplyr")
		table <- dplyr::bind_rows(res) # much much fasta ! 
		# dim(table); head(table)		
		rm(res)
		
		# Compute: mean TOTAL and frequency of TOTAL < 0, for each cell id
		l <- length(unique(species))
		ddf <- data.frame(table %>% group_by(id) %>% 
				summarise(x = unique(x), y = unique(y), 
					mean = mean(TOTAL), f = (length(TOTAL[TOTAL < 0]))/l ) 
		) # eo ddf
		# head(ddf); summary(ddf) # Ok, ready for mapping
		
		# Load world coastline etc.
		library("maps")
		pac <- map_data(map = "world2")
		ddf$x2 <- ddf$x # not compulsory, but just in case you want to keep the original vector of longitudes
		ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360
		
		map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = ddf) + 
				scale_fill_gradient2(name = "Mean MESS", high = "#313695", mid = "white", low = "#a50026") + 
				geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
				coord_quickmap() + 
				scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
         	   		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
				scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
	     	   		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
				theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
					panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

		
		map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = f), data = ddf) + 
				scale_fill_viridis(name = "Frequency of\nMESS < 0", discrete = F) + 
				geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
				coord_quickmap() + 
				scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
          	  		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
				scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
	      	  		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
				theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
					panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

		
		# Compute, for each variable, hwo many times they appear as driving predictor where TOTAl < 0
		vars <- table[table$TOTAL < 0,]
		# Plot density of obs per var
		plot <- ggplot(vars, aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
					scale_fill_brewer(name = "Basin", palette = "Spectral") + 
					geom_density(position = "stack", alpha = 0.5) + ylab("Density") + xlab("") + 
					theme_linedraw()
		#
		setwd(zoo.wd)
		ggsave(plot = plot, filename = paste("plot_MoD_","annual","_",p,".jpg", sep = ""), dpi = 300, height = 4, width = 6)
		ggsave(plot = map1, filename = paste("map_avg_MESS_","annual","_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
		ggsave(plot = map2, filename = paste("map_freq_","annual","_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
		
		# Clean 
		rm(plot, map1, map2, vars, ddf, table); gc()
		setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
	
} # eo for loop
	
# Phyto now
setwd(paste(WD,"/phytoplankton_15_01_19/total_background/species_data/niche.modelling/", sep = ""))
phyto.wd <- getwd()
for(p in pools) {
	
		# Go to pool's dir
		message(paste("Plotting annual MESS for pool ", p, sep = ""))
		setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
		species <- dir()
		
		require("parallel")
		res <- mclapply(species, function(sp) {
					d <- get(load(sp))
					d$id <- paste(d$x, d$y, sep = "_")
					return(d)
				}, mc.cores = 27
		) # eo mclapply
		# Rbind
		# table <- do.call(rbind, res)
		# ??rbind.fill; ??bind_rows
		require("dplyr")
		table <- dplyr::bind_rows(res) # much much fasta ! 
		# dim(table); head(table)		
		rm(res)
		
		# Compute: mean TOTAL and frequency of TOTAL < 0, for each cell id
		l <- length(unique(species))
		ddf <- data.frame(table %>% group_by(id) %>% 
				summarise(x = unique(x), y = unique(y), 
					mean = mean(TOTAL), f = (length(TOTAL[TOTAL < 0]))/l ) 
		) # eo ddf
		# head(ddf); summary(ddf) # Ok, ready for mapping
		
		# Load world coastline etc.
		library("maps")
		pac <- map_data(map ="world2")
		ddf$x2 <- ddf$x # not compulsory, but just in case you want to keep the original vector of longitudes
		ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360
		
		map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = mean), data = ddf) + 
				scale_fill_gradient2(name = "Mean MESS", high = "#313695", mid = "white", low = "#a50026") + 
				geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
				coord_quickmap() + 
				scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
         	   		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
				scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
	     	   		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
				theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
					panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

		
		map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = f), data = ddf) + 
				scale_fill_viridis(name = "Frequency of\nMESS < 0", discrete = F) + 
				geom_polygon(aes(x = long, y = lat, group = group), data = pac, fill = "grey65", colour = "black", size = 0.3) + 
				coord_quickmap() + 
				scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360), 
          	  		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
				scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90), 
	      	  		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) + 
				theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
					panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 

		
		# Compute, for each variable, hwo many times they appear as driving predictor where TOTAl < 0
		vars <- table[table$TOTAL < 0,]
		# Plot density of obs per var
		plot <- ggplot(vars, aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
					scale_fill_brewer(name = "Basin", palette = "Spectral") + 
					geom_density(position = "stack", alpha = 0.5) + ylab("Density") + xlab("") + 
					theme_linedraw()
		#
		setwd(phyto.wd)
		ggsave(plot = plot, filename = paste("plot_MoD_","annual","_",p,".jpg", sep = ""), dpi = 300, height = 4, width = 6)
		ggsave(plot = map1, filename = paste("map_avg_MESS_","annual","_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
		ggsave(plot = map2, filename = paste("map_freq_","annual","_",p,".jpg", sep = ""), dpi = 300, height = 5, width = 7)
		
		# Clean 
		rm(plot, map1, map2, vars, ddf, table); gc()
		setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
	
} # eo for loop



# ---------------------------------------------------------------

### 21/05/2019: From the annual estimates of MESS frequency, identify those grid cells where f(MESS<0) is > 5% 
### (i.e. those cells where extrapolation occurs more than 5% of the time across each species and each month) 
### Save those cells, and remover them from the pool's average annual maps of diversity (prior to computing the ensemble diversity map)

setwd(paste(WD,"/species_data_v9v3.1/total_background/niche.modelling/", sep = ""))
zoo.wd <- getwd()
for(p in pools) {
	
		# Go to pool's dir
		message(paste("Extracting cells with f(MESS<0) > 5% for pool ", p, sep = ""))
		setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
		species <- dir()
		
		require("parallel")
		res <- mclapply(species, function(sp) {
					d <- get(load(sp))
					d$id <- paste(d$x, d$y, sep = "_")
					return(d)
				}, mc.cores = 27
		) # eo mclapply
		require("dplyr")
		table <- dplyr::bind_rows(res) 
		rm(res)
		
		# Compute: mean TOTAL and frequency of TOTAL < 0, for each cell id
		l <- length(unique(species))
		ddf <- data.frame(table %>% group_by(id) %>% 
				summarise(x = unique(x), y = unique(y), 
					mean = mean(TOTAL), f = (length(TOTAL[TOTAL < 0]))/l ) 
		) # eo ddf
		ddf$x2 <- ddf$x # not compulsory, but just in case you want to keep the original vector of longitudes
		ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360
	
		# Compute, for each variable, hwo many times they appear as driving predictor where TOTAl < 0
		setwd(WD)
		message(paste("Saving data for pool ", p, sep = ""))
		save(ddf, file = paste("mess_values_zoo_",p,".Rdata", sep = "") )
		
		# Clean 
		rm(ddf, table); gc()
		setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
	
} # eo for loop
	
	
# Phyto now
setwd(paste(WD,"/phytoplankton_15_01_19/total_background/species_data/niche.modelling/", sep = ""))
phyto.wd <- getwd()
for(p in pools) {
	
		# Go to pool's dir
		message(paste("Extracting cells with f(MESS<0) > 5% for pool ", p, sep = ""))
		setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
		species <- dir()
		
		require("parallel")
		res <- mclapply(species, function(sp) {
					d <- get(load(sp))
					d$id <- paste(d$x, d$y, sep = "_")
					return(d)
				}, mc.cores = 27
		) # eo mclapply
		require("dplyr")
		table <- dplyr::bind_rows(res) 
		rm(res)
		
		# Compute: mean TOTAL and frequency of TOTAL < 0, for each cell id
		l <- length(unique(species))
		ddf <- data.frame(table %>% group_by(id) %>% 
				summarise(x = unique(x), y = unique(y), 
				mean = mean(TOTAL), f = (length(TOTAL[TOTAL < 0]))/l ) 
		) # eo ddf
		
		# Load world coastline etc.
		ddf$x2 <- ddf$x # not compulsory, but just in case you want to keep the original vector of longitudes
		ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360
				
		# Compute, for each variable, hwo many times they appear as driving predictor where TOTAl < 0
		setwd(WD)
		message(paste("Saving data for pool ", p, sep = ""))
		save(ddf, file = paste("mess_values_phyto_",p,".Rdata", sep = "") )
		
		# Clean 
		rm(ddf, table); gc()
		setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
	
} # eo for loop

# Check one file out
data <- get(load("mess_values_phyto_p1.Rdata"))
dim(data)
head(data)
summary(data)
# nrow(data[data$f > 0.05,])

# ---------------------------------------------------------------


