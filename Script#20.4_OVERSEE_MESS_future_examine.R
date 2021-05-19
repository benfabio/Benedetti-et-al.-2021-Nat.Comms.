
##### 02/09/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- retrieve the MESS data for each species, for each pool, and for the 2 kingdoms separately for RCP 8.5 GFDL-ESM2M
#	- compute the frequency of extrapolation (average TOTAL) and identify the main env variables causing the extrapolation
 
### Last update : 03/09/2019

# --------------------------------------------------------------------------------------------------------------------------------

library("reshape2")
library("tidyverse")
library("viridis")
library("RColorBrewer")
library("modEvA")
library("raster")
library("maptools")
library("marmap")

world2 <- map_data(map = "world2")

# --------------------------------------------------------------------------------------------------------------------------------

# Master directory 
WD <- getwd()
# Vector of months
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")

# For testing
#p <- "p1"
#m <- "Apr"

### For phytoplankton first
setwd(paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future_rcp85/", sep = ""))
phyto.wd <- getwd()
# dir()
res <- lapply(pools, function(p) {
            
            message(paste("", sep = ""))
            message(paste("Retrieving MESS data for p || ",p, sep = ""))
            message(paste("", sep = ""))
            setwd(paste(phyto.wd,"/","MESS_",p,"/", sep = ""))
            
            require("parallel")
            mess <- mclapply(dir(), function(f) {
                        # Retrieve the species name and the mont from the filename
                        # f <- dir()[1]
                        chars <- data.frame(do.call(rbind, strsplit(f, "_")))
                        sp <- paste(chars$X1, chars$X2, sep = "_")
                        message(paste("Retrieving MESS data for ",sp, sep = ""))
                        m <- str_replace_all(as.character(chars$X3), ".Rdata", "")
                        d <- get(load(f))
                        d$species <- sp
                        d$month <- m 
                        # Return
                        return(d)
                }, mc.cores = 25 
            ) # eo mclapply
            
            # Rbind
            require("dplyr")
            messy <- bind_rows(mess)
            # head(messy)
            rm(mess); gc()
            messy$pool <- p 
            messy$id <- paste(messy$x, messy$y, sep = "_")
            
            # Return
            return(messy)
        
        } # eo FUN
    
) # eo lapply 
# Rbind
mess.phyto <- bind_rows(res)
dim(mess.phyto); head(mess.phyto); length(unique(mess.phyto$species))
rm(res);gc()

### For zooplankton second
setwd(paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future_rcp85/", sep = ""))
zoo.wd <- getwd()
res <- lapply(pools, function(p) {
            
            message(paste("", sep = ""))
            message(paste("Retrieving MESS data for p || ",p, sep = ""))
            message(paste("", sep = ""))
            setwd(paste(zoo.wd,"/","MESS_",p,"/", sep = ""))
            
            require("parallel")
            mess <- mclapply(dir(), function(f) {
                
                        # Retrieve the species name and the mont from the filename
                        # f <- dir()[6062]
                        chars <- data.frame(do.call(rbind, strsplit(f, "_")))
                        if(length(chars) == 3) {
                            sp <- paste(chars$X1, chars$X2, sep = "_")
                            message(paste("Retrieving MESS data for ",sp, sep = ""))
                            m <- str_replace_all(as.character(chars$X3), ".Rdata", "")
                            d <- get(load(f))
                            d$species <- sp
                            d$month <- m 
                        } else if (length(chars) == 4) { 
                            sp <- paste(chars$X1, chars$X2, chars$X3, sep = "_")
                            m <- str_replace_all(as.character(chars$X4), ".Rdata", "")
                            d <- get(load(f))
                            d$species <- sp
                            d$month <- m 
                        } # eo if else loop
                        
                        # Return
                        return(d)
                        
                }, mc.cores = 25 
            ) # eo mclapply
            
            # Rbind
            require("dplyr")
            messy <- bind_rows(mess)
            # head(messy)
            rm(mess); gc()
            messy$pool <- p 
            messy$id <- paste(messy$x, messy$y, sep = "_")
            
            # Return
            return(messy)
        
        } # eo FUN
    
) # eo lapply 
# Rbind
mess.zoo <- bind_rows(res)
dim(mess.zoo); head(mess.zoo); unique(mess.zoo$month)
rm(res);gc()
# Add cell id 


### OK, for phytoplankton and zooplankton, compute frequency of extrapolation 
d.phyto <- data.frame(mess.phyto %>% 
    group_by(id,month) %>% 
    summarise(x = unique(x), y = unique(y), mean = mean(TOTAL,na.rm=T), freq = length(TOTAL[TOTAL<0]) ) 
) # eo ddf
# summary(d.phyto)
# head(d.phyto)
# Compute % frequency by dividing by number of spp * n pools
d.phyto$freq <- (d.phyto$freq)/(length(unique(mess.phyto$species))*4)

# Same with zoo
d.zoo <- data.frame(mess.zoo %>% 
    group_by(id,month) %>% 
    summarise(x = unique(x), y = unique(y), mean = mean(TOTAL,na.rm=T), freq = length(TOTAL[TOTAL<0]) ) 
) # eo ddf
# summary(d.zoo)
# head(d.zoo)
d.zoo$freq <- (d.zoo$freq)/(length(unique(mess.zoo$species))*4)

### Now, from the monthly average TOTAL & freq, compute the annual averages
a.phyto <- data.frame(d.phyto %>% 
    group_by(id) %>% 
    summarise(x = unique(x), y = unique(y), mean = mean(mean, na.rm = T), freq = mean(freq, na.rm = T) ) 
) # eo ddf

# For zoo now
a.zoo <- data.frame(d.zoo %>% 
    group_by(id) %>% 
    summarise(x = unique(x), y = unique(y), mean = mean(mean, na.rm = T), freq = mean(freq, na.rm = T) ) 
) # eo ddf
# 
dim(a.phyto); dim(a.zoo)
summary(a.phyto)
summary(a.zoo)

setwd(WD)

# Map Average total
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean), data = a.phyto) +
    geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = mean), data = a.phyto) +
    scale_fill_gradient2(name = "Mean MESS\n(annual)", low = "#3288bd", high = "#d53e4f", mid = "white") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save map
ggsave(plot = map, filename = paste("map_mean_MESS_phyto_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean), data = a.zoo) +
    geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = mean), data = a.zoo) +
    scale_fill_gradient2(name = "Mean MESS\n(annual)", low = "#3288bd", high = "#d53e4f", mid = "white") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save map
ggsave(plot = map, filename = paste("map_mean_MESS_zoo_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)


### Map frequencies
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = freq), data = a.phyto) +
    geom_contour(colour = "grey50", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = freq), data = a.phyto) +
    scale_fill_distiller(name = "Annual frequency\n(MESS<0)", palette = "YlOrRd", direction = 1) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save map
ggsave(plot = map, filename = paste("map_mean_freq_phyto_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = freq), data = a.zoo) +
    geom_contour(colour = "grey50", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = freq), data = a.zoo) +
    scale_fill_distiller(name = "Annual frequency\n(MESS<0)", palette = "YlOrRd", direction = 1) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Save map
ggsave(plot = map, filename = paste("map_mean_freq_zoo_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)


### Very nice, and plot distribution of variables contributing to MESS < 0 (from mess.phyto & mess.zoo)
plot <- ggplot(data = mess.phyto[which(mess.phyto$TOTAL < 0),], aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
			scale_fill_brewer(name = "Variable", palette = "Spectral") + 
			geom_density(position = "stack", alpha = 0.5) + ylab("Density") + xlab("") + 
			theme_classic()
ggsave(plot = plot, filename = paste("distrib_vars_MESS_phyto_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)
            
#
plot <- ggplot(data = mess.zoo[which(mess.zoo$TOTAL < 0),], aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
			scale_fill_brewer(name = "Variable", palette = "Spectral") + 
			geom_density(position = "stack", alpha = 0.5) + ylab("Density") + xlab("") + 
			theme_classic()
ggsave(plot = plot, filename = paste("distrib_vars_MESS_zoo_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)


### 03/09/19: Tally the counts per variable
summary(factor(mess.phyto[which(mess.phyto$TOTAL < 0),c("MoD")]))
summary(factor(mess.zoo[which(mess.zoo$TOTAL < 0),c("MoD")]))

# Examine %
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "logNO3"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# 0.519%
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "SST"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# 0.389%
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "logChl"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# 0.042%
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "dSST"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# 0.019%
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "logSiO2"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# 0.023%
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "Sistar"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# ~0
nrow(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "Nstar"),]) / nrow( mess.phyto[which(mess.phyto$TOTAL < 0),] )
# ~0


### Zoo
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "logNO3"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# 0.488%
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "SST"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# 0.244%
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "logChl"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# 0.105%
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "dSST"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# 0.021%
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "logSiO2"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# 0.064%
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "Sistar"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# ~0
nrow(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "dO2"),]) / nrow( mess.zoo[which(mess.zoo$TOTAL < 0),] )
# 0.073%

### Summarize latitudinal distrib of these extrapolation events
summary( abs(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "logNO3"),c("y")]) ) #  9.50   20.50  34.50  
summary( abs(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "SST"),c("y")]) ) #     4.50   9.50   14.50   
summary( abs(mess.phyto[which(mess.phyto$TOTAL < 0 & mess.phyto$MoD == "logChl"),c("y")]) ) #  17.50  25.50  29.50   

summary( abs(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "logNO3"),c("y")]) ) # 7.50  17.50  32.50
summary( abs(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "SST"),c("y")]) )    # 3.50  8.50   14.50
summary( abs(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "logChl"),c("y")]) ) # 16.50 23.50  29.50
summary( abs(mess.zoo[which(mess.zoo$TOTAL < 0 & mess.zoo$MoD == "dO2"),c("y")]) )    # 12.50  15.5  19.50

### Show this with boxplot
plot <- ggplot(data = mess.phyto[which(mess.phyto$TOTAL < 0),]) + 
    geom_boxplot(aes(x = factor(MoD), y = abs(y), fill = factor(MoD)), colour = "black") + 
    scale_fill_brewer(name = "Variable", palette = "Spectral") + 
    ylab("Absolute Latitude") + xlab("Variables") + theme_classic()
            
ggsave(plot = plot, filename = paste("boxplot_distrib_lat_MESS_phyto_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)

plot <- ggplot(data = mess.zoo[which(mess.zoo$TOTAL < 0),]) + 
    geom_boxplot(aes(x = factor(MoD), y = abs(y), fill = factor(MoD)), colour = "black") + 
    scale_fill_brewer(name = "Variable", palette = "Spectral") + 
    ylab("Absolute Latitude") + xlab("Variables") + theme_classic()
            
ggsave(plot = plot, filename = paste("boxplot_distrib_lat_MESS_zoo_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)





