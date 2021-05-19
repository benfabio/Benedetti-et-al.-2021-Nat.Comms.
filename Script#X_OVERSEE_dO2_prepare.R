
# ---------------------------------------------------------------------------------------

library("ncdf4")
library("raster")
library("rgdal")
library("sp")
library("tidyverse")
library("stringr")
library("reshape2")
library("viridis")

# ----------------------------------------------------------------------------------------

### Grey coastline
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
  	  	theme(
    		panel.background = element_rect(fill = "white"),  # background
    		legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70")
  		 )
) 
quartz()
ggplot() + coast + theme_bw()


### Load some of the ncdf layers as raster and 
ras <- raster::raster("woa13_all_o_monthly.nc", level = 24, band = 1) # band = Jan
plot(ras)

ras2 <- raster::raster("woa13_all_o_monthly.nc", level = 24, band = 7) # # band = Jul
#ras2
quartz()
plot(ras2)

### OK, now extract the dO2 climatologies @175m depth for each month, map with coastline and save as data.frame
bands <- c(1:12)
res <- lapply(bands, function(b) {
			ras <- raster::raster("woa13_all_o_monthly.nc", level = 24, band = b)
			ddf <- as.data.frame(ras, xy = TRUE)
			colnames(ddf) <- c("x","y","dO2")
			if(b == 1) {
				return(ddf)
			} else {
				return(data.frame(dO2 = ddf[,"dO2"]))
			} # eo if else loop
	
		} # eo fun
) # eo lapply
table <- do.call(cbind, res)
dim(table)
colnames(table)[c(3:14)] <- c("dO2_Jan","dO2_Feb","dO2_Mar","dO2_Apr","dO2_May","dO2_Jun","dO2_Jul","dO2_Aug",
								"dO2_Sep","dO2_Oct","dO2_Nov","dO2_Dec")

# Save as .txt
write.table(table, file = "dO2_stack_1d.txt", sep = ";")

# And now as raster stack
coordinates(table) <- ~ x + y
gridded(table) <- TRUE
# coerce to raster
raster <- stack(table)
raster
plot(raster)
# Save
save(raster, file = "dO2_stack_1d.Rdata")

# And plot each monthly clim and arrange in panel 
m_table <- melt(table, id.vars = c("x","y"))
colnames(m_table) <- c("x","y","var","value")
plot <- ggplot(data = m_table) + geom_raster(aes(x = x, y = y, fill = value) ) + coast + scale_fill_viridis() + facet_wrap(~ var, scales = "free") 
ggsave(plot = plot, filename = paste("map_","dO2",".pdf", sep = ""), dpi = 300, width = 35, height = 15)


# ----------------------------------------------------------------------------------------

### 27/09/19: Examine seasonality of dissolved o2 at 175m depth

world2 <- map_data("world2") # for coastline

### Load monthly do2 clims
clim <- get(load("dO2_stack_1d.Rdata"))
# class(clim) # is raster stack
clim <- as.data.frame(clim, xy = T)
str(clim)
### Add cell ids: 1° and 4°
clim$id1 <- paste(clim$x, clim$y, sep = "_")
# To round to the nearest integer
mround <- function(x,base){
        base*round(x/base)
} # eo fun
# mround(163,5) 
head(mround(clim$x,4)) # OK
clim$id4 <- paste(mround(clim$x,4), mround(clim$y,4), sep = "_")
clim$x4 <- mround(clim$x, 4)
clim$y4 <- mround(clim$y, 4)

### Using melt and dplyr, compute min, max and range across the year within each 1d or 4d cells
require("dplyr")
require("reshape2")
mclim <- melt(clim, id.vars = c("x","y","x4","y4","id1","id4"))
head(mclim); summary(mclim)
ddf1 <- data.frame(na.omit(mclim) %>% group_by(id1) %>% summarize(x = unique(x), y = unique(y), min = min(value,na.rm=T), max = max(value,na.rm=T) ) ) # eo ddf
summary(ddf1)
### 4d
ddf4 <- data.frame(na.omit(mclim) %>% group_by(id4) %>% summarize(x = unique(x4), y = unique(y4), min = min(value,na.rm=T), max = max(value,na.rm=T) ) ) # eo ddf
summary(ddf4)

### Compute annual range 
head(ddf1); head(ddf4)
ddf1$range <- (ddf1$max) - (ddf1$min)
ddf4$range <- (ddf4$max) - (ddf4$min)

### Map results
# Rotate x 
ddf1$x2 <- ddf1$x 
ddf1[ddf1$x < 0 ,"x2"] <- (ddf1[ddf1$x < 0 ,"x"]) + 360
ddf4$x2 <- ddf4$x 
ddf4[ddf4$x < 0 ,"x2"] <- (ddf4[ddf4$x < 0 ,"x"]) + 360

#quartz()
m1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = range), data = ddf1) +
	scale_fill_viridis(name = "dO2\nannual range") +
    #geom_contour(colour = "grey75", binwidth = 1, size = 0.25, aes(x = x2, y = y, z = range), data = ddf1) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + ggtitle("Baseline")
#
#quartz()
m2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = range), data = ddf4) +
	scale_fill_viridis(name = "dO2\nannual range") +
    geom_contour(colour = "grey75", binwidth = 1, size = 0.25, aes(x = x2, y = y, z = range), data = ddf4) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + ggtitle("Baseline")
#
setwd("/Users/fabiobenedetti/Desktop/")
ggsave(plot = m1, filename = "map_dO2_annual_range_1d.pdf", width = 7, height = 3, dpi = 300)
ggsave(plot = m2, filename = "map_dO2_annual_range_4d.pdf", width = 7, height = 3, dpi = 300)


### Zonal plot showing the monthly variations with boxplots/violin plots, from mclim
head(mclim)
# Order months per number
mclim$month <- NA
mclim[mclim$variable == "dO2_Jan","month"] <- 1
mclim[mclim$variable == "dO2_Feb","month"] <- 2
mclim[mclim$variable == "dO2_Mar","month"] <- 3
mclim[mclim$variable == "dO2_Apr","month"] <- 4
mclim[mclim$variable == "dO2_May","month"] <- 5
mclim[mclim$variable == "dO2_Jun","month"] <- 6
mclim[mclim$variable == "dO2_Jul","month"] <- 7
mclim[mclim$variable == "dO2_Aug","month"] <- 8
mclim[mclim$variable == "dO2_Sep","month"] <- 9
mclim[mclim$variable == "dO2_Oct","month"] <- 10
mclim[mclim$variable == "dO2_Nov","month"] <- 11
mclim[mclim$variable == "dO2_Dec","month"] <- 12

#quartz()
plot <- ggplot(mclim, aes(x=factor(month), y=value)) + geom_violin(colour = "black", fill = "grey65") + 
    geom_boxplot(colour = "black", fill = "white", width = .15) + xlab("Months") + ylab("dO2 at 175m\n(ml/l)") +  
    theme_classic() + ggtitle("Monthly variations in dO2 - Baseline")
    
ggsave(plot = plot, filename = "plot_dO2_monthly_ranges_1d.pdf", width = 7, height = 5, dpi = 300)

### And examine future monthly variations in O2
setwd("/Users/fabiobenedetti/Desktop/~work/PostDocs/ETHZ/OVERSEE/data/GFDL-ESM2SM/future_monthly_clims/diff/rcp85/")

res <- lapply(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), function(m) {
			d <- read.table(paste("clim_2100-2071_rcp85_diff_",m,"_GFDL-ESM2M_24_07_19_v3.txt",sep = ""), sep = "\t", h = T)
			if(m == "Jan") {
			    return(d[,c("id","x","y","dO2")])
			} else {
			    return( data.frame(dO2 = d[,c("dO2")]) )
			}
		} # eo fun
) # eo lapply
fut.clims <- do.call(cbind, res)
rm(res)
summary(fut.clims)
colnames(fut.clims)[c(4:15)] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# melt and compute range & stuff
mclim.fut <- melt(fut.clims, id.vars = c("x","y","id"))
head(mclim.fut); summary(mclim.fut)
ddf.fut <- data.frame(na.omit(mclim.fut) %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), min = min(value,na.rm=T), max = max(value,na.rm=T) ) ) # eo ddf
summary(ddf.fut)
ddf.fut$range <- (ddf.fut$max) - (ddf.fut$min)

# Rotate x 
ddf.fut$x2 <- ddf.fut$x 
ddf.fut[ddf.fut$x < 0 ,"x2"] <- (ddf.fut[ddf.fut$x < 0 ,"x"]) + 360

# map
#quartz()
map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = range), data = ddf.fut) +
	scale_fill_viridis(name = "dO2\nannual range") +
    #geom_contour(colour = "grey75", binwidth = 1, size = 0.25, aes(x = x2, y = y, z = range), data = ddf1) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + ggtitle("2100 - RCP8.5")
        
ggsave(plot = map, filename = "map_dO2_annual_range_1d_rcp85.pdf", width = 7, height = 3, dpi = 300)

### And examine future monthly variations
head(mclim.fut)
# Order months per number
mclim.fut$month <- NA
mclim.fut[mclim.fut$variable == "Jan","month"] <- 1
mclim.fut[mclim.fut$variable == "Feb","month"] <- 2
mclim.fut[mclim.fut$variable == "Mar","month"] <- 3
mclim.fut[mclim.fut$variable == "Apr","month"] <- 4
mclim.fut[mclim.fut$variable == "May","month"] <- 5
mclim.fut[mclim.fut$variable == "Jun","month"] <- 6
mclim.fut[mclim.fut$variable == "Jul","month"] <- 7
mclim.fut[mclim.fut$variable == "Aug","month"] <- 8
mclim.fut[mclim.fut$variable == "Sep","month"] <- 9
mclim.fut[mclim.fut$variable == "Oct","month"] <- 10
mclim.fut[mclim.fut$variable == "Nov","month"] <- 11
mclim.fut[mclim.fut$variable == "Dec","month"] <- 12

#quartz()
plot <- ggplot(mclim.fut, aes(x=factor(month), y=value)) + geom_violin(colour = "black", fill = "grey65") + 
    geom_boxplot(colour = "black", fill = "white", width = .15) + xlab("Months") + ylab("dO2 at 175m\n(ml/l)") +  
    theme_classic() + ggtitle("Monthly variations in dO2 - 2100 RCP85")
    
ggsave(plot = plot, filename = "plot_dO2_monthly_ranges_1d_rcp85.pdf", width = 7, height = 5, dpi = 300)
### Ok, save the plots you want and summarize the in a .ppt


### 27/09/19: Examine monthly vaitaions in doO2 from the model outputs purely (no delta)
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/")
dir() # "clims_future_2071-2100_o2_rcp85_GFDL-ESM2M.txt" 
model <- read.table("clims_future_2071-2100_o2_rcp85_GFDL-ESM2M.txt", h = T, sep = "\t")
dim(model); str(model); summary(na.omit(model)) # not the right unit yet
# Melt
m.model <- melt(model, id.vars = c("id","x","y"))
#head(m.model)
colnames(m.model)[c(4:5)] <- c("month","value")
# summary(m.model)

# Multiply by 10e5 (convert M to µM)
m.model$value <- (m.model$value)*1.029
m.model$value <- (m.model$value)/44.661e-6
summary(m.model)
# ok
m.model$month2 <- NA
m.model[m.model$month == "Jan","month2"] <- 1
m.model[m.model$month == "Feb","month2"] <- 2
m.model[m.model$month == "Mar","month2"] <- 3
m.model[m.model$month == "Apr","month2"] <- 4
m.model[m.model$month == "May","month2"] <- 5
m.model[m.model$month == "Jun","month2"] <- 6
m.model[m.model$month == "Jul","month2"] <- 7
m.model[m.model$month == "Aug","month2"] <- 8
m.model[m.model$month == "Sep","month2"] <- 9
m.model[m.model$month == "Oct","month2"] <- 10
m.model[m.model$month == "Nov","month2"] <- 11
m.model[m.model$month == "Dec","month2"] <- 12

plot <- ggplot(m.model, aes(x=factor(month2), y=value)) + geom_violin(colour = "black", fill = "grey65") + 
    geom_boxplot(colour = "black", fill = "white", width = .15) + xlab("Months") + ylab("dO2 at 175m\n(ml/l)") +  
    scale_y_continuous(limits = c(0,10)) + theme_classic() + ggtitle("Monthly variations in dO2 - 2100 RCP85 (no delta)")

setwd("/net/kryo/work/fabioben/OVERSEE/data/")    
ggsave(plot = plot, filename = "plot_dO2_monthly_ranges_1d_rcp85_nodelta.pdf", width = 7, height = 5, dpi = 300)

### And map annual range 
ddf.mod <- data.frame(na.omit(m.model) %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), min = min(value,na.rm=T), max = max(value,na.rm=T) ) ) # eo ddf
ddf.mod$range <- (ddf.mod$max) - (ddf.mod$min)
summary(ddf.mod)

# Map
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = range), data = ddf.mod) +
	scale_fill_viridis(name = "dO2\nannual range") +
    geom_contour(colour = "grey75", binwidth = 1, size = 0.25, aes(x = x, y = y, z = range), data = ddf.mod) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + ggtitle("2100 - RCP8.5 (model only)")
        
ggsave(plot = map, filename = "map_dO2_annual_range_1d_rcp85_nodelta.pdf", width = 7, height = 3, dpi = 300)


### And same with baseline (model only)
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/")
dir() # "clims_baseline_1971-2000_o2_rcp26_GFDL-ESM2M.txt" 
model.base <- read.table("clims_baseline_1971-2000_o2_rcp85_GFDL-ESM2M.txt", h = T, sep = "\t")
m.model.base <- melt(model.base, id.vars = c("id","x","y"))
# head(m.model.base)
colnames(m.model.base)[c(4:5)] <- c("month","value")
# summary(m.model.base)

# Multiply by 10e5 (convert M to µM)
m.model.base$value <- (m.model.base$value)*1.029
m.model.base$value <- (m.model.base$value)/44.661e-6
summary(m.model.base)
# ok
m.model.base$month2 <- NA
m.model.base[m.model.base$month == "Jan","month2"] <- 1
m.model.base[m.model.base$month == "Feb","month2"] <- 2
m.model.base[m.model.base$month == "Mar","month2"] <- 3
m.model.base[m.model.base$month == "Apr","month2"] <- 4
m.model.base[m.model.base$month == "May","month2"] <- 5
m.model.base[m.model.base$month == "Jun","month2"] <- 6
m.model.base[m.model.base$month == "Jul","month2"] <- 7
m.model.base[m.model.base$month == "Aug","month2"] <- 8
m.model.base[m.model.base$month == "Sep","month2"] <- 9
m.model.base[m.model.base$month == "Oct","month2"] <- 10
m.model.base[m.model.base$month == "Nov","month2"] <- 11
m.model.base[m.model.base$month == "Dec","month2"] <- 12

plot <- ggplot(m.model.base, aes(x=factor(month2), y=value)) + geom_violin(colour = "black", fill = "grey65") + 
    geom_boxplot(colour = "black", fill = "white", width = .15) + xlab("Months") + ylab("dO2 at 175m\n(ml/l)") +  
    scale_y_continuous(limits = c(0,10)) + theme_classic() + ggtitle("Monthly variations in dO2 - Baseline (model only)")

setwd("/net/kryo/work/fabioben/OVERSEE/data/")    
ggsave(plot = plot, filename = "plot_dO2_monthly_ranges_baseline_GFDL.pdf", width = 7, height = 5, dpi = 300)


### And map annual range 
ddf.mod <- data.frame(na.omit(m.model.base) %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), min = min(value,na.rm=T), max = max(value,na.rm=T) ) ) # eo ddf
ddf.mod$range <- (ddf.mod$max) - (ddf.mod$min)
summary(ddf.mod)

# Map
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = range), data = ddf.mod) +
	scale_fill_viridis(name = "dO2\nannual range") +
    geom_contour(colour = "grey75", binwidth = 1, size = 0.25, aes(x = x, y = y, z = range), data = ddf.mod) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + ggtitle("Baseline (model only)")
        
ggsave(plot = map, filename = "map_dO2_annual_range_baseline_GFDL.pdf", width = 7, height = 3, dpi = 300)



### 02/09/19: Same as you did for PAR/RSDS, combine the model-based (no delta) baseline and future climatologies of dO2
### and examine change in seasonality
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM/")
dir() # "clims_future_2071-2100_o2_rcp85_GFDL-ESM2M.txt" 
model.base <- read.table("clims_baseline_1971-2000_dO2_rcp85_GFDL-ESM2M.txt", h = T, sep = "\t")
dim(model.base); str(model.base); summary(na.omit(model.base)) # not the right unit yet

model.base <- melt(model.base, id.vars = c("id","x","y"))
colnames(model.base)[c(4:5)] <- c("month","value")

# summary(model.base)
# ok
model.base$month2 <- NA
model.base[model.base$month == "Jan","month2"] <- 1
model.base[model.base$month == "Feb","month2"] <- 2
model.base[model.base$month == "Mar","month2"] <- 3
model.base[model.base$month == "Apr","month2"] <- 4
model.base[model.base$month == "May","month2"] <- 5
model.base[model.base$month == "Jun","month2"] <- 6
model.base[model.base$month == "Jul","month2"] <- 7
model.base[model.base$month == "Aug","month2"] <- 8
model.base[model.base$month == "Sep","month2"] <- 9
model.base[model.base$month == "Oct","month2"] <- 10
model.base[model.base$month == "Nov","month2"] <- 11
model.base[model.base$month == "Dec","month2"] <- 12

### And same for future
model.fut <- read.table("clims_future_2071-2100_o2_rcp85_GFDL-ESM2M.txt", h = T, sep = "\t")
dim(model.fut); str(model.fut); summary(na.omit(model.fut)) # not the right unit yet

model.fut <- melt(model.fut, id.vars = c("id","x","y"))
colnames(model.fut)[c(4:5)] <- c("month","value")
model.fut$value <- (model.fut$value)*1.029
model.fut$value <- (model.fut$value)/44.661e-6
# summary(model.fut)
# ok
model.fut$month2 <- NA
model.fut[model.fut$month == "Jan","month2"] <- 1
model.fut[model.fut$month == "Feb","month2"] <- 2
model.fut[model.fut$month == "Mar","month2"] <- 3
model.fut[model.fut$month == "Apr","month2"] <- 4
model.fut[model.fut$month == "May","month2"] <- 5
model.fut[model.fut$month == "Jun","month2"] <- 6
model.fut[model.fut$month == "Jul","month2"] <- 7
model.fut[model.fut$month == "Aug","month2"] <- 8
model.fut[model.fut$month == "Sep","month2"] <- 9
model.fut[model.fut$month == "Oct","month2"] <- 10
model.fut[model.fut$month == "Nov","month2"] <- 11
model.fut[model.fut$month == "Dec","month2"] <- 12

# Rbind after informing time period and plot
model.fut$period <- "Future"
model.base$period <- "Baseline"

ddf <- rbind(model.base, model.fut)
head(ddf)

### Plot monthly distributions of RSDS
setwd("/net/kryo/work/fabioben/OVERSEE/data/")
plot <- ggplot(ddf, aes(x= factor(month2), y=value, fill = factor(period))) + geom_violin(colour = "black", position = "dodge", alpha = 0.5) + 
    geom_boxplot(colour = "black", width = .25, position = "dodge", notch = T) + 
    scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
    xlab("Months") + ylab("O2 at 175m (µl/l)") + theme_classic() + ggtitle("Baseline and future monthly variations in dO2 - GFDL-TOPAZ")

ggsave(plot = plot, filename = "plot_o2_monthly_ranges_GFDL-TOPAZ.pdf", width = 12, height = 6, dpi = 300)
# Nice, do the same for a particular lat band

plot <- ggplot(ddf[ddf$y < 30 & ddf$y > -30,], aes(x= factor(month2), y=value, fill = factor(period))) + geom_violin(colour = "black", position = "dodge", alpha = 0.5) + 
    geom_boxplot(colour = "black", width = .25, position = "dodge", notch = T) + 
    scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
    xlab("Months") + ylab("O2 at 175m (µl/l)") + theme_classic() + ggtitle("Baseline and future monthly variations in dO2 - GFDL-TOPAZ")

ggsave(plot = plot, filename = "plot_o2_monthly_ranges_GFDL-TOPAZ_tropics.pdf", width = 12, height = 6, dpi = 300)

plot <- ggplot(ddf[ddf$y < 60 & ddf$y > 35,], aes(x= factor(month2), y=value, fill = factor(period))) + geom_violin(colour = "black", position = "dodge", alpha = 0.5) + 
    geom_boxplot(colour = "black", width = .25, position = "dodge", notch = T) + 
    scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
    xlab("Months") + ylab("O2 at 175m (µl/l)") + theme_classic() + ggtitle("Baseline and future monthly variations in dO2 - GFDL-TOPAZ")

ggsave(plot = plot, filename = "plot_o2_monthly_ranges_GFDL-TOPAZ_temperate.pdf", width = 12, height = 6, dpi = 300)

plot <- ggplot(ddf[ddf$y < 90 & ddf$y > 60,], aes(x= factor(month2), y=value, fill = factor(period))) + geom_violin(colour = "black", position = "dodge", alpha = 0.5) + 
    geom_boxplot(colour = "black", width = .25, position = "dodge", notch = T) + 
    scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
    xlab("Months") + ylab("O2 at 175m (µl/l)") + theme_classic() + ggtitle("Baseline and future monthly variations in dO2 - GFDL-TOPAZ")

ggsave(plot = plot, filename = "plot_o2_monthly_ranges_GFDL-TOPAZ_poles.pdf", width = 12, height = 6, dpi = 300)

plot <- ggplot(ddf[ddf$y < -10 & ddf$y > -20,], aes(x= factor(month2), y=value, fill = factor(period))) + geom_violin(colour = "black", position = "dodge", alpha = 0.5) + 
    geom_boxplot(colour = "black", width = .25, position = "dodge", notch = T) + 
    scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
    xlab("Months") + ylab("O2 at 175m (µl/l)") + theme_classic() + ggtitle("Baseline and future monthly variations in dO2 - GFDL-TOPAZ")

ggsave(plot = plot, filename = "plot_o2_monthly_ranges_GFDL-TOPAZ_ETP.pdf", width = 12, height = 6, dpi = 300)



