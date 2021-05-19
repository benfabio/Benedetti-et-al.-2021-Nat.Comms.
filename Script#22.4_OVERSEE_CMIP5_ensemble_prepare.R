
##### 11/11/2019 - ETHZ - Fabio Benedetti © UP Group+ IBP+ ETH Zürich
##### Script for : 
#	- extracting the CMIP5 models outputs from the netCDF files prepared by C. Laufkötter
#   - preparing the baseline (to be defined) and future (2081-2100) climatologies from the monthly climatologies
#   - derive the difference (delta) to apply to the in situ climatologies
#   - plot and map outputs along the way 
#   - do the same for O2 and O2 from the various sorts of netCDF...
 
### Last update: 18/11/2019

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("ncdf4")
library("raster")
library("reshape2")
library("scales")
library("maps")
library("RColorBrewer")
library("viridis")

# Coastline
world2 <- map_data("world2")

# WD
WD <- getwd()

# Vector of months names
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# --------------------------------------------------------------------------------------------------------------------------------

##### 1°) Start with ALL MODELS by LOTTE (SST, Si, NO3, PO4, CHL (nano+diat))
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/oxygen_files_separate_esm")

### 1.A) Extract and format SST monthly outputs 
# sst1.nc
nc <- nc_open("sst1.nc")
nc # units: days since 2006-01-01 00:00:00
names(nc$var) # DATA_pisces + DATA_CNRM # ok
# And with raster ? 
ras <- raster::stack("sst1.nc", varname = "DATA_pisces") # By default --> DATA_pisces
ras
# dimensions  : 180, 360, 64800, 1068  (nrow, ncol, ncell, nlayers)
# nlayers ? = 2880/12 = 89 years --> one layer = one month
r <- as.data.frame(ras, xy = T)
dim(r) # 64800  1070
head(r[,c(1:10)])
summary(r[,c(1:10)]) # x = 0-360°, y = -90°/+90°
colnames(r) # 3 to 1068 = 1068 months and years
# Ok, let's go :-) Extract and make data.frames out of the .nc files

### Choose ESM
ras <- raster::stack("sst1.nc", varname = "DATA_pisces")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract SST data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
# dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
# cols[c(1:240)] # for baseline
# cols[c(829:length(cols))] # for future conditions
sst_baseline <- r[,c(1,2,3:242)]
sst_fut <- r[,c(1,2,831:length(r))]
# dim(sst_baseline); dim(sst_fut)
head(sst_baseline); head(sst_fut)

# Add an id 
sst_baseline$id <- factor(paste(sst_baseline$x, sst_baseline$y, sep = "_"))
sst_fut$id <- factor(paste(sst_fut$x, sst_fut$y, sep = "_"))
head(sst_baseline$id); head(sst_fut$id)

# Melt, and no need to convert to SST !
m_sst_base <- melt(sst_baseline, id.vars = c("id","x","y"))
m_sst_fut <- melt(sst_fut, id.vars = c("id","x","y"))
colnames(m_sst_base)[c(4,5)] <- c("month","SST")
colnames(m_sst_fut)[c(4,5)] <- c("month","SST")
# dim(m_sst_base); dim(m_sst_fut)
# head(m_sst_base); head(m_sst_fut)
summary(m_sst_base)
summary(m_sst_fut) # +4°C in max temperatures, ok guys
# head(m_sst_base)
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_sst_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_sst_base$month), "_") ))[,1]
m_sst_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_sst_fut$month), "_") ))[,1]
# head(m_sst_base); head(m_sst_fut)

# Clean some stuff
rm(sst_fut, sst_baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_sst_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SST = mean(SST) ) )
clims_fut <- data.frame( m_sst_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SST = mean(SST) ) )
# dim(clims_base); dim(clims_fut)
# head(clims_base); head(clims_fut)
### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "SST")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "SST")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_sst_rcp85_CNRM-PISCES_2031-2100.Rdata")
save(clims_base, file = "clims_mon_sst_rcp85_IPSL-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_sst_rcp85_IPSL-PISCES_2100-2081.Rdata")
#
# # maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta SST (°C)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_sst_rcp85_CNRM-PISCES_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop



# sst2.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Monthly forecast files by Lotte (ALL MODELS)")
nc <- nc_open("sst2.nc")
nc # units: days since 2008-01-01 00:00:00
names(nc$var) # DATA_GFDL + DATA_PELAGOS + DATA_nemuro # ok, skip PELAGOS though

### Choose ESL
ras <- raster::stack("sst2.nc", varname = "DATA_nemuro")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract SST data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
sst_baseline <- r[,c(1,2,3:242)]
sst_fut <- r[,c(1,2,831:length(r))]
# dim(sst_baseline); dim(sst_fut)
head(sst_baseline); head(sst_fut)

# Add an id 
sst_baseline$id <- factor(paste(sst_baseline$x, sst_baseline$y, sep = "_"))
sst_fut$id <- factor(paste(sst_fut$x, sst_fut$y, sep = "_"))
head(sst_baseline$id); head(sst_fut$id)

# Melt, and no need to convert to SST !
m_sst_base <- melt(sst_baseline, id.vars = c("id","x","y"))
m_sst_fut <- melt(sst_fut, id.vars = c("id","x","y"))
colnames(m_sst_base)[c(4,5)] <- c("month","SST")
colnames(m_sst_fut)[c(4,5)] <- c("month","SST")
summary(m_sst_base)
summary(m_sst_fut) # +4°C in max temperatures, ok guys
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_sst_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_sst_base$month), "_") ))[,1]
m_sst_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_sst_fut$month), "_") ))[,1]
head(m_sst_base); head(m_sst_fut)

# Clean some stuff
rm(sst_fut, sst_baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_sst_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SST = mean(SST) ) )
clims_fut <- data.frame( m_sst_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SST = mean(SST) ) )
# dim(clims_base); dim(clims_fut)
# head(clims_base); head(clims_fut)
### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "SST")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "SST")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
save(deltas, file = "clims_mon_diff_sst_rcp85_MRI-NEMURO_2031-2100.Rdata")
# And monthly means
save(clims_base, file = "clims_mon_sst_rcp85_MRI-NEMURO_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_sst_rcp85_MRI-NEMURO_2100-2081.Rdata")

# # maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta SST (°C)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_sst_rcp85_MRI-NEMURO_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


# sst3.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Monthly forecast files by Lotte (ALL MODELS)")
nc <- nc_open("sst3.nc")
nc # units: days since 2008-01-01 00:00:00
names(nc$var) # DATA_Hadgem + DATA_BEC # ok, skip DATA_Hadgem though
### Choose ESL
ras <- raster::stack("sst3.nc", varname = "DATA_BEC")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract SST data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
sst_baseline <- r[,c(1,2,3:242)]
sst_fut <- r[,c(1,2,831:length(r))]
# dim(sst_baseline); dim(sst_fut)
head(sst_baseline); head(sst_fut)

# Add an id 
sst_baseline$id <- factor(paste(sst_baseline$x, sst_baseline$y, sep = "_"))
sst_fut$id <- factor(paste(sst_fut$x, sst_fut$y, sep = "_"))
head(sst_baseline$id); head(sst_fut$id)

# Melt, and no need to convert to SST !
m_sst_base <- melt(sst_baseline, id.vars = c("id","x","y"))
m_sst_fut <- melt(sst_fut, id.vars = c("id","x","y"))
colnames(m_sst_base)[c(4,5)] <- c("month","SST")
colnames(m_sst_fut)[c(4,5)] <- c("month","SST")
summary(m_sst_base)
summary(m_sst_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_sst_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_sst_base$month), "_") ))[,1]
m_sst_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_sst_fut$month), "_") ))[,1]
head(m_sst_base); head(m_sst_fut)

# Clean some stuff
rm(sst_fut, sst_baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_sst_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SST = mean(SST) ) )
clims_fut <- data.frame( m_sst_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SST = mean(SST) ) )
### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "SST")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "SST")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
save(deltas, file = "clims_mon_diff_sst_rcp85_CESM-BEC_2031-2100.Rdata")
save(clims_base, file = "clims_mon_sst_rcp85_CESM-BEC_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_sst_rcp85_CESM-BEC_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta SST (°C)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_sst_rcp85_CESM-BEC_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop
#


### 1.B) Extract and format NO3 monthly outputs ----------------------------------------------------------------------
# no31.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Monthly forecast files by Lotte (ALL MODELS)")
nc <- nc_open("no31.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_pisces + DATA_CNRM # ok, use both

### Choose ESM
ras <- raster::stack("no31.nc", varname = "DATA_pisces")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract SST data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","NO3")
colnames(m_fut)[c(4,5)] <- c("month","NO3")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), NO3 = mean(NO3) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), NO3 = mean(NO3) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$NO3)*1000)); summary(log1p((clims_fut$NO3)*1000))
clims_base$NO3 <- (clims_base$NO3)*1000
clims_fut$NO3 <- (clims_fut$NO3)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "NO3")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "NO3")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
save(deltas, file = "clims_mon_diff_no3_rcp85_CNRM-PISCES_2031-2100.Rdata")
save(clims_base, file = "clims_mon_no3_rcp85_IPSL-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_no3_rcp85_IPSL-PISCES_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta NO3 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_no3_rcp85_CNRM-PISCES_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop



# no32.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("no32.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_GFDL + DATA_nemuro # ok, use both
### Choose ESM
ras <- raster::stack("no32.nc", varname = "DATA_nemuro")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract SST data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","NO3")
colnames(m_fut)[c(4,5)] <- c("month","NO3")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), NO3 = mean(NO3) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), NO3 = mean(NO3) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$NO3)*1000)); summary(log1p((clims_fut$NO3)*1000))
clims_base$NO3 <- (clims_base$NO3)*1000
clims_fut$NO3 <- (clims_fut$NO3)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "NO3")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "NO3")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
save(clims_base, file = "clims_mon_no3_rcp85_MRI-NEMURO_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_no3_rcp85_MRI-NEMURO_2100-2081.Rdata")
save(deltas, file = "clims_mon_diff_no3_rcp85_MRI-NEMURO_2031-2100.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta NO3 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_no3_rcp85_MRI-NEMURO_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


# no33.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("no33.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_Hadgem + DATA_BEC # ok, use both
### Choose ESM
ras <- raster::stack("no33.nc", varname = "DATA_BEC")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract SST data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","NO3")
colnames(m_fut)[c(4,5)] <- c("month","NO3")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), NO3 = mean(NO3) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), NO3 = mean(NO3) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$NO3)*1000)); summary(log1p((clims_fut$NO3)*1000))
clims_base$NO3 <- (clims_base$NO3)*1000
clims_fut$NO3 <- (clims_fut$NO3)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "NO3")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "NO3")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
save(clims_base, file = "clims_mon_no3_rcp85_CESM-BEC_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_no3_rcp85_CESM-BEC_2100-2081.Rdata")
#save(deltas, file = "clims_mon_diff_no3_rcp85_CESM-BEC_2031-2100.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta NO3 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_no3_rcp85_CESM-BEC_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop



### 1.C) Extract and format PO4 monthly outputs  ----------------------------------------------------------------------
# po41.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("po41.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_pisces + DATA_CNRM # ok, use both
### Choose ESM
ras <- raster::stack("po41.nc", varname = "DATA_CNRM")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","PO4")
colnames(m_fut)[c(4,5)] <- c("month","PO4")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), PO4 = mean(PO4) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), PO4 = mean(PO4) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$PO4)*1000)); summary(log1p((clims_fut$PO4)*1000))
clims_base$PO4 <- (clims_base$PO4)*1000
clims_fut$PO4 <- (clims_fut$PO4)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "PO4")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "PO4")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_po4_rcp85_CNRM-PISCES_2031-2100.Rdata")
save(clims_base, file = "clims_mon_po4_rcp85_CNRM-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_po4_rcp85_CNRM-PISCES_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta PO4 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_po4_rcp85_CNRM-PISCES_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


# po42.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("po42.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_pisces + DATA_CNRM # ok, use both
### Choose ESM
ras <- raster::stack("po42.nc", varname = "DATA_GFDL")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","PO4")
colnames(m_fut)[c(4,5)] <- c("month","PO4")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), PO4 = mean(PO4) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), PO4 = mean(PO4) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$PO4)*1000)); summary(log1p((clims_fut$PO4)*1000))
clims_base$PO4 <- (clims_base$PO4)*1000
clims_fut$PO4 <- (clims_fut$PO4)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "PO4")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "PO4")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_po4_rcp85_GFDL-TOPAZ_2031-2100.Rdata")
save(clims_base, file = "clims_mon_po4_rcp85_GFDL-TOPAZ_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_po4_rcp85_GFDL-TOPAZ_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta PO4 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_po4_rcp85_GFDL-TOPAZ_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


# po43.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("po43.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_Hadgem + DATA_BEC # ok, use DATA_BEC
### Choose ESM
ras <- raster::stack("po43.nc", varname = "DATA_BEC")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","PO4")
colnames(m_fut)[c(4,5)] <- c("month","PO4")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), PO4 = mean(PO4) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), PO4 = mean(PO4) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$PO4)*1000)); summary(log1p((clims_fut$PO4)*1000))
clims_base$PO4 <- (clims_base$PO4)*1000
clims_fut$PO4 <- (clims_fut$PO4)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "PO4")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "PO4")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_po4_rcp85_CESM-BEC_2031-2100.Rdata")
save(clims_base, file = "clims_mon_po4_rcp85_CESM-BEC_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_po4_rcp85_CESM-BEC_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta PO4 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_po4_rcp85_CESM-BEC_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


### 12/11/19: For MRI-NEMURO, compute monthly climatologies and diff of PO4 based on NO3 and the Redfield C:N:P ratio
### (106:16:1) --> divide average NO3 by 16 to get average PO4
clims_base_no3 <- get(load("clims_mon_no3_rcp85_MRI-NEMURO_2031-2012.Rdata"))
clims_fut_no3 <- get(load("clims_mon_no3_rcp85_MRI-NEMURO_2100-2081.Rdata"))
dim(clims_base_no3); dim(clims_fut_no3)
# Convert to PO4
clims_base_po4 <- (clims_base_no3[,c(4:15)])/16
clims_fut_po4 <- (clims_fut_no3[,c(4:15)])/16
summary(clims_base_po4)
summary(clims_fut_po4)
summary(log1p(clims_base_po4$Jan))
summary(log1p(clims_fut_po4$Oct))

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base_po4$id, x = clims_base_po4$x, y = clims_base_po4$y, 
	Jan = (clims_fut_po4$Jan) - (clims_base_po4$Jan), Feb = (clims_fut_po4$Feb) - (clims_base_po4$Feb), 
    Mar = (clims_fut_po4$Mar) - (clims_base_po4$Mar), 
	Apr = (clims_fut_po4$Apr) - (clims_base_po4$Apr), May = (clims_fut_po4$May) - (clims_base_po4$May), 
    Jun = (clims_fut_po4$Jun) - (clims_base_po4$Jun),
	Jul = (clims_fut_po4$Jul) - (clims_base_po4$Jul), Aug = (clims_fut_po4$Aug) - (clims_base_po4$Aug), 
    Sep = (clims_fut_po4$Sep) - (clims_base_po4$Sep), 
	Oct = (clims_fut_po4$Oct) - (clims_base_po4$Oct), Nov = (clims_fut_po4$Nov) - (clims_base_po4$Nov), 
    Dec = (clims_fut_po4$Dec) - (clims_base_po4$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
save(deltas, file = "clims_mon_diff_po4_rcp85_MRI-NEMURO_2031-2100.Rdata")
save(clims_base_po4, file = "clims_mon_po4_rcp85_MRI-NEMURO_2031-2012.Rdata")
save(clims_fut_po4, file = "clims_mon_po4_rcp85_MRI-NEMURO_2100-2081.Rdata")



### 1.D) Extract and format SiO2 monthly outputs  ----------------------------------------------------------------------
# si1.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("si1.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_pisces + DATA_CNRM # ok, use both
### Choose ESM
ras <- raster::stack("si1.nc", varname = "DATA_CNRM")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","SiO2")
colnames(m_fut)[c(4,5)] <- c("month","SiO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SiO2 = mean(SiO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SiO2 = mean(SiO2) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$SiO2)*1000)); summary(log1p((clims_fut$SiO2)*1000))
clims_base$SiO2 <- (clims_base$SiO2)*1000
clims_fut$SiO2 <- (clims_fut$SiO2)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "SiO2")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "SiO2")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_sio2_rcp85_CNRM-PISCES_2031-2100.Rdata")
save(clims_base, file = "clims_mon_sio2_rcp85_CNRM-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_sio2_rcp85_CNRM-PISCES_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta SiO2 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_sio2_rcp85_CNRM-PISCES_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


# si2.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("si2.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # "DATA_GFDL" + "DATA_nemuro" # ok, use both
### Choose ESM
ras <- raster::stack("si2.nc", varname = "DATA_nemuro")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","SiO2")
colnames(m_fut)[c(4,5)] <- c("month","SiO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SiO2 = mean(SiO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SiO2 = mean(SiO2) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$SiO2)*1000)); summary(log1p((clims_fut$SiO2)*1000))
clims_base$SiO2 <- (clims_base$SiO2)*1000
clims_fut$SiO2 <- (clims_fut$SiO2)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "SiO2")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "SiO2")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_sio2_rcp85_MRI-NEMURO_2031-2100.Rdata")
save(clims_base, file = "clims_mon_sio2_rcp85_MRI-NEMURO_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_sio2_rcp85_MRI-NEMURO_2100-2081.Rdata")

# # maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta SiO2 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_sio2_rcp85_MRI-NEMURO_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


# si3.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
nc <- nc_open("si3.nc")
nc  # units = mol m-3 --> will need to convert to µmol
names(nc$var) # DATA_BEC # ok, use both
### Choose ESM
ras <- raster::stack("si3.nc", varname = "DATA_BEC")
ras

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,3:242)]
fut <- r[,c(1,2,831:length(r))]
# dim(baseline); dim(fut)
head(baseline); head(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","SiO2")
colnames(m_fut)[c(4,5)] <- c("month","SiO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SiO2 = mean(SiO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), SiO2 = mean(SiO2) ) )
head(clims_base); head(clims_fut)
# Convert to µmol --> *1000
# summary(log1p((clims_base$SiO2)*1000)); summary(log1p((clims_fut$SiO2)*1000))
clims_base$SiO2 <- (clims_base$SiO2)*1000
clims_fut$SiO2 <- (clims_fut$SiO2)*1000

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "SiO2")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "SiO2")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
dir()
#save(deltas, file = "clims_mon_diff_sio2_rcp85_CESM-BEC_2031-2100.Rdata")
save(clims_base, file = "clims_mon_sio2_rcp85_CESM-BEC_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_sio2_rcp85_CESM-BEC_2100-2081.Rdata")

# maps
# for(momo in months) {
#
#         d <- deltas[,c("x","y",momo)]
#         colnames(d)[c(3)] <- "delta"
#
#         map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
#                  scale_fill_gradient2(name = "Delta SiO2 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#         ggsave(plot = map, filename = paste("map_delta_sio2_rcp85_CESM-BEC_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#         rm(d);gc()
# } # eo for loop


### 12/11/19: log1p() transform nutrients climatologies
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
ESM <- c("CNRM-PISCES","IPSL-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")
# m <- "MRI-NEMURO" # for testing for loop below
# log transform in a for loop 
for(m in ESM) {

        ### A) logNO3 ------------------------------------------------------------------------------------------------------------------
        message(paste("log transforming NO3 for ", m, sep = "")) 
        base <- get(load(paste("clims_mon_no3_rcp85_",m,"_2031-2012.Rdata", sep = "")))
        fut <- get(load(paste("clims_mon_no3_rcp85_",m,"_2100-2081.Rdata", sep = "")))
        # log transform
        # summary(log1p(base[,c(4:15)])); summary(log1p(fut[,c(4:15)]))
        base[,c(4:15)] <- log1p(base[,c(4:15)])
        fut[,c(4:15)] <- log1p(fut[,c(4:15)])
        # summary(base); summary(fut)
        
        # Compute delta between future and baseline
		delta <- data.frame(id = base$id, x = base$x, y = base$y, 
			    Jan = (fut$Jan) - (base$Jan), Feb = (fut$Feb) - (base$Feb), 
                Mar = (fut$Mar) - (base$Mar), Apr = (fut$Apr) - (base$Apr), 
                May = (fut$May) - (base$May), Jun = (fut$Jun) - (base$Jun),
			    Jul = (fut$Jul) - (base$Jul), Aug = (fut$Aug) - (base$Aug), 
                Sep = (fut$Sep) - (base$Sep), Oct = (fut$Oct) - (base$Oct), 
                Nov = (fut$Nov) - (base$Nov), Dec = (fut$Dec) - (base$Dec)
		) # eo ddf
        # summary(delta)
        
        # Map and save
        save(x = base, file = paste("clims_mon_logno3_rcp85_",m,"_2031-2012.Rdata", sep = "") )
        save(x = fut, file = paste("clims_mon_logno3_rcp85_",m,"_2100-2081.Rdata", sep = "") )
        save(x = delta, file = paste("clims_mon_diff_logno3_rcp85_",m,"_2031-2100.Rdata", sep = "") )
        
        ### Print maps
        months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
        # For momo in months
        for(momo in months) {
            
            d1 <- delta[,c("x","y",momo)]
            colnames(d1)[c(3)] <- "delta"
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d1) ) +
                scale_fill_gradient2(name = "Delta logNO3 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
            ggsave(plot = map, filename = paste("map_delta_logno3_rcp85_",m,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
                 
        } # eo for loop
        
        # ### B) logPO4 ------------------------------------------------------------------------------------------------------------------
#         message(paste("log transforming PO4 for ", m, sep = ""))
#         base <- get(load(paste("clims_mon_po4_rcp85_",m,"_2031-2012.Rdata", sep = "")))
#         fut <- get(load(paste("clims_mon_po4_rcp85_",m,"_2100-2081.Rdata", sep = "")))
#         # log transform
#         # summary(log1p(base[,c(4:15)])); summary(log1p(fut[,c(4:15)]))
#         base[,c(4:15)] <- log1p(base[,c(4:15)])
#         fut[,c(4:15)] <- log1p(fut[,c(4:15)])
#         # summary(base); summary(fut)
#
#         # Compute delta between future and baseline
#         delta <- data.frame(id = base$id, x = base$x, y = base$y,
#                 Jan = (fut$Jan) - (base$Jan), Feb = (fut$Feb) - (base$Feb),
#                 Mar = (fut$Mar) - (base$Mar), Apr = (fut$Apr) - (base$Apr),
#                 May = (fut$May) - (base$May), Jun = (fut$Jun) - (base$Jun),
#                 Jul = (fut$Jul) - (base$Jul), Aug = (fut$Aug) - (base$Aug),
#                 Sep = (fut$Sep) - (base$Sep), Oct = (fut$Oct) - (base$Oct),
#                 Nov = (fut$Nov) - (base$Nov), Dec = (fut$Dec) - (base$Dec)
#         ) # eo ddf
#         # summary(delta)
#
#         # Map and save
#         save(x = base, file = paste("clims_mon_logpo4_rcp85_",m,"_2031-2012.Rdata", sep = "") )
#         save(x = fut, file = paste("clims_mon_logpo4_rcp85_",m,"_2100-2081.Rdata", sep = "") )
#         save(x = delta, file = paste("clims_mon_diff_logpo4_rcp85_",m,"_2031-2100.Rdata", sep = "") )
#
#         ### Print maps
#         # For momo in months
#         for(momo in months) {
#
#             d1 <- delta[,c("x","y",momo)]
#             colnames(d1)[c(3)] <- "delta"
#
#             map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d1) ) +
#                 scale_fill_gradient2(name = "Delta logPO4 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
#                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                     labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                     labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#             ggsave(plot = map, filename = paste("map_delta_logpo4_rcp85_",m,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
#
#         } # eo for loop
 
 
        ### C) logSiO2 ------------------------------------------------------------------------------------------------------------------
        message(paste("log transforming SiO2 for ", m, sep = "")) 
        base <- get(load(paste("clims_mon_sio2_rcp85_",m,"_2031-2012.Rdata", sep = "")))
        fut <- get(load(paste("clims_mon_sio2_rcp85_",m,"_2100-2081.Rdata", sep = "")))
        # log transform
        # summary(log1p(base[,c(4:15)])); summary(log1p(fut[,c(4:15)]))
        base[,c(4:15)] <- log1p(base[,c(4:15)])
        fut[,c(4:15)] <- log1p(fut[,c(4:15)])
        # summary(base); summary(fut)
        
        # Compute delta between future and baseline
		delta <- data.frame(id = base$id, x = base$x, y = base$y, 
			    Jan = (fut$Jan) - (base$Jan), Feb = (fut$Feb) - (base$Feb), 
                Mar = (fut$Mar) - (base$Mar), Apr = (fut$Apr) - (base$Apr), 
                May = (fut$May) - (base$May), Jun = (fut$Jun) - (base$Jun),
			    Jul = (fut$Jul) - (base$Jul), Aug = (fut$Aug) - (base$Aug), 
                Sep = (fut$Sep) - (base$Sep), Oct = (fut$Oct) - (base$Oct), 
                Nov = (fut$Nov) - (base$Nov), Dec = (fut$Dec) - (base$Dec)
		) # eo ddf
        # summary(delta)
        
        # Map and save
        save(x = base, file = paste("clims_mon_logsio2_rcp85_",m,"_2031-2012.Rdata", sep = "") )
        save(x = fut, file = paste("clims_mon_logsio2_rcp85_",m,"_2100-2081.Rdata", sep = "") )
        save(x = delta, file = paste("clims_mon_diff_logsio2_rcp85_",m,"_2031-2100.Rdata", sep = "") )
        
        ### Print maps
        # For momo in months
        for(momo in months) {
            
            d1 <- delta[,c("x","y",momo)]
            colnames(d1)[c(3)] <- "delta"
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d1) ) +
                scale_fill_gradient2(name = "Delta logSiO2 (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
            ggsave(plot = map, filename = paste("map_delta_logsio2_rcp85_",m,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
                 
        } # eo for loop

} # eo for m in ESM

### Check distrib of logPO4, logSiO2
d <- get(load("clims_mon_logsio2_rcp85_MRI-NEMURO_2100-2081.Rdata"))
summary(d) # OK



### 1.E) Extract and format CHL (nano+diat) monthly outputs  ----------------------------------------------------------------------
# diatomsCHL1.nc + nanoCHL1.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Monthly forecast files by Lotte (ALL MODELS)")
nc <- nc_open("diatomsCHL1.nc")
nc  # units =  kg.m-3 --> will need to convert to  mg.m-3 --> x1000000
names(nc$var) # DIAT_pisces + DIAT_CNRM # ok, use both
### Choose ESM
rasD <- raster::stack("diatomsCHL1.nc", varname = "DIAT_CNRM")
### Nano now
nc <- nc_open("nanoCHL1.nc")
nc  # units =  kg.m-3 --> will need to convert to  mg.m-3 --> x1000000
names(nc$var) # DIAT_pisces + DIAT_CNRM # ok, use both
rasN <- raster::stack("nanoCHL1.nc", varname = "DIAT_CNRM")

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
# Ok, seems like a good start
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
# paste(m,l, sep = "_") # yeah
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
rD <- as.data.frame(rasD, xy = T)
dim(rD)
colnames(rD)[c(3:length(rD))] <- cols
head(rD[,c(1:25)])
summary(rD[,c(1:10)])
# Same for Nanophyto
rN <- as.data.frame(rasN, xy = T)
dim(rN)
colnames(rN)[c(3:length(rN))] <- cols
head(rN[,c(1:25)])
summary(rN[,c(1:10)])

# Ok, select the columns you are not interested in
baselineD <- rD[,c(1,2,3:242)]
futD <- rD[,c(1,2,831:length(rD))]
# Add an id 
baselineD$id <- factor(paste(baselineD$x, baselineD$y, sep = "_"))
futD$id <- factor(paste(futD$x, futD$y, sep = "_"))
# Do the same for Nano
baselineN <- rN[,c(1,2,3:242)]
futN <- rN[,c(1,2,831:length(rN))]
baselineN$id <- factor(paste(baselineN$x, baselineN$y, sep = "_"))
futN$id <- factor(paste(futN$x, futN$y, sep = "_"))

# Melt, and no need to convert to SST !
m_baseD <- melt(baselineD, id.vars = c("id","x","y")) ; m_futD <- melt(futD, id.vars = c("id","x","y"))
colnames(m_baseD)[c(4,5)] <- c("month","Diat") ; colnames(m_futD)[c(4,5)] <- c("month","Diat")
m_baseN <- melt(baselineN, id.vars = c("id","x","y")) ; m_futN <- melt(futN, id.vars = c("id","x","y"))
colnames(m_baseN)[c(4,5)] <- c("month","Nano") ; colnames(m_futN)[c(4,5)] <- c("month","Nano")

# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_baseD$m <- data.frame(do.call(rbind, strsplit(as.character(m_baseD$month), "_") ))[,1]
m_futD$m <- data.frame(do.call(rbind, strsplit(as.character(m_futD$month), "_") ))[,1]
m_baseN$m <- data.frame(do.call(rbind, strsplit(as.character(m_baseN$month), "_") ))[,1]
m_futN$m <- data.frame(do.call(rbind, strsplit(as.character(m_futN$month), "_") ))[,1]

### Ok, add them to get total chl
m_base <- data.frame(id = m_baseD$id, x = m_baseD$x, y = m_baseD$y, CHL = (m_baseD$Diat)+(m_baseN$Nano), m = m_baseD$m)
m_fut <- data.frame(id = m_futD$id, x = m_futD$x, y = m_futD$y, CHL = (m_futD$Diat)+(m_futN$Nano), m = m_futD$m)

# Clean some stuff
rm(futD, futN, baselineN, baselineD, rD, rN, rasD, rasN, m_futD, m_futN, m_baseD, m_baseN); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), CHL = mean(CHL) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), CHL = mean(CHL) ) )
head(clims_base); head(clims_fut)
# Convert from kg to to mg --> *1000000
# summary(log10((clims_base$CHL)*1000000)); summary(log10((clims_fut$CHL)*1000000))
clims_base$CHL <- (clims_base$CHL)*1000000
clims_fut$CHL <- (clims_fut$CHL)*1000000

### And log10 transform
clims_base$logChl <- log10(clims_base$CHL)
clims_fut$logChl <- log10(clims_fut$CHL)

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "logChl")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "logChl")
dim(clims_base); dim(clims_fut)
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)
# nrow(na.omit(clims_base[clims_base$Nov < -3,])) # only 286 points though

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(deltas, file = "clims_mon_diff_logchl_rcp85_CNRM-PISCES_2031-2100.Rdata")
save(clims_base, file = "clims_mon_logchl_rcp85_CNRM-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_logchl_rcp85_CNRM-PISCES_2100-2081.Rdata")
# And print maps
for(momo in months) {

    d <- deltas[,c("x","y",momo)]
    colnames(d)[c(3)] <- "delta"

    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
        scale_fill_gradient2(name = "Delta CHL log(mg/m3)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    ggsave(plot = map, filename = paste("map_delta_logchl_rcp85_","CNRM-PISCES","_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
    rm(d);gc()
} # eo for loop



# diatomsCHL2.nc + nanoCHL2.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Monthly forecast files by Lotte (ALL MODELS)")
nc <- nc_open("diatomsCHL2.nc")
nc  # units =  kg.m-3 --> will need to convert to  mg.m-3 --> x1000000
names(nc$var) # DIAT_GFDL + DIAT_nemuro
### Choose ESM
rasD <- raster::stack("diatomsCHL2.nc", varname = "DIAT_GFDL")
### Nano now
nc <- nc_open("nanoCHL2.nc")
nc  # units =  kg.m-3 --> will need to convert to  mg.m-3 --> x1000000
names(nc$var) # DIAT_GFDL + DIAT_nemuro
rasN <- raster::stack("nanoCHL2.nc", varname = "DIAT_GFDL")

length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
rD <- as.data.frame(rasD, xy = T)
colnames(rD)[c(3:length(rD))] <- cols
rN <- as.data.frame(rasN, xy = T)
colnames(rN)[c(3:length(rN))] <- cols

# Ok, select the columns you are not interested in
baselineD <- rD[,c(1,2,3:242)]
futD <- rD[,c(1,2,831:length(rD))]
# Add an id 
baselineD$id <- factor(paste(baselineD$x, baselineD$y, sep = "_"))
futD$id <- factor(paste(futD$x, futD$y, sep = "_"))
# Do the same for Nano
baselineN <- rN[,c(1,2,3:242)]
futN <- rN[,c(1,2,831:length(rN))]
baselineN$id <- factor(paste(baselineN$x, baselineN$y, sep = "_"))
futN$id <- factor(paste(futN$x, futN$y, sep = "_"))

# Melt, and no need to convert to SST !
m_baseD <- melt(baselineD, id.vars = c("id","x","y")) ; m_futD <- melt(futD, id.vars = c("id","x","y"))
colnames(m_baseD)[c(4,5)] <- c("month","Diat") ; colnames(m_futD)[c(4,5)] <- c("month","Diat")
m_baseN <- melt(baselineN, id.vars = c("id","x","y")) ; m_futN <- melt(futN, id.vars = c("id","x","y"))
colnames(m_baseN)[c(4,5)] <- c("month","Nano") ; colnames(m_futN)[c(4,5)] <- c("month","Nano")

# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_baseD$m <- data.frame(do.call(rbind, strsplit(as.character(m_baseD$month), "_") ))[,1]
m_futD$m <- data.frame(do.call(rbind, strsplit(as.character(m_futD$month), "_") ))[,1]
m_baseN$m <- data.frame(do.call(rbind, strsplit(as.character(m_baseN$month), "_") ))[,1]
m_futN$m <- data.frame(do.call(rbind, strsplit(as.character(m_futN$month), "_") ))[,1]

### Ok, add them to get total chl
m_base <- data.frame(id = m_baseD$id, x = m_baseD$x, y = m_baseD$y, CHL = (m_baseD$Diat)+(m_baseN$Nano), m = m_baseD$m)
m_fut <- data.frame(id = m_futD$id, x = m_futD$x, y = m_futD$y, CHL = (m_futD$Diat)+(m_futN$Nano), m = m_futD$m)

# Clean some stuff
rm(futD, futN, baselineN, baselineD, rD, rN, rasD, rasN, m_futD, m_futN, m_baseD, m_baseN); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), CHL = mean(CHL) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), CHL = mean(CHL) ) )
head(clims_base); head(clims_fut)
# Convert from kg to to mg --> *1000000
# summary(log10((clims_base$CHL)*1000000)); summary(log10((clims_fut$CHL)*1000000))
### WATCHOUT: FOR GFDL-TOPAZ UNITS ARE ALREADY mg of Chl-a
# clims_base$CHL <- (clims_base$CHL)*1000000
# clims_fut$CHL <- (clims_fut$CHL)*1000000

### And log10 transform
clims_base$logChl <- log10(clims_base$CHL)
clims_fut$logChl <- log10(clims_fut$CHL)

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "logChl")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "logChl")
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)
# nrow(na.omit(clims_base[clims_base$Nov < -3,])) # only 286 points though

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(deltas, file = "clims_mon_diff_logchl_rcp85_GFDL-TOPAZ_2031-2100.Rdata")
save(clims_base, file = "clims_mon_logchl_rcp85_GFDL-TOPAZ_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_logchl_rcp85_GFDL-TOPAZ_2100-2081.Rdata")
# And print maps
for(momo in months) {

    d <- deltas[,c("x","y",momo)]
    colnames(d)[c(3)] <- "delta"

    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
        scale_fill_gradient2(name = "Delta CHL log(mg/m3)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    ggsave(plot = map, filename = paste("map_delta_logchl_rcp85_","GFDL-TOPAZ","_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
    rm(d);gc()
} # eo for loop



# diatomsCHL3.nc + nanoCHL3.nc
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Monthly forecast files by Lotte (ALL MODELS)")
nc <- nc_open("diatomsCHL3.nc")
nc  # units =  kg.m-3 --> will need to convert to  mg.m-3 --> x1000000
names(nc$var) # "DIAT_Hadgem" "DIAT_BEC"   
rasD <- raster::stack("diatomsCHL3.nc", varname = "DIAT_BEC")

### Nano now
nc <- nc_open("nanoCHL3.nc")
nc  # units =  kg.m-3 --> will need to convert to  mg.m-3 --> x1000000
names(nc$var) # "DIAT_Hadgem" "DIAT_BEC"   
rasN <- raster::stack("nanoCHL3.nc", varname = "DIAT_BEC")

length(rep(months, times = 89)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(2012:2100, each=12))
m <- rep(months, times = 89)
l <- rep(2012:2100, each= 12)
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
rD <- as.data.frame(rasD, xy = T)
colnames(rD)[c(3:length(rD))] <- cols
rN <- as.data.frame(rasN, xy = T)
colnames(rN)[c(3:length(rN))] <- cols

# Ok, select the columns you are not interested in
baselineD <- rD[,c(1,2,3:242)]
futD <- rD[,c(1,2,831:length(rD))]
# Add an id 
baselineD$id <- factor(paste(baselineD$x, baselineD$y, sep = "_"))
futD$id <- factor(paste(futD$x, futD$y, sep = "_"))
# Do the same for Nano
baselineN <- rN[,c(1,2,3:242)]
futN <- rN[,c(1,2,831:length(rN))]
baselineN$id <- factor(paste(baselineN$x, baselineN$y, sep = "_"))
futN$id <- factor(paste(futN$x, futN$y, sep = "_"))

# Melt, and no need to convert to SST !
m_baseD <- melt(baselineD, id.vars = c("id","x","y")) ; m_futD <- melt(futD, id.vars = c("id","x","y"))
colnames(m_baseD)[c(4,5)] <- c("month","Diat") ; colnames(m_futD)[c(4,5)] <- c("month","Diat")
m_baseN <- melt(baselineN, id.vars = c("id","x","y")) ; m_futN <- melt(futN, id.vars = c("id","x","y"))
colnames(m_baseN)[c(4,5)] <- c("month","Nano") ; colnames(m_futN)[c(4,5)] <- c("month","Nano")

# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_baseD$m <- data.frame(do.call(rbind, strsplit(as.character(m_baseD$month), "_") ))[,1]
m_futD$m <- data.frame(do.call(rbind, strsplit(as.character(m_futD$month), "_") ))[,1]
m_baseN$m <- data.frame(do.call(rbind, strsplit(as.character(m_baseN$month), "_") ))[,1]
m_futN$m <- data.frame(do.call(rbind, strsplit(as.character(m_futN$month), "_") ))[,1]

### Ok, add them to get total chl
m_base <- data.frame(id = m_baseD$id, x = m_baseD$x, y = m_baseD$y, CHL = (m_baseD$Diat)+(m_baseN$Nano), m = m_baseD$m)
m_fut <- data.frame(id = m_futD$id, x = m_futD$x, y = m_futD$y, CHL = (m_futD$Diat)+(m_futN$Nano), m = m_futD$m)

# Clean some stuff
rm(futD, futN, baselineN, baselineD, rD, rN, rasD, rasN, m_futD, m_futN, m_baseD, m_baseN); gc()

# Compute the 12 monthly clims
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), CHL = mean(CHL) ) )
clims_fut <- data.frame( m_fut %>% group_by(id,m) %>% summarise(x = unique(x), y = unique(y), CHL = mean(CHL) ) )
head(clims_base); head(clims_fut)
# Convert from kg to to mg --> *1000000
# summary(log10((clims_base$CHL))); summary(log10((clims_fut$CHL)))
#clims_base$CHL <- (clims_base$CHL)*1000000
#clims_fut$CHL <- (clims_fut$CHL)*1000000
### CESM-BEC already in mg/m3 !!!!

### And log10 transform
clims_base$logChl <- log10(clims_base$CHL)
clims_fut$logChl <- log10(clims_fut$CHL)

### dcast to put the 12 months and coordinates as columns
clims_base <- dcast(data = clims_base, id + x + y ~ m, value.var = "logChl")
clims_fut <- dcast(data = clims_fut, id + x + y ~ m, value.var = "logChl")
head(clims_base); head(clims_fut)
summary(clims_base); summary(clims_fut)
# nrow(na.omit(clims_base[clims_base$Nov < -3,])) # only 286 points though

# Compute the 12 monthly deltas SST and map
deltas <- data.frame(id = clims_base$id, x = clims_base$x, y = clims_base$y, 
	Jan = (clims_fut$Jan) - (clims_base$Jan), Feb = (clims_fut$Feb) - (clims_base$Feb), Mar = (clims_fut$Mar) - (clims_base$Mar), 
	Apr = (clims_fut$Apr) - (clims_base$Apr), May = (clims_fut$May) - (clims_base$May), Jun = (clims_fut$Jun) - (clims_base$Jun),
	Jul = (clims_fut$Jul) - (clims_base$Jul), Aug = (clims_fut$Aug) - (clims_base$Aug), Sep = (clims_fut$Sep) - (clims_base$Sep), 
	Oct = (clims_fut$Oct) - (clims_base$Oct), Nov = (clims_fut$Nov) - (clims_base$Nov), Dec = (clims_fut$Dec) - (clims_base$Dec)
) # eo ddf
summary(deltas)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(deltas, file = "clims_mon_diff_logchl_rcp85_CESM-BEC_2031-2100.Rdata")
save(clims_base, file = "clims_mon_logchl_rcp85_CESM-BEC_2031-2012.Rdata")
save(clims_fut, file = "clims_mon_logchl_rcp85_CESM-BEC_2100-2081.Rdata")

# And print maps
for(momo in months) {

    d <- deltas[,c("x","y",momo)]
    colnames(d)[c(3)] <- "delta"

    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
        scale_fill_gradient2(name = "Delta CHL log(mg/m3)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    ggsave(plot = map, filename = paste("map_delta_logchl_rcp85_","CESM-BEC","_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
    rm(d);gc()
} # eo for loop


# -----------------------------------------------------------------

##### 2°) Derive the 3 secondary variables: dSST from SST climatologies, N* and Si*

### 2.A) Compute annual SST range (dSST = max - min) in a for loop
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
ESM <- c("CNRM-PISCES","IPSL-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")
# m <- "CNRM-PISCES" # for testing
for(m in ESM) {
    
        message(paste("Computing dSST for ", m, sep = "")) 
        base <- get(load(paste("clims_mon_sst_rcp85_",m,"_2031-2012.Rdata", sep = "")))
        fut <- get(load(paste("clims_mon_sst_rcp85_",m,"_2100-2081.Rdata", sep = "")))
		require("matrixStats")
		base$dSST <- (rowMaxs(as.matrix(base[,c(4:15)]), na.rm = F)) - (rowMins(as.matrix(base[,c(4:15)]), na.rm = F))
        fut$dSST <- (rowMaxs(as.matrix(fut[,c(4:15)]), na.rm = F)) - (rowMins(as.matrix(fut[,c(4:15)]), na.rm = F))
        delta <- data.frame(id = fut$id, x = fut$x, y = fut$y, dSST = (fut$dSST) - (base$dSST) )
        # summary(delta$dSST)
        base2 <- base[,c(1:3,16)]
        fut2 <- fut[,c(1:3,16)]
        # dim(base2); dim(fut2)
        
        save(x = base2, file = paste("clims_ann_dsst_rcp85_",m,"_2031-2012.Rdata", sep = "") )
        save(x = fut2, file = paste("clims_ann_dsst_rcp85_",m,"_2100-2081.Rdata", sep = "") )
        save(x = delta, file = paste("clims_ann_diff_dsst_rcp85_",m,"_2031-2100.Rdata", sep = "") )
        
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dSST), data = na.omit(delta) ) +
            scale_fill_gradient2(name = "Delta dSST (°C)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
        ggsave(plot = map, filename = paste("map_delta_dSST_rcp85_",m,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
        
} # eo for loop


### 2.B) Compute monthly N* (no3 - 16*po4) and Si* (sio2 - no3)
m <- "CNRM-PISCES" # for testing
for(m in ESM) {
    
        message(paste("Computing N* and Si* for ", m, sep = "")) 
        base_no3 <- get(load(paste("clims_mon_no3_rcp85_",m,"_2031-2012.Rdata", sep = "")))
        fut_no3 <- get(load(paste("clims_mon_no3_rcp85_",m,"_2100-2081.Rdata", sep = "")))
        base_po4 <- get(load(paste("clims_mon_po4_rcp85_",m,"_2031-2012.Rdata", sep = "")))
        fut_po4 <- get(load(paste("clims_mon_po4_rcp85_",m,"_2100-2081.Rdata", sep = "")))
        base_sio2 <- get(load(paste("clims_mon_sio2_rcp85_",m,"_2031-2012.Rdata", sep = "")))
        fut_sio2 <- get(load(paste("clims_mon_sio2_rcp85_",m,"_2100-2081.Rdata", sep = "")))
		
		### Compute N* : no3 - (16*(po4))
		base_nstar <- data.frame(id = base_no3$id, x = base_no3$x, y = base_no3$y,
			Jan = base_no3$Jan - (16*(base_po4$Jan)), Feb = base_no3$Feb - (16*(base_po4$Feb)), Mar = base_no3$Mar - (16*(base_po4$Mar)), 
			Apr = base_no3$Apr - (16*(base_po4$Apr)), May = base_no3$May - (16*(base_po4$May)), Jun = base_no3$Jun - (16*(base_po4$Jun)),
			Jul = base_no3$Jul - (16*(base_po4$Jul)), Aug = base_no3$Aug - (16*(base_po4$Aug)), Sep = base_no3$Sep - (16*(base_po4$Sep)), 
			Oct = base_no3$Oct - (16*(base_po4$Oct)), Nov = base_no3$Nov - (16*(base_po4$Nov)), Dec = base_no3$Jan - (16*(base_po4$Jan))
		) # eo ddf
	
		fut_nstar <- data.frame(id = fut_no3$id, x = fut_no3$x, y = fut_no3$y,
			Jan = fut_no3$Jan - (16*(fut_po4$Jan)), Feb = fut_no3$Feb - (16*(fut_po4$Feb)), Mar = fut_no3$Mar - (16*(fut_po4$Mar)), 
			Apr = fut_no3$Apr - (16*(fut_po4$Apr)), May = fut_no3$May - (16*(fut_po4$May)), Jun = fut_no3$Jun - (16*(fut_po4$Jun)),
			Jul = fut_no3$Jul - (16*(fut_po4$Jul)), Aug = fut_no3$Aug - (16*(fut_po4$Aug)), Sep = fut_no3$Sep - (16*(fut_po4$Sep)), 
			Oct = fut_no3$Oct - (16*(fut_po4$Oct)), Nov = fut_no3$Nov - (16*(fut_po4$Nov)), Dec = fut_no3$Jan - (16*(fut_po4$Jan))
		) # eo ddf
	
		# Compute delta N*
		deltas_nstar <- data.frame(id = fut$id, x = fut$x, y = fut$y, 
			Jan = (fut_nstar$Jan) - (base_nstar$Jan), Feb = (fut_nstar$Feb) - (base_nstar$Feb), Mar = (fut_nstar$Mar) - (base_nstar$Mar), 
			Apr = (fut_nstar$Apr) - (base_nstar$Apr), May = (fut_nstar$May) - (base_nstar$May), Jun = (fut_nstar$Jun) - (base_nstar$Jun),
			Jul = (fut_nstar$Jul) - (base_nstar$Jul), Aug = (fut_nstar$Aug) - (base_nstar$Aug), Sep = (fut_nstar$Sep) - (base_nstar$Sep), 
			Oct = (fut_nstar$Oct) - (base_nstar$Oct), Nov = (fut_nstar$Nov) - (base_nstar$Nov), Dec = (fut_nstar$Dec) - (base_nstar$Dec)
		) # eo ddf
        
        ### Compute Si* : no3 - (16*(po4))
		base_sistar <- data.frame(id = base_no3$id, x = base_no3$x, y = base_no3$y,
			Jan = (base_sio2$Jan)-(base_no3$Jan), Feb = (base_sio2$Feb)-(base_no3$Feb), Mar = (base_sio2$Mar)-(base_no3$Mar), 
			Apr = (base_sio2$Apr)-(base_no3$Apr), May = (base_sio2$May)-(base_no3$May), Jun = (base_sio2$Jun)-(base_no3$Jun),
			Jul = (base_sio2$Jul)-(base_no3$Jul), Aug = (base_sio2$Aug)-(base_no3$Aug), Sep = (base_sio2$Sep)-(base_no3$Sep), 
			Oct = (base_sio2$Oct)-(base_no3$Oct), Nov = (base_sio2$Nov)-(base_no3$Nov), Dec = (base_sio2$Dec)-(base_no3$Dec)
		) # eo ddf
	
		fut_sistar <- data.frame(id = fut_no3$id, x = fut_no3$x, y = fut_no3$y,
			Jan = (fut_sio2$Jan)-(fut_no3$Jan), Feb = (fut_sio2$Feb)-(fut_no3$Feb), Mar = (fut_sio2$Mar)-(fut_no3$Mar), 
			Apr = (fut_sio2$Apr)-(fut_no3$Apr), May = (fut_sio2$May)-(fut_no3$May), Jun = (fut_sio2$Jun)-(fut_no3$Jun),
			Jul = (fut_sio2$Jul)-(fut_no3$Jul), Aug = (fut_sio2$Aug)-(fut_no3$Aug), Sep = (fut_sio2$Sep)-(fut_no3$Sep), 
			Oct = (fut_sio2$Oct)-(fut_no3$Oct), Nov = (fut_sio2$Nov)-(fut_no3$Nov), Dec = (fut_sio2$Dec)-(fut_no3$Dec)
		) # eo ddf
        
		# Compute delta Si*
		deltas_sistar <- data.frame(id = fut$id, x = fut$x, y = fut$y, 
			Jan = (fut_sistar$Jan) - (base_sistar$Jan), Feb = (fut_sistar$Feb) - (base_sistar$Feb), Mar = (fut_sistar$Mar) - (base_sistar$Mar), 
			Apr = (fut_sistar$Apr) - (base_sistar$Apr), May = (fut_sistar$May) - (base_sistar$May), Jun = (fut_sistar$Jun) - (base_sistar$Jun),
			Jul = (fut_sistar$Jul) - (base_sistar$Jul), Aug = (fut_sistar$Aug) - (base_sistar$Aug), Sep = (fut_sistar$Sep) - (base_sistar$Sep), 
			Oct = (fut_sistar$Oct) - (base_sistar$Oct), Nov = (fut_sistar$Nov) - (base_sistar$Nov), Dec = (fut_sistar$Dec) - (base_sistar$Dec)
		) # eo ddf
        
        # Save N* files
        save(x = base_nstar, file = paste("clims_mon_nstar_rcp85_",m,"_2031-2012.Rdata", sep = "") )
        save(x = fut_nstar, file = paste("clims_mon_nstar_rcp85_",m,"_2100-2081.Rdata", sep = "") )
        save(x = deltas_nstar, file = paste("clims_mon_diff_nstar_rcp85_",m,"_2031-2100.Rdata", sep = "") )
        # Save Si* files
        save(x = base_sistar, file = paste("clims_mon_sistar_rcp85_",m,"_2031-2012.Rdata", sep = "") )
        save(x = fut_sistar, file = paste("clims_mon_sistar_rcp85_",m,"_2100-2081.Rdata", sep = "") )
        save(x = deltas_sistar, file = paste("clims_mon_diff_sistar_rcp85_",m,"_2031-2100.Rdata", sep = "") )
        
        ### Print maps
        months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
        for(momo in months) {
            
            d1 <- deltas_nstar[,c("x","y",momo)]
            colnames(d1)[c(3)] <- "delta"
            d2 <- deltas_sistar[,c("x","y",momo)]
            colnames(d2)[c(3)] <- "delta"
            
            map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d1) ) +
                scale_fill_gradient2(name = "Delta N* (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
            ggsave(plot = map1, filename = paste("map_delta_nstar_rcp85_",m,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
            
            map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d2) ) +
                scale_fill_gradient2(name = "Delta Si* (µM)\n(2100-2031)", low = "blue", high = "red", mid = "white") +
                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
            ggsave(plot = map2, filename = paste("map_delta_sistar_rcp85_",m,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
                 
        } # eo for loop
        
} # eo for loop

### Ok, check the distribution of Si* and N*
dir()
d <- get(load("clims_mon_sistar_rcp85_MRI-NEMURO_2031-2012.Rdata"))
summary(d)


# -----------------------------------------------------------------

##### 3°) Extract yearly O2 outputs from light.nc for each ESM, get baseline and future monthly means, and derive seasonality based on in situ monthly climatologies

### 3.A) Compute mean annual climatologies of O2 for each of the 5 ESMs
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Irradiance files (light.nc by Lotte)")
nc <- nc_open("light.nc")
nc 
# O2 unit = W/m^2
# time --> Size:89 89 years

names(nc$var) # First 5 levels to be used ! You may extract this in a for loop
ESM <- names(nc$var)[1:5]; ESM

### In a for loop, extract the yearly outputs, compute annual means over the two periods of interest and save annual mean
# esm <- ESM[1] # for testing
for(esm in ESM) {
    
        ### Choose ESM
        setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/irradiance_data_light.nc_byLotte")
        message(paste("Preparing annual climatologies for ", esm, sep = ""))
        # ras <- raster::stack("light.nc")
        ras <- raster::stack("light.nc", varname = esm)
        # Extract data for RCP8.5 and provide colnames
        r <- as.data.frame(ras, xy = T)
        # dim(r); summary(r) --> good unit, good range of values 
        rm(ras); gc()

        # Ok, now add the years as colnames: vector of c(2012:2100); length(c(2012:2100))
        colnames(r)[c(3:length(r))] <- c(2012:2100)

        # Now, compute annual mean over same period of interest as the other variables (2012-2031 & 2081-2100)
        baseline <- r[,c(1,2,3:22)]
        fut <- r[,c(1,2,72:length(r))]

        # Add an id 
        baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
        fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
        # head(baseline$id); head(fut$id)
        ### IDs don't seem to match previous IDS...you can still compute annual averages and then re-grid though
        # Compute the 12 monthly clims after melting
        require("dplyr","reshape2")
        m_base <- melt(baseline, id.vars = c("id","x","y"))
        m_fut <- melt(fut, id.vars = c("id","x","y"))
        colnames(m_base)[c(4,5)] <- c("year","O2") # it's called O2 because I was lazy, but it is PAR no worries
        colnames(m_fut)[c(4,5)] <- c("year","O2")

        # Clean some stuff
        rm(fut, baseline, r); gc()

        # Compute the annual clim
        clims_base <- data.frame( m_base %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), O2 = mean(O2) ) )
        clims_fut <- data.frame( m_fut %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), O2 = mean(O2) ) )
        # summary(clims_base); summary(clims_fut)
    
        # OK, before saving, change x and y so they match the others 
        ### !!! ONLY FOR DATA_pisces ¡¡¡
        if( esm == "DATA_pisces") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
                clims_base$x <- (clims_base$x) - 0.5
                clims_base$y <- (clims_base$y) - 90.5
                clims_base$id <- paste(clims_base$x, clims_base$y, sep = "_")
                clims_base <- clims_base[order(clims_base$id),]
                # And same for clim_gut
                clims_fut$x <- (clims_fut$x) - 0.5
                clims_fut$y <- (clims_fut$y) - 90.5
                clims_fut$id <- paste(clims_fut$x, clims_fut$y, sep = "_")
                clims_fut <- clims_fut[order(clims_fut$id),]
                # head(clims_base); head(clims_fut); summary(clims_base); summary(clims_fut) 
        } # eo if loop

        ### Save and print annual maps
        # Choose ESM name based on 'esm' 
        if( esm == "DATA_nemuro" ) {
                name <- "MRI-NEMURO"
        } else if( esm == "DATA_BEC" ) {
                name <- "CESM-BEC"
        } else if( esm == "DATA_pisces" ) {
                name <- "IPSL-PISCES"
        } else if( esm == "DATA_GFDL" ) {
                name <- "GFDL-TOPAZ"
        } else if( esm == "DATA_CNRM" ) {
                name <- "CNRM-PISCES"
        } # eo else if loop
        
        message(paste("Saving annual climatologies for ", name, sep = ""))
        setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
        
        # Save climatologies
        save(x = clims_base, file = paste("clims_ann_par_rcp85_",name,"_2031-2012.Rdata", sep = "") )
        save(x = clims_fut, file = paste("clims_ann_par_rcp85_",name,"_2100-2081.Rdata", sep = "") )
        
        # Maps
        world2 <- map_data("world2")
        map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = O2), data = na.omit(clims_base) ) +
            scale_fill_viridis(name = "O2 (W/m2)\n(2031-2012)" ) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
        map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = O2), data = na.omit(clims_fut) ) +
            scale_fill_viridis(name = "O2 (W/m2)\n(2100-2081)" ) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                
        ggsave(plot = map1, filename = paste("map_ann_par_rcp85_",name,"_2031-2012.pdf", sep = ""), height = 3, width = 6, dpi = 300)
        ggsave(plot = map2, filename = paste("map_ann_par_rcp85_",name,"_2100-2081.pdf", sep = ""), height = 3, width = 6, dpi = 300)
    
} # eo for loop


### 3.B) Extract anomalies of monthly O2 climatologies to annual mean and derive baseline and future modelled monthly O2 for each ESM
# 3.B.1) Get in situ monthly means, compute annual O2 average and derive anomalies
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
# m <- "apr"
res <- lapply(months, function(m) {
            d <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")[,c("x","y","O2")]
            d$id <- paste(d$x, d$y, sep = "_")
            d <- d[order(d$id),]
            d$month <- m
            return(d)
    } # eo FUN
) # eo lapply
# Cbind
ddf <- do.call(rbind, res)
head(ddf)
# Dcast to put monthly means as columns
clims <- dcast(data = ddf, id + x + y ~ month, value.var = "O2")
dim(clims); summary(clims)

# Compute annual average 
require("matrixStats")
clims$ann <- rowMeans(as.matrix(clims[,c(4:15)]), na.rm = T)
summary(clims$ann) # summary(clims$Ann)
### Change colnames and compute difference to annual mean
colnames(clims)[c(4:16)] <- c("Apr","Aug","Dec","Feb","Jan","Jul","Jun","Mar","May","Nov","Oct","Sep","Ann")
anoms <- data.frame(id = clims$id, x = clims$x, y = clims$y, 
            Jan = (clims$Jan)-(clims$Ann), Feb = (clims$Feb)-(clims$Ann), Mar = (clims$Mar)-(clims$Ann),
            Apr = (clims$Apr)-(clims$Ann), May = (clims$May)-(clims$Ann), Jun = (clims$Jun)-(clims$Ann),
            Jul = (clims$Jul)-(clims$Ann), Aug = (clims$Aug)-(clims$Ann), Sep = (clims$Sep)-(clims$Ann),
            Oct = (clims$Oct)-(clims$Ann), Nov = (clims$Nov)-(clims$Ann), Dec = (clims$Dec)-(clims$Ann)
) # eo ddf
summary(anoms)
# Rotate x coordinates
anoms$x2 <- anoms$x 
anoms[anoms$x < 0 ,"x2"] <- (anoms[anoms$x < 0 ,"x"]) + 360
# And modify cell id accordingly
anoms$id <- paste(anoms$x2, anoms$y, sep = "_")

### Print maps of monthly anomalies to the annual mean 
months <- c("Apr","Aug","Dec","Feb","Jan","Jul","Jun","Mar","May","Nov","Oct","Sep")
for(momo in months) {

    d <- anoms[,c("x2","y",momo)]
    colnames(d)[c(3)] <- "delta"

    map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = delta), data = d) +
        scale_fill_gradient2(name = "O2 anomaly\n(W/m2)", low = "#3288bd", high = "#d53e4f", mid = "white") +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    ggsave(plot = map, filename = paste("map_mon_anom_O2_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
    rm(d);gc()
} # eo for loop

### Ok nice, save anoms before deriving baseline and future seasonal cycle
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Irradiance files (light.nc by Lotte)")
anoms <- anoms[order(anoms$id),]
save(x = anoms, file = "clims_mon+ann_anoms_par_insitu.Rdata")


# 3.B.2) Using monthly anoms, derive monthly climatologies from baseline and future annual clims (looping over ESM)
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
ESM <- c("IPSL-PISCES","CNRM-PISCES","GFDL-TOPAZ","MRI-NEMURO","CESM-BEC")
# esm <- "CNRM-PISCES" # for testing
# For loop
for(esm in ESM) {
    
        # Get annual data
        message(paste("Preparing baseline and future monthly O2 for ", esm, sep = ""))
        message(paste("", sep = ""))
        clim_base <- get(load(paste("clims_ann_par_rcp85_",esm,"_2031-2012.Rdata", sep = "")))
        clim_fut <- get(load(paste("clims_ann_par_rcp85_",esm,"_2100-2081.Rdata", sep = "")))
        # head(clim_base); head(clim_fut); head(anoms)
        # tail(clim_base); tail(clim_fut); tail(anoms)
        
        # Derive monthly clims from anoms for baseline
        clim_base$Jan <- (anoms$Jan) + (clim_base$O2)
        clim_base$Feb <- (anoms$Feb) + (clim_base$O2)
        clim_base$Mar <- (anoms$Mar) + (clim_base$O2)
        clim_base$Apr <- (anoms$Apr) + (clim_base$O2)
        clim_base$May <- (anoms$May) + (clim_base$O2)
        clim_base$Jun <- (anoms$Jun) + (clim_base$O2)
        clim_base$Jul <- (anoms$Jul) + (clim_base$O2)
        clim_base$Aug <- (anoms$Aug) + (clim_base$O2)
        clim_base$Sep <- (anoms$Sep) + (clim_base$O2)
        clim_base$Oct <- (anoms$Oct) + (clim_base$O2)
        clim_base$Nov <- (anoms$Nov) + (clim_base$O2)
        clim_base$Dec <- (anoms$Dec) + (clim_base$O2)
        # and now for future
        clim_fut$Jan <- (anoms$Jan) + (clim_fut$O2)
        clim_fut$Feb <- (anoms$Feb) + (clim_fut$O2)
        clim_fut$Mar <- (anoms$Mar) + (clim_fut$O2)
        clim_fut$Apr <- (anoms$Apr) + (clim_fut$O2)
        clim_fut$May <- (anoms$May) + (clim_fut$O2)
        clim_fut$Jun <- (anoms$Jun) + (clim_fut$O2)
        clim_fut$Jul <- (anoms$Jul) + (clim_fut$O2)
        clim_fut$Aug <- (anoms$Aug) + (clim_fut$O2)
        clim_fut$Sep <- (anoms$Sep) + (clim_fut$O2)
        clim_fut$Oct <- (anoms$Oct) + (clim_fut$O2)
        clim_fut$Nov <- (anoms$Nov) + (clim_fut$O2)
        clim_fut$Dec <- (anoms$Dec) + (clim_fut$O2)
        
        # summary(clim_base); summary(clim_fut)
        # nrow(na.omit(clim_base[clim_base$Aug < 0,]))
        # Replace all negative values to zeroes
        clim_base[,c(5:16)][clim_base[,c(5:16)] < 0 & !is.na(clim_base[,c(5:16)])] <- 0
        clim_fut[,c(5:16)][clim_fut[,c(5:16)] < 0 & !is.na(clim_fut[,c(5:16)])] <- 0
        # summary(clim_base); summary(clim_fut)
        colnames(clim_base)[4] <- "Annual"
        colnames(clim_fut)[4] <- "Annual"
        # dim(clim_base); dim(clim_fut)
        # Compute delta between present and future
        # Compute the 12 monthly deltas SST and map
        deltas <- data.frame(id = clim_base$id, x = clim_base$x, y = clim_base$y, 
        	Jan = (clim_fut$Jan) - (clim_base$Jan), Feb = (clim_fut$Feb) - (clim_base$Feb), Mar = (clim_fut$Mar) - (clim_base$Mar), 
        	Apr = (clim_fut$Apr) - (clim_base$Apr), May = (clim_fut$May) - (clim_base$May), Jun = (clim_fut$Jun) - (clim_base$Jun),
        	Jul = (clim_fut$Jul) - (clim_base$Jul), Aug = (clim_fut$Aug) - (clim_base$Aug), Sep = (clim_fut$Sep) - (clim_base$Sep), 
        	Oct = (clim_fut$Oct) - (clim_base$Oct), Nov = (clim_fut$Nov) - (clim_base$Nov), Dec = (clim_fut$Dec) - (clim_base$Dec)
        ) # eo ddf
        # summary(deltas)

        ### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
        setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
        save(deltas, file = paste("clims_mon_diff_par_rcp85_",esm,"_2031-2100.Rdata", sep = "") )
        save(clim_base, file = paste("clims_mon_par_rcp85_",esm,"_2031-2012.Rdata", sep = "") )
        save(clim_fut, file = paste("clims_mon_par_rcp85_",esm,"_2100-2081.Rdata", sep = "") )

        # And print maps
        for(momo in months) {

            d <- deltas[,c("x","y",momo)]
            colnames(d)[c(3)] <- "delta"

            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
                scale_fill_gradient2(name = "Delta O2 (W/m2)\n(2100-2031)", low = "#3288bd", high = "#d53e4f", mid = "white") +
                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

            ggsave(plot = map, filename = paste("map_delta_par_rcp85_",esm,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
            rm(d);gc()
            
        } # eo for loop
    
} # eo for loop


# -----------------------------------------------------------------


##### 4°) Extract dO2 @ 175m outputs for each ESM, get baseline and future annual means like for O2 above
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Oxygen files (annual, separate models)")
Wd <- getwd()
###

### 4.A) Start with the easiest because already monthly and in one netCDF file: GFDL-TOPAZ
setwd(paste(Wd,"/","GFDL-TOPAZ (monthly)", sep = ""))
nc <- nc_open("o2_175m_1861_2100_rcp85.nc_regrid.nc")
nc 
#
ras <- raster::stack("o2_175m_1861_2100_rcp85.nc_regrid.nc", varname = "O2_REGRID")
ras # 2880 time layers

months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
length(rep(months, times = 240)) # ok, now add the years: vector of same length but with years repeating 12 times
length(rep(1861:2100, each=12))
# Ok
m <- rep(months, times = 89)
l <- rep(1861:2100, each= 12)
cols <- paste(m,l, sep = "_")

# Extract data for RCP8.5 and provide colnames
r <- as.data.frame(ras, xy = T)
dim(r)
colnames(r)[c(3:length(r))] <- cols
head(r[,c(1:25)])
summary(r[,c(1:10)])
# Ok, select the columns you are not interested in
baseline <- r[,c(1,2,1815:2054)]
fut <- r[,c(1,2,2643:length(r))]
dim(baseline); dim(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
head(baseline$id); head(fut$id)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","dO2")
colnames(m_fut)[c(4,5)] <- c("month","dO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, r, ras); gc()

# Compute the annual climatology
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
head(clims_base); head(clims_fut); dim(clims_base); dim(clims_fut)
# Convert to ml/l
# O2 (mol/kg OR mol/L): multiply by 1.029 (seawater density) and divide by 44.661e-6 (conversion from moles to L by the molar volume of O2)
clims_base$dO2 <- (clims_base$dO2)*1.029
clims_fut$dO2 <- (clims_fut$dO2)*1.029
clims_base$dO2 <- (clims_base$dO2)/44.661e-6
clims_fut$dO2 <- (clims_fut$dO2)/44.661e-6
# summary(clims_base); summary(clims_fut)
# head(clims_base); head(clims_fut)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(clims_base, file = "clims_ann_o2_rcp85_GFDL-TOPAZ_2031-2012.Rdata")
save(clims_fut, file = "clims_ann_o2_rcp85_GFDL-TOPAZ_2100-2081.Rdata")

# Maps
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = na.omit(clims_base) ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2012-2031)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_GFDL-TOPAZ_2012-2031.pdf", sep = ""), height = 3, width = 6, dpi = 300)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = na.omit(clims_fut) ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2081-2100)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_GFDL-TOPAZ_2081-2100.pdf", sep = ""), height = 3, width = 6, dpi = 300)


### -----------------------------------------------------------------

### 4.2) CESM-BEC
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Oxygen files (annual, separate models)")
Wd <- getwd()
setwd(paste(Wd,"/","CESM-BEC (monthly)", sep = ""))
dir()
# There are 89 netCDF files, one for each year ranging from 2012 to 2100, with monthly output
nc <- nc_open("CESM1_monthly_O2_2012.1x1d.nc")
nc
names(nc$var) # "O2" but not depth
# z_t  Size:15
# long_name: depth from surface to midpoint of layer
# units: centimeters !!!
# positive: down
# valid_min: 500 --> 5m
# valid_max: 537500 --> 5375m 

# Need to find the level corresponding to 175m depth
z <- ncvar_get(nc,"z_t"); z # --> need the second depth layer (150m because it's the midpoint between 100m and 200m depth)
### You mught want to try with the third level too though

### With mclapply by years, extract each yearly netCDF and combine them into one large data.frame
years <- c(2012:2100) # For the files
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") # for the colnames
# yy <- 2014 # for testing
require("parallel")
res <- mclapply(years, function(yy) {
    
            message(paste("Preparing oxygen climatologies for ",yy, sep = ""))
            ras <- raster::brick(paste("CESM1_monthly_O2_",yy,".1x1d.nc", sep = ""), varname = "O2", level = 15) #
            r <- as.data.frame(ras, xy = T)
            rm(ras); gc()
            # summary(r); head(r); dim(r)
            # Correct coordinates
            r$x <- (r$x) - 0.5
            r$y <- (r$y) - 90.5
            # Provide cell id and order
            r$id <- paste(r$x, r$y, sep = "_")
            r <- r[order(r$id),]
            # Change colnames
            colnames(r)[c(3:14)] <- paste(months, yy, sep = "_")
            # Return
            if(yy == 2012) {
                return( r[,c("id","x","y",colnames(r)[c(3:14)])] )
            } else {
                return( r[,c(colnames(r)[c(3:14)])] )
            } # eo if else loop
        
    }, mc.cores = 20
    
) # eo mclapply
# Cbind
require("dplyr")
ddf <- dplyr::bind_cols(res)
dim(ddf)
colnames(ddf)
# OK, now compute baseline and future like for GFDL-TOPAZ

# Ok, select the columns you are not interested in
baseline <- ddf[,c(1,2,3,4:243)]
fut <- ddf[,c(1,2,3,832:length(ddf))]
dim(baseline); dim(fut)
# Add an id 
head(baseline); head(fut)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","dO2")
colnames(m_fut)[c(4,5)] <- c("month","dO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean some stuff
rm(fut, baseline, ddf); gc()

# Compute the annual climatology
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
head(clims_base); head(clims_fut); dim(clims_base); dim(clims_fut)
### Convert mol/m^3 to ml/l
# O2 (mol/kg OR mol/L): multiply by 1.029 (seawater density) and divide by 44.661e-6 (conversion from moles to L by the molar volume of O2)
clims_base$dO2 <- (clims_base$dO2)*0.001
clims_fut$dO2 <- (clims_fut$dO2)*0.001
clims_base$dO2 <- (clims_base$dO2)*1.029
clims_fut$dO2 <- (clims_fut$dO2)*1.029
clims_base$dO2 <- (clims_base$dO2)/44.661e-6
clims_fut$dO2 <- (clims_fut$dO2)/44.661e-6
summary(clims_base); summary(clims_fut)
head(clims_base); head(clims_fut)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(clims_base, file = "clims_ann_o2_rcp85_CESM-BEC_2031-2012.Rdata")
save(clims_fut, file = "clims_ann_o2_rcp85_CESM-BEC_2100-2081.Rdata")

# Maps
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = na.omit(clims_base) ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2012-2031)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_CESM-BEC_2012-2031.pdf", sep = ""), height = 3, width = 6, dpi = 300)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = na.omit(clims_fut) ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2081-2100)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_CESM-BEC_2081-2100.pdf", sep = ""), height = 3, width = 6, dpi = 300)



### -----------------------------------------------------------------

### 4.3) MRI-NEMURO
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/Oxygen files (annual, separate models)")
Wd <- getwd()
setwd(paste(Wd,"/","MRI-NEMURO (monthly)", sep = ""))
# dir()
# There are 93 netCDF files, one for each year ranging from 2008 to 2100. Monthly ouput? Yes
nc <- nc_open("MRICOM_O2_2010.nc")
nc
names(nc$var) # o2

# Need to find the level corresponding to 175m depth
z <- ncvar_get(nc, "level"); z
# OK, need level 14 ! using raster::brick()
# test
yy <- 2010
ras <- raster::brick(paste("MRICOM_O2_",yy,".nc", sep = ""), varname = "O2", level = 14) #
# ras
r <- as.data.frame(ras, xy = T)
rm(ras); gc()
# summary(r); head(r); dim(r)
# Provide cell id and order
r$id <- paste(r$x, r$y, sep = "_")
r <- r[order(r$id),]
# Change colnames
colnames(r)[c(3:14)] <- paste(months, yy, sep = "_")

rm(ras); gc()

### With mclapply by years, extract each yearly netCDF and combine them into one large data.frame
years <- c(2008:2100) # For the files
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") # for the colnames
# yy <- 2014 # for testing
require("parallel")
res <- mclapply(years, function(yy) {
    
            message(paste("Preparing oxygen climatologies for ",yy, sep = ""))
            ras <- raster::brick(paste("MRICOM_O2_",yy,".nc", sep = ""), varname = "O2", level = 14) #
            r <- as.data.frame(ras, xy = T)
            rm(ras); gc()
            # summary(r); head(r); dim(r)
            # Provide cell id and order
            r$id <- paste(r$x, r$y, sep = "_")
            r <- r[order(r$id),]
            # Change colnames
            colnames(r)[c(3:14)] <- paste(months, yy, sep = "_")
            # Return
            if(yy == 2008) {
                return( r[,c("id","x","y",colnames(r)[c(3:14)])] )
            } else {
                return( r[,c(colnames(r)[c(3:14)])] )
            } # eo if else loop
        
    }, mc.cores = 23
    
) # eo mclapply
# Cbind
require("dplyr")
ddf <- dplyr::bind_cols(res)
dim(ddf)
colnames(ddf)

# Ok, select the columns you are not interested in
baseline <- ddf[,c(1,2,3,52:291)]
fut <- ddf[,c(1,2,3,880:length(ddf))]
dim(baseline); dim(fut)
# Add an id 
head(baseline); head(fut)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("month","dO2")
colnames(m_fut)[c(4,5)] <- c("month","dO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$m <- data.frame(do.call(rbind, strsplit(as.character(m_base$month), "_") ))[,1]
m_fut$m <- data.frame(do.call(rbind, strsplit(as.character(m_fut$month), "_") ))[,1]
head(m_base); head(m_fut)

# Clean stuff
rm(fut, baseline, ddf); gc()

# Compute the annual climatology
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
head(clims_base); head(clims_fut); dim(clims_base); dim(clims_fut)
### Convert mol/m^3 to ml/l
# O2 (mol/kg OR mol/L): multiply by 1.029 (seawater density) and divide by 44.661e-6 (conversion from moles to L by the molar volume of O2)
clims_base$dO2 <- (clims_base$dO2)*0.001
clims_fut$dO2 <- (clims_fut$dO2)*0.001
clims_base$dO2 <- (clims_base$dO2)*1.029
clims_fut$dO2 <- (clims_fut$dO2)*1.029
clims_base$dO2 <- (clims_base$dO2)/44.661e-6
clims_fut$dO2 <- (clims_fut$dO2)/44.661e-6
summary(clims_base); summary(clims_fut)
head(clims_base); head(clims_fut)

### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(clims_base, file = "clims_ann_o2_rcp85_MRI-NEMURO_2031-2012.Rdata")
save(clims_fut, file = "clims_ann_o2_rcp85_MRI-NEMURO_2100-2081.Rdata")

# Maps
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = na.omit(clims_base) ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2012-2031)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_MRI-NEMURO_2012-2031.pdf", sep = ""), height = 3, width = 6, dpi = 300)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = na.omit(clims_fut) ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2081-2100)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_MRI-NEMURO_2081-2100.pdf", sep = ""), height = 3, width = 6, dpi = 300)



### -----------------------------------------------------------------

setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/oxygen_files_separate_esm")
Wd <- getwd()
setwd(paste(Wd,"/","CNRM-PISCES", sep = ""))
# There are 10 netCDF files, one for each decade ranging from 2006 to 2100
nc <- nc_open("o2_Oyr_CNRM-CM5_rcp85_r1i1p1_2006-2015_regrid.nc")
nc
names(nc$var) # o2

# Need to find the level corresponding to 175m depth
z <- ncvar_get(nc, "lev_bnds"); z
# OK, need level 16 ! using raster::brick()

### With mclapply, return all annual values at level 16 for each period
periods <- c("2006-2015","2016-2025","2026-2035","2036-2045","2046-2055","2056-2065","2066-2075","2076-2085","2086-2095","2096-2100")
require("parallel")
# p <- "2016-2025"
res <- mclapply(periods, function(p) {
            
            # Useless message
            message(paste("Extracting O2 yearly outputs for ",p, sep = ""))
            ras <- raster::brick(paste("o2_Oyr_CNRM-CM5_rcp85_r1i1p1_",p,"_regrid.nc", sep = ""), varname = "o2", level = 16)
            r <- as.data.frame(ras, xy = T)
            # summary(r)
            #map <- ggplot() + geom_point(aes(x = x, y = y), data = r[which(r$X2016.07.02 == 0),]) +
                #scale_colour_viridis(name = "Mean O2 (ml/l)" ) + coord_quickmap() + theme_bw()
        
            #ggsave(plot = map, filename = paste("test.pdf", sep = ""), dpi = 300, width = 6, height = 3)

            # Change colnames depdning on period
            if(p == "2006-2015") {
                colnames(r)[c(3:length(r))] <- c(2006:2015)
            } else if(p == "2016-2025") {
                colnames(r)[c(3:length(r))] <- c(2016:2025)
            } else if(p == "2026-2035") {
                colnames(r)[c(3:length(r))] <- c(2026:2035)
            } else if(p == "2036-2045") {
                colnames(r)[c(3:length(r))] <- c(2036:2045)
            } else if(p == "2046-2055") {
                colnames(r)[c(3:length(r))] <- c(2046:2055)
            } else if(p == "2056-2065") {
                colnames(r)[c(3:length(r))] <- c(2056:2065)
            } else if(p == "2066-2075") {
                colnames(r)[c(3:length(r))] <- c(2066:2075)
            } else if(p == "2076-2085") {
                colnames(r)[c(3:length(r))] <- c(2076:2085)
            } else if(p == "2086-2095") {
                colnames(r)[c(3:length(r))] <- c(2086:2095)
            } else if(p == "2096-2100") {
                colnames(r)[c(3:length(r))] <- c(2096:2100)
            } # eo if else loop
            
            ### 15/11/19: land cells are marked by values == 0 --> replace by NAs
            #require("dplyr")
            #r[,c(3:length(r))] <- na_if(r[,c(3:length(r))], 0)
            
            # Return
            if(p == "2006-2015") {
                return(r)
            } else {
                return(r[,c(3:length(r))])
            } # eo if else loop
            
    }, mc.cores = 15
    
) # eo mclappy
require("dplyr")
ddf <- dplyr::bind_cols(res)
dim(ddf)
colnames(ddf)
# summary(ddf)
# ddf[1:1000,"x"]
# Add an id 
ddf$id <- paste(ddf$x, ddf$y, sep = "_")

# Ok, select the columns you are not interested in
baseline <- ddf[,c(1,2,98,9:28)]
fut <- ddf[,c(1,2,98,78:97)]
dim(baseline); dim(fut)
head(baseline); head(fut)

# Melt, and no need to convert to SST !
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("year","dO2")
colnames(m_fut)[c(4,5)] <- c("year","dO2")
summary(m_base)
summary(m_fut) 
# Add a column with just months to compute the 12 monthly climatologies with dplyr 
m_base$year <- data.frame(do.call(rbind, strsplit(as.character(m_base$year), "_") ))[,1]
m_fut$year <- data.frame(do.call(rbind, strsplit(as.character(m_fut$year), "_") ))[,1]
head(m_base); head(m_fut)

# Clean stuff
rm(fut, baseline, ddf); gc()

# Compute the annual climatology
require("dplyr")
clims_base <- data.frame( m_base %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2, na.rm=T) ) )
clims_fut <- data.frame( m_fut %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2, na.rm=T) ) )
head(clims_base); head(clims_fut); dim(clims_base); dim(clims_fut)
### Convert mol/m^3 to ml/l
# O2 (mol/kg OR mol/L): multiply by 1.029 (seawater density) and divide by 44.661e-6 (conversion from moles to L by the molar volume of O2)
clims_base$dO2 <- (clims_base$dO2)*0.001
clims_fut$dO2 <- (clims_fut$dO2)*0.001
clims_base$dO2 <- (clims_base$dO2)*1.029
clims_fut$dO2 <- (clims_fut$dO2)*1.029
clims_base$dO2 <- (clims_base$dO2)/44.661e-6
clims_fut$dO2 <- (clims_fut$dO2)/44.661e-6
summary(clims_base); summary(clims_fut)
head(clims_base); head(clims_fut)

### Ok, final step is to identify land grid cells and make sure the coords fit the other monthly climatologies
# To do so simply use a standard SST climatology from "byLotte"
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_forecast_byLotte")
grid <- raster::raster("sst1.nc", varname = "DATA_pisces")
test <- as.data.frame(grid, xy = T)
summary(test) 
# where test$DATA_pisces == NaN --> land
unique(test$x)
unique(clims_base$x)
# Correct x 
clims_base$x <- (clims_base$x)+0.5
clims_fut$x <- (clims_fut$x)+0.5
clims_base$id <- factor(paste(clims_base$x, clims_base$y, sep = "_"))
clims_fut$id <- factor(paste(clims_fut$x, clims_fut$y, sep = "_"))
clims_base <- clims_base[order(clims_base$id),]
clims_fut <- clims_fut[order(clims_fut$id),]
# Do the same with test
test$id <- factor(paste(test$x, test$y, sep = "_"))
test <- test[order(test$id),]
setdiff(test$id, clims_fut$id) # ok
setdiff(test$id, clims_base$id) # ok 
# Identify land cells and attribute NaNs
test$land <- NA
test[is.na(test$DATA_pisces),"land"] <- "yes"
test[!is.na(test$DATA_pisces),"land"] <- "no"
# Use this and ids levels to convert land o2 to NaN
landcells <- test[test$land == "yes","id"] # length(landcells)
# unique(clims_base[clims_base$id %in% landcells,"dO2"])
clims_base[clims_base$id %in% landcells,"dO2"] <- NA
clims_fut[clims_fut$id %in% landcells,"dO2"] <- NA
clims_base[clims_base$dO2 == 0.000000 & !is.na(clims_base$dO2),"dO2"] <- NA
clims_fut[clims_fut$dO2 == 0.0 & !is.na(clims_fut$dO2),"dO2"] <- NA
summary(clims_base)
summary(clims_fut)

# Maps
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = clims_base ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2012-2031)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_CNRM-PISCES_2012-2031.pdf", sep = ""), height = 3, width = 6, dpi = 300)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = clims_fut ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2081-2100)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_CNRM-PISCES_2081-2100.pdf", sep = ""), height = 3, width = 6, dpi = 300)

### OK, save 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(clims_base, file = "clims_ann_o2_rcp85_CNRM-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_ann_o2_rcp85_CNRM-PISCES_2100-2081.Rdata")



### -----------------------------------------------------------------

setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/oxygen_files_separate_esm")
Wd <- getwd()

### 4.X) IPSL-PISCES
setwd(paste(Wd,"/","IPSL-PISCES", sep = ""))
dir()
nc <- nc_open("o2_Oyr_IPSL-CM5A-LR_rcp85_r1i1p1_2006-2105_regridded.nc")
nc
# O2 units = mol m-3
# Time bins = 100 years (Jul 2006 -> Jul 2006)
# 31 depth levels
#names(nc$var) # o2
z <- ncvar_get(nc, "lev_bnds"); z
# OK, need level 16 ! using raster::brick()

setwd(paste(Wd,"/","IPSL-PISCES", sep = ""))
ras <- raster::brick("o2_Oyr_IPSL-CM5A-LR_rcp85_r1i1p1_2006-2105_regridded.nc", varname = "o2", level = 16)
ras # non conform x*y coordinates but OK for now. 100 time bands
r <- as.data.frame(ras, xy = T)
dim(r)
# Ok, now add the years as colnames: vector of c(2006:2105); length(c(2006:2105))
colnames(r)[c(3:length(r))] <- c(2006:2105)
# Now, compute annual mean over same period of interest as the other variables (2012-2031 & 2081-2100)
baseline <- r[,c(1,2,9:28)]
fut <- r[,c(1,2,78:97)]
dim(baseline); dim(fut)
# Add an id 
baseline$id <- factor(paste(baseline$x, baseline$y, sep = "_"))
fut$id <- factor(paste(fut$x, fut$y, sep = "_"))
# head(baseline$id); head(fut$id); tail(baseline$id); tail(fut$id)
### IDs don't seem to match previous IDS...you can still compute annual averages and then re-grid though
# Compute the 12 monthly clims after melting
require("dplyr","reshape2")
m_base <- melt(baseline, id.vars = c("id","x","y"))
m_fut <- melt(fut, id.vars = c("id","x","y"))
colnames(m_base)[c(4,5)] <- c("year","dO2")
colnames(m_fut)[c(4,5)] <- c("year","dO2")
# Clean some stuff
rm(fut, baseline, r); gc()
# Compute annual mean
clims_base <- data.frame( m_base %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
clims_fut <- data.frame( m_fut %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), dO2 = mean(dO2) ) )
summary(clims_base); summary(clims_fut)
dim(clims_base); dim(clims_fut)
# head(clims_base$id); head(clims_fut$id); tail(clims_base$id); tail(clims_fut$id)

### Now, need to convert mol/m3 to ml/l. You know how to convert from mol/L to ml/l so just start converting mol/m3 to mol/L
clims_base$dO2 <- (clims_base$dO2)*0.001
clims_fut$dO2 <- (clims_fut$dO2)*0.001
# convert mol/L to ml/l
clims_base$dO2 <- (clims_base$dO2)*1.029
clims_fut$dO2 <- (clims_fut$dO2)*1.029
clims_base$dO2 <- (clims_base$dO2)/44.661e-6
clims_fut$dO2 <- (clims_fut$dO2)/44.661e-6
#
summary(clims_base); summary(clims_fut)
dim(clims_base); dim(clims_fut)

### Humm...might need to roatate the x coords a lil bit, compare the x values to the other netCDF outputs
#setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_forecast_byLotte")
#grid <- raster::raster("sst1.nc", varname = "DATA_pisces")
#test <- as.data.frame(grid, xy = T)[,c("x","y")]
#summary(test$x) # Reference
#summary(clims_base$x) # ahhh missing 0.5 ;)
#unique(test$x)
#unique(clims_base$x)
# OK, add 0.5, give new ID, and re-order priori to saving !
# Check latitudes too
#unique(test$y)
#unique(clims_base$y)
# Lats are fine
clims_base$x <- (clims_base$x)+0.5
clims_fut$x <- (clims_fut$x)+0.5
clims_base$id <- factor(paste(clims_base$x, clims_base$y, sep = "_"))
clims_fut$id <- factor(paste(clims_fut$x, clims_fut$y, sep = "_"))
clims_base <- clims_base[order(clims_base$id),]
clims_fut <- clims_fut[order(clims_fut$id),]

# Maps
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = clims_base ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2012-2031)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_IPSL-PISCES_2012-2031.pdf", sep = ""), height = 3, width = 6, dpi = 300)

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = dO2), data = clims_fut ) +
        scale_fill_viridis(name = "Mean O2 (ml/l)\n(2081-2100)" ) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_ann_o2_rcp85_IPSL-PISCES_2081-2100.pdf", sep = ""), height = 3, width = 6, dpi = 300)


### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
save(clims_base, file = "clims_ann_o2_rcp85_IPSL-PISCES_2031-2012.Rdata")
save(clims_fut, file = "clims_ann_o2_rcp85_IPSL-PISCES_2100-2081.Rdata")


### -----------------------------------------------------------------

### 15/11/19: While waiting for help to disentangle the grid of CNRM-PISCES and IPSL-PISCES, compute the monthly deltas of dO2 for the 3 other ESMs using the anoamies to the annual average from the in situ climatologies :)

### A) First, like for O2, get in situ climatologies of O2 @ 175m depth
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
# m <- "apr"
res <- lapply(months, function(m) {
            d <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")[,c("x","y","dO2")]
            d$id <- paste(d$x, d$y, sep = "_")
            d <- d[order(d$id),]
            d$month <- m
            return(d)
    } # eo FUN
) # eo lapply
# Cbind
ddf <- do.call(rbind, res)
head(ddf)
# Dcast to put monthly means as columns
clims <- dcast(data = ddf, id + x + y ~ month, value.var = "dO2")
dim(clims); summary(clims)

# Compute annual average 
require("matrixStats")
clims$ann <- rowMeans(as.matrix(clims[,c(4:15)]), na.rm = T)
summary(clims$ann) # summary(clims$Ann)
### Change colnames and compute difference to annual mean
colnames(clims)[c(4:16)] <- c("Apr","Aug","Dec","Feb","Jan","Jul","Jun","Mar","May","Nov","Oct","Sep","Ann")
anoms <- data.frame(id = clims$id, x = clims$x, y = clims$y, 
            Jan = (clims$Jan)-(clims$Ann), Feb = (clims$Feb)-(clims$Ann), Mar = (clims$Mar)-(clims$Ann),
            Apr = (clims$Apr)-(clims$Ann), May = (clims$May)-(clims$Ann), Jun = (clims$Jun)-(clims$Ann),
            Jul = (clims$Jul)-(clims$Ann), Aug = (clims$Aug)-(clims$Ann), Sep = (clims$Sep)-(clims$Ann),
            Oct = (clims$Oct)-(clims$Ann), Nov = (clims$Nov)-(clims$Ann), Dec = (clims$Dec)-(clims$Ann)
) # eo ddf
# Rotate x coordinates
anoms$x2 <- anoms$x 
anoms[anoms$x < 0 ,"x2"] <- (anoms[anoms$x < 0 ,"x"]) + 360
# And modify cell id accordingly
anoms$id <- paste(anoms$x2, anoms$y, sep = "_")
summary(anoms)

### Print maps of monthly anomalies to the annual mean 
#months <- c("Apr","Aug","Dec","Feb","Jan","Jul","Jun","Mar","May","Nov","Oct","Sep")
#for(momo in months) {

    d <- anoms[,c("x2","y",momo)]
    colnames(d)[c(3)] <- "delta"

    map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = delta), data = d) +
        scale_fill_gradient2(name = "O2 anomaly\n(ml/l)", low = "#3288bd", high = "#d53e4f", mid = "white") +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    ggsave(plot = map, filename = paste("map_mon_anom_dO2_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
    rm(d);gc()
} # eo for loop

### Ok nice, save anoms before deriving baseline and future seasonal cycle
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
#anoms <- anoms[order(anoms$id),]
# save(x = anoms, file = "clims_mon+ann_anoms_o2_insitu.Rdata")
anoms <- get(load("clims_mon+ann_anoms_o2_insitu.Rdata"))
dim(anoms); summary(anoms)

### B) Get model-based dO2 climatologies and derive the monthly variations from the observed monthly anomalies
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
ESMs <- c("IPSL-PISCES","CNRM-PISCES","CESM-BEC") # "GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")
# For every ESM, get the annual baseline and future climatologies, and derive monthly ones from 'anoms'
#esm <- "ISPL-PISCES" ; esm <- "CESM-BEC"
for(esm in ESMs) {
    
        # Get annual data
        message(paste("Preparing baseline and future monthly O2 for ", esm, sep = ""))
        message(paste("", sep = ""))
        clim_base <- get(load(paste("clims_ann_o2_rcp85_",esm,"_2031-2012.Rdata", sep = "")))
        clim_fut <- get(load(paste("clims_ann_o2_rcp85_",esm,"_2100-2081.Rdata", sep = "")))
        # head(clim_base); head(clim_fut); head(anoms)
        # tail(clim_base); tail(clim_fut); tail(anoms)
        # dim(clim_base); dim(anoms)# 
        
        # Derive monthly clims from anoms for baseline
        clim_base$Jan <- (anoms$Jan) + (clim_base$dO2)
        clim_base$Feb <- (anoms$Feb) + (clim_base$dO2)
        clim_base$Mar <- (anoms$Mar) + (clim_base$dO2)
        clim_base$Apr <- (anoms$Apr) + (clim_base$dO2)
        clim_base$May <- (anoms$May) + (clim_base$dO2)
        clim_base$Jun <- (anoms$Jun) + (clim_base$dO2)
        clim_base$Jul <- (anoms$Jul) + (clim_base$dO2)
        clim_base$Aug <- (anoms$Aug) + (clim_base$dO2)
        clim_base$Sep <- (anoms$Sep) + (clim_base$dO2)
        clim_base$Oct <- (anoms$Oct) + (clim_base$dO2)
        clim_base$Nov <- (anoms$Nov) + (clim_base$dO2)
        clim_base$Dec <- (anoms$Dec) + (clim_base$dO2)
        # and now for future
        clim_fut$Jan <- (anoms$Jan) + (clim_fut$dO2)
        clim_fut$Feb <- (anoms$Feb) + (clim_fut$dO2)
        clim_fut$Mar <- (anoms$Mar) + (clim_fut$dO2)
        clim_fut$Apr <- (anoms$Apr) + (clim_fut$dO2)
        clim_fut$May <- (anoms$May) + (clim_fut$dO2)
        clim_fut$Jun <- (anoms$Jun) + (clim_fut$dO2)
        clim_fut$Jul <- (anoms$Jul) + (clim_fut$dO2)
        clim_fut$Aug <- (anoms$Aug) + (clim_fut$dO2)
        clim_fut$Sep <- (anoms$Sep) + (clim_fut$dO2)
        clim_fut$Oct <- (anoms$Oct) + (clim_fut$dO2)
        clim_fut$Nov <- (anoms$Nov) + (clim_fut$dO2)
        clim_fut$Dec <- (anoms$Dec) + (clim_fut$dO2)
        
        # summary(clim_base); summary(clim_fut)
        # nrow(na.omit(clim_base[clim_base$Aug < 0,]))
        # Replace all negative and zero values values to the lowest value that is still above zero
        val1 <- min(clim_base[,c(5:16)][clim_base[,c(5:16)] > 0], na.rm = T)
        val2 <- min(clim_fut[,c(5:16)][clim_fut[,c(5:16)] > 0], na.rm = T)
        clim_base[,c(5:16)][clim_base[,c(5:16)] <= 0 & !is.na(clim_base[,c(5:16)])] <- val1
        clim_fut[,c(5:16)][clim_fut[,c(5:16)] <= 0 & !is.na(clim_fut[,c(5:16)])] <- val2
        # summary(clim_base); summary(clim_fut)
        colnames(clim_base)[4] <- "Annual"
        colnames(clim_fut)[4] <- "Annual"
        # dim(clim_base); dim(clim_fut)
        # Compute delta between present and future
        # Compute the 12 monthly deltas SST and map
        # summary(clim_base$Apr) ; summary(clim_fut$Apr)
        deltas <- data.frame(id = clim_base$id, x = clim_base$x, y = clim_base$y, 
        	Jan = (clim_fut$Jan) - (clim_base$Jan), Feb = (clim_fut$Feb) - (clim_base$Feb), Mar = (clim_fut$Mar) - (clim_base$Mar), 
        	Apr = (clim_fut$Apr) - (clim_base$Apr), May = (clim_fut$May) - (clim_base$May), Jun = (clim_fut$Jun) - (clim_base$Jun),
        	Jul = (clim_fut$Jul) - (clim_base$Jul), Aug = (clim_fut$Aug) - (clim_base$Aug), Sep = (clim_fut$Sep) - (clim_base$Sep), 
        	Oct = (clim_fut$Oct) - (clim_base$Oct), Nov = (clim_fut$Nov) - (clim_base$Nov), Dec = (clim_fut$Dec) - (clim_base$Dec)
        ) # eo ddf
        # summary(deltas)
        # deltas[c(1:200),c(4:6)]

        ### Save on your dir and plot from your computer because something is obviously wrong with the graphics device 
        setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims/")
        save(deltas, file = paste("clims_mon_diff_o2_rcp85_",esm,"_2031-2100.Rdata", sep = "") )
        save(clim_base, file = paste("clims_mon_o2_rcp85_",esm,"_2031-2012.Rdata", sep = "") )
        save(clim_fut, file = paste("clims_mon_o2_rcp85_",esm,"_2100-2081.Rdata", sep = "") )

        # And print maps
        for(momo in months) {

            d <- deltas[,c("x","y",momo)]
            colnames(d)[c(3)] <- "delta"

            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = delta), data = na.omit(d) ) +
                scale_fill_gradient2(name = "Delta O2 (ml/l)\n(2100-2031)", low = "#3288bd", high = "#d53e4f", mid = "white") +
                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

            ggsave(plot = map, filename = paste("map_delta_o2_rcp85_",esm,"_",momo,".pdf", sep = ""), height = 3, width = 6, dpi = 300)
            rm(d);gc()
            
        } # eo for loop
    
} # eo for loop



### ----------------------------------------------------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------------------------------------------------


