
# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("RColorBrewer")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("geosphere")
library("rgdal")
library("raster")
library("betapart")

world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------

### Set the working directories, vectors etc.
WD <- getwd()
rcp <- "rcp85"
# Vectors
SDMs <- c('GAM','GLM','RF','ANN')
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
# Vector of earth system models
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Load individual projections in perc change in SR for phyto and zoo
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("table_ann_changes_",dir())]#; files
# Separate the ones for richness (no "beta.div") from those with beta.div
files2 <- files[grep("beta.div", files)]#; files2
files1 <- files[!(files %in% files2)]; files1

require("parallel")
res <- mclapply(files1, function(f) {
            # Useless message
            message(paste("Loading results from ",f, sep = ""))
            t <- get(load(f))
            # Get attributes from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            t$ESM <- terms[4,1] 
            t$SDM <- terms[5,1] 
            p <- terms[6,1]
            t$pool <- str_replace(as.character(p),".Rdata","")
            # Return
            return(t)
    }, mc.cores = 25
) # eo mclapply
# Rbind
table <- dplyr::bind_rows(res)
rm(res); gc()

### Computing ensembles (all, SDM/ESM/pool)
ens <- data.frame(table %>%
        group_by(id) %>%
        summarize(x = unique(x), y = unique(y), 
            rich_tot = mean(perc_rich_tot, na.rm = T), sd_rich_tot = sd(perc_rich_tot, na.rm = T),
            rich_phyto = mean(perc_rich_phyto, na.rm = T), sd_rich_phyto = sd(perc_rich_phyto, na.rm = T),
            rich_zoo = mean(perc_rich_zoo, na.rm = T), sd_rich_zoo = sd(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
# summary(ens)


### 2°) First, map stdev (Figs. A and B of the overall panel)
mapA <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd_rich_phyto), data = ens[!is.na(ens$sd_rich_phyto),]) +
    geom_contour(colour = "grey25", binwidth = 33, size = 0.25, aes(x = x, y = y, z = sd_rich_phyto),
                    data = ens[!is.na(ens$sd_rich_phyto),]) +
 	scale_fill_distiller(name = "Standard\ndeviation", palette = "YlOrRd", direction = 1) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapB <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd_rich_zoo), data = ens[!is.na(ens$sd_rich_zoo),]) +
    geom_contour(colour = "grey25", binwidth = 33, size = 0.25, aes(x = x, y = y, z = sd_rich_zoo),
                    data = ens[!is.na(ens$sd_rich_zoo),]) +
 	scale_fill_distiller(name = "Standard\ndeviation", palette = "YlOrRd", direction = 1) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

require("ggpubr")
ggarrange(mapA, mapB, ncol = 2, nrow = 1, labels = c("A1","A2"))    
      
    
### 3°) The 8 maps (LETTERS[3:10]) for the each SDM type 
# Per SDM (for facet)
ens.SDM <- data.frame(table %>%
        group_by(id,SDM) %>%
        summarize(x = unique(x), y = unique(y), 
            rich_tot = mean(perc_rich_tot, na.rm = T),
            rich_phyto = mean(perc_rich_phyto, na.rm = T),
            rich_zoo = mean(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
summary(ens.SDM) 
min <- floor(min(ens.SDM[,c("rich_phyto","rich_zoo")], na.rm = T) )
min

mapC <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "GLM"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_phyto >= 55 & ens.SDM$SDM == "GLM"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "GLM"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapD <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "GLM"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_zoo >= 55 & ens.SDM$SDM == "GLM"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "GLM"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapE <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "GAM"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_phyto >= 55 & ens.SDM$SDM == "GAM"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "GAM"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapF <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "GAM"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_zoo >= 55 & ens.SDM$SDM == "GAM"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "GAM"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapG <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "ANN"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_phyto >= 55 & ens.SDM$SDM == "ANN"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "ANN"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapH <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "ANN"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_zoo >= 55 & ens.SDM$SDM == "ANN"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "ANN"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )       
#
mapI <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "RF"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_phyto >= 55 & ens.SDM$SDM == "RF"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.SDM[which(ens.SDM$rich_phyto < 55 & ens.SDM$SDM == "RF"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapJ <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "RF"),]) +
    geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_zoo >= 55 & ens.SDM$SDM == "RF"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.SDM[which(ens.SDM$rich_zoo < 55 & ens.SDM$SDM == "RF"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

### Test panel
labels <- paste(expand.grid(c(1,2), LETTERS[1:10])$Var2, expand.grid(c(1,2), LETTERS[1:10])$Var1, sep = "")
panel1 <- ggarrange(mapA, mapB, mapC, mapD, mapE, mapF, mapG, mapH, mapI, mapJ, ncol = 2, nrow = 5, labels = labels)    
# Cool

### Same, with ESM
ens.ESM <- data.frame(table %>%
        group_by(id,ESM) %>%
        summarize(x = unique(x), y = unique(y), 
        rich_tot = mean(perc_rich_tot, na.rm = T),
        rich_phyto = mean(perc_rich_phyto, na.rm = T),
        rich_zoo = mean(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
head(ens.ESM)

min <- floor(min(ens.ESM[,c("rich_phyto","rich_zoo")], na.rm = T) )
min

unique(ens.ESM$ESM)

mapC <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "CESM-BEC"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_phyto >= 55 & ens.ESM$ESM == "CESM-BEC"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "CESM-BEC"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapD <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "CESM-BEC"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_zoo >= 55 & ens.ESM$ESM == "CESM-BEC"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "CESM-BEC"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapE <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "CNRM-PISCES"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_phyto >= 55 & ens.ESM$ESM == "CNRM-PISCES"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "CNRM-PISCES"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapF <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "CNRM-PISCES"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_zoo >= 55 & ens.ESM$ESM == "CNRM-PISCES"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "CNRM-PISCES"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapG <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "GFDL-TOPAZ"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_phyto >= 55 & ens.ESM$ESM == "GFDL-TOPAZ"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "GFDL-TOPAZ"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapH <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "GFDL-TOPAZ"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_zoo >= 55 & ens.ESM$ESM == "GFDL-TOPAZ"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "GFDL-TOPAZ"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )       
#
mapI <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "IPSL-PISCES"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_phyto >= 55 & ens.ESM$ESM == "IPSL-PISCES"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "IPSL-PISCES"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapJ <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "IPSL-PISCES"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_zoo >= 55 & ens.ESM$ESM == "IPSL-PISCES"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "IPSL-PISCES"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapK <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
        data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "MRI-NEMURO"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_phyto >= 55 & ens.ESM$ESM == "MRI-NEMURO"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
            data = ens.ESM[which(ens.ESM$rich_phyto < 55 & ens.ESM$ESM == "MRI-NEMURO"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

mapL <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "MRI-NEMURO"),]) +
    geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_zoo >= 55 & ens.ESM$ESM == "MRI-NEMURO"),], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
            data = ens.ESM[which(ens.ESM$rich_zoo < 55 & ens.ESM$ESM == "MRI-NEMURO"),] ) +
 	scale_fill_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,55)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )   
              
### Arrange in second panel
labels <- paste(expand.grid(c(1,2), LETTERS[6:10])$Var2, expand.grid(c(1,2), LETTERS[6:10])$Var1, sep = "")
panel2 <- ggarrange(mapC, mapD, mapE, mapF, mapG, mapH, mapI, mapJ, mapK, mapL, ncol = 2, nrow = 5, labels = labels)    

### Save both panels
setwd("/net/kryo/work/fabioben/OVERSEE/data/")
ggsave(plot = panel1, filename = "panel_maps_uncertainty_A-E.jpg", dpi = 300, width = 10, height = 10)
ggsave(plot = panel2, filename = "panel_maps_uncertainty_F-J.jpg", dpi = 300, width = 10, height = 10)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
 