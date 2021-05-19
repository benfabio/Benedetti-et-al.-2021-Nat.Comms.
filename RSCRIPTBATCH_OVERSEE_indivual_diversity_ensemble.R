
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

### 1°) Examine ensembles of baseline annual species richness
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("div_baseline_",dir())]; files

require("parallel")
# f <- files[1]
res <- mclapply(files, function(f) {
            # Useless message
            message(paste("Loading results from ",f, sep = ""))
            t <- read.table(f, sep = "\t", h = T) # dim(t); colnames(t)
            # Get attributes from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            group <- terms[2,1] 
            sdm <- terms[7,1] 
            p <- terms[6,1]
            m <- terms[8,1]
            mm <- str_replace(as.character(m),".txt","")
            ddf <- data.frame(id = t$cell_id, x = t$x, y = t$y, 
                         month = mm, group = group,
                         sdm = sdm, pool = p, rich = t$rich
            ) # eo ddf
            # Return
            return(ddf)
    }, mc.cores = 35
) # eo mclapply
# Rbind
table.base <- dplyr::bind_rows(res)
dim(table.base)
head(table.base)
rm(res); gc()
table.base <- table.base[order(table.base$id),]

### Computing ensembles (all, SDM/ESM/pool)
ens <- data.frame(table.base %>% group_by(id,group) %>%
        summarize(x = unique(x), y = unique(y), mean.rich = mean(rich, na.rm = T), sd.rich = sd(rich, na.rm = T) ) 
) # ddf
summary(ens)
head(ens)

# Dcast to separate between 
dens <- dcast(ens, id + x + y ~ group, value.var = c("mean.rich") )
dim(dens)
summary(dens)
dens$tot <- (dens$phyto)+(dens$zoo)
# Define limts of the color scale
#min <- floor(min(c(dens$phyto,dens$zoo), na.rm = T))
#max <- floor(max(c(dens$phyto,dens$zoo), na.rm = T))

### 23/06/2020: As Niki suggested, rather plot % of species modelled 
dens$phyto2 <- (dens$phyto)/336
dens$zoo2 <- (dens$zoo)/524
dens$tot2 <- (dens$tot)/860
summary(dens)
min <- 0
max <- 0.7
# min;max
map1 <- ggplot() + geom_tile(aes(x = x, y = y, fill = tot2), data = na.omit(dens)) +
	#scale_fill_viridis(name = "Mean annual\nrichness", begin = 0, end = 1, limits = c(min,max)) + 
    scale_fill_viridis(name = "", rescaler = function(x, to = c(0,1), from = NULL) {
        ifelse(x < 0.4, scales::rescale(x, to = to, from = c(min(x, na.rm = T), 0.4)), 1) } ) +
    geom_contour(colour = "grey75", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = tot2), data = na.omit(dens)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

map2 <- ggplot() + geom_tile(aes(x = x, y = y, fill = phyto2), data = na.omit(dens)) +
	#scale_fill_viridis(name = "Mean annual\nrichness", limits = c(min,max) ) +
    scale_fill_viridis(name = "", rescaler = function(x, to = c(0,1), from = NULL) {
        ifelse(x < 0.4, scales::rescale(x, to = to, from = c(min(x, na.rm = T), 0.4)), 1) } ) +
    geom_contour(colour = "grey75", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = phyto2), data = na.omit(dens)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

map3 <- ggplot() + geom_tile(aes(x = x, y = y, fill = zoo2), data = na.omit(dens)) +
	#scale_fill_viridis(name = "Mean annual\nrichness", limits = c(min,max) ) +
    scale_fill_viridis(name = "", rescaler = function(x, to = c(0,1), from = NULL) {
        ifelse(x < 0.4, scales::rescale(x, to = to, from = c(min(x, na.rm = T), 0.4)), 1) } ) +
    geom_contour(colour = "grey75", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = zoo2), data = na.omit(dens)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

setwd(WD)

ggsave(plot = map2, filename = paste("map_ensemble_annual_rich_phyto_baseline.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map3, filename = paste("map_ensemble_annual_rich_zoo_baseline.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map1, filename = paste("map_ensemble_annual_rich_tot_baseline.jpg", sep = ""), dpi = 300, width = 7, height = 5)	

### 09/08/2020: and the true postscript versions for Niki
grDevices::cairo_ps(filename = "map_ensemble_annual_rich_phyto_baseline.ps", width = 7, height = 5, fallback_resolution = 300)
print(map2)
dev.off()
grDevices::cairo_ps(filename = "map_ensemble_annual_rich_zoo_baseline.ps", width = 7, height = 5, fallback_resolution = 300)
print(map3)
dev.off()
grDevices::cairo_ps(filename = "map_ensemble_annual_rich_tot_baseline.ps", width = 7, height = 5, fallback_resolution = 300)
print(map1)
dev.off()


### 23/06/2020: And plot zonal patterns, one per SDM for baseline richness
head(table.base)
zonal <- data.frame(table.base %>%
             group_by(group,sdm,y) %>%
             summarise(rich = mean(rich,na.rm=T) )
) # eo ddf
summary(zonal)

plot2 <- ggplot() + 
   geom_path(aes(y = y, x = rich/336), data = na.omit(zonal[zonal$group == "phyto" & zonal$sdm == "GLM",]), colour = "black", linetype = "solid") +
   geom_path(aes(y = y, x = rich/336), data = na.omit(zonal[zonal$group == "phyto" & zonal$sdm == "GAM",]), colour = "black", linetype = "longdash") +
   geom_path(aes(y = y, x = rich/336), data = na.omit(zonal[zonal$group == "phyto" & zonal$sdm == "ANN",]), colour = "black", linetype = "dashed") +
   geom_path(aes(y = y, x = (rich/336)*3.2), data = na.omit(zonal[zonal$group == "phyto" & zonal$sdm == "RF",]), colour = "grey50", linetype = "solid") +
   scale_x_continuous(name = "", limits = c(0,0.55), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6)) + 
   scale_y_continuous(position = "left", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
   theme_classic() 

plot3 <- ggplot() + 
   geom_path(aes(y = y, x = rich/524), data = na.omit(zonal[zonal$group == "zoo" & zonal$sdm == "GLM",]), colour = "black", linetype = "solid") +
   geom_path(aes(y = y, x = rich/524), data = na.omit(zonal[zonal$group == "zoo" & zonal$sdm == "GAM",]), colour = "black", linetype = "longdash") +
   geom_path(aes(y = y, x = rich/524), data = na.omit(zonal[zonal$group == "zoo" & zonal$sdm == "ANN",]), colour = "black", linetype = "dashed") +
   geom_path(aes(y = y, x = (rich/524)*3.2), data = na.omit(zonal[zonal$group == "zoo" & zonal$sdm == "RF",]), colour = "grey50", linetype = "solid") +
   scale_y_continuous(position = "left", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
   scale_x_continuous(name = "", limits = c(0,0.55), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6)) + 
   theme_classic() 

# Dcast to separate between 
dzon <- dcast(zonal, sdm + y ~ group, value.var = c("rich") )
dim(dzon) ; summary(dzon)
dzon$tot <- (dzon$phyto)+(dzon$zoo)
   
plot1 <- ggplot() + 
   geom_path(aes(y = y, x = tot/860), data = na.omit(dzon[dzon$sdm == "GLM",]), colour = "black", linetype = "solid") +
   geom_path(aes(y = y, x = tot/860), data = na.omit(dzon[dzon$sdm == "GAM",]), colour = "black", linetype = "longdash") +
   geom_path(aes(y = y, x = tot/860), data = na.omit(dzon[dzon$sdm == "ANN",]), colour = "black", linetype = "dashed") +
   geom_path(aes(y = y, x = (tot/860)*3.2), data = na.omit(dzon[dzon$sdm == "RF",]), colour = "grey50", linetype = "solid") +
   scale_y_continuous(position = "left", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
   scale_x_continuous(name = "", limits = c(0,0.55), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6)) + 
   theme_classic() 


# zonal <- data.frame(dens %>%
#             group_by(y) %>%
#             summarise(mean_plankton = mean(tot, na.rm=T), sd_plankton = sd(tot, na.rm=T),
#             mean_phyto = mean(phyto, na.rm=T), sd_phyto = sd(phyto, na.rm=T),
#             mean_zoo = mean(zoo, na.rm=T), sd_zoo = sd(zoo, na.rm=T) )
# ) # eo ddf
# plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_plankton - sd_plankton, xmax = mean_plankton + sd_plankton),
#                 fill = "grey70", data = na.omit(zonal)) +
#    geom_path(aes(y = y, x = mean_plankton), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(0,290)) +
#    xlab("Mean annual richness") + scale_y_continuous(position = "left", name = "Latitude (°)") +
#    theme_classic()
#
# plot2 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_phyto - sd_phyto, xmax = mean_phyto + sd_phyto),
#                 fill = "grey70", data = na.omit(zonal)) +
#    geom_path(aes(y = y, x = mean_phyto), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(0,200)) +
#    xlab("Mean annual richness") + scale_y_continuous(position = "left", name = "Latitude (°)") +
#    theme_classic()
# #
# plot3 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_zoo - sd_zoo, xmax = mean_zoo + sd_zoo),
#                 fill = "grey70", data = na.omit(zonal)) +
#    geom_path(aes(y = y, x = mean_zoo), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(0,200)) +
#    xlab("Mean annual richness") + scale_y_continuous(position = "left", name = "Latitude (°)") +
#    theme_classic()
#
ggsave(plot = plot1, filename = "plot_zonal_ann_plankton_rich_baseline_ensemble.jpg", dpi = 300, width = 3, height = 4)
ggsave(plot = plot2, filename = "plot_zonal_ann_phyto_rich_baseline_ensemble.jpg", dpi = 300, width = 3, height = 4)
ggsave(plot = plot3, filename = "plot_zonal_ann_zoo_rich_baseline_ensemble.jpg", dpi = 300, width = 3, height = 4)

grDevices::cairo_ps(filename = "plot_zonal_ann_plankton_rich_baseline_ensemble.ps", width = 3, height = 4, fallback_resolution = 300)
print(plot1)
dev.off()
grDevices::cairo_ps(filename = "plot_zonal_ann_phyto_rich_baseline_ensemble.ps", width = 3, height = 4, fallback_resolution = 300)
print(plot2)
dev.off()
grDevices::cairo_ps(filename = "plot_zonal_ann_zoo_rich_baseline_ensemble.ps", width = 3, height = 4, fallback_resolution = 300)
print(plot3)
dev.off()


### 20/04/2020: Use grid arrange or ggExtra to arrange the plots in panels that will make them 
library("ggpubr")
# ?ggarrange
# ggarrange(plot1,map1, labels = c("A","B"), ncol = 2, nrow = 1, widths = c(1,2.7), align = "v")
# ggarrange(plot2,map2, labels = c("C","D"), ncol = 2, nrow = 1, widths = c(1,2.7), align = "v")
# ggarrange(plot3,map3, labels = c("E","F"), ncol = 2, nrow = 1, widths = c(1,2.7), align = "v")
#ggarrange(plot1,map1,plot2,map2,plot3,map3, labels = c("A","B","C","D","E","F"), ncol = 2, nrow = 3, widths = c(1,2.7) )

setwd("/net/kryo/work/fabioben/OVERSEE/data")
panel <- ggarrange(plot1,map1,plot2,map2,plot3,map3, labels = c("A","B","C","D","E","F"), ncol = 2, nrow = 3, widths = c(1,2.7) )
ggsave(plot = panel, filename = "panel_ensemble_baseline_rich_Fig1A-F.pdf", dpi = 300, width = 8, height = 8)

grDevices::cairo_ps(filename = "panel_ensemble_baseline_rich_Fig1A-F.ps", width = 8, height = 8, fallback_resolution = 300)
print(panel)
dev.off()

### And save one with as .pdf with the color bar
map <- ggplot() + geom_tile(aes(x = x, y = y, fill = zoo2), data = na.omit(dens)) +
    scale_fill_viridis(name = "", rescaler = function(x, to = c(0,1), from = NULL) {
        ifelse(x < 0.4, scales::rescale(x, to = to, from = c(min(x, na.rm = T), 0.4)), 1) }, limits = c(0,0.7) ) +
    geom_contour(colour = "grey75", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = zoo2), data = na.omit(dens)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

ggsave(plot = map, filename = paste("map4colorbar1.pdf", sep = ""), dpi = 300, width = 8, height = 5)	


### And check ensemble per SDMs
ens.SDM <- data.frame(table %>% group_by(id,group,sdm) %>%
        summarize(x = unique(x), y = unique(y), mean.rich = mean(rich, na.rm = T), sd.rich = sd(rich, na.rm = T) ) 
) # ddf

# dim(ens.SDM[which(ens.SDM$sdm == "GAM" & dens$group == "phyto"),])
min1 <- min(ens.SDM[which(ens.SDM$group == "phyto"),"mean.rich"], na.rm = T)
max1 <- max(ens.SDM[which(ens.SDM$group == "phyto"),"mean.rich"], na.rm = T)
min2 <- min(ens.SDM[which(ens.SDM$group == "zoo"),"mean.rich"], na.rm = T)
max2 <- max(ens.SDM[which(ens.SDM$group == "zoo"),"mean.rich"], na.rm = T)

for(sdm in unique(ens.SDM$sdm) ) {
    
    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich),
                data = ens.SDM[which(ens.SDM$group == "phyto" & ens.SDM$sdm == sdm),] ) +
    	scale_fill_viridis(name = "Annual richness", limits = c(min1,max1)) +
        geom_contour(colour = "grey75", binwidth = 50, size = 0.25, aes(x = x, y = y, z = mean.rich),
                data = ens.SDM[which(ens.SDM$group == "phyto" & ens.SDM$group == sdm),] ) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    #
    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich),
                data = ens.SDM[which(ens.SDM$group == "zoo" & ens.SDM$sdm == sdm),] ) +
    	scale_fill_viridis(name = "Annual richness", limits = c(min2,max2)) +
        geom_contour(colour = "grey75", binwidth = 50, size = 0.25, aes(x = x, y = y, z = mean.rich),
                data = ens.SDM[which(ens.SDM$group == "zoo" & ens.SDM$group == sdm),] ) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    
    ggsave(plot = map1, filename = paste("map_ensemble_annual_rich_phyto_baseline_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	
    ggsave(plot = map2, filename = paste("map_ensemble_annual_rich_zoo_baseline_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	
            
} # eo sdm in SDMs

### Put all 4 SDM-scale zonal on the same plot? (with facet)
head(ens.SDM)
zonal <- data.frame(ens.SDM %>% 
        group_by(y,sdm,group) %>% 
        summarise(mean = mean(mean.rich, na.rm=T), sd = sd(mean.rich, na.rm=T) 
    )    
) # eo ddf
summary(zonal)

plot1 <- ggplot(data = zonal[zonal$group == "phyto",]) +
   geom_ribbon(aes(x = y, ymin = mean - sd, ymax = mean + sd), fill = "grey70" ) +
   geom_line(aes(x = y, y = mean), colour = "black" ) + 
   ylab("Mean annual species richness\n(Phytoplankton)") + xlab("Latitude (°)") + 
   theme_classic() + coord_flip() + facet_wrap(~factor(sdm), ncol = 2, scales = "fixed")

plot2 <- ggplot(data = zonal[zonal$group == "zoo",]) +
   geom_ribbon(aes(x = y, ymin = mean - sd, ymax = mean + sd), fill = "grey70" ) +
   geom_line(aes(x = y, y = mean), colour = "black" ) + 
   ylab("Mean annual species richness\n(Zooplankton)") + xlab("Latitude (°)") + 
   theme_classic() + coord_flip() + facet_wrap(~factor(sdm), ncol = 2, scales = "fixed")
#
ggsave(plot = plot1, filename = "plot_zonal_ann_phyto_rich_baseline_ensemble_SDMs.jpg", dpi = 300, width = 6, height = 6)
ggsave(plot = plot2, filename = "plot_zonal_ann_zoo_rich_baseline_ensemble_SDMs.jpg", dpi = 300, width = 6, height = 6)


### And check ensemble per pools
ens.pool <- data.frame(table %>% group_by(id,group,pool) %>%
        summarize(x = unique(x), y = unique(y), mean.rich = mean(rich, na.rm = T), sd.rich = sd(rich, na.rm = T) ) 
) # ddf
# summary(ens.pool)
min1 <- min(ens.pool[which(ens.pool$group == "phyto"),"mean.rich"], na.rm = T)
max1 <- max(ens.pool[which(ens.pool$group == "phyto"),"mean.rich"], na.rm = T)
min2 <- min(ens.pool[which(ens.pool$group == "zoo"),"mean.rich"], na.rm = T)
max2 <- max(ens.pool[which(ens.pool$group == "zoo"),"mean.rich"], na.rm = T)

for(p in unique(ens.pool$pool) ) {
    
    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich),
                data = ens.pool[which(ens.pool$group == "phyto" & ens.pool$pool == p),] ) +
    	scale_fill_viridis(name = "Annual richness", limits = c(min1,max1)) +
        geom_contour(colour = "grey75", binwidth = 50, size = 0.25, aes(x = x, y = y, z = mean.rich),
                data = ens.pool[which(ens.pool$group == "phyto" & ens.pool$pool == p),] ) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    #
    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich),
                data = ens.pool[which(ens.pool$group == "zoo" & ens.pool$pool == p),] ) +
    	scale_fill_viridis(name = "Annual richness", limits = c(min2,max2)) +
        geom_contour(colour = "grey75", binwidth = 50, size = 0.25, aes(x = x, y = y, z = mean.rich),
                data = ens.pool[which(ens.pool$group == "zoo" & ens.pool$pool == p),] ) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    
    ggsave(plot = map1, filename = paste("map_ensemble_annual_rich_phyto_baseline_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	
    ggsave(plot = map2, filename = paste("map_ensemble_annual_rich_zoo_baseline_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	
            
} # eo p in pools


# ---------------------------------------------------------

### 2°) Examine ensembles of changes in annual species richness
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("table_ann_changes_",dir())]; files
# Separate the ones for richness (no "beta.div") from those with beta.div
files2 <- files[grep("beta.div", files)]; files2
files1 <- files[!(files %in% files2)]; files1

require("parallel")
# f <- files1[1]
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
    }, mc.cores = 30
) # eo mclapply
# Rbind
table.diff <- dplyr::bind_rows(res)
dim(table.diff)
head(table.diff)
rm(res); gc()


### 26/06/2020: To describe the results based on ALL 80 ensemble members rather than from the spatial variability of the ensemble projection:
# - compute spatial averages for ALL 80 combinations in table, and THEN compute median (±IQR). Focus on those regions:
# - Southern Ocean (< -60°)
# - tropical band (-30°/+30°)
# - temperate latitudes (abs(ens$y) >= 40 & abs(ens$y))
# - arctic ocean: y > 70

### First, add an id for the ensemble member:
table.diff$member <- factor(paste(table.diff$ESM, table.diff$SDM, table.diff$pool, sep = "_"))
unique(table.diff$member) # 80 values as should be
# And a vector of region
table.diff$region <- NA
table.diff[which(table.diff$y <= -60),"region"] <- "SO"
table.diff[which(table.diff$y >= 70),"region"] <- "AO"
table.diff[which(abs(table.diff$y) <= 30),"region"] <- "Tropics"
table.diff[which(abs(table.diff$y) >= 40 & abs(table.diff$y) <= 55),"region"] <- "Temperate"
unique(table.diff$region) 

### 07/08/2020: For Niki, compute and map regions of model members agreement according to varying thresholds
ddf <- data.frame(table.diff %>% group_by(id) %>% summarize( 
                x = unique(x), y = unique(y), 
                n_posi_tot = sum(perc_rich_tot > 0)/80,
                n_nega_tot = sum(perc_rich_tot < 0)/80,
                n_posi_phyto = sum(perc_rich_phyto > 0)/80,
                n_nega_phyto = sum(perc_rich_phyto < 0)/80,
                n_posi_zoo = sum(perc_rich_zoo > 0)/80,
                n_nega_zoo = sum(perc_rich_zoo < 0)/80
        ) # eo summarize 
) # eo ddf
head(ddf) ; summary(ddf)

map_agreement_90_tot <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_tot >= 0.9,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_agreement_85_tot <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_tot >= 0.85,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_agreement_75_tot <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_tot >= 0.75,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
              
#
map_agreement_90_phyto <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_phyto >= 0.9,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_agreement_85_phyto <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_phyto >= 0.85,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_agreement_75_phyto <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_phyto >= 0.75,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#  
map_agreement_90_zoo <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_zoo >= 0.9,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_agreement_85_zoo <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_zoo >= 0.85,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_agreement_75_zoo <- ggplot() + geom_tile(aes(x = x, y = y), data = ddf[ddf$n_posi_zoo >= 0.75,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()

setwd(WD)              
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_plankton_90.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_90_tot)
dev.off()
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_plankton_85.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_85_tot)
dev.off()
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_plankton_75.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_75_tot)
dev.off()

grDevices::cairo_ps(filename = "map_model_member_agreement_diff_phytoplankton_90.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_90_phyto)
dev.off()
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_phytoplankton_85.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_85_phyto)
dev.off()
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_phytoplankton_75.eps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_75_phyto)
dev.off()

grDevices::cairo_ps(filename = "map_model_member_agreement_diff_zooplankton_90.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_90_zoo)
dev.off()
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_zooplankton_85.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_85_zoo)
dev.off()
grDevices::cairo_ps(filename = "map_model_member_agreement_diff_zooplankton_75.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_agreement_75_zoo)
dev.off()


### Use dplyr to compute median response per region and member and sign agreement
ddf <- data.frame(na.omit(table.diff) %>% group_by(member,region) %>% summarize( 
                perc_tot = median(perc_rich_tot,na.rm=T), 
                perc_phyto = median(perc_rich_phyto,na.rm=T),
                perc_zoo = median(perc_rich_zoo,na.rm=T)
        ) # eo summarize 
) # eo ddf
head(ddf) 

### Derive sign agreement:
ddf$sign_tot <- NA
ddf$sign_phyto <- NA
ddf$sign_zoo <- NA
for(i in c(1:nrow(ddf))) {
        
        message(paste(i, sep = ""))
        
        ### Use 3 if else loops to attribute sign of change
        if( ddf[i,"perc_tot"] > 0) {
                ddf[i,"sign_tot"] <- "increase"
        } else {
            ddf[i,"sign_tot"] <- "decrease"
        }
        
        # 
        if( ddf[i,"perc_phyto"] > 0) {
                ddf[i,"sign_phyto"] <- "increase"
        } else {
                ddf[i,"sign_phyto"] <- "decrease"
        }
        
        # 
        if( ddf[i,"perc_zoo"] > 0) {
                ddf[i,"sign_zoo"] <- "increase"
        } else {
               ddf[i,"sign_zoo"] <- "decrease"
        }
    
} # eo for loop

### Display % agreement
for(reg in unique(ddf$region)) {
    
        tot_posi <- ( (nrow(ddf[ddf$region == reg & ddf$sign_tot == "increase",]))/80 )*100
        tot_nega <- ( (nrow(ddf[ddf$region == reg & ddf$sign_tot == "decrease",]))/80 )*100
        message(paste("% agreement on total plankton increase in SR in the ",reg," is ",tot_posi, "%", sep = ""))
        message(paste("% agreement on total plankton decrease in SR in the ",reg," is ",tot_nega, "%", sep = ""))
        message(paste("", sep = ""))
        
        # Same with phyto and zoo
        phyto_posi <- ( (nrow(ddf[ddf$region == reg & ddf$sign_phyto == "increase",]))/80 )*100
        phyto_nega <- ( (nrow(ddf[ddf$region == reg & ddf$sign_phyto == "decrease",]))/80 )*100
        message(paste("% agreement in total phytoplankton increase in SR in the ",reg," is ",phyto_posi, "%", sep = ""))
        message(paste("% agreement on total phytoplankton decrease in SR in the ",reg," is ",phyto_nega, "%", sep = ""))
        message(paste("", sep = ""))
        
        zoo_posi <- ( (nrow(ddf[ddf$region == reg & ddf$sign_zoo == "increase",]))/80 )*100
        zoo_nega <- ( (nrow(ddf[ddf$region == reg & ddf$sign_zoo == "decrease",]))/80 )*100
        message(paste("% agreement in total zooplankton increase in SR in the ",reg," is ",zoo_posi, "%", sep = ""))
        message(paste("% agreement on total zooplankton decrease in SR in the ",reg," is ",zoo_nega, "%", sep = ""))
        message(paste("", sep = ""))
        
}

### And compute median ±IQR of regional averages for each resp var + agreement %
ddf2 <- data.frame(na.omit(ddf) %>% group_by(region) %>% summarize( 
                med_tot = median(perc_tot, na.rm=T), q25_t = quantile(perc_tot, na.rm=T)[2], q75_t = quantile(perc_tot, na.rm=T)[5],
                med_phyto = median(perc_phyto, na.rm=T), q25_p = quantile(perc_phyto, na.rm=T)[2], q75_p = quantile(perc_phyto, na.rm=T)[5],
                med_zoo = median(perc_zoo, na.rm=T), q25_z = quantile(perc_zoo, na.rm=T)[2], q75_z = quantile(perc_zoo, na.rm=T)[5]
        ) # eo summarize 
) # eo ddf
ddf2


### And the global numbers
ddf_glob <- data.frame(na.omit(table.diff) %>% group_by(member) %>% summarize( 
                perc_tot = median(perc_rich_tot,na.rm=T), 
                perc_phyto = median(perc_rich_phyto,na.rm=T),
                perc_zoo = median(perc_rich_zoo,na.rm=T)
        ) # eo summarize 
) # eo ddf
head(ddf_glob) 
# Display global reponses
median(ddf_glob$perc_tot) ; IQR(ddf_glob$perc_tot) ; quantile(ddf_glob$perc_tot)
median(ddf_glob$perc_phyto) ; IQR(ddf_glob$perc_phyto) ; quantile(ddf_glob$perc_phyto)
median(ddf_glob$perc_zoo) ; IQR(ddf_glob$perc_zoo) ; quantile(ddf_glob$perc_zoo)
(nrow(ddf_glob[which(ddf_glob$perc_tot > 0),])/80)*100
(nrow(ddf_glob[which(ddf_glob$perc_phyto > 0),])/80)*100
(nrow(ddf_glob[which(ddf_glob$perc_zoo > 0),])/80)*100


### Compute global zonal averages across all 80 ensemble members for zoo and phyto responses and plot panel
glob_zonal <- data.frame(table.diff %>% group_by(member,y) %>% summarize( 
                perc_tot = mean(perc_rich_tot,na.rm=T), 
                perc_phyto = mean(perc_rich_phyto,na.rm=T),
                perc_zoo = mean(perc_rich_zoo,na.rm=T)
        ) # eo summarize 
) # eo ddf
head(glob_zonal) 
# Change labels 
glob_zonal$member2 <- rep(c(1:80), each = 179, len = 14320)

summary(glob_zonal)

panel <- ggplot() + geom_path(aes(y = y, x = perc_phyto), data = na.omit(glob_zonal), colour = "black", linetype = "solid") +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), position = "right") +
   scale_x_continuous(name = "", limits = c(-40,200), breaks = c(-40,-20,0,50,100,200)) + 
   geom_vline(xintercept = 0, linetype = "dotted") +
   theme_classic() + facet_wrap(~factor(member2), nrow = 40, ncol = 20, scales = "fixed")
#
setwd(WD)
ggsave(plot = panel, filename = "ensemble_members_∆SRphyto_zonal.jpg", dpi = 300, width = 30, height = 10)

panel <- ggplot() + geom_path(aes(y = y, x = perc_zoo), data = na.omit(glob_zonal), colour = "black", linetype = "solid") +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), position = "right") +
   scale_x_continuous(name = "", limits = c(-45,170), breaks = c(-45,-20,0,25,50,100,170)) + 
   geom_vline(xintercept = 0, linetype = "dotted") +
   theme_classic() + facet_wrap(~factor(member2), nrow = 40, ncol = 20, scales = "fixed")
#
ggsave(plot = panel, filename = "ensemble_members_∆SRzoo_zonal.jpg", dpi = 300, width = 30, height = 10)

### Examine the members that highly diverge for phyto (13:16, 29:32)
glob_zonal[which(glob_zonal$member2 %in% c(13:16)),"member"]
glob_zonal[which(glob_zonal$member2 %in% c(29:32)),"member"]
glob_zonal[which(glob_zonal$member2 %in% c(77:80)),"member"]


### 26/06/20: Perform non parametric tests of variance to test how signif the changes in SR are, on global and regional scales
### First, need to combine table.base with table.diff --> already exists in "/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections"
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# Example of data.frame to use
files <- dir()[grep("fut_SR_",dir())]
f <- files[1]
res <- mclapply(files, function(f) {
            
            # Useless message
            message(paste("Loading results from ",f, sep = ""))
            t <- get(load(f))
            # Get attributes from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            sdm <- terms[5,1] 
            esm <- terms[4,1]
            p <- terms[6,1]
            pool <- str_replace(p,".Rdata","")
            t$SDM <- sdm
            t$ESM <- esm
            t$pool <- pool
            t <- t[order(t$id),]
            
            # Return
            return(t)
            
    }, mc.cores = 30
    
) # eo mclapply
# Rbind
table <- dplyr::bind_rows(res)
rm(res); gc()
colnames(table)[6] <- "rich_tot"
dim(table) ; head(table)

### First, add an ensemble member IR
table$member <- paste(table$SDM, table$ESM, table$pool, sep = "_")
### Add region like above
table$region <- NA
table[which(table$y <= -60),"region"] <- "SO"
table[which(table$y >= 70),"region"] <- "AO"
table[which(abs(table$y) <= 30),"region"] <- "Tropics"
table[which(abs(table$y) >= 40 & abs(table$y) <= 55),"region"] <- "Temperate"

### Compute medians of base and future SR per region and member
ens.med <- data.frame(table %>% group_by(region,member) %>% summarize(
            base_tot = median(rich_tot,na.rm=T), 
            base_phyto = median(rich_phyto,na.rm=T),
            base_zoo = median(rich_zoo,na.rm=T),
            fut_tot = median(fut_rich_tot,na.rm=T), 
            fut_phyto = median(fut_rich_phyto,na.rm=T),
            fut_zoo = median(fut_rich_zoo,na.rm=T)
    ) # eo summarize
) # eo ddf
summary(ens.med)
head(ens.med)

### Perform wilcox test
wilcox.test(na.omit(ens.med[ens.med$region == "AO","base_zoo"]), na.omit(ens.med[ens.med$region == "AO","fut_zoo"]), paired = T)
wilcox.test(na.omit(ens.med[ens.med$region == "SO","base_zoo"]), na.omit(ens.med[ens.med$region == "SO","fut_zoo"]), paired = T)
wilcox.test(na.omit(ens.med[ens.med$region == "Temperate","base_zoo"]), na.omit(ens.med[ens.med$region == "Temperate","fut_zoo"]), paired = T)
wilcox.test(na.omit(ens.med[ens.med$region == "Tropics","base_zoo"]), na.omit(ens.med[ens.med$region == "Tropics","fut_zoo"]), paired = T)

### Same as above but global
glob.med <- data.frame(table %>% group_by(member) %>% summarize(
            base_tot = median(rich_tot,na.rm=T), 
            base_phyto = median(rich_phyto,na.rm=T),
            base_zoo = median(rich_zoo,na.rm=T),
            fut_tot = median(fut_rich_tot,na.rm=T), 
            fut_phyto = median(fut_rich_phyto,na.rm=T),
            fut_zoo = median(fut_rich_zoo,na.rm=T)
    ) # eo summarize
) # eo ddf
head(glob.med)

### Perform wilcox test
wilcox.test(na.omit(glob.med[,"base_tot"]), na.omit(glob.med[,"fut_tot"]), paired = T)
wilcox.test(na.omit(glob.med[,"base_phyto"]), na.omit(glob.med[,"fut_phyto"]), paired = T)
wilcox.test(na.omit(glob.med[,"base_zoo"]), na.omit(glob.med[,"fut_zoo"]), paired = T)


### 23/04/2020: For describing results in the manuscript
# summary(ens)
# IQR(ens$rich_phyto, na.rm = T)
# # colnames(ens)
# summary(ens[abs(ens$y) < 30,])
# IQR(ens[abs(ens$y) < 30,"rich_zoo"], na.rm = T)
#
# summary(ens[abs(ens$y) >= 40 & abs(ens$y) <= 55,])
# IQR(ens[abs(ens$y) >= 40 & abs(ens$y) <= 55,"rich_tot"], na.rm = T)
# IQR(ens[abs(ens$y) >= 40 & abs(ens$y) <= 55,"rich_phyto"], na.rm = T)
# IQR(ens[abs(ens$y) >= 40 & abs(ens$y) <= 55,"rich_zoo"], na.rm = T)
#
# # For the Southern Ocean
# summary(ens[ens$y < -60,])
# IQR(ens[ens$y < -60,"rich_zoo"], na.rm = T)
#
# # For the phyto in the arctic
# summary(ens[ens$y > 75,])
# IQR(ens[ens$y > 75,"rich_phyto"], na.rm = T)


### 01/04/2020: Compute ratio of average to sd for signifiance test?
ens$robust_tot <- (ens$rich_tot)/(ens$sd_rich_tot)
ens$robust_phyto <- (ens$rich_phyto)/(ens$sd_rich_phyto)
ens$robust_zoo <- (ens$rich_zoo)/(ens$sd_rich_zoo)

ggplot() + geom_point(aes(x = rich_tot, y = sd_rich_tot), data = ens, colour = "grey50", alpha = 0.5) + 
    geom_abline(intercept = 1, linetype = "longdash") + geom_vline(xintercept = 0, linetype = "dashed") + 
    xlab("Mean plankton richness difference (%)") + ylab("Standard deviation") + theme_classic()
#    
ggplot() + geom_point(aes(x = rich_tot, y = robust_tot), data = ens, colour = "grey50", alpha = 0.5) + 
    geom_hline(yintercept = 1, linetype = "dashed") + xlab("Mean phytoplankton richness difference (%)") + 
    ylab("Ensemble projection robustness\n(ensemble mean:inter-model std)") + theme_classic()
    
ggplot() + geom_point(aes(x = rich_phyto, y = robust_phyto), data = ens, colour = "grey50", alpha = 0.5) + 
    geom_hline(yintercept = 1, linetype = "dashed") + xlab("Mean phytoplankton richness difference (%)") + 
    ylab("Ensemble projection robustness\n(ensemble mean:inter-model std)") + theme_classic()
#
ggplot() + geom_point(aes(x = rich_zoo, y = robust_zoo), data = ens, colour = "grey50", alpha = 0.5) + 
    geom_hline(yintercept = 1, linetype = "dashed") + xlab("Mean zooplankton richness difference (%)") + 
    ylab("Ensemble projection robustness\n(ensemble mean:inter-model std)") + theme_classic()


### Make maps ! 
### Compute global zonal averages across all 80 ensemble members for zoo and phyto responses and plot panel
ens <- data.frame(table.diff %>% group_by(id) %>% summarize( 
                x = unique(x), y = unique(y), 
                rich_tot = mean(perc_rich_tot,na.rm=T), 
                rich_phyto = mean(perc_rich_phyto,na.rm=T),
                rich_zoo = mean(perc_rich_zoo,na.rm=T)
        ) # eo summarize 
) # eo ddf
summary(ens) 

min <- floor(min(ens[,c("rich_tot","rich_phyto","rich_zoo")], na.rm = T) ) ; min

map1 <- ggplot() + geom_tile(aes(x = x, y = y, fill = rich_tot), data = ens[ens$rich_tot < 50,]) +
    geom_tile(aes(x = x, y = y), data = ens[ens$rich_tot >= 50,], fill = "#b2182b") +
    #geom_point(aes(x = x, y = y), data = ens[ens$robust_tot > 1,], colour = "grey30", alpha = 0.1, size = 0.05, shape = 4) + 
    #geom_contour(colour = "grey25", breaks = 1, size = 0.25, aes(x = x, y = y, z = robust_tot), data = ens) +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_tot), data = ens[ens$rich_tot < 50,] ) +
 	scale_fill_gradient2(name = "Richness\ndifference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

map2 <- ggplot() + geom_tile(aes(x = x, y = y, fill = rich_phyto), data = ens[ens$rich_phyto < 50,]) +
    geom_tile(aes(x = x, y = y), data = ens[ens$rich_phyto >= 50,], fill = "#b2182b") +
    #geom_point(aes(x = x, y = y), data = ens[ens$robust_phyto > 1,], colour = "grey30", alpha = 0.1, size = 0.05, shape = 4) + 
    #geom_contour(colour = "grey25", breaks = 1, size = 0.25, aes(x = x, y = y, z = robust_phyto), data = ens ) +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), data = ens[ens$rich_phyto < 50,] ) +
 	scale_fill_gradient2(name = "Richness\ndifference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

map3 <- ggplot() + geom_tile(aes(x = x, y = y, fill = rich_zoo), data = ens[ens$rich_zoo < 50,]) +
    geom_tile(aes(x = x, y = y), data = ens[ens$rich_zoo >= 50,], fill = "#b2182b") +
    #geom_point(aes(x = x, y = y), data = ens[ens$robust_zoo > 1,], colour = "grey30", alpha = 0.1, size = 0.05, shape = 4) + 
    #geom_contour(colour = "grey25", breaks = 1, size = 0.25, aes(x = x, y = y, z = robust_zoo), data = ens) +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), data = ens[ens$rich_zoo < 50,] ) +
 	scale_fill_gradient2(name = "Richness\ndifference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
              labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")
   
setwd(WD)       
ggsave(plot = map1, filename = "map_ensemble_annual_perc_rich_tot_robust.ps", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_ensemble_annual_perc_rich_phyto_robust.ps", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_ensemble_annual_perc_rich_zoo_robust.ps", dpi = 300, width = 7, height = 5)

grDevices::cairo_ps(filename = "map_ensemble_annual_perc_rich_tot_robust.ps", width = 7, height = 5, fallback_resolution = 300)
print(map1)
dev.off()
grDevices::cairo_ps(filename = "map_ensemble_annual_perc_rich_phyto_robust.ps", width = 7, height = 5, fallback_resolution = 300)
print(map2)
dev.off()
grDevices::cairo_ps(filename = "map_ensemble_annual_perc_rich_zoo_robust.ps", width = 7, height = 5, fallback_resolution = 300)
print(map3)
dev.off()


### Also standard deviation to highlight uncertainties
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd_rich_tot), data = ens[!is.na(ens$sd_rich_tot),]) +
    geom_contour(colour = "#fff7bc", binwidth = 25, size = 0.25, aes(x = x, y = y, z = sd_rich_tot), data = ens[!is.na(ens$sd_rich_tot),]) +
 	scale_fill_viridis(name = "Uncertainty\n(stdev)", option = "cividis") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd_rich_phyto), data = ens[!is.na(ens$sd_rich_phyto),]) +
    geom_contour(colour = "#fff7bc", binwidth = 25, size = 0.25, aes(x = x, y = y, z = sd_rich_phyto), data = ens[!is.na(ens$sd_rich_phyto),]) +
 	scale_fill_viridis(name = "Uncertainty\n(stdev)", option = "cividis") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd_rich_zoo), data = ens[!is.na(ens$sd_rich_zoo),]) +
    geom_contour(colour = "#fff7bc", binwidth = 25, size = 0.25, aes(x = x, y = y, z = sd_rich_zoo), data = ens[!is.na(ens$sd_rich_zoo),]) +
 	scale_fill_viridis(name = "Uncertainty\n(stdev)", option = "cividis") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
setwd(WD)       
ggsave(plot = map1, filename = paste("map_ensemble_annual_stdev_rich_tot.jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = paste("map_ensemble_annual_stdev_rich_phyto.jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = paste("map_ensemble_annual_stdev_rich_zoo.jpg", sep = ""), dpi = 300, width = 7, height = 5)
        

### 23/06/2020: And plot zonal patterns by separating ESMs
head(table)

zonal <- data.frame(table.diff %>% group_by(y,ESM) %>%
             summarise(perc_tot = mean(perc_rich_tot, na.rm = T),
             perc_phyto = mean(perc_rich_phyto, na.rm = T),
             perc_zoo = mean(perc_rich_zoo, na.rm = T) )
) # eo ddf
summary(zonal)
unique(zonal$ESM)

plot2 <- ggplot() + 
   geom_path(aes(y = y, x = perc_phyto), data = na.omit(zonal[zonal$ESM == "MRI-NEMURO",]), colour = "black", linetype = "solid") +
   geom_path(aes(y = y, x = perc_phyto), data = na.omit(zonal[zonal$ESM == "IPSL-PISCES",]), colour = "black", linetype = "longdash") +
   geom_path(aes(y = y, x = perc_phyto), data = na.omit(zonal[zonal$ESM == "GFDL-TOPAZ",]), colour = "black", linetype = "dashed") +
   geom_path(aes(y = y, x = perc_phyto), data = na.omit(zonal[zonal$ESM == "CNRM-PISCES",]), colour = "grey50", linetype = "solid") +
   geom_path(aes(y = y, x = perc_phyto), data = na.omit(zonal[zonal$ESM == "CESM-BEC",]), colour = "grey50", linetype = "dashed") +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), position = "right") +
   scale_x_continuous(name = "", limits = c(-25,95), breaks = c(-25,0,25,50,75,95)) + 
   geom_vline(xintercept = 0, linetype = "dotted") +
   theme_classic() 

plot3 <- ggplot() + 
   geom_path(aes(y = y, x = perc_zoo), data = na.omit(zonal[zonal$ESM == "MRI-NEMURO",]), colour = "black", linetype = "solid") +
   geom_path(aes(y = y, x = perc_zoo), data = na.omit(zonal[zonal$ESM == "IPSL-PISCES",]), colour = "black", linetype = "longdash") +
   geom_path(aes(y = y, x = perc_zoo), data = na.omit(zonal[zonal$ESM == "GFDL-TOPAZ",]), colour = "black", linetype = "dashed") +
   geom_path(aes(y = y, x = perc_zoo), data = na.omit(zonal[zonal$ESM == "CNRM-PISCES",]), colour = "grey50", linetype = "solid") +
   geom_path(aes(y = y, x = perc_zoo), data = na.omit(zonal[zonal$ESM == "CESM-BEC",]), colour = "grey50", linetype = "dashed") +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), position = "right") +
   scale_x_continuous(name = "", limits = c(-25,95), breaks = c(-25,0,25,50,75,95)) + 
   geom_vline(xintercept = 0, linetype = "dotted") +
   theme_classic() 

plot1 <- ggplot() + 
   geom_path(aes(y = y, x = perc_tot), data = na.omit(zonal[zonal$ESM == "MRI-NEMURO",]), colour = "black", linetype = "solid") +
   geom_path(aes(y = y, x = perc_tot), data = na.omit(zonal[zonal$ESM == "IPSL-PISCES",]), colour = "black", linetype = "longdash") +
   geom_path(aes(y = y, x = perc_tot), data = na.omit(zonal[zonal$ESM == "GFDL-TOPAZ",]), colour = "black", linetype = "dashed") +
   geom_path(aes(y = y, x = perc_tot), data = na.omit(zonal[zonal$ESM == "CNRM-PISCES",]), colour = "grey50", linetype = "solid") +
   geom_path(aes(y = y, x = perc_tot), data = na.omit(zonal[zonal$ESM == "CESM-BEC",]), colour = "grey50", linetype = "dashed") +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), position = "right") +
   scale_x_continuous(name = "", limits = c(-25,95), breaks = c(-25,0,25,50,75,95)) + 
   geom_vline(xintercept = 0, linetype = "dotted") +
   theme_classic() 

grDevices::cairo_ps(filename = "plot_zonal_ann_plankton_perc_ensemble.ps", width = 3, height = 4, fallback_resolution = 300)
print(plot1)
dev.off()
grDevices::cairo_ps(filename = "plot_zonal_ann_phyto_perc_ensemble.ps", width = 3, height = 4, fallback_resolution = 300)
print(plot2)
dev.off()
grDevices::cairo_ps(filename = "plot_zonal_ann_zoo_perc_ensemble.ps", width = 3, height = 4, fallback_resolution = 300)
print(plot3)
dev.off()


# summary(ens)
# zonal <- data.frame(ens %>%
#         group_by(y) %>%
#         summarise(mean_plankton = mean(rich_tot, na.rm=T), sd_plankton = sd(rich_tot, na.rm=T),
#         mean_phyto = mean(rich_phyto, na.rm=T), sd_phyto = sd(rich_phyto, na.rm=T),
#         mean_zoo = mean(rich_zoo, na.rm=T), sd_zoo = sd(rich_zoo, na.rm=T) )
# ) # eo ddf
# summary(zonal)
#
# # Define maxi and mins
# min <- -20
# max <- 70
#
# plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_plankton - sd_plankton, xmax = mean_plankton + sd_plankton), fill = "grey70", data = na.omit(zonal)) +
#    geom_path(aes(y = y, x = mean_plankton), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(min,max)) +
#    xlab("Mean richness difference (%)") + scale_y_continuous(position = "right", name = "Latitude (°)") +
#    geom_vline(xintercept = 0, linetype = "dashed") + theme_classic()
#
# plot2 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_phyto - sd_phyto, xmax = mean_phyto + sd_phyto), fill = "grey70", data = na.omit(zonal)) +
#    geom_path(aes(y = y, x = mean_phyto), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(min,max)) +
#    xlab("Mean richness difference (%)") + scale_y_continuous(position = "right", name = "Latitude (°)") +
#    geom_vline(xintercept = 0, linetype = "dashed") + theme_classic()
# #
# plot3 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_zoo - sd_zoo, xmax = mean_zoo + sd_zoo), fill = "grey70", data = na.omit(zonal)) +
#    geom_path(aes(y = y, x = mean_zoo), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(min,max)) +
#    xlab("Mean richness difference (%)") + scale_y_continuous(position = "right", name = "Latitude (°)") +
#    geom_vline(xintercept = 0, linetype = "dashed") + theme_classic()
#
# #
# ggsave(plot = plot1, filename = "plot_zonal_ann_plankton_perc_ensemble.jpg", dpi = 300, width = 3, height = 4)
# ggsave(plot = plot2, filename = "plot_zonal_ann_phyto_perc_ensemble.jpg", dpi = 300, width = 3, height = 4)
# ggsave(plot = plot3, filename = "plot_zonal_ann_zoo_perc_ensemble.jpg", dpi = 300, width = 3, height = 4)



### 23/06/20: Use ggarrange and save panels
library("ggpubr")
setwd("/net/kryo/work/fabioben/OVERSEE/data")
ggarrange(map1, plot1, map2, plot2, map3, plot3, labels = c("G","H","I","J","K","L"), ncol = 2, nrow = 3, widths = c(2.7,1) )
panel <- ggarrange(map1, plot1, map2, plot2, map3, plot3, labels = c("G","H","I","J","K","L"), ncol = 2, nrow = 3, widths = c(2.7,1) )
ggsave(plot = panel, filename = "panel_ensemble_future_perc_Fig1G-L.eps", dpi = 300, width = 8, height = 8)

grDevices::cairo_ps(filename = "panel_ensemble_future_perc_Fig1G-L.ps", width = 8, height = 8, fallback_resolution = 300)
print(panel)
dev.off()

### And save a map for colorbar
map4colorbar <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_tot), data = ens[ens$rich_tot < 50,]) +
    geom_raster(aes(x = x, y = y), data = ens[ens$rich_tot >= 50,], fill = "#b2182b") +
    #geom_point(aes(x = x, y = y), data = ens[ens$robust_tot > 1,], colour = "grey30", alpha = 0.1, size = 0.05, shape = 4) + 
    #geom_contour(colour = "grey25", breaks = 1, size = 0.25, aes(x = x, y = y, z = robust_tot), data = ens) +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_tot), data = ens[ens$rich_tot < 50,] ) +
 	scale_fill_gradient2(name = "Richness\ndifference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")
  
# ggsave(plot = map4colorbar, filename = paste("map4colorbar2.pdf", sep = ""), dpi = 300, width = 8, height = 5)

grDevices::cairo_ps(filename = "map4colorbar2.eps", width = 8, height = 5, fallback_resolution = 300)
print(map4colorbar)
dev.off()

### Per SDM (for facet)
ens.SDM <- data.frame(table %>%
        group_by(id,SDM) %>%
        summarize(x = unique(x), y = unique(y), 
            rich_tot = mean(perc_rich_tot, na.rm = T),
            rich_phyto = mean(perc_rich_phyto, na.rm = T),
            rich_zoo = mean(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
summary(ens.SDM) # dim(ens.SDM)

min <- floor(min(ens.SDM[,c("rich_tot","rich_phyto","rich_zoo")], na.rm = T) )
min

for(sdm in unique(ens.SDM$SDM) ) {
    
    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_tot), 
                data = ens.SDM[which(ens.SDM$rich_tot < 50 & ens.SDM$SDM == sdm),]) +
        geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_tot >= 50 & ens.SDM$SDM == sdm),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_tot), 
                data = ens.SDM[which(ens.SDM$rich_tot < 50 & ens.SDM$SDM == sdm),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
            data = ens.SDM[which(ens.SDM$rich_phyto < 50 & ens.SDM$SDM == sdm),]) +
        geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_phyto >= 50 & ens.SDM$SDM == sdm),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
                data = ens.SDM[which(ens.SDM$rich_phyto < 50 & ens.SDM$SDM == sdm),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
                data = ens.SDM[which(ens.SDM$rich_zoo < 50 & ens.SDM$SDM == sdm),]) +
        geom_raster(aes(x = x, y = y), data = ens.SDM[which(ens.SDM$rich_zoo >= 50 & ens.SDM$SDM == sdm),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
                data = ens.SDM[which(ens.SDM$rich_zoo < 50 & ens.SDM$SDM == sdm),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
    ggsave(plot = map1, filename = paste("map_ensemble_annual_perc_rich_tot_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    ggsave(plot = map2, filename = paste("map_ensemble_annual_perc_rich_phyto_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    ggsave(plot = map3, filename = paste("map_ensemble_annual_perc_rich_zoo_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                
} # eo sdm in SDM


# Per ESM (for facet)
ens.ESM <- data.frame(table %>%
        group_by(id,ESM) %>%
        summarize(x = unique(x), y = unique(y), 
        rich_tot = mean(perc_rich_tot, na.rm = T),
        rich_phyto = mean(perc_rich_phyto, na.rm = T),
        rich_zoo = mean(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
head(ens.ESM)

min <- floor(min(ens.ESM[,c("rich_tot","rich_phyto","rich_zoo")], na.rm = T) )
min

for(esm in unique(ens.ESM$ESM) ) {
    
    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_tot), 
                data = ens.ESM[which(ens.ESM$rich_tot < 50 & ens.ESM$ESM == esm),]) +
        geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_tot >= 50 & ens.ESM$ESM == esm),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_tot), 
                data = ens.ESM[which(ens.ESM$rich_tot < 50 & ens.ESM$ESM == esm),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
            data = ens.ESM[which(ens.ESM$rich_phyto < 50 & ens.ESM$ESM == esm),]) +
        geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_phyto >= 50 & ens.ESM$ESM == esm),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
                data = ens.ESM[which(ens.ESM$rich_phyto < 50 & ens.ESM$ESM == esm),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
                data = ens.ESM[which(ens.ESM$rich_zoo < 50 & ens.ESM$ESM == esm),]) +
        geom_raster(aes(x = x, y = y), data = ens.ESM[which(ens.ESM$rich_zoo >= 50 & ens.ESM$ESM == esm),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
                data = ens.ESM[which(ens.ESM$rich_zoo < 50 & ens.ESM$ESM == esm),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
    ggsave(plot = map1, filename = paste("map_ensemble_annual_perc_rich_tot_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    ggsave(plot = map2, filename = paste("map_ensemble_annual_perc_rich_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    ggsave(plot = map3, filename = paste("map_ensemble_annual_perc_rich_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                
} # eo sdm in SDM

# Per pred pool (for facet)
ens.pool <- data.frame(table %>%
        group_by(id,pool) %>%
        summarize(x = unique(x), y = unique(y), 
        rich_tot = mean(perc_rich_tot, na.rm = T),
        rich_phyto = mean(perc_rich_phyto, na.rm = T),
        rich_zoo = mean(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
summary(ens.pool)

min <- floor(min(ens.pool[,c("rich_tot","rich_phyto","rich_zoo")], na.rm = T) )
min

for(p in unique(ens.pool$pool) ) {
    
    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_tot), 
                data = ens.pool[which(ens.pool$rich_tot < 50 & ens.pool$pool == p),]) +
        geom_raster(aes(x = x, y = y), data = ens.pool[which(ens.pool$rich_tot >= 50 & ens.pool$pool == p),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_tot), 
                data = ens.pool[which(ens.pool$rich_tot < 50 & ens.pool$pool == p),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto),
            data = ens.pool[which(ens.pool$rich_phyto < 50 & ens.pool$pool == p),]) +
        geom_raster(aes(x = x, y = y), data = ens.pool[which(ens.pool$rich_phyto >= 50 & ens.pool$pool == p),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_phyto), 
                data = ens.pool[which(ens.pool$rich_phyto < 50 & ens.pool$pool == p),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

    map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo),
                data = ens.pool[which(ens.pool$rich_zoo < 50 & ens.pool$pool == p),]) +
        geom_raster(aes(x = x, y = y), data = ens.pool[which(ens.pool$rich_zoo >= 50 & ens.pool$pool == p),], fill = "#b2182b") +
        geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = rich_zoo), 
                data = ens.pool[which(ens.pool$rich_zoo < 50 & ens.pool$pool == p),] ) +
     	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
    ggsave(plot = map1, filename = paste("map_ensemble_annual_perc_rich_tot_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    ggsave(plot = map2, filename = paste("map_ensemble_annual_perc_rich_phyto_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    ggsave(plot = map3, filename = paste("map_ensemble_annual_perc_rich_zoo_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                
} # eo sdm in SDM

### Clearly the RF models can't really be compared to the other...which is problematic for the turn-over metrics because one cannot rely on the same thresholds...

# ----------------------------------------------------------

### 3°) Re-define the thresholds to be used for each SDM type ! choose a range of 10 thresholds that bets match rich patterns based on HSI
setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
files.phyto <- dir()[grep("table_ann_compo_phyto_baseline",dir())]; files.phyto
files.zoo <- dir()[grep("table_ann_compo_zoo_baseline",dir())]; files.zoo
### For each file, compute anual richness, and examine for each SDM which thresholds work best to match the richness pattern
# f <- files.phyto[1]
res.phyto <- mclapply(files.phyto, function(f) {
                
                # Read data
                message(paste("Loading ",f, sep = ""))
                table <- get(load(f))
                terms <- do.call(cbind, strsplit(as.character(f),"_"))
                sdm <- terms[6,1] 
                pool <- terms[7,1] 
                pool <- str_replace(as.character(pool),".Rdata","")
                # head(table)
                table$rich <- rowSums(as.matrix(table[,c(4:length(table))]))
                
                ### For various thresholds, binarily transform the HSI and re-compute richness, assess r2
                message(paste("Thresholding", sep = ""))
                message(paste("",sep = ""))
                thresholds <- seq(from = 0.1, to = 0.25, by = 0.01)
                # t <- 0.5
                thresh <- lapply(thresholds, function(t) { 
                
                            message(paste("Using t = ",t, sep = ""))
                            table2 <- table
                            # Convert to P/1 using rbinom and initial probability
                            for(sp in colnames(table2)[c(4:c(length(table2)-1))] ) {
                                table2[,c(sp)][table2[,c(sp)] > t] <- 1
                                table2[,c(sp)][table2[,c(sp)] <= t] <- 0
                            } # eo for loop
                            # Compute rich again from 1/0
                            table2$rich2 <- rowSums(as.matrix(table2[c(4:c(length(table2)-1))] ))
                            # summary(table2[,c("rich","rich2")])
                            lm <- lm(rich2 ~ rich, data = table2)
                            # str(summary(lm))
                            r2 <- summary(lm)$adj.r.squared
                            rho <- cor(table2$rich2, table2$rich, method = "spearman")
                            # Return
                            ddf <- data.frame(t = t, r2 = r2, cor = rho)
                
                    } # eo FUN 
                ) # eo lapply
                # Rbind
                ddf <- bind_rows(thresh)
                ddf$SDM <- sdm
                ddf$pool <- pool
                # summary(ddf)
                # Rturn 
                return(ddf)
            }, mc.cores = 12
            
) # eo mclapply
# Rbind
ddf.phyto <- bind_rows(res.phyto)
head(ddf.phyto); dim(ddf.phyto)
summary(ddf.phyto)

### Same for zoo
res.zoo <- mclapply(files.zoo, function(f) {
                
                # Read data
                message(paste("Loading ",f, sep = ""))
                table <- get(load(f))
                terms <- do.call(cbind, strsplit(as.character(f),"_"))
                sdm <- terms[6,1] 
                pool <- terms[7,1] 
                pool <- str_replace(as.character(pool),".Rdata","")
                # head(table)
                table$rich <- rowSums(as.matrix(table[,c(4:length(table))]))
                
                ### For various thresholds, binarily transform the HSI and re-compute richness, assess r2
                message(paste("Thresholding", sep = ""))
                message(paste("",sep = ""))
                thresholds <- seq(from = 0.1, to = 0.25, by = 0.01)
                # t <- 0.5
                thresh <- lapply(thresholds, function(t) { 
                
                            message(paste("Using t = ",t, sep = ""))
                            table2 <- table
                            # Convert to P/1 using rbinom and initial probability
                            for(sp in colnames(table2)[c(4:c(length(table2)-1))] ) {
                                table2[,c(sp)][table2[,c(sp)] > t] <- 1
                                table2[,c(sp)][table2[,c(sp)] <= t] <- 0
                            } # eo for loop
                            # Compute rich again from 1/0
                            table2$rich2 <- rowSums(as.matrix(table2[c(4:c(length(table2)-1))] ))
                            # summary(table2[,c("rich","rich2")])
                            lm <- lm(rich2 ~ rich, data = table2)
                            # str(summary(lm))
                            r2 <- summary(lm)$adj.r.squared
                            rho <- cor(table2$rich2, table2$rich, method = "spearman")
                            # Return
                            ddf <- data.frame(t = t, r2 = r2, cor = rho)
                
                    } # eo FUN 
                ) # eo lapply
                
                # Rbind
                ddf <- bind_rows(thresh)
                ddf$SDM <- sdm
                ddf$pool <- pool
                # summary(ddf)
                # Rturn 
                return(ddf)
            }, mc.cores = 12
            
) # eo mclapply
# Rbind
ddf.zoo <- bind_rows(res.zoo)
head(ddf.zoo); dim(ddf.zoo)
summary(ddf.zoo)

# Make room
rm(res.phyto,res.zoo); gc()

# Examine distrbution of r2 and correlation coeff across thresholds and SDMs (facet_grid)
library("RColorBrewer")
library("wesanderson")
library("ggthemes")

plot1 <- ggplot(data = ddf.phyto, aes(x = factor(t), y = r2, fill = factor(SDM))) + geom_violin(colour = "black") +
    scale_fill_manual(name = "ESM", values = wes_palettes$Zissou1) + 
    xlab("") + ylab("Adjusted R2") + scale_y_continuous(limits = c(0,1)) + theme_classic() +
    facet_wrap(~ factor(ddf.phyto$SDM), ncol = 2)

plot2 <- ggplot(data = ddf.zoo, aes(x = factor(t), y = r2, fill = factor(SDM))) + geom_violin(colour = "black") +
    scale_fill_manual(name = "ESM", values = wes_palettes$Zissou1) + 
    xlab("") + ylab("Adjusted R2") + scale_y_continuous(limits = c(0,1)) + theme_classic() +
    facet_wrap(~ factor(ddf.zoo$SDM), ncol = 2)
#
plot3 <- ggplot(data = ddf.phyto, aes(x = factor(t), y = cor, fill = factor(SDM))) + geom_violin(colour = "black") +
    scale_fill_manual(name = "ESM", values = wes_palettes$Zissou1) + 
    xlab("") + ylab("Spearman's correlation coefficient") + theme_classic() +
    facet_wrap(~ factor(ddf.phyto$SDM), ncol = 2)

plot4 <- ggplot(data = ddf.zoo, aes(x = factor(t), y = cor, fill = factor(SDM))) + geom_violin(colour = "black") +
    scale_fill_manual(name = "ESM", values = wes_palettes$Zissou1) + 
    xlab("") + ylab("Spearman's correlation coefficient") + theme_classic() +
    facet_wrap(~ factor(ddf.zoo$SDM), ncol = 2)

# Save
ggsave(plot = plot1, filename = "plot_thresholds_r2_phytov3.jpg", dpi = 300, width = 10, height = 10)
ggsave(plot = plot2, filename = "plot_thresholds_r2_zoov3.jpg", dpi = 300, width = 10, height = 10)
ggsave(plot = plot3, filename = "plot_thresholds_cor_phytov3.jpg", dpi = 300, width = 10, height = 10)
ggsave(plot = plot4, filename = "plot_thresholds_cor_zoov3.jpg", dpi = 300, width = 10, height = 10)

### Ok, So thresholds to be used for computing beta.div changes
# GLM/GAM/ANN: 0.25-0.4, by 0.01 ; seq(from = 0.25, to = 0.39, by = 0.01) 
# RF: seq(from = 0.10, to = 0.24, by = 0.01) 


# ----------------------------------------------------------

### 4°) Examine ensembles of changes in annual species composition
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("table_ann_changes_",dir())]; files
# Separate the ones for richness (no "beta.div") from those with beta.div
files2 <- files[grep("beta.div", files)]; files2
files1 <- files[!(files %in% files2)]; files1

require("parallel")
# f <- files2[1]
res <- mclapply(files2, function(f) {
            # Useless message
            message(paste("Loading results from ",f, sep = ""))
            t <- get(load(f))
            # Get attributes from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            t$ESM <- terms[6,1] 
            t$SDM <- terms[7,1] 
            p <- terms[8,1]
            t$pool <- str_replace(as.character(p),".Rdata","")
            # Return
            return(t)
    }, mc.cores = 30
) # eo mclapply
# Rbind
table <- dplyr::bind_rows(res)
dim(table)
head(table)
rm(res); gc()
# unique(table$jac)
colnames(table)
summary(table)

### Computing ensembles (all, SDM/ESM/pool)
ens <- data.frame(table %>%
        group_by(id) %>%
        summarize(x = unique(x), y = unique(y), 
            jac = mean(jac, na.rm = T),
            jne = mean(jne, na.rm = T),
            jtu = mean(jtu, na.rm = T),
            jac_phyto = mean(jac_phyto, na.rm = T),
            jne_phyto = mean(jne_phyto, na.rm = T),
            jtu_phyto = mean(jtu_phyto, na.rm = T),
            jac_zoo = mean(jac_zoo, na.rm = T),
            jne_zoo = mean(jne_zoo, na.rm = T),
            jtu_zoo = mean(jtu_zoo, na.rm = T)
        ) 
) # ddf
summary(ens)

### Draw zonal patterns
zonal <- data.frame(ens %>% 
            group_by(y) %>% 
            summarise(mean_jtu_plankton = mean(jtu, na.rm=T), sd_plankton = sd(jtu, na.rm=T), 
            mean_jtu_phyto = mean(jtu_phyto, na.rm=T), sd_phyto = sd(jtu_phyto, na.rm=T),
            mean_jtu_zoo = mean(jtu_zoo, na.rm=T), sd_zoo = sd(jtu_zoo, na.rm=T) )    
) # eo ddf
summary(zonal)

# Define maxi and mins 
min <- 0
max <- 0.65

plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_jtu_plankton - sd_plankton, xmax = mean_jtu_plankton + sd_plankton), fill = "grey75", data = na.omit(zonal)) +
   geom_path(aes(y = y, x = mean_jtu_plankton), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(min,max)) + 
   xlab("Mean annual species turn-over") + 
   scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
   theme_classic() 

plot2 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_jtu_phyto - sd_phyto, xmax = mean_jtu_phyto + sd_phyto), fill = "grey75", data = na.omit(zonal)) +
   geom_path(aes(y = y, x = mean_jtu_phyto), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(min,max)) + 
   xlab("Mean annual species turn-over") +
   scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
   theme_classic() 
#
plot3 <- ggplot() + geom_ribbon(aes(y = y, xmin = mean_jtu_zoo - sd_zoo, xmax = mean_jtu_zoo + sd_zoo), fill = "grey75", data = na.omit(zonal)) +
   geom_path(aes(y = y, x = mean_jtu_zoo), data = na.omit(zonal), colour = "black" ) + scale_x_continuous(limits = c(min,max)) + 
   xlab("Mean annual species turn-over") + 
   scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
   theme_classic() 

# Save
setwd(WD)
ggsave(plot = plot1, filename = "plot_zonal_ann_plankton_perc_ensemble.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot2, filename = "plot_zonal_ann_phyto_perc_ensemble.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot3, filename = "plot_zonal_ann_zoo_perc_ensemble.pdf", dpi = 300, width = 3, height = 4)


### Maps
# Annual plankton jac
map_jac <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ens) +
	scale_fill_viridis(name = "Mean annual\nJaccard", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

# Annual plankton nestedness
map_jne <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = ens) +
	scale_fill_viridis(name = "Mean annual\nNestedness", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

# Annual plankton turn-over
map_jtu <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ens) +
	scale_fill_viridis(name = "Mean annual\nTurn-over", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
              labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

### Phyto
map_jac_phyto <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = ens) +
	scale_fill_viridis(name = "Mean annual\nJaccard", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

# Annual Phyto nestedness
map_jne_phyto <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = ens) +
	scale_fill_viridis(name = "Mean annual\nNestedness", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
              labels = c("GM","60°E","120°E","180°","120°W","60°W","GM") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

# Annual Phyto turn-over
map_jtu_phyto <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = ens) +
	scale_fill_viridis(name = "Mean annual\nTurn-over", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

### Zoo
map_jac_zoo <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = ens) +
	scale_fill_viridis(name = "Mean annual\nJaccard", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
              labels = c("GM","60°E","120°E","180°","120°W","60°W","GM") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

# Annual Zoo nestedness
map_jne_zoo <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = ens) +
	scale_fill_viridis(name = "Mean annual\nNestedness", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")

# Annual Zoo turn-over
map_jtu_zoo <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = ens) +
	scale_fill_viridis(name = "Mean annual\nTurn-over", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "grey70", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("GM","60°E","120°E","180°","120°W","60°W","GM") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none")

 
ggsave(plot = map_jac, filename = paste("map_annual_jac_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jne, filename = paste("map_annual_jne_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jtu, filename = paste("map_annual_jtu_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jac_phyto, filename = paste("map_annual_jac_phyto_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jne_phyto, filename = paste("map_annual_jne_phyto_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jtu_phyto, filename = paste("map_annual_jtu_phyto_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jac_zoo, filename = paste("map_annual_jac_zoo_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)		
ggsave(plot = map_jne_zoo, filename = paste("map_annual_jne_zoo_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	
ggsave(plot = map_jtu_zoo, filename = paste("map_annual_jtu_zoo_","ensemble",".pdf", sep = ""), dpi = 300, width = 7, height = 5)	


### Use ggarrange to show plankton turn-over
ggarrange(map_jtu_phyto, plot2, map_jtu_zoo, plot3, map_jtu, plot1, labels = c("A","B","C","D","E","F"), ncol = 2, nrow = 3, widths = c(2.75,1))
panel <- ggarrange(map_jtu_phyto, plot2, map_jtu_zoo, plot3, map_jtu, plot1, labels = c("a","b","c","d","e","f"), ncol = 2, nrow = 3, widths = c(2.75,1))
setwd("/net/kryo/work/fabioben/OVERSEE/data")
ggsave(plot = panel, filename = "panel_ensemble_future_jtu_Fig3A-F.pdf", dpi = 300, width = 8.5, height = 7.5)

### And just save one map with color palette at bottom for panel to be cleaner
mapwithpalette <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = ens) +
	scale_fill_viridis(name = "Mean annual\nTurn-over", limits = c(0,0.85)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom")
#
ggsave(plot = mapwithpalette, filename = "mapwithpalette.pdf", dpi = 300, width = 8, height = 5)


### 22/04/2020: Compute beta ratio: JTU/Jac to examine the contribution of SST to dissimilarity
colnames(ens)
ens$jtu_ratio <- (ens$jtu)/(ens$jac)
ens$jtu_ratio_phyto <- (ens$jtu_phyto)/(ens$jac_phyto)
ens$jtu_ratio_zoo <- (ens$jtu_zoo)/(ens$jac_zoo)
summary(ens)

summary(ens[abs(ens$y) < 30,"jtu"])
summary(ens[abs(ens$y) > 60,"jtu"])

### Make some maps
ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_ratio), data = ens) +
	scale_fill_viridis(name = "ßratio", limits = c(0.1,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_ratio), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right")
#
ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_ratio_phyto), data = ens) +
	scale_fill_viridis(name = "ßratio", limits = c(0.1,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_ratio_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right")
# 
ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_ratio_zoo), data = ens) +
	scale_fill_viridis(name = "ßratio", limits = c(0.1,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_ratio_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right")
        
### For the manuscript: median (±IQR) turn-over
### Not normally distributed so examine median rather than mean
median(ens$jtu, na.rm = T) ; IQR(ens$jtu, na.rm = T)
median(ens$jtu_phyto, na.rm = T) ; IQR(ens$jtu_phyto, na.rm = T)
median(ens$jtu_zoo, na.rm = T) ; IQR(ens$jtu_zoo, na.rm = T)          

# ANd to examine extratropical STT rates 
ens2 <- ens[ens$y < 30,]
median(ens2$jtu, na.rm = T) ; IQR(ens2$jtu, na.rm = T)
median(ens2$jtu_phyto, na.rm = T) ; IQR(ens2$jtu_phyto, na.rm = T)
median(ens2$jtu_zoo, na.rm = T) ; IQR(ens2$jtu_zoo, na.rm = T) 

ens2 <- ens[ens$y > 60,]
median(ens2$jtu, na.rm = T) ; IQR(ens2$jtu, na.rm = T)
median(ens2$jtu_phyto, na.rm = T) ; IQR(ens2$jtu_phyto, na.rm = T)
median(ens2$jtu_zoo, na.rm = T) ; IQR(ens2$jtu_zoo, na.rm = T) 

mean(ens$jtu_phyto, na.rm = T)
mean(ens$jtu_zoo, na.rm = T)

### Perform paired wilcoxon tests to check if there diff between the 2 trophic levels
# STT
subset <- na.omit(ens[,c("id","jtu_phyto","jtu_zoo")])
m <- melt(subset, id.vars = "id")
summary(m)
colnames(m) <- c("id","Group","STT")
wilcox.test(STT ~ factor(Group), data = m, paired = T) 
# V = 147890000, p-value < 2.2e-16 ; alternative hypothesis: true location shift is not equal to 0
### --> signif. variation
kruskal.test(STT ~ factor(Group), data = m)
# p-value < 2.2e-16

# Richness
subset <- na.omit(ens[,c("id","jac_phyto","jac_zoo")])
m <- melt(subset, id.vars = "id")
colnames(m) <- c("id","Group","Jac")
wilcox.test(Jac ~ factor(Group), data = m, paired = T) 
# V = 453110000, p-value < 2.2e-16
### 
kruskal.test(Jac ~ factor(Group), data = m)
# Kruskal-Wallis chi-squared = 4304.5, p-value < 2.2e-16

### And examine just fr the N Hemisphere
subset <- na.omit(ens[ens$y > 60,c("id","jtu_phyto","jtu_zoo")])
summary(subset)
m <- melt(subset, id.vars = "id")
colnames(m) <- c("id","Group","STT")
wilcox.test(STT ~ factor(Group), data = m, paired = T) 
# V = 2970900, p-value < 2.2e-16
### --> signif. variation
kruskal.test(STT ~ factor(Group), data = m)
# 


### Ensemble per SDM ---------------------------------------------------------------
ens.SDM <- data.frame(table %>%
        group_by(id,SDM) %>%
        summarize(x = unique(x), y = unique(y), 
            jac = mean(jac, na.rm = T),
            jne = mean(jne, na.rm = T),
            jtu = mean(jtu, na.rm = T),
            jac_phyto = mean(jac_phyto, na.rm = T),
            jne_phyto = mean(jne_phyto, na.rm = T),
            jtu_phyto = mean(jtu_phyto, na.rm = T),
            jac_zoo = mean(jac_zoo, na.rm = T),
            jne_zoo = mean(jne_zoo, na.rm = T),
            jtu_zoo = mean(jtu_zoo, na.rm = T)
        ) 
) # ddf
summary(ens.SDM) # dim(ens.SDM)

### Maps
for(sdm in unique(ens.SDM$SDM) ) {

    message(paste("Mapping changes in community for ",sdm, sep= ""))
    
    # Annual plankton jac
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual plankton nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual plankton turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    ### Phyto
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_phyto_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Phyto nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_phyto_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Phyto turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_phyto_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    ### Zoo
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_zoo_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)		

    # Annual Zoo nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_zoo_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Zoo turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = ens.SDM[ens.SDM$SDM == sdm,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = ens.SDM[ens.SDM$SDM == sdm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_zoo_",sdm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	    
    
} # eo for loop - sdm in SDM


### Ensemble per ESM ---------------------------------------------------------------
ens.ESM <- data.frame(table %>%
        group_by(id,ESM) %>%
        summarize(x = unique(x), y = unique(y), 
            jac = mean(jac, na.rm = T),
            jne = mean(jne, na.rm = T),
            jtu = mean(jtu, na.rm = T),
            jac_phyto = mean(jac_phyto, na.rm = T),
            jne_phyto = mean(jne_phyto, na.rm = T),
            jtu_phyto = mean(jtu_phyto, na.rm = T),
            jac_zoo = mean(jac_zoo, na.rm = T),
            jne_zoo = mean(jne_zoo, na.rm = T),
            jtu_zoo = mean(jtu_zoo, na.rm = T)
        ) 
) # ddf
summary(ens.ESM)

### Maps
for(esm in unique(ens.ESM$ESM) ) {

    message(paste("Mapping changes in community for ",esm, sep= ""))
    
    # Annual plankton jac
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual plankton nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual plankton turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    ### Phyto
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Phyto nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Phyto turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    ### Zoo
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)		

    # Annual Zoo nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Zoo turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = ens.ESM[ens.ESM$ESM == esm,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = ens.ESM[ens.ESM$ESM == esm,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	    
    
} # eo for loop - esm in ESM


### Ensemble per pred pool ---------------------------------------------------------------
ens.pool <- data.frame(table %>%
        group_by(id,pool) %>%
        summarize(x = unique(x), y = unique(y), 
            jac = mean(jac, na.rm = T),
            jne = mean(jne, na.rm = T),
            jtu = mean(jtu, na.rm = T),
            jac_phyto = mean(jac_phyto, na.rm = T),
            jne_phyto = mean(jne_phyto, na.rm = T),
            jtu_phyto = mean(jtu_phyto, na.rm = T),
            jac_zoo = mean(jac_zoo, na.rm = T),
            jne_zoo = mean(jne_zoo, na.rm = T),
            jtu_zoo = mean(jtu_zoo, na.rm = T)
        ) 
) # ddf
summary(ens.pool)

### Maps
for(p in unique(ens.pool$pool) ) {

    message(paste("Mapping changes in community for ",p, sep= ""))
    
    # Annual plankton jac
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual plankton nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual plankton turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    ### Phyto
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_phyto_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Phyto nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_phyto_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Phyto turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_phyto_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    ### Zoo
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nJaccard", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jac_zoo_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)		

    # Annual Zoo nestedness
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nNestedness", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jne_zoo_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

    # Annual Zoo turn-over
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = ens.pool[ens.pool $pool == p,]) +
    	scale_fill_viridis(name = "Annual\nTurn-over", limits = c(0,1)) +
        geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = ens.pool[ens.pool $pool == p,]) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    # 
    ggsave(plot = map, filename = paste("map_annual_jtu_zoo_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	    
    
} # eo for loop - esm in ESM



# ----------------------------------------------------------

### 5°) Quantify variance due to pool vs ESM vs SDM
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("table_ann_changes_",dir())]; files
# Separate the ones for richness (no "beta.div") from those with beta.div
files2 <- files[grep("beta.div", files)]; files2
files1 <- files[!(files %in% files2)]; files1
require("parallel")
# f <- files1[1]
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
dim(table)
head(table)
rm(res); gc()
# summary(table)
fit <- aov(diff_rich_tot ~ ESM*SDM*pool, data = table)
str(summary(fit))
ssq <- summary(fit)[[1]][,2]
ssq/sum(ssq)

### Or explore variance with violin plots
library("wesanderson")
library("ggthemes")

plot1 <- ggplot(aes(x = factor(SDM), y = perc_rich_tot, fill = factor(SDM)), data = table) + geom_violin(colour = "black") + 
        geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
        scale_fill_manual(name = "", values = wes_palettes$Zissou1) + 
        scale_y_continuous(limits = c(-80,150)) + geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("SDM") + ylab("Annual % change in plankton richness") + theme_classic()
#
plot2 <- ggplot(aes(x = factor(ESM), y = perc_rich_tot, fill = factor(ESM)), data = table) + geom_violin(colour = "black") + 
        geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
        scale_fill_manual(name = "", values = wes_palettes$Zissou1) + 
        scale_y_continuous(limits = c(-80,150)) + geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("ESM") + ylab("Annual % change in plankton richness") + theme_classic()
#
plot3 <- ggplot(aes(x = factor(pool), y = perc_rich_tot, fill = factor(pool)), data = table) + geom_violin(colour = "black") + 
        geom_boxplot(colour = "black", fill = "white", width = 0.1) + 
        scale_fill_manual(name = "", values = wes_palettes$Zissou1) + 
        scale_y_continuous(limits = c(-80,150)) + geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("Pool") + ylab("Annual % change in plankton richness") + theme_classic()
#
ggsave(plot = plot1, filename = "plot_violins_ann_perc_tot_SDM.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = plot2, filename = "plot_violins_ann_perc_tot_ESM.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = plot3, filename = "plot_violins_ann_perc_tot_pool.jpg", dpi = 300, width = 7, height = 5)

            
# ----------------------------------------------------------

### 6°) Examine covariance between changes in diversity (richness/composition)

### A) Retrieve changes in richness 
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("table_ann_changes_",dir())]; files
# Separate the ones for richness (no "beta.div") from those with beta.div
files2 <- files[grep("beta.div", files)]; files2
files1 <- files[!(files %in% files2)]; files1

require("parallel")
# f <- files1[1]
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
table1 <- dplyr::bind_rows(res)
dim(table1)
head(table1)
rm(res); gc()

### B) Retrieve changes in composition 
res <- mclapply(files2, function(f) {
            # Useless message
            message(paste("Loading results from ",f, sep = ""))
            t <- get(load(f))
            # Get attributes from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            t$ESM <- terms[6,1] 
            t$SDM <- terms[7,1] 
            p <- terms[8,1]
            t$pool <- str_replace(as.character(p),".Rdata","")
            # Return
            return(t)
    }, mc.cores = 25
) # eo mclapply
# Rbind
table2 <- dplyr::bind_rows(res)
dim(table2)
head(table2)
rm(res); gc()


### 23/03/2020: Examine covariance between changes in richness and turn-over
colnames(table1); colnames(table2)
dim(table1); dim(table2)
### Compute mean ensembles
ens1 <- data.frame(table1 %>%
        group_by(id) %>%
        summarize(x = unique(x), y = unique(y), 
            rich_phyto = mean(perc_rich_phyto, na.rm = T),
            rich_zoo = mean(perc_rich_zoo, na.rm = T)
        ) 
) # ddf
summary(ens1)

ens2 <- data.frame(table2 %>%
        group_by(id) %>%
        summarize(x = unique(x), y = unique(y), 
            jtu_phyto = mean(jtu_phyto, na.rm = T),
            jtu_zoo = mean(jtu_zoo, na.rm = T)
        ) 
) # ddf
summary(ens2)
# dim(ens1); dim(ens2)
# Combine 
commons <- intersect(unique(ens1$id), unique(ens2$id)) # length(commons)
ddf <- cbind(ens1[which(ens1$id %in% commons),], ens2[which(ens2$id %in% commons),c("jtu_phyto","jtu_zoo")])
summary(ddf); dim(ddf)
# Plots JTU vs. %rich
p1 <- ggplot(data = ddf) + geom_point(aes(x = rich_phyto, y = jtu_phyto, colour = jtu_phyto), alpha = 0.5) + 
    scale_colour_viridis(name = "", option = "viridis") + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    xlab("Mean annual %∆SR") + ylab("Mean annual species turn-over") +
    scale_x_continuous(limits = c(-35,119)) + scale_y_continuous(limits = c(0,0.85)) + 
    theme_classic()
# Zooplankton
p2 <- ggplot(data = ddf) + geom_point(aes(x = rich_zoo, y = jtu_zoo, colour = jtu_zoo), alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    scale_colour_viridis(name = "", option = "viridis") + 
    xlab("Mean annual %∆SR") + ylab("Mean annual species turn-over") +
    scale_x_continuous(limits = c(-35,119)) + scale_y_continuous(limits = c(0,0.85)) + 
    theme_classic()
#
ggsave(plot = p1, filename = "plot_perc.vs.jtu_phyto_annual_ensembles.jpg", dpi = 300, width = 5, height = 4.5)
ggsave(plot = p2, filename = "plot_perc.vs.jtu_zoo_annual_ensembles.jpg", dpi = 300, width = 5, height = 4.5)


### 21/04/2020: Combine those with the 2 zonal plots and the 2 maps of turn-over in a panel


# ----------------------------------------------------------

### D) Retrieve diff climatologies 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
files <- dir()[grep("mon_diff",dir())]; files
# d <- get(load(files[5]))
f <- files[5]
res <- mclapply(files, function(f) {
            # Useless message
            message(paste("Loading diffs from ",f, sep = ""))
            d <- get(load(f))
            # Get attributes from filename
            terms <- do.call(cbind, strsplit(as.character(f),"_")) # terms
            ESM <- terms[6,1] 
            var <- terms[4,1] 
            # Compute mean annual diff
            d$Annual <- rowMeans(as.matrix(d[,c(4:15)])) # summary(d$Annual)
            # Return
            return( data.frame(id = d$id, x = d$x, y = d$y, mean = d$Annual, ESM = ESM, var = var) )
    }, mc.cores = 30
) # eo mclapply
# Rbind
diffs <- dplyr::bind_rows(res)
dim(diffs)
head(diffs)
rm(res); gc()
unique(diffs$var) # missing range sst because annual and not monthly
unique(diffs$ESM)

### Merge with table1/table2 which contains the mean annual change in richness/ composition
# Per ESM (for facet)
ens.ESM <- data.frame(table1 %>%
        group_by(id,ESM) %>%
        summarize(x = unique(x), y = unique(y), 
            rich_tot = mean(diff_rich_tot, na.rm = T),
            rich_phyto = mean(diff_rich_phyto, na.rm = T),
            rich_zoo = mean(diff_rich_zoo, na.rm = T)
        ) 
) # ddf
# head(ens.ESM) ; head(diffs)

jtu <- data.frame(table2 %>%
        group_by(id,ESM) %>%
        summarize(x = unique(x), y = unique(y), 
            jtu_phyto = mean(jtu_phyto, na.rm = T),
            jtu_zoo = mean(jtu_zoo, na.rm = T)
        ) 
) # ddf
summary(jtu)

#var.phyto <- c("sst","logchl","logno3","logsio2","nstar","sistar")
#var.zoo <- c("sst","o2","logchl","logno3","logsio2","nstar","sistar")

### In for loops, plot biplot of %∆SR and diffs in variables
# esm <- "CESM-BEC"

setwd("/net/kryo/work/fabioben/OVERSEE/data/")

for(esm in unique(ens.ESM$ESM)) {
    
        # Plotting for...
        message(paste("Making plots for ",esm, sep = ""))
        ens2 <- ens.ESM[which(ens.ESM$ESM == esm),]
        diff <- diffs[which(diffs$ESM == esm),]
        commons <- intersect(unique(ens2$id), unique(diff$id)) # length(commons)
        # head(ens2); head(diff)
        ens3 <- ens2[which(ens2$id %in% commons),]
        
        for(v in var.phyto) {
            
            # v <- "sst"
            message(paste("Plotting relationship between changes in phytoplankton diversity and ",v, sep = ""))
            diff2 <- diff[which(diff$id %in% commons & diff$var == v),]
            diff2 <- diff2[order(diff2$id),]
            ens3 <- ens3[order(ens3$id),]
            ens3[,v] <- diff2$mean # head(ens3)
            # Add domain to examine per lat band
            ens3$domain <- NA
            ens3[which(abs(ens3$y) <= 30 & abs(ens3$y) ),"domain"] <- "Tropical"
            ens3[which(abs(ens3$y) > 30 & abs(ens3$y) <= 60),"domain"] <- "Temperate"
            ens3[which(abs(ens3$y) > 30 & abs(ens3$y) > 60),"domain"] <- "Polar"
            ens3$inter1 <- 0
            ens3$inter2 <- 0
            # levels(factor(ens3$domain))
            p <- ggplot() + geom_point(aes(x = get(v), y = rich_phyto, colour = factor(domain)), alpha = 0.5, data = ens3) + 
                    scale_colour_manual(name="", values = c("#2166ac","#4393c3","#92c5de")) + 
                    geom_hline(aes(yintercept = unique(inter1)), linetype = "dashed", data = ens3) + 
                    geom_vline(aes(xintercept = unique(inter2)), linetype = "dashed", data = ens3) +
                    xlab(v) + ylab("Mean annual %∆SR") + theme_classic() +
                    facet_wrap(~ factor(ens3$domain), ncol = 3)
                    
            ggsave(plot = p, filename = paste("plot_annual_perc_phyto_vs_",v,"_",esm,".jpg", sep = ""), dpi = 300, width = 7.5, height = 3)
            
        } # eo 2nd for loop
        
        for(v in var.zoo) {
            
            message(paste("Plotting relationship between changes in zooplankton diversity and ",v, sep = ""))
            diff2 <- diff[which(diff$id %in% commons & diff$var == v),]
            diff2 <- diff2[order(diff2$id),]
            ens3 <- ens3[order(ens3$id),]
            ens3[,v] <- diff2$mean # head(ens3)
            # Add domain to examine per lat band
            ens3$domain <- NA
            ens3[which(abs(ens3$y) <= 30 & abs(ens3$y) ),"domain"] <- "Tropical"
            ens3[which(abs(ens3$y) > 30 & abs(ens3$y) <= 60),"domain"] <- "Temperate"
            ens3[which(abs(ens3$y) > 30 & abs(ens3$y) > 60),"domain"] <- "Polar"
            ens3$inter1 <- 0
            ens3$inter2 <- 0
            
            p <- ggplot() + geom_point(aes(x = get(v), y = rich_zoo, colour = factor(domain)), alpha = 0.5, data = ens3) + 
                    scale_colour_manual(name="", values = c("#2166ac","#4393c3","#92c5de")) +
                    geom_hline(aes(yintercept = unique(inter1)), linetype = "dashed", data = ens3) + 
                    geom_vline(aes(xintercept = unique(inter2)), linetype = "dashed", data = ens3) +
                    xlab(v) + ylab("Mean annual %∆SR") + theme_classic() +
                    facet_wrap(~ factor(ens3$domain), ncol = 3)
                    
            ggsave(plot = p, filename = paste("plot_annual_perc_zoo_vs_",v,"_",esm,".jpg", sep = ""), dpi = 300, width = 7.5, height = 3)
           
        } # eo 3rd for loop
        
         message(paste("", sep = ""))
    
} # eo 1st for loop



### 24/03/2020: Examine covariation through a PCA (how original, I know!):
# - in a for loop, per ESM, combine diff with jtu and perc rich for phyto and zoo, also combine annual delta SST range
# - perform PCA (scale.unit = T), save plots as nice as possible
# - map 4 PCs
library("FactoMineR")
vars <- c("sst","o2","logchl","logno3","logsio2","nstar","sistar")
# for dSST: clims_ann_diff_dsst_rcp85_CESM-BEC_2031-2100.Rdata
esm <- "IPSL-PISCES"

for(esm in unique(ens.ESM$ESM)) {
    
        message(paste("Performing PCA for ",esm, sep = ""))
        ens2 <- ens.ESM[which(ens.ESM$ESM == esm),]
        diff <- diffs[which(diffs$ESM == esm),]
        commons <- intersect(unique(ens2$id), unique(diff$id)) # length(commons)
        # Select
        ens3 <- ens2[which(ens2$id %in% commons),]
        ens3 <- ens3[order(ens3$id),]
        # Provide vars
        for(v in vars) {
            diff2 <- diff[which(diff$id %in% commons & diff$var == v),]
            diff2 <- diff2[order(diff2$id),]
            ens3[,v] <- diff2[diff2$var == v,"mean"]
        }
        
        # Provide dSST diff
        setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
        dsst <- get(load(paste("clims_ann_diff_dsst_rcp85_",esm,"_2031-2100.Rdata", sep = "")))
        dsst <- dsst[order(dsst$id),]
        dsst2 <- dsst[dsst$id %in% commons,] # dim(dsst2)
        ens3$dsst <- dsst2$dSST
        rm(dsst,diff2); gc()
        setwd("/net/kryo/work/fabioben/OVERSEE/data/")
        # Add JTU estimates
        jtu2 <- jtu[jtu$id %in% commons & jtu$ESM == esm,]
        jtu2 <- jtu2[order(jtu2$id),] # dim(jtu2)
        commons2 <- intersect(unique(ens3$id), unique(jtu2$id)) # length(commons2)
        ens4 <- ens3[ens3$id %in% commons2,] # dim(ens4); colnames(ens4)
        ens4$jtu_phyto <- jtu2$jtu_phyto
        ens4$jtu_zoo <- jtu2$jtu_zoo
        # Perform PCA
        # colnames(ens4)
        data4pca <- na.omit(ens4[,c(1:4,6:17)])
        # dim(data4pca); summary(data4pca); colnames(data4pca)
        pca <- PCA(X = data4pca[,c(4:length(data4pca))], ncp = 5, scale.unit = T, graph = F, quanti.sup = c(2,3,12,13))
        # summary(pca) ; str(pca)
        # Provide for data4pca for mapping
        data4pca$PC1 <- pca$ind$coord[,1]
        data4pca$PC2 <- pca$ind$coord[,2]
        data4pca$PC3 <- pca$ind$coord[,3]
        data4pca$PC4 <- pca$ind$coord[,4]
        pc1 <- ceiling(pca$eig[1,2])
        pc2 <- ceiling(pca$eig[2,2])
        pc3 <- ceiling(pca$eig[3,2])
        pc4 <- ceiling(pca$eig[4,2])
        
        # Save maps
        message(paste("Plotting results of the PCA for ",esm, sep = ""))
        p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC1), data = data4pca) +
            geom_contour(colour = "grey60", binwidth = 1.5, size = 0.25, aes(x = x, y = y, z = PC1), data = data4pca) +
 	        scale_fill_gradient2(name = paste("PC1 (",pc1,"%)", sep=""), low = "#3288bd", high = "#d53e4f", mid = "white") +
 	        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                  labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		          labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                  
        #
        p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC2), data = data4pca) +
            geom_contour(colour = "grey60", binwidth = 1.5, size = 0.25, aes(x = x, y = y, z = PC2), data = data4pca) +
         	scale_fill_gradient2(name = paste("PC2 (",pc2,"%)", sep=""), low = "#3288bd", high = "#d53e4f", mid = "white") +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC3), data = data4pca) +
            geom_contour(colour = "grey60", binwidth = 1.5, size = 0.25, aes(x = x, y = y, z = PC3), data = data4pca) +
         	scale_fill_gradient2(name = paste("PC3 (",pc3,"%)", sep=""), low = "#3288bd", high = "#d53e4f", mid = "white") +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC4), data = data4pca) +
            geom_contour(colour = "grey60", binwidth = 1.5, size = 0.25, aes(x = x, y = y, z = PC4), data = data4pca) +
         	scale_fill_gradient2(name = paste("PC4 (",pc4,"%)", sep=""), low = "#3288bd", high = "#d53e4f", mid = "white") +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                
        ggsave(plot = p1, filename = paste("map_","PC1","_",esm,"v2.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
        ggsave(plot = p2, filename = paste("map_","PC2","_",esm,"v2.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
        ggsave(plot = p3, filename = paste("map_","PC3","_",esm,"v2.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
        ggsave(plot = p4, filename = paste("map_","PC4","_",esm,"v2.jpg", sep = ""), dpi = 300, width = 7, height = 5)	
        
        # Save PCA plot (PC1 vs PC2 and PC3 vs PC4)
        # Extract data from a FactoMineR::PCA object
        require("ggrepel"); require("wesanderson"); require("ggthemes")
        augment.PCA <- function(x, dims=c(1,2,3,4), which="col") {
          .get <- function(x, element, dims) {
            y <- as.data.frame(x[[element]]$coord[,dims])
            if (nrow(y) == 0) {
              y <- NULL
            } else {
              y$type <- element
            }
            return(y)
          }
          if (which == "col") {
            y <- rbind(.get(x, "var", dims), .get(x, "quanti.sup", dims))
          } else {
            y <- rbind(.get(x, "ind", dims), .get(x, "quali.sup", dims))
          }
          y$var <- row.names(y)
          row.names(y) <- NULL
          return(y)
        }
        pcad <- augment.PCA(pca)
        # class(pcad); str(pcad) ; pcad$var ; pcad$type
        # Rename vars and types of vars for colours in plot
        pcad$var <- c("Lat","SST","O2","log(Chl)","log(NO3)","log(SiOH4)","N*","Si*","rSST","%∆Phyto","%∆Zoo","T-O Phyto","T-O Zoo")
        pcad$type <- c("geo","env","env","env","env","env","env","env","env","bio","bio","bio","bio")
        # levels(factor(pcad$type))
        # Make plots
        plot1 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
          annotate(geom="path", colour="black", x=cos(seq(0,2*pi,length.out=100)), y=sin(seq(0,2*pi,length.out=100))) +
          geom_segment(aes(x=0, xend=Dim.1, y=0, yend=Dim.2, colour=factor(type)), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
          scale_colour_manual(name = "", values = c("#F21A00","#3B9AB2","#E1AF00") ) + 
          geom_text_repel(aes(x=Dim.1, y=Dim.2, colour=type, label=var), 
                  data=filter(pcad,(Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
          xlab(paste("PC1 (",pc1,"%)",sep="")) + ylab(paste("PC2 (",pc2,"%)",sep="")) + theme_bw()
        
        plot2 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
          annotate(geom="path", colour="black", x=cos(seq(0,2*pi,length.out=100)), y=sin(seq(0,2*pi,length.out=100))) +
          geom_segment(aes(x=0, xend=Dim.3, y=0, yend=Dim.4, colour=factor(type)), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
          scale_colour_manual(name = "", values = c("#F21A00","#3B9AB2","#E1AF00") ) + 
          geom_text_repel(aes(x=Dim.3, y=Dim.4, colour=type, label=var),
                  data=filter(pcad,(Dim.3^2+Dim.4^2) > 0.2^2), segment.alpha=0.5) +
          xlab(paste("PC3 (",pc3,"%)",sep="")) + ylab(paste("PC4 (",pc4,"%)",sep="")) + theme_bw()
          
         ggsave(plot = plot1, filename = paste("plot_","PC1xPC2","_",esm,"v2.jpg", sep = ""), width = 5, height = 5)
         ggsave(plot = plot2, filename = paste("plot_","PC3xPC4","_",esm,"v2.jpg", sep = ""), width = 5, height = 5)
        
         rm(pcad, data4pca, ens4, ens3); gc()
        
} # eo for loop - esm in ESMs
    



# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


