
##### 19/08/19 - ETHZ - Fabio Benedetti © UP Group+ IBP+ ETH Zürich
##### Script to prepare the figures for the mansucript, as follows:
#	- Figure 1: panel of contemporary annual diversity estimates (total plankton, phyto-, zooplankton) 
#   - Figure 2: anomalies of log(S) to eV for phyto- and zooplankton
#   - Figure 3: heatmap of rank correlations (or slopes of lm) between diversity estimates (total, phyto-, zooplankton + anoms to eV) 
#       and selected predictors
#   - Figure 4: panel of future diversity changes 
 
### Last update: 30/10/19

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("raster")
library("rgdal")
library("sp")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("vegan")
library("FactoMineR")
library("cmocean")

# Coastline
world2 <- map_data("world2")

WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Figure 1: panel of contemporary annual diversity estimates (total plankton, phyto-, zooplankton) 
setwd(paste(WD,"/tables_composition","/", sep=""))
phyto.base <- read.table("table_phyto_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
zoo.base <- read.table("table_zoo_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
base <- cbind(phyto.base,zoo.base[,c(4:length(zoo.base))])
dim(base)
colnames(base) # phyto = 4:341  ; zoo = 342:length(base)
base$x2 <- base$x 
base[base$x < 0 ,"x2"] <- (base[base$x < 0 ,"x"]) + 360
### Compute species richness for total plankton, phyto- and zooplankton 
base$rich_plankton <- rowSums(as.matrix(base[,c(4:865)]))
base$rich_phyto <- rowSums(as.matrix(base[,c(4:341)]))
base$rich_zoo <- rowSums(as.matrix(base[,c(342:865)]))
summary(base)
setwd(WD)

### Bin all diversity indices
# 1) HSI_plankton
base$rich_plankton_bin <- factor(cut_interval(base$rich_plankton,10))
levels(base$rich_plankton_bin)
levels <- str_replace_all(levels(base$rich_plankton_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(base$rich_plankton_bin) <- levels

# 2) HSI_phyto
base$rich_phyto_bin <- factor(cut_interval(base$rich_phyto,10))
levels(base$rich_phyto_bin)
levels <- str_replace_all(levels(base$rich_phyto_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(base$rich_phyto_bin) <- levels

# 3) HSI_zoo
base$rich_zoo_bin <- factor(cut_interval(base$rich_zoo,10))
levels(base$rich_zoo_bin)
levels <- str_replace_all(levels(base$rich_zoo_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(base$rich_zoo_bin) <- levels

map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_plankton_bin)), data = base) +
 	scale_fill_viridis(name = "", na.value = "white", discrete = T) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_phyto_bin)), data = base) +
 	scale_fill_viridis(name = "", na.value = "white", discrete = T) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_zoo_bin)), data = base) +
 	scale_fill_viridis(name = "", na.value = "white", discrete = T) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
ggsave(plot = map1, filename = "map_Fig1A_plankton_rich_baseline.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = map2, filename = "map_Fig1C_phyto_rich_baseline.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = map3, filename = "map_Fig1E_zoo_rich_baseline.pdf", dpi = 300, width = 7, height = 4)


### Make zonal plots too
library("dplyr")
zonal <- data.frame(base %>% 
            group_by(y) %>% 
            summarise(mean_plankton = mean(rich_plankton,na.rm=T), sd_plankton = sd(rich_plankton,na.rm=T), 
            mean_phyto = mean(rich_phyto,na.rm=T), sd_phyto = sd(rich_phyto,na.rm=T),
            mean_zoo = mean(rich_zoo,na.rm=T), sd_zoo = sd(rich_zoo,na.rm=T) )    
) # eo ddf
summary(zonal)

plot1 <- ggplot() + geom_ribbon(aes(x = y, ymin = mean_plankton - sd_plankton, ymax = mean_plankton + sd_plankton), 
            fill = "grey70", data = zonal) +
            geom_line(aes(x = y, y = mean_plankton), data = zonal, colour = "black" ) + scale_y_continuous(limits = c(0,300)) + 
		    ylab("Mean annual species richness\n(Plankton)") + xlab("Latitude (°)") + 
		    theme_classic() + coord_flip()
#
plot2 <- ggplot() + geom_ribbon(aes(x = y, ymin = mean_phyto - sd_phyto, ymax = mean_phyto + sd_phyto), 
            fill = "grey70", data = zonal) +
            geom_line(aes(x = y, y = mean_phyto), data = zonal, colour = "black" ) + scale_y_continuous(limits = c(0,200)) + 
		    ylab("Mean annual species richness\n(Phytoplankton)") + xlab("Latitude (°)") + 
		    theme_classic() + coord_flip()
#
plot3 <- ggplot() + geom_ribbon(aes(x = y, ymin = mean_zoo - sd_zoo, ymax = mean_zoo + sd_zoo),
            fill = "grey70", data = zonal) +
            geom_line(aes(x = y, y = mean_zoo), data = zonal, colour = "black" ) + scale_y_continuous(limits = c(0,200)) + 
		    ylab("Mean annual species richness\n(Zooplankton)") + xlab("Latitude (°)") + 
		    theme_classic() + coord_flip()
#
ggsave(plot = plot1, filename = "plot_Fig1B_plankton_rich_baseline.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot2, filename = "plot_Fig1D_phyto_rich_baseline.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot3, filename = "plot_Fig1F_zoo_rich_baseline.pdf", dpi = 300, width = 3, height = 4)


# -----------------------------------------------------------

### 2°) Figure 2: panel of anomalies of log(S) to eV, for phyto- and zooplankton, with:
#   - biplot with ln(S) ~ eV x2
#   - map of anomaly of ln(S) to eV (lm)

# First, get annual clims of env predictors
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
files <- dir()[grep("21_02_19",dir())]
clims <- lapply(files, function(f) {
				d <- read.table(f, h = T, sep = ";")
				return(d)
		} # eo FUN
) # eo lapply
clims <- dplyr::bind_rows(clims)
# And use dplyr to compute annual climatologies
clims$x2 <- clims$x 
clims[clims$x < 0 ,"x2"] <- (clims[clims$x < 0 ,"x"]) + 360
clims$id <- factor(paste(clims$x2, clims$y, sep = "_"))
aclim <- data.frame(clims %>% group_by(id) %>% 
				summarise(x = unique(x2), y = unique(y), Bathy = mean(Bathy,na.rm=T),
				SST = mean(SST,na.rm=T), SSS = mean(SSS,na.rm=T), dSST = mean(deltaT,na.rm=T), Wind = mean(Wind,na.rm=T),
				logEKE = mean(logEKE,na.rm=T), MLD = mean(MLD1,na.rm=T), PAR = mean(PAR,na.rm=T), O2 = mean(dO2,na.rm=T), 
				NO3 = mean(logNO3,na.rm=T), SiO2 = mean(logSiO2,na.rm=T), Nstar = mean(Nstar,na.rm=T), Sistar = mean(Sistar,na.rm=T),
				Chla = mean(logChl,na.rm=T)
			) # eo summarise	
) # eo ddf
# Exclude non opean ocean cells: SSS > 30 and Bathy < -175
aclim <- aclim[aclim$SSS >= 20,]
aclim <- aclim[aclim$Bathy <= -175,]
# Check
aclim <- aclim[!is.na(aclim$id),]
aclim <- aclim[order(aclim$id),]

setwd(WD)

### Combine with 'base'
base$cell_id <- factor(paste(base$x2, base$y, sep = "_"))
base <- base[order(base$cell_id),]
head(base[,c(1:3)]);head(aclim[,c(1:3)])
# dim(base); dim(aclim)
# colnames(base)
ddf <- cbind(aclim, base[which(base$cell_id %in% unique(aclim$id)),c(867:869)])
# colnames(ddf)
# summary(ddf)

### Add K and compute eV
ddf$absT <- ddf$SST + 273.15
ddf$kT <- 1 / (ddf$absT*-8.6173303*10^(-5) )

### Compute div anomalies to SST based on a linear regression
lm <- lm(log(rich_phyto) ~ kT, data = ddf, na.action = na.exclude)
summary(lm) # Adjusted R-squared = 0.7738; p-value: < 2.2e-16
ddf$phyto_anom <- residuals(lm)
# And anoms between rich and SST
lm2 <- lm(rich_phyto ~ SST, data = ddf, na.action = na.exclude)
summary(lm2) # Adjusted R-squared = 0.664; p-value: < 2.2e-16
ddf$phyto_anom2 <- residuals(lm2)

lm <- lm(log(rich_zoo) ~ kT, data = ddf, na.action = na.exclude)
summary(lm) # Adjusted R-squared = 0.849; p-value: < 2.2e-16
ddf$zoo_anom <- residuals(lm)
# And anoms between rich and SST
lm2 <- lm(rich_zoo ~ SST, data = ddf, na.action = na.exclude)
summary(lm2) # Adjusted R-squared = 0.831; p-value: < 2.2e-16
ddf$zoo_anom2 <- residuals(lm2)

summary(ddf)

#
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = phyto_anom), data = ddf[!is.na(ddf$phyto_anom),] ) +
 	scale_fill_gradient2(name = "Anomalies to\nlinear model",
    low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-0.54,1.35)) +
    geom_contour(colour = "grey60", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = phyto_anom), 
            data = ddf[!is.na(ddf$phyto_anom),] ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = zoo_anom), data = ddf[!is.na(ddf$zoo_anom),]) +
 	scale_fill_gradient2(name = "Anomalies to\nlinear model",
    low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-0.54,1.35)) +
    geom_contour(colour = "grey60", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = zoo_anom), 
            data = ddf[!is.na(ddf$zoo_anom),] ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
setwd(WD)
ggsave(plot = map4, filename = "map_Fig2B_phyto_anoms.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = map5, filename = "map_Fig2D_zoo_anoms.pdf", dpi = 300, width = 7, height = 4)


### And biplots
summary(log(ddf$rich_phyto)) # 3-5.5
summary(log(ddf$rich_zoo)) # 3.5-5.5
plot4 <- ggplot() + geom_point(aes(x = kT, y = log(rich_phyto)), data = ddf, colour = "grey70") +
            geom_smooth(aes(x = kT, y = log(rich_phyto)), data = ddf, colour = "black", method = "lm") + 
            scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness\n(ln)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot5 <- ggplot() + geom_point(aes(x = kT, y = log(rich_zoo)), data = ddf, colour = "grey70") +
            geom_smooth(aes(x = kT, y = log(rich_zoo)), data = ddf, colour = "black", method = "lm") + 
            scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness\n(ln)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
ggsave(plot = plot4, filename = "plot_Fig2A_phyto_anoms.pdf", dpi = 300, width = 4.5, height = 4)
ggsave(plot = plot5, filename = "plot_Fig2C_zoo_anoms.pdf", dpi = 300, width = 4.5, height = 4)


### 27/08/2019: Also map and plot anomalies to SST
mp1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = phyto_anom2), data = ddf[!is.na(ddf$phyto_anom2),] ) +
 	scale_fill_gradient2(name = "Anomalies to\nSST",
    low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-56,162)) +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = phyto_anom2), 
            data = ddf[!is.na(ddf$phyto_anom2),] ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mp2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = zoo_anom2), data = ddf[!is.na(ddf$zoo_anom2),] ) +
 	scale_fill_gradient2(name = "Anomalies to\nSST",
    low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-56,162)) +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = zoo_anom2), 
            data = ddf[!is.na(ddf$zoo_anom2),] ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
setwd(WD)
ggsave(plot = mp1, filename = "map_Fig2B_phyto_anoms2.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = mp2, filename = "map_Fig2D_zoo_anoms2.pdf", dpi = 300, width = 7, height = 4)


### 05/09/19: Graphs of ln(S)~eV not convincing for phytoplankton (too different looking from Righetti et al. 2019): 
### fit 2nd and 3r order polynomials and plot

### For phytoplankton
lm.2 <- lm(log(rich_phyto) ~ kT + I(kT^2), data = ddf, na.action = na.exclude)
lm.3 <- lm(log(rich_phyto) ~ kT + I(kT^2) + I(kT^3), data = ddf, na.action = na.exclude)
summary(lm.2) # Adjusted R-squared = 0.825; p-value: < 2.2e-16
summary(lm.3) # Adjusted R-squared = 0.8371; p-value: < 2.2e-16
ddf$phyto_anom_2d <- residuals(lm.2)
ddf$phyto_anom_3d <- residuals(lm.3)
ddf$phyto_2d <- predict(lm.2)
ddf$phyto_3d <- predict(lm.3)

### For zoo
lm.2 <- lm(log(rich_zoo) ~ kT + I(kT^2), data = ddf, na.action = na.exclude)
lm.3 <- lm(log(rich_zoo) ~ kT + I(kT^2) + I(kT^3), data = ddf, na.action = na.exclude)
summary(lm.2) # Adjusted R-squared = 0.849; p-value: < 2.2e-16
summary(lm.3) # Adjusted R-squared = 0.9705; p-value: < 2.2e-16
ddf$zoo_anom_2d <- residuals(lm.2)
ddf$zoo_anom_3d <- residuals(lm.3)
ddf$zoo_2d <- predict(lm.2)
ddf$zoo_3d <- predict(lm.3)

### Make plots
plot1 <- ggplot() + geom_point(aes(x = kT, y = log(rich_phyto)), data = ddf, colour = "grey70") +
            geom_smooth(aes(x = kT, y = log(rich_phyto)), data = ddf, colour = "black", method = "loess") + 
            geom_line(aes(x = kT, y = phyto_2d), data = ddf, colour = "#b2df8a") +
            geom_line(aes(x = kT, y = phyto_3d), data = ddf, colour = "#33a02c") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness\n(log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot2 <- ggplot() + geom_point(aes(x = kT, y = log(rich_zoo)), data = ddf, colour = "grey70") +
            geom_smooth(aes(x = kT, y = log(rich_zoo)), data = ddf, colour = "black", method = "loess") + 
            geom_line(aes(x = kT, y = zoo_2d), data = ddf, colour = "#a6cee3") +
            geom_line(aes(x = kT, y = zoo_3d), data = ddf, colour = "#1f78b4") + 
            scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness\n(log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
ggsave(plot = plot1, filename = "plot_phyto_anoms_withloess.pdf", dpi = 300, width = 6, height = 7)
ggsave(plot = plot2, filename = "plot_zoo_anoms_withloess.pdf", dpi = 300, width = 6, height = 7)


### Summarize slopes across the regimes defined by Righetti et al. 2019 :
# Above 23°C: 0.2939  (0.26 above 24°C)
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST >= 23,]))
# Above 19°C:  0.37812
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST >= 19,]))
# 19°C-11°C: 0.3886
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST < 19 & ddf$SST > 11,]))
# Below 11°C: 0.0054
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST <= 11,]))
# Below 9°C: -0.027
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST <= 9,]))
# Below 5°C: -0.127
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST <= 5,]))
# Below 7°C: -0.0599
summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST <= 7,]))


### Same with zoo
# Above 23°C: -0.2654 ! (-0.3412 above 24°C)
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST >= 23,]))
# Above 19°C: -0.0810 
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST >= 19,]))
# 19°C-11°C: 0.7569
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 19 & ddf$SST > 11,]))
# 17°C-11°C: 0.6877
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 17 & ddf$SST > 11,]))
# Below 11°C: 0.0222
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 11,]))
# Below 9°C: -0.0785
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 9,]))
# Below 7°C: -0.1775
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 7,]))
# Below 5°C: -0.3172
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 5,]))

### For phytoplankton, display the slopes from lm objects of varying SST thresholds so you identify the point where it goes below 0.32
for(t in seq(from = -1, to = 30, by = 1)) {
     lm <- lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST <= t,])
     message(paste("Below t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
     message(paste("Below t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
     message(paste("", sep = ""))
} 

# Same with zoo
for(t in seq(from = -1, to = 30, by = 1)) {
      lm <- lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= t,])
      message(paste("Below t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Below t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 

summary(lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST <= 9,]))
summary(lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST <= 9,]))


### And when above the threshold t 
for(t in seq(from = -1, to = 30, by = 1)) {
     lm <- lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST >= t,])
     message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
     message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
     message(paste("", sep = ""))
} 

# Same with zoo
for(t in seq(from = -1, to = 30, by = 1)) {
      lm <- lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST >= t,])
      message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 

### And from 11 to an upper threshold t 
for(t in seq(from = 12, to = 30, by = 1)) {
     lm <- lm(log(rich_phyto) ~ kT, data = ddf[ddf$SST < t & ddf$SST > 11,] )
     message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
     message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
     message(paste("", sep = ""))
} 

# Same with zoo
for(t in seq(from = 12, to = 30, by = 0.5)) {
      lm <- lm(log(rich_zoo) ~ kT, data = ddf[ddf$SST < t & ddf$SST > 11,] )
      message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 


### Niki & Meike want to see the heatmap of points density per bins of both kT and log(S) illustrated with a heatmap

lm3.phyto <- lm(log(rich_phyto) ~ kT + I(kT^2) + I(kT^3), data = ddf, na.action = na.exclude)
lm3.zoo <- lm(log(rich_zoo) ~ kT + I(kT^2) + I(kT^3), data = ddf, na.action = na.exclude)
ddf$phyto_3d <- predict(lm3.phyto)
ddf$zoo_3d <- predict(lm3.zoo)

ddf$kT_bin <- factor(cut_interval(ddf$kT,15))
levels(ddf$kT_bin)
levels <- str_replace_all(levels(ddf$kT_bin), ",", ":")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(ddf$kT_bin) <- levels

ddf$logphyto_bin <- factor(cut_interval(log(ddf$rich_phyto),10))
levels(ddf$logphyto_bin)
levels <- str_replace_all(levels(ddf$logphyto_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(ddf$logphyto_bin) <- levels

ddf$logzoo_bin <- factor(cut_interval(log(ddf$rich_zoo),15))
levels(ddf$logzoo_bin)
levels <- str_replace_all(levels(ddf$logzoo_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(ddf$logzoo_bin) <- levels

# Use dplyr to count n obs per binxbin 
require("viridis")
dens <- ddf %>% group_by(kT_bin,logphyto_bin) %>% summarise(n = n())
summary(data.frame(dens))
# Plot heatmap with geom_raster
heat1 <- ggplot() + geom_tile(aes(x = factor(kT_bin), y = factor(logphyto_bin), fill = log(n)), data = na.omit(dens)) + 
        scale_fill_viridis(name = "Density\n(log)", option = "A") + xlab("Thermal energy (1/kT)") + 
        ylab("Phytoplankton richness\n(log)") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#
ggsave(plot = heat1, filename = "heatmap_dens_richxkT_phyto.pdf", width = 9, height = 4, dpi = 300)

###
dens <- ddf %>% group_by(kT_bin,logzoo_bin) %>% summarise(n = n())
summary(data.frame(dens))
# Plot heatmap with geom_raster
heat2 <- ggplot() + geom_tile(aes(x = factor(kT_bin), y = factor(logzoo_bin), fill = log(n)), data = na.omit(dens)) + 
        scale_fill_viridis(name = "Density\n(log)", option = "A") + xlab("Thermal energy (1/kT)") + 
        ylab("Zooplankton richness\n(log)") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
ggsave(plot = heat2, filename = "heatmap_dens_richxkT_zoo.pdf", width = 9, height = 4, dpi = 300)




### UPDATE 24/09/19: need finer bins and overlay fit

lm3.phyto <- lm(log(rich_phyto) ~ kT + I(kT^2) + I(kT^3), data = ddf, na.action = na.exclude)
lm3.zoo <- lm(log(rich_zoo) ~ kT + I(kT^2) + I(kT^3), data = ddf, na.action = na.exclude)
ddf$phyto_3d <- predict(lm3.phyto)
ddf$zoo_3d <- predict(lm3.zoo)

ddf$kT_bin <- factor(cut_interval(ddf$kT,40))
levels(ddf$kT_bin)
levels <- str_replace_all(levels(ddf$kT_bin), ",", ":")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(ddf$kT_bin) <- levels

ddf$logphyto_bin <- factor(cut_interval(log(ddf$rich_phyto),20))
levels(ddf$logphyto_bin)
levels <- str_replace_all(levels(ddf$logphyto_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(ddf$logphyto_bin) <- levels

ddf$logzoo_bin <- factor(cut_interval(log(ddf$rich_zoo),20))
levels(ddf$logzoo_bin)
levels <- str_replace_all(levels(ddf$logzoo_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(ddf$logzoo_bin) <- levels

### Try adding a stat_density_2d(aes(fill = stat(level)), geom = "polygon")
# plot1 <-  <- ggplot() + stat_density_2d(aes(x = kT, y = phyto_3d), data = ddf, colour = "grey65") +
#             geom_line(aes(x = kT, y = phyto_3d), data = ddf, colour = "black") +
#             scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness\n(log)") +
#             xlab("Thermal energy (1/kT)") + theme_classic()
# #
# plot2 <- ggplot() + stat_density_2d(aes(x = kT, y = zoo_3d), data = ddf, colour = "grey65") +
#             geom_line(aes(x = kT, y = zoo_3d), data = ddf, colour = "black") +
#             scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness\n(log)") +
#             xlab("Thermal energy (1/kT)") + theme_classic()
#
ggsave(plot = plot1, filename = "plot_phyto_fit.pdf", dpi = 300, width = 4, height = 4)
ggsave(plot = plot2, filename = "plot_zoo_fit.pdf", dpi = 300, width = 4, height = 4)

#    stat_density_2d(aes(x=eruptions, y=waiting, color=geyser_types, fill=geyser_types), geom="polygon", bins=4, alpha=.1) +

### With polygon for density
plot1 <- ggplot() + stat_density_2d(data = na.omit(ddf), aes(x = kT, y = log(rich_phyto) ), 
                fill = "#66c2a5", fill = "#66c2a5", geom = "polygon", bins = 6, alpha = .3) +
            geom_line(aes(x = kT, y = phyto_3d), data = na.omit(ddf), colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness\n(log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot2 <- ggplot() + stat_density_2d(data = na.omit(ddf), aes(x = kT, y = log(rich_zoo) ), 
                fill = "#d53e4f", fill = "#d53e4f", geom = "polygon", bins = 6, alpha = .3) +
            geom_line(aes(x = kT, y = zoo_3d), data = na.omit(ddf), colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness\n(log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
ggsave(plot = plot1, filename = "plot_phyto_dens_polygon_fit.pdf", dpi = 300, width = 5, height = 4)
ggsave(plot = plot2, filename = "plot_zoo_dens_polygon_fit.pdf", dpi = 300, width = 5, height = 4)


### And by adding geom_point below stat_density_2d ?
plot1 <- ggplot() + geom_point(data = na.omit(ddf), aes(x = kT, y = log(rich_phyto) ), colour = "grey75", size = 1, alpha = .2) + 
            stat_density_2d(data = na.omit(ddf), aes(x = kT, y = log(rich_phyto) ), 
                colour = "#3288bd", fill = "#3288bd", geom = "polygon", bins = 6, alpha = .2) +
            geom_line(aes(x = kT, y = phyto_3d), data = na.omit(ddf), colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness (log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot2 <- ggplot() + geom_point(data = na.omit(ddf), aes(x = kT, y = log(rich_zoo) ), colour = "grey75", size = 1, alpha = .2) + 
            stat_density_2d(data = na.omit(ddf), aes(x = kT, y = log(rich_zoo) ), 
                colour = "#d53e4f", fill = "#d53e4f", geom = "polygon", bins = 6, alpha = .2) +
            geom_line(aes(x = kT, y = zoo_3d), data = na.omit(ddf), colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness (log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
ggsave(plot = plot1, filename = "plot_phyto_dens_points_polygon_fit.pdf", dpi = 300, width = 5, height = 5)
ggsave(plot = plot2, filename = "plot_zoo_dens_points_polygon_fit.pdf", dpi = 300, width = 5, height = 5)


### Values returned by stat_density_2d ?? # https://stackoverflow.com/questions/36157437/r-extract-coordinates-plotted-by-ggplot2stat-density-2d
# head( ggplot_build(plot1)$data[[2]] )
summary( ggplot_build(plot1)$data[[2]] ) # phyto
dim( ggplot_build(plot1)$data[[2]] ) 
summary( ggplot_build(plot2)$data[[2]] ) # zoo

str(ggplot_build(plot1)$data)
summary( ggplot_build(plot1)$data[[3]] )
head(ggplot_build(plot1)$data[[3]])

### NOTE: levels are computed using MASS::kde2d()
# ?kde2d 
# https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
# https://stackoverflow.com/questions/32206623/what-does-level-mean-in-ggplotstat-density2d

### And add a dahsed dline representing the expected slopes and intercepts
# To find the intercept = https://stackoverflow.com/questions/33292969/linear-regression-with-specified-slope
ddf2 <- na.omit(ddf)
slopeP <- 0.32
interP <- mean(log(ddf2$rich_phyto) - slopeP * ddf2$kT )
slopeZ <- 0.65
interZ <- mean(log(ddf2$rich_zoo) - slopeZ * ddf2$kT )
# Derive expected model
ddf2$phyto.mte <- slopeP*(ddf2$kT) + interP
ddf2$zoo.mte <- slopeZ*(ddf2$kT) + interZ
summary(ddf2)

### Compare R2 of MTE models to 3rd polynomial fit's
summary(lm(log(ddf2$rich_phyto) ~ ddf2$phyto.mte, data = ddf2)) # R-squared:  0.7843
summary(lm(log(ddf2$rich_zoo) ~ ddf2$zoo.mte, data = ddf2)) # R-squared:  0.853 
summary(lm3.phyto) # 0.8371 
summary(lm3.zoo) # 0.9705


### And by adding geom_point below stat_density_2d ?
plot1 <- ggplot() + geom_point(data = ddf2, aes(x = kT, y = log(rich_phyto) ), colour = "grey75", size = 1, alpha = .2) + 
            stat_density_2d(data = ddf2, aes(x = kT, y = log(rich_phyto) ), 
                colour = "#3288bd", fill = "#3288bd", geom = "polygon", bins = 6, alpha = .2) +
            geom_line(aes(x = kT, y = phyto.mte), data = ddf2, colour = "black", linetype = "dashed") +    
            geom_line(aes(x = kT, y = phyto_3d), data = ddf2, colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness (log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot2 <- ggplot() + geom_point(data = ddf2, aes(x = kT, y = log(rich_zoo) ), colour = "grey75", size = 1, alpha = .2) + 
            stat_density_2d(data = ddf2, aes(x = kT, y = log(rich_zoo) ), 
                colour = "#d53e4f", fill = "#d53e4f", geom = "polygon", bins = 6, alpha = .2) +
            geom_line(aes(x = kT, y = zoo.mte), data = ddf2, colour = "black", linetype = "dashed") +        
            geom_line(aes(x = kT, y = zoo_3d), data = ddf2, colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness (log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()

ggsave(plot = plot1, filename = "plot_phyto_dens_points_polygon_fit+expectedMTE.pdf", dpi = 300, width = 5, height = 5)
ggsave(plot = plot2, filename = "plot_zoo_dens_points_polygon_fit+expectedMTE.pdf", dpi = 300, width = 5, height = 5)


### ANd yet another possibility is to collide histograms next to the dotplot (using ggExtra)
require("ggExtra")

### add histgrams 
plot1 <- ggplot() + geom_point(data = ddf2, aes(x = kT, y = log(rich_phyto) ), colour = "grey75", size = 1, alpha = .2) + 
            stat_density_2d(data = ddf2, aes(x = kT, y = log(rich_phyto) ), 
                colour = "#3288bd", fill = "#3288bd", geom = "polygon", bins = 6, alpha = .2) +
            geom_line(aes(x = kT, y = phyto.mte), data = ddf2, colour = "black", linetype = "dashed") +    
            geom_line(aes(x = kT, y = phyto_3d), data = ddf2, colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Phytoplankton species richness (log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot2 <- ggplot() + geom_point(data = ddf2, aes(x = kT, y = log(rich_zoo) ), colour = "grey75", size = 1, alpha = .2) + 
            stat_density_2d(data = ddf2, aes(x = kT, y = log(rich_zoo) ), 
                colour = "#d53e4f", fill = "#d53e4f", geom = "polygon", bins = 6, alpha = .2) +
            geom_line(aes(x = kT, y = zoo.mte), data = ddf2, colour = "black", linetype = "dashed") +        
            geom_line(aes(x = kT, y = zoo_3d), data = ddf2, colour = "black") +
            scale_y_continuous(limits = c(3,5.5)) + ylab("Zooplankton species richness (log)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()

p1 <- ggExtra::ggMarginal(plot1, type = "histogram", data = ddf2, x = kT, y = log(rich_phyto), colour = "black", fill = "white")
p2 <- ggExtra::ggMarginal(plot2, type = "histogram",data = ddf2, x = kT, y = log(rich_zoo), colour = "black", fill = "white")

setwd(WD)
ggsave(plot = p1, filename = "plot_phyto_dens_points_polygon_fit+expectedMTE+hist.pdf", dpi = 300, width = 7, height = 7)
ggsave(plot = p2, filename = "plot_zoo_dens_points_polygon_fit+expectedMTE+hist.pdf", dpi = 300, width = 7, height = 7)

### But missing the scales on the density histograms...add them manually 
p1 <- ggplot(ddf2, aes(x=kT)) + geom_histogram(colour = "black", fill = "white") + 
        xlab("Thermal energy (1/kT)") + ylab("Density") + scale_x_continuous(limits = c(-42.77,-38.33)) + 
        theme_classic()
p2 <- ggplot(ddf2, aes(x=log(rich_phyto))) + geom_histogram(colour = "black", fill = "white") + 
        xlab("Phytoplankton species richness (log)") + ylab("Density") +  scale_x_continuous(limits = c(3,5.5)) + 
        theme_classic() + coord_flip()
p3 <- ggplot(ddf2, aes(x=log(rich_zoo))) + geom_histogram(colour = "black", fill = "white") + 
        xlab("Zooplankton species richness (log)") + ylab("Density") + scale_x_continuous(limits = c(3,5.5)) +
         theme_classic() + coord_flip()

ggsave(plot = p1, filename = "plot_dens_points_kT.pdf", dpi = 300, width = 7, height = 3)
ggsave(plot = p2, filename = "plot_dens_points_logphyto.pdf", dpi = 300, width = 3, height = 7)
ggsave(plot = p3, filename = "plot_dens_points_logzoo.pdf", dpi = 300, width = 3, height = 7)


# -----------------------------------------------------------

### 3°) Figure 3: Heatmap of rank correlation and/or linear slopes between the 5 diversity metrics and the env predictors

### A°) Spearman's rank correlations
library("corrplot")
colnames(ddf)
# Compute p-val matrix
p.mat <- cor.mtest(na.omit(ddf[,c(5:17,22,18,19,20,23:26)]))$p
# Correlation matrix
M <- cor(na.omit(ddf[,c(5:17,22,18,19,20,23:26)]), method = "spearman")
head(M)
# Provide names to p.mat's dimensions
rownames(p.mat) <- rownames(M)
colnames(p.mat) <- colnames(M)
head(p.mat)

# Need to melt M so you have 3 columns: env variables, diversity indices and corr values
m <- reshape2::melt(data = M, value.name = "cor", id.vars = colnames(ddf)[c(5:17,22)] )
# remove the line where Var2 %in% colnames(aclim)[6:17]
m <- m[!(m$Var2 %in% colnames(ddf)[c(5:17,22)]),]
m <- m[(m$Var1 %in% colnames(ddf)[c(5:17,22)]),]
colnames(m) <- c("env","div","corr")

# create a new variable from incidence
m$corrfactor <- cut(m$corr, 
		breaks = c(-1,-0.75,-0.5,-0.25,-0.10,-0.01,0.01,0.1,0.25,0.5,0.75,1),
		labels = c("-1.0","-0.75","-0.50","-0.25","-0.1","0","0.1","0.25","0.50","0.75","1.0") 
)
m$corrfactor <- factor(as.character(m$corrfactor),levels=rev(levels(m$corrfactor)))
# Change the levels of the factors according to what you want to display on the plot
unique(m$env)
unique(m$div)
# For re-naming
levels(m$env)[levels(m$env) == "logEKE"] <- "EKE"
levels(m$env)[levels(m$env) == "NO3"] <- "Nitrates"
levels(m$env)[levels(m$env) == "SiO2"] <- "Silicates"
levels(m$env)[levels(m$env) == "Nstar"] <- "N*"
levels(m$env)[levels(m$env) == "Sistar"] <- "Si*"
levels(m$env)[levels(m$env) == "Chla"] <- "Chlorophyll"
levels(m$env)[levels(m$env) == "kT"] <- "kT"
levels(m$div)[levels(m$div) == "rich_plankton"] <- "Plankton richness"
levels(m$div)[levels(m$div) == "rich_phyto"] <- "Phytoplankton richness"
levels(m$div)[levels(m$div) == "rich_zoo"] <- "Zooplankton richness"
levels(m$div)[levels(m$div) == "phyto_anom"] <- "Phytoplankton anomalies (log)"
levels(m$div)[levels(m$div) == "zoo_anom"] <- "Zooplankton anomalies (log)"
levels(m$div)[levels(m$div) == "phyto_anom2"] <- "Phytoplankton anomalies"
levels(m$div)[levels(m$div) == "zoo_anom2"] <- "Zooplankton anomalies"


### Do the same for p.mat
p <- melt(data = p.mat, value.name = "pval", id.vars = colnames(ddf)[c(5:17,22)] )
head(p)
# remove the line where Var2 %in% colnames(aclim)[6:17]
p <- p[!(p$Var2 %in% colnames(ddf)[c(5:17,22)]),]
p <- p[(p$Var1 %in% colnames(ddf)[c(5:17,22)]),]
colnames(p) <- c("env","div","pval")
# Change the levels of the factors according to what you want to display on the plot
unique(p$env)
unique(p$div)

### Change p$pval levels to "*" etc.
summary(p$pval)
p[p$pval > 0.01,]
# 318 dSST rich_phyto 0.0106
# 428   O2  zoo_anom2 0.1251

#define a colour for fonts
textcol <- "grey30"
# Diverging color palette
#levels(m$corrfactor) ; #unique(m$corrfactor) 11 levels but only 10 are realized (no 0)
#values <- c("#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
#quartz()
h1 <- ggplot(m,aes(x = env, y = div, fill = corr)) + geom_tile() + 
    geom_tile(colour = "white", size = 0.25, show_guide = F) +
    scale_fill_gradient2("Rank correlation", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-1,1)) +
	labs(x = "",y = "",title = "Rank correlations") + 
    scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
	geom_vline(aes(xintercept = 36), size = 3.4, alpha = 0.24) + coord_fixed() +
	theme_grey(base_size = 10) + geom_text(label = round(m$corr,3) ) +
	theme(	
		legend.title=element_blank(),
		#remove legend margin
		legend.margin = grid::unit(0,"cm"),
		#change legend text properties
		legend.text=element_text(colour=textcol,size=7,face="bold"),
		#change legend key height
		legend.key.height=grid::unit(0.8,"cm"),
		#set a slim legend
		legend.key.width=grid::unit(0.2,"cm"),
		#set x axis text size and colour
		axis.text.x=element_text(size=10,colour=textcol),
		#set y axis text colour and adjust vertical justification
		axis.text.y=element_text(vjust = 0.2,colour=textcol),
		#change axis ticks thickness
		axis.ticks=element_line(size=0.4),
		#change title font, size, colour and justification
		plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
		#remove plot background
		plot.background=element_blank(),
		#remove plot border
		panel.border=element_blank()
) 

setwd(WD)
ggsave(plot = h1, filename = "heat_Fig.3_rank_corrv2.pdf", dpi = 300, width = 12, height = 6)


### B°) Slopes
# https://stackoverflow.com/questions/51953709/fast-pairwise-simple-linear-regression-between-variables-in-a-data-frame
# pairwise_simple_lm <- function(dat) {
#
#       # matrix and its dimension (n: numbeta.ser of data; p: numbeta.ser of variables)
#       dat <- as.matrix(dat)
#       n <- nrow(dat)
#       p <- ncol(dat)
#       # variable summary: mean, (unscaled) covariance and (unscaled) variance
#       m <- colMeans(dat)
#       V <- crossprod(dat) - tcrossprod(m * sqrt(n))
#       d <- diag(V)
#       # R-squared (explained variance) and its complement
#       R2 <- (V ^ 2) * tcrossprod(1 / d)
#       R2_complement <- 1 - R2
#       R2_complement[seq.int(from = 1, by = p + 1, length = p)] <- 0
#       # slope and intercept
#       beta <- V * rep(1 / d, each = p) # slope
#       alpha <- m - beta * rep(m, each = p) # intercept
#       # residual sum of squares and standard error
#       RSS <- R2_complement * d
#       sig <- sqrt(RSS * (1 / (n - 2)))
#       # statistics for slope
#       beta.se <- sig * rep(1 / sqrt(d), each = p)
#       beta.tv <- beta / beta.se
#       beta.pv <- 2 * pt(abs(beta.tv), n - 2, lower.tail = FALSE)
#       # F-statistic and p-value
#       Ff.fv <- (n - 2) * R2 / R2_complement
#       Ff.pv <- pf(Ff.fv, 1, n - 2, lower.tail = FALSE)
#       # export
#       data.frame(LHS = rep(colnames(dat), times = p),
#              RHS = rep(colnames(dat), each = p),
#              alpha = c(alpha),
#              beta = c(beta),
#              beta.se = c(beta.se),
#              beta.tv = c(beta.tv),
#              beta.pv = c(beta.pv),
#              sig = c(sig),
#              R2 = c(R2),
#              Ff.fv = c(Ff.fv),
#              Ff.pv = c(Ff.pv),
#              stringsAsFactors = FALSE
#      ) # eo ddf
#
# } # eo FUN
#
# lm.res <- pairwise_simple_lm(dat = na.omit(ddf[,c(5:17,22,18,19,20,23,24)]))
# head(lm.res)
# summary(lm.res)
# str(lm.res)
# #
# lm.res <- lm.res[which(lm.res$LHS %in% colnames(ddf)[c(5:17,22)]),]
# lm.res <- lm.res[which(lm.res$RHS %in% c("rich_plankton","rich_phyto","rich_zoo","phyto_anom","zoo_anom")),]
# # OK
# summary(lm.res)
#
# ### Colour should be beta (slope, both positive and negative) and the label should be the R2
# # lm.res$beta
# # round(lm.res$R2,2)
# # For re-naming
# lm.res$RHS <- factor(lm.res$RHS)
# lm.res$LHS <- factor(lm.res$LHS)
# levels(lm.res$LHS)
# levels(lm.res$RHS)
# # Change levels
# levels(lm.res$LHS)[levels(lm.res$LHS) == "logEKE"] <- "EKE"
# levels(lm.res$LHS)[levels(lm.res$LHS) == "NO3"] <- "Nitrates"
# levels(lm.res$LHS)[levels(lm.res$LHS) == "SiO2"] <- "Silicates"
# levels(lm.res$LHS)[levels(lm.res$LHS) == "Nstar"] <- "N*"
# levels(lm.res$LHS)[levels(lm.res$LHS) == "Sistar"] <- "Si*"
# levels(lm.res$LHS)[levels(lm.res$LHS) == "Chla"] <- "Chlorophyll"
# levels(lm.res$LHS)[levels(lm.res$LHS) == "kT"] <- "kT"
# levels(lm.res$RHS)[levels(lm.res$RHS) == "rich_plankton"] <- "Plankton richness"
# levels(lm.res$RHS)[levels(lm.res$RHS) == "rich_phyto"] <- "Phytoplankton richness"
# levels(lm.res$RHS)[levels(lm.res$RHS) == "rich_zoo"] <- "Zooplankton richness"
# levels(lm.res$RHS)[levels(lm.res$RHS) == "phyto_anom"] <- "Phytoplankton anomalies"
# levels(lm.res$RHS)[levels(lm.res$RHS) == "zoo_anom"] <- "Zooplankton anomalies"
#
#
# h2 <- ggplot(lm.res, aes(x=LHS,y=RHS,fill=beta) ) + geom_tile() + geom_tile(colour="white",size=0.25, show_guide=FALSE) +
#     labs(x="",y="",title="Adjusted R2 and slope of linear models") +
#     scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
#     scale_fill_gradient2(name = "Slope", low = "#3288bd", high = "#d53e4f", mid = "white") +
#     geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24) + coord_fixed() +
#     theme_grey(base_size=10) + geom_text(label = round(lm.res$R2,3) ) +
#     theme(
#         legend.title=element_blank(),
#         legend.margin = grid::unit(0,"cm"),
#         legend.text=element_text(colour=textcol,size=7,face="bold"),
#         legend.key.height=grid::unit(0.8,"cm"),
#         legend.key.width=grid::unit(0.2,"cm"),
#         axis.text.x=element_text(size=10,colour=textcol),
#         axis.text.y=element_text(vjust = 0.2,colour=textcol),
#         axis.ticks=element_line(size=0.4),
#         plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
#         plot.background=element_blank(),
#         panel.border=element_blank()
# )
#
# setwd(WD)
# ggsave(plot = h2, filename = "heat_Fig.3_lm.pdf", dpi = 300, width = 12, height = 6)

### 27/08/19: Make bivariate plots to assess global PDR
plot <- ggplot() + geom_point(aes(x = Chla, y = rich_phyto, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Phytoplankton species richness") + 
     scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_richxchla_phyto_HSI.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = Chla, y = rich_zoo, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Zooplankton species richness") + 
      scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_richxchla_zoo_HSI.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = rich_phyto, y = rich_zoo, colour = Chla, size = abs(y)), data = ddf, alpha = 0.5) + 
         scale_colour_viridis(name = "Phytoplankton biomass\nlog(mgC.m3)") + 
         xlab("Phytoplankton species richness") + ylab("Zooplankton species richness") + 
         theme_classic()
ggsave(plot = plot, filename = "plot_phytoxzooxchla_rich_HSI.pdf", dpi = 300, width = 8, height = 6)

### 27/08/19: And with a fit
#summary(lm(rich_phyto ~ Chla, data = ddf)) # -15.60; 0.036
#summary(lm(rich_zoo ~ Chla, data = ddf)) # -90.38; 0.383
#summary(lm(phyto_anom ~ Chla, data = ddf)) # 0.287;  0.209
#summary(lm(zoo_anom ~ Chla, data = ddf)) # -0.079;  0.017

plot <- ggplot() + geom_point(aes(x = Chla, y = rich_phyto, colour = abs(y)), data = ddf, alpha = 0.6) + 
            geom_smooth(aes(x = Chla, y = rich_phyto, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
            xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Phytoplankton species richness") + 
            scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
            scale_y_continuous(limits = c(20,235)) + theme_classic()
            
ggsave(plot = plot, filename = "plot_richxchla_phyto_HSI_fit.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = Chla, y = rich_zoo, colour = abs(y)), data = ddf, alpha = 0.6) + 
             geom_smooth(aes(x = Chla, y = rich_zoo, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
             xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Zooplankton species richness") + 
             scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
             scale_y_continuous(limits = c(20,235)) + theme_classic()
             
ggsave(plot = plot, filename = "plot_richxchla_zoo_HSI_fit.pdf", dpi = 300, width = 6, height = 4)

plot <- ggplot() + geom_point(aes(x = Chla, y = phyto_anom, colour = abs(y)), data = ddf, alpha = 0.6) + 
            geom_smooth(aes(x = Chla, y = phyto_anom, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
            xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Phytoplankton richness anomalies") + 
            scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
            scale_y_continuous(limits = c(-0.54,1.35)) + geom_hline(yintercept = 0, linetype = "dashed") + 
            theme_classic()
            
ggsave(plot = plot, filename = "plot_anomsxchla_phyto_HSI_fit.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = Chla, y = zoo_anom, colour = abs(y)), data = ddf, alpha = 0.6) + 
             geom_smooth(aes(x = Chla, y = zoo_anom, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
             xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Zooplankton richness anomalies") + 
             scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
             scale_y_continuous(limits = c(-0.54,1.35)) + geom_hline(yintercept = 0, linetype = "dashed") + 
             theme_classic()
             
ggsave(plot = plot, filename = "plot_anomsxchla_zoo_HSI_fit.pdf", dpi = 300, width = 6, height = 4)


### 27/08/19: And with Nutrients availability? 
summary(lm(rich_phyto ~ NO3, data = ddf)) # -15.263; 0.382
summary(lm(rich_zoo ~ NO3, data = ddf)) # -38.137; 0.769
summary(lm(phyto_anom ~ NO3, data = ddf)) # 0.0137; 0.0052
summary(lm(zoo_anom ~ NO3, data = ddf)) # -0.0477; 0.0718

summary(lm(rich_phyto ~ SiO2, data = ddf)) # -12.9874; 0.198
summary(lm(rich_zoo ~ SiO2, data = ddf)) # -37.266; 0.526
summary(lm(phyto_anom ~ SiO2, data = ddf)) # 0.080227; 0.1296
summary(lm(zoo_anom ~ SiO2, data = ddf)) # -0.009277; 0.00191


plot <- ggplot() + geom_point(aes(x = NO3, y = rich_phyto, colour = abs(y)), data = ddf, alpha = 0.6) + 
            geom_smooth(aes(x = NO3, y = rich_phyto, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
            xlab("Nitrates concentration - log(µg/L)") + ylab("Phytoplankton species richness") + 
            scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
            scale_y_continuous(limits = c(20,235)) + theme_classic()
            
ggsave(plot = plot, filename = "plot_richxno3_phyto_HSI_fit.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = NO3, y = rich_zoo, colour = abs(y)), data = ddf, alpha = 0.6) + 
             geom_smooth(aes(x = NO3, y = rich_zoo, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
             xlab("Nitrates concentration - log(µg/L)") + ylab("Zooplankton species richness") + 
             scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
             scale_y_continuous(limits = c(20,235)) + theme_classic()
             
ggsave(plot = plot, filename = "plot_richxno3_zoo_HSI_fit.pdf", dpi = 300, width = 6, height = 4)

plot <- ggplot() + geom_point(aes(x = NO3, y = phyto_anom, colour = abs(y)), data = ddf, alpha = 0.6) + 
            geom_smooth(aes(x = NO3, y = phyto_anom, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
            xlab("Nitrates concentration - log(µg/L)") + ylab("Phytoplankton richness anomalies") + 
            scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
            scale_y_continuous(limits = c(-0.54,1.35)) + geom_hline(yintercept = 0, linetype = "dashed") + 
            theme_classic()
            
ggsave(plot = plot, filename = "plot_anomsxno3_phyto_HSI_fit.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = NO3, y = zoo_anom, colour = abs(y)), data = ddf, alpha = 0.6) + 
             geom_smooth(aes(x = NO3, y = zoo_anom, colour = abs(y)), data = ddf, method = "lm", colour = "black") + 
             xlab("Nitrates concentration - log(µg/L)") + ylab("Zooplankton richness anomalies") + 
             scale_colour_distiller(name = "Latitude", palette = "RdYlBu", direction = 1) + 
             scale_y_continuous(limits = c(-0.54,1.35)) + geom_hline(yintercept = 0, linetype = "dashed") + 
             theme_classic()
             
ggsave(plot = plot, filename = "plot_anomsxno3_zoo_HSI_fit.pdf", dpi = 300, width = 6, height = 4)



# Anomalies between S and SST
plot <- ggplot() + geom_point(aes(x = Chla, y = phyto_anom2, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Phytoplankton richness anomalies") + 
     scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_anoms2xchla_phyto.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = Chla, y = zoo_anom2, colour = abs(y)), data = ddf) + 
      xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Zooplankton richness anomalies") + 
      scale_colour_viridis(name = "Latitude", option = "A") + 
      theme_classic()
ggsave(plot = plot, filename = "plot_anoms2xchla_zoo.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = phyto_anom2, y = zoo_anom2, colour = Chla, size = abs(y)), data = ddf, alpha = 0.5) + 
         scale_colour_viridis(name = "Phytoplankton biomass\nlog(mgC.m3)") + 
         xlab("Phytoplankton richness anomalies") + ylab("Zooplankton richness anomalies") + 
         theme_classic()
ggsave(plot = plot, filename = "plot_phytoxzooxchla_anoms2.pdf", dpi = 300, width = 8, height = 6)



### 20/08/19: Find the best looking bivariate plots to hilight PDRs
# After exploring various possibilities, I think we need to do so per bioregions
# Perform PCA or CA
library("FactoMineR")
colnames(base)
pca <- PCA(X = base[,c(4:865)], graph = F, ncp = 5, scale.unit = F)
# Use cell's coords along CA1:CA4 to compute euclidean distances and cluster
summary(pca)
#                        Dim.1   Dim.2   Dim.3   Dim.4   Dim.5 
# Variance              16.711   3.254   1.973   1.404   1.064 
# % of var.             58.969  11.484   6.963   4.956   3.755
# Cumulative % of var.  58.969  70.452  77.415  82.371  86.126
# str(pca)
coords <- data.frame(pca$ind$coord)
colnames(coords) <- c("PC1","PC2","PC3","PC4","PC5")
# Compute euclidean distance using the CA scores
dist <- dist(coords[,c(1:5)], "euclidean")
# Perform HAC with Ward's algorithm
fit1 <- stats::hclust(dist, "ward.D2")
fit2 <- stats::hclust(dist, "average")
# Save dendrograms
#pdf(paste("dendro_ward_annual_plankton_comp_PCA_HSI_PC5.pdf", sep = ""), width = 15, height = 10)
#plot(fit1)
#dev.off()
#pdf(paste("dendro_UGPMA_annual_plankton_comp_PCA_HSI_PC5.pdf", sep = ""), width = 15, height = 10)
#plot(fit2)
#dev.off()

# Based on Ward's
base$k2_ward <- cutree(fit1, 2)
base$k3_ward <- cutree(fit1, 3)
base$k4_ward <- cutree(fit1, 4)
base$k5_ward <- cutree(fit1, 5)
base$k6_ward <- cutree(fit1, 6)
# Based on Ward's
base$k2_ugpma <- cutree(fit2, 2)
base$k3_ugpma <- cutree(fit2, 3)
base$k4_ugpma <- cutree(fit2, 4)
base$k5_ugpma <- cutree(fit2, 5)
base$k6_ugpma <- cutree(fit2, 6)

# Check their distrib
# colnames(base)
for(col in c(870:length(base))) {
    # Name
    name <- colnames(base)[col]
    meta <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(base[,col])), data = base) +
     	scale_fill_manual(name = "Metacommunity", values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
    setwd(WD)
    ggsave(plot = meta, filename = paste("map_bioregions_",name,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
    
} # eo for loop

### Save vectors of metacommunities !!!
# write.table(x = base[,c(1:3,866,870:879)], file = "metacommunities_global_plankton_annual_HSI_27_08_19.txt")

### Ok, stick to k3_ward (Tropical.Temperate/Polar) or k4_ugpma (has upwellings)
head(base[,c(1:3)])
head(ddf[,c(1:3)])
ddf$k4_ugpma <- factor(base[which(base$cell_id %in% unique(ddf$id)),c("k4_ugpma")])
ddf$k3_ward <- factor(base[which(base$cell_id %in% unique(ddf$id)),c("k3_ward")])
# Rename levels
levels(ddf$k4_ugpma)[levels(ddf$k4_ugpma) == "1"] <- "Tropics"
levels(ddf$k4_ugpma)[levels(ddf$k4_ugpma) == "2"] <- "Upwellings"
levels(ddf$k4_ugpma)[levels(ddf$k4_ugpma) == "3"] <- "Transition zones"
levels(ddf$k4_ugpma)[levels(ddf$k4_ugpma) == "4"] <- "High latitudes"

levels(ddf$k3_ward)[levels(ddf$k3_ward) == "1"] <- "Tropical"
levels(ddf$k3_ward)[levels(ddf$k3_ward) == "2"] <- "Temperate"
levels(ddf$k3_ward)[levels(ddf$k3_ward) == "3"] <- "Polar"

### WATCH OUT: you have some cells classified as "transition zones in the Artic" (y > 80) --> NaN
ddf[which(ddf$k4_ugpma == "Transition zones" & ddf$y > 80),"k4_ugpma"] <- NA
ddf[which(ddf$k3_ward == "Temperate" & ddf$y > 80),"k3_ward"] <- NA



### Plot PDR per metacommunity
# plot1 <- ggplot() + geom_point(aes(x = rich_phyto, y = rich_zoo, colour = Chla), data = ddf) +
#            scale_colour_viridis(name = "Phytoplankton biomass") +
#            xlab("Phytoplankton richness") + ylab("Zooplankton richness") +
#            theme_bw() + facet_wrap(~ factor(ddf$k4_ugpma), ncol = 2, nrow = 2)
# 
# plot2 <- ggplot() + geom_point(aes(x = rich_phyto, y = rich_zoo, colour = Chla), data = ddf) +
#            scale_colour_viridis(name = "Phytoplankton biomass") +
#            xlab("Phytoplankton richness") + ylab("Zooplankton richness") +
#            theme_bw() + facet_wrap(~ factor(ddf$k3_ward), ncol = 2, nrow = 2)
# 
# ggsave(plot = plot1, filename = "Fig.3B_ugpma_PDR_rich.pdf", dpi = 300, width = 9, height = 3)
# ggsave(plot = plot2, filename = "Fig.3B_ward_PDR_rich.pdf", dpi = 300, width = 9, height = 3)
#
#
# ### Plot PDR per metacommunity based on richness anomalies
# plot1 <- ggplot() + geom_point(aes(x = phyto_anom, y = zoo_anom, colour = Chla), data = ddf) +
#            scale_colour_viridis(name = "Phytoplankton biomass") +
#            xlab("Phytoplankton richness anomalies") + ylab("Zooplankton richness anomalies") +
#            theme_bw() + facet_wrap(~ factor(ddf$k4_ugpma), ncol = 2, nrow = 2)
# #
# plot2 <- ggplot() + geom_point(aes(x = phyto_anom, y = zoo_anom, colour = Chla), data = ddf) +
#            scale_colour_viridis(name = "Phytoplankton biomass") +
#            xlab("Phytoplankton richness anomalies") + ylab("Zooplankton richness anomalies") +
#            theme_bw() + facet_wrap(~ factor(ddf$k3_ward), ncol = 2, nrow = 2)
# 
# ggsave(plot = plot1, filename = "Fig.3B_ugpma_PDR_anoms.pdf", dpi = 300, width = 10, height = 5)
# ggsave(plot = plot2, filename = "Fig.3B_ward_PDR_anoms.pdf", dpi = 300, width = 10, height = 5)


### And facet per kingdom also
colnames(ddf)
ddf2 <- ddf[,c(2,3,5:17,19,20,27,28)]
colnames(ddf2)
m2 <- melt(data = ddf2, id.vars = colnames(ddf2)[c(1:15,18:length(ddf2))] )
colnames(m2)[18] <- c("div.var")
# Correct factor levels
levels(m2$div.var)[levels(m2$div.var) == "rich_phyto"] <- "Phytoplankton"
levels(m2$div.var)[levels(m2$div.var) == "rich_zoo"] <- "Zooplankton"

# Based on k4_ugpma 
p1 <- ggplot( m2[!is.na(m2$k4_ugpma),] ) + geom_point(aes(x = Chla, y = value, colour = factor(k4_ugpma)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(k4_ugpma) ) 

ggsave(plot = p1, filename = "Fig.3B_ugpma_PDR_rich_fit.pdf", dpi = 300, width = 10, height = 4)
# Based on k3_ward
p2 <- ggplot(m2[!is.na(m2$k3_ward),]) + geom_point(aes(x = Chla, y = value, colour = factor(k3_ward)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(k3_ward) ) 

ggsave(plot = p2, filename = "Fig.3B_ward_PDR_rich_fit.pdf", dpi = 300, width = 10, height = 4)
 

### Ok, nice, do the same with anomalies
colnames(ddf)
ddf2 <- ddf[,c(2,3,5:17,23,25,27,28)]
colnames(ddf2)
m2 <- melt(data = ddf2, id.vars = colnames(ddf2)[c(1:15,18:length(ddf2))] )
colnames(m2)[18] <- c("div.var")
# Correct factor levels
levels(m2$div.var)[levels(m2$div.var) == "phyto_anom"] <- "Phytoplankton"
levels(m2$div.var)[levels(m2$div.var) == "zoo_anom"] <- "Zooplankton"

# Based on k4_ugpma
p1 <- ggplot(m2[!is.na(m2$k4_ugpma),]) + geom_point(aes(x = Chla, y = value, colour = factor(k4_ugpma)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness anomalies") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(k4_ugpma) ) 
#
ggsave(plot = p1, filename = "Fig.3B_ugpma_PDR_anoms_fit.pdf", dpi = 300, width = 10, height = 4)
# Based on k3_ward
p2 <- ggplot( m2[!is.na(m2$k3_ward),] ) + geom_point(aes(x = Chla, y = value, colour = factor(k3_ward)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness anomalies") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(k3_ward) ) 
#
ggsave(plot = p2, filename = "Fig.3B_ward_PDR_anoms_fit.pdf", dpi = 300, width = 10, height = 4)

### And perform linear models - Adjusted R2 values (all p-values < 2.2e-16)
# k4_ugpma
summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k4_ugpma == "Tropics",])) # 0.4762
summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k4_ugpma == "Upwellings",])) # 0.6567
summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k4_ugpma == "Transition zones",])) # 0.14
summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k4_ugpma == "High latitudes",])) # 0.227

summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k4_ugpma == "Tropics",])) # 0.004; p-value: 8.917e-16
summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k4_ugpma == "Upwellings",])) # 0.4573
summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k4_ugpma == "Transition zones",])) # 0.27
summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k4_ugpma == "High latitudes",])) # n.s
# k3_ward
# summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k3_ward == "Tropical",])) # 0.6292
# summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k3_ward == "Temperate",])) #  0.0108
# summary(lm(rich_phyto ~ Chla, data = ddf[ddf$k3_ward == "Polar",])) # 0.215
#
# summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k3_ward == "Tropical",])) # 0.038
# summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k3_ward == "Temperate",])) # 0.445
# summary(lm(rich_zoo ~ Chla, data = ddf[ddf$k3_ward == "Polar",])) # 0.044

### Anomalies
# k4_ugpma
summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k4_ugpma == "Tropics",])) # 0.5483
summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k4_ugpma == "Upwellings",])) #  0.48
summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k4_ugpma == "Transition zones",])) # 0.448 
summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k4_ugpma == "High latitudes",])) # 0.297 

summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k4_ugpma == "Tropics",])) #  0.0588
summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k4_ugpma == "Upwellings",])) # 0.203
summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k4_ugpma == "Transition zones",])) # 0.0758
summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k4_ugpma == "High latitudes",])) # 0.05
# k3_ward
# summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k3_ward == "Tropical",])) # 0.6703
# summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k3_ward == "Temperate",])) # 0.1572
# summary(lm(phyto_anom ~ Chla, data = ddf[ddf$k3_ward == "Polar",])) # 0.251
#
# summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k3_ward == "Tropical",])) # 0.0019; p-value: 9.967e-07
# summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k3_ward == "Temperate",])) #  0.042
# summary(lm(zoo_anom ~ Chla, data = ddf[ddf$k3_ward == "Polar",])) #  0.081


### And zoo div ~ phyto div
summary(lm(rich_zoo ~ rich_phyto, data = ddf)) # Adjusted R-squared: 0.525; p-value: < 2.2e-16 # because both governed by SST
summary(lm(zoo_anom ~ phyto_anom, data = ddf)) # Adjusted R-squared: 0.0102; p-value: < 2.2e-16
### And anoms2
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf)) # Adjusted R-squared: 0.007; p-value: < 2.2e-16

#
p <- ggplot(ddf[!is.na(ddf$k4_ugpma),]) + geom_point(aes(x = rich_phyto, y = rich_zoo, colour = factor(k4_ugpma)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = rich_phyto, y = rich_zoo), colour = "black", method = "lm") + 
	xlab("Phytoplankton species richness") + ylab("Zooplankton species richness") + 
	theme_bw() + facet_grid(. ~ factor(k4_ugpma) )  
ggsave(plot = p, filename = "Fig.3C_ugpma_phytoxzoo_rich_fit.pdf", dpi = 300, width = 10, height = 2.5)
#
summary(lm(rich_zoo ~ rich_phyto, data = ddf[ddf$k4_ugpma == "Tropics",])) # -0.284687; 0.1065; p-value < 2.2e-16
summary(lm(rich_zoo ~ rich_phyto, data = ddf[ddf$k4_ugpma == "Upwellings",])) # -0.42675; 0.581; p-value < 2.2e-16
summary(lm(rich_zoo ~ rich_phyto, data = ddf[ddf$k4_ugpma == "Transition zones",])) # n.s.
summary(lm(rich_zoo ~ rich_phyto, data = ddf[ddf$k4_ugpma == "High latitudes",])) #  1.45623; 0.324; p-value < 2.2e-16
#
p <- ggplot(ddf[!is.na(ddf$k4_ugpma),]) + geom_point(aes(x = phyto_anom, y = zoo_anom, colour = factor(k4_ugpma)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = phyto_anom, y = zoo_anom), colour = "black", method = "lm") + 
	ylab("Zooplankton richness anomalies") + xlab("Phytoplankton richness anomalies") + 
	theme_bw() + facet_grid(. ~ factor(k4_ugpma) ) 
ggsave(plot = p, filename = "Fig.3C_ugpma_phytoxzoo_anoms_fit.pdf", dpi = 300, width = 10, height = 2.5)
#
summary(lm(zoo_anom ~ phyto_anom, data = ddf[ddf$k4_ugpma == "Tropics",])) # -0.446; 0.121; p-value < 2.2e-16
summary(lm(zoo_anom ~ phyto_anom, data = ddf[ddf$k4_ugpma == "Upwellings",])) # n.s. 
summary(lm(zoo_anom ~ phyto_anom, data = ddf[ddf$k4_ugpma == "Transition zones",])) # -0.324; 0.294; p-value < 2.2e-16
summary(lm(zoo_anom ~ phyto_anom, data = ddf[ddf$k4_ugpma == "High latitudes",])) # 0.5971; 0.354; p-value < 2.2e-16


### 27/08/19: And examine PDRs with anoms2
colnames(ddf)
ddf2 <- ddf[,c(2,3,5:17,24,26,27,28)]
colnames(ddf2)
m2 <- melt(data = ddf2, id.vars = colnames(ddf2)[c(1:15,18:length(ddf2))] )
colnames(m2)[18] <- c("div.var")
# Correct factor levels
levels(m2$div.var)[levels(m2$div.var) == "phyto_anom2"] <- "Phytoplankton"
levels(m2$div.var)[levels(m2$div.var) == "zoo_anom2"] <- "Zooplankton"


# Based on k4_ugpma
p <- ggplot(m2[!is.na(m2$k4_ugpma),]) + geom_point(aes(x = Chla, y = value, colour = factor(k4_ugpma)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness anomalies") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(k4_ugpma) ) 
#
ggsave(plot = p, filename = "Fig.3B_ugpma_PDR_anoms2_fit.pdf", dpi = 300, width = 10, height = 4)

# Based on k3_ward
p <- ggplot( m2[!is.na(m2$k3_ward),] ) + geom_point(aes(x = Chla, y = value, colour = factor(k3_ward)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness anomalies") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(k3_ward) ) 
#
ggsave(plot = p, filename = "Fig.3B_ward_PDR_anoms2_fit.pdf", dpi = 300, width = 10, height = 4)

# Cross-plot
p <- ggplot(ddf[!is.na(ddf$k4_ugpma),]) + geom_point(aes(x = phyto_anom2, y = zoo_anom2, colour = factor(k4_ugpma)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = phyto_anom2, y = zoo_anom2), colour = "black", method = "lm") + 
	ylab("Zooplankton richness anomalies") + xlab("Phytoplankton richness anomalies") + 
	theme_bw() + facet_grid(. ~ factor(k4_ugpma), scales = "free_x" ) 
ggsave(plot = p, filename = "Fig.3C_ugpma_phytoxzoo_anoms2_fit.pdf", dpi = 300, width = 10, height = 2.5)

p <- ggplot(ddf[!is.na(ddf$k3_ward),]) + geom_point(aes(x = phyto_anom2, y = zoo_anom2, colour = factor(k3_ward)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = phyto_anom2, y = zoo_anom2), colour = "black", method = "lm") + 
	ylab("Zooplankton richness anomalies") + xlab("Phytoplankton richness anomalies") + 
	theme_bw() + facet_grid(. ~ factor(k3_ward), scales = "free_x" ) 
ggsave(plot = p, filename = "Fig.3C_ward_phytoxzoo_anoms2_fit.pdf", dpi = 300, width = 10, height = 2.5)


### Summarize linear models
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k4_ugpma == "Tropics",])) # -0.653; 0.12
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k4_ugpma == "Upwellings",])) # -0.269; 0.137
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k4_ugpma == "Transition zones",])) # -0.539; 0.233
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k4_ugpma == "High latitudes",])) # 1.130; 0.514

summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k3_ward == "Tropical",])) # n.s
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k3_ward == "Temperate",])) # -0.391; 0.059
summary(lm(zoo_anom2 ~ phyto_anom2, data = ddf[ddf$k3_ward == "Polar",])) # 1.289; 0.583



# -----------------------------------------------------------

### 4°) Panel of climate change impacts: panel of contemporary annual diversity estimates (total plankton, phyto-, zooplankton) 
phyto.base <- read.table("table_phyto_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
phyto.fut <- read.table("table_phyto_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
zoo.base <- read.table("table_zoo_annual_composition_baseline_14_08.txt", h = T, sep = "\t")
zoo.fut <- read.table("table_zoo_annual_composition_2100-2071_rcp85_14_08.txt", h = T, sep = "\t")
base <- cbind(phyto.base,zoo.base[,c(4:length(zoo.base))])
fut <- cbind(phyto.fut,zoo.fut[,c(4:length(zoo.fut))])
base <- na.omit(base)
fut <- na.omit(fut)
base$x2 <- base$x 
base[base$x < 0 ,"x2"] <- (base[base$x < 0 ,"x"]) + 360
base$cell_id <- factor(paste(base$x2, base$y, sep = "_"))
fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
base <- base[order(base$cell_id),]
fut <- fut[order(fut$cell_id),]
#head(base$cell_id); head(fut$cell_id)
base <- base[base$cell_id %in% unique(fut$cell_id),]
fut <- fut[fut$cell_id %in% unique(base$cell_id),]
fut$x2 <- fut$x
# Add a factor specifying the time period
base$period <- factor("baseline")
fut$period <- factor("future")

### Compute baseline and future richness diff and associated % 
base$rich_plankton <- rowSums(as.matrix(base[,c(4:865)]))
base$rich_phyto <- rowSums(as.matrix(base[,c(4:341)]))
base$rich_zoo <- rowSums(as.matrix(base[,c(342:865)]))
fut$rich_plankton <- rowSums(as.matrix(fut[,c(4:865)]))
fut$rich_phyto <- rowSums(as.matrix(fut[,c(4:341)]))
fut$rich_zoo <- rowSums(as.matrix(fut[,c(342:865)]))
# Differences
base$diff_plankton <- (fut$rich_plankton) - (base$rich_plankton)
base$diff_phyto <- (fut$rich_phyto) - (base$rich_phyto)
base$diff_zoo <- (fut$rich_zoo) - (base$rich_zoo)
# % diff
base$perc_plankton <- ((base$diff_plankton)/(base$rich_plankton))*100
base$perc_phyto <- ((base$diff_phyto)/(base$rich_phyto))*100
base$perc_zoo <- ((base$diff_zoo)/(base$rich_zoo))*100
# Check
summary(base[,c(868:length(base))])

### Maps
map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc_plankton), data = na.omit(base[base$perc_plankton < 50,])) +
    geom_raster(aes(x = x2, y = y), data = na.omit(base[base$perc_plankton >= 50,]), fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x2, y = y, z = perc_plankton), 
            data = na.omit(base[base$perc_plankton < 50,]) ) +
 	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-26,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc_phyto), data = na.omit(base[base$perc_phyto < 50,])) +
    geom_raster(aes(x = x2, y = y), data = na.omit(base[base$perc_phyto >= 50,]), fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x2, y = y, z = perc_phyto), 
            data = na.omit(base[base$perc_phyto < 50,]) ) +
 	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-45,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc_zoo), data = na.omit(base[base$perc_zoo < 50,])) +
    geom_raster(aes(x = x2, y = y), data = na.omit(base[base$perc_zoo >= 50,]), fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x2, y = y, z = perc_zoo), 
            data = na.omit(base[base$perc_zoo < 50,]) ) +
 	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-26,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )	
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_Fig.4A_perc_plankton_2100-2000_rcp85_v3.pdf", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_Fig.4C_perc_phyto_2100-2000_rcp85_v3.pdf", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_Fig.4E_perc_zoo_2100-2000_rcp85_v3.pdf", dpi = 300, width = 7, height = 5)


# ### For testing other diverging gradients
# map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc_zoo), data = na.omit(base[base$perc_zoo < 50,])) +
#      geom_raster(aes(x = x2, y = y), data = na.omit(base[base$perc_zoo >= 50,]), fill = "#b2182b") +
#      geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x2, y = y, z = perc_zoo),
#              data = na.omit(base[base$perc_zoo < 50,]) ) +
#      scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-26,50)) +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#      coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#      scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#      theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
# ggsave(plot = map, filename = "map_Fig.4_test2.pdf", dpi = 300, width = 7, height = 5)


### Make zonal plots too
library("dplyr")
zonal <- data.frame(base %>% 
            group_by(y) %>% 
            summarise(mean_plankton = mean(perc_plankton,na.rm=T), sd_plankton = sd(perc_plankton,na.rm=T), 
            mean_phyto = mean(perc_phyto,na.rm=T), sd_phyto = sd(perc_phyto,na.rm=T),
            mean_zoo = mean(perc_zoo,na.rm=T), sd_zoo = sd(perc_zoo,na.rm=T) )    
) # eo ddf
summary(zonal)

### Plots
plot1 <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed") + 
            geom_ribbon(aes(x = y, ymin = mean_plankton - sd_plankton, ymax = mean_plankton + sd_plankton), 
            fill = "grey70", data = zonal, alpha = 0.66) +
            geom_line(aes(x = y, y = mean_plankton), data = zonal, colour = "black" ) + scale_y_continuous(limits = c(-25,65)) + 
		    ylab("Difference in species richness\n(Plankton)") + xlab("Latitude (°)") + 
		    theme_classic() + coord_flip()
#
plot2 <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed") + 
            geom_ribbon(aes(x = y, ymin = mean_phyto - sd_phyto, ymax = mean_phyto + sd_phyto), 
            fill = "grey70", data = zonal, alpha = 0.66) +
            geom_line(aes(x = y, y = mean_phyto), data = zonal, colour = "black" ) + scale_y_continuous(limits = c(-25,66)) + 
		    ylab("Difference in species richness\n(Phytoplankton)") + xlab("Latitude (°)") + 
		    theme_classic() + coord_flip()
#
plot3 <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed") + 
            geom_ribbon(aes(x = y, ymin = mean_zoo - sd_zoo, ymax = mean_zoo + sd_zoo),
            fill = "grey70", data = zonal, alpha = 0.66) +
            geom_line(aes(x = y, y = mean_zoo), data = zonal, colour = "black" ) + scale_y_continuous(limits = c(-25,65)) + 
		    ylab("Difference in species richness\n(Zooplankton)") + xlab("Latitude (°)") + 
		    theme_classic() + coord_flip()
#
ggsave(plot = plot1, filename = "plot_Fig.4B_plankton_rich_baseline_binom.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot2, filename = "plot_Fig.4D_phyto_rich_baseline_binom.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot3, filename = "plot_Fig.4F_zoo_rich_baseline_binom.pdf", dpi = 300, width = 3, height = 4)


### 28/08/19: Make a heatmap or rank correlations between delta env and % changes in richness to investigate drivers of change
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/GFDL-ESM2SM")
# dir()
vars <- c("dSST","SST","dO2","logNO3","logSiO2","Nstar","Sistar","logChl")
### Load all and compute annual climatology along the way, except for dSST
# v <- "SST"
# v <- "dSST"
clims <- lapply(vars, function(v) {            
            message(paste("Loading ", v, sep = ""))           
            d <- read.table(paste("clims_deltas_diff_2100-2000_",v,"_rcp85_GFDL-ESM2M.txt", sep = ""), sep = "\t", h = T)
            if(v == "dSST") {
                colnames(d)[4] <- v
                return(d)
            } else {
                d$diff <- rowMeans(as.matrix(d[,c(4:length(d))]))
                colnames(d)[length(d)] <- v
                # summary(d)
                return( data.frame(v = d[,length(d)]) )
            } # eo if else loop   
    } # eo FUN
) # eo FUN
# Rbind
clim <- do.call(cbind,clims)
dim(clim); summary(clim)
colnames(clim)[c(5:11)] <- c("SST","dO2","logNO3","logSiO2","Nstar","Sistar","logChl")
rm(clims);gc()

dim(clim)
dim(base)

### Combine div changes of plankton/phyto/zoo with those in a ddf
ddf <- cbind( base[,c("cell_id","x2","y","perc_plankton","perc_phyto","perc_zoo")], 
            clim[clim$id %in% base$cell_id,c("SST","dSST","dO2","logNO3","logSiO2","Nstar","Sistar","logChl")]
) # eo ddf 
summary(ddf)

### A°) Spearman's rank correlations
library("corrplot")
colnames(ddf)
# Compute p-val matrix
p.mat <- cor.mtest(na.omit(ddf[,c(4:14)]))$p
# Correlation matrix
M <- cor(na.omit(ddf[,c(4:14)]), method = "spearman")
head(M)
# Provide names to p.mat's dimensions
rownames(p.mat) <- rownames(M)
colnames(p.mat) <- colnames(M)
head(p.mat)

# Need to melt M so you have 3 columns: env variables, diversity indices and corr values
m <- reshape2::melt(data = M, value.name = "cor", id.vars = colnames(ddf)[c(7:14)] )
# remove the line where Var2 %in% colnames(aclim)[6:17]
m <- m[!(m$Var2 %in% colnames(ddf)[c(7:14)]),]
m <- m[(m$Var1 %in% colnames(ddf)[c(7:14)]),]
colnames(m) <- c("env","div","corr")

# create a new variable from incidence
m$corrfactor <- cut(m$corr, 
		breaks = c(-1,-0.75,-0.5,-0.25,-0.10,-0.01,0.01,0.1,0.25,0.5,0.75,1),
		labels = c("-1.0","-0.75","-0.50","-0.25","-0.1","0","0.1","0.25","0.50","0.75","1.0") 
)
m$corrfactor <- factor(as.character(m$corrfactor),levels=rev(levels(m$corrfactor)))
# Change the levels of the factors according to what you want to display on the plot
unique(m$env)
unique(m$div)
# For re-naming
levels(m$env)[levels(m$env) == "logNO3"] <- "Nitrates"
levels(m$env)[levels(m$env) == "logSiO2"] <- "Silicates"
levels(m$env)[levels(m$env) == "Nstar"] <- "N*"
levels(m$env)[levels(m$env) == "Sistar"] <- "Si*"
levels(m$env)[levels(m$env) == "logChl"] <- "Chlorophyll"
levels(m$div)[levels(m$div) == "perc_plankton"] <- "Plankton"
levels(m$div)[levels(m$div) == "perc_phyto"] <- "Phytoplankton"
levels(m$div)[levels(m$div) == "perc_zoo"] <- "Zooplankton"


### Do the same for p.mat
p <- melt(data = p.mat, value.name = "pval", id.vars = colnames(ddf)[c(7:14)] )
head(p)
# remove the line where Var2 %in% colnames(aclim)[6:17]
p <- p[!(p$Var2 %in% colnames(ddf)[c(7:14)]),]
p <- p[(p$Var1 %in% colnames(ddf)[c(7:14)]),]
colnames(p) <- c("env","div","pval")
### Change p$pval levels to "*" etc.
summary(p$pval)
p[p$pval > 0.01,]  #0

#define a colour for fonts
textcol <- "grey30"

h1 <- ggplot(m,aes(x = env, y = div, fill = corr)) + geom_tile() + 
    geom_tile(colour = "white", size = 0.25, show_guide = F) +
    scale_fill_gradient2("Rank correlation", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-0.8,0.8)) +
	labs(x = "",y = "",title = "Rank correlations") + 
    scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
	geom_vline(aes(xintercept = 36), size = 3.4, alpha = 0.24) + coord_fixed() +
	theme_grey(base_size = 10) + geom_text(label = round(m$corr,3) ) +
	theme(	
		legend.title=element_blank(),
		#remove legend margin
		legend.margin = grid::unit(0,"cm"),
		#change legend text properties
		legend.text=element_text(colour=textcol,size=7,face="bold"),
		#change legend key height
		legend.key.height=grid::unit(0.8,"cm"),
		#set a slim legend
		legend.key.width=grid::unit(0.2,"cm"),
		#set x axis text size and colour
		axis.text.x=element_text(size=10,colour=textcol),
		#set y axis text colour and adjust vertical justification
		axis.text.y=element_text(vjust = 0.2,colour=textcol),
		#change axis ticks thickness
		axis.ticks=element_line(size=0.4),
		#change title font, size, colour and justification
		plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
		#remove plot background
		plot.background=element_blank(),
		#remove plot border
		panel.border=element_blank()
) 

setwd(WD)
ggsave(plot = h1, filename = "heat_Fig.4_rank_corrv.pdf", dpi = 300, width = 8, height = 4)

### 30/08/19: To help you interpret this heatmap: map the deltas of these vars and plot the bivariate relationships with 
### phyto- and zooplankton diversity
colnames(ddf)[c(4:14)] <- c("Plankton","Phytoplankton","Zooplankton","SST","dSST","O2","NO3","SiOH4","N*","Si*","Chl")
vars <- c("SST","dSST","O2","NO3","SiOH4","N*","Si*","Chl")
### For each v in vars, map deltas and plot relationships with div differences (%)
for(v in vars) {
    
        # Make map
        if(v == "Si*") {
            binw <- 3  
        } else {
            binw <- 1   
        }
        map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = ddf[,c(v)]), data = ddf ) +
            geom_contour(colour = "grey60", binwidth = binw, size = 0.25, aes(x = x2, y = y, z = ddf[,c(v)]), data = ddf) +
 	        scale_fill_gradient2(name = v, low = "#3288bd", high = "#d53e4f", mid = "white") +
 	        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	        coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	        scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		          labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        # Save map
        ggsave(plot = map, filename = paste("map_delta_2100-2000_rcp85_",v,".pdf", sep = ""), dpi = 300, width = 6, height = 4)
        
        # Make plots
        plot1 <- ggplot() + geom_point(aes(x = ddf[,c(v)], y = Phytoplankton, fill = abs(y)), data = ddf, pch = 21, colour = "black", alpha = 0.5) + 
                scale_fill_distiller(palette = "RdYlBu", name = "Latitude") + 
                geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
                xlab(v) + ylab("Richness difference\n(Phytoplankton)") + theme_bw()
        
        plot2 <- ggplot() + geom_point(aes(x = ddf[,c(v)], y = Zooplankton, fill = abs(y)), data = ddf, pch = 21, colour = "black", alpha = 0.5) + 
                scale_fill_distiller(palette = "RdYlBu", name = "Latitude") + 
                geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
                xlab(v) + ylab("Richness difference\n(Zooplankton)") + theme_bw()
                
        # Save plots
        ggsave(plot = plot1, filename = paste("plot_",v,"xdiff_phyto_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)
        ggsave(plot = plot2, filename = paste("plot_",v,"xdiff_zoo_2100-2000_rcp85.pdf", sep = ""), dpi = 300, width = 6, height = 4)
    
} # eo v in vars

### And examine lm between richness changes and vars changes
# SST
summary(lm(Phytoplankton ~ SST, data = ddf)) # R-squared: 0.4376; 8.098
summary(lm(Zooplankton ~ SST, data = ddf)) # R-squared: 0.009259; -0.946
# dSST
summary(lm(Phytoplankton ~ dSST, data = ddf)) #  0.012; -1.64
summary(lm(Zooplankton ~ dSST, data = ddf)) #  0.006; -0.950
# dO2
summary(lm(Phytoplankton ~ O2, data = ddf)) # - don't care really
summary(lm(Zooplankton ~ O2, data = ddf)) # 0.061; -8.289
# Nitrates
summary(lm(Phytoplankton ~ NO3, data = ddf)) # 0.062; 15.92
summary(lm(Zooplankton ~ NO3, data = ddf)) # 0.029; -8.810
# Silicates
summary(lm(Phytoplankton ~ SiOH4, data = ddf)) # 0.095; 15.23
summary(lm(Zooplankton ~ SiOH4, data = ddf)) # 0.1988; -17.630 
# N*
summary(lm(Phytoplankton ~ ddf[,c("N*")], data = ddf)) # 0.0523; 6.647
summary(lm(Zooplankton ~ ddf[,c("N*")], data = ddf)) # 0.0067; -1.904
# Si*
summary(lm(Phytoplankton ~ ddf[,c("Si*")], data = ddf)) #  0.0021; 0.226
summary(lm(Zooplankton ~ ddf[,c("Si*")], data = ddf)) # 0.034; -0.714
# Chlorophyll
summary(lm(Phytoplankton ~ Chl, data = ddf)) # 0.0053; -12.24
summary(lm(Zooplankton ~ Chl, data = ddf)) # 0.094; 41.465

### From linear models and spearman's correlation coefficients, it is hard to see the effect of STS on zooplankton richness changes.
### That is likely due to the fact that that the initial relationship is NON UNIMODAL but rather trimodal (Fig. 2):
### decreases from lowest SST to low SST, then increases linearly until 23°C. So, could decrease in tropics because future SST > 23°C
### but still increase at high latitudes vecause of northwards shift due to species tracking their niche optima

# Split the data in 3 domains (Trop/Temp/Polar) and examine zooplankton div changes ~ SST in them
ddf$domain <- NA
ddf[abs(ddf$y) < 30,c("domain")] <- "Tropical"
ddf[abs(ddf$y) >= 30 & abs(ddf$y) <= 60,c("domain")] <- "Temperate"
ddf[abs(ddf$y) > 60,c("domain")] <- "Polar"
levels(factor(ddf$domain)); summary(factor(ddf$domain))

plot2 <- ggplot() + geom_point(aes(x = SiOH4, y = Zooplankton, fill = factor(domain)), data = ddf, pch = 21, colour = "black", alpha = 0.5) + 
        scale_fill_manual(name = "Latitude", values = c("#3288bd","#abd9e9","#d53e4f")) + 
        geom_smooth(aes(x = SiOH4, y = Zooplankton), data = ddf, colour = "black", method = "lm") + 
        #geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
        xlab("SiOH4") + ylab("Richness difference\n(Zooplankton)") + theme_classic() + 
        facet_wrap(.~factor(ddf$domain), scales = "free")
#
ggsave(plot = plot2, filename = paste("plot_","SiOH4","xdiff_zoo_2100-2000_rcp85_facet.pdf", sep = ""), dpi = 300, width = 8, height = 3)


 



