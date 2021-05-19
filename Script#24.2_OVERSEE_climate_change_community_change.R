
##### 22/08/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	    - Examining changes in alpha and beta diversity based on binary data (thresholded HSI), for zoo/phyto/total
#       - Examine metadommunity distribution and pairwise rank correlations
#       - Examine climate change impacts on alpha and beta.diversity
 
### Last update: 22/08/2019

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("dplyr")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")
library("viridis")
library("scales")
library("maps")
library("betapart")
library("cmocean")

world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------


### 1°) Set the working directories, vectors etc. to retrieve the species' mean thresholds
WD <- getwd()
setwd( paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future_rcp85/", sep = "") )
zoo.wd <- getwd()
setwd( paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future_rcp85/", sep = "") )
phyto.wd <- getwd()
setwd(WD)

# Vector of SDMs
SDMs <- c('GLM','GAM','RF','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
# p <- "p4"

### 2°) Plot distrbution of models'TSS values for zoo and phyto separately, using boxplots per SDM and facet per pool
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(zoo.wd,"/eval_scores_",p,"/", sep = "") )
					zoo_spp <- str_replace_all(dir(), "eval_scores_", "")
					zoo_spp <- str_replace_all(zoo_spp, ".Rdata", "")
					scores <- lapply(zoo_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- dplyr::bind_rows(scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}
) # eo lapply
# Rbind
table_scores_zoo <- do.call(rbind,res)
table_scores_zoo$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_zoo)), pattern = "_", n = 2)))[,1])
table_scores_zoo$kingdom <- "Zooplankton"
#dim(table_scores_zoo); head(table_scores_zoo)
#summary(table_scores_zoo)
rm(res)
### And for phyto
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(phyto.wd,"/eval_scores_",p,"/", sep = "") )
					phyto_spp <- str_replace_all(dir(), "eval_scores_", "")
					phyto_spp <- str_replace_all(phyto_spp, ".Rdata", "")
					scores <- lapply(phyto_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind, scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}
) # eo lapply
# Rbind
table_scores_phyto <- do.call(rbind,res)
table_scores_phyto$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_phyto)), pattern = "_", n = 2)))[,1])
table_scores_phyto$kingdom <- "Phytoplankton"
#dim(table_scores_phyto); head(table_scores_phyto)
#summary(table_scores_phyto)
rm(res)

# Compute the species' mean/median probability thresholds
require("dplyr")
zoo.thresh <- data.frame(na.omit(table_scores_zoo) %>%
  		group_by(species) %>%
  		summarise(avg_TSS = mean(TSS), cutoff = median(Cutoff_TSS) )
) # eo ddf
zoo.thresh[scZ$avg_TSS <= 0.3,] # none, only probel with ONE RUNE of S.magnus
zoo_spp <- unique(zoo.thresh$species)
zoo_spp <- zoo_spp[!(zoo_spp %in% c("Spinocalanus_magnus"))]

# For phytoplankton
phyto.thresh <- data.frame(table_scores_phyto %>%
  		group_by(species) %>%
  		summarise(avg_TSS = mean(TSS), cutoff = median(Cutoff_TSS) )
) # eo ddf
phyto_spp <- unique(phyto.thresh$species)
phyto_spp <- phyto_spp[!(phyto_spp %in% c("Actiniscus_pentasterias"))]

summary(phyto.thresh)
summary(zoo.thresh)

setwd(WD)


# --------------------------------------------------------------------------------------------------------------------------------


### 2°) Load baseline and future community tables
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
head(base$cell_id); head(fut$cell_id)
base <- base[base$cell_id %in% unique(fut$cell_id),]
fut <- fut[fut$cell_id %in% unique(base$cell_id),]
fut$x2 <- fut$x
# Add a factor specifying the time period
base$period <- factor("baseline")
fut$period <- factor("future")
#length(setdiff(base$cell_id, fut$cell_id)) # 0
#length(setdiff(fut$cell_id, base$cell_id)) # 0

### Use the species' mean thresholds above to convert probabilities to binary distribution
# Make species names match
colnames(base)[c(4:865)] <- gsub("[.]","",as.character(colnames(base)[c(4:865)])) 
colnames(fut)[c(4:865)] <- gsub("[.]","",as.character(colnames(fut)[c(4:865)])) 
zoo_spp <- gsub("\\(|\\)", "", as.character(zoo_spp))
zoo.thresh$species <- gsub("\\(|\\)", "", as.character(zoo.thresh$species))
#length(setdiff(zoo_spp, colnames(base)[c(342:865)])) # 0
#length(setdiff(zoo_spp, colnames(fut)[c(342:865)])) # 0
# Convert probabilities to binary distribution
fut2 <- fut
base2 <- base

# Converting phytoplankton probabilities
# sp <- "Ceratium_breve"
for(sp in colnames(base)[c(4:341)] ) {
        message(paste("Converting probabilities for ",sp, sep = ""))
        # Get species' cutoff
        cut <- phyto.thresh[phyto.thresh$species == sp,"cutoff"]/1000
        base2[which(base2[,c(sp)] >= cut),c(sp)] <- 1
        base2[which(base2[,c(sp)] < cut),c(sp)] <- 0
        fut2[which(fut2[,c(sp)] >= cut),c(sp)] <- 1
        fut2[which(fut2[,c(sp)] < cut),c(sp)] <- 0
} # eo for loop

# Converting zooplankton probabilities
# sp <- "Vogtia_glabra"
for(sp in colnames(base)[c(342:865)] ) {
        message(paste("Converting probabilities for ",sp, sep = ""))
         # Get species' threshold
         cut <- zoo.thresh[zoo.thresh$species == sp,"cutoff"]/1000
         base2[which(base2[,c(sp)] >= cut),c(sp)] <- 1
         base2[which(base2[,c(sp)] < cut),c(sp)] <- 0
         fut2[which(fut2[,c(sp)] >= cut),c(sp)] <- 1
         fut2[which(fut2[,c(sp)] < cut),c(sp)] <- 0   
} # eo for loop

# Check
summary(base2)
summary(fut2)

### Compute baseline, and future, species richness (total, phyto- and zooplankton)
base2$rich_plankton <- rowSums(as.matrix(base2[,c(4:865)]))
base2$rich_phyto <- rowSums(as.matrix(base2[,c(4:341)]))
base2$rich_zoo <- rowSums(as.matrix(base2[,c(342:865)]))
fut2$rich_plankton <- rowSums(as.matrix(fut2[,c(4:865)]))
fut2$rich_phyto <- rowSums(as.matrix(fut2[,c(4:341)]))
fut2$rich_zoo <- rowSums(as.matrix(fut2[,c(342:865)]))
#
summary(base2[,c(866:length(base2))])
summary(fut2[,c(866:length(fut2))])

### Map richness patterns for baseline
base2$rich_plankton_bin <- factor(cut_interval(base2$rich_plankton, 10))
levels(base2$rich_plankton_bin)
levels <- str_replace_all(levels(base2$rich_plankton_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(base2$rich_plankton_bin) <- levels

base2$rich_phyto_bin <- factor(cut_interval(base2$rich_phyto, 10))
levels(base2$rich_phyto_bin)
levels <- str_replace_all(levels(base2$rich_phyto_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(base2$rich_phyto_bin) <- levels

base2$rich_zoo_bin <- factor(cut_interval(base2$rich_zoo, 10))
levels(base2$rich_zoo_bin)
levels <- str_replace_all(levels(base2$rich_zoo_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(base2$rich_zoo_bin) <- levels

map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_plankton_bin)), data = na.omit(base2)) +
 	scale_fill_viridis(name = "", discrete = T ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_phyto_bin)), data = na.omit(base2)) +
 	scale_fill_viridis(name = "", discrete = T ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_zoo_bin)), data = na.omit(base2)) +
 	scale_fill_viridis(name = "", discrete = T ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )		
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_rich_binom_plankton_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_rich_binom_phyto_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_rich_binom_zoo_baseline.jpg", dpi = 300, width = 7, height = 5)

### Map ln(Richness) too
map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = log1p(rich_plankton)), data = na.omit(base2)) +
 	scale_fill_viridis(name = "", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = log1p(rich_phyto)), data = na.omit(base2)) +
 	scale_fill_viridis(name = "", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = log1p(rich_zoo)), data = na.omit(base2)) +
 	scale_fill_viridis(name = "", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )		
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_logrich_binom_plankton_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_logrich_binom_phyto_baseline.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_logrich_binom_zoo_baseline.jpg", dpi = 300, width = 7, height = 5)


### Make zonal plots too
library("dplyr")
zonal <- data.frame(base2 %>% 
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
ggsave(plot = plot1, filename = "plot_zonal_plankton_rich_baseline_binom.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot2, filename = "plot_zonal_phyto_rich_baseline_binom.pdf", dpi = 300, width = 3, height = 4)
ggsave(plot = plot3, filename = "plot_zonal_zoo_rich_baseline_binom.pdf", dpi = 300, width = 3, height = 4)



### Same, but for future
fut2$rich_plankton_bin <- factor(cut_interval(fut2$rich_plankton, 10))
levels(fut2$rich_plankton_bin)
levels <- str_replace_all(levels(fut2$rich_plankton_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(fut2$rich_plankton_bin) <- levels

fut2$rich_phyto_bin <- factor(cut_interval(fut2$rich_phyto, 10))
levels(fut2$rich_phyto_bin)
levels <- str_replace_all(levels(fut2$rich_phyto_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(fut2$rich_phyto_bin) <- levels

fut2$rich_zoo_bin <- factor(cut_interval(fut2$rich_zoo, 10))
levels(fut2$rich_zoo_bin)
levels <- str_replace_all(levels(fut2$rich_zoo_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels(fut2$rich_zoo_bin) <- levels

map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_plankton_bin)), data = na.omit(fut2)) +
 	scale_fill_viridis(name = "", discrete = T ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_phyto_bin)), data = na.omit(fut2)) +
 	scale_fill_viridis(name = "", discrete = T ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(rich_zoo_bin)), data = na.omit(fut2)) +
 	scale_fill_viridis(name = "", discrete = T ) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )		
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_rich_binom_plankton_2100_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_rich_binom_phyto_2100_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_rich_binom_zoo_2100_rcp85.jpg", dpi = 300, width = 7, height = 5)

### Map ln(Richness) too
map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = log1p(rich_plankton)), data = na.omit(fut2)) +
 	scale_fill_viridis(name = "", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = log1p(rich_phyto)), data = na.omit(fut2)) +
 	scale_fill_viridis(name = "", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = log1p(rich_zoo)), data = na.omit(fut2)) +
 	scale_fill_viridis(name = "", discrete = F) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )		
# Save
setwd(WD)
ggsave(plot = map1, filename = "map_logrich_binom_plankton_2100_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_logrich_binom_phyto_2100_rcp85.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_logrich_binom_zoo_2100_rcp85.jpg", dpi = 300, width = 7, height = 5)


# --------------------------------------------------------------------------------------------------------------------------------


### 3°) Load contemporary annual clims and combine with diversity estimates
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

### Combine with 'base'
base2$cell_id <- factor(paste(base2$x2, base$y, sep = "_"))
base2 <- base2[order(base2$cell_id),]
head(base2[,c(1:3)]); head(aclim[,c(1:3)])
# dim(base2); dim(aclim)
colnames(base2)
ddf <- cbind(aclim[which(aclim$id %in% unique(base2$cell_id)),], base2[which(base2$cell_id %in% unique(aclim$id)),c(868:870)] ) 
dim(ddf)
summary(ddf)

### Compute metabolic energy availability and diversity anomalies
# Add K and compute eV
ddf$absT <- ddf$SST + 273.15
ddf$kT <- 1 / (ddf$absT*-8.6173303*10^(-5) )

# Compute div anomalies to SST based on a linear regression
lm <- lm(log1p(rich_phyto) ~ kT, data = ddf, na.action = na.exclude)
summary(lm) # Adjusted R-squared:  0.2495 
#ddf$phyto_pred <- predict(lm)
ddf$phyto_anom <- residuals(lm)

lm <- lm(log1p(rich_zoo) ~ kT, data = ddf, na.action = na.exclude)
summary(lm) # Adjusted R-squared:  0.5385 
ddf$zoo_anom <- residuals(lm)

### Map and plot anomalies
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = phyto_anom), data = ddf) +
 	scale_fill_gradient2(name = "Anomalies to\nlinear model",
    low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-3.3,2.4)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = zoo_anom), data = ddf) +
 	scale_fill_gradient2(name = "Anomalies to\nlinear model",
    low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(-3.3,2.4)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
setwd(WD)
ggsave(plot = map4, filename = "map_phyto_rich_binom_anoms.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = map5, filename = "map_zoo_rich_binom_anoms.pdf", dpi = 300, width = 7, height = 4)

### And biplots
summary(log1p(ddf$rich_phyto)) # 0-5.6
summary(log1p(ddf$rich_zoo)) # 1.6-5.5
plot4 <- ggplot() + geom_point(aes(x = kT, y = log1p(rich_phyto)), data = ddf, colour = "grey70") +
            geom_smooth(aes(x = kT, y = log1p(rich_phyto)), data = ddf, colour = "black", method = "lm") + 
            scale_y_continuous(limits = c(0,5.6)) + ylab("Phytoplankton species richness\n(ln)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
plot5 <- ggplot() + geom_point(aes(x = kT, y = log1p(rich_zoo)), data = ddf, colour = "grey70") +
            geom_smooth(aes(x = kT, y = log1p(rich_zoo)), data = ddf, colour = "black", method = "lm") + 
            scale_y_continuous(limits = c(0,5.6)) + ylab("Zooplankton species richness\n(ln)") + 
            xlab("Thermal energy (1/kT)") + theme_classic()
#
ggsave(plot = plot4, filename = "plot_phyto_rich_binom_anoms.pdf", dpi = 300, width = 4.5, height = 4)
ggsave(plot = plot5, filename = "plot_zoo_rich_binom_anoms.pdf", dpi = 300, width = 4.5, height = 4)


### Plot heatmap of rank correlations
library("corrplot")
colnames(ddf)
# Compute p-val matrix
p.mat <- cor.mtest(na.omit(ddf[,c(5:17,22,18,19,20,23,24)]))$p
# Correlation matrix
M <- cor(na.omit(ddf[,c(5:17,22,18,19,20,23,24)]), method = "spearman")
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
levels(m$div)[levels(m$div) == "phyto_anom"] <- "Phytoplankton anomalies"
levels(m$div)[levels(m$div) == "zoo_anom"] <- "Zooplankton anomalies"

### Do the same for p.mat
p <- melt(data = p.mat, value.name = "pval", id.vars = colnames(ddf)[c(5:17,22)] )
head(p)
# remove the line where Var2 %in% colnames(aclim)[6:17]
p <- p[!(p$X1 %in% colnames(ddf)[c(5:17,22)]),]
p <- p[(p$X2 %in% colnames(ddf)[c(5:17,22)]),]
colnames(p) <- c("env","div","pval")
# Change the levels of the factors according to what you want to display on the plot
unique(p$env)
unique(p$div)
# For re-naming
levels(p$env)[levels(p$env) == "logEKE"] <- "EKE"
levels(p$env)[levels(p$env) == "NO3"] <- "Nitrates"
levels(p$env)[levels(p$env) == "SiO2"] <- "Silicates"
levels(p$env)[levels(p$env) == "Nstar"] <- "N*"
levels(p$env)[levels(p$env) == "Sistar"] <- "Si*"
levels(p$env)[levels(p$env) == "Chla"] <- "Chlorophyll"
levels(p$env)[levels(p$env) == "kT"] <- "kT"
levels(p$div)[levels(p$div) == "rich_plankton"] <- "Plankton richness"
levels(p$div)[levels(p$div) == "rich_phyto"] <- "Phytoplankton richness"
levels(p$div)[levels(p$div) == "rich_zoo"] <- "Zooplankton richness"
levels(p$div)[levels(p$div) == "phyto_anom"] <- "Phytoplankton anomalies"
levels(p$div)[levels(p$div) == "zoo_anom"] <- "Zooplankton anomalies"

### Change p$pval levels to "*" etc.
summary(p$pval)
# Humm a bit doubtful that ALL rank corr have a p-value < 0.0125794

#define a colour for fonts
textcol <- "grey30"
# Diverging color palette
levels(m$corrfactor) ; #unique(m$corrfactor) 11 levels but only 10 are realized (no 0)
#values <- c("#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
#quartz()
h1 <- ggplot(m,aes(x = env, y = div, fill = corrfactor)) + geom_tile() + 
    geom_tile(colour = "white", size = 0.25, show_guide = F) +
	labs(x = "",y = "",title = "Rank correlations") + scale_fill_manual(values = cmocean("curl")(15)[4:13], na.value = "grey90") +
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
ggsave(plot = h1, filename = "heatmap_rank_corr_binoms.pdf", dpi = 300, width = 12, height = 6)

### Examine how zoo and phyto diversity relate to Phytoplankton biomass
plot <- ggplot() + geom_point(aes(x = Chla, y = rich_phyto, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Phytoplankton species richness") + 
     scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_phytoxchla_rich_binom.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = Chla, y = rich_zoo, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Zooplankton species richness") + 
      scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_zooxchla_rich_binom.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = rich_phyto, y = rich_zoo, colour = Chla, size = abs(y)), data = ddf, alpha = 0.5) + 
         scale_colour_viridis(name = "Phytoplankton biomass\nlog(mgC.m3)") + 
         xlab("Phytoplankton species richness") + ylab("Zooplankton species richness") + 
         theme_classic()
ggsave(plot = plot, filename = "plot_phytoxzooxchla_rich_binom.pdf", dpi = 300, width = 8, height = 6)


### Nice ! Same with anomalies of richness
plot <- ggplot() + geom_point(aes(x = Chla, y = phyto_anom, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Phytoplankton richness anomalies") + 
     scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_phytoxchla_rich_binom_anoms.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = Chla, y = zoo_anom, colour = abs(y)), data = ddf) + 
     xlab("Phytoplankton biomass - log(mgC.m3)") + ylab("Zooplankton richness anomalies") + 
      scale_colour_viridis(name = "Latitude", option = "A") + 
     theme_classic()
ggsave(plot = plot, filename = "plot_zooxchla_rich_binom_anoms.pdf", dpi = 300, width = 6, height = 4)
#
plot <- ggplot() + geom_point(aes(x = phyto_anom, y = zoo_anom, colour = Chla, size = abs(y)), data = ddf, alpha = 0.5) + 
         scale_colour_viridis(name = "Phytoplankton biomass\nlog(mgC.m3)") + 
         xlab("Phytoplankton richness anomalies") + ylab("Zooplankton richness anomalies") + 
         theme_classic()
ggsave(plot = plot, filename = "plot_phytoxzooxchla_rich_binom_anoms.pdf", dpi = 300, width = 8, height = 6)


# --------------------------------------------------------------------------------------------------------------------------------

### 4°) Now, examine those PDR per plankton metacommunities

### Option A) Simply separate into Tropical/ Temperate and Polar metacomm
ddf$domain <- NA
ddf[which(abs(ddf$y) < 30),"domain"] <- "Tropical (<30°)"
ddf[which(abs(ddf$y) >= 30 & abs(ddf$y) < 60),"domain"] <- "Temperate (30°-60°)"
ddf[which(abs(ddf$y) >= 60),"domain"] <- "Polar (>60°)"
# levels(factor(zoo.fut$domain))
colnames(ddf)
ddf2 <- ddf[,c(2,3,5:17,19,20,25)]
colnames(ddf2)
m2 <- melt(data = ddf2, id.vars = colnames(ddf2)[c(1:15,18)])
colnames(m2)[17] <- c("div.var")
# Correct factor levels
levels(m2$div.var)[levels(m2$div.var) == "rich_phyto"] <- "Phytoplankton"
levels(m2$div.var)[levels(m2$div.var) == "rich_zoo"] <- "Zooplankton"

# Based on k4_ugpma 
p1 <- ggplot( m2[!is.na(m2$domain),] ) + geom_point(aes(x = Chla, y = value, colour = factor(domain)), alpha = 0.5 ) + 
    scale_colour_manual(name = "Metacommunity",
        values = c("#1f78b4","#a6cee3","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
	geom_smooth(aes(x = Chla, y = value), colour = "black", method = "lm") + 
	ylab("Species richness") + xlab("Phytoplankton biomass") + 
	theme_bw() + facet_grid(factor(div.var) ~ factor(domain) ) 
    
ggsave(plot = p1, filename = "plot_PDR_rich_domain_binom_fit.pdf", dpi = 300, width = 10, height = 4)

# Perform PCA or CA
library("FactoMineR")
colnames(base2)

### Check if a species as only zeroes
#colSums <- colSums(as.matrix(base2[,c(4:865)]))
#summary(colSums)
#colSums[colSums < 5] # Metridia_venusta; 618

pca <- PCA(X = base[,c(4:865)], graph = F, ncp = 5, scale.unit = F)
# Use cell's coords along CA1:CA4 to compute euclidean distances and cluster
summary(pca)
#                        Dim.1   Dim.2   Dim.3   Dim.4   Dim.5
# Variance              19.531   3.809   2.302   1.640   1.240 
# % of var.             59.143  11.535   6.970   4.967   3.755
# Cumulative % of var.  59.143  70.678  77.648  82.614  86.370                    
# str(pca)
coords <- data.frame(pca$ind$coord)
colnames(coords) <- c("CA1","CA2","CA3","CA4","CA5") #,"CA6","CA7","CA8","CA9","CA10")
# Compute euclidean distance using the CA scores
dist <- dist(coords[,c(1:5)], "euclidean")
# Perform HAC with Ward's algorithm
fit1 <- stats::hclust(dist, "ward.D2")
fit2 <- stats::hclust(dist, "average")

# Save plot of dendrogram
pdf(paste("dendro_ward_annual_plankton_comp_PCA.pdf", sep = ""), width = 15, height = 10)
plot(fit1)
dev.off()
pdf(paste("dendro_UGPMA_annual_plankton_comp_PCA.pdf", sep = ""), width = 15, height = 10)
plot(fit2)
dev.off()

# Based on Ward's
#base2$k2_ward <- cutree(fit1, 2)
base2$k3_ward <- cutree(fit1, 3)
base2$k4_ward <- cutree(fit1, 4)
base2$k5_ward <- cutree(fit1, 5)
#base2$k6_ward <- cutree(fit1, 6)
# Based on Ward's
#base2$k2_ugpma <- cutree(fit2, 2)
base2$k3_ugpma <- cutree(fit2, 3)
base2$k4_ugpma <- cutree(fit2, 4)
base2$k5_ugpma <- cutree(fit2, 5)
#base2$k6_ugpma <- cutree(fit2, 6)

### Ok, stick to k3_ward (Tropical.Temperate/Polar) or k4_ugpma (has upwellings)
head(base2[,c(1:3)])
head(ddf[,c(1:3)])
ddf$k4_ugpma <- factor(base[which(base2$cell_id %in% unique(ddf$id)),c("k4_ugpma")])
ddf$k3_ward <- factor(base[which(base2$cell_id %in% unique(ddf$id)),c("k3_ward")])
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

bioregions <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(k3_ward)), data = base2) +
 	scale_fill_manual(name = "Metacommunity", values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

setwd(WD)
ggsave(plot = bioregions, filename = "map_metacomm_k3_ward.pdf", dpi = 300, width = 6, height = 3)






