
# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")
library("viridis")
library("scales")
library("maps")
library("betapart")
library("cmocean")
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

### In a mclapply, compute changes in annual alpha+beta diversity for each combination of: SDM x pool x ESM
### (maps for annual 80 combinations)
# esm <- "CESM-BEC"
# sdm <- "GAM"
# p <- "p2"

diversity.changer <- function(esm = ESMs) {
    
                setwd(WD)
                message(paste("", sep = ""))
                message(paste("Computing changes in diversity for ",esm, sep = ""))
                message(paste("", sep = ""))
                 
                for(sdm in SDMs) {
                    
                        message(paste("Based on ",sdm, sep = ""))
                    
                        for(p in pools) {
                            
                            message(paste("and pool = ",p, sep = ""))
                            setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
                            require("dplyr")
                            
                            # Load monthly baseline projections 
                            basies.zoo <- lapply(months, function(m) {
                                        b <- read.table(paste("table_zoo_mon_composition+div_baseline_",p,"_",sdm,"_",m,".txt", sep=""), h = T, sep = "\t")
                                        #b$month <- m
                                        return(b)
                                } # eo FUN
                            ) # eo lapply - basies.zoo
                            # Rbind
                            base.zoo <- bind_rows(basies.zoo)
                            
                            basies.phyto <- lapply(months, function(m) {
                                        b <- read.table(paste("table_phyto_mon_composition+div_baseline_",p,"_",sdm,"_",m,".txt", sep=""), h = T, sep = "\t")
                                        #b$month <- m
                                        return(b) 
                                } # eo FUN
                            ) # eo lapply - basies.phyto
                            # Rbind
                            base.phyto <- bind_rows(basies.phyto)
                             
                            # Load monthly future projections
                            future.zoo <- lapply(months, function(m) {
                                        b <- read.table(paste("table_zoo_mon_composition+div_2100-2000_",esm,"_",p,"_",sdm,"_",m,".txt", sep=""), h = T, sep = "\t")
                                        #b$month <- m
                                        return(b) 
                                } # eo FUN
                            ) # eo lapply - future.zoo
                            # Rbind
                            fut.zoo <- bind_rows(future.zoo)
                            
                            future.phyto <- lapply(months, function(m) {
                                        b <- read.table(paste("table_phyto_mon_composition+div_2100-2000_",esm,"_",p,"_",sdm,"_",m,".txt", sep=""), h = T, sep = "\t")
                                        #b$month <- m
                                        return(b) 
                                } # eo FUN
                            ) # eo lapply - future.phyto
                            # Rbind
                            fut.phyto <- bind_rows(future.phyto)
                            
                            # dim(base.zoo); dim(base.phyto); dim(fut.zoo); dim(fut.phyto)
                            rm(future.phyto, future.zoo, basies.phyto, basies.zoo) ; gc()
                            
                            ### Compute annual average diversity (H' and rich) and mean % change in SR
                            message(paste("Calculating changes in richness", sep = ""))
                            ann.base.phyto <- data.frame(base.phyto %>% group_by(cell_id) %>%
                                    summarize(x = unique(x), y = unique(y), rich = mean(rich, na.rm = T)) 
                            )
                            
                            ann.base.zoo <- data.frame(base.zoo %>% group_by(cell_id) %>%
                                    summarize(x = unique(x), y = unique(y), rich = mean(rich, na.rm = T)) 
                            )
                            
                            ann.fut.phyto <- data.frame(fut.phyto %>% group_by(cell_id) %>%
                                    summarize(x = unique(x), y = unique(y), rich = mean(rich, na.rm = T)) 
                            )
                            
                            ann.fut.zoo <- data.frame(fut.zoo %>% group_by(cell_id) %>%
                                    summarize(x = unique(x), y = unique(y), rich = mean(rich, na.rm = T)) 
                            )
                            
                            # dim(ann.base.phyto); dim(ann.base.zoo); dim(ann.fut.phyto); dim(ann.fut.zoo)
                            base <- data.frame(id = ann.base.phyto$cell_id, x = ann.base.phyto$x, y = ann.base.phyto$y,
                                        rich_phyto = ann.base.phyto$rich, rich_zoo = ann.base.zoo$rich) 
                                        
                            fut <- data.frame(id = ann.fut.zoo$cell_id, x = ann.fut.zoo$x, y = ann.fut.zoo$y,
                                        rich_phyto = ann.fut.phyto$rich, rich_zoo = ann.fut.zoo$rich) 
                            
                            base$rich <- (base$rich_phyto)+(base$rich_zoo)
                            fut$rich <- (fut$rich_phyto)+(fut$rich_zoo)
                            
                            # # Maps, first choose appropriate limits for scale
 #                            min <- floor(min(c(base$rich_phyto,base$rich_zoo), na.rm = T))
 #                            max <- floor(max(c(base$rich_phyto,base$rich_zoo), na.rm = T))
 #
 #                            map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_phyto), data = base) +
 #                                scale_fill_viridis(name = "Annual richness", limits = c(min,max) ) +
 #                                geom_contour(colour = "grey75", binwidth = 50, size = 0.25, aes(x = x, y = y, z = rich_phyto), data = base) +
 #                                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 #                                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
 #                                           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 #                                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 #                                           labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
 #                                   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 #
 #                            map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich_zoo), data = base) +
 #                                scale_fill_viridis(name = "Annual richness", limits = c(min,max) ) +
 #                                geom_contour(colour = "grey75", binwidth = 50, size = 0.25, aes(x = x, y = y, z = rich_zoo), data = base) +
 #                                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 #                                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
 #                                           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 #                                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 #                                           labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
 #                                   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 #
 #                            map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = base) +
 #                                scale_fill_viridis(name = "Annual richness") + geom_contour(colour = "grey75", binwidth = 50,
 #                                            size = 0.25, aes(x = x, y = y, z = rich), data = base) +
 #                                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 #                                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
 #                                           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 #                                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 #                                           labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
 #                                   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 #
 #                            setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/maps")
 #                            ggsave(plot = map1, filename = paste("map_annual_rich_phyto_baseline_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
 #                            ggsave(plot = map2, filename = paste("map_annual_rich_zoo_baseline_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
 #                            ggsave(plot = map3, filename = paste("map_annual_rich_tot_baseline_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
 #                    
                            # Calculate % diff rich for all 3 groups
                            base <- base[order(base$id),]
                            fut <- fut[order(fut$id),]
                            
                            base$fut_rich_tot <- fut[which(fut$id %in% base$id),c("rich")]
                            base$fut_rich_phyto <- fut[which(fut$id %in% base$id),c("rich_phyto")]
                            base$fut_rich_zoo <- fut[which(fut$id %in% base$id),c("rich_zoo")]
                            
                            ### Save base (you'll need it to compute Mann-Whitney tests)
                            setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/")
                            save(base, file = paste("table_base+fut_SR_",esm,"_",sdm,"_",p,".Rdata", sep = ""))             
      
                            # changes <- data.frame(id = base$id, x = base$x, y = base$y,
 #                                diff_rich_tot = (base$fut_rich_tot)-(base$rich),
 #                                diff_rich_phyto = (base$fut_rich_phyto)-(base$rich_phyto),
 #                                diff_rich_zoo = (base$fut_rich_zoo)-(base$rich_zoo)
 #                            ) # eo ddf
 #                            # summary(changes)
 #                            # Add percentage changes
 #                            changes$perc_rich_tot <- ((changes$diff_rich_tot)/base$rich)*100
 #                            changes$perc_rich_phyto <- ((changes$diff_rich_phyto)/base$rich_phyto)*100
 #                            changes$perc_rich_zoo <- ((changes$diff_rich_zoo)/base$rich_zoo)*100
 #
 #                            min <- floor(min(changes[,c("perc_rich_tot","perc_rich_phyto","perc_rich_zoo")], na.rm = T) )
 #
 #                            map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_tot), data = changes[changes$perc_rich_tot < 50,]) +
 #                                geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_tot >= 50,], fill = "#b2182b") +
 #                                geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_tot),
 #                                        data = changes[changes$perc_rich_tot < 50,] ) +
 #                                 scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 #                                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 #                                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
 #                                           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 #                                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 #                                           labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
 #                                   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 #
 #                            map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_phyto), data = changes[changes$perc_rich_phyto < 50,]) +
 #                                geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_phyto >= 50,], fill = "#b2182b") +
 #                                geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_phyto),
 #                                        data = changes[changes$perc_rich_phyto < 50,] ) +
 #                                 scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 #                                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 #                                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
 #                                           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 #                                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 #                                           labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
 #                                   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 #
 #                            map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_zoo), data = changes[changes$perc_rich_zoo < 50,]) +
 #                                geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_zoo >= 50,], fill = "#b2182b") +
 #                                geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_zoo),
 #                                        data = changes[changes$perc_rich_zoo < 50,] ) +
 #                                 scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 #                                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 #                                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
 #                                           labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 #                                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 #                                           labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
 #                                   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 #
 #                            ggsave(plot = map1, filename = paste("map_annual_perc_rich_tot_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
 #                            ggsave(plot = map2, filename = paste("map_annual_perc_rich_phyto_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
 #                            ggsave(plot = map3, filename = paste("map_annual_perc_rich_zoo_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
 #
 #                            ### Save changes in richness (you'll compute ensemble from those)
 #                            setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/")
 #                            save(changes, file = paste("table_ann_changes_",esm,"_",sdm,"_",p,".Rdata", sep = ""))
 #                            # Make room for change sin beta div
 #                            rm(ann.base.phyto,ann.base.zoo,ann.fut.phyto,ann.fut.zoo,base,fut,map1,map2,map3,changes); gc()
 #
 #                            ### Computing changes in beta diversity
 #                            message(paste("Computing species mean annual composition", sep = ""))
 #                            # Compute mean annual species' HSI (melt+dcast) Remove NAs?
 #                            # dim(base.zoo); dim(base.phyto); dim(fut.zoo); dim(fut.phyto)
 #                            base.phyto <- na.omit(base.phyto)
 #                            base.zoo <- na.omit(base.zoo)
 #                            fut.phyto <- na.omit(fut.phyto)
 #                            fut.zoo <- na.omit(fut.zoo)
 #                            # Melt
 #                            m.base.phyto <- melt(base.phyto[,c(1:c(length(base.phyto)-2))], id.vars = c("cell_id","x","y") )
 #                            m.base.zoo <- melt(base.zoo[,c(1:c(length(base.zoo)-2))], id.vars = c("cell_id","x","y") )
 #                            m.fut.phyto <- melt(fut.phyto[,c(1:c(length(fut.phyto)-2))], id.vars = c("cell_id","x","y") )
 #                            m.fut.zoo <- melt(fut.zoo[,c(1:c(length(fut.zoo)-2))], id.vars = c("cell_id","x","y") )
 #                            colnames(m.base.zoo)[c(4,5)] <- c("species","HSI")
 #                            colnames(m.base.phyto)[c(4,5)] <- c("species","HSI")
 #                            colnames(m.fut.phyto)[c(4,5)] <- c("species","HSI")
 #                            colnames(m.fut.zoo)[c(4,5)] <- c("species","HSI")
 #                            # Dcast with averaging as a FUN
 #                            d.base.zoo <- dcast(m.base.zoo, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
 #                            d.base.phyto <- dcast(m.base.phyto, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
 #                            d.fut.zoo <- dcast(m.fut.zoo, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
 #                            d.fut.phyto <- dcast(m.fut.phyto, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
 #                            # dim(d.base.zoo); dim(d.base.phyto); dim(d.fut.zoo); dim(d.fut.phyto)
 #
 #                            # Re-order to make sure they can be combined
 #                            d.base.zoo <- d.base.zoo[order(d.base.zoo$cell_id),]
 #                            d.base.phyto <- d.base.phyto[order(d.base.phyto$cell_id),]
 #                            d.fut.zoo <- d.fut.zoo[order(d.fut.zoo$cell_id),]
 #                            d.fut.phyto <- d.fut.phyto[order(d.fut.phyto$cell_id),]
 #                            # Make room
 #                            rm(m.base.phyto,m.base.zoo,m.fut.phyto,m.fut.zoo); gc()
 #
 #                            # Save annual comp tables too (useful for species' range shifts)
 #                            save(d.base.zoo, file = paste("table_ann_compo_zoo_baseline_",sdm,"_",p,".Rdata", sep = ""))
 #                            save(d.base.phyto, file = paste("table_ann_compo_phyto_baseline_",sdm,"_",p,".Rdata", sep = ""))
 #                            save(d.fut.zoo, file = paste("table_ann_compo_zoo_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = ""))
 #                            save(d.fut.phyto, file = paste("table_ann_compo_phyto_2100-2000_",esm,"_",sdm,"_",p,".Rdata", sep = ""))
 #                         
                           #  # Bind 'em
                           #  base <- cbind(d.base.phyto[which(d.base.phyto$cell_id %in% unique(d.base.zoo$cell_id)),],
                           #              d.base.zoo[which(d.base.zoo$cell_id %in% unique(d.base.phyto$cell_id)),c(4:length(d.base.zoo))]
                           #          )
                           #  fut <- cbind(d.fut.phyto[which(d.fut.phyto$cell_id %in% unique(d.fut.zoo$cell_id)),],
                           #              d.fut.zoo[which(d.fut.zoo$cell_id %in% unique(d.fut.phyto$cell_id)),c(4:length(d.fut.zoo))]
                           #          )
                           #  # Remove points in the colnames
                           #  colnames(base)[c(4:865)] <- gsub("[.]","",as.character(colnames(base)[c(4:865)]))
                           #  colnames(fut)[c(4:865)] <- gsub("[.]","",as.character(colnames(fut)[c(4:865)]))
                           #  # dim(base); dim(fut)
                           #  base <- na.omit(base)
                           #  fut <- na.omit(fut)
                           #  base <- base[which(base$cell_id %in% unique(fut$cell_id)),]
                           #  fut <- fut[which(fut$cell_id %in% unique(base$cell_id)),]
                           #
                           #  # Calculating beta div changes for various thresholds
                           #  thresholds <- seq(from = 0.35, to = 0.45, by = 0.01)
                           #  # t <- 0.35
                           #  message(paste("Calculating changes in beta div", sep = ""))
                           #  thresh <- lapply(thresholds, function(t) {
                           #
                           #              base2 <- base
                           #              fut2 <- fut
                           #              # Get months names based on couples
                           #              message(paste("Using t = ",t, sep = ""))
                           #              # Convert to P/1 using rbinom and initial probability
                           #              for(sp in colnames(base2)[c(4:865)] ) {
                           #                      base2[,c(sp)][base2[,c(sp)] > t] <- 1
                           #                      base2[,c(sp)][base2[,c(sp)] <= t] <- 0
                           #                      # Future 2 now
                           #                      fut2[,c(sp)][fut2[,c(sp)] > t] <- 1
                           #                      fut2[,c(sp)][fut2[,c(sp)] <= t] <- 0
                           #              } # eo for loop
                           #              # And compute beta.div changes
                           #              M = nrow(base2[,c(4:865)])
                           #              N = ncol(base2[,c(4:865)])
                           #
                           #              ### With paralelling :
                           #              require("doParallel")
                           #              require("plyr")
                           #              registerDoParallel(cores = 15)
                           #              # Compute beta.div changes ofr all plankton
                           #              div <- data.frame()
                           #              d <- cbind(base2[,c(4:865)], fut2[,c(4:865)])
                           #              d$bit <- cut(1:M, 50, labels = F)
                           #              div <- ddply(d, ~ bit, function(x) {
                           #                           beta.temp(x[,1:N], x[,(N+1):(2*N)],"jaccard")
                           #                      },.parallel = T
                           #              ) # eo ddply
                           #              # summary(div)
                           #
                           #              # For phytoplankton
                           #              M = nrow(base2[,c(4:341)])
                           #              N = ncol(base2[,c(4:341)])
                           #              div.phyto <- data.frame()
                           #              d <- cbind(base2[,c(4:341)], fut2[,c(4:341)])
                           #              d$bit <- cut(1:M, 50, labels = F)
                           #              # head(d)
                           #              div.phyto <- ddply(d, ~ bit, function(x) {
                           #                           beta.temp(x[,1:N], x[,(N+1):(2*N)],"jaccard")
                           #                      },.parallel = T
                           #              ) # eo ddply
                           #
                           #              # For zooplankton
                           #              M = nrow(base2[,c(342:865)])
                           #              N = ncol(base2[,c(342:865)])
                           #              div.zoo <- data.frame()
                           #              d <- cbind(base2[,c(342:865)], fut2[,c(342:865)])
                           #              d$bit <- cut(1:M, 50, labels = F)
                           #              div.zoo <- ddply(d, ~ bit, function(x) {
                           #                           beta.temp(x[,1:N], x[,(N+1):(2*N)],"jaccard")
                           #                      },.parallel = T
                           #              ) # eo ddply
                           #
                           #              if( sum(is.na(div$beta.jtu)) >= 1 ) {
                           #                  # summary(div.phyto)
                           #                  div[which(is.na(div$beta.jac)),c("beta.jac","beta.jne","beta.jtu")] <- 0
                           #                  div[which(is.na(div$beta.jne) & !is.na(div$beta.jac)),c("beta.jne")] <- 1
                           #                  div[which(is.na(div$beta.jtu) & !is.na(div$beta.jac)),c("beta.jtu")] <- 1
                           #              }
                           #
                           #              if( sum(is.na(div.phyto$beta.jtu)) >= 1 ) {
                           #                  # summary(div.phyto)
                           #                  div.phyto[which(is.na(div.phyto$beta.jac)),c("beta.jac","beta.jne","beta.jtu")] <- 0
                           #                  div.phyto[which(is.na(div.phyto$beta.jne) & !is.na(div.phyto$beta.jac)),c("beta.jne")] <- 1
                           #                  div.phyto[which(is.na(div.phyto$beta.jtu) & !is.na(div.phyto$beta.jac)),c("beta.jtu")] <- 1
                           #              }
                           #
                           #              if( sum(is.na(div.zoo$beta.jtu)) >= 1 ) {
                           #                  # summary(div.phyto)
                           #                  div.zoo[which(is.na(div.zoo$beta.jac)),c("beta.jac","beta.jne","beta.jtu")] <- 0
                           #                  div.zoo[which(is.na(div.zoo$beta.jne) & !is.na(div.zoo$beta.jac)),c("beta.jne")] <- 1
                           #                  div.zoo[which(is.na(div.zoo$beta.jtu) & !is.na(div.zoo$beta.jac)),c("beta.jtu")] <- 1
                           #              }
                           #
                           #              # Gather all in a data.frame and return
                           #              div$x <- base2$x
                           #              div$y <- base2$y
                           #              div$id <- factor(paste(div$x, div$y, sep = "_"))
                           #              # And add zoo and phyto components of beta.div
                           #              div$jac.phyto <- div.phyto$beta.jac
                           #              div$jne.phyto <- div.phyto$beta.jne
                           #              div$jtu.phyto <- div.phyto$beta.jtu
                           #              div$jac.zoo <- div.zoo$beta.jac
                           #              div$jne.zoo <- div.zoo$beta.jne
                           #              div$jtu.zoo <- div.zoo$beta.jtu
                           #              # Add the couple
                           #              div$thresh <- factor(t)
                           #              # Return(div)
                           #              rm(d); gc()
                           #              return(div)
                           #
                           #      } # eo 2nd FUN
                           #
                           #  )  # eo 2nd lapply - annual
                           #  # Rbind
                           #  detach("package:plyr", unload = T)
                           #  require("dplyr")
                           #  table <- dplyr::bind_rows(thresh)
                           #  rm(thresh); gc()
                           #  # Compute average changes across thresholds
                           #  div <- data.frame(table %>%
                           #          group_by(id) %>%
                           #          summarise(x = unique(x), y = unique(y),
                           #          jac = mean(beta.jac,na.rm=T), jne = mean(beta.jne,na.rm=T), jtu = mean(beta.jtu,na.rm=T),
                           #          jac_phyto = mean(jac.phyto,na.rm=T), jne_phyto = mean(jne.phyto,na.rm=T), jtu_phyto = mean(jtu.phyto,na.rm=T),
                           #          jac_zoo = mean(jac.zoo,na.rm=T), jne_zoo = mean(jne.zoo,na.rm=T), jtu_zoo = mean(jtu.zoo,na.rm=T),
                           #          diff_all = mean(diff_all,na.rm=T), diff_phyto = mean(diff_phyto,na.rm=T), diff_zoo = mean(diff_zoo,na.rm=T)
                           #          ) # eo summarise
                           #  ) # eo ddf
                           #
                           #  message(paste("Saving changes in beta div", sep = ""))
                           #  rm(table); gc()
                           #  setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/")
                           #  save(div, file = paste("table_ann_changes_beta.div_",esm,"_",sdm,"_",p,".Rdata", sep = ""))
                           #
                           #  setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/maps")
                           #  # Annual plankton jac
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = div) +
                           #      scale_fill_viridis(name = "Annual Jaccard", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jac_tot_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  # Annual plankton nestedness
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = div) +
                           #      scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jne_tot_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  # Annual plankton turn-over
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = div) +
                           #      scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jtu_tot_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  ### Phyto
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = div) +
                           #      scale_fill_viridis(name = "Annual Jaccard index", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jac_phyto_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  # Annual Phyto nestedness
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = div) +
                           #      scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jne_phyto_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  # Annual Phyto turn-over
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = div) +
                           #      scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jtu_phyto_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #
                           #  ### Zoo
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = div) +
                           #      scale_fill_viridis(name = "Annual Jaccard index", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           # ggsave(plot = map, filename = paste("map_annual_jac_zoo_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = div) +
                           #      scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jne_zoo_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  # Annual Zoo turn-over
                           #  map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = div) +
                           #      scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
                           #      geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = div) +
                           #       geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                           #       coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                           #                 labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                           #       scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                           #                 labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                           #         theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           #           panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                           #
                           #  ggsave(plot = map, filename = paste("map_annual_jne_zoo_",esm,"_",sdm,"_",p,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
                           #
                           #  ### Clean stuff and exit loop
                           #  rm(base, fut, d.base.phyto, d.base.zoo, d.fut.phyto, d.fut.zoo, map)
                           #  gc()
                           
                            setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/")
                            
                        } # eo for loop - sdm in SDMs
                    
                } # eo for loop - sdm in SDMs
    
} # eo FUN - diversity.changer

require("parallel")
mclapply(X = ESMs, FUN = diversity.changer, mc.cores = 5)


# --------------------------------------------------------------------------------------------------------------------------------


