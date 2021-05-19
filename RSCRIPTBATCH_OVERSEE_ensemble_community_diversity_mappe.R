
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
library("ncdf4")

world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------

### Set the working directories, vectors etc.
WD <- getwd()
setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")

months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

### 1°) First, map all monthly richness and shannon-weiver indices values (baseline and future)
# files1 <- dir()[grep("_compo_baseline_", dir())]; files1 # baseline
# files2 <- dir()[grep("_compo_2100-2000", dir())]; files2 # baseline
# Examine some txt files and check
# Test
# f <- "table_phyto_mean_compo_baseline_apr.txt"
# f <- "table_phyto_mean_compo_2100-2000_MRI-NEMURO_rcp85_apr.txt"
# diversity.mapper <- function(f = files1) {
#
#         if( f %in% files1 ) {
#
#             # Useless message
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#             group <- do.call(cbind,strsplit(as.character(f), split = "_"))[2,1]
#             if(group == "zoo") { group2 <- "Zooplankton" }
#             if(group == "phyto") { group2 <- "Phytoplankton" }
#             m <- do.call(cbind,strsplit(as.character(f), split = "_"))[6,1]
#             m <- str_replace(m,".txt","") # m
#             message(paste("Mapping baseline SR and H' for ",m,  sep = ""))
#             table <- read.table(f, sep = "\t")
#             # summary(table[,c("x","y")])
#             # Maps!
#             map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = na.omit(table)) +
#              scale_fill_viridis(name = paste(group2,"\n diversity (SR)",sep="")) +
#              geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#              coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#              scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#              theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#
#             map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = H), data = na.omit(table)) +
#                 scale_fill_viridis(name = paste(group2,"\n diversity (H)",sep="")) +
#                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                              panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#             ### Save
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps")
#             ggsave(plot = map1, filename = paste("map_rich_",group,"_base_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
#             ggsave(plot = map2, filename = paste("map_shannon_",group,"_base_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#
#         } else if ( f %in% files2 ) {
#
#             # Useless message
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#             group <- do.call(cbind,strsplit(as.character(f), split = "_"))[2,1]
#             if(group == "zoo") { group2 <- "Zooplankton" }
#             if(group == "phyto") { group2 <- "Phytoplankton" }
#             esm <- do.call(cbind,strsplit(as.character(f), split = "_"))[6,1]
#             m <- do.call(cbind,strsplit(as.character(f), split = "_"))[8,1]
#             m <- str_replace(m,".txt","") # m
#             message(paste("Mapping future SR and H' for ",m," under ",esm,  sep = ""))
#             table <- read.table(f, sep = "\t")
#
#             # Maps!
#             map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = na.omit(table)) +
#              scale_fill_viridis(name = paste(group2,"\n diversity (SR)",sep="")) +
#              geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#              coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#              scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#              theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#
#             map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = H), data = na.omit(table)) +
#                 scale_fill_viridis(name = paste(group2,"\n diversity (H)",sep="")) +
#                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                 coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#                 scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#                 theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                              panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#             ### Save
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps")
#             ggsave(plot = map1, filename = paste("map_rich_",group,"_2100-2000_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
#             ggsave(plot = map2, filename = paste("map_shannon_",group,"_2100-2000_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#
#         } # eo if else loop
#
# }
#
# require("parallel")
# mclapply(files2, FUN = diversity.mapper, mc.cores = 25)


# ---------------------------------------------------------------

### 09/03/2020: 2°) Same but on the annual scale: compute mean annual diversity (rich and H') and map

### First, for the baseline period
# With lapply, extract each monthly output and compute mean
# m <- "apr"
# require("parallel")
# res <- mclapply(months, function(m) {
#
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#             # Extract monthly files
#             message(paste("Extracting baseline monthly diversity indices for ",m, sep = ""))
#             table.phyto <- read.table(paste("table_phyto_mean_compo_baseline_",m,".txt",sep=""), sep = "\t")
#             table.zoo <- read.table(paste("table_zoo_mean_compo_baseline_",m,".txt",sep=""), sep = "\t")
#             # dim(table.phyto) ; dim(table.zoo)
#             ddf <- data.frame(id = table.phyto$cell_id, x = table.phyto$x, y = table.phyto$y, month = m,
#                 rich_phyto = table.phyto$rich, rich_zoo = table.zoo$rich, Hphyto = table.phyto$H, Hzoo = table.zoo$H)
#             # summary(ddf)
#             return(ddf)
#
#         }, mc.cores = 20
#
# ) # eo lapply
# # Rbind
# table <- dplyr::bind_rows(res)
# dim(table) ; head(table)
# ddf <- data.frame(table %>% group_by(id) %>% summarize(x = unique(x), y = unique(y),
#                 rich_phyto = mean(rich_phyto, na.rm=T), rich_zoo = mean(rich_zoo, na.rm=T),
#                 Hphyto = mean(Hphyto, na.rm=T), Hzoo = mean(Hzoo, na.rm=T) )
# ) # eo ddf
# summary(ddf)
#
# # Maps!
# map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = Hphyto), data = na.omit(ddf)) +
#  scale_fill_viridis(name = "Annual\nPhytoplankton\ndiversity (H')") +
#  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#  scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#  theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#
# map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = Hzoo), data = na.omit(ddf)) +
#     scale_fill_viridis(name = "Annual\nZooplankton\ndiversity (H')") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#     theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                  panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
# ### Save
# setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps")
# ggsave(plot = map1, filename = "map_annual_shannon_phyto_baseline.jpg", dpi = 300, width = 7, height = 5)
# ggsave(plot = map2, filename = "map_annual_shannon_zoo_baseline.jpg", dpi = 300, width = 7, height = 5)
# setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")


### Also, compute the species' mean ANNUAL HSI for baeline and future (per ESM) and save tables
# require("parallel")
# res <- mclapply(months, function(m) {
#
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#             # Extract monthly files
#             message(paste("Extracting plankton baseline monthly composition for ",m, sep = ""))
#             table.phyto <- read.table(paste("table_phyto_mean_compo_baseline_",m,".txt",sep=""), sep = "\t")
#             table.zoo <- read.table(paste("table_zoo_mean_compo_baseline_",m,".txt",sep=""), sep = "\t")
#             # Melt after removing the div metrics from the columns
#             m.table.phyto <- reshape2::melt(data = table.phyto[,c(1:341)], id.vars = colnames(table.phyto)[c(1:3)])
#             colnames(m.table.phyto)[c(4:5)] <- c("species","HSI")
#             m.table.zoo <- reshape2::melt(data = table.zoo[,c(1:527)], id.vars = colnames(table.zoo)[c(1:3)])
#             colnames(m.table.zoo)[c(4:5)] <- c("species","HSI")
#             m.table.phyto$group <- "Phytoplankton"
#             m.table.zoo$group <- "Zooplankton"
#             # rbind them, add month, and return
#             m.table <- rbind(m.table.phyto, m.table.zoo)
#             m.table$month <- m
#             return(m.table)
#
#         }, mc.cores = 20
#
# ) # eo lapply
# table <- dplyr::bind_rows(res)
# dim(table) ; head(table)
# rm(res); gc()
# # Compute mean annual species' level HSI
# ddf <- data.frame(table %>% group_by(cell_id,group,species) %>% summarize(x = unique(x), y = unique(y), mean = mean(HSI, na.rm = T)))
# head(ddf); dim(ddf)
# # Dcast and save
# d_ddf <- dcast(ddf[,c(1,3:6)], cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean")
# head(d_ddf); dim(d_ddf)
#
# # And map annual total diversity metrics
# require("vegan")
# d_ddf$rich <- rowSums(as.matrix(d_ddf[,c(4:length(d_ddf))]))
# d_ddf$H <- diversity(d_ddf[,c(4:(length(d_ddf)-1))], index = "shannon")
# summary(d_ddf)
#
# # Maps!
# map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = na.omit(d_ddf)) +
#  scale_fill_viridis(name = "Annual\nPlankton\ndiversity (H')") +
#  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#  scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#  theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#
# map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = H), data = na.omit(d_ddf)) +
#     scale_fill_viridis(name = "Annual\nPlankton\ndiversity (H')") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
#                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#                        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#     theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                  panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
# ### Save
# setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps")
# ggsave(plot = map1, filename = "map_annual_rich_plankton_baseline.jpg", dpi = 300, width = 7, height = 5)
# ggsave(plot = map2, filename = "map_annual_shannon_plankton_baseline.jpg", dpi = 300, width = 7, height = 5)
# # Save
# setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
# write.table(d_ddf, file = "table_annual_mean_compo_plankton_baseline.txt", sep = "\t")

### Same for ESMs (future annual mean HSI)
require("parallel")
# Vector of earth system models
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
rcp <- "rcp85"

# res <- mclapply(ESMs, function(esm) {
#
#             setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
#             res2 <- lapply(months, function(m) {
#                         # Extract monthly files
#                         message(paste("Extracting phytoplankton baseline monthly composition for ",esm,", for ",m, sep = ""))
#                         table.phyto <- read.table(paste("table_phyto_mon_composition_2100-2000_",esm,"_rcp85_",m,".txt",sep=""), sep = "\t")
#                         # Melt after removing the div metrics from the columns
#                         m.table.phyto <- reshape2::melt(data = table.phyto, id.vars = colnames(table.phyto)[c(1:4)] )
#                         colnames(m.table.phyto)[c(5:6)] <- c("species","HSI")
#                         m.table.phyto$month <- m
#                         # head(m.table.phyto); dim(m.table.phyto)
#                         return(m.table.phyto)
#                     } # eo FUN
#             ) # eo lapply
#             # Rbind
#             table <- dplyr::bind_rows(res2)
#             dim(table) ; head(table)
#             rm(res2); gc()
#
#             # Compute mean annual species' level HSI
#             message(paste("Computing phytoplankton species' mean annual HSI, for ",esm, sep = ""))
#             annual <- data.frame(table %>%
#                     group_by(cell_id,species) %>%
#                     summarize(x = unique(x), y = unique(y), mean = mean(HSI,na.rm=T))
#             ) # eo ddf
#             # head(annual); dim(annual)
#             # Dcast
#             d_annual <- dcast(annual, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean")
#             # head(d_annual); dim(d_annual)
#             # Save
#             message(paste("Saving table of phytoplankton species' mean annual HSI, for ",esm, sep = ""))
#             write.table(d_annual, file = paste("table_annual_mean_compo_phyto_2100-2000_rcp85_",esm,".txt", sep = ""), sep = "\t")
#             #rm(annual, table); gc()
#
#         }, mc.cores = 5
#
# ) # eo mclapply
# Check those vibes
# test <- read.table("table_annual_mean_compo_phyto_2100-2000_rcp85_GFDL-TOPAZ.txt", h = T, sep = "\t")
# dim(test)

### And now zooplankton
# esm <- "IPSL-PISCES"
# months <- c("jan","feb","mar")

res <- mclapply(ESMs, function(esm) {
            
            setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
            res2 <- lapply(months, function(m) {
                        # Extract monthly files
                        message(paste("Extracting zooplankton baseline monthly composition for ",esm,", for ",m, sep = ""))
                        table.zoo <- read.table(paste("table_zoo_mon_composition_2100-2000_",esm,"_rcp85_",m,".txt",sep=""), sep = "\t")
                        # Melt after removing the div metrics from the columns
                        m.table.zoo <- reshape2::melt(data = table.zoo, id.vars = colnames(table.zoo)[c(1:4)])
                        colnames(m.table.zoo)[c(5:6)] <- c("species","HSI")
                        m.table.zoo$month <- m
                        # head(m.table.zoo); dim(m.table.zoo)
                        return(m.table.zoo)
                    } # eo FUN
            ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res2)
            # dim(table) ; head(table)
            rm(res2); gc()
            
            # Compute mean annual species' level HSI
            message(paste("Computing zooplankton species' mean annual HSI, for ",esm, sep = ""))
            annual <- data.frame(table %>%
                    group_by(cell_id,species) %>%
                    summarize(x = unique(x), y = unique(y), mean = mean(HSI,na.rm=T))
            ) # eo ddf
            # head(annual); dim(annual)
            # Dcast
            d_annual <- dcast(annual, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "mean")
            # head(d_annual); dim(d_annual)
            # Save
            message(paste("Saving table of zooplankton species' mean annual HSI, for ",esm, sep = ""))
            write.table(d_annual, file = paste("table_annual_mean_compo_zoo_2100-2000_rcp85_",esm,".txt", sep = ""), sep = "\t")
            #rm(annual, table); gc()
            
        }, mc.cores = 5
    
) # eo mclapply


# ---------------------------------------------------------------

### 3°) Map the monthly and annual differences in richness and H'
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
# esm <- "GFDL-TOPAZ"
for(esm in ESMs) {

         setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
         message(paste("Loading annual community projections for ",esm, sep = ""))
         base <- read.table("table_annual_mean_compo_plankton_baseline.txt", h = T, sep = "\t")
         # remove last 2 cols
         base <- base[,c(1:(length(base)-2))]
         fut.phyto <- read.table(paste("table_annual_mean_compo_phyto_2100-2000_rcp85_",esm,".txt", sep = ""), h = T, sep = "\t")
         fut.zoo <- read.table(paste("table_annual_mean_compo_zoo_2100-2000_rcp85_",esm,".txt", sep = ""), h = T, sep = "\t")
         # dim(fut.phyto) ; dim(fut.zoo) ; colnames(fut.phyto) ; colnames(fut.zoo)
         base <- base[order(base$cell_id),]
         fut.phyto <- fut.phyto[order(fut.phyto$cell_id),]
         fut.zoo <- fut.zoo[order(fut.zoo$cell_id),]
         # Merge (cbind) fut.phyto and fut.zoo for full plankton diversity estimates
         fut <- cbind(fut.phyto, fut.zoo[,c(4:length(fut.zoo))]) # dim(fut)
         # 4:341 -> phyto; 342:865 -> zoopl
         
         # Derive multiple diversity metrics for each group (total/phyto/zoo)
         base$rich_total <- rowSums(as.matrix(base[,c(4:length(base))]))
         base$rich_phyto <- rowSums(as.matrix(base[,c(4:341)]))
         base$rich_zoo <- rowSums(as.matrix(base[,c(342:865)]))
         require("vegan")
         base$Htotal <- diversity(base[,c(4:865)], index = "shannon")
         base$Hphyto <- diversity(base[,c(4:341)], index = "shannon")
         base$Hzoo <- diversity(base[,c(342:865)], index = "shannon")
         
         # Same but for fut
         fut$rich_total <- rowSums(as.matrix(fut[,c(4:length(fut))]))
         fut$rich_phyto <- rowSums(as.matrix(fut[,c(4:341)]))
         fut$rich_zoo <- rowSums(as.matrix(fut[,c(342:865)]))
         fut$Htotal <- diversity(fut[,c(4:865)], index = "shannon")
         fut$Hphyto <- diversity(fut[,c(4:341)], index = "shannon")
         fut$Hzoo <- diversity(fut[,c(342:865)], index = "shannon")
         
         # And derive % changes ! head(fut[,c(1:4)]); head(base[,c(1:4)])
         message(paste("Calculating changes in diversity for ",esm, sep = ""))
         base$fut_rich_tot <- fut[which(fut$cell_id %in% base$cell_id),c("rich_total")]
         base$fut_rich_phyto <- fut[which(fut$cell_id %in% base$cell_id),c("rich_phyto")]
         base$fut_rich_zoo <- fut[which(fut$cell_id %in% base$cell_id),c("rich_zoo")]
         base$fut_Htotal <- fut[which(fut$cell_id %in% base$cell_id),c("Htotal")]
         base$fut_Hphyto <- fut[which(fut$cell_id %in% base$cell_id),c("Hphyto")]
         base$fut_Hzoo <- fut[which(fut$cell_id %in% base$cell_id),c("Hzoo")]
      
         changes <- data.frame(id = base$cell_id, x = base$x, y = base$y, 
             diff_rich_tot = (base$fut_rich_tot)-(base$rich_total), 
             diff_rich_phyto = (base$fut_rich_phyto)-(base$rich_phyto), 
             diff_rich_zoo = (base$fut_rich_zoo)-(base$rich_zoo),
             diff_Htot = (base$fut_Htotal)-(base$Htotal),
             diff_Hphyto = (base$fut_Hphyto)-(base$Hphyto),
             diff_Hzoo = (base$fut_Hzoo)-(base$Hzoo) 
         ) # eo ddf 
         # summary(changes)
         # Add percentage changes
         changes$perc_rich_tot <- ((changes$diff_rich_tot)/base$rich_total)*100
         changes$perc_rich_phyto <- ((changes$diff_rich_phyto)/base$rich_phyto)*100
         changes$perc_rich_zoo <- ((changes$diff_rich_zoo)/base$rich_zoo)*100
         changes$perc_Htot <- ((changes$diff_Htot)/base$Htotal)*100
         changes$perc_Hphyto <- ((changes$diff_Hphyto)/base$Hphyto)*100
         changes$perc_Hzoo <- ((changes$diff_Hzoo)/base$Hzoo)*100
         
         # Mapping time with cobtours and diverging color gradients
         message(paste("Mapping annual changes in diversity (richness and H') for ",esm, sep = ""))
         
         # Define the minimum limts for the richness diff limits of the 
         min <- floor(min(changes[,c("perc_rich_tot","perc_rich_phyto","perc_rich_zoo")], na.rm = T) )
         min2 <- floor(min(changes[,c("perc_Htot","perc_Hphyto","perc_Hzoo")], na.rm = T) )
         max2 <- abs(min2)
         
         map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_tot), data = changes[changes$perc_rich_tot < 50,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_tot >= 50,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_tot), 
                     data = changes[changes$perc_rich_tot < 50,] ) +
          	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_phyto), data = changes[changes$perc_rich_phyto < 50,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_phyto >= 50,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_phyto), 
                     data = changes[changes$perc_rich_phyto < 50,] ) +
          	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_zoo), data = changes[changes$perc_rich_zoo < 50,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_zoo >= 50,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_zoo), 
                     data = changes[changes$perc_rich_zoo < 50,] ) +
          	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Htot), data = changes[changes$perc_Htot < max2,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_Htot >= max2,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Htot), 
                     data = changes[changes$perc_Htot < max2,] ) +
          	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hphyto), data = changes[changes$perc_Hphyto < max2,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_Hphyto >= max2,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hphyto), 
                     data = changes[changes$perc_Hphyto < max2,] ) +
          	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hzoo), data = changes[changes$perc_Hzoo < max2,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_Hzoo >= max2,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hzoo), 
                     data = changes[changes$perc_Hzoo < max2,] ) +
          	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         
         # Save maps as .jpg for now
         setwd(paste(getwd(),"/","maps", sep = ""))
         ggsave(plot = map1, filename = paste("map_annual_perc_rich_tot_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map2, filename = paste("map_annual_perc_rich_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map3, filename = paste("map_annual_perc_rich_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map4, filename = paste("map_annual_perc_Htot_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map5, filename = paste("map_annual_perc_Hphyto_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map6, filename = paste("map_annual_perc_Hzoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         
         # Saving estimates
         setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
         message(paste("Saving diversity changes estimates for ",esm, sep = ""))
         message(paste("", sep = ""))
         write.table(changes, file = paste("table_annual_changes_",esm,".txt",sep=""), sep = "\t")   
         
} # eo for loop

### And monthly
# esm <- "GFDL-TOPAZ"
# m <- "apr"
for(esm in ESMs) {

     for(m in months) {
         
         setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
         message(paste("Loading monthly community projections for ",m," || ",esm, sep = ""))
         base.phyto <- read.table(paste("table_phyto_mean_compo_baseline_",m,".txt", sep = ""),)
         base.zoo <- read.table(paste("table_zoo_mean_compo_baseline_",m,".txt", sep = ""),)
         fut.phyto <- read.table(paste("table_phyto_mean_compo_2100-2000_",esm,"_rcp85_",m,".txt", sep = ""), h = T, sep = "\t")
         fut.zoo <- read.table(paste("table_zoo_mean_compo_2100-2000_",esm,"_rcp85_",m,".txt", sep = ""), h = T, sep = "\t")
         # colnames(base.phyto); colnames(base.zoo); colnames(fut.phyto); colnames(fut.zoo)
         # dim(base.phyto); dim(base.zoo); dim(fut.phyto); dim(fut.zoo)
         base.phyto <- base.phyto[order(base.phyto$cell_id),]; base.zoo <- base.zoo[order(base.zoo$cell_id),]
         fut.phyto <- fut.phyto[order(fut.phyto$cell_id),]; fut.zoo <- fut.zoo[order(fut.zoo$cell_id),]
         # Merge for total plankton diversity changes
         base <- cbind(base.phyto[,c(1:341)], base.zoo[,c(4:527)])
         fut <- cbind(fut.phyto[,c(1:341)], fut.zoo[,c(4:527)])
         
         # Derive multiple diversity metrics for each group (total/phyto/zoo)
         require("vegan")
         base$rich_total <- rowSums(as.matrix(base[,c(4:length(base))]))
         base$rich_phyto <- rowSums(as.matrix(base[,c(4:341)]))
         base$rich_zoo <- rowSums(as.matrix(base[,c(342:865)]))
         base$Htotal <- diversity(base[,c(4:865)], index = "shannon")
         base$Hphyto <- diversity(base[,c(4:341)], index = "shannon")
         base$Hzoo <- diversity(base[,c(342:865)], index = "shannon")
         
         # Same but for fut
         fut$rich_total <- rowSums(as.matrix(fut[,c(4:length(fut))]))
         fut$rich_phyto <- rowSums(as.matrix(fut[,c(4:341)]))
         fut$rich_zoo <- rowSums(as.matrix(fut[,c(342:865)]))
         fut$Htotal <- diversity(fut[,c(4:865)], index = "shannon")
         fut$Hphyto <- diversity(fut[,c(4:341)], index = "shannon")
         fut$Hzoo <- diversity(fut[,c(342:865)], index = "shannon")
         
         # And derive % changes ! head(fut[,c(1:4)]); head(base[,c(1:4)])
         message(paste("Calculating changes in diversity for ",m," || ",esm, sep = ""))
         base$fut_rich_tot <- fut[which(fut$cell_id %in% base$cell_id),c("rich_total")]
         base$fut_rich_phyto <- fut[which(fut$cell_id %in% base$cell_id),c("rich_phyto")]
         base$fut_rich_zoo <- fut[which(fut$cell_id %in% base$cell_id),c("rich_zoo")]
         base$fut_Htotal <- fut[which(fut$cell_id %in% base$cell_id),c("Htotal")]
         base$fut_Hphyto <- fut[which(fut$cell_id %in% base$cell_id),c("Hphyto")]
         base$fut_Hzoo <- fut[which(fut$cell_id %in% base$cell_id),c("Hzoo")]
      
         changes <- data.frame(id = base$cell_id, x = base$x, y = base$y, 
             diff_rich_tot = (base$fut_rich_tot)-(base$rich_total), 
             diff_rich_phyto = (base$fut_rich_phyto)-(base$rich_phyto), 
             diff_rich_zoo = (base$fut_rich_zoo)-(base$rich_zoo),
             diff_Htot = (base$fut_Htotal)-(base$Htotal),
             diff_Hphyto = (base$fut_Hphyto)-(base$Hphyto),
             diff_Hzoo = (base$fut_Hzoo)-(base$Hzoo) 
         ) # eo ddf 
         # summary(changes)
         # Add percentage changes
         changes$perc_rich_tot <- ((changes$diff_rich_tot)/base$rich_total)*100
         changes$perc_rich_phyto <- ((changes$diff_rich_phyto)/base$rich_phyto)*100
         changes$perc_rich_zoo <- ((changes$diff_rich_zoo)/base$rich_zoo)*100
         changes$perc_Htot <- ((changes$diff_Htot)/base$Htotal)*100
         changes$perc_Hphyto <- ((changes$diff_Hphyto)/base$Hphyto)*100
         changes$perc_Hzoo <- ((changes$diff_Hzoo)/base$Hzoo)*100
         
         # Mapping time with cobtours and diverging color gradients
         message(paste("Mapping annual changes in diversity (richness and H') for ",m," || ",esm, sep = ""))
         
         # Define the minimum limts for the richness diff limits of the 
         min <- floor(min(changes[,c("perc_rich_tot","perc_rich_phyto","perc_rich_zoo")], na.rm = T) )
         min2 <- floor(min(changes[,c("perc_Htot","perc_Hphyto","perc_Hzoo")], na.rm = T) )
         max2 <- abs(min2)
         # min; min2; max2
         
         map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_tot), data = changes[changes$perc_rich_tot < 50,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_tot >= 50,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_tot), 
                     data = changes[changes$perc_rich_tot < 50,] ) +
          	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_phyto), data = changes[changes$perc_rich_phyto < 50,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_phyto >= 50,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_phyto), 
                     data = changes[changes$perc_rich_phyto < 50,] ) +
          	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_zoo), data = changes[changes$perc_rich_zoo < 50,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_rich_zoo >= 50,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_zoo), 
                     data = changes[changes$perc_rich_zoo < 50,] ) +
          	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Htot), data = changes[changes$perc_Htot < max2,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_Htot >= max2,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Htot), 
                     data = changes[changes$perc_Htot < max2,] ) +
          	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hphyto), data = changes[changes$perc_Hphyto < max2,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_Hphyto >= max2,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hphyto), 
                     data = changes[changes$perc_Hphyto < max2,] ) +
          	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         #
         map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hzoo), data = changes[changes$perc_Hzoo < max2,]) +
             geom_raster(aes(x = x, y = y), data = changes[changes$perc_Hzoo >= max2,], fill = "#b2182b") +
             geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hzoo), 
                     data = changes[changes$perc_Hzoo < max2,] ) +
          	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
          	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
          	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
          	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
          		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
         
         # Save maps as .jpg for now
         setwd(paste(getwd(),"/","maps", sep = ""))
         ggsave(plot = map1, filename = paste("map_mon_perc_rich_tot_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map2, filename = paste("map_mon_perc_rich_phyto_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map3, filename = paste("map_mon_perc_rich_zoo_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map4, filename = paste("map_mon_perc_Htot_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map5, filename = paste("map_mon_perc_Hphyto_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         ggsave(plot = map6, filename = paste("map_mon_perc_Hzoo_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
         
         # Saving estimates
         setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
         message(paste("Saving diversity changes estimates for ",m," || ",esm, sep = ""))
         message(paste("", sep = ""))
         write.table(changes, file = paste("table_mon_changes_",esm,"_",m,".txt",sep=""), sep = "\t")   

     } # eo for loop - m in months

} # eo for loop - esm in ESMs


# ---------------------------------------------------------------

### 3°) Map average across ESMs (ensemble projections)

### A) Annual
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
res <- mclapply(ESMs, function(esm) {
            # Read the changes table
            d <- read.table(paste("table_annual_changes_",esm,".txt", sep = ""), sep = "\t", h = T)
            d$ESM <- esm
            return(d)
    }, mc.cores = 5
) # eo mclapply - esm
# Rbind
ens <- dplyr::bind_rows(res)
dim(ens); str(ens); summary(ens); unique(ens$ESM)
rm(res); gc()
# Compute ensemble projections
ensemble <- data.frame(ens %>%
    group_by(id) %>%
    summarize(x = unique(x), y = unique(y), perc_rich_tot = mean(perc_rich_tot,na.rm=T), 
        perc_rich_phyto = mean(perc_rich_phyto,na.rm=T), perc_rich_zoo = mean(perc_rich_zoo,na.rm=T), 
        perc_Htot = mean(perc_Htot,na.rm=T), perc_Hphyto = mean(perc_Hphyto,na.rm=T), perc_Hzoo = mean(perc_Hzoo,na.rm=T)
    )
) # eo ddf
head(ensemble); summary(ensemble)

# Define the minimum limts for the richness diff limits of the 
min <- floor(min(ensemble[,c("perc_rich_tot","perc_rich_phyto","perc_rich_zoo")], na.rm = T) )
min2 <- floor(min(ensemble[,c("perc_Htot","perc_Hphyto","perc_Hzoo")], na.rm = T) )
max2 <- abs(min2); min; min2; max2

map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_tot), data = ensemble[ensemble$perc_rich_tot < 50,]) +
    geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_rich_tot >= 50,], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_tot), 
            data = ensemble[ensemble$perc_rich_tot < 50,] ) +
 	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_phyto), data = ensemble[ensemble$perc_rich_phyto < 50,]) +
    geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_rich_phyto >= 50,], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_phyto), 
            data = ensemble[ensemble$perc_rich_phyto < 50,] ) +
 	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_zoo), data = ensemble[ensemble$perc_rich_zoo < 50,]) +
    geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_rich_zoo >= 50,], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_zoo), 
            data = ensemble[ensemble$perc_rich_zoo < 50,] ) +
 	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Htot), data = ensemble[ensemble$perc_Htot < max2,]) +
    geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_Htot >= max2,], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Htot), 
            data = ensemble[ensemble$perc_Htot < max2,] ) +
 	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hphyto), data = ensemble[ensemble$perc_Hphyto < max2,]) +
    geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_Hphyto >= max2,], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hphyto), 
            data = ensemble[ensemble$perc_Hphyto < max2,] ) +
 	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hzoo), data = ensemble[ensemble$perc_Hzoo < max2,]) +
    geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_Hzoo >= max2,], fill = "#b2182b") +
    geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hzoo), 
            data = ensemble[ensemble$perc_Hzoo < max2,] ) +
 	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

# Save maps as .jpg for now
ggsave(plot = map1, filename = paste("map_annual_perc_rich_tot_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = paste("map_annual_perc_rich_phyto_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = paste("map_annual_perc_rich_zoo_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map4, filename = paste("map_annual_perc_Htot_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map5, filename = paste("map_annual_perc_Hphyto_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)
ggsave(plot = map6, filename = paste("map_annual_perc_Hzoo_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)


### B) Monthly (12*6 maps)
for(m in months) {
    
        message(paste("Mapping ensemble diversity projections for ",m, sep = ""))
        setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/")
        ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
        res <- mclapply(ESMs, function(esm) {
                    # Read the changes table
                    d <- read.table(paste("table_mon_changes_",esm,"_",m,".txt", sep = ""), sep = "\t", h = T)
                    d$ESM <- esm
                    return(d)
            }, mc.cores = 5
        ) # eo mclapply - esm
        # Rbind
        ens <- dplyr::bind_rows(res)
        # dim(ens); str(ens); summary(ens); unique(ens$ESM)
        rm(res); gc()
        # Compute ensemble projections
        ensemble <- data.frame(ens %>%
            group_by(id) %>%
            summarize(x = unique(x), y = unique(y), perc_rich_tot = mean(perc_rich_tot,na.rm=T), 
                perc_rich_phyto = mean(perc_rich_phyto,na.rm=T), perc_rich_zoo = mean(perc_rich_zoo,na.rm=T), 
                perc_Htot = mean(perc_Htot,na.rm=T), perc_Hphyto = mean(perc_Hphyto,na.rm=T), perc_Hzoo = mean(perc_Hzoo,na.rm=T)
            )
        ) # eo ddf
        # head(ensemble); summary(ensemble)

        # Define the minimum limts for the richness diff limits of the 
        min <- floor(min(ensemble[,c("perc_rich_tot","perc_rich_phyto","perc_rich_zoo")], na.rm = T) )
        min2 <- floor(min(ensemble[,c("perc_Htot","perc_Hphyto","perc_Hzoo")], na.rm = T) )
        max2 <- abs(min2); min; min2; max2

        map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_tot), data = ensemble[ensemble$perc_rich_tot < 50,]) +
            geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_rich_tot >= 50,], fill = "#b2182b") +
            geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_tot), 
                    data = ensemble[ensemble$perc_rich_tot < 50,] ) +
         	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_phyto), data = ensemble[ensemble$perc_rich_phyto < 50,]) +
            geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_rich_phyto >= 50,], fill = "#b2182b") +
            geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_phyto), 
                    data = ensemble[ensemble$perc_rich_phyto < 50,] ) +
         	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_rich_zoo), data = ensemble[ensemble$perc_rich_zoo < 50,]) +
            geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_rich_zoo >= 50,], fill = "#b2182b") +
            geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc_rich_zoo), 
                    data = ensemble[ensemble$perc_rich_zoo < 50,] ) +
         	scale_fill_gradient2(name = "Richness difference\n(%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min,50)) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Htot), data = ensemble[ensemble$perc_Htot < max2,]) +
            geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_Htot >= max2,], fill = "#b2182b") +
            geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Htot), 
                    data = ensemble[ensemble$perc_Htot < max2,] ) +
         	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hphyto), data = ensemble[ensemble$perc_Hphyto < max2,]) +
            geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_Hphyto >= max2,], fill = "#b2182b") +
            geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hphyto), 
                    data = ensemble[ensemble$perc_Hphyto < max2,] ) +
         	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        #
        map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc_Hzoo), data = ensemble[ensemble$perc_Hzoo < max2,]) +
            geom_raster(aes(x = x, y = y), data = ensemble[ensemble$perc_Hzoo >= max2,], fill = "#b2182b") +
            geom_contour(colour = "grey60", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = perc_Hzoo), 
                    data = ensemble[ensemble$perc_Hzoo < max2,] ) +
         	scale_fill_gradient2(name = "H' difference (%)", low = "#3288bd", high = "#d53e4f", mid = "white", limits = c(min2,max2)) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

        # Save maps as .jpg for now
        setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps/")
        ggsave(plot = map1, filename = paste("map_mon_perc_rich_tot_","ensemble_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
        ggsave(plot = map2, filename = paste("map_mon_perc_rich_phyto_","ensemble_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
        ggsave(plot = map3, filename = paste("map_mon_perc_rich_zoo_","ensemble_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
        ggsave(plot = map4, filename = paste("map_mon_perc_Htot_","ensemble_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
        ggsave(plot = map5, filename = paste("map_mon_perc_Hphyto_","ensemble_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
        ggsave(plot = map6, filename = paste("map_mon_perc_Hzoo_","ensemble_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 5)
    
    
} # eo for loop - m in months 



# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------