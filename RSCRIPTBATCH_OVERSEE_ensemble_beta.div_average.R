
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

setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
# months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
thresholds <- seq(from = 0.35, to = 0.45, by = 0.01) 
# t <- 0.41
# esm <- "CESM-BEC"

res <- lapply(ESMs, function(esm) {
    
                message(paste("", sep = ""))
                message(paste("Converting HSI to 1/0 for ",esm, sep = ""))
                message(paste("", sep = ""))
                if(esm %in% c("CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")) {
                        detach("package:dplyr", unload = T)
                } # eo if loop
                
                annual <- lapply(thresholds, function(t) {
                    
                            # Get months names based on couples
                            message(paste("Using t = ",t, sep = ""))
                            # Load the baseline composition for these two months, for both zoo and phyto
                            base <- read.table(paste("table_annual_mean_compo_plankton_baseline.txt", sep = ""), sep = "\t")
                            # colnames(base)
                            base <- base[,c(1:(length(base)-2))]
                            phyto.fut <- read.table(paste("table_annual_mean_compo_phyto_2100-2000_rcp85_",esm,".txt", sep = ""), sep = "\t")
                            zoo.fut <- read.table(paste("table_annual_mean_compo_zoo_2100-2000_rcp85_",esm,".txt", sep = ""), sep = "\t")
                            # dim(base); dim(phyto.fut); dim(zoo.fut)
                            # Cbind to have all plankton
                            fut <- cbind(phyto.fut, zoo.fut[,c(4:length(zoo.fut))])
                            rm(phyto.fut,zoo.fut);gc()
                            # Remove points in the colnames
                            colnames(base)[c(4:865)] <- gsub("[.]","",as.character(colnames(base)[c(4:865)])) 
                            colnames(fut)[c(4:865)] <- gsub("[.]","",as.character(colnames(fut)[c(4:865)])) 
                            base <- na.omit(base)
                            fut <- na.omit(fut)
                            # dim(base); dim(fut)
                            # Compute richness based on HSI
                            base$rich_plankton <- rowSums(as.matrix(base[,c(4:865)]))
                            base$rich_phyto <- rowSums(as.matrix(base[,c(4:341)]))
                            base$rich_zoo <- rowSums(as.matrix(base[,c(342:865)]))
                            fut$rich_plankton <- rowSums(as.matrix(fut[,c(4:865)]))
                            fut$rich_phyto <- rowSums(as.matrix(fut[,c(4:341)]))
                            fut$rich_zoo <- rowSums(as.matrix(fut[,c(342:865)]))
                            # Uneven because monthly climatologies -> make them even
                            base <- base[order(base$cell_id),]
                            fut <- fut[order(fut$cell_id),]
                            base <- base[which(base$cell_id %in% unique(fut$cell_id)),]
                            fut <- fut[which(fut$cell_id %in% unique(base$cell_id)),]
                            # summary(base$rich_plankton); summary(fut$rich_plankton)
                            
                            # Convert to P/1 using rbinom and initial probability
                            for(sp in colnames(base)[c(4:865)] ) {
                                    message(paste("Converting probabilities for ",sp, sep = ""))
                                    base[,c(sp)][base[,c(sp)] > t] <- 1
                                    base[,c(sp)][base[,c(sp)] <= t] <- 0
                                    # Future 2 now
                                    fut[,c(sp)][fut[,c(sp)] > t] <- 1
                                    fut[,c(sp)][fut[,c(sp)] <= t] <- 0
                            } # eo for loop
                            # head(base); head(fut)
                    
                            # And compute beta.div changes 
                            M = nrow(base[,c(4:865)])
                            N = ncol(base[,c(4:865)])
                            
                            ### With paralelling : 
                            require("doParallel")
                            require("plyr")
                            registerDoParallel(cores = 25)
                            # Compute beta.div changes ofr all plankton
                            div <- data.frame()
                            d <- cbind(base[,c(4:865)], fut[,c(4:865)]) 
                            d$bit <- cut(1:M, 50, labels = F)
                            div <- ddply(d, ~ bit, function(x) {
                            	 	    beta.temp(x[,1:N], x[,(N+1):(2*N)],"jaccard")
                                    },.parallel = T
                            ) # eo ddply
                            # summary(div)
            
                            # For phytoplankton
                            M = nrow(base[,c(4:341)])
                            N = ncol(base[,c(4:341)])
                            div.phyto <- data.frame()
                            d <- cbind(base[,c(4:341)], fut[,c(4:341)]) 
                            d$bit <- cut(1:M, 50, labels = F)
                            # head(d)
                            div.phyto <- ddply(d, ~ bit, function(x) {
                            	 	    beta.temp(x[,1:N], x[,(N+1):(2*N)],"jaccard")
                                    },.parallel = T
                            ) # eo ddply
                            
                            # For zooplankton
                            M = nrow(base[,c(342:865)])
                            N = ncol(base[,c(342:865)])
                            div.zoo <- data.frame()
                            d <- cbind(base[,c(342:865)], fut[,c(342:865)]) 
                            d$bit <- cut(1:M, 50, labels = F)
                            div.zoo <- ddply(d, ~ bit, function(x) {
                            	 	    beta.temp(x[,1:N], x[,(N+1):(2*N)],"jaccard")
                                    },.parallel = T
                            ) # eo ddply
                            
                            if( sum(is.na(div$beta.jtu)) >= 1 ) {
                                # summary(div.phyto)
                                div[which(is.na(div$beta.jac)),c("beta.jac","beta.jne","beta.jtu")] <- 0
                                div[which(is.na(div$beta.jne) & !is.na(div$beta.jac)),c("beta.jne")] <- 1
                                div[which(is.na(div$beta.jtu) & !is.na(div$beta.jac)),c("beta.jtu")] <- 1
                            } 
                            
                            if( sum(is.na(div.phyto$beta.jtu)) >= 1 ) {
                                # summary(div.phyto)
                                div.phyto[which(is.na(div.phyto$beta.jac)),c("beta.jac","beta.jne","beta.jtu")] <- 0
                                div.phyto[which(is.na(div.phyto$beta.jne) & !is.na(div.phyto$beta.jac)),c("beta.jne")] <- 1
                                div.phyto[which(is.na(div.phyto$beta.jtu) & !is.na(div.phyto$beta.jac)),c("beta.jtu")] <- 1
                            } 
                            
                            if( sum(is.na(div.zoo$beta.jtu)) >= 1 ) {
                                # summary(div.phyto)
                                div.zoo[which(is.na(div.zoo$beta.jac)),c("beta.jac","beta.jne","beta.jtu")] <- 0
                                div.zoo[which(is.na(div.zoo$beta.jne) & !is.na(div.zoo$beta.jac)),c("beta.jne")] <- 1
                                div.zoo[which(is.na(div.zoo$beta.jtu) & !is.na(div.zoo$beta.jac)),c("beta.jtu")] <- 1
                            } 
                            
                            # Gather all in a data.frame and return
                            div$x <- base$x
                            div$y <- base$y
                            div$id <- factor(paste(div$x, div$y, sep = "_"))
                            div$rich_all_base <- base$rich_plankton
                            div$rich_phyto_base <- base$rich_phyto
                            div$rich_zoo_base <- base$rich_zoo
                            div$rich_all_fut <- fut$rich_plankton
                            div$rich_phyto_fut <- fut$rich_phyto
                            div$rich_zoo_fut <- fut$rich_zoo
                            # Compute differences
                            div$diff_all <- (div$rich_all_fut) - (div$rich_all_base)
                            div$diff_phyto <- (div$rich_phyto_fut) - (div$rich_phyto_base)
                            div$diff_zoo <- (div$rich_zoo_fut) - (div$rich_zoo_base)
                            # And add zoo and phyto components of beta.div
                            div$jac.phyto <- div.phyto$beta.jac
                            div$jne.phyto <- div.phyto$beta.jne
                            div$jtu.phyto <- div.phyto$beta.jtu
            
                            div$jac.zoo <- div.zoo$beta.jac
                            div$jne.zoo <- div.zoo$beta.jne
                            div$jtu.zoo <- div.zoo$beta.jtu

                            # div[is.na(div$jac.phyto),c("x","y","jac.phyto")]

                            # Add the couple
                            div$thresh <- factor(t)
                            # Return(div)
                            rm(d,base.m1,base.m2); gc()
                            return(div)
                    
                    } # eo 2nd FUN
                    
                ) # eo 2nd lapply - annual
                # Rbind
                detach("package:plyr", unload = T)
                require("dplyr")
                table <- bind_rows(annual) 
                
                # Compute average changes across thresholds
                div <- data.frame(table %>% 
                        group_by(id) %>% 
                        summarise(x = unique(x), y = unique(y),
                        rich_all_base = mean(rich_all_base,na.rm=T), rich_phyto_base = mean(rich_phyto_base,na.rm=T), rich_zoo_base = mean(rich_zoo_base,na.rm=T),
                        rich_all_fut = mean(rich_all_fut,na.rm=T), rich_phyto_fut = mean(rich_phyto_fut,na.rm=T), rich_zoo_fut = mean(rich_zoo_fut,na.rm=T),
                        jac = mean(beta.jac,na.rm=T), jne = mean(beta.jne,na.rm=T), jtu = mean(beta.jtu,na.rm=T), 
                        jac_phyto = mean(jac.phyto,na.rm=T), jne_phyto = mean(jne.phyto,na.rm=T), jtu_phyto = mean(jtu.phyto,na.rm=T), 
                        jac_zoo = mean(jac.zoo,na.rm=T), jne_zoo = mean(jne.zoo,na.rm=T), jtu_zoo = mean(jtu.zoo,na.rm=T), 
                        diff_all = mean(diff_all,na.rm=T), diff_phyto = mean(diff_phyto,na.rm=T), diff_zoo = mean(diff_zoo,na.rm=T)
                        ) # eo summarise
                ) # eo ddf
                # summary(div)
                div$ESM <- esm
                # Return div
                rm(table); gc()
                return(div)
    
    } # eo 1st FUN
    
) # eo 1st lapply - thresholds
# Rbind 
require("dplyr")
ddf <- bind_rows(res) 
dim(ddf)
head(ddf); unique(ddf$ESM)
summary(ddf)

# Save
write.table(ddf, file = "table_annual_comm_beta.div_ensemble.txt", sep = "\t")
dim(ddf[ddf$ESM == "GFDL-TOPAZ",]); dim(ddf[ddf$ESM == "IPSL-PISCES",]); dim(ddf[ddf$ESM == "MRI-NEMURO",])

### Map changes in beta.div for each ESM
for(esm in ESMs) { 
    
        # Mapping time
        setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps")
        message(paste("Mapping annual beta div changes for ",esm, sep = ""))
        d <- ddf[ddf$ESM == esm,]
        
        require("viridis")
        # Annual plankton jac
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = d) +
        	scale_fill_viridis(name = "Annual Jaccard", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = d) +
        	scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = d) +
        	scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = d) +
        	scale_fill_viridis(name = "Annual Jaccard index", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = d) +
        	scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = d) +
        	scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = d) +
        	scale_fill_viridis(name = "Annual Jaccard index", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = d) +
        	scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = d) +
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
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = d) +
        	scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
            geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = d) +
         	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
         	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                       labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
         	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
         		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
           	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        # 
        ggsave(plot = map, filename = paste("map_annual_jtu_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

        setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")
        
} # 

### And map ensemble projections
ens <- data.frame(ddf %>% 
        group_by(id) %>% 
        summarize(x = unique(x), y = unique(y), 
            jac = mean(jac,na.rm=T), jne = mean(jne,na.rm=T), jtu = mean(jtu,na.rm=T), 
            jac_phyto = mean(jac_phyto,na.rm=T), jne_phyto = mean(jne_phyto,na.rm=T), jtu_phyto = mean(jtu_phyto,na.rm=T),
            jac_zoo = mean(jac_zoo,na.rm=T), jne_zoo = mean(jne_zoo,na.rm=T), jtu_zoo = mean(jtu_zoo,na.rm=T)
        ) 
) # eo ddf - ens
summary(ens)


setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/maps")
# Annual plankton jac
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ens) +
	scale_fill_viridis(name = "Annual Jaccard", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jac_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

# Annual plankton nestedness
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = ens) +
	scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jne_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

# Annual plankton turn-over
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ens) +
	scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jtu_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

### Phyto
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_phyto), data = ens) +
	scale_fill_viridis(name = "Annual Jaccard index", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jac_phyto_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

# Annual Phyto nestedness
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_phyto), data = ens) +
	scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jne_phyto_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

# Annual Phyto turn-over
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_phyto), data = ens) +
	scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_phyto), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jtu_phyto_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	


### Zoo
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac_zoo), data = ens) +
	scale_fill_viridis(name = "Annual Jaccard index", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jac_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jac_zoo_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)		

# Annual Zoo nestedness
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne_zoo), data = ens) +
	scale_fill_viridis(name = "Annual Nestedness", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jne_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jne_zoo_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	

# Annual Zoo turn-over
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu_zoo), data = ens) +
	scale_fill_viridis(name = "Annual Turn-over", limits = c(0,1)) +
    geom_contour(colour = "grey75", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = jtu_zoo), data = ens) +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# 
ggsave(plot = map, filename = paste("map_annual_jtu_zoo_","ensemble",".jpg", sep = ""), dpi = 300, width = 7, height = 5)	


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

