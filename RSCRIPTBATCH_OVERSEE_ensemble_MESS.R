
# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("reshape2")
library("tidyverse")
library("viridis")
library("RColorBrewer")
library("modEvA")
library("maps")
library("parallel")

world2 <- map_data("world2")

firstup <- function(x) {
  substr(x,1,1) <- toupper(substr(x,1,1))
  x
} # FUN to make first letter capital

WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------


# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
# vector of ESMs
ESM <- c("CESM-BEC","GFDL-TOPAZ","IPSL-PISCES","CNRM-PISCES","MRI-NEMURO")

# m <- "nov"
# p <- "p1"
# esm <- "IPSL-PISCES"

for(esm in ESM) {
    
        message(paste("Computing MESS for ",esm, sep = ""))
        
        for(p in pools) {
            
            message(paste("     based on pool ",p, sep = ""))
        	
            # Set pool of variables (zooplankton)
        	if(p == "p1") {
        		vars <- c("SST","dSST","dO2","logNO3","logChl")
        	} else if (p == "p2") {
        		vars <- c("SST","dSST","dO2","logSiO2","logChl")
        	} else if (p == "p3") {
        		vars <- c("SST","dSST","dO2","logSiO2","logChl","Nstar")
        	} else if (p == "p4") {
        		vars <- c("SST","dSST","dO2","logNO3","logChl","Sistar")
        	} # eo if else loop
            
            for(m in months) {
                
                    message(paste("             for month = ",m, sep = ""))
                    message(paste("", sep = ""))
                    
                    # get baseline clim
            		setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
            		envT1 <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
                    colnames(envT1)[7] <- "dSST"
                    # dim(envT1); summary(envT1)
                    # Exclude SSS < 25 & bathy > 175
                    envT1 <- envT1[-which(envT1$SSS < 20),]
                    envT1 <- envT1[-which(envT1$Bathy > -175),]
                    envT1$x2 <- envT1$x 
                    envT1[envT1$x < 0 ,"x2"] <- (envT1[envT1$x < 0 ,"x"]) + 360
                    envT1$id <- paste(envT1$x2, envT1$y, sep = "_")
                    envT1 <- envT1[order(envT1$id),]
                    
            		# get future monthly clim
                    mm <- firstup(m) # future clims need a capital first letter

                    setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims/",esm,"/", sep = ""))
            		envT2 <- read.table(paste("clims_mon_",mm,"_",esm,"_rcp85_base+2100-2031.txt", sep = ""), h = T, sep = "\t")
                    # colnames(envT1); colnames(envT2)
                    # dim(envT2); summary(envT2)

            		# Perform MESS
            		V <- na.omit(envT1[,vars])
            		P <- na.omit(envT2[,vars])
                    # summary(V)
                    # summary(P)
            		# Need to split the "P matrix" (projection data) into sub groups of cells and then apply MESS wiht mclapply
            		P$bit <- cut(1:nrow(P), 20, labels = F)
            		# b <- 1 # for testing
            		require("parallel")
            		res <- mclapply(unique(P$bit), function(b) {
            					message(paste("Doing bit || ", b, sep = ""))
            					mess <- modEvA::MESS(V = V, P = P[P$bit == b,vars])
            					return(mess)	
            				}, mc.cores = 30
            		) # eo lapply
                    require("dplyr")
            		mess <- bind_rows(res)
            		rm(res,V,P)
		
            		# Provide coords to mess
            		mess$x <- na.omit(envT2[,c("x","y",vars)])[,"x"]
            		mess$y <- na.omit(envT2[,c("x","y",vars)])[,"y"]
            		mess <- mess[,c("x","y","TOTAL","MoD")]
                    # summary(mess)
		
            		# Go to proper dir() and plot MESS map and distrbution of MoD < 0 (if any)
            		setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/MESS_tables_31_03_20/plots/", sep = ""))
                    # Map MESS values
                    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = TOTAL), data = mess) +
                     	scale_fill_gradient2(name = "MESS", low = "#c51b7d", high = "#7fbc41", mid = "white") +
                        geom_contour(colour = "grey40", binwidth = 10, size = 0.25, aes(x = x, y = y, z = TOTAL), data = mess) +
                     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
                    # Plot distrib of extrap variables
                    ggsave(plot = map, filename = paste("map_MESS_zoo_env_2100-2000_rcp85_",esm,"_",p,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 4)
                        
                    plot <- ggplot(data = mess[which(mess$TOTAL < 0),], aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
                    			scale_fill_brewer(name = "Variable", palette = "Paired") + 
                    			geom_density(position = "stack") + ylab("Density") + xlab("") + 
                    			theme_classic()
                                
                    ggsave(plot = plot, filename = paste("distrib_vars_MESS_zoo_2100-2000_rcp85_",esm,"_",p,"_",m,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
            		
                    ### Save the data in specific directory
                    setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/MESS_tables_31_03_20", sep = ""))
                    save(mess, file = paste("table_MESS_zoo_2100-2000_rcp85_",esm,"_",p,"_",m,".Rdata", sep = ""))
            		rm(mess,V,P,envT1,envT2); gc()
                    setwd(WD)

            } # eo m in months
            
        } # eo for loop - p in pools    
    
} # eo for loop - esm in ESM


### Same, but for phytoplankton vars !

for(esm in ESM) {
    
        message(paste("Computing MESS for ",esm, sep = ""))
        
        for(p in pools) {
            
            message(paste("     based on pool ",p, sep = ""))
        	
            # Set pool of variables (zooplankton)
        	if(p == "p1") {
        		vars <- c("SST","dSST","PAR","logNO3","logChl","Nstar")
        	} else if (p == "p2") {
        		vars <- c("SST","dSST","PAR","logSiO2","logChl","Nstar")
        	} else if (p == "p3") {
        		vars <- c("SST","dSST","PAR","logNO3","logChl","Nstar","Sistar")
        	} else if (p == "p4") {
        		vars <- c("SST","dSST","logChl","PAR","logNO3")
        	} # eo else if loop
            
            for(m in months) {
                
                    message(paste("             for month = ",m, sep = ""))
                    message(paste("", sep = ""))
                    
                    # get baseline clim
            		setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
            		envT1 <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
                    colnames(envT1)[7] <- "dSST"
                    # dim(envT1)
                    # Exclude SSS < 25 & bathy > 175
                    envT1 <- envT1[-which(envT1$SSS < 20),]
                    envT1 <- envT1[-which(envT1$Bathy > -175),]
                    envT1$x2 <- envT1$x 
                    envT1[envT1$x < 0 ,"x2"] <- (envT1[envT1$x < 0 ,"x"]) + 360
                    envT1$id <- paste(envT1$x2, envT1$y, sep = "_")
                    envT1 <- envT1[order(envT1$id),]
                    
            		# get future monthly clim
                    mm <- firstup(m) # future clims need a capital first letter

                    setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims/",esm,"/", sep = ""))
                    # clims_mon_Dec_CESM-BEC_rcp85_base+2100-2031.txt
            		envT2 <- read.table(paste("clims_mon_",mm,"_",esm,"_rcp85_base+2100-2031.txt", sep = ""), h = T, sep = "\t")
                    # colnames(envT1); colnames(envT2)

            		# Perform MESS
            		V <- na.omit(envT1[,vars])
            		P <- na.omit(envT2[,vars])
                  	# Need to split the "P matrix" (projection data) into sub groups of cells and then apply MESS wiht mclapply
            		P$bit <- cut(1:nrow(P), 20, labels = F)
            		# b <- 1 # for testing
            		require("parallel")
            		res <- mclapply(unique(P$bit), function(b) {
            					message(paste("Doing bit || ", b, sep = ""))
            					mess <- modEvA::MESS(V = V, P = P[P$bit == b,vars])
            					return(mess)	
            				}, mc.cores = 30
            		) # eo lapply
                    require("dplyr")
            		mess <- bind_rows(res)
            		rm(res,V,P)
		
            		# Provide coords to mess
            		mess$x <- na.omit(envT2[,c("x","y",vars)])[,"x"]
            		mess$y <- na.omit(envT2[,c("x","y",vars)])[,"y"]
            		mess <- mess[,c("x","y","TOTAL","MoD")]
                    # summary(mess)
		
            		# Go to proper dir() and plot MESS map and distrbution of MoD < 0 (if any)
            		setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/MESS_tables_31_03_20/plots/", sep = ""))
                    # Map MESS values
                    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = TOTAL), data = mess) +
                     	scale_fill_gradient2(name = "MESS", low = "#c51b7d", high = "#7fbc41", mid = "white") +
                        geom_contour(colour = "grey40", binwidth = 10, size = 0.25, aes(x = x, y = y, z = TOTAL), data = mess) +
                     	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                     	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                                   labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                     	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                     		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                    # 
                    ggsave(plot = map, filename = paste("map_MESS_phyto_2100-2000_rcp85_",esm,"_",p,"_",m,".jpg", sep = ""), dpi = 300, width = 7, height = 4)
                        
                    plot <- ggplot(data = mess[which(mess$TOTAL < 0),], aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
                    			scale_fill_brewer(name = "Variable", palette = "Paired") + 
                    			geom_density(position = "stack") + ylab("Density") + xlab("") + 
                    			theme_classic()
                                
                    ggsave(plot = plot, filename = paste("distrib_vars_MESS_phyto_2100-2000_rcp85_",esm,"_",p,"_",m,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
        
            		### Save the data in specific directory
                    setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/MESS_tables_31_03_20", sep = ""))
                    save(mess, file = paste("table_MESS_phyto_2100-2000_rcp85_",esm,"_",p,"_",m,".Rdata", sep = ""))
                    rm(mess,V,P,envT1,envT2); gc()
                    setwd(WD)

            } # eo m in months
            
        } # eo for loop - p in pools    
    
} # eo for loop - esm in ESM


# --------------------------------------------------

### Get all results from the code above and compyte mean ensemble MESS values and (map).
setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/MESS_tables_31_03_20", sep = ""))
files.phyto <- dir()[grep("_phyto_",dir())]
files.zoo <- dir()[grep("_zoo_",dir())]

### A) Get all phyto results
# f <- files.phyto[1]
res.phyto <- mclapply(files.phyto, function(f) {
            # Useless message
            message(paste("Loading ",f, sep = ""))
            d <- get(load(f))
            # get terms
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            esm <- terms[6,1]
            p <- terms[7,1]
            month <- terms[8,1]
            month <- str_replace(as.character(month),".Rdata","")
            d$ESM <- esm
            d$pool <- p
            d$month <- month
            # Return
            return(d)
    }, mc.cores = 30
) # eo mclapply - f in files.phyto
# Rbind
table <- bind_rows(res.phyto)
dim(table); summary(table)
rm(res.phyto); gc()

# First examine overall MoD
ggplot(data = table[which(table$TOTAL < 0),], aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
    scale_fill_brewer(name = "Variable", palette = "Paired") + 
    geom_density(position = "stack") + ylab("Density") + xlab("") + 
    theme_classic() + facet_wrap(~factor(month), ncol = 4, scales = "free")
    
# and face per ESM                         
ggplot(data = table[which(table$TOTAL < 0),], aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
    scale_fill_brewer(name = "Variable", palette = "Paired") + 
    geom_density(position = "stack") + ylab("Density") + xlab("") + 
    theme_classic() + facet_wrap(~factor(ESM), ncol = 3, scales = "free")

# Variability across months seems higher than variations across ESMs
setwd(paste("/net/kryo/work/fabioben/OVERSEE/data", sep = ""))
for(esm in unique(table$ESM)) {
    
        plot <- ggplot(data = table[which(table$TOTAL < 0 & table$ESM == esm),],
                        aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
            scale_fill_brewer(name = "Variable", palette = "Paired") + 
            geom_density(position = "stack") + ylab("Density") + xlab("") + 
            theme_classic() + facet_wrap(~factor(month), ncol = 4, scales = "free")
        # save 
        ggsave(plot = plot, filename = paste("plot_distrib_MoD_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 17, height = 7)
        
} # eo for loop - esm in unique(table$ESM)
    
### Compute ensemble 
table$id <- paste(table$x, table$y, sep = "_")
ens.phy <- data.frame(table %>% group_by(id, ESM) %>% summarize(x = unique(x), y = unique(y), mess = mean(TOTAL,na.rm=T)) ) # eo ddf
dim(ens.phy); summary(ens.phy)

for(esm in unique(ens.phy$ESM)) {
    
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mess), data = ens.phy[ens.phy$ESM == esm,]) +
            scale_fill_gradient2(name = "MESS", low = "#c51b7d", high = "#7fbc41", mid = "white") +
            geom_contour(colour = "grey40", binwidth = 25, size = 0.25, aes(x = x, y = y, z = mess), data = ens.phy[ens.phy$ESM == esm,]) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
        # save map
        ggsave(plot = map, filename = paste("map_ann_MESS_phyto_",esm,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
        
} # eo for loop - esm in unique(table$ESM)

### And per month?
ens2.phy <- data.frame(table %>% group_by(id,ESM,month) %>% summarize(x = unique(x), y = unique(y), mess = mean(TOTAL,na.rm=T)) ) # eo ddf
summary(ens2.phy)

for(esm in unique(ens2.phy$ESM)) {
    
    for(m in unique(ens2.phy$month)) {
    
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mess), data = ens2.phy[ens2.phy$ESM == esm & ens2.phy$month == m,]) +
            scale_fill_gradient2(name = "MESS", low = "#c51b7d", high = "#7fbc41", mid = "white") +
            geom_contour(colour = "grey40", binwidth = 20, size = 0.25, aes(x = x, y = y, z = mess),
                    data = ens2.phy[ens2.phy$ESM == esm & ens2.phy$month == m,]) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
        # save map
        ggsave(plot = map, filename = paste("map_mon_MESS_phyto_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
    
    } # eo for loop - m in months
    
} # eo for loop - esm in unique(table$ESM)



### B) Get zoo results
setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/MESS_tables_31_03_20", sep = ""))
res.zoo <- mclapply(files.zoo, function(f) {
            # Useless message
            message(paste("Loading ",f, sep = ""))
            d <- get(load(f))
            # get terms
            terms <- do.call(cbind, strsplit(as.character(f),"_"))
            esm <- terms[6,1]
            p <- terms[7,1]
            month <- terms[8,1]
            month <- str_replace(as.character(month),".Rdata","")
            d$ESM <- esm
            d$pool <- p
            d$month <- month
            # Return
            return(d)
    }, mc.cores = 25
) # eo mclapply - f in files.phyto
# Rbind
table <- bind_rows(res.zoo)
dim(table); summary(table)
rm(res.zoo); gc()

setwd(paste("/net/kryo/work/fabioben/OVERSEE/data", sep = ""))

# Variability across months seems higher than variations across ESMs
setwd(paste("/net/kryo/work/fabioben/OVERSEE/data", sep = ""))
for(esm in unique(table$ESM)) {
    
        plot <- ggplot(data = table[which(table$TOTAL < 0 & table$ESM == esm),],
                        aes(factor(MoD), stat(count), fill = factor(MoD)) ) + 
            scale_fill_brewer(name = "Variable", palette = "Paired") + 
            geom_density(position = "stack") + ylab("Density") + xlab("") + 
            theme_classic() + facet_wrap(~factor(month), ncol = 4, scales = "free")
        # save 
        ggsave(plot = plot, filename = paste("plot_distrib_MoD_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 17, height = 7)
        
} # eo for loop - esm in unique(table$ESM)
    
### Compute ensemble 
table$id <- paste(table$x, table$y, sep = "_")
ens.zoo <- data.frame(table %>% group_by(id, ESM) %>% summarize(x = unique(x), y = unique(y), mess = mean(TOTAL,na.rm=T)) ) # eo ddf
summary(ens.zoo) ; dim(ens.zoo)

for(esm in unique(ens.zoo$ESM)) {
    
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mess), data = ens.zoo[ens.zoo$ESM == esm,]) +
            scale_fill_gradient2(name = "MESS", low = "#c51b7d", high = "#7fbc41", mid = "white") +
            geom_contour(colour = "grey40", binwidth = 20, size = 0.25, aes(x = x, y = y, z = mess), data = ens.zoo[ens.zoo$ESM == esm,]) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
        # save map
        ggsave(plot = map, filename = paste("map_ann_MESS_zoo_",esm,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
        
} # eo for loop - esm in unique(table$ESM)

### And per month?
ens2.zoo <- data.frame(table %>% group_by(id,ESM,month) %>% summarize(x = unique(x), y = unique(y), mess = mean(TOTAL,na.rm=T)) ) # eo ddf
summary(ens2.zoo)

for(esm in unique(ens2.zoo$ESM)) {
    
    for(m in unique(ens2.zoo$month)) {
    
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mess), data = ens2.zoo[ens2.zoo$ESM == esm & ens2.zoo$month == m,]) +
            scale_fill_gradient2(name = "MESS", low = "#c51b7d", high = "#7fbc41", mid = "white") +
            geom_contour(colour = "grey40", binwidth = 20, size = 0.25, aes(x = x, y = y, z = mess),
                    data = ens2.zoo[ens2.zoo$ESM == esm & ens2.zoo$month == m,]) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                    labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
            scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                    labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
        
        # save map
        ggsave(plot = map, filename = paste("map_mon_MESS_zoo_",esm,"_",m,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
    
    } # eo for loop - m in months
    
} # eo for loop - esm in unique(table$ESM)


### 29/04/2020: Make a 10 maps (trophic level*5 ESMs) panel to show mean annual negative MESS values 
unique(ens.phy$ESM)
unique(ens.zoo$ESM)
# summary(ens.zoo$mess) ; summary(ens.phy$mess)

### To test the discrete colorscale
ens.phy$discr <- NA
ens.phy[which(ens.phy$mess < 0 & ens.phy$mess >= -2.5),"discr"] <- "0-2.5"
ens.phy[which(ens.phy$mess < -2.5 & ens.phy$mess >= -5),"discr"] <- "2.5-5"
ens.phy[which(ens.phy$mess < -5 & ens.phy$mess >= -10),"discr"] <- "5-10"
ens.phy[which(ens.phy$mess < -10 & ens.phy$mess >= -15),"discr"] <- "10-15"
ens.phy[which(ens.phy$mess < -15 & ens.phy$mess >= -30),"discr"] <- "15-30"
ens.phy[which(ens.phy$mess < -30),"discr"] <- ">30"
ens.zoo$discr <- NA
ens.zoo[which(ens.zoo$mess < 0 & ens.zoo$mess >= -2.5),"discr"] <- "0-2.5"
ens.zoo[which(ens.zoo$mess < -2.5 & ens.zoo$mess >= -5),"discr"] <- "2.5-5"
ens.zoo[which(ens.zoo$mess < -5 & ens.zoo$mess >= -10),"discr"] <- "5-10"
ens.zoo[which(ens.zoo$mess < -10 & ens.zoo$mess >= -15),"discr"] <- "10-15"
ens.zoo[which(ens.zoo$mess < -15 & ens.zoo$mess >= -30),"discr"] <- "15-30"
ens.zoo[which(ens.zoo$mess < -30),"discr"] <- ">30"

cls <- c("0-2.5"="#fed976","2.5-5"="#feb24c","5-10"="#fd8d3c","10-15"="#fc4e2a","15-30"="#e31a1c",">30"="#bd0026")
ord <- c("0-2.5","2.5-5","5-10","10-15","15-30",">30")

# summary(ens.phy[ens.phy$ESM == "CESM-BEC" & ens.phy$mess < 0,])
mapA <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.phy[ens.phy$ESM == "CESM-BEC" & ens.phy$mess < 0,]) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "CESM-BEC" & ens.phy$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(CESM-BEC)", palette = "YlOrRd", direction = 1) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapB <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.zoo[ens.zoo$ESM == "CESM-BEC" & ens.zoo$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(CESM-BEC)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "CESM-BEC" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapC <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.phy[ens.phy$ESM == "CNRM-PISCES" & ens.phy$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(CNRM-PISCES)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "CNRM-PISCES" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapD <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.zoo[ens.zoo$ESM == "CNRM-PISCES" & ens.zoo$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(CNRM-PISCES)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "CNRM-PISCES" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapE <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.phy[ens.phy$ESM == "GFDL-TOPAZ" & ens.phy$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(GFDL-TOPAZ)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "GFDL-TOPAZ" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapF <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.zoo[ens.zoo$ESM == "GFDL-TOPAZ" & ens.zoo$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(GFDL-TOPAZ)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "GFDL-TOPAZ" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapG <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.phy[ens.phy$ESM == "IPSL-PISCES" & ens.phy$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(IPSL-PISCES)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "IPSL-PISCES" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapH <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.zoo[ens.zoo$ESM == "IPSL-PISCES" & ens.zoo$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(IPSL-PISCES)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "IPSL-PISCES" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapI <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.phy[ens.phy$ESM == "MRI-NEMURO" & ens.phy$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(MRI-NEMURO)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "MRI-NEMURO" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
mapJ <- ggplot() + #geom_tile(aes(x = x, y = y, fill = abs(mess)), data = ens.zoo[ens.zoo$ESM == "MRI-NEMURO" & ens.zoo$mess < 0,]) +
    #scale_fill_distiller(name = "Mean annual\ndissimilarity\n(MRI-NEMURO)", palette = "YlOrRd", direction = 1) +
    geom_tile(aes(x = x, y = y, fill = factor(discr)), data = ens.phy[ens.phy$ESM == "MRI-NEMURO" & ens.phy$mess < 0,]) +
    scale_fill_manual(name = "Mean annual\nenvironmental\ndissimilarity", values = cls, breaks = ord) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey95", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
            labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
            labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "grey60"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
            
### Combine in panel
library("ggpubr")
ggarrange(mapA, mapB, mapC, mapC, mapD, mapE, mapF, mapG, mapH, mapI, mapJ, ncol = 2, nrow = 5, labels = LETTERS[1:10])
panel <- ggarrange(mapA, mapB, mapC, mapD, mapE, mapF, mapG, mapH, mapI, mapJ, ncol = 2, nrow = 5, labels = LETTERS[1:10])
ggsave(plot = panel, filename = "panel_ann_MESS_ESMs_discrete.jpg", dpi = 300, height = 13, width = 11)


### 07/08/2020: Map regions of mean MESS < 0 for phyto, zoo and tot
dim(ens.phy) ; dim(ens.zoo)
ens.phy2 <- data.frame(ens.phy %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), mess = mean(mess, na.rm=T)) ) # eo ddf
ens.zoo2 <- data.frame(ens.zoo %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), mess = mean(mess, na.rm=T)) ) # eo ddf
dim(ens.phy2) ; dim(ens.zoo2)

map_mean_annual_mess_phyto <- ggplot() + geom_tile(aes(x = x, y = y), data = ens.phy2[ens.phy2$mess < 0,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_mean_annual_mess_zoo <- ggplot() + geom_tile(aes(x = x, y = y), data = ens.zoo2[ens.zoo2$mess < 0,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
map_mean_annual_mess_plankton <- ggplot() +
    geom_tile(aes(x = x, y = y), data = ens.phy2[ens.phy2$mess < 0,], fill = "red") +
    geom_tile(aes(x = x, y = y), data = ens.zoo2[ens.zoo2$mess < 0,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
# Save as .eps
grDevices::cairo_ps(filename = "map_mean_annual_mess_phyto.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_mean_annual_mess_phyto)
dev.off()
grDevices::cairo_ps(filename = "map_mean_annual_mess_zoo.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_mean_annual_mess_zoo)
dev.off()
grDevices::cairo_ps(filename = "map_mean_annual_mess_plankton.ps", width = 7, height = 5, fallback_resolution = 300)
print(map_mean_annual_mess_plankton)
dev.off()


### test .ps or geom_tile instead

test.map1 <- ggplot() + geom_raster(aes(x = x, y = y), data = ens.phy2[ens.phy2$mess < 0,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
#
test.map2 <- ggplot() + geom_tile(aes(x = x, y = y), data = ens.phy2[ens.phy2$mess < 0,], fill = "red") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
               labels = c("180°W","120°W","60°W","0°E","60°E","120°E","180°E") ) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		      labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N") ) + theme_void()
              
grDevices::cairo_ps(filename = "test_map1.ps", width = 7, height = 5, fallback_resolution = 300)
print(test.map1)
dev.off()
              
grDevices::cairo_ps(filename = "test_map2.ps", width = 7, height = 5, fallback_resolution = 300)
print(test.map2)
dev.off()
              

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------