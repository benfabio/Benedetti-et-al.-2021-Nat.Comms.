
##### 19/11/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, D-USYS, ETH Zürich
##### Script for : 
#		- extracting the mon_diff climatologies you built from the CMIP5 netCDF of Charlotte
#		- apply difference change to in situ climatologies to create the future fields
#		- variables to compute from: SST, dSST, logChl, logNO3, logSiO2, Si*, N*, dO2, PAR

### Last update: 15/01/2021

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("ncdf4")
library("raster")
library("sp")
library("reshape2")
library("scales")
library("maps")
library("cmocean")
library("RColorBrewer")
library("viridis")

# Coastline
world2 <- map_data("world2")

CapStr <- function(y) {
          c <- strsplit(y, " ")[[1]]
          paste(toupper(substring(c, 1,1)), substring(c, 2), sep = "", collapse = " ")
} # eo fun
capitalize_str <- function(charcter_string) {
      sapply(charcter_string, CapStr)
}

# --------------------------------------------------------------------------------------------------------------------------------

### In a singlt for loop (per ESMs), create monthly climatologies
# months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
vars <- c("sst","dsst","logchl","o2","logsio2","logno3","nstar","sistar","par")
ESMs <- c("IPSL-PISCES","CNRM-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")
rcp <- "rcp85"

### First, get all monthly in situ baseline clims
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d")
res <- lapply(months, function(m) {
            message(paste(m, sep = ""))
            cl <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
            # dim(cl); head(cl)
            cl$id <- factor(paste(cl$x, cl$y, sep = "_"))
            cl <- cl[order(cl$id),]
            # Add month 
            cl$month <- capitalize_str(m)
            return(cl)
        } # eo fun
) # eo lapply
# Rbind
clims_obs <- dplyr::bind_rows(res)
# dim(clims_obs); summary(clims_obs)
colnames(clims_obs)[7] <- "dSST"
# colnames(clims_obs)
# unique(clims_obs$month)

# Need to rotate x coords and adjust cell ids
clims_obs$x2 <- clims_obs$x 
clims_obs[clims_obs$x < 0 ,"x2"] <- (clims_obs[clims_obs$x < 0 ,"x"]) + 360
clims_obs$id <- factor(paste(clims_obs$x2, cl$y, sep = "_"))
clims_obs <- clims_obs[order(clims_obs$id),]

# For testing
m <- "Apr" ; esm <- "GFDL-TOPAZ"

# For every ESM, load the monthly diff climatologies and add them to the observed baseline climatologies, and plot the distrb of each variablexmonthsxESM
for(esm in ESMs) {
 
        message(paste("Preparing future monthly climatologies for ", esm, sep = ""))
        setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
        # Load all monthly climatologies of diff 
        diff_sst <- get(load( paste("clims_mon_diff_","sst","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_dsst <- get(load( paste("clims_ann_diff_","dsst","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_logchl <- get(load( paste("clims_mon_diff_","logchl","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_o2 <- get(load( paste("clims_mon_diff_","o2","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_logsio2 <- get(load( paste("clims_mon_diff_","logsio2","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_logno3 <- get(load( paste("clims_mon_diff_","logno3","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_nstar <- get(load( paste("clims_mon_diff_","nstar","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_sistar <- get(load( paste("clims_mon_diff_","sistar","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        diff_par <- get(load( paste("clims_mon_diff_","par","_",rcp,"_",esm,"_2031-2100.Rdata", sep = "") ))
        # summary(diff_par)
        
        # Create new monthly climatology based on obs insitu + diff
        months2 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
        # m <- "Oct"
        for(m in months2) {
            
            # Subset 'clims_obs'
            message(paste("Preparing future monthly climatologies for ", m, sep = ""))
            mon_subset <- clims_obs[clims_obs$month == m,]
            # Create new monthly clim
            ddf_mon <- data.frame(id = mon_subset$id, x = mon_subset$x2, y = mon_subset$y,
                        SST = (mon_subset$SST) + (diff_sst[,m]), dSST = (mon_subset$dSST) + (diff_dsst[,"dSST"]), 
                        PAR = (mon_subset$PAR) + (diff_par[,m]), dO2 = (mon_subset$dO2) + (diff_o2[,m]),
                        logChl = (mon_subset$logChl) + (diff_logchl[,m]), logNO3 = (mon_subset$logNO3) + (diff_logno3[,m]), 
                        logSiO2 = (mon_subset$logSiO2) + (diff_logsio2[,m]), Nstar = (mon_subset$Nstar) + (diff_nstar[,m]), 
                        Sistar = (mon_subset$Sistar) + (diff_sistar[,m])   
            ) # eo ddf
            # summary(ddf_mon)
            
            ### Need to correct the negative values that might have been created for: dSST, PAR, dO2, logNO3, logSiO2
            if( nrow(ddf_mon[ddf_mon$dSST < 0 & !is.na(ddf_mon$dSST),]) > 0 ) {
                
                n <- nrow(ddf_mon[ddf_mon$dSST < 0 & !is.na(ddf_mon$dSST),])
                value <- min(ddf_mon[,c("dSST")][ddf_mon[,c("dSST")] > 0], na.rm = T) # value
                message(paste("Need to replace ", n, " dSST values < 0 by ", value, sep = ""))
                ddf_mon[ddf_mon$dSST < 0 & !is.na(ddf_mon$dSST),"dSST"] <- value
                
            } # eo if loop for dSST
            
            if( nrow(ddf_mon[ddf_mon$PAR < 0 & !is.na(ddf_mon$PAR),]) > 0 ) {
                
                n <- nrow(ddf_mon[ddf_mon$PAR < 0 & !is.na(ddf_mon$PAR),])
                value <- min(ddf_mon[,c("PAR")][ddf_mon[,c("PAR")] > 0], na.rm = T) # value
                message(paste("Need to replace ", n, " PAR values < 0 by ", value, sep = ""))
                ddf_mon[ddf_mon$PAR < 0 & !is.na(ddf_mon$PAR),"PAR"] <- value
                
            } # eo if loop for PAR
            
            if( nrow(ddf_mon[ddf_mon$dO2 < 0 & !is.na(ddf_mon$dO2),]) > 0 ) {
                
                n <- nrow(ddf_mon[ddf_mon$dO2 < 0 & !is.na(ddf_mon$dO2),])
                value <- min(ddf_mon[,c("dO2")][ddf_mon[,c("dO2")] > 0], na.rm = T) # value
                message(paste("Need to replace ", n, " dO2 values < 0 by ", value, sep = ""))
                ddf_mon[ddf_mon$dO2 < 0 & !is.na(ddf_mon$dO2),"dO2"] <- value
                
            } # eo if loop for dO2
            
            if( nrow(ddf_mon[ddf_mon$logNO3 < 0 & !is.na(ddf_mon$logNO3),]) > 0 ) {
                
                n <- nrow(ddf_mon[ddf_mon$logNO3 < 0 & !is.na(ddf_mon$logNO3),])
                value <- min(ddf_mon[,c("logNO3")][ddf_mon[,c("logNO3")] > 0], na.rm = T) # value
                message(paste("Need to replace ", n, " logNO3 values < 0 by ", value, sep = ""))
                ddf_mon[ddf_mon$logNO3 < 0 & !is.na(ddf_mon$logNO3),"logNO3"] <- value
                
            } # eo if loop for logNO3
            
            if( nrow(ddf_mon[ddf_mon$logSiO2 < 0 & !is.na(ddf_mon$logSiO2),]) > 0 ) {
                
                n <- nrow(ddf_mon[ddf_mon$logSiO2 < 0 & !is.na(ddf_mon$logSiO2),])
                value <- min(ddf_mon[,c("logSiO2")][ddf_mon[,c("logSiO2")] > 0], na.rm = T) # value
                message(paste("Need to replace ", n, " logSiO2 values < 0 by ", value, sep = ""))
                ddf_mon[ddf_mon$logSiO2 < 0 & !is.na(ddf_mon$logSiO2),"logSiO2"] <- value
                
            } # eo if loop for logSiO2
            
            ### And finally, plot change in distrib between monthly obs and monthly future based on obs + delta
    		m.f.clim <- melt(ddf_mon, id.vars = c("id","x","y"))
    		vars <- c("SST","dSST","logChl","dO2","logSiO2","logNO3","Nstar","Sistar","PAR")
    		mon_subset2 <- mon_subset[,c("id","x2","y",vars)]
    		m.obs <- melt(mon_subset2, id.vars = c("id","x2","y"))
    		# setdiff(unique(m.obs$id), unique(m.f.clim$id))
    		# setdiff(unique(m.obs$variable), unique(m.f.clim$variable))
    		# setdiff(unique(m.obs$x2), unique(m.f.clim$x))
    		m.obs$group <- "baseline"
    		m.f.clim$group <- "future"
    		colnames(m.obs)[2] <- "x"
    		# Rbind
    		data2rbind <- rbind(m.obs, m.f.clim)
    		# dim(data2rbind) ; head(data2rbind) ; str()
    		# Plot distrbution of "value", facet per variable and colour per "group"
    		plot <- ggplot(data2rbind, aes(x = value, fill = factor(group))) + geom_histogram(alpha=.5, position="identity") + 
    				scale_fill_manual(name = m, labels = c("Baseline","Future"), values = c("#2166ac","#b2182b")) +
    				theme_light() + facet_wrap(factor(data2rbind$variable), ncol = 3, nrow = 3, scales = "free")
		
    		ggsave(plot = plot, filename = paste("plot_shift_distrib_",m,"_",esm,".pdf", sep = ""), height = 10, width = 15, dpi = 300)	
	
            rm(data2rbind, plot, m.f.clim, m.obs ,mon_subset2) ; gc()
            
            ### Return if ddf_mon presents right dimensions (64800x12)
            if( length(ddf_mon) == 12 & nrow(ddf_mon) == 64800 ) {
            
                message(paste("Saving ddf_mon for ", m, " for ESM == ", esm,  sep = ""))
                setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/")
                write.table(x = ddf_mon, file = paste("clims_mon_",m,"_",esm,"_rcp85_base+2100-2031.txt", sep = ""), sep = "\t")
            
            } else {
            
                message(paste("ERROR !!!!!! ddf_mon has wrong dimensions !!!!!!", sep = ""))
                
            }
            
            # Clean stuff and continue for loop
            gc()
            
            setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
            message(paste("", sep = ""))
            
        } # eo for loop - m in months
    
} # eo for loop - esm in ESMs
    
    
# ------------------------------------------------------------ 

### 19/11/19: Check the final future conditions by mapping 

setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims")
wd <- getwd()
esm <- "GFDL-TOPAZ"
for(esm in ESMs) {
        
        # Useless message
        message(paste("", sep = ""))
        message(paste("Mapping future fields for ", esm, sep = ""))
      
        setwd( paste(wd,"/",esm,"/", sep = "") )
        # Read
        months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
        res <- lapply(months, function(m) {
                        cl <- read.table(paste("clims_mon_",m,"_",esm,"_rcp85_base+2100-2031.txt", sep = ""), sep = "\t", h = T)
                        cl$month <- m 
                        return(cl)
            } # eo FUN
        ) # eo lapply
        clims_fut <- dplyr::bind_rows(res)
        # dim(clims_fut) ; summary(clims_fut)
        rm(res); gc()
        
        # Need to melt and then map per variable and facet per month !
        molten <- melt(clims_fut, id.vars = c("id","x","y","month"))
        # colnames(molten) ; dim(molten)
        # vv <- "dSST"
        for(vv in vars) {
            
            message(paste("Mapping future fields of ", vv, sep = ""))
            sub <- molten[molten$variable == vv,]
            # dim(sub) 12*64800 rows
            # Find a way to automatically define the levels of the isopleths for the geom_contour
            # quantile(sub$value, na.rm = T)
            if( vv == "dSST" ) {
                
                map1 <- ggplot(data = sub[sub$month == "Apr",]) + geom_raster(aes(x = x, y = y, fill = value) ) +
                    geom_contour(colour = "grey60", binwidth = 5, size = 0.25, aes(x = x, y = y, z = value) ) +
                    scale_fill_viridis(name = vv) + geom_polygon(aes(x = long, y = lat, group = group), 
                                data = world2, fill = "grey70", colour = "black", size = 0.3) +
                 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
                 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + theme_bw()
                        
                ggsave(plot = map1, filename = paste("maps_mon_",vv,"_rcp85_",esm,"_base+2100-2031.pdf", sep = ""), height = 6, width = 3, dpi = 300)        
                    
            } else {
                
                bin <- round( abs(quantile(sub$value, na.rm = T)[5]) - abs(quantile(sub$value, na.rm = T)[4]) ) # bin
                map2 <- ggplot(data = sub) + geom_raster(aes(x = x, y = y, fill = value) ) +
                    geom_contour(colour = "grey60", binwidth = bin, size = 0.25, aes(x = x, y = y, z = value) ) +
                    scale_fill_viridis(name = vv) + coord_quickmap() + theme_bw() + 
                    facet_wrap(~ factor(month), ncol = 3) 
                    
                ggsave(plot = map2, filename = paste("maps_mon_",vv,"_rcp85_",esm,"_base+2100-2031.pdf", sep = ""), height = 16, width = 20, dpi = 300)    
                    
            } # eo if else loop 

        } # eo for loop - v in vars
        
        message(paste("", sep = ""))
    
} # eo for loop - for esm in ESMs


# ------------------------------------------------------------ 

### 15/01/21: Comparing the ranges of of future changes in environmental predictors: sst, dsst VS. logsio2, logno3, nstar, sistar 

setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
wd <- getwd()
dir()

### For each ESM and each v in vars: 
# - load baseline climatology
# - load diff climatology
# - compute and return % change
require("parallel")

vars <- c("sst","logsio2","logno3","nstar","sistar","o2","dsst")
ESMs <- c("IPSL-PISCES","CNRM-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")

# For testing: 
v <- "sst"
esm <- "IPSL-PISCES"

res.vars <- mclapply(vars, function(v) {
    
            message(paste(v, sep = ""))
            message(paste("", sep = ""))
            message(paste("", sep = ""))
            
            res.ESM <- lapply(ESMs, function(esm) {
                
                        message(paste(esm, sep = ""))
                        
                        if(v %in% c("sst","logsio2","logno3","nstar","sistar","o2")) {
                            
                            # Load baseline
                            base <- get(load(paste("clims_mon_",v,"_rcp85_",esm,"_2100-2081.Rdata", sep = "")))
                            # Load diff
                            diff <- get(load(paste("clims_mon_diff_",v,"_rcp85_",esm,"_2031-2100.Rdata", sep = "")))
                            # dim(base); dim(diff)
                            # head(base$id); head(diff$id)
                        
                            if(v == "o2" & length(base == 16)) {
                                # Drop 'Annual col' in 
                                base <- base[,c(1:3,5:length(base))]
                            }
                        
                            # Melt both and cbind
                            m.base <- melt(base, id.vars = c("id","x","y"))
                            colnames(m.base)[c(4:5)] <- c('month','base')
                            m.diff <- melt(diff, id.vars = c("id","x","y"))
                            colnames(m.diff)[c(4:5)] <- c('month','diff')
                            # dim(m.diff); dim(m.diff)
                        
                            # Cbind
                            m.base$diff <- m.diff$diff
                            rm(diff,m.diff) ; gc()
                            # Compute %
                            m.base$perc <- (m.base$diff/m.base$base)*100
                            # quantile(m.base[abs(m.base$perc) < 1000,"perc"], na.rm = T, probs = seq(0,1,0.05))
                            # Narrow down between 5th & 95th percentiles
                            low.bound <- quantile(m.base[abs(m.base$perc) < 1000,"perc"], na.rm = T, probs = seq(0,1,0.05))[2]
                            up.bound <- quantile(m.base[abs(m.base$perc) < 1000,"perc"], na.rm = T, probs = seq(0,1,0.05))[20]
                        
                            m.base <- m.base[m.base$perc > low.bound,]
                            m.base <- m.base[m.base$perc < up.bound,]
                        
                            # For each cell, compute mean annual %
                            clim <- data.frame(m.base %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), mean = mean(perc, na.rm = T)) )
                            summary(clim)
                            clim$ESM <- esm
                        
                            return(clim)
                            
                        } else if(v == "dsst") {
                            
                            # Load baseline
                            base <- get(load(paste("clims_ann_",v,"_rcp85_",esm,"_2100-2081.Rdata", sep = "")))
                            # Load diff
                            diff <- get(load(paste("clims_ann_diff_",v,"_rcp85_",esm,"_2031-2100.Rdata", sep = "")))
                            # dim(base); dim(diff)
                            # head(base$id); head(diff$id)
                        
                            # Melt both and cbind
                            m.base <- melt(base, id.vars = c("id","x","y"))
                            colnames(m.base)[c(4:5)] <- c('month','base')
                            m.diff <- melt(diff, id.vars = c("id","x","y"))
                            colnames(m.diff)[c(4:5)] <- c('month','diff')
                            # dim(m.diff); dim(m.diff)
                        
                            # Cbind
                            m.base$diff <- m.diff$diff
                            rm(diff,m.diff) ; gc()
                            # Compute %
                            m.base$perc <- (m.base$diff/m.base$base)*100
                            # quantile(m.base[abs(m.base$perc) < 1000,"perc"], na.rm = T, probs = seq(0,1,0.05))
                            # Narrow down between 5th & 95th percentiles
                            low.bound <- quantile(m.base[abs(m.base$perc) < 1000,"perc"], na.rm = T, probs = seq(0,1,0.05))[2]
                            up.bound <- quantile(m.base[abs(m.base$perc) < 1000,"perc"], na.rm = T, probs = seq(0,1,0.05))[20]
                        
                            m.base <- m.base[m.base$perc > low.bound,]
                            m.base <- m.base[m.base$perc < up.bound,]
                        
                            # For each cell, compute mean annual %
                            clim <- data.frame(m.base %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), mean = mean(perc, na.rm = T)) )
                            # summary(clim)
                            clim$ESM <- esm
                        
                            return(clim)
                            
                        } # eo if else loop 
                
                } # eo lapply
                
            ) # eo lapply
            
            # Rbind
            table <- dplyr::bind_rows(res.ESM)
            # dim(table) ; summary(table)
            rm(res.ESM)
            table$var <- v
            
            return(table)
    
    }, mc.cores = length(vars) 
    
) # eo mclapply
# Rbind and plot distributions across vars and ESMs (facet)
ddf <- dplyr::bind_rows(res.vars)
head(ddf) ; dim(ddf)
summary(ddf) 
rm(res.vars) ; gc()

# Adjust some labels
ddf$var <- factor(ddf$var)
levels(factor(ddf$var))

levels(ddf$var)[levels(ddf$var) == "sst"] <- "SST"
levels(ddf$var)[levels(ddf$var) == "sistar"] <- "Si*"
levels(ddf$var)[levels(ddf$var) == "nstar"] <- "N*"
levels(ddf$var)[levels(ddf$var) == "logsio2"] <- "Silicates"
levels(ddf$var)[levels(ddf$var) == "logno3"] <- "Nitrates"
levels(ddf$var)[levels(ddf$var) == "o2"] <- "Oxygen (175m)"
levels(ddf$var)[levels(ddf$var) == "dsst"] <- "SST range"

data.frame(ddf %>% group_by(var) %>% summarize(mean = mean(mean, na.rm = T)))
data.frame(ddf %>% group_by(var) %>% summarize(med = median(mean, na.rm = T)))


### Plot
ggplot(data = ddf[abs(ddf$mean) < 100,], aes(x = var, y = mean)) + geom_boxplot(colour = "black", fill = "grey65") +
    theme_bw() + facet_wrap(.~factor(ESM), scales = "free_y") + geom_hline(yintercept = 0, linetype = "dashed")

### Plot
plot1 <- ggplot(data = na.omit(ddf[abs(ddf$mean) < 100,]), aes(x = var, y = mean)) + 
    geom_boxplot(colour = "black", fill = "grey65") + xlab("") +
    ylab("Mean annual difference (%)\n(2081-2100 minus 2012-2031)") +
    theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave(plot = plot, filename = "plot_distrib_mean_ann_perc_predictors_2100-2031_all_ESMs.jpg", dpi = 300, width = , height = )

### Split per lat bands: tropics/temper/polar
ddf$domain <- NA
ddf[which(abs(ddf$y) < 30),'domain'] <- "Tropical"
ddf[which(abs(ddf$y) > 30),'domain'] <- "Temperate"
ddf[which(abs(ddf$y) > 60),'domain'] <- "Polar"

plot2 <- ggplot(data = na.omit(ddf[abs(ddf$mean) < 100,]), aes(x = var, y = mean)) + 
    geom_boxplot(colour = "black", fill = "grey65") + xlab("") +
    ylab("Mean annual difference (%)\n(2081-2100 minus 2012-2031)") +
    theme_bw() + facet_wrap(.~factor(domain), scales = "free_y") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave(plot = plot, filename = "plot_distrib_mean_ann_perc_predictors_2100-2031_all_ESMs_domains.jpg", dpi = 300, width = 9, height = 5)

library("ggpubr")
panel <- ggarrange(plot1, plot2, labels = c("A","B"), align = "hv", ncol = 2, nrow = 1, widths = c(1,3))
setwd(wd)
ggsave(plot = panel, filename = "plot_distrib_mean_ann_perc_predictors_2100-2031_all_ESMs.jpg", dpi = 300, width = 15, height = 4.5)


### Test
library("PMCMR")
kruskal.test(mean ~ factor(var), data = na.omit(ddf[abs(ddf$mean) < 100,]))
# Kruskal-Wallis chi-squared = 354543, df = 6, p-value < 2.2e-16

kruskal.test(mean ~ factor(var), data = na.omit(ddf[abs(ddf$mean) < 100 & ddf$domain == "Tropical",]))
# Kruskal-Wallis chi-squared = 92212, df = 6, p-value < 2.2e-16
kruskal.test(mean ~ factor(var), data = na.omit(ddf[abs(ddf$mean) < 100 & ddf$domain == "Temperate",]))
# Kruskal-Wallis chi-squared = 186067, df = 6, p-value < 2.2e-16
kruskal.test(mean ~ factor(var), data = na.omit(ddf[abs(ddf$mean) < 100 & ddf$domain == "Polar",]))
# Kruskal-Wallis chi-squared = 106925, df = 6, p-value < 2.2e-16

posthoc.kruskal.dunn.test(mean ~ factor(var), p.adjust = "bonf", data = na.omit(ddf[abs(ddf$mean) < 100 & ddf$domain == "Temperate",]))





