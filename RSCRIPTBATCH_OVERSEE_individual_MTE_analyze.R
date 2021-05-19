
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
library("parallel")

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

### 1°) Get ensembles of baseline annual species richness
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# format: "table_ann_changes_",esm,"_",sdm,"_",p,".Rdata"
files <- dir()[grep("div_baseline_",dir())]; files
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
    }, mc.cores = 30
) # eo mclapply
# Rbind
table <- dplyr::bind_rows(res)
dim(table)
head(table)
rm(res); gc()
table <- table[order(table$id),]

### Computing ensembles (all, SDM/ESM/pool)
ens <- data.frame(table %>% group_by(id,group) %>%
        summarize(x = unique(x), y = unique(y), mean.rich = mean(rich, na.rm = T) ) 
) # ddf
summary(ens)
head(ens)
# Dcast to separate between 
dens_base <- dcast(ens, id + x + y ~ group, value.var = c("mean.rich") )
dim(dens_base); head(dens_base)


### 2°) Get ensembles of changes in annual species richness
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
dim(table); head(table)
rm(res); gc()

### Computing ensembles (all, SDM/ESM/pool)
ens_fut <- data.frame(table %>% group_by(id) %>% 
        summarize(x = unique(x), y = unique(y), diff_phyto = mean(diff_rich_phyto, na.rm = T), diff_zoo = mean(diff_rich_zoo, na.rm = T) ) 
) # ddf
summary(ens_fut)
# colnames(ens)

### Combine "ens_fut" with "dens_base" to get fut richness (not diff)
ens_fut <- ens_fut[order(ens_fut$id),]
dens_base <- dens_base[order(dens_base$id),]
commons <- intersect(unique(ens_fut$id), unique(dens_base$id)) # length(commons)
ens_fut <- ens_fut[ens_fut$id %in% commons,]
dens_base <- dens_base[dens_base$id %in% commons,]
dim(ens_fut) ; dim(dens_base)
head(ens_fut); head(dens_base)
dens_base$diff_phyto <- ens_fut$diff_phyto
dens_base$diff_zoo <- ens_fut$diff_zoo
colnames(dens_base)[c(4,5)] <- c("base_phyto","base_zoo")

# Derive future richness: 
dens_base$fut_phyto <- (dens_base$base_phyto)+(dens_base$diff_phyto)
dens_base$fut_zoo <- (dens_base$base_zoo)+(dens_base$diff_zoo)
summary(dens_base)
# OK all set

# ------------------------------------------------------------

### 3°) Get baseline and future SST monthly cliamtologies 

### A) Baseline
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
		summarise(x = unique(x2), y = unique(y), Bathy = mean(Bathy,na.rm=T), SST = mean(SST,na.rm=T), SSS = mean(SSS,na.rm=T)
	) # eo summarise	
) # eo ddf
# Exclude non opean ocean cells: SSS > 30 and Bathy < -175
aclim <- aclim[aclim$SSS >= 20,]
aclim <- aclim[aclim$Bathy <= -175,]
# Check
aclim <- aclim[!is.na(aclim$id),]
aclim <- aclim[order(aclim$id),] # ok
base_sst <- aclim[,c("id","x","y","SST")]
rm(aclim,clims,files) ; gc()


### B) Future 
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims")
# One directory per ESM...
# esm <- "MRI-NEMURO"
res <- lapply(ESMs, function(esm) {
            
            message(paste("Loading future SST for ",esm, sep = ""))
            setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims/",esm, sep = ""))
            files <- dir() # f <- files[1]
            monthlies <- lapply(files, function(f) {
                        d <- read.table(f, h = T, sep = "\t")
                        return(d[,c("id","x","y","SST")])
                }
            ) # eo lapply
            # Rbind
            table <- bind_rows(monthlies)          
            # dim(table); head(table)      
            # Compute annual mean SST
            mean <- data.frame(table %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), SST = mean(SST,na.rm=T)))
            #dim(mean)
            mean <- mean[order(mean$id),]
            mean$ESM <- esm
            return(mean)
    
    } # eo FUN
    
) # eo lapply
# Rbind
table <- bind_rows(res)
head(table); dim(table)
# Compute ensemble mean
fut_sst <- data.frame(table %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), SST = mean(SST, na.rm = T)))

### Check if base_sst and fut_sst are ok
# summary(base_sst)
# summary(fut_sst)
# dim(base_sst) ; dim(fut_sst)
commons <- intersect(unique(base_sst$id), unique(fut_sst$id)) # length(commons)
base_sst <- base_sst[base_sst$id %in% commons,]
fut_sst <- fut_sst[fut_sst$id %in% commons,]
colnames(base_sst)[4] <- "SST_base" 
colnames(fut_sst)[4] <- "SST_fut"
# Provide fut SST to "base_sst"
base_sst$SST_fut <- fut_sst$SST_fut
summary(base_sst)


# ------------------------------------------------------------

### 4°) Merge "dens_base" and "base_sst" to examine contemporary and future MTE predictions
commons <- intersect(unique(base_sst$id), unique(dens_base$id)) # length(commons)
base_sst <- base_sst[base_sst$id %in% commons,]
dens_base <- dens_base[dens_base$id %in% commons,]
dim(base_sst) ; dim(dens_base)
dens_base$SST_base <- base_sst$SST_base
dens_base$SST_fut <- base_sst$SST_fut

### Ok, convert to eV
# base
dens_base$absT <- dens_base$SST_base + 273.15
dens_base$kT_base <- 1 / (dens_base$absT*-8.6173303*10^(-5) )
# fut
dens_base$absT <- dens_base$SST_fut + 273.15
dens_base$kT_fut <- 1 / (dens_base$absT*-8.6173303*10^(-5) )
summary(dens_base)
#Chekc log scale for richness (y scales)
summary(log(dens_base$base_phyto)) # 3-5.5
summary(log(dens_base$base_zoo)) # 3.5-5.5
summary(log(dens_base$fut_phyto)) # 2.9-5.5
summary(log(dens_base$fut_zoo)) # 3.5-5.3
# limits = c(2.9,5.5)

# For x axis
summary(dens_base$kT_base) # -42.8 - -38.3
summary(dens_base$kT_fut) # -42.71 - -37.9
# limits = c(-42.8,-37.9)

### Save objetc so you don't have to re-do the codes above anymore
setwd(WD)
# save(dens_base, file = "table_4MTE_06_04_20.Rdata")
dens_base <- get(load("table_4MTE_06_04_20.Rdata"))

### Check relationshps between ln(SR) ~ eV
# ggplot() + geom_point(aes(x = kT_base, y = log(base_phyto)), data = dens_base, colour = "grey70") +
#             geom_smooth(aes(x = kT_base, y = log(base_phyto)), data = dens_base, colour = "black", method = "lm") +
#             scale_y_continuous(limits = c(2.9,5.5)) + ylab("Phytoplankton species richness (ln)") +
#             xlab("Thermal energy (1/kT)") + theme_classic()
# #
# ggplot() + geom_point(aes(x = kT_base, y = log(base_zoo)), data = dens_base, colour = "grey70") +
#             geom_smooth(aes(x = kT_base, y = log(base_zoo)), data = dens_base, colour = "black", method = "lm") +
#             scale_y_continuous(limits = c(2.9,5.5)) + ylab("Zooplankton species richness (ln)") +
#             xlab("Thermal energy (1/kT)") + theme_classic()
head(dens_base)            

### And when above the threshold t 
for(t in seq(from = -1, to = 30, by = 0.5)) {
     lm <- lm(log(base_phyto) ~ kT_base, data = dens_base[dens_base$SST_base >= t,])
     message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
     message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
     message(paste("", sep = ""))
} 

# Same with zoo
for(t in seq(from = -1, to = 30, by = 0.5)) {
      lm <- lm(log(base_zoo) ~ kT_base, data = dens_base[dens_base$SST_base >= t,])
      message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 

### And based on future SST ?
for(t in seq(from = -1, to = 33, by = 0.1)) {
      lm <- lm(log(fut_phyto) ~ kT_fut, data = ddf[ddf$SST_fut >= t,])
      message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 

for(t in seq(from = -1, to = 33, by = 0.1)) {
      lm <- lm(log(fut_zoo) ~ kT_fut, data = ddf[ddf$SST_fut >= t,])
      message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 


### And from 11°C to an upper threshold t 
for(t in seq(from = 12, to = 30, by = 0.5)) {
     lm <- lm(log(base_phyto) ~ kT_base, data = dens_base[dens_base$SST_base < t & dens_base$SST_base > 11,] )
     message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
     message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
     message(paste("", sep = ""))
} 

# Same with zoo
for(t in seq(from = 12, to = 30, by = 0.5)) {
      lm <- lm(log(base_zoo) ~ kT_base, data = dens_base[dens_base$SST_base < t & dens_base$SST_base > 11,] )
      message(paste("Above t = ",t,"°C || Slope to kT = ", round(lm$coefficients[2],4), sep = ""))
      message(paste("Above t = ",t,"°C || r2 to kT = ", round(summary(lm)$adj.r.squared,4), sep = ""))
      message(paste("", sep = ""))
} 

### Linear fits that best match MTE predictions
# - for phyto: from 22°C (-39.32 eV) to max SST
# - for zoo: from 11°C (-40.84 eV) to 20°C (-39.58 eV)

### To find matches between eV and °C
summary(dens_base[which(dens_base$SST_fut > 29.4 & dens_base$SST_fut < 29.6),"kT_fut"])

### To find SST/kT ranges for tropics in the future ?
summary(ddf[which(ddf$y > -25 & ddf$y < 25),c("kT_base","SST_base","kT_fut","SST_fut")])


# ------------------------------------------------------------

library("tidyverse")
library("RColorBrewer")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("parallel")

### To plot MTE predictions over data points and fitted polynomial
WD <- getwd()
setwd(WD)
dens_base <- get(load("table_4MTE_06_04_20.Rdata"))

ddf <- na.omit(dens_base)
slopeP <- 0.32
slopeZ <- 0.65
interP <- mean(log(ddf$base_phyto) - slopeP * ddf$kT_base)
interZ <- mean(log(ddf$base_zoo) - slopeZ * ddf$kT_base)
# Derive expected model
ddf$phyto.mte.base <- slopeP*(ddf$kT_base) + interP
ddf$zoo.mte.base <- slopeZ*(ddf$kT_base) + interZ

# Same based on fut
interP <- mean(log(ddf$fut_phyto) - slopeP * ddf$kT_fut)
interZ <- mean(log(ddf$fut_zoo) - slopeZ * ddf$kT_fut)
ddf$phyto.mte.fut <- slopeP*(ddf$kT_fut) + interP
ddf$zoo.mte.fut <- slopeZ*(ddf$kT_fut) + interZ
# summary(ddf)

### Add 3rd order polynomials
lm3.phyto <- lm(log(base_phyto) ~ kT_base + I(kT_base^2) + I(kT_base^3), data = ddf, na.action = na.exclude)
summary(lm3.phyto) # R-squared:  0.8396 
lm3.zoo <- lm(log(base_zoo) ~ kT_base + I(kT_base^2) + I(kT_base^3), data = ddf, na.action = na.exclude)
summary(lm3.zoo) # R-squared:  0.9704 
ddf$poly_phyto_base <- predict(lm3.phyto)
ddf$poly_zoo_base <- predict(lm3.zoo)

### And ft polynom for future rich as well
lm3.phyto <- lm(log(fut_phyto) ~ kT_fut + I(kT_fut^2) + I(kT_fut^3), data = ddf, na.action = na.exclude)
summary(lm3.phyto) # R-squared: 0.8782
lm3.zoo <- lm(log(fut_zoo) ~ kT_fut + I(kT_fut^2) + I(kT_fut^3), data = ddf, na.action = na.exclude)
summary(lm3.zoo) # R-squared: 0.9638
ddf$poly_phyto_fut <- predict(lm3.phyto)
ddf$poly_zoo_fut <- predict(lm3.zoo)

summary(ddf)

### Plot MTE predictions versus polynom, try to add secondary x axis for SST
# Find the proper labels to put as secondary SST axis
summary(dens_base[which(dens_base$kT_fut > -38.05),"SST_base"])

### Phytopl
plot1 <- ggplot() + #geom_vline(xintercept = -38.9, linetype = "dotted") + geom_vline(xintercept = -38.1, linetype = "dotted") + 
    geom_point(data = ddf, aes(x = kT_base, y = log(base_phyto) ), colour = "grey75", size = 1, alpha = .2) + 
    stat_density_2d(data = ddf, aes(x = kT_base, y = log(base_phyto) ), 
                    colour = "#4d9221", fill = "#4d9221", geom = "polygon", bins = 6, alpha = .2) +
    geom_line(aes(x = kT_base, y = phyto.mte.base), data = ddf, colour = "black", linetype = "dashed") +    
    geom_line(aes(x = kT_base, y = poly_phyto_base), data = ddf, colour = "black") +
    scale_y_continuous(limits = c(2.9,5.5), sec.axis = dup_axis(name="")) + 
    scale_x_continuous(limits = c(-42.8,-37.9), sec.axis = sec_axis(~.*1, name = "Mean anual SST (°C)",
                    breaks = c(-42,-41,-40,-39,-38), labels = c("3.1","9.9","17.0","24.4","29.4"))) + 
    ylab("Mean annual phytoplankton species richness (log)") + 
    xlab("Mean annual thermal energy (1/eV)") + theme_classic()
#
plot2 <- ggplot() + geom_vline(xintercept = -38.9, linetype = "dotted") + geom_vline(xintercept = -38.1, linetype = "dotted") + 
    geom_point(data = ddf, aes(x = kT_fut, y = log(fut_phyto) ), colour = "grey75", size = 1, alpha = .2) + 
    stat_density_2d(data = ddf, aes(x = kT_fut, y = log(fut_phyto) ), 
                colour = "#4d9221", fill = "#4d9221", geom = "polygon", bins = 6, alpha = .2) +
    geom_line(aes(x = kT_fut, y = phyto.mte.fut), data = ddf, colour = "black", linetype = "dashed") +    
    geom_line(aes(x = kT_fut, y = poly_phyto_fut), data = ddf, colour = "black") +
    scale_y_continuous(limits = c(2.9,5.5), sec.axis = dup_axis(name="")) + 
    scale_x_continuous(limits = c(-42.8,-37.9), sec.axis = sec_axis(~.*1, name = "Mean anual SST (°C)",
                    breaks = c(-42,-41,-40,-39,-38), labels = c("3.1","9.9","17.0","24.4","29.4"))) + 
    ylab("Mean annual phytoplankton species richness (log)") + 
    xlab("Mean annual thermal energy (1/eV)") + theme_classic()

### Zoo            
plot3 <- ggplot() + #geom_vline(xintercept = -38.9, linetype = "dotted") + geom_vline(xintercept = -38.1, linetype = "dotted") + 
    geom_point(data = ddf, aes(x = kT_base, y = log(base_zoo) ), colour = "grey75", size = 1, alpha = .2) + 
    stat_density_2d(data = ddf, aes(x = kT_base, y = log(base_zoo) ), 
                colour = "#c51b7d", fill = "#c51b7d", geom = "polygon", bins = 6, alpha = .2) +
    geom_line(aes(x = kT_base, y = zoo.mte.base), data = ddf, colour = "black", linetype = "dashed") +    
    geom_line(aes(x = kT_base, y = poly_zoo_base), data = ddf, colour = "black") +
    scale_y_continuous(limits = c(2.9,5.5), sec.axis = dup_axis(name="")) + 
    scale_x_continuous(limits = c(-42.8,-37.9), sec.axis = sec_axis(~.*1, name = "Mean anual SST (°C)",
                    breaks = c(-42,-41,-40,-39,-38), labels = c("3.1","9.9","17.0","24.4","29.4"))) + 
    ylab("Mean annual zooplankton species richness (log)") + 
    xlab("Mean annual thermal energy (1/eV)") + theme_classic()
#
plot4 <- ggplot() + geom_vline(xintercept = -38.9, linetype = "dotted") + geom_vline(xintercept = -38.1, linetype = "dotted") + 
    geom_point(data = ddf, aes(x = kT_fut, y = log(fut_zoo) ), colour = "grey75", size = 1, alpha = .2) + 
    stat_density_2d(data = ddf, aes(x = kT_fut, y = log(fut_zoo) ), 
                colour = "#c51b7d", fill = "#c51b7d", geom = "polygon", bins = 6, alpha = .2) +
    geom_line(aes(x = kT_fut, y = zoo.mte.fut), data = ddf, colour = "black", linetype = "dashed") +    
    geom_line(aes(x = kT_fut, y = poly_zoo_fut), data = ddf, colour = "black") +
    scale_y_continuous(limits = c(2.9,5.5), sec.axis = dup_axis(name="")) + 
    scale_x_continuous(limits = c(-42.8,-37.9), sec.axis = sec_axis(~.*1, name = "Mean anual SST (°C)",
                    breaks = c(-42,-41,-40,-39,-38), labels = c("3.1","9.9","17.0","24.4","29.4"))) + 
    ylab("Mean annual zooplankton species richness (log)") + 
    xlab("Mean annual thermal energy (1/eV)") + theme_classic()

### Identify the eV range that correspond to a decrease in zoo SR but increase in phyto SR
#summary(ddf[which(ddf$diff_zoo < 0 & ddf$diff_phyto > 0 & ddf$SST_base > 15),c("SST_base","SST_fut","kT_base","kT_fut")])

### Idea, for future relationships, omit contours and color points as a fun of diff
plot5 <- ggplot() + geom_vline(xintercept = -38.9, linetype = "dotted") + geom_vline(xintercept = -38.1, linetype = "dotted") + 
        geom_point(data = ddf, aes(x = kT_fut, y = log(fut_phyto), colour = diff_phyto), size = 1, alpha = .8) +
        scale_colour_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", guide = F) +
        geom_line(aes(x = kT_fut, y = phyto.mte.fut), data = ddf, colour = "black", linetype = "dashed") +    
        geom_line(aes(x = kT_fut, y = poly_phyto_fut), data = ddf, colour = "black") +
        scale_y_continuous(limits = c(2.9,5.5), sec.axis = dup_axis(name="")) + 
        scale_x_continuous(limits = c(-42.8,-37.9), sec.axis = sec_axis(~.*1, name = "SST (°C) - Future",
                    breaks = c(-42,-41,-40,-39,-38), labels = c("3.1","9.9","17.0","24.4","29.4"))) + 
        ylab("Mean annual phytoplankton species richness (log)") + 
        xlab("Thermal energy (1/eV) - Future") + theme_classic() + theme(legend.position = "bottom")
            
plot6 <- ggplot() + geom_vline(xintercept = -38.9, linetype = "dotted") + geom_vline(xintercept = -38.1, linetype = "dotted") + 
        geom_point(data = ddf, aes(x = kT_fut, y = log(fut_zoo), colour = diff_zoo), size = 1, alpha = .8) +
        scale_colour_gradient2(name = "Richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white", guide = F) +
        geom_line(aes(x = kT_fut, y = zoo.mte.fut), data = ddf, colour = "black", linetype = "dashed") +    
        geom_line(aes(x = kT_fut, y = poly_zoo_fut), data = ddf, colour = "black") +
        scale_y_continuous(limits = c(2.9,5.5), sec.axis = dup_axis(name="")) + 
        scale_x_continuous(limits = c(-42.8,-37.9), sec.axis = sec_axis(~.*1, name = "SST (°C) - Future",
                    breaks = c(-42,-41,-40,-39,-38), labels = c("3.1","9.9","17.0","24.4","29.4"))) + 
        ylab("Mean annual zooplankton species richness (log)") + 
        xlab("Thermal energy (1/eV) - Future") + theme_classic() + theme(legend.position = "bottom")

### Save plots
ggsave(plot = plot1, filename = "plot_rich_phyto_MTE_base.jpg", dpi = 300, width = 5, height = 5)
ggsave(plot = plot2, filename = "plot_rich_phyto_MTE_fut.jpg", dpi = 300, width = 5, height = 5)
ggsave(plot = plot3, filename = "plot_rich_zoo_MTE_base.jpg", dpi = 300, width = 5, height = 5)
ggsave(plot = plot4, filename = "plot_rich_zoo_MTE_fut.jpg", dpi = 300, width = 5, height = 5)
# And diff
ggsave(plot = plot5, filename = "plot_diff_phyto_MTE_v2.jpg", dpi = 300, width = 5, height = 5)
ggsave(plot = plot6, filename = "plot_diff_zoo_MTE_v2.jpg", dpi = 300, width = 5, height = 5)

### Arrange in panel for MS
require("ggpubr")
# 1st panel: MTE plots for baseline
ggarrange(plot1, plot3, labels = c("A","B"), ncol = 2, nrow = 1, align = "h", common.legend = T)
panel <- ggarrange(plot1, plot3, labels = c("A","B"), ncol = 2, nrow = 1, align = "h", common.legend = T)
ggsave(plot = panel, filename = "Fig.2A-B.jpg", dpi = 300, width = 8, height = 4)

# 2nd panel: plots for baseline + future
panel2 <- ggarrange(plot1, plot3, plot5, plot6, labels = c("A","B","C","D"), ncol = 2, nrow = 2, common.legend = F)
ggsave(plot = panel2, filename = "Fig.2v2.A-D.jpg", dpi = 300, width = 8, height = 8)


### 07/04/20: Detecting turnpoints
infl <- c(F, diff(diff(ddf[order(ddf$kT_fut, decreasing = F),"poly_zoo_fut"])>0) !=0 )
infl[infl == T]
# 2 iflx points as showed on the curve on the plot; 5373 22472 and 
match(TRUE,infl)

fit <- ddf[order(ddf$kT_fut, decreasing = F),c("kT_fut","SST_fut")]
fit[c(5373,22472),]

### Same for phyto?
infl <- c(F, diff(diff(ddf[order(ddf$kT_fut, decreasing = F),"poly_phyto_fut"])>0) !=0 )
infl[infl == T] # 1



