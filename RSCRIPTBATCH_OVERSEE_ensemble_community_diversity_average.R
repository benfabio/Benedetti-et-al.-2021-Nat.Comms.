
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

world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Set the working directories, vectors etc.
WD <- getwd()
setwd("/net/hydro/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85")

months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of pools
pools <- c("p1","p2","p3","p4")
# Vector of earth system models
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
rcp <- "rcp85"

# ### Make sure all tables have similar dimensions
# for(f in dir()) {
#
#     message(paste("", sep = ""))
#     t <- read.table(f, sep = "\t")
#     message(paste("Dimensions for ",f," are: ", sep = ""))
#     message(dim(t))
#     message(paste("", sep = ""))
#
# } # eo for loop

### Examine some outputs for phyto & zoo
# test.table1 <- read.table("table_zoo_mon_composition_2100-2000_MRI-NEMURO_rcp85_nov.txt", sep = "\t")
# dim(test.table1); head(test.table1[,c(1:10)])
#
# test.table2 <- read.table("table_zoo_mon_composition_2100-2000_CESM-BEC_rcp85_may.txt", sep = "\t")
# dim(test.table2); head(test.table2[,c(1:10)])

### With mclapply, extract the monthly composition tables and average across all pools of env predictors, save mean tables and derive diversity estimates

# f <- "table_zoo_mon_composition_baseline_jul.txt"
community.averager <- function(f = files) {
    
            # Useless message
            group <- do.call(cbind,strsplit(as.character(f), split = "_"))[2,1]
            #esm <- do.call(cbind,strsplit(as.character(f), split = "_"))[6,1]
            m <- do.call(cbind,strsplit(as.character(f), split = "_"))[6,1]
            m <- str_replace(m,".txt","") # m
            message(paste("Averaging species' HSI across all 4 pools for the baseline period, for ",m,  sep = ""))
            
            # Load data table 
            table <- read.table(f, sep = "\t")
            # dim(table); head(table[,c(1:10)]) # summary(table[,c(1:10)])
            # melt and use dplyr to compute mean HSI per species before dcasting
            mtable <- reshape2::melt(data = table, id.vars = colnames(table)[c(1:4)])
            rm(table); gc()
            # head(mtable)
            require("dplyr")
            ddf <- data.frame(mtable %>% group_by(cell_id, variable) %>% summarize(x = unique(x), y = unique(y), mean = mean(value)) ) # eo ddf
            # dim(ddf) ; head(ddf)
            rm(mtable); gc()
            # Dcast
            d_ddf <- dcast(ddf, cell_id + x + y ~ variable, fun.aggregate = mean, na.rm = T, value.var = "mean")
            # dim(d_ddf); head(d_ddf)
            # Add richness and Shannon index values while you are at it 
            d_ddf$rich <- rowSums(as.matrix(d_ddf[,c(4:length(d_ddf))]))
            require("vegan")
            d_ddf$H <- diversity(d_ddf[,c(4:(length(d_ddf)-1))], index = "shannon")
            # summary(d_ddf[,c(520:length(d_ddf))])
                        
            message(paste("Saving table and diversity metric for the baseline, for ",m, sep = ""))
            write.table(d_ddf, file = paste("table_",group,"_mean_compo_baseline","_",m,".txt", sep = ""), sep = "\t")
            rm(d_ddf); gc()
            message(paste(" ",sep = ""))
            message(paste(" ",sep = ""))
            
} # eo FUN - community.averager

files <- dir()[grep("mon_composition_baseline",dir())] ; files
require("parallel")
mclapply(files, FUN = community.averager, mc.cores = 25)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------



