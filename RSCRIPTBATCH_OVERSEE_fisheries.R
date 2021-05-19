
# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("raster")
library("viridis")
library("scales")
library("maps")
library("cmocean")
library("RColorBrewer")
library("ggthemes")
library("parallel")

WD <- getwd()
world <- map_data("world")
world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------

### A°) Check catch data files (Ind and NInd)
# setwd(paste(WD,"/","catch_data", sep = ""))
# dir()
# catches <- read.csv("CatchNInd1985_1989.csv", sep = ",", h = T)
# dim(catches) # 9,183,467
# head(catches)
# length(unique(catches$ID)) # 11,478 || so some rows do have same ID
# length(unique(catches$Cell)) # 138,355
# unique(catches$Cell)[1:1000]

### Describe file: 
# ID = seem to be the same ID as in the master index file below, use it to trace back the data
# Cell = cell id to be matched with the code in codes.xlsx ! 
# Reported = large scale + small fisheries catch rate (tons per km^2 in year)
# IUU = Rate of illegal and otherwise unreported catch (same unit)
# Discards = Rate of associated discards at sea (same unit)

# ----------------------------------------------------------

### B°) Check master index an codes files

### NOTE:   if file = IndexInd.csv --> INDUSTRIAL FISHING
###         if file = IndexNInd.csv --> NON INDUSTRIAL FISHING
setwd(WD)
IndexNInd <- read.csv("IndexNInd.csv", sep = ",", h = T)
dim(IndexNInd) # 246189
head(IndexNInd)
length(unique(IndexNInd$ID)) # 246,189 -> one ID per row...maybe the ID in the catch data are to be traced back in this index ID

IndexInd <- read.csv("IndexInd.csv", sep = ",", h = T)
dim(IndexInd) # 42,963
head(IndexInd)
length(unique(IndexInd$ID)) # 42,963 -> one ID per row...maybe the ID in the catch data are to be traced back in this index ID

# And codes files 
cells <- read.csv("codes_cells_Watson2017.csv", sep = ";", h = T)
#head(cells) ; dim(cells) # 149,928 cells
gears <- read.csv("codes_gear_Watson2017.csv", sep = ";", h = T)
countries <- read.csv("codes_countries_Watson2017.csv", sep = ";", h = T)
taxa <- read.csv("codes_taxa_Watson2017.csv", sep = ";", h = T)
# unique(taxa$TaxonName) # 1344 levels
# unique(taxa$CommonName) # 1322

### Fom one ind$ID, you can get:
# IYear -> year
# CNumber -> country number so country in the codes.xlsx
# Taxonkey -> type of animal, varying resolution but soemtimes species, also in code.xlsx
# Gear --> gear used to fish (also in codes.xlsx)
# NumCells --> I guess number of cell fitting the ID (so common year/country/taxon etc)
# Reported = large scale (LSF) + small fisheries (SSF) catch rate (tons per km^2 in year)
### NOTE: LSF >>> SSF so mainly LSF contrib
# IUUTotal = Rate of illegal and otherwise unreported catch (same unit)
# Discards = Rate of associated discards at sea (same unit)

### OK so combining the industrial and/or the non industrial data tabmes (.csv) with codes.xlsx and the master index files, 
### you should be able to get for each ID a 
# Year
# Long/ Lat
# Gear
# Taxonkey
# Issuing country
# Reported catch rate in tons per km2 per year


# ----------------------------------------------------------

### For both industrial and non industrial data, concatenate the data based on the index and the code files
setwd(paste(WD,"/","catch_data", sep = ""))
files.indus <- dir()[grep("CatchInd", dir())]
files.nonindus <- dir()[grep("CatchNInd", dir())]
# Don't really care about the time period because you'll get the year of the data from the "ind" file


### Perform a lapply for each file within which you'll perform a mclapply with like 40 cores to retrieve the data from each ID
# f <- files.indus[6] # for testing
#f <- "CatchInd1995_1999.csv"

lapply(files.indus, function(f) {

            # Message
            message(paste("", sep = ""))
            message(paste("Preparing catch data for ", f, sep = ""))
            message(paste("", sep = ""))
            catch <- read.csv(f, h = T, ",")
            # head(catch) ; dim(catch)
            ### Change of strategy: split the catch data.frame into 50 bits and perform the mclapply on those bits so it dos not swap
            ### Save every bit along the way 
            catch$bit <- cut(1:nrow(catch), 50, labels = F)
            # unique(catch$bit)
            # b <- 9
            for(b in unique(catch$bit)) {
                
                    message(paste("for bit = ", b, " -------------------------------------------- ", sep = ""))
                    catch2 <- catch[which(catch$bit == b),]
                    IDs <- unique(catch2$ID) # id event, from 'IndexInd' object
                    commons <- intersect(unique(catch2$ID), unique(IndexInd$ID))
                    # length(commons); length(IDs)
                    
                    require("parallel")
                    res <- mclapply(IDs, function(i) {
                
                                message(paste("# ", which(IDs == i), " - ",i," || ", round((which(IDs == i)/length(IDs)),4)*100, "%", sep = ""))
                                # head( catch[catch$ID == i,] ) # dim(catch[catch$ID == i,])        
                                catch3 <- catch2[catch2$ID == i,]
                                # dim(catch2)
                        
                                ### Supply metadata (long/lat/year/taxa etc.)
                                # Supply Long & Lat from 'cells' using the 'Cell' column for catch2
                                catch3$Long <- NA; catch3$Lat <- NA; catch3$area <- NA
                                for(c in unique(catch3$Cell)) {
                                        # c <- unique(catch2$Cell)[3] # for testing
                                        # cells[cells$Cell == c,]
                                        catch3[which(catch3$Cell == c),"Long"] <- cells[cells$Cell == c,"Lon"]
                                        catch3[which(catch3$Cell == c),"Lat"] <- cells[cells$Cell == c,"Lat"]
                                        catch3[which(catch3$Cell == c),"area"] <- cells[cells$Cell == c,"Areasqkm"]
                                } # eo for 
                            
                                # Supply Year from IndexInd
                                if( isTRUE(i %in% IndexInd$ID) ) {
                            
                                    catch3$Gear <- NA
                                    catch3$TaxonName <- NA
                                    catch3$CommonName <- NA
                                    catch3$Year <- NA
                            
                                    catch3$Year <- unique(IndexInd[IndexInd$ID == i,"IYear"])[1]
                                    Taxonkey <-  unique(IndexInd[IndexInd$ID == i,"Taxonkey"])[1]
                                    Gearkey <-  unique(IndexInd[IndexInd$ID == i,"Gear"])[1]
                            
                                    catch3$Gear <-  unique(gears[gears$Gear == Gearkey,"VBDesc"])[1]
                                    catch3$TaxonName <-  unique(taxa[taxa$TaxonKey == Taxonkey,"TaxonName"])[1]
                                    catch3$CommonName <-  unique(taxa[taxa$TaxonKey == Taxonkey,"CommonName"])[1]
                            
                                } else {
                            
                                    catch3$Gear <- NA
                                    catch3$TaxonName <- NA
                                    catch3$CommonName <- NA
                                    catch3$Year <- NA
                            
                                } # if else loop
                                # Check
                                return(catch3)
                        
                        }, mc.cores = 43
            
                    ) # eo mclapply - i in IDs
                    # Rbind
                    #t <- dplyr::bind_rows(res)
                    t <- do.call(rbind, res)
                    #str(t); head(t); dim(t)
                    # Find the problematic obs
                    #dim( t[t$Lat == unique(t$Lat)[grep("Error",unique(t$Lat))],] )
                    # t[t$Lat == unique(t$Lat)[grep("Error",unique(t$Lat))],"ID"] 
                    tt <- t[-which(t$Lat == unique(t$Lat)[grep("Error",unique(t$Lat))]),]
                    tt$ID <- as.integer(tt$ID)
                    tt$Cell <- as.integer(tt$Cell)
                    tt$Reported <- as.numeric(tt$Reported)
                    tt$Discards <- as.numeric(tt$Discards)
                    tt$IUU <- as.numeric(tt$IUU)
                    tt$bit <- as.integer(tt$bit)
                    tt$Long <- as.numeric(tt$Long)
                    tt$Lat <- as.numeric(tt$Lat)
                    tt$area <- as.numeric(tt$area)

                    rm(res) ; gc()
                    # comms <- intersect(unique(t$ID), unique(catch$ID)) # length(comms)
                    # Who's missing? 
                    # unique(catch$ID)[!(unique(catch$ID) %in% comms)]
                    
                    ### Save 
                    filename <- paste(str_replace(f,".csv",""),"_bit_",b,".Rdata", sep = "")
                    message(paste("Saving catch data for ", filename, sep = ""))
                    save(tt, file = filename)
                
            } # eo for loop - b in unique(catch$bit)

        } # EO FUN

) # eo 

# ----------------------------------------------------------

### Same but with NON industrial
setwd(paste(WD,"/","catch_data", sep = ""))
# f <- res.nonindus[5]
lapply(files.nonindus, function(f) {

            # Message
            message(paste("", sep = ""))
            message(paste("Preparing catch data for ", f, sep = ""))
            message(paste("", sep = ""))
            catch <- read.csv(f, h = T, ",")
            # head(catch) ; dim(catch)
            ### Change of strategy: split the catch data.frame into 50 bits and perform the mclapply on those bits so it dos not swap
            ### Save every bit along the way 
            catch$bit <- cut(1:nrow(catch), 50, labels = F)
            # unique(catch$bit)
            
            for(b in unique(catch$bit)) {
                
                    message(paste("for bit = ", b, " -------------------------------------------- ", sep = ""))
                    catch2 <- catch[which(catch$bit == b),]
                    IDs <- unique(catch2$ID) # id event, from 'IndexInd' object
                    commons <- intersect(unique(catch2$ID), unique(IndexInd$ID))
                    # length(commons); length(IDs)
                    
                    require("parallel")
                    # i <- IDs[100]
                    # i <- commons[1]
                    # i <- 8343
                    res <- mclapply(IDs, function(i) {
                
                                message(paste("# ", which(IDs == i), " - ",i," || ", round((which(IDs == i)/length(IDs)),4)*100, "%", sep = ""))
                                # head( catch[catch$ID == i,] ) # dim(catch[catch$ID == i,])        
                                catch3 <- catch2[catch2$ID == i,]
                                # dim(catch2)
                        
                                ### Supply metadata (long/lat/year/taxa etc.)
                                # Supply Long & Lat from 'cells' using the 'Cell' column for catch2
                                catch3$Long <- NA; catch3$Lat <- NA; catch3$area <- NA
                                for(c in unique(catch3$Cell)) {
                                        # c <- unique(catch2$Cell)[3] # for testing
                                        # cells[cells$Cell == c,]
                                        catch3[which(catch3$Cell == c),"Long"] <- cells[cells$Cell == c,"Lon"]
                                        catch3[which(catch3$Cell == c),"Lat"] <- cells[cells$Cell == c,"Lat"]
                                        catch3[which(catch3$Cell == c),"area"] <- cells[cells$Cell == c,"Areasqkm"]
                                } # eo for 
                        
                                # Supply Year from IndexInd
                                if( isTRUE(i %in% IndexInd$ID) ) {
                            
                                    catch3$Gear <- NA
                                    catch3$TaxonName <- NA
                                    catch3$CommonName <- NA
                                    catch3$Year <- NA
                            
                                    catch3$Year <- unique(IndexInd[IndexInd$ID == i,"IYear"])[1]
                                    Taxonkey <-  unique(IndexInd[IndexInd$ID == i,"Taxonkey"])[1]
                                    Gearkey <-  unique(IndexInd[IndexInd$ID == i,"Gear"])[1]
                            
                                    catch3$Gear <-  unique(gears[gears$Gear == Gearkey,"VBDesc"])[1]
                                    catch3$TaxonName <-  unique(taxa[taxa$TaxonKey == Taxonkey,"TaxonName"])[1]
                                    catch3$CommonName <-  unique(taxa[taxa$TaxonKey == Taxonkey,"CommonName"])[1]
                            
                                } else {
                            
                                    catch3$Gear <- NA
                                    catch3$TaxonName <- NA
                                    catch3$CommonName <- NA
                                    catch3$Year <- NA
                            
                                } # if else loop
                                # Check
                                return(catch3)
                        
                        }, mc.cores = 43
            
                    ) # eo mclapply - i in IDs
                    # Rbind
                    t <- do.call(rbind, res)
                    rm(res) ; gc()
                    
                    ### Save 
                    filename <- paste(str_replace(f,".csv",""),"_bit_",b,".Rdata", sep = "")
                    message(paste("Saving catch data for ", filename, sep = ""))
                    save(t, file = filename)
                
            } # eo for loop - b in unique(catch$bit)

        } # EO FUN

) # eo 

# ----------------------------------------------------------

### 14/04/2020: Gather all data and compute annual climatologies
setwd(paste(WD,"/","catch_data", sep = ""))
files.indus <- dir()[grep("CatchInd", dir())] # files.indus
files.indus <- files.indus[grep(".Rdata", files.indus)]
files.nonindus <- dir()[grep("CatchNInd", dir())] # files.indus
files.nonindus <- files.nonindus[grep(".Rdata", files.nonindus)]

# Use mclapply
#f <- "CatchInd1995_1999_bit_10.Rdata"
res.indus <- mclapply(files.indus, function(f) {
                d <- get(load(f))
                # str(d)
                # d$Long <- as.numeric(d$Long)
                message(paste(f, sep = ""))
                if( is.numeric(d$Long) ) { 
                    message("Numerics are numerics - OK")
                } else {
                    message(" ------------------------------------------ NOT NUMERIC ------------------------------------------ ")
                }
                d$Year <- factor(d$Year)
                return(d)
    }, mc.cores = 40
) # eo mclapply - res.indus
# Rbind
indus <- dplyr::bind_rows(res.indus)
dim(indus); head(indus)
# 236,993,519 observations
str(indus)
rm(res.indus) ; gc()

### For non indus 
res.nonindus <- mclapply(files.nonindus, function(f) {
                d <- get(load(f))
                # str(d)
                #d$Long <- as.numeric(d$Long)
                message(paste(f, sep = ""))
                if( is.numeric(d$Long) ) { 
                    message("Numerics are numerics - OK")
                } else {
                    message(" ------------------------------------------ NOT NUMERIC ------------------------------------------ ")
                }
                d$Year <- factor(d$Year)
                return(d)
            }, mc.cores = 40
) # eo mclapply - res.indus
# Rbind
non <- dplyr::bind_rows(res.nonindus)
dim(non) # 30,412,706 observations
head(non)
str(non)
summary(non)
rm(res.nonindus) ; gc()

indus$Year <- as.numeric(indus$Year)
non$Year <- as.numeric(non$Year)

### Combine non industrial and industrial
all_catches <- rbind(indus, non)
dim(all_catches) # 267,406,225 observations
str(all_catches)

### Save as master table so you don't have to re-load all the other files "_bit_.Rdata"
save(all_catches, file = "table_all_catches_1985-2015.Rdata")

# ----------------------------------------------------------

library("tidyverse")
library("reshape2")
library("raster")
library("viridis")
library("scales")
library("maps")
library("cmocean")
library("RColorBrewer")
library("ggthemes")
library("parallel")

WD <- getwd()
world <- map_data("world")
world2 <- map_data("world2")

### 15/04/2020: Making sure "table_all_catches_1985-2015.Rdata" match the maps showed by Reg Watson in his papers
setwd(paste(WD,"/","catch_data", sep = ""))
data <- get(load("table_all_catches_1985-2015.Rdata"))
dim(data) # 267406225 obs
str(data)

### Check Years
summary(factor(data$Year)) # 91% of data don't have year ...

### Assess distrbution of catch da
#plot <- ggplot(data, aes(y = Reported)) + geom_boxplot(fill = "grey50", colour = "black") + theme_classic()
#ggsave(plot = plot, filename = "boxplot_distrb_reported.jpg", dpi = 300, height = 4, width = 4)
#plot <- ggplot(data, aes(y = IUU)) + geom_boxplot(fill = "grey50", colour = "black") + theme_classic()
#ggsave(plot = plot, filename = "boxplot_distrb_unreported.jpg", dpi = 300, height = 4, width = 4)
### Plotting takes an extremely long time because of the nb of data...do so after computing sums and averages

### Compute ensemble mean after creating a cell id (unless the current one already works?)
data$id <- factor(paste(data$Long, data$Lat, sep = "_"))
length(unique(data$id)) # 146,830 cells, which corresponds to 98% of all marine cells in the 'cells' data frame
length(unique(data$Cell)) # ok makes sense

require("dplyr")
ddf <- data.frame(data %>%
    group_by(id) %>%
    summarize(x = unique(Long), y = unique(Lat),
            mean_report = mean(Reported, na.rm = T),
            mean_iuu = mean(IUU, na.rm = T),
            sum_report = sum(Reported),
            sum_iuu = sum(IUU) 
        )
) # eo ddf
# Recompute another mean by dividing the "sum" col by the number of years considered in the time period 1985:2015 -> 31
ddf$mean2_report <- ddf$sum_report/31
ddf$mean2_iuu <- ddf$mean_iuu/31
# And combine in reported cacthes and iuu catches in tot: 
ddf$tot_mean <- (ddf$mean_report)+(ddf$mean_iuu)
ddf$tot_mean2 <- (ddf$mean2_report)+(ddf$mean2_iuu)
ddf$tot_sum <- (ddf$sum_report)+(ddf$sum_iuu)
# Check 
summary(ddf)
### Save clim
save(ddf, file = "clim_fisheries_all_1985-2015.Rdata")

### Assess distrbution of catch data (reported + unreported)
ggplot(ddf[ddf$tot_mean < 100,], aes(x = tot_mean)) + geom_histogram(colour = "black", fill = "white") + 
    geom_vline(aes(xintercept = mean(tot_mean, na.rm = T)), color = "red", linetype = "dashed")
#
ggplot(ddf[ddf$tot_mean2 < 500,], aes(x = tot_mean2)) + geom_histogram(colour = "black", fill = "white") + 
    geom_vline(aes(xintercept = mean(tot_mean2, na.rm = T)), color = "red", linetype = "dashed")
#       
ggplot(ddf[ddf$tot_sum < 500,], aes(x = tot_sum)) + geom_histogram(colour = "black", fill = "white") + 
    geom_vline(aes(xintercept = mean(tot_sum, na.rm = T)), color = "red", linetype = "dashed")

### Apply log transformation to assess which wold be better suited for our analyses        
ggplot(ddf, aes(x = log1p(tot_mean))) + geom_histogram(colour = "black", fill = "white") + 
    geom_vline(aes(xintercept = mean(log1p(tot_mean), na.rm = T)), color = "red", linetype = "dashed")
#
ggplot(ddf, aes(x = log1p(tot_mean2))) + geom_histogram(colour = "black", fill = "white") + 
    geom_vline(aes(xintercept = mean(log1p(tot_mean2), na.rm = T)), color = "red", linetype = "dashed")
#       
ggplot(ddf, aes(x = log1p(tot_sum))) + geom_histogram(colour = "black", fill = "white") + 
    geom_vline(aes(xintercept = mean(log1p(tot_sum), na.rm = T)), color = "red", linetype = "dashed")
    

### Check correlation between the two averages    
cor(ddf$tot_mean, ddf$tot_mean2, method = "spearman") # 0.87 ok 
summary(lm(tot_mean2 ~ tot_mean, data = ddf)) # Adjusted R-squared: 0.3928 (not great)

### Create a discrete scale that matches the maps of Watson et al. (Spectral color palette with 7 categories)
summary(ddf$tot)
summary(ddf$tot2) # looks better for the tons/year scale 

### For tons/year (Link & Watson, 2019)
# - 0-2 
# - 2-10
# - 10-20
# - 20-40
# - 40-60
# - 60-100
# - >100

ddf$tot_mean_discr <- NA
ddf[which(ddf$tot_mean < 2),"tot_mean_discr"] <- "0-2"
ddf[which(ddf$tot_mean >= 2 & ddf$tot_mean < 10),"tot_mean_discr"] <- "2-10"
ddf[which(ddf$tot_mean >= 10 & ddf$tot_mean < 20),"tot_mean_discr"] <- "10-20"
ddf[which(ddf$tot_mean >= 20 & ddf$tot_mean < 40),"tot_mean_discr"] <- "20-40"
ddf[which(ddf$tot_mean >= 40 & ddf$tot_mean < 60),"tot_mean_discr"] <- "40-60"
ddf[which(ddf$tot_mean >= 60 & ddf$tot_mean < 100),"tot_mean_discr"] <- "60-100"
ddf[which(ddf$tot_mean >= 100),"tot_mean_discr"] <- ">100"

ddf$tot_mean2_discr <- NA
ddf[which(ddf$tot_mean2 < 2),"tot_mean2_discr"] <- "0-2"
ddf[which(ddf$tot_mean2 >= 2 & ddf$tot_mean2 < 10),"tot_mean2_discr"] <- "2-10"
ddf[which(ddf$tot_mean2 >= 10 & ddf$tot_mean2 < 20),"tot_mean2_discr"] <- "10-20"
ddf[which(ddf$tot_mean2 >= 20 & ddf$tot_mean2 < 40),"tot_mean2_discr"] <- "20-40"
ddf[which(ddf$tot_mean2 >= 40 & ddf$tot_mean2 < 60),"tot_mean2_discr"] <- "40-60"
ddf[which(ddf$tot_mean2 >= 60 & ddf$tot_mean2 < 100),"tot_mean2_discr"] <- "60-100"
ddf[which(ddf$tot_mean2 >= 100),"tot_mean2_discr"] <- ">100"

ddf$tot_sum_discr <- NA
ddf[which(ddf$tot_sum < 2),"tot_sum_discr"] <- "0-2"
ddf[which(ddf$tot_sum >= 2 & ddf$tot_sum < 10),"tot_sum_discr"] <- "2-10"
ddf[which(ddf$tot_sum >= 10 & ddf$tot_sum < 20),"tot_sum_discr"] <- "10-20"
ddf[which(ddf$tot_sum >= 20 & ddf$tot_sum < 40),"tot_sum_discr"] <- "20-40"
ddf[which(ddf$tot_sum >= 40 & ddf$tot_sum < 60),"tot_sum_discr"] <- "40-60"
ddf[which(ddf$tot_sum >= 60 & ddf$tot_sum < 100),"tot_sum_discr"] <- "60-100"
ddf[which(ddf$tot_sum >= 100),"tot_sum_discr"] <- ">100"

summary(factor(ddf$tot_mean_discr))
summary(factor(ddf$tot_mean2_discr))
summary(factor(ddf$tot_sum_discr))

### Map just to make sure the patterns make sense
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean_discr)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean2_discr)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_sum_discr)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

# ggsave
ggsave(plot = map1, filename = "map_tot_mean_discr_1985-2015.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map2, filename = "map_tot_mean2_discr_1985-2015.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map3, filename = "map_tot_sum_discr_1985-2015.jpg", dpi = 300, height = 5, width = 7)


### Continuous colorbar with logged data
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(tot_mean)), data = ddf) +
   	scale_fill_viridis(name = "Annual\nmean catches\nlog(t/km2.yr)", limits = c(0,4)) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(tot_mean2)), data = ddf) +
   	scale_fill_viridis(name = "Annual\nmean catches\nlog(t/km2.yr)", limits = c(0,6.5)) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(tot_sum)), data = ddf) +
   	scale_fill_viridis(name = "Total\nmean catches\nlog(t/km2.yr)") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

# ggsave
ggsave(plot = map1, filename = "map_tot_mean_log_1985-2015.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map2, filename = "map_tot_mean2_log_1985-2015.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map3, filename = "map_tot_sum_log_1985-2015.jpg", dpi = 300, height = 5, width = 7)

### And check a continuous color palette for non logged mean2
ggplot() + geom_raster(aes(x = x, y = y, fill = tot_mean2), data = ddf) +
   	scale_fill_viridis(name = "Annual\nmean catches\n(t/km2.yr)", limits = c(0,500)) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    




# ----------------------------------------------------------

### Re-compute a climatology but based on more recent records like in Watson's papers to assess if you're biased by the lower catches that occurred over 1985-1990

### DO the same for 2010-2014 too

setwd(paste(WD,"/","catch_data", sep = ""))
files.indus <- dir()[grep("CatchInd", dir())] # files.indus
files.indus <- files.indus[grep(".Rdata", files.indus)]
# Restrict to 2010-2014 (files 101 to 350)
files.indus <- files.indus[251:300]

files.nonindus <- dir()[grep("CatchNInd", dir())] # files.indus
files.nonindus <- files.nonindus[grep(".Rdata", files.nonindus)]
# Restrict to 2010-2014 (files 101 to 350)
files.nonindus <- files.nonindus[251:300]

# Use mclapply
res.indus <- mclapply(files.indus, function(f) {
                d <- get(load(f))
                # str(d)
                # d$Long <- as.numeric(d$Long)
                message(paste(f, sep = ""))
                if( is.numeric(d$Long) ) { 
                    message("Numerics are numerics - OK")
                } else {
                    message(" ------------------------------------------ NOT NUMERIC ------------------------------------------ ")
                }
                d$Year <- factor(d$Year)
                return(d)
    }, mc.cores = 40
) # eo mclapply - res.indus
# Rbind
indus <- dplyr::bind_rows(res.indus)
rm(res.indus) ; gc()
dim(indus); head(indus)
# 45,754,091 observations
str(indus)
summary(indus)
summary(factor(indus$Year))

### For non indus 
res.nonindus <- mclapply(files.nonindus, function(f) {
                d <- get(load(f))
                # str(d)
                #d$Long <- as.numeric(d$Long)
                message(paste(f, sep = ""))
                if( is.numeric(d$Long) ) { 
                    message("Numerics are numerics - OK")
                } else {
                    message(" ------------------------------------------ NOT NUMERIC ------------------------------------------ ")
                }
                d$Year <- factor(d$Year)
                return(d)
            }, mc.cores = 40
) # eo mclapply - res.indus
# Rbind
non <- dplyr::bind_rows(res.nonindus)
rm(res.nonindus) ; gc()
dim(non) # 21,717,528 observations
head(non)
str(non)
summary(non)

indus$Year <- as.numeric(indus$Year)
non$Year <- as.numeric(non$Year)

### Combine non industrial and industrial
all_catches <- rbind(indus, non)
dim(all_catches) # 198,384,380 observations
str(all_catches)

### Save as master table so you don't have to re-load all the other files "_bit_.Rdata"
save(all_catches, file = "table_all_catches_2010-2014.Rdata")

### Compute ensemble mean after creating a cell id (unless the current one already works?)
all_catches$id <- factor(paste(all_catches$Long, all_catches$Lat, sep = "_"))
# length(unique(all_catches$id)) # 146,830 cells, which corresponds to 98% of all marine cells in the 'cells' data frame
# length(unique(all_catches$Cell)) # ok makes sense

require("dplyr")
ddf <- data.frame(all_catches %>%
    group_by(id) %>%
    summarize(x = unique(Long), y = unique(Lat),
        mean_report = mean(Reported, na.rm = T),
        mean_iuu = mean(IUU, na.rm = T),
        sum_report = sum(Reported),
        sum_iuu = sum(IUU) 
        )
) # eo ddf
# Check 
summary(ddf)

# ANd re-compute another mean by dividing the "sum" col by the number of years considered in the time period 1995:2015 -> 21
ddf$mean2_report <- ddf$sum_report/4
ddf$mean2_iuu <- ddf$mean_iuu/4

# And combine in reported cacthes and iuu catches in tot: 
ddf$tot_mean <- (ddf$mean_report)+(ddf$mean_iuu)
ddf$tot_mean2 <- (ddf$mean2_report)+(ddf$mean2_iuu)
ddf$tot_sum <- (ddf$sum_report)+(ddf$sum_iuu)

summary(ddf)

### Save clim file
save(ddf, file = "clim_fisheries_all_2010-2014.Rdata")

### Create a discrete color palette that matches Watson's maps
ddf$tot_mean_discr <- NA
ddf[which(ddf$tot_mean < 2),"tot_mean_discr"] <- "0-2"
ddf[which(ddf$tot_mean >= 2 & ddf$tot_mean < 10),"tot_mean_discr"] <- "2-10"
ddf[which(ddf$tot_mean >= 10 & ddf$tot_mean < 20),"tot_mean_discr"] <- "10-20"
ddf[which(ddf$tot_mean >= 20 & ddf$tot_mean < 40),"tot_mean_discr"] <- "20-40"
ddf[which(ddf$tot_mean >= 40 & ddf$tot_mean < 60),"tot_mean_discr"] <- "40-60"
ddf[which(ddf$tot_mean >= 60 & ddf$tot_mean < 100),"tot_mean_discr"] <- "60-100"
ddf[which(ddf$tot_mean >= 100),"tot_mean_discr"] <- ">100"

ddf$tot_mean2_discr <- NA
ddf[which(ddf$tot_mean2 < 2),"tot_mean2_discr"] <- "0-2"
ddf[which(ddf$tot_mean2 >= 2 & ddf$tot_mean2 < 10),"tot_mean2_discr"] <- "2-10"
ddf[which(ddf$tot_mean2 >= 10 & ddf$tot_mean2 < 20),"tot_mean2_discr"] <- "10-20"
ddf[which(ddf$tot_mean2 >= 20 & ddf$tot_mean2 < 40),"tot_mean2_discr"] <- "20-40"
ddf[which(ddf$tot_mean2 >= 40 & ddf$tot_mean2 < 60),"tot_mean2_discr"] <- "40-60"
ddf[which(ddf$tot_mean2 >= 60 & ddf$tot_mean2 < 100),"tot_mean2_discr"] <- "60-100"
ddf[which(ddf$tot_mean2 >= 100),"tot_mean2_discr"] <- ">100"

ddf$tot_sum_discr <- NA
ddf[which(ddf$tot_sum < 2),"tot_sum_discr"] <- "0-2"
ddf[which(ddf$tot_sum >= 2 & ddf$tot_sum < 10),"tot_sum_discr"] <- "2-10"
ddf[which(ddf$tot_sum >= 10 & ddf$tot_sum < 20),"tot_sum_discr"] <- "10-20"
ddf[which(ddf$tot_sum >= 20 & ddf$tot_sum < 40),"tot_sum_discr"] <- "20-40"
ddf[which(ddf$tot_sum >= 40 & ddf$tot_sum < 60),"tot_sum_discr"] <- "40-60"
ddf[which(ddf$tot_sum >= 60 & ddf$tot_sum < 100),"tot_sum_discr"] <- "60-100"
ddf[which(ddf$tot_sum >= 100),"tot_sum_discr"] <- ">100"

summary(factor(ddf$tot_mean_discr))
summary(factor(ddf$tot_mean2_discr))
summary(factor(ddf$tot_sum_discr))

### Map just to make sure the patterns make sense
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean_discr)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean2_discr)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_sum_discr)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

# ggsave
ggsave(plot = map1, filename = "map_tot_mean_discr_2010-2014.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map2, filename = "map_tot_mean2_discr_2010-2014.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map3, filename = "map_tot_sum_discr_2010-2014.jpg", dpi = 300, height = 5, width = 7)

### Try replotting tot_mean_discr in kg/km2.year (Watson, 2017)
ddf$tot_mean_kg <- (ddf$tot_mean)*1000
ddf$tot_mean_discr_kg <- NA
ddf[which(ddf$tot_mean_kg < 0.1),"tot_mean_discr_kg"] <- "0-0.1"
ddf[which(ddf$tot_mean_kg >= 0.1 & ddf$tot_mean_kg < 1),"tot_mean_discr_kg"] <- "0.1-1"
ddf[which(ddf$tot_mean_kg >= 1 & ddf$tot_mean_kg < 2),"tot_mean_discr_kg"] <- "1-2"
ddf[which(ddf$tot_mean_kg >= 2 & ddf$tot_mean_kg < 10),"tot_mean_discr_kg"] <- "2-10"
ddf[which(ddf$tot_mean_kg >= 10 & ddf$tot_mean_kg < 50),"tot_mean_discr_kg"] <- "10-50"
ddf[which(ddf$tot_mean_kg >= 50 & ddf$tot_mean_kg < 200),"tot_mean_discr_kg"] <- "50-200"
ddf[which(ddf$tot_mean_kg >= 200),"tot_mean_discr_kg"] <- ">200"

ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean_discr_kg)), data = ddf) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-0.1","0.1-1","1-2","2-10","10-50","50-200",">200"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    
            
### That does not seem to be it...
### Reduce the spatial resilution to 1°x1° and recompute climatology. Presently you have a 0.5°x0.5°
unique(round(all_catches$Long/0.5)*0.5)
unique(round(all_catches$Lat/0.5)*0.5)
all_catches$id2 <- factor(paste(round(all_catches$Long/0.5)*0.5, round(all_catches$Lat/0.5)*0.5, sep = "_"))
length(unique(all_catches$id2)) # 38450 so 1°x1°...

require("dplyr")
ddf2 <- data.frame(all_catches %>%
    group_by(id2) %>%
    summarize(x = unique(round(Long/0.5)*0.5), y = unique(round(Lat/0.5)*0.5),
        mean_report = mean(Reported, na.rm = T),
        mean_iuu = mean(IUU, na.rm = T),
        sum_report = sum(Reported),
        sum_iuu = sum(IUU) 
        )
) # eo ddf
# Check 
summary(ddf2)

# Recompute another mean by dividing the "sum" col by the number of years considered in the time period 1995:2015 -> 21
ddf2$mean2_report <- ddf2$sum_report/4
ddf2$mean2_iuu <- ddf2$mean_iuu/4

# And combine in reported cacthes and iuu catches in tot: 
ddf2$tot_mean <- (ddf2$mean_report)+(ddf2$mean_iuu)
ddf2$tot_mean2 <- (ddf2$mean2_report)+(ddf2$mean2_iuu)
ddf2$tot_sum <- (ddf2$sum_report)+(ddf2$sum_iuu)

### Create a discrete color palette that matches Watson's maps
ddf2$tot_mean_discr <- NA
ddf2[which(ddf2$tot_mean < 2),"tot_mean_discr"] <- "0-2"
ddf2[which(ddf2$tot_mean >= 2 & ddf2$tot_mean < 10),"tot_mean_discr"] <- "2-10"
ddf2[which(ddf2$tot_mean >= 10 & ddf2$tot_mean < 20),"tot_mean_discr"] <- "10-20"
ddf2[which(ddf2$tot_mean >= 20 & ddf2$tot_mean < 40),"tot_mean_discr"] <- "20-40"
ddf2[which(ddf2$tot_mean >= 40 & ddf2$tot_mean < 60),"tot_mean_discr"] <- "40-60"
ddf2[which(ddf2$tot_mean >= 60 & ddf2$tot_mean < 100),"tot_mean_discr"] <- "60-100"
ddf2[which(ddf2$tot_mean >= 100),"tot_mean_discr"] <- ">100"

ddf2$tot_mean2_discr <- NA
ddf2[which(ddf2$tot_mean2 < 2),"tot_mean2_discr"] <- "0-2"
ddf2[which(ddf2$tot_mean2 >= 2 & ddf2$tot_mean2 < 10),"tot_mean2_discr"] <- "2-10"
ddf2[which(ddf2$tot_mean2 >= 10 & ddf2$tot_mean2 < 20),"tot_mean2_discr"] <- "10-20"
ddf2[which(ddf2$tot_mean2 >= 20 & ddf2$tot_mean2 < 40),"tot_mean2_discr"] <- "20-40"
ddf2[which(ddf2$tot_mean2 >= 40 & ddf2$tot_mean2 < 60),"tot_mean2_discr"] <- "40-60"
ddf2[which(ddf2$tot_mean2 >= 60 & ddf2$tot_mean2 < 100),"tot_mean2_discr"] <- "60-100"
ddf2[which(ddf2$tot_mean2 >= 100),"tot_mean2_discr"] <- ">100"

ddf2$tot_sum_discr <- NA
ddf2[which(ddf2$tot_sum < 2),"tot_sum_discr"] <- "0-2"
ddf2[which(ddf2$tot_sum >= 2 & ddf2$tot_sum < 10),"tot_sum_discr"] <- "2-10"
ddf2[which(ddf2$tot_sum >= 10 & ddf2$tot_sum < 20),"tot_sum_discr"] <- "10-20"
ddf2[which(ddf2$tot_sum >= 20 & ddf2$tot_sum < 40),"tot_sum_discr"] <- "20-40"
ddf2[which(ddf2$tot_sum >= 40 & ddf2$tot_sum < 60),"tot_sum_discr"] <- "40-60"
ddf2[which(ddf2$tot_sum >= 60 & ddf2$tot_sum < 100),"tot_sum_discr"] <- "60-100"
ddf2[which(ddf2$tot_sum >= 100),"tot_sum_discr"] <- ">100"

summary(factor(ddf2$tot_mean_discr))
summary(factor(ddf2$tot_mean2_discr))
summary(factor(ddf2$tot_sum_discr))

### Map just to make sure the patterns make sense
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean_discr)), data = ddf2) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_mean2_discr)), data = ddf2) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(tot_sum_discr)), data = ddf2) +
   	scale_fill_manual(name = "Annual total\nmean catches\n(ton/km2/year)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    

# ggsave
ggsave(plot = map1, filename = "map_tot_mean_discr_2010-2014_1d.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map2, filename = "map_tot_mean2_discr_2010-2014_1d.jpg", dpi = 300, height = 5, width = 7)
ggsave(plot = map3, filename = "map_tot_sum_discr_2010-2014_1d.jpg", dpi = 300, height = 5, width = 7)

### And with continuous palette?
ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(tot_mean2)), data = ddf2) +
   	scale_fill_distiller(name = "Annual mean\n catches\nlog(t/km2.yr)", palette = "Spectral") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )    
#
ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(tot_sum)), data = ddf2) +
   	scale_fill_distiller(name = "Total mean\n catches\nlog(t/km2.yr)", palette = "Spectral") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )  
            
            

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

