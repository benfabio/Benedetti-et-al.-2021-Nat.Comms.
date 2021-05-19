
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

### First, load the codes containing the taxa names, gears and cells coordinates

cells <- read.csv("codes_cells_Watson17.csv", sep = ";", h = T, dec = ",")
gears <- read.csv("codes_gears_Watson17.csv", sep = ";", h = T)
taxa <- read.csv("codes_taxa_Watson17.csv", sep = ";", h = T)
# Check their str
# str(cells) ; summary(cells) ; dim(cells)
# str(gears)
# str(taxa)
# unique(taxa$TaxonName) # 1356 levels
# unique(taxa$CommonName) # 1331

# ----------------------------------------------------------

### For both industrial and non industrial data, concatenate the data based on the index and the code files
setwd("/net/kryo/work/updata/fisheries_watson2017")
files <- dir()[grep("Catch",dir())][9:13]

### Perform a lapply for each file within which you'll perform a mclapply with like 40 cores to retrieve the data from each ID
# f <- files[5] # for testing

lapply(files, function(f) {

            # Message
            setwd("/net/kryo/work/updata/fisheries_watson2017")
            message(paste("", sep = ""))
            message(paste("Preparing catch data for ", f, sep = ""))
            message(paste("", sep = ""))
            catch <- read.csv(f, h = T, ",")
            # head(catch) ; dim(catch) ; str(catch)
            ### Change of strategy: split the catch data.frame into 50 bits and perform the mclapply on those bits so it dos not swap
            ### Save every bit along the way 
            catch$bit <- cut(1:nrow(catch), 43, labels = F)
            # b <- 20
            require("parallel")
            res <- mclapply(unique(catch$bit), function(b) {
                
                        message(paste("for bit = ", b, " -------------------------------------------- ", sep = ""))
                        catch3 <- catch[which(catch$bit == b),]
                        # dim(catch3)
                        
                        # Supply Long & Lat from 'cells' using the 'Cell' column for catch2
                        catch3$Long <- NA ; catch3$Lat <- NA ; catch3$area <- NA
                        catch3$GearUsed <- NA ; catch3$TaxonName <- NA ; catch3$CommonName <- NA
                        # provide ina  for loop based on unique values (because several rows of 'catch3' can have the same codes)
                        message(paste("providing long & lat", sep = ""))
                        # Find cells keys that are common to catch3 AND cells$Cell (index of cell key)
                        commons <- intersect(unique(catch3$Cell), unique(cells$Cell))
                        
                        for(c in commons) {
                            
                            # c <- unique(catch3$Cell)[4]
                            message(paste("# ", which(commons == c), " - ",c," || ", round((which(commons == c)/length(commons)),4)*100, "%", sep = ""))   
                            
                            if( isTRUE(c %in% cells$Cell) ) {
                                catch3[catch3$Cell == c,"Long"] <- unique(cells[cells$Cell == c,"Long"])
                                catch3[catch3$Cell == c,"Lat"] <- unique(cells[cells$Cell == c,"Lat"])
                                catch3[catch3$Cell == c,"area"] <- unique(cells[cells$Cell == c,"area"])
                            } else {
                                message(paste("Cell key was not macthed in cells' index",sep=""))
                            }
                            
                        } # eo for loop - c in Cell
                        
                        # Same with gear code 
                        message(paste("providing fishing gear used", sep = ""))
                        for(g in unique(catch3$Gear) ) {
                            
                            message(paste("# ", which(unique(catch3$Gear) == g), " - ",g," || ", round((which(unique(catch3$Gear) == g)/length(unique(catch3$Gear))),4)*100, "%", sep = ""))   
                            
                            if( isTRUE(g %in% gears$Gear) ) {
                                catch3[catch3$Gear == g,"GearUsed"] <- unique(gears[gears$Gear == g,"FleetGearName"])
                            } else {
                                message(paste("Gear key was not macthed in gears' index",sep=""))
                            }
                            
                        } # eo for loop - g in Gear
                        
                        # And again same with taxonkeys
                        message(paste("providing taxon name", sep = ""))
                        for(t in unique(catch3$Taxonkey) ) {
                            
                            message(paste("# ", which(unique(catch3$Taxonkey) == t), " - ",t," || ", round((which(unique(catch3$Taxonkey) == t)/length(unique(catch3$Taxonkey))),4)*100, "%", sep = ""))   
                            
                            if( isTRUE(t %in% taxa$Taxonkey) ) {
                                catch3[catch3$Taxonkey == t,"TaxonName"] <- as.character(unique(taxa[taxa$Taxonkey == t,"TaxonName"]))
                                catch3[catch3$Taxonkey == t,"CommonName"] <- as.character(unique(taxa[taxa$Taxonkey == t,"CommonName"]))    
                            } else {
                                message(paste("Taxon key was not macthed in taxa' index",sep=""))
                            }
                            
                        } # eo for loop - t in Taxonkey
                        # Check
                        # str(catch3) ; summary(catch3)
                        # head( catch3[is.na(catch3$Lat),] )
                        # Return
                        return(catch3)
                
                }, mc.cores = 43

            ) # eo mclapply - b in bits
            # Rbind
            t <- dplyr::bind_rows(res)
            # str(t); head(t); dim(t)
            rm(res) ; gc()
            
            ### Save 
            setwd(paste(WD,"/","catch_data", sep = ""))
            filename <- paste(str_replace(f,".csv",""),"_treated_16_04_10",".Rdata", sep = "") # filename
            message(paste("Saving catch data for ", filename, sep = ""))
            save(t, file = filename)
            rm(t) ; gc()
            setwd("/net/kryo/work/updata/fisheries_watson2017")

        } # EO FUN

) # eo 


# --------------------------------------------------------------------------------------------------------------------------------

### 17/04/2020: When the code above is finished, examine results and compute climatologies.
### To do so: 
# - rbind all data
# - sum all catches (reported and unreported separately) within each cell and for each year
# - compute mean catches within each cell based on all years (90-19: 29 years)
setwd(paste(WD,"/","catch_data", sep = ""))
# dir()
# f <- "Catch2010_2014_treated_16_04_10.Rdata"

require("parallel")
res <- mclapply(dir()[grep("20_04_20",dir())], function(f) {
            d <- get(load(f))
            # summary(d)
            d <- d[!is.na(d$Lat),]
            return(d)
        }, mc.cores = 10
) # eo lapply
# Rbind
data <- dplyr::bind_rows(res)
dim(data) # 261 893 626
str(data)
head(data)
rm(res); gc()

total <- data.frame(data %>% group_by(Cell,IYear,TaxonType) %>%
            summarize(x = unique(Long), y = unique(Lat),
                Reported = sum(c(ReportedIND,ReportedNIND)),
                IUU = sum(c(IUUIND,IUUNIND)),
                Discards = sum(c(DiscardsIND,DiscardsNIND))
        )
) # eo ddf
dim(total) # 17042792 -> 132097 cells * 29 years
# length(unique(total$Cell)) # 147632 cells 
str(total)
head(total)

### Use if else loop to change the TaxonType (Pelagic vs. Demersal)
unique(total$TaxonType)
total$Type <- NA
total[which(total$TaxonType %in% unique(total$TaxonType)[c(3,4,7,10,14,16:18,20,22,23,27,29)]),"Type"] <- "Pelagic"
total[which(total$TaxonType %in% unique(total$TaxonType)[c(1:2,5:6,8:9,11:13,15,19,21,24:26,28)]),"Type"] <- "Demersal"
unique(total$Type)
summary(factor(total$Type)) # 61% pelagic fisheries, 38% demersal


### Save so you don't have to re-load the data above
#save(total, file = "table_total_catches_1990-2019v5.Rdata")
total <- get(load("table_total_catches_1990-2019v5.Rdata"))

# Melt to plot the time serie splot with geom_area
ts <- data.frame(total %>% group_by(IYear,Type) %>%
            summarize(Reported = sum(Reported),
                IUU = sum(IUU), Discards = sum(Discards)
        )
) # eo ddf
head(ts) ; dim(ts)
ts$Total <- (ts$Reported)+(ts$IUU)

# Convert
require("lubridate")
ts$Year <- lubridate::ymd(ts$IYear, truncated = 2L)

### Plot time series using geom_area()
ggplot(ts, aes(x = Year, y = Total/1000000, fill = factor(Type))) +
    geom_area(color = "black", size = 0.2, alpha = 0.8, position = 'stack') +
    scale_fill_manual(name = "", values = c("#F21A00","#3B9AB2")) + 
    scale_y_continuous(limits = c(0,140)) + scale_x_date(date_breaks = "3 years", date_labels = "%Y") + 
    labs(x = "Year", y = "Total fisheries catches (Mt)") + theme_minimal()
# Very nice, matches Watson's plot
rm(m); gc()

### Compute annual average from total now
clim <- data.frame(total %>% group_by(Cell) %>%
            summarize(x = unique(x), y = unique(y),
                Reported = mean(Reported, na.rm = T),
                Unreported = mean(IUU, na.rm = T),
                Discarded = mean(Discards, na.rm = T))
        )
) # eo ddf
head(clim) ; dim(clim)
clim$Total <- (clim$Reported)+(clim$Unreported)
summary(clim)

# clim_pel <- clim[clim$Type == "Pelagic",]
# clim_dem <- clim[clim$Type == "Demersal",]

### Mapping time ! 
# Logged otherwise can't see anything
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim) +
   	scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)", option = "viridis") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )  
            
# # Same, for demersals
# map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim) +
#        scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)", option = "viridis") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(60,120,180,-180,-120,-60,0),
#               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
#     scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
#               labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
#     theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#             panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

### Discrete color palette as in Link & Watson
clim$Total2 <- NA
clim[which(clim$Total < 2),"Total2"] <- "0-2"
clim[which(clim$Total >= 2 & clim$Total < 10),"Total2"] <- "2-10"
clim[which(clim$Total >= 10 & clim$Total < 20),"Total2"] <- "10-20"
clim[which(clim$Total >= 20 & clim$Total < 40),"Total2"] <- "20-40"
clim[which(clim$Total >= 40 & clim$Total < 60),"Total2"] <- "40-60"
clim[which(clim$Total >= 60 & clim$Total < 100),"Total2"] <- "60-100"
clim[which(clim$Total >= 100),"Total2"] <- ">100"
# Check nb per category
summary(factor(clim$Total2))           

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(Total2)), data = clim) +
   	scale_fill_manual(name = "Mean annual\ncatches\n(t/km2.yr)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )             

# Awesome, save maps 
ggsave(plot = map1, filename = "map_mean_ann_pelagics_catches_logged_1990-2019.jpg", dpi = 300, width = 7, height = 5) 
ggsave(plot = map2, filename = "map_mean_ann_demersals_catches_logged_1990-2019.jpg", dpi = 300, width = 7, height = 5) 
ggsave(plot = map3, filename = "map_mean_ann_pelagics_catches_discr_1990-2019.jpg", dpi = 300, width = 7, height = 5) 
            
# Save clim
save(clim_pel, file = "clim_pelagic_catches_1990-2019.Rdata")    
save(clim_dem, file = "clim_demersal_catches_1990-2019.Rdata") 



### 23/04/2020: Do the same but only for small surface pelagics. Restrict 'total' to plankton-feeding fisheries 
### But also do it for various size classes of pelagics
unique(total$TaxonType) # restrict to "pelagic <30 cm
small_pel <- total[total$TaxonType == "pelagic <30 cm                                    ",]
med_pel <- total[total$TaxonType == "pelagic 30 - 90 cm                                ",]
large_pel <- total[total$TaxonType == "pelagic >=90 cm                                   ",]
all_pel <- rbind(small_pel, med_pel, large_pel)

dim(small_pel) ; dim(med_pel); dim(large_pel) ; dim(all_pel)
# 2378847, 1995666, 3135447
#head(med_pel)
# save(small_pel, file = "table_total_catches_pelagics30-90cm_1990-2019v5.Rdata")
# save(med_pel, file = "table_total_catches_pelagics30-90cm_1990-2019v5.Rdata")
# save(large_pel, file = "table_total_catches_pelagics>90cm_1990-2019v5.Rdata")

# Clim of small pelagics
clim_small_pel <- data.frame(small_pel %>% group_by(Cell) %>%
            summarize(x = unique(x), y = unique(y),
                Reported = mean(Reported, na.rm = T),
                Unreported = mean(IUU, na.rm = T),
                Discarded = mean(Discards, na.rm = T)
        )
) # eo ddf
clim_small_pel$Total <- (clim_small_pel$Reported)+(clim_small_pel$Unreported)

# Clim of medium pelagics
clim_med_pel <- data.frame(med_pel %>% group_by(Cell) %>%
            summarize(x = unique(x), y = unique(y),
                Reported = mean(Reported, na.rm = T),
                Unreported = mean(IUU, na.rm = T),
                Discarded = mean(Discards, na.rm = T)
        )
) # eo ddf
clim_med_pel$Total <- (clim_med_pel$Reported)+(clim_med_pel$Unreported)

# Clim of large pelagics
clim_large_pel <- data.frame(large_pel %>% group_by(Cell) %>%
            summarize(x = unique(x), y = unique(y),
                Reported = mean(Reported, na.rm = T),
                Unreported = mean(IUU, na.rm = T),
                Discarded = mean(Discards, na.rm = T)
        )
) # eo ddf
clim_large_pel$Total <- (clim_large_pel$Reported)+(clim_large_pel$Unreported)

# Clim of large pelagics
clim_all_pel <- data.frame(all_pel %>% group_by(Cell) %>%
            summarize(x = unique(x), y = unique(y),
                Reported = mean(Reported, na.rm = T),
                Unreported = mean(IUU, na.rm = T),
                Discarded = mean(Discards, na.rm = T)
        )
) # eo ddf
clim_all_pel$Total <- (clim_all_pel$Reported)+(clim_all_pel$Unreported)


# summary(clim_small_pel); summary(clim_med_pel); summary(clim_large_pel)
# summary(log1p(clim_small_pel$Total)); summary(log1p(clim_med_pel$Total)); summary(log1p(clim_large_pel$Total))
# summary(log1p(clim_all_pel$Total))
### maps their logged catches with same scale (0-11.8)
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim_small_pel) +
  	scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)\n(<30cm)", option = "magma", limits = c(0,11.8)) +
   geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
   coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
             labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
	      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
   		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )  
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim) +
  	scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)\n(30-90cm)", option = "magma", limits = c(0,11.8)) +
   geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
   coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
             labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
	      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
   		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )  
#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim_med_pel) +
  	scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)\n(>90cm)", option = "magma", limits = c(0,11.8)) +
   geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
   coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
             labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
   scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
	      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
   		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )  
        
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim_all_pel) +
    scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)\n(all)", option = "magma", limits = c(0,11.8)) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
        labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )  
                               
###
save(clim_small_pel, file = "clim_pelagics<30cm_catches_1990-2019.Rdata") 
save(clim_med_pel, file = "clim_pelagics30-90cm_catches_1990-2019.Rdata") 
save(clim_large_pel, file = "clim_pelagics>90cm_catches_1990-2019.Rdata") 
save(clim_all_pel, file = "clim_all_pelagics_catches_1990-2019.Rdata")
# Save maps
ggsave(plot = map1, filename = "map_ann_total_catches_pelagics<30cm_1990-2019.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map2, filename = "map_ann_total_catches_pelagics30-90cm_1990-2019.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map3, filename = "map_ann_total_catches_pelagics>90cm_1990-2019.jpg", dpi = 300, width = 7, height = 5)
ggsave(plot = map4, filename = "map_ann_total_catches_all_pelagics_1990-2019.jpg", dpi = 300, width = 7, height = 5)



# ----------------------------------------------------------

### 21/04/2020: Complete the annual climatologies of pelagics fisheries catch rates:
# - Degrade resolution from 0.5° to 1° cells.
# - Fill in gaps with grid.expand
# - Make sure it foolows the exact same grid as WOA (1°x1°)
# - Discard land cells using a SST product from WOA 
clim_pel <- get(load("clim_pelagics<30cm_catches_1990-2019.Rdata"))
### First, get the standard cell grid
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors")
ras <- raster("woa13_all_o_monthly.nc") 
ras
#plot(ras)
grid <- as.data.frame(ras, xy = T)
dim(grid) # ok
unique(grid$x)
unique(grid$y)
# Ok, make clim_pel match with these coordinates...try interpolation
clim_pel$x1 <- round(clim_pel$x/0.5)*0.5
clim_pel$y1 <- round(clim_pel$y/0.5)*0.5
clim_pel$id <- factor(paste(clim_pel$x1, clim_pel$y1, sep = "_"))
### Compute mean within these new cells
clim1d <- data.frame(clim_pel %>% group_by(id) %>% summarize(x = unique(x1), y = unique(y1), Total = mean(Total, na.rm =T)) ) # eo ddf
summary(clim1d)
head(grid) # need to re-order clim1d maybe
clim1d <- clim1d[order(clim1d$x, decreasing = F),]
# OK, try to interp this

x <- seq(from = -180, to = 180, by = 1)
y <- seq(from = -90, to = 90, by = 1)
d1 <- expand.grid(x = x, y = y)
d1 <- data.frame(d1)
class(d1)
unique(d1$x)
unique(d1$y)
d1$id <- factor(paste(d1$x, d1$y, sep = "_"))
d1$Total <- as.numeric(NA)
str(d1)

### Make them follow the same ORDER
clim1d <- clim1d[order(clim1d$id),]
d1 <- d1[order(d1$id),]
# Provie id as rownames
rownames(d1) <- d1$id
rownames(clim1d) <- clim1d$id
# Check
head(clim1d)
head(d1)
# Define common grid cells
commons <- intersect(unique(clim1d$id), unique(d1$id)) # length(commons)
# dim(d1[commons,])
d1[commons,"Total"] <- clim1d$Total
summary(d1)
# Replace NA by zeroes
d1$Total[is.na(d1$Total)] <- 0
# Cool, adjust x and y in d1 to match grid
dim(d1) #  65341 
dim(grid) # reference -> 64800
# Need to remove 541 cells
# Remove 0.5 to every positve long and lat
# Add 0.5 to every negative long and lat
d1$x2 <- NA
d1$y2 <- NA
d1[which(d1$x < 0),"x2"] <- (d1[which(d1$x < 0),"x"])+0.5
d1[which(d1$y < 0),"y2"] <- (d1[which(d1$y < 0),"y"])+0.5
d1[which(d1$x >= 0),"x2"] <- (d1[which(d1$x >= 0),"x"])-0.5
d1[which(d1$y >= 0),"y2"] <- (d1[which(d1$y >= 0),"y"])-0.5
summary(d1)
# Recompute an anevarge based on new id
d1$id2 <- factor(paste(d1$x2, d1$y2, sep = "_"))
length(unique(d1$id2)) # 64800 nice 

clim <- data.frame(d1 %>% group_by(id2) %>% summarize(x = unique(x2), y = unique(y2), Total = mean(Total, na.rm = T)) ) # eo ddf
summary(clim)
# Quickmap
ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim) +
   	scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)") + 
    geom_polygon(aes(x = long, y = lat, group = group), data = world, 
        fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme_bw()
# looks ok

### Filter land cells with the 'grid', add an ID to grid
grid$id <- factor(paste(grid$x, grid$y, sep = "_"))
commons <- intersect(unique(grid$id), unique(clim$id2)) # length(commons)
grid <- grid[order(grid$id),]
colnames(grid)[3] <- "O2"
head(grid)
head(clim)
clim$o2 <- grid$O2
summary(clim)
# Convert Total to NA based on O2

clim[is.na(clim$o2),"Total"] <- NA
summary(clim)

ggplot() + geom_raster(aes(x = x, y = y, fill = log1p(Total)), data = clim) +
   	scale_fill_viridis(name = "Mean annual\ncatches\nlog(t/km2.yr)") + 
    geom_polygon(aes(x = long, y = lat, group = group), data = world,
            fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme_bw()

### Nice, save 
clim2save <- clim[,c(1:4)]
head(clim2save)
clim2save$logged <- log1p(clim2save$Total)
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Fisheries_Watson2017/catch_data")
save(clim2save, file = "clim_pelagic<30cm_catches_1990-2019_1d.Rdata")    

### Check with discrete colorscale
clim2save$Total2 <- NA
clim2save[which(clim2save$Total < 2),"Total2"] <- "0-2"
clim2save[which(clim2save$Total >= 2 & clim2save$Total < 10),"Total2"] <- "2-10"
clim2save[which(clim2save$Total >= 10 & clim2save$Total < 20),"Total2"] <- "10-20"
clim2save[which(clim2save$Total >= 20 & clim2save$Total < 40),"Total2"] <- "20-40"
clim2save[which(clim2save$Total >= 40 & clim2save$Total < 60),"Total2"] <- "40-60"
clim2save[which(clim2save$Total >= 60 & clim2save$Total < 100),"Total2"] <- "60-100"
clim2save[which(clim2save$Total >= 100),"Total2"] <- ">100"
# Check nb per category
summary(factor(clim2save$Total2))           

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(Total2)), data = clim2save) +
   	scale_fill_manual(name = "Mean annual\ncatch rates\n(t/km2.yr)",
            breaks = c("0-2","2-10","10-20","20-40","40-60","60-100",">100"),
            values = c("#3288bd","#66c2a5","#abdda4","#fee08b","#fdae61","#f46d43","#d53e4f") ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(60,120,180,-180,-120,-60,0),
              labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )             

ggsave(plot = map3, filename = "map_mean_ann_pelagics_catches_discr_1990-2019_1d.jpg", dpi = 300, width = 7, height = 5) 


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

