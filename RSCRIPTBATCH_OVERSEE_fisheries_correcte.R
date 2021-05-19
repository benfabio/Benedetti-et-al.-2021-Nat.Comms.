
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

cells <- read.csv("codes_cells_Watson17.csv", sep = ";", h = T, dec = ",")
gears <- read.csv("codes_gears_Watson17.csv", sep = ";", h = T)
taxa <- read.csv("codes_taxa_Watson17.csv", sep = ";", h = T)

### For both industrial and non industrial data, concatenate the data based on the index and the code files
setwd("/net/kryo/work/updata/fisheries_watson2017")
files <- dir()[grep("Catch",dir())][9:14]

### Perform a lapply for each file within which you'll perform a mclapply with like 40 cores to retrieve the data from each ID
# f <- files[6] # for testing

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
            catch$bit <- cut(1:nrow(catch), 42, labels = F)
            # b <- 20
            require("parallel")
            res <- mclapply(unique(catch$bit), function(b) {
                
                        message(paste("for bit = ", b, " -------------------------------------------- ", sep = ""))
                        catch3 <- catch[which(catch$bit == b),]
                        # dim(catch3)
                        
                        # Supply Long & Lat from 'cells' using the 'Cell' column for catch2
                        catch3$Long <- NA ; catch3$Lat <- NA ; catch3$area <- NA
                        catch3$GearUsed <- NA ; catch3$TaxonName <- NA ; catch3$CommonName <- NA ; catch3$TaxonType <- NA
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
                        message(paste("providing taxon name and descript", sep = ""))
                        for(t in unique(catch3$Taxonkey) ) {
                            
                            message(paste("# ", which(unique(catch3$Taxonkey) == t), " - ",t," || ", round((which(unique(catch3$Taxonkey) == t)/length(unique(catch3$Taxonkey))),4)*100, "%", sep = ""))   
                            
                            if( isTRUE(t %in% taxa$Taxonkey) ) {
                                catch3[catch3$Taxonkey == t,"TaxonName"] <- as.character(unique(taxa[taxa$Taxonkey == t,"TaxonName"]))
                                catch3[catch3$Taxonkey == t,"CommonName"] <- as.character(unique(taxa[taxa$Taxonkey == t,"CommonName"]))   
                                catch3[catch3$Taxonkey == t,"TaxonType"] <- as.character(unique(taxa[taxa$Taxonkey == t,"Descript"]))   
                            } else {
                                message(paste("Taxon key was not macthed in taxa' index",sep=""))
                            }
                            
                        } # eo for loop - t in Taxonkey
                        # Check
                        # str(catch3) ; summary(catch3)
                        # Return
                        return(catch3)
                
                }, mc.cores = 42

            ) # eo mclapply - b in bits
            # Rbind
            t <- dplyr::bind_rows(res)
            # str(t); head(t); dim(t)
            rm(res) ; gc()
            
            ### Save 
            setwd(paste(WD,"/","catch_data", sep = ""))
            filename <- paste(str_replace(f,".csv",""),"_treated_20_04_20",".Rdata", sep = "") # filename
            message(paste("Saving catch data for ", filename, sep = ""))
            save(t, file = filename)
            rm(t) ; gc()
            setwd("/net/kryo/work/updata/fisheries_watson2017")

        } # EO FUN

) # eo 


#########################################################################################################

