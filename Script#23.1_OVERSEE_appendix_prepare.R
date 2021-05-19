
##### 30/07/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for preparing the Appendices of the Supplementary Material:
#	x Table S1: workflow table summarizing the steps involved in implemenying the zoo dataset
#	x Document S2: list of characters used to clean the fossil/sedimentary records
#	x Figure S3: distribution of sampling effort in space and time (like the SI of Righetti et al., 2019)
#	x Document S4: preliminary tests to evaluate the difference between total-background data and target-background data
#	x Document S5: distribution of sampling effort and richness along global annual mean pCO2 gradients (and SSS for zoo)
#	x Document S6: distribution of predictors rankings and explanatory power of GLMs and RF models
#	x Figure S7: distribution of species TSS scores, per predictors pool
#	- Document S8: robustness of diversiry patterns to the omission of Chl-a from the predictors
#	- Document S9: rarefaction analysis to demonstrate that the dip of zooplankton diversity in the tropics is not driven by sampling biases

# Note: x = done; - = to do


### Last update: 20/03/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("tidyverse")
library("stringr")
library("reshape2")
library("biomod2")
library("viridis")
library("RColorBrewer")
library("scales")
library("maps")
library("ggsci")
library("wesanderson")
library("ggthemes")
library("viridis")

world2 <- map_data(map = "world2")

# --------------------------------------------------------------------------------------------------------------------------------

WD <- getwd()

### 1°) Prepare Figure S3: distribution of sampling effort in space and time (like the SI of Righetti et al., 2019)

### Start from directory : /net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1")
# Load all obs
res <- lapply(dir()[grep("30_04_19",dir())], function(f) {
			# f <- "Foraminifera_30_04_19.Rdata"
			d <- get(load(f))
			# colnames(d)
			return(d[,c("species","genus","family","order","class","phylum","month","year","xbin_1d","ybin_1d","source")])
	} # eo FUN
) # eo lapply
ddf <- dplyr::bind_rows(res) 
dim(ddf);str(ddf)
rm(res);gc()
setwd(WD)
head(ddf)

# Add cell id
ddf$id <- factor(paste(ddf$xbin, ddf$ybin, sep = "_"))
### Compute and sampling effort 
ddf$x2 <- ddf$xbin_1d 
ddf[ddf$xbin_1d < 0 ,"x2"] <- (ddf[ddf$xbin_1d < 0 ,"xbin_1d"]) + 360
# head(ddf)
require("dplyr")
effort <- data.frame(ddf %>% group_by(id) %>% summarize(x = unique(x2), y = unique(ybin_1d), n = n(), rich = length(unique(species))) ) # eo ddf
dim(effort)
summary(effort)
effort$logn <- log(effort$n)

### Make bins for nicer color palette
effort$logn_bin <- factor(cut_interval(effort$logn,9))
levels(effort$logn_bin)
levels <- str_replace_all(levels(effort$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort$logn_bin) <- levels
levels(effort$logn_bin)[levels(effort$logn_bin) == "0-1.06"] <- "0-1.0"
levels(effort$logn_bin)[levels(effort$logn_bin) == "1.06-2.12"] <- "1.0-2.0"
levels(effort$logn_bin)[levels(effort$logn_bin) == "2.12-3.19"] <- "2.0-3.2"
levels(effort$logn_bin)[levels(effort$logn_bin) == "3.19-4.25"] <- "3.2-4.2"
levels(effort$logn_bin)[levels(effort$logn_bin) == "4.25-5.31"] <- "4.2-5.3"
levels(effort$logn_bin)[levels(effort$logn_bin) == "5.31-6.37"] <- "5.3-6.4"
levels(effort$logn_bin)[levels(effort$logn_bin) == "6.37-7.44"] <- "6.4-7.4"
levels(effort$logn_bin)[levels(effort$logn_bin) == "7.44-8.5"] <- "7.4-8.5"
levels(effort$logn_bin)[levels(effort$logn_bin) == "8.5-9.56"] <- "8.5-9.5"

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(logn_bin)), data = effort) + 
			scale_fill_viridis(name = "Sampling effort\n(log)", discrete = T, option = "inferno") +
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
           		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

setwd("/net/hydro/work/fabioben/OVERSEE/data/")
ggsave(plot = map, filename = "map_effort_logn_appendix3_zoo.jpg", dpi = 300, width = 7, height = 4)

### Nice, now plot density distribution per month and 1° lat bin
effort2 <- data.frame(ddf %>% group_by(month,ybin_1d) %>% summarize(n = n(), rich = length(unique(species))) )
summary(effort2)
effort2$logn <- log(effort2$n)
effort2$logn_bin <- factor(cut_interval(effort2$logn,9))
levels(effort$logn_bin)
levels <- str_replace_all(levels(effort2$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort2$logn_bin) <- levels
# Rename better
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "0-0.967"] <- "0-1.0"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "0.967-1.93"] <- "1.0-2.0"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "1.93-2.9"] <- "2.0-3.0"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "2.9-3.87"] <- "3.0-3.9"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "3.87-4.83"] <- "3.9-4.8"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "4.83-5.8"] <- "4.8-5.8"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "5.8-6.77"] <- "5.8-6.8"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "6.77-7.73"] <- "6.8-7.7"
levels(effort2$logn_bin)[levels(effort2$logn_bin) == "7.73-8.7"] <- "7.7-8.7"

plot <- ggplot() + geom_raster(aes(x = factor(ybin_1d), y = factor(month), fill = factor(logn_bin)), data = effort2) + 
			scale_fill_viridis(name = "Sampling effort\n(log)", discrete = T, option = "inferno") +
			ylab("Month") + xlab("Latitude (°)") + theme_light() + 
			theme(axis.text.x = element_text(angle=90, hjust=1, size = 5)) 
#
ggsave(plot = plot, filename = "plot_effort_months_lat_logn_appendix3.jpg", dpi = 300, width = 13, height = 4)

### Gut, and now just latitudinal density distribution
# Lets make 2° or 5° bins first
ddf$ybin_5d <- ceiling(ddf$ybin_1d/5)*5
#summary(ddf$ybin_5d)
unique(ddf$ybin_5d)
effort3 <- data.frame(ddf %>% group_by(ybin_1d) %>% summarize(n = n(), rich = length(unique(species))) )
summary(effort3)

# Make plot
plot <- ggplot(effort3) + geom_col(aes(x = factor(ybin_1d), y = n), fill = inferno(9)[3], color = "black") + 
		ylab("Sampling effort") + xlab("Latitude (°)") + theme_light() +
		theme(axis.text.x = element_text(angle=90)) #+ theme(axis.text.y = element_text(angle=90)) + 
		#coord_flip()
#
ggsave(plot = plot, filename = "plot_effort_lat_n_appendix3_v3.jpg", dpi = 300, width = 13, height = 4)


### Nice, and finally: same as effort2 but per...years!
effort4 <- data.frame(ddf %>% group_by(year,ybin_5d) %>% summarize(n = n(), rich = length(unique(species))) )
summary(effort4)
effort4$logn <- log(effort4$n)
effort4$logn_bin <- factor(cut_interval(effort4$logn,9))
levels(effort4$logn_bin)
levels <- str_replace_all(levels(effort4$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort4$logn_bin) <- levels
# Rename better
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "0-0.973"] <- "0-1.0"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "0.973-1.95"] <- "1.0-2.0"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "1.95-2.92"] <- "2.0-3.0"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "2.92-3.89"] <- "3.0-4.0"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "3.89-4.87"] <- "4.0-4.9"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "4.87-5.84"] <- "4.9-5.8"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "5.84-6.81"] <- "5.8-6.8"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "6.81-7.79"] <- "6.8-7.8"
levels(effort4$logn_bin)[levels(effort4$logn_bin) == "7.79-8.76"] <- "7.8-8.7"

plot <- ggplot() + geom_raster(aes(x = factor(ybin_5d), y = factor(year), fill = factor(logn_bin)), data = effort4) + 
			scale_fill_viridis(name = "Sampling effort\n(log)", discrete = T, option = "inferno") +
			ylab("Year") + xlab("Latitude (5°)") + theme_light() + 
			theme(axis.text.x = element_text(angle=90, hjust=1, size = 5)) 
#
ggsave(plot = plot, filename = "plot_effort_years_lat_logn_appendix3.jpg", dpi = 300, width = 7, height = 14)


### 20/03/2020: Same as above but for phytoplankton
setwd("/net/hydro/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
files <- dir()[grep("data_total_",dir())]
# f <- files[11]
res <- lapply(files, function(f) {
			d <- read.table(f, h = T, sep = ";")
			return(d[,c("species","genus","family","order","class","phylum","month","av.year","cell_id","x","y")])
	} # eo FUN
) # eo lapply
ddf <- dplyr::bind_rows(res) 
dim(ddf);str(ddf)
rm(res);gc()
head(ddf)
unique(ddf$y)

### Compute and sampling effort 
require("maps")
world2 <- map_data(map = "world2")
ddf$x2 <- ddf$x 
ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0,"x"]) + 360
# head(ddf)
require("dplyr")
effort <- data.frame(ddf %>% group_by(cell_id) %>% summarize(x = unique(x2), y = unique(y), n = n(), rich = length(unique(species))) ) # eo ddf
dim(effort)
summary(effort)
effort$logn <- log(effort$n)

### Make bins for nicer color palette
effort$logn_bin <- factor(cut_interval(effort$logn,9))
levels(effort$logn_bin)
levels <- str_replace_all(levels(effort$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort$logn_bin) <- levels
levels(effort$logn_bin)[levels(effort$logn_bin) == "4.36-4.78"] <- "4.3-4.8"
levels(effort$logn_bin)[levels(effort$logn_bin) == "4.78-5.2"] <- "4.8-5.2"
levels(effort$logn_bin)[levels(effort$logn_bin) == "5.2-5.62"] <- "5.2-5.6"
levels(effort$logn_bin)[levels(effort$logn_bin) == "5.62-6.03"] <- "5.6-6.1"
levels(effort$logn_bin)[levels(effort$logn_bin) == "6.03-6.45"] <- "6.1-6.5"
levels(effort$logn_bin)[levels(effort$logn_bin) == "6.45-6.87"] <- "6.5-6.9"
levels(effort$logn_bin)[levels(effort$logn_bin) == "6.87-7.29"] <- "6.9-7.3"
levels(effort$logn_bin)[levels(effort$logn_bin) == "7.29-7.71"] <- "7.3-7.7"
levels(effort$logn_bin)[levels(effort$logn_bin) == "7.71-8.13"] <- "7.7-8.2"

library("ggsci")
library("wesanderson")
library("ggthemes")
library("viridis")

map <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(logn_bin)), data = effort) + 
			scale_fill_viridis(name = "Sampling effort\n(log)", discrete = T, option = "inferno") +
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
           		labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

setwd("/net/hydro/work/fabioben/OVERSEE/data/")
ggsave(plot = map, filename = "map_effort_logn_appendix3_phyto.jpg", dpi = 300, width = 7, height = 4)




# -------------------------------------------------------------


### 2°) Difference between total-background data and target-background data ? 
### Some material is already on your computer (TSS distributions per groups), but you do need to replot the diversity patterns (total and Chaetognatha)

require("maps")
world2 <- map_data(map = "world2")

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/group_bckgrnd/niche.modelling")
group.wd <- getwd()
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/niche.modelling")
tot.wd <- getwd()

# Get species list
setwd(tot.wd)
species <- dir()
# Get spp scores to assess whether a model has failed (NA values in scores table)
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/eval_scores/")
files <- dir()
spp <- str_replace_all(files,"TSS_","")
spp <- str_replace_all(spp,".Rdata","")
scores <- mclapply(spp, function(sp) {
 				message(paste("Loading scores for ",sp, sep = ""))
 				s <- get(load(paste("TSS_",sp,".Rdata", sep = "")))
 				s$species <- sp
 				return(s)
 			}, mc.cores = 15
) # eo lapply
scores_tbl <- dplyr::bind_rows(scores)
rm(scores) ; gc()
# Add an ENM column
summary(scores_tbl)
scores <- data.frame(scores_tbl %>% group_by(species) %>% summarize(TSS = mean(TSS)))
nrow(scores[scores$TSS < 0.3,]) # 19 spp with group background; 18 Chaetognatha and 1 Pteropod (Limacina retroversa)
# ZERO spp with mean TSS < 0.3 using total background
# scores[scores$TSS < 0.3,] # 

### Choose a season/month
s <- "winter"

# Retreive coordinates 
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
if(s == "spring") {
	m <- "apr"	
} else if(s == "fall") {
	m <- "oct"
} else if(s == "summer") {
	m <- "aug"
} else if (s == "winter") {
	m <- "jan"
} # eo else if loop
# Get corresponding clims
coords <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
coords$id <- factor(paste(coords$x, coords$y, sep = "_"))

# Get spp probabilities
require("parallel")
#sp <- "Oncaea.media"
#sp <- "Thetys.vagina"
probas <- mclapply(species, function(sp) {
	
			# Got to species dir
			setwd( paste(tot.wd,"/",sp,"/", sep = "") )
			message(paste("Loading projections for ", sp, "  ================================", sep = ""))
		
			# Need to modify sp when there are 2 names and add brackets
			if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3 ) {
				# Then add brackets around the second piece
				sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
						strsplit(sp, ".", fixed = TRUE)[[1]][2],").", 
						strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )		 
			}
		
			# If the 4 seasonal projections are done
			if( sum(grepl("proj_projection_", dir())) == 4 ) {
			
				### Load projections for each SDM
				setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",s, sep = ""),"/", sep = "") )
				d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",s,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
				# GAM
				resModelGam <- d[,"GAM",,]
				resModelGam <- apply(resModelGam, 1, mean, na.rm = F) 
				resModelGam <- (resModelGam/1000) 
				# GLM
				resModelGlm <- d[,"GLM",,]
				resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F) 
				resModelGlm <- (resModelGlm/1000) 
				# RF
				resModelRF <- d[,"RF",,]
				resModelRF <- apply(resModelRF, 1, mean, na.rm = F) 
				resModelRF <- (resModelRF/1000) 
				# ANN
				resModelANN <- d[,"ANN",,]
				resModelANN <- apply(resModelANN, 1, mean, na.rm = F) 
				resModelANN <- (resModelANN/1000) 

				# Return
				return( data.frame(cell_id = paste(coords$x, coords$y, sep = "_"), x = coords$x, y = coords$y, species = sp, 
						GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
			
			} else {
			
				message(paste("Skipping because no projection (yet)", sp, "  ================================", sep = ""))
			
			} # eo if else loop

		}, mc.cores = 20
) # eo lapply
# Rbind
tbl <- dplyr::bind_rows(probas)
rm(probas); gc()

# Compute total species richness by computing species average HSI (average across SDMs) and then the sum of these average
tbl$mean_HSI <- rowMeans( as.matrix(tbl[,c(5:length(tbl))]) )

# rm cells 
coords <- coords[-which(coords$SSS < 20),] 
coords <- coords[-which(coords$Bathy > -175),] 
cells2keep <- unique(coords$id)
tbl <- tbl[which(tbl$cell_id %in% cells2keep),]
dim(tbl)


### Calculate sum of HSI
require("dplyr")
rich <- data.frame(tbl %>%
			group_by(cell_id) %>%
			summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI) )
) # eo ddf
summary(rich)
# And flip x coordinates to match coastline
rich$x2 <- rich$x 
rich[rich$x < 0 ,"x2"] <- (rich[rich$x < 0 ,"x"]) + 360

# Go to monthly plot.dir
setwd(WD)
require("viridis")
require("ggplot2")
map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = rich), data = rich) + 
			scale_fill_viridis(name = "Species\nrichness") + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
	
ggsave(plot = map, filename = paste("map_group_rich_",s,".jpg", sep = ""), dpi = 300, height = 3, width = 6)

### To plot diversity patterns per groups 
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/occurence_data_groups/v9/v9v8v5.1v3.1/")
example <- get(load(dir()[1]))
names <- colnames(example)[c(1:8,10:13)]
rm(example)
files <- dir()[grep("_30_04_19",dir())]
# Get all obs data
matched <- mclapply(files, function(f) {
 				message(paste(f, sep = ""))
				data <- get(load(f))
 				return(data[,names])
 			}, mc.cores = 15
 ) # eo lapply
match <- dplyr::bind_rows(matched)
rm(matched) ; gc()
all_spp_names <- unique(match$species)
# length(all_spp_names)
counts <- data.frame(match %>%
 		group_by(species) %>%
 		summarise(n = n(), fam = unique(family)[1], ord = unique(order)[1], class = unique(class)[1], phylum = unique(phylum)[1] )
)
# Define groupings for per group plots (maps and eval scores)
counts[counts$class == "Maxillopoda","class"] <- "Hexanauplia"
counts$groupings <- NA
counts[counts$class == "Hexanauplia","groupings"] <- "Copepoda"
counts[counts$phylum == "Mollusca","groupings"] <- "Pteropoda"
counts[counts$phylum == "Chaetognatha","groupings"] <- "Chaetognatha"
counts[counts$phylum == "Cnidaria","groupings"] <- "Jellyfish"
counts[counts$phylum == "Ctenophora","groupings"] <- "Jellyfish"
counts[counts$phylum == "Foraminifera","groupings"] <- "Foraminifera"
counts[counts$phylum == "Chordata","groupings"] <- "Chordata"
counts[counts$phylum == "Annelida","groupings"] <- "Annelida"
counts[counts$class == "Ostracoda","groupings"] <- "Other_arthropoda"
counts[counts$class == "Branchiopoda","groupings"] <- "Other_arthropoda"
counts[counts$class == "Malacostraca","groupings"] <- "Malacostraca"


# Provide groupings to tbl
tbl$phylum <- NA
tbl$class <- NA
tbl$ord <- NA
tbl$groupings <- NA
tbl$species <- gsub("\\.", "_", tbl$species)
# Provide class and groupings to tbl # - takes some time since 'tbl' has many many rows
for(sp in unique(tbl$species) ) {
		# Get corresponding phyl name from counts and provide to tbl
		tbl[tbl$species == sp,"phylum"] <- unique(as.character(counts[counts$species == sp,"phylum"]))
		tbl[tbl$species == sp,"class"] <- unique(as.character(counts[counts$species == sp,"class"]))
		tbl[tbl$species == sp,"ord"] <- unique(as.character(counts[counts$species == sp,"ord"]))
		tbl[tbl$species == sp,"groupings"] <- unique(as.character(counts[counts$species == sp,"groupings"]))
} # eo for loop - sp in unique(tbl$species)

# And plot richness patterns
# g <- "Jellyfish"
for(g in unique(tbl$groupings) ) {
		message(paste("Mapping global richness of ", g, sep = ""))
		#nsp <- length(unique(tbl[tbl$groupings == g,"species"]))
		tbl2 <- tbl[tbl$groupings == g,]
		require("dplyr")
		rich2 <- data.frame(tbl2 %>% 
				group_by(cell_id) %>%
				summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI) )
		) # eo ddf
		rich2$x2 <- rich2$x 
		rich2[rich2$x < 0 ,"x2"] <- (rich2[rich2$x < 0 ,"x"]) + 360
		setwd(WD)
		map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = rich), data = rich2) + 
					scale_fill_viridis(name = "Species\nrichness") + 
					geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
					coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
		          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
					scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
			      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
					theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
						panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

		ggsave(plot = map, filename = paste("map_rich_total_",m,"_",g,".jpg", sep = ""), dpi = 300, height = 3, width = 6)	
} # eo for loop - g in unique(tbl$groupings)






# -------------------------------------------------------------

### 3°) Map of global annual pCO2 and SSS for appendix 5 (rest of material already on computer? )
require("maps")
world2 <- map_data(map = "world2")
# Get all 12 monthly clims
files <- dir()[grep("_21_02_19",dir())]
clims <- lapply(files, function(f) {
				d <- read.table(f, h = T, sep = ";")
				d$id <- factor(paste(d$x,d$y,sep = "_"))
				return(d[,c("id","x","y","SSS","pCO2")])
		} # eo FUN
) # lapply
clims <- dplyr::bind_rows(clims)
summary(clims)
# Compute annuals climatologies of SS, and pCO2
annual <- data.frame(clims %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), SSS = mean(SSS,na.rm=T), pCO2 = mean(pCO2,na.rm=T)) ) 
summary(annual)
# And flip x coordinates to match coastline
annual$x2 <- annual$x 
annual[annual$x < 0 ,"x2"] <- (annual[annual$x < 0 ,"x"]) + 360
#quartz()
map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = pCO2), data = annual) + 
			#scale_fill_viridis(name = "Annual pCO2\n(µatm)") + 
			scale_fill_distiller(name = "Annual pCO2", palette = "RdYlBu", direction = -1) + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#
#quartz()
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = SSS), data = annual) + 
			#scale_fill_viridis(name = "Annual SSS", limits = c(30,41), option = "A") + 
			scale_fill_distiller(name = "Annual SSS", limits = c(29,41), palette = "RdYlBu", direction = -1) + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#	
ggsave(plot = map1, filename = paste("map_annual_pCO2.jpg", sep = ""), dpi = 300, height = 3, width = 6)
ggsave(plot = map2, filename = paste("map_annual_SSS.jpg", sep = ""), dpi = 300, height = 3, width = 6)				


### 02/08/2019: Re-map by binning the palette
# For SSS ONLY BETWEEN 25 and 41
annual2 <- annual[annual$SSS >= 30,c("id","x","y","SSS","x2")]
annual2$SSS_bin <- factor(cut_interval(annual2$SSS,9))
levels(annual2$SSS_bin)
levels <- str_replace_all(levels(annual2$SSS_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(annual2$SSS_bin) <- levels
map1 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(SSS_bin)), data = na.omit(annual2) ) + 
			scale_fill_brewer(name = "Annual SSS", palette = "RdYlBu", direction = -1) + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map1, filename = paste("map_annual_SSS_bins.jpg", sep = ""), dpi = 300, height = 3, width = 6)


### And for pCO2
annual3 <- annual[,c("id","x","y","pCO2","x2")]
annual3$pCO2_bin <- factor(cut_interval(annual3$pCO2,9))
levels(annual3$pCO2_bin)
levels <- str_replace_all(levels(annual3$pCO2_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(annual3$pCO2_bin) <- levels
map2 <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(pCO2_bin) ), data = na.omit(annual3) ) + 
			scale_fill_brewer(name = "Annual pCO2", palette = "RdYlBu", direction = -1) + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
#	
ggsave(plot = map2, filename = paste("map_annual_pCO2_bin.jpg", sep = ""), dpi = 300, height = 3, width = 6)				



# -------------------------------------------------------------


### 4°) Zooplankton species richness patterns and SDMs' TSS distrbution for the 3 predictors set that highlihhted the impact of SSS
# - with SSS: test 10 - SST, SSS, dSST, dO2, logChl, logNO3
# - with Lon instead of SSS: test 11 - SST, Lon, dSST, dO2, logChl, logNO3
# - without SSS or Lon: test 12 - SST, dSST, dO2, logChl, logNO3

# For the maps, as always
require("maps")
world2 <- map_data(map = "world2")

### You already have the TSS distribution plots for all tests and all SDMs/groups, just need to re-make some maps like above: 
# Set the directories
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/niche.modelling10_05_02")
wd.t10 <- getwd()
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/niche.modelling11_15_03")
wd.t11 <- getwd()
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_v9data_for_tests/species_data_v9v3.1/total_bckgrnd/niche.modelling12_18_03")
wd.t12 <- getwd()

### Vector of tests
tests <- c("t10","t11","t12")

### Vector of seasons/months
months <- c("apr","jul","oct","jan")

# For testing the double for loop below: 
t <- "t10"
m <- "apr"

### For plotting diversity patterns in a for loop
for(t in tests) {
	
		message(paste("Making maps of richness for test ", t, "  ================================", sep = ""))
		if(t == "t10") {
			tot.wd <- wd.t10
		} else if (t == "t11") {
			tot.wd <- wd.t11
		} else if(t == "t12") {
			tot.wd <- wd.t12
		} # eo else if loop
	
		for(m in months) {
		
			message(paste("Making maps of richness for ", m, "  ================================", sep = ""))
			setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
			coords <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
			coords$id <- factor(paste(coords$x, coords$y, sep = "_"))
			# rm cells outside open ocean realm
			coords <- coords[-which(coords$SSS < 20),] 
			coords <- coords[-which(coords$Bathy > -175),] 
		
			setwd(tot.wd)
			species <- dir()
			# sp <- "Agetus.flaccus"
			require("parallel")
			probas <- mclapply(species, function(sp) {
	
						# Got to species dir
						setwd( paste(tot.wd,"/",sp,"/", sep = "") )
						message(paste("Loading projections for ", sp, "  ================================", sep = ""))
		
						# Need to modify sp when there are 2 names and add brackets
						if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3 ) {
							# Then add brackets around the second piece
							sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
									strsplit(sp, ".", fixed = TRUE)[[1]][2],").", 
									strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )		 
						}
		
						# If the 4 seasonal projections are done
						if( sum(grepl("proj_projection_", dir())) == 4 ) {
			
							### Load projections for each SDM
							setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
							d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
							# GAM
							resModelGam <- d[,"GAM",,]
							resModelGam <- apply(resModelGam, 1, mean, na.rm = F) 
							resModelGam <- (resModelGam/1000) 
							# GLM
							resModelGlm <- d[,"GLM",,]
							resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F) 
							resModelGlm <- (resModelGlm/1000) 
							# RF
							resModelRF <- d[,"RF",,]
							resModelRF <- apply(resModelRF, 1, mean, na.rm = F) 
							resModelRF <- (resModelRF/1000) 
							# ANN
							resModelANN <- d[,"ANN",,]
							resModelANN <- apply(resModelANN, 1, mean, na.rm = F) 
							resModelANN <- (resModelANN/1000) 

							# Return
							return( data.frame(cell_id = coords$id, x = coords$x, y = coords$y, species = sp, 
									GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
			
						} else {
			
							message(paste("Skipping because no projection (yet)", sp, "  ================================", sep = ""))
			
						} # eo if else loop

					}, mc.cores = 20
			) # eo lapply
			
			# Rbind
			tbl <- dplyr::bind_rows(probas)
			rm(probas); gc()
			# Compute total species richness by computing species average HSI (average across SDMs) and then the sum of these average
			tbl$mean_HSI <- rowMeans( as.matrix(tbl[,c(5:length(tbl))]) )

			#cells2keep <- unique(coords$id)
			#tbl <- tbl[which(tbl$cell_id %in% cells2keep),]

			### Calculate sum of HSI
			require("dplyr")
			rich <- data.frame(tbl %>%
						group_by(cell_id) %>%
						summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI) )
			) # eo ddf
			# And flip x coordinates to match coastline
			rich$x2 <- rich$x 
			rich[rich$x < 0 ,"x2"] <- (rich[rich$x < 0 ,"x"]) + 360
			rm(tbl);gc()

			# Save map as .jpg in WD
			setwd(WD)
			require("viridis")
			require("ggplot2")
			map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = rich), data = rich) + 
						scale_fill_viridis(name = "Species\nrichness") + 
						geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
						coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
			          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
						scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
				      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
						theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
							panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
	
			ggsave(plot = map, filename = paste("map_total_rich_",m,"_",t,".jpg", sep = ""), dpi = 300, height = 3, width = 6)
		
	
		} # EO second for loop - m in months
	
	
} # EO first for loop - t in test


### 14/08/19: Ok, and then compute the difference between t12 and the two others to highlight the bias towards higher SSS 
m <- "apr"

divs <- lapply(tests, function(t) {
	
			message(paste("Extracting global zooplankton richness for test ", t, sep = ""))
			if(t == "t10") {
				tot.wd <- wd.t10
			} else if (t == "t11") {
				tot.wd <- wd.t11
			} else if(t == "t12") {
				tot.wd <- wd.t12
			} # eo else if loop

			message(paste("Making maps of richness for ", m, "  ================================", sep = ""))
			setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
			coords <- read.table(paste("glob_stack_month_",m,"_21_02_19.txt", sep = ""), h = T, sep = ";")
			coords$id <- factor(paste(coords$x, coords$y, sep = "_"))
			# rm cells outside open ocean realm
			coords <- coords[-which(coords$SSS < 20),] 
			coords <- coords[-which(coords$Bathy > -175),] 
		
			setwd(tot.wd)
			species <- dir()
			# sp <- "Agetus.flaccus"
			require("parallel")
			probas <- mclapply(species, function(sp) {
	
						# Got to species dir
						setwd( paste(tot.wd,"/",sp,"/", sep = "") )
						message(paste("Loading projections for ", sp, "  ================================", sep = ""))
		
						# Need to modify sp when there are 2 names and add brackets
						if( length(strsplit(sp, ".", fixed = TRUE)[[1]]) == 3 ) {
							# Then add brackets around the second piece
							sp <- paste(strsplit(sp, ".", fixed = TRUE)[[1]][1],".(",
									strsplit(sp, ".", fixed = TRUE)[[1]][2],").", 
									strsplit(sp, ".", fixed = TRUE)[[1]][3], sep = "" )		 
						}
		
						# If the 4 seasonal projections are done
						if( sum(grepl("proj_projection_", dir())) == 4 ) {
			
							### Load projections for each SDM
							setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
							d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
							# GAM
							resModelGam <- d[,"GAM",,]
							resModelGam <- apply(resModelGam, 1, mean, na.rm = F) 
							resModelGam <- (resModelGam/1000) 
							# GLM
							resModelGlm <- d[,"GLM",,]
							resModelGlm <- apply(resModelGlm, 1, mean, na.rm = F) 
							resModelGlm <- (resModelGlm/1000) 
							# RF
							resModelRF <- d[,"RF",,]
							resModelRF <- apply(resModelRF, 1, mean, na.rm = F) 
							resModelRF <- (resModelRF/1000) 
							# ANN
							resModelANN <- d[,"ANN",,]
							resModelANN <- apply(resModelANN, 1, mean, na.rm = F) 
							resModelANN <- (resModelANN/1000) 

							# Return
							return( data.frame(cell_id = coords$id, x = coords$x, y = coords$y, species = sp, 
									GLM = resModelGlm, GAM = resModelGam, RF = resModelRF, ANN = resModelANN ) )
			
						} else {
			
							message(paste("Skipping because no projection (yet)", sp, "  ================================", sep = ""))
			
						} # eo if else loop

					}, mc.cores = 25
					
			) # eo lapply
			
			# Rbind
			tbl <- dplyr::bind_rows(probas)
			rm(probas); gc()
			# Compute total species richness by computing species average HSI (average across SDMs) and then the sum of these average
			tbl$mean_HSI <- rowMeans( as.matrix(tbl[,c(5:length(tbl))]) )
			# Calculate sum of HSI
			require("dplyr")
			rich <- data.frame(tbl %>%
						group_by(cell_id) %>%
						summarise(x = unique(x), y = unique(y), rich = sum(mean_HSI) )
			) # eo ddf
			# Add test info
			rich$test <- t
			# Return
			return(rich)
	
		} # EO first for loop - t in test
		
) # eo lapply		
# Examine 'divs' and rbind
str(divs)
div <- dplyr::bind_cols(divs)
dim(div) # 3 times 41553 cells
head(div) 
rm(divs);gc()
# Compute difference in diversity: t10-t12 and t11-t12 (difference to no SSS nor Long)
colnames(div)
div <- div[,c("cell_id","x","y","rich","rich1","rich2")] # rich = test10, rich1 = test11 & rich2 = test12
div$diff1 <- (div$rich) - (div$rich2) # t10-t12
div$diff2 <- (div$rich1) - (div$rich2) # t11-t12
summary(div)

# And flip x coordinates to match coastline and map
div$x2 <- div$x 
div[div$x < 0 ,"x2"] <- (div[div$x < 0 ,"x"]) + 360

setwd(WD)
require("viridis")
require("ggplot2")
map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = diff1), data = div) + 
			scale_fill_gradient2(name = "Species richness\ndifference", low = "#3288bd", high = "#d53e4f", mid = "white",
			limits = c(-50,20)) + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_total_rich_diff_",m,"_","t10-t12",".jpg", sep = ""), dpi = 300, height = 3, width = 6)

map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = diff2), data = div) + 
			scale_fill_gradient2(name = "Species richness\ndifference",low = "#3288bd", high = "#d53e4f", mid = "white",
			limits = c(-50,20)) + 
			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
			coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
          	 	labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
			scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
	      		labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
			theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
				panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

ggsave(plot = map, filename = paste("map_total_rich_diff_",m,"_","t11-t12",".jpg", sep = ""), dpi = 300, height = 3, width = 6)



# -------------------------------------------------------------


### 5°) For appendix 6: distribution of predictors ranks for the two tests (GLMs and RF) for total zooplankton and its sub-groups
### Also find a way to show changes in mean R2 + text to explain how you finally chose your predictors

### A) Phytoplankton: distribution of R2 and ranks per variable pool and then by groups
scores <- read.table("poolz_skillz_table_phyto_tot_26_03_19.txt", h = T, sep = ";")
#dim(scores)
head(scores)
unique(scores$grouping)
unique(scores$pool)
scores <- scores[scores$grouping %in% c("haptophyta","bacillariophyceae","dinoflagellata"),]
# Add a capital at the beginning of the group level
levels(scores$grouping)[levels(scores$grouping) == "haptophyta"] <- "Haptophyta"
levels(scores$grouping)[levels(scores$grouping) == "bacillariophyceae"] <- "Bacillariophyceae"
levels(scores$grouping)[levels(scores$grouping) == "dinoflagellata"] <- "Dinoflagellata"

# For each pool, display the set of variables
for(p in unique(scores$pool)) {
		vars <- unique(scores[scores$pool == p,"var"])
		message(paste("Vars in pool ",p," are : ", sep = ""))
		message(paste(vars, sep = ""))
		message(paste("---------------------------", sep = ""))
} # eo for loop

# To extend color palette
colourCount = length(unique(scores$pool))
getPalette = colorRampPalette(brewer.pal(9,"Set3"))

### Distribution of models R2 per pools of variables
#quartz()
plot1 <- ggplot(scores[scores$test == "test1",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (GLM)") + scale_y_continuous(limits = c(-0.1,1))
#
#quartz()
plot2 <- ggplot(scores[scores$test == "test2",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (RF)") + scale_y_continuous(limits = c(-0.1,1))

#quartz()
plot3 <- ggplot(scores[scores$test == "test1",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (GLM)") + scale_y_continuous(limits = c(-0.1,1)) + 
			facet_grid(factor(grouping) ~ .)
#quartz()
plot4 <- ggplot(scores[scores$test == "test2",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (RF)") + scale_y_continuous(limits = c(-0.1,1)) + 
			facet_grid(factor(grouping) ~ .)
# Save plots 
ggsave(plot = plot1, filename = "boxplot_R2_pools_test1_phyto.jpg", dpi = 300, height = 4, width = 4)
ggsave(plot = plot2, filename = "boxplot_R2_pools_test2_phyto.jpg", dpi = 300, height = 4, width = 4)
ggsave(plot = plot3, filename = "boxplot_R2_pools_test1_phyto_groups.jpg", dpi = 300, height = 8, width = 4)
ggsave(plot = plot4, filename = "boxplot_R2_pools_test2_phyto_groups.jpg", dpi = 300, height = 8, width = 4)


### Distribution of variables' normalized ranks
#quartz()
plot1 <- ggplot(scores[scores$test == "test1",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = T) + 
			scale_fill_brewer(name = "", palette = "Set3") + theme_bw() +
			xlab("Variables") + ylab("Variable ranking (GLM)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) 
#quartz()
plot2 <- ggplot(scores[scores$test == "test2",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = T) + 
			scale_fill_brewer(name = "", palette = "Set3") + theme_bw() +
			xlab("Variables") + ylab("Variable ranking (RF)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) 
# Save
ggsave(plot = plot1, filename = "boxplot_ranks_test1_phyto.jpg", dpi = 300, height = 4, width = 6)
ggsave(plot = plot2, filename = "boxplot_ranks_test2_phyto.jpg", dpi = 300, height = 4, width = 6)

#quartz()
plot1 <- ggplot(scores[scores$test == "test1",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = F) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("Variable") + ylab("Variable ranking (GLM)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) + facet_grid(factor(grouping) ~ factor(pool) )
#
#quartz()
plot2 <- ggplot(scores[scores$test == "test2",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = F) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("Variable") + ylab("Variable ranking (RF)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) + facet_grid(factor(grouping) ~ factor(pool) )
# Save
ggsave(plot = plot1, filename = "boxplot_ranks_test1_phyto_groups.jpg", dpi = 300, height = 8, width = 25)
ggsave(plot = plot2, filename = "boxplot_ranks_test2_phyto_groups.jpg", dpi = 300, height = 8, width = 25)


####################################  And now for Zooplankton
### B) Zooplankton plot distribution of R2 and ranks per variable pool and then by groups
scores <- read.table("poolz_skillz_table_zoo_tot_26_03_19.txt", h = T, sep = ";")
unique(scores$grouping)
unique(scores$pool)
# Filter groups of interest
scores <- scores[scores$grouping %in% c("Jellyfish","Copepoda","Chaetognatha","Malacostraca","Pteropoda","Foraminifera","Chordata"),]

# For each pool, display the set of variables
for(p in unique(scores$pool)) {
		vars <- unique(scores[scores$pool == p,"var"])
		message(paste("Vars in pool ",p," are : ", sep = ""))
		message(paste(vars, sep = ""))
		message(paste("---------------------------", sep = ""))
} # eo for loop

colourCount = length(unique(scores$pool))
getPalette = colorRampPalette(brewer.pal(9,"Set3"))

### Distribution of models R2 per pools of variables
#quartz()
plot1 <- ggplot(scores[scores$test == "test1",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (GLM)") + scale_y_continuous(limits = c(-0.1,1))
#
#quartz()
plot2 <- ggplot(scores[scores$test == "test2",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (RF)") + scale_y_continuous(limits = c(-0.1,1))

#quartz()
plot3 <- ggplot(scores[scores$test == "test1",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (GLM)") + scale_y_continuous(limits = c(-0.1,1)) + 
			facet_grid(factor(grouping) ~ .)
#quartz()
plot4 <- ggplot(scores[scores$test == "test2",], aes(x = factor(pool), y = R2, fill = factor(pool))) + geom_boxplot(notch = T) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("") + ylab("Adjusted R2 (RF)") + scale_y_continuous(limits = c(-0.1,1)) + 
			facet_grid(factor(grouping) ~ .)
# Save plots 
ggsave(plot = plot1, filename = "boxplot_R2_pools_test1_zoo.jpg", dpi = 300, height = 4, width = 4)
ggsave(plot = plot2, filename = "boxplot_R2_pools_test2_zoo.jpg", dpi = 300, height = 4, width = 4)
ggsave(plot = plot3, filename = "boxplot_R2_pools_test1_zoo_groups.jpg", dpi = 300, height = 12, width = 4)
ggsave(plot = plot4, filename = "boxplot_R2_pools_test2_zoo_groups.jpg", dpi = 300, height = 12, width = 4)


### Distribution of variables' normalized ranks
#quartz()
plot1 <- ggplot(scores[scores$test == "test1",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = T) + 
			scale_fill_brewer(name = "", palette = "Set3") + theme_bw() +
			xlab("Variables") + ylab("Variable ranking (GLM)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) 
#quartz()
plot2 <- ggplot(scores[scores$test == "test2",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = T) + 
			scale_fill_brewer(name = "", palette = "Set3") + theme_bw() +
			xlab("Variables") + ylab("Variable ranking (RF)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) 
# Save
ggsave(plot = plot1, filename = "boxplot_ranks_test1_zoo.jpg", dpi = 300, height = 4, width = 6)
ggsave(plot = plot2, filename = "boxplot_ranks_test2_zoo.jpg", dpi = 300, height = 4, width = 6)

#quartz()
plot1 <- ggplot(scores[scores$test == "test1",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = F) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("Variable") + ylab("Variable ranking (GLM)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) + facet_grid(factor(grouping) ~ factor(pool) )
#
#quartz()
plot2 <- ggplot(scores[scores$test == "test2",], aes(x = factor(var), y = norm, fill = factor(var))) + geom_boxplot(notch = F) + 
			scale_fill_manual(name = "", values = getPalette(colourCount)) + theme_bw() +
			xlab("Variable") + ylab("Variable ranking (RF)") + scale_y_continuous(limits = c(0,1)) +
			theme(axis.text.x = element_text(angle = 90)) + facet_grid(factor(grouping) ~ factor(pool) )
# Save
ggsave(plot = plot1, filename = "boxplot_ranks_test1_zoo_groups.jpg", dpi = 300, height = 8, width = 25)
ggsave(plot = plot2, filename = "boxplot_ranks_test2_zoo_groups.jpg", dpi = 300, height = 8, width = 25)



# -------------------------------------------------------------

### 6°) For appendix 7: distribution of TSS scores

WD <- getwd()
setwd( paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future/", sep = "") )
zoo.wd <- getwd()
setwd( paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future/", sep = "") )
phyto.wd <- getwd()
setwd(WD)

# Vector of SDMs
SDMs <- c('GLM','GAM','RF','ANN')
# Vector of pools
pools <- c("p1","p2","p3","p4")
# p <- "p4"

### 2°) Plot distrbution of models'TSS values for zoo and phyto separately, using boxplots per SDM and facet per pool
setwd(zoo.wd)
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(zoo.wd,"/eval_scores_",p,"/", sep = "") )
					zoo_spp <- str_replace_all(dir(), "eval_scores_", "")
					zoo_spp <- str_replace_all(zoo_spp, ".Rdata", "")
					scores <- lapply(zoo_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind, scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}
) # eo lapply
# Rbind
table_scores_zoo <- do.call(rbind,res)
head(table_scores_zoo) # No SDM indicator...
table_scores_zoo$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_zoo)), pattern = "_", n = 2)))[,1])
table_scores_zoo$kingdom <- "Zooplankton"
dim(table_scores_zoo); head(table_scores_zoo)
summary(table_scores_zoo)
rm(res)
#table_scores_zoo[is.na(table_scores_zoo$TSS),]


### And for phyto
setwd(phyto.wd)
res <- lapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(phyto.wd,"/eval_scores_",p,"/", sep = "") )
					phyto_spp <- str_replace_all(dir(), "eval_scores_", "")
					phyto_spp <- str_replace_all(phyto_spp, ".Rdata", "")
					scores <- lapply(phyto_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind, scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}
) # eo lapply
# Rbind
table_scores_phyto <- do.call(rbind,res)
table_scores_phyto$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_phyto)), pattern = "_", n = 2)))[,1])
table_scores_phyto$kingdom <- "Phytoplankton"
dim(table_scores_phyto); head(table_scores_phyto)
summary(table_scores_phyto)
rm(res)
# table_scores_phyto[is.na(table_scores_phyto$TSS),]

### Join both tables and make plot
str(table_scores_phyto)
str(table_scores_zoo)
colnames(table_scores_phyto); colnames(table_scores_zoo)
scores <- rbind(table_scores_phyto,table_scores_zoo)
summary(scores)
rm(table_scores_phyto,table_scores_zoo);gc()

### Make plots (box plots per SDM, facet per kingdomxpools)
# Change name of pools to Set 1, Set 2 etc.
scores$pool <- factor(scores$pool)
levels(scores$pool)[levels(scores$pool) == "p1"] <- "Set 1"
levels(scores$pool)[levels(scores$pool) == "p2"] <- "Set 2"
levels(scores$pool)[levels(scores$pool) == "p3"] <- "Set 3"
levels(scores$pool)[levels(scores$pool) == "p4"] <- "Set 4"
levels(scores$pool) # ok

# Make plot
setwd(WD)
bp <- ggplot(scores, aes(x=factor(pool), y= TSS) ) + geom_boxplot(aes(fill = factor(kingdom))) + 
	 	scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + theme_light() + xlab("") + 
		ylab("True Skill Statistics (TSS)") + scale_y_continuous(limits = c(0,1)) + 
		facet_grid(factor(SDM) ~ factor(kingdom) ) 
# save
ggsave(plot = bp, filename = "boxplots_TSS_Appendix7.jpg", dpi = 300, width = 7, height = 13)

bp <- ggplot(scores, aes(x=factor(pool), y= ROC) ) + geom_boxplot(aes(fill = factor(kingdom))) + 
	 	scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + theme_light() + xlab("") + 
		ylab("Area Under the Curve (AUC)") + scale_y_continuous(limits = c(0,1)) + 
		facet_grid(factor(SDM) ~ factor(kingdom) ) 
# save
ggsave(plot = bp, filename = "boxplots_AUC_Appendix7.jpg", dpi = 300, width = 7, height = 13)



