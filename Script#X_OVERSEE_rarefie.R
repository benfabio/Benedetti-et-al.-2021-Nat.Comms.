

##### 16/04/2019 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- loading the v3.1v5.1 dataset of zooplankton occurrences (fitted with env)
#	- examine patterns of sampling effort and diversity effort (= sampled richness) per bands of 2° or 5° latitude
#	- find a way to rarefy the bands until you obtain a rarefied dataset with no dip in sampled diversity near the equator 
#	- separate data set into species darta and re-draw psAbs from it
#	- test ENMs (pool of vars: SST, dSST, logChl, logNO3, dO2)
 
### Last update : 07/11/2019

# -------------------------------------------------------------------------------------------------------------------

library("dplyr")
library("stringr")
library("reshape2")
library("tidyverse")
library("viridis")
library("RColorBrewer")
library("scales")

# -------------------------------------------------------------------------------------------------------------------

# Main dir() - Desktop
WD <- getwd()

# Go load coastline
setwd("/Users/fabiobenedetti/Desktop/~work/PostDocs/ETHZ/OVERSEE/data/")
cl <- read.csv("world_coast.csv", h = T)
coast <- list(
 	   # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	 	geom_polygon(aes(x = lon, y = lat), data = cl, fill = "grey55"),
  	  	geom_path(aes(x = lon, y = lat), data = cl, colour = "black", linetype = 1),
  	  	# appropriate projection
  	  	coord_quickmap(),
  	  	# remove extra space around the coast
  	  	scale_x_continuous(name = "Longitude", 
                     breaks = c(-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180), 
                     labels = c("-180°E","-150°E","-120°E","-90°E","-60°E","-30°E","0°E","30°E","60°E","90°E","120°E","150°E","180°E"), 
                     expand = c(0,0)), 
  
  		scale_y_continuous(name = "Latitude", 
                     breaks = c(-90,-60,-30,0,30,60,90), 
                     labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"),
                     expand = c(0,0)),
  		# dark gray background for the panel and legend
  	  	theme(
    		panel.background = element_rect(fill = "white"),  # background
    		legend.key = element_rect(fill = "grey50"),
    		panel.grid.major = element_line(colour = "grey70")
  		 )
) 
setwd(WD)

### Load the data, restrict to obs == 1
data <- get(load("zoo_data_fitted.Rdata"))
data <- data[data$obs == 1,]

# Define 2° and then 5° lat bands (like Menegotto & Rangel 2018) 
data$y_2d <- ceiling(data$y/2) * 2
data$y_5d <- ceiling(data$y/5) * 5
head(data)

# Plot occurrences density per band and observed richness 
quartz()
ggplot(data = data, aes(y_5d, stat(count))) + geom_density(position = "stack", fill = "grey55", clour = "black") + 
	ylab("Density") + xlab("Latitude (5°)") + theme_bw() 

# Histogram version
quartz()
ggplot(data, aes(x = y_5d)) + geom_histogram(colour="black", fill="grey65") + xlab("Latitude (5°)") + ylab("Count") + theme_bw() 
quartz()
ggplot(data, aes(x = y_2d)) + geom_histogram(colour="black", fill="grey65") + xlab("Latitude (2°)") + ylab("Count") + theme_bw() 
# Density v2
quartz()
p <- ggplot() + geom_density(aes(x = y_5d), data = data, color="black", fill = "grey55") + 
    xlab("Latitude") + ylab("Density of points") + theme_bw()

ggsave(plot = p, filename = "plot_density_sampling_effort_5d.jpg", dpi = 300, width = 5, height = 3)

# Plot total richness per band 
table2d <- data.frame( data %>% group_by(y_2d) %>% summarise(n = n(), rich = length(unique(species))) )
table5d <- data.frame( data %>% group_by(y_5d) %>% summarise(n = n(), rich = length(unique(species))) )
# Plot relationships
quartz()
ggplot(table2d, aes(x = n, y = rich)) + geom_point(colour = "black") + xlab("Effort") + ylab("Sampled richness")
# --> different rates of species accumulation between tropics and extratropics

quartz()
ggplot(table2d, aes(x = y_2d, y = rich)) + geom_point(colour = "black") + xlab("Latitude (2°)") + ylab("Sampled richness")
quartz()
ggplot(table5d, aes(x = y_5d, y = rich)) + geom_point(colour = "black") + xlab("Latitude (5°)") + ylab("Sampled richness")
# Peaks in sampled diversity = -25°N & + 35-40°N with a ckear dip between -6 and -16°N. 


### 17/04/2019: the 'dip' in sampled richness corresponds to a dip in sampling effort with n oscillating between 500 and 1000
### First, you should assess the minimum rarefaction sizes (size argument in the sample() function) that allows you to obtain 
### a flater pattern in diversity: 
#  - explore the distribution of sampled richness with varying rarefaction sizes (200,300,400,500,576)
#  - repeat process 30 times
#  - plot latitudinal patterns of rarefied sampled richness with boxplots and with facet_wrap() ~ rarefaction factor
#  - NOTE: at y_2d %in% c(-76,-74,90) you may not have enough occurrences to properly perform this test without  
#  - Do it for 5° bands too, with higher sizes (since more obs per band): c(500,1000,1250,1500,1750,2000)

sizes <- c(200,300,400,500,576)

data2rarefy <- data[-which(data$y_2d %in% c(-76,-74,90)),]

# s <- 500
# b <- "-8"
require("parallel")
rarefied <- mclapply(unique(data2rarefy$y_2d), function(b) {
 				
				# Select the bin you're interested in
				message(paste("rarefying obs for band || ", b, sep = ""))
 				d <- data2rarefy[which(data2rarefy$y_2d == b),]
 				
				# Randomly re-sample s times (s in sizes), and do so 30 times 
				rar <- lapply(sizes, function(s) {
							
							# Random re-sampling - 30 times
							rarar <- lapply(c(1:30), function(r) {
										dd <- d[sample(1:nrow(d), size = s),]
										dd$run <- r
										return(dd)
									}
							) # eo 3rd lapply
							
							# Rbind and return
							new_d <- do.call(rbind, rarar)
							new_d$size_rar <- s
							rm(rarar)
							return(new_d)
						} 
				) # eo 2nd lapply
				
				# Rbind
				new_d <- do.call(rbind, rar)
				rm(rar)
				# length( unique(new_d$species) )
				# length( unique( d$species) )
				# Return new_d
 				return(new_d)
 		   } # eo 1st FUN
		   , mc.cores = 2  
		   
) # eo 1st lapply

# Examine results
# str(rarefied)
# length(rarefied)
new_d <- do.call(rbind, rarefied)
dim(new_d)
str(new_d)
rm(rarefied)

### Plot distrbution of sampled diversity per band and per rarefaction size
require("dplyr")
div <- data.frame( new_d %>% group_by(y_2d, size_rar, run) %>% summarise(rich = length(unique(species)), n = n()) ) 
head(div)
summary(div)
# Check if n obs do not vary per band 
ggplot(div, aes(x = factor(y_2d), y = n)) + 
	geom_boxplot(fill = "grey50", colour = "black") + 
	xlab("Latitude (2°)") + ylab("Count") + theme_bw() + 
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	facet_wrap(factor(size_rar) ~ ., ncol = 2)	
# gut gut
	
#quartz()
plot <- ggplot(div, aes(x = factor(y_2d), y = rich)) + 
	geom_boxplot(fill = "grey50", colour = "black") + 
	xlab("Latitude (2°)") + ylab("Sampled richness") + theme_bw() + 
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	facet_wrap(factor(size_rar) ~ ., ncol = 2)	
# Nice, do the same with 5° bands
ggsave(plot = plot, filename = "plot_distrib_div_rar_2d.jpg", dpi = 300, width = 13, height = 10)

div$raw_rich <- NA # fill w/ for loop
for(b in unique(div$y_2d)) {
		div[div$y_2d == b,"raw_rich"] <- table2d[table2d$y_2d == b,"rich"]
} # eo for loop

# And compute percentage
div$perc <- (div$rich/div$raw_rich)*100

# And plot
plot <- ggplot(div, aes(x = factor(y_2d), y = perc)) + 
	geom_boxplot(fill = "grey50", colour = "black") + 
	xlab("Latitude (5°)") + ylab("% sampled richness") + theme_bw() + 
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	facet_wrap(factor(size_rar) ~ ., ncol = 2)	
#
ggsave(plot = plot, filename = "plot_distrib_perc_rar_2d.jpg", dpi = 300, width = 13, height = 10)



### ------------------------------------------------------------

sizes <- c(500,1000,1250,1500,1750,2000)

data2rarefy <- data[-which(data$y_5d %in% c(-75,90)),]

# s <- 500
# b <- "-8"
require("parallel")
rarefied <- mclapply(unique(data2rarefy$y_5d), function(b) {
				# Select the bin you're interested in
				message(paste("rarefying obs for band || ", b, sep = ""))
 				d <- data2rarefy[which(data2rarefy$y_5d == b),]
				# Randomly re-sample s times (s in sizes), and do so 30 times 
				rar <- lapply(sizes, function(s) {
							# Random re-sampling - 30 times
							rarar <- lapply(c(1:30), function(r) {
										dd <- d[sample(1:nrow(d), size = s),]
										dd$run <- r
										return(dd)
									}
							) # eo 3rd lapply
							# Rbind and return
							new_d <- do.call(rbind, rarar)
							new_d$size_rar <- s
							rm(rarar)
							return(new_d)
						} 
				) # eo 2nd lapply
				# Rbind
				new_d <- do.call(rbind, rar)
				rm(rar)
 				return(new_d)
 		   } # eo 1st FUN
		   , mc.cores = 2  
) # eo 1st lapply
# Examine results
new_d <- do.call(rbind, rarefied)
dim(new_d)
str(new_d)
rm(rarefied)

### Plot distrbution of sampled diversity per band and per rarefaction size
require("dplyr")
div <- data.frame( new_d %>% group_by(y_5d, size_rar, run) %>% summarise(rich = length(unique(species)), n = n()) ) 
head(div)
summary(div)

#quartz()
plot <- ggplot(div, aes(x = factor(y_5d), y = rich)) + 
	geom_boxplot(fill = "grey50", colour = "black") + 
	xlab("Latitude (5°)") + ylab("Sampled richness") + theme_bw() + 
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	facet_wrap(factor(size_rar) ~ ., ncol = 2)	
#
ggsave(plot = plot, filename = "plot_distrib_div_rar_5d.jpg", dpi = 300, width = 13, height = 10)

### Also check the % of richness covered by the rarefied data compared to the unrarefied one (table5d)
div$raw_rich <- NA # fill w/ for loop
for(b in unique(div$y_5d)) {
		div[div$y_5d == b,"raw_rich"] <- table5d[table5d$y_5d == b,"rich"]
} # eo for loop

# And compute percentage
div$perc <- (div$rich/div$raw_rich)*100

# And plot
plot <- ggplot(div, aes(x = factor(y_5d), y = perc)) + 
	geom_boxplot(fill = "grey50", colour = "black") + 
	xlab("Latitude (5°)") + ylab("% sampled richness") + theme_bw() + 
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	facet_wrap(factor(size_rar) ~ ., ncol = 2)	
#
ggsave(plot = plot, filename = "plot_distrib_perc_rar_5d.jpg", dpi = 300, width = 13, height = 10)

# And plot density to make sure sampling effort is even :p
ggplot(new_d, aes(x = factor(y_2d) )) + 
	geom_density(fill = "grey50", colour = "black", position = "stack") + 
	xlab("Latitude (2°)") + ylab("Density") + theme_bw() + 
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	facet_wrap(factor(size_rar) ~ ., ncol = 2)	
	
# summary(factor( new_d[new_d$run == 1 & new_d$size_rar == 200,] ))

### To help choose between rarefaction levels: map sampling effort per size_rar ? 
# length(unique(new_d[new_d$size_rar == 2000 & new_d$run == 3,"species"]))
quartz()
ggplot() + geom_point(aes(x = x, y = y), data = new_d[new_d$size_rar == 2000 & new_d$run == 3,]) + 
	coast + coord_quickmap() + theme_bw()


### OK, same with 1°x1° data, map effort and richness in the rarefied dataset
rar_ddf <- new_d[new_d$size_rar == 2000 & new_d$run == 25,]
colnames(rar_ddf); colnames(data[which(data$y_5d %in% c(-75,90)),])
# And add the data from the bands you omitted earlier 
rar_ddf <- rbind( data[which(data$y_5d %in% c(-75,90)),] , rar_ddf[,c(1:39)] )
dim(rar_ddf)
head(rar_ddf)

rar_ddf$id <- factor( paste(round(rar_ddf$x, .1), round(rar_ddf$y, .1), sep = "_") )

effort <- data.frame( rar_ddf %>% group_by(id) %>% summarise(x = unique(round(x,.1)), y = unique(round(y,.1)), n = n(), rich = length(unique(species))) )
summary(effort)
# dim(effort)

### Map rarefied sampling effort
quartz()
ggplot() + geom_raster(aes(x = x, y = y, fill = log(n)), data = effort) + 
	scale_fill_viridis(name = "Effort - logged") + coast + coord_quickmap() + theme_bw()
#
quartz()
ggplot() + geom_point(aes(x = y, y = n), data = effort) + xlab("Latitude") + ylab("Effort") + theme_bw()

### Map rarefied diversity
quartz()
ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = effort) + 
	scale_fill_viridis(name = "Richness") + coast + coord_quickmap() + theme_bw()

quartz()
ggplot() + geom_point(aes(x = y, y = rich), data = effort) + xlab("Latitude") + ylab("Richness") + theme_bw()

### OK fine, one last step before saving and draxing psAbs: look at n counts per species 
counts <- data.frame( rar_ddf %>% group_by(species) %>% summarise(n = n() ) )
counts[order(counts$n),]
nrow(counts[counts$n > 50,]) # 257 species to model

write.table(rar_ddf, file = "rarefied_data_zooplankton_17_04_19.txt", sep = ";")

ddf <- read.table("rarefied_data_zooplankton_17_04_19.txt", h = T, sep = ";")
# head(ddf)
effort <- data.frame( ddf %>% group_by(id) %>% summarise(x = unique(round(x,.1)), y = unique(round(y,.1)), n = n(), rich = length(unique(species))) )
summary(effort)
# Save plot of sampling effort
#p <- ggplot() + geom_point(aes(x = y, y = n), data = effort) + xlab("Latitude") + ylab("Rarefied sampling effort") + theme_bw()
p <- ggplot() + geom_density(aes(x = y_5d), data = ddf, color="black", fill="grey55") + xlab("Latitude") + ylab("Density of points") + theme_bw()
ggsave(plot = p, filename = "plot_density_sampling_effort_rar.jpg", dpi = 300, width = 5, height = 3)

#quartz()
#ggplot() + geom_point(aes(x = y, y = rich), data = effort) + xlab("Latitude") + ylab("Richness") + theme_bw()



