
##### 23/11/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for preparing the Revisions of the manuscript submitted to Nat Comms.:

#	x Gather baseline and future changes in diversity (ensemble) for plankton groups (all, phyto-/zoo- and PFTs you can define) + services proxies
#   x Examine correlations in baseline conditions between services ans groups SR, at 'bioregions' scales and glibal and lat bands
#   x Plot correlation networks with 'corrr'
#   x PLot and quantify correlations between phyto div and nutrients fileds --> test nutrients use efficiency?
#   x Compute ensemble map of copepod median size (for each pool*SDM and then pool*SDM*ESM) and relate to services again (global/regional)
#   x Compute ensemble map of Diatoms and Dinos median cell diameter (biovolume?)
#   x Compute mean/median diatom bioviolume based on he MARDEAT data (Leblanc et al., 2012)

#   x For revisions #2: check longitudinal patterns in LDGs and changes in SR
#   - For revisions #2: analyze overlap between present projections and those of Ibarbalz et al. (2019)

# Note: x = done; - = to do

### Last update: 18/05/21

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("tidyverse")
library("stringr")
library("reshape2")
library("RColorBrewer")
library("scales")
library("maps")
library("ggthemes")
library("viridis")

world2 <- map_data(map = "world2")
WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Gather plankton groups diversity projections with clusters of changes and proxy variables for ES
setwd("/net/kryo/work/fabioben/OVERSEE/data")
data <- get(load("table_ann_changes_diversity+services_ensemble_klusters_21_04_20.Rdata"))
# dim(data) ; colnames(data) ; str(data)
# Ok, this object has: the regions and the proxies of ES and the large plankton diversity patterns

### Then, need to retrieve the group level patterns (those used in Fig. S1 for instance --> Script called groups ordinate.R)
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology")
ddf <- get(load("table_longformat_group_rich_29_01_20.Rdata"))
# Dcast per kingdom and then per group
# case 1: kingdoms
cast1 <- dcast(data = ddf[,c(1,2,4:6)], cell_id + x + y ~ kingdom, fun.aggregate = sum, value.var = "rich_base")
dim(cast1) ; head(cast1) ; summary(cast1)
# case 2: main groups
cast2 <- dcast(data = ddf[,c(1,3,4:6)], cell_id + x + y ~ group, fun.aggregate = sum, value.var = "rich_base")
head(cast2) ; summary(cast2)
# Re-provide kingdoms to cast2
cast2$Phytoplankton <- cast1$Phytoplankton
cast2$Zooplankton <- cast1$Zooplankton
colnames(cast2)
rm(cast1) ; gc()

### And do the same based on differences in richness ! 
ddf$diff <- (ddf$rich_fut)-(ddf$rich_base)
# summary(ddf$diff)
# Dcast
cast_fut1 <- dcast(data = ddf[,c(1,2,4,5,8)], cell_id + x + y ~ kingdom, fun.aggregate = sum, value.var = "diff")
# case 2: main groups
cast_fut2 <- dcast(data = ddf[,c(1,3,4,5,8)], cell_id + x + y ~ group, fun.aggregate = sum, value.var = "diff")
# Re-provide kingdoms to cast2
cast_fut2$Phytoplankton <- cast_fut1$Phytoplankton
cast_fut2$Zooplankton <- cast_fut1$Zooplankton
colnames(cast_fut2)

### Quick summary: 
# - 'data' has the ES proxiyes and the regions
# - 'cast2' has the baseline groups' diversity
# - 'cast_fut2' has the fuetur groups' diversity

# Make sure they follow same order
data <- data[order(data$id),]
cast2 <- cast2[order(cast2$cell_id),]
cast_fut2 <- cast_fut2[order(cast_fut2$cell_id),]
dim(data) ; dim(cast2) ; dim(cast_fut2)
# length(intersect(cast2$cell_id, cast_fut2$cell_id)) # cast2 & cast_fut2 have the same ids, so they're already on the same format
length(intersect(data$id, cast_fut2$cell_id))
commons <- intersect(data$id, cast_fut2$cell_id)
data2 <- data[data$id %in% commons,]
cast_base <- cast2[cast2$cell_id %in% commons,]
cast_fut <- cast_fut2[cast_fut2$cell_id %in% commons,]

# Put them in the same ddf (data2)
data.revi <- data2[,c("id","x","y","pca_euclid_pam_k6","jac","jtu_phyto","jtu_zoo","AllNorm","OceanNorm","FPOCex","NPP","e","slope2","catches")]
data.revi$phyto_base <- cast_base$Phytoplankton
data.revi$zoo_base <- cast_base$Zooplankton
data.revi$phyto_fut <- cast_fut$Phytoplankton
data.revi$zoo_fut <- cast_fut$Zooplankton
# And groups:
groups <- colnames(cast_fut)[c(5:7,9,11,13:16,18)] ; groups
data.revi[,paste(groups,"base", sep = "_")] <- cast_base[,groups]
data.revi[,paste(groups,"fut", sep = "_")] <- cast_fut[,groups]
# Check
# head(data.revi) ; summary(data.revi)


# ------------------------------------------------------------

### 2°) Examine relationships between the groups' diversity pattern (log maybe) and the ES proxy variables
### A. Overall plankton groups (phyto- and zooplankton) baseline diversity

ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = e)) +
     theme_classic() + ylab("e ratio") + xlab("Phytoplankton SR") + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")

ggplot(data = data.revi) + geom_point(aes(x = log(Bacillariophyceae_base), y = e)) +
   theme_classic() + ylab("NPP (log)") + xlab("Diatom SR") + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")

# Plot in a for loop to go faster (across ES variables)
vars <- colnames(data.revi)[c(9:14)] ; vars
# v <- "catches"
for(v in vars) {
    
    message(paste("Plotting B-ES relationship for : ",v, sep = ""))
    
    # For NPP, FPOC and fish catches, log transform
    if(v %in% c("FPOCex","NPP")) {
        
        p1 <- ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Phytoplankton SR") +
            theme_classic() + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")
    
        p2 <- ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Zooplankton SR") +
            theme_classic() + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")
        
    } else if (v == "catches") {
        
        p1 <- ggplot(data = data.revi[data.revi$catches > 5,]) + geom_point(aes(x = log(phyto_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Phytoplankton SR") +
            theme_classic() + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")
    
        p2 <- ggplot(data = data.revi[data.revi$catches > 5,]) + geom_point(aes(x = log(zoo_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Zooplankton SR") +
            theme_classic() + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")
        
    } else  {
        
        p1 <- ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = get(v)), alpha = .5) +
            ylab(paste(v, sep = "")) + xlab("Phytoplankton SR") +
            theme_classic() + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")
    
        p2 <- ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = get(v)), alpha = .5) +
            ylab(paste(v, sep = "")) + xlab("Zooplankton SR") +
            theme_classic() + facet_wrap(.~factor(pca_euclid_pam_k6), scales = "free")
        
    } # eo if else loop
        
    # ggsave
    setwd(WD)
    ggsave(plot = p1, filename = paste("plot_B-ES_revisions_",v,"_phyto_base_clusters.jpg", sep = ""), dpi = 300, width = 10, height = 5)
    ggsave(plot = p2, filename = paste("plot_B-ES_revisions_",v,"_zoo_base_clusters.jpg", sep = ""), dpi = 300, width = 10, height = 5)
        
} # eo for loop - v in vars

### Assess strength of linear relationships - report on table
summary(lm(log(FPOCex) ~ log(zoo_base), data = data.revi[data.revi$pca_euclid_pam_k6 == 6,]))

summary(lm(log(catches) ~ log(zoo_base), data = data.revi[data.revi$catches > 5 & data.revi$pca_euclid_pam_k6 == 6,]))

summary(lm(slope2 ~ log(zoo_base), data = data.revi[data.revi$pca_euclid_pam_k6 == 6,]))


### Also check global relationships ? Or across latitudinal bands

# 1. Phyto vs NPP/ FPOC/ catches
summary(lm(log(NPP) ~ log(phyto_base), data = data.revi)) # Adjusted R-squared:  0.4575 ! (positive)
summary(lm(log(FPOCex) ~ log(phyto_base), data = data.revi)) # no, R2 too small
summary(lm(log(catches) ~ log(phyto_base), data = data.revi[data.revi$catches > 5,])) # no, R2 too small

ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = log(NPP)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(phyto_base), y = log(NPP)), method = "lm", colour = "#e31a1c") +
    ylab("NPP (log)") + xlab("Phytoplankton SR (log)") + theme_classic()

# 2. Phyto vs. Ocean norm div/ e ratio/ slope2
summary(lm(OceanNorm ~ log(phyto_base), data = data.revi)) # R-squared: 0.6274 
summary(lm(e ~ log(phyto_base), data = data.revi)) # Adjusted R-squared:  0.5249 
summary(lm(slope2 ~ log(phyto_base), data = data.revi)) # no, r2 too weak

ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = OceanNorm), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(phyto_base), y = OceanNorm), method = "lm", colour = "#e31a1c") +
    ylab("Oceanic diversity (normalized)") + xlab("Phytoplankton SR (log)") + theme_classic()
#
ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = e), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(phyto_base), y = e), method = "lm", colour = "#e31a1c") +
    ylab("e ratio") + xlab("Phytoplankton SR (log)") + theme_classic()


# 3. Zooplankton vs NPP/ FPOC/ fish catches
summary(lm(log(NPP) ~ log(zoo_base), data = data.revi)) # R-squared:  0.321 
summary(lm(log(FPOCex) ~ log(zoo_base), data = data.revi)) # no, r2 too small
summary(lm(log(catches) ~ log(zoo_base), data = data.revi[data.revi$catches > 5,])) # R-squared:  0.1244 

ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = log(NPP)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(zoo_base), y = log(NPP)), method = "lm", colour = "#e31a1c") +
    ylab("NPP (log)") + xlab("Zooplankton SR (log)") + theme_classic()
#
ggplot(data = data.revi[data.revi$catches > 5,]) + geom_point(aes(x = log(zoo_base), y = log(catches)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(zoo_base), y = log(catches)), method = "lm", colour = "#e31a1c") +
    ylab("Catches (log)") + xlab("Zooplankton SR (log)") + theme_classic()


# 4. Zooplankton vs. Ocean norm div/ e ratio/ slope2
summary(lm(OceanNorm ~ log(zoo_base), data = data.revi)) # R-squared:  0.7245 
summary(lm(e ~ log(zoo_base), data = data.revi)) # R-squared:  0.5773 
summary(lm(slope2 ~ log(zoo_base), data = data.revi)) # R-squared:  0.2036 

ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = OceanNorm), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(zoo_base), y = OceanNorm), method = "lm", colour = "#e31a1c") +
    ylab("Oceanic diversity (normalized)") + xlab("Zooplankton SR (log)") + theme_classic()
#
ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = e), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(zoo_base), y = e), method = "lm", colour = "#e31a1c") +
    ylab("e ratio") + xlab("Zooplankton SR (log)") + theme_classic()
#
ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = slope2), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(zoo_base), y = slope2), method = "lm", colour = "#e31a1c") +
    ylab("PSD") + xlab("Zooplankton SR (log)") + theme_classic()


### And across latitudinal bands/ climate regime (Tropical/ Temperate/ Polar)
# Add them to data.revi
data.revi$band <- NA
data.revi[abs(data.revi$y) < 30,"band"] <- "Tropical"
data.revi[abs(data.revi$y) >= 30,"band"] <- "Temperate"
data.revi[abs(data.revi$y) > 60,"band"] <- "Polar"
summary(factor(data.revi$band))

for(v in vars) {
    
    message(paste("Plotting B-ES relationship for : ",v, sep = ""))
    
    # For NPP, FPOC and fish catches, log transform
    if(v %in% c("FPOCex","NPP","catches")) {
        
        p1 <- ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Phytoplankton SR") +
            theme_classic() + facet_wrap(.~factor(band), scales = "free")
    
        p2 <- ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Zooplankton SR") +
            theme_classic() + facet_wrap(.~factor(band), scales = "free")
        
    } else if(v == "catches") {
        
        p1 <- ggplot(data = data.revi[data.revi$catches > 5,]) + geom_point(aes(x = log(phyto_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Phytoplankton SR") +
            theme_classic() + facet_wrap(.~factor(band), scales = "free")
    
        p2 <- ggplot(data = data.revi[data.revi$catches > 5,]) + geom_point(aes(x = log(zoo_base), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Zooplankton SR") +
            theme_classic() + facet_wrap(.~factor(band), scales = "free")
        
    } else {
        
        p1 <- ggplot(data = data.revi) + geom_point(aes(x = log(phyto_base), y = get(v)), alpha = .5) +
            ylab(paste(v, sep = "")) + xlab("Phytoplankton SR") +
            theme_classic() + facet_wrap(.~factor(band), scales = "free")
    
        p2 <- ggplot(data = data.revi) + geom_point(aes(x = log(zoo_base), y = get(v)), alpha = .5) +
            ylab(paste(v, sep = "")) + xlab("Zooplankton SR") +
            theme_classic() + facet_wrap(.~factor(band), scales = "free")
        
    } # eo if else loop
        
    # Ggsave
    setwd(WD)
    ggsave(plot = p1, filename = paste("plot_B-ES_revisions_",v,"_phyto_base_bands.jpg", sep = ""), dpi = 300, width = 10, height = 3)
    ggsave(plot = p2, filename = paste("plot_B-ES_revisions_",v,"_zoo_base_bands.jpg", sep = ""), dpi = 300, width = 10, height = 3)
        
} # eo for loop - v in vars


### Assess strength of linear relationships - report on table
summary(lm(log(FPOCex) ~ log(zoo_base), data = data.revi[data.revi$band == "Polar",])) 

summary(lm(log(catches) ~ log(zoo_base), data = data.revi[data.revi$catches > 5 & data.revi$band == "Tropical",])) 

summary(lm(slope2 ~ log(zoo_base), data = data.revi[data.revi$band == "Polar",])) 
# cor(data.revi[data.revi$band == "Polar","e"], log(data.revi[data.revi$band == "Polar","zoo_base"]), method = "spearman")

### And re-plot the distribution of future changes in SR across the 3 lat bands for phyto- and zooplankton & help interpret results above
ggplot(data = data.revi, aes(x = factor(band), y = phyto_fut)) + geom_boxplot(colour = "black", fill = "grey70") +
    geom_hline(yintercept = 0, linetype = "dashed") + xlab("") + ylab("Difference in Phytoplankton richness") +
    theme_classic() 

ggplot(data = data.revi, aes(x = factor(band), y = zoo_fut)) + geom_boxplot(colour = "black", fill = "grey70") +
    geom_hline(yintercept = 0, linetype = "dashed") + xlab("") + ylab("Difference in Zooplankton richness") +
    theme_classic() 

# --> phyto SR increases by >15% in tropics and > 5% in temperate (but high turn-over) ; no increase in poles
# --> zoo SR decreases in the tropics, increases slightly in temperate regimes (turn-over) and does not really change in poles 



### B. Smaller plankton groups baseline diversity  ------------------------------------------------------------
colnames(data.revi)

summary(lm(slope2 ~ log(Chordata_base), data = data.revi))

summary(lm(log(catches) ~ log(Chordata_base), data = data.revi[data.revi$catches > 5,])) 

ggplot(data = data.revi[data.revi$catches > 5,]) + geom_point(aes(x = log(Copepoda_base), y = log(catches)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(Copepoda_base), y = log(catches)), method = "lm", colour = "#e31a1c") +
    ylab("Fish catches (log)") + xlab("Richness (log)") + theme_classic()
    
### Melt data.revi to put groups' baseline SR as column and plot B-ES relationships per cluster with faceting (clusterxgroup)
m.data.revi <- melt(data.revi[,c(1:14,19:28,39)], id.vars = c("id","x","y","band","pca_euclid_pam_k6","jac","jtu_phyto","jtu_zoo","AllNorm","OceanNorm","FPOCex","NPP","e","slope2","catches"))
head(m.data.revi)
summary(m.data.revi)
colnames(m.data.revi)[c(16,17)] <- c("group","rich")
unique(m.data.revi$group) # remove the "base"
m.data.revi$group <- str_replace_all(m.data.revi$group, "_base", "")

vars <- colnames(data.revi)[c(9:14)] ; vars
v <- "NPP"

for(v in vars) {
    
    message(paste("Plotting B-ES relationships across groups & clusters of severity for : ",v, sep = ""))
    
    # For NPP, FPOC and fish catches, log transform
    if(v %in% c("FPOCex","NPP","catches")) {
        
        p1 <- ggplot(data = m.data.revi) + geom_point(aes(x = log(rich), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = "")) + xlab("Richness (log)") +
            theme_bw() + facet_grid(factor(group) ~ factor(pca_euclid_pam_k6), scales = "free")
        
    } else if(v == "catches") {
        
        p1 <- ggplot(data = m.data.revi[m.data.revi$catches > 5,]) + geom_point(aes(x = log(rich), y = log(get(v))), alpha = .5) +
            ylab(paste(v," (log)", sep = ""))  + xlab("Richness (log)") +
            theme_bw() + facet_grid(factor(group) ~ factor(pca_euclid_pam_k6), scales = "free")
        
    } else {
        
        p1 <- ggplot(data = m.data.revi) + geom_point(aes(x = log(rich), y = get(v)), alpha = .5) +
            ylab(paste(v, sep = "")) + xlab("Richness (log)") +
            theme_bw() + facet_grid(factor(group) ~ factor(pca_euclid_pam_k6), scales = "free")
        
    } # eo if else loop
        
    # Ggsave
    setwd(WD)
    ggsave(plot = p1, filename = paste("plot_B-ES_revisions_",v,"_groups_clusters.jpg", sep = ""), dpi = 300, width = 15, height = 13)
        
} # eo for loop - v in vars

### Assess lm() and report on dedicated excel sheet
colnames(data.revi)
# summary(lm(slope2 ~ log(Chordata_base), data = data.revi[data.revi$pca_euclid_pam_k6 == 6,]))
summary(lm(log(catches) ~ log(Chordata_base), data = data.revi[data.revi$catches > 5 & data.revi$pca_euclid_pam_k6 == 6,])) 


### And plot the distribution of %∆SR for groups and clusters to help you interpret your results in terms of changes in ES
colnames(data.revi)
m.data.revi <- melt(data.revi[,c(1:4,19:38)], id.vars = c("id","x","y","pca_euclid_pam_k6"))
#head(m.data.revi)
#summary(m.data.revi)
colnames(m.data.revi)[c(5,6)] <- c("group","rich")
head(do.call(rbind,strsplit(as.character(m.data.revi$group), split = "_")))
m.data.revi$period <- do.call(rbind,strsplit(as.character(m.data.revi$group), split = "_"))[,2]
m.data.revi$group <- str_replace_all(m.data.revi$group, "_base", "")
m.data.revi$group <- str_replace_all(m.data.revi$group, "_fut", "")
unique(m.data.revi$group) # 

# Dcast to put base and fut in separate columns
d.data.revi <- dcast(m.data.revi, formula = id+x+y+pca_euclid_pam_k6+group ~ period, value.var = "rich")
#head(d.data.revi)

### Derive % changes 
perc <- data.frame(d.data.revi %>% group_by(id, group) %>% summarize(x = unique(x), y = unique(y), pca_euclid_pam_k6 = unique(pca_euclid_pam_k6), perc = ((fut/base))*100 ) ) 
summary(perc)

### !!! Clusters number do not match those in Fig. 4 yet --> change
perc$clusters <- perc$pca_euclid_pam_k6
perc[perc$pca_euclid_pam_k6 == 1,"clusters"] <- 4
perc[perc$pca_euclid_pam_k6 == 2,"clusters"] <- 3
perc[perc$pca_euclid_pam_k6 == 3,"clusters"] <- 6
perc[perc$pca_euclid_pam_k6 == 4,"clusters"] <- 2
perc[perc$pca_euclid_pam_k6 == 6,"clusters"] <- 1

### Re-label some groups for figures consistency
perc$group2 <- perc$group
unique(perc$group) ; unique(perc$group2)
perc[perc$group == "Bacillariophyceae","group2"] <- "Diatoms"
perc[perc$group == "Chaetognatha","group2"] <- "Chaetognaths"
perc[perc$group == "Chordata","group2"] <- "Chordates"
perc[perc$group == "Copepoda","group2"] <- "Copepods"
perc[perc$group == "Dinoflagellata","group2"] <- "Dinoflagellates"
perc[perc$group == "Pteropoda","group2"] <- "Pteropods"

perc$group2 <- factor(perc$group2, levels = c("Diatoms","Dinoflagellates","Haptophyta",
                                    "Copepods","Malacostraca","Jellyfish","Chordates",
                                    "Chaetognaths","Pteropods","Foraminifera"))

plot <- ggplot(data = perc[perc$perc < 100,], aes(x = factor(clusters), y = perc, fill = factor(clusters))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") + 
    ylab(paste("% difference in mean annual richness", sep = ""))  + xlab("") + theme_bw() +
    facet_wrap(.~factor(group2), scales = "free_y", nrow = 3) + geom_hline(yintercept = 0, linetype = "dashed")
#
ggsave(plot = plot, filename = "plot_distrib_perc_groups_clusters_severity.jpg", dpi = 300, width = 10, height = 9)

# And global plot
plot <- ggplot(data = perc[perc$perc < 100,], aes(x = group2, y = perc, fill = factor(group2))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Groups", palette = "RdYlBu") + 
    ylab(paste("% difference in mean annual richness", sep = ""))  + xlab("") + theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = plot, filename = "plot_distrib_perc_groups_global.jpg", dpi = 300, width = 6, height = 6)

# And plot perc per cluster 1 group
plot <- ggplot(data = perc[perc$perc < 100,], aes(x = factor(clusters), y = perc, fill = factor(clusters))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Groups", palette = "RdYlBu") + 
    ylab(paste("% difference in mean annual richness", sep = ""))  + xlab("") + theme_bw() +
    facet_wrap(.~factor(group2), scales = "free_y", nrow = 3) + geom_hline(yintercept = 0, linetype = "dashed")
    
ggsave(plot = plot, filename = "plot_distrib_perc_groups_clusters.jpg", dpi = 300, width = 10, height = 8)


# ------------------------------------------------------------

### 26/11/2020: Use 'corrr' to plot corr coeff (spearman's) in a network fashion.
### Go on pesonal computer for that because I couldn't install 'corrr' on kryo somehow.

library("tidyverse")
library("reshape2")
library("viridis")
library("RColorBrewer")
library("scales")
library("ggthemes")
library("corrr")
library("corrgram")
library("igraph")
library("ggraph")

WD <- getwd()

### Draw the correlation coeff heatmaps for each net
get_lower_tri <- function(cormat){
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
}
get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
}
reorder_cormat <- function(cormat){
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <- cormat[hc$order, hc$order]
        return(cormat)
}


data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
# colnames(data.revi)
names <- colnames(data.revi)[c(1,9:16,19,23,25,22,27,26,21,20,28,24)] ; names
mydata <- data.revi[,names]
# Change some colnames etc.
colnames(mydata)[2] <- "Oceanic\nbiodiversity"
colnames(mydata)[5] <- "e ratio"
colnames(mydata)[6] <- "PSI"
colnames(mydata)[7] <- "Fish catches"
colnames(mydata)[c(8:19)] <- str_replace_all(colnames(mydata)[c(8:19)] , "_base", "")
colnames(mydata)[c(8:11,13,16:18)] <- c("Phytoplankton","Zooplankton","Diatoms","Dinoflagellates","Copepods","Chordates","Chaetognaths","Pteropods")
head(mydata)

### Combine with annual climatologies of nutrients fields
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d")
files <- dir()[grep("_18_09_20",dir())] ; files
clims <- lapply(files, function(f) {
            d <- read.table(f, h = T, sep = ";")
            d$x2 <- d$x
            d[d$x < 0 ,"x2"] <- (d[d$x < 0 ,"x"]) + 360
            d$id2 <- paste(d$x2, d$y, sep = "_")
            d <- d[,c("id2","x2","y","SST","logChl","logNO3","logSiO2")]
            return(d)
    } # eo 
)
clims <- bind_rows(clims)

### Derive mean annual climatologies
ann.clim <- data.frame(clims %>% group_by(id2) %>%
                summarize(x = unique(x2), y = unique(y), SST = mean(SST, na.rm = T), 
                NO3 = mean(logNO3, na.rm = T), SiO2 = mean(logSiO2, na.rm = T), 
                Chlorophyll = mean(logChl, na.rm = T))
) # eo ddf
# summary(ann.clim)
setwd(WD)

# OK, re-order, find common cells and cbind
ann.clim <- ann.clim[order(ann.clim$id2),]
mydata <- mydata[order(mydata$id),]
commons <- intersect(ann.clim$id2, mydata$id)
mydata <- mydata[mydata$id %in% commons,]
ann.clim <- ann.clim[ann.clim$id2 %in% commons,]
mydata[,c("SST","NO3","SiO2","Chlorophyll")] <- ann.clim[,c("SST","NO3","SiO2","Chlorophyll")]

### Transform the columns needed
mydata$FPOCex <- log(mydata$FPOCex)
mydata$NPP <- log(mydata$NPP)
mydata[,"Fish catches"] <- log1p(mydata[,"Fish catches"])
# mydata[,c(7:length(mydata))] <- log1p(mydata[,c(7:length(mydata))])
head(mydata)

### Try some basic corr heatmap first 
cormat <- round(cor(na.omit( mydata[,c(2:7,20:23,8:19)] ), method = "spearman"),2)
#head(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = T)
colnames(melted_cormat) <- c("V1","V2","rho")
melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
groups <- colnames(mydata)[c(8:19)]
vars <- colnames(mydata)[c(2:7,20:23)]
colnames(melted_cormat) <- c("ES","Group","rho")
#melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]

max.val <- max(abs(melted_cormat$rho))

# looks good to go
plot <- ggplot(melted_cormat, aes(factor(Group), factor(ES), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0,
            limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group), factor(ES), label = rho), color = "black", size = 3)
#
ggsave(plot = plot, filename = "heatmap.jpg", dpi = 300, width = 8.5, height = 8.5)


### Re-supply clusters of severity
data.revi2 <- data.revi[order(data.revi$id),]
data.revi2 <- data.revi2[data.revi2$id %in% commons,]
mydata$clusters <- data.revi2$clusters  
head(mydata)

### Same as above but facet plot per clusters ! (use lapply before to retunr all cormats)
res <- lapply(unique(mydata$clusters), function(i) {
    
            message(paste("Computing corr matrox for cluster  ",i, sep = ""))
            
            ### Try some basic corr heatmap first 
            cormat <- round(cor(na.omit( mydata[mydata$clusters == i,c(2:7,20:23,8:19)] ), method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = T)
            colnames(melted_cormat) <- c("V1","V2","rho")
            melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
            groups <- colnames(mydata)[c(8:19)]
            vars <- colnames(mydata)[c(2:7,20:23)]
            colnames(melted_cormat) <- c("ES","Group","rho")
            melted_cormat$cluster <- i
            return(melted_cormat)
    
    } # eo lapply 
    
) # eo lapply - i in clusters
table.corr <- dplyr::bind_rows(res) ; rm(res) ; gc()
max.val <- max(abs(table.corr$rho)) ; max.val

plot <- ggplot(table.corr, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(cluster), ncol = 2)

ggsave(plot = plot, filename = "corr_heat_B-ES_clusters_groups_facet.jpg", dpi = 300, width = 15, height = 20)


### Try network_plot()
# ?network_plot
# Output a network plot of a correlation data frame in which variables that are more highly correlated appear closer together
# and are joined by stronger paths. Paths are also colored by their sign (blue for positive and red for negative).
# The proximity of the points are determined using multidimensional clustering.

# Uses correlate() # ?correlate
#data.frame(correlate(na.omit(mydata), method = "spearman", diagonal = 1))
#na.omit(mydata[mydata,]) %>% correlate(method = "spearman") %>% network_plot(min_cor = 0.5, colours = c("skyblue1", "white", "indianred2"))
#na.omit(mydata[mydata$clusters == 4,c(1:16)]) %>% correlate(method = "spearman") %>% network_plot(min_cor = 0.25, colours = c("#3288bd", "white", "#d53e4f"))

### In a for loop, plot and save correlation network and correlation heatmap
# for(i in unique(mydata$clusters)) {
#
#     message(paste("Making plots for cluster #",i, sep = ""))
#
#     cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:16)]), method = "spearman"),2)
#     upper_tri <- get_upper_tri(cormat)
#     melted_cormat <- melt(upper_tri, na.rm = F)
#     colnames(melted_cormat) <- c("V1","V2","rho")
#     melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
#     groups <- colnames(mydata)[c(7:16)]
#     vars <- c("Biodiversity","NPP","FPOCex","e ratio","PSD","catches")
#     melted_cormat <- melted_cormat[melted_cormat$V1 %in% vars,]
#     melted_cormat <- melted_cormat[melted_cormat$V2 %in% groups,]
#     colnames(melted_cormat) <- c("ES","Group","rho")
#     max.val <- max(abs(melted_cormat$rho))
#     cor.heat <- ggplot(melted_cormat, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
#         scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
#                 midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
#         xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#             axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#         geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3)
#
#     cor.net <- na.omit(mydata[mydata$clusters == i,c(1:16)]) %>%
#                 correlate(method = "spearman") %>%
#                 network_plot(min_cor = 0.5, colours = c("#3288bd","white","#d53e4f"))
#     # Save plots
#     ggsave(plot = cor.heat, filename = paste("corr_heat_B-ES_cluster#",i,".pdf", sep = ""), dpi = 300, width = 6, height = 9)
#     ggsave(plot = cor.net, filename = paste("corr_net_B-ES_cluster#",i,".pdf", sep = ""), dpi = 300, width = 13, height = 13)
#
# } # eo for loop - c in clusters
#
# # ### And draw global scale plots
# # cormat <- round(cor(na.omit(mydata[,c(1:16)]), method = "spearman"),2)
# # upper_tri <- get_upper_tri(cormat)
# # melted_cormat <- melt(upper_tri, na.rm = F)
# # colnames(melted_cormat) <- c("V1","V2","rho")
# # melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
# # groups <- colnames(mydata)[c(7:16)]
# # vars <- c("Biodiversity","NPP","FPOCex","e ratio","PSD","catches")
# # melted_cormat <- melted_cormat[melted_cormat$V1 %in% vars,]
# # melted_cormat <- melted_cormat[melted_cormat$V2 %in% groups,]
# # colnames(melted_cormat) <- c("ES","Group","rho")
# # max.val <- max(abs(melted_cormat$rho))
# # cor.heat <- ggplot(melted_cormat, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
# #     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
# #             midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
# #     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
# #         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
# #     geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3)
# #
# # cor.net <- na.omit(mydata[,c(1:16)]) %>%
# #             correlate(method = "spearman") %>%
# #             network_plot(min_cor = 0.25, colours = c("#3288bd","white","#d53e4f"))
# # # Save plots
# # ggsave(plot = cor.heat, filename = paste("corr_heat_B-ES_global.pdf", sep = ""), dpi = 300, width = 6, height = 9)
# # ggsave(plot = cor.net, filename = paste("corr_net_B-ES_global.pdf", sep = ""), dpi = 300, width = 13, height = 13)
#
#
# ### Same as above but only with Phyto- and Zooplankton SR
# data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
# names <- colnames(data.revi)[c(9:16)] ; names
# mydata <- data.revi[,names]
# # Change some colnames etc.
# colnames(mydata)[1] <- "Oceanic diversity"
# colnames(mydata)[4] <- "e ratio"
# colnames(mydata)[5] <- "PSD"
# colnames(mydata)[c(7:8)] <- c("Phytoplankton","Zooplankton")
#
# # Transform the columns needed
# mydata$FPOCex <- log(mydata$FPOCex)
# mydata$NPP <- log(mydata$NPP)
# mydata$catches <- log1p(mydata$catches)
# mydata[,c(7:length(mydata))] <- log1p(mydata[,c(7:length(mydata))])
#
# ### Plot global correlations coeff
# cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
# upper_tri <- get_upper_tri(cormat)
# melted_cormat <- melt(upper_tri, na.rm = F)
# colnames(melted_cormat) <- c("V1","V2","rho")
# melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
# groups <- c("Phytoplankton","Zooplankton")
# vars <- c("Oceanic diversity","NPP","FPOCex","e ratio","PSD","catches")
# melted_cormat <- melted_cormat[melted_cormat$V1 %in% vars,]
# melted_cormat <- melted_cormat[melted_cormat$V2 %in% groups,]
# colnames(melted_cormat) <- c("ES","Group","rho")
# # melted_cormat
# max.val <- max(abs(melted_cormat$rho))
# cor.heat <- ggplot(melted_cormat, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
#             midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3)
#
# cor.net <- na.omit(mydata) %>%
#             correlate(method = "spearman") %>%
#             network_plot(min_cor = 0.25, colours = c("#3288bd","white","#d53e4f"))
# # Save plots
# ggsave(plot = cor.heat, filename = paste("corr_heat_B-ES_global_nogroups.jpg", sep = ""), dpi = 300, width = 6, height = 9)
# ggsave(plot = cor.net, filename = paste("corr_net_B-ES_global_nogroups.jpg", sep = ""), dpi = 300, width = 8, height = 8)
#
# ### And now per clusters of severity
# mydata$clusters <- data.revi$clusters
# for(i in unique(mydata$clusters)) {
#
#     message(paste("Making plots for cluster #",i, sep = ""))
#
#     cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:8)]), method = "spearman"),2)
#     upper_tri <- get_upper_tri(cormat)
#     melted_cormat <- melt(upper_tri, na.rm = F)
#     colnames(melted_cormat) <- c("V1","V2","rho")
#     melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
#     groups <- c("Phytoplankton","Zooplankton")
#     vars <- c("Oceanic diversity","NPP","FPOCex","e ratio","PSD","catches")
#     melted_cormat <- melted_cormat[melted_cormat$V1 %in% vars,]
#     melted_cormat <- melted_cormat[melted_cormat$V2 %in% groups,]
#     colnames(melted_cormat) <- c("ES","Group","rho")
#     max.val <- max(abs(melted_cormat$rho))
#     cor.heat <- ggplot(melted_cormat, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
#         scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
#                 midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
#         xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#             axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#         geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3)
#
#     cor.net <- na.omit(mydata[mydata$clusters == i,c(1:8)]) %>%
#                 correlate(method = "spearman") %>%
#                 network_plot(min_cor = 0.25, colours = c("#3288bd","white","#d53e4f"))
#     # Save plots
#     ggsave(plot = cor.heat, filename = paste("corr_heat_B-ES_cluster#",i,"_nogroups.pdf", sep = ""), dpi = 300, width = 6, height = 9)
#     ggsave(plot = cor.net, filename = paste("corr_net_B-ES_cluster#",i,"_nogroups.pdf", sep = ""), dpi = 300, width = 8, height = 8)
#
# } # eo for loop - c in clusters

### Same as above but facet plot per clusters ! (use lapply before to retunr all cormats)
res <- lapply(unique(mydata$clusters), function(i) {
    
            message(paste("Computing corr matrox for cluster  ",i, sep = ""))
            
            cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:8)]), method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = F)
            colnames(melted_cormat) <- c("V1","V2","rho")
            melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
            groups <- c("Phytoplankton","Zooplankton")
            vars <- c("Oceanic diversity","NPP","FPOCex","e ratio","PSD","catches")
            melted_cormat <- melted_cormat[melted_cormat$V1 %in% vars,]
            melted_cormat <- melted_cormat[melted_cormat$V2 %in% groups,]
            colnames(melted_cormat) <- c("ES","Group","rho")
            # max.val <- max(abs(melted_cormat$rho))
            melted_cormat$cluster <- i
            return(melted_cormat)
    
    } # eo lapply 
    
) # eo lapply - i in clusters
# Rbind
table.corr <- dplyr::bind_rows(res)
#dim(table.corr) ; head(table.corr)
max.val <- max(abs(table.corr$rho)) ; max.val

plot <- ggplot(table.corr, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(cluster), ncol = 2)

ggsave(plot = plot, filename = "corr_heat_B-ES_clusters_nogroups_facet.pdf", dpi = 300, width = 8, height = 8)
 
 
### Same with PFTs:
data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
names <- colnames(data.revi)[c(9:14,19:28)] ; names
mydata <- data.revi[,names]
# Change some colnames etc.
colnames(mydata)[1] <- "Biodiversity"
colnames(mydata)[4] <- "e ratio"
colnames(mydata)[5] <- "PSD"
colnames(mydata)[c(7:16)] <- str_replace_all(colnames(mydata)[c(7:16)] , "_base", "")
colnames(mydata)[14] <- "Cnidaria"
# Transform the columns needed
mydata$FPOCex <- log(mydata$FPOCex)
mydata$NPP <- log(mydata$NPP)
mydata$catches <- log1p(mydata$catches)
mydata[,c(7:length(mydata))] <- log1p(mydata[,c(7:length(mydata))])
mydata$clusters <- data.revi$clusters  

res <- lapply(unique(mydata$clusters), function(i) {
    
            message(paste("Computing corr matrix for cluster  ",i, sep = ""))
            
            cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = F)
            colnames(melted_cormat) <- c("V1","V2","rho")
            melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
            groups <- colnames(mydata)[c(7:16)]
            vars <- c("Biodiversity","NPP","FPOCex","e ratio","PSD","catches")
            melted_cormat <- melted_cormat[melted_cormat$V1 %in% vars,]
            melted_cormat <- melted_cormat[melted_cormat$V2 %in% groups,]
            colnames(melted_cormat) <- c("ES","Group","rho")
            melted_cormat$cluster <- i
            
            return(melted_cormat)
    
    } # eo lapply 
    
) # eo lapply - i in clusters
# Rbind
table.corr <- dplyr::bind_rows(res)
#dim(table.corr) ; head(table.corr)
max.val <- max(abs(table.corr$rho)) ; max.val

plot <- ggplot(table.corr, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(cluster), ncol = 3)

ggsave(plot = plot, filename = "corr_heat_B-ES_clusters_groups_facet.pdf", dpi = 300, width = 10, height = 10)


# ------------------------------------------------------------

### 10/02/21: Plot correlations of total phytoplankton and the 3 phytoplankton PFTs and nutrients variables: NO3, SiO2, 
get_lower_tri <- function(cormat) {
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
}
get_upper_tri <- function(cormat) {
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
}
reorder_cormat <- function(cormat) {
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <- cormat[hc$order, hc$order]
        return(cormat)
}

data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t") ; colnames(data.revi)
# Change some colnames etc.
colnames(data.revi)[c(14:28)] <- str_replace_all(colnames(data.revi)[c(14:28)] , "_base", "")
colnames(data.revi)[c(15)] <- "Phytoplankton"
colnames(data.revi)[c(16)] <- "Zooplankton"
colnames(data.revi)[c(19)] <- "Diatoms"
colnames(data.revi)[c(23)] <- "Dinoflagellates"
colnames(data.revi)[c(25)] <- "Haptophytes"

### Combine with annual climatologies of nutrients fields
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d")
files <- dir()[grep("_18_09_20",dir())] ; files
clims <- lapply(files, function(f) {
            d <- read.table(f, h = T, sep = ";")
            # Rotate longitudes and add new ID
            d$x2 <- d$x
            d[d$x < 0 ,"x2"] <- (d[d$x < 0 ,"x"]) + 360
            d$id2 <- paste(d$x2, d$y, sep = "_")
            d <- d[,c("id2","x2","y","logNO3","logSiO2","Nstar","Sistar","Pstar")]
            return(d)
    } # eo 
)
clims <- bind_rows(clims)
#head(clims) ; summary(clims)
setwd(WD)
### Derive mean annual climatologies
ann.clim <- data.frame(clims %>% group_by(id2) %>%
                summarize(x = unique(x2), y = unique(y),
                NO3 = mean(logNO3, na.rm = T), SiO2 = mean(logSiO2, na.rm = T),
                Nstar = mean(Nstar, na.rm = T), Sistar = mean(Sistar, na.rm = T),
                Pstar = mean(Pstar, na.rm = T))
) # eo ddf
summary(ann.clim)

# ggplot(aes(x = Nstar, y = Pstar), data = ann.clim) + geom_point() + theme_classic() +

# OK, re-order, find common cells and cbind
ann.clim <- ann.clim[order(ann.clim$id),]
data.revi <- data.revi[order(data.revi$id),]
commons <- intersect(ann.clim$id, data.revi$id)
data.revi <- data.revi[data.revi$id %in% commons,]
ann.clim <- ann.clim[ann.clim$id2 %in% commons,]
data.revi[,c("NO3","SiO2","N*","Si*","P*")] <- ann.clim[,c("NO3","SiO2","Nstar","Sistar","Pstar")]
### Add absolute Si* and N*
data.revi[,"Si* (abs)"] <- abs(data.revi[,"Si*"])
data.revi[,"N* (abs)"] <- abs(data.revi[,"N*"])
data.revi[,"P* (abs)"] <- abs(data.revi[,"P*"])
summary(data.revi) ; dim(data.revi)

### Global scale corr heatmap
#colnames(data.revi)
cormat <- round(cor(na.omit(data.revi[,c("Phytoplankton","Diatoms","Dinoflagellates","Haptophytes",
                                        "NO3","SiO2","N*","Si*","P*","Si* (abs)","N* (abs)","P* (abs)")]),
                                        method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = T)
colnames(melted_cormat) <- c("V1","V2","rho")
melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
colnames(melted_cormat) <- c("Group1","Group2","rho")
melted_cormat <- melted_cormat[!(melted_cormat$Group1 == melted_cormat$Group2),]
melted_cormat <- melted_cormat[melted_cormat$Group1 %in% c("Phytoplankton","Diatoms","Dinoflagellates","Haptophytes"),]

p <- ggplot(melted_cormat, aes(factor(Group1), factor(Group2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(-0.9,0.9), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group1), factor(Group2), label = rho), color = "black", size = 3)

setwd(WD)
ggsave(plot = p, filename = "heatmap_corr_ann_rich_base_phytoxnutrients.jpg", dpi = 300, width = 4, height = 6)

### Not super obvious...per latitudinal bands?
corrs.bands <- lapply(unique(data.revi$band), function(b) {
    
                    cormat <- round(cor(na.omit(data.revi[data.revi$band == b,c("Phytoplankton","Diatoms","Dinoflagellates",
                                "Haptophytes","NO3","SiO2","N*","Si*","Si* (abs)","N* (abs)")]), method = "spearman"),2)
                    upper_tri <- get_upper_tri(cormat)
                    melted_cormat <- melt(upper_tri, na.rm = T)
                    colnames(melted_cormat) <- c("V1","V2","rho")
                    melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
                    colnames(melted_cormat) <- c("Group1","Group2","rho")
                    melted_cormat <- melted_cormat[!(melted_cormat$Group1 == melted_cormat$Group2),]
                    melted_cormat <- melted_cormat[melted_cormat$Group1 %in% c("Phytoplankton","Diatoms","Dinoflagellates","Haptophytes"),]
                    melted_cormat$band <- b
                    return(melted_cormat)
                       
    } # eo 
)
corrs.bands <- bind_rows(corrs.bands)
corrs.bands
max(corrs.bands$rho) # 0.88
p <- ggplot(corrs.bands, aes(factor(Group1), factor(Group2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(-0.9,0.9), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group1), factor(Group2), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(band))

setwd(WD)
ggsave(plot = p, filename = "heatmap_corr_ann_rich_base_phytoxnutrients_bands.jpg", dpi = 300, width = 7, height = 6)


### Check global relationships between Si* and diversity
m.data.revi <- melt(data.revi[,c("Phytoplankton","Diatoms","Dinoflagellates","Haptophytes","Si*")], id.vars = c("Si*"))
#head(m.data.revi)
colnames(m.data.revi)[1] <- "Si"

p1 <- ggplot(data = m.data.revi, aes(y = Si, x = value)) +
    geom_point(colour = "grey55", alpha = 0.5) +
    #geom_smooth(colour = "black", method = "lm", se=TRUE, formula = y ~ poly(x,3) ) + 
    geom_smooth(colour = "black", method = "lm", se=TRUE) + 
    ylab('Si* (µM)') + xlab('Species richness') + theme_classic() +
    facet_wrap(~factor(variable), scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed")
    
# and with absolute Si*+ 3rd order polynms
p2 <- ggplot(data = m.data.revi, aes(y = abs(Si), x = (value))) +
    geom_point(colour = "grey55", alpha = 0.5) +
    # geom_smooth(colour = "black", method = "lm", se=TRUE, formula = y ~ poly(x,2) ) + 
    geom_smooth(colour = "black", method = "gam") + 
    ylab('Absolute Si* (µM)') + xlab('Species richness') + theme_classic() +
    facet_wrap(~factor(variable), scales = "free") 
    
### Save the 2 plots above: 
setwd(WD)
ggsave(plot = p1, filename = "plot_ann_rich_base_phytogroupsxSi*_facet_lmfit.jpg", dpi = 300, width = 6, height = 6)
ggsave(plot = p2, filename = "plot_ann_rich_base_phytogroupsxabsSi*_facet.jpg", dpi = 300, width = 6, height = 6)


### Evaluate correlation between ABSOLUTE excess (or deficiency) of Silicates relative to NO3
cor(data.revi[,"Phytoplankton"], abs(data.revi[,"Si*"]), method = "spearman")
# -0.481 for looged SR ; same with unlogged

summary(lm(abs(get('Si*'))~ log(Phytoplankton), data = data.revi))
# Residual standard error: 8.603 on 35021 degrees of freedom
#Multiple R-squared:  0.203,	Adjusted R-squared:  0.2029 
#F-statistic:  8918 on 1 and 35021 DF,  p-value: < 2.2e-16
### Polynomial: 
summary(lm(abs(get('Si*')) ~ log(Phytoplankton) + I(log(Phytoplankton)^2) + I(log(Phytoplankton)^3), data = data.revi))
# or:
summary(lm( abs(get('Si*')) ~ poly(log(Phytoplankton),3), data = data.revi)) # signif decrease
summary(lm( abs(get('Si*')) ~ poly(log(Diatoms),3), data = data.revi)) # signif decrease
summary(lm( abs(get('Si*')) ~ poly(log(Dinoflagellates),3), data = data.revi)) # signif decrease
summary(lm( abs(get('Si*')) ~ poly(log(Haptophytes),3), data = data.revi)) # signif decrease

cor(data.revi[,"Diatoms"], abs(data.revi[,"Si*"]), method = "spearman")
# -0.257 for logged ; same
cor(data.revi[,"Dinoflagellates"], abs(data.revi[,"Si*"]), method = "spearman")
# -0.579 ; same
cor(data.revi[,"Haptophytes"], abs(data.revi[,"Si*"]), method = "spearman")
# -0.539 ; same

summary(lm(abs(get('Si*'))~ (Diatoms), data = data.revi))
# ns
summary(lm(abs(get('Si*'))~ log(Dinoflagellates), data = data.revi))
# Multiple R-squared:  0.3202,	Adjusted R-squared:  0.3202 
# F-statistic: 1.649e+04 on 1 and 35021 DF,  p-value: < 2.2e-16
summary(lm(abs(get('Si*'))~ log(Haptophytes), data = data.revi))
# Multiple R-squared:  0.3214,	Adjusted R-squared:  0.3214 
# F-statistic: 1.659e+04 on 1 and 35021 DF,  p-value: < 2.2e-16



### Check global relationships between N* and diversity
m.data.revi <- melt(data.revi[,c("Phytoplankton","Diatoms","Dinoflagellates","Haptophytes","N*")], id.vars = c("N*"))
#head(m.data.revi)

p1 <- ggplot(data = m.data.revi, aes(y = get("N*"), x = value)) +
    geom_point(colour = "grey55", alpha = 0.25) +
    geom_smooth(colour = "black", method = "lm", se = T) + 
    ylab('N* (µM)') + xlab('Species richness') + theme_classic() +
    facet_wrap(~factor(variable), scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed")

# and with absolute N*
p2 <- ggplot(data = m.data.revi, aes(y = abs(get('N*')), x = value)) +
    geom_point(colour = "grey55", alpha = 0.5) +
    #geom_smooth(colour = "black", method = "lm") + 
    ylab('Absolute N* (µM)') + xlab('Species richness') + theme_classic() +
    facet_wrap(~factor(variable), scales = "free") 

### Save the 2 plots above: 
setwd(WD)
ggsave(plot = p1, filename = "plot_ann_rich_base_phytogroupsxN*_facet_lmfit.jpg", dpi = 300, width = 6, height = 6)
ggsave(plot = p2, filename = "plot_ann_rich_base_phytogroupsxabsN*_facet.jpg", dpi = 300, width = 6, height = 6)


### Check correlations?
cor(data.revi[,"Phytoplankton"], data.revi[,"N*"], method = "spearman")
# -0.1686 ; +0.167 for unlooged and not absolute
cor((data.revi[,"Diatoms"]), data.revi[,"N*"], method = "spearman")
# -0.2453 ; +0.25 for unlooged and not absolute
cor((data.revi[,"Dinoflagellates"]), data.revi[,"N*"], method = "spearman")
# -0.0769 ; 
cor((data.revi[,"Haptophytes"]), data.revi[,"N*"], method = "spearman")
# -0.179 ; 

summary(lm(get('N*') ~ Phytoplankton, data = data.revi))
summary(lm(get('N*') ~ Diatoms, data = data.revi))
# Sign decrease but vety weak R2

### And quick panel plot of geom_points: groups x N*/Si*
# m.data.revi <- melt(data.revi[,c("band","Phytoplankton","Diatoms","Dinoflagellates","Haptophytes","Si*")], id.vars = c("band","Si*"))
# head(m.data.revi)
#
# ggplot(data = m.data.revi, aes(y = abs(get('Si*')), x = log(value))) +
#     geom_point(colour = "grey55", alpha = 0.5) + geom_smooth(colour = "black") +
#     ylab('Absolute Si* (µM)') + xlab('Species richness (log)') + theme_classic() +
#     facet_wrap(factor(variable) ~ factor(band), scales = "free") +
#     geom_hline(yintercept = 0, linetype = "dashed")
#
# ### Same but with
# m.data.revi <- melt(data.revi[,c("band","Phytoplankton","Diatoms","Dinoflagellates","Haptophytes","N*")], id.vars = c("band","N*"))
# head(m.data.revi)
#
# ggplot(data = m.data.revi, aes(y = abs(get('N*')), x = log(value))) +
#     geom_point(colour = "grey55", alpha = 0.5) + geom_smooth(colour = "black") +
#     ylab('Absolute N* (µM)') + xlab('Species richness (log)') + theme_classic() +
#     facet_wrap(factor(variable) ~ factor(band), scales = "free") +
#     geom_hline(yintercept = 0, linetype = "dashed")


### Check global relationships between N* and diversity
m.data.revi <- melt(data.revi[,c("Phytoplankton","Diatoms","Dinoflagellates","Haptophytes","P*")], id.vars = c("P*"))
#head(m.data.revi)

p1 <- ggplot(data = m.data.revi, aes(y = get("P*"), x = value)) +
    geom_point(colour = "grey55", alpha = 0.25) +
    geom_smooth(colour = "black", method = "lm", se = T) + 
    ylab('P* (µM)') + xlab('Species richness') + theme_classic() +
    facet_wrap(~factor(variable), scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed")

setwd(WD)
ggsave(plot = p1, filename = "plot_ann_rich_base_phytogroupsxP*_facet_lmfit.jpg", dpi = 300, width = 6, height = 6)

cor(data.revi[,"Phytoplankton"], abs(data.revi[,"P*"]), method = "spearman")
# -0.17
cor(data.revi[,"Diatoms"], abs(data.revi[,"P*"]), method = "spearman")
# -0.25
cor(data.revi[,"Dinoflagellates"], abs(data.revi[,"P*"]), method = "spearman")
# ns
cor(data.revi[,"Haptophytes"], abs(data.revi[,"P*"]), method = "spearman")
# -0.18


ggplot(data = m.data.revi, aes(y = get("P*"), x = value)) +
    geom_point(colour = "grey55", alpha = 0.25) +
    geom_smooth(colour = "black", method = "lm", se=TRUE, formula = y ~ poly(x,3) ) + 
    ylab('P* (µM)') + xlab('Species richness') + theme_classic() +
    facet_wrap(~factor(variable), scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed")

# ------------------------------------------------------------

### 08/02/21: Comments from Meike Vogt: look at plankton groups' pairwise correlations to say stuff about...trophic interactions?
data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
names <- colnames(data.revi)[c(9:16,19,23,25,22,27,26,21,20,28,24)] ; names
mydata <- data.revi[,names]
# Change some colnames etc.
colnames(mydata)[c(7:18)] <- str_replace_all(colnames(mydata)[c(7:18)] , "_base", "")
colnames(mydata)[c(7:10,12,15:17)] <- c("Phytoplankton","Zooplankton","Diatoms","Dinoflagellates","Copepods","Chordates","Chaetognaths","Pteropods")

# Transform the columns needed
mydata$clusters <- data.revi$clusters  

### Global scale corr heatmap
cormat <- round(cor(na.omit(mydata[,c(7:18)]), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = T)
colnames(melted_cormat) <- c("V1","V2","rho")
melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
colnames(melted_cormat) <- c("Group1","Group2","rho")
max.val <- max(abs(melted_cormat$rho)) ; max.val

p <- ggplot(melted_cormat, aes(factor(Group1), factor(Group2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group1), factor(Group2), label = rho), color = "black", size = 3)

ggsave(plot = p, filename = "corr_heat_groupsxgroups_ann_base_rich_global_facet.jpg", dpi = 300, width = 5.5, height = 5.5)


### And now the same but per clusters/ regions
# i <- 2
res <- lapply(unique(mydata$clusters), function(i) {
    
            message(paste("Computing corr matrix for cluster  ",i, sep = ""))
            
            cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(7:18)]), method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = T)
            colnames(melted_cormat) <- c("V1","V2","rho")
            melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
            colnames(melted_cormat) <- c("Group1","Group2","rho")
            melted_cormat$cluster <- i
            
            return(melted_cormat)
    
    } # eo lapply 
    
) # eo lapply - i in clusters
# Rbind
table.corr <- dplyr::bind_rows(res)
max.val <- max(abs(table.corr$rho)) ; max.val
plot <- ggplot(table.corr, aes(factor(Group1), factor(Group2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group1), factor(Group2), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(cluster), ncol = 3)

ggsave(plot = plot, filename = "corr_heat_groupsxgroups_ann_rich_clusters_facet.jpg", dpi = 300, width = 12, height = 10)


### Same, but based on ∆SR of the PFTs. Global and then per regions
names <- colnames(data.revi)[c(17,18,29,33,35,32,37,36,31,30,38,34)] ; names
mydata <- data.revi[,names]
# Change some colnames etc.
colnames(mydata) <- str_replace_all(colnames(mydata),"_fut","")
colnames(mydata)[c(1,2,3,4,6,9:11)] <- c("Phytoplankton","Zooplankton","Diatoms","Dinoflagellates","Copepods","Chordates","Chaetognaths","Pteropods")
mydata$clusters <- data.revi$clusters  
 
### Global scale corr heatmap
cormat <- round(cor(na.omit(mydata[,c(1:12)]), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = T)
colnames(melted_cormat) <- c("V1","V2","rho")
melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
colnames(melted_cormat) <- c("Group1","Group2","rho")
max.val <- max(abs(melted_cormat$rho)) ; max.val
p <- ggplot(melted_cormat, aes(factor(Group1), factor(Group2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group1), factor(Group2), label = rho), color = "black", size = 3)
#
ggsave(plot = p, filename = "corr_heat_groupsxgroups_ann_diff_rich_global.jpg", dpi = 300, width = 5.5, height = 5.5)


### And now the same but per clusters/ regions
# i <- 4
res <- lapply(unique(mydata$clusters), function(i) {
    
            message(paste("Computing corr matrix for cluster  ",i, sep = ""))
            
            cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:12)]), method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = T)
            colnames(melted_cormat) <- c("V1","V2","rho")
            melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
            colnames(melted_cormat) <- c("Group1","Group2","rho")
            melted_cormat$cluster <- i
            
            return(melted_cormat)
    
    } # eo lapply 
    
) # eo lapply - i in clusters
# Rbind
table.corr <- dplyr::bind_rows(res)
#table.corr
max.val <- max(abs(table.corr$rho)) ; max.val
plot <- ggplot(table.corr, aes(factor(Group1), factor(Group2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(Group1), factor(Group2), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(cluster), ncol = 3)

ggsave(plot = plot, filename = "corr_heat_groupsxgroups_ann_diff_rich_clusters_facet.jp", dpi = 300, width = 12, height = 10)


# --------------------------------------------------------------------------------------------------------------------------------

### 08/02/21: Add panel (10 plots) of zonal LDG per PFT to be added to the Appendices (Fig. S1) in order to strengthen SDM predictions

data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
names <- colnames(data.revi)[c(19:28)] ; names
colnames(data.revi)[c(19:28)] <- str_replace_all(names,"_base","")
colnames(data.revi)[26] <- "Cnidaria"
colnames(data.revi)[c(19:28)] <- c("Diatoms","Chaetognaths","Chordates","Copepods","Dinoflagellates","Foraminifera",
                                    "Haptophytes","Jellyfish","Krill","Pteropods")
groups <- colnames(data.revi)[c(19:28)] ; groups
# g <- "Copepods"

# myplots <- list()  # new empty list

# for(g in groups) {
#
#     max <- max(data.revi[,g])
#     zonal <- data.frame(data.revi %>% group_by(y) %>% summarize(avg_sr = mean(get(g)/max, na.rm = T), sd = sd(get(g)/max, na.rm = T)))
#
#     z.plot <- ggplot() + geom_ribbon(aes(y = y, xmin = avg_sr-sd, xmax = avg_sr+sd), fill = "grey75", data = zonal) +
#              geom_path(aes(y = y, x = avg_sr), data = zonal, colour = "black", linetype = "solid") +
#              scale_x_continuous(name = paste(g," richness", sep = "")) +
#              scale_y_continuous(position = "right", name = "", limits = c(-90,90),
#                  breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
#              theme_classic()
#
#     myplots[[g]] <- z.plot  # add each plot into plot list
# } # eo for loop - g in groups
# #multiplot(plotlist = myplots, cols = 4)
#
# require("ggpubr")
# panel <- ggarrange(myplots[["Diatoms"]],myplots[["Dinoflagellates"]],myplots[["Haptophytes"]],
#     myplots[["Copepods"]],myplots[["Krill"]],myplots[["Jellyfish"]],
#     myplots[["Pteropods"]],myplots[["Chaetognaths"]],myplots[["Chordates"]],
#     myplots[["Foraminifera"]], labels = LETTERS[1:10], align = "hv", ncol = 3, nrow = 4)
# 
# ggsave(plot = panel, filename = "panel_zonal_plots_ann_rich_base_ens_groups.jpg", dpi = 300, width = 10, height = 15)

### Same but with maps of normalized richness
# myplots <- list()  # new empty list

for(g in groups) {
   
    message(paste("Mapping ",g," baseline richness", sep = ""))
   
    subset <- data.revi[,c("x","y",g)]
    # max <- max(subset[,g]) 
    ### For plotting %SR instead of SR, divide by number of modelled species per PFG
    if(g == "Diatoms") {
        max <- 154
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Dinoflagellates") {
        max <- 154
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Haptophytes") {
        max <- 24
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Copepods") {
        max <- 272
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Krill") {
        max <- 73
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Chaetognaths") {
        max <- 27
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Jellyfish") {
        max <- 69
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Chordates") {
        max <- 17
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Foraminifera") {
        max <- 26
        subset[,"perc"] <- subset[,g]/ max
    } else if(g == "Pteropods") {
        max <- 22
        subset[,"perc"] <- subset[,g]/ max
    } 
   
    zonal <- data.frame(subset %>% group_by(y) %>% summarize(avg_perc = mean(perc, na.rm = T), sd = sd(perc, na.rm = T)) )  
    
    # Map and Zonal plot
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = subset) +
        geom_contour(colour = "grey60", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = perc), data = subset) +
        scale_fill_viridis(name = paste("% Species richness\n","(",max," species)", sep = ""),) +
        geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
        coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
        scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
        scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

    z.plot <- ggplot() + geom_ribbon(aes(y = y, xmin = avg_perc-sd, xmax = avg_perc+sd), fill = "grey75", data = zonal) +
         geom_path(aes(y = y, x = avg_perc), data = zonal, colour = "black", linetype = "solid") +
         scale_x_continuous(name = "") + scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
             breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
         theme_classic()

    # Assemble in panel
    require("ggpubr")
    require("grDevices")
    panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")   
    
    setwd(WD)
    ggsave(plot = panel, filename = paste("map+zonal_ann_perc_base_ens_",g,".pdf", sep = ""), dpi = 300, width = 8.25, height = 3)    
           
} # eo for loop - g in groups



# --------------------------------------------------------------------------------------------------------------------------------

### 3°) 30/11/2020: Use copepod maximum body lengths from Razouls et al. 

### First, get the dataset of copepod species body length
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology")
sizes <- read.csv("copepod_sizes_razouls.csv", h = T, sep = ";", dec = ",")
dim(sizes) ; str(sizes) ; summary(sizes)
# Remove brackets in species names
sizes$sp.name <- gsub("\\(|\\)", "", sizes$sp.name)
specieswithsizes <- unique(sizes[sizes$n >= 3,"sp.name"])
specieswithsizes

### Then, fetch mean annual baseline (and then future) copepod community composition
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
files.base <- dir()[grep("ann_compo_zoo_baseline",dir())]
files.fut <- dir()[grep("ann_compo_zoo_2100-2000",dir())] # for later

### For each files.base, retrieve copepod community composition and derive estiate of annual size structure 
# f <- files.base[8] # for testing function below
require("parallel")
res <- mclapply(files.base, function(f) {
            
            # Message
            message(paste(f, sep = ""))
            comm <- get(load(f)) # dim(comm) ; colnames(comm)
            # Change colnames (remove brackets or points)
            colnames(comm) <-  gsub("\\.|\\.", "", colnames(comm))
            # Colnames of copeopods should match the 'specieswithsizes' vector
            common.spp <- intersect(specieswithsizes,colnames(comm))
            # length(common.spp) # 249 species, 47.5% of zooplankton species; 29% of ALL species modelled (860)
            subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
            rm(comm) ; gc()
            # Melt, to put species names as vector
            m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
            colnames(m.sub)[c(4,5)] <- c("species","HSI")
            m.sub$size <- NA
            for(sp in unique(m.sub$species)) {
                s <- sizes[sizes$sp.name == sp,"max"]
                m.sub[m.sub$species == sp,"size"] <- s
            } # eo for loop
            # Use summaize to derive estimate of HSI weighted median body length (size structure)
            # require("Hmisc")  !!! https://stackoverflow.com/questions/33807624/understanding-ddply-error-message 
            require("matrixStats") ; require("dplyr")
            # weightedMedian(x = m.sub$size, w = m.sub$HSI)
            med.size <- data.frame(m.sub %>% group_by(cell_id) %>%
                summarize(x = unique(x), y = unique(y),
                med.size = weightedMedian(x = size, w = HSI), 
                var.size = weightedVar(x = size, w = HSI)
            )) 
            
            # head(med.size) ; summary(med.size)
            # Quick map to check
             # ggplot() + geom_raster(aes(x = x, y = y, fill = var.size), data = med.size) +
#                  geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = var.size), data = med.size) +
#                  scale_fill_viridis(name = "Copepod\nsize variance\n(mm)") +
#                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#                  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#                  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#                  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

            return(med.size)
    
        }, mc.cores = length(files.base)
        
) # eo mclapply - f in files.base
# Rbind
table.med.size <- bind_rows(res) ; rm(res) ; gc()
# head(table.med.size) ; dim(table.med.size) ; summary(table.med.size)
# Derive ensemble estimate of median size structure
med.size <- data.frame(table.med.size %>% group_by(cell_id) %>%
    summarize(x = unique(x), y = unique(y),
        med.size = median(med.size, na.rm = T),
        med.size.var = median(var.size, na.rm = T)
))

summary(med.size)

map.size.base <- ggplot() + geom_raster(aes(x = x, y = y, fill = med.size), data = med.size) +
    geom_contour(colour = "grey60", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = med.size), data = med.size) +
    scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map.var.size.base <- ggplot() + geom_raster(aes(x = x, y = y, fill = med.size.var), data = med.size) +
    geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = med.size.var), data = med.size) +
    scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
setwd(WD)
ggsave(plot = map.size.base, filename = "map_ens_ann_med_copepod_size_baseline.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map.var.size.base, filename = "map_ens_ann_med_variance_copepod_size_baseline.pdf", width = 7, height = 4, dpi = 300)


### Combine 'med.size' with table_data_for_revisions_26_11_20
setwd(WD)
data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
# head(data.revi)
commons <- intersect(data.revi$id, med.size$cell_id) ; length(commons)
data.revi <- data.revi[order(data.revi$id),]
med.size <- med.size[order(med.size$cell_id),]
data.revi <- data.revi[data.revi$id %in% commons,]
med.size <- med.size[med.size$cell_id %in% commons,]
data.revi$med.size <- med.size$med.size
data.revi$med.size.var <- med.size$med.size.var

# Re-number clusters of severity propoerly
data.revi$clusters <- data.revi$pca_euclid_pam_k6
data.revi[data.revi$pca_euclid_pam_k6 == 1,"clusters"] <- 4
data.revi[data.revi$pca_euclid_pam_k6 == 2,"clusters"] <- 3
data.revi[data.revi$pca_euclid_pam_k6 == 3,"clusters"] <- 6
data.revi[data.revi$pca_euclid_pam_k6 == 4,"clusters"] <- 2
data.revi[data.revi$pca_euclid_pam_k6 == 6,"clusters"] <- 1


### OK, couple of things to do from here:
# - plot distribution of med.size per clusters (boxplots)
# - plot correlation heatmap bewteen medi.size and the other variables (log SR and ES)
# - plot some bivariate plot for the strong and interesting corr found above (with POC, e ratio or PSD)

plot <- ggplot(data = data.revi, aes(x = factor(clusters), y = med.size, fill = factor(clusters))) +
   geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") + 
   ylab(paste("Copepod median body size (mm)", sep = ""))  + xlab("") + theme_classic()
 
ggsave(plot = plot, filename = "plot_ens_ann_med_copepod_size_clusters_facet.pdf", width = 6, height = 4.5, dpi = 300)
# Clear gradients: 5 > 1 > 2 > 6 > 3-4
plot <- ggplot(data = data.revi, aes(x = factor(clusters), y = med.size.var, fill = factor(clusters))) +
   geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") + 
   ylab(paste("Copepod median body size variance (mm)", sep = ""))  + xlab("") + theme_classic()
 
ggsave(plot = plot, filename = "plot_ens_ann_med.var_copepod_size_clusters_facet.pdf", width = 6, height = 4.5, dpi = 300)


# Corr coef heatmap to analyze trend between copepod size str indices and relevant variables: 
get_lower_tri <- function(cormat){
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
}
get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
}
reorder_cormat <- function(cormat){
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <- cormat[hc$order, hc$order]
        return(cormat)
}

### !! Add mean annual SiO2, NO3, N* and Si*
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d")
files <- dir()[grep("_18_09_20",dir())] ; files
clims <- lapply(files, function(f) {
            d <- read.table(f, h = T, sep = ";")
            # Rotate longitudes and add new ID
            d$x2 <- d$x
            d[d$x < 0 ,"x2"] <- (d[d$x < 0 ,"x"]) + 360
            d$id2 <- paste(d$x2, d$y, sep = "_")
            d <- d[,c("id2","x2","y","SST","logChl")]
            return(d)
    } # eo 
)
clims <- bind_rows(clims)

setwd(WD)
### Derive mean annual climatologies
ann.clim <- data.frame(clims %>% group_by(id2) %>%
                summarize(x = unique(x2), y = unique(y), SST = mean(SST, na.rm = T), logChl = mean(logChl, na.rm = T))
) # eo ddf
summary(ann.clim)

# OK, re-order, find common cells and cbind
ann.clim <- ann.clim[order(ann.clim$id),]
data.revi <- data.revi[order(data.revi$id),]
commons <- intersect(ann.clim$id, data.revi$id)
data.revi <- data.revi[data.revi$id %in% commons,]
ann.clim <- ann.clim[ann.clim$id2 %in% commons,]
data.revi[,c("SST","logChl")] <- ann.clim[,c("SST","logChl")]

### Plot heatmpa of global correlations between copepod SR, size strcuture and ecosystem functions that can be linked to copepod size
names <- colnames(data.revi)[c(10:13,43:44,22,41,42)] ; names
mydata <- data.revi[,names]
# Change some colnames etc.
#colnames(mydata)[1] <- "Biodiversity"
colnames(mydata)[3] <- "e ratio"
colnames(mydata)[4] <- "PSI"
colnames(mydata)[6] <- "Chlorophyll\n(log)"
colnames(mydata)[7] <- 'Copepod\nSR'
#colnames(mydata)[14] <- "Cnidaria"
colnames(mydata)[8] <- "Copepod median\nbody size"
colnames(mydata)[9] <- "Copepod median\nbody size variance"
# Transform the columns needed
mydata$FPOCex <- log(mydata$FPOCex)
mydata$NPP <- log(mydata$NPP)
# mydata$catches <- log1p(mydata$catches)
# mydata[,c(7:16)] <- log1p(mydata[,c(7:16)])

# Corr matrix
cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = T)
colnames(melted_cormat) <- c("V1","V2","rho")
melted_cormat
melted_cormat <- melted_cormat[melted_cormat$V2 %in% c("Copepod\nSR","Copepod median\nbody size","Copepod median\nbody size variance"),]
melted_cormat <- melted_cormat[!(melted_cormat$V2 == melted_cormat$V1),]
#melted_cormat <- melted_cormat[melted_cormat$V2 == "Copepod median size variance",]
colnames(melted_cormat) <- c("ES","Group","rho")
max.val <- max(abs(melted_cormat$rho)) ; max.val

heatmap <- ggplot(melted_cormat, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + coord_fixed() +
    geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3)

setwd(WD)
ggsave(plot = heatmap, filename = "heatmap_ens_ann_copepod_size_struct.vs.ES.pdf", width = 10, height = 4, dpi = 300)




### Examine bivariate plots
# With e ratio
summary(lm(med.size ~ e, data = data.revi))
ggplot(data = data.revi) + geom_point(aes(x = med.size, y = e), alpha = .5, colour = "grey50") +
    ylab("e ratio") + xlab("Copepod median size (mm)") + theme_classic() 
# With NPP?
summary(lm(med.size ~ log(NPP), data = data.revi))
ggplot(data = data.revi) + geom_point(aes(x = med.size, y = log(NPP)), alpha = .5, colour = "grey50") +
    ylab("NPP (log)") + xlab("Copepod median size (mm)") + theme_classic() 
# With FPOCex? Nope
# summary(lm(med.size ~ (FPOCex), data = data.revi))
# With slope? Nope
#summary(lm(med.size ~ slope2, data = data.revi))
# Catches
summary(lm(med.size ~ log1p(catches), data = data.revi))
ggplot(data = data.revi) + geom_point(aes(x = med.size, y = log1p(catches) ), alpha = .5, colour = "grey50") +
    ylab("Fish catches (log)") + xlab("Copepod median size (mm)") + theme_classic() 

# Meh, makes no sense honestly...analyze correlations cluster-wise
names <- colnames(data.revi)[c(9:14,19:28,40,41)] ; names
mydata <- data.revi[,names]
# Change some colnames etc.
colnames(mydata)[1] <- "Biodiversity"
colnames(mydata)[4] <- "e ratio"
colnames(mydata)[5] <- "PSD"
colnames(mydata)[c(7:16)] <- str_replace_all(colnames(mydata)[c(7:16)] , "_base", "")
colnames(mydata)[14] <- "Cnidaria"
colnames(mydata)[18] <- "Copepod size"
# Transform the columns needed
mydata$FPOCex <- log(mydata$FPOCex)
mydata$NPP <- log(mydata$NPP)
mydata$catches <- log1p(mydata$catches)
mydata[,c(7:16)] <- log1p(mydata[,c(7:16)])

res <- lapply(unique(mydata$clusters), function(i) {
    
            message(paste("Computing corr matrix for cluster  ",i, sep = ""))
            
            cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:16,18)]), method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = F)
            colnames(melted_cormat) <- c("V1","V2","rho")
            melted_cormat <- melted_cormat[melted_cormat$V2 == "Copepod size",]
            melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
            melted_cormat$cluster <- i
            
            return(melted_cormat)
    
    } # eo lapply 
    
) # eo lapply - i in clusters
# Rbind
table.corr <- dplyr::bind_rows(res)
max.val <- max(abs(table.corr$rho)) ; max.val

# Make sure 'V1' follows the order you want though
table.corr$V1 <- factor(table.corr$V1, levels = unique(table.corr$V1))

plot <- ggplot(table.corr, aes(factor(V1), factor(V2), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(V1), factor(V2), label = rho), color = "black", size = 3) +
    facet_wrap(.~factor(cluster), ncol = 2)

ggsave(plot = plot, filename = "corr_heat_B-ES_size_clusters_groups_facet.pdf", dpi = 300, width = 15, height = 7)
 
 
 
### OK, next step: evaluate changes in copepod med size per clusters between future and baseline
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
files.fut 

### For each files.base, retrieve copepod community composition and derive estiate of annual size structure 
f <- files.fut[1] # for testing function below
res <- mclapply(files.fut, function(f) {
            
            # Message
            message(paste(f, sep = ""))
            comm <- get(load(f)) # dim(comm) ; colnames(comm)
            # Change colnames (remove brackets or points)
            colnames(comm) <-  gsub("\\.|\\.", "", colnames(comm))
            # Colnames of copeopods should match the 'specieswithsizes' vector
            common.spp <- intersect(specieswithsizes,colnames(comm))
            # length(common.spp) # 249 species, 47.5% of zooplankton species; 29% of ALL species modelled (860)
            subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
            rm(comm) ; gc()
            # Melt, to put species names as vector
            m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
            colnames(m.sub)[c(4,5)] <- c("species","HSI")
            m.sub$size <- NA
            for(sp in unique(m.sub$species)) {
                s <- sizes[sizes$sp.name == sp,"max"]
                m.sub[m.sub$species == sp,"size"] <- s
            } # eo for loop
            require("matrixStats") ; require("dplyr")
            # weightedMedian(x = m.sub$size, w = m.sub$HSI)
            med.size <- data.frame(m.sub %>% group_by(cell_id) %>%
                    summarize(x = unique(x), y = unique(y),
                        med.size = weightedMedian(x = size, w = HSI),
                        var.size = weightedVar(x = size, w = HSI)
            ))
            # head(med.size) ; summary(med.size)
            # Return with ESM, SDM etc.
            terms <- do.call(cbind,strsplit(f, split = "_"))
            med.size$SDM <- terms[7,1]
            med.size$ESM <- terms[6,1]
            med.size$pool <- str_replace_all(terms[8,1],".Rdata","")
            
            return(med.size)
    
        }, mc.cores = length(files.base)
        
) # eo mclapply - f in files.base
# Rbind
table.med.size <- bind_rows(res) ; rm(res) ; gc()
head(table.med.size) ; dim(table.med.size) ; summary(table.med.size)

# Derive ensemble estimate of median size structure across all ESM and SDM
med.size.fut <- data.frame(table.med.size %>%
            group_by(cell_id) %>%
            summarize(x = unique(x), y = unique(y), 
            med.size = median(med.size, na.rm = T),
            med.size.var = median(var.size, na.rm = T) 
) )

summary(med.size.fut)

ggplot() + geom_raster(aes(x = x, y = y, fill = med.size.var), data = med.size.fut) +
    geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = med.size.var), data = med.size.fut) +
    scale_fill_viridis(name = "Copepod median\nsize variance (mm)") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### Derive % difference in copepod median size
# dim(med.size) ; dim(med.size.fut)
commons <- intersect(data.revi$id, med.size.fut$cell_id) ; length(commons)
data.revi <- data.revi[order(data.revi$id),]
med.size.fut <- med.size.fut[order(med.size.fut$cell_id),]
data.revi <- data.revi[data.revi$id %in% commons,]
med.size.fut <- med.size.fut[med.size.fut$cell_id %in% commons,]
data.revi$med.size.fut <- med.size.fut$med.size
data.revi$med.size.var.fut <- med.size.fut$med.size.var

### Derive % change
data.revi$diff.size <- (data.revi$med.size.fut)-(data.revi$med.size)
data.revi$perc.size <- (data.revi$diff.size/(data.revi$med.size))*100
data.revi$diff.size.var <- (data.revi$med.size.var.fut)-(data.revi$med.size.var)
data.revi$perc.size.var <- (data.revi$diff.size.var/(data.revi$med.size.var))*100

summary(data.revi)

cor(data.revi$perc.size, data.revi$perc.size.var, method = "spearman") # 0.490

map.perc.med.size <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc.size), data = data.revi) +
    geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = perc.size), data = data.revi) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map.perc.med.size.var <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc.size.var), data = data.revi) +
    geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = perc.size.var), data = data.revi) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### And plot distrbition % difference in copepod median 
plot1 <- ggplot(data = data.revi, aes(x = factor(clusters), y = perc.size, fill = factor(clusters))) +
   geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") + 
   ylab(paste("% Difference in copepod median body size (mm)", sep = "")) + xlab("") + theme_classic() + 
   geom_hline(yintercept = 0, linetype = "dashed") 
#
plot2 <- ggplot(data = data.revi, aes(x = factor(clusters), y = perc.size.var, fill = factor(clusters))) +
   geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") + 
   ylab(paste("% Difference in copepod median body size variance (mm)", sep = "")) + xlab("") + theme_classic() + 
   geom_hline(yintercept = 0, linetype = "dashed") 


# --> decrease in zooplankton size structure in clusters 1 & 2: linked to decrease in HSI + poleward migrations of warm water communities

setwd(WD)
ggsave(plot = map.perc.med.size, filename = "map_ens_perc_med_copepod_size_2100-2000.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map.perc.med.size.var, filename = "map_ens_perc_med.var_copepod_size_2100-2000.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = plot1, filename = "plot_ens_perc_ann_med_copepod_size_clusters.pdf", width = 6, height = 4.5, dpi = 300)
ggsave(plot = plot2, filename = "plot_ens_perc_ann_med.var_copepod_size_clusters.pdf", width = 6, height = 4.5, dpi = 300)


# --------------------------------------------------------------------------------------------------------------------------------


### 08/02/2021: Use phytoplankton mean cell diameters from Damiano

# ### First, get the dataset of phyto species cell diamaters from Damiano
# setwd("/net/kryo/work/fabioben/OVERSEE/data/biology")
# sizes <- read.csv("names_of_taxa_data_Phytoplankton_size_united_ESSD.csv", h = T, sep = ",")
# dim(sizes) ; str(sizes) ; summary(sizes)
# # ADD underscores to taxon name
# sizes$taxon <- gsub(" ", "_", sizes$taxon )
# length(unique(sizes$taxon)) # 1'716
# # Subset modelled species
# sizes <- sizes[sizes$modelled == "YES",]
# dim(sizes)
# # Column of interest --> idealized_diameter_in_micrometers_mean_of_sets
# # OR: biovolume_in_micrometers_3_weighted_mean_of_sets ?
#
# ### Then, fetch mean annual baseline (and then future) copepod community composition
# setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# files.base <- dir()[grep("ann_compo_phyto_baseline",dir())]
# files.fut <- dir()[grep("ann_compo_phyto_2100-2000",dir())] # for later
#
# ### For each files.base, retrieve phytoplankton community composition and estimate of annual cell dimater structure
# require("parallel")
#
# # f <- files.base[8] # for testing function below
# res <- mclapply(files.base, function(f) {
#
#             # Message
#             message(paste(f, sep = ""))
#             comm <- get(load(f)) # dim(comm) ; colnames(comm)
#             # Change colnames (remove brackets or points)
#             colnames(comm) <-  gsub("\\.|\\.", "", colnames(comm))
#             # Colnames of copeopods should match the 'specieswithsizes' vector
#             common.spp <- intersect(sizes$taxon,colnames(comm))
#             # species in the colnames of 'comm' that could have a estimate of cell diameter
#
#             ### Check % of species without NA in their diameter estimate
#             #sub.sizes <- sizes %>% drop_na(idealized_diameter_in_micrometers_mean_of_sets)
#             #dim(sub.sizes[sub.sizes$taxon %in% common.spp,]) # 239 ! --> 71ù of all phyto species!
#
#             subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
#             rm(comm) ; gc()
#             # Melt, to put species names as vector
#             m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
#             colnames(m.sub)[c(4,5)] <- c("species","HSI")
#             m.sub$diameter <- NA
#             for(sp in unique(m.sub$species)) {
#                 s <- sizes[sizes$taxon == sp,"idealized_diameter_in_micrometers_mean_of_sets"]
#                 m.sub[m.sub$species == sp,"diameter"] <- s
#             } # eo for loop
#             # Check
#             # summary(m.sub)
#             # Use summaize to derive estimate of HSI weighted median body length (size structure)
#             #require("Hmisc")  !!! https://stackoverflow.com/questions/33807624/understanding-ddply-error-message
#             require("matrixStats") ; require("dplyr")
#             # weightedMedian(x = m.sub$size, w = m.sub$HSI)
#             mean.size <- data.frame(m.sub %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y),
#                             mean.diameter = weightedMean(x = diameter, w = HSI, na.rm = T)) )
#             # head(med.size) ; summary(med.size)
#             # Quick map to check
#             # ggplot() + geom_raster(aes(x = x, y = y, fill = med.diameter), data = med.size) +
# #                  geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = med.diameter), data = med.size) +
# #                  scale_fill_viridis(name = "Phytoplankton\nmedian diameter\n(µm)") +
# #                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
# #                  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
# #                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
# #                  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
# #                  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
#             return(mean.size)
#
#         }, mc.cores = length(files.base)
#
# ) # eo mclapply - f in files.base
# # Rbind
# table.mean.size <- bind_rows(res) ; rm(res) ; gc()
# # head(table.mean.size) ; dim(table.mean.size) ; summary(table.mean.size)
# # Derive ensemble estimate of median size structure
# mean.size <- data.frame(table.mean.size %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), ens.diameter = mean(mean.diameter, na.rm = T)))
# # summary(mean.size)
# map.size.base <- ggplot() + geom_raster(aes(x = x, y = y, fill = ens.diameter), data = mean.size) +
#     geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = ens.diameter), data = mean.size) +
#     scale_fill_viridis(name = "Phytoplankton\nmean diameter\n(µm)") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
# #
# setwd(WD)
# ggsave(plot = map.size.base, filename = "map_ens_ann_mean_phytoplankton_cell_diameter_baseline.pdf", width = 7, height = 4, dpi = 300)
#
#
# ### Combine 'med.size' with table_data_for_revisions_26_11_20
# setwd(WD)
# data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
# # head(data.revi)
# commons <- intersect(data.revi$id, mean.size$cell_id) ; length(commons)
# data.revi <- data.revi[order(data.revi$id),]
# mean.size <- mean.size[order(mean.size$cell_id),]
# data.revi <- data.revi[data.revi$id %in% commons,]
# mean.size <- mean.size[mean.size$cell_id %in% commons,]
# data.revi$size.phyto <- mean.size$ens.diameter
#
# # Re-number clusters of severity propoerly
# data.revi$clusters <- data.revi$pca_euclid_pam_k6
# data.revi[data.revi$pca_euclid_pam_k6 == 1,"clusters"] <- 4
# data.revi[data.revi$pca_euclid_pam_k6 == 2,"clusters"] <- 3
# data.revi[data.revi$pca_euclid_pam_k6 == 3,"clusters"] <- 6
# data.revi[data.revi$pca_euclid_pam_k6 == 4,"clusters"] <- 2
# data.revi[data.revi$pca_euclid_pam_k6 == 6,"clusters"] <- 1
#
# ### OK, couple of things to do from here:
# # - plot distribution of size.phyto per clusters (boxplots)
# # - plot correlation heatmap bewteen size.phyto and the other variables (log SR and ES)
# # - plot some bivariate plot for the strong and interesting corr found above (with POC, e ratio or PSD)
#
# plot <- ggplot(data = data.revi, aes(x = factor(clusters), y = size.phyto, fill = factor(clusters))) +
#    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") +
#    ylab(paste("Phytoplankton cell diameter (µm)", sep = "")) + xlab("") + theme_classic()
#
# ggsave(plot = plot, filename = "plot_ens_ann_mean_phyto_cell_diamater_clusters_facet.pdf", width = 6, height = 6, dpi = 300)
#
#
# # Corr coef heatmap
# get_lower_tri <- function(cormat){
#       cormat[upper.tri(cormat)] <- NA
#       return(cormat)
# }
# get_upper_tri <- function(cormat){
#       cormat[lower.tri(cormat)]<- NA
#       return(cormat)
# }
# reorder_cormat <- function(cormat){
#         dd <- as.dist((1-cormat)/2)
#         hc <- hclust(dd)
#         cormat <- cormat[hc$order, hc$order]
#         return(cormat)
# }
#
# names <- colnames(data.revi)[c(9:14,19:28,41)] ; names
# mydata <- data.revi[,names]
# # Change some colnames etc.
# colnames(mydata)[1] <- "Biodiversity"
# colnames(mydata)[4] <- "e ratio"
# colnames(mydata)[5] <- "PSD"
# colnames(mydata)[c(7:16)] <- str_replace_all(colnames(mydata)[c(7:16)] , "_base", "")
# #colnames(mydata)[14] <- "Cnidaria"
# colnames(mydata)[17] <- "Phytoplankton\nmean diameter"
# # Transform the columns needed
# mydata$FPOCex <- log(mydata$FPOCex)
# mydata$NPP <- log(mydata$NPP)
# mydata$catches <- log1p(mydata$catches)
# mydata[,c(7:16)] <- log1p(mydata[,c(7:16)])
# # Corr matrix
# cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
# upper_tri <- get_upper_tri(cormat)
# melted_cormat <- melt(upper_tri, na.rm = T)
# colnames(melted_cormat) <- c("V1","V2","rho")
# melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
# melted_cormat <- melted_cormat[melted_cormat$V2 == "Phytoplankton\nmean diameter",]
# colnames(melted_cormat) <- c("ES","Group","rho")
# max.val <- max(abs(melted_cormat$rho)) ; max.val
#
# ggplot(melted_cormat, aes(factor(Group), factor(ES), fill = rho)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
#             midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
#     xlab("") + ylab("") + theme_minimal() + coord_fixed() +
#     geom_text(aes(factor(Group), factor(ES), label = rho), color = "black", size = 3)
#
# ### Examine bivariate plots
# # With e ratio
# summary(lm(size.phyto ~ e, data = data.revi))
# ggplot(data = data.revi) + geom_point(aes(x = size.phyto, y = e), alpha = .5, colour = "grey50") +
#     geom_smooth(aes(x = size.phyto, y = e), colour = "red") +
#     ylab("e ratio") + xlab("Phytoplankton median cell diameter (µm)") + theme_classic()
# # With NPP? NO
# #summary(lm(size.phyto ~ log(NPP), data = data.revi))
# # With FPOCex? YES
# summary(lm(size.phyto ~ log(FPOCex), data = data.revi))
# ggplot(data = data.revi) + geom_point(aes(x = size.phyto, y = log(FPOCex)), alpha = .5, colour = "grey50") +
#     geom_smooth(aes(x = size.phyto, y = log(FPOCex)), colour = "red") +
#     ylab("FPOCex (log)") + xlab("Phytoplankton median cell diameter (µm)") + theme_classic()
# # With log FPOCex > 2???
# summary(lm(size.phyto ~ log(FPOCex), data = data.revi[which(log(data.revi$FPOCex)>2.8),]))
#
# # With slope? YES
# summary(lm(size.phyto ~ slope2, data = data.revi))
# # Multiple R-squared:  0.2831,    Adjusted R-squared:  0.2831
# # F-statistic: 1.382e+04 on 1 and 35000 DF,  p-value: < 2.2e-16
#
# ### Polynomial tests:
# summary(lm(size.phyto ~ slope2 + I(slope2^2) + I(slope2^3), data = data.revi))
# summary(lm(size.phyto ~ slope2 + I(slope2^2), data = data.revi))
#
# ggplot(data = data.revi) + geom_point(aes(x = size.phyto, y = slope2, colour = y), alpha = .5) +
#      geom_smooth(aes(x = size.phyto, y = slope2), colour = "black", method = "lm") +
#      scale_colour_distiller(name = "Latitude", palette = "RdYlBu") +
#      ylab("Plankton size (slope of PSD)") + xlab("Phytoplankton median cell diameter (µm)") + theme_classic()
# #
# ggplot(data = data.revi) + geom_point(aes(x = size.phyto, y = slope2), alpha = .5, colour = "grey50") +
#      geom_smooth(aes(x = size.phyto, y = slope2), colour = "red",
#          formula = size.phyto ~ slope2 + I(slope2^2) + I(slope2^3)) +
#      ylab("Plankton size (slope of PSD)") + xlab("Phytoplankton mean cell diameter (µm)") + theme_classic()
#
# ### What if we exclude the poles
# summary(lm(size.phyto ~ slope2, data = data.revi[abs(data.revi$y)<70,]))
#
#
# ### --> Works well with FPOCex, e ratio and PSD
#
# ### Analyze correlations cluster-wise
# names <- colnames(data.revi)[c(9:14,19:28,40,41)] ; names
# mydata <- data.revi[,names]
# # Change some colnames etc.
# colnames(mydata)[1] <- "Biodiversity"
# colnames(mydata)[4] <- "e ratio"
# colnames(mydata)[5] <- "PSD"
# colnames(mydata)[c(7:16)] <- str_replace_all(colnames(mydata)[c(7:16)] , "_base", "")
# colnames(mydata)[14] <- "Cnidaria"
# colnames(mydata)[18] <- "Phytoplankton diameter"
# # Transform the columns needed
# mydata$FPOCex <- log(mydata$FPOCex)
# mydata$NPP <- log(mydata$NPP)
# mydata$catches <- log1p(mydata$catches)
# mydata[,c(7:16)] <- log1p(mydata[,c(7:16)])
#
# res <- lapply(unique(mydata$clusters), function(i) {
#
#             message(paste("Computing corr matrix for cluster  ",i, sep = ""))
#
#             cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:16,18)]), method = "spearman"),2)
#             upper_tri <- get_upper_tri(cormat)
#             melted_cormat <- melt(upper_tri, na.rm = F)
#             colnames(melted_cormat) <- c("V1","V2","rho")
#             melted_cormat <- melted_cormat[melted_cormat$V2 == "Phytoplankton diameter",]
#             melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
#             melted_cormat$cluster <- i
#
#             return(melted_cormat)
#
#     } # eo lapply
#
# ) # eo lapply - i in clusters
# # Rbind
# table.corr <- dplyr::bind_rows(res)
# max.val <- max(abs(table.corr$rho)) ; max.val
#
# # Make sure 'V1' follows the order you want though
# table.corr$V1 <- factor(table.corr$V1, levels = unique(table.corr$V1))
#
# plot <- ggplot(table.corr, aes(factor(V1), factor(V2), fill = rho)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
#             midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(V1), factor(V2), label = rho), color = "black", size = 3) +
#     facet_wrap(.~factor(cluster), ncol = 2)
#
# setwd(WD)
# ggsave(plot = plot, filename = "corr_heat_mean_phyto_diameter_clusters_groups_facet.pdf", dpi = 300, width = 15, height = 7)
#
#
# ### OK, next step: evaluate changes in copepod med size per clusters between future and baseline
# setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
# files.fut
#
# ### For each files.base, retrieve copepod community composition and derive estiate of annual size structure
# #f <- files.fut[1] # for testing function below
# res <- mclapply(files.fut, function(f) {
#
#             # Message
#             message(paste(f, sep = ""))
#             comm <- get(load(f)) # dim(comm) ; colnames(comm)
#             # Colnames should match the 'taxon' vector
#             common.spp <- intersect(sizes$taxon,colnames(comm))
#             # length(common.spp)
#             subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
#             rm(comm) ; gc()
#             # Melt, to put species names as vector
#             m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
#             colnames(m.sub)[c(4,5)] <- c("species","HSI")
#             m.sub$size <- NA
#             for(sp in unique(m.sub$species)) {
#                 s <- sizes[sizes$taxon == sp,"idealized_diameter_in_micrometers_mean_of_sets"]
#                 m.sub[m.sub$species == sp,"size"] <- s
#             } # eo for loop
#             require("matrixStats") ; require("dplyr")
#             # weightedMedian(x = m.sub$size, w = m.sub$HSI)
#             mean.size <- data.frame(m.sub %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y),
#                             mean.size = weightedMean(x = size, w = HSI, na.rm = TRUE)))
#             # head(med.size) ; summary(med.size)
#             # Return with ESM, SDM etc.
#             terms <- do.call(cbind,strsplit(f, split = "_"))
#             mean.size$SDM <- terms[7,1]
#             mean.size$ESM <- terms[6,1]
#             mean.size$pool <- str_replace_all(terms[8,1],".Rdata","")
#
#             return(mean.size)
#
#         }, mc.cores = 40
#
# ) # eo mclapply - f in files.base
# # Rbind
# table.mean.size <- bind_rows(res) ; rm(res) ; gc()
# head(table.mean.size) ; dim(table.mean.size) ; summary(table.mean.size)
#
# # Derive ensemble estimate of median size structure across all ESM and SDM
# mean.size.fut <- data.frame(table.mean.size %>% group_by(cell_id) %>%
#         summarize(x = unique(x), y = unique(y), mean.size = mean(mean.size, na.rm = T)))
# summary(mean.size.fut)
#
# ggplot() + geom_raster(aes(x = x, y = y, fill = mean.size), data = mean.size.fut) +
#     geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = mean.size), data = mean.size.fut) +
#     scale_fill_viridis(name = "Phytoplankton mean\ncell diameter\n(µm)") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# ### Derive % difference in copepod median size
# # dim(mean.size) ; dim(mean.size.fut)
# commons <- intersect(data.revi$id, mean.size.fut$cell_id) ; length(commons)
# data.revi <- data.revi[order(data.revi$id),]
# mean.size.fut <- mean.size.fut[order(mean.size.fut$cell_id),]
# data.revi <- data.revi[data.revi$id %in% commons,]
# mean.size.fut <- mean.size.fut[mean.size.fut$cell_id %in% commons,]
# data.revi$size.phyto.fut <- mean.size.fut$mean.size
#
# ### Derive % change
# data.revi$diff.size.phyto <- (data.revi$size.phyto.fut)-(data.revi$size.phyto)
# data.revi$perc.size.phyto <- (data.revi$diff.size.phyto/(data.revi$size.phyto))*100
# summary(data.revi)
#
# map.perc.mean.size <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc.size.phyto), data = data.revi) +
#     geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = perc.size.phyto), data = data.revi) +
#     scale_fill_gradient2(name = "% Difference in\nphytoplankton median cell\ndiameter size (µm)", mid = "white", high = "#d53e4f", low = "#3288bd") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# ### And plot distrbition % difference in copepod median
# plot <- ggplot(data = data.revi, aes(x = factor(clusters), y = perc.size.phyto, fill = factor(clusters))) +
#    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") +
#    ylab(paste("% difference in phytoplankton median cell diameter (µm)", sep = "")) + xlab("") + theme_classic() +
#    geom_hline(yintercept = 0, linetype = "dashed")
#
# setwd(WD)
# ggsave(plot = map.perc.mean.size, filename = "map_ens_perc_mean_phyto_cell_diameter_2100-2000.pdf", width = 7, height = 4, dpi = 300)
# ggsave(plot = plot, filename = "plot_ens_perc_ann_mean_phyto_cell_diameter_clusters_facet.pdf", width = 6, height = 6, dpi = 300)




### 12/02/2021: Re-do the phytoplankton species size analysis but better: focus on Diatoms (like you focused on Copepods) and use the mean estimates of species-level cell volume, area, S/V ratio and carbon content from Leblanc et al. (2012)

library("raster")
library("tidyverse")
library("stringr")
library("reshape2")
library("RColorBrewer")
library("scales")
library("maps")
library("ggthemes")
library("viridis")

world2 <- map_data(map = "world2")
WD <- getwd()

setwd("/net/kryo/work/fabioben/OVERSEE/data/")
sizes <- read.table("table_traits_diatoms_species_modelled_11_02_21.txt", h = T, sep = "\t", dec = ",")
# dim(sizes) ; str(sizes) ; summary(sizes)
length(unique(sizes$Species)) # 156
# summary(factor(sizes$Level.of.meas))
# ggplot(aes(x = log10(Volume_mean), y = A.V_mean), data = sizes[sizes$A.V_mean < 4,]) +
#     geom_point(colour = "black") +
#     # geom_smooth(aes(x = log10(Volume_mean), y = A.V_mean),
# #         colour = "red", se = T, method = "lm", formula = y ~ x + I(x^2)) +
#     xlab("Mean cell violume log10(µm3)") + ylab("Mean cell S/V (µm-1)") +
#     theme_classic()

### Weird: S/V ratios above 10
sizes[order(sizes$A.V_mean , decreasing = T),]
sizes <- sizes[sizes$A.V_mean < 5,]
sizes <- sizes %>% drop_na(A.V_mean)
#dim(sizes)

### Then, fetch mean annual baseline (and then future) copepod community composition
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
files.base <- dir()[grep("ann_compo_phyto_baseline",dir())]
files.fut <- dir()[grep("ann_compo_phyto_2100-2000",dir())] # for later

### For each files.base, retrieve phytoplankton community composition and estimate of annual cell dimater structure 
require("parallel")

# f <- files.base[8] # for testing function below
res <- mclapply(files.base, function(f) {
            
            # Message
            message(paste(f, sep = ""))
            comm <- get(load(f)) # dim(comm) ; colnames(comm)
            # Change colnames (remove brackets or points)
            colnames(comm) <-  gsub("\\.|\\.", "", colnames(comm))
            # Colnames of copeopods should match the 'specieswithsizes' vector
            common.spp <- intersect(sizes$Species, colnames(comm))
            # % of species with measurements relative to all phyto species modelled?
            # length(common.spp) / length(colnames(comm)[c(4:length(comm))])
            ### --> 45%
                        
            subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
            rm(comm) ; gc()
            # Melt, to put species names as vector
            m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
            colnames(m.sub)[c(4,5)] <- c("species","HSI")
            m.sub$volume <- NA
            m.sub$SV <- NA
            m.sub$carbon <- NA
            for(sp in unique(m.sub$species)) {
                m.sub[m.sub$species == sp,"volume"] <- sizes[sizes$Species == sp,"Volume_mean"]
                m.sub[m.sub$species == sp,"SV"] <- sizes[sizes$Species == sp,"A.V_mean"]
                m.sub[m.sub$species == sp,"carbon"] <- sizes[sizes$Species == sp,"Carbon_mean"]
            } # eo for loop
            # summary(m.sub)
            require("matrixStats") ; require("dplyr")
            # weightedMedian(x = m.sub$size, w = m.sub$HSI)
            comm.stats <- data.frame(m.sub %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y),
                            comm.vol = weightedMedian(x = volume, w = HSI, na.rm = T), 
                            comm.vol.var = weightedVar(x = volume, w = HSI, na.rm = T),
                            comm.SV = weightedMedian(x = SV, w = HSI, na.rm = T), 
                            comm.SV.var = weightedVar(x = SV, w = HSI, na.rm = T),
                            comm.carb = weightedMedian(x = carbon, w = HSI, na.rm = T), 
                            comm.carb.var = weightedVar(x = carbon, w = HSI, na.rm = T)
                    ) 
            )
                     
            return(comm.stats)
    
        }, mc.cores = length(files.base)
        
) # eo mclapply - f in files.base
# Rbind
comm.stats <- bind_rows(res) ; rm(res) ; gc()
head(comm.stats) ; dim(comm.stats) #; summary(comm.stats)
 
# Derive ensemble estimate of median size structure
mean.comm.stats <- data.frame(comm.stats %>% group_by(cell_id)
                %>% summarize(x = unique(x), y = unique(y),
                  comm.vol = median(comm.vol, na.rm = T), 
                  comm.vol.range = median(comm.vol.var, na.rm = T),
                  comm.SV = median(comm.SV, na.rm = T), 
                  comm.SV.range = median(comm.SV.var, na.rm = T),
                  comm.carb = median(comm.carb, na.rm = T),
                  comm.carb.range = median(comm.carb.var, na.rm = T)
                  ) # eo summarize
) # eo ddf
summary(mean.comm.stats)

map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.vol/10000), data = mean.comm.stats) +
    geom_contour(colour = "grey60", binwidth = 2, size = 0.25, aes(x = x, y = y, z = comm.vol/10000), data = mean.comm.stats) +
    scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.SV), data = mean.comm.stats) +
    geom_contour(colour = "grey60", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = comm.SV), data = mean.comm.stats) +
     scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.carb), data = mean.comm.stats) +
    geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = comm.carb), data = mean.comm.stats) +
     scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(comm.vol.range)), data = mean.comm.stats) +
    geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = log10(comm.vol.range)), data = mean.comm.stats) +
    scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.SV.range), data = mean.comm.stats) +
    geom_contour(colour = "grey60", binwidth = 0.20, size = 0.25, aes(x = x, y = y, z = comm.SV.range), data = mean.comm.stats) +
    scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.carb.range), data = mean.comm.stats) +
    geom_contour(colour = "grey60", binwidth = 100, size = 0.25, aes(x = x, y = y, z = comm.carb.range), data = mean.comm.stats) +
    scale_fill_viridis(name = "") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

setwd(WD)
ggsave(plot = map1, filename = "map_ens_ann_med_diatoms_cell_vol_baseline_10^5.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map2, filename = "map_ens_ann_med_diatoms_cell_SV_baseline.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map3, filename = "map_ens_ann_med_diatoms_cell_carbon_baseline.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map4, filename = "map_ens_ann_med_diatoms_cell_vol_range_log10_baseline.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map5, filename = "map_ens_ann_med_diatoms_cell_SV_range_baseline.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map6, filename = "map_ens_ann_med_diatoms_cell_carbon_range_baseline.pdf", width = 7, height = 4, dpi = 300)


### Combine 'med.size' with table_data_for_revisions_26_11_20
setwd("/net/kryo/work/fabioben/OVERSEE/data/")
data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
commons <- intersect(data.revi$id, mean.comm.stats$cell_id) ; length(commons)
data.revi <- data.revi[order(data.revi$id),]
mean.comm.stats <- mean.comm.stats[order(mean.comm.stats$cell_id),]
data.revi <- data.revi[data.revi$id %in% commons,]
mean.comm.stats <- mean.comm.stats[mean.comm.stats$cell_id %in% commons,]
# dim(mean.comm.stats) ; dim(data.revi)

# Cbind the 4 variables of Diatom community size structure you just estimated
data.revi[,c("Vol","SV","Carbon","Vol.range","SV.range","Carbon.range")] <- mean.comm.stats[,c("comm.vol","comm.SV","comm.carb","comm.vol.range","comm.SV.range","comm.carb.range")]

# Re-number clusters of severity propoerly
data.revi$clusters <- data.revi$pca_euclid_pam_k6
data.revi[data.revi$pca_euclid_pam_k6 == 1,"clusters"] <- 4
data.revi[data.revi$pca_euclid_pam_k6 == 2,"clusters"] <- 3
data.revi[data.revi$pca_euclid_pam_k6 == 3,"clusters"] <- 6
data.revi[data.revi$pca_euclid_pam_k6 == 4,"clusters"] <- 2
data.revi[data.revi$pca_euclid_pam_k6 == 6,"clusters"] <- 1

### OK, couple of things to do from here:
# - plot distribution of size.phyto per clusters (boxplots)
# - plot correlation heatmap bewteen size.phyto and the other variables (log SR and ES)
# - plot some bivariate plot for the strong and interesting corr found above (with POC, e ratio or PSD)

ggplot(data = data.revi, aes(x = factor(clusters), y = log(Vol), fill = factor(clusters))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") +
    ylab(paste("Diatom median cell volume log(µm3)", sep = "")) + xlab("") + theme_classic()
 
ggplot(data = data.revi, aes(x = factor(clusters), y = log(Vol.range), fill = factor(clusters))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") +
    ylab(paste("Diatom median cell volume range log(µm3)", sep = "")) + xlab("") + theme_classic()
 
ggplot(data = data.revi, aes(x = factor(clusters), y = get("S/V"), fill = factor(clusters))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") +
    ylab(paste("Diatom median cell S/V (µm-1)", sep = "")) + xlab("") + theme_classic()
 
ggplot(data = data.revi, aes(x = factor(clusters), y = Carbon, fill = factor(clusters))) +
    geom_boxplot(colour = "black") + scale_fill_brewer(name = "Clusters", palette = "RdYlBu") +
    ylab(paste("Diatom median cell mg C/m3", sep = "")) + xlab("") + theme_classic()

#  If ever nneeded, plots 3 and 4 are those that make most sense

### Corr coeff heatmap
get_lower_tri <- function(cormat){
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
}
get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
}
reorder_cormat <- function(cormat){
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <- cormat[hc$order, hc$order]
        return(cormat)
}

# names <- colnames(data.revi)[c(9:14,19:28,41:44)] ; names

### !! Add mean annual SiO2, NO3, N* and Si*
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d")
files <- dir()[grep("_18_09_20",dir())] ; files
clims <- lapply(files, function(f) {
            d <- read.table(f, h = T, sep = ";")
            # Rotate longitudes and add new ID
            d$x2 <- d$x
            d[d$x < 0 ,"x2"] <- (d[d$x < 0 ,"x"]) + 360
            d$id2 <- paste(d$x2, d$y, sep = "_")
            d <- d[,c("id2","x2","y","SST","logNO3","logSiO2","Nstar","Sistar","Pstar","MLD","PAR")]
            return(d)
    } # eo 
)
clims <- bind_rows(clims)

setwd(WD)
### Derive mean annual climatologies
ann.clim <- data.frame(clims %>% group_by(id2) %>%
                summarize(x = unique(x2), y = unique(y), SST = mean(SST, na.rm = T),
                NO3 = mean(logNO3, na.rm = T), SiO2 = mean(logSiO2, na.rm = T),
                Nstar = mean(Nstar, na.rm = T), Sistar = mean(Sistar, na.rm = T),
                Pstar = mean(Pstar, na.rm = T), MLD = mean(MLD, na.rm = T), PAR = mean(PAR, na.rm = T))
) # eo ddf
summary(ann.clim)

# OK, re-order, find common cells and cbind
ann.clim <- ann.clim[order(ann.clim$id),]
data.revi <- data.revi[order(data.revi$id),]
commons <- intersect(ann.clim$id, data.revi$id)
data.revi <- data.revi[data.revi$id %in% commons,]
ann.clim <- ann.clim[ann.clim$id2 %in% commons,]
data.revi[,c("SST","NO3","SiO2")] <- ann.clim[,c("SST","NO3","SiO2")]
### Add absolute Si* and N*
# data.revi[,"Si* (abs)"] <- abs(data.revi[,"Si*"])
# data.revi[,"N* (abs)"] <- abs(data.revi[,"N*"])
# summary(data.revi) #; dim(data.revi)

### Plot 2 zonal plots: median Volume and PSD --> need to check if lat patterns match observations
zonal <- data.frame(data.revi %>% group_by(y) %>%
            summarize(Vol.avg = mean(Vol, na.rm = T)/10000, Vol.sd = sd(Vol, na.rm=T)/10000,
             PSD = mean(slope2, na.rm = T), PSD.sd = sd(slope2, na.rm=T)
) )
summary(zonal)

p1 <- ggplot() + geom_ribbon(aes(y = y, xmin = Vol.avg-Vol.sd, xmax = Vol.avg+Vol.sd), fill = "grey75", data = zonal) +
    geom_path(aes(y = y, x = Vol.avg), data = zonal, colour = "black", linetype = "solid") +
    scale_x_continuous(name = paste('Mean (+sd) Diatom community\n median cell volume (10^5 µm3)', sep="")) +
    scale_y_continuous(position = "right", name = "", limits = c(-90,90),
    breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
    theme_classic() + geom_hline(yintercept = -60, linetype = "dashed") + 
    geom_hline(yintercept = 60, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dotted")
    
p2 <- ggplot() + geom_ribbon(aes(y = y, xmin = PSD-PSD.sd, xmax = PSD+PSD.sd), fill = "grey75", data = zonal) +
    geom_path(aes(y = y, x = PSD), data = zonal, colour = "black", linetype = "solid") +
    scale_x_continuous(name = paste('Mean (+sd) Particles Size Index\n(Kostadinov et al. 2009)', sep="")) +
    scale_y_continuous(position = "right", name = "", limits = c(-90,90),
    breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
    theme_classic() + geom_hline(yintercept = -60, linetype = "dashed") + 
    geom_hline(yintercept = 60, linetype = "dashed")+ geom_hline(yintercept = 0, linetype = "dotted")
#
require("ggpubr")
# ggarrange(p1,p2,ncol = 2, nrow = 1, align = "hv", labels = c("A","B"))
setwd("/net/kryo/work/fabioben/OVERSEE/data/")
panel <- ggarrange(p1, p2, ncol = 2, nrow = 1, align = "hv")
#  panel
ggsave(plot = panel, filenam = "panel_zonal_ann_base_diatoms_vol+PSD.jpg", dpi = 300, width = 6, height = 5)


### To re-check regional correlations
data.revi2 <- data.revi %>% drop_na(Vol, slope2)
# Global:
cor(log(data.revi2[,"Vol"]), data.revi2[,"slope2"], method = "spearman") # 0.26
cor(data.revi2[,"Vol"]/10000, data.revi2[,"slope2"], method = "spearman") # 0.26
# Without ANY poles
cor(log(data.revi2[abs(data.revi2$y) < 60,"Vol"]), data.revi2[abs(data.revi2$y) < 60,"slope2"], method = "spearman") # 0.5357
# ONLY with poles
cor(log(data.revi2[abs(data.revi2$y) > 60,"Vol"]), data.revi2[abs(data.revi2$y) > 60,"slope2"], method = "spearman") # 0.48
# In the Arctic Ocean:
cor(log(data.revi2[data.revi2$y > 60,"Vol"]), data.revi2[data.revi2$y > 60,"slope2"], method = "spearman") # 0.30
cor(log(data.revi2[data.revi2$y > 70,"Vol"]), data.revi2[data.revi2$y > 70,"slope2"], method = "spearman") # 0.28
cor(log(data.revi2[data.revi2$y > 80,"Vol"]), data.revi2[data.revi2$y > 80,"slope2"], method = "spearman") # NS
# In the Southern Ocean
cor(log(data.revi2[data.revi2$y < -60,"Vol"]), data.revi2[data.revi2$y < -60,"slope2"], method = "spearman") # 0.12
cor(log(data.revi2[data.revi2$y < -65,"Vol"]), data.revi2[data.revi2$y < -65,"slope2"], method = "spearman") # 0.47
cor(log(data.revi2[data.revi2$y < -70,"Vol"]), data.revi2[data.revi2$y < -70,"slope2"], method = "spearman") # 0.61
# Southern Hemisphere
cor(log(data.revi2[data.revi2$y > 0,"Vol"]), data.revi2[data.revi2$y > 0,"slope2"], method = "spearman") # 0.31
# N Hemisphre
cor(log(data.revi2[data.revi2$y < 0,"Vol"]), data.revi2[data.revi2$y < 0,"slope2"], method = "spearman") # 0.24

summary(lm(slope2 ~ log(Vol), data = data.revi[abs(data.revi$y) < 70,]))
# Adjusted R-squared:  0.14
summary(lm(slope2 ~ log(Vol), data = data.revi[abs(data.revi$y) >= 70,]))
# 0.13

# lat.thresh <- 70
# plot <- ggplot() + geom_point(aes(x = log(Vol), y = slope2), alpha = .3, data = data.revi[abs(data.revi$y) < lat.thresh,], colour = "#d6604d") +
#      geom_point(aes(x = log(Vol), y = slope2), alpha = .3, data = data.revi[abs(data.revi$y) >= lat.thresh,], colour = "#4393c3") +
#      geom_smooth(aes(x = log(Vol), y = slope2), colour = "black", se = T,
#         method = "lm", data = data.revi[abs(data.revi$y) < lat.thresh,], linetype = "dashed") +
#      geom_smooth(aes(x = log(Vol), y = slope2), colour = "black", se = T,
#         method = "lm", data = data.revi[abs(data.revi$y) >= lat.thresh,], linetype = "longdash") +
#     geom_smooth(aes(x = log(Vol), y = slope2), colour = "black", se = T, method = "lm", data = data.revi) +
#      ylab("Particles size index") + xlab("Diatoms median cell volume log(µm3)") + theme_classic()
#
# ggsave(plot = plot, filename = "plot_ann_base_med_diatoms_volxPSD_fits.jpg", dpi = 300, width = 5, height = 4.5)
#


### Heatmap of rank correlations
names <- colnames(data.revi)[c(10:13,47:49,19,41:46)] ; names

# For global analysis
#mydata <- data.revi[,names]

# For filtering out high latitudes
# mydata <- data.revi[data.revi$y > -60,names]
mydata <- data.revi[data.revi$y < 63,names]
# Change some colnames etc.
colnames(mydata)
#colnames(mydata)[1] <- "Biodiversity"
colnames(mydata)[3] <- "e ratio"
colnames(mydata)[4] <- "PSI"
colnames(mydata)[8] <- "Diatoms SR"
colnames(mydata)[9] <- "Diatoms median\nvolume"
colnames(mydata)[10] <- "Diatoms median\nS/V ratio"
colnames(mydata)[11] <- "Diatoms median\nC content"
colnames(mydata)[12] <- "Diatoms median\nvolume diversity"
colnames(mydata)[13] <- "Diatoms median\nS/V ratio diversity"
colnames(mydata)[14] <- "Diatoms median\nC content diversity"

stats <- colnames(mydata)[c(8:14)] ; stats
env <- colnames(mydata)[c(1:7)]

# Transform the columns needed
mydata$FPOCex <- log(mydata$FPOCex)
mydata$NPP <- log(mydata$NPP)
# Log transform cell volume and volume range
mydata[,"Diatoms median\nvolume"] <- mydata[,"Diatoms median\nvolume"]/10000
mydata[,"Diatoms median\nvolume diversity"] <- log10(mydata[,"Diatoms median\nvolume diversity"])

# Corr matrix
cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = T)
#colnames(melted_cormat) <- c("V1","V2","rho")
### Separate Diatom community size structure indices from other covariates
 # diatoms <- stats
 # env <- colnames(mydata)[!(stats)]
# melted_cormat <- melted_cormat[which(melted_cormat$V1 %in% diatoms & melted_cormat$V2 %in% env),]
colnames(melted_cormat) <- c("ES","Group","rho")
melted_cormat <- melted_cormat[!(melted_cormat$Group == melted_cormat$ES),]
#melted_cormat <- melted_cormat[(melted_cormat$ES %in% stats),]
#melted_cormat <- melted_cormat[(melted_cormat$Group %in% c(env,diatoms)),]

max.val <- max(abs(melted_cormat$rho)) ; max.val

heatmap <- ggplot(melted_cormat, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
            midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
    xlab("") + ylab("") + theme_minimal() + coord_fixed() +
    geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#
setwd("/net/kryo/work/fabioben/OVERSEE/data/")
ggsave(plot = heatmap, filename = "heatmap_corr_ann_base_ens_diatoms_volsxenv_nopoles.jpg", dpi = 300, width = 7, height = 7)

### Examine bivariate plots
# With SST
summary(lm(log(Vol) ~ SST, data = data.revi))
ggplot(data = data.revi) + geom_point(aes(x = SST, y = log(Vol)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = SST, y = log(Vol)), colour = "red", se = T, method = "gam") + 
    xlab("SST (°C)") + ylab("Diatoms median cell volume log(µm3)") + theme_classic() 
#
ggplot(data = data.revi) + geom_point(aes(x = SST, y = get("S/V")), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = SST, y = get("S/V")), colour = "red", se = T, method = "gam") + 
    xlab("SST (°C)") + ylab("Diatoms median cell S/V (µm-1)") + theme_classic() 

# With e ratio? No

# With NPP? Yes
summary(lm(log(NPP) ~ log(Vol), data = data.revi)) # Adjusted R-squared: 0.2466
ggplot(data = data.revi) + geom_point(aes(x = log(Vol), y = log(NPP)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(Vol), y = log(NPP)), colour = "red", se = T, method = "lm") + 
    ylab("NPP (log)") + xlab("Diatoms median cell volume log(µm3)") + theme_classic() 

# With FPOCex? YES
summary(lm(log(FPOCex) ~ log(Vol), data = data.revi)) # Adjusted R-squared: 0.2343 
ggplot(data = data.revi) + geom_point(aes(x = log(Vol), y = log(FPOCex)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(Vol), y = log(FPOCex)), colour = "red", se = T, method = "lm") + 
    ylab("FPOCex (log)") + xlab("Diatoms median cell volume log(µm3)") + theme_classic() 

# With slope of PSD?
summary(lm(slope2 ~ log(Vol), data = data.revi))
ggplot(data = data.revi) + geom_point(aes(x = log(Vol), y = slope2, colour = y), alpha = .5) +
     geom_smooth(aes(x = log(Vol), y = slope2), colour = "black", se = T, method = "lm", formula = y ~ x+I(x^2)) +
     scale_colour_distiller(palette = "RdYlBu", name = "Latitude") + 
     ylab("Particles size index") + xlab("Diatoms median cell volume log(µm3)") + theme_classic()

# Without poles
summary(lm(slope2 ~ log(Vol), data = data.revi[abs(data.revi$y)<60,] )) # 0.2313 

ggplot(data = data.revi) + geom_point(aes(x = log(Vol), y = slope2, colour = y), alpha = .5) +
     geom_smooth(aes(x = log(Vol), y = slope2), colour = "black", se = T, method = "lm") +
     scale_colour_distiller(palette = "RdYlBu", name = "Latitude") + 
     ylab("Particles size index") + xlab("Diatoms median cell volume log(µm3)") + theme_classic()
#
ggplot(data = data.revi[abs(data.revi$y) < 65,]) + geom_point(aes(x = log(Vol), y = slope2, colour = y), alpha = .5) +
     geom_smooth(aes(x = log(Vol), y = slope2), colour = "black", se = T, method = "lm") +
     scale_colour_distiller(palette = "RdYlBu", name = "Latitude") + 
     ylab("Particles size index") + xlab("Diatoms median cell volume log(µm3)") + theme_classic()

### Si* vs. volume and volume range
summary(lm( get("Si* (abs)") ~ log(Vol), data = data.revi)) # Adjusted R-squared:  0.19 
#summary(lm(get("Si*") ~ log(Vol), data = data.revi[abs(data.revi$y)<70,]))
ggplot(data = data.revi) + geom_point(aes(x = log(Vol), y = get("Si* (abs)")), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(Vol), y = get("Si* (abs)")), colour = "black", se = T, method = "gam") + 
    ylab("Absolute Si*") + xlab("Diatoms median cell volume log(µm3)") + theme_classic() 
#
ggplot(data = data.revi) + geom_point(aes(x = log(Vol.range), y = get("Si* (abs)")), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(Vol.range), y = get("Si* (abs)")), colour = "black", se = T, method = "gam") + 
    ylab("Absolute Si*") + xlab("Diatoms median cell volume range, log(µm3)") + theme_classic() 
#
ggplot(data = data.revi) + geom_point(aes(x = get("S/V"), y = get("Si* (abs)")), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = get("S/V"), y = get("Si* (abs)")), colour = "black", se = T, method = "gam") + 
    ylab("Absolute Si*") + xlab("Diatoms median S/V (µm-1)") + theme_classic() 
#
ggplot(data = data.revi) + geom_point(aes(x = Bacillariophyceae_base, y = get("Si* (abs)")), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = Bacillariophyceae_base, y = get("Si* (abs)")), colour = "black", se = T, method = "gam") + 
    ylab("Absolute Si*") + xlab("Diatoms species richness") + theme_classic() 
#
ggplot(data = data.revi) + geom_point(aes(x = Bacillariophyceae_base, y = log(Vol)), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = Bacillariophyceae_base, y = log(Vol)), colour = "black", se = T, method = "lm") + 
    ylab("Diatoms median cell volume (µm3)") + xlab("Diatoms species richness") + theme_classic() 
#
ggplot(data =  data.revi[abs(data.revi$y) < 60,]) +
    geom_point(aes(x = log(Bacillariophyceae_base), y = get("S/V")), alpha = .5, colour = "grey50") +
    geom_smooth(aes(x = log(Bacillariophyceae_base), y = get("S/V")), colour = "black", se = T, method = "lm") + 
    ylab("Diatoms median cell S/V") + xlab("Diatoms species richness") + theme_classic() 


# ### Check the position of mean cell volume and mean S/V ratio in env niche space like the C-S-R model or Margalef's mandala
# ggplot(data = data.revi, aes(x= PAR, y= SiO2, colour = log(Vol))) +
#      geom_point(alpha = .5) + scale_colour_viridis(name = "Diatom median\ncell volume\n(µm3)") +
#      xlab("PAR") + ylab("SiO2 (µM)") + theme_classic() + scale_x_reverse()
#
# ggplot(data = data.revi, aes(x= PAR, y= SiO2, colour = get("S/V") )) +
#      geom_point(alpha = .5) + scale_colour_viridis(name = "Diatom median\nS/V (µm-1)") +
#      xlab("PAR") + ylab("SiO2 (µM)") + theme_classic() + scale_x_reverse()

# ### Analyze correlations cluster-wise
# names <- colnames(data.revi)[c(9:14,19:28,40:44)] ; names
# mydata <- data.revi[,names]
# colnames(mydata)[1] <- "Biodiversity"
# colnames(mydata)[4] <- "e ratio"
# colnames(mydata)[5] <- "PSD"
# colnames(mydata)[c(7:16)] <- str_replace_all(colnames(mydata)[c(7:16)] , "_base", "")
# colnames(mydata)[18] <- "Diatoms mean\ncell volume"
# colnames(mydata)[19] <- "Diatoms mean\ncell volume div."
# colnames(mydata)[20] <- "Diatoms mean\nS/V ratio"
# colnames(mydata)[21] <- "Diatoms mean\nCarbon content"
# stats <- c("Diatoms mean\ncell volume","Diatoms mean\ncell volume div.","Diatoms mean\nS/V ratio","Diatoms mean\nCarbon content")
# # Transform the columns needed
# mydata$FPOCex <- log(mydata$FPOCex)
# mydata$NPP <- log(mydata$NPP)
# mydata$catches <- log1p(mydata$catches)
# mydata[,c(7:16)] <- log1p(mydata[,c(7:16)])
# # Log transform cell volume and volume range
# mydata[,"Diatoms mean\ncell volume"] <- log(mydata[,"Diatoms mean\ncell volume"])
# mydata[,"Diatoms mean\ncell volume div."] <- log(mydata[,"Diatoms mean\ncell volume div."])
#
# res <- lapply(unique(mydata$clusters), function(i) {
#
#             message(paste("Computing corr matrix for cluster  ",i, sep = ""))
#
#             cormat <- round(cor(na.omit(mydata[mydata$clusters == i,c(1:16,18)]), method = "spearman"),2)
#
#             upper_tri <- get_upper_tri(cormat)
#             melted_cormat <- melt(upper_tri, na.rm = T)
#             colnames(melted_cormat) <- c("V1","V2","rho")
#             melted_cormat <- melted_cormat[!(melted_cormat$V1 == melted_cormat$V2),]
#             melted_cormat <- melted_cormat[!(melted_cormat$V1 %in% stats),]
#             melted_cormat <- melted_cormat[melted_cormat$V2 %in% stats,]
#             colnames(melted_cormat) <- c("ES","Group","rho")
#
#             melted_cormat$cluster <- i
#
#             return(melted_cormat)
#
#     } # eo lapply
#
# ) # eo lapply - i in clusters
# # Rbind
# table.corr <- dplyr::bind_rows(res)
# max.val <- max(abs(table.corr$rho)) ; max.val
#
# # Make sure 'V1' follows the order you want though
# table.corr$ES <- factor(table.corr$ES, levels = unique(table.corr$ES))
#
# ggplot(table.corr, aes(factor(ES), factor(Group), fill = rho)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white",
#             midpoint = 0, limits = c(max.val*-1,max.val), name = "Spearman's\ncorrelation\ncoefficient") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(ES), factor(Group), label = rho), color = "black", size = 3) +
#     facet_wrap(.~factor(cluster), ncol = 2)
#
# setwd(WD)
# ggsave(plot = plot, filename = "corr_heat_mean_phyto_diameter_clusters_groups_facet.pdf", dpi = 300, width = 15, height = 7)


### OK, next step: evaluate changes in copepod med size per clusters between future and baseline
setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
files.fut 

### For each files.base, retrieve copepod community composition and derive estiate of annual size structure 
#f <- files.fut[1] # for testing function below
res <- mclapply(files.fut, function(f) {
            
            # Message
            message(paste(f, sep = ""))
            comm <- get(load(f)) # dim(comm) ; colnames(comm)
            # Colnames should match the 'taxon' vector
            common.spp <- intersect(sizes$Species,colnames(comm))
            # length(common.spp) 
            subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
            rm(comm) ; gc()
            
            m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
            colnames(m.sub)[c(4,5)] <- c("species","HSI")
            m.sub$volume <- NA
            m.sub$SV <- NA
            m.sub$carbon <- NA
            for(sp in unique(m.sub$species)) {
                m.sub[m.sub$species == sp,"volume"] <- sizes[sizes$Species == sp,"Volume_mean"]
                m.sub[m.sub$species == sp,"SV"] <- sizes[sizes$Species == sp,"A.V_mean"]
                m.sub[m.sub$species == sp,"carbon"] <- sizes[sizes$Species == sp,"Carbon_mean"]
            } # eo for loop

            require("matrixStats") ; require("dplyr")
            # weightedMedian(x = m.sub$size, w = m.sub$HSI)
            comm.stats <- data.frame(m.sub %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y),
                            comm.vol = weightedMedian(x = volume, w = HSI, na.rm = T), 
                            comm.vol.var = weightedVar(x = volume, w = HSI, na.rm = T),
                            comm.SV = weightedMedian(x = SV, w = HSI, na.rm = T), 
                            comm.SV.var = weightedVar(x = SV, w = HSI, na.rm = T),
                            comm.carb = weightedMedian(x = carbon, w = HSI, na.rm = T), 
                            comm.carb.var = weightedVar(x = carbon, w = HSI, na.rm = T)
                    ) 
            )
            
            # Return with ESM, SDM etc.
            terms <- do.call(cbind,strsplit(f, split = "_"))
            comm.stats$SDM <- terms[7,1]
            comm.stats$ESM <- terms[6,1]
            comm.stats$pool <- str_replace_all(terms[8,1],".Rdata","")
            
            return(comm.stats)
    
        }, mc.cores = 40
        
) # eo mclapply - f in files.base
# Rbind
comm.stats.fut <- bind_rows(res) ; rm(res) ; gc()
head(comm.stats.fut) ; dim(comm.stats.fut)
# summary(comm.stats.fut)

# Derive ensemble estimate of median size structure
mean.comm.stats.fut <- data.frame(comm.stats.fut %>% group_by(cell_id) %>%
                summarize(x = unique(x), y = unique(y),
                  comm.vol = median(comm.vol, na.rm = T), 
                  comm.vol.range = median(comm.vol.var, na.rm = T),
                  comm.SV = median(comm.SV, na.rm = T), 
                  comm.SV.range = median(comm.SV.var, na.rm = T),
                  comm.carb = median(comm.carb, na.rm = T),
                  comm.carb.range = median(comm.carb.var, na.rm = T)
                ) # eo summarize
) # eo ddf

summary(mean.comm.stats$comm.SV)
summary(mean.comm.stats.fut$comm.SV)

### Maps
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(comm.vol)), data = mean.comm.stats.fut) +
    geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = log(comm.vol)), data = mean.comm.stats.fut) +
    scale_fill_viridis(name = "Diatoms median\ncell volume\nlog(µm3)") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.SV), data = mean.comm.stats.fut) +
    geom_contour(colour = "grey60", binwidth = 0.2, size = 0.25, aes(x = x, y = y, z = comm.SV), data = mean.comm.stats.fut) +
    scale_fill_viridis(name = "Diatoms median\nS/V (µm-1)") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

#
map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = comm.carb), data = mean.comm.stats.fut) +
    geom_contour(colour = "grey60", binwidth = 50, size = 0.25, aes(x = x, y = y, z = comm.carb), data = mean.comm.stats.fut) +
    scale_fill_viridis(name = "Diatoms median\ncell mgC/m3") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(comm.vol.range)), data = mean.comm.stats.fut) +
    geom_contour(colour = "grey60", binwidth = 1, size = 0.25, aes(x = x, y = y, z = log(comm.vol.range)), data = mean.comm.stats.fut) +
    scale_fill_viridis(name = "Diatoms median\nvolume range\nlog(µm3)") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


setwd(WD)
ggsave(plot = map1, filename = "map_ens_ann_med_diatoms_cell_vol_2100-2012.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map2, filename = "map_ens_ann_med_diatoms_cell_SV_2100-2012.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map3, filename = "map_ens_ann_med_diatoms_cell_carbon_2100-2012.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map4, filename = "map_ens_ann_med_diatoms_cell_vol_range_2100-2012.pdf", width = 7, height = 4, dpi = 300)

### Derive % difference in copepod median size
# dim(mean.size) ; dim(mean.size.fut)
commons <- intersect(data.revi$id, mean.comm.stats.fut$cell_id) ; length(commons)
data.revi <- data.revi[order(data.revi$id),]
mean.comm.stats.fut <- mean.comm.stats.fut[order(mean.comm.stats.fut$cell_id),]
data.revi <- data.revi[data.revi$id %in% commons,]
mean.comm.stats.fut <- mean.comm.stats.fut[mean.comm.stats.fut$cell_id %in% commons,]

# Cbind the 4 variables of Diatom community size structure you just estimated
colnames(mean.comm.stats.fut)
data.revi[,c("Vol.fut","SV.fut","Carbon.fut","Vol.range.fut","SV.range.fut","Carbon.range.fut")] <- mean.comm.stats.fut[,c("comm.vol","comm.SV","comm.carb","comm.vol.range","comm.SV.range","comm.carb.range")]

### Derive % changes
# Vol
data.revi$diff.vol <- (data.revi$Vol.fut)-(data.revi$Vol)
data.revi$perc.vol <- (data.revi$diff.vol/(data.revi$Vol))*100
# S/V
data.revi$diff.SV <- (data.revi$SV.fut)-(data.revi$SV)
data.revi$perc.SV <- (data.revi$diff.SV/(data.revi$SV))*100
# Carbon
data.revi$diff.carb <- (data.revi$Carbon.fut)-(data.revi$Carbon)
data.revi$perc.carb <- (data.revi$diff.carb/(data.revi$Carbon))*100

# And their diff in ranges
data.revi$diff.vol.range <- (data.revi$Vol.range.fut)-(data.revi$Vol.range)
data.revi$perc.vol.range <- (data.revi$diff.vol.range/(data.revi$Vol.range))*100

data.revi$diff.SV.range <- (data.revi$SV.range.fut)-(data.revi$SV.range)
data.revi$perc.SV.range <- (data.revi$diff.SV.range/(data.revi$SV.range))*100

data.revi$diff.Carbon.range <- (data.revi$Carbon.range.fut)-(data.revi$Carbon.range)
data.revi$perc.Carbon.range <- (data.revi$diff.Carbon.range/(data.revi$Carbon.range))*100

summary(data.revi)

map.perc.volume <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc.vol),
            data = data.revi[data.revi$perc.vol < 50,]) +
    geom_raster(aes(x = x, y = y), data = data.revi[data.revi$perc.vol > 50,], fill = "#d53e4f") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc.vol), data = data.revi) +
    geom_point(aes(x = x, y = y), data = data.revi[abs(data.revi$y) > 60,], colour = "grey55", alpha = 0.07) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    geom_hline(yintercept = 60, linetype = "dashed") + geom_hline(yintercept = -60, linetype = "dashed")
#
map.perc.SV <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc.SV), data = data.revi) +
    geom_contour(colour = "grey60", binwidth = 5, size = 0.25, aes(x = x, y = y, z = perc.SV), data = data.revi) +
    geom_point(aes(x = x, y = y), data = data.revi[abs(data.revi$y) > 60,], colour = "grey55", alpha = 0.07) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    geom_hline(yintercept = 60, linetype = "dashed") + geom_hline(yintercept = -60, linetype = "dashed")
#
map.perc.carb <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc.carb), data = data.revi) +
    geom_contour(colour = "grey60", binwidth = 10, size = 0.25, aes(x = x, y = y, z = perc.carb), data = data.revi) +
    geom_point(aes(x = x, y = y), data = data.revi[abs(data.revi$y) > 60,], colour = "grey55", alpha = 0.07) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    geom_hline(yintercept = 60, linetype = "dashed") + geom_hline(yintercept = -60, linetype = "dashed")

setwd("/net/kryo/work/fabioben/OVERSEE/data/")
ggsave(plot = map.perc.volume, filename = "map_ens_perc_med_diatoms_volume_2100-2000_masked.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map.perc.SV, filename = "map_ens_perc_med_diatoms_SV_2100-2000_masked.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map.perc.carb, filename = "map_ens_perc_med_diatoms_carbon_2100-2000_masked.pdf", width = 7, height = 4, dpi = 300)


### !!! Maps of ∆SR in Copepod SR & Diatoms
data.revi$Copepods.perc <- (data.revi$Copepoda_fut/(data.revi$Copepoda_base))*100
data.revi$Diatoms.perc <- (data.revi$Bacillariophyceae_fut/(data.revi$Bacillariophyceae_base))*100
summary(data.revi)

map.perc.copepoda <- ggplot() + geom_raster(aes(x = x, y = y, fill = Copepods.perc), data = data.revi[data.revi$Copepods.perc <= 40,]) +
    geom_raster(aes(x = x, y = y), data = data.revi[data.revi$Copepods.perc > 40,], fill = "#d53e4f") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = Copepods.perc), data = data.revi) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) 
#
map.perc.diatoms <- ggplot() + geom_raster(aes(x = x, y = y, fill = Diatoms.perc), data = data.revi[data.revi$Diatoms.perc <= 50,]) +
    geom_raster(aes(x = x, y = y), data = data.revi[data.revi$Diatoms.perc > 50,], fill = "#d53e4f") +
    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = Diatoms.perc), data = data.revi) +
    scale_fill_gradient2(name = "", mid = "white", high = "#d53e4f", low = "#3288bd") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) 

ggsave(plot = map.perc.copepoda, filename = "map_ens_perc_rich_Copepods_2100-2000_masked.pdf", width = 7, height = 4, dpi = 300)
ggsave(plot = map.perc.diatoms, filename = "map_ens_perc_rich_Diatoms_2100-2000_masked.pdf", width = 7, height = 4, dpi = 300)



# --------------------------------------------------------------------------------------------------------------------------------

### 18/05/21: Checking for longitudinal patterns in our LDGs and future changes in SR
data.revi <- read.table("table_data_for_revisions_26_11_20.txt", h = T, sep = "\t")
colnames(data.revi)
summary(data.revi)

### First: baseline phyto
zonal.phyto.base <- data.frame(data.revi %>% group_by(x) %>% summarize(avg_perc = mean(phyto_base, na.rm = T), sd = sd(phyto_base, na.rm = T)) )  
summary(zonal.phyto.base)

# To compare to lat variations
zonal2.phyto.base <- data.frame(data.revi %>% group_by(y) %>% summarize(avg_perc = mean(phyto_base, na.rm = T), sd = sd(phyto_base, na.rm = T)) )  
#quantile(zonal.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[1]
#quantile(zonal2.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal2.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[1]
max(zonal.phyto.base$avg_perc) - min(zonal.phyto.base$avg_perc)
max(zonal2.phyto.base$avg_perc) - min(zonal2.phyto.base$avg_perc)
# Ratio of: 2.01 

# Map zonal longitudinal plot
plota <- ggplot() + geom_point(aes(x = x, y = phyto_base), data = data.revi, alpha = .33, colour = "grey75") +
    geom_ribbon(aes(x = x, ymin = avg_perc-sd, ymax = avg_perc+sd), fill = "#d53e4f", data = zonal.phyto.base, alpha = .33) +
     geom_path(aes(x = x, y = avg_perc), data = zonal.phyto.base, colour = "black", linetype = "solid") +
     scale_y_continuous(name = "Mean annual phytoplanton species richness\n(contemporary ocean)") +
     scale_x_continuous(name = "Longitude", limits = c(0,360),
             breaks = c(0,60,120,180,240,300,360),
             labels = c("GM","60°E","120°E","180°","120°W","60°W","GM")) +
     theme_classic() 



### Second: baseline zoo
zonal.zoo.base <- data.frame(data.revi %>% group_by(x) %>% summarize(avg_perc = mean(zoo_base, na.rm = T), sd = sd(zoo_base, na.rm = T)) )  
summary(zonal.zoo.base)

# Map zonal longitudinal plot
plotb <- ggplot() + geom_point(aes(x = x, y = zoo_base), data = data.revi, alpha = .33, colour = "grey75") +
    geom_ribbon(aes(x = x, ymin = avg_perc-sd, ymax = avg_perc+sd), fill = "#d53e4f", data = zonal.zoo.base, alpha = .33) +
     geom_path(aes(x = x, y = avg_perc), data = zonal.zoo.base, colour = "black", linetype = "solid") +
     scale_y_continuous(name = "Mean annual zooplankton species richness\n(contemporary ocean)") +
     scale_x_continuous(name = "Longitude", limits = c(0,360),
             breaks = c(0,60,120,180,240,300,360),
             labels = c("GM","60°E","120°E","180°","120°W","60°W","GM")) +
     theme_classic() 

# To compare to lat variations
zonal2.zoo.base <- data.frame(data.revi %>% group_by(y) %>% summarize(avg_perc = mean(zoo_base, na.rm = T), sd = sd(zoo_base, na.rm = T)) )  
#quantile(zonal.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[1]
#quantile(zonal2.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal2.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[1]
max(zonal.zoo.base$avg_perc) - min(zonal.zoo.base$avg_perc)
max(zonal2.zoo.base$avg_perc) - min(zonal2.zoo.base$avg_perc)
# Ratio of: 2.16 

### Derive percentages ∆SR in phyto- and zooplankton SR from the diff and the base values
data.revi$perc.phyto <- (data.revi$phyto_fut/(data.revi$phyto_base))*100
data.revi$perc.zoo <- (data.revi$zoo_fut/(data.revi$zoo_base))*100
summary(data.revi)

### Third: %∆SR phyto
zonal.phyto.perc <- data.frame(data.revi %>% group_by(x) %>% summarize(avg_perc = mean(perc.phyto, na.rm = T), sd = sd(perc.phyto, na.rm = T)) )  
summary(zonal.phyto.perc)
# To compare to lat variations
zonal2.phyto.perc <- data.frame(data.revi %>% group_by(y) %>% summarize(avg_perc = mean(perc.phyto, na.rm = T), sd = sd(perc.phyto, na.rm = T)) )  
#quantile(zonal.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[1]
#quantile(zonal2.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal2.phyto.base$avg_perc, probs = seq(0, 1, 0.1))[1]
max(zonal.phyto.perc$avg_perc) - min(zonal.phyto.perc$avg_perc)
max(zonal2.phyto.perc$avg_perc) - min(zonal2.phyto.perc$avg_perc)
# Ratio of 4.49


# Map zonal longitudinal plot
plotc <- ggplot() + geom_point(aes(x = x, y = perc.phyto), data = data.revi, alpha = .33, colour = "grey75") +
    geom_ribbon(aes(x = x, ymin = avg_perc-sd, ymax = avg_perc+sd), fill = "#d53e4f", data = zonal.phyto.perc, alpha = .33) +
     geom_path(aes(x = x, y = avg_perc), data = zonal.phyto.perc, colour = "black", linetype = "solid") +
     geom_hline(yintercept = 0, linetype = "dashed") + 
     scale_y_continuous(name = "% difference in mean annual phytoplanton\nspecies richness (future - baseline)") +
     scale_x_continuous(name = "Longitude", limits = c(0,360),
             breaks = c(0,60,120,180,240,300,360),
             labels = c("GM","60°E","120°E","180°","120°W","60°W","GM")) +
     theme_classic() 

### Fourth: % ∆SR zoo
zonal.zoo.perc <- data.frame(data.revi %>% group_by(x) %>% summarize(avg_perc = mean(perc.zoo, na.rm = T), sd = sd(perc.zoo, na.rm = T)) )  
summary(zonal.zoo.perc)

# Map zonal longitudinal plot
plotd <- ggplot() + geom_point(aes(x = x, y = perc.zoo), data = data.revi, alpha = .33, colour = "grey75") +
    geom_ribbon(aes(x = x, ymin = avg_perc-sd, ymax = avg_perc+sd), fill = "#d53e4f", data = zonal.zoo.perc, alpha = .33) +
     geom_path(aes(x = x, y = avg_perc), data = zonal.zoo.perc, colour = "black", linetype = "solid") +
     geom_hline(yintercept = 0, linetype = "dashed") + 
     scale_y_continuous(name = "% difference in mean annual zooplankton\nspecies richness (future - baseline)") +
     scale_x_continuous(name = "Longitude", limits = c(0,360),
             breaks = c(0,60,120,180,240,300,360),
             labels = c("GM","60°E","120°E","180°","120°W","60°W","GM")) +
     theme_classic() 

# To compare to lat variations
zonal2.zoo.perc <- data.frame(data.revi %>% group_by(y) %>% summarize(avg_perc = mean(perc.zoo, na.rm = T), sd = sd(perc.zoo, na.rm = T)) )  
#quantile(zonal.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[1]
#quantile(zonal2.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[10] - quantile(zonal2.zoo.base$avg_perc, probs = seq(0, 1, 0.1))[1]
max(zonal.zoo.perc$avg_perc) - min(zonal.zoo.perc$avg_perc)
max(zonal2.zoo.perc$avg_perc) - min(zonal2.zoo.perc$avg_perc)
# Ratio of 2.42

### Save zonal plots in a panel
require("ggpubr") # ?ggarrange
panel <- ggarrange(plota, plotb, plotc, plotd, ncol = 2, nrow = 2, align = "hv", labels = c("a","b","c","d"))
panel
# Save panel as .jpg
ggsave(plot = panel, filename = "panel_long_gradients_mean_ann_rich+diff_phyto+zoo_reviews.jpg", dpi = 300, width = 12, height = 8)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------