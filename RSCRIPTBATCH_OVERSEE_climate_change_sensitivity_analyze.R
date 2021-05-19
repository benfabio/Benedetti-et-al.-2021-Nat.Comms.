
# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("cmocean")
library("RColorBrewer")
library("wesanderson")
library("ggthemes")

world2 <- map_data("world2")

# --------------------------------------------------------------------------------------------------------------------------------

### Set the working directories, vectors etc.
WD <- getwd()

setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
# Vector of month names
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# Vector of ESM
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
# Vector of vars to examine
vars <- c("sst","logno3","logsio2","o2","logchl")

### Examine distribution (violins + boxplot) of the delta (diff) for each and ESM
# v <- "sst"
# esm <- "CESM-BEC" 
library("parallel")
res <- mclapply(vars, function(v) {
            
            res2 <- lapply(ESMs, function(esm) {
                
                        # Useless message
                        message(paste("Retrieving delta clims for ",v," and for ",esm, sep = ""))
                        d <- get(load( paste("clims_mon_diff_",v,"_rcp85_",esm,"_2031-2100.Rdata", sep = "") ))
                        # dim(d); colnames(d)
                        md <- melt(d, id.vars = c("id","x","y"))
                        colnames(md)[c(4:5)] <- c("month","diff")
                        md$ESM <- esm
                        return(md)
                
                    } # eo FUN
                
            ) # eo lapply - ESM
            # Rbind
            table <- dplyr::bind_rows(res2)
            # dim(table) ; summary(table) 
            table$var <- v
            rm(res2); gc()
            return(table)
    
    }, mc.cores = 5
    
) # eo mclapply - v in vars
# Rbind
table <- dplyr::bind_rows(res)
dim(table) ; summary(table) 
unique(table$var)

setwd("/net/kryo/work/fabioben/OVERSEE/data/")

# quartz()
plot <- ggplot(data = na.omit(table[table$var == "sst",]), aes(x = factor(ESM), y = diff, fill = factor(ESM))) + 
    scale_fill_manual(name = "ESM", values = wes_palettes$Zissou1) + 
    geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", notch = F, width = .2) + 
    xlab("") + ylab("Modelled delta NO3 log(µg/l)") + theme_bw() + facet_wrap(~ month, ncol = 4, scales = "free_y") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    scale_y_continuous(limits = c(-0.3,0.3))
#
ggsave(plot = plot, filename = "plot_sensitivity_logno3_ESMs.jpg", dpi = 300, width = 14, height = 15)


### 09/04/2020: Following up on Gesa' question: check sensitivity to CC between NH and SH
head(table)
table$domain <- NA ; table[table$y >= 0,"domain"] <- "NH" ; table[table$y <= 0,"domain"] <- "SH"

plot <- ggplot(data = na.omit(table[table$var == "sst",]), aes(x = factor(domain), y = diff, fill = factor(domain))) + 
    scale_fill_manual(name = "Hemisphere", values = c("#3B9AB2","#F21A00")) + 
    geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", notch = F, width = .2) + 
    xlab("") + ylab("Modelled delta SST (°C)") + theme_bw() + facet_grid(factor(month) ~ factor(ESM), scales = "free_y") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    scale_y_continuous(limits = c(-3,10))
    
ggsave(plot = plot, filename = "plot_sensitivity_sst_ESMsxmonth.jpg", dpi = 300, width = 10, height = 30)



# --------------------------------------------------------------------------------------------------------------------------------

