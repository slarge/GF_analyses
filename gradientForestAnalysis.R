rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: Comparative approach for identifying ecosystem thresholds.
# AUTHOR: Scott Large 2015
# REVIEWED BY:
# VERSION: 0.1
#

######################
# CHANGES/ ADDITIONS #
######################
# Need to add: 


# Done:

############
# PACKAGES #
############
# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# install.packages("extendedForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
library(extendedForest)
library(reshape2)
# require(ggplot2)
##################
# Set up folders #
##################
setwd("~/git/slarge/gradientForest")
dump.dir <- "data/GF/"
figure.dir <- "output/"

######################
# Load any functions #
######################
source("gradForestFunctions.R")
#
#############
# Load data #
#############
set.seed(627)
cc <- read.csv("data/CCIEA-RPW.csv")
#
# Partion based on area... initally for coastwide
ccALL <- cc[cc$Coastwide == 1,]
ccALL$year <- as.numeric(ccALL$year)
#
ccALL$timeseries <- gsub("\\(", "", ccALL$timeseries)
ccALL$timeseries <- gsub("\\)", "", ccALL$timeseries)
ccALL$timeseries <- gsub(" ", "_", ccALL$timeseries)
ccALL$timeseries <- gsub("_-_", "_", ccALL$timeseries)
#
# Wide df with columns as variables and select only the columns with data
dat.full <- dcast(ccALL, year ~ timeseries, value.var = "value")
str(dat.full)
dat.full$year[!is.na(dat.full[colnames(dat.full) %in% ind.name])]


dat.sc <- dat.full[,apply(dat.full, 2, function(x)!all(is.na(x)))]
# if any column has less than 10 years, omit.
cuts <- 10
len.list <- sapply(dat.sc, function(x) length(na.contiguous(x)))
keep.list <- names(len.list[len.list >= cuts])
dat.kl <- dat.sc[, keep.list] 
dat <- dat.kl[apply(dat.kl, 1, function(x)!any(is.na(x))),]
#
colnames(dat)
# Indicator and Driver names
ind.name <- c("GF_spp_richness_coastwide", "GF-Simp_coastwide", "GF_MTL_coastwide", "Scav_ratio")
dri.name <- as.character(colnames(dat[!colnames(dat) %in% c(ind.name, "Ocean-based_pollution") ]))
dri.name <- dri.name[!dri.name == "year"]
dat <- dat[colnames(dat) %in% c(ind.name, dri.name)]
#
colnames(dat)
# Maximum level of splits
lev <- floor(log2(nrow(dat) * 0.368/2))
#
## GF analysis ##
gf <-   gradientForest(data = dat, 
                       predictor.vars = dri.name, 
                       response.vars = ind.name,
                       ntree = 10, 
                       transform = NULL,
                       maxLevel = lev,
                       corr.threshold = 0.5, 
                       compact = T,
                       trace = T)

save(gf, file = paste0(gfDump, "NEUSgf_v01.RDATA"))
imp.vars <- names(importance(gf))

# Collect PCA info
Trns_grid <- predict(gf, dat[, imp.vars])
row.names(Trns_grid) <- dat[,"YEAR"]
PCs <- prcomp(Trns_grid[, imp.vars])

save(PCs, file = paste0(gfDump, "NEUSpcaDAT_v01.RDATA"))

