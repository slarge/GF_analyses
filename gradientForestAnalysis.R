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
source("pcaBiplotFunction.R")
#
#############
# Load data #
#############
set.seed(627)
cc <- read.csv("data/CCIEA-RPW.csv")
#
# Subset area... initally for coastwide
ccALL <- cc[cc$Coastwide == 1,]
ccALL$year <- as.numeric(ccALL$year)
#
ccALL$timeseries <- gsub("\\(", "", ccALL$timeseries)
ccALL$timeseries <- gsub("\\)", "", ccALL$timeseries)
ccALL$timeseries <- gsub(" ", "_", ccALL$timeseries)
ccALL$timeseries <- gsub("_-_", "_", ccALL$timeseries)
#
# Wide df with columns as variables
dat.full <- dcast(ccALL, year ~ timeseries, value.var = "value")
#
######################################################################### 
# NOTE: It gets tricky here, b/c response variables are only 2003-2012, #
# whereas most pressure variables have longer timeseries or timeseries  # 
# that don't overlap.                                                   #
#########################################################################
#
# First cull the data to the years where we have response variables (may need to change
# for different area subsets)
dat.full <- dat.full[dat.full$year %in% c(2003:2012),]
#
# Many of the pressure variables do not have 2012 data:
## Option 1- Impute:
dat <- na.roughfix(dat.full)
#
## Option 2- Get rid of the columns without full time series:
#
# dat.sc <- dat.full[,apply(dat.full, 2, function(x)!all(is.na(x)))]
#
# Somewhat roundabout way of doing this, but can also be used to find the best
# contiguous set of data.
# if any column has less than XX years, omit.
#
# cuts <- 10
# len.list <- sapply(dat.sc, function(x) length(na.contiguous(x)))
# keep.list <- names(len.list[len.list >= cuts])
# dat.kl <- dat.sc[, keep.list] 
# dat <- dat.kl[apply(dat.kl, 1, function(x)!any(is.na(x))),]
#
## Option 3- Simply use as is (what I initially tried, but it might not run...)
# dat <- dat.full

#
# Indicator and Driver names
ind.name <- c("GF_spp_richness_coastwide", "GF-Simp_coastwide", "GF_MTL_coastwide", "Scav_ratio")
dri.name <- as.character(colnames(dat[!colnames(dat) %in% c(ind.name)]))
dri.name <- dri.name[!dri.name == "year"]

dat <- dat[colnames(dat) %in% c(ind.name, dri.name)]
#
# Maximum level of splits
lev <- floor(log2(nrow(dat) * 0.368/2))
#
## GF analysis ##
gf <-   gradientForest(data = dat, 
                       predictor.vars = dri.name, 
                       response.vars = ind.name,
                       ntree = 10000, 
                       transform = NULL,
                       maxLevel = lev,
                       corr.threshold = 0.5, 
                       compact = F,
                       trace = T)
#
save(gf, file = paste0(dump.dir, "CCE-gf_v01.RDATA"))
#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(gf, type = "Weighted", sort = TRUE))
#
indicatorIMP <- gf$imp.rsq
indicatorNA <- colnames(indicatorIMP)[apply(indicatorIMP, 2, function(x)all(is.na(x)))]
indicatorIMP <- indicatorIMP[apply(indicatorIMP, 2, function(x)!all(is.na(x))),]
#
indicatorDF <- data.frame("VARIABLE" = colnames(indicatorIMP),
                 "TYPE" = "INDICATOR", melt(apply(indicatorIMP, 2, mean, na.rm = T),
                                            variable.name = "VARIABLE",
                                            value.name = "MEAN"),
                 melt(apply(indicatorIMP, 2, min, na.rm = T),
                      variable.name = "VARIABLE", 
                      value.name = "MIN"),
                 melt(apply(indicatorIMP, 2, max, na.rm = T),
                      variable.name = "VARIABLE", 
                      value.name = "MAX"))
#
pressureIMP <- gf$imp.rsq
pressureNA <- row.names(pressureIMP)[apply(pressureIMP, 1, function(x)all(is.na(x)))]
pressureIMP <- pressureIMP[apply(pressureIMP, 1, function(x)!all(is.na(x))),]
pressureDF <- data.frame("VARIABLE" = row.names(pressureIMP),
                 "TYPE" = "PRESSURE", melt(apply(pressureIMP, 1, mean, na.rm = T),
                                           variable.name = "VARIABLE",
                                           value.name = "MEAN"),
                 melt(apply(pressureIMP, 1, min, na.rm = T),
                      variable.name = "VARIABLE", 
                      value.name = "MIN"),
                 melt(apply(pressureIMP, 1, max, na.rm = T),
                      variable.name = "VARIABLE", 
                      value.name = "MAX"))
#
indicatorDF[,c(3:5)] <- apply(indicatorDF[,c(3:5)], 2, signif, 2)
pressureDF[,c(3:5)] <- apply(pressureDF[,c(3:5)], 2, signif, 2)
# indicatorDF$VARIABLE <- factor(indicatorDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))
# pressureDF$VARIABLE <- factor(pressureDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))

pressureDF <- rbind(pressureDF, cbind(VARIABLE = pressureNA,TYPE = "PRESSURE", MEAN = "NA", MIN = "NA", MAX = "NA"))
#
modPerformance <- rbind(pressureDF, indicatorDF)
##
# write.csv(modPerformance, file = paste0(figure.dir, "CCE-gfModelPerformance_v001.csv"),
#           row.names = F)
#########
# PLOTS #
#########
#
imp.vars <- names(importance(gf)[importance(gf) > 0])
#
# Overall Importance
png(file = paste0(figure.dir, "CCE-overallImportance_v01.png"), width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(gf,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = paste0(figure.dir, "CCE-splitRatio_v01.png"), width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(gf,
        imp.vars = imp.vars,
#           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = paste0(figure.dir, "CCE-splitImportance_v01.png"), width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(gf,
             imp.vars = imp.vars,
#           imp.vars.names = #c(ADD NAMES HERE),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = paste0(figure.dir, "CCE-splitData_v01.png"), width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(gf,
       imp.vars = imp.vars,
#           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.png"), width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(gf, 
          imp.vars = imp.vars, 
#           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 6,
#           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = paste0(figure.dir, "CCE-cumImportance_v01.png"), width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(gf, 
          imp.vars = imp.vars, 
#           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 6,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(gf, dat[, imp.vars])
row.names(Trns_grid) <- c(2003:2012)
PCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(paste0(figure.dir, "CCE-PCA_v01.png"), PCbiplot(PCs),
       width = 83, height = 83, units = "mm", scale = 2)
#