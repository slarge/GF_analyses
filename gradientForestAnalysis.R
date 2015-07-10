install.packages("extendedForest", repos="http://R-Forge.R-project.org")
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

library(gradientForest)

cc <- read.csv("data/CCIEA-RPW.csv")

str(cc)
