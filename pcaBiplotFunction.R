PCbiplot <- function(PC, x="PC1", y="PC2") {
  require(ggplot2)
  require(grid)
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x]))))
  
  datapc <- transform(datapc,
                      v1 = 1.25 * mult * (get(x)),
                      v2 = 1.25 * mult * (get(y)))
  
#   limits <- rbind(xlim = c(min(data$PC1) + min(data$PC1) * .15 + min(datapc$v1) * 1,
#                            max(data$PC1) + max(data$PC1) * .15 + max(datapc$v1)* 1),
#                   ylim = c(min(data$PC2) + min(data$PC2) * .15 + min(datapc$v2)* 1,
#                            max(data$PC2) + max(data$PC2) * .15 + max(datapc$v2)* 1))
  
  pcplot <- ggplot(data, 
                   aes(x=PC1, y=PC2)) + 
    geom_hline(yintercept = 0, 
               size=.5,
               col = "grey70") + 
    geom_vline(xintercept = 0, 
               size=.5,
               col = "grey70") +
    geom_path(alpha = .2) +
    geom_segment(data=datapc,
                 aes(x=0, y=0, xend=v1, yend=v2), 
                 arrow=arrow(length=unit(0.1,"cm")), 
                 alpha=0.6, color="RED") +
    geom_text(alpha = .9, 
              size = 2, 
              aes(label = obsnames)) +
   
    geom_text(data = datapc, 
              aes(x = v1, y = v2, label = varnames),
              size = 2.75, 
              vjust = 1,
              alpha = 1,
              color="BLACK") +
   
    labs(list(colour = "Relative\n effect", 
              shape = "Ecosystem\n driver",
#               title = PC$ECO,
              x = paste0("PC1 (", round(summary(PC)$importance[2,1] * 100, 2), "%)"), 
              y = paste0("PC2 (", round(summary(PC)$importance[2,2] * 100, 2), "%)"))) +
  coord_equal(ratio = 1) +
#     coord_equal(ratio = 1, ylim = limits[2,], xlim = limits[1,]) +
#     scale_x_continuous(limits = limits[1,])+
#     scale_y_continuous(limits = limits[2,]) +
    theme_bw()
  pcplot
}