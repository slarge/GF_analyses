whiten <- function (dens, lambda = 0.9) 
{
  dens$y <- lambda * dens$y + (1 - lambda)/diff(range(dens$x))
  dens
}


spRatio <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topright", 
          mfrow = n2mfrow(length(imp.vars)), 
#           mfrow = length(imp.vars),
          mar = c(4.5, 1.5, 0.5,4.5), 
          omi = c(0.1, 0.25, 0.1, 0.1), bin = F, nbin = 101, 
          leg.panel = 1, barwidth = 1, cex.legend = 0.8, line.ylab = 1.5, 
          ...) 
{
  if (is.null(imp.vars)) 
    imp.vars <- imp.var.names <- names(sort(rowMeans(obj$imp.rsq), 
                                            decreasing = T))[1:2]
  is.binned <- function(obj) {
    compact <- obj$call$compact
    if (is.null(compact)) 
      FALSE
    else eval(compact)
  }
  normalize.histogram <- function(ci, integral = 1, bin = F, 
                                  nbin = 101) {
    if (bin) {
      brks <- seq(min(ci$x), max(ci$x), len = nbin)
      xx <- cut(ci$x, breaks = brks, inc = T)
      yy <- tapply(ci$y, xx, sum)
      yy[is.na(yy)] <- 0
      ci <- list(x = 0.5 * (brks[-1] + brks[-nbin]), y = yy)
    }
    dx <- min(diff(ci$x))
    Id <- sum(ci$y * dx)
    ci$y <- ci$y/Id * integral
    ci
  }
  normalize.density <- function(d, integral = 1, integrate = T) {
    Id <- if (integrate) 
      integrate.density(d)
    else 1
    d$y <- d$y/Id * integral
    d
  }
  integrate.density <- function(d) {
    integrate(approxfun(d, rule = 2), lower = min(d$x), 
              upper = max(d$x))$value
  }
  scale.density <- function(d, scale = 1/mean(d$y)) {
    d$y <- d$y * scale
    d
  }
  par(mfrow = mfrow)
  nice.names <- structure(as.list(imp.vars.names), names = imp.vars)
  for (i in imp.vars) {
    imp <- importance(obj)[i]
    resA <- obj$res[obj$res$var == i, ]
    splits <- resA$split
    w <- pmax(resA$improve.norm, 0)
    X <- na.omit(obj$X[, i])
    rX <- range(X)
    dX <- diff(rX)
    dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                    to = rX[2])
    if ((dX/dImp$bw) > 50) 
      dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                      to = rX[2], bw = dX/50)
    dImpNorm <- normalize.density(dImp, imp, integrate = T)
    dObs <- density(X, from = rX[1], to = rX[2])
    if ((dX/dObs$bw) > 50) 
      dObs <- density(X, from = rX[1], to = rX[2], bw = dX/50)
    dObs <- whiten(dObs, lambda = 0.9)
    dObsNorm <- normalize.density(dObs, imp, integrate = T)
    ci <- cumimp(obj, i, standardize = F)
    ci$y <- diff(c(0, ci$y))
    ci <- normalize.histogram(ci, imp, bin = bin | !is.binned(obj), 
                              nbin = nbin)
    dStd <- dImp
    dStd$y <- dImp$y/dObs$y
    dStdNorm <- try(normalize.density(dStd, imp, integrate = T))
    if (class(dStdNorm) == "try-error") 
      dStdNorm <- normalize.histogram(dStd, imp)
    plot(ci, type = "h", col = "grey60", xlim = range(splits), 
         lwd = barwidth, ylim = c(0, max(dImpNorm$y, dObsNorm$y, 
                                         dStdNorm$y) * 1.1), lend = 2, xlab = nice.names[[i]], 
         ylab = "", ...)
#     lines(dImpNorm, col = "black", lwd = 2)
#     lines(dObsNorm, col = "red", lwd = 2)
    lines(dStdNorm, col = "blue", lwd = 2)
    abline(h = mean(dStdNorm$y)/mean(dStd$y), lty = 2, col = "blue")
#     if (i == imp.vars[leg.panel]) 
#       legend(leg.posn, legend = c("Density of splits", 
#                                   "Density of data", "Ratio of densities", "Ratio=1"), 
#              lty = c(1, 1, 1, 2), col = c("black", "red", 
#                                           "blue", "blue"), cex = cex.legend, bty = "n", 
#              lwd = 1)
  }
  mtext("Ratio of Densities", side = 2, line = line.ylab, outer = T)
}


spImportance <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topright", 
                     mfrow = n2mfrow(length(imp.vars)), mar = c(4.5, 1.5, 0.5, 
                                                                4.5), omi = c(0.1, 0.25, 0.1, 0.1), bin = F, nbin = 101, 
                     leg.panel = 1, barwidth = 1, cex.legend = 0.8, line.ylab = 1.5, 
                     ...) 
{
  if (is.null(imp.vars)) 
    imp.vars <- imp.var.names <- names(sort(rowMeans(obj$imp.rsq), 
                                            decreasing = T))[1:2]
  is.binned <- function(obj) {
    compact <- obj$call$compact
    if (is.null(compact)) 
      FALSE
    else eval(compact)
  }
  normalize.histogram <- function(ci, integral = 1, bin = F, 
                                  nbin = 101) {
    if (bin) {
      brks <- seq(min(ci$x), max(ci$x), len = nbin)
      xx <- cut(ci$x, breaks = brks, inc = T)
      yy <- tapply(ci$y, xx, sum)
      yy[is.na(yy)] <- 0
      ci <- list(x = 0.5 * (brks[-1] + brks[-nbin]), y = yy)
    }
    dx <- min(diff(ci$x))
    Id <- sum(ci$y * dx)
    ci$y <- ci$y/Id * integral
    ci
  }
  normalize.density <- function(d, integral = 1, integrate = T) {
    Id <- if (integrate) 
      integrate.density(d)
    else 1
    d$y <- d$y/Id * integral
    d
  }
  integrate.density <- function(d) {
    integrate(approxfun(d, rule = 2), lower = min(d$x), 
              upper = max(d$x))$value
  }
  scale.density <- function(d, scale = 1/mean(d$y)) {
    d$y <- d$y * scale
    d
  }
  par(mfrow = mfrow)
  nice.names <- structure(as.list(imp.vars.names), names = imp.vars)
  for (i in imp.vars) {
    imp <- importance(obj)[i]
    resA <- obj$res[obj$res$var == i, ]
    splits <- resA$split
    w <- pmax(resA$improve.norm, 0)
    X <- na.omit(obj$X[, i])
    rX <- range(X)
    dX <- diff(rX)
    dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                    to = rX[2])
    if ((dX/dImp$bw) > 50) 
      dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                      to = rX[2], bw = dX/50)
    dImpNorm <- normalize.density(dImp, imp, integrate = T)
    dObs <- density(X, from = rX[1], to = rX[2])
    if ((dX/dObs$bw) > 50) 
      dObs <- density(X, from = rX[1], to = rX[2], bw = dX/50)
    dObs <- whiten(dObs, lambda = 0.9)
    dObsNorm <- normalize.density(dObs, imp, integrate = T)
    ci <- cumimp(obj, i, standardize = F)
    ci$y <- diff(c(0, ci$y))
    ci <- normalize.histogram(ci, imp, bin = bin | !is.binned(obj), 
                              nbin = nbin)
    dStd <- dImp
    dStd$y <- dImp$y/dObs$y
    dStdNorm <- try(normalize.density(dStd, imp, integrate = T))
    if (class(dStdNorm) == "try-error") 
      dStdNorm <- normalize.histogram(dStd, imp)
    plot(ci, type = "h", col = "grey60", xlim = range(splits), 
         lwd = barwidth, ylim = c(0, max(dImpNorm$y, dObsNorm$y, 
                                         dStdNorm$y) * 1.1), lend = 2, xlab = nice.names[[i]], 
         ylab = "", ...)
        lines(dImpNorm, col = "black", lwd = 2)
    #     lines(dObsNorm, col = "red", lwd = 2)
#     lines(dStdNorm, col = "blue", lwd = 2)
#     abline(h = mean(dStdNorm$y)/mean(dStd$y), lty = 2, col = "blue")
    #     if (i == imp.vars[leg.panel]) 
    #       legend(leg.posn, legend = c("Density of splits", 
    #                                   "Density of data", "Ratio of densities", "Ratio=1"), 
    #              lty = c(1, 1, 1, 2), col = c("black", "red", 
    #                                           "blue", "blue"), cex = cex.legend, bty = "n", 
    #              lwd = 1)
  }
  mtext("Density of Splits", side = 2, line = line.ylab, outer = T)
}


spData <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topright", 
                          mfrow = n2mfrow(length(imp.vars)), mar = c(4.5, 1.5, 0.5, 
                                                                     4.5), omi = c(0.1, 0.25, 0.1, 0.1), bin = F, nbin = 101, 
                          leg.panel = 1, barwidth = 1, cex.legend = 0.8, line.ylab = 1.5, 
                          ...) 
{
  if (is.null(imp.vars)) 
    imp.vars <- imp.var.names <- names(sort(rowMeans(obj$imp.rsq), 
                                            decreasing = T))[1:2]
  is.binned <- function(obj) {
    compact <- obj$call$compact
    if (is.null(compact)) 
      FALSE
    else eval(compact)
  }
  normalize.histogram <- function(ci, integral = 1, bin = F, 
                                  nbin = 101) {
    if (bin) {
      brks <- seq(min(ci$x), max(ci$x), len = nbin)
      xx <- cut(ci$x, breaks = brks, inc = T)
      yy <- tapply(ci$y, xx, sum)
      yy[is.na(yy)] <- 0
      ci <- list(x = 0.5 * (brks[-1] + brks[-nbin]), y = yy)
    }
    dx <- min(diff(ci$x))
    Id <- sum(ci$y * dx)
    ci$y <- ci$y/Id * integral
    ci
  }
  normalize.density <- function(d, integral = 1, integrate = T) {
    Id <- if (integrate) 
      integrate.density(d)
    else 1
    d$y <- d$y/Id * integral
    d
  }
  integrate.density <- function(d) {
    integrate(approxfun(d, rule = 2), lower = min(d$x), 
              upper = max(d$x))$value
  }
  scale.density <- function(d, scale = 1/mean(d$y)) {
    d$y <- d$y * scale
    d
  }
  par(mfrow = mfrow)
  nice.names <- structure(as.list(imp.vars.names), names = imp.vars)
  for (i in imp.vars) {
    imp <- importance(obj)[i]
    resA <- obj$res[obj$res$var == i, ]
    splits <- resA$split
    w <- pmax(resA$improve.norm, 0)
    X <- na.omit(obj$X[, i])
    rX <- range(X)
    dX <- diff(rX)
    dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                    to = rX[2])
    if ((dX/dImp$bw) > 50) 
      dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                      to = rX[2], bw = dX/50)
    dImpNorm <- normalize.density(dImp, imp, integrate = T)
    dObs <- density(X, from = rX[1], to = rX[2])
    if ((dX/dObs$bw) > 50) 
      dObs <- density(X, from = rX[1], to = rX[2], bw = dX/50)
    dObs <- whiten(dObs, lambda = 0.9)
    dObsNorm <- normalize.density(dObs, imp, integrate = T)
    ci <- cumimp(obj, i, standardize = F)
    ci$y <- diff(c(0, ci$y))
    ci <- normalize.histogram(ci, imp, bin = bin | !is.binned(obj), 
                              nbin = nbin)
    dStd <- dImp
    dStd$y <- dImp$y/dObs$y
    dStdNorm <- try(normalize.density(dStd, imp, integrate = T))
    if (class(dStdNorm) == "try-error") 
      dStdNorm <- normalize.histogram(dStd, imp)
    plot(ci, type = "h", col = "grey60", xlim = range(splits), 
         lwd = barwidth, ylim = c(0, max(dImpNorm$y, dObsNorm$y, 
                                         dStdNorm$y) * 1.1), lend = 2, xlab = nice.names[[i]], 
         ylab = "", ...)
    #     lines(dImpNorm, col = "black", lwd = 2)
        lines(dObsNorm, col = "red", lwd = 2)
    #     lines(dStdNorm, col = "blue", lwd = 2)
    #     abline(h = mean(dStdNorm$y)/mean(dStd$y), lty = 2, col = "blue")
    #     if (i == imp.vars[leg.panel]) 
    #       legend(leg.posn, legend = c("Density of splits", 
    #                                   "Density of data", "Ratio of densities", "Ratio=1"), 
    #              lty = c(1, 1, 1, 2), col = c("black", "red", 
    #                                           "blue", "blue"), cex = cex.legend, bty = "n", 
    #              lwd = 1)
  }
  mtext("Density of Data", side = 2, line = line.ylab, outer = T)
}


spCumPlot <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topleft", 
          leg.nspecies = 10, legend = TRUE, mfrow = rev(n2mfrow(length(imp.vars) * 
                                                                  (show.species + show.overall))), show.species = TRUE, 
          show.overall = TRUE, mar = c(0, 2.1, 1.1, 0), omi = c(0.75, 
                                                                0.75, 0.1, 0.1), common.scale = F, line.ylab = 1, cex.legend = 0.75, 
          ...) 
{
  if (is.null(imp.vars)) 
    imp.vars <- imp.var.names <- names(importance(obj))[1:2]
  par(mfrow = mfrow)
  require(RColorBrewer)

  plotPalette <- colorRampPalette(brewer.pal(length(levels(obj$res.u$spec)), "Dark2"))
  cols <- plotPalette(length(levels(obj$res.u$spec)))

#   cols <- rainbow(length(levels(obj$res.u$spec)))
  names(cols) <- levels(obj$res.u$spec)
  xaxt <- if (show.overall) 
    "n"
  else "s"
  if (show.species) {
    for (varX in imp.vars) {
      CU <- cumimp(obj, varX, "Species")
      xlim <- range(sapply(CU, "[[", "x"), na.rm = T)
      ylim <- range(sapply(CU, "[[", "y"), na.rm = T)
      plot(xlim, ylim, type = "n", xlab = if (show.overall) 
        ""
        else imp.vars.names[imp.vars == varX], ylab = "", 
        xaxt = xaxt, ...)
      for (species in names(CU)) {
        isub <- seq(1, length(CU[[species]]$x), len = pmin(500, 
                                                           length(CU[[species]]$x)))
        lines(CU[[species]]$x[isub], CU[[species]]$y[isub], 
              type = "s", col = cols[species], lwd = 1.5)
      }
      no.species <- length(names(cols))
      imp.sp <- sapply(CU, function(cu) max(cu$y))
      best <- order(-imp.sp)[1:min(leg.nspecies, length(imp.sp))]
      if (legend) 
        legend(x = leg.posn, legend = names(cols)[best], 
               pch = rep(16, no.species)[best], col = cols[best], 
               bty = "n", cex = cex.legend, pt.lwd = 2)
    }
  }
  if (show.overall) {
    for (varX in imp.vars) {
      CU <- cumimp(obj, varX)
      ymax <- max(CU$y)
      if (varX == imp.vars[1]) 
        ymax1 <- ymax
      isub <- seq(1, length(CU$x), len = pmin(500, length(CU$x)))
      plot(CU$x[isub], CU$y[isub], type = "s", ylab = "", 
           xlab = imp.vars.names[imp.vars == varX], ylim = c(0, 
                                                             if (common.scale) ymax1 else ymax), ...)
    }
  }
  mtext("Cumulative importance", side = 2, line = line.ylab, 
        outer = TRUE)
}

impPlot <- function (obj, ..., las = 1, cex.axis = 0.7, cex.names = cex.axis, 
          horiz = TRUE) 
  {
#     imp.a <- importance(obj, "Accuracy")
    imp.w <- importance(obj, "Weighted")
#     o.a <- order(imp.a)
    o.w <- order(imp.w)
#     barplot(imp.a[o.a], names = names(imp.a[o.a]), horiz = horiz, 
#             main = "Accuracy importance", las = las, cex.axis = cex.axis, 
#             cex.names = cex.names, ...)
    barplot(imp.w[o.w], names = names(imp.w[o.w]), horiz = horiz, 
            las = 1, 
#             main = expression(paste(R^2, " weighted importance")),
            main = "",
            las = las, cex.axis = cex.axis, cex.names = cex.names, 
            ...)
  }

       
combinedCumPlot <- function (obj, weight = c("uniform", "species", "rsq.total", 
                                             "rsq.mean", "site", "site.species", "site.rsq.total", "site.rsq.mean")[3], 
                              use.diff = FALSE, prednames = names(obj$X)[-1], show.weights = FALSE, 
                              show.gears = TRUE, sort = TRUE, ...) 
#   obj <- f2
#     weight = "rsq.total"
{
  if ((nw <- length(weight)) > 1) 
    show.gears <- show.weights <- FALSE
  gearnames <- names(obj$dens)[-1]
  CU <- if (!show.gears) 
    list()
  else lapply(prednames, function(predictor) do.call("rbind", 
                                                     lapply(obj$gears[[predictor]], function(gear) {
                                                       cu <- cumimp(obj, predictor, gear = gear)
                                                       dens <- density(obj, predictor, gear = gear, gridded = T)
                                                       if (use.diff) 
                                                         cu$y <- diff(c(0, cu$y))
                                                       data.frame(predictor = predictor, gear = gear, value = cu$x, 
                                                                  CU = cu$y, dens = dens$y)
                                                     })))
  CU <- c(CU, lapply(prednames, function(predictor) do.call("rbind", 
                                                            lapply(weight, function(w) {
                                                              cu <- cumimp(obj, predictor, weight = w)
                                                              if (use.diff) cu$y <- diff(c(0, cu$y))
                                                              data.frame(predictor = predictor, gear = paste("combined", 
                                                                                                             w, sep = "."), value = cu$x, CU = cu$y, dens = -1)
                                                            }))))
  CU <- do.call("rbind", CU)
  imp <- importance(obj)
  o <- order(-imp)
  CU$predictor <- ordered(CU$predictor, levels = if (sort) 
    names(sort(-imp))
    else prednames)

#   CU$predictor <- as.character(CU$predictor)
#   CU$predictor[CU$predictor == "EXP" &
#      CU$gear %in% c("ESS", "WSS")] <- "EXP_CA"
#   CU$predictor <- as.factor(CU$predictor)
#   
  
  colPal <- c("GB" = "#E41A1C", "GOM" = "#377EB8", "MAB" = "#4DAF4A", "SS" = "#984EA3", 
              "ESS" = "#FF7F00", "WSS" = "#A65628", "combined.rsq.total" =  "#000000")
  sizPal <- c(rep(1,6), 1.75)
  alphPal <- c(rep(1,6), .5)
  names(sizPal) <- names(colPal)  

  cuPlot <- ggplot(CU, aes(x = value, y = CU)) +
    geom_line(aes(color = gear, size = gear, alpha = gear)) +
    scale_color_manual(values = colPal) +
    scale_size_manual(values = sizPal, guide = "none") +
    scale_alpha_manual(values = alphPal, guide = "none") +
    facet_wrap(~predictor, scales = "free_x")+
    labs(list(color = "Ecosystem",
              x = "Predictor value", 
              y = "Cumulative Importance")) +
    theme_bw()
  cuPlot
}


shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else
    out<-x
  out
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}
