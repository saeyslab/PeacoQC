channel_data <- ff@exprs[,channel]
channel_data <- channel_data[rep(1:5000,1000)]

microbenchmark::microbenchmark(FindThemPeaks(channel_data), times = 10)

profvis::profvis({FindThemPeaks(channel_data)})



microbenchmark::microbenchmark(
  {peaks <- lapply(breaks, function(x){FindThemPeaks(ff@exprs[x,channel])})},
  {e <- ff@exprs[,channel]; peaks <- lapply(breaks, function(x){FindThemPeaks(e[x])})},
  times = 10)


install.packages("snow")
library(snow)
z=vector('list',4)
z=1:12
system.time({peaks <- lapply(breaks, function(x){FindThemPeaks(ff@exprs[x,channel])})})

FindThemPeaks_NonPar <-  function(ff, channel, breaks){
  peaks <- lapply(breaks, function(x) {FindThemPeaks(ff@exprs[x,channel])})
}

FindThemPeaks_Par <- function(ff, channel, breaks){
  e <- ff@exprs[,channel]
  cl<-makeCluster(4, type="SOCK")
  clusterExport(cl,"e")
  peaks <- clusterApply(cl, breaks, function(x) {PeacoQC:::FindThemPeaks(e[x])})
  stopCluster(cl)
}

microbenchmark::microbenchmark(
  {FindThemPeaks_NonPar(ff, channel, breaks)},
  {FindThemPeaks_Par(ff, channel, breaks)},
  times = 10)



profvis::profvis({DetermineAllPeaks(ff, channel, breaks)})


test1 <- DetermineAllPeaks(ff, channel, breaks)
test2 <- DetermineAllPeaks(ff, channel, breaks)

FindThemPeaks <- function (channel_data)
{
  n <- which(!is.na(channel_data))
  if (length(n) < 3) {
    return(NA)
  }

  dens <- stats::density(channel_data[which(!is.na(channel_data))], adjust = 1)
  dens <- stats::smooth.spline(dens$x, dens$y, spar = 0.6)
  dens$y[which(dens$y < 0)] <- 0

  if (all(is.na(dens)))
    return(NA)

  # w <- 1
  # peaks <- c()
  # for (i in 1:(length(dens$y) - 2)) {
  #   if (dens$y[i + w] > dens$y[(i + w + 1):(i + 2 * w)] && dens$y[i +
  #       w] > dens$y[i:(i + w - 1)] && dens$y[i + w] > 1/3 *
  #       max(dens$y)) {
  #     peaks <- c(peaks, dens$x[i + w])
  #   }
  # }
  n <- length(dens$y)
  selection <- (dens$y[2:(n-1)] > dens$y[1:(n-2)]) &
    (dens$y[2:(n-1)] > dens$y[3:n] ) &
    (dens$y[2:(n-1)] > (1/3 * max(dens$y)))
  peaks <- dens$x[-1][selection]

  if (all(is.na(peaks))) {
    # warning("No peaks could be found, returning the maximum value of density.")
    peaks <- dens$x[which.max(dens$y)]
  }
  if (any(peaks < 0) ){
    return(NA)
  }

  return(peaks)
}





microbenchmark::microbenchmark(
  DetermineAllPeaks(ff, channel, breaks),
  DetermineAllPeaks2(ff, channel, breaks), times = 10)
