

# ------------------------- Calculate CV --------------------------------------

CV_peaks <- function(ff, kept_cells, channels) {
    ff_filtered <- ff[kept_cells, ]
    peaks <- PeacoQC(ff_filtered, determine_good_cells = FALSE, channels = channels)$AllPeaks
    cv <- apply(peaks, 2, function(x) {
        (stats::sd(x)/mean(x)) * 100
    })
}


# ----------------------------- Determine all peaks for one channel -----------

DetermineAllPeaks <- function(ff, channel, breaks){

  full_channel_peaks <- FindThemPeaks(ff@exprs[,channel])

  if (length(full_channel_peaks)<1){
    return(NA)
  }


  peak_results <- list()

  minimum <- min(ff@exprs[,channel])
  maximum <- max(ff@exprs[,channel])
  range <- abs(minimum) + abs(maximum)

  peaks <- lapply(breaks, function(x){FindThemPeaks(ff@exprs[x,channel])})

  # peaks <- lapply(splits, with, Peaks)

  names(peaks) <- seq_len(length(peaks))

  peak_frame <- plyr::ldply(peaks, cbind)

  colnames(peak_frame) <- c("Bin", "Peak")


  # Remove NA values (bins where no peaks were found)

  if (length(which(is.na(peak_frame$Peak))) > 0){

    peak_frame <- peak_frame[!is.na(peak_frame[,2]),]
    peaks <-peaks[!is.na(peaks)]

  }

  lengths <- sapply(peaks, length)


  nr_peaks <- unique(lengths)
  nr_peaks_in_bins <- tabulate(match(lengths, nr_peaks))


  total_full_length_peaks <- length(full_channel_peaks)


  most_occuring_peaks <- max(nr_peaks[which(tabulate(match(lengths, nr_peaks))
    > (0.2)*length(lengths))])

  ind_bins_nr_peaks <- lengths == most_occuring_peaks
  limited_nr_peaks <- do.call(rbind, peaks[ind_bins_nr_peaks])
  medians_to_use <- apply(limited_nr_peaks, 2, stats::median)



  # Trying to find the medians for the bins with the smaller amount of peaks
  #and using these as initial cluster centra kmeans

  if (length(medians_to_use) > 1) {

    # peak_clustering <- kmeans(peak_frame[,2], centers = peaks_medians)
    peak_clustering <- stats::kmeans(peak_frame[,2], centers = medians_to_use)

    cluster_medians <- peak_clustering$centers

    names(cluster_medians) <- rownames(peak_clustering$centers)
    peak_frame <- cbind(peak_frame, "Cluster"= factor(peak_clustering$cluster))

    final_medians <- c()

    to_use_clusters <- names(cluster_medians)[sapply(full_channel_peaks,
      function(x)which.min(abs(x-cluster_medians)))]

    peak_frame <- peak_frame[peak_frame$Cluster %in% to_use_clusters,]

    final_medians <- medians_to_use[which(names(cluster_medians) %in% to_use_clusters)]


  } else { # If only one peak was found per bin, no kmeans has to happen
    peak_frame <- cbind(peak_frame, "Cluster"= rep("1", dim(peak_frame)[1]))
    final_medians <- stats::median(peak_frame$Peak)

  }


  bin_list <- lapply(unique(peak_frame$Bin),
    DuplicatePeaks,
    peak_frame,
    final_medians)

  peak_frame <- do.call(rbind, bin_list)
  rownames(peak_frame) <- c(1:dim(peak_frame)[1])


  #Removes clusters that are too small to be used for isolation forest (< 1/3 of the data)

  peak_frame <- TooSmallClusters(peak_frame, peak_frame$Cluster)

  return(peak_frame)

}


FindThemPeaks <- function (channel_data)
{
  n <- which(!is.na(channel_data))
  if (length(n) < 3) {
    return(list(Peaks = NA, P.ind = 0, P.h = 0))
  }

  dens <- stats::density(channel_data[which(!is.na(channel_data))], adjust = 1)
  dens <- stats::smooth.spline(dens$x, dens$y, spar = 0.6)
  dens$y[which(dens$y < 0)] <- 0

  if (all(is.na(dens)))
    return(NA)


  peaks.ind <- c()
  for (i in 1:(length(dens$y) - 2)) {
    if (dens$y[i + 1] > dens$y[(i + 2)] && dens$y[i +
        1] > dens$y[i] && dens$y[i + 1] > 0.1 *
        max(dens$y)) {
       peaks.ind <- c(peaks.ind, i + 1)
    }
  }
  if (all(is.na(peaks.ind))) {
    # warning("No peaks could be found, returning the maximum value of density.")
    peaks.ind <- which.max(dens$y)
  }

  return(dens$x[peaks.ind])
}







  # # Eerste afgeleide
  # first_derivative <- paste(as.character(sign(diff(dens$y))), collapse = "")
  # first_derivative <- gsub("1", "+", gsub("-1", "-", first_derivative))
  # peak_attribute <- "[+]{1,}[-]{1,}"
  # possible_peaks <- gregexpr(peak_attribute, first_derivative)[[1]]
  #
  # if (possible_peaks[1] < 0)
  #   return(NA)
  # x1 <- possible_peaks
  # x2 <- possible_peaks + attr(possible_peaks, "match.length")
  # attributes(x1) <- NULL
  # attributes(x2) <- NULL
  # n <- length(x1)
  # peak_value <- peak_position <- numeric(n)
  # for (i in 1:n) {
  #   peak_position[i] <- which.max(dens$y[x1[i]:x2[i]]) + x1[i] - 1
  #   peak_value[i] <- dens$y[peak_position[i]]
  # }
  # minpeakheight <- max(dens$y) * 0.1
  #
  # inds <- which(peak_value >= minpeakheight & peak_value - pmax(dens$y[x1], dens$y[x2]) >= 0)
  # if (length(peak_value) == 0)
  #   return(NA)
  #
  # peaks <- channel_data[peak_position[inds]]
  #
  # return(peaks)
# }


# ------------------------------- Remove duplicate peaks ----------------------

DuplicatePeaks <- function(i, data, medians ){
  to_use <- data[data$Bin == i,]
  duplicates <- as.numeric(to_use$Cluster[ duplicated(to_use$Cluster)])
  for (duplex in duplicates){
    ind <- which(to_use$Cluster == duplex)

    diff <- abs(to_use$Peak[to_use$Cluster == duplex] - medians[duplex])
    to_use <- to_use[-ind[which.max(diff)],]
  }

  to_use <- to_use[order(to_use$Cluster),]


  return(to_use)
}

# ------------------------ Find the medians per cluster -----------------------

clusterMedian = function(i, data, clusters) {
  ind = (clusters == i)
  mediancluster <-  stats::median(data[ind])
  names(mediancluster) <- i
  return(mediancluster)
}


# ------------------------- Remove too small clusters -------------------------

TooSmallClusters <- function(data, clusters ){
  for (cluster in unique(clusters)){
    if (length(which(clusters == cluster)) < length(unique(data$Bin))*(2/3)){
      ind <-  which(clusters == cluster)
      data <- data[-ind,]

    }
  }
  return(data)
}


# ---------------------------- Function to remove outliers based on mad -----



MADOutliers <- function(peak, MAD) {


  kernel <- stats::ksmooth(seq_along(peak),
    peak,
    x.points = seq_along(peak),
    bandwidth = 50)

  median_peak <- stats::median(kernel$y, na.rm=TRUE)
  mad_peak <- stats::mad(kernel$y)

  upper_interval <- stats::median(median_peak, na.rm=TRUE)+MAD*(mad_peak)
  lower_interval <- stats::median(median_peak, na.rm=TRUE)-MAD*(mad_peak)

  outliers <- ifelse(kernel$y > upper_interval, TRUE,
    ifelse(kernel$y < lower_interval,
      TRUE, FALSE))
  return(outliers)}

MADOutlierMethod<- function(peak_frame, MAD) {
  to_remove_bins_df <- apply(peak_frame, 2, MADOutliers, MAD)

  return(to_remove_bins_df)
}


# -------------------------------- Binning ------------------------------------

SplitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)

  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

SplitWithOverlapMids <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)

  mids <- sapply(1:length(starts), function(i) starts[i] + ceiling(overlap/2))
  mids <- mids[mids < max(vec)]
  return(mids)
}


# ----------------------- Isolation Tree ------------------------------------



isolationTreeSD <- function(x, max_depth = as.integer(ceiling(log2(nrow(x)))),
  gain_limit = 0.55){

  res <- data.frame(id = 1,
    left_child = NA,
    right_child = NA,
    gain = NA,
    split_column = NA,
    split_value = NA,
    depth = 0,
    to_split = TRUE,
    path_length = NA )
  selection <- matrix(TRUE, ncol = nrow(x))


  nodes_to_split <- which(res$to_split)

  while(length(nodes_to_split) > 0){

    node <- nodes_to_split[1]

    rows <- which(selection[node, ])

    if(length(rows) > 3 & res[node, "depth"] < max_depth){
      gain_max <- gain_limit
      split_value <- NA
      split_col <- NA

      col <- 1
      for(col in seq_len(ncol(x))){
        x_col <- sort(x[rows, col])
        base_sd <- stats::sd(x_col)
        gain_max_col <- 0
        split_v <- NA
        mean_val <- c()
        for(i in c(1:(length(x_col)-1))){

          sd_1 <- stats::sd(x_col[c(1:i)])
          sd_2 <- stats::sd(x_col[c((i+1):length(x_col))])
          mean_val <- c(mean_val, mean(x_col[c(1:i)]))


          if (i == 1){
            sd_1 <- 0
          } else if (i == length(x_col) - 1){
            sd_2 <- 0
          }
          gain <- (base_sd - mean(c(sd_1,sd_2)))/base_sd

          if (is.na(gain) == TRUE){
            next()
          }
          if(gain >= gain_max_col){
            gain_max_col <- gain
            split_v <- x_col[i]
          }
        }
        if(gain_max_col > gain_max){
          gain_max <- gain_max_col
          split_value <- split_v
          split_col <- col
        }
      }

      if(!is.na(split_value)) {
        left_node <- selection[node, ] & x[ ,split_col] <= split_value
        right_node <- selection[node, ] & x[ , split_col] > split_value

        if(length(which(left_node)) == length(rows) ||
            length(which(right_node)) == length(rows)){
          res[node, "to_split"] <- FALSE
          nodes_to_split <- which(res$to_split==TRUE)
          n_datapoints <- length(which(selection[,node]))
          res[node, "path_length"] <- avgPL(n_datapoints) + res[node, "depth"]


          next()
        }

        max_node <- nrow(res)
        left_id <- max_node + 1
        right_id <- max_node + 2
        res[node, c("left_child", "right_child", "gain",
          "split_column","split_value","to_split")] <- c(
            left_id, right_id, gain_max,
            colnames(x)[split_col], split_value, FALSE)
        res <- rbind(res,
          data.frame(id = c(left_id, right_id),
            left_child = NA,
            right_child = NA,
            gain = NA,
            split_column = NA,
            split_value = NA,
            depth = res[node,"depth"]+1,
            to_split = TRUE,
            path_length = NA))


        selection <- rbind(selection, left_node, right_node)




      } else {
        res[node, "to_split"] <- FALSE
        n_datapoints <- length(which(selection[node,]))
        res[node, "path_length"] <- avgPL(n_datapoints) + res[node, "depth"]


      }

    } else {
      res[node, "to_split"] <- FALSE
      n_datapoints <- length(which(selection[node,]))
      res[node, "path_length"] <- avgPL(n_datapoints) + res[node, "depth"]


    }

    nodes_to_split <- which(res$to_split==TRUE)

  }

  res$n_datapoints = rowSums(selection)
  res$anomaly_score = 2^(-(res$path_length)/(avgPL(sum(res$n_datapoints))))

#
#   if (plot == TRUE){
#     relations <- data.frame(from = c(res$id, res$id),
#       to = c(res$left_child, res$right_child))
#
#     relations <- stats::na.omit(relations)
#
#     g_1 <- graph_from_data_frame(relations, directed = TRUE, vertices = res)
#     plot.igraph(g_1,
#       layout = layout_as_tree(g_1),
#       vertex.label = rowSums(selection))
#   }


  return(list(res = res,
    selection = selection))
}



# ----------------------------- avgPL -----------------------------------------

avgPL <- function(n_datapoints){

  if (n_datapoints -1 == 0){
    AVG <- 0
  } else {
    AVG <- 2*(log(n_datapoints - 1) +  0.5772156649) -
      (2*(n_datapoints -1))/(n_datapoints)
  }

  return (AVG)
}


# --------------------------- append column to ff -----------------------------

Append_GoodCells <- function(ff, bad_cells){

  pd <- flowWorkspace::pData(flowCore::parameters(ff))
  new_pd_name <-  ncol(ff@exprs) +1
  new_pd_name <- paste0("$P", new_pd_name)
  new_pd <- data.frame(name = "GoodCells", desc = "PeacoQC_Cells", range = 1,
    minRange = 0, maxRange = 1)
  rownames(new_pd) <- new_pd_name
  pd <- rbind(pd, new_pd)
  ff@exprs <- cbind(ff@exprs, "GoodCells" = as.numeric(bad_cells))
  flowWorkspace::pData(flowCore::parameters(ff)) <- pd
  ff

}

