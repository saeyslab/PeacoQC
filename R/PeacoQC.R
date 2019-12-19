

PeacoQC <- function(ff,
                    channels,
                    # You could set this to "MAD" or "IT" you want to only use that outlier detection method
                    determine_good_cells = TRUE,
                    plot = FALSE,
                    output_directory = NULL,
                    MAD = 6,
                    events_per_bin = 2000,
                    IT_limit = 0.55,
                    consecutive_bins = 5,
                    plot_tree = F
){

  # Searching for the name of the ff
  if (length(ff@description$FILENAME)>0){
    message(paste0("Starting quality control analysis for ",
                   basename(ff@description$FILENAME)))
  }

  # Timing everything
  start_time <- Sys.time()

  # If no channels are given we will use all channels except for the time
  if (missing(channels)){
    time_index <- grep("Time",colnames(ff@exprs))
    channels <- colnames(ff@exprs)[-time_index]

  } else if (is.numeric(channels)){ #If numbers are given we select the correct channel names
    channels <- colnames(ff@exprs)[channels]
  }

  results <- list()

  ff_orig <- ff

  # Split the ff up in bins (overlapping)
  different_breaks <- splitWithOverlap(c(1:nrow(ff)),
                                       events_per_bin,
                                       ceiling(events_per_bin/2))

  names(different_breaks) <- seq_along(different_breaks)


  # If not enough bins are made, at least 100 should be present
  if (length(different_breaks) < 100){
    events_per_bin <- ceiling(nrow(ff)/100) *2
    different_breaks <- splitWithOverlap(c(1:nrow(ff)),
                                         events_per_bin,
                                         ceiling(events_per_bin/2))

    names(different_breaks) <- seq_along(different_breaks)

  }

  # This is purely for plotting purposes (Change to plotPeacoQC function?)
  mid_breaks_channel <- splitWithOverlapMids(c(1:nrow(ff)),
                                             events_per_bin,
                                             ceiling(events_per_bin/2))
  m <-  nrow(ff@exprs)%%events_per_bin


  # if not zero, you have something to add to your sequence
  if(m != 0) {

    mid_breaks_channel = c(mid_breaks_channel, nrow(ff@exprs))
  }



  results$Breaks <- different_breaks

  results$MidBreaks <- mid_breaks_channel


  # Check if there is an increasing or decreasing trent in the channels
  weird_channel_increasing <- NULL
  weird_channel_decreasing <- NULL


  for (channel in channels){

    channel_medians <- sapply(different_breaks,
                              function(x){median(ff@exprs[x,channel])})

    smoothed <- ksmooth(seq_along(channel_medians),
                        channel_medians,
                        x.points = seq_along(channel_medians),
                        bandwidth = 50)

    increasing <- cummax(smoothed$y) == smoothed$y
    decreasing <- cummin(smoothed$y) == smoothed$y

    if (length(which(increasing)) > (1/2)*length(increasing)){
      warning(paste0("There seems to be an increasing trent in channel ",
                     FlowSOM::get_markers(ff, channel),
                     "for file ",basename(ff@description$FILENAME),
                     ". Please inspect this before doing any further analysis"))
      plot(channel_medians, main = FlowSOM::get_markers(ff, channel))
      weird_channel_increasing <- c(weird_channel_increasing, channel)
    } else if (length(which(decreasing)) > (1/2)*length(decreasing)){
      warning(paste0("There seems to be a decreasing trent in channel ",
                     FlowSOM::get_markers(ff, channel),
                     "for file ",basename(ff@description$FILENAME),
                     ". Please inspect this before doing any further analysis"))
      plot(channel_medians, main = FlowSOM::get_markers(ff, channel))
      weird_channel_decreasing <- c(weird_channel_decreasing, channel)
    }
  }

  results$WeirdChannels <- list("Increasing"= weird_channel_increasing,
                                "Decreasing" = weird_channel_decreasing)

  peak_values_bin  <- list()


  # ------ Determine all peaks found in all channels and store them ---------
  i <- 0
  message("Calculating peaks")
  pb <-  utils::txtProgressBar(min = 0, max = length(channels), initial = 0,
                               char = "+", style = 3, width = 50)
  setTxtProgressBar(pb,i)

  channel <- channels[2]

  for (channel in channels){
    i <- i +1
    peak_frame <- DetermineAllPeaks(ff, channel, different_breaks)
    if (class(peak_frame) == "logical"){
      setTxtProgressBar(pb,i)
      next()
    }

    peak_values <- list()

    for (peak in unique(peak_frame$Cluster)){
      peak_ind <- peak_frame$Cluster == peak
      peak_data <- peak_frame[peak_ind,]
      peak_vector <- rep(median(peak_data$Peak),
                         length(different_breaks))

      peak_vector[as.numeric(peak_data$Bin)] <- peak_data$Peak

      names(peak_vector) <- seq_len(length(peak_vector))

      peak_values[[peak]] <- peak_vector

    }

    peak_values_bin[[channel]] <- peak_values

    results[[channel]] <- peak_frame

    setTxtProgressBar(pb,i)
  }

  close(pb)

  peak_values <- unlist(peak_values_bin, recursive = FALSE)

  all_peaks <- do.call(cbind, peak_values)

  all_peaks <- all_peaks %>% as.data.frame()

  rownames(all_peaks) <- names(peak_values[[1]])


  outlier_bins <- rep(TRUE, nrow(all_peaks))



  # ------------------------ Isolation Tree  ----------------------------------

  if (determine_good_cells == TRUE || determine_good_cells == "IT"){

    if (length(which(outlier_bins)) < 3){
      stop(paste0("There is an issue with file ",basename(ff@description$FILENAME),
                  ". There are no good cells left to perform quality control."))
    }

    data_for_trees <- all_peaks[outlier_bins,]

    tree <- isolationTreeSD(data_for_trees, gain_limit = IT_limit, plot = plot_tree)
    results$IT <- tree

    scores_to_use <- na.omit(tree$res[,c("n_datapoints", "anomaly_score")])
    node_to_keep <- rownames(scores_to_use)[which.max(scores_to_use$n_datapoints)]

    bins_isolated <- tree$selection[as.numeric(node_to_keep),]
    names(bins_isolated) <- which(outlier_bins == TRUE)

    outlier_bins[outlier_bins] <- bins_isolated

    removed_cells_IT <- different_breaks[names(bins_isolated)
                                         [which(bins_isolated == FALSE)]] %>%
      unlist
    removed_cells_IT <- removed_cells_IT[!duplicated(removed_cells_IT)]

    perc_IT <- length(removed_cells_IT)/nrow(ff_orig)
    perc_IT <- perc_IT * 100

    message(paste0("IT analysis removed ",
                   paste0(round(perc_IT,2),
                          " % of the measurements" )))

    results$ITPercentage <- perc_IT

  }


  # ------------------------ Outliers based on mad --------------------------

  # Bins that were selected based on their mad


  if (determine_good_cells == TRUE || determine_good_cells == "MAD"){

    if (length(which(outlier_bins)) < 3){
      stop(paste0("There is an issue with file ",basename(ff@description$FILENAME),
                  ". There are no good cells left to perform quality control."))
    }

    perc_not_outliers <- length(which(outlier_bins == TRUE))/length(outlier_bins)


    outlier_bins_df <- Flag.Outliers(all_peaks[outlier_bins,],MAD)

    outlier_bins_MAD <- apply(outlier_bins_df, 1, any)

    names(outlier_bins_MAD) <- which(outlier_bins == TRUE)

    contribution_MAD <- apply(outlier_bins_df, 2,
                              function(x){length(which(x == TRUE))/
                                  (length(x)*perc_not_outliers)})

    contribution_MAD <- round(contribution_MAD*100,2)

    mad_cells <- rep(TRUE, nrow(ff_orig))

    removed_cells_MAD <- different_breaks[names(outlier_bins_MAD)
                                          [which(outlier_bins_MAD == TRUE)]] %>%
      unlist
    removed_cells_MAD <- removed_cells_MAD[!duplicated(removed_cells_MAD)]

    perc_outliers <- length(removed_cells_MAD)/nrow(ff_orig)

    perc_outliers <- perc_outliers * 100

    message(paste0("MAD analysis removed ",
                   paste0(round(perc_outliers,2),
                          " % of the measurements" )))

    results$ContributionMad <- contribution_MAD

    results$MADPercentage <- perc_outliers

    mad_cells[removed_cells_MAD] <- FALSE

    results$OutlierMads <- mad_cells

    outlier_bins[outlier_bins] <- !outlier_bins_MAD

  }



  # ------------------------- indicate bad cells ----------------------------


  if (determine_good_cells %in% c(TRUE, "IT", "MAD")){

    # Determine which cells should also be removed because they fall inbetween FALSE regions
    outlier_bins_new <- inverse.rle(within.list(rle(outlier_bins),
                                                values[lengths<consecutive_bins]
                                                <- FALSE))

    consecutive_cells <- rep(TRUE, nrow(ff_orig))

    consecutive_bins_to_remove <- !(outlier_bins_new == outlier_bins)

    names(consecutive_bins_to_remove) <- seq_along(outlier_bins_new)

    removed_cells_consec <- different_breaks[names(consecutive_bins_to_remove)
                                             [which(consecutive_bins_to_remove == TRUE)]] %>%
      unlist
    removed_cells_consec <- removed_cells_consec[!duplicated(removed_cells_consec)]


    consecutive_cells[removed_cells_consec] <- FALSE

    results$ConsecutiveCells <- consecutive_cells


    outlier_bins <- outlier_bins_new

    to_remove_bins <- which(outlier_bins == FALSE)


    removed_cells <- different_breaks[to_remove_bins] %>% unlist

    removed_cells <- removed_cells[!duplicated(removed_cells)]


    # Summary of entire analysis
    bad_cells <- rep(TRUE, nrow(ff_orig))

    bad_cells[removed_cells] <- FALSE

    ff_orig <- Append_GoodCells(ff_orig,  bad_cells)


    percentage_removed <- (length(which(bad_cells == FALSE))/length(bad_cells))*100

    message(paste0("The algorithm removed ",
                   round(percentage_removed,2),
                   "% of the measurements"))

    results$remove_bins <- to_remove_bins
    results$GoodCells <- bad_cells
    results$PercentageRemoved <- percentage_removed
    if(PercentageRemoved > 70){
      warning(paste0("There was ", round(percentage_removed,3),
                     "% of the measurements removed for file ",
                     basename(ff@description$FILENAME)))
    }
  }

  # ----------------------- Final results ------------------------------------

  results$AllPeaks <- all_peaks
  results$FinalFF <- ff_orig
  end_time <- Sys.time()

  results$Time <- end_time - start_time

  #---------------- Does the file needs to be plotted? -----------------------

  if (plot == TRUE){
    message("Plotting the results")
    PlotPeacoQC(ff = results$FinalFF,
                show_peaks = results,
                output_directory = output_directory,
                channels = channels,
                acc_score = paste0(round(results$PercentageRemoved, 3),"% of the data was removed."))
  }

  return(results)

}


# ---------------------- PeacoQC: detection of peaks --------------------------


DetermineAllPeaks <- function(ff, channel, different_breaks){

  full_channel_peaks <- find_them_peaks(ff@exprs[,channel])$Peaks

  if (length(full_channel_peaks)<1){
    # warning(paste0("There seems to be an issue in channel ", channel,
    #                ". No peaks were found."))
    # plot(density(ff@exprs[,channel]), main = paste0("Density: ", channel))
    return(NA)
  }


  peak_results <- list()

  marker <- FlowSOM::get_markers(ff, channel)

  minimum <- min(ff@exprs[,channel])
  maximum <- max(ff@exprs[,channel])
  range <- abs(minimum) + abs(maximum)
  set.seed(1)

  splits <- lapply(different_breaks, function(x){find_them_peaks(ff@exprs[x,channel])})

  peaks <- lapply(splits, with, Peaks)

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
  medians_to_use <- apply(limited_nr_peaks, 2, median)



  # Trying to find the medians for the bins with the smaller amount of peaks
  #and using these as initial cluster centra kmeans

  if (length(medians_to_use) > 1) {

    # peak_clustering <- kmeans(peak_frame[,2], centers = peaks_medians)
    peak_clustering <- kmeans(peak_frame[,2], centers = medians_to_use)

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
    final_medians <- median(peak_frame$Peak)

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


# ----------------------------- Find them peaks -------------------------------


find_them_peaks <- function (obj)
{
  x <- obj
  n <- which(!is.na(x))
  if (length(n) < 3) {
    return(list(Peaks = NA, P.ind = 0, P.h = 0))
  }

  data <- x

  dens <- density(data[which(!is.na(data))], adjust = 1)
  dens <- smooth.spline(dens$x, dens$y, spar = 0.6)
  dens$y[which(dens$y < 0)] <- 0

  all.peaks <- actual_find_them_peaks(dens$y)

  peaks <- list("Heights" = all.peaks[,1], "Peaks" = dens$x[all.peaks[,2]])

  return(peaks)
}


actual_find_them_peaks <- function (x, nups = 1, ndowns = nups,
                                    minpeakheight = 1/4)
{
  # Eerste afgeleide
  xc <- paste(as.character(sign(diff(x))), collapse = "")
  xc <- gsub("1", "+", gsub("-1", "-", xc))
  peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
  rc <- gregexpr(peakpat, xc)[[1]]

  if (rc[1] < 0)
    return(NULL)
  x1 <- rc
  x2 <- rc + attr(rc, "match.length")
  attributes(x1) <- NULL
  attributes(x2) <- NULL
  n <- length(x1)
  xv <- xp <- numeric(n)
  for (i in 1:n) {
    xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
    xv[i] <- x[xp[i]]
  }
  minpeakheight <- max(x) * minpeakheight

  inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >=
                  0)
  X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
  if (length(X) == 0)
    return(c())
  return(X)
}



# ------------------------------- Remove duplicate peaks ----------------------

DuplicatePeaks <- function(i, data, medians ){
  to_use <- data[data$Bin == i,]
  duplicates <- to_use$Cluster[ duplicated(to_use$Cluster)] %>% as.numeric()
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
  mediancluster <-  median(data[ind])
  names(mediancluster) <- i
  return(mediancluster)
}


# ------------------------- Remove too small clusters -------------------------

#if the cluster is smaller than 1/10th of the binnr, remove it!

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



FlagOutliers <- function(peak, MAD) {


  kernel <- ksmooth(seq_along(peak), peak, x.points = seq_along(peak), bandwidth = 50)
  median_peak <- median(kernel$y, na.rm=TRUE)
  mad_peak <- mad(kernel$y)

  upper_interval <- median(median_peak, na.rm=TRUE)+MAD*(mad_peak)
  lower_interval <- median(median_peak, na.rm=TRUE)-MAD*(mad_peak)

  outliers <- ifelse(kernel$y > upper_interval, TRUE, ifelse(kernel$y < lower_interval, TRUE, FALSE))
  return(outliers)}

Flag.Outliers<- function(peak_frame, MAD) {
  to_remove_bins_df <- apply(peak_frame, 2, FlagOutliers, MAD)

  return(to_remove_bins_df)
}


# -------------------------------- Binning ------------------------------------

splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)

  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

splitWithOverlapMids <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)

  mids <- sapply(1:length(starts), function(i) starts[i] + ceiling(overlap/2))
  mids <- mids[mids < max(vec)]
  return(mids)
}


# ----------------------- Isolation Tree ------------------------------------



isolationTreeSD <- function(x, max_depth = as.integer(ceiling(log2(nrow(x)))),
                            gain_limit = 0.55, plot = FALSE){

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
        base_sd <- sd(x_col)
        gain_max_col <- 0
        split_v <- NA
        mean_val <- c()
        for(i in c(1:(length(x_col)-1))){

          sd_1 <- sd(x_col[c(1:i)])
          sd_2 <- sd(x_col[c((i+1):length(x_col))])
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


  if (plot == TRUE){
    relations <- data.frame(from = c(res$id, res$id),
                            to = c(res$left_child, res$right_child))

    relations <- na.omit(relations)

    g_1 <- graph_from_data_frame(relations, directed = TRUE, vertices = res)
    plot.igraph(g_1,
                layout = layout_as_tree(g_1),
                vertex.label = rowSums(selection))
  }


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

  pd <- pData(parameters(ff))
  new_pd_name <-  ncol(ff@exprs) +1
  new_pd_name <- paste0("$P", new_pd_name)
  new_pd <- data.frame(name = "GoodCells", desc = "PeacoQC_Cells", range = 1,
                       minRange = 0, maxRange = 1)
  rownames(new_pd) <- new_pd_name
  pd <- rbind(pd, new_pd)
  ff@exprs <- cbind(ff@exprs, "GoodCells" = as.numeric(bad_cells))
  pData(parameters(ff)) <- pd
  ff

}
