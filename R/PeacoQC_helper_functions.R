# ---------------------------- Make nice warning messages ---------------------

StrMessage <- function(x, prefix=" ", initial=""){
    strwrap(x, prefix=prefix, initial=initial)
}

# ----------------------------- Determine all peaks  --------------------------

DeterminePeaksAllChannels <- function(ff, channels, breaks,
                                        remove_zeros, results,
                                      peak_removal,min_nr_bins_peakdetection){

    peak_values_bin  <- list()

    # ------ Determine all peaks found in all channels and store them ---------
    i <- 0
    message("Calculating peaks")
    pb <- utils::txtProgressBar(min=0,
                                max=length(channels),
                                initial=0,
                                char="+",
                                style=3,
                                width=50,
                                file=stderr())
    utils::setTxtProgressBar(pb, i)

    for (channel in channels){
        i <- i +1

        channel_data <- flowCore::exprs(ff)[,channel]
        peak_frame <- DetermineAllPeaks(channel_data, breaks, remove_zeros,
                                        peak_removal,min_nr_bins_peakdetection)
        if (is(peak_frame, "logical")){
            utils::setTxtProgressBar(pb, i)
            next()
        }
        peak_values_bin[[channel]] <- ExtractPeakValues(peak_frame, breaks)
        results[[channel]] <- peak_frame
        utils::setTxtProgressBar(pb, i)
    }

    close(pb)


    all_peak_values <- unlist(peak_values_bin, recursive=FALSE)
    all_peaks <- as.data.frame(do.call(cbind, all_peak_values))
    rownames(all_peaks) <- names(all_peak_values[[1]])

    return(list("all_peaks"=all_peaks, "results"=results))
}


DetermineAllPeaks <- function(channel_data, breaks, remove_zeros, peak_removal,
                              min_nr_bins_peakdetection){

    full_channel_peaks <- FindThemPeaks(channel_data,
                                        remove_zeros,
                                        peak_removal)

    if(all(is.na(full_channel_peaks))) return(NA)

    peak_results <- list()

    minimum <- min(channel_data)
    maximum <- max(channel_data)
    range <- abs(minimum) + abs(maximum)


    channel_breaks <- lapply(breaks, function(x)channel_data[x])


    peaks <- list()

    for (channel_break in names(channel_breaks)){
        result_peak <- FindThemPeaks(channel_breaks[[channel_break]],
                                     remove_zeros, peak_removal)
        peaks[[channel_break]] <- result_peak
    }

    peak_frame <- plyr::ldply(peaks, cbind)

    colnames(peak_frame) <- c("Bin", "Peak")
    if (length(which(is.na(peak_frame$Peak))) > 0){

        peak_frame <- peak_frame[!is.na(peak_frame[, 2]), ]
        peaks <-peaks[!is.na(peaks)]

    }

    lengths <- vapply(peaks, length, FUN.VALUE=numeric(1))
    nr_peaks <- unique(lengths)
    nr_peaks_in_bins <- tabulate(match(lengths, nr_peaks))

    most_occuring_peaks <- max(nr_peaks[which(nr_peaks_in_bins
                                              > (min_nr_bins_peakdetection/100)*length(lengths))])

    ind_bins_nr_peaks <- lengths == most_occuring_peaks
    limited_nr_peaks <- do.call(rbind, peaks[ind_bins_nr_peaks])
    medians_to_use <- apply(limited_nr_peaks, 2, stats::median)

    if (length(medians_to_use) > 1) {

      names(medians_to_use) <- c(1:length(medians_to_use))

      peak_frame <- cbind(peak_frame,
                          "Cluster"= factor(sapply(peak_frame[,2],
                                            function(x){which.min(abs(x-medians_to_use))})))
      final_medians <- medians_to_use

    } else { # If only one peak was found per bin, no kmeans has to happen
        peak_frame <- cbind(peak_frame, "Cluster"=rep("1",
                                                      dim(peak_frame)[1]))
        final_medians <- stats::median(peak_frame$Peak)
    }

    bin_list <- lapply(unique(peak_frame$Bin),
                       DuplicatePeaks,
                       peak_frame,
                       final_medians)

    peak_frame <- do.call(rbind, bin_list)
    rownames(peak_frame) <- seq_len(dim(peak_frame)[1])

    peak_frame <- TooSmallClusters(peak_frame, peak_frame$Cluster)

    return(peak_frame)
}


FindThemPeaks <- function (channel_data, remove_zeros, peak_removal = (1/3))
{


    if (remove_zeros){
        # Remove all zeros before calculating densities
        channel_data <- channel_data[channel_data !=0]
        }

    if (length(channel_data) < 3) {
        return(NA)
    }

    dens <- stats::density(channel_data[!is.na(channel_data)], adjust=1)

    if (all(is.na(dens)))
        return(NA)

    n <- length(dens$y)
    selection <- (dens$y[2:(n-1)] > dens$y[seq_len(n-2)]) &
        (dens$y[2:(n-1)] > dens$y[3:n] ) &
        (dens$y[2:(n-1)] > (peak_removal * max(dens$y)))
    peaks <- dens$x[-1][selection]

    if (all(is.na(peaks))) {
        peaks <- dens$x[which.max(dens$y)]
    }

    return(peaks)
}


# ------------------------------- Remove duplicate peaks ----------------------

DuplicatePeaks <- function(i, data, medians ){
    to_use <- data[data$Bin == i, ]
    duplicates <- as.numeric(to_use$Cluster[ duplicated(to_use$Cluster)])
    if(length(duplicates)>0){

        for (duplex in duplicates){
            ind <- which(to_use$Cluster == duplex)

            diff <- abs(to_use$Peak[to_use$Cluster == duplex] - medians[duplex])
            to_use <- to_use[-ind[which.max(diff)], ]
        }

        to_use <- to_use[order(to_use$Cluster), ]
    }

    return(to_use)
}

# ------------------------ Find the medians per cluster -----------------------

clusterMedian=function(i, data, clusters) {
    ind=(clusters == i)
    mediancluster <-  stats::median(data[ind])
    names(mediancluster) <- i
    return(mediancluster)
}


# ------------------------- Remove too small clusters -------------------------

TooSmallClusters <- function(data, clusters ){
  nr_bins <- max(as.numeric(data$Bin))
  to_remove_clusters <- names(table(clusters))[which(table(clusters) < (1/2)*nr_bins)]
  data <- data[!data$Cluster %in% to_remove_clusters,]

    return(data)
}


# ------------------ Function to remove outliers based on mad -----------------



MADOutliers <- function(peak, MAD) {


    kernel <- stats::smooth.spline(seq_along(peak), peak, spar=0.5)


    median_peak <- stats::median(kernel$y, na.rm=TRUE)
    mad_peak <- stats::mad(kernel$y)

    upper_interval <- stats::median(median_peak, na.rm=TRUE)+MAD*(mad_peak)
    lower_interval <- stats::median(median_peak, na.rm=TRUE)-MAD*(mad_peak)

    outliers <- ifelse(kernel$y > upper_interval, TRUE,
                        ifelse(kernel$y < lower_interval,
                                TRUE, FALSE))
    return(outliers)}

MADOutlierMethod<- function(peaks, outlier_bins, MAD, breaks, nr_cells) {

    peak_frame <- peaks[outlier_bins, , drop = FALSE]

    names_bins <- which(outlier_bins)

    to_remove_bins_df <- apply(peak_frame, 2, MADOutliers, MAD)
    outlier_bins_MAD <- apply(to_remove_bins_df, 1, any)
    contribution_MAD <- apply(to_remove_bins_df, 2, function(x){
        names(x) <- names_bins
        removed <- RemovedBins(breaks, x, nr_cells)
        round((length(removed$cell_ids)/nr_cells)*100, 2)
    })
    names(outlier_bins_MAD) <- which(outlier_bins)

    return(list("MAD_bins"=outlier_bins_MAD,
                "Contribution_MAD"=contribution_MAD))
}


# -------------------------------- Binning ------------------------------------

SplitWithOverlap <- function(vec, seg.length, overlap) {
    starts=seq(1, length(vec), by=seg.length-overlap)
    ends  =starts + seg.length - 1
    ends[ends > length(vec)]=length(vec)

    lapply(seq_along(starts), function(i) vec[starts[i]:ends[i]])
}

SplitWithOverlapMids <- function(vec, seg.length, overlap) {
    starts=seq(1, length(vec), by=seg.length-overlap)
    ends  =starts + seg.length - 1
    ends[ends > length(vec)]=length(vec)

    mids <- vapply(seq_along(starts),
                    function(i) starts[i] + ceiling(overlap/2),
                    FUN.VALUE=numeric(1))
    mids <- mids[mids < max(vec)]
    return(mids)
}


# ----------------------- Isolation Tree ------------------------------------



isolationTreeSD <- function(x, max_depth=as.integer(ceiling(log2(nrow(x)))),
                            gain_limit){

    res <- data.frame(id=1,
                        left_child=NA,
                        right_child=NA,
                        gain=NA,
                        split_column=NA,
                        split_value=NA,
                        depth=0,
                        to_split=TRUE,
                        path_length=NA )
    selection <- matrix(TRUE, ncol=nrow(x))


    nodes_to_split <- base::which(res$to_split)

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
                # x_col <- x[rows, col]


                base_sd <- stats::sd(x_col)
                gain_max_col <- 0
                split_v <- NA
                mean_val <- c()
                for(i in seq_len((length(x_col)-1))){

                    sd_1 <- stats::sd(x_col[seq_len(i)])
                    sd_2 <- stats::sd(x_col[c((i+1):length(x_col))])
                    mean_val <- c(mean_val, mean(x_col[seq_len(i)]))


                    if (i == 1){
                        sd_1 <- 0
                    } else if (i == length(x_col) - 1){
                        sd_2 <- 0
                    }
                    gain <- (base_sd - mean(c(sd_1, sd_2)))/base_sd

                    if (is.na(gain)){
                        next()
                    }
                    if(gain >=gain_max_col){
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
                left_node <- selection[node, ] & x[ , split_col] <= split_value
                right_node <- selection[node, ] & x[ , split_col] > split_value

                if(length(which(left_node)) == length(rows) ||
                length(which(right_node)) == length(rows)){
                    res[node, "to_split"] <- FALSE
                    nodes_to_split <- which(res$to_split)
                    n_datapoints <- length(which(selection[, node]))
                    res[node, "path_length"] <- avgPL(n_datapoints) +
                        res[node, "depth"]
                    next()
                }

                max_node <- nrow(res)
                left_id <- max_node + 1
                right_id <- max_node + 2

                res[node, "left_child"] <- left_id
                res[node, "right_child"] <- right_id
                res[node, "gain"] <- gain_max
                res[node, "split_column"] <- colnames(x)[split_col]
                res[node, "split_value"] <- split_value
                res[node, "to_split"] <- FALSE

                gain_limit <- gain_max
                res <- rbind(res,
                            data.frame(id=c(left_id, right_id),
                                        left_child=NA,
                                        right_child=NA,
                                        gain=NA,
                                        split_column=NA,
                                        split_value=NA,
                                        depth=res[node, "depth"]+1,
                                        to_split=TRUE,
                                        path_length=NA))

                selection <- rbind(selection, left_node, right_node)


            } else {
                res[node, "to_split"] <- FALSE
                n_datapoints <- length(which(selection[node, ]))
                res[node, "path_length"] <- avgPL(n_datapoints) +
                    res[node, "depth"]
            }

        } else {
            res[node, "to_split"] <- FALSE
            n_datapoints <- length(which(selection[node, ]))
            res[node, "path_length"] <- avgPL(n_datapoints) + res[node, "depth"]


        }

        nodes_to_split <- base::which(res$to_split)

    }

    res$n_datapoints=rowSums(selection)
    res$anomaly_score=2^(-(res$path_length)/(avgPL(sum(res$n_datapoints))))


    scores_to_use <- stats::na.omit(res[,
                                        c("n_datapoints", "anomaly_score")])
    nodes_to_keep <- rownames(scores_to_use)[
        which.max(scores_to_use$n_datapoints)]


    outlier_bins <- selection[as.numeric(nodes_to_keep), ]
    names(outlier_bins) <- seq_along(outlier_bins)


    return(list("res"=res,
                "outlier_bins"=outlier_bins))
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

AppendCellID <- function(new_ff, cell_id){

    matrix_cell_id <- matrix(data = cell_id, ncol = 1,
                             dimnames = list(c(), list("Original_ID")))
    new_ff <- flowCore::fr_append_cols(new_ff, matrix_cell_id)
    return(new_ff)

}


# ------------------------- Make Breaks ---------------------------------------


MakeBreaks <- function(events_per_bin, nr_events){
    breaks <- SplitWithOverlap(seq_len(nr_events),
                                events_per_bin,
                                ceiling(events_per_bin/2))

    names(breaks) <- seq_along(breaks)

    return(list("breaks"=breaks, "events_per_bin"=events_per_bin))
}


FindEventsPerBin <- function(remove_zeros, ff,
                             channels,
                             min_cells = 150,
                             max_bins = 500,
                             step = 500
                             ){

  nr_events <- nrow(ff)

  if (remove_zeros){
    max_bins_mass <- min(apply(flowCore::exprs(ff)[,channels, drop = FALSE], 2,
                          function(x)sum(x != 0)))/min_cells

    if (max_bins_mass < max_bins){
      max_bins <- max_bins_mass
    }

  }


  max_cells <- ceiling((nr_events/max_bins)*2)
  max_cells <- ((max_cells%/%step)*step) + step

  events_per_bin <- max(min_cells, max_cells)

  return(events_per_bin)

    return(events_per_bin)
}


# ------------------ Find Increasing or decreasing channels -------------------

FindIncreasingDecreasingChannels <- function(breaks, ff, channels, plot,
                                             filename){

    weird_channel_decreasing <- list()
    weird_channel_increasing <- list()

    for (channel in channels){

        channel_data <- flowCore::exprs(ff)[,channel]

        channel_medians <- vapply(breaks,
                                    function(x){
                                        stats::median(channel_data[x])},
                                    FUN.VALUE=numeric(1))

        smoothed <- stats::ksmooth(seq_along(channel_medians),
                                    channel_medians,
                                    x.points=seq_along(channel_medians),
                                    bandwidth=50)

        increasing <- cummax(smoothed$y) == smoothed$y
        decreasing <- cummin(smoothed$y) == smoothed$y

        marker_name <- flowCore::getChannelMarker(ff, channel)$desc
        if(is.na(marker_name))
            marker_name <- channel

        if (length(which(increasing)) > (3/4)*length(increasing)){
            weird_channel_increasing <- c(weird_channel_increasing, channel)
        } else if (length(which(decreasing)) > (3/4)*length(decreasing)){
            weird_channel_decreasing <- c(weird_channel_decreasing, channel)
        }

    }

    if (length(weird_channel_decreasing) > 0 |
        length(weird_channel_increasing) > 0) {
        warning(StrMessage(c("There seems to be an increasing or decreasing
            trend in a channel ",
            " for ", basename(filename),
            ". Please inspect this in the overview figure.")))
        if (plot != FALSE){
            plot <- TRUE
        }
    }

    changing_channel <- "No increasing or decreasing effect"
    if (length(weird_channel_decreasing)>0){
        changing_channel <- "Decreasing channel"
    }
    if (length(weird_channel_increasing)>0){
        changing_channel <- "Increasing channel"
    }
    if (length(weird_channel_decreasing) >0 &
        length(weird_channel_increasing)>0){
        changing_channel <- "Increasing and decreasing channel"
    }


    return(list("Increasing"=weird_channel_increasing,
                "Decreasing"=weird_channel_decreasing,
                "Changing_channel"=changing_channel,
                "plot" = plot))
}


# ------------------------- Extract peak values --------------------------------

ExtractPeakValues <- function(peak_frame, breaks){

    peak_values <- list()

    for (peak in unique(peak_frame$Cluster)){
        peak_ind <- peak_frame$Cluster == peak
        peak_data <- peak_frame[peak_ind, ]
        peak_vector <- rep(stats::median(peak_data$Peak),
                            length(breaks))

        peak_vector[as.numeric(peak_data$Bin)] <- peak_data$Peak

        names(peak_vector) <- seq_len(length(peak_vector))

        peak_values[[peak]] <- peak_vector

    }

    return(peak_values)
}


# ---------------------- Removed cells ----------------------------------------

RemovedCells <- function(breaks, outlier_bins){

    removed_cells <- unlist(breaks[names(outlier_bins)
                                    [which(outlier_bins)]])
    removed_cells <- removed_cells[!duplicated(removed_cells)]

    return(removed_cells)

}

# ------------------------------ Removed bins ---------------------------------

RemovedBins <- function(breaks, outlier_bins, nr_cells){
    removed_cells <- RemovedCells(breaks, outlier_bins)

    # Summary of entire analysis
    bad_cells <- rep(TRUE, nr_cells)
    bad_cells[removed_cells] <- FALSE

    return(list("cells"=bad_cells,
                "cell_ids"=removed_cells))
}

# ------------------------- Check input PeacoQCSS ----------------------------
CheckInputSignalStability <- function(ff, channels, determine_good_cells, plot,
                                      save_fcs, output_directory, report,
                                      time_channel_parameter){
    if(!is(ff, "flowFrame") | is.null(ff))
        stop("ff should be a flowFrame.")
    if(is.null(output_directory) & save_fcs)
        warning("Since the output directory is NULL,
            no fcs files will be stored.")
    if(is.null(output_directory) & report)
        warning("Since the output directory is NULL,
            no report will be generated.")
    if(is.null(output_directory) & plot)
        warning("Since the output directory is NULL,
            no plots will be generated.")
    if(!(determine_good_cells %in% c("all", "IT", "MAD", FALSE)))
        stop(StrMessage("The parameter determine_good_cells should be of
            following values: all, IT or MAD."))
    if(!is.numeric(channels) & !all(channels%in% colnames(flowCore::exprs(ff)))|
    is.null(channels))
        stop(StrMessage("Make sure that you use indices or the colnames in the
            expression matrix in the flowframe to indicate which channels you
            want to use."))

    # Check for time channel
    time_channel <- grep(time_channel_parameter, colnames(flowCore::exprs(ff)),
                            ignore.case=TRUE)
    if (any(diff(flowCore::exprs(ff)[, time_channel]) < 0))
        warning(StrMessage("There is an inconsistancy in the time channel.
            It seems that not all the cells are ordered according to time
            in the flowframe."))


    if (nrow(ff) < 500){
        warning("There are less than 500 cells available in the flowFrame.
                This could result in unsufficient cells that are neccessary
                for the IT or the MAD analysis.")
    }

}

# -------------------------------- Remove short regions ------------------------
RemoveShortRegions <- function(ff,
                                outlier_bins,
                                consecutive_bins,
                                breaks,
                                results){
    # Determine which cells should also be removed because they
    # fall inbetween FALSE regions

    nr_cells <- nrow(ff)
    outlier_bins_new <- inverse.rle(within.list(rle(outlier_bins),
                                    values[lengths<consecutive_bins] <- FALSE))
    consecutive_bins_to_remove <- !(outlier_bins_new == outlier_bins)
    names(consecutive_bins_to_remove) <- seq_along(outlier_bins_new)

    consecutive_cells <- RemovedBins(breaks,
                                        consecutive_bins_to_remove,
                                        nr_cells)
    names(outlier_bins_new) <- names(outlier_bins)
    outlier_bins <- outlier_bins_new
    bad_cells <- RemovedBins(breaks, !outlier_bins, nr_cells)

    results$ConsecutiveCells <- consecutive_cells$cells
    results$ConsecutiveCellsPercentage <- (length(consecutive_cells$cell_ids)/
                                              nr_cells)*100
    results$GoodCells <- bad_cells$cells
    results$PercentageRemoved <- (length(bad_cells$cell_ids)/
                                    nr_cells)*100
    results$FinalFF <- ff[results$GoodCells, ]
    return(results)
}


