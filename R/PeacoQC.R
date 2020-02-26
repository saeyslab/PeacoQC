#' @title Remove margin events
#'
#' @description \code{RemoveMargins} will remove margin events from the
#' flowframe based on the internal description of the fcs file.
#'
#' @usage RemoveMargins(ff, channels,
#' channel_specifications = NULL, output = "frame")
#'
#' @param ff A flowframe
#' @param channels The channel indices or channel names that have to be checked
#' for margin events
#' @param channel_specifications A list of lists with parameter specifications
#' for certain channels. This parameter should only be used if the values in
#' the internal parameters description is too strict or wrong for a number or
#' all channels. This should be one list per channel with first a minRange and
#' then a maxRange value. This list should have the channel name found back in
#' \code{colnames(ff@exprs)}. If a channel is not listed in this parameter, its
#' default internal values will be used.
#' @param output If set to "full", a list with the filtered flowframe and the
#' indices of the margin event is returned. If set to "frame", only the
#' filtered flowframe is returned.
#'
#' @examples
#' # Read in raw data
#' fileName <- system.file("extdata", "111.fcs", package = "PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels where the margin events should be removed
#' channels <- c(1,3,5:14,18,21)
#'
#' # Remove margins
#'
#' ff_cleaned <- RemoveMargins(ff, channels)
#'
#' # If an internal value is wrong for a channels (e.g. FSC-A)
#'
#' channel_specifications <- list("FSC-A" = c(-111, 262144))
#' ff_cleaned <- RemoveMargins(
#'     ff,
#'     channels,
#'     channel_specifications = channel_specifications)
#' @return This function returns either a filtered flowframe when the
#' \code{output} parameter is set to "frame" or a list containing the filtered
#' flowframe and a TRUE/FALSE list indicating the margin events.
#' @importFrom flowWorkspace pData
#' @importFrom methods is
#' @export

RemoveMargins <- function(
    ff,
    channels,
    channel_specifications = NULL,
    output = "frame") {

    if (!is(ff, "flowFrame"))
        stop("ff should be a flowframe.")
    if (!is(channel_specifications, "list") &
            !is.null(channel_specifications))
        stop("channel_specifications should be a list of lists.")
    if(!all(lengths(channel_specifications) == 2) &
            !is.null(channel_specifications))
        stop(StrMessage("channel_specifications should be a list of lists.
            Every list should have the channel name and should contain a
            minRange and maxRange value."))
    if(is.null(names(channel_specifications)) &
            !is.null(channel_specifications))
        stop(StrMessage("channel_specifications should be a list of named lists.
            Make sure that the names correspend with the channel names."))
    if(!all(names(channel_specifications) %in% colnames(ff@exprs)) &
            !is.null(channel_specifications))
        stop(StrMessage("channel_specifications should be a list of named lists.
            Make sure that the names correspend with the channel names."))
    if(!is.numeric(channels) & !all(channels%in% colnames(ff@exprs)) |
            is.null(channels))
        stop(StrMessage("Make sure that you use indices or the colnames in the
            expression matrix in the flowframe to indicate which channels you
            want to use."))

    meta <- flowWorkspace::pData(ff@parameters)
    rownames(meta) <- meta[, "name"]

    if(!is.null(channel_specifications)){
        meta[names(channel_specifications),
            c("minRange", "maxRange")] <- do.call(rbind, channel_specifications)
    }

    # Initialize variables
    selection <- rep(TRUE, times = dim(ff)[1])
    e <- ff@exprs

    if(is.numeric(channels)){
        channels <- colnames(ff@exprs)[channels]
    }
        # Make selection
    for (d in channels) {

        selection <- selection &
            e[, d] > max(min(meta[d, "minRange"], 0), min(e[, d])) &
            e[, d] < min(meta[d, "maxRange"], max(e[, d]))
    }


    if (length(which(selection == FALSE))/length(selection) > 0.1) {
        warning(StrMessage(paste0("More then ",
            round(length(which(selection == FALSE))/length(selection) * 100, 2),
            "% is considered as a margin event in file ",
            basename(ff@description$FILENAME), ". This should be verified.")))
    }
    if (output == "full"){
        return(
            list("flowframe" = ff[selection,],
                "indices_margins" = which(selection == FALSE)))
    } else if (output == "frame"){
        return(ff[selection,])
    }
}
#' @title Peak-based detection of high quality cytometry data
#'
#' @description \code{PeacoQC} will determine peaks on the channels in the
#' flowframe. Then it will remove anomalies caused by e.g. clogs, changes in
#' speed etc. by using an IsolationTree and/or the MAD method.
#'
#' @usage
#' PeacoQCSignalStability(ff,channels, determine_good_cells = "all",
#'         plot = TRUE, save_fcs = TRUE, output_directory = ".",
#'         name_directory = "PeacoQC_results", report = TRUE,
#'         events_per_bin = 2000, MAD = 5, IT_limit = 0.55,
#'         consecutive_bins = 5, remove_zeros = FALSE, suffix_fcs = "_QC", ...)
#'
#' @param ff A flowframe or the location of an fcs file
#' @param channels Indices or names of the channels in the flowframe on which
#' peaks have to be determined.
#' @param determine_good_cells If set to FALSE, the algorithm will only
#' determine peaks. If it is set to "all", the bad measurements will be
#' filtered out based on the MAD and IT analysis. It can also be put to "MAD"
#' or "IT" to only use one method of filtering.
#' @param plot If set to TRUE, the \code{PlotPeacoQC} function is run to make
#' an overview plot of the deleted measurements.
#' @param save_fcs If set to TRUE, the cleaned fcs file will be saved in the
#' \code{output_directory} as: filename_QC.fcs.
#' @param output_directory Directory where a new folder will be created that
#' consists of the generated fcs files, plots and report. If set to NULL,
#' nothing will be stored.
#' @param name_directory Name of folder that will be generated in
#' \code{output_directory}.
#' @param report Overview text report that is generated after PeacoQC is run.
#' If set to FALSE, no report will be generated.
#' @param events_per_bin Number of events that are put in one bin.
#' Default is 2000.
#' @param MAD The MAD parameter. Default is 5. If this is increased, the
#' algorithm becomes less strict.
#' @param IT_limit The IsolationTree parameter. Default is 0.55. If this is
#' increased, the algorithm becomes less strict.
#' @param consecutive_bins If 'good' bins are located between bins that are
#' removed, they will also be marked as 'bad'. The default is 5.
#' @param remove_zeros If this is set to TRUE, the zero values will be removed
#' before the peak detection step. They will not be indicated as 'bad' value.
#' This is recommended when cleaning mass cytometry data.
#' @param suffix_fcs The suffix given to the new fcs files.
#' @param ... Options to pass on to the \code{PlotPeacoQC} function
#' (display_cells, manual_cells, prefix)
#'
#' @return This function returns a \code{list} with a number of items. It will
#' include "FinalFF" where the transformed, compensated and cleaned flowframe is
#' stored. It also contains the starting parameters and the information
#' necessary to give to \code{PlotPeacoQC} if the two functions are run
#' seperatly. The GoodCells list is also given where 'good' measurements are
#' indicated as TRUE and the to be removed measurements as FALSE.
#'
#' @examples
#' # General pipeline for preprocessing and quality control with PeacoQC
#'
#' # Read in compensated and transformed data
#' fileName <- system.file("extdata", "111_Comp_Trans.fcs", package = "PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels where the margin events should be removed
#' # and on which the quality control should be done
#' channels <- c(1,3,5:14,18,21)
#'
#' #Run PeacoQC
#' PeacoQC_res <- PeacoQCSignalStability(ff, channels,
#'                     determine_good_cells = "all",
#'                     plot = TRUE, save_fcs = TRUE)
#'
#' @importFrom methods is
#'
#' @export
#'

PeacoQCSignalStability <- function(ff,
    channels,
    determine_good_cells = "all",
    plot = TRUE,
    save_fcs = TRUE,
    output_directory = ".",
    name_directory = "PeacoQC_results",
    report = TRUE,
    events_per_bin = 2000,
    MAD = 5,
    IT_limit = 0.55,
    consecutive_bins = 5,
    remove_zeros = FALSE,
    suffix_fcs = "_QC",
    ...
){

    if(!is(ff, "flowFrame") | is.null(ff))
        stop("ff should be a flowFrame.")
    if(is.null(output_directory) & save_fcs == TRUE)
        warning("Since the output directory is NULL,
            no fcs files will be stored.")
    if(is.null(output_directory) & report == TRUE)
        warning("Since the output directory is NULL,
            no report will be generated.")
    if(is.null(output_directory) & plot == TRUE)
        warning("Since the output directory is NULL,
            no plots will be generated.")
    if(!(determine_good_cells %in% c("all", "IT", "MAD", FALSE)))
        stop(StrMessage("The parameter determine_good_cells should be of
            following values: all, IT or MAD."))
    if(!is.numeric(channels) & !all(channels%in% colnames(ff@exprs)) |
            is.null(channels))
        stop(StrMessage("Make sure that you use indices or the colnames in the
            expression matrix in the flowframe to indicate which channels you
            want to use."))

    # Check for time channel
    time_channel <- grep("time", colnames(ff@exprs), ignore.case = TRUE)
    if (any(diff(ff@exprs[,time_channel]) < 0) == TRUE)
        warning(StrMessage("There is an inconsistancy in the Time channel.
            It seems that not al the cells are ordered according to time
            in the flowframe."))

    # Make a new directory where all results will be stored
    if(!is.null(output_directory)){
        suppressWarnings(dir.create(output_directory))
        storing_directory <- file.path(output_directory, name_directory)
        suppressWarnings(dir.create(storing_directory))
        if (save_fcs == TRUE){
            fcs_directory <- file.path(storing_directory, "fcs_files")
            suppressWarnings(dir.create(fcs_directory))
        }
        if (report == TRUE & determine_good_cells %in% c("all", "IT", "MAD")){
            # Make initial percentages in case only one analysis is used
            perc_IT <- "Not_used"
            perc_MAD <- "Not_used"

            report_location <- file.path(storing_directory,
                paste0("PeacoQC_report",".txt"))
            if (!file.exists(report_location)) {
                utils::write.table(t(c("Filename", "Nr. Measurements",
                    "% Full analysis", "Analysis by", "% MAD analysis",
                    "% IT analysis", "MAD","IT limit", "Consecutive bins",
                    "events_per_bin", "Increasing/Decreasing channel")),
                    report_location, sep = "\t", row.names = FALSE,
                    quote = FALSE, col.names = FALSE)
            }
        }
    } else { # If the directory is NULL, no things can be saved in any location
        plot <- FALSE
        save_fcs <- FALSE
        report <- FALSE
    }

    # Searching for the name of the ff
    if (length(ff@description$FILENAME)>0){
        message(paste0("Starting quality control analysis for ",
            basename(ff@description$FILENAME)))
    }

    # Make an empty list for the eventual results
    results <- list()
    results$Analysis <- determine_good_cells

    # Make sure that channels only consist out of the colnames of ff@exprs
    results$Channels <- channels
    if (is.numeric(channels)){
        channels <- colnames(ff@exprs)[channels]
    }

    # Split the ff up in bins (overlapping)
    breaks <- SplitWithOverlap(seq_len(nrow(ff)),
        events_per_bin,
        ceiling(events_per_bin/2))

    names(breaks) <- seq_along(breaks)


    # If not enough bins are made, at least 100 should be present
    if (length(breaks) < 100){
        events_per_bin <- ceiling(nrow(ff)/100) *2
        breaks <- SplitWithOverlap(seq_len(nrow(ff)),
            events_per_bin,
            ceiling(events_per_bin/2))

        names(breaks) <- seq_along(breaks)

    }

    results$EventsPerBin <- events_per_bin

    # Check if there is an increasing or decreasing trent in the channels
    weird_channel_increasing <- NULL
    weird_channel_decreasing <- NULL


    for (channel in channels){

        channel_medians <- vapply(breaks,
            function(x){stats::median(ff@exprs[x,channel])},
            FUN.VALUE = numeric(1))

        smoothed <- stats::ksmooth(seq_along(channel_medians),
            channel_medians,
            x.points = seq_along(channel_medians),
            bandwidth = 50)

        increasing <- cummax(smoothed$y) == smoothed$y
        decreasing <- cummin(smoothed$y) == smoothed$y

        marker_name <- flowCore::getChannelMarker(ff, channel)$desc
        if(is.na(marker_name))
            marker_name <- channel

        if (length(which(increasing)) > (3/4)*length(increasing)){
            warning(StrMessage(paste0(
                "There seems to be an increasing trent in channel ",
                marker_name,
                " for ",basename(ff@description$FILENAME),
                ". Please inspect this before doing any further analysis.")),
                plot(channel_medians, main = marker_name))
            weird_channel_increasing <- c(weird_channel_increasing, channel)
        } else if (length(which(decreasing)) > (3/4)*length(decreasing)){
            warning(StrMessage(paste0(
                "There seems to be a decreasing trent in channel ",
                marker_name,
                " for ",basename(ff@description$FILENAME),
                ". Please inspect this before doing any further analysis.")),
                plot(channel_medians, main = marker_name))
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
        char = "+", style = 3, width = 50, file = stderr())
    utils::setTxtProgressBar(pb,i)

    channel <- channels[10]


    for (channel in channels){
        i <- i +1

        peak_frame <- DetermineAllPeaks(ff, channel, breaks, remove_zeros)

        if (is(peak_frame, "logical")){
            utils::setTxtProgressBar(pb,i)
            next()
        }

        peak_values <- list()

        for (peak in unique(peak_frame$Cluster)){
            peak_ind <- peak_frame$Cluster == peak
            peak_data <- peak_frame[peak_ind,]
            peak_vector <- rep(stats::median(peak_data$Peak),
                length(breaks))

            peak_vector[as.numeric(peak_data$Bin)] <- peak_data$Peak

            names(peak_vector) <- seq_len(length(peak_vector))

            peak_values[[peak]] <- peak_vector

        }

        peak_values_bin[[channel]] <- peak_values

        results[[channel]] <- peak_frame

        utils::setTxtProgressBar(pb,i)
    }

    close(pb)

    peak_values <- unlist(peak_values_bin, recursive = FALSE)

    all_peaks <- do.call(cbind, peak_values)

    all_peaks <- as.data.frame(all_peaks)

    rownames(all_peaks) <- names(peak_values[[1]])


    outlier_bins <- rep(TRUE, nrow(all_peaks))



    # ------------------------ Isolation Tree  --------------------------------

    if (determine_good_cells == "all" || determine_good_cells == "IT"){

        if (length(which(outlier_bins)) < 3){
            stop(StrMessage(paste0("There is an issue with file ",
                basename(ff@description$FILENAME),
                ". There are no good cells left to perform quality control.")))
        }

        data_for_trees <- all_peaks[outlier_bins,]

        tree <- isolationTreeSD(data_for_trees, gain_limit = IT_limit)
        results$IT <- tree

        scores_to_use <- stats::na.omit(tree$res[,
            c("n_datapoints", "anomaly_score")])
        node_to_keep <- rownames(scores_to_use)[
            which.max(scores_to_use$n_datapoints)]

        IT_cells <- rep(TRUE, nrow(ff))


        bins_isolated <- tree$selection[as.numeric(node_to_keep),]
        names(bins_isolated) <- which(outlier_bins == TRUE)

        outlier_bins[outlier_bins] <- bins_isolated

        removed_cells_IT <- unlist(breaks[names(bins_isolated)
            [which(bins_isolated == FALSE)]])
        removed_cells_IT <- removed_cells_IT[!duplicated(removed_cells_IT)]

        IT_cells[removed_cells_IT] <- FALSE

        results$OutlierIT <- IT_cells


        perc_IT <- length(removed_cells_IT)/nrow(ff)
        perc_IT <- perc_IT * 100

        message(paste0("IT analysis removed ",
            paste0(round(perc_IT,2),
                " % of the measurements" )))

        results$ITPercentage <- perc_IT



    }


    # ------------------------ Outliers based on mad --------------------------

    # Bins that were selected based on their mad


    if (determine_good_cells == "all" || determine_good_cells == "MAD"){

        if (length(which(outlier_bins)) < 3){
            stop(StrMessage(paste0("There is an issue with file ",
                basename(ff@description$FILENAME),
                ". There are no good cells left to perform quality control.")))
        }

        perc_not_outliers <- length(which(outlier_bins == TRUE))/
            length(outlier_bins)


        outlier_bins_df <- MADOutlierMethod(all_peaks[outlier_bins,],MAD)

        outlier_bins_MAD <- apply(outlier_bins_df, 1, any)

        names(outlier_bins_MAD) <- which(outlier_bins == TRUE)

        contribution_MAD <- apply(outlier_bins_df, 2,
            function(x){length(which(x == TRUE))/
                    (length(x)*perc_not_outliers)})

        contribution_MAD <- round(contribution_MAD*100,2)

        mad_cells <- rep(TRUE, nrow(ff))

        removed_cells_MAD <- unlist(breaks[names(outlier_bins_MAD)
            [which(outlier_bins_MAD == TRUE)]])
        removed_cells_MAD <- removed_cells_MAD[!duplicated(removed_cells_MAD)]

        perc_MAD <- length(removed_cells_MAD)/nrow(ff)

        perc_MAD <- perc_MAD * 100

        message(paste0("MAD analysis removed ",
            paste0(round(perc_MAD,2),
                " % of the measurements" )))

        results$ContributionMad <- contribution_MAD

        results$MADPercentage <- perc_MAD

        mad_cells[removed_cells_MAD] <- FALSE

        results$OutlierMads <- mad_cells

        outlier_bins[outlier_bins] <- !outlier_bins_MAD

    }



    # ------------------------- indicate bad cells ----------------------------


    if (determine_good_cells %in% c("all", "IT", "MAD")){

        # Determine which cells should also be removed because they
        # fall inbetween FALSE regions
        outlier_bins_new <- inverse.rle(within.list(rle(outlier_bins),
            values[lengths<consecutive_bins]
            <- FALSE))

        consecutive_cells <- rep(TRUE, nrow(ff))

        consecutive_bins_to_remove <- !(outlier_bins_new == outlier_bins)

        names(consecutive_bins_to_remove) <- seq_along(outlier_bins_new)

        removed_cells_consec <- unlist(breaks[
            names(consecutive_bins_to_remove)
            [which(consecutive_bins_to_remove == TRUE)]])
        removed_cells_consec <- removed_cells_consec[
            !duplicated(removed_cells_consec)]


        consecutive_cells[removed_cells_consec] <- FALSE

        results$ConsecutiveCells <- consecutive_cells


        outlier_bins <- outlier_bins_new

        to_remove_bins <- which(outlier_bins == FALSE)


        removed_cells <- unlist(breaks[to_remove_bins])

        removed_cells <- removed_cells[!duplicated(removed_cells)]


        # Summary of entire analysis
        bad_cells <- rep(TRUE, nrow(ff))

        bad_cells[removed_cells] <- FALSE

        ff_orig <- Append_GoodCells(ff,  bad_cells)


        percentage_removed <- (length(which(bad_cells == FALSE))/
                length(bad_cells))*100

        message(paste0("The algorithm removed ",
            round(percentage_removed,2),
            "% of the measurements"))

        results$remove_bins <- to_remove_bins
        results$GoodCells <- bad_cells
        results$PercentageRemoved <- percentage_removed
        results$FinalFF <- ff[results$GoodCells,]

        if(percentage_removed > 70){
            warning(StrMessage(paste0("There was ", round(percentage_removed,3),
                "% of the measurements removed for file ",
                basename(ff@description$FILENAME))))
        }


        # -----------------  Does the file need to be saved in an fcs? ---------
        if (save_fcs == TRUE & !is.null(output_directory)){
            message("Saving fcs file")
            flowCore::write.FCS(ff[results$GoodCells,],file.path(fcs_directory,
                paste0(sub(".fcs",
                    "",
                    basename(ff@description$FILENAME)),
                    paste0(suffix_fcs,".fcs"))))
        }


        # ----------------- Does an overview file need to be generated? --------
        changing_channel <- "No increasing or decreasing effect"
        if (report == TRUE & !is.null(output_directory)){
            if (!is.null(weird_channel_decreasing)){
                changing_channel <- "Decreasing channel"
            }
            if (!is.null(weird_channel_increasing)){
                changing_channel <- "Increasing channel"
            }
            if (!is.null(weird_channel_decreasing) &
                    !is.null(weird_channel_increasing)){
                changing_channel <- "Increasing and decreasing channel"
            }


            utils::write.table(t(c(
                basename(ff@description$FILENAME),
                nrow(ff),
                percentage_removed,
                determine_good_cells,perc_MAD,
                perc_IT, MAD, IT_limit, consecutive_bins,
                events_per_bin, changing_channel)),
                report_location, sep = "\t", append = TRUE, row.names = FALSE,
                quote = FALSE, col.names = FALSE)
        }
    }


    #---------------- Does the file need to be plotted? ------------------------

    if (plot == TRUE & !is.null(output_directory)){
        if(determine_good_cells %in% c("all", "IT", "MAD")){
            title_FR = paste0(round(results$PercentageRemoved, 3),
                "% of the data was removed.")
        } else {
            title_FR = ""
        }

        message("Plotting the results")
        PlotPeacoQC(ff = ff,
            display_peaks = results,
            output_directory = storing_directory,
            channels = results$Channels,
            title_FR = title_FR,
            ...)
    }

    return(results)

}
#' @title Visualise deleted cells of PeacoQC
#'
#' @description \code{PlotPeacoQC} will generate a png file with on overview of
#' the flow rate and the different selected channels. These will be annotated
#' based on the measurements that were removed by PeacoQC. It is also possible
#' to only display the quantiles and median or only the measurements without
#' any annotation.
#'
#' @usage
#' PlotPeacoQC(ff, channels, output_directory = ".", display_cells = 5000,
#'             manual_cells = NULL, title_FR = NULL, display_peaks = TRUE,
#'             prefix = "PeacoQC_", time_unit = 100, ...)
#'
#' @param ff A flowframe
#' @param channels Indices of names of the channels in the flowframe that have
#' to be displayed
#' @param output_directory Directory where the plots should be generated. Set
#' to NULL if no plots need to be generated
#' @param display_cells The number of measurements that should be displayed.
#' (The number of dots that are displayed for every channel)
#' @param manual_cells Give a vector (TRUE/FALSE) with annotations for each cell
#' to compare the automated QC with.
#' @param title_FR The title that has to be displayed above the flow rate figure
#' @param display_peaks If the result of \code{PeacoQC} is given, all the
#' quality control results will be visualised. If set to TRUE: \code{PeacoQC}
#' will be run and only the peaks will be displayed without any quality control.
#' If set to FALSE, no peaks will be displayed and only the events will be
#' displayed.
#' @param prefix The prefix that will be given to the generated png file
#' @param time_unit The number of time units grouped together for visualising
#' event rate. The default is set to 100, resulting in events per second for
#' most flow datasets. Suggested to adapt for mass cytometry data.
#' @param ... Arguments to be given to \code{PeacoQC} if \code{display_peaks}
#' is set to TRUE.
#'
#' @return This function returns nothing but generates a png file in the
#' output_directory
#'
#' @import ggplot2
#' @importFrom methods is
#'
#' @examples
#'
#' ## Plotting the results of PeacoQC
#'
#' # Read in transformed and compensated data
#' fileName <- system.file("extdata", "111_Comp_Trans.fcs", package = "PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels on which the quality control should be done and the
#' # plots should be made
#' channels <- c(1,3,5:14,18,21)
#'
#' # Run PeacoQC
#' PeacoQC_res <- PeacoQCSignalStability(ff,
#'     channels,
#'     determine_good_cells = "all",
#'     plot = FALSE,
#'     save_fcs = TRUE)
#'
#' # Run PlotPeacoQC
#' PlotPeacoQC(ff, channels, display_peaks = PeacoQC_res)
#'
#' ## Plot only the peaks (No quality control)
#' PlotPeacoQC(ff, channels, display_peaks = TRUE)
#'
#' ## Plot only the dots of the file
#' PlotPeacoQC(ff, channels, display_peaks = FALSE)
#'
#' @export
PlotPeacoQC <- function(ff,
    channels,
    output_directory = ".",
    display_cells = 5000,
    manual_cells = NULL,
    title_FR = NULL,
    display_peaks = TRUE,
    prefix = "PeacoQC_",
    time_unit = 100,
    ...) {

    requireNamespace("ggplot2")


    if(!is(ff, "flowFrame") | is.null(ff))
        stop("ff should be a flowFrame.")
    if(!is.numeric(channels) & !all(channels%in% colnames(ff@exprs)) |
            is.null(channels))
        stop(StrMessage("Make sure that you use indices or the colnames in the
            expression matrix in the flowframe to indicate which channels you
            want to use."))
    if(is.null(output_directory))
        stop("There should be a path given to the output_directory parameter.")


    # Check for time channel
    time_channel <- grep("time", colnames(ff@exprs), ignore.case = TRUE)
    if (any(diff(ff@exprs[,time_channel]) < 0) == TRUE)
        warning(StrMessage("There is an inconsistancy in the Time channel.
            It seems that not al the cells are ordered according
            to time in the flowframe."))


    # Make a new directory where all results will be stored
    if(!is.null(output_directory)){

        plot_directory <- file.path(output_directory, "PeacoQC_plots")
        suppressWarnings(dir.create(plot_directory))
    }
    # Plot acc score if one is listed in given arguments
    if (!is.null(title_FR)) {
        scores_time <- title_FR
    } else {
        scores_time <- ""
    }

    # If display_peaks == TRUE, the peaks should be calculated by using PeacoQC
    if (is(display_peaks, "logical")){
        if(display_peaks == TRUE){
            message("Running PeacoQC to determine peaks")
            peaks <- PeacoQCSignalStability(ff,
                channels,
                determine_good_cells = FALSE,
                output_directory = NULL, plot = FALSE, save_fcs = FALSE,
                report = FALSE,
                ...)
        } else{peaks <- FALSE}
    } else { peaks <- display_peaks}

    filename <- basename(ff@description$FILENAME)

    # Name to put on plotfile
    name <- sub(".fcs", "", filename)

    n_channels <- length(channels)
    if (is.numeric(channels)){
        channels <- colnames(ff@exprs)[channels]
    }

    # Determining grid to plot
    n_row <- floor(sqrt(n_channels + 2))
    n_col <- ceiling((n_channels + 2)/n_row)

    # Subset of points that will be plotted
    if (nrow(ff) >= 50000) {
        subset_timeplot <- sort(sample(seq_len(nrow(ff)), 50000))
    } else {
        subset_timeplot <- seq_len(nrow(ff))
    }

    if (display_cells > nrow(ff)) {
        message(StrMessage("There are less then the number of display cells
            available. Setting the number of display cells to the number of
            measurements."))
        display_cells <- nrow(ff)
    }

    subset_signalplot <- sort(sample(seq_len(nrow(ff)), display_cells))



    # Calculate backgroundvalues for points on plot, rectangle block and
    # CV values after automated qc

    if (is(display_peaks, "list")){

        # Make blocks for automated gated algorithms to display on plot
        run_length <- rle(peaks$GoodCells)

        full_QC_vector <- ifelse(peaks$GoodCells == TRUE, TRUE, FALSE)


        consecutive_cells <- which(peaks$ConsecutiveCells == FALSE)
        if (length(consecutive_cells) > 0) {
            full_QC_vector[!peaks$ConsecutiveCells] <- "consecutive"
            run_length <- rle(full_QC_vector)

        }

        mad_cells <- which(peaks$OutlierMads == FALSE)
        if (length(mad_cells) > 0) {
            full_QC_vector[!peaks$OutlierMads] <- "mad"
            run_length <- rle(full_QC_vector)

        }

        fill_blocks <- ifelse(run_length$values == TRUE, "Good values",
            ifelse(run_length$values == "consecutive", "In consecutive bins",
                ifelse(run_length$values ==
                        FALSE, "IT", "MAD")))

        x_min <- c(1, cumsum(run_length$lengths)[-length(run_length$lengths)])
        x_max <- cumsum(run_length$lengths)

        overview_blocks <- data.frame(x_min = x_min,
            x_max = x_max,
            y_min = -Inf,
            y_max = Inf,
            fill_blocks = fill_blocks)

        x_min <- ff@exprs[, "Time"][c(
            1,
            cumsum(run_length$lengths)[-length(run_length$lengths)])]
        x_max <- ff@exprs[, "Time"][cumsum(run_length$lengths)]

        overview_blocks_background <- data.frame(x_min = x_min,
            x_max = x_max,
            y_min = -Inf,
            y_max = Inf,
            fill_blocks = fill_blocks)

    }

    # Make blocks for manual gates to display on signalplot + cv score
    # after manual gating flowrate
    if (!is.null(manual_cells)) {
        run_length_man <- rle(manual_cells)


        fill_blocks_man <- ifelse(
            run_length_man$values,
            "deepskyblue1",
            "indianred1")

        x_min_man <- c(1,
            cumsum(run_length_man$lengths)[-length(run_length_man$lengths)])
        x_max_man <- cumsum(run_length_man$lengths)

        manual_blocks <- data.frame(x_min = x_min_man,
            x_max = x_max_man,
            y_min = -Inf,
            y_max = Inf,
            fill_blocks = fill_blocks_man)
    }


    # Initiate empty lists for all plots
    plot_list <- list()


    if (length(time_channel) > 0){
        # Calculating time breaks
        h <- graphics::hist(ff@exprs[subset_timeplot, "Time"],
            breaks = seq(min(ff@exprs[,"Time"]), max(ff@exprs[, "Time"]) +
                    time_unit, by = time_unit),
            plot = FALSE)

        idcs <- findInterval(ff@exprs[subset_timeplot, "Time"], h$breaks)


        # Build time plot (FlowRate)

        p_time <- ggplot() + theme_bw()

        p_time <- p_time + theme(panel.grid = element_blank())


        if (is(display_peaks, "list")){

            p_time <- p_time + geom_rect(data = overview_blocks_background,
                mapping = aes(xmin = x_min,
                    xmax = x_max,
                    ymin = y_min,
                    ymax = y_max,
                    fill = fill_blocks),
                alpha = 0.4,
                show.legend = TRUE) + scale_fill_manual(name = "",
                    values = c(IT = "indianred1",
                        MAD = "mediumpurple1",
                        `In consecutive bins` = "plum1",
                        `Good Values` = "white"),
                    guide = guide_legend(override.aes = list(alpha = 0.4)))
            p_time <- p_time + theme(legend.key = element_rect(colour = "snow4"))

        }

        p_time <- p_time +
            geom_point(aes(x = h$mids, y = h$counts)) +
            ggtitle(scores_time) +
            xlab("Time") + ylab("Nr of cells per second") +
            theme(plot.title = element_text(size = 10))



        # Save time plot in plot list
        plot_list[["Time"]] <- p_time

    }


    # -------------------------------- Individual channels --------------------

    channel <- channels[2]

    if(is(peaks, "list")){
        events_per_bin <- peaks$EventsPerBin

        mid_breaks <- SplitWithOverlapMids(seq_len(nrow(ff)),
            events_per_bin,
            ceiling(events_per_bin/2))
        m <-  nrow(ff@exprs)%%events_per_bin


        # if not zero, you have something to add to your sequence
        if(m != 0) {

            mid_breaks = c(mid_breaks, nrow(ff@exprs))
        }
    }


    for (channel in channels) {

        marker <- flowCore::getChannelMarker(ff, channel)$desc
        if(is.na(marker))
            marker <- channel

        minimum <- min(ff@exprs[, channel])
        maximum <- max(ff@exprs[, channel])
        range <- abs(minimum) + abs(maximum)

        if (is(display_peaks, "list")){

            # Show contributions of every channel in MAD and IF

            contribution_MAD <- sum(peaks$ContributionMad[grep(channel,
                names(peaks$ContributionMad))])
            contribution_IT <- ifelse(length(grep(channel,
                peaks$IT$res$split_column)) > 0, "+", "/")

            contributions <- paste0(marker, "\n", "IT: ", contribution_IT,
                " MAD: ", contribution_MAD, "%")

            if (length(grep(channel, peaks$WeirdChannels$Increasing)) > 0) {
                contributions <- paste0(contributions, "\n",
                    "WARNING: Increasing channel.")
            } else if (length(grep(channel,
                peaks$WeirdChannels$Decreasing))> 0){
                contributions <- paste0(contributions, "\n",
                    "WARNING: Decreasing channel.")
            }
        } else{contributions <- marker}


        flowdata <- data.frame(Cells = subset_signalplot,
            channel = ff@exprs[subset_signalplot, channel])


        # Initial plot
        p <- ggplot() +
            ylab("Value") +
            xlab("Cells") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0),
                panel.grid = element_blank())

        if (is(display_peaks, "list")) {

            p <- p + geom_rect(data = overview_blocks,
                mapping = aes(xmin = x_min,
                    xmax = x_max,
                    ymin = y_min,
                    ymax = y_max,
                    fill = fill_blocks),
                alpha = 0.4,
                show.legend = TRUE) +
                scale_fill_manual(name = "",
                    values = c(IT = "indianred1",
                        MAD = "mediumpurple1",
                        `In consecutive bins` = "plum1",
                        `Good Values` = "white"),
                    guide = guide_legend(override.aes = list(alpha = 0.4)))
            p <- p + theme(legend.key = element_rect(colour = "snow4"))

        }

        p <- p + geom_point(data = flowdata,
            aes(x = Cells, y = channel),
            size = 0.3,
            col = "snow4")

        if(is(peaks, "list")){

            peak_frame <- peaks[[channel]]

            if (is(peak_frame, "data.frame")) {

                peak_frame$Bin <- as.numeric(mid_breaks)[as.numeric(
                    peak_frame$Bin)]

                colours <- paste0("grey",
                    seq_len(20))[seq_len(
                        max(as.numeric(peak_frame[,"Cluster"])))]

                p <- p + geom_line(data = peak_frame, aes(x = Bin,
                    y = Peak,
                    color = Cluster),
                    size = 1,
                    show.legend = FALSE) +
                    scale_color_manual(values = colours)
            } else {
                contributions <- paste0(contributions, " No peak was found.")
            }
        }


        p <- p + ggtitle(contributions)


        # Add manual blocks to plot
        if (!is.null(manual_cells)) {

            if (length(grep("SC", channel)) == 0) {

                manual_blocks$y_min <- minimum - 1
                manual_blocks$y_max <- minimum - 1.25
            } else {
                manual_blocks$y_min <- minimum - 45000
                manual_blocks$y_max <- minimum - 58000
            }
            p <- p + geom_rect(data = manual_blocks,
                mapping = aes(xmin = x_min,
                    xmax = x_max,
                    ymin = y_min,
                    ymax = y_max),
                fill = manual_blocks$fill_blocks,
                alpha = 1)

            p <- p + theme(legend.position = "none")

        }

        p <- p + theme(plot.subtitle = element_text(size = 10, color = "black"))

        # plot(p)
        plot_list[[channel]] <- p

    }


    # To make nice plots
    g <- ggplotGrob(plot_list[[1]] + theme(legend.position = "bottom"))$grobs


    # Use original argument to see if there was a list given or not
    if (!(is(display_peaks, "logical"))){
        legend <- g[[which(vapply(g,
                function(x) x$name,
                FUN.VALUE = character(1)) == "guide-box")]]
        lheight <- sum(legend$height)
        plots <- do.call(gridExtra::arrangeGrob,
            c(lapply(plot_list, function(x)x + theme(legend.position = "none")),
                nrow = n_row,
                ncol = n_col))
        new <- gridExtra::arrangeGrob(
            plots,
            legend,
            heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight))
    } else {
        plots <- do.call(gridExtra::arrangeGrob, c(lapply(plot_list,
            function(x) x), nrow = n_row, ncol = n_col))
        new <- gridExtra::arrangeGrob(plots)}

    ggsave(paste0(plot_directory, "/", prefix, name, ".png"),
        new, width = n_col * 5, height = n_row * 3, limitsize = FALSE)

}

#' @title Make overview heatmap of quality control analysis
#'
#' @description \code{PeacoQCHeatmap} will make a heatmap to display all the
#' results generated by \code{PeacoQC}. It will include the percentages of
#' measurements that are removed in total, by the IT method and by the MAD
#' method. It will also show the parameters that were used during the
#' quality control.
#'
#' @usage
#' PeacoQCHeatmap(report_location, show_values = TRUE, show_row_names = TRUE,
#' latest_tests = FALSE, title = "PeacoQC report", ...)
#'
#' @param report_location The path to the PeacoQC report generated by
#' \code{PeacoQC}.
#' @param show_values If set to TRUE, the percentages of removed values
#' will be displayed on the heatmap.
#' @param show_row_names If set to FALSE, the filenames will not be displayed
#' on the heatmap.
#' @param latest_tests If this is set to TRUE, only the latest quality control
#' run will be displayed in the heatmap.
#' @param title The title that should be given to the heatmap.
#' @param ... Extra parameters to be given to the \code{Heatmap} function
#' (eg. row_split)
#' @return This function returns nothing but generates a heatmap that can be
#' saved as pdf or png
#'
#' @import ComplexHeatmap
#' @importFrom utils read.delim write.table
#' @importFrom grDevices png dev.off
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.text gpar
#'
#' @examples
#'
#' # Find path to PeacoQC report
#' location <- system.file("extdata", "PeacoQC_report.txt", package = "PeacoQC")
#'
#' # Make heatmap overview of quality control run
#' PeacoQCHeatmap(report_location = location)
#'
#' # Make heatmap with only the runs of the last test
#' PeacoQCHeatmap(report_location = location, latest_tests = TRUE)
#'
#' # Make heatmap with row annotation
#' PeacoQCHeatmap(report_location = location,
#'     row_split = c("r1", "r2", rep("r3",2), rep("r4", 16)))
#'
#' @importFrom grDevices colorRampPalette
#' @export

PeacoQCHeatmap <- function(
    report_location,
    show_values = TRUE,
    show_row_names = TRUE,
    latest_tests = FALSE,
    title = "PeacoQC report",
    ...){

    if (!file.exists(report_location))
        stop(StrMessage("The path specified in the report_location parameter
            is wrong or incomplete."))

    if(show_row_names == FALSE & latest_tests == FALSE)
        warning(StrMessage("If there are duplicates in the report file,
            they will be displayed on the heatmap without their filename."))

    report_table <- utils::read.delim(report_location,
        check.names = FALSE,
        stringsAsFactors = FALSE)
    report_table[report_table == "Not_used"] <- NA



    if (latest_tests == TRUE){
        report_table <- report_table[!duplicated(report_table$Filename,
            fromLast = TRUE),]
        rownames(report_table) <- report_table$Filename
    } else {

        # Find the duplicates and rename them

        table_names <- report_table$Filename
        y <- table(table_names) > 1
        tmp <- !duplicated(table_names) & (table_names %in% names(which(y)))
        unique_table_names <- make.unique(table_names)

        unique_table_names[tmp] <- paste0(unique_table_names[tmp], ".0")
        unique_table_names <- sub("(.*.fcs)(\\.)(.*)",
            "\\1_\\3", unique_table_names)

        rownames(report_table)  <- unique_table_names
    }
    annotation_frame <- data.frame(
        "Consecutive bins" = factor(report_table$`Consecutive bins`),
        "IT limit" = factor(report_table$`IT limit`),
        "MAD" = factor(report_table$MAD),check.names = FALSE)

    rownames(annotation_frame) <- rownames(report_table)

    t1 <- colorRampPalette(c("#8D99AE", "#2B2D42"))
    col_cons <- t1(length(unique(annotation_frame$`Consecutive bins`)))
    t2 <- colorRampPalette(c("#EBB9DF", "#7D1D3F"))
    col_MAD <- t2(length(unique(annotation_frame$MAD)))
    t3 <- colorRampPalette(c("#B2CEDE", "#AD7A99"))
    col_IT <- t3(length(unique(annotation_frame$`IT limit`)))
    names(col_cons) <- unique(annotation_frame$`Consecutive bins`)
    names(col_MAD) <- unique(annotation_frame$MAD)
    names(col_IT) <- unique(annotation_frame$`IT limit`)


    analysis <- report_table$`Analysis by`

    if("all" %in% analysis | all(c("IT", "MAD") %in% analysis)){
        ha <- rowAnnotation(df = annotation_frame,
            col = list("Consecutive bins" = col_cons,
                "MAD" = col_MAD,
                "IT limit" = col_IT))

    } else if(length(unique(analysis)) == 1 & unique(analysis) == "IT"){
        ha <- rowAnnotation(
            "Consecutive bins" = annotation_frame$`Consecutive bins`,
            "IT limit" = annotation_frame$`IT limit`,
            col = list("Consecutive bins" = col_cons,
                "IT limit" = col_IT))
    }else if(length(unique(analysis)) == 1 & unique(analysis) == "MAD"){
        ha <- rowAnnotation(
            "Consecutive bins" = annotation_frame$`Consecutive bins`,
            "MAD" = annotation_frame$MAD,
            col = list("Consecutive bins" = col_cons,
                "MAD" = col_MAD))

    }

    col_incr_decr_channel <- c()

    if("No increasing or decreasing effect" %in%
            report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
            "No increasing or decreasing effect" = "#26C485")
    }
    if ("Increasing channel" %in% report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
            "Increasing channel" = "#AF3800")
    }
    if ("Decreasing channel" %in% report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
            "Decreasing channel" = "#721817")
    }
    if ("Increasing and decreasing channel" %in%
            report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
            "Increasing and decreasing channel" = "#A50104")
    }


    annotation_right <- rowAnnotation(
        "Incr/Decr" = as.factor(report_table[,
            "Increasing/Decreasing channel"]),
        col = list("Incr/Decr" = col_incr_decr_channel))




    report_matrix <- data.matrix(report_table[,c(3,5,6)])

    if(show_values == TRUE){
        cell_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.1f", report_matrix[i, j]), x, y,
                gp = gpar(fontsize = 10))
        }
    } else{
        cell_fun = NULL
    }


    ph <- Heatmap(report_matrix,
        cell_fun = cell_fun,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_title = title,
        column_title_gp = grid::gpar(fontface = "bold"),
        col = circlize::colorRamp2(c(0,20,100), c("#EBEBD3", "#FFD151","red")),
        left_annotation = ha,
        right_annotation = annotation_right,
        heatmap_legend_param = list(direction = "horizontal"),
        name = "Removed percentage",
        na_col = "Grey",
        show_row_names = show_row_names,
        ...)

    draw(ph, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")


}
#' @title Full PeacoQC pre-processing pipeline
#'
#' @description Method to run general pre-processing and quality control
#' workflow. The method will first remove margins and compensate and transform
#' the flowframe. This is followed by the \code{PeacoQCSignalStability}
#' function that will check the data on anomalies that occured during
#' measurement.
#'
#' @usage
#'
#' PeacoQC(ff, channels, remove_margins = TRUE, compensation_matrix,
#'         transformation_list, channel_specifications = NULL,
#'         determine_good_cells = "all", plot = TRUE, save_fcs = TRUE,
#'         output_directory = ".", name_directory = "PeacoQC_results",
#'         report = TRUE, events_per_bin = 2000, MAD = 5, IT_limit = 0.55,
#'         consecutive_bins = 5, remove_zeros = FALSE, suffix_fcs = "_QC",  ...)
#'
#' @param ff A flowframe or the location of an fcs file
#' @param channels Indices or namers of the channels in the flowframe on which
#' peaks have to be determined.
#' @param remove_margins If set to FALSE, the margins will not be removed.
#' @param channel_specifications A list of lists with parameter specifications
#' for certain channels. This parameter should only be used if the values in
#' the internal parameters description is too strict or wrong for a number or
#' all channels. This should be one list per channel with first a minRange and
#' then a maxRange value. This list should have the channel name found back in
#' \code{colnames(ff@exprs)}. If a channel is not listed in this parameter, its
#' default internal values will be used.
#' @param compensation_matrix The compensation matrix that will be used by the
#' flowCore function compensate.
#' @param transformation_list The transformation list for all the channels
#' that should be transformed.
#' @param determine_good_cells If set to FALSE, the algorithm will only
#' determine peaks. If it is set to "all", the bad measurements will be
#' filtered out based on the MAD and IT analysis. It can also be put to "MAD"
#' or "IT" to only use one method of filtering.
#' @param plot If set to TRUE, the \code{PlotPeacoQC} function is run to make
#' an overview plot of the deleted measurements.
#' @param save_fcs If set to TRUE, the compensated, transformed and cleaned fcs
#' file will be saved in the \code{output_directory} as: filename_QC.fcs.
#' @param output_directory Directory where a new folder will be created that
#' consists of the generated fcs files, plots and report. If set to NULL,
#' nothing will be stored.
#' @param name_directory Name of folder that will be generated in
#' \code{output_directory}.
#' @param report Overview text report that is generated after PeacoQC is run.
#' If set to FALSE, no report will be generated.
#' @param events_per_bin Number of events that are put in one bin.
#' Default is 2000.
#' @param MAD The MAD parameter. Default is 5. If this is increased, the
#' algorithm becomes less strict.
#' @param IT_limit The IsolationTree parameter. Default is 0.55. If this is
#' increased, the algorithm becomes less strict.
#' @param consecutive_bins If 'good' bins are located between bins that are
#' removed, they will also be marked as 'bad'. The default is 5.
#' @param remove_zeros If this is set to TRUE, the zero values will be removed
#' before the peak detection step. They will not be indicated as 'bad' value.
#' This is recommended when cleaning mass cytometry data.
#' @param suffix_fcs The suffix given to the new fcs files.
#' @param ... Options to pass on to the \code{PlotPeacoQC} function
#' (display_cells, manual_cells, prefix)
#'
#' @return This function returns a \code{list} with a number of items. It will
#' include "FinalFF" where the transformed, compensated and cleaned flowframe is
#' stored. It also contains the starting parameters and the information
#' necessary to give to \code{PlotPeacoQC} if the two functions are run
#' seperatly. The GoodCells list is also given where 'good' measurements are
#' indicated as TRUE and the to be removed measurements as FALSE.
#'
#' @examples
#' # General pipeline for preprocessing and quality control with PeacoQC
#'
#' # Read in raw data
#' fileName <- system.file("extdata", "111.fcs", package = "PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels where the margin events should be removed
#' # and on which the quality control should be done
#' channels <- c(1,3,5:14,18,21)
#'
#' # Compensation matrix (is most of the time stored in the flowframe as
#' # flowCore::description(ff)$SPILL or flowCore::description(ff)$SPILLOVER)
#' compensation_matrix <- flowCore::description(ff)$SPILL
#'
#' # Store the transformation list
#' transformation_list <- flowCore::estimateLogicle(ff,
#'                                 colnames(compensation_matrix))
#'
#' #Run PeacoQC
#' PeacoQC_res <- PeacoQC(ff = ff,
#'                         channels = channels,
#'                         compensation_matrix = compensation_matrix,
#'                         transformation_list = transformation_list)
#'
#' @importFrom flowCore compensate transform
#'
#' @export

PeacoQC <- function(
    ff,
    channels,
    remove_margins = TRUE,
    compensation_matrix,
    transformation_list,
    channel_specifications = NULL,
    determine_good_cells = "all",
    plot = TRUE,
    save_fcs = TRUE,
    output_directory = ".",
    name_directory = "PeacoQC_results",
    report = TRUE,
    events_per_bin = 2000,
    MAD = 5,
    IT_limit = 0.55,
    consecutive_bins = 5,
    remove_zeros = FALSE,
    suffix_fcs = "_QC",
    ...
) {


    if (is.null(compensation_matrix) & remove_margins == TRUE)
        warning(StrMessage("If the data is already compensated, it could be that
            the RemoveMargins function will not work correctly."))

    if(remove_margins == TRUE){
        #Remove margins
        ff <- RemoveMargins(
            ff = ff,
            channels = channels,
            channel_specifications = channel_specifications)
    }


    if (!is.null(compensation_matrix)){

        #Compensation
        ff <- flowCore::compensate(ff, compensation_matrix)
    }


    if (!is.null(transformation_list)){
        #Transformation
        ff <- flowCore::transform(ff, transformation_list)
    }

    #PeacoQCSignalStability
    results_peacoQC <- PeacoQCSignalStability(
        ff = ff,
        channels = channels,
        determine_good_cells = determine_good_cells,
        plot = plot,
        save_fcs = save_fcs,
        output_directory = output_directory,
        name_directory = name_directory,
        report = report,
        events_per_bin = events_per_bin,
        MAD = MAD,
        IT_limit = IT_limit,
        consecutive_bins = consecutive_bins,
        remove_zeros = remove_zeros,
        suffix_fcs = suffix_fcs,
        ...)

    return(results_peacoQC)

}











