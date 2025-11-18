#' @title Remove margin events of flow cytometry data
#'
#' @description \code{RemoveMargins} will remove margin events from the
#' flowframe based on the internal description of the fcs file.
#'
#' @usage RemoveMargins(ff, channels,
#' channel_specifications=NULL, output="frame")
#'
#' @param ff A flowframe that contains flow cytometry data.
#' @param channels The channel indices or channel names that have to be checked
#' for margin events
#' @param channel_specifications A list of vectors with parameter specifications
#' for certain channels. This parameter should only be used if the values in
#' the internal parameters description is too strict or wrong for a number or
#' all channels. This should be one list per channel with first a minRange and
#' then a maxRange value in a vector. This list should have the channel name
#' found back in \code{colnames(flowCore::exprs(ff))}. If a channel is not
#' listed in this parameter, its default internal values will be used. The
#' default of this parameter is NULL.
#' @param remove_min The channel indices or channel names that have to be checked
#' for minimum margins. The default of this parameter is channels.
#' @param remove_max The channel indices or channel names that have to be checked
#' for maximum margins. The default of this parameter is channels.
#' @param output If set to "full", a list with the filtered flowframe and the
#' indices of the margin event is returned. If set to "frame", only the
#' filtered flowframe is returned. The default is "frame".
#'
#' @examples
#' # Read in raw data
#' fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels where the margin events should be removed
#' channels <- c(1, 3, 5:14, 18, 21)
#'
#' # Remove margins
#'
#' ff_cleaned <- RemoveMargins(ff, channels)
#'
#' # If an internal value is wrong for a channels (e.g. FSC-A)
#'
#' channel_specifications <- list("FSC-A"=c(-111, 262144),
#'                                "SSC-A"=c(-111, 262144))
#' ff_cleaned <- RemoveMargins(
#'     ff,
#'     channels,
#'     channel_specifications=channel_specifications)
#' @return This function returns either a filtered flowframe when the
#' \code{output} parameter is set to "frame" or a list containing the filtered
#' flowframe and a TRUE/FALSE list indicating the margin events. An extra column
#' named "Original_ID" is added to the flowframe where the cells are given their
#'  original cell id.
#' @importFrom flowCore pData
#' @importFrom methods is
#' @export

RemoveMargins <- function(
        ff,
        channels,
        channel_specifications=NULL,
        output="frame",
        remove_min = channels,
        remove_max = channels) {

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
            Make sure that the names correspond with the channel names."))
    if(!all(names(channel_specifications) %in% colnames(flowCore::exprs(ff))) &
       !is.null(channel_specifications))
        stop(StrMessage("channel_specifications should be a list of named lists.
            Make sure that the names correspond with the channel names."))
    if(!is.numeric(channels) & !all(channels%in% colnames(flowCore::exprs(ff)))|
       is.null(channels))
        stop(StrMessage("Make sure that you use indices or the colnames in the
            expression matrix in the flowframe to indicate which channels you
            want to use."))

    meta <- flowCore::pData(flowCore::parameters(ff))
    rownames(meta) <- meta[, "name"]

    if(!is.null(channel_specifications)){
        meta[names(channel_specifications),
             c("minRange", "maxRange")] <- do.call(rbind, channel_specifications)
    }

    # Initialize variables
    selection <- rep(TRUE, times=dim(ff)[1])
    e <- flowCore::exprs(ff)

    if(is.numeric(channels)){
        channels <- colnames(flowCore::exprs(ff))[channels]
    }
    # Make selection
    margin_matrix <- matrix(nrow = length(channels), ncol = 2, dimnames = list(channels, c("min", "max")))
    for (d in channels) {

      if(d %in% remove_min){
        min_margin_ev <- e[, d] <= max(min(meta[d, "minRange"], 0), min(e[, d]))
        margin_matrix[d, "min"] <- sum(min_margin_ev)
        selection <- selection & !min_margin_ev
      }
      if(d %in% remove_max){
        max_margin_ev <- e[, d] > min(meta[d, "maxRange"], max(e[, d]))
        margin_matrix[d, "max"] <- sum(max_margin_ev)
        selection <- selection & !max_margin_ev
            
      }
      
    }

    if ((length(flowCore::keyword(ff)$FILENAME) > 0) &&
        !is.na(flowCore::keyword(ff)$FILENAME)) {
        filename <- basename(flowCore::keyword(ff)$FILENAME)
    } else if ((length(flowCore::keyword(ff)[["$FIL"]]) > 0) &&
               !is.na(flowCore::keyword(ff)[["$FIL"]])) {
        filename <- basename(flowCore::keyword(ff)[["$FIL"]])
    }

    if (length(which(selection == FALSE))/length(selection) > 0.1) {
        warning(StrMessage(c("More than ",
                             round(length(which(selection == FALSE))/length(selection) * 100, 2),
                             "% is considered as a margin event in file ",
                             filename,
                             ". This should be verified.")))
    }

    new_ff <- ff[selection, ]
    attr(new_ff, "margin_matrix") <- margin_matrix
    if (!("Original_ID" %in% colnames(flowCore::exprs(new_ff)))){
        new_ff <- AppendCellID(new_ff, which(selection))}
    if (output == "full"){
        return(
            list("flowframe"=new_ff,
                 "indices_margins"=which(selection == FALSE)))
    } else if (output == "frame"){
        return(new_ff)
    }
}

#' @title Remove doublet events from flow cytometry data
#'
#' @description \code{RemoveDoublets} will remove doublet events from the
#' flowframe based on two channels.
#'
#' @usage RemoveDoublets(ff, channel1="FSC-A", channel2="FSC-H", nmad=4,
#' verbose=FALSE, output="frame")
#'
#' @param ff A flowframe that contains flow cytometry data.
#' @param channel1 The first channels that will be used to determine the
#' doublet events. Default is "FSC-A"
#' @param channel2 The second channels that will be used to determine the
#' doublet events. Default is "FSC-H"
#' @param nmad Bandwidth above the ratio allowed (cells are kept if their
#' ratio is smaller than the median ratio + \code{nmad} times the median
#' absolute deviation of the ratios). Default is 4.
#' @param verbose If set to TRUE, the median ratio and width will be printed.
#' Default is FALSE.
#' @param output If set to "full", a list with the filtered flowframe and the
#' indices of the doublet event is returned. If set to "frame", only the
#' filtered flowframe is returned. The default is "frame".
#' @param b Parameter to manually shift the accepted ratio.
#'
#' @examples
#' # Read in data
#' fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Remove doublets
#' ff_cleaned <- RemoveDoublets(ff)
#' @return This function returns either a filtered flowframe when the
#' \code{output} parameter is set to "frame" or a list containing the filtered
#' flowframe and a TRUE/FALSE list indicating the margin events. An extra column
#' named "Original_ID" is added to the flowframe where the cells are given their
#' original cell id.
#' @importFrom methods is
#' @export

RemoveDoublets <- function(ff,
                           channel1="FSC-A",
                           channel2="FSC-H",
                           nmad=4,
                           verbose=FALSE,
                           output="frame",
                           b=0){

    if (!is(ff, "flowFrame"))
        stop("ff should be a flowframe.")

    # Calculate the ratios
    ratio <- flowCore::exprs(ff)[,channel1] /
        (1e-10+ flowCore::exprs(ff)[,channel2] + b)

    # Define the region that is accepted
    r <- stats::median(ratio)
    r_m <- stats::mad(ratio)
    if(verbose) message(paste0("Median ratio: ", r,", width: ", nmad*r_m))

    # Make selection
    selection <- ratio < r+nmad*r_m

    new_ff <- ff[selection, ]

    if (!("Original_ID" %in% colnames(flowCore::exprs(new_ff)))){
        new_ff <- AppendCellID(new_ff, which(selection))}

    if (output == "full"){
        return(
            list("flowframe"=new_ff,
                 "indices_doublets"=which(selection == FALSE)))
    } else if (output == "frame"){
        return(new_ff)
    }
}

#' @title Peak-based detection of high quality cytometry data
#'
#' @description \code{PeacoQC} will determine peaks on the channels in the
#' flowframe. Then it will remove anomalies caused by e.g. clogs, changes in
#' speed etc. by using an IsolationTree and/or the MAD method.
#'
#' @usage
#' PeacoQC(ff, channels, determine_good_cells="all",
#'         plot=20, save_fcs=TRUE, output_directory=".",
#'         name_directory="PeacoQC_results", report=TRUE,
#'         events_per_bin=FindEventsPerBin(remove_zeros, ff, channels,
#'         min_cells, max_bins, step), min_cells=150, max_bins=500, step=500,
#'         MAD=6, IT_limit=0.6, consecutive_bins=5, remove_zeros=FALSE,
#'         suffix_fcs="_QC", force_IT=150, peak_removal = (1/3),
#'         min_nr_bins_peakdetection = 10, time_channel_parameter = "Time",
#'          ...)
#'
#' @param ff A flowframe or the location of an fcs file. Make sure that the
#' flowframe is compensated and transformed. If it is mass cytometry data, only
#' a transformation is necessary.
#' @param channels Indices or names of the channels in the flowframe on which
#' peaks have to be determined.
#' @param determine_good_cells If set to FALSE, the algorithm will only
#' determine peaks. If it is set to "all", the bad measurements will be
#' filtered out based on the MAD and IT analysis. It can also be put to "MAD"
#' or "IT" to only use one method of filtering.
#' @param min_cells The minimum amount of cells (nonzero values) that should be
#' present in one bin. Lowering this parameter can affect the robustness of the
#' peak detection. Default is 150.
#' @param max_bins The maximum number of bins that can be used in the cleaning
#' process. If this value is lowered, larger bins will be made. Default is 500.
#' @param step The step in events_per_bin to which the parameter is reduced to.
#' Default is 500.
#' @param plot When PeacoQC removes more than the specified percentage, an
#' overview plot will be made of all the selected channels and the deleted
#' measurements. If set to TRUE, the \code{PlotPeacoQC} function is
#' run to make an overview plot of the deleted measurements, even when
#' nothing is removed. Default is set to 20. If an increasing or decreasing
#' trend is found, a figure will also be made except if plot is set to FALSE.
#' @param save_fcs If set to TRUE, the cleaned fcs file will be saved in the
#' \code{output_directory} as: filename_QC.fcs. The _QC name can be altered with
#' the \code{suffix_fcs} parameter. An extra column named "Original_ID" is added
#' to this fcs file where the cells are given their original cell id.
#' Default is TRUE.
#' @param output_directory Directory where a new folder will be created that
#' consists of the generated fcs files, plots and report. If set to NULL,
#' nothing will be stored.The default folder is the working directory.
#' @param name_directory Name of folder that will be generated in
#' \code{output_directory}. The default is "PeacoQC_results".
#' @param report Overview text report that is generated after PeacoQC is run.
#' If set to FALSE, no report will be generated. The default is TRUE.
#' @param events_per_bin Number of events that are put in one bin.
#' Default is calculated based on the rows in \code{ff}
#' @param MAD The MAD parameter. Default is 6. If this is increased, the
#' algorithm becomes less strict.
#' @param IT_limit The IsolationTree parameter. Default is 0.55. If this is
#' increased, the algorithm becomes less strict.
#' @param consecutive_bins If 'good' bins are located between bins that are
#' removed, they will also be marked as 'bad'. The default is 5.
#' @param remove_zeros If this is set to TRUE, the zero values will be removed
#' before the peak detection step. They will not be indicated as 'bad' value.
#' This is recommended when cleaning mass cytometry data. Default is FALSE.
#' @param suffix_fcs The suffix given to the new fcs files. Default is "_QC".
#' @param force_IT If the number of determined bins is less than this number,
#' the IT analysis will not be performed. Default is 150 bins.
#' @param peak_removal During the peak detection step, peaks are only kept if
#' they are \code{peak_removal} percentage of the maximum height peak. Default is
#' 1/3
#' @param min_nr_bins_peakdetection The percentage of number of bins in which
#' the maximum number of peaks has to be present. Default is 10.
#' @param time_channel_parameter Name of the time channel in ff if present.
#' Default is "Time".
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
#' # Read in raw fcs file
#' fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels where the margin events should be removed
#' # and on which the quality control should be done
#' channels <- c(1, 3, 5:14, 18, 21)
#'
#' ff <- RemoveMargins(ff=ff, channels=channels, output="frame")
#'
#' # Compensate and transform the data
#'
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
#' ff <- flowCore::transform(ff,
#'                             flowCore::estimateLogicle(ff,
#'                             colnames(flowCore::keyword(ff)$SPILL)))
#' #Run PeacoQC
#' PeacoQC_res <- PeacoQC(ff, channels,
#'                         determine_good_cells="all",
#'                         save_fcs=TRUE)
#'
#' @importFrom methods is
#'
#' @export
#'
PeacoQC <- function(ff,
                    channels,
                    determine_good_cells="all",
                    plot=20,
                    save_fcs=TRUE,
                    output_directory=".",
                    name_directory="PeacoQC_results",
                    report=TRUE,
                    events_per_bin=FindEventsPerBin(remove_zeros, ff,
                                                    channels, min_cells,
                                                    max_bins, step),
                    min_cells = 150,
                    max_bins = 500,
                    step = 500,
                    MAD=6,
                    IT_limit=0.6,
                    consecutive_bins=5,
                    remove_zeros=FALSE,
                    suffix_fcs="_QC",
                    force_IT = 150,
                    peak_removal = (1/3),
                    min_nr_bins_peakdetection = 10,
                    time_channel_parameter = "Time",
                    ...
){

    CheckInputSignalStability(ff, channels, determine_good_cells, plot,
                              save_fcs, output_directory, report,
                              time_channel_parameter)

    # Make a new directory where all results will be stored
    if(!is.null(output_directory)){
        suppressWarnings(dir.create(output_directory, recursive = TRUE))
        storing_directory <- file.path(output_directory, name_directory)
        suppressWarnings(dir.create(storing_directory))
        if (save_fcs){
            fcs_directory <- file.path(storing_directory, "fcs_files")
            suppressWarnings(dir.create(fcs_directory))
        }
        if (report & determine_good_cells %in% c("all", "IT", "MAD")){
            # Make initial percentages in case only one analysis is used
            perc_IT <- "Not_used"
            perc_MAD <- "Not_used"

            report_location <- file.path(storing_directory,
                                         paste0("PeacoQC_report", ".txt"))
            if (!file.exists(report_location)) {
                utils::write.table(t(c("Filename",
                                       "Nr. Measurements before cleaning",
                                       "Nr. Measurements after cleaning",
                                       "% Full analysis", "Analysis by",
                                       "% IT analysis", "% MAD analysis",
                                       "% Consecutive cells",
                                       "MAD", "IT limit", "Consecutive bins",
                                       "Events per bin",
                                       "Increasing/Decreasing channel")),
                                   report_location, sep="\t",
                                   row.names=FALSE,
                                   quote=FALSE, col.names=FALSE)
            }
        }
    } else { # If the directory is NULL, no things can be saved in any location
        plot <- FALSE
        save_fcs <- FALSE
        report <- FALSE
    }

    # Searching for the name of the ff

    if ((length(flowCore::keyword(ff)$FILENAME) > 0) &&
        !is.na(flowCore::keyword(ff)$FILENAME)) {
        filename <- basename(flowCore::keyword(ff)$FILENAME)
    } else if ((length(flowCore::keyword(ff)[["$FIL"]]) > 0) &&
               !is.na(flowCore::keyword(ff)[["$FIL"]])) {
        filename <- basename(flowCore::keyword(ff)[["$FIL"]])
    }
    if (length(filename)>0){
        message("Starting quality control analysis for ",
                basename(filename))
    }

    # Make an empty list for the eventual results
    results <- list()
    results$Analysis <- determine_good_cells

    results$Channels <- channels
    if (is.numeric(channels)){
        channels <- colnames(flowCore::exprs(ff))[channels]}

    # Make the breaks for the entire flowframe
    res_breaks <- MakeBreaks(events_per_bin, nrow(ff))
    breaks <- res_breaks$breaks
    results$nr_bins <- length(res_breaks$breaks)

    results$EventsPerBin <- res_breaks$events_per_bin

    # Check if there is an increasing or decreasing trend in the channels

    list_weird_channels <- FindIncreasingDecreasingChannels(breaks,
                                                            ff, channels,
                                                            plot,
                                                            filename)
    plot <- list_weird_channels$plot
    results$WeirdChannels <- list_weird_channels[1:3]

    all_peaks_res <- DeterminePeaksAllChannels(ff, channels,
                                               breaks, remove_zeros, results,
                                               peak_removal, min_nr_bins_peakdetection)
    results <- all_peaks_res$results
    outlier_bins <- rep(TRUE, nrow(all_peaks_res$all_peaks))
    names(outlier_bins) <- seq_along(outlier_bins)

    # ------------------------ Isolation Tree  --------------------------------

    if (determine_good_cells == "all" || determine_good_cells == "IT"){
        if (length(res_breaks$breaks) >= force_IT ){
            IT_res <- isolationTreeSD(x = all_peaks_res$all_peaks,
                                      gain_limit=IT_limit)
            IT_cells <- RemovedBins(breaks, !IT_res$outlier_bins, nrow(ff))
            outlier_bins <- IT_res$outlier_bins
            results$IT <- IT_res$res
            results$OutlierIT <- IT_cells$cells
            results$ITPercentage <- (length(IT_cells$cell_ids)/nrow(ff))* 100
            message("IT analysis removed ", round(results$ITPercentage, 2),
                    "% of the measurements" )
        } else {
            warning(StrMessage("There are not enough bins for a robust isolation tree
                    analysis."))
            results$ITPercentage <- NA
        }
    }


    # ------------------------ Outliers based on mad --------------------------

    if (determine_good_cells == "all" || determine_good_cells == "MAD"){

        MAD_results <- MADOutlierMethod(all_peaks_res$all_peaks, outlier_bins,
                                        MAD, breaks, nrow(ff))
        mad_cells <- RemovedBins(breaks, MAD_results$MAD_bins, nrow(ff))
        results$OutlierMads <- mad_cells$cells
        results$ContributionMad <- MAD_results$Contribution_MAD
        results$MADPercentage <- (length(mad_cells$cell_ids)/nrow(ff))*100

        outlier_bins[outlier_bins] <- !MAD_results$MAD_bins

        message("MAD analysis removed ", round(results$MADPercentage, 2),
                "% of the measurements")
    }

    # ------------------------- indicate bad cells ----------------------------

    if (determine_good_cells %in% c("all", "IT", "MAD")){

        results <- RemoveShortRegions(ff, outlier_bins, consecutive_bins,
                                      breaks, results)

        message("The algorithm removed ",
                round(results$PercentageRemoved, 2),
                "% of the measurements")

        if(results$PercentageRemoved > 70){
            warning(StrMessage("More than 70% was removed from file ",
                               basename(filename)))
        }

        if(plot != FALSE &
           (results$PercentageRemoved >= plot | plot %in% c(TRUE, "all"))){
            plot <- TRUE
        } else {plot <- FALSE}
        new_ff <- ff[results$GoodCells, ]

        if (!("Original_ID" %in% colnames(new_ff@exprs))){
            new_ff <- AppendCellID(new_ff, which(results$GoodCells))
            results$FinalFF <- new_ff
        }

        # -----------------  Does the file need to be saved in an fcs? ---------
        if (save_fcs & !is.null(output_directory)){
            message("Saving fcs file")

            flowCore::write.FCS(new_ff,
                                file.path(fcs_directory,
                                          paste0(sub(".fcs",
                                                     "",
                                                     basename(filename)),
                                                 paste0(suffix_fcs, ".fcs"))))
        }

        # ----------------- Does an overview file need to be generated? --------
        if (report & !is.null(output_directory)){
            utils::write.table(t(c(
                basename(filename),
                nrow(ff),
                nrow(new_ff),
                results$PercentageRemoved,
                determine_good_cells, results$ITPercentage,
                results$MADPercentage, results$ConsecutiveCellsPercentage,
                MAD, IT_limit, consecutive_bins,
                events_per_bin, list_weird_channels$Changing_channel)),
                report_location, sep="\t", append=TRUE, row.names=FALSE,
                quote=FALSE, col.names=FALSE)
        }
    }

    #---------------- Does the file need to be plotted? ------------------------

    if(plot){
        if (!is.null(output_directory)){
            if(determine_good_cells %in% c("all", "IT", "MAD")){
                title_FR <- paste0(round(results$PercentageRemoved, 3),
                                   "% of the data was removed.")
            } else {
                title_FR <- ""
            }
            message("Plotting the results")
            PlotPeacoQC(ff=ff,
                        display_peaks=results,
                        output_directory=storing_directory,
                        channels=results$Channels,
                        title_FR=title_FR,
                        ...)
        }
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
#' PlotPeacoQC(ff, channels, output_directory=".", display_cells=2000,
#'             manual_cells=NULL, title_FR=NULL, display_peaks=TRUE,
#'             prefix="PeacoQC_", time_unit=100, time_channel_parameter="Time",
#'              ...)
#'
#' @param ff A flowframe
#' @param channels Indices of names of the channels in the flowframe that have
#' to be displayed
#' @param output_directory Directory where the plots should be generated. Set
#' to NULL if no plots need to be generated. The default is the working
#' directory.
#' @param display_cells The number of measurements that should be displayed.
#' (The number of dots that are displayed for every channel) The default is
#' 2000.
#' @param manual_cells Give a vector (TRUE/FALSE) with annotations for each cell
#' to compare the automated QC with. The default is NULL.
#' @param title_FR The title that has to be displayed above the flow rate
#' figure. Default is NULL.
#' @param display_peaks If the result of \code{PeacoQC} is given, all the
#' quality control results will be visualised. If set to TRUE: \code{PeacoQC}
#' will be run and only the peaks will be displayed without any quality control.
#' If set to FALSE, no peaks will be displayed and only the events will be
#' displayed. Default is TRUE.
#' @param prefix The prefix that will be given to the generated png file.
#' Default is "PeacoQC_".
#' @param time_unit The number of time units grouped together for visualising
#' event rate. The default is read from the $`TIMESTEP` keyword in the flowframe. 
#' If this keyword is not found, time_unit is set to 10.000. If the time plot has stacked lines, 
#' change this value.
#' @param time_channel_parameter Name of the time channel in ff if present.
#' Default is "Time".
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
#' fileName <- system.file("extdata", "111_Comp_Trans.fcs", package="PeacoQC")
#' ff <- flowCore::read.FCS(fileName)
#'
#' # Define channels on which the quality control should be done and the
#' # plots should be made
#' channels <- c(1, 3, 5:14, 18, 21)
#'
#' # Run PeacoQC
#' PeacoQC_res <- PeacoQC(ff,
#'     channels,
#'     determine_good_cells="all",
#'     plot=FALSE,
#'     save_fcs=TRUE)
#'
#' # Run PlotPeacoQC
#' PlotPeacoQC(ff, channels, display_peaks=PeacoQC_res)
#'
#' ## Plot only the peaks (No quality control)
#' PlotPeacoQC(ff, channels, display_peaks=TRUE)
#'
#' ## Plot only the dots of the file
#' PlotPeacoQC(ff, channels, display_peaks=FALSE)
#'
#' @export
PlotPeacoQC <- function(ff,
                        channels,
                        output_directory=".",
                        display_cells=2000,
                        manual_cells=NULL,
                        title_FR=NULL,
                        display_peaks=TRUE,
                        prefix="PeacoQC_",
                        time_unit=GetTimeUnit(ff),
                        time_channel_parameter="Time",
                        ...) {

    requireNamespace("ggplot2")


    input_res <- CheckInputPlot(ff, channels, output_directory,
                                display_cells, display_peaks, ...)
    display_cells <- input_res$display_cells
    peaks <- input_res$peaks
    plot_directory <- input_res$plot_directory


    # Initiate empty lists for all plots
    plot_list <- list()

    # Check for time channel
    if (!is.null(time_channel_parameter)){
        if (all(!grepl(time_channel_parameter,
                  colnames(flowCore::exprs(ff)),
                  ignore.case=TRUE))){
            time_channel <- NULL
        } else{
        time_channel <- grep(time_channel_parameter,
                             colnames(flowCore::exprs(ff)),
                             ignore.case=TRUE)}
    } else(time_channel <-  NULL)

    if (is.numeric(channels)){
        channels <- colnames(flowCore::exprs(ff))[channels]
    }

    # Plot acc score if one is listed in given arguments
    if (!is.null(title_FR)) {
        scores_time <- title_FR
    } else {
        scores_time <- ""
    }

    if ((length(flowCore::keyword(ff)$FILENAME) > 0) &&
        !is.na(flowCore::keyword(ff)$FILENAME)) {
        filename <- basename(flowCore::keyword(ff)$FILENAME)
    } else if ((length(flowCore::keyword(ff)[["$FIL"]]) > 0) &&
               !is.na(flowCore::keyword(ff)[["$FIL"]])) {
        filename <- basename(flowCore::keyword(ff)[["$FIL"]])
    }
    # Name to put on plotfile
    name <- sub(".fcs", "", filename)


    if (is(display_peaks, "list") && display_peaks$Analysis != FALSE){
        blocks <- MakeOverviewBlocks(ff, peaks, time_channel)
    } else {
        blocks <- NULL
    }


    if (!is.null(manual_cells)){
        manual_blocks <-  MakeManualBlocks(manual_cells)
    } else {
        manual_blocks <- NULL
    }

    if (!is.null(time_channel)){
        if (is(display_peaks, "list") && display_peaks$Analysis != FALSE){
            p_time <- BuildTimePlot(ff, blocks$overview_blocks_time,
                                    scores_time, time_unit,
                                    time_id = time_channel)
        } else{ p_time <- BuildTimePlot(ff, scores_time=scores_time,
                                        time_unit=time_unit,
                                        time_id = time_channel) }

        plot_list[["Time"]] <- p_time

    }

    plot_list <- BuildChannelPlots(channels, peaks, display_peaks,
                                   display_cells, manual_cells, manual_blocks,
                                   ff, blocks, plot_list)

    MakeNicePlots(display_peaks, plot_list, channels, plot_directory,
                  prefix, name)
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
#' PeacoQCHeatmap(report_location, show_values=TRUE, show_row_names=TRUE,
#' latest_tests=FALSE, title="PeacoQC report", ...)
#'
#' @param report_location The path to the PeacoQC report generated by
#' \code{PeacoQC}.
#' @param show_values If set to TRUE, the percentages of removed values
#' will be displayed on the heatmap. Default is TRUE.
#' @param show_row_names If set to FALSE, the filenames will not be displayed
#' on the heatmap. Default is TRUE.
#' @param latest_tests If this is set to TRUE, only the latest quality control
#' run will be displayed in the heatmap. Default is FALSE.
#' @param title The title that should be given to the heatmap. Default is
#' "PeacoQC_report".
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
#' location <- system.file("extdata", "PeacoQC_report.txt", package="PeacoQC")
#'
#' # Make heatmap overview of quality control run
#' PeacoQCHeatmap(report_location=location)
#'
#' # Make heatmap with only the runs of the last test
#' PeacoQCHeatmap(report_location=location, latest_tests=TRUE)
#'
#' # Make heatmap with row annotation
#'PeacoQCHeatmap(report_location=location,
#'               row_split=c(rep("r1",7), rep("r2", 55)))
#'
#' @importFrom grDevices colorRampPalette
#' @export

PeacoQCHeatmap <- function(
        report_location,
        show_values=TRUE,
        show_row_names=TRUE,
        latest_tests=FALSE,
        title="PeacoQC report",
        ...){

    if (!file.exists(report_location))
        stop(StrMessage("The path specified in the report_location parameter
            is wrong or incomplete."))

    if(show_row_names == FALSE & latest_tests == FALSE)
        warning(StrMessage("If there are duplicates in the report file,
            they will be displayed on the heatmap without their filename."))

    report_table <- utils::read.delim(report_location,
                                      check.names=FALSE,
                                      stringsAsFactors=FALSE)
    report_table[report_table == "Not_used"] <- NA

    if (latest_tests){
        report_table <- report_table[!duplicated(report_table$Filename,
                                                 fromLast=TRUE), ]
        rownames(report_table) <- report_table$Filename
    } else {

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
        "Consecutive bins"=factor(report_table$`Consecutive bins`),
        "IT limit"=factor(report_table$`IT limit`),
        "MAD"=factor(report_table$MAD),
        "Events per bin" = report_table$`Events per bin`,
        check.names=FALSE)


    rownames(annotation_frame) <- rownames(report_table)

    t1 <- colorRampPalette(c("#8D99AE", "#2B2D42"))
    col_cons <- t1(length(unique(annotation_frame$`Consecutive bins`)))
    t2 <- colorRampPalette(c("#EBB9DF", "#7D1D3F"))
    col_MAD <- t2(length(unique(annotation_frame$MAD)))
    t3 <- colorRampPalette(c("#B2CEDE", "#AD7A99"))
    col_IT <- t3(length(unique(annotation_frame$`IT limit`)))
    col_events <- colorRamp2(c(0, max(annotation_frame$`Events per bin`)/2,
                               max(annotation_frame$`Events per bin`)),
                             c("#5AAA95", "white", "#474973"))



    col_events <- colorRamp2(c(0, max(annotation_frame$`Events per bin`)),
                             c("#5AAA95", "#474973"))
    # t4 <- colorRampPalette(c("#5AAA95", "#474973"))
    # col_events <- t4(length(unique(annotation_frame$`Events per bin`)))

    names(col_cons) <- unique(annotation_frame$`Consecutive bins`)
    names(col_MAD) <- unique(annotation_frame$MAD)
    names(col_IT) <- unique(annotation_frame$`IT limit`)
    # names(col_events) <- unique(annotation_frame$`Events per bin`)


    analysis <- report_table$`Analysis by`

    if("all" %in% analysis | all(c("IT", "MAD") %in% analysis)){
        ha <- rowAnnotation(df=annotation_frame,
                            col=list("Consecutive bins"=col_cons,
                                     "MAD"=col_MAD,
                                     "IT limit"=col_IT,
                                     "Events per bin"=col_events),
                            annotation_legend_param =
                                list("IT limit" = list(ncol = 1),
                                     "MAD" = list(ncol = 1),
                                     "Consecutive bins" = list(ncol = 1),
                                     "Events per bin" = list(direction = "horizontal")))

    } else if(length(unique(analysis)) == 1 & unique(analysis) == "IT"){
        ha <- rowAnnotation(
            "Consecutive bins"=annotation_frame$`Consecutive bins`,
            "IT limit"=annotation_frame$`IT limit`,
            col=list("Consecutive bins"=col_cons,
                     "IT limit"=col_IT,
                     "Events per bin"=col_events),
            annotation_legend_param = list(
                "IT limit" = list(ncol = 1),
                "Consecutive bins" = list(ncol = 1),
                "Events per bin" = list(direction = "horizontal")))
    } else if(length(unique(analysis)) == 1 & unique(analysis) == "MAD"){
        ha <- rowAnnotation(
            "Consecutive bins"=annotation_frame$`Consecutive bins`,
            "MAD"=annotation_frame$MAD,
            col=list("Consecutive bins"=col_cons,
                     "MAD"=col_MAD,
                     "Events per bin"=col_events),
            annotation_legend_param = list(
                "MAD" = list(ncol = 1),
                "Consecutive bins" = list(ncol = 1),
                "Events per bin" = list(direction = "horizontal")))

    }

    col_incr_decr_channel <- c()

    if("No increasing or decreasing effect" %in%
       report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
                                   "No increasing or decreasing effect"="#26C485")
    }
    if ("Increasing channel" %in% report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
                                   "Increasing channel"="#AF3800")
    }
    if ("Decreasing channel" %in% report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
                                   "Decreasing channel"="#721817")
    }
    if ("Increasing and decreasing channel" %in%
        report_table$`Increasing/Decreasing channel`){
        col_incr_decr_channel <- c(col_incr_decr_channel,
                                   "Increasing and decreasing channel"="#A50104")
    }

    annotation_right <- rowAnnotation(
        "Incr/Decr"=as.factor(report_table[,
                                           "Increasing/Decreasing channel"]),
        col=list("Incr/Decr"=col_incr_decr_channel))

    report_matrix <- data.matrix(report_table[, c(4, 6, 7, 8)])

    if(show_values){
        cell_fun=function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.1f", report_matrix[i, j]), x, y,
                      gp=gpar(fontsize=10))
        }
    } else{
        cell_fun=NULL
    }

    ph <- Heatmap(report_matrix,
                  cell_fun=cell_fun,
                  cluster_columns=FALSE,
                  cluster_rows=FALSE,
                  column_title=title,
                  column_title_gp=grid::gpar(fontface="bold"),
                  col=circlize::colorRamp2(c(0, 20, 100), c("#EBEBD3", "#FFD151", "red")),
                  left_annotation=ha,
                  right_annotation=annotation_right,
                  heatmap_legend_param=list(direction="horizontal"),
                  name="Removed percentage",
                  na_col="Grey",
                  show_row_names=show_row_names,
                  ...)

    draw(ph, annotation_legend_side="bottom", heatmap_legend_side="bottom")
}

