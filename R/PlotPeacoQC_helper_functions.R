MakeOverviewBlocks <- function(ff, peaks, time_channel){
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

    fill_blocks <- ifelse(run_length$values == TRUE, "Good Values",
                        ifelse(run_length$values == "consecutive",
                                "In consecutive bins",
                                ifelse(run_length$values == FALSE, "IT",
                                        "MAD")))

    x_min <- c(1, cumsum(run_length$lengths)[-length(run_length$lengths)])
    x_max <- cumsum(run_length$lengths)

    overview_blocks <- data.frame(x_min=x_min,
                                    x_max=x_max,
                                    y_min=-Inf,
                                    y_max=Inf,
                                    fill_blocks=fill_blocks)

    if(!is.null(time_channel)){
    overview_blocks_time <- NULL
    if (length(time_channel) > 0){

        x_min <- flowCore::exprs(ff)[, time_channel][c(1,
            cumsum(run_length$lengths)[-length(run_length$lengths)])]
        x_max <- flowCore::exprs(ff)[, time_channel][cumsum(run_length$lengths)]

        overview_blocks_time <- data.frame(x_min=x_min,
                                            x_max=x_max,
                                            y_min=-Inf,
                                            y_max=Inf,
                                            fill_blocks=fill_blocks)
    }} else{overview_blocks_time <- NULL}

    return(list("Overview_blocks"=overview_blocks,
                "overview_blocks_time"=overview_blocks_time))
}

MakeManualBlocks <- function(manual_cells){
    run_length_man <- rle(manual_cells)


    fill_blocks_man <- ifelse(
        run_length_man$values,
        "deepskyblue1",
        "indianred1")

    x_min_man <- c(1,
                    cumsum(run_length_man$lengths)
                    [-length(run_length_man$lengths)])
    x_max_man <- cumsum(run_length_man$lengths)

    manual_blocks <- data.frame(x_min=x_min_man,
                                x_max=x_max_man,
                                y_min=-Inf,
                                y_max=Inf,
                                fill_blocks=fill_blocks_man)
    return(manual_blocks)
}


BuildTimePlot <- function(ff, blocks=NULL, scores_time, time_unit=100, time_id){

    h <- graphics::hist(flowCore::exprs(ff)[ seq_len(nrow(ff)), time_id],
                        breaks=seq(min(flowCore::exprs(ff)[, time_id]),
                                        max(flowCore::exprs(ff)[, time_id]) +
                                        time_unit, by=time_unit),
                        plot=FALSE)

    idcs <- findInterval(flowCore::exprs(ff)[ seq_len(nrow(ff)), time_id], h$breaks)

    p_time <- ggplot() + theme_bw()

    p_time <- p_time + theme(panel.grid=element_blank())

    if (!is.null(blocks)){
        p_time <- p_time + geom_rect(data=blocks,
                                        mapping=aes(xmin=x_min,
                                                    xmax=x_max,
                                                    ymin=y_min,
                                                    ymax=y_max,
                                                    fill=fill_blocks),
                                        alpha=0.4,
                                        show.legend=TRUE) +
                            scale_fill_manual(name="",
                            values=c(IT="indianred1",
                                        MAD="mediumpurple1",
                                        `In consecutive bins`="plum1",
                                        `Good Values`="white"),
                            guide=guide_legend(override.aes=
                                                list(alpha=0.4)))
        p_time <- p_time + theme(legend.key=element_rect(colour="snow4"))

    }

    p_time <- p_time +
        geom_point(aes(x=h$mids, y=h$counts)) +
        ggtitle(scores_time) +
        xlab("Time") + ylab("Nr of cells per second") +
        theme(plot.title=element_text(size=10))

    return(p_time)

}


MakeContributions <- function(peaks, channel, marker){

    contributions <- marker

    contribution_MAD_channel <- max(peaks$ContributionMad
                            [grep(channel,
                                    names(peaks$ContributionMad))])
    if (exists("IT", peaks)){
        contribution_IT_channel <- ifelse(length(grep(channel,
                                        peaks$IT$split_column)) > 0, "+", "/")
    } else{ contribution_IT_channel <- "/"}

    if (contribution_MAD_channel > 0 || contribution_IT_channel != "/" ){
        contributions <- paste0(contributions, "\n")

        if (contribution_IT_channel != "/"){
            contributions <- paste0(contributions, "IT: +")}
        if (contribution_MAD_channel > 0 & contribution_IT_channel != "/"){
            contributions <- paste0(contributions, " ")
        }
        if (contribution_MAD_channel > 0){
            contributions <- paste0(contributions, "MAD: ",
                                    contribution_MAD_channel, "%")
            }
    }

    if (length(grep(channel, peaks$WeirdChannels$Increasing)) > 0) {
        contributions <- paste0(contributions, "\n",
                                "WARNING: Increasing channel.")
    } else if (length(grep(channel,
                            peaks$WeirdChannels$Decreasing))> 0){
        contributions <- paste0(contributions, "\n",
                                "WARNING: Decreasing channel.")

    }

    return(contributions)

}

BuildManualPartPlot <- function(p, channel, manual_blocks){

    if (length(grep("SC", channel)) == 0) {
        manual_blocks$y_min <- minimum - 1
        manual_blocks$y_max <- minimum - 1.25
    } else {
        manual_blocks$y_min <- minimum - 45000
        manual_blocks$y_max <- minimum - 58000
    }
    p <- p + geom_rect(data=manual_blocks,
                        mapping=aes(xmin=x_min,
                                    xmax=x_max,
                                    ymin=y_min,
                                    ymax=y_max),
                        fill=manual_blocks$fill_blocks,
                        alpha=1)

    p <- p + theme(legend.position="none")

    return(p)
}


BuildPeaksPlot <- function(p, peaks, channel, mid_breaks, contributions){
    peak_frame <- peaks[[channel]]

    if (is(peak_frame, "data.frame")) {

        peak_frame$Bin <- as.numeric(mid_breaks)[as.numeric(
            peak_frame$Bin)]

        colours <- paste0("grey",
                            seq_len(20))[seq_len(length(
                                unique(peak_frame[,"Cluster"])))]

        p <- p + geom_line(data=peak_frame, aes(x=Bin,
                                                y=Peak,
                                                color=Cluster),
                            linewidth=1,
                            show.legend=FALSE) +
            scale_color_manual(values=colours)
    } else {
        contributions <- paste0(contributions, " No peak was found.")
    }

    return(list("plot"=p, "contributions"=contributions))
}


BuildBackgroundQCPlot <- function(p, overview_blocks){

    p <- p + geom_rect(data=overview_blocks,
                        mapping=aes(xmin=x_min,
                                    xmax=x_max,
                                    ymin=y_min,
                                    ymax=y_max,
                                    fill=fill_blocks),
                        alpha=0.4,
                        show.legend=TRUE) +
        scale_fill_manual(name="",
                            values=c(IT="indianred1",
                                    MAD="mediumpurple1",
                                    `In consecutive bins`="plum1",
                                    `Good Values`="white"),
                            guide=guide_legend(override.aes=
                                                list(alpha=0.4)))
    p <- p + theme(legend.key=element_rect(colour="snow4"))

    return(p)
}

MakeMidBreaks <- function(peaks, ncells){
    events_per_bin <- peaks$EventsPerBin

    mid_breaks <- SplitWithOverlapMids(seq_len(ncells),
                                        events_per_bin,
                                        ceiling(events_per_bin/2))
    m <-  ncells%%events_per_bin

    if(m !=0) {
        mid_breaks=c(mid_breaks, ncells)
    }

    return(mid_breaks)
}

BuildChannelPlots <- function(channels, peaks, display_peaks, display_cells,
                                manual_cells, manual_blocks,  ff, blocks,
                                plot_list){
    if(is(peaks, "list")){
        mid_breaks <- MakeMidBreaks(peaks, nrow(ff))
    }

    subset_signalplot <- sort(sample(seq_len(nrow(ff)), display_cells))

    for (channel in channels) {

        marker <- flowCore::getChannelMarker(ff, channel)$desc
        if(is.na(marker))
            marker <- channel

        minimum <- min(flowCore::exprs(ff)[, channel])
        maximum <- max(flowCore::exprs(ff)[, channel])
        range <- abs(minimum) + abs(maximum)

        if (is(display_peaks, "list") && display_peaks$Analysis != FALSE){
            contributions <- MakeContributions(peaks, channel, marker)
        } else{contributions <- marker}


        flowdata <- data.frame(Cells=subset_signalplot,
                                channel=flowCore::exprs(ff)
                                [subset_signalplot, channel])

        p <- ggplot() +
            ylab("Value") +
            xlab("Cells") +
            theme_bw() +
            theme(plot.title=element_text(hjust=0),
                    panel.grid=element_blank())

        if (is(display_peaks, "list") && display_peaks$Analysis != FALSE) {
            p <- BuildBackgroundQCPlot(p, blocks$Overview_blocks)
        }

        p <- p + geom_point(data=flowdata,
                            aes(x=Cells, y=channel),
                            size=0.3,
                            col="snow4")

        if(is(peaks, "list")){
            peaks_plot <- BuildPeaksPlot(p, peaks, channel,
                                            mid_breaks, contributions)
            p <- peaks_plot$plot
            contributions <- peaks_plot$contributions
        }

        p <- p + ggtitle(contributions)

        if (!is.null(manual_cells)) {
            p <- BuildManualPartPlot(p, manual_blocks, channel)
        }

        p <- p + theme(plot.subtitle=element_text(size=10, color="black"))

        plot_list[[channel]] <- p
    }

    return(plot_list)
}


CheckInputPlot <- function(ff, channels, output_directory,
                            display_cells, display_peaks, ...){
    if(!is(ff, "flowFrame") | is.null(ff))
        stop("ff should be a flowFrame.")
    if(!is.numeric(channels) & !all(channels%in% colnames(flowCore::exprs(ff)))|
    is.null(channels))
        stop(StrMessage("Make sure that you use indices or the colnames in the
            expression matrix in the flowframe to indicate which channels you
            want to use."))
    if(is.null(output_directory))
        stop("There should be a path given to the output_directory parameter.")

    # Make a new directory where all results will be stored
    if(!is.null(output_directory)){

        plot_directory <- file.path(output_directory, "PeacoQC_plots")
        suppressWarnings(dir.create(plot_directory))
    }

    if (display_cells > nrow(ff)) {
        message(StrMessage("There are less then the number of display cells
            available. Setting the number of display cells to the number of
            measurements."))
        display_cells <- nrow(ff)
    }

    # If display_peaks == TRUE, the peaks should be calculated by using PeacoQC
    if (is(display_peaks, "logical")){
        if(display_peaks){
            message("Running PeacoQC to determine peaks")
            peaks <- PeacoQC(ff,
                            channels,
                            determine_good_cells=FALSE,
                            output_directory=NULL,
                            plot=FALSE,
                            save_fcs=FALSE,
                            report=FALSE,
                            ...)
        } else{peaks <- FALSE}
    } else { peaks <- display_peaks}


    return(list("display_cells"=display_cells,
                "peaks"=peaks,
                "plot_directory"=plot_directory))
}

MakeNicePlots <- function(display_peaks, plot_list, channels, plot_directory,
                            prefix, name){
    n_channels <- length(channels)

    n_row <- floor(sqrt(n_channels + 2))
    n_col <- ceiling((n_channels + 2)/n_row)

    g <- ggplotGrob(plot_list[[1]] + theme(legend.position="bottom"))$grobs

    if (!(is(display_peaks, "logical")) && display_peaks$Analysis != FALSE){
        legend <- g[[which(vapply(g,
                                    function(x) x$name,
                                    FUN.VALUE=character(1)) == "guide-box")]]
        lheight <- sum(legend$height)
        plots <- do.call(gridExtra::arrangeGrob,
                            c(lapply(plot_list, function(x)x +
                                        theme(legend.position="none")),
                            nrow=n_row,
                            ncol=n_col))
        final_plots <- gridExtra::arrangeGrob(
            plots,
            legend,
            heights=grid::unit.c(grid::unit(1, "npc") - lheight, lheight))
    } else {
        plots <- do.call(gridExtra::arrangeGrob,
                            c(lapply(plot_list,
                                        function(x) x),
                                nrow=n_row, ncol=n_col))
        final_plots <- gridExtra::arrangeGrob(plots)}

    ggsave(paste0(plot_directory, "/", prefix, name, ".png"),
            final_plots, width=n_col * 5,
            height=n_row * 3, limitsize=FALSE)
}

GetTimeUnit <- function(ff){
  time_unit <- 10000
  
  if(!is.null(ff@description$`$TIMESTEP`)){
    unit <- as.numeric(ff@description$`$TIMESTEP`)
    time_unit<- 1/unit
  } else {
    warning('`$TIMESTEP` value not found in the FlowFrame, running PlotPeacoQC with deafult `time_unit`=10000',
            call.=FALSE)
  }
  return(time_unit)
}