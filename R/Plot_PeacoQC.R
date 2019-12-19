
# --------------------- Plot the peaks of the PeacoQC function -----------------

PlotPeacoQC <- function(ff, output_directory,
                        channels,
                        display_cells = 5000,
                        manual_cells = NULL,
                        acc_score = NULL,
                        show_peaks = TRUE,
                        prefix = "PeacoQC_") {


    # Make folder

    suppressWarnings(dir.create(output_directory))

    # Plot acc score if one is listed in given arguments
    if (!is.null(acc_score)) {
        scores_time <- acc_score
    } else {
        scores_time <- ""
    }

    filename <- basename(ff@description$FILENAME)

    # Name to put on plotfile
    name <- sub(".fcs", "", filename)

    n_channels <- length(channels)
    if (missing(channels)) {
        time_index <- grep("Time", colnames(ff@exprs))
        channels <- colnames(ff@exprs)[-time_index]

    } else if (is.numeric(channels)) {
        channels <- colnames(ff@exprs)[channels]
    }


    # Determining grid to plot
    n_row <- floor(sqrt(n_channels + 2))
    n_col <- ceiling((n_channels + 2)/n_row)
    set.seed(1)

    # Subset of points that will be plotted
    if (nrow(ff) >= 50000) {
        subset_timeplot <- sort(sample(seq_len(nrow(ff)), 50000))
    } else {
        subset_timeplot <- seq_len(nrow(ff))
    }

    if (display_cells > nrow(ff)) {
        print("There are less then the number of display cells available.
              Setting the number of display cells to the number of measurements.")
        display_cells <- nrow(ff)
    }

    subset_signalplot <- sort(sample(seq_len(nrow(ff)), display_cells))

    # Calculating time breaks
    breaks <- cut(ff@exprs[, "Time"],
                  breaks = seq(0, max(ff@exprs[, "Time"]) + 100, by = 100),
                  labels = FALSE)
    mid_breaks <- seq(50, max(ff@exprs[, "Time"]) + 50, by = 100)


    h <- hist(ff@exprs[subset_timeplot, "Time"],
              breaks = seq(0, max(ff@exprs[, "Time"]) + 100, by = 100),
              plot = FALSE)

    idcs <- findInterval(ff@exprs[subset_timeplot, "Time"], h$breaks)

    h_2 <- hist(ff@exprs[, "Time"],
                breaks = seq(0, max(ff@exprs[, "Time"]) + 100, by = 100),
                plot = FALSE)

    idcs_2 <- findInterval(ff@exprs[, "Time"], h$breaks)




    # Calculate backgroundvalues for points on plot, rectangle block and CV values after automated qc

    if (!is.null(show_peaks$GoodCells)) {


        # Make blocks for automated gated algorithms to display on plot
        run_length <- rle(show_peaks$GoodCells)

        full_QC_vector <- ifelse(show_peaks$GoodCells == TRUE, TRUE, FALSE)


        consecutive_cells <- which(show_peaks$ConsecutiveCells == FALSE)
        if (length(consecutive_cells) > 0) {
            full_QC_vector[!show_peaks$ConsecutiveCells] <- "consecutive"
            run_length <- rle(full_QC_vector)

        }

        mad_cells <- which(show_peaks$OutlierMads == FALSE)
        if (length(mad_cells) > 0) {
            full_QC_vector[!show_peaks$OutlierMads] <- "mad"
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

        x_min <- ff@exprs[, "Time"][c(1, cumsum(run_length$lengths)[-length(run_length$lengths)])]
        x_max <- ff@exprs[, "Time"][cumsum(run_length$lengths)]

        overview_blocks_background <- data.frame(x_min = x_min,
                                                 x_max = x_max,
                                                 y_min = -Inf,
                                                 y_max = Inf,
                                                 fill_blocks = fill_blocks)




    }

    # Make blocks for manual gates to display on signalplot + cv score after manual gating flowrate
    if (!is.null(manual_cells)) {
        run_length_man <- rle(manual_cells)


        fill_blocks_man <- ifelse(run_length_man$values, "deepskyblue1", "indianred1")

        x_min_man <- c(1, cumsum(run_length_man$lengths)[-length(run_length_man$lengths)])
        x_max_man <- cumsum(run_length_man$lengths)

        manual_blocks <- data.frame(x_min = x_min_man,
                                    x_max = x_max_man,
                                    y_min = -Inf,
                                    y_max = Inf,
                                    fill_blocks = fill_blocks_man)
    }


    # Initiate empty lists for all plots
    plot_list <- list()



    # Build time plot (FlowRate)

    p_time <- ggplot() + theme_bw()

    p_time <- p_time + theme(panel.grid = element_blank())



    if (!is.null(show_peaks$GoodCells)) {

        p_time <- p_time + geom_rect(data = overview_blocks_background,
                                     mapping = aes(xmin = x_min,
                                                   xmax = x_max,
                                                   ymin = y_min,
                                                   ymax = y_max,
                                                   fill = fill_blocks),
                                     alpha = 0.4,
            show.legend = T) + scale_fill_manual(name = "",
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



    # -------------------------------- Individual channels -----------------------


    indices_matched <- show_peaks$MidBreaks

    channel <- channels[2]


    for (channel in channels) {

        marker <- FlowSOM::get_markers(ff, channel)

        minimum <- min(ff@exprs[, channel])
        maximum <- max(ff@exprs[, channel])
        range <- abs(minimum) + abs(maximum)


        if (!is.null(show_peaks$GoodCells)) {

            # Show contributions of every channel in MAD and IF

            contribution_MAD <- sum(show_peaks$ContributionMad[grep(channel,
                                                                    names(show_peaks$ContributionMad))])
            contribution_IT <- ifelse(length(grep(channel, show_peaks$IT$res$split_column)) > 0, "+", "/")

            contributions <- paste0(marker, "\n", "IT: ", contribution_IT, " MAD: ", contribution_MAD, "%")

        } else (contributions <- marker)


        if (length(grep(channel, show_peaks$WeirdChannels$Increasing)) > 0) {
            contributions <- paste0(contributions, "\n", "WARNING: Increasing channel.")
        } else if (length(grep(channel, show_peaks$WeirdChannels$Decreasing)) > 0) {
            contributions <- paste0(contributions, "\n", "WARNING: Decreasing channel.")
        }

        flowdata <- data.frame(Cells = subset_signalplot, channel = ff@exprs[subset_signalplot, channel])


        # Initial plot
        p <- ggplot() + ylab("Value") + xlab("Cells") + theme_bw() + theme(plot.title = element_text(hjust = 0), panel.grid = element_blank())


        if (!is.null(show_peaks$GoodCells)) {

            p <- p + geom_rect(data = overview_blocks, mapping = aes(xmin = x_min,
                                                                     xmax = x_max,
                                                                     ymin = y_min,
                                                                     ymax = y_max,
                                                                     fill = fill_blocks),
                               alpha = 0.4,
                               show.legend = T) +
                scale_fill_manual(name = "",
                                  values = c(IT = "indianred1",
                                             MAD = "mediumpurple1",
                                             `In consecutive bins` = "plum1",
                                             `Good Values` = "white"),
                                  guide = guide_legend(override.aes = list(alpha = 0.4)))
            p <- p + theme(legend.key = element_rect(colour = "snow4"))
            p <- p + geom_point(data = flowdata,
                                aes(x = Cells, y = channel),
                                col = "snow4", size = 0.3)

        } else {
            p <- p + geom_point(data = flowdata,
                                aes(x = Cells, y = channel),
                                size = 0.3,
                                col = "snow4")
        }

        peak_frame <- show_peaks[[channel]]

        if (class(peak_frame) == "data.frame") {

            length_bins <- sapply(show_peaks$Breaks, length)

            peak_frame$Bin <- as.numeric(indices_matched)[as.numeric(peak_frame$Bin)]

            colors <- colorRampPalette(c("grey0", "grey15"))

            p <- p + geom_line(data = peak_frame, aes(x = Bin,
                                                      y = Peak,
                                                      color = Cluster),
                               size = 1,
                               show.legend = F) +
              scale_color_manual(values = colors(length(unique(peak_frame$Cluster))))
        } else {
            contributions <- paste0(contributions, " No peak was found!")
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
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    plots <- do.call(arrangeGrob, lapply(plot_list, function(x) x + theme(legend.position = "none")))
    new <- arrangeGrob(plots, legend, heights = unit.c(unit(1, "npc") - lheight, lheight))


    ggsave(paste0(output_directory, "/", prefix, name, ".png"),
           new, width = n_col * 5, height = n_row * 3, limitsize = F)

}
