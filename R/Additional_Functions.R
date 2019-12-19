# ------------------- Remove margins -------------------------------------------

# Function to remove events on the margins
removeMargins <- function(flowFrame, dimensions) {
    # Look up the accepted ranges for the dimensions
    
    meta <- pData(flowFrame@parameters)
    rownames(meta) <- meta[, "name"]
    
    
    
    # Initialize variables
    selection <- rep(TRUE, times = dim(flowFrame)[1])
    e <- flowFrame@exprs
    
    
    # Make selection
    for (d in dimensions) {
        ## first: the min between the minRange and 0 (minRange could be wrong), than the max between this and the min found in the expression matrix
        
        selection <- selection & e[, d] > max(min(meta[d, "minRange"], 0), min(e[, d])) & e[, d] < min(meta[d, "maxRange"], max(e[, d]))
    }
    
    
    if (length(which(selection == FALSE))/length(selection) > 0.1) {
        warning(paste0("More then ", round(length(which(selection == FALSE))/length(selection) * 100, 2), "% is considered as margin event. This should be verified."))
    }
    
    return(selection)
}


# ------------------------- Calculate CV --------------------------------------

CV_peaks <- function(ff, kept_cells) {
    ff_filtered <- ff[kept_cells, ]
    peaks <- PeacoQC(ff_filtered, determine_good_cells = FALSE, channels = channels)$AllPeaks
    cv <- apply(peaks, 2, function(x) {
        (sd(x)/mean(x)) * 100
    })
}
