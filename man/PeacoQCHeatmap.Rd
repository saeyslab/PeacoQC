% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeacoQC.R
\name{PeacoQCHeatmap}
\alias{PeacoQCHeatmap}
\title{Make overview heatmap of quality control analysis}
\usage{
PeacoQCHeatmap(report_location, show_values=TRUE, show_row_names=TRUE,
latest_tests=FALSE, title="PeacoQC report", ...)
}
\arguments{
\item{report_location}{The path to the PeacoQC report generated by
\code{PeacoQC}.}

\item{show_values}{If set to TRUE, the percentages of removed values
will be displayed on the heatmap. Default is TRUE.}

\item{show_row_names}{If set to FALSE, the filenames will not be displayed
on the heatmap. Default is TRUE.}

\item{latest_tests}{If this is set to TRUE, only the latest quality control
run will be displayed in the heatmap. Default is FALSE.}

\item{title}{The title that should be given to the heatmap. Default is
"PeacoQC_report".}

\item{...}{Extra parameters to be given to the \code{Heatmap} function
(eg. row_split)}
}
\value{
This function returns nothing but generates a heatmap that can be
saved as pdf or png
}
\description{
\code{PeacoQCHeatmap} will make a heatmap to display all the
results generated by \code{PeacoQC}. It will include the percentages of
measurements that are removed in total, by the IT method and by the MAD
method. It will also show the parameters that were used during the
quality control.
}
\examples{

# Find path to PeacoQC report
location <- system.file("extdata", "PeacoQC_report.txt", package="PeacoQC")

# Make heatmap overview of quality control run
PeacoQCHeatmap(report_location=location)

# Make heatmap with only the runs of the last test
PeacoQCHeatmap(report_location=location, latest_tests=TRUE)

# Make heatmap with row annotation
PeacoQCHeatmap(report_location=location,
              row_split=c(rep("r1",7), rep("r2", 55)))

}
