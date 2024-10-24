% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Extractor functions MA.R
\name{process_array}
\alias{process_array}
\title{Exception handling for using slabR}
\usage{
process_array(array, scheme, count, collapse_days)
}
\arguments{
\item{scheme}{a summary scheme as expected by `array_w_mean()`.}

\item{count}{an optional vector of steps to subset along. See ncdf4 documentation.}

\item{collapse_days}{a logical switch setting whether to extract all days within a single netcdf file, or to take the mean across days.}

\item{file}{the name of a netcdf file containing the title variables.}
}
\value{
A matrix with a column of group averages per title variable. If multiple days have been extracted these run sequentially within a column.
}
\description{
This function catches the special cases when arrays are not 3D, which causes errors for Armadillo in SlabR.
}
\details{
The function checks how many depth and time steps there are in an array of NEMO-ERSEM output. The function
then ensures that slabR is passed a 3D array, either by collapsing or applying the summary scheme over the time dimension,
or if there are too few dimensions (because a depth layer is being extracted) a simple join is performed to the scheme.
}
\seealso{
Other NEMO-ERSEM variable extractors: 
\code{\link{NEMO_ERSEM}()},
\code{\link{extractors_1D}},
\code{\link{extractors_MA}},
\code{\link{extractors_slabR}},
\code{\link{get_spatial}()}
}
\concept{NEMO-ERSEM variable extractors}
