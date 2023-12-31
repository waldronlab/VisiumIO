% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TENxSpatialCSV.R
\docType{class}
\name{TENxSpatialCSV-class}
\alias{TENxSpatialCSV-class}
\alias{.TENxSpatialCSV}
\alias{TENxSpatialCSV}
\alias{import,TENxSpatialCSV,ANY,ANY-method}
\title{Represent and import spatial CSV data from 10X Genomics}
\usage{
TENxSpatialCSV(resource, colnames = .TISSUE_POS_COLS)

\S4method{import}{TENxSpatialCSV,ANY,ANY}(con, format, text, ...)
}
\arguments{
\item{resource}{character(1) The path to the file}

\item{colnames}{\code{character()} A vector specifying the column names of the
CSV, defaults to \code{c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")}.}

\item{con}{The connection from which data is loaded or to which data is
saved. If this is a \code{character} vector, it is assumed to be a file name and
a corresponding file connection is created and then closed after exporting
the object. If it is a \linkS4class{BiocFile} derivative, the data is loaded from
or saved to the underlying resource. If missing, the function will return
the output as a character vector, rather than writing to a connection.}

\item{format}{The format of the output. If missing and \code{con} is a
file name, the format is derived from the file extension. This argument is
unnecessary when \code{con} is a derivative of \linkS4class{BiocFile}.}

\item{text}{If \code{con} is missing, this can be a character vector directly
providing the string data to import.}

\item{...}{Additional inputs to the low level class generator functions}
}
\value{
TENxSpatialCSV: An object of class \link{TENxSpatialCSV}

import-method: A \code{DataFrame} object containing the data from the CSV
file
}
\description{
\code{TENxSpatialCSV} is a class to represent and import spatial CSV
files with specific column names. It is a composed class of
\link[TENxIO:TENxFile]{TENxIO::TENxFile} and contains additional slots for the column names and
whether the CSV is a list-type of file.
}
\details{
Typically, the user will not create an object of this class directly
but rather use the \code{\link[=TENxVisium]{TENxVisium()}} constructor function to create an object
of this class in the background. The column names are set to the default
values of \code{c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")}. The column names can be
changed by specifying the \code{colnames} argument in the constructor function.
}
\section{Slots}{

\describe{
\item{\code{isList}}{\code{logical(1)} A scalar specifying whether the CSV is a list-type
of file}

\item{\code{colnames}}{\code{character()} A vector specifying the column names of the CSV}
}}

\examples{
sample_dir <- system.file(
    file.path("extdata", "10xVisium", "section1"),
    package = "SpatialExperiment"
)
spatial_dir <- Filter(
  function(x) endsWith(x, "spatial"), list.dirs(sample_dir)
)
csvresource <- file.path(spatial_dir, "tissue_positions_list.csv")
TENxSpatialCSV(csvresource)
head(import(TENxSpatialCSV(csvresource)), 4)

}
