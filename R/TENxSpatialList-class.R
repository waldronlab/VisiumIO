#' @docType class
#'
#' @title A class to represent and import spatial Visium data
#'
#' @description This class is a composed class of [TENxFileList], which can
#'   contain a list of [TENxFile] objects, and a [TENxSpatialList] object. It is
#'   meant to handle spatial Visium data from 10X Genomics.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class.
#'
#' @inheritParams TENxVisium
#'
#' @slot images `character()` The image names to use with `grep` and include in
#'   the list of files.
#'
#' @slot scaleJSON `character(1)` The file name of the scale factors JSON file,
#'   defaults to 'scalefactors_json.json'.
#'
#' @slot tissuePos `character(1)` The file name of the tissue positions file;
#'   typically a `.parquet` or `.csv` file.
#'
#' @slot sampleId `character(1)` A scalar specifying the sample identifier.
#'
#' @slot binSize The bin size of the images to import. The default slot value is
#'   `character()`. It typically corresponds to the directory name
#'   `square_000um` where `000` is the bin value.
#'
#' @exportClass TENxSpatialList
.TENxSpatialList <- setClass(
    "TENxSpatialList",
    contains = "TENxFileList",
    slots = c(
        images = "character",
        scaleJSON = "character",
        tissuePos = "character",
        sampleId = "character",
        binSize = "character"
    )
)

.check_file_pattern <- function(obj, pattern) {
    fname <- switch(
        pattern,
        "tissue_positions.*" = "tissue positions",
        "scalefactors.*\\.json$" = "scalefactor JSON"
    )
    if (!any(grepl(pattern, names(obj))))
        paste0("The '", fname, "' file was not found")
}

.check_file <- function(obj, filename) {
    if (!filename %in% names(obj))
        paste0("The '", filename, "' file was not found")
}

.validTENxSpatialList <- function(object) {
    c(
        .check_file_pattern(object, "tissue_positions.*"),
        .check_file_pattern(object, "scalefactors.*\\.json$"),
        .check_file(object, object@scaleJSON)
    )
}

S4Vectors::setValidity2("TENxSpatialList", .validTENxSpatialList)

.SCALE_JSON_FILE <- "scalefactors_json.json"

#' @rdname TENxSpatialList-class
#'
#' @inheritParams TENxVisium
#'
#' @param bin_size `character(1)` The bin size of the images to import. The
#'   default is `008`. It corresponds to the directory name `square_000um` where
#'   `000` is the bin value.
#'
#' @importFrom BiocIO decompress
#'
#' @returns A `SpatialExperiment` object
#'
#' @examples
#' spatial_dir <- system.file(
#'     file.path("extdata", "10xVisium", "section1", "outs", "spatial"),
#'     package = "SpatialExperiment"
#' )
#'
#' TENxSpatialList(resources = spatial_dir, images = "lowres")
#'
#' TENxSpatialList(resources = spatial_dir, images = "lowres") |>
#'     metadata() |> lapply(names)
#'
#' @export
TENxSpatialList <- function(
    resources,
    sample_id = "sample01",
    images = c("lowres", "hires", "detected", "aligned", "aligned_fiducials"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*",
    bin_size = character(0L),
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    if (!is(resources, "TENxFileList"))
        resources <- TENxFileList(resources, ...)
    if (resources@compressed)
        resources <- decompress(con = resources)
    tissuePos <- grep(tissuePattern, names(resources), value = TRUE)
    if (!length(tissuePos))
        stop("No tissue positions file found with pattern: ", tissuePattern)

    if (missing(bin_size) && any(grepl("square_\\d{3}", path(resources)))) {
        bin_size <- unique(
            gsub(".*?square_(\\d{3}).*", "\\1", path(resources))
        )
        if (!identical(length(bin_size), 1L))
            stop("Multiple 'bin_size' values found in the directory")
    }

    .TENxSpatialList(
        resources, images = images, scaleJSON = jsonFile,
        tissuePos = tissuePos, sampleId = sample_id,
        binSize = bin_size
    )
}

#' @describeIn TENxSpatialList Import a `TENxSpatialList` object
#'
#' @param ... Parameters to pass to the format-specific method.
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocIO FileForFormat
#'
#' @exportMethod import
setMethod("import", "TENxSpatialList", function(con, format, text, ...) {
    jsonFile <- con@scaleJSON
    sampid <- con@sampleId
    sfs <- jsonlite::fromJSON(txt = path(con)[jsonFile])

    DFs <- lapply(con@images, function(image) {
        .getImgRow(con = con, sampleId = sampid, image = image, scaleFx = sfs)
    })
    fff <- FileForFormat(
        path(con)[con@tissuePos],
        prefix = "TENxSpatial", suffix = NULL
    )
    ffcolData <- import(fff)
    if (length(con@binSize))
        ffcolData[["bin_size"]] <- con@binSize
    list(
        imgData = DataFrame(
            do.call(rbind, DFs)
        ),
        colData = ffcolData
    )
})

.getImgRow <- function(con, sampleId, image, scaleFx) {
    scfactor <- NA_integer_
    fileNames <- names(con)
    filePaths <- path(con)
    imgFile <- grep(image, fileNames, value = TRUE)
    imgPath <- filePaths[endsWith(filePaths, imgFile)]
    if (!length(imgPath))
        stop(
            "The '", image, "' image was not found in the list of file names.",
            call. = FALSE
        )
    spi <- SpatialExperiment::SpatialImage(imgPath)
    scaleName <- grep(image, names(scaleFx), value = TRUE)
    if (length(scaleName))
        scfactor <- unlist(scaleFx[scaleName])

    S4Vectors::DataFrame(
        sample_id = sampleId,
        image_id = image,
        data = I(list(spi)),
        scaleFactor = scfactor
    )
}
