#' @include TENxVisiumList-class.R
setClassUnion("TENxGeoJSON_OR_NULL", c("TENxGeoJSON", "NULL"))

#' @include TENxParquet-class.R
setClassUnion("TENxMappingParquet_OR_NULL", c("TENxParquet", "NULL"))

setClassUnion("character_OR_NULL", c("character", "NULL"))

#' @docType class
#'
#' @title A class to represent and import multiple Visium HD samples
#'
#' @description This class contains a `SimpleList` of [TENxVisiumHD] objects
#'   each corresponding to one sample. The provided `spacerangerOut` folder
#'   should contain a `binned_outputs` folder where multiple `bin_size`
#'   subfolders are present, e.g., `square_002um`.
#'
#' @details Typically, the user will provide a path to a directory containing
#'  the output of the `spaceranger count` command. The `spaceranger count`
#'  command outputs a folder containing the "raw" or "filtered"
#'  `()_feature_bc_matrix`.
#'
#' @inheritParams TENxVisiumList-class
#'
#' @slot mapping A [TENxMappingParquet] object or `NULL` containing the barcode
#'   mapping data.
#'
#' @return A [SpatialExperiment][SpatialExperiment::SpatialExperiment-class]
#'   object
#'
#' @exportClass TENxVisiumHD
.TENxVisiumHD <- setClass(
    Class = "TENxVisiumHD",
    contains = "TENxVisium",
    slots = c(
        binSize = "character_OR_NULL",
        cellseg = "logical",
        geojson = "TENxGeoJSON_OR_NULL",
        boundary = "character",
        mapping = "TENxMappingParquet_OR_NULL"
    )
)

.getSpatialPath <- function(path, bin_size = NULL, type = c("bc", "cell")) {
    type <- match.arg(type)
    squaref <- ""

    if (identical(type, "bc")) {
        out_path <- file.path(path, "outs")
        if (dir.exists(out_path))
            squaref <- "outs"
        if (!is.null(bin_size)) {
            path <- file.path(path, "binned_outputs")
            squaref <- paste0("square_", bin_size, "um")
        }
    }

    spatial_out <- file.path(path, squaref, "spatial")

    if (!dir.exists(spatial_out))
        stop("The 'spatial' directory was not found.")

    spatial_out
}

.filter_sort_mtx_files <- function(namesvec) {
    files <- .FEATURE_BC_MATRIX_FILES
    names(files) <- files
    res <- lapply(files, function(file) {
        namesvec[startsWith(namesvec, file)]
    })
    unlist(res)
}

.exclude_mtx_files <- function(filelist) {
    files <- .FEATURE_BC_MATRIX_FILES
    names(files) <- files
    res <- lapply(files, function(file) {
        startsWith(names(filelist), file)
    })
    filelist[!Reduce(`|`, res)]
}

.exclude_h5_files <- function(filelist) {
    filelist[tools::file_ext(names(filelist)) != "h5"]
}

.check_filter_mtx <- function(filelist) {
    afiles <- .filter_sort_mtx_files(names(filelist))
    if (!identical(names(afiles), .FEATURE_BC_MATRIX_FILES))
        stop(
            "'TENxFileList' does not contain the expected files:\n  ",
            .FEATURE_BC_MATRIX_FILES_PRINT
        )
    filelist[afiles]
}

.find_convert_resources <- function(
    path,
    processing,
    format,
    bin_size = NULL,
    type = c("bc", "cell")
) {
    type <- match.arg(type)

    if (!is(path, "TENxFileList")) {
        squaref <-
            .getSpatialPath(path = path, bin_size = bin_size, type = type) |>
                dirname()

        path <- .find_file_or_dir(
            reldir = squaref,
            processing = processing,
            format = format,
            type = type
        )

        fdirname <- paste0(processing, "_feature_", type, "_matrix")
        fdirpath <- file.path(squaref, fdirname)
        spatialpath <- file.path(squaref, "spatial")

        if (
            (identical(format, "mtx") && !all(dir.exists(fdirpath))) ||
            !all(dir.exists(spatialpath))
        ) {
            input_arg <- if (identical(type, "bc"))
                "'spacerangerOut'"
            else
                "'segmented_outputs'"

            stop(
                "The 'spatial' or '", fdirname, "' directory was not found.",
                "\n  Verify ", input_arg, " and 'processing' inputs.",
                call. = FALSE
            )
        }
    } else {
        if (identical(format, "h5"))
            path <- .filter_h5_files(path, processing, format, type)
    }
    if (identical(format, "mtx"))
        path <- .check_filter_mtx(path)
    path
}

.find_convert_spatial <- function(path, bin_size, type, ...) {
    if (!is(path, "TENxFileList")) {
        path <- .getSpatialPath(path = path, bin_size = bin_size, type = type)
    } else {
        path <- .exclude_mtx_files(path)
        path <- .exclude_h5_files(path)
    }
    TENxSpatialList(path, ...)
}

.OUTS_FOLDERS <- c("outs", "binned_outputs", "segmented_outputs")

.find_convert_maps <- function(path, pattern) {
    if (is(path, "TENxFileList") || is(path, "TENxFile"))
        path <- dirname(path(path)) |> unique()

    while (any(.OUTS_FOLDERS %in% strsplit(path, .Platform$file.sep)[[1L]])) {
        path <- dirname(path)
        .find_convert_maps(path, pattern)
    }

    mapfile <- list.files(
        path = path, pattern = pattern, full.names = TRUE, recursive = TRUE
    )
    if (!isScalarCharacter(mapfile) || !file.exists(mapfile))
        NULL
    else
        TENxParquet(mapfile)
}

#' @rdname TENxVisiumHD-class
#'
#' @inheritParams TENxVisium
#' @inheritParams TENxVisiumList
#'
#' @param segmented_outputs `character(1)` The path to the `segmented_outputs`
#'   directory
#'
#' @param boundary `character(1)` The type of segmentation boundary to use. The
#'   options are `"cell_segmentations"` (default) or `"nucleus_segmentations"`.
#'
#' @param bin_size `character(1)` The bin size of the images to import. The
#'   default is `008`. It corresponds to the directory name `square_000um` where
#'   `000` is the bin value.
#'
#' @examples
#'
#' vdir <- system.file(
#'     "extdata", package = "VisiumIO", mustWork = TRUE
#' )
#'
#' ## with spacerangerOut folder
#' TENxVisiumHD(spacerangerOut = vdir, bin_size = "002", images = "lowres")
#'
#' TENxVisiumHD(spacerangerOut = vdir, bin_size = "002", images = "lowres") |>
#'     import()
#'
#' ## indicate h5 format
#' TENxVisiumHD(
#'     spacerangerOut = vdir, bin_size = "002",
#'     images = "lowres", format = "h5"
#' )
#'
#' TENxVisiumHD(
#'     spacerangerOut = vdir, bin_size = "002",
#'     images = "lowres", format = "h5"
#' ) |>
#'     import()
#'
#' ## use resources and spatialResource arguments as file paths
#' TENxVisiumHD(
#'     resources = file.path(
#'         vdir, "binned_outputs", "square_002um",
#'         "filtered_feature_bc_matrix.h5"
#'     ),
#'     spatialResource = file.path(
#'         vdir, "binned_outputs", "square_002um",
#'         "spatial"
#'     ),
#'     bin_size = "002", processing = "filtered",
#'     images = "lowres", format = "h5"
#' ) |>
#'     import()
#'
#' ## provide the spatialResource argument as a TENxFileList
#' TENxVisiumHD(
#'     resources = file.path(
#'         vdir, "binned_outputs", "square_002um",
#'         "filtered_feature_bc_matrix.h5"
#'     ),
#'     spatialResource = TENxFileList(
#'         file.path(
#'             vdir, "binned_outputs", "square_002um",
#'             "spatial"
#'         )
#'     ),
#'     bin_size = "002", images = "lowres", format = "h5"
#' ) |>
#'     import()
#'
#' ## with segmented_outputs folder
#' seg_outs <- system.file(
#'     "extdata", "segmented_outputs", package = "VisiumIO", mustWork = TRUE
#' )
#' TENxVisiumHD(
#'     segmented_outputs = seg_outs,
#'     format = "h5",
#'     images = "lowres"
#' ) |>
#'     import()
#' @export
TENxVisiumHD <- function(
    resources,
    spatialResource,
    spacerangerOut,
    segmented_outputs,
    sample_id = "sample01",
    processing = c("filtered", "raw"),
    format = c("mtx", "h5"),
    images = c("lowres", "hires", "detected", "aligned_fiducials"),
    bin_size = c("008", "016", "002"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions\\.parquet",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    mappingPattern = "barcode_mappings\\.parquet",
    boundary = c("cell_segmentations", "nucleus_segmentations"),
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    processing <- match.arg(processing)
    bin_size <-
        if (missing(bin_size)) NULL else match.arg(bin_size)
    format <- match.arg(format)
    boundary <- match.arg(boundary)
    cellseg <- FALSE
    geojson <- NULL

    if (!missing(resources) && !missing(spatialResource)) {
        stopifnot(
            (isScalarCharacter(resources) && file.exists(resources)) ||
                is(resources, "TENxFileList_OR_TENxH5"),
            (isScalarCharacter(spatialResource) &&
                file.exists(spatialResource)) ||
                    is(spatialResource, "TENxFileList")
        )
        if (
            !is(resources, "TENxFileList_OR_TENxH5") &&
            identical(tools::file_ext(resources), "h5")
        )
            resources <- TENxH5(resources, ranges = NA_character_)
        else if (is.character(resources))
            resources <- TENxFileList(resources, ...)

        if (!is(spatialResource, "TENxSpatialList"))
            spatialResource <- TENxSpatialList(
                resources = spatialResource,
                sample_id = sample_id,
                images = images,
                jsonFile = jsonFile,
                tissuePattern = tissuePattern
            )
    } else {
        if (!missing(segmented_outputs)) {
            stopifnot(
                isScalarCharacter(segmented_outputs),
                dir.exists(segmented_outputs)
            )
            tissuePattern <- bin_size <- NULL
            geojson <- TENxGeoJSON(
                file.path(
                    segmented_outputs, paste0(boundary, ".geojson")
                )
            )
            data_folder <- segmented_outputs
            cellseg <- TRUE
            type <- "cell"
        } else if (!missing(spacerangerOut)) {
            stopifnot(
                isScalarCharacter(spacerangerOut),
                dir.exists(spacerangerOut)
            )
            data_folder <- spacerangerOut
            type <- "bc"
        }
        resources <- .find_convert_resources(
            path = data_folder,
            processing = processing,
            format = format,
            bin_size = bin_size,
            type = type
        )
        spatialResource <- .find_convert_spatial(
            path = data_folder,
            bin_size = bin_size,
            type = type,
            sample_id = sample_id,
            images = images,
            jsonFile = jsonFile,
            tissuePattern = tissuePattern
        )
    }
    mapping <- .find_convert_maps(
        path = resources, pattern = mappingPattern
    )

    txv <- TENxVisium(
        resources = resources,
        spatialResource = spatialResource,
        sample_id = sample_id,
        processing = processing,
        format = format,
        images = images,
        jsonFile = jsonFile,
        tissuePattern = tissuePattern,
        spatialCoordsNames = spatialCoordsNames,
        mappingPattern = mappingPattern,
        ...
    )

    .TENxVisiumHD(
        txv,
        cellseg = cellseg,
        boundary = boundary,
        geojson = geojson,
        binSize = bin_size,
        mapping = mapping
    )
}

# import TENxVisiumHD method ----------------------------------------------

#' @describeIn TENxVisiumHD-class Import Visium HD data from multiple bin sizes
#'
#' @inheritParams TENxVisiumList
#'
#' @author E. Y. Dong, M. Ramos
#'
#' @exportMethod import
setMethod("import", "TENxVisiumHD", function(con, format, text, ...) {
    .import_fun <-
        if (con@cellseg)
            .import_cellseg
        else if (!is.null(con@binSize))
            .import_binsize

    .import_fun(con)
})

.import_cellseg <- function(con) {
    checkInstalled("sf")

    sce <- import(con@resources)
    slist <- import(con@spatialList)
    img <- slist[["imgData"]]

    geo_data <- import(con@geojson)
    geo_data <- .repair_cell_ids(geo_data, colnames(sce))
    centroids <- sf::st_centroid(geo_data)

    sce <- .add_map_to_sce(sce, con)
    keep_idx <- colnames(sce) %in% centroids[["cell_id"]]
    sce <- sce[, keep_idx]

    coords <- sf::st_coordinates(centroids)
    colnames(coords) <- con@coordNames
    rownames(coords) <- centroids[["cell_id"]]

    SpatialExperiment(
        assays = list(counts = assay(sce)),
        rowData = rowData(sce),
        mainExpName = mainExpName(sce),
        altExps = altExps(sce),
        sample_id = con@sampleId,
        colData = colData(sce),
        spatialCoords = coords,
        imgData = img,
        metadata = list(
            resouces = metadata(sce),
            spatialList = metadata(con@spatialList),
            cellseg = geo_data,
            boundary = con@boundary
        )
    )
}

.add_map_to_sce <- function(sce, con) {
    bin_size <- con@binSize
    hasMap <- !is.null(con@mapping)
    if (!hasMap)
        return(sce)

    map <- import(con@mapping)
    if (!is.null(bin_size)) {
        binCol <- grepv(bin_size, names(map), TRUE)
        matches <- match(colnames(sce), map[[binCol]])
        if (length(binCol) && length(matches) && !all(is.na(matches))) {
            matched_map <-
                map[matches, , drop = FALSE]
            colData(sce) <- cbind(colData(sce), matched_map)
        }
    } else {
        map <- map[map[["cell_id"]] %in% colnames(sce), ]
        map <- map[map[["in_cell"]], ]
        map <- split(map, map[["cell_id"]])
        idx <- match(colnames(sce), names(map))
        sce$map <- map[idx]
    }
    sce
}

.import_binsize <- function(con) {
    sce <- import(con@resources)
    slist <- import(con@spatialList)
    img <- slist[["imgData"]]
    sce <- .add_map_to_sce(sce, con)

    SpatialExperiment(
        assays = list(counts = assay(sce)),
        rowData = rowData(sce),
        mainExpName = mainExpName(sce),
        altExps = altExps(sce),
        sample_id = con@sampleId,
        colData = colData(sce),
        imgData = img,
        metadata = list(
            resouces = metadata(sce),
            spatialList = metadata(con@spatialList)
        )
    )
}
