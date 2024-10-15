.compare_list_bcodes <- function(bcode_list) {
    inRes <- do.call(
        `%in%`, unname(lapply(bcode_list, unlist, use.names = FALSE))
    )
    comp1 <- names(bcode_list)[1L]
    res_frame <- cbind.data.frame(bcode_list[[comp1]], inRes)
    colname1 <- paste0(comp1, "_", "barcodes")
    colname2 <- paste0("in_", names(bcode_list)[2L])
    names(res_frame) <- c(colname1, colname2)
    res_frame
}

.get_processing_type <- function(filepath) {
    filename <- basename(filepath)
    splitname <- strsplit(filename, "_")[[1L]]
    haskey <- splitname %in% c("raw", "filtered")
    if (!identical(sum(haskey), 1L))
        stop("'processing' type could not be determined from the file name.")
    splitname[haskey]
}

#' Compare barcodes between raw and filtered data
#'
#' @description This function compares the barcodes between raw and filtered
#'  data **depending** on the order of `processing`. Typically, the "raw"
#'  barcodes are compared to the "filtered" ones. The presence of raw
#'  barcodes in the filtered data are marked as `TRUE` in the resulting
#'  `data.frame`.
#'
#' @param from_resource `character(1)` The path to the resource file whose
#'   barcodes are used as the basis of the comparison; typically, the "raw"
#'   feature barcodes are used.
#'
#' @param to_resource `character(1)` The path to the resource file whose
#'   barcodes are compared to the `from_resource`; typically, the "filtered"
#'   feature barcodes.
#'
#' @param spacerangerOut `character(1)` A scalar vector specifying the path
#'   to the `space ranger` output directory.
#'
#' @param processing `character(2)` A vector of length 2 that corresponds to the
#'   processing type. The processing types are typically "raw" and "filtered".
#'   These are the prefixes of the folder names `raw_feature_bc_matrix` and
#'   `filtered_feature_bc_matrix`. The order of the vector determines the
#'   comparison. By default, `processing = c("raw", "filtered")`, which means
#'   barcodes in the raw data are compared to the filtered data.
#'
#' @return A `data.frame` with barcodes of the first element in the `processing`
#'   data type as the first column and a logical vector indicating whether the
#'   barcodes are found in the second element in `processing`. For example,
#'   if processing is `c("raw", "filtered")`, then the first column will be the
#'   barcodes in the `raw` data and the second column will be a logical vector
#'   indicating whether the barcodes are found in the `filtered` data.
#'
#' @examples
#' if (interactive()) {
#'     compareBarcodes(
#'         from_resource =
#'             "./V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz",
#'         to_resource =
#'             "./V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz",
#'     ) |>
#'     head()
#' }
#' @export
compareBarcodes <- function(
    from_resource, to_resource, spacerangerOut,
    processing = c("raw", "filtered"),
    ...
) {
    stopifnot(
        all(c("raw", "filtered") %in% processing),
        identical(length(processing), 2L)
    )
    res <- structure(vector("list", length = 2L), .Names = processing)
    if (!missing(spacerangerOut)) {
        if (isScalarCharacter(spacerangerOut))
            stopifnot(dir.exists(spacerangerOut))
        for (process in processing)
            res[[process]] <-
                .find_convert_resources(spacerangerOut, process, ...)
    } else {
        if (missing(from_resource) || missing(to_resource))
            stop("Both *_resource arguments must be provided")
        stopifnot(
            isScalarCharacter(from_resource), isScalarCharacter(to_resource),
            file.exists(from_resource), file.exists(to_resource)
        )
        resources <- c(from_resource, to_resource)
        rnames <- vapply(resources, .get_processing_type, character(1))
        names(resources) <- rnames
        message(
            "Comparing 'processing' types from ", rnames[1], " to ", rnames[2]
        )
        for (process in processing) {
            txfl <- TENxFileList(resources[process])
            if (txfl@compressed)
                res[[process]] <- decompress(con = txfl)
            else
                res[[process]] <- txfl
        }
    }

    bcode_files <- lapply(
        res, function(x) {
            TENxFile(path(x)[startsWith(names(x), "barcode")])
        }
    )
    bcode_frames <- lapply(bcode_files, import)
    .compare_list_bcodes(bcode_frames)
}
