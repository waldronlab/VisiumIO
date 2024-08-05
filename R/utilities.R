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

#' Compare barcodes between raw and filtered data
#'
#' @description This function compares the barcodes between raw and filtered
#'  data **depending** on the order of `processing`. Typically, the "raw"
#'  barcodes are compared to the "filtered" ones. The presence of raw
#'  barcodes in the filtered data are marked as `TRUE` in the resulting
#'  `data.frame`.
#'
#' @param raw_filt_resources `character(2)` A vector of length 2, where the
#'   first element is the path to the **raw** resources and the second element
#'   is the path to the **filtered** resources. If the resources are already
#'   converted to `TENxFileList` objects, then this argument can be a list of
#'   length 2, where the first element is the raw resources and the second
#'   element is the filtered resources.
#'
#' @param spacerangerOut `character(1)` A scalar vector specifying the path
#'   to the `spaceranger` output directory.
#'
#' @param processing `character(2)` A vector of length 2 that corresponds to the
#'   processing type. The processing types are typically "raw" and "filtered".
#'   These are the prefixes of the folder names `raw_feature_bc_matrix` and
#'   `filtered_feature_bc_matrix`. The order of the vector determines the
#'   comparison. For example, if `processing = c("raw", "filtered")`, then the
#'   barcodes in the raw data are compared to the filtered data. Default is
#'   `c("raw", "filtered")`.
#'
#' @param outs `character(1)` A single string specifying the name of the
#'   `spaceranger` output directory. Default is "outs".
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
#'         raw_filt_resources = c(
#'             "~/data/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz",
#'             "~/data/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"
#'         ),
#'         processing = c("raw", "filtered"),
#'         outs = "outs"
#'     ) |>
#'     head()
#' }
#' @export
compareBarcodes <- function(
    raw_filt_resources, spacerangerOut,
    processing = c("raw", "filtered"), outs = "outs", ...
) {
    if (!identical(length(processing), 2L))
        stop("Length of 'processing' must be 2, e.g., c('raw', 'filtered')")
    res <- structure(vector("list", 2L), .Names = processing)
    if (!missing(spacerangerOut)) {
        if (isScalarCharacter(spacerangerOut))
            stopifnot(dir.exists(spacerangerOut))
        for (process in processing) {
            res[[process]] <-
                .find_convert_resources(spacerangerOut, process, ...)
        }
    } else {
        if (
            !identical(length(raw_filt_resources), 2L) &&
            isCharacter(raw_filt_resources)
        )
            stop("'raw_filt_resources' must be a character vector of length 2")
        keywordInFile <- mapply(
            function(file, proc) {
                grepl(pattern = proc, x = file, fixed = TRUE)
            }, raw_filt_resources, processing
        )
        if (!all(keywordInFile))
            stop("No 'processing' match in 'raw_filt_resources' file names")
        for (process in processing) {
            file <- grep(process, raw_filt_resources, value = TRUE)
            txfl <- TENxFileList(file)
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
