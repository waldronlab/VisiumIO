library(rhdf5)
orig <-
    "~/data/Visium_HD_Mouse_Brain_binned_outputs/binned_outputs/square_002um/"
setwd(orig)

R.utils::gunzip(
    "filtered_feature_bc_matrix/matrix.mtx.gz",
    "filtered_feature_bc_matrix/test_matrix.mtx",
    remove = FALSE
)
## reads as dgTMatrix
mm <- readMM("filtered_feature_bc_matrix/test_matrix.mtx")
writeMM(mm[1:10, 1:10], file = "~/gh/VisiumIO/inst/extdata/VisiumHD/binned_outputs/square_002um/filtered_feature_bc_matrix/matrix.mtx")
## file.remove("filtered_feature_bc_matrix/test_matrix.mtx")

setwd("~/gh/VisiumIO/inst/extdata/binned_outputs/square_002um/")
readMM("filtered_feature_bc_matrix/matrix.mtx.gz")

# features.tsv.gz ---------------------------------------------------------

setwd(orig)
feats <- read.delim(
    "filtered_feature_bc_matrix/features.tsv.gz",
    header = FALSE, sep = "\t", nrows = 10
)
pkg_data <- "~/gh/VisiumIO/inst/extdata/binned_outputs"
write.table(
    feats, file =
        file.path(
            pkg_data, "square_002um/filtered_feature_bc_matrix/features.tsv.gz"
        ),
    row.names = FALSE, sep = "\t", col.names = FALSE
)
# read.delim(
#     file.path(
#           pkg_data, "square_002um/filtered_feature_bc_matrix/features.tsv.gz"
#     ),
#     header = FALSE, sep = "\t"
# )

# barcodes.tsv.gz ---------------------------------------------------------

bcodes <- read.delim(
    "filtered_feature_bc_matrix/barcodes.tsv.gz",
    header = FALSE, sep = "\t", nrows = 10
)
write.table(
    bcodes,
    file = file.path(
        pkg_data,
        "square_002um/filtered_feature_bc_matrix/barcodes.tsv.gz"
    ),
    row.names = FALSE, sep = "\t", col.names = FALSE
)
# read.delim(
#     file.path(
#         pkg_data,
#         "square_002um/filtered_feature_bc_matrix/barcodes.tsv.gz"
#     ),
#     header = FALSE, sep = "\t"
# )

# subset MTX format -------------------------------------------------------

dir.create(
    "~/gh/VisiumIO/inst/extdata/VisiumHD/binned_outputs/square_002um",
    recursive = TRUE
)
file.copy("test_write_matrix.mtx", "~/gh/VisiumIO/inst/extdata/VisiumHD/")

setwd("~/gh/VisiumIO/inst/extdata/VisiumHD/binned_outputs/")
h5f <- file.path(orig, "filtered_feature_bc_matrix.h5")
stopifnot(
    file.exists(h5f)
)
h5new <- "square_002um/filtered_feature_bc_matrix.h5"
newmat <- HDF5Array::TENxMatrix(h5f, "matrix")[1:10, 1:10]
HDF5Array::writeTENxMatrix(newmat, h5new, group = "matrix", verbose = TRUE)

tkeys <- h5read(h5f, "/matrix/features/_all_tag_keys")
myidx <- list(1:10)
h1 <- h5read(h5f, "/matrix/features/feature_type", index = myidx)
h2 <- h5read(h5f, "/matrix/features/genome", index = myidx)
h3 <- h5read(h5f, "/matrix/features/id", index = myidx)
h4 <- h5read(h5f, "/matrix/features/name", index = myidx)
h5 <- h5read(h5f, "/matrix/features/target_sets", index = myidx)

h5createGroup(h5new, "/matrix/features/")
h5write(tkeys, h5new, "/matrix/features/_all_tag_keys")
h5write(h1, h5new, "/matrix/features/feature_type")
h5write(h2, h5new, "/matrix/features/genome")
h5write(h3, h5new, "/matrix/features/id")
h5write(h4, h5new, "/matrix/features/name")
h5write(h5, h5new, "/matrix/features/target_sets")

TENxH5(h5new)
file.info(h5new)$size

park <- arrow::read_parquet(
    "~/data/Visium_HD_Mouse_Brain_binned_outputs/binned_outputs/square_002um/spatial/tissue_positions.parquet"
)
ss <- TENxFileList(
    "inst/extdata/VisiumHD/binned_outputs/square_002um/filtered_feature_bc_matrix/"
) |> import()
writepq <-
    park[park$barcode %in% colnames(ss),]

arrow::write_parquet(
    writepq,
    sink = "inst/extdata/VisiumHD/binned_outputs/square_002um/spatial/tissue_positions.parquet"
)
