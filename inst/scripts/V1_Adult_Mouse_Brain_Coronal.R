## URLs of source files
# https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_1
# https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_2

library(rhdf5)

# section 1 ---------------------------------------------------------------
h5f1 <-
    "~/data/V1_Adult_Mouse_Brain_Coronal_Section_1_raw_feature_bc_matrix.h5"
allrows1 <- h5read(h5f1, "/matrix/features/id")

# section 2 ---------------------------------------------------------------
h5f2 <-
    "~/data/V1_Adult_Mouse_Brain_Coronal_Section_2_raw_feature_bc_matrix.h5"
allrows2 <- h5read(h5f2, "/matrix/features/id")

## common in both
commonrows <- intersect(allrows1, allrows2)
## use first 50
commonrows <- commonrows[1:50]

h5new1 <-
    "~/data/V1_Adult_Mouse_Brain_Coronal_Section_1_raw_feature_bc_matrix_new.h5"
# file.remove(h5new1)

newmat1 <- HDF5Array::TENxMatrix(h5f1, "matrix")

library(VisiumIO)
## example object from SpatialExperiment
sample_dir <- system.file(
    file.path("extdata", "10xVisium", "section1"),
    package = "SpatialExperiment"
)

## using spacerangerOut folder
tv <- TENxVisium(
    spacerangerOut = sample_dir, processing = "raw", images = "lowres"
)

s1 <- import(tv)

scols <- na.omit(match(colnames(s1), colnames(newmat1)))
newmat1 <- newmat1[, scols]

identical(colnames(newmat1), colnames(s1))

tkeys <- h5read(h5f1, "/matrix/features/_all_tag_keys")
myidx <- list(na.omit(match(rownames(s1), rownames(newmat1))))

newmat1 <- newmat1[myidx[[1]], ]
## some are missing
identical(rownames(newmat1), rownames(s1))
all(rownames(newmat1) %in% rownames(s1))
setdiff(rownames(s1), rownames(newmat1))
#' [1] "ENSMUSG00000104217"

HDF5Array::writeTENxMatrix(newmat1, h5new1, group = "matrix", verbose = TRUE)

h5ls(h5f1)
h1 <- h5read(h5f1, "/matrix/features/feature_type", index = myidx)
h2 <- h5read(h5f1, "/matrix/features/genome", index = myidx)
h3 <- h5read(h5f1, "/matrix/features/id", index = myidx)
h4 <- h5read(h5f1, "/matrix/features/name", index = myidx)

h5createGroup(h5new1, "/matrix/features/")
h5write(tkeys, h5new1, "/matrix/features/_all_tag_keys")
h5write(h1, h5new1, "/matrix/features/feature_type")
h5write(h2, h5new1, "/matrix/features/genome")
h5write(h3, h5new1, "/matrix/features/id")
h5write(h4, h5new1, "/matrix/features/name")

TENxH5(h5new1, ranges = NA_character_) |> import()

dir.create(sec1 <- "~/gh/TENxIO/inst/extdata/10xVisium/section1/outs", recursive = TRUE)
file.copy(h5new1, file.path(sec1, gsub("_new", "", basename(h5new1))))

# section 2 ---------------------------------------------------------------

h5new2 <-
    "~/data/V1_Adult_Mouse_Brain_Coronal_Section_2_raw_feature_bc_matrix_new.h5"
# file.remove(h5new2)

newmat2 <- HDF5Array::TENxMatrix(h5f2, "matrix")

scols <- match(colnames(s1), colnames(newmat2))
newmat2 <- newmat2[, scols]

## all colnames found
identical(colnames(newmat2), colnames(s1))

tkeys <- h5read(h5f2, "/matrix/features/_all_tag_keys")
myidx <- list(na.omit(match(rownames(s1), rownames(newmat2))))

newmat2 <- newmat2[myidx[[1]], ]

## some are missing
identical(rownames(newmat2), rownames(s1))
all(rownames(newmat2) %in% rownames(s1))
setdiff(rownames(s1), rownames(newmat2))
#' [1] "ENSMUSG00000104217"

HDF5Array::writeTENxMatrix(newmat2, h5new2, group = "matrix", verbose = TRUE)

h5ls(h5f2)
h1 <- h5read(h5f2, "/matrix/features/feature_type", index = myidx)
h2 <- h5read(h5f2, "/matrix/features/genome", index = myidx)
h3 <- h5read(h5f2, "/matrix/features/id", index = myidx)
h4 <- h5read(h5f2, "/matrix/features/name", index = myidx)

h5createGroup(h5new2, "/matrix/features/")
h5write(tkeys, h5new2, "/matrix/features/_all_tag_keys")
h5write(h1, h5new2, "/matrix/features/feature_type")
h5write(h2, h5new2, "/matrix/features/genome")
h5write(h3, h5new2, "/matrix/features/id")
h5write(h4, h5new2, "/matrix/features/name")

TENxH5(h5new2, ranges = NA_character_) |> import()

dir.create(sec2 <- "~/gh/TENxIO/inst/extdata/10xVisium/section2/outs", recursive = TRUE)
file.copy(h5new2, file.path(sec2, gsub("_new", "", basename(h5new2))))

# SpatialExperiment files -------------------------------------------------

se_dir <- system.file(
    file.path("extdata", "10xVisium"),
    package = "SpatialExperiment"
)

file.copy(se_dir, "~/gh/VisiumIO/inst/extdata")
