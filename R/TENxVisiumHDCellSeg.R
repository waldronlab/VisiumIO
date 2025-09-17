#' Flip the Y-axis of cell or nucleus segmentations to align with H&E image. 
#'
#' @param sf a sf object read from a `.geojson` file.
#' @param type "POINT" for cell centroid, or "POLYGON" for cell segmentation 
#' mask. Default is "POINT".
#' @param img_height total length along the Y axis of the image. Obtained by 
#' reading in `hires` or `lowre`s `.png` under `/spatial` folder with 
#' `magick::image_read()`.
#' @param scalef scaling factor from `/spatial/scalefactors_json.json` file
#'
#' @returns a sf object with Y-axis of the points or polygons flipped.
#'
#' @examples
#' geo_data_flipped <- flip_sf_Y(sf = geo_data, type = "POLYGON", 
#'                               img_height = 3886, scalef = 0.079)
flip_sf_Y <- function(sf, type = "POINT", img_height, scalef){
  st_geometry(sf) <- st_sfc( 
    lapply(st_geometry(sf), function(geom) {
      coords <- st_coordinates(geom)
      coords[, 2] <- img_height/scalef - coords[, 2] 
      if(type == "POINT"){
        st_point(coords)
      }else{ # type == "POLYGON"
        st_polygon(list(matrix(coords[, 1:2], ncol = 2)))
      }
    }),
    crs = st_crs(sf)
  )
  
  return(sf)
}



#' Reader for Visium HD cell segmented output.
#'
#' @param spacerangerOut path to unzipped Visium HD data download with `/segmented_outputs`
#' @param res "hires" or "lowres" for what .png to read in. 
#'
#' @returns a `SpatialFeatureExperiment` object.
#'
#' @examples
#' vhdsfe <- readVisiumHDCellSeg(spacerangerOut = "~/Desktop", res = "hires")
TENxVisiumHDCellSeg <- function(
    spacerangerOut,
    res
){
  if(res == "hires"){
    png_name = "tissue_hires_image.png"
    scalef_name = "tissue_hires_scalef"
  }else{ # res = "lowres"
    png_name = "tissue_lowres_image.png"
    scalef_name = "tissue_lowres_scalef"
  }
  
  ## Coord ---------------------------------------------------------------
  # Read the GeoJSON file of cell segmentation
  geo_data <- st_read(file.path(spacerangerOut, "segmented_outputs", 
                                "cell_segmentations.geojson"))
  scalef <- rjson::fromJSON(file = file.path(spacerangerOut, "segmented_outputs", "spatial", 
                                             "scalefactors_json.json"))
  image <- image_read(file.path(spacerangerOut, "segmented_outputs", 
                                "spatial", png_name))
  info <- image_info(image)
  
  # Get centroid
  st_crs(geo_data) <- NA
  rownames(geo_data) <- geo_data$cell_id
  centroids <- st_centroid(geo_data)
  centroids$cell_id <- as.character(centroids$cell_id)
  
  
  ## Countmat ------------------------------------------------------------
  countmat_file <- file.path(spacerangerOut, "segmented_outputs", 
                             "filtered_feature_cell_matrix.h5")
  vhdcellsce <- DropletUtils::read10xCounts(countmat_file, col.names = TRUE)
  
  colnames(vhdcellsce) <- sub("cellid_0*([0-9]+)-.*", "\\1", vhdcellsce$Barcode)
  
  rownames(vhdcellsce) <- rowData(vhdcellsce)$Symbol
  rownames(vhdcellsce) <- make.unique(rownames(vhdcellsce))
  
  
  ## Matching coordinates with count matrix ------------------------------
  # Subset to common cells
  common_cells <- intersect(centroids$cell_id, colnames(vhdcellsce))
  
  geo_data <- geo_data[geo_data$cell_id %in% common_cells, ]
  centroids <- centroids[centroids$cell_id %in% common_cells, ]
  vhdcellsce <- vhdcellsce[, colnames(vhdcellsce) %in% common_cells]
  
  # # Sanity check on ordering of cell ids
  # all(as.character(centroids$cell_id) == colnames(vhdcellsce))
  
  # Extract coordinates
  coords <- st_coordinates(centroids)
  colnames(coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  
  
  ## Construct SPE -------------------------------------------------------
  vhd <- SpatialExperiment(
    assays = list(counts = as(counts(vhdcellsce), "dgCMatrix")),
    rowData = rowData(vhdcellsce),
    colData = cbind(colData(vhdcellsce), coords),
    spatialCoordsNames = colnames(coords),
    scaleFactors = scalef$tissue_hires_scalef, 
    imageSources = file.path(spacerangerOut, "segmented_outputs", "spatial", png_name), 
    loadImage = TRUE
  )
  imgData(vhd)$image_id <- res
  
  ## Add polygons to metadata()
  metadata(vhd)$cellSeg <- geo_data
  
  
  # Flip Y in for metadata() ---------------------------------------------
  # Flip cellSeg Y
  metadata(vhd)$updatecellSeg <- 
    flip_sf_Y(geo_data, type = "POLYGON", img_height = info$height, 
              scalef = imgData(vhd)$scaleFactor) 
  
  # Add high/low res coords to `colData()` -----------------------------------
  vhd[[paste0("pxl_col_in_", res)]] <- 
    spatialCoords(vhd)[, "pxl_col_in_fullres"] * imgData(vhd)$scaleFactor
  vhd[[paste0("pxl_row_in_", res)]] <- 
    spatialCoords(vhd)[, "pxl_row_in_fullres"] * imgData(vhd)$scaleFactor
  
  return(vhd)
  
}