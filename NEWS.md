## Changes in version 1.6.2

### New features

* Added the `st_invert_y` helper function (from @estellad, #15) to invert
  y-coordinates in `sf` objects.
* Added support for `TENxGeoJSON` format to import cell segmentation data
  from the `segmented_outputs` folder. This allows working with segmentation
  data as `sf` objects within `SpatialExperiment`
* `tissuePos` and `binSize` can be optional parameters for `TENxVisiumHD`
  imports, esp. for `segmented_outputs` inputs 
  
## Bug fixes and minor improvements

* Fixed issue when using `sample_id` parameter in `TENxVisiumHD`
  (@michaelplynch, #17)
* Allow `spacerangerOut` inputs in `TENxVisium` to have an `outs` folder and optionally
  check for the `spatial` subfolder within it (@estellad, #18).
* Folder inputs to arguments such as `spacerangerOut` are not required to
contain an `outs/` directory. These inputs should be used directly, if
available.

## Changes in version 1.4.0

### New features

* Added support for CytAssist images (@ZheFrench, #8)
* Enabled use of alternative readers via `VisiumIO.csvreader` option 
* Included metadata in `SpatialExperiment` outputs

### Bug fixes

* Resolved issue where `*_feature_bc_matrix` folder checks were incorrectly
  triggered for non-mtx formats (@estellad, #4)
* Fixed missing package anchor in `@param` documentation entries

## Changes in version 1.2.0

### New features

* Support VisiumHD file formats including `parquet`.
* Include `format` argument for `h5` file imports (default remains `mtx`).
* Add example data for Visium and VisiumHD imports.

## Changes in version 1.0.0

* Package released in Bioconductor!
