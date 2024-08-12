validate_soil_layers <- function(soil_layers) {
  # scoping for soil_layers
  valid_soil <- c("bdod", "cec", "cfvo", "clay", "nitrogen", 
                  "ocd", "phh2o", "sand", "silt", "soc")
  
  if (is.null(soil_layers)) {
    stop("`soil_layers` parameter cannot be NULL")
  }
  
  if (length(soil_layers) == 0) {
    stop("`soil_layers` cannot be empty.")
  }
  
  if (!all(sapply(soil_layers, is.character))) {
    stop("All elements in `soil_layers` must be character string.")
  }
  
  invalid_soils <- setdiff(soil_layers, valid_soil)
  
  if (length(invalid_soils) > 0) {
    if(length(invalid_soils) == length(soil_layers)) {
      stop("None of the values in soil_layers are valid. Valid values are: ",
           paste(valid_soil, collapse = ", "), ".")
    } else {
      message("The following values in soil_layers are invalid and will be ignored: ", 
              paste(invalid_soils, collapse = ", "))
    }
  }
  
  # Return valid soil layers
  valid_soil_layers <- unique(soil_layers[soil_layers %in% valid_soil])
  
  return(valid_soil_layers)
}

get_soil_data <- function(valid_soil_layers, bounding_box) {
  # Error handling for `valid_soil_layers`
  if (is.null(valid_soil_layers)) {
    stop("`valid_soil_layers` parameter cannot be NULL")
  }
  
  if (length(valid_soil_layers) == 0) {
    stop("`valid_soil_layers` cannot be empty.")
  }
  
  if (!all(sapply(valid_soil_layers, is.character))) {
    stop("All elements in `valid_soil_layers` must be character string.")
  }
  
  # To download soil raster files for an AOI, the AOI layer needs to be reprojected. 
  # As the ISRIC layers use the Homolosine projection, we need to reproject the SF object:
  igh <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  
  sg_url <- "/vsicurl/https://files.isric.org/soilgrids/latest/data/"
  
  # Set soil depths param: 
  depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")
  
  # Set prediction quanitles:
  prediction_quantiles <- c("Q0.05", "Q0.5", "mean", "Q0.95", "uncertainty")
  
  # Create all combinations: 
  combinations <- expand.grid(depths, prediction_quantiles)
  
  # Paste the combinations together
  soil_comb <- paste(combinations$Var1, combinations$Var2, sep="_")
  
  # Initialize file_list variable
  file_list <- NULL
  
  for (i in 1:length(valid_soil_layers)) {
    for (s in 1:length(soil_comb)) {
      soil_datos <- paste0(valid_soil_layers[i], "/", valid_soil_layers[i], "_", soil_comb[s], ".vrt")
      soil_file <- tempfile(pattern=basename(soil_datos), tmpdir=tempdir(), fileext=".tif")
      
      # Check if temporary file can be created
      if (file.exists(soil_file)) {
        file.remove(soil_file)
        message("Overwriting existing file: ", soil_file)
      }
      
      # Try to excute the gdal-utils function and handle errors
      tryCatch({
        sf::gdal_utils(
          util="translate",
          source=paste0(sg_url, soil_datos),
          destination=soil_file,
          options=c("-projwin", bounding_box, "-projwin_srs", igh, "-q"),
          quiet=TRUE
        )
        
        file_list <- c(file_list, soil_file)
      }, error=function(e) {
        # Notify the error:
        message(sprintf("Failed to process: %s. Error: %s", soil_datos, e$message))
      })
    }
  }
  
  return(file_list)
}

get_bbox_in_crs <- function(sf_aoi, radius_src=NULL, crs_val) {
  
  # Add a 50m buffer to account for differences with reprojection. 
  sf_aoi <- sf::st_buffer(sf_aoi, dist=50)
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(sf_aoi, max_area = area_limit)
  
  # Transform aoi to target CRS
  sf_aoi_crs <- sf::st_transform(sf_aoi, crs_val)
  
  # Return bounding box of area of interest: 
  bbox_crs <- sf::st_bbox(sf_aoi_crs)
  
  # Use the bbox data to define the boundary box limits as used by GDAL. 
  # ul means upper left
  # lr means lower right
  ulx <- bbox_crs$xmin
  uly <- bbox_crs$ymax
  lrx = bbox_crs$xmax
  lry <- bbox_crs$ymin
  bb_crs <- c(ulx, uly, lrx, lry)
  
  return(bb_crs)
}

rast_reproject <- function(raster_layer, sf_layer) {
  # Get sf crs
  sf_crs_wkt <- sf::st_crs(sf_layer)$wkt
  
  # Reproject the raster
  rast_project <- terra::project(raster_layer, sf_crs_wkt)
  
  return(rast_project)
}