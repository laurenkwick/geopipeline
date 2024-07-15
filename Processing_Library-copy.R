# library(devtools)
# create_package("C:/Users/laure/Documents/rrr/geopipeline")

###############################
###############################
##### Dependent Libraries #####
###############################
###############################


# library(rsi)
library(sf)
library(terra)
library(fst)
library(units)
library(pak)
pak::pak("Permian-Global-Research/rsi")

####################################
####################################
##### Private Helper Functions #####
####################################
####################################

# Description: Apply a buffer to sf point geometries and return an sf polygon object
# Arguments:
  # boundary_layer: sf object with POINT geometry type.
  # radius: numeric value specifying distance of buffer to apply.
# Value: sf object with polygon geometry (if radius > 0). sf object with point geometry (if radius NULL)

point_buffer <- function(boundary_layer, radius=NULL, max_radius) {
  matching_layer <- boundary_layer

  # Check that radius argument is not null and is numerical
  if (!is.null(radius) && !is.numeric(radius)) {
    stop(deparse(substitute(radius)), " must be numeric or NULL when ", substitute(boundary_layer), " has POINT geometry type. \n")
  }

  if (is.null(radius)) {
    message("No radius defined. Will not apply a buffer to the point data. \n")
    return(matching_layer)
  }

  if (radius > max_radius) {
    warning(deparse(substitute(radius)), " exceeds ", max_radius, " m. Processing time and performance may be impacted. \n",
            immediate. = TRUE)
  }

  # Add buffer to point geometry and return new sf polygon layer
  matching_layer <- sf::st_buffer(boundary_layer, dist=radius)

  return(matching_layer)
}

# Description: Validate user-specified start and end dates for searching over STAC data.
# Arguments:
  # start_date - character representing date in format "YYYY-MM-DD"
  # end_date - character representing date in format "YYYY-MM-DD"
# Value: character vector of length 2 containing start date and end date to be used in get_stac_data search.

valid_dates <- function(start_date, end_date){
  # Try to convert the input strings into dates
  start_date_converted <- try(as.Date(start_date, format = "%Y-%m-%d"), silent=TRUE)
  end_date_converted <- try(as.Date(end_date, format = "%Y-%m-%d"), silent=TRUE)

  # Check if the conversion was successful
  if (inherits(start_date_converted, "try-error") || is.na(start_date_converted)) {
    stop(deparse(substitute(start_date)), " is not a real date or is not in the format YYYY-MM-DD. \n")
  }
  if (inherits(end_date_converted, "try-error") || is.na(end_date_converted)) {
    stop(deparse(substitute(end_date)), " is not a real date or is not in the format YYYY-MM-DD. \n")
  }

  # Check if the start date is later than the end date
  if (start_date_converted > end_date_converted) {
    stop(deparse(substitute(start_date)), " is later than ", substitute(end_date),". \n")
  }

  return(c(start_date, end_date))
}

# Description: Apply Sentinel-2 DN to SR (surface reflectance) conversion.
# Arguments:
  # raster - terra SpatRaster object
# Value: character, File location where scaled GeoTIFF is stored.

s2_scale <- function(raster) {
  # Check if the "SCL" layer exists
  if ("SCL" %in% names(raster)) {
    # Identify the "SCL" layer index
    scl_layer_index <- which(names(raster) == "SCL")
  } else {
    scl_layer_index <- NULL
  }

  # Create a copy of the raster to store the scaled layers
  scaled_raster <- raster

  # Loop through all layers and apply scaling factor, except for the "SCL" layer
  for (i in 1:terra::nlyr(raster)) {
    if (!is.null(scl_layer_index) && i != scl_layer_index) {
      scaled_raster[[i]] <- raster[[i]] / 10000
    } else if (!is.null(scl_layer_index) && i == scl_layer_index) {
      scaled_raster[[i]] <- raster[[i]]
    } else {
      scaled_raster[[i]] <- raster[[i]] / 10000
    }
  }

  output_filename <- tempfile(pattern="SR_image", tmpdir=tempdir(), fileext=".tif")
  w <- terra::writeRaster(scaled_raster, filename = output_filename,
                          datatype="FLT4S", filetype="GTiff",
                          gdal=c("COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS", "PREDICTOR=2"),
                          tempdir=tempdir(), NAflag=NA, verbose=FALSE)

  return(output_filename)
}

# Description: Generate data frame of indices to calculate.
# Arguments:
  # band_names: character of length > 0. Band names as required by the Awesome Spectral Indices (ASI) repository.
  # index_application: character specifying application domains to include in calculations. This value is
    # based on the application_domain column from the ASI repository. At the time of writing
    # this, available application domains for Sentinel-2 data are “vegetation”, “water”,
    # “burn”, “soil”, “urban”, “snow”, “kernel”.
  # index_names: character specifying index names to calculate. This value is based on the short_name column from
    # the ASI repository.
# Value: data frame object containing list of indices that match user-specified search criteria and data set.

indices_df <- function(band_names, index_application=NULL, index_names=NULL) {

  # Return a data frame of spectral indices that match to Sentinel2 sensor and bands
  index_df <- rsi::spectral_indices() |>
    rsi::filter_bands(band_names)

  # Return error message if the dataframe contains no data
  if (nrow(index_df)==0) {
    message("No indices found that match the bands in dataset. Will not calculate indices. \n")
    return(NULL)
  }

  # Initialize a dataframe which will store the filtered indices
  selection_df <- NULL

  # Loop through the user-entered list of application domains and subset index dataframe
  if (!is.null(index_application)) {
    for (a in index_application) {
      if (!a %in% index_df$application_domain) {
        message(a, " is not a valid applciation domain. Proceeding without. \n")
      }
      selection_df <- rbind(selection_df, subset(index_df, application_domain==a))
    }
  }

  # Loop through the user-entered list of index names and subset index dataframe.
  if (!is.null(index_names)) {
    for (n in index_names) {
      if (!n %in% index_df$short_name) {
        message(n, " is not a valid index name. Proceeding without. \n")
      }

      if (!n %in% selection_df$short_name){
        selection_df <- rbind(selection_df, subset(index_df, short_name==n))
      }
    }
  }

  # Stop if no match to user-specified indices.
  if (nrow(selection_df)==0) {
    message("No indices were found that match ", deparse(substitute(index_application)), " and/or ",
            deparse(substitute(index_names))," criteria. Will not calculate indices. \n")
    return(NULL)
  }

  message(nrow(selection_df), " indices were found based on criteria. \n")

  return(selection_df)
}

# Description: Extract spectral and index values for point geometries with a specified buffer. Returns
  # a data.frame object with column corresponding to spectral and index values and rows corresponding to
  # pixels.
# Arguments:
  # raster - terra SpatRaster object to extract values from
  # buffer_layer - sf object with polygon geometry type.
  # point_layer - sf object with point geometry type. buffer_layer was derived from this layer, via the
    # point_buffer function.
# Value:
  # data frame containing spectral and index values.

extract_vals_df <- function(raster, buffer_layer, col_ID=NULL) {

  # Update row.names to use user-defined column ID if populated
  if (!is.null(col_ID)) {
    row.names(buffer_layer) <- buffer_layer[[col_ID]]
  }

  # Convert sf data.frame into a terra SpatVector object.
  layer_spatVect <- terra::vect(buffer_layer)

  # Extract values from the raster for the set of locations.
  aoi_df <- terra::extract(raster, layer_spatVect)

  # Create a mapping from polygon IDs to original row names
  id_to_rowname <- setNames(row.names(buffer_layer), 1:nrow(buffer_layer))

  # Use original row names in output data frame:
  aoi_df$ID <- id_to_rowname[aoi_df$ID]

  # Print warning message of number of rows with miss values.
  # Exclude the geom column when checking for complete cases.
  if (sum(!complete.cases(aoi_df)) > 0) {
    message("There are ", sum(!complete.cases(aoi_df))," out of ", nrow(aoi_df), " rows with NA values. \n")
  }

  return(aoi_df)
}

valid_col_ID <- function(aoi, column_ID=NULL) {

  # Message for no column_ID provided.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi))) && (is.null(column_ID))) {
    message("No unique column ID provided. Generating default unique IDs using `row.names()` \n")
    return(NULL)
  }

  # Message for column_ID provided, but polygon geometry type. s2_process generates a geotiff for polygon
  # geometry, rather than a data.frame, so this argument will not be used.
  if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi))) && (!is.null(column_ID))) {
    message("Data frame not created for polygon geometry. Ignoring ", deparse(substitute(column_ID)), ". \n")
    return(NULL)
  }

  if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi))) && (is.null(column_ID))) {
    return(NULL)
  }

  # Check if the column_ID exists in the sf data.frame
  if (!column_ID %in% names(aoi)) {
    message("Column ", column_ID, " does not exist in the data frame. Generating default unique IDs using `row.names()` \n")
    return(NULL)
  }

  # Check if the values in the column are unique
  column_vals <- aoi[[column_ID]]
  if (length(unique(column_vals)) != length(column_vals)) {
    message("Column ", column_ID, " does not contain unique values. Generating default unique IDs using `row.names()` \n")
    return(NULL)
  }

  # If both criteria above are met, return TRUE. We will use the user-defined column_ID
  return(column_ID)
}

s1_ascending_query <- function(bbox,
                               stac_source,
                               collection,
                               start_date,
                               end_date,
                               limit,
                               ...) {
  #`bbox` is guaranteed to be in 4326 already
  geometry <- rstac::cql2_bbox_as_geojson(bbox)
  #`start_date` and `end_date` will be processed and in RFC-3339 formats
  datetime <- rstac::cql2_interval(start_date, end_date)

  request <- rstac::ext_filter(
    rstac::stac(stac_source),
    collection == {{ "sentinel-1-rtc" }} &&
      t_intersects(datetime, {{ datetime }}) &&
      s_intersects(geometry, {{ geometry }}) &&
      `sat:orbit_state` == "ascending"
  )
  rstac::items_fetch(rstac::post_request(request))
}

s1_descending_query <- function(bbox,
                                stac_source,
                                collection,
                                start_date,
                                end_date,
                                limit,
                                ...) {
  #`bbox` is guaranteed to be in 4326 already
  geometry <- rstac::cql2_bbox_as_geojson(bbox)
  #`start_date` and `end_date` will be processed and in RFC-3339 formats
  datetime <- rstac::cql2_interval(start_date, end_date)

  request <- rstac::ext_filter(
    rstac::stac(stac_source),
    collection == {{ "sentinel-1-rtc" }} &&
      t_intersects(datetime, {{ datetime }}) &&
      s_intersects(geometry, {{ geometry }}) &&
      `sat:orbit_state` == "descending"
  )
  rstac::items_fetch(rstac::post_request(request))
}

s1_dB <- function(raster) {
  # Create a copy of the raster to store the transformed values
  transformed_raster <- raster

  # Loop through all layers and apply the logarithmic transformation
  for (i in 1:terra::nlyr(raster)) {
    transformed_raster[[i]] <- 10*log10(raster[[i]])
    names(transformed_raster[[i]]) <- paste0(names(raster[[i]]), "_dB")
  }

  # Combine the two rasters so that we have backscatter data and dB
  combined <- c(raster, transformed_raster)

  output_filename <- tempfile(pattern="dB_image", tmpdir=tempdir(), fileext=".tif")
  w <- terra::writeRaster(combined, filename = output_filename,
                          datatype="FLT4S", filetype="GTiff",
                          gdal=c("COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS", "PREDICTOR=2"),
                          tempdir=tempdir(), NAflag=NA, verbose=FALSE)

  return(output_filename)
}

check_bbox_area <- function(sf_object, max_area) {

  # Grab the bounding box of the AOI layer
  bbox_aoi <- st_bbox(sf_object)

  # Convert bounding box to an sfc object
  boundary <- st_as_sfc(bbox_aoi)

  # Calculate the area of the object
  area <- st_area(boundary)
  area_sqm <- set_units(area, m^2)

  # Check if the area exceeds the max_area
  if (area_sqm > set_units(max_area, m^2)) {
    warning("The area of the AOI bounding box exceeds ", max_area, " sqm. Processing time and performance may be impacted. \n",
            immediate. = TRUE)
  }

}

check_date_in_range <- function(start_dt, end_dt, target_yr) {
  # Convert the input strings to Date objects
  start_date <- as.Date(start_dt)
  end_date <- as.Date(end_dt)
  
  # Create the start and end dates for the target year
  target_start <- as.Date(paste0(target_yr,"-01-01"))
  target_end <- as.Date(paste0(target_yr, "-12-31"))
  
  # Check if the date ranges overlap and return the target search range
  if (!(end_date < target_start || start_date > target_end)) {
    return(c(target_start, target_end))
  }
  
  # Otherwise, return NULL
  return(NULL)
}

validate_soil_layers <- function(soil_layers) {
  # scoping for soil_layers
  valid_soil <- c("bdod", "cec", "cfvo", "clay", "nitrogen", 
                  "ocd", "ocs", "phh2o", "sand", "silt", "soc")
  
  if (is.null(soil_layers)) {
    stop("`soil_layers` parameter cannot be NULL")
  }
  
  invalid_soils <- NULL
  
  for (i in 1:length(soil_layers)) {
    if (!soil_layers[i] %in% valid_soil) {
      invalid_soils <- c(invalid_soils, soil_layers[i])
    }
  }
  
  if (length(invalid_soils) > 0) {
    if(length(invalid_soils) == length(soil_layers)) {
      stop("None of the values in soil_layers are valid. Valid values are: ",
           paste(valid_soil, collapse = ", "), ".")
    } else {
      warning("The following values in soil_layers are invalid and will be ignored: ", 
              paste(invalid_soils, collapse = ", "))
    }
  }
  
  # Return valid soil layers
  valid_soil_layers <- soil_layers[soil_layers %in% valid_soil]
  
  return(unique(valid_soil_layers))
}

get_soil_data <- function(valid_soil_layers, bounding_box) {
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
      }
      
      # Try to excute the gdal-utils function and handle errors
      tryCatch({
        sf::gdal_utils(
          util="translate",
          source=paste0(sg_url, soil_datos),
          destination=soil_file,
          options=c("-projwin", bounding_box, "-projwin_srs", igh),
          quiet=FALSE
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


#################################
#################################
##### User-Facing Functions #####
#################################
#################################

# Description: User-facing processing function. This function takes in a suite of parameters and performs
  # one of two possible tasks: (1) Serialization of a data.frame object containing spectral data and indices
  # for point geometry data. (2) Generation of a GeoTIFF containing spectral data and indices for a region
  # of interest.
# Arguments:
  # aoi_layer: sf object. Must have POINT or POLYGON geometry. Dictates the output of the function:
    # POINT geometry returns a data frame. POLYGON geometry returns a GeoTIFF.
  # radius: numeric. Required when aoi_layer has POINT geometry. This is the distance used for buffering.
  # start_dt, end_dt: character. The first and last date, respectively, of imagery to include. Must be
    # "YYYY-MM-DD" format. end_dt must occur after start_dt.
  # composite_method: character. The name of the function used to combine downloaded images into a single
    # composite. Must be one of "sum", "median", "mean", "min", or "max". Default is "median"
  # resample_method: character. Resampling method to use. See available resampling methods in gdalwarp help
    # text. Default is "bilinear".
  # app_domains: character specifying application domains to include in calculations. This value is
    # based on the application_domain column from the ASI repository. At the time of writing
    # this, available application domains for Sentinel-2 data are “vegetation”, “water”,
    # “burn”, “soil”, “urban”, “snow”, “kernel”.
  # idx_names: character specifying index names to calculate. This value is based on the short_name column from
    # the ASI repository.
  # file_path: character. The output file path and name to write either the output .tif or .fst file. Do
    # not include the file extension in this variable, it will be added automatically.
# Value: character vector of length 1 or greater. Each element is either an .fst or .tif file.

s2_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL,
                       composite_method=c("median", "mean", "sum", "min", "max", NULL),
                       resample_method=c("bilinear", "average", "rms", "nearest", "gauss", "cube", "cubicspline", "lanczos", "average_magphase", "mode"),
                       keep_SCL=FALSE, apply_mask=TRUE, app_domains=NULL, idx_names=NULL, file_path) {

  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }

  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }

  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)

  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)

  # scoping for composite_method.
  valid_composite <- c("sum", "mean", "median", "min", "max")
  if(!is.null(composite_method[1])){
    if (!composite_method[1] %in% valid_composite) {
      stop(deparse(substitute(composite_method[1])), " must be one of ", paste(valid_composite, collapse=", "), ". \n")
    }
  }

  if (is.null(composite_method[1]) && ((as.Date(dates[2])-as.Date(dates[1]))) > 100) {
    warning("`composite_method` is `NULL` and date range is greater than 100 days. Several files will be created which may impact processing time and performance. \n",
            immediate. = TRUE)
  }

  # scoping for resample_method
  valid_resample <- c("nearest", "average", "rms", "bilinear", "gauss", "cube", "cubicspline",
                      "lanczos", "average_magphase", "mode")
  if (!resample_method[1] %in% valid_resample) {
    stop(deparse(substitute(resample_method[1])), " must be one of ", paste(valid_resample, collapse=", "), ". \n")
  }

  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer

  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }

  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }

  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)

  # Set rsi band mapping object
  s2_asset <- rsi::sentinel2_band_mapping$aws_v1

  # Define band mapping vector.
  s2_vals <- c("B","A","G","N","N2","R","RE1","RE2","RE3","S1","S2")
  names(s2_vals) <- c("blue", "coastal", "green", "nir", "nir08", "red", "rededge1", "rededge2", "rededge3", "swir16", "swir22")

  # Initialize a variable to store indices
  df_indices <- NULL

  # Set value for Sentinel-2 band names and return a dataframe of indices to use IF user provides
  # application domain or specific index names to search over. Otherwise, proceed with only spectral data.
  if (!is.null(app_domains) | !is.null(idx_names)) {
    df_indices <- indices_df(band_names=s2_vals, index_application=app_domains, index_names=idx_names)
  }

  if (keep_SCL == TRUE && !is.null(composite_method)) {
    message("`keep_SCL` ignored when `composite_method` is not `NULL`. \n")
  }

  # Add SCL band if composite method is NULL and keep_SCL is TRUE
  if (keep_SCL == TRUE && is.null(composite_method)) {
    s2_vals <- c(s2_vals, scl="SCL")
  }

  # Initialize pixel masking variables
  mask_layer <- NULL
  mask_fxn <- NULL

  # Define mask band and function if apply_mask = TRUE
  if (apply_mask == TRUE) {
    mask_layer <- attr(s2_asset, "mask_band")
    mask_fxn <- attr(s2_asset, "mask_function")
  }

  # Initialize image variable
  image <- NULL

  tryCatch({
    # Retrieve raster data
    # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
    image <- rsi::get_sentinel2_imagery(
      aoi = sf::st_buffer(aoi_search, dist=50),
      start_date = dates[1],
      end_date = dates[2],
      asset_names = s2_vals,
      stac_source = attr(s2_asset, "stac_source"),
      collection = attr(s2_asset, "collection_name"),
      query_function = attr(s2_asset, "query_function"),
      download_function = attr(s2_asset, "download_function"),
      sign_function = attr(s2_asset, "sign_function"),
      mask_band = mask_layer,
      mask_function = mask_fxn,
      output_filename = tempfile(pattern="S2img", tmpdir=tempdir(), fileext=".tif"),
      composite_function = composite_method[1],
      gdalwarp_options = c("-r", resample_method, "-multi", "-overwrite", "-co",
                           "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
    )
  }, error = function(e) {
    message("No items were found for this combination of AOI and date range for Sentinel-2. \n")
  })

  if (is.null(image)) {
    stop("No Sentinel-2 items were found. \n")
  }

  file_list <- NULL

  for (i in 1:length(image)) {
    # Convert image file into raster
    DN_raster <- terra::rast(image[i])

    # Scale the data
    SR_image <- s2_scale(DN_raster)

    # Add indices
    if (!is.null(df_indices)) {
      if (nrow(df_indices) > 0) {
        indices_image <- rsi::calculate_indices(raster = SR_image, indices = df_indices,
                                                output_filename = tempfile(pattern="index_image", tmpdir=tempdir(), fileext=".tif"))
        SR_image <- rsi::stack_rasters(c(SR_image,indices_image), output_filename=paste0(file_path,"_",i,".tif"))
      }
    }

    # Initialize dataframe for extracting pixel values
    vals_df <- NULL

    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      SR_raster <- terra::rast(SR_image)
      vals_df <- extract_vals_df(SR_raster, aoi_search, use_uniqueID)

      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }

    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", SR_image)
      file_list <- append(file_list,SR_image)
    }

  }

  return(file_list)
}

s1_image <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL,
                     composite_method=c("median", "mean", "sum", "min", "max", NULL),
                     resample_method=c("bilinear", "average", "rms", "nearest", "gauss", "cube", "cubicspline", "lanczos", "average_magphase", "mode"),
                     calc_dB=FALSE, app_domains=NULL, idx_names=NULL, file_path) {

  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), "must be an sf object. \n")
  }

  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of", paste(gtypes, collapse=", "),". \n")
  }

  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)

  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)

  # scoping for composite_method.
  valid_composite <- c("sum", "mean", "median", "min", "max")
  if(!is.null(composite_method[1])){
    if (!composite_method[1] %in% valid_composite) {
      stop(deparse(substitute(composite_method[1])), " must be one of ", paste(valid_composite, collapse=", "),". \n")
    }
  }

  if (is.null(composite_method[1]) && ((as.Date(dates[2])-as.Date(dates[1]))) > 100) {
    warning("`composite_method` is `NULL` and date range is greater than 100 days. Several files will be created which may impact processing time and performance. \n",
            immediate. = TRUE)
  }

  # scoping for resample_method
  valid_resample <- c("nearest", "average", "rms", "bilinear", "gauss", "cube", "cubicspline",
                      "lanczos", "average_magphase", "mode")
  if (!resample_method[1] %in% valid_resample) {
    stop(deparse(substitute(resample_method[1])), " must be one of ", paste(valid_resample, collapse=", "),". \n")
  }

  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer

  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }

  # Print a message if user give radius but sf object has polygon geometry.
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }

  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)

  # Initialize a variable to store indices
  df_indices <- NULL

  # Set the collection to Sentinel 1 RTC
  rtc_collection <- "sentinel-1-rtc"

  # Initialize image variables
  image_asc <- NULL
  image_des <- NULL

  # Try to get imagery from ascending orbit
  tryCatch({
    image_asc <- rsi::get_sentinel1_imagery(
      aoi = sf::st_buffer(aoi_search, dist=50),
      start_date = dates[1],
      end_date = dates[2],
      collection = rtc_collection,
      output_filename = tempfile(pattern="S1img_asc", tmpdir=tempdir(), fileext=".tif"),
      query_function = s1_ascending_query,
      composite_function = composite_method[1],
      gdalwarp_options = c("-r", resample_method[1], "-multi", "-overwrite", "-co",
                           "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
    )
  }, error = function(e) {
    message("No items were found for this combination of AOI and date range for Sentinel-1 ascending orbit. \n")
  })

  # Try to get imagery from descending orbit
  tryCatch({
    image_des <- rsi::get_sentinel1_imagery(
      aoi = sf::st_buffer(aoi_search, dist=50),
      start_date = dates[1],
      end_date = dates[2],
      collection = rtc_collection,
      output_filename = tempfile(pattern="S1img_des", tmpdir=tempdir(), fileext=".tif"),
      query_function = s1_descending_query,
      composite_function = composite_method[1],
      gdalwarp_options = c("-r", resample_method[1], "-multi", "-overwrite", "-co",
                           "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
    )
  }, error = function(e) {
    message("No items were found for this combination of AOI and date range for Sentinel-1 descending orbit. \n")
  })

  # Check if both images are NULL
  if(is.null(image_asc) && is.null(image_des)) {
    stop("No Sentinel-1 items were found for both ascending and descending orbit for this combination of AOI and date range. \n")
  }

  images_combined <- c(image_asc, image_des)
  file_list <- NULL

  for (i in 1:length(images_combined)) {

    # Set value for Sentinel-2 band names and return a dataframe of indices to use IF user provides
    # application domain or specific index names to search over. Otherwise, proceed with only spectral data.
    out_raster <- terra::rast(images_combined[i])

    if (!is.null(app_domains) | !is.null(idx_names)) {
      s1_vals <- names(out_raster)
      df_indices <- indices_df(band_names=s1_vals, index_application=app_domains, index_names=idx_names)
    }

    if (calc_dB==TRUE) {
      images_combined[i] <- s1_dB(out_raster)
    }

    # Add indices
    if (!is.null(df_indices)) {
      if (nrow(df_indices) > 0) {
        indices_image <- rsi::calculate_indices(raster = out_raster, indices = df_indices,
                                                output_filename = tempfile(pattern="index_image", tmpdir=tempdir(), fileext=".tif"))

        images_combined[i] <- rsi::stack_rasters(c(images_combined[i],indices_image), output_filename=paste0(file_path,"_",i,".tif"))
      }
    }

    # Initialize dataframe for extracting pixel values
    vals_df <- NULL

    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      radar_raster <- terra::rast(images_combined[i])
      vals_df <- extract_vals_df(radar_raster, aoi_search, use_uniqueID)

      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }

    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", images_combined[i])
      file_list <- append(file_list,images_combined[i])
    }

  }

  return(file_list)
}

dem_process <- function(aoi_layer, radius=NULL, uniqueID=NULL, file_path) {

  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }

  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }

  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)


  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer

  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }

  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }

  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)

  # Initialize image variable
  image <- NULL
  
  # Set default image file name: 
  image_file <- tempfile(pattern="dem_img", tmpdir=tempdir(), fileext=".tif")
  
  # Override image file name for polygon geometry
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    image_file <- paste0(file_path, ".tif")
  }

  tryCatch({
    # Retrieve raster data
    # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
    image <- rsi::get_dem(
      aoi = sf::st_buffer(aoi_search, dist=50),
      output_filename = image_file
    )
  }, error = function(e) {
    message("No DEM items were found.")
  })

  if (is.null(image)) {
    stop("No DEM items were found.")
  }

  file_list <- NULL

  for (i in 1:length(image)) {
    # Convert image file into raster
    DEM_raster <- terra::rast(image[i])

    vals_df <- NULL

    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      vals_df <- extract_vals_df(DEM_raster, aoi_search, use_uniqueID)

      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }

    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", image[i])
      file_list <- append(file_list,image[i])
    }

  }

  return(file_list)
}


esa_worldcover_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Initialize image variable
  image <- NULL
  
  # Set up sign function to sign the item URLs for PC
  esa_coll <- "esa-worldcover"
  esa_assets <- "map"
  pc_stac_source <- "https://planetarycomputer.microsoft.com/api/stac/v1/"
  
  for(i in seq(2020, 2021, 1)) {
    # Check if query dates contain the ESA Worldcover Temporal extent
    search_dates <- check_date_in_range(dates[1], dates[2], as.character(i))
    if (is.null(search_dates)) {
      next
    }
    
    # Set default image file name: 
    image_file <- tempfile(pattern=paste0("esa_img", as.character(i)), tmpdir=tempdir(), fileext=".tif")
    
    # Override image file name for polygon geometry
    if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      image_file <- paste0(file_path, "_", as.character(i), ".tif")
    }
    
    tryCatch({
      # Retrieve raster data
      # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
      image <- c(image, rsi::get_stac_data(
        aoi = sf::st_buffer(aoi_search, dist=50),
        start_date = as.character(search_dates[1]),
        end_date = as.character(search_dates[2]),
        pixel_x_size = 10,
        pixel_y_size = 10,
        asset_names = esa_assets,
        stac_source = pc_stac_source,
        collection = esa_coll,
        sign_function = rsi::sign_planetary_computer,
        rescale_bands=FALSE,
        download_function = rsi::rsi_download_rasters,
        composite_function = "max",
        gdalwarp_options = c("-r", "near", "-multi", "-overwrite", "-co",
                             "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS"),
        output_filename = image_file
      ))
    }, error = function(e) {
      message("No ESA WorldCover items were found for ", as.character(i))
    })
    
  }
  
  if (is.null(image)) {
    stop("No ESA WorldCover items were found.")
  }
  
  file_list <- NULL
  
  for (i in 1:length(image)) {
    # Convert image file into raster
    ESA_raster <- terra::rast(image[i])
    
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      vals_df <- extract_vals_df(ESA_raster, aoi_search, use_uniqueID)
      
      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", image[i])
      file_list <- append(file_list,image[i])
    }
    
  }
  
  return(file_list)
}

alos_fnf_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Initialize image variable
  image <- NULL
  
  # Set up sign function to sign the item URLs for PC
  alos_coll <- "alos-fnf-mosaic"
  alos_assets <- "C"
  pc_stac_source <- "https://planetarycomputer.microsoft.com/api/stac/v1/"
  
  for (i in seq(2015, 2020, 1)) {
    # Check if query dates contain the ALOS FNF temporal extent
    search_dates <- check_date_in_range(dates[1], dates[2], as.character(i))
    if (is.null(search_dates)) {
      next
    }
    
    # Set default image file name
    image_file <- tempfile(pattern=paste0("alos_img", as.character(i)), tmpdir=tempdir(), fileext=".tif")
    
    # Override image file name for polygon geometry
    if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      image_file <- paste0(file_path, "_", as.character(i), ".tif")
    }
    
    tryCatch({
      # Retrieve raster data
      # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
      image <- c(image, rsi::get_stac_data(
        aoi = sf::st_buffer(aoi_search, dist=50),
        start_date = as.character(search_dates[1]),
        end_date = as.character(search_dates[2]),
        pixel_x_size = 25,
        pixel_y_size = 25,
        asset_names = alos_assets,
        stac_source = pc_stac_source,
        collection = alos_coll,
        sign_function = rsi::sign_planetary_computer,
        rescale_bands=FALSE,
        download_function = rsi::rsi_download_rasters,
        composite_function = "max",
        gdalwarp_options = c("-r", "near", "-multi", "-overwrite", "-co",
                             "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS"),
        output_filename = image_file
      ))
    }, error = function(e) {
      message("No ALOS FNF items were found for ", as.character(i))
    })
  }
  
  if (is.null(image)) {
    stop("No ALOS FNF items were found.")
  }
  
  file_list <- NULL
  
  for (i in 1:length(image)) {
    # Convert image file into raster
    FNF_raster <- terra::rast(image[i])
    
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      vals_df <- extract_vals_df(FNF_raster, aoi_search, use_uniqueID)
      
      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", image[i])
      file_list <- append(file_list,image[i])
    }
    
  }
  
  return(file_list)
}

soil_process <- function(aoi_layer, radius=NULL, uniqueID=NULL, soil_layers) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # Set up aoi search parameter to use in API call
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Get bounding box in Homolosine projection
  # As the ISRIC layers use the Homolosine projection, we need to reproject the SF object:
  igh <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  bbox_igh <- get_bbox_in_crs(sf_aoi=aoi_search, radius_src=radius, crs_val=igh)
  
  # Get valid soil layers to use in GDAL translate utility
  valid_soil_srch <- validate_soil_layers(soil_layers)
  
  # Get list of GeoTIFF soil files: 
  soil_files <- get_soil_data(valid_soil_layers=valid_soil_srch, bounding_box=bbox_igh)
  
  # Convert output GeoTIFFs into terra SpatRaster object
  soil_rast <- terra::rast(soil_files)
  
  # Reproject the raster to match the original sf object projection
  soil_rast <- rast_reproject(raster_layer=soil_rast, sf_layer=aoi_layer)
  
  # Extract values into dataframe for point geometry. 
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    #Intialize dataframe variable: 
    vals_df <- NULL
    
    # Convert to raster and extract pixel values
    vals_df <- extract_vals_df(soil_rast, aoi_search, use_uniqueID)
    
    # Write to fst
    file_nm <- "testing_soil_df.fst"
    fst::write_fst(vals_df, path=file_nm)
    message("File saved to: ", file_nm)
  }
  
  # Save as .tif for polygon geometry. 
  if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    file_nm <- "testing_soil_img.tif"
    writeRaster(soil_rast, filename=file_nm)
    message("File saved to: ", file_nm)
  }
  
  return(file_nm)
  
}

###################
###################
##### Testing #####
###################
###################

# Set working directory
setwd("~/ProcessingModule")

# Create sf objects for testing:
sf_point <- st_read("BraulioCarrillo.gpkg",layer="bioacoustic_plots")
sf_poly <- st_read("BraulioCarrillo.gpkg", layer="study_area")

sf_poly_lg <- st_read("BC_testc.gpkg", layer="large_aoi")

##############
# SENTINEL-2 #
##############

# Run module for point geometry:
s2_t1 <- s2_process(sf_point, radius=100, start_dt="2023-01-01", end_dt="2023-01-15", uniqueID="Id",keep_SCL=TRUE,
                    apply_mask = FALSE, idx_names=c("NDVI","ndkas","NDMI"), file_path = tempfile(pattern="s2_df_", tmpdir=tempdir()))
s2_t1_scl <- s2_process(sf_point, radius=50, start_dt="2023-01-01", end_dt="2023-01-15", uniqueID="Id",keep_SCL=TRUE,
                        apply_mask = FALSE, composite_method=NULL, idx_names="NDVI", file_path = tempfile(pattern="s2_df_", tmpdir=tempdir()))

# Run module for polygon geometry:
s2_t2 <- s2_process(sf_poly, start_dt="2023-01-01", end_dt="2023-01-15", idx_names="NDVI", file_path = tempfile(pattern="s2_tif_", tmpdir=tempdir()))

# Read in the fst file
s2_df_1 <- fst::read.fst(s2_t1[1])

s2_df_scl_1 <- fst::read.fst(s2_t1_scl[1])
s2_df_scl_2 <- fst::read.fst(s2_t1_scl[2])

# Generate SpatRasters
s2_rast_1 <- terra::rast(s2_t2[1])

##############
# SENTINEL-1 #
##############

# Run module for point geometry:
s1_t1 <- s1_image(sf_point, radius=100, start_dt="2023-01-01", end_dt="2023-01-20", uniqueID="Id", composite_method="median",
                  calc_dB=TRUE, file_path = tempfile(pattern="s1_df_", tmpdir=tempdir()))

# Run module for polygon geometry:
s1_t2 <- s1_image(sf_poly, start_dt="2023-01-01", end_dt="2023-01-15", calc_dB=TRUE,
                  idx_names="NDPolI", file_path = tempfile(pattern="s1_tif_", tmpdir=tempdir()))

s1_t2_lg <- s1_image(sf_poly_lg, start_dt="2023-01-01", end_dt="2023-01-15", calc_dB=TRUE, composite_method=NULL,
                     idx_names="NDPolI", file_path = tempfile(pattern="s1_tif_", tmpdir=tempdir()))

s1_lg1 <- terra::rast(s1_t2_lg[1])
s1_lg2 <- terra::rast(s1_t2_lg[2])
s1_lg3 <- terra::rast(s1_t2_lg[3])


# Read in the fst file
s1_df_1 <- fst::read.fst(s1_t1[1])
s1_df_2 <- fst::read.fst(s1_t1[2])

# Generate SpatRasters
s1_rast_1 <- terra::rast(s1_t2[1])
s1_rast_2 <- terra::rast(s1_t2[2])

##################
# COPERNICUS DEM #
##################

# Run module for point geometry:
dem_t1 <- dem_process(sf_point, radius=100, uniqueID="Id", file_path = tempfile(pattern="dem_df_", tmpdir=tempdir()))

# Run module for polygon geometry: 
dem_t2 <- dem_process(sf_poly, file_path = tempfile(pattern="dem_tif_", tmpdir=tempdir()))

# Read in the fst file
dem_df <- fst::read.fst(dem_t1[1])

# Generate SpatRasters
dem_rast <- terra::rast(dem_t2)

##################
# ESA WORLDCOVER #
##################

# Run module for point geometry: 
esa_t1 <- esa_worldcover_process(sf_point, radius=100, start_dt="2019-01-01", end_dt="2020-01-01", uniqueID="Id",
                                 file_path=tempfile(pattern="esa_df_", tmpdir=tempdir()))

# Run module for polygon geometry: 
esa_t2 <- esa_worldcover_process(sf_poly, radius=200, start_dt="2019-01-01", end_dt="2022-01-01", uniqueID="Id",
                                 file_path=tempfile(pattern="esa_tif_", tmpdir=tempdir()))

# Read in the fst file:
esa_df1 <- fst::read.fst(esa_t1[1])

# Generate SpatRasters
esa_rast1 <- terra::rast(esa_t2[1])
esa_rast2 <- terra::rast(esa_t2[2])

############################
# ALOS FOREST / NON-FOREST #
############################

# Run module for point geometry: 
fnf_t1 <- alos_fnf_process(sf_point, radius=100, start_dt="2015-01-01", 
                           end_dt="2018-01-01", uniqueID="Id", file_path = tempfile(pattern="fnf_df_", tmpdir=tempdir()))

# Run module for polygon geometry:
fnf_t2 <- alos_fnf_process(sf_poly, start_dt="2015-01-01", end_dt="2018-01-01", file_path = tempfile(pattern="fnf_tif_", tmpdir=tempdir()))

# Read in the fst files: 
fnf_df1 <- fst::read.fst(fnf_t1[1])
fnf_df2 <- fst::read.fst(fnf_t1[2])

# Generate SpatRasters
fnf_rast1 <- terra::rast(fnf_t2[1])
fnf_rast2 <- terra::rast(fnf_t2[2])
