

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

point_buffer <- function(boundary_layer, radius=NULL, max_dist) {
  matching_layer <- boundary_layer

  # Check that radius argument is not null and is numerical
  if (!is.null(radius) && !is.numeric(radius)) {
    stop(deparse(substitute(radius)), " must be numeric or NULL when ", substitute(boundary_layer), " has POINT geometry type. \n")
  }

  if (is.null(radius)) {
    message("No radius defined. Will not apply a buffer to the point data. \n")
    return(matching_layer)
  }

  if (radius > max_dist) {
    warning(deparse(substitute(radius)), " exceeds ", max_dist, " m. Processing time and performance may be impacted. \n",
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

s2_scale <- function(raster, file_pth=NULL, counter=NULL) {
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
  
  if (!is.null(file_pth)) {
    output_filename <- paste0(file_pth,"_",counter, ".tif")
  }
  
  if (is.null(file_pth)) {
    output_filename <- tempfile(pattern="sentinel2_sr_image", tmpdir=tempdir(), fileext=".tif")
  }

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
  
  if (is.null(index_application) && is.null(index_names)) {
    return(NULL)
  }
  
  # Check if "nir09" is in the band_names vector and remove it if present
  if ("nir09" %in% names(band_names)) {
    band_names <- band_names[names(band_names) != "nir09"]
  }
  
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

extract_vals_df <- function(raster, buffer_layer, col_ID=NULL, file_name=NULL, add_dates=FALSE) {
  
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
  
  if (isTRUE(add_dates)) {
    colnames(aoi_df)[-1] <- paste(names(file_name), colnames(aoi_df)[-1], sep=" ")
  }
  
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

find_first_date <- function(start_date, 
                            end_date,
                            stac_search_fxn,
                            aoi_search,
                            band_mapping,
                            layers_search) {
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by="day")
  #return(dates)
  for (i in 1:length(dates)) {
    date_start <- as.character(dates[i])
    date_end <- as.character(dates[i]+1)
    
    # Call get_stac_data for the specific data
    result <- try(stac_search_fxn(
      aoi = aoi_search,
      start_date = date_start,
      end_date = date_end,
      asset_names = layers_search,
      stac_source = attr(band_mapping, "stac_source"),
      collection = attr(band_mapping, "collection_name"),
      query_function = attr(band_mapping, "query_function"),
      download_function = attr(band_mapping, "download_function"),
      sign_function = attr(band_mapping, "sign_function"),
      mask_band = attr(band_mapping, "mask_band"),
      mask_function = attr(band_mapping, "mask_function"),
      output_filename = tempfile(pattern="date_test", tmpdir=tempdir(), fileext=".tif")
    ), silent = TRUE)
    
    # Check if the result is successful
    if (!inherits(result,"try-error") && length(result) > 0) {
      return(as.character(dates[i]))
    }
  }
  return(NULL)
}


validate_geometry <- function(aoi_layer) {
  valid_types <- c("POLYGON", "POINT")
  if (!inherits(aoi_layer, "sf") || !any(valid_types %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object with ", paste(valid_types, collapse=", "), " geometry.\n")
  }
}

validate_method <- function(method, valid_methods, param_name) {
  if (!is.null(method) && !method[1] %in% valid_methods) {
    stop(deparse(substitute(method)), " for ", param_name, " must be one of ", paste(valid_methods, collapes=", "), ".\n")
  }
}


prepare_aoi <- function(aoi_layer, radius=NULL) {
  aoi_search <- aoi_layer
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    aoi_search <- point_buffer(aoi_layer, radius, max_dist=1000)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  check_bbox_area(aoi_search, max_area=5000000000)
  return(aoi_search)
}

select_layers <- function(layer_list=NULL, band_mapping) {
  # Initialize assets
  layer_vals <- band_mapping
  
  # Check if user selected specific layers
  if (!is.null(layer_list)) {
    if (!all(sapply(layer_list, is.character))) {
      stop("All elements in ", deparse(substitute(layer_list))," must be character strings.")
    }
    
    valid_layers <- c("blue", "coastal", "green", "nir", "nir08", "red", "rededge1", "rededge2", "rededge3", "swir16", "swir22")
    layer_match <- intersect(layer_list, valid_layers)
    
    if (length(layer_match) > 0) {
      layer_vals <- layer_vals[layer_match]
      no_match_layer <- setdiff(layer_list, layer_match)
      if (length(no_match_layer) > 0) {
        message("The following values in ", deparse(substitute(layer_list)), " are invalid and will be ignored: ", 
                paste(no_match_layer, collapse = ", "),"\n")
      }
    }
  }
  return(layer_vals)
}

add_scl_layer <- function(keep_mask=FALSE, composite, layers) {
  
  if (keep_mask == TRUE && !is.null(composite)) {
    message(deparse(substitute(keep_mask))," ignored when ", deparse(substitute(composite)), " is not `NULL`. \n")
  }
  
  # Add SCL band if composite method is NULL and keep_mask is TRUE
  if (keep_mask == TRUE && is.null(composite)) {
    layers <- c(layers, scl="SCL")
  }
  
  return(layers)
}

process_image_data <- function(stac_search_fxn,
                               aoi,
                               date_range,
                               band_mapping,
                               search_layers=band_mapping,
                               mask_name=NULL,
                               mask_function=NULL,
                               composite="median",
                               resample="bilinear",
                               max_date_range=NULL,
                               date_interval=NULL) {
  
  
  # If composite_method is not NULL:
  output_files <- NULL
  
  if (!is.null(composite)) {
    tryCatch({
      # Retrieve raster data
      # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
      output_files <- stac_search_fxn(
        aoi = sf::st_buffer(aoi, dist=50),
        start_date = date_range[1],
        end_date = date_range[2],
        asset_names = search_layers,
        stac_source = attr(band_mapping, "stac_source"),
        collection = attr(band_mapping, "collection_name"),
        query_function = attr(band_mapping, "query_function"),
        download_function = attr(band_mapping, "download_function"),
        sign_function = attr(band_mapping, "sign_function"),
        mask_band = mask_name,
        mask_function = mask_function,
        output_filename = tempfile(pattern=attr(band_mapping, "collection_name"), tmpdir=tempdir(), fileext=".tif"),
        composite_function = composite,
        gdalwarp_options = c("-r", resample, "-multi", "-overwrite", "-co",
                             "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
      )
    }, error = function(e) {
      message("Error encountered: ", e$message)
      message("No items were found for this combination of AOI and date range for ", attr(band_mapping, "collection_name"), ". \n")
    })
    if (!is.null(output_files)) {
      names(output_files) <- paste0(as.character(date_range[1]),"_",as.character(date_range[2]))
    }
    
  }
  
  if (is.null(composite)) {
    # Warning for long date range
    if ((as.Date(date_range[2])-as.Date(date_range[1])) > max_date_range) { 
      warning(deparse(substitute(composite))," is `NULL` and date range is greater than ", max_date_range, " days. Several files will be created which may impact processing time and performance. \n",
              immediate. = TRUE)
    }
    
    first_date <- find_first_date(start_date = date_range[1], 
                                  end_date = date_range[2],
                                  stac_search_fxn = stac_search_fxn,
                                  aoi_search = aoi,
                                  band_mapping = band_mapping,
                                  layers_search = search_layers)
    
    if (!is.null(first_date)) {
      date_sequence <- seq.Date(as.Date(first_date), as.Date(date_range[2]), by = paste(as.character(date_interval),"days"))
      tryCatch({
        # Download data for each date in the sequence
        output_files <- vapply(date_sequence,
                               function(date) {
                                 tryCatch({
                                   stac_search_fxn(
                                     aoi = sf::st_buffer(aoi, dist=50),
                                     start_date = as.character(date),
                                     end_date = as.character(date+(date_interval-1)),
                                     asset_names = search_layers,
                                     stac_source = attr(band_mapping, "stac_source"),
                                     collection = attr(band_mapping, "collection_name"),
                                     query_function = attr(band_mapping, "query_function"),
                                     download_function = attr(band_mapping, "download_function"),
                                     sign_function = attr(band_mapping, "sign_function"),
                                     mask_band = mask_name,
                                     mask_function = mask_function,
                                     output_filename = file.path(tempdir(), glue::glue("{date}.tif")),
                                     gdalwarp_options = c("-r", resample, "-multi", "-overwrite", "-co",
                                                          "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
                                   )
                                 }, error = function(e) {
                                   message("No data found for ", date)
                                   return("")
                                 })
                            
                               }, character(1)
        )
      }, error = function(e) {
        message("Error encountered: ", e$message)
        message("No items were found for this combination of AOI and date range for ", attr(band_mapping, "collection_name"), ". \n")
      })
      if (!is.null(output_files)) {
        names(output_files) <- as.character(date_sequence)
      }
    }
  }
  
  if (is.null(output_files)) {
    stop("No ", attr(band_mapping, "collection_name"), " items were found. \n")
  }
  
  return(output_files)
}

add_indices <- function(raster, file_name, count, indices_df) {
  
  # If no indices, scale data
  if (is.null(indices_df)) {
    out_image <- s2_scale(raster, file_pth=file_name, counter=count)
  }
  
  # Add indices
  if (!is.null(indices_df)) {
    out_image <- s2_scale(raster)
    if (nrow(indices_df) > 0) {
      indices_image <- rsi::calculate_indices(raster=raster, indices=indices_df, output_filename=tempfile(pattern="sentinel2_indices", tmpdir=tempdir(), fileext=".tif"))
      out_image <- rsi::stack_rasters(c(out_image, indices_image), output_filename=paste0(file_name,"_",count,".tif"))
    }
  }
  
  return(out_image)
}

df_out <- function(num_images, file_name, values_df, count) {
  # Write as fst file if only one image 
  if (num_images==1) {
    file_nm <- paste0(file_name, ".fst")
    fst::write_fst(values_df, path=file_nm)
    message("File saved to: ", file_nm)
    file_list <- file_nm
  }
  
  # Combine data frames for separate dates and write as fst file
  if (num_images > 1) {
    if (count==1) {
      combined_df <- values_df
    }
    else {
      combined_df <- cbind(combined_df, vals_df[-1])
    }
  }
  
  if (count==num_images && !is.null(combined_df)) {
    file_nm <- paste0(file_name, ".fst")
    fst::write_fst(combined_df, path=file_nm)
    message("File saved to: ", file_nm)
    file_list <- file_nm
  }
  
  return(file_list)
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

s1_dB <- function(raster, file_pth=NULL, counter=NULL) {
  # Create a copy of the raster to store the transformed values
  transformed_raster <- raster

  # Loop through all layers and apply the logarithmic transformation
  for (i in 1:terra::nlyr(raster)) {
    transformed_raster[[i]] <- 10*log10(raster[[i]])
    names(transformed_raster[[i]]) <- paste0(names(raster[[i]]), "_dB")
  }

  # Combine the two rasters so that we have backscatter data and dB
  combined <- c(raster, transformed_raster)
  
  if (!is.null(file_pth)) {
    output_filename <- paste0(file_pth,"_",counter, ".tif")
  }
  
  if (is.null(file_pth)) {
    output_filename <- tempfile(pattern="sentinel1_db_image", tmpdir=tempdir(), fileext=".tif")
  }
  
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
  area_sqm <- units::set_units(area, m^2)

  # Check if the area exceeds the max_area
  if (area_sqm > units::set_units(max_area, m^2)) {
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


