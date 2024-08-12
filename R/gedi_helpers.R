gedi_finder <- function(product) {
  # Define the base CMR granule search url, including the ORNL_CLOUD provider name and max page size
  cmr <- "https://cmr.earthdata.nasa.gov/search/granules.json?pretty=true&provider=ORNL_CLOUD&page_size=2000&concept_id="
  
  # Set up list where key is GEDI shortname and value is CMR concept ID. 
  # For now, we only have GEDI L3, but will potentially add GEDI L4B
  concept_ids <- list('GEDI03' = 'C2153683336-ORNL_CLOUD')
  
  # Check if the product is valid
  if (!product %in% names(concept_ids)) {
    stop("Invalid product name. Available products are: ", paste(names(concept_ids), collapse=","))
  }
  
  # CMR uses pagination for queries with more features returned than the page size. 
  page <- 1
  
  # Initialize a list to store and append granule links
  granules <- list() 
  
  # Send GET request to CMR granule search endpoint w/ product concept ID & page number
  cmr_response <- httr::GET(sprintf('%s%s&pageNum=%s', cmr, concept_ids[[product]], page))
  
  # Verify if the request submission was successful
  if (cmr_response$status_code==200) {
    
    # Send GET request to CMR granule search endpoint w/ product concept ID & page number
    cmr_url <- sprintf("%s%s&pageNum=%s", cmr, concept_ids[[product]],page)
    cmr_response <- content(GET(cmr_url))$feed$entry
    
    # If 2000 features are returned, move to the next page and submit another request, and append to the response
    while(length(cmr_response) %% 2000 == 0){
      page <- page + 1
      cmr_url <- sprintf("%s%s&pageNum=%s", cmr, concept_ids[[product]],page)
      cmr_response <- c(cmr_response, content(GET(cmr_url))$feed$entry)
    }
    
    # CMR returns more info than just the Data Pool links, below use for loop to grab each DP link, and add to list
    for (i in 1:length(cmr_response)) {
      granules[[i]] <- cmr_response[[i]]$links[[1]]$href
    }
    
    # Return the list of links
    return(granules)
  } else {
    # If the request did not complete successfully, print out the response from CMR
    print(content(GET(sprintf("%s%s&pageNum=%s", cmr, concept_ids[[product]],page)))$errors)
  }
}

gedi_get_latest_data <- function(granules) {
  # Error checks
  if (is.null(granules)) {
    stop("`granules` input is NULL. Provide a valid list of granules.")
  }
  
  if (!is.list(granules) || length(granules) == 0) {
    stop("`granules` must be a non-empty list.")
  }
  
  if (length(granules) < 5) {
    stop("`granules` list has fewer than 5 elements.")
  }
  
  # Subset output to last 5 granules (most recent datasets)
  granules_5 <- granules[(length(granules)-4):length(granules)]
  
  # Transform from list to character vector
  granules_chr <- as.character(granules_5)
  
  return(granules_chr)
}

gedi_raster <- function(file_list) {
  # Transform .tif files into one raster
  gedi_rast <- tryCatch({
    terra::rast(file_list)
  }, error = function(e) {
    stop("An error occured while creating the raster: ", e$message)
  })
  
  return(gedi_rast)
}

# Define function to download data. 
gedi_download <- function(files, netrc_path) {
  # Check for valid `files` input
  if (is.null(files) || length(files) ==0) {
    stop("The files input is NULL or empty. Provide a list of file URLs.")
  }
  
  # Check for valid `netrc_path` input
  if (!file.exists(netrc_path)) {
    stop(sprintf("The netrc file at %s does not exist. Provide a valid path.", netrc_path))
  }
  
  # Initialize the file list
  file_list <- NULL
  
  # Create a temporary directory
  temp_dir <- tempdir()
  
  # Loop through each file URL
  for (i in 1:length(files)) {
    # Keep original file name
    filename <- tail(strsplit(files[i], '/')[[1]], n=1) 
    
    # Define the full path for the file in the temp directory
    full_path <- file.path(temp_dir, filename)
    
    # Write file to disk (authenticating with netrc) using current directory/filename
    response <- tryCatch(
      httr::GET(files[i], write_disk(full_path, overwrite=TRUE),
                config(netrc=TRUE, netrc_file=netrc_path), set_cookies("LC"="cookies")),
      error = function(e) {
        warning(sprintf("Error occurred while downloading %s: %s", files[i], e$message))
        return(NULL)
      }
    )
    
    # Check if the response is NULL (error occurred)
    if (is.null(response)) {
      next
    }
    
    # Check to see if file downloaded correctly
    if (response$status_code == 200) {
      file_list <- c(file_list, full_path)
    } else {
      message(sprintf("%s not downloaded. Verify that your username and password are correct in %s. Status code: %s", 
                      filename, netrc_path, response$status_code))
    }
  }
  
  return(file_list)
}

gedi_reproject <- function(raster, aoi_layer) {
  # Check if inputs are valid
  if (!inherits(raster, "SpatRaster")) {
    stop("The input raster is not a valid SpatRaster object.")
  }
  
  if (!inherits(aoi_layer, "sf")) {
    stop("the input aoi_layer is not a valid sf object.")
  }
  
  # Check if the CRS is defined for both inputs:
  if (is.null(terra::crs(raster))) {
    stop("The raster object does not have a defined CRS.")
  }
  
  if (is.null(sf::st_crs(aoi_layer))) {
    stop("The aoi_layer object does not have a defined CRS.")
  }
  
  tryCatch({
    # Reproject AOI to match raster data CRS
    aoi_reproject <- sf::st_transform(aoi_layer, terra::crs(raster))
    
    # Add small buffer to account for distances with reprojection. 
    aoi_reproject_buffer <- sf::st_buffer(aoi_reproject, dist=50)
    
    # Crop the raster using the reprojected aoi
    rast_crop <- terra::crop(raster, aoi_reproject_buffer)
    
    # Get the crs of the org sf object 
    aoi_crs <- sf::st_crs(aoi_layer)$wkt
    
    # Reproject the raster to match the original aoi layer
    rast_crop_proj <- terra::project(rast_crop, aoi_crs)
    
    return(rast_crop_proj)
  }, error = function(e) {
    stop("An error occurred during the reprojection process: ", e$message)
  })
}