
#### slabR ####

#' Exception handling for using slabR
#'
#' This function catches the special cases when arrays are not 3D, which causes errors for Armadillo in SlabR.
#'
#' @details The function checks how many depth and time steps there are in an array of NEMO-ERSEM output. The function
#' then ensures that slabR is passed a 3D array, either by collapsing or applying the summary scheme over the time dimension,
#' or if there are too few dimensions (because a depth layer is being extracted) a simple join is performed to the scheme.
#'
#' @param file the name of a netcdf file containing the title variables.
#' @param scheme a summary scheme as expected by `array_w_mean()`.
#' @param count an optional vector of steps to subset along. See ncdf4 documentation.
#' @param collapse_days a logical switch setting whether to extract all days within a single netcdf file, or to take the mean across days.
#' @return A matrix with a column of group averages per title variable. If multiple days have been extracted these run sequentially within a column.
#' @family NEMO-ERSEM variable extractors
#' @export
NULL

process_array <- function(array, scheme = scheme, count = count, collapse_days = TRUE) {

  # If there are multiple depths and a single time step then slabR is fine as dim4 is dropped
  if(count[3] != 1 & count[4] == 1) {

    processed <- array_w_mean(array, scheme)                            # Summarise diffusivity according to scheme

  }

  # If there is only one depth layer and multiple time steps slabR breaks as we have matrices of data
  # This is a special case that crops up when trying to interpolate depth layers, but finding that layer already exists
  if(count[3] == 1 & count[4] != 1) {

    # Either Collect all time steps
    if(isFALSE(collapse_days)) {

      processed <- reshape2::melt(array, varnames=c('x', 'y', "t"), value.name = "processed") %>%
        right_join(scheme) %>%
        `$`(processed)
    }

    # Or Average across time steps
    if(isTRUE(collapse_days)) {

      processed <- reshape2::melt(rowMeans(na.rm = T, array, dims = 2), varnames=c('x', 'y'), value.name = "processed") %>%
        right_join(scheme) %>%
        `$`(processed)
    }

  }

  # If there are multiple time steps and depths slabr breaks as we have more than 3 dimensions
  if(count[3] != 1 & count[4] != 1) {

    # Either Collect all time steps by using apply on the 4th dimension
    if(isFALSE(collapse_days)) {

      processed <- as.vector(apply(array, 4, array_w_mean, scheme))   # Summarise diffusivity according to scheme

    }

    # Or Average across time steps
    if(isTRUE(collapse_days)) {

      processed <- array_w_mean(rowMeans(na.rm = T, array, dims = 3), scheme) # Summarise diffusivity according to scheme

    }
  }

  return(processed)

}

#' Extract summaries from Mission Atlantic NEMO-ERSEM arrays using slabR
#'
#' These functions read in target variables from NEMO-ERSEM model outputs and return weighted averages
#' according to a summary scheme.
#'
#' @details Each variable of interest in a netcdf file is imported, only reading within an x/y window specified
#' with `start` and `count`. The values are then passed to `array_w_mean()` to summarise according to a scheme.
#'
#' As the Mission Atlantic data files contain multiple time steps, these functions calculate the mean monthly
#' array BEFORE using `array_w_mean()` for the spatial operation.
#'
#' Each model output file type contains a different variable, extracted by the relevant function variant:
#'
#' | File type &nbsp; &nbsp; &nbsp; &nbsp;| Variables    |
#' |--------------|-------------|
#' | thetao_con   | Temperature. Celsius|
#' | so_abs       | Salinity.|
#' | uo           | Zonal currents. m/s|
#' | vo           | Meridional currents. m/s|
#' | wo           | Vertical velocity.|
#' | difvho       | Vertical diffusivity. m2/s|
#' | R1_n         | Dissolved organic nitrogen. mmolN.m-3|
#' | O2_o         | Oxygen. mmol.m-3|
#' | N4_n         | Ammonium. mmolN.m-3|
#' | N3_n         | Nitrate. mmolN.m-3|
#' | RP_n_result  | Detritus (Particulate organic nitrogen). mmolN.m-3|
#' | B1_n         | Bacterial nitrogen. mmolN.m-3|
#' | P1_n         | Diatom nitrogen. mmolN.m-3|
#' | P234_n_result| Other phytoplankton nitrogen. mmolN.m-3 |
#'
#' Some function variants have different arguments:
#'
#' The Mission Atlantic dataset was interpolated onto a regular lat-lon and depth grid. There is no need for any grid_w arguments in contrast to the nemomedusR package.
#'
#' @md
#' @param path the path to the NEMO-ERSEM model outputs.
#' @param file the name of a netcdf file containing the title variables.
#' @param scheme a summary scheme as expected by `array_w_mean()`.
#' @param start an optional vector of indices to start subsetting at. See ncdf4 documentation.
#' @param count an optional vector of steps to subset along. See ncdf4 documentation.
#' @param collapse_days a logical switch setting whether to extract all days within a single netcdf file, or to take the mean across days.
#' @param ... soaks up unused function arguments passed by the wrapper functions handling file architecture.
#' @return A matrix with a column of group averages per title variable for a single day.
#' @family NEMO-ERSEM variable extractors
#' @name extractors_MA
NULL

#' @rdname extractors_MA
#' @export
get_so_abs_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "so_abs", start, count)                # Extract an array of salinities
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Salinity = array_w_mean(rowMeans(na.rm = T,nc_saline, dims = 3), scheme))                                # Summarise salinity according to scheme

  return(all)
}

#' @rdname extractors_MA
#' @export
get_thetao_con_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_temp <- ncdf4::ncvar_get(nc_raw, "thetao_con", start, count)              # Extract an array of temperatures
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Temperature = array_w_mean(rowMeans(na.rm = T,nc_temp, dims = 3), scheme))                               # Summarise temperature according to scheme

    return(all)
}

#' @rdname extractors_MA
#' @export
get_R1_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_DON <- ncdf4::ncvar_get(nc_raw, "R1_n", start, count)                       # Extract an array for the variable
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
    DON = array_w_mean(rowMeans(na.rm = T,nc_DON, dims = 3), scheme))                                    # summarise Dissolved Organic Nitrogen
    return(all)
}

#' @rdname extractors_MA
#' @export
get_O2_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                     # Open up a netcdf file to see it's raw contents (var names)
  nc_O2 <- ncdf4::ncvar_get(nc_raw, "O2_o", start, count)                          # Extract an array for the variable
  ncdf4::nc_close(nc_raw)                                                          # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
    O2 = array_w_mean(rowMeans(na.rm = T,nc_O2, dims = 3), scheme))                                              # summarise dissolved oxygen according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_N4_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_NH4 <- ncdf4::ncvar_get(nc_raw, "N4_n", start, count)
  ncdf4::nc_close(nc_raw)                                                          # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
    NH4 = array_w_mean(rowMeans(na.rm = T,nc_NH4, dims = 3), scheme))                                            # summarise Dissolved ammonium according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_N3_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_NO3 <- ncdf4::ncvar_get(nc_raw, "N3_n", start, count)                       # Extract an array for the variable
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    NO3 = array_w_mean(rowMeans(na.rm = T,nc_NO3, dims = 3), scheme))                                        # summarise dissolved nitrate according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_RP_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
    nc_DET <- ncdf4::ncvar_get(nc_raw, "RP_n_result", start, count)
  ncdf4::nc_close(nc_raw)                                                          # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
        Detritus = array_w_mean(rowMeans(na.rm = T,nc_DET, dims = 3), scheme))                                   # summarise Detritus according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_B1_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_B1 <- ncdf4::ncvar_get(nc_raw, "B1_n", start, count)                        # Extract an array for the variable
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    Bacteria = array_w_mean(rowMeans(na.rm = T,nc_B1, dims = 3), scheme))                                    # summarise Bacterial nitrogen according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_P1_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_PHD <- ncdf4::ncvar_get(nc_raw, "P1_n", start, count)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    Diatoms = array_w_mean(rowMeans(na.rm = T,nc_PHD, dims = 3), scheme))                                    # summarise Diatoms according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_P234_n_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_PHO <- ncdf4::ncvar_get(nc_raw, "P234_n_result", start, count)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    Other_phytoplankton = array_w_mean(rowMeans(na.rm = T,nc_PHO, dims = 3), scheme))                        # summarise Phytoplankton according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_difvho_slabR   <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), collapse_days = TRUE, ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_dif <- ncdf4::ncvar_get(nc_raw, "difvho", start, count)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  # all <- cbind(                                                                # Bind as matrix
  #   Vertical_diffusivity = array_w_mean(rowMeans(na.rm = T,nc_dif, dims = 2), scheme))                       # Summarise diffusivity according to scheme

  # Vertical_diffusivity <- reshape2::melt(rowMeans(na.rm = T,nc_dif, dims = 2), varnames=c('x', 'y'), value.name = "Vertical_diffusivity") %>%
  #   right_join(scheme)

  all <- cbind(                                                                # Bind as matrix
    Vertical_diffusivity = process_array(nc_dif, collapse_days = collapse_days))

  return(all)
}

#' @rdname extractors_MA
#' @export
get_wo_slabR   <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "wo", start, count)                       # Extract an array for the variable
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  # all <- cbind(                                                                # Bind as matrix
  #   Vertical_velocity = array_w_mean(rowMeans(na.rm = T,nc_vel, dims = 2), scheme))                          # Summarise vertical velocity according to scheme

  # Vertical_velocity <- reshape2::melt(rowMeans(na.rm = T,nc_vel, dims = 2), varnames=c('x', 'y'), value.name = "Vertical_velocity") %>%
  #   right_join(scheme)

  # all <- cbind(                                                                # Bind as matrix
  #   Vertical_velocity = Vertical_velocity$Vertical_velocity)

  all <- cbind(                                                                # Bind as matrix
    Vertical_velocity = process_array(nc_vel, collapse_days = collapse_days))

  return(all)
}

#' @rdname extractors_MA
#' @export
get_vo_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vo", start, count)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(Meridional = array_w_mean(rowMeans(na.rm = T,nc_merid, dims = 3), scheme))                  # summarise meridional currents according to scheme
  return(all)
}

#' @rdname extractors_MA
#' @export
get_uo_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "uo", start, count) # Pull zonal current
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(Zonal = array_w_mean(rowMeans(na.rm = T,nc_zonal, dims = 3), scheme))                       # summarise zonal currents according to scheme
  return(all)
}

