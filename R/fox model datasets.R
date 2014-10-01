#------------------------------------------------------------
#' Fox habitat suitability map of Tasmania.
#'
#' A matrix representing fox habitat suitability on a 3 km 
#' raster suitable for input into \code{foxsim}.  Cell values
#' are coded (1) Suitable habitat, (0) unsuitable habitat and
#' (-1) water
#' 
#' 
#' @format A matrix of size 114 x 106 created from a raster file
#' by calling \code{as.matrix(raster)}\cr
#' @name habitatmat
NULL
#------------------------------------------------------------
#' Road network for Tasmania.
#'
#' A vector representing the main road network for Tasmania 
#' on a 3 km raster suitable for input into \code{foxsim}.  Cell values
#' are coded (1) road presence, (0) road absence
#' 
#' 
#' @format A vector of size 12084 created from a raster file
#' by calling \code{as.double(as.matrix(raster))}\cr
#' The dimensions of the original raster are stored in the 
#' \code{dims} attibute (accessed by \code{attributes(roadvec)}
#' @name roadmat
NULL
#-------------------------------------------------------------
#' List of introduction cells for foxes in Tasmania.
#'
#' A list of length 2 representing the possible introduction points
#' (cell references) for foxes in Tasmania based on possible introductions
#' through the ports at Burnie, Devonport and Hobart as well as rumoured 
#' introductions in Longford, Oatlands and St Helens.  Each element is a matrix
#' of cell references (row and column number) of cells within 3 km of the Port
#' locations and 30 km of the rumoured introductions.  
#' 
#' @format A list of length 2 with elements consisting of 
#' a list of length 3 containing a matrix of cell locations, one
#' for each possible introduction point.
#' @name incpoints
NULL
#-------------------------------------------------------------
#' Fox carcass discoveries in Tasmania.
#'
#' A list of length 2 representing the number of fox carcasses
#' discovered in Tasmania since 2001 either through road kills
#' or hunter kills.  
#' 
#' @format A list of length 2 with the following elements\cr
#' xr - number of road killed carcasses each year.\cr 
#' xs - number of hunter killed carcasses each year.\cr
#' Years are numbered from 2001.
#' 
#' @name carcass_obs
NULL