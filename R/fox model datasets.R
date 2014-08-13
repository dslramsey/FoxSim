#------------------------------------------------------------
#' Fox habitat suitability map of Tasmania.
#'
#' A vector representing fox habitat suitability on a 3 km 
#' raster suitable for input into \code{foxsim}.  Cell values
#' are coded (1) Suitable habitat, (0) unsuitable habitat and
#' (-1) water
#' 
#' 
#' @format A vector of size 12084 created from a raster file
#' by calling \code{as.double(as.matrix(raster))}
#' @name habitatvec
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
#' by calling \code{as.double(as.matrix(raster))}
#' @name roadvec
NULL
#-------------------------------------------------------------
#' List of introduction cells for foxes in Tasmania.
#'
#' A list of length 3 representing the possible introduction points
#' (cellnumbers) for foxes in Tasmania based on rumoured introductions 
#' in Longford, Oatlands and St Helens.  Cell numbers are all habitat
#' cells within 20 km of each location.#'  
#' 
#' @format A list of length 3 with each element consisting of 
#' a vector of cell numbers representing possible introduction
#' points.
#' @name incpoints
NULL