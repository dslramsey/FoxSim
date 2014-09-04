\name{matchspatial}
\alias{matchspatial}
\title{Spatial matching of cell locations}
\usage{
matchspatial(n, m, locs, x, ncells, state)
}

\arguments{
\item{n}{number of rows in the spatial map \code{x}}

\item{m}{number of columns in the spatial map \code{x}}

\item{locs}{matrix of row  and column numbers of cell locations for matching}

\item{x}{vector of size \code{n} * \code{m} with values coded
0 or \code{state} representing cells that will be spatially
matched against the locations in \code{locs}}

\item{ncells}{the radius for a successful spatial match given as the number 
of cells}

\item{state}{the cell value indicating which cells will be used for matching}

}
\value{
a logical vector of size \code{n} * \code{m} with values set to 
\code{TRUE} for cell locations in \code{x} that are within \code{ncells}
distance of locations given in \code{locs} (i.e. within
\code{ncells} of \code{x[locs(_,0) + n * locs(_,1)]}}
\description{
\code{matchspatial} performs spatial matching of cells in the vector \code{x} 
against locations given in the matrix \code{locs}. Used for spatial matching 
of simulated observations with observed locations.  Coded in C++ for speed.
}