\name{foxsim}
\alias{foxsim}
\title{Tasmanian Fox Cellular Automata}
\usage{
foxsim(nr, nc, ksize, x, roads, incpoints, Kern, parms)
}

\arguments{
\item{nr}{number of rows in the habitat map \code{x}}

\item{nc}{number of columns in the habitat map \code{x}}

\item{ksize}{dimension of the (square) dispersal kernel matrix}

\item{x}{habitat map vector of size \code{nr} * \code{nc}. Values for each
cell are coded 1 - suitable habitat, 0 - unsuitable habitat or (-1) - water}

\item{roads}{road map vector of size \code{nr} * \code{nc}. Values for each
cell are coded 1 - road or 0 - no road}

\item{incpoints}{list containing a vector of cells for possible fox introductions.
The length of the list indicates the number of simultaneous introductions}

\item{Kern}{List of possible dispersal kernels as a square matrix}

\item{parms}{list of parameters of the population and observation processes
including\cr
\code{syear}: introduction year,\cr
\code{eyear}: end year,\cr
\code{psurv}: cell survial probability,\cr
\code{Ryear}: cell reproduction rate,\cr
\code{proad}: probability of road kill,\cr
\code{pshot}: probability of hunter kill}
}
\value{
a list with the following elements
\code{xr} simulated number of road killed carcasses for each year,
\code{xs} simulated number of hunter killed caracsses for each year,
\code{pop} simulated population size each year (number of occupied cells)
}
\description{
\code{foxsim} simulates a cellular automata of an invading fox population
in Tasmania and is designed to be fitted to observations of fox carcasses.
The model is fully coded in C++ for speed.
}