% Generated by roxygen2 (4.0.2): do not edit by hand
\name{PMC.update}
\alias{PMC.update}
\title{Sequential Population Monte Carlo Sampler - Updater}
\usage{
PMC.update(post, N, x0, SeqTol, priors, Data, parallel = FALSE,
  ncores = NULL, logfile = NULL, save.post = NULL)
}
\arguments{
\item{post}{output from a single tolerance from \code{PMS.sampler}}

\item{N}{size of the population of parameters to be estimated}

\item{x0}{observations}

\item{SeqTol}{the sequence of (increasingly smaller) tolerances}

\item{priors}{list of prior distribution parameters}

\item{Data}{list containing other required inputs into \code{foxsim}}

\item{parallel}{logical indicating use of parallel computations}

\item{ncores}{integer indicating number of cores to use}

\item{logfile}{character indicating the filename to write out intermediate
results from call to parallel functions (i.e. parSapply)}

\item{save.post}{character indicating filename to save intermediate
tolerance results}
}
\value{
a list with elements equal to the number of tolerances
for each tolerance level a list with the following elements
\code{theta} the posterior distribution with dimensions N x npar
\code{xr} accepted simulated road killed carcasses
\code{xs} accepted simulated hunter killed carcasses
\code{pop} accepted population size trajectory for each year
(number of occupied cells)
\code{weights} weights matrix for this tolerance
\code{elapsed} elapsed time (in seconds).
}
\description{
\code{PMC.update} Updates an existing posterior distribution \code{post}
produced by \code{PMC.sampler} using the specified tolerance.  All other
parameters are as for \code{PMC.sampler}.
}
\seealso{
\code{\link{PMC.sampler}}
}

