\name{FoxSim-package}
\alias{FoxSim-package}
\docType{package}
\title{
FoxSim - Tasmanian Fox Cellular Automata
}
\description{
A cellular automata to simulate fox (Vulpes vulpes) invasion
in Tasmania coded in C++.  In addition, the package includes
functions to fit the model to observations of carcasses using
sequential Monte Carlo ABC methods.
}
\details{
\tabular{ll}{
Package: \tab FoxSim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-24\cr
License: \tab Free as in free beer\cr
  }
}
\author{
Dave Ramsey - Arthur Rylah Institute

Maintainer: Dave Ramsey <david.ramsey@depi.vic.gov.au>
}
\references{
Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009) Adaptive approximate Bayesian computation. Biometrika,(96),pp 983-990. 
}
\examples{
data(habitatmat)
data(roadmat)
data(incpoints)
data(carcass_obs)
data(carcass_locs)
data(spotlocs)

kern1<- kernel2D(eps=3, kfun="gamma", shape=1, scale=3, max.disp=30)
kern2<- kernel2D(eps=3, kfun="uniform", sigma=30, max.disp=30)

kern<- list(kern1, kern2)

MyData<- list(habitat.mat=habitatmat,road.mat=roadmat,ipoints=incpoints, kern.list=kern, syear=1995, eyear=2014, carcass.road=carcass_locs$xr, carcass.shot=carcass_locs$xs, spotlocs=spotlocs, ncell=10)


priors<-list(list("uniform",-0.5,2.49),
             list("uniform",-0.5,2.49),
             list("uniform",-0.5,11.49),
             list("uniform",-0.5,6.49),
             list("beta",1,1),
             list("beta",17,83),
             list("beta",17,83),
             list("beta",10,90),
             list("uniform",0.1,2))

Tol<- c(6,5)

#fit<- PMC.sampler(N=10, carcass_obs, Tol, priors, MyData, parallel=FALSE)
}
