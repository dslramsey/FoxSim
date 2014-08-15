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
data(habitatvec)
data(roadvec)
data(incpoints)
data(carcass_obs)
nd<- attributes(habitatvec)$dims

kern1<- kernel2D(eps=3, kfun="gamma", shape=1, scale=3, max.disp=30)
kern2<- kernel2D(eps=3, kfun="uniform", sigma=30, max.disp=30)

kern<- list(kern1, kern2)
kdim<- dim(kern1)[1] 

MyData<- list(nr=nd[1],nc=nd[2],kdim=kdim, hab.vec=habitatvec, road.vec=roadvec,  
              ipoints=incpoints, kern.list=kern, eyear=2014)

priors<- list(list("uniform",-3.49,3.49),
              list("beta",1,1),
              list("beta",17,83),
              list("beta",17,83),
              list("uniform",0.1,2))

Tol<- c(6,5)

fit<- PMC.sampler(N=10, carcass_obs, Tol, priors, MyData, parallel=FALSE)
}
