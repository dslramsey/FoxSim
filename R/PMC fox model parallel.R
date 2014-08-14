#' Sequential Population Monte Carlo Sampler
#'
#' \code{PMC.sampler} performs parameter estimation using sequential Monte Carlo
#' within an ABC framework.  A population of parameters are fitted to a sequence
#' of increasingly smaller tolerances to approximate the posterior distribution
#' of interest.
#' 
#' @param N size of the population of parameters to be estimated
#' @param x0 observations
#' @param SeqTol the sequence of (increasingly smaller) tolerances
#' @param priors list of prior distribution parameters
#' @param Data list containing other required inputs into \code{foxsim}
#' @param parallel logical indicating use of parallel computations
#' @param ncores integer indicating number of cores to use
#' @param logfile character indicating the filename to write out intermediate
#' results from call to parallel functions (i.e. parSapply)
#' @param save.post character indicating filename to save intermediate
#' tolerance results
#' 
#' @return a list with elements equal to the number of tolerances
#' for each tolerance level a list with the following elements
#' \code{theta} the posterior distribution with dimensions N x npar
#' \code{xr} accepted simulated road killed carcasses
#' \code{xs} accepted simulated hunter killed carcasses
#' \code{pop} accepted population size trajectory for each year 
#' (number of occupied cells)
#' \code{weights} weights matrix for this tolerance 
#' \code{elapsed} elapsed time (in seconds)
#' @export
#' 
PMC.sampler<- function(N, x0, SeqTol, priors, Data, parallel=F, ncores=NULL, logfile=NULL, save.post=NULL){
  post.list<- list()    
  if(parallel & !is.null(ncores)) {
    library(parallel)
    cl<- makeCluster(ncores,outfile=logfile)
    clusterEvalQ(cl, {
      library(Rcpp)
      library(RcppArmadillo)
      library(FoxSim)
    })
    #clusterExport(cl, varlist=c("ModelABC")) # export functions
  }
  start = Sys.time()
  # Use rejection sampling for first pass
  if(parallel)
    xx<- parSapply(cl, 1:N, ABC.reject, x0, SeqTol[1],priors, Data)
  else
    xx<- sapply(1:N, ABC.reject, x0, SeqTol[1], priors, Data)
  
    theta<- do.call('rbind',xx[1,])
    xr.sim<- t(sapply(xx[2,],function(x) x[(length(x)-(Data$eyear-2001)):length(x)])) #from 2001:end.year
    xs.sim<- t(sapply(xx[3,],function(x) x[(length(x)-(Data$eyear-2001)):length(x)])) #from 2001:end.year
    pop<- xx[4,]
  
    duration = difftime(Sys.time(), start, units = "secs")
    cat("Completed particles for tol ",SeqTol[1],"\n")
    cat("Time elapsed:", display.time(duration),"\n")
  
    w<- rep(1/N,N) 
    
    post.list[[paste0("Tol",SeqTol[1])]]<- list(theta=theta, xr=xr.sim, xs=xs.sim, pop=pop, weights=w, elapsed=duration)
  
    if(!is.null(save.post)) save(post.list,file=paste0(save.post,"tol",SeqTol[1],".RData"))
    # Weights are uniform
  
  # Now select a weighted sample for each tolerance updating weights and proposal at each iteration
  # proposals are multivariate normal
  
  for(j in 2:length(SeqTol)) {
      start = Sys.time()
      thetaold<- theta
      wold<- w
      VarCov = as.matrix(2 * cov.wt(thetaold, wold)$cov) ## Adjust std of proposal
      
      if(parallel)
        xx<- parSapply(cl, 1:N, ABC.weighted, thetaold, x0, SeqTol[j], VarCov, wold, priors, Data)
      else
        xx<- sapply(1:N, ABC.weighted, thetaold, x0, SeqTol[j], VarCov, wold, priors, Data)
      
      theta<- do.call('rbind',xx[1,])
      xr.sim<- t(sapply(xx[2,],function(x) x[(length(x)-(Data$eyear-2001)):length(x)])) #from 2001:end.year
      xs.sim<- t(sapply(xx[3,],function(x) x[(length(x)-(Data$eyear-2001)):length(x)])) #from 2001:end.year
      pop<- xx[4,]
      
      w<- calc.weights(theta, thetaold, wold, VarCov, priors)
      
      duration = difftime(Sys.time(), start, units = "secs")
      
      post.list[[paste0("Tol",SeqTol[j])]]<- list(theta=theta, xr=xr.sim, xs=xs.sim, pop=pop, weights=w, elapsed=duration)
      if(!is.null(save.post)) save(post.list,file=paste0(save.post,"tol",SeqTol[j],".RData"))      
      cat("Completed current tolerance ",SeqTol[j],"\n")
      cat("Time elapsed:", display.time(duration),"\n")
  }  
if(parallel) stopCluster(cl)                   
post.list

}
#-----------------------------------------------------------------------
#' Sequential Population Monte Carlo Sampler - Updater
#'
#' \code{PMC.update} Updates an existing posterior distribution \code{post}
#' produced by \code{PMC.sampler} using the specified tolerance.  All other
#' parameters are as for \code{PMC.sampler}.
#' 
#' @param post output from a single tolerance from \code{PMS.sampler}
#' @inheritParams PMC.sampler
#' 
#' @return a list with elements equal to the number of tolerances
#' for each tolerance level a list with the following elements
#' \code{theta} the posterior distribution with dimensions N x npar
#' \code{xr} accepted simulated road killed carcasses
#' \code{xs} accepted simulated hunter killed carcasses
#' \code{pop} accepted population size trajectory for each year 
#' (number of occupied cells)
#' \code{weights} weights matrix for this tolerance 
#' \code{elapsed} elapsed time (in seconds).
#' 
#' @seealso \code{\link{PMC.sampler}}
#' @export
#' 
PMC.update<- function(post, N, x0, SeqTol, priors, Data, parallel=F, ncores=NULL, logfile=NULL, save.post=NULL){
  post.list<- list()    
  if(parallel & !is.null(ncores)) {
    library(parallel)
    cl<- makeCluster(ncores,outfile=logfile)
    clusterEvalQ(cl, {
      library(Rcpp)
      library(RcppArmadillo)
      library(FoxSim)
    })
    #clusterExport(cl, varlist=c("ModelABC")) # export functions
  }
  # Update a weighted sample for each tolerance updating weights and proposal at each iteration
  # proposals are multivariate normal
  
  theta<- post$theta
  w<- post$weights
  
  for(j in 1:length(SeqTol)) {
    start = Sys.time()
    thetaold<- theta
    wold<- w
    VarCov = as.matrix(2 * cov.wt(thetaold, wold)$cov) ## Adjust std of proposal
    
    if(parallel)
      xx<- parSapply(cl, 1:N, ABC.weighted, thetaold, x0, SeqTol[j], VarCov, wold, priors, Data)
    else
      xx<- sapply(1:N, ABC.weighted, thetaold, x0, SeqTol[j], VarCov, wold, priors, Data)
    
    theta<- do.call('rbind',xx[1,])
    xr.sim<- t(sapply(xx[2,],function(x) x[(length(x)-(Data$eyear-2001)):length(x)]))
    xs.sim<- t(sapply(xx[3,],function(x) x[(length(x)-(Data$eyear-2001)):length(x)]))
    pop<- xx[4,]
    
    w<- calc.weights(theta, thetaold, wold, VarCov, priors)
    
    duration = difftime(Sys.time(), start, units = "secs")
    
    post.list[[paste0("Tol",SeqTol[j])]]<- list(theta=theta, xr=xr.sim, xs=xs.sim, pop=pop, weights=w, elapsed=duration)
    if(!is.null(save.post)) save(post.list,file=paste0(save.post,"tol",SeqTol[j],".RData"))      
    cat("Completed current tolerance ",SeqTol[j],"\n")
    cat("Time elapsed:", display.time(duration),"\n")
  }  
  if(parallel) stopCluster(cl)                   
  post.list
  
}
#------------------------------------------------------------------
#' Append PMC model objects
#'
#' \code{PMC.append} appends posterior samples from 2 or more objects
#' produced by \code{PMC.sampler}.  Only posterior samples with the same
#' tolerance are appended and this needs to be specified using \code{TolName}. 
#' 
#' @param x A list with at least 2 elements containing \code{PMC.sampler}
#' objects.  Each element of \code{x} must contain samples with at least 
#' one similar tolerance value.
#' @param TolName Character value indicating the tolerance value used for
#' appending samples.
#' 
#' @return A list with elements of \code{x} appended.\cr
#' \code{theta} - population parameter matrix.\cr  
#' \code{xr} - simulated road killed carcasses.\cr
#' \code{xs} - simulated hunter killed carcasses.\cr
#' \code{pop} - population size trajectory for each year 
#' (number of occupied cells).
#' 
#' @examples
#' PMC.append(list(fit1, fit2), "Tol2")
#'  
#' @seealso \code{\link{PMC.sampler}}
#' 
#' @export
#' 
PMC.append<- function(x, TolName) {
  # append posterior samples from objects with the same tolerance 
  # x is a list of objects produced by PMC.sampler  
  # TolName is the specific element (tolerance) of each object
  # to be appended
  if(!is.list(x)) stop("x is not a list")
  if(length(x) < 2) stop("Need at least 2 objects to append")
  n<- length(x)
  tmp<- x[[1]][[TolName]]
  theta<- tmp$theta
  xr<- tmp$xr
  xs<- tmp$xs
  pop<- tmp$pop
  for(i in 2:n) {
    tmp<- x[[i]][[TolName]]
    if(is.null(tmp)) stop("This tolerance does not occur in each object")
    theta<- rbind(theta, tmp$theta)
    xr<- rbind(xr, tmp$xr)
    xs<- rbind(xs, tmp$xs)
    pop<- c(pop, tmp$pop)
  }
  list(theta=theta,xr=xr,xs=xs,pop=pop)
}
#-----------------------------------------------------------------------
#' Model parser for PMC
#'
#' \code{ModelABC} is a convenience function that parses parameters
#' and ancillary data required for \code{foxsim}. This function is used
#' by \code{PMC.sampler} and is not normally called by the user
#' 
#' @param parm Vector of parameter proposals produced by \code{propose.theta}
#' for use in parameter estimation via \code{PMS.sampler}.
#' @param Data List of ancillary data required by \code{foxsim}.
#' 
#' @return a list with elements 
#' \code{xr} - simulated road killed carcasses,
#' \code{xs} - simulated hunter killed carcasses,
#' \code{pop} - population size trajectory for each year 
#' (number of occupied cells)
#'  
#' @seealso \code{\link{propose.theta}}
#' @export
#' 
ModelABC<- function(parm, Data) {
# Called by PMC.sampler  
  syear<- round(parm[1] + 1998)
  eyear<- Data$eyear 
  Parms<- list(syear=syear,eyear=eyear,psurv=parm[2],proad=parm[3],pshot=parm[4],Ryear=parm[5])
  # Fox cellular automata C++ function from library(FoxSim) 
  mod<- foxsim(Data$nr, Data$nc, Data$kdim, Data$hab.vec, Data$road.vec, Data$ipoints, Data$kern.list, Parms)
  
  list(xr=mod[[1]],xs=mod[[2]],pop=mod[[3]])
  
}
#-------------------------------------------------------------------------
#' Dispersal kernel 
#'
#' \code{Kernel2D} produces a square dispersal kernel matrix for use in
#' \code{foxsim}. Values are the probability of dispersal from the
#' origin cell (centre cell) to each cell in the matrix.  
#' 
#' @param eps cellsize in km used for distance calculations
#' @param max.disp Maximum dispersal distance (km). The resulting
#' matrix has dimensions \code{ceiling(max.disp/eps)}.
#' @param kfun Character vector indicating the kernel used to calculate the cell
#' probabilities. Possible values are \code{"circle"}, \code{"Gauss"},
#' \code{"uniform"}, or \code{"gamma"}.
#' @param sigma Parameter for the Gaussian, circle and uniform kernels i.e. for
#' \code{kfun="Gauss"}, \code{kfun="circle"} or \code{kfun="uniform"}.
#' @param shape,scale Shape and scale parameters for the gamma distribution
#' \code{kfun="gamma"}.
#' 
#' @return a square matrix with the dispersal probabilities for each cell
#' (normalised to sum to 1.0). 
#'
#' @examples
#' kernel2D(kfun="gamma") 
#' @export
#' 
kernel2D<- function(eps=3, max.disp=30, kfun=c("circle","Gauss","uniform","gamma"),sigma=5,shape=1,scale=3) {
  # 2D kernel density function fully vectorised
  calc.dist<- function(x,y){
    return(sqrt((x^2 + y^2))) 
  }  
  kfun<- match.arg(kfun)
  mx<-  ceiling(max.disp/eps)
  gxy <- -mx:mx * eps
  z <- outer(gxy, gxy, FUN=calc.dist) 
  switch(kfun,
         Gauss = {zprob<- dnorm(z, 0, sigma)},
         uniform = {zprob<- dunif(z,0,sigma)},
         gamma = {zprob<- dgamma(z,shape,1/scale)},
         circle = {zprob<- ifelse(z<=sigma,1,0)}
  )
  zprob[is.infinite(zprob)]<- 1
  znorm<- zprob/sum(zprob)
  znorm
}
#-------------------------------------------------------------------------------
#' Distance measure for ABC estimation
#'
#' \code{distm} calculates the sum of Euclidean distances between observed
#' and simulated data. Used by \code{PMC.sampler}.
#' 
#' @param xstar Vector of simulated observations.
#' @param x0 Vector of observed observations.
#' 
#' @return Sum of the Euclidean distance measure
#'  
#' @seealso \code{\link{PMC.sampler}}
#' @export
#' 
distm<- function(xstar, x0) { 
  sum(abs(xstar - x0))
}
#--------------------------------------------------------------------------------
#' Helper function required for ABC estimation
#'
#' \code{pre.pad} pads out the observation vector with zeros
#' so that it has the same length as the simulated observations.
#' Used by \code{distm}.
#' 
#' @param x Vector of simulated observations.
#' @param obs Vector of observed observations.
#' 
#' @return \code{obs} with zero observations added to the beginning 
#' of the vector so it is the same length as \code{x}.  This facilitates
#' distance measure calculation using \code{distm}.
#'  
#' @seealso \code{\link{PMC.sampler}} and \code{\link{distm}}
#' @export
#' 
pre.pad<- function(x, obs) {
  nzeros<- length(x) - length(obs)
  c(rep(0,nzeros),obs)
}
#--------------------------------------------------------------------------------
#' Rejection sampler for ABC estimation
#'
#' \code{ABC.reject} is a ABC rejection sampler used by 
#' \code{PMC.sampler}. Rejection sampling is performed for the 
#' first specified tolerance only and samples parameters from the prior
#' distribution.
#'  
#' @param i parameter index used by \code{sapply}
#' @param x0 List of observation vectors
#' @param tol Tolerance level (Euclidean distance)
#' @param priors List specifying parameters of the prior
#' distributions.
#' @param Data Ancillary data required by \code{foxsim}.
#' 
#' @return a list with elements
#' \code{theta} Vector of accepted parameters
#' \code{xr} Vector of accepted road killed carcasses
#' \code{xs} Vector of accepted hunter killed carcasses
#' \code{pop} Vector of accepted population size for each year 
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{ABC.weighted}}
#' @export
#' 
ABC.reject<- function(i, x0, tol, priors, Data) {
# ABC rejection algorithm sampling from priors 
  found<- FALSE
    while(!found){
      thetac=sample.priors(priors)
      xc=ModelABC(thetac, Data)
      xr0<- pre.pad(xc$xr,x0$xr)
      xs0<- pre.pad(xc$xs,x0$xs)
      if(distm(xc$xr,xr0) < tol & distm(xc$xs,xs0) < 1){
        found<- TRUE
        theta<- thetac
        xr<- xc$xr
        xs<- xc$xs
        pop<- xc$pop
      }
    }
  cat("completed particle ",i," for seq ",tol,"\n")
  list(theta,xr,xs,pop)
  }  
#------------------------------------------------------------------------
#' Weighted rejection sampler for ABC estimation
#'
#' \code{ABC.weighted} is a ABC rejection sampler used by 
#' \code{PMC.sampler} that uses weighted proposals.  Weighted
#' sampling is used for all tolerance levels after the first.
#' Weights are calcuated from the previous population.
#'  
#' @param i parameter index used by \code{sapply}
#' @param parms Parameter population from the previous tolerance step
#' @param x0 List of observation vectors
#' @param tol Tolerance level (Euclidean distance)
#' @param VarCov Variance/covariance matrix for multivariate
#' proposals used by \code{propose.theta}.
#' @param w Weights matrix used to select a particular particle
#' from the previous population as the basis for proposal.  Used
#' by \code{pick.particle}.
#' @param priors List specifying parameters of the prior
#' distributions.
#' @param Data Ancillary data required by \code{foxsim}.
#' 
#' @return a list with elements
#' \code{theta} Vector of accepted parameters
#' \code{xr} Vector of accepted road killed carcasses
#' \code{xs} Vector of accepted hunter killed carcasses
#' \code{pop} Vector of accepted population size for each year 
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{ABC.reject}}, 
#'  \code{\link{propose.theta}},  \code{\link{pick.particle}},
#'  \code{\link{calc.weights}}
#' @export
ABC.weighted<- function(i, parms, x0, tol, VarCov, w, priors, Data) {
# ABC proposal algorithm sampling from weighted posterior 
  found<- FALSE
  while(!found){
    thetao=pick.particle(parms, w)             
    thetac=propose.theta(thetao, VarCov, priors)     
    xc=ModelABC(thetac, Data) 
    xr0<- pre.pad(xc$xr,x0$xr)
    xs0<- pre.pad(xc$xs,x0$xs)
    if(distm(xc$xr,xr0) < tol & distm(xc$xs,xs0) < 1){
      found<- TRUE
      theta<- thetac
      xr<- xc$xr
      xs<- xc$xs
      pop<- xc$pop
    } 
  }
  cat("completed particle ",i," for seq ",tol,"\n")
  list(theta,xr,xs,pop)
}
#--------------------------------------------------------------------------------
#' Weights calculation for sequential Monte Carlo sampling
#'
#' \code{calc.weights} calculates the weights for each particle in
#' the population of parameters. The weights algorithm follows Beaumont et al. 
#' (2010) and is implemented here for multvariate normal proposals. Used in
#' \code{PMC.sampler}.
#'  
#' @param theta Current parameter population
#' @param thetaold Previous parameter population
#' @param w Vector of weights from previous population
#' @param VarCov Variance/covariance matrix used for multivariate
#' proposals.
#' @param priors List specifying parameters of the prior
#' distributions.
#' 
#' @return a vector of weights for the current population 
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{propose.theta}},
#'   \code{\link{pick.particle}}
#' @export
calc.weights<- function(theta, thetaold, w, VarCov, priors){
# calculate weights for posterior. Kernel is multivariate normal
  imat<- solve(VarCov)
  mvnorm<- function(mu, ivcov) {
    Z<- t(mu)   
    exp(-0.5*apply((ivcov %*% Z) * Z, 2, sum))
  }
  num<- apply((apply(theta, 1, prior.density, priors=priors)),2,prod)
  den<- apply(theta, 1, function(x) {sum(w * mvnorm(t(t(thetaold)-x),imat))})
  wstar<- num/den            
  wup<- wstar/sum(wstar)
  wup
}
#--------------------------------------------
#' Check for prior support
#'
#' \code{has.support} checks its arguments have non-zero prior
#' support.
#'  
#' @param parm Vector of parameter values
#' @param priors List specifying parameters of the prior
#' distributions.
#' 
#' @return logical indicating whether all elements of \code{parm} 
#' have non-zero prior support
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{propose.theta}}, 
#' \code{\link{prior.density}}
#' @export
has.support <- function(parm, prior) {
# test for prior support  
  test<- prior.density(parm, prior)
  if(any(test==0))
    return(FALSE)
  else
    return(TRUE)
}
#---------------------------------------------------------------------------------
#' Multivariate parameter proposals
#'
#' \code{propose.theta} proposes new parameter values based on 
#' perturbing the input vector \code{theta} using a multivariate 
#' normal distribution.  Multivariate normal proposals are 
#' implemented using a Cholsky decomposition of \code{vcov}.
#'  
#' @param theta Vector of parameter values to perturb
#' @param vcov Variance/covariance matrix for proposals
#' @param priors List specifying parameters of the prior
#' distributions.
#' 
#' @return Vector of proposed parameters the same length as
#' \code{theta} 
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{has.support}}
#' @export
propose.theta<- function(theta, vcov, priors) {
  # multivariate normal proposals using Cholsky
  n<- length(theta)
  test<- FALSE
  while(!test) {
    prop<- as.vector(theta + rbind(rnorm(n)) %*% chol(vcov))
    if(has.support(prop,priors)) test<- TRUE
  }
  prop
  }
#-------------------------------------------------------------------------
#' Weighted particle selection
#'
#' \code{pick.particle} selects a single particle (parameter
#' vector) from the population \code{param} at random using
#' probabilities based on \code{w}. Hence, particles with 
#' higher weights have higher chance of selection.  Selected 
#' particles are used as the basis for new proposals using
#' \code{propose.theta}.
#'  
#' @param param Population parameter matrix
#' @param vcov Variance/covariance matrix for proposals
#' @param priors List specifying parameters of the prior
#' distributions.
#' 
#' @return Vector of proposed parameters the same length as
#' \code{theta} 
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{has.support}}
#' @export
pick.particle <- function(param, w) {
  weight_cum = cumsum(w/sum(w))
  pos = 1:length(w)
  p = min(pos[weight_cum > runif(1)])
  param[p,]
}
#---------------------------------------------------------------------------------
#' Sample parameters from prior distributions
#'
#' \code{sample.priors} samples a parameter vector from the
#' supplied prior distribution specified in \code{priors}
#'  
#' @param priors List specifying parameters of the prior
#' distributions. The length of the list should be the same
#' as the number of parameters.  Each element of the list is
#' itself a list of length 3 indicating the name of the prior
#' distribution and its parameters e.g: \code{list("name", par1, par2)}.
#' 
#' \code{name}: Distribution name.  Currently supported distributions are \code{"uniform"}
#'  \code{"normal"}, \code{"lognormal"}, \code{"gamma"} and
#'  \code{"beta"}.\cr
#'  \code{par1, par2}: Parameters appropriate for the specified 
#'  distribution.
#' 
#' @return Vector of random parameters from the prior distributions
#' 
#' @seealso \code{\link{prior.density}}, \code{\link{has.support}}
#' @export
sample.priors<- function(priors) {
  n<- length(priors)
  param<- rep(NA,n)
  for(i in 1:n) {
    name<- priors[[i]][[1]]
    switch(EXPR = name, 
           uniform = {param[i]<- runif(1, priors[[i]][[2]],priors[[i]][[3]])}, 
           normal = {param[i]<- rnorm(1, priors[[i]][[2]],priors[[i]][[3]])}, 
           lognormal = {param[i]<- rlnorm(1, priors[[i]][[2]],priors[[i]][[3]])}, 
           gamma = {param[i]<- rgamma(1, priors[[i]][[2]],priors[[i]][[3]])},
           beta = {param[i]<- rbeta(1, priors[[i]][[2]],priors[[i]][[3]])})
    if(!name %in% c("uniform","normal","lognormal","gamma","beta"))
      stop("Prior distribution not recognised")
  }
  param
}
#---------------------------------------------------------------------------------
#' prior density of parameter vector
#'
#' \code{prior.density} evaluates the density of the supplied 
#' parameter values in \code{val} under the prior distribution.
#'  
#' @param val Parameter values to be evaluated
#'  @param priors List specifying parameters of the prior
#' distributions. The length of the list should be the same
#' as the number of parameters.  Each element of the list is
#' itself a list of length 3 indicating the name of the prior
#' distribution and its parameters e.g: \code{list("name", par1, par2)}.
#' 
#' \code{name}: Distribution name.  Currently supported distributions are \code{"uniform"}
#'  \code{"normal"}, \code{"lognormal"}, \code{"gamma"} and
#'  \code{"beta"}.\cr
#'  \code{par1, par2}: Parameters appropriate for the specified 
#'  distribution.
#'  @param logdens Logical indicating whether to calculate the log
#'  of the density estimate.
#' 
#' @return Vector of density estimates under the prior distributions
#' 
#' @seealso \code{\link{sample.priors}}, \code{\link{has.support}}
#' @export
prior.density<- function(val, priors, logdens=FALSE) {
  n<- length(priors)
  param<- rep(NA,n)
  if(length(val) != length(priors)) stop("parameter length mismatch")
  for(i in 1:n) {
  name<- priors[[i]][[1]]
  switch(EXPR = name, 
         uniform = {param[i]<- dunif(val[i], priors[[i]][[2]],priors[[i]][[3]],log=logdens)}, 
         normal = {param[i]<- dnorm(val[i], priors[[i]][[2]],priors[[i]][[3]],log=logdens)}, 
         lognormal = {param[i]<- dlnorm(val[i], priors[[i]][[2]],priors[[i]][[3]],log=logdens)}, 
         gamma = {param[i]<- dgamma(val[i], priors[[i]][[2]],prior[[i]][[3]],log=logdens)},
         beta = {param[i]<- dbeta(val[i], priors[[i]][[2]],priors[[i]][[3]],log=logdens)})
  if(!name %in% c("uniform","normal","lognormal","gamma","beta"))
    stop("Prior distribution not recognised")
  }
  param
}
#-----------------------------------------------------------------------------------
#' Display a pretty elapsed time
#'
#' \code{display.time} takes a \code{difftime} argument 
#' in seconds and calculates the elapsed time
#' in days, hours, minutes and seconds.
#' 
#' @param z Value in seconds produced by \code{difftime}
#'   
#' @return Character vector indicating the elapsed time
#' in days, hours, minutes and seconds
#' 
#' @seealso \code{\link{PMC.sampler}}
#' @export
display.time<- function(z) {
  # Display elapsed time in D:M:H:S format
  # z is a difftime object in seconds
  z<- as.numeric(z)
  d<- z %/% 86400 # days
  if(d > 0){
    h <- (z - 86400) %/% 60 %/% 60  # hours
    m <- (z - 86400) %/% 60 %% 60 # mins 
    s <- (z - 86400) %% 60 %% 60 # secs 
    if(d==1)
      paste0(d," day ",h," hours ",m," mins ",round(s)," secs")
    else
      paste0(d," days ",h," hours ",m," mins ",round(s)," secs")
  }
  else {
    h <- z %/% 60 %/% 60  # hours
    m <- z %/% 60 %% 60 # mins 
    s <- z %% 60 %% 60 # secs
    paste0(h," hours ",m," mins ",round(s)," secs")
  }
}