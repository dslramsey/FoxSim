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
#' @param kern type of proposals for peturbing variables. Choices are
#' "multi" or "uni" indicating multivariate or componentwise repectively
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
PMC.sampler<- function(N, x0, SeqTol, priors, Data, kern = c("multi","uni"), parallel=FALSE, ncores=NULL, logfile=NULL, save.post=NULL){
  seed<- floor(runif(1, 1, 10000))
  if(parallel & !is.null(ncores)) {
    library(parallel)
    cl<- makePSOCKcluster(ncores,outfile=logfile)
    on.exit(stopCluster(cl))
    clusterSetRNGStream(cl, seed)
    clusterEvalQ(cl, {
      library(Rcpp)
      library(RcppArmadillo)
      library(FoxSim)
    })
  }
  
  kern<- match.arg(kern, c("multi","uni"))
  start = Sys.time()
  cTol<- SeqTol$ctol
  sTol<- SeqTol$stol
  nseq<- length(cTol)
  
  # Use rejection sampling for first pass
  if(parallel)
    xx<- clusterApplyLB(cl, 1:N, ABC.reject, x0, cTol[1], sTol[1], priors, Data)
  else
    xx<- lapply(1:N, ABC.reject, x0, cTol[1], sTol[1], priors, Data)
  
    theta<- t(sapply(xx, function(x) rbind(x[[1]])))
    xr.sim<- t(sapply(xx, function(x) rbind(x[[2]])))
    xs.sim<- t(sapply(xx, function(x) rbind(x[[3]])))
    xspot.sim<- t(sapply(xx, function(x) rbind(x[[4]])))
    xscats.sim<- t(sapply(xx, function(x) rbind(x[[5]])))
    fpscats.sim<- t(sapply(xx, function(x) rbind(x[[6]])))
    poploc<- lapply(xx, function(x) x[[7]])
    pop<- lapply(poploc, function(x) sapply(x, function(y) sum(y==2)))
     
    duration = difftime(Sys.time(), start, units = "secs")
    cat("Completed particles for tol ",cTol[1], " and ",sTol[1],"\n")
    cat("Time elapsed:", display.time(duration),"\n")
    
    w<- rep(1/N,N) 
    
    post<- list(theta=theta, xr=xr.sim, xs=xs.sim, xspot=xspot.sim, xscats=xscats.sim, fpscats=fpscats.sim, 
                pop=pop, weights=w, cTol=cTol[1], sTol=sTol[1], elapsed=duration)
  
    if(!is.null(save.post)) saveRDS(post,file=paste0(save.post,"tol",1,".rds"))
    if(!is.null(save.post)) saveRDS(poploc, file=paste0(save.post,"pop",1,".rds"))   
    # Weights are uniform
   rm(xx,poploc,xr.sim,xs.sim,xspot.sim,xscats.sim,fpscats.sim)
  # Now select a weighted sample for each tolerance updating weights and proposal at each iteration
  # proposals are multivariate normal
  if(nseq > 1) {
  for(j in 2:nseq) {
      start = Sys.time()
      thetaold<- theta
      wold<- w
      if(identical(kern, "multi"))
        VarCov = as.matrix(2 * cov.wt(thetaold, wold)$cov) ## multivariate proposals
      else
        VarCov = as.matrix(2 * diag(var.wt(thetaold, wold))) ## componentwise proposals
      
      if(parallel)
        xx<- clusterApplyLB(cl, 1:N, ABC.weighted, thetaold, x0, cTol[j], sTol[j], VarCov, wold, priors, Data)
      else
        xx<- lapply(1:N, ABC.weighted, thetaold, x0, cTol[j], sTol[j], VarCov, wold, priors, Data)
      
      theta<- t(sapply(xx, function(x) rbind(x[[1]])))
      xr.sim<- t(sapply(xx, function(x) rbind(x[[2]])))
      xs.sim<- t(sapply(xx, function(x) rbind(x[[3]])))
      xspot.sim<- t(sapply(xx, function(x) rbind(x[[4]])))
      xscats.sim<- t(sapply(xx, function(x) rbind(x[[5]])))
      fpscats.sim<- t(sapply(xx, function(x) rbind(x[[6]])))
      poploc<- lapply(xx, function(x) x[[7]])
      pop<- lapply(poploc, function(x) sapply(x, function(y) sum(y==2)))
      w<- calc.weights(theta, thetaold, wold, VarCov, priors)
      
      duration = difftime(Sys.time(), start, units = "secs")
      
      post<- list(theta=theta, xr=xr.sim, xs=xs.sim, xspot=xspot.sim, xscats=xscats.sim, 
                  fpscats=fpscats.sim, pop=pop, weights=w, cTol=cTol[j], sTol=sTol[j], elapsed=duration)
      
      if(!is.null(save.post)) saveRDS(post,file=paste0(save.post,"tol",j,".rds")) 
      if(!is.null(save.post)) saveRDS(poploc, file=paste0(save.post,"pop",j,".rds"))
      
      cat("Completed current tolerance ",cTol[j]," and ",sTol[j],"\n")
      cat("Time elapsed:", display.time(duration),"\n")
      rm(xx,poploc,xr.sim,xs.sim,xspot.sim,xscats.sim) 
    }
  }
          
post

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
PMC.update<- function(post, N, x0, SeqTol, priors, Data, kern = c("multi","uni"), parallel=FALSE, ncores=NULL, logfile=NULL, save.post=NULL){
  post.list<- list()    
  seed<- floor(runif(1, 1, 10000))
  if(parallel & !is.null(ncores)) {
    library(parallel)
    cl<- makePSOCKcluster(ncores,outfile=logfile)
    on.exit(stopCluster(cl))
    clusterSetRNGStream(cl, seed)
    clusterEvalQ(cl, {
      library(Rcpp)
      library(RcppArmadillo)
      library(FoxSim)
    })
  }
  
  # Update a weighted sample for each tolerance updating weights and proposal at each iteration
  # proposals are multivariate normal
  kern<- match.arg(kern, c("multi","uni"))
  
  theta<- post$theta
  w<- post$weights

  cTol<- SeqTol$ctol
  sTol<- SeqTol$stol
  nseq<- length(cTol)
  
  for(j in 1:nseq) {
    start = Sys.time()
    thetaold<- theta
    wold<- w
    if(identical(kern, "multi"))
      VarCov = as.matrix(2 * cov.wt(thetaold, wold)$cov) ## multivariate proposals
    else
      VarCov = as.matrix(2 * diag(var.wt(thetaold, wold))) ## componentwise proposals
    
    if(parallel)
      xx<- clusterApplyLB(cl, 1:N, ABC.weighted, thetaold, x0, cTol[j], sTol[j], VarCov, wold, priors, Data)
    else
      xx<- lapply(1:N, ABC.weighted, thetaold, x0, cTol[j], sTol[j], VarCov, wold, priors, Data)
    
    theta<- t(sapply(xx, function(x) rbind(x[[1]])))
    xr.sim<- t(sapply(xx, function(x) rbind(x[[2]])))
    xs.sim<- t(sapply(xx, function(x) rbind(x[[3]])))
    xspot.sim<- t(sapply(xx, function(x) rbind(x[[4]])))
    xscats.sim<- t(sapply(xx, function(x) rbind(x[[5]])))
    fpscats.sim<- t(sapply(xx, function(x) rbind(x[[6]])))
    poploc<- lapply(xx, function(x) x[[7]])
    pop<- lapply(poploc, function(x) sapply(x, function(y) sum(y==2)))
    
    w<- calc.weights(theta, thetaold, wold, VarCov, priors)
    
    duration = difftime(Sys.time(), start, units = "secs")
    
    post<- list(theta=theta, xr=xr.sim, xs=xs.sim, xspot=xspot.sim, xscats=xscats.sim, 
                fpscats=fpscats.sim, pop=pop, weights=w, ctol=cTol[j], stol=sTol[j], elapsed=duration)
    
    if(!is.null(save.post)) saveRDS(post,file=paste0(save.post,"tol",j,".rds"))  
    if(!is.null(save.post)) saveRDS(poploc, file=paste0(save.post,"pop",j,".rds"))
    cat("Completed current tolerance cTol= ",cTol[j]," and sTol= ",sTol[j],"\n")
    cat("Time elapsed:", display.time(duration),"\n")
    rm(poploc,xx) # cleanup
  }  
  if(parallel) stopCluster(cl)                   
  post
  
}
#---------------------------------------------------------------------
#' Append PMC model objects saved in rds files
#'
#' \code{PMC.appendFiles} appends posterior samples from 2 or more 
#' \code{PMC.sampler} objects saved to rds binary files.  Only 
#' posterior samples with the same tolerance are appended and this 
#' needs to be specified using \code{TolName}. 
#' 
#' @param x A character vector of length at least 2 containing the relative
#' or full paths and filenames for all saved files.  files should be binary 
#' files saved from \code{PMC.sampler} (by setting \code{save.post}) and all 
#' files should contain at least one similar tolerance value.
#' @param TolName Character value indicating the tolerance value used for
#' appending samples.
#' @param fileout character value of the filename to write the appended object 
#' to a file. Default is \code{NULL} which does not write anything.
#' 
#' @return A list with elements of \code{x} appended.\cr
#' \code{theta} - population parameter matrix.\cr  
#' \code{xr} - simulated road killed carcasses.\cr
#' \code{xs} - simulated hunter killed carcasses.\cr
#' \code{xspot} - simulated spotlight detections. \cr
#' \code{xscats} - simulated scat detections. \cr
#' \code{pop} - population size trajectory for each year 
#' (number of occupied cells).
#' 
#' @examples
#' filenm<- paste0("post",1:3,".rds") # vector of filenames
#' PMC.appendFiles(filenm, "Tol2")
#'  
#' @seealso \code{\link{PMC.sampler}}
#' 
#' @export
#' 
PMC.appendFiles<- function(x, fileout=NULL) {
  # append posterior samples from saved rds files with the same tolerance 
  # x is a character vector of filenames  
  if(!is.character(x)) stop("x is not a character vector")
  if(length(x) < 2) stop("Need at least 2 files to append")
  n<- length(x)
  tmp<- readRDS(x[1])
  tol1<- c(tmp$ctol,tmp$stol)
  theta<- tmp$theta
  xr<- tmp$xr
  xs<- tmp$xs
  xspot<- tmp$xspot
  xscats<- tmp$xscats
  fpscats<- tmp$fpscats
  pop<- tmp$pop
  for(i in 2:n) {
    tmp<- readRDS(x[i])
    tol2<- c(tmp$ctol,tmp$stol)
    if(!all.equal(tol1,tol2)) stop(paste0("File ",x[i]," has a different tolerance"))
    theta<- rbind(theta, tmp$theta)
    xr<- rbind(xr, tmp$xr)
    xs<- rbind(xs, tmp$xs)
    xspot<- rbind(xspot, tmp$xspot)
    xscats<- rbind(xscats, tmp$xscats)
    fpscats<- rbind(fpscats,tmp$fpscats)
    pop<- c(pop, tmp$pop)
  }
  temp<- list(theta=theta,xr=xr,xs=xs,xspot=xspot,xscats=xscats,fpscats=fpscats,pop=pop)
  if(!is.null(fileout)) saveRDS(temp,file=fileout)
  temp
}
#-----------------------------------------------------------------------
#' Combine posterior samples of spatial maps of fox occupancy saved in rds files
#'
#' \code{Make.Raster} combines posterior samples from the saved file(s) of 
#' spatial maps of fox occupancy produced by \code{PMC.sampler} into 
#' a single raster file. Cell values are the occupancy probability calculated
#' by averaging the presence/absence values for each cell over all the posterior
#' samples.  The occupancy probabilities can be optionally smoothed using 
#' a Gaussian window with bandwidth \code{bw}.    
#' 
#' @param X A character vector containing the relative or full paths
#' and filenames for all saved files.  files should be binary 
#' files saved from \code{PMC.sampler} (by setting \code{save.post}) and 
#' should have filenames \code{pop#.rds} where \code{#} refers to the unique
#' file identifer.
#' @param rast raster habitat file that was used for simulating fox occupancy 
#' in \code{PMC.sampler}
#' @param years vector of years used for simulating fox occupancy
#' @param bw bandwidth (in meters) used for smoothing the final occupancy map
#' 
#' @return A list of raster files of the same length as \code{years}.
#' Each raster in the list has the same dimensions as \code{rast} 
#' containing the estimates of the fox occupancy probability for each cell.  
#' 
#' @examples
#' #filenm<- paste0("pop",1:3,".rds") # vector of filenames
#' #Make.Raster(filenm, habitat, 1995:2014, 3000)
#'  
#' @seealso \code{\link{PMC.sampler}}
#' 
#' @export
#' 
Make.Raster<- function(X, rast, years, bw=NULL) {
  if(!is.character(X)) stop("X must be a character vector of filenames")
  N<- length(X)
  nd<- dim(rast)[c(1,2)]
  tot.years<- length(years)
  if(!is.null(bw)) wmat<- focalWeight(rast, bw, "Gauss")
  
  occ<- vector("list",tot.years)
  count.rast<- rep(0,tot.years)
  for(j in 1:N) {
    cat("doing filename ",X[j]," file ",j, " of ",N,"\n")
    sim.list<- readRDS(X[j])
    n<- length(sim.list)
    
    for(i in 1:n) {    
      nyears<- length(sim.list[[i]]) 
      syear<- tot.years - nyears + 1  
      
      for(j in 1:nyears) {
        tmprast<- raster(rast)    
        r<- sim.list[[i]][[j]]
        r<- ifelse(r==2,1,0)
        dim(r)<- c(nd[1],nd[2])
        values(tmprast)<- r
        if(is.null(occ[[syear]])) occ[[syear]]<- tmprast
        else occ[[syear]]<- overlay(occ[[syear]],tmprast,fun="sum")
        count.rast[syear]<- count.rast[syear] + 1
        syear<- syear + 1
      }
    }
  }
  
  for(i in 1:tot.years) {
    if(is.null(occ[[i]])) {
      occ[[i]]<- raster(rast)
      values(occ[[i]])<- NA
      next
    }
    else{
      occ[[i]]<- occ[[i]]/count.rast[i]
      occ[[i]]<- focal(occ[[i]],w=wmat)
    }
  }
  # tidy up
  rm(sim.list)
  list(occ=occ,counts=count.rast)
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
#' \code{xr} - matrix of simulated road killed carcasses,
#' \code{xs} - matrix of simulated hunter killed carcasses,
#' \code{pop} - matrix of fox occupancy,
#' \code{xspot} - number of fox detections from spotlight surveys 
#'  
#' @seealso \code{\link{propose.theta}}
#' @export
#' 
ModelABC<- function(parm, Data) {
# Called by PMC.sampler
  nintro<- Data$nintro
  iend<- 2*nintro
  syear<- Data$syear
  eyear<- Data$eyear 
  pintro<- round(parm[1:nintro])
  yintro<- round(parm[(nintro+1):iend] + syear)
  nyears<- length(syear:eyear)
  if(sum(Data$nfpyrs) > 0) {
    pFP.parms<- parm[(iend+8):(iend+8+sum(Data$nfpyrs)-1)]
    pFP<- as.list(rep(0, nyears))
    pFP[Data$nfpyrs]<- pFP.parms
  } else pFP<- as.list(rep(0, nyears))
  Parms<- list(pintro=pintro,yintro=yintro,syear=syear,eyear=eyear,psurv=parm[iend+1],Ryear=parm[iend+2],proad=parm[iend+3],pshot=parm[iend+4],pbait=parm[iend+7],nintro=nintro)
  # Fox cellular automata C++ function from library(FoxSim) 
  mod<- foxsim(Data$habitat.mat, Data$road.mat, Data$ipoints, Data$kern.list, Data$baitlocs, Parms)
  xspot<- sapply(mod[[3]], function(x) spotlight.survey(x, Data$spotlocs, parm[iend+5], 0.2, 3))
  xscat<- mapply(scat.survey, mod[[3]], Data$scatsearch, pFP, MoreArgs=list(drate=parm[iend+6],parms=Data$scat.pars),SIMPLIFY=FALSE)
  list(years=syear:eyear,xr=mod[[1]],xs=mod[[2]],pop=mod[[3]],xspot=xspot,xscat=xscat)
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
#------------------------------------------------------------------------------
#' Release points 
#'
#' \code{release.points} returns cell coordinates of \code{habitat}
#' that are within \code{max.dist} of \code{intro.points}. Coordinates
#' are given as row and column number corresponding to \code{habitat}.
#' 
#' @param intro.points Geographic coordinates of introduction points.
#' points must lie with the bounding box of \code{habitat}
#' @param max.dist Maximum distance (m) from \code{intro.points} to include as
#' a potential release cell. Only cells coded as suitable habitat (1) are
#' included.
#' @param habitat Raster map of habitat with suitable habitat coded as 1
#' 
#' @return a list of size \code{nrow(intro.points)} with cell coordinates. 
#'
#' @examples
#' release.points(intro, fox.suit, 10000) # suitable cells within 10km from intro. 
#' @export
#' 
release.points<- function(intro.points, habitat, max.dist) {
  # returns cell coordinates (rn, cn) of release points within
  # max.dist of intro.points
  xy<- coordinates(intro.points)
  n<- nrow(xy)
  release.list<- list()
  for(i in 1:n) {
    tmp<- distanceFromPoints(habitat, xy[i,])
    cells<-  Which(tmp <= max.dist & habitat == 1, cells=T)
    release.list[[i]]<- rowColFromCell(habitat,cells)
  }
  release.list
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
  # pad observations with zeros
  nzeros<- length(x) - length(obs)
  c(rep(0,nzeros),obs)
}
#-----------------------------------------------------------------
#' Spatial Matching of cell locations
#'
#' \code{match.locations} performs spatial matching of cells in the vector
#'  \code{x} against locations given in the matrix \code{locs}. Used for 
#'  spatial matching of simulated observations with observed locations. 
#'  Wrapper for the C++ function \code{matchspatial}.
#'  
#' @param x Habitat matrix of cell locations that will be spatially matched 
#' against the locations in \code{locs}.
#' @param locs Matrix specifying the row and column number of
#' cell locations to match against.
#' @param ncell The spatial tolerance for specifying a match given
#' as the radial distance in cell units.
#' @param Val integer representing the cell value used for matching
#' 
#' @return a logical matrix with values 
#' set to \code{TRUE} for cell locations in \code{x} that are within
#'  \code{ncells} distance of locations given in \code{locs} 
#'  (i.e. within \code{ncells} of \code{x[locs[,1] + nr * locs[,2]]} 
#' 
#' @seealso \code{\link{ABC.reject}}, \code{\link{ABC.weighted}}
#' @export
match.locations<- function(x, locs, ncell, Val) {
  zz<- matchspatial(locs, x, ncell, Val)
  zz
}
#-----------------------------------------------------------------
#' Expected scats calculations
#'
#' \code{expected.scats} calculates the number of scats that are
#' expected to be available for detection within the home range of
#' a single fox. This is calculated as the equilibrium value 
#' between scat production and degradation.
#' Used by \code{scat.survey}.
#'  
#' @param pr scat production rate given as No. of scats/day 
#' @param dr scat degradation rate vector where \code{dr[1]}
#' is the (log) mean degradation rate per day and \code{dr[2]} is
#' the (log) standard deviation.
#' @param lf estimate of the proportion of scats that are 
#' deposited on linear features.
#' 
#' @return the (integer) number of scats 
#' 
#' @seealso \code{\link{scat.survey}}
#' @export
expected.scats<- function(pr, dr, lf) {
  b<- exp(rnorm(1,dr[1],dr[2])) # degradation rate
  a<- rpois(1,pr) # scat production per fox per day
  exp.scats<- -a*exp(-b*Inf)/b - (-a/b) # arbitrary long time 
  as.integer(exp.scats * lf) # this is expected scats per day
}
#-----------------------------------------------------------------
#' Scat survey monitoring observation model
#'
#' \code{scat.survey} performs detection of foxes from scat
#' surveys.  Simulated scat monitoring occurs on given cells
#' with the probability of detection given as a exponential
#' function of search distance with a given base detection rate. 
#' If the cell is occupied by a fox, then the number of scats
#' available for detection is given by \code{expected.scats}.
#'  
#' @param occ fox occupancy map matrix 
#' produced by \code{PMC.sampler}.
#' @param scatlocs matrix containing locations of cells subject to
#' scat searches with the first two columns containing the cell 
#' coordinates (row and column number) with the third column 
#' containing the transect distance (km).
#' @param drate base scat detection probability (per km).
#' @param parms parameters for calculating the expected number
#' of scats available for detection (see \code{expected.scats}).
#' 
#' @return matrix of the same dimension as \code{occ} with 1 
#' indicating cells with detected scats and 0 otherwise.
#' 
#' @seealso \code{\link{expected.scats}}, \code{\link{ModelABC}}
#' @export
scat.survey<- function(occ, scatlocs, fprate, drate, parms) {
  # Scat observation process
  # occ is occupancy status, scatlocs is the coordinates of the search effort
  # drate is the (log) per unit detection rate, parms are the parameters of
  # scat generation process
  atrisk<- which(occ[cbind(scatlocs[,1],scatlocs[,2])] == 2) #occupied sites 
  frisk<- which(occ[cbind(scatlocs[,1],scatlocs[,2])] != 2) #unoccupied sites
  n<- length(atrisk)
  nfp<- length(frisk)
  occ[,]<- 0  # set value of all cells to zero
  if(n > 0) { 
    exp.scats<- replicate(n, expected.scats(pr=parms$pr, dr=parms$dr, lf=parms$lf))
    slength<- scatlocs[atrisk,3]
    prob<- 1-exp(-exp(drate)*slength)
    dscats<- rbinom(n, exp.scats, prob)
    nonzero<- dscats > 0
    scoord<- matrix(scatlocs[atrisk[nonzero],c(1,2)],ncol=2)
    #occ[scoord[,1],scoord[,2]]<- dscats[nonzero]
    occ[scoord]<- 1
  }
  if(nfp > 0) {
    fpscats<- rbinom(nfp, 1, fprate)
    nonzero<- fpscats > 0
    fcoord<- matrix(scatlocs[frisk[nonzero],c(1,2)],ncol=2)
    #occ[scoord[,1],scoord[,2]]<- dscats[nonzero]
    occ[fcoord]<- 2
  }
  occ
}
#-----------------------------------------------------------------
#' Spotlight monitoring observation model
#'
#' \code{spotlight.survey} performs detection of foxes from spotlight
#' surveys.  Simulated spotlight monitoring occurs on given cells
#' with the probability of detection simulated using a random
#' Poisson process with a given base detection rate.  A distance to first
#' detection is then simulated and compared to the transect length 
#' in the cell to determine realised detections.
#'  
#' @param occ fox occupancy map matrix 
#' produced by \code{PMC.sampler}.
#' @param spotlocs matrix containing cell locations of spotlight
#' transects with the first two columns containing the cell coordinates
#' (row and column number) with the third column containing the distance of
#' the spotlight transect in the cell (km).
#' @param prob base spotlight detection probability (per km).
#' @param strip The spatlight stripwidth
#' @param cellsize the size of the cells in \code{occ} (km).
#' 
#' @return the total number of simulated detections
#' 
#' @seealso \code{\link{ModelABC}}
#' @export
spotlight.survey<- function(occ, spotlocs, prob, strip, cellsize){
  # spotlight survey detection probability on cells  
  drate<- -log(1-prob)  # detection rate per spotlight km
  atrisk<- which(occ[cbind(spotlocs[,1],spotlocs[,2])] == 2)  #occupied cells on spotlight transects
  pfox<- spotlocs[atrisk,3]*strip/cellsize^2 #probability of fox in transect
  trate<- spotlocs[atrisk,3] # length of spotlight transect in km 
  ttdet<- -log(runif(length(atrisk)))/drate  # distance to next event
  dtime<- ((ttdet <= trate) & (runif(length(atrisk)) <= pfox)) #fox in transect and detected 
  sum(dtime)
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
ABC.reject<- function(i, x0, ctol, stol, priors, Data) {
# ABC rejection algorithm sampling from priors 
  found<- FALSE
  ins<- Data$nintro
  ine<- 2*ins
  if(Data$Exact) extol<- 1 else extol<- ctol
    while(!found){
      thetac=sample.priors(priors)
      iord<- order(thetac[(ins+1):ine])
      thetac[1:ins]<- thetac[iord]
      thetac[(ins+1):ine]<- thetac[(iord+ins)]
      xc=ModelABC(thetac, Data)
      xr.sim<- sapply(xc$xr, sum)
      xs.sim<- sapply(xc$xs, sum)
      xspot.sim<- sum(xc$xspot)
      TPscat<- sapply(xc$xscat, function(x) sum(x==1)) # scats on occupied sites
      FPscat<- sapply(xc$xscat, function(x) sum(x==2)) # scats on unoccupied sites
      xscat.sim<- TPscat + FPscat #total is sum of true positive and false positive
      xr0<- pre.pad(xr.sim,x0$xr)
      xs0<- pre.pad(xs.sim,x0$xs)
      xscat0<- pre.pad(xscat.sim,x0$scat)
     
      if(distm(xr.sim,xr0) < ctol & distm(xs.sim,xs0) < extol & xspot.sim < extol 
         & distm(xscat.sim,xscat0) < stol){
        if(Data$MatchSpatial) {
          xr.sim<- lapply(xc$xr, match.locations,Data$carcass.road,Data$ncellC,1)
          xs.sim<- lapply(xc$xs, match.locations,Data$carcass.shot,Data$ncellC,1)
          xr.sim<- sapply(xr.sim, sum)
          xs.sim<- sapply(xs.sim, sum)
          TPscat<- mapply(match.locations, xc$xscat, Data$scats, 
                          MoreArgs = list(ncell=Data$ncellS,Val=1),SIMPLIFY=FALSE)  
          FPscat<- mapply(match.locations, xc$xscat, Data$scats, 
                          MoreArgs = list(ncell=Data$ncellS,Val=2),SIMPLIFY=FALSE) 
          xscat.sim<- sapply(TPscat, sum) + sapply(FPscat, sum)
          x1=do.call('c',mapply(function(x,y) as.numeric(x[y]), TPscat, Data$scats))
          x2=do.call('c',mapply(function(x,y) as.numeric(x[y]), FPscat, Data$scats))
          x1[which(x2==1)]<- 2
        } else x1<- rep(0, sum(xscat0)) # no spatial matching
        if(distm(xr.sim,xr0) < ctol & distm(xs.sim,xs0) < extol & distm(xscat.sim,xscat0) < stol){
          found<- TRUE
          theta<- thetac
          xr<- xr.sim
          xs<- xs.sim
          xspot<- xc$xspot
          xscats<- xscat.sim
          fpscats<- x1
          pop<- xc$pop
        }
      }
    }
  cat("completed particle ",i," for ctol ",ctol," and stol ",stol,"\n")
  list(theta,xr,xs,xspot,xscats,fpscats,pop)
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
ABC.weighted<- function(i, parms, x0, ctol, stol, VarCov, w, priors, Data) {
# ABC proposal algorithm sampling from weighted posterior 
  found<- FALSE
  ins<- Data$nintro # for ordering intoduction years and places
  ine<- 2*ins
  if(Data$Exact) extol<- 1 else extol<- ctol
  while(!found){
    thetao=pick.particle(parms, w)             
    thetac=propose.theta(thetao, VarCov, priors)  
    iord<- order(thetac[(ins+1):ine])
    thetac[1:ins]<- thetac[iord]
    thetac[(ins+1):ine]<- thetac[(iord+ins)]
    xc=ModelABC(thetac, Data) 
    xr.sim<- sapply(xc$xr, sum)
    xs.sim<- sapply(xc$xs, sum)
    xspot.sim<- sum(xc$xspot)
    TPscat<- sapply(xc$xscat, function(x) sum(x==1)) # scats on occupied sites
    FPscat<- sapply(xc$xscat, function(x) sum(x==2)) # scats on unoccupied sites
    xscat.sim<- TPscat + FPscat #total is sum of true positive and false positive
    xr0<- pre.pad(xr.sim,x0$xr)
    xs0<- pre.pad(xs.sim,x0$xs)
    xscat0<- pre.pad(xscat.sim,x0$scat)
    
    if(distm(xr.sim,xr0) < ctol & distm(xs.sim,xs0) < extol & xspot.sim < extol 
       & distm(xscat.sim,xscat0) < stol){
      if(Data$MatchSpatial) {
        xr.sim<- lapply(xc$xr, match.locations,Data$carcass.road,Data$ncellC,1)
        xs.sim<- lapply(xc$xs, match.locations,Data$carcass.shot,Data$ncellC,1)
        xr.sim<- sapply(xr.sim, sum)
        xs.sim<- sapply(xs.sim, sum)
        TPscat<- mapply(match.locations, xc$xscat, Data$scats, 
                        MoreArgs = list(ncell=Data$ncellS,Val=1),SIMPLIFY=FALSE)  
        FPscat<- mapply(match.locations, xc$xscat, Data$scats, 
                        MoreArgs = list(ncell=Data$ncellS,Val=2),SIMPLIFY=FALSE) 
        xscat.sim<- sapply(TPscat, sum) + sapply(FPscat, sum)
        
        x1=do.call('c',mapply(function(x,y) as.numeric(x[y]), TPscat, Data$scats))
        x2=do.call('c',mapply(function(x,y) as.numeric(x[y]), FPscat, Data$scats))
        x1[which(x2==1)]<- 2
      } else x1<- rep(0, sum(xscat0)) # no spatial matching
      if(distm(xr.sim,xr0) < ctol & distm(xs.sim,xs0) < extol & distm(xscat.sim,xscat0) < stol){
        found<- TRUE
        theta<- thetac
        xr<- xr.sim
        xs<- xs.sim
        xspot<- xc$xspot
        xscats<- xscat.sim
        fpscats<- x1
        pop<- xc$pop
      }
    }
  }
  cat("completed particle ",i," for ctol ",ctol," and stol ",stol,"\n")
  list(theta,xr,xs,xspot,xscats,fpscats,pop)
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
  imat<- qr.solve(VarCov)
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
#--------------------------------------------------------------------------------
#' Weighted variance calculation
#'
#' \code{var.wt} calculates the weighted variance for each column
#' of a matrix and is used for specifying a componenwise peturbation 
#' kernel.  Used in \code{PMC.sampler}.
#'  
#' @param x matrix 
#' @param w Vector of weights of length \code{nrow(x)}
#' #' 
#' @return a vector of weighted variances, one for each column of 
#' \code{x} 
#' 
#' @seealso \code{\link{PMC.sampler}}, \code{\link{propose.theta}},
#'   \code{\link{pick.particle}}
#' @export
var.wt<- function(x, w){
  # weighted variance function applied colwise
  wvar<- function(z, w) {
    n<- length(z)
    centre<- sum(z * w)
    xsqr<- (z - centre)^2
    return(sum(w * xsqr) * (n/(n-1)))
  }
  apply(x, 2, function(x) wvar(x, w=w))
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
#' implemented using a Cholski decomposition of \code{vcov}.
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
  # multivariate normal proposals using Cholski
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
    h <- (z - d*86400) %/% 60 %/% 60  # hours
    m <- (z - d*86400) %/% 60 %% 60 # mins 
    s <- (z - d*86400) %% 60 %% 60 # secs 
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