## ------------------------------------------------------------------------
library(FoxSim, quietly=T)
data(habitatvec)
data(roadvec)


## ------------------------------------------------------------------------
nd<- attributes(habitatvec)$dims
nd

## ------------------------------------------------------------------------
data(carcass_obs)
carcass_obs

## ------------------------------------------------------------------------
# local neighbourhood dispersal kernel
kern1<- kernel2D(eps=3, kfun="gamma", shape=1, scale=3, max.disp=30)

# Vagrant (searching) dispersal kernel. equal dispersal probability  
# in all directions out to 30 km.
kern2<- kernel2D(eps=3, kfun="circle", sigma=30, max.disp=30)

kern<- list(kern1, kern2)

kdim<- dim(kern1)[1] # save dimensions of kernel


## ----,fig.width=10,fig.height=8------------------------------------------
par(mfrow=c(1,2),pty="s")
x<- y<- seq(-30, 30, 3)
image(x, y, kern1, col=terrain.colors(50),main="local dispersal kernel")
image(x, y, kern2, col=terrain.colors(50),main="vagrant dispersal kernel")


## ------------------------------------------------------------------------
data(incpoints)
incpoints

## ------------------------------------------------------------------------
Parms<- list(syear=1998,eyear=2013,psurv=0.7,Ryear=1.1,proad=0.2,pshot=0.1)

# Now do a single run of the model

foxsim(nd[1], nd[2], kdim, habitatvec, roadvec, incpoints, kern, Parms)
       

## ------------------------------------------------------------------------
MyData<- list(nr=nd[1],nc=nd[2],kdim=kdim, hab.vec=habitatvec, road.vec=roadvec,  
              ipoints=incpoints, kern.list=kern, eyear=2014)


## ------------------------------------------------------------------------
priors<- list(list("uniform",-3.49,3.49),
              list("beta",1,1),
              list("beta",17,83),
              list("beta",17,83),
              list("uniform",0.1,2))

## ------------------------------------------------------------------------
Tol<- c(5,4)

fit<- PMC.sampler(N=10, carcass_obs, Tol, priors, MyData, parallel=F)
theta<- fit$Tol4$theta  # parameter estimates for final tolerance
theta[,1]<- round(theta[1] + 1998) # transform syear to integer year of intro.
theta


