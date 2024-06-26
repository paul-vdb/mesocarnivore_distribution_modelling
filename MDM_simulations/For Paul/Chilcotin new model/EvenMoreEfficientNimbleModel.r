require(nimble)

## Way more complicated new version:
fastSCRFunction <- nimbleFunction(
  setup = function(M, xlim, ylim, dx, dy, traps_camera, traps_hair, maxdist){
    # z_ <- numeric(M)
    # sigma_ <- 1
    # lam0_id <- 0.1 
    # lam0_noid <- 0.1
    
    xrange <- diff(xlim)
    yrange <- diff(ylim)
    dx <- xrange / ceiling(xrange/dx)
    dy <- yrange / ceiling(yrange/dy)

    maxd2 <- maxdist^2
    xgrd <- seq(xlim[1], xlim[2], by = dx)
    nx <- length(xgrd)
    ygrd <- seq(ylim[1], ylim[2], by = dy)
    ny <- length(ygrd)

    grd <- as.matrix(expand.grid(xgrd, ygrd))
    grd.indx <- cbind(floor((grd[,1]-xlim[1])/xrange*nx) + 1 ,floor((grd[,2]-ylim[1])/yrange*ny) + 1)
    
    grd_hair <- t(apply(grd, 1, FUN = function(x){(x[1] - traps_hair[,1])^2 + (x[2] - traps_hair[,2])^2}))
    grd_hair <- (grd_hair <= maxd2)*1

    grd_camera <- t(apply(grd, 1, FUN = function(x){(x[1] - traps_camera[,1])^2 + (x[2] - traps_camera[,2])^2}))
    grd_camera <- (grd_camera <= maxd2)*1

    J_camera <- nrow(traps_camera)
    J_hair <- nrow(traps_hair)

    d2_hair <- matrix(0, nrow = M, ncol = J_hair)
    d2_camera <- matrix(0, nrow = M, ncol = J_camera)

    p_hair <- matrix(0, nrow = M, ncol = J_hair)
    p_camera <- matrix(0, nrow = M, ncol = J_camera)
  },
  run = function(){},
  methods = list(
    calcGridCell = function(s = double(1), z = double(), id = double()){
      returnType(double())
      if(z == 0)
        return(0)

      cell <- c(floor((s[1]-xlim[1])/xrange*nx), floor((s[2]-ylim[1])/yrange*ny)) + 1
      sindex <- cell[1] + (cell[2]-1)*nx  ## This is the order of expand.grid.
      if(sindex < 1 | sindex > (nx*ny)) stop("Sindex error")
      
      if(z != 0){
        indcam <- which(grd_camera[sindex, ] == 1)
        if( dim(indcam)[1] > 0 ){
          d2_camera[id,indcam] <<- (traps_camera[indcam,1] - s[1])^2 + (traps_camera[indcam,2] - s[2])^2
        }        
        indhair <- which(grd_hair[sindex, ] == 1)
        if( dim(indhair)[1] > 0 ){
          d2_hair[id,indhair] <<- (traps_hair[indhair,1] - s[1])^2 + (traps_hair[indhair,2] - s[2])^2
        }
      }
      return(sindex)
    },
    calcProbDet = function(z = double(), sindex = double(), sigma = double(), id = double()){
      returnType(double(0))

      pcam <- numeric(value = 0, length = J_camera)
      phair <- numeric(value = 0, length = J_hair)
      
      if(z == 1){
        indcam <- which(grd_camera[sindex, ] == 1)
        if( dim(indcam)[1] > 0 ){
          pcam[indcam] <- exp(-0.5*d2_camera[id, indcam]/sigma^2)
        }        
        indhair <- which(grd_hair[sindex, ] == 1)
        if( dim(indhair)[1] > 0 ){
          phair[indhair] <- exp(-0.5*d2_hair[id, indhair]/sigma^2)
        }
      }
      p_hair[id,] <<- phair
      p_camera[id,] <<- pcam
      return(0)
    },
    dTrapBinom_NoID = function(x = double(), lam0 = double(), z = double(1), 
        nocc = double(), j = double(), nothing = double(1), log = integer(0, default = 0)){
      returnType(double())
      ind <- which(p_camera[,j] > 0)
      if(dim(ind)[1] != 0){
        miss <- 1-lam0*p_camera[ind, j]
        logpmiss = sum(log(miss))
        ans <- (nocc-x)*logpmiss + x*log(1-exp(logpmiss))
      }else{
        ans <- 0
      }
      if(log) 
        return(ans)
      else
        return(exp(ans))
    },
    dTrapBinom = function(x = double(1), lam0 = double(), z = double(), id = double(),
        nocc = double(1), detfn = double(0, default = 1), nothing = double(), log = integer(0, default = 0)){ #, detfn = character(0, default = "halfnormal")){
      returnType(double())
      if(z == 0){
        ans <- 0
      }else{
        ppos <- which(p_hair[id,] > 0)  ## Setting detections far away as hard zeros.
        if(dim(ppos)[1] != 0){
          if(detfn == 1) ## "halfnormal"
            pobs <- p_hair[id, ppos]*lam0
          if(detfn == 2) ## hazard halfnormal
            pobs <- 1-exp(-lam0*pobs)  ## This one is consistent with poisson halfnormal model.  
          if(sum(x) == 0) {
            ans <- sum(nocc[ppos]*log(1-pobs))
          }else{
            ans <- sum( x[ppos]*log(pobs) + (nocc[ppos]-x[ppos])*log(1-pobs) )
          }
        }else{
          ans <- 0
        }
      }
      if(log) 
        return(ans)
      else
        return(exp(ans))
    }
  )  
)

## Combined SCR and OCC model built for speed.
scr_occ_model_maxdist <- nimbleCode({
  psi ~ dbeta(1, 1)    # M*psi = E[N(1)]
  p0.S ~ dbeta(1, 1)
  p0.O ~ dbeta(1, 1)
  
  sigma ~ dgamma(6, 4) # followed burgar et. al. 2018 for calculations, 5 to 300km2
  N <- sum(z[1:M])     # Number of females
    
  for(i in 1:M) {
    z[i] ~ dbern(psi)

    s[i,1] ~ dunif(xlims[1], xlims[2])
    s[i,2] ~ dunif(ylims[1], ylims[2])

    sindex[i] <- fastSCR$calcGridCell(s=s[i,1:2], z=z[i], id=i)
    p[i] <- fastSCR$calcProbDet(z=z[i], sindex=sindex[i], sigma=sigma, id=i)
    
    ## p0 goes in here as it didn't impact that other part.
    y[i,1:J.s] ~ fastSCR$dTrapBinom(lam0=p0.S, nocc=nocc.s[1:J.s], z=z[i], detfn=1, id=i, nothing = p[i])
  }#m
  
  for(j in 1:J.o) { #for every trap
      O[j] ~ fastSCR$dTrapBinom_NoID(lam0=p0.O, nocc=nocc.o[j], j=j, z=z[1:M], nothing = p[1:M]) 
    }#j.o
  }#m
)


## Combined SCR and OCC model built for speed.
scr_occ_model <- nimbleCode({
  psi ~ dbeta(1, 1)    # M*psi = E[N(1)]
  p0.S ~ dbeta(1, 1)
  p0.O ~ dbeta(1, 1)
  
  sigma ~ dgamma(6, 4) # followed burgar et. al. 2018 for calculations, 5 to 300km2
  N <- sum(z[1:M])     # Number of females
    
  for(i in 1:M) {
    z[i] ~ dbern(psi)

    s[i,1] ~ dunif(xlims[1], xlims[2])
    s[i,2] ~ dunif(ylims[1], ylims[2])

    d2.s[i,1:J.s] <- calcDistSquared(s = s[i, 1:2], traps = X.s[1:J.s, 1:2], z = z[i])
    d2.o[i,1:J.o] <- calcDistSquared(s = s[i, 1:2], traps = X.o[1:J.o, 1:2], z = z[i])

    ## This has purposefully left out p0 for computation speed.
    p.S[i, 1:J.s] <- calcTrapDetection_HN(distsSquared = d2.s[i,1:J.s], sigma = sigma, z = z[i])
    p.O[i, 1:J.o] <- calcTrapDetection_HN(distsSquared = d2.o[i,1:J.o], sigma = sigma, z = z[i])
    
    ## p0 goes in here as it didn't impact that other part.
    y[i,1:J.s] ~ dTrapBinom(p = p.S[i,1:J.s], lam0 = p0.S, nocc = nocc.s[1:J.s], z = z[i], detfn = 1)
  }#m
  
  for(j in 1:J.o) { #for every trap
      O[j] ~ dTrapBinom_NoID(p = p.O[1:M, j], lam0 = p0.O, nocc = nocc.o[j], z = z[1:M]) 
    }#j.o
  }#m
)

## This model input has counts instead of binomials
scr_count_model <- nimbleCode({
  psi ~ dbeta(1, 1)    # M*psi = E[N(1)]
  lam ~ unif(0, 100)  ## removed p.O + p.S
  
  sigma ~ dgamma(6, 4) # followed burgar et. al. 2018 for calculations, 5 to 300km2
  N <- sum(z[])     # Number of females
    
  for(i in 1:M) {
    z[i] ~ dbern(psi)

    s[i,1] ~ dunif(xlims[1], xlims[2])
    s[i,2] ~ dunif(ylims[1], ylims[2])

    d2.s[i,1:J.s] <- calcDistSquared(s = s[i, 1:2], traps = X.s[1:J.s, 1:2], z = z[i])
    d2.o[i,1:J.o] <- calcDistSquared(s = s[i, 1:2], traps = X.o[1:J.o, 1:2], z = z[i])

    ## This has purposefully left out p0 for computation speed.
    p.S[i, 1:J.s] <- calcTrapDetection_HN(distsSquared = d2.s[i,1:J.s], sigma = sigma, z = z[i])
    p.O[i, 1:J.o] <- calcTrapDetection_HN(distsSquared = d2.o[i,1:J.o], sigma = sigma, z = z[i])
    
    ## p0 goes in here as it didn't impact that other part.
    y[i,1:J.s] ~ dTrapBinom(p = p.S[i,1:J.s], lam0 = lam, nocc = nocc.s[1:J.s], z= z[i], detfn = 2)
  }#m
  
  ## This is Chandler and Royle 2013
  for(j in 1:J.o) { #for every trap
      O[j] ~ dTrapPoisson_NoID(p = p.O[1:M, j], lam0 = lam, T = nocc.o[j], z = z[1:M]) 
    }#j.o
  }#m
)



