require(nimble)

## Way more complicated new version:
fastSCRFunction <- nimbleFunction(
  setup = function(M, xlim, ylim, dx, dy, traps_camera, traps_hair, maxdist){
    # z_ <- numeric(M)
    # s_indx <- numeric(M)
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
  },
  run = function(){},
  methods = list(
    calcGridCell = function(s = double(1), z = double()){
      returnType(double())
      if(z == 0)
        return(0)

      cell <- c(floor((s[1]-xlim[1])/xrange*nx), floor((s[2]-ylim[1])/yrange*ny)) + 1
      s_index <- cell[1] + (cell[2]-1)*nx  ## This is the order of expand.grid.
      return(s_index)
    },
    calcDistSquared = function(s = double(1), device = character(0, "camera"), z = double(), sindex = double() ){
      returnType(double(1))
        
      if(device == "camera")
        dist2 <- numeric(value = Inf, length = J_camera)
      else
        dist2 <- numeric(value = Inf, length = J_hair)
      
      if(z == 0) return(dist2)
      
      if(device == "camera"){
        ind <- which(grd_camera[sindex, ] == 1)
        dist2[ind] <- (traps_camera[ind,1] - s[1])^2 + (traps_camera[ind,2] - s[2])^2
      }else{
        ind <- which(grd_hair[sindex, ] == 1)
        dist2[ind] <- (traps_hair[ind,1] - s[1])^2 + (traps_hair[ind,2] - s[2])^2    
      }
      return(dist2)
    },
    calcTrapDetection_HN = function(distsSquared = double(1), sigma = double(), z = double(), sindex = double(), device = character(0, default = "camera")){ #, detfn = character(0, default = "halfnormal")){
      returnType(double(1))
      
      if(device == "camera")
        ans <- numeric(value = 0, length = J_camera)
      else
        ans <- numeric(value = 0, length = J_hair)
      if(z == 0)
        return(ans)
   
      if(device == "camera"){
        ind <- which(grd_camera[sindex, ] == 1)
        ans[ind] <- exp(-0.5*distsSquared[ind]/sigma^2)
      }else{
        ind <- which(grd_hair[sindex, ] == 1)
        ans[ind] <- exp(-0.5*distsSquared[ind]/sigma^2)
      }
      return(ans)
    },
    dTrapPoisson_NoID = function(x = double(), p = double(1), lam0 = double(), T = double(), z = double(1), log = integer(0, default = 0)){ #, detfn = character(0, default = "halfnormal")){
      returnType(double())
      H <- lam0*sum(p[z == 1])
      ans <- x*log(T*H) - T*H
      if(log) 
        return(ans)
      else
        return(exp(ans))
    },
    dTrapBinom_NoID = function(x = double(), p = double(1), lam0 = double(), z = double(1), 
        nocc = double(), j = double(), sindices = double(1), log = integer(0, default = 0)){
      returnType(double())
      ind <- which(grd_camera[sindices[z==1],j] == 1)
      miss <- 1-lam0*p[ind]
      miss <- miss[miss > 0]
      logpmiss = sum(log(miss))
      ans <- (nocc-x)*logpmiss + x*log(1-exp(logpmiss))
      if(log) 
        return(ans)
      else
        return(exp(ans))
    },
    dTrapBinom = function(x = double(1), p = double(1), lam0 = double(), z = double(), 
        nocc = double(1), detfn = double(0, default = 1), log = integer(0, default = 0)){ #, detfn = character(0, default = "halfnormal")){
      returnType(double())
      if(z == 0){
        ans <- 0
      }else{
        ppos <- which(p > 0)  ## Setting detections far away as hard zeros.
        if(detfn == 1) ## "halfnormal"
          pobs <- p[ppos]*lam0
        if(detfn == 2) ## hazard halfnormal
          pobs <- 1-exp(-lam0*pobs)  ## This one is consistent with poisson halfnormal model.
          
        if(sum(x) == 0) {
          ans <- sum(nocc[ppos]*log(1-pobs))
        }else{
          ans <- sum( x[ppos]*log(pobs) + (nocc[ppos]-x[ppos])*log(1-pobs) )
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
  N <- sum(z[])     # Number of females
    
  for(i in 1:M) {
    z[i] ~ dbern(psi)

    s[i,1] ~ dunif(xlims[1], xlims[2])
    s[i,2] ~ dunif(ylims[1], ylims[2])

    sindex[i] <- fastSCR$calcGridCell(s=s[i,1:2], z=z[i])

    d2.s[i,1:J.s] <- fastSCR$calcDistSquared(s = s[i, 1:2], device = "hair", z = z[i], sindex=sindex[i])
    d2.o[i,1:J.o] <- fastSCR$calcDistSquared(s = s[i, 1:2], device = "camera", z = z[i], sindex=sindex[i])

    ## This has purposefully left out p0 for computation speed.
    p.S[i, 1:J.s] <- fastSCR$calcTrapDetection_HN(distsSquared=d2.s[i,1:J.s], sigma=sigma, sindex=sindex[i], z=z[i], device="hair")
    p.O[i, 1:J.o] <- fastSCR$calcTrapDetection_HN(distsSquared=d2.o[i,1:J.o], sigma=sigma, sindex=sindex[i], z=z[i], device="camera")
    
    ## p0 goes in here as it didn't impact that other part.
    y[i,1:J.s] ~ fastSCR$dTrapBinom(p=p.S[i,1:J.s], lam0=p0.S, nocc=nocc.s[1:J.s], z=z[i], detfn=1)
  }#m
  
  for(j in 1:J.o) { #for every trap
      O[j] ~ fastSCR$dTrapBinom_NoID(p=p.O[1:M, j], lam0=p0.O, nocc=nocc.o[j], j=j, sindices=sindex[1:M], z=z[1:M]) 
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




## ARCHIVE FOR LATER:
# calcDistSquared <- nimbleFunction(
  # run = function(s = double(1), traps = double(2), z = double()){
    # returnType(double(1))
    # if(z == 0)
      # return(numeric(value = 0, length = dim(traps)[1]))

      # dist2 <- (traps[,1] - s[1])^2 + (traps[,2] - s[2])^2
      # return(dist2)
  # }
# )

# calcTrapDetection_HN = nimbleFunction(
  # run = function(distsSquared = double(1), sigma = double(), z = double()){ #, detfn = character(0, default = "halfnormal")){
    # returnType(double(1))
    # if(z == 0)
      # return(numeric(value = 0, length = dim(distsSquared)[1]))
    ## if(detfn == "halfnormal")
    # return(exp(-0.5*distsSquared/sigma^2))
  # }
# )

# dTrapPoisson_NoID = nimbleFunction(
  # run = function(x = double(), p = double(1), lam0 = double(), T = double(), z = double(1), log = integer(0, default = 0)){ #, detfn = character(0, default = "halfnormal")){
    # returnType(double())
    # H <- lam0*sum(p[z == 1])
    # ans <- x*log(T*H) - T*H
    # if(log) 
      # return(ans)
    # else
      # return(exp(ans))
  # }
# )

# dTrapBinom_NoID = nimbleFunction(
  # run = function(x = double(), p = double(1), lam0 = double(), z = double(1), 
      # nocc = double(), log = integer(0, default = 0)){
    # returnType(double())
    # logpmiss = sum(log(1-lam0*p[z==1]))
    # ans <- (nocc-x)*logpmiss + x*log(1-exp(logpmiss))
    # if(log) 
      # return(ans)
    # else
      # return(exp(ans))
  # }
# )

# dTrapBinom = nimbleFunction(
  # run = function(x = double(1), p = double(1), lam0 = double(), z = double(), 
      # nocc = double(1), detfn = double(0, default = 1), log = integer(0, default = 0)){ #, detfn = character(0, default = "halfnormal")){
    # returnType(double())
    # if(z == 0){
      # ans <- 0
    # }else{
      # ppos <- which(p > 0)
      # if(detfn == 1) ## "halfnormal"
      # pobs <- p[ppos]*lam0
      # if(detfn == 2) ## hazard halfnormal
        # pobs <- 1-exp(-lam0*pobs)  ## This one is consistent with poisson halfnormal model.
        
      # if(sum(x) == 0) {
        # ans <- sum(nocc[ppos]*log(1-pobs))
      # }else{
        # ans <- sum( x[ppos]*log(pobs) + (nocc[ppos]-x[ppos])*log(1-pobs) )
      # }
    # }
    # if(log) 
      # return(ans)
    # else
      # return(exp(ans))
  # }
# )
