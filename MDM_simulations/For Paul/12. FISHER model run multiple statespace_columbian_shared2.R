##########################
### fit marginal closed model for different state spaces in the same model, Columbian
##########################

#Objects in environment

#Area.list <-  list with the area of each state space in sqkm
#bin.o_num <-  array of fisher records in each camera, each row is a state space (1-7) and each column is the number of fisher detections in that camera. Max number of cameras 188.

#Y.s.binom_num <- array with number of individual fisher captured at each hair trap. each row is a unique fisher and each column is a hair trap. 

#J<-  vector with number of  sampling sites per state (hair snaggers and cameras)
#Traps_list_num <- 3D array with hair trap coordinates per state space. dim (346, 2, 7)
#cam_list_num <- 3D array with camera coordinates per state space. dim (188, 2, 7)
#xlims <- limits of x coordinates for each state space, each row is a state space. 
#ylims <- limits of y coordinates for each state space, each row is a state space.

#nAnimals <- maximum number of identified fishers per state space

#animal.state <- state space that each individual fisher belonged to. 

#nstates <- number of states in the model, 7. 


#setwd("C:/LocalR/data_collation/Chilcotin new model")

load("full_data.RData")
start.time <- Sys.time()

library(rjags)
library(jagsUI)
library(abind)

nstates <- length(Area.list)
Area <- unlist(Area.list) # in km2
M<-c(500, 500, 300, 500, 500, 300, 1500)

init_simple <- function() {
  
  zi <- list()
  sii <- list()  
  si <-  list()

  for(state in 1:7){
    for (i in 1:M[state]) {
      
      zi[[state]] <- matrix(0L, M[state], 1)
      zi[[state]][1:M[state]] <- 1  #  zi[1:M)] <- 1 give 1's to all individuals
      si[[state]] <- cbind(runif(M[state], xlims[state, 1], xlims[state, 2]),
                       runif(M[state], ylims[state, 1], ylims[state, 2]))
      
      max_rows_s<- max(sapply(si, nrow))
      max_cols_s <- max(sapply(si, ncol))
      # Pad each array to the maximum dimensions
      padded_arrays_s_list <- lapply(si, function(mat) {
        padded_mat<- matrix(0, nrow = max_rows_s, ncol = max_cols_s)
        padded_mat[1:nrow(mat), 1:ncol(mat)] <- mat
        return(padded_mat)
      })
      
    max_rows_z<- max(sapply(zi, nrow))
    max_cols_z<- max(sapply(zi, ncol))

      # Pad each array to the maximum dimensions
      padded_arrays_z_list <- lapply(zi, function(mat) {
        padded_mat <- matrix(0, nrow = max_rows_z, ncol = max_cols_z)
        padded_mat[1:nrow(mat), 1:ncol(mat)] <- mat
        return(padded_mat)
      })
      
    # Combine the padded arrays into a 3D matrix
    s_num <-  abind(padded_arrays_s_list, along = 3)
    z_num <-  abind(padded_arrays_z_list, along = 2)
    
  }
 
  } 
  list(z= z_num, s= s_num, p0.S=0.1, p0.O=0.1,
       sigma=3)
}

pars <- c("N","psi","p0.S","p0.O","sigma","N", "Density")


  name.i <- "Multiple_state_space_columbian"
  out.i <- paste("out.",  "_", sep = "")

  data <- list(
    y= Y.s.binom_num2,
    O = bin.o_num, 
    M = M,
    X= X, 
    J= J,
    trapType= trapType,
    xlims = xlims,
    ylims = ylims,
    nAnimals= sum(nAnimals),
    NOCC= NOCC_num,# number of occasion per cameras
    animal.state = animal.state,
    nstates= nstates,
    Area = Area
  )
  out <-
    jags(
      "margSingle_IM_fisher_Example_States - Copy_cindy.JAG",
      data = data,
      inits = init_simple,
      parallel = FALSE, n.cores= 1,
      n.chains = 3,
      n.burnin = 10,
      n.adapt = 80,
      n.iter = 100,
      parameters.to.save = pars
    )
  assign(out.i, out)
  save(list = out.i, file = paste(out.i, "test_fisher_multi_state_columbian.Rdata", sep = ""))
  rm(name.i, out.i, out)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


