#graphs IM cariboo simulations###
library(tidyverse)
library(ggplot2)
library(jagsUI)
library(ggplot2)
library(grid)
library(gridExtra)
setwd("C:/LocalR/mesocarnivore_distribution_modelling/MDM_simulations")

# Check if M is enough by looking at N.ever alive output
# Change p.o and p.s init parameters to 0.05
# simulate with fisher parameters and lower p.o and p.s

#1. Select the model 

 #M= 500
#priors not set.
cariboo.list <- list()
cariboo.list[[1]] <- out.fisher_ICM_1
cariboo.list[[2]]  <- out.fisher_ICM_2
cariboo.list[[3]]  <- out.fisher_ICM_3
cariboo.list[[4]]  <- out.fisher_ICM_4
cariboo.list[[5]]  <- out.fisher_ICM_5
#cariboo.list[[6]]  <- out.fisher_ICM_chilcotin_newN__1


Abundance <- list()
N.df <- list()

for(i in 1:length(cariboo.list)) {    
  N.df[[i]] <- as.data.frame(cariboo.list[[i]]$sims.list$N)
  colnames(N.df[[i]]) <- "N"
  Abundance[[i]] <- N.df[[i]] %>% pivot_longer(cols = N)
  
}
plot1<- ggplot(Abundance[[1]], aes(x = name, y = value))+ 
  geom_violin()
plot2<- ggplot(Abundance[[2]], aes(x = name, y = value))+ 
  geom_violin()
plot3<- ggplot(Abundance[[3]], aes(x = name, y = value))+ 
  geom_violin()
plot4<- ggplot(Abundance[[4]], aes(x = name, y = value))+ 
  geom_violin()
plot5<- ggplot(Abundance[[5]], aes(x = name, y = value))+ 
  geom_violin()
# plot6<- ggplot(Abundance[[6]], aes(x = name, y = value))+ 
#   geom_violin()

plot <- list(plot1, plot2, plot3, plot4, plot5)#, plot6)

n <- length(plot)
grid.arrange(grobs= plot, ncol=2)


m2 <- out.fisher_ICM_chilcotin_newN__1 # 5k iterations
model2 <- open.cariboo_test_IOM_1 ## 50k iterations
model3 <- open.chilcotin_test_IOM_chilcotin_1 # 50k iterations

print(cariboo_m1,3)
plot(model)
plot(model2)
plot(model3)

#2. View model output table 
jags.View(model)

#3. rough plot of parameters 

trace.list <- list()
trace.list[[1]] <- traceplot(cariboo.list[[1]] ,param= c("p0.S", "N", "p0.O", "psi"))
trace.list[[2]] <- traceplot(cariboo.list[[2]] ,param= c("p0.S", "N", "p0.O", "psi"))
trace.list[[3]] <- traceplot(cariboo.list[[3]] ,param= c("p0.S", "N", "p0.O", "psi"))
trace.list[[4]] <- traceplot(cariboo.list[[4]] ,param= c("p0.S", "N", "p0.O", "psi"))
trace.list[[5]] <- traceplot(cariboo.list[[5]] ,param= c("p0.S", "N", "p0.O", "psi"))

grid.arrange(grobs= trace.list, ncol=2)

traceplot(m2,param= c("p0.S", "N", "p0.O"))
traceplot(model,param= "deviance")
traceplot(model2,param= "deviance")

whiskerplot(m2, parameters = c("N"))
whiskerplot(model3, parameters = c("N"))
whiskerplot(model, parameters = c("p0.S"))
whiskerplot(model2, parameters = c("p0.S"))
whiskerplot(model, parameters = c("p0.O"))
whiskerplot(model2, parameters = c("p0.O"))

#3. create a violin plot 

N.df <- as.data.frame(m2$sims.list$N)
colnames(N.df) <- c("N1")
Abundance <- N.df %>% pivot_longer(cols = c(N1)) 

ggplot(Abundance, aes(x = name, y = value))+ geom_violin()

# 4. Start here with new models
m2 <- out.fisher_ICM_chilcotin_newN__1

model.list <- list()
model.list[[1]] <- out.fisher_ICM_chilcotin_newN__1
model.list[[2]] <- out.fisher_ICM_chilcotin_newN__2
model.list[[3]] <- m2
Abundance <- list()
N.df <- list()

for(i in 1:length(model.list)) {    
N.df[[i]] <- as.data.frame(model.list[[i]]$sims.list$N)
colnames(N.df[[i]]) <- "N"
Abundance[[i]] <- N.df[[i]] %>% pivot_longer(cols = N)

}
plot1<- ggplot(Abundance[[1]], aes(x = name, y = value))+ 
  geom_violin() + geom_hline(yintercept=500, linetype="dashed", 
                             color = "red", size=2)
plot2<- ggplot(Abundance[[2]], aes(x = name, y = value))+ 
  geom_violin() + geom_hline(yintercept=500, linetype="dashed", 
                               color = "red", size=2)
plot3<- ggplot(Abundance[[3]], aes(x = name, y = value))+ 
  geom_violin()+ geom_hline(yintercept=500, linetype="dashed", 
                            color = "red", size=2)
plot4<- ggplot(Abundance[[4]], aes(x = name, y = value))+ 
  geom_violin()

plot <- list(plot1, plot2, plot3)#, plot4)

n <- length(plot)
grid.arrange(grobs= plot, ncol=3)
