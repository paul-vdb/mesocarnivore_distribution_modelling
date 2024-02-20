#graphs IM cariboo simulations###
library(tidyverse)
library(ggplot2)
library(jagsUI)
setwd("C:/LocalR/mesocarnivore_distribution_modelling/MDM_simulations")

# Check if M is enough by looking at N.ever alive output
# Change p.o and p.s init parameters to 0.05
# simulate with fisher parameters and lower p.o and p.s


#1. Select the model 
model <- open.cariboo_test_IOM_1
print(model,3)
plot(model)
#2. View model output table 
jags.View(model)

#3. rough plot of parameters 

traceplot(model,param= c("p0.S", "N", "p0.O"))
traceplot(model,param= "deviance")
traceplot(model,param= "deviance")

whiskerplot(model, parameters = c("N"))
whiskerplot(model, parameters = c("p0.S"))
whiskerplot(model, parameters = c("p0.O"))

#3. create a violin plot 

N.df <- as.data.frame(model$sims.list$N)
colnames(N.df) <- c("N1", "N2", "N3", "N4")
test <- N.df %>% pivot_longer(cols = c(N1, N2, N3, N4)) 

ggplot(test, aes(x = name, y = value))+ geom_violin()
       



