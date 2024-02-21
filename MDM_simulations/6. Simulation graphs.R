#graphs IM cariboo simulations###
library(tidyverse)
library(ggplot2)
library(jagsUI)
setwd("C:/LocalR/mesocarnivore_distribution_modelling/MDM_simulations/Bear parameters")

# Check if M is enough by looking at N.ever alive output
# Change p.o and p.s init parameters to 0.05
# simulate with fisher parameters and lower p.o and p.s


#1. Select the model 
model <- out.fisher_IM_1
model <- open.cariboo_test_IOM_1 # 5k iterations
model2 <- open.cariboo_test_IOM_1 ## 50k iterations
model3 <- open.chilcotin_test_IOM_chilcotin_1 # 50k iterations
print(model,3)
plot(model)
plot(model2)
plot(model3)

#2. View model output table 
jags.View(model)

#3. rough plot of parameters 

traceplot(model,param= c("p0.S", "N", "p0.O"))
traceplot(model2,param= c("p0.S", "N", "p0.O"))
traceplot(model,param= "deviance")
traceplot(model2,param= "deviance")

whiskerplot(model, parameters = c("N"))
whiskerplot(model3, parameters = c("N"))
whiskerplot(model, parameters = c("p0.S"))
whiskerplot(model2, parameters = c("p0.S"))
whiskerplot(model, parameters = c("p0.O"))
whiskerplot(model2, parameters = c("p0.O"))

#3. create a violin plot 

N.df <- as.data.frame(model3$sims.list$N)
colnames(N.df) <- c("N1", "N2", "N3", "N4")
Abundance <- N.df %>% pivot_longer(cols = c(N1, N2, N3, N4)) 

ggplot(Abundance, aes(x = name, y = value))+ geom_violin()
       
N.df <- as.data.frame(model$sims.list$N)
colnames(N.df) <- "N"
Abundance <- N.df %>% pivot_longer(cols = N) 

ggplot(Abundance, aes(x = name, y = value))+ geom_violin()



