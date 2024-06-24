#graphs IM cariboo simulations###
library(tidyverse)
library(ggplot2)
library(jagsUI)
library(ggplot2)
library(grid)
library(gridExtra)
library(miceadds)

setwd("C:/LocalR/mesocarnivore_distribution_modelling/MDM_simulations")

#Chilcotin 3k-5k simulations ####

#1. Read files with simulations

file_list3k <- list.files("VM_simulations_chilcotin",  full.names="TRUE")
file_list5k <- list.files("VM_simulation_chilcotin_5k",  full.names="TRUE")

#file_list <- list.files("VM_simulations_cariboo",  full.names="TRUE")
file_list <- list.files("VM_simulations_omineca",  full.names="TRUE")

# this loops read the 1 object saved in each simulations, no more objects can be saved in that Rdata file or it will get one at random 

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

#list of all simulations for chilcotin

chilcotin.sim3k <- lapply(file_list3k, load_obj)
chilcotin.sim5k <- lapply(file_list5k, load_obj)
# Check if M is enough by looking at N.ever alive output
# Change p.o and p.s init parameters to 0.05
# simulate with fisher parameters and lower p.o and p.s

#1. Select the model 

#M= 500

# cariboo.list[[3]]  <- out.fisher_ICM_3
# cariboo.list[[4]]  <- out.fisher_ICM_4
# cariboo.list[[5]]  <- out.fisher_ICM_5
# #cariboo.list[[6]]  <- out.fisher_ICM_chilcotin_newN__1

# chilcotin.sim[[19]] <- out.fisher_ICM_chilcotin_newN_C_1

Abundance <- list()
N.df <- list()
plots <- list()

for(i in 1:length(chilcotin.sim)) {    
  N.df[[i]] <- as.data.frame(chilcotin.sim[[i]]$sims.list$N)
  colnames(N.df[[i]]) <- "N"
  #Abundance[[i]] <- N.df[[i]] %>% pivot_longer(cols = N)
  }

# to create only one graph with all simulations 

# N.df_all <- purrr::map_df(N.df, data.frame, .id = 'name')
#  N.df_all$group <- paste("3k", N.df_all$name)

## run above again for N.df_all2
# N.df_all2 <- purrr::map_df(N.df, data.frame, .id = 'name') ## comparison 
# N.df_all2$group <- paste("C", N.df_all2$name)

N.df_all3 <- purrr::map_df(N.df, data.frame, .id = 'name') ## comparison 
N.df_all3$group <- paste("O", N.df_all3$name)

N.df_all <- bind_rows(N.df_all, N.df_all2, N.df_all3)

chilcotin_3k_simulations <- ggplot(N.df_all, aes(x = group, y= N))+ geom_violin() + geom_hline(yintercept=500, linetype="dashed",  color = "red", size=1) + ylab("Pop. estimate")+ xlab("Simulations") + theme( axis.text.x=element_blank())

chilcotin_3k_simulations


## plot psi

psi <- list()
psi.df <- list()
plots_psi <- list()

for(i in 1:length(chilcotin.sim)) {    
  psi.df[[i]] <- as.data.frame(chilcotin.sim[[i]]$sims.list$psi)
  colnames(psi.df[[i]]) <- "psi"
#  psi[[i]] <- psi.df[[i]] %>% pivot_longer(cols = psi)
#  plots_psi[[i]] <- ggplot(psi[[i]], aes(x = reorder(name, sort(as.numeric(name))), y = value))+ geom_violin() + geom_hline(yintercept=0.3, linetype="dashed",  color = "red", size=1)+ xlab("") + ylab("Psi") #+ theme( axis.text.x=element_blank())
}

#grid.arrange(grobs= plots_psi, ncol=6)

psi_all <- purrr::map_df(psi.df, data.frame, .id = 'name')
psi_all_plot <- ggplot(psi_all, aes(x = reorder(name, sort(as.numeric(name))), y= psi))+ geom_violin() + geom_hline(yintercept=0.3, linetype="dashed",  color = "red", size=1) + ylab("Psi")+ xlab("Simulations")
# + theme( axis.text.x=element_blank())
psi_all_plot

## plot sigma 

sigma <- list()
sigma.df <- list()
plots_sigma<- list()

for(i in 1:length(chilcotin.sim)) {    
  sigma.df[[i]] <- as.data.frame(chilcotin.sim[[i]]$sims.list$sigma)
  colnames(sigma.df[[i]]) <- "sigma"
#  sigma[[i]] <- sigma.df[[i]] %>% pivot_longer(cols = sigma)
#  plots_sigma[[i]] <- ggplot(sigma[[i]], aes(x = reorder(name, sort(as.numeric(name))), y = value))+ geom_violin() + geom_hline(yintercept= 3, linetype="dashed",  color = "red", size=1)+ xlab("") + ylab("sigma") #+ theme( axis.text.x=element_blank())
}

#grid.arrange(grobs= plots_sigma, ncol=6)

sigma_all <- purrr::map_df(sigma.df, data.frame, .id = 'name')
sigma_all_plot <- ggplot(sigma_all, aes(x = reorder(name, sort(as.numeric(name))), y= sigma))+ geom_violin() + geom_hline(yintercept=3, linetype="dashed",  color = "red", size=1) + ylab("sigma")+ xlab("Simulations") #+ theme( axis.text.x=element_blank())
sigma_all_plot


### Columbian simulation plots ####

columbian.list <- list()
columbian.list[[1]] <- out.fisher_IMC_columbian_1
columbian.list[[2]]  <- out.fisher_IMC_columbianB_1


Abundance <- list()
N.df <- list()
plots <- list()

for(i in 1:length(columbian.list)) {    
  N.df[[i]] <- as.data.frame(columbian.list[[i]]$sims.list$N)
  colnames(N.df[[i]]) <- "N"
  Abundance[[i]] <- N.df[[i]] %>% pivot_longer(cols = N)
  plots[[i]] <- ggplot(Abundance[[i]], aes(x = name, y = value))+     geom_violin() + geom_hline(yintercept=500, linetype="dashed",  color = "red", size=1)+ xlab("") + ylab("Pop. estimate") + theme( axis.text.x=element_blank())
}

grid.arrange(grobs= plots, ncol=2)

#2. View model output table 
jags.View(model)

#3. rough plot of parameters 

trace.list3k <- list()


for(i in 1:length(chilcotin.sim3k)) {    

  trace.list3k[[i]] <-traceplot(chilcotin.sim3k[[i]] ,param = "N")} 

#lapply(trace.list3k,function(x){ggsave(file=paste(x,"pdf",sep="."),get(x))})

#c("p0.S", "N", "p0.O", "psi"))

trace.list5k <- list()

for(i in 1:length(chilcotin.sim5k)) {    
  
  trace.list5k[[i]] <-traceplot(chilcotin.sim5k[[i]] ,param= "N")
}




trace.list1 <- traceplot(out.fisher_ICM_chilcotin_newN__1 ,param= c("p0.S", "N", "p0.O", "psi"))
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
  geom_violin() + ylim(0, 1000) geom_hline(yintercept=500, linetype="dashed", 
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


plot1


