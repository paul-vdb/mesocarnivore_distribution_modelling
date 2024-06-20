
## . clean data for DNA analysis 
library(terra)
library(ggplot2)
library(sf)
library(readr)
library(farver)
library(ggplot2)
library(rmapshaper)
library(stars)
library(tidyterra)
library(dplyr)
library(rmapshaper)
# X and Y data already obtained on Simulations. 

# Install packages that are not yet installed and load them using Pacman
pacman::p_load(ggplot2, tidyr, dplyr, readxl, knitr, lubridate, tidyverse, lme4, lattice, mgcv, pacman, sf)

#1. Read compiled DNA and Camera files ####

setwd("I:/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Mesocarnivore_Monitoring_Program/2.Data/Mesocarnivores DB/1. Master Data")

DNA_data <- read_csv("DNA_data_MDB_03-06.csv",col_types = cols(Start_date = col_date(format = "%Y-%m-%d"), End_date = col_date(format = "%Y-%m-%d"))) %>% 
  drop_na(c(Latitude_DD, Longitude_DD)) %>%  
  sf::st_as_sf(coords= c('Longitude_DD', 'Latitude_DD'), crs= 4326, remove = FALSE ) %>%
  st_transform( crs=3005) %>% 
  dplyr::mutate(UTM_X = sf::st_coordinates(.)[,1], UTM_Y = sf::st_coordinates(.)[,2])


DNA_data$DATA_TYPE <- "DNA"

# 2. Read Area of interest- Fisher Columbian Population #### 

setwd("C:/LocalR")
meso_grid <- st_read("BC_meso_grid.shp")
#grid_columbian <- st_read("BC_meso_grid_columbian.shp")
#  columbian_area <- st_union(grid_columbian, by_feature = FALSE) %>% ms_simplify(., keep = 0.01, keep_shapes = FALSE)

# 3. add grid information ####

DNA_grid <- st_intersection(DNA_data, meso_grid)

#4. Filter by population, omineca ####

setwd("C:/LocalR")
#setwd("C:/Users/cindyhurtado/OneDrive - Government of BC/VM")

subpopulations <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Subpopulations")

subpop <- ms_simplify(subpopulations, keep = 0.001,
                      keep_shapes = FALSE)

chilcotin <- subpop |> dplyr::filter(Subpop == "Chilcotin")

chilcotin <- st_buffer( chilcotin, 1000)
DNA_grid2 <- st_intersection( DNA_grid, chilcotin)


DNA_grid_columbian <- DNA_grid2 %>% subset(Project_name== "Enterprise" | Project_name == "Chilcotin" | Project_name == "Bonaparte" | Project_name == "williston" | Project_name == "Skeena" | Project_name == "Omineca-Lower Chilako" | Project_name == "JPRF") # studies from 2018 onwards

#cameras_grid2 # contains all cameras in columbian area, no need to filter any more. 

# DNA_columbian_3k<-  DNA_grid_columbian %>% group_by(Project_name, Station_ID, .keep_all = TRUE)
DNA_columbian_3k <- DNA_grid_columbian

# The hope is to make this into a list of identical tables.

head(DNA_columbian_3k)

# remove animals IDS that are not from Fisher, not filter as we need all stations 

DNA_columbian_3k <- DNA_columbian_3k %>% 
  mutate(
    Animal_ID = case_when(
      Species == "Pekania pennanti" ~ Animal_ID,
      Species != ~ "Pekania pennanti" ~ NA
    )
  )

library(data.table)
DNA_columbian_3k <- DNA_columbian_3k %>% 
  dplyr::group_by(Project_name, Station_ID) %>% 
  arrange(Project_name, Station_ID, End_date) %>% 
  dplyr::mutate(occasion = rleid(End_date)) %>%
  drop_na(occasion)

#3. Extract coordinates in Albers 
coord.scale <- 1000

X.s <- DNA_columbian_3k %>% select(Station_ID, UTM_X, UTM_Y) %>%  st_drop_geometry()

X.s <-X.s[!duplicated(X.s), ]  %>%  ungroup() %>% select(UTM_X,UTM_Y)
X.s <- X.s/coord.scale
X.s <- as.matrix(X.s, dim = c(dim(X.s)[1], 2))
X.s <- unname(X.s) # removed this to see if it helps with error

#Create Y.s array

# This script is a little exploration of what Cindy wants for the 'DNA_columbian_3k' table.

# We need a row for unique animal IDs, columns will be unique station IDs filled by a YES (1) or NO (0)
# for detections. Note that this is for ALL occasions and thus has some duplicated data in it.
# It's a preliminary check of the code!
# DNA_columbian_3k %>% 
#   dplyr::ungroup() %>% 
#   sf::st_drop_geometry() %>% 
#   dplyr::filter(occasion == 1, Animal_ID == '1011') %>% 
#   dplyr::filter(!duplicated(Animal_ID)) %>% 
#   dplyr::filter(!is.na(Animal_ID)) %>% 
#   dplyr::select(Animal_ID, Station_ID, Species) %>% 
#   tidyr::pivot_wider(names_from = Station_ID, values_from = Species,names_prefix = 'st_') %>% 
#   # Depending on the species that fills in each cell, let's replace Pekania p. with 1, and everything 
#   # else we'll replace with 0.
#   dplyr::mutate(across(-Animal_ID, \(x) dplyr::case_when(
#     is.na(x) ~ 0,
#     x == 'Pekania pennanti' ~ 1,
#     T ~ 0)))

# This is the full 'framework' of the matrix - it has ALL stations and animal IDs.
matrix_frame = DNA_columbian_3k %>% 
  ungroup() %>%
  sf::st_drop_geometry() %>% 
  dplyr::mutate(st_filler = 'goop') %>% 
  dplyr::select(Station_ID, Animal_ID, st_filler) %>%
  dplyr::distinct() %>% 
  tidyr::pivot_wider(names_from = Station_ID, values_from = st_filler, names_prefix = 'st_') %>% 
  dplyr::filter(!duplicated(Animal_ID)) %>%
  dplyr::filter(!is.na(Animal_ID))

n.ind.En <- n_distinct(na.omit(matrix_frame$Animal_ID)) 

# tbl_to_join = matrix_frame %>% 
#   dplyr::select(-all_of(stations_already_in_this_occasion)) %>% 
#   dplyr::mutate(across(-Animal_ID, \(x) x = NA))

# For each 'occasion', make a table showing 1's and 0's by unique animal ID and stations as columns.
tbls_by_occasion_list = 1:max(DNA_columbian_3k$occasion) %>% 
  map( ~ {
    dat_for_occ = DNA_columbian_3k %>% 
      dplyr::filter(occasion == .x) %>% 
      dplyr::ungroup() %>% 
      sf::st_drop_geometry() %>% 
      dplyr::filter(!duplicated(Animal_ID)) %>% 
      dplyr::filter(!is.na(Animal_ID)) %>% 
      dplyr::select(Animal_ID, Station_ID, Species) %>% 
      tidyr::pivot_wider(names_from = Station_ID, values_from = Species,names_prefix = 'st_') %>% 
      # Depending on the species that fills in each cell, let's replace Pekania p. with 1, and everything 
      # else we'll replace with 0.
      dplyr::mutate(dplyr::across(-Animal_ID, \(x) dplyr::case_when(
        is.na(x) ~ 0,
        x == 'Pekania pennanti' ~ 1,
        T ~ 0)))
    
    stations_already_in_this_occasion = names(dat_for_occ[,-1])
    
    tbl_to_join = matrix_frame %>% 
      dplyr::select(-all_of(stations_already_in_this_occasion)) %>% 
      dplyr::mutate(dplyr::across(-Animal_ID, \(x) x = 0))
    
    dat_for_occ %>% 
      dplyr::left_join(tbl_to_join) %>% 
      dplyr::bind_rows(tbl_to_join[!tbl_to_join$Animal_ID %in% dat_for_occ$Animal_ID,]) %>% 
      dplyr::mutate(dplyr::across(-Animal_ID, \(x) ifelse(is.na(x), 0, x)))
  })

names(tbls_by_occasion_list) <- paste0('occ_',c(1:max(DNA_columbian_3k$occasion)))

tbls_by_occasion_list$occ_1
# tbls_by_occasion_list$occ_8
tbls_by_occasion_list$occ_9

# Let's reorder the stations for this output based on the order in the 'matrix_frame'
# we made; the same thing for the order of rows by Animal_ID.
all_dat = tbls_by_occasion_list %>% 
  dplyr::bind_rows(.id = 'occasion') %>% 
  dplyr::arrange(occasion,Animal_ID) %>% 
  dplyr::select(occasion, Animal_ID, names(matrix_frame)[-1])

# Please use all_dat going forward!

k.En <- max(DNA_columbian_3k$occasion)
J.s.En <-length(colnames(matrix_frame))-1 # not counting Animal ID
## create Y.s 
Y.s <- array(0, c(n.ind.En, J.s.En, k.En)) 

all_dat %>% 
  tidyr::pivot_longer(cols = -c(occasion, Animal_ID)) %>%  
  dplyr::group_by(Animal_ID) %>% 
  dplyr::summarise(total_value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(total_value != 0)
# All 110 unique Animal_IDs have been detected at least once across the occasions and stations.

m_occ1 <- all_dat %>% filter(occasion== "occ_1") %>% select(-c(Animal_ID, occasion))
m_occ2 <- all_dat %>% filter(occasion== "occ_2") %>% select(-c(Animal_ID, occasion))
m_occ3 <- all_dat %>% filter(occasion== "occ_3") %>% select(-c(Animal_ID, occasion))
m_occ4 <- all_dat %>% filter(occasion== "occ_4") %>% select(-c(Animal_ID, occasion))
#m_occ5 <- all_dat %>% filter(occasion== "occ_5") %>% select(-c(Animal_ID, occasion))
# m_occ6 <- all_dat %>% filter(occasion== "occ_6") %>% select(-c(Animal_ID, occasion))
# m_occ7 <- all_dat %>% filter(occasion== "occ_7") %>% select(-c(Animal_ID, occasion))
# m_occ8 <- all_dat %>% filter(occasion== "occ_8") %>% select(-c(Animal_ID, occasion))
# m_occ9 <- all_dat %>% filter(occasion== "occ_9") %>% select(-c(Animal_ID, occasion))

Y.s[,,1] <- array(data=unlist(m_occ1), dim= c(n.ind.En, J.s.En))# remove Animal_ID
Y.s[,,2] <- array(data = unlist(m_occ2), dim= c(n.ind.En, J.s.En))
Y.s[,,3] <- array(data = unlist(m_occ3), dim= c(n.ind.En, J.s.En))
Y.s[,,4] <- array(data = unlist(m_occ4), dim= c(n.ind.En, J.s.En))
#Y.s[,,5] <- array(data = unlist(m_occ5), dim= c(n.ind.En, J.s.En))
# Y.s[,,6] <- array(data = unlist(m_occ6), dim= c(110, 1312))
# Y.s[,,7] <- array(data = unlist(m_occ7), dim= c(110, 1312))
# Y.s[,,8] <- array(data = unlist(m_occ8), dim= c(110, 1312))
# Y.s[,,9] <- array(data = unlist(m_occ9), dim= c(110, 1312))
# # 
# class(Y.s) <- "numeric"

# Y.s[is.na(Y.s)] <- 0

Y.s.enterprise <- Y.s


#Do we need to make all NA to 0? NA shows trap not active

# Session1_final[is.na(Session1_final)] <- 0
# Session2_final[is.na(Session2_final)] <- 0
# Session3_final[is.na(Session3_final)] <- 0


# prueba <- colSums(Session1_final[,-1 ])
# prueba2 <- colSums(Session2_final[,-1 ])
# prueba3 <- colSums(Session3_final[,-1 ])


## GET OPERABILITY MATRIX

Unique_traps <-  as.data.frame(unique(DNA_columbian_3k$Station_ID))
colnames(Unique_traps) <- "Station_ID"

Enterprise_matrix2 <- DNA_columbian_3k

session1_trap <- Enterprise_matrix2 %>% 
  ungroup() %>% 
  st_drop_geometry()%>% 
  select(Station_ID, occasion)%>% 
  subset(occasion==1) %>%
  dplyr::rename(Occ_1= occasion)

session2_trap <- Enterprise_matrix2 %>% 
  ungroup() %>% 
  st_drop_geometry()%>% 
  select(Station_ID, occasion)%>% 
  subset(occasion==2) %>%
  dplyr::rename(Occ_2= occasion)

session2_trap$Occ_2 <- 1

session3_trap <- Enterprise_matrix2 %>% 
  ungroup() %>% 
  st_drop_geometry()%>% 
  select(Station_ID, occasion)%>% 
  subset(occasion==3) %>%
  dplyr::rename(Occ_3= occasion)

session3_trap$Occ_3 <- 1

session4_trap <- Enterprise_matrix2 %>% 
  ungroup() %>% 
  st_drop_geometry()%>% 
  select(Station_ID, occasion)%>% 
  subset(occasion==4) %>%
  dplyr::rename(Occ_4= occasion)

session4_trap$Occ_4 <- 1

# session5_trap <- Enterprise_matrix2 %>% 
#   ungroup() %>% 
#   st_drop_geometry()%>% 
#   select(Station_ID, occasion)%>% 
#   subset(occasion==5) %>%
#   dplyr::rename(Occ_5= occasion)
# 
# session5_trap$Occ_5 <- 1


O.s_df <- dplyr::left_join(Unique_traps, session1_trap) %>% 
  dplyr::left_join(session2_trap) %>% 
  dplyr::left_join(session3_trap, by= "Station_ID")%>% 
  dplyr::left_join(session4_trap) #%>% 
  #dplyr::left_join(session5_trap, by= "Station_ID") #%>% replace(is.na(.), 0) 


O.s <- as.matrix(O.s_df[,-1]) ## operability matrix
O.s[is.na(O.s)] <- 0

# 6. Create O.o array for cameras ####

##6.1 read file from script 3. camera data prep/ camera grid short. ####
setwd("I:/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Mesocarnivore_Monitoring_Program/2.Data/Mesocarnivores DB/1. Master Data")

cameras_grid_short <- read_csv("columbian_cam_deployment_05-13-2024.csv")%>% drop_na(c(Latitude_DD, Longitude_DD)) %>%  dplyr::rename(End_Deployment_date2= New_End_Deployment_date) %>%  sf::st_as_sf(., coords= c('Longitude_DD', 'Latitude_DD'), crs= 4326, remove = FALSE ) %>% st_transform( crs=3005)%>% dplyr::mutate(UTM_X = sf::st_coordinates(.)[,1], UTM_Y = sf::st_coordinates(.)[,2])

#2.4 restrict by region 

cam_to_add <- cameras_grid_short %>% filter(Station_ID == "RC-05") # remove later
cameras_grid_short <- st_intersection( cameras_grid_short, chilcotin)
cameras_grid_short <- bind_rows(cameras_grid_short, cam_to_add)
# cameras_grid$Start_Deployment_date <- ymd(cameras_grid$Start_Deployment_date)
# cameras_grid$End_Deployment_date <- ymd(cameras_grid$End_Deployment_date)
# cameras_grid$Cam_days <- interval(cameras_grid$Start_Deployment_date, cameras_grid$End_Deployment_date)/ddays(1)
# # this does not include the first day, different from operation matrix

#remove any deployments without end dates

# tmp <- cameras_grid[is.na(cameras_grid$End_Deployment_date)==F,]

#List by project to reduce NAs and stack days

# tmp.list <- split(tmp, f= tmp$Project_name)

# wrapper <- function(df) {
#   as.data.frame %>% list %>% return
# }
# cameras_grid_list <-  cameras_grid %>% group_by(Project_name) %>%
#   do(res = wrapper(.)) 

# Create an empty list to store our months, all projects have at least one individual identified, so empty list will always be filled with something. 

#says daily but changed days to months :)

# length(tmp.list)
# 
# daily_lookup <- list()
# daily_lookup2 <- list()
# daily_lookup3 <- list()
# daily_lookup4 <- list()
# daily_lookup5 <- list()
# daily_lookup6 <- list()
# daily_lookup7 <- list()

#project with one camera is causing problems, it did not have any fisher detectoins so I removed it.

cameras_grid_short <- filter(cameras_grid_short, Project_name != "6115")

cam.list <- split(cameras_grid_short, f= cameras_grid_short$Project_name)

CampOp.list <- list()

#project with one camera is causing problems, it did not have any fisher detectoins so I removed it.

# only for chilcotin

#cam.list[[1]][1,3] <- "RC-05" 

setwd("I:/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Mesocarnivore_Monitoring_Program/2.Data/Mesocarnivores DB/1. Master Data")
for (i in 1:length(cam.list)){
  
  CampOp.list[[i]] <- cameraOperation(
    CTtable= cam.list[[i]],
    stationCol = "Station_ID",
    setupCol= "Start_Deployment_date",
    retrievalCol= "End_Deployment_date2",
    hasProblems = FALSE,
    dateFormat = "ymd",
    occasionStartTime = 0
    #writecsv = FALSE
  )
}



# 
# #Cam operation for project 1
# 
# for(i in 1:nrow(tmp.list[[1]])) {
#   if(ymd(tmp.list[[1]]$Start_Deployment_date[i])!=ymd(tmp.list[[1]]$End_Deployment_date[i]))
#   {
#     daily_lookup[[i]] <- data.frame("date"=seq(ymd(tmp.list[[1]]$Start_Deployment_date[i]), ymd(tmp.list[[1]]$End_Deployment_date[i]),  by= "months"), "Station_ID"=tmp.list[[1]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup1<- bind_rows(daily_lookup)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup1 <- row_lookup1[duplicated(row_lookup1)==F,]
# 
# # create operation matrix
# library(reshape2)
# row_lookup1$active <- 1
# CamOp1 <- dcast(row_lookup1, Station_ID ~ date, value.var = "active")
# # names1 <- c(paste0("V", 1:ncol(CamOp1)))
# # colnames(CamOp1) <- names1 # dont change the names yet cos you need this for a native join 
# # 
# #Cam operation for project 2
# 
# for(i in 1:nrow(tmp.list[[2]])) {
#   if(ymd(tmp.list[[2]]$Start_Deployment_date[i])!=ymd(tmp.list[[2]]$End_Deployment_date[i]))
#   {
#     daily_lookup2[[i]] <- data.frame("date"=seq(ymd(tmp.list[[2]]$Start_Deployment_date[i]), ymd(tmp.list[[2]]$End_Deployment_date[i]),  by= "months"), "Station_ID"=tmp.list[[2]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup2<- bind_rows(daily_lookup2)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup2 <- row_lookup2[duplicated(row_lookup2)==F,]
# 
# # create operation matrix
# library(reshape2)
# row_lookup2$active <- 1
# CamOp2 <- dcast(row_lookup2, Station_ID ~ date, value.var = "active")
# # names2 <- c(paste0("V", 1:ncol(CamOp2)))
# # colnames(CamOp2) <- names2
# 
# #Cam operation for project 3
# 
# for(i in 1:nrow(tmp.list[[3]])) {
#   if(ymd(tmp.list[[3]]$Start_Deployment_date[i])!=ymd(tmp.list[[3]]$End_Deployment_date[i]))
#   {
#     daily_lookup3[[i]] <- data.frame("date"=seq(ymd(tmp.list[[3]]$Start_Deployment_date[i]), ymd(tmp.list[[3]]$End_Deployment_date[i]),  by= "month"), "Station_ID"=tmp.list[[3]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup3<- bind_rows(daily_lookup3)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup3 <- row_lookup3[duplicated(row_lookup3)==F,]
# 
# # create operation matrix
# row_lookup3$active <- 1
# CamOp3 <- dcast(row_lookup3, Station_ID ~ date, value.var = "active")
# # names3 <- c(paste0("V", 1:ncol(CamOp3)))
# # colnames(CamOp3) <- names3
# 
# #Cam operation for project 4
# 
# for(i in 1:nrow(tmp.list[[4]])) {
#   if(ymd(tmp.list[[4]]$Start_Deployment_date[i])!=ymd(tmp.list[[4]]$End_Deployment_date[i]))
#   {
#     daily_lookup4[[i]] <- data.frame("date"=seq(ymd(tmp.list[[4]]$Start_Deployment_date[i]), ymd(tmp.list[[4]]$End_Deployment_date[i]),  by= "months"), "Station_ID"=tmp.list[[4]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup4<- bind_rows(daily_lookup4)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup4 <- row_lookup4[duplicated(row_lookup4)==F,]
# 
# # create operation matrix
# row_lookup4$active <- 1
# CamOp4 <- dcast(row_lookup4, Station_ID ~ date, value.var = "active")
# # names4 <- c(paste0("V", 1:ncol(CamOp4)))
# # colnames(CamOp4) <- names4
# 
# #Cam operation for project 5
# 
# for(i in 1:nrow(tmp.list[[5]])) {
#   if(ymd(tmp.list[[5]]$Start_Deployment_date[i])!=ymd(tmp.list[[5]]$End_Deployment_date[i]))
#   {
#     daily_lookup5[[i]] <- data.frame("date"=seq(ymd(tmp.list[[5]]$Start_Deployment_date[i]), ymd(tmp.list[[5]]$End_Deployment_date[i]),  by= "months"), "Station_ID"=tmp.list[[5]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup5<- bind_rows(daily_lookup5)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup5 <- row_lookup5[duplicated(row_lookup5)==F,]
# 
# # create operation matrix
# row_lookup5$active <- 1
# CamOp5 <- dcast(row_lookup5, Station_ID ~ date, value.var = "active")
# # names5 <- c(paste0("V", 1:ncol(CamOp5)))
# # colnames(CamOp5) <- names5
# 
# #Cam operation for project 6
# 
# for(i in 1:nrow(tmp.list[[6]])) {
#   if(ymd(tmp.list[[6]]$Start_Deployment_date[i])!=ymd(tmp.list[[6]]$End_Deployment_date[i]))
#   {
#     daily_lookup6[[i]] <- data.frame("date"=seq(ymd(tmp.list[[6]]$Start_Deployment_date[i]), ymd(tmp.list[[6]]$End_Deployment_date[i]),  by= "months"), "Station_ID"=tmp.list[[6]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup6<- bind_rows(daily_lookup6)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup6 <- row_lookup6[duplicated(row_lookup6)==F,]
# 
# # create operation matrix
# library(reshape2)
# row_lookup6$active <- 1
# CamOp6 <- dcast(row_lookup6, Station_ID ~ date, value.var = "active")
# names6 <- c(paste0("V", 1:ncol(CamOp6)))
# colnames(CamOp6) <- names6

# #Cam operation for project 7
# 
# for(i in 1:nrow(tmp.list[[7]])) {
#   if(ymd(tmp.list[[7]]$Start_Deployment_date[i])!=ymd(tmp.list[[7]]$End_Deployment_date[i]))
#   {
#     daily_lookup7[[i]] <- data.frame("date"=seq(ymd(tmp.list[[7]]$Start_Deployment_date[i]), ymd(tmp.list[[7]]$End_Deployment_date[i]),  by= "days"), "Station_ID"=tmp.list[[7]]$Station_ID[i])}
# }
# 
# # Merge the lists into a dataframe
# row_lookup7<- bind_rows(daily_lookup7)
# 
# # Remove duplicates - when start and end days are the same for successive deployments
# row_lookup7 <- row_lookup7[duplicated(row_lookup7)==F,]
# 
# # create operation matrix
# row_lookup7$active <- 1
# CamOp7 <- dcast(row_lookup7, Station_ID ~ date, value.var = "active")
# # names7 <- c(paste0("V", 1:ncol(CamOp7)))
# # colnames(CamOp7) <- names7
# 
# # list all matrices
# 
# CamOpT <- list(CamOp1, CamOp2, CamOp3, CamOp4, CamOp5, CamOp6)
# 
# 
# for(i in 1:length(CamOpT)){
# 
#   CamOpT[[i]][CamOpT[[i]]==1] <- 0
# }

# 9.2. Create occupancy matrix = O.o ###

## 6.2. Read in species data  ####

setwd("I:/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Mesocarnivore_Monitoring_Program/2.Data/Mesocarnivores DB/1. Master Data")

sp_cam <- read_csv("fisher_detection_cam_05-13-2024.csv", col_types = cols(Date_Time = col_datetime()))
sp_cam <- as.data.frame(sp_cam)
# Convert DateTime to character to run next script

sp_cam$Date_Time2 <- as.character(sp_cam$Date_Time)

sp_cam <- sp_cam %>% drop_na(Date_Time2) 

## 6.3. Create table with independent events ####
out_dir<- "I:/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Mesocarnivore_Monitoring_Program/2.Data/Mesocarnivores DB/1. Master Data"

Inp_eventscam1 <-  filterRecordTable(
  sp_cam,
  minDeltaTime = 24*60,
  deltaTimeComparedTo= "lastIndependentRecord",
  speciesCol = "Species",
  stationCol= "Station_ID",
  camerasIndependent = FALSE,
  recordDateTimeCol = "Date_Time2",
  recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
  removeDuplicateRecords = TRUE,
  # exclude= c("Unknown", "Other", "Setup"),
  timeZone= "US/Pacific",
  writecsv = TRUE,
  outDir= out_dir
)

#remove events outsite of the target area omineca

Inp_eventscam <- Inp_eventscam1 %>% filter(Project_name %in% c("Fisher Exclusion Box Program", "Interior BC Moose & Predator Camera Study", "Itcha Ilgachuz", "South Chilcotin Mountains Camera Sampling")) %>% 
  filter(Station_ID %in% c("RC-05","RC-06","PGS-01","PGS-02",     "PGS-03",     "PGS-04",     "PGS-05",     "PGS-06",     "PGS-07","PGS-09", "PGS-11" ,    "PGS-12",     "PGS-14",     "PGS-15",     "PGS-16",     "PGS-17", "PGS-18","PGS-19","PGS-20",  "PGS-25", "PGS-26", "PGS-28",     "PGS-29",     "PGS-30",     "PGS-31","PGS-32","PGS-34","PGS-35","PGS-38",     "PGS-39",     "PGS-40",     "PGS-41",     "PGS-42","PGS-43",  "PGS-44",     "PGS-45",     "PGS-46",     "PGS-47"  ,   "PGS-48"  ,   "ITCHA1-01",  "ITCHA1-07", "ITCHA1-25",  "ITCHA2-02",  "ITCHA2-03",  "ITCHA2-04",  "ITCHA2-05",  "ITCHA2-06",  "ITCHA2-07",  "ITCHA2-08", "ITCHA2-09",  "ITCHA2-10",  "ITCHA2-11" , "ITCHA2-12" , "ITCHA2-15",  "ITCHA2-16",  "ITCHA2-17"  ,"ITCHA2-25" ,"ITCHA2-26",  "ITCHA2-27",  "ITCHA2-27B", "ITCHA2-28",  "ITCHA2-29"  ,"ITCHA2-30",  "ITCHA3-01" , "ITCHA3-07" ,"ITCHA3-24",  "ITCHA3-30",  "ITCHA4-01",  "ITCHA4-02" , "ITCHA4-03",  "ITCHA4-04" , "ITCHA4-05",  "ITCHA4-06" ,"ITCHA4-07",  "ITCHA4-08",  "ITCHA4-09" , "ITCHA4-10",  "ITCHA4-11" , "ITCHA4-12" , "ITCHA4-13" , "ITCHA4-14" ,"ITCHA4-15",  "ITCHA4-16",  "ITCHA4-17",  "ITCHA4-18",  "ITCHA4-19",  "ITCHA4-20",  "ITCHA4-21",  "ITCHA4-22" ,"ITCHA4-22B", "ITCHA4-23",  "ITCHA4-24", "ITCHA4-25",  "ITCHA4-25B", "ITCHA4-26" , "ITCHA4-27",  "ITCHA4-28","ITCHA4-29",  "ITCHA4-30",  "ITCHA5-01", "ITCHA5-01B", "ITCHA5-02",  "ITCHA5-03",  "ITCHA5-04",  "ITCHA5-05", "ITCHA5-06",  "ITCHA5-07",  "ITCHA5-08", "ITCHA5-09",  "ITCHA5-10",  "ITCHA5-11" , "ITCHA5-12",  "ITCHA5-12B","ITCHA5-13",  "ITCHA5-13B", "ITCHA5-14",  "ITCHA5-14B", "ITCHA5-15",  "ITCHA5-16",  "ITCHA5-17" , "ITCHA5-18","ITCHA5-19",  "ITCHA5-20",  "ITCHA5-21" , "ITCHA5-22",  "ITCHA5-23" , "ITCHA5-24" , "ITCHA5-25",  "ITCHA5-26","ITCHA5-27",  "ITCHA5-28",  "ITCHA5-29",  "ITCHA5-30",  "ITCHA6-01",  "ITCHA6-02",  "ITCHA6-03",  "ITCHA6-04","ITCHA6-05",  "ITCHA6-06",  "ITCHA6-07",  "ITCHA6-08",  "ITCHA6-09",  "ITCHA6-10",  "ITCHA6-11",  "ITCHA6-12","ITCHA6-13",  "ITCHA6-14",  "ITCHA6-15",  "ITCHA6-16",  "ITCHA6-17",  "ITCHA6-17B", "ITCHA6-18",  "ITCHA6-19","ITCHA6-20" , "ITCHA6-21",  "ITCHA6-22",  "ITCHA6-23",  "ITCHA6-24",  "ITCHA6-24B", "ITCHA6-25",  "ITCHA6-26","ITCHA6-27",  "ITCHA6-28",  "ITCHA6-29",  "ITCHA6-30",  "SCE5",       "SCE6" ,      "SCE7" ,      "SCF5" ,"SCF5_2",     "SCF6",       "SCF7",       "SCF8",       "SCF9",       "SCG4" ,      "SCG9",       "SCH10" ,"SCH10_2",    "SCH11",      "SCH11_2",    "SCH5",       "SCH9",       "SCI10" ,     "SCI11"  ,    "SCI3"  ,"SCI4"   ,    "SCI8" ,      "SCI9" ,      "SCJ10" ,     "SCJ10a" ,    "SCJ10a_2" ,  "SCJ4" ,      "SCJ9"  ,"SCK10"   ,   "SCK10_2"  ,  "SCK3"  ,     "SCK4"    ,   "SCK5",       "SCK6" ,      "SCK7" ,      "SCK8"      ,"SCK9" ,      "SCL10"    ,  "SCL3"    ,   "SCL6"  ,     "SCL6_2",     "SCL7" ,      "SCL8"  ,     "SCL8_2"   ,"SCL9",       "SCM9"  ,     "SCM9_2")) # added RC-delete to generate empty matrix, need to remove it later


Inp_eventscam.list <- split(Inp_eventscam, f= Inp_eventscam$Project_name)
DetHist <- list()

#detection history only for one site with detected fisher. 


for(i in 1:length(Inp_eventscam.list)){
  DetHist[[i]] <- detectionHistory(recordTable = Inp_eventscam.list[[i]],
                                   camOp = CampOp.list[[i]],
                                   stationCol = "Station_ID",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "Date_Time",
                                   species = "Pekania pennanti",
                                   occasionLength = 30,
                                   day1 = "station",
                                   datesAsOccasionNames = FALSE,
                                   includeEffort = FALSE,
                                   timeZone = "America/Vancouver")
}



# make one table

df <- ldply(DetHist, data.frame) # this doesnt add the rownmaes so need to make sure stations are the same on X.o when mapped. 

#remove cam_to_add from detHist

dim(df)
df <- df[-1,] #removes added station
# 
# ## filter by species
# 
# fisher_cam <- subset(Inp_eventscam, Species== "Pekania pennanti")
# 
# # generate a list making each project an element of the list. 
# 
# tmp_cam <- split(Inp_eventscam, f= Inp_eventscam$Project_name)
# length(tmp_cam)  ## all detections
# 
# tmp_fisher_cam <- split(fisher_cam, f= fisher_cam$Project_name) ## only fisher detections
# length(tmp_fisher_cam) # some projects had no fisher detections, make those 0 for all active days
# 
# 
# # generate the detection matrix of fisher camera detections
# 
# fisher_cam_matrix <- list()
# 
# for (i in 1:length(tmp_fisher_cam)){
#   
#   fisher_cam_matrix[[i]] <- tmp_fisher_cam[[i]]  %>% select(Station_ID, Date, Count) %>%  pivot_wider(id_cols = Station_ID, names_from = Date, values_from = Count)
#   
# }
# 
# 
# # populate empty with fisher table
# 
# for(i in 1:length(CamOpT)){
#   CamOpT[[i]] <- CamOpT[[i]] %>% remove_rownames %>% column_to_rownames(var="Station_ID")
# }
# 
# for(i in 1:length(fisher_cam_matrix)){
#   fisher_cam_matrix[[i]] <- fisher_cam_matrix[[i]] %>% remove_rownames %>% column_to_rownames(var="Station_ID")
# }
# 
# #we create a new matrix as big as CamOpT with the values of m.2 in it
# 
# # loop didnt work because of problem with one of the matrices having different dimensions. Try again. 
# # 
# # mres <- list()
# # cam_matrix <- list()
# # 
# # for( i in 1:length(fisher_cam_matrix)){
# #   for(j in 3:length(CamOpT)){
# #     
# #     mres[[i]]<-CamOpT[[j]]
# #     mres[[i]][rownames(fisher_cam_matrix[[i]]),colnames(fisher_cam_matrix[[i]])]<-fisher_cam_matrix[[i]]
# #     #Then we use pmax
# #     cam_matrix[[i]] <- pmax(CamOpT[[j]],mres[[i]],na.rm=FALSE)
# #   }
# #   
# # }
# 
# mres2<-CamOpT[[2]]
# mres2[rownames(fisher_cam_matrix[[1]]),colnames(fisher_cam_matrix[[1]])]<-fisher_cam_matrix[[1]]
# #Then we use pmax
# cam_matrix_2 <- pmax(CamOpT[[2]],mres2,na.rm=FALSE)
# names2 <- c(paste0("V", 1:ncol(cam_matrix_2)))
# colnames(cam_matrix_2) <- names2
# 
# mres3<-CamOpT[[3]]
# mres3[rownames(fisher_cam_matrix[[2]]),colnames(fisher_cam_matrix[[2]])]<-fisher_cam_matrix[[2]]
# #Then we use pmax
# cam_matrix_3 <- pmax(CamOpT[[3]],mres3,na.rm=FALSE)
# names3 <- c(paste0("V", 1:ncol(cam_matrix_3)))
# colnames(cam_matrix_3) <- names3
# 
# mres4<-CamOpT[[4]]
# mres4[rownames(fisher_cam_matrix[[3]]),colnames(fisher_cam_matrix[[3]])]<-fisher_cam_matrix[[3]]
# #Then we use pmax
# cam_matrix_4 <- pmax(CamOpT[[4]],mres4,na.rm=FALSE)
# names4 <- c(paste0("V", 1:ncol(cam_matrix_4)))
# colnames(cam_matrix_4) <- names4
# 
# mres5<-CamOpT[[5]]
# mres5[rownames(fisher_cam_matrix[[4]]),colnames(fisher_cam_matrix[[4]])]<-fisher_cam_matrix[[4]]
# #Then we use pmax
# cam_matrix_5 <- pmax(CamOpT[[5]],mres5,na.rm=FALSE)  ## one camera with a fisher detectoin was removed from the locations, based on the map?
# names5 <- c(paste0("V", 1:ncol(cam_matrix_5)))
# colnames(cam_matrix_5) <- names5
# 
# mres6<-CamOpT[[6]]
# mres6[rownames(fisher_cam_matrix[[5]]),colnames(fisher_cam_matrix[[5]])]<-fisher_cam_matrix[[5]]
# #Then we use pmax
# cam_matrix_6 <- pmax(CamOpT[[6]],mres6,na.rm=FALSE)
# names6 <- c(paste0("V", 1:ncol(cam_matrix_6)))
# colnames(cam_matrix_6) <- names6
# 
# ## Remove the names of all columns and combine all matrices 
# 
# names1 <- c(paste0("V", 1:ncol(CamOpT[[1]])))
# colnames(CamOpT[[1]]) <- names1
# 
# names2 <- c(paste0("V", 1:ncol(CamOpT[[2]])))
# colnames(mres2) <- names2
# 
# names3 <- c(paste0("V", 1:ncol(CamOpT[[3]])))
# colnames(mres3) <- names3
# 
# names4 <- c(paste0("V", 1:ncol(CamOpT[[4]])))
# colnames(mres4) <- names4
# 
# names5 <- c(paste0("V", 1:ncol(CamOpT[[5]])))
# colnames(mres5) <- names5
# 
# names6 <- c(paste0("V", 1:ncol(CamOpT[[6]])))
# colnames(mres6) <- names6


### merge all camera operation matrices
library(plyr)
O.o2 <- as.matrix(df[]) # will be operability matrix
O.o <- as.matrix(df[]) 
dim(O.o2)


O.o2[O.o2==0] <- 1
O.o2 <- replace(O.o2, is.na(O.o2), 0)### Operability matrix

#O.o <- replace(O.o, is.na(O.o), 0) ## detection matrix


# Add an extra dimension of size 1
O <- array(0, c(dim(O.o)[1], dim(O.o)[2], T)) 

#namesEnterprise <- c(paste0("V", 1:J.s.En))
colnames(O[,,1]) <-NULL

O[,,1] <- as.matrix(O.o)

## add coordinates for camera
cameras_grid_short <- cameras_grid_short %>% filter(Station_ID != "RC-05")# removed added station

X.o <- cameras_grid_short %>% select(Station_ID, UTM_X, UTM_Y) %>%  st_drop_geometry()

X.o <-X.o[!duplicated(X.o), ]  %>%  ungroup() %>% select(UTM_X,UTM_Y)
X.o <- X.o/coord.scale
X.o <- as.matrix(X.o, dim = c(dim(X.s)[1], 2))
X.o <- unname(X.o)

# Define limits of the state-space and add it to the plot

traps.scale <- rbind(X.o,X.s)

library(scales)
#coordinates <- select(Ssites_cariboo2, Long, Lat) ## scaled coordinates
plot(traps.scale, xlab='V1', ylab='V2', frame=FALSE, las=1, pch=10, col='#002A64FF', asp=1)
summary(traps.scale)
xlim <- c(min(traps.scale[,1])+2, max(traps.scale[,1])+2)
ylim <- c(min(traps.scale[,2])+2, max(traps.scale[,2])+2)


## number of occassion per trap 

n.occ.s1 <- DNA_columbian_3k %>% 
          dplyr::group_by(Station_ID) %>% 
          dplyr::summarise(N_occ= max(occasion)) %>% 
          st_drop_geometry()
n.occ.s <- n.occ.s1[,2]
n.occ.s <- as.numeric(n.occ.s$N_occ)

n.occ.o <- rowSums(O.o2, na.rm=TRUE)



## other variables needed for the model 
M <- 1500
T=1
psi <- 0.3

s <- array(NA, c(M, 2, T)) # empty array to fill with activity centers, 300 ind, sampled across 4 periods, 2 columns for coordinates. 
z <- a <- matrix(NA, M, T) # empty matrix for population membership
s[,,1] <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])) # random activity centers
z[,1] <- rbinom(M, 1, psi) # create first year's pop with M from psi.
##
dim(s)

dim(Y.s)
# We'll need to put the things we want exported from inside the function
# into a 'return' statement. Ideally, we should name each of the exports.
#return(c(s,Z,Y.s = 'Y.s'))
#}

# Run the function!
#results = process_data()

