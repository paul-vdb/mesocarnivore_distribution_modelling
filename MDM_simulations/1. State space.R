### Defining state space for Columbian populations ###

library(terra)
library(ggplot2)
library(sf)
library(readr)
library(farver)
library(ggplot2)
library(rmapshaper)
library(stars)
library(tidyterra)
require(sf)

#1. Read and plot layers ####

setwd("C:/LocalR")
populations <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Fisher_population_units")
barrier <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Geographic_barrier")
subpopulations <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Subpopulations")
subpop_bound <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Subpopulation_boundaries")
Interdigit <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Interdigitation_zone")
boreal <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Boreal_population")
columbian <- sf::st_read("BC_Fisher_populations_2024.gdb", layer = "Columbian_population")
grid_columbian <- st_read("BC_meso_grid_columbian.shp")
columbian_area <- st_union(grid_columbian, by_feature = FALSE)
columbian_area <- ms_simplify(columbian_area, keep = 0.01,
                      keep_shapes = FALSE)
#plot

bc_bound = bcmaps::bc_bound()
regs = bcmaps::nr_regions()

plot1 = ggplot() +
  geom_sf(data = bc_bound) +
  geom_sf(data = populations, aes(fill= Pop_unit)) +
  geom_sf(data = barrier, color = 'black') 
plot1


subpop_plot = ggplot() +
  geom_sf(data = bc_bound) +
  geom_sf(data = subpopulations, aes(fill= Subpop))+
  geom_sf(data = subpop_bound, color = 'black')
subpop_plot

plot3 = ggplot() +
  geom_sf(data = bc_bound) +
  geom_sf(data = Interdigit, color= 'yellow')
plot3

#2. simplify shape file of  subpopulations ####

library(rmapshaper)
subpop <- ms_simplify(subpopulations, keep = 0.001,
                                keep_shapes = FALSE)

subpop_plot = ggplot() +
  geom_sf(data = bc_bound) +
  geom_sf(data = subpop, aes(fill= Subpop))+
  geom_sf(data = subpop_bound, color = 'black')+
  coord_sf(crs = 3005) #3005 albers
subpop_plot

#3. Select desired population or Remove boreal population

subpop2 <- subpop |> dplyr::filter(Subpop == "Cariboo")
#subpop2 <- subpop |> dplyr::filter(Subpop != "Boreal")

#3. convert to raster, 3km (res=3000) is too small for area ####

template = rast(vect(subpop),res=3000)

subpop_raster <- rasterize(vect(subpop2), field= "Subpop", template)
plot(subpop_raster)

rasterplot = ggplot() +
  geom_spatraster(data = subpop_raster)+
  #geom_sf(data = bc_bound, color= 'black', fill= NA) +
  geom_sf(data = columbian_area, fill= NA)+
  #geom_sf(data = grid_columbian, fill= NA)+
  coord_sf(crs = 3005)  
rasterplot 

#4. Clip extent to columbian area, can also dissolve 3 subpopulations and clip extent with that 

subpop_raster_cp <- terra::mask(subpop_raster, columbian_area)
subpop_raster_cp <- trim(subpop_raster_cp)
subpop_raster_cp 
writeRaster(subpop_raster_cp,'Cariboo_3k.tif')

rasterplot2 = ggplot() +
  geom_spatraster(data = subpop_raster_cp)+
  #geom_sf(data = bc_bound, color= 'black', fill= NA) +
  geom_sf(data = columbian_area, fill= NA)+
  #geom_sf(data = grid_columbian, fill= NA)+
  geom_sf(data = Camera_data_sf, color = 'lightgreen') +
  geom_sf(data = DNA_data_sf, color= 'lightblue', alpha=0.1) +
  geom_sf(data = MDM_academics_sf, color= 'lightgreen', alpha=0.1) +
  #geom_sf(data = Academic_sp_fisher_sf, color= 'black', alpha=0.1)+
  #geom_sf(data = CAM_sp_map_fisher_sf, color= 'blue', alpha=0.1)+
  geom_sf(data = df_fisher_sf, color= 'pink') +
  coord_sf(crs = 3005)  
rasterplot2 


