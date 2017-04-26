rm(list=ls())
library(ncdf4)
library(abind)
library(fields)

options(scipen=999)

data <- nc_open("data/gefs_temp")

init.time <- ncvar_get(data, "intTime")
valid.time <- ncvar_get(data, "intValidTime")

#Temperature is in kelvin
temp <- ncvar_get(data, "Temperature_height_above_ground") - 273.15
lat <- data$dim$lat$vals
lon <- data$dim$lon$vals - 360

lon.lat <- as.matrix(expand.grid(lon,lat))

obs <- cbind(c(temp[,,1,1]), lon.lat)
quilt.plot(lon.lat, c(temp[,,1,1]))
US(add=TRUE)