rm(list=ls())
library(ncdf4)
library(abind)
library(fields)

options(scipen=999)

data <- nc_open("data/air.mon.mean.v401.nc")
air <- ncvar_get(data, "air")
lat <- data$dim$lat$vals
lon <- data$dim$lon$vals - 360

top = 81
left = 470
right = 587
bottom = 132
start = 1140
end = 1380

#Get data for the continental US for the past 20 years
us.air <- air[left:right, top:bottom, start:end]
lon.lat <- as.matrix(expand.grid(lon[left:right],lat[top:bottom]))

obs <- cbind(c(us.air[,,1]), lon.lat)
quilt.plot(lon.lat, c(us.air[,,10]))
US(add=TRUE)


# save(testWind, file = "GEFS_NREL_Wind_Data.RData")
# save(testDirWind, file = "GEFS_NREL_Wind_Dir_Data.RData")
#nc_close(data)
#rm(data)



