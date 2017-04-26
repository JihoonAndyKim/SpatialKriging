rm(list=ls())
library(ncdf4)
library(abind)
library(fields)
library(mapproj)
library(geoR)
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
plot(lon.lat)
quilt.plot(lon.lat, c(temp[,,1,1]))
US(add=TRUE)

temp1 <- mapproject(x=lon.lat[,1],y=lon.lat[,2],projection="sinusoidal")
lon.lat <- cbind(temp1$x,temp1$y)
mod <- lm(c(temp[,,1,1]) ~ lon.lat[,1] + lon.lat[,2]) # does both lat/lon
summary(mod)

z <- mod$residuals
lonvals <- seq(min(lon.lat[,1]),max(lon.lat[,1]),length.out=200)
latvals <- seq(min(lon.lat[,2]),max(lon.lat[,2]),length.out=200)

pred.grd <- as.matrix(expand.grid(lonvals,latvals))
quilt.plot(lon.lat, z)

breaks <- seq(0,max(rdist(lon.lat)),length.out=50)

v.bc <- variog(coords=lon.lat,data=z,estimator.type="classical",breaks=breaks,
               bin.cloud=TRUE) # binned

plot(v.bc,main="binned",pch=19,cex=0.5)
vfit <- variofit(v.bc,cov.model="exponential",weights="cressie")
lines(vfit)

a = vfit$cov.pars[2]
sigma = sqrt(vfit$cov.pars[1])
tausq = vfit$nugget

## Matrices for kriging
dist0.mat <- rdist(lon.lat,pred.grd)
Sigma0 <- sigma^2 * exp(-dist0.mat/a)

dist.mat <- rdist(lon.lat)
Sigma <- sigma^2 * exp(-dist.mat/a) + tausq

## Simple kriging predictor
z.hat <- t(Sigma0) %*% solve(Sigma) %*% z

#Let's take a look at the kriged data
c <- mod$coefficients
quilt.plot(pred.grd, c[1] + c[2]*pred.grd[,1] + c[3]*pred.grd[,2] + z.hat)
