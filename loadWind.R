rm(list=ls())
library(ncdf4)
library(abind)
library(fields)
library(mapproj)
library(geoR)
options(scipen=999)

u <- nc_open("data/u")
v <- nc_open("data/v")

u.init.time <- ncvar_get(u, "intTime")
v.init.time <- ncvar_get(v, "intTime")

#Temperature is in kelvin
u.wind <- ncvar_get(u, "U-component_of_wind_height_above_ground_80m")
v.wind <- ncvar_get(v, "V-component_of_wind_height_above_ground_80m")
absoluteWind <- sqrt(u.wind^2 + v.wind^2)

lat <- u$dim$lat$vals
lon <- u$dim$lon$vals - 360

lon.lat <- as.matrix(expand.grid(lon,lat))

obs <- cbind(c(absoluteWind), lon.lat)
plot(lon.lat)
quilt.plot(lon.lat, c(absoluteWind))
US(add=TRUE)

resolution <- 250

temp1 <- mapproject(x=lon.lat[,1],y=lon.lat[,2],projection="mollweide")
lon.lat <- cbind(temp1$x,temp1$y)
mod <- lm(c(absoluteWind) ~ lon.lat[,1] + lon.lat[,2]) # does both lat/lon
summary(mod)

z <- mod$residuals


lonvals <- seq(min(lon.lat[,1]),max(lon.lat[,1]),length.out=resolution)
latvals <- seq(min(lon.lat[,2]),max(lon.lat[,2]),length.out=resolution)

pred.grd <- as.matrix(expand.grid(lonvals,latvals))
quilt.plot(lon.lat, z)

breaks <- seq(0,max(rdist(lon.lat)),length.out=50)

v.bc <- variog(coords=lon.lat,data=z,estimator.type="classical",breaks=breaks,
               bin.cloud=TRUE) # binned

plot(v.bc,main="binned",pch=19,cex=0.5)
#seems that cauchy is the best
vfit <- variofit(v.bc,cov.model="spherical",weights="cressie")
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
pred <- c[1] + c[2]*pred.grd[,1] + c[3]*pred.grd[,2] + z.hat
quilt.plot(pred.grd[,1], pred.grd[,2], pred, nx = resolution, ny = resolution)


library(rgdal)
reverse <- function(x, y, z) {
  my.points = data.frame(x=x * 18040096/2 , y=y * 9020048)
  my.points = SpatialPoints(my.points, CRS('+proj=moll'))
  my.points = as.data.frame(spTransform(my.points, CRS('+proj=longlat')))
  my.points$x = my.points$x + 264
  my.points$x[my.points$x > 180] = my.points$x[my.points$x > 180] - 360
  return(my.points)
}
new.loc <- reverse(pred.grd[,1], pred.grd[,2])


quilt.plot(new.loc, pred, nx = resolution, ny = resolution)
US(add=TRUE)

#NREL Datatower, 38842
#nrel <- which(new.loc[,1] < -105 & new.loc[,1] > -105.5 & new.loc[,2] > 39.75 & new.loc[,2] < 40)
