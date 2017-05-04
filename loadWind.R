rm(list=ls())
library(ncdf4)
library(abind)
library(fields)
library(mapproj)
library(geoR)
options(scipen=999)

resolution <- 250
u <- nc_open("data/u.nc")
v <- nc_open("data/v.nc")

u.init.time <- ncvar_get(u, "intTime")
v.init.time <- ncvar_get(v, "intTime")

#Temperature is in kelvin
u.wind <- ncvar_get(u, "U-component_of_wind_height_above_ground")
v.wind <- ncvar_get(v, "V-component_of_wind_height_above_ground")
absoluteWind <- sqrt(u.wind^2 + v.wind^2)

lat <- u$dim$lat$vals
lon <- u$dim$lon$vals - 360

lon.lat <- as.matrix(expand.grid(lon,lat))

obs <- cbind(c(absoluteWind), lon.lat)
quilt.plot(lon.lat, c(absoluteWind))
US(add=TRUE)

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
vfit <- variofit(v.bc,cov.model="cauchy",weights="cressie")
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


quilt.plot(new.loc, pred, nx = resolution, ny = resolution, zlim = c(0,15))
US(add=TRUE)



#---------------------------------------------------------
#------------------ OBSERVATION DATA ---------------------
#---------------------------------------------------------

#NREL Datatower, 38842
#nrel <- which(new.loc[,1] < -105 & new.loc[,1] > -105.5 & new.loc[,2] > 39.75 & new.loc[,2] < 40)
obs <- read.csv("data/newobs.csv")
obs <- obs[which(obs[,4] != 'M'),]
obs <- obs[, c(3,2,4)]
top = 51 # north lat
left = -127 # west long
right = -65 # east long
bottom =  23 # south lat
#convert factors to numeric and as well as converting knots to m/s
obs[,3] <- as.numeric(as.character(obs[,3])) * 0.514444

obs.index <- which(obs[,1] > right | obs[,1] < left | obs[,2] > top | obs[,2] < bottom)
obs <- obs[-obs.index,]
quilt.plot(obs[,1], obs[,2], obs[,3])
US(add=TRUE)

obs.lon.lat <- cbind(obs[,1], obs[,2])
temp1 <- mapproject(x=obs.lon.lat[,1],y=obs.lon.lat[,2],projection="mollweide")
obs.lon.lat <- cbind(temp1$x,temp1$y)
obs.mod <- lm(obs[,3] ~ obs.lon.lat[,1] + obs.lon.lat[,2]) # does both lat/lon
summary(obs.mod)

obs.z <- obs.mod$residuals
quilt.plot(obs.lon.lat, obs.z)

breaks <- seq(0,max(rdist(obs.lon.lat)),length.out=50)

v.bc <- variog(coords=obs.lon.lat,data=obs.z,estimator.type="classical",breaks=breaks,
               bin.cloud=TRUE) # binned

plot(v.bc,main="binned",pch=19,cex=0.5)
#seems that exp is the best
vfit <- variofit(v.bc,cov.model="exponential",weights="cressie")
lines(vfit)

a = vfit$cov.pars[2]
sigma = sqrt(vfit$cov.pars[1])
tausq = vfit$nugget

## Matrices for kriging
dist0.mat <- rdist(obs.lon.lat,pred.grd)
Sigma0 <- sigma^2 * exp(-dist0.mat/a)

dist.mat <- rdist(obs.lon.lat)
Sigma <- sigma^2 * exp(-dist.mat/a) + tausq

## Simple kriging predictor
obs.z.hat <- t(Sigma0) %*% solve(Sigma) %*% obs.z

#Let's take a look at the kriged data
c <- obs.mod$coefficients
obs.pred <- c[1] + c[2]*pred.grd[,1] + c[3]*pred.grd[,2] + obs.z.hat
quilt.plot(pred.grd[,1], pred.grd[,2], obs.pred, nx = resolution, ny = resolution)


quilt.plot(new.loc, obs.pred, nx = resolution, ny = resolution)
US(add=TRUE)

par(mfrow=c(1,2))
quilt.plot(new.loc, pred, nx = resolution, ny = resolution)
US(add=TRUE)
quilt.plot(new.loc, obs.pred, nx = resolution, ny = resolution)
US(add=TRUE)

MSE = sum((pred-obs.pred)^2)


#---------------------------------------------------------
#-------------- GLS with Ordinary Kriging ----------------
#---------------------------------------------------------


#NREL Datatower, 38842
#nrel <- which(new.loc[,1] < -105 & new.loc[,1] > -105.5 & new.loc[,2] > 39.75 & new.loc[,2] < 40)
obs <- read.csv("data/newobs.csv")
obs <- obs[which(obs[,4] != 'M'),]
obs <- obs[, c(3,2,4)]
top = 51 # north lat
left = -127 # west long
right = -65 # east long
bottom =  23 # south lat
#convert factors to numeric and as well as converting knots to m/s
obs[,3] <- as.numeric(as.character(obs[,3])) * 0.514444

obs.index <- which(obs[,1] > right | obs[,1] < left | obs[,2] > top | obs[,2] < bottom)
obs <- obs[-obs.index,]
z.gls <- obs[,3]
quilt.plot(obs[,1], obs[,2], obs[,3])
US(add=TRUE)

obs.lon.lat <- cbind(obs[,1], obs[,2])
temp1 <- mapproject(x=obs.lon.lat[,1],y=obs.lon.lat[,2],projection="mollweide")
obs.lon.lat <- cbind(temp1$x,temp1$y)
obs.mod <- lm(obs[,3] ~ obs.lon.lat[,1] + obs.lon.lat[,2]) # does both lat/lon
summary(obs.mod)

obs.z <- obs.mod$residuals
quilt.plot(obs.lon.lat, obs.z)

breaks <- seq(0,max(rdist(obs.lon.lat)),length.out=50)

v.bc <- variog(coords=obs.lon.lat,data=obs.z,estimator.type="classical",breaks=breaks,
               bin.cloud=TRUE) # binned

plot(v.bc,main="binned",pch=19,cex=0.5)
#seems that exp is the best
vfit <- variofit(v.bc,cov.model="exponential",weights="cressie")
lines(vfit)

a = vfit$cov.pars[2]
sigma = sqrt(vfit$cov.pars[1])
tausq = vfit$nugget

## Matrices for kriging
dist0.mat <- rdist(obs.lon.lat,pred.grd)
Sigma0 <- sigma^2 * exp(-dist0.mat/a)

dist.mat <- rdist(obs.lon.lat)
Sigma <- sigma^2 * exp(-dist.mat/a) + tausq


#Beta GLS

Sigma <- sigma^2 * exp(-rdist(obs.lon.lat)/a)
diag(Sigma) <- diag(Sigma) + tausq
Sigma.i <- solve(Sigma)
X <- cbind(1,obs.lon.lat) # design matrix of regression
beta.gls <- solve(t(X) %*% Sigma.i %*% X) %*% t(X) %*% Sigma.i %*% z.gls

v.bc.gls <- variog(coords=obs.lon.lat,data=z.gls,estimator.type="classical",breaks=breaks,
                   bin.cloud=TRUE) # binned

plot(v.bc.gls,main="binned",pch=19,cex=0.5)
#seems that exp is the best
vfit.gls <- variofit(v.bc.gls,cov.model="exponential",weights="cressie")
lines(vfit.gls)

a = vfit.gls$cov.pars[2]
sigma = sqrt(vfit.gls$cov.pars[1])
tausq = vfit.gls$nugget

## Matrices for kriging
dist0.mat <- rdist(obs.lon.lat,pred.grd)
Sigma0 <- sigma^2 * exp(-dist0.mat/a)

dist.mat <- rdist(obs.lon.lat)
Sigma <- sigma^2 * exp(-dist.mat/a) + tausq

## Ordinary kriging predictor
Sigma.i <- solve(Sigma)
ones <- rep(1,length(z.gls))
ones.mat <- matrix(1,nr=dim(Sigma0)[1],nc=dim(Sigma0)[2])
coef <- c(1 - (t(ones) %*% Sigma.i %*% Sigma0)) / c(t(ones) %*% Sigma.i %*% ones)
obs.z.hat.ok <- ( t(Sigma0) + coef * t(ones.mat) ) %*% Sigma.i %*% z.gls

#Let's take a look at the kriged data
c <- obs.mod$coefficients
obs.pred <- c[1] + c[2]*pred.grd[,1] + c[3]*pred.grd[,2] + obs.z.hat.ok
quilt.plot(pred.grd[,1], pred.grd[,2], obs.pred, nx = resolution, ny = resolution)


quilt.plot(new.loc, obs.pred, nx = resolution, ny = resolution)
US(add=TRUE)

par(mfrow=c(1,2))
quilt.plot(new.loc, pred, nx = resolution, ny = resolution)
US(add=TRUE)
quilt.plot(new.loc, obs.pred, nx = resolution, ny = resolution)
US(add=TRUE)

MSE = sum((pred-obs.pred)^2)