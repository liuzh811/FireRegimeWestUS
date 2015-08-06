##################################################################################
# Spatial synchrony analysis and Kennel density plot MTBS fires in the Western United States
##################################################################################

# Load libraries
library(ncf)          # Functions for spatial synchrony analysis
library(sp)           # classes and methods for spatial data
library(rgdal)        # functions for spatial data input/output
library(raster)
library("maptools")
library("maps")
library("rgeos") 
library("shapefiles")

setwd(".\\data")

##################################################################################
# Spatial synchrony analysis
##################################################################################
# Read in MTBS fire points 
firep.sp <- readOGR(dsn = ".\\data", layer = "MTBSJoshwest3Final")
fishnet.nm = c("WestGrid_20","WestGrid_50","WestGrid_100","WestGrid_150")

for(i in 1:4){
# read into fishnet
wgrid100.sp <- readOGR(dsn = ".\\data", layer = fishnet.nm[i])
wgrid100.sp@data$Id = 1:nrow(wgrid100.sp) #add by zhihua 1/4/14

# Map the fire points and grid cells
plot(firep.sp, pch=1, cex=0.1, col="red")
plot(wgrid100.sp, add=TRUE)

# Identify the grid cell associated with each fire point
gridindex.df <- over(firep.sp, wgrid100.sp)

# Add the grid cell identifier to the MTBS fire dataset
firep2.df <- data.frame(firep.sp@data, gridindex.df)

# Spatial synchrony analysis begins here
# for fire frequency
attach(firep2.df)
ffsum.df <- aggregate(AreaHa, by=list(gridcell=Id, year=YEAR_), FUN="length")
detach(firep2.df)
ffsum2.df <- with(ffsum.df, ffsum.df[year > 0,])
ffsum3.df <- reshape(ffsum2.df, timevar="year", idvar="gridcell", direction="wide")
ffsum3.df[is.na(ffsum3.df)] <- 0
gridff.sp <- merge(wgrid100.sp, ffsum3.df, by.x="Id", by.y="gridcell", all.x=FALSE)
plot(gridff.sp, col="blue")
gridcoords <- coordinates(gridff.sp)
ffdata.df <- gridff.sp@data[,2:28]
ff.sncf <- Sncf(gridcoords[,1], gridcoords[,2], ffdata.df, na.rm=T, resamp=1000)

# for median fire size
attach(firep2.df)
fssum.df <- aggregate(AreaHa, by=list(gridcell=Id, year=YEAR_), FUN="mean")
detach(firep2.df)
fssum2.df <- with(fssum.df, fssum.df[year > 0,])
fssum3.df <- reshape(fssum2.df, timevar="year", idvar="gridcell", direction="wide")
fssum3.df[is.na(fssum3.df)] <- 0
gridfs.sp <- merge(wgrid100.sp, fssum3.df, by.x="Id", by.y="gridcell", all.x=FALSE)
plot(gridfs.sp, col="blue")
gridcoords <- coordinates(gridfs.sp)
fsdata.df <- gridfs.sp@data[,2:28]
fs.sncf <- Sncf(gridcoords[,1], gridcoords[,2], fsdata.df, na.rm=T, resamp=1000)

# for fire severity
attach(firep2.df)
fireper.df <- aggregate(FSHPRO, by=list(gridcell=Id, year=YEAR_), FUN="median")
detach(firep2.df)
fireper3.df <- reshape(fireper.df, timevar="year", idvar="gridcell", direction="wide")
fireper3.df[is.na(fireper3.df)] <- 0.1
gridper.sp <- merge(wgrid100.sp, fireper3.df, by.x="Id", by.y="gridcell", all.x=FALSE)
plot(gridper.sp, col="blue")
gridcoords <- coordinates(gridper.sp)
perdata.df <- gridper.sp@data[,2:28]

per.sncf <- Sncf(gridcoords[,1], gridcoords[,2], perdata.df, na.rm=T, resamp=1000)

# plot results
output.figure.names = paste(".data\\spline_correlogram_", fishnet.nm[i], ".png", sep = "")
png(output.figure.names, height = 2000, width = 1200, res = 300, units = "px")
par(mfrow=c(3,1),mar=c(0,4,0,0),oma=c(1,1,0,0))

##plot the results for fire frequency 
matplot(as.numeric(ff.sncf$boot$boot.summary$predicted$x)/1000, as.numeric(ff.sncf$boot$boot.summary$predicted$y[2,]), cex.axis = 1.5, cex.lab=1.5, ylab = "", xlab = "", xaxt='n', type = "l",  ylim = c(-0.5, 0.5)) ##5% percentile
matplot(as.numeric(ff.sncf$boot$boot.summary$predicted$x/1000), as.numeric(ff.sncf$boot$boot.summary$predicted$y[10,]), type = "l", add = TRUE) ## 95% percentile

polygon(c(as.numeric(ff.sncf$boot$boot.summary$predicted$x/1000), rev(as.numeric(ff.sncf$boot$boot.summary$predicted$x/1000))), 
		c(as.numeric(ff.sncf$boot$boot.summary$predicted$y[10,]), rev(as.numeric(ff.sncf$boot$boot.summary$predicted$y[2,]))),
        col = "grey90", border = NA)
matplot(as.numeric(ff.sncf$boot$boot.summary$predicted$x/1000), as.numeric(ff.sncf$boot$boot.summary$predicted$y[6,]), type = "l", add = TRUE) ##mean 
axis(side = 1,  lwd.ticks=1.5, labels = NA) #add x axis ticks
abline(h = 0, lty = 2)
text(100, 0.45, "a)", cex = 1.5)

##plot the results for mean fire size 
matplot(as.numeric(fs.sncf$boot$boot.summary$predicted$x/1000), as.numeric(fs.sncf$boot$boot.summary$predicted$y[2,]), cex.axis = 1.5, cex.lab=1.5,ylab = "", xlab = "",xaxt='n', type = "l", xaxt='n', ylim = c(-0.5, 0.5)) ##5% percentile
matplot(as.numeric(fs.sncf$boot$boot.summary$predicted$x/1000), as.numeric(fs.sncf$boot$boot.summary$predicted$y[10,]), type = "l", add = TRUE) ## 95% percentile

polygon(c(as.numeric(fs.sncf$boot$boot.summary$predicted$x/1000), rev(as.numeric(fs.sncf$boot$boot.summary$predicted$x/1000))), 
		c(as.numeric(fs.sncf$boot$boot.summary$predicted$y[10,]), rev(as.numeric(fs.sncf$boot$boot.summary$predicted$y[2,]))),
        col = "grey90", border = NA)
matplot(as.numeric(fs.sncf$boot$boot.summary$predicted$x/1000), as.numeric(fs.sncf$boot$boot.summary$predicted$y[6,]), type = "l", add = TRUE) ##mean 
axis(side = 1,  lwd.ticks=1.5, labels = NA) #add x axis ticks
abline(h = 0, lty = 2)
text(100, 0.45, "b)", cex = 1.5)

##plot the results for high percent high severiy
par(mar=c(4, 4, 0, 0))
matplot(as.numeric(per.sncf$boot$boot.summary$predicted$x/1000), as.numeric(per.sncf$boot$boot.summary$predicted$y[2,]), cex.axis = 1.5,cex.lab=1.5,ylab = "", xlab = "",type = "l",  ylim = c(-0.5, 0.5)) ##5% percentile
matplot(as.numeric(per.sncf$boot$boot.summary$predicted$x/1000), as.numeric(per.sncf$boot$boot.summary$predicted$y[10,]), type = "l", add = TRUE) ## 95% percentile

polygon(c(as.numeric(per.sncf$boot$boot.summary$predicted$x/1000), rev(as.numeric(per.sncf$boot$boot.summary$predicted$x/1000))), 
		c(as.numeric(per.sncf$boot$boot.summary$predicted$y[10,]), rev(as.numeric(per.sncf$boot$boot.summary$predicted$y[2,]))),
        col = "grey90", border = NA)
matplot(as.numeric(per.sncf$boot$boot.summary$predicted$x)/1000, as.numeric(per.sncf$boot$boot.summary$predicted$y[6,]), type = "l", add = TRUE) ##mean 
axis(side = 1,  lwd.ticks=1.5, labels = NA, tcl = -0.2) #add x axis ticks
abline(h = 0, lty = 2)
text(100, 0.45, "c)", cex = 1.5)

mtext("Geographic distance (km)", side = 1, cex = 1.5,outer=TRUE,padj = -0.5,adj = 0.5) #http://stat.ethz.ch/R-manual/R-devel/library/graphics/html/mtext.html
mtext("Covariance", side = 2, cex = 1.5,outer=TRUE,padj = 0.5,adj = 0.5) 

dev.off()

print(paste("Finish Calculating", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

##################################################################################
# Smoothed Maps
##################################################################################
library(spatstat)
# read into boundary
b = readShapeSpatial(".data\\Bailey_west_dis.shp")
b1 = as(b, "SpatialPolygons")
b2 = as(b1, "owin")
b2.im = as.im(b2)

firep.sp = firep.sp[b,]

temp.ppp <- as(firep.sp, "ppp")
temp.ppp <- temp.ppp[b2]

# ppp object with no marks
firep.ppp <- unmark(temp.ppp)

# Marked ppp object - fire size
firesize.ppp <- firep.ppp
marks(firesize.ppp) <- log10(firep.sp$AreaHa)

# Marked ppp object - fire severity
firesev.ppp <- firep.ppp
marks(firesev.ppp) <- log10(firep.sp$FSHPRO + 1)

# Density plot 
firep.den = density(firep.ppp, sigma=50000)
firep.den2 = firep.den
firep.den2$v = firep.den$v*10000000000/27  # number of fire per million ha per year

# Smoothed fire size plot
firesize.den = markmean(firesize.ppp, sigma=100000)
firesize.den2 = firesize.den
firesize.den2$v = 10^(firesize.den$v)

# Smoothed fire severity plot
firesev.den = markmean(firesev.ppp, weights = firep.sp$AreaHa, sigma=50000) #weighted by fire size, 50000 seems to be an optimal radius
firesev.den2 = firesev.den
firesev.den2$v = 10^(firesev.den$v)

#plot results
png(".data\\SmoothMap.png",height = 3000, width = 2500, res = 300, units = "px")
par(mfrow=c(2,2),mar=c(0,0,1,1))

plot(firep.den2,main = "a): Fire Occurrence (# Fires/Million ha*year)", terrain.colors(10),
     ribside="right",ribsep=0.03,ribwid=0.05,ribn=2048,ribscale=1.5,ribargs=list(cex.axis=1.5))
	 
firesize.den3= firesize.den2
firesize.den3$v= log10(firesize.den2$v)
plot(firesize.den2,main = "b): Fire Size (ha)", terrain.colors(10),
     ribside="right",ribsep=0.03,ribwid=0.05,ribn=2048,ribscale=1.5,ribargs=list(cex.axis=1.5))

plot(firesev.den2, main = "c): Percent High Severity (%)", terrain.colors(10),
     ribside="right",ribsep=0.03,ribwid=0.05,ribn=2048,ribscale=1.5,ribargs=list(cex.axis=1.5))

dev.off()

