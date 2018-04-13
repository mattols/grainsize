#
#
# Plots to test CASI Spectra and lookup methods
# Matt Olson - UofU - 04/12/2018
#
# # # # # # # #

require(hsdar)
require(rgdal)
require(raster)
require(caTools)

# LOAD REQUIRED DATA
setwd("~/data/ASO")

### READ IN SNICAR DATA
tbl_snicar = read.csv("tables/LookupTable2.csv")
tbl_snicar[,1] = tbl_snicar[,1]*1000

### READ IN CASI DATA
# SENATOR BECK
r_sb = brick("CASI/20170220/174008/CASI_2017_02_20_174008_atm_ort") #SBB
# obtain a list of wavelengths
aso_wv_72 = as.numeric(unlist(lapply(strsplit(names(r_sb),"bsq\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))

# SAGEHEN
r_sh = brick("CASI/sagehen_mosaic/CASI_2016_04_17_mosaic_atm_corr_ort.bsq")
# obtain a list of wavelengths
aso_wv_67 = aso_wv_72[c(6:72)]

# Choose REGION
r = r_sh # r_sb
aso_wv = aso_wv_67 # aso_wv_72

# create new dataframe with ASO wavelength properties
aso_df = data.frame("center"=aso_wv*1000,'fwhm'=3.5)
aso_wv = aso_wv * 1000 # convert to nm
names(r) = aso_wv

# SPECTRAL RESAMPLING (all at once)
s_mat = speclib(as.matrix(tbl_snicar[,2:ncol(tbl_snicar)]),
                tbl_snicar[,1])
snicar_resample <- spectralResampling(s_mat,aso_df)

##############################################################################
##############################################################################
# SMART PLOT
# SHOWS BEST ESTIMATE BASED ON INTEGRAL AND SLOPE FOR A GIVEN RANGE
# SOLVING FOR SINGLE CASI PIXEL VALUE
require(data.table)

# Funcitons
create_lookup <- function(snicar_resample,band_range){
  # Computes table based on range and snicar data
  shortb = snicar_resample[c(1:1471),band_range[1]]@spectra@spectra_ma
  longb = snicar_resample[c(1:1471),tail(band_range,1)]@spectra@spectra_ma
  slope = (longb - shortb)/(length(band_range) - 1)
  cont = array(dim=c(length(band_range),length(shortb)))
  for(i in 0:(length(band_range) - 1)){
    cont[i+1,] = (i)*slope + shortb
  }
  integrand = array(dim=c(length(band_range),length(shortb)))
  for (i in 1:length(band_range)){
    integrand[i,] <- (cont[i,] - snicar_resample[c(1:1471),band_range[i]]@spectra@spectra_ma )/cont[i,]
  } 
  
  # Generate lookup table
  integral = colSums(integrand)
  grain_size = seq(30,1500)
  lookup = data.frame(integral,slope, grain_size)
  # ADD NA VALUES FOR EXTREMES
  # lookup <- rbind(c(lookup[1,1] - abs(lookup[1,1] - lookup[2,1]),
  #                   lookup[1,2] - abs(lookup[1,2] - lookup[2,2]),NA),lookup) # anything too small (not likely snow)
  # lookup <- rbind(c(lookup[1471,1] + abs(lookup[1471,1] - lookup[1470,1]),
  #                   lookup[1471,2] + abs(lookup[1471,2] - lookup[1470,2]),NA),lookup) # or anything too large
  return(lookup)
}

lookup_values <- function(casi_integral,casi_slope,snicar_spectral,b_range){
  # look up the best match for integral and slope
  lookup <- create_lookup(snicar_spectral,b_range)
  DT <- data.table(lookup)
  setkey(DT,integral)
  grain_integral <- DT[J(casi_integral), roll='nearest']
  setkey(DT,slope)
  grain_slope <- DT[J(casi_slope), roll='nearest']
  return(c(grain_integral$grain_size,grain_slope$grain_size))
  # return(matrix(c(grain_integral$integral, grain_slope$slope,
  #                 grain_integral$grain_size,grain_slope$grain_size),ncol=2))
  # # returns a matrix row 1 is integral, row 2 is slope
  # # grain estimation for each (column 2)
}


##### # # # # # # # # # # # # #
# PLOT CASI SPECTRA VS MODELED

casi_test <- function(test_casi_full,band_range,pixel){
  par(mfrow=c(3,2))
  par(oma=c(3.5,3.5,1,3))
  par(mar=c(0.5,0.5,0,0))
  
  rng = band_range
  
  # plot CASI
  pix = pixel            # POSSIBLE VARIABLE -> ITERATE THROUGH
  y_lim <- range(test_casi_full[pix,]@spectra@spectra_ma)+ 
    (rep(diff(range(test_casi_full[pix,]@spectra@spectra_ma))/30,2)*c(-1,1))
  #y_lim[1] = 60 # FIX!
  y_lim2 <- range(test_casi_full[pix,rng]@spectra@spectra_ma) + 
    (rep(diff(range(test_casi_full[pix,rng]@spectra@spectra_ma))/10,2)*c(-1,1))
  plot(test_casi_full[pix,],ylim=y_lim,xaxt='n') # plot some pixel
  legend('bottom',c("       CASI",paste(as.character(aso_wv[c(rng[1],tail(rng,1))]),collapse='--')),cex=1.2,bty='n')
  polygon(c(aso_wv[rng[1]],aso_wv[tail(rng,1)]+10,aso_wv[tail(rng,1)]+10,aso_wv[rng[1]]),
          c(y_lim2[2],y_lim2[2],y_lim2[1],y_lim2[1]),col= rgb(255, 0, 0, 5, maxColorValue=255),border='firebrick')
  plot(test_casi_full[pix,rng],xaxt='n',yaxt='n')
  #points(aso_wv[rng],rep(y_lim2[2]-10,length(rng)))
  #legend('top',legend='band locations',cex=1,bty='n')
  axis(4)
  box(col='firebrick')
  
  # CASI CALCULATE SLOPE AND INTEGRATION
  shortb = test_casi_full[pix,rng[1]]@spectra@spectra_ma
  longb = test_casi_full[pix,tail(rng,1)]@spectra@spectra_ma
  casi_slope = (longb - shortb)/(length(rng)-1)
  cont <- c(casi_slope)*seq(0,length(rng)-1) + c(shortb)
  integrand <- (cont - test_casi_full[pix,rng]@spectra@spectra_ma )/cont
  casi_integral = sum(integrand)
  
  # Plot slope and Integral
  lines(test_casi_full@wavelength[rng],cont, col='firebrick',lty=2)
  legend('left',c(paste("Slope:",as.character(round(casi_slope,4))),
                  paste("Integral:",as.character(round(casi_integral,4)))),bty='n')
  
  
  results <- lookup_values(casi_integral = casi_integral,casi_slope = c(casi_slope),
                           snicar_spectral = snicar_resample,b_range = rng)
  
  # plot SNICAR for integration
  gs = results[1]-29 # some grain size
  y_lim <- range(snicar_resample[gs,]@spectra@spectra_ma)+ 
    (rep(diff(range(snicar_resample[gs,]@spectra@spectra_ma))/20,2)*c(-1,1))
  y_lim2 <- range(snicar_resample[gs,rng]@spectra@spectra_ma) + 
    (rep(diff(range(snicar_resample[gs,rng]@spectra@spectra_ma))/20,2)*c(-1,1))
  plot(snicar_resample[gs,],type='l',ylim=y_lim,xaxt='n')
  legend('bottomleft',c("SNICAR","Integration method",paste("grain size =",gs+29)),cex=1.2,bty='n')
  #legend('bottomleft',paste("grain size =",gs+29),bty='n')
  polygon(c(aso_wv[rng[1]],aso_wv[tail(rng,1)]+10,aso_wv[tail(rng,1)]+10,aso_wv[rng[1]]),
          c(y_lim2[2],y_lim2[2],y_lim2[1],y_lim2[1]),col= rgb(255, 0, 0, 5, maxColorValue=255),border='firebrick')
  plot(snicar_resample[gs,rng],type='l',yaxt='n',xaxt='n')
  #points(aso_wv[rng],rep(y_lim2[2],length(rng)))
  axis(4)
  box(col='firebrick')
  
  # SNICAR CALCULATE SLOPE AND INTEGRATION
  shortb = snicar_resample[gs,rng[1]]@spectra@spectra_ma
  longb = snicar_resample[gs,tail(rng,1)]@spectra@spectra_ma
  slope = (longb - shortb)/(length(rng)-1)
  cont <- c(slope)*seq(0,length(rng)-1) + c(shortb)
  integrand <- (cont - snicar_resample[gs,rng]@spectra@spectra_ma )/cont
  integral = sum(integrand)
  
  # Plot slope and Integral
  lines(snicar_resample@wavelength[rng],cont, col='firebrick',lty=2)
  legend('bottomleft',c(paste("CASI integral:",round(casi_integral,4)),
                        paste("closest match:",as.character(round(integral,4))),
                        paste("SNICAR grain estimation:",results[1])),bty='n')
  ### SLOPE
  # plot SNICAR for slope
  gs = results[2]-29 # some grain size
  y_lim <- range(snicar_resample[gs,]@spectra@spectra_ma)+ 
    (rep(diff(range(snicar_resample[gs,]@spectra@spectra_ma))/20,2)*c(-1,1))
  y_lim2 <- range(snicar_resample[gs,rng]@spectra@spectra_ma) + 
    (rep(diff(range(snicar_resample[gs,rng]@spectra@spectra_ma))/20,2)*c(-1,1))
  plot(snicar_resample[gs,],type='l',ylim=y_lim)
  legend('bottomleft',c("SNICAR","Slope method", paste("grain size =",gs+29)),cex=1.2,bty='n')
  #legend('bottomleft',paste("grain size =",gs+29),bty='n')
  polygon(c(aso_wv[rng[1]],aso_wv[tail(rng,1)]+10,aso_wv[tail(rng,1)]+10,aso_wv[rng[1]]),
          c(y_lim2[2],y_lim2[2],y_lim[1],y_lim[1]),col= rgb(255, 0, 0, 5, maxColorValue=255),border='firebrick')
  plot(snicar_resample[gs,rng],type='l',yaxt='n')
  points(aso_wv[rng],rep(y_lim2[2]-0.01,length(rng)))
  axis(4)
  box(col='firebrick')
  
  # SNICAR CALCULATE SLOPE AND INTEGRATION
  shortb = snicar_resample[gs,rng[1]]@spectra@spectra_ma
  longb = snicar_resample[gs,tail(rng,1)]@spectra@spectra_ma
  slope = (longb - shortb)/(length(rng)-1)
  cont <- c(slope)*seq(0,length(rng)-1) + c(shortb)
  integrand <- (cont - snicar_resample[gs,rng]@spectra@spectra_ma )/cont
  integral = sum(integrand)
  
  # Plot slope and Integral
  lines(snicar_resample@wavelength[rng],cont, col='firebrick',lty=2)
  legend('bottomleft',c(paste("CASI slope:",round(casi_slope,4)),
                        paste("closest match:",as.character(round(slope,4))),
                        paste("SNICAR grain estimation:",results[2])),bty='n')
  
  
  # par(mar=c(4.5,4.5,2.0,1))
  mtext(text="Wavelength (nm)",side=1,line=2,outer=TRUE)
  mtext(text="Reflectivity",side=2,line=2,outer=TRUE)
  par(oma=c(0,0,0,0))
  par(mar=c(4.5,4.5,2.0,1))
}


# ...
# RUN FROM HERE

# select other casi image
# SENATOR BECK
# ROI
roi_l <- crop(r, extent(c(xmin=260100, xmax=260800, ymin=4198600, ymax=4199600)))
roi_s <- crop(r, extent(c(xmin=260300, xmax=260370, ymin=4198890, ymax=4198950)))

# Choose REGION
r = r_sh # r_sb
aso_wv = aso_wv_67 # aso_wv_72

# SAGEHEN
# find ROI small
par(mfrow=c(1,1))
plot(r[[10]])
bb_l <- c(732500,733500,4367500,4369000)
abline(v=bb_l[1],lty=2) #xmin
abline(v=bb_l[2],lty=2) #xmax
abline(h=bb_l[3],lty=2) #ymin
abline(h=bb_l[4],lty=2) #ymax
bb <- c(733000,733300,4367900,4368100) #x_min,x_max,y_min,y_max
abline(v=bb[1]) #xmin
abline(v=bb[2]) #xmax
abline(h=bb[3]) #ymin
abline(h=bb[4]) #ymax

# ROI
roi_l <- crop(r, extent(c(xmin=bb_l[1], xmax=bb_l[2], ymin=bb_l[3], ymax=bb_l[4])))
roi_s <- crop(r, extent(c(xmin=bb[1], xmax=bb[2], ymin=bb[3], ymax=bb[4])))

# ASSIGN TEST ROI 
test_casi_full = speclib(spectra=roi_s[],wavelength=aso_wv,
                         fwhm=rep(3.5,length(aso_wv)))

par(mfrow=c(1,2))
plot(test_casi_full[3300,])
# rectangular bounding box
e = extent(roi_s)
p <- as(e, 'SpatialPolygons')
plotRGB(roi_l, r=23, g=19, b=8, stretch="lin")
plot(p, lwd=2,border='red',add=T)
# RESET AFTER RGB IMAGE
par(mar=par()$mar)

# PLOT ALL SPECTRA IN ROI_S
par(mfrow=c(1,2))
adj = 0
for(i in 1:(2*ncol(roi_s))){
  if (i==1){
    plot(test_casi_full[i+adj,],type='l',ylim=c(0,max(roi_s@data@max)),main=paste((2*ncol(roi_s)),'spectral pixels'))
  } else if((i %% 10) == 0){ # only plot every10nth value
    # do nothing
    lines(NA,NA)
  }else{
    lines(test_casi_full@wavelength, test_casi_full@spectra@spectra_ma[i+adj,c(1:length(aso_wv))],col='black',lty=2)
  }
}
e = extent(roi_s)
p <- as(e, 'SpatialPolygons')
plotRGB(roi_l, r=23, g=19, b=8, stretch="lin")
plot(p, lwd=2,border='red',add=T)
# RESET AFTER RGB IMAGE
par(mar=par()$mar)


# # # #  # #
## # # ## # #
# call function
range_of_interest = c(57:67) # SB 62:72 SH 57:67
paste(range(aso_wv[range_of_interest]),collapse="--")
casi_test(test_casi_full, band_range = range_of_interest, pixel = 10)

# randomly select a few pixels
pixels <- sample(1:ncell(roi_s), 10)
range_of_interest = c(56:65) # SB 62:72 SH 57:67
paste(range(aso_wv[range_of_interest]),collapse="--")
for (i in pixels){
  casi_test(test_casi_full, band_range = range_of_interest, pixel = i)
  print(paste("pixel number:", i, "for", match(i,pixels),"of", length(pixels)))
  title(paste("cellnumber:",i))
}
# end

nlookup <- create_lookup(snicar_resample,range_of_interest)
head(nlookup)
tail(nlookup)
nlookup[1420:1422,]

# could add pixel variable
##############################################################################
##############################################################################
##############################################################################






































# ...




















##############################################################################
par(mar=c(4.5,4.5,2.0,1))

#PLOT SPECTRA AT SMALLER TEST SITE
# Test Location (No NA's)
# rn_l <- crop(r, extent(c(xmin=260100, xmax=260800, ymin=4198500, ymax=4199900)))
# rn <- crop(r, extent(c(xmin=260500, xmax=260600, ymin=4199250, ymax=4199400)))
# sites
rn_l <- crop(r, extent(c(xmin=260100, xmax=260800, ymin=4198600, ymax=4199600)))
rn <- crop(r, extent(c(xmin=260515, xmax=260560, ymin=4199290, ymax=4199350)))
rn2 <- crop(r, extent(c(xmin=260300, xmax=260370, ymin=4198890, ymax=4198950)))
test_casi_full = speclib(spectra=rn2[],wavelength=aso_wv, fwhm=rep(3.5,length(aso_wv)))

par(mfrow=c(1,2))
for(i in 1:ncell(rn2)){
  if (i==1){
    plot(test_casi_full[i,],type='l',ylim=c(0,150),main=paste(ncell(rn2),'spectral pixels'))
  }else{
    lines(test_casi_full@wavelength, test_casi_full@spectra@spectra_ma[i,c(1:72)],col='black',lty=2)
  }
}

# rectangular bounding box
e = extent(rn2)
p <- as(e, 'SpatialPolygons')
plotRGB(rn_l, r=23, g=19, b=8, stretch="lin")
plot(p, lwd=2,border='red',add=T)

# RESET AFTER RGB IMAGE
par(mar=par()$mar)


##############################################################################
##############################################################################
##############################################################################
# old
# 2 X 2 plot CASI VS SNICAR
par(mfrow=c(2,2))
par(oma=c(3.5,3.5,1,3))
par(mar=c(0.5,0.5,0,0))

rng = c(59:72)

# plot CASI
pix = 10
y_lim = c(60,100)
plot(test_casi_full[pix,],ylim=y_lim+15,xaxt='n') # plot some pixel
legend('bottom',"CASI",cex=1.2,bty='n')
polygon(c(aso_wv[rng[1]],aso_wv[tail(rng,1)]+10,aso_wv[tail(rng,1)]+10,aso_wv[rng[1]]),
        c(110,110,75,75),col= rgb(255, 0, 0, 5, maxColorValue=255),border='firebrick')
plot(test_casi_full[pix,rng],xaxt='n',yaxt='n')
points(aso_wv[rng],rep(y_lim[2],length(rng)))
legend('top',legend='band locations',cex=1,bty='n')

axis(4)
box(col='firebrick')

# CASI CALCULATE SLOPE AND INTEGRATION
shortb = test_casi_full[pix,rng[1]]@spectra@spectra_ma
longb = test_casi_full[pix,tail(rng,1)]@spectra@spectra_ma
slope = (longb - shortb)/(length(rng)-1)
cont <- c(slope)*seq(0,length(rng)-1) + c(shortb)
integrand <- (cont - test_casi_full[pix,rng]@spectra@spectra_ma )/cont
integral = sum(integrand)

# Plot slope and Integral
lines(test_casi_full@wavelength[rng],cont, col='firebrick',lty=2)
legend('left',c(paste("Slope:",as.character(round(slope,4))),
                paste("Integral:",as.character(round(integral,4)))),bty='n')

# plot SNICAR
gs = 1471 # some grain size
y_lim <- c(min(snicar_resample[gs,]@spectra@spectra_ma)-0.1,max(snicar_resample[gs,]@spectra@spectra_ma)+0.1)
plot(snicar_resample[gs,],type='l',ylim=y_lim)
legend('bottom',c("SNICAR",paste("grain size =",gs+29)),cex=1.2,bty='n')
#legend('bottomleft',paste("grain size =",gs+29),bty='n')
polygon(c(aso_wv[rng[1]],aso_wv[tail(rng,1)]+10,aso_wv[tail(rng,1)]+10,aso_wv[rng[1]]),
        c(y_lim[2]-.2,y_lim[2]-.2,y_lim[1],y_lim[1]),col= rgb(255, 0, 0, 5, maxColorValue=255),border='firebrick')
plot(snicar_resample[gs,rng],type='l',yaxt='n')
points(aso_wv[rng],rep(y_lim[2],length(rng)))
axis(4)
box(col='firebrick')

# SNICAR CALCULATE SLOPE AND INTEGRATION
shortb = snicar_resample[gs,rng[1]]@spectra@spectra_ma
longb = snicar_resample[gs,tail(rng,1)]@spectra@spectra_ma
slope = (longb - shortb)/(length(rng)-1)
cont <- c(slope)*seq(0,length(rng)-1) + c(shortb)
integrand <- (cont - snicar_resample[gs,rng]@spectra@spectra_ma )/cont
integral = sum(integrand)

# Plot slope and Integral
lines(snicar_resample@wavelength[rng],cont, col='firebrick',lty=2)
legend('left',c(paste("Slope:",as.character(round(slope,4))),
                paste("Integral:",as.character(round(integral,4)))),bty='n')

# par(oma=c(0,0,0,0))
# par(mar=c(4.5,4.5,2.0,1))
mtext(text="Wavelength (nm)",side=1,line=2,outer=TRUE)
mtext(text="Reflectivity",side=2,line=2,outer=TRUE)







