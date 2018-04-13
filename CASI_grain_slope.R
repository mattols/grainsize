#
# Grain Size Estimation basedon slope
# Nolin and Dozier method developed for CASI
#
# Matt Olson - University of Utah 04/2018
# # # # # # # # # # # #

require(hsdar)
require(rgdal)
require(raster)
require(caTools)

setwd("~/data/ASO")

### READ IN SNICAR DATA
tbl_snicar = read.csv("tables/LookupTable2.csv")
tbl_snicar[,1] = tbl_snicar[,1]*1000

### READ IN CASI DATA
r = brick("CASI/20170220/174008/CASI_2017_02_20_174008_atm_ort")
# obtain a list of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"bsq\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
# create new dataframe with ASO wavelength properties
aso_df = data.frame("center"=aso_wv*1000,'fwhm'=3.5)
aso_wv = aso_wv * 1000 # convert to nm
names(r) = aso_wv
stk = r[[c(61:70)]] # or (59) 61 to 70

# SPECTRAL RESAMPLING (all at once)
s_mat = speclib(as.matrix(tbl_snicar[,2:ncol(tbl_snicar)]),
                tbl_snicar[,1])
snicar_resample <- spectralResampling(s_mat,aso_df)

# Compute endpoints for each grain size
shortb = snicar_resample[c(1:1471),61]@spectra@spectra_ma
longb = snicar_resample[c(1:1471),70]@spectra@spectra_ma

# rise over run (11 bands)
slope = (longb - shortb)/9 #11 or 9

###PLOTS
par(mfrow=c(3,5))
for(i in seq(100,1500,100)){
  contin <- slope[i]*seq(0,9,1) + shortb[i]
  plot(snicar_resample[i,c(61:70)],type='l')
  lines(snicar_resample@wavelength[61:70],contin,col='red',lty=2)
  legend('bottomleft',paste(as.character(round(slope[i],4))),bty='n')
  legend('topright',paste(as.character(i)),bty='n')
}

contin <- slope[200]*seq(0,9,1) + shortb[200]
plot(snicar_resample[200,c(61:70)],type='l')
lines(snicar_resample@wavelength[61:70],contin,col='red',lty=2)

# Generate lookup table
grain_size = seq(30,1500)
lookup = data.frame(slope,grain_size)
lookup2 <- rbind(c(-0.0062,NA),lookup) # anything too small (not likely snow)
lookup2 <- rbind(lookup2, c(-0.022287,NA)) # or anything too large
# arbitrary numbers chosen
head(lookup2)

# END RESULT == LOOKUP table

##########################################################
### ALL AT ONCE

# Create grain size raster
gz = raster(r)
spec_vals <- stk[] 
n_bands <- dim(spec_vals)[2]

# Create continuum 
shortb <- spec_vals[,1]
longb <- spec_vals[,n_bands]
slopes = (longb - shortb)/(n_bands-1)

# Lookup table (longest time) IMPROVE!!!!!
strt = Sys.time()
DT <- data.table(lookup2)
setkey(DT,slope)
dt_size <- DT[J(slopes), roll='nearest']
gz <- make_raster(dt_size$grain_size,gz)
print(Sys.time() - strt)

##############################################################

casi_mat = speclib(spectra=stk[],wavelength=aso_wv[c(61:70)], fwhm=rep(3.5,length(aso_wv[c(61:70)])))

stkt = rn[[c(61:70)]]
test_casi = speclib(spectra=stkt[],wavelength=aso_wv[c(61:70)], fwhm=rep(3.5,length(aso_wv[c(61:70)])))

par(mfrow=c(3,3))
for(i in seq(10800,11001)){
  contin <- slopes[i]*seq(0,9,1) + shortb[i]
  plot(test_casi[i,],type='l')
  lines(casi_mat@wavelength,contin,col='red',lty=2)
  legend('bottomleft',paste(as.character(round(slope[i],4))),bty='n')
  legend('topright',paste(as.character(i)),bty='n')
}


