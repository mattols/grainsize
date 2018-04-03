#
# Grain Size Estimation 
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
r = brick("CASI/class_image/CASI_2017_02_21_ort2_mosaic_specalb")
# obtain a list of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
# create new dataframe with ASO wavelength properties
aso_df = data.frame("center"=aso_wv*1000,'fwhm'=3.5)

grain_size = 30
snicar_tbl = tbl_snicar
s = speclib(snicar_tbl[,grain_size-28],snicar_tbl[,1])

# SPECTRAL RESAMPLING (all at once)
s_mat = speclib(as.matrix(snicar_tbl[,2:ncol(snicar_tbl)]),
                snicar_tbl[,1])
snicar_resample <- spectralResampling(s_mat,aso_df)

# plot 30 and 1500 micro meter grain sizes
plot(snicar_resample[1],ylim=c(0.3,1))
lines((snicar_resample@wavelength),
      (snicar_resample[300,]@spectra@spectra_ma),
      col='forestgreen')
lines((snicar_resample@wavelength),
      (snicar_resample[1471,]@spectra@spectra_ma),
      col='blue')
abline(v=1030,col='red',lty=2)

#### Integrated Depth
# change dimentions first >> 
ss <- snicar_resample[,c(62:72)] #tail(aso_df[,1],11) ???
# Continuum

# Compute endpoints for each grain size
shortb = snicar_resample[c(1:1471),62]@spectra@spectra_ma
longb = snicar_resample[c(1:1471),72]@spectra@spectra_ma

# Plot start and end of feature for each grain size
plot(snicar_resample[c(1:1471),62]@spectra@spectra_ma,ylim=c(0.3,1),type='l')
lines(seq(1,1471),snicar_resample[c(1:1471),72]@spectra@spectra_ma,col='green')

# rise over run (11 bands)
slope = (longb - shortb)/10

#CONTINUUM NOT CORRECT !!!!!
cont = array(dim=c(11,length(shortb)))
for(i in 1:11){
  cont[i,] = (i-1)*slope + shortb
}
#test
plot(snicar_resample[1000,c(62:72)],type='l')
lines(snicar_resample@wavelength[62:72],cont[,1000],col='red',lty=2)

# Create integrand
integrand = array(dim=c(11,length(shortb)))
for (i in 1:11){
  integrand[i,] <- (cont[i,] - snicar_resample[c(1:1471),i+61]@spectra@spectra_ma )/cont[i,]
} 


# END RESULT == LOOKUP table



