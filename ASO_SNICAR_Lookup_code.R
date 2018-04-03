#
# Calculating Band Depth of the Ice Absorption Feature
# Through Spectral Resampling and integrating between Continuum and
# Spectral Reflectance Values
#
# Matt Olson - University of Utah 03/2018
# # # # # # # # # # # #

require(hsdar)
require(rgdal)
require(raster)

setwd("~/data/ASO_CASI_SNOWX")

### READ IN SNICAR DATA
tbl_snicar = read.csv("LookupTable2.csv")
tbl_snicar[,1] = tbl_snicar[,1]*1000

### READ IN ASO DATA
r = brick("CASI_2017_02_21_ort2_mosaic_specalb")
# obtain a list of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
# create new dataframe with ASO wavelength properties
aso_df = data.frame("center"=aso_wv*1000,'fwhm'=3.5)

# Function for single grain size
grain_band_depth <- function(grain_size, snicar_tbl, aso_df){
  # create a spectral object from table
  s = speclib(snicar_tbl[,grain_size-28],snicar_tbl[,1])
  
  # SPECTRAL RESAMPLING
  new_spec <- spectralResampling(s,aso_df)
  
  # plot results
  par(mfrow=c(3,1))
  par(mar=c(4.5,4.5,2,2))
  plot(new_spec, xlab = "Full Wavelength (nm)")
  abline(v=950,col='red', lty=2)
  legend('bottomleft',paste(grain_size,'microns'),bty='n',cex=1)
  legend('topright',c("~950 nm","1040 nm"),col='red',bty='n',cex=0.8)
  
  ## SOLVE FOR INTEGRATED AREA FROM CONTINUUM
  # Continuum is linear fit between. 0.95 microns to 1.04
  sp_tr <- transformSpeclib(new_spec[,c(62:72)],method='ch',out='bd')
  
  # continuum line
  s_line <- transformSpeclib(new_spec[,c(62:72)],method='ch',out='raw')
  plot(s_line,ispec=1)
  legend('topright',"Continuum",lty=2,bty='n',cex=1)
  
  # extract values
  d = data.frame('wv'= sp_tr@wavelength,'val'=sp_tr@spectra@spectra_ma[1,])
  
  # Integrate
  bd = sum(d$val[which(d$wv >= 950)])
  
  # Add plot
  # Observe band depth from Continuum
  plot(sp_tr)
  legend('topleft',paste('Integration:',round(bd,4)),bty='n',cex=1)
  # THIS is as band depth not as a percentage!!!
  
  return(bd)
}

####
# create final Lookup table
look_df <- data.frame(matrix(nrow=length(seq(30,1500)),ncol=2))
gr_sizes = seq(30,1500)
for (i in 1:length(gr_sizes)){
  print(paste('...band depth for grain size', gr_sizes[i], 'of', tail(gr_sizes,1)))
  out = grain_band_depth(gr_sizes[i], tbl, aso)
  look_df[i,1] = out # band depth
  look_df[i,2] = gr_sizes[i]
}
names(look_df) = c('band_depth','grain_size')

# PLOT BD ACROSS GRAIN SIZE
par(mfrow=c(1,1))
plot(look_df$band_depth~look_df$grain_size,type='l',
     xlab='Grain size (microns)', ylab="Band Depth")

# SAVE FILE
write.csv(look_df,file="lookup.csv", row.names=FALSE)

# Negative values
look_df2 = look_df
look_df2[,1] = look_df2[,1]*-1
write.csv(look_df2,file="lookup2.csv", row.names=FALSE)
