#
#
# Spectral Resampling
#
# # # # # # # # # # # #

require(hsdar)
require(rgdal)
require(raster)

setwd("~/data/ASO_CASI_SNOWX")

r = brick("CASI_2017_02_21_ort2_mosaic_specalb")
r
s = speclib("CASI_2017_02_21_ort2_mosaic_specalb") # doesnt work
s = speclib(r)

# obtain a list of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"\\.{3}"),
                         function(x) substr(x[2],1,nchar(x[2])-1))))

# ?
aso1 = speclib(tbl[,2:ncol(tbl)],aso_wv*1000,fwhm=3.5)

# Read in as Header
map = readGDAL("CASI_2017_02_21_ort2_mosaic_specalb")
summary(map)
image(map)


# spectral object
aso_spec = speclib(map, aso_wv)

#
####
# readin SNICAR data
tbl = read.csv("LookupTable2.csv")
tbl[,1] = tbl[,1]*1000
# create a spectral object
s = speclib(tbl[,2:ncol(tbl)],tbl[,1])
s = speclib(tbl[,2],tbl[,1])
w = wavelength(s)
fw = fwhm(s)


new_spec <- spectralResampling(s,aso_spec)



### WORKED!!!
#####
# SNICAR DATA
tbl = read.csv("LookupTable2.csv")
tbl[,1] = tbl[,1]*1000

# create a spectral object
s = speclib(tbl[,2:ncol(tbl)],tbl[,1]) ## NEED TO DO THIS FOR ALL WAVELENGTHS
grain_size = 100
s = speclib(tbl[,grain_size],tbl[,1])

#####
# ASO DATA
setwd("~/data/ASO_CASI_SNOWX")
r = brick("CASI_2017_02_21_ort2_mosaic_specalb")
# obtain a list of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
aso = data.frame("center"=aso_wv*1000,'fwhm'=3.5)

# SPECTRAL RESAMPLING
new_spec <- spectralResampling(s,aso)

# plot results
par(mfrow=c(3,1))
par(mar=c(4.5,4.5,2,2))
plot(new_spec, xlab = "Full Wavelength (nm)")
abline(v=950,col='red', lty=2)
legend('bottomleft',paste(grain_size,'microns'),bty='n',cex=1)
legend('topright',c("~950 nm","1040 nm"),col='red',bty='n',cex=0.8)


########
# SOLVE FOR INTEGRATED AREA FROM CONTINUUM
# TransformSpeclib
# Continuum is linear fit between. 0.95 microns to 1.04
#sp_tr <- transformSpeclib(new_spec,method='ch',out='ratio') # convex or segmented hull?
#sp_tr <- transformSpeclib(new_spec,method='ch',out='bd')
sp_tr <- transformSpeclib(new_spec[,c(62:72)],method='ch',out='bd')

# continuum line
s_line <- transformSpeclib(new_spec[,c(62:72)],method='ch',out='raw')
plot(s_line,ispec=1)
legend('topright',"Continuum",lty=2,bty='n',cex=1)

# extract values
d =data.frame('wv'= sp_tr@wavelength,'val'=sp_tr@spectra@spectra_ma[1,])

# Integrate
bd = sum(d$val[which(d$wv >= 950)])

# Add plot
# Observe bband depth from Continuum
plot(sp_tr)
legend('topleft',paste('Integration:',round(bd,4)),bty='n',cex=1)
# THIS is as band depth not as a percentage!!!

## Need to figure out how to add all grain sizes to initial snicar spectral object (s)
## how to subset data by wavelength so continuum is only between given points




