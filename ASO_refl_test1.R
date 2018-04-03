#
# 
# Utilizing SNICAR lookup table to determine grain size feature
#
# Matt Olson - University of Utah 03/2018
# # # # # # # # # # # #

require(hsdar)
require(rgdal)
require(raster)
require(caTools)

setwd("~/data/ASO")

### READ IN LOOKUP TABLE
lookup = read.csv("lookup.csv")

# READ AS ENVI FILE
filename = "CASI_2016_04_17_mosaic_atm_corr_ort"
r = read.ENVI(paste(filename, ".bsq", sep=""), headerfile=paste(filename, ".hdr", sep=""))

### READ IN ASO DATA
r1 = brick("CASI_2016_04_17_mosaic_atm_corr_ort.bsq")
aso_wv = as.numeric(unlist(lapply(strsplit(names(r1),"\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
aso_wv = aso_wv * 1000 # convert to nm
names(r1) = aso_wv

# Convert to Spectral
s = speclib(spectra=r1,wavelength=aso_wv, fwhm=rep(3.5,length(aso_wv)))

s2 = speclib(spectra=r1[5500005][c(57:67)],wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))

stk = r1[[c(57:67)]]
s3 = speclib(spectra=stk,wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))

# Function for single grain size
grain_lookup <- function(aso_spec, lookup_table){
  ## SOLVE FOR INTEGRATED AREA FROM CONTINUUM
  # Continuum is linear fit between. 0.95 microns to 1.04
  sp_tr <- transformSpeclib(aso_spec[,c(57:67)],method='ch',out='bd')
  
  # continuum line
  s_line <- transformSpeclib(aso_spec[,c(57:67)],method='ch',out='raw')
  
  ## PLOT
  par(mfrow=c(2,1))
  par(mar=c(4.5,4.5,2,2))
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
  
  # Use lookup table with calculated band depth
  ## Needs improvement!
  gr_size = lookup_table[which(abs(lookup_table[,1]-bd)==min(abs(lookup_table[,1]-bd))),2]
  
  return(gr_size)
}

grain_lookup(s,lookup)



#### DIFFERENT
grain_lookup <- function(rast_brick, lookup_table){
  grain_list = list()
  for (i in 1:ncell(rast_brick)){
    if (is.na(rast_brick[i][c(57:67)])[1]){
      print(paste(".> skip",i,'of', ncell(rast_brick)))
      grain_list = c(grain_list, NA)
    } else{
      print(paste("..calculating spectra for cell", i, 'of', ncell(rast_brick)))
      
      spec_range = speclib(rast_brick[i][c(57:67)],aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
      sp_tr <- transformSpeclib(spec_range,
                                method='ch',out='bd')
      
      # Calculate band depth
      bd = sum(sp_tr@spectra@spectra_ma[1,])
      
      # Lookup table
      gr_size = lookup_table[which(abs(lookup_table[,1]-bd)==min(abs(lookup_table[,1]-bd))),2]
      grain_list = c(grain_list, gr_size)
    }
  }
  return(grain_list)
}

grain_lookup(r1,lookup)


####################################
### ALL AT ONCE



stk = r1[[c(57:67)]]
s3 = speclib(spectra=brick(stk),wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
strt = Sys.time()
complete_stk <- stk[complete.cases(getValues(stk[[1]]))]
end = Sys.time() - strt 
print(end)
s4 = speclib(spectra=complete_stk,wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))

strt = Sys.time()
sp_tr <- transformSpeclib(s4[c(105050:105054),],method='ch',out='bd')
end = Sys.time() - strt 
print(end)

s_line <- transformSpeclib(s4[c(105050:105054),],method='ch',out='raw')

par(mfrow=c(4,1))
for (i in 1:4){
  
  bd = round(sum(sp_tr[i,]@spectra@spectra_ma),4)
  gr_size = lookup[which(abs(lookup[,1]-bd)==min(abs(lookup[,1]-bd))),2]
  plot(s_line,ispec=i)
  legend('bottomleft',c(paste('Integration:',bd),paste('Grain Size:',gr_size)),bty='n',cex=1)

}


## Not run:
supported_functions <- hsdar_parallel()
supported_functions


## Example for Linux and other systems where doMC is available
## Load library
library(doMC)
## Register number of workers
registerDoMC(3)
## Transform speclib using 3 cores
bd <- transformSpeclib(spectral_data)

## Close the cluster (important to get rid of processes)
closeCluster(cl)



par(mfrow=c(1,1))
plot(r1[[2]])

# Crop to smaller extent
p <- crop(stk,extent(c(xmin=734000,xmax=736000,ymin=4365700,ymax=4368000)))
p <- crop(stk,extent(c(xmin=734000,xmax=735000,ymin=4365700,ymax=4366500)))
p <- crop(stk,extent(c(xmin=734000,xmax=734300,ymin=4365700,ymax=4366000)))
p <- crop(stk,extent(c(xmin=734000,xmax=734100,ymin=4365700,ymax=4365800)))

spec_p = speclib(spectra=as.matrix(p),wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))

# Single
strt = Sys.time()
# p <- crop(stk,extent(c(xmin=734000,xmax=736000, # 35 mins
#                        ymin=4365700,ymax=4367000)))
p <- crop(stk,extent(c(xmin=734000,xmax=736000, #16 minutes
                       ymin=4365700,ymax=4366600)))
p
spec_p = speclib(spectra=as.matrix(p),wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
sp_tr <- transformSpeclib(spec_p,method='ch',out='bd')
dim(sp_tr@spectra@spectra_ma)
#Integrate
band_depth = rowSums(sp_tr@spectra@spectra_ma)
print(paste('...for',ncell(p),'pixels:'))
end = Sys.time() - strt 
print(end)

pix = c(1122,27722,88911,199800,289044)
mins = c((0.61/60),(31.27/60),3.95,16.74,35.43)
plot(pix, mins,
     xlab="pixel number",ylab='minutes to complete')

# determine the trend and extrapolate to 8 million
exponential.model <- lm(log(mins)~ pix)
exponential.model <- lm(mins ~ pix)
pixvalues <- seq(0, 8000000, 1)
Counts.exponential2 <- exp(predict(exponential.model,list(pix=pixvalues)))
plot(pix, mins,pch=16,xlim=c(0,800000),ylim=c(0,100))
lines(pixvalues, Counts.exponential2,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")
(tail(Counts.exponential2,1)/60)/24 #this is how many days it will take

# PARALLEL
strt = Sys.time()
## Not run:
supported_functions <- hsdar_parallel()
supported_functions


## Example for Linux and other systems where doMC is available
## Load library
library(doMC)
## Register number of workers
registerDoMC(3)

# Perform transformations
sp_tr <- transformSpeclib(spec_p,method='ch',out='bd')

## Close the cluster (important to get rid of processes)
closeCluster(cl)

end = Sys.time() - strt 
print(end)

# improves process from 31 sec to 27 sec
#
