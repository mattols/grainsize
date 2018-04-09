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

### DATA

# READ IN SNICAR DATA
tbl_snicar = read.csv("tables/LookupTable2.csv")
tbl_snicar[,1] = tbl_snicar[,1]*1000

# READ IN CASI DATA
r = brick("CASI/20170220/174008/CASI_2017_02_20_174008_atm_ort")
# List of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"bsq\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
# New dataframe with ASO wavelength properties
aso_df = data.frame("center"=aso_wv*1000,'fwhm'=3.5)
aso_wv = aso_wv * 1000 # convert to nm
names(r) = aso_wv # rename

# select bands used for integration
stk = r[[c(62:72)]]

### LOOKUP TABLE ###

# SPECTRAL RESAMPLING (all at once)
s_mat = speclib(as.matrix(tbl_snicar[,2:ncol(tbl_snicar)]),
                tbl_snicar[,1])
snicar_resample <- spectralResampling(s_mat,aso_df)

# plot 30 and 1500 micro meter grain sizes ****
plot(snicar_resample[1],ylim=c(0.3,1))
lines((snicar_resample@wavelength),
      (snicar_resample[300,]@spectra@spectra_ma),
      col='forestgreen')
lines((snicar_resample@wavelength),
      (snicar_resample[1471,]@spectra@spectra_ma),
      col='blue')
abline(v=1030,col='red',lty=2)

#### Integrated Depth (NOT USED)
ss <- snicar_resample[,c(62:72)] #tail(aso_df[,1],11) ???

# Compute endpoints for each grain size
shortb = snicar_resample[c(1:1471),62]@spectra@spectra_ma
longb = snicar_resample[c(1:1471),72]@spectra@spectra_ma

# Plot start and end of feature for each grain size ****
plot(snicar_resample[c(1:1471),62]@spectra@spectra_ma,ylim=c(0.3,1),type='l')
lines(seq(1,1471),snicar_resample[c(1:1471),72]@spectra@spectra_ma,col='green')

# rise over run (11 bands)
slope = (longb - shortb)/10

#CONTINUUM 
cont = array(dim=c(11,length(shortb)))
for(i in 1:11){
  cont[i,] = (i-1)*slope + shortb
}

#test ****PLOT
plot(snicar_resample[200,c(62:72)],type='l')
lines(snicar_resample@wavelength[62:72],cont[,200],col='red',lty=2)
integrand[,200]

# Create integrand
integrand = array(dim=c(11,length(shortb)))
for (i in 1:11){
  integrand[i,] <- (cont[i,] - snicar_resample[c(1:1471),i+61]@spectra@spectra_ma )/cont[i,]
} 

# Generate lookup table
integral = colSums(integrand)
grain_size = seq(30,1500)
lookup = data.frame(integral,grain_size)
lookup2 <- rbind(c(0.09,NA),lookup) # anything too small (not likely snow)
lookup2 <- rbind(lookup2, c(0.7,NA)) # or anything too large
# arbitrary numbers chosen
head(lookup2)

# END RESULT == LOOKUP table

##########################################################################################################
## # # # # ### ### ## ## #   # # # # # # # # #
# CREATE

### casi_mat = speclib(r,
#                 aso_wv)
####
#### Integrated Depth
# CONVERT to spectral (only essential bands)
#### casi_mat2 <- casi_mat[,c(62:72)] #tail(aso_df[,1],11) ???
# casi_mat = speclib(spectra=brick(stk),wavelength=aso_wv[c(62:72)], fwhm=rep(3.5,length(aso_wv[c(62:72)])))

# DO BY ROW
#####################
strt = Sys.time()
gz <- raster(r)

for (j in 1:nrow(stk)){
  print(paste("...solving for column", j, "of", ncol(stk)))
  cl <- stk[,j] # specify row
  # CONTINUUM 
  shortb <- stk[,j][,1]
  longb <- stk[,j][,11]
  slope = (longb - shortb)/10
  cont = array(dim=c(11,length(shortb)))
  for(i in 1:11){
    cont[i,] = (i-1)*slope + shortb
  }
  # Create integrand
  integrand = array(dim=c(11,length(shortb)))
  for (i in 1:11){
    integrand[i,] <- (cont[i,] - cl[,i] )/cont[i,]
  } 
  # use lookup table
  integral = colSums(integrand)
  gr_list = list()
  for (i in 1:length(integral)){
    gr_size = lookup2[which(abs(lookup2[,1]-integral[i])==min(abs(lookup2[,1]-integral[i]))),2] # still doesn't work for negatives???
    gr_list = c(gr_list,gr_size)
  }  
  end = Sys.time() - strt 
  print(end)
  gz[,j] <- unlist(gr_list)
}
end = Sys.time() - strt 
print(end)

require(mailR)
sender <- "robot.slave7@gmail.com"
password = 'iamarobotslave'
user.name = 'robot.slave7@gmail.com'

recipients <- c("olsonmeu@gmail.com") # Replace with one or more valid addresses
email <- send.mail(from = sender,
                   to = recipients,
                   subject= paste("CASI code complete"),
                   body = paste("Completed in", as.character(round(end,2)),'hours'),
                   smtp = list(host.name = "smtp.gmail.com", port = 587, user.name = sender, passwd = password, ssl = TRUE),
                   # smtp = list(host.name = "aspmx.l.google.com", port = 25),
                   authenticate = TRUE,
                   send = TRUE)


writeRaster(gz, filename = "gz.grd",
            format="raster")
names(gz) <- "grain_size"
r2 = gz
filename = "CASI_2017_02_20_174008_snow_grain.grd"
r2 <- writeRaster(r2, filename, overwrite=TRUE)
hdr(r2, format="ENVI")

# read in mask (snow == 1)
mask1 = brick("CASI/20170220/174008/CASI_2017_02_20_174008_classified_ort")
gz_mask <- raster::mask(gz, mask1, filename="", inverse=FALSE,
     maskvalue=c(0,2,3,4))
r2 = gz_mask
filename = "CASI_2017_02_20_174008_snow_grain_mask.grd"
r2 <- writeRaster(r2, filename, overwrite=TRUE)
hdr(r2, format="ENVI")

# plot RGB and grain size (blue)
par(mfrow=c(1,2))
par(oma=c(0,0,0,0))
plotRGB(r, r=23, g=19, b=8, stretch="lin")
plot(gz_mask,col=blues9,add=T)
plotRGB(r, r=23, g=19, b=8, stretch="lin")

# plot RGB and grain size (rainbow)
par(mfrow=c(1,2))
par(oma=c(0,0,0,0))
plotRGB(r, r=23, g=19, b=8, stretch="lin")
plot(gz_mask,col=rev(rainbow(8,alpha=0.5)),add=T)
plotRGB(r, r=23, g=19, b=8, stretch="lin")

# Overlay
par(mfrow=c(1,1))
plotRGB(r, r=23, g=19, b=8, stretch="lin")
plot(gz_mask,col=rev(rainbow(8,alpha=0.5)),add=T)

# hist
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,2,0))
hist(gz_mask[[1]], col='deepskyblue2',
     main= "CASI_2017_02_20_174008",
     xlab="Grain size (micro m)")



# Show histogram distribution of grain size by aspect
# show grain size vs. elevation also









# stop() # NEXT IS IMPROVED METHOD
###
#
#
#
####
##########################################################
### ALL AT ONCE (VERY FAST)

# Transform list to raster
make_raster <-function(list,reference_raster){
  # Creates raster from a list 
  raster(
    matrix(list,nrow=nrow(reference_raster),
           ncol=ncol(reference_raster), byrow=TRUE),
    xmn=reference_raster@extent@xmin, xmx=reference_raster@extent@xmax,
    ymn=reference_raster@extent@ymin, ymx=reference_raster@extent@ymax, 
    crs=crs(reference_raster))
}

# Select bands in rasterbrick
stk = r[[c(63:72)]]

# Create grain size raster
gz = raster(r)
spec_vals <- stk[] 
n_bands <- dim(spec_vals)[2]

# Create continuum 
shortb <- spec_vals[,1]
longb <- spec_vals[,n_bands]
slope = (longb - shortb)/(n_bands-1)
cont = array(dim=c(n_bands, length(shortb)))
for(i in 1:n_bands){
  print(i)
  cont[i,] = (i-1)*slope + shortb
}

# Create integrand
integrand = array(dim=c(n_bands, length(shortb)))
for (i in 1:n_bands){
  integrand[i,] <- (cont[i,] - spec_vals[,i] )/cont[i,]
} 

# Lookup table
strt = Sys.time()
integrals = colSums(integrand)
DT <- data.table(lookup2)
setkey(DT,integral)
dt_size <- DT[J(integrals), roll='nearest']
gz <- make_raster(dt_size$grain_size,gz)
print(Sys.time() - strt)
##############################################################





#TEST shape of modeled and measures spectral profiles
# Snicar
for (i in seq(1,1471,30)){
  if (i == 1){
    plot(snicar_resample[i][,c(62:72)],ylim=c(0.3,0.93))
  } else{
    lines((snicar_resample@wavelength[c(62:72)]),
          (snicar_resample[i,]@spectra@spectra_ma[c(62:72)]),
          col='black',lty=2)
  }
}
# CASI???


