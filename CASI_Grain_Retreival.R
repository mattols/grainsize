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

#### Integrated Depth
# change dimentions first >> 
ss <- snicar_resample[,c(62:72)] #tail(aso_df[,1],11) ???
# Continuum

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

#test ****
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
lookup2 <- rbind(c(0.08,NA),lookup) # anything too small (not likely snow)
lookup2 <- rbind(lookup2, c(0.71,NA)) # or anything too large
# arbitrary numbers chosen
head(lookup)

# END RESULT == LOOKUP table

##########################################################################################################
## # # # # ### ### ## ## #   # # # # # # # # #
# CREATE
### READ IN CASI DATA
r = brick("CASI/20170220/174008/CASI_2017_02_20_174008_atm_ort")
# obtain a list of wavelengths
aso_wv = as.numeric(unlist(lapply(strsplit(names(r),"bsq\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
# create new dataframe with ASO wavelength properties
aso_df = data.frame("center"=aso_wv*1000,'fwhm'=3.5)
aso_wv = aso_wv * 1000 # convert to nm
names(r) = aso_wv

# casi_mat = speclib(r,
#                 aso_wv)
####
#### Integrated Depth
# CONVERT to spectral (only essential bands)
# casi_mat = speclib(spectra=brick(stk),wavelength=aso_wv[c(62:72)], fwhm=rep(3.5,length(aso_wv[c(62:72)])))
# casi_mat2 <- casi_mat[,c(62:72)] #tail(aso_df[,1],11) ???

# DO BY ROW
#####################
strt = Sys.time()
gz <- stk[[1]]
gz[] <- NA

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



