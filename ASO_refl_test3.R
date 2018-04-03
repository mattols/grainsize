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

setwd("C:/Users/u1037042/Documents/MattOlson/projects/ASO")

### READ IN LOOKUP TABLE
lookup = read.csv("lookup.csv")
lookup2 <- rbind(c(0,NA),lookup)
lookup2 <- rbind(c(-0.001,0),lookup2)
# # READ AS ENVI FILE
filename = "CASI_2016_04_17_mosaic_atm_corr_ort"
r = read.ENVI(paste(filename, ".bsq", sep=""), headerfile=paste(filename, ".hdr", sep=""))
header = read.ENVI(filename, headerfile=paste(filename, ".hdr", sep=""))
### READ IN ASO DATA
r1 = brick("CASI_2016_04_17_mosaic_atm_corr_ort.bsq")
aso_wv = as.numeric(unlist(lapply(strsplit(names(r1),"\\.{3}"),
                                  function(x) substr(x[2],1,nchar(x[2])-1))))
aso_wv = aso_wv * 1000 # convert to nm
names(r1) = aso_wv

# PLOT IN COLOR (http://neondataskills.org/R/Multi-Band-Rasters-In-R/)
par(mfrow=c(1,1))
plotRGB(r1, r=18, g=14, b=3, stretch="lin") # stretch = "hist" 

# NDSI of image
# ndsi <- (r1)

# # Convert to Spectral
# s = speclib(spectra=r1,wavelength=aso_wv, fwhm=rep(3.5,length(aso_wv)))

# CONVERT to spectral (only essential bands)
stk = r1[[c(57:67)]]
s3 = speclib(spectra=brick(stk),wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))


## BAND DEPTH
strt = Sys.time()
sp_tr <- transformSpeclib(s3,method='ch',out='bd')
end = Sys.time() - strt 
print(end)


# IN CHUNKS (BREAK UP)
p <- crop(stk,extent(c(xmin=734000,xmax=735000,ymin=4365700,ymax=4366500)))

spec_p = speclib(spectra=as.matrix(p),wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))



###################################################################################
###################################################################################
###################################################################################

###################################################################################
### TESTING TESTING TESTING
# Test (FAKE DATA)
m = matrix(c(1,2,3,4,5, NA),ncol=3)
ml = matrix(c(2.0,2.1,2.2,2.3,2.4,2.5),ncol=3)
mm = matrix(c(m/ml,m/(ml+.1),m/(ml+.2),m/(ml+.3),m/(ml+.45),m/(ml+.55),m/(ml+.66),m/(ml+.7),
              m/(ml+.8),m/(ml+.9),m/(ml+1)),ncol=11)
mm[is.na(mm[])] <- 1
st = speclib(spectra=mm,wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
sp_t <- transformSpeclib(st,method='ch',out='bd')
s_line <- transformSpeclib(st,method='ch',out='raw')
# Observe results
par(mfrow=c(3,2))
for (i in 1:6){
  bd = round(sum(sp_t[i,]@spectra@spectra_ma),4)
  gr_size = lookup2[which(abs(lookup2[,1]-bd)==min(abs(lookup2[,1]-bd))),2]
  plot(s_line,ispec=i)
  legend('bottomleft',c(paste('Integration:',bd),paste('Grain Size:',gr_size)),bty='n',cex=1)
  legend('topright',as.character(i))
}


## MORE TESTING (ACTUAL DATA)
p <- crop(stk,extent(c(xmin=734000,xmax=734020,ymin=4365700,ymax=4365720)))
p_mat <- as.matrix(p)
st = speclib(spectra=p_mat,wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
st

##### One row at a time

# blank matrix
gz <- stk[[1]]
gz[] <- NA

# TIME TO DO ONE ROW
strt = Sys.time()
mat <- stk[1,]
mat[is.na(mat[])] <- 1
st = speclib(spectra=mat,wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
sp_t <- transformSpeclib(st,method='ch',out='bd')
#s_line <- transformSpeclib(st,method='ch',out='raw')
end = Sys.time() - strt 
print(end)

# place back
bd = rowSums(sp_t@spectra@spectra_ma)
gr_list = list()
for (i in 1:length(bd)){
  gr_size = lookup2[which(abs(lookup2[,1]-bd[i])==min(abs(lookup2[,1]-bd[i]))),2] # still doesn't work for negatives???
  gr_list = c(gr_list,gr_size)
}
#gr_size = lapply(bd, function(x) lookup2[which(abs(lookup2[,1]-bd[x])==min(abs(lookup2[,1]-bd[x]))),2]) # still doesn't work for negatives???

gz[1,] <- unlist(gr_list)



# ALL AT ONCE?
mat <- as.matrix(stk)
mat[is.na(mat[])] <- 0






# DO BY ROW
#####################
strt = Sys.time()
gz <- stk[[1]]
gz[] <- NA

for (j in 1:nrow(stk)){
  print(paste("...solving for row", j, "of", nrow(stk)))
  mat <- stk[j,]
  mat[is.na(mat[])] <- 1
  st = speclib(spectra=mat,wavelength=aso_wv[c(57:67)], fwhm=rep(3.5,length(aso_wv[c(57:67)])))
  sp_t <- transformSpeclib(st,method='ch',out='bd')
  # place back
  bd = rowSums(sp_t@spectra@spectra_ma)
  gr_list = list()
  for (i in 1:length(bd)){
    gr_size = lookup2[which(abs(lookup2[,1]-bd[i])==min(abs(lookup2[,1]-bd[i]))),2] # still doesn't work for negatives???
    gr_list = c(gr_list,gr_size)
  }  
  end = Sys.time() - strt 
  print(end)
  gz[j,] <- unlist(gr_list)
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
                   subject= paste("ASO Grain Size code completed"),
                   body = paste("Completed in", as.character(round(end,2)),'hours'),
                   smtp = list(host.name = "smtp.gmail.com", port = 587, user.name = sender, passwd = password, ssl = TRUE),
                   # smtp = list(host.name = "aspmx.l.google.com", port = 25),
                   authenticate = TRUE,
                   send = TRUE)


writeRaster(gz, filename = "gz.grd",
            format="raster")

plot(gz)
summary(gz)

names(gz) <- "grain_size"
r2 = gz
filename = "CASI_2016_04_17_snow_grain.grd"
r2 <- writeRaster(r2, filename, overwrite=TRUE)
hdr(r2, format="ENVI")

# filename = "CASI_2016_04_17_snow_grain"
# write.ENVI (gz, filename, interleave = "bsq")
# filename = "CASI_2016_04_17_snow_grain.bsq"
# grain_header <- writeRaster(gz, filename, overwrite=TRUE)
# hdr(grain_header, format="ENVI") 


 









# Threshold = 0.7? or lower 0.6

# NDVI
nir = 127       # 0.9 micro m
vis = 53      # 0.6
ndvi = (nir-vis)/(nir+vis)
ndvi


















