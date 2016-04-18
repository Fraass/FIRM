### Setup file for Fraass & Lowery Foraminiferal Isotope Reproducibility Model (FIRM) ###
  # This file needs to be run prior to running the model.
  # This file will install the required package (ncdf) to download the Levitus dataset, the raw data
  # that the model runs on.
  # It also does a little big of transformation of that dataset, simply to be used more simply in the model.
  # This _does_ need to have access to the internet to download the files which are 
  # on the NOAA FTP site. 
  # After this file is run, FIRM (v.0.93).R should run fine.




## Downloading and loading neccessary packages and files
  install.packages('ncdf4') #download from whichever cran mirror you'd like
  library('ncdf4') #loading package for opening Levitus salinity & temperature profiles
  
  #temperature profiles
  temp.profile<-'ftp://ftp.cdc.noaa.gov/Datasets/nodc.woa98/temperat/seasonal/otemp.anal1deg.nc'
  download.file(temp.profile, "otemp.anal1deg.nc", method = "auto",
              quiet = FALSE, mode="wb", cacheOK = TRUE)
  #salinity profiles
  salt.profile<-'ftp://ftp.cdc.noaa.gov/Datasets/nodc.woa98/salinity/seasonal/salt.anal1deg.nc'
  download.file(salt.profile, "salt.anal1deg.nc", method = "auto",
              quiet = FALSE, mode="wb", cacheOK = TRUE)

## Transforming dataset
  #extracting temp profile
  temp.nc <- nc_open("otemp.anal1deg.nc")
  temp.nc

  t.lon<-ncvar_get(temp.nc,'lon')
  t.lat<-ncvar_get(temp.nc,'lat')
  t.depth<-ncvar_get(temp.nc,'level')
  t.time<-ncvar_get(temp.nc,'time')
  temp.temp<-ncvar_get(temp.nc,'otemp')

  #extracting salt profile
  salt.nc <- nc_open("salt.anal1deg.nc")
  salt.nc

  s.lon<-ncvar_get(salt.nc,'lon')
  s.lat<-ncvar_get(salt.nc,'lat')
  s.depth<-ncvar_get(salt.nc,'level')
  s.time<-ncvar_get(salt.nc,'time')
  salt.salt<-ncvar_get(salt.nc,'salt')

  #unloading temp and salt files 
  #(values used are still loaded, this just removes the large unprocessed datafile)
  rm(salt.nc)
  rm(temp.nc)

## Diagnostic Plot [ If a plot doesn't appear, there is an issue above ]
  #extract temperature/salinity profile (Again, diagnostic plot)
  lon<-155 #longitude (not longitude currently, but cell in levitus file)
  lat<-50 #latitude (not latitude currently, but cell in levitus file)
  time<-3 
  temp1<-temp.temp[lon,lat,,time]
  salt1<-salt.salt[lon,lat,,time]
  plot(temp1,
     t.depth,
     type='b',
     ylim=c(5500,0),
     pch=16,
     lwd=2,
     ylab="Depth (m)",
     xlab='Temperature (degC)'
  )
  lines(salt1/max(salt1)*400-385, #scale for salinity is wrong. Just checking if the data is there.
      s.depth,
      type='b',
      col='red',
      pch=16,
      lwd=2)

