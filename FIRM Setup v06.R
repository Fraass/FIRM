### Setup file for Fraass & Lowery Planktic Isotope Model ###
# This file needs to be run prior to running the model.
# This file will install the required package (ncdf) to download the Levitus dataset, the raw data
# that the model runs on.
# It also does a little big of transformation of that dataset, simply to be used more simply in the model.
# This _does_ need to have access to the internet to download the files which are 
# on the NOAA FTP site. 
# After this file is run, Planktic Isotope Model (v.#.#).R should run fine.


"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y

#reading in salinity Data

sal.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/salinity/csv/decav/1.00/woa13_decav_s13mn01v2.csv.gz'
tempfile()->tmp
download.file(sal.profile,tmp)
gzfile(tmp,'r+')
sal.profile.win<-read.csv(gzfile(tmp),skip=1,header=T)
names(sal.profile.win)[1:3]<-c('lat','lon','X0')

sal.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/salinity/csv/decav/1.00/woa13_decav_s14mn01v2.csv.gz'
tempfile()->tmp
download.file(sal.profile,tmp)
gzfile(tmp,'r+')
sal.profile.spr<-read.csv(gzfile(tmp),skip=1,header=T)
names(sal.profile.spr)[1:3]<-c('lat','lon','X0')

sal.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/salinity/csv/decav/1.00/woa13_decav_s14mn01v2.csv.gz'
tempfile()->tmp
download.file(sal.profile,tmp)
gzfile(tmp,'r+')
sal.profile.sum<-read.csv(gzfile(tmp),skip=1,header=T)
names(sal.profile.sum)[1:3]<-c('lat','lon','X0')

sal.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/salinity/csv/decav/1.00/woa13_decav_s14mn01v2.csv.gz'
tempfile()->tmp
download.file(sal.profile,tmp)
gzfile(tmp,'r+')
sal.profile.fal<-read.csv(gzfile(tmp),skip=1,header=T)
names(sal.profile.fal)[1:3]<-c('lat','lon','X0')


#reading in temperature data
temp.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t13mn01v2.csv.gz'
tempfile()->tmp
download.file(temp.profile,tmp)
gzfile(tmp,'r+')
temp.profile.win<-read.csv(gzfile(tmp),skip=1,header=T)
names(temp.profile.win)[1:3]<-c('lat','lon','X0')

temp.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t14mn01v2.csv.gz'
tempfile()->tmp
download.file(temp.profile,tmp)
gzfile(tmp,'r+')
temp.profile.spr<-read.csv(gzfile(tmp),skip=1,header=T)
names(temp.profile.spr)[1:3]<-c('lat','lon','X0')

temp.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t15mn01v2.csv.gz'
tempfile()->tmp
download.file(temp.profile,tmp)
gzfile(tmp,'r+')
temp.profile.sum<-read.csv(gzfile(tmp),skip=1,header=T)
names(temp.profile.sum)[1:3]<-c('lat','lon','X0')

temp.profile<-'http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t16mn01v2.csv.gz'
tempfile()->tmp
download.file(temp.profile,tmp)
gzfile(tmp,'r+')
temp.profile.fal<-read.csv(gzfile(tmp),skip=1,header=T)
names(temp.profile.fal)[1:3]<-c('lat','lon','X0')



#Creating 4-D array
  #1st D: longitude
  #2D: latitude
  #3D: depth
  #4D: season (1, NHWinter, 3 NHSummer)

install.packages('abind')
install.packages('reshape2')
library(abind) #to help build the 4-D array
library(reshape2) #needed to quickly generate a grid of temp/salinity data at lat-long coords
At.lon<-unique(temp.profile.win$lon) #all possible latitude coordinates in WOA
At.lat<-unique(temp.profile.win$lat) #all possible longitude coordinates in WOA
seq(min(At.lon),max(At.lon),by=1)->At.lon #fixing as the above misses a few gridcells
seq(-89.5,89.5,by=1)->At.lat #fixing as the above misses a few gridcells

convertWOAdataforFIRM<-function(season.data){
  #output is the 3d Array
  
  output<-dcast(season.data[,1:3],lon~lat,value.var="X0") #matrix of lat/lon for SST
  
  rownames(output)<-output$lon #moving longitude to the 'rowname' attribute
  output[,2:length(output[1,])]->output #removing the first column
  
  
  #coercing the matrix to have the correct number of rows and columns (from jlhoward http://stackoverflow.com/questions/20335637/insert-rows-into-the-middle-of-a-column)
  df       <- data.frame(id=rownames(output),value=output)
  all.rows <- data.frame(id=At.lon,value=0)
  output[which(all.rows$id %in% df$id),]$value <-  df$value #rows corrected
  
  missing<-At.lat %w/o% colnames(output) #finds missing columns
 
  if(length(missing) != 0){
    
    for(o in 1:length(missing)){
      if(length(which(as.numeric(colnames(output)) > missing[o])) == 0){
        output<-cbind(output,
                              matrix(NA,nrow =length(output[,1])))
        names(output)[ncol(output)]<-missing[o]
      }
      if(length(which(as.numeric(colnames(output)) > missing[o])) != 0){
        target<-min(which(as.numeric(colnames(output)) > missing[o]))
        if(target == 1){
          #if target ==1, then 1:(target-1) is nonsense (couple of lines below)
          landing<-matrix(NA,nrow =length(output[,1]))
          output<-cbind(landing,
                                output[,target:length(output[1,]),drop=F])
          names(output)[1]<-missing[o]
        }
        if(target != 1){
          landing<-cbind(output[,1:(target-1),drop=F],
                         matrix(NA,nrow =length(output[,1])))
          
          names(landing)[ncol(landing)]<-missing[o]
          
          output<-cbind(landing,
                                output[,target:length(output[1,]),drop=F])
        }
      } #inserts missing column into correct position
    }
    
  }
  
  #Now do for each depth, adding them as the array's third dimension
  
  for(i in 4:{length(colnames(season.data))}){ #first two columns in season.data
    #are lat,lon, and we've already done the 0 mbsl
    holding.matrix<-dcast(season.data[,c(1:2,i)],
                          lon~lat,
                          value.var=colnames(season.data)[i]
    )
    
    
    
    rownames(holding.matrix)<-holding.matrix$lon #moving longitude to the 'rowname' attribute
    holding.matrix[,2:length(holding.matrix[1,])]->holding.matrix #removing the first column
    
    
    #coercing the matrix to have the correct number of rows and columns (from jlhoward http://stackoverflow.com/questions/20335637/insert-rows-into-the-middle-of-a-column)
    df       <- data.frame(id=rownames(holding.matrix),value=holding.matrix)
    all.rows <- data.frame(id=At.lon,value=0)
    holding.matrix[which(all.rows$id %in% df$id),]$value <-  df$value #rows corrected
    
    missing<-At.lat %w/o% colnames(holding.matrix) #finds missing columns
    
    
    
    if(length(missing) != 0){
      
      for(o in 1:length(missing)){
        if(length(which(as.numeric(colnames(holding.matrix)) > missing[o])) == 0){
          holding.matrix<-cbind(holding.matrix,
                        matrix(NA,nrow =length(holding.matrix[,1])))
          names(holding.matrix)[ncol(holding.matrix)]<-missing[o]
        }
        if(length(which(as.numeric(colnames(holding.matrix)) > missing[o])) != 0){
          target<-min(which(as.numeric(colnames(holding.matrix)) > missing[o]))
          if(target == 1){
            #if target ==1, then 1:(target-1) is nonsense (couple of lines below)
            landing<-matrix(NA,nrow =length(holding.matrix[,1]))
            holding.matrix<-cbind(landing,
                                  holding.matrix[,target:length(holding.matrix[1,]),drop=F])
            names(holding.matrix)[1]<-missing[o]
          }
          if(target != 1){
            landing<-cbind(holding.matrix[,1:(target-1),drop=F],
                           matrix(NA,nrow =length(holding.matrix[,1])))
            
            names(landing)[ncol(landing)]<-missing[o]
            
            holding.matrix<-cbind(landing,
                                  holding.matrix[,target:length(holding.matrix[1,]),drop=F])
          }
        } #inserts missing column into correct position
      }
     } 
    
   
    
    
    
    
    abind(output,
          holding.matrix,
          along=3)->output
  }
    
  
  return(output)
}


#pulling together seasonal data
temp.temp.win<-convertWOAdataforFIRM(temp.profile.win)
rm(temp.profile.win)
temp.temp.spr<-convertWOAdataforFIRM(temp.profile.spr)
rm(temp.profile.spr)
temp.temp.sum<-convertWOAdataforFIRM(temp.profile.sum)
rm(temp.profile.sum)
temp.temp.fal<-convertWOAdataforFIRM(temp.profile.fal)
rm(temp.profile.fal)

#generating temperature 4D array
abind(temp.temp.win,
      temp.temp.spr,
      temp.temp.sum,
      temp.temp.fal,
      along=4)->temp.temp
rm(temp.temp.win)
rm(temp.temp.spr)
rm(temp.temp.sum)
rm(temp.temp.fal)


#pulling together seasonal salinity data
sal.temp.win<-convertWOAdataforFIRM(sal.profile.win)
rm(sal.profile.win)
sal.temp.spr<-convertWOAdataforFIRM(sal.profile.spr)
rm(sal.profile.spr)
sal.temp.sum<-convertWOAdataforFIRM(sal.profile.sum)
rm(sal.profile.sum)
sal.temp.fal<-convertWOAdataforFIRM(sal.profile.fal)
rm(sal.profile.fal)

#pulling together salinity 4D array
abind(sal.temp.win,
      sal.temp.spr,
      sal.temp.sum,
      sal.temp.fal,
      along=4)->salt.salt
rm(sal.temp.win)
rm(sal.temp.spr)
rm(sal.temp.sum)
rm(sal.temp.fal)

#Indexing objects
  #temperature
  t.lon<-At.lon
  t.lat<-At.lat
  #getting depths
    a<-colnames(temp.profile.win)[3:length(colnames(temp.profile.win))]
    for(i in 1:length(a)){
      substring(a[i],2)->a[i]
    }
    as.numeric(a)->a
    t.depth<-a
    rm(a)

  #Salinity (follows temperature for the most part.) 
    #IF you are rewriting this to use a different underlying dataset with
    #unique resolutions in temperature and salinity, this must be changed!
  s.lon<-At.lon
  s.lat<-At.lat
  #getting depths
    a<-colnames(sal.profile.win)[3:length(colnames(sal.profile.win))]
    for(i in 1:length(a)){
      substring(a[i],2)->a[i]
    }
    as.numeric(a)->a
    s.depth<-a
    rm(a)


## Diagnostic Plot [ If a plot doesn't appear, there is an issue above ]
  #extract temperature/salinity profile (Again, diagnostic plot)
  lon<-155 #longitude (not longitude currently, but cell in WOA file)
  lat<-50 #latitude (not latitude currently, but cell in WOA file)
  time<-3 
  temp1<-temp.temp[lon,lat,,time]
  salt1<-salt.salt[lon,lat,,time]
  par(mfcol=c(1,2))
  plot(temp1,
       t.depth,
       type='b',
       ylim=c(5500,0),
       pch=16,
       lwd=2,
       ylab="Depth (m)",
       xlab='Temperature (degC)'
  )
  plot(salt1, #scale for salinity is wrong. Just checking if the data is there.
        s.depth,
        type='b',
        col='red',
       ylab="Salinity",
        ylim=c(5500,0),
        pch=16,
        lwd=2)

