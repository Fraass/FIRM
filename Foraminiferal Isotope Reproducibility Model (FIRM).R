##Planktic Foram Isotope Reproducability##
#A. J. Fraass(1,2) & C. M. Lowery(1,3),  In Revision
#1-University of Massachussetts - Amherst, Department of Geosciences
#2-Present Affiliation: Smithsonian Institution, National Museum of Natural History, Department of Paleobiology
#2-Present Affiliation: University of Texas, Institute for Geophysics


########Number of Repititions######
10000->reps



### Variables ###
## Range of possible depths ##
  depthrange<-c(0,250) #Upper and lower bound for depths in meters

##Location Choice##
  lon<-161 #Cell value E/W within WOA13 dataset)
  lat<-88 #Cell value N/S within WOA13 dataset)


#### NOTE: Find nearest lat.long superceeds the above command ###

###Find the nearest lat./long. coordinates
  find.loc<-'yes' #'yes' or 'no'
  find.lat<-1.4
  find.lon<-157.3
    #160.5410000,2.4330,      Site 803
    #80.5910000,-61.5776667,  Site 744
    #-89.683,-1.216           Koutavas et al. 2006

##  Distribution for depths ##
  dist<-"Even" 
    #possible: 
    #Normal (this uses mean=mean(depthrange), sd={mean(depthrange)-depthrange[1]}/3)
    #Even (even chance for all depths in range), 

##  Number of Individuals ##
  ind<-33

##  Choose d18O equation ##
  equation<-'spero'
    #possible: 
    ##'erez&luz' Cultured Globigerinoides sacculifer (Mixed Layer) 1983
    ##'mielke' Cultured Globorotalia menardii (Thermocline) 2001 [in Spero et al2003]
    ##'spero' Cultured Globigerinoides sacculifer (mixed layer-High light) [in Spero et al2003]
    ##'kim&oneil' Synthetic Calcite 1997 [reformulated by Bemis et al 1998]
    ##'lynch-stieglitz' Insitu Cibicidoides & Planulina [arraged by Cramer et al 2011]
    ##'bemis' Cultured Orbulina universa (high light) [Bemis et al. 1998]

## Mass varies? ##
  mvar<-33
    #Allow mass of individuals to vary. 
    # 0 means all bugs contribute exactly the same amount of gas
    # 10 means bugs could vary up to 10%

##  Diagenesis?  ##
  diag.ind<-0 #Number of individuals diagenetically altered 
  diag<-50 #0-100% of original signal destroyed
    #0% - pure, undiagenetically altered
    #100% - pure d18O bottom water
    #round(1/10*ind) #~10% of individuals diagentically altered

## 'Vital Effect' value ##
  err<-0.146
    #Here using the standard error from the Bemis et al 1998 [low light]
    #That error is ~0.7degC ~= .146permil 

##  Misidentification  ##
  misid<-0 #0-100% chance of including a misidentified bug
  maxmisd<-500 #max depth possible for misid. bug (must be within site lower depth)
    #Suggested Values:
    #undergrad - 10
    #masters - 1
    #PhD - .1
    #faculty - .01
    #geochemist - 100 

#Seasonality
  time<-c(1,2,3,4) #Seasons to choose from (1,2,3,4)
    #Example syntax: time<-c(2,3,4) OR time<-c(1)
  weights<-c(10,20,50,20) #chance that an individual grew in each season (integer, 1 to 100 , should sum to 100)
    #Example syntax: weights<-c(25,50,25)
    #parameter weights is not used a single season model run

#Machine Error
  macer<-0.07 
    #Suggested value: 1sigma error = 0.07 
  

# Tweakable Variables for calculation speed ###
  int<-.1 
    #possible depth point every 10 cm = 0.1 m
    #this has not been tested extensively. 
    #10 cm runs fast enough for current model configuration, and is only here for future expansion.
    #Tweak at own risk.






#### Model calculations are below. ###
  # This model has been explicitly designed to run without any function() calls, and is essentially
  # exposed to allow the user to see, with hopefully minimal R knowledge, to understand the calculations.
  # Those wishing to modify this are more than welcome to, but please at least cite the original publication.
  # Better yet, contact Chris and I, we plan to continue augmenting this for various purposes, and welcome collaboration.

if(find.loc == 'yes'){ #This finds a specific latitude and longitude if using the "location choice" parameter
  #dist matrix 
  a<-dist(c(find.lon,s.lon))
  b<-which(a == min(dist(c(find.lon,s.lon))))
  lon<-s.lon[b[1]]
  lon<-which(s.lon == lon)
  #dist matrix
  a<-dist(c(find.lat,s.lat))
  b<-which(a == min(dist(c(find.lat,s.lat))))
  lat<-s.lat[b[1]]
  lat<-which(s.lat == lat)
  rm(a);rm(b)}

#Plot for Temperature and Salinity ####
  rm(temp1) #clearing a previous run
  rm(salt1)
  temp1<-temp.temp[lon,lat,,time] #temp.temp is the levitus dataset
  if(length(time) > 1){ #this shows up frequently. time is the seasonality parameter, and has to be treated differently if there is a single season or multiple seasons.
   temp1<-temp.temp[lon,lat,,1:4] #temp.temp is the levitus dataset
   min(which(is.na(temp1[,1] == F)))->a
   a-1->a
   temp1<-temp1[1:a,1:4] #in case there isn't data to the base of the Levitus profile
   salt1<-salt.salt[lon,lat,,1:4] #salt.salt is from the levitus dataset
   salt1<-salt1[1:a,1:4]}
  if(length(time) == 1){
    temp1<-temp.temp[lon,lat,,time]
    min(which(is.na(temp1 == F)))->a
    a-1->a
    temp1<-temp1[1:a]
    salt1<-salt.salt[lon,lat,,time]
    salt1<-salt1[1:a]
  }

  par(mfcol=c(1,3))
  ### 
  plot(0,
     0,
     type='n',
     ylim=c(5500,0),
     xlim=c(min(temp1,na.rm=T),max(temp1,na.rm=T)),
     pch=16,
     lwd=2,
     ylab="Depth (m)",
     xlab='Temperature (degC)'
  )
  for(i in time){
    lines(c(temp1[,i],rep(NA,times= length(t.depth) - length(temp1[,i]))),
          t.depth,
          lwd=2)
    points(c(temp1[,i],rep(NA,times= length(t.depth) - length(temp1[,i]))),
           t.depth)
  }
  
  plot(0, 
     0,
     type='p',
     col='red',
     ylim=c(5500,0),
     xlim=c(min(salt1,na.rm=T),max(salt1,na.rm=T)),
     ylab="Depth (m)",
     xlab='Salinity PSU')
  for(i in time){
    lines(c(salt1[,i],rep(NA,times= length(s.depth) - length(salt1[,i]))),
          s.depth,
          lwd=2,
          col='red')
    points(c(salt1[,i],rep(NA,times= length(s.depth) - length(salt1[,i]))),
           s.depth,
           col='red')
  }
  #Interpolation for Temperature & Salinity from Levitus '98 Profiles
  temp<-rep(NA,times={max(t.depth)-min(t.depth)}/int+1)  #because 0 must be a point
  salt<-temp
  temp.depth<-seq(min(t.depth),max(t.depth),by=int) #adding a depth scale
  salt.depth<-seq(min(s.depth),max(s.depth),by=int)



  ##### Calculation of d18O values ###
  if(length(time) == 1){   #d18O profile for a single season
      temp1[1]->temp[1]
      salt1[1]->salt[1]
    #linear interpolation
    for(i in 2:length(t.depth)){
      {t.depth[i]-t.depth[i-1]}/int->n.p #number of points 
      {temp1[i]-temp1[i-1]}/n.p->dif 
      for(o in 1:n.p){
       temp[o+t.depth[i-1]/.1]<-temp1[i-1]+dif*{o}
      }
    } 
    temp1[length(temp1)]->temp[length(temp)]
    lines(temp,temp.depth)
    #interpolating for salinity
    for(i in 2:length(s.depth)){
      {s.depth[i]-s.depth[i-1]}/int->n.p #number of points 
      {salt1[i]-salt1[i-1]}/n.p->dif
      for(o in 1:n.p){
       salt[o+s.depth[i-1]/.1]<-salt1[i-1]+dif*{o}
      }
    }
  salt1[length(salt1)]->salt[length(salt)]
  lines(salt/max(salt)*400-385, 
        salt.depth,col="red")
  #Calculating d18Osw from salinity profile [Equation from Conroy et al. 2014]####
  d18O.sw<-NA
  d18O.sw[1:which(salt.depth == 77.5)]<- -7.82 + 0.23*salt[1:which(salt.depth == 77.5)] #0-77.5m
  d18O.sw[
    {which(salt.depth == 77.5)+1}:length(salt)
    ]<- -11.83 + 0.35*salt[
      {which(salt.depth == 77.5)+1}:length(salt)] #77.4m-Seafloor
  }

  if(length(time) > 1){ #calculation of d18O for multiple seasons
    temp<-matrix(ncol=4,nrow={max(t.depth)-min(t.depth)}/int+1)  #because 0 must be a point
    salt<-temp
    temp.depth<-seq(min(t.depth),max(t.depth),by=int) #adding a depth scale
    salt.depth<-seq(min(s.depth),max(s.depth),by=int)
    temp1[1,]->temp[1,]
    salt1[1,]->salt[1,]
    # interpolating for temperature
    for(p in 1:4){
      for(i in 2:length(temp1[,1])){
        {t.depth[i]-t.depth[i-1]}/int->n.p #number of points 
        {temp1[i,p]-temp1[i-1,p]}/n.p->dif 
        for(o in 1:n.p){
          temp[o+t.depth[i-1]/int,p]<-temp1[i-1,p]+dif*{o}
          }
        } 
      }
    temp1[length(temp1[,1]),]->temp[length(temp[,1]),]
    plot(0,0,
       ylim=c(max(salt.depth),0),
       xlim=c(max(temp,na.rm=T),min(temp,na.rm=T)),
      ylab='temperature',
       type='n'
    )
    for(o in 1:4){
      lines(temp[,o],temp.depth[1:length(temp[,o])])
    }
    #interpolating for salinity
    for(p in 1:4){
      for(i in 2:length(salt1[,1])){
        {s.depth[i]-s.depth[i-1]}/int->n.p #number of points 
        {salt1[i,p]-salt1[i-1,p]}/n.p->dif
        for(o in 1:n.p){
          salt[o+s.depth[i-1]/int,p]<-salt1[i-1,p]+dif*{o}
          }
        }
      }
    salt1[length(salt1[,1]),]->salt[length(salt[,1]),]
    for(o in 1:4){
      lines(salt[,o]/max(salt[,o])*400-385,
          salt.depth[1:length(salt[,o])],col='red')
    }
    plot(0,0,
       ylim=c(max(salt.depth),0),
       xlim=c(max(salt,na.rm=T),min(salt,na.rm=T)),
       type='n',ylab='salinity'
    )
    for(o in 1:4){
      lines(salt[,o],salt.depth[1:length(salt[,o])],col='red')
    }
  
  #Calculating d18Osw from salinity profile [Equation from Conroy et al. 2014]####
  d18O.sw<-matrix(nrow=length(salt[,1]),ncol=4)
  d18O.sw[1:which(salt.depth == 77.5),]<- -7.82 + 0.23*salt[1:which(salt.depth == 77.5),] #0-77.5m
  d18O.sw[
    {which(salt.depth == 77.5)+1}:length(salt[,1]),
    ]<- -7.82 + 0.23*salt[
      {which(salt.depth == 77.5)+1}:length(salt[,1]),] #77.4m-Seafloor
  
}

  if(length(time) == 1){ #diagnostic plots
    plot(d18O.sw,
       salt.depth,
       type='l',
       ylim=c(max(salt.depth),0),
       pch=16,
       lwd=2,
       ylab="Depth (m)",
       xlab='d18Osw'
    )}
  if(length(time) > 1){#diagnostic plots
    plot(d18O.sw[,1],
       salt.depth[1:length(d18O.sw[,1])],
       type='l',
       ylim=c(max(salt.depth),0),
       pch=16,
       lwd=2,
       ylab="Depth (m)",
       xlab='d18Osw'
    )
    for(o in 2:4){#diagnostic plots
      lines(d18O.sw[,o],
          salt.depth[1:length(d18O.sw[,1])])}
  }



# Repetition code starts here ####
    d18Ostore<-matrix(NA,ncol=reps,nrow=ind) #creating object to store entire synthetic dataset, including each individual 'specimen'
    d18Otot.store<-rep(NA,times=reps) #creating object to store synthetic mass spec output
  temp.store<-matrix(NA,ncol=reps,nrow=ind)
  d18Owd.store<-matrix(NA,ncol=reps,nrow=ind)
  for(A in 1:reps){
    #Sampling of Depths ###
      if(dist == "Normal"){
        rnorm(n=ind,
              mean=mean(depthrange), 
              sd={mean(depthrange)-depthrange[1]}/3
        )->depth.samples
        depth.samples[which(depth.samples < 0)]<-0
        #y <- dnorm(temp.depth,
        #           mean=mean(depthrange),
        #           sd={mean(depthrange)-depthrange[1]}/3)
        #lines(max(salt1)-y/{max(salt1)-min(salt1)}*10,temp.depth)
        #rm(y)
      }
      if(dist == "Even"){
        seq(from=temp.depth[which(temp.depth == min(depthrange))],
            to=temp.depth[which(temp.depth == max(depthrange))],
            by=int
        )->d
        sample(d,
               size=ind,
               replace=TRUE
        )->depth.samples
        rect(max(salt1),
             depthrange[1],
             max(salt1)-.05,
             depthrange[2],
             col='black')
      }
    
    #Seasonality sampling ###
      if(length(time) == 1){season<-"NO"}
      if(length(time) == 2){
        if(ind <= 100){
          season<-c(rep(time[1],times=weights[1]),
                  rep(time[2],times=weights[2]))}
        if(ind > 100){
          season<-c(rep(time[1],times=weights[1]*ind/100),
                  rep(time[2],times=weights[2]*ind/100))
        }
        
      }
      if(length(time) == 3){
        if(ind <= 100){
        season<-c(rep(time[1],times=weights[1]),
                  rep(time[2],times=weights[2]),
                  rep(time[3],times=weights[3]))
        }
        if(ind > 100){
          season<-c(rep(time[1],times=weights[1]*ind/100),
                    rep(time[2],times=weights[2]*ind/100),
                    rep(time[3],times=weights[3]*ind/100))
        }
          
        } #My apologies to anybody who got this far in the code. I got really lazy here.
      if(length(time) == 4){
        if(ind <= 100){
        season<-c(rep(time[1],times=weights[1]),
                  rep(time[2],times=weights[2]),
                  rep(time[3],times=weights[3]),
                  rep(time[4],times=weights[4]))
        }
        if(ind > 100){
          season<-c(rep(time[1],times=weights[1]*ind/100),
                    rep(time[2],times=weights[2]*ind/100),
                    rep(time[3],times=weights[3]*ind/100),
                    rep(time[4],times=weights[4]*ind/100))
        }
      }
      if(length(time) > 1){sample(season,size=ind)->season.sample}
    
    #Missidentification (The Undergrad Lab assistant Correction) [  ] ####
      if(misid != 0){
        if(misid >=1){misseq<-seq(0,100,by=1)}
        if(misid >=.1 & misid <1){misseq<-seq(0,100,by=.1)}
        if(misid <.1){misseq<-seq(0,100,by=.01)}
        sample(misseq,size=ind)->mis
        for(i in 1:ind){
          if(mis[i] <= misid){
            sample(seq(0,maxmisd,by=int),size=1)->depth.samples[i]} #finding individuals to be replaced
          if(season != "NO"){
            sample(1:4,size=1)->season.sample[i] #if there is no seasonality included, here we grab a random bug from a random depth during a random season
            }
          }
      }
      temp.s<-NA
      salt.s<-NA
      d18O.wd<-NA
      if(season[1] != "NO"){
        for(w in 1:ind){
          temp.s[w]<-temp[depth.samples[w]*10+1,season.sample[w]]
          #salt.s[w]<-salt[depth.samples[w]*10+1,season.sample[w]]
          d18O.wd[w]<-d18O.sw[depth.samples[w]*10+1,season.sample[w]]
        }
      }
      if(season == "NO"){
        temp.s<-temp[depth.samples*10+1] #THIS NEEDS TO BE FIXED FOR CHANGING THE int PARAMETER
        #salt.s<-salt[depth.samples*10+1]
        d18O.wd<-d18O.sw[depth.samples*10+1]
      }
    
    #Calculating pure oxygen isotope signal for each individual ####
      if(equation == "erez&luz"){
        d18O<-{
          #{d18O.wd+-0.22}-{100/449}*{temp.s-17} #-0.22 is to convert VSMOW to VPDB
          #d18O.wd+{{89939-5000*temp.s}/22450}
          d18O.wd-2/3*sqrt(75*temp.s+11494)+{11267/150}
        }  
      }
      if(equation == "kim&oneil"){
        d18O<-{
          #d18O.wd-{20/91}*temp.s+ {4951/1300} #-0.27 is to convert VSMOW to VPDB
          {1/900}*{900*d18O.wd-100*sqrt(2)*sqrt(450*temp.s+19667)+22957}
        }  
      }
      if(equation == "mielke"){
        d18O<-{
          #{d18O.wd+-0.27}-{5/257}*{10*temp.s- 149} #-0.27 is to convert VSMOW to VPDB
          d18O.wd-{100*temp.s/513}+{135149/51300}
        }  
      }
      if(equation == "spero"){
        d18O<-{
          #{d18O.wd+-0.27}-100*{temp.s-12}/557 #-0.27 is to convert VSMOW to VPDB
          d18O.wd-{100*temp.s/557}+{104961/55700}
        }  
      }
      if(equation == "lynch-stieglitz"){
        d18O<-{
          #{d18O.wd+-0.27}+5/238*{161-10*temp.s} #-0.27 is to convert VSMOW to VPDB
          d18O.wd-{20*temp.s/87}+{9157/2900}
        }  
      }
      if(equation == "bemis"){
        d18O<-{
          #{d18O.wd+-0.27}+5/238*{161-10*temp.s} #-0.27 is to convert VSMOW to VPDB
          d18O.wd-{5*temp.s/24}+{3401/1200}
        }  
      }      
    
    #Modifying pure oxygen isotope signal with 'vital effect' error value [ See notes in in parameter section ] ####
      if(err != 0){
        for(i in 1:length(d18O)){
          #sample(seq(from=-err, to=err, by=0.001),size=1)->q
          rnorm(n=1,mean=0,sd=err)->q #modifying due to machine error
          d18O[i]<-d18O[i]+q
          rm(q)
        }}
    
    #Carbonate Ion / pH [ Not used for Planktics ] ####

    #Individuals of differing masses ####
      contrib<-rep(NA,times=ind)
      seq(100,100-mvar,by=-.01)->t # getting the variations in mass
      for(i in 1:ind){
        contrib[i]<-sample(t,size=1) # randomizing the mass variations
      }
      rm(t)
      contrib<-contrib/sum(contrib)*100 # creating an object with the mass contributions to final synthetic mass spec value
        
    #Diagenesis [ Only uses Kim&oneil equation ] ####
      temp1[length(temp1)]->tb
      d18O.wd[length(d18O.wd)]->tsw
      d18Ob<-{tsw+-0.27}-{20/91}*tb+ {46/13} #finding bottom water calcite
      #-0.27 is to convert VSMOW to VPDP
      rm(tb);rm(tsw)
      d18O-d18Ob->difO #difference between sampled values and the bottom inorganic calcite
      -difO*{diag/100}+d18O->d18Od
    

    #Calculating single value for isotopes ####
    #Diagenetic alterations
      d18O.a<-d18O
      if(diag.ind > 0){
        sample(1:ind,size=diag.ind)->t
        d18O.a[t]<-d18Od[t]
        rm(t)
      }
      if(mvar == 0){
        mean(d18O.a)->t.d18O # no mass variability
      }
      if(mvar != 0){
        sum(d18O.a*contrib)/100->t.d18O # modifying the contributions from each individual for the final synthetic mass value
      }
  
    #storing Data ####
      d18Ostore[,A]<-d18O.a
      d18Otot.store[A]<-t.d18O
      temp.store[,A]<-temp.s
      d18Owd.store[,A]<-d18O.wd
      #print(A) # Check progress in the model
}
    #Machine Error
      if(macer != 0){
        rnorm(n=reps,mean=0,sd=macer)->a #modifying due to machine error
        a+d18Otot.store->d18Otot.store
      }

#Summary Plots#####
  par(mfcol=c(4,1))
  hist(temp.store
       ,col='red'
       ,border=F
       ,xlab='Temperature'
       ,main="Histogram of Individual Values"
       ,xlim=c(min(temp.store),max(temp.store))) #for no diag
  #,xlim=c(d18Ob,max(d18O))) #for lots of diag
  abline(v=mean(temp.store),lwd=3)
  hist(d18Owd.store
       ,col='purple'
       ,border=F
       ,xlab='d18Osw'
       ,main="Histogram of Individual Values"
       ,xlim=c(min(d18Owd.store),max(d18Owd.store))) #for no diag
  #,xlim=c(d18Ob,max(d18O))) #for lots of diag
  abline(v=mean(d18Owd.store),lwd=3)
  hist(d18Ostore
       ,col='forestgreen'
       ,border=F
       ,xlab='d18O'
       ,main="Histogram of Individual Values"
       ,xlim=c(min(d18Ostore),max(d18Ostore))) #for no diag
  #,xlim=c(d18Ob,max(d18O))) #for lots of diag
  abline(v=mean(d18Ostore),lwd=3)
  hist(d18Otot.store
       ,col='blue'
       ,xlab='d18O'
       ,border=F
       ,main=paste(paste('Histogram of', ind,sep=' '),'Individuals',sep=' ')
       ,xlim=c(min(d18Ostore),max(d18Ostore))) #for no diag
  #,xlim=c(d18Ob,max(d18O)) #for lots of diag   
  sort(d18Otot.store)[c(.025*reps,.975*reps)]#95%C.I.
  abline(v=sort(d18Otot.store)[c(.025*reps,.975*reps)],
         col='red',
         
         lwd=2)
  abline(v=mean(d18Otot.store),lwd=3)
  
  
  #Calculating various results
    {
      print(paste('High Bound (95%) = ',
                  round(
                    sort(
                      d18Otot.store
                      )[.025*reps]
                    ,digits=2)
                  )
            )
      print(paste('Low Bound (95%)  = ',
                  round(
                    sort(
                      d18Otot.store
                      )[.975*reps]
                    ,digits=2
                    )
                  )
            )
      print(paste('Mean Value       = ',
                  round(
                    mean(
                    d18Otot.store
                             )
                        ,digits=2
                        )
                  )
            )
      print(paste('+/- 95%          = ',
                  round(
                    abs(
                      {sort(
                        d18Otot.store
                        )[c(.025*reps)]-sort(
                          d18Otot.store
                          )[c(.975*reps)]}/2)
                    ,digits=2
                    )
                  )
            )
      print(paste('SD               = ',
                  sd(d18Otot.store)
      )
      )
    }
