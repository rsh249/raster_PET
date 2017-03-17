##############
###R script for calculating Thornthwaite-Mather Potential Evapotranspiration###
##This was written by:
#####Robert S. Harbert#########
#####rharbert@amnh.org#########
##Direct all questions here####
##T-M PET is a quick and dirty estimate$
## of potential evapotranspiration based on
## Thornthwaite (1948). Relying only on
## mean monthly temperatures and daylength
## this can be calculated quickly and easily
## from Worldclim (www.worldclim.org -- Hijmans etal. 2005)
## raster files -- tmean -- and a daylength formula.
## While not a substitute for more sophisticated
## methods of PET estimation, Thornthwaite's method
## is reasonable at coarse spatial scales and
## long periods of time (i.e., monthly or annual) 
## [see: Lu et al. 2005; doi:10.1111/j.1752-1688.2005.tb03759.x]

##############


###Software notes###
# Note that this script tries to keep all raster files in RAM
# As such, it is recommended to run this script on a machine with at least 16GB of RAM.
#
####


library(raster)
rasterOptions(maxmemory=200000000);
#ext = extent(c(-150, -40, 10, 60)) ## Set to study area if desired
#get mean monthly temperatures and precipitation grids. 
tmean <- get_worldclim(varset='tmean', res=10) ## Update spatial resolution if needed (options 30, 2.5, 10 corresponding to 30 arcseconds, 2.5 arcminutes, and 10 arcminutes respectively)
tmean = stack(tmean[[c(1,5:12)]], tmean[[2:4]])
#tmean = brick(tmean)
if(max(maxValue(na.omit(tmean)))>100){
  tmean = tmean/10;#necessary if using WorldClim temperature grids 
}
prec <- get_worldclim(varset='prec', res=10);
prec = stack(prec[[c(1,5:12)]], prec[[2:4]])
prec = brick(prec)

#tmean = crop(tmean, ext);
#prec = crop(prec, ext);
#Matrix to raster
 # Calculate daylength (avg) for the month
 library(geosphere)
 days = tmean[[1]];
 for(i in 1:12){
   v = extent(days)@ymax
   grid <- matrix(nrow=nrow(days), ncol=ncol(days))
   for(r in 1:nrow(days)){
     v=v-res(days)[1]
     #print(v)
     J = (i*30)-15; #Julian Day of year
     dl = daylength(v, J)
     #print(dl);
     for(c in 1:ncol(days)){
         grid[r,] = dl;
     }
       #print(dl); print(J); print(v);
   }
   ras <- raster(grid)
   e <- extent(days)
   extent(ras) <- e
   days <- stack(days, ras)
   rm(ras)
   rm(grid)
   print(i)
 }

 days = days[[-1]]
 days = brick(days)
 days = mask(days, tmean[[1]])
 names(days) = c('Jan', 'Feb', 'Mar', "Apr", 'May', 'Jun', 'Jul', "Aug", 'Sep', 'Oct', 'Nov', 'Dec')

 #plot(days)


#writeRaster(days, file='dl')
#dl <- stack('dl.gri')
#dl <- mask(dl, tmean[[1]])

##Thorwaite PET

#Calculate I -- Annual "heat index"
I_hold = tmean[[1]]> 5000 ; # should be 0 everywhere so does not effect summing at the end

m <- nlayers(tmean)
for (i in 1:m){
  tm =tmean[[i]];
  tm1 = tm>=0; #But only where tmean > 0, elsewhere PET is 0?
  tm = tm*tm1
  I = ((tm)/5)
  I=I**1.514
  I_hold = stack(I_hold, I)
  
}
I_hold = I_hold[[-1]]
I_hold = brick(I_hold);
names(I_hold) = names(days)
Ival = sum(I_hold);

alpha = ((0.000000675)*Ival**3) - ((0.0000771)*Ival**2) + ((0.01792)*Ival) + 0.49239

#N array of number of days each month
N <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

PET <- function(n, daylength, monthd, TM, alpha, I) { #n = iteration (month) # 'a' = daylength grid, 'b' = N -- array of days in month, 'c' tmean grid, 'd' alpha grid
  mdays = (daylength/12) * (monthd/30)
  #plot(mdays)
  P = 16 * mdays *  (((10 * TM)/I)**alpha); #in mm/month
  tm1 = TM>=0;
  P = P*tm1;
  return(P)
}
potET <- I_hold[[1]]
for (i in 1:12) {
  est <- PET(i, daylength = days[[i]], monthd = N[i], TM = tmean[[i]], alpha = alpha, I = Ival)
  potET <- stack(potET, est)
}
potET=potET[[-1]]
potET[is.na(potET)] = 0;
potET = mask(potET, tmean[[1]])
