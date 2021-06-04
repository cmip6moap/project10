###################################################
# Code By Chloe Brimicombe and Claudia Di Napoli  #
# University of Reading for the CMIP6 Hackathon   #
# June 2021                                       #
###################################################


#import statments
require(ncdf4)
require(raster)
require(rasterVis)
require(rworldmap)
require(map)
require(maptools)
require(lattice)
require(latticeExtra)
require(RColorBrewer)
require(animation)
require(ggplot2)

#functions

# removes the years to keep only 2020-2100
#remove -1 for bcc model
limitraster <- function(raster){
  raster <- raster[[61:(nlayers(raster))-1]]
  return(raster)
}

#mean over the year
yearly_mean<- function(month){
  month <- crop(month,c(-180,180,-60,90))
  nmonth<- nlayers(month)
  monthmean <- month[[1:(nmonth/12)]]
  i=1
  y=1
  while(i<= nlayers(month)){
    monthmean[[y]] <- calc(month[[i:(i+11)]], function(x){mean(x, na.rm=T)})
    i<- i+12
    y<- y+1
  }
  return(monthmean)
}
#max over the year
yearly_max<- function(month){
  month <- crop(month,c(-180,180,-60,90)) #crops out Antarctica
  nmonth<- nlayers(month)
  monthmax <- month[[1:(nmonth/12)]]
  i=1
  y=1
  while(i<= nlayers(month)){
    monthmax[[y]] <- calc(month[[i:(i+11)]], function(x){max(x, na.rm=T)})
    i<- i+12
    y<- y+1
  }
  return(monthmax)
}
#decadal mean
decadal_mean<- function(year){
  year <- year[[1:(nlayers(year)-1)]]
  nyear<- nlayers(year)
  yearmean <- year[[1:(nyear/10)]]
  i=1
  y=1
  while(i<= nlayers(year)){
    yearmean[[y]] <- calc(year[[i:(i+9)]], function(x){mean(x)})
    i<- i+10
    y<- y+1
  }
  return(yearmean)
}
#decadal max
decadal_max<- function(year){
  year <- year[[1:(nlayers(year)-1)]]
  nyear<- nlayers(year)
  yearmax <- year[[1:(nyear/10)]]
  i=1
  y=1
  while(i<= nlayers(year)){
    yearmax[[y]] <- calc(year[[i:(i+9)]], function(x){max(x)})
    i<- i+10
    y<- y+1
  }
  return(yearmax)
}


#read in and limit
bcc126 <- limitraster(brick("bcc126.nc"))
bcc245 <- limitraster(brick("bcc245.nc"))
bcc585 <- limitraster(brick("bcc585.nc"))

hadgem126 <- limitraster(brick("hadgem126.nc"))
hadgem245 <- limitraster(brick("hadgem245.nc"))
hadgem585 <- limitraster(brick("hadgem585.nc"))
nlayers(bcc245)
cmcc126 <- limitraster(brick("cmcc126.nc"))
cmcc245 <- limitraster(brick("cmcc245.nc"))
cmcc585 <- limitraster(brick("cmcc585.nc"))

#model monthly mean
p126 <- (bcc126 + hadgem126 + cmcc126) / 3
p245 <- (bcc245 + hadgem245 + cmcc245) / 3
p585 <- (bcc585 + hadgem585 + cmcc585) / 3

#annual mean
ybcc126 <- yearly_mean(bcc126)
ybcc245 <- yearly_mean(bcc245)
ybcc585 <- yearly_mean(bcc585)

yhadgem126 <- yearly_mean(hadgem126)
yhadgem245 <- yearly_mean(hadgem245)
yhadgem585 <- yearly_mean(hadgem585)

ycmcc126 <- yearly_mean(cmcc126)
ycmcc245 <- yearly_mean(cmcc245)
ycmcc585 <- yearly_mean(cmcc585)

#annual max
ybcc126m <- yearly_max(bcc126)
ybcc245m <- yearly_max(bcc245)
ybcc585m <- yearly_max(bcc585)

yhadgem126m <- yearly_max(hadgem126)
yhadgem245m <- yearly_max(hadgem245)
yhadgem585m <- yearly_max(hadgem585)

ycmcc126m <- yearly_max(cmcc126)
ycmcc245m <- yearly_max(cmcc245)
ycmcc585m <- yearly_max(cmcc585)


# it masks the cells under 39 to NAs
p245un <- py245m -273.15
p126un <- py126m -273.15
p585un <- py585m -237.15

p126un[p126un[] < 39 ] = NA 
p245un[p245un[] < 39 ] = NA 
p585un[p585un[] < 39 ] = NA 

#masking to land surface area proportions for RCP4.5 and RCP2.6
world<-getMap('world')
py245unhab <- py245m - 273.15 > 39
py245unhabtotal <- py245m - 273.15 > -100
py245unhabland <- mask(py245unhab,world)
py245unhabtotalmask <- mask(py245unhabtotal,world)
pysum245unhab <- cellStats(py245unhab,sum)
pysum245unhabtotal <- cellStats(py245unhabtotalmask,sum)
prop245 <- (pysum245unhab/pysum245unhabtotal)*100

py126unhab <- py126m - 273.15 > 39
py126unhabtotal <- py126m - 273.15 > -100
py126unhabland <- mask(py126unhab,world)
py126unhabtotalmask <- mask(py126unhabtotal,world)
pysum126unhab <- cellStats(py126unhab,sum)
pysum126unhabtotal <- cellStats(py126unhabtotalmask,sum)
prop126 <- (pysum126unhab/pysum126unhabtotal)*100

dfproj <- data.frame()
#plotting yearly max area above 39C
plot(c(2020:2100),prop245,type="o",pch=16,col="darkorchid4",xlab="Years",ylab="Percentage of Land Surface(%)")
lines(c(2020:2100),prop126,type="o",pch=15,col="cadetblue",xlab="Years",ylab="Percentage of Land Surface(%)",axes=F)
legend(2020, 5, legend=c("RCP 4.5", "RCP 2.6"),
       col=c("darkorchid4", "cadetblue"), 
       pch=c(16,15),lty=1:1, cex=0.8,box.lty = 0)
plot(p126un[[1]])

cp126ym <- cellStats(p126un,mean)
cp245ym <-cellStats(p245un,mean)
cp585ym <- cellStats(p585un,mean)

#all of the models plot together as above
years<- rep.int(2020:2100,1)
length(years)
modelsnamesy <- c(rep('mean126',81),rep('mean245',81),rep('mean585',81))

dfyearunmean <- data.frame(years=c(years,years,years),
                           models=c(cp126ym,cp245ym,cp585ym),
                           modelname= modelsnamesy)
plot(cp585ym)
ggplot(data=dfyearunmean,aes(x=years,y=models,col=modelname))+geom_line()

#model annual mean

py126m <- (ybcc126m + yhadgem126m + ycmcc126m) / 3
py245m <- (ybcc245m + yhadgem245m + ycmcc245m) / 3
py585m <- (ybcc585m + yhadgem585m + ycmcc585m) / 3

#


#decadal mean
dbcc126 <- decadal_mean(bcc126)
dbcc245 <- decadal_mean(bcc245)
dbcc585 <- decadal_mean(bcc585)

dhadgem126 <- decadal_mean(hadgem126)
dhadgem245 <- decadal_mean(hadgem245)
dhadgem585 <- decadal_mean(hadgem585)

dcmcc126 <- decadal_mean(cmcc126)
dcmcc245 <- decadal_mean(cmcc245)
dcmcc585 <- decadal_mean(cmcc585)

pd126 <- decadal_mean(py126)
pd245 <- decadal_mean(py245)
pd585 <- decadal_mean(py585)

pd126m <- decadal_max(py126)
pd245m <- decadal_max(py245)
pd585m <- decadal_max(py585)

dbcc126m <- decadal_max(ybcc126m)
dbcc245m <- decadal_max(ybcc245m)
dbcc585m <- decadal_max(ybcc585m)

dhadgem126m <- decadal_max(yhadgem126m)
dhadgem245m <- decadal_max(yhadgem245m)
dhadgem585m <- decadal_max(yhadgem585m)

dcmcc126m <- decadal_max(ycmcc126m)
dcmcc245m <- decadal_max(ycmcc245m)
dcmcc585m <- decadal_max(ycmcc585m)

dbcc126mm <- decadal_mean(ybcc126m)
dbcc245mm <- decadal_mean(ybcc245m)
dbcc585mm <- decadal_mean(ybcc585m)

dhadgem126mm <- decadal_mean(yhadgem126m)
dhadgem245mm <- decadal_mean(yhadgem245m)
dhadgem585mm <- decadal_mean(yhadgem585m)

dcmcc126mm <- decadal_mean(ycmcc126m)
dcmcc245mm <- decadal_mean(ycmcc245m)
dcmcc585mm <- decadal_mean(ycmcc585m)

#model average 
dmyp126 <- (ybcc126m + yhadgem126m + ycmcc126m) / 3
dmyp245 <- (ybcc245m + yhadgem245m + ycmcc245m) / 3
dmyp585 <- (ybcc585m + yhadgem585m + ycmcc585m) / 3

#model average of the decadal mean of the yearly max
dmmp126 <- (dcmcc126mm + dhadgem126mm + dbcc126mm) / 3
dmmp245 <- (dcmcc245mm + dhadgem245mm + dbcc245mm) / 3
dmmp585 <- (dcmcc585mm + dhadgem585mm + dbcc585mm) / 3

plot(pd126-273.15)

# above uninhabitable
py126un <- py126-273.15 > 38.8
py245un <- py245-273.15 > 38.8
py585un <- py585-273.15 > 38.8

#cellstats max
cbcc126m <- cellStats(bcc126,max)
chadgem126m <- cellStats(hadgem126,max)
ccmcc126 <- cellStats(cmcc126,max)
cp126m <- cellStats(p126,max)

cbcc245m <- cellStats(bcc245,max)
chadgem245m <- cellStats(hadgem245,max)
ccmcc245 <- cellStats(cmcc245,max)
cp245m <- cellStats(p245,max)

cbcc585m <- cellStats(bcc585,max)
chadgem585m <- cellStats(hadgem585,max)
ccmcc585 <- cellStats(cmcc585,max)
cp585m <- cellStats(p585,max)




#cellstats mean

cbcc126me <- cellStats(bcc126,mean)
chadgem126me <- cellStats(hadgem126,mean)
ccmcc126me <- cellStats(cmcc126,mean)
cp126me <- cellStats(p126,mean)

cbcc245me <- cellStats(bcc245,mean)
chadgem245me <- cellStats(hadgem245,mean)
ccmcc245me <- cellStats(cmcc245,mean)
cp245me <- cellStats(p245,mean)

cbcc585me <- cellStats(bcc585,mean)
chadgem585me <- cellStats(hadgem585,mean)
ccmcc585me <- cellStats(cmcc585,mean)
cp585me <- cellStats(p585,mean)


#plotting
cutpts<-c(-100,-40,-27,-13,0,9,26,32,38,46,80)
extreme_cs <-rgb(8/255,48/255,107/255,1)
verystrong_cs<- rgb(8/255,81/255,156/255,1)
strong_cs<-rgb(33/255,113/255,181/255,1)
moderate_cs <-rgb(66/255,146/255,198/255,1)
slight_cs <-rgb(158/255,202/255,225/255,1)
nothermalstress <-rgb(217/255,240/255,163/255,1)
my_palette <-c(extreme_cs,extreme_cs,verystrong_cs,strong_cs,moderate_cs,slight_cs,nothermalstress,'darkorange','orangered','red3','red4','darkorchid4')
col2hex(my_palette)
my_palette <- rev(my_palette)
#grid::grid.raster(my_palette, int=F)
#myColorkey <- list(at=cutpts, labels=list(at=cutpts,
                                         # labels=cutpts))
color2hex(my_palette)

plot(dmmp245)
names(dmmp245)<- c(".2020s",".2030s",".2040s",".2050s",".2060s",
                   ".2070s",".2080s",".2090s")
names(dmmp126)<- c(".2020s",".2030s",".2040s",".2050s",".2060s",
                   ".2070s",".2080s",".2090s")
names(dmmp585)<- c(".2020s",".2030s",".2040s",".2050s",".2060s",
                   ".2070s",".2080s",".2090s")

d1 <- dmmp245[[1]]
plt1 <- rasterVis::levelplot(dmmp585-273.15, 
                             main = "RCP 8.5 Decadal Mean of Max",
                             layout=c(2,4),
                             xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                             col.regions=my_palette,at=cutpts) #cutcol)
plt1 + latticeExtra::layer(sp.lines(countries, col="black", lwd=1))


#plotting line graph

modelnamem <- c(rep("mean",length(cp126m)),rep("cbcc",length(cp126m)),
               rep("cmcc",length(cp126m)),rep("chadgem",length(cp126m)))

modelname <- c(rep("max",length(cp126m)),rep("cbcc",length(cp126m)),
                rep("cmcc",length(cp126m)),rep("chadgem",length(cp126m)))


index<- rep.int(1:12,length(cp126m)/12)
count <- rep.int(1,length(cp126m)/12)
years<- rep.int(2020:2100,12)
yearso <- sort(years)
plot(cp126m-273.15,col="black",size="2",type="o")

df126 <- data.frame(years = c(yearso,yearso,yearso,yearso),
                    index= c(index,index,index,index),
                    models =c(cp126m,cbcc126m,ccmcc126,chadgem126m)-273.15,
                    modelname = modelname,
                    count=count)
df126mm <- data.frame(years = c(yearso),
                    index= c(index),
                    models =c(cp126m)-273.15,
                    modelname = rep("max",length(cp126m)),
                    count=count)

df245 <- data.frame(years = c(yearso,yearso,yearso,yearso),
                    index= c(index,index,index,index),
                    models =c(cp245m,cbcc245m,ccmcc245,chadgem245m)-273.15,
                    modelname = modelname)

df585 <- data.frame(years = c(yearso,yearso,yearso,yearso),
                    index= c(index,index,index,index),
                    models =c(cp585m,cbcc585m,ccmcc585,chadgem585m)-273.15,
                    modelname = modelname)

df126m <- data.frame(years = c(yearso,yearso,yearso,yearso),
                    index= c(index,index,index,index),
                    models =c(cp126me,cbcc126me,ccmcc126me,chadgem126me)-273.15,
                    modelname = modelnamem)

df245m <- data.frame(years = c(yearso,yearso,yearso,yearso),
                    index= c(index,index,index,index),
                    models =c(cp245me,cbcc245me,ccmcc245me,chadgem245me)-273.15,
                    modelname = modelnamem)

df585m <- data.frame(years = c(yearso,yearso,yearso,yearso),
                    index= c(index,index,index,index),
                    models =c(cp585me,cbcc585me,ccmcc585me,chadgem585me)-273.15,
                    modelname = modelnamem)


ggplot(data=df126)+geom_point(aes(x=years,y=models,col=modelname),size=2)+
  theme_bw()+geom_hline(yintercept=39,size=1)+xlab("Years")+
  ylab("UTCI (°C)")+theme(text = element_text(size = 20))+
  facet_wrap(~index)

ggplot(data=df245)+geom_point(aes(x=years,y=models,col=modelname),size=2)+
  theme_bw()+geom_hline(yintercept=38.8,size=1)+xlab("Years")+
  ylab("UTCI (°C)")+theme(text = element_text(size = 20))+
  facet_wrap(~index)

ggplot(data=df585)+geom_point(aes(x=years,y=models,col=modelname),size=2)+
  theme_bw()+geom_hline(yintercept=38.8,size=1)+xlab("Years")+
  ylab("UTCI (°C)")+theme(text = element_text(size = 20))+
  facet_wrap(~index)

ggplot(data=df126m)+geom_point(aes(x=years,y=models,col=modelnamem),size=2)+
  theme_bw()+geom_hline(yintercept=38.8,size=1)+xlab("Years")+
  ylab("UTCI (°C)")+theme(text = element_text(size = 20))+
  facet_wrap(~index)

ggplot(data=df245m)+geom_point(aes(x=years,y=models,col=modelnamem),size=2)+
  theme_bw()+geom_hline(yintercept=38.8,size=1)+xlab("Years")+
  ylab("UTCI (°C)")+theme(text = element_text(size = 20))+
  facet_wrap(~index)

ggplot(data=df585m)+geom_point(aes(x=years,y=models,col=modelnamem),size=2)+
  theme_bw()+geom_hline(yintercept=38.8,size=1)+xlab("Years")+
  ylab("UTCI (°C)")+theme(text = element_text(size = 20))+
  facet_wrap(~index)

#animation
pd254mun = pd245m > 38.8


saveGIF({
  for(i in c(1:80)){
    year <- 2019
    l <- levelplot(py245[[i]]-273.15,main=paste0("RCP4.5 Projection for Year ",year+i),at=cutpts, col.regions=(my_palette), margin=FALSE)
    l<- l + latticeExtra::layer(sp.lines(countries, col="black", lwd=1))
    plot(l)
  }
}, interval=0.2, movie.name="pd245long.gif")
