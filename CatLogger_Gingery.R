#################################################################################################################
#################################################################################################################
##
##
##                            Gingery CatLogger Code  
##
#################################################################################################################
#################################################################################################################


##############################################################################################
##
##                                   Collar LE 
##              
##                The following code does the following:
##          (1)  Calculates the minimum dis. between catloggers to a true point (Trimble) 
##                      This calculate a SD of error associated with the cat logger
##
##############################################################################################

#Packages
library(geosphere)
library(distances)

cat <- read.csv("CatLocs_TrimbleLocs.csv", header= TRUE)

long1 <- cat$clong
lat1 <- cat$clat
long2 <- cat$tlong
lat2 <- cat$tlat

# Calculate distance (m) between two points
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

cat$error <- earth.dist(long1, lat1, long2, lat2)
#Average linear error (31.58)
mean(cat$error)
# SD
sd(cat$error)

# Write .csv with error column for future use/ reference
write.csv(cat, file = "CatLocs_TrimbleLocs.csv", na="NA")

###############################################################################################
#
#                                           Fix Success Rate
#
#       We calculated Fix success rate by evaluating the number
#          of recorded positions that would be expected for each
#         neonate (givn start and end dates with hourly locations expected) and 
#         the number of locations that were actually recorded. This can
#          be repeated using the file "all.csv"
# 
#
###############################################################################################

#Packages
library(rgdal)

#Note:  Fix success rate calculated by evaluating 
#         the number of recorded positions that would be expected for each
#         neonate (given start and end dates with hourly locations expected) and 
#         the number of locations that were actually recorded. This can be repeated
#         using the file "all.csv"

allcats  <- read.csv("all.csv")
allcats$CatID <- as.factor(allcats$Collar)#make CatID a factor

fsr <- subset(allcats, Collar== "66362" ) # to determine how many observed locations we got on each neonate
# We obtained data from 20 collars

# CollarID : FSR : At least 1 week of data? (For use in land cover comparison)
#1 6275 : 0.88  : Yes
#2 6293 : 0.93 : Yes
#3 6301 : 0.76 : Eliminated from all analyses
#4 6593 : 0.89 : No
#5 6597 : 0.95 : No
#6 6602 : 0.92 : Yes
#7 6603 : 0.99 : Yes
#8 6619 : 0.92 : Yes
#9 6620 : 0.98 : Yes
#10 6621 : 0.94: Yes
#11 6623 : 0.67 : Eliminated from all analyses
#12 6631 : 0.93 : Yes
#13 6636 : 0.97: No
#14 6640 : 1 : No
#15 6643 : 0.85 : No
#16 6646 : 0.97 : Yes
#17 6652 : 1: No
#18 6655 : 0.87 : Yes
#19 6656 : 0.57 : Eliminated from all analyses
#20 66362 : 0.90 : No 
 
#transform coordinates from WGS84 to UTM18 so BBMM has meters to work with.
#create new spatial data frame
coords <- data.frame(x = allcats$Longitude, y = allcats$Latitude)
current <-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"

#Create spatial points data frame, projection defined same as above
coords.spdf <- SpatialPointsDataFrame(coords= coords, data=allcats, proj4string = CRS(current))

#projecting points into UTM18N
pa.crs <- CRS("+proj=utm +zone=18 +datum=WGS84")
cat.UTM <- spTransform(coords.spdf, CRS=pa.crs)

#creates new csv with necessary UTM coordinates as columns x and y
#The following csv also won't contain the 3 collars with FSR < 75%
write.csv(cat.UTM, file = "AllUseableCollars.csv", na="NA")


###############################################################################################
#
#                          Creating Traditional Home Range Representations
#
#                                 MCPs, buffers, maternal core ranges
#
#                           10 neonates were used to create MCPs and buffers:
#
#                    6275, 6293, 6602, 6603, 6619, 6620, 6621, 6631, 6646, 6655
#
###############################################################################################

#############                         Minimum convex polygons                               ###


#Selected 14 Random locations for each of the 10 neonates during the first week of life in ArcGIS
# The file WkOld_RandomPointsForMCP represents the 14 random locations per neonate
cats <-read.csv("WkOld_RandomPointsForMCP.csv",header=T) 

#Subsetting neonates
x = subset(cats, cats$Collar == "6619")
xy <- x[c(6:7)]

library(rgdal) #writeOGR
library(adehabitatHR) #mcp function
#library(maptools) # writePolyshape

#Creates class Spatial Points for 
xy <- SpatialPoints(xy)
#head(xysp)
plot(xy)
proj4string(xy) <- CRS("+proj=utm +zone=18 +ellps=WGS84")

## estimates the MCP
cp <- mcp(xy, percent=95)#(95% is the default)
plot(cp)
plot(xy, add=TRUE)
## MCP size
as.data.frame(cp)

#Writing the shapefile with writeOGR to include proj info
writeOGR(cp, dsn=".", layer ="6619MCP", driver = "ESRI Shapefile")

# MCPs were mapped in ArcGIS and MCPs in the southern unit were merged into one shapefile
#       and MCPs in the northern unit were merged in a second shapefile.


#############                         Buffers                                            ###

# Completed in ArcGIS

#############                         Maternal core ranges                               ###

# Brownian bridge movement models for maternal ranges
#
#     Four VIT caught neonates were monitored with Catloggers.
#
#     Code below creates a BBMM contour at a 50% and 95% contour for moms with a VIT monitored
#     neonate. Time frame = first 7 days post-capture. 
#

#Maternal BBMMs
library(BBMM)
#Making each collar recognizable
cats <-read.csv("VITmoms.csv",header=T) 
str(cats)
cats$CatID <- as.factor(cats$CollarID)#make CatID a factor
#Date/Time format
cats$Date2 <- as.Date(cats$UTC_Date, "%m/%d/%Y")
cats$No <- gsub('-', "",cats$Date2)
cats$NoT <- gsub(":", "", cats$UTC_Time)
cats$datetime <- paste(cats$No, cats$NoT) #Combines Two Columns

##MAKE SURE DATE IS THE SAME FORMAT for R to read DATE
cats$NewDate<-as.POSIXct(strptime(paste(cats$Date2, cats$UTC_Time), "%Y-%m-%d %H:%M:%S"))
cats<- cats[order(cats$CatID, cats$NewDate),]
timediff<- diff(cats$NewDate)  
cats <- cats[-1,] # remove first entry without any difference
cats$timelag <-as.numeric(abs(timediff))

# Created two BBMM contours (95% and 50% contour) for each mom , one mom at a time. 

#Desired date is specific to the first 7 days that their neonate was also monitored. The exact same dates were used
# for their neonates 7 day BBMM occurance distribution.
# To determine dates: use the file 'Catloggers" that lists the meta-data for each collar including the day
#     monitoring began for each neonate and who their mother was.
cat <- subset(cats, cats$CatID == "17181" & cats$timelag > 50 & cats$Date2 < as.Date("2017-06-05"))  

#cat<-subset(cats, cats$Date2 > as.Date("2017-05-9") & cats$Date2 < as.Date("2017-06-12") )  # another way to subset by dates
cat <- cat[-1,] #Removes the first record that has a wrong timelag
cat$CatID <- factor(cat$CatID)

BBMM = brownian.bridge(x=cat$Easting, y=cat$Northing, time.lag = cat$timelag,
                       location.error=10, cell.size=20)
#bbmm.summary(BBMM)
#I created the maternal home ranges at a 50% and 95% contour
contours = bbmm.contour(BBMM, levels= c(50, 95),
                        locations= NULL, plot=TRUE)
cont <- 50 #sets the contour that you want to plot later
# Create a shapefile with contour lines
library(maptools)
library(raster)
library(rgeos)
library(sp)

out <- data.frame(x = BBMM$x, y = BBMM$y, z = BBMM$probability)

# Make sure the data is properly projected
out.raster <- rasterFromXYZ(out, crs = CRS("+proj=utm +zone=18 +datum=NAD83"),
                            digits = 4) 

raster.contour <- rasterToContour(out.raster, levels = contours$Z) 

#checking
plot(raster.contour)
raster.contour <- spChFIDs(raster.contour,
                           paste(cont, "% Contour Line", sep="")) 

library(rgdal) # can't be loaded before CRS call in line above
writeOGR(obj = raster.contour, dsn=".", layer="Cont50_17181_1wk", driver="ESRI Shapefile")

# Home ranges were then viewed, mapped, and sized in ArcGIS
# VIT mom-fawn figure created in ARG GIS

###############################################################################################
#
#                          Creating BBMMs for neonates at 7-days post-monitoring                   #
#                       (same 10 neonates that were used to create MCPs and buffers)
#
#           The same 10 neonates were used here as were used to create MCPs and buffers:
#
#                    6275, 6293, 6602, 6603, 6619, 6620, 6621, 6631, 6646, 6655
#
###############################################################################################

#Packages
library(BBMM)

#Making each collar recognizable
cats <-read.csv("AllUseableCollars.csv",header=T) 
attach(cats)
cats$CatID <- as.factor(cats$Collar)#make CatID a factor

#Date/Time Format
# classifies Date, changes format
cats$Date2 <- as.Date(cats$Date, "%m/%d/%Y")
# gets rid of the hyphens from the classified date collumn
cats$No <- gsub('-', "",cats$Date2)
# gets rid of the ":" between hour min sex in the time collumn
cats$NoT <- gsub(":", "", cats$Time)
cats$datetime <- paste(cats$No, cats$NoT) #Combines Two Columns
cats$DT <- as.POSIXct(strptime(cats$datetime, format = '%Y %m %d %H %M'))
cats <- cats[order(cats$CatID, cats$DT),]
timediff<- diff(cats$DT)*60  
cats <- cats[-1,] # remove first entry without any difference
cats$timelag <-as.numeric(abs(timediff))

#Subset neonate and first 7 days of monitoring
cat <- subset(cats, cats$CatID == "6631")  
cat <-subset(cat, cat$Date2 < as.Date("2017-05-29"))  
cat <- cat[-1,] #Removes the first record that has a wrong timelag
cat$CatID <- factor(cat$CatID)

BBMM = brownian.bridge(x=cat$x, cat$y, time.lag=cats$timelag,
                       location.error=40, cell.size=10)
#bbmm.summary(BBMM)
#Viewing the points in  R with their contour
contours = bbmm.contour(BBMM, levels= c(95),
                        locations= NULL, plot= TRUE)

x <- points(x = cat$x, y = cat$y) # Points
cont <- 95 #sets the contour that you want to plot later

# Create a shapefile with contour lines
library(maptools)
library(raster)
library(rgeos)
library(sp)

out <- data.frame(x = BBMM$x, y = BBMM$y, z = BBMM$probability)

# Make sure the data is properly projected
out.raster <- rasterFromXYZ(out, crs = CRS("+proj=utm +zone=18 +datum=NAD83"),
                            digits = 5, res = c(10,10)) 

raster.contour <- rasterToContour(out.raster, levels = contours$Z) 

#checking
plot(raster.contour)

raster.contour <- spChFIDs(raster.contour,
                           paste(cont, "% Contour Line", sep="")) 
library(rgdal) # can't be loaded before CRS call in line above

writeOGR(obj = raster.contour, dsn=".", layer="BBMM_6631_1wk", driver="ESRI Shapefile")

# Home range polygon shapefiles were created in ArcGIS using the contours created here
#  and BBMMs in each study area (SS and NS) southern 
# unit were merged into individual shapefiles. 


######################################################################################################
#
#                Evaluating percent forest cover in each HR estimation (BBMMs, MCPs, buffers) 
#   
#   Step A:    Determining Percent forest cover in the BBMM contours
#
#
#   Step B:    Determining percent forest cover in MCPs
#
#
#   Step C:    Determining percent forest cover in Buffers
#
#
######################################################################################################


####################################################### 
###                                                   #
###                                                   #
###    A : Determining Percent forest in BBMMS       #
###                                                   #
###                                                   #
#######################################################

###                       What is showen is for the BBMMS in the Southern Study Area          ######
###                               However, it was also modified and used for the              ######
###                               northern study (layers were just replaced)                  ######

### Bring in Northern BBMMS  ##
library(raster) #raster function
library(rgeos) #readOGR
library(rgdal)

t.buff <- readOGR (dsn=".", layer= "South_BBMM")
south <- readOGR(dsn="." , layer = "SouthBlock")
#new.crs <-CRS("+proj=utm +zone=18 +ellps=WGS84")
new.crs <- CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
#proj4string(t.buff)
#proj4string(south)
t.buff <- spTransform(t.buff, CRS=new.crs)
south <- spTransform(south, CRS=new.crs)

plot(south)
plot(t.buff, add=TRUE)

#      IMPORTING CDL LAYER , RECLASSIFYING , AND CUTTING       #
#   Reclassifying CDL   
CDL <- raster("CDL_2015_42.tif")
## 0 is other (water and stuff)  , 1 is forest, and 2 is agriculture
mfor <- c(-Inf, 0.5, 0, 0.6, 6.5, 2, 6.6, 9.5, 0 , 9.6, 14.5, 2, 14.6, 20.5, 0, 20.6, 39.5, 2, 39.6, 40.5, 0, 40.6, 61.5, 2, 61.6, 62.5, 0, 62.6, 63.5, 1, 63.6, 65.5, 0, 65.6, 72.5, 2, 72.6, 73.5, 0, 73.6, 77.5, 2, 77.6, 140.5, 0, 140.6, 143.5, 1, 143.6, 175.5, 0, 175.6, 176.5, 2, 176.6, 203.5, 0, 203.6, 214.5, 2, 214.6, 215.5, 0, 216.5, 227.5, 2, 227.6, 228.5, 0, 228.6, 250.5, 2, 250.6, 253.5, 2, 253.6, 254.5, 2, 254.6, Inf, 0)
rclfor <- matrix(mfor, ncol=3, byrow=TRUE)
rcfor <- reclassify(CDL, rclfor)
plot(rcfor)

### Cutting CDL raster 
Albers.crs <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
south2 <- spTransform(south, CRS = Albers.crs)
foragclip2 <- crop(rcfor, south2)
plot(foragclip2)
plot(south2, add=TRUE)
buff.st <- spTransform(t.buff, CRS = Albers.crs)
plot(buff.st, add=TRUE)

#      CALCULATING LAND COVER PROPORTIONS             #

library(plyr) #ddply function
library(SDMTools) # ClassStat  function
library(maptools) #readShapeSpatial function

indatacs2 <- buff.st
innamescs2 <- unique(buff.st@data$Collar)
outnamescs2 <- innamescs2
# set up output table
#output <- as.data.frame(matrix(0,nrow=length(innames),ncol=38))
#makes shapefile for each county and makes a list in a txt file
# begin loop to create separate county shapefiles
for (i in 1:length(innamescs2)){
  data <- indatacs2[which(indatacs2$Collar==innamescs2[i]),]
  if(dim(data)[1] != 0){
    writePolyShape(data,fn=paste(outnamescs2[i],sep="/"),factor2char=TRUE)
    write.table(innamescs2, "Listcs2.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)
  }
}
#reading in the table of shapefile names
#Read in a list of shapefiles files from above
Listshpscs2<-read.table("Listcs2.txt",sep="\t",header=F)
#colnames(Listshps) <- c("id")
Listshpscs2
#take shapefile and overlay on raster
shape2 <- function(Listshpscs2) {
  file <- as.character(Listshpscs2[1,])
  shp <- readShapeSpatial(file)
  mask <- mask(foragclip2,shp)
  ### Calculate the Class statistics in each county
  cl.data <- ClassStat(mask)
}
resultscs2 <- ddply(Listshpscs2, 1, shape2)
resultscs2
write.csv(resultscs2, file="South_BBMM.csv")

####################################################### 
###                                                   #
###                                                   #
###    B: Determining Percent forest in MCPs        #
###                                                   #
###                                                   #
#######################################################

###                       What is showen is for the MCPs in the Southern Study Area          ######
###                               However, it was also modified and used for the              ######
###                               Northern study (layers were just replaced)                  ######

### Bring in MCPs  
library(raster) #raster function
library(rgeos) #readOGR
library(rgdal)

t.buff <- readOGR (dsn=".", layer= "South_MCPs")
south <- readOGR(dsn="." , layer = "SouthBlock")
#new.crs <-CRS("+proj=utm +zone=18 +ellps=WGS84")
new.crs <- CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
#proj4string(t.buff)
#proj4string(south)
t.buff <- spTransform(t.buff, CRS=new.crs)
south <- spTransform(south, CRS=new.crs)
plot(south)
plot(t.buff, add=TRUE)

#      IMPORTING CDL LAYER , RECLASSIFYING , AND CUTTING       

#   Reclassifying CDL          
CDL <- raster("CDL_2015_42.tif")
## 0 is other (water and stuff)  , 1 is forest, and 2 is agriculture
mfor <- c(-Inf, 0.5, 0, 0.6, 6.5, 2, 6.6, 9.5, 0 , 9.6, 14.5, 2, 14.6, 20.5, 0, 20.6, 39.5, 2, 39.6, 40.5, 0, 40.6, 61.5, 2, 61.6, 62.5, 0, 62.6, 63.5, 1, 63.6, 65.5, 0, 65.6, 72.5, 2, 72.6, 73.5, 0, 73.6, 77.5, 2, 77.6, 140.5, 0, 140.6, 143.5, 1, 143.6, 175.5, 0, 175.6, 176.5, 2, 176.6, 203.5, 0, 203.6, 214.5, 2, 214.6, 215.5, 0, 216.5, 227.5, 2, 227.6, 228.5, 0, 228.6, 250.5, 2, 250.6, 253.5, 2, 253.6, 254.5, 2, 254.6, Inf, 0)
rclfor <- matrix(mfor, ncol=3, byrow=TRUE)
rcfor <- reclassify(CDL, rclfor)
plot(rcfor)

### Cutting CDL raster
Albers.crs <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
south2 <- spTransform(south, CRS = Albers.crs)

foragclip2 <- crop(rcfor, south2)
plot(foragclip2)
plot(south2, add=TRUE)
buff.st <- spTransform(t.buff, CRS = Albers.crs)
plot(buff.st, add=TRUE)

#      CALCULATING LAND COVER PROPORTIONS             
library(plyr) #ddply function
library(SDMTools) # ClassStat  function
library(maptools) #readShapeSpatial function

indatacs2 <- buff.st
innamescs2 <- unique(buff.st@data$Collar)
outnamescs2 <- innamescs2
# set up output table
#output <- as.data.frame(matrix(0,nrow=length(innames),ncol=38))
#makes shapefile for each county and makes a list in a txt file
# begin loop to create separate county shapefiles
for (i in 1:length(innamescs2)){
  data <- indatacs2[which(indatacs2$Collar==innamescs2[i]),]
  if(dim(data)[1] != 0){
    writePolyShape(data,fn=paste(outnamescs2[i],sep="/"),factor2char=TRUE)
    write.table(innamescs2, "Listcs2.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)
  }
}
#reading in the table of shapefile names
#Read in a list of shapefiles files from above
Listshpscs2<-read.table("Listcs2.txt",sep="\t",header=F)
#colnames(Listshps) <- c("id")
Listshpscs2
#take shapefile and overlay on raster
shape2 <- function(Listshpscs2) {
  file <- as.character(Listshpscs2[1,])
  shp <- readShapeSpatial(file)
  mask <- mask(foragclip2,shp)
  ### Calculate the Class statistics in each county
  cl.data <- ClassStat(mask)
}
resultscs2 <- ddply(Listshpscs2, 1, shape2)
resultscs2
write.csv(resultscs2, file="South_MCP.csv")

####################################################### 
###
###
###    C2 : Determining Percent forest in Buffers     #
###
###
#######################################################


###                       What is showen is for the Buffers in the Southern Study Area          ######
###                               However, it was also modified and used for the              ######
###                               northern study (layers were just replaced)                  ######


### Bring in  Buffers
library(raster) #raster function
library(rgeos) #readOGR
library(rgdal)

t.buff <- readOGR (dsn=".", layer= "South_CapBuffs")
south <- readOGR(dsn="." , layer = "SouthBlock")
#new.crs <-CRS("+proj=utm +zone=18 +ellps=WGS84")
new.crs <- CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
#proj4string(t.buff)
#proj4string(south)
t.buff <- spTransform(t.buff, CRS=new.crs)
south <- spTransform(south, CRS=new.crs)
plot(south)
plot(t.buff, add=TRUE)

#      IMPORTING CDL LAYER , RECLASSIFYING , AND CUTTING    

#   Reclassifying CDL     

CDL <- raster("CDL_2015_42.tif")
## 0 is other (water and stuff)  , 1 is forest, and 2 is agriculture
mfor <- c(-Inf, 0.5, 0, 0.6, 6.5, 2, 6.6, 9.5, 0 , 9.6, 14.5, 2, 14.6, 20.5, 0, 20.6, 39.5, 2, 39.6, 40.5, 0, 40.6, 61.5, 2, 61.6, 62.5, 0, 62.6, 63.5, 1, 63.6, 65.5, 0, 65.6, 72.5, 2, 72.6, 73.5, 0, 73.6, 77.5, 2, 77.6, 140.5, 0, 140.6, 143.5, 1, 143.6, 175.5, 0, 175.6, 176.5, 2, 176.6, 203.5, 0, 203.6, 214.5, 2, 214.6, 215.5, 0, 216.5, 227.5, 2, 227.6, 228.5, 0, 228.6, 250.5, 2, 250.6, 253.5, 2, 253.6, 254.5, 2, 254.6, Inf, 0)
rclfor <- matrix(mfor, ncol=3, byrow=TRUE)
rcfor <- reclassify(CDL, rclfor)
plot(rcfor)

### Cutting CDL raster for Northern study area Buffers
Albers.crs <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
south2 <- spTransform(south, CRS = Albers.crs)
foragclip2 <- crop(rcfor, south2)
plot(foragclip2)
plot(south2, add=TRUE)
buff.st <- spTransform(t.buff, CRS = Albers.crs)
plot(buff.st, add=TRUE)

#      CALCULATING LAND COVER PROPORTIONS                                                                  # 
library(plyr) #ddply function
library(SDMTools) # ClassStat  function
library(maptools) #readShapeSpatial function

indatacs2 <- buff.st
innamescs2 <- unique(buff.st@data$Collar)
outnamescs2 <- innamescs2
# set up output table
#output <- as.data.frame(matrix(0,nrow=length(innames),ncol=38))
#makes shapefile for each county and makes a list in a txt file
# begin loop to create separate county shapefiles
for (i in 1:length(innamescs2)){
  data <- indatacs2[which(indatacs2$Collar==innamescs2[i]),]
  if(dim(data)[1] != 0){
    writePolyShape(data,fn=paste(outnamescs2[i],sep="/"),factor2char=TRUE)
    write.table(innamescs2, "Listcs2.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)
  }
}
#reading in the table of shapefile names
#Read in a list of shapefiles files from above
Listshpscs2<-read.table("Listcs2.txt",sep="\t",header=F)
#colnames(Listshps) <- c("id")
Listshpscs2
#take shapefile and overlay on raster
shape2 <- function(Listshpscs2) {
  file <- as.character(Listshpscs2[1,])
  shp <- readShapeSpatial(file)
  mask <- mask(foragclip2,shp)
  ### Calculate the Class statistics in each county
  cl.data <- ClassStat(mask)
}
resultscs2 <- ddply(Listshpscs2, 1, shape2)
resultscs2
write.csv(resultscs2, file="South_Buffer.csv")

###############################################################################
#
#                 Distance Moved from Capture Location
#
#               Euclidean distance between each position
#     and the capture location for every day.
#     
#
#     We used data from 15 neonates.
#         - neonates must have at least 2 days of location data
#
#           Neonates used: 6293, 6275, 6621, 6631, 6640 , 
#           6646,6602, 6603,6619, 6620, 6655, 6597, 6636, 6643, 6652
#
################################################################################

setwd("C:/Users/tessg/OneDrive - The Pennsylvania State University/Tess' Working Folder/Cat Collars/AA_DatNewNew/Distance")
data1 <- read.csv("AllUseableCollars.csv", header= TRUE)

#Assigning collar id's as sequential numbers
data1$Collar <- as.factor(data1$Collar)
data2<- transform(data1$Collar, id=match(data1$Collar, unique(data1$Collar)))
data1$group <- as.factor(data2$id)

#The following code takes the difference between the first x and y 
#       (the capture location) and the subsequent x and y positions

#Did this one collar at a time 
data <- subset(data1, data1$Collar == 6597)

#One way
for ( i in 2:(nrow(data))) {
  #   data$Dist[i]  
  data$Xdist[i] <- data$x[i] - data$x[1]
  data$Ydist[i] <- data$y[i] - data$y[1]
  data$Dist[i]   <- as.integer(sqrt(data$Xdist[i]^2+data$Ydist[i]^2))
} 

#A second way
data$Ydist  <- data$y[1] - data$y[1]
for (i in 2:(nrow(data))) {
  
  data$Ydist[i]  <- data$y[i] - data$y[1]
}

data$Xdist  <- data$x[1] - data$x[1]
for (i in 2:(nrow(data))) {
  
  data$Xdist[i]  <- data$x[i] - data$x[1]
}

data$Dist[1]  <- as.integer(sqrt(data$Xdist[1]^2+data$Ydist[1]^2))
for (i in 2:(nrow(data))) {
  
  data$Dist[i]  <- as.integer(sqrt(data$Xdist[i]^2+data$Ydist[i]^2))
}

#Write one .csv for each neonate
write.csv(data, file = "Dist6597.csv", na="NA")

#getting 1 .csv that has all neonates and the distance moved each hour

files <- list.files(pattern = '\\.csv')
tables <- lapply(files, read.csv, header = TRUE)

#following line won't work if number of columns isn't exact in every sheet
combined.df <- do.call(rbind , tables)
write.csv(combined.df,"all.distance.csv", 
          row.names=TRUE,
          col.names=TRUE)

# Using the .csv created above:
#Calculating the mean distance each DAY using
#   the distance moved per hour above 
#Did this for each individual fawn
cats <-read.csv("all.distance.csv",header=T) 
cats$CatID <- as.factor(cats$Collar)#make CatID a factor
# classifies Date as a date 
cats$Date2 <- as.Date(cats$Date, "%m/%d/%Y")
#I ran this line and the next line 1 time for each fawn
cat1 <- subset(cats, cats$CatID == "6652")
#x = unique(as.Date(cat1$Date2))

library(plyr)

y = ddply(cat1, .(Date2), summarize, 
          mean_Dist = mean(Dist),
          sd = sd(Dist))
y <- y[-1,] # remove first entry without any difference
y$sex <- (y$sex ="female")
y$day <- (seq.int(nrow(y)))
y$Collar <- (y$Collar = "6652")
write.csv(y, file = "6652dailymean.csv")

#Creating a .csv that details the mean distance each day
files <- list.files(pattern = '\\.csv')
tables <- lapply(files, read.csv, header = TRUE)
combined.df <- do.call(rbind , tables)
write.csv(combined.df,"indiv.means.csv", 
          row.names=TRUE,
          col.names=TRUE)
###############################################################################
#
#                 Home range size (BBMM 95% contour) over 30 days
#
#     We used data from 11 neonates.
#         - neonates must have at least 4 days of location data
#
#   Neonates used: 6275, 6293, 6597, 6602, 6603, 6619, 6620, 6621, 6631, 6646, 6655
#
################################################################################

#Packages
library(BBMM)
#Making each collar recognizable
cats <-read.csv("AllUseableCollars.csv",header=T) 
attach(cats)
cats$CatID <- as.factor(cats$Collar)#make CatID a factor
#Getting the date and time into the right format
#Modified from David Walter Spatial Ecology code 
# classifies Date as a date and changes it's format
cats$Date2 <- as.Date(cats$Date, "%m/%d/%Y")
# gets rid of the hyphens from the classified date collumn
cats$No <- gsub('-', "",cats$Date2)
# gets rid of the ":" between hour min sex in the time collumn
cats$NoT <- gsub(":", "", cats$Time)
cats$datetime <- paste(cats$No, cats$NoT) #Combines Two Columns
cats$DT <- as.POSIXct(strptime(cats$datetime, format = '%Y %m %d %H %M'))
cats <- cats[order(cats$CatID, cats$DT),]
timediff<- diff(cats$DT)*60  #In seconds
cats <- cats[-1,] # remove first entry without any difference
cats$timelag <-as.numeric(abs(timediff))

#  Each BBMM was run individually.
#        For example, so I ran the following BBMM code approximately 15 times for each neonate
#        to get a BBMM at 2 days, a BBMM at 4 days, a BBMM at 6 days.. until 30 days subsetting out
#         particular dates that are specific to that neonate and the date it began to be monitored. 
cat <- subset(cats, cats$CatID == "6619")  #Subsetting 1 individual
cat <-subset(cat, cat$Date2 < as.Date("2017-06-08"))  # Filtering needed dates
cat <- cat[-1,] #Removes the first record 
cat$CatID <- factor(cat$CatID)
#BBMM
BBMM = brownian.bridge(x=cat$x, cat$y, time.lag=cats$timelag,
                       location.error=40, cell.size=10)
#Viewing the points in  R with their contour
contours = bbmm.contour(BBMM, levels= c(95),
                        locations= NULL, plot= TRUE)

x <- points(x = cat$x, y = cat$y) # Points 
cont <- 95 #sets the contour
# Create a shapefile with contour lines
library(maptools)
library(raster)
library(rgeos)
library(sp)

out <- data.frame(x = BBMM$x, y = BBMM$y, z = BBMM$probability)
# Projection
out.raster <- rasterFromXYZ(out, crs = CRS("+proj=utm +zone=18 +datum=NAD83"),
                            digits = 5, res = c(10,10)) 
raster.contour <- rasterToContour(out.raster, levels = contours$Z) 
#checking
plot(raster.contour)
raster.contour <- spChFIDs(raster.contour,
                           paste(cont, "% Contour Line", sep="")) 
library(rgdal) # can't be loaded before CRS call in line above
writeOGR(obj = raster.contour, dsn=".", layer="BBMM_6619_2_day", driver="ESRI Shapefile")

# Home ranges were mapped and sized in ArcGIS where resutling BBMM sizes were recorded
#    in an excel sheet named "HrSizesDuringFirstMonth" used to model HR sizes and plot


######################################################################################################################
#
#                           Temporal Trends over 30 days
#  
#                          Movement distance and HR size 
#                      
#     Step 1:      2 Random effect models where individual intercepts and slope vary
#
#     Step 2:      Plot- Figure 1 in Gingery et al. 
#       
######################################################################################################################


########################          Step 1             ################################
#Movement Model
dist <-read.csv("indiv.means.csv",header=T) 
dist = subset(dist, dist$day < 31)

a <- dist$day 
b <- dist$mean_Dist
ID <- dist$Collar
y <- dist$mean_Dist
x <- dist$day

library(lme4)
library(effects)
#Random effect model with varying intercept and slop
m1 <- lmer(data=dist, y ~ x + (1+ x|ID), REML = F)
summary(m1)
#Indiv Intercepts- movement 
coef(m1)$ID

#Home range size model
move <-read.csv("HrSizesDuringFirstMonth.csv",header=T)
a.2 <- move$Day 
b.2 <- move$Hrsize_ha
ID.2 <- move$Collar
y.2 <- move$Hrsize_ha
x.2 <- move$Day

library(lme4)
library(effects)
#Random effect model with varying intercept and slop
m1.2 <- lmer(data=move, y.2 ~ x.2 + (1+ x.2|ID.2), REML = F)
summary(m1.2)
#Neonate intercept- HR
coef(m1.2)$ID.2

##########################              Step 2 - figure      ###################

#Setting the upper and lower limits - movement
CollarID <- as.numeric(dist$Collar, dist)
hr <- as.numeric(dist$mean_Dist, dist)
#takes the mean of HRSize by each collar id
l2ptp <- by(hr , INDICES= CollarID, FUN= mean, na.rm=T)
l2ptp <- l2ptp[!is.na(l2ptp)]
x <- l2ptp # rename for plotting
x <- as.numeric(l2ptp)
M2.coeftp <- coef(m1)

#Creates my x value for movement graph which is age in two day intervals
xtp <-seq(min(dist$day),max(dist$day), length=7)
# takes the intercept, slope for the overall random effect model and times it by day
predtp <- fixef(m1)["(Intercept)"] + fixef(m1)["x"] * xtp

#slopes- movement
tphupred <- matrix(NA, nrow=length(xtp), ncol=12)
for(j in 1:12){
  tphupred[,j] <- M2.coeftp[[1]][j,1] + M2.coeftp[[1]][j,2]*xtp
}

#Setting the upper and lower limits - HR
#Makes collar id and HR numeric
CollarID.2 <- as.numeric(move$Collar, move)
hr.2 <- as.numeric(move$Hrsize_ha, move)
#takes the mean of HRSize by each collar id
l2ptp.2 <- by(hr.2 , INDICES= CollarID.2, FUN= mean, na.rm=T)
l2ptp.2 <- l2ptp.2[!is.na(l2ptp.2)]
x.2 <- l2ptp.2 # rename for plotting
x.2 <- as.numeric(l2ptp.2)
M2.coeftp.2 <- coef(m1.2)

#Creates my x value for HR graph - age in two day intervals
xtp.2 <-seq(min(move$Day),max(move$Day), length=15)
# takes the intercept, slope for the overall random effect model and times it by day
predtp.2 <- fixef(m1.2)["(Intercept)"] + fixef(m1.2)["x.2"] * xtp.2

# HU-specific slopes- HR
tphupred.2 <- matrix(NA, nrow=length(xtp.2), ncol=12)
for(j in 1:12){
  tphupred.2[,j] <- M2.coeftp.2[[1]][j,1] + M2.coeftp.2[[1]][j,2]*xtp.2
}

#starting graphing parameters 
res <- 6
name_figure <- "Test_2.png"
png(filename = name_figure, height = 900*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)
size.labels = 1
size.text = 1
axissize <- 1
x.label = 'Days after capture'
nf <- layout(matrix(c(1:2),nrow=2,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1.5,0.80,1.5,1), oma=c(2,2,0.3,0.1) )

#Movement Plot
plot(b ~ a,data=dist, axes=F,pch=16, ylim=c(0,1000), xlim=c(0,30))
axis(side=1,cex.axis=size.text,pos=0, tck=-0.01, mgp=c(0,0.5,0),at=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30),labels=T )
axis(side=2,cex.axis=size.text,pos=0, las=1, mgp=c(0,0.5,0),tck=-0.01,at=c( 0, 200, 400, 600, 800, 1000), labels=c( 0, 200, 400, 600, 800, 1000))
text(2 ,1000, label= "A", font = 2, cex = 1.7 )
# Add random slopes- Movement
for(i in 1:12){
  points(xtp, tphupred[,i], lwd=1, col='grey', lty=1, type='l')
}
# Average effect
points(xtp, predtp, lwd=4, col='black', lty=1, type='l')

#HR Plot
plot(b.2 ~ a.2,data=move, axes=F,xlim = c(0,30), ylim=c(0,60), pch=16)
axis(side=1,cex.axis=size.text,pos= 0,  tck=-0.01, mgp=c(0,0.5,0),at=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30),labels=T )
axis(side=2,cex.axis=size.text,pos=0, las=1, mgp=c(0,0.5,0),tck=-0.01,at=c(0, 10, 20, 30, 40, 50, 60), labels=c( 0, 10, 20, 30, 40, 50, 60))

# Add random slopes- HR
for(i in 1:12){
  points(xtp.2, tphupred.2[,i], lwd=1, col='grey', lty=1, type='l')
}
# Average effect
points(xtp.2, predtp.2, lwd=4, col='black', lty=1, type='l')

# Labels and legend

mtext(expression(paste('Days after capture')), line = 0.1, side = 1, cex = size.text, outer=T)
mtext(expression(paste('Average distance from capture location (m)')), line = 0.1, side = 2, cex = size.text, outer=T, at=0.78)
mtext(expression(paste('Home range size (ha)')), line = 0.1, side = 2, cex = size.text, outer=T, at=0.3)
legend(12, 55, legend=c("Mean relationship", "Individual neonates"), col=c("black", "grey"), lty= 1:1, lwd= 4:1, box.lty=0)
text(2 ,60, label= "B", font = 2, cex = 1.7 )

# Ending Graph
par(def.par)
dev.off()

################################################################################################################
#                                                                                                              #
#                                                 End                                                          #
#                                                                                                              #
################################################################################################################



