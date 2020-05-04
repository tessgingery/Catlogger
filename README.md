# Catlogger
R Code used for several analyses in a manuscript detailing white-tailed deer neonate home ranges

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
