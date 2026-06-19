requiredPackages <- c(
  #GENERAL USE LIBRARIES --------#
  "here", # Library for reproducible workflow
  "rstudioapi",  # Library for reproducible workflow
  "maptools", #plotting world map
  "ggplot2", #for plotting
  
  #Download presence data--------#
  "robis", # Specific library to get the occurrence data
  "rgbif",# Specific library to get the occurrence data
  "CoordinateCleaner", #to remove outlier
  "rgdal", # to work with Spatial data
  "sf", # to work with spatial data (shapefiles)
  "data.table", #for reading data,
  "dplyr", #for reading data,
  "tidyr", #for reading data
  "marmap", #bathymetry getNOAA.bathy remotes::install_github("ericpante/marmap")
  
  #Create pseudo-absence data--------#
  "tidyverse", 
  "scales",
  "ggridges",
  "maps",     # some basic country maps
  "mapdata",   # higher resolution maps
  "mapproj",
  "mapplots",   # ICES rectangles
  "gridExtra",
  "lubridate",
  "raster" # to work with Spatial data
)

#funcion de instalacion de paquetes de datos
install_load_function <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

install_load_function(requiredPackages)

###DataDownload

# url where FAO shapfile is stored
url<-"https://www.fao.org/fishery/geoserver/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=fifao:FAO_MAJOR&outputFormat=SHAPE-ZIP"

# Download file
download.file(url,"gam_modelling/data/FAO_AREA.zip",mode="wb")

# Unzip downloaded file
unzip("gam_modelling/data/FAO_AREA.zip",
      exdir="gam_modelling/spatial")

# Load FAO (spatial multipolygon)
FAO<- st_read(file.path("gam_modelling/spatial", "FAO_MAJOR.shp"))

# Select Atlantic Ocean FAO Area 
FAO_Pa <- FAO[FAO$OCEAN=="Pacific",]

# Transform to spatial polygons dataframe
study_area<- sf:::as_Spatial(FAO_Pa)
save(study_area, file=file.path("gam_modelling/spatial",
                                file="study_area.RData"))
plot(study_area)
rm(FAO,FAO_Pa)

#Get data from GBIF
mydata.gbif<-occ_data(scientificName="Lontra felina (Molina, 1782)", hasCoordinate = TRUE, limit=100000)$data
save(mydata.gbif,file=file.path("gam_modelling/occurrences",
                                file="mydata_gbif.RData"))
load(here::here ("gam_modelling/occurrences", "mydata_gbif.RData"))

# Check names for GBIF data
names(mydata.gbif)

# Select columns of interest
mydata.gbif <- mydata.gbif %>%
  dplyr::select("acceptedScientificName",
                "verbatimLocality",
                "decimalLongitude",
                "decimalLatitude",
                "year",
                "month",
                "day",
                "eventDate")
mydata.gbif <- mydata.gbif %>% 
  dplyr::mutate(occurrenceStatus=1) %>%
  dplyr::mutate(DEP=NA) %>%
  dplyr::rename(Lugar= "verbatimLocality")

#Upload data collected from the field work
campo_nutrias<- read.csv("gam_modelling/occurrences/capa_final_wgs84.csv", header = T,sep=";", row.names=NULL)
save(campo_nutrias,file = file.path("gam_modelling/occurrences",
                                    file="campo.RData"))
load(here::here ("gam_modelling/occurrences", "campo.RData"))

#Join data from campo and GBIF 
campo_nutrias <- campo_nutrias %>% 
  dplyr::rename(DEP = "DEP..")
mydata.fus<-rbind(campo_nutrias,mydata.gbif)

# Remove unused files
rm(mydata.gbif, campo_nutrias)

# Give date format to eventDate and fill out month and date_year columns
mydata.fus$eventDate <- as.Date(mydata.fus$eventDate)
mydata.fus$date_year <- as.numeric(mydata.fus$year)
mydata.fus$month <- as.numeric(mydata.fus$month)

out.dist <- cc_outl(x=mydata.fus,
                    lon = "decimalLongitude", lat = "decimalLatitude",
                    species = "acceptedScientificName",
                    method="distance", tdi=1000, # distance method with tdi=1000km
                    thinning=T, thinning_res=0.5,
                    value="flagged") 

# Remove outliers from the data
mydata.fus <- mydata.fus[out.dist, ]

# First create a vector containing longitude, latitude and event date information
date <- cbind(mydata.fus$decimalLongitude,mydata.fus$decimalLatitude,mydata.fus$eventDate)

# Remove the duplicated records
mydata.fus<-mydata.fus[!duplicated(date),] 

# Remove unused files
rm(date)

# Assign coordinate format and projection to be able to use FAO Pacific as a mask
dat <- data.frame(cbind(mydata.fus$decimalLongitude,mydata.fus$decimalLatitude))
ptos<-as.data.table(dat,keep.columnnames=TRUE)

coordinates(ptos) <- ~ X1 + X2

# Assign projection
proj4string(ptos) <-proj4string(study_area)

# Select only occurrences from FAO Atlantic
match2<-data.frame(subset(mydata.fus,!is.na(over(ptos, study_area)[,1])))

# Extract the FAO area of each point
match3<-data.frame(subset(over(ptos, study_area), !is.na(over(ptos, study_area)[,1])))

# Create data frame with area, name, long, lat and year 
df0<-cbind(F_AREA=match3$F_AREA,match2)[,c("F_AREA","acceptedScientificName","decimalLongitude","decimalLatitude","year","occurrenceStatus")]

# Rename some columns
names(df0)[3:5]<-c("LON","LAT","YEAR")

#database construction for Maxent
occu_maxent <- df0%>%
  dplyr::select("acceptedScientificName",
                "LON",
                "LAT")

write.csv(occu_maxent, file= "gam_modelling/occurrences/occu_maxent.csv")

# Add bathymetry from NOAA
install.packages("ncdf4")
install.packages("rgdal")
install.packages("marmap")
library (ncdf4)
library(marmap)
options(timeout = 1000)
bathy <- marmap::getNOAA.bathy(lon1=-100,lon2=-70,lat1=-90,lat2=0, resolution = 1, keep=TRUE,
                               antimeridian=FALSE, path= "gam_modelling/spatial/bath")
save(bathy,file=file.path("gam_modelling/spatial/bath",
                          file="bathy.RData"))
load("gam_modelling/spatial/bath/bathy.RData")

df0$bathymetry <- get.depth(bathy, df0[,c("LON","LAT")], locator=F)$depth

# Remove unused files
rm(mydata.fus, match2, match3, dat, ptos)

#Create psudo-absence data ###

# Remove points in land
df0<-subset(df0,bathymetry<0)

# Select only years from 2000 to 2024
df0<- subset(df0, YEAR<=2024 & YEAR>=2000)

# Convert to spatial point data frame
df<-df0 ; coordinates(df)<- ~LON+LAT
crs(df)<-crs(study_area)

# Convert to sf
df.sf<-st_as_sf(df)
study_area <- crop(study_area, extent(-85,-50,
                                      -60,0))
study_area.sf<-st_union(st_as_sf(study_area))
st_write(study_area.sf, "gam_modelling/spatial/study_area.shp") #crea un shp en la carpeta destino

#plot data in SF
ggplot(study_area.sf) + 
  geom_sf() + 
  geom_sf(data=st_union(df.sf),
          size=1,
          alpha=0.5)

# Basic ggplot
global <- map_data("worldHires")

p0 <- ggplot() + 
  annotation_map(map=global, fill="grey")+
  geom_sf(data=study_area.sf,fill=5)

print(p0)

# Function to find your UTM. 
lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}

(EPSG_2_UTM <- lonlat2UTM(c(mean(df$LON), mean(df$LAT))))

# Transform study_area and data points to UTMs (in m)
aux <- st_transform(study_area.sf, EPSG_2_UTM)
df.sf.utm <- st_transform(df.sf, EPSG_2_UTM)

# Generate the pseudo-absence data frame
pseudo <- matrix(data=NA, nrow=dim(df0)[1], ncol=dim(df0)[2])
pseudo <- data.frame(pseudo)
names(pseudo) <- names(df0)

# Set the seed
set.seed(1)

# Sample from the defined area
rp.sf <- st_sample(aux, size=dim(df.sf.utm)[1], type="random") # randomly sample points

# Transform to lat and lon and extract coordinates as data.frame
rp.sf <- st_transform(rp.sf, 4326)
rp <- as.data.frame(st_coordinates(rp.sf)) 
pseudo$LON <- rp$X
pseudo$LAT <- rp$Y

# Complete the rest of columns
pseudo$acceptedScientificName <- df0$acceptedScientificName
pseudo$occurrenceStatus  <- 0

##grafico presencia ausencia pseudo
p0 +
  geom_sf(data=rp.sf, col=6, shape=4,size=0.5)+
  geom_sf(data=df.sf.utm, col=1, alpha=0.8,size=0.5)+
  ggtitle(unique(df$scientificName))

# Zoom
p0 +
  geom_sf(data=rp.sf, col=6, shape=4,size=1)+
  geom_sf(data=df.sf.utm, col=1, alpha=0.8,size=0.5)+
  coord_sf(xlim=c(-80,-70), ylim=c(-45,-10))+
  ggtitle(unique(df$acceptedScientificName))

# Join the two data sets 
PAdata <- rbind(df0, pseudo)[,c("acceptedScientificName","LON","LAT","YEAR","occurrenceStatus")]

# Save the final dataset of occurrence and pseudo-absence points
save(list=c("PAdata"),file=file.path("gam_modelling/outpout_modelling",file="PAdata.RData"))




