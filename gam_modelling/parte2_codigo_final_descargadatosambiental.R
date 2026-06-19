requiredPackages <- c(
  #GENERAL USE LIBRARIES --------#
  "here", # Library for reproducible workflow
  "rstudioapi",  # Library for reproducible workflow
  "maptools", #plotting world map
  "ggplot2", #for plotting
  "knitr",  # format tables
  "kableExtra", # format tables
  "raster", # to work with spatial data
  "dplyr",  
  #DOWNLOAD FROM PUBLIC REPOSITORIES --------#
  "sdmpredictors" #to access Bio-ORACLE dataset
)
install_load_function <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

install_load_function(requiredPackages)

#####seleccion de variables de las fuentes de datos Bioracle, WorlClim, 
bioracle_variables <- list_layers("Bio-ORACLE")


kable(bioracle_variables)%>% 
  kable_styling("striped") %>% 
  scroll_box(height="600px", width = "100%")

target <- c("BO22_damean",
            "BO22_parmean",
            "BO22_ph",
            "BO_chlomean",
            "BO_dissox",
            "BO_salinity",
            "BO_sstmean",
            "BO_bathymean",
            "BO_curvelmean_bdmean",
            "BO_curvelmean_bdmin",
            "BO2_curvelmean_ss",
            "BO2_dissoxmean_ss",
            "BO2_ironmean_ss",
            "BO2_salinitymean_ss")

# Extrat details from the list
myvars <- bioracle_variables %>% 
  dplyr::filter (bioracle_variables$layer_code %in% target)

myvars$name

# Download layers
myBioracle.layers <- load_layers(c("BO22_damean",
                                   "BO22_parmean",
                                   "BO22_ph",
                                   "BO_chlomean",
                                   "BO_dissox",
                                   "BO_salinity",
                                   "BO_sstmean",
                                   "BO_bathymean",
                                   "BO2_curvelmean_ss",
                                   "BO2_dissoxmean_ss",
                                   "BO2_ironmean_ss",
                                   "BO2_salinitymean_ss"),datadir = "data/spatial/Bioracle_variables") 

save (list = "myBioracle.layers",
      file = "data/spatial/myBioracle.layers.RData")

#code to crop the raster for a specific area or study area
load("data/spatial/study_area.RData")
mylayers <- crop(myBioracle.layers, extent(-85,-50,
                                           -60,0))
plot(mylayers)

#create raster layers for the variables used for the model 
writeRaster(mylayers, filename="gam_modelling/env/mylayers.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
