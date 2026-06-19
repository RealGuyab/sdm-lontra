#cargar paquetes de trabajo
requiredPackages <- c(
  "here", 
  "rstudioapi", 
  "ggplot2", 
  "tidyverse", 
  "rgdal", 
  "raster", 
  "maps", 
  "RColorBrewer", 
  "scam", 
  "ggpubr"
)

install_load_function <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[,
                                                 "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

install_load_function(requiredPackages)

#Set Plot parameter
theme_set(theme_bw(base_size = 16))


#Preparando data ambiental
load(here::here("gam_modelling", "spatial", "study_area.RData"))
mylayers <- stack(here::here("gam_modelling", "env",
                             "mylayers.tif"))
env_dataframe <- raster::as.data.frame(mylayers,
                                       xy = TRUE)
#renombrar capas 
names(env_dataframe) <- c("x", "y", "BO22_damean",
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
                                     "BO2_salinitymean_ss")

#Establcer proyecciones
# Load SC-GAM model
load(here::here("gam_modelling", "outpout_modelling", "selected_model.Rdata"))

# Predicting
predict <- predict(model, newdata = env_dataframe,
                   type = "response", se.fit = T)

env_dataframe$fit <- predict$fit
env_dataframe$se.fit <- predict$se.fit

save(env_dataframe, file = "gam_modelling/results/projection.Rdata")

#Mapeo
# Load PA data
load(here::here ("gam_modelling/outpout_modelling/PAdata_with_env7_final.RData"))


proj_map <-ggplot()+
  geom_raster(data=subset(env_dataframe),
              aes(x,y,fill=fit)) +
  scale_fill_gradient2(low="blue", 
                       mid="orange",
                       high="red",
                       midpoint = 0.5,
                       limits = c(0,1)) +
  ggtitle("Occurrence probabilty Lontra Felina")+ 
  geom_point(data=subset(data,occurrenceStatus==1),
             aes(LON,LAT),
             col=1,
             size=0.3) +
  theme_pubclean(base_size = 14)+
  theme(panel.background = element_blank(),
        plot.title = element_text(face = "italic"), 
        #text = element_text(size = 14), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position="right") +
  labs(y="latitude", x = "longitude")

print(proj_map)
dir()
#Export as raster image
ggsave(filename = "Lontra_felina_img.tif",
       plot = proj_map, device = "tiff", path = ("gam_modelling/results"),
       height = 22, width = 30,
       units = "cm", dpi = 300)


#Export as raster tif
install.packages("terra")
library("terra")
pred_raster <- terra::rast(mylayers, nlyrs = 1)
pred_raster[] <- env_dataframe$fit
crs(pred_raster) <- "EPSG:4326"
writeRaster(pred_raster, here::here("gam_modelling/results/Lontra_felina.tif"),
            overwrite = TRUE)
