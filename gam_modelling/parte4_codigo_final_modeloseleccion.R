##Model Fit
requiredPackages <- c(
  #GENERAL USE LIBRARIES --------#
  "here", # Library for reproducible workflow
  "rstudioapi",  # Library for reproducible workflow
  "stringr",
  "RColorBrewer",  
  "ggplot2",
  "dplyr",
  
  #SPATIAL DATA --------#
  "rgdal",
  "fields",
  "maps" ,
  "raster",
  
  #MODEL FIT --------#
  "scam",
  "plotmo",
  "SDMTools",
  "pkgbuild",
  "dismo"
)

install_load_function <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

install_load_function(requiredPackages)

find_rtools()
install.packages("remotes")
remotes::install_version("SDMTools", version = "1.1-221.2")

# General settings for ggplot (black-white background, larger base_size)
theme_set(theme_bw(base_size = 16))

#Configuracion del directorio de trabajo y carga de datos de Psudoausencia y datos ambientales
load(here::here ("gam_modelling/outpout_modelling/PAdata_with_env7_final.Rdata"))

###Univariable Models
##BO2_curvelmean_ss
model_curvelmean_ss <- scam (occurrenceStatus ~  s(BO2_curvelmean_ss, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_curvelmean_ss)
plotmo(model_curvelmean_ss, level = 0.95, pt.col=8)

##BO2_ironmean_ss
model_iron <- scam (occurrenceStatus ~  s(BO2_ironmean_ss, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_iron)
plotmo(model_iron, level = 0.95, pt.col=8)

##BO22_bathy
model_bathy <- scam (occurrenceStatus ~  s(BO_bathymean, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_bathy)
plotmo(model_bathy,level = 0.95, pt.col=8)

##BO22_parmean
model_parmean <- scam (occurrenceStatus ~  s(BO22_parmean, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_parmean)
plotmo(model_parmean,level = 0.95, pt.col=8)

##BO_salinity
model_salinity <- scam (occurrenceStatus ~  s(BO_salinity, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_salinity)
plotmo(model_salinity,level = 0.95, pt.col=8)

##BO_sstmean
model_sstmean <- scam (occurrenceStatus ~  s(BO_sstmean, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_sstmean)
plotmo(model_sstmean, level = 0.95, pt.col=8)

##BO22_ph
model_ph <- scam (occurrenceStatus ~  s(BO22_ph, k=8,bs="cv"), family=binomial(link="logit"), data=data)
summary(model_ph)
plotmo(model_ph, level = 0.95, pt.col=8)

####Modelo Multivariado
model_large <- scam (occurrenceStatus ~  
                s(BO22_parmean, k=8,bs="cv") +
               s(BO2_curvelmean_ss, k=8,bs="cv") +
               s(BO2_ironmean_ss, k=8,bs="cv") +
               s(BO_bathymean, k=8,bs="cv") +
               s(BO_salinity, k=8,bs="cv") +
               s(BO_sstmean, k=8,bs="cv") +
               s(BO22_ph, k=8,bs="cv"),
               family=binomial(link="logit"), data=data)
summary(model_large)
plotmo(model_large,level = 0.95, pt.col=8)


####Multivariate Model without less relevant variables
model <- scam (occurrenceStatus ~  
                 s(BO22_parmean, k=8,bs="cv") +
                 #s(BO2_curvelmean_ss, k=8,bs="cv") +
                 #s(BO2_ironmean_ss, k=8,bs="cv") +
                 s(BO_bathymean, k=8,bs="cv") +
                 #s(BO_salinity, k=8,bs="cv") +
                 s(BO_sstmean, k=8,bs="cv"),
               #s(BO22_ph, k=8,bs="cv"),
               family=binomial(link="logit"), data=data)
summary(model)
plotmo(model,level = 0.95, pt.col=8)

#analisis del modelo
old.par <- par(mfrow = c(2, 2))
scam.check(model)

install.packages("pROC")
library(pROC)
roc_curve <- roc(data$occurrenceStatus, fitted(model))
plot(roc_curve, col="blue")
auc(roc_curve)

save(list="model", file="gam_modelling/outpout_modelling/selected_model.RData")

