###Preparation final dataset 
requiredPackages <- c(
  #GENERAL USE LIBRARIES --------#
  "here", # Library for reproducible workflow
  "rstudioapi",  # Library for reproducible workflow
  
  #EXTRACT ENVIRONMENTAL DATA AND PLOTS
  "sp", # spatial data
  "raster", #spatial data
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggcorrplot",
  
  #CORRELATION ANALYSIS
  "GGally", #correlation analysis
  "HH" #calculate VIF
)

install_load_function <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

install_load_function(requiredPackages)

#PRESENCE ABSENCE (PA) DATA FOR MODELLING

# Load presence-absence data
load("gam_modelling/outpout_modelling/PAdata.RData")

# Load environmental rasters
mylayers<-stack("gam_modelling/env/mylayers.tif")

#Extract values for every presence-absence data in raster stack
raster_ex <- raster::extract(x=mylayers, y=PAdata[,c("LON","LAT")], method="bilinear", na.rm=TRUE, df=T) 

colnames(raster_ex)[-1]<-c("BO22_damean",
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

head(raster_ex)

##Merge environment data and resence absence data 
data <- cbind(PAdata, raster_ex)
dim(data)
str(data)
head(data)
summary(data)
#errase the data NA 
data <- data %>% 
  dplyr::select (-YEAR) %>% #we remove year column because pseudoabsences miss this info
  na.omit()
save(list="data", file="gam_modelling/outpout_modelling/PA_DATA_with_env.RData")

###Explratory Plots 

tmp <- data[, c("BO22_damean",
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
                "BO2_salinitymean_ss")]
tmp <- pivot_longer(data=tmp, cols=everything()) 

#bloxplot exploratory
ggplot(data=tmp, aes(x=name, y=value)) + 
  geom_boxplot()+
  facet_wrap(~name, scales="free")

#plots
ggplot(data=tmp, aes(x=name, y=value)) + 
  geom_violin(fill="red", alpha=0.3)+
  geom_boxplot(width=0.1)+
  facet_wrap(~name, scales="free")

###Exploratory plots
tmp <- data[, c("BO22_damean",
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
                "BO2_salinitymean_ss",
                "occurrenceStatus")]

tmp <- pivot_longer(data=tmp, cols=!occurrenceStatus) 

ggplot(data=tmp, aes(x=factor(occurrenceStatus), y=value, fill=factor(occurrenceStatus), group=factor(occurrenceStatus))) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(fill="white", width=0.1)+
  facet_wrap(~name, scales="free")+
  theme(legend.position = "bottom",legend.background = element_rect(fill = "white", colour = NA))

ggplot(data=tmp, aes(x=value, fill=factor(occurrenceStatus), group=factor(occurrenceStatus))) + 
  geom_density(lwd=1, alpha=0.3)+
  facet_wrap(~name, scales="free")+
  theme(legend.position = "bottom",legend.background = element_rect(fill = "white", colour = NA))

###Correlation analysis between variables
tmp <- data[, c("BO22_damean",
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
                "BO2_salinitymean_ss")]

ggpairs(tmp)

mat <- cor(tmp, use="complete.obs") 
p.mat <- cor_pmat(tmp)

ggcorrplot(mat, type = "lower", lab=T, p.mat = p.mat)

###VIF analysis 

# Select variables for VIF calculation
library(dplyr)
v.table <- data %>% 
  dplyr::select (BO22_damean,
                 BO22_parmean,
                 BO22_ph,
                 BO_chlomean,
                 BO_dissox,
                 BO_salinity,
                 BO_sstmean,
                 BO_bathymean,
                 BO2_curvelmean_ss,
                 BO2_dissoxmean_ss, #eliminado por alta correlacion
                 BO2_ironmean_ss,
                 BO2_salinitymean_ss)

# Get VIF results
library(HH)
out.vif <- HH::vif(v.table)
sort(out.vif)

###First removal of correlated variables
v.table <- v.table %>% 
  dplyr::select (-BO2_dissoxmean_ss)
# Get new VIF results
out.vif <- vif(v.table)
sort(out.vif)

###Second removal of correlated variables
v.table <- v.table %>% 
  dplyr::select (-BO22_damean)
# Get new VIF results
out.vif <- vif(v.table)
sort(out.vif)

###Third removal of correlated variables
v.table <- v.table %>% 
  dplyr::select (-BO_dissox)
# Get new VIF results
out.vif <- vif(v.table)
sort(out.vif)

###Fourth removal of correlated variables
v.table <- v.table %>% 
  dplyr::select (-BO2_salinitymean_ss)
# Get new VIF results
out.vif <- vif(v.table)
sort(out.vif)

###Fifth removal of correlated variables
v.table <- v.table %>% 
  dplyr::select (-BO_chlomean)
# Get new VIF results
out.vif <- vif(v.table)
sort(out.vif)

##Removal of high VIF variables from the previous database
data <- data %>% dplyr::select (-BO2_dissoxmean_ss,-BO22_damean,-BO_dissox,-BO2_salinitymean_ss,
                                -BO_chlomean)
#se guarda la data para analisis de modelo
save(list="data", file="gam_modelling/outpout_modelling/PAdata_with_env7_final.RData")

