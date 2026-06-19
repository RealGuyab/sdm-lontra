#Load Libraries
install.packages("Rtools")
install.packages("raster")
install.packages("rgdal")
install.packages("SDMPlay")
install.packages("rJava")
install.packages("caret")
install.packages("dismo")
install.packages("gbm")
install.packages("pROC")  # para AUC y ROC

# si es necesario instalar vcd para kappa:
install.packages("vcd")

library(raster)
library(rgdal)
library(dismo)
library(SDMPlay)
#library(rJava)
library(caret)
library(pROC)
library(vcd)

#3# leer los raster configurados
bathy   <- raster("brt_modelling/variables/recorte_bo_bathymean_lonlat_tif_tif.asc")
sstmean <- raster("brt_modelling/variables/recorte_bo_sstmean_lonlat_tif_tif.asc")
parmean <- raster("brt_modelling/variables/recorte_present_surface_par_mean_bov2_2_tif_tif.asc")

layers <- stack(bathy, sstmean, parmean)
names(layers) <- c("bathy",
                   "sstmean",
                   "parmean")
summary(layers)
par(mfrow=c(1,1))
plot(layers)

#4# leer las ocurrencias (oc)
oc <- read.csv("occurence_data/occu_maxent.csv")
library(sf)
oc_sf <- st_as_sf(oc, coords = c("LON", "LAT"), crs = 4326)
st_crs(oc_sf)
plot(subset(layers, 1))
points(oc_sf, pch = 21, cex = 1, bg="red")

#5# remover variables innecesarias
rm(bathy, sstmean, parmean)

# Definir semillas para reproducibilidad
global_seeds <- list(
  simple_model = 100,
  replicas = c(101, 102, 103, 104, 105, 106, 107, 108, 109, 110)
)

#6# background sampling reproducible
names(oc)
set.seed(global_seeds$simple_model)
dots_data <- SDMPlay:::SDMtab(oc[,c(2,3)], layers, unique.data = TRUE, same = FALSE, background.nb = 200)
background <- dots_data[dots_data$id == 0, 2:3]

#7# BRT simple (single run) reproducible
set.seed(global_seeds$simple_model)
dots_data_sp1 <- SDMPlay:::SDMtab(oc[,c("LON","LAT")], layers, unique.data = TRUE, same = FALSE,
                                  background.nb = 2 * nrow(oc))
brt_sp1_1 <- SDMPlay:::compute.brt(x = dots_data_sp1, proj.predictors = layers,
                                   tc = 2, lr = 0.005, bf = 0.75, n.trees = 100)

#8# Evaluacion y AUC para ese modelo simple
evaluation_sp1 <- SDMPlay:::SDMeval(brt_sp1_1)
# obtener ROC/AUC
preds1 <- predict(brt_sp1_1$response, dots_data_sp1[,c("bathy","sstmean","parmean")], 
                  n.trees = brt_sp1_1$response$gbm.call$best.trees, type = "response")
obs1 <- ifelse(dots_data_sp1$id == 1, 1, 0)
roc1 <- roc(obs1, preds1)
cat("AUC modelo simple:", auc(roc1), "\n")

#9# Configurar 10 réplicas de BRT con semillas para reproducibilidad
# Creamos plantillas de datos y modelos
datos_list <- vector("list", 10)
modelos_list <- vector("list", 10)
for(i in 1:10) {
  set.seed(global_seeds$replicas[i])
  # generar nuevo sample reproducible
  datos_list[[i]] <- SDMPlay:::SDMtab(oc[,c("LON","LAT")], layers, unique.data = TRUE, same = FALSE,
                                      background.nb = 2 * nrow(oc))
  set.seed(global_seeds$replicas[i])
  modelos_list[[i]] <- SDMPlay:::compute.brt(x = datos_list[[i]], proj.predictors = layers,
                                             tc = 2, lr = 0.005, bf = 0.75, n.trees = 100)
}

# Funcion para extraer metrics incluyendo Kappa
tiene_metricas <- function(brt_obj, data_sample, umbral) {
  preds <- predict(brt_obj$response, data_sample[,c("bathy","sstmean","parmean")], 
                   n.trees = brt_obj$response$n.trees, type = "response")
  obs <- factor(ifelse(data_sample$id == 1, "Presencia", "Ausencia"))
  bin <- factor(ifelse(preds > umbral, "Presencia", "Ausencia"))
  cm <- confusionMatrix(bin, obs)
  return(c(
    Accuracy    = cm$overall['Accuracy'],
    Kappa       = cm$overall['Kappa'],
    Sensitivity = cm$byClass['Sensitivity'],
    Specificity = cm$byClass['Specificity'],
    Precision   = cm$byClass['Precision'],
    F1_Score    = cm$byClass['F1']
  ))
}

# Umbral óptimo del primer modelo
umbral_opt <- evaluation_sp1$MaxSSS

# Data frame para resultados
df_results <- data.frame(
  Modelo = paste0("BRT_", 1:10),
  Accuracy = NA, Kappa = NA, Sensitivity = NA,
  Specificity = NA, Precision = NA, F1_Score = NA,
  stringsAsFactors = FALSE
)

# Llenar resultados
for(i in 1:10) {
  met <- tiene_metricas(modelos_list[[i]], datos_list[[i]], umbral_opt)
  df_results[i, 2:7] <- met
}

# Mostrar tabla de resultados
print(df_results)

# Calculo de Kappa promedio
kappa_promedio <- mean(df_results$Kappa, na.rm = TRUE)
cat(kappa_promedio)

#10# Raster promedio de las 10 predicciones
to_stack <- stack(lapply(modelos_list, function(x) x$raster.prediction))
br_mean <- calc(to_stack, fun = mean)
plot(br_mean, main = "Mean BRT (10 réplicas)")
writeRaster(br_mean, filename = "brt_modelling/results/BRT_model_promedio-3.asc", format = "ascii", overwrite = TRUE)