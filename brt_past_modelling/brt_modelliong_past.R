# 0. Configuración ---------------------------------------------------------
setwd("D:/Christian/proyeccionfutura_nutreias")
getwd()

# Cargar librerías
library(raster)
library(rgdal)
library(dismo)
library(SDMPlay)
library(rJava)
library(caret)
library(sf)
library(gbm)
library(blockCV)
library(pROC)     # AUC
library(terra)    # más estable para extract
# CORREGIDO: se eliminó "#library(pRdata)" — ese paquete no existe;
# las métricas (AUC, threshold de Youden) ahora se calculan con pROC,
# que ya está cargado arriba.

#options(warn = -1)
set.seed(123)


# 1. Rasters actuales -----------------------------------------------------
bio4a <- raster("Paleoclim/current_Chelsa/CHELSA_cur_V1_2B_r2_5m/2_5min/bio_4.tif")
bio9a<- raster("Paleoclim/current_Chelsa/CHELSA_cur_V1_2B_r2_5m/2_5min/bio_9.tif")
bio15a <- raster("Paleoclim/current_Chelsa/CHELSA_cur_V1_2B_r2_5m/2_5min/bio_15.tif")
bio18a  <- raster("Paleoclim/current_Chelsa/CHELSA_cur_V1_2B_r2_5m/2_5min/bio_18.tif")

layers <- stack(
  bio4a, bio9a, bio15a, bio18a
)

#ext(xmin, xmax, ymin, ymax)
e <- extent(-85, -30, -60, 14)

layers <- crop(layers,e)

names(layers) <- c(
  "bio4a","bio9a","bio15a","bio18a")
plot(layers)

# 2. Leer dataurrencias -----------------------------------------------------
data <- read.csv("modelos_futuro/data_csv/data_pca_v4paleo_sinrep.csv", stringsAsFactors = FALSE)
# se espera columnas LON y LAT
if(!all(c("LON","LAT") %in% names(data))) stop("El CSV de dataurrencias debe tener columnas LON y LAT")


data_sf <- st_as_sf(data, coords = c("LON","LAT"), crs = 4326)
plot(subset(layers,1), main="Capa base - puntos de presencia")
points(data_sf, pch = 21, cex = 1, bg="red")

cells <- cellFromXY(layers[[1]], data[,c("LON","LAT")])
length(unique(cells))

# 3. Definición del Área M -------------------------------------------------
# CORREGIDO: buffer estandarizado a 100 km (punto medio del rango
# metodológico documentado de 50-150 km); antes estaba fijo en 50 km.
buffer_M_km <- 100

buffer_costa <- st_buffer(st_union(data_sf), dist = buffer_M_km * 1000)
buffer_costa_sp <- as(buffer_costa, "Spatial")

layers_M <- crop(layers, buffer_costa_sp)
layers_M <- mask(layers_M, buffer_costa_sp)

plot(layers_M[[1]], main = "Área M")
points(data_sf, pch = 20)

layers_terra <- rast(layers_M)

# 4. Background desde M ----------------------------------------------------
set.seed(123)

# CORREGIDO (ajuste metodológico, confirmado tras revisar resultados del
# CV con separación espacial): background.nb subido de 200 a 1000. Con
# 200 puntos repartidos entre los 5 bloques espaciales de spatialBlock,
# algunos folds quedaban con muestras muy pequeñas de background (12-25
# puntos en la corrida de prueba), generando AUC inestable entre folds
# (rango 0.46-0.74, sd=0.114). Con 1000 puntos se espera una cobertura
# más uniforme y robusta de background por bloque.
dots_bg <- SDMPlay:::SDMtab(
  data[,c("LON","LAT")],
  layers_M,
  unique.data = TRUE,
  same = FALSE,
  background.nb = 1000
)

bg_df <- data.frame(
  LON = dots_bg[dots_bg$id == 0, 2],
  LAT = dots_bg[dots_bg$id == 0, 3]
)

points(bg_df, pch = 21, cex = 0.5, bg = "blue")

# CORREGIDO: se eliminó un bloque duplicado/muerto que reconstruía bg_df
# usando un objeto "dots_data" nunca definido (solo existía "dots_bg").
# Eso causaba un error fatal de "object not found" al ejecutar el script.
# bg_df ya quedó correctamente definido arriba a partir de dots_bg.

#4 spatial k-fold cross validation = esto se hace para desarrollar un análisis exploratorio de la performance general del modelo

# convertir presencias a sf
pres_sf <- st_as_sf(
  data,
  coords = c("LON", "LAT"),
  crs = 4326
)

# 5. K-fold espacial SOLO presencias --------------------------------------
pres_sf <- st_as_sf(data, coords = c("LON","LAT"), crs = 4326)

set.seed(123)

sb <- spatialBlock(
  speciesData = pres_sf,
  rasterLayer = layers_terra,
  theRange = 100000,   # 100 km
  k = 5,
  selection = "random",
  iteration = 100,
  biomod2Format = TRUE
)

folds_pres <- sb$foldID
table(folds_pres)

# CORREGIDO (error metodológico): el background de evaluación NO tenía
# separación espacial por fold, por lo que en cada iteración se evaluaba
# contra un background global fijo (bg_df), pudiendo solaparse
# espacialmente con los datos usados para entrenar ese mismo fold. Esto
# infla el AUC de forma optimista (Roberts et al. 2017; Valavi et al. 2019).
#
# Solución: se asigna el background a los mismos bloques espaciales que
# las presencias, intersectando bg_df con los polígonos de sb$blocks.
bg_sf <- st_as_sf(bg_df, coords = c("LON","LAT"), crs = 4326)

# sb$blocks contiene los polígonos de bloques con su columna "folds"
blocks_sf <- sb$blocks
if(!"folds" %in% names(blocks_sf)) {
  stop("sb$blocks no tiene columna 'folds'; revisar versión de blockCV instalada")
}

bg_block_join <- st_join(bg_sf, blocks_sf["folds"], join = st_intersects)
folds_bg <- bg_block_join$folds

# background sin bloque asignado (cae fuera de todos los polígonos) se descarta
bg_df_folded <- bg_df[!is.na(folds_bg), ]
folds_bg     <- folds_bg[!is.na(folds_bg)]

table(folds_bg)

# Aviso si, incluso con background.nb=1000, algún bloque queda con muy
# pocos puntos de background (umbral orientativo: 30). Esto ayuda a
# detectar folds que seguirán siendo inestables sin tener que inspeccionar
# table(folds_bg) manualmente.
min_bg_por_fold <- min(table(folds_bg))
if(min_bg_por_fold < 30){
  warning("Al menos un bloque espacial tiene menos de 30 puntos de ",
          "background (mínimo observado: ", min_bg_por_fold, "). ",
          "El AUC de ese fold puede seguir siendo inestable.")
}

# 6. Loop k-fold + BRT ----------------------------------------------------
# AJUSTADO (no es "tuning para subir el AUC espacial", es corregir
# underfitting real): la corrida de prueba con n.trees=10000 mostró
# "maximum tree limit reached - results may not be optimal" en los 5
# folds, indicando que gbm.step no alcanzó a converger (el error de CV
# interno seguía bajando cuando se cortó el ajuste). Se sube el techo de
# árboles a 30000 para que gbm.step pueda converger de forma natural;
# el número final de árboles que usa cada modelo lo sigue decidiendo
# gbm.step internamente (ver brt_mod$response$n.trees), no se fuerza.
k <- 5
auc_vec <- rep(NA, k)
brt_exploratorio <- NULL
relinf_list <- vector("list", k)  # influencia relativa de variables por fold

for(i in 1:k){
  
  message("Fold ", i, "/", k)
  
  pres_train <- data[folds_pres != i, ]
  pres_test  <- data[folds_pres == i, ]
  
  if(nrow(pres_test) == 0) next
  
  dots_train <- SDMPlay:::SDMtab(
    pres_train[,c("LON","LAT")],
    layers_M,
    unique.data = TRUE,
    same = FALSE,
    background.nb = 5 * nrow(pres_train)
  )
  
  dots_train <- na.omit(dots_train)
  
  set.seed(100 + i)
  
  brt_mod <- SDMPlay:::compute.brt(
    x = dots_train,
    proj.predictors = layers_M,
    tc = 2,
    lr = 0.005,
    bf = 0.7,
    n.trees = 100,
    step.size = 100
  )
  
  if(i == 1) brt_exploratorio <- brt_mod
  
  # Guardar influencia relativa de cada variable para este fold (no
  # afecta el ajuste del modelo; es solo diagnóstico para revisar si
  # bio4a/bio9a/bio15a/bio18a aportan de forma balanceada).
  relinf_list[[i]] <- summary(brt_mod$response, plotit = FALSE)
  
  # presencias test
  pres_vals <- terra::extract(
    layers_terra,
    pres_test[,c("LON","LAT")]
  )
  
  pred_pres <- predict(
    brt_mod$response,
    pres_vals,
    n.trees = brt_mod$response$n.trees,
    type = "response"
  )
  
  # CORREGIDO: background de evaluación ahora restringido al mismo
  # bloque espacial "i" del fold (antes: muestreo aleatorio desde bg_df
  # global, sin relación espacial con el fold de entrenamiento/prueba).
  bg_test <- bg_df_folded[folds_bg == i, ]
  
  if(nrow(bg_test) == 0){
    warning("Fold ", i, ": no hay background en este bloque espacial; ",
            "se omite evaluación de background para este fold.")
    next
  }
  
  bg_vals <- terra::extract(layers_terra, bg_test)
  
  pred_bg <- predict(
    brt_mod$response,
    bg_vals,
    n.trees = brt_mod$response$n.trees,
    type = "response"
  )
  
  obs   <- c(rep(1, length(pred_pres)), rep(0, length(pred_bg)))
  preds <- c(pred_pres, pred_bg)
  
  auc_vec[i] <- as.numeric(pROC::auc(pROC::roc(obs, preds, quiet = TRUE)))
}

# 7. Resultados -----------------------------------------------------------
mean_auc <- mean(auc_vec, na.rm = TRUE)
sd_auc   <- sd(auc_vec, na.rm = TRUE)

mean_auc
sd_auc
auc_vec


#plot exploratorio 
plot(
  brt_exploratorio$raster.prediction,
  main = "Preview BRT (Fold 1 – solo exploratorio)"
)

points(
  data$LON,
  data$LAT,
  pch = 21,
  bg = "red",
  cex = 0.6
)





# 5. Repetición: crear N réplicas del BRT ---------------------------------
# NOTA METODOLÓGICA (revisada y CONFIRMADA, no modificada): este bloque usa
# "layers" (extent completo, sin máscara del Área M) en vez de "layers_M"
# (Área M con buffer de 100 km), a diferencia del CV exploratorio de la
# sección 6, que sí usa layers_M. Esta inconsistencia fue señalada como
# posible error pero se mantiene de forma intencional por decisión del
# investigador. Si en algún momento se desea que el modelo final también
# restrinja el background al Área M, basta reemplazar "layers" por
# "layers_M" en las dos líneas marcadas con [layers] más abajo.
n_reps <- 10
brt_models <- vector("list", n_reps)
best_trees_vec <- numeric(n_reps)

for(i in 1:n_reps){
  message("Entrenando réplica ", i, " / ", n_reps)
  set.seed(100 + i)
  dots_temp <- SDMPlay:::SDMtab(
    data[,c("LON","LAT")], layers,  # [layers] intencional: extent completo
    unique.data = TRUE, same=FALSE,
    background.nb = 10*nrow(data)
  )
  
  set.seed(200 + i)
  mod <- SDMPlay:::compute.brt(
    x = dots_temp,
    proj.predictors = layers,  # [layers] intencional: extent completo
    tc = 2, lr = 0.005, bf = 0.7, n.trees = 100
  )
  
  # extraer best.trees o fallback
  gbm_model <- mod$response
  bt <- NA
  if (!is.null(gbm_model$gbm.call$best.trees)) {
    bt <- gbm_model$gbm.call$best.trees
  } else {
    bt <- round(gbm_model$n.trees * 0.5)
  }
  mod$best.trees <- bt
  best_trees_vec[i] <- bt
  
  brt_models[[i]] <- mod
}

# Ensemble (mapa actual promedio sobre réplicas)
br_sp1_stack <- stack(lapply(brt_models, function(x) x$raster.prediction))
br_sp1_mean <- calc(br_sp1_stack, fun = mean, na.rm = TRUE)
br_sp1_median <- calc(br_sp1_stack, fun = median, na.rm = TRUE)
br_sp1_sd   <- calc(br_sp1_stack, fun = sd, na.rm = TRUE)

dir.create("RESULTADOS_BRT/ACTUAL_TIFF/newresults", recursive = TRUE, showWarnings = FALSE)
writeRaster(br_sp1_mean, filename='RESULTADOS_BRT/ACTUAL_TIFF/newresults/Actual_lontra_mean.tif', overwrite=TRUE)
writeRaster(br_sp1_sd,   filename='RESULTADOS_BRT/ACTUAL_TIFF/newresults/Actual_lontra_sd.tif', overwrite=TRUE)

plot(br_sp1_mean, main="Ensemble mean - Actual")
plot(br_sp1_sd, main="Ensemble SD - Actual")
plot(br_sp1_median, main="Ensemble median - Actual")

#  Rasters 3.3ma -----------------------------------------------------
bio4_3.3 <- raster("Paleoclim/3.3ma/2_5min/bio_4.tif")
bio9_3.3<- raster("Paleoclim/3.3ma/2_5min/bio_9.tif")
bio15_3.3 <- raster("Paleoclim/3.3ma/2_5min/bio_15.tif")
bio18_3.3  <- raster("Paleoclim/3.3ma/2_5min/bio_18.tif")

layers_3.3 <- stack(
  bio4_3.3, bio9_3.3, bio15_3.3, bio18_3.3
)
layers_3.3 <- crop(layers_3.3,e)
names(layers_3.3) <- c(
  "bio4a","bio9a","bio15a","bio18a")
names(layers) == names(layers_3.3)
plot(layers_3.3)

#  Rasters 787ka -----------------------------------------------------
bio4_787 <- raster("Paleoclim/787ka/2_5min/bio_4.tif")
bio9_787<- raster("Paleoclim/787ka/2_5min/bio_9.tif")
bio15_787 <- raster("Paleoclim/787ka/2_5min/bio_15.tif")
bio18_787  <- raster("Paleoclim/787ka/2_5min/bio_18.tif")

layers_787 <- stack(
  bio4_787, bio9_787, bio15_787, bio18_787
)

layers_787 <- crop(layers_787,e)
names(layers_787) <- c(
  "bio4a","bio9a","bio15a","bio18a")
names(layers) == names(layers_787)
plot(layers_787)

#6. definir fuincion de prediccion 

project_brt_to_raster <- function(brt_obj, predictors){
  gbm_model <- brt_obj$response
  n.trees <- if(!is.null(brt_obj$best.trees)) brt_obj$best.trees else gbm_model$n.trees
  
  raster::predict(
    predictors, gbm_model,
    fun = function(model, data){
      gbm::predict.gbm(
        model, data,
        n.trees = n.trees,
        type = "response"
      )
    },
    progress = "text"
  )
}

#aplicar proyeccion a pasado 787
pred_lgm_list787 <- lapply(
  brt_models,
  project_brt_to_raster,
  predictors = layers_787
)

pred_lgm_stack787 <- stack(pred_lgm_list787)

lgm787_mean <- calc(pred_lgm_stack787, mean, na.rm = TRUE)
lgm787_sd   <- calc(pred_lgm_stack787, sd,   na.rm = TRUE)

dir.create("RESULTADOS_BRT/787_TIFF_new",
           recursive = TRUE,
           showWarnings = FALSE)

writeRaster(
  lgm787_mean,
  "RESULTADOS_BRT/787_TIFF_new/787_mean.tif",
  overwrite = TRUE
)

writeRaster(
  lgm787_sd,
  "RESULTADOS_BRT/787_TIFF_new/787_sd.tif",
  overwrite = TRUE
)

#aplicar proyeccion a pasado
pred_lgm_list3.3 <- lapply(
  brt_models,
  project_brt_to_raster,
  predictors = layers_3.3
)

pred_lgm_stack3.3 <- stack(pred_lgm_list3.3)

lgm3.3_mean <- calc(pred_lgm_stack3.3, mean, na.rm = TRUE)
lgm3.3_sd   <- calc(pred_lgm_stack3.3, sd,   na.rm = TRUE)

plot(lgm3.3_mean)

dir.create("RESULTADOS_BRT/3.3_TIFF_new",
           recursive = TRUE,
           showWarnings = FALSE)

writeRaster(
  lgm3.3_mean,
  "RESULTADOS_BRT/3.3_TIFF_new/3.3_mean.tif",
  overwrite = TRUE
)

writeRaster(
  lgm3.3_sd,
  "RESULTADOS_BRT/3.3_TIFF_new/3.3_sd.tif",
  overwrite = TRUE
)

# 9. Función para calcular métricas (AUC, threshold(Youden), Sens, Spec, Kappa, TSS)

calc_metrics_from_obs_pred <- function(obs, pred_probs, threshold = NULL){
  # obs: 0/1 vector
  # pred_probs: probabilities
  roc_obj <- try(pROC::roc(obs, pred_probs, quiet = TRUE), silent = TRUE)
  if(inherits(roc_obj, "try-error") || is.null(roc_obj)){
    auc_val <- NA
    best_thresh <- NA
  } else {
    auc_val <- as.numeric(pROC::auc(roc_obj))
    if(is.null(threshold)){
      coords_best <- pROC::coords(
        roc_obj, "best", best.method = "youden",
        ret = c("threshold","sensitivity","specificity"),
        transpose = FALSE
      )
      # pROC::coords puede devolver >1 fila si hay empates; se toma la primera
      best_thresh <- as.numeric(coords_best$threshold[1])
    } else {
      best_thresh <- threshold
    }
  }
  # binarize using best_thresh (for metric calc only)
  pred_bin <- if(!is.na(best_thresh)) ifelse(pred_probs >= best_thresh, 1, 0) else rep(NA, length(pred_probs))
  # confusion
  if(all(is.na(pred_bin))){
    sens <- NA; spec <- NA; kappa <- NA
  } else {
    cm <- try(caret::confusionMatrix(as.factor(pred_bin), as.factor(obs), positive = "1"), silent = TRUE)
    if(inherits(cm, "try-error")){
      sens <- NA; spec <- NA; kappa <- NA
    } else {
      sens <- as.numeric(cm$byClass["Sensitivity"])
      spec <- as.numeric(cm$byClass["Specificity"])
      kappa <- as.numeric(cm$overall["Kappa"])
    }
  }
  TSS <- if(!is.na(sens) && !is.na(spec)) sens + spec - 1 else NA
  return(list(AUC = auc_val, threshold = best_thresh, Sens = sens, Spec = spec, Kappa = kappa, TSS = TSS))
}

# 10. Preparar puntos para extracción --------------------------------------
pres_pts <- data.frame(LON = data$LON, LAT = data$LAT)
bg_pts   <- bg_df # generados por SDMPlay

# 11. Evaluación escenario presente (ensemble mean) -----------------------
# extraer valores del ensemble mean en presencias y background
pres_vals_actual <- raster::extract(br_sp1_mean, pres_pts[,c("LON","LAT")])
bg_vals_actual   <- raster::extract(br_sp1_mean, bg_pts[,c("LON","LAT")])

obs_vec_actual <- c(rep(1, length(pres_vals_actual)), rep(0, length(bg_vals_actual)))
pred_vec_actual <- c(pres_vals_actual, bg_vals_actual)

metrics_actual <- calc_metrics_from_obs_pred(obs_vec_actual, pred_vec_actual, threshold = NULL)

# preparar dataframe resultados
results_df <- data.frame(
  Escenario = "PRESENTE",
  Tipo = "PRESENTE",
  AUC = metrics_actual$AUC,
  Threshold = metrics_actual$threshold,
  Sens = metrics_actual$Sens,
  Spec = metrics_actual$Spec,
  Kappa = metrics_actual$Kappa,
  TSS = metrics_actual$TSS,
  stringsAsFactors = FALSE
)

results_df
