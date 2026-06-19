#instalando librerias necesarias
requiredPackages <- c("here", "rstudioapi",
                      "stringr", "RColorBrewer", "ggplot2",
                      "dplyr", "tidyverse", "R.utils", "ggpubr",
                      "hrbrthemes", "fields", "maps", "raster",
                      "scam", "plotmo", "pkgbuild", "dismo",
                      "SDMTools","caret")
install_load_function <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[,
                                                 "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

install_load_function(requiredPackages)

#installing SDMTools manually
# find_rtools()
# install.packages('remotes')
# remotes::install_version('SDMTools',
# version = '1.1-221.2')
library(SDMTools)

#Preliminar metrics
library(caret)

# Load Data
load(here::here("gam_modelling/outpout_modelling/selected_model.Rdata"))
data <- model$model
scgam.pred <- predict(model, newdata = data, type = "response")
data$scgam.pred <- as.vector(scgam.pred)

obs <- data$occurrenceStatus
predSCGAM_P <- data$scgam.pred

# Optimice Threshold
myoptim <- optim.thresh(obs, predSCGAM_P)
myThreshold <- median(as.numeric(myoptim[["max.sensitivity+specificity"]]))
myoptim

pred_bin <- ifelse(predSCGAM_P >= myThreshold, 1, 0)

#convert into factors
obs_f <- factor(obs, levels = c(1, 0))        
pred_f <- factor(pred_bin, levels = c(1, 0)) 

#confusion Matrix and metrics
conf <- confusionMatrix(pred_f, obs_f)
conf

kappa_value <- conf$overall["Kappa"]
accuracy_value <- conf$overall["Accuracy"]
sensitivity <- conf$byClass["Sensitivity"]
specificity <- conf$byClass["Specificity"]
precision   <- conf$byClass["Pos Pred Value"]
npv         <- conf$byClass["Neg Pred Value"]
f1_score    <- 2 * ((precision * sensitivity) / (precision + sensitivity))

results <- list(
  Threshold   = myThreshold,
  Kappa       = kappa_value,
  Accuracy    = accuracy_value,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision   = precision,
  NPV         = npv,
  F1_Score    = f1_score
)

print(results)

###Cross validation K-fold

# Number of groups
set.seed(123)
k <- 5

# Generate groups
groups <- kfold(data, k, by = data$occurrencestatus)

# Initialise the confusion matrix and
# the accuracy table:
myCM <- NULL
myACC <- NULL

# Get the formula of the selected model
formula <- summary(model)[["formula"]]

# Get the smoothing parameters of the
# selected model
sp <- model$sp

# Loop for each group k
for (j in 1:k) {
  # Preparation of Training Sites
  p_Training <- data[groups != j, ]
  
  # Model fit
  selected_model.sp.j <- scam(formula,
                              family = binomial(link = "logit"),
                              data = p_Training, sp = c(sp))
  
  # Predict Model
  p_validacion <- data[groups == j, ]
  
  model.sp.j.pred <- predict(selected_model.sp.j,
                                      newdata = p_validacion, type = "response")
  p_validacion$Pred <- model.sp.j.pred
  
  # Confussion matrix and accuracy
  # table for fold j
  obs <- p_validacion$occurrenceStatus
  predSCGAM <- p_validacion$Pred
  myCM <- rbind(myCM, as.numeric(confusion.matrix(obs,
                                                  predSCGAM, threshold = myThreshold)))
  myACC <- rbind(myACC, accuracy(obs, predSCGAM,
                                 threshold = myThreshold))
}

# Mean values across k-folds
validation_summary <- cbind(Threshold = myThreshold,
                            mean_AUC = mean(myACC$AUC), mean_Omision = mean(myACC$omission.rate),
                            mean_sensitivity = mean(myACC$sensitivity),
                            mean_specificity = mean(myACC$specificity),
                            mean_Prop.Corr = mean(myACC$prop.correct))

validation_summary

save(validation_summary, file = here::here("Nutrias/modelo/validation_summary.RData"))


