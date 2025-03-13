# ---
# title: "Non-linear Relationship Exploration (NTEMS Data) - Negative Binomial GLMM with Offset"
# author: "Your Name"
# created: "2024-05-09"
# description: "Exploring non-linear relationships with negative binomial GLMMs and AIC comparison, including offset"
# ---

# Load packages ----
library(lme4)
library(MuMIn)
library(fitdistrplus)
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(segmented) 

# Load data
NTEMS <- read.csv("Output/Tabular Data/NTEMS_w_clustering.csv")

# Convert site clustering random effects to factors
NTEMS$spatial_group_1000 <- as.factor(NTEMS$spatial_group_1000)
NTEMS$spatial_group_2000 <- as.factor(NTEMS$spatial_group_2000)
NTEMS$spatial_group_3000 <- as.factor(NTEMS$spatial_group_3000)
NTEMS$spatial_group_4000 <- as.factor(NTEMS$spatial_group_4000)
NTEMS$spatial_group_5000 <- as.factor(NTEMS$spatial_group_5000)
NTEMS$spatial_group_6000 <- as.factor(NTEMS$spatial_group_6000)
NTEMS$spatial_group_7000 <- as.factor(NTEMS$spatial_group_7000)
NTEMS$spatial_group_8000 <- as.factor(NTEMS$spatial_group_8000)
NTEMS$spatial_group_9000 <- as.factor(NTEMS$spatial_group_9000)
NTEMS$spatial_group_10000 <- as.factor(NTEMS$spatial_group_10000)
NTEMS$spatial_group_11000 <- as.factor(NTEMS$spatial_group_11000)
NTEMS$spatial_group_12000 <- as.factor(NTEMS$spatial_group_12000)
NTEMS$spatial_group_13000 <- as.factor(NTEMS$spatial_group_13000)
NTEMS$spatial_group_14000 <- as.factor(NTEMS$spatial_group_14000)
NTEMS$spatial_group_15000 <- as.factor(NTEMS$spatial_group_15000)


# Random effect models
model_lin1 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_1000), data = NTEMS)
model_lin2 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_2000), data = NTEMS)
model_lin3 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_3000), data = NTEMS)
model_lin4 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_4000), data = NTEMS)
model_lin5 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)
model_lin6 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_6000), data = NTEMS)
model_lin7 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_7000), data = NTEMS)
model_lin8 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_8000), data = NTEMS)
model_lin9 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_9000), data = NTEMS)
model_lin10 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)
model_lin11 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_11000), data = NTEMS)
model_lin12 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_12000), data = NTEMS)
model_lin13 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_13000), data = NTEMS)
model_lin14 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_14000), data = NTEMS)
model_lin15 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)




###### Compare random effects using AIC

# Create a list of your models
model_list <- list(
  model_lin1, model_lin2, model_lin3, model_lin4, model_lin5,
  model_lin6, model_lin7, model_lin8, model_lin9, model_lin10, model_lin11,
  model_lin12, model_lin13, model_lin14, model_lin15
)

# Extract AIC values and create a data frame
aic_values <- data.frame(
  model_name = paste0("model_lin", 1:15),
  aic = sapply(model_list, AIC)
)

# Find the model with the lowest AIC
best_model <- aic_values[which.min(aic_values$aic), ]

# Print the results
print(aic_values)
print(paste("The model with the lowest AIC is:", best_model$model_name))
print(paste("AIC:", best_model$aic))

# If you want to know the spatial group associated with the best model
best_model_number <- as.numeric(gsub("model_lin", "", best_model$model_name))

# Adjust the spatial group calculation based on your model definitions
spatial_group_names <- paste0("spatial_group_", seq(1000, 15000, by = 1000))
best_spatial_group <- spatial_group_names[best_model_number]

print(paste("Spatial group associated with lowest AIC is:", best_spatial_group))













################# BTNW

### Prop_con_1

# Linear model
model_lin <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_1), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(BTNW ~ poly(prop_con_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(BTNW ~ poly(prop_con_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Fourth model
model_fourth <- glmer.nb(BTNW ~ poly(prop_con_1, 4, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(BTNW ~ log(prop_con_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(BTNW ~ exp(prop_con_1) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(BTNW ~ I(1 / (1 + exp(-prop_con_1))) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(BTNW ~ I(prop_con_1 / (1 + prop_con_1)) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glmer.nb(BTNW ~ prop_con_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise), AIC(model_fourth))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

######## Cubic is best






### Clumpy_1


# Linear model
model_lin <- glmer.nb(BTNW ~ clumpy_1 + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(BTNW ~ poly(clumpy_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(BTNW ~ poly(clumpy_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Fourth model
model_fourth <- glmer.nb(BTNW ~ poly(clumpy_1, 4, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(BTNW ~ log(clumpy_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(BTNW ~ exp(clumpy_1) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(BTNW ~ I(1 / (1 + exp(-clumpy_1))) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(BTNW ~ I(clumpy_1 / (1 + clumpy_1)) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$clumpy_1)
NTEMS$above_thresh <- as.numeric(NTEMS$clumpy_1 > threshold)
model_piecewise <- glmer.nb(BTNW ~ clumpy_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


Model      AIC
5      Logarithmic 4642.886
4           Fourth 4895.820
6      Exponential 4897.215
1           Linear 4897.276
7         Logistic 4897.328
8 Michaelis-Menten 4897.628
2        Quadratic 4899.251
9        Piecewise 4899.435
3            Cubic 4900.496

#  is best




### Age 

# Linear model
model_lin <- glmer.nb(BTNW ~ age_mn_1 + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(BTNW ~ poly(age_mn_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(BTNW ~ poly(age_mn_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Fourth model
model_fourth <- glmer.nb(BTNW ~ poly(age_mn_1, 4, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(BTNW ~ log(age_mn_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(BTNW ~ exp(age_mn_1) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(BTNW ~ I(1 / (1 + exp(-age_mn_1))) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(BTNW ~ I(clumpy_1 / (1 + age_mn_1)) + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$age_mn_1)
NTEMS$above_thresh <- as.numeric(NTEMS$age_mn_1 > threshold)
model_piecewise <- glmer.nb(BTNW ~ age_mn_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group_5000), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

########  is best















#################### TEWA



### Prop_con_1

# Linear model
model_lin <- glmer.nb(TEWA ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(TEWA ~ poly(prop_con_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(TEWA ~ poly(prop_con_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Fourth model
model_fourth <- glmer.nb(TEWA ~ poly(clumpy_1, 4, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(TEWA ~ log(prop_con_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(TEWA ~ exp(prop_con_1) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(TEWA ~ I(1 / (1 + exp(-prop_con_1))) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(TEWA ~ I(prop_con_1 / (1 + prop_con_1)) + offset(log(survey_effort)) + (1|spatial_group_15000), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glmer.nb(TEWA ~ prop_con_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


### Quadratic is best










### Clumpy_1

# Linear model
model_lin <- glmer.nb(TEWA ~ clumpy_1 + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(TEWA ~ poly(clumpy_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(TEWA ~ poly(clumpy_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(TEWA ~ log(clumpy_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(TEWA ~ exp(clumpy_1) + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(TEWA ~ I(1 / (1 + exp(-clumpy_1))) + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(TEWA ~ I(clumpy_1 / (1 + clumpy_1)) + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$clumpy_1)
NTEMS$above_thresh <- as.numeric(NTEMS$clumpy_1 > threshold)
model_piecewise <- glmer.nb(TEWA ~ clumpy_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group_10000), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

# Log is best




### Age 

# Linear model
model_lin <- glmer.nb(TEWA ~ age_mn_1 + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(TEWA ~ poly(age_mn_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(TEWA ~ poly(age_mn_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(TEWA ~ log(age_mn_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(TEWA ~ exp(age_mn_1) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(TEWA ~ I(1 / (1 + exp(-age_mn_1))) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(TEWA ~ I(clumpy_1 / (1 + age_mn_1)) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$age_mn_1)
NTEMS$above_thresh <- as.numeric(NTEMS$age_mn_1 > threshold)
model_piecewise <- glmer.nb(TEWA ~ age_mn_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ]) # Sort models by AIC (lower is better)

########  is best












################# BBWA

### Prop_con_1

# Linear model
model_lin <- glmer.nb(BBWA ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group_1), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(BBWA ~ poly(prop_con_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(BBWA ~ poly(prop_con_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(BBWA ~ log(prop_con_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(BBWA ~ exp(prop_con_1) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(BBWA ~ I(1 / (1 + exp(-prop_con_1))) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(BBWA ~ I(prop_con_1 / (1 + prop_con_1)) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glmer.nb(BBWA ~ prop_con_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

######## Cubic is best






### Clumpy_1

# Linear model
model_lin <- glmer.nb(BBWA ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(BBWA ~ poly(prop_con_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(BBWA ~ poly(prop_con_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(BBWA ~ log(prop_con_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(BBWA ~ exp(prop_con_1) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(BBWA ~ I(1 / (1 + exp(-prop_con_1))) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(BBWA ~ I(prop_con_1 / (1 + prop_con_1)) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glmer.nb(BBWA ~ prop_con_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

# Cubic is best




### Age 

# Linear model
model_lin <- glmer.nb(BBWA ~ age_mn_1 + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Quadratic model
model_quad <- glmer.nb(BBWA ~ poly(age_mn_1, 2, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Cubic model
model_cubic <- glmer.nb(BBWA ~ poly(age_mn_1, 3, raw = TRUE) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Logarithmic model (adding a small constant to avoid log(0))
model_log <- glmer.nb(BBWA ~ log(age_mn_1 + 0.001) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Exponential model
model_exp <- glmer.nb(BBWA ~ exp(age_mn_1) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glmer.nb(BBWA ~ I(1 / (1 + exp(-age_mn_1))) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Michaelis-Menten (Saturating) model
model_mm <- glmer.nb(BBWA ~ I(clumpy_1 / (1 + age_mn_1)) + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$age_mn_1)
NTEMS$above_thresh <- as.numeric(NTEMS$age_mn_1 > threshold)
model_piecewise <- glmer.nb(BBWA ~ age_mn_1 * above_thresh + offset(log(survey_effort)) + (1|spatial_group), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Logarithmic", "Exponential", "Logistic", "Michaelis-Menten", "Piecewise"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_log),
          AIC(model_exp), AIC(model_logistic), AIC(model_mm), AIC(model_piecewise))
)

print(aic_values[order(aic_values$AIC), ]) # Sort models by AIC (lower is better)

########  is best
















































# --- MODEL SPECIFICATION (EXAMPLE - CHANGE THESE FOR EACH MODEL) ---
response <- "BTNW"
predictor <- "prop_con_"
scale <- "1"
predictor_var <- paste0(predictor, scale)
# --- END MODEL SPECIFICATION ---

if (predictor_var %in% names(NTEMS)) {
  # Linear model
  nl1 <- glmer.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("Linear Model AIC:")
  print(AIC(nl1))
  
  # Squared model
  nl2 <- glmer.nb(BTNW ~ I(prop_con_1^2) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("Squared Model AIC:")
  print(AIC(nl2))
  
  # Cubic model
  nl3 <- glmer.nb(BTNW ~ I(prop_con_1^3) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("Cubic Model AIC:")
  print(AIC(nl3))
  
  # ^4 model
  nl4 <- glmer.nb(BTNW ~ I(prop_con_1^4) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("^4 Model AIC:")
  print(AIC(nl4))
  
  # ^5 model
  nl5 <- glmer.nb(BTNW ~ I(prop_con_1^5) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("^5 Model AIC:")
  print(AIC(nl5))
  
  # Log model
  nl6 <- glmer.nb(BTNW ~ log(prop_con_1 + 1) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("Log Model AIC:")
  print(AIC(nl6))
  
  # Square root model
  nl7 <- glmer.nb(BTNW ~ sqrt(prop_con_1) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("Square Root Model AIC:")
  print(AIC(nl7))
  
  # 1/n model
  nl8 <- glmer.nb(BTNW ~ I(1/(prop_con_1 + 0.0001)) + offset(log(survey_effort)) + (1|spatial_group),
                  data = NTEMS, na.action = na.exclude)
  print("1/n Model AIC:")
  print(AIC(nl8))
  
  # Piecewise regression (segmented)
  nl9_formula <- BTNW ~ prop_con_1 + offset(log(survey_effort)) + (1|spatial_group)
  nl9_data <- NTEMS[!is.na(NTEMS$BTNW) & !is.na(NTEMS$prop_con_1) & !is.na(NTEMS$survey_effort), ]
  nl9_glmer <- glmer.nb(nl9_formula, data = nl9_data)
  tryCatch({
    nl9_segmented <- segmented(nl9_glmer, seg.Z = ~prop_con_1, control = seg.control(display = FALSE))
    print("Piecewise Model AIC:")
    print(AIC(nl9_segmented))
  }, error = function(e) {
    print("Piecewise Model Failed to converge.")
  })
  
  # Plot if squared model is best
  if (AIC(nl2) == min(AIC(nl1), AIC(nl2), AIC(nl3), AIC(nl4), AIC(nl5), AIC(nl6), AIC(nl7), AIC(nl8))) {
    sc_1 <- attr(NTEMS$prop_con_1, 'scaled:scale')
    ce_1 <- attr(NTEMS$prop_con_1, 'scaled:center')
    
    dd2 <- ggpredict(nl2, terms = "prop_con_1[all]") %>%
      as_tibble() %>%
      mutate(x = x * sc_1 + ce_1)
    
    peak_x <- dd2$x[which.max(dd2$predicted)]
    
    plot <- ggplot(dd2) +
      aes(x, predicted, ymin = conf.low, ymax = conf.high) +
      geom_ribbon(alpha = 0.3) +
      geom_line() +
      labs(x = "prop_con_1", y = "BTNW") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white"),
            axis.line = element_line(colour = "black")) +
      geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")
    
    ggsave(paste0("3_output/figures/", response, "_", predictor_var, ".png"), plot,
           width = 6, height = 4, dpi = 300)
  }
} else {
  print(paste("Predictor variable", predictor_var, "not found in data."))
}








