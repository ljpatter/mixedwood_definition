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
library(DHARMa)
library(MASS) 
library(dplyr) 

# Load data
NTEMS <- read.csv("Output/Tabular Data/NTEMS_w_clustering.csv")
LULC  <- read.csv("Output/Tabular Data/LULC_w_clustering.csv")

# Convert site clustering random effects to factors
NTEMS[paste0("cluster_k_", 3:10)] <- lapply(NTEMS[paste0("cluster_k_", 3:10)], as.factor)



########### Compare clustering approaches


# Create a list of cluster variables
cluster_vars <- c("cluster_k_3", "cluster_k_4", "cluster_k_5", "cluster_k_6", "cluster_k_7", "cluster_k_8", "cluster_k_9", "cluster_k_10")

# Create a list to store model results
model_results <- list()

# Loop through each cluster variable and fit a GLM
for (cluster_var in cluster_vars) {
  # Create a formula for the GLM
  formula <- as.formula(paste("BTNW ~ prop_con_1 + (1 |", cluster_var, ")"))
  
  # Fit the GLM with Poisson distribution
  model <- glmer(formula, data = NTEMS, family = poisson(link = "log"))
  
  # Store the model results
  model_results[[cluster_var]] <- model
}

# Calculate AIC for each model
aic_values <- sapply(model_results, AIC)

# Print AIC values
print(aic_values)

# Find the model with the lowest AIC
best_model <- model_results[[which.min(aic_values)]]
best_cluster_var <- names(aic_values)[which.min(aic_values)]

# Print the best model and cluster variable
cat("Best model (lowest AIC):", best_cluster_var, "\n")
print(summary(best_model))





#### Assess LULC vs NTEMS for clumpy and prop_con across all spatial scales

# Assuming NTEMS and LULC data frames are loaded and glm.nb is available

# Run models
model1 <- glm.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)), data = NTEMS)
model2 <- glm.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)), data = LULC)
model3 <- glm.nb(BTNW ~ clumpy_1 + offset(log(survey_effort)), data = NTEMS)
model4 <- glm.nb(BTNW ~ clumpy_1 + offset(log(survey_effort)), data = LULC)
model5 <- glm.nb(BTNW ~ prop_con_2 + offset(log(survey_effort)), data = NTEMS)
model6 <- glm.nb(BTNW ~ prop_con_2 + offset(log(survey_effort)), data = LULC)
model7 <- glm.nb(BTNW ~ prop_con_3 + offset(log(survey_effort)), data = NTEMS)
model8 <- glm.nb(BTNW ~ prop_con_3 + offset(log(survey_effort)), data = LULC)
model9 <- glm.nb(BTNW ~ clumpy_2 + offset(log(survey_effort)), data = NTEMS)
model10 <- glm.nb(BTNW ~ clumpy_2 + offset(log(survey_effort)), data = LULC)
model11 <- glm.nb(BTNW ~ clumpy_3 + offset(log(survey_effort)), data = NTEMS)
model12 <- glm.nb(BTNW ~ clumpy_3 + offset(log(survey_effort)), data = LULC)

# Store AIC values
aic_values <- c(
  `NTEMS, prop_con_1` = AIC(model1),
  `LULC, prop_con_1` = AIC(model2),
  `NTEMS, clumpy_1` = AIC(model3),
  `LULC, clumpy_1` = AIC(model4),
  `NTEMS, prop_con_2` = AIC(model5),
  `LULC, prop_con_2` = AIC(model6),
  `NTEMS, prop_con_3` = AIC(model7),
  `LULC, prop_con_3` = AIC(model8),
  `NTEMS, clumpy_2` = AIC(model9),
  `LULC, clumpy_2` = AIC(model10),
  `NTEMS, clumpy_3` = AIC(model11),
  `LULC, clumpy_3` = AIC(model12)
)

# Print the AIC values with names
print(aic_values)

###### NTEMS provides a better fit for prop_con_1 and clumpy_1 across all spatial scales.










################# BTNW
### Prop_con_1

# Linear model
model_lin <- glm.nb(BTNW ~ prop_con_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(BTNW ~ prop_con_1 + I(prop_con_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes linear and quadratic terms)
model_cubic <- glm.nb(BTNW ~ prop_con_1 + I(prop_con_1^2) + I(prop_con_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(BTNW ~ prop_con_1 + I(prop_con_1^2) + I(prop_con_1^3) + I(prop_con_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(BTNW ~ exp(prop_con_1) + offset(log(survey_effort)), data = NTEMS)

# Sigmoidal (Logistic) model
model_logistic <- glm.nb(BTNW ~ I(1 / (1 + exp(-prop_con_1))) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glm.nb(BTNW ~ prop_con_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Square Root model
model_sqrt <- glm.nb(BTNW ~ sqrt(prop_con_1) + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(BTNW ~ I(prop_con_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Square Root", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_sqrt), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

Model      AIC
3       Cubic 2927.879
5 Exponential 2928.662
4      Fourth 2929.712
1      Linear 2929.717
2   Quadratic 2931.003
6   Piecewise 2932.903
7 Square Root 2933.862
8  Cubic Root 2937.320

# Cubic is best




### Clumpy_1

# Linear model
model_lin <- glm.nb(BTNW ~ clumpy_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(BTNW ~ clumpy_1 + I(clumpy_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes lower-order terms)
model_cubic <- glm.nb(BTNW ~ clumpy_1 + I(clumpy_1^2) + I(clumpy_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(BTNW ~ clumpy_1 + I(clumpy_1^2) + I(clumpy_1^3) + I(clumpy_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(BTNW ~ exp(clumpy_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$clumpy_1)
NTEMS$above_thresh <- as.numeric(NTEMS$clumpy_1 > threshold)
model_piecewise <- glm.nb(BTNW ~ clumpy_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(BTNW ~ I(clumpy_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


Model      AIC
7  Cubic Root 2870.265
1      Linear 3000.184
5 Exponential 3001.234
2   Quadratic 3002.040
3       Cubic 3003.037
6   Piecewise 3003.849
4      Fourth 3004.123


# Cubic root is best







### Age

# Linear model
model_lin <- glm.nb(BTNW ~ age_mn_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(BTNW ~ age_mn_1 + I(age_mn_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes all lower-order terms)
model_cubic <- glm.nb(BTNW ~ age_mn_1 + I(age_mn_1^2) + I(age_mn_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(BTNW ~ age_mn_1 + I(age_mn_1^2) + I(age_mn_1^3) + I(age_mn_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(BTNW ~ exp(age_mn_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$age_mn_1)
NTEMS$above_thresh <- as.numeric(NTEMS$age_mn_1 > threshold)
model_piecewise <- glm.nb(BTNW ~ age_mn_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Square Root model
model_sqrt <- glm.nb(BTNW ~ sqrt(age_mn_1) + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(BTNW ~ I(age_mn_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Square Root", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_sqrt), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


Model      AIC
3       Cubic 2848.705
4      Fourth 2850.566
2   Quadratic 2854.883
6   Piecewise 2860.807
8  Cubic Root 3000.305
5 Exponential 3002.376
7 Square Root 3002.999
1      Linear 3009.431

# Cubic is best











#################### TEWA


### Prop_con_1

# Linear model
model_lin <- glm.nb(TEWA ~ prop_con_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(TEWA ~ prop_con_1 + I(prop_con_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes lower-order terms)
model_cubic <- glm.nb(TEWA ~ prop_con_1 + I(prop_con_1^2) + I(prop_con_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(TEWA ~ prop_con_1 + I(prop_con_1^2) + I(prop_con_1^3) + I(prop_con_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(TEWA ~ exp(prop_con_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glm.nb(TEWA ~ prop_con_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Square Root model
model_sqrt <- glm.nb(TEWA ~ sqrt(prop_con_1) + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(TEWA ~ I(prop_con_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Square Root", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_sqrt), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

Model      AIC
1      Linear 17384.58
4      Fourth 17385.99
5 Exponential 17386.21
2   Quadratic 17386.54
6   Piecewise 17387.95
3       Cubic 17388.46
7 Square Root 17388.47
8  Cubic Root 17392.43

# Linear is best




### Clumpy_1

# Linear model
model_lin <- glm.nb(TEWA ~ clumpy_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(TEWA ~ clumpy_1 + I(clumpy_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes lower-order terms)
model_cubic <- glm.nb(TEWA ~ clumpy_1 + I(clumpy_1^2) + I(clumpy_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(TEWA ~ clumpy_1 + I(clumpy_1^2) + I(clumpy_1^3) + I(clumpy_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(TEWA ~ exp(clumpy_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$clumpy_1)
NTEMS$above_thresh <- as.numeric(NTEMS$clumpy_1 > threshold)
model_piecewise <- glm.nb(TEWA ~ clumpy_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(TEWA ~ I(clumpy_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

Model      AIC
7  Cubic Root 16955.95
4      Fourth 17454.78
6   Piecewise 17457.51
5 Exponential 17458.92
2   Quadratic 17458.96
3       Cubic 17460.66
1      Linear 17461.05

# Cubic root is best




####### Age

# Linear model
model_lin <- glm.nb(TEWA ~ age_mn_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(TEWA ~ age_mn_1 + I(age_mn_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes all lower-order terms)
model_cubic <- glm.nb(TEWA ~ age_mn_1 + I(age_mn_1^2) + I(age_mn_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(TEWA ~ age_mn_1 + I(age_mn_1^2) + I(age_mn_1^3) + I(age_mn_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(TEWA ~ exp(age_mn_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$age_mn_1)
NTEMS$above_thresh <- as.numeric(NTEMS$age_mn_1 > threshold)
model_piecewise <- glm.nb(TEWA ~ age_mn_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Square Root model
model_sqrt <- glm.nb(TEWA ~ sqrt(age_mn_1) + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(TEWA ~ I(age_mn_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Square Root", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_sqrt), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


Model      AIC
4      Fourth 17152.68
3       Cubic 17169.53
2   Quadratic 17201.27
6   Piecewise 17220.37
1      Linear 17302.26
7 Square Root 17338.14
8  Cubic Root 17350.72
5 Exponential 17437.13

# Fourth is best









#################### BBWA

### Prop_con_1

# Linear model
model_lin <- glm.nb(BBWA ~ prop_con_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(BBWA ~ prop_con_1 + I(prop_con_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes lower-order terms)
model_cubic <- glm.nb(BBWA ~ prop_con_1 + I(prop_con_1^2) + I(prop_con_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(BBWA ~ prop_con_1 + I(prop_con_1^2) + I(prop_con_1^3) + I(prop_con_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(BBWA ~ exp(prop_con_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$prop_con_1)
NTEMS$above_thresh <- as.numeric(NTEMS$prop_con_1 > threshold)
model_piecewise <- glm.nb(BBWA ~ prop_con_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Square Root model
model_sqrt <- glm.nb(BBWA ~ sqrt(prop_con_1) + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(BBWA ~ I(prop_con_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Square Root", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_sqrt), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


> print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)
Model      AIC
6   Piecewise 4462.789
2   Quadratic 4464.419
8  Cubic Root 4465.605
3       Cubic 4465.838
7 Square Root 4466.325
4      Fourth 4467.614
1      Linear 4469.171
5 Exponential 4471.106

# Quadratic is best





### Clumpy_1

# Linear model
model_lin <- glm.nb(BBWA ~ clumpy_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(BBWA ~ clumpy_1 + I(clumpy_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes lower-order terms)
model_cubic <- glm.nb(BBWA ~ clumpy_1 + I(clumpy_1^2) + I(clumpy_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(BBWA ~ clumpy_1 + I(clumpy_1^2) + I(clumpy_1^3) + I(clumpy_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(BBWA ~ exp(clumpy_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$clumpy_1)
NTEMS$above_thresh <- as.numeric(NTEMS$clumpy_1 > threshold)
model_piecewise <- glm.nb(BBWA ~ clumpy_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(BBWA ~ I(clumpy_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)

Model      AIC
7  Cubic Root 4331.157
4      Fourth 4436.374
6   Piecewise 4440.431
2   Quadratic 4444.148
3       Cubic 4446.112
5 Exponential 4452.676
1      Linear 4458.490

# Cubic root is best








### Age

# Linear model
model_lin <- glm.nb(BBWA ~ age_mn_1 + offset(log(survey_effort)), data = NTEMS)

# Quadratic model (includes linear term)
model_quad <- glm.nb(BBWA ~ age_mn_1 + I(age_mn_1^2) + offset(log(survey_effort)), data = NTEMS)

# Cubic model (includes all lower-order terms)
model_cubic <- glm.nb(BBWA ~ age_mn_1 + I(age_mn_1^2) + I(age_mn_1^3) + offset(log(survey_effort)), data = NTEMS)

# Fourth model (includes all lower-order terms)
model_fourth <- glm.nb(BBWA ~ age_mn_1 + I(age_mn_1^2) + I(age_mn_1^3) + I(age_mn_1^4) + offset(log(survey_effort)), data = NTEMS)

# Exponential model
model_exp <- glm.nb(BBWA ~ exp(age_mn_1) + offset(log(survey_effort)), data = NTEMS)

# Piecewise model (Threshold at median value)
threshold <- median(NTEMS$age_mn_1)
NTEMS$above_thresh <- as.numeric(NTEMS$age_mn_1 > threshold)
model_piecewise <- glm.nb(BBWA ~ age_mn_1 * above_thresh + offset(log(survey_effort)), data = NTEMS)

# Square Root model
model_sqrt <- glm.nb(BBWA ~ sqrt(age_mn_1) + offset(log(survey_effort)), data = NTEMS)

# Cubic Root model
model_cuberoot <- glm.nb(BBWA ~ I(age_mn_1^(1/3)) + offset(log(survey_effort)), data = NTEMS)

# Compare AIC values
aic_values <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Fourth", "Exponential", "Piecewise", "Square Root", "Cubic Root"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_fourth),
          AIC(model_exp), AIC(model_piecewise), AIC(model_sqrt), AIC(model_cuberoot))
)

print(aic_values[order(aic_values$AIC), ])  # Sort models by AIC (lower is better)


Model      AIC
2   Quadratic 4325.576
3       Cubic 4327.171
4      Fourth 4329.166
6   Piecewise 4344.216
5 Exponential 4459.558
8  Cubic Root 4461.444
7 Square Root 4464.903
1      Linear 4473.057

# Quadratic is best







# Best models

BTNW (Black-throated Green Warbler)
Prop_con_1 → Cubic
Clumpy_1 → Cubic Root
Age → Cubic


TEWA (Tennessee Warbler)
Prop_con_1 → Linear
Clumpy_1 → Cubic Root
Age → Fourth


BBWA (Bay-breasted Warbler)
Prop_con_1 → Quadratic
Clumpy_1 → Cubic Root
Age → Quadratic


