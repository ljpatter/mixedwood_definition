# ---
# title: "Mixed regression models"
# author: "Brendan Casey"
# created: ""
# description: "Building mixed regression models"
# ---

#Setup ----

##Load packages----
library(lme4)
library(MuMIn)
library(fitdistrplus) #checking dsistribution of response variable
library(tidyverse)
library(ggplot2)
library(ggeffects)


# Load data
NTEMS <- read.csv("Output/Tabular Data/NTEMS_w_clustering.csv")


##Import data----
# load("0_data/manual/bd_cov_scaled_2024-05-03.rData")
load("0_data/manual/bd_cov_scaled_grouped_2024-05-07.rData")
d<-bd_cov_scaled_grouped%>%
  mutate(spatial_group=as.factor(spatial_group))


##////////////////////////////////////////////////////////////////

#Shannon diversity ----

##Evaluate nonlinear effects----
### Fractional mixedwood----
nl1<-lme4::lmer(shan ~prcB_1000_mean   + (1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(shan ~poly(prcB_1000_mean,2)   + (1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(shan ~poly(prcB_1000_mean,3)   + (1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(shan ~poly(prcB_1000_mean,4)   + (1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$prcB_1000_mean, 'scaled:scale')
ce_1 <- attr(d$prcB_1000_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

# Find the x-coordinate corresponding to the peak of the trend line
peak_x <- dd2$x[which.max(dd2$predicted)]

# Create the plot
plot <- ggplot(dd2) +
  aes(x, predicted, ymin = conf.low, ymax = conf.high) +
  # Add ribbon for confidence interval
  geom_ribbon(alpha = 0.3) +  
  # Add line for the predicted values
  geom_line() +  
  # Add axis labels
  labs(x = "Fractional Broadleaf", y = "Shannon diversity" ) +  
  # Apply minimal theme
  theme_minimal() +  
  theme(# Remove major gridlines
    panel.grid.major = element_blank(),  
    # Remove minor gridlines
    panel.grid.minor = element_blank(), 
    # Set panel background color to white
    panel.background = element_rect(fill = "white"), 
    # Set axis line color to black
    axis.line = element_line(colour = "black")) + 
  # Add vertical line
  geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")  

# Save the plot
ggsave("3_output/figures/shan_fractional_broadleaf.png", plot, 
       width = 6, height = 4, dpi = 300)


#^2

### Forest age----
nl1<-lme4::lmer(shan ~age_mean_500   + (1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(shan ~poly(age_mean_500,2)   + (1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(shan ~poly(age_mean_500,3)   + (1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(shan ~poly(age_mean_500,4)   + (1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$age_mean_500, 'scaled:scale')
ce_1 <- attr(d$age_mean_500, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("age_mean_500 [all]"))%>%
  as_tibble()%>%
  mutate(x = x* sc_1 + ce_1)

ggplot(dd2)+
  aes(x, predicted, ymin = conf.low, ymax = conf.high)+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = TRUE)

#^2

### Height upper----
nl1<-lme4::lmer(shan ~height_500_mean   + (1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(shan ~poly(height_500_mean,2)   + (1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(shan ~poly(height_500_mean,3)   + (1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(shan ~poly(height_500_mean,4)   + (1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$height_500_mean, 'scaled:scale')
ce_1 <- attr(d$height_500_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("height_500_mean [all]"))%>%
  as_tibble()%>%
  mutate(x = x* sc_1 + ce_1)

ggplot(dd2)+
  aes(x, predicted, ymin = conf.low, ymax = conf.high)+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = TRUE)

#^2
##////////////////////////////////////////////////////////////////

## 150----
shan_m_1_150<-lme4::lmer(shan ~  poly(age_mean_150,2, raw=TRUE) 
                         * poly(prcB_150_mean,2, raw=TRUE) 
                         + poly(height_150_mean,2, raw=TRUE)
                         + (1|spatial_group/location)
                         + offset(log(survey_effort)),
                         data=d,  na.action = na.fail)


shan_dredge_150<-MuMIn::dredge(shan_m_1_150,
                               fixed = c("offset(log(survey_effort))",
                                         "(1|spatial_group/location)"),                    
                               extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))


shan_m_2_150 <- lme4::lmer(shan ~ poly(prcB_150_mean,2, raw=TRUE) 
                           + poly(height_150_mean,2, raw=TRUE)
                           + (1|spatial_group/location)
                           + offset(log(survey_effort)),
                           data=d,  na.action = na.fail)

shan_m_2_150_boot <- bootMer(shan_m_2_150, function(x) fixef(x), nsim = 100)


## 500----
shan_m_1_500<-lme4::lmer(shan ~  poly(age_mean_500,2, raw=TRUE) 
                         * poly(prcB_500_mean,2, raw=TRUE) 
                         + poly(height_500_mean,2, raw=TRUE)
                         + (1|spatial_group/location)
                         + offset(log(survey_effort)),
                         data=d,  na.action = na.fail)


shan_dredge_500<-MuMIn::dredge(shan_m_1_500,
                               fixed = c("offset(log(survey_effort))",
                                         "(1|spatial_group/location)"),                    
                               extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

shan_m_2_500 <- lme4::lmer(shan ~  poly(age_mean_500,2, raw=TRUE) 
                           * poly(prcB_500_mean,2, raw=TRUE) 
                           + poly(height_500_mean,2, raw=TRUE)
                           + (1|spatial_group/location)
                           + offset(log(survey_effort)),
                           data=d,  na.action = na.fail)

shan_m_2_500_boot <- bootMer(shan_m_2_500, function(x) fixef(x), nsim = 100)

## 1000----
shan_m_1_1000<-lme4::lmer(shan ~  poly(age_mean_1000,2, raw=TRUE) 
                          * poly(prcB_1000_mean,2, raw=TRUE) 
                          + poly(height_1000_mean,2, raw=TRUE)
                          + (1|spatial_group/location)
                          + offset(log(survey_effort)),
                          data=d,  na.action = na.fail)


shan_dredge_1000<-MuMIn::dredge(shan_m_1_500,
                                fixed = c("offset(log(survey_effort))",
                                          "(1|spatial_group/location)"),                    
                                extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

shan_m_2_1000 <- lme4::lmer(shan ~  poly(age_mean_1000,2, raw=TRUE) 
                            * poly(prcB_1000_mean,2, raw=TRUE) 
                            + poly(height_1000_mean,2, raw=TRUE)
                            + (1|spatial_group/location)
                            + offset(log(survey_effort)),
                            data=d,  na.action = na.fail)

shan_m_2_1000_boot <- bootMer(shan_m_2_1000, function(x) fixef(x), nsim = 100)

## Save----
save(shan_dredge_1000, shan_dredge_150, shan_dredge_500, 
     shan_m_1_1000, shan_m_1_150, shan_m_1_500, shan_m_2_1000, 
     shan_m_2_1000_boot, shan_m_2_150, shan_m_2_150_boot, 
     shan_m_2_500, shan_m_2_500_boot, 
     file = "3_output/models/shan_models.RData")

##////////////////////////////////////////////////////////////////
#Richness ----

##Evaluate nonlinear effects----
### Fractional mixedwood
nl1<-lme4::glmer(richness ~prcB_1000_mean + (1|location),
                 data=d, family=poisson,  na.action = na.exclude) 

nl2<-lme4::glmer(richness ~poly(prcB_500_mean,2) + (1|location), 
                 data=d, family=poisson, na.action = na.exclude, 
                 control=glmerControl(optCtrl=list(maxfun=1e9))) 

nl3<-lme4::glmer(richness ~poly(prcB_1000_mean,3) + (1|location), 
                 data=d, family=poisson,  na.action = na.exclude, 
                 control=glmerControl(optCtrl=list(maxfun=1e9))) 

nl4<-lme4::glmer(richness ~poly(prcB_1000_mean,4) + (1|location), 
                 data=d,family=poisson,  na.action = na.exclude, 
                 control=glmerControl(optCtrl=list(maxfun=1e9))) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$prcB_500_mean, 'scaled:scale')
ce_1 <- attr(d$prcB_500_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("prcB_500_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

dd3<-ggpredict(nl3, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

dd4<-ggpredict(nl4, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)
# Find the x-coordinate corresponding to the peak of the trend line
peak_x <- dd2$x[which.max(dd2$predicted)]

# Create the plot
plot <- ggplot(dd2) +
  aes(x, predicted, ymin = conf.low, ymax = conf.high) +
  # Add ribbon for confidence interval
  geom_ribbon(alpha = 0.3) +  
  # Add line for the predicted values
  geom_line() +  
  # Add axis labels
  labs(x = "Fractional Broadleaf", y = "Richness") +  
  # Apply minimal theme
  theme_minimal() +  
  theme(# Remove major gridlines
    panel.grid.major = element_blank(),  
    # Remove minor gridlines
    panel.grid.minor = element_blank(), 
    # Set panel background color to white
    panel.background = element_rect(fill = "white"), 
    # Set axis line color to black
    axis.line = element_line(colour = "black")) + 
  # Add vertical line
  geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")  

# Save the plot
ggsave("3_output/figures/richness_fractional_broadleaf.png", plot, width = 6, 
       height = 4, dpi = 300)

##////////////////////////////////////////////////////////////////

## 150----
richness_m_1_150<-lme4::glmer(richness ~  poly(age_mean_150,2, raw=TRUE) 
                              * poly(prcB_150_mean,2, raw=TRUE) 
                              + poly(height_150_mean,2, raw=TRUE)
                              + (1|spatial_group/location)
                              + offset(log(survey_effort)),
                              data=d, family=poisson, na.action = na.fail,
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_dredge_150<-MuMIn::dredge(richness_m_1_150,
                                   fixed = c("offset(log(survey_effort))",
                                             "(1|spatial_group/location)"),                    
                                   extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

richness_m_2_150 <-lme4::glmer(richness ~  poly(age_mean_150,2, raw=TRUE) 
                               * poly(prcB_150_mean,2, raw=TRUE) 
                               + poly(height_150_mean,2, raw=TRUE)
                               + (1|spatial_group/location)
                               + offset(log(survey_effort)),
                               data=d, family=poisson, na.action = na.fail)

richness_m_2_150_boot <- bootMer(richness_m_2_150, function(x) fixef(x), nsim = 100)

## 500----
richness_m_1_500<-lme4::glmer(richness ~  poly(age_mean_500,2, raw=TRUE) 
                              * poly(prcB_500_mean,2, raw=TRUE) 
                              + poly(height_500_mean,2, raw=TRUE)
                              + (1|spatial_group/location)
                              + offset(log(survey_effort)),
                              data=d, family=poisson, na.action = na.fail,
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_dredge_500<-MuMIn::dredge(richness_m_1_500,
                                   fixed = c("offset(log(survey_effort))",
                                             "(1|spatial_group/location)"),                    
                                   extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

richness_m_2_500 <- lme4::glmer(richness ~  poly(age_mean_500,2, raw=TRUE) 
                                + poly(prcB_500_mean,2, raw=TRUE) 
                                + poly(height_500_mean,2, raw=TRUE)
                                + (1|spatial_group/location)
                                + offset(log(survey_effort)),
                                data=d, family=poisson, na.action = na.fail,
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_m_2_500_boot <- bootMer(richness_m_2_500, function(x) fixef(x), nsim = 100)

## 1000----
richness_m_1_1000<-lme4::glmer(richness ~  poly(age_mean_1000,2, raw=TRUE) 
                               * poly(prcB_1000_mean,2, raw=TRUE) 
                               + poly(height_1000_mean,2, raw=TRUE)
                               + (1|spatial_group/location)
                               + offset(log(survey_effort)),
                               data=d, family=poisson, na.action = na.fail,
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_dredge_1000<-MuMIn::dredge(richness_m_1_500,
                                    fixed = c("offset(log(survey_effort))",
                                              "(1|spatial_group/location)"),                    
                                    extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

richness_m_2_1000 <- lme4::glmer(richness ~  poly(age_mean_1000,2, raw=TRUE) 
                                 * poly(prcB_1000_mean,2, raw=TRUE) 
                                 + poly(height_1000_mean,2, raw=TRUE)
                                 + (1|spatial_group/location)
                                 + offset(log(survey_effort)),
                                 data=d, family=poisson, na.action = na.fail)	
richness_m_2_1000_boot <- bootMer(richness_m_2_1000, function(x) fixef(x), nsim = 100)

## Save----
save(richness_dredge_1000, richness_dredge_150, richness_dredge_500, 
     richness_m_1_1000, richness_m_1_150, richness_m_1_500, 
     richness_m_2_1000, richness_m_2_1000_boot, richness_m_2_150, 
     richness_m_2_150_boot, richness_m_2_500, richness_m_2_500_boot,
     file = "3_output/models/richness_models.RData")

##////////////////////////////////////////////////////////////////
# Evenness ----

##Evaluate nonlinear effects----
### Fractional mixedwood
nl1<-lme4::lmer(pilou_even ~prcB_1000_mean  + (1|spatial_group)+(1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(pilou_even ~poly(prcB_1000_mean,2)+ (1|spatial_group)+(1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(pilou_even ~poly(prcB_1000_mean,3)+ (1|spatial_group)+(1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(pilou_even ~poly(prcB_1000_mean,4)+ (1|spatial_group)+(1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$prcB_1000_mean, 'scaled:scale')
ce_1 <- attr(d$prcB_1000_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl1, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

dd3<-ggpredict(nl3, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

# Find the x-coordinate corresponding to the peak of the trend line
peak_x <- dd2$x[which.min(dd2$predicted)]

# Create the plot
plot <- ggplot(dd2) +
  aes(x, predicted, ymin = conf.low, ymax = conf.high) +
  # Add ribbon for confidence interval
  geom_ribbon(alpha = 0.3) +  
  # Add line for the predicted values
  geom_line() +  
  # Add axis labels
  labs(x = "Fractional Broadleaf", y = "Pielou's evenness") +  
  # Apply minimal theme
  theme_minimal() +  
  theme(# Remove major gridlines
    panel.grid.major = element_blank(),  
    # Remove minor gridlines
    panel.grid.minor = element_blank(), 
    # Set panel background color to white
    panel.background = element_rect(fill = "white"), 
    # Set axis line color to black
    axis.line = element_line(colour = "black")) + 
  # Add vertical line
  geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")  

# Save the plot
ggsave("3_output/figures/pilou_even_fractional_broadleaf.png", plot, width = 6, 
       height = 4, dpi = 300)

##////////////////////////////////////////////////////////////////

## 150----
pilou_even_m_1_150<-lme4::lmer(pilou_even ~  poly(age_mean_150,2, raw=TRUE) 
                               * poly(prcB_150_mean,2, raw=TRUE) 
                               + poly(height_150_mean,2, raw=TRUE)
                               + (1|spatial_group/location)
                               + offset(log(survey_effort)),
                               data=d,  na.action = na.fail,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_dredge_150<-MuMIn::dredge(pilou_even_m_1_150,
                                     fixed = c("offset(log(survey_effort))",
                                               "(1|spatial_group/location)"),                    
                                     extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

pilou_even_m_2_150 <- lme4::lmer(pilou_even ~  poly(age_mean_150,2, raw=TRUE) 
                                 + prcB_150_mean
                                 + poly(height_150_mean,2, raw=TRUE)
                                 + (1|location)
                                 + offset(log(survey_effort)),
                                 data=d,  na.action = na.fail,
                                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

pilou_even_m_2_150_boot <- bootMer(pilou_even_m_2_150, function(x) fixef(x), nsim = 100)

## 500----
pilou_even_m_1_500<-lme4::lmer(pilou_even ~  poly(age_mean_500,2, raw=TRUE) 
                               + poly(prcB_500_mean,2, raw=TRUE) 
                               + poly(height_500_mean,2, raw=TRUE)
                               + (1|spatial_group)
                               + offset(log(survey_effort)),
                               data=d,  na.action = na.fail,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_dredge_500<-MuMIn::dredge(pilou_even_m_1_500,
                                     fixed = c("offset(log(survey_effort))",
                                               "(1|spatial_group/location)"),                    
                                     extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

pilou_even_m_2_500 <- lme4::lmer(pilou_even ~  poly(age_mean_500,2, raw=TRUE) 
                                 + prcB_500_mean
                                 + poly(height_500_mean,2, raw=TRUE)
                                 + (1|location)
                                 + offset(log(survey_effort)),
                                 data=d,  na.action = na.fail,
                                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_m_2_500_boot <- bootMer(pilou_even_m_2_500, function(x) fixef(x), nsim = 100)

## 1000----
pilou_even_m_1_1000<-lme4::lmer(pilou_even ~  poly(age_mean_1000,2, raw=TRUE) 
                                * poly(prcB_1000_mean,2, raw=TRUE) 
                                + poly(height_1000_mean,2, raw=TRUE)
                                + (1|spatial_group)
                                + offset(log(survey_effort)),
                                data=d,  na.action = na.fail,
                                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_dredge_1000<-MuMIn::dredge(pilou_even_m_1_500,
                                      fixed = c("offset(log(survey_effort))",
                                                "(1|spatial_group/location)"),                    
                                      extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

pilou_even_m_2_1000 <- lme4::lmer(pilou_even ~  poly(age_mean_1000,2, raw=TRUE) 
                                  + prcB_1000_mean
                                  + poly(height_1000_mean,2, raw=TRUE)
                                  + (1|location)
                                  + offset(log(survey_effort)),
                                  data=d,  na.action = na.fail,
                                  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))



pilou_even_m_2_1000_boot <- bootMer(pilou_even_m_2_1000, function(x) fixef(x), nsim = 100)


## Save----
# save(pilou_even_dredge_1000, pilou_even_dredge_150, 
#      pilou_even_dredge_500, pilou_even_m_1_1000, 
#      pilou_even_m_1_150, pilou_even_m_1_500, 
#      pilou_even_m_2_1000, pilou_even_m_2_1000_boot, 
#      pilou_even_m_2_150, pilou_even_m_2_150_boot, 
#      pilou_even_m_2_500, pilou_even_m_2_500_boot,
#      file = "3_output/models/pilou_even_models.RData")

save(   pilou_even_m_2_1000,  
        pilou_even_m_2_150, 
        pilou_even_m_2_500, 
        # pilou_even_m_2_150_boot, 
        # pilou_even_m_2_500_boot,
        # pilou_even_m_2_1000_boot,
        # pilou_even_dredge_1000, pilou_even_dredge_150, 
        # pilou_even_dredge_500, pilou_even_m_1_1000, 
        # pilou_even_m_1_150, pilou_even_m_1_500, 
        
        file = "3_output/models/pilou_even_models.RData")

##////////////////////////////////////////////////////////////////



library(caret)
# Define the control parameters for cross-validation
ctrl <- trainControl(method = "cv", number = 3)

# Perform cross-validation
cv_results <- train(nl2, method = "lmer", trControl = ctrl)

cv_results <- bootMer(nl2, FUN = fixef, nsim = 3)  # Use appropriate number of nsim

boot_gm1 <- bootMer(gm1, FUN = function(x) fixef(x), nsim = 100)

pilou_even_m_2_150_boot <- bootMer(pilou_even_m_2_150, function(x) fixef(x), nsim = 100)
