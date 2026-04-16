#### COX IPW FOR TTE ANALYSIS####
## By Charlie Cunniffe
## 14/04/26


#################################################################################
###    Fundamentals of Confounding Adjustment                                 ### 
###    Hands-on Session 5 - ESTIMATING THE EFFECT OF A POINT INTERVENTION     ###
###                         ON A SURVIVAL OUTCOME                             ### 
###    CAUSALab                                                               ### 
#################################################################################


# Load required packages
if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if (!require("data.table")) install.packages("data.table"); library(data.table)
if (!require("survival")) install.packages("survival"); library(survival)
if (!require("survminer")) install.packages("survminer"); library(survminer)
if (!require("boot")) install.packages("boot"); library(boot)
if (!require("MatchIt")) install.packages("MatchIt"); library(MatchIt)
if (!require("geepack")) install.packages("geepack"); library(geepack)
library(mice)
# Load necessary libraries
library(survival)
library(survey)
library(dplyr)
library(car)
set.seed(2602)

imp_test <- complete(imp2, 30)

########IPW COX ########
# Step 2: Define a function to fit the Cox model with IPW for each imputed dataset
#imputed_data <- imp_test
fit_cox_with_ipw <- function(imputed_data) {
  imputed_data <- as.data.frame(imputed_data) %>% filter(age<80)
  # Estimate the propensity score model for treatment assignment
  ps_model <- glm(treatconcurrent ~ age + #I(age*age) +
                    time_since_study_start_days + #I(time_since_study_start_days*time_since_study_start_days) +
                    #season +
                    imd_average_decile +
                    ecog +
                    smoking_history +
                    tumour_side_left +
                    tumour_size +
                    node_size +
                    cardiac_comor +
                    resp_comor +
                    #renal_comor +
                    other_comor# +
                  #alcohol +
                  #substance_abuse, 
                  ,family = binomial(), data = imputed_data)
  
  #print(vif(ps_model))
  ps_scores <- predict(ps_model, type = "response")
  
  # Compute the inverse probability weights
  ip_weights <- ifelse(imputed_data$treatconcurrent == 1, 1 / ps_scores, 1 / (1 - ps_scores))
  
  # Create a weighted survival design using the survey package
  svy_data <- svydesign(ids = ~1, data = imputed_data, weights = ~ip_weights)
  
  # Fit the Cox model using the weighted data
  cox_model <- svycoxph(Surv(time_days, deceased) ~ treatconcurrent 
                        + age + #I(age*age) +
                          time_since_study_start_days + #I(time_since_study_start_days*time_since_study_start_days) +
                          #season +
                          imd_average_decile +
                          ecog +
                          smoking_history +
                          tumour_side_left +
                          tumour_size +
                          node_size +
                          cardiac_comor +
                          resp_comor +
                          #renal_comor +
                          other_comor# +
                        #alcohol +
                        #substance_abuse
                        , design = svy_data, data = imputed_data)
  #cox_model<-coxph(Surv(time_days, deceased)~treatconcurrent, weights = ip_weights, data = imputed_data )
  
  return(cox_model)
}
# Step 3: Apply the function to each imputed dataset (correct approach)
# Use the 'complete' function to extract imputed datasets from 'imp2'
cox_results <- lapply(1:imp2$m, function(i) {
  # Extract the i-th imputed dataset
  imputed_data <- complete(imp2, i)
  
  # Apply the function to the imputed dataset
  fit_cox_with_ipw(imputed_data)
})

# Now, cox_results is a list of Cox models for each imputed dataset
fit_cox_with_ipw(imp_test)
# Step 4: Pool the results from the multiple imputations
# Use the 'mitools' or 'mice' package to pool results
# Assuming 'cox_results' contains a list of coxph objects:

pooled_results <- pool(cox_results)

# Display the pooled results
summary(pooled_results)
