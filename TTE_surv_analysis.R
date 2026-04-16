#### SURVIVAL ANALYSIS FOR TTE ANALYSIS####
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


#######
imp_test <- imp_test %>% filter(age<80 & age>45)

#################################################################################
###                                                                           ###
###     DATA ANALYSIS EXERCISE 1                                              ###
###                                                                           ###
#################################################################################

# Note: many of the data manipulation procedures can be accomplished using 
# tidyverse or data.table. Where possible, both approaches are provided below, 
# and will generate equivalent results. However, data.table will generally be 
# more computationally efficient with larger datasets.    

################################
### a) Creating survtime variable that represents the duration (in months) of follow-up 
### and then expanding the dataset using the survtime variable
imp_test$survtime <- round(imp_test$time_days / 30, digits = 0)

imp_test.surv <- uncount(imp_test, weights=survtime, .remove=F)

################################
### b) Creating variables for time
imp_test.surv <- imp_test.surv %>% group_by(patient_id) %>% dplyr::mutate(time=row_number()-1) %>% ungroup()
# Creating variable for timesq
imp_test.surv$timesq <- imp_test.surv$time^2

################################
### c) Creating event variable
imp_test.surv$event <- ifelse(
  imp_test.surv$time==imp_test.surv$survtime-1 & imp_test.surv$deceased==1, 1, 0)



#################################################################################
###                                                                           ###
###     DATA ANALYSIS EXERCISE 2                                              ###
###                                                                           ###
#################################################################################

################################
### a) Fitting a pooled logistic regression model with time, treatment and 
###    product term between treatment and time

fit.pool1 <- glm(
  formula = event==1 ~ treatconcurrent + time + treatconcurrent*time,
  family = binomial(link = 'logit'),
  data = imp_test.surv)

summary(fit.pool1)

################################
### b) Fitting a pooled logistic regression model with time (linear and 
###    quadratic terms), treatment and product terms between treatment and time  

fit.pool2 <- glm(
  formula = event==1 ~ treatconcurrent + time + timesq + treatconcurrent*time + treatconcurrent*timesq,
  family = binomial(link = 'logit'),
  data = imp_test.surv)
summary(fit.pool2)

################################
### c) Computing risks at each time point of follow-up

# Create datasets to store results
# Include all time points under each treatment level
results0 <- data.frame(treatconcurrent=0, time=seq(0,119), timesq=(seq(0,119))^2)
results1 <- data.frame(treatconcurrent=1, time=seq(0,119), timesq=(seq(0,119))^2)

# Obtain predicted hazards from pooled logistic regression model
results0$hazard0 <- predict(fit.pool2, results0, type="response")
results1$hazard1 <- predict(fit.pool2, results1, type="response")

# Estimate survival probabilities from hazards
# S(t) = cumulative product of (1 - h(t))
results0$surv0 <- cumprod(1-results0$hazard0)
results1$surv1 <- cumprod(1-results1$hazard1)

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
results0$risk0 <- 1 - results0$surv0
results1$risk1 <- 1 - results1$surv1

#View(results0)
#View(results1)

################################
### d) Constructing risk curves

# Combine results for each treatment group into a single dataset
results.combined <- merge(results0, results1, by=c("time", "timesq"))

# Create a new "time" variable to reflect the fact that risks start at 0
# and are estimated at the end of each time interval
results.combined$time_updated <- results.combined$time + 1
results.combined <- results.combined %>% add_row(time_updated=0, risk0=0, risk1=0) %>% arrange(time_updated)

# Creating plot
ggplot(results.combined, aes(x=time_updated)) + 
  geom_line(aes(y = risk0, colour = "0")) + 
  geom_line(aes(y = risk1, colour = "1")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
  ylab("Cumulative Incidence") + 
  labs(colour="Concurrent treatment:") +
  ggtitle("Unadjusted cummulative risk curve") +
  theme_bw() + 
  theme(legend.position="bottom")

ggplot(results.combined, aes(x=time_updated)) + 
  geom_line(aes(y = 1-risk0, colour = "0")) + 
  geom_line(aes(y = 1-risk1, colour = "1")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
  ylab("Cumulative Incidence") + 
  labs(colour="Concurrent Treatment:") +
  ggtitle("Unadjusted inverse cummulative risk curve")+
  theme_bw() + 
  theme(legend.position="bottom")

survplot <- ggsurvplot(survfit(Surv(survtime, deceased)~treatconcurrent, data = imp_test), data = imp_test,
           ggtheme = theme_bw(), conf.int = TRUE)
survplot$plot + #ggplot(results.combined, aes(x=time_updated)) + 
  geom_line(data = results.combined, aes(y = 1-results.combined$risk0, x= results.combined$time_updated, colour = "0")) + 
  geom_line(data = results.combined, aes(y = 1-results.combined$risk1, x= results.combined$time_updated, colour = "1")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
  ylab("Cumulative Incidence") + 
  labs(colour="Concurrent Treatment:") +
  ggtitle("Unadjusted inverse cummulative risk curve")+
  theme_bw() + 
  theme(legend.position="bottom")

################################
### e) Use the risks you calculated in (c) to obtain the risk at end of follow-up

risk1 <- results1[results1$time==119,]$risk1
risk0 <- results0[results0$time==119,]$risk0
risk1
risk0

################################
### f) Use the risks to estimate the effect on the risk ratio and risk
### difference scales

risk_dif_10yr <- risk1 - risk0
risk_ratio_10yr <- risk1 / risk0

################################
### i) Optional: pooled logistic regression and compare results against a Cox 
###    proportional hazards model

# Fitting pooled logistic regression model with treatment and time (linear 
# and quadratic terms)
fit.pool3 <- glm(
  formula = event==1 ~ treatconcurrent + time + timesq,
  family = binomial(link = 'logit'),
  data = imp_test.surv)

summary(fit.pool3)
exp(fit.pool3$coefficients)

# Cox proportional hazards model for comparison
fit.cox <- coxph(Surv(survtime, deceased) ~ treatconcurrent,
                 data = imp_test)
summary(fit.cox)



#################################################################################
###                                                                           ###
###     DATA ANALYSIS EXERCISE 3                                              ###
###                                                                           ###
#################################################################################

################################
### a) Fitting PLR

plr.model <- glm(
  event ~ treatconcurrent + 
    age + I(age*age) +
    time_since_study_start_days + I(time_since_study_start_days*time_since_study_start_days) +
    season +
    imd_average_decile +
    ecog +
    smoking_history +
    tumour_side_left +
    tumour_size +
    node_size +
    cardiac_comor +
    resp_comor +
    #renal_comor +
    other_comor +
    time + I(time*time) + 
    I(treatconcurrent*time) + I(treatconcurrent*time*time), 
  data=imp_test.surv, 
  family=binomial())

summary(plr.model)


################################
### b) Calculating risks at each time point, standardizing and standardized risks

# Creating dataset with all time points for each individual under each 
# treatment level

# Note: everyone will have 120 rows of data, regardless of whether or not they
# developed the outcome in the original dataset. These datasets will be used 
# to store and calculate final results

# Had everyone been untreated
surv.results0 <- uncount(imp_test, weights=120, .remove=F)
surv.results0 <- surv.results0 %>% group_by(patient_id) %>% mutate(time=row_number()-1) %>% ungroup() 
surv.results0$treatconcurrent <- 0

# Had everyone been treated
surv.results1 <- surv.results0
surv.results1$treatconcurrent <- 1

# Calculating risks from hazards

# Hazards based on predicted probabilities from PLR
surv.results0$hazard0 <- predict(plr.model, newdata=surv.results0, type="response")  
surv.results1$hazard1 <- predict(plr.model, newdata=surv.results1, type="response")  

# Survival from cumulative product of (1-hazard) for each individual
surv.results0 <- arrange(surv.results0, patient_id, time)
surv.results1 <- arrange(surv.results1, patient_id, time)

surv.results0 <- surv.results0 %>% group_by(patient_id) %>% mutate(surv0 = cumprod(1-hazard0)) 
surv.results1 <- surv.results1 %>% group_by(patient_id) %>% mutate(surv1 = cumprod(1-hazard1))

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
surv.results0$risk0 <- 1 - surv.results0$surv0
surv.results1$risk1 <- 1 - surv.results1$surv1

# Standardization by averaging
# Note: we calculate the averages stratified by time point


surv.results0 <- surv.results0 %>%
  group_by(time) %>%
  summarize(
    treatconcurrent = first(treatconcurrent),  # assuming treatconcurrent is constant within each time point
    surv0 = mean(surv0, na.rm = TRUE),        # calculate mean survival probability
    risk0 = mean(risk0, na.rm = TRUE)         # calculate mean risk
  )

surv.results1 <- surv.results1 %>%
  group_by(time) %>%
  summarize(
    treatconcurrent = first(treatconcurrent),  # assuming treatconcurrent is constant within each time point
    surv1 = mean(surv1, na.rm = TRUE),        # calculate mean survival probability
    risk1 = mean(risk1, na.rm = TRUE)         # calculate mean risk
  )




################################
### c) Risks at end of follow-up
risk1 <- surv.results1[surv.results1$time==119,]$risk1
risk0 <- surv.results0[surv.results0$time==119,]$risk0
risk1
risk0

################################
### d) Risk ratio and risk difference
risk_ratio_10yr_adj <- risk1 / risk0
risk_dif_10yr_adj <- risk1 - risk0


### !!!!!!! here
# 1. Compute the Cumulative Hazard for each group
surv.results0$cum_hazard0 <- -log(surv.results0$surv0)
surv.results1$cum_hazard1 <- -log(surv.results1$surv1)

# 2. Compute the Hazard Ratio (HR) at each time point
surv.graph <- merge(surv.results0, surv.results1, by = "time")
surv.graph$HR <- surv.graph$cum_hazard1 / surv.graph$cum_hazard0

# 3. Calculate the Average Hazard Ratio (e.g., for all time points)
average_HR <- mean(surv.graph$HR, na.rm = TRUE)




################################
### f)	Use bootstrapping to compute a valid 95% confidence interval for the causal risk ratio
#imp_test <- imp_test %>% mutate(across(where(~length(unique(.))<13), ~as.factor(.)))
survival.std.boot <- function(data, indices) {
  d <- data[indices,]
  
  # Expanding dataset in person-time format
  d$count <- as.numeric(rownames(d))
  d.gf <- data.frame(count=rep(d$count, times=d$survtime))
  d.gf <- merge(d.gf, d, by="count")  
  d.gf$time <- sequence(rle(d.gf$count)$lengths)-1
  d.gf$event <- ifelse(d.gf$time==d.gf$survtime-1 & 
                         d.gf$deceased==1, 1, 0)
  d.gf$timesq <- d.gf$time^2  
  
  # Fitting PLR with covariates
  boot.model <- glm(
    event ~ treatconcurrent + 
      age + I(age*age) +
      time_since_study_start_days + I(time_since_study_start_days*time_since_study_start_days) +
      season +
      imd_average_decile +
      ecog +
      smoking_history +
      tumour_side_left +
      tumour_size +
      node_size +
      cardiac_comor +
      resp_comor +
      #renal_comor +
      other_comor +
      time + I(time*time) + 
      I(treatconcurrent*time) + I(treatconcurrent*time*time), 
    data=d.gf, 
    family=binomial())
  
  # Creating  dataset with all time points for each individual under 
  # each treatment level
  d.gf.pred <- data.frame(count=rep(d$count, times=120))
  d.gf.pred <- merge(d.gf.pred, d, by="count")    
  d.gf.pred$time <- rep(seq(0, 119), nrow(d))
  d.gf.pred$timesq <- d.gf.pred$time^2
  
  gf.treatconcurrent0 <- gf.treatconcurrent1 <- d.gf.pred
  gf.treatconcurrent0$treatconcurrent <- 0
  gf.treatconcurrent1$treatconcurrent <- 1
  
  # Assigning of estimated hazard to each person-month */
  gf.treatconcurrent0$hazard0 <- predict(boot.model, gf.treatconcurrent0, type="response")
  gf.treatconcurrent1$hazard1 <- predict(boot.model, gf.treatconcurrent1, type="response")  
  
  # Computing risks from survival for each person-month
  gf.treatconcurrent0.surv <- gf.treatconcurrent0 %>% group_by(count) %>% mutate(risk0 = (1-cumprod(1-hazard0)))
  gf.treatconcurrent1.surv <- gf.treatconcurrent1 %>% group_by(count) %>% mutate(risk1 = (1-cumprod(1-hazard1)))
  
  # Calculating average risk for each time point
  gf.surv0 <- aggregate(gf.treatconcurrent0.surv, by=list(gf.treatconcurrent0.surv$time), FUN=mean)[c("treatconcurrent", "time", "risk0")]
  gf.surv1 <- aggregate(gf.treatconcurrent1.surv, by=list(gf.treatconcurrent1.surv$time), FUN=mean)[c("treatconcurrent", "time", "risk1")]
  
  # Data management to output risks, risk differences and risk ratios
  gf.graph <- merge(gf.surv0, gf.surv1, by=c("time"))
  gf.graph$riskdiff <- gf.graph$risk1-gf.graph$risk0
  gf.graph$riskratio <- gf.graph$risk1/gf.graph$risk0
  
  gf.graph <- gf.graph[order(gf.graph$time),]
  
  risk1 <- gf.graph[gf.graph$time==119,]$risk1
  risk0 <- gf.graph[gf.graph$time==119,]$risk0
  
  risk_dif_10yr_adj <- risk1 - risk0
  risk_ratio_10yr_adj <- risk1 / risk0
  
  return(c(risk0, risk1, risk_dif_10yr_adj, risk_ratio_10yr_adj)) 
  
}

# survival.std.boot <- function(data, indices) { 
#   repeat {
#     d <- data[indices, ]  # Create the resample
#     
#     # Check if any factor variable has less than 2 unique levels
#     factor_vars <- sapply(d, is.factor)
#     problematic_factors <- sapply(d[, factor_vars, drop = FALSE], function(x) length(unique(x)) < 2)
#     
#     # If any factor variable has fewer than 2 levels, resample again
#     if (all(!problematic_factors)) {
#       break  # Exit the loop if all factors have at least 2 levels
#     }
#   }
#   
#   # Expanding dataset in person-time format
#   d$count <- as.numeric(rownames(d))
#   d.gf <- data.frame(count=rep(d$count, times=d$survtime))
#   d.gf <- merge(d.gf, d, by="count")
#   d.gf$time <- sequence(rle(d.gf$count)$lengths)-1
#   d.gf$event <- ifelse(d.gf$time == d.gf$survtime-1 & d.gf$deceased == 1, 1, 0)
#   d.gf$timesq <- d.gf$time^2  
#   
#   # Fitting PLR with covariates
#   boot.model <- glm(
#     event ~ treatconcurrent + age + I(age*age) + time_since_study_start_days + 
#       I(time_since_study_start_days*time_since_study_start_days) + season +
#       imd_average_decile + ecog + smoking_history + tumour_side_left + 
#       tumour_size + node_size + cardiac_comor + resp_comor + renal_comor + 
#       other_comor + time + I(time*time) + 
#       I(treatconcurrent*time) + I(treatconcurrent*time*time), 
#     data=d.gf, 
#     family=binomial()
#   )
#   
#   # Creating dataset with all time points for each individual under 
#   # each treatment level
#   d.gf.pred <- data.frame(count=rep(d$count, times=120))
#   d.gf.pred <- merge(d.gf.pred, d, by="count")
#   d.gf.pred$time <- rep(seq(0, 119), nrow(d))
#   d.gf.pred$timesq <- d.gf.pred$time^2
#   
#   gf.treatconcurrent0 <- gf.treatconcurrent1 <- d.gf.pred
#   gf.treatconcurrent0$treatconcurrent <- 0
#   gf.treatconcurrent1$treatconcurrent <- 1
#   
#   # Assigning of estimated hazard to each person-month 
#   gf.treatconcurrent0$hazard0 <- predict(boot.model, gf.treatconcurrent0, type="response")
#   gf.treatconcurrent1$hazard1 <- predict(boot.model, gf.treatconcurrent1, type="response")  
#   
#   # Computing risks from survival for each person-month
#   gf.treatconcurrent0.surv <- gf.treatconcurrent0 %>% group_by(count) %>% mutate(risk0 = (1 - cumprod(1 - hazard0)))
#   gf.treatconcurrent1.surv <- gf.treatconcurrent1 %>% group_by(count) %>% mutate(risk1 = (1 - cumprod(1 - hazard1)))
#   
#   # Calculating average risk for each time point
#   gf.surv0 <- aggregate(gf.treatconcurrent0.surv, by=list(gf.treatconcurrent0.surv$time), FUN=mean)[c("treatconcurrent", "time", "risk0")]
#   gf.surv1 <- aggregate(gf.treatconcurrent1.surv, by=list(gf.treatconcurrent1.surv$time), FUN=mean)[c("treatconcurrent", "time", "risk1")]
#   
#   # Data management to output risks, risk differences and risk ratios
#   gf.graph <- merge(gf.surv0, gf.surv1, by=c("time"))
#   gf.graph$riskdiff <- gf.graph$risk1 - gf.graph$risk0
#   gf.graph$riskratio <- gf.graph$risk1 / gf.graph$risk0
#   
#   gf.graph <- gf.graph[order(gf.graph$time),]
#   
#   risk1 <- gf.graph[gf.graph$time == 119,]$risk1
#   risk0 <- gf.graph[gf.graph$time == 119,]$risk0
#   risk_dif_10yr_adj <- risk1 - risk0
#   risk_ratio_10yr_adj <- risk1 / risk0
#   
#   return(c(risk0, risk1, risk_dif_10yr_adj, risk_ratio_10yr_adj))
# }
# 5 bootstrap samples for demonstration purposes 
set.seed(5435305)
survival.std.results <- boot(data=imp_test, statistic=survival.std.boot, R=50)

boot.ci(survival.std.results,
        conf = 0.95,
        type = "perc",
        index = 4)


################################
### g) Construct marginal cumulative incidence (risk) curves for all-cause mortality, by treatment group.

# Data processing before plotting
surv.graph <- merge(surv.results0, surv.results1, by=c("time"))
surv.graph$time_updated <- surv.graph$time + 1
surv.graph <- surv.graph %>% add_row(time_updated = 0, risk0 = 0, risk1 = 0) %>% arrange(time_updated)

# Plot
ggplot(surv.graph, aes(x=time_updated, y=risk)) + 
  geom_line(aes(y = risk0, colour = "0")) + 
  geom_line(aes(y = risk1, colour = "1")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
  ylab("Cumulative Incidence") + 
  labs(colour="Concurrent Treatment:") +
  theme_bw() + 
  theme(legend.position="bottom")


#################################################################################
###                                                                           ###
###     DATA ANALYSIS EXERCISE 4                                              ###
###                                                                           ###
#################################################################################


################################
### a) Fit a non-saturated logistic model to estimate the denominator for 
###   nonstabilized IP weights use MSM with nonstabilized weights 

ipw.model <- glm(
  treatconcurrent ~ age + I(age*age) + time_since_study_start_days + 
    I(time_since_study_start_days*time_since_study_start_days) + season +
    imd_average_decile + ecog + smoking_history + tumour_side_left + 
    tumour_size + node_size + cardiac_comor + resp_comor + renal_comor + 
    other_comor  + time + I(time*time) + 
    I(treatconcurrent*time) + I(treatconcurrent*time*time), 
  family=binomial(), 
  data=imp_test.surv)

################################
### b) Use the predicted probabilities from this model to estimate the 
### nonstabilized IP weights for each individual

imp_test.surv$ipw.denom <- predict(ipw.model, imp_test.surv, type="response")

imp_test.surv$w_a <- ifelse(imp_test.surv$treatconcurrent==1, 1/imp_test.surv$ipw.denom, 1/(1-imp_test.surv$ipw.denom))


################################
### c) Fit an IP-weighted pooled logistic regression model using the 
### nonstabilized IP weights

# Converting data from wide to long form (note: both tidyverse and data.table 
# versions are provided below. Running either set of scripts will generate the
# same long format of the imp_test.surv dataset)


# Expand using tidyverse
imp_test.surv <- uncount(imp_test, weights=survtime, .remove=F) # Expand
imp_test.surv <- imp_test.surv %>% group_by(patient_id) %>% mutate(time=row_number()-1) %>% ungroup() # Time variable
imp_test.surv$timesq <- imp_test.surv$time^2 # Time squared variable
imp_test.surv$event <- ifelse( # Event variable
  imp_test.surv$time==imp_test.surv$survtime-1 & imp_test.surv$deceased==1, 1, 0)

# MSM with nonstabilized weights
options(warn=-1) # Need to suppress warning or else geeglm will encounter
# error due to non-integer number of successes as a result
# of weights

msm.w <- glm(
  event ~ treatconcurrent + time + timesq + treatconcurrent*time + treatconcurrent*timesq,
  family=binomial(), 
  weight=w_a,
  data=imp_test.surv)

summary(msm.w)

################################
### d) Use estimates from the pooled logistic model to estimate the risk 
###  at each month of follow-up for each treatment group

# Create a dataset to store results
# Include all time points under each treatment level
results0.w <- data.frame(treatconcurrent=0, time=seq(0,119), timesq=(seq(0,119))^2)
results1.w <- data.frame(treatconcurrent=1, time=seq(0,119), timesq=(seq(0,119))^2)

# Obtain predicted hazards from pooled logistic regression model
results0.w$hazard0 <- predict(msm.w, results0.w, type="response")
results1.w$hazard1 <- predict(msm.w, results1.w, type="response")

# Estimate survival probabilities from hazards
# S(t) = cumulative product of (1 - h(t))
results0.w$surv0 <- cumprod(1-results0.w$hazard0)
results1.w$surv1 <- cumprod(1-results1.w$hazard1)

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
results0.w$risk0 <- 1 - results0.w$surv0
results1.w$risk1 <- 1 - results1.w$surv1

# Risks
risk0.w <- results0.w[results0.w$time==119,]$risk0
risk1.w <- results1.w[results1.w$time==119,]$risk1

risk1.w
risk0.w

################################
### e) Estimating risk difference and risk ratio

risk1.w - risk0.w
risk1.w / risk0.w


################################
### f) Use bootstrapping to compute a valid 95% confidence interval for the risk ratio

survival.ipw.boot <- function(data, indices) {
  d <- data[indices,]
  
  # Fit model for the denominator
  boot.num.model <- glm(
    treatconcurrent ~ age + I(age*age) + time_since_study_start_days + 
      I(time_since_study_start_days*time_since_study_start_days) + season +
      imd_average_decile + ecog + smoking_history + tumour_side_left + 
      tumour_size + node_size + cardiac_comor + resp_comor + renal_comor + 
      other_comor + alcohol + substance_abuse + time + I(time*time) + 
      I(treatconcurrent*time) + I(treatconcurrent*time*time), 
    family=binomial(), 
    data=d)
  
  d$ipw.denom <- predict(boot.num.model, d, type="response")
  
  # Expanding dataset in person-time format
  boot.surv <- d[rep(seq(nrow(d)), times=d$survtime,)] # Expand
  boot.surv[, time:=seq(1:.N)-1, by=patient_id] # Time variable
  boot.surv$timesq <- boot.surv$time^2 # Time squared variable
  boot.surv$event <- ifelse( # Event variable
    boot.surv$time==boot.surv$survtime-1 & boot.surv$deceased==1, 1, 0)  
  
  boot.surv$w_a <- ifelse(boot.surv$treatconcurrent==1, 1/boot.surv$ipw.denom, 
                          1/(1-boot.surv$ipw.denom))
  
  
  # MSM with stabilized weights
  options(warn=-1) # Need to suppress warning or else geeglm will encounter
  # error due to non-integer number of successes as a result
  # of weights
  
  boot.msm.w <- glm(
    event ~ treatconcurrent + time + timesq + treatconcurrent*time * treatconcurrent*timesq,
    family=binomial(), weight=w_a,
    data=boot.surv)
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(treatconcurrent=0, time=seq(0,119), timesq=(seq(0,119))^2)
  results1 <- data.frame(treatconcurrent=1, time=seq(0,119), timesq=(seq(0,119))^2)
  
  # Obtain predicted hazards from pooled logistic regression model
  results0$hazard0 <- predict(boot.msm.w, results0, type="response")
  results1$hazard1 <- predict(boot.msm.w, results1, type="response")
  
  # Estimate survival probabilities from hazards
  # S(t) = cumulative product of (1 - h(t))
  results0$surv0 <- cumprod(1-results0$hazard0)
  results1$surv1 <- cumprod(1-results1$hazard1)
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  
  # Risks, risk differences and risk ratios
  risk0 <- results0[results0$time==119,]$risk0
  risk1 <- results1[results1$time==119,]$risk1
  
  rd <- risk1 - risk0
  rr <- risk1 / risk0
  
  return(c(risk1, risk0, rd, rr)) 
}

# 100 bootstrap samples for demonstration purposes 
set.seed(5435305)
survival.ipw.results <- boot(data=imp_test.surv, statistic=survival.ipw.boot, R=100)

boot.ci(survival.ipw.results,
        conf = 0.95,
        type = "perc",
        index = 4)


################################
### h) Fit a saturated logistic model for the numerator and use the predicted 
### probabilities from this model to estimate the stabilized IP weights 
### for each individual

treatconcurrent.model <- glm(
  treatconcurrent ~ 1, 
  family=binomial(), 
  data=imp_test)

imp_test.surv$ipw.num <- predict(treatconcurrent.model, imp_test.surv, type="response")

imp_test.surv$sw_a <- ifelse(imp_test.surv$treatconcurrent==1, imp_test.surv$ipw.num/imp_test.surv$ipw.denom, 
                             (1-imp_test.surv$ipw.num)/(1-imp_test.surv$ipw.denom))


################################
### i) Fit an IP weighted pooled logistic model using stabilized weights 

# MSM with stabilized weights
options(warn=-1) # Need to suppress warning or else geeglm will encounter
# error due to non-integer number of successes as a result
# of weights

msm.sw <- glm(
  event ~ treatconcurrent + time + timesq + treatconcurrent*time + treatconcurrent*timesq,
  family=binomial(), weight=sw_a,
  data=imp_test.surv)

summary(msm.sw)


################################
### j) Use estimates from the pooled logistic model to estimate the risk 
###  at each month of follow-up for each treatment group

# Create a dataset to store results
# Include all time points under each treatment level
results0.sw <- data.frame(treatconcurrent=0, time=seq(0,119), timesq=(seq(0,119))^2)
results1.sw <- data.frame(treatconcurrent=1, time=seq(0,119), timesq=(seq(0,119))^2)

# Obtain predicted hazards from pooled logistic regression model
results0.sw$hazard0 <- predict(msm.sw, results0.sw, type="response")
results1.sw$hazard1 <- predict(msm.sw, results1.sw, type="response")

# Estimate survival probabilities from hazards
# S(t) = cumulative product of (1 - h(t))
results0.sw$surv0 <- cumprod(1-results0.sw$hazard0)
results1.sw$surv1 <- cumprod(1-results1.sw$hazard1)

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
results0.sw$risk0 <- 1 - results0.sw$surv0
results1.sw$risk1 <- 1 - results1.sw$surv1

# Risks, risk differences and risk ratios
risk0.sw <- results0.sw[results0.sw$time==119,]$risk0
risk1.sw <- results1.sw[results1.sw$time==119,]$risk1

risk1.sw
risk0.sw


################################
### k) Estimating risk difference and risk ratio 

risk1.sw - risk0.sw
risk1.sw / risk0.sw


################################
### l) Use bootstrapping to compute a valid 95% confidence interval for the risk ratio

survival.ipw.sw.boot <- function(data, indices) {
  d <- data[indices,]
  
  # Fit model for the denominator
  ipw.model <- glm(
    treatconcurrent ~ age + I(age*age) + time_since_study_start_days + 
      I(time_since_study_start_days*time_since_study_start_days) + season +
      imd_average_decile + ecog + smoking_history + tumour_side_left + 
      tumour_size + node_size + cardiac_comor + resp_comor + renal_comor + 
      other_comor + alcohol + substance_abuse + time + I(time*time) + 
      I(treatconcurrent*time) + I(treatconcurrent*time*time), 
    family=binomial(), 
    data=d)
  
  d$ipw.denom <- predict(ipw.model, d, type="response")
  
  
  # Expanding dataset in person-time format
  boot.surv <- d[rep(seq(nrow(d)), times=d$survtime,)] # Expand
  boot.surv[, time:=seq(1:.N)-1, by=patient_id] # Time variable
  boot.surv$timesq <- boot.surv$time^2 # Time squared variable
  boot.surv$event <- ifelse( # Event variable
    boot.surv$time==boot.surv$survtime-1 & boot.surv$deceased==1, 1, 0)  
  
  # Fit model for the numerator 
  treatconcurrent.model <- glm(
    treatconcurrent ~ 1, 
    family=binomial(), 
    data=d)
  
  boot.surv$ipw.num <- predict(treatconcurrent.model, boot.surv, type="response")
  
  boot.surv$sw_a <- ifelse(boot.surv$treatconcurrent==1, boot.surv$ipw.num/boot.surv$ipw.denom, 
                           (1-boot.surv$ipw.num)/(1-boot.surv$ipw.denom))
  
  # MSM with stabilized weights
  options(warn=-1) # Need to suppress warning or else geeglm will encounter
  # error due to non-integer number of successes as a result
  # of weights
  
  boot.msm.sw <- glm(
    event ~ treatconcurrent + time + timesq + treatconcurrent*time * treatconcurrent*timesq,
    family=binomial(), weight=sw_a,
    data=boot.surv)
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0.sw <- data.frame(treatconcurrent=0, time=seq(0,119), timesq=(seq(0,119))^2)
  results1.sw <- data.frame(treatconcurrent=1, time=seq(0,119), timesq=(seq(0,119))^2)
  
  # Obtain predicted hazards from pooled logistic regression model
  results0.sw$hazard0 <- predict(boot.msm.sw, results0.sw, type="response")
  results1.sw$hazard1 <- predict(boot.msm.sw, results1.sw, type="response")
  
  # Estimate survival probabilities from hazards
  # S(t) = cumulative product of (1 - h(t))
  results0.sw$surv0 <- cumprod(1-results0.sw$hazard0)
  results1.sw$surv1 <- cumprod(1-results1.sw$hazard1)
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  results0.sw$risk0 <- 1 - results0.sw$surv0
  results1.sw$risk1 <- 1 - results1.sw$surv1
  
  # Risks, risk differences and risk ratios
  risk0.sw <- results0.sw[results0.sw$time==119,]$risk0
  risk1.sw <- results1.sw[results1.sw$time==119,]$risk1
  
  rd <- risk1.sw - risk0.sw
  rr <- risk1.sw / risk0.sw
  
  return(c(risk1.sw, risk0.sw, rd, rr)) 
}

# 100 bootstrap samples for demonstration purposes 
set.seed(5435305)
survival.ipw.sw.results <- boot(data=imp_test.surv, statistic=survival.ipw.sw.boot, R=100)

# Confidence interval for the risk ratio
boot.ci(survival.ipw.sw.results,
        conf = 0.95,
        type = "perc",
        index = 4)


################################
### o) Construct marginal cumulative incidence (risk) curves for all-cause mortality, by treatment group.

# Combine results for each treatment group into a single dataset
results.combined <- merge(results0.sw, results1.sw, by=c("time", "timesq"))

# Create a new "time" variable to reflect the fact that risks start at 0
# and are estimated at the end of each time interval
results.combined$time_updated <- results.combined$time + 1
results.combined <- results.combined %>% add_row(time_updated=0, risk0=0, risk1=0) %>% arrange(time_updated)  

# plot
ggplot(results.combined, aes(x=time_updated, y=surv)) + 
  geom_line(aes(y = risk0, colour = "0")) + 
  geom_line(aes(y = risk1, colour = "1")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
  ylab("Cumulative Incidence") + 
  labs(colour="Concurrent Treatment:") +
  theme_bw() + 
  theme(legend.position="bottom")




#################################################################################
###                                                                           ###
###     OPTIONAL: DATA ANALYSIS EXERCISE 5                                    ###
###                                                                           ###
#################################################################################


################################
### a) Using the nearest neighbor approach, create dataset in which all treated 
###   are matched to untreated individuals.

set.seed(1234)

# Run the matching algorithm based on nearest neighbor approach
nearest.match <- matchit(treatconcurrent ~ age + I(age*age) + time_since_study_start_days + 
                           I(time_since_study_start_days*time_since_study_start_days) + season +
                           imd_average_decile + ecog + smoking_history + tumour_side_left + 
                           tumour_size + node_size + cardiac_comor + resp_comor + renal_comor + 
                           other_comor + alcohol + substance_abuse + time + I(time*time) + 
                           I(treatconcurrent*time) + I(treatconcurrent*time*time), 
                         method = "nearest", distance = "glm", link="logit", discard = "treat",
                         caliper=0.2,
                         data = imp_test)


# Extract the data 
imp_test.matched <- match.data(nearest.match)

################################
### b) Modify the matched sample from wide to person-time format

imp_test.matched$survtime <- ifelse(
  imp_test.matched$deceased==0, 120,
  (imp_test.matched$yrdth-83)*12+imp_test.matched$modth) # yrdth ranges from 83 to 92

# Expand using tidyverse
imp_test.matched.surv <- uncount(imp_test.matched, weights=survtime, .remove=F) # Expand
imp_test.matched.surv <- imp_test.matched.surv %>% group_by(patient_id) %>% mutate(time=row_number()-1) %>% ungroup() # Time variable
imp_test.matched.surv$timesq <- imp_test.matched.surv$time^2 # Time squared variable
imp_test.matched.surv$event <- ifelse( # Event variable
  imp_test.matched.surv$time==imp_test.matched.surv$survtime-1 & imp_test.matched.surv$deceased==1, 1, 0)


################################
### c) Fitting a pooled logistic regression model with time (linear and 
###    quadratic terms), treatment and product terms between treatment and time  

fit.pool.match <- glm(
  formula = event==1 ~ treatconcurrent + time + timesq + treatconcurrent*time + treatconcurrent*timesq,
  family = binomial(link = 'logit'),
  data = imp_test.matched.surv)

summary(fit.pool.match)

################################
### d) Use estimates from the pooled logistic model to estimate the risk 
###  at each month of follow-up for each treatment group

# Create datasets to store results
# Include all time points under each treatment level
results0_match <- data.frame(treatconcurrent=0, time=seq(0,119), timesq=(seq(0,119))^2)
results1_match <- data.frame(treatconcurrent=1, time=seq(0,119), timesq=(seq(0,119))^2)

# Obtain predicted hazards from pooled logistic regression model
results0_match$hazard0 <- predict(fit.pool.match, results0_match, type="response")
results1_match$hazard1 <- predict(fit.pool.match, results1_match, type="response")

# Estimate survival probabilities from hazards
# S(t) = cumulative product of (1 - h(t))
results0_match$surv0 <- cumprod(1-results0_match$hazard0)
results1_match$surv1 <- cumprod(1-results1_match$hazard1)

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
results0_match$risk0 <- 1 - results0_match$surv0
results1_match$risk1 <- 1 - results1_match$surv1


risk1_match <- results1_match[results1_match$time==119,]$risk1
risk0_match <- results0_match[results0_match$time==119,]$risk0

risk1_match
risk0_match

################################
### e) Use the risks to estimate the effect on the risk ratio and risk 
### difference scales

risk1_match - risk0_match
risk1_match / risk0_match

################################
### g) Constructing risk curves

# Combine results for each treatment group into a single dataset
results.match.combined <- merge(results0_match, results1_match, by=c("time", "timesq"))

# Create a new "time" variable to reflect the fact that risks start at 0
# and are estimated at the end of each time interval
results.match.combined$time_updated <- results.match.combined$time + 1
results.match.combined <- results.match.combined %>% add_row(time_updated=0, risk0=0, risk1=0) %>% arrange(time_updated)

# Creating plot
ggplot(results.match.combined, aes(x=time_updated)) + 
  geom_line(aes(y = risk0, colour = "0")) + 
  geom_line(aes(y = risk1, colour = "1")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 120), breaks=seq(0,120,12)) +
  ylab("Cumulative Incidence") + 
  labs(colour="Concurrent Treatment:") +
  theme_bw() + 
  theme(legend.position="bottom")


