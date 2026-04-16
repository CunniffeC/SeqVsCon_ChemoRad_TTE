#### MULTIPLE IMPUTATION FOR TTE ANALYSIS####
## By Charlie Cunniffe
## 14/04/26


###########################
##### RUN THIRD ##########
##########################




##### LIBRARIES ####
library(mice)
library(survival)
library(rio) ## for xlsx files
library(howManyImputations)
library(ggmice)
set.seed(2602)

#TTE_final <- read.csv("data_location/data.csv")

##### CHECK MISSINGNESS PATTERNS ####
summary(TTE_final) ## provides expected ranges and number missing for each factor

plot_pattern(TTE_final, rotate = TRUE)

##### REMOVE TOO MISSING VARIABLES #####
TTE_final <- TTE_final %>% select(-kras, - braf, -met, -fgfr3, -pi3_k)
plot_pattern(TTE_final, rotate = TRUE)

TTE_final_mutations <- TTE_final
TTE_final <- TTE_final %>% select(-egfr, -pdl1_score, -alk)
plot_pattern(TTE_final, rotate = TRUE)

##### IMPUTATION #####
m=30 ## number of imputed datasets to make (default is 5 but more is recommended)
maxit = 30 ## number of iterations (default is 5 but more is recommended)

## make an empty mids object 
imp0 <- mice(TTE_final, maxit = 0)

## get the predictor matrix of the correct shape for your data
pred <- imp0$predictorMatrix 
## rows are being predicted by the columns containing 1's
## default pred is all 1's except for the diagonal
pred[c("patient_id","date_seen"),] <- 0
pred[,c("patient_id", "date_seen")] <- 0## should at least set both the column and row for patient ID to 0's
pred ## print to check its right

## get the methods matrix
meth <- imp0$meth
## specify the imputation method for each variable 
## defaults:
## pmm, predictive mean matching (numeric data) 
## logreg, logistic regression imputation (binary data, factor with 2 levels) 
## polyreg, polytomous regression imputation for unordered categorical data (factor > 2 levels) 
## polr, proportional odds model for (ordered, > 2 levels))

# COMPLETE IMPUTATION
imp2 <- mice(TTE_final, m=m, maxit = maxit, pred=pred, meth=meth)

## CHECK IMPUTED VALUES
plot(x=imp2, y=c(), layout = c(2,5))
## plots the mean or sd across the cohort on y, line for each imputed data set, x is the itterations
## want it to look stable, across all (or at least the last group of) iterations
imp2.dat <- complete(imp2, "long", include = TRUE)

imp2.dat <- imp2.dat %>% 
  mutate(imputed = .imp > 0,
         imputed = factor(imputed,
                          levels = c(F,T),
                          labels = c("Observed", "Imputed")))

smoking_history_imp <- prop.table(table(imp2.dat$smoking_history,
                 imp2.dat$imputed),
           margin = 2)
smoking_history_imp

imd_average_decile_imp <- prop.table(table(imp2.dat$imd_average_decile,
                                        imp2.dat$imputed),
                                  margin = 2)
imd_average_decile_imp
ecog_imp <- prop.table(table(imp2.dat$ecog,
                                        imp2.dat$imputed),
                                  margin = 2)
ecog_imp
tumour_side_left_imp <- prop.table(table(imp2.dat$tumour_side_left,
                                        imp2.dat$imputed),
                                  margin = 2)
tumour_side_left_imp
tumour_size_imp <- prop.table(table(imp2.dat$tumour_size,
                 imp2.dat$imputed),
           margin = 2)
tumour_size_imp 
node_size_imp <- prop.table(table(imp2.dat$node_size,
                 imp2.dat$imputed),
           margin = 2)
node_size_imp 
cardiac_comor_imp <- prop.table(table(imp2.dat$cardiac_comor,
                 imp2.dat$imputed),
           margin = 2)
cardiac_comor_imp 
resp_comor_imp <- prop.table(table(imp2.dat$resp_comor,
                                        imp2.dat$imputed),
                                  margin = 2)
resp_comor_imp
renal_comor_imp <- prop.table(table(imp2.dat$renal_comor,
                                        imp2.dat$imputed),
                                  margin = 2)
renal_comor_imp 
other_comor_imp <- prop.table(table(imp2.dat$other_comor,
                                        imp2.dat$imputed),
                                  margin = 2)
other_comor_imp 
alcohol_imp <- prop.table(table(imp2.dat$alcohol,
                                        imp2.dat$imputed),
                                  margin = 2)
alcohol_imp 
substance_abuse_imp <- prop.table(table(imp2.dat$substance_abuse,
                                        imp2.dat$imputed),
                                  margin = 2)
substance_abuse_imp 

prelim_fit <- with(imp2, coxph(Surv(time_days, deceased)~ 
                                 treatconcurrent + 
                                 age + 
                                 time_since_study_start_days +
                                 # season +
                                 imd_average_decile +
                                 ecog +
                                 smoking_history +
                                 tumour_side_left +
                                 tumour_size +
                                 node_size +
                                 cardiac_comor +
                                 resp_comor +
                                 # renal_comor +
                                 other_comor 
                                 ))
results.prelim_fit <- summary(pool(prelim_fit), conf.int=T)
results.prelim_fit

