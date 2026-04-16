#### TABLE ONE AND POSITIVITY CHECKS TTE ANALYSIS ####
## By Charlie Cunniffe
## 15/04/26



###########################
##### RUN SECOND ##########
##########################




###### Libraries ######
#library(plyr)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(data.table)
library(DescTools)
library(mice)
library(anytime)
library(irr)
library(tableone)
library(gtsummary)
library(gt)
library(flextable)
library(officer)
library(survival)
library(paletteer)
my_colors <- paletteer::paletteer_d("lisa::FridaKahlo")


#### data ####
#TTE_final <- read.csv("E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/manchester4_update.csv")
#write.csv(TTE_final$pid, "E:/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/final_patient_id_list.csv", row.names = FALSE)
TTE_final <- TTE_final %>%
  mutate_if(is.integer, as.numeric)%>%
  mutate_if(is.logical, as.numeric)%>%
  mutate(y= ifelse(deceased == 0, - time_days, time_days))%>%
  #rename(tx = treatconcurrent,
  #       Patient_ID = X)%>%
  mutate(Treatment = as.factor(ifelse(treatconcurrent == 1, "Concurrent", "Sequential")))


####summary table ####
tab_data <- TTE_final %>%
  select(-patient_id, -date_seen) %>%
  mutate(agerange = ifelse(is.na(age), "Unknown",
                           ifelse(age < 60, "< 60",
                                  ifelse(age < 65, "60 - 64",
                                         ifelse(age < 70, "65 - 70", ">=70"))))) %>%
  mutate(time_5yr = ifelse( time_days >= 5*365, "Yes", ifelse(deceased == 0, "Censored", "No"))) %>%
  select(-deceased, -treatconcurrent, -y)
tab_sum <- list(
  time_5yr ~ "5 year survival",
  time_days ~ "Survival (days)",
  sex_male ~ "Male Sex",
  age ~ "Median age (years)",
  agerange ~ "Age range (years)",
  imd_average_decile ~ "Index of multiple deprivation (decile)",
  ecog ~ "Performance status (ECOG)",
  #Treatment ~ "Treatment",
  smoking_history ~ "Smoking History",
  tumour_side_left ~ "Left sided tumour",
  tumour_size ~ "Tumour volume (T stage)",
  node_size ~ "Nodal volume (N stage)",
  resp_comor ~ "Pulmonary comorbidity (ACE)",
  cardiac_comor ~ "Cardiac comorbidity (ACE)",
  renal_comor ~ "Renal comorbidity (ACE)",
  other_comor ~ "Other comorbidity (ACE)",
  alcohol ~ "Alcohol abuse (ACE)",
  substance_abuse ~ "Illicit drug abuse (ACE)",
  season ~ "Month treatment began",
  time_since_study_start_days ~ "Days since study start"
)
#get table one
tbl.gts.one <- tab_data %>%
  tbl_summary(by = Treatment, label = tab_sum)
tbl.gts.one


#### full cohort ####
#save
tbl.gts.one.tib <- as_tibble(tbl.gts.one)
write.csv(tbl.gts.one.tib, "/mnt/data/charliec/ChemoRad_TTE/table1_full.csv", row.names = FALSE)
tabone <- tbl.gts.one %>%
  as_flex_table()%>%
  autofit()%>%
  flextable::save_as_docx(path="/mnt/data/charliec/ChemoRad_TTE/table1_full.docx")


age_hist <- ggplot(data = TTE_final, aes(x=age, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Age (years)")+
  labs(title = "Age distribution in Manchester by received treatment \n")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
age_hist

p_score_mod <- glm(Treatment ~ stageIIIBC + age + histo, data = tab_data, family = binomial)
tab_data$p_score <- predict(p_score_mod, newdata = tab_data, type = "response")
p_score_hist <- ggplot(data = tab_data, aes(x=p_score, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Propensity score")+
  labs(title = "Propensity score distribution in Manchester by received treatment \n")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
p_score_hist


#### RCT cohort ####
sel_tab_data <- tab_data %>%
  select(-p_score) %>%
  filter( age <= 75 & ecog <2)
#get table Treatment
sel.tbl.gts.Treatment <- sel_tab_data %>%
  tbl_summary(by = Treatment, label = tab_sum, statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n}"))
sel.tbl.gts.Treatment
#save
sel.tbl.gts.Treatment.tib <- as_tibble(sel.tbl.gts.Treatment)
write.csv(sel.tbl.gts.Treatment.tib, "E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/tab_Treatment_age75ecog01.csv", row.names = FALSE)
sel.tabTreatment <- sel.tbl.gts.Treatment %>%
  as_flex_table()%>%
  autofit()%>%
  flextable::save_as_docx(path="E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/tab_Treatment_age75ecog01.docx")
age_hist_age75ecog01 <- ggplot(data = sel_tab_data, aes(x=age, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Age (years)")+
  labs(title = "Age distribution in Manchester by received treatment, \nRCT adjacent (age<75 & ECOG<2)")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
age_hist_age75ecog01
p_score_mod <- glm(Treatment ~ stageIIIBC + age + histo, data = sel_tab_data, family = binomial)
sel_tab_data$p_score <- predict(p_score_mod, newdata = sel_tab_data, type = "response")
p_score_hist_age75ecog01 <- ggplot(data = sel_tab_data, aes(x=p_score, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Propensity score")+
  labs(title = "Propensity score distribution in Manchester by received treatment, \nRCT adjacent (age<75 & ECOG<2)")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
p_score_hist_age75ecog01
#### overlap cohort ####
sel2_tab_data <- tab_data %>%
  select(-p_score) %>%
  filter( age <= 80)
#get table Treatment
sel2.tbl.gts.Treatment <- sel2_tab_data %>%
  tbl_summary(by = Treatment, label = tab_sum, statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n}"))
sel2.tbl.gts.Treatment
#save
sel2.tbl.gts.Treatment.tib <- as_tibble(sel2.tbl.gts.Treatment)
write.csv(sel2.tbl.gts.Treatment.tib, "E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/tab_Treatment_age80.csv", row.names = FALSE)
sel2.tabTreatment <- sel2.tbl.gts.Treatment %>%
  as_flex_table()%>%
  autofit()%>%
  flextable::save_as_docx(path="E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/tab_Treatment_age80.docx")
age_hist_age80 <- ggplot(data = sel2_tab_data, aes(x=age, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Age (years)")+
  labs(title = "Age distribution in Manchester by received treatment, \ncomplete overlap (age<80)")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
age_hist_age80
p_score_mod <- glm(Treatment ~ stageIIIBC + age + histo, data = sel2_tab_data, family = binomial)
sel2_tab_data$p_score <- predict(p_score_mod, newdata = sel2_tab_data, type = "response")
p_score_hist_age80 <- ggplot(data = sel2_tab_data, aes(x=p_score, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Propensity score")+
  labs(title = "Propensity score distribution in Manchester by received treatment, \ncomplete overlap (age<80)")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
p_score_hist_age80
#### RCT overlap cohort ####
sel3_tab_data <- tab_data %>%
  select(-p_score) %>%
  filter( age <= 75 & ecog <2 & age>= 45)
#get table Treatment
sel3.tbl.gts.Treatment <- sel3_tab_data %>%
  tbl_summary(by = Treatment, label = tab_sum, statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n}"))
sel3.tbl.gts.Treatment
#save
sel3.tbl.gts.Treatment.tib <- as_tibble(sel3.tbl.gts.Treatment)
write.csv(sel3.tbl.gts.Treatment.tib, "E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/tab_Treatment_age45-75ecog01.csv", row.names = FALSE)
sel3.tabTreatment <- sel3.tbl.gts.Treatment %>%
  as_flex_table()%>%
  autofit()%>%
  flextable::save_as_docx(path="E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/protectcode/peppercorn_Rcodes_280425/data/tab_Treatment_age45-75ecog01.docx")
age_hist_RCT <- ggplot(data = sel3_tab_data, aes(x=age, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Age (years)")+
  labs(title = "Age distribution in Manchester by received treatment,\nRCT adjacent complete overlap (45<age<80 & ECOG<2 )")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
age_hist_RCT
p_score_mod <- glm(Treatment ~ stageIIIBC + age + histo, data = sel3_tab_data, family = binomial)
sel3_tab_data$p_score <- predict(p_score_mod, newdata = sel3_tab_data, type = "response")
p_score_hist_RCT <- ggplot(data = sel3_tab_data, aes(x=p_score, fill = Treatment, colour = Treatment)) + 
  geom_histogram(position = 'dodge', alpha = 0.5) +
  xlab("Propensity score")+
  labs(title = "Propensity score distribution in Manchester by received treatment,\nRCT adjacent complete overlap (45<age<80 & ECOG<2 )")+
  scale_fill_paletteer_d("khroma::highcontrast") +
  scale_colour_paletteer_d("khroma::highcontrast")+
  theme_light()+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "right")
p_score_hist_RCT