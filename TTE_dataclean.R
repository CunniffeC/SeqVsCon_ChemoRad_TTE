##### libraries #####
library(tidyverse)
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

##### functions #####
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
my.first <- function(x) ifelse( !all(is.na(x)), first(x, na.rm=T), NA)
my.closest <- function(x, y) which.min(abs(x-y))




final_date <-  as.Date("2023-07-01", format ="%Y-%m-%d" ) # as.numeric?

set.seed(2602)


##### import data #####
dat_protect <- read.csv("/mnt/data/charliec/ChemoRad_TTE/Data/manchester4_update.csv")%>%
  rename(patient_id = pid)%>%
  rename(date_seen = fup_start)%>%
  mutate(date_seen = anydate(date_seen))
id_dates <- dat_protect %>% select(patient_id, date_seen)
dat_DS1 <- read.csv("/mnt/data/charliec/ChemoRad_TTE/Data/TTE_DS_data.csv") %>% mutate(date_seen = anydate(date_seen))
dat_DS <- left_join(id_dates, dat_DS1, by = c("patient_id", "date_seen"))
dat_patient <- read.csv("/mnt/data/charliec/ChemoRad_TTE/Data/TTE_patient_data.csv") %>% select(patient_id, imd_average_decile)

TTE_final <- dat_protect %>%
  select(patient_id, date_seen, ecog, deceased, time,  age, treatconcurrent) #weight,

TTE_final <- left_join(TTE_final, dat_patient, by="patient_id")

#smoking status, side 
clean_DS <- dat_DS %>% 
  rename(tumour_side  = side)%>%
  select(patient_id, smoking_history, tumour_side) %>%
  distinct()
TTE_final <- left_join(TTE_final, clean_DS, by = "patient_id")

#extract T and N stage
dat_stage <- dat_DS %>%
  select(patient_id, clinical_stage, path_stage) %>%
  distinct() %>%
  mutate(clinical_stage = ifelse(is.na(clinical_stage), path_stage, clinical_stage)) %>%
  select(-path_stage)

dat_stage$stage <- gsub("\\,.*","", dat_stage$clinical_stage)
dat_stage$stage <- gsub("\\M.*","", dat_stage$stage)
dat_stage$stage <- gsub("p", "", dat_stage$stage)
dat_stage$stage <- gsub(" ", "", dat_stage$stage)
dat_stage$stage <- gsub("T", "", dat_stage$stage)
dat_stage <- dat_stage %>%
  mutate(tumour_size = gsub("\\N.*","",stage))%>%
  mutate(node_size = gsub(".*N", "",stage)) %>%
  mutate(stage = ifelse(stage == "", NA, ifelse(stage == " ", NA, stage)))%>%
  drop_na(stage)

clean_stage <- dat_stage %>%
  select(patient_id, tumour_size, node_size)%>%
  distinct()

TTE_final <- left_join(TTE_final, clean_stage, by = "patient_id", relationship = "many-to-many") %>%
  distinct()

#break up comorbidities
cardiac_terms <- c(
  "Myocardial Infarct",
  "Angina / Coronary Artery Disease",
  "Congestive Heart Failure (CHF)",
  "Arrhythmias",
  "Hypertension",
  "Venous Disease",
  "Peripheral Arterial Disease"
)
gastro_terms <- c(
  "Hepatic",
  "Stomach / Intestine",
  "Pancreas"
)
neuro_terms <- c(
  "Stroke",
  "Dementia",
  "Paralysis",
  "neuromuscular"
)
malig_terms <- c(
  "Solid Tumor including melanoma",
  "Leukemia and Myeloma",
  "Lymphoma"
)


dat_comor <- dat_DS %>%
  select(patient_id, acecomorbidities)%>%
  drop_na(acecomorbidities)%>%
  distinct() %>%
  mutate(overall_comor = ifelse(grepl("Overall", acecomorbidities, fixed = T), parse_number(acecomorbidities), NA)) %>%
  mutate(alcohol = ifelse(grepl("Alcohol", acecomorbidities, fixed = T), parse_number(gsub(".*Alcohol", "", acecomorbidities)),NA))%>%
  mutate(resp_comor = ifelse(grepl("Respiratory", acecomorbidities, fixed = T), parse_number(gsub(".*Respiratory", "", acecomorbidities)),NA))%>%
  mutate(substance_abuse = ifelse(grepl("Drugs", acecomorbidities, fixed = T), parse_number(gsub(".*Drugs", "", acecomorbidities)),NA))%>%
  mutate(renal_comor = ifelse(grepl("renal", acecomorbidities, fixed = T), parse_number(gsub(".*renal", "", acecomorbidities)),NA))%>%
  mutate(endocrine_comor = ifelse(grepl("Diabetes", acecomorbidities, fixed = T), parse_number(gsub(".*Diabetes", "", acecomorbidities)),NA))%>%
  mutate(psych_comor = ifelse(grepl("Psychiatric", acecomorbidities, fixed = T), parse_number(gsub(".*Psychiatric", "", acecomorbidities)),NA))%>%
  mutate(rheum_comor = ifelse(grepl("Rheumatologic", acecomorbidities, fixed = T), parse_number(gsub(".*Rheumatologic", "", acecomorbidities)),NA))%>%
  mutate(AIDS_comor = ifelse(grepl("AIDS", acecomorbidities, fixed = T), parse_number(gsub(".*AIDS", "", acecomorbidities)),NA))%>%
  mutate(obese_comor = ifelse(grepl("Obesity", acecomorbidities, fixed = T), parse_number(gsub(".*Obesity", "", acecomorbidities)),NA))%>%
  #cardiac comor
  rowwise() %>%
  mutate(
    cardiac_grades = list(
      unlist(lapply(cardiac_terms, function(term) {
        if (str_detect(acecomorbidities, fixed(term))) {
          str_extract(acecomorbidities, paste0(term, " grade \\d")) %>%
            parse_number()
        } else {
          NA
        }
      }))
    ),
    cardiac_grades = list(na.omit(cardiac_grades)),
    cardiac_comor = case_when(
      length(cardiac_grades) == 0 ~ NA_real_,
      sum(cardiac_grades == 2) >= 2 ~ 3,
      TRUE ~ max(cardiac_grades)
    )
  ) %>%
  ungroup()%>%
  #gastro comor
  rowwise() %>%
  mutate(
    gastro_grades = list(
      unlist(lapply(gastro_terms, function(term) {
        if (str_detect(acecomorbidities, fixed(term))) {
          str_extract(acecomorbidities, paste0(term, " grade \\d")) %>%
            parse_number()
        } else {
          NA
        }
      }))
    ),
    gastro_grades = list(na.omit(gastro_grades)),
    gastro_comor = case_when(
      length(gastro_grades) == 0 ~ NA_real_,
      sum(gastro_grades == 2) >= 2 ~ 3,
      TRUE ~ max(gastro_grades)
    )
  ) %>%
  ungroup()%>%
  #neuro comor
  rowwise() %>%
  mutate(
    neuro_grades = list(
      unlist(lapply(neuro_terms, function(term) {
        if (str_detect(acecomorbidities, fixed(term))) {
          str_extract(acecomorbidities, paste0(term, " grade \\d")) %>%
            parse_number()
        } else {
          NA
        }
      }))
    ),
    neuro_grades = list(na.omit(neuro_grades)),
    neuro_comor = case_when(
      length(neuro_grades) == 0 ~ NA_real_,
      sum(neuro_grades == 2) >= 2 ~ 3,
      TRUE ~ max(neuro_grades)
    )
  ) %>%
  ungroup()%>%
  #malig comor
  rowwise() %>%
  mutate(
    malig_grades = list(
      unlist(lapply(malig_terms, function(term) {
        if (str_detect(acecomorbidities, fixed(term))) {
          str_extract(acecomorbidities, paste0(term, " grade \\d")) %>%
            parse_number()
        } else {
          NA
        }
      }))
    ),
    malig_grades = list(na.omit(malig_grades)),
    malig_comor = case_when(
      length(malig_grades) == 0 ~ NA_real_,
      sum(malig_grades == 2) >= 2 ~ 3,
      TRUE ~ max(malig_grades)
    )
  ) %>%
  ungroup()%>%
  rowwise() %>%
  mutate(
    n_twos = sum(c_across(
      ends_with("_comor") & 
        !all_of(c("cardiac_comor", "resp_comor", "renal_comor", "alcohol", "substance_abuse", "overall_comor"))
    ) == 2, na.rm = TRUE),
    
    max_val = max(c_across(
      ends_with("_comor") & 
        !all_of(c("cardiac_comor", "resp_comor", "renal_comor", "alcohol", "substance_abuse", "overall_comor"))
    ), na.rm = TRUE),
    
    other_comor = case_when(
      is.infinite(max_val) ~ NA_real_,
      n_twos >= 2 ~ 3,
      TRUE ~ max_val
    )
  ) %>%
  ungroup() %>%
  select(-n_twos, -max_val)



#make season

#check mutationrates 
