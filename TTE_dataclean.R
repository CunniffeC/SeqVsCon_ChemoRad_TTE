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
dat_comor <- dat_DS %>%
  select(patient_id, acecomorbidities)%>%
  drop_na(acecomorbidities)%>%
  distinct() %>%
  mutate(overall_comor = ifelse(grepl("Overall", acecomorbidities, fixed = T), parse_number(acecomorbidities), NA)) %>%
  mutate(alcohol = ifelse(grepl("Alcohol", acecomorbidities, fixed = T), parse_number(gsub(".*Alcohol", "", acecomorbidities)),NA))%>%
  mutate(respiratory_comor = ifelse(grepl("Respiratory", acecomorbidities, fixed = T), parse_number(gsub(".*Respiratory", "", acecomorbidities)),NA))%>%
  mutate(substance_abuse = ifelse(grepl("Drugs", acecomorbidities, fixed = T), parse_number(gsub(".*Drugs", "", acecomorbidities)),NA))%>%
  mutate(renal_comor = ifelse(grepl("renal", acecomorbidities, fixed = T), parse_number(gsub(".*renal", "", acecomorbidities)),NA))#%>%
#make season

#check mutationrates 
