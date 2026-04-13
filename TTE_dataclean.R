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
patient <- read.csv("E:/Theragnostics/Queries/Charlie/PROTECT_validation_R/data/PROTECT_patient5.csv")




