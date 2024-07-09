# Data preparation
# 
# The purpose of this script is to load the group summary raw data, process it,
# and save it as Rds file in the processed-data folder.
#
# Created by Isaac Chen 2024.05.06

# Clear the global environment
rm(list = ls())

## ---- load_packages ----
library(here)
library(dplyr)
library(stringr)
library(skimr)

## ---- load_rawdata ----
data_path <- here("data", "raw_data", "Group_summary_table_2024.04.10.csv")
raw_data <- read.csv(data_path)

## ---- recode_variables ----
# Trim white spaces to avoid multiple same categories
# Re-categorize etiology as TBI vs Non-TBI
# Re-categorize race as White vs Black vs Other
processed_data <- raw_data |>
  mutate(
    across(where(is.character), str_trim), 
    Transfer_Status = ifelse(
      Transfer_Status == "Unplanned", "Unplanned", "None/Planned"
    ),
    Screening_Status_24hr = ifelse(
      Delirium_24hours == 1, "Ever Positive", "Never Positive"
    ),
    Etiology_TBIOther = if_else(Etiology == "TBI", "TBI", "Non-TBI")
  ) 

## ---- select_variables ----
processed_data <- processed_data |> 
  select(
    Subject_ID, Excluded, Screening_Status_24hr, Delirium_24hours,
    Transfer_Status, Etiology_TBIOther, Sex, Race, Compliance_24hours, 
    AgeAtFirstAdmission, AcuteLOS, RehabLOS, CAPD_24hrScore,
    WeeFIM_CogDFQ_admission, WeeFIM_CogDFQ_discharge, 
  )

## ---- rename_variables ----
# Transform quantitative variables (all columns after age) to numeric data type
# Rename variables to comply with tidyverse style guide and clarity
processed_data <- processed_data |>
  mutate(across(AgeAtFirstAdmission:WeeFIM_CogDFQ_discharge, as.numeric)) |>
  rename(
    Acute_LOS = AcuteLOS,
    Rehab_LOS = RehabLOS,
    Screening_Status_24hr_code = Delirium_24hours,
    Compliance_24hr = Compliance_24hours,
    Age_Adm = AgeAtFirstAdmission,
    WeeFIM_CogDFQ_Adm= WeeFIM_CogDFQ_admission,
    WeeFIM_CogDFQ_ultDC = WeeFIM_CogDFQ_discharge,
  )

## ---- code_factors ----
processed_data <- processed_data |> 
  mutate(
    Screening_Status_24hr = factor(
      Screening_Status_24hr,levels = c("Never Positive", "Ever Positive")
    ),
    Transfer_Status = factor(
      Transfer_Status, levels = c("None/Planned", "Unplanned")
    ),
    Etiology_TBIOther = factor(Etiology_TBIOther, levels = c("Non-TBI", "TBI")),
    Sex = factor(Sex, levels = c("Male", "Female")),
    Race = factor(
      Race, levels = c("CA", "AA", "LAT", "AS" ,"OT")
    ),
  )

## ---- data_summary ----
# Check if everything looks good
skim(processed_data)

## ---- save_data --------
# location to save file
processed_data_path <- here("data","processed_data","CAPD_prognosis.rds")
saveRDS(processed_data, file = processed_data_path)


