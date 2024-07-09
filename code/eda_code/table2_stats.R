# Table 2 Patient Outcomes
#
# This script aims to produce all the values that appear in Table 1.
# 
# Created by Isaac Chen 2024.05.20

# Clear all the variables in the environment
rm(list = ls())

## ---- load_packages ----
lib_list = c("here", "dplyr", "forcats")
lapply(lib_list, library, character.only = TRUE)

## ---- load_rawdata ----
data_path <- here("data", "processed_data", "CAPD_prognosis.rds")
CAPD_prog <- readRDS(data_path)

## ---- filter_excluded ----
# Filter out patients who are excluded from analysis
CAPD_prog <- CAPD_prog |>
  filter(Excluded == "")

## ---- binary_predictor ----
# Create binary predictor variable and encode it as a factor
CAPD_prog <- CAPD_prog |> 
  mutate(
    Age_Adm_12cutoff = ifelse(Age_Adm >= 12, ">=12 yrs", "<12 yrs"), 
    Acute_LOS_30cutoff = ifelse(Acute_LOS > 30, ">30 Days", "<=30 Days"),
    WeeFIM_CogDFQ_70cutoff = ifelse(
      WeeFIM_CogDFQ_ultDC >= 70, "Low Assistance", "Medium-High Assistance"
    )
  ) |>
  mutate(
    Age_Adm_12cutoff = factor(
      Age_Adm_12cutoff, 
      levels = c("<12 yrs", ">=12 yrs")
    ),
    Acute_LOS_30cutoff = factor(
      Acute_LOS_30cutoff, 
      levels = c("<=30 Days", ">30 Days")
    ),
    WeeFIM_CogDFQ_70cutoff = factor(
      WeeFIM_CogDFQ_70cutoff, 
      levels = c("Low Assistance", "Medium-High Assistance")
    )
  )

# Table 2 stats
# Entire cohort
table(CAPD_prog$Transfer_Status, CAPD_prog$Screening_Status_24hr)
chisq.test(table(CAPD_prog$Transfer_Status, CAPD_prog$Screening_Status_24hr), correct = F)

# Only children with WeeFIM at discharge
CAPD_prog <- CAPD_prog |> 
  filter(!is.na(WeeFIM_CogDFQ_ultDC))

table(CAPD_prog$WeeFIM_CogDFQ_70cutoff, CAPD_prog$Screening_Status_24hr)
fisher.test(table(CAPD_prog$WeeFIM_CogDFQ_70cutoff, CAPD_prog$Screening_Status_24hr))

tapply(CAPD_prog$WeeFIM_CogDFQ_ultDC, CAPD_prog$Screening_Status_24hr, length)
tapply(CAPD_prog$WeeFIM_CogDFQ_ultDC, CAPD_prog$Screening_Status_24hr, 
       FUN = function(x) quantile(x, 0.25))
tapply(CAPD_prog$WeeFIM_CogDFQ_ultDC, CAPD_prog$Screening_Status_24hr, 
       FUN = function(x) quantile(x, 0.50))
tapply(CAPD_prog$WeeFIM_CogDFQ_ultDC, CAPD_prog$Screening_Status_24hr, 
       FUN = function(x) quantile(x, 0.75))

wilcox.test(CAPD_prog$WeeFIM_CogDFQ_ultDC ~ CAPD_prog$Screening_Status_24hr,
            correct = F)





