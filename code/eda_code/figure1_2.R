# Plot of Continuous Predictors by Screening Status to determine Cutoffs
# 
# This script aims to plot all the continuous predictors by screening status
# as a visual means to decide what is the best cutoff to use for the logistic
# regression analysis.
# 
# Created by Isaac Chen 2024.05.21

# Clear all the variables in the environment
rm(list = ls())

## ---- load_packages ----
lib_list = c("here", "dplyr", "forcats", "ggplot2", "ggdist", "gghalves")
lapply(lib_list, library, character.only = TRUE)

## ---- make_dir ----
fig_dir <- here("results", "figures", "continuous_predictors_cutoff")
dir.create(fig_dir)

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
  ) |>
  mutate(
    Age_Adm_12cutoff = factor(
      Age_Adm_12cutoff, levels = c("<12 yrs", ">=12 yrs")
    ),
    Acute_LOS_30cutoff = factor(
      Acute_LOS_30cutoff, levels = c("<=30 Days", ">30 Days")
    )
  ) 

## ---- age_screening_plot ----
# Age at Admission by Screening Status
fig_age_screening <- CAPD_prog |> 
  ggplot(aes(x = Screening_Status_24hr, y = Age_Adm)) +
  stat_halfeye(
    aes(fill = Screening_Status_24hr),
    alpha = 0.5,
    adjust = 0.5,
    width = 0.5,
    .width = 0,
    justification = -0.3,
    trim = FALSE
  ) +
  geom_boxplot(
    aes(color = Screening_Status_24hr),
    width = 0.12,
    outlier.shape = NA
  ) +
  geom_half_point(
    aes(color = Screening_Status_24hr),
    side = "l",
    range_scale = 0.5,
    size = 0.8
  ) +
  geom_hline(yintercept = 12, col = "red", linetype = "dashed") +
  labs(y = "Age at Admission (years)", title = "Sample Size (n = 113)") +
  scale_color_manual(values = c("#66C2A5", "#FC8D62")) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  scale_y_continuous(breaks = seq(0, 20, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = c("#66C2A5", "#FC8D62")),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  coord_flip()

fig_age_screening

## ---- acuteLOS_screening_plot ----
# Acute care LOS by Screening Status
fig_AcuteLOS_screening <- CAPD_prog |> 
  ggplot(aes(x = Screening_Status_24hr, y = Acute_LOS)) +
  stat_halfeye(
    aes(fill = Screening_Status_24hr),
    alpha = 0.5,
    adjust = 0.5,
    width = 0.5,
    .width = 0,
    justification = -0.3,
    trim = FALSE
  ) +
  geom_boxplot(
    aes(color = Screening_Status_24hr),
    width = 0.12,
    outlier.shape = NA
  ) +
  geom_half_point(
    aes(color = Screening_Status_24hr),
    side = "l",
    range_scale = 0.5,
    size = 0.8
  ) +
  geom_hline(yintercept = 30, col = "red", linetype = "dashed") +
  labs(y = "Time From Acute Care Admission to Rehabilitaion Admission (days)",
       title = "Sample Size (n = 113)") +
  scale_color_manual(values = c("#66C2A5", "#FC8D62")) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  scale_y_continuous(breaks = seq(0, 150, 15)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = c("#66C2A5", "#FC8D62")),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  coord_flip()

fig_AcuteLOS_screening

## ---- save_fig ----
fig_path <- file.path(fig_dir, "fig_rain_AgeScreening_12cutoff.png")
ggsave(fig_path, plot = fig_age_screening, width = 8, height = 4)

fig_path <- file.path(fig_dir, "fig_rain_TTAScreening_30cutoff.png")
ggsave(fig_path, plot = fig_AcuteLOS_screening, width = 8, height = 4)

