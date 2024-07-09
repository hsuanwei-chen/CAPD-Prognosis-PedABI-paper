# Logistic Regression for Predicting WeeFIM Cog DFQ using 70 cutoff
# 
# This analysis script aims to evaluate whether screening status upon 24 hours
# of admission is predictive of WeeFIM Cog DFQ using 70 cutoff and offers useful 
# information in addition to demographic variables.
# 
# Created by Isaac Chen 2024.05.08

# Clear all the variables in the environment
rm(list = ls())

## ---- load_packages ----
lib_list = c("here", "dplyr", "forcats", "rms", "pROC", "gtsummary", "ggplot2", 
             "ggpubr")
lapply(lib_list, library, character.only = TRUE)

## ---- load_rawdata ----
data_path <- here("data", "processed_data", "CAPD_prognosis.rds")
CAPD_prog <- readRDS(data_path)

## ---- filter_excluded ----
# Filter out patients who are excluded from analysis and are missing WeeFIMs
CAPD_prog <- CAPD_prog |>
  filter(Excluded == "") |> 
  filter(!is.na(WeeFIM_CogDFQ_Adm) & !is.na(WeeFIM_CogDFQ_ultDC))

## ---- binary_outcome ----
# Create binary predictor variable and encode it as a factor
CAPD_prog <- CAPD_prog |>
  mutate(
    WeeFIM_CogDFQ_70cutoff = ifelse(
      WeeFIM_CogDFQ_ultDC >= 70, "Low Assistance", "Medium-High Assistance"
    )
  ) |> 
  mutate(
    WeeFIM_CogDFQ_70cutoff = factor(
      WeeFIM_CogDFQ_70cutoff, 
      levels = c("Low Assistance", "Medium-High Assistance")
    ),
  )

########## TABLE
## ---- tbl_contingency ----
# 2x2 table between screening status and transfer status
tbl_contingency <- table(
  fct_rev(CAPD_prog$Screening_Status_24hr), 
  fct_rev(CAPD_prog$WeeFIM_CogDFQ_70cutoff)
)
tbl_contingency

########### PRE-TEST MODEL
## ---- pre-test_SelectedModel_summary ----
# Create model with all variables included
pre_full_mod_rms <- 
  lrm(
    WeeFIM_CogDFQ_70cutoff ~ Age_Adm + Acute_LOS + Sex + Etiology_TBIOther + 
      Race, data = CAPD_prog, x= TRUE, y = TRUE
  )

# Variable selection via backward elimination
fastbw(pre_full_mod_rms, type = "individual")

# Reduced model
pre <- glm(
  WeeFIM_CogDFQ_70cutoff ~ Age_Adm + Acute_LOS, 
  data = CAPD_prog, family = "binomial"
)
pre_rms <- lrm(
  WeeFIM_CogDFQ_70cutoff ~ Age_Adm + Acute_LOS, 
  data = CAPD_prog, x= TRUE, y = TRUE
)

# Model summary
tbl_regression(pre, exponentiate = TRUE)

## ---- pre-test_SelectedModel_performance ----
# Statistical performance
pre_rms

# Discrimination Slope
tbl_pre <- data.frame(
  group = CAPD_prog$WeeFIM_CogDFQ_70cutoff,
  predicted = pre$fitted.values
)
pre_slope <- round(diff(tapply(tbl_pre$predicted, tbl_pre$group, mean)), 2)

## ---- pre-test_SelectedModel_IntVal ----
# Internal validation by bootstrapping
# C = Dx/2 + 0.5
set.seed(123)
validate(pre_full_mod_rms, B = 500, bw = TRUE, type = "individual")

########## POST-TEST MODEL
## ---- post-test_SelectedModel_summary ----
# Reduced model + screening status
post <- glm(
  WeeFIM_CogDFQ_70cutoff ~ Age_Adm + Acute_LOS + Screening_Status_24hr, 
  data = CAPD_prog, family = "binomial"
)
post_rms <- lrm(
  WeeFIM_CogDFQ_70cutoff ~ Age_Adm + Acute_LOS + Screening_Status_24hr, 
  data = CAPD_prog, x = TRUE, y = TRUE
)

# Model Summary
tbl_regression(post, exponentiate = TRUE)

## ---- post-test_SelectedModel_performance ----
# Statistical performance
post_rms

# Discrimination Slope
tbl_post <- data.frame(
  group = CAPD_prog$WeeFIM_CogDFQ_70cutoff,
  predicted = post$fitted.values
)
post_slope <- round(diff(tapply(tbl_post$predicted, tbl_post$group, mean)), 2)

## ---- post-test_SelectedModel_IntVal ----
# Internal validation by bootstrapping
# C = Dx/2 + 0.5
validate(post_rms, B = 500)

## ---- between_model_performance ----
anova(pre, post, test = "Chisq")

########## CHECK ASSUMPTIONS
## ---- post-test_AssumptionCheck ----
post_aug <- augment(post)
post_aug <- post_aug |>
  mutate(
    .predicted = post$fitted.values,
    .dfbeta = dfbeta(post),
    .lrg.resid = .std.resid > 2 | .std.resid < -2
  )

# Check the range of residuals and the number of cases with large residuals
range(post_aug$.std.resid)
sum(post_aug$.lrg.resid)

# Take a closer look at leverage and DFbeta for those cases
post_aug |>
  filter(.lrg.resid == T) |>
  select(.std.resid, .hat, .dfbeta)

# Check for multicolinearity
vif(post_rms)

########### VISUALIZATIONS
## ---- fig_box_ProbUnplanned_screening----
# Box plot of probability by transfer status for pre-test model
fig_box_pre <- tbl_pre |> 
  ggplot(aes(x = group, y = predicted)) +
  geom_boxplot() + 
  geom_point(
    data = aggregate(predicted ~ group, data = tbl_pre, mean),
    aes(x = group, y = predicted),
    col = "red"
  ) + 
  labs(
    title = "Pre-test Model",
    subtitle = paste0("Discrimination Slope = ", pre_slope),
    y = "Predicted Probabilites without CAPD"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme(axis.title.x = element_blank())

# Box plot of probability by transfer status for post-test model
fig_box_post <- tbl_post |> 
  ggplot(aes(x = group, y = predicted)) +
  geom_boxplot() + 
  geom_point(data = aggregate(predicted ~ group, data = tbl_post, mean),
             aes(x = group, y = predicted),
             col = "red") + 
  labs(title = "Post-test Model",
       subtitle = paste0("Discrimination Slope = ", post_slope),
       y = "Predicted Probabilites with CAPD") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1., 0.1)) +
  theme(axis.title.x = element_blank())

# Combine both plot into one figure
fig_box_prepost <- ggarrange(fig_box_pre, fig_box_post, labels = c("A", "B"))
fig_box_prepost

# Save figure
fig_path <- here(
  "results", "figures", "WeeFIM_70cutoff", "fig_box_PrePostProbabilites.png"
)
ggsave(fig_path, plot = fig_box_prepost, width = 5, height = 3)

## ---- fig_scatter_ProbUnplanned_screening ----
# Probability scatter plot
tbl_prepost <- data.frame(
  group = CAPD_prog$WeeFIM_CogDFQ_70cutoff,
  pre_predicted = pre$fitted.values,
  post_predicted = post$fitted.values
)

# Scatter plot of pre-test vs post-test probabilities
fig_scatter_prepost <- tbl_prepost|> 
  ggplot(aes(x = pre_predicted, y = post_predicted, col = fct_rev(group))) +
  geom_point() +
  geom_abline(linetype = "dashed", col = "grey") +
  labs(
    x = "Predicted Risks without Delirium Screening (Model 1)", 
    y = "Predicted Risks with Delirium Screening (Model 2)",
  ) + 
  scale_color_discrete(labels = c("Medium-High Assistance", "Low Assistance")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_bw() + 
  theme(legend.position = c(.8,.15),
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
fig_scatter_prepost

## ---- fig_hist_ProbUnplanned_screening----
# Probability back-to-back histogram
fig_hist_prepost <- tbl_prepost|> 
  ggplot() +
  geom_histogram(
    aes(x = pre_predicted, y = -after_stat(count)), 
    binwidth = 0.05, 
    color = "black", 
    fill= "#404080",
    boundary = 0
  ) +
  annotate(
    "text", x = 0.9, 
    y = -12, 
    label = "Model 1\n (Without Screenig)", 
    color="#404080", 
    size = 4
  ) +
  geom_histogram(
    aes(x = post_predicted, y = after_stat(count)), 
    binwidth = 0.05, 
    color = "black", 
    fill = "#69b3a2", 
    boundary = 0
  ) +
  annotate(
    "text",
    x = 0.9, 
    y = 15, 
    label = "Model 2\n (With Screening)", 
    color="#69b3a2", 
    size = 4
  ) +
  labs(
    x = "Predicted Risks of Medium-High Assistance", 
    y = "Count"
  ) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(-20, 45, 10)) + 
  theme_bw() + 
  coord_flip()
fig_hist_prepost

# Combine both plot into one figure
fig_scatter_hist_prepost <- ggarrange(fig_scatter_prepost, fig_hist_prepost, 
                                      labels = c("A", "B"))
fig_scatter_hist_prepost

# Save figure
fig_path <- here(
  "results", "figures", "WeeFIM_70cutoff", "fig_scatter_hist_PrePostProbabilites.png"
)
ggsave(fig_path, plot = fig_scatter_hist_prepost)
