###############################################
##### DOSE RESPONSE MODEL TEMPLATE SCRIPT #####
###############################################

# Clean the environment
rm(list = ls())

# Load packages
pacman::p_load(tidyverse,
               rstan,
               flextable,
               rio,
               janitor,
               RColorBrewer,
               loo)

# Load data
df <- import("") # add your own data here or import sample data below

# df <- import("example_data.csv", drop = 1)

head(df) # Look at the top few rows of the data

# clean data
## We need the following columns in the data, all cleaned:
### ~ Insecticide = name of insecticide
### ~ Concentration = dose of insecticide
### ~ Test_n = Number of mosquitoes tested
### ~ Test_mort_n = Number of mosquitoes that died
### ~ Test_mort_perc = Percentage mortality (=100*number died/number tested)
### ~ rlookup = A unique name for saving files
### ~ OPTIONAL: Covariate = A variable to test difference in dose response (site, year, species etc.)

df <- df %>% # This is how we would clean the example dataset
  rename(Concentration = Dose,
         Test_mort_perc = Mortality_perc,
         Test_n = Subjects,
         Test_mort_n = Responded) %>%
  mutate(rlookup = "example")

#######################
##### FIXED MODEL #####
#######################

df_fixed <- df

# Source the model code
source("dose_response_fixed.R")

# Set up for Stan model
options(mc.cores=4)
rstan::rstan_options(auto_write = TRUE)

# Run the model: this can take some time
dose_response_fixed(df_fixed)

## Model checks
model_checks <- read.csv("outputs/csv/modelcheck_fixed_example.csv")

# For a good model, we want 0 not converged and 0 divergent iterations:
model_checks

### Extract dose-response data
data_fixed <- read.csv("outputs/csv/df1_s_fixed_example.csv")

### Clean columns
data_fixed <- data_fixed %>%
  mutate(Test_mort_perc=ifelse(is.na(Test_mort_perc),Mortality_perc,Test_mort_perc))

#### LOAD LC DATA FOR SUMMARY TABLE
LC50_fixed <- read.csv("outputs/csv/fixed_summ_LC50_example.csv")
LC99_fixed <- read.csv("outputs/csv/fixed_summ_LC99_example.csv")

# Combine LC50 and LC99 data together
LC_summary_fixed <- bind_rows(LC50_fixed, LC99_fixed) %>%
  select(-X, -LC_mean) %>%
  pivot_wider(names_from = LC_value, values_from = c(LC_median:LC_upper), names_vary = "slowest") %>%
  mutate(across(LC_median_50:LC_upper_99, ~ round(.x,3))) %>%
  mutate(Insecticide = Insecticide,
         LC50 = LC_median_50,
         `LC50 (95%CI)` = paste0(LC_lower_50,"-",LC_upper_50),
         LC99 = LC_median_99,
         `LC99 (95%CI)` = paste0(LC_lower_99,"-",LC_upper_99),
         .keep = "none")

# Save or make summary table for LC50 and LC99
write.csv(LC_summary_fixed, "outputs/csv/LC_summary_table_fixed.csv")

LC_summary_fixed %>%
  flextable() %>%
  bold(part = "header")

## PLOT
ggplot(data = data_fixed,
       aes(x = Concentration)) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.border = element_rect(fill = NA),
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.5),
        legend.box="vertical") +
  geom_point(data = filter(data_fixed, is.na(dat)),
             aes(y = Test_mort_perc,
                 color = Insecticide)) +
  geom_ribbon(data = filter(data_fixed, !is.na(dat)),
              aes(ymin=lower, ymax=upper,
                  fill = Insecticide), alpha=0.2) +
  geom_line(data = filter(data_fixed, !is.na(dat)),
            aes(y=Test_mort_perc,
                color = Insecticide)) +
  ylab("Mosquito mortality (%)") +
  xlab("Insecticide concentration (%)") +
  guides(color = guide_legend(title = "Insecticide"),
         fill = guide_legend(title = "Insecticide")) +
  scale_x_sqrt()

ggsave("outputs/plots/fixed_model_dose_response.png", bg = "white",
       width = 8, height = 8)

################################
##### LINEAR MODEL BY YEAR #####
################################

# Data cleaning
df_linear <- df %>% # For the linear model, we need to select the covariate, in this case Year
  mutate(Covariate = year)

# Source the function for the linear model
source("dose_response_linear.R")

# Run the model: this can take some time
dose_response_linear(df_linear)

# Model checks
model_checks <- read.csv("outputs/csv/modelcheck_linear_example.csv")

# For a good model, we want 0 not converged and 0 divergent iterations:
model_checks

# Extract dose-response data
data_linear <- read.csv("outputs/csv/df1_s_linear_example.csv")

# Clean columns
data_linear <- data_linear %>%
  mutate(Test_mort_perc=ifelse(is.na(Test_mort_perc),Mortality_perc,Test_mort_perc))

# LOAD LC DATA FOR SUMMARY TABLE
LC50_linear <- read.csv("outputs/csv/linear_summ_LC50_example.csv")
LC99_linear <- read.csv("outputs/csv/linear_summ_LC99_example.csv")

# Combine LC50 and LC99 data together
LC_summary_linear <- bind_rows(LC50_linear, LC99_linear) %>%
  select(-X, -LC_mean) %>%
  pivot_wider(names_from = LC_value, values_from = c(LC_median:LC_upper), names_vary = "slowest") %>%
  mutate(across(LC_median_50:LC_upper_99, ~ round(.x,3))) %>%
  mutate(Year = Covariate,
         Insecticide = Insecticide,
         LC50 = LC_median_50,
         `LC50 (95%CI)` = paste0(LC_lower_50,"-",LC_upper_50),
         LC99 = LC_median_99,
         `LC99 (95%CI)` = paste0(LC_lower_99,"-",LC_upper_99),
         .keep = "none")

# Save or make summary table for LC50 and LC99
write.csv(LC_summary_linear, "outputs/csv/LC_summary_table_linear.csv")

LC_summary_linear %>%
  flextable() %>%
  bold(part = "header")

# PLOT
ggplot(data = data_linear,
       aes(x = Concentration)) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.border = element_rect(fill = NA),
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.5),
        legend.box="vertical") +
  geom_point(data = filter(data_linear, is.na(dat)),
             aes(y = Test_mort_perc,
                 color = factor(Covariate))) +
  geom_ribbon(data = filter(data_linear, !is.na(dat)),
              aes(ymin=lower, ymax=upper,
                  group = factor(Covariate),
                  fill = factor(Covariate)), alpha=0.2) +
  geom_line(data = filter(data_linear, !is.na(dat)),
            aes(y=Test_mort_perc,
                group = factor(Covariate),
                color = factor(Covariate))) +
  ylab("Mosquito mortality (%)") +
  xlab("Insecticide concentration (%)") +
  guides(color = guide_legend(title = "Year"),
         fill = guide_legend(title = "Year")) +
  scale_x_sqrt()

ggsave("outputs/plots/linear_model_dose_response.png", bg = "white",
       width = 8, height = 8)

####################################
##### INDIVIDUAL MODEL BY TIME #####
####################################

# For this model we want to assess the difference in insecticide resistance betwee
# different strains of mosquitoes, so we need to change the covariate variable to Strain
df_individual <- df %>%
  mutate(Covariate = Strain)

# Source the function for the individual model
source("dose_response_individual.R")

# Run the model: this can take some time
dose_response_individual(df_individual)

## Model checks
model_checks <- read.csv("outputs/csv/modelcheck_individual_example.csv")

# For a good model, we want 0 not converged and 0 divergent iterations:
model_checks

### Extract dose-response data
data_individual <- read.csv("outputs/csv/df1_s_individual_example.csv")

### Clean columns
data_individual <- data_individual %>%
  mutate(Test_mort_perc=ifelse(is.na(Test_mort_perc),Mortality_perc,Test_mort_perc))

#### LOAD LC DATA FOR SUMMARY TABLE
LC50_individual <- read.csv("outputs/csv/individual_summ_LC50_example.csv")
LC99_individual <- read.csv("outputs/csv/individual_summ_LC99_example.csv")

# Combine LC50 and LC99 data together
LC_summary_individual <- bind_rows(LC50_individual, LC99_individual) %>%
  select(-X, -LC_mean) %>%
  pivot_wider(names_from = LC_value, values_from = c(LC_median:LC_upper), names_vary = "slowest") %>%
  mutate(across(LC_median_50:LC_upper_99, ~ round(.x,3))) %>%
  mutate(Strain = Covariate,
         Insecticide = Insecticide,
         LC50 = LC_median_50,
         `LC50 (95%CI)` = paste0(LC_lower_50,"-",LC_upper_50),
         LC99 = LC_median_99,
         `LC99 (95%CI)` = paste0(LC_lower_99,"-",LC_upper_99),
         .keep = "none")

# Save or make summary table for LC50 and LC99
write.csv(LC_summary_individual, "outputs/csv/LC_summary_table_individual.csv")

LC_summary_individual %>%
  flextable() %>%
  bold(part = "header")

## PLOT
ggplot(data = data_individual,
       aes(x = Concentration)) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.border = element_rect(fill = NA),
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.5),
        legend.box="vertical") +
  geom_point(data = filter(data_individual, is.na(dat)),
             aes(y = Test_mort_perc,
                 color = factor(Covariate))) +
  geom_ribbon(data = filter(data_individual, !is.na(dat)),
              aes(ymin=lower, ymax=upper,
                  group = factor(Covariate),
                  fill = factor(Covariate)), alpha=0.2) +
  geom_line(data = filter(data_individual, !is.na(dat)),
            aes(y=Test_mort_perc,
                group = factor(Covariate),
                color = factor(Covariate))) +
  ylab("Mosquito mortality (%)") +
  xlab("Insecticide concentration (%)") +
  guides(color = guide_legend(title = "Strain"),
         fill = guide_legend(title = "Strain")) +
  scale_x_sqrt()

ggsave("outputs/plots/linear_model_dose_individual.png", bg = "white",
       width = 8, height = 8)

############################
##### MODEL COMPARISON #####
############################

# Model comparison uses leave-one-out cross-validation (LOO-CV)
# This essentially takes one data point from each model and removes it,
# fits a model with the rest of the data and then sees how well the removed data
# point fits the model. It then repeats this for all data points one at a time.
# The version in the LOO package approximates this rather than actually computing it

# The higher the elpd the better the model fit

# Fixed vs linear to see if insecticide resistance different over time
LOO_fixed <- readRDS("outputs/LOO/fixed_example.rds")
LOO_linear <- readRDS("outputs/LOO/linear_example.rds")

# Next we directly compare the two models, if there is a big difference in elpd between the models
# and this difference is much bigger than the standard error between them, then this suggests the model with the higher elpd is the best

LOO_compare_year <- as.data.frame(loo_compare(list("fixed" = LOO_fixed,
                                                    "linear" = LOO_linear)))

# To formally test  we can compute a p-value using a z-score and the standard normal distribution
# to see if the difference in elpd is significantly bigger than the difference in standard error

LOO_compare_year$p_value <- pnorm(LOO_compare_year$elpd_diff/LOO_compare_year$se_diff, lower.tail = TRUE) # = p-value
LOO_compare_year

# And we can repeat the process for for individual model to see if there is a difference between Strains in resistance
LOO_individual <- readRDS("outputs/LOO/individual_example.rds")

LOO_compare_strain <- as.data.frame(loo_compare(list("fixed" = LOO_fixed,
                                                   "individual" = LOO_individual)))

LOO_compare_strain$p_value <- pnorm(LOO_compare_strain$elpd_diff/LOO_compare_strain$se_diff, lower.tail = TRUE) # = p-value
LOO_compare_strain
