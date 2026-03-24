# Load required libraries
library(dplyr)
library(tidyr)

# Source the generate.counterfactor function
source("generateCounterfactor.R")

# Load obsdata5
obsdata5 <- readRDS("obsdata5.rds")

# Extract baseline (time=0) for each subject
baselines <- obsdata5 %>%
  filter(follow_up == 0) %>%
  select(id, X1, X2, X3, X4, age) %>%
  rename(X1_0 = X1, X2_0 = X2, age_0 = age)


# Truth for pp
J <- 120  # number of months
beta.adherence <- 100
beta0.outcome <- -5
beta.trt.effect <- -0.8
beta0.censoring <- -100

counterfactual_values <- c(0, 1)
scenario_results <- vector("list", length(counterfactual_values))

for (scenario_index in seq_along(counterfactual_values)) {
  counterfactual.A <- counterfactual_values[scenario_index]
  counterfactor_list <- vector("list", nrow(baselines))

  for (i in seq_len(nrow(baselines))) {
    subject_baseline <- baselines[i, , drop = FALSE]

    counterfactor_list[[i]] <- generate.counterfactor(
      J = J,
      baseline = subject_baseline,
      A0 = counterfactual.A,
      beta.adherence = beta.adherence,
      beta0.outcome = beta0.outcome,
      beta.trt.effect = beta.trt.effect,
      beta0.censoring = beta0.censoring,
      seed = NULL
    )

    if (i %% 10 == 0) {
      cat(sprintf("Processed %d/%d subjects for counterfactual.A = %d\n", i, nrow(baselines), counterfactual.A))
    }
  }

  scenario_results[[scenario_index]] <- do.call(rbind, counterfactor_list) %>%
    mutate(counterfactual.A = counterfactual.A)
}

# Combine all counterfactual trajectories across both scenarios
counterfactor_data <- bind_rows(scenario_results)

# Calculate the observed event probability at each time point within each counterfactual scenario.
# The denominator is the number of subjects still present in the generated data at that time point.
event_probability_by_time <- counterfactor_data %>%
  group_by(counterfactual.A, time) %>%
  summarise(prob_outcome = mean(Y, na.rm = TRUE), .groups = "drop") %>%
  group_by(counterfactual.A) %>%
  mutate(
    survival_prob = cumprod(1 - prob_outcome), # P(alive up to current time)
    survival_prob_lag = lag(survival_prob, default = 1), # P(alive up to last time)
    death_prob = survival_prob_lag * prob_outcome, # P(dead at current time)
    cumulative_incidence = cumsum(death_prob)
  ) %>%
  ungroup()

cumulative_incidence_wide <- event_probability_by_time %>%
  select(time, counterfactual.A, cumulative_incidence) %>%
  pivot_wider(
    names_from = counterfactual.A,
    values_from = cumulative_incidence,
    names_prefix = "cumulative_incidence_trt"
  ) %>% mutate(risk_difference = cumulative_incidence_trt1 - cumulative_incidence_trt0,
    risk_ratio = cumulative_incidence_trt1 / cumulative_incidence_trt0)
risk.truth.pp <- cumulative_incidence_wide


# Truth for ITT
J <- 120  # Number of months
beta.adherence <- 6
beta0.outcome <- -5
beta.trt.effect <- -0.8
beta0.censoring <- -3
counterfactual_values <- c(0, 1)
scenario_results <- vector("list", length(counterfactual_values))

for (scenario_index in seq_along(counterfactual_values)) {
  counterfactual.A <- counterfactual_values[scenario_index]
  counterfactor_list <- vector("list", nrow(baselines))

  for (i in seq_len(nrow(baselines))) {
    subject_baseline <- baselines[i, , drop = FALSE]

    counterfactor_list[[i]] <- generate.counterfactor(
      J = J,
      baseline = subject_baseline,
      A0 = counterfactual.A,
      beta.adherence = beta.adherence,
      beta0.outcome = beta0.outcome,
      beta.trt.effect = beta.trt.effect,
      beta0.censoring = beta0.censoring,
      seed = NULL
    )

    if (i %% 10 == 0) {
      cat(sprintf("Processed %d/%d subjects for counterfactual.A = %d\n", i, nrow(baselines), counterfactual.A))
    }
  }

  scenario_results[[scenario_index]] <- do.call(rbind, counterfactor_list) %>%
    mutate(counterfactual.A = counterfactual.A)
}

# Combine all counterfactual trajectories across both scenarios
counterfactor_data <- bind_rows(scenario_results)

# Calculate the observed event probability at each time point within each counterfactual scenario.
# The denominator is the number of subjects still present in the generated data at that time point.
event_probability_by_time <- counterfactor_data %>%
  group_by(counterfactual.A, time) %>%
  summarise(prob_outcome = mean(Y, na.rm = TRUE), .groups = "drop") %>%
  group_by(counterfactual.A) %>%
  mutate(
    survival_prob = cumprod(1 - prob_outcome), # P(alive up to current time)
    survival_prob_lag = lag(survival_prob, default = 1), # P(alive up to last time)
    death_prob = survival_prob_lag * prob_outcome, # P(dead at current time)
    cumulative_incidence = cumsum(death_prob)
  ) %>%
  ungroup()

cumulative_incidence_wide <- event_probability_by_time %>%
  select(time, counterfactual.A, cumulative_incidence) %>%
  pivot_wider(
    names_from = counterfactual.A,
    values_from = cumulative_incidence,
    names_prefix = "cumulative_incidence_trt"
  ) %>% mutate(risk_difference = cumulative_incidence_trt1 - cumulative_incidence_trt0,
    risk_ratio = cumulative_incidence_trt1 / cumulative_incidence_trt0)
risk.truth.ITT <- cumulative_incidence_wide

saveRDS(list(PP.truth = risk.truth.pp, ITT.truth = risk.truth.ITT), "risk_truth.rds")
