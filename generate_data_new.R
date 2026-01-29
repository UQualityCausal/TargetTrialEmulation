# We developed the function below to simulate an observational dataset  
# Using simulated data offers key advantages: 
# it avoids data security concerns and allows us to know the true treatment effects, 
# enabling us to evaluate the performance of different analytic strategies, including TTE and causal inference methods. 

library(dplyr)
library(lubridate)

# Input Parameters:
# n: Number of patients
# J: Number of follow-up time points (in months)
# beta0.trt0: Intercept for initial treatment assignment
# beta.adherence: Coefficient for treatment adherence
# beta0.outcome: Intercept for outcome model
# beta.trt.effect: Effect of treatment on outcome
# beta0.censoring: Intercept for censoring model
# deviation_prob: Probability that a subject will deviate from treatment

generate.obs.data <- function(n, J, beta0.trt0 = 0, beta.adherence = 1, 
                              beta0.outcome = -6.5, beta.trt.effect = -1.2, 
                              beta0.censoring = -3.5, deviation_prob = 0.1) {
  
  # Initialize baseline data
  baseline <- data.frame(
    id = 1:n, # Subject identifier
    X1_0 = rbinom(n, 1, 0.5), # X1: Time-varying categorical variable
    X2_0 = rnorm(n, 0, 1), # X2: Time-varying numeric variable
    Age_0 = rnorm(n, 35, 12), # Subject Age
    X3 = rbinom(n, 1, 0.5), # X3: Fixed categorical variable
    X4 = rnorm(n, 0, 1) # X4: Fixed numeric variable 
  ) %>%
    rowwise() %>%
    mutate(
      logit_pi0 = beta0.trt0 + 0.5* X1_0 + 0.5 * X2_0 - 0.2 * X3 + X4 - 0.3 * (Age_0-35)/12,
      A_0 = rbinom(1, 1, plogis(logit_pi0)) # Use initial treatment probability plogis(logit_pi0) to determine treatment A
    ) %>%
    select(-logit_pi0) %>% # Remove intermediate variable logit_pi0
    ungroup()
  
  # Expand data to long format: one row per person per time
  long_data <- baseline %>%
    slice(rep(1:n(), each = J)) %>% # Each patient has J (time points) rows
    mutate(Time = rep(0:(J-1), times = n)) %>% # Add a time variable that repeats the time variable for each patients
    group_by(id) %>%
    mutate(
      Age = Age_0 + Time / 12, # Unit for Time is month
      Age_s = (Age - 35) / 12, # Standardized Age
      A = if_else(Time == 0, A_0, NA),
      Y = if_else(Time == 0, 0, NA), # Y=0: negative outcome
      C = if_else(Time == 0, 0, NA), # C=0: not censored
      E = if_else(Time == 0, 1, NA), # E=1: eligible
      X1 = if_else(Time == 0, X1_0, NA),
      X2 = if_else(Time == 0, X2_0, NA)
    ) %>% ungroup()
  
  # Simulate time-varying covariates, treatment, outcome, censoring
  long_data <- long_data %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x # Refer to the subset of rows for one subject
      X3 <- df$X3[1]
      X4 <- df$X4[1]
      
      for (j in 2:nrow(df)) {
        
        # Find last non-NA treatment up to previous time
        A_values_to_prev <- df$A[1:(j-1)]
        A_prev <- tail(A_values_to_prev[!is.na(A_values_to_prev)], 1)
        
        # Update time-varying covariates
        df$X1[j] <- rbinom(1, 1, plogis(-A_prev)) # plogis(-0) = 0.5, plogis(-1) = 0.2689414
        df$X2[j] <- rnorm(1, mean = -0.3 * A_prev, sd = 1)
        
        if (rbinom(1, 1, deviation_prob) == 1) {
          df$A[j] <- NA  # Treatment deviation
        } else {
          # Treatment model: probability of adhering to treatment A at time j is a function of treatment adherence at time j-1 and covariates at time j
          logit_pi <- beta.adherence * A_prev + 0.5 * df$X1[j] + 
            0.5 * df$X2[j] - 0.2 * X3 + X4 - 0.3 * df$Age_s[j]
          df$A[j] <- rbinom(1, 1, plogis(logit_pi))
        }
        
        # Find last non-NA treatment up to current time
        A_values_to_current <- df$A[1:j]
        A_current <- tail(A_values_to_current[!is.na(A_values_to_current)], 1)
        
        # Outcome model (e.g., death) 
        if (df$Y[j - 1] == 0) {
          logit_lambda <- beta0.outcome + beta.trt.effect * A_current + 0.5 * df$X1[j] + 
            0.5 * df$X2[j] + X3 + X4 + 0.5 * df$Age_s[j]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        }
        else{ 
          df$Y[j] <- 1 # Once outcome is positive, remain positive
        }
        
        # Censoring model
        if (df$C[j - 1] == 0) { 
          logit_q <- beta0.censoring - A_prev - 0.5 * df$X1[j] + 
            0.5 * df$X2[j] - 0.2 * X3 + 0.2 * X4 - df$Age_s[j]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } 
        else{ 
          df$C[j] <- 1 # Once censored, remain censored
        }
      }
      df
    }) %>%
    ungroup()
  
  # Drop data after censoring or event
  t0 <- long_data %>% 
    filter(C == 1 | Y == 1) %>% 
    group_by(id) %>% 
    summarise(t0 = min(Time))
  
  long_data <- long_data %>% 
    left_join(t0, by = "id") %>% 
    filter(is.na(t0) | Time <= t0) %>% 
    select(-t0)
  
  return(long_data)
}

# Extend data with pre-treatment period
extend.data <- function(data) {
  months.prev <- sample(2:5, 1) # Add 2-5 periods of pre-treatment history
  data.prev <- data[rep(1, months.prev), ] %>% # Replicate the first row months.prev times to keep fixed variables
    mutate(
      Time = -(1:months.prev), # Negative time for pre-treatment period
      A = NA, # No treatment
      Y = 0, # No outcome event
      C = 0, # Not censored
      X1 = rbinom(n(), 1, 0.5), # Random Time-varying covariates
      X2 = rnorm(n(), 0, 1), # Random Time-varying covariates
      Age = data$Age[1] - (1:months.prev)/12 # Age counts backwards
    )
  
  data.extended <- rbind(data.prev, data) %>% 
    arrange(id, Time)
  
  return(data.extended)
}

# Input Parameters:
# n: Number of patients
# J: Number of time points
# n_ineligible: Number of subjects who never initiate treatment
generate.final.data <- function(n, J, n_ineligible = 100) {
  
  simdata <- generate.obs.data(n = n, J = J)
  
  # Extend data for all subjects to add pre-treatment period
  simdata.extended <- simdata %>% 
    mutate(id0 = id) %>% 
    group_by(id0) %>%
    group_modify(~extend.data(.x)) %>%
    ungroup() %>%
    select(-id0)
  
  # Add ineligible subjects who never initiate either treatment
  # These subjects will be filtered out during data processing
  simdata.no.trt <- generate.obs.data(n = n_ineligible, J = J)
  simdata.no.trt <- simdata.no.trt %>% 
    mutate(A = NA, id = id + max(simdata.extended$id)) 
  
  # Combine eligible and ineligible subjects
  simdata.final <- rbind(simdata.extended, simdata.no.trt) %>%
    arrange(id, Time)
  
  start_year <- sample(2014:2023, 1)
  start_date <- as.Date(paste(start_year, "01", "01", sep = "-"))
  
  simdata.final <- simdata.final %>%
    group_by(id) %>%
    mutate(Time = format(start_date + months(Time), "%Y-%m")) %>%
    ungroup()
  
  return(simdata.final)
}

# Set seed for reproducibility
set.seed(412)

simdata0 <- generate.final.data(n = 1000, J = 80)
obsdata <- simdata0 %>% select(id, Time, X1, X2, X3, X4, Age, A, Y, C) 

attr(obsdata$id, "label") <- "Patient ID"
attr(obsdata$Time, "label") <- "Time index for longitudinal records"
attr(obsdata$X1, "label") <- "Non-ACEI or ARB antihypertensive medication use over time"
attr(obsdata$X2, "label") <- "Standardized systolic blood pressure over time"
attr(obsdata$X3, "label") <- "Biological sex (F/M)"
attr(obsdata$X4, "label") <- "Standardized diastolic blood pressure at baseline"
attr(obsdata$Age, "label") <- "Age over time (years)"
attr(obsdata$A, "label") <- "Treatment indicator of ARB use over time (ACEI=0)"
attr(obsdata$Y, "label") <- "Event indicator of cardiovascular disease"
attr(obsdata$C, "label") <- "Indicator of early dropout"

saveRDS(obsdata, "obsdata.rds")