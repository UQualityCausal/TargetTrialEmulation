# We developed the function below to simulate an observational dataset  
# Using simulated data offers key advantages: 
# it avoids data security concerns and allows us to know the true treatment effects, 
# enabling us to evaluate the performance of different analytic strategies, including TTE and causal inference methods. 


library(dplyr)

# Input Parameters:
# n: Number of patients
# J: Number of follow-up time points (in years)
# beta0.trt0: Intercept for initial treatment assignment
# beta.adherence: Coefficient for treatment adherence
# beta0.outcome: Intercept for outcome model
# beta.trt.effect: Effect of treatment on outcome
# beta0.censoring: Intercept for censoring model

generate.obs.data <- function(n, J, beta0.trt0, beta.adherence, beta0.outcome, beta.trt.effect, beta0.censoring) {
  
  # Create baseline data
  baseline <- data.frame(
    id = 1:n, 
    X1_0 = rbinom(n, 1, 0.5), 
    X2_0 = rnorm(n, 0, 1), 
    age_0 = rnorm(n, 35, 12),
    X3 = rbinom(n, 1, 0.5),
    X4 = rnorm(n, 0, 1)) %>% 
    rowwise() %>%
    mutate(logit_pi0 = beta0.trt0 + 0.5 * X1_0 + 0.5 * X2_0 - 0.2 * X3 + X4 - 0.03 * age_0,
           A0 = rbinom(1, 1, plogis(logit_pi0))) %>%
    select(-logit_pi0) %>% 
    ungroup()
  
  # Expand baseline data to long format
  long_data <- baseline %>%
    slice(rep(1:n(), each = J)) %>% 
    mutate(time = rep(0:(J - 1), times = n)) %>% # Add a follow-up time variable for each patient
    group_by(id) %>%
    mutate(
      age = age_0 + time, 
      age_s = (age - 35) / 12, # Standardized age
      A_prev = NA,
      Y_prev = NA,
      C_prev = NA,
      A=if_else(time == 0, A0, NA),
      Y = if_else(time == 0, 0, NA),
      C = if_else(time == 0, 0, NA),
      eligible = 1,
      X1 = if_else(time == 0, X1_0, NA),
      X2 = if_else(time == 0, X2_0, NA)) %>%
    ungroup()
  
  # Simulate time-varying covariates, treatment, outcome, censoring, and eligibility
  long_data <- long_data %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x # Refer to the subset of rows for one patient
      for (j in 2:nrow(df)) {
        
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j - 1])) # plogis(-0) = 0.5, plogis(-1) = 0.2689414
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j - 1], sd = 1) 
        df$A_prev[j] <- df$A[j - 1]
        df$Y_prev[j] <- df$Y[j - 1]
        df$C_prev[j] <- df$C[j - 1]
        
        # ?? Treatment or treatment adherence: 
        # The probability of adhering to treatment A at time j is a function of treatment adherence at time j-1 and covariates at time j
        logit_pi <- beta.adherence* (df$A_prev[j]-0.5 ) + 0.5 * df$X1[j] + 0.5 * df$X2[j] - 0.2 * df$X3[j] + df$X4[j] - 0.3 * df$age_s[j]
        df$A[j] <- rbinom(1, 1, plogis(logit_pi)) 
        
        # Eligibility criteria: at least 18 years old, have not received treatment, and have not experienced the event of interest
        df$eligible[j] <- as.integer(df$age[j] >= 18 & all(df$A[1:j - 1] == 0) & all(df$Y[1:j - 1] == 0))
        
        # Outcome
        if (df$Y_prev[j] == 0) {
          logit_lambda <- beta0.outcome + beta.trt.effect * df$A[j] + 0.5 * df$X1[j] + 0.5 * df$X2[j] + df$X3[j] + df$X4[j] + 0.5 * df$age_s[j]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        } else {
          df$Y[j] <- 1  # once outcome is positive, remain positive
        }
        
        # Censoring
        if (df$C_prev[j] == 0){
          logit_q <- beta0.censoring - df$A_prev[j] - 0.5 * df$X1[j] + 0.5 * df$X2[j] - 0.2 * df$X3[j] + 0.2 * df$X4[j] - df$age_s[j]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } else {
          df$C[j] <- 1  # once censored, remain censored
        }
      }
      df
    }) %>%
    ungroup()
  
  # Drop data after censoring or event
  t0 <- long_data %>%
    filter(C==1 | Y==1) %>%
    group_by(id) %>%
    summarise(t0=min(time))
  
  long_data <- long_data %>%
    left_join(t0) %>%
    filter(is.na(t0)| time<=t0) %>%
    select(-t0)
  
  # Drop data before first eligible
  long_data <- long_data %>%
    group_by(id) %>%
    mutate(first_eligible = min(time[eligible == 1], na.rm = TRUE)) %>%
    filter(time >= first_eligible) %>%
    ungroup() %>%
    select(-first_eligible)
  
  # Ensure all variables are numeric
  long_data <- long_data %>%
    mutate(across(c(id, time, A, Y, C, eligible, X1, X2, X3, X4, age_s), as.numeric))
  
  long_data
}

# Set seed for reproducibility
set.seed(412)
simdata0 <- generate.obs.data(n=1000, J=20, beta0.trt0=1, beta.adherence=2, 
                              beta0.outcome=-5, beta.trt.effect=-1.2, beta0.censoring=-1)
obsdata <- simdata0 %>%
  select(id, time, X1, X2, X3, X4, age, A, Y, C, age_s)

attr(obsdata$age, "label") <- "Age over time (years)"
attr(obsdata$X3, "label") <- "Biological sex (F/M)"
attr(obsdata$X4, "label") <- "Standardized systolic blood pressure at baseline"
attr(obsdata$id, "label") <- "Patient ID"
attr(obsdata$time, "label") <- "Time index for longitudinal records"
attr(obsdata$age_s, "label") <- "Standardized age over time (years)"
attr(obsdata$X1, "label") <- "Non-ACEI or ARB antihypertensive medication use over time"
attr(obsdata$X2, "label") <- "Standardized systolic blood pressure over time"
attr(obsdata$A, "label") <- "Treatment indicator of ARB use over time (ACEI=0)"
attr(obsdata$Y, "label") <- "Event indicator of cardiovascular disease"
attr(obsdata$C, "label") <- "Indicator of early dropout"

saveRDS(obsdata, "obsdata.rds")
