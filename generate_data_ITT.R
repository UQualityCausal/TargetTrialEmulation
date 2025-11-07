library(dplyr)

generate.obs.data <- function(n, J) {
  ###### Input Parameters: 
  # n: Number of subjects
  # J: Number of time points
  ###### 
  
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
      logit_pi0 = 0.5* X1_0 + 0.5 * X2_0 - 0.2 * X3 + X4 - 0.3 * (Age_0-35)/12,
      A_0 = rbinom(1, 1, plogis(logit_pi0)) # Use initial treatment probability plogis(logit_pi0) to determine treatment A
    ) %>%
    select(-logit_pi0) %>% # Remove intermediate variable logit_pi0
    ungroup()
  
  # Expand to long format: one row per person per time
  long_data <- baseline %>%
    slice(rep(1:n(), each = J)) %>% # Each patient has J (time points) rows
    mutate(Time = rep(0:(J-1), times = n)) %>% # Add a time variable that repeats the time points sequence n times (for n patients)
    group_by(id) %>%
    mutate(
      Age = Age_0 + Time, # Unit for Time is year
      Age_s = (Age - 35) / 12, # Standardized Age
      A = if_else(Time == 0, A_0, NA),
      Y = if_else(Time == 0, 0, NA), # Y=0: negative outcome
      C = if_else(Time == 0, 0, NA), # C=0: not censored
      E = if_else(Time == 0, 1, NA), # E=1: eligible
      X1 = if_else(Time == 0, X1_0, NA),
      X2 = if_else(Time == 0, X2_0, NA)
    ) %>% ungroup()
  
  # Simulate time-varying data
  long_data <- long_data %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x # Refer to the subset of rows for one subject
      X3 <- df$X3[1]
      X4 <- df$X4[1]
      
      for (j in 2:nrow(df)) {
        
        # Update time-varying covariates
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j - 1])) # plogis(-0) = 0.5, plogis(-1) = 0.2689414
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j - 1], sd = 1)
        
        # Treatment model
        logit_pi <- df$A[j - 1] + 0.5 * df$X1[j] + 
          0.5 * df$X2[j] - 0.2 * X3 + X4 - 0.3 * df$Age_s[j]
        df$A[j] <- rbinom(1, 1, plogis(logit_pi))
        
        # Outcome model (e.g., death) 
        if (df$Y[j - 1] == 0) {
          logit_lambda <- -5 - 1.2 * df$A[j] + 0.5 * df$X1[j] + 
            0.5 * df$X2[j] + X3 + X4 + 0.5 * df$Age_s[j]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        }
        else{ 
          df$Y[j] <- 1 # Once outcome is positive, remain positive
        }
        
        # Censoring model
        if (df$C[j - 1] == 0) { 
          logit_q <- -1 - df$A[j - 1] - 0.5 * df$X1[j] + 
            0.5 * df$X2[j] - 0.2 * X3 + 0.2 * X4 - df$Age_s[j]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } 
        else{ 
          df$C[j] <- 1 # Once censored, remain censored
        }
        
        # # Eligibility criteria: at least 18 years old, have not experienced outcome of interest/censoring, have not received treatment in the past two years
        # prev_Y_zero <- all(df$Y[1:(j-1)] == 0, na.rm = TRUE)
        # prev_C_zero <- all(df$C[1:(j-1)] == 0, na.rm = TRUE)
        # no_recent_A <- (df$A[j-1] + ifelse(j > 2, df$A[j-2], 0)) == 0
        # df$E[j] <- as.integer(df$Age[j] >= 18 & prev_Y_zero & 
        #                         prev_C_zero & no_recent_A)
        
      }
      df
    }) %>%
    ungroup()
  
  # Remove rows after censoring or outcome event happened
  t0 <- long_data %>% filter(C==1 | Y==1) %>% 
    group_by(id) %>% summarise(t0=min(Time))
  long_data <- long_data %>% left_join(t0, by = "id") %>% 
    filter(is.na(t0) | Time<=t0) %>% select(-t0)
  
  # # Remove data before first eligible time
  # long_data <- long_data %>% group_by(id) %>%
  #   mutate(first_eligible = ifelse(any(E == 1, na.rm = TRUE), 
  #                                  min(Time[E == 1], na.rm = TRUE), NA)) %>%
  #   filter(!is.na(first_eligible) & Time >= first_eligible) %>%
  #   ungroup() %>% select(-first_eligible)
  
  return(long_data)
}

# Extend data with pre-treatment period
extend.data <- function(data) {
  yrs.prev <- sample(2:5, 1) # Add 2-5 years of pre-treatment history
  data.prev <- data[rep(1, yrs.prev), ] %>% # Replicate the first row yrs.prev times to keep fixed variables
    mutate(
      Time = -(1:yrs.prev), # Negative time for pre-treatment period
      A = NA, # No treatment
      Y = 0, # No outcome event
      C = 0, # Not censored
      X1 = rbinom(n(), 1, 0.5), # Random Time-varying covariates
      X2 = rnorm(n(), 0, 1), # Random Time-varying covariates
      Age = data$Age[1] - (1:yrs.prev) # Age counts backwards
    )
  
  data.extended <- rbind(data.prev, data) %>% 
    arrange(id, Time)
  
  return(data.extended)
}

generate.final.data <- function(n, J, n_ineligible = 100) {
  ###### Input Parameters: 
  # n: Number of eligible subjects
  # J: Number of time points
  # n_ineligible: Number of subjects who never initiate treatment
  # A=0/1: Two active treatment arms
  ###### 
  
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
  simdata.final <- rbind(simdata.extended, simdata.no.trt)
  
  return(simdata.final)
}
