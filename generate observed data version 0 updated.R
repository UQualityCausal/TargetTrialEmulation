# simulate two-treatment longitudinal/clinical-record data
# Author: (revised)
# Purpose: Provide two functions:
#   - generate.obs.data(): simulate person-month longitudinal data with a single treatment indicator A
#   - generate.2trt.obs.data(): build on generate.obs.data() to create a data set where some subjects can initial one of two treatments and 
#                               deviate from the baseline treatment to the other treatment or no treatment. Some subjects never initial either of the two treatments.
#

library(dplyr)
library(lubridate)

#' Simulate person-month longitudinal data (single treatment indicator)
#'
#' This function simulates longitudinal person-month records for `n` subjects up to `J` months.
#' It simulates baseline covariates, time-varying covariates, treatment A (binary), outcome Y and censoring C.
#' Once an outcome (Y==1) or censoring (C==1) occurs, later rows for that subject are dropped.
#'
#' @param n integer: number of subjects
#' @param J integer: number of months to simulate per subject (time = 0:(J-1))
#' @param beta0.trt0 numeric: intercept in the treatment model at baseline (logit scale)
#' @param beta.adherence numeric: coefficient for adhereing the previous treatment in the treatment model (logit scale)
#' @param beta0.outcome numeric: intercept in the outcome model (logit scale)
#' @param beta.trt.effect numeric: treatment effect on outcome (logit scale)
#' @param beta0.censoring numeric: intercept in the censoring model (logit scale)
#' @param seed optional integer: if provided, sets seed for reproducibility
#' @return data.frame in long (person-month) format with columns:
#'   id, time, X1, X2, X3, X4, age, A, Y, C
#' @examples
#' sim <- generate.obs.data(n = 100, J = 120,
#'                          beta0.trt0 = 0, beta.adherence = 10,
#'                          beta0.outcome = -5, beta.trt.effect = -1.2,
#'                          beta0.censoring = -3, seed = 123)
generate.obs.data <- function(n, J,
                              beta0.trt0, beta.adherence,
                              beta0.outcome, beta.trt.effect,
                              beta0.censoring,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Baseline demographics: vectorized creation
  mean.age <- 35
  sd.age <- 12
  
  baseline <- data.frame(
    id = seq_len(n),
    X1_0 = rbinom(n, 1, 0.5),
    X2_0 = rnorm(n, 0, 1),
    age_0 = rnorm(n, mean.age, sd.age),
    X3 = rbinom(n, 1, 0.5),         # binary baseline covariate (e.g., sex)
    X4 = rnorm(n, 0, 1)            # baseline continuous covariate
  ) %>%
    mutate(
      # baseline treatment assignment (vectorized)
      A0 = rbinom(n(), 1, plogis(beta0.trt0))
    )
  
  # Expand to person-month long format
  long_data <- baseline %>%
    slice(rep(1:n(), each = J)) %>%
    group_by(id) %>%
    arrange(id, .by_group = TRUE) %>%     # ensure rows are ordered
    mutate(time = row_number() - 1) %>%   # safe creation of time index 0,...,J-1
    ungroup() %>%
    mutate(
      age = age_0 + time / 12,                     # age in years; time unit = months
      ages = (age - mean.age) / sd.age,            # standardized age
      # initialize history columns
      A = if_else(time == 0, A0, NA_integer_),
      Y = if_else(time == 0, 0L, NA_integer_),
      C = if_else(time == 0, 0L, NA_integer_),
      X1 = if_else(time == 0, X1_0, NA_real_),
      X2 = if_else(time == 0, X2_0, NA_real_),
      A_prev = NA_integer_,
      Y_prev = NA_integer_,
      C_prev = NA_integer_
    ) %>%
    ungroup()
  
  # Simulate forward month-by-month within each id
  long_data <- long_data %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x
      for (j in seq_len(nrow(df))) {
        if (j == 1) {
          # j==1 are baseline rows already initialized
          next
        }
        # set history
        df$A_prev[j] <- df$A[j - 1]
        df$Y_prev[j] <- df$Y[j - 1]
        df$C_prev[j] <- df$C[j - 1]
        
        # time-varying covariates: depend on previous treatment
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j - 1]))                 # binary covariate
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j - 1], sd = 1)       # continuous covariate
        
        # Treatment allocation at month j (depends on previous A and covariates)
        logit_pi <- beta.adherence * (df$A_prev[j] - 0.5) +
          0.5 * df$X1[j] + 0.5 * df$X2[j] - 0.2 * df$X3[j] +
          df$X4[j] - 0.3 * df$ages[j]
        df$A[j] <- rbinom(1, 1, plogis(logit_pi))
        
        # Outcome: once Y==1 it stays 1
        if (df$Y_prev[j] == 0) {
          logit_lambda <- beta0.outcome + beta.trt.effect * df$A[j] +
            0.5 * df$X1[j] + 0.5 * df$X2[j] + df$X3[j] + df$X4[j] + 0.5 * df$ages[j]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        } else {
          df$Y[j] <- 1L
        }
        
        # Censoring: once C==1 it stays 1
        if (df$C_prev[j] == 0) {
          logit_q <- beta0.censoring - df$A_prev[j] - 0.5 * df$X1[j] +
            0.5 * df$X2[j] - 0.2 * df$X3[j] + 0.2 * df$X4[j] - df$ages[j]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } else {
          df$C[j] <- 1L
        }
      }
      df
    }) %>%
    ungroup()
  
  # Drop records after first event (Y==1) or first censoring (C==1)
  first_event_time <- long_data %>%
    filter(C == 1 | Y == 1) %>%
    group_by(id) %>%
    summarise(first_time = min(time), .groups = "drop")
  
  long_data <- long_data %>%
    left_join(first_event_time, by = "id") %>%
    filter(is.na(first_time) | time <= first_time) %>%
    select(-first_time)
  
  # Ensure numeric types and select desired columns
  long_data <- long_data %>%
    mutate(across(c(id, time, A, Y, C, X1, X2, X3, X4, age), as.numeric)) %>%
    select(id, time, X1, X2, X3, X4, age, A, Y, C)
  
  return(long_data)
}


#' Create a two-treatment observational dataset (derived from generate.obs.data)
#'
#' This function:
#' 1) reinterprets A in generate.obs.data as one of two active treatments (A==0 vs A==1)
#' 2) randomly sets A = NA after deviation times to simulate periods of "neither treatment"
#' 3) prepends a randomized pre-treatment period (>= 36 months) during which A = NA, Y = 0, C = 0
#' 4) appends an extra set of subjects who never initiate treatment (A = NA)
#'
#' The goal is to produce realistic-looking clinical-record longitudinal data for
#' demonstrating the data-processing required to emulate a target trial.
#'
#' @param n integer: base number of subjects to simulate (per call to generate.obs.data)
#' @param J integer: months of post-index follow-up used in generate.obs.data
#' @param beta0.trt0 numeric
#' @param beta.adherence numeric
#' @param beta0.outcome numeric
#' @param beta.trt.effect numeric
#' @param beta0.censoring numeric
#' @param seed optional integer: if provided, sets seed (for reproducibility)
#' @return data.frame with extended pre-treatment months, treatment gaps (NA), and extra no-treatment subjects
#' @examples
#' sim2 <- generate.2trt.obs.data(n = 100, J = 120,
#'                                beta0.trt0 = 0, beta.adherence = 10,
#'                                beta0.outcome = -5, beta.trt.effect = -1.2,
#'                                beta0.censoring = -1, seed = 123)
generate.2trt.obs.data <- function(n, J,
                                   beta0.trt0, beta.adherence,
                                   beta0.outcome, beta.trt.effect,
                                   beta0.censoring,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # 1) Generate base data (A=0 vs A=1 represent two treatments)
  simdata <- generate.obs.data(n = n, J = J,
                               beta0.trt0 = beta0.trt0,
                               beta.adherence = beta.adherence,
                               beta0.outcome = beta0.outcome,
                               beta.trt.effect = beta.trt.effect,
                               beta0.censoring = beta0.censoring,
                               seed = NULL) # seed handled above
  
  # 2) Introduce deviations: after the subject first deviates from assigned baseline treatment,
  #    they have a 20% chance each month to take neither treatment (A = NA)
  assigned_trt <- simdata %>%
    filter(time == 0) %>%
    select(id, assigned_trt = A)
  
  simdata <- simdata %>%
    left_join(assigned_trt, by = "id") %>%
    group_by(id) %>%
    mutate(
      deviate.time = ifelse(any(A != assigned_trt & !is.na(A)),
                            min(time[A != assigned_trt & !is.na(A)]),
                            Inf)
    ) %>%
    ungroup() %>%
    mutate(
      A = ifelse(time >= deviate.time & runif(n()) < 0.20, NA, A)   # 20% chance to set to NA post-deviation
    ) %>%
    select(-assigned_trt, -deviate.time)
  
  # 3) Prepend a pre-treatment period (random length >= 36 months) for each subject.
  #    In the pre-period: A = NA, Y = 0, C = 0. Covariates X1/X2 randomized (since not used later).
  extend.data <- function(data_row) {
    # months.prev sampled per subject (36 to 60)
    months.prev <- sample(36:60, 1)
    
    # create pre-period rows by sampling rows from the subject and updating fields
    # Note: sample with replacement to create plausible covariate variety in the pre-period
    data_prev <- data_row %>%
      slice(sample(x = seq_len(nrow(data_row)), size = months.prev, replace = TRUE)) %>%
      mutate(
        time = -seq(from = months.prev, to = 1),
        A = NA_real_,
        Y = 0,
        C = 0,
        X1 = rbinom(n(), 1, 0.5),
        X2 = rnorm(n(), 0, 1),
        # set age for pre-period by subtracting months
        age = data_row$age[1] - seq(from = months.prev, to = 1) / 12
      )
    
    combined <- bind_rows(data_prev, data_row) %>%
      arrange(id, time)
    
    return(combined)
  }
  
  # Apply extension per subject
  simdata.extended <- simdata %>%
    mutate(.id_tmp = id) %>%
    group_by(.id_tmp) %>%
    group_modify(~ extend.data(.x)) %>%
    ungroup() %>%
    select(-.id_tmp)
  
  # 4) Add subjects who never initiate any treatment (A = NA for all months).
  simdata.no.trt <- generate.obs.data(n = n, J = J,
                                      beta0.trt0 = beta0.trt0,
                                      beta.adherence = beta.adherence,
                                      beta0.outcome = beta0.outcome,
                                      beta.trt.effect = beta.trt.effect,
                                      beta0.censoring = beta0.censoring,
                                      seed = NULL) %>%
    mutate(A = NA_real_) %>%
    arrange(id, time) %>%
    mutate(id = id + max(simdata.extended$id))  # shift ids so they don't overlap
  
  simdata.final <- bind_rows(simdata.extended, simdata.no.trt)
  
  # Return a data.frame sorted by id and time
  simdata.final <- simdata.final %>% arrange(id, time)
  
  return(simdata.final)
}


# -------------------------
# Example usage / checks
# -------------------------
# small quick-check (uncomment to run)
sim1 <- generate.obs.data(n = 100, J = 120, beta0.trt0 = 0, beta.adherence = 10,
                          beta0.outcome = -5, beta.trt.effect = -1.2,
                          beta0.censoring = -3, seed = 123)

sim2 <- generate.2trt.obs.data(n = 100, J = 120, beta0.trt0 = 0, beta.adherence = 6,
                               beta0.outcome = -5, beta.trt.effect = -1.2,
                               beta0.censoring = -1, seed = 123)
#
# head(sim2)

# -------------------------
# Decoration: map time index to real calendar dates and set variable labels
# -------------------------
# If you want to map time 0 for each subject to a random calendar month between 2000-01-01 and 2005-12-01:
# (uncomment and run after generating `simdata`)
simdata <- generate.2trt.obs.data(n = 100, J = 120, beta0.trt0 = 0, beta.adherence = 6,
                                  beta0.outcome = -5, beta.trt.effect = -1.2,
                                  beta0.censoring = -1, seed = 123)
time0_dates <- sample(seq(as.Date("2000-01-01"), as.Date("2005-12-01"), by = "month"),
                      length(unique(simdata$id)), replace = TRUE)
time0_df <- data.frame(id = unique(simdata$id), time0 = time0_dates)
obsdata <- simdata %>%
  left_join(time0_df, by = "id") %>%
  mutate(time = time0 %m+%  months(time)) %>%
  select(-time0)

# Set simple labels as attributes (use lowercase 'age')
attr(obsdata$id, "label")    <- "Patient ID"
attr(obsdata$time, "label")  <- "Time index for longitudinal records (months)"
attr(obsdata$X1, "label")    <- "Non-ACEI or ARB antihypertensive medication use over time"
attr(obsdata$X2, "label")    <- "Standardized systolic blood pressure over time"
attr(obsdata$X3, "label")    <- "Biological sex (binary indicator)"
attr(obsdata$X4, "label")    <- "Standardized diastolic blood pressure at baseline"
attr(obsdata$age, "label")   <- "Age over time (years)"
attr(obsdata$A, "label")     <- "Treatment indicator over time (ARB = 1, ACEI = 0; NA = neither)"
attr(obsdata$Y, "label")     <- "Event indicator of cardiovascular disease"
attr(obsdata$C, "label")     <- "Indicator of early dropout / censoring"
