# We developed the function below to simulate an observational dataset  
# Using simulated data offers key advantages: 
# it avoids data security concerns and allows us to know the true treatment effects, 
# enabling us to evaluate the performance of different analytic strategies, including TTE and causal inference methods. 
# The code includes two main functions:
#   - generate.obs.data(): simulate person-month longitudinal data with a binary treatment indicator A where 1 means under active treatment and 0 means under control treatment.
#   - generate.2trt.obs.data(): simulate person-month longitudinal data with a treatment indicator A where 1 means under one treatment, 0 means under another treatment and NA means no treatment.
#                               Some subjects initiate one treatment then discontinue or switch to the other treatment at some time point, and some subjects never initiate either treatment.


library(dplyr)
library(lubridate)

#' Create a one-treatment observational dataset. A subject may under treatment (A==1) or control (A==0).
#'
#' This function simulates longitudinal person-month records for `n` subjects up to `J` months.
#' The simulated data includes baseline covariates, time-varying covariates, treatment A (binary), outcome Y, and censoring C.
#' Once an outcome (Y==1) or censoring (C==1) occurs, rows after that are dropped.
#'
#' @param n integer: number of subjects
#' @param J integer: number of months to simulate per subject (time = 0:(J-1))
#' @param beta0.trt0 numeric: intercept in the treatment model at baseline (logit scale)
#' @param beta.adherence numeric: impact of the received treatment at j-1 on receiving treatment at j (logit scale)
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
  mean.age <- 50
  sd.age <- 12
  
  baseline <- data.frame(id = seq_len(n), 
                         X1_0 = rbinom(n, 1, 0.5),            # baseline value of a time-varying binary covariate (e.g., medication use)
                         X2_0 = rnorm(n, 0, 1),               # baseline value of a time-varying continuous covariate (e.g., SBP)
                         age_0 = rnorm(n, mean.age, sd.age),  # baseline value of a time-varying age
                         X3 = rbinom(n, 1, 0.5),              # baseline binary covariate (e.g., sex)
                         X4 = rnorm(n, 0, 1)) %>%             # baseline continuous covariate (e.g., BMI)
    mutate(A0 = rbinom(n(), 1, plogis(beta0.trt0+0.5 * X1_0 + 0.5 * X2_0 - 0.2 * X3 +
                                        X4 - 0.3 * (age_0 - mean.age) / sd.age)))           # baseline treatment assignment 
  
  # Expand to person-month long format
  long_data <- baseline %>%
    slice(rep(1:n(), each = J)) %>%
    group_by(id) %>%
    arrange(id, .by_group = TRUE) %>%     # ensure rows are ordered
    mutate(time = row_number() - 1) %>%   # safe creation of time index 0,...,J-1
    mutate(age = age_0 + time / 12,                     # age in years; time unit = months
           ages = (age - mean.age) / sd.age,            # standardized age
      
           # initialize history columns
           A =  NA_integer_,
           Y =  NA_integer_,
           C =  NA_integer_,
           X1 = NA_real_,
           X2 =  NA_real_) %>%
    ungroup()
  
  # Simulate forward month-by-month within each id
  long_data <- long_data %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x
      for (j in seq_len(nrow(df))) { 
        if (j == 1) {
          # j==1 are baseline rows 
          df$C[j] <- 0L
          df$Y[j] <- 0L
          df$X1[j] <- df$X1_0[j]
          df$X2[j] <- df$X2_0[j]
          next
        }
        if (j == 2) {
          # j==2 start follow up
        df$A[j] <- df$A0[j]
        df$C[j] <- 0L
        df$Y[j] <- 0L
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j]+2*(df$X1[j-1]-0.5)))  # binary covariate, the coefficient before (df$X1[j-1]-0.5) determins ho
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j]+0.5*df$X2[j-1], sd = 1) 
          next
        }

        # Censoring: Censoring depends on previous treatment and previous covariates. once C==1 it stays 1
        if (df$C[j-1]== 0) {
          logit_q <- beta0.censoring - df$A[j-1] - 0.5 * df$X1[j-1] +
            0.5 * df$X2[j-1] - 0.2 * df$X3[j-1] + 0.2 * df$X4[j-1] - df$ages[j-1]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } else {
          df$C[j] <- 1L
        }
        
        # Treatment allocation at month j (depends on previous A and previous covariates)
        logit_pi <- beta.adherence * (df$A[j-1]- 0.5) +
          0.5 * df$X1[j-1] + 0.5 * df$X2[j-1] - 0.2 * df$X3[j-1] +
          df$X4[j-1] - 0.3 * df$ages[j-1]
        df$A[j] <- rbinom(1, 1, plogis(logit_pi))
        
        # Outcome: depends on current A and previous covariates.  once Y==1 it stays 1
        if (df$Y[j-1]== 0) {
          logit_lambda <- beta0.outcome + beta.trt.effect * df$A[j] +
            0.5 * df$X1[j-1] + 0.5 * df$X2[j-1] + df$X3[j-1] + df$X4[j-1] + 0.5 * df$ages[j-1]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        } else {
          df$Y[j] <- 1L
        }
        
        # Time-varying covariates: depend on previous treatment   
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j]+2*(df$X1[j-1]-0.5)))                 # binary covariate, the coefficient before (df$X1[j-1]-0.5) determines how much the current value depends on the previous value
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j]+0.5*df$X2[j-1], sd = 1)        # continuous covariate. This is problematic. Crystal, I need to discuss with you about this.  

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
  
  # make sure the variables are available in the temporal order of Censoring -> /*Adherence ->*/ Death -> outcome of interest -> time-varying covariates  
  long_data$A[long_data$C==1] <- NA
  long_data$Y[long_data$C==1] <- NA
  long_data$X1[long_data$C==1] <- NA
  long_data$X2[long_data$C==1] <- NA
  long_data$X1[long_data$Y==1] <- NA
  long_data$X2[long_data$Y==1] <- NA

  return(long_data)
}


# Add a pre-treatment period (random length >= 36 months) for each subject.
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

#' Create a two-treatment observational dataset (derived from generate.obs.data)
#'
#' This function:
#' 1) reinterprets A in generate.obs.data as one of two active treatments (A==0 vs A==1)
#' 2) randomly selects 20% subjects and sets their A = NA after deviation to simulate scenarios of medication discontinuation
#' 3) adds a pre-treatment period (36-60 months) during which neither of the study treatment is used 
#' 4) adds a set of subjects who never initiate either of the study treatment
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
  
  # 3) Add a pre-treatment period (random length >= 36 months) for each subject.
  #    In the pre-period: A = NA, Y = 0, C = 0. Covariates X1/X2 randomized (since not used later).
  
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
  # Combine all data to make the data to be more like the real clinical data:some some subjects initial one of the two treatments, some subjects never initiate either treatment. 
  # And those initial a treatment can discontinue or switch to the other treatment at some time point.
  simdata.final <- bind_rows(simdata.extended, simdata.no.trt)
  
  # Return a data.frame sorted by id and time
  simdata.final <- simdata.final %>% arrange(id, time)
  
  return(simdata.final)
}

#Now we use the function to generate 2 sample data set: 
# 1)Placebo-controlled design. Patients either receive active treatment or control treatment
 simdata <- generate.obs.data(n = 10000, J = 60,
                          beta0.trt0 = 0, beta.adherence = 10,
                          beta0.outcome = -5, beta.trt.effect = -1.2,
                          beta0.censoring = -3, seed = 123)
 
 # Apply extension per subject  
 simdata.extended <- simdata %>%
   mutate(.id_tmp = id) %>%
   group_by(.id_tmp) %>%
   group_modify(~ extend.data(.x)) %>%
   ungroup() %>%
   select(-.id_tmp)
 simdata.extended <- simdata.extended %>%mutate(A=ifelse(!is.na(A) ,A,0))%>% # if A is NA, treat it as no treatment
   filter(time>=-24) # only keep data within 2 years before index date
 
 # -------------------------
 # Decoration: to make the data more like the real clinical data we map time index to real calendar dates and set variable labels
 # -------------------------
 # 1) map time 0 for each subject to a random calendar month between 2000-01-01 and 2005-12-01:
 
 time0_dates <- sample(seq(as.Date("2000-01-01"), as.Date("2005-12-01"), by = "month"),
                       length(unique(simdata.extended$id )), replace = TRUE)
 time0_df <- data.frame(id = unique(simdata.extended$id ), time0 = time0_dates)
 obsdata  <- simdata.extended %>%
   left_join(time0_df, by = "id") %>%
   mutate(time = time0 %m+%  months(time)) %>%
   select(-time0)
 
 # 2) Set simple labels as attributes for each variable
 attr(obsdata$id, "label")    <- "Patient ID"
 attr(obsdata$time, "label")  <- "Time index for longitudinal records (months)"
 attr(obsdata$X1, "label")    <- "Non-ACEI or ARB antihypertensive medication use over time"
 attr(obsdata$X2, "label")    <- "Standardized systolic blood pressure over time"
 attr(obsdata$X3, "label")    <- "Biological sex (M=1, F=0)"
 attr(obsdata$X4, "label")    <- "Standardized diastolic blood pressure at baseline"
 attr(obsdata$age, "label")   <- "Age over time (years)"
 attr(obsdata$A, "label")     <- "Treatment indicator over time (ARB = 1, control = 0 )"
 attr(obsdata$Y, "label")     <- "Event indicator of cardiovascular disease"
 attr(obsdata$C, "label")     <- "Indicator of early dropout / censoring"
 saveRDS(obsdata, "obsdata1trt.rds")
 
 
# 2)Two-treatment design. Patients either receive one of the two active treatments or neither treatment
simdata <- generate.2trt.obs.data(n=10000, J=36, beta0.trt0 = 0, beta.adherence = 6,
                                  beta0.outcome = -5, beta.trt.effect = -1.2,
                                  beta0.censoring = -1, seed = 123)
# -------------------------
# Decoration: to make the data more like the real clinical data we map time index to real calendar dates and set variable labels
# -------------------------
# 1) map time 0 for each subject to a random calendar month between 2000-01-01 and 2005-12-01:

time0_dates <- sample(seq(as.Date("2000-01-01"), as.Date("2005-12-01"), by = "month"),
                      length(unique(simdata$id)), replace = TRUE)
time0_df <- data.frame(id = unique(simdata$id), time0 = time0_dates)
obsdata  <- simdata %>%
  left_join(time0_df, by = "id") %>%
  mutate(time = time0 %m+%  months(time)) %>%
  select(-time0)

# 2) Set simple labels as attributes for each variable
attr(obsdata$id, "label")    <- "Patient ID"
attr(obsdata$time, "label")  <- "Time index for longitudinal records (months)"
attr(obsdata$X1, "label")    <- "Non-ACEI or ARB antihypertensive medication use over time"
attr(obsdata$X2, "label")    <- "Standardized systolic blood pressure over time"
attr(obsdata$X3, "label")    <- "Biological sex (M=1, F=0)"
attr(obsdata$X4, "label")    <- "Standardized diastolic blood pressure at baseline"
attr(obsdata$age, "label")   <- "Age over time (years)"
attr(obsdata$A, "label")     <- "Treatment indicator over time (ARB = 1, ACEI = 0; NA = neither)"
attr(obsdata$Y, "label")     <- "Event indicator of cardiovascular disease"
attr(obsdata$C, "label")     <- "Indicator of early dropout / censoring"
saveRDS(obsdata, "obsdata2trt.rds")
