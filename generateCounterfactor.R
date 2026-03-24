#' Create a one-treatment observational dataset. A subject may under treatment (A==1) or control (A==0).
#'
#' This function simulates longitudinal person-month records for `n` subjects up to `J` months.
#' The simulated data includes baseline covariates, time-varying covariates, treatment A (binary), outcome Y, and censoring C.
#' Once an outcome (Y==1) or censoring (C==1) occurs, rows after that are dropped.
#'
#' @param J integer: number of months to simulate per subject (time = 0:(J-1))
#' @param beta.adherence numeric: impact of the received treatment at j-1 on receiving treatment at j (logit scale)
#' @param beta0.outcome numeric: intercept in the outcome model (logit scale)
#' @param beta.trt.effect numeric: treatment effect on outcome (logit scale)
#' @param beta0.censoring numeric: intercept in the censoring model (logit scale)
#' @param baseline data.frame: a single-row data frame with baseline characteristics for one subject.
#'   Should contain columns: id, X1_0, X2_0, age_0, X3, X4.
#' @param A0 numeric: baseline treatment assignment (0 or 1) for the subject
#' @param seed optional integer: if provided, sets seed for reproducibility
#' @return data.frame in long (person-month) format with columns:
#'   id, time, X1, X2, X3, X4, age, A, Y, C
#' @examples
#' baseline <- data.frame(id = 1, X1_0 = 1, X2_0 = 0.5, age_0 = 55, X3 = 0, X4 = -0.3)
#' sim <- generate.counterfactor(J = 120, baseline = baseline, A0 = 1,
#'   beta.adherence = 10,
#'   beta0.outcome = -5, beta.trt.effect = -1.2,
#'   beta0.censoring = -3, seed = 123)
library(dplyr)

# Suppress NSE notes for dplyr verbs in script-style code.
A <- Y <- C <- X1 <- X2 <- X3 <- X4 <- NULL
X1_0 <- X2_0 <- A0 <- NULL
age <- age_0 <- ages <- NULL
id <- time <- first_time <- NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "A", "Y", "C", "X1", "X2", "X3", "X4",
    "X1_0", "X2_0", "A0",
    "age", "age_0", "ages",
    "id", "time", "first_time"
  ))
}

generate.counterfactor <- function(J, baseline, A0,
                                   beta.adherence,
                                   beta0.outcome, beta.trt.effect,
                                   beta0.censoring,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Baseline demographics
  mean.age <- 50
  sd.age <- 12

  # Ensure baseline has exactly 1 row and add A0
  baseline <- baseline[1, , drop = FALSE] %>%
    mutate(A0 = A0)

  # Expand to person-month long format (for 1 subject)
  long_data <- baseline %>%
    slice(rep(1, each = J)) %>%
    group_by(id) %>%
    arrange(id, .by_group = TRUE) %>%     # Ensure rows are ordered
    mutate(time = row_number() - 1) %>%   # Creation of time index 0,...,J-1
    mutate(age = .data$age_0 + .data$time / 12, # Age in years; time unit = months
      ages = (.data$age - mean.age) / sd.age, # Standardized age

      # Initialize history columns
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
          # j==1 means baseline (time=0)
          df$C[j] <- 0L
          df$Y[j] <- 0L
          df$X1[j] <- df$X1_0[j]
          df$X2[j] <- df$X2_0[j]
          next
        }
        if (j == 2) {
          # j==2 means follow-up month 1 (time=1)
          df$A[j] <- df$A0[j]
          df$C[j] <- 0L
          df$Y[j] <- 0L
          df$X1[j] <- rbinom(1, 1, plogis(-df$A[j] + 2 * (df$X1[j - 1] - 0.5)))  
          df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j] + 0.5 * df$X2[j - 1], sd = 1)
          next
        }

        # Censoring: Censoring depends on previous treatment and previous covariates. once C==1 it stays 1
        if (df$C[j - 1] == 0) {
          logit_q <- beta0.censoring - df$A[j - 1] - 0.5 * df$X1[j - 1] +
            0.5 * df$X2[j - 1] - 0.2 * df$X3[j - 1] + 0.2 * df$X4[j - 1] - df$ages[j - 1]
          df$C[j] <- rbinom(1, 1, plogis(logit_q))
        } else {
          df$C[j] <- 1L
        }

        # Treatment allocation at month j (depends on previous A and previous covariates)
        logit_pi <- beta.adherence * (df$A[j - 1] - 0.5) +
          0.5 * df$X1[j - 1] + 0.5 * df$X2[j - 1] - 0.2 * df$X3[j - 1] +
          0.5 * df$X4[j - 1] - 0.3 * df$ages[j - 1]
        df$A[j] <- rbinom(1, 1, plogis(logit_pi))

        # Outcome: depends on current A and previous covariates.  once Y==1 it stays 1
        if (df$Y[j - 1] == 0) {
          logit_lambda <- beta0.outcome + beta.trt.effect * df$A[j] +
            0.5 * df$X1[j - 1] + 0.5 * df$X2[j - 1] + df$X3[j - 1] - 0.5 * df$X4[j - 1] + 0.5 * df$ages[j - 1]
          df$Y[j] <- rbinom(1, 1, plogis(logit_lambda))
        } else {
          df$Y[j] <- 1L
        }

        # Time-varying covariates: depend on previous treatment
        df$X1[j] <- rbinom(1, 1, plogis(-df$A[j] + 2 * (df$X1[j - 1] - 0.5))) # Binary covariate, the coefficient before (df$X1[j-1]-0.5) determines how much the current value depends on the previous value
        df$X2[j] <- rnorm(1, mean = -0.3 * df$A[j] + 0.5 * df$X2[j - 1], sd = 1) # Continuous covariate.  

      }
      df
    }) %>%
    ungroup()

  # Drop records after first event (Y==1) or first censoring (C==1)
  first_event_time <- long_data %>%
    filter(.data$C == 1 | .data$Y == 1) %>%
    group_by(id) %>%
    summarise(first_time = min(time), .groups = "drop")

  long_data <- long_data %>%
    left_join(first_event_time, by = "id") %>%
    filter(is.na(.data$first_time) | .data$time <= .data$first_time) %>%
    select(-any_of("first_time"))

  # Ensure numeric types and select desired columns
  numeric_cols <- c("id", "time", "A", "Y", "C", "X1", "X2", "X3", "X4", "age")
  output_cols <- c("id", "time", "X1", "X2", "X3", "X4", "age", "A", "Y", "C")
  long_data <- long_data %>%
    mutate(across(all_of(numeric_cols), as.numeric)) %>%
    select(all_of(output_cols))

  # Make sure the variables are available in the temporal order of Censoring -> Adherence -> Death -> Outcome of interest -> Time-varying covariates
  long_data$A[long_data$C == 1] <- NA
  long_data$Y[long_data$C == 1] <- NA
  long_data$X1[long_data$C == 1] <- NA
  long_data$X2[long_data$C == 1] <- NA
  long_data$X1[long_data$Y == 1] <- NA
  long_data$X2[long_data$Y == 1] <- NA

  return(long_data)
}
