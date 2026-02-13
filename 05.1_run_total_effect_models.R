# Purpose: Run the MR component of Graph MR to derive the total causal effects θ₁.

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages
library("tidyverse")
library("vroom")
library("parallel")

# How many cores do you need?
ncores <- 6

# Start array indexing
task_id_string <- Sys.getenv("PBS_ARRAY_INDEX")
k <- as.numeric(task_id_string)

# Arguments: method_type, dataset (any label, e.g. proteomics, mydata)
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
dataset     <- if (length(args) >= 2L) args[2] else Sys.getenv("GRAPH_MR_DATASET")
if (!nzchar(dataset)) stop("Provide dataset name as second argument or set GRAPH_MR_DATASET.")
base_dir <- paste0(method_type, "_", dataset)

# Read preprocessed exposures (generic layout)
exposures_path <- paste0(data_dir, base_dir, "/preprocessed/exposures.tsv")
if (!file.exists(exposures_path)) stop("Preprocessed exposures not found: ", exposures_path)
input_df <- vroom(exposures_path)

# Define the pattern based on method_type
pattern <- switch(
  method_type,
  cftlmr = "cftlmr_scores\\.csv",
  cfmr   = "cfmr_scores\\.csv",
  stop("Unknown method_type: ", method_type)
)

# Use the pattern in list.files
input_s_list <- list.files(
  path = paste0(results_dir, base_dir, "/genetic_scores"),
  pattern = pattern,
  full.names = TRUE,
  recursive = TRUE
)

# Dedicate one PBS array element to each element in the input_s_list (list of genetic scores S) 
input_s_list <- input_s_list[k]

# Read age and sex for MR (generic: from pipeline layout)
covariates_path <- paste0(data_dir, base_dir, "/gwas_input/covariates_age_sex.txt")
if (!file.exists(covariates_path)) stop("Covariates not found: ", covariates_path, ". Run preprocessing (see docs/INPUT_SPEC.md).")
age_sex <- vroom(covariates_path) %>%
  dplyr::select(barcode, AGE, SEX)
  
# Set the output path for the total effects
output <- paste0(results_dir, base_dir, "/total_fx/")
dir.create(output, showWarnings = F, recursive = T)



# Start script ----
writeLines(c(" ", " ", "Started script 05.1 - run total fx MR.", " ", " "))

# Run MR for each predicted exposure (X̂) and outcome (Y) pair in turn 
# Estimates total causal effect θ₁ of predicted exposure X̂ on outcome Y using two-sample 2SLS
run_total_effects_mr <- function(input_s_list) {
  mclapply(input_s_list, mc.cores = ncores, FUN = function(s_i){
    # Extract the exposure name string 
    exposure_base_name <- basename(dirname(s_i))
    exposure_s <- paste0(exposure_base_name, "_scores")  # Variable name for exposure genetic score S
    # Read exposure genetic scores S and rename the S column
    s_i <- vroom(s_i, delim = ",", show_col_types = F) %>% 
      rename(!!exposure_s := gs)
    # Inner join the fold-specific exposure genetic score S and all phenotype levels
    regression_table <- inner_join(s_i, input_df, by = "barcode") %>% 
      inner_join(age_sex, by = "barcode") %>%
      as.data.frame() %>% 
      mutate(barcode = as.character(barcode))  %>%
      distinct(barcode, .keep_all = TRUE) %>% 
      filter(!grepl("-9", barcode)) %>%
      column_to_rownames("barcode") %>% 
      rename_with(~ gsub(" ", "_", .))
    
    # Identify the exposure phenotype column X that matches the genetic score S
    # The exposure phenotype X should have the same name as exposure_base_name
    exposure_phenotype <- NULL
    if (exposure_base_name %in% colnames(regression_table)) {
      exposure_phenotype <- exposure_base_name
    }     
    if (is.null(exposure_phenotype)) {
      warning(paste("Could not find exposure phenotype for", exposure_base_name))
    }
    # Total effects estimation: Two-sample 2SLS to estimate θ₁ (total causal effect of X on Y)
    # First stage: X ~ S (estimate β_X, the S--X association, only where X is observed)
    # Second stage: Y ~ X̂ (estimate θ₁, where X̂ is predicted from S for all with Y and S)
    regress_yvars_on_xs <- function(regression_table, exposure_phenotype, exposure_s, adjust_for_age_sex = FALSE) {
      # Create a vector of variable names excluding the genetic score S and covariates (AGE, SEX)
      # We want ALL variables as outcomes Y, including the exposure phenotype X itself (cartesian product)
      y_vars <- setdiff(colnames(regression_table), c(exposure_s, "AGE", "SEX"))  # All variables except S and covariates
      
      # Helper function to create an NA-filled result row for a given outcome
      create_na_result <- function(y_i) {
        data.frame(Response = y_i,
                   Exposure = ifelse(is.null(exposure_phenotype), NA_character_, exposure_phenotype),
                   Causal_estimate = NA_real_,
                   P_value = NA_real_,
                   F_stat = NA_real_,
                   Model = NA_character_,
                   N_first_stage = NA_integer_,
                   N_second_stage = NA_integer_,
                   R_squared = NA_real_,
                   stringsAsFactors = FALSE)
      }
      
      # If exposure phenotype X not found, return NA-filled rows for all outcomes
      if (is.null(exposure_phenotype) || !exposure_phenotype %in% colnames(regression_table)) {
        warning(paste("Exposure phenotype X", exposure_phenotype, "not found in regression table. Returning NA-filled results."))
        return(do.call(rbind, lapply(y_vars, create_na_result)))
      }
      
      # Apply two-sample 2SLS to each outcome variable Y
      iv_results <- lapply(y_vars, function(y_i) {
        yvec <- regression_table[[y_i]]  # Outcome Y
        xvec <- regression_table[[exposure_phenotype]]  # Exposure X
        svec <- regression_table[[exposure_s]]  # Genetic score S
        
        tryCatch({
          # Two-sample 2SLS
          # First stage: X ~ S (or X ~ S + AGE + SEX if adjusting)
          if (adjust_for_age_sex) {
            first_stage_complete <- !is.na(xvec) & !is.na(svec) & !is.na(regression_table$AGE) & !is.na(regression_table$SEX)
            first_stage_formula <- paste0(exposure_phenotype, " ~ ", exposure_s, " + AGE + SEX")
          } else {
            first_stage_complete <- !is.na(xvec) & !is.na(svec)
            first_stage_formula <- paste0(exposure_phenotype, " ~ ", exposure_s)
          }
          n1 <- sum(first_stage_complete)
          
          # Check if we have enough data for first stage
          if (n1 < 3) {
            stop("Insufficient data for first stage regression (n < 3)")
          }
          
          first_stage_lm <- lm(first_stage_formula, data = regression_table[first_stage_complete, ])
          
          # Safely extract F-statistic
          # For unadjusted models: overall F-statistic (tests S)
          # For adjusted models: conditional F-statistic for S given AGE and SEX (using t²)
          first_stage_f_stat <- tryCatch({
            if (adjust_for_age_sex) {
              # Get conditional F-statistic for S using t-statistic from full model
              # F = t² for single coefficient, and t-statistic from full model already conditions on AGE and SEX
              coef_summary <- summary(first_stage_lm)$coefficients
              if (!is.null(coef_summary) && exposure_s %in% rownames(coef_summary)) {
                t_stat <- coef_summary[exposure_s, "t value"]
                if (!is.na(t_stat) && is.finite(t_stat)) {
                  t_stat^2  # F = t² for single coefficient
                } else {
                  NA
                }
              } else {
                NA
              }
            } else {
              # For unadjusted model, overall F-statistic is the F-statistic for S
              fstat <- summary(first_stage_lm)$fstatistic
              if (!is.null(fstat) && length(fstat) > 0 && is.finite(fstat[1])) {
                fstat[1]
              } else {
                NA
              }
            }
          }, error = function(e) {
            NA
          })
          
          # Check if first stage model is valid
          if (is.na(first_stage_f_stat) || !is.finite(first_stage_f_stat)) {
            stop("First stage model failed or F-statistic is invalid")
          }
          
          # Get first-stage coefficient β_X and its variance Var(β_X)
          beta_X <- tryCatch({
            coef(first_stage_lm)[exposure_s]
          }, error = function(e) {
            stop("Could not extract first stage coefficient")
          })
          
          if (is.na(beta_X) || !is.finite(beta_X)) {
            stop("First stage coefficient is NA or invalid")
          }
          
          first_stage_vcov <- tryCatch({
            vcov(first_stage_lm)
          }, error = function(e) {
            stop("Could not compute first stage variance-covariance matrix")
          })
          
          var_beta_X <- tryCatch({
            first_stage_vcov[exposure_s, exposure_s]
          }, error = function(e) {
            stop("Could not extract variance of first stage coefficient")
          })
          
          # Predict X̂ = S * β̂_X for all individuals who have Y and S (even if X is missing)
          # Second stage complete cases: need Y and S (and AGE/SEX if adjusting)
          if (adjust_for_age_sex) {
            second_stage_complete <- !is.na(yvec) & !is.na(svec) & !is.na(regression_table$AGE) & !is.na(regression_table$SEX)
          } else {
            second_stage_complete <- !is.na(yvec) & !is.na(svec)
          }
          n2 <- sum(second_stage_complete)
          
          # Check if we have enough data for second stage
          if (n2 < 3) {
            stop("Insufficient data for second stage regression (n < 3)")
          }
          
          xhat <- tryCatch({
            predict(first_stage_lm, newdata = regression_table[second_stage_complete, ])
          }, error = function(e) {
            stop("Could not predict X̂ from first stage model")
          })
          
          if (any(is.na(xhat)) || any(!is.finite(xhat))) {
            stop("Predicted X̂ contains NA or infinite values")
          }
          
          # Determine if Y is binary (using second stage complete cases)
          unique_vals <- unique(na.omit(yvec[second_stage_complete]))
          is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
          
          if (is_binary) {
            # Second stage: logistic regression Y ~ X̂ (or Y ~ X̂ + AGE + SEX if adjusting)
            # Note: Standard errors from GLM do not fully account for first-stage uncertainty
            # in X̂, but this is standard practice in MR with binary outcomes
            
            second_stage_data <- data.frame(
              y = yvec[second_stage_complete],
              xhat = xhat
            )
            if (adjust_for_age_sex) {
              second_stage_data$AGE <- regression_table$AGE[second_stage_complete]
              second_stage_data$SEX <- regression_table$SEX[second_stage_complete]
            }
            
            # Capture warnings from glm() and check for convergence issues
            glm_warnings <- NULL
            glm_model <- NULL
            model_failed <- FALSE
            
            # Fit the model while capturing and suppressing warnings
            withCallingHandlers({
              if (adjust_for_age_sex) {
                glm_model <- glm(y ~ xhat + AGE + SEX, data = second_stage_data, family = "binomial")
              } else {
                glm_model <- glm(y ~ xhat, data = second_stage_data, family = "binomial")
              }
            }, warning = function(w) {
              glm_warnings <<- c(glm_warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            })
            
            # Check if model converged and if problematic warnings occurred
            if (!is.null(glm_model)) {
              # Check convergence
              if (!glm_model$converged) {
                model_failed <- TRUE
              }
              # Check for specific problematic warnings
              if (!is.null(glm_warnings)) {
                problematic_warnings <- grepl("fitted probabilities numerically 0 or 1 occurred", glm_warnings, ignore.case = TRUE) |
                                       grepl("algorithm did not converge", glm_warnings, ignore.case = TRUE)
                if (any(problematic_warnings)) {
                  model_failed <- TRUE
                }
              }
            } else {
              model_failed <- TRUE
            }
            
            # If model failed due to warnings or non-convergence, return NA-filled results
            if (model_failed) {
              result <- data.frame(Response = y_i,  # Outcome Y
                                   Exposure = exposure_phenotype,  # Exposure X
                                   Causal_estimate = NA_real_,  # θ̂₁: NA due to convergence/separation issues
                                   P_value = NA_real_,
                                   F_stat = first_stage_f_stat,  # First-stage F-statistic for S (conditional on covariates if adjusted)
                                   Model = "2S2SLS_logistic",
                                   N_first_stage = n1,  # Sample size for first stage (X and S observed)
                                   N_second_stage = n2,  # Sample size for second stage (Y and S observed)
                                   R_squared = NA_real_)  # McFadden's pseudo-R² (partial if adjusted)
              return(result)
            }
            
            # If model is valid, proceed with extracting results
            theta_1_hat <- coef(glm_model)["xhat"]  # θ̂₁: total causal effect (log-odds ratio)
            
            # Extract standard error and p-value from GLM
            # Note: This variance does not account for uncertainty in first-stage estimation
            # of β_X, which is a known limitation of two-stage logistic regression
            glm_summary <- summary(glm_model)
            se_glm <- glm_summary$coefficients["xhat", "Std. Error"]
            p_glm <- glm_summary$coefficients["xhat", "Pr(>|z|)"]
            
            # McFadden's pseudo R-squared
            # For adjusted models: compute partial pseudo R² for xhat (using covariate-only model as baseline)
            # For unadjusted models: compute standard pseudo R² (using intercept-only model as baseline)
            mcfadden_r2 <- tryCatch({
              if (adjust_for_age_sex) {
                # Partial McFadden R²: baseline is Y ~ AGE + SEX (not intercept-only)
                # This isolates the contribution of xhat after accounting for covariates
                baseline_model <- withCallingHandlers({
                  glm(y ~ AGE + SEX, data = second_stage_data, family = "binomial")
                }, warning = function(w) {
                  invokeRestart("muffleWarning")
                })
                # Partial McFadden R² = 1 - (log-likelihood of full model / log-likelihood of covariate-only model)
                1 - (logLik(glm_model)[1] / logLik(baseline_model)[1])
              } else {
                # Standard McFadden R²: baseline is intercept-only model
                null_model <- withCallingHandlers({
                  glm(y ~ 1, data = second_stage_data, family = "binomial")
                }, warning = function(w) {
                  invokeRestart("muffleWarning")
                })
                # McFadden R² = 1 - (log-likelihood of fitted model / log-likelihood of null model)
                1 - (logLik(glm_model)[1] / logLik(null_model)[1])
              }
            }, error = function(e) NA)
            
            # Build result data.frame with only relevant columns
            result <- data.frame(Response = y_i,  # Outcome Y
                                 Exposure = exposure_phenotype,  # Exposure X
                                 Causal_estimate = theta_1_hat,  # θ̂₁: total causal effect (log-odds ratio)
                                 P_value = p_glm,
                                 F_stat = first_stage_f_stat,  # First-stage F-statistic for S (conditional on covariates if adjusted)
                                 Model = "2S2SLS_logistic",
                                 N_first_stage = n1,  # Sample size for first stage (X and S observed)
                                 N_second_stage = n2,  # Sample size for second stage (Y and S observed)
                                 R_squared = mcfadden_r2)  # McFadden's pseudo-R² (partial if adjusted)
            
            return(result)
          } else {
            # Second stage: linear regression Y ~ X̂ (or Y ~ X̂ + AGE + SEX if adjusting)
            second_stage_data <- data.frame(
              y = yvec[second_stage_complete],
              xhat = xhat,
              s = svec[second_stage_complete]
            )
            if (adjust_for_age_sex) {
              second_stage_data$AGE <- regression_table$AGE[second_stage_complete]
              second_stage_data$SEX <- regression_table$SEX[second_stage_complete]
              second_stage_lm <- lm(y ~ xhat + AGE + SEX, data = second_stage_data)
            } else {
              second_stage_lm <- lm(y ~ xhat, data = second_stage_data)
            }
            theta_1_hat <- coef(second_stage_lm)["xhat"]  # θ̂₁: total causal effect of X on Y
            beta_all <- coef(second_stage_lm)  # All coefficients from second stage
            
            # Calculate variance using Inoue & Solon (2010) formula for two-sample 2SLS
            # This accounts for uncertainty in both first and second stages
            
            # Step 1: Build X1_hat matrix (predicted X for second stage sample)
            if (adjust_for_age_sex) {
              X1_hat <- model.matrix(~ xhat + AGE + SEX, data = second_stage_data)
            } else {
              X1_hat <- model.matrix(~ xhat, data = second_stage_data)
            }
            y1 <- yvec[second_stage_complete]
            
            # Step 2: Calculate sigma_2 (second stage residual variance)
            k_p <- ncol(X1_hat)  # Number of parameters
            pred_y1 <- X1_hat %*% beta_all
            eps <- y1 - pred_y1
            sigma_2 <- as.numeric((t(eps) %*% eps) / (n2 - k_p))
            
            # Step 3: Build X2 matrix (actual X for first stage sample)
            if (adjust_for_age_sex) {
              X2 <- model.matrix(~ xvec[first_stage_complete] + 
                                 regression_table$AGE[first_stage_complete] + 
                                 regression_table$SEX[first_stage_complete])
            } else {
              X2 <- model.matrix(~ xvec[first_stage_complete])
            }
            
            # Step 4: Build pred_X2 (X2 with X replaced by predicted values from first stage)
            pred_X2 <- X2
            pred_x_first_stage <- predict(first_stage_lm, newdata = regression_table[first_stage_complete, ])
            pred_X2[, 2] <- pred_x_first_stage  # Replace X column (column 2) with predicted values
            
            # Step 5: Calculate sigma_nu (first stage residual variance matrix)
            if (adjust_for_age_sex) {
              Z2 <- model.matrix(~ svec[first_stage_complete] + 
                                 regression_table$AGE[first_stage_complete] + 
                                 regression_table$SEX[first_stage_complete])
            } else {
              Z2 <- model.matrix(~ svec[first_stage_complete])
            }
            k_q <- ncol(Z2)  # Number of instruments
            eps_1s <- X2 - pred_X2
            sigma_nu <- (t(eps_1s) %*% eps_1s) / (n1 - k_q)
            
            # Step 6: Calculate sigma_f = sigma_2 + (n2/n1) * beta' * sigma_nu * beta
            beta_vec <- matrix(beta_all, ncol = 1)  # Convert to column vector
            sigma_f <- sigma_2 + (n2 / n1) * as.numeric(t(beta_vec) %*% sigma_nu %*% beta_vec)
            
            # Step 7: Calculate Var = sigma_f * (X1_hat' * X1_hat)^(-1)
            X1_hat_t_X1_hat <- t(X1_hat) %*% X1_hat
            Var_beta_ts2sls <- sigma_f * solve(X1_hat_t_X1_hat)
            
            # Step 8: Extract SE for xhat coefficient (row/column 2, since column 1 is intercept)
            two_sample_se <- sqrt(Var_beta_ts2sls[2, 2])
            
            # Calculate p-value using t-distribution
            t_stat <- theta_1_hat / two_sample_se
            robust_p <- 2 * pt(abs(t_stat), df = n2 - k_p, lower.tail = FALSE)
            
            # Build result data.frame with only relevant columns
            # For adjusted models: compute partial R² for xhat (isolates genetic score's contribution)
            # For unadjusted models: compute standard R²
            r_squared <- tryCatch({
              if (adjust_for_age_sex) {
                # Partial R² for xhat after accounting for AGE and SEX
                # Partial R² = (R²_full - R²_reduced) / (1 - R²_reduced)
                reduced_model <- lm(y ~ AGE + SEX, data = second_stage_data)
                r2_full <- summary(second_stage_lm)$r.squared
                r2_reduced <- summary(reduced_model)$r.squared
                if (!is.na(r2_full) && !is.na(r2_reduced) && r2_reduced < 1) {
                  (r2_full - r2_reduced) / (1 - r2_reduced)
                } else {
                  NA
                }
              } else {
                # Standard R² for unadjusted model
                summary(second_stage_lm)$r.squared
              }
            }, error = function(e) NA)
            
            result <- data.frame(Response = y_i,  # Outcome Y
                                 Exposure = exposure_phenotype,  # Exposure X
                                 Causal_estimate = theta_1_hat,  # θ̂₁: total causal effect
                                 P_value = robust_p,
                                 F_stat = first_stage_f_stat,  # First-stage F-statistic for S (conditional on covariates if adjusted)
                                 Model = "2S2SLS",
                                 N_first_stage = n1,  # Sample size for first stage (X and S observed)
                                 N_second_stage = n2,  # Sample size for second stage (Y and S observed)
                                 R_squared = r_squared)  # R² for xhat (partial R² if adjusted)
            
            return(result)
          }
        }, error = function(e) {
          # Handle the error and return minimal data.frame
          # Try to compute sample sizes even if models fail
          if (adjust_for_age_sex) {
            first_stage_complete <- !is.na(xvec) & !is.na(svec) & !is.na(regression_table$AGE) & !is.na(regression_table$SEX)
            second_stage_complete <- !is.na(yvec) & !is.na(svec) & !is.na(regression_table$AGE) & !is.na(regression_table$SEX)
            first_stage_formula_temp <- paste0(exposure_phenotype, " ~ ", exposure_s, " + AGE + SEX")
          } else {
            first_stage_complete <- !is.na(xvec) & !is.na(svec)
            second_stage_complete <- !is.na(yvec) & !is.na(svec)
            first_stage_formula_temp <- paste0(exposure_phenotype, " ~ ", exposure_s)
          }
          first_stage_f_stat <- NA_real_
          
          # Try to get first stage F-stat if possible
          if (sum(first_stage_complete) >= 3) {
            tryCatch({
              first_stage_lm_temp <- lm(first_stage_formula_temp, 
                                        data = regression_table[first_stage_complete, ])
              if (adjust_for_age_sex) {
                # Get conditional F-statistic for S using t-statistic (F = t²)
                coef_summary_temp <- summary(first_stage_lm_temp)$coefficients
                if (!is.null(coef_summary_temp) && exposure_s %in% rownames(coef_summary_temp)) {
                  t_stat_temp <- coef_summary_temp[exposure_s, "t value"]
                  if (!is.na(t_stat_temp) && is.finite(t_stat_temp)) {
                    first_stage_f_stat <- t_stat_temp^2
                  }
                }
              } else {
                # For unadjusted model, overall F-statistic is the F-statistic for S
                fstat <- summary(first_stage_lm_temp)$fstatistic
                if (!is.null(fstat) && length(fstat) > 0 && is.finite(fstat[1])) {
                  first_stage_f_stat <- fstat[1]
                }
              }
            }, error = function(e2) {})
          }
          
          return(data.frame(Response = y_i,  # Outcome Y
                            Exposure = exposure_phenotype,  # Exposure X
                            Causal_estimate = NA_real_,  # θ̂₁: total causal effect (NA if model failed)
                            P_value = NA_real_,
                            F_stat = first_stage_f_stat,  # First-stage F-statistic for S (conditional on covariates if adjusted)
                            Model = NA_character_,
                            N_first_stage = sum(first_stage_complete),  # Sample size for first stage
                            N_second_stage = sum(second_stage_complete),  # Sample size for second stage
                            R_squared = NA_real_))  # R-squared (NA if model failed)
        })
      })
      
      # Filter out NULL or invalid results before binding
      iv_results <- iv_results[!sapply(iv_results, is.null)]
      iv_results <- iv_results[sapply(iv_results, function(x) is.data.frame(x) && nrow(x) > 0)]
      
      # If no valid results, return NA-filled rows for all outcomes
      if (length(iv_results) == 0) {
        return(do.call(rbind, lapply(y_vars, create_na_result)))
      }
      
      # Combine valid results
      combined_results <- do.call(rbind, iv_results)
      
      # Ensure we have a row for every outcome variable (fill missing ones with NA)
      missing_outcomes <- setdiff(y_vars, combined_results$Response)
      if (length(missing_outcomes) > 0) {
        na_rows <- do.call(rbind, lapply(missing_outcomes, create_na_result))
        combined_results <- rbind(combined_results, na_rows)
      }
      
      return(combined_results)
    }
    
    # Run unadjusted models
    regression_results_unadjusted <- regress_yvars_on_xs(regression_table, exposure_phenotype, exposure_s, adjust_for_age_sex = FALSE)
    
    # Run age/sex-adjusted models
    regression_results_adjusted <- regress_yvars_on_xs(regression_table, exposure_phenotype, exposure_s, adjust_for_age_sex = TRUE)
    
    # Always write the results files, even if all values are NA
    # This ensures we get one output file per exposure score for each model type
    vroom_write(regression_results_unadjusted, file = paste0(output, "/", exposure_s, "_total_fx.csv"), delim = ",", append = F)
    vroom_write(regression_results_adjusted, file = paste0(output, "/", exposure_s, "_total_fx_agesex.csv"), delim = ",", append = F)
    
    # Print a statement for debugging purposes
    n_valid_unadj <- sum(!is.na(regression_results_unadjusted$Causal_estimate))
    n_valid_adj <- sum(!is.na(regression_results_adjusted$Causal_estimate))
    n_total <- nrow(regression_results_unadjusted)
    cat(c("We conducted MR with this exposure genetic score S: ", exposure_s, 
          " (unadjusted: ", n_valid_unadj, "/", n_total, " valid results; adjusted: ", n_valid_adj, "/", n_total, " valid results)\n"))
  }
  )
}

total_effects_mr_results <- run_total_effects_mr(input_s_list)

writeLines(c(" ", " ", "Finished script 05.1 - finished running total fx MR", " ", " "))
