library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(gtsummary)


cancer_perception_df <- read.csv("CNA Catchment Area Data - 06.26.26.csv") %>% 
  dplyr::select(X, Cancerbeliefs_1, Cancerbeliefs_2, Cancerbeliefs_4) %>%
  mutate(
    across(
      c(Cancerbeliefs_1, Cancerbeliefs_2, Cancerbeliefs_4),
      ~ case_when(
        . == "Strongly disagree" ~ 0,
        . == "Disagree" ~ 1,
        . == "Neither agree nor disagree" ~ 2,
        . == "Agree" ~ 3,
        . == "Strongly agree" ~ 4,
        TRUE ~ NA_real_
      )
    )
  )

cancer_dem_df <- read.csv("Clean and Recategorized CNA Demographics Data - 06.11.26 .csv") %>%
  dplyr::select(X, Age_Group, Education, Race.Ethnicity, Sex, Income, Insurance.Coverage, Employment, Marital.Status, MedHist_cancerdx, Area_Type)

cna_df <- merge(cancer_perception_df, cancer_dem_df, by = 'X')

tbl_summary(cna_df)

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)

recode_cancer_belief <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  case_when(
    x == "Strongly disagree" ~ 0,
    x == "Disagree" ~ 1,
    x == "Neither agree nor disagree" ~ 2,
    x == "Agree" ~ 3,
    x == "Strongly agree" ~ 4,
    TRUE ~ NA_real_
  )
}

cna_df <- cna_df %>%
  mutate(
    across(
      c(Cancerbeliefs_1, Cancerbeliefs_2, Cancerbeliefs_4),
      recode_cancer_belief
    )
  )

cna <- cna_df %>%
  select(
    X,
    Cancerbeliefs_1, Cancerbeliefs_2, Cancerbeliefs_4,
    Area_Type,
    Age_Group, Education, Race.Ethnicity, Sex, Income,
    Insurance.Coverage, Employment, Marital.Status
  ) %>%
  mutate(
    A = case_when(
      Area_Type == "Urban" ~ 0L,
      Area_Type == "Rural" ~ 1L,
      TRUE ~ NA_integer_
    )
  )

############################################################
# GROUP-SPECIFIC ACS MARGINS
############################################################

acs_prop_urban <- list(
  Age_Group = c(
    "18 to 34 years" = 25.213726,
    "35 to 44 years" = 13.937665,
    "45 to 64 years" = 30.251810,
    "65+ years"      = 30.596798
  ),
  Education = c(
    "Bachelor's Degree or Higher"               = 30.199767,
    "Less than High School to High School Grad" = 36.996966,
    "Post-HS / Some College"                    = 32.803267
  ),
  Employment = c(
    "Employed"   = 51.322268,
    "Others"     = 46.224392,
    "Unemployed" = 2.453340
  ),
  Income = c(
    "$0-$19,999"      = 12.005602,
    "$100,000+"       = 34.026895,
    "$20,000-$49,999" = 23.140250,
    "$50,000-$74,999" = 17.248741,
    "$75,000-$99,999" = 13.578512
  ),
  Insurance.Coverage = c(
    "No"  = 9.773532,
    "Yes" = 90.226468
  ),
  Marital.Status = c(
    "Living alone"                     = 52.401499,
    "Married or living with a partner" = 47.598501
  ),
  Race.Ethnicity = c(
    "All Other Races"    = 14.068061,
    "Hispanic"           = 10.822895,
    "Non-Hispanic Black" = 11.149668,
    "Non-Hispanic White" = 63.959375
  ),
  Sex = c(
    "Female" = 51.777023,
    "Male"   = 48.222977
  )
)

acs_prop_rural <- list(
  Age_Group = c(
    "18 to 34 years" = 24.577741,
    "35 to 44 years" = 16.030568,
    "45 to 64 years" = 32.501265,
    "65+ years"      = 26.890426
  ),
  Education = c(
    "Bachelor's Degree or Higher"               = 13.781333,
    "Less than High School to High School Grad" = 57.194259,
    "Post-HS / Some College"                    = 29.024408
  ),
  Employment = c(
    "Employed"   = 44.611567,
    "Others"     = 52.537612,
    "Unemployed" = 2.850821
  ),
  Income = c(
    "$0-$19,999"      = 17.395816,
    "$100,000+"       = 21.714967,
    "$20,000-$49,999" = 29.586525,
    "$50,000-$74,999" = 19.003919,
    "$75,000-$99,999" = 12.298773
  ),
  Insurance.Coverage = c(
    "No"  = 14.162332,
    "Yes" = 85.837668
  ),
  Marital.Status = c(
    "Living alone"                     = 58.770167,
    "Married or living with a partner" = 41.229833
  ),
  Race.Ethnicity = c(
    "All Other Races"    = 10.221594,
    "Hispanic"           = 7.191337,
    "Non-Hispanic Black" = 15.331070,
    "Non-Hispanic White" = 67.255999
  ),
  Sex = c(
    "Female" = 47.223494,
    "Male"   = 52.776506
  )
)

acs_prop_urban <- lapply(acs_prop_urban, function(x) x / sum(x))
acs_prop_rural <- lapply(acs_prop_rural, function(x) x / sum(x))

margin_vars <- names(acs_prop_urban)

for (v in margin_vars) {
  lev <- union(
    names(acs_prop_urban[[v]]),
    names(acs_prop_rural[[v]])
  )
  cna[[v]] <- factor(cna[[v]], levels = lev)
}

n_before <- nrow(cna)

cna <- cna %>%
  filter(
    !is.na(A),
    if_all(all_of(margin_vars), ~ !is.na(.x))
  )

cat("\nOriginal CNA N = ", n_before, "\n", sep = "")
cat("CNA N used for group-specific ET = ", nrow(cna), "\n", sep = "")
cat("Dropped = ", n_before - nrow(cna), "\n", sep = "")
print(table(cna$Area_Type, cna$A))

############################################################
# BUILD CALIBRATION SYSTEM
############################################################

build_calibration_system <- function(data, margins) {

  H <- matrix(
    1,
    nrow = nrow(data),
    ncol = 1,
    dimnames = list(NULL, "(Intercept)")
  )

  target_mean <- c(`(Intercept)` = 1)

  column_info <- tibble(
    matrix_column = "(Intercept)",
    Variable = "Intercept",
    Category = "Population total",
    ACS_target = 1
  )

  for (v in names(margins)) {

    probs <- margins[[v]]
    categories <- names(probs)

    for (category in categories[-1]) {

      column_name <- paste0(v, "__", category)
      indicator <- as.numeric(data[[v]] == category)

      H <- cbind(H, indicator)
      colnames(H)[ncol(H)] <- column_name

      target_mean[column_name] <- probs[category]

      column_info <- bind_rows(
        column_info,
        tibble(
          matrix_column = column_name,
          Variable = v,
          Category = category,
          ACS_target = as.numeric(probs[category])
        )
      )
    }
  }

  storage.mode(H) <- "double"

  list(
    H = H,
    target_mean = target_mean[colnames(H)],
    column_info = column_info
  )
}

scale_calibration_system <- function(H, target_mean) {

  H <- as.matrix(H)
  target_mean <- as.numeric(target_mean)
  names(target_mean) <- colnames(H)

  scale_vec <- rep(1, ncol(H))
  names(scale_vec) <- colnames(H)

  for (j in seq_len(ncol(H))) {
    if (all(abs(H[, j] - 1) < 1e-12)) {
      scale_vec[j] <- 1
    } else {
      s <- sqrt(mean(H[, j]^2))
      if (!is.finite(s) || s < 1e-8) s <- 1
      scale_vec[j] <- s
    }
  }

  list(
    H = sweep(H, 2, scale_vec, FUN = "/"),
    target_mean = target_mean / scale_vec,
    scale = scale_vec
  )
}

fit_et_weights <- function(
    H,
    target_mean,
    maxit = 5000,
    balance_tol = 1e-6
) {

  H <- as.matrix(H)
  n <- nrow(H)
  p <- ncol(H)

  target_mean <- as.numeric(target_mean)
  target_total <- n * target_mean

  dual_obj <- function(lambda) {
    eta <- as.vector(H %*% lambda)
    if (any(!is.finite(eta)) || max(eta) > 700) return(1e100)
    w <- exp(eta)
    if (any(!is.finite(w))) return(1e100)
    sum(w) - sum(lambda * target_total)
  }

  dual_grad <- function(lambda) {
    w <- exp(as.vector(H %*% lambda))
    as.vector(crossprod(H, w) - target_total)
  }

  fit <- optim(
    par = rep(0, p),
    fn = dual_obj,
    gr = dual_grad,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-12)
  )

  lambda_hat <- fit$par
  w <- exp(as.vector(H %*% lambda_hat))

  weighted_mean <- as.vector(crossprod(H, w) / sum(w))
  balance_error <- max(abs(weighted_mean - target_mean))

  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("Entropy tilting produced invalid weights.")
  }

  if (!is.finite(balance_error) || balance_error > balance_tol) {
    stop(
      paste0(
        "ET balance tolerance not met. Max error = ",
        signif(balance_error, 6)
      )
    )
  }

  list(
    weights = w,
    lambda = lambda_hat,
    convergence = fit$convergence,
    balance_error = balance_error
  )
}

fit_group_et <- function(cna, a, margins) {

  dat <- cna %>% filter(A == a)

  cal <- build_calibration_system(
    data = dat,
    margins = margins
  )

  scaled <- scale_calibration_system(
    H = cal$H,
    target_mean = cal$target_mean
  )

  fit <- fit_et_weights(
    H = scaled$H,
    target_mean = scaled$target_mean
  )

  dat$ET_weight <- fit$weights

  list(
    data = dat,
    H_original = cal$H,
    target_mean_original = cal$target_mean,
    column_info = cal$column_info,
    H = scaled$H,
    target_mean = scaled$target_mean,
    fit = fit
  )
}

############################################################
# FIT URBAN AND RURAL SEPARATELY
############################################################

fit_urban <- fit_group_et(
  cna = cna,
  a = 0L,
  margins = acs_prop_urban
)

fit_rural <- fit_group_et(
  cna = cna,
  a = 1L,
  margins = acs_prop_rural
)

cna_weighted <- bind_rows(
  fit_urban$data,
  fit_rural$data
)

############################################################
# WEIGHT DIAGNOSTICS
############################################################

weight_diagnostics <- cna_weighted %>%
  group_by(Area_Type, A) %>%
  summarise(
    N = n(),
    Min = min(ET_weight),
    Q01 = quantile(ET_weight, 0.01),
    Q05 = quantile(ET_weight, 0.05),
    Median = median(ET_weight),
    Mean = mean(ET_weight),
    Q95 = quantile(ET_weight, 0.95),
    Q99 = quantile(ET_weight, 0.99),
    Max = max(ET_weight),
    ESS = sum(ET_weight)^2 / sum(ET_weight^2),
    .groups = "drop"
  )

weight_diagnostics

############################################################
# BALANCE CHECK
############################################################

make_balance_table <- function(fit_obj, group_label) {

  weighted_target_original <- as.vector(
    crossprod(
      fit_obj$H_original,
      fit_obj$data$ET_weight
    ) /
      sum(fit_obj$data$ET_weight)
  )

  unweighted_original <- colMeans(fit_obj$H_original)

  fit_obj$column_info %>%
    mutate(
      Group = group_label,
      CNA_unweighted = unweighted_original[
        match(matrix_column, colnames(fit_obj$H_original))
      ],
      CNA_ET = weighted_target_original[
        match(matrix_column, colnames(fit_obj$H_original))
      ],
      Difference = CNA_ET - ACS_target
    )
}

balance_table <- bind_rows(
  make_balance_table(fit_urban, "Urban"),
  make_balance_table(fit_rural, "Rural")
)

balance_table_percent <- balance_table %>%
  filter(Variable != "Intercept") %>%
  mutate(
    across(
      c(ACS_target, CNA_unweighted, CNA_ET, Difference),
      ~ 100 * .x
    )
  )

balance_table_percent

############################################################
# GROUP-SPECIFIC EE/SANDWICH
############################################################

et_mean_ee_group <- function(
    y,
    H,
    target_mean,
    w
) {

  H <- as.matrix(H)
  n <- nrow(H)
  p <- ncol(H)

  observed <- !is.na(y)
  y_work <- ifelse(observed, y, 0)
  M <- as.numeric(observed)

  mu_hat <- sum(M * w * y_work) /
    sum(M * w)

  U_lambda <- sweep(
    H * as.numeric(w),
    2,
    target_mean,
    FUN = "-"
  )

  U_mu <- M * w * (y_work - mu_hat)
  U <- cbind(U_lambda, U_mu)

  A_mat <- matrix(
    0,
    nrow = p + 1,
    ncol = p + 1
  )

  A_mat[1:p, 1:p] <-
    crossprod(
      H,
      H * as.numeric(w)
    ) / n

  dUmu_dlambda <-
    H *
    as.numeric(
      M * w * (y_work - mu_hat)
    )

  A_mat[p + 1, 1:p] <-
    colMeans(dUmu_dlambda)

  A_mat[p + 1, p + 1] <-
    -mean(M * w)

  Meat <- crossprod(U) / n
  A_inv <- solve(A_mat)

  V <- A_inv %*%
    Meat %*%
    t(A_inv) /
    n

  var_mu <- V[p + 1, p + 1]
  se_mu <- sqrt(pmax(var_mu, 0))

  tibble(
    N = sum(observed),
    Estimate = mu_hat,
    Variance = var_mu,
    SE = se_mu,
    LCL = mu_hat - 1.96 * se_mu,
    UCL = mu_hat + 1.96 * se_mu
  )
}

############################################################
# URBAN / RURAL RESULTS
############################################################

outcomes <- c(
  "Cancerbeliefs_1",
  "Cancerbeliefs_2",
  "Cancerbeliefs_4"
)

final_results <- map_dfr(
  outcomes,
  function(outcome) {

    urban_res <- et_mean_ee_group(
      y = fit_urban$data[[outcome]],
      H = fit_urban$H,
      target_mean = fit_urban$target_mean,
      w = fit_urban$data$ET_weight
    )

    rural_res <- et_mean_ee_group(
      y = fit_rural$data[[outcome]],
      H = fit_rural$H,
      target_mean = fit_rural$target_mean,
      w = fit_rural$data$ET_weight
    )

    delta <- rural_res$Estimate - urban_res$Estimate

    # Separate group-specific estimating systems use disjoint
    # observations, so the covariance is zero and variances add.
    var_delta <- rural_res$Variance + urban_res$Variance
    se_delta <- sqrt(pmax(var_delta, 0))

    tibble(
      Outcome = outcome,

      Urban_N = urban_res$N,
      Urban_Mean = urban_res$Estimate,
      Urban_SE = urban_res$SE,
      Urban_LCL = urban_res$LCL,
      Urban_UCL = urban_res$UCL,

      Rural_N = rural_res$N,
      Rural_Mean = rural_res$Estimate,
      Rural_SE = rural_res$SE,
      Rural_LCL = rural_res$LCL,
      Rural_UCL = rural_res$UCL,

      Delta_Rural_minus_Urban = delta,
      Delta_SE = se_delta,
      Delta_LCL = delta - 1.96 * se_delta,
      Delta_UCL = delta + 1.96 * se_delta,

      Z = delta / se_delta,
      P_value = 2 * pnorm(-abs(delta / se_delta)),
      Calibration = "Separate Urban/Rural ACS margins"
    )
  }
)

final_results_display <- final_results %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  )

final_results_display

fit_linear_weights <- function(
    H,
    target_mean,
    balance_tol = 1e-8
) {

  H <- as.matrix(H)
  n <- nrow(H)
  p <- ncol(H)

  target_mean <- as.numeric(target_mean)
  target_total <- n * target_mean

  # Solve
  #   H^T(1 + H lambda) = n * target_mean
  # for lambda.
  calibration_matrix <- crossprod(H)
  calibration_rhs <- target_total - colSums(H)

  qr_fit <- qr(calibration_matrix, tol = 1e-10)

  if (qr_fit$rank < p) {
    stop(
      paste0(
        "Linear calibration system is rank deficient: rank = ",
        qr_fit$rank,
        " but p = ",
        p,
        ". Check for empty or redundant calibration categories."
      )
    )
  }

  lambda_hat <- qr.coef(qr_fit, calibration_rhs)
  w <- 1 + as.vector(H %*% lambda_hat)

  if (any(!is.finite(w))) {
    stop("Linear calibration produced non-finite weights.")
  }

  # The intercept in H constrains sum(w) to equal n.
  weighted_mean <- as.vector(crossprod(H, w) / sum(w))
  balance_error <- max(abs(weighted_mean - target_mean))

  if (!is.finite(balance_error) || balance_error > balance_tol) {
    stop(
      paste0(
        "Linear calibration balance tolerance not met. Max error = ",
        signif(balance_error, 6)
      )
    )
  }

  list(
    weights = w,
    lambda = lambda_hat,
    rank = qr_fit$rank,
    balance_error = balance_error,
    n_negative = sum(w < 0),
    n_zero = sum(w == 0)
  )
}

fit_group_linear <- function(cna, a, margins) {

  dat <- cna %>% filter(A == a)

  cal <- build_calibration_system(
    data = dat,
    margins = margins
  )

  # Use the same column scaling as ET for numerical stability.
  scaled <- scale_calibration_system(
    H = cal$H,
    target_mean = cal$target_mean
  )

  fit <- fit_linear_weights(
    H = scaled$H,
    target_mean = scaled$target_mean
  )

  dat$Linear_weight <- fit$weights

  list(
    data = dat,
    H_original = cal$H,
    target_mean_original = cal$target_mean,
    column_info = cal$column_info,
    H = scaled$H,
    target_mean = scaled$target_mean,
    fit = fit
  )
}

############################################################
# FIT LINEAR CALIBRATION WITHIN URBAN AND RURAL
############################################################

fit_urban_linear <- fit_group_linear(
  cna = cna,
  a = 0L,
  margins = acs_prop_urban
)

fit_rural_linear <- fit_group_linear(
  cna = cna,
  a = 1L,
  margins = acs_prop_rural
)

cna_linear_weighted <- bind_rows(
  fit_urban_linear$data,
  fit_rural_linear$data
)

############################################################
# LINEAR-WEIGHT DIAGNOSTICS
############################################################

linear_weight_diagnostics <- cna_linear_weighted %>%
  group_by(Area_Type, A) %>%
  summarise(
    N = n(),
    Min = min(Linear_weight),
    Q01 = quantile(Linear_weight, 0.01),
    Q05 = quantile(Linear_weight, 0.05),
    Median = median(Linear_weight),
    Mean = mean(Linear_weight),
    Q95 = quantile(Linear_weight, 0.95),
    Q99 = quantile(Linear_weight, 0.99),
    Max = max(Linear_weight),
    N_negative = sum(Linear_weight < 0),
    Percent_negative = 100 * mean(Linear_weight < 0),
    N_nonpositive = sum(Linear_weight <= 0),
    Percent_nonpositive = 100 * mean(Linear_weight <= 0),
    ESS = sum(Linear_weight)^2 / sum(Linear_weight^2),
    .groups = "drop"
  )

linear_weight_diagnostics

# Records receiving negative linear-calibration weights, ordered
# from the most negative weight upward. X is retained as the CNA
# record identifier.
negative_linear_weights <- cna_linear_weighted %>%
  filter(Linear_weight < 0) %>%
  arrange(Linear_weight)

negative_linear_weights

if (nrow(negative_linear_weights) > 0) {
  cat(
    "\nLinear calibration produced ",
    nrow(negative_linear_weights),
    " negative weights (",
    round(100 * nrow(negative_linear_weights) / nrow(cna_linear_weighted), 2),
    "% overall).\n",
    sep = ""
  )
} else {
  message(
    "Linear calibration produced no negative weights in these data under the specified margins."
  )
}

# The dashed red line marks zero; bars to its left are negative
# calibration weights.
linear_weight_histogram <- ggplot(
  cna_linear_weighted,
  aes(x = Linear_weight, fill = Area_Type)
) +
  geom_histogram(
    bins = 40,
    color = "white",
    alpha = 0.80
  ) +
  geom_vline(
    xintercept = 0,
    color = "red3",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ Area_Type,
    scales = "free_y"
  ) +
  labs(
    title = "Distribution of group-specific linear calibration weights",
    subtitle = "Weights below zero demonstrate the possibility of negative linear calibration weights",
    x = "Linear calibration weight",
    y = "Number of CNA respondents",
    fill = "Area type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

linear_weight_histogram

############################################################
# LINEAR-CALIBRATION BALANCE CHECK
############################################################

make_linear_balance_table <- function(fit_obj, group_label) {

  weighted_target_original <- as.vector(
    crossprod(
      fit_obj$H_original,
      fit_obj$data$Linear_weight
    ) /
      sum(fit_obj$data$Linear_weight)
  )

  unweighted_original <- colMeans(fit_obj$H_original)

  fit_obj$column_info %>%
    mutate(
      Group = group_label,
      CNA_unweighted = unweighted_original[
        match(matrix_column, colnames(fit_obj$H_original))
      ],
      CNA_linear = weighted_target_original[
        match(matrix_column, colnames(fit_obj$H_original))
      ],
      Difference = CNA_linear - ACS_target
    )
}

linear_balance_table <- bind_rows(
  make_linear_balance_table(fit_urban_linear, "Urban"),
  make_linear_balance_table(fit_rural_linear, "Rural")
)

linear_balance_table_percent <- linear_balance_table %>%
  filter(Variable != "Intercept") %>%
  mutate(
    across(
      c(ACS_target, CNA_unweighted, CNA_linear, Difference),
      ~ 100 * .x
    )
  )

linear_balance_table_percent

############################################################
# LINEAR-CALIBRATION EE/SANDWICH
############################################################

linear_mean_ee_group <- function(
    y,
    H,
    target_mean,
    w
) {

  H <- as.matrix(H)
  n <- nrow(H)
  p <- ncol(H)

  observed <- !is.na(y)
  y_work <- ifelse(observed, y, 0)
  M <- as.numeric(observed)

  mu_hat <- sum(M * w * y_work) /
    sum(M * w)

  # Calibration estimating functions:
  #   U_lambda,i = w_i H_i - target_mean.
  U_lambda <- sweep(
    H * as.numeric(w),
    2,
    target_mean,
    FUN = "-"
  )

  # Weighted-mean estimating function:
  #   U_mu,i = M_i w_i (Y_i - mu).
  U_mu <- M * w * (y_work - mu_hat)
  U <- cbind(U_lambda, U_mu)

  A_mat <- matrix(
    0,
    nrow = p + 1,
    ncol = p + 1
  )

  # Because dw_i/dlambda = H_i for linear calibration,
  # dU_lambda/dlambda = H_i H_i^T.
  A_mat[1:p, 1:p] <- crossprod(H) / n

  # dU_mu/dlambda = M_i (Y_i - mu) H_i.
  dUmu_dlambda <-
    H *
    as.numeric(
      M * (y_work - mu_hat)
    )

  A_mat[p + 1, 1:p] <- colMeans(dUmu_dlambda)
  A_mat[p + 1, p + 1] <- -mean(M * w)

  Meat <- crossprod(U) / n
  A_inv <- solve(A_mat)

  V <- A_inv %*%
    Meat %*%
    t(A_inv) /
    n

  var_mu <- V[p + 1, p + 1]
  se_mu <- sqrt(pmax(var_mu, 0))

  tibble(
    N = sum(observed),
    Estimate = mu_hat,
    Variance = var_mu,
    SE = se_mu,
    LCL = mu_hat - 1.96 * se_mu,
    UCL = mu_hat + 1.96 * se_mu
  )
}

############################################################
# LINEAR-CALIBRATION URBAN / RURAL RESULTS
############################################################

linear_results <- map_dfr(
  outcomes,
  function(outcome) {

    urban_res <- linear_mean_ee_group(
      y = fit_urban_linear$data[[outcome]],
      H = fit_urban_linear$H,
      target_mean = fit_urban_linear$target_mean,
      w = fit_urban_linear$data$Linear_weight
    )

    rural_res <- linear_mean_ee_group(
      y = fit_rural_linear$data[[outcome]],
      H = fit_rural_linear$H,
      target_mean = fit_rural_linear$target_mean,
      w = fit_rural_linear$data$Linear_weight
    )

    delta <- rural_res$Estimate - urban_res$Estimate

    # Urban and Rural observations are disjoint, so variances add.
    var_delta <- rural_res$Variance + urban_res$Variance
    se_delta <- sqrt(pmax(var_delta, 0))

    tibble(
      Outcome = outcome,

      Urban_N = urban_res$N,
      Urban_Mean = urban_res$Estimate,
      Urban_SE = urban_res$SE,
      Urban_LCL = urban_res$LCL,
      Urban_UCL = urban_res$UCL,

      Rural_N = rural_res$N,
      Rural_Mean = rural_res$Estimate,
      Rural_SE = rural_res$SE,
      Rural_LCL = rural_res$LCL,
      Rural_UCL = rural_res$UCL,

      Delta_Rural_minus_Urban = delta,
      Delta_SE = se_delta,
      Delta_LCL = delta - 1.96 * se_delta,
      Delta_UCL = delta + 1.96 * se_delta,

      Z = delta / se_delta,
      P_value = 2 * pnorm(-abs(delta / se_delta)),
      Calibration = "Separate Urban/Rural linear calibration to ACS margins"
    )
  }
)

linear_results_display <- linear_results %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  )

linear_results_display

unweighted_area_one_outcome <- function(
    data,
    outcome
) {

  dat <- data %>%
    filter(!is.na(.data[[outcome]]))

  urban <- dat %>%
    filter(A == 0L)

  rural <- dat %>%
    filter(A == 1L)

  n_u <- nrow(urban)
  n_r <- nrow(rural)

  mu_u <- mean(urban[[outcome]])
  mu_r <- mean(rural[[outcome]])

  var_u <- var(urban[[outcome]]) / n_u
  var_r <- var(rural[[outcome]]) / n_r

  se_u <- sqrt(var_u)
  se_r <- sqrt(var_r)

  delta <- mu_r - mu_u

  # Urban and Rural observations are disjoint, so variances add.
  var_delta <- var_u + var_r
  se_delta <- sqrt(var_delta)

  z <- delta / se_delta
  p_value <- 2 * pnorm(-abs(z))

  tibble(
    Outcome = outcome,

    Urban_N = n_u,
    Urban_Mean = mu_u,
    Urban_SE = se_u,
    Urban_LCL = mu_u - 1.96 * se_u,
    Urban_UCL = mu_u + 1.96 * se_u,

    Rural_N = n_r,
    Rural_Mean = mu_r,
    Rural_SE = se_r,
    Rural_LCL = mu_r - 1.96 * se_r,
    Rural_UCL = mu_r + 1.96 * se_r,

    Delta_Rural_minus_Urban = delta,
    Delta_SE = se_delta,
    Delta_LCL = delta - 1.96 * se_delta,
    Delta_UCL = delta + 1.96 * se_delta,

    Z = z,
    P_value = p_value,
    Method = "Unweighted"
  )
}

unweighted_area_results <- map_dfr(
  outcomes,
  ~ unweighted_area_one_outcome(
    data = cna,
    outcome = .x
  )
)

unweighted_area_results_display <- unweighted_area_results %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )

unweighted_area_results_display

weighted_area_results_for_comparison <- final_results %>%
  mutate(
    Method = "Group-specific ET"
  ) %>%
  select(
    Outcome,
    Method,
    Urban_N,
    Urban_Mean,
    Urban_SE,
    Urban_LCL,
    Urban_UCL,
    Rural_N,
    Rural_Mean,
    Rural_SE,
    Rural_LCL,
    Rural_UCL,
    Delta_Rural_minus_Urban,
    Delta_SE,
    Delta_LCL,
    Delta_UCL,
    Z,
    P_value
  )

linear_area_results_for_comparison <- linear_results %>%
  mutate(
    Method = "Group-specific linear calibration"
  ) %>%
  select(
    Outcome,
    Method,
    Urban_N,
    Urban_Mean,
    Urban_SE,
    Urban_LCL,
    Urban_UCL,
    Rural_N,
    Rural_Mean,
    Rural_SE,
    Rural_LCL,
    Rural_UCL,
    Delta_Rural_minus_Urban,
    Delta_SE,
    Delta_LCL,
    Delta_UCL,
    Z,
    P_value
  )

unweighted_area_results_for_comparison <- unweighted_area_results %>%
  select(
    Outcome,
    Method,
    Urban_N,
    Urban_Mean,
    Urban_SE,
    Urban_LCL,
    Urban_UCL,
    Rural_N,
    Rural_Mean,
    Rural_SE,
    Rural_LCL,
    Rural_UCL,
    Delta_Rural_minus_Urban,
    Delta_SE,
    Delta_LCL,
    Delta_UCL,
    Z,
    P_value
  )

area_results_comparison <- bind_rows(
  unweighted_area_results_for_comparison,
  linear_area_results_for_comparison,
  weighted_area_results_for_comparison
) %>%
  arrange(
    Outcome,
    factor(
      Method,
      levels = c(
        "Unweighted",
        "Group-specific linear calibration",
        "Group-specific ET"
      )
    )
  )

area_results_comparison_display <- area_results_comparison %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )

area_results_comparison_display

area_target_prop <- c(
  Urban = 0.911,
  Rural = 0.089
)

overall_results <- final_results %>%
  transmute(
    Outcome,
    Estimate =
      area_target_prop["Urban"] * Urban_Mean +
      area_target_prop["Rural"] * Rural_Mean,
    Variance =
      area_target_prop["Urban"]^2 * Urban_SE^2 +
      area_target_prop["Rural"]^2 * Rural_SE^2
  ) %>%
  mutate(
    SE = sqrt(Variance),
    LCL = Estimate - 1.96 * SE,
    UCL = Estimate + 1.96 * SE,
    Calibration =
      "Separate Urban/Rural ACS margins combined by target Area_Type proportion"
  )

overall_results

linear_overall_results <- linear_results %>%
  transmute(
    Outcome,
    Estimate =
      area_target_prop["Urban"] * Urban_Mean +
      area_target_prop["Rural"] * Rural_Mean,
    Variance =
      area_target_prop["Urban"]^2 * Urban_SE^2 +
      area_target_prop["Rural"]^2 * Rural_SE^2
  ) %>%
  mutate(
    SE = sqrt(Variance),
    LCL = Estimate - 1.96 * SE,
    UCL = Estimate + 1.96 * SE,
    Calibration =
      "Separate Urban/Rural linear calibration combined by target Area_Type proportion"
  )

linear_overall_results

############################################################
# SAVE OUTPUTS
############################################################

write.csv(
  cna_weighted,
  "CNA_cancerbeliefs_group_specific_ET_weighted.csv",
  row.names = FALSE
)

write.csv(
  balance_table,
  "CNA_group_specific_ET_balance_check.csv",
  row.names = FALSE
)

write.csv(
  overall_results,
  "CNA_group_specific_ET_cancerbeliefs_overall_results.csv",
  row.names = FALSE
)

write.csv(
  final_results,
  "CNA_group_specific_ET_cancerbeliefs_urban_rural_results.csv",
  row.names = FALSE
)


write.csv(
  unweighted_area_results,
  "CNA_unweighted_cancerbeliefs_urban_rural_results.csv",
  row.names = FALSE
)

write.csv(
  area_results_comparison,
  "CNA_unweighted_linear_and_ET_urban_rural_results.csv",
  row.names = FALSE
)

write.csv(
  cna_linear_weighted,
  "CNA_cancerbeliefs_group_specific_linear_weighted.csv",
  row.names = FALSE
)

write.csv(
  linear_weight_diagnostics,
  "CNA_group_specific_linear_weight_diagnostics.csv",
  row.names = FALSE
)

write.csv(
  negative_linear_weights,
  "CNA_group_specific_negative_linear_weights.csv",
  row.names = FALSE
)

ggsave(
  filename = "CNA_group_specific_linear_weight_histogram.png",
  plot = linear_weight_histogram,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)

write.csv(
  linear_balance_table,
  "CNA_group_specific_linear_balance_check.csv",
  row.names = FALSE
)

write.csv(
  linear_overall_results,
  "CNA_group_specific_linear_cancerbeliefs_overall_results.csv",
  row.names = FALSE
)

write.csv(
  linear_results,
  "CNA_group_specific_linear_cancerbeliefs_urban_rural_results.csv",
  row.names = FALSE
)
