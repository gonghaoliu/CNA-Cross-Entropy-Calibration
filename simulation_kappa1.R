############################################################

library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(ggplot2)

############################################################
# 1. Helpers
############################################################

expit <- function(x) {
  1 / (1 + exp(-x))
}


############################################################
# Selection predictors
############################################################

get_selection_predictors <- function(
    data,
    selection_model = c("nonlinear", "linear")
) {

  selection_model <- match.arg(selection_model)

  if (selection_model == "nonlinear") {

    H <- cbind(
      X1 = data$X1,
      X2 = data$X2,
      X4 = data$X4,
      X1_X2_plus_sin_X1 =
        data$X1 * data$X2 +
        sin(data$X1)
    )

  } else {

    H <- cbind(
      X1 = data$X1,
      X2 = data$X2,
      X4 = data$X4,
      X5 = data$X5
    )
  }

  H <- as.matrix(H)
  storage.mode(H) <- "double"

  H
}


############################################################
# Numerical scaling helper
############################################################

scale_calibration_system <- function(H, target) {

  H <- as.matrix(H)
  target <- as.numeric(target)

  p <- ncol(H)

  scale_vec <- rep(1, p)

  for (j in seq_len(p)) {

    # Keep intercept unchanged
    if (
      all(
        abs(
          H[, j] - 1
        ) < 1e-12
      )
    ) {
      scale_vec[j] <- 1
    } else {

      s <- sqrt(
        mean(
          H[, j]^2
        )
      )

      if (
        !is.finite(s) ||
        s < 1e-8
      ) {
        s <- 1
      }

      scale_vec[j] <- s
    }
  }

  H_scaled <-
    sweep(
      H,
      2,
      scale_vec,
      FUN = "/"
    )

  target_scaled <-
    target /
    scale_vec

  colnames(H_scaled) <- colnames(H)

  list(
    H = H_scaled,
    target = target_scaled,
    scale = scale_vec
  )
}


############################################################
# 2. Calibration basis
############################################################

get_basis <- function(
    data,
    basis = c("H0", "HS", "HY", "HF"),
    selection_model = c("nonlinear", "linear"),
    intercept = TRUE
) {

  basis <- match.arg(basis)
  selection_model <- match.arg(selection_model)

  X1 <- data$X1
  X2 <- data$X2
  X3 <- if ("X3" %in% names(data)) data$X3 else NULL
  X4 <- if ("X4" %in% names(data)) data$X4 else NULL

  if (basis == "H0") {

    H <- cbind(
      X1 = X1,
      X2 = X2
    )

  } else if (basis == "HS") {

    H <- get_selection_predictors(
      data = data,
      selection_model = selection_model
    )

  } else if (basis == "HY") {

    H <- cbind(
      X1 = X1,
      X2 = X2,
      X3 = X3,
      X4 = X4,
      X4_sq = X4^2
    )

  } else if (basis == "HF") {

    H <- cbind(
      get_selection_predictors(
        data = data,
        selection_model = selection_model
      ),
      X3 = X3,
      X4_sq = X4^2
    )
  }

  if (intercept) {
    H <- cbind(
      `(Intercept)` = 1,
      H
    )
  }

  H <- as.matrix(H)
  storage.mode(H) <- "double"

  H
}


############################################################
# 3. Generate target population
############################################################

generate_population <- function(
    N = 100000,
    pA = 0.40,
    
    X1_mean = 50,
    X1_sd = 9,
    
    # X2 model
    gamma20 = -1.2,
    gamma21 = 2.4,
    
    # X3 model
    # Rescaled so that the outcome and Delta remain on a
    # moderate numerical scale. X3 | A is still normal and
    # depends only on A.
    X3_mean0 = 0,
    X3_A_effect = 1,
    X3_sd = 1,

    # X4 model: shared selection-outcome variable
    X4_min = 0,
    X4_max = 1,

    # X5 model: selection-only variable for the linear design
    X5_min = 0,
    X5_max = 1,
    
    # Outcome model
    beta0 = 0,
    betaA = 2.0,
    beta1 = 0.05,
    beta2 = 0.50,
    beta3 = 2,

    # X4 has the same linear coefficient in both A groups.
    # The A-by-X4^2 term below creates the group-specific nonlinear
    # component and prevents its bias from canceling in Delta.
    beta4 = 2,

    # Quadratic outcome term absent from HS; set to 0 for no curvature.
    betaA_X4_sq = 6,
    
    sigma_y = 1
) {
  
  A <- rbinom(
    N,
    size = 1,
    prob = pA
  )
  
  
  ##########################################################
  # X1
  ##########################################################
  
  X1 <- rnorm(
    N,
    mean = X1_mean,
    sd = X1_sd
  )
  
  ##########################################################
  # X2
  ##########################################################
  
  p_X2 <- expit(
    gamma20 +
      gamma21 * A
  )
  
  X2 <- rbinom(
    N,
    size = 1,
    prob = p_X2
  )
  
  
  ##########################################################
  # X4
  # Shared variable: it enters both selection designs and the
  # outcome model. Its omission makes H0 genuinely biased.
  ##########################################################

  X4 <- runif(
    N,
    min = X4_min,
    max = X4_max
  )


  ##########################################################
  # X5
  # Selection-only variable used in the linear design. It is
  # excluded from the outcome model so HY remains outcome-
  # correct but selection-incomplete.
  ##########################################################

  X5 <- runif(
    N,
    min = X5_min,
    max = X5_max
  )


  ##########################################################
  # X3
  # Outcome-only variable. It is normal within each A group
  # and depends only on A.
  ##########################################################

  X3 <- rnorm(
    N,
    mean =
      X3_mean0 +
      X3_A_effect * A,
    sd = X3_sd
  )
  
  
  ##########################################################
  # Outcome model
  ##########################################################
  
  m_true <-
    beta0 +
    betaA * A +
    beta1 * X1 +
    beta2 * X2 +
    beta3 * X3 +
    beta4 * X4 +
    betaA_X4_sq * A * X4^2
  
  
  Y <-
    m_true +
    rnorm(
      N,
      mean = 0,
      sd = sigma_y
    )
  
  
  ##########################################################
  # Return full population
  ##########################################################
  
  tibble(
    id = seq_len(N),
    A = A,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    X5 = X5,
    m_true = m_true,
    Y = Y
  )
}


############################################################
# 3B. Analytic superpopulation truth for the default DGP
############################################################

get_superpopulation_truth <- function(
    gamma20 = -1.2,
    gamma21 = 2.4,
    X1_mean = 50,
    X3_mean0 = 0,
    X3_A_effect = 1,
    X4_min = 0,
    X4_max = 1,
    beta0 = 0,
    betaA = 2,
    beta1 = 0.05,
    beta2 = 0.50,
    beta3 = 2,
    beta4 = 2,
    betaA_X4_sq = 6
) {
  mean_X4 <- (X4_min + X4_max) / 2
  mean_X4_sq <- (X4_min^2 + X4_min * X4_max + X4_max^2) / 3

  mu <- vapply(c(0, 1), function(a) {
    mean_X2 <- expit(gamma20 + gamma21 * a)
    mean_X3 <- X3_mean0 + X3_A_effect * a

    beta0 + betaA * a + beta1 * X1_mean + beta2 * mean_X2 +
      beta3 * mean_X3 + beta4 * mean_X4 +
      betaA_X4_sq * a * mean_X4_sq
  }, numeric(1))

  tibble(A = c(0, 1), true_mu = mu)
}


############################################################
# 4. Generate CNA / nonprobability sample
############################################################

generate_sample <- function(
    pop,
    expected_n = 400,
    selection_model = c("nonlinear", "linear"),
    kappa = 1,
    lambda0_base = NULL,
    lambda1_base = NULL
) {

  selection_model <- match.arg(selection_model)

  if (is.null(lambda0_base)) {
    lambda0_base <-
      if (selection_model == "nonlinear") {
        c(
          X1 = 0.001,
          X2 = 0.03,
          X4 = 1.0,
          X1_X2_plus_sin_X1 = 0.02
        )
      } else {
        c(
          X1 = 0.001,
          X2 = 0.03,
          X4 = 1.0,
          X5 = 0.75
        )
      }
  }

  if (is.null(lambda1_base)) {
    lambda1_base <-
      if (selection_model == "nonlinear") {
        c(
          X1 = 0.001,
          X2 = 0.03,
          X4 = 1.5,
          X1_X2_plus_sin_X1 = 0.035
        )
      } else {
        c(
          X1 = 0.001,
          X2 = 0.03,
          X4 = 1.5,
          X5 = 1.0
        )
      }
  }

  N <- nrow(pop)

  if (
    !is.numeric(expected_n) ||
    length(expected_n) != 1 ||
    !is.finite(expected_n) ||
    expected_n <= 0 ||
    expected_n >= N
  ) {
    stop("expected_n must be one number strictly between 0 and N.")
  }
  
  pop$pi_R <- NA_real_
  
  
  ##########################################################
  # Selection generated separately within A = 0 and A = 1
  ##########################################################
  
  for (a in c(0, 1)) {
    
    idx <-
      which(
        pop$A == a
      )
    
    pop_a <-
      pop[
        idx,
        ,
        drop = FALSE
      ]
    
    N_a <-
      nrow(pop_a)
    
    
    ########################################################
    # Expected sample size within A group
    ########################################################
    
    n_expected_a <-
      expected_n *
      N_a / N
    
    target_fraction <-
      n_expected_a /
      N_a
    
    
    ########################################################
    # A-specific selection coefficients
    ########################################################
    
    lambda_base <-
      if (a == 0) {
        
        lambda0_base
        
      } else {
        
        lambda1_base
        
      }
    
    ########################################################
    # True HS selection basis
    ########################################################
    
    HS_no_intercept <- get_selection_predictors(
      data = pop_a,
      selection_model = selection_model
    )

    required_names <- colnames(HS_no_intercept)

    if (
      is.null(names(lambda_base)) ||
      !setequal(names(lambda_base), required_names)
    ) {
      stop(
        "The names of lambda0_base and lambda1_base must match: ",
        paste(required_names, collapse = ", ")
      )
    }

    lambda_star <-
      kappa *
      lambda_base[required_names]
    
    
    ########################################################
    # Selection linear predictor excluding intercept
    ########################################################
    
    selection_part <-
      -as.vector(
        HS_no_intercept %*%
          lambda_star
      )
    
    
    ########################################################
    # Exact exponential selection model
    ########################################################

    selection_max <-
      max(selection_part)

    log_mean_exp_selection <-
      selection_max +
      log(
        mean(
          exp(
            selection_part -
              selection_max
          )
        )
      )

    alpha_a <-
      log(target_fraction) -
      log_mean_exp_selection

    log_pi_a <-
      alpha_a +
      selection_part

    pi_a <-
      exp(log_pi_a)
    
    
    if (
      any(!is.finite(pi_a)) ||
      any(pi_a <= 0) ||
      any(pi_a >= 1)
    ) {
      
      stop(
        paste0(
          "Invalid exponential selection probability; ",
          "reduce kappa or the selection coefficients."
        )
      )
    }
    
    
    pop$pi_R[idx] <-
      pi_a
  }
  
  
  ##########################################################
  # Draw selected sample
  ##########################################################
  
  pop$R <- rbinom(
    nrow(pop),
    size = 1,
    prob = pop$pi_R
  )
  
  
  ##########################################################
  # Observed CNA sample
  ##########################################################
  
  cna <-
    pop %>%
    filter(R == 1) %>%
    select(
      id,
      A,
      X1,
      X2,
      X3,
      X4,
      X5,
      Y
    )
  
  
  ##########################################################
  # Return
  ##########################################################
  
  list(
    pop = pop,
    cna = cna,
    selection_model = selection_model,
    expected_n = expected_n
  )
}


############################################################
# 5. Population marginal summaries ONLY
############################################################

make_marginal_targets <- function(pop) {

  overall <- pop %>%
    group_by(A) %>%
    summarise(
      N_a = n(),
      prop_X2_1 = mean(X2 == 1),
      mean_sin_X1 = mean(sin(X1)),
      mean_X3 = mean(X3),
      mean_X4 = mean(X4),
      mean_X4_sq = mean(X4^2),
      mean_X5 = mean(X5),
      .groups = "drop"
    )

  x1_by_x2 <- pop %>%
    group_by(A, X2) %>%
    summarise(
      N_cell = n(),
      prop_cell =
        n() /
        sum(pop$A == first(A)),
      mean_X1 = mean(X1),
      .groups = "drop"
    )

  list(
    overall = overall,
    X1_by_X2 = x1_by_x2
  )
}


############################################################
# 6. Recover population totals from marginal summaries
############################################################

get_marginal_components <- function(
    margins,
    a
) {

  m <- margins$overall %>%
    filter(A == a)

  m0 <- margins$X1_by_X2 %>%
    filter(A == a, X2 == 0)

  m1 <- margins$X1_by_X2 %>%
    filter(A == a, X2 == 1)

  if (
    nrow(m) != 1 ||
    nrow(m0) != 1 ||
    nrow(m1) != 1
  ) {
    stop(
      "Marginal summaries incomplete for A = ",
      a
    )
  }

  N_a <- m$N_a

  p1 <- m$prop_X2_1
  p0 <- 1 - p1

  mean_X1 <-
    p0 * m0$mean_X1 +
    p1 * m1$mean_X1

  mean_X1_X2 <-
    p1 *
    m1$mean_X1

  mean_X1_X2_plus_sin_X1 <-
    mean_X1_X2 +
    m$mean_sin_X1

  list(
    N_a = N_a,
    mean_X1 = mean_X1,
    mean_X2 = p1,
    mean_X1_X2_plus_sin_X1 =
      mean_X1_X2_plus_sin_X1,
    mean_X3 = m$mean_X3,
    mean_X4 = m$mean_X4,
    mean_X4_sq = m$mean_X4_sq,
    mean_X5 = m$mean_X5
  )
}


get_target_totals <- function(
    margins,
    a,
    basis = c("H0", "HS", "HY", "HF"),
    selection_model = c("nonlinear", "linear")
) {

  basis <- match.arg(basis)
  selection_model <- match.arg(selection_model)

  z <- get_marginal_components(
    margins = margins,
    a = a
  )

  if (basis == "H0") {

    tq <- c(
      `(Intercept)` = z$N_a,
      X1 = z$N_a * z$mean_X1,
      X2 = z$N_a * z$mean_X2
    )

  } else if (basis == "HS") {

    if (selection_model == "nonlinear") {

      tq <- c(
        `(Intercept)` = z$N_a,
        X1 = z$N_a * z$mean_X1,
        X2 = z$N_a * z$mean_X2,
        X4 = z$N_a * z$mean_X4,
        X1_X2_plus_sin_X1 =
          z$N_a * z$mean_X1_X2_plus_sin_X1
      )

    } else {

      tq <- c(
        `(Intercept)` = z$N_a,
        X1 = z$N_a * z$mean_X1,
        X2 = z$N_a * z$mean_X2,
        X4 = z$N_a * z$mean_X4,
        X5 = z$N_a * z$mean_X5
      )
    }

  } else if (basis == "HY") {

    tq <- c(
      `(Intercept)` = z$N_a,
      X1 = z$N_a * z$mean_X1,
      X2 = z$N_a * z$mean_X2,
      X3 = z$N_a * z$mean_X3,
      X4 = z$N_a * z$mean_X4,
      X4_sq = z$N_a * z$mean_X4_sq
    )

  } else {

    if (selection_model == "nonlinear") {

      tq <- c(
        `(Intercept)` = z$N_a,
        X1 = z$N_a * z$mean_X1,
        X2 = z$N_a * z$mean_X2,
        X4 = z$N_a * z$mean_X4,
        X1_X2_plus_sin_X1 =
          z$N_a * z$mean_X1_X2_plus_sin_X1,
        X3 = z$N_a * z$mean_X3,
        X4_sq = z$N_a * z$mean_X4_sq
      )

    } else {

      tq <- c(
        `(Intercept)` = z$N_a,
        X1 = z$N_a * z$mean_X1,
        X2 = z$N_a * z$mean_X2,
        X4 = z$N_a * z$mean_X4,
        X5 = z$N_a * z$mean_X5,
        X3 = z$N_a * z$mean_X3,
        X4_sq = z$N_a * z$mean_X4_sq
      )
    }
  }

  as.numeric(tq) |>
    setNames(names(tq))
}


############################################################
# 7. Entropy calibration + EE variance
############################################################

fit_entropy_ee_group <- function(
    cna,
    margins,
    a,
    basis = c("H0", "HS", "HY", "HF"),
    selection_model = c("nonlinear", "linear"),
    maxit = 2000,
    balance_tol = 1e-6
) {

  basis <- match.arg(basis)
  selection_model <- match.arg(selection_model)

  cna_a <- cna %>% filter(A == a)
  n_a <- nrow(cna_a)

  fail <- function(reason) {
    list(success = FALSE, reason = reason)
  }

  if (n_a < 5) {
    return(fail("n_a < 5"))
  }

  ##########################################################
  # Manuscript notation:
  # H(X) = {1, H_tilde(X)^T}^T.
  # The ET exponent contains H_tilde only. Normalization
  # handles the intercept, so no intercept parameter is fitted.
  ##########################################################

  H_original <- get_basis(
    cna_a,
    basis = basis,
    selection_model = selection_model,
    intercept = FALSE
  )

  N_a <- get_marginal_components(margins, a)$N_a

  target_total_H <- get_target_totals(
    margins = margins,
    a = a,
    basis = basis,
    selection_model = selection_model
  )

  # Drop the intercept total and convert totals to benchmark means.
  target_mean_original <- as.numeric(target_total_H[-1]) / N_a
  names(target_mean_original) <- names(target_total_H)[-1]

  scale_vec <- apply(H_original, 2, sd)
  scale_vec[!is.finite(scale_vec) | scale_vec < 1e-8] <- 1

  H <- sweep(H_original, 2, target_mean_original, FUN = "-")
  H <- sweep(H, 2, scale_vec, FUN = "/")
  target_mean <- rep(0, ncol(H))
  names(target_mean) <- colnames(H)
  k <- ncol(H)

  # Manuscript: d_i=1 and q_ai=d_i/sum(d_i)=1/n_a.
  d_i <- rep(1, n_a)
  q_ai <- d_i / sum(d_i)

  dual_obj <- function(lambda) {
    eta <- as.vector(H %*% lambda)
    if (any(!is.finite(eta))) return(1e100)
    eta_max <- max(eta)
    log_denom <- eta_max + log(sum(d_i * exp(eta - eta_max)))
    log_denom - sum(lambda * target_mean)
  }

  dual_grad <- function(lambda) {
    eta <- as.vector(H %*% lambda)
    if (any(!is.finite(eta))) return(rep(1e100, k))
    eta_max <- max(eta)
    u_stable <- d_i * exp(eta - eta_max)
    p_norm <- u_stable / sum(u_stable)
    as.vector(crossprod(H, p_norm) - target_mean)
  }

  fit_lambda <- tryCatch(
    optim(
      par = rep(0, k),
      fn = dual_obj,
      gr = dual_grad,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-10)
    ),
    error = function(e) NULL
  )

  if (is.null(fit_lambda)) {
    return(fail("ET optimization error"))
  }

  lambda_hat <- fit_lambda$par

  newton_converged <- FALSE

  for (newton_iter in seq_len(50)) {
    eta_newton <- as.vector(H %*% lambda_hat)

    if (any(!is.finite(eta_newton))) {
      break
    }

    eta_newton_max <- max(eta_newton)
    u_newton <- d_i * exp(eta_newton - eta_newton_max)
    p_newton <- u_newton / sum(u_newton)

    gradient_newton <-
      as.vector(crossprod(H, p_newton) - target_mean)

    if (max(abs(gradient_newton)) < 1e-12) {
      newton_converged <- TRUE
      break
    }

    H_bar_newton <- as.vector(crossprod(H, p_newton))
    H_deviation <- sweep(H, 2, H_bar_newton, FUN = "-")
    hessian_newton <-
      crossprod(H_deviation, H_deviation * as.numeric(p_newton))

    newton_direction <- tryCatch(
      solve(
        hessian_newton + diag(1e-12, k),
        gradient_newton
      ),
      error = function(e) NULL
    )

    if (is.null(newton_direction) ||
        any(!is.finite(newton_direction))) {
      break
    }

    current_moment_error <- max(abs(gradient_newton))
    step_size <- 1
    accepted_step <- FALSE

    for (line_search_iter in seq_len(30)) {
      candidate_lambda <-
        lambda_hat - step_size * newton_direction
      candidate_eta <- as.vector(H %*% candidate_lambda)
      candidate_eta_max <- max(candidate_eta)
      candidate_u <-
        d_i * exp(candidate_eta - candidate_eta_max)
      candidate_p <- candidate_u / sum(candidate_u)
      candidate_gradient <-
        as.vector(crossprod(H, candidate_p) - target_mean)
      candidate_moment_error <-
        max(abs(candidate_gradient))

      # Use the calibration residual itself for line search.
      # This avoids another objective-change stopping rule when
      # the dual objective is already flat to machine precision.
      if (is.finite(candidate_moment_error) &&
          candidate_moment_error < current_moment_error) {
        lambda_hat <- candidate_lambda
        accepted_step <- TRUE
        break
      }

      step_size <- step_size / 2
    }

    if (!accepted_step) {
      break
    }
  }

  newton_converged <-
    max(abs(dual_grad(lambda_hat))) < 1e-12

  eta_hat <- as.vector(H %*% lambda_hat)

  if (any(!is.finite(eta_hat))) {
    return(fail("nonfinite ET linear predictor"))
  }

  eta_max <- max(eta_hat)
  u_stable <- d_i * exp(eta_hat - eta_max)
  p_norm <- as.vector(u_stable / sum(u_stable))

  if (any(!is.finite(p_norm)) || any(p_norm <= 0)) {
    return(fail("invalid normalized ET weights"))
  }

  weighted_mean_original <- colSums(H_original * p_norm)
  balance_error <- max(abs(weighted_mean_original - target_mean_original))
  normalization_error <- abs(sum(p_norm) - 1)

  if (!is.finite(balance_error) || balance_error > balance_tol ||
      normalization_error > balance_tol) {
    return(fail(paste0(
      "ET balance error = ", signif(balance_error, 4),
      "; normalization error = ", signif(normalization_error, 4),
      "; optim convergence = ", fit_lambda$convergence,
      "; Newton convergence = ", newton_converged,
      "; scaled moment error = ",
      signif(max(abs(dual_grad(lambda_hat))), 4)
    )))
  }

  Y <- cna_a$Y
  mu_hat <- sum(p_norm * Y)

  u <- u_stable
  H_centered <- sweep(H, 2, target_mean, FUN = "-")

  U_lambda <- H_centered * as.numeric(u)
  U_mu <- u * (Y - mu_hat)
  U <- cbind(U_lambda, U_mu)

  A_mat <- matrix(0, nrow = k + 1, ncol = k + 1)

  # Derivative of u_i(H_i-Hbar) with respect to lambda.
  A_mat[1:k, 1:k] <-
    crossprod(H_centered, H * as.numeric(u)) / n_a

  # Derivative of u_i(Y_i-mu) with respect to lambda and mu.
  A_mat[k + 1, 1:k] <-
    colMeans(H * as.numeric(u * (Y - mu_hat)))
  A_mat[k + 1, k + 1] <- -mean(u)

  Meat <- crossprod(U) / n_a
  A_inv <- tryCatch(solve(A_mat), error = function(e) NULL)

  if (is.null(A_inv)) {
    return(fail("singular ET joint-EE bread"))
  }

  sandwich <- A_inv %*% Meat %*% t(A_inv) / n_a
  var_mu <- sandwich[k + 1, k + 1]

  if (!is.finite(var_mu) || var_mu < 0) {
    return(fail("invalid ET variance"))
  }

  list(
    success = TRUE,
    reason = "ok",
    mu = mu_hat,
    se = sqrt(var_mu),
    ess = 1 / sum(p_norm^2),
    weights = p_norm,
    q_ai = q_ai,
    balance_error = balance_error,
    normalization_error = normalization_error,
    optim_convergence = fit_lambda$convergence,
    newton_convergence = newton_converged,
    scaled_moment_error = max(abs(dual_grad(lambda_hat)))
  )
}

fit_entropy_ee <- function(
    cna,
    margins,
    basis,
    selection_model = c("nonlinear", "linear")
) {

  selection_model <- match.arg(selection_model)

  fit0 <-
    fit_entropy_ee_group(
      cna = cna,
      margins = margins,
      a = 0,
      basis = basis,
      selection_model = selection_model
    )

  fit1 <-
    fit_entropy_ee_group(
      cna = cna,
      margins = margins,
      a = 1,
      basis = basis,
      selection_model = selection_model
    )

  if (
    !fit0$success ||
    !fit1$success
  ) {
    return(
      tibble(success = FALSE)
    )
  }

  tibble(
    success = TRUE,
    mu0 = fit0$mu,
    mu1 = fit1$mu,
    se_mu0 = fit0$se,
    se_mu1 = fit1$se,
    ess0 = fit0$ess,
    ess1 = fit1$ess
  )
}


############################################################
# 8. Linear calibration + EE variance
############################################################

fit_linear_group <- function(
    cna,
    margins,
    a,
    basis = "HS",
    selection_model = c("nonlinear", "linear"),
    balance_tol = 1e-6
) {

  selection_model <- match.arg(selection_model)

  cna_a <-
    cna %>%
    filter(A == a)

  n_a <- nrow(cna_a)

  N_a <-
    get_marginal_components(
      margins,
      a
    )$N_a

  if (n_a < 5) {
    return(
      list(success = FALSE)
    )
  }

  Qs <- get_basis(
    cna_a,
    basis = basis,
    selection_model = selection_model,
    intercept = TRUE
  )

  tq <- get_target_totals(
    margins = margins,
    a = a,
    basis = basis,
    selection_model = selection_model
  )

  Q_original <- Qs
  target_original <- as.numeric(tq) / N_a
  scaled <- scale_calibration_system(Q_original, target_original)
  Qs <- scaled$H
  tq <- scaled$target
  d <- rep(1 / n_a, n_a)

  p <- ncol(Qs)

  M <-
    crossprod(
      Qs,
      Qs *
        as.numeric(d)
    )

  rhs <-
    tq -
    as.vector(
      crossprod(
        Qs,
        d
      )
    )

  eta_hat <- tryCatch(
    solve(
      M,
      rhs
    ),
    error = function(e) NULL
  )

  if (is.null(eta_hat)) {
    return(
      list(success = FALSE)
    )
  }

  w <-
    as.vector(
      d *
        (
          1 +
            Qs %*%
            eta_hat
        )
    )

  if (any(!is.finite(w))) {
    return(
      list(success = FALSE)
    )
  }

  negative_weight <-
    any(w < 0)

  weighted_total <-
    colSums(
      Qs *
        w
    )

  balance_error <-
    max(
      abs(
        weighted_total -
          tq
      ) /
        pmax(
          1,
          abs(tq)
        )
    )

  if (
    !is.finite(balance_error) ||
    balance_error >
    balance_tol
  ) {
    return(
      list(success = FALSE)
    )
  }

  Y <- cna_a$Y

  if (
    abs(sum(w)) <
    1e-10
  ) {
    return(
      list(success = FALSE)
    )
  }

  # Manuscript: mu_hat = sum_i p_i Y_i; intercept balance gives sum(p)=1.
  mu_hat <- sum(w * Y)

  # EE variance
  U_eta <-
    sweep(
      Qs *
        as.numeric(w),
      2,
      tq / n_a,
      FUN = "-"
    )

  U_mu <-
    w *
    (
      Y -
        mu_hat
    )

  U <- cbind(
    U_eta,
    U_mu
  )

  A_mat <- matrix(
    0,
    nrow = p + 1,
    ncol = p + 1
  )

  A_mat[
    1:p,
    1:p
  ] <-
    crossprod(
      Qs,
      Qs *
        as.numeric(d)
    ) /
    n_a

  dUmu_deta <-
    Qs *
    as.numeric(
      d *
        (
          Y -
            mu_hat
        )
    )

  A_mat[
    p + 1,
    1:p
  ] <-
    colMeans(
      dUmu_deta
    )

  A_mat[
    p + 1,
    p + 1
  ] <-
    -mean(w)

  Meat <-
    crossprod(U) /
    n_a

  A_inv <- tryCatch(
    solve(A_mat),
    error = function(e) NULL
  )

  if (is.null(A_inv)) {
    return(
      list(success = FALSE)
    )
  }

  sandwich <-
    A_inv %*%
    Meat %*%
    t(A_inv) /
    n_a

  var_mu <-
    sandwich[
      p + 1,
      p + 1
    ]

  if (
    !is.finite(var_mu) ||
    var_mu < 0
  ) {
    return(
      list(success = FALSE)
    )
  }

  ess <-
    sum(w)^2 /
    sum(w^2)

  list(
    success = TRUE,
    mu = mu_hat,
    se = sqrt(var_mu),
    ess = ess,
    negative_weight = negative_weight,
    weights = w,
    normalized_weights = w,
    expansion_weights = N_a * w,
    balance_error = balance_error,
    normalization_error = abs(sum(w) - 1)
  )
}


fit_linear <- function(
    cna,
    margins,
    basis = "HS",
    selection_model = c("nonlinear", "linear")
) {

  selection_model <- match.arg(selection_model)

  fit0 <-
    fit_linear_group(
      cna = cna,
      margins = margins,
      a = 0,
      basis = basis,
      selection_model = selection_model
    )

  fit1 <-
    fit_linear_group(
      cna = cna,
      margins = margins,
      a = 1,
      basis = basis,
      selection_model = selection_model
    )

  if (
    !fit0$success ||
    !fit1$success
  ) {
    return(
      tibble(success = FALSE)
    )
  }

  tibble(
    success = TRUE,
    mu0 = fit0$mu,
    mu1 = fit1$mu,
    se_mu0 = fit0$se,
    se_mu1 = fit1$se,
    ess0 = fit0$ess,
    ess1 = fit1$ess,
    negative_weight =
      fit0$negative_weight |
      fit1$negative_weight
  )
}


############################################################
# 9. Conventional joint-cell logistic IPW (NOT calibration)
############################################################

make_joint_cells <- function(pop, selection_model = c("nonlinear", "linear")) {
  selection_model <- match.arg(selection_model)
  # X3 is not used by either IPW specification.
  vars <- if (selection_model == "nonlinear") c("X1", "X2", "X4") else
    c("X1", "X2", "X4", "X5")
  pop %>%
    group_by(across(all_of(c("A", vars)))) %>%
    summarise(N_cell = n(), n_cell = sum(R), .groups = "drop")
}

fit_logit_ipw_group <- function(
    cna, joint_cells, a, basis = c("HS", "H0"),
    selection_model = c("nonlinear", "linear")) {
  basis <- match.arg(basis)
  selection_model <- match.arg(selection_model)
  fail <- function(reason) list(success = FALSE, reason = reason)
  ca <- cna[cna$A == a, , drop = FALSE]
  cells <- joint_cells[joint_cells$A == a, , drop = FALSE]
  if (nrow(ca) < 5 || nrow(cells) == 0) return(fail("insufficient observations"))
  if (any(cells$n_cell < 0 | cells$n_cell > cells$N_cell) ||
      sum(cells$n_cell) != nrow(ca)) return(fail("inconsistent cell counts"))
  # Within H0, collapse the richer joint table to the fitted profile.
  if (basis == "H0") {
    cells <- cells %>% group_by(A, X1, X2) %>%
      summarise(N_cell = sum(N_cell), n_cell = sum(n_cell), .groups = "drop")
  }
  H <- get_basis(cells, basis, selection_model)
  Hs <- get_basis(ca, basis, selection_model)
  # Center/scale predictors for stable logistic likelihood fitting.
  center <- colSums(H * cells$N_cell) / sum(cells$N_cell)
  center[1] <- 0
  Z <- sweep(H, 2, center, "-")
  scale <- sqrt(colSums(Z^2 * cells$N_cell) / sum(cells$N_cell))
  scale[1] <- 1
  if (any(!is.finite(scale) | scale < 1e-12)) return(fail("constant predictor"))
  Z <- sweep(Z, 2, scale, "/")
  Zs <- sweep(sweep(Hs, 2, center, "-"), 2, scale, "/")
  fit <- tryCatch(glm.fit(
    x = Z, y = cbind(cells$n_cell, cells$N_cell - cells$n_cell),
    family = binomial(link = "logit"),
    control = glm.control(epsilon = 1e-10, maxit = 100)
  ), error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged) ||
      any(!is.finite(fit$coefficients))) return(fail("logistic fit failed"))
  pc <- plogis(drop(Z %*% fit$coefficients))
  ps <- plogis(drop(Zs %*% fit$coefficients))
  if (any(ps <= 0 | ps >= 1) || any(pc <= 0 | pc >= 1))
    return(fail("boundary propensity"))
  w <- 1 / ps
  mu <- sum(w * ca$Y) / sum(w)
  u <- w * (ca$Y - mu)
  k <- ncol(Z)
  # Joint empirical sandwich: individual Bernoulli score contributions
  # reconstructed from cell counts, NOT treating each cell as one person.
  # Selected: [(1-pi)H, w(Y-mu)]; nonselected: [-pi H, 0].
  bread <- matrix(0, k + 1, k + 1)
  bread[1:k, 1:k] <- -crossprod(Z, Z * (cells$N_cell * pc * (1-pc)))
  bread[k+1, 1:k] <- -colSums(Zs * (u * (1-ps)))
  bread[k+1, k+1] <- -sum(w)
  meat <- matrix(0, k + 1, k + 1)
  meat[1:k, 1:k] <- crossprod(Z, Z *
    (cells$n_cell * (1-pc)^2 + (cells$N_cell-cells$n_cell) * pc^2))
  meat[1:k, k+1] <- colSums(Zs * (u * (1-ps)))
  meat[k+1, 1:k] <- meat[1:k, k+1]
  meat[k+1, k+1] <- sum(u^2)
  inv <- tryCatch(solve(bread), error = function(e) NULL)
  if (is.null(inv)) return(fail("singular joint IPW bread"))
  V <- inv %*% meat %*% t(inv)
  var_mu <- V[k+1, k+1]
  if (!is.finite(var_mu) || var_mu < 0) return(fail("invalid IPW variance"))
  # This diagnostic need NOT be zero: IPW is not calibrated.
  target_mean <- colSums(H * cells$N_cell) / sum(cells$N_cell)
  imbalance <- colSums(Hs * w) / sum(w) - target_mean
  list(success = TRUE, reason = "ok", mu = mu, se = sqrt(var_mu),
       ess = sum(w)^2 / sum(w^2), weights = w,
       max_abs_mean_imbalance = max(abs(imbalance)),
       score_error = max(abs(crossprod(Z, cells$n_cell-cells$N_cell*pc))) /
         sum(cells$N_cell), n_cells = nrow(cells))
}

fit_logit_ipw <- function(cna, joint_cells, basis = c("HS", "H0"),
                          selection_model = c("nonlinear", "linear")) {
  basis <- match.arg(basis)
  selection_model <- match.arg(selection_model)
  f0 <- fit_logit_ipw_group(cna, joint_cells, 0, basis, selection_model)
  f1 <- fit_logit_ipw_group(cna, joint_cells, 1, basis, selection_model)
  if (!isTRUE(f0$success) || !isTRUE(f1$success)) return(tibble(success = FALSE))
  tibble(success = TRUE, mu0 = f0$mu, mu1 = f1$mu,
         se_mu0 = f0$se, se_mu1 = f1$se, ess0 = f0$ess, ess1 = f1$ess,
         ipw_imbalance0 = f0$max_abs_mean_imbalance,
         ipw_imbalance1 = f1$max_abs_mean_imbalance)
}

############################################################
# 9B. Optional diagnostic helper
############################################################

diagnose_one_dataset <- function(cna, margins, joint_cells, selection_model) {
  list(
    ET_HS_A0 = fit_entropy_ee_group(cna, margins, 0, "HS", selection_model),
    ET_HS_A1 = fit_entropy_ee_group(cna, margins, 1, "HS", selection_model),
    ET_HY_A0 = fit_entropy_ee_group(cna, margins, 0, "HY", selection_model),
    ET_HY_A1 = fit_entropy_ee_group(cna, margins, 1, "HY", selection_model),
    ET_HF_A0 = fit_entropy_ee_group(cna, margins, 0, "HF", selection_model),
    ET_HF_A1 = fit_entropy_ee_group(cna, margins, 1, "HF", selection_model),
    Linear_HS_A0 = fit_linear_group(cna, margins, 0, "HS", selection_model),
    Linear_HS_A1 = fit_linear_group(cna, margins, 1, "HS", selection_model),
    IPW_HS_A0 = fit_logit_ipw_group(cna, joint_cells, 0, "HS", selection_model),
    IPW_HS_A1 = fit_logit_ipw_group(cna, joint_cells, 1, "HS", selection_model),
    IPW_H0_A0 = fit_logit_ipw_group(cna, joint_cells, 0, "H0", selection_model),
    IPW_H0_A1 = fit_logit_ipw_group(cna, joint_cells, 1, "H0", selection_model)
  )
}

run_all_scenarios <- function(
    cna,
    margins,
    joint_cells,
    selection_model = c("nonlinear", "linear")
) {

  selection_model <- match.arg(selection_model)

  # S1: ET HS
  S1 <- tryCatch(
    fit_entropy_ee(
      cna = cna,
      margins = margins,
      basis = "HS",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S1) ||
    !isTRUE(S1$success[1])
  ) {
    return(NULL)
  }

  S1 <- S1 %>%
    mutate(
      scenario = "S1",
      method = "ET_HS"
    )

  # S2: ET HY
  S2 <- tryCatch(
    fit_entropy_ee(
      cna = cna,
      margins = margins,
      basis = "HY",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S2) ||
    !isTRUE(S2$success[1])
  ) {
    return(NULL)
  }

  S2 <- S2 %>%
    mutate(
      scenario = "S2",
      method = "ET_HY"
    )

  # S3: ET HF
  S3 <- tryCatch(
    fit_entropy_ee(
      cna = cna,
      margins = margins,
      basis = "HF",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S3) ||
    !isTRUE(S3$success[1])
  ) {
    return(NULL)
  }

  S3 <- S3 %>%
    mutate(
      scenario = "S3",
      method = "ET_HF"
    )

  # S4: Linear HS
  S4 <- tryCatch(
    fit_linear(
      cna = cna,
      margins = margins,
      basis = "HS",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S4) ||
    !isTRUE(S4$success[1])
  ) {
    return(NULL)
  }

  S4 <- S4 %>%
    mutate(
      scenario = "S4",
      method = "Linear_HS"
    )

  # S5: joint-cell binomial likelihood IPW using HS.
  S5 <- tryCatch(
    fit_logit_ipw(
      cna = cna,
      joint_cells = joint_cells,
      basis = "HS",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S5) ||
    !isTRUE(S5$success[1])
  ) {
    return(NULL)
  }

  S5 <- S5 %>%
    mutate(
      scenario = "S5",
      method = "IPW_logit_HS"
    )

  # S6a: ET H0
  S6_ET <- tryCatch(
    fit_entropy_ee(
      cna = cna,
      margins = margins,
      basis = "H0",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S6_ET) ||
    !isTRUE(S6_ET$success[1])
  ) {
    return(NULL)
  }

  S6_ET <- S6_ET %>%
    mutate(
      scenario = "S6a",
      method = "ET_H0"
    )

  # S6b: joint-cell binomial likelihood IPW using H0.
  S6_IPW <- tryCatch(
    fit_logit_ipw(
      cna = cna,
      joint_cells = joint_cells,
      basis = "H0",
      selection_model = selection_model
    ),
    error = function(e) NULL
  )

  if (
    is.null(S6_IPW) ||
    !isTRUE(S6_IPW$success[1])
  ) {
    return(NULL)
  }

  S6_IPW <- S6_IPW %>%
    mutate(
      scenario = "S6b",
      method = "IPW_logit_H0"
    )

  bind_rows(
    S1,
    S2,
    S3,
    S4,
    S5,
    S6_ET,
    S6_IPW
  )
}


############################################################
# 11. Main simulation
############################################################

run_simulation <- function(
    n_success = 1000,
    N = 100000,
    expected_n = 400,
    selection_model = c("nonlinear", "linear"),
    kappa = 1,
    lambda0_base = NULL,
    lambda1_base = NULL,
    seed = 20260814
) {

  selection_model <- match.arg(selection_model)

  set.seed(seed)

  result_storage <-
    vector(
      "list",
      n_success
    )

  linear_weight_storage <-
    vector(
      "list",
      n_success
    )

  success_count <- 0
  attempt_count <- 0

  while (
    success_count <
    n_success
  ) {

    attempt_count <- attempt_count + 1
    if (attempt_count > 5 * n_success) {
      stop("Too many failed attempts: ", success_count, "/", attempt_count - 1,
           ". Inspect diagnose_one_dataset(); do not hide failure rates.")
    }

    # Progress diagnostics: shows whether R is still moving
    # even when a replication is rejected before success.
    if (attempt_count %% 10 == 0) {
      message(
        paste0(
          "Attempt: ",
          attempt_count,
          " | Successful: ",
          success_count,
          "/",
          n_success,
          " | Selection: ",
          selection_model,
          " | Expected n: ",
          expected_n
        )
      )
    }

    ########################################################
    # DGP population
    ########################################################

    pop0 <-
      generate_population(
        N = N
      )

    ########################################################
    # Draw selected sample
    ########################################################

    sample_object <- tryCatch(
      generate_sample(
        pop = pop0,
        expected_n = expected_n,
        selection_model = selection_model,
        kappa = kappa,
        lambda0_base = lambda0_base,
        lambda1_base = lambda1_base
      ),
      error = function(e) NULL
    )

    if (is.null(sample_object)) {
      if (attempt_count %% 10 == 0) {
        message("  Last attempt failed during sample generation.")
      }
      next
    }

    # Full pop stays inside the simulation loop only
    pop <-
      sample_object$pop

    # Estimators receive this selected sample
    cna <-
      sample_object$cna

    if (
      length(
        unique(cna$A)
      ) < 2
    ) {
      if (attempt_count %% 10 == 0) {
        message("  Last attempt failed because one A group was absent.")
      }
      next
    }

    truth <- get_superpopulation_truth()
    true_mu0 <- truth$true_mu[truth$A == 0]
    true_mu1 <- truth$true_mu[truth$A == 1]

    finite_truth <- pop %>%
      group_by(A) %>%
      summarise(finite_mu = mean(Y), .groups = "drop")
    finite_mu0 <- finite_truth$finite_mu[finite_truth$A == 0]
    finite_mu1 <- finite_truth$finite_mu[finite_truth$A == 1]

    margins <-
      make_marginal_targets(
        pop
      )

    joint_cells <- make_joint_cells(pop, selection_model)

    fits <- tryCatch(
      run_all_scenarios(
        cna = cna,
        margins = margins,
        joint_cells = joint_cells,
        selection_model = selection_model
      ),
      error = function(e) NULL
    )

    if (is.null(fits)) {
      if (attempt_count %% 10 == 0) {
        diagnostic_fits <- tryCatch(
          diagnose_one_dataset(cna, margins, joint_cells, selection_model),
          error = function(e) NULL
        )

        if (is.null(diagnostic_fits)) {
          message("  Last attempt failed; diagnostic refit also failed.")
        } else {
          failed_names <- names(diagnostic_fits)[
            !vapply(diagnostic_fits, function(x) isTRUE(x$success), logical(1))
          ]
          failed_reasons <- vapply(
            diagnostic_fits[failed_names],
            function(x) if (is.null(x$reason)) "unspecified" else x$reason,
            character(1)
          )
          message(
            "  Failed methods: ",
            paste0(failed_names, " [", failed_reasons, "]", collapse = "; ")
          )
        }
      }
      next
    }

    linear0 <- tryCatch(
      fit_linear_group(
        cna = cna,
        margins = margins,
        a = 0,
        basis = "HS",
        selection_model = selection_model
      ),
      error = function(e) NULL
    )

    linear1 <- tryCatch(
      fit_linear_group(
        cna = cna,
        margins = margins,
        a = 1,
        basis = "HS",
        selection_model = selection_model
      ),
      error = function(e) NULL
    )

    if (
      is.null(linear0) ||
      is.null(linear1) ||
      !isTRUE(linear0$success) ||
      !isTRUE(linear1$success)
    ) {
      if (attempt_count %% 10 == 0) {
        message("  Last attempt failed while saving Linear_HS weights.")
      }
      next
    }

    success_count <-
      success_count + 1

    fits <- fits %>%
      mutate(
        selection_model = selection_model,
        expected_n = expected_n,
        population_N = N,
        replicate = success_count,
        observed_n = nrow(cna),
        observed_n0 = sum(cna$A == 0),
        observed_n1 = sum(cna$A == 1),
        true_mu0 = true_mu0,
        true_mu1 = true_mu1,
        finite_mu0 = finite_mu0,
        finite_mu1 = finite_mu1
      )

    result_storage[[success_count]] <-
      fits

    linear_weight_storage[[success_count]] <-
      bind_rows(
        tibble(
          selection_model = selection_model,
          expected_n = expected_n,
          population_N = N,
          replicate = success_count,
          A = 0,
          weight =
            as.numeric(
              linear0$weights
            )
        ),
        tibble(
          selection_model = selection_model,
          expected_n = expected_n,
          population_N = N,
          replicate = success_count,
          A = 1,
          weight =
            as.numeric(
              linear1$weights
            )
        )
      ) %>%
      mutate(
        negative_weight =
          weight < 0
      )

    if (
      success_count %% 10 == 0
    ) {
      message(
        paste0(
          "Successful replications: ",
          success_count,
          "/",
          n_success,
          " | Total attempts: ",
          attempt_count,
          " | Selection: ",
          selection_model,
          " | Expected n: ",
          expected_n
        )
      )
    }
  }

  list(
    sim_results =
      bind_rows(
        result_storage
      ),
    linear_weight_results =
      bind_rows(
        linear_weight_storage
      )
  )
}


############################################################
# 12. Base selection coefficients for both designs
############################################################

lambda0_nonlinear <- c(
  X1 = 0.001,
  X2 = 0.03,
  X4 = 1.0,
  X1_X2_plus_sin_X1 = 0.02
)

lambda1_nonlinear <- c(
  X1 = 0.001,
  X2 = 0.03,
  X4 = 1.5,
  X1_X2_plus_sin_X1 = 0.035
)

# In the linear design, X4 is the shared predictor and X5 is
# the selection-only predictor. Both are Uniform(0, 1).
lambda0_linear <- c(
  X1 = 0.001,
  X2 = 0.03,
  X4 = 1.0,
  X5 = 0.75
)

lambda1_linear <- c(
  X1 = 0.001,
  X2 = 0.03,
  X4 = 1.5,
  X5 = 1.0
)


############################################################
# 12B. OPTIONAL QUICK DIAGNOSTIC BEFORE 1000 REPLICATIONS
############################################################

if (!isTRUE(getOption("simulation.functions_only", FALSE))) {
simulation_settings <- tribble(
  ~selection_model, ~expected_n, ~seed,
  "nonlinear",      400,         20260814,
  "nonlinear",      800,         20260815,
  "linear",         400,         20260816,
  "linear",         800,         20260817
)

simulation_output_list <-
  pmap(
    simulation_settings,
    function(selection_model, expected_n, seed) {

      if (selection_model == "nonlinear") {
        lambda0_current <- lambda0_nonlinear
        lambda1_current <- lambda1_nonlinear
      } else {
        lambda0_current <- lambda0_linear
        lambda1_current <- lambda1_linear
      }

      message(
        paste0(
          "Starting selection model = ",
          selection_model,
          ", expected n = ",
          expected_n
        )
      )

      run_simulation(
        n_success = getOption("simulation.n_success", 1000),
        N = 100000,
        expected_n = expected_n,
        selection_model = selection_model,
        kappa = 1,
        lambda0_base = lambda0_current,
        lambda1_base = lambda1_current,
        seed = seed
      )
    }
  )

simulation_output <- list(
  sim_results =
    map_dfr(
      simulation_output_list,
      "sim_results"
    ),
  linear_weight_results =
    map_dfr(
      simulation_output_list,
      "linear_weight_results"
    )
)


############################################################
# 14. Extract results
############################################################

sim_results <-
  simulation_output$sim_results

linear_weight_results <-
  simulation_output$linear_weight_results


############################################################
# 15. Add replication-specific errors and 95% CI coverage
############################################################

sim_results <-
  sim_results %>%
  mutate(

    true_delta =
      true_mu1 -
      true_mu0,

    delta =
      mu1 -
      mu0,

    # Estimation errors use the fixed analytic superpopulation truth.
    error_mu0 =
      mu0 -
      true_mu0,

    error_mu1 =
      mu1 -
      true_mu1,

    error_delta =
      delta -
      true_delta,

    # mu0 and mu1 estimated in separate A groups
    se_delta =
      sqrt(
        se_mu1^2 +
        se_mu0^2
      ),

    delta_lcl =
      delta -
      1.96 * se_delta,

    delta_ucl =
      delta +
      1.96 * se_delta,

    delta_covered =
      as.integer(
        true_delta >= delta_lcl &
        true_delta <= delta_ucl
      )
  )


############################################################
# 16. Verify common successful replications
############################################################

sim_results %>%
  count(
    selection_model,
    expected_n,
    scenario,
    method
  )


############################################################
# 17. Verify true means identical across scenarios
############################################################

truth_check <-
  sim_results %>%
  group_by(
    selection_model,
    expected_n,
    replicate
  ) %>%
  summarise(
    n_mu0 =
      n_distinct(
        true_mu0
      ),
    n_mu1 =
      n_distinct(
        true_mu1
      ),
    n_delta =
      n_distinct(
        true_delta
      ),
    .groups = "drop"
  )

table(
  truth_check$n_mu0
)

table(
  truth_check$n_mu1
)

table(
  truth_check$n_delta
)


############################################################
# 18. FINAL SUMMARY
############################################################

simulation_summary <-
  sim_results %>%
  group_by(
    selection_model,
    expected_n,
    scenario,
    method
  ) %>%
  summarise(

    true_mu0 =
      mean(
        true_mu0
      ),

    true_mu1 =
      mean(
        true_mu1
      ),

    # Correct empirical SE: SD of estimation errors, not
    # SD of the raw estimates when truth changes by replicate
    emp_se_mu0 =
      sd(
        error_mu0
      ),

    emp_se_mu1 =
      sd(
        error_mu1
      ),

    avg_est_se_mu0 =
      mean(
        se_mu0
      ),

    avg_est_se_mu1 =
      mean(
        se_mu1
      ),

    mean_mu0 =
      mean(
        mu0
      ),

    mean_mu1 =
      mean(
        mu1
      ),

    true_delta =
      mean(
        true_delta
      ),

    mean_delta =
      mean(
        delta
      ),

    emp_se_delta =
      sd(
        error_delta
      ),

    avg_est_se_delta =
      mean(
        se_delta
      ),

    n_covered_delta =
      sum(
        delta_covered
      ),

    n_sim =
      n(),

    coverage_rate_delta =
      mean(
        delta_covered
      ),

    .groups = "drop"
  ) %>%
  rename(
    mu0 = mean_mu0,
    mu1 = mean_mu1,
    delta = mean_delta
  ) %>%
  select(
    selection_model,
    expected_n,
    scenario,
    method,

    true_mu0,
    true_mu1,

    mu0,
    mu1,

    emp_se_mu0,
    avg_est_se_mu0,

    emp_se_mu1,
    avg_est_se_mu1,

    true_delta,
    delta,

    emp_se_delta,
    avg_est_se_delta,

    n_covered_delta,
    n_sim,
    coverage_rate_delta
  )

simulation_summary


############################################################
# 19. Bias check
############################################################

bias_check <-
  simulation_summary %>%
  mutate(

    bias_mu0 =
      mu0 -
      true_mu0,

    bias_mu1 =
      mu1 -
      true_mu1,

    bias_delta =
      delta -
      true_delta
  ) %>%
  select(
    selection_model,
    expected_n,
    scenario,
    method,
    bias_mu0,
    bias_mu1,
    bias_delta
  )

bias_check


############################################################
# 20. Linear_HS weight diagnostics
############################################################

linear_weight_summary <-
  linear_weight_results %>%
  group_by(
    selection_model,
    expected_n
  ) %>%
  summarise(
    n_weights = n(),
    n_negative = sum(weight < 0),
    proportion_negative = mean(weight < 0),
    min_weight = min(weight),
    q01_weight = quantile(weight, 0.01),
    q05_weight = quantile(weight, 0.05),
    median_weight = median(weight),
    mean_weight = mean(weight),
    q95_weight = quantile(weight, 0.95),
    q99_weight = quantile(weight, 0.99),
    max_weight = max(weight),
    .groups = "drop"
  )

linear_weight_summary


############################################################
# 21. Identify the replication with the most negative
#     Linear_HS individual weight within each design
############################################################

linear_negative_by_rep <-
  linear_weight_results %>%
  group_by(
    selection_model,
    expected_n,
    replicate
  ) %>%
  summarise(
    any_negative =
      any(weight < 0),
    n_negative =
      sum(weight < 0),
    min_weight =
      min(weight),
    .groups = "drop"
  )

linear_negative_rep_summary <-
  linear_negative_by_rep %>%
  group_by(
    selection_model,
    expected_n
  ) %>%
  summarise(
    n_replications = n(),
    n_replications_with_negative =
      sum(any_negative),
    proportion_replications_with_negative =
      mean(any_negative),
    overall_most_negative_weight =
      min(min_weight),
    .groups = "drop"
  )

linear_negative_rep_summary


# Replication containing the most negative weight within each
# selection-model-by-sample-size design
worst_linear_reps <-
  linear_negative_by_rep %>%
  group_by(
    selection_model,
    expected_n
  ) %>%
  slice_min(
    order_by = min_weight,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

worst_linear_reps


# All Linear_HS weights from the selected replication in each
# design
worst_linear_weights <-
  linear_weight_results %>%
  semi_join(
    worst_linear_reps %>%
      select(
        selection_model,
        expected_n,
        replicate
      ),
    by = c(
      "selection_model",
      "expected_n",
      "replicate"
    )
  )

worst_linear_weights %>%
  group_by(
    selection_model,
    expected_n,
    replicate
  ) %>%
  summarise(
    n_weights = n(),
    n_negative = sum(weight < 0),
    proportion_negative = mean(weight < 0),
    min_weight = min(weight),
    max_weight = max(weight),
    mean_weight = mean(weight),
    median_weight = median(weight),
    .groups = "drop"
  )


############################################################
# 22. Plot the worst replication within each design
############################################################

linear_weight_plot <-
  ggplot(
    worst_linear_weights,
    aes(
      x = weight
    )
  ) +
  geom_histogram(
    bins = 50,
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_grid(
    selection_model ~ expected_n,
    scales = "free_y",
    labeller = label_both
  ) +
  labs(
    title =
      "Distribution of Linear Calibration Weights",
    subtitle =
      "Worst replication selected separately within each design",
    x = "Normalized linear calibration weight p",
    y = "Frequency"
  ) +
  theme_classic()

linear_weight_plot


############################################################
# 23. Optional same worst replications, separated by A
############################################################

linear_weight_plot_by_A <-
  ggplot(
    worst_linear_weights,
    aes(
      x = weight
    )
  ) +
  geom_histogram(
    bins = 50,
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_grid(
    selection_model + expected_n ~ A,
    scales = "free_y",
    labeller = label_both
  ) +
  labs(
    title =
      "Distribution of Linear Calibration Weights by A Group",
    subtitle =
      "Worst replication selected separately within each design",
    x = "Normalized linear calibration weight p",
    y = "Frequency"
  ) +
  theme_classic()

linear_weight_plot_by_A


############################################################
# 24. Save outputs
############################################################

write.csv(
  sim_results,
  "simulation_replication_results.csv",
  row.names = FALSE
)

write.csv(
  simulation_summary,
  "simulation_summary.csv",
  row.names = FALSE
)

write.csv(
  linear_weight_results,
  "linear_HS_all_weights.csv",
  row.names = FALSE
)

write.csv(
  linear_negative_by_rep,
  "linear_HS_negative_weight_by_replication.csv",
  row.names = FALSE
)

write.csv(
  worst_linear_weights,
  "linear_HS_worst_replication_weights.csv",
  row.names = FALSE
)

ggsave(
  filename =
    "linear_HS_worst_replication_weight_distribution.png",
  plot =
    linear_weight_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename =
    "linear_HS_worst_replication_weight_distribution_by_A.png",
  plot =
    linear_weight_plot_by_A,
  width = 8,
  height = 6,
  dpi = 300
)

} # end optional simulation execution
