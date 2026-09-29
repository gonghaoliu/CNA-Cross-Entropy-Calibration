############################################################

library(dplyr)
library(tibble)
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

    # Estimation basis only: keep the true selection predictors in
    # get_selection_predictors() and generate_sample() unchanged.
    H <- if (selection_model == "nonlinear") {
      cbind(X1 = X1, X2 = X2, X4 = X4)
    } else {
      get_selection_predictors(data, selection_model)
    }

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
    # Outcome-only variable, dependent on A but not the other predictors.
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
  
  p_X2 <- expit(gamma20 + gamma21 * A)
  
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
  # Outcome-only variable, dependent on A but not the other predictors.
  ##########################################################

  X3 <- rnorm(
    N,
    mean = X3_mean0 + X3_A_effect * A,
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
# 4. Generate CNA / nonprobability sample
############################################################

generate_sample <- function(
    pop,
    expected_n = 400,
    selection_model = c("nonlinear", "linear"),
    kappa = 3,
    lambda0_base = NULL,
    lambda1_base = NULL
) {

  selection_model <- match.arg(selection_model)

  # X4 is the shared selection-outcome variable. The combined
  # nonlinear term is selection-only in the nonlinear design,
  # and X5 is selection-only in the linear design.
  if (is.null(lambda0_base)) {
    lambda0_base <-
      if (selection_model == "nonlinear") {
        c(
          X1 = 0.001,
          X2 = -1.0,
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
          X2 = -1.75,
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
    
    n_expected_a <-
      expected_n *
      N_a / N
    
    target_fraction <-
      n_expected_a /
      N_a
    
    lambda_base <-
      if (a == 0) {
        
        lambda0_base
        
      } else {
        
        lambda1_base
        
      }
    
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
        X4 = z$N_a * z$mean_X4
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
  et_se_method <- match.arg(
    getOption("simulation.et_se_method", "HC3"),
    c("HC3", "HC0")
  )

  cna_a <- cna %>% filter(A == a)
  n_a <- nrow(cna_a)

  fail <- function(reason) {
    list(success = FALSE, reason = reason)
  }

  if (n_a < 5) {
    return(fail("n_a < 5"))
  }

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
  var_mu_hc0 <- sandwich[k + 1, k + 1]

  if (!is.finite(var_mu_hc0) || var_mu_hc0 < 0) {
    return(fail("invalid ET variance"))
  }

  H_reg <- cbind(`(Intercept)` = 1, H)
  H_weighted <- H_reg * as.numeric(p_norm)
  gram_inv <- tryCatch(
    solve(crossprod(H_reg, H_weighted)),
    error = function(e) NULL
  )
  if (is.null(gram_inv)) {
    return(fail("singular ET HC3 weighted-regression matrix"))
  }
  beta_weighted <-
    as.vector(gram_inv %*% crossprod(H_reg, p_norm * Y))
  residual <- Y - as.vector(H_reg %*% beta_weighted)
  leverage <-
    p_norm * rowSums((H_reg %*% gram_inv) * H_reg)
  leverage_denom <- 1 - leverage
  if (any(!is.finite(leverage_denom)) ||
      any(leverage_denom <= 1e-8)) {
    return(fail("ET HC3 leverage too close to one"))
  }
  var_mu_hc3 <-
    sum((p_norm * residual / leverage_denom)^2)
  if (!is.finite(var_mu_hc3) || var_mu_hc3 < 0) {
    return(fail("invalid ET HC3 variance"))
  }

  list(
    success = TRUE,
    reason = "ok",
    mu = mu_hat,
    se = sqrt(if (et_se_method == "HC3") var_mu_hc3 else var_mu_hc0),
    se_hc0 = sqrt(var_mu_hc0),
    se_hc3 = sqrt(var_mu_hc3),
    max_leverage = max(leverage),
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
    se_mu0_hc0 = fit0$se_hc0,
    se_mu1_hc0 = fit1$se_hc0,
    se_mu0_hc3 = fit0$se_hc3,
    se_mu1_hc3 = fit1$se_hc3,
    max_leverage0 = fit0$max_leverage,
    max_leverage1 = fit1$max_leverage,
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
# 9. Optional diagnostic helper
############################################################

diagnose_one_dataset <- function(cna, margins, selection_model) {
  list(
    ET_HS_A0 = fit_entropy_ee_group(cna, margins, 0, "HS", selection_model),
    ET_HS_A1 = fit_entropy_ee_group(cna, margins, 1, "HS", selection_model),
    ET_HY_A0 = fit_entropy_ee_group(cna, margins, 0, "HY", selection_model),
    ET_HY_A1 = fit_entropy_ee_group(cna, margins, 1, "HY", selection_model),
    ET_HF_A0 = fit_entropy_ee_group(cna, margins, 0, "HF", selection_model),
    ET_HF_A1 = fit_entropy_ee_group(cna, margins, 1, "HF", selection_model),
    Linear_HS_A0 = fit_linear_group(cna, margins, 0, "HS", selection_model),
    Linear_HS_A1 = fit_linear_group(cna, margins, 1, "HS", selection_model),
    ET_H0_A0 = fit_entropy_ee_group(cna, margins, 0, "H0", selection_model),
    ET_H0_A1 = fit_entropy_ee_group(cna, margins, 1, "H0", selection_model)
  )
}

run_all_scenarios <- function(
    cna,
    margins,
    selection_model = c("nonlinear", "linear")
) {
  selection_model <- match.arg(selection_model)
  specs <- tibble(
    scenario = c("S1", "S4"),
    method = c("ET_HS", "Linear_HS"),
    basis = c("HS", "HS")
  )
  results <- vector("list", nrow(specs))
  status <- vector("list", nrow(specs))

  for (j in seq_len(nrow(specs))) {
    fit <- tryCatch(
      if (specs$method[[j]] == "Linear_HS") {
        fit_linear(cna, margins, "HS", selection_model)
      } else {
        fit_entropy_ee(cna, margins, specs$basis[[j]], selection_model)
      },
      error = function(e) e
    )
    success <- !inherits(fit, "error") && !is.null(fit) &&
      nrow(fit) == 1 && isTRUE(fit$success[[1]])

    status[[j]] <- tibble(
      scenario = specs$scenario[[j]],
      method = specs$method[[j]],
      success = success,
      reason = if (success) "ok" else if (inherits(fit, "error")) {
        conditionMessage(fit)
      } else {
        "estimator fit failed"
      }
    )
    if (success) {
      results[[j]] <- fit %>% mutate(
        scenario = specs$scenario[[j]],
        method = specs$method[[j]]
      )
    }
  }

  list(results = bind_rows(results), status = bind_rows(status))
}


############################################################
# 11. Main simulation: 1000 accepted populations per setting
############################################################

run_simulation <- function(
    n_success = 1000,
    N = 100000,
    simulation_settings = tibble(
      selection_model = c("nonlinear", "nonlinear", "linear", "linear"),
      expected_n = c(400, 800, 400, 800)
    ),
    kappa = 3,
    seed = 20260814,
    max_attempts = 100L * n_success
) {
  if (length(n_success) != 1L || !is.finite(n_success) ||
      n_success < 1 || n_success != floor(n_success)) {
    stop("n_success must be a positive integer")
  }
  if (length(max_attempts) != 1L || !is.finite(max_attempts) ||
      max_attempts < n_success || max_attempts != floor(max_attempts)) {
    stop("max_attempts must be an integer at least n_success")
  }
  n_settings <- nrow(simulation_settings)
  if (n_settings < 1L) stop("simulation_settings must contain at least one design")
  result_storage <- vector("list", n_success * n_settings)
  linear_weight_storage <- vector("list", n_success * n_settings)
  population_reference_storage <- vector("list", n_success * n_settings)
  attempt_storage <- vector("list", n_settings)

  for (s in seq_len(n_settings)) {
    design <- simulation_settings$selection_model[[s]]
    expected_n <- simulation_settings$expected_n[[s]]
    lambda0 <- if (design == "nonlinear") lambda0_nonlinear else lambda0_linear
    lambda1 <- if (design == "nonlinear") lambda1_nonlinear else lambda1_linear
    accepted <- 0L
    attempts <- 0L
    setting_attempts <- vector("list", max_attempts)
    message("Starting ", design, ", expected n = ", expected_n,
            ": ", n_success, " successful populations; N = ", N)

    while (accepted < n_success) {
      attempts <- attempts + 1L
      if (attempts > max_attempts) {
        stop("Only ", accepted, " replications succeeded for ", design,
             ", expected n = ", expected_n, " in ", max_attempts, " attempts")
      }
      set.seed(seed + (s - 1L) * max_attempts + attempts - 1L)
      pop <- generate_population(N = N)
      pop_mu0 <- mean(pop$Y[pop$A == 0])
      pop_mu1 <- mean(pop$Y[pop$A == 1])
      margins <- make_marginal_targets(pop)
      if (attempts == 1L) message("Attempt 1: population generated; fitting methods")

      sample_object <- tryCatch(
        generate_sample(pop, expected_n, design, kappa, lambda0, lambda1),
        error = function(e) e
      )
      failure <- NULL
      if (inherits(sample_object, "error")) {
        failure <- conditionMessage(sample_object)
      } else {
        cna <- sample_object$cna
        rm(sample_object)
        if (length(unique(cna$A)) < 2) {
          failure <- "one A group absent"
        } else {
          fitted <- run_all_scenarios(cna, margins, design)
          if (!all(fitted$status$success)) {
            failure <- paste(
              fitted$status$method[!fitted$status$success],
              fitted$status$reason[!fitted$status$success], sep = ": ",
              collapse = "; "
            )
          }
        }
      }

      if (!is.null(failure)) {
        setting_attempts[[attempts]] <- tibble(
          selection_model = design, expected_n = expected_n,
          attempt = attempts, accepted = FALSE, reason = failure
        )
        message(design, ", n = ", expected_n, ", attempt ", attempts,
                " rejected; accepted ", accepted, "/", n_success,
                "; reason: ", failure)
        next
      }

      accepted <- accepted + 1L
      idx <- (s - 1L) * n_success + accepted
      population_reference_storage[[idx]] <- tibble(
        selection_model = design, expected_n = expected_n,
        replicate = accepted, attempt = attempts, population_N = N,
        pop_mu0 = pop_mu0, pop_mu1 = pop_mu1
      )
      setting_attempts[[attempts]] <- tibble(
        selection_model = design, expected_n = expected_n,
        attempt = attempts, accepted = TRUE, reason = "ok"
      )
      result_storage[[idx]] <- fitted$results %>% mutate(
        selection_model = design, expected_n = expected_n,
        population_N = N, replicate = accepted,
        observed_n = nrow(cna),
        observed_n0 = sum(cna$A == 0),
        observed_n1 = sum(cna$A == 1)
      )

      linear0 <- tryCatch(fit_linear_group(cna, margins, 0, "HS", design),
                          error = function(e) NULL)
      linear1 <- tryCatch(fit_linear_group(cna, margins, 1, "HS", design),
                          error = function(e) NULL)
      if (!is.null(linear0) && !is.null(linear1) &&
          isTRUE(linear0$success) && isTRUE(linear1$success)) {
        linear_weight_storage[[idx]] <- bind_rows(
          tibble(selection_model = design, expected_n = expected_n,
                 population_N = N, replicate = accepted, A = 0,
                 weight = as.numeric(linear0$weights)),
          tibble(selection_model = design, expected_n = expected_n,
                 population_N = N, replicate = accepted, A = 1,
                 weight = as.numeric(linear1$weights))
        ) %>% mutate(negative_weight = weight < 0)
      }
      message(design, ", n = ", expected_n, ", attempt ", attempts,
              " accepted; accepted ", accepted, "/", n_success)
    }
    attempt_storage[[s]] <- bind_rows(setting_attempts[seq_len(attempts)])
    message("Finished ", design, ", n = ", expected_n,
            ": ", accepted, " accepted in ", attempts, " attempts")
  }
  message("Finished: ", n_success * n_settings, " accepted populations total")

  list(
    sim_results = bind_rows(result_storage),
    linear_weight_results = bind_rows(linear_weight_storage),
    population_reference_results = bind_rows(population_reference_storage),
    attempt_log = bind_rows(attempt_storage)
  )
}


############################################################
# 12. Base selection coefficients for both designs
############################################################

# The nonlinear X2 coefficients offset X1*X2 at X1 around 50.
lambda0_nonlinear <- c(
  X1 = 0.001,
  X2 = -1.0,
  X4 = 1.0,
  X1_X2_plus_sin_X1 = 0.02
)

lambda1_nonlinear <- c(
  X1 = 0.001,
  X2 = -1.75,
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
kappa_current <- getOption("simulation.kappa", 3)
if (!is.numeric(kappa_current) || length(kappa_current) != 1 ||
    !is.finite(kappa_current) || kappa_current <= 0) {
  stop("simulation.kappa must be one positive finite number")
}
population_N_current <- getOption("simulation.population_N", 100000)
et_se_method_current <- match.arg(
  getOption("simulation.et_se_method", "HC3"),
  c("HC3", "HC0")
)
if (!is.numeric(population_N_current) || length(population_N_current) != 1 ||
    !is.finite(population_N_current) || population_N_current < 2 ||
    population_N_current != floor(population_N_current)) {
  stop("simulation.population_N must be an integer greater than 1")
}
run_tag <- paste0("kappa", gsub("\\.", "_", as.character(kappa_current)),
                  "_N", population_N_current,
                  "_nonlinearHS_X1_X2_X4_S1_S4_ETse", et_se_method_current)
simulation_settings <- tribble(
  ~selection_model, ~expected_n,
  "nonlinear",      400,
  "nonlinear",      800
)

simulation_output <- run_simulation(
  n_success = getOption("simulation.n_success", 1000),
  N = population_N_current,
  simulation_settings = simulation_settings,
  kappa = kappa_current,
  seed = 20260814,
  max_attempts = getOption("simulation.max_attempts",
                           100L * getOption("simulation.n_success", 1000))
)


############################################################
# 14. Extract results
############################################################

sim_results <-
  simulation_output$sim_results

linear_weight_results <-
  simulation_output$linear_weight_results

population_reference_results <-
  simulation_output$population_reference_results

attempt_log <- simulation_output$attempt_log

mc_truth <- tibble(
  true_mu0 = mean(population_reference_results$pop_mu0),
  true_mu1 = mean(population_reference_results$pop_mu1),
  n_reference_populations = nrow(population_reference_results)
) %>%
  mutate(true_delta = true_mu1 - true_mu0)

mc_truth


############################################################
# 15. Add replication-specific errors and 95% CI coverage
############################################################

sim_results <-
  sim_results %>%
  mutate(
    true_mu0 = mc_truth$true_mu0[[1]],
    true_mu1 = mc_truth$true_mu1[[1]],

    true_delta =
      true_mu1 -
      true_mu0,

    delta =
      mu1 -
      mu0,

    # All replications use the same Monte Carlo averaged truth.
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

    # Keep the previous uncorrected ET sandwich interval as a diagnostic.
    # It is NA for linear calibration, which has separate variance code.
    se_delta_hc0 =
      sqrt(se_mu1_hc0^2 + se_mu0_hc0^2),
    delta_covered_hc0 =
      as.integer(abs(error_delta) <= 1.96 * se_delta_hc0),

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
# 16. Verify 1000 accepted replications per setting and method
############################################################

replication_counts <- sim_results %>%
  count(
    selection_model,
    expected_n,
    scenario,
    method
  )

reference_counts <- population_reference_results %>%
  count(selection_model, expected_n)

stopifnot(
  nrow(population_reference_results) ==
    getOption("simulation.n_success", 1000) * nrow(simulation_settings),
  nrow(reference_counts) == nrow(simulation_settings),
  all(reference_counts$n == getOption("simulation.n_success", 1000)),
  nrow(replication_counts) == nrow(simulation_settings) * 2L,
  all(replication_counts$n == getOption("simulation.n_success", 1000))
)

replication_counts

attempt_summary <- attempt_log %>%
  group_by(selection_model, expected_n) %>%
  summarise(n_attempted = n(), n_accepted = sum(accepted),
            n_rejected = sum(!accepted), .groups = "drop")

print(attempt_summary)


############################################################
# 17. Verify true means identical across all settings and scenarios
############################################################

truth_check <-
  sim_results %>%
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
      )
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

    # The reference truth is fixed, so SD(errors) = SD(estimates).
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

    avg_est_se_delta_hc0 =
      if (all(is.na(se_delta_hc0))) NA_real_
      else mean(se_delta_hc0, na.rm = TRUE),

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

    coverage_rate_delta_hc0 =
      if (all(is.na(delta_covered_hc0))) NA_real_
      else mean(delta_covered_hc0, na.rm = TRUE),

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
    coverage_rate_delta,
    avg_est_se_delta_hc0,
    coverage_rate_delta_hc0
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
  paste0("simulation_replication_results_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  simulation_summary,
  paste0("simulation_summary_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  mc_truth,
  paste0("simulation_mc_truth_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  population_reference_results,
  paste0("simulation_population_reference_means_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  attempt_log,
  paste0("simulation_attempt_log_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  linear_weight_results,
  paste0("linear_HS_all_weights_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  linear_negative_by_rep,
  paste0("linear_HS_negative_weight_by_replication_", run_tag, ".csv"),
  row.names = FALSE
)

write.csv(
  worst_linear_weights,
  paste0("linear_HS_worst_replication_weights_", run_tag, ".csv"),
  row.names = FALSE
)

ggsave(
  filename =
    paste0("linear_HS_worst_replication_weight_distribution_", run_tag, ".png"),
  plot =
    linear_weight_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename =
    paste0("linear_HS_worst_replication_weight_distribution_by_A_", run_tag, ".png"),
  plot =
    linear_weight_plot_by_A,
  width = 8,
  height = 6,
  dpi = 300
)

} # end optional simulation execution
