library(survey)

set.seed(20260925)

N <- 100000L
expected_n <- 400L
n_mc <- 1000L
cell_vars <- c("A", paste0("X", 1:5))

# Rows correspond to A = 0 and A = 1; covariates are independent given A.
p_x <- rbind(
  c(X1 = .40, X2 = .35, X3 = .45, X4 = .40, X5 = .50),
  c(X1 = .55, X2 = .55, X3 = .60, X4 = .55, X5 = .65)
)

# Coefficients of X1, X2, X4 in the true logistic selection model.
selection_beta <- rbind(
  c(X1 = -.30, X2 = -.50, X4 = -.70),
  c(X1 = -.60, X2 = -.80, X4 = -1.10)
)

# The nonlinear setting adds gamma[A + 1] * X1 * X4.
# IPW fits only A-specific main effects, omitting X1:X4.
selection_gamma <- list(linear = c(0, 0), nonlinear = c(.30, 1.70))

expit <- function(z) plogis(z)

# Normalized ET with d_i = 1 and no intercept in the exponent.
# The standard error projects the outcome residual on the calibrated basis.
et_mean <- function(sample_a, cells_a) {
  basis <- paste0("X", 1:5)
  target <- colSums(as.matrix(cells_a[basis]) * cells_a$N) /
    sum(cells_a$N)
  H0 <- as.matrix(sample_a[basis])
  scales <- apply(H0, 2, sd)
  scales[!is.finite(scales) | scales < 1e-8] <- 1
  H <- sweep(sweep(H0, 2, target, "-"), 2, scales, "/")

  weights_at <- function(lambda) {
    eta <- drop(H %*% lambda)
    z <- exp(eta - max(eta))
    z / sum(z)
  }
  objective <- function(lambda) {
    eta <- drop(H %*% lambda)
    max(eta) + log(sum(exp(eta - max(eta))))
  }
  gradient <- function(lambda) drop(crossprod(H, weights_at(lambda)))

  fit <- optim(rep(0, ncol(H)), objective, gradient,
               method = "BFGS", control = list(maxit = 2000, reltol = 1e-11))
  lambda <- fit$par

  for (iteration in seq_len(50L)) {
    w <- weights_at(lambda)
    g <- drop(crossprod(H, w))
    if (max(abs(g)) < 1e-9) break
    centered <- sweep(H, 2, g, "-")
    curvature <- crossprod(centered, centered * as.numeric(w))
    direction <- solve(curvature + diag(1e-12, ncol(H)), g)
    step <- 1
    improved <- FALSE
    for (j in seq_len(30L)) {
      candidate <- lambda - step * drop(direction)
      if (all(is.finite(candidate)) &&
          max(abs(gradient(candidate))) < max(abs(g))) {
        lambda <- candidate
        improved <- TRUE
        break
      }
      step <- step / 2
    }
    if (!improved) break
  }

  w <- weights_at(lambda)
  if (max(abs(drop(crossprod(H, w)))) > 1e-7)
    stop("ET calibration did not balance")
  mu <- sum(w * sample_a$Y)
  residual <- sample_a$Y - mu
  J <- crossprod(H, H * as.numeric(w))
  slope <- solve(J, crossprod(H, w * residual))
  influence <- w * (residual - drop(H %*% slope))
  c(mean = mu, se = sqrt(sum(influence^2)))
}

one_replication <- function(setting, replication, attempt) {
  cells <- expand.grid(A = 0:1, X1 = 0:1, X2 = 0:1,
                       X3 = 0:1, X4 = 0:1, X5 = 0:1)
  cells$cell <- seq_len(nrow(cells))
  cells$N <- 0L

  n_A1 <- rbinom(1L, N, .40)
  for (a in 0:1) {
    rows <- which(cells$A == a)
    x <- as.matrix(cells[rows, paste0("X", 1:5)])
    prob <- apply(x, 1L, function(row)
      prod(ifelse(row == 1, p_x[a + 1L, ], 1 - p_x[a + 1L, ])))
    cells$N[rows] <- drop(rmultinom(
      1L, if (a == 1L) n_A1 else N - n_A1, prob = prob
    ))
  }
  stopifnot(sum(cells$N) == N, all(cells$N > 0))

  cells$mu <- with(cells,
    2 * A + .6 * X1 + X2 + .8 * X3 + 2 * X4 +
      .7 * X5 + 6 * A * X4)

  cells$pi <- NA_real_
  for (a in 0:1) {
    rows <- which(cells$A == a)
    Z <- as.matrix(cells[rows, c("X1", "X2", "X4")])
    eta <- drop(Z %*% selection_beta[a + 1L, ]) +
      selection_gamma[[setting]][a + 1L] *
      cells$X1[rows] * cells$X4[rows]
    # The realized finite population has E(n | population) = expected_n.
    alpha <- uniroot(function(intercept)
      weighted.mean(expit(intercept + eta), cells$N[rows]) -
        expected_n / N, interval = c(-20, 20))$root
    cells$pi[rows] <- expit(alpha + eta)
  }
  cells$n <- rbinom(nrow(cells), cells$N, cells$pi)

  # Full poststratification needs a sampled individual in every joint cell.
  # Reject early, before generating outcomes or fitting the other methods.
  missing_cells <- sum(cells$N > 0L & cells$n == 0L)
  if (missing_cells > 0L)
    return(data.frame(empty_poststrata = missing_cells))

  sample <- cells[rep(cells$cell, cells$n),
                  c("cell", cell_vars, "mu"), drop = FALSE]
  if (nrow(sample) == 0L || length(unique(sample$A)) < 2L)
    stop("A group has no sampled individuals")
  sample$Y <- rnorm(nrow(sample), mean = sample$mu, sd = 1)

  # Generate the unsampled outcome sum to obtain the finite-population truth.
  # Estimators below receive only cells' N and the sampled records.
  sampled_sum <- numeric(nrow(cells))
  by_cell <- rowsum(sample$Y, sample$cell, reorder = FALSE)
  sampled_sum[as.integer(rownames(by_cell))] <- by_cell[, 1L]
  unsampled_error <- rnorm(nrow(cells),
                           sd = sqrt(cells$N - cells$n))
  cells$sum_Y_pop <- cells$N * cells$mu +
    sampled_sum - cells$n * cells$mu + unsampled_error
  truth_a <- vapply(0:1, function(a) {
    j <- cells$A == a
    sum(cells$sum_Y_pop[j]) / sum(cells$N[j])
  }, numeric(1))
  delta_true <- truth_a[2L] - truth_a[1L]

  # Poststratify on all 64 A x X1 x X2 x X3 x X4 x X5 cells.
  ps_vars <- cell_vars
  pop_ps <- cells[, c(ps_vars, "N"), drop = FALSE]
  sample_ps <- sample
  for (v in ps_vars) {
    pop_ps[[v]] <- factor(pop_ps[[v]], levels = 0:1)
    sample_ps[[v]] <- factor(sample_ps[[v]], levels = 0:1)
  }
  pop_table <- xtabs(N ~ A + X1 + X2 + X3 + X4 + X5, data = pop_ps)
  post_result <- tryCatch({
    design <- survey::svydesign(ids = ~1, weights = ~1,
                                data = sample_ps)
    post <- survey::postStratify(
      design, strata = ~ A + X1 + X2 + X3 + X4 + X5,
      population = pop_table
    )
    means <- survey::svyby(~Y, ~A, post, survey::svymean)
    delta <- unname(means$Y[as.character(means$A) == "1"] -
                    means$Y[as.character(means$A) == "0"])
    # Disjoint group-specific poststrata have zero cross-group covariance.
    se <- sqrt(sum(as.numeric(survey::SE(means))^2))
    c(delta = delta, se = se)
  }, error = function(e) c(delta = NA_real_, se = NA_real_))

  # Grouped conventional IPW from all 64 joint population/sample counts.
  # Omit the true X1:X4 term from the A-specific IPW propensity model.
  ipw_result <- tryCatch({
    fit <- glm(cbind(n, N - n) ~ A * (X1 + X2 + X3 + X4 + X5),
               family = binomial(link = "logit"), data = cells)
    if (!isTRUE(fit$converged)) stop("IPW did not converge")
    phat <- as.numeric(fitted(fit))
    if (any(!is.finite(phat) | phat <= 0 | phat >= 1))
      stop("Invalid IPW fitted probabilities")
    w <- 1 / phat[sample$cell]
    is1 <- sample$A == 1L
    mu1 <- weighted.mean(sample$Y[is1], w[is1])
    mu0 <- weighted.mean(sample$Y[!is1], w[!is1])
    # Joint estimating-equation linearization of the grouped binomial
    # score and the two weighted group means, with known population totals.
    Z <- model.matrix(fit)
    Z_sample <- Z[sample$cell, , drop = FALSE]
    phat_sample <- phat[sample$cell]
    direct <- numeric(nrow(sample))
    direct[is1] <- w[is1] * (sample$Y[is1] - mu1) / sum(w[is1])
    direct[!is1] <- -w[!is1] * (sample$Y[!is1] - mu0) /
      sum(w[!is1])
    J <- crossprod(Z, Z * as.numeric(cells$N * phat * (1 - phat)))
    derivative <- crossprod(Z_sample, direct * (1 - phat_sample))
    correction <- solve(J, derivative)
    score_adjustment <- drop(Z %*% correction)
    selected_influence <- direct -
      (1 - phat_sample) * score_adjustment[sample$cell]
    unsampled_influence <- phat * score_adjustment
    se <- sqrt(sum(selected_influence^2) +
               sum((cells$N - cells$n) * unsampled_influence^2))
    c(delta = mu1 - mu0, se = se)
  }, error = function(e) c(delta = NA_real_, se = NA_real_))

  # ET balances population means of X1 through X5 within each A group.
  et_result <- tryCatch({
    a1 <- et_mean(sample[sample$A == 1L, , drop = FALSE],
                  cells[cells$A == 1L, , drop = FALSE])
    a0 <- et_mean(sample[sample$A == 0L, , drop = FALSE],
                  cells[cells$A == 0L, , drop = FALSE])
    c(delta = unname(a1["mean"] - a0["mean"]),
      se = sqrt(unname(a1["se"]^2 + a0["se"]^2)))
  }, error = function(e) c(delta = NA_real_, se = NA_real_))

  data.frame(setting, replication, attempt, observed_n = nrow(sample),
             empty_poststrata = missing_cells, delta_true,
             delta_post = unname(post_result["delta"]),
             delta_ipw = unname(ipw_result["delta"]),
             delta_et = unname(et_result["delta"]),
             se_post = unname(post_result["se"]),
             se_ipw = unname(ipw_result["se"]),
             se_et = unname(et_result["se"]), row.names = NULL)
}

results <- vector("list", 2L * n_mc)
attempts_by_setting <- vector("list", 2L)
position <- 0L
for (s in seq_along(selection_gamma)) {
  setting <- names(selection_gamma)[s]
  accepted <- attempted <- rejected_empty <- rejected_fit <- 0L
  while (accepted < n_mc) {
    attempted <- attempted + 1L
    if (attempted > 250L * n_mc)
      stop("Could not obtain ", n_mc,
           " successful replications with all 64 poststrata observed for ",
           setting)
    row <- one_replication(setting, accepted + 1L, attempted)
    if (row$empty_poststrata > 0L) {
      rejected_empty <- rejected_empty + 1L
      next
    }
    estimates <- unlist(row[c("delta_true", "delta_post", "delta_ipw",
                              "delta_et", "se_post", "se_ipw", "se_et")],
                        use.names = FALSE)
    if (!all(is.finite(estimates)) || any(estimates[5:7] <= 0)) {
      rejected_fit <- rejected_fit + 1L
      next
    }
    accepted <- accepted + 1L
    position <- position + 1L
    results[[position]] <- row
    if (accepted %% 100L == 0L)
      message(setting, ": ", accepted, "/", n_mc,
              " accepted (", attempted, " attempts)")
  }
  attempts_by_setting[[s]] <- data.frame(
    setting, attempted, accepted, rejected_empty, rejected_fit,
    acceptance_rate = accepted / attempted
  )
}
mc_results <- do.call(rbind, results)
attempt_summary <- do.call(rbind, attempts_by_setting)


summary_rows <- list()
position <- 0L
for (setting in c("linear", "nonlinear")) {
  subset <- mc_results[mc_results$setting == setting, , drop = FALSE]
  draw_counts <- attempt_summary[attempt_summary$setting == setting, ]
  for (method in c("post", "ipw", "et")) {
    estimates <- subset[[paste0("delta_", method)]]
    se <- subset[[paste0("se_", method)]]
    error <- estimates - subset$delta_true
    position <- position + 1L
    summary_rows[[position]] <- data.frame(
      setting, method, N, expected_n, n_mc,
      attempted = draw_counts$attempted,
      rejected_empty = draw_counts$rejected_empty,
      rejected_fit = draw_counts$rejected_fit,
      acceptance_rate = draw_counts$acceptance_rate,
      mean_observed_n = mean(subset$observed_n),
      truth = mean(subset$delta_true),
      estimate = mean(estimates),
      bias = mean(error),
      ESE = sd(estimates),
      TSE = mean(se),
      SER = mean(se) / sd(estimates),
      coverage = mean(abs(error) <= 1.96 * se),
      RMSE = sqrt(mean(error^2)),
      MCSE_bias = sd(error) / sqrt(n_mc)
    )
  }
}
mc_summary <- do.call(rbind, summary_rows)
print(mc_summary, row.names = FALSE)

table_summary <- mc_summary[, c("setting", "expected_n", "method", "truth",
                                "estimate", "bias", "ESE", "TSE", "SER",
                                "coverage")]
table_summary <- table_summary[order(match(table_summary$setting,
                                           c("linear", "nonlinear")),
                                     match(table_summary$method,
                                           c("post", "ipw", "et"))), ]
table_summary$spec <- c("1", "2", "3", "1", "2", "3")
table_summary$method <- rep(c("Poststratification", "IPW", "ET (X1-X5)"),
                            times = 2L)
table_summary <- table_summary[, c("setting", "expected_n", "spec", "method",
                                   "truth", "estimate", "bias", "ESE", "TSE",
                                   "SER", "coverage")]
names(table_summary) <- c("Selection", "E(n)", "Spec.", "Estimator", "Truth",
                          "Estimate", "Bias", "ESE", "TSE", "SER", "Coverage")
print(table_summary, row.names = FALSE, digits = 4)

write.csv(mc_results, "joint64cells_post64_X1X4_omitted_binomial_ipw_all5_replications.csv",
          row.names = FALSE)
write.csv(mc_summary, "joint64cells_post64_X1X4_omitted_binomial_ipw_all5_summary.csv",
          row.names = FALSE)
write.csv(table_summary, "joint64cells_post64_X1X4_omitted_binomial_ipw_all5_table.csv",
          row.names = FALSE)
write.csv(attempt_summary, "joint64cells_post64_X1X4_omitted_binomial_ipw_all5_attempts.csv",
          row.names = FALSE)
