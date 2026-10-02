cross_entropy_weights <- function(
    data,
    benchmarks,
    group = NULL,
    base_weights = NULL,
    normalize = c("mean1", "sum1"),
    na_action = c("error", "omit"),
    weight_name = "ET_weight",
    maxit = 5000,
    balance_tol = 1e-6
) {
  normalize <- match.arg(normalize)
  na_action <- match.arg(na_action)
  if (!is.data.frame(data) || nrow(data) < 1L) {
    stop("data must be a nonempty data frame.", call. = FALSE)
  }
  if (anyDuplicated(names(data))) {
    stop("Sample column names must be unique.", call. = FALSE)
  }
  if (!is.character(weight_name) || length(weight_name) != 1L ||
      is.na(weight_name) || !nzchar(weight_name) || weight_name %in% names(data)) {
    stop("weight_name must name a new output column.", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
      maxit < 1 || maxit != floor(maxit)) {
    stop("maxit must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(balance_tol) || length(balance_tol) != 1L ||
      !is.finite(balance_tol) || balance_tol <= 0) {
    stop("balance_tol must be a positive finite number.", call. = FALSE)
  }
  if (!is.null(group) && (!is.character(group) || length(group) != 1L ||
      is.na(group) || !group %in% names(data) ||
      group %in% c("variable", "category", "target"))) {
    stop("group must name a sample column other than variable/category/target.",
         call. = FALSE)
  }
  if (!is.null(group) && (!is.atomic(data[[group]]) || is.matrix(data[[group]]))) {
    stop("The grouping column must be a vector.", call. = FALSE)
  }

  if (!is.data.frame(benchmarks) && is.list(benchmarks)) {
    if (!is.null(group)) {
      stop("Grouped calibration requires a benchmark data frame.", call. = FALSE)
    }
    if (length(benchmarks) == 0L || is.null(names(benchmarks)) ||
        anyNA(names(benchmarks)) || any(!nzchar(names(benchmarks))) ||
        anyDuplicated(names(benchmarks))) {
      stop("The benchmark list must have unique variable names.", call. = FALSE)
    }
    benchmarks <- do.call(rbind, lapply(names(benchmarks), function(v) {
      z <- benchmarks[[v]]
      if (!is.numeric(z) || is.complex(z) || length(z) == 0L) {
        stop("Each benchmark must be a numeric mean or named proportions: ",
             v, call. = FALSE)
      }
      if (is.null(names(z))) {
        if (length(z) != 1L) {
          stop("An unnamed benchmark must be one numeric mean: ", v,
               call. = FALSE)
        }
        category <- NA_character_
      } else {
        if (anyNA(names(z)) || any(!nzchar(names(z)))) {
          stop("Categorical benchmarks need category names: ", v,
               call. = FALSE)
        }
        category <- names(z)
      }
      data.frame(variable = v, category = category, target = as.numeric(z),
                 stringsAsFactors = FALSE)
    }))
  }
  if (!is.data.frame(benchmarks) || nrow(benchmarks) < 1L ||
      anyDuplicated(names(benchmarks)) ||
      !all(c("variable", "target") %in% names(benchmarks))) {
    stop("benchmarks must have variable and target columns.", call. = FALSE)
  }
  if (!is.numeric(benchmarks$target) || is.complex(benchmarks$target) ||
      any(!is.finite(benchmarks$target))) {
    stop("Benchmark targets must be finite numbers.", call. = FALSE)
  }
  if (!"category" %in% names(benchmarks)) benchmarks$category <- NA_character_
  benchmarks$variable <- as.character(benchmarks$variable)
  benchmarks$category <- as.character(benchmarks$category)
  if (anyNA(benchmarks$variable) || any(!nzchar(benchmarks$variable)) ||
      !all(benchmarks$variable %in% names(data))) {
    stop("Every benchmark variable must name a sample column.", call. = FALSE)
  }
  if (!is.null(group) && (!group %in% names(benchmarks) ||
      anyNA(benchmarks[[group]]))) {
    stop("The benchmark table must contain nonmissing group labels in ", group,
         ".", call. = FALSE)
  }
  if (!is.null(group) && group %in% benchmarks$variable) {
    stop("The grouping column cannot also be a calibration variable.",
         call. = FALSE)
  }

  n <- nrow(data)
  d <- base_weights
  if (is.null(d)) d <- rep(1, n)
  if (is.character(d) && length(d) == 1L && !is.na(d) && d %in% names(data)) {
    d <- data[[d]]
  }
  if (!is.numeric(d) || is.complex(d) || length(d) != n ||
      any(!is.finite(d)) || any(d <= 0)) {
    stop("base_weights must supply one positive finite weight per sample row.",
         call. = FALSE)
  }

  variables <- unique(benchmarks$variable)
  missing <- rep(FALSE, n)
  for (v in variables) {
    if (!is.atomic(data[[v]]) || is.matrix(data[[v]])) {
      stop("Calibration variables must be vectors: ", v, call. = FALSE)
    }
    missing <- missing | is.na(data[[v]])
    if (any(is.na(benchmarks$category[benchmarks$variable == v]))) {
      if (!is.numeric(data[[v]]) || is.complex(data[[v]])) {
        stop("A numeric-mean benchmark requires a numeric sample column: ", v,
             call. = FALSE)
      }
      missing <- missing | !is.finite(data[[v]])
    }
  }
  if (!is.null(group)) missing <- missing | is.na(data[[group]])
  if (any(missing) && na_action == "error") {
    stop("Missing or nonfinite calibration/group values in ", sum(missing),
         " rows. Clean them or set na_action = 'omit'.", call. = FALSE)
  }
  included <- which(!missing)
  if (length(included) == 0L) {
    stop("No complete sample rows remain.", call. = FALSE)
  }
  sample_group <- if (is.null(group)) rep("All", n) else as.character(data[[group]])
  benchmark_group <- if (is.null(group)) rep("All", nrow(benchmarks)) else
    as.character(benchmarks[[group]])
  groups <- unique(sample_group[included])
  if (!setequal(groups, unique(benchmark_group))) {
    stop("Sample groups and benchmark groups must match after exclusions.",
         call. = FALSE)
  }

  fit_one_group <- function(rows, b, label) {
    x <- data[rows, , drop = FALSE]
    H <- matrix(numeric(0), nrow = length(rows), ncol = 0L)
    targets <- numeric(0)
    column_labels <- character(0)
    report <- vector("list", length(unique(b$variable)))
    for (j in seq_along(unique(b$variable))) {
      v <- unique(b$variable)[j]
      bv <- b[b$variable == v, , drop = FALSE]
      if (all(is.na(bv$category))) {
        if (nrow(bv) != 1L) stop("Duplicate numeric benchmark for ", v,
                                    call. = FALSE)
        h <- matrix(x[[v]], ncol = 1L)
        t <- bv$target
        lo <- min(h)
        hi <- max(h)
        if (lo < hi && (t <= lo || t >= hi)) {
          stop("For positive weights, the mean target for ", v,
               " must be inside its sample range (group ", label, ").",
               call. = FALSE)
        }
        labels <- v
      } else {
        if (anyNA(bv$category) || anyDuplicated(bv$category)) {
          stop("Use either one mean or unique category proportions for ", v,
               call. = FALSE)
        }
        if (any(bv$target < 0 | bv$target > 1) ||
            abs(sum(bv$target) - 1) > 1e-8) {
          stop("Category targets for ", v,
               " must be proportions from 0 to 1 summing to 1.", call. = FALSE)
        }
        values <- as.character(x[[v]])
        if (!all(values %in% bv$category)) {
          stop("Sample category labels are missing from benchmarks for ", v,
               " (group ", label, ").", call. = FALSE)
        }
        h <- vapply(bv$category, function(z) as.numeric(values == z),
                    numeric(nrow(x)))
        h <- matrix(h, nrow = nrow(x), ncol = nrow(bv))
        t <- bv$target
        observed <- colSums(h)
        if (any(observed == 0 & t > 0) || any(observed > 0 & t == 0)) {
          stop("Category support is incompatible with positive weights for ",
               v, " (group ", label, ").", call. = FALSE)
        }
        labels <- paste(v, bv$category, sep = "::")
      }
      H <- cbind(H, h)
      targets <- c(targets, t)
      column_labels <- c(column_labels, labels)
      report[[j]] <- data.frame(group = label, variable = v,
                               category = bv$category, target = t,
                               stringsAsFactors = FALSE)
    }
    colnames(H) <- make.unique(column_labels)
    report <- do.call(rbind, report)
    sample_mean <- colMeans(H)
    scales <- apply(H, 2L, sd)
    scales[!is.finite(scales) | scales < 1e-8] <- 1
    D <- sweep(sweep(H, 2L, sample_mean, "-"), 2L, scales, "/")
    rank_fit <- qr(D, tol = 1e-10)
    k <- rank_fit$rank
    keep <- if (k > 0L) rank_fit$pivot[seq_len(k)] else integer(0)
    target_shift <- (targets - sample_mean) / scales
    implied <- rep(0, ncol(H))
    if (k > 0L) {
      relations <- qr.coef(qr(D[, keep, drop = FALSE], tol = 1e-10), D)
      implied <- as.vector(target_shift[keep] %*% relations)
    }
    if (any(abs((implied - target_shift) * scales) > balance_tol)) {
      stop("Benchmarks conflict with constant or redundant sample variables",
           " (group ", label, ").", call. = FALSE)
    }

    # Remove redundant constraints; verify every original benchmark below.
    G <- sweep(sweep(H[, keep, drop = FALSE], 2L, targets[keep], "-"),
               2L, scales[keep], "/")
    log_base <- log(d[rows])
    log_base <- log_base - max(log_base)
    evaluate <- function(lambda) {
      z <- as.vector(G %*% lambda) + log_base
      if (any(!is.finite(z))) return(NULL)
      zmax <- max(z)
      u <- exp(z - zmax)
      p <- u / sum(u)
      list(value = zmax + log(sum(u)), p = p,
           gradient = as.vector(crossprod(G, p)))
    }
    lambda <- rep(0, k)
    convergence <- 0L
    newton_steps <- 0L
    if (k > 0L) {
      fit <- tryCatch(optim(
        par = lambda,
        fn = function(z) { e <- evaluate(z); if (is.null(e)) 1e100 else e$value },
        gr = function(z) { e <- evaluate(z); if (is.null(e)) rep(1e100, k) else e$gradient },
        method = "BFGS", control = list(maxit = maxit, reltol = 1e-12)
      ), error = function(e) stop("ET optimization failed (group ", label,
                                  "): ", conditionMessage(e), call. = FALSE))
      lambda <- fit$par
      convergence <- fit$convergence
      for (iteration in seq_len(50L)) {
        e <- evaluate(lambda)
        if (is.null(e) || max(abs(e$gradient)) < 1e-12) break
        deviation <- sweep(G, 2L, e$gradient, "-")
        hessian <- crossprod(deviation, deviation * e$p)
        direction <- tryCatch(solve(hessian + diag(1e-12, k), e$gradient),
                              error = function(e) NULL)
        if (is.null(direction) || any(!is.finite(direction))) break
        accepted <- FALSE
        for (step in 0:29) {
          candidate <- lambda - direction * 2^(-step)
          trial <- evaluate(candidate)
          if (!is.null(trial) && max(abs(trial$gradient)) < max(abs(e$gradient))) {
            lambda <- candidate
            newton_steps <- newton_steps + 1L
            accepted <- TRUE
            break
          }
        }
        if (!accepted) break
      }
    }
    e <- evaluate(lambda)
    if (is.null(e) || any(!is.finite(e$p)) || any(e$p <= 0)) {
      stop("Invalid ET weights (group ", label, ").", call. = FALSE)
    }
    report$unweighted <- sample_mean
    report$weighted <- as.vector(crossprod(H, e$p))
    report$difference <- report$weighted - report$target
    error <- max(abs(report$difference))
    if (!is.finite(error) || error > balance_tol) {
      stop("ET could not match benchmarks in group ", label,
           "; maximum error = ", signif(error, 6),
           ". Check joint feasibility and sample support.", call. = FALSE)
    }
    weights <- e$p * if (normalize == "mean1") length(rows) else 1
    coefficients <- setNames(rep(0, ncol(H)), colnames(H))
    coefficients[keep] <- lambda / scales[keep]
    list(weights = weights, balance = report, coefficients = coefficients,
         diagnostics = data.frame(
           group = label, n = length(rows), weight_sum = sum(weights),
           min_weight = min(weights), max_weight = max(weights),
           ess = 1 / sum(e$p^2), balance_error = error,
           independent_constraints = k, optim_convergence = convergence,
           newton_steps = newton_steps,
           stringsAsFactors = FALSE
         ))
  }

  weights <- rep(NA_real_, n)
  fits <- setNames(vector("list", length(groups)), groups)
  for (i in seq_along(groups)) {
    g <- groups[i]
    rows <- included[sample_group[included] == g]
    fits[[i]] <- fit_one_group(rows, benchmarks[benchmark_group == g, , drop = FALSE], g)
    weights[rows] <- fits[[i]]$weights
  }
  output <- data
  output[[weight_name]] <- weights
  balance <- do.call(rbind, lapply(fits, function(z) z$balance))
  diagnostics <- do.call(rbind, lapply(fits, function(z) z$diagnostics))
  rownames(balance) <- NULL
  rownames(diagnostics) <- NULL
  list(data = output, weights = weights, balance = balance,
       diagnostics = diagnostics,
       coefficients = lapply(fits, function(z) z$coefficients),
       excluded_rows = which(missing), normalize = normalize,
       group_variable = group, call = match.call())
}

# Example: 

sample_data <- data.frame(
    Age = c(25, 35, 45, 55, 65, 40),
    Sex = c("Female", "Male", "Female", "Male", "Female", "Male"),
    Outcome = c(2, 3, 1, 4, 2, 3)
    )

benchmarks <- list(Age = 45, Sex = c(Female = 0.55, Male = 0.45))

fit <- cross_entropy_weights(sample_data, benchmarks)

fit$data
fit$balance
fit$diagnostics
weighted.mean(fit$data$Outcome, fit$data$ET_weight)
  
