### Australia-wide comparison: standard CCEMG versus ABC whole-sum DWB, l_N = 1.
### Section 5 wheat-yield model; all 119 regions are estimated as one panel.
### JCK is deliberately excluded.

rm(list = ls())

required_packages <- c("MASS", "readxl", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()

variables <- c("FD", "DDlow", "DDmod", "DDhigh", "RF", "FD_RF", "DDhigh_RF")
variable_labels <- c(FD = "FD", DDlow = "DDlow", DDmod = "DDmod",
                     DDhigh = "DDhigh", RF = "RF", FD_RF = "FD x RF",
                     DDhigh_RF = "DDhigh x RF")
sample_label <- "Australia"
abc_method <- "ABC + whole-sum DWB"

read_integer <- function(name, default, minimum = 0L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < minimum) default else value
}

read_number <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || !is.finite(value)) default else value
}

matrix_pinv <- function(x) MASS::ginv(x, tol = 1e-10)

fit_cce <- function(Y, X, rows = seq_len(nrow(Y))) {
  Y <- Y[rows, , drop = FALSE]
  X <- X[rows, , , drop = FALSE]
  Tn <- nrow(Y)
  N <- ncol(Y)
  k <- dim(X)[3L]
  y_bar <- rowMeans(Y)
  x_bar <- apply(X, c(1L, 3L), mean)
  controls <- cbind(1, y_bar, x_bar)

  betas <- matrix(NA_real_, N, k)
  x_tilde <- array(NA_real_, c(Tn, N, k))
  residuals <- matrix(NA_real_, Tn, N)
  gram <- array(NA_real_, c(k, k, N))
  intercept_y <- numeric(N)
  proxy_y <- matrix(NA_real_, N, k + 1L)
  intercept_x <- matrix(NA_real_, N, k)
  proxy_x <- array(NA_real_, c(N, k + 1L, k))
  residual_se <- numeric(N)
  joint_rank <- integer(N)

  control_inverse <- matrix_pinv(crossprod(controls))
  for (i in seq_len(N)) {
    Xi <- X[, i, , drop = FALSE][, 1L, ]
    yi <- Y[, i]
    coefficient_x <- control_inverse %*% crossprod(controls, Xi)
    coefficient_y <- control_inverse %*% crossprod(controls, yi)
    Xi_res <- Xi - controls %*% coefficient_x
    yi_res <- yi - controls %*% coefficient_y
    beta_i <- matrix_pinv(crossprod(Xi_res)) %*% crossprod(Xi_res, yi_res)
    nuisance_y <- control_inverse %*% crossprod(controls, yi - Xi %*% beta_i)
    error_i <- as.numeric(yi - Xi %*% beta_i - controls %*% nuisance_y)
    rank_i <- qr(cbind(controls, Xi), tol = 1e-10)$rank
    degrees_freedom <- max(1L, Tn - rank_i)

    betas[i, ] <- beta_i
    x_tilde[, i, ] <- Xi_res
    residuals[, i] <- error_i
    gram[, , i] <- crossprod(Xi_res) / Tn
    intercept_y[i] <- nuisance_y[1L]
    proxy_y[i, ] <- nuisance_y[-1L]
    intercept_x[i, ] <- coefficient_x[1L, ]
    proxy_x[i, , ] <- coefficient_x[-1L, , drop = FALSE]
    residual_se[i] <- sqrt(sum(error_i^2) / degrees_freedom)
    joint_rank[i] <- rank_i
  }

  list(betas = betas, x_tilde = x_tilde, residuals = residuals,
       gram = gram, intercept_y = intercept_y, proxy_y = proxy_y,
       intercept_x = intercept_x, proxy_x = proxy_x,
       residual_se = residual_se, y_bar = y_bar, x_bar = x_bar,
       joint_rank = joint_rank, N = N, T = Tn, k = k)
}

analytical_bias_correction <- function(fit, Y, X, factor_count = 2L) {
  N <- fit$N
  Tn <- fit$T
  k <- fit$k
  if (factor_count < 1L || factor_count > min(k + 1L, Tn - 1L)) {
    stop("FACTOR_COUNT is outside the admissible range.")
  }

  averages <- scale(cbind(fit$y_bar, fit$x_bar), center = TRUE, scale = FALSE)
  factors <- svd(averages, nu = factor_count, nv = 0L)$u[, seq_len(factor_count),
                                                        drop = FALSE] * sqrt(Tn)
  factor_controls <- cbind(1, factors)
  factor_inverse <- matrix_pinv(crossprod(factor_controls))

  p_units <- array(NA_real_, c(factor_count + 1L, k + 2L, N))
  pi_units <- matrix(NA_real_, N, factor_count + 1L)
  Pi_units <- array(NA_real_, c(factor_count + 1L, k, N))

  for (i in seq_len(N)) {
    Xi <- X[, i, , drop = FALSE][, 1L, ]
    yi <- Y[, i]
    beta_i <- fit$betas[i, ]
    coefficient_x <- factor_inverse %*% crossprod(factor_controls, Xi)
    coefficient_y <- factor_inverse %*%
      crossprod(factor_controls, yi - Xi %*% beta_i)
    Ai <- coefficient_x[1L, ]
    Gamma_i <- coefficient_x[-1L, , drop = FALSE]
    ai <- coefficient_y[1L]
    gamma_i <- coefficient_y[-1L]
    Bi <- c(ai + sum(Ai * beta_i), Ai)
    Ci <- cbind(gamma_i + Gamma_i %*% beta_i, Gamma_i)
    p_units[, , i] <- rbind(c(1, Bi), cbind(0, Ci))
    pi_units[i, ] <- c(ai, gamma_i)
    Pi_units[, , i] <- rbind(Ai, Gamma_i)
  }

  pbar <- apply(p_units, c(1L, 2L), mean)
  pbar_plus <- matrix_pinv(pbar)
  b0_unit <- matrix(NA_real_, N, k)
  b1_unit <- matrix(NA_real_, N, k)
  for (i in seq_len(N)) {
    selector <- cbind(0, fit$betas[i, ], diag(k))
    sigma_vector <- c(0, fit$residual_se[i], rep(1, k))
    inverse_gram <- matrix_pinv(fit$gram[, , i])
    b0_unit[i, ] <- selector %*% pbar_plus %*% pi_units[i, ]
    b1_unit[i, ] <- inverse_gram %*% t(Pi_units[, , i]) %*%
      t(pbar_plus) %*% sigma_vector
  }

  b0_hat <- sqrt(Tn / N) * colMeans(b0_unit)
  b1_hat <- sqrt(Tn / N) * colMeans(b1_unit)
  beta_mg <- colMeans(fit$betas)
  beta_abc <- beta_mg - (b0_hat + b1_hat) / sqrt(N * Tn)
  unit_corrected_slopes <- fit$betas - (b0_unit + b1_unit) / N
  equivalent <- colMeans(unit_corrected_slopes)
  if (!isTRUE(all.equal(beta_abc, equivalent, tolerance = 1e-8))) {
    stop("ABC scaling identity failed.")
  }
  list(beta_abc = beta_abc, b0_hat = b0_hat, b1_hat = b1_hat,
       b0_unit = b0_unit, b1_unit = b1_unit,
       unit_corrected_slopes = unit_corrected_slopes, pbar = pbar)
}

standard_ccemg_covariance <- function(betas) {
  N <- nrow(betas)
  centered <- sweep(betas, 2L, colMeans(betas), "-")
  crossprod(centered) / (N * (N - 1L))
}

whole_sum_dwb_l1 <- function(fit, abc, bootstrap_reps, confidence) {
  N <- fit$N
  Tn <- fit$T
  if (N <= 1L) stop("The l_N=1 covariance requires at least two units.")

  centered <- sweep(abc$unit_corrected_slopes, 2L, abc$beta_abc, "-")
  unit_influence <- sqrt(Tn) * sqrt(N / (N - 1L)) * centered
  multipliers <- matrix(rnorm(N * bootstrap_reps), nrow = N,
                        ncol = bootstrap_reps)
  statistic_draws <- t(multipliers) %*% unit_influence / sqrt(N)

  conditional_covariance <- crossprod(unit_influence) / N
  estimator_covariance <- conditional_covariance / (N * Tn)
  corrected_slope_covariance <- crossprod(centered) / (N * (N - 1L))
  if (!isTRUE(all.equal(estimator_covariance, corrected_slope_covariance,
                        tolerance = 1e-10))) {
    stop("The l_N=1 covariance identity failed.")
  }
  bootstrap_covariance <- stats::cov(statistic_draws) / (N * Tn)

  alpha <- 1 - confidence
  q_low <- apply(statistic_draws, 2L, stats::quantile, probs = alpha / 2,
                 names = FALSE, type = 8L)
  q_high <- apply(statistic_draws, 2L, stats::quantile, probs = 1 - alpha / 2,
                  names = FALSE, type = 8L)
  scale_nt <- sqrt(N * Tn)
  lower <- abc$beta_abc - q_high / scale_nt
  upper <- abc$beta_abc - q_low / scale_nt
  se <- sqrt(pmax(diag(estimator_covariance), 0))
  list(lower = lower, upper = upper, se = se, q_low = q_low, q_high = q_high,
       covariance = estimator_covariance,
       conditional_covariance = conditional_covariance,
       bootstrap_covariance = bootstrap_covariance,
       unit_influence = unit_influence, bandwidth = 1L,
       finite_sample_correction = TRUE,
       covariance_identity_max_abs_error = max(abs(
         estimator_covariance - corrected_slope_covariance)),
       bootstrap_covariance_max_abs_difference = max(abs(
         bootstrap_covariance - estimator_covariance)))
}

load_australia_panel <- function(data) {
  subset <- data
  years <- sort(unique(subset$Year))
  units <- sort(unique(subset$Id))
  Tn <- length(years)
  N <- length(units)
  k <- length(variables)
  Y <- matrix(NA_real_, Tn, N)
  X <- array(NA_real_, c(Tn, N, k))
  regions <- character(N)
  source_states <- character(N)
  for (i in seq_along(units)) {
    unit_data <- subset[subset$Id == units[i], , drop = FALSE]
    unit_data <- unit_data[order(unit_data$Year), , drop = FALSE]
    if (nrow(unit_data) != Tn || !identical(as.integer(unit_data$Year), as.integer(years))) {
      stop("Australia: panel is not balanced for Id ", units[i])
    }
    Y[, i] <- unit_data$ln_yields
    X[, i, ] <- as.matrix(unit_data[, variables])
    regions[i] <- as.character(unit_data$Region[1L])
    source_states[i] <- as.character(unit_data$State_norm[1L])
  }
  if (anyNA(Y) || anyNA(X)) stop("Australia: missing model values.")
  list(Y = Y, X = X, years = years, units = units, regions = regions,
       source_states = source_states)
}

estimate_australia <- function(data, bootstrap_reps, confidence,
                               factor_count, seed) {
  set.seed(seed)
  panel <- load_australia_panel(data)
  fit <- fit_cce(panel$Y, panel$X)
  beta_mg <- colMeans(fit$betas)
  abc <- analytical_bias_correction(fit, panel$Y, panel$X, factor_count)
  standard_covariance <- standard_ccemg_covariance(fit$betas)
  standard_se <- sqrt(pmax(diag(standard_covariance), 0))
  z_value <- stats::qnorm((1 + confidence) / 2)
  standard_lower <- beta_mg - z_value * standard_se
  standard_upper <- beta_mg + z_value * standard_se
  bootstrap <- whole_sum_dwb_l1(fit, abc, bootstrap_reps, confidence)

  method <- rep(c("CCEMG", abc_method), times = length(variables))
  estimate <- as.vector(rbind(beta_mg, abc$beta_abc))
  se <- as.vector(rbind(standard_se, bootstrap$se))
  lower <- as.vector(rbind(standard_lower, bootstrap$lower))
  upper <- as.vector(rbind(standard_upper, bootstrap$upper))
  result <- data.frame(
    Sample = sample_label,
    Variable = rep(variables, each = 2L),
    Variable_label = rep(unname(variable_labels[variables]), each = 2L),
    Method = method,
    Estimate = estimate,
    Standard_error = se,
    CI_lower = lower,
    CI_upper = upper,
    CI_length = upper - lower,
    p_value_normal = 2 * stats::pnorm(-abs(estimate / se)),
    N = fit$N, T = fit$T, Confidence_level = confidence,
    Bootstrap_reps = ifelse(method == "CCEMG", NA, bootstrap_reps),
    DWB_bandwidth = ifelse(method == "CCEMG", NA, bootstrap$bandwidth),
    stringsAsFactors = FALSE
  )
  state_counts <- table(panel$source_states)
  diagnostic <- data.frame(
    Sample = sample_label, N = fit$N, T = fit$T,
    Start_year = min(panel$years), End_year = max(panel$years),
    State_counts = paste(paste0(names(state_counts), ":", as.integer(state_counts)),
                         collapse = " | "),
    Joint_rank_min = min(fit$joint_rank), Joint_rank_max = max(fit$joint_rank),
    Pbar_rank = qr(abc$pbar, tol = 1e-10)$rank,
    Pbar_condition_number = kappa(abc$pbar), ABC_factor_count = factor_count,
    DWB_bandwidth = bootstrap$bandwidth,
    DWB_dimension = "cross_section",
    DWB_finite_sample_correction = bootstrap$finite_sample_correction,
    Conditional_covariance_min_eigenvalue = min(eigen(
      bootstrap$conditional_covariance, symmetric = TRUE, only.values = TRUE)$values),
    Conditional_covariance_max_eigenvalue = max(eigen(
      bootstrap$conditional_covariance, symmetric = TRUE, only.values = TRUE)$values),
    Covariance_identity_max_abs_error = bootstrap$covariance_identity_max_abs_error,
    Bootstrap_covariance_max_abs_difference =
      bootstrap$bootstrap_covariance_max_abs_difference
  )
  slope_detail <- do.call(rbind, lapply(seq_len(fit$N), function(i) {
    data.frame(Sample = sample_label, Source_state = panel$source_states[i],
               Region = panel$regions[i], Variable = variables,
               Individual_CCE_slope = fit$betas[i, ],
               ABC_b0_unit_component = abc$b0_unit[i, ],
               ABC_b1_unit_component = abc$b1_unit[i, ],
               Individual_ABC_corrected_slope = abc$unit_corrected_slopes[i, ],
               DWB_centered_contribution = bootstrap$unit_influence[i, ])
  }))
  list(result = result, diagnostic = diagnostic, slopes = slope_detail)
}

plot_australia <- function(results, output_file) {
  results$Method_y <- ifelse(results$Method == "CCEMG", 1, 0)
  results$Variable_label <- factor(results$Variable_label,
                                   levels = unname(variable_labels[variables]))
  colors <- c("CCEMG" = "#D55E00", "ABC + whole-sum DWB" = "#0072B2")
  plot <- ggplot2::ggplot(results) +
    ggplot2::geom_rect(ggplot2::aes(xmin = CI_lower, xmax = CI_upper,
                                    ymin = Method_y - 0.18, ymax = Method_y + 0.18,
                                    fill = Method), alpha = 0.22) +
    ggplot2::geom_point(ggplot2::aes(x = Estimate, y = Method_y, color = Method), size = 2.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::facet_wrap(~Variable_label, scales = "free_x", ncol = 1) +
    ggplot2::scale_y_continuous(breaks = c(0, 1), labels = c("ABC + DWB", "CCEMG")) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = "Australia: CCEMG and ABC whole-sum DWB mean effects",
                  x = "Estimated mean effect with 95% confidence interval", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank())
  ggplot2::ggsave(output_file, plot, width = 10.5, height = 12.5, dpi = 220,
                  bg = "white")
}

data_file <- Sys.getenv("DATA_FILE", unset = file.path(script_dir, "data", "Aus_data_with_agg_new.xlsx"))
output_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path(script_dir, "results"))
bootstrap_reps <- read_integer("BOOT_REPS", 9999L, 99L)
confidence <- read_number("CONFIDENCE_LEVEL", 0.95)
factor_count <- read_integer("FACTOR_COUNT", 2L, 1L)
seed <- read_integer("EMPIRICAL_SEED", 20261001L, 1L)
if (confidence <= 0 || confidence >= 1) stop("CONFIDENCE_LEVEL must lie in (0,1).")

data <- as.data.frame(readxl::read_excel(data_file))
required_columns <- c("Id", "Region", "Year", "State", "ln_yields", variables)
if (!all(required_columns %in% names(data))) {
  stop("Missing required columns: ", paste(setdiff(required_columns, names(data)), collapse = ", "))
}
data$State_norm <- toupper(trimws(as.character(data$State)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fit <- estimate_australia(data, bootstrap_reps, confidence, factor_count, seed)
results <- fit$result
diagnostics <- fit$diagnostic
slopes <- fit$slopes
utils::write.csv(results, file.path(output_dir, "empirical_ccemg_abc_results.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(output_dir, "empirical_diagnostics.csv"), row.names = FALSE)
utils::write.csv(slopes, file.path(output_dir, "individual_cce_slopes_and_abc_components.csv"), row.names = FALSE)
lengths <- reshape(results[, c("Sample", "Variable", "Variable_label", "Method", "CI_length")],
                   idvar = c("Sample", "Variable", "Variable_label"),
                   timevar = "Method", direction = "wide")
utils::write.csv(lengths, file.path(output_dir, "confidence_interval_lengths.csv"), row.names = FALSE)
plot_australia(results, file.path(output_dir,
                                  "Australia_mean_effects_with_shaded_CI.png"))
print(results)
cat("\nResults written to:", normalizePath(output_dir), "\n")
