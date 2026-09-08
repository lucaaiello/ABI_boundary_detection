rm(list = ls())

options(stringsAsFactors = FALSE)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2026082201)

library(Rcpp)
library(sf)
library(tictoc)

if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Package 'coda' is required to compute MCMC ESS and MCSE diagnostics.")
}

find_real_data_dir <- function(candidate_dir) {
  candidate_dir <- normalizePath(candidate_dir, winslash = "/", mustWork = FALSE)
  candidates <- unique(c(
    candidate_dir,
    file.path(candidate_dir, "Real Data Analysis"),
    dirname(candidate_dir),
    file.path(dirname(candidate_dir), "Real Data Analysis")
  ))

  for (candidate in candidates) {
    if (dir.exists(file.path(candidate, "Data"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NULL
}

get_project_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])

  if (length(path) > 0) {
    script_dir <- dirname(normalizePath(path[1], winslash = "/", mustWork = TRUE))
    real_data_dir <- find_real_data_dir(script_dir)
    if (!is.null(real_data_dir)) {
      return(real_data_dir)
    }
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
    real_data_dir <- find_real_data_dir(script_dir)
    if (!is.null(real_data_dir)) {
      return(real_data_dir)
    }
  }

  real_data_dir <- find_real_data_dir(getwd())
  if (!is.null(real_data_dir)) {
    return(real_data_dir)
  }

  stop("Could not locate the 'Real Data Analysis' directory from working directory: ", getwd())
}

posterior_summary_row <- function(samples, parameter) {
  q <- quantile(samples, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.frame(
    parameter = parameter,
    median = unname(q[2]),
    lower_95 = unname(q[1]),
    upper_95 = unname(q[3]),
    stringsAsFactors = FALSE
  )
}

running_mean <- function(x) {
  cumsum(x) / seq_along(x)
}

estimate_chain_precision <- function(samples, parameter) {
  samples <- as.numeric(samples)
  samples <- samples[is.finite(samples)]
  n_draws <- length(samples)

  if (n_draws < 3L || !is.finite(var(samples)) || var(samples) <= 0) {
    warning("Cannot estimate ESS and MCSE for ", parameter, ": insufficient variation in the retained draws.")
    return(list(
      n_draws = n_draws,
      spectral_variance = NA_real_,
      effective_sample_size = NA_real_,
      method = NA_character_
    ))
  }

  spectral_fit <- tryCatch(
    coda::spectrum0.ar(coda::mcmc(samples)),
    error = function(error) NULL
  )
  method <- "coda::spectrum0.ar"

  if (is.null(spectral_fit) || !is.finite(spectral_fit$spec[[1]]) || spectral_fit$spec[[1]] <= 0) {
    spectral_fit <- tryCatch(
      coda::spectrum0(coda::mcmc(samples)),
      error = function(error) NULL
    )
    method <- "coda::spectrum0"
  }

  if (is.null(spectral_fit) || !is.finite(spectral_fit$spec[[1]]) || spectral_fit$spec[[1]] <= 0) {
    warning("Cannot estimate ESS and MCSE for ", parameter, ": spectral variance estimation failed.")
    return(list(
      n_draws = n_draws,
      spectral_variance = NA_real_,
      effective_sample_size = NA_real_,
      method = NA_character_
    ))
  }

  spectral_variance <- as.numeric(spectral_fit$spec[[1]])
  effective_sample_size <- n_draws * var(samples) / spectral_variance

  list(
    n_draws = n_draws,
    spectral_variance = spectral_variance,
    effective_sample_size = effective_sample_size,
    method = method
  )
}

summarize_mcmc_precision <- function(chain_list, parameter) {
  clean_chains <- lapply(chain_list, function(samples) {
    samples <- as.numeric(samples)
    samples[is.finite(samples)]
  })
  chain_stats <- lapply(clean_chains, estimate_chain_precision, parameter = parameter)
  n_draws_total <- sum(vapply(chain_stats, function(x) x$n_draws, numeric(1)))
  pooled_samples <- unlist(clean_chains, use.names = FALSE)
  posterior_sd <- if (length(pooled_samples) > 1L) sd(pooled_samples) else NA_real_

  valid <- vapply(
    chain_stats,
    function(x) is.finite(x$spectral_variance) && is.finite(x$effective_sample_size),
    logical(1)
  )

  if (n_draws_total == 0L || !all(valid)) {
    effective_sample_size <- NA_real_
    spectral_variance <- NA_real_
    mcse_mean <- NA_real_
  } else {
    effective_sample_size <- sum(vapply(chain_stats, function(x) x$effective_sample_size, numeric(1)))
    weighted_spectral_sum <- sum(vapply(
      chain_stats,
      function(x) x$n_draws * x$spectral_variance,
      numeric(1)
    ))
    spectral_variance <- weighted_spectral_sum / n_draws_total
    mcse_mean <- sqrt(weighted_spectral_sum) / n_draws_total
  }

  data.frame(
    parameter = parameter,
    n_chains = length(chain_list),
    n_draws_total = n_draws_total,
    posterior_sd = posterior_sd,
    spectral_variance = spectral_variance,
    effective_sample_size = effective_sample_size,
    ess_per_draw = effective_sample_size / n_draws_total,
    mcse_mean = mcse_mean,
    mcse_over_posterior_sd = mcse_mean / posterior_sd,
    method = paste(unique(na.omit(vapply(chain_stats, function(x) x$method, character(1)))), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

standardize_population <- function(x, variable_name = "boundary covariate") {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(variable_name, " must contain only finite values.")
  }

  centered <- x - mean(x)
  sd_population <- sqrt(mean(centered^2))
  if (!is.finite(sd_population) || sd_population <= 0) {
    stop(variable_name, " must be non-constant.")
  }

  centered / sd_population
}

save_mcmc_diagnostics <- function(posterior_samples, output_dir, dataset_name) {
  parameter_names <- setdiff(names(posterior_samples), "draw")

  if (length(parameter_names) == 0) {
    return(invisible(NULL))
  }

  pdf(
    file.path(output_dir, paste0(dataset_name, "_mcmc_diagnostics.pdf")),
    width = 10,
    height = 9
  )
  on.exit(dev.off(), add = TRUE)

  for (parameter_name in parameter_names) {
    samples <- as.numeric(posterior_samples[[parameter_name]])
    draw_index <- posterior_samples$draw
    sample_mean <- mean(samples, na.rm = TRUE)

    par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

    plot(
      draw_index,
      samples,
      type = "l",
      xlab = "Saved draw",
      ylab = parameter_name,
      main = paste(dataset_name, "-", parameter_name, "traceplot")
    )
    abline(h = sample_mean, col = "red", lty = 2)

    acf(
      samples,
      lag.max = min(200, length(samples) - 1),
      main = paste(dataset_name, "-", parameter_name, "autocorrelation")
    )

    plot(
      draw_index,
      running_mean(samples),
      type = "l",
      xlab = "Saved draw",
      ylab = paste("Running mean of", parameter_name),
      main = paste(dataset_name, "-", parameter_name, "running mean")
    )
    abline(h = sample_mean, col = "red", lty = 2)
  }

  invisible(NULL)
}

build_edge_table <- function(W, Z_mat, W_border_prob, W_posterior, area_ids) {
  edge_idx <- which(upper.tri(W) & (W == 1), arr.ind = TRUE)

  data.frame(
    from_index = edge_idx[, 1],
    to_index = edge_idx[, 2],
    from_id = area_ids[edge_idx[, 1]],
    to_id = area_ids[edge_idx[, 2]],
    edge_dissimilarity = Z_mat[edge_idx],
    boundary_prob = W_border_prob[edge_idx],
    w_posterior = W_posterior[edge_idx],
    boundary_median = as.integer(W_posterior[edge_idx] < 0.5),
    stringsAsFactors = FALSE
  )
}

compute_eta_scaling <- function(W, Z_mat, threshold = log(2.0)) {
  edge_mask <- upper.tri(W) & (W == 1)
  z_edges <- Z_mat[edge_mask]

  if (length(z_edges) == 0) {
    stop("Cannot compute eta scaling because the adjacency matrix has no edges.")
  }

  z_median_adjacent <- median(z_edges)
  eta_max <- threshold / (z_median_adjacent + 1e-8)

  list(
    z_median_adjacent = z_median_adjacent,
    eta_max = eta_max
  )
}

validate_application_inputs <- function(
  gdf,
  W,
  id_col,
  observed_col,
  expected_col,
  covariate_col,
  gpkg_name,
  adj_csv_name
) {
  required_columns <- c(id_col, observed_col, expected_col, covariate_col)
  missing_columns <- setdiff(required_columns, names(gdf))
  if (length(missing_columns) > 0) {
    stop(
      "Missing column(s) in ", gpkg_name, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  n_areas <- nrow(gdf)
  if (!identical(dim(W), c(n_areas, n_areas))) {
    stop(
      "Adjacency matrix in ", adj_csv_name, " has dimensions ",
      paste(dim(W), collapse = " x "), " but the spatial data contain ",
      n_areas, " areas."
    )
  }
  if (any(!is.finite(W))) {
    stop("Adjacency matrix contains non-finite values: ", adj_csv_name)
  }
  if (any(!(W %in% c(0, 1)))) {
    stop("Adjacency matrix must be binary: ", adj_csv_name)
  }
  if (!all(W == t(W))) {
    stop("Adjacency matrix must be symmetric: ", adj_csv_name)
  }
  if (any(diag(W) != 0)) {
    stop("Adjacency matrix diagonal must be zero: ", adj_csv_name)
  }

  observed <- as.numeric(gdf[[observed_col]])
  expected <- as.numeric(gdf[[expected_col]])
  covariate <- as.numeric(gdf[[covariate_col]])
  if (any(!is.finite(observed)) || any(observed < 0)) {
    stop("Observed counts must be finite and non-negative in ", gpkg_name)
  }
  if (any(!is.finite(expected)) || any(expected <= 0)) {
    stop("Expected counts must be finite and strictly positive in ", gpkg_name)
  }
  if (any(!is.finite(covariate)) || sd(covariate) <= 0) {
    stop("Boundary covariate must be finite and non-constant in ", gpkg_name)
  }

  invisible(NULL)
}

extract_acceptance_summary <- function(accept_vec) {
  accept_vec <- as.numeric(accept_vec)

  safe_rate <- function(num_idx, den_idx) {
    if (length(accept_vec) < den_idx || is.na(accept_vec[den_idx]) || accept_vec[den_idx] <= 0) {
      return(NA_real_)
    }
    accept_vec[num_idx] / accept_vec[den_idx]
  }

  data.frame(
    parameter = c("beta_0", "phi", "eta", "rho"),
    accepted = c(
      if (length(accept_vec) >= 1) accept_vec[1] else NA_real_,
      if (length(accept_vec) >= 3) accept_vec[3] else NA_real_,
      if (length(accept_vec) >= 5) accept_vec[5] else NA_real_,
      if (length(accept_vec) >= 7) accept_vec[7] else NA_real_
    ),
    attempted = c(
      if (length(accept_vec) >= 2) accept_vec[2] else NA_real_,
      if (length(accept_vec) >= 4) accept_vec[4] else NA_real_,
      if (length(accept_vec) >= 6) accept_vec[6] else NA_real_,
      if (length(accept_vec) >= 8) accept_vec[8] else NA_real_
    ),
    acceptance_rate = c(
      safe_rate(1, 2),
      safe_rate(3, 4),
      safe_rate(5, 6),
      safe_rate(7, 8)
    ),
    stringsAsFactors = FALSE
  )
}

run_dagarbayes_application <- function(
  dataset_name,
  gpkg_name,
  adj_csv_name,
  id_col,
  observed_col,
  expected_col,
  covariate_col,
  data_dir,
  results_dir,
  n_iter = 300000L,
  burnin = 100000L,
  thin = 20L,
  n_adapt = 100000L,
  threshold = log(2.0),
  verbose = TRUE,
  seed = NULL
) {
  message("Running DAGARBayes for ", dataset_name, "...")

  if (!is.null(seed)) {
    set.seed(seed)
    message(dataset_name, ": using R seed ", seed)
  }

  output_dir <- file.path(results_dir, dataset_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gpkg_path <- file.path(data_dir, gpkg_name)
  adj_path <- file.path(data_dir, adj_csv_name)

  if (!file.exists(gpkg_path)) {
    stop("Could not find spatial data file: ", gpkg_path)
  }
  if (!file.exists(adj_path)) {
    stop("Could not find adjacency matrix: ", adj_path)
  }

  gdf <- st_read(gpkg_path, quiet = TRUE)
  W <- as.matrix(read.csv(adj_path, check.names = FALSE))
  storage.mode(W) <- "double"

  validate_application_inputs(
    gdf = gdf,
    W = W,
    id_col = id_col,
    observed_col = observed_col,
    expected_col = expected_col,
    covariate_col = covariate_col,
    gpkg_name = gpkg_name,
    adj_csv_name = adj_csv_name
  )

  area_ids <- as.character(gdf[[id_col]])
  covariate_raw <- as.numeric(gdf[[covariate_col]])
  covariate_std <- standardize_population(covariate_raw)
  Z_mat <- as.matrix(dist(covariate_std, diag = TRUE, upper = TRUE))
  eta_scaling <- compute_eta_scaling(W, Z_mat, threshold)
  alpha_max <- eta_scaling$eta_max

  message(
    dataset_name, ": N = ", nrow(gdf),
    ", edges = ", sum(W[upper.tri(W)] == 1),
    ", median adjacent dissimilarity = ",
    format(eta_scaling$z_median_adjacent, digits = 6),
    ", eta upper bound M = ", format(alpha_max, digits = 6)
  )

  gdf$covariate_std <- covariate_std

  fit_start <- proc.time()[["elapsed"]]
  fit <- dagar_poisson_boundary_mwg2(
    y_r = as.numeric(gdf[[observed_col]]),
    e_r = as.numeric(gdf[[expected_col]]),
    A_r = W,
    Z_r = Z_mat,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    n_adapt = n_adapt,
    alpha_max = alpha_max,
    threshold = threshold,
    save_phi = FALSE,
    save_fitted = TRUE,
    save_loglike = FALSE,
    verbose = verbose
  )
  elapsed_sec <- proc.time()[["elapsed"]] - fit_start

  beta0_samples <- as.numeric(fit$beta0)
  sigma2_samples <- as.numeric(fit$sigma2_w)
  eta_samples <- as.numeric(fit$eta)
  rho_samples <- as.numeric(fit$rho)

  summary_df <- do.call(
    rbind,
    list(
      posterior_summary_row(beta0_samples, "beta_0"),
      posterior_summary_row(sigma2_samples, "sigma2_w"),
      posterior_summary_row(eta_samples, "eta"),
      posterior_summary_row(rho_samples, "rho")
    )
  )
  rownames(summary_df) <- NULL

  posterior_samples <- data.frame(
    draw = seq_along(beta0_samples),
    beta_0 = beta0_samples,
    sigma2_w = sigma2_samples,
    eta = eta_samples,
    rho = rho_samples,
    stringsAsFactors = FALSE
  )

  precision_df <- do.call(
    rbind,
    lapply(
      setdiff(names(posterior_samples), "draw"),
      function(parameter) summarize_mcmc_precision(list(posterior_samples[[parameter]]), parameter)
    )
  )
  precision_df <- cbind(dataset_name = dataset_name, precision_df, stringsAsFactors = FALSE)

  save_mcmc_diagnostics(posterior_samples, output_dir, dataset_name)

  area_table <- data.frame(
    index = seq_len(nrow(gdf)),
    area_id = area_ids,
    observed = as.numeric(gdf[[observed_col]]),
    expected = as.numeric(gdf[[expected_col]]),
    covariate_raw = covariate_raw,
    covariate_std = covariate_std,
    fitted_values = as.numeric(fit$fitted_values),
    risk_dagarbayes = as.numeric(fit$fitted_values) / as.numeric(gdf[[expected_col]]),
    stringsAsFactors = FALSE
  )

  W_border_prob <- as.matrix(fit$localised_structure$W.border.prob)
  W_posterior <- as.matrix(fit$localised_structure$W.posterior)
  edge_table <- build_edge_table(W, Z_mat, W_border_prob, W_posterior, area_ids)

  acceptance_df <- extract_acceptance_summary(fit$accept)
  runtime_df <- data.frame(
    dataset_name = dataset_name,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    n_adapt = n_adapt,
    n_saved_draws = nrow(posterior_samples),
    elapsed_sec = elapsed_sec,
    seconds_per_saved_draw = elapsed_sec / max(1L, nrow(posterior_samples)),
    seconds_per_1000_saved_draws = 1000.0 * elapsed_sec / max(1L, nrow(posterior_samples)),
    seed = seed,
    z_median_adjacent = eta_scaling$z_median_adjacent,
    M_real = alpha_max,
    alpha_max = alpha_max,
    threshold = threshold,
    stringsAsFactors = FALSE
  )

  config_df <- data.frame(
    dataset_name = dataset_name,
    gpkg_path = file.path("Data", gpkg_name),
    adjacency_path = file.path("Data", adj_csv_name),
    N = nrow(gdf),
    edge_count = sum(W[upper.tri(W)] == 1),
    mean_neighbors = mean(rowSums(W)),
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    n_adapt = n_adapt,
    seed = seed,
    z_median_adjacent = eta_scaling$z_median_adjacent,
    M_real = alpha_max,
    alpha_max = alpha_max,
    threshold = threshold,
    verbose = verbose,
    stringsAsFactors = FALSE
  )

  write.csv(summary_df, file.path(output_dir, paste0(dataset_name, "_posterior_summary.csv")), row.names = FALSE)
  write.csv(posterior_samples, file.path(output_dir, paste0(dataset_name, "_posterior_samples.csv")), row.names = FALSE)
  write.csv(precision_df, file.path(output_dir, paste0(dataset_name, "_mcmc_precision.csv")), row.names = FALSE)
  write.csv(area_table, file.path(output_dir, paste0(dataset_name, "_area_table.csv")), row.names = FALSE)
  write.csv(edge_table, file.path(output_dir, paste0(dataset_name, "_edge_table.csv")), row.names = FALSE)
  write.csv(W, file.path(output_dir, paste0(dataset_name, "_adjacency_matrix.csv")), row.names = FALSE)
  write.csv(Z_mat, file.path(output_dir, paste0(dataset_name, "_dissimilarity_matrix.csv")), row.names = FALSE)
  write.csv(W_border_prob, file.path(output_dir, paste0(dataset_name, "_W_border_prob.csv")), row.names = FALSE)
  write.csv(W_posterior, file.path(output_dir, paste0(dataset_name, "_W_posterior.csv")), row.names = FALSE)
  write.csv(as.matrix(fit$A_filtered_last), file.path(output_dir, paste0(dataset_name, "_A_filtered_last.csv")), row.names = FALSE)
  write.csv(acceptance_df, file.path(output_dir, paste0(dataset_name, "_acceptance_summary.csv")), row.names = FALSE)
  write.csv(runtime_df, file.path(output_dir, paste0(dataset_name, "_runtime_summary.csv")), row.names = FALSE)
  write.csv(config_df, file.path(output_dir, paste0(dataset_name, "_config.csv")), row.names = FALSE)

  st_write(gdf, file.path(output_dir, paste0(dataset_name, "_data.gpkg")), delete_dsn = TRUE, quiet = TRUE)
  saveRDS(fit, file.path(output_dir, paste0(dataset_name, "_fit.rds")))

  capture.output(
    {
      cat("Posterior summary\n")
      print(summary_df)
      cat("\nAcceptance summary\n")
      print(acceptance_df)
      cat("\nMCMC precision diagnostics\n")
      print(precision_df)
      cat("\nRuntime summary\n")
      print(runtime_df)
      cat("\nFit structure\n")
      str(fit, max.level = 1)
    },
    file = file.path(output_dir, paste0(dataset_name, "_model_summary.txt"))
  )

  message(
    dataset_name, ": saved posterior summaries, samples, fitted risks, edge-level boundary outputs, and diagnostics to ",
    output_dir
  )

  invisible(
    list(
      fit = fit,
      summary = summary_df,
      posterior_samples = posterior_samples,
      mcmc_precision = precision_df,
      areas = area_table,
      edges = edge_table,
      acceptance = acceptance_df,
      runtime = runtime_df
    )
  )
}

env_flag <- function(name, default = TRUE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

project_dir <- get_project_dir()
data_dir <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "setup_and_diagnostics", "DAGARBayes")
sampler_path <- file.path(project_dir, "dagar_poisson_boundary_mwg.cpp")

if (!dir.exists(data_dir)) {
  stop("Could not find Data directory under: ", project_dir)
}

if (!file.exists(sampler_path)) {
  stop("Could not find sampler source at: ", sampler_path)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

Rcpp::sourceCpp(sampler_path)

run_glasgow_california <- env_flag("RUN_GLASGOW_CALIFORNIA", TRUE)
run_south_korea <- env_flag("RUN_SOUTH_KOREA", TRUE)

if (run_glasgow_california) {
  tic()
  glasgow_fit <- run_dagarbayes_application(
    dataset_name = "glasgow",
    gpkg_name = "respiratory_data_glasgow.gpkg",
    adj_csv_name = "adjacency_matrix_glasgow.csv",
    id_col = "IZ",
    observed_col = "observed",
    expected_col = "expected",
    covariate_col = "incomedep",
    data_dir = data_dir,
    results_dir = results_dir,
    seed = 2026082201L
  )
  toc()

  california_fit <- run_dagarbayes_application(
    dataset_name = "california",
    gpkg_name = "respiratory_data_california.gpkg",
    adj_csv_name = "adjacency_matrix_california.csv",
    id_col = "county",
    observed_col = "lung_O_count",
    expected_col = "lung_E_count",
    covariate_col = "smoking",
    data_dir = data_dir,
    results_dir = results_dir,
    seed = 2026082203L
  )
  toc()
} else {
  message("Skipping Glasgow and California DAGARBayes runs because RUN_GLASGOW_CALIFORNIA is false.")
}

if (run_south_korea) {
  tic()
  south_korea_fit <- run_dagarbayes_application(
    dataset_name = "south_korea",
    gpkg_name = file.path("South_Korea", "mortality_data_south_korea.gpkg"),
    adj_csv_name = file.path("South_Korea", "adjacency_matrix_south_korea.csv"),
    id_col = "area_id",
    observed_col = "observed_lung_cancer",
    expected_col = "expected_lung_cancer",
    covariate_col = "smoking_pct",
    data_dir = data_dir,
    results_dir = results_dir,
    seed = 2026082204L
  )
  toc()
} else {
  message("Skipping South Korea DAGARBayes run because RUN_SOUTH_KOREA is false.")
}

message("Finished requested DAGARBayes runs.")
