rm(list = ls())

options(stringsAsFactors = FALSE)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2026081201)

library(CARBayes)
library(sf)

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

extract_beta0_samples <- function(chain) {
  beta_samples <- chain$samples$beta

  if (is.matrix(beta_samples) || is.data.frame(beta_samples)) {
    return(as.numeric(beta_samples[, 1]))
  }

  as.numeric(beta_samples)
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

run_carbayes_application <- function(
  dataset_name,
  gpkg_name,
  adj_csv_name,
  id_col,
  observed_col,
  expected_col,
  covariate_col,
  data_dir,
  results_dir,
  burnin = 100000,
  n.sample = 300000,
  thin = 20,
  seed = NULL
) {
  message("Running CARBayes for ", dataset_name, "...")

  if (!is.null(seed)) {
    set.seed(seed)
    message(dataset_name, ": using R seed ", seed)
  }

  output_dir <- file.path(results_dir, dataset_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gpkg_path <- file.path(data_dir, gpkg_name)
  adj_path <- file.path(data_dir, adj_csv_name)

  gdf <- st_read(gpkg_path, quiet = TRUE)
  W <- as.matrix(read.csv(adj_path, check.names = FALSE))

  if (!(id_col %in% names(gdf))) {
    stop("Missing id column '", id_col, "' in ", gpkg_name)
  }

  area_ids <- as.character(gdf[[id_col]])
  covariate_raw <- as.numeric(gdf[[covariate_col]])
  covariate_std <- standardize_population(covariate_raw)
  Z_mat <- as.matrix(dist(covariate_std, diag = TRUE, upper = TRUE))

  gdf$covariate_std <- covariate_std

  formula_obj <- as.formula(sprintf("%s ~ offset(log(%s))", observed_col, expected_col))
  Z_list <- list(Z_covariate = Z_mat)
  names(Z_list) <- paste0("Z_", covariate_col)

  chain <- S.CARdissimilarity(
    formula = formula_obj,
    data = gdf,
    family = "poisson",
    W = W,
    Z = Z_list,
    W.binary = TRUE,
    burnin = burnin,
    n.sample = n.sample,
    thin = thin
  )

  beta0_samples <- extract_beta0_samples(chain)
  alpha_samples <- as.numeric(chain$samples$alpha)

  summary_rows <- list(
    posterior_summary_row(beta0_samples, "beta_0"),
    posterior_summary_row(alpha_samples, "alpha")
  )

  if ("tau2" %in% names(chain$samples)) {
    summary_rows[[length(summary_rows) + 1]] <- posterior_summary_row(
      as.numeric(chain$samples$tau2),
      "tau2"
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  posterior_samples <- data.frame(
    draw = seq_along(alpha_samples),
    beta_0 = beta0_samples,
    alpha = alpha_samples,
    stringsAsFactors = FALSE
  )

  if ("tau2" %in% names(chain$samples)) {
    posterior_samples$tau2 <- as.numeric(chain$samples$tau2)
  }

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
    fitted_values = as.numeric(chain$fitted.values),
    risk_carbayes = as.numeric(chain$fitted.values) / as.numeric(gdf[[expected_col]]),
    stringsAsFactors = FALSE
  )

  W_border_prob <- as.matrix(chain$localised.structure$W.border.prob)
  W_posterior <- as.matrix(chain$localised.structure$W.posterior)
  edge_table <- build_edge_table(W, Z_mat, W_border_prob, W_posterior, area_ids)

  write.csv(summary_df, file.path(output_dir, paste0(dataset_name, "_posterior_summary.csv")), row.names = FALSE)
  write.csv(posterior_samples, file.path(output_dir, paste0(dataset_name, "_posterior_samples.csv")), row.names = FALSE)
  write.csv(precision_df, file.path(output_dir, paste0(dataset_name, "_mcmc_precision.csv")), row.names = FALSE)
  write.csv(area_table, file.path(output_dir, paste0(dataset_name, "_area_table.csv")), row.names = FALSE)
  write.csv(edge_table, file.path(output_dir, paste0(dataset_name, "_edge_table.csv")), row.names = FALSE)
  write.csv(W, file.path(output_dir, paste0(dataset_name, "_adjacency_matrix.csv")), row.names = FALSE)
  write.csv(Z_mat, file.path(output_dir, paste0(dataset_name, "_dissimilarity_matrix.csv")), row.names = FALSE)
  write.csv(W_border_prob, file.path(output_dir, paste0(dataset_name, "_W_border_prob.csv")), row.names = FALSE)
  write.csv(W_posterior, file.path(output_dir, paste0(dataset_name, "_W_posterior.csv")), row.names = FALSE)

  st_write(gdf, file.path(output_dir, paste0(dataset_name, "_data.gpkg")), delete_dsn = TRUE, quiet = TRUE)
  saveRDS(chain, file.path(output_dir, paste0(dataset_name, "_chain.rds")))

  capture.output(
    {
      print(chain)
      cat("\nMCMC precision diagnostics\n")
      print(precision_df)
    },
    file = file.path(output_dir, paste0(dataset_name, "_model_summary.txt"))
  )

  message(
    dataset_name, ": saved posterior summaries, samples, fitted risks, and edge-level boundary outputs to ",
    output_dir
  )

  invisible(
    list(
      chain = chain,
      summary = summary_df,
      mcmc_precision = precision_df,
      areas = area_table,
      edges = edge_table
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
results_dir <- file.path(project_dir, "setup_and_diagnostics")

if (!dir.exists(data_dir)) {
  stop("Could not find Data directory under: ", project_dir)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(tictoc)

run_glasgow_california <- env_flag("RUN_GLASGOW_CALIFORNIA", TRUE)
run_south_korea <- env_flag("RUN_SOUTH_KOREA", TRUE)

if (run_glasgow_california) {
  tic()

  glasgow_fit <- run_carbayes_application(
    dataset_name = "glasgow",
    gpkg_name = "respiratory_data_glasgow.gpkg",
    adj_csv_name = "adjacency_matrix_glasgow.csv",
    id_col = "IZ",
    observed_col = "observed",
    expected_col = "expected",
    covariate_col = "incomedep",
    data_dir = data_dir,
    results_dir = results_dir,
    seed = 2026081201
  )

  toc()

  california_fit <- run_carbayes_application(
    dataset_name = "california",
    gpkg_name = "respiratory_data_california.gpkg",
    adj_csv_name = "adjacency_matrix_california.csv",
    id_col = "county",
    observed_col = "lung_O_count",
    expected_col = "lung_E_count",
    covariate_col = "smoking",
    data_dir = data_dir,
    results_dir = results_dir,
    seed = 2026081203
  )

  toc()
} else {
  message("Skipping Glasgow and California CARBayes runs because RUN_GLASGOW_CALIFORNIA is false.")
}

if (run_south_korea) {
  tic()

  south_korea_fit <- run_carbayes_application(
    dataset_name = "south_korea",
    gpkg_name = file.path("South_Korea", "mortality_data_south_korea.gpkg"),
    adj_csv_name = file.path("South_Korea", "adjacency_matrix_south_korea.csv"),
    id_col = "area_id",
    observed_col = "observed_lung_cancer",
    expected_col = "expected_lung_cancer",
    covariate_col = "smoking_pct",
    data_dir = data_dir,
    results_dir = results_dir,
    seed = 2026081204
  )

  toc()
} else {
  message("Skipping South Korea CARBayes run because RUN_SOUTH_KOREA is false.")
}

message("Finished requested CARBayes runs.")
