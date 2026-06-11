rm(list = ls())

options(stringsAsFactors = FALSE)

library(Rcpp)
library(sf)
library(tictoc)

get_project_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])

  if (length(path) > 0) {
    return(dirname(normalizePath(path[1])))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }

  normalizePath(getwd())
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

compute_alpha_max <- function(Z_mat) {
  z_positive <- Z_mat[upper.tri(Z_mat) & (Z_mat > 0)]
  if (length(z_positive) == 0) {
    return(1.0)
  }
  -log(0.5) / median(z_positive)
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
  verbose = TRUE
) {
  message("Running DAGARBayes for ", dataset_name, "...")

  output_dir <- file.path(results_dir, dataset_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gpkg_path <- file.path(data_dir, gpkg_name)
  adj_path <- file.path(data_dir, adj_csv_name)

  gdf <- st_read(gpkg_path, quiet = TRUE)
  W <- as.matrix(read.csv(adj_path, check.names = FALSE))
  storage.mode(W) <- "double"

  if (!(id_col %in% names(gdf))) {
    stop("Missing id column '", id_col, "' in ", gpkg_name)
  }

  area_ids <- as.character(gdf[[id_col]])
  covariate_raw <- as.numeric(gdf[[covariate_col]])
  covariate_std <- as.numeric(scale(covariate_raw))
  Z_mat <- as.matrix(dist(covariate_std, diag = TRUE, upper = TRUE))
  alpha_max <- compute_alpha_max(Z_mat)

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
    alpha_max = alpha_max,
    threshold = threshold,
    stringsAsFactors = FALSE
  )

  config_df <- data.frame(
    dataset_name = dataset_name,
    gpkg_path = gpkg_path,
    adjacency_path = adj_path,
    N = nrow(gdf),
    edge_count = sum(W[upper.tri(W)] == 1),
    mean_neighbors = mean(rowSums(W)),
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    n_adapt = n_adapt,
    alpha_max = alpha_max,
    threshold = threshold,
    verbose = verbose,
    stringsAsFactors = FALSE
  )

  write.csv(summary_df, file.path(output_dir, paste0(dataset_name, "_posterior_summary.csv")), row.names = FALSE)
  write.csv(posterior_samples, file.path(output_dir, paste0(dataset_name, "_posterior_samples.csv")), row.names = FALSE)
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
      areas = area_table,
      edges = edge_table,
      acceptance = acceptance_df,
      runtime = runtime_df
    )
  )
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

Rcpp::sourceCpp("dagar_poisson_boundary_mwg.cpp")

tic()
glasgow_fit <- run_dagarbayes_application(
  dataset_name = "glasgow",
  gpkg_name = "respiratory_data.gpkg",
  adj_csv_name = "adjacency_matrix.csv",
  id_col = "IZ",
  observed_col = "observed",
  expected_col = "expected",
  covariate_col = "incomedep",
  data_dir = data_dir,
  results_dir = results_dir
)
toc()

tic()
california_fit <- run_dagarbayes_application(
  dataset_name = "california",
  gpkg_name = "respiratory_data_california.gpkg",
  adj_csv_name = "adjacency_matrix_california.csv",
  id_col = "county",
  observed_col = "lung_O_count",
  expected_col = "lung_E_count",
  covariate_col = "smoking",
  data_dir = data_dir,
  results_dir = results_dir
)
toc()

message("Finished DAGARBayes runs for Glasgow and California.")
