rm(list = ls())

options(stringsAsFactors = FALSE)

library(CARBayes)
library(sf)

setwd("C:/Dati/Lavori/Aiello_Banerjee_2026/Code/ABI_poisson_regression/Real Data Analysis")

# get_project_dir <- function() {
#   args <- commandArgs(trailingOnly = FALSE)
#   file_arg <- "--file="
#   path <- sub(file_arg, "", args[grep(file_arg, args)])

#   if (length(path) > 0) {
#     return(dirname(normalizePath(path[1])))
#   }

#   if (!is.null(sys.frames()[[1]]$ofile)) {
#     return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
#   }

#   normalizePath(getwd())
# }

get_project_dir <- function() {
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
  thin = 20
) {
  message("Running CARBayes for ", dataset_name, "...")

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
  covariate_std <- as.numeric(scale(covariate_raw))
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
  write.csv(area_table, file.path(output_dir, paste0(dataset_name, "_area_table.csv")), row.names = FALSE)
  write.csv(edge_table, file.path(output_dir, paste0(dataset_name, "_edge_table.csv")), row.names = FALSE)
  write.csv(W, file.path(output_dir, paste0(dataset_name, "_adjacency_matrix.csv")), row.names = FALSE)
  write.csv(Z_mat, file.path(output_dir, paste0(dataset_name, "_dissimilarity_matrix.csv")), row.names = FALSE)
  write.csv(W_border_prob, file.path(output_dir, paste0(dataset_name, "_W_border_prob.csv")), row.names = FALSE)
  write.csv(W_posterior, file.path(output_dir, paste0(dataset_name, "_W_posterior.csv")), row.names = FALSE)

  st_write(gdf, file.path(output_dir, paste0(dataset_name, "_data.gpkg")), delete_dsn = TRUE, quiet = TRUE)
  saveRDS(chain, file.path(output_dir, paste0(dataset_name, "_chain.rds")))

  capture.output(print(chain), file = file.path(output_dir, paste0(dataset_name, "_model_summary.txt")))

  message(
    dataset_name, ": saved posterior summaries, samples, fitted risks, and edge-level boundary outputs to ",
    output_dir
  )

  invisible(
    list(
      chain = chain,
      summary = summary_df,
      areas = area_table,
      edges = edge_table
    )
  )
}

project_dir <- get_project_dir()
data_dir <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "setup_and_diagnostics")

if (!dir.exists(data_dir)) {
  stop("Could not find Data directory under: ", project_dir)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(tictoc)

tic()

glasgow_fit <- run_carbayes_application(
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

california_fit <- run_carbayes_application(
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

message("Finished CARBayes runs for Glasgow and California.")
