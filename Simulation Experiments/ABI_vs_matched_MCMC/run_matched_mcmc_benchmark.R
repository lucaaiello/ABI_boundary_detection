script_config <- list(
  benchmark_dir = NULL,
  results_dir = NULL,
  dataset_ids = NULL,
  max_datasets = NULL,
  n_iter = 20000L,
  burnin = 10000L,
  thin = 1L,
  n_adapt = 10000L,
  seed = 123L,
  threshold = log(2.0),
  save_draws = TRUE,
  verbose_sampler = FALSE,
  use_null_makevars = FALSE
)


parse_args <- function(args) {
  output <- list()
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (!startsWith(argument, "--")) {
      stop("Unexpected argument: ", argument)
    }
    key <- substring(argument, 3L)
    if (index == length(args) || startsWith(args[[index + 1L]], "--")) {
      output[[key]] <- TRUE
      index <- index + 1L
    } else {
      output[[key]] <- args[[index + 1L]]
      index <- index + 2L
    }
  }
  output
}


pick_value <- function(arguments, config, key) {
  if (!is.null(arguments[[key]])) {
    return(arguments[[key]])
  }
  config[[gsub("-", "_", key, fixed = TRUE)]]
}


as_flag <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(default)
  }
  if (is.logical(value)) {
    return(isTRUE(value))
  }
  normalized <- tolower(trimws(as.character(value)))
  if (normalized %in% c("1", "true", "yes", "y")) {
    return(TRUE)
  }
  if (normalized %in% c("0", "false", "no", "n")) {
    return(FALSE)
  }
  stop("Cannot parse logical value: ", value)
}


as_int <- function(value, default = NULL) {
  if (is.null(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    stop("Cannot parse integer value: ", value)
  }
  parsed
}


as_num <- function(value, default = NULL) {
  if (is.null(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) {
    stop("Cannot parse numeric value: ", value)
  }
  parsed
}


get_script_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_args, value = TRUE)
  if (length(file_argument) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]))))
  }
  normalizePath(getwd())
}


find_repository_root <- function(start_dir, max_up = 8L) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (step in seq_len(max_up + 1L)) {
    checkpoint <- file.path(current, "Training", "Checkpoints", "poisson_dagar.keras")
    experiment <- file.path(
      current,
      "Simulation Experiments",
      "ABI_vs_matched_MCMC",
      "dagar_poisson_boundary_matched_mwg.cpp"
    )
    if (file.exists(checkpoint) && file.exists(experiment)) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop("Could not locate the ABI_poisson_regression repository root.")
}


portable_project_path <- function(path, project_root) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  prefix <- paste0(normalized_root, "/")
  if (startsWith(normalized_path, prefix)) {
    return(substring(normalized_path, nchar(prefix) + 1L))
  }
  normalized_path
}


parse_dataset_ids <- function(value) {
  if (is.null(value) || identical(trimws(as.character(value)), "")) {
    return(NULL)
  }
  ids <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  ids[nzchar(ids)]
}


summarize_draws <- function(samples, truth, dataset_id, parameter) {
  quantiles <- as.numeric(quantile(samples, c(0.025, 0.5, 0.975), names = FALSE))
  posterior_mean <- mean(samples)
  data.frame(
    dataset_id = dataset_id,
    parameter = parameter,
    posterior_mean = posterior_mean,
    posterior_sd = sd(samples),
    posterior_median = quantiles[[2L]],
    lower_95 = quantiles[[1L]],
    upper_95 = quantiles[[3L]],
    truth = truth,
    bias_mean = posterior_mean - truth,
    abs_error_mean = abs(posterior_mean - truth),
    covered_95 = as.integer(quantiles[[1L]] <= truth && truth <= quantiles[[3L]]),
    stringsAsFactors = FALSE
  )
}


estimate_chain_precision <- function(samples, parameter) {
  samples <- as.numeric(samples)
  samples <- samples[is.finite(samples)]
  n_draws <- length(samples)
  posterior_sd <- if (n_draws > 1L) sd(samples) else NA_real_
  if (n_draws < 3L || !is.finite(posterior_sd) || posterior_sd <= 0) {
    return(data.frame(
      parameter = parameter,
      n_draws = n_draws,
      posterior_sd = posterior_sd,
      effective_sample_size = NA_real_,
      ess_per_draw = NA_real_,
      mcse_mean = NA_real_,
      mcse_over_posterior_sd = NA_real_,
      method = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  spectral_fit <- tryCatch(coda::spectrum0.ar(coda::mcmc(samples)), error = function(error) NULL)
  if (is.null(spectral_fit) || !is.finite(spectral_fit$spec) || spectral_fit$spec <= 0) {
    spectral_variance <- stats::var(samples)
    method <- "iid_fallback"
  } else {
    spectral_variance <- as.numeric(spectral_fit$spec)
    method <- "coda::spectrum0.ar"
  }
  effective_sample_size <- n_draws * stats::var(samples) / spectral_variance
  mcse_mean <- sqrt(spectral_variance / n_draws)
  data.frame(
    parameter = parameter,
    n_draws = n_draws,
    posterior_sd = posterior_sd,
    effective_sample_size = effective_sample_size,
    ess_per_draw = effective_sample_size / n_draws,
    mcse_mean = mcse_mean,
    mcse_over_posterior_sd = mcse_mean / posterior_sd,
    method = method,
    stringsAsFactors = FALSE
  )
}


compute_auc <- function(probability, truth) {
  positive <- which(truth == 1)
  negative <- which(truth == 0)
  if (length(positive) == 0L || length(negative) == 0L) {
    return(NA_real_)
  }
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[positive]) - length(positive) * (length(positive) + 1) / 2) /
    (length(positive) * length(negative))
}


compute_average_precision <- function(probability, truth) {
  positive_total <- sum(truth == 1)
  if (positive_total == 0L) {
    return(NA_real_)
  }
  order_index <- order(probability, decreasing = TRUE)
  ordered_truth <- truth[order_index]
  true_positive <- cumsum(ordered_truth == 1)
  precision <- true_positive / seq_along(true_positive)
  sum(precision[ordered_truth == 1]) / positive_total
}


compute_boundary_probabilities <- function(eta_draws, edge_z, threshold) {
  draw_matrix <- outer(
    eta_draws,
    edge_z,
    FUN = function(eta, dissimilarity) as.numeric(dissimilarity * eta > threshold)
  )
  colMeans(draw_matrix)
}


write_gzip_csv <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(data, connection, row.names = FALSE)
}


arguments <- parse_args(commandArgs(trailingOnly = TRUE))
script_dir <- get_script_dir()
project_root <- find_repository_root(script_dir)

benchmark_dir <- pick_value(arguments, script_config, "benchmark-dir")
if (is.null(benchmark_dir)) {
  benchmark_dir <- file.path(script_dir, "datasets", "benchmark_bank_seed123_n100")
}
benchmark_dir <- normalizePath(benchmark_dir, mustWork = TRUE)

results_dir <- pick_value(arguments, script_config, "results-dir")
if (is.null(results_dir)) {
  results_dir <- file.path(benchmark_dir, "matched_mcmc_results_all100")
}

dataset_ids <- parse_dataset_ids(pick_value(arguments, script_config, "dataset-ids"))
max_datasets <- as_int(pick_value(arguments, script_config, "max-datasets"))
n_iter <- as_int(pick_value(arguments, script_config, "n-iter"), 20000L)
burnin <- as_int(pick_value(arguments, script_config, "burnin"), 10000L)
thin <- as_int(pick_value(arguments, script_config, "thin"), 1L)
n_adapt <- as_int(pick_value(arguments, script_config, "n-adapt"), 10000L)
seed <- as_int(pick_value(arguments, script_config, "seed"), 123L)
threshold <- as_num(pick_value(arguments, script_config, "threshold"), log(2.0))
save_draws <- as_flag(pick_value(arguments, script_config, "save-draws"), TRUE)
verbose_sampler <- as_flag(pick_value(arguments, script_config, "verbose-sampler"), FALSE)

if (n_iter <= burnin || burnin < 0L || thin <= 0L || n_adapt < 0L) {
  stop("Require n_iter > burnin >= 0, thin > 0, and n_adapt >= 0.")
}
if (!is.null(max_datasets) && max_datasets <= 0L) {
  stop("max_datasets must be positive when provided.")
}

if (as_flag(pick_value(arguments, script_config, "use-null-makevars"), FALSE)) {
  Sys.setenv(R_MAKEVARS_USER = "NUL")
}
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required.")
}
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Package 'coda' is required for ESS and MCSE diagnostics.")
}

manifest <- read.csv(file.path(benchmark_dir, "benchmark_manifest.csv"), check.names = FALSE)
r_inputs_root <- file.path(benchmark_dir, "r_inputs")
if (!dir.exists(r_inputs_root)) {
  stop("Missing R input bundles. Run prepare_benchmark_bank.py first.")
}
if (!is.null(dataset_ids)) {
  missing_ids <- setdiff(dataset_ids, manifest$dataset_id)
  if (length(missing_ids) > 0L) {
    stop("Unknown dataset IDs: ", paste(missing_ids, collapse = ", "))
  }
  manifest <- manifest[manifest$dataset_id %in% dataset_ids, , drop = FALSE]
}
if (!is.null(max_datasets)) {
  manifest <- head(manifest, max_datasets)
}
if (nrow(manifest) == 0L) {
  stop("No datasets selected.")
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
per_dataset_dir <- file.path(results_dir, "per_dataset")
dir.create(per_dataset_dir, recursive = TRUE, showWarnings = FALSE)

Rcpp::sourceCpp(file.path(script_dir, "dagar_poisson_boundary_matched_mwg.cpp"))

config <- data.frame(
  benchmark_dir = portable_project_path(benchmark_dir, project_root),
  results_dir = portable_project_path(results_dir, project_root),
  dataset_ids = if (is.null(dataset_ids)) "" else paste(dataset_ids, collapse = ","),
  max_datasets = if (is.null(max_datasets)) NA_integer_ else max_datasets,
  n_iter = n_iter,
  burnin = burnin,
  thin = thin,
  n_adapt = n_adapt,
  chains = 1L,
  seed = seed,
  threshold = threshold,
  beta0_prior_mean = 0.0,
  beta0_prior_sd = 0.5,
  sigma2_w_prior = "HalfNormal(scale=0.5) on sigma2_w",
  eta_raw_prior = "Uniform(0,1)",
  rho_prior = "Uniform(0,1)",
  isolate_repair = TRUE,
  latent_construction = "u~DAGAR; w=u-mean(u)",
  save_draws = save_draws,
  stringsAsFactors = FALSE
)
write.csv(config, file.path(results_dir, "matched_mcmc_config.csv"), row.names = FALSE)

all_parameter_summaries <- list()
all_chain_diagnostics <- list()
all_acceptance <- list()
all_runtime <- list()
all_edge_metrics <- list()

for (row_index in seq_len(nrow(manifest))) {
  row <- manifest[row_index, , drop = FALSE]
  dataset_id <- as.character(row$dataset_id[[1L]])
  dataset_dir <- file.path(r_inputs_root, dataset_id)
  node_table <- read.csv(file.path(dataset_dir, "node_table.csv"), check.names = FALSE)
  edge_table <- read.csv(file.path(dataset_dir, "edge_table.csv"), check.names = FALSE)
  adjacency <- as.matrix(read.csv(file.path(dataset_dir, "A.csv"), header = FALSE, check.names = FALSE))
  dissimilarity <- as.matrix(read.csv(file.path(dataset_dir, "Z.csv"), header = FALSE, check.names = FALSE))
  ordering <- as.integer(node_table$ordering_1based)
  m_value <- as.numeric(row$M[[1L]])

  message(sprintf(
    "[%d/%d] Running matched single-chain MCMC for %s (N=%d, edges=%d)",
    row_index,
    nrow(manifest),
    dataset_id,
    row$N[[1L]],
    row$edge_count[[1L]]
  ))
  set.seed(seed + row_index * 1000L)
  start_time <- proc.time()[["elapsed"]]
  fit <- dagar_poisson_boundary_matched_mwg(
    y_r = as.numeric(node_table$y),
    e_r = as.numeric(node_table$e),
    A_r = adjacency,
    Z_r = dissimilarity,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    n_adapt = n_adapt,
    ordering_r = ordering,
    eta_max = m_value,
    beta0_prior_mean = 0.0,
    beta0_prior_sd = 0.5,
    sigma2_w_prior_sd = 0.5,
    threshold = threshold,
    save_loglike = FALSE,
    verbose = verbose_sampler
  )
  elapsed <- proc.time()[["elapsed"]] - start_time

  eta_draws <- as.numeric(fit$eta)
  draws <- data.frame(
    draw = seq_along(fit$beta0),
    beta0 = as.numeric(fit$beta0),
    sigma2_w = as.numeric(fit$sigma2_w),
    eta_raw = eta_draws / m_value,
    eta = eta_draws,
    rho = as.numeric(fit$rho),
    stringsAsFactors = FALSE
  )
  truths <- c(
    beta0 = as.numeric(row$beta0_true[[1L]]),
    sigma2_w = as.numeric(row$sigma2_w_true[[1L]]),
    eta_raw = as.numeric(row$eta_raw_true[[1L]]),
    eta = as.numeric(row$eta_true[[1L]]),
    rho = as.numeric(row$rho_true[[1L]])
  )
  parameter_summaries <- do.call(
    rbind,
    lapply(names(truths), function(parameter) {
      summarize_draws(draws[[parameter]], truths[[parameter]], dataset_id, parameter)
    })
  )
  chain_diagnostics <- do.call(
    rbind,
    lapply(names(truths), function(parameter) {
      diagnostic <- estimate_chain_precision(draws[[parameter]], parameter)
      diagnostic$dataset_id <- dataset_id
      diagnostic[, c("dataset_id", setdiff(names(diagnostic), "dataset_id")), drop = FALSE]
    })
  )
  acceptance <- data.frame(
    dataset_id = dataset_id,
    parameter = names(fit$acceptance),
    acceptance_rate = as.numeric(fit$acceptance),
    stringsAsFactors = FALSE
  )
  runtime <- data.frame(
    dataset_id = dataset_id,
    chain = 1L,
    elapsed_sec = elapsed,
    n_saved_draws = nrow(draws),
    seconds_per_saved_draw = elapsed / max(1L, nrow(draws)),
    seconds_per_1000_saved_draws = 1000.0 * elapsed / max(1L, nrow(draws)),
    stringsAsFactors = FALSE
  )

  boundary_probability <- compute_boundary_probabilities(eta_draws, edge_table$edge_z, threshold)
  boundary_median <- as.integer(boundary_probability > 0.5)
  boundary_truth <- as.integer(edge_table$boundary_true)
  boundary_count_draws <- vapply(
    eta_draws,
    function(eta) sum(edge_table$edge_z * eta > threshold),
    numeric(1)
  )
  boundary_count_interval <- as.numeric(quantile(boundary_count_draws, c(0.025, 0.975)))
  edge_probabilities <- edge_table
  edge_probabilities$dataset_id <- dataset_id
  edge_probabilities$boundary_prob_mcmc <- boundary_probability
  edge_probabilities$boundary_median_mcmc <- boundary_median
  edge_metrics <- data.frame(
    dataset_id = dataset_id,
    edge_count = nrow(edge_table),
    true_boundary_count = sum(boundary_truth),
    posterior_boundary_count_mpm = sum(boundary_median),
    auroc = compute_auc(boundary_probability, boundary_truth),
    average_precision = compute_average_precision(boundary_probability, boundary_truth),
    brier = mean((boundary_probability - boundary_truth) ^ 2),
    sensitivity_mpm = if (sum(boundary_truth == 1) > 0L) {
      mean(boundary_median[boundary_truth == 1] == 1)
    } else {
      NA_real_
    },
    specificity_mpm = if (sum(boundary_truth == 0) > 0L) {
      mean(boundary_median[boundary_truth == 0] == 0)
    } else {
      NA_real_
    },
    boundary_count_mean_draws = mean(boundary_count_draws),
    boundary_count_lower_95 = boundary_count_interval[[1L]],
    boundary_count_upper_95 = boundary_count_interval[[2L]],
    boundary_count_truth_in_95 = as.integer(
      boundary_count_interval[[1L]] <= sum(boundary_truth) &&
        sum(boundary_truth) <= boundary_count_interval[[2L]]
    ),
    stringsAsFactors = FALSE
  )

  dataset_output <- file.path(per_dataset_dir, dataset_id)
  dir.create(dataset_output, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    parameter_summaries,
    file.path(dataset_output, paste0(dataset_id, "_parameter_summary.csv")),
    row.names = FALSE
  )
  write.csv(
    edge_probabilities,
    file.path(dataset_output, paste0(dataset_id, "_edge_probabilities.csv")),
    row.names = FALSE
  )
  write.csv(
    edge_metrics,
    file.path(dataset_output, paste0(dataset_id, "_edge_metrics.csv")),
    row.names = FALSE
  )
  if (save_draws) {
    write_gzip_csv(
      draws,
      file.path(dataset_output, paste0(dataset_id, "_posterior_draws.csv.gz"))
    )
  }

  all_parameter_summaries[[row_index]] <- parameter_summaries
  all_chain_diagnostics[[row_index]] <- chain_diagnostics
  all_acceptance[[row_index]] <- acceptance
  all_runtime[[row_index]] <- runtime
  all_edge_metrics[[row_index]] <- edge_metrics
}

write.csv(
  do.call(rbind, all_parameter_summaries),
  file.path(results_dir, "combined_parameter_summaries.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, all_chain_diagnostics),
  file.path(results_dir, "combined_chain_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, all_acceptance),
  file.path(results_dir, "combined_acceptance.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, all_runtime),
  file.path(results_dir, "combined_runtime.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, all_edge_metrics),
  file.path(results_dir, "combined_edge_metrics.csv"),
  row.names = FALSE
)

message("\nSaved matched MCMC benchmark outputs to: ", results_dir)
