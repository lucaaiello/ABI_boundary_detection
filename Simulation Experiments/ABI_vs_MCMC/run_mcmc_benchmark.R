rm(list = ls())

options(stringsAsFactors = FALSE)

script_config <- list(
  benchmark_dir = "C:/Dati/Lavori/Aiello_Banerjee_2026/Code/ABI_poisson_regression/Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100",
  results_dir = "C:/Dati/Lavori/Aiello_Banerjee_2026/Code/ABI_poisson_regression/Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100/mcmc_results_all100",
  dataset_ids = NULL,
  max_datasets = NULL,
  n_iter = 20000L,
  burnin = 10000L,
  thin = 1L,
  n_adapt = 10000L,
  chains = 1L,
  seed = 123L,
  threshold = log(2.0),
  save_draws = TRUE,
  verbose_sampler = FALSE,
  use_null_makevars = FALSE
)


parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop("Unexpected argument: ", arg)
    }
    key <- sub("^--", "", arg)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[parts[1]]] <- parts[2]
      i <- i + 1
      next
    }
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- "TRUE"
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

get_script_dir <- function() {
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

find_existing_parent <- function(start_dir, target_file, max_up = 6L) {
  if (is.null(start_dir)) {
    return(NULL)
  }
  current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (step in seq_len(max_up + 1L)) {
    candidate <- file.path(current, target_file)
    if (file.exists(candidate)) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  NULL
}

resolve_project_root <- function(script_dir, benchmark_dir) {
  target_file <- "ABI_poisson_regression.code-workspace"
  candidates <- c(
    find_existing_parent(script_dir, target_file),
    find_existing_parent(benchmark_dir, target_file),
    find_existing_parent(getwd(), target_file)
  )
  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) == 0) {
    stop(
      "Could not locate ", target_file,
      ". Checked upward from the script directory, benchmark directory, and current working directory."
    )
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

as_int <- function(x, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  as.integer(x)
}

as_num <- function(x, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  as.numeric(x)
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "1", "yes", "y")
}

pick_value <- function(args, config, key, default = NULL) {
  cli_value <- args[[key]]
  if (!is.null(cli_value)) {
    return(cli_value)
  }
  config_key <- gsub("-", "_", key, fixed = TRUE)
  config_value <- config[[config_key]]
  if (!is.null(config_value)) {
    return(config_value)
  }
  default
}

as_dataset_ids <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  if (length(x) == 1) {
    pieces <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  } else {
    pieces <- as.character(x)
  }
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces)]
  if (length(pieces) == 0) {
    return(NULL)
  }
  pieces
}

summarize_draws <- function(samples, truth, dataset_id, parameter) {
  q <- as.numeric(quantile(samples, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  mean_val <- mean(samples)
  data.frame(
    dataset_id = dataset_id,
    parameter = parameter,
    posterior_mean = mean_val,
    posterior_sd = sd(samples),
    posterior_median = q[2],
    lower_95 = q[1],
    upper_95 = q[3],
    truth = truth,
    bias_mean = mean_val - truth,
    abs_error_mean = abs(mean_val - truth),
    covered_95 = as.integer(q[1] <= truth && truth <= q[3]),
    stringsAsFactors = FALSE
  )
}

compute_rhat <- function(chain_list) {
  if (length(chain_list) < 2) {
    return(NA_real_)
  }
  n <- min(vapply(chain_list, length, integer(1)))
  if (n < 2) {
    return(NA_real_)
  }
  chains <- lapply(chain_list, function(x) x[seq_len(n)])
  m <- length(chains)
  chain_means <- vapply(chains, mean, numeric(1))
  chain_vars <- vapply(chains, var, numeric(1))
  w <- mean(chain_vars)
  if (!is.finite(w) || w <= 0) {
    return(NA_real_)
  }
  b <- n * var(chain_means)
  var_hat <- ((n - 1) / n) * w + (b / n)
  sqrt(var_hat / w)
}

compute_auc <- function(prob, truth) {
  pos <- which(truth == 1)
  neg <- which(truth == 0)
  if (length(pos) == 0 || length(neg) == 0) {
    return(NA_real_)
  }
  ranks <- rank(prob, ties.method = "average")
  (sum(ranks[pos]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

compute_average_precision <- function(prob, truth) {
  pos_total <- sum(truth == 1)
  if (pos_total == 0) {
    return(NA_real_)
  }
  ord <- order(prob, decreasing = TRUE)
  truth_ord <- truth[ord]
  tp <- cumsum(truth_ord == 1)
  precision <- tp / seq_along(tp)
  sum(precision[truth_ord == 1]) / pos_total
}

compute_brier <- function(prob, truth) {
  mean((prob - truth) ^ 2)
}

compute_boundary_probabilities <- function(eta_draws, edge_z, threshold) {
  draw_mat <- outer(eta_draws, edge_z, FUN = function(eta, z) as.numeric(z * eta > threshold))
  colMeans(draw_mat)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

benchmark_dir <- pick_value(args, script_config, "benchmark-dir")
if (is.null(benchmark_dir)) {
  stop("Please provide --benchmark-dir, or set script_config$benchmark_dir at the top of this file.")
}
benchmark_dir <- normalizePath(benchmark_dir, mustWork = TRUE)
project_root <- resolve_project_root(get_script_dir(), benchmark_dir)

results_dir <- pick_value(args, script_config, "results-dir")
if (is.null(results_dir)) {
  results_dir <- file.path(benchmark_dir, "mcmc_results")
}
results_dir <- normalizePath(results_dir, mustWork = FALSE)

if (as_flag(pick_value(args, script_config, "use-null-makevars"), default = FALSE)) {
  Sys.setenv(R_MAKEVARS_USER = "NUL")
}

n_iter <- as_int(pick_value(args, script_config, "n-iter"), default = 4000L)
burnin <- as_int(pick_value(args, script_config, "burnin"), default = 2000L)
thin <- as_int(pick_value(args, script_config, "thin"), default = 2L)
n_adapt <- as_int(pick_value(args, script_config, "n-adapt"), default = 2000L)
chains <- as_int(pick_value(args, script_config, "chains"), default = 4L)
seed <- as_int(pick_value(args, script_config, "seed"), default = 123L)
threshold <- as_num(pick_value(args, script_config, "threshold"), default = log(2.0))
save_draws <- as_flag(pick_value(args, script_config, "save-draws"), default = TRUE)
verbose_sampler <- as_flag(pick_value(args, script_config, "verbose-sampler"), default = FALSE)

if (n_iter <= 0 || burnin < 0 || thin <= 0 || n_adapt < 0 || chains <= 0) {
  stop("Invalid MCMC settings.")
}
if (n_iter <= burnin) {
  stop("n_iter must be larger than burnin so that posterior draws are saved.")
}

manifest_path <- file.path(benchmark_dir, "benchmark_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Could not find manifest: ", manifest_path)
}
manifest <- read.csv(manifest_path, check.names = FALSE)
r_inputs_root <- file.path(benchmark_dir, "r_inputs")
if (!dir.exists(r_inputs_root)) {
  stop(
    "This benchmark bank does not contain the required 'r_inputs/' directory.\n",
    "It was likely exported with an older version of export_benchmark_datasets.py.\n",
    "Please rerun the exporter on this bank, for example:\n",
    "python ABI_vs_MCMC/export_benchmark_datasets.py --num-datasets ",
    nrow(manifest),
    " --seed ",
    seed,
    " --output-dir \"",
    benchmark_dir,
    "\""
  )
}

dataset_ids <- as_dataset_ids(pick_value(args, script_config, "dataset-ids"))
if (!is.null(dataset_ids)) {
  manifest <- manifest[manifest$dataset_id %in% dataset_ids, , drop = FALSE]
}

max_datasets <- as_int(pick_value(args, script_config, "max-datasets"), default = NULL)
if (!is.null(max_datasets) && nrow(manifest) > max_datasets) {
  manifest <- manifest[seq_len(max_datasets), , drop = FALSE]
}

if (nrow(manifest) == 0) {
  stop("No datasets selected.")
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
per_dataset_dir <- file.path(results_dir, "per_dataset")
dir.create(per_dataset_dir, recursive = TRUE, showWarnings = FALSE)

config_df <- data.frame(
  benchmark_dir = benchmark_dir,
  results_dir = results_dir,
  dataset_ids = if (is.null(dataset_ids)) "" else paste(dataset_ids, collapse = ","),
  max_datasets = if (is.null(max_datasets)) NA_integer_ else max_datasets,
  n_iter = n_iter,
  burnin = burnin,
  thin = thin,
  n_adapt = n_adapt,
  chains = chains,
  seed = seed,
  threshold = threshold,
  save_draws = save_draws,
  stringsAsFactors = FALSE
)
write.csv(config_df, file.path(results_dir, "mcmc_config.csv"), row.names = FALSE)

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required.")
}

Rcpp::sourceCpp(file.path(project_root, "dagar_poisson_boundary_mwg.cpp"))

all_parameter_summaries <- list()
all_chain_diagnostics <- list()
all_acceptance <- list()
all_runtime <- list()
all_runtime_summary <- list()
all_edge_metrics <- list()

for (row_idx in seq_len(nrow(manifest))) {
  row <- manifest[row_idx, , drop = FALSE]
  dataset_id <- row$dataset_id[[1]]
  dataset_dir <- file.path(r_inputs_root, dataset_id)
  if (!dir.exists(dataset_dir)) {
    stop(
      "Missing R input bundle for ", dataset_id, ": ", dataset_dir, "\n",
      "Please rerun the exporter so that r_inputs/ is created for every dataset."
    )
  }

  node_table <- read.csv(file.path(dataset_dir, "node_table.csv"), check.names = FALSE)
  edge_table <- read.csv(file.path(dataset_dir, "edge_table.csv"), check.names = FALSE)
  a_mat <- as.matrix(read.csv(file.path(dataset_dir, "A.csv"), header = FALSE, check.names = FALSE))
  z_mat <- as.matrix(read.csv(file.path(dataset_dir, "Z.csv"), header = FALSE, check.names = FALSE))

  y <- as.numeric(node_table$y)
  e <- as.numeric(node_table$e)
  m_val <- as.numeric(row$M[[1]])

  chain_draws <- vector("list", chains)
  acceptance_rows <- vector("list", chains)
  runtime_rows <- vector("list", chains)

  message(sprintf("[%d/%d] Running CARBayes-style DAGAR MCMC for %s (N=%d, edges=%d)", row_idx, nrow(manifest), dataset_id, row$N[[1]], row$edge_count[[1]]))
  dataset_start_time <- proc.time()[["elapsed"]]

  for (chain_idx in seq_len(chains)) {
    chain_seed <- seed + row_idx * 1000L + chain_idx
    set.seed(chain_seed)

    start_time <- proc.time()[["elapsed"]]
    fit <- dagar_poisson_boundary_mwg2(
      y_r = y,
      e_r = e,
      A_r = a_mat,
      Z_r = z_mat,
      n_iter = n_iter,
      burnin = burnin,
      thin = thin,
      n_adapt = n_adapt,
      alpha_max = m_val,
      save_phi = FALSE,
      save_fitted = FALSE,
      save_loglike = TRUE,
      verbose = verbose_sampler
    )
    elapsed <- proc.time()[["elapsed"]] - start_time

    chain_loglik <- if (!is.null(fit$samples_loglike)) {
      rowSums(as.matrix(fit$samples_loglike))
    } else {
      rep(NA_real_, length(fit$beta0))
    }

    eta_draws <- as.numeric(fit$eta)
    eta_raw_draws <- eta_draws / m_val

    draw_df <- data.frame(
      dataset_id = dataset_id,
      chain = chain_idx,
      draw = seq_along(fit$beta0),
      beta0 = as.numeric(fit$beta0),
      sigma2_w = as.numeric(fit$sigma2_w),
      eta_raw = eta_raw_draws,
      eta = eta_draws,
      rho = as.numeric(fit$rho),
      loglik = as.numeric(chain_loglik),
      stringsAsFactors = FALSE
    )
    chain_draws[[chain_idx]] <- draw_df

    accept_vec <- as.numeric(fit$accept)
    beta_accept_rate <- if (length(accept_vec) >= 2 && accept_vec[2] > 0) accept_vec[1] / accept_vec[2] else NA_real_
    phi_accept_rate <- if (length(accept_vec) >= 4 && accept_vec[4] > 0) accept_vec[3] / accept_vec[4] else NA_real_
    eta_accept_rate <- if (length(accept_vec) >= 6 && accept_vec[6] > 0) accept_vec[5] / accept_vec[6] else NA_real_
    rho_accept_rate <- if (length(accept_vec) >= 8 && accept_vec[8] > 0) accept_vec[7] / accept_vec[8] else NA_real_

    acceptance_rows[[chain_idx]] <- data.frame(
      dataset_id = dataset_id,
      chain = chain_idx,
      beta0_accept = beta_accept_rate,
      sigma2_w_accept = NA_real_,
      eta_accept = eta_accept_rate,
      rho_accept = rho_accept_rate,
      phi_accept = phi_accept_rate,
      w_accept_mean = NA_real_,
      w_accept_median = NA_real_,
      stringsAsFactors = FALSE
    )

    runtime_rows[[chain_idx]] <- data.frame(
      dataset_id = dataset_id,
      chain = chain_idx,
      elapsed_sec = elapsed,
      n_saved_draws = nrow(draw_df),
      stringsAsFactors = FALSE
    )
  }

  draws_df <- do.call(rbind, chain_draws)
  acceptance_df <- do.call(rbind, acceptance_rows)
  runtime_df <- do.call(rbind, runtime_rows)
  dataset_wall_elapsed_sec <- proc.time()[["elapsed"]] - dataset_start_time
  chain_elapsed_sum_sec <- sum(runtime_df$elapsed_sec)
  runtime_summary_df <- data.frame(
    dataset_id = dataset_id,
    n_chains = chains,
    n_saved_draws_total = nrow(draws_df),
    n_saved_draws_per_chain = nrow(chain_draws[[1]]),
    chain_elapsed_sum_sec = chain_elapsed_sum_sec,
    chain_elapsed_mean_sec = mean(runtime_df$elapsed_sec),
    chain_elapsed_min_sec = min(runtime_df$elapsed_sec),
    chain_elapsed_max_sec = max(runtime_df$elapsed_sec),
    dataset_wall_elapsed_sec = dataset_wall_elapsed_sec,
    postprocess_overhead_sec = dataset_wall_elapsed_sec - chain_elapsed_sum_sec,
    seconds_per_saved_draw = dataset_wall_elapsed_sec / max(1L, nrow(draws_df)),
    seconds_per_1000_saved_draws = 1000.0 * dataset_wall_elapsed_sec / max(1L, nrow(draws_df)),
    stringsAsFactors = FALSE
  )

  parameter_summaries <- do.call(
    rbind,
    list(
      summarize_draws(draws_df$beta0, as.numeric(row$beta0_true[[1]]), dataset_id, "beta0"),
      summarize_draws(draws_df$sigma2_w, as.numeric(row$sigma2_w_true[[1]]), dataset_id, "sigma2_w"),
      summarize_draws(draws_df$eta_raw, as.numeric(row$eta_raw_true[[1]]), dataset_id, "eta_raw"),
      summarize_draws(draws_df$eta, as.numeric(row$eta_true[[1]]), dataset_id, "eta"),
      summarize_draws(draws_df$rho, as.numeric(row$rho_true[[1]]), dataset_id, "rho")
    )
  )

  chain_diag_df <- data.frame(
    dataset_id = dataset_id,
    parameter = c("beta0", "sigma2_w", "eta_raw", "eta", "rho"),
    rhat = c(
      compute_rhat(lapply(chain_draws, function(x) x$beta0)),
      compute_rhat(lapply(chain_draws, function(x) x$sigma2_w)),
      compute_rhat(lapply(chain_draws, function(x) x$eta_raw)),
      compute_rhat(lapply(chain_draws, function(x) x$eta)),
      compute_rhat(lapply(chain_draws, function(x) x$rho))
    ),
    n_chains = chains,
    n_draws_per_chain = nrow(chain_draws[[1]]),
    stringsAsFactors = FALSE
  )

  boundary_prob <- compute_boundary_probabilities(draws_df$eta, edge_table$edge_z, threshold = threshold)
  edge_prob_df <- edge_table
  edge_prob_df$dataset_id <- dataset_id
  edge_prob_df$boundary_prob_mcmc <- boundary_prob
  edge_prob_df$boundary_median_mcmc <- as.integer(boundary_prob > 0.5)

  truth_boundary <- as.numeric(edge_prob_df$boundary_true)
  boundary_count_draws <- vapply(draws_df$eta, function(eta) sum(edge_table$edge_z * eta > threshold), numeric(1))
  boundary_count_q <- as.numeric(quantile(boundary_count_draws, c(0.025, 0.975)))
  edge_metrics_df <- data.frame(
    dataset_id = dataset_id,
    edge_count = nrow(edge_prob_df),
    true_boundary_count = sum(truth_boundary == 1),
    posterior_boundary_count_mpm = sum(edge_prob_df$boundary_median_mcmc == 1),
    auroc = compute_auc(edge_prob_df$boundary_prob_mcmc, truth_boundary),
    average_precision = compute_average_precision(edge_prob_df$boundary_prob_mcmc, truth_boundary),
    brier = compute_brier(edge_prob_df$boundary_prob_mcmc, truth_boundary),
    sensitivity_mpm = if (sum(truth_boundary == 1) > 0) {
      mean(edge_prob_df$boundary_median_mcmc[truth_boundary == 1] == 1)
    } else {
      NA_real_
    },
    specificity_mpm = if (sum(truth_boundary == 0) > 0) {
      mean(edge_prob_df$boundary_median_mcmc[truth_boundary == 0] == 0)
    } else {
      NA_real_
    },
    boundary_count_mean_draws = mean(boundary_count_draws),
    boundary_count_lower_95 = boundary_count_q[1],
    boundary_count_upper_95 = boundary_count_q[2],
    boundary_count_truth_in_95 = as.integer(
      sum(truth_boundary == 1) >= boundary_count_q[1] &&
      sum(truth_boundary == 1) <= boundary_count_q[2]
    ),
    stringsAsFactors = FALSE
  )

  dataset_out_dir <- file.path(per_dataset_dir, dataset_id)
  dir.create(dataset_out_dir, recursive = TRUE, showWarnings = FALSE)

  if (save_draws) {
    write.csv(draws_df, file.path(dataset_out_dir, paste0(dataset_id, "_posterior_draws.csv")), row.names = FALSE)
  }
  write.csv(parameter_summaries, file.path(dataset_out_dir, paste0(dataset_id, "_parameter_summary.csv")), row.names = FALSE)
  write.csv(chain_diag_df, file.path(dataset_out_dir, paste0(dataset_id, "_chain_diagnostics.csv")), row.names = FALSE)
  write.csv(acceptance_df, file.path(dataset_out_dir, paste0(dataset_id, "_acceptance.csv")), row.names = FALSE)
  write.csv(runtime_df, file.path(dataset_out_dir, paste0(dataset_id, "_runtime.csv")), row.names = FALSE)
  write.csv(runtime_summary_df, file.path(dataset_out_dir, paste0(dataset_id, "_runtime_summary.csv")), row.names = FALSE)
  write.csv(edge_prob_df, file.path(dataset_out_dir, paste0(dataset_id, "_edge_probabilities.csv")), row.names = FALSE)
  write.csv(edge_metrics_df, file.path(dataset_out_dir, paste0(dataset_id, "_edge_metrics.csv")), row.names = FALSE)

  all_parameter_summaries[[length(all_parameter_summaries) + 1L]] <- parameter_summaries
  all_chain_diagnostics[[length(all_chain_diagnostics) + 1L]] <- chain_diag_df
  all_acceptance[[length(all_acceptance) + 1L]] <- acceptance_df
  all_runtime[[length(all_runtime) + 1L]] <- runtime_df
  all_runtime_summary[[length(all_runtime_summary) + 1L]] <- runtime_summary_df
  all_edge_metrics[[length(all_edge_metrics) + 1L]] <- edge_metrics_df
}

write.csv(do.call(rbind, all_parameter_summaries), file.path(results_dir, "combined_parameter_summaries.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_chain_diagnostics), file.path(results_dir, "combined_chain_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_acceptance), file.path(results_dir, "combined_acceptance.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_runtime), file.path(results_dir, "combined_runtime.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_runtime_summary), file.path(results_dir, "combined_runtime_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_edge_metrics), file.path(results_dir, "combined_edge_metrics.csv"), row.names = FALSE)

message("\nSaved MCMC benchmark outputs to: ", results_dir)
