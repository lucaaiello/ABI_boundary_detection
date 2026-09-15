#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

struct DagarStructure {
  std::vector<double> lambda;
  std::vector<double> b_row;
  std::vector<std::vector<int> > preds;
  std::vector<std::vector<int> > succs;
  double sum_log_lambda;
};

inline int flat_index(const int i, const int j, const int n) {
  return i * n + j;
}

inline double inv_logit(const double x) {
  if (x >= 0.0) {
    const double e = std::exp(-x);
    return 1.0 / (1.0 + e);
  }
  const double e = std::exp(x);
  return e / (1.0 + e);
}

inline double clip_open_interval(const double x, const double lower, const double upper) {
  const double eps = 1e-12;
  return std::min(std::max(x, lower + eps), upper - eps);
}

inline bool accept_log_mh(const double log_ratio) {
  return std::log(R::runif(0.0, 1.0)) < log_ratio;
}

inline double log_half_normal_kernel(const double x, const double scale) {
  if (!R_FINITE(x) || x <= 0.0 || !R_FINITE(scale) || scale <= 0.0) {
    return R_NegInf;
  }
  const double standardized = x / scale;
  return -0.5 * standardized * standardized;
}

inline double truncnorm_log_density(
    const double x,
    const double mean,
    const double sd,
    const double lower,
    const double upper
) {
  if (x <= lower || x >= upper || sd <= 0.0) {
    return R_NegInf;
  }
  const double z_lower = (lower - mean) / sd;
  const double z_upper = (upper - mean) / sd;
  const double norm_const =
    R::pnorm5(z_upper, 0.0, 1.0, 1, 0) - R::pnorm5(z_lower, 0.0, 1.0, 1, 0);
  if (!R_FINITE(norm_const) || norm_const <= 0.0) {
    return R_NegInf;
  }
  return R::dnorm4(x, mean, sd, 1) - std::log(norm_const);
}

inline double rtruncnorm_one(
    const double mean,
    const double sd,
    const double lower,
    const double upper
) {
  const double cdf_lower = R::pnorm5((lower - mean) / sd, 0.0, 1.0, 1, 0);
  const double cdf_upper = R::pnorm5((upper - mean) / sd, 0.0, 1.0, 1, 0);
  const double u = R::runif(cdf_lower, cdf_upper);
  return mean + sd * R::qnorm5(u, 0.0, 1.0, 1, 0);
}

inline double tune_sd(
    const double current,
    const int accepted,
    const int attempted,
    const double lower_target,
    const double upper_target,
    const double upper_bound = R_PosInf
) {
  if (attempted == 0) {
    return current;
  }
  const double rate = static_cast<double>(accepted) / static_cast<double>(attempted);
  double updated = current;
  if (rate < lower_target) {
    updated *= 0.8;
  } else if (rate > upper_target) {
    updated *= 1.2;
  }
  updated = std::max(updated, 1e-8);
  if (R_FINITE(upper_bound)) {
    updated = std::min(updated, upper_bound);
  }
  return updated;
}

void validate_square_matrix(const NumericMatrix& matrix, const char* name) {
  if (matrix.nrow() != matrix.ncol()) {
    stop("%s must be square.", name);
  }
}

std::vector<int> parse_ordering(
    const int n,
    const Nullable<IntegerVector>& ordering_r
) {
  std::vector<int> ordering(n);
  if (ordering_r.isNull()) {
    std::iota(ordering.begin(), ordering.end(), 0);
    return ordering;
  }

  const IntegerVector input(ordering_r);
  if (input.size() != n) {
    stop("ordering must have length n.");
  }
  std::vector<int> seen(n, 0);
  for (int position = 0; position < n; ++position) {
    const int index = input[position] - 1;
    if (index < 0 || index >= n || seen[index]) {
      stop("ordering must be a 1-based permutation of 1:n.");
    }
    seen[index] = 1;
    ordering[position] = index;
  }
  return ordering;
}

void build_original_graph(
    const NumericMatrix& adjacency,
    std::vector<unsigned char>& graph,
    std::vector<std::pair<int, int> >& edges
) {
  const int n = adjacency.nrow();
  graph.assign(n * n, 0);
  edges.clear();
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (adjacency(i, j) > 0.5 || adjacency(j, i) > 0.5) {
        graph[flat_index(i, j, n)] = 1;
        graph[flat_index(j, i, n)] = 1;
        edges.push_back(std::make_pair(i, j));
      }
    }
  }
}

void build_filtered_graph_with_repair(
    const int n,
    const NumericMatrix& dissimilarity,
    const std::vector<unsigned char>& original_graph,
    const std::vector<std::pair<int, int> >& edges,
    const double eta,
    const double threshold,
    std::vector<unsigned char>& filtered_graph
) {
  filtered_graph.assign(n * n, 0);
  std::vector<int> degree(n, 0);
  for (std::size_t edge = 0; edge < edges.size(); ++edge) {
    const int i = edges[edge].first;
    const int j = edges[edge].second;
    if (dissimilarity(i, j) * eta <= threshold) {
      filtered_graph[flat_index(i, j, n)] = 1;
      filtered_graph[flat_index(j, i, n)] = 1;
      ++degree[i];
      ++degree[j];
    }
  }

  // Match the simulator: scan nodes in index order and reconnect an isolate
  // to its least-dissimilar original neighbor.
  for (int i = 0; i < n; ++i) {
    if (degree[i] != 0) {
      continue;
    }
    int best_neighbor = -1;
    double best_dissimilarity = std::numeric_limits<double>::infinity();
    for (int j = 0; j < n; ++j) {
      if (!original_graph[flat_index(i, j, n)]) {
        continue;
      }
      const double candidate = dissimilarity(i, j);
      if (candidate < best_dissimilarity) {
        best_dissimilarity = candidate;
        best_neighbor = j;
      }
    }
    if (best_neighbor >= 0) {
      filtered_graph[flat_index(i, best_neighbor, n)] = 1;
      filtered_graph[flat_index(best_neighbor, i, n)] = 1;
      ++degree[i];
      ++degree[best_neighbor];
    }
  }
}

DagarStructure build_dagar_structure(
    const std::vector<unsigned char>& graph,
    const std::vector<int>& ordering,
    const double rho
) {
  const int n = static_cast<int>(ordering.size());
  const double rho_squared = rho * rho;
  DagarStructure output;
  output.lambda.assign(n, 0.0);
  output.b_row.assign(n, 0.0);
  output.preds.assign(n, std::vector<int>());
  output.succs.assign(n, std::vector<int>());
  output.sum_log_lambda = 0.0;

  for (int position = 0; position < n; ++position) {
    const int i = ordering[position];
    for (int previous = 0; previous < position; ++previous) {
      const int j = ordering[previous];
      if (graph[flat_index(i, j, n)]) {
        output.preds[i].push_back(j);
      }
    }
    const int predecessor_count = static_cast<int>(output.preds[i].size());
    const double denominator = 1.0 + std::max(predecessor_count - 1, 0) * rho_squared;
    output.b_row[i] = predecessor_count > 0 ? rho / denominator : 0.0;
    output.lambda[i] = denominator / (1.0 - rho_squared);
    output.sum_log_lambda += std::log(output.lambda[i]);
    for (std::size_t index = 0; index < output.preds[i].size(); ++index) {
      output.succs[output.preds[i][index]].push_back(i);
    }
  }
  return output;
}

double compute_quadratic_form(
    const std::vector<double>& field,
    const DagarStructure& dagar,
    std::vector<double>& residuals
) {
  const int n = static_cast<int>(field.size());
  residuals.assign(n, 0.0);
  double quadratic = 0.0;
  for (int i = 0; i < n; ++i) {
    double predecessor_sum = 0.0;
    for (std::size_t index = 0; index < dagar.preds[i].size(); ++index) {
      predecessor_sum += field[dagar.preds[i][index]];
    }
    residuals[i] = field[i] - dagar.b_row[i] * predecessor_sum;
    quadratic += dagar.lambda[i] * residuals[i] * residuals[i];
  }
  return quadratic;
}

inline double clipped_poisson_log_kernel(
    const double log_rate,
    const double y,
    const double log_lambda_lower,
    const double log_lambda_upper,
    double& lambda,
    double& clipped_log_rate
) {
  clipped_log_rate = std::min(std::max(log_rate, log_lambda_lower), log_lambda_upper);
  lambda = std::exp(clipped_log_rate);
  return y * clipped_log_rate - lambda;
}

double compute_loglikelihood(
    const double beta0,
    const std::vector<double>& field,
    const double field_mean,
    const std::vector<double>& log_exposure,
    const std::vector<double>& y,
    const double log_lambda_lower,
    const double log_lambda_upper,
    std::vector<double>& lambda,
    std::vector<double>& clipped_log_rate
) {
  const int n = static_cast<int>(field.size());
  double loglikelihood = 0.0;
  for (int i = 0; i < n; ++i) {
    const double log_rate = log_exposure[i] + beta0 + field[i] - field_mean;
    loglikelihood += clipped_poisson_log_kernel(
      log_rate,
      y[i],
      log_lambda_lower,
      log_lambda_upper,
      lambda[i],
      clipped_log_rate[i]
    );
  }
  return loglikelihood;
}

NumericMatrix graph_to_matrix(const std::vector<unsigned char>& graph, const int n) {
  NumericMatrix output(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      output(i, j) = static_cast<double>(graph[flat_index(i, j, n)]);
    }
  }
  return output;
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List dagar_poisson_boundary_matched_mwg(
    const NumericVector& y_r,
    const NumericVector& e_r,
    const NumericMatrix& A_r,
    const NumericMatrix& Z_r,
    const int n_iter,
    const int burnin = 10000,
    const int thin = 1,
    const int n_adapt = 10000,
    Nullable<IntegerVector> ordering_r = R_NilValue,
    const double beta0_init = NA_REAL,
    const double sigma2_w_init = 0.25,
    const double eta_init = NA_REAL,
    const double eta_max = NA_REAL,
    const double rho_init = 0.5,
    const double proposal_sd_beta0_init = 0.01,
    const double proposal_sd_field_init = 0.10,
    const double proposal_sd_log_sigma2_init = 0.10,
    const double proposal_sd_eta_init = NA_REAL,
    const double proposal_sd_logit_rho_init = 0.10,
    const double beta0_prior_mean = 0.0,
    const double beta0_prior_sd = 0.5,
    const double sigma2_w_prior_sd = 0.5,
    const double threshold = 0.6931471805599453,
    const double lambda_lower = 1e-2,
    const double lambda_upper = 1e6,
    const bool save_loglike = false,
    const bool verbose = false
) {
  if (n_iter <= burnin || burnin < 0 || thin <= 0 || n_adapt < 0) {
    stop("Require n_iter > burnin >= 0, thin > 0, and n_adapt >= 0.");
  }
  if (!R_FINITE(eta_max) || eta_max <= 0.0) {
    stop("eta_max must be finite and positive.");
  }
  if (beta0_prior_sd <= 0.0 || sigma2_w_prior_sd <= 0.0) {
    stop("Prior scales must be positive.");
  }
  if (lambda_lower <= 0.0 || lambda_upper <= lambda_lower) {
    stop("Require 0 < lambda_lower < lambda_upper.");
  }
  validate_square_matrix(A_r, "A");
  validate_square_matrix(Z_r, "Z");

  const int n = y_r.size();
  if (n < 2 || e_r.size() != n || A_r.nrow() != n || Z_r.nrow() != n) {
    stop("The data vectors and graph matrices have incompatible dimensions.");
  }

  std::vector<double> y(n);
  std::vector<double> log_exposure(n);
  double sum_y = 0.0;
  double sum_e = 0.0;
  for (int i = 0; i < n; ++i) {
    if (y_r[i] < 0.0 || e_r[i] <= 0.0) {
      stop("Counts must be non-negative and exposures must be positive.");
    }
    y[i] = y_r[i];
    log_exposure[i] = std::log(e_r[i]);
    sum_y += y_r[i];
    sum_e += e_r[i];
  }

  const std::vector<int> ordering = parse_ordering(n, ordering_r);
  std::vector<unsigned char> original_graph;
  std::vector<std::pair<int, int> > undirected_edges;
  build_original_graph(A_r, original_graph, undirected_edges);
  if (undirected_edges.empty()) {
    stop("The adjacency matrix contains no edges.");
  }

  double beta0 = NumericVector::is_na(beta0_init)
    ? std::log((sum_y + 0.5) / (sum_e + 0.5))
    : beta0_init;
  double sigma2_w = sigma2_w_init > 0.0 ? sigma2_w_init : 0.25;
  double log_sigma2_w = std::log(sigma2_w);
  double eta = NumericVector::is_na(eta_init)
    ? R::runif(0.0, eta_max / 3.0)
    : clip_open_interval(eta_init, 0.0, eta_max);
  double rho = clip_open_interval(rho_init, 0.0, 1.0);
  double logit_rho = std::log(rho) - std::log(1.0 - rho);

  std::vector<double> field(n, 0.0);
  double field_mean = 0.0;
  std::vector<unsigned char> filtered_graph;
  build_filtered_graph_with_repair(
    n, Z_r, original_graph, undirected_edges, eta, threshold, filtered_graph
  );
  DagarStructure dagar = build_dagar_structure(filtered_graph, ordering, rho);
  std::vector<double> residuals;
  double quadratic = compute_quadratic_form(field, dagar, residuals);

  const double log_lambda_lower = std::log(lambda_lower);
  const double log_lambda_upper = std::log(lambda_upper);
  std::vector<double> lambda(n, 0.0);
  std::vector<double> clipped_log_rate(n, 0.0);
  double loglikelihood = compute_loglikelihood(
    beta0,
    field,
    field_mean,
    log_exposure,
    y,
    log_lambda_lower,
    log_lambda_upper,
    lambda,
    clipped_log_rate
  );

  double proposal_sd_beta0 = std::max(proposal_sd_beta0_init, 1e-8);
  double proposal_sd_field = std::max(proposal_sd_field_init, 1e-8);
  double proposal_sd_log_sigma2 = std::max(proposal_sd_log_sigma2_init, 1e-8);
  double proposal_sd_eta = NumericVector::is_na(proposal_sd_eta_init)
    ? 0.02 * eta_max
    : std::max(proposal_sd_eta_init, 1e-8);
  double proposal_sd_logit_rho = std::max(proposal_sd_logit_rho_init, 1e-8);

  int beta_accept = 0, beta_attempt = 0, beta_window_accept = 0, beta_window_attempt = 0;
  int field_accept = 0, field_attempt = 0, field_window_accept = 0, field_window_attempt = 0;
  int sigma_accept = 0, sigma_attempt = 0, sigma_window_accept = 0, sigma_window_attempt = 0;
  int eta_accept = 0, eta_attempt = 0, eta_window_accept = 0, eta_window_attempt = 0;
  int rho_accept = 0, rho_attempt = 0, rho_window_accept = 0, rho_window_attempt = 0;

  const int n_keep = (n_iter - burnin) / thin;
  NumericVector samples_beta0(n_keep);
  NumericVector samples_sigma2_w(n_keep);
  NumericVector samples_eta(n_keep);
  NumericVector samples_rho(n_keep);
  NumericVector samples_loglike(save_loglike ? n_keep : 0);
  int save_index = 0;

  for (int iteration = 1; iteration <= n_iter; ++iteration) {
    if (iteration % 250 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Intercept-only move. This changes every Poisson mean but not the field prior.
    {
      ++beta_attempt;
      ++beta_window_attempt;
      const double proposal = beta0 + proposal_sd_beta0 * R::rnorm(0.0, 1.0);
      std::vector<double> lambda_proposal(n, 0.0);
      std::vector<double> clipped_proposal(n, 0.0);
      const double loglikelihood_proposal = compute_loglikelihood(
        proposal,
        field,
        field_mean,
        log_exposure,
        y,
        log_lambda_lower,
        log_lambda_upper,
        lambda_proposal,
        clipped_proposal
      );
      const double prior_difference = -0.5 * (
        std::pow((proposal - beta0_prior_mean) / beta0_prior_sd, 2.0) -
        std::pow((beta0 - beta0_prior_mean) / beta0_prior_sd, 2.0)
      );
      if (accept_log_mh(loglikelihood_proposal - loglikelihood + prior_difference)) {
        beta0 = proposal;
        loglikelihood = loglikelihood_proposal;
        lambda.swap(lambda_proposal);
        clipped_log_rate.swap(clipped_proposal);
        ++beta_accept;
        ++beta_window_accept;
      }
    }

    // Joint field/intercept moves leave all but one centered linear predictor unchanged.
    for (int i = 0; i < n; ++i) {
      ++field_attempt;
      ++field_window_attempt;
      const double delta = proposal_sd_field * R::rnorm(0.0, 1.0);
      const double beta0_proposal = beta0 + delta / static_cast<double>(n);
      const double field_mean_proposal = field_mean + delta / static_cast<double>(n);
      const double field_proposal_i = field[i] + delta;

      double lambda_proposal_i = 0.0;
      double clipped_proposal_i = 0.0;
      const double proposed_component = clipped_poisson_log_kernel(
        log_exposure[i] + beta0_proposal + field_proposal_i - field_mean_proposal,
        y[i],
        log_lambda_lower,
        log_lambda_upper,
        lambda_proposal_i,
        clipped_proposal_i
      );
      const double current_component = y[i] * clipped_log_rate[i] - lambda[i];

      double delta_quadratic = dagar.lambda[i] * (
        std::pow(residuals[i] + delta, 2.0) - std::pow(residuals[i], 2.0)
      );
      for (std::size_t index = 0; index < dagar.succs[i].size(); ++index) {
        const int successor = dagar.succs[i][index];
        const double proposed_residual = residuals[successor] - dagar.b_row[successor] * delta;
        delta_quadratic += dagar.lambda[successor] * (
          proposed_residual * proposed_residual - residuals[successor] * residuals[successor]
        );
      }
      const double beta_prior_difference = -0.5 * (
        std::pow((beta0_proposal - beta0_prior_mean) / beta0_prior_sd, 2.0) -
        std::pow((beta0 - beta0_prior_mean) / beta0_prior_sd, 2.0)
      );
      const double log_ratio = proposed_component - current_component
        - 0.5 * delta_quadratic / sigma2_w + beta_prior_difference;

      if (accept_log_mh(log_ratio)) {
        field[i] = field_proposal_i;
        field_mean = field_mean_proposal;
        beta0 = beta0_proposal;
        residuals[i] += delta;
        for (std::size_t index = 0; index < dagar.succs[i].size(); ++index) {
          const int successor = dagar.succs[i][index];
          residuals[successor] -= dagar.b_row[successor] * delta;
        }
        quadratic += delta_quadratic;
        loglikelihood += proposed_component - current_component;
        lambda[i] = lambda_proposal_i;
        clipped_log_rate[i] = clipped_proposal_i;
        ++field_accept;
        ++field_window_accept;
      }
    }

    // Log-scale Metropolis update for the half-normal prior on sigma2_w itself.
    {
      ++sigma_attempt;
      ++sigma_window_attempt;
      const double log_proposal = log_sigma2_w + proposal_sd_log_sigma2 * R::rnorm(0.0, 1.0);
      const double proposal = std::exp(log_proposal);
      const double current_target =
        -0.5 * n * log_sigma2_w - 0.5 * quadratic / sigma2_w +
        log_half_normal_kernel(sigma2_w, sigma2_w_prior_sd) + log_sigma2_w;
      const double proposed_target =
        -0.5 * n * log_proposal - 0.5 * quadratic / proposal +
        log_half_normal_kernel(proposal, sigma2_w_prior_sd) + log_proposal;
      if (accept_log_mh(proposed_target - current_target)) {
        sigma2_w = proposal;
        log_sigma2_w = log_proposal;
        ++sigma_accept;
        ++sigma_window_accept;
      }
    }

    // Uniform eta prior; only the repaired graph-dependent DAGAR density changes.
    {
      ++eta_attempt;
      ++eta_window_attempt;
      const double proposal = rtruncnorm_one(eta, proposal_sd_eta, 0.0, eta_max);
      std::vector<unsigned char> graph_proposal;
      build_filtered_graph_with_repair(
        n, Z_r, original_graph, undirected_edges, proposal, threshold, graph_proposal
      );
      const DagarStructure dagar_proposal = build_dagar_structure(graph_proposal, ordering, rho);
      std::vector<double> residuals_proposal;
      const double quadratic_proposal = compute_quadratic_form(field, dagar_proposal, residuals_proposal);
      const double density_difference =
        0.5 * (dagar_proposal.sum_log_lambda - dagar.sum_log_lambda) -
        0.5 * (quadratic_proposal - quadratic) / sigma2_w;
      const double proposal_correction =
        truncnorm_log_density(eta, proposal, proposal_sd_eta, 0.0, eta_max) -
        truncnorm_log_density(proposal, eta, proposal_sd_eta, 0.0, eta_max);
      if (accept_log_mh(density_difference + proposal_correction)) {
        eta = proposal;
        filtered_graph.swap(graph_proposal);
        dagar = dagar_proposal;
        residuals.swap(residuals_proposal);
        quadratic = quadratic_proposal;
        ++eta_accept;
        ++eta_window_accept;
      }
    }

    // Logit-scale proposal with the Jacobian required by Uniform(0, 1) on rho.
    {
      ++rho_attempt;
      ++rho_window_attempt;
      const double logit_proposal = logit_rho + proposal_sd_logit_rho * R::rnorm(0.0, 1.0);
      const double proposal = clip_open_interval(inv_logit(logit_proposal), 0.0, 1.0);
      const DagarStructure dagar_proposal = build_dagar_structure(filtered_graph, ordering, proposal);
      std::vector<double> residuals_proposal;
      const double quadratic_proposal = compute_quadratic_form(field, dagar_proposal, residuals_proposal);
      const double density_difference =
        0.5 * (dagar_proposal.sum_log_lambda - dagar.sum_log_lambda) -
        0.5 * (quadratic_proposal - quadratic) / sigma2_w;
      const double jacobian_difference =
        std::log(proposal) + std::log(1.0 - proposal) -
        std::log(rho) - std::log(1.0 - rho);
      if (accept_log_mh(density_difference + jacobian_difference)) {
        rho = proposal;
        logit_rho = logit_proposal;
        dagar = dagar_proposal;
        residuals.swap(residuals_proposal);
        quadratic = quadratic_proposal;
        ++rho_accept;
        ++rho_window_accept;
      }
    }

    if (iteration % 100 == 0 && iteration <= n_adapt && iteration < burnin) {
      proposal_sd_beta0 = tune_sd(
        proposal_sd_beta0, beta_window_accept, beta_window_attempt, 0.30, 0.40
      );
      proposal_sd_field = tune_sd(
        proposal_sd_field, field_window_accept, field_window_attempt, 0.40, 0.50
      );
      proposal_sd_log_sigma2 = tune_sd(
        proposal_sd_log_sigma2, sigma_window_accept, sigma_window_attempt, 0.30, 0.40
      );
      proposal_sd_eta = tune_sd(
        proposal_sd_eta, eta_window_accept, eta_window_attempt, 0.40, 0.50, eta_max / 4.0
      );
      proposal_sd_logit_rho = tune_sd(
        proposal_sd_logit_rho, rho_window_accept, rho_window_attempt, 0.30, 0.40
      );
      beta_window_accept = beta_window_attempt = 0;
      field_window_accept = field_window_attempt = 0;
      sigma_window_accept = sigma_window_attempt = 0;
      eta_window_accept = eta_window_attempt = 0;
      rho_window_accept = rho_window_attempt = 0;
    }

    if (verbose && iteration % 1000 == 0) {
      Rcout << "Iteration " << iteration << " / " << n_iter
            << " | beta0=" << beta0
            << ", sigma2_w=" << sigma2_w
            << ", eta=" << eta
            << ", rho=" << rho << "\n";
    }

    if (iteration > burnin && (iteration - burnin) % thin == 0) {
      samples_beta0[save_index] = beta0;
      samples_sigma2_w[save_index] = sigma2_w;
      samples_eta[save_index] = eta;
      samples_rho[save_index] = rho;
      if (save_loglike) {
        samples_loglike[save_index] = loglikelihood;
      }
      ++save_index;
    }
  }

  IntegerVector ordering_output(n);
  for (int i = 0; i < n; ++i) {
    ordering_output[i] = ordering[i] + 1;
  }

  NumericVector acceptance = NumericVector::create(
    _["beta0"] = static_cast<double>(beta_accept) / std::max(beta_attempt, 1),
    _["field"] = static_cast<double>(field_accept) / std::max(field_attempt, 1),
    _["sigma2_w"] = static_cast<double>(sigma_accept) / std::max(sigma_attempt, 1),
    _["eta"] = static_cast<double>(eta_accept) / std::max(eta_attempt, 1),
    _["rho"] = static_cast<double>(rho_accept) / std::max(rho_attempt, 1)
  );

  SEXP loglike_output = R_NilValue;
  if (save_loglike) {
    loglike_output = samples_loglike;
  }
  const double centered_field_mean =
    std::accumulate(field.begin(), field.end(), 0.0) / static_cast<double>(n) - field_mean;

  return List::create(
    _["beta0"] = samples_beta0,
    _["sigma2_w"] = samples_sigma2_w,
    _["eta"] = samples_eta,
    _["rho"] = samples_rho,
    _["loglike"] = loglike_output,
    _["acceptance"] = acceptance,
    _["ordering"] = ordering_output,
    _["A_filtered_last"] = graph_to_matrix(filtered_graph, n),
    _["field_mean_last"] = field_mean,
    _["centered_field_mean_last"] = centered_field_mean,
    _["proposal_sds_final"] = List::create(
      _["beta0"] = proposal_sd_beta0,
      _["field"] = proposal_sd_field,
      _["log_sigma2_w"] = proposal_sd_log_sigma2,
      _["eta"] = proposal_sd_eta,
      _["logit_rho"] = proposal_sd_logit_rho
    )
  );
}
