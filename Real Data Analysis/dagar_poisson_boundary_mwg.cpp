#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

inline int flat_index(const int i, const int j, const int n) {
  return i * n + j;
}

inline bool accept_log_mh(const double log_alpha) {
  return R_FINITE(log_alpha) && (std::log(R::runif(0.0, 1.0)) < log_alpha);
}

inline double inv_logit(const double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

inline double clip_open_interval(
    const double x,
    const double lower,
    const double upper,
    const double eps = 1e-10
) {
  const double lo = lower + eps;
  const double hi = upper - eps;
  if (x <= lo) {
    return lo;
  }
  if (x >= hi) {
    return hi;
  }
  return x;
}

inline double safe_exp_clip(
    const double log_rate,
    const double lambda_lower,
    const double lambda_upper,
    const double log_lambda_lower,
    const double log_lambda_upper,
    double& log_lambda_out
) {
  if (log_rate <= log_lambda_lower) {
    log_lambda_out = log_lambda_lower;
    return lambda_lower;
  }
  if (log_rate >= log_lambda_upper) {
    log_lambda_out = log_lambda_upper;
    return lambda_upper;
  }
  log_lambda_out = log_rate;
  return std::exp(log_rate);
}

inline double truncnorm_log_density(
    const double x,
    const double mean,
    const double sd,
    const double lower,
    const double upper
) {
  if (sd <= 0.0 || x < lower || x > upper) {
    return R_NegInf;
  }

  const double z_lower = (lower - mean) / sd;
  const double z_upper = (upper - mean) / sd;
  const double cdf_lower = R::pnorm5(z_lower, 0.0, 1.0, 1, 0);
  const double cdf_upper = R::pnorm5(z_upper, 0.0, 1.0, 1, 0);
  const double norm_const = cdf_upper - cdf_lower;
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
  const double z_lower = (lower - mean) / sd;
  const double z_upper = (upper - mean) / sd;
  const double cdf_lower = R::pnorm5(z_lower, 0.0, 1.0, 1, 0);
  const double cdf_upper = R::pnorm5(z_upper, 0.0, 1.0, 1, 0);
  const double u = R::runif(cdf_lower, cdf_upper);
  const double z = R::qnorm5(u, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

inline double tune_sd_scalar(
    const double sd,
    const int accept_count,
    const int attempt_count,
    const double min_rate,
    const double max_rate
) {
  if (attempt_count <= 0) {
    return sd;
  }
  const double rate = 100.0 * static_cast<double>(accept_count) /
    static_cast<double>(attempt_count);
  if (rate > max_rate) {
    return sd + 0.1 * sd;
  }
  if (rate < min_rate) {
    return sd - 0.1 * sd;
  }
  return sd;
}

inline double tune_sd_scalar_bounded(
    const double sd,
    const int accept_count,
    const int attempt_count,
    const double min_rate,
    const double max_rate,
    const double sd_max
) {
  double out = tune_sd_scalar(sd, accept_count, attempt_count, min_rate, max_rate);
  if (out > sd_max) {
    out = sd_max;
  }
  return out;
}

struct DagarStructure {
  std::vector<double> lambda;
  std::vector<double> b_row;
  std::vector< std::vector<int> > preds;
  std::vector< std::vector<int> > succs;
  double sum_log_lambda;
};

void validate_square_matrix(const NumericMatrix& mat, const std::string& name) {
  if (mat.nrow() != mat.ncol()) {
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

  const IntegerVector ord_in(ordering_r);
  if (ord_in.size() != n) {
    stop("ordering must have length n.");
  }

  std::vector<int> seen(n, 0);
  for (int pos = 0; pos < n; ++pos) {
    const int idx = ord_in[pos] - 1;
    if (idx < 0 || idx >= n) {
      stop("ordering must be a 1-based permutation of 1:n.");
    }
    if (seen[idx]) {
      stop("ordering contains repeated indices.");
    }
    seen[idx] = 1;
    ordering[pos] = idx;
  }
  return ordering;
}

void build_original_graph_info(
    const NumericMatrix& A,
    std::vector<unsigned char>& original_graph,
    std::vector< std::pair<int, int> >& undirected_edges
) {
  const int n = A.nrow();
  original_graph.assign(n * n, 0);
  undirected_edges.clear();

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const bool linked = (A(i, j) > 0.5) || (A(j, i) > 0.5);
      if (!linked) {
        continue;
      }
      original_graph[flat_index(i, j, n)] = 1;
      original_graph[flat_index(j, i, n)] = 1;
      undirected_edges.push_back(std::make_pair(i, j));
    }
  }
}

void build_filtered_graph_no_repair(
    const int n,
    const NumericMatrix& Z,
    const std::vector< std::pair<int, int> >& undirected_edges,
    const double alpha,
    const double threshold,
    std::vector<unsigned char>& filtered_graph
) {
  filtered_graph.assign(n * n, 0);

  for (std::size_t e = 0; e < undirected_edges.size(); ++e) {
    const int i = undirected_edges[e].first;
    const int j = undirected_edges[e].second;
    if (Z(i, j) * alpha <= threshold) {
      filtered_graph[flat_index(i, j, n)] = 1;
      filtered_graph[flat_index(j, i, n)] = 1;
    }
  }
}

DagarStructure build_dagar_structure(
    const std::vector<unsigned char>& graph,
    const std::vector<int>& ordering,
    const double rho
) {
  const int n = static_cast<int>(ordering.size());
  const double rho2 = rho * rho;

  DagarStructure out;
  out.lambda.assign(n, 0.0);
  out.b_row.assign(n, 0.0);
  out.preds.assign(n, std::vector<int>());
  out.succs.assign(n, std::vector<int>());
  out.sum_log_lambda = 0.0;

  for (int pos = 0; pos < n; ++pos) {
    const int i = ordering[pos];
    std::vector<int>& preds_i = out.preds[i];

    for (int q = 0; q < pos; ++q) {
      const int j = ordering[q];
      if (graph[flat_index(i, j, n)]) {
        preds_i.push_back(j);
      }
    }

    const int n_lt = static_cast<int>(preds_i.size());
    const double denom = 1.0 + std::max(n_lt - 1, 0) * rho2;
    const double b_val = (n_lt > 0) ? (rho / denom) : 0.0;
    const double lambda_i = denom / (1.0 - rho2);

    out.b_row[i] = b_val;
    out.lambda[i] = lambda_i;
    out.sum_log_lambda += std::log(lambda_i);

    for (std::size_t idx = 0; idx < preds_i.size(); ++idx) {
      out.succs[preds_i[idx]].push_back(i);
    }
  }

  return out;
}

double compute_quadratic_form(
    const std::vector<double>& phi,
    const DagarStructure& dagar,
    std::vector<double>& residuals
) {
  const int n = static_cast<int>(phi.size());
  residuals.assign(n, 0.0);
  double Q = 0.0;

  for (int i = 0; i < n; ++i) {
    double pred_sum = 0.0;
    const std::vector<int>& preds_i = dagar.preds[i];
    for (std::size_t idx = 0; idx < preds_i.size(); ++idx) {
      pred_sum += phi[preds_i[idx]];
    }
    residuals[i] = phi[i] - dagar.b_row[i] * pred_sum;
    Q += dagar.lambda[i] * residuals[i] * residuals[i];
  }

  return Q;
}

double compute_loglik_full(
    const double beta0,
    const std::vector<double>& phi,
    const std::vector<double>& log_e,
    const std::vector<double>& y,
    const double lambda_lower,
    const double lambda_upper,
    const double log_lambda_lower,
    const double log_lambda_upper,
    std::vector<double>& lambda_out,
    std::vector<double>& log_lambda_out
) {
  const int n = static_cast<int>(phi.size());
  double loglik = 0.0;

  for (int i = 0; i < n; ++i) {
    const double log_rate = log_e[i] + beta0 + phi[i];
    lambda_out[i] = safe_exp_clip(
      log_rate,
      lambda_lower,
      lambda_upper,
      log_lambda_lower,
      log_lambda_upper,
      log_lambda_out[i]
    );
    loglik += y[i] * log_lambda_out[i] - lambda_out[i];
  }

  return loglik;
}

inline double compute_loglik_single_move(
    const int i,
    const double phi_prop_i,
    const double beta0,
    const std::vector<double>& phi,
    const std::vector<double>& log_e,
    const std::vector<double>& y,
    const double lambda_lower,
    const double lambda_upper,
    const double log_lambda_lower,
    const double log_lambda_upper,
    double& lambda_prop_i,
    double& log_lambda_prop_i
) {
  const double log_rate_curr = log_e[i] + beta0 + phi[i];
  double log_lambda_curr_i;
  const double lambda_curr_i = safe_exp_clip(
    log_rate_curr,
    lambda_lower,
    lambda_upper,
    log_lambda_lower,
    log_lambda_upper,
    log_lambda_curr_i
  );

  const double log_rate_prop = log_e[i] + beta0 + phi_prop_i;
  lambda_prop_i = safe_exp_clip(
    log_rate_prop,
    lambda_lower,
    lambda_upper,
    log_lambda_lower,
    log_lambda_upper,
    log_lambda_prop_i
  );

  return y[i] * (log_lambda_prop_i - log_lambda_curr_i) -
    (lambda_prop_i - lambda_curr_i);
}

NumericMatrix graph_to_matrix(const std::vector<unsigned char>& graph, const int n) {
  NumericMatrix out(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out(i, j) = static_cast<double>(graph[flat_index(i, j, n)]);
    }
  }
  return out;
}

void compute_boundary_summaries(
    const NumericMatrix& Z,
    const std::vector< std::pair<int, int> >& undirected_edges,
    const NumericVector& alpha_draws,
    const double threshold,
    NumericMatrix& W_posterior,
    NumericMatrix& W_border_prob
) {
  const int n = Z.nrow();
  const int n_draws = alpha_draws.size();

  W_posterior = NumericMatrix(n, n);
  W_border_prob = NumericMatrix(n, n);
  std::fill(W_posterior.begin(), W_posterior.end(), NA_REAL);
  std::fill(W_border_prob.begin(), W_border_prob.end(), NA_REAL);

  for (std::size_t e = 0; e < undirected_edges.size(); ++e) {
    const int i = undirected_edges[e].first;
    const int j = undirected_edges[e].second;
    const double zij = Z(i, j);

    int n_smooth = 0;
    for (int s = 0; s < n_draws; ++s) {
      if (zij * alpha_draws[s] <= threshold) {
        ++n_smooth;
      }
    }

    const double border_prob = 1.0 - static_cast<double>(n_smooth) /
      static_cast<double>(n_draws);
    const double posterior_keep = (2 * n_smooth >= n_draws) ? 1.0 : 0.0;

    W_border_prob(i, j) = border_prob;
    W_border_prob(j, i) = border_prob;
    W_posterior(i, j) = posterior_keep;
    W_posterior(j, i) = posterior_keep;
  }
}

} // namespace

// [[Rcpp::export]]
Rcpp::List dagar_poisson_boundary_mwg2(
    const NumericVector& y_r,
    const NumericVector& e_r,
    const NumericMatrix& A_r,
    const NumericMatrix& Z_r,
    const int n_iter,
    const int burnin = 5000,
    const int thin = 1,
    const int n_adapt = 5000,
    Nullable<IntegerVector> ordering_r = R_NilValue,
    Nullable<NumericVector> phi_init_r = R_NilValue,
    const double beta0_init = NA_REAL,
    const double tau2_init = 0.25,
    const double alpha_init = NA_REAL,
    const double alpha_max = NA_REAL,
    const double rho_init = 0.5,
    const double proposal_sd_beta0_init = 0.01,
    const double proposal_sd_phi_init = 0.10,
    const double proposal_sd_alpha_init = NA_REAL,
    const double proposal_sd_logit_rho_init = 0.10,
    const double beta0_prior_mean = 0.0,
    const double beta0_prior_var = 1e5,
    const double tau2_prior_shape = 1.0,
    const double tau2_prior_scale = 0.01,
    const double threshold = 0.6931471805599453,
    const double lambda_lower = 1e-2,
    const double lambda_upper = 1e6,
    const bool save_phi = true,
    const bool save_fitted = true,
    const bool save_loglike = true,
    const bool verbose = false
) {
  if (n_iter <= 0) {
    stop("n_iter must be positive.");
  }
  if (burnin < 0) {
    stop("burnin must be non-negative.");
  }
  if (thin <= 0) {
    stop("thin must be positive.");
  }
  if (n_adapt < 0) {
    stop("n_adapt must be non-negative.");
  }
  if (lambda_lower <= 0.0 || lambda_upper <= lambda_lower) {
    stop("lambda bounds must satisfy 0 < lambda_lower < lambda_upper.");
  }
  if (beta0_prior_var <= 0.0) {
    stop("beta0_prior_var must be positive.");
  }
  if (tau2_prior_shape <= 0.0 || tau2_prior_scale <= 0.0) {
    stop("tau2 prior shape and scale must be positive.");
  }
  if (rho_init <= 0.0 || rho_init >= 1.0) {
    stop("rho_init must lie strictly between 0 and 1.");
  }

  validate_square_matrix(A_r, "A");
  validate_square_matrix(Z_r, "Z");

  const int n = y_r.size();
  if (e_r.size() != n) {
    stop("y and e must have the same length.");
  }
  if (A_r.nrow() != n || Z_r.nrow() != n) {
    stop("A and Z must have dimensions n x n, where n = length(y).");
  }

  std::vector<double> y(n);
  std::vector<double> log_e(n);
  double sum_y = 0.0;
  double sum_e = 0.0;
  for (int i = 0; i < n; ++i) {
    if (y_r[i] < 0.0) {
      stop("y must contain non-negative counts.");
    }
    if (e_r[i] <= 0.0) {
      stop("e must contain strictly positive exposures.");
    }
    y[i] = y_r[i];
    log_e[i] = std::log(e_r[i]);
    sum_y += y[i];
    sum_e += e_r[i];
  }

  double alpha_max_use = alpha_max;
  if (NumericVector::is_na(alpha_max_use)) {
    std::vector<double> z_nonzero;
    z_nonzero.reserve(n * n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double zij = Z_r(i, j);
        if (zij > 0.0) {
          z_nonzero.push_back(zij);
        }
      }
    }
    if (z_nonzero.empty()) {
      stop("Could not compute alpha_max because Z has no positive off-diagonal values.");
    }
    std::sort(z_nonzero.begin(), z_nonzero.end());
    const int m = static_cast<int>(z_nonzero.size());
    const double z_crit = (m % 2 == 1)
      ? z_nonzero[m / 2]
      : 0.5 * (z_nonzero[m / 2 - 1] + z_nonzero[m / 2]);
    alpha_max_use = -std::log(0.5) / z_crit;
  }
  if (alpha_max_use <= 0.0) {
    stop("alpha_max must be positive.");
  }

  double max_z = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (Z_r(i, j) > max_z) {
        max_z = Z_r(i, j);
      }
    }
  }
  const double alpha_threshold = (max_z > 0.0) ? (-std::log(0.5) / max_z) : NA_REAL;

  std::vector<int> ordering = parse_ordering(n, ordering_r);

  std::vector<unsigned char> original_graph;
  std::vector< std::pair<int, int> > undirected_edges;
  build_original_graph_info(A_r, original_graph, undirected_edges);

  double beta0 = NumericVector::is_na(beta0_init)
    ? std::log((sum_y + 0.5) / (sum_e + 0.5))
    : beta0_init;

  double tau2 = (tau2_init > 0.0) ? tau2_init : 0.25;
  double alpha = NumericVector::is_na(alpha_init)
    ? R::runif(0.0, alpha_max_use / 3.0)
    : clip_open_interval(alpha_init, 0.0, alpha_max_use);
  double rho = rho_init;
  double logit_rho = std::log(rho) - std::log(1.0 - rho);

  std::vector<double> phi(n, 0.0);
  if (!phi_init_r.isNull()) {
    const NumericVector phi_in(phi_init_r);
    if (phi_in.size() != n) {
      stop("phi_init must have length n.");
    }
    for (int i = 0; i < n; ++i) {
      phi[i] = phi_in[i];
    }
  }

  std::vector<unsigned char> current_graph;
  build_filtered_graph_no_repair(n, Z_r, undirected_edges, alpha, threshold, current_graph);

  DagarStructure dagar = build_dagar_structure(current_graph, ordering, rho);
  std::vector<double> residuals;
  double Q = compute_quadratic_form(phi, dagar, residuals);

  const double log_lambda_lower = std::log(lambda_lower);
  const double log_lambda_upper = std::log(lambda_upper);

  std::vector<double> lambda_curr(n, 0.0);
  std::vector<double> log_lambda_curr(n, 0.0);
  double loglik_curr = compute_loglik_full(
    beta0,
    phi,
    log_e,
    y,
    lambda_lower,
    lambda_upper,
    log_lambda_lower,
    log_lambda_upper,
    lambda_curr,
    log_lambda_curr
  );

  double proposal_sd_beta0 = std::max(proposal_sd_beta0_init, 1e-8);
  double proposal_sd_phi = std::max(proposal_sd_phi_init, 1e-8);
  double proposal_sd_alpha = NumericVector::is_na(proposal_sd_alpha_init)
    ? (0.02 * alpha_max_use)
    : std::max(proposal_sd_alpha_init, 1e-8);
  double proposal_sd_logit_rho = std::max(proposal_sd_logit_rho_init, 1e-8);

  int beta_accept_total = 0;
  int beta_attempt_total = 0;
  int phi_accept_total = 0;
  int phi_attempt_total = 0;
  int alpha_accept_total = 0;
  int alpha_attempt_total = 0;
  int rho_accept_total = 0;
  int rho_attempt_total = 0;

  int beta_accept_window = 0;
  int beta_attempt_window = 0;
  int phi_accept_window = 0;
  int phi_attempt_window = 0;
  int alpha_accept_window = 0;
  int alpha_attempt_window = 0;
  int rho_accept_window = 0;
  int rho_attempt_window = 0;

  const int n_keep = (n_iter > burnin) ? ((n_iter - burnin + thin - 1) / thin) : 0;
  NumericMatrix samples_beta(n_keep, 1);
  NumericMatrix samples_tau2(n_keep, 1);
  NumericMatrix samples_alpha(n_keep, 1);
  NumericMatrix samples_rho(n_keep, 1);
  NumericMatrix samples_phi;
  NumericMatrix samples_fitted;
  NumericMatrix samples_loglike;
  if (save_phi) {
    samples_phi = NumericMatrix(n_keep, n);
  }
  if (save_fitted) {
    samples_fitted = NumericMatrix(n_keep, n);
  }
  if (save_loglike) {
    samples_loglike = NumericMatrix(n_keep, n);
  }

  int save_idx = 0;

  for (int iter = 1; iter <= n_iter; ++iter) {
    if ((iter % 250) == 0) {
      Rcpp::checkUserInterrupt();
    }

    // beta0 update
    {
      ++beta_attempt_total;
      ++beta_attempt_window;

      const double beta0_prop = beta0 + proposal_sd_beta0 * R::rnorm(0.0, 1.0);
      std::vector<double> lambda_prop(n, 0.0);
      std::vector<double> log_lambda_prop(n, 0.0);
      const double loglik_prop = compute_loglik_full(
        beta0_prop,
        phi,
        log_e,
        y,
        lambda_lower,
        lambda_upper,
        log_lambda_lower,
        log_lambda_upper,
        lambda_prop,
        log_lambda_prop
      );

      const double log_prior_diff =
        -0.5 * (
          (beta0_prop - beta0_prior_mean) * (beta0_prop - beta0_prior_mean) -
          (beta0 - beta0_prior_mean) * (beta0 - beta0_prior_mean)
        ) / beta0_prior_var;

      const double log_alpha_mh = (loglik_prop - loglik_curr) + log_prior_diff;
      if (accept_log_mh(log_alpha_mh)) {
        beta0 = beta0_prop;
        loglik_curr = loglik_prop;
        lambda_curr.swap(lambda_prop);
        log_lambda_curr.swap(log_lambda_prop);
        ++beta_accept_total;
        ++beta_accept_window;
      }
    }

    // phi sweep
    for (int i = 0; i < n; ++i) {
      ++phi_attempt_total;
      ++phi_attempt_window;

      const double delta = proposal_sd_phi * R::rnorm(0.0, 1.0);
      const double phi_prop_i = phi[i] + delta;

      double lambda_prop_i;
      double log_lambda_prop_i;
      const double loglik_diff = compute_loglik_single_move(
        i,
        phi_prop_i,
        beta0,
        phi,
        log_e,
        y,
        lambda_lower,
        lambda_upper,
        log_lambda_lower,
        log_lambda_upper,
        lambda_prop_i,
        log_lambda_prop_i
      );

      double delta_Q = dagar.lambda[i] *
        ((residuals[i] + delta) * (residuals[i] + delta) - residuals[i] * residuals[i]);

      const std::vector<int>& succs_i = dagar.succs[i];
      for (std::size_t idx = 0; idx < succs_i.size(); ++idx) {
        const int k = succs_i[idx];
        const double rk_prop = residuals[k] - dagar.b_row[k] * delta;
        delta_Q += dagar.lambda[k] * (rk_prop * rk_prop - residuals[k] * residuals[k]);
      }

      const double log_prior_diff = -0.5 * delta_Q / tau2;
      const double log_alpha_mh = loglik_diff + log_prior_diff;

      if (accept_log_mh(log_alpha_mh)) {
        phi[i] = phi_prop_i;
        residuals[i] += delta;
        for (std::size_t idx = 0; idx < succs_i.size(); ++idx) {
          const int k = succs_i[idx];
          residuals[k] -= dagar.b_row[k] * delta;
        }
        Q += delta_Q;
        loglik_curr += loglik_diff;
        lambda_curr[i] = lambda_prop_i;
        log_lambda_curr[i] = log_lambda_prop_i;
        ++phi_accept_total;
        ++phi_accept_window;
      }
    }

    // match CARBayes identifiability step
    {
      const double mean_phi = std::accumulate(phi.begin(), phi.end(), 0.0) /
        static_cast<double>(n);
      for (int i = 0; i < n; ++i) {
        phi[i] -= mean_phi;
      }
      Q = compute_quadratic_form(phi, dagar, residuals);
      loglik_curr = compute_loglik_full(
        beta0,
        phi,
        log_e,
        y,
        lambda_lower,
        lambda_upper,
        log_lambda_lower,
        log_lambda_upper,
        lambda_curr,
        log_lambda_curr
      );
    }

    // tau2 Gibbs update under IG prior
    {
      const double shape_post = tau2_prior_shape + 0.5 * static_cast<double>(n);
      const double scale_post = tau2_prior_scale + 0.5 * Q;
      tau2 = 1.0 / R::rgamma(shape_post, 1.0 / scale_post);
    }

    // alpha update using CARBayes-style truncated-normal proposal
    {
      ++alpha_attempt_total;
      ++alpha_attempt_window;

      const double alpha_prop = rtruncnorm_one(alpha, proposal_sd_alpha, 0.0, alpha_max_use);
      std::vector<unsigned char> graph_prop;
      build_filtered_graph_no_repair(n, Z_r, undirected_edges, alpha_prop, threshold, graph_prop);

      const DagarStructure dagar_prop = build_dagar_structure(graph_prop, ordering, rho);
      std::vector<double> residuals_prop;
      const double Q_prop = compute_quadratic_form(phi, dagar_prop, residuals_prop);

      const double logprob_curr = 0.5 * dagar.sum_log_lambda - 0.5 * Q / tau2;
      const double logprob_prop = 0.5 * dagar_prop.sum_log_lambda - 0.5 * Q_prop / tau2;
      const double log_hastings =
        truncnorm_log_density(alpha, alpha_prop, proposal_sd_alpha, 0.0, alpha_max_use) -
        truncnorm_log_density(alpha_prop, alpha, proposal_sd_alpha, 0.0, alpha_max_use);

      const double log_alpha_mh = (logprob_prop - logprob_curr) + log_hastings;
      if (accept_log_mh(log_alpha_mh)) {
        alpha = alpha_prop;
        current_graph.swap(graph_prop);
        dagar = dagar_prop;
        residuals.swap(residuals_prop);
        Q = Q_prop;
        ++alpha_accept_total;
        ++alpha_accept_window;
      }
    }

    // rho update on logit scale with ABI-style Uniform(0,1) prior
    {
      ++rho_attempt_total;
      ++rho_attempt_window;

      const double logit_rho_prop = logit_rho + proposal_sd_logit_rho * R::rnorm(0.0, 1.0);
      const double rho_prop = clip_open_interval(inv_logit(logit_rho_prop), 0.0, 1.0);

      const DagarStructure dagar_prop = build_dagar_structure(current_graph, ordering, rho_prop);
      std::vector<double> residuals_prop;
      const double Q_prop = compute_quadratic_form(phi, dagar_prop, residuals_prop);

      const double logprob_curr = 0.5 * dagar.sum_log_lambda - 0.5 * Q / tau2;
      const double logprob_prop = 0.5 * dagar_prop.sum_log_lambda - 0.5 * Q_prop / tau2;
      const double log_jacobian =
        (std::log(rho_prop) + std::log(1.0 - rho_prop)) -
        (std::log(rho) + std::log(1.0 - rho));

      const double log_alpha_mh = (logprob_prop - logprob_curr) + log_jacobian;
      if (accept_log_mh(log_alpha_mh)) {
        rho = rho_prop;
        logit_rho = logit_rho_prop;
        dagar = dagar_prop;
        residuals.swap(residuals_prop);
        Q = Q_prop;
        ++rho_accept_total;
        ++rho_accept_window;
      }
    }

    if ((iter % 100) == 0 && iter < burnin) {
      proposal_sd_beta0 = tune_sd_scalar(
        proposal_sd_beta0,
        beta_accept_window,
        beta_attempt_window,
        30.0,
        40.0
      );
      proposal_sd_phi = tune_sd_scalar(
        proposal_sd_phi,
        phi_accept_window,
        phi_attempt_window,
        40.0,
        50.0
      );
      proposal_sd_alpha = tune_sd_scalar_bounded(
        proposal_sd_alpha,
        alpha_accept_window,
        alpha_attempt_window,
        40.0,
        50.0,
        alpha_max_use / 4.0
      );
      proposal_sd_logit_rho = tune_sd_scalar(
        proposal_sd_logit_rho,
        rho_accept_window,
        rho_attempt_window,
        40.0,
        50.0
      );

      beta_accept_window = beta_attempt_window = 0;
      phi_accept_window = phi_attempt_window = 0;
      alpha_accept_window = alpha_attempt_window = 0;
      rho_accept_window = rho_attempt_window = 0;
    }

    if (verbose && (iter % 1000 == 0)) {
      Rcout << "Iteration " << iter
            << " / " << n_iter
            << " | beta0=" << beta0
            << ", tau2=" << tau2
            << ", alpha=" << alpha
            << ", rho=" << rho
            << "\n";
    }

    if (iter > burnin && ((iter - burnin) % thin == 0)) {
      samples_beta(save_idx, 0) = beta0;
      samples_tau2(save_idx, 0) = tau2;
      samples_alpha(save_idx, 0) = alpha;
      samples_rho(save_idx, 0) = rho;

      if (save_phi) {
        for (int i = 0; i < n; ++i) {
          samples_phi(save_idx, i) = phi[i];
        }
      }
      if (save_fitted) {
        for (int i = 0; i < n; ++i) {
          samples_fitted(save_idx, i) = lambda_curr[i];
        }
      }
      if (save_loglike) {
        for (int i = 0; i < n; ++i) {
          samples_loglike(save_idx, i) = y[i] * log_lambda_curr[i] - lambda_curr[i];
        }
      }
      ++save_idx;
    }
  }

  NumericVector fitted_values(n, NA_REAL);
  if (save_fitted && n_keep > 0) {
    for (int i = 0; i < n; ++i) {
      double accum = 0.0;
      for (int s = 0; s < n_keep; ++s) {
        accum += samples_fitted(s, i);
      }
      fitted_values[i] = accum / static_cast<double>(n_keep);
    }
  }

  NumericMatrix W_posterior;
  NumericMatrix W_border_prob;
  if (n_keep > 0) {
    NumericVector alpha_draws_only(n_keep);
    for (int s = 0; s < n_keep; ++s) {
      alpha_draws_only[s] = samples_alpha(s, 0);
    }
    compute_boundary_summaries(
      Z_r,
      undirected_edges,
      alpha_draws_only,
      threshold,
      W_posterior,
      W_border_prob
    );
  }

  IntegerVector ordering_out(n);
  for (int i = 0; i < n; ++i) {
    ordering_out[i] = ordering[i] + 1;
  }

  SEXP samples_phi_out = R_NilValue;
  SEXP samples_fitted_out = R_NilValue;
  SEXP samples_loglike_out = R_NilValue;
  if (save_phi) {
    samples_phi_out = samples_phi;
  }
  if (save_fitted) {
    samples_fitted_out = samples_fitted;
  }
  if (save_loglike) {
    samples_loglike_out = samples_loglike;
  }

  NumericVector beta0_draws(n_keep);
  NumericVector tau2_draws(n_keep);
  NumericVector alpha_draws(n_keep);
  NumericVector rho_draws(n_keep);
  for (int s = 0; s < n_keep; ++s) {
    beta0_draws[s] = samples_beta(s, 0);
    tau2_draws[s] = samples_tau2(s, 0);
    alpha_draws[s] = samples_alpha(s, 0);
    rho_draws[s] = samples_rho(s, 0);
  }

  List out = List::create(
    _["samples_beta"] = samples_beta,
    _["samples_tau2"] = samples_tau2,
    _["samples_alpha"] = samples_alpha,
    _["samples_rho"] = samples_rho,
    _["samples_phi"] = samples_phi_out,
    _["samples_fitted"] = samples_fitted_out,
    _["samples_loglike"] = samples_loglike_out,
    _["fitted_values"] = fitted_values,
    _["localised_structure"] = List::create(
      _["W.posterior"] = W_posterior,
      _["W.border.prob"] = W_border_prob
    ),
    _["accept"] = NumericVector::create(
      static_cast<double>(beta_accept_total),
      static_cast<double>(beta_attempt_total),
      static_cast<double>(phi_accept_total),
      static_cast<double>(phi_attempt_total),
      static_cast<double>(alpha_accept_total),
      static_cast<double>(alpha_attempt_total),
      static_cast<double>(rho_accept_total),
      static_cast<double>(rho_attempt_total)
    ),
    _["beta0"] = beta0_draws,
    _["sigma2_w"] = tau2_draws,
    _["eta"] = alpha_draws,
    _["alpha"] = alpha_draws,
    _["rho"] = rho_draws,
    _["alpha_max"] = alpha_max_use,
    _["alpha_threshold"] = alpha_threshold,
    _["A_filtered_last"] = graph_to_matrix(current_graph, n),
    _["ordering"] = ordering_out,
    _["proposal_sds_final"] = List::create(
      _["beta0"] = proposal_sd_beta0,
      _["phi"] = proposal_sd_phi,
      _["alpha"] = proposal_sd_alpha,
      _["logit_rho"] = proposal_sd_logit_rho
    )
  );

  return out;
}
