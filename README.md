# Amortized Bayesian Inference for Spatial Boundary Detection

This repository contains the code and selected outputs for an amortized
Bayesian inference workflow for Poisson areal count data with spatial boundary
detection. The main implementation trains a reusable neural posterior
approximator for a localized DAGAR model and evaluates it through simulation
experiments, model-matched MCMC benchmarks, ablation studies, and real-data
applications to the Glasgow respiratory disease data and the California lung
cancer data.

The short version: if you want to understand what produced a result in the
paper, start from the map below. The repository is organized around the same
story as the manuscript:

1. Train one amortized posterior approximator.
2. Validate it on simulated areal graphs.
3. Benchmark it against model-matched MCMC.
4. Check which engineered summaries matter through ablations.
5. Deploy the trained approximator on real data and compare with CARBayes and
   DAGARBayes reference analyses.

## Project at a Glance

The statistical target is posterior inference for

```text
(beta_0, sigma_w^2, eta, rho)
```

where `beta_0` is the global log-risk level, `sigma_w^2` is the latent spatial
variance, `eta` controls covariate-driven boundary formation, and `rho` controls
residual DAGAR spatial dependence. Boundary probabilities are derived from the
posterior draws by applying the edge-deletion mechanism to each neighboring
pair.

The main trained model is the ABI-DAGAR network stored at:

```text
Training/Checkpoints/poisson_dagar.keras
```

The core training notebook is:

```text
Training/ABI_poisson_regression_DAGAR.ipynb
```

## Repository Layout

```text
ABI_poisson_regression/
|-- Training/
|   |-- ABI_poisson_regression_DAGAR.ipynb
|   |-- training_history.csv
|   |-- Checkpoints/
|   |   `-- poisson_dagar.keras
|   `-- Images/
|       `-- history_plot.png
|
|-- Simulation Experiments/
|   |-- simulation_more_detailed.ipynb
|   |-- Images/
|   |-- ABI_vs_MCMC/
|   `-- Ablation_experiments/
|
|-- Real Data Analysis/
|   |-- Data/
|   |-- CARBayes_california_glasgow.R
|   |-- DAGARBayes_california_glasgow.R
|   |-- dagar_poisson_boundary_mwg.cpp
|   |-- CARBayes_vs_ABI_real_data.ipynb
|   |-- DAGARBayes_vs_ABI_real_data.ipynb
|   |-- setup_and_diagnostics/
|   |-- results_ABI_vs_CARBayes/
|   `-- results_ABI_vs_DAGARBayes/
|
`-- .gitignore
```

Some exploratory files and heavier intermediate artifacts are intentionally
ignored. The public-facing tree keeps the files needed to reproduce the
manuscript results and the selected outputs used in the text.

## Main Workflow

### 1. Train the ABI-DAGAR posterior approximator

Entry point:

```text
Training/ABI_poisson_regression_DAGAR.ipynb
```

This notebook defines the simulator, summary construction, BayesFlow workflow,
summary network, and inference network. The trained approximator uses a
SetTransformer summary network and a spline CouplingFlow inference network.

Primary outputs:

```text
Training/Checkpoints/poisson_dagar.keras
Training/training_history.csv
Training/Images/history_plot.png
```

The manuscript reports the one-time training stage and then treats the saved
checkpoint as the reusable posterior approximator for all downstream analyses.

### 2. Run the main simulation validation

Entry point:

```text
Simulation Experiments/simulation_more_detailed.ipynb
```

This notebook evaluates the trained ABI-DAGAR network on held-out simulated
datasets with varying graph sizes and topologies. It produces parameter
recovery, posterior calibration, boundary-probability diagnostics, decision-rule
diagnostics, posterior predictive checks, and runtime summaries.

Primary figure outputs:

```text
Simulation Experiments/Images/history_plot.png
Simulation Experiments/Images/poisson_dagar_recovery.png
Simulation Experiments/Images/poisson_dagar_calibration_histogram.png
Simulation Experiments/Images/poisson_dagar_calibration_ecdf.png
Simulation Experiments/Images/parameter_coverage_curve.png
Simulation Experiments/Images/boundary_sens_spec_curve.png
Simulation Experiments/Images/boundary_mpm_histograms.png
Simulation Experiments/Images/boundary_probability_quality.png
Simulation Experiments/Images/predictive_checks.png
Simulation Experiments/Images/graph_size_posterior_zscores.png
Simulation Experiments/Images/regime_specific_parameter_errors.png
```

Selected table outputs:

```text
Simulation Experiments/boundary_mpm_results.csv
Simulation Experiments/boundary_mpm_summary.csv
```

### 3. Compare ABI with model-matched MCMC

Folder:

```text
Simulation Experiments/ABI_vs_MCMC/
```

This folder contains the model-matched simulation benchmark used to compare the
trained ABI-DAGAR posterior approximation with a DAGAR MCMC sampler on the same
synthetic datasets.

Key files:

```text
Simulation Experiments/ABI_vs_MCMC/export_benchmark_datasets.py
Simulation Experiments/ABI_vs_MCMC/benchmark_utils.py
Simulation Experiments/ABI_vs_MCMC/dagar_poisson_boundary_mwg.cpp
Simulation Experiments/ABI_vs_MCMC/run_abi_benchmark.ipynb
Simulation Experiments/ABI_vs_MCMC/run_mcmc_benchmark.R
Simulation Experiments/ABI_vs_MCMC/compare_abi_vs_mcmc.ipynb
```

Recommended order if regenerating from scratch:

1. `export_benchmark_datasets.py` creates the fixed benchmark bank.
2. `run_abi_benchmark.ipynb` applies the trained ABI checkpoint to that bank.
3. `run_mcmc_benchmark.R` runs the model-matched DAGAR MCMC benchmark.
4. `compare_abi_vs_mcmc.ipynb` compares parameter recovery, boundary recovery,
   uncertainty, runtime, and break-even behavior.

The comparison notebook reads from:

```text
Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100/
|-- abi_results_all100/
|-- mcmc_results_all100/
`-- datasets/
```

It writes manuscript-ready outputs to:

```text
Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100/comparison_abi_vs_mcmc/
```

Important plot outputs include:

```text
comparison_abi_vs_mcmc/plots/parameter_recovery_truth_scatter.png
comparison_abi_vs_mcmc/plots/parameter_recovery_agreement_scatter.png
comparison_abi_vs_mcmc/plots/parameter_bias_interval_boxplots.png
comparison_abi_vs_mcmc/plots/boundary_metric_bars.png
comparison_abi_vs_mcmc/plots/runtime_comparison.png
```

### 4. Run the summary-statistic ablation study

Folder:

```text
Simulation Experiments/Ablation_experiments/
```

The ablation study evaluates whether each group of engineered summaries
contributes to posterior recovery and boundary detection. The baseline uses the
full summary representation. The ablations remove one block at a time or retain
only the core observation features.

Training notebooks:

```text
Simulation Experiments/Ablation_experiments/Training/00_baseline/train_baseline.ipynb
Simulation Experiments/Ablation_experiments/Training/01_ablate_core_observation_features/train_ablate_core_observation_features.ipynb
Simulation Experiments/Ablation_experiments/Training/02_ablate_graph_topology_features/train_ablate_graph_topology_features.ipynb
Simulation Experiments/Ablation_experiments/Training/03_ablate_dissimilarity_features/train_ablate_dissimilarity_features.ipynb
Simulation Experiments/Ablation_experiments/Training/04_ablate_local_spatial_features/train_ablate_local_spatial_features.ipynb
Simulation Experiments/Ablation_experiments/Training/05_ablate_global_graph_features/train_ablate_global_graph_features.ipynb
Simulation Experiments/Ablation_experiments/Training/06_core_observation_only/train_core_observation_only.ipynb
```

Comparison notebook:

```text
Simulation Experiments/Ablation_experiments/Results/compare_ablation_networks.ipynb
```

This notebook evaluates all available ablation checkpoints on the same shared
validation bank, including the 4050-dataset validation design with roughly 50
datasets for each graph size `N = 40, ..., 120`.

Main output folder:

```text
Simulation Experiments/Ablation_experiments/Results/comparison_outputs/
```

Important plots:

```text
comparison_outputs/Images/ablation_overall_summary.png
comparison_outputs/Images/ablation_recovery_heatmaps.png
comparison_outputs/Images/ablation_delta_vs_baseline_heatmaps.png
comparison_outputs/Images/ablation_boundary_heatmaps.png
comparison_outputs/Images/ablation_boundary_delta_vs_baseline_heatmaps.png
comparison_outputs/Images/ablation_training_histories.png
```

Important tables:

```text
comparison_outputs/Tables/overall_ablation_summary.csv
comparison_outputs/Tables/parameter_recovery_summary_display_parameters.csv
comparison_outputs/Tables/parameter_recovery_vs_baseline.csv
comparison_outputs/Tables/boundary_metric_summary.csv
comparison_outputs/Tables/boundary_metric_vs_baseline.csv
```

### 5. Real-data deployment and benchmark comparison

Folder:

```text
Real Data Analysis/
```

The real-data analysis applies the trained ABI-DAGAR network to the Glasgow and
California datasets and compares the resulting boundary conclusions with two
reference analyses:

1. CARBayes localized smoothing through `S.CARdissimilarity`.
2. A model-matched DAGARBayes sampler implemented with Rcpp.

Shared empirical inputs:

```text
Real Data Analysis/Data/respiratory_data.gpkg
Real Data Analysis/Data/respiratory_data_california.gpkg
Real Data Analysis/Data/adjacency_matrix.csv
Real Data Analysis/Data/adjacency_matrix_california.csv
```

#### CARBayes reference analysis

Producer script:

```text
Real Data Analysis/CARBayes_california_glasgow.R
```

This script runs the CARBayes benchmark for both datasets and writes the
exports consumed by the comparison notebook.

Required exported inputs:

```text
Real Data Analysis/setup_and_diagnostics/glasgow/glasgow_posterior_summary.csv
Real Data Analysis/setup_and_diagnostics/glasgow/glasgow_posterior_samples.csv
Real Data Analysis/setup_and_diagnostics/glasgow/glasgow_area_table.csv
Real Data Analysis/setup_and_diagnostics/glasgow/glasgow_edge_table.csv
Real Data Analysis/setup_and_diagnostics/california/california_posterior_summary.csv
Real Data Analysis/setup_and_diagnostics/california/california_posterior_samples.csv
Real Data Analysis/setup_and_diagnostics/california/california_area_table.csv
Real Data Analysis/setup_and_diagnostics/california/california_edge_table.csv
```

Comparison notebook:

```text
Real Data Analysis/CARBayes_vs_ABI_real_data.ipynb
```

Output folder:

```text
Real Data Analysis/results_ABI_vs_CARBayes/
```

Important plots:

```text
results_ABI_vs_CARBayes/glasgow/glasgow_boundary_agreement.png
results_ABI_vs_CARBayes/california/california_boundary_agreement.png
results_ABI_vs_CARBayes/glasgow/glasgow_edge_probability_vs_dissimilarity.png
results_ABI_vs_CARBayes/california/california_edge_probability_vs_dissimilarity.png
results_ABI_vs_CARBayes/glasgow/glasgow_fitted_risk_surface.png
results_ABI_vs_CARBayes/california/california_fitted_risk_surface.png
results_ABI_vs_CARBayes/glasgow/glasgow_post_pred_check.png
results_ABI_vs_CARBayes/california/california_post_pred_check.png
```

#### DAGARBayes reference analysis

Producer files:

```text
Real Data Analysis/DAGARBayes_california_glasgow.R
Real Data Analysis/dagar_poisson_boundary_mwg.cpp
```

The R script compiles the C++ sampler with `Rcpp::sourceCpp`, runs the
DAGARBayes benchmark for both datasets, and writes the exports consumed by the
comparison notebook.

Required exported inputs:

```text
Real Data Analysis/setup_and_diagnostics/DAGARBayes/glasgow/glasgow_posterior_summary.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/glasgow/glasgow_posterior_samples.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/glasgow/glasgow_area_table.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/glasgow/glasgow_edge_table.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/california/california_posterior_summary.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/california/california_posterior_samples.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/california/california_area_table.csv
Real Data Analysis/setup_and_diagnostics/DAGARBayes/california/california_edge_table.csv
```

Comparison notebook:

```text
Real Data Analysis/DAGARBayes_vs_ABI_real_data.ipynb
```

Output folder:

```text
Real Data Analysis/results_ABI_vs_DAGARBayes/
```

Important plots:

```text
results_ABI_vs_DAGARBayes/glasgow/glasgow_boundary_agreement_dagarbayes.png
results_ABI_vs_DAGARBayes/california/california_boundary_agreement_dagarbayes.png
results_ABI_vs_DAGARBayes/glasgow/glasgow_edge_probability_vs_dissimilarity.png
results_ABI_vs_DAGARBayes/california/california_edge_probability_vs_dissimilarity.png
results_ABI_vs_DAGARBayes/glasgow/glasgow_fitted_risk_surface.png
results_ABI_vs_DAGARBayes/california/california_fitted_risk_surface.png
results_ABI_vs_DAGARBayes/glasgow/glasgow_post_pred_check.png
results_ABI_vs_DAGARBayes/california/california_post_pred_check.png
```

## Environment

The ABI workflow was developed in Python with BayesFlow and Keras on the
TensorFlow backend. The recorded environment used:

```text
Python 3.10
BayesFlow 2.0.8
Keras 3.12.1
TensorFlow 2.21.0
NumPy
SciPy
Pandas
Matplotlib
GeoPandas
```

The R benchmark scripts use:

```text
R
CARBayes
sf
Rcpp
tictoc
```

The DAGARBayes and ABI-vs-MCMC benchmarks require a working Rcpp/C++ toolchain.

## Reproducibility Notes

Run notebooks from their own analysis folder unless a notebook explicitly
discovers the project root. Several configuration cells define paths near the
top of the notebook; if you clone the repository somewhere else, check those
cells first, especially paths pointing to the saved checkpoint.

The public repository keeps selected outputs used by the manuscript while
hiding exploratory notebooks, local IDE files, CAR-ABI side experiments, raw
heavy MCMC artifacts, and other intermediate files. The `.gitignore` file is
part of the reproducibility setup and documents what is intentionally excluded.

Folders with spaces, such as `Simulation Experiments` and `Real Data Analysis`,
should be quoted when used from a shell.

## Suggested Reading Order

If you are new to the repository, read/run in this order:

1. `Training/ABI_poisson_regression_DAGAR.ipynb`
2. `Simulation Experiments/simulation_more_detailed.ipynb`
3. `Simulation Experiments/ABI_vs_MCMC/compare_abi_vs_mcmc.ipynb`
4. `Simulation Experiments/Ablation_experiments/Results/compare_ablation_networks.ipynb`
5. `Real Data Analysis/CARBayes_vs_ABI_real_data.ipynb`
6. `Real Data Analysis/DAGARBayes_vs_ABI_real_data.ipynb`

If the saved outputs are already present, the comparison notebooks can be read
or rerun directly without regenerating every upstream MCMC or training artifact.

