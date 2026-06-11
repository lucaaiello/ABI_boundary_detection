## Reproducing the manuscript and supplementary results

This repository contains the code used to reproduce the results reported in the manuscript and Supplementary Materials for *Amortized Bayesian Boundary Detection on Varying-Size Graphs*.

The workflow has three main blocks:

1. **Simulation experiments** validating ABI-DAGAR under the proposed Poisson-DAGAR boundary-detection model. These include parameter recovery, posterior calibration, boundary-probability diagnostics, posterior predictive checks, the model-matched MCMC-DAGAR benchmark, computation summaries, and the ablation study of summary statistics.
2. **Real-data analyses** applying the trained ABI-DAGAR approximation to the Glasgow respiratory disease and California lung cancer datasets, with comparisons against `CARBayes` and a model-matched MCMC-DAGAR sampler.
3. **Supplementary diagnostics and benchmark analyses** extending the main validation with additional calibration plots, posterior predictive checks, ablation summaries, runtime summaries, and real-data model-matched MCMC-DAGAR comparisons.

Figures are stored next to the analysis block that generated them, not in a single global image folder. Numerical summaries for manuscript and supplementary tables are stored in the CSV outputs listed below.

### Path shortcuts

To keep the README readable, the output tables below use these path shortcuts.

| Shortcut | Expands to |
| --- | --- |
| `$TRAIN` | `Training` |
| `$SIM` | `Simulation Experiments` |
| `$SIM_IMG` | `Simulation Experiments/Images` |
| `$SIM_TAB` | `Simulation Experiments/Tables` |
| `$MCMC` | `Simulation Experiments/ABI_vs_MCMC` |
| `$MCMC_BANK` | `Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100` |
| `$MCMC_CMP` | `Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100/comparison_abi_vs_mcmc` |
| `$MCMC_PLOTS` | `Simulation Experiments/ABI_vs_MCMC/datasets/benchmark_bank_seed123_n100/comparison_abi_vs_mcmc/plots` |
| `$ABL` | `Simulation Experiments/Ablation_experiments` |
| `$ABL_IMG` | `Simulation Experiments/Ablation_experiments/Results/comparison_outputs/Images` |
| `$ABL_TAB` | `Simulation Experiments/Ablation_experiments/Results/comparison_outputs/Tables` |
| `$REAL` | `Real Data Analysis` |
| `$CAR_SETUP` | `Real Data Analysis/setup_and_diagnostics` |
| `$DAGAR_SETUP` | `Real Data Analysis/setup_and_diagnostics/DAGARBayes` |
| `$REAL_CAR` | `Real Data Analysis/results_ABI_vs_CARBayes` |
| `$REAL_DAGAR` | `Real Data Analysis/results_ABI_vs_DAGARBayes` |

### Recommended order of execution

Run the analysis blocks in the order below if reproducing the full workflow from scratch.

### 1. Train the ABI-DAGAR amortized posterior approximation

Run:

| File | Role |
| --- | --- |
| `$TRAIN/ABI_poisson_regression_DAGAR.ipynb` | Trains the SetTransformer summary network and conditional normalizing flow for ABI-DAGAR. |

Main outputs:

| Folder | File | Used in |
| --- | --- | --- |
| `$TRAIN/Images` | `history_plot.png` | Main manuscript Fig. 1 |
| `$TRAIN` | `training_history.csv` | Training-loss values and training diagnostics |
| `$TRAIN/Checkpoints` | `poisson_dagar.keras` | Trained ABI-DAGAR network used by downstream analyses |

### 2. Run the ABI-DAGAR simulation experiments

Run:

| File | Role |
| --- | --- |
| `$SIM/simulation_more_detailed.ipynb` | Performs held-out validation, posterior calibration, boundary-probability diagnostics, posterior predictive checks, and runtime summaries. |

Main manuscript figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$SIM_IMG` | `poisson_dagar_recovery.png` | Main manuscript Fig. 2 |
| `$SIM_IMG` | `poisson_dagar_calibration_histogram.png` | Main manuscript Fig. 3 |
| `$SIM_IMG` | `boundary_sens_spec_curve.png` | Main manuscript Fig. 4 |
| `$SIM_IMG` | `boundary_mpm_histograms.png` | Main manuscript Fig. 5 |

Supplementary figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$SIM_IMG` | `parameter_coverage_curve.png` | Supplement Fig. S1 |
| `$SIM_IMG` | `poisson_dagar_calibration_ecdf.png` | Supplement Fig. S2 |
| `$SIM_IMG` | `graph_size_posterior_zscores.png` | Supplement Fig. S3 |
| `$SIM_IMG` | `regime_specific_parameter_errors.png` | Supplement Fig. S4 |
| `$SIM_IMG` | `boundary_probability_quality.png` | Supplement Fig. S5 |
| `$SIM_IMG` | `predictive_checks.png` | Supplement Fig. S6 |

Numerical summary sources:

| Result | Source files |
| --- | --- |
| Main simulation parameter-recovery table | `$SIM_TAB/parameter_recovery_summary.csv`; `$SIM_TAB/parameter_recovery_detail.csv` |
| Additional boundary-probability table | `$SIM_TAB/boundary_probability_summary.csv`; `$SIM_TAB/boundary_probability_detail.csv`; `$SIM/boundary_mpm_summary.csv`; `$SIM/boundary_mpm_results.csv` |
| Posterior predictive table | `$SIM_TAB/predictive_check_summary.csv`; `$SIM_TAB/predictive_check_detail.csv` |
| Simulation computation table | `$SIM_TAB/computation_summary.csv`; `$SIM_TAB/hardware_summary.csv`; `$TRAIN/training_history.csv` |

### 3. Run the simulation benchmark against MCMC-DAGAR

Run:

| File | Role |
| --- | --- |
| `$MCMC/export_benchmark_datasets.py` | Exports the fixed benchmark bank of simulated datasets. |
| `$MCMC/run_abi_benchmark.ipynb` | Runs ABI-DAGAR on the benchmark datasets. |
| `$MCMC/run_mcmc_benchmark.R` | Runs the model-matched MCMC-DAGAR sampler on the same datasets. |
| `$MCMC/compare_abi_vs_mcmc.ipynb` | Produces the ABI-DAGAR versus MCMC-DAGAR comparison tables and plots. |

Key input/output folders:

| Folder | Contents |
| --- | --- |
| `$MCMC_BANK/datasets` | Fixed simulated datasets used by both methods. |
| `$MCMC_BANK/abi_results_all100` | ABI-DAGAR posterior summaries and edge probabilities. |
| `$MCMC_BANK/mcmc_results_all100` | MCMC-DAGAR posterior summaries and edge probabilities. |
| `$MCMC_CMP` | Comparison CSV summaries. |
| `$MCMC_PLOTS` | Comparison figures. |

Main manuscript figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$MCMC_PLOTS` | `parameter_recovery_bars.png` | Main manuscript Fig. 6 |
| `$MCMC_PLOTS` | `parameter_recovery_truth_scatter.png` | Main manuscript Fig. 7 |
| `$MCMC_PLOTS` | `parameter_recovery_agreement_scatter.png` | Main manuscript Fig. 7 |
| `$MCMC_PLOTS` | `parameter_bias_interval_boxplots.png` | Main manuscript Fig. 7 |

Supplementary figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$MCMC_PLOTS` | `boundary_metric_bars.png` | Supplement Fig. S7 |
| `$MCMC_PLOTS` | `runtime_comparison.png` | Supplement Fig. S7 |

Numerical summary sources:

| Result | Source files |
| --- | --- |
| Parameter-level ABI-DAGAR versus MCMC-DAGAR comparison | `$MCMC_CMP/parameter_method_summary.csv`; `$MCMC_CMP/parameter_agreement_summary.csv`; `$MCMC_CMP/parameter_pairwise_by_dataset.csv` |
| Boundary-probability and decision-rule comparison | `$MCMC_CMP/edge_metric_method_summary.csv`; `$MCMC_CMP/edge_agreement_summary.csv`; `$MCMC_CMP/edge_agreement_by_dataset.csv` |
| Runtime and break-even summaries | `$MCMC_CMP/runtime_summary.csv`; `$MCMC_CMP/runtime_comparison_by_dataset.csv`; `$MCMC_CMP/break_even_summary.csv`; `$MCMC_CMP/break_even_scenarios.csv` |

### 4. Run the ablation study of summary statistics

Run the ablation training notebooks and then the comparison notebook:

| File pattern | Role |
| --- | --- |
| `$ABL/Training/*/train_*.ipynb` | Trains the baseline and ablated networks using shared training and validation banks. |
| `$ABL/Results/compare_ablation_networks.ipynb` | Evaluates all available ablation networks on the shared validation set and creates comparison outputs. |

Main manuscript figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$ABL_IMG` | `ablation_delta_vs_baseline_heatmaps.png` | Main manuscript Fig. 8 |
| `$ABL_IMG` | `ablation_boundary_delta_vs_baseline_heatmaps.png` | Main manuscript Fig. 9 |

Supplementary figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$ABL_IMG` | `ablation_recovery_heatmaps.png` | Supplement Fig. S8 |
| `$ABL_IMG` | `ablation_boundary_heatmaps.png` | Supplement Fig. S9 |
| `$ABL_IMG` | `ablation_training_histories.png` | Supplement Fig. S10 |

Additional ablation outputs:

| Folder | File | Contents |
| --- | --- | --- |
| `$ABL_IMG` | `ablation_overall_summary.png` | Compact overall ablation summary. |
| `$ABL_TAB` | `overall_ablation_summary.csv` | Aggregate recovery metrics by summary representation. |
| `$ABL_TAB` | `parameter_recovery_summary_display_parameters.csv` | Parameter-specific recovery metrics for displayed parameters. |
| `$ABL_TAB` | `parameter_recovery_vs_baseline.csv` | Parameter-recovery differences versus the full-summary baseline. |
| `$ABL_TAB` | `boundary_metric_summary.csv` | Boundary-detection metrics by summary representation. |
| `$ABL_TAB` | `boundary_metric_vs_baseline.csv` | Boundary-detection differences versus the full-summary baseline. |
| `$ABL_TAB` | `training_history_summary.csv` | Final training and validation losses by ablation run. |

### 5. Run the real-data ABI-DAGAR and `CARBayes` analyses

Run:

| File | Role |
| --- | --- |
| `$REAL/CARBayes_california_glasgow.R` | Fits the localized `CARBayes` benchmark for Glasgow and California and saves posterior/edge diagnostics. |
| `$REAL/CARBayes_vs_ABI_real_data.ipynb` | Applies ABI-DAGAR to the real datasets and compares ABI-DAGAR with `CARBayes`. |

Input/result folders:

| Folder | Contents |
| --- | --- |
| `$REAL/Data` | Glasgow and California spatial data and adjacency files. |
| `$CAR_SETUP/glasgow` and `$CAR_SETUP/california` | `CARBayes` posterior samples, posterior summaries, edge tables, and boundary probabilities. |
| `$REAL_CAR/glasgow` and `$REAL_CAR/california` | ABI-DAGAR outputs and ABI-DAGAR versus `CARBayes` comparison outputs. |

Main manuscript figure:

| Folder | File | Used in |
| --- | --- | --- |
| `$REAL_CAR/glasgow` | `glasgow_boundary_agreement.png` | Main manuscript Fig. 10 |
| `$REAL_CAR/california` | `california_boundary_agreement.png` | Main manuscript Fig. 10 |

Supplementary figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$REAL_CAR/glasgow` | `glasgow_post_pred_check.png` | Supplement Fig. S11 |
| `$REAL_CAR/california` | `california_post_pred_check.png` | Supplement Fig. S11 |
| `$REAL_CAR/glasgow` | `glasgow_risk_comparison.png` | Supplement Fig. S12 |
| `$REAL_CAR/california` | `california_risk_comparison.png` | Supplement Fig. S12 |
| `$REAL_CAR/glasgow` | `glasgow_edge_probability_vs_dissimilarity.png` | Supplement Fig. S13 |
| `$REAL_CAR/california` | `california_edge_probability_vs_dissimilarity.png` | Supplement Fig. S13 |

Numerical summary sources:

| Result | Source files |
| --- | --- |
| Main posterior-summary comparison table | `$REAL_CAR/glasgow/parameter_summary_comparison.csv`; `$REAL_CAR/california/parameter_summary_comparison.csv`; benchmark inputs in `$CAR_SETUP/glasgow` and `$CAR_SETUP/california` |
| Main boundary-agreement table | `$REAL_CAR/combined_edge_metrics.csv`; `$REAL_CAR/glasgow/edge_metrics.csv`; `$REAL_CAR/california/edge_metrics.csv` |
| Supplement fitted-risk comparison table | `$REAL_CAR/combined_risk_metrics.csv`; `$REAL_CAR/glasgow/risk_metrics.csv`; `$REAL_CAR/california/risk_metrics.csv`; per-edge/per-area values in `risk_comparison.csv` |
| Supplement boundary-probability table | `$REAL_CAR/combined_edge_metrics.csv`; `$REAL_CAR/glasgow/edge_comparison.csv`; `$REAL_CAR/california/edge_comparison.csv` |
| Supplement real-data runtime table | Timing information printed by `$REAL/CARBayes_california_glasgow.R` and `$REAL/CARBayes_vs_ABI_real_data.ipynb` |

### 6. Run the real-data model-matched MCMC-DAGAR comparison

Run:

| File | Role |
| --- | --- |
| `$REAL/DAGARBayes_california_glasgow.R` | Fits the model-matched MCMC-DAGAR sampler for Glasgow and California. |
| `$REAL/dagar_poisson_boundary_mwg.cpp` | C++ sampler used by the MCMC-DAGAR real-data script. |
| `$REAL/DAGARBayes_vs_ABI_real_data.ipynb` | Compares ABI-DAGAR with MCMC-DAGAR on the real datasets. |

Input/result folders:

| Folder | Contents |
| --- | --- |
| `$DAGAR_SETUP/glasgow` and `$DAGAR_SETUP/california` | MCMC-DAGAR posterior samples, posterior summaries, acceptance summaries, edge tables, and boundary probabilities. |
| `$REAL_DAGAR/glasgow` and `$REAL_DAGAR/california` | ABI-DAGAR versus MCMC-DAGAR comparison outputs. |

Supplementary figures:

| Folder | File | Used in |
| --- | --- | --- |
| `$REAL_DAGAR/glasgow` | `glasgow_boundary_agreement_dagarbayes.png` | Supplement Fig. S14 |
| `$REAL_DAGAR/california` | `california_boundary_agreement_dagarbayes.png` | Supplement Fig. S14 |

Numerical summary sources:

| Result | Source files |
| --- | --- |
| Supplement ABI-DAGAR versus MCMC-DAGAR posterior-summary table | `$REAL_DAGAR/glasgow/parameter_summary_comparison.csv`; `$REAL_DAGAR/california/parameter_summary_comparison.csv`; MCMC inputs in `$DAGAR_SETUP/glasgow` and `$DAGAR_SETUP/california` |
| Supplement ABI-DAGAR versus MCMC-DAGAR agreement table | `$REAL_DAGAR/combined_edge_metrics.csv`; `$REAL_DAGAR/glasgow/edge_metrics.csv`; `$REAL_DAGAR/california/edge_metrics.csv` |

### Main manuscript output index

| Item | Folder | File(s) | Contents |
| --- | --- | --- | --- |
| Fig. 1 | `$TRAIN/Images` | `history_plot.png` | Training loss over 100 epochs for ABI-DAGAR. |
| Fig. 2 | `$SIM_IMG` | `poisson_dagar_recovery.png` | Parameter recovery on 200 held-out simulated datasets. |
| Fig. 3 | `$SIM_IMG` | `poisson_dagar_calibration_histogram.png` | Simulation-based calibration rank histograms. |
| Fig. 4 | `$SIM_IMG` | `boundary_sens_spec_curve.png` | Sensitivity and specificity as a function of the number of selected boundaries. |
| Fig. 5 | `$SIM_IMG` | `boundary_mpm_histograms.png` | Sensitivity and specificity under the median-probability rule. |
| Fig. 6 | `$MCMC_PLOTS` | `parameter_recovery_bars.png` | ABI-DAGAR versus MCMC-DAGAR mean absolute error and empirical coverage. |
| Fig. 7 | `$MCMC_PLOTS` | `parameter_recovery_truth_scatter.png`; `parameter_recovery_agreement_scatter.png`; `parameter_bias_interval_boxplots.png` | Detailed ABI-DAGAR versus MCMC-DAGAR parameter-level comparison. |
| Fig. 8 | `$ABL_IMG` | `ablation_delta_vs_baseline_heatmaps.png` | Scalar recovery differences relative to the full-summary baseline. |
| Fig. 9 | `$ABL_IMG` | `ablation_boundary_delta_vs_baseline_heatmaps.png` | Boundary-detection differences relative to the full-summary baseline. |
| Fig. 10 | `$REAL_CAR/glasgow`; `$REAL_CAR/california` | `glasgow_boundary_agreement.png`; `california_boundary_agreement.png` | Boundary agreement between ABI-DAGAR and `CARBayes`. |

### Main manuscript table sources

| Item | Produced by | Numerical source files |
| --- | --- | --- |
| Table 1 | `$SIM/simulation_more_detailed.ipynb` | `$SIM_TAB/parameter_recovery_summary.csv`; `$SIM_TAB/parameter_recovery_detail.csv` |
| Table 2 | `$REAL/CARBayes_california_glasgow.R`; `$REAL/CARBayes_vs_ABI_real_data.ipynb` | `$REAL_CAR/glasgow/parameter_summary_comparison.csv`; `$REAL_CAR/california/parameter_summary_comparison.csv`; benchmark inputs in `$CAR_SETUP/glasgow` and `$CAR_SETUP/california` |
| Table 3 | `$REAL/CARBayes_vs_ABI_real_data.ipynb` | `$REAL_CAR/combined_edge_metrics.csv`; `$REAL_CAR/glasgow/edge_metrics.csv`; `$REAL_CAR/california/edge_metrics.csv` |

### Supplementary output index

| Item | Folder | File(s) | Contents |
| --- | --- | --- | --- |
| Fig. S1 | `$SIM_IMG` | `parameter_coverage_curve.png` | Empirical coverage versus nominal coverage. |
| Fig. S2 | `$SIM_IMG` | `poisson_dagar_calibration_ecdf.png` | Simulation-based calibration ECDF-difference plots. |
| Fig. S3 | `$SIM_IMG` | `graph_size_posterior_zscores.png` | Posterior z-scores stratified by graph-size bin. |
| Fig. S4 | `$SIM_IMG` | `regime_specific_parameter_errors.png` | Mean absolute recovery error by true rho, true eta, graph size, and edge count. |
| Fig. S5 | `$SIM_IMG` | `boundary_probability_quality.png` | Boundary-probability reliability and dissimilarity-profile diagnostics. |
| Fig. S6 | `$SIM_IMG` | `predictive_checks.png` | Posterior predictive diagnostics on held-out simulated datasets. |
| Fig. S7 | `$MCMC_PLOTS` | `boundary_metric_bars.png`; `runtime_comparison.png` | Additional ABI-DAGAR versus MCMC-DAGAR boundary and runtime diagnostics. |
| Fig. S8 | `$ABL_IMG` | `ablation_recovery_heatmaps.png` | Parameter-specific recovery metrics for the ablation study. |
| Fig. S9 | `$ABL_IMG` | `ablation_boundary_heatmaps.png` | Boundary-detection diagnostics across summary representations. |
| Fig. S10 | `$ABL_IMG` | `ablation_training_histories.png` | Training and validation loss trajectories for the ablation study. |
| Fig. S11 | `$REAL_CAR/glasgow`; `$REAL_CAR/california` | `glasgow_post_pred_check.png`; `california_post_pred_check.png` | Posterior predictive checks for the real-data applications. |
| Fig. S12 | `$REAL_CAR/glasgow`; `$REAL_CAR/california` | `glasgow_risk_comparison.png`; `california_risk_comparison.png` | Fitted-risk comparison between ABI-DAGAR and `CARBayes`. |
| Fig. S13 | `$REAL_CAR/glasgow`; `$REAL_CAR/california` | `glasgow_edge_probability_vs_dissimilarity.png`; `california_edge_probability_vs_dissimilarity.png` | Posterior boundary probability versus standardized edge dissimilarity. |
| Fig. S14 | `$REAL_DAGAR/glasgow`; `$REAL_DAGAR/california` | `glasgow_boundary_agreement_dagarbayes.png`; `california_boundary_agreement_dagarbayes.png` | Boundary agreement between ABI-DAGAR and model-matched MCMC-DAGAR. |

### Supplementary table sources

| Item | Produced by | Numerical source files |
| --- | --- | --- |
| Table S1 | `$SIM/simulation_more_detailed.ipynb` | `$SIM_TAB/boundary_probability_summary.csv`; `$SIM_TAB/boundary_probability_detail.csv`; `$SIM/boundary_mpm_summary.csv`; `$SIM/boundary_mpm_results.csv` |
| Table S2 | `$SIM/simulation_more_detailed.ipynb` | `$SIM_TAB/predictive_check_summary.csv`; `$SIM_TAB/predictive_check_detail.csv` |
| Table S3 | `$SIM/simulation_more_detailed.ipynb`; `$TRAIN/ABI_poisson_regression_DAGAR.ipynb` | `$SIM_TAB/computation_summary.csv`; `$SIM_TAB/hardware_summary.csv`; `$TRAIN/training_history.csv` |
| Table S4 | `$REAL/CARBayes_vs_ABI_real_data.ipynb` | `$REAL_CAR/combined_risk_metrics.csv`; `$REAL_CAR/glasgow/risk_metrics.csv`; `$REAL_CAR/california/risk_metrics.csv` |
| Table S5 | `$REAL/CARBayes_vs_ABI_real_data.ipynb` | `$REAL_CAR/combined_edge_metrics.csv`; `$REAL_CAR/glasgow/edge_metrics.csv`; `$REAL_CAR/california/edge_metrics.csv` |
| Table S6 | `$REAL/DAGARBayes_california_glasgow.R`; `$REAL/DAGARBayes_vs_ABI_real_data.ipynb` | `$REAL_DAGAR/glasgow/parameter_summary_comparison.csv`; `$REAL_DAGAR/california/parameter_summary_comparison.csv`; MCMC inputs in `$DAGAR_SETUP/glasgow` and `$DAGAR_SETUP/california` |
| Table S7 | `$REAL/DAGARBayes_vs_ABI_real_data.ipynb` | `$REAL_DAGAR/combined_edge_metrics.csv`; `$REAL_DAGAR/glasgow/edge_metrics.csv`; `$REAL_DAGAR/california/edge_metrics.csv` |
| Table S8 | `$REAL/CARBayes_california_glasgow.R`; `$REAL/CARBayes_vs_ABI_real_data.ipynb` | Runtime values printed by the corresponding R script and notebook cells |
| Table S9 | `$TRAIN/ABI_poisson_regression_DAGAR.ipynb` | Network architecture, software versions, training settings, and `$TRAIN/training_history.csv` |

### Finding outputs quickly

From the repository root, useful searches are:

```bash
rg --files | rg "poisson_dagar_recovery.png|parameter_recovery_summary.csv"
rg --files | rg "boundary_metric_bars.png|edge_metric_method_summary.csv"
rg --files | rg "glasgow_boundary_agreement|parameter_summary_comparison.csv"
```

For a full rerun, start with the training notebook, then run the simulation validation, benchmark comparison, ablation comparison, and real-data notebooks in the order listed above.
