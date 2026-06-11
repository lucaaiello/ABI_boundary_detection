## Reproducing the manuscript results

This repository contains the code used to reproduce all results reported in the manuscript and Supplementary Materials for *Amortized Bayesian Boundary Detection on Varying-Size Graphs*. The analyses follow the same order as the paper: first train the ABI-DAGAR amortized posterior approximation, then validate it on held-out simulated datasets, then compare it with model-matched MCMC-DAGAR and `CARBayes`, and finally run the additional ablation and ABI-CAR sensitivity analyses.

The final figures used in the manuscript and Supplementary Materials are stored in the `Images/` folder. The numerical values reported in tables are contained in the corresponding LaTeX source files, `main.tex` and `supplementary.tex`, and are produced by the simulation, benchmark, posterior predictive, ablation, and real-data analysis code in this repository.

### Recommended order of execution

The code should be run in the following order.

#### 1. Train the ABI-DAGAR amortized posterior approximation

Run the ABI-DAGAR training code first. This step simulates training datasets from the Poisson-DAGAR boundary-detection model and trains the SetTransformer summary network and conditional normalizing flow.

This step produces the trained amortized posterior approximation and the training-history output used in the manuscript.

Main output:

| Output                    | Used in                |
| ------------------------- | ---------------------- |
| `Images/history_plot.png` | Main manuscript Fig. 1 |

#### 2. Run held-out ABI-DAGAR simulation validation

After training, run the held-out simulation validation. This step applies the trained amortized posterior approximation to 200 held-out simulated datasets and computes parameter recovery, posterior calibration, posterior boundary probabilities, posterior predictive checks, and computational summaries.

Main manuscript outputs from this step:

| Output                                            | Used in                 |
| ------------------------------------------------- | ----------------------- |
| `Images/poisson_dagar_recovery.png`               | Main manuscript Fig. 2  |
| `Images/poisson_dagar_calibration_histogram.png`  | Main manuscript Fig. 3  |
| `Images/boundary_sens_spec_curve.png`             | Main manuscript Fig. 4  |
| `Images/boundary_mpm_histograms.png`              | Main manuscript Fig. 5  |
| Table values for `tab:sim_recovery` in `main.tex` | Main manuscript Table 1 |

Supplementary outputs from this step:

| Output                                                         | Used in             |
| -------------------------------------------------------------- | ------------------- |
| `Images/parameter_coverage_curve.png`                          | Supplement Fig. S1  |
| `Images/poisson_dagar_calibration_ecdf.png`                    | Supplement Fig. S2  |
| `Images/graph_size_posterior_zscores.png`                      | Supplement Fig. S3  |
| `Images/regime_specific_parameter_errors.png`                  | Supplement Fig. S4  |
| `Images/boundary_probability_quality.png`                      | Supplement Fig. S5  |
| `Images/predictive_checks.png`                                 | Supplement Fig. S6  |
| Table values for `tab:sim_boundary_sup` in `supplementary.tex` | Supplement Table S1 |
| Table values for `tab:sim_ppc_sup` in `supplementary.tex`      | Supplement Table S2 |
| Table values for `tab:sim_computation` in `supplementary.tex`  | Supplement Table S3 |

#### 3. Run the model-matched MCMC-DAGAR simulation benchmark

Next, run the MCMC-DAGAR benchmark on the additional bank of 100 held-out simulated datasets. This benchmark compares ABI-DAGAR with a dataset-specific MCMC sampler targeting the same thresholded Poisson-DAGAR model.

Main manuscript outputs from this step:

| Output                                            | Used in                |
| ------------------------------------------------- | ---------------------- |
| `Images/parameter_recovery_bars.png`              | Main manuscript Fig. 6 |
| `Images/parameter_recovery_truth_scatter.png`     | Main manuscript Fig. 7 |
| `Images/parameter_recovery_agreement_scatter.png` | Main manuscript Fig. 7 |
| `Images/parameter_bias_interval_boxplots.png`     | Main manuscript Fig. 7 |

Supplementary outputs from this step:

| Output                            | Used in            |
| --------------------------------- | ------------------ |
| `Images/boundary_metric_bars.png` | Supplement Fig. S7 |
| `Images/runtime_comparison.png`   | Supplement Fig. S7 |

#### 4. Run the summary-statistic ablation study

The ablation study should be run after the baseline ABI-DAGAR validation, because all ablated representations are compared with the full-summary baseline. For each summary representation, a separate amortized posterior approximator is retrained from scratch using the same network architecture and training protocol, and all runs are evaluated on the same validation datasets.

Outputs from this step:

| Output                                   | Used in             |
| ---------------------------------------- | ------------------- |
| `Images/ablation_recovery_heatmaps.png`  | Supplement Fig. S8  |
| `Images/ablation_boundary_heatmaps.png`  | Supplement Fig. S9  |
| `Images/ablation_training_histories.png` | Supplement Fig. S10 |

#### 5. Run the real-data ABI-DAGAR analyses and `CARBayes` benchmark

After validating the trained ABI-DAGAR amortizer on simulated datasets, apply the same trained network to the Glasgow respiratory disease and California lung cancer datasets. Then run the localized `CARBayes` benchmark on the same datasets using the same standardized dissimilarity covariates.

Main manuscript outputs from this step:

| Output                                                             | Used in                                             |
| ------------------------------------------------------------------ | --------------------------------------------------- |
| `Images/glasgow_boundary_agreement.png`                            | Main manuscript real-data boundary agreement figure |
| `Images/california_boundary_agreement.png`                         | Main manuscript real-data boundary agreement figure |
| Table values for `tab:carbayes-parameter-comparison` in `main.tex` | Main manuscript Table 2                             |
| Table values for `tab:carbayes-agreement` in `main.tex`            | Main manuscript Table 3                             |

Supplementary outputs from this step:

| Output                                                               | Used in             |
| -------------------------------------------------------------------- | ------------------- |
| `Images/glasgow_post_pred_check.png`                                 | Supplement Fig. S11 |
| `Images/california_post_pred_check.png`                              | Supplement Fig. S11 |
| `Images/glasgow_risk_comparison.png`                                 | Supplement Fig. S12 |
| `Images/california_risk_comparison.png`                              | Supplement Fig. S12 |
| `Images/glasgow_edge_probability_vs_dissimilarity.png`               | Supplement Fig. S13 |
| `Images/california_edge_probability_vs_dissimilarity.png`            | Supplement Fig. S13 |
| Table values for `tab:realdata_risk_extra` in `supplementary.tex`    | Supplement Table S4 |
| Table values for `tab:realdata_edge_extra` in `supplementary.tex`    | Supplement Table S5 |
| Table values for `tab:realdata_runtime_extra` in `supplementary.tex` | Supplement Table S8 |

#### 6. Run the real-data model-matched MCMC-DAGAR benchmark

The real-data MCMC-DAGAR benchmark compares ABI-DAGAR with a dataset-specific MCMC implementation of the same DAGAR boundary model. This step is separate from the external `CARBayes` benchmark and is used to distinguish model-class differences from amortization error.

Outputs from this step:

| Output                                                                 | Used in             |
| ---------------------------------------------------------------------- | ------------------- |
| `Images/glasgow_boundary_agreement_dagarbayes.png`                     | Supplement Fig. S14 |
| `Images/california_boundary_agreement_dagarbayes.png`                  | Supplement Fig. S14 |
| Table values for MCMC-DAGAR posterior summaries in `supplementary.tex` | Supplement Table S6 |
| Table values for `tab:mcmc_dagar_agreement` in `supplementary.tex`     | Supplement Table S7 |

#### 7. Run the ABI-CAR sensitivity analysis

The ABI-CAR analysis is an additional sensitivity experiment. It repeats the amortized Bayesian workflow under a localized Leroux CAR prior matching the spatial-dependence specification used by Lee and Mitchell (2012) and the `CARBayes` benchmark. This analysis serves as a positive-control experiment: when the amortized model is aligned with the localized CAR prior, the resulting boundary conclusions closely reproduce the `CARBayes` benchmark.

Outputs from this step:

| Output                                                                             | Used in              |
| ---------------------------------------------------------------------------------- | -------------------- |
| `Images/poisson_car_recovery.png`                                                  | Supplement Fig. S15  |
| `Images/parameter_coverage_curve_car.png`                                          | Supplement Fig. S15  |
| `Images/glasgow_boundary_agreement_car.png`                                        | Supplement Fig. S16  |
| `Images/california_boundary_agreement_car.png`                                     | Supplement Fig. S16  |
| Table values for `tab:supp_abicar_recovery` in `supplementary.tex`                 | Supplement Table S10 |
| Table values for ABI-CAR real-data posterior summaries in `supplementary.tex`      | Supplement Table S11 |
| Table values for ABI-CAR boundary and fitted-risk agreement in `supplementary.tex` | Supplement Table S12 |

#### 8. Compile the manuscript and Supplementary Materials

After the figures and numerical summaries have been generated, compile `main.tex` and `supplementary.tex`. The manuscript and supplement read the final figures from the `Images/` folder. The table values are included directly in the LaTeX source files.

Because the main manuscript and supplement use cross-references between files, compile both documents enough times for cross-references to resolve correctly.

### Output index

The following tables provide a direct index between the manuscript outputs and the corresponding files in the repository.

#### Main manuscript figures

| Manuscript item                     | File(s) in repository                                                                                                                           | Contents                                                                                            |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| Fig. 1                              | `Images/history_plot.png`                                                                                                                       | Training loss over 100 epochs for the amortized posterior approximation.                            |
| Fig. 2                              | `Images/poisson_dagar_recovery.png`                                                                                                             | Parameter recovery on 200 held-out simulated datasets.                                              |
| Fig. 3                              | `Images/poisson_dagar_calibration_histogram.png`                                                                                                | Simulation-based calibration rank histograms for the four model parameters.                         |
| Fig. 4                              | `Images/boundary_sens_spec_curve.png`                                                                                                           | Sensitivity and specificity as functions of the number of selected boundaries.                      |
| Fig. 5                              | `Images/boundary_mpm_histograms.png`                                                                                                            | Sensitivity and specificity under the median-probability rule.                                      |
| Fig. 6                              | `Images/parameter_recovery_bars.png`                                                                                                            | ABI-DAGAR versus model-matched MCMC-DAGAR: mean absolute error and empirical 95% interval coverage. |
| Fig. 7                              | `Images/parameter_recovery_truth_scatter.png`; `Images/parameter_recovery_agreement_scatter.png`; `Images/parameter_bias_interval_boxplots.png` | Detailed parameter-level comparison between ABI-DAGAR and model-matched MCMC-DAGAR.                 |
| Real-data boundary agreement figure | `Images/glasgow_boundary_agreement.png`; `Images/california_boundary_agreement.png`                                                             | Boundary agreement between ABI-DAGAR and `CARBayes` in the Glasgow and California applications.     |

#### Main manuscript tables

| Manuscript item | Location in source                                          | Contents                                                                                                                   |
| --------------- | ----------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| Table 1         | `main.tex`, table label `tab:sim_recovery`                  | Parameter recovery on 200 held-out simulated datasets: bias, RMSE, Pearson correlation, (R^2), and empirical 95% coverage. |
| Table 2         | `main.tex`, table label `tab:carbayes-parameter-comparison` | Posterior medians and 95% credible intervals for ABI-DAGAR and `CARBayes` in the Glasgow and California applications.      |
| Table 3         | `main.tex`, table label `tab:carbayes-agreement`            | Edge-level and fitted-risk agreement between ABI-DAGAR and `CARBayes` under the median-probability rule.                   |

#### Supplementary figures

| Supplement item | File(s) in repository                                                                                             | Contents                                                                                                                 |
| --------------- | ----------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| Fig. S1         | `Images/parameter_coverage_curve.png`                                                                             | Empirical coverage versus nominal coverage for the amortized posterior approximation.                                    |
| Fig. S2         | `Images/poisson_dagar_calibration_ecdf.png`                                                                       | Simulation-based calibration ECDF-difference plots.                                                                      |
| Fig. S3         | `Images/graph_size_posterior_zscores.png`                                                                         | Posterior (z)-scores stratified by graph-size bin.                                                                       |
| Fig. S4         | `Images/regime_specific_parameter_errors.png`                                                                     | Mean absolute recovery error by true (\rho), true (\eta), graph size, and edge count.                                    |
| Fig. S5         | `Images/boundary_probability_quality.png`                                                                         | Boundary-probability reliability and dissimilarity-profile diagnostics.                                                  |
| Fig. S6         | `Images/predictive_checks.png`                                                                                    | Posterior predictive diagnostics for counts, Moran-type summaries, edge contrasts, and the implied number of boundaries. |
| Fig. S7         | `Images/boundary_metric_bars.png`; `Images/runtime_comparison.png`                                                | Additional ABI-DAGAR and MCMC-DAGAR benchmark diagnostics and runtime comparison.                                        |
| Fig. S8         | `Images/ablation_recovery_heatmaps.png`                                                                           | Parameter-specific recovery metrics for the ablation study.                                                              |
| Fig. S9         | `Images/ablation_boundary_heatmaps.png`                                                                           | Boundary-detection diagnostics across summary representations.                                                           |
| Fig. S10        | `Images/ablation_training_histories.png`                                                                          | Training and validation loss trajectories for the baseline and ablated summary representations.                          |
| Fig. S11        | `Images/glasgow_post_pred_check.png`; `Images/california_post_pred_check.png`                                     | Posterior predictive checks for the Glasgow and California real-data applications.                                       |
| Fig. S12        | `Images/glasgow_risk_comparison.png`; `Images/california_risk_comparison.png`                                     | Fitted-risk comparison between ABI-DAGAR and `CARBayes`.                                                                 |
| Fig. S13        | `Images/glasgow_edge_probability_vs_dissimilarity.png`; `Images/california_edge_probability_vs_dissimilarity.png` | Posterior boundary probability versus standardized edge dissimilarity.                                                   |
| Fig. S14        | `Images/glasgow_boundary_agreement_dagarbayes.png`; `Images/california_boundary_agreement_dagarbayes.png`         | Boundary agreement between ABI-DAGAR and the model-matched MCMC-DAGAR benchmark.                                         |
| Fig. S15        | `Images/poisson_car_recovery.png`; `Images/parameter_coverage_curve_car.png`                                      | ABI-CAR simulation diagnostics.                                                                                          |
| Fig. S16        | `Images/glasgow_boundary_agreement_car.png`; `Images/california_boundary_agreement_car.png`                       | Boundary agreement between ABI-CAR and `CARBayes`.                                                                       |

#### Supplementary tables

| Supplement item | Location in source                                                                                | Contents                                                                                 |
| --------------- | ------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| Table S1        | `supplementary.tex`, table label `tab:sim_boundary_sup`                                           | Additional boundary-probability summaries on held-out simulated datasets.                |
| Table S2        | `supplementary.tex`, table label `tab:sim_ppc_sup`                                                | Additional posterior predictive summaries on held-out simulated datasets.                |
| Table S3        | `supplementary.tex`, table label `tab:sim_computation`                                            | Computational details for the simulation study.                                          |
| Table S4        | `supplementary.tex`, table label `tab:realdata_risk_extra`                                        | Additional fitted-risk comparison metrics between ABI-DAGAR and `CARBayes`.              |
| Table S5        | `supplementary.tex`, table label `tab:realdata_edge_extra`                                        | Additional boundary-probability comparison metrics between ABI-DAGAR and `CARBayes`.     |
| Table S6        | `supplementary.tex`, subsection `Additional comparison with a model-matched MCMC-DAGAR benchmark` | Posterior medians and 95% credible intervals for ABI-DAGAR and model-matched MCMC-DAGAR. |
| Table S7        | `supplementary.tex`, table label `tab:mcmc_dagar_agreement`                                       | Agreement between ABI-DAGAR and model-matched MCMC-DAGAR.                                |
| Table S8        | `supplementary.tex`, table label `tab:realdata_runtime_extra`                                     | Real-data posterior sampling runtimes for ABI-DAGAR and `CARBayes`.                      |
| Table S9        | `supplementary.tex`, table label `tab:supp_network_impl`                                          | Neural-network and training configuration for the ABI-DAGAR implementation.              |
| Table S10       | `supplementary.tex`, table label `tab:supp_abicar_recovery`                                       | ABI-CAR parameter recovery on 200 held-out simulated datasets.                           |
| Table S11       | `supplementary.tex`, subsection `Real-data comparison with CARBayes`                              | Posterior medians and 95% credible intervals for ABI-CAR and `CARBayes`.                 |
| Table S12       | `supplementary.tex`, subsection `Real-data comparison with CARBayes`                              | Boundary and fitted-risk agreement between ABI-CAR and `CARBayes`.                       |

### Locating outputs from the command line

The final figure files can be found directly in the `Images/` folder. The numerical table entries can be located in the LaTeX source by searching for the table labels listed above, for example:

```bash
grep -n "tab:sim_recovery" main.tex
grep -n "tab:carbayes-agreement" main.tex
grep -n "tab:sim_boundary_sup" supplementary.tex
grep -n "tab:supp_network_impl" supplementary.tex
```

Similarly, figure files can be located from the LaTeX source by searching for the corresponding image name or figure label, for example:

```bash
grep -n "poisson_dagar_recovery.png" main.tex
grep -n "boundary_probability_quality.png" supplementary.tex
grep -n "glasgow_boundary_agreement.png" main.tex supplementary.tex
```

The repository is organized so that the manuscript and Supplementary Materials can be checked against the exact output files used for all reported figures and tables.
