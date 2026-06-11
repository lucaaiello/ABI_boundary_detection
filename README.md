## Reproducing the manuscript and supplementary results

This repository contains the code used to reproduce the results reported in the manuscript and Supplementary Materials for *Amortized Bayesian Boundary Detection on Varying-Size Graphs*. The workflow is organized around three main blocks:

1. **Simulation experiments**, used to validate ABI-DAGAR under the proposed generative model. These include parameter recovery, posterior calibration, boundary-probability diagnostics, posterior predictive checks, the simulation benchmark against model-matched MCMC-DAGAR, computational summaries, and the ablation study of summary statistics.

2. **Real-data analyses**, where the trained ABI-DAGAR approximation and the localized `CARBayes` benchmark are applied to the Glasgow respiratory disease and California lung cancer datasets and then compared at the parameter, boundary-probability, boundary-set, fitted-risk, and runtime levels.

3. **Supplementary sensitivity and diagnostic analyses**, including additional simulation diagnostics, additional real-data diagnostics, the real-data model-matched MCMC-DAGAR comparison, and the ABI-CAR sensitivity analysis.

The final figures used in the manuscript and Supplementary Materials are stored in the `Images/` folder. The numerical values reported in tables are contained in the corresponding LaTeX source files, `main.tex` and `supplementary.tex`, and are produced by the simulation, benchmark, posterior predictive, ablation, and real-data analysis code in this repository.

### Recommended order of execution

The analyses should be run in the following order.

### 1. Train the ABI-DAGAR amortized posterior approximation

The first step is to train the ABI-DAGAR amortized posterior approximation. This step simulates training datasets from the Poisson-DAGAR boundary-detection model and trains the SetTransformer summary network together with the conditional normalizing flow.

Main output:

| Output                    | Used in                |
| ------------------------- | ---------------------- |
| `Images/history_plot.png` | Main manuscript Fig. 1 |

This trained ABI-DAGAR approximation is then reused in the simulation validation and in the Glasgow and California real-data analyses.

### 2. Run the ABI-DAGAR simulation experiments

The simulation experiments evaluate the trained ABI-DAGAR approximation on held-out simulated datasets generated from the proposed Poisson-DAGAR boundary-detection model. This block contains parameter recovery, posterior calibration, posterior boundary probabilities, decision-rule diagnostics, posterior predictive checks, and computational summaries.

Main manuscript outputs:

| Output                                            | Used in                 |
| ------------------------------------------------- | ----------------------- |
| `Images/poisson_dagar_recovery.png`               | Main manuscript Fig. 2  |
| `Images/poisson_dagar_calibration_histogram.png`  | Main manuscript Fig. 3  |
| `Images/boundary_sens_spec_curve.png`             | Main manuscript Fig. 4  |
| `Images/boundary_mpm_histograms.png`              | Main manuscript Fig. 5  |
| Table values for `tab:sim_recovery` in `main.tex` | Main manuscript Table 1 |

Supplementary outputs:

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

### 3. Run the simulation benchmark against MCMC-DAGAR

The simulation benchmark compares ABI-DAGAR with a model-matched MCMC implementation of the same thresholded Poisson-DAGAR model, denoted MCMC-DAGAR. This is still part of the simulation validation block, because both methods are applied to held-out simulated datasets with known generating values.

Main manuscript outputs:

| Output                                            | Used in                |
| ------------------------------------------------- | ---------------------- |
| `Images/parameter_recovery_bars.png`              | Main manuscript Fig. 6 |
| `Images/parameter_recovery_truth_scatter.png`     | Main manuscript Fig. 7 |
| `Images/parameter_recovery_agreement_scatter.png` | Main manuscript Fig. 7 |
| `Images/parameter_bias_interval_boxplots.png`     | Main manuscript Fig. 7 |

Supplementary outputs:

| Output                            | Used in            |
| --------------------------------- | ------------------ |
| `Images/boundary_metric_bars.png` | Supplement Fig. S7 |
| `Images/runtime_comparison.png`   | Supplement Fig. S7 |

### 4. Run the ablation study of summary statistics

The ablation study is part of the simulation experiments. It evaluates whether the engineered graph-aware summaries provide useful information for posterior inference. For each ablated representation, a separate amortized posterior approximator is retrained from scratch using the same architecture and training protocol. The resulting models are evaluated on the same validation datasets and compared with the full-summary baseline.

Outputs:

| Output                                                    | Used in                                                                             |
| --------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `Images/ablation_delta_vs_baseline_heatmaps.png`          | Main manuscript ablation figure, if the ablation study is included in the main text |
| `Images/ablation_boundary_delta_vs_baseline_heatmaps.png` | Main manuscript ablation figure, if the ablation study is included in the main text |
| `Images/ablation_recovery_heatmaps.png`                   | Supplement Fig. S8                                                                  |
| `Images/ablation_boundary_heatmaps.png`                   | Supplement Fig. S9                                                                  |
| `Images/ablation_training_histories.png`                  | Supplement Fig. S10                                                                 |

If the ablation study is placed entirely in the Supplementary Materials, the first two files above may be omitted from the final manuscript build or kept only as auxiliary output files.

### 5. Run the real-data ABI-DAGAR and `CARBayes` analyses

The real-data analyses apply the trained ABI-DAGAR approximation to the Glasgow respiratory disease and California lung cancer datasets. The same datasets are also analyzed using the localized smoothing model implemented in `CARBayes`. The main real-data comparison is therefore between ABI-DAGAR and `CARBayes`.

This block produces posterior summaries, selected boundary sets, fitted-risk comparisons, boundary-probability comparisons, and runtime summaries for the Glasgow and California applications.

Main manuscript outputs:

| Output                                                             | Used in                                             |
| ------------------------------------------------------------------ | --------------------------------------------------- |
| `Images/glasgow_boundary_agreement.png`                            | Main manuscript real-data boundary agreement figure |
| `Images/california_boundary_agreement.png`                         | Main manuscript real-data boundary agreement figure |
| Table values for `tab:carbayes-parameter-comparison` in `main.tex` | Main manuscript Table 2                             |
| Table values for `tab:carbayes-agreement` in `main.tex`            | Main manuscript Table 3                             |

Supplementary outputs:

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

### 6. Run the real-data MCMC-DAGAR comparison

The Supplementary Materials also contain a real-data comparison between ABI-DAGAR and a model-matched MCMC-DAGAR implementation. This is separate from the main `CARBayes` real-data comparison. Its purpose is to distinguish model-class differences from amortization error: `CARBayes` is an external localized-CAR benchmark, whereas MCMC-DAGAR targets the same DAGAR boundary model as ABI-DAGAR.

Supplementary outputs:

| Output                                                                               | Used in             |
| ------------------------------------------------------------------------------------ | ------------------- |
| `Images/glasgow_boundary_agreement_dagarbayes.png`                                   | Supplement Fig. S14 |
| `Images/california_boundary_agreement_dagarbayes.png`                                | Supplement Fig. S14 |
| Table values for ABI-DAGAR and MCMC-DAGAR posterior summaries in `supplementary.tex` | Supplement Table S6 |
| Table values for `tab:mcmc_dagar_agreement` in `supplementary.tex`                   | Supplement Table S7 |

If desired, the files containing `dagarbayes` in the name can be renamed to `mcmc_dagar` for consistency with the terminology used in the manuscript, for example:

```text
Images/glasgow_boundary_agreement_mcmc_dagar.png
Images/california_boundary_agreement_mcmc_dagar.png
```

The corresponding file names should then also be updated in `supplementary.tex`.

### 7. Run the ABI-CAR sensitivity analysis

The ABI-CAR analysis is an additional sensitivity experiment reported in the Supplementary Materials. It repeats the amortized Bayesian workflow under a localized Leroux CAR prior matching the spatial-dependence specification used by Lee and Mitchell (2012) and the `CARBayes` benchmark. This analysis serves as a positive-control experiment: when the amortized model is aligned with the localized CAR prior, the resulting boundary conclusions closely reproduce the `CARBayes` benchmark.

Supplementary outputs:

| Output                                                                             | Used in              |
| ---------------------------------------------------------------------------------- | -------------------- |
| `Images/poisson_car_recovery.png`                                                  | Supplement Fig. S15  |
| `Images/parameter_coverage_curve_car.png`                                          | Supplement Fig. S15  |
| `Images/glasgow_boundary_agreement_car.png`                                        | Supplement Fig. S16  |
| `Images/california_boundary_agreement_car.png`                                     | Supplement Fig. S16  |
| Table values for `tab:supp_abicar_recovery` in `supplementary.tex`                 | Supplement Table S10 |
| Table values for ABI-CAR and `CARBayes` posterior summaries in `supplementary.tex` | Supplement Table S11 |
| Table values for ABI-CAR boundary and fitted-risk agreement in `supplementary.tex` | Supplement Table S12 |

### 8. Compile the manuscript and Supplementary Materials

After the figures and numerical summaries have been generated, compile `main.tex` and `supplementary.tex`. The manuscript and supplement read the final figures from the `Images/` folder. The table values are included directly in the LaTeX source files.

Because the main manuscript and supplement use cross-references between files, compile both documents enough times for cross-references to resolve correctly.

### Output index

The following index maps the manuscript and supplementary outputs to the corresponding figure files and table locations.

#### Main manuscript figures

| Manuscript item                     | File(s) in repository                                                                                                                           | Contents                                                                                      |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| Fig. 1                              | `Images/history_plot.png`                                                                                                                       | Training loss over 100 epochs for the amortized posterior approximation.                      |
| Fig. 2                              | `Images/poisson_dagar_recovery.png`                                                                                                             | Parameter recovery on 200 held-out simulated datasets.                                        |
| Fig. 3                              | `Images/poisson_dagar_calibration_histogram.png`                                                                                                | Simulation-based calibration rank histograms.                                                 |
| Fig. 4                              | `Images/boundary_sens_spec_curve.png`                                                                                                           | Ranking diagnostic for posterior boundary probabilities.                                      |
| Fig. 5                              | `Images/boundary_mpm_histograms.png`                                                                                                            | Sensitivity and specificity under the median-probability rule.                                |
| Fig. 6                              | `Images/parameter_recovery_bars.png`                                                                                                            | ABI-DAGAR versus MCMC-DAGAR simulation benchmark: mean absolute error and empirical coverage. |
| Fig. 7                              | `Images/parameter_recovery_truth_scatter.png`; `Images/parameter_recovery_agreement_scatter.png`; `Images/parameter_bias_interval_boxplots.png` | Detailed ABI-DAGAR versus MCMC-DAGAR parameter-level comparison.                              |
| Real-data boundary agreement figure | `Images/glasgow_boundary_agreement.png`; `Images/california_boundary_agreement.png`                                                             | Boundary agreement between ABI-DAGAR and `CARBayes`.                                          |

#### Main manuscript tables

| Manuscript item | Location in source                                    | Contents                                                                                |
| --------------- | ----------------------------------------------------- | --------------------------------------------------------------------------------------- |
| Table 1         | `main.tex`, label `tab:sim_recovery`                  | Simulation parameter recovery: bias, RMSE, correlation, (R^2), and empirical coverage.  |
| Table 2         | `main.tex`, label `tab:carbayes-parameter-comparison` | Posterior medians and 95% credible intervals for ABI-DAGAR and `CARBayes`.              |
| Table 3         | `main.tex`, label `tab:carbayes-agreement`            | Real-data agreement between ABI-DAGAR and `CARBayes` under the median-probability rule. |

#### Supplementary figures

| Supplement item | File(s) in repository                                                                                             | Contents                                                                              |
| --------------- | ----------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| Fig. S1         | `Images/parameter_coverage_curve.png`                                                                             | Empirical coverage versus nominal coverage.                                           |
| Fig. S2         | `Images/poisson_dagar_calibration_ecdf.png`                                                                       | Simulation-based calibration ECDF-difference plots.                                   |
| Fig. S3         | `Images/graph_size_posterior_zscores.png`                                                                         | Posterior (z)-scores stratified by graph-size bin.                                    |
| Fig. S4         | `Images/regime_specific_parameter_errors.png`                                                                     | Mean absolute recovery error by true (\rho), true (\eta), graph size, and edge count. |
| Fig. S5         | `Images/boundary_probability_quality.png`                                                                         | Boundary-probability reliability and dissimilarity-profile diagnostics.               |
| Fig. S6         | `Images/predictive_checks.png`                                                                                    | Posterior predictive diagnostics on held-out simulated datasets.                      |
| Fig. S7         | `Images/boundary_metric_bars.png`; `Images/runtime_comparison.png`                                                | Additional simulation MCMC-DAGAR benchmark diagnostics and runtime comparison.        |
| Fig. S8         | `Images/ablation_recovery_heatmaps.png`                                                                           | Parameter-specific recovery metrics for the ablation study.                           |
| Fig. S9         | `Images/ablation_boundary_heatmaps.png`                                                                           | Boundary-detection diagnostics across summary representations.                        |
| Fig. S10        | `Images/ablation_training_histories.png`                                                                          | Training and validation loss trajectories for the ablation study.                     |
| Fig. S11        | `Images/glasgow_post_pred_check.png`; `Images/california_post_pred_check.png`                                     | Posterior predictive checks for the real-data applications.                           |
| Fig. S12        | `Images/glasgow_risk_comparison.png`; `Images/california_risk_comparison.png`                                     | Fitted-risk comparison between ABI-DAGAR and `CARBayes`.                              |
| Fig. S13        | `Images/glasgow_edge_probability_vs_dissimilarity.png`; `Images/california_edge_probability_vs_dissimilarity.png` | Posterior boundary probability versus standardized edge dissimilarity.                |
| Fig. S14        | `Images/glasgow_boundary_agreement_dagarbayes.png`; `Images/california_boundary_agreement_dagarbayes.png`         | Real-data boundary agreement between ABI-DAGAR and MCMC-DAGAR.                        |
| Fig. S15        | `Images/poisson_car_recovery.png`; `Images/parameter_coverage_curve_car.png`                                      | ABI-CAR simulation diagnostics.                                                       |
| Fig. S16        | `Images/glasgow_boundary_agreement_car.png`; `Images/california_boundary_agreement_car.png`                       | Boundary agreement between ABI-CAR and `CARBayes`.                                    |

#### Supplementary tables

| Supplement item | Location in source                                      | Contents                                                                                       |
| --------------- | ------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| Table S1        | `supplementary.tex`, label `tab:sim_boundary_sup`       | Additional simulation boundary-probability summaries.                                          |
| Table S2        | `supplementary.tex`, label `tab:sim_ppc_sup`            | Additional simulation posterior predictive summaries.                                          |
| Table S3        | `supplementary.tex`, label `tab:sim_computation`        | Simulation computational details.                                                              |
| Table S4        | `supplementary.tex`, label `tab:realdata_risk_extra`    | Additional real-data fitted-risk comparison metrics between ABI-DAGAR and `CARBayes`.          |
| Table S5        | `supplementary.tex`, label `tab:realdata_edge_extra`    | Additional real-data boundary-probability comparison metrics between ABI-DAGAR and `CARBayes`. |
| Table S6        | `supplementary.tex`, real-data MCMC-DAGAR subsection    | Posterior medians and 95% credible intervals for ABI-DAGAR and MCMC-DAGAR.                     |
| Table S7        | `supplementary.tex`, label `tab:mcmc_dagar_agreement`   | Real-data agreement between ABI-DAGAR and MCMC-DAGAR.                                          |
| Table S8        | `supplementary.tex`, label `tab:realdata_runtime_extra` | Real-data posterior sampling runtimes for ABI-DAGAR and `CARBayes`.                            |
| Table S9        | `supplementary.tex`, label `tab:supp_network_impl`      | Neural-network and training configuration for ABI-DAGAR.                                       |
| Table S10       | `supplementary.tex`, label `tab:supp_abicar_recovery`   | ABI-CAR parameter recovery on held-out simulated datasets.                                     |
| Table S11       | `supplementary.tex`, ABI-CAR real-data subsection       | Posterior medians and 95% credible intervals for ABI-CAR and `CARBayes`.                       |
| Table S12       | `supplementary.tex`, ABI-CAR real-data subsection       | Boundary and fitted-risk agreement between ABI-CAR and `CARBayes`.                             |

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
