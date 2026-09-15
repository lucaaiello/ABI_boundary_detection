# ABI-DAGAR versus prior-matched MCMC-DAGAR benchmark

This folder implements a paired comparison between the trained ABI-DAGAR
posterior approximator and a single-chain MCMC sampler targeting the simulator's
posterior. The two methods use the same fixed bank of 100 held-out datasets, the
same Poisson likelihood, DAGAR construction, deterministic isolate repair,
latent-field centering, boundary threshold, identity ordering, and scalar
priors.

The fixed bank of 100 held-out datasets is included in this folder. The
preparation script converts that bank into the ignored R-readable matrix
bundles needed by the MCMC runner; it can also import a compatible external
bank through `--source-dir`.

## Matched model

The scalar priors are

- `beta0 ~ Normal(0, 0.5^2)`;
- `sigma2_w ~ HalfNormal(0.5)`, with the half-normal placed on `sigma2_w`;
- `eta_raw ~ Uniform(0, 1)` and `eta = M * eta_raw`;
- `rho ~ Uniform(0, 1)`.

For each proposed `eta`, the filtered graph is thresholded and repaired using
the simulator's deterministic isolate rule. The latent field is represented as

`u ~ Normal(0, sigma2_w * Q^{-1})`, `w = u - mean(u)`,

and `w` enters the Poisson log mean. This is the generative construction used by
the ABI training simulator.

## Workflow

Run the components from the repository root in this order:

1. `python "Simulation Experiments/ABI_vs_matched_MCMC/prepare_benchmark_bank.py"`
2. Execute `Simulation Experiments/ABI_vs_matched_MCMC/run_abi_benchmark.ipynb`.
3. `Rscript "Simulation Experiments/ABI_vs_matched_MCMC/run_matched_mcmc_benchmark.R"`
4. Execute `Simulation Experiments/ABI_vs_matched_MCMC/compare_abi_vs_matched_mcmc.ipynb`.

Both inference runners default to all 100 datasets and 10,000 retained draws.
The R runner uses one chain by design. For a smoke test, pass `--max-datasets 1
--n-iter 1000 --burnin 500 --n-adapt 500` to the R runner and set
`ABI_MATCHED_MAX_DATASETS=1` and `ABI_MATCHED_NUM_SAMPLES=500` in the
environment before executing the ABI notebook. The comparison notebook accepts
`ABI_MATCHED_ABI_RESULTS_DIR`, `ABI_MATCHED_MCMC_RESULTS_DIR`, and
`ABI_MATCHED_COMPARISON_DIR` overrides for testing alternate result folders.

## Files

| File | Role |
| --- | --- |
| `prepare_benchmark_bank.py` | Copies the fixed held-out inputs and rebuilds R-readable bundles |
| `matched_abi_runner.py` | Shared implementation used by the ABI notebook |
| `dagar_poisson_boundary_matched_mwg.cpp` | Prior- and construction-matched Rcpp sampler |
| `run_matched_mcmc_benchmark.R` | Runs the single-chain matched MCMC benchmark |
| `run_abi_benchmark.ipynb` | Applies the trained ABI checkpoint and saves posterior draws |
| `matched_comparison.py` | Comparison, MCSE-controlled error decomposition, table, and plotting helpers |
| `compare_abi_vs_matched_mcmc.ipynb` | Compares error components, marginal posteriors, recovery, boundaries, and runtime |

Generated data and results are placed below
`datasets/benchmark_bank_seed123_n100/`. The versioned reference outputs
include the summaries and figures used in the manuscript. Regenerable R input
bundles and gzip-compressed posterior draws are intentionally excluded from
version control.
