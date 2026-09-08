# ABI-DAGAR versus MCMC-DAGAR benchmark

This folder implements the descriptive simulation benchmark between the trained
ABI-DAGAR posterior approximator and a closely aligned Metropolis-within-Gibbs
DAGAR sampler. Both methods operate on the same fixed bank of 100 simulated
datasets. The likelihood, DAGAR construction, threshold mechanism, and identity
ordering agree; prior differences are explicit in the R runner and manuscript
supplement.

## Workflow

Run the components in this order from the repository root:

1. `python "Simulation Experiments/ABI_vs_MCMC/export_benchmark_datasets.py" --num-datasets 100 --seed 123`
2. Run `Simulation Experiments/ABI_vs_MCMC/run_abi_benchmark.ipynb` top to bottom.
3. Run `Rscript "Simulation Experiments/ABI_vs_MCMC/run_mcmc_benchmark.R"` with the production settings recorded in `mcmc_config.csv`.
4. Run `Simulation Experiments/ABI_vs_MCMC/compare_abi_vs_mcmc.ipynb` top to bottom.

The submitted result bank is:

`datasets/benchmark_bank_seed123_n100/`

## Files

| File | Role |
| --- | --- |
| `benchmark_utils.py` | Shared simulator and serialization helpers |
| `export_benchmark_datasets.py` | Creates the fixed synthetic benchmark bank |
| `dagar_poisson_boundary_mwg.cpp` | Rcpp Metropolis-within-Gibbs sampler |
| `run_abi_benchmark.ipynb` | Applies the trained ABI checkpoint to the bank |
| `run_mcmc_benchmark.R` | Runs the MCMC-DAGAR comparator and computes ESS/MCSE diagnostics |
| `compare_abi_vs_mcmc.ipynb` | Produces reported comparison tables and PNG figures |

## Retained results

The repository retains all 100 synthetic `.npz` inputs, combined parameter,
boundary, timing, acceptance, and chain-diagnostic CSVs, and per-dataset edge
probabilities required by the comparison notebook. Large posterior-draw files,
R-friendly duplicate input bundles, and other per-dataset intermediates are not
versioned because the supplied runners recreate them.

Reported tables are under `comparison_abi_vs_mcmc/`; reported plots are under
`comparison_abi_vs_mcmc/plots/`.
