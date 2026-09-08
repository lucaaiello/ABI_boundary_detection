# Standalone manuscript figures

These scripts generate manuscript figures that are not direct outputs of an
analysis notebook. Run all three from the repository root with:

```powershell
python "Figure Generation/generate_standalone_manuscript_figures.py"
```

The generators are deterministic and write PNG files directly to
`Figure Generation/Images/`.

| Generator | Output | Data dependencies |
|---|---|---|
| `generate_abi_workflow.py` | `abi_train_validate_deploy_workflow.png` | None; conceptual diagram |
| `generate_boundary_probability_residual_contrast.py` | `boundary_probability_residual_contrast.png` | `results_ABI_vs_CARBayes/*/edge_comparison.csv` and `results_exploratory_analysis/*_harmonized_data.csv` |
| `generate_deployment_support_training_distribution.py` | `deployment_support_training_distribution.png` | Original training input design, empirical harmonized data, and adjacency matrices |

The deployment-support generator uses seed `20260907` and saves its complete
3,000-configuration metric bank, together with the empirical markers, to
`Real Data Analysis/results_exploratory_analysis/deployment_support_metrics_seed20260907_n3000.csv`.

After running the exploratory and real-data notebooks, standardize and stage
the reported map panels with:

```powershell
python "Figure Generation/refresh_manuscript_map_assets.py"
```

This preserves each map's aspect ratio, sets the color bar to 75% of the map
height, and copies the final PNG panels to `Figure Generation/Images/`.
