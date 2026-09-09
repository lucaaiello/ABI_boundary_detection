# Data access and redistribution

The repository contains the analysis-ready areal inputs used by the three case
studies. No individual-level health data are used by the analyses.

The repository's MIT License applies only to original code developed for this
project. It does not apply to third-party datasets or override their respective
access, licensing, citation, or redistribution terms.

## Greater Glasgow

The respiratory-disease data are distributed through the R package
`CARBayesdata` as `respiratorydata`. The repository includes the analysis-ready
polygon file and adjacency matrix used by the notebooks. See the package
documentation for the original variables and geographic metadata.

## California

The county-level lung-cancer incidence data originate from SEER*Stat. The
analysis-ready GeoPackage is included at:

`Real Data Analysis/Data/respiratory_data_california.gpkg`

It contains one row per California county and the following fields:

| Field | Description |
| --- | --- |
| `county` | County identifier used by the adjacency matrix |
| `lung_O_count` | Observed lung-cancer incidence count |
| `lung_E_count` | Expected lung-cancer incidence count |
| `smoking` | Adult smoking prevalence |
| `lung_standard_ratio` | Observed-to-expected incidence ratio |
| `geometry` | County polygon geometry |

The county adjacency matrix, analysis code, posterior summaries, edge-level
outputs, aggregate comparison metrics, and reported figures are also included.
The original source and current access conditions are documented by the
[official SEER data page](https://seer.cancer.gov/data/access.html).

## South Korea

The South Korean input pipeline uses public KOSIS mortality and population
tables, the KDCA Community Health Status Indicators workbook, and the May 2019
municipal boundary layer. Run:

```powershell
python "Real Data Analysis/Data/South_Korea/download_south_korea_data.py"
python "Real Data Analysis/Data/South_Korea/prepare_south_korea_data.py"
```

The folder-level README documents the source tables, identifier crosswalk, graph
restriction, and generated files. SHA-256 checksums are retained for all source
archives.

## Responsible use

Users are responsible for complying with each source provider's current access,
citation, and redistribution terms.
