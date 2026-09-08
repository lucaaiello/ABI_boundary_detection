# South Korea application data

This directory contains the reproducible input preparation for the South Korean
real-data application of the shared ABI-DAGAR network.

## Primary application

- Spatial units: 241 municipalities in the largest queen-contiguous component of
  the 2019 South Korean municipal boundary layer.
- Outcome: deaths from malignant neoplasms of the trachea, bronchus, and lung
  (ICD-10 C33-C34), pooled over 2015-2019.
- Expected counts: indirect age-sex standardization using annual national
  cause-specific rates and annual municipal midyear populations.
- Boundary covariate: age-standardized current smoking prevalence, averaged over
  2015-2019 and standardized across the 241 retained municipalities.

The processed GeoPackage also retains suicide mortality, circulatory-disease
mortality, and National Basic Livelihood Security recipient prevalence as
documented alternative variables. The reported application uses lung-cancer
mortality and current-smoking prevalence.

## Public sources

- Municipal cause-specific mortality: KOSIS table `DT_1B34E13`.
  <https://kosis.kr/statisticsList/mass/mass_list_e.jsp?org_id=101&tbl_id=DT_1B34E13&vw_cd=MT_ZTITLE>
- Municipal midyear population by sex and five-year age: KOSIS table
  `DT_1B040M5`.
  <https://kosis.kr/statisticsList/mass/mass_list_e.jsp?list_id=A_7&org_id=101&tbl_id=DT_1B040M5>
- Cause-specific mortality by province, sex, and age: KOSIS table `DT_1B34E11`.
  <https://kosis.kr/statisticsList/mass/mass_list_e.jsp?org_id=101&tbl_id=DT_1B34E11>
- Korean Community Health Status Indicators database, version 1.7: KDCA.
  <https://chs.kdca.go.kr/chs/recsRoom/dataBaseMain.do>
- 2019 municipal boundary archive (`TL_SCCO_SIG`): GIS Developer mirror of the
  Korean administrative boundary file.
  <http://www.gisdeveloper.co.kr/download/admin_shp/SIG_201905.zip>

## Reproduction

Run the following commands from the repository root:

```powershell
python "Real Data Analysis/Data/South_Korea/download_south_korea_data.py"
python "Real Data Analysis/Data/South_Korea/prepare_south_korea_data.py"
```

The preparation script validates all required cells and stops on unmatched
geographies. It explicitly handles the 2018 rename of Incheon Nam-gu to
Michuhol-gu and propagates parent-city KDCA values to constituent non-autonomous
districts only when the source database reports the indicator at parent-city
level.

## Main outputs

- `mortality_data_south_korea.gpkg`: analysis-ready 241-area GeoPackage.
- `adjacency_matrix_south_korea.csv`: 241 by 241 queen-contiguity matrix.
- `south_korea_kosis_crosswalk.csv`: boundary-to-KOSIS geography crosswalk.
- `excluded_island_municipalities.csv`: nine units outside the largest connected
  component.
- `south_korea_data_manifest.json`: topology, covariate, and outcome checks.
- `raw/checksums.csv`: checksums for the downloaded archives.
