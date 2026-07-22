# Data

Raw spectral and reference soil analysis data are **not stored in this GitHub
repository**. They are archived on Zenodo (see badge in the top-level
`README.md`) to keep the code repository lightweight.

Download the archived dataset from Zenodo and place the files here before
running `scripts/01_calibration_pipeline_v3.R`:

| File | Description |
|---|---|
| `asd_vnir_raw_all.csv` | Raw ASD FieldSpec VNIR-NIR reflectance spectra (350–2500 nm), one row per sample |
| `isc_nir_raw_all.csv` | Raw NIRVascan (Inno-Spectra NIR-S-G1) reflectance spectra (900–1700 nm), one row per sample |
| `db_soil_analysis_all.csv` | Reference (wet-chemistry) soil analysis results: Total Carbon (%), Total Nitrogen (%), pH (H2O), plant-available P (Double Lactate) |

Zenodo DOI: `10.5281/zenodo.XXXXXXX` (placeholder — update after first deposit)
