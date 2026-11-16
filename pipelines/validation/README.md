# External validation of EII v2

This suite compares EII v2 and its three components with independent global layers and
with the published EII of UNEP-WCMC (Hill et al. 2026, *Ecological Indicators* 190, 115314).
It reports numbers only: correlations, area summaries and agreement statistics, with no
pass/fail thresholds. The results and their interpretation are on the documentation page
[Validation](https://landler-io.github.io/ecosystem-integrity-index/methodology/validation/).

## What is validated

1. **Cross-indicator correlations.** The cross-indicator validation in the supplement of
   Hill et al. (2026, section A2) is repeated for our product, with the same layers, the same
   Eckert IV (ESRI:54012) 1 km grid and 500,000 equal-area random land points (seed 42):
   - EII against the human footprint (WCS Human Influence Index 2013), the Biodiversity
     Habitat Index (BHI v2, habitat and species, 2015) and the Mean Species Abundance
     (GLOBIO4, 2015), Pearson r and Spearman rho with p-values (their Fig A4);
   - correlations between the three components (Fig A3);
   - ecoregion means against Beyer et al. (2020) ecoregion intactness (Fig A5);
   - the distribution of EII and of structural integrity over human-influence classes from
     anthromes, the human footprint and Low Impact Areas, and over the three Locke et al.
     (2019) conditions, with the probability that a least-influenced km² scores higher than a
     most-influenced one (AUC; Fig A6);
   - functional integrity against nature's contributions to people (Chaplin-Kramer et al.
     2019: pollination, nutrient retention, nitrogen retention, NDR) and soil organic carbon
     (GSOCmap; Fig A8);
   - compositional integrity against BHI, MSA and the defaunation index (Benítez-López et
     al. 2019; Fig A9).
2. **Product comparison with the WCMC EII.** Our EII for 2021, the year nearest to the
   latest WCMC map, against the published `EII_2020_1km.tif` on the WCMC grid (EPSG:4326,
   ~1 km), area-weighted: Pearson r, bias (ours minus WCMC), MAE, RMSE and difference
   quantiles on the shared land mask, per WWF biome, per RESOLVE ecoregion and per Sayre et
   al. (2020) World Terrestrial Ecosystem. Our 2025 map is compared in the same way.

The paper's own values are stored next to ours in the tables (`paper_r`).

## Inputs

All bulk data lives under the data root (`EII_VALIDATION_DATA_ROOT`, default
`./validation_data`). External layers keep their own licences (several are
non-commercial); they are not redistributed with this repository. The tables in
[`outputs/`](outputs) contain derived per-unit means of two of them: the UNEP-WCMC EII
(Hill et al. 2026; CC BY 4.0) per biome, ecoregion and World Terrestrial Ecosystem, and the
ecoregion intactness of Beyer et al. (2020; CC BY 4.0) per ecoregion.

| Layer | Source | Path under `inputs/layers/` |
|---|---|---|
| WCMC EII V1 (2000-2020) | `run_validation.py fetch-wcmc` (public S3 `EII_V1.zip`, CC BY 4.0) | `../wcmc/EII_<year>_1km.tif` |
| Human footprint 2013 (HII) | [WCS Human Footprint](https://www.wcshumanfootprint.org/data-access) | `hfi_2013/hii_2013-01-01.tif` |
| Anthromes v2, 2000 | [Harvard Dataverse doi:10.7910/DVN/6KKAWS](https://doi.org/10.7910/DVN/6KKAWS) | `anthromes/anthro2_a2000.tif` |
| Low Impact Areas | [Zenodo 4943288](https://zenodo.org/records/4943288) (CC0) | `low_impact_areas/{Very_Low,Low}_Impact.tif` |
| BHI v2 habitat / species | [CSIRO doi:10.25919/tt2t-h452](https://doi.org/10.25919/tt2t-h452) | `bhi/data/BHI_v2_Proportion_of_{Habitat,Species}/BILBI_P_BHIv2_*_2015.tif` |
| GLOBIO4 MSA 2015 | [globio.info](https://www.globio.info/globio-data-downloads) (CC BY 4.0) | `msa/TerrestrialMSA_2015_World.tif` |
| Ecoregion intactness 2013 | [Zenodo 15850563](https://zenodo.org/records/15850563) (CC BY 4.0) | `ecoregion_intactness/Intactness_2013.tif` |
| NCP layers (Chaplin-Kramer et al. 2019) | [Zenodo 21249617](https://zenodo.org/records/21249617) and linked parts | `chaplin-kramer/extracted/*.tif` |
| GSOCmap v1.5.0 | [FAO GSOCmap](https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/global-soil-organic-carbon-map-gsocmap/en/) | `gsoc/GSOCSEQ.GSOCMAP1-5-0.tif` |
| Defaunation index | [figshare 7012598](https://figshare.com/articles/dataset/Defaunation_index_maps/7012598) (CC BY 4.0) | `defaunation/DefInd.tif` |
| Copernicus LC100 2015 crop and built-up cover | [Zenodo 3939038](https://zenodo.org/records/3939038) | `cgls_lc100_2015/*CoverFraction-layer_EPSG-4326.tif` |
| World Terrestrial Ecosystems 2020 | [USGS doi:10.5066/P9DO61LP](https://doi.org/10.5066/P9DO61LP) | `sayre_wte/unpacked/commondata/raster_data/WorldEcosystem.tif` |

A raster of RESOLVE Ecoregions 2017 ecoregion ids (`EII_VALIDATION_ECOREGIONS`, default
`<data root>/inputs/layers/ecoregions/resolve_ecoregions.tif`) provides ecoregion and biome
ids; `pipelines/preprocessing/reference/eco_id_to_biome.json` maps ecoregions to WWF biomes. Exact file names, scale factors and
nodata handling are in `E_LAYERS` / `W_LAYERS` of `run_validation.py`. Two layers of the
paper are not included: the Red List of Ecosystems for South Africa (not obtainable) and
the per-taxon defaunation maps.

Our product is read from `EII_VALIDATION_EII_ROOT` (default `./eii_v2`), the tiles written
by the production run (`<year>/<tile>.tif`, bands EII, functional, structural, compositional
first), i.e. the GeoTIFFs behind the Earth Engine collection
`projects/landler-open-data/assets/products/eii/v2_0/release/global`.

## Running

From the repository root, with the `[pipelines]` extra installed
(`pip install -e ".[pipelines]"`; the point sampling needs `pyproj`):

```bash
export EII_VALIDATION_DATA_ROOT=/path/to/validation_data
export EII_VALIDATION_EII_ROOT=/path/to/eii_v2
export EII_VALIDATION_ECOREGIONS=/path/to/resolve_ecoregions.tif
python pipelines/validation/run_validation.py fetch-wcmc    # once, ~2.4 GB download
python pipelines/validation/run_validation.py list          # warp jobs, done / todo
python pipelines/validation/run_validation.py all --workers 8
```

`all` runs four stages, which can also be called one by one:

| Stage | Does | Time (12 workers) |
|---|---|---|
| `warp` | external layers to the Eckert IV 1 km grid (once, shared by runs); our EII, F, S, C per year to Eckert IV and our EII to the WCMC grid | ~2 min per year; external layers once, about an hour |
| `points` | 500,000 equal-area land points sampled from every Eckert IV raster | ~2 min |
| `analyse` | correlation, area and WCMC comparison tables | ~5 min |
| `figures` | PNG figures from the tables | ~1 min |

Warping is memory-bounded and parallel (`warping.py`, rasterio only, no GDAL command-line
tools needed); finished rasters are skipped on rerun.

## Outputs

Per run, under `EII_VALIDATION_RUN_DIR` (default `<data root>/runs/<name of the EII root>`):

```
our_eii/        VRT mosaics of the product tiles per year
e_grid/         our EII, F, S, C per year on the Eckert IV 1 km grid
w_grid/         our EII per year on the WCMC grid
points_e.csv.gz the 500k point table (all layers and our bands)
outputs/        run_manifest.json, tables/*.csv, figures/*.png
validation.log
```

A copy of the summary tables and figures of the current product is kept in
[`outputs/`](outputs) in this directory. Main tables:

| Table | Content |
|---|---|
| `a4_pairs.csv`, `a4_mean_abs_r.csv` | EII (ours and WCMC) against human footprint, BHI and MSA; mean absolute r |
| `a3_components.csv` | correlations between components and share of the limiting component |
| `a5_ecoregion_corr.csv` | ecoregion means against ecoregion intactness |
| `a6_summary.csv` | EII and S per human-influence class and Locke condition, AUC |
| `a8_functional.csv`, `a9_compositional.csv` | F against NCP and SOC; C against BHI, MSA and retained fauna (defaunation index, both orientations) |
| `wcmc_global.csv` | r, bias, MAE, RMSE, difference quantiles and mask overlap against WCMC 2020 |
| `wcmc_biome.csv`, `wcmc_ecoregion.csv`, `wcmc_sayre.csv` | the same per biome, ecoregion, ecosystem |
| `wcmc_biome_own_mask.csv`, `wcmc_sayre_summary.csv` | biome and ecosystem means on each product's own mask |
| `wcmc_points_corr.csv` | Pearson and Spearman against WCMC on the equal-area points |

## Code

| Module | Role |
|---|---|
| `run_validation.py` | CLI, paths, layer registry, warp jobs |
| `warping.py` | parallel blockwise warping onto both grids |
| `grid.py` | the Eckert IV 1 km grid |
| `analysis.py` | point sampling, correlations, area summaries, WCMC comparison |
| `metrics.py` | Pearson/Spearman with p-values, Locke et al. conditions |
| `figures.py` | figures |
| `report_figures.py` | the three figures of the methodology report (`docs/report/figures/`; its land-cover sample is cached there as `landcover_sample_<year>.csv.gz`) |
| `wcmc.py` | download and checksum verification of the WCMC EII |

Offline tests: `tests/test_pipelines/test_validation.py`.
