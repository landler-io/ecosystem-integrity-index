# Local Modulation (Natural Capital)

**Local Modulation** adjusts the base Ecosystem Integrity Index (EII) using plot-level **Natural Capital (NC)** performance. It adds three dimensions, namely biodiversity, soil, and water as KPIs, combines them into an NC score, and shifts EII by up to ±0.05. This keeps a shared structure (the same three dimensions everywhere) while allowing each ecosystem to use the most relevant KPIs within those dimensions, for example, habitat fraction for biodiversity in forests or water table depth for peatlands. In the follwing we describe the approach, follwed up by an implemented example of how it works in practice.

## Modulation Formula

The modulated EII is:

$$ EII_{mod} = EII + (NC - 0.5) \times R $$

where \(R\) is the modulation range (default 0.1). The result is clamped to \([0, 1]\).

- **NC = 0** → EII decreases by 0.05
- **NC = 0.5** → EII unchanged
- **NC = 1** → EII increases by 0.05

The NC score is a weighted average of three KPIs (default: 1/3 each):

$$ NC = w_b \cdot KPI_{biodiv} + w_s \cdot KPI_{soil} + w_w \cdot KPI_{water} $$

Custom weights are supported; they must sum to 1.0.

---

## Client Functions

All functions are in `eii.client` and require `ee.Initialize()`.

### `get_modulated_eii`

Retrieves EII with Natural Capital modulation applied.

**Signature:**

```python
get_modulated_eii(
    geometry,
    kpis=None,
    kpi_layers=None,
    compute_default_kpis=False,
    aggregation_method="min_fuzzy_logic",
    compute_mode="precomputed",
    include_components=True,
    kpi_weights=None,
    modulation_range=0.1,
    biodiversity_max_threshold=0.5,
    scale=300,
    stats=None,
    percentiles=None,
    output_format="dict",
)
```

**KPI input (exactly one required):**

| Option | Description |
|--------|-------------|
| `kpis` | Pre-computed dict `{"biodiversity": float, "soil": float, "water": float}` (0–1). Fastest. |
| `kpi_layers` | User-provided `ee.Image` dict with keys `"biodiversity"`, `"soil"`, `"water"` (0–1). Flexible. |
| `compute_default_kpis=True` | Compute KPIs from built-in GEE datasets. Demo/fallback. |

**Returns:** Dict (or GeoDataFrame if `output_format="geodataframe"`) with `values` containing:

- `eii` — Base EII statistics
- `eii_modulated` — Modulated EII
- `nc_score` — Natural Capital score
- `biodiversity_kpi`, `soil_kpi`, `water_kpi` — Individual KPI statistics

---

### `get_default_kpis`

Computes default KPIs for an AOI and returns polygon-level means (0–1).

**Signature:**

```python
get_default_kpis(
    geometry,
    biodiversity_max_threshold=0.5,
    scale=100,
)
```

**Returns:** `{"biodiversity": float, "soil": float, "water": float}`

---

### `get_kpi_layers`

Returns default KPI raster layers (normalized 0–1) for an AOI.

**Signature:**

```python
get_kpi_layers(
    geometry,
    biodiversity_max_threshold=0.5,
)
```

**Returns:** `{"biodiversity": ee.Image, "soil": ee.Image, "water": ee.Image}`

Use for visualization, custom aggregation, or as `kpi_layers` in `get_modulated_eii`.

---

### `get_nc_score`

Computes only the Natural Capital score (0–1) from KPI inputs.

**Signature:**

```python
get_nc_score(
    geometry,
    kpis=None,
    kpi_layers=None,
    compute_default_kpis=False,
    kpi_weights=None,
    biodiversity_max_threshold=0.5,
    scale=100,
)
```

**Returns:** `float` (0–1). KPI input rules match `get_modulated_eii`: exactly one of `kpis`, `kpi_layers`, or `compute_default_kpis=True`.

---

## KPI Definitions and Data Sources

The KPI layer choices are illustrative and not part of the core EII methodology. They are included to demonstrate how local modulation can work end-to-end, and can be replaced with domain-appropriate KPIs once standards or project-specific definitions are agreed. The defaults below mirror common natural-capital stock proxies in the agrifood sector:

* Biodiversity: percentage of natural habitats
* Soil: soil organic carbon
* Water: soil available water capacity

These should not be interpreted as authoritative KPIs; they are toy examples for demonstrating the workflow.

### Biodiversity KPI

**Concept:** Fraction of natural and semi-natural land cover within the AOI, normalized between 0 and a configurable upper bound (default 50%).

**Data:** ESA WorldCover 10 m (v200). Natural classes: tree cover, shrubland, grassland, bare/sparse vegetation, snow/ice, permanent water, herbaceous wetland, mangroves, moss/lichen. Cropland and built-up are excluded.

**Processing:** Binary natural mask → aggregate to 100 m (1 ha) → fraction → linear normalization:

$$\text{KPI}_{biodiv} = \text{clamp}\left( \frac{f - f_{min}}{f_{max} - f_{min}}, 0, 1 \right)$$

Default: \(f_{min}=0\), \(f_{max}=0.5\) (50% natural = score 1). Configurable via `biodiversity_max_threshold`.

**Resolution:** 100 m.

---

### Soil KPI

**Concept:** Soil organic carbon (SOC) relative to a climate-zone reference, with a minimum ratio floor. Values near or above the reference (well-managed or pristine) map to 1.

**Data:**

- **Actual SOC:** SoilGrids 250 m (`projects/soilgrids-isric/soc_mean`), g/kg (0–30 cm depth by default).
- **Climate zones:** GloH2O Köppen-Geiger v3 (Beck et al. 2023), `projects/landler-open-data/assets/datasets/climatezones/gloh2o-koeppen-v3`, band `b1`.
- **Reference SOC:** Hardcoded table by Köppen class (e.g. Af 45, BWh 8, Cfb 45, ET 150 g/kg). From literature (Jobbágy & Jackson 2000, Post et al. 1982, etc.). Unmapped zones use `DEFAULT_SOC_REFERENCE` (30 g/kg).

**Processing:**

1. Ratio \(r = \text{SOC} / \text{SOC}_{ref}\).
2. Apply minimum threshold \(\tau\) (default 0.1):

$$\text{KPI}_{soil} = \text{clamp}\left( \frac{r - \tau}{1 - \tau}, 0, 1 \right)$$

- \(r \le \tau\) → 0; \(r \ge 1\) → 1 (capped).

**Resolution:** 250 m.

---

### Water KPI

**Concept:** Available water capacity (AWC) relative to a texture-specific range. AWC is the difference between field capacity (33 kPa) and permanent wilting point (1500 kPa). The KPI measures how close the pixel is to its theoretical optimum for its texture.

**Data:** SoilGrids 250 m: `sand_mean`, `clay_mean`, `soc_mean`. AWC is **not** taken from a water-content asset; it is derived via the **Saxton & Rawls (2006)** pedotransfer function from sand, clay, and SOC. SOC is converted to soil organic matter (SOM) as \(\text{SOM} = \text{SOC} \times 1.9\).

**Processing:**

1. **AWC (vol%):** Saxton & Rawls PTF → \(\theta_{33}\), \(\theta_{1500}\) (m³/m³) → \(\text{AWC} = (\theta_{33} - \theta_{1500}) \times 100\).
2. **Texture-specific bounds:**
   - **Max AWC:** Parabolic function of sand/clay around an optimal loam (~35% sand, ~25% clay), clamped to [8, 28] vol%.
   - **Min AWC:** Texture-specific lower bound (heuristic), clamped to [0, 12] vol%.
3. **Normalization:**

$$\text{KPI}_{water} = \text{clamp}\left( \frac{\text{AWC} - \text{AWC}_{min}}{\text{AWC}_{max} - \text{AWC}_{min}}, 0, 1 \right)$$

Values near the texture-specific maximum map to 1.

**Resolution:** 250 m.

---

## Usage Examples

**Modulated EII with default KPIs (GEE):**

```python
import ee
from eii.client import get_modulated_eii

ee.Initialize()
aoi = ee.Geometry.Rectangle([10.5, 47.5, 11.0, 48.0])

result = get_modulated_eii(aoi, compute_default_kpis=True)
print(f"Base EII:      {result['values']['eii']['mean']:.3f}")
print(f"Modulated EII: {result['values']['eii_modulated']['mean']:.3f}")
print(f"NC score:      {result['values']['nc_score']['mean']:.3f}")
```

**Pre-computed KPIs:**

```python
result = get_modulated_eii(aoi, kpis={"biodiversity": 0.45, "soil": 0.62, "water": 0.38})
```

**KPI rasters for mapping or custom aggregation:**

```python
from eii.client import get_kpi_layers

layers = get_kpi_layers(aoi)
# layers["biodiversity"], layers["soil"], layers["water"] are ee.Image (0–1)

result = get_modulated_eii(aoi, kpi_layers=layers, output_format="geodataframe")
```

**Only the NC score:**

```python
from eii.client import get_nc_score

nc = get_nc_score(aoi, compute_default_kpis=True)
```

---

## Data Assets Summary

| KPI | Dataset | GEE asset / source | Resolution |
|-----|---------|--------------------|------------|
| Biodiversity | ESA WorldCover v200 | `ESA/WorldCover/v200` | 10 m → 100 m aggregated |
| Soil (SOC) | SoilGrids 250 m | `projects/soilgrids-isric/soc_mean` | 250 m |
| Soil (climate) | GloH2O Köppen v3 | `projects/landler-open-data/assets/datasets/climatezones/gloh2o-koeppen-v3` | ~1 km |
| Water (sand, clay, SOC) | SoilGrids 250 m | `projects/soilgrids-isric/sand_mean`, `clay_mean`, `soc_mean` | 250 m |

---

## References

- Beck et al. (2023). High-resolution (1 km) Köppen-Geiger maps for 1901–2099. *Scientific Data* 10, 724.
- Jobbágy & Jackson (2000). The vertical distribution of soil organic carbon and its relation to climate and vegetation. *Ecol. Appl.*
- Poggio et al. (2021). SoilGrids 2.0. *SOIL* 7, 217–240.
- Saxton & Rawls (2006). Soil water characteristic estimates by texture and organic matter for hydrologic solutions. *Soil Sci. Soc. Am. J.*
