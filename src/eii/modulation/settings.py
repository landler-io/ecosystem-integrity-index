"""
Settings for Natural Capital modulation of EII.

Contains configuration for computing the three KPI dimensions
(biodiversity, soil, water) and combining them into a Natural Capital score.
"""

# ---------------------------------------------------------------------------
# Modulation Parameters
# ---------------------------------------------------------------------------

# Total modulation range: NC score shifts EII by ±(MODULATION_RANGE/2)
# Default: 0.1 means NC=0 -> -0.05, NC=1 -> +0.05
MODULATION_RANGE = 0.1

# Default KPI weights (must sum to 1.0)
DEFAULT_KPI_WEIGHTS = {
    "biodiversity": 1 / 3,
    "soil": 1 / 3,
    "water": 1 / 3,
}

# ---------------------------------------------------------------------------
# Biodiversity KPI Settings (ESA WorldCover)
# ---------------------------------------------------------------------------

WORLDCOVER_ASSET = "ESA/WorldCover/v200"

# WorldCover class codes considered "natural" or "semi-natural"
# 10=Tree cover, 20=Shrubland, 30=Grassland, 60=Bare/sparse,
# 70=Snow/ice, 80=Water, 90=Wetland, 95=Mangroves, 100=Moss/lichen
NATURAL_LANDCOVER_CLASSES = [10, 20, 30, 60, 70, 80, 90, 95, 100]

# Non-natural classes (for reference): 40=Cropland, 50=Built-up

# Resolution for aggregating 10m pixels to calculate fraction
BIODIVERSITY_AGGREGATION_SCALE = 100  # meters (1 hectare)

# Normalization thresholds for biodiversity KPI
BIODIVERSITY_MIN_THRESHOLD = 0.0  # 0% natural -> score 0
BIODIVERSITY_MAX_THRESHOLD = 0.50  # 50% natural -> score 1 (configurable)

# ---------------------------------------------------------------------------
# Soil KPI Settings (SoilGrids + Köppen Climate Zones)
# ---------------------------------------------------------------------------

# SoilGrids 2.0 assets are published per property in the community catalog.
# See https://gee-community-catalog.org/projects/isric/
SOILGRIDS_SOC_ASSET = "projects/soilgrids-isric/soc_mean"
SOILGRIDS_WV0033_ASSET = "projects/soilgrids-isric/wv0033_mean"
SOILGRIDS_WV1500_ASSET = "projects/soilgrids-isric/wv1500_mean"
SOILGRIDS_SAND_ASSET = "projects/soilgrids-isric/sand_mean"
SOILGRIDS_CLAY_ASSET = "projects/soilgrids-isric/clay_mean"
# GloH2O Koppen-Geiger v3 (Beck et al. 2023) uploaded to EE.
# Legend codes are stored in pipelines/preprocessing/climatezones/legend.txt.
KOPPEN_ASSET = "projects/landler-open-data/assets/datasets/climatezones/gloh2o-koeppen-v3"

# Default soil depth for SOC extraction
SOIL_DEPTH = "0-30cm"

# Reference SOC values by Köppen climate zone (g/kg, 0-30cm depth)
# Values represent typical SOC in undisturbed/natural soils per climate zone
# Derived from literature synthesis (Jobbágy & Jackson 2000, post et al. 1982, etc.)
SOC_REFERENCE_BY_CLIMATE = {
    # Tropical climates (A)
    "Af": 45,  # Tropical rainforest
    "Am": 40,  # Tropical monsoon
    "Aw": 25,  # Tropical savanna
    # Arid climates (B)
    "BWh": 8,  # Hot desert
    "BWk": 10,  # Cold desert
    "BSh": 15,  # Hot steppe
    "BSk": 20,  # Cold steppe
    # Temperate climates (C)
    "Csa": 25,  # Mediterranean hot summer
    "Csb": 35,  # Mediterranean warm summer
    "Csc": 40,  # Mediterranean cold summer
    "Cwa": 30,  # Humid subtropical dry winter
    "Cwb": 35,  # Subtropical highland dry winter
    "Cwc": 40,  # Subpolar oceanic dry winter
    "Cfa": 35,  # Humid subtropical
    "Cfb": 45,  # Oceanic
    "Cfc": 60,  # Subpolar oceanic
    # Continental climates (D)
    "Dsa": 35,  # Mediterranean-influenced hot summer continental
    "Dsb": 40,  # Mediterranean-influenced warm summer continental
    "Dsc": 50,  # Mediterranean-influenced subarctic
    "Dsd": 60,  # Mediterranean-influenced extremely cold subarctic
    "Dwa": 35,  # Monsoon-influenced hot summer continental
    "Dwb": 45,  # Monsoon-influenced warm summer continental
    "Dwc": 60,  # Monsoon-influenced subarctic
    "Dwd": 80,  # Monsoon-influenced extremely cold subarctic
    "Dfa": 40,  # Hot summer continental
    "Dfb": 50,  # Warm summer continental
    "Dfc": 80,  # Subarctic
    "Dfd": 100,  # Extremely cold subarctic
    # Polar climates (E)
    "ET": 150,  # Tundra
    "EF": 50,  # Ice cap (limited soil)
}

# Default SOC reference for unmapped climate zones
DEFAULT_SOC_REFERENCE = 30  # g/kg

# Minimum SOC ratio threshold for normalization (ratio floor).
# Ratios <= threshold map to 0, ratio >= 1 maps to 1.
SOIL_MIN_THRESHOLD = 0.1

# Köppen-Geiger integer codes mapping to climate zone names
# Based on Beck et al. (2023) legend (pipelines/preprocessing/climatezones/legend.txt)
KOPPEN_CODE_TO_NAME = {
    1: "Af",
    2: "Am",
    3: "Aw",
    4: "BWh",
    5: "BWk",
    6: "BSh",
    7: "BSk",
    8: "Csa",
    9: "Csb",
    10: "Csc",
    11: "Cwa",
    12: "Cwb",
    13: "Cwc",
    14: "Cfa",
    15: "Cfb",
    16: "Cfc",
    17: "Dsa",
    18: "Dsb",
    19: "Dsc",
    20: "Dsd",
    21: "Dwa",
    22: "Dwb",
    23: "Dwc",
    24: "Dwd",
    25: "Dfa",
    26: "Dfb",
    27: "Dfc",
    28: "Dfd",
    29: "ET",
    30: "EF",
}

# ---------------------------------------------------------------------------
# Water KPI Settings (Saxton & Rawls AWC from SoilGrids)
# ---------------------------------------------------------------------------

# AWC is derived via Saxton & Rawls (2006) pedotransfer functions from
# sand, clay, and SOC (converted to SOM via som = soc * 1.9).

# Texture-based AWC bounds (vol%)
# AWC varies with soil texture; loamy soils have highest AWC
# Sandy and clayey extremes have lower AWC
AWC_BASE_MAX = 25.0  # Maximum AWC for ideal loam (vol%)
AWC_MIN_BOUND = 8.0  # Minimum possible max AWC (sandy soils)
AWC_MAX_BOUND = 28.0  # Maximum possible max AWC (silt loam)

# Texture-based minimum AWC parameters (heuristic lower bound)
AWC_BASE_MIN = 5.0  # Minimum AWC for ideal loam (vol%)
AWC_MIN_MIN_BOUND = 0.0  # Minimum possible min AWC
AWC_MIN_MAX_BOUND = 12.0  # Maximum possible min AWC

# Texture penalty coefficients for AWC normalization
# Penalty increases as texture deviates from optimal (loam ~35% sand, ~25% clay)
AWC_OPTIMAL_SAND = 35.0  # Optimal sand content (%)
AWC_OPTIMAL_CLAY = 25.0  # Optimal clay content (%)
AWC_SAND_PENALTY = 0.005  # Penalty coefficient for sand deviation
AWC_CLAY_PENALTY = 0.008  # Penalty coefficient for clay deviation
