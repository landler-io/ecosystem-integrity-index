"""
Constants of the v2 EII calculation.

The calculation runs on local rasters (numpy); ``pipelines/production/compute_eii_v2.py`` is
the global tile driver.
"""

from typing import Literal

NPP_STORED_SCALE = 1e-4
"""Multiply stored CLMS NPP by this factor to get gC m-2 yr-1.

The annual pipeline retains CLMS DMP digital numbers with an unapplied 0.01
scale, then multiplies by 100 for integer storage. Stored annual sums are
therefore physical NPP times 10,000. The same factor converts the stored
intra-annual daily standard deviation to gC m-2 day-1.
"""

OBSERVED_NPP_YEARS = (2014, 2025)
"""First and last annual image of CLMS NPP (DMP V2.0.2, ``pipelines/preprocessing/npp``)."""
OBSERVED_NPP_FIRST_YEAR = 2019
"""First annual image EII v2 uses, the first OLCI year. DMP V2.0.2 changes from PROBA-V to OLCI
at 2019 with no overlap and a break of up to 25 % by biome, so earlier years are not
comparable."""
OBSERVED_NPP_WINDOW_YEARS = 3
"""Each reference year averages the 3 annual sums ending in it, from the first OLCI year on."""
OBSERVED_NPP_MIN_VALID_YEARS = 2
"""Pixels with fewer valid (non-negative) annual sums in the window are masked."""
REFERENCE_YEARS = tuple(range(2021, 2026))
"""Published annual EII years. 2020 is not published because its window holds only two OLCI
years (2019-2020). ``observed_npp_window`` accepts 2020 for analyses."""
POTENTIAL_NPP_RESPONSE_YEARS = (2019, 2025)
"""The potential-NPP training response is the mean of these annual sums (inclusive)."""
POTENTIAL_NPP_MIN_VALID_YEARS = 5
"""Of the 7 years behind the potential-NPP training response."""

POTENTIAL_NPP_FLOOR = 50.0
"""gC m-2 yr-1. ``deviation_tolerance_integrity`` divides deviations by max(potential, floor), so
deserts are scored on an absolute scale (relative error is about 100 % below it)."""
TOLERANCE_K = 1.0
"""Deficit tolerance multiplier. Tolerance = K x sqrt(u^2 + cv^2): u is the relative
half-width of the 68 % potential-NPP interval on the side of the deviation, cv the
modelled natural inter-annual CV (optional)."""
TOLERANCE_K_SURPLUS = 1.5
"""Surplus tolerance multiplier. Surplus produces most false alarms on intact land,
so its tolerance is wider than the deficit one."""
SURPLUS_ZERO = 2.0
"""Relative surplus ``(A - P) / max(P, floor)`` at which functional integrity reaches 0, i.e.
three times the potential. A deficit reaches 0 at a relative deviation of 1 (total loss)."""
RECENTRE_MIN_REFERENCE_PIXELS = 200
"""Biomes with this many reference pixels or fewer keep a recentring factor of 1."""
RECENTRE_REFERENCE_HMI_MAX = 0.05
RECENTRE_REFERENCE_BII_MIN = 0.95
"""Default recentring reference: HMI below the maximum and BII at or above the minimum,
taken from the inputs rather than from the run's structural and compositional scores."""
RECENTRE_RATIO_BINS = 400
RECENTRE_RATIO_MAX = 4.0
"""Bins of the observed/potential ratio histogram of ``recentring_sums`` (on [0, max), plus one
overflow bin); the per-bin sums add up over tiles."""
N_BIOMES = 15
"""RESOLVE biome numbers 1-14; 0 is unknown and is never recentred."""
BIOME_NAMES = {
    1: "Tropical moist broadleaf",
    2: "Tropical dry broadleaf",
    3: "Tropical conifer",
    4: "Temperate broadleaf",
    5: "Temperate conifer",
    6: "Boreal",
    7: "Savanna",
    8: "Temperate grassland",
    9: "Flooded grassland",
    10: "Montane grassland",
    11: "Tundra",
    12: "Mediterranean",
    13: "Desert and xeric",
    14: "Mangrove",
}
"""Short names of the RESOLVE biome numbers."""

LANDSCAPE_DECAY_M = 5000.0
"""Exponential kernel scale, 1/beta with beta = 0.2 per km (Beyer et al. 2020)."""
LANDSCAPE_RADIUS_M = 26500.0
"""Kernel truncation, the 99.5 % quantile of the exponential distance (Beyer et al. 2020)."""
LANDSCAPE_OPENING_M = 300.0
"""Grey opening of habitat quality in the landscape term only: patches narrower than this ground
disk (a 3x3 cross at 300 m) add no landscape quality."""
QUALITY_MIDPOINT = 0.4
"""HMI at which the logistic habitat-quality curve crosses one half, calibrated against
WorldCover 2021 natural land cover by ``pipelines/training/calibrate_structural.py``."""
QUALITY_STEEPNESS = 10.0
"""Slope of the logistic quality curve per unit HMI (set, not fitted)."""

METERS_PER_DEGREE = 111319.49079327357
"""Equatorial metres per degree, as Earth Engine converts kernel radii on EPSG:4326."""

WETLAND_FLAG_CLASSES = {
    1: (*range(8, 22), 29, 30, 31, 32),  # vegetated inland and coastal wetland, saltmarsh
    2: tuple(range(22, 28)),  # peatland (arctic/boreal, temperate, tropical)
    3: (28,),  # mangrove
}
"""``wetland_flag`` code -> GLWD v2 classes; open water (1-7) and paddy rice (33) stay 0."""
WETLAND_FLAG_MIN_COVER_PCT = {1: 51, 2: 33, 3: 51}
"""Minimum GLWD cover (percent, inclusive) for each flag code."""

COASTAL_SETTING_KM = 50
"""Coastal predictors are 0 farther inland."""
MAX_TIDAL_RANGE_M = 16.0
"""Bay of Fundy; larger EOT20 values are estuary artefacts (St Lawrence)."""

AggregationMethod = Literal["minimum", "product", "min_fuzzy_logic", "geometric_mean"]
DEFAULT_AGGREGATION_METHOD: AggregationMethod = "min_fuzzy_logic"

# FI change term. Changing any of these values changes the published product.
FI_CHANGE_SERIES = {"TV": ("T", "V"), "TNc": ("T", "Nc"), "Nc": ("Nc",)}
"""Series -> sensors: T MODIS Terra MOD13Q1 061, V VIIRS Suomi NPP VNP13A1 002, Nc VIIRS
NOAA-20 VJ113A1 002 chained to Suomi NPP. ``series_for`` picks the one of a reference year."""
FI_CHANGE_K = {"TV": 1.0, "TNc": 1.0, "Nc": 1.0}
"""k_c per series: FI_change = 1 where 1 - r <= k_c h."""
FI_CHANGE_ZERO_R = 0.0
"""FI_change reaches 0 at this r (linear from the tolerance edge)."""
FI_CHANGE_H_FLOOR = 0.02
"""Floor of the tolerance h predicted by the stored LightGBM models."""
FI_CHANGE_BASELINES = {"T": (2009, 2018), "V": (2012, 2018), "Nc": (2012, 2018)}
"""Fixed baseline years per sensor (inclusive); r divides by the median over them."""
FI_CHANGE_MIN_BASELINE_YEARS = 6
FI_CHANGE_WINDOW_YEARS = 3
FI_CHANGE_MIN_WINDOW_YEARS = 2
"""r = mean of the 3-year window ending in the reference year (>= 2 valid years) divided by the
median of the corrected baseline (>= 6 valid years)."""
FI_CHANGE_CHAIN_YEARS = (2018, 2020)
FI_CHANGE_CHAIN_MIN_YEARS = 2
FI_CHANGE_CHAIN_SWITCH = 2021
"""NOAA-20 chain: f = median of annual V / N over 2018-2020 (>= 2 valid years); Nc = V before
2021 and N x f from 2021."""
EVI_COMPOSITES_PER_YEAR = 23
"""16-day composites on DOY 1 + 16 k; a sensor year is complete with all of them."""
EVI_MIN_COMPOSITES = 15
"""Valid or snow composites a pixel needs for an annual mean."""
EVI_MAX_GAP = 2
"""Gaps of up to this many composites are interpolated linearly."""
EVI_SNOW_PERCENTILE = 5
"""Snow composites take the pixel's 5th percentile of valid EVI over the sensor's record."""
EVI_SNOW_CODE = -20000
EVI_INVALID_CODE = -30000
"""Encoded composites: EVI x 1e4 where valid, else one of these codes."""
FI_CHANGE_CONSTANTS_SHA256 = {
    "drift_factors.csv": "0b58e318849bcbe05996fde5613f3f9992924666267793173fa7eb3a1c9fe9df",
    "h_model_TV.txt": "d2527506b5c0090a57ddb292bf05d13b83ad2e453180e80794e2f9eaece4bf44",
    "h_model_TNc.txt": "6fdbd899754c072cdd80e3a090fcd777f0ebfc8975a6f7816674c7c533fd4231",
    "h_model_Nc.txt": "fed699c403e704228bfe35fa62505e191b53a6d7ee911794436475836567c422",
    "k_c.json": "2648c961d452c927ad2d3b9bddd0c0abba14af0b60022a28f2461ad0297393bf",
}
"""Stored constants in ``eii/compute/fi_change_constants/``, derived by
``pipelines/training/derive_fi_change_constants.py`` and checked at load."""
