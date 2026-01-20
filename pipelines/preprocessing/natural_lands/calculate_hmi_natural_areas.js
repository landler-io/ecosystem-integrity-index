/**
 * HMI-Based Natural Areas Identification
 *
 * Identifies natural areas using HMI v3 time-series threshold, similar to wilderness definition.
 * This captures smaller natural features (river beds, floodplains, small patches) that aren't
 * formally protected but have very low human modification across the entire time period (1990-2020).
 *
 * Processing steps:
 * 1. Load all years of HMI v3 time-series (1990-2020, threat_code "AA")
 * 2. Take pixelwise maximum across all years (worst-case modification over time)
 * 3. Fill missing values using 5x5 pixel kernel with max-aggregation
 * 4. Buffer high-modification areas (HMI > threshold) by specified pixels
 * 5. Identify low-modification areas outside buffered zones
 * 6. Filter by minimum area and convert to polygons
 */

// ============================================================================
// CONFIGURATION PARAMETERS
// ============================================================================

// HMI thresholds
var HMI_THRESHOLD = 0.1;              // Maximum HMI value to consider "natural" (lower = stricter)
var HIGH_MOD_THRESHOLD = 0.4;         // HMI threshold for high-modification areas to buffer
var BUFFER_PIXELS = 10;                // Number of pixels to buffer high-modification areas (10 pixels = 3km at 300m resolution)

// Area filtering
var MIN_AREA_KM2 = 1;                // Minimum area in km² to include (filters out tiny fragments)
var SIMPLIFY_TOLERANCE = 1000;      // Geometry simplification tolerance in meters

// AOI for testing (set to null for global computation)
// Option 1: Use a rectangle (uncomment and modify coordinates)

// Option 2: Use a polygon from an asset (uncomment and modify path)

// Option 3: Global computation (no AOI constraint)
var AOI = (typeof geometry !== 'undefined') ? geometry : null;

// Controls for immediate single-AOI computation
var RUN_SINGLE_AOI = true;  // Set to true/false explicitly as needed

// Controls for batch exporting via tiling
var ENABLE_GRID_EXPORT = false;   // Set to true to start Drive exports
var GRID_TILE_SIZE_DEGREES = 50;  // Tile size in degrees (e.g., 5, 10)
var DRIVE_FILE_FORMAT = "GeoJSON";
var START_TILE_INDEX = 0;         // Start index within grid list
var END_TILE_INDEX = null;        // Set to a number to limit tiles per run, or null for all

// Raster export settings
var EXPORT_RASTER = false;         // Export binary raster mask in addition to vectors
var RASTER_SCALE = 300;           // Export resolution in meters (match HMI resolution)

// ============================================================================
// DYNAMIC EXPORT NAMING (encodes parameters for variant comparison)
// ============================================================================

/**
 * Format a number for use in folder/file names.
 * Replaces decimal points with 'p' (e.g., 0.01 -> "0p01")
 */
function formatParamForName(value) {
  return String(value).replace(".", "p");
}

/**
 * Generate parameter suffix for folder names.
 * Format: hmi{threshold}_high{highMod}_buf{pixels}_area{minArea}
 * Example: hmi0p01_high0p10_buf30_area5
 */
function getParamSuffix() {
  return "hmi" + formatParamForName(HMI_THRESHOLD) +
         "_high" + formatParamForName(HIGH_MOD_THRESHOLD) +
         "_buf" + BUFFER_PIXELS +
         "_area" + MIN_AREA_KM2;
}

// Dynamic folder and prefix names encoding parameters
var PARAM_SUFFIX = getParamSuffix();
var DRIVE_EXPORT_FOLDER = "HMI_Natural_Areas_" + PARAM_SUFFIX;
var RASTER_EXPORT_FOLDER = "HMI_Natural_Areas_Raster_" + PARAM_SUFFIX;
var EXPORT_DESCRIPTION_PREFIX = "hmi_natural_" + PARAM_SUFFIX + "_";
var RASTER_EXPORT_PREFIX = "hmi_mask_" + PARAM_SUFFIX + "_";

// Log the computed names for verification
print("=== Export Configuration ===");
print("Parameters: HMI=" + HMI_THRESHOLD + ", HIGH_MOD=" + HIGH_MOD_THRESHOLD +
      ", BUFFER=" + BUFFER_PIXELS + "px, MIN_AREA=" + MIN_AREA_KM2 + "km²");
print("Vector folder: " + DRIVE_EXPORT_FOLDER);
print("Raster folder: " + RASTER_EXPORT_FOLDER);

// Optional: preview grid tiles directly in the map (small ranges only)
var ENABLE_TILE_PREVIEW = false;
var PREVIEW_TILE_START = 1;
var PREVIEW_TILE_END = 3;  // number of tiles to preview (exclusive)

// ============================================================================
// MAIN FUNCTION
// ============================================================================

function getHmiBasedNaturalAreas(hmiThreshold, minAreaKm2, highModThreshold, bufferPixels, simplifyTolerance, aoi) {
  // Log start
  if (aoi) {
    print("Processing HMI v3 time-series data with threshold " + hmiThreshold + " for AOI");
  } else {
    print("Loading HMI v3 time-series data with threshold " + hmiThreshold + " (global)");
  }

  var processingGeometry = null;
  var outputGeometry = null;
  if (aoi) {
    var bufferDistanceMeters = bufferPixels * 300;  // approximate 300m per pixel
    processingGeometry = aoi.buffer(bufferDistanceMeters);
    outputGeometry = aoi;
  }

  // Load HMI v3 time-series collection (1990-2020, change-consistent series)
  var hmCollection = ee.ImageCollection("projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M");
  // Filter to threat_code "AA" (all threats combined) and all years
  var hmiSeries = hmCollection.filter(ee.Filter.eq("threat_code", "AA"));

  // If AOI is provided, filter collection to buffered AOI bounds for efficiency
  if (processingGeometry) {
    hmiSeries = hmiSeries.filterBounds(processingGeometry);
  }

  // Standardize band names across all images before aggregation
  // Map over collection to select first band and rename to "hmi" for consistency
  var standardizeBands = function(img) {
    return img.select([0]).rename("hmi");
  };

  hmiSeries = hmiSeries.map(standardizeBands);

  // Step 1: Temporal aggregation - take pixelwise maximum across all years
  // This captures worst-case modification over the entire time period
  var hmiMax = hmiSeries.max();

  // Clip to buffered AOI if provided
  if (processingGeometry) {
    hmiMax = hmiMax.clip(processingGeometry);
  }

  // Step 2: Spatial interpolation - fill missing values using 5x5 pixel kernel with max-aggregation
  // Create 5x5 kernel for focal max (2 pixels radius = 5x5 kernel)
  var kernel5x5 = ee.Kernel.square({radius: 2, units: "pixels"});
  var hmiFilled = hmiMax.focalMax({kernel: kernel5x5});

  // Step 3: Spatial distance filter - buffer high-modification areas
  // Identify high-modification areas
  var highModMask = hmiFilled.gt(highModThreshold);

  // Buffer high-modification areas by specified number of pixels
  // Use focal max with circular kernel to create buffer
  var bufferKernel = ee.Kernel.circle({radius: bufferPixels, units: "pixels"});
  var highModBuffered = highModMask.focalMax({kernel: bufferKernel});

  // Step 4: Create natural areas mask
  // Areas with HMI below threshold AND outside buffered high-modification zones
  var lowHmiMask = hmiFilled.lt(hmiThreshold);
  var outsideBufferMask = highModBuffered.eq(0);
  var naturalMask = lowHmiMask.and(outsideBufferMask);

  // Step 5: Convert to polygons
  var reduceGeometry = null;
  if (processingGeometry) {
    // Use buffered AOI bounds for reduceToVectors to constrain computation
    reduceGeometry = processingGeometry.bounds();
  } else {
    reduceGeometry = null;
  }

  var naturalAreas = naturalMask.reduceToVectors({
    geometry: reduceGeometry,  // Constrain to AOI if provided, otherwise global
    scale: 300,                 // Match HMI resolution
    geometryType: "polygon",
    eightConnected: false,
    maxPixels: 1e12
  });

  // Keep only polygons where the mask label equals 1 (natural areas)
  naturalAreas = naturalAreas.filter(ee.Filter.eq("label", 1));

  // If AOI was provided, clip results back to original AOI to remove buffered outskirts
  if (outputGeometry) {
    naturalAreas = naturalAreas
      .map(function(feature) {
        return feature.intersection(outputGeometry, ee.ErrorMargin(1));
      })
      .filter(ee.Filter.bounds(outputGeometry));
  }

// Step 6: Compute area property before filtering
naturalAreas = naturalAreas.map(function(feature) {
  var areaM2 = feature.geometry().area({maxError: 1});
  return feature.set("area_m2", areaM2);
});

// Step 7: Filter by minimum area
  var minAreaM2 = minAreaKm2 * 1e6;
naturalAreas = naturalAreas.filter(ee.Filter.gt("area_m2", minAreaM2));

// Step 8: Simplify geometries for efficiency
  naturalAreas = naturalAreas.map(function(feature) {
    return feature.simplify({maxError: simplifyTolerance});
  });

// Step 9: Add metadata
  var addMetadata = function(feature) {
    var areaKm2 = feature.geometry().area({maxError: 1}).divide(1e6);
    return feature.set({
      "source": "HMI_threshold",
      "hmi_threshold": hmiThreshold,
      "high_mod_threshold": highModThreshold,
      "buffer_pixels": bufferPixels,
      "area_km2": areaKm2
    });
  };

  naturalAreas = naturalAreas.map(addMetadata);

  return naturalAreas;
}

/**
 * Returns the natural areas as a binary raster mask (for sampling workflows).
 * This is more efficient than vectorizing when combining with other raster layers.
 */
function getNaturalAreasMask(hmiThreshold, highModThreshold, bufferPixels, aoi) {
  var processingGeometry = null;
  var outputGeometry = null;
  if (aoi) {
    var bufferDistanceMeters = bufferPixels * 300;
    processingGeometry = aoi.buffer(bufferDistanceMeters);
    outputGeometry = aoi;
  }

  // Load HMI v3 time-series collection
  var hmCollection = ee.ImageCollection("projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M");
  var hmiSeries = hmCollection.filter(ee.Filter.eq("threat_code", "AA"));

  if (processingGeometry) {
    hmiSeries = hmiSeries.filterBounds(processingGeometry);
  }

  // Standardize bands
  var standardizeBands = function(img) {
    return img.select([0]).rename("hmi");
  };
  hmiSeries = hmiSeries.map(standardizeBands);

  // Temporal max
  var hmiMax = hmiSeries.max();
  if (processingGeometry) {
    hmiMax = hmiMax.clip(processingGeometry);
  }

  // Gap-fill with focal max
  var kernel5x5 = ee.Kernel.square({radius: 2, units: "pixels"});
  var hmiFilled = hmiMax.focalMax({kernel: kernel5x5});

  // Buffer high-modification areas
  var highModMask = hmiFilled.gt(highModThreshold);
  var bufferKernel = ee.Kernel.circle({radius: bufferPixels, units: "pixels"});
  var highModBuffered = highModMask.focalMax({kernel: bufferKernel});

  // Create natural areas mask
  var lowHmiMask = hmiFilled.lt(hmiThreshold);
  var outsideBufferMask = highModBuffered.eq(0);
  var naturalMask = lowHmiMask.and(outsideBufferMask).rename("natural_mask");

  // Clip to output geometry if provided
  if (outputGeometry) {
    naturalMask = naturalMask.clip(outputGeometry);
  }

  // Set mask metadata as properties
  naturalMask = naturalMask
    .set("hmi_threshold", hmiThreshold)
    .set("high_mod_threshold", highModThreshold)
    .set("buffer_pixels", bufferPixels);

  return naturalMask;
}

function sanitizeForDescription(text) {
  return text.replace(/[^A-Za-z0-9_]+/g, "_");
}

function createGridMetadata(tileSizeDegrees) {
  var tiles = [];
  for (var lat = -90; lat < 90; lat += tileSizeDegrees) {
    for (var lon = -180; lon < 180; lon += tileSizeDegrees) {
      var minLat = lat;
      var minLon = lon;
      var maxLat = lat + tileSizeDegrees;
      var maxLon = lon + tileSizeDegrees;
      var geom = ee.Geometry.Rectangle([minLon, minLat, maxLon, maxLat], null, false);
      var tileId = "lat_" + minLat.toFixed(2) + "_lon_" + minLon.toFixed(2);
      tiles.push({
        id: tileId,
        safeId: sanitizeForDescription(tileId),
        geometry: geom,
        feature: ee.Feature(geom, { tile_id: tileId }),
      });
    }
  }
  return tiles;
}

var GRID_METADATA = createGridMetadata(GRID_TILE_SIZE_DEGREES);
var GRID_FEATURE_COLLECTION = ee.FeatureCollection(
  GRID_METADATA.map(function(tile) {
    return tile.feature;
  })
);

var gridImage = ee.Image().paint(GRID_FEATURE_COLLECTION, 1, 1);
Map.addLayer(gridImage, {palette: ['black'], opacity: 0.8}, 'Grid Cells');


// ============================================================================
// EXECUTION
// ============================================================================

var hmiNaturalAreas = ee.FeatureCollection([]);
if (RUN_SINGLE_AOI && AOI) {
  hmiNaturalAreas = getHmiBasedNaturalAreas(
    HMI_THRESHOLD,
    MIN_AREA_KM2,
    HIGH_MOD_THRESHOLD,
    BUFFER_PIXELS,
    SIMPLIFY_TOLERANCE,
    AOI
  );

  print("Natural areas FeatureCollection:", hmiNaturalAreas);

  var totalArea = hmiNaturalAreas.aggregate_sum("area_km2");
  print("Total area (km²):", totalArea);
} else {
  print("RUN_SINGLE_AOI disabled or AOI null. Using grid-based workflow.");
}

// ============================================================================
// VISUALIZATION
// ============================================================================

// Load HMI for visualization
var hmCollection = ee.ImageCollection("projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M");
var hmi = hmCollection.filter(ee.Filter.eq("threat_code", "AA")).first().select([0]);

// Create visualization parameters
var hmiVis = {
  min: 0,
  max: 1,
  palette: ["green", "yellow", "orange", "red"]
};

// Create visualization masks using temporal max across all years (same as processing)
var hmiSeriesVis = hmCollection.filter(ee.Filter.eq("threat_code", "AA"))
  .map(function(img) { return img.select([0]).rename("hmi"); });
var hmiMaxVis = hmiSeriesVis.max();

// Low-HMI mask: all pixels below threshold (before buffering)
var lowHmiMaskVis = hmiMaxVis.lt(HMI_THRESHOLD).selfMask();

// High-HMI mask: all pixels above high-mod threshold (triggers buffer)
var highHmiMaskVis = hmiMaxVis.gt(HIGH_MOD_THRESHOLD).selfMask();

// Ensure AOI is a FeatureCollection (even if it is a Geometry)
var aoiParams = ee.FeatureCollection(AOI);

var polystyle = {
  color: 'FF0000',           // Outline color (red)
  width: 2,                  // Outline width (pixels)
  fillColor: '00000000'      // Transparent fill (hex ARGB)
};

// Apply the style function to the FeatureCollection
var styledParams = aoiParams.style(polystyle);


// Add to map
if (RUN_SINGLE_AOI && AOI) {
  Map.centerObject(AOI, 10);
  Map.addLayer(styledParams, {}, 'AOI Outline');
  Map.addLayer(hmi, hmiVis, "HMI (for reference)", false);
  Map.addLayer(lowHmiMaskVis, {palette: ['lightgreen']}, "Low HMI Pixels (< " + HMI_THRESHOLD + ")", false);
  Map.addLayer(highHmiMaskVis, {palette: ['red']}, "High HMI Pixels (> " + HIGH_MOD_THRESHOLD + ")", false);
  Map.addLayer(hmiNaturalAreas, {color: "blue"}, "HMI Natural Areas (filtered)", true);
} else {
  Map.setCenter(0, 0, 2);
  Map.addLayer(hmi, hmiVis, "HMI (for reference)", false);
  print("Single-AOI visualization disabled. Enable RUN_SINGLE_AOI and set AOI to inspect results on the map.");
}

// ============================================================================
// GRID EXPORT WORKFLOW (BATCH EXPORTS TO DRIVE)
// ============================================================================

if (ENABLE_GRID_EXPORT) {
  var totalTiles = GRID_METADATA.length;
  var startIndex = START_TILE_INDEX || 0;
  var endIndex =
    END_TILE_INDEX !== null && END_TILE_INDEX !== undefined
      ? Math.min(END_TILE_INDEX, totalTiles)
      : totalTiles;

  for (var idx = startIndex; idx < endIndex; idx++) {
    var tileMeta = GRID_METADATA[idx];
    if (!tileMeta) {
      continue;
    }
    var tileGeometry = tileMeta.geometry;
    if (AOI) {
      tileGeometry = tileGeometry.intersection(AOI, ee.ErrorMargin(1));
    }

    var tileNaturalAreas = getHmiBasedNaturalAreas(
      HMI_THRESHOLD,
      MIN_AREA_KM2,
      HIGH_MOD_THRESHOLD,
      BUFFER_PIXELS,
      SIMPLIFY_TOLERANCE,
      tileGeometry
    );

    var description = EXPORT_DESCRIPTION_PREFIX + tileMeta.safeId;

    // Export vector polygons
    Export.table.toDrive({
      collection: tileNaturalAreas,
      description: description,
      folder: DRIVE_EXPORT_FOLDER,
      fileNamePrefix: description,
      fileFormat: DRIVE_FILE_FORMAT,
    });

    // Export raster mask
    if (EXPORT_RASTER) {
      var tileMask = getNaturalAreasMask(
        HMI_THRESHOLD,
        HIGH_MOD_THRESHOLD,
        BUFFER_PIXELS,
        tileGeometry
      );

      var rasterDescription = RASTER_EXPORT_PREFIX + tileMeta.safeId;

      Export.image.toDrive({
        image: tileMask.toInt8(),
        description: rasterDescription,
        folder: RASTER_EXPORT_FOLDER,
        fileNamePrefix: rasterDescription,
        region: tileGeometry,
        scale: RASTER_SCALE,
        maxPixels: 1e12,
        fileFormat: "GeoTIFF"
      });
    }

  }

  if (ENABLE_TILE_PREVIEW) {
    var previewStart = PREVIEW_TILE_START || startIndex;
    var previewEnd =
      PREVIEW_TILE_END || Math.min(previewStart + 1, GRID_METADATA.length);
    var previewCollection = ee.FeatureCollection([]);

    for (var p = previewStart; p < previewEnd; p++) {
      var previewMeta = GRID_METADATA[p];
      if (!previewMeta) {
        continue;
      }
      var previewGeom = previewMeta.geometry;
      if (AOI) {
        previewGeom = previewGeom.intersection(AOI, ee.ErrorMargin(1));
      }
      var previewAreas = getHmiBasedNaturalAreas(
        HMI_THRESHOLD,
        MIN_AREA_KM2,
        HIGH_MOD_THRESHOLD,
        BUFFER_PIXELS,
        SIMPLIFY_TOLERANCE,
        previewGeom
      );
      previewCollection = previewCollection.merge(previewAreas);
    }

    Map.addLayer(previewCollection, { color: "cyan" }, "Preview Natural Areas");
  }
} else {
  print("Grid export disabled. Set ENABLE_GRID_EXPORT = true to start Drive exports.");
}
