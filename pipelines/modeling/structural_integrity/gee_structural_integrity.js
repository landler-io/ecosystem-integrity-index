/**
 * Structural Integrity Calculator (Quality-Weighted Core Area)
 *
 * This script calculates structural integrity using quality-weighted core area:
 * - Habitat mask from HMI (< threshold = natural/semi-natural)
 * - Erosion by edge depth to identify core habitat
 * - Weight core pixels by habitat quality class
 * - Calculate weighted core proportion within neighborhood
 *
 * Quality Classes:
 *   HMI < 0.1: Pristine (weight 4)
 *   HMI 0.1-0.2: Low-impact (weight 3)
 *   HMI 0.2-0.3: Moderate (weight 2)
 *   HMI 0.3-0.4: Semi-natural (weight 1)
 *   HMI >= 0.4: Modified (weight 0)
 *
 * Score interpretation:
 *   1.0 = All pristine core habitat
 *   0.25 = All semi-natural core habitat
 *   0.0 = No core habitat (fragmented or modified)
 *
 * Usage:
 *   1. Define your AOI below (or draw on the map and import as 'geometry')
 *   2. Set EXPORT_TO_DRIVE = true if you want to export
 *   3. Run the script
 */

// =============================================================================
// CONFIGURATION
// =============================================================================

var INTERACTIVE_MODE = false;  // Set to false for global export (disables map/stats)

// Area of Interest
var AOI;
if (INTERACTIVE_MODE) {
  // Interactive mode: Use smaller region/drawing for testing
  AOI = ee.Geometry.Rectangle([-10, 35, 5, 45]); // Example: Iberian Peninsula
  // AOI = geometry; // Uncomment if you drew a geometry on the map
} else {
  // Batch mode: Global export (matching grid boundaries: -60 to 80 latitude)
  AOI = ee.Geometry.Rectangle([-180, -60, 180, 80], 'EPSG:4326', false);
}

// Core area parameters
var EDGE_DEPTH_M = 300;       // Edge effect penetration depth (meters)
var NEIGHBORHOOD_M = 5000;    // Landscape analysis radius (meters)
var HMI_THRESHOLD = 0.4;      // Values below this = natural/semi-natural habitat
var SCALE_M = 300;            // Output resolution (meters)

// Quality class thresholds and weights
var QUALITY_THRESHOLDS = [0.1, 0.2, 0.3, 0.4];  // Pristine, Low-impact, Moderate, Semi-natural
var QUALITY_WEIGHTS = [4, 3, 2, 1];              // Higher weight = higher quality
var MAX_WEIGHT = 4;                              // For normalization

// Export settings
var EXPORT_TO_DRIVE = true;  // Set to true to export to Drive
var DRIVE_FOLDER = 'structural_integrity_outputs';
var EXPORT_DESCRIPTION = 'structural_integrity_corearea';

// Visualization
var VIS_PARAMS = {
  min: 0,
  max: 1,
  palette: ['#d73027', '#fc8d59', '#fee08b', '#d9ef8b', '#91cf60', '#1a9850']
};

// =============================================================================
// DATA SOURCES
// =============================================================================

var HMI_ASSET = 'projects/sat-io/open-datasets/GHM/HM_2022_300M';

/**
 * Load Global Human Modification Index (2022, all threats combined)
 */
function loadHMI() {
  return ee.ImageCollection(HMI_ASSET)
    .filter(ee.Filter.eq('threat_code', 'AA'))
    .first();
}

// =============================================================================
// CORE AREA CALCULATION
// =============================================================================

/**
 * Create quality class image from HMI.
 * Higher values = higher habitat quality.
 *
 * @param {ee.Image} hmi - Human Modification Index image
 * @returns {ee.Image} Quality class (0-4)
 */
function createQualityClass(hmi) {
  // Start with 0 (modified/non-habitat)
  var quality = ee.Image(0);

  // Assign quality classes based on HMI thresholds
  // Order matters: most restrictive (pristine) last to override
  quality = quality.where(hmi.lt(0.4), 1);  // Semi-natural
  quality = quality.where(hmi.lt(0.3), 2);  // Moderate
  quality = quality.where(hmi.lt(0.2), 3);  // Low-impact
  quality = quality.where(hmi.lt(0.1), 4);  // Pristine

  return quality.rename('quality_class');
}

/**
 * Calculate structural integrity using quality-weighted core area.
 *
 * Core area is habitat that survives erosion by the edge depth. Each core
 * pixel is weighted by its habitat quality class (pristine = 4, semi-natural = 1).
 * The final score reflects both fragmentation AND habitat quality.
 *
 * @param {ee.Geometry} aoi - Area of interest (optional, for clipping)
 * @param {number} edgeDepthM - Edge effect penetration depth in meters
 * @param {number} neighborhoodM - Landscape analysis radius in meters
 * @param {number} scaleM - Output scale in meters for reprojection
 * @returns {ee.Image} Structural integrity score (0-1)
 */
function calculateStructuralIntegrity(aoi, edgeDepthM, neighborhoodM, scaleM) {
  var hmi = loadHMI();

  // Binary habitat classification (HMI < threshold = natural/semi-natural)
  var habitatBinary = hmi.lt(HMI_THRESHOLD).rename('habitat');

  // Find core habitat (survives erosion by edge depth)
  // focal_min with binary mask effectively erodes the habitat patches
  var coreHabitat = habitatBinary.focal_min({
    radius: edgeDepthM,
    kernelType: 'circle',
    units: 'meters'
  }).rename('core');

  // Create quality class image
  var qualityClass = createQualityClass(hmi);

  // Weight core pixels by their quality class
  // Non-core pixels (edges and non-habitat) get 0, core pixels get their quality weight
  // Do NOT mask - we want 0s to contribute to the neighborhood mean so that
  // non-habitat pixels still get influenced by nearby core habitat
  var weightedCore = coreHabitat.multiply(qualityClass);

  // Calculate weighted core proportion within neighborhood
  // Divide by MAX_WEIGHT to normalize to 0-1 range
  var structuralIntegrity = weightedCore.reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.circle({
      radius: neighborhoodM,
      units: 'meters'
    })
  })
  .divide(MAX_WEIGHT)
  .unmask(0)
  .reproject({crs: 'EPSG:4326', scale: scaleM})
  .rename('structural_integrity');

  // Clip to AOI if provided
  if (aoi) {
    structuralIntegrity = structuralIntegrity.clip(aoi);
  }

  return structuralIntegrity;
}

// =============================================================================
// MAIN EXECUTION
// =============================================================================

print('=== Structural Integrity Calculator ===');
print('Edge Depth:', EDGE_DEPTH_M, 'm');
print('Neighborhood:', NEIGHBORHOOD_M, 'm');
print('HMI Threshold:', HMI_THRESHOLD);
print('Scale:', SCALE_M, 'm');

// Calculate structural integrity
var structuralIntegrity = calculateStructuralIntegrity(
  AOI,
  EDGE_DEPTH_M,
  NEIGHBORHOOD_M,
  SCALE_M
);

// Load supporting layers for visualization
var hmi = loadHMI();
var habitatBinary = hmi.lt(HMI_THRESHOLD);
var qualityClass = createQualityClass(hmi);
var coreHabitat = habitatBinary.focal_min({
  radius: EDGE_DEPTH_M,
  kernelType: 'circle',
  units: 'meters'
});
var weightedCore = coreHabitat.multiply(qualityClass);

// =============================================================================
// VISUALIZATION
// =============================================================================

if (INTERACTIVE_MODE) {
  Map.centerObject(AOI);

  // Add layers to map
  Map.addLayer(hmi.clip(AOI), {min: 0, max: 1, palette: ['green', 'yellow', 'red']},
    'Human Modification Index', false);

  Map.addLayer(habitatBinary.clip(AOI).selfMask(), {palette: ['#228b22']},
    'Habitat (HMI < ' + HMI_THRESHOLD + ')', false);

  Map.addLayer(coreHabitat.clip(AOI).selfMask(), {palette: ['#006400']},
    'Core Habitat (after ' + EDGE_DEPTH_M + 'm erosion)', false);

  Map.addLayer(qualityClass.clip(AOI).selfMask(),
    {min: 1, max: 4, palette: ['#fdae61', '#fee08b', '#a6d96a', '#1a9850']},
    'Quality Class (1=Semi-natural, 4=Pristine)', false);

  Map.addLayer(weightedCore.clip(AOI).selfMask(),
    {min: 1, max: 4, palette: ['#fdae61', '#fee08b', '#a6d96a', '#1a9850']},
    'Quality-Weighted Core', false);

  Map.addLayer(structuralIntegrity, VIS_PARAMS,
    'Structural Integrity (Quality-Weighted Core)');

  // Add legend
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '8px 15px'
    }
  });

  var legendTitle = ui.Label({
    value: 'Structural Integrity',
    style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0'}
  });
  legend.add(legendTitle);

  var palette = VIS_PARAMS.palette;
  var labels = ['0 (Fragmented)', '0.2', '0.4', '0.6', '0.8', '1 (Intact)'];

  for (var i = 0; i < palette.length; i++) {
    var colorBox = ui.Label({
      style: {
        backgroundColor: palette[i],
        padding: '8px',
        margin: '0 4px 4px 0'
      }
    });
    var description = ui.Label({
      value: labels[i],
      style: {margin: '0 0 4px 6px'}
    });
    var row = ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.Flow('horizontal')
    });
    legend.add(row);
  }

  Map.add(legend);
}

// =============================================================================
// STATISTICS
// =============================================================================

if (INTERACTIVE_MODE) {
  // Calculate area-weighted statistics for the AOI
  // NOTE: For global runs, we increase the scale to avoid timeout/memory errors in the browser.
  // Uses 10km scale for quick stats, but keeps 300m for export.
  var stats = structuralIntegrity.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.stdDev(), '', true)
      .combine(ee.Reducer.percentile([10, 25, 50, 75, 90]), '', true),
    geometry: AOI,
    // scale: SCALE_M,
    scale: 10000,
    maxPixels: 1e9
  });

  print('=== AOI Statistics ===');
  print('Mean:', stats.get('structural_integrity_mean'));
  print('Std Dev:', stats.get('structural_integrity_stdDev'));
  print('10th percentile:', stats.get('structural_integrity_p10'));
  print('25th percentile:', stats.get('structural_integrity_p25'));
  print('Median:', stats.get('structural_integrity_p50'));
  print('75th percentile:', stats.get('structural_integrity_p75'));
  print('90th percentile:', stats.get('structural_integrity_p90'));
}

// =============================================================================
// EXPORT (Optional)
// =============================================================================

if (EXPORT_TO_DRIVE) {
  // Add metadata
  var exportImage = structuralIntegrity.set({
    'edge_depth_m': EDGE_DEPTH_M,
    'neighborhood_m': NEIGHBORHOOD_M,
    'hmi_threshold': HMI_THRESHOLD,
    'scale_m': SCALE_M,
    'method': 'quality_weighted_core_area',
    'quality_thresholds': QUALITY_THRESHOLDS.join(','),
    'quality_weights': QUALITY_WEIGHTS.join(','),
    'hmi_source': HMI_ASSET,
    'export_date': ee.Date(Date.now()).format('YYYY-MM-dd')
  });

  Export.image.toDrive({
    image: exportImage,
    description: EXPORT_DESCRIPTION,
    folder: DRIVE_FOLDER,
    region: AOI,
    scale: SCALE_M,
    maxPixels: 1e13,
    crs: 'EPSG:4326',
    formatOptions: {
      cloudOptimized: true,
      compression: 'DEFLATE'
    },
    skipEmptyTiles: true
  });

  print('=== Export Task Created ===');
  print('Check the Tasks tab to start the export.');
  print('Drive Folder:', DRIVE_FOLDER);
} else {
  print('=== Export Disabled ===');
  print('Set EXPORT_TO_DRIVE = true to enable export.');
}
