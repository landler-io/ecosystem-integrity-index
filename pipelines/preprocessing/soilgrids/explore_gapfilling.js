// Compare original and gap-filled SoilGrids sand content (15-30cm)

var original = ee.Image("projects/soilgrids-isric/sand_mean")
    .select("sand_15-30cm_mean");

var gapfilled = ee.Image("projects/landler-open-data/assets/datasets/soilgrids/sand/sand_15-30cm_mean_gapfilled");

var vis = {
  min: 0,
  max: 800,
  palette: ["#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"]
};

Map.addLayer(original, vis, "Sand (original)");
Map.addLayer(gapfilled, vis, "Sand (gap-filled)");

var diff = gapfilled.subtract(original);
Map.addLayer(diff.selfMask(), {min: -100, max: 100, palette: ["blue", "white", "red"]}, "Difference", false);

var originalMask = original.mask();
var filledOnly = gapfilled.updateMask(originalMask.not());
Map.addLayer(filledOnly, vis, "Filled pixels only", false);

Map.setCenter(10, 50, 5);
