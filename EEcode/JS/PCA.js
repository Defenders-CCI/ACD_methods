
// This helper function returns a list of new band names.
function getNewBandNames(prefix, bands) {
  //var seq = ee.List.sequence(1, bands.length());
  //return seq.map(function(b) {
  return bands.map(function(band){
    return ee.String(prefix).cat(band);
  });
 //   return ee.String(prefix).cat(ee.Number(b).int());
  //});
}

function center(img, scale, region){
  var mn = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1000000000
  });
  return img.subtract(mn.toImage(img.bandNames()));
}

function getPrincipalComponents(centered, scale, region) {
  // Collapse the bands of the image into a 1D array per pixel.
  var arrays = centered.toArray();

  // Compute the covariance of the bands within the region.
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });

  // Get the 'array' covariance result and cast to an array.
  // This represents the band-to-band covariance within the region.
  var covarArray = ee.Array(covar.get('array'));

  // Perform an eigen analysis and slice apart the values and vectors.
  var eigens = covarArray.eigen();

  // This is a P-length vector of Eigenvalues.
  var eigenValues = eigens.slice(1, 0, 1);
  //print(eigenValues);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigens.slice(1, 1);
  //print(eigenVectors);
  // Convert the array image to 2D arrays for matrix computations.
  var arrayImage = arrays.toArray(1);

  // Left multiply the image array by the matrix of eigenvectors.
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd', centered.bandNames())]);

  // Turn the PCs into a P-band image, normalized by SD.
  return principalComponents
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc', centered.bandNames())])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
}

exports.getNewBandNames = getNewBandNames;
exports.center = center;
exports.getPCs = getPrincipalComponents;

exports.doc = 'PCA contains functions to compute principal component scores for multi band images.'+
  '\n getNewBandNames(prefix, bands):'+
  '\n   @param {ee.String} prefix String to append on to each element'+
  '\n   @param {ee.List} bands Band names to be concatenated with prefix'+
  '\n   @return {ee.List} List of modified band names.'+
  '\n center(img):'+
  '\n   @param {ee.Image} centered Image'+
  '\n   @param {Integer} scale Scale for reduction computations in meters (e.g. 30)'+
  '\n   @param {ee.Geometry} region Bounding geometry.'+
  '\n   @return {ee.Image} Image with mean-value centered bands.'+
  '\n getPCs(centered, scale, region):'+
  '\n   @param {ee.Image} centered Mean centered n-band image from which to extract principal components'+
  '\n   @param {Integer} scale Scale for reduction computations in meters (e.g. 30)'+
  '\n   @param {ee.Geometry} region Bounding geometry.'+
  '\n   @return {ee.Image} N-band image containing principal component scores.';


// EXAMPLE: Get the PCs at the specified scale and in the specified region
/*
var single = ee.Image(S2.filterBounds(geometry).limit(1, 'CLOUDY_PIXEL_PERCENTAGE', true).first()).clip(geometry);
var med = S2.filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 0.3).median().clip(geometry);
var centered = center(single);
var medcenter = center(med);
Map.addLayer(centered, {bands:['B4', 'B3', 'B2']}, 'centered');
Map.addLayer(medcenter, {bands:['B4', 'B3', 'B2']}, 'single');
var pcImage = getPrincipalComponents(centered.select(['B1', 'B2', 'B3', 'B4', 'B6', 'B8', 'B11', 'B12']), 30, geometry);
var pcMed = getPrincipalComponents(medcenter.select(['B1', 'B2', 'B3', 'B4', 'B6', 'B8', 'B11', 'B12']), 30, geometry);
Map.addLayer(pcImage, {}, 'PCAs');
Map.addLayer(pcMed, {}, 'medPCAs');
*/
