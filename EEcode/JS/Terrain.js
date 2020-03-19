//Function calculating hillshade for an image using the image's solar angle properties
function shade(img, elev){
  var az = img.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var solel = ee.Number(90).subtract(img.get('MEAN_SOLAR_ZENITH_ANGLE'));
  var shd = ee.Terrain.hillshade(elev, az, solel);
  return shd;
}

function illuminate(img, elev){
  var az = ee.Image.constant(img.get('MEAN_SOLAR_AZIMUTH_ANGLE'));
  var solel = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'));
  var aspect = ee.Terrain.aspect(elev);//.clip(geometry);
  var slope = ee.Terrain.slope(elev);//.clip(geometry);
  var illum = slope.cos().multiply(solel.cos()).add(slope.sin().multiply(solel.sin()).multiply(az.subtract(aspect).cos()));
  return illum.rename('illumination');
}

//Function calculating and correcting hillshade for each image in a collection.
//Return Image Collection
function emp_stat(imgCol, bands, aoi, elev){
  var otherbands = ee.Image(imgCol.first()).bandNames().removeAll(bands);
  var names = ee.Image(imgCol.first()).bandNames();
  var illumCol = imgCol.map(function(image){
    var illumination = illuminate(image, elev.clip(aoi));
    return image.addBands(ee.Image(1)).addBands(illumination);
  });

  var regbands = ee.List(['constant', 'illumination']).cat(bands);
  var coeffs = illumCol.select(regbands).reduce(ee.Reducer.linearRegression({
    numX: 2,
    numY: ee.List(bands).length()})
    ).select('coefficients');
  var betas = coeffs
//Coefficients is 2d array (numX by numY). Return only the slope coefficients
  .arraySlice(0, 1, 2)
//Slope coefficients are stil 2d (1xnumY).  Transform to 1d array of length numY
  .arrayProject([1])
//Transform the array image into multiband image
  .arrayFlatten([bands]);

  var int = coeffs.arraySlice(0, 0, 1).arrayProject([1]).arrayFlatten([bands]);

  var pavg = illumCol.select(bands).reduce(ee.Reducer.mean());

  var corrected = illumCol.map(function(image){
    //Predicted reflectance based on regression equation and observed hillshade
    var pred = image.select('illumination').multiply(betas.select(bands)).add(int.select(bands));

    //Multiply difference by the slope(s) and subtract from current value
    var correction = image.select(bands).subtract(pred).add(pavg);//.copyProperties(image);
    return image.select(otherbands).addBands(correction);
  });

  return corrected;
}

//Function calculating and correcting hillshade for each image in a collection.
//Return Image Collection
function c_correct(imgCol, bands, aoi, elev){
  var otherbands = ee.Image(imgCol.first()).bandNames().removeAll(bands);
  //var names = ee.Image(imgCol.first()).bandNames();
  var illumCol = imgCol.map(function(image){
    var illumination = illuminate(image, elev.clip(aoi));
    return image.addBands(ee.Image(1)).addBands(illumination);
  });

  var regbands = ee.List(['constant', 'illumination']).cat(bands);
  var coeffs = illumCol.select(regbands).reduce(ee.Reducer.linearRegression({
    numX: 2,
    numY: ee.List(bands).length()})
    ).select('coefficients');
  var betas = coeffs
//Coefficients is 2d array (numX by numY). Return only the slope coefficients
  .arraySlice(0, 1, 2)
//Slope coefficients are stil 2d (1xnumY).  Transform to 1d array of length numY
  .arrayProject([1])
//Transform the array image into multiband image
  .arrayFlatten([bands]);

  var int = coeffs.arraySlice(0, 0, 1).arrayProject([1]).arrayFlatten([bands]);

  var c = int.divide(betas);

  var mnshade = illumCol.select('illumination').reduce(ee.Reducer.mean());

  var corrected = illumCol.map(function(image){
    //cosine of zenith
    var zen = ee.Image.constant(image.get('MEAN_SOLAR_ZENITH_ANGLE')).cos();
    var num = zen.add(c).divide(image.select('illumination').add(c));
    //Difference between shading in current image, and mean value (i.e. brightness) at pixel
    //var diffshade = image.select('hillshade').subtract(mnshade);
    //Multiply difference by the slope(s) and subtract from current value
    var correction = image.select(bands).multiply(num);//.copyProperties(image);
    return image.select(otherbands.cat(['illumination'])).addBands(correction);
  });

  return corrected;
}

exports.shade = shade;
exports.c_correct = c_correct;
exports.emp_stat = emp_stat;
exports.illumiate = illuminate;

exports.doc = 'Terrain contains two functions, "shade" and "terrainCorrect2". These calculate the hillhade'+
'\n of a single image, and correct reflectance values for all images in a collection based on hillshade'+
'\n shade(img, elev):'+
'\n   @param {ee.Image} img Image with solar azimuth and zenith properties'+
'\n   @param {ee.Image} elev Digital elevation model (e.g. SRTM)'+
'\n   @return {ee.Image} Output from ee.Terrain.hillshade'+
'\n'+
'\n illuminate(img, elev): calculates illumination as cos incidence angle'+
'\n   @param {ee.Image} img Image with solar azimuth and zenith properties'+
'\n   @param {ee.Image} elev Digital elevation model (e.g. SRTM)'+
'\n   @return {ee.Image} Single band image with per pixel illumination values [-1:1]'+
'\n'+
'\n c_correct(imgCol, bands, aoi, elev): computes C-correction (Telleit 1982)'+
'\n   @param {ee.ImageCollection} imgCol Image collection over which to apply terrain shadow correction.'+
'\n   @param {ee.List} bands List of band names to be corrected'+
'\n   @param {ee.Geometry} aoi Bounding geometry'+
'\n   @param {ee.Image} elev Digital elevation model (e.g. SRTM)'+
'\n   @return {ee.ImageCollection} Image collection contianing selected bands, corrected for terrain shadow.'+
'\n em_stat(imgCol, bands, aoi, elev): computes empirical-statistical correction (Meyer et al. 1993)'+
'\n   @param {ee.ImageCollection} imgCol Image collection over which to apply terrain shadow correction.'+
'\n   @param {ee.List} bands List of band names to be corrected'+
'\n   @param {ee.Geometry} aoi Bounding geometry'+
'\n   @param {ee.Image} elev Digital elevation model (e.g. SRTM)'+
'\n   @return {ee.ImageCollection} Image collection contianing selected bands, corrected for terrain shadow.';

//EXAMPLE
/*
var corrected = c_correct(S2.filterBounds(geometry), ['B2', 'B3', 'B4', 'B8', 'B10', 'B11', 'B12'], geometry, SRTM);
var doi = '2017-01-01';
var today = ee.Date(Date.now());
var after = ee.Image(corrected.filterDate(doi, today).limit(1, 'CLOUDY_PIXEL_PERCENTAGE').first()).clip(geometry);
var before = ee.Image(corrected.filterDate('2014-01-01', doi).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 0.6).limit(1, 'system:time_start', false).first()).clip(geometry);

Map.addLayer(before, {}, 'before');
Map.addLayer(after, {}, 'after');
*/
