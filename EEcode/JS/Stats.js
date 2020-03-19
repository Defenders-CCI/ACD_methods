function normalize(img){
  var bands = img.bandNames();
  var mean = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    scale: 300
  }).toImage()
  var sd = img.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    scale: 300
  }).toImage(bands)

  return img.subtract(mean).divide(sd);
}

exports.normalize = normalize;

//Function calulating CDF probability of chi-square statistic
//param chi: single band image of scores from standard chi-squared distribution
//param df: degrees of freedom (usually # of bands)
//Return: single band image of probabilities
exports.chi_p = function(chi, df){
  var cdf =  ee.Image(chi.divide(2)).gammainc(ee.Number(df).divide(2));
  return cdf.rename(['p']);
};

exports.gamma_p = function(stat, df){
  var shape = ee.Image(1);
  var scale = ee.Image(df);
  var denom = shape.gamma();
  var num = shape.gammainc(stat.divide(scale));
  return num.divide(denom).rename(['p']);
};

//Function calulating normal probability of z-score
// using normal CDF approximation: p = 1/1+exp(-1.65451*z)
//param z: single band image of z-scores from standard normal distribution
//Return: single band image of probabilities
exports.norm_p = function(z){
  return ee.Image.constant(1).subtract(z.multiply(-1.65451).exp().add(1).pow(-1));
};

exports.norm_pdf = function(x, u, s){
  x = ee.Image(x);
  u = ee.Image(u);
  s = ee.Image(s);
  var pow = x.subtract(u).pow(2).divide(s.pow(2).multiply(2)).multiply(-1);
  var first = s.pow(2).multiply(Math.PI).multiply(2).sqrt().pow(-1);
  return pow.exp().multiply(first);
};

exports.ND = function(img, NIR, R, G, SWIR1, SWIR2){
  var NBR = img.normalizedDifference([NIR, SWIR2]).rename(["nbr"]);
  var NDSI = img.normalizedDifference([G, SWIR1]).divide(img.select([NIR])).rename(["ndsi"]);
  //var NNDSI = norm(NDSI, ["ndsi"], region).rename(["nndsi"]);
  var NDWI = img.normalizedDifference([G, NIR]).rename(['ndwi']);
  var NDVI = img.normalizedDifference([NIR, R]).rename(["ndvi"]);

  return img.addBands(NDVI)
  .addBands(NDSI)
  .addBands(NBR)
  .addBands(NDWI);
};

//Function converting multiband image into single band image of LDA scores
//param img: multiband image (ee.Image)
//param int: intercept parameter from LDA analysis (numeric)
//param xbands: string list of n band names (list)
//param coefficients: numeric list of length n containing LDA coefficients (list)
exports.ldaScore = function(img, int, xbands, coefficients){
  var bands = img.select(xbands);

  var coeffs = ee.Dictionary.fromLists(xbands, coefficients).toImage(xbands);

  var score = bands.multiply(coeffs).addBands(ee.Image(int)).reduce(ee.Reducer.sum());

  return score;
  };

//Function calculating PDF of bivariate normal distribution
function bvNorm(imgt, imgs){
  imgt = imgt.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']);
  imgs = imgs.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']);
  var corr = imgt.bandNames().map(function(band){
    var img = imgt.select([band]).addBands(imgs.select([band]));
    var corr = img.reduceRegion({
      reducer: ee.Reducer.pearsonsCorrelation(),
      geometry: aoi,
      scale: scl,
      maxPixels: max
    }).toImage(['correlation']).rename([band]);
    return corr;
  });

  corr = ee.Image.cat(corr.get(0), corr.get(1), corr.get(2), corr.get(3), corr.get(4), corr.get(5));
  print(corr);
  var mns = standardize(imgs, aoi, scl, 1);
  var sds = imgs.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: scl,
    maxPixels: max
  }).toImage(imgs.bandNames());

  var mnt = ee.Image(0);
  var sdt = ee.Image(1);
  //var corr = var1.addBands(var2).reduceRegion({
  //  reducer: ee.Reducer.pearsonsCorrelation(),
  //  geometry: aoi,
  //  scale: scl,
  //  maxPixels: max
  //}).toImage(['correlation']);

  var corr2 = ee.Image(1).subtract(corr.pow(2));

  var dt = imgt.subtract(mnt);

  var ds = imgs.subtract(mns);

  var zt = dt.pow(2).divide(sdt.pow(2));

  var z2 = dt.multiply(ds).multiply(corr).multiply(2).divide(sdt.multiply(sds));

  var zs = ds.pow(2).divide(sds.pow(2));

  var z = zt.subtract(z2).add(zs);

  var exp = z.divide(corr2.multiply(-2)).exp();

  var first = corr2.sqrt().multiply(sdt).multiply(sds).multiply(Math.PI).multiply(2).pow(-1);

  var p = first.multiply(exp);
  return p;
}


exports.doc = 'Stats contains helpful functions for calculating statistical distributions'+
'\n ND(img, NIR, R, G, SWIR1, SWIR2): calculates normalized difference metrics (NDVI, NDSI, NDWI, NBR)'+
'\n   @param {ee.Image} img Image containing near ifrared, red, green, and two shorwavve infrared bands.'+
'\n   @param {ee.String} NIR Band name containing near-infrared reflectance.'+
'\n   @param {ee.String} R Band name containing red reflectance.'+
'\n   ...'+
'\n   @return {ee.Image} Image containing original bands plus 4 ND metrics'+
'\n ldaScore(img, int, xbands, coefficients): calculates lda discriminant score given list of coefficients'+
'\n   @param {ee.Image} img Image with bands corresponding to lda variables'+
'\n   @param {Numeric} int Intercept for LDA function'+
'\n   @param {ee.List} xbands Band names containing variables used in LDA.'+
'\n   @param {ee.List} coefficients Numeric list containing LDA coefficients'+
'\n   @return {ee.Image} Single band image containing computed LDA score per pixel'+
'\n norm_pdf(x, u, s): compute density of normal distribution'+
'\n   @param {Numeric} x observed value(s) for variable x'+
'\n   @param {Numeric} u mean value(s) for variable x'+
'\n   @param {Numeric} s standard deviation(s) of variable x'+
'\n   @return {ee.Image} Single band image containing computed LDA score per pixel';
