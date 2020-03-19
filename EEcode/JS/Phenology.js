function norm_pdf(x, u, s){
  x = ee.Image(x);
  u = ee.Image(u);
  s = ee.Image(s);
  var pow = x.subtract(u).pow(2).divide(s.pow(2).multiply(2)).multiply(-1);
  var first = s.pow(2).multiply(Math.PI).multiply(2).sqrt().pow(-1);
  return pow.exp().multiply(first);
}

var addHarmonics = function(img, freq){
  var date = ee.Date(img.get('system:time_start'));
  var time= date.difference(ee.Date('2000-01-01'), 'year').multiply(freq);
  var sin = time.multiply(2 * Math.PI).sin();
  var cos = time.multiply(2 * Math.PI).cos();
  var independent = ee.Image([1, time, sin, cos]).double().rename(['offset', 'time', 'sin', 'cos']);
  return img.addBands(independent);
};

exports.addHarmonics = addHarmonics;

var harmonic_predict = function(img, model, responses){
  var xvars = img.select(['offset', 'time', 'sin', 'cos']).toArray().toArray(1);
  var coeffs = model.select('coefficients').toArray().matrixTranspose();
  var predicted = coeffs.matrixMultiply(xvars).arrayProject([0]).arrayFlatten([responses]);
  var diff = predicted.subtract(img.select(responses)).pow(2);
  return(diff);
};

exports.harmonic_predict = function(img, model, responses){
  var xvars = img.select(['offset', 'time', 'sin', 'cos']).toArray().toArray(1);
  var coeffs = model.select('coefficients').toArray().matrixTranspose();
  var predicted = coeffs.matrixMultiply(xvars).arrayProject([0]).arrayFlatten([responses]);
  var diff = predicted.subtract(img.select(responses)).pow(2);
  return(diff);
};

exports.outliers = function(imgCol, mon1, mon2, responses, freq){
  var today = ee.Date(Date.now());
  imgCol = imgCol.map(function(img){
    return addHarmonics(img, freq);
    });
  var historic = imgCol.filterDate(imgCol.aggregate_min('system:time_start'), mon1);
  var monitor = imgCol.filterDate(mon1, mon2);
  var vars = ee.List(['offset', 'time', 'sin', 'cos']).cat(ee.List(responses));
  var model = historic.select(vars)
  .reduce(ee.Reducer.linearRegression({
    numX: 4,
    numY: ee.List(responses).length(),
    }));
  var residuals = model.select('residuals')
  .arrayProject([0])
  .arrayFlatten([responses]);
  var diffs = monitor.map(function(img){
    return harmonic_predict(img, model, responses);
  });
  return diffs.mean().sqrt().divide(residuals);
};

function mvKern(imgCol, date1, date2, bands, ht){
  var names = bands.map(function(string){
    return ee.String(string).cat('_p');
  });

  var dated = imgCol.map(function(img){
    var date = ee.Date(img.get('system:time_start'));
    var day = ee.Image(date.getRelative('day', 'year')).rename('day');
    return img.addBands(day);
  });

  var test = dated.filterDate(date1, date2);
  var train = dated.filterDate('2000-01-01', date1);
  //Define bandwidth for bands base don rule of thumb estimator
  var n = train.select(bands).reduce(ee.Reducer.count());
  var hb = train.select(bands).reduce(ee.Reducer.stdDev()).multiply(1.06).multiply(n.pow(-0.2));

  var ps = test.map(function(testimg){
    var ks = train.map(function(trainimg){
      var timep = norm_pdf(testimg.select('day'), trainimg.select('day'), ht);
      var bp = norm_pdf(testimg.select(bands), trainimg.select(bands), hb);
      var k = timep.multiply(bp).divide(ht.add(hb));
      return k.rename(names);
    });
    return ks.sum().divide(train.size());
  });

  return ps;
}
exports.mvKern = mvKern;

exports.doc = 'Phenology is a collection of functions to calculate harmonic trends across time series.'+
  '\n addHarmonics(img, freq):'+
  '\n   @param {ee.Image} img The image to which variables "time", "cos", "sin", and "offset" will be added'+
  '\n   @param {Int} freq The number of annual cycles to model'+
  '\n   @return {ee.Image} The original image with variables "time", "cos", "sin", and "offset"'+
  '\n'+
  '\n harmonic_predict(img, model, responses):'+
  '\n   @param {ee.Image} img The image for which predictions based on a harmonic model will be estimaed'+
  '\n   @param {ee.Image} model Result from linear regression reducer containing model coefficients and residuals'+
  '\n   @param {ee.List} repsonses band names for the response variables'+
  '\n   @return {ee.Image} squared difference between observed response and harmonic model predictions.'+
  '\n outliers(imgCol, mon1, mon2, responses, freq):'+
  '\n   @param {ee.ImageCollection} imgCol The image collection to which a harmonic model will be fit and estimated'+
  '\n   @param {ee.String} mon1 Start date dividing imgCol into model fitting and monitoring images'+
  '\n   @param {ee.String} mon2 End date of monitoring period'+
  '\n   @param {ee.List} repsonses band names for the response variables'+
  '\n   @param {Int} freq The number of annual cycles to model'+
  '\n   @return {ee.Image} Image RMSE of pixel response values for monitoring period, standardized by residuals.'+
  '\n mvKern(imgCol, date1, date2, bands, ht): implements a bivariate kernel density estimator'+
  '\n   @param {ee.ImageCollection} imgCol: The image collection from which kernel density function is estimated'+
  '\n   @param {ee.String} mon1: Start date dividing imgCol into model fitting and monitoring images'+
  '\n   @param {ee.String} mon2: End date of monitoring period'+
  '\n   @param {ee.List} bands: band names for the response variables'+
  '\n   @param {ee.Image} ht: Temporal bandwidth (days).  Default = 30'+
  '\n   @return {ee.Image} Image containing p-value bands for each input band.';

