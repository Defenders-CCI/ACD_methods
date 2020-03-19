
//Function calulating CDF probability of chi-square statistic
//param chi: single band image of scores from standard chi-squared distribution
//param df: degrees of freedom (usually # of bands)
//Return: single band image of probabilities
function chi_p(chi, df){
  var cdf = ee.Image(chi.divide(2)).gammainc(ee.Number(df).divide(2));
  return cdf.rename(['p']);
}

//function paramaterizing a gamma distribution from data, and returning p-values of observations
//param chi: image of chi-square statistics, or other values following a gamma distribution
//param aoi: geometry object defining area of analisys
//param scl: scale (m) at which reduction operations are performed
//return: image with single band ['p'] containing p-values of original image values
function gamma(chi, aoi, scl){
  var mode = chi.reduceRegion({
    reducer: ee.Reducer.mode(),
    geometry: aoi,
    scale: scl,
    maxPixels: 1e13
  }).toImage();
  var sd = chi.reduceRegion({
    reducer: ee.Reducer.sampleStdDev(),
    geometry: aoi,
    scale: scl,
    maxPixels: 1e13
  }).toImage();
  var rate = (((((sd.pow(2)).multiply(4)).add(mode.pow(2))).sqrt()).add(mode)).divide(sd.pow(2).multiply(2));
  //shape calculated from observed mode
  var shape = mode.multiply(rate).add(1);
  //shape calculated assuming mode = 0 (i.e. no change)
  //var shape = ee.Image.constant(1);
  var p = //ee.Image.constant(1).subtract(
    shape.gammainc(chi.multiply(rate)).rename(['p']);
  return p;
}

//Function computing the chi square statistic and p-value sum of squared band values
//param img: multi-band image
//param aoi: geometry object defining area of analisys
//param scl: scale (m) at which reduction operations are performed
//Return: image with bands ['chi', 'p'] equal to chi-square statistic and right-tail p-value
function chisq(img, aoi, scl){
  var k = img.bandNames().length();
  var k2 = k.divide(2);
  //return dictionary of standard deviations for each input band
  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    geometry: aoi,
    scale: 100,
    maxPixels: 1e13,
    tileScale: 4
  });
  //var sd = img.reduceRegion({
  //  reducer: ee.Reducer.stdDev(),
  //  geometry: aoi,
  //  scale: scl,
  //  maxPixels: 100000000000}).toImage(img.bandNames());
  //var mean = img.reduceRegion({
  //  reducer: ee.Reducer.mean(),
  //  geometry: aoi,
  //  scale: scl,
  //  maxPixels: 100000000000}).toImage(img.bandNames());
  //calculte chi-square statistic as sum(i:p) of (MADij/stdDev(i))^2
  var means = img.bandNames().map(function(band){
    return ee.String(band).cat('_mean');
  });

  var sds = img.bandNames().map(function(band){
    return ee.String(band).cat('_stdDev');
  });
  var mean = stats.toImage(means);
  var sd = stats.toImage(sds);
  var chi = img.subtract(mean).divide(sd).pow(2).reduce(ee.Reducer.sum()).rename(['chi']);
  //var chi = (((img.subtract(mean)).divide(sd)).pow(2)).reduce(ee.Reducer.sum()).rename(['chi']);
  //var p = ee.Image.constant(1).subtract(chi.divide(2).gammainc(k2)).rename(['p']);
  var p = chi_p(chi, k).multiply(-1).add(1).max(0.001);
  return chi.addBands(p);
}

//Function to compute correlation matrix between two sets of variables
//param set1: image with n bands equal to first set of variables
//param set2: image with m bands equal to second set of variables
//param aoi: geometry object defining area of analisys
//Return: n x m array containing correlation coefficients between variables in set1 and set2
function corrmat(set1, set2, aoi, weights){
  var emptylist = ee.List([]);
  var xbands = set1.bandNames();
  var ybands = set2.bandNames();
  var xs = xbands.iterate(function(xband, list){
    var ys = ybands.iterate(function(yband, list){
      var selected = set1.select([xband]).addBands(set2.select([yband]));

      var cor = selected.reduceRegion({
        reducer: ee.Reducer.pearsonsCorrelation(),
        geometry: aoi,
        scale: 100,
        maxPixels: 1e13,
        tileScale: 4
        }).get('correlation');

      return ee.List(list).add(cor);
    }, emptylist);

    return ee.List(list).add(ys);
  }, emptylist);
  return ee.Array(xs);
}

//Function to compute correlation matrix between two sets of variables
//param {ee.Image} set1: N-band image equal to first set of variables
//param {ee.Image} set2: M-band image equal to second set of variables
//param {ee.Geometry} aoi: Area of analisys
//param {ee.Image} weights Single band image containing per pixel weights
//Return: n x m array containing correlation coefficients between variables in set1 and set2
function wCorrmat(set1, set2, aoi, weights){
  weights = weights.clip(aoi);

  var b1 = weights.bandNames().get(0)
  var n = set1.bandNames().length();

  var img = set1.addBands(set2);

  var weightsImage = img.multiply(ee.Image.constant(0).add(weights))

  var N = img.bandNames().length();

  var xbands = ee.List(img.bandNames().slice(0, n));

  var ybands = ee.List(img.bandNames().slice(n, N));

  var sumWeights = ee.Number(weights.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: aoi,
    scale: 10,
    tileScale: 4,
    maxPixels: 1e13}).get(b1));

  var nPixels = ee.Number(weights.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: aoi,
    scale: 10,
    tileScale: 4,
    maxPixels: 1e13}).get(b1));

  var means = img.addBands(weightsImage)
    .reduceRegion({
      reducer: ee.Reducer.mean().repeat(N).splitWeights(),
      geometry: aoi,
      scale: 100,
      maxPixels: 1e13
    }).toArray().project([1]);

  var centered = img.toArray().subtract(means);

  var weighted = img.multiply(weights.sqrt());

  var xs = xbands.map(function(xbnd){
    var ys = ybands.map(function(ybnd){

      var selected = weighted.select([xbnd, ybnd])//.addBands(weighted.select([ybnd]));

      var cor = selected.reduceRegion({
        reducer: ee.Reducer.pearsonsCorrelation(),
        geometry: aoi,
        scale: 100,
        maxPixels: 1e13,
        tileScale: 4
        }).get('correlation');

      return cor;
        });
    return ys;
  });
  var corw = ee.Array(xs).multiply(nPixels).divide(sumWeights);
  return corw;
}

//Function to perform canonical correlation analysis using singular value deocmposition
//on a correlation matrix between two variable sets.
//param cor: ee.Array n x m correlation matrix
//param before: image with n bands
//param after: image with m bands
//param length: minimum of n & m
//Return image with j = min(n, m) bands ["V1", "V2",..."Vj", "p"], corresponding to differences aX - bY across j canonical variates
function canon(cor, before, after, length){
  var labels = ee.List.sequence(1, length).map(function(item){
    return ee.String("V").cat(ee.Number(item).toInt().format());
  });
  var decomp = cor.matrixSingularValueDecomposition();
  //create n * min(n, m) matrix of canonical covariates a
  var U = ee.Array(decomp.get("U"));
  //create m * min(n, m) matrix of canonical covariates b
  var V = ee.Array(decomp.get("V"));
  //get diagonal elements of SVD equivalent to CCA correlations
  var S = ee.Array(decomp.get("S")).matrixDiagonal();
  //turn images into array images with 1 x nbands matrices at each pixel
  var before2D = before.toArray().toArray(1).arrayTranspose();
  var after2D = after.toArray().toArray(1).arrayTranspose();
  var a = before2D.matrixMultiply(ee.Image(U)).arrayProject([1]).arrayFlatten([labels]);
  var b = after2D.matrixMultiply(ee.Image(V)).arrayProject([1]).arrayFlatten([labels]);
  return ee.Dictionary({'img': a.subtract(b), 'rhos': S});
}

//Function to conduct MAD transformation between two images
// using single value decomposition on correlation matrix
//param before: before image
//param after: after image
//param aoi: geometry object defining area of analysis
//Return image with j + 2 bands ["V1", "V2", ..."Vj", "chi", "p"]
function mad(before, after, aoi){
  var length = before.bandNames().length().min(after.bandNames().length());
  var corMat = corrmat(before, after, aoi);
  var cca = canon(corMat, before, after, length);
  var cvs = ee.Image(cca.get('img'));
  var ccachi = chisq(cvs, aoi, 100);
  var rhos = cca.get('rhos');
  return cvs.addBands(ccachi);
  //return ee.Dictionary({'img': cvs.addBands(ccachi), 'rhos':rhos});
}

//Function to conduct weighted MAD transformation between two images
// using single value decomposition on correlation matrix
//param before (ee.Image): before image
//param after (ee.Image): after image
//param aoi (ee.Geometry): area of analysis
//Return ee.Dictionary 'img' j + 2 bands ["V1", "V2", ..."Vj", "chi", "p"]
function madw(before, after, aoi, weights){
  var length = before.bandNames().length().min(after.bandNames().length());
  var corMat = wCorrmat(before, after, aoi, weights);
  var cca = canon(corMat, before, after, length);
  var cvs = ee.Image(cca.get('img'));
  var ccachi = chisq(cvs, aoi, 100);
  var rhos = ee.Array(cca.get('rhos'));
  //net++
  return ee.Dictionary({
    'img': cvs.addBands(ccachi),
    'rhos':rhos,
    'before': before,
    'after': after,
    'aoi': aoi
  });
}

// Function running weighted MAD using previous MAD output.  For use with iterate
// param current (): placeholder for use in GEE iterate
// param prev(ee.Dictionary): output from madw containing
//                            before, after, and mad output images
// returns (ee.Dictionary):
function iw_mad(current, prev){
  prev = ee.Dictionary(prev);
  var aoi = ee.Geometry(prev.get('aoi'));
  var rhos = ee.Array(prev.get('rhos'));
  var before = ee.Image(prev.get('before'));
  var after = ee.Image(prev.get('after'));
  var length = before.bandNames().length().min(after.bandNames().length())
  var output = ee.Image(prev.get('img'));
  var weights = output.select(['p']);
  var corMat = wCorrmat(before, after, aoi, weights);
  var cca = canon(corMat, before, after, length);
  var cvs = ee.Image(cca.get('img'));
  var ccachi = chisq(cvs, aoi, 100);
  var updates = ee.Array(cca.get('rhos'));
  var net = updates.sort()
    .subtract(rhos).abs()
    .reduce(ee.Reducer.max(), [0])
    .get([0,0]);
  var converged = net.lte(0.01)
  var result = ee.Dictionary({
    'img': cvs.addBands(ccachi),
    'rhos': updates,
    'before': before,
    'after': after,
    'aoi': aoi,
    'iter': current
  });
  return ee.Algorithms.If(converged, prev, result)
}

//Function to perform iterative reweighting and MAD transformations
//param niter: number of iterations
//Returns output same as mad()
function mad_iw(before, after, aoi, niter){
  var first = madw(before, after, aoi, ee.Image(1));
  var rhos = ee.Array(first.get('rhos'));
  var output = ee.Image(first.get('img'))
  var weights = output.select(['p']);
  var net = 10;
  var iter = 1;
  while (iter < niter & net > 0.01) {
    var update = madw(before, after, aoi, weights)
    output = ee.Image(update.get(['img']))
    weights = output.select(['p'])
    var updates = ee.Array(update.get('rhos'))
    net = updates.sort().subtract(rhos).abs()
    .reduce(ee.Reducer.max(), [0]).get([0,0]);
    rhos = updates;
    iter++
    }
  return update.set('iter', iter);
}

//Function running MAD analysis by converting a covariance matrix to
//correlation matrix and performing singular value decomposition
//param before: before image
//param after: after image
//param aoi: geometry defining area of analysis
//param scale: scale of spatial reductions
function MAD(before, after, aoi, scale){
  var beforeArray = before.toArray()
  var afterArray = after.toArray()
  var img = before.addBands(after).toArray()
  var N = before.bandNames().size()
  var labels = ee.List.sequence(1, N).map(function(item){
    return ee.String("V").cat(ee.Number(item).toInt().format());
  });
  var covar = img.reduceRegion({
      reducer: ee.Reducer.covariance(),
      geometry: aoi,
      scale : scale,
      maxPixels: 1e13
  }).toArray().slice(0, 0, N).slice(1, N, N.multiply(2))

  var diag = covar.matrixDiagonal().matrixToDiag()
  var inv_diag = diag.matrixInverse()
  var cormat = inv_diag.matrixMultiply(covar).matrixMultiply(inv_diag)
  var decomp = cormat.matrixSingularValueDecomposition()
  //create n * min(n, m) matrix of canonical covariates a
  var U = ee.Array(decomp.get("U"))
  //create m * min(n, m) matrix of canonical covariates b
  var V = ee.Array(decomp.get("V"))
  var S = ee.Array(decomp.get("S"))

  S = S.matrixDiagonal().project([0])

  var before2D = beforeArray.toArray(1).arrayTranspose();
  var after2D = afterArray.toArray(1).arrayTranspose();
  var a = before2D.matrixMultiply(ee.Image(U)).arrayProject([1]).arrayFlatten([labels]);
  var b = after2D.matrixMultiply(ee.Image(V)).arrayProject([1]).arrayFlatten([labels]);
  var diff = a.subtract(b)
  return(diff)
}

//Function running MAD analysis by converting a covariance matrix to
//correlation matrix and performing singular value decomposition using weighted inputs
//param before: before image
//param after: after image
//param weights: single band image of pixel weights
//param aoi: geometry defining area of analysis
//param scale: scale of spatial reductions
function wMAD(before, after, weights, aoi, scale){
  var beforeArray = before.toArray()
  var afterArray = after.toArray()
  var img = before.addBands(after)
  var N = img.bandNames().size()
  var imgweights = img.multiply(0).add(weights)
  var means = img.addBands(imgweights).reduceRegion({
    reducer: ee.Reducer.mean().repeat(N).splitWeights(),
    geometry: aoi,
    scale: scale,
    tileScale: 12,
    maxPixels: 1e13
  })
  .toArray()
  .project([1])

  var sumweights = weights.reduceRegion({
    reducer: ee.Reducer.count().combine(ee.Reducer.sum(), '', true),
    geometry: aoi,
    scale: 10,
    tileScale: 12,
    maxPixels: 1e13
  })

  print(sumweights)
  var sum = sumweights.get('p_sum')
  var npix = sumweights.get('p_count')

  var centered = img.toArray().subtract(means)

  var covar = centered.multiply(weights.sqrt()).reduceRegion({
      reducer: ee.Reducer.centeredCovariance(),
      geometry: aoi,
      scale : scale,
      maxPixels: 1e13
  }).toArray().slice(0, 0, N.divide(2)).slice(1, N.divide(2), N)

  var wcovar = covar.multiply(npix).divide(sum)

  var diag = wcovar.matrixDiagonal().matrixToDiag()
  var inv_diag = diag.matrixInverse()
  var cormat = inv_diag.matrixMultiply(covar).matrixMultiply(inv_diag)
  var decomp = cormat.matrixSingularValueDecomposition()
  //create n * min(n, m) matrix of canonical covariates a
  var U = ee.Array(decomp.get("U"))
  //create m * min(n, m) matrix of canonical covariates b
  var V = ee.Array(decomp.get("V"))
  var S = ee.Array(decomp.get("S"))

  S = S.matrixDiagonal().project([0])

  var before2D = beforeArray.toArray(1).arrayTranspose();
  var after2D = afterArray.toArray(1).arrayTranspose();
  var a = before2D.matrixMultiply(ee.Image(U)).arrayProject([1]).arrayFlatten([['V1', 'V2', 'V3', 'V4']]);
  var b = after2D.matrixMultiply(ee.Image(V)).arrayProject([1]).arrayFlatten([['V1', 'V2', 'V3', 'V4']]);
  var diff = a.subtract(b)
  return(diff)
}


function runMAD(before, after, aoi, ag){
  var DEM = ee.Image("USGS/SRTMGL1_003");
  var CDL = ee.Image("USDA/NASS/CDL/2017");
  var demMask = DEM.select('elevation').lte(3500);
  var agMask = CDL.select('cultivated').eq(1);
  before = before.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']);
  after = after.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']);

  var recent = after.median().clip(aoi);
  var past = before.median().clip(aoi);

  recent = ee.Image(
    ee.Algorithms.If(ag == 'yes', recent.updateMask(agMask.and(demMask)), recent.updateMask(demMask))
    );
  past = ee.Image(
    ee.Algorithms.If(ag == 'yes', past.updateMask(agMask.and(demMask)), past.updateMask(demMask))
    );

  var mad = ee.List.sequence(1, 10).iterate(iw_mad, madw(past, recent, aoi, ee.Image(1)));
  var madout = ee.Image(ee.Dictionary(mad).get('img'))
  return madout;
}

exports.runMAD = runMAD;
exports.mad = mad;
exports.madw = madw;
exports.wCorrmat = wCorrmat
exports.iw_mad = iw_mad;
exports.mad_iw = mad_iw;

exports.doc = 'MAD contains functions for running the multivariate alteration detection algorithm'+
'\n runMAD(before, after, aoi, ag): calculates canonical correspondence variates and chi-squre score'+
'\n   @param {ee.ImageCollection} before Image Collection encompassing historic conditions.'+
'\n   @param {ee.ImageCollection} after Image Collection in which to detect change.'+
'\n   @param {ee.Geometry} aoi Bounding geometry for change detection.'+
'\n   @param {String} ag indicator if agriculture should be masked {null, "yes"}'+
'\n   @return {ee.Image} image containing chi-score and p-values.';

//EXAMPLE
/*
LS8 = LS8.filterBounds(geometry).filterMetadata('CLOUD_COVER', 'less_than', 10);

var before = LS8.filterDate('2016-06-01', '2016-09-01');
var after = LS8. filterDate('2017-06-01', '2017-09-01');


before = before.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
after = after.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7']);

var recent = after.median().clip(geometry);
var past = before.median().clip(geometry);

Map.addLayer(past, {bands:['B4', 'B3', 'B2'], min:300, max:3000, gamma:1}, 'before');
Map.addLayer(recent, {bands:['B4', 'B3', 'B2'], min:300, max:3000, gamma:1}, 'after');

var cca = mad(past, recent, geometry, ee.Image(1));
Map.addLayer(cca, {}, 'mad');
//var iw = mad_iw(before, after, geometry, 20);
*/






