function LStoS2(img, sensor){
  if(sensor == 'LS7'){
    return img.rename(['B1', 'B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'PC', 'Cirrus', 'TIR1', 'TIR2', 'QA']);
  }else if(sensor == 'LS8'){
    return img.rename(['B1', 'B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'T1', 'T2', 'sr_aerosol', 'pixel_qa', 'radsat_qa']);
  }
  }

exports.LStoS2 = LStoS2;

// This helper function returns a list of new band names.
function getNewBandNames(suffix, bands) {
  //var seq = ee.List.sequence(1, bands.length());
  //return seq.map(function(b) {
  return bands.map(function(band){
    return ee.String(band).cat(suffix);
  });
 //   return ee.String(prefix).cat(ee.Number(b).int());
  //});
}

//Function to add date info to an image.
//param img: image
//Returns image with ['date'] band holding date of image collection
function addDate(img){
  var dated = img.addBands(img.metadata('system:time_start', 'date'));
  return dated;
}

//Function converting multiband image into single band image of LDA scores
//param img: multiband image (ee.Image)
//param int: intercept parameter from LDA analysis (numeric)
//param xbands: string list of n band names (list)
//param coefficients: numeric list of length n containing LDA coefficients (list)
function ldaScore(img, int, xbands, coefficients){
  var bands = img.select(xbands);

  var coeffs = ee.Dictionary.fromLists(xbands, coefficients).toImage(xbands);

  var score = bands.multiply(coeffs).addBands(ee.Image(int)).reduce(ee.Reducer.sum());

  return score;
  }

exports.ldaScore = ldaScore;

//Function calulating CDF probability of chi-square statistic
//param chi: single band image of scores from standard chi-squared distribution
//param df: degrees of freedom (usually # of bands)
//Return: single band image of probabilities
/*
function chi_p(chi, df){
  var shape = ee.Image(df).divide(2);
  var denom = shape.gamma();
  var num = shape.gammainc(chi.divide(2));
  return (num.divide(denom)).rename(['p']);
}
*/

function chi_p(chi, df){
var cdf = ee.Image(chi.divide(2)).gammainc(ee.Number(df).divide(2));
return cdf.rename(['p']);
}

function gamma_p(stat, df){
  var shape = ee.Image(1);
  var scale = ee.Image(df);
  var denom = shape.gamma();
  var num = shape.gammainc(stat.divide(scale));
  return num.divide(denom).rename(['p']);
}

//Function calulating normal probability of z-score
// using normal CDF approximation: p = 1/1+exp(-1.65451*z)
//param z: single band image of z-scores from standard normal distribution
//Return: single band image of probabilities
function norm_p(z){
  return ee.Image.constant(1).subtract(z.multiply(-1.65451).exp().add(1).pow(-1));
}

//Function to calculate normalized difference metrics and add bands to collections
// param img: multispectral image
// param NIR: band name corresponding to near-infrared spectrum reflectance
// param R: band name corresponding to red spectrum reflectance
// ...
//Return: image with bands ['ndvi', 'nbr', 'ndsi']
function ND(img, NIR, R, G, SWIR1, SWIR2){
  var NBR = img.normalizedDifference([NIR, SWIR2]).rename(["nbr"]);
  var NDSI = img.normalizedDifference([G, SWIR1]).divide(img.select([NIR])).rename(["ndsi"]);
  //var NNDSI = norm(NDSI, ["ndsi"], region).rename(["nndsi"]);
  var NDWI = img.normalizedDifference([G, NIR]).rename(['ndwi']);
  var NDVI = img.normalizedDifference([NIR, R]).rename(["ndvi"]);
  //var NNDVI = norm(NDVI, ["ndvi"], region).rename(["nndvi"]);

  return img.addBands(NDVI)
  .addBands(NDSI)
  .addBands(NBR)
  .addBands(NDWI);
  //.addBands(NNDVI)
  //;
}

//CHANGE METRICS
//Function calculating 'change vector'
//Return: image with ['cv'] band
//param b: before image
//param a: after image
//param bnds: list of band names (typically R, G, B, NIR, SWIR1, SWIR2)
function RCV(b, a, bnds, aoi){
  var diff = b.select(bnds).subtract(a.select(bnds));
  var maxab = b.select(bnds).max(a.select(bnds)).pow(2);
  var stat = diff.divide(maxab);
  var diff_stats = diff.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    geometry: aoi,
    scale: 100,
    maxPixels: 1e13,
    tileScale: 6
  });

  var means = bnds.map(function(band){
    return ee.String(band).cat('_mean');
  });

  var sds = bnds.map(function(band){
    return ee.String(band).cat('_stdDev');
  })
  var mean = diff_stats.toImage(means);
  var sd = diff_stats.toImage(sds);
  var cv = diff.subtract(mean).divide(sd).pow(2).reduce(ee.Reducer.sum()).rename(['cv']);
  var rcv = stat.subtract(mean.divide(maxab)).divide(sd.divide(maxab)).reduce(ee.Reducer.sum()).rename(['rcvmax']);
  return cv.addBands(rcv);
}

function CV(b, a, bnds, aoi){
  var diff = b.select(bnds).subtract(a.select(bnds));
  var diff_sd = diff.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 6
  }).toImage(bnds);
  var diff_mn = diff.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 6
  }).toImage(bnds);
return(
  diff.subtract(diff_mn).divide(diff_sd).pow(2)
.reduce(ee.Reducer.sum())
.rename(['cv'])
);
}

exports.CV = CV;

//Function calculating 'relative change vector maximum'
//Return: image with ['rcvmax'] band
//v b: before image
//param a: after image
//param bnds: list of band names (typically R, G, B, NIR, SWIR1, SWIR2)
function rcvmax(b, a, bnds, aoi){
  var diff = b.select(bnds).subtract(a.select(bnds));
  var maxab = b.select(bnds).max(a.select(bnds)).pow(2);
  var stat = diff.divide(maxab);
  var diff_sd = diff.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 6
  }).toImage(bnds).divide(maxab);
  var diff_mn = diff.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 6
  }).toImage(bnds).divide(maxab);

  return(
    //stat.reduce(ee.Reducer.sum()).rename(['rcvmax'])
    stat.subtract(diff_mn).divide(diff_sd).reduce(ee.Reducer.sum()).rename(['rcvmax'])
    );
}

//Function calculating changes in common band values between two images
//Return: image with bands equal to common bands in b and a
//param b: before image
//param a: after image
//param bnds: list of band names to calculate differences between before and after
function d(b, a, bnds){
  return(
    b.select(bnds).subtract(a.select(bnds))
    );
}

//Function calculating z-score and p-value for all bands in a change image.
//param change: multiband image representing change in pixels values btwn time a and b
//param aoi: geometry restricting analysis
//param scl: numeric scale at which to sample image to calculate global mean/sd
//Return: image with 2(n) bands of original image storing z-scores and p-values
function calc_zp(change, aoi, scl){
  var norm = change.select(['ndvi', 'ndsi', 'nbr', 'ndwi', 'rcvmax']);
  var chi = change.select(['cv']).rename(['cv_z']);
  var norm_bands = norm.bandNames();
  var chi_bands = ['cv'];
  var cat_mean = norm_bands.map(function(string){
    string = ee.String(string).cat('_mode');
    return string;
  });
  var cat_sd = norm_bands.map(function(string){
    string = ee.String(string).cat('_stdDev');
    return string;
  });
  var cat_p = norm_bands.map(function(string){
    string = ee.String(string).cat('_p');
    return string;
  });
  var cat_z = norm_bands.map(function(string){
    string = ee.String(string).cat('_z');
    return string;
  });
  var stats = change.select(norm_bands).reduceRegion({
    reducer: ee.Reducer.mode().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    //kernel: ee.Kernel.fixed(150, 150)
    geometry: aoi,
    scale: scl,
    maxPixels: 1e13,
    tileScale: 6
  });//.rename(norm_bands, cat_mean);
  /*
  var sd = change.select(norm_bands).reduceRegion({
    reducer: ee.Reducer.stdDev(),
    //kernel: ee.Kernel.fixed(150, 150)
    geometry: aoi,
    scale: scl,
    maxPixels: 100000000000,
    tileScale: 6
  }).rename(norm_bands, cat_sd);
  */
  //var stats = mean.combine(sd);
  var img_z = norm.subtract(stats.toImage(cat_mean))
  .divide(stats.toImage(cat_sd)).rename(cat_z);
  var np = norm_p(img_z.abs()).multiply(2).rename(cat_p);
  var cp = chi_p(chi, 6).multiply(-1).add(1).rename(['cv_p']);
  return chi.addBands(img_z).addBands(np).addBands(cp);
}

//Function to iteratively reweight pixels of a chane image as per Nielsen 2007.
// Currently takes a number of iterations as argument, but want to
// update to use a single response metric threshold (e.g. total change in p)
//param change: multiband image representing change in pixels values btwn time a and b
//param niter: number of iterations to perform
//Return: image with 2(n) original image bands storing z-score and p-values from
// the final iteration.
function iw(change, aoi, niter){
  var bands = change.bandNames();
  var cat_p = bands.map(function(string){
    string = ee.String(string).cat('_p');
    return string;
  });
  var net = 1;
  var zs = calc_zp(change, aoi, 100);

  while (net <= niter){
    var dp = zs.select(cat_p).max(0.001).multiply(change).rename(bands);
    zs = calc_zp(dp, aoi, 100);
    net++;
  }
  return zs;
}
exports.iw = iw;

//Function to add area to feature attributes.  Will be replaced by shape metrics function
//param ft: feature
function sz(ft){var area = ft.area(5);
  return ft.set({'area':area});
}

//Function to run complete iteratively reweighted national landcover change analysis.
//param aoi:
//param doi: date after which change will be detected
//param cvz: z-score threshold on change vector metric
//param rcvz: z-score threshold on relative change vector metric
//param ndviz: " " ndvi
//param ndsiz: " " ndsi
//param nbrz: " " nbr
//param size: minimum size (ac) of change to be reported
exports.runIW = function(before, after, aoi, ag){
  var demMask = DEM.select('elevation').lte(3500);
  var agMask = CDL.select('cultivated').eq(1);
  //var projdate = ee.Date(doi);
  //var today = ee.Date(Date.now());
  //var prior = ee.Date.fromYMD(projdate.get('year').subtract(1), projdate.get('month'), projdate.get('day'));
  var rgbn = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
  var time1 = after;//imgCol.filterDate(projdate, today).filterBounds(aoi);
  var time2 = before;//imgCol.filterDate(prior, projdate).filterBounds(aoi);

  var recent = time1.median().clip(aoi);//.qualityMosaic('date').clip(aoi);
  var past = time2.median().clip(aoi);
  //recent = recent.register(past, 100);
  /*
  if (ee.List(before.get('tags')).contains('l8sr')){
    recent = LStoS2(recent, "LS8");
    past = LStoS2(past, "LS8");
  }
  */

  recent = ee.Image(
    ee.Algorithms.If(ag == 'yes', recent.updateMask(agMask.and(demMask)), recent.updateMask(demMask))
    );
  past = ee.Image(
    ee.Algorithms.If(ag == 'yes', past.updateMask(agMask.and(demMask)), past.updateMask(demMask))
    );

  var now = ND(recent, 'B8', 'B4', 'B3', 'B11', 'B12');
  var old = ND(past, 'B8', 'B4', 'B3', 'B11', 'B12');
  //Map.addLayer(old, {bands:['B4', 'B3', 'B2'], min:250, max:3000, gamma:1}, 'old');
  //Map.addLayer(now, {bands:['B4', 'B3', 'B2'], min:250, max:3000, gamma:1}, 'now');
//CREATE IMAGE WITH BANDS FOR CHANGE METRICS CV, RCV, NDVI, NBR, NDSI
//calculate cv and rcvmax from before and after images
  //var cv = RCV(old, now, rgbn, aoi);
  var cv = CV(old, now, rgbn, aoi);
//calculate rcv from before and after images
  var rcv = rcvmax(old, now, rgbn, aoi);

//calculate combined normalized difference metrics from before and after images
  var diff = d(old, now, ['ndvi', 'ndsi', 'ndwi', 'nbr']);

//combine cv, rcv, and normalized difference images into single image
  var change = cv.addBands(diff).addBands(rcv);

  //var zchange = zp(change, aoi, 30);

  var iwchange = iw(change, aoi, 10);

  return iwchange;
};

function changePolys(size, intercept, lda, cvz, rcvz, ndviz, ndsiz, ndwiz, nbrz){
  var ac = ee.Number(size).multiply(4047);
  var ldascore = ldaScore(iwchange, intercept, ['cv_z', 'rcvmax_z', 'ndvi_z', 'ndsi_z', 'ndwi_z', 'nbr_z'], [cvz, rcvz, ndviz, ndsiz, ndwiz, nbrz]);
  //NOTE: This is where we will utilizethe 'ldaScore' function to convert multiple
  //thresholds to a single score
  var selected = ldascore.lte(lda).focal_min(1, 'square', 'pixels').focal_max(1, 'square', 'pixels');

  var polys = selected.reduceToVectors({
      geometry: aoi,
      scale: 10,
      eightConnected: true,
      maxPixels: 1e13
      });

  //NOTE: This is where we will apply the shape metrics in the future...for now just size
  //var change_poly_met = change_poly.filter(ee.Filter.eq('label', 1)).map(metrics);
  //var change_polys = polys.filterMetadata("label", "equals", 1).map(metrics);
  var change_polys = polys.filter(ee.Filter.eq("label", 1)).map(sz);

  var big_polys = change_polys.filter(ee.Filter.gte("area", ac));
  //var count = big_polys.aggregate_count('label');

  //var indicator = (count.getInfo() >=1 )? true : false;

  return big_polys;
}

exports.doc = 'IW contains functions for running the iteratively re-weighted land cover change detection algorithm'+
'\n runIW(before, after, aoi): calculates standardized change metrics (CV, RCVmax, NDVI, NDSI, NDWI, NBR)'+
'\n   @param {ee.ImageCollection} before Image Collection for encompassing historic conditions.'+
'\n   @param {ee.ImageCollection} after Image Collection in which to detect change.'+
'\n   @param {ee.Geometry} aoi Bounding geometry for change detection.'+
'\n   @param {String} sensor indicator if the image colleciton is Landsat or Sentinel. Default Sentinel'+
'\n   @return {ee.Image} 12-band image containing z-score and p-values for 5 change metrics.'+
'\n LStoS2(img): renames bands of LS8 Tier 1 TOA reflectance to match those of S2'+
'\n   @param {ee.Image} imgCol Image Collection for change detection.'+
'\n   @return {ee.Image} LS8 TOA image with renamed bands';
