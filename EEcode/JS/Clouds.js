var Stats = require('users/defendersofwildlifeGIS/DOW_ST_Collab:Stats');
// Convert processed sentinel toa reflectance to raw values, and extract azimuth/zenith metadata
function sentinel2toa(img) {
  var toa = img.select(['B1','B2','B3','B4', 'B5', 'B6','B7','B8', 'B8A','B9','B10', 'B11','B12'])
     .divide(10000)
     .set('solar_azimuth',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
     .set('solar_zenith', img.get('MEAN_SOLAR_ZENITH_ANGLE'))
     .set('viewing_azimuth', img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8'))
     .set('viewing_zenith', img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8'))
     .set('CLOUDY_PIXEL_PERCENTAGE', img.get('CLOUDY_PIXEL_PERCENTAGE'));
  return toa.addBands(img.select(['QA60']));
}

exports.sentinel2toa = sentinel2toa;

function rescale(img, exp, thresholds) {
    return img.expression(exp, {img: img})
        .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
  }

//Algorithm to compute liklihood of water
//Builds on logic from Google cloudScore algorithm
// developed by Ian Housmann
function waterScore(img){
  img = sentinel2toa(img);
      // Compute several indicators of water and take the minimum of them.
      var score = ee.Image(1.0);

      //Set up some params
      var darkBands = ['B3','B4','B8','B11','B12'];//,'nir','swir1','swir2'];
      var brightBand = 'B2';
      var shadowSumBands = ['B8','B11','B12'];
      //Water tends to be dark
      var sum = img.select(shadowSumBands).reduce(ee.Reducer.sum());
      var sum = rescale(sum,'img',[0.35,0.2]).clamp(0,1)
      score = score.min(sum);

      //It also tends to be relatively bright in the blue band
      var mean = img.select(darkBands).reduce(ee.Reducer.mean());
      var std = img.select(darkBands).reduce(ee.Reducer.stdDev());
      var z = (img.select([brightBand]).subtract(std)).divide(mean);
      z = rescale(z,'img',[0,1]).clamp(0,1);
      score = score.min(z);


      // // Water is at or above freezing
      // score = score.min(rescale(img, 'img.temp', [273, 275]));


      // // Water is nigh in ndsi (aka mndwi)
      var ndsi = img.normalizedDifference(['B3', 'B11']);
      ndsi = rescale(ndsi, 'img', [0.3, 0.8]);


      score = score.min(ndsi);

      return score.clamp(0,1).rename('waterScore');

      }
exports.waterScore = waterScore;

basicQA = function(img){
  var qa = img.select('QA60').int16();

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  var dated = img.updateMask(mask);
  return dated;
};

exports.basicQA = basicQA

// Function to cloud mask from the Fmask band of Landsat 8 SR data.
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = ee.Number(2).pow(3).int();
  var cloudsBitMask = ee.Number(2).pow(4).int();
  var snowBitMask = ee.Number(2).pow(5).int();
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
      .and(qa.bitwiseAnd(snowBitMask).eq(0));

  // Return the masked image, scaled to [0, 1].
  return image.updateMask(mask);//.divide(10000);
}

exports.Fmask = maskL8sr;

function cloudBands(img){
  var ndmi = img.normalizedDifference(['B8','B11']).rename('ndmi');
  var ndsi = img.normalizedDifference(['B3', 'B11']).rename('ndsi');
  var cirrus = img.select(['B1', 'B10']).reduce(ee.Reducer.sum()).rename('cirrus');
  var vis = img.select(['B4', 'B3', 'B2']).reduce(ee.Reducer.sum()).rename('vis');
  return img.addBands(ndmi).addBands(ndsi).addBands(cirrus).addBands(vis);
}

exports.cloudBands = cloudBands;

function darkC (img, R, G, B){
  var R = img.select(R);
  var G = img.select(G);
  var B = img.select(B);
  var maxRB = R.max(B);
  var maxGB = G.max(B);
  var maxRG = R.max(G);
  var C1 = G.divide(maxRB).atan().rename('C1');
  var C2 = R.divide(maxGB).atan().rename('C2');
  var C3 = B.divide(maxRG).atan().rename('C3');
  return img.addBands(C1).addBands(C2).addBands(C3);
}

exports.darkC = darkC;

//Algorithm to compute liklihood of clouds in Sentinel-2 iamagery
//Builds on logic from Google cloudScore algorithm for Landsat
// developed by Ian Housmann ian.housmann@gmail.com
function sentinelCloudScore(img) {
  var im = sentinel2toa(img);
  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1);

  // Clouds are reasonably bright in the blue and cirrus bands.
  score = score.min(rescale(im, 'img.B2', [0.1, 0.5]));
  score = score.min(rescale(im, 'img.B1', [0.1, 0.3]));
  score = score.min(rescale(im, 'img.B1 + img.B10', [0.15, 0.2]));

  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(im, 'img.B4 + img.B3 + img.B2', [0.2, 0.8]));

  //Clouds are moist
  var ndmi = im.normalizedDifference(['B8','B11']);
  score=score.min(rescale(ndmi, 'img', [-0.1, 0.1]));

  // However, clouds are not snow.
  var ndsi = im.normalizedDifference(['B3', 'B11']);
  score=score.min(rescale(ndsi, 'img', [0.8, 0.6]));

  score = score.multiply(100).byte();

  return img.addBands(score.rename('cloudScore'));
}

exports.sentinelCloudScore = sentinelCloudScore;
/*function(img){
  img = sentinel2toa(img);
  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1);

  // Clouds are reasonably bright in the blue and cirrus bands.
  score = score.min(rescale(img, 'img.B2', [0.1, 0.5]));
  score = score.min(rescale(img, 'img.B1', [0.1, 0.3]));
  score = score.min(rescale(img, 'img.B1 + img.B10', [0.15, 0.2]));

  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.B4 + img.B3 + img.B2', [0.2, 0.8]));

  //Clouds are moist
  var ndmi = img.normalizedDifference(['B8','B11']);
  score=score.min(rescale(ndmi, 'img', [-0.1, 0.1]));

  // However, clouds are not snow.
  var ndsi = img.normalizedDifference(['B3', 'B11']);
  score=score.min(rescale(ndsi, 'img', [0.8, 0.6]));

  score = score.multiply(100).byte();

  return img.addBands(score.rename('cloudScore'));
};
*/
///////////////////////////////////////////////////////////

//Function for finding dark outliers in time series
//Masks pixels that are dark, and dark outliers

function TSp(c, bands){
  //Caclulate mean and standard deviation of band values across time series
  var stdDev = c.select(bands).reduce(ee.Reducer.stdDev());
  var mean = c.select(bands).reduce(ee.Reducer.median());

  //For each image in time series..
  c = c.map(function(img){
    //create new band names with '_z' suffix
    var cat_z = bands.map(function(string){
      return ee.String(string).cat('_z');
    });
    //calculate z-score for each specified band in each image
    var z_img = img.select(bands).subtract(mean).divide(stdDev).rename(cat_z);
    return img.addBands(z_img);
  });

  return c;//.select(bandNames);
}

exports.TSp = TSp;
/*
function(c, bands){
  //Caclulate mean and standard deviation of band values across time series
  var stdDev = c.select(bands).reduce(ee.Reducer.stdDev());
  var mean = c.select(bands).reduce(ee.Reducer.median());

  //For each image in time series..
  c = c.map(function(img){
    //create new band names with '_z' suffix
    var cat_z = bands.map(function(string){
      return ee.String(string).cat('_z');
    });
    //calculate z-score for each specified band in each image
    var z_img = img.select(bands).subtract(mean).divide(stdDev).rename(cat_z);
    return img.addBands(z_img);
  });

  return c;//.select(bandNames);
};
*/
function ESAcloud(toa) {
  // author: Nick Clinton

  var qa = toa.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);

  // clear if both flags set to zero.
  var clear = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));

  var cloud = clear.eq(0)

  return cloud
}

function shadowMask(toa,cloud){
  // Author: Gennadii Donchyts
  // License: Apache 2.0

  // solar geometry (radians)
  var azimuth =ee.Number(toa.get('solar_azimuth')).multiply(Math.PI).divide(180.0).add(ee.Number(0.5).multiply(Math.PI));
  var zenith  =ee.Number(0.5).multiply(Math.PI ).subtract(ee.Number(toa.get('solar_zenith')).multiply(Math.PI).divide(180.0));

  // find where cloud shadows should be based on solar geometry
  var nominalScale = cloud.projection().nominalScale();
  var cloudHeights = ee.List.sequence(200,10000,500);
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    var shadowVector = zenith.tan().multiply(cloudHeight);
    var x = azimuth.cos().multiply(shadowVector).divide(nominalScale).round();
    var y = azimuth.sin().multiply(shadowVector).divide(nominalScale).round();
    return cloud.changeProj(cloud.projection(), cloud.projection().translate(x, y));
  });
  var potentialShadow = ee.ImageCollection.fromImages(shadows).max();

  // shadows are not clouds
  potentialShadow = potentialShadow.and(cloud.not());

  // (modified by Sam Murphy) dark pixel detection
  var darkPixels = toa.normalizedDifference(['B2', 'B12']).gt(0.25).rename(['dark_pixels']);

  // shadows are dark
  var shadow = potentialShadow.and(darkPixels).rename('shadows');

  return shadow
}

function cloud_and_shadow_mask(img) {

  var toa = sentinel2toa(img);

  var cloud = sentinelCloudScore(img).select('cloudScore').gte(20);

  var shadow = shadowMask(toa,cloud);

  var mask = cloud.or(shadow).eq(0);

  return img.updateMask(mask);

}

exports.cloud_and_shadow_mask = cloud_and_shadow_mask;

//right now this is taking the difference of cloud_z and shadow_z...way too slow..go back to adding 1s and 0s
//param clouds: cloud mask image with clouds scored as '1', and all other pixels '0'
//param shadows: shadow mask image with shadows scored as '1' and all other pixels '0'
// NOTE: Extent must exceed that of clouds by max of cloudHeights
//Return: Estimate of cloud height in image (number)

function projectShadows (img, aoi, cloudband, cloudthresh, darkband, darkthresh, cloudHeights, yMult){
//function projectShadows(cloudMask, img, darkthresh, cloudHeights, yMult){
  var max = cloudHeights.reduce(ee.Reducer.max());

  if(yMult === undefined || yMult === null){
    yMult = ee.Algorithms.If(ee.Algorithms.IsEqual(img.select([3]).projection(), ee.Projection("EPSG:4326")),1,-1);
  }
  var meanAzimuth = img.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var meanZenith = img.get('MEAN_SOLAR_ZENITH_ANGLE');
  var viewAzimuth = img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8');
  var viewZenith = img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8');

  //Get scale of image for 'B8' will be 10
  var nominalScale = img.select(cloudband).projection().nominalScale();
  print(nominalScale);
  //Find where cloud shadows should be based on solar geometry
  //Convert to radians
  var azR = ee.Number(meanAzimuth).add(180).multiply(Math.PI).divide(180.0);
  var zenR  = ee.Number(meanZenith).multiply(Math.PI).divide(180.0);

  var vazR = ee.Number(viewAzimuth).add(180).multiply(Math.PI).divide(180.0);
  var vzenR = ee.Number(viewZenith).multiply(Math.PI).divide(180.0);

  //Clip an image extending beyond border of target image to match all clouds to their shadows
  //These have the same projection and scale
  var image = img.clip(aoi);
  var cloudMask = image.select(cloudband).gte(cloudthresh).focal_min(2).focal_max(3);
  var buff = aoi.buffer(zenR.tan().multiply(max));
  var buff_image = img.clip(buff);

  //Find dark pixels
  //This image has same projection and scale as 'image' and 'buff_image'
  //var darkPixels = buff_image.select(darkband).reduce(ee.Reducer.sum()).round();
  var darkPixels = buff_image.select(darkband).reduce(ee.Reducer.sum()).lt(darkthresh)
    .focal_min(2).focal_max(3);
  //.gte(1);
  //print(cloudMask.max(darkPixels).bandNames());
  //Find the shadows
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);

    var cloudCastedDistance = vzenR.tan().multiply(cloudHeight);//Distance cloud is projected onto image
    var shadowCastedDistance = zenR.tan().multiply(cloudHeight);//Distance shadow is cast

    var cldx = vazR.sin().multiply(cloudCastedDistance).multiply(yMult).divide(nominalScale);//X distance of cloud
    var cldy = vazR.cos().multiply(cloudCastedDistance).divide(nominalScale);//Y distance of cloud
    //Vertical location of cloud above ground
    var cloudpos = cloudMask.changeProj(cloudMask.select(cloudband).projection(), cloudMask.select(cloudband).projection().translate(cldx, cldy));

    var x = azR.sin().multiply(shadowCastedDistance).divide(nominalScale);//X distance of shadow
    var y = azR.cos().multiply(shadowCastedDistance).divide(nominalScale).multiply(yMult);//Y distance of shadow
    //print(x,y);

    var shift = cloudpos.changeProj(cloudpos.select(cloudband).projection(), cloudpos.select(cloudband).projection().translate(x, y));
    shift = shift.updateMask(cloudMask.not());
    var total = shift.max(darkPixels).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: buff,
      scale: 10
    }).get(cloudband);

    return total;
  });

  //return the position of the maximum overlap score
  var maxhgt = shadows.indexOf(shadows.reduce(ee.Reducer.min()));
  //return the most likely cloud height
  var height_est = cloudHeights.get(maxhgt);

  var darkScale = darkPixels.projection().nominalScale();
  var cloudCastedDistance = vzenR.tan().multiply(height_est);//Distance cloud is projected onto image
  var shadowCastedDistance = zenR.tan().multiply(height_est);//Distance chadow is cast from cloud
  //Reverse the signs of x and y variables to move shadows 'back' towards clouds
  var x = azR.sin().multiply(shadowCastedDistance).divide(darkScale).multiply(yMult);//X distance of shadow
  var y = azR.cos().multiply(shadowCastedDistance).divide(darkScale);//Y distance of shadow
  var cldx = vazR.sin().multiply(cloudCastedDistance).divide(darkScale);//X distance of cloud
  var cldy = vazR.cos().multiply(cloudCastedDistance).divide(darkScale).multiply(yMult);//Y distance of cloud


  //Create shadow mask
  //var shadowMask = img.changeProj(img.select(['cloudScore_z']).projection(), img.select(['cloudScore_z']).projection().translate(x.add(cldx), y.add(cldy)));

  return ee.Dictionary.fromLists(['x', 'y', 'height'], [x.add(cldx), y.add(cldy), height_est]);
}

//Boolean whether projected clouds intersect shadows
//param x: x distance to project clouds to align with shadows
//param y: y distance to project clouds to align with shadows
//param clouds: single-band image indicating clouds[1] vs. no-clouds[0]
//param shadows: single-band image indicating areas of shadow[1] vs. no-shadow[0]
function comboMask(x, y, shadows, clouds){
  //x = ee.Number(x);
  //y = ee.Number(y);
  var xFilter = ee.Filter.intersects({
    leftField: '.geo',
    rightField: '.geo',
    });

  var cloudshift = clouds.changeProj(clouds.projection(), clouds.projection().translate(x, y));

  var cloud_polys = cloudshift.reduceToVectors({
    eightConnected: true,
    labelProperty: 'cloud'
  }).filterMetadata('cloud', 'equals', 1);

  var shadow_polys = shadows.reduceToVectors({
    eightConnected: true,
    labelProperty: 'shadow'
  }).filterMetadata('shadow', 'equals', 1);

  var shadowmask = shadows.changeProj(shadows.projection(), shadows.projection().translate(ee.Number(x).multiply(-1), ee.Number(y).multiply(-1)));

  cloud_polys = ee.Join.saveAll('matches').apply(cloud_polys, shadow_polys, xFilter);

  var keep_clouds = cloud_polys.reduceToImage({
    properties: ['cloud'],
    reducer: ee.Reducer.first()
  });

  var cloudmask = keep_clouds.changeProj(keep_clouds.projection(), keep_clouds.projection().translate(ee.Number(x).multiply(-1), ee.Number(y).multiply(-1)));

  var mask = cloudmask.or(shadowmask).or(shadows);

  return mask;
}

exports.doc = "A collection of funcitons for identifying and masking clouds"+
'\n basicQA(img): mask clouds using Sentinels built in QA bands'+
  '\n @param {ee.Image} img: Sentinel 2 image'+
  '\n @return {ee.Image} original image with pixels identfied as cloudy or cirrus masked'+
'\n sentinelCloudScore(img): computes a cloud probability score [0-100] and adds as a band to input image'+
  '\n @param {ee.Image} img: Sentinel 2 image'+
  '\n @return {ee.Image} original image with "cloudScore" band added'+
'\n waterScore(img): computes a water probability score [0-1] and adds as a band to input image'+
  '\n @param {ee.Image} img: Sentinel 2 image'+
  '\n @return {ee.Image} image with 1 "waterScore" band';

//EXAMPLE
/*
print(Stats.doc);
var scored = S2.filterBounds(geometry).map(function(image){
  var shadow = darkC(image, "B4", "B3", "B2");
  var cloud = sentinelCloudScore(image);
  return image.addBands(cloud.select('cloudScore')).addBands(shadow.select('C1', 'C2', 'C3'));
});


var zs = TSp(scored, ['C1', 'C2', 'C3', 'B8', 'B11', 'B12']);

var ldas = zs.map(function(image){
  var lda = Stats.ldaScore(image, -0.6772855, ['B11_z', 'B12_z', 'B8_z', 'C1_z', 'C2_z', 'C3_z'], [-1.73655, 1.482, 0.1155737, -0.2842, -0.1455765, 0.53579]);
  return image.addBands(lda.rename('shadowLDA'));
});

var test = ee.Image(ldas.filterDate('2017-06-01', '2018-01-01').filterMetadata("CLOUDY_PIXEL_PERCENTAGE", "less_than", 25).sort("CLOUDY_PIXEL_PERCENTAGE", true).first());
var cloud = test.select('cloudScore');
var shadow = test.select(['B8_z', 'B11_z', 'B12_z', 'C1_z', 'C2_z', 'C3_z', 'shadowLDA']);

//Map.addLayer(test, {bands: ['B4', 'B3', 'B2'], min: 1000, max: 3000}, 'test');
//Map.addLayer(cloud.clip(geometry), {bands: ['cloudScore'], min: 15, max: 50}, 'cloud');
//Map.addLayer(shadow.clip(geometry), {bands: ['C1_z', 'C2_z', 'C3_z']}, 'shadow');

*/
