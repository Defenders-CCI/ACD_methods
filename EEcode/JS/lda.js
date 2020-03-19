//Function converting multiband image into single band image of LDA scores
//param img: multiband image (ee.Image)
//param int: intercept parameter from LDA analysis (numeric)
//param xbands: string list of n band names (list)
//param coefficients: numeric list of length n containing LDA coefficients (list)
exports.ldaScore = function (img, int, xbands, coefficients){
  var bands = img.select(xbands);
  
  var coeffs = ee.Dictionary.fromLists(xbands, coefficients).toImage(xbands);
  
  var score = bands.multiply(coeffs).addBands(ee.Image(int)).reduce(ee.Reducer.sum());

  return score;
  };