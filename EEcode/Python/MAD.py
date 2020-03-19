import ee
from stats import chi_p

ee.Initialize()
#function paramaterizing a gamma distribution from data, and returning p-values of observations
#param chi: image of chi-square statistics, or other values following a gamma distribution
#param aoi: geometry object defining area of analisys
#param scl: scale (m) at which reduction operations are performed
#return: image with single band ['p'] containing p-values of original image values
def gamma(chi, aoi, scl):
    """
    Parameterize a gamma distribution from data and return p-values
    
    Parameters:
        chi (ee.Image): image containing chi-square statistics or other values
        following a gamma distribution
        aoi (ee.Geometry): area of interest
        scl (int): scale (m) at which reduction operations are performed
    Returns:
        ee.Image: single band ['p']
    """
    mode = chi.reduceRegion(
      reducer = ee.Reducer.mode(),
      geometry = aoi,
      scale = scl,
      maxPixels = 1e13
    ).toImage()
    sd = chi.reduceRegion(
      reducer = ee.Reducer.sampleStdDev(),
      geometry = aoi,
      scale = scl,
      maxPixels = 1e13
    ).toImage()
    rate = (((((sd.pow(2)).multiply(4)).add(mode.pow(2))).sqrt()).add(mode)).divide(sd.pow(2).multiply(2))
    #shape calculated from observed mode
    shape = mode.multiply(rate).add(1)
    #shape calculated assuming mode = 0 (i.e. no change)
    #shape = ee.Image.constant(1)
    p = shape.gammainc(chi.multiply(rate)).rename(['p'])
    return p

def rnm_u(band):
    return ee.String(band).cat('_mean')

def rnm_s(band):
    return ee.String(band).cat('_stdDev')

def chisq(img, aoi, scl):
    """
    Compute the chi square statistic and p-value from an image as the sum
    of squared differences in pixels from band means
    
    Parameters:
        img (ee.Image): multiband image
        aoi (ee.Geometry): area of interest
        scl (int): scale at which reduction operations are performed
        
    Returns:
        ee.Image: bands ['chi', 'p']
        
    """
    k = img.bandNames().length()
    k2 = k.divide(2)
  #return dictionary of standard deviations for each input band
    stats = img.reduceRegion(
            reducer = ee.Reducer.mean().combine(
                    reducer2 = ee.Reducer.stdDev(), sharedInputs = true),
            geometry = aoi,
            scale = scl,
            maxPixels = 1e13,
            tileScale = 4
            )

  #calculte chi-square statistic as sum(i:p) of (MADij/stdDev(i))^2
    means = img.bandNames().map(rnm_u)
  
    sds = img.bandNames().map(rnm_s)
  
    mean = stats.toImage(means)
    sd = stats.toImage(sds)
    chi = img.subtract(mean).divide(sd).pow(2).reduce(ee.Reducer.sum()).rename(['chi'])
  #chi = (((img.subtract(mean)).divide(sd)).pow(2)).reduce(ee.Reducer.sum()).rename(['chi'])
  #p = ee.Image.constant(1).subtract(chi.divide(2).gammainc(k2)).rename(['p'])
    p = chi_p(chi, k).multiply(-1).add(1)
    return chi.addBands(p)

def correlate(x, y, aoi):
    selected = x.addBands(y)
    cor = selected.reduceRegion(
            reducer = ee.Reducer.pearsonsCorrelation(),
            geometry = aoi,
            scale = 100,
            maxPixels = 1e13,
            tileScale = 4
            ).get('correlation')
    return cor

def corrmat(set1, set2, aoi):
    """
    Compute a correlation matrix between two sets of variables
    
    Parameters:
        set1 (ee.Image): image with n bands equal to first set of variables
        set2 (ee.Image): image with m bands equal to second set of variables
        aoi (ee.Geometry): area of interest
        
    Returns:
        ee.Array: n x m array containing correlation coefficients between
        variables in set1 and set2
    """
    emptylist = ee.List([])
    xbands = set1.bandNames()
    ybands = set2.bandNames()
    ylen = len(ybands)
    cors = [correlate(set1.select(x), set2.select(y), aoi) for x in xbands for y in ybands]
    nestedcors = [cors[i:i+ylen] for i in range(0, len(cors), ylen)]
    return ee.Array(nestedcors)
# =============================================================================
#   old javascript way of creating a correlation matrix
#   xs = xbands.iterate(function(xband, list){
#     ys = ybands.iterate(function(yband, list){
#       selected = set1.select([xband]).addBands(set2.select([yband]))
#       
#       cor = selected.reduceRegion(
#         reducer = ee.Reducer.pearsonsCorrelation(),
#         geometry = aoi,
#         scale = 100,
#         maxPixels = 1e13,
#         tileScale = 4
#         ).get('correlation')
#         
#       return ee.List(list).add(cor)   
#     }, emptylist)
# 
#     return ee.List(list).add(ys)
#   }, emptylist)
# =============================================================================



#Function to compute correlation matrix between two sets of variables 
#param {ee.Image} set1: N-band image equal to first set of variables
#param {ee.Image} set2: M-band image equal to second set of variables
#param {ee.Geometry} aoi: Area of analisys
#param {ee.Image} weights Single band image containing per pixel weights
#Return: n x m array containing correlation coefficients between variables in set1 and set2
def wCorrmat(set1, set2, aoi, weights):
    """
    Compute a correlation matrix between two sets of variables
    
    Parameters:
        set1 (ee.Image): image with n bands equal to first set of variables
        set2 (ee.Image): image with m bands equal to second set of variables
        aoi (ee.Geometry): area of interest
        weights (ee.Image): single band iage containing per pixel weights
        
    Returns:
        ee.Array: n x m array containing correlation coefficients between
    """
    
    weights = weights.clip(aoi)
    
    b1 = weights.bandNames().get(0)
    
    n = set1.bandNames().length()
    
    img = set1.addBands(set2)
    
    weightsImage = img.multiply(ee.Image.constant(0).add(weights)) 
  
    N = img.bandNames().length()
  
    xbands = ee.List(img.bandNames().slice(0, n))

    ybands = ee.List(img.bandNames().slice(n, N))

    sumWeights = ee.Number(weights.reduceRegion(
        reducer = ee.Reducer.sum(),
        geometry = aoi,
        scale = 10,
        tileScale = 4,
        maxPixels = 1e13).get(b1))
    
    nPixels = ee.Number(weights.reduceRegion(
        reducer = ee.Reducer.count(),
        geometry = aoi,
        scale = 10,
        tileScale = 4,
        maxPixels = 1e13).get(b1))

    means = img.addBands(weightsImage).reduceRegion(
            reducer = ee.Reducer.mean().repeat(N).splitWeights(),
            geometry = aoi,
            scale = 100,
            maxPixels = 1e13
            ).toArray().project([1])

    centered = img.toArray().subtract(means)
  
    weighted = img.multiply(weights.sqrt())
  
    cors = [correlate(set1.select(x), set2.select(y), aoi) for x in xbands for y in ybands]
    nestedcors = [cors[i:i+ylen] for i in range(0, len(cors), ylen)]
    corw = ee.Array(nestedcors).multiply(nPixels).divide(sumWeights)
    return corw
# =============================================================================
#   xs = xbands.map(function(xbnd){
#     ys = ybands.map(function(ybnd){
# 
#       selected = weighted.select([xbnd, ybnd])#.addBands(weighted.select([ybnd]))
#       
#       cor = selected.reduceRegion(
#         reducer = ee.Reducer.pearsonsCorrelation(),
#         geometry = aoi,
#         scale = 100,
#         maxPixels = 1e13,
#         tileScale = 4
#         ).get('correlation')
#       
#       return cor
#         })
#     return ys
#   })
# =============================================================================


def canon(cor, before, after, length):
    """
    Perform canonical correlation analysis using SVD on a correlation matrix
    Parameters:
        cor (ee.Array): n x m correlation matrix
        before (ee.Image): image with n bands
        after (ee.Image): image with m bands
        length (int): minimum of [n,m]
    Returns:
        image with j = min(n,m) bands ['V1', 'V2',...'Vj', 'chi', 'p'] corresponding to differences aX - bY across j canonical variates
  """
    labels = ee.List.sequence(1, length).map(lambda item: ee.String("V").cat(ee.Number(item).toInt().format()))
    decomp = cor.matrixSingularValueDecomposition()
      #create n * min(n, m) matrix of canonical covariates a
    U = ee.Array(decomp.get("U"))
      #create m * min(n, m) matrix of canonical covariates b
    V = ee.Array(decomp.get("V"))
      #get diagonal elements of SVD equivalent to CCA correlations
    S = ee.Array(decomp.get("S")).matrixDiagonal()
      #turn images into array images with 1 x nbands matrices at each pixel
    before2D = before.toArray().toArray(1).arrayTranspose() 
    after2D = after.toArray().toArray(1).arrayTranspose()
    a = before2D.matrixMultiply(ee.Image(U)).arrayProject([1]).arrayFlatten([labels])
    b = after2D.matrixMultiply(ee.Image(V)).arrayProject([1]).arrayFlatten([labels])
    return ee.Dictionary({'img': a.subtract(b), 'rhos': S})

def mad(before, after, aoi):
    """
    Conduct the MAD transformation between two images
    
    Parameters:
        before (ee.Image): n band image 
        after (ee.Image): m band image 
        aoi (ee.Geometry): area of interest
    Returns:
        ee.Image: with j = min(n, m) bands ['V1', 'V2',...'Vj', 'chi', 'p]
    """
    length = before.bandNames().length().min(after.bandNames().length())
    corMat = corrmat(before, after, aoi)
    cca = canon(corMat, before, after, length)
    cvs = ee.Image(cca.get('img'))
    ccachi = chisq(cvs, aoi, 100)
    rhos = cca.get('rhos')
    return cvs.addBands(ccachi)
  #return ee.Dictionary({'img': cvs.addBands(ccachi), 'rhos':rhos})


#Function to conduct weighted MAD transformation between two images
# using single value decomposition on correlation matrix
#param before (ee.Image): before image
#param after (ee.Image): after image
#param aoi (ee.Geometry): area of analysis
#Return ee.Dictionary 'img' j + 2 bands ["V1", "V2", ..."Vj", "chi", "p"] 
def madw(before, after, aoi, weights):
  length = before.bandNames().length().min(after.bandNames().length())
  corMat = wCorrmat(before, after, aoi, weights)
  cca = canon(corMat, before, after, length)
  cvs = ee.Image(cca.get('img'))
  ccachi = chisq(cvs, aoi, 100)
  rhos = ee.Array(cca.get('rhos'))
  #net++
  return ee.Dictionary({
    'img': cvs.addBands(ccachi),
    'rhos':rhos,
    'before': before,
    'after': after,
    'aoi': aoi
  })


# Function running weighted MAD using previous MAD output.  For use with iterate
# param current (): placeholder for use in GEE iterate
# param prev(ee.Dictionary): output from madw containing
#                            before, after, and mad output images
# returns (ee.Dictionary): 
def iw_mad(current, prev):
  prev = ee.Dictionary(prev)
  aoi = ee.Geometry(prev.get('aoi'))
  rhos = ee.Array(prev.get('rhos'))
  before = ee.Image(prev.get('before'))
  after = ee.Image(prev.get('after'))
  length = before.bandNames().length().min(after.bandNames().length())
  output = ee.Image(prev.get('img'))
  weights = output.select(['p'])
  corMat = wCorrmat(before, after, aoi, weights)
  cca = canon(corMat, before, after, length)
  cvs = ee.Image(cca.get('img'))
  ccachi = chisq(cvs, aoi, 100)
  updates = ee.Array(cca.get('rhos'))
  net = updates.sort()
    .subtract(rhos).abs()
    .reduce(ee.Reducer.max(), [0])
    .get([0,0])
  converged = net.lte(0.01)
  result = ee.Dictionary({
    'img': cvs.addBands(ccachi),
    'rhos': updates,
    'before': before,
    'after': after,
    'aoi': aoi
  })
  return ee.Algorithms.If(converged, prev, result)


#Function to perform iterative reweighting and MAD transformations
#param niter: number of iterations
#Returns output same as mad()
def mad_iw(before, after, aoi):
  cca = mad(before, after, aoi, ee.Image(1))
  net = rho.subtract(cca.get('rho'))

  while net > 0.001:
    ps = cca.select(['p']).sqrt()

    sumWeights = ee.Number(ps.reduceRegion(
      reducer = ee.Reducer.sum(),
      geoetry = aoi,
      scale = 10,
      maxPixels = 1e13).get('p'))
   

# =============================================================================
#     nPixels = ps.reduceRegion(
#       reducer = ee.Reducer.count(),
#       geometry = aoi,
#       scale = 10,
#       maxPixels = 1e13).get('p')
# =============================================================================

    weights = ee.Image.constant(count).divide(ee.Image.constant(sumweight)).multiply(ps)
    
    #This version will apply new weights to previously re-weighted images
    #before = before.multiply(weights)
    #after = after.multiply(weights)
    #cca = mad(before, after, geometry, 30)
    #This version will apply new weights to original images each iteration
    b = before#.multiply(ps)
    a = after#.multiply(ps)
    cca = mad(b, a, aoi, ps)
    net++
  
  return cca

def runMAD(before, after, aoi, ag):
    """
    Run the MAD analysis. Currently uses a single mad iteration
    
    Parameters:
        before (ee.ImageCollection): image collection representing baseline state
        after (ee.ImageCollection): image collection representing current state
        aoi (ee.Geometry): are of interest
        ag ('Yes/No'): mask agricultural areas?
    Returns:
        ee.Image
    """
    DEM = ee.Image("USGS/SRTMGL1_003")
    CDL = ee.Image("USDA/NASS/CDL/2017")
    demMask = DEM.select('elevation').lte(3500)
    agMask = CDL.select('cultivated').eq(1)
    before = before.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'])
    after = after.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'])
    
    recent = after.median().clip(aoi)
    past = before.median().clip(aoi)
      
    recent = ee.Image(
      ee.Algorithms.If(ag == 'yes', recent.updateMask(agMask.and(demMask)), recent.updateMask(demMask))
      )
    past = ee.Image(
      ee.Algorithms.If(ag == 'yes', past.updateMask(agMask.and(demMask)), past.updateMask(demMask))
      )
      
    mad = mad(past, recent, aoi, 100)
    chi = chisq(mad, aoi, 100)
    #  mad = ee.List.sequence(1, 10).iterate(iw_mad, madw(past, recent, aoi, ee.Image(1)))
    #  madout = ee.Image(ee.Dictionary(mad).get('img'))
    #  return madout
    return mad.addBands(chi)

#EXAMPLE
# =============================================================================
# 
# LS8 = LS8.filterBounds(geometry).filterMetadata('CLOUD_COVER', 'less_than', 10)
# 
# before = LS8.filterDate('2016-06-01', '2016-09-01')
# after = LS8. filterDate('2017-06-01', '2017-09-01')
# 
# 
# before = before.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'])
# after = after.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'])
# 
# recent = after.median().clip(geometry)
# past = before.median().clip(geometry)
#   
# Map.addLayer(past, {bands:['B4', 'B3', 'B2'], min:300, max:3000, gamma:1}, 'before')
# Map.addLayer(recent, {bands:['B4', 'B3', 'B2'], min:300, max:3000, gamma:1}, 'after')
# 
# cca = mad(past, recent, geometry, ee.Image(1))
# Map.addLayer(cca, {}, 'mad')
# =============================================================================









