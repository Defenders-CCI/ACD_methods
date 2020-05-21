# Import the Earth Engine Python Package
import ee
import stats

# Initialize Earth Engine
ee.Initialize()

# Import the Sentinel-2 image collection
S2 = ee.ImageCollection("COPERNICUS/S2")
CDL = ee.ImageCollection("USDA/NASS/CDL/2017")
DEM = ee.ImageCollection("USGS/SRTMGL1_003")

def ND(img, NIR, R, G, SWIR1, SWIR2):
    """
    Calcuate multiple normalized difference metrics
    
    Parameters:
        img (ee.Image): multispectral image
        NIR (string): name of the near infrared band
        R (string): name of the red band
        G (string): ...
        SWIR1 (string):...
        SWIR2 (string):...
        
    Returns:
        ee.Image: image with bands ['ndvi', 'nbr', 'ndwi', 'ndsi']
    """
    NBR = img.normalizedDifference([NIR, SWIR2]).rename(["nbr"])
    NDSI = img.normalizedDifference([G, SWIR1]).divide(img.select([NIR])).rename(["ndsi"])
    NDWI = img.normalizedDifference([G, NIR]).rename(["ndwi"])
    NDVI = img.normalizedDifference([NIR, R]).rename(["ndvi"])
    return img.addBands(NDVI).addBands(NDSI).addBands(NBR).addBands(NDWI)

def CV(b, a, bnds, aoi):
    """
    Calculate the change vector between two images
    
    Parameters:
        b (ee.Image): 'before' image
        a (ee.Image): 'after' image
        bnds (ee.List<str>): band names included in calculation
        aoi (ee.Geometry): area of interest
        
    Returns:
        ee.Image: image with single change vector ['cv'] band
    """
    
    diff = b.select(bnds).subtract(a.select(bnds))
    diff_sd = diff.reduceRegion(
        reducer = ee.Reducer.stdDev(),
        geometry = aoi,
        scale = 30,
        maxPixels = 1e13,
        tileScale = 6).toImage(bnds)
    diff_mn = diff.reduceRegion(
        reducer= ee.Reducer.mean(),
        geometry= aoi,
        scale= 30,
        maxPixels= 1e13,
        tileScale = 6).toImage(bnds)
    return (diff.subtract(diff_mn).divide(diff_sd).pow(2).reduce(ee.Reducer.sum()).rename(['cv']))


def rcvmax(b, a, bnds, aoi):
    """
    Calculate the relative change vector max metric between two images
    
    Parameters:
        b (ee.Image): 'before' image
        a (ee.Image): 'after' image
        bnds (ee.List<str>): band names included in calculation
        aoi (ee.Geometry): area of interest
        
    Returns:
        ee.Image: image with single change vector ['rcvmax'] band
    """
    diff = b.select(bnds).subtract(a.select(bnds))
    # bands = diff.bandNames()
    # print('diff bands:',  bands.getInfo())
    maxab = b.select(bnds).max(a.select(bnds)).pow(2)
    # bands = maxab.bandNames()
    # print('maxab bands:',  bands.getInfo())
    stat = diff.divide(maxab)
    # bands = stat.bandNames()
    # print('stat bands:',  bands.getInfo())
    diff_sd = diff.reduceRegion(
        reducer = ee.Reducer.stdDev(),
        geometry= aoi,
        scale= 30,
        maxPixels= 1e13,
        tileScale=6).toImage(bnds).divide(maxab)
    # bands = diff_sd.bandNames()
    # print('diff_sd bands:',  bands.getInfo())
    
    diff_mn = diff.reduceRegion(
        reducer= ee.Reducer.mean(),
        geometry= aoi,
        scale= 30,
        maxPixels= 1e13,
        tileScale=6).toImage(bnds).divide(maxab)
    # bands = diff_mn.bandNames()
    # print('diff_mn bands:',  bands.getInfo())

    return (stat.subtract(diff_mn).divide(diff_sd).reduce(ee.Reducer.sum()).rename(['rcvmax']))

def d(b, a, bnds):
    """
    Calculate difference between two images
    
    Parameters:
        b (ee.Image): 'before' image
        a (ee.Image): 'after' image
        bnds (ee.List<str>): band names included in calculation
        
    Returns:
        ee.Image: difference image with original bands
    """    
    return b.select(bnds).subtract(a.select(bnds))

# Function(s) to rename bands.
    
def rnm_u(instring):
    newstring = ee.String(instring).cat('_mode')
    return(newstring)

def rnm_s(instring):
    newstring = ee.String(instring).cat('_stdDev')
    return(newstring)

def rnm_z(instring):
    newstring = ee.String(instring).cat('_z')
    return(newstring)

def rnm_p(instring):
    newstring = ee.String(instring).cat('_p')
    return(newstring)

def calc_zp(change, aoi, scl):
    """
    Calculate the z-score and p-values from normal and chi-squared distributions
    
    Parameters:
        change (ee.Image): n-band image containing change metrics
        aoi (ee.Geometry): area of interest
        scl (int): scale at which to sample images for statistic computations
    
    Returns:
        ee.Image: image with 2n bands containing z-score and p-values
    """
    
    norm = change.select(['ndvi', 'ndsi', 'nbr', 'ndwi', 'rcvmax'])
    chi = change.select(['cv']).rename(['cv_z'])
    norm_bands = norm.bandNames()

    cat_mean = norm_bands.map(rnm_u)
    cat_sd = norm_bands.map(rnm_s)
    cat_p = norm_bands.map(rnm_p)
    cat_z = norm_bands.map(rnm_z)

    # mean = change.select(norm_bands).reduceRegion(
    #     reducer=ee.Reducer.mode(),
    #     geometry=aoi,
    #     scale=scl,
    #     maxPixels=1e13).rename(norm_bands, cat_mean)
    #
    # sd = change.select(norm_bands).reduceRegion(
    #     reducer=ee.Reducer.stdDev(),
    #     geometry=aoi,
    #     scale=scl,
    #     maxPixels=1e13).rename(norm_bands, cat_sd)
    #mystats = mean.combine(sd)
    mystats = change.select(norm_bands).reduceRegion(
        reducer=ee.Reducer.mode().combine(
            reducer2=ee.Reducer.stdDev(),
            sharedInputs=True
            ),
        geometry=aoi,
        scale=scl,
        maxPixels=1e13,
        tileScale=6
    ) # //.rename(norm_bands, cat_mean);

    img_z = norm.subtract(mystats.toImage(cat_mean)).divide(mystats.toImage(cat_sd)).rename(cat_z)
    np = stats.norm_p(img_z.abs()).multiply(2).rename(cat_p)
    #np = ee.Image.constant(1).subtract(img_z.abs().multiply(-1.65451).exp().add(1).pow(-1)).multiply(2).rename(ee.String().cat('_p'))

    cp = stats.chi_p(chi, 6).multiply(-1).add(1).rename(['cv_p'])
    return chi.addBands(img_z).addBands(np).addBands(cp)

def iw(change, aoi, niter):
    """
    Iteratively reweight the pixels of an image 
    
    Parameters:
        change (ee.Image): n-band image of change metrics
        aoi (ee.Geometry): area of interest
        niter (int): number of reweighting iterations
    
    Returns:
        ee.Image: z-score image output of calc_zp()
    """
    bands = change.bandNames()
    cat_p = bands.map(rnm_p)
    net = 1
    zs = calc_zp(change, aoi, 300)
    while net <= niter:
        dp = zs.select(cat_p).max(0.001).multiply(change).rename(bands)
        zs = calc_zp(dp, aoi, 30)
        net += 1
    return zs

def runIW(before, after, aoi, ag):
    """
    Run the complete iteratively weighted change analysis
    
    Parameters:
        before (ee.ImageCollection): images representing the reference landscape
        after (ee.ImageCollection): images representing the after condition
        aoi: (ee.Geometry): area of interest
        ag ('yes/no'): mask agricultural areas using Cultivated Lands Dataset?
        
    Returns:
        ee.Image: z-score image output of iw()
        
    """
    CDL = ee.Image("USDA/NASS/CDL/2017")
    DEM = ee.Image("USGS/SRTMGL1_003")

    demMask = DEM.select(['elevation']).lte(3500)
    agMask = CDL.select(['cultivated']).eq(1)
    rgbn = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']

    time1 = after
    time2 = before
    recent = time1.median().clip(aoi)
    past = time2.median().clip(aoi)
    recent = ee.Image(ee.Algorithms.If(ag == 'yes', recent.updateMask(agMask.And(demMask)), recent.updateMask(demMask)))
    past = ee.Image(ee.Algorithms.If(ag == 'yes', past.updateMask(agMask.And(demMask)), past.updateMask(demMask)))
    now = ND(recent, 'B8', 'B4', 'B3', 'B11', 'B12')
    old = ND(past, 'B8', 'B4', 'B3', 'B11', 'B12')

    # bands = now.bandNames()
    # list = bands.getInfo()
    # print('now bands:', bands.getInfo())


    # CREATE IMAGE WITH BANDS FOR CHANGE METRICS CV, RCV, NDVI, NBR, NDSI
    # Calculate cv from before and after images
    cv = CV(old, now, rgbn, aoi)

    # Calculate rcv from before and after images
    rcv = rcvmax(old, now, rgbn, aoi)
    #bands = rcv.bandNames()
    #list = bands.getInfo()
    #print('rcv bands:',  rcv.bandNames().getInfo())

    # Calculate combined normalized difference metrics from before and after images
    diff = d(old, now, ['ndvi', 'ndsi', 'ndwi', 'nbr'])

    #bands = diff.bandNames()
    #list = bands.getInfo()
    #print('diff bands:',  bands.getInfo())

    # Combine cv, rcv, and normalized difference images into single image
    change = cv.addBands(diff).addBands(rcv)

    #bands = change.bandNames()
    #list = bands.getInfo()
    #print('change bands:',  bands.getInfo())
    # zchange not used, but still need to call zp
    #zchange = calc_zp(change, aoi, 30)

    iwchange = iw(change, aoi, 10)
    #bands = iwchange.bandNames()
    #list = bands.getInfo()
    #print('iwchange bands:',  bands.getInfo())


    return iwchange
