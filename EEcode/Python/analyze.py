import ee

from datetime import datetime
import clouds
import iw
import MAD_mc
import stats
import terrain
import sys
import os

# Initialize Earth Engine
ee.Initialize()

CDL = ee.Image("USDA/NASS/CDL/2017")
S2 = ee.ImageCollection("COPERNICUS/S2")
DEM = ee.Image("USGS/SRTMGL1_003")

def mask(img):
    scored = clouds.basicQA(img)
    scored = clouds.sentinelCloudScore(scored)
    waterMask = clouds.waterScore(img).select('waterScore').lte(0.5)
    shadowMask = img.select('B11').gt(900)
    return scored.updateMask(scored.select('cloudScore').lte(15).And(shadowMask).And(waterMask))

def sz(ft):
    area = ft.area(5)
    return ft.set({'area': area})

def displaySize(size):
    print('size:', size)

#def analyze_iw(aoi, doi, cvz, nbrz, ndsiz, ndviz, ndwiz, rcvz, size, intercept, lda):
def analyze_iw(aoi, doi, dictionary, size):
    """
    Function that pre-processes sentinel-2 imagery and runs the LCC change detection algorithm
    
    Parameters:
        aoi(ee.Geometry): area of interest
        doi(ee.Date): date of interest
        cvz(float): LDA coefficeint associated with the cv metric
        nbrz (float): LDA coefficient associated with the NBR metric
        ndsiz (float): LDA coefficient associated with the NDSI metric
        ndviz (float): ...
        ndwiz (float): ...
        rcvz (float): LDA coefficient associated with the rcv_max metric
        size (float): minimum size (ac) of changes to output
        intercept (float): LDA intercept coefficient
        lda (float): LDA score threshold defining change/no-change
        
    Returns:
        ee.Geometry: ?
    """
    try:
        sq_meters = ee.Number(size).multiply(4047)
        projdate = ee.Date(doi);
        today = projdate.advance(3, 'month');
    
        proj_dt = str(datetime.fromtimestamp(int(projdate.getInfo()['value']) / 1e3))[:10]
        print('proj_dt:', proj_dt)
        prior = ee.Date.fromYMD(projdate.get('year').subtract(1), projdate.get('month'), projdate.get('day'))
        prior_dt = str(datetime.fromtimestamp(int(prior.getInfo()['value']) / 1e3))[:10]
        print('prior_dt:', prior_dt)

        rgbn = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']

        masked = S2.filterDate(prior, today).filterBounds(aoi).map(mask)
        corrected = terrain.c_correct(masked, rgbn, aoi, DEM)

        after = corrected.filterDate(projdate, today)
        count = after.size()
        #print('after size:', count.getInfo())
        reference = after.sort('system:time_start', False)
        time0 = ee.Image(reference.first()).get('system:time_start')
        recent_date = str(datetime.fromtimestamp(int(time0.getInfo()) / 1e3))[:10]

        before = corrected.filterDate(prior, projdate)
        count = before.size()
        print('before size:', count.getInfo())
        reference = before.sort('system:time_start', False)
        time0 = reference.first().get('system:time_start')
        past_date = str(datetime.fromtimestamp(int(time0.getInfo()) / 1e3))[:10]
 
        # run the IW algorithm between the before and after collections within the user defined AOI.
        # by default, ag fields are masked by 'yes'
        iwout = iw.runIW(before, after, aoi, 'yes').clip(aoi)

        # calculate LDA score to discriminate change/no-change pixels in iwout.  Requires thresholds from habitat dictionary        
        scored = stats.ldaScore(
                iwout,
                ['cv_z', 'rcvmax_z', 'ndvi_z', 'ndsi_z', 'ndwi_z', 'nbr_z'],
                dictionary)
        
#        scored = stats.ldaScore(iwout, 0 ['cv_z', 'rcvmax_z', 'ndvi_z', 'ndsi_z', 'ndwi_z', 'nbr_z'],
#                                [cvz, rcvz, ndviz, ndsiz, ndwiz, nbrz]).clip(aoi)

        # create a binary [0, 1] image representing change and no-change pixels.  Erode and dilate changed areas
        selected = scored.gte(dictionary['lda'])\
        .focal_min(1, 'square', 'pixels')\
        .focal_max(1, 'square', 'pixels')

        # mask image to retain only pixels equal to '1'
        selected = selected.updateMask(selected)
        #maskSelected = selected.updateMask(selected.eq(0))
        # mask out no-change areas (i.e. selected = 0) here.  Creates fewer polygons which should save memory
        # selected = selected.updateMask(selected.eq(1))
        #print('selected is a ', type(selected))

        scale = 10
        tileScale = 6
        
        # convert binary image to polygons.  Note: this creates polygons for both contiguous change and contiguous no-change areas
        polys = selected.reduceToVectors(
            geometry=aoi,
            scale=scale,
            tileScale=tileScale,
            eightConnected=True,
            bestEffort=True,
            maxPixels=1e13)

        #print('polys is a ', type(polys))
        count = polys.size()
        print(count.getInfo())
        #print('polys size:', count.getInfo(displaySize))

        # return only polygons corresponding to change pixels
        polys = polys.map(sz)

        # filter out change polygons smaller than the user defined minimum area
        polys = polys.filter(ee.Filter.gte('area', sq_meters))

        # indicator = True


#        task_keyed = ee.batch.Export.table.toCloudStorage(
#            collection=polys,
#            description=cd_id,
#            fileFormat= 'GeoJSON',
#            bucket='alerts-storage',
#            path='geojson/' + cd_id
#        )
#        task_keyed.start()

        return "OK", past_date, recent_date, polys
    except Exception as error:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print ("")
        print ("*******************************")
        print ("Unexpected error in analyze.py")
        print (exc_type, fname, exc_tb.tb_lineno)
        #print("sys.exc_info:", sys.exc_info()[0])
        print ("Error:", error)
        print ("*******************************")
        print ("")
        return "error"

def analyze_mad(aoi, doi, size, niters):
    iterations = niters
    # service_account = config.EE_ACCOUNT
    # credentials = ee.ServiceAccountCredentials(service_account, config.EE_PRIVATE_KEY_FILE)
    # ee.Initialize(credentials)
    sq_meters = ee.Number(size).multiply(4047)
    projdate = ee.Date(doi);
    today = projdate.advance(3, 'month');
    #projdate = ee.Date(doi)
    #today = ee.Date(datetime.now())
    #today = ee.Date('2018-11-01')
    
    prior = projdate.advance(-1, 'year')
    
    rgbn = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
    
    nbands = len(rgbn)
    print(nbands)
    masked = S2.filterDate(prior, today).filterBounds(aoi).map(mask)
    
    corrected = terrain.c_correct(masked, rgbn, aoi, DEM)

    after = corrected.filterDate(projdate, today)
    image2 = after.median().select(rgbn)
    
    before = corrected.filterDate(prior, projdate)
    image1 = before.median().select(rgbn)
    
    image = image1.addBands(image2).clip(aoi)
    
    npix = image.select(0).reduceRegion(
            reducer = ee.Reducer.count(),
            geometry = aoi,
            scale = 10,
            maxPixels = 1e13,
            tileScale = 4).get('B2')
    
    inputlist = ee.List.sequence(1, iterations)
        
    first = ee.Dictionary({
            'done':ee.Number(0),
            'image': image,
            'allrhos': [ee.List.sequence(1, nbands)],
            'chi2': ee.Image.constant(1),
            'MAD': ee.Image.constant(0),
            'size': npix})
    output = ee.Dictionary(inputlist.iterate(MAD_mc.imad, first))
    return output

