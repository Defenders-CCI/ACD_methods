import ee

# Initialize Earth Engine
#ee.Initialize()


# Function calculating hillshade for an image using the image's solar angle properties
def shade(img, elev):
    az = img.get('MEAN_SOLAR_AZIMUTH_ANGLE')
    solel = ee.Number(90).subtract(img.get('MEAN_SOLAR_ZENITH_ANGLE'))
    shd = ee.Terrain.hillshade(elev, az, solel)
    return shd


def illuminate(img, elev):
    az = ee.Image.constant(img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
    solel = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'))
    aspect = ee.Terrain.aspect(elev)
    slope = ee.Terrain.slope(elev)
    illum = slope.cos().multiply(solel.cos()).add(slope.sin().multiply(solel.sin()).multiply(az.subtract(aspect).cos()))
    return illum.rename(['illumination'])

# def illumImg(image):
#     return image.addBands(ee.Image(1)).addBands(illuminate(image, elev.clip(aoi)))

# Function calculating and correcting hillshade for each image in a collection.
# Return Image Collection
def c_correct(imgCol, bands, aoi, elev):
    print('Running c_correct algorithm')
    otherbands = ee.Image(imgCol.first()).bandNames().removeAll(bands)
    #print('otherbands:', otherbands)
    illumCol = imgCol.map(lambda image: image.addBands(ee.Image(1)).addBands(illuminate(image, elev.clip(aoi))))
    #illumCol.map(illumImg)

    regbands = ee.List(['constant', 'illumination']).cat(bands)
    coeffs = illumCol.select(regbands).reduce(ee.Reducer.linearRegression(
        numX = 2,
        numY = ee.List(bands).length())
    ).select('coefficients')
    # Coefficients is 2d array (numX by numY). Return only the slope coefficients
    # Slope coefficients are stil 2d (1xnumY).  Transform to 1d array of length numY
    # Transform the array image into multiband image
    betas = coeffs.arraySlice(0, 1, 2).arrayProject([1]).arrayFlatten([bands])

    intercept = coeffs.arraySlice(0, 0, 1).arrayProject([1]).arrayFlatten([bands])

    c = intercept.divide(betas)

    # not used
    #mnshade = illumCol.select('illumination').reduce(ee.Reducer.mean())

    def correct(image):
        zen = ee.Image.constant(image.get('MEAN_SOLAR_ZENITH_ANGLE')).cos()
        num = zen.add(c).divide(image.select('illumination').add(c))
        correction = image.select(bands).multiply(num)
        return image.select(otherbands.cat(['illumination'])).addBands(correction)

    corrected = illumCol.map(correct)

    return corrected
