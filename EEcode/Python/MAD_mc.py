# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:56:04 2019

@author: MEvans
"""

import ee

ee.Initialize()
#util = require('users/mortcanty/changedetection:utilities')
# ****************************
# Modules for IR-MAD algorithm
# ****************************

def chi2cdf(chi2,df):
  #Chi square cumulative distribution function '''
  return ee.Image(chi2).divide(2).gammainc(ee.Number(df).divide(2))

def covarw(image, weights, nPix):
# Return the weighted centered image and its weighted covariance matrix''' 
#var image = image1.addBands(image2)
#var weights = ee.Image.constant(1)
    image = ee.Image(image) 
    geometry = image.geometry() 
    bandNames = image.bandNames()
    N = bandNames.length() 
    scale = 10
    #var scale = image.select(0).projection().nominalScale()
    weightsImage = image.multiply(ee.Image.constant(0)).add(weights)
    means = image.addBands(weightsImage) \
         .reduceRegion(
           reducer = ee.Reducer.mean().repeat(N).splitWeights(),
           geometry = geometry,
           scale = 100,
           maxPixels = 1e13,
           tileScale = 8) \
         .toArray() \
         .project([1])
    centered = image.toArray().subtract(means)
    B1 = centered.bandNames().get(0)
    b1 = weights.bandNames().get(0)
    nPixels = nPix
    #nPixels = ee.Number(centered.reduceRegion(ee.Reducer.count(),geometry,scale,null,null,false,1e9).get(B1)) 
    #nPixels = ee.Number(weights.reduceRegion(ee.Reducer.count(),geometry,scale,None,None,False,1e9).get(b1)) 
    sumWeights = ee.Number(weights.reduceRegion(ee.Reducer.sum(),geometry,scale,None,None,False,1e9).get(b1))
    covw1 = centered.multiply(weights.sqrt()) \
                   .toArray() \
                   .reduceRegion(
                     reducer = ee.Reducer.centeredCovariance(),
                     geometry = geometry,
                     scale = 100,
                     maxPixels = 1e13,
                     tileScale = 8) \
                   .get('array')
    covw = ee.Array(covw1).multiply(nPixels).divide(sumWeights)
    return ee.Dictionary({'centeredimage':centered.arrayFlatten([bandNames]), 'covw': covw})

def geneiv(C1,B1):
# Generalized eigenproblem C*X = lambda*B*X
    C = ee.Array(C1)
    B = ee.Array(B1)  
#  Li = choldc(B)^-1
    Li = ee.Array(B.matrixCholeskyDecomposition().get('L')).matrixInverse()
#  solve symmetric eigenproblem Li*C*Li^T*x = lambda*x
    Xa = Li.matrixMultiply(C) \
               .matrixMultiply(Li.matrixTranspose()) \
               .eigen()
#  eigenvalues as a row vector
    lambdas = Xa.slice(1,0,1).matrixTranspose()
#  eigenvectors as columns
    X = Xa.slice(1,1).matrixTranspose()  
#  generalized eigenvectors as columns, Li^T*X
    eigenvecs = Li.matrixTranspose().matrixMultiply(X)
    return ee.Dictionary({'lambdas':lambdas, 'eigenvecs':eigenvecs})

def radcal (current,prev):
# iterator function for orthogonal regression and interactive radiometric normalization 
    k = ee.Number(current)
    prev = ee.Dictionary(prev)    
#  image is concatenation of reference and target    
    image = ee.Image(prev.get('image'))
    ncmask = ee.Image(prev.get('ncmask'))
    nbands = ee.Number(prev.get('nbands'))
    geometry = ee.Geometry(prev.get('geometry'))
    coeffs = ee.List(prev.get('coeffs'))
    normalized = ee.Image(prev.get('normalized'))
    scale = image.select(0).projection().nominalScale()
#  orthoregress reference onto target  
    image1 = image.clip(geometry).select(k.add(nbands),k).updateMask(ncmask).rename(['x','y'])
    means = image1.reduceRegion(ee.Reducer.mean(),geometry,scale,None,None,False,1e9) \
                  .toArray() \
                  .project([0])              
    Xm = means.get([0])
    Ym = means.get([1])
    S = ee.Array(image1.toArray() 
                           .reduceRegion(ee.Reducer.covariance(),geometry,scale,None,None,False,1e9) 
                           .get('array'))     
#  Pearson correlation     
    R = S.get([0,1]).divide(S.get([0,0]).multiply(S.get([1,1])).sqrt())
    eivs = S.eigen()
    e1 = eivs.get([0,1])
    e2 = eivs.get([0,2])
#  slope and intercept    
    b = e2.divide(e1)
    a = Ym.subtract(b.multiply(Xm))
    coeffs = coeffs.add(ee.List([b,a,R]))
#  normalize kth band in target    
    normalized = normalized.addBands(image.select(k.add(nbands)).multiply(b).add(a))
    return ee.Dictionary({'image':image,'ncmask':ncmask,'nbands':nbands,'geometry':geometry,'coeffs':coeffs,'normalized':normalized})    


def imad1(current,prev):
# Iteratively re-weighted MAD
    image = ee.Image(ee.Dictionary(prev).get('image'))
    chi2 = ee.Image(ee.Dictionary(prev).get('chi2'))
    allrhos = ee.List(ee.Dictionary(prev).get('allrhos'))
    nPix = ee.Dictionary(prev).get('size')
    region = image.geometry()
    nBands = image.bandNames().length().divide(2)
    weights = chi2cdf(chi2,nBands).subtract(1).multiply(-1)
    result = covarw(image,weights,nPix) 
    centeredImage = ee.Image(result.get('centeredimage'))
    #covarArray is 2n covar matrix of bands in image 1 and image 2
    covarArray = ee.Array(result.get('covw'))
    bNames = centeredImage.bandNames()
    bNames1 = bNames.slice(0,nBands)
    bNames2 = bNames.slice(nBands)
    centeredImage1 = centeredImage.select(bNames1)
    centeredImage2 = centeredImage.select(bNames2)
    s11 = covarArray.slice(0,0,nBands).slice(1,0,nBands)
    s22 = covarArray.slice(0,nBands).slice(1,nBands)
    s12 = covarArray.slice(0,0,nBands).slice(1,nBands)
    s21 = covarArray.slice(0,nBands).slice(1,0,nBands)      
    c1 = s12.matrixMultiply(s22.matrixInverse()).matrixMultiply(s21)
    b1 = s11
    c2 = s21.matrixMultiply(s11.matrixInverse()).matrixMultiply(s12)
    b2 = s22
#  solution of generalized eigenproblems 
    result = geneiv(c1,b1)
    lambdas = ee.Array(result.get('lambdas'))
    A = ee.Array(result.get('eigenvecs'))
    result = geneiv(c2,b2)
    B = ee.Array(result.get('eigenvecs'))
    rhos = lambdas.sqrt().project(ee.List([1]))
#  sort in increasing order
    keys = ee.List.sequence(nBands,1,-1)
    A = A.sort([keys])
    B = B.sort([keys]) 
    rhos = rhos.sort(keys)
#  test for convergence    
    lastrhos = ee.Array(allrhos.get(-1))
    done = rhos.subtract(lastrhos) \
                   .abs() \
                   .reduce(ee.Reducer.max(),ee.List([0])) \
                   .lt(ee.Number(0.001)) \
                   .toList() \
                   .get(0)      
    allrhos = allrhos.cat([rhos.toList()]) 
#  MAD variances    
    sigma2s = rhos.subtract(1).multiply(-2).toList()
    sigma2s = ee.Image.constant(sigma2s)
#  ensure sum of positive correlations between X and U is positive
    tmp = s11.matrixDiagonal().sqrt()
    ones = tmp.multiply(0).add(1)
    tmp = ones.divide(tmp).matrixToDiag()
    s = tmp.matrixMultiply(s11).matrixMultiply(A).reduce(ee.Reducer.sum(),[0]).transpose()
    A = A.matrixMultiply(s.divide(s.abs()).matrixToDiag())
#  ensure positive correlation
    tmp = A.transpose().matrixMultiply(s12).matrixMultiply(B).matrixDiagonal()
    tmp = tmp.divide(tmp.abs()).matrixToDiag()
    B = B.matrixMultiply(tmp)        
#  canonical and MAD variates 
    centeredImage1Array = centeredImage1.toArray().toArray(1)
    centeredImage2Array = centeredImage2.toArray().toArray(1)
    U = ee.Image(A.transpose()).matrixMultiply(centeredImage1Array) \
                      .arrayProject([0]) \
                      .arrayFlatten([bNames1]) 
    V = ee.Image(B.transpose()).matrixMultiply(centeredImage2Array) \
                      .arrayProject([0]) \
                      .arrayFlatten([bNames2])   
    MAD = U.subtract(V)
#  chi square image
    chi2 = MAD.pow(2) \
              .divide(sigma2s) \
              .reduce(ee.Reducer.sum()) \
              .clip(region)               
    return ee.Dictionary({'done':done,'image':image,'allrhos':allrhos,'chi2':chi2,'MAD':MAD, 'size': nPix})    


# the algorithm
def imad(current,prev):
    done =  ee.Number(ee.Dictionary(prev).get('done')) 
    return ee.Algorithms.If(done,prev,imad1(current,prev))


