# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.meas.multifit as measMult
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipses
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import sys
import lsst.pex.policy as pexPol

import numpy
import lsst.afw.display.ds9 as ds9
import eups
import math

def main():
    numpy.set_printoptions(threshold=numpy.nan)
    i =0
    centroid = afwGeom.makePointD(45,45)
    psf = afwDet.createPsf("DoubleGaussian", 9, 9, 1.0)
    affine = afwGeom.AffineTransform()

    flux = 35.0
    axes = geomEllipses.Axes(30,25,0).transform(affine)
    print >> sys.stderr, axes

    logShear = geomEllipses.LogShear(axes)

    sersicIndex = 1.0
    try:
        cache = measMult.Cache.load("/home/dubcovsky/multifit/cache_10x4000")
    except:
        pol = pexPol.Policy()
        cache = measMult.makeRobustSersicCache(pol)
        cache.save("robustCache")
    
    measMult.SersicMorphology.setSersicCache(cache)

    model = measMult.createSersicModel(flux, centroid, logShear, sersicIndex)
    #model = measMult.createExponentialModel(flux, centroid, logShear)
    #model = measMult.createPointSourceModel(flux, centroid)

    fp = model.computeProjectionFootprint(psf, affine)
    box = fp.getBBox()
    bigBox = afwImage.BBox(box.getLLC(), box.getWidth()+120, box.getHeight()+120);
    bigBox.shift(-60, -60);
    bigFP = afwDet.Footprint(bigBox);

    proj = model.makeProjection(psf, affine, bigFP)

    modelImage = afwImage.MaskedImageF(bigBox.getWidth(), bigBox.getHeight())
    modelImage.setXY0(bigBox.getX0(), bigBox.getY0())
    
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros_like(imageVector)
    measMult.expandImageF(bigFP, modelImage, imageVector, varianceVector)

    analyticDerivative = proj.computeNonlinearParameterDerivative()
    nonlinearParameters = model.getNonlinearParameters()
    
    eps = 2.2e-16

    
    columnVectors = []
    for n in range(model.getNonlinearParameterSize()):
        h = math.sqrt(eps) * nonlinearParameters[n]
        nonlinearParameters[n] += h        
        print  >> sys.stderr, "n: %g"%n
        print  >> sys.stderr, "h: %g"%h
        print >> sys.stderr, "params plus h: %s"%nonlinearParameters
        model.setNonlinearParameters(nonlinearParameters)
        print >> sys.stderr, model.getNonlinearParameters()

        plus = numpy.copy(proj.computeModelImage())

        nonlinearParameters[n] -= 2*h

        print >> sys.stderr, "params minus h: %s"%nonlinearParameters
        model.setNonlinearParameters(nonlinearParameters)
        minus = numpy.copy(proj.computeModelImage())

        partial = (plus - minus)
        partial = partial / (2*h)
        
        var = numpy.zeros_like(partial.base[0])
        fpBox = fp.getBBox()
        derivativeImage = afwImage.MaskedImageD(
                bigBox.getWidth(), bigBox.getHeight()
        )
        derivativeImage.setXY0(bigBox.getLLC())
        measMult.expandImageD(bigFP, derivativeImage, numpy.array(analyticDerivative)[n, : ], var)
        ds9.mtv(derivativeImage, frame=n)

        columnVectors.append(partial)
        nonlinearParameters[n] += h
   
    ds9.mtv(modelImage, frame = model.getNonlinearParameterSize()+1)
    
    numericDerivative = numpy.concatenate(columnVectors)
    diff = analyticDerivative - numericDerivative
    factor = analyticDerivative / numericDerivative

    print >> sys.stderr, "analytic:\n%s"% analyticDerivative
    print >> sys.stderr, "numeric:\n%s"% numericDerivative
    print >> sys.stderr, "diff:\n%s"%diff
    print >> sys.stderr, "factor:\n%s"%factor

    
if __name__== "__main__":
    main()