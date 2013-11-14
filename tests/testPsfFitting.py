#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import unittest
import os
import numpy
try:
    import scipy.integrate
except ImportError:
    scipy = None

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.meas.multifit

numpy.random.seed(500)

class PsfFittingTestCase(lsst.utils.tests.TestCase):

    def testMultiShapeletPsfModel(self):
        orders = [2, 3, 4]
        model = lsst.meas.multifit.makeMultiShapeletPsfModel(orders)
        self.assertEqual(model.getNonlinearDim(), 5*len(orders))
        self.assertEqual(model.getAmplitudeDim(), sum(lsst.shapelet.computeSize(order) for order in orders))
        self.assertEqual(model.getFixedDim(), 0)
        self.assertEqual(model.getBasisCount(), len(orders))
        for order, component in zip(orders, lsst.meas.multifit.MultiModel.cast(model).getComponents()):
            self.assertEqual(component.getNonlinearDim(), 5)
            self.assertEqual(component.getAmplitudeDim(), lsst.shapelet.computeSize(order))
            self.assertEqual(component.getFixedDim(), 0)
            self.assertEqual(component.getBasisCount(), 1)

    def testSingleGaussianFit(self):
        orders = [0]
        model = lsst.meas.multifit.makeMultiShapeletPsfModel(orders)
        nonlinear = numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar)
        amplitudes = numpy.ones(model.getAmplitudeDim(), dtype=lsst.meas.multifit.Scalar)
        fixed = numpy.zeros(model.getFixedDim(), dtype=lsst.meas.multifit.Scalar)
        ellipses = model.makeEllipseVector()
        ellipses[0].setCore(lsst.afw.geom.ellipses.Axes(4.0, 4.0, 0.0))
        model.readEllipses(ellipses, nonlinear, fixed)
        msf = model.makeShapeletFunction(nonlinear, amplitudes, fixed)
        image = numpy.zeros((25, 25), dtype=float)
        xy0 = lsst.afw.geom.Point2I(-12, -12)
        msf.evaluate().addToImage(image, xy0)
        image = image.astype(lsst.meas.multifit.Pixel)
        likelihood = lsst.meas.multifit.MultiShapeletPsfLikelihood(image, xy0, model, fixed)
        self.assertEqual(likelihood.getDataDim(), image.size)
        self.assertEqual(likelihood.getAmplitudeDim(), model.getAmplitudeDim())
        self.assertEqual(likelihood.getNonlinearDim(), model.getNonlinearDim())
        self.assertEqual(likelihood.getFixedDim(), model.getFixedDim())
        self.assertClose(likelihood.getData(), image.flat)
        matrix = numpy.zeros((likelihood.getAmplitudeDim(), likelihood.getDataDim()),
                             dtype=lsst.meas.multifit.Pixel).transpose()
        likelihood.computeModelMatrix(matrix, nonlinear)
        self.assertClose(image, matrix.reshape(image.shape), rtol=1E-8, atol=1E-8, plotOnFailure=True)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFittingTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
