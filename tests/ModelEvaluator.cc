// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ModelEvaluator

#include "boost/pointer_cast.hpp"
#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/deprecated.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/ModelFactory.h"

namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;
namespace image = lsst::afw::image;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;

image::Wcs::Ptr makeWcs(geom::PointD const & crVal) {
    geom::PointD crPix= geom::makePointD(5312.4, 9609);
    Eigen::Matrix2d cdMatrix;
    cdMatrix << 0.00005194, 0.0, 0.0, -0.00005194;
    return boost::make_shared<image::Wcs>(crVal, crPix, cdMatrix);
}

image::Exposure<float>::Ptr makeCharExp(
    multifit::Model::Ptr model, geom::PointD const & crVal
) {
    PTR(image::Wcs) wcs = makeWcs(crVal);

    PTR(det::Psf) psf = det::createPsf("DoubleGaussian", 23, 23, 2);
    det::Footprint::Ptr fp(model->computeProjectionFootprint(psf, wcs));
    image::BBox bbox = fp->getBBox();

    image::MaskedImage<float> mi(
        bbox.getWidth(), 
        bbox.getHeight() 
    );
    mi.setXY0(bbox.getX0(), bbox.getY0());
    multifit::ModelProjection::Ptr projection(model->makeProjection(psf, wcs, fp));
    ndarray::Array<multifit::Pixel const, 1, 1> modelImage(projection->computeModelImage());
    ndarray::Array<multifit::Pixel, 1 ,1> variance(ndarray::allocate(ndarray::makeVector(fp->getNpix())));
    variance = 0.25;

    multifit::expandImage(*fp, mi, modelImage, variance);
    image::Exposure<float>::Ptr exp(new image::Exposure<float>(mi, *wcs));
    exp->setPsf(psf);

    return exp;
}

BOOST_AUTO_TEST_CASE(PsModel) {
    double flux = 1;
    geom::Point2D pixel = geom::makePointD(45,45);

    geom::PointD crVal=geom::makePointD(150.11883, 2.20639);
    PTR(image::Wcs) wcs0 = makeWcs(crVal);   

    lsst::afw::coord::Coord::Ptr coord = wcs0->pixelToSky(pixel);

    multifit::Model::Ptr model = 
        multifit::ModelFactory::createPointSourceModel(
            flux, 
            coord->getPosition(lsst::afw::coord::DEGREES)
        );

    std::list<image::Exposure<float>::Ptr> exposureList;
    exposureList.push_back(makeCharExp(model, crVal));
    
    crVal = geom::makePointD(150.11863, 2.20583);
    exposureList.push_back(makeCharExp(model, crVal));

    crVal= geom::makePointD(150.11917, 2.20639);
    exposureList.push_back(makeCharExp(model, crVal));

    multifit::ModelEvaluator evaluator(model);
    evaluator.setExposureList<float, image::MaskPixel, image::VariancePixel>(exposureList);

    BOOST_CHECK_EQUAL(evaluator.getNProjections(), 3);    
    BOOST_CHECK(evaluator.getNPixels() > 0);

    ndarray::Array<multifit::Pixel const, 1, 1> img;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> modelImage;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> lpd, npd;

    ndarray::shallow(img) = evaluator.getDataVector();
    modelImage = evaluator.computeModelImage();
    lpd = evaluator.computeLinearParameterDerivative();
    npd = evaluator.computeNonlinearParameterDerivative();
    
    //test for nan's in matrices
    for (int i = 0; i < evaluator.getNPixels(); ++i){
        BOOST_CHECK_EQUAL(img[i], img[i]);
        BOOST_CHECK_EQUAL(modelImage[i], modelImage[i]);

        for (int j = 0; j < evaluator.getLinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(lpd(i, j), lpd(i,j));

        for (int j = 0; j < evaluator.getNonlinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(npd(i,j), npd(i,j));
    }
}
BOOST_AUTO_TEST_CASE(SersicModel) {
    //define the ellipse parameters in pixel coordinates
    double flux = 1;
    geom::Point2D pixel = geom::makePointD(45,45);
    geom::ellipses::Axes axes(3,5,0);

    geom::PointD crVal=geom::makePointD(150.11883, 2.20639);
    PTR(image::Wcs) wcs0 = makeWcs(crVal);   

    lsst::afw::coord::Coord::Ptr coord = wcs0->pixelToSky(pixel);

    multifit::Cache::ConstPtr cache;

    try {
        cache = multifit::Cache::load("testCache", "Sersic", false);
    } catch (...){
        lsst::pex::policy::Policy pol;
        cache = multifit::makeSersicCache(pol);
    }
    multifit::components::SersicMorphology::setSersicCache(cache);

    //transform the ellipse parameters to be in sky coordinates
    geom::AffineTransform transform = wcs0->linearizePixelToSky(
        coord->getPosition(lsst::afw::coord::DEGREES)
    );
    axes.transform(transform).inPlace();

    multifit::Model::Ptr model = multifit::ModelFactory::createSersicModel(
        flux, 
        coord->getPosition(lsst::afw::coord::DEGREES),
        axes, 
        1.0
    );

    std::list<image::Exposure<float>::Ptr> exposureList;
    exposureList.push_back(makeCharExp(model, crVal));
    
    crVal = geom::makePointD(150.11863, 2.20583);
    exposureList.push_back(makeCharExp(model, crVal));

    crVal= geom::makePointD(150.11917, 2.20639);
    exposureList.push_back(makeCharExp(model, crVal));

    multifit::ModelEvaluator evaluator(model);
    evaluator.setExposureList<float, image::MaskPixel, image::VariancePixel>(exposureList);

    BOOST_CHECK_EQUAL(evaluator.getNProjections(), 3);    
    BOOST_CHECK(evaluator.getNPixels() > 0);
    ndarray::Array<multifit::Pixel const, 1, 1> img;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> modelImage;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> lpd, npd;

    ndarray::shallow(img) = evaluator.getDataVector();
    modelImage = evaluator.computeModelImage();
    lpd = evaluator.computeLinearParameterDerivative();
    npd = evaluator.computeNonlinearParameterDerivative();
    
    //test for nan's in matrices
    for (int i = 0; i < evaluator.getNPixels(); ++i){
        BOOST_CHECK_EQUAL(img[i], img[i]);
        BOOST_CHECK_EQUAL(modelImage[i], modelImage[i]);

        for (int j = 0; j < evaluator.getLinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(lpd(i, j), lpd(i,j));

        for (int j = 0; j < evaluator.getNonlinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(npd(i,j), npd(i,j));
    }

    
}