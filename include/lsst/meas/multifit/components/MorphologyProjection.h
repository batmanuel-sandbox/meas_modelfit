// -*- lsst-c++ -*-

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
 
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H

#include <boost/shared_ptr.hpp>

#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModelProjection;

namespace components {

class Morphology;

/**
 * A projection of a Morphology object.
 *
 * A MorphologyProjection should only exist as a data member of a 
 * ComponentModel.
 */
class MorphologyProjection : private boost::noncopyable {
public:
    typedef boost::shared_ptr<MorphologyProjection> Ptr;
    typedef boost::shared_ptr<MorphologyProjection const> ConstPtr;
    
    typedef Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> ParameterJacobianMatrix;
    typedef Eigen::Matrix<Pixel,Eigen::Dynamic, 4> TransformJacobianMatrix;
    typedef boost::shared_ptr<ParameterJacobianMatrix const> ParameterJacobianMatrixPtr;
    typedef boost::shared_ptr<TransformJacobianMatrix const> TransformJacobianMatrixPtr;

    virtual ~MorphologyProjection() {}

    /// Return the Morphology this is a projection of.
    boost::shared_ptr<Morphology const> getMorphology() const { 
        return _morphology; 
    }

    /// Return the kernel dimensions this MorphologyProjection expects.
    lsst::afw::geom::Extent2I const & getKernelDimensions() const { 
        return _kernelDimensions; 
    }

    /**
     * AffineTransform that transforms this projection to the global coordinate 
     * space.
     */
    lsst::afw::geom::AffineTransform::ConstPtr getTransform() const { 
        return _transform; 
    }

    /**
     *  Return the matrix that maps the output of 
     *  ComponentModelProjection::computeProjectedParameterDerivative to the 
     *  morphology block of the nonlinear parameter derivative.
     */
    virtual ParameterJacobianMatrixPtr computeProjectedParameterJacobian() const = 0;

    /**
     *  Return the matrix that deprojects the output of 
     *  ComponentModelProjection::computeProjectedParameterDerivative to the 
     *  morphology terms of the WCS parameter derivative.
     */
    virtual TransformJacobianMatrixPtr computeTransformParameterJacobian() const = 0;

    int const getLinearParameterSize() const; 
    int const getNonlinearParameterSize() const; 
protected:

    /**
     * Handle a change in the linear parameters, as propogated by the owning
     * ComponentModelProjection.
     */
    virtual void _handleLinearParameterChange() {}

    /**
     * Handle a change in the nonlinear morphology parameters, as propogated 
     * by the owning ComponentModelProjection.
     */
    virtual void _handleNonlinearParameterChange() {}

    /**
     * Construct a MorphologyProjection.
     */
    MorphologyProjection(
        boost::shared_ptr<Morphology const> const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions, 
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) : _morphology(morphology), 
        _kernelDimensions(kernelDimensions), 
        _transform(transform) {}

private:
    friend class Morphology;
    friend class lsst::meas::multifit::ComponentModelProjection;

    boost::shared_ptr<Morphology const> _morphology;
    lsst::afw::geom::Extent2I _kernelDimensions;
    lsst::afw::geom::AffineTransform::ConstPtr _transform;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H
