// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_psf_h_INCLUDED
#define LSST_MEAS_MULTIFIT_psf_h_INCLUDED

#include "boost/scoped_ptr.hpp"

#include "lsst/pex/config.h"

#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/Likelihood.h"
#include "lsst/meas/multifit/optimizer.h"

namespace lsst { namespace meas { namespace multifit {

/// Control object used to define one piece of multishapelet fit to a PSF model; see PsfFitter.
class PsfFitterComponentControl {
public:

    PsfFitterComponentControl(int order_=0, double radiusFactor_=1.0) :
        order(order_), positionPriorSigma(0.0), ellipticityPriorSigma(0.0),
        radiusFactor(radiusFactor_), radiusPriorSigma(0.0)
    {}

    LSST_CONTROL_FIELD(
        order, int,
        "shapelet order for this component; negative to disable this component completely"
    );
    LSST_CONTROL_FIELD(
        positionPriorSigma, double,
        "sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, "
        "relative to the center of the PSF image"
    );
    LSST_CONTROL_FIELD(
        ellipticityPriorSigma, double,
        "sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta"
    );
    LSST_CONTROL_FIELD(
        radiusFactor, double,
        "Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either "
        "the second-moments radius of the PSF image (in an initial fit), or the radius of the primary "
        "component in a previous fit.  Ignored if the previous fit included this component (as then we "
        "can just use that radius)."
    );
    LSST_CONTROL_FIELD(
        radiusPriorSigma, double,
        "sigma in a Gaussian prior on ln(radius/fiducialRadius)"
    );

};

/// Control object used to configure a multishapelet fit to a PSF model; see PsfFitter.
class PsfFitterControl {
public:

    PsfFitterControl() : inner(-1, 0.5), primary(0, 1.0), wings(0, 2.0), outer(-1, 4.0) {}

    LSST_NESTED_CONTROL_FIELD(
        inner, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Innermost shapelet expansion, used to fit PSFs with very sharp cores"
    );

    LSST_NESTED_CONTROL_FIELD(
        primary, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Primary shapelet expansion, typically used to fit the bulk of the PSF "
    );

    LSST_NESTED_CONTROL_FIELD(
        wings, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Wing shapelet expansion (between primary and outer), typically used to fit the wings of the PSF"
    );

    LSST_NESTED_CONTROL_FIELD(
        outer, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Outermost shapelet expansion, used to fit PSFs with very broad wings"
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.multifit.multifitLib, OptimizerControl,
        "Configuration of the optimizer used to do the fitting"
    );

};

/**
 *  @brief Class for fitting multishapelet models to PSF images
 *
 *  This class fits up to four shapelet expansions simultaneously to a PSF image, with the relative radii
 *  and number of shapelet coefficients for each expansion separately configurable.  These expansions are
 *  also named; this allows us to map different fits with some expansions disabled to each other, in order
 *  to first fit an approximate model and follow this up with a more complete model, using the approximate
 *  model as a starting point.
 *
 *  The configuration also defines a simple Bayesian prior for the fit, defined using simple independent
 *  Gaussians for the ellipse parameters of each component.  The priors can be disabled by setting their
 *  width (xxPriorSigma in the control object) to infinity, and those parameters can be held fixed at
 *  their input values by setting the prior width to zero.  The priors are always centered at the input
 *  value, meaning that it may be more appropriate to think of the priors as a form of regularization,
 *  rather than a rigorous prior.  In fact, it's impossible to use a prior here rigorously without a
 *  noise model for the PSF image, which is something the LSST Psf class doesn't provide, and here is
 *  just provided as a constant noise sigma to be provided by the user (who generally just has to chose
 *  a small number arbitrarily).  Decreasing the noise sigma will of course decrease the effect of the
 *  priors (and vice versa).  In any case, having some sort of regularization is probably a good idea,
 *  as this is a very high-dimensional fit.
 */
class PsfFitter {
public:

    /// Initialize the fitter class with the given control object.
    PsfFitter(PsfFitterControl const & ctrl);

    /**
     *  Return the Model object that corresponds to the configuration.
     *
     *  In addition to the shapelet coefficients (stored in the "amplitudes" array), this Model
     *  stores all the initial ellipse parameters in the "fixed" array, as these are used to
     *  define the center of the prior; the "nonlinear" parameters are the free-to-vary ellipse
     *  parameters minus the corresponding initial values.
     */
    PTR(Model) getModel() const { return _model; }

    /**
     *  Return the Prior object that corresponds to the configuration.
     *
     *  This Prior class only supports evaluate() and evaluateDerivatives(), reflecting the fact
     *  that we only intend to use it with a Optimizer, not a Sampler.
     */
    PTR(Prior) getPrior() const { return _prior; }

    /**
     *  Adapt a differently-configured previous fit to be used as an starting point for this PsfFitter.
     *
     *  @param[in] previousFit     The return value of apply() from a differently-configured
     *                             instance of PsfFitter.
     *  @param[in] previousModel   The Model associated with the PsfFitter used to create previousFit.
     *
     *  @return a new MultiShapelet function that may be passed directly to apply().  When possible,
     *  the ellipse and shapelet coefficeints will be copied from previousFit; higher-order coefficients
     *  will be set to zero, and any components used in this but unused in the previous fit will have their
     *  ellipses set relative to the previous fit's "primary" component.
     */
    shapelet::MultiShapeletFunction adapt(
        shapelet::MultiShapeletFunction const & previousFit,
        PTR(Model) previousModel
    ) const;

    /**
     *  Perform an initial fit to a PSF image.
     *
     *  @param[in]  image       The image to fit, typically the result of Psf::computeKernelImage().  The
     *                          image's xy0 should be set such that the center of the PSF is at (0,0).
     *  @param[in]  noiseSigma  An estimate of the noise in the image.  As LSST PSF images are generally
     *                          assumed to be noise-free, this is really just a fiddle-factor for the user.
     *  @param[in]  moments     Second moments of the PSF, typically result of Psf::computeShape() or running
     *                          some other adaptive moments code on the PSF image.  This will be used to
     *                          set the initial ellipses of the multishapelet model.
     */
    shapelet::MultiShapeletFunction apply(
        afw::image::Image<Pixel> const & image,
        Scalar noiseSigma,
        afw::geom::ellipses::Quadrupole const & moments
    ) const;

    /**
     *  Perform a fit to a PSF image, using a previous fit as a starting point
     *
     *  @param[in]  image       The image to fit, typically the result of Psf::computeKernelImage().  The
     *                          image's xy0 should be set such that the center of the PSF is at (0,0).
     *  @param[in]  noiseSigma  An estimate of the noise in the image.  As LSST PSF images are generally
     *                          assumed to be noise-free, this is really just a fiddle-factor for the user.
     *  @param[in]  initial     The result of a previous call to apply(), using an identically-configured
     *                          PsfFitter instance.  To use a result from a differently-configured PsfFitter,
     *                          use adapt().
     */
    shapelet::MultiShapeletFunction apply(
        afw::image::Image<Pixel> const & image,
        Scalar noiseSigma,
        shapelet::MultiShapeletFunction const & initial
    ) const;

private:
    PsfFitterControl _ctrl;
    PTR(Model) _model;
    PTR(Prior) _prior;
};

/**
 *  Likelihood object used to fit multishapelet models to PSF model images; mostly for internal use
 *  by PsfFitter.
 */
class MultiShapeletPsfLikelihood : public Likelihood {
public:

    MultiShapeletPsfLikelihood(
        ndarray::Array<Pixel const,2,1> const & image,
        afw::geom::Point2I const & xy0,
        PTR(Model) model,
        Scalar sigma,
        ndarray::Array<Scalar const,1,1> const & fixed
    );

    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        bool doApplyWeights=true
    ) const;

    virtual ~MultiShapeletPsfLikelihood();

private:
    class Impl;
    boost::scoped_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_psf_h_INCLUDED
