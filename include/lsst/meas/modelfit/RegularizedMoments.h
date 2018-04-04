// -*- lsst-c++ -*-
/*
 * This file is part of package_name.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_MEAS_MODELFIT_REGULARIZEDMOMENTS_H
#define LSST_MEAS_MODELFIT_REGULARIZEDMOMENTS_H

#include <Eigen/Dense>

namespace lsst {
namespace meas {
namespace modelfit {

/**
 *  @brief This class is a model for calculating observed shape moments and their derivatives
 *
 *  This class models the biased weighted moments of an astronomical source given a vector of
 *  moments intrinsic to an object, and a vector of moments representing a weighting function.
 *
 *  The class is constructed with the moments of weight function, and the moments of the object
 *  are set by calling the at method. This allows a single instance to be reused to evaluate
 *  multiple positions in object moment space. This structure is computationally more efficient
 *  as well since the weight function would remain constant for many evaluation of image moments 
 *
 *  This class is designed to be used as an objective function in an optimizer which will find the
 *  maximally likely image moments given moments measured from an image. As such the class
 *  provides two methods in addition to the at method already mentioned. One method computes the
 *  biased weighted moments given the specified weight and "true" moments, and the other computes
 *  the derivative of the objective function along each of the image moment variables.
 */
class MomentsModel {
public:
    using Moments = Eigen::Matrix<double, 6, 1>;
    using Jacobian = Eigen::Matrix<double, 6, 6>;
    using FirstMoment = Eigen::Matrix<double, 2, 1>;
    using SecondMoment = Eigen::Matrix<double, 2, 2>;

    /**
     * Construct a MomentsModel object with a given vector of weight function moments
     *
     * @param[in]   W       An eigen vector of the moments of the elliptical Gaussian weighting function.
     *                      Most likely the first element (the zeroth moment) should be one, representing
     *                      a normalized weight function.
     */
    MomentsModel(Moments const & W): W(W) {
    }

    /**
     *  Sets the moments corresponding to the astronomical object's intrinsic moments
     *
     *  When this function is called, it sets the astronomical object's intrinsic moments for all subsequent
     *  calculations. Additionally after setting the moment, this function calculates and caches the value
     *  of the biased weighted moments given the weights and intrinsic moments.
     *
     *  @param[in]  Q       An eigen vector of the intrinsic moments of the astronomical object
     */
    void at(Moments const & Q);

    /**
     *  Returns the value of the biased weighted moments given the input weight and object's intrinsic moments
     *
     *  @param[out]         Eigen vector of the model biased weighted moments given the input weight and
                            intrinsic moments
     */
    Moments computeValues();

    /**
     *  Computes and returns the gradient of the model biased weighted moments along each of the intrinsic
     *  moment dimensions.
     *
     *  @param[out]         A six by six Eigen matrix, with the biased weighted moments down the column, and
     *                      the derivative of the column moment with respect to the intrinsic moments along
     *                      the rows.
     */
    Jacobian computeJacobian();

private:
    // Calculates the Model evaluated at the input intrinsic moments given the specified weighted moments
    void makeValue();

    // Storage variables
    Moments W, Q, value;

    double constant;

    FirstMoment alpha;

    SecondMoment beta;
};

class ShapePrior {
public:
    using Moments = MomentsModel::Moments;
    using SecondMoment = MomentsModel::SecondMoment;

    virtual void at(SecondMoment second, double flux) = 0;

    virtual double computeLogProbability() const = 0;

    virtual Moments computeLogDerivative() const = 0;
};


class ClassificationPrior {
public:
    virtual double pStar(double flux) const = 0;
};

class StarGalaxyPrior : ShapePrior {
public:
    void at(SecondMoment second, double flux) override;

    double computeLogProbability() const override;

    double Moments computeLogDerivative() const override;

    StarGalaxyPrior(std::unique_ptr<ShapePrior> galaxyPrior, std::unique_ptr<ClassificationPrior> classify,
                    SecondMoment psfMoments, double psfFuzz):_galaxy(std::move(galaxyPrior)), 
                                                             _classification(std::move(classify)),
                                                             _psf(psfMoments), _psfFuzz(psfFuzz){
    };


private:
        std::unique_ptr<ShapePrior> _galaxy;
        std::unique_ptr<ClassificationPrior _classification;
        SecondMoment _psf;
        double _psfFuzz;
};

// Tests for classes in anonymous name spaces
bool testConstant(double tol=1e-6);

bool testAlphaX(double tol=1e-6);
bool testAlphaY(double tol=1e-6);

bool testBetaX(double tol=1e-6);
bool testBetaY(double tol=1e-6);
bool testBetaXY(double tol=1e-6);
}}} // Close namespace lsst::meas::modelfit

#endif // LSST_MEAS_MODELFIT_REGULARIZEDMOMENTS_H

