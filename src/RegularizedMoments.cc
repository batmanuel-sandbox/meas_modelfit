// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018  AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include <cmath>
#include <exception>
#include <string>
#include <tuple>

#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/modelfit/RegularizedMoments.h"

namespace lsst {
namespace meas {
namespace modelfit {

namespace {
/*
 *  Anonymous namespace to contain code that calculates all the intermediate products used in
 *  calculating the biased waited moments model value and derivative. See corresponding document
 *  for explanation of parameter names.
 */

using Moments = MomentsModel::Moments; 
using FirstMoment = MomentsModel::FirstMoment;
using SecondMoment = MomentsModel::SecondMoment;

// Builds moments with specific values used in the tests for the anonymous namesapce functionality
std::tuple<Moments, Moments> buildTestMoments(){
    Moments Q, W;
    Q << 6, 4, 3, 2, 1, 4;
    W << 2, 4, 3.1, 2.5, 1.2, 3.7;
    return std::make_tuple(Q, W);
}

struct AlphaX {
    static double computeValue(Moments const & Q, Moments const & W) {
        double valueTop = std::pow(Q[4],2)*W[1]+(std::pow(W[4],2)-W[3]*W[5])*Q[1]-(Q[2]*W[4]-W[2]*W[4]+
                          W[1]*W[5])*Q[3]+(Q[2]*W[3]-W[2]*W[3]+Q[1]*W[4]+W[1]*W[4])*Q[4]-(Q[3]*W[1]+
                          Q[1]*W[3])*Q[5];
        double valueBottom = std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-
                             W[3]*W[5];
        return valueTop/valueBottom;
    }

    static Moments computeGradient(Moments const & Q, Moments const & W, double muBottom, double sigBottom) {
        Moments vec;
        vec(0, 0) = 0;
        vec(1, 0) = -1*(Q[5]*W[3]-Q[4]*W[4]-std::pow(W[4],2)+W[3]*W[5])/muBottom;

        vec(2, 0) = (Q[4]*W[3]-Q[3]*W[4])/muBottom;

        double vec3Top = W[2]*std::pow(W[4],3)+W[1]*W[3]*std::pow(W[5],2)-(Q[2]*W[4]-W[2]*W[4])*
                         std::pow(Q[4],2)-(Q[1]*W[3]-W[1]*W[3])*std::pow(Q[5],2)+(std::pow(W[4],2)*W[5]-
                         W[3]*std::pow(W[5],2))*Q[1]-(std::pow(W[4],3)-W[3]*W[4]*W[5])*Q[2]+
                         (2*W[2]*std::pow(W[4],2)+Q[1]*W[4]*W[5]-
                         (2*std::pow(W[4],2)-W[3]*W[5])*Q[2]-(W[2]*W[3]+W[1]*W[4])*W[5])*Q[4]+
                         (Q[2]*W[3]*W[4]-W[2]*W[3]*W[4]-W[1]*std::pow(W[4],2)+2*W[1]*W[3]*W[5]+
                         (std::pow(W[4],2)-2*W[3]*W[5])*Q[1]+(Q[2]*W[3]-W[2]*W[3]+Q[1]*W[4]-
                         W[1]*W[4])*Q[4])*Q[5]-(W[2]*W[3]*W[4]+W[1]*std::pow(W[4],2))*W[5];

        vec(3, 0) = vec3Top/sigBottom;

        double vec4Top = W[2]*W[3]*std::pow(W[4],2)-W[1]*std::pow(W[4],3)+(Q[2]*W[3]-W[2]*W[3]+
                         Q[1]*W[4]-W[1]*W[4])*std::pow(Q[4],2)+(std::pow(W[4],3)-W[3]*W[4]*W[5])*Q[1]-
                         (W[3]*std::pow(W[4],2)-std::pow(W[3],2)*W[5])*Q[2]+
                         (2*W[2]*std::pow(W[4],2)+Q[1]*W[4]*W[5]-(2*std::pow(W[4],2)-W[3]*W[5])*Q[2]-
                         (W[2]*W[3]+W[1]*W[4])*W[5])*Q[3]-2*(W[1]*std::pow(W[4],2)-W[1]*W[3]*W[5]-
                         (std::pow(W[4],2)-W[3]*W[5])*Q[1]+(Q[2]*W[4]-W[2]*W[4])*Q[3])*Q[4]+
                         (Q[2]*std::pow(W[3],2)-W[2]*std::pow(W[3],2)-Q[1]*W[3]*W[4]+W[1]*W[3]*W[4]+
                         (Q[2]*W[3]-W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[3]-2*(Q[1]*W[3]-
                         W[1]*W[3])*Q[4])*Q[5]-(W[2]*std::pow(W[3],2)-W[1]*W[3]*W[4])*W[5];

        vec(4, 0) = -1*vec4Top/sigBottom;

        double vec5Top = (Q[2]*W[4]-W[2]*W[4])*std::pow(Q[3],2)+(Q[1]*W[3]-W[1]*W[3])*std::pow(Q[4],2)+
                         (Q[2]*W[3]*W[4]-W[2]*W[3]*W[4]-Q[1]*std::pow(W[4],2)+W[1]*std::pow(W[4],2))*Q[3]-
                         (Q[2]*std::pow(W[3],2)-W[2]*std::pow(W[3],2)-Q[1]*W[3]*W[4]+W[1]*W[3]*W[4]+
                         (Q[2]*W[3]-W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[3])*Q[4];

        vec(5, 0) = -1*vec5Top/sigBottom;
        return vec;
    }
};

struct AlphaY {
    static double computeValue(Moments const & Q, Moments const & W) {
        double valueTop = std::pow(Q[4],2)*W[2]-Q[2]*Q[3]*W[5]+(std::pow(W[4],2)-W[3]*W[5])*Q[2]+
                          (Q[2]*W[4]+W[2]*W[4]+Q[1]*W[5]-W[1]*W[5])*Q[4]-(Q[3]*W[2]+W[2]*W[3]+Q[1]*W[4]-
                          W[1]*W[4])*Q[5];

        double valueBottom = std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-
                             W[3]*W[5];

        return valueTop/valueBottom;
    }

    static Moments computeGradient(Moments const & Q, Moments const & W, double muBottom, double sigBottom) {
        Moments vec;
        vec(0, 0) = 0;

        vec(1, 0) = (Q[4]*W[5] - Q[5]*W[4])/muBottom;

        vec(2, 0) = (Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5])/muBottom;

        double vec3Top = (Q[2]*W[5]-W[2]*W[5])*std::pow(Q[4],2)+(Q[1]*W[4]-W[1]*W[4])*std::pow(Q[5],2)+
                         (Q[2]*W[4]*W[5]-W[2]*W[4]*W[5]-Q[1]*std::pow(W[5],2)+W[1]*std::pow(W[5],2))*Q[4]-
                         (Q[2]*std::pow(W[4],2)-W[2]*std::pow(W[4],2)-Q[1]*W[4]*W[5]+W[1]*W[4]*W[5]+
                         (Q[2]*W[4]-W[2]*W[4]+Q[1]*W[5]-W[1]*W[5])*Q[4])*Q[5];

        vec(3, 0) = -1*vec3Top/sigBottom;

        double vec4Top = W[2]*std::pow(W[4],3)+W[1]*W[3]*std::pow(W[5],2)-(Q[2]*W[4]-W[2]*W[4]+Q[1]*W[5]-
                         W[1]*W[5])*std::pow(Q[4],2)+(std::pow(W[4],2)*W[5]-W[3]*std::pow(W[5],2))*Q[1]-
                         (std::pow(W[4],3)-W[3]*W[4]*W[5])*Q[2]+
                         (Q[2]*W[4]*W[5]-W[2]*W[4]*W[5]-Q[1]*std::pow(W[5],2)+W[1]*std::pow(W[5],2))*Q[3]+
                         2*(W[2]*std::pow(W[4],2)-W[2]*W[3]*W[5]-(std::pow(W[4],2)-W[3]*W[5])*Q[2]+
                         (Q[2]*W[5]-W[2]*W[5])*Q[3])*Q[4]-(Q[2]*W[3]*W[4]-W[2]*W[3]*W[4]+
                         2*W[1]*std::pow(W[4],2)-W[1]*W[3]*W[5]-(2*std::pow(W[4],2)-W[3]*W[5])*Q[1]+
                         (Q[2]*W[4]-W[2]*W[4]+Q[1]*W[5]-W[1]*W[5])*Q[3]-2*(Q[1]*W[4]-W[1]*W[4])*Q[4])*Q[5]-
                         (W[2]*W[3]*W[4]+W[1]*std::pow(W[4],2))*W[5];

        vec(4, 0) = vec4Top/sigBottom;

        double vec5Top = W[2]*W[3]*std::pow(W[4],2)-W[1]*std::pow(W[4],3)+(Q[2]*W[5]-W[2]*W[5])*
                         std::pow(Q[3],2)+(Q[1]*W[4]-W[1]*W[4])*std::pow(Q[4],2)+(std::pow(W[4],3)-
                         W[3]*W[4]*W[5])*Q[1]-(W[3]*std::pow(W[4],2)-std::pow(W[3],2)*W[5])*Q[2]+
                         (W[2]*std::pow(W[4],2)-Q[1]*W[4]*W[5]-(std::pow(W[4],2)-2*W[3]*W[5])*Q[2]-
                         (2*W[2]*W[3]-W[1]*W[4])*W[5])*Q[3]-(Q[2]*W[3]*W[4]-W[2]*W[3]*W[4]+2*W[1]*
                         std::pow(W[4],2)-W[1]*W[3]*W[5]-(2*std::pow(W[4],2)-W[3]*W[5])*Q[1]+
                         (Q[2]*W[4]-W[2]*W[4]+Q[1]*W[5]-W[1]*W[5])*Q[3])*Q[4]-(W[2]*std::pow(W[3],2)-
                         W[1]*W[3]*W[4])*W[5];

        vec(5, 0) = -1*vec5Top/sigBottom;
        return vec;
    }
};

struct BetaX {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (std::pow(Q[4],2)*W[3]-Q[3]*Q[5]*W[3]+(std::pow(W[4],2)-W[3]*W[5])*Q[3])/
               (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]);
    }

    static Moments computeGradient(Moments const & Q, Moments const & W, double sigBottom) {
        Moments vec;

        vec(0, 0) = 0;
        vec(1, 0) = 0;
        vec(2, 0) = 0;

        double vec3Top = std::pow(Q[5],2)*std::pow(W[3],2)+std::pow(Q[4],2)*std::pow(W[4],2)+
                         std::pow(W[4],4)-2*W[3]*std::pow(W[4],2)*W[5]+std::pow(W[3],2)*std::pow(W[5],2)+
                         2*(std::pow(W[4],3)-W[3]*W[4]*W[5])*Q[4]-2*(Q[4]*W[3]*W[4]+W[3]*std::pow(W[4],2)-
                         std::pow(W[3],2)*W[5])*Q[5];

        vec(3, 0) = vec3Top/sigBottom;

        double vec4Top = 2*(std::pow(Q[4],2)*W[3]*W[4]-(std::pow(W[4],3)-W[3]*W[4]*W[5])*Q[3]-
                         (Q[3]*std::pow(W[4],2)-W[3]*std::pow(W[4],2)+std::pow(W[3],2)*W[5])*Q[4]-
                         (Q[4]*std::pow(W[3],2)-Q[3]*W[3]*W[4])*Q[5]);

        vec(4, 0) = vec4Top/sigBottom;

        double vec5Top = std::pow(Q[4],2)*std::pow(W[3],2)-2*Q[3]*Q[4]*W[3]*W[4]+
                         std::pow(Q[3],2)*std::pow(W[4],2);

        vec(5, 0) = vec5Top/sigBottom;

        return vec;
    }
};

struct BetaXY {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (std::pow(Q[4],2)*W[4]-Q[3]*Q[5]*W[4]+(std::pow(W[4],2)-W[3]*W[5])*Q[4])/
               (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]);
    }

    static Moments computeGradient(Moments const & Q, Moments const & W, double sigBottom) {
        Moments vec;

        vec(0, 0) = 0;
        vec(1, 0) = 0;
        vec(2, 0) = 0;

        double vec3Top = std::pow(Q[5],2)*W[3]*W[4]+std::pow(Q[4],2)*W[4]*W[5]+(std::pow(W[4],2)*W[5]-
                         W[3]*std::pow(W[5],2))*Q[4]-(std::pow(W[4],3)-W[3]*W[4]*W[5]+(std::pow(W[4],2)+
                         W[3]*W[5])*Q[4])*Q[5];

        vec(3, 0) = vec3Top/sigBottom;

        double vec4Top = std::pow(W[4],4)-2*W[3]*std::pow(W[4],2)*W[5]+std::pow(W[3],2)*std::pow(W[5],2)+
                         (std::pow(W[4],2)+W[3]*W[5])*std::pow(Q[4],2)-(std::pow(W[4],2)*W[5]-
                         W[3]*std::pow(W[5],2))*Q[3]+2*(std::pow(W[4],3)-Q[3]*W[4]*W[5]-W[3]*W[4]*W[5])*Q[4]-
                         (2*Q[4]*W[3]*W[4]+W[3]*std::pow(W[4],2)-std::pow(W[3],2)*W[5]-(std::pow(W[4],2)+
                         W[3]*W[5])*Q[3])*Q[5];

        vec(4, 0) = vec4Top/sigBottom;

        double vec5Top = std::pow(Q[4],2)*W[3]*W[4]+std::pow(Q[3],2)*W[4]*W[5]-(std::pow(W[4],3)-
                         W[3]*W[4]*W[5])*Q[3]+(W[3]*std::pow(W[4],2)-std::pow(W[3],2)*W[5]-
                         (std::pow(W[4],2)+W[3]*W[5])*Q[3])*Q[4];

        vec(5, 0) = vec5Top/sigBottom;

        return vec;
    }
};

struct BetaY {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (std::pow(Q[4],2)*W[5]+(std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5])*Q[5])/
            (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]);
    }

    static Moments computeGradient(Moments const & Q, Moments const & W, double sigBottom) {
        Moments vec;

        vec(0, 0) = 0;
        vec(1, 0) = 0;
        vec(2, 0) = 0;

        double vec3Top = std::pow(Q[5],2)*std::pow(W[4],2)-2*Q[4]*Q[5]*W[4]*W[5]+std::pow(Q[4],2)*
                         std::pow(W[5],2);

        vec(3, 0) = vec3Top/sigBottom;

        double vec4Top = 2*(std::pow(Q[4],2)*W[4]*W[5]+(std::pow(W[4],2)*W[5]-Q[3]*std::pow(W[5],2)-
                         W[3]*std::pow(W[5],2))*Q[4]-(Q[4]*std::pow(W[4],2)+std::pow(W[4],3)-Q[3]*W[4]*W[5]-
                         W[3]*W[4]*W[5])*Q[5]);

        vec(4, 0) = vec4Top/sigBottom;

        double vec5Top = std::pow(Q[4],2)*std::pow(W[4],2)+std::pow(W[4],4)-2*W[3]*std::pow(W[4],2)*W[5]+
                         std::pow(Q[3],2)*std::pow(W[5],2)+std::pow(W[3],2)*std::pow(W[5],2)-2*
                         (std::pow(W[4],2)*W[5]-W[3]*std::pow(W[5],2))*Q[3]+2*(std::pow(W[4],3)-
                         Q[3]*W[4]*W[5]-W[3]*W[4]*W[5])*Q[4];

        vec(5, 0) = vec5Top/sigBottom;

        return vec;
    }
};

struct Constant {
    static double computeValue(Moments Q, Moments W) {
        double expTop = std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+std::pow(Q[1],2)*
                        W[5]+std::pow(W[1],2)*W[5]+2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-
                        W[1]*W[4])*Q[2]+(std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-2*
                        ((Q[1]-W[1])*Q[2]-Q[1]*W[2]+W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-2*Q[1]*W[1]+
                        std::pow(W[1],2))*Q[5];
        
        double expBottom = 2*(std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-
                           W[3]*W[5]);

        double bottom = 2*pi*std::sqrt((-1*std::pow((Q[4]+W[4]),2)+(Q[3]+W[3])*(Q[5]+W[5])));

        return std::exp(expTop/expBottom)/bottom;

    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;

        vec(0,0) = 0;

        double vec1Top = ((Q[2]-W[2])*Q[4]-(Q[1]-W[1])*Q[5]+
                         Q[2]*W[4]-W[2]*W[4]-Q[1]*W[5]+W[1]*W[5])*
                         std::exp((std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+
                         std::pow(Q[1],2)*W[5]+std::pow(W[1],2)*W[5]+
                         2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[2]+
                         (std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-2*((Q[1]-W[1])*Q[2]-Q[1]*W[2]+
                         W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-2*Q[1]*W[1]+std::pow(W[1],2))*Q[5])/
                         (2*(std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-
                         W[3]*W[5])));

        double vec12Bottom = std::sqrt(-4*std::pow(pi,2)*std::pow((Q[4]+W[4]),2)+
                             4*std::pow(pi,2)*(Q[3]+W[3])*(Q[5]+W[5]))*(std::pow(Q[4],2)-(Q[3]+W[3])*
                             Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]);

        vec(1,0) = -1*vec1Top/vec12Bottom;

        double vec2Top = ((Q[2]-W[2])*Q[3]-(Q[1]-W[1])*Q[4]+
                         Q[2]*W[3]-W[2]*W[3]-Q[1]*W[4]+W[1]*W[4])*
                         std::exp((std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+
                         std::pow(Q[1],2)*W[5]+std::pow(W[1],2)*W[5]+
                         2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[2]+
                         (std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-2*((Q[1]-W[1])*Q[2]-Q[1]*W[2]+
                         W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-2*Q[1]*W[1]+std::pow(W[1],2))*Q[5])/
                         (2*(std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-
                         Q[3]*W[5]-W[3]*W[5])));

        vec(2,0) = vec2Top/vec12Bottom; 

        double bottomComp = 4*std::pow(pi,2)*((Q[3]+W[3])*(Q[5]+W[5]) - std::pow((Q[4]+W[4]),2));
        double vec345FirstBottom = std::pow(bottomComp ,1.5);
        double vec345SecondBottom = std::sqrt(bottomComp);

        double vec345ExpTop = std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+
                              std::pow(Q[1],2)*W[5]+std::pow(W[1],2)*W[5]+
                              2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[2]+
                              (std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-2*((Q[1]-W[1])*Q[2]-
                              Q[1]*W[2]+W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-2*Q[1]*W[1]+
                              std::pow(W[1],2))*Q[5];

        double vec345ExpBottom = 2*(std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-
                                 Q[3]*W[5]-W[3]*W[5]);

        double vec345Exp = std::exp(vec345ExpTop/vec345ExpBottom);


        double vec3FirstTop = 2*std::pow(pi,2)*(Q[5]+W[5])*vec345Exp;

        double vec3SecondTopFirst = (std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))/
                                    (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-
                                    Q[3]*W[5]-W[3]*W[5]);

        double vec3SecondTopSecondTop = (std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+
                                        std::pow(Q[1],2)*W[5]+std::pow(W[1],2)*W[5]+
                                        2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[2]+
                                        (std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-
                                        2*((Q[1]-W[1])*Q[2]-Q[1]*W[2]+
                                        W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-2*Q[1]*W[1]+
                                        std::pow(W[1],2))*Q[5])*(Q[5]+W[5]);

        double vec3SecondTopSecondBottom = std::pow((std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+
                                           std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]),2);

        vec(3,0) = -1* vec3FirstTop/vec345FirstBottom +
                   (vec3SecondTopFirst + vec3SecondTopSecondTop/vec3SecondTopSecondBottom*vec345Exp)/
                   (2*vec345SecondBottom);

        double vec4FirstTop = 4*std::pow(pi,2)*(Q[4]+W[4])*vec345Exp;

        
        double vec4SecondTopFirst = ((Q[1]-W[1])*Q[2]-Q[1]*W[2]+W[1]*W[2])/
                                    (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-
                                    Q[3]*W[5]-W[3]*W[5]);

        double vec4SecondTopSecondTop = (std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+
                                        std::pow(Q[1],2)*W[5]+std::pow(W[1],2)*W[5]+
                                        2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-W[1]*W[4])*Q[2]+
                                        (std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-
                                        2*((Q[1]-W[1])*Q[2]-Q[1]*W[2]+W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-
                                        2*Q[1]*W[1]+std::pow(W[1],2))*Q[5])*(Q[4]+W[4]);

        double vec4SecondTopSecondBottom = std::pow((std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+
                                           std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]),2);


        vec(4, 0) = vec4FirstTop/vec345FirstBottom - 
                    (vec4SecondTopFirst + vec4SecondTopSecondTop/vec4SecondTopSecondBottom*vec345Exp)/
                    vec345SecondBottom;

        double vec5FirstTop = 2*std::pow(pi,2)*(Q[3]+W[3])*vec345Exp;

        double vec5SecondTopFirst = (std::pow(Q[1],2)-2*Q[1]*W[1]+std::pow(W[1],2))/
                                    (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-
                                    Q[3]*W[5]-W[3]*W[5]);

        double vec5SecondTopSecondTop = (std::pow(Q[2],2)*W[3]+std::pow(W[2],2)*W[3]-2*W[1]*W[2]*W[4]+
                                        std::pow(Q[1],2)*W[5]+std::pow(W[1],2)*W[5]
                                        +2*(W[2]*W[4]-W[1]*W[5])*Q[1]-2*(W[2]*W[3]+Q[1]*W[4]-
                                        W[1]*W[4])*Q[2]+(std::pow(Q[2],2)-2*Q[2]*W[2]+std::pow(W[2],2))*Q[3]-
                                        2*((Q[1]-W[1])*Q[2]-Q[1]*W[2]+W[1]*W[2])*Q[4]+(std::pow(Q[1],2)-
                                        2*Q[1]*W[1]+std::pow(W[1],2))*Q[5])*(Q[3]+W[3]);

        double vec5SecondTopSecondBottom = std::pow((std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+
                                           std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]),2);

        vec(5, 0) = -1*vec5FirstTop/vec345FirstBottom +
                    (vec5SecondTopFirst+(vec5SecondTopSecondTop/vec5SecondTopSecondBottom))*vec345Exp/
                    (2*vec345SecondBottom);

        return vec;

    }
private:
    static double constexpr pi = lsst::afw::geom::PI;

};

// Creates the alpha first moment vector from its sub components
MomentsModel::FirstMoment makeAlpha(Moments const & Q, Moments const & W) {
    double x = AlphaX::computeValue(Q, W);
    double y = AlphaY::computeValue(Q, W);

    return MomentsModel::FirstMoment(x, y);
};

// Calculates quantities used in calculating the alpha derivatives. These quantities are used in multiple
// places so are calculated once and passed where needed
auto makeAlphaGradDivisors(Moments const & Q, Moments const & W) {
    double muBottom = (std::pow(Q[4],2)-(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5]);
    double sigBottom = std::pow(Q[4],4)+4*std::pow(Q[4],3)*W[4]+std::pow(W[4],4)-2*W[3]*std::pow(W[4],2)*W[5]+
                       std::pow(Q[3],2)*std::pow(W[5],2)+std::pow(W[3],2)*std::pow(W[5],2)+
                       2*(3*std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5])*std::pow(Q[4],2)+(std::pow(Q[3],2)+
                       2*Q[3]*W[3]+std::pow(W[3],2))*std::pow(Q[5],2)-2*(std::pow(W[4],2)*W[5]-W[3]*
                       std::pow(W[5],2))*Q[3]+4*(std::pow(W[4],3)-Q[3]*W[4]*W[5]-W[3]*W[4]*W[5])*Q[4]-
                       2*((Q[3]+W[3])*std::pow(Q[4],2)+W[3]*std::pow(W[4],2)-std::pow(Q[3],2)*W[5]-
                       std::pow(W[3],2)*W[5]+(std::pow(W[4],2)-2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+
                       W[3]*W[4])*Q[4])*Q[5];
    return std::make_tuple(muBottom, sigBottom);
}

// Calculates and returns the derivatives for the alpha components.
auto makeAlphaGrad(Moments const & Q, Moments const & W) {
    double muBottom, sigBottom;
    std::tie(muBottom, sigBottom) = makeAlphaGradDivisors(Q, W);

    Moments x = AlphaX::computeGradient(Q, W, muBottom, sigBottom);
    Moments y = AlphaY::computeGradient(Q, W, muBottom, sigBottom);

    return std::make_tuple(x, y);
}

// Calculates and assembles the matrix of the second moment from its sub components
MomentsModel::SecondMoment makeBeta(Moments const & Q, Moments const & W) {
    double x = BetaX::computeValue(Q, W);
    double xy = BetaXY::computeValue(Q, W);
    double y = BetaY::computeValue(Q, W);

    MomentsModel::SecondMoment beta;
    
    beta(0, 0) = x;
    beta(1, 0) = xy;
    beta(0, 1) = xy;
    beta(1, 1) = y;

    return beta;
}

// Calculates quantities used in calculating the beta derivatives. These quantities are used in multiple
// places so are calculated once and passed where needed
double makeBetaGradDivisor(Moments const & Q, Moments const & W) {
     double bottom = std::pow(Q[4],4)+4*std::pow(Q[4],3)*W[4]+std::pow(W[4],4)-2*W[3]*std::pow(W[4],2)*W[5]+
                        std::pow(Q[3],2)*std::pow(W[5],2)+std::pow(W[3],2)*std::pow(W[5],2)+
                        2*(3*std::pow(W[4],2)-Q[3]*W[5]-W[3]*W[5])*std::pow(Q[4],2)+(std::pow(Q[3],2)+
                        2*Q[3]*W[3]+std::pow(W[3],2))*std::pow(Q[5],2)-2*(std::pow(W[4],2)*W[5]-
                        W[3]*std::pow(W[5],2))*Q[3]+4*(std::pow(W[4],3)-Q[3]*W[4]*W[5]-W[3]*W[4]*W[5])*Q[4]-
                        2*((Q[3]+W[3])*std::pow(Q[4],2)+W[3]*std::pow(W[4],2)-std::pow(Q[3],2)*W[5]-
                        std::pow(W[3],2)*W[5]+(std::pow(W[4],2)-2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*
                        Q[4])*Q[5];
    return bottom;
}

// Calculates and returns the derivatives for the beta components
auto makeBetaGrad(Moments const & Q, Moments const & W) {
    double bottom = makeBetaGradDivisor(Q, W);
    return std::make_tuple(BetaX::computeGradient(Q, W, bottom), BetaXY::computeGradient(Q, W, bottom),
                           BetaY::computeGradient(Q, W, bottom));
}

// Tests if two Moments vectors are approximately equal to some tolerance
bool approxEqual(Moments const & first, Moments const & second, double tol=1e-6){
    for (size_t i = 0; i < 6; ++i){
        if (abs(first(i, 0) - second(i, 0)) > tol) {
            return false;
        }
    }
    return true;
}

} //end anonymous

// Test the various calculations done by objects in the anonymous namespace

bool testAlphaX(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = AlphaX::computeValue(Q, W);
    double muBottom, sigBottom;
    std::tie(muBottom, sigBottom) = makeAlphaGradDivisors(Q, W);
    Moments firstRes = AlphaX::computeGradient(Q, W, muBottom, sigBottom);
    double zeroTruth = 4.00033545790003;
    Moments firstTruth;
    firstTruth << 0,
                  0.557195571955720,
                  -0.00335457900033546,
                  -0.00411214444247764,
                  0.00843596158202444,
                  -0.0000506394012127103;
    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstRes, firstTruth, tol);
}

bool testAlphaY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = AlphaY::computeValue(Q, W);
    double muBottom, sigBottom;
    std::tie(muBottom, sigBottom) = makeAlphaGradDivisors(Q, W);
    Moments firstRes = AlphaY::computeGradient(Q, W, muBottom, sigBottom);
    double zeroTruth = 3.05300234820530;
    Moments firstTruth;
    firstTruth << 0,
                  0.0369003690036900,
                  0.469976517946998,
                  -0.000272327446521693,
                  -0.00291142797372289,
                  0.00709458010990103;

    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstRes, firstTruth, tol);
}

bool testBetaX(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaX::computeValue(Q, W);
    double bottom = makeBetaGradDivisor(Q, W);
    Moments firstRes = BetaX::computeGradient(Q, W, bottom);
    double zeroTruth = 1.11103656491110;
    Moments firstTruth;
    firstTruth << 0,
                  0,
                  0,
                  0.310466905407061,
                  -0.00373831312952513,
                  0.0000112532002694914;

    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstRes, firstTruth, tol);
}

bool testBetaXY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaXY::computeValue(Q, W);
    double bottom = makeBetaGradDivisor(Q, W);
    Moments firstRes = BetaXY::computeGradient(Q, W, bottom);
    double zeroTruth = 0.543777255954378;
    Moments firstTruth;
    firstTruth << 0,
                  0,
                  0,
                  0.0205607222123882,
                  0.261745049520270,
                  -0.00157657335775577;

    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstRes, firstTruth, tol);
}

bool testBetaY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaY::computeValue(Q, W);
    double bottom = makeBetaGradDivisor(Q, W);
    Moments firstRes = BetaY::computeGradient(Q, W, bottom);
    double zeroTruth = 1.91680644079168;
    Moments firstTruth;
    firstTruth << 0,
                  0,
                  0,
                  0.00136163723260849,
                  0.0346846138706271,
                  0.220877927421585;

    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstRes, firstTruth, tol);
}


bool testConstant(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = Constant::computeValue(Q, W);
    Moments firstRes = Constant::computeGradient(Q, W);
    double zeroTruth = 0.0915084542604366/afw::geom::PI;
    Moments firstTruth;
    firstTruth << 0,
                  -0.000675339145833491/afw::geom::PI,
                  0.00138137552556848/afw::geom::PI,
                  -0.0118159430257175/afw::geom::PI,
                  0.00674319680500958/afw::geom::PI,
                  -0.00689645127785071/afw::geom::PI;

    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstRes, firstTruth, tol);
}

// Sets the location in moment space to evaluate
void MomentsModel::at(Moments const & inputQ) {
   Q = Moments(inputQ); 
   makeValue();
}

// Calculates the value of the weighted biased moments
void MomentsModel::makeValue() {
    alpha = makeAlpha(Q, W);
    beta = makeBeta(Q, W);
    constant = Constant::computeValue(Q, W);

    double zero = Q[0]*constant;
    FirstMoment one = zero*alpha;
    SecondMoment two = zero*(beta + alpha*alpha.transpose());

    value << zero, one[0], one[1], two(0, 0), two(0, 1), two(1, 1);
}

// Returns the value of the weighted biased moments calculation
Moments MomentsModel::computeValues() {
    return value;
}

// Calculates and returns the gradients of the function for each of the moments along
// each of the intrinsic moment axes
MomentsModel::Jacobian MomentsModel::computeJacobian() {
    // Calculate the matrices that will be needed for the rest of this
    // function
    Moments constantGrad = Constant::computeGradient(Q, W);

    Moments alphaXGrad, alphaYGrad;
    std::tie(alphaXGrad, alphaYGrad) = makeAlphaGrad(Q, W);

    Moments betaXGrad, betaXYGrad, betaYGrad;
    std::tie(betaXGrad, betaXYGrad, betaYGrad) = makeBetaGrad(Q, W);

    // Calculate the gradient along the first moment
    Moments zerothGrad;

    // the first term in the equation is only needed for the derivative
    // with respect to Q_0, otherwise the term always is zero
    zerothGrad = Q[0]*constantGrad;
    zerothGrad[0] += constant;

    // Calculate the gradient along the two components of the first moment
    Moments firstX, firstY;

    firstX = alpha[0]*zerothGrad + value[0]*alphaXGrad;
    firstY = alpha[1]*zerothGrad + value[0]*alphaYGrad;

    // Calculate the gradient along each of the components of the second
    // moment
    SecondMoment modBeta = beta + alpha*alpha.transpose();

    Moments secondX, secondXY, secondY;

    secondX = modBeta(0, 0)*zerothGrad + value[0]*(betaXGrad + 2*alpha(0, 0)*alphaXGrad);
    secondXY = modBeta(0, 1)*zerothGrad + 
               value[0]*(betaXYGrad + 2*(alphaXGrad*alpha(1, 0) +alphaYGrad*alpha(0, 0)));
    secondY = modBeta(1, 1)*zerothGrad + value[0]*(betaYGrad + 2*alpha(1, 0)*alphaYGrad);

    // Build the result and return it
    Jacobian result;

    result.row(0) = zerothGrad.transpose();
    result.row(1) = firstX.transpose();
    result.row(2) = firstY.transpose();
    result.row(3) = secondX.transpose();
    result.row(4) = secondXY.transpose();
    result.row(5) = secondY.transpose();

    return result;
}

// anonymous namespace fore prior functions
namespace {
    auto convertShear(Quadrupole const & second) {
        lsst::afw::geom::ellipse::Separable<ConformalShear, TraceRadius> shear;

        eigen::Matrix<double, 3, 3> dConformalShear = shear.dAssign(second);

        double e = std::sqrt(std::pow(shear.getE1(), 2) + std::pow(shear.getE2(), 2));

        eigen::Matrix<double, 3, 4> dIsotropy;
        dIsotropy.zero(3, 4);
        dIsotropy(0, 0) = 1; // Carry through the flux term unchanged
        dIsotropy(1, 1) = shear.getE1()/e;
        dIsotropy(1, 2) = shear.getE2()/e;
        dIsotropy(2, 3) = 1; // Carry through the radial term unchanged

        return std::make_tuple(e, shear.getTraceRadius(), dIsotropy, dConformalShear);
    }

    auto convertBoxCox(double flux, double e, double radius, ShapePrior::Triplet lambdaVec) {
        double bcFlux = (std::pow(flux, lambdaVec[0]) - 1)/lambdaVec[0];
        double bcE = (std::pow(e, lambdaVec[1]) - 1)/lambdaVec[1];
        double bcRadius = (std::pow(radius, lambdaVec[2]) - 1)/lambdaVec[2];

        eigen::Matrix<double, 3, 3> dBoxcox;
        dBoxCox.zero(3, 3);
        dBoxCox(0, 0) << std::pow(flux, lambdaVec[0]-1);
        dBoxCox(1, 1) << std::pow(e, lambdaVec[1]-1);
        dBoxCox(2, 2) << std::pow(radius, lambdaVec[2]-1);

        return std::make_tuple(bcFlux, bcE, bcRadius, dBoxcox);
    }
} // end anonymous namespace for prior functions

void GalaxyPrior::setParameters(Quadrupole const & second, double flux) {
    // Convert to the sky coordinates
    Quadrupole::Transformer skyQuadTransform(second.transform(transformation.geometric.getLinear()));

    Quadrupole skyQuad(skyQuadTransform.copy());

    Quadrupole::Transformer::DerivativeMatrix dSky = skyQuadTransform.d();

    double radius, shape;
    eigen::Matrix<double, 3, 3> dConformalShear;
    eigen::Matrix<double, 3, 4> dIsotropy;

    // convert to shier coordinates
    std::tie(shape, radius, dIsotropy, dConformalShear) = convertShear(skyQuad);

    double bcFlux, bcRadius, bcShape; 
    eigen::Matrix<double, 3, 3> dBoxCox;

    // convert to box cox
    std::tie(bcFlux, bcRadius, bcShape, dBoxCox) = convertBoxCox(flux, shape, radius, boxCoxParams);

    convertedParameters = Tripplet(bcFlux, bcRadius, bcShape);
    logProbability = GMM.evaluate(3, convertedParameters);

    Tripplet probGrad;
    GMM.evaluateDerivatives(convertedParameters, probGrad);

    logDerivative = probGrad*dBoxCox*dIsotropy*dConformalShear*dSky;
}

double computeLogProbability() const {
   return logProbability; 
}

PriorGrad computeLogDerivative() const {
    return logDerivative;
}

void GalaxyPrior::GalaxyPrior(Quadrupole psf, const std::shared_ptr<geom::SkyWcs const> wcs,
                              const std::shared_ptr<image::calib const> calib,
                              afw::geom::Point2D const & location): wcs(wcs), calib(calib),
                              boxCoxParams({-0.29366909334755653, -1.9465686350090095, 0.2299077796501515}){

    eigen::Matrix<double, 5, 1> weights;
    eigen::Matrix<double, 5, 3> means;
    std::vector<eigen::Matrix<double, 3, 3>> covariance(5);

    weights << 0.19396571, 0.15883761, 0.20203353, 0.1760483 , 0.26911486;

    means << 1.41521881, -1.74300431, -1.03472421,
             2.18321609, -1.25418083, -1.17840354,
             0.75704646, -2.53491701, -1.24165501,
             0.70224685, -3.19682058, -1.61864423,
             1.58472117, -2.76518943, -1.72620849;

    covariance(0) << 0.21332199, -0.0006038 , -0.04732416,
                     -0.0006038 ,  0.24309651,  0.11458489,
                    -0.04732416,  0.11458489,  0.16618354;

    covariance(1) << 0.12946727,  0.07342582,  0.01472728,
                     0.07342582,  0.26472774,  0.12345566,
                     0.01472728,  0.12345566,  0.2841772;

    covariance(2) << 0.13257331,  0.06555505, -0.06795207,
                     0.06555505,  0.26217165,  0.02080632,
                     -0.06795207,  0.02080632,  0.17185127;

    covaraince(3) << 0.19290655,  0.17001021, -0.09368398,
                     0.17001021,  0.35668521, -0.06725099,
                     -0.09368398, -0.06725099,  0.16826779;

    covariance(4) << 0.20690505,  0.11141595, -0.02040214,
                     0.11141595,  0.21465827,  0.09367842,
                    -0.02040214,  0.09367842,  0.17228232;

    std::vector<Mixture::Component> ComponentList(weights.size());
    
    for (std::size_t i = 0; i < weights.size(); ++i){
        ComponentList.emplace_back(weights[i], means.row(i), covariance[i]);
    }

    GMM = std::make_unique<Mixture>(weights.size(), ComponentList);
    fluxProjection = GMM.project(Parameter::Flux);

    UnitSystem source(wcs, calib);
    UnitSystem destination(wcs->pixelToSky(location), magZero);

    localTransfrom = LocalUnitTransform(location, source, destination);

    logDerivative.zero(4);
}


}}} // Close lsst::meas::modelfit
