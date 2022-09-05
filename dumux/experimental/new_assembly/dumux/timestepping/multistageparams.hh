// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TimeStepping
 * \brief Parameters for a stage in a multi-stage time integration step.
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_PARAMS_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_PARAMS_HH

#include <vector>
#include <cmath>
#include <memory>

namespace Dumux {

//! forward declaration
template<class Scalar>
class MultiStageMethod;

//! Data object for the parameters of a given stage
template<class Scalar>
class MultiStageParams
{
    struct Params {
        Scalar alpha, betaDt, timeAtStage, dtFraction;
        bool skipTemporal, skipSpatial;
    };
public:
    //! Extract params for stage i from method m
    MultiStageParams(const MultiStageMethod<Scalar>& m, std::size_t i, const Scalar t, const Scalar dt)
    : size_(i+1)
    {
        params_.resize(size_);
        for (std::size_t k = 0; k < size_; ++k)
        {
            auto& p = params_[k];
            p.alpha = m.temporalWeight(i, k);
            p.betaDt = m.spatialWeight(i, k)*dt;
            p.timeAtStage = t + m.timeStepWeight(k)*dt;
            p.dtFraction = m.timeStepWeight(k);

            using std::abs;
            p.skipTemporal = (abs(p.alpha) < 1e-6);
            p.skipSpatial = (abs(p.betaDt) < 1e-6);
        }
    }

    std::size_t size() const
    { return size_; }

    //! weights of the temporal operator residual (\f$ \alpha_{ik} \f$)
    Scalar temporalWeight(std::size_t k) const
    { return params_[k].alpha; }

    //! weights of the spatial operator residual (\f$ \beta_{ik} \f$)
    Scalar spatialWeight(std::size_t k) const
    { return params_[k].betaDt; }

    //! the time at which we have to evaluate the operators
    Scalar timeAtStage(std::size_t k) const
    { return params_[k].timeAtStage; }

    //! the fraction of a time step corresponding to the k-th stage
    Scalar timeStepFraction(std::size_t k) const
    { return params_[k].dtFraction; }

    //! If \f$ \alpha_{ik} = 0\f$
    bool skipTemporal(std::size_t k) const
    { return params_[k].skipTemporal; }

    //! If \f$ \beta_{ik} = 0\f$
    bool skipSpatial(std::size_t k) const
    { return params_[k].skipSpatial; }

private:
    std::size_t size_;
    std::vector<Params> params_;
};

} // namespace Dumux

#endif
