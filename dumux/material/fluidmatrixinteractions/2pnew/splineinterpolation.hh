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
 * \ingroup Fluidmatrixinteractions
 * \brief Replaces the actual curve by a spline in the common interval
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_SPLINE_INTERPOLATION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_SPLINE_INTERPOLATION_HH

#include <vector>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/cubicspline.hh>

#include <dumux/material/fluidmatrixinteractions/2pnew/regularization.hh>

namespace Dumux {
namespace FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The default regularization policy
 */
template<class Scalar>
class TwoPSplineInterpolation : public TwoPRegularization<Scalar>
{
    using ParentType = TwoPRegularization<Scalar>;
public:
    using Return = typename ParentType::Return;

    //! The parameters
    template<class S>
    struct Params
    {
        std::size_t numPoints;
        Scalar minSwThreshold, maxSwThreshold;
    };

    //! Initialize the spline
    template<class MaterialLaw>
    void init(const MaterialLaw* m, const std::string& paramGroup)
    {
        const auto numPoints = getParam<std::size_t>(paramGroup + ".Spline.NumPoints", 10);
        minSwThreshold_ = getParam<Scalar>(paramGroup + ".Spline.MinSw", 0.1);
        maxSwThreshold_ = getParam<Scalar>(paramGroup + ".Spline.MaxSw", 1.0);

        initSpline_(m, numPoints);
    }

    template<class MaterialLaw, class BaseParams, class EffToAbsParams>
    void init(const MaterialLaw* m, const BaseParams& bp, const EffToAbsParams& etap, const Params<Scalar>& p)
    {
        const auto numPoints = p.numPoints;
        minSwThreshold_ = p.minSwThreshold;
        maxSwThreshold_ = p.maxSwThreshold;

        initSpline_(m, numPoints);
    }

    /*!
     * \brief The regularized capillary pressure-saturation curve
     */
    Return pc(const Scalar sw) const
    {
        if (sw < minSwThreshold_) return {};
        if (sw > maxSwThreshold_) return {};

        return splinePcSw_.eval(sw);
    }

private:
    template<class MaterialLaw>
    void initSpline_(const MaterialLaw* m, const std::size_t numPoints)
    {
        const double delta = (maxSwThreshold_-minSwThreshold_)/static_cast<double>(numPoints-1);
        std::vector<double> sw(numPoints);
        for (int i = 0; i < sw.size(); ++i)
            sw[i] = minSwThreshold_ + i*delta;

        auto pc = sw;
        std::transform(sw.begin(), sw.end(), pc.begin(), [&](auto sw) { return m->template pc<false>(sw); });

        splinePcSw_.updatePoints(sw, pc);

        std::cout << "Created spline interpolation with " << numPoints << " control points." << std::endl;
    }

    Scalar minSwThreshold_, maxSwThreshold_;
    CubicSpline<Scalar> splinePcSw_;
};

} // end namespace FluidMatrix
} // end namespace Dumux

#endif
