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
 * \brief Specification of the material parameters
 *       for the spline law constitutive relations.
 */
#ifndef DUMUX_SPLINE_LAW_PARAMS_HH
#define DUMUX_SPLINE_LAW_PARAMS_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/monotonecubicspline.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *       for the spline law constitutive relations.
 * \see SplineLaw
 */
template <class ScalarT>
class SplineLawParams
{
public:
    using Scalar = ScalarT;

    SplineLawParams()
    {
        auto swData = getParam<std::vector<double>>("SplineLaw.SwData");
        auto pcData = getParam<std::vector<double>>("SplineLaw.PcData");
        auto krwData = getParam<std::vector<double>>("SplineLaw.KrwData");
        auto krnData = getParam<std::vector<double>>("SplineLaw.KrnData");

        update(swData, pcData, krwData, krnData);
    }

    SplineLawParams(const std::vector<double>& swData,
                    const std::vector<double>& pcData,
                    const std::vector<double>& krwData,
                    const std::vector<double>& krnData)
    {
        update(swData, pcData, krwData, krnData);
    }

    void update(const std::vector<double>& swData,
                const std::vector<double>& pcData,
                const std::vector<double>& krwData,
                const std::vector<double>& krnData)
    {
        pcSpline_.updatePoints(swData, pcData);
        krwSpline_.updatePoints(swData, krwData);
        krnSpline_.updatePoints(swData, krnData);
    }

    /*!
     * \brief Equality comparison with another set of params
     */
    template<class OtherParams>
    bool operator== (const OtherParams& otherParams) const
    {
        DUNE_THROW(Dune::NotImplemented, "comparison of spline law parameters");
    }

    const MonotoneCubicSpline<Scalar>& pcSpline() const
    { return pcSpline_; }

    const MonotoneCubicSpline<Scalar>& krwSpline() const
    { return krwSpline_; }

    const MonotoneCubicSpline<Scalar>& krnSpline() const
    { return krnSpline_; }

private:
    MonotoneCubicSpline<Scalar> pcSpline_;
    MonotoneCubicSpline<Scalar> krwSpline_;
    MonotoneCubicSpline<Scalar> krnSpline_;
};
} // namespace Dumux

#endif
