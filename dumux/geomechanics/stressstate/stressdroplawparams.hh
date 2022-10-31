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
 * \ingroup Geomechanics
 * \brief Stress drop law checks if the shear failure occurs and reduces
 *  the shear stress directly from the failure stress state.
 *  The failure stress state maybe rotated.
 */
#ifndef DUMUX_STRESS_STATE_STRESS_DROP_LAW_PARAMS_HH
#define DUMUX_STRESS_STATE_STRESS_DROP_LAW_PARAMS_HH

#include <cmath>
#include <algorithm>
#include <dune/common/exceptions.hh>
namespace Dumux {

/*!
 * \brief Stress Drop Law Params
 * per default, the slope is negative and
 * cohesion is negative.
 * \tparam Scalar
 */
template <class Scalar, class Line>
class StressDropLawParams{
public:
    /*!
     * \brief Construct a new Stress Drop Law Params object
     *
     * \param internalFrictionAngle angle in degree
     * \param cohesion
     */
    StressDropLawParams(const Scalar& internalFrictionAngle,
                        const Scalar& cohesion,
                        const Scalar& stressDrop)
    : internalFrictionAngle_{internalFrictionAngle}
    , cohesion_{cohesion}
    , stressDrop_(stressDrop)
    {
        if(internalFrictionAngle != std::clamp(internalFrictionAngle, -90.0, 0.0)){
            DUNE_THROW(Dune::RangeError,"internal friction angle should be between -90 and 0 degree.");
        }

        if(cohesion < 0.0) {
            DUNE_THROW(Dune::RangeError,"cohesion should not be smaller than 0.");
        }

        if(stressDrop < 0.0) {
            DUNE_THROW(Dune::RangeError,"stress drop value should be positiv.");
        }

        using std::tan;
        envelopeCurveSlope_ = tan(internalFrictionAngle * M_PI / 180);
    }

    Scalar cohesion() const
    { return cohesion_;}

    Scalar envelopeCurveSlope() const
    { return envelopeCurveSlope_; }

    Scalar stressDrop() const
    { return stressDrop_; }

    Line envelopeLine(const bool& inUpperHalfPlane = true) const
    {
        const Scalar sign = inUpperHalfPlane ? 1 : -1;
        return Line(sign*envelopeCurveSlope_, sign*cohesion_);
    }
private:
    Scalar internalFrictionAngle_;
    Scalar cohesion_;
    Scalar envelopeCurveSlope_;
    Scalar stressDrop_;
}; // end class StressDropLawParams

} // end namespace dumux
#endif
