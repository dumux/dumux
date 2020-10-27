// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup FreeflowModels
 * \brief The available free flow slip conditions in Dumux
 */
#ifndef DUMUX_FREEFLOW_SLIPCONDITIONS_HH
#define DUMUX_FREEFLOW_SLIPCONDITIONS_HH

#include <string>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
    * \brief The available free flow slip conditions in Dumux
    * \ingroup FreeflowModels
    */
enum class SlipCondition
{
    BJS, BJ, ER
};

/**
 * \brief return the name of the Turbulence Model
 */
std::string slipConditionToString(SlipCondition condition)
{
    switch (condition)
    {
        case SlipCondition::BJS: return "Beavers-Joseph-Saffman";
        case SlipCondition::BJ:  return "Beavers-Joseph";
        case SlipCondition::ER:  return "Eggenweiler-Rybak";
        default:                 return "Beavers-Joseph-Saffman";
    }
}

/**
 * \brief return the tpye of slip condition
 */
SlipCondition slipCondition()
{
    static const std::string condition = getParam<std::string>("Problem.SlipCondition", "BJS");
    if (condition == "BJS")
        return SlipCondition::BJS;
    else if (condition == "BJ")
        return SlipCondition::BJ;
    else if (condition == "ER")
        return SlipCondition::ER;
    else
        DUNE_THROW(Dune::InvalidStateException, condition + " is not a valid slip condition");
}

template<class GridGeometry>
class SlipVelocity
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
      };

public:
    template<class Problem, class Scalar>
    static Scalar velocity(const Problem& problem,
                           const Element& element,
                           const SubControlVolume& scv,
                           const SubControlVolumeFace& ownScvf,
                           const SubControlVolumeFace& faceOnPorousBoundary,
                           const Scalar velocitySelf, //vel perpendicular to tangential vel
                           const Scalar tangentialVelocityGradient) //dv/dx (=0)
    {
        using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
        static const SlipCondition condition_ = slipCondition();

        // create a unit normal vector oriented in positive coordinate direction
        auto orientation = ownScvf.unitOuterNormal();
        orientation[ownScvf.directionIndex()] = 1.0;

        Scalar factor(0.0);
        if (condition_ == SlipCondition::ER)
            factor = -1.0 / problem.epsInterface(faceOnPorousBoundary) / problem.factorNTangential(faceOnPorousBoundary);
        else
            factor = problem.betaBJ(element, faceOnPorousBoundary, orientation); //beta = alpha/sqrt(K)

        const Scalar distanceNormalToBoundary = (faceOnPorousBoundary.center() - scv.center()).two_norm();
        const auto porousMediumTerm = (condition_ == SlipCondition::BJS) ? VelocityVector(0.0)
                                                                         : problem.porousMediumTerm(element, faceOnPorousBoundary);
        return (tangentialVelocityGradient*distanceNormalToBoundary
                + porousMediumTerm * orientation * factor * distanceNormalToBoundary
                + velocitySelf) / (factor*distanceNormalToBoundary + 1.0);
    }
};

} // end namespace Dumux

#endif
