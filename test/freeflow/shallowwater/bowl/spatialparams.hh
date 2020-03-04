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
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the bowl problem.
 */
#ifndef DUMUX_BOWL_SPATIAL_PARAMETERS_HH
#define DUMUX_BOWL_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the bowl test.
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class BowlSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         BowlSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = BowlSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    BowlSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        const Scalar bowlDepthAtCenter =  getParam<Scalar>("Problem.BowlDepthAtCenter");
        const Scalar bowlParaboloidRadius =  getParam<Scalar>("Problem.BowlParaboloidRadius");
        initFrictionLaw();
        bedSurface_.assign(this->gridGeometry().gridView().size(0), 0.0);

        //compute the bedSurface for all elements
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& globalPos = element.geometry().center();
            bedSurface_[eIdx] = bowlDepthAtCenter *(((globalPos[0] * globalPos[0]) +
                                (globalPos[1] * globalPos[1])) /
                                (bowlParaboloidRadius * bowlParaboloidRadius) - 1.0) + 1.0;
        }

    }

    /*!
     * \brief Initialize the FrictionLaw
     */
    void initFrictionLaw()
    {
        frictionLaw_ = std::make_unique<FrictionLawNoFriction<VolumeVariables>>();
    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity() const
    {
        return gravity_;
    }

    /*! \brief Get the frictionLaw.
    *
    * Get the frictionLaw, which already includes the friction value.
    *
    * \return frictionLaw
    */

    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    {
        return *frictionLaw_;
    }

    /*! \brief Define the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return The bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return bedSurface_[eIdx];
    }

private:
    std::vector<Scalar> bedSurface_;
    Scalar gravity_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
