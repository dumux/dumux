// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the bowl problem.
 */
#ifndef DUMUX_BOWL_SPATIAL_PARAMETERS_HH
#define DUMUX_BOWL_SPATIAL_PARAMETERS_HH

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/spatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the bowl test.
 * \note Analytical solution from Thacker (1981) section 4
 *       William Thacker, "Some exact solutions to the nonlinear shallow-water wave equations", Journal
 *       of Fluid Mechanics, 107:499–508, 1981
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class BowlSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, BowlSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = BowlSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
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
        frictionLaw_ = std::make_unique<FrictionLawNoFriction<VolumeVariables>>();
        bedSurface_.assign(this->gridGeometry().gridView().size(0), 0.0);

        // compute the bedSurface for all elements
        const Scalar D0 = getParam<Scalar>("Problem.BowlDepthAtCenter");
        const Scalar L = getParam<Scalar>("Problem.BowlParaboloidRadius");
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& globalPos = element.geometry().center();
            // see Thacker (1981) Eq. (13) + shift bowl so that midpoint is a z=0
            // the shift is arbitrary and doesn't affect the analytical solution
            const auto radiusSquared = globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1];
            const auto D = D0*(1.0 - radiusSquared/(L*L));
            bedSurface_[eIdx] = D0 - D;
        }
    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
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
    { return *frictionLaw_; }

    /*! \brief Define the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return The bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    { return bedSurface_[scv.elementIndex()]; }

private:
    std::vector<Scalar> bedSurface_;
    Scalar gravity_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
