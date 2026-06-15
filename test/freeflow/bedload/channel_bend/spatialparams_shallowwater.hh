// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \copydoc Dumux::ChannelBendTestSpatialParamsShallowWater
 */
#ifndef DUMUX_CHANNEL_BEND_SPATIAL_PARAMETERS_SHALLOWWATER_HH
#define DUMUX_CHANNEL_BEND_SPATIAL_PARAMETERS_SHALLOWWATER_HH

#include <dumux/freeflow/spatialparams.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>

namespace Dumux {

/*!
 * \ingroup BedloadTests
 * \brief The spatial parameter class for the channel bend test (shallow water part)
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class ChannelBendTestSpatialParamsShallowWater
: public FreeFlowSpatialParams<GridGeometry, Scalar, ChannelBendTestSpatialParamsShallowWater<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = ChannelBendTestSpatialParamsShallowWater<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    ChannelBendTestSpatialParamsShallowWater(std::shared_ptr<const GridGeometry> gridGeometry,
                                             std::vector<Scalar>& initialBedSurface)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity", 9.81);
        // transform Manning friction coefficient to Nikuradse equivalent sand roughness and use it as roughness height
        Scalar equivalentSandRoughness = pow(25.68/(1.0/manningFrictionCoefficient_), 6);
        frictionLaw_ = std::make_unique<FrictionLawManning<VolumeVariables>>(gravity_,
                                                                             manningFrictionCoefficient_,
                                                                             equivalentSandRoughness);
        bedSurface_ = initialBedSurface;
    }

    /*! \brief Get the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return the bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return bedSurface_[eIdx];
    }

    /*! \brief Get the friction law.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return friction law
    */
    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    {
        return *frictionLaw_;
    }

    /*! \brief Get the gravitation.
    *
    * \param globalPos The global position
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Get the water density.
    *
    * \return water density
    */
    Scalar waterDensity() const
    {
        return waterDensity_;
    }

    /*! \brief Update the coupled variables.
    *
    * Update the bed surface.
    *
    * \param bedSurface The bed surface
    */
    void updateCoupledVariables(std::vector<Scalar> bedSurface)
    {
        bedSurface_ = bedSurface;
    }
private:
    std::vector<Scalar> bedSurface_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
    Scalar gravity_;
    Scalar manningFrictionCoefficient_ = 0.0167;  // The same value like in the bedload part [s/m^(1/3)]
    Scalar waterDensity_ = 1000.0;  // [kg/m^3]
};

} // end namespace Dumux

#endif
