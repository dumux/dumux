// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the rough channel problem with limited
 *        Nikuradse friction
 */
#ifndef DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH
#define DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>
#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the test with limited Nikuradse
 * friction.
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class RoughChannelSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar,
                               RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    RoughChannelSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        initFrictionLaw_();
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
        // todo depends on index e.g. eIdx = scv.elementIndex();
        return 10.0 - element.geometry().center()[0] * bedSlope_;
    }

private:
    void initFrictionLaw_()
    {
        // Get the roughnessHeight which will be used to limit the friction law
        const Scalar roughnessHeight = getParam<Scalar>("Problem.RoughnessHeight", 0.0);
        const Scalar ks = getParam<Scalar>("Problem.Ks"); // equivalent sand roughness
        frictionLaw_ = std::make_unique<FrictionLawNikuradse<VolumeVariables>>(ks, roughnessHeight);
    }

    Scalar gravity_;
    Scalar bedSlope_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
