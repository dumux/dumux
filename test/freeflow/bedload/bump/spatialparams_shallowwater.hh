// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \copydoc Dumux::BumpTestSpatialParamsShallowWater
 */
#ifndef DUMUX_BUMP_SPATIAL_PARAMETERS_SHALLOWWATER_HH
#define DUMUX_BUMP_SPATIAL_PARAMETERS_SHALLOWWATER_HH

#include <dumux/freeflow/spatialparams.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>

namespace Dumux {

/*!
 * \ingroup BedloadTests
 * \brief The spatial parameter class for the bump test (shallow water part)
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class BumpTestSpatialParamsShallowWater
: public FreeFlowSpatialParams<GridGeometry, Scalar, BumpTestSpatialParamsShallowWater<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = BumpTestSpatialParamsShallowWater<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    BumpTestSpatialParamsShallowWater(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity", 9.81);
        frictionLaw_ = std::make_unique<FrictionLawNoFriction<VolumeVariables>>();

        int nElems = this->gridGeometry().gridView().size(0);
        // bedSurface_ must be filled to use updateCoupledVariables()
        bedSurface_.assign(nElems, 0.0);
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
    * \return gravity constant
    */
    Scalar gravity() const
    {
        return gravity_;
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
    Scalar waterDensity_ = 1000.0;  // [kg/m^3]
};

} // end namespace Dumux

#endif
