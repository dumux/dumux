// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem.
 */

#ifndef DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem.
 */
template<class GridGeometry, class Scalar>
class TracerTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         TracerTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                                       TracerTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        porosity_ = getParam<Scalar>("Problem.Porosity");
        velocity_ = getParam<Scalar>("Problem.FrontVelocity")*porosity_;
    }

    /*!
     * \brief if we are in the water
     */
    bool isWater(const GlobalPosition& globalPos) const
    {
        const auto minCoord = this->gridGeometry().bBoxMin()[dimWorld-1];
        const auto coord = globalPos[dimWorld-1]-minCoord;
        const auto delta = this->gridGeometry().bBoxMax()[dimWorld-1]-minCoord;
        if (coord < 0.2*delta + eps_ || coord > 0.8*delta + eps_ || (coord > 0.4*delta - eps_ && coord < 0.6*delta + eps_))
            return false;
        else
            return true;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    //! Fluid properties that are spatial parameters in the tracer model
    //! They can possibly vary with space but are usually constants
    //! saturation field
    Scalar saturation(const Element &element,
                      const SubControlVolume& scv) const
    { return isWater(scv.center()) ? 1.0 : 0.0; }

    //! Fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    //! Fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 0.018; }

    //! Velocity field
    GlobalPosition velocity(const SubControlVolumeFace& scvf) const
    {
        GlobalPosition vel(0.0);
        vel[0] = isWater(scvf.center()) ? velocity_ : 0.0;
        return vel;
    }

    //! Velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return velocity(scvf) * scvf.unitOuterNormal() * scvf.area()
               * elemVolVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
    }

private:
    Scalar velocity_, porosity_;
    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
