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

#include <dumux/porousmediumflow/properties.hh>
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
        alphaL_ = getParam<Scalar>("Problem.AlphaL", 0.0);
        alphaT_ = getParam<Scalar>("Problem.AlphaT", 0.0);
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    //! Fluid properties that are spatial parameters in the tracer model
    //! They can possibly vary with space but are usually constants

    //! Fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    //! Fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 18.0; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 18.0; }

    //! Velocity field
    GlobalPosition velocity(const SubControlVolumeFace& scvf) const
    {
        GlobalPosition vel(1e-5);
        const auto globalPos = scvf.ipGlobal();
        const auto& x = globalPos[0];
        const auto& y = globalPos[1];

        vel[0] *= x*x * (1.0 - x)*(1.0 - x) * (2.0*y - 6.0*y*y + 4.0*y*y*y);
        vel[1] *= -1.0*y*y * (1.0 - y)*(1.0 - y) * (2.0*x - 6.0*x*x + 4.0*x*x*x);

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

    /*!
     * \brief Defines the dispersion tensor \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    std::array<Scalar, 2> dispersionAlphas(const GlobalPosition& globalPos,
                                           const int phaseIdx = 0,
                                           const int compIdx = 0) const
    { return { alphaL_, alphaT_ }; }

private:
    Scalar alphaL_;
    Scalar alphaT_;
};

} // end namespace Dumux

#endif
