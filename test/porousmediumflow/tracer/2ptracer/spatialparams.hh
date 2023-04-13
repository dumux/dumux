// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
#ifndef DUMUX_TWOP_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TWOP_TRACER_TEST_SPATIAL_PARAMS_HH

#include <iostream>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
template<class GridGeometry, class Scalar>
class TwoPTracerTestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           TwoPTracerTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                         TwoPTracerTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TwoPTracerTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    //! Fluid properties that are spatial parameters in the tracer model
    //! They can possible vary with space but are usually constants

    //! Fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return density_[this->gridGeometry().elementMapper().index(element)];; }

    void setDensity(const std::vector<Scalar>& d)
    { density_ = d; }

    //! fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 0.018; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 0.018; }

    //! Velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    { return volumeFlux_[scvf.index()]; }

    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

    //! saturation from twoPProblem
    Scalar saturation(const Element &element,
                      const SubControlVolume& scv) const
    { return saturation_[scv.dofIndex()]; }

    void setSaturation(const std::vector<Scalar>& s)
    { saturation_ = s; }

private:
    std::vector<Scalar> volumeFlux_;
    std::vector<Scalar> density_;
    std::vector<Scalar> saturation_;
};

} // end namespace Dumux

#endif
