// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The spatial parameters class for the soil problem
 */
#ifndef DUMUX_TEST_ROOT_SOIL_BENCHMARK_SOIL_SPATIALPARAMS_HH
#define DUMUX_TEST_ROOT_SOIL_BENCHMARK_SOIL_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The spatial parameters class for the soil problem
 */
template<class GridGeometry, class Scalar>
class SoilSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, SoilSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = SoilSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;
    // export material law type
    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

    SoilSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSw_("Soil.SpatialParams")
    {
        // perm and poro
        permeability_ = getParam<Scalar>("Soil.SpatialParams.Permeability");
        porosity_ = getParam<Scalar>("Soil.SpatialParams.Porosity");
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSw_);  }

    /*!
     * \brief Returns the temperature within the domain in [K].
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 10.0; }

private:
    const PcKrSwCurve pcKrSw_;
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
