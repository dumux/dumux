// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of the spatial parameters for the 1pnc problems.
 */

#ifndef DUMUX_1PNC_TEST_SPATIAL_PARAMS_HH
#define DUMUX_1PNC_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/porousmediumflow/fvspatialparamsnonequilibrium.hh>
#include <dumux/material/fluidmatrixinteractions/1pia/fluidsolidinterfacialareashiwang.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Definition of the spatial parameters for the 1pnc test problems.
 */
template<class GridGeometry, class Scalar>
class OnePNCNonequilibriumTestSpatialParams
: public FVPorousMediumFlowSpatialParamsNonEquilibrium<GridGeometry, Scalar,
                                       OnePNCNonequilibriumTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsNonEquilibrium<GridGeometry, Scalar,
                                           OnePNCNonequilibriumTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;
    using FluidSolidInterfacialAreaFormulation = FluidSolidInterfacialAreaShiWang<Scalar>;

    OnePNCNonequilibriumTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeability_ = 1e-11;
        porosity_ = 0.4;

        characteristicLength_ = 5e-4;
        factorEnergyTransfer_ = 0.5;
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Returns the characteristic length for the mass transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar characteristicLengthAtPos(const  GlobalPosition & globalPos) const
    { return characteristicLength_ ; }

    /*!
     * \brief Return the pre factor the the energy transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    { return factorEnergyTransfer_; }

private:
    Scalar permeability_;
    Scalar porosity_;

    Scalar factorEnergyTransfer_ ;
    Scalar characteristicLength_ ;
};

} // end namespace Dumux

#endif
