// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCMinTests
 * \brief Definition of the spatial parameters for the thermochemistry
 *        problem which uses the non-insothermal 1pncmin model.
 */

#ifndef DUMUX_THERMOCHEM_SPATIAL_PARAMS_HH
#define DUMUX_THERMOCHEM_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCMinTests
 * \brief Definition of the spatial parameters for the thermochemistry
 *        problem which uses the non-insothermal 1pncmin model.
 */
template<class GridGeometry, class Scalar>
class ThermoChemSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         ThermoChemSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ThisType = ThermoChemSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    enum { dimWorld=GridView::dimensionworld };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;

    ThermoChemSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Returns the intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     *
     *  Solution dependent permeability function.
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    { return 8.53e-12; }

private:

   Scalar eps_;
};

} // end namespace Dumux

#endif
