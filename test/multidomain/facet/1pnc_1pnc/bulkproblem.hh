// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief The problem for the bulk domain in the single-phase two-components facet coupling test
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_BULKPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_BULKPROBLEM_HH

#include <dune/alugrid/grid.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePNCBulkProblem;

namespace Properties {
// create the type tag nodes
NEW_TYPE_TAG(OnePNCBulk, INHERITS_FROM(OnePNC));
NEW_TYPE_TAG(OnePNCBulkTpfa, INHERITS_FROM(OnePNCBulk, CCTpfaFacetCouplingModel));
NEW_TYPE_TAG(OnePNCBulkBox, INHERITS_FROM(OnePNCBulk, BoxFacetCouplingModel));

// Set the grid type
SET_TYPE_PROP(OnePNCBulk, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
// Set the problem type
SET_TYPE_PROP(OnePNCBulk, Problem, OnePNCBulkProblem<TypeTag>);
// set the spatial params
SET_TYPE_PROP(OnePNCBulk,
              SpatialParams,
              OnePSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), typename GET_PROP_TYPE(TypeTag, Scalar)>);
// the fluid system
SET_TYPE_PROP(OnePNCBulk,
              FluidSystem,
              FluidSystems::OnePAdapter< FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), FluidSystems::H2ON2DefaultPolicy<true>> >);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePNCBulk, UseMoles, true);
} // end namespace Properties
/*!
 * \ingroup OnePTests
 * \brief Test problem for the incompressible one-phase model
 *        with coupling across the bulk grid facets
 */
template<class TypeTag>
class OnePNCBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    // copy some indices for convenience
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),
    };

public:
    OnePNCBulkProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                      std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                      const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , boundaryMoleFraction_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundaryMoleFraction"))
    {}

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if ( (globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + 1e-6 && globalPos[1] < 0.0)
             || (globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - 1e-6 && globalPos[1] > 0.0) )
            values.setAllDirichlet();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        if (globalPos[0] < 1e-6)
            values[N2Idx] = boundaryMoleFraction_;
        return values;
    }

    //! evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables({1.0e5, 0.0}); }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar boundaryMoleFraction_;
};

} // end namespace Dumux

#endif
