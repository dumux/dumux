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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief The problem for the bulk domain in the single-phase facet coupling test
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULKPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULKPROBLEM_HH

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePBulkProblem;

namespace Properties {
// create the type tag nodes
NEW_TYPE_TAG(OnePBulk, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePBulkTpfa, INHERITS_FROM(OnePBulk, CCTpfaFacetCouplingModel));

// Set the grid type
SET_TYPE_PROP(OnePBulk, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
// Set the problem type
SET_TYPE_PROP(OnePBulk, Problem, OnePBulkProblem<TypeTag>);
// set the spatial params
SET_TYPE_PROP(OnePBulk, SpatialParams, OnePSpatialParams<TypeTag>);

// the fluid system
SET_PROP(OnePBulk, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

// Disable caching (for testing purposes)
SET_BOOL_PROP(OnePBulk, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(OnePBulk, EnableGridFluxVariablesCache, false);
SET_BOOL_PROP(OnePBulk, EnableFVGridGeometryCache, false);

} // end namespace Properties
/*!
 * \ingroup OnePTests
 * \brief Test problem for the incompressible one-phase model
 *        with coupling across the bulk grid facets
 */
template<class TypeTag>
class OnePBulkProblem : public PorousMediumFlowProblem<TypeTag>
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

public:
    OnePBulkProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , useAnalyticProblem_(getParam<bool>("Problem.UseAnalyticProblem"))
    , useInteriorDirichletBCs_(getParam<bool>("Problem.UseInteriorDirichletBCs"))
    , lowDimPermeability_(getParam<Scalar>("LowDim.SpatialParams.Permeability"))
    , aperture_(getParam<Scalar>("Problem.LowDimAperture"))
    {}

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();
        if (!useAnalyticProblem_)
        {
            const auto xMax = this->fvGridGeometry().bBoxMax()[0];
            const auto xMin = this->fvGridGeometry().bBoxMin()[0];
            const auto x = globalPos[0];
            if (x < xMin + 1e-6 || x > xMax - 1e-6)
                values.setAllNeumann();
        }
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        if (!useInteriorDirichletBCs_) values.setAllNeumann();
        else values.setAllDirichlet();
        return values;
    }

    //! Evaluate the source term at a given position
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        if (useAnalyticProblem_)
        {
            using std::cos;
            using std::cosh;
            Scalar u = (1.0 - lowDimPermeability_)*cos(globalPos[0])*cosh(aperture_/2);
            return NumEqVector(u);
        }

        return NumEqVector(0.0);
    }

    //! evaluates the exact solution at a given position.
    Scalar exact(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::cosh;
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return lowDimPermeability_*cos(x)*cosh(y) + (1.0 - lowDimPermeability_)*cos(x)*cosh(aperture_/2);
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        if (useAnalyticProblem_)
            return PrimaryVariables(exact(globalPos));
        else
        {
            auto values = initialAtPos(globalPos);
            if (globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + 1e-6)
                values[0] *= 2.0;
            return values;
        }
    }

    //! evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); }

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
    bool useAnalyticProblem_;
    bool useInteriorDirichletBCs_;
    Scalar lowDimPermeability_;
    Scalar aperture_;
};

} // end namespace Dumux

#endif
