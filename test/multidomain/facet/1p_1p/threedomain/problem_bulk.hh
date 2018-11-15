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
 * \brief The problem for the bulk domain in the single-phase facet coupling test
 */
#ifndef DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_BULKPROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_BULKPROBLEM_HH

#include <dune/alugrid/grid.hh>

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
// Create new type tags
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
SET_TYPE_PROP(OnePBulk, Grid, Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>);
// Set the problem type
SET_TYPE_PROP(OnePBulk, Problem, OnePBulkProblem<TypeTag>);
// set the spatial params
SET_PROP(OnePBulk, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = OnePSpatialParams<FVGridGeometry, Scalar>;
};

// the fluid system
SET_PROP(OnePBulk, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

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

    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    //! The constructor
    OnePBulkProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManagerPtr,
                    const std::string& paramGroup = "Bulk")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Specifies the kind of boundary condition on a given boundary position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + 1e-6 || globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - 1e-6)
            values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Specifies which kind of interior boundary condition should be
     *        used for which equation on a given sub-control volume face
     *        that couples to a facet element.
     *
     * \param element The finite element the scvf is embedded in
     * \param scvf The sub-control volume face
     */
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluate the Dirichlet boundary conditions at a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        const auto y = globalPos[1];
        const auto yMin = this->fvGridGeometry().bBoxMin()[1];
        const auto yMax = this->fvGridGeometry().bBoxMax()[1];

        return PrimaryVariables( {2.0 - (y-yMin)/(yMax-yMin)} );
    }

    //! Evaluat the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); }

    //! Return the temperature in \f$\mathrm{[K]}\f$ in the domain
    Scalar temperature() const
    { return 283.15; /*10°C*/ }

    //! Return const reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
};

} // end namespace Dumux

#endif
