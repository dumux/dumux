// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FacetTests
 * \brief The problem for the bulk domain in the single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULKPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULKPROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief Test problem for the incompressible one-phase model
 *        with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    OnePBulkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" +
                         getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[1] < 1e-6 || globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
            values.setAllDirichlet();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! Evaluates the Neumann boundary condition for a boundary segment.
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const auto g = this->spatialParams().gravity(globalPos)[dimWorld-1];
        const auto h = 1.0 - (3.0-globalPos[dimWorld-1])*g;
        return PrimaryVariables(h);
    }

    //! Returns the temperature in \f$\mathrm{[K]}\f$ in the domain.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
