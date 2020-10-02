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
 * \ingroup BoundaryTests
 * \brief The free-flow sub-problem of coupled FreeFlow/Darcy convergence test
 */

#ifndef DUMUX_FREEFLOW_SUBPROBLEM_HH
#define DUMUX_FREEFLOW_SUBPROBLEM_HH

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "testcase.hh"

namespace Dumux {
template <class TypeTag>
class FreeFlowSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct FreeFlowOneP { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowOneP> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOneP> { using type = Dumux::FreeFlowSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief The Stokes sub-problem of coupled Stokes-Darcy convergence test
 */
template <class TypeTag>
class FreeFlowSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    //! export the Indices
    using Indices = typename ModelTraits::Indices;

    FreeFlowSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                       std::shared_ptr<CouplingManager> couplingManager,
                       const TestCase testCase,
                       const std::string& name)
    : ParentType(gridGeometry, "FreeFlow")
    , couplingManager_(couplingManager)
    , testCase_(testCase)
    {
        problemName_ = name + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

   /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        switch (testCase_)
        {
            case TestCase::ShiueExampleOne:
                return rhsShiueEtAlExampleOne_(globalPos);
            case TestCase::ShiueExampleTwo:
                return rhsShiueEtAlExampleTwo_(globalPos);
            case TestCase::Rybak:
                return rhsRybak_(globalPos);
            case TestCase::Schneider:
                return rhsSchneiderEtAl_(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
        }
        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
              for the Beavers-Joseph-Saffman boundary condition
     */
    auto permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element, scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        switch (testCase_)
        {
            case TestCase::ShiueExampleOne:
                return analyticalSolutionShiueEtAlExampleOne_(globalPos);
            case TestCase::ShiueExampleTwo:
                return analyticalSolutionShiueEtAlExampleTwo_(globalPos);
            case TestCase::Rybak:
                return analyticalSolutionRybak_(globalPos);
            case TestCase::Schneider:
                return analyticalSolutionSchneiderEtAl_(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

private:

    // see Rybak et al., 2015: "Multirate time integration for coupled saturated/unsaturated porous medium and free flow systems"
    PrimaryVariables analyticalSolutionRybak_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::sin; using std::cos;
        sol[Indices::velocityXIdx] = -cos(M_PI*x)*sin(M_PI*y);
        sol[Indices::velocityYIdx] = sin(M_PI*x)*cos(M_PI*y);
        sol[Indices::pressureIdx] = 0.5*y*sin(M_PI*x);
        return sol;
    }

    // see Rybak et al., 2015: "Multirate time integration for coupled saturated/unsaturated porous medium and free flow systems"
    NumEqVector rhsRybak_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        NumEqVector source(0.0);
        using std::sin; using std::cos;
        source[Indices::momentumXBalanceIdx] = 0.5*M_PI*y*cos(M_PI*x) - 2*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y);
        source[Indices::momentumYBalanceIdx] = 0.5*sin(M_PI*x) + 2*M_PI*M_PI*sin(M_PI*x)*cos(M_PI*y);
        return source;
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    PrimaryVariables analyticalSolutionShiueEtAlExampleOne_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::exp; using std::sin; using std::cos;
        sol[Indices::velocityXIdx] = -1/M_PI * exp(y) * sin(M_PI*x);
        sol[Indices::velocityYIdx] = (exp(y) - exp(1)) * cos(M_PI*x);
        sol[Indices::pressureIdx] = 2*exp(y) * cos(M_PI*x);
        return sol;
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    NumEqVector rhsShiueEtAlExampleOne_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        NumEqVector source(0.0);
        source[Indices::momentumXBalanceIdx] = exp(y)*sin(M_PI*x) * (1/M_PI -3*M_PI);
        source[Indices::momentumYBalanceIdx] = cos(M_PI*x) * (M_PI*M_PI*(exp(y)- exp(1)) + exp(y));
        return source;
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    PrimaryVariables analyticalSolutionShiueEtAlExampleTwo_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        sol[Indices::velocityXIdx] = (y-1.0)*(y-1.0) + x*(y-1.0) + 3.0*x - 1.0;
        sol[Indices::velocityYIdx] = x*(x-1.0) - 0.5*(y-1.0)*(y-1.0) - 3.0*y + 1.0;
        sol[Indices::pressureIdx] = 2.0*x + y - 1.0;
        return sol;
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    NumEqVector rhsShiueEtAlExampleTwo_(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    // see Schneider et al., 2019: "Coupling staggered-grid and MPFA finite volume methods for
    // free flow/porous-medium flow problems"
    PrimaryVariables analyticalSolutionSchneiderEtAl_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::sin;
        static constexpr Scalar omega = M_PI;
        const Scalar sinOmegaX = sin(omega*x);

        sol[Indices::velocityXIdx] = y;
        sol[Indices::velocityYIdx] = -y*sinOmegaX;
        sol[Indices::pressureIdx] = -y*y*sinOmegaX*sinOmegaX;
        return sol;
    }

    // see Schneider et al., 2019: "Coupling staggered-grid and MPFA finite volume methods for
    // free flow/porous-medium flow problems"
    NumEqVector rhsSchneiderEtAl_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        static constexpr Scalar omega = M_PI;
        const Scalar sinOmegaX = sin(omega*x);
        const Scalar cosOmegaX = cos(omega*x);

        NumEqVector source(0.0);
        source[Indices::conti0EqIdx] = -sinOmegaX;
        source[Indices::momentumXBalanceIdx] = -2*omega*y*y*sinOmegaX*cosOmegaX
                                               -2*y*sinOmegaX + omega*cosOmegaX;
        source[Indices::momentumYBalanceIdx] = -omega*y*y*cosOmegaX - omega*omega*y*sinOmegaX;
        return source;
    }

    static constexpr Scalar eps_ = 1e-7;
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManager_;
    TestCase testCase_;
};
} // end namespace Dumux

#endif
