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
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"
#include "testcase.hh"

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct DarcyOneP { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOneP> { using type = Dumux::DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOneP> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOneP>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test
 */
template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr auto velocityXIdx = 0;
    static constexpr auto velocityYIdx = 1;
    static constexpr auto pressureIdx = 2;

public:
    //! export the Indices
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    const TestCase testCase,
                    const std::string& name)
    : ParentType(gridGeometry, spatialParams, "Darcy")
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
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();
        else
            values.setAllDirichlet();

        return values;
    }

        /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const auto p = analyticalSolution(scvf.center())[pressureIdx];
        return PrimaryVariables(p);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
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
     * \brief Evaluates the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Element &element) const
    {
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    auto analyticalSolution(const GlobalPosition& globalPos) const
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

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:

    // see Rybak et al., 2015: "Multirate time integration for coupled saturated/unsaturated porous medium and free flow systems"
    Dune::FieldVector<Scalar, 3> analyticalSolutionRybak_(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::exp; using std::sin; using std::cos;
        sol[velocityXIdx] = -0.5*M_PI*y*y*cos(M_PI*x);
        sol[velocityYIdx] = -1.0*y*sin(M_PI*x);
        sol[pressureIdx] = 0.5*y*y*sin(M_PI*x);
        return sol;
    }

    // see Rybak et al., 2015: "Multirate time integration for coupled saturated/unsaturated porous medium and free flow systems"
    NumEqVector rhsRybak_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::sin;
        return NumEqVector((0.5*M_PI*y*M_PI*y - 1)*sin(M_PI*x));
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    Dune::FieldVector<Scalar, 3> analyticalSolutionShiueEtAlExampleOne_(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::exp; using std::sin; using std::cos;
        sol[pressureIdx] = (exp(y) - y*exp(1)) * cos(M_PI*x);
        sol[velocityXIdx] = M_PI*(-exp(1)*y + exp(y))*sin(M_PI*x);
        sol[velocityYIdx] = (exp(1) - exp(y))*cos(M_PI*x);

        return sol;
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    NumEqVector rhsShiueEtAlExampleOne_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        return NumEqVector(cos(M_PI*x) * (M_PI*M_PI*(exp(y) - y*exp(1)) - exp(y)));
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    Dune::FieldVector<Scalar, 3> analyticalSolutionShiueEtAlExampleTwo_(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        sol[pressureIdx] = x*(1.0-x)*(y-1.0) + (y-1.0)*(y-1.0)*(y-1.0)/3.0 + 2.0*x + 2.0*y + 4.0;
        sol[velocityXIdx] = x*(y - 1.0) + (x - 1.0)*(y - 1.0) - 2.0;
        sol[velocityYIdx] = x*(x - 1.0) - 1.0*(y - 1.0)*(y - 1.0) - 2.0;

        return sol;
    }

    // see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
    NumEqVector rhsShiueEtAlExampleTwo_(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    // see Schneider et al., 2019: "Coupling staggered-grid and MPFA finite volume methods for
    // free flow/porous-medium flow problems"
    Dune::FieldVector<Scalar, 3> analyticalSolutionSchneiderEtAl_(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        static constexpr Scalar omega = M_PI;
        static constexpr Scalar c = 0.0;
        using std::exp; using std::sin; using std::cos;
        const Scalar sinOmegaX = sin(omega*x);
        const Scalar cosOmegaX = cos(omega*x);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(y+1);

        sol[pressureIdx] = (expYPlusOne + 2 - expTwo)*sinOmegaX;
        sol[velocityXIdx] = c/(2*omega)*expYPlusOne*sinOmegaX*sinOmegaX
                            -omega*(expYPlusOne + 2 - expTwo)*cosOmegaX;
        sol[velocityYIdx] = (0.5*c*(expYPlusOne + 2 - expTwo)*cosOmegaX
                            -(c*cosOmegaX + 1)*exp(y-1))*sinOmegaX;

        return sol;
    }

    // see Schneider et al., 2019: "Coupling staggered-grid and MPFA finite volume methods for
    // free flow/porous-medium flow problems (with c = 0)"
    NumEqVector rhsSchneiderEtAl_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        static constexpr Scalar omega = M_PI;
        static constexpr Scalar c = 0.0;
        const Scalar cosOmegaX = cos(omega*x);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(y+1);

        const Scalar result = (-(c*cosOmegaX + 1)*exp(y - 1)
                                             + 1.5*c*expYPlusOne*cosOmegaX
                                             + omega*omega*(expYPlusOne - expTwo + 2))
                                             *sin(omega*x);
        return NumEqVector(result);
    }

    static constexpr Scalar eps_ = 1e-7;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::string problemName_;
    TestCase testCase_;
};
} // end namespace Dumux

#endif
