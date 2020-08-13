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
 * \brief The Stokes sub-problem of coupled Stokes-Darcy BJ/NewIC test
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

namespace Dumux {
template <class TypeTag>
class StokesSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct StokesOneP { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief The Stokes sub-problem of coupled Stokes-Darcy convergence test
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
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

    using VelocityVector = typename Element::Geometry::GlobalCoordinate;


    enum class TestCase
    {
        BJSymmetrized, NewICNonSymmetrized
    };

public:
    //! export the Indices
    using Indices = typename ModelTraits::Indices;

    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Stokes"), couplingManager_(couplingManager)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        const auto testCaseInput = getParamFromGroup<std::string>(this->paramGroup(), "Problem.TestCase", "BJSymmetrized");
        if (testCaseInput == "BJSymmetrized")
            testCase_ = TestCase::BJSymmetrized;
        else if (testCaseInput == "NewICNonSymmetrized")
            testCase_ = TestCase::NewICNonSymmetrized;
        else
            DUNE_THROW(Dune::InvalidStateException, testCaseInput + " is not a valid test case");
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
            case TestCase::BJSymmetrized:
                return rhsBJSymmetrized_(globalPos);
            case TestCase::NewICNonSymmetrized:
                return rhsNewICNonSymmetrized_(globalPos);
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

        const auto& globalPos = scvf.dofPosition();

        if(onUpperBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        // left/right wall
        if (onRightBoundary_(globalPos) || (onLeftBoundary_(globalPos)))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        if (couplingManager().isCoupledEntity(CouplingManager::freeFlowIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setSlipCondition(Indices::momentumXBalanceIdx);
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

        if(couplingManager().isCoupledEntity(CouplingManager::freeFlowIdx, scvf))
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
    * \brief Returns the velocity in the porous medium (which is 0 by default according to Saffmann).
    */
    VelocityVector porousMediumVelocity(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().porousMediumVelocity(element, scvf);
    }



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
              for the Beavers-Joseph(-Saffman) boundary condition
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element, scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::porousMediumIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }


    //for new IC, input parameters that stay in the tangential coupling condition
    /*!
     * \brief Returns the scale separation parameter epsilon required as input parameter for the
              new coupling conditions
     */
    Scalar epsInterface(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::porousMediumIdx).spatialParams().epsInterfaceAtPos(scvf.center());
    }

    /*!
     * \brief Returns the boundary layer constant N_1_bl required as input parameter for the
              new coupling condition for the tangential component
     */
    Scalar factorNTangential(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::porousMediumIdx).spatialParams().factorNTangentialAtPos(scvf.center());
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
            case TestCase::BJSymmetrized:
                return analyticalSolutionBJSymmetrized_(globalPos);
            case TestCase::NewICNonSymmetrized:
                return analyticalSolutionNewICNonSymmetrized_(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

private:

    // see exact solution for BJ-IC with symmetrized stress tensor (by Elissa Eggenweiler)
    PrimaryVariables analyticalSolutionBJSymmetrized_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::sin; using std::cos;
        sol[Indices::velocityXIdx] = (y-1.0)*(y-1.0) + x*(y-1.0) - x - 7.0 + 2.0*x*x + 2.0*y*y;
        sol[Indices::velocityYIdx] = (x-1.0)*x - 0.5*(y-1.0)*(y-1.0) - 3.0*y - 3.0 -4.0*y*(x-1.0);
        sol[Indices::pressureIdx] = -8.0*x + 1.0*y + 7.0 + x*x;
        return sol;
    }

    // see exact solution for BJ-IC with symmetrized stress tensor (by Elissa Eggenweiler)
    NumEqVector rhsBJSymmetrized_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        NumEqVector source(0.0);
        using std::sin; using std::cos;
        source[Indices::momentumXBalanceIdx] = 2.0*x - 18.0;
        source[Indices::momentumYBalanceIdx] = 0.0;
        return source;
    }

    // see exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    PrimaryVariables analyticalSolutionNewICNonSymmetrized_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::exp; using std::sin; using std::cos;
        sol[Indices::velocityXIdx] = (y-1.0)*(y-1.0) + x*(y-1.0) + 3.0*x + 1.5;
        sol[Indices::velocityYIdx] = 0.5*x*(x-1.0) - 0.5*(y-1.0)*(y-1.0) - 3.0*y + 2.0;
        sol[Indices::pressureIdx] = 2.0*x + y - 4.0;
        return sol;
    }

    // see exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    NumEqVector rhsNewICNonSymmetrized_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        NumEqVector source(0.0);
        source[Indices::momentumXBalanceIdx] = -2.0;
        source[Indices::momentumYBalanceIdx] = 1.0;
        return source;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    static constexpr Scalar eps_ = 1e-7;
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManager_;
    TestCase testCase_;
};
} // end namespace Dumux

#endif // DUMUX_STOKES_SUBPROBLEM_HH
