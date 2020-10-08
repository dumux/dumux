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
 * \brief The Darcy sub-problem of coupled Stokes-Darcy BJ/NewIC test
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct DarcyOneP { using InheritsFrom = std::tuple<OneP, BoxModel>; };
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
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy BJ/NewIC test
 */
template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
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

    enum class TestCase
    {
        BJSymmetrized, NewICNonSymmetrized
    };

public:
    //! export the Indices
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Darcy"), couplingManager_(couplingManager)
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
      * \param scv The boundary sub control volume
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        //corners will not be handled as coupled
        if(onLeftBoundary_(scv.dofPosition()) || onRightBoundary_(scv.dofPosition()))
        {
            values.setAllDirichlet();
            return values;
        }

        auto fvGeometry = localView(this->gridGeometry());
        fvGeometry.bindElement(element);
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (couplingManager().isCoupledEntity(CouplingManager::porousMediumIdx, element, scvf))
            {
                values.setAllCouplingNeumann();
            }
        }
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        const auto p = analyticalSolution(globalPos)[pressureIdx];
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

        if (couplingManager().isCoupledEntity(CouplingManager::porousMediumIdx, element, scvf))
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

        return values;
    }

    // \}

    /*!
    * \name Volume terms
    */
    // \{
    /*!
    * \brief Evaluate the source term for all phases within a given
    *        sub-control-volume.
    *
    * \param element The element for which the source term is set
    * \param fvGeomentry The fvGeometry
    * \param elemVolVars The element volume variables
    * \param scv The subcontrolvolume
    */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        GlobalPosition globalPos = scv.geometry().center();
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
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
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
            case TestCase::BJSymmetrized:
                return analyticalSolutionBJSymmetrized_(globalPos);
            case TestCase::NewICNonSymmetrized:
                return analyticalSolutionNewICNonSymmetrized_(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:

    // see exact solution for BJ-IC with symmetrized stress tensor (by Elissa Eggenweiler)
    Dune::FieldVector<Scalar, 3> analyticalSolutionBJSymmetrized_(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::exp; using std::sin; using std::cos;
        sol[velocityXIdx] = 2.0*x*x + 2.0*x*y - 6.0*x - 2.0*y*y + 3.0*y - 5.0;
        sol[velocityYIdx] = 1.0*x*x + 4.0*x*y - 9.0*x - 1.0*y*y - 1.0;
        sol[pressureIdx] = -x*(x-1.0)*(y-1.0) + 1.0/3.0*(y-1.0)*(y-1.0)*(y-1.0) + 4.0 + 2.0*y + x*x;
        return sol;
    }

    // see exact solution for BJ-IC with symmetrized stress tensor (by Elissa Eggenweiler)
    NumEqVector rhsBJSymmetrized_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        using std::sin;
        return NumEqVector(8.0*x - 6.0);
    }


    // see exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    Dune::FieldVector<Scalar, 3> analyticalSolutionNewICNonSymmetrized_(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::exp; using std::sin; using std::cos;
        sol[velocityXIdx] = x*(y-1.0) + (x-1.0)*(y-1.0) - 3.0;
        sol[velocityYIdx] = x*(x-1.0) - (y-1.0)*(y-1.0) - 2.0;
        sol[pressureIdx] = x*(1.0-x)*(y-1.0) + 1.0/3.0*(y-1.0)*(y-1.0)*(y-1.0) + 3.0*x + 2.0*y + 1.0;
        return sol;
    }

    // see exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    NumEqVector rhsNewICNonSymmetrized_(const GlobalPosition& globalPos) const
    {
        return NumEqVector(0.0);
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;  }
    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;  }
    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;  }
    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;  }

    static constexpr Scalar eps_ = 1e-7;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::string problemName_;
    TestCase testCase_;
};
} // end namespace Dumux

#endif //DUMUX_DARCY_SUBPROBLEM_HH
