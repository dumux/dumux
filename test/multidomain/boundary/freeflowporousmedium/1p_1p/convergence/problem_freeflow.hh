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

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dune/common/fmatrix.hh>

#include "testcase.hh"
#include "analyticalsolutions.hh"

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief The Stokes sub-problem of coupled Stokes-Darcy convergence test
 */
template <class TypeTag>
class FreeFlowSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using NumEqVector = typename ParentType::NumEqVector;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class BC {
        dirichlet, stress, mixed
    };

public:
    //! export the Indices
    using Indices = typename ModelTraits::Indices;

    FreeFlowSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                       std::shared_ptr<CouplingManager> couplingManager,
                       const DarcyStokesTestCase testCase,
                       const std::string& name)
    : ParentType(gridGeometry, couplingManager, "FreeFlow")
    , couplingManager_(couplingManager)
    , testCase_(testCase)
    {
        problemName_ = name + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        auto bc = getParamFromGroup<std::string>(this->paramGroup(), "Problem.BoundaryConditions", "Dirichlet");
        if (bc == "Dirichlet")
            boundaryConditions_ = BC::dirichlet;
        else if (bc == "Stress")
            boundaryConditions_ = BC::stress;
        else if (bc == "Mixed")
            boundaryConditions_ = BC::mixed;
        else
            DUNE_THROW(Dune::Exception, "Wrong BC type choose: Dirichlet, Stress or Mixed");

        std::cout << "Free flow domain: Using " << bc << " boundary conditions" << std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    const std::string& name() const
    { return problemName_; }

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
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
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            // TODO: is isFrontal really needed here?
            if (couplingManager().isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, scvf) && scvf.isFrontal())
            {
                values.setCouplingNeumann(Indices::momentumYBalanceIdx);
                values.setCouplingNeumann(Indices::momentumXBalanceIdx);
            }
            else
            {
                if (boundaryConditions_ == BC::dirichlet)
                    values.setAllDirichlet();
                else if (boundaryConditions_ == BC::stress)
                    values.setAllNeumann();
                else
                {
                    if (onLeftBoundary_(scvf.center()))
                        values.setAllNeumann();
                    else
                        values.setAllDirichlet();
                }
            }
        }
        else
        {
            if (couplingManager().isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex, scvf))
                values.setAllCouplingNeumann();
            else
            {
                if (boundaryConditions_ == BC::dirichlet)
                    values.setAllDirichlet();
                else if (boundaryConditions_ == BC::stress)
                    values.setAllNeumann();
                else
                {
                    if (onLeftBoundary_(scvf.center()))
                        values.setAllNeumann();
                    else
                        values.setAllDirichlet();
                }
            }
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;

        // momentum boundary conditions
        if constexpr (ParentType::isMomentumProblem())
        {
            if (couplingManager().isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, scvf))
            {
                values += couplingManager().momentumCouplingCondition(
                    CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex,
                    fvGeometry, scvf, elemVolVars
                );

                values += FluxHelper::slipVelocityMomentumFlux(
                    *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache
                );
            }
            else
            {
                const auto stress = stressTensorAtPos(globalPos);
                const auto n = scvf.unitOuterNormal();
                stress.mv(n, values);
            }
        }

        // mass boundary conditions
        else
        {
            if (couplingManager().isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex, scvf))
            {
                values = couplingManager().massCouplingCondition(
                    CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex,
                    fvGeometry, scvf, elemVolVars
                );
            }
            else
            {
                using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
                values = FluxHelper::scalarOutflowFlux(
                    *this, element, fvGeometry, scvf, elemVolVars
                );
            }
        }

        return values;
    }

    /*!
     * \brief Returns true if the scvf lies on a porous slip boundary
     */
    bool onSlipBoundary(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        return scvf.boundary()
            && couplingManager().isCoupled(
                CouplingManager::freeFlowMomentumIndex,
                CouplingManager::porousMediumIndex,
                scvf
            );
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
     * \brief Returns the sources within the domain.
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        const auto select = [&](const auto& rhs) -> NumEqVector {
            if constexpr (ParentType::isMomentumProblem())
                return { rhs[0], rhs[1] };
            else
                return { rhs[2] };
        };

        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                return select(ShiueOne::stokesRHS(globalPos));
            case DarcyStokesTestCase::ShiueExampleTwo:
                return select(ShiueTwo::stokesRHS(globalPos));
            case DarcyStokesTestCase::Rybak:
                return select(Rybak::stokesRHS(globalPos));
            case DarcyStokesTestCase::Schneider:
                return select(Schneider::stokesRHS(globalPos));
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
              for the Beavers-Joseph-Saffman boundary condition
     */
    auto permeability(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    { return couplingManager().darcyPermeability(fvGeometry, scvf); }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    { return couplingManager().problem(CouplingManager::porousMediumIndex).spatialParams().beaversJosephCoeffAtPos(scvf.ipGlobal()); }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     * Returns vector with entries: (velocity-x | velocity-y | pressure)
     */
    Dune::FieldVector<Scalar, 3> fullAnalyticalSolution(const GlobalPosition& globalPos) const
    {
        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                return ShiueOne::stokes(globalPos);
            case DarcyStokesTestCase::ShiueExampleTwo:
                return ShiueTwo::stokes(globalPos);
            case DarcyStokesTestCase::Rybak:
                return Rybak::stokes(globalPos);
            case DarcyStokesTestCase::Schneider:
                return Schneider::stokes(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     * Returns analytical solution depending on the type of sub-problem
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        using namespace Solution::DarcyStokes;
        const auto sol = fullAnalyticalSolution(globalPos);
        if constexpr (ParentType::isMomentumProblem())
            return { sol[0], sol[1] };
        else
            return { sol[2] };
    }

    // \}

private:
    /*!
     * \brief Returns the stress tensor within the domain.
     * \param globalPos The global position
     */
    Dune::FieldMatrix<double, 2, 2> stressTensorAtPos(const GlobalPosition &globalPos) const
    {
        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                return ShiueOne::stokesStress(globalPos);
            case DarcyStokesTestCase::ShiueExampleTwo:
                return ShiueTwo::stokesStress(globalPos);
            case DarcyStokesTestCase::Rybak:
                return Rybak::stokesStress(globalPos);
            case DarcyStokesTestCase::Schneider:
                return Schneider::stokesStress(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    static constexpr Scalar eps_ = 1e-7;
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManager_;
    DarcyStokesTestCase testCase_;
    BC boundaryConditions_;
};
} // end namespace Dumux

#endif
