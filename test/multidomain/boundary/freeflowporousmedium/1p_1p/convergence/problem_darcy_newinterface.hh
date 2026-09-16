// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test stated in terms of the new interfaces
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_NEWINTERFACE_HH
#define DUMUX_DARCY_SUBPROBLEM_NEWINTERFACE_HH

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

#include "testcase.hh"
#include "analyticalsolutions.hh"

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test stated in terms of the new interfaces
 *
 * The analytic solution, the source term and the boundary values are those of DarcySubProblem.
 * Dirichlet values reach the assembler as constraints, which the problem collects once over the
 * boundary faces that carry no flux condition.
 */
template <class TypeTag>
class DarcySubProblemNewInterface : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;
    using BoundaryFluxes = typename ParentType::Traits::NumEqVector;

    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class BC {
        dirichlet, neumann, mixed
    };

public:
    //! export the Indices
    using Indices = typename ModelTraits::Indices;

    DarcySubProblemNewInterface(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                std::shared_ptr<CouplingManager> couplingManager,
                                std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                                const DarcyStokesTestCase testCase,
                                const std::string& name)
    : ParentType(gridDiscretization, spatialParams, "Darcy")
    , couplingManager_(couplingManager)
    , testCase_(testCase)
    {
        problemName_ = name + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        auto bc = getParamFromGroup<std::string>(this->paramGroup(), "Problem.BoundaryConditions", "Dirichlet");
        if (bc == "Dirichlet")
            boundaryConditions_ = BC::dirichlet;
        else if (bc == "Neumann")
            boundaryConditions_ = BC::neumann;
        else if (bc == "Mixed")
            boundaryConditions_ = BC::mixed;
        else
            DUNE_THROW(Dune::Exception, "Wrong BC type choose: Dirichlet, Neumann or Mixed");

        std::cout << "Porous medium domain: Using " << bc << " boundary conditions" << std::endl;
    }

    //! Call this function after the problem and the coupling manager have been initialized
    void setConstraints()
    { appendDirichletConstraints_(); }

    /*!
     * \name Problem parameters
     */
    // \{

    const std::string& name() const
    { return problemName_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary face.
     *
     * \param elemDisc The element discretization
     * \param boundaryFace The boundary face
     */
    BoundaryTypes boundaryTypes(const ElementDiscretization& elemDisc,
                                const typename ElementDiscretization::BoundaryFace& boundaryFace) const
    {
        BoundaryTypes values;

        if (isFluxBoundary_(elemDisc, boundaryFace))
            values.setAllFluxBoundary();

        return values;
    }

    /*!
     * \brief Evaluates the boundary flux at a given interpolation point.
     *
     * \param elemDisc The element discretization
     * \param elemVars All variables for the element
     * \param faceIpData Face interpolation point data
     */
    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    {
        BoundaryFluxes values(0.0);

        const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
        if (couplingManager().isCoupled(CouplingManager::porousMediumIndex,
                                        CouplingManager::freeFlowMassIndex,
                                        elemDisc, elemDisc.intersectionIndex(scvf)))
            values[Indices::conti0EqIdx] = 1.0/scvf.area() * couplingManager().massCouplingCondition(
                CouplingManager::porousMediumIndex, CouplingManager::freeFlowMassIndex,
                elemDisc, scvf, elemVars
            );
        else
        {
            const auto sol = fullAnalyticalSolution(faceIpData.global());
            const auto n = faceIpData.unitOuterNormal();
            auto v = n; v[0] = sol[0]; v[1] = sol[1];
            values[Indices::conti0EqIdx] = v*n;
        }

        return values;
    }

    //! The Dirichlet values the assembler constrains the solution to
    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                return ShiueOne::darcyRHS(globalPos);
            case DarcyStokesTestCase::ShiueExampleTwo:
                return ShiueTwo::darcyRHS(globalPos);
            case DarcyStokesTestCase::Rybak:
                return Rybak::darcyRHS(globalPos);
            case DarcyStokesTestCase::Schneider:
                return Schneider::darcyRHS(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    /*!
     * \brief Evaluates the initial value for a position.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& pos) const
    { return PrimaryVariables(0.0); }

    // \}

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
                return ShiueOne::darcy(globalPos);
            case DarcyStokesTestCase::ShiueExampleTwo:
                return ShiueTwo::darcy(globalPos);
            case DarcyStokesTestCase::Rybak:
                return Rybak::darcy(globalPos);
            case DarcyStokesTestCase::Schneider:
                return Schneider::darcy(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        const auto sol = fullAnalyticalSolution(globalPos);
        return { sol[2] };
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    template<class BoundaryFace>
    bool isFluxBoundary_(const ElementDiscretization& elemDisc,
                         const BoundaryFace& boundaryFace) const
    {
        if (couplingManager().isCoupled(CouplingManager::porousMediumIndex,
                                        CouplingManager::freeFlowMassIndex,
                                        elemDisc,
                                        boundaryFace.intersectionIndex()))
            return true;
        else if (boundaryConditions_ == BC::dirichlet)
            return false;
        else if (boundaryConditions_ == BC::neumann)
            return true;
        else
            return onLeftBoundary_(boundaryFace.center());
    }

    //! Every local dof on a boundary face carrying no flux condition is constrained to the analytical pressure
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);

            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if (isFluxBoundary_(elemDisc, boundaryFace))
                    continue;

                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto globalPos = ipData(elemDisc, localDof).global();

                    ConstraintInfo info;
                    info.setAll();

                    ConstraintValues values(analyticalSolution(globalPos));
                    constraints_.push_back(
                        DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()}
                    );
                }
            }
        }
    }

    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridDiscretization().bBoxMin()[0] + eps_; }

    static constexpr Scalar eps_ = 1e-7;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::string problemName_;
    DarcyStokesTestCase testCase_;
    BC boundaryConditions_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
