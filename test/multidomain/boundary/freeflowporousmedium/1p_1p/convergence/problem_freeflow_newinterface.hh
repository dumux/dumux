// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The free-flow sub-problem of coupled FreeFlow/Darcy convergence test
 */

#ifndef DUMUX_FREEFLOW_SUBPROBLEM_NEWINTERFACE_HH
#define DUMUX_FREEFLOW_SUBPROBLEM_NEWINTERFACE_HH

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/common/constraintinfo.hh>

#include "testcase.hh"
#include "analyticalsolutions.hh"

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief The Stokes sub-problem of coupled Stokes-Darcy convergence test
 */
template <class TypeTag, class BaseProblem>
class FreeFlowSubProblemNewInterface : public BaseProblem
{
    using ParentType = BaseProblem;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<typename GridGeometry::GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class BC {
        dirichlet, stress, mixed
    };

public:
    //! export the Indices
    using Indices = typename ModelTraits::Indices;

    FreeFlowSubProblemNewInterface(std::shared_ptr<const GridGeometry> gridGeometry,
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

    // Call this function after problem and coupling manager init(...) has been called
    void setConstraints()
    {
        if constexpr (ParentType::isMomentumProblem())
        {
            if(!(boundaryConditions_ == BC::stress))
                appendDirichletConstraints_();
        }
    }

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

        if constexpr (ParentType::isMomentumProblem())
        {
            if(isMomentumFluxBoundary_(elemDisc, boundaryFace))
                values.setAllFluxBoundary();
        }
        else
            values.setAllFluxBoundary();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions at a given interpolation point.
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
        const auto& globalPos = faceIpData.global();

        // momentum boundary conditions
        if constexpr (ParentType::isMomentumProblem())
        {
            // We need to take care here: If the interpolation point lies both on the coupling interface and the
            // domain boundary (i.e., the leftmost and rightmost points of the interface), do not evaluate the coupling or slip condition
            // because they would require data from the boundary not available in this case (velocities for evaluating gradients,
            // those would only be available for Dirichlet BCs). Instead, directly use the given Neumann boundary stress.
            // TODO: Maybe couplingManager().isCoupled(...) could return false for these scvfs.
            //
            //                       | <-- left Neumann boundary
            //                       |
            // interpolation point --> o##### <-- scvf lying on coupling interface with interpolation point touching left Neumann boundary
            //
            if (isCoupled_(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, elemDisc, faceIpData)
                && globalPos[0] > this->gridGeometry().bBoxMin()[0] + eps_
                && globalPos[0] < this->gridGeometry().bBoxMax()[0] - eps_)
            {
                // normal momentum coupling condition
                const auto pmPressure = couplingManager().pmPressure(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex,
                                                                     elemDisc, elemVars, faceIpData);
                values += (pmPressure - this->referencePressure()) * faceIpData.unitOuterNormal();

                // tangential momentum coupling condition (slip condition)
                const auto& element = elemDisc.element();
                const auto elemSol = elementSolution(element, elemVars, elemDisc);
                const auto v = evalSolution(element, element.geometry(), elemDisc.gridGeometry(), elemSol, globalPos);
                // Slip is in x orientation
                GlobalPosition t(0.0);
                t[0] = 1.0;
                auto vt = v*t;
                values.axpy((this->effectiveViscosity(element, elemDisc, faceIpData)*this->betaBJ(elemDisc, faceIpData, t))*vt, t);
            }
            else
            {
                const auto stress = stressTensorAtPos(globalPos);
                const auto n = faceIpData.unitOuterNormal();
                stress.mv(n, values);
            }
        }
        // mass boundary conditions
        else
        {
            const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
            if (isCoupled_(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex, elemDisc, faceIpData))
            {
                values = 1.0/scvf.area() * couplingManager().massCouplingCondition(
                    CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex,
                    elemDisc, scvf, elemVars
                );
            }
            else
            {
                values =  this->velocity(elemDisc, faceIpData)
                    * elemVars[scvf.insideScvIdx()].density() * faceIpData.unitOuterNormal();
            }
        }

        return values;
    }

    /*!
     * \brief Return  Dirichlet boundary constraints and internal constraints.
     */
    const auto& constraints() const
    { return constraints_; }

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
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        const auto select = [&](const auto& rhs) -> Sources {
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
    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    { return InitialValues(0.0); }

    /*!
     * \brief Returns the beta value which is the alpha value divided by the square root of the (scalar-valued) interface permeability.
     */
    template<class IpData>
    Scalar betaBJ(const ElementDiscretization& elemDisc, const IpData& ipData, const GlobalPosition& tangentialVector) const
    {
        const auto& K = couplingManager().darcyPermeability(elemDisc, ipData);
        const auto alphaBJ = couplingManager().problem(CouplingManager::porousMediumIndex).spatialParams().beaversJosephCoeffAtPos(ipData.global());
        const auto interfacePermeability = vtmv(tangentialVector, K, tangentialVector);
        using std::sqrt;
        return alphaBJ / sqrt(interfacePermeability);
    }

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
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        using namespace Solution::DarcyStokes;
        const auto sol = fullAnalyticalSolution(globalPos);
        if constexpr (ParentType::isMomentumProblem())
            return { sol[0], sol[1] };
        else
            return { sol[2] };
    }

    /*!
     * \brief Returns the gradient of the analytical solution of the problem at a given position.
     * \param globalPos The global position
     */
    Dune::FieldVector<GlobalPosition, DirichletValues::dimension> gradAnalyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        Dune::FieldVector<GlobalPosition, DirichletValues::dimension> values;
        Dune::FieldMatrix<Scalar, 3, 2> grad(0.0);

        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                grad = ShiueOne::stokesGrad(globalPos);
                break;
            case DarcyStokesTestCase::ShiueExampleTwo:
                grad = ShiueTwo::stokesGrad(globalPos);
                break;
            case DarcyStokesTestCase::Rybak:
                grad = Rybak::stokesGrad(globalPos);
                break;
            case DarcyStokesTestCase::Schneider:
                grad = Schneider::stokesGrad(globalPos);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx][0] = grad[0][0];
            values[Indices::velocityXIdx][1] = grad[0][1];
            values[Indices::velocityYIdx][0] = grad[1][0];
            values[Indices::velocityYIdx][1] = grad[1][1];
        }
        else
        {
            values[Indices::pressureIdx][0] = grad[2][0];
            values[Indices::pressureIdx][1] = grad[2][1];
        }

        return values;
    }

    // \}

private:
    template<std::size_t i, std::size_t j, class FaceIpData>
    bool isCoupled_(Dune::index_constant<i> domainI,
                    Dune::index_constant<j> domainJ,
                    const ElementDiscretization& elemDisc,
                    const FaceIpData& faceIpData) const
    {
        if constexpr (requires { faceIpData.scvfIndex(); })
        {
            const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
            return couplingManager_->isCoupled(domainI, domainJ, elemDisc, elemDisc.intersectionIndex(scvf));
        }
        else if constexpr (requires { faceIpData.boundaryFaceIndex(); })
            return couplingManager_->isCoupled(domainI, domainJ, elemDisc, elemDisc.boundaryFace(faceIpData.boundaryFaceIndex()).intersectionIndex());
        else
            DUNE_THROW(Dune::InvalidStateException, "FaceIpData must provide either scvfIndex() or intersectionIndex()");
    }

    template<class BoundaryFace>
    bool isMomentumFluxBoundary_(const ElementDiscretization& elemDisc,
                                 const BoundaryFace& boundaryFace) const
    {
        if(couplingManager().isCoupled(CouplingManager::freeFlowMomentumIndex,
                                       CouplingManager::porousMediumIndex,
                                       elemDisc,
                                       boundaryFace.intersectionIndex()))
        {
            return true;
        }
        else if (boundaryConditions_ == BC::dirichlet)
            return false;
        else if (boundaryConditions_ == BC::stress)
            return true;
        else
        {
            if (onLeftBoundary_(boundaryFace.center()))
                return true;
            else
                return false;
        }
    }

    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);

            for(const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if(!isMomentumFluxBoundary_(elemDisc, boundaryFace))
                {
                    for(const auto& localDof : localDofs(elemDisc, boundaryFace))
                    {
                        const auto& globalPos = ipData(elemDisc, localDof).global();
                        ConstraintInfo info;
                        info.setAll();

                        DirichletValues dirichletValues(this->analyticalSolution(globalPos));
                        constraints_.push_back(DirichletConstraintData{std::move(info), std::move(dirichletValues), localDof.dofIndex()});
                    }
                }
            }
        }
    }

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
    std::vector<DirichletConstraintData> constraints_;
};
} // end namespace Dumux

#endif
