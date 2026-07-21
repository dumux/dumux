// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Free-flow coupling manager for pure FE (Galerkin) mass-momentum coupling.
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_FE_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_FE_HH

#include <memory>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>

#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include "typetraits.hh"

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief Coupling manager for pure-FE Navier-Stokes mass-momentum coupling.
 * \note SubDomain<0> = momentum (velocity, e.g. PQ2FEModel), SubDomain<1> = mass (pressure, e.g. PQ1FEModel)
 */
template<class Traits>
class FEFreeFlowCouplingManager
: public CouplingManager<Traits>
{
    using ParentType = CouplingManager<Traits>;
public:
    static constexpr auto freeFlowMomentumIndex = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIndex = typename Traits::template SubDomain<1>::Index();

    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using GridVariables = typename Traits::template SubDomain<id>::GridVariables;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;

    using Scalar = typename Traits::Scalar;
    using SolutionVector = typename Traits::SolutionVector;
    using CouplingStencilType = std::vector<std::size_t>;
    using GridVariablesTuple = typename Traits::template TupleOfSharedPtr<GridVariables>;

    using GlobalPosition = typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate;
    using VelocityVector = GlobalPosition;
    static constexpr int dim = GridView<freeFlowMomentumIndex>::dimension;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

    using MomentumDiscretizationMethod = typename GridGeometry<freeFlowMomentumIndex>::DiscretizationMethod;
    using MassDiscretizationMethod = typename GridGeometry<freeFlowMassIndex>::DiscretizationMethod;

public:

    //! index of the pressure primary variable in the mass subdomain
    static constexpr auto pressureIdx =
        GetPropType<SubDomainTypeTag<freeFlowMassIndex>, Properties::ModelTraits>::Indices::pressureIdx;

    /*!
     * \name Setup
     */
    // \{

    //! use as regular coupling manager (stationary)
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        gridVariables_ = gridVariables;
        this->updateSolution(curSol);
        computeCouplingStencils_();
    }

    //! use as regular coupling manager (transient)
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol)
    {
        init(momentumProblem, massProblem, std::forward<GridVariablesTuple>(gridVariables), curSol);
        prevSol_ = &prevSol;
        isTransient_ = true;
    }

    //! use as binary coupling manager in a multi-model context
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const typename ParentType::SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        gridVariables_ = gridVariables;
        this->attachSolution(curSol);
        computeCouplingStencils_();
    }

    // \}

    /*!
     * \name Quantities the momentum subdomain reads from the mass (pressure) subdomain
     */
    // \{

    //! Interpolated pressure of the mass subdomain at an interpolation point on the momentum element
    template<class IpData>
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const IpData& ipData,
                    const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        const auto& gg = this->problem(freeFlowMassIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMassIndex]
                                                   :  this->curSol(freeFlowMassIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        return evalSolutionAtLocalPos(element, element.geometry(), gg, elemSol, ipData.local())[pressureIdx];
    }

    //! Density evaluated at an interpolation point (constant one-phase fluid for the Stokes stage)
    template<class IpData>
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const IpData& ipData,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->problem(freeFlowMassIndex).densityAtPos(ipData.global());
    }

    //! Effective viscosity evaluated at an interpolation point
    template<class IpData>
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const IpData& ipData,
                              const bool considerPreviousTimeStep = false) const
    {
        return this->problem(freeFlowMassIndex).effectiveViscosityAtPos(ipData.global());
    }

    // \}

    /*!
     * \name Quantities the mass (continuity) subdomain reads from the momentum subdomain
     */
    // \{

    //! Interpolated velocity of the momentum subdomain at an interpolation point on the mass element
    template<class IpData>
    VelocityVector velocity(const FVElementGeometry<freeFlowMassIndex>& fvGeometry,
                            const IpData& ipData,
                            const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        const auto& element = fvGeometry.element();
        const auto& gg = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMomentumIndex]
                                                   :  this->curSol(freeFlowMomentumIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        const auto priVars = evalSolutionAtLocalPos(element, element.geometry(), gg, elemSol, ipData.local());
        VelocityVector v(0.0);
        for (int i = 0; i < dim; ++i)
            v[i] = priVars[i];
        return v;
    }

    //! Interpolated momentum velocity at the element center
    VelocityVector elementVelocity(const FVElementGeometry<freeFlowMassIndex>& fvGeometry) const
    {
        const auto& element = fvGeometry.element();
        const auto& gg = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto elemSol = elementSolution(element, this->curSol(freeFlowMomentumIndex), gg);
        const auto center = referenceElement(element).position(0, 0);
        const auto priVars = evalSolutionAtLocalPos(element, element.geometry(), gg, elemSol, center);
        VelocityVector v(0.0);
        for (int i = 0; i < dim; ++i)
            v[i] = priVars[i];
        return v;
    }

    //! Divergence of the momentum velocity at an interpolation point on the mass element
    template<class IpData>
    Scalar velocityDivergence(const FVElementGeometry<freeFlowMassIndex>& fvGeometry,
                              const IpData& ipData,
                              const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        const auto& element = fvGeometry.element();
        const auto& gg = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMomentumIndex]
                                                   :  this->curSol(freeFlowMomentumIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        const auto grads = evalGradientsAtLocalPos(element, element.geometry(), gg, elemSol, ipData.local());
        Scalar divV = 0.0;
        for (int i = 0; i < dim; ++i)
            divV += grads[i][i]; // grads[eqIdx][dir] = d(priVar_eqIdx)/d(x_dir)
        return divV;
    }

    // \}

    /*!
     * \name Coupling stencils
     */
    // \{

    //! mass element residual depends on all momentum dofs of that element
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMassIndex> domainI,
                                               const Element<freeFlowMassIndex>& elementI,
                                               Dune::index_constant<freeFlowMomentumIndex> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMassIndex).gridGeometry().elementMapper().index(elementI);
        return massToMomentumStencils_[eIdx];
    }

    //! momentum element residual depends on all mass (pressure) dofs of that element
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                               const Element<freeFlowMomentumIndex>& elementI,
                                               Dune::index_constant<freeFlowMassIndex> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMomentumIndex).gridGeometry().elementMapper().index(elementI);
        return momentumToMassStencils_[eIdx];
    }

    // \}

    /*!
     * \name Coupling context
     */
    // \{

    //! update the solution when a coupled dof changes (no cached cross-domain volvars for pure FE)
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    // \}

    /*!
     * \name Multithreaded assembly
     */
    // \{

    void computeColorsForAssembly()
    {
        elementSets_ = computeColoring(this->problem(freeFlowMomentumIndex).gridGeometry()).sets;
    }

    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        const auto& grid = this->problem(freeFlowMomentumIndex).gridGeometry().gridView().grid();
        for (const auto& elements : elementSets_)
            Dumux::parallelFor(elements.size(), [&](const std::size_t k)
            {
                const auto element = grid.entity(elements[k]);
                assembleElement(element);
            });
    }

    // \}

private:
    void computeCouplingStencils_()
    {
        const auto& momentumGG = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto& massGG = this->problem(freeFlowMassIndex).gridGeometry();
        auto momentumFvGeometry = localView(momentumGG);
        auto massFvGeometry = localView(massGG);

        massToMomentumStencils_.clear();
        massToMomentumStencils_.resize(massGG.gridView().size(0));
        momentumToMassStencils_.clear();
        momentumToMassStencils_.resize(momentumGG.gridView().size(0));

        for (const auto& element : elements(momentumGG.gridView()))
        {
            momentumFvGeometry.bindElement(element);
            massFvGeometry.bindElement(element);
            const auto eIdx = momentumFvGeometry.elementIndex();

            for (const auto& localDof : localDofs(momentumFvGeometry))
                massToMomentumStencils_[eIdx].push_back(localDof.dofIndex());

            for (const auto& localDof : localDofs(massFvGeometry))
                momentumToMassStencils_[eIdx].push_back(localDof.dofIndex());
        }
    }

    std::vector<CouplingStencilType> momentumToMassStencils_;
    std::vector<CouplingStencilType> massToMomentumStencils_;

    GridVariablesTuple gridVariables_;
    const SolutionVector* prevSol_ = nullptr;
    bool isTransient_ = false;
    std::deque<std::vector<ElementSeed<freeFlowMomentumIndex>>> elementSets_;
};

} // end namespace Dumux

#endif
