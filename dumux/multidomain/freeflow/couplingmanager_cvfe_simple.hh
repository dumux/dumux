// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Freeflow coupling managers (Navier-Stokes mass-momentum coupling)
 * Simpler implementation of the coupling manager
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_CVFE_SIMPLE_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_CVFE_SIMPLE_HH

#include <memory>
#include <tuple>
#include <vector>
#include <deque>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include "typetraits.hh"

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for free flow systems
 * \note coupling manager for control volume finite element schemes
 * Simpler implementation of the coupling manager
 */
template<class Traits>
class CVFEFreeFlowCouplingManagerSimple
: public CouplingManager<Traits>
{
    using ParentType = CouplingManager<Traits>;
public:
    static constexpr auto freeFlowMomentumIndex = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIndex = typename Traits::template SubDomain<1>::Index();

    // this can be used if the coupling manager is used inside a meta-coupling manager (e.g. multi-binary)
    // to manager the solution vector storage outside this class
    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridVariables = typename Traits::template SubDomain<id>::GridVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using GridFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using VolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::VolumeVariables>;

    using Scalar = typename Traits::Scalar;
    using SolutionVector = typename Traits::SolutionVector;

    using CouplingStencilType = std::vector<std::size_t>;

    using GridVariablesTuple = typename Traits::template TupleOfSharedPtr<GridVariables>;

    using FluidSystem = typename VolumeVariables<freeFlowMassIndex>::FluidSystem;

    using VelocityVector = typename SubControlVolumeFace<freeFlowMassIndex>::GlobalPosition;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

    static_assert(std::is_same_v<VelocityVector, typename SubControlVolumeFace<freeFlowMomentumIndex>::GlobalPosition>);

    using MomentumDiscretizationMethod = typename GridGeometry<freeFlowMomentumIndex>::DiscretizationMethod;
    using MassDiscretizationMethod = typename GridGeometry<freeFlowMassIndex>::DiscretizationMethod;

public:

    static constexpr auto pressureIdx = VolumeVariables<freeFlowMassIndex>::Indices::pressureIdx;

    //! use as regular coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol)
    {
        init(momentumProblem, massProblem, curSol);
    }

    //! use as regular coupling manager in a transient setting
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

    //! use as regular coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              const SolutionVector& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        this->updateSolution(curSol);

        computeCouplingStencils_();
    }

    //! use as sub coupling manager in a meta coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              const typename ParentType::SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        this->attachSolution(curSol);

        computeCouplingStencils_();
    }

    /*!
     * \brief Returns the pressure at a given sub control volume face
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                    const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        const auto& gg = this->problem(freeFlowMassIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMassIndex]
                                                   :  this->curSol(freeFlowMassIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        return evalSolution(element, element.geometry(), gg, elemSol, scvf.ipGlobal())[pressureIdx];
    }

    /*!
     * \brief Returns the pressure at a given sub control volume
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolume<freeFlowMomentumIndex>& scv,
                    const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        const auto& gg = this->problem(freeFlowMassIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMassIndex]
                                                   :  this->curSol(freeFlowMassIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        return evalSolution(element, element.geometry(), gg, elemSol, scv.dofPosition())[pressureIdx];
    }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    VelocityVector faceVelocity(const Element<freeFlowMassIndex>& element,
                                const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        const auto geometry = element.geometry();
        const auto& gg = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto elemSol = elementSolution(element, this->curSol(freeFlowMomentumIndex), gg);

        return evalSolution(
            element,
            geometry,
            gg, elemSol,
            scvf.ipGlobal()
        );
    }

    /*!
     * \brief Returns the velocity at the element center.
     */
    VelocityVector elementVelocity(const FVElementGeometry<freeFlowMassIndex>& fvGeometry) const
    {
        const auto& element = fvGeometry.element();
        const auto geometry = element.geometry();
        const auto& gg = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto elemSol = elementSolution(element, this->curSol(freeFlowMomentumIndex), gg);

        return evalSolution(
            element,
            geometry,
            gg, elemSol,
            geometry.center()
        );
    }

    template<std::size_t j>
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                               const Element<freeFlowMomentumIndex>& elementI,
                                               const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                               Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMassIndex> domainI,
                                              const Element<freeFlowMassIndex>& elementI,
                                              Dune::index_constant<freeFlowMomentumIndex> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMassIndex).gridGeometry().elementMapper().index(elementI);
        return massAndEnergyToMomentumStencils_[eIdx];
    }

    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                               const Element<freeFlowMomentumIndex>& elementI,
                                               Dune::index_constant<freeFlowMassIndex> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMomentumIndex).gridGeometry().elementMapper().index(elementI);
        return momentumToMassAndEnergyStencils_[eIdx];
    }

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

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    {
        if constexpr (MomentumDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            // use coloring of the mass discretization for both domains
            // the diamond coloring is a subset (minimum amount of colors) of cctpfa/box coloring
            elementSets_ = computeColoring(this->problem(freeFlowMassIndex).gridGeometry()).sets;
        }
        else
        {
            // use coloring of the momentum discretization for both domains
            elementSets_ = computeColoring(this->problem(freeFlowMomentumIndex).gridGeometry()).sets;
        }
    }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainId the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel
        // for this we have to color the elements so that we don't get
        // race conditions when writing into the global matrix
        // each color can be assembled using multiple threads
        const auto& grid = this->problem(freeFlowMassIndex).gridGeometry().gridView().grid();
        for (const auto& elements : elementSets_)
        {
            Dumux::parallelFor(elements.size(), [&](const std::size_t eIdx)
            {
                const auto element = grid.entity(elements[eIdx]);
                assembleElement(element);
            });
        }
    }

private:

    void computeCouplingStencils_()
    {
        const auto& momentumGridGeometry = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto& massGridGeometry = this->problem(freeFlowMassIndex).gridGeometry();
        auto momentumFvGeometry = localView(momentumGridGeometry);
        auto massFvGeometry = localView(massGridGeometry);

        massAndEnergyToMomentumStencils_.clear();
        massAndEnergyToMomentumStencils_.resize(massGridGeometry.gridView().size(0));

        momentumToMassAndEnergyStencils_.clear();
        momentumToMassAndEnergyStencils_.resize(momentumGridGeometry.gridView().size(0));

        assert(massAndEnergyToMomentumStencils_.size() == momentumToMassAndEnergyStencils_.size());

        for (const auto& element : elements(momentumGridGeometry.gridView()))
        {
            momentumFvGeometry.bindElement(element);
            massFvGeometry.bindElement(element);
            const auto eIdx = momentumFvGeometry.elementIndex();

            for (const auto& scv : scvs(momentumFvGeometry))
                massAndEnergyToMomentumStencils_[eIdx].push_back(scv.dofIndex());

            for (const auto& scv : scvs(massFvGeometry))
                momentumToMassAndEnergyStencils_[eIdx].push_back(scv.dofIndex());
        }
    }

    CouplingStencilType emptyStencil_;
    std::vector<CouplingStencilType> momentumToMassAndEnergyStencils_;
    std::vector<CouplingStencilType> massAndEnergyToMomentumStencils_;

    const SolutionVector* prevSol_;
    bool isTransient_;

    std::deque<std::vector<ElementSeed<freeFlowMomentumIndex>>> elementSets_;
};

namespace Detail {

// declaration (specialize for different discretization types)
template<class Traits, class DiscretizationMethod = typename Detail::MomentumDiscretizationMethod<Traits>::type>
struct CouplingManagerSupportsMultithreadedAssemblyCVFESelector;

template<class Traits, class D>
struct CouplingManagerSupportsMultithreadedAssemblyCVFESelector<Traits, DiscretizationMethods::CVFE<D>>
{ using type = std::true_type; };

} // end namespace Detail

//! whether we support multithreaded assembly
template<class T>
struct CouplingManagerSupportsMultithreadedAssembly<CVFEFreeFlowCouplingManagerSimple<T>>
: public Detail::CouplingManagerSupportsMultithreadedAssemblyCVFESelector<T>::type
{};

} // end namespace Dumux

#endif
