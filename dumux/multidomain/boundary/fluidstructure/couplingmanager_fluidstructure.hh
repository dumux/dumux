// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Coupling manager for the Biot model
 */
#ifndef DUMUX_BIOT_COUPLINGMANAGER_HH
#define DUMUX_BIOT_COUPLINGMANAGER_HH

#include <memory>
#include <tuple>
#include <vector>
#include <array>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

template<class MDTraits>
class BiotCouplingManager
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    // the sub domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using Indices
        = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>::VolumeVariables::Indices;

    template<std::size_t id> using CouplingStencil = std::vector<GridIndex<id>>;

    using GlobalPosition = typename Element<0>::Geometry::GlobalCoordinate;

public:
    static constexpr auto mechanicsIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto flowIdx = typename MDTraits::template SubDomain<1>::Index();

    //! export traits
    using MultiDomainTraits = MDTraits;

    BiotCouplingManager() = default;

    BiotCouplingManager(std::shared_ptr<GridGeometry<mechanicsIdx>> mechanicsGG,
                        std::shared_ptr<GridGeometry<flowIdx>> flowGG)
    {
        const auto& mechanicsGridGeometry = *mechanicsGG;
        const auto& flowGridGeometry = *flowGG;
        auto mechanicsFVGeometry = localView(mechanicsGridGeometry);
        auto flowFVGeometry = localView(flowGridGeometry);

        flowToMechStencils_.clear();
        flowToMechStencils_.resize(flowGridGeometry.gridView().size(0));

        mechToFlowStencils_.clear();
        mechToFlowStencils_.resize(mechanicsGridGeometry.gridView().size(0));

        assert(flowToMechStencils_.size() == mechToFlowStencils_.size());

        for (const auto& element : elements(mechanicsGridGeometry.gridView()))
        {
            mechanicsFVGeometry.bindElement(element);
            flowFVGeometry.bindElement(element);
            const auto eIdx = mechanicsFVGeometry.elementIndex();

            for (const auto& scv : scvs(mechanicsFVGeometry))
                flowToMechStencils_[eIdx].push_back(scv.dofIndex());

            for (const auto& scv : scvs(flowFVGeometry))
                mechToFlowStencils_[eIdx].push_back(scv.dofIndex());
        }
    }

    void init(std::shared_ptr<Problem<mechanicsIdx>> mechanicsProblem,
              std::shared_ptr<Problem<flowIdx>> flowProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        this->setSubProblems(std::make_tuple(mechanicsProblem, flowProblem));
    }

    template<std::size_t i, std::size_t j>
    const std::vector<std::size_t>& couplingStencil(Dune::index_constant<i> domainI,
                                                    const Element<i>& element,
                                                    Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if constexpr (domainI == mechanicsIdx)
            return mechToFlowStencils_[eIdx];
        else
            return flowToMechStencils_[eIdx];
    }

    auto divU(typename GridGeometry<flowIdx>::LocalView const& fvGeometry,
              typename GridGeometry<flowIdx>::SubControlVolume const& scv) const
    {
        const auto& gg = this->problem(mechanicsIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(mechanicsIdx), gg);

        const auto gradU = evalGradients(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            scv.center()
        );

        double divU = 0.0;
        for (int i = 0; i < gradU.size(); ++i)
            divU += gradU[i][i];
        return divU;
    }

    auto displacement(typename GridGeometry<flowIdx>::LocalView const& fvGeometry,
                      typename GridGeometry<flowIdx>::SubControlVolumeFace const& scvf) const
    {
        const auto& gg = this->problem(mechanicsIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(mechanicsIdx), gg);
        return evalSolution(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            scvf.ipGlobal()
        );
    }

    auto totalPressure(typename GridGeometry<mechanicsIdx>::LocalView const& fvGeometry,
                       typename GridGeometry<mechanicsIdx>::SubControlVolumeFace const& scvf) const
    {
        const auto& gg = this->problem(flowIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(flowIdx), gg);
        return evalSolution(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            scvf.ipGlobal()
        )[1];
    }

    /*!
     * \brief the solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    template<std::size_t i>
    auto& curSol(Dune::index_constant<i> domainIdx)
    { return ParentType::curSol(domainIdx); }

    /*!
     * \brief the solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    template<std::size_t i>
    const auto& curSol(Dune::index_constant<i> domainIdx) const
    { return ParentType::curSol(domainIdx); }

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    { elementSets_ = computeColoring(this->problem(flowIdx).gridGeometry()).sets; }

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
            DUNE_THROW(Dune::InvalidStateException,
                "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel
        // for this we have to color the elements so that we don't get
        // race conditions when writing into the global matrix or modifying grid variable caches
        // each color can be assembled using multiple threads
        const auto& grid = this->problem(domainId).gridGeometry().gridView().grid();
        for (const auto& elements : elementSets_)
        {
            Dumux::parallelFor(elements.size(), [&](const std::size_t n)
            {
                const auto element = grid.entity(elements[n]);
                assembleElement(element);
            });
        }
    }

protected:
    using ParentType::curSol;

private:
    //! coloring for multithreaded assembly
    std::deque<std::vector<ElementSeed<mechanicsIdx>>> elementSets_;

    std::vector<std::vector<std::size_t>> mechToFlowStencils_;
    std::vector<std::vector<std::size_t>> flowToMechStencils_;
};

// we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<BiotCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
