// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup KirchhoffLovePlate
 * \brief Coupling manager for the Kirchhoff-Love model
 */
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_COUPLINGMANAGER_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_COUPLINGMANAGER_HH

#include <memory>
#include <tuple>
#include <vector>
#include <array>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

template<class MDTraits>
class KirchhoffLovePlateCouplingManager
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
    static constexpr auto rotationIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto deformationIdx = typename MDTraits::template SubDomain<1>::Index();

    //! export traits
    using MultiDomainTraits = MDTraits;

    KirchhoffLovePlateCouplingManager() = default;

    KirchhoffLovePlateCouplingManager(std::shared_ptr<GridGeometry<rotationIdx>> rotationsGG,
                                      std::shared_ptr<GridGeometry<deformationIdx>> deformationGG)
    {
        const auto& rotationsGridGeometry = *rotationsGG;
        const auto& deformationGridGeometry = *deformationGG;
        auto rotGeo = localView(rotationsGridGeometry);
        auto defGeo = localView(deformationGridGeometry);

        defToRotStencils_.clear();
        defToRotStencils_.resize(deformationGridGeometry.gridView().size(0));

        rotToDefStencils_.clear();
        rotToDefStencils_.resize(rotationsGridGeometry.gridView().size(0));

        assert(defToRotStencils_.size() == rotToDefStencils_.size());

        for (const auto& element : elements(rotationsGridGeometry.gridView()))
        {
            rotGeo.bindElement(element);
            defGeo.bindElement(element);
            const auto eIdx = rotGeo.elementIndex();

            for (const auto& localDof : localDofs(rotGeo))
                defToRotStencils_[eIdx].push_back(localDof.dofIndex());

            for (const auto& localDof : localDofs(defGeo))
                rotToDefStencils_[eIdx].push_back(localDof.dofIndex());
        }
    }

    void init(std::shared_ptr<Problem<rotationIdx>> momentumProblem,
              std::shared_ptr<Problem<deformationIdx>> massProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
    }

    template<std::size_t i, std::size_t j>
    const CouplingStencil<j>& couplingStencil(Dune::index_constant<i> domainI,
                                              const Element<i>& element,
                                              Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if constexpr (domainI == rotationIdx)
            return rotToDefStencils_[eIdx];
        else
            return defToRotStencils_[eIdx];
    }

    auto rotation(typename GridGeometry<deformationIdx>::LocalView const& fvGeometry,
                  typename GridGeometry<deformationIdx>::SubControlVolumeFace const& scvf) const
    {
        const auto& gg = this->problem(rotationIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(rotationIdx), gg);
        return evalSolution(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            scvf.ipGlobal()
        );
    }

    auto deformationAndPotentials(typename GridGeometry<rotationIdx>::LocalView const& fvGeometry,
                                  typename GridGeometry<rotationIdx>::SubControlVolumeFace const& scvf) const
    {
        const auto& gg = this->problem(deformationIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(deformationIdx), gg);
        return evalSolution(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            scvf.ipGlobal()
        );
    }

    auto shearCurlPotentialIdx() const
    { return Indices<deformationIdx>::shearCurlPotentialIdx; }

    auto shearGradPotentialIdx() const
    { return Indices<deformationIdx>::shearGradPotentialIdx; }

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
    { elementSets_ = computeColoring(this->problem(deformationIdx).gridGeometry()).sets; }

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
    std::deque<std::vector<ElementSeed<deformationIdx>>> elementSets_;

    std::vector<CouplingStencil<deformationIdx>> rotToDefStencils_;
    std::vector<CouplingStencil<rotationIdx>> defToRotStencils_;
};

// we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<KirchhoffLovePlateCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
