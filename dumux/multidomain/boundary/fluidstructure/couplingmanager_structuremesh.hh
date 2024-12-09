// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FLUIDSTRUCTURE_COUPLINGMANAGER_STRUCTURE_MESH_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FLUIDSTRUCTURE_COUPLINGMANAGER_STRUCTURE_MESH_HH

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
class StructureMeshCouplingManager
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
    static constexpr auto structureIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto meshMotionIdx = typename MDTraits::template SubDomain<1>::Index();

    //! export traits
    using MultiDomainTraits = MDTraits;

    StructureMeshCouplingManager() = default;

    StructureMeshCouplingManager(std::shared_ptr<GridGeometry<structureIdx>> structureGG,
                                 std::shared_ptr<GridGeometry<meshMotionIdx>> meshMotionGG)
    {
        const auto& structureGridGeometry = *structureGG;
        const auto& meshMotionGridGeometry = *meshMotionGG;

        meshToStructureStencils_.clear();
        meshToStructureStencils_.resize(meshMotionGridGeometry.gridView().size(0));
        isCoupledMeshMotionDof_.clear();
        isCoupledMeshMotionDof_.resize(meshMotionGridGeometry.gridView().size(0), false);
        meshToStructureElements_.clear();
        meshToStructureElements_.resize(meshMotionGridGeometry.gridView().size(0));

        structureToMeshMotionStencils_.clear();
        structureToMeshMotionStencils_.resize(structureGridGeometry.gridView().size(0));
        isCoupledStructureDof_.clear();
        isCoupledStructureDof_.resize(structureGridGeometry.gridView().size(0), false);
        structureToMeshMotionElements_.clear();
        structureToMeshMotionElements_.resize(structureGridGeometry.gridView().size(0));

        auto structureFVGeometry = localView(structureGridGeometry);
        for (const auto& element : elements(structureGridGeometry.gridView()))
        {
            structureFVGeometry.bindElement(element);
            const auto eIdx = structureFVGeometry.elementIndex();

            if (!structureFVGeometry.hasBoundaryScvf())
                continue;

            for (const auto& scv : scvs(structureFVGeometry))
            {
                const auto indices = intersectingEntities(scv.dofPosition(), meshMotionGridGeometry.boundingBoxTree());

                // skip if no intersection was found
                if (indices.empty())
                    continue;

                isCoupledStructureDof_[scv.dofIndex()] = true;

                for (const auto idx : indices)
                {
                    structureToMeshMotionElements_[eIdx].push_back(idx);
                    meshToStructureElements_[idx].push_back(eIdx);
                    const auto meshMotionElement = meshMotionGridGeometry.element(idx);

                    const auto meshMotionFVGeometry = localView(meshMotionGridGeometry).bindElement(meshMotionElement);
                    for (const auto& coupledScv : scvs(meshMotionFVGeometry))
                        structureToMeshMotionStencils_[eIdx].push_back(coupledScv.dofIndex());
                }
            }
        }

        auto meshMotionFVGeometry = localView(meshMotionGridGeometry);
        for (const auto& element : elements(meshMotionGridGeometry.gridView()))
        {
            meshMotionFVGeometry.bindElement(element);
            const auto eIdx = meshMotionFVGeometry.elementIndex();

            if (!meshMotionFVGeometry.hasBoundaryScvf())
                continue;

            for (const auto& scv : scvs(meshMotionFVGeometry))
            {
                const auto indices = intersectingEntities(scv.dofPosition(), structureGridGeometry.boundingBoxTree());

                // skip if no intersection was found
                if (indices.empty())
                    continue;

                isCoupledMeshMotionDof_[scv.dofIndex()] = true;

                for (const auto idx : indices)
                {
                    meshToStructureElements_[eIdx].push_back(idx);
                    structureToMeshMotionElements_[idx].push_back(eIdx);
                    const auto structureElement = structureGridGeometry.element(idx);

                    const auto structureFVGeometry = localView(structureGridGeometry).bindElement(structureElement);
                    for (const auto& coupledScv : scvs(structureFVGeometry))
                        meshToStructureStencils_[eIdx].push_back(coupledScv.dofIndex());
                }
            }
        }

        // make stencils unique
        for (std::size_t eIdx = 0; eIdx < meshToStructureStencils_.size(); ++eIdx)
        {
            std::sort(meshToStructureStencils_[eIdx].begin(), meshToStructureStencils_[eIdx].end());
            meshToStructureStencils_[eIdx].erase(
                std::unique(meshToStructureStencils_[eIdx].begin(), meshToStructureStencils_[eIdx].end()),
                meshToStructureStencils_[eIdx].end()
            );

            std::sort(meshToStructureElements_[eIdx].begin(), meshToStructureElements_[eIdx].end());
            meshToStructureElements_[eIdx].erase(
                std::unique(meshToStructureElements_[eIdx].begin(), meshToStructureElements_[eIdx].end()),
                meshToStructureElements_[eIdx].end()
            );
        }

        for (std::size_t eIdx = 0; eIdx < structureToMeshMotionStencils_.size(); ++eIdx)
        {
            std::sort(structureToMeshMotionStencils_[eIdx].begin(), structureToMeshMotionStencils_[eIdx].end());
            structureToMeshMotionStencils_[eIdx].erase(
                std::unique(structureToMeshMotionStencils_[eIdx].begin(), structureToMeshMotionStencils_[eIdx].end()),
                structureToMeshMotionStencils_[eIdx].end()
            );

            std::sort(structureToMeshMotionElements_[eIdx].begin(), structureToMeshMotionElements_[eIdx].end());
            structureToMeshMotionElements_[eIdx].erase(
                std::unique(structureToMeshMotionElements_[eIdx].begin(), structureToMeshMotionElements_[eIdx].end()),
                structureToMeshMotionElements_[eIdx].end()
            );
        }
    }

    void init(std::shared_ptr<Problem<structureIdx>> structureProblem,
              std::shared_ptr<Problem<meshMotionIdx>> meshMotionProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        this->setSubProblems(std::make_tuple(structureProblem, meshMotionProblem));
    }

    template<std::size_t i, std::size_t j>
    const std::vector<std::size_t>& couplingStencil(Dune::index_constant<i> domainI,
                                                    const Element<i>& element,
                                                    Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if constexpr (domainI == structureIdx)
            return emptyStencil_;
        else
            return meshToStructureStencils_[eIdx];
    }

    auto displacement(typename GridGeometry<meshMotionIdx>::LocalView::Element const& element,
                      typename GridGeometry<meshMotionIdx>::SubControlVolume const& scv) const
    {
        if (!isCoupledMeshMotionDof_[scv.dofIndex()])
            DUNE_THROW(Dune::InvalidStateException, "This dof is not coupled!");

        const auto& gg = this->problem(structureIdx).gridGeometry();
        const auto& ggMM = this->problem(meshMotionIdx).gridGeometry();
        const auto meshEIdx = ggMM.elementMapper().index(element);
        if (meshToStructureElements_[meshEIdx].empty())
            DUNE_THROW(Dune::InvalidStateException, "This mesh element (" << meshEIdx << ") is not coupled!");

        const auto structureEIdx = meshToStructureElements_[meshEIdx][0];
        // std::cout << "structureEIdx: " << structureEIdx << std::endl;
        // std::cout << gg.gridView().size(0) << std::endl;
        const auto structureElement = gg.element(structureEIdx);
        const auto elemSol = elementSolution(structureElement, curSol(structureIdx), gg);

        return evalSolution(
            structureElement,
            structureElement.geometry(),
            gg, elemSol,
            scv.dofPosition()
        );
    }

    template<std::size_t id>
    bool isCoupledDof(std::size_t dofIdx, Dune::index_constant<id> domainIdx) const
    {
        if constexpr (id == structureIdx)
            return isCoupledStructureDof_[dofIdx];
        else
            return isCoupledMeshMotionDof_[dofIdx];
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
    { elementSets_ = computeColoring(this->problem(meshMotionIdx).gridGeometry()).sets; }

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
    std::deque<std::vector<ElementSeed<structureIdx>>> elementSets_;

    std::vector<bool> isCoupledStructureDof_, isCoupledMeshMotionDof_;
    std::vector<std::vector<std::size_t>> structureToMeshMotionStencils_, meshToStructureStencils_;
    std::vector<std::vector<std::size_t>> structureToMeshMotionElements_, meshToStructureElements_;

    std::vector<std::size_t> emptyStencil_;
};

// we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<StructureMeshCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
