// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FoepplVonKarmanPlate
 * \brief Coupling manager for the Föppl-von Kármán model
 */
#ifndef DUMUX_FOEPPL_VON_KARMAN_PLATE_COUPLINGMANAGER_HH
#define DUMUX_FOEPPL_VON_KARMAN_PLATE_COUPLINGMANAGER_HH

#include <memory>
#include <tuple>
#include <vector>
#include <deque>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Coupling manager for the three sub-problems of the Föppl-von Kármán plate
 *
 * All three sub-problems are discretized on the same grid, so the coupling is
 * element-local: the coupling stencil of an element is the set of local degrees
 * of freedom of the coupled sub-problem on that same element.
 *
 * The coupling is directed: the rotations need the potentials \f$ (\varphi, \psi) \f$,
 * the deformation needs the rotations \f$ \boldsymbol{\theta} \f$ and the in-plane
 * displacement gradient \f$ \nabla\mathbf{u} \f$, and the in-plane displacements
 * need the slope \f$ \nabla w \f$. Rotations and in-plane displacements do not
 * appear in each other's equations.
 */
template<class MDTraits>
class FoepplVonKarmanPlateCouplingManager
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

public:
    static constexpr auto rotationIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto deformationIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto inPlaneIdx = typename MDTraits::template SubDomain<2>::Index();

    //! export traits
    using MultiDomainTraits = MDTraits;

    FoepplVonKarmanPlateCouplingManager() = default;

    FoepplVonKarmanPlateCouplingManager(std::shared_ptr<GridGeometry<rotationIdx>> rotationGG,
                                        std::shared_ptr<GridGeometry<deformationIdx>> deformationGG,
                                        std::shared_ptr<GridGeometry<inPlaneIdx>> inPlaneGG)
    {
        const auto numElements = rotationGG->gridView().size(0);
        std::get<rotationIdx>(stencils_).assign(numElements, CouplingStencil<rotationIdx>{});
        std::get<deformationIdx>(stencils_).assign(numElements, CouplingStencil<deformationIdx>{});
        std::get<inPlaneIdx>(stencils_).assign(numElements, CouplingStencil<inPlaneIdx>{});

        auto rotGeo = localView(*rotationGG);
        auto defGeo = localView(*deformationGG);
        auto inPlaneGeo = localView(*inPlaneGG);

        for (const auto& element : elements(rotationGG->gridView()))
        {
            rotGeo.bindElement(element);
            defGeo.bindElement(element);
            inPlaneGeo.bindElement(element);
            const auto eIdx = rotGeo.elementIndex();

            for (const auto& localDof : localDofs(rotGeo))
                std::get<rotationIdx>(stencils_)[eIdx].push_back(localDof.dofIndex());

            for (const auto& localDof : localDofs(defGeo))
                std::get<deformationIdx>(stencils_)[eIdx].push_back(localDof.dofIndex());

            for (const auto& localDof : localDofs(inPlaneGeo))
                std::get<inPlaneIdx>(stencils_)[eIdx].push_back(localDof.dofIndex());
        }
    }

    void init(std::shared_ptr<Problem<rotationIdx>> rotationProblem,
              std::shared_ptr<Problem<deformationIdx>> deformationProblem,
              std::shared_ptr<Problem<inPlaneIdx>> inPlaneProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        this->setSubProblems(std::make_tuple(rotationProblem, deformationProblem, inPlaneProblem));
    }

    template<std::size_t i, std::size_t j>
    const CouplingStencil<j>& couplingStencil(Dune::index_constant<i> domainI,
                                              const Element<i>& element,
                                              Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        // the moment equilibrium does not contain the in-plane state and the
        // in-plane equilibrium does not contain the rotations
        if constexpr ((i == rotationIdx && j == inPlaneIdx) || (i == inPlaneIdx && j == rotationIdx))
            return std::get<j>(emptyStencils_);
        else
        {
            const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
            return std::get<j>(stencils_)[eIdx];
        }
    }

    //! The rotation vector \f$ \boldsymbol{\theta} \f$ at a point of an element
    template<class FVElementGeometry, class GlobalPosition>
    auto rotation(const FVElementGeometry& fvGeometry, const GlobalPosition& globalPos) const
    { return evalAt_(rotationIdx, fvGeometry, position_(globalPos)); }

    //! The deformation and potentials \f$ (\varphi, w, \psi) \f$ at a point of an element
    template<class FVElementGeometry, class GlobalPosition>
    auto deformationAndPotentials(const FVElementGeometry& fvGeometry, const GlobalPosition& globalPos) const
    { return evalAt_(deformationIdx, fvGeometry, position_(globalPos)); }

    //! The in-plane displacement gradient \f$ \nabla\mathbf{u} \f$ at a point of an element
    template<class FVElementGeometry, class GlobalPosition>
    auto inPlaneDisplacementGradient(const FVElementGeometry& fvGeometry, const GlobalPosition& globalPos) const
    { return evalGradientsAt_(inPlaneIdx, fvGeometry, position_(globalPos)); }

    //! The slope \f$ \nabla w \f$ at a point of an element
    template<class FVElementGeometry, class GlobalPosition>
    auto verticalDeformationGradient(const FVElementGeometry& fvGeometry, const GlobalPosition& globalPos) const
    { return evalGradientsAt_(deformationIdx, fvGeometry, position_(globalPos))[Indices<deformationIdx>::verticalDeformationIdx]; }

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

private:
    //! Accept either a position or anything that carries one, so both a sub control volume
    //! face and a quadrature point can name where the coupled state is read
    template<class Position>
    static decltype(auto) position_(const Position& position)
    {
        if constexpr (requires { position.ipGlobal(); })
            return position.ipGlobal();
        else if constexpr (requires { position.global(); })
            return position.global();
        else
            return (position);
    }

    template<std::size_t i, class FVElementGeometry, class GlobalPosition>
    auto evalAt_(Dune::index_constant<i> domainIdx,
                 const FVElementGeometry& fvGeometry,
                 const GlobalPosition& globalPos) const
    {
        const auto& gg = this->problem(domainIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), ParentType::curSol(domainIdx), gg);
        return evalSolution(
            fvGeometry.element(), fvGeometry.element().geometry(), gg, elemSol, globalPos
        );
    }

    template<std::size_t i, class FVElementGeometry, class GlobalPosition>
    auto evalGradientsAt_(Dune::index_constant<i> domainIdx,
                          const FVElementGeometry& fvGeometry,
                          const GlobalPosition& globalPos) const
    {
        const auto& gg = this->problem(domainIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), ParentType::curSol(domainIdx), gg);
        return evalGradients(
            fvGeometry.element(), fvGeometry.element().geometry(), gg, elemSol, globalPos
        );
    }

    //! coloring for multithreaded assembly
    std::deque<std::vector<ElementSeed<deformationIdx>>> elementSets_;

    //! for each domain, the local degrees of freedom on each element
    std::tuple<std::vector<CouplingStencil<rotationIdx>>,
               std::vector<CouplingStencil<deformationIdx>>,
               std::vector<CouplingStencil<inPlaneIdx>>> stencils_;

    std::tuple<CouplingStencil<rotationIdx>,
               CouplingStencil<deformationIdx>,
               CouplingStencil<inPlaneIdx>> emptyStencils_;
};

template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<FoepplVonKarmanPlateCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
