// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HyperelasticVolIso
 * \brief Coupling manager for the mixed u–p hyperelastic model.
 *
 * Both subdomains (momentum and pressure) share the same grid.
 *
 * Services provided:
 * - `pressure(elemDisc, ipData)`:
 *   evaluates \f$ p_s \f$ at an interpolation point from the stored pressure solution.
 *   Called by the momentum residual.
 * - `deformationGradient(elemDisc, ipData)`:
 *   evaluates \f$ \mathbf{F} \f$ at an interpolation point from the stored displacement solution.
 *   Called by the pressure residual.
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_COUPLING_MANAGER_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_COUPLING_MANAGER_HH

#include <deque>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>

#include <dumux/assembly/coloring.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/assemblytraits.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/parallel/parallel_for.hh>

namespace Dumux {

/*!
 * \ingroup HyperelasticVolIso
 * \brief Coupling manager for the mixed u–p hyperelastic model.
 *
 * \tparam MDTraits MultiDomainTraits with subdomain 0 = momentum, subdomain 1 = pressure.
 */
template<class MDTraits>
class HyperelasticVolIsoCouplingManager : public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridDiscretization = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using ElementDiscretization = typename GridDiscretization<id>::LocalView;
    template<std::size_t id> using GridView = typename GridDiscretization<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    static constexpr std::size_t momIdx_ = 0;
    static constexpr std::size_t presIdx_ = 1;

    using GlobalPosition = typename Element<0>::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element<0>::Geometry::LocalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, GridView<0>::dimension, GridView<0>::dimension>;
    static constexpr int dim = GridView<0>::dimension;

public:
    static constexpr auto momentumIdx = Dune::index_constant<momIdx_>();
    static constexpr auto pressureIdx = Dune::index_constant<presIdx_>();

    using CouplingStencilType = std::vector<std::size_t>;

    HyperelasticVolIsoCouplingManager(
        std::shared_ptr<GridDiscretization<momIdx_>> momGridDisc,
        std::shared_ptr<GridDiscretization<presIdx_>> presGridDisc)
    {
        buildStencils_(*momGridDisc, *presGridDisc);
    }

    void init(std::shared_ptr<GetPropType<SubDomainTypeTag<momIdx_>, Properties::Problem>> momProblem,
              std::shared_ptr<GetPropType<SubDomainTypeTag<presIdx_>, Properties::Problem>> presProblem,
              const SolutionVector& curSol)
    {
        ParentType::updateSolution(curSol);
        this->setSubProblems(std::make_tuple(momProblem, presProblem));
    }

    template<std::size_t i, std::size_t j>
    const CouplingStencilType& couplingStencil(Dune::index_constant<i>,
                                               const Element<j>& element,
                                               Dune::index_constant<j>) const
    {
        static_assert(i != j);
        const auto eIdx = this->problem(Dune::index_constant<j>{})
                              .gridDiscretization().elementMapper().index(element);
        if constexpr (i == momIdx_) return stencilsMomToPres_[eIdx];
        else                        return stencilsPresTOMom_[eIdx];
    }

    //! Evaluate \f$ p_s \f$ at an interpolation point (P1 interpolation from pressure solution).
    template<Concept::IpData IpData>
    Scalar pressure(const ElementDiscretization<momIdx_>& elemDisc,
                    const IpData& ipData) const
    {
        const auto& gridDisc = this->problem(pressureIdx).gridDiscretization();
        const auto elemSol = elementSolution(elemDisc.element(), this->curSol(pressureIdx), gridDisc);
        return evalSolutionAtLocalPos(elemDisc.element(), elemDisc.element().geometry(),
                                      gridDisc, elemSol, ipData.local())[0];
    }

    //! Evaluate \f$ \mathbf{F} \f$ at an interpolation point (displacement shape functions).
    template<Concept::IpData IpData>
    Tensor deformationGradient(const ElementDiscretization<presIdx_>& elemDisc,
                               const IpData& ipData) const
    {
        const auto& gridDisc = this->problem(momentumIdx).gridDiscretization();
        const auto elemDiscMom = localView(gridDisc).bindElement(elemDisc.element());
        return computeF_(elemDisc.element(), ipData.local(), elemDiscMom);
    }

    void computeColorsForAssembly()
    {
        elementSets_ = computeColoring(
            this->problem(momentumIdx).gridDiscretization()).sets;
    }

    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i>, AssembleElementFunc&& f) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException,
                "Call computeColorsForAssembly before parallel assembly!");
        const auto& grid = this->problem(Dune::index_constant<i>{})
                               .gridDiscretization().gridView().grid();
        for (const auto& elems : elementSets_)
            Dumux::parallelFor(elems.size(), [&](std::size_t n)
            { f(grid.entity(elems[n])); });
    }

    template<std::size_t i>
    auto& curSol(Dune::index_constant<i> id) { return ParentType::curSol(id); }
    template<std::size_t i>
    const auto& curSol(Dune::index_constant<i> id) const { return ParentType::curSol(id); }

private:
    Tensor computeF_(const Element<presIdx_>& element,
                     const LocalPosition& ipLocal,
                     const ElementDiscretization<momIdx_>& elemDiscMom) const
    {
        const auto& geom = element.geometry();
        const auto& localBasis = elemDiscMom.feLocalBasis();
        using ShapeJac = typename std::decay_t<decltype(localBasis)>::Traits::JacobianType;
        std::vector<ShapeJac> shapeJac;
        localBasis.evaluateJacobian(ipLocal, shapeJac);
        const auto jacInvT = geom.jacobianInverseTransposed(ipLocal);

        Tensor F(0.0);
        for (const auto& localDof : localDofs(elemDiscMom))
        {
            GlobalPosition gradN;
            jacInvT.mv(shapeJac[localDof.index()][0], gradN);
            for (int d = 0; d < dim; ++d)
                F[d].axpy(this->curSol(momentumIdx)[localDof.dofIndex()][d], gradN);
        }
        for (int d = 0; d < dim; ++d) F[d][d] += 1.0;
        return F;
    }

    void buildStencils_(const GridDiscretization<momIdx_>& momGridDisc,
                        const GridDiscretization<presIdx_>& presGridDisc)
    {
        const std::size_t n = presGridDisc.gridView().size(0);
        stencilsMomToPres_.resize(n);
        stencilsPresTOMom_.resize(n);

        auto elemDiscMom = localView(momGridDisc);
        auto elemDiscPres = localView(presGridDisc);

        for (const auto& element : elements(presGridDisc.gridView()))
        {
            elemDiscMom.bindElement(element);
            elemDiscPres.bindElement(element);
            const std::size_t eIdx = presGridDisc.elementMapper().index(element);

            for (const auto& localDof : localDofs(elemDiscPres))
                stencilsMomToPres_[eIdx].push_back(localDof.dofIndex());
            // Use localDofs (not scvs) to include non-SCV DOFs (e.g. PQ2 edge midpoints)
            for (const auto& localDof : localDofs(elemDiscMom))
                stencilsPresTOMom_[eIdx].push_back(localDof.dofIndex());
        }
    }

    std::deque<std::vector<ElementSeed<momIdx_>>> elementSets_;
    std::vector<CouplingStencilType> stencilsMomToPres_, stencilsPresTOMom_;
};

template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<HyperelasticVolIsoCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux
#endif
