// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElasticLargeDef
 * \brief Coupling manager for the large-deformation poroelastic mixed
 *        \f$ \mathbf{d} \f$–\f$ p_s \f$–\f$ p_f \f$ model.
 */
#ifndef DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_COUPLING_MANAGER_HH
#define DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_COUPLING_MANAGER_HH

#include <array>
#include <deque>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/assembly/coloring.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/assemblytraits.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/parallel/parallel_for.hh>

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Coupling manager for the large-deformation poroelastic model.
 *
 * All three subdomains share the same grid.
 *
 * \tparam MDTraits MultiDomainTraits with the subdomains momentum (0), solid pressure (1) and fluid pressure (2)
 */
template<class MDTraits>
class PoroElasticLargeDefCouplingManager : public CouplingManager<MDTraits>
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
    static constexpr std::size_t spIdx_ = 1;
    static constexpr std::size_t fpIdx_ = 2;

    using GlobalPosition = typename Element<0>::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element<0>::Geometry::LocalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, GridView<0>::dimension, GridView<0>::dimension>;
    static constexpr int dim = GridView<0>::dimension;

public:
    static constexpr auto momentumIdx = Dune::index_constant<momIdx_>{};
    static constexpr auto solidPressureIdx = Dune::index_constant<spIdx_>{};
    static constexpr auto fluidPressureIdx = Dune::index_constant<fpIdx_>{};

    using CouplingStencilType = std::vector<std::size_t>;

    PoroElasticLargeDefCouplingManager(
        std::shared_ptr<GridDiscretization<momIdx_>> momGridDisc,
        std::shared_ptr<GridDiscretization<spIdx_>> spGridDisc,
        std::shared_ptr<GridDiscretization<fpIdx_>> fpGridDisc)
    {
        buildStencils_(*momGridDisc, *spGridDisc, *fpGridDisc);
    }

    void init(
        std::shared_ptr<GetPropType<SubDomainTypeTag<momIdx_>, Properties::Problem>> momProblem,
        std::shared_ptr<GetPropType<SubDomainTypeTag<spIdx_>, Properties::Problem>> spProblem,
        std::shared_ptr<GetPropType<SubDomainTypeTag<fpIdx_>, Properties::Problem>> fpProblem,
        const SolutionVector& curSol)
    {
        ParentType::updateSolution(curSol);
        this->setSubProblems(std::make_tuple(momProblem, spProblem, fpProblem));
    }

    /*!
     * \brief Evaluate the residual of domain i w.r.t. a perturbation in domain j.
     *
     * This is the full stage residual, since the fluid storage depends on the displacement through \f$J\f$.
     */
    template<std::size_t i, class LocalAssemblerI, std::size_t j>
    typename LocalAssemblerI::LocalResidual::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<i>,
                         const LocalAssemblerI& localAssemblerI,
                         Dune::index_constant<j>, std::size_t)
    {
        const auto [tWeight, sWeight] = localAssemblerI.assembler().currentStageWeights();
        return localAssemblerI.evalLocalResidualForStage(localAssemblerI.curElemVars(), tWeight, sWeight);
    }

    template<std::size_t i, std::size_t j>
    const CouplingStencilType& couplingStencil(Dune::index_constant<i>,
                                               const Element<j>& element,
                                               Dune::index_constant<j>) const
    {
        static_assert(i != j);
        const auto eIdx = this->problem(Dune::index_constant<j>{})
                              .gridDiscretization().elementMapper().index(element);
        return stencils_[j][i][eIdx];
    }

    template<Concept::IpData IpData>
    Scalar solidPressure(const ElementDiscretization<momIdx_>& elemDisc,
                         const IpData& ipData) const
    { return pressure_(elemDisc.element(), ipData, this->curSol(solidPressureIdx)); }

    template<Concept::IpData IpData>
    Scalar fluidPressure(const ElementDiscretization<momIdx_>& elemDisc,
                         const IpData& ipData) const
    { return pressure_(elemDisc.element(), ipData, this->curSol(fluidPressureIdx)); }

    template<Concept::IpData IpData>
    Tensor deformationGradient(const ElementDiscretization<fpIdx_>& elemDisc,
                               const IpData& ipData) const
    {
        const auto& momGridDisc = this->problem(momentumIdx).gridDiscretization();
        const auto& element = elemDisc.element();
        const auto& momSol = this->curSol(momentumIdx);

        // the basis gradients only depend on the reference geometry and are cached per element and
        // integration point; points are compared in local coordinates, so the tolerance scales with
        // the element size
        thread_local std::size_t cachedEIdx = std::numeric_limits<std::size_t>::max();
        thread_local std::vector<IpGradients_> ipCache;

        const auto eIdx = momGridDisc.elementMapper().index(element);
        if (eIdx != cachedEIdx) { ipCache.clear(); cachedEIdx = eIdx; }

        const auto& ipLocal = ipData.local();
        const IpGradients_* grads = nullptr;
        for (const auto& entry : ipCache)
            if ((entry.ipLocal - ipLocal).two_norm2() < 1e-20) { grads = &entry; break; }

        if (grads == nullptr)
        {
            const auto elemDiscMom = localView(momGridDisc).bindElement(element);
            const auto& localBasis = elemDiscMom.feLocalBasis();
            using ShapeJac = typename std::decay_t<decltype(localBasis)>::Traits::JacobianType;
            thread_local std::vector<ShapeJac> shapeJac;
            localBasis.evaluateJacobian(ipLocal, shapeJac);
            const auto jacInvT = element.geometry().jacobianInverseTransposed(ipLocal);

            IpGradients_ entry;
            entry.ipLocal = ipLocal;
            for (const auto& localDof : localDofs(elemDiscMom))
            {
                GlobalPosition gradN;
                jacInvT.mv(shapeJac[localDof.index()][0], gradN);
                entry.dofGrads.push_back({localDof.dofIndex(), gradN});
            }
            ipCache.push_back(std::move(entry));
            grads = &ipCache.back();
        }

        Tensor F(0.0);
        for (const auto& dg : grads->dofGrads)
            for (int d = 0; d < dim; ++d)
                F[d].axpy(momSol[dg.dofIndex][d], dg.gradN);
        for (int d = 0; d < dim; ++d) F[d][d] += 1.0;
        return F;
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
    const SolutionVector& curSol() const { return ParentType::curSol(); }

private:
    struct DofGrad_ { std::size_t dofIndex; GlobalPosition gradN; };
    static constexpr std::size_t maxMomLocalDofs_ = 32;
    struct IpGradients_
    {
        LocalPosition ipLocal;
        Dune::ReservedVector<DofGrad_, maxMomLocalDofs_> dofGrads;
    };

    struct DofShape_ { std::size_t dofIndex; Scalar shapeValue; };
    static constexpr std::size_t maxPresLocalDofs_ = 8;
    struct IpShapes_
    {
        LocalPosition ipLocal;
        Dune::ReservedVector<DofShape_, maxPresLocalDofs_> dofShapes;
    };

    // solid and fluid pressure are both interpolated with the basis of the fluid pressure
    // discretization, so both need the same discretization
    template<class IpData, class SolVec>
    Scalar pressure_(const Element<fpIdx_>& element,
                     const IpData& ipData,
                     const SolVec& sol) const
    {
        thread_local std::size_t cachedEIdx = std::numeric_limits<std::size_t>::max();
        thread_local std::vector<IpShapes_> ipCache;

        const auto& presGridDisc = this->problem(fluidPressureIdx).gridDiscretization();
        const auto eIdx = presGridDisc.elementMapper().index(element);
        if (eIdx != cachedEIdx) { ipCache.clear(); cachedEIdx = eIdx; }

        const auto& ipLocal = ipData.local();
        const IpShapes_* shapes = nullptr;
        for (const auto& entry : ipCache)
            if ((entry.ipLocal - ipLocal).two_norm2() < 1e-20) { shapes = &entry; break; }

        if (shapes == nullptr)
        {
            const auto elemDiscPres = localView(presGridDisc).bindElement(element);
            const auto& localBasis = elemDiscPres.feLocalBasis();
            using RangeType = typename std::decay_t<decltype(localBasis)>::Traits::RangeType;
            thread_local std::vector<RangeType> shapeVals;
            localBasis.evaluateFunction(ipLocal, shapeVals);

            IpShapes_ entry;
            entry.ipLocal = ipLocal;
            for (const auto& localDof : localDofs(elemDiscPres))
                entry.dofShapes.push_back({localDof.dofIndex(), shapeVals[localDof.index()][0]});
            ipCache.push_back(std::move(entry));
            shapes = &ipCache.back();
        }

        Scalar p = 0.0;
        for (const auto& ds : shapes->dofShapes)
            p += sol[ds.dofIndex][0] * ds.shapeValue;
        return p;
    }

    void buildStencils_(const GridDiscretization<momIdx_>& momGridDisc,
                        const GridDiscretization<spIdx_>& spGridDisc,
                        const GridDiscretization<fpIdx_>& fpGridDisc)
    {
        const std::size_t n = fpGridDisc.gridView().size(0);
        for (auto& arr : stencils_)
            for (auto& s : arr)
            {
                s.clear();
                s.resize(n);
            }

        auto elemDiscMom = localView(momGridDisc);
        auto elemDiscSp = localView(spGridDisc);

        for (const auto& element : elements(fpGridDisc.gridView()))
        {
            elemDiscMom.bindElement(element);
            elemDiscSp.bindElement(element);
            const std::size_t eIdx = fpGridDisc.elementMapper().index(element);

            for (const auto& localDof : localDofs(elemDiscMom))
            {
                stencils_[momIdx_][spIdx_][eIdx].push_back(localDof.dofIndex());
                stencils_[momIdx_][fpIdx_][eIdx].push_back(localDof.dofIndex());
            }
            for (const auto& localDof : localDofs(elemDiscSp))
            {
                stencils_[spIdx_][momIdx_][eIdx].push_back(localDof.dofIndex());
                stencils_[fpIdx_][momIdx_][eIdx].push_back(localDof.dofIndex());
            }
        }
    }

    //! stencils_[j][i][eIdx]: the dofs of domain j that the residual of element eIdx of domain i depends on
    std::array<std::array<std::vector<CouplingStencilType>, 3>, 3> stencils_;

    std::deque<std::vector<ElementSeed<momIdx_>>> elementSets_;
};

template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<PoroElasticLargeDefCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux
#endif
