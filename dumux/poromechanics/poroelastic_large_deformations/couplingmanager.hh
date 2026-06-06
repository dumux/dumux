// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElasticLargeDef
 * \brief Coupling manager for the large-deformation poroelastic mixed u–p_s–p_f model.
 *
 * All three subdomains share the same grid.
 *
 * Services provided:
 * - `solidPressureAtPoint(fvGeom, ip)`:   interpolates \f$p_s\f$ at a global point.
 * - `fluidPressureAtPoint(fvGeom, ip)`:   interpolates \f$p_f\f$ at a global point.
 * - `deformationGradientAtPoint(fvGeom, ip)`:
 *     evaluates \f$\mathbf{F}\f$ from the current displacement field.
 *     The assembler keeps `curSol` up-to-date (via `updateGridVariables`) so
 *     this always reflects the correct Newton iterate or stage state.
 */
#ifndef DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_COUPLING_MANAGER_HH
#define DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_COUPLING_MANAGER_HH

#include <vector>
#include <limits>
#include <dune/common/fmatrix.hh>
#include <dune/common/reservedvector.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Coupling manager for the large-deformation poroelastic model.
 *
 * \tparam MDTraits MultiDomainTraits with
 *   - subdomain 0 = momentum     (displacement)
 *   - subdomain 1 = solid pressure
 *   - subdomain 2 = fluid pressure
 */
template<class MDTraits>
class PoroElasticLargeDefCouplingManager : public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    static constexpr std::size_t momIdx_  = 0;
    static constexpr std::size_t spIdx_   = 1;
    static constexpr std::size_t fpIdx_   = 2;

    using GlobalPosition = typename Element<0>::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, GridView<0>::dimension, GridView<0>::dimension>;
    static constexpr int dim = GridView<0>::dimension;

public:
    static constexpr auto momentumIdx      = Dune::index_constant<momIdx_>{};
    static constexpr auto solidPressureIdx = Dune::index_constant<spIdx_>{};
    static constexpr auto fluidPressureIdx = Dune::index_constant<fpIdx_>{};

    using CouplingStencilType = std::vector<std::size_t>;

    PoroElasticLargeDefCouplingManager(
        std::shared_ptr<GridGeometry<momIdx_>> momGG,
        std::shared_ptr<GridGeometry<spIdx_>>  spGG,
        std::shared_ptr<GridGeometry<fpIdx_>>  fpGG)
    {
        buildStencils_(*momGG, *spGG, *fpGG);
    }

    void init(
        std::shared_ptr<GetPropType<SubDomainTypeTag<momIdx_>, Properties::Problem>> momProblem,
        std::shared_ptr<GetPropType<SubDomainTypeTag<spIdx_>,  Properties::Problem>>  spProblem,
        std::shared_ptr<GetPropType<SubDomainTypeTag<fpIdx_>,  Properties::Problem>>  fpProblem,
        const SolutionVector& curSol)
    {
        ParentType::updateSolution(curSol);
        this->setSubProblems(std::make_tuple(momProblem, spProblem, fpProblem));
    }

    /*!
     * \brief Evaluate the residual of domain i w.r.t. a perturbation in domain j.
     *
     * Uses the full stage-weighted residual (storage + spatial) so that the cross-domain
     * Jacobian blocks include the temporal coupling — in particular, dJ/dd for the fluid
     * pressure storage term.  The stage weights are retrieved from the assembler view so
     * this works transparently with any multi-stage method.
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
                              .gridGeometry().elementMapper().index(element);
        return stencils_[j][i][eIdx];
    }

    //! Evaluate \f$p_s\f$ at a global position (P1 interpolation from solid pressure solution).
    Scalar solidPressureAtPoint(const typename GridGeometry<momIdx_>::LocalView& fvGeom,
                                const GlobalPosition& ip) const
    { return pressureAtPoint_(fvGeom.element(), ip, this->curSol(solidPressureIdx)); }

    //! Evaluate \f$p_f\f$ at a global position (P1 interpolation from fluid pressure solution).
    Scalar fluidPressureAtPoint(const typename GridGeometry<momIdx_>::LocalView& fvGeom,
                                const GlobalPosition& ip) const
    { return pressureAtPoint_(fvGeom.element(), ip, this->curSol(fluidPressureIdx)); }

    /*!
     * \brief Evaluate the deformation gradient \f$\mathbf{F}\f$ at a global position.
     *
     * Always uses the coupling manager's current solution (`curSol`), which the
     * assembler keeps up-to-date via `updateGridVariables` at each Newton iteration
     * and at each multi-stage evaluation.
     *
     * \param fvGeom  Local FV geometry view bound to the element of interest
     * \param ip      Global position at which to evaluate \f$\mathbf{F}\f$
     */
    Tensor deformationGradientAtPoint(
        const typename GridGeometry<fpIdx_>::LocalView& fvGeom,
        const GlobalPosition& ip) const
    {
        const auto& momGG = this->problem(momentumIdx).gridGeometry();
        const auto& element = fvGeom.element();
        const auto& momSol = this->curSol(momentumIdx);

        // Element-wise cache of the reference-configuration basis gradients gradN
        // (= grad_X N) per momentum local dof, keyed by the queried integration
        // point. gradN depends only on the (undeformed) element geometry, so it is
        // constant for an element across the whole solve — only the displacement
        // *values* change when the numeric-differentiation deflects a dof. Caching
        // it removes the dominant cost (the local-basis Jacobian evaluation, which
        // was otherwise recomputed on every perturbation). thread_local keeps this
        // correct under the coloured multithreaded assembly (each thread processes
        // an element to completion before moving on).
        thread_local std::size_t cachedEIdx = std::numeric_limits<std::size_t>::max();
        thread_local std::vector<IpGradients_> ipCache;

        const auto eIdx = momGG.elementMapper().index(element);
        if (eIdx != cachedEIdx) { ipCache.clear(); cachedEIdx = eIdx; }

        const IpGradients_* grads = nullptr;
        for (const auto& entry : ipCache)
            if ((entry.ip - ip).two_norm2() < 1e-20) { grads = &entry; break; }

        if (grads == nullptr)
        {
            const auto fvGeomMom = localView(momGG).bindElement(element);
            const auto& geom = element.geometry();
            const auto ipLocal = geom.local(ip);
            const auto& localBasis = fvGeomMom.feLocalBasis();
            using ShapeJac = typename std::decay_t<decltype(localBasis)>::Traits::JacobianType;
            thread_local std::vector<ShapeJac> shapeJac;
            localBasis.evaluateJacobian(ipLocal, shapeJac);
            const auto jacInvT = geom.jacobianInverseTransposed(ipLocal);

            IpGradients_ entry;
            entry.ip = ip;
            for (const auto& localDof : localDofs(fvGeomMom))
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
            this->problem(momentumIdx).gridGeometry()).sets;
    }

    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i>, AssembleElementFunc&& f) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException,
                "Call computeColorsForAssembly before parallel assembly!");
        const auto& grid = this->problem(Dune::index_constant<i>{})
                               .gridGeometry().gridView().grid();
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
    //! Cached reference-config basis gradient of one momentum local dof.
    struct DofGrad_ { std::size_t dofIndex; GlobalPosition gradN; };
    //! All momentum dof gradients at one integration point (key = ip).
    static constexpr std::size_t maxMomLocalDofs_ = 32;
    struct IpGradients_
    {
        GlobalPosition ip;
        Dune::ReservedVector<DofGrad_, maxMomLocalDofs_> dofGrads;
    };

    //! Cached P1 (Box) shape value of one pressure local dof.
    struct DofShape_ { std::size_t dofIndex; Scalar shapeValue; };
    static constexpr std::size_t maxPresLocalDofs_ = 8;
    struct IpShapes_
    {
        GlobalPosition ip;
        Dune::ReservedVector<DofShape_, maxPresLocalDofs_> dofShapes;
    };

    /*!
     * \brief Interpolate a P1 (Box) pressure field at a global point.
     *
     * The solid- and fluid-pressure subdomains are both Box/P1 on the same grid, so
     * they share the same vertex dof indices and shape functions. The shape values
     * \f$N_i(\mathbf{x})\f$ are geometry-only, so they are cached per (element, ip) —
     * exactly like the deformation-gradient basis gradients — and only the (cheap)
     * dot product with the current solution is recomputed. This keeps the numeric
     * differentiation correct (the perturbed dof value flows through `sol`) while
     * avoiding the repeated shape-function evaluation in `evalSolution`.
     */
    template<class SolVec>
    Scalar pressureAtPoint_(const Element<fpIdx_>& element,
                            const GlobalPosition& ip,
                            const SolVec& sol) const
    {
        thread_local std::size_t cachedEIdx = std::numeric_limits<std::size_t>::max();
        thread_local std::vector<IpShapes_> ipCache;

        const auto& presGG = this->problem(fluidPressureIdx).gridGeometry();
        const auto eIdx = presGG.elementMapper().index(element);
        if (eIdx != cachedEIdx) { ipCache.clear(); cachedEIdx = eIdx; }

        const IpShapes_* shapes = nullptr;
        for (const auto& entry : ipCache)
            if ((entry.ip - ip).two_norm2() < 1e-20) { shapes = &entry; break; }

        if (shapes == nullptr)
        {
            const auto fvGeomPres = localView(presGG).bindElement(element);
            const auto& geom = element.geometry();
            const auto ipLocal = geom.local(ip);
            const auto& localBasis = fvGeomPres.feLocalBasis();
            using RangeType = typename std::decay_t<decltype(localBasis)>::Traits::RangeType;
            thread_local std::vector<RangeType> shapeVals;
            localBasis.evaluateFunction(ipLocal, shapeVals);

            IpShapes_ entry;
            entry.ip = ip;
            for (const auto& localDof : localDofs(fvGeomPres))
                entry.dofShapes.push_back({localDof.dofIndex(), shapeVals[localDof.index()][0]});
            ipCache.push_back(std::move(entry));
            shapes = &ipCache.back();
        }

        Scalar p = 0.0;
        for (const auto& ds : shapes->dofShapes)
            p += sol[ds.dofIndex][0] * ds.shapeValue;
        return p;
    }

    void buildStencils_(const GridGeometry<momIdx_>& momGG,
                        const GridGeometry<spIdx_>&  spGG,
                        const GridGeometry<fpIdx_>&  fpGG)
    {
        // stencils_[domain j][domain i][element index] = DOF indices in domain i
        // that influence element in domain j
        for (auto& arr : stencils_) for (auto& s : arr) s.clear();

        const std::size_t n = fpGG.gridView().size(0);
        for (auto& arr : stencils_) for (auto& s : arr) s.resize(n);

        auto fvMom  = localView(momGG);
        auto fvSp   = localView(spGG);
        auto fvFp   = localView(fpGG);

        for (const auto& element : elements(fpGG.gridView()))
        {
            fvMom.bindElement(element);
            fvSp.bindElement(element);
            fvFp.bindElement(element);
            const std::size_t eIdx = fpGG.elementMapper().index(element);

            // momentum ← solid pressure (p_s appears in stress via p_s constraint)
            for (const auto& localDof : localDofs(fvMom))
                stencils_[momIdx_][spIdx_][eIdx].push_back(localDof.dofIndex());
            // momentum ← fluid pressure
            for (const auto& localDof : localDofs(fvMom))
                stencils_[momIdx_][fpIdx_][eIdx].push_back(localDof.dofIndex());
            // solid pressure ← momentum (p_s^eq depends on J)
            for (const auto& localDof : localDofs(fvSp))
                stencils_[spIdx_][momIdx_][eIdx].push_back(localDof.dofIndex());
            // fluid pressure ← momentum (J in storage and flux)
            for (const auto& localDof : localDofs(fvFp))
                stencils_[fpIdx_][momIdx_][eIdx].push_back(localDof.dofIndex());
        }
    }

    // stencils_[domainI][domainJ][elementIndex] = DOF indices in J affecting element in I
    std::array<std::array<std::vector<CouplingStencilType>, 3>, 3> stencils_;

    std::deque<std::vector<ElementSeed<momIdx_>>> elementSets_;
};

template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<PoroElasticLargeDefCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux
#endif
