// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Assembly
 * \brief A local operator wrapper for multi-stage time stepping schemes
 */
#ifndef DUMUX_EXPERIMENTAL_MULTISTAGE_FV_LOCAL_OPERATOR_HH
#define DUMUX_EXPERIMENTAL_MULTISTAGE_FV_LOCAL_OPERATOR_HH

#include <cmath>
#include <dumux/discretization/extrusion.hh>

namespace Dumux::Experimental {

namespace Detail {

//! Lightweight proxies that scale every scalar `+=` written through
//! `proxy[i][j][r][c] += d` by a fixed weight `w` (so the underlying matrix
//! receives `+= w*d`). Used to apply a multistage spatial/temporal weight to
//! hand-coded analytic Jacobian contributions without modifying the (validated)
//! derivative math, which writes directly into the global matrix.
template<class Scalar>
class ScaledScalarRef
{
public:
    ScaledScalarRef(Scalar& ref, double w) : ref_(ref), w_(w) {}
    template<class T>
    ScaledScalarRef& operator+=(const T& d) { ref_ += w_*d; return *this; }
    template<class T>
    ScaledScalarRef& operator-=(const T& d) { ref_ -= w_*d; return *this; }
private:
    Scalar& ref_;
    double w_;
};

template<class BlockRow>
class ScaledBlockRowProxy
{
public:
    ScaledBlockRowProxy(BlockRow& row, double w) : row_(row), w_(w) {}
    auto operator[](std::size_t c) { return ScaledScalarRef<std::decay_t<decltype(row_[c])>>(row_[c], w_); }
private:
    BlockRow& row_;
    double w_;
};

template<class Block>
class ScaledBlockProxy
{
public:
    ScaledBlockProxy(Block& block, double w) : block_(block), w_(w) {}
    auto operator[](std::size_t r) { return ScaledBlockRowProxy<std::decay_t<decltype(block_[r])>>(block_[r], w_); }
private:
    Block& block_;
    double w_;
};

template<class Row>
class ScaledRowProxy
{
public:
    ScaledRowProxy(Row& row, double w) : row_(row), w_(w) {}
    auto operator[](std::size_t j) { return ScaledBlockProxy<std::decay_t<decltype(row_[j])>>(row_[j], w_); }
private:
    Row& row_;
    double w_;
};

template<class Matrix>
class ScaledMatrixProxy
{
public:
    ScaledMatrixProxy(Matrix& A, double w) : A_(A), w_(w) {}
    auto operator[](std::size_t i) { return ScaledRowProxy<std::decay_t<decltype(A_[i])>>(A_[i], w_); }
private:
    Matrix& A_;
    double w_;
};

} // end namespace Detail

template<class LocalOperator>
class MultiStageFVLocalOperator
{
    using ElementOperatorResultVector = typename LocalOperator::ElementResidualVector;
public:
    // compatibility
    using ElementResidualVector = ElementOperatorResultVector;

    MultiStageFVLocalOperator(const LocalOperator& op)
    : op_(op)
    , spatialWeight_(1.0)
    , temporalWeight_(1.0)
    {}

    // discretization-agnostic interface (apart from FV)
    template<class FVGeometry, class ElemVolVars>
    ElementOperatorResultVector evalStorage(
        const FVGeometry& fvGeometry,
        const ElemVolVars& elemVolVars
    ) const {
        ElementOperatorResultVector result(fvGeometry.numScv());

        if (std::abs(temporalWeight_) > 1e-6)
        {
            for (const auto& scv : scvs(fvGeometry))
                result[scv.localDofIndex()] +=
                    op_.computeStorage(op_.problem(), scv, elemVolVars[scv])
                    * elemVolVars[scv].extrusionFactor()
                    * Extrusion_t<typename FVGeometry::GridGeometry>::volume(fvGeometry, scv)
                    * temporalWeight_;
        }

        return result;
    }

    // discretization-agnostic interface (apart from FV)
    template<class FVGeometry, class ElemVolVars, class ElemFluxVars, class ElemBCTypes>
    ElementOperatorResultVector evalFluxAndSource(
        const typename FVGeometry::Element&, // not needed, here for compatibility
        const FVGeometry& fvGeometry,
        const ElemVolVars& elemVolVars,
        const ElemFluxVars& elemFluxVarsCache,
        const ElemBCTypes& bcTypes
    ) const {
        ElementOperatorResultVector result(fvGeometry.numScv());
        if (std::abs(spatialWeight_) > 1e-6)
        {
            result = op_.evalFluxAndSource(fvGeometry.element(), fvGeometry, elemVolVars, elemFluxVarsCache, bcTypes);
            for (auto& r : result)
                r *= spatialWeight_;
        }

        return result;
    }

    // interface allowing for optimization when computing the cell-centered finite volume Jacobian
    template<class Problem, class FVGeometry, class ElemVolVars, class ElemFluxVars>
    auto evalFlux(
        const Problem&, // not needed
        const typename FVGeometry::Element& element, // can be neighbor
        const FVGeometry& fvGeometry,
        const ElemVolVars& elemVolVars,
        const ElemFluxVars& elemFluxVarsCache,
        const typename FVGeometry::SubControlVolumeFace& scvf
    ) const {
        using NumEqVector = std::decay_t<decltype(op_.evalFlux(op_.problem(), element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf))>;
        NumEqVector result(0.0);
        if (std::abs(spatialWeight_) > 1e-6)
        {
            result = op_.evalFlux(op_.problem(), element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
            result *= spatialWeight_;
        }
        return result;
    }

    void spatialWeight(double w) { spatialWeight_ = w; }
    double spatialWeight() const { return spatialWeight_; }

    void temporalWeight(double w) { temporalWeight_ = w; }
    double temporalWeight() const { return temporalWeight_; }

    const auto& problem() const
    { return op_.problem(); }

    /*!
     * \name Analytic-Jacobian forwarders applying the multistage stage weights.
     *
     * The per-stage weighted residual is  tWeight*storage + sWeight*(flux+source),
     * so the assembled Jacobian block is
     *   A = tWeight * d(volume*computeStorage)/du  +  sWeight * [ d(flux)/du + d(source)/du ].
     * The wrapped operator provides RAW (un-time-scaled) storage derivatives and the
     * (un-scaled) flux/source derivatives; here we attach the temporal/spatial weights.
     */
    // \{

    //! storage block (diagonal): += temporalWeight * volume * d(computeStorage)/du
    template<class PartialDerivativeMatrix, class Problem, class Element, class FVGeometry, class VolVars, class SubControlVolume>
    void addStorageDerivatives(PartialDerivativeMatrix& A_ii,
                               const Problem& problem,
                               const Element& element,
                               const FVGeometry& fvGeometry,
                               const VolVars& curVolVars,
                               const SubControlVolume& scv) const
    {
        if (std::abs(temporalWeight_) < 1e-6)
            return;
        std::decay_t<PartialDerivativeMatrix> raw(0.0);
        op_.addStorageDerivativesRaw(raw, problem, element, fvGeometry, curVolVars, scv);
        A_ii.axpy(temporalWeight_, raw);
    }

    //! source block (diagonal): += spatialWeight * d(source)/du
    template<class PartialDerivativeMatrix, class Problem, class Element, class FVGeometry, class VolVars, class SubControlVolume>
    void addSourceDerivatives(PartialDerivativeMatrix& A_ii,
                              const Problem& problem,
                              const Element& element,
                              const FVGeometry& fvGeometry,
                              const VolVars& curVolVars,
                              const SubControlVolume& scv) const
    {
        if (std::abs(spatialWeight_) < 1e-6)
            return;
        std::decay_t<PartialDerivativeMatrix> raw(0.0);
        op_.addSourceDerivatives(raw, problem, element, fvGeometry, curVolVars, scv);
        A_ii.axpy(spatialWeight_, raw);
    }

    //! flux block (whole element stencil): += spatialWeight * d(flux)/du
    template<class JacobianMatrix, class Problem, class Element, class FVGeometry, class ElemVolVars, class ElemFluxVarsCache, class SubControlVolumeFace>
    void addFluxDerivatives(JacobianMatrix& A,
                            const Problem& problem,
                            const Element& element,
                            const FVGeometry& fvGeometry,
                            const ElemVolVars& curElemVolVars,
                            const ElemFluxVarsCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
    {
        if (std::abs(spatialWeight_) < 1e-6)
            return;
        Detail::ScaledMatrixProxy<JacobianMatrix> scaledA(A, spatialWeight_);
        op_.addFluxDerivatives(scaledA, problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);
    }

    //! Robin/Neumann flux block (diagonal of inside dof): += spatialWeight * d(flux)/du
    template<class PartialDerivativeMatrices, class Problem, class Element, class FVGeometry, class ElemVolVars, class ElemFluxVarsCache, class SubControlVolumeFace>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& A_i,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVGeometry& fvGeometry,
                                 const ElemVolVars& curElemVolVars,
                                 const ElemFluxVarsCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {
        if (std::abs(spatialWeight_) < 1e-6)
            return;
        Detail::ScaledRowProxy<PartialDerivativeMatrices> scaledA_i(A_i, spatialWeight_);
        op_.addRobinFluxDerivatives(scaledA_i, problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);
    }

    // \}

    // some old interface (TODO: get rid of this)
    // (stationary is also the wrong word by the way, systems with zero
    // time derivative are in a steady-state but not necessarily stationary/not-moving at all)
    // can we decide this somehow from the temporal weight?
    bool isStationary() const
    { return false; }

private:
    LocalOperator op_;
    double spatialWeight_, temporalWeight_; // TODO: get correct type
};

} // end namespace Dumux::Experimental

#endif
