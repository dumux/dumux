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
 * \brief A local operator wrapper for the implicit solves of time integration methods
 */
#ifndef DUMUX_EXPERIMENTAL_MULTISTAGE_FV_LOCAL_OPERATOR_HH
#define DUMUX_EXPERIMENTAL_MULTISTAGE_FV_LOCAL_OPERATOR_HH

#include <cmath>
#include <type_traits>

#include <dumux/discretization/extrusion.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \brief Weights the operators of a finite volume local residual
 *
 * The storage term is `computeStorage(problem, scv, volVars)` and, for problems that enable
 * inertia (see ImplicitSolveParams), the inertia term is `computeInertia(problem, element, scv, volVars)`,
 * both integrated over the sub-control volumes.
 */
template<class LocalOperator>
class MultiStageFVLocalOperator
{
    using ElementOperatorResultVector = typename LocalOperator::ElementResidualVector;
    using OperatorProblem = std::decay_t<decltype(std::declval<const LocalOperator&>().problem())>;
    static constexpr bool enableInertia = Detail::TimeStepping::enablesInertia<OperatorProblem>();

public:
    // compatibility
    using ElementResidualVector = ElementOperatorResultVector;

    MultiStageFVLocalOperator(const LocalOperator& op)
    : op_(op)
    , spatialWeight_(1.0)
    , temporalWeight_(1.0)
    , inertiaWeight_(1.0)
    {}

    //! The weighted terms with time derivatives (storage and inertia)
    template<class FVGeometry, class ElemVolVars>
    ElementOperatorResultVector evalStorage(
        const FVGeometry& fvGeometry,
        const ElemVolVars& elemVolVars
    ) const {
        ElementOperatorResultVector result(fvGeometry.numScv());

        if (!skipTemporal_)
            addStorage_(result, fvGeometry, elemVolVars, temporalWeight_);

        if constexpr (enableInertia)
            if (!skipInertia_)
                addInertia_(result, fvGeometry, elemVolVars, inertiaWeight_);

        return result;
    }

    //! The unweighted storage operator
    template<class FVGeometry, class ElemVolVars>
    ElementOperatorResultVector evalStorageOperator(
        const FVGeometry& fvGeometry,
        const ElemVolVars& elemVolVars
    ) const {
        ElementOperatorResultVector result(fvGeometry.numScv());
        addStorage_(result, fvGeometry, elemVolVars, 1.0);
        return result;
    }

    //! The unweighted inertia operator (zero if the problem does not enable inertia)
    template<class FVGeometry, class ElemVolVars>
    ElementOperatorResultVector evalInertiaOperator(
        const FVGeometry& fvGeometry,
        const ElemVolVars& elemVolVars
    ) const {
        ElementOperatorResultVector result(fvGeometry.numScv());
        if constexpr (enableInertia)
            addInertia_(result, fvGeometry, elemVolVars, 1.0);
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
        if (!skipSpatial_)
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
        if (!skipSpatial_)
        {
            result = op_.evalFlux(op_.problem(), element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
            result *= spatialWeight_;
        }
        return result;
    }

    //! set the weight of the spatial operator and whether it is structurally zero
    void spatialWeight(double w, bool skip = false) { spatialWeight_ = w; skipSpatial_ = skip; }
    double spatialWeight() const { return spatialWeight_; }
    bool skipSpatial() const { return skipSpatial_; }

    //! set the weight of the temporal operator and whether it is structurally zero
    void temporalWeight(double w, bool skip = false) { temporalWeight_ = w; skipTemporal_ = skip; }
    double temporalWeight() const { return temporalWeight_; }
    bool skipTemporal() const { return skipTemporal_; }

    //! set the weight of the inertia operator and whether it is structurally zero
    void inertiaWeight(double w, bool skip = false) { inertiaWeight_ = w; skipInertia_ = skip; }
    double inertiaWeight() const { return inertiaWeight_; }
    bool skipInertia() const { return skipInertia_; }

    const auto& problem() const
    { return op_.problem(); }

    // some old interface (TODO: get rid of this)
    // (stationary is also the wrong word by the way, systems with zero
    // time derivative are in a steady-state but not necessarily stationary/not-moving at all)
    // can we decide this somehow from the temporal weight?
    bool isStationary() const
    { return false; }

private:
    template<class FVGeometry, class ElemVolVars>
    void addStorage_(ElementOperatorResultVector& result,
                     const FVGeometry& fvGeometry,
                     const ElemVolVars& elemVolVars,
                     const double weight) const
    {
        for (const auto& scv : scvs(fvGeometry))
            result[scv.localDofIndex()] +=
                op_.computeStorage(op_.problem(), scv, elemVolVars[scv])
                * elemVolVars[scv].extrusionFactor()
                * Extrusion_t<typename FVGeometry::GridGeometry>::volume(fvGeometry, scv)
                * weight;
    }

    template<class FVGeometry, class ElemVolVars>
    void addInertia_(ElementOperatorResultVector& result,
                     const FVGeometry& fvGeometry,
                     const ElemVolVars& elemVolVars,
                     const double weight) const
    {
        static_assert(requires { op_.computeInertia(op_.problem(), fvGeometry.element(), *scvs(fvGeometry).begin(), elemVolVars[*scvs(fvGeometry).begin()]); },
            "The problem enables inertia, so the local residual has to implement computeInertia(problem, element, scv, volVars)");

        for (const auto& scv : scvs(fvGeometry))
            result[scv.localDofIndex()] +=
                op_.computeInertia(op_.problem(), fvGeometry.element(), scv, elemVolVars[scv])
                * elemVolVars[scv].extrusionFactor()
                * Extrusion_t<typename FVGeometry::GridGeometry>::volume(fvGeometry, scv)
                * weight;
    }

    LocalOperator op_;
    double spatialWeight_, temporalWeight_, inertiaWeight_; // TODO: get correct type
    bool skipSpatial_ = false, skipTemporal_ = false, skipInertia_ = false;
};

} // end namespace Dumux::Experimental

#endif
