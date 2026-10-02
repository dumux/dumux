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

        if (!skipTemporal_)
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

    const auto& problem() const
    { return op_.problem(); }

    // some old interface (TODO: get rid of this)
    // (stationary is also the wrong word by the way, systems with zero
    // time derivative are in a steady-state but not necessarily stationary/not-moving at all)
    // can we decide this somehow from the temporal weight?
    bool isStationary() const
    { return false; }

private:
    LocalOperator op_;
    double spatialWeight_, temporalWeight_; // TODO: get correct type
    bool skipSpatial_ = false, skipTemporal_ = false;
};

} // end namespace Dumux::Experimental

#endif
