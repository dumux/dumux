// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldLevelSetModel
 * \brief Interpolation of the phase field (and its gradient) at an
 *        interpolation point of a hybrid CVFE element.
 */
#ifndef DUMUX_PHASEFIELD_LEVEL_SET_FLUX_HH
#define DUMUX_PHASEFIELD_LEVEL_SET_FLUX_HH

#include <algorithm>
#include <type_traits>
#include <dune/common/fvector.hh>

#include <dumux/phasefield/common/interpolationcontext.hh>
#include <dumux/phasefield/levelset/indices.hh>

namespace Dumux {

/*!
 * \ingroup PhaseFieldLevelSetModel
 * \brief Interpolates the phase field and its gradient at an interpolation
 *        point. A thin, named-accessor wrapper around the shared
 *        `Dumux::PhaseField::InterpolationContext`.
 * \note Shared by the control-volume (`PhaseFieldLevelSetLocalResidual::fluxIntegral`)
 *       and finite-element (`PhaseFieldLevelSetFELocalResidualTerms::addFluxAndSourceTerms`)
 *       residual paths, so the interpolation logic only has to be
 *       implemented (and kept correct) once.
 */
template<class FVElementGeometry, class ElementVariables, class IpCache>
class PhaseFieldLevelSetFluxFunctionContext
{
    using GlobalPosition = std::decay_t<decltype(std::declval<IpCache>().gradN(0u))>;
    using Scalar = typename GlobalPosition::value_type;
    using PrimaryVariables = Dune::FieldVector<Scalar, 1>;
    using Context = PhaseField::InterpolationContext<PrimaryVariables, FVElementGeometry, ElementVariables, IpCache>;

public:
    PhaseFieldLevelSetFluxFunctionContext(const FVElementGeometry& fvGeometry,
                                            const ElementVariables& elemVars,
                                            const IpCache& ipCache)
    : context_(fvGeometry, elemVars, ipCache)
    {}

    Scalar phaseField() const
    { return context_.value(PhaseFieldLevelSetIndices::phaseFieldIdx); }

    const GlobalPosition& gradPhaseField() const
    { return context_.gradient(PhaseFieldLevelSetIndices::phaseFieldIdx); }

    /*!
     * \brief The (regularized) interface unit normal n = grad(phi)/|grad(phi)|.
     *
     * Away from the interface, grad(phi) vanishes, so the gradient magnitude
     * is bounded away from zero by a small numerical regularization constant
     * to avoid a 0/0 division (the resulting normal direction is meaningless
     * there, but the compression term it multiplies, phi(1-phi), also
     * vanishes away from the interface, so this has no effect on the result).
     */
    GlobalPosition normal() const
    {
        static constexpr Scalar regularization = 1e-6;
        using std::max;
        auto n = gradPhaseField();
        n /= max(n.two_norm(), regularization);
        return n;
    }

private:
    Context context_;
};

} // end namespace Dumux

#endif
