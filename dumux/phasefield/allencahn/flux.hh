// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AllenCahnModel
 * \brief Interpolation of the phase field (and its gradient) at an
 *        interpolation point of a hybrid CVFE element.
 */
#ifndef DUMUX_PHASEFIELD_ALLENCAHN_FLUX_HH
#define DUMUX_PHASEFIELD_ALLENCAHN_FLUX_HH

#include <type_traits>
#include <dune/common/fvector.hh>

#include <dumux/phasefield/common/interpolationcontext.hh>
#include <dumux/phasefield/allencahn/indices.hh>

namespace Dumux {

/*!
 * \ingroup AllenCahnModel
 * \brief Interpolates the phase field and its gradient at an interpolation
 *        point. A thin, named-accessor wrapper around the shared
 *        `Dumux::PhaseField::InterpolationContext`.
 * \note Shared by the control-volume (`AllenCahnLocalResidual::fluxIntegral`)
 *       and finite-element (`AllenCahnFELocalResidualTerms::addFluxAndSourceTerms`)
 *       residual paths, and by the test problem's source term, so the
 *       interpolation logic only has to be implemented (and kept correct) once.
 */
template<class FVElementGeometry, class ElementVariables, class IpCache>
class AllenCahnFluxFunctionContext
{
    using GlobalPosition = std::decay_t<decltype(std::declval<IpCache>().gradN(0u))>;
    using Scalar = typename GlobalPosition::value_type;
    using PrimaryVariables = Dune::FieldVector<Scalar, 1>;
    using Context = PhaseField::InterpolationContext<PrimaryVariables, FVElementGeometry, ElementVariables, IpCache>;

public:
    AllenCahnFluxFunctionContext(const FVElementGeometry& fvGeometry,
                                 const ElementVariables& elemVars,
                                 const IpCache& ipCache)
    : context_(fvGeometry, elemVars, ipCache)
    {}

    Scalar phaseField() const
    { return context_.value(AllenCahnIndices::phaseFieldIdx); }

    const GlobalPosition& gradPhaseField() const
    { return context_.gradient(AllenCahnIndices::phaseFieldIdx); }

private:
    Context context_;
};

} // end namespace Dumux

#endif
