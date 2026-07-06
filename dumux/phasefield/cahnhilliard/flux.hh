// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CahnHilliardModel
 * \brief Interpolation of concentration and chemical potential (and their
 *        gradients) at an interpolation point of a hybrid CVFE element.
 */
#ifndef DUMUX_PHASEFIELD_CAHNHILLIARD_FLUX_HH
#define DUMUX_PHASEFIELD_CAHNHILLIARD_FLUX_HH

#include <type_traits>
#include <dune/common/fvector.hh>

#include <dumux/phasefield/common/interpolationcontext.hh>
#include <dumux/phasefield/cahnhilliard/indices.hh>

namespace Dumux {

/*!
 * \ingroup CahnHilliardModel
 * \brief Interpolates concentration, chemical potential, and their gradients
 *        at an interpolation point. A thin, named-accessor wrapper around the
 *        shared `Dumux::PhaseField::InterpolationContext`.
 * \note Shared by the control-volume (`CahnHilliardLocalResidual::fluxIntegral`)
 *       and finite-element (`CahnHilliardFELocalResidualTerms::addFluxAndSourceTerms`)
 *       residual paths, and by the test problem's source term, so the
 *       interpolation logic only has to be implemented (and kept correct) once.
 */
template<class FVElementGeometry, class ElementVariables, class IpCache>
class CahnHilliardFluxFunctionContext
{
    using GlobalPosition = std::decay_t<decltype(std::declval<IpCache>().gradN(0u))>;
    using Scalar = typename GlobalPosition::value_type;
    using PrimaryVariables = Dune::FieldVector<Scalar, 2>;
    using Context = PhaseField::InterpolationContext<PrimaryVariables, FVElementGeometry, ElementVariables, IpCache>;

public:
    CahnHilliardFluxFunctionContext(const FVElementGeometry& fvGeometry,
                                    const ElementVariables& elemVars,
                                    const IpCache& ipCache)
    : context_(fvGeometry, elemVars, ipCache)
    {}

    Scalar concentration() const
    { return context_.value(CahnHilliardIndices::concentrationIdx); }

    Scalar chemicalPotential() const
    { return context_.value(CahnHilliardIndices::chemicalPotentialIdx); }

    const GlobalPosition& gradConcentration() const
    { return context_.gradient(CahnHilliardIndices::concentrationIdx); }

    const GlobalPosition& gradChemicalPotential() const
    { return context_.gradient(CahnHilliardIndices::chemicalPotentialIdx); }

private:
    Context context_;
};

} // end namespace Dumux

#endif
