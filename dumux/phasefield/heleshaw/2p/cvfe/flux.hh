// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawModel
 * \brief Interpolation of pressure, phase field, and chemical potential (and
 *        their gradients) at an interpolation point of a hybrid CVFE element.
 */
#ifndef DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_FLUX_HH
#define DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_FLUX_HH

#include <type_traits>
#include <dune/common/fvector.hh>

#include <dumux/phasefield/common/interpolationcontext.hh>
#include <dumux/phasefield/heleshaw/2p/cvfe/indices.hh>

namespace Dumux {

/*!
 * \ingroup HeleShawModel
 * \brief Interpolates pressure, phase field, chemical potential, and their
 *        gradients at an interpolation point. A thin, named-accessor wrapper
 *        around the shared `Dumux::PhaseField::InterpolationContext`.
 */
template<class FVElementGeometry, class ElementVariables, class IpCache>
class HeleShawTwoPCVFEFluxFunctionContext
{
    using GlobalPosition = std::decay_t<decltype(std::declval<IpCache>().gradN(0u))>;
    using Scalar = typename GlobalPosition::value_type;
    using PrimaryVariables = Dune::FieldVector<Scalar, 3>;
    using Context = PhaseField::InterpolationContext<PrimaryVariables, FVElementGeometry, ElementVariables, IpCache>;

public:
    HeleShawTwoPCVFEFluxFunctionContext(const FVElementGeometry& fvGeometry,
                                        const ElementVariables& elemVars,
                                        const IpCache& ipCache)
    : context_(fvGeometry, elemVars, ipCache)
    {}

    Scalar pressure() const
    { return context_.value(HeleShawTwoPCVFEIndices::pressureIdx); }

    Scalar phaseField() const
    { return context_.value(HeleShawTwoPCVFEIndices::phaseFieldIdx); }

    Scalar chemicalPotential() const
    { return context_.value(HeleShawTwoPCVFEIndices::chemPotIdx); }

    const GlobalPosition& gradPressure() const
    { return context_.gradient(HeleShawTwoPCVFEIndices::pressureIdx); }

    const GlobalPosition& gradPhaseField() const
    { return context_.gradient(HeleShawTwoPCVFEIndices::phaseFieldIdx); }

    const GlobalPosition& gradChemicalPotential() const
    { return context_.gradient(HeleShawTwoPCVFEIndices::chemPotIdx); }

private:
    Context context_;
};

} // end namespace Dumux

#endif
