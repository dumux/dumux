// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief Interpolation of the primary variables (and their gradients) at an
 *        interpolation point of a hybrid CVFE element, shared across
 *        phase-field models.
 */
#ifndef DUMUX_PHASEFIELD_COMMON_INTERPOLATION_CONTEXT_HH
#define DUMUX_PHASEFIELD_COMMON_INTERPOLATION_CONTEXT_HH

#include <array>
#include <type_traits>

#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::PhaseField {

/*!
 * \ingroup PhaseFieldModels
 * \brief Interpolates the primary variables and their gradients at an
 *        interpolation point from the local dof values and the shape
 *        function/gradient data cached at that point.
 * \tparam PrimaryVariables The primary variable vector type (fixes the number of equations)
 */
template<class PrimaryVariables, class FVElementGeometry, class ElementVariables, class IpCache>
class InterpolationContext
{
    using GlobalPosition = std::decay_t<decltype(std::declval<IpCache>().gradN(0u))>;
    using Scalar = typename PrimaryVariables::value_type;
    static constexpr int numEq = PrimaryVariables::dimension;

public:
    InterpolationContext(const FVElementGeometry& fvGeometry,
                        const ElementVariables& elemVars,
                        const IpCache& ipCache)
    : value_(0.0)
    {
        gradient_.fill(GlobalPosition(0.0));

        const auto& shapeValues = ipCache.shapeValues();
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& vars = elemVars[localDof];
            const auto& gradN = ipCache.gradN(localDof.index());
            const auto shapeValue = shapeValues[localDof.index()][0];

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                value_[eqIdx] += shapeValue*vars.priVars()[eqIdx];
                gradient_[eqIdx].axpy(vars.priVars()[eqIdx], gradN);
            }
        }
    }

    //! the interpolated value of primary variable `eqIdx`
    Scalar value(int eqIdx) const
    { return value_[eqIdx]; }

    //! the interpolated gradient of primary variable `eqIdx`
    const GlobalPosition& gradient(int eqIdx) const
    { return gradient_[eqIdx]; }

    //! the interpolated primary variable vector
    const PrimaryVariables& values() const
    { return value_; }

private:
    PrimaryVariables value_;
    std::array<GlobalPosition, numEq> gradient_;
};

} // end namespace Dumux::PhaseField

#endif
