// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Variables deflection policy
 */

#ifndef DUMUX_DISCRETIZATION_CVFE_VARIABLES_DEFLECTION_POLICY_HH
#define DUMUX_DISCRETIZATION_CVFE_VARIABLES_DEFLECTION_POLICY_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Detail::CVFE {

template<class MutableVariablesView, class FVElementGeometry>
class VariablesDeflectionPolicy
{
    using Variables = typename MutableVariablesView::Variables;
    static constexpr int maxNumLocalDofs = Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>();

public:
    VariablesDeflectionPolicy(MutableVariablesView elementVariables,
                              const FVElementGeometry& fvGeometry,
                              bool deflectAllVariables)
    : elementVariables_(elementVariables)
    , fvGeometry_(fvGeometry)
    , deflectAll_(deflectAllVariables)
    {
        if (deflectAll_)
            for (const auto& localDof : localDofs(fvGeometry))
                origVariables_.push_back(elementVariables_[localDof]);
    }

    template<class LocalDof>
    void store(const LocalDof& localDof)
    {
        if (!deflectAll_)
        {
            origVariables_.clear();
            origVariables_.push_back(elementVariables_[localDof]);
        }
    }

    template<class ElementSolution, class LocalDof, class Problem>
    void update(const ElementSolution& elemSol,
                const LocalDof& localDof,
                const Problem& problem)
    {
        if (deflectAll_)
            for (const auto& localDofI : localDofs(fvGeometry_))
                elementVariables_[localDofI].update(elemSol, problem, fvGeometry_, localDofI);
        else
        {
            elementVariables_[localDof].update(elemSol, problem, fvGeometry_, localDof);
        }
    }

    template<class LocalDof>
    void restore(const LocalDof& localDof)
    {
        if (!deflectAll_)
            elementVariables_[localDof] = origVariables_[0];
        else
        {
            int idx = 0;
            for (const auto& localDofI : localDofs(fvGeometry_))
            {
                elementVariables_[localDofI] = origVariables_[idx];
                idx++;
            }
        }
    }

private:
    MutableVariablesView elementVariables_;
    const FVElementGeometry& fvGeometry_;
    const bool deflectAll_;
    Dune::ReservedVector<Variables, maxNumLocalDofs> origVariables_;
};

} // end namespace Dumux::Detail::CVFE

#endif
