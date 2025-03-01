// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Deflection helper to deflect and restore grid variables.
 */

#ifndef DUMUX_VARIABLES_DEFLECTION_HELPER_HH
#define DUMUX_VARIABLES_DEFLECTION_HELPER_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Deflection helper to deflect and restore grid variables.
 */
template<class GridVariablesCache, class FVElementGeometry>
class VariablesDeflectionHelper
{
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Variables = typename ElementVariables::Variables;

    static constexpr int maxNumLocalDofs = Dumux::Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>();

public:
    VariablesDeflectionHelper(GridVariablesCache& gridVariablesCache,
                              ElementVariables& elementVariables,
                              const FVElementGeometry& fvGeometry,
                              bool deflectAllVariables)
    : gridVariablesCache_(gridVariablesCache)
    , elementVariables_(elementVariables)
    , fvGeometry_(fvGeometry)
    , deflectAll_(deflectAllVariables)
    {
        if (deflectAll_)
            for (const auto& localDof : localDofs(fvGeometry))
                origVariables_.push_back(accessor_(localDof));
    }

    template<class LocalDof>
    void setCurrent(const LocalDof& localDof)
    {
        if (!deflectAll_)
        {
            origVariables_.clear();
            origVariables_.push_back(accessor_(localDof));
        }
    }

    template<class ElementSolution, class LocalDof, class Problem>
    void deflect(const ElementSolution& elemSol,
                 const LocalDof& localDof,
                 const Problem& problem)
    {
        if (deflectAll_)
            for (const auto& localDofI : localDofs(fvGeometry_))
                accessor_(localDofI).update(elemSol, problem, fvGeometry_, localDofI);
        else
        {
            accessor_(localDof).update(elemSol, problem, fvGeometry_, localDof);
        }
    }

    template<class LocalDof>
    void restore(const LocalDof& localDof)
    {
        if (!deflectAll_)
            accessor_(localDof) = origVariables_[0];
        else
        {
            int idx = 0;
            for (const auto& localDofI : localDofs(fvGeometry_))
            {
                accessor_(localDofI) = origVariables_[idx];
                idx++;
            }
        }
    }

private:
    template<class LocalDof>
    Variables& accessor_(const LocalDof& localDof)
    {
        if constexpr (GridVariablesCache::cachingEnabled)
            return gridVariablesCache_.volVars(localDof);
        else
            return elementVariables_[localDof];
    }

    GridVariablesCache& gridVariablesCache_;
    ElementVariables& elementVariables_;
    const FVElementGeometry& fvGeometry_;
    const bool deflectAll_;
    Dune::ReservedVector<Variables, maxNumLocalDofs> origVariables_;
};

} // end namespace Dumux

#endif
