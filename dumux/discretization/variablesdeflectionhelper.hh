// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

#include <dumux/common/typetraits/localdofs.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#ifndef DOXYGEN

namespace Dumux {

template<class GridLocalVariables, class FVElementGeometry>
class VariablesDeflectionHelper
{
    using ElementLocalVariables = typename GridLocalVariables::LocalView;
    using LocalVariables = typename ElementLocalVariables::LocalVariables;

    static constexpr int maxNumLocalDofs = Detail::maxNumLocalDofs<FVElementGeometry>();

public:
    VariablesDeflectionHelper(GridLocalVariables& gridLocalVariables,
                              ElementLocalVariables& elementLocalVariables,
                              const FVElementGeometry& fvGeometry,
                              bool deflectAllVariables)
    : gridLocalVariables_(gridLocalVariables)
    , elementLocalVariables_(elementLocalVariables)
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
    LocalVariables& accessor_(const LocalDof& localDof)
    {
        if constexpr (GridLocalVariables::cachingEnabled)
            return gridLocalVariables_.localVars(localDof);
        else
            return elementLocalVariables_[localDof];
    }

    GridLocalVariables& gridLocalVariables_;
    ElementLocalVariables& elementLocalVariables_;
    const FVElementGeometry& fvGeometry_;
    const bool deflectAll_;
    Dune::ReservedVector<LocalVariables, maxNumLocalDofs> origVariables_;
};

} // end namespace Dumux

#endif // DOXYGEN

#endif
