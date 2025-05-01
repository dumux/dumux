// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Volume variables deflection helper.
 */

#ifndef DUMUX_ASSEMBLY_VOLVARS_DEFLECTION_HELPER_HH
#define DUMUX_ASSEMBLY_VOLVARS_DEFLECTION_HELPER_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#ifndef DOXYGEN

namespace Dumux::Detail {

template<class GridVolumeVariables, class FVElementGeometry>
class VolVarsDeflectionHelper
{
    using ElementVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename ElementVariables::VolumeVariables;

    static constexpr int maxNumLocalDofs = Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>();

public:
    VolVarsDeflectionHelper(GridVolumeVariables& gridVolumeVariables,
                            ElementVariables& elementVariables,
                            const FVElementGeometry& fvGeometry,
                            bool deflectAllVariables)
    : gridVolumeVariables_(gridVolumeVariables)
    , elementVariables_(elementVariables)
    , fvGeometry_(fvGeometry)
    , deflectAll_(deflectAllVariables)
    {
        if (deflectAll_)
            for (const auto& scv : scvs(fvGeometry))
                origVolVars_.push_back(accessor_(scv));
    }

    template<class LocalDof>
    void setCurrent(const LocalDof& localDof)
    {
        if (!deflectAll_)
        {
            origVolVars_.clear();
            for(const auto& scv : scvs(fvGeometry_, localDof))
                origVolVars_.push_back(accessor_(scv));
        }
    }

    template<class ElementSolution, class LocalDof, class Problem>
    void deflect(const ElementSolution& elemSol,
                 const LocalDof& localDof,
                 const Problem& problem)
    {
        if (deflectAll_)
            for (const auto& scv : scvs(fvGeometry_))
                accessor_(scv).update(elemSol, problem, fvGeometry_.element(), scv);
        else
        {
            for(const auto& scv : scvs(fvGeometry_, localDof))
                accessor_(scv).update(elemSol, problem, fvGeometry_.element(), scv);
        }
    }

    template<class LocalDof>
    void restore(const LocalDof& localDof)
    {
        if (!deflectAll_)
        {
            unsigned int idx = 0;
            for(const auto& scv : scvs(fvGeometry_, localDof))
            {
                accessor_(scv) = origVolVars_[idx];
                idx++;
            }
        }
        else
            for (const auto& scv : scvs(fvGeometry_))
                accessor_(scv) = origVolVars_[scv.indexInElement()];
    }

private:
    template<class SCV>
    VolumeVariables& accessor_(const SCV& scv)
    {
        if constexpr (GridVolumeVariables::cachingEnabled)
            return gridVolumeVariables_.volVars(scv);
        else
            return elementVariables_[scv];
    }

    GridVolumeVariables& gridVolumeVariables_;
    ElementVariables& elementVariables_;
    const FVElementGeometry& fvGeometry_;
    const bool deflectAll_;
    Dune::ReservedVector<VolumeVariables, maxNumLocalDofs> origVolVars_;
};

} // end namespace Dumux::Detail

#endif // DOXYGEN

#endif
