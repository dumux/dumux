// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Variables deflection policy
 */

#ifndef DUMUX_ASSEMBLY_CVFE_VOLVARS_DEFLECTION_POLICY_HH
#define DUMUX_ASSEMBLY_CVFE_VOLVARS_DEFLECTION_POLICY_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#ifndef DOXYGEN

namespace Dumux::Detail::CVFE {

template<class MutableVariablesView, class FVElementGeometry>
class VolVarsDeflectionPolicy
{
    using VolumeVariables = typename MutableVariablesView::VolumeVariables;
    static constexpr int maxNumLocalDofs = Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>();

public:
    VolVarsDeflectionPolicy(MutableVariablesView elementVariables,
                            const FVElementGeometry& fvGeometry,
                            bool deflectAllVariables)
    : elementVariables_(elementVariables)
    , fvGeometry_(fvGeometry)
    , deflectAll_(deflectAllVariables)
    {
        if (deflectAll_)
            for (const auto& scv : scvs(fvGeometry))
                origVolVars_.push_back(elementVariables_[scv]);
    }

    template<class LocalDof>
    void setCurrent(const LocalDof& localDof)
    {
        if (!deflectAll_)
        {
            origVolVars_.clear();
            for(const auto& scv : scvs(fvGeometry_, localDof))
                origVolVars_.push_back(elementVariables_[scv]);
        }
    }

    template<class ElementSolution, class LocalDof, class Problem>
    void deflect(const ElementSolution& elemSol,
                 const LocalDof& localDof,
                 const Problem& problem)
    {
        if (deflectAll_)
            for (const auto& scv : scvs(fvGeometry_))
                elementVariables_[scv].update(elemSol, problem, fvGeometry_.element(), scv);
        else
        {
            for(const auto& scv : scvs(fvGeometry_, localDof))
                elementVariables_[scv].update(elemSol, problem, fvGeometry_.element(), scv);
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
                elementVariables_[scv] = origVolVars_[idx];
                idx++;
            }
        }
        else
            for (const auto& scv : scvs(fvGeometry_))
                elementVariables_[scv] = origVolVars_[scv.indexInElement()];
    }

private:
    MutableVariablesView elementVariables_;
    const FVElementGeometry& fvGeometry_;
    const bool deflectAll_;
    Dune::ReservedVector<VolumeVariables, maxNumLocalDofs> origVolVars_;
};

template<typename ElemVars, typename FVG>
using DefinesDeflectionPolicyType = typename ElemVars::template DeflectionPolicy<FVG>;

template<typename ElemVars, typename FVG>
constexpr inline bool definesVariablesDeflectionPolicy()
{ return Dune::Std::is_detected<DefinesDeflectionPolicyType, ElemVars, FVG>::value; }

template<class GridVarsCache, class ElemVars, class FVG>
auto makeVariablesDeflectionPolicy(GridVarsCache& gridVarsCache, ElemVars& elemVars, const FVG& fvg, bool deflectAllVolVars)
{
    if constexpr (definesVariablesDeflectionPolicy<ElemVars, FVG>())
    {
        using DeflectionPolicy = typename ElemVars::template DeflectionPolicy<FVG>;
        return DeflectionPolicy(elemVars.asMutableView(gridVarsCache), fvg, deflectAllVolVars);
    }
    else
    {
        using DeflectionPolicy = Detail::CVFE::VolVarsDeflectionPolicy<typename GridVarsCache::MutableLocalView, FVG>;
        return DeflectionPolicy(elemVars.asMutableView(gridVarsCache), fvg, deflectAllVolVars);
    }
};

} // end namespace Dumux::Detail::CVFE

#endif // DOXYGEN

#endif
