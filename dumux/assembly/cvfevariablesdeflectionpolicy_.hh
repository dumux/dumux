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

#ifndef DUMUX_ASSEMBLY_CVFE_VARIABLES_DEFLECTION_POLICY_HH
#define DUMUX_ASSEMBLY_CVFE_VARIABLES_DEFLECTION_POLICY_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/volvarscontext.hh>

#ifndef DOXYGEN

namespace Dumux::Detail::CVFE {

template<class VariablesContext, class FVElementGeometry>
class VariablesDeflectionPolicy
{
    using Variables = typename VariablesContext::Variables;

    static constexpr int maxNumLocalDofs = Dumux::Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>();

public:
    VariablesDeflectionPolicy(VariablesContext& context,
                              const FVElementGeometry& fvGeometry)
    : context_(context)
    , fvGeometry_(fvGeometry)
    {
        if(context_.dependsOnAllElementDofs())
            for(const auto& variable : context_.variables())
                origVariables_.push_back(variable);
    }

    template<class LocalDof>
    void store(const LocalDof& localDof)
    {
        if (!context_.dependsOnAllElementDofs())
        {
            origVariables_.clear();
            for(const auto& variable : context_.variables(localDof))
                origVariables_.push_back(variable);
        }
    }

    template<class ElementSolution, class LocalDof, class Problem>
    void update(const ElementSolution& elemSol,
                const LocalDof& localDof,
                const Problem& problem)
    {
        context_.update(elemSol, problem, localDof);
    }

    template<class LocalDof>
    void restore(const LocalDof& localDof)
    {
        unsigned int idx = 0;
        for(auto& variable : context_.variables(localDof))
        {
            variable = origVariables_[idx];
            idx++;
        }
    }

private:
    VariablesContext& context_;
    const FVElementGeometry& fvGeometry_;
    Dune::ReservedVector<Variables, maxNumLocalDofs> origVariables_;
};

template<typename ElemVars, typename FVG>
using DefinesVariablesContextType = typename ElemVars::template VariablesContext<FVG>;

template<typename ElemVars, typename FVG>
constexpr inline bool definesVariablesContext()
{ return Dune::Std::is_detected<DefinesVariablesContextType, ElemVars, FVG>::value; }

template<class GridVarsCache, class ElemVars, class FVG>
auto makeVariablesContext(GridVarsCache& gridVarsCache, ElemVars& elemVars, const FVG& fvg, bool localDependency)
{
    if constexpr (definesVariablesContext<ElemVars, FVG>())
    {
        using VariablesContext = typename ElemVars::template VariablesContext<FVG>;
        return VariablesContext(elemVars.asMutableView(gridVarsCache), fvg, localDependency);
    }
    else
    {
        using VariablesContext = Detail::CVFE::LocalDofVolVarsContext<typename GridVarsCache::MutableLocalView, FVG>;
        return VariablesContext(elemVars.asMutableView(gridVarsCache), fvg, localDependency);
    }
};

} // end namespace Dumux::Detail::CVFE

#endif // DOXYGEN

#endif
