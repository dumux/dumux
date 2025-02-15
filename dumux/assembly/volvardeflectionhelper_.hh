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

#include <dumux/common/typetraits/localdofs.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#ifndef DOXYGEN

namespace Dumux::Detail {

template<class VolVarAccessor, class FVElementGeometry>
class VolVarsDeflectionHelper
{
    static constexpr int maxNumLocalDofs = Detail::maxNumLocalDofs<FVElementGeometry>();

    using SubControlVolume = typename FVElementGeometry::GridGeometry::SubControlVolume;
    using VolumeVariablesRef = std::invoke_result_t<VolVarAccessor, const SubControlVolume&>;
    using VolumeVariables = std::decay_t<VolumeVariablesRef>;
    static_assert(std::is_lvalue_reference_v<VolumeVariablesRef>
                    && !std::is_const_v<std::remove_reference_t<VolumeVariablesRef>>);

public:
    VolVarsDeflectionHelper(VolVarAccessor&& accessor,
                            const FVElementGeometry& fvGeometry,
                            bool deflectAllVolVars)
    : deflectAll_(deflectAllVolVars)
    , accessor_(std::move(accessor))
    , fvGeometry_(fvGeometry)
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

    template<class ElementSolution,
                class LocalDof,
                class Problem>
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
    const bool deflectAll_;
    VolVarAccessor accessor_;
    const FVElementGeometry& fvGeometry_;
    Dune::ReservedVector<VolumeVariables, maxNumLocalDofs> origVolVars_;
};

} // end namespace Dumux::Detail

#endif // DOXYGEN

#endif
