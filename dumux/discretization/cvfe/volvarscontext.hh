// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Context for grid volume variables.
 */

#ifndef DUMUX_DISCRETIZATION_CVFE_LOCALDOF_VOLVARS_CONTEXT_HH
#define DUMUX_DISCRETIZATION_CVFE_LOCALDOF_VOLVARS_CONTEXT_HH

#include <ranges>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Detail::CVFE {

/*!
 * \ingroup Discretization
 * \brief Context for grid volume variables.
 */
template<class MutableVariablesView, class FVElementGeometry>
class LocalDofVolVarsContext
{
    using LocalDof = Dumux::Detail::LocalDofs::LocalDof_t<FVElementGeometry>;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using LocalIndexType = typename SubControlVolume::Traits::LocalIndexType;
public:
    using Variables = typename MutableVariablesView::VolumeVariables;

    LocalDofVolVarsContext(MutableVariablesView elementVariables,
                           const FVElementGeometry& fvGeometry,
                           bool dependsOnAllElementDofs = false)
    : elementVariables_(elementVariables)
    , fvGeometry_(fvGeometry)
    , dependsOnAllElementDofs_(dependsOnAllElementDofs)
    { }

    //! Return all variables (reference) related to the element
    std::ranges::range auto variables()
    {
        return makeVariablesRange_(0u, fvGeometry_.numScv(), [&](LocalIndexType i) -> const SubControlVolume& {return fvGeometry_.scv(i);});
    }

    //! Return all variables (reference) related to a localDof
    std::ranges::range auto variables(const LocalDof& localDof)
    {
        if (!dependsOnAllElementDofs())
        {
            auto&& scvRange = scvs(fvGeometry_, localDof);
            auto numScvsForLocalDof = scvRange.size();
            return makeVariablesRange_(0u, numScvsForLocalDof, [scvs = std::move(scvRange)](LocalIndexType i) -> const SubControlVolume& {return scvs[i];});
        }
        else
            return variables();
    }

    //! Update all variables that depend on localDof
    template<class ElementSolution, class Problem>
    void updateVariables(const ElementSolution& elemSol,
                         const Problem& problem,
                         const LocalDof& localDof)
    {
        if (!dependsOnAllElementDofs())
            for(const auto& scv : scvs(fvGeometry_, localDof))
                elementVariables_[scv].update(elemSol, problem, fvGeometry_.element(), scv);
        else
            for(const auto& scv : scvs(fvGeometry_))
                elementVariables_[scv].update(elemSol, problem, fvGeometry_.element(), scv);
    }

    bool dependsOnAllElementDofs() const
    { return dependsOnAllElementDofs_;}

private:
    std::ranges::range auto makeVariablesRange_(LocalIndexType first, LocalIndexType last, std::function<const SubControlVolume&(LocalIndexType)>&& getScv)
    {
        return std::views::all(std::views::iota(first, last)) |
               std::views::transform([this, scv = std::move(getScv)](const LocalIndexType& i) -> Variables& {
                   return this->elementVariables_[scv(i)];
               });
    }

    MutableVariablesView elementVariables_;
    const FVElementGeometry& fvGeometry_;
    bool dependsOnAllElementDofs_;
};

} // end namespace Dumux::Detail::CVFE

#endif
