// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Context for grid variables.
 */

#ifndef DUMUX_DISCRETIZATION_CVFE_LOCALDOF_VARIABLES_CONTEXT_HH
#define DUMUX_DISCRETIZATION_CVFE_LOCALDOF_VARIABLES_CONTEXT_HH

#include <ranges>

#include <dumux/common/typetraits/localdofs_.hh>

namespace Dumux::Detail::CVFE {

/*!
 * \ingroup Discretization
 * \brief Context for grid variables.
 */
template<class MutableVariablesView, class FVElementGeometry>
class LocalDofVariablesContext
{
    using LocalDof = Dumux::Detail::LocalDofs::LocalDof_t<FVElementGeometry>;

    using LocalIndexType = typename LocalDof::LocalIndexType;
public:
    using Variables = typename MutableVariablesView::Variables;

    LocalDofVariablesContext(MutableVariablesView elementVariables,
                             const FVElementGeometry& fvGeometry,
                             bool dependsOnAllElementDofs = false)
    : elementVariables_(elementVariables)
    , fvGeometry_(fvGeometry)
    , dependsOnAllElementDofs_(dependsOnAllElementDofs)
    { }

    //! Return all variables (reference) related to the element
    std::ranges::range auto variables()
    {
        auto&& localDofRange = localDofs(fvGeometry_);
        auto numLocalDofs = localDofRange.size();
        return makeVariablesRange_(0u, numLocalDofs, [localDofs = std::move(localDofRange)](LocalIndexType i){return localDofs[i];});
    }

    //! Return all variables (reference) related to a localDof
    std::ranges::range auto variables(const LocalDof& localDof)
    {
        if(!dependsOnAllElementDofs())
            return makeVariablesRange_(0u, 1u, [&](LocalIndexType i){return localDof;});
        else
            return variables();
    }

    //! Update all variables that depend on localDof
    template<class ElementSolution, class Problem>
    void updateVariables(const ElementSolution& elemSol,
                         const Problem& problem,
                         const LocalDof& localDof)
    {
        if(!dependsOnAllElementDofs())
            elementVariables_[localDof].update(elemSol, problem, fvGeometry_, localDof);
        else
            for (const auto& localDofI : localDofs(fvGeometry_))
                elementVariables_[localDofI].update(elemSol, problem, fvGeometry_, localDofI);
    }

    bool dependsOnAllElementDofs() const
    { return dependsOnAllElementDofs_;}

private:
    std::ranges::range auto makeVariablesRange_(LocalIndexType first, LocalIndexType last, std::function<LocalDof(LocalIndexType)>&& getLocalDof)
    {
        return std::views::all(std::views::iota(first, last)) |
               std::views::transform([this, localDof = std::move(getLocalDof)](const LocalIndexType& i) -> Variables& {
                   return this->elementVariables_[localDof(i)];
               });
    }

    MutableVariablesView elementVariables_;
    const FVElementGeometry& fvGeometry_;
    bool dependsOnAllElementDofs_;
};

} // end namespace Dumux::Detail::CVFE

#endif
