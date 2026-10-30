// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief An element cache that is updated on its first access after binding
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_LAZY_ELEMENT_CACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_LAZY_ELEMENT_CACHE_HH

#include <utility>

namespace Dumux::Detail::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief An element cache that is updated on its first access after binding
 * \note Accessing the cache of a bound object is not thread-safe.
 * \tparam Cache the cache, providing `update(problem, element, elemDisc, elemVars)`
 * \tparam Problem the problem
 * \tparam ElementVariables the variables of the element the cache is computed from
 */
template<class Cache, class Problem, class ElementVariables>
class LazyElementCache
{
    using Update = void(*)(Cache&, const void*, const Problem&, const ElementVariables&);

public:
    /*!
     * \brief Bind to an element discretization, the cache is updated on the next access
     * \note The element discretization has to stay bound to the same element while the cache is accessed.
     */
    template<class ElementDiscretization>
    void bind(const ElementDiscretization& elemDisc)
    {
        elemDisc_ = &elemDisc;
        update_ = [](Cache& cache, const void* elemDiscPtr, const Problem& problem, const ElementVariables& elemVars)
        {
            const auto& elemDisc = *static_cast<const ElementDiscretization*>(elemDiscPtr);
            cache.update(problem, elemDisc.element(), elemDisc, elemVars);
        };
        isUpToDate_ = false;
    }

    //! The cache, updated if this is the first access after binding
    const Cache& get(const Problem& problem, const ElementVariables& elemVars) const
    {
        if (!isUpToDate_)
        {
            update_(cache_, elemDisc_, problem, elemVars);
            isUpToDate_ = true;
        }
        return cache_;
    }

    //! The cache, updated if this is the first access after binding
    Cache& get(const Problem& problem, const ElementVariables& elemVars)
    {
        std::as_const(*this).get(problem, elemVars);
        return cache_;
    }

private:
    mutable Cache cache_;
    mutable bool isUpToDate_ = false;
    const void* elemDisc_ = nullptr;
    Update update_ = nullptr;
};

} // end namespace Dumux::Detail::CVFE

#endif
