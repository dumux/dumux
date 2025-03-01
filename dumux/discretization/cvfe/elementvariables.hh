// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The element variables class
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_ELEMENT_VARIABLES_HH
#define DUMUX_DISCRETIZATION_CVFE_ELEMENT_VARIABLES_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/elementsolution.hh>

#include "variablesdeflectionpolicy.hh"

namespace Dumux::Detail::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief The (stencil) element variables class for control-volume finite element
 * \note The class is specialized for versions with and without caching
 * \tparam GVC the grid variables cache type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVC, bool cachingEnabled>
class CVFEElementVariables;

/*!
 * \ingroup CVFEDiscretization
 * \brief The (stencil) element variables class for control-volume finite element with caching
 * \note the variables are stored for the whole grid view in the corresponding GridVariables class
 */
template<class GVC>
class CVFEElementVariables<GVC, /*cachingEnabled*/true>
{
    class MutableVariablesView
    {
    public:
        MutableVariablesView(GVC& gridCache)
        : gridCache_(gridCache) {}

        using Variables = typename GVC::VolumeVariables;

        template<class LocalDof>
        Variables& operator [](const LocalDof& localDof) const
        { return gridCache_.volVars(localDof.elementIndex(), localDof.index()); }
    private:
        GVC& gridCache_;
    };

public:
    //! export type of the grid variables
    using GridVariablesCache = GVC;
    //! TODO: Delete after new variables concept has been implemented
    //! this is currently needed to support old interface
    using GridVolumeVariables = GridVariablesCache;

    //! export type of the mutable version of the view
    using MutableView = MutableVariablesView;

    //! export type of the variables
    using Variables = typename GridVariablesCache::VolumeVariables;
    //! TODO: Delete after new variables concept has been implemented
    //! this is currently needed to support old interface
    using VolumeVariables = Variables;

    //! export type of deflection policy
    template<class FVElementGeometry>
    using DeflectionPolicy = VariablesDeflectionPolicy<MutableView, FVElementGeometry>;

    //! Constructor
    CVFEElementVariables(const GridVariablesCache& gridVariablesCache)
    : gridVariablesCachePtr_(&gridVariablesCache) {}

    const Variables& operator [](std::size_t localDofIdx) const
    { return gridVolVars().volVars(eIdx_, localDofIdx); }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    const Variables& operator [](const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Dumux::Detail::LocalDofs::isLocalDofType<ScvOrLocalDof>())
            return gridVolVars().volVars(eIdx_, scvOrLocalDof.index());
        else
            return gridVolVars().volVars(eIdx_, scvOrLocalDof.indexInElement());
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                              const FVElementGeometry& fvGeometry,
                              const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // For compatibility reasons with the case of not storing the variables.
    // function to be called before assembling an element, preparing the variables within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        bindElement(element, fvGeometry, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                     const FVElementGeometry& fvGeometry,
                                     const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // function to prepare the variables within the element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    {
        eIdx_ = fvGeometry.gridGeometry().elementMapper().index(element);
    }

    //! The grid variables cache object we are a restriction of
    //! TODO: Rename after new variables concept has been implemented
    //! this is currently needed to support old interface
    const GridVariablesCache& gridVolVars() const
    { return *gridVariablesCachePtr_; }

    /*!
    * \brief return a local view on variables that is always mutable, regardless of the caching policy
    * \pre bind has to be called before, otherwise using the view will result in undefined behavior
    */
   MutableView asMutableView(GridVariablesCache& gridVars)
   { return { gridVars }; }

private:
    const GridVariablesCache* gridVariablesCachePtr_;
    std::size_t eIdx_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) element variables class for control-volume finite element without caching
 */
template<class GVC>
class CVFEElementVariables<GVC, /*cachingEnabled*/false>
{
    using ThisType = CVFEElementVariables<GVC, /*cachingEnabled*/false>;

    class MutableVariablesView
    {
    public:
        MutableVariablesView(ThisType& view)
        : view_(view) {}

        using Variables = typename GVC::VolumeVariables;

        template<class LocalDof>
        Variables& operator [](const LocalDof& localDof) const
        { return view_[localDof]; }
    private:
        ThisType& view_;
    };

public:
    //! export type of the grid variables
    using GridVariablesCache = GVC;
    //! TODO: Delete after new variables concept has been implemented
    //! this is currently needed to support old interface
    using GridVolumeVariables = GridVariablesCache;

    //! export type of the mutable version of the view
    using MutableView = MutableVariablesView;

    //! export type of the variables
    using Variables = typename GridVariablesCache::VolumeVariables;
    //! TODO: Delete after new variables concept has been implemented
    //! this is currently needed to support old interface
    using VolumeVariables = Variables;

    //! export type of deflection policy
    template<class FVElementGeometry>
    using DeflectionPolicy = VariablesDeflectionPolicy<MutableView, FVElementGeometry>;

    //! Constructor
    CVFEElementVariables(const GridVariablesCache& gridVarsCache)
    : gridVariablesCachePtr_(&gridVarsCache) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                              const FVElementGeometry& fvGeometry,
                              const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // specialization for control-volume finite element, simply forwards to the bindElement method
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        bindElement(element, fvGeometry, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                          const FVElementGeometry& fvGeometry,
                                          const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // specialization for control-volume finite element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    {
        // get the solution at the dofs of the element
        auto elemSol = elementSolution(element, sol, fvGeometry.gridGeometry());

        // resize variables to the required size
        variables_.resize(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));

        // update cv related dofs where there exists a localDof
        for (const auto& localDof : localDofs(fvGeometry))
            variables_[localDof.index()].update(elemSol, gridVolVars().problem(), fvGeometry, localDof);
    }

    const Variables& operator [](std::size_t localIdx) const
    { return variables_[localIdx]; }

    Variables& operator [](std::size_t localIdx)
    { return variables_[localIdx]; }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    const Variables& operator [](const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Dumux::Detail::LocalDofs::isLocalDofType<ScvOrLocalDof>())
            return variables_[scvOrLocalDof.index()];
        else
            return variables_[scvOrLocalDof.indexInElement()];
    }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    Variables& operator [](const ScvOrLocalDof& scvOrLocalDof)
    {
        if constexpr (Dumux::Detail::LocalDofs::isLocalDofType<ScvOrLocalDof>())
            return variables_[scvOrLocalDof.index()];
        else
            return variables_[scvOrLocalDof.indexInElement()];
    }

    //! The grid variables cache object we are a restriction of
    //! TODO: Rename after new variables concept has been implemented
    //! this is currently needed to support old interface
    const GridVariablesCache& gridVolVars() const
    { return *gridVariablesCachePtr_; }

    /*!
    * \brief return a local view on variables that is always mutable, regardless of the caching policy
    * \pre bind has to be called before, otherwise using the view will result in undefined behavior
    */
   MutableView asMutableView(GridVariablesCache&)
   { return { *this }; }

private:
    const GridVariablesCache* gridVariablesCachePtr_;
    std::vector<Variables> variables_;
};

} // end namespace Dumux

#endif
