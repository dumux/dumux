// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEDiscretization
 * \brief The element variables class
 */
#ifndef DUMUX_DISCRETIZATION_FE_ELEMENT_VARIABLES_HH
#define DUMUX_DISCRETIZATION_FE_ELEMENT_VARIABLES_HH

#include <array>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/concepts/localdofs_.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/discretization/cvfe/elementinterpolationpointdatacache.hh>
#include <dumux/discretization/cvfe/lazyelementcache.hh>
#include <dumux/discretization/cvfe/variablesdeflectionpolicy.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup FEDiscretization
 * \brief The (stencil) element variables class for finite element schemes
 * \note The class is specialized for versions with and without caching
 * \tparam GVC the grid variables cache type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVC, bool cachingEnabled>
class FEElementVariables;

/*!
 * \ingroup FEDiscretization
 * \brief The (stencil) element variables class for finite element schemes with caching
 * \note the variables are stored for the whole grid view in the corresponding GridVariables class
 */
template<class GVC>
class FEElementVariables<GVC, /*cachingEnabled*/true>
{
    class MutableVariablesView
    {
    public:
        MutableVariablesView(GVC& gridCache)
        : gridCache_(gridCache) {}

        using Variables = typename GVC::Variables;

        template<Dumux::Concept::LocalDof LocalDof>
        Variables& operator [](const LocalDof& localDof) const
        { return gridCache_.variables(localDof); }
    private:
        GVC& gridCache_;
    };

    class MutableVariablesViewWithIpCacheAccess
    {
    public:
        MutableVariablesViewWithIpCacheAccess(GVC& gridCache)
        : gridCache_(gridCache)
        {}

        using Variables = typename GVC::Variables;

        template<Dumux::Concept::LocalDof LocalDof>
        Variables& operator [](const LocalDof& localDof) const
        { return gridCache_.variables(localDof); }

        auto& cache(std::size_t eIdx) const
        { return gridCache_.cache(eIdx); }

    private:
        GVC& gridCache_;
    };

public:
    //! export type of the grid variables
    using GridVariablesCache = GVC;

    //! export type of the mutable version of the view
    using MutableView = std::conditional_t<
        GridVariablesCache::InterpolationPointData::isSolDependent,
        MutableVariablesViewWithIpCacheAccess,
        MutableVariablesView
    >;

    //! export type of the variables
    using Variables = typename GridVariablesCache::Variables;

    //! export interpolation point data
    using InterpolationPointData = typename GridVariablesCache::InterpolationPointData;

    //! export type of deflection policy
    template<class ElementDiscretization>
    using DeflectionPolicy = std::conditional_t<
        InterpolationPointData::isSolDependent,
        Dumux::Detail::CVFE::VariablesDeflectionPolicyWithIpCacheUpdate<MutableView, ElementDiscretization>,
        Dumux::Detail::CVFE::VariablesDeflectionPolicy<MutableView, ElementDiscretization>
    >;

    //! Constructor
    FEElementVariables(const GridVariablesCache& gridVariablesCache)
    : gridVariablesCachePtr_(&gridVariablesCache) {}

    const Variables& operator [](std::size_t localDofIdx) const
    { return gridVariablesCache().variables(eIdx_, localDofIdx); }

    template<Dumux::Concept::LocalDof LocalDof>
    const Variables& operator [](const LocalDof& localDof) const
    { return gridVariablesCache().variables(localDof); }

    template<class IpData>
    requires requires (const IpData& ipData) { ipData.localDofIndex(); }
    const Variables& operator [](const IpData& ipData) const
    { return gridVariablesCache().variables(eIdx_, ipData.localDofIndex()); }

    template<Concept::BoundaryFaceQpIpData IpData>
    friend const InterpolationPointData& cache(const FEElementVariables& elemVars,
                                               const IpData& ipData)
    { return elemVars.ipDataCache().boundaryFaceCache(ipData.boundaryFaceIndex(), ipData.qpIndex()); }

    template<Concept::QIpData IpData>
    friend const InterpolationPointData& cache(const FEElementVariables& elemVars,
                                               const IpData& ipData)
    { return elemVars.ipDataCache().elementCache(ipData.qpIndex()); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementDiscretization, class SolutionVector>
    FEElementVariables bind(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
                            const ElementDiscretization& elemDisc,
                            const SolutionVector& sol)  &&
    {
        this->bindElement(element, elemDisc, sol);
        return std::move(*this);
    }

    // For compatibility reasons with the case of not storing the variables.
    // function to be called before assembling an element, preparing the variables within the stencil
    template<class ElementDiscretization, class SolutionVector>
    void bind(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
              const ElementDiscretization& elemDisc,
              const SolutionVector& sol) &
    {
        bindElement(element, elemDisc, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementDiscretization, class SolutionVector>
    FEElementVariables bindElement(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
                                   const ElementDiscretization& elemDisc,
                                   const SolutionVector& sol)  &&
    {
        this->bindElement(element, elemDisc, sol);
        return std::move(*this);
    }

    // function to prepare the variables within the element
    template<class ElementDiscretization, class SolutionVector>
    void bindElement(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
                     const ElementDiscretization& elemDisc,
                     const SolutionVector& sol) &
    {
        const auto& gridDiscretization = elemDisc.gridDiscretization();
        eIdx_ = gridDiscretization.elementMapper().index(element);
        if constexpr (!GridVariablesCache::interpolationPointDataCachingEnabled)
            ipDataCache_.bind(elemDisc);
    }

    //! The grid variables cache object we are a restriction of
    const GridVariablesCache& gridVariablesCache() const
    { return *gridVariablesCachePtr_; }

    /*!
    * \brief return a local view on variables that is always mutable, regardless of the caching policy
    * \pre bind has to be called before, otherwise using the view will result in undefined behavior
    */
   MutableView asMutableView(GridVariablesCache& gridVars)
   { return { gridVars }; }

private:
    using ElementInterpolationPointDataCache = Dumux::Detail::CVFE::ElementAndBoundaryInterpolationPointDataCache<InterpolationPointData>;

    const ElementInterpolationPointDataCache& ipDataCache() const
    {
        if constexpr (GridVariablesCache::interpolationPointDataCachingEnabled)
            return gridVariablesCache().cache(eIdx_);
        else
            return ipDataCache_.get(gridVariablesCache().problem(), gridVariablesCache().localDofVariables(eIdx_));
    }

    const GridVariablesCache* gridVariablesCachePtr_;
    std::size_t eIdx_;
    Dumux::Detail::CVFE::LazyElementCache<ElementInterpolationPointDataCache, typename GridVariablesCache::Problem, std::vector<Variables>> ipDataCache_;
};

/*!
 * \ingroup FEDiscretization
 * \brief The local (stencil) element variables class for finite element schemes without caching
 */
template<class GVC>
class FEElementVariables<GVC, /*cachingEnabled*/false>
{
    using ThisType = FEElementVariables<GVC, /*cachingEnabled*/false>;
    class MutableVariablesView
    {
    public:
        MutableVariablesView(ThisType& view)
        : view_(view) {}

        using Variables = typename GVC::Variables;

        template<Dumux::Concept::LocalDof LocalDof>
        Variables& operator [](const LocalDof& localDof) const
        { return view_[localDof]; }
    private:
        ThisType& view_;
    };

    class MutableVariablesViewWithIpCacheAccess
    {
    public:
        MutableVariablesViewWithIpCacheAccess(ThisType& view)
        : view_(view)
        {}

        using Variables = typename GVC::Variables;

        template<Dumux::Concept::LocalDof LocalDof>
        Variables& operator [](const LocalDof& localDof) const
        { return view_[localDof]; }

        auto& cache(std::size_t) const
        { return view_.ipDataCache_.get(view_.gridVariablesCache().problem(), view_.variables_); }

    private:
        ThisType& view_;
    };

public:
    //! export type of the grid variables
    using GridVariablesCache = GVC;

    //! export type of the mutable version of the view
    using MutableView = std::conditional_t<
        GridVariablesCache::InterpolationPointData::isSolDependent,
        MutableVariablesViewWithIpCacheAccess,
        MutableVariablesView
    >;

    //! export type of the variables
    using Variables = typename GridVariablesCache::Variables;

    //! export interpolation point data
    using InterpolationPointData = typename GridVariablesCache::InterpolationPointData;

    //! export type of deflection policy
    template<class ElementDiscretization>
    using DeflectionPolicy = std::conditional_t<
        InterpolationPointData::isSolDependent,
        Dumux::Detail::CVFE::VariablesDeflectionPolicyWithIpCacheUpdate<MutableView, ElementDiscretization>,
        Dumux::Detail::CVFE::VariablesDeflectionPolicy<MutableView, ElementDiscretization>
    >;

    //! Constructor
    FEElementVariables(const GridVariablesCache& gridVarsCache)
    : gridVariablesCachePtr_(&gridVarsCache)
    {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementDiscretization, class SolutionVector>
    FEElementVariables bind(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
                            const ElementDiscretization& elemDisc,
                            const SolutionVector& sol)  &&
    {
        this->bindElement(element, elemDisc, sol);
        return std::move(*this);
    }

    // specialization for finite element schemes, simply forwards to the bindElement method
    template<class ElementDiscretization, class SolutionVector>
    void bind(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
              const ElementDiscretization& elemDisc,
              const SolutionVector& sol) &
    {
        bindElement(element, elemDisc, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementDiscretization, class SolutionVector>
    FEElementVariables bindElement(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
                                   const ElementDiscretization& elemDisc,
                                   const SolutionVector& sol)  &&
    {
        this->bindElement(element, elemDisc, sol);
        return std::move(*this);
    }

    // specialization for finite element schemes
    template<class ElementDiscretization, class SolutionVector>
    void bindElement(const typename ElementDiscretization::GridDiscretization::GridView::template Codim<0>::Entity& element,
                     const ElementDiscretization& elemDisc,
                     const SolutionVector& sol) &
    {
        // get the solution at the dofs of the element
        const auto& gridDiscretization = elemDisc.gridDiscretization();
        auto elemSol = elementSolution(element, sol, gridDiscretization);

        // resize variables to the required size
        variables_.resize(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));

        // update variables related to localDofs
        for (const auto& localDof : localDofs(elemDisc))
            variables_[localDof.index()].update(elemSol, gridVariablesCache().problem(), elemDisc, ipData(elemDisc, localDof));

        if constexpr (GridVariablesCache::interpolationPointDataCachingEnabled)
            eIdx_ = gridDiscretization.elementMapper().index(element);
        else
        {
            ipDataCache_.bind(elemDisc);
            if constexpr (InterpolationPointData::isSolDependent)
                ipDataCache_.get(gridVariablesCache().problem(), variables_);
        }
    }
    const Variables& operator [](std::size_t localIdx) const
    { return variables_[localIdx]; }

    Variables& operator [](std::size_t localIdx)
    { return variables_[localIdx]; }

    template<Dumux::Concept::LocalDof LocalDof>
    const Variables& operator [](const LocalDof& localDof) const
    { return variables_[localDof.index()]; }

    template<Dumux::Concept::LocalDof LocalDof>
    Variables& operator [](const LocalDof& localDof)
    { return variables_[localDof.index()]; }

    template<class IpData>
    requires requires (const IpData& ipData) { ipData.localDofIndex(); }
    const Variables& operator [](const IpData& ipData) const
    { return variables_[ipData.localDofIndex()]; }

    template<class IpData>
    requires requires (const IpData& ipData) { ipData.localDofIndex(); }
    Variables& operator [](const IpData& ipData)
    { return variables_[ipData.localDofIndex()]; }

    template<Concept::BoundaryFaceQpIpData IpData>
    friend const InterpolationPointData& cache(const FEElementVariables& elemVars,
                                               const IpData& ipData)
    { return elemVars.ipDataCache().boundaryFaceCache(ipData.boundaryFaceIndex(), ipData.qpIndex()); }

    template<Concept::QIpData IpData>
    friend const InterpolationPointData& cache(const FEElementVariables& elemVars,
                                               const IpData& ipData)
    { return elemVars.ipDataCache().elementCache(ipData.qpIndex()); }

    //! The grid variables cache object we are a restriction of
    const GridVariablesCache& gridVariablesCache() const
    { return *gridVariablesCachePtr_; }

    /*!
    * \brief return a local view on variables that is always mutable, regardless of the caching policy
    * \pre bind has to be called before, otherwise using the view will result in undefined behavior
    */
   MutableView asMutableView(GridVariablesCache&)
   { return { *this }; }

private:
    using ElementInterpolationPointDataCache = Dumux::Detail::CVFE::ElementAndBoundaryInterpolationPointDataCache<InterpolationPointData>;

    const ElementInterpolationPointDataCache& ipDataCache() const
    {
        if constexpr (GridVariablesCache::interpolationPointDataCachingEnabled)
            return gridVariablesCache().cache(eIdx_);
        else
            return ipDataCache_.get(gridVariablesCache().problem(), variables_);
    }

    const GridVariablesCache* gridVariablesCachePtr_;
    std::vector<Variables> variables_;
    std::size_t eIdx_ = 0;
    Dumux::Detail::CVFE::LazyElementCache<ElementInterpolationPointDataCache, typename GridVariablesCache::Problem, std::vector<Variables>> ipDataCache_;
};

} // end namespace Dumux::Experimental

#endif
