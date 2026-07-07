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
#include <memory>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/concepts/localdofs_.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

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
    { return elemVars.gridVariablesCache().boundaryFaceCache(elemVars.eIdx_, ipData.boundaryFaceIndex(), ipData.qpIndex()); }

    template<Concept::QIpData IpData>
    friend const InterpolationPointData& cache(const FEElementVariables& elemVars,
                                               const IpData& ipData)
    { return elemVars.gridVariablesCache().elementCache(elemVars.eIdx_, ipData.qpIndex()); }

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
        eIdx_ = elemDisc.gridGeometry().elementMapper().index(element);
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
    const GridVariablesCache* gridVariablesCachePtr_;
    std::size_t eIdx_;
};

/*!
 * \ingroup FEDiscretization
 * \brief The local (stencil) element variables class for finite element schemes without caching
 */
template<class GVC>
class FEElementVariables<GVC, /*cachingEnabled*/false>
{
    using ThisType = FEElementVariables<GVC, /*cachingEnabled*/false>;
    using GridDiscretization = std::decay_t<decltype(std::declval<GVC>().problem().gridDiscretization())>;
    using GridView = typename GridDiscretization::GridView;

    //!< maximum number of boundary faces per element, here assumed to be the number of faces of a dim-dimensional hypercube
    static constexpr std::size_t maxNumBoundaryFaces = GridView::dimension << 1;

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
        { return *view_.ipDataCache_; }

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
    , ipDataCache_(std::make_shared<InterpolationPointDataCache>())
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
        auto elemSol = elementSolution(element, sol, elemDisc.gridDiscretization());

        // resize variables to the required size
        variables_.resize(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));

        // update variables related to localDofs
        for (const auto& localDof : localDofs(elemDisc))
            variables_[localDof.index()].update(elemSol, gridVariablesCache().problem(), elemDisc, ipData(elemDisc, localDof));

        if constexpr (InterpolationPointData::isSolDependent)
        {
            auto newIpDataCache = std::make_shared<InterpolationPointDataCache>(*ipDataCache_);
            newIpDataCache->update(gridVariablesCache().problem(), element, elemDisc, variables_);
            ipDataCache_ = std::move(newIpDataCache);
        }
        else
            ipDataCache_->update(gridVariablesCache().problem(), element, elemDisc, variables_);
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
    { return elemVars.ipDataCache_->boundaryFaceCache(ipData.boundaryFaceIndex(), ipData.qpIndex()); }

    template<Concept::QIpData IpData>
    friend const InterpolationPointData& cache(const FEElementVariables& elemVars,
                                               const IpData& ipData)
    { return elemVars.ipDataCache_->elementCache(ipData.qpIndex()); }

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
    class InterpolationPointDataCache
    {
    public:
        InterpolationPointDataCache()
        {}

        template<class Problem, class ElementDiscretization, class ElementVariables>
        void update(const Problem& problem,
                    const typename ElementDiscretization::Element& element,
                    const ElementDiscretization& elemDisc,
                    const ElementVariables& elemVars)
        {
            updateElementCache_(problem, element, elemDisc, elemVars);
        }

        // access operator
        const InterpolationPointData& elementCache(std::size_t qpIdx) const
        { return elementCache_[qpIdx]; }

        // access operator
        InterpolationPointData& elementCache(std::size_t qpIdx)
        { return elementCache_[qpIdx]; }

        // access operator
        const InterpolationPointData& boundaryFaceCache(std::size_t bfIdx, std::size_t qpIdx) const
        { return boundaryFaceCache_[bfIdx][qpIdx]; }

        // access operator
        InterpolationPointData& boundaryFaceCache(std::size_t bfIdx, std::size_t qpIdx)
        { return boundaryFaceCache_[bfIdx][qpIdx]; }

    private:
        template<class Problem, class ElementDiscretization, class ElementVariables>
        void updateElementCache_(const Problem& problem,
                                 const typename ElementDiscretization::Element& element,
                                 const ElementDiscretization& elemDisc,
                                 const ElementVariables& elemVars)
        {
            const auto elemQuadRule = Dumux::CVFE::quadratureRule(elemDisc, element);
            elementCache_.resize(std::ranges::size(elemQuadRule));
            for (const auto& qpData : elemQuadRule)
                elementCache_[qpData.ipData().qpIndex()].update(problem, element, elemDisc, elemVars, qpData.ipData());

            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                auto& bfCache = boundaryFaceCache_[boundaryFace.index()];
                const auto quadRule = Dumux::CVFE::quadratureRule(elemDisc, boundaryFace);
                bfCache.resize(std::ranges::size(quadRule));
                for (const auto& qpData : quadRule)
                    bfCache[qpData.ipData().qpIndex()].update(problem,
                                                              element,
                                                              elemDisc,
                                                              elemVars,
                                                              qpData.ipData());
            }
        }

        std::vector<InterpolationPointData> elementCache_;
        std::array<std::vector<InterpolationPointData>, maxNumBoundaryFaces> boundaryFaceCache_;
    };

    const GridVariablesCache* gridVariablesCachePtr_;
    std::vector<Variables> variables_;
    std::shared_ptr<InterpolationPointDataCache> ipDataCache_;
};

} // end namespace Dumux::Experimental

#endif
