// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEDiscretization
 * \brief Element data cache
 */
#ifndef DUMUX_DISCRETIZATION_FE_ELEMENT_DATACACHE_HH
#define DUMUX_DISCRETIZATION_FE_ELEMENT_DATACACHE_HH

#include <array>
#include <cstddef>
#include <optional>
#include <ranges>
#include <type_traits>
#include <vector>
#include <utility>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup FEDiscretization
 * \brief The data caches for an element when using FE discretizations
 * \note The class is specialized for caching enabled/disabled and for different quadrature rules
 *       If grid caching is enabled the data caches are stored for the whole gridview in the corresponding
 *       GridDataCache which is memory intensive but faster. For caching disabled the
 *       data caches are locally computed for each element whenever needed.
 */
template<class GDC, bool cachingEnabled>
class FEElementDataCache;

/*!
 * \ingroup FEDiscretization
 * \brief The data caches for an element with caching enabled (for general quadrature rules)
 */
template<class GDC>
class FEElementDataCache<GDC, true>
{
public:
    //! export the type of the grid data cache
    using GridDataCache = GDC;

    //! export the type of the local data cache
    using DataCache = typename GDC::DataCache;

    FEElementDataCache(const GridDataCache& global)
    : gridDataCachePtr_(&global) {}

    // Function is called prior to residual calculations on the element.
    // We assume the ElementGeometryView to be bound at this point
    template<class ElementGeometryView, class ElementVariables>
    void bind(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
              const ElementGeometryView& eGeometryView,
              const ElementVariables& elemVars) &
    { bindElement(element, eGeometryView, elemVars); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementGeometryView, class ElementVariables>
    FEElementDataCache bind(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const ElementGeometryView& eGeometryView,
                            const ElementVariables& elemVars) &&
    {
        this->bind(element, eGeometryView, elemVars);
        return std::move(*this);
    }

    template<class ElementGeometryView, class ElementVariables>
    void bindElement(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const ElementGeometryView& eGeometryView,
                     const ElementVariables& elemVars) &
    { eIdx_ = eGeometryView.gridGeometry().elementMapper().index(element); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementGeometryView, class ElementVariables>
    FEElementDataCache bindElement(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
                                   const ElementGeometryView& eGeometryView,
                                   const ElementVariables& elemVars) &&
    {
        this->bindElement(element, eGeometryView, elemVars);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class ElementGeometryView, class ElementVariables>
    void update(const typename ElementGeometryView::Element& element,
                const ElementGeometryView& eGeometryView,
                const ElementVariables& elemVars) {}

    // access cache for given interpolation point data of a boundary intersection quadrature point
    template<Concept::IntersectionQpIpData IpData>
    const DataCache& operator [](const IpData& ipData) const
    {
        return gridDataCache().boundaryIntersectionCache(eIdx_, ipData.intersectionIndex(), ipData.qpIndex());
    }

    // access cache for given interpolation point data of an element quadrature point
    template<Concept::QIpData IpData>
    const DataCache& operator [](const IpData& ipData) const
    {
        return gridDataCache().elementCache(eIdx_, ipData.qpIndex());
    }

    //! The global object we are a restriction of
    const GridDataCache& gridDataCache() const
    {  return *gridDataCachePtr_; }

private:
    const GridDataCache* gridDataCachePtr_;
    std::size_t eIdx_; //!< currently bound element
};

/*!
 * \ingroup FEDiscretization
 * \brief The data caches for an element with caching disabled (general quadrature)
 */
template<class GDC>
class FEElementDataCache<GDC, false>
{
    using GridGeometry = std::decay_t<decltype(std::declval<GDC>().problem().gridGeometry())>;
    using GridView = typename GridGeometry::GridView;

    //!< maximum number of boundary intersections per element, here assumed to be the number of faces of a dim-dimensional hypercube
    static constexpr std::size_t maxNumIntersections = GridView::dimension << 1;

public:
    //! export the type of the grid data cache
    using GridDataCache = GDC;

    //! export the type of the local data cache
    using DataCache = typename GDC::DataCache;

    FEElementDataCache(const GridDataCache& global)
    : gridDataCachePtr_(&global) {}

    // Function is called prior to residual calculations on the element
    // We assume the ElementGeometryView to be bound at this point
    template<class ElementGeometryView, class ElementVariables>
    void bind(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
              const ElementGeometryView& eGeometryView,
              const ElementVariables& elemVars) &
    {
        bindElement(element, eGeometryView, elemVars);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementGeometryView, class ElementVariables>
    FEElementDataCache bind(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const ElementGeometryView& eGeometryView,
                            const ElementVariables& elemVars) &&
    {
        this->bind(element, eGeometryView, elemVars);
        return std::move(*this);
    }

    template<class ElementGeometryView, class ElementVariables>
    void bindElement(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const ElementGeometryView& eGeometryView,
                     const ElementVariables& elemVars) &
    {
        // Resize element cache based on element quadrature rule
        const auto elemQuadRule = CVFE::quadratureRule(eGeometryView, element);
        elementCache_.resize(std::ranges::size(elemQuadRule));
        for (const auto& qpData : elemQuadRule)
            elementCache_[qpData.ipData().qpIndex()].update(gridDataCache().problem(), element, eGeometryView, elemVars, qpData.ipData());

        // Set boundary intersection cache entries and reset non-boundary ones
        for (const auto& intersection : intersections(eGeometryView.gridGeometry().gridView(), element))
        {
            const auto iIdx = intersection.indexInInside();
            if (intersection.boundary())
            {
                auto& boundaryCache = boundaryIntersectionCache_[iIdx];
                boundaryCache.emplace();
                const auto quadRule = CVFE::quadratureRule(eGeometryView, intersection);
                boundaryCache->resize(quadRule.size());
                for (const auto& qpData : quadRule)
                    (*boundaryCache)[qpData.ipData().qpIndex()].update(gridDataCache().problem(),
                                                                       element,
                                                                       eGeometryView,
                                                                       elemVars,
                                                                       qpData.ipData());
            }
            else
                boundaryIntersectionCache_[iIdx].reset();
        }
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class ElementGeometryView, class ElementVariables>
    FEElementDataCache bindElement(const typename ElementGeometryView::GridGeometry::GridView::template Codim<0>::Entity& element,
                                   const ElementGeometryView& eGeometryView,
                                   const ElementVariables& elemVars) &&
    {
        this->bindElement(element, eGeometryView, elemVars);
        return std::move(*this);
    }

    /*!
     * \brief Update the caches if the variables have changed and the cache is solution-dependent
     * \note Results in undefined behaviour if called before bind() or with a different element
     */
    template<class ElementGeometryView, class ElementVariables>
    void update(const typename ElementGeometryView::Element& element,
                const ElementGeometryView& eGeometryView,
                const ElementVariables& elemVars)
    {
        if constexpr (DataCache::isSolDependent)
            bindElement(element, eGeometryView, elemVars);
    }

    // access cache for a given interpolation point data of an intersection quadrature point
    template<Concept::IntersectionQpIpData IpData>
    const DataCache& operator [](const IpData& ipData) const
    { return (*boundaryIntersectionCache_[ipData.intersectionIndex()])[ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an intersection quadrature point
    template<Concept::IntersectionQpIpData IpData>
    DataCache& operator [](const IpData& ipData)
    { return (*boundaryIntersectionCache_[ipData.intersectionIndex()])[ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an element quadrature point
    template<Concept::QIpData IpData>
    const DataCache& operator [](const IpData& ipData) const
    { return elementCache_[ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an element quadrature point
    template<Concept::QIpData IpData>
    DataCache& operator [](const IpData& ipData)
    { return elementCache_[ipData.qpIndex()]; }

    //! The global object we are a restriction of
    const GridDataCache& gridDataCache() const
    {  return *gridDataCachePtr_; }

private:
    const GridDataCache* gridDataCachePtr_;
    std::vector<DataCache> elementCache_;
    std::array<std::optional<std::vector<DataCache>>, maxNumIntersections> boundaryIntersectionCache_;
};

} // end namespace Dumux

#endif
