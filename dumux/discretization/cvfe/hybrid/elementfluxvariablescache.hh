// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_HYBRID_CVFE_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_HYBRID_CVFE_ELEMENT_FLUXVARSCACHE_HH

#include <array>
#include <cassert>
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
 * \ingroup CVFEDiscretization
 * \brief The flux variables caches for an element when using hybrid CVFE discretizations
 * \note The class is specialized for caching enabled/disabled and for different quadrature rules
 * If grid caching is enabled the flux caches are stored for the whole gridview in the corresponding
 * GridFluxVariablesCache which is memory intensive but faster. For caching disabled the
 * flux caches are locally computed for each element whenever needed.
 */
template<class GFVC, bool cachingEnabled>
class HybridCVFEElementFluxVariablesCache;

/*!
 * \ingroup CVFEDiscretization
 * \brief The flux variables caches for an element with caching enabled (for general quadrature rules)
 */
template<class GFVC>
class HybridCVFEElementFluxVariablesCache<GFVC, true>
{
public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    HybridCVFEElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    // Function is called prior to flux calculations on the element.
    // We assume the FVGeometries to be bound at this point
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) &
    { bindElement(element, fvGeometry, elemVolVars); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    HybridCVFEElementFluxVariablesCache bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) &
    { eIdx_ = fvGeometry.gridGeometry().elementMapper().index(element); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    HybridCVFEElementFluxVariablesCache bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) &
    { bindElement(element, fvGeometry, elemVolVars); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    HybridCVFEElementFluxVariablesCache bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}

    // access cache for given interpolation point data of an scvf quadrature point
    template<Concept::ScvfQpIpData IpData>
    const FluxVariablesCache& operator [](const IpData& ipData) const
    {
        return gridFluxVarsCache().cache(eIdx_, ipData.scvfIndex(), ipData.qpIndex());
    }

    // access cache for given interpolation point data of a boundary intersection quadrature point
    template<Concept::IntersectionQpIpData IpData>
    const FluxVariablesCache& operator [](const IpData& ipData) const
    {
        return gridFluxVarsCache().boundaryIntersectionCache(eIdx_, ipData.intersectionIndex(), ipData.qpIndex());
    }

    // access cache for given interpolation point data of an element quadrature point
    template<Concept::QIpData IpData>
    const FluxVariablesCache& operator [](const IpData& ipData) const
    {
        return gridFluxVarsCache().elementCache(eIdx_, ipData.qpIndex());
    }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
    std::size_t eIdx_; //!< currently bound element
};

/*!
 * \ingroup CVFEDiscretization
 * \brief The flux variables caches for an element with caching disabled (general quadrature)
 */
template<class GFVC>
class HybridCVFEElementFluxVariablesCache<GFVC, false>
{
    using GridGeometry = std::decay_t<decltype(std::declval<GFVC>().problem().gridGeometry())>;
    using GridView = typename GridGeometry::GridView;

    //!< maximum number of boundary intersections per element, here assumed to be the number of faces of a dim-dimensional hypercube
    static constexpr std::size_t maxNumIntersections = GridView::dimension << 1;

public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    HybridCVFEElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    // Function is called prior to flux calculations on the element
    // We assume the FVGeometries to be bound at this point
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) &
    {
        bindElement(element, fvGeometry, elemVolVars);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    HybridCVFEElementFluxVariablesCache bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) &
    {
        // temporary resizing of the cache
        fluxVarsCache_.resize(fvGeometry.numScvf());
        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto quadRule = CVFE::quadratureRule(fvGeometry, scvf);
            fluxVarsCache_[scvf.index()].resize(std::ranges::size(quadRule));
            for (const auto& qpData : quadRule)
                fluxVarsCache_[scvf.index()][qpData.ipData().qpIndex()].update(
                    gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, qpData.ipData()
                );
        }

        // Resize element cache based on element quadrature rule
        const auto elemQuadRule = CVFE::quadratureRule(fvGeometry, element);
        elementCache_.resize(std::ranges::size(elemQuadRule));
        for (const auto& qpData : elemQuadRule)
            elementCache_[qpData.ipData().qpIndex()].update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, qpData.ipData());

        // Set boundary intersection cache entries and reset non-boundary ones
        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            const auto iIdx = intersection.indexInInside();
            if (intersection.boundary())
            {
                auto& boundaryCache = boundaryIntersectionCache_[iIdx];
                boundaryCache.emplace();
                const auto quadRule = CVFE::quadratureRule(fvGeometry, intersection);
                boundaryCache->resize(quadRule.size());
                for (const auto& qpData : quadRule)
                    (*boundaryCache)[qpData.ipData().qpIndex()].update(gridFluxVarsCache().problem(),
                                                                       element,
                                                                       fvGeometry,
                                                                       elemVolVars,
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
    template<class FVElementGeometry, class ElementVolumeVariables>
    HybridCVFEElementFluxVariablesCache bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) &
    {
        fluxVarsCache_.resize(fvGeometry.numScvf());
        const auto quadRule = CVFE::quadratureRule(fvGeometry, scvf);
        fluxVarsCache_[scvf.index()].resize(std::ranges::size(quadRule));
        for (const auto& qpData : quadRule)
            fluxVarsCache_[scvf.index()][qpData.ipData().qpIndex()].update(
                gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, qpData.ipData()
            );
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    HybridCVFEElementFluxVariablesCache bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    /*!
     * \brief Update the caches if the volume variables have changed and the cache is solution-dependent
     * \note Results in undefined behaviour if called before bind() or with a different element
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        if constexpr (FluxVariablesCache::isSolDependent)
            bindElement(element, fvGeometry, elemVolVars);
    }

    // access cache for a given interpolation point data of an scvf quadrature point
    template<Concept::ScvfQpIpData IpData>
    const FluxVariablesCache& operator [](const IpData& ipData) const
    { return fluxVarsCache_[ipData.scvfIndex()][ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an scvf quadrature point
    template<Concept::ScvfQpIpData IpData>
    FluxVariablesCache& operator [](const IpData& ipData)
    { return fluxVarsCache_[ipData.scvfIndex()][ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an intersection quadrature point
    template<Concept::IntersectionQpIpData IpData>
    const FluxVariablesCache& operator [](const IpData& ipData) const
    { return (*boundaryIntersectionCache_[ipData.intersectionIndex()])[ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an intersection quadrature point
    template<Concept::IntersectionQpIpData IpData>
    FluxVariablesCache& operator [](const IpData& ipData)
    { return (*boundaryIntersectionCache_[ipData.intersectionIndex()])[ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an element quadrature point
    template<Concept::QIpData IpData>
    const FluxVariablesCache& operator [](const IpData& ipData) const
    { return elementCache_[ipData.qpIndex()]; }

    // access cache for a given interpolation point data of an element quadrature point
    template<Concept::QIpData IpData>
    FluxVariablesCache& operator [](const IpData& ipData)
    { return elementCache_[ipData.qpIndex()]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_;
    std::vector<FluxVariablesCache> elementCache_;
    std::array<std::optional<std::vector<FluxVariablesCache>>, maxNumIntersections> boundaryIntersectionCache_;
};

} // end namespace Dumux

#endif
