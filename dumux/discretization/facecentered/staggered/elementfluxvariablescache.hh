// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_ELEMENT_FLUXVARSCACHE_HH

#include <cstddef>
#include <vector>
#include <utility>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief The flux variables caches for an element
 * \note The class is specialized for a version with and without caching
 * If grid caching is enabled the flux caches are stored for the whole gridview in the corresponding
 * GridFluxVariablesCache which is memory intensive but faster. For caching disabled the
 * flux caches are locally computed for each element whenever needed.
 */
template<class GFVC, bool cachingEnabled>
class FaceCenteredStaggeredElementFluxVariablesCache;

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief The flux variables caches for an element with caching enabled
 */
template<class GFVC>
class FaceCenteredStaggeredElementFluxVariablesCache<GFVC, true>
{
public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    FaceCenteredStaggeredElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    FaceCenteredStaggeredElementFluxVariablesCache bind(const typename FVElementGeometry::Element& element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    // Function is called by the BoxLocalJacobian prior to flux calculations on the element.
    // We assume the FVGeometries to be bound at this point
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::Element& element,
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
    FaceCenteredStaggeredElementFluxVariablesCache bindElement(const typename FVElementGeometry::Element& element,
                                                               const FVElementGeometry& fvGeometry,
                                                               const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) &
    {
        eIdx_ = fvGeometry.gridGeometry().elementMapper().index(element);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    FaceCenteredStaggeredElementFluxVariablesCache bindScvf(const typename FVElementGeometry::Element& element,
                                                            const FVElementGeometry& fvGeometry,
                                                            const ElementVolumeVariables& elemVolVars,
                                                            const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) &
    {
        bindElement(element, fvGeometry, elemVolVars);
    }

    // access operator
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return gridFluxVarsCache()[scvf]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
    std::size_t eIdx_; //!< currently bound element
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief The flux variables caches for an element with caching disabled
 */
template<class GFVC>
class FaceCenteredStaggeredElementFluxVariablesCache<GFVC, false>
{
public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    FaceCenteredStaggeredElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    FaceCenteredStaggeredElementFluxVariablesCache bind(const typename FVElementGeometry::Element& element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    // Function is called by the BoxLocalJacobian prior to flux calculations on the element.
    // We assume the FVGeometries to be bound at this point
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::Element& element,
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
    FaceCenteredStaggeredElementFluxVariablesCache bindElement(const typename FVElementGeometry::Element& element,
                                                               const FVElementGeometry& fvGeometry,
                                                               const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) &
    {
        // temporary resizing of the cache
        fluxVarsCache_.resize(fvGeometry.numScvf());
        // for (auto&& scvf : scvfs(fvGeometry))
            // (*this)[scvf].update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf); TODO
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    FaceCenteredStaggeredElementFluxVariablesCache bindScvf(const typename FVElementGeometry::Element& element,
                                                            const FVElementGeometry& fvGeometry,
                                                            const ElementVolumeVariables& elemVolVars,
                                                            const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) &
    {
        fluxVarsCache_.resize(fvGeometry.numScvf());
        // (*this)[scvf].update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf); TODO
    }

    // access operator
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    // access operator
    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
};

} // end namespace Dumux

#endif
