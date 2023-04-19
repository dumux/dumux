// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondFVElementGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_FV_ELEMENT_GEOMETRY_HH

#include <type_traits>
#include <optional>

#include <dune/common/reservedvector.hh>
#include <dune/common/iteratorrange.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/facecentered/diamond/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief Element-wise grid geometry (local view)
 */
template<class GG, bool cachingEnabled>
class FaceCenteredDiamondFVElementGeometry;

/*!
 * \ingroup DiamondDiscretization
 * \brief Element-wise grid geometry (local view)
 */
template<class GG>
class FaceCenteredDiamondFVElementGeometry<GG, /*cachingEnabled*/true>
{
    using ThisType = FaceCenteredDiamondFVElementGeometry<GG, /*cachingEnabled*/true>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;

public:
    //! export type of subcontrol volume face
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GG;

    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 2*GridView::dimension;

    FaceCenteredDiamondFVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache)
    {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    { return ggCache_->scvs(eIdx_)[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    { return ggCache_->scvf(eIdx_)[scvfIdx]; }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto
    scvs(const FaceCenteredDiamondFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto
    scvfs(const FaceCenteredDiamondFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvfs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(element().type()).localBasis();
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    {
        return ggCache_->scvs(eIdx_).size();
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    {
        return ggCache_->scvfs(eIdx_).size();
    }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return ggCache_->hasBoundaryScvf(eIdx_); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    FaceCenteredDiamondFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    void bind(const Element& element) &
    {
        this->bindElement(element);
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    FaceCenteredDiamondFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Bind only element-local
    void bindElement(const Element& element) &
    {
        element_ = element;
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return ggCache_->gridGeometry(); }

    //! The bound element index
    std::size_t elementIndex() const
    { return eIdx_; }

    //! Geometry of a sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        return {
            SubControlVolume::Traits::geometryType(geo.type()),
            GeometryHelper(geo).getScvCorners(scv.indexInElement())
        };
    }

    //! Geometry of a sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        if (scvf.boundary())
        {
            // use the information that each boundary scvf corresponds to one scv constructed around the same facet
            const auto localFacetIndex = scvf.insideScvIdx();
            return {
                referenceElement(geo).type(localFacetIndex, 1),
                GeometryHelper(geo).getBoundaryScvfCorners(localFacetIndex)
            };
        }
        else
        {
            return {
                SubControlVolumeFace::Traits::interiorGeometryType(geo.type()),
                GeometryHelper(geo).getScvfCorners(scvf.index())
            };
        }
    }

private:
    std::optional<Element> element_;
    GridIndexType eIdx_;
    const GGCache* ggCache_;
};

} // end namespace Dumux

#endif
