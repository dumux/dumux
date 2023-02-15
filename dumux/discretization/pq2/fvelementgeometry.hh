// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief Base class for the local finite volume geometry for the pq2 method
 *        This builds up the sub control volumes and sub control volume faces
 *        for an element.
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ2_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

#include <dumux/discretization/pq2/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ2Discretization
 * \brief Base class for the finite volume geometry vector for pq2 models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class PQ2FVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class PQ2FVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = PQ2GeometryHelper<GridView, typename GG::SubControlVolume, typename GG::SubControlVolumeFace>;
public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (here for cubes)
    // ToDo get this from GG
    static constexpr std::size_t maxNumElementScvs = (1<<dim)+ dim*(1<<(dim-1));

    //! Constructor
    PQ2FVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache)
    {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return ggCache_->scvs(eIdx_)[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return ggCache_->scvfs(eIdx_)[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const PQ2FVElementGeometry& fvGeometry)
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
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const PQ2FVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvfs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(element_->type()).localBasis();
    }

    //! Get a local finite element basis
    const auto& feLocalCoefficients() const
    {
        return gridGeometry().feCache().get(element_->type()).localCoefficients();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return ggCache_->scvs(eIdx_).size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return ggCache_->scvfs(eIdx_).size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    PQ2FVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! this function is for compatibility reasons with cc methods
    //! The pq2 stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    PQ2FVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element) &
    {
        element_ = element;
        // cache element index
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

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return ggCache_->hasBoundaryScvf(eIdx_); }

    //! The bound element index
    std::size_t elementIndex() const
    { return eIdx_; }

    //! Geometry of a sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        const GeometryHelper helper(geo);
        const auto idx =  helper.localKeyToLocalScvIndex(this->feLocalCoefficients().localKey(scv.indexInElement()));
        return {
            helper.getScvGeometryType(idx),
            helper.getScvCorners(idx)
        };
    }

    //! Geometry of a sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        if (scvf.boundary())
        {
            GeometryHelper helper(geo);
            const auto localScvfIdx = scvf.index() - helper.numInteriorScvf();
            const auto [localFacetIndex, isScvfLocalIdx]
                = ggCache_->scvfBoundaryGeometryKeys(eIdx_)[localScvfIdx];
            return {
                helper.getBoundaryScvfGeometryType(localFacetIndex, localScvfIdx),
                helper.getBoundaryScvfCorners(localFacetIndex, isScvfLocalIdx)
            };
        }
        else
        {
            GeometryHelper helper(geo);
            return {
                helper.getInteriorScvfGeometryType(scvf.index()),
                helper.getScvfCorners(scvf.index())
            };
        }
    }

private:
    const GGCache* ggCache_;
    GridIndexType eIdx_;

    std::optional<Element> element_;
};

} // end namespace Dumux

#endif
