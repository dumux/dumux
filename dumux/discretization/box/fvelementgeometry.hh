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
 * \ingroup BoxDiscretization
 * \brief Base class for the local finite volume geometry for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for an element.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BOX_FV_ELEMENT_GEOMETRY_HH

#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class BoxFVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class BoxFVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementScvs = (1<<dim);

    //! Constructor
    BoxFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return gridGeometry().scvs(eIdx_)[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return gridGeometry().scvfs(eIdx_)[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvs(fvGeometry.eIdx_).begin(), g.scvs(fvGeometry.eIdx_).end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvfs(fvGeometry.eIdx_).begin(), g.scvfs(fvGeometry.eIdx_).end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(elemGeometryType_).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return gridGeometry().scvs(eIdx_).size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return gridGeometry().scvfs(eIdx_).size();
    }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        elemGeometryType_ = element.type();
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(eIdx_); }

private:
    Dune::GeometryType elemGeometryType_;
    const GridGeometry* gridGeometryPtr_;

    GridIndexType eIdx_;
};

//! specialization in case the FVElementGeometries are not stored
template<class GG>
class BoxFVElementGeometry<GG, false>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;

    using GeometryHelper = BoxGeometryHelper<GridView, dim,
                                             typename GG::SubControlVolume,
                                             typename GG::SubControlVolumeFace>;
public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementScvs = (1<<dim);

    //! Constructor
    BoxFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvs_.begin(), fvGeometry.scvs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(elemGeometryType_).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        eIdx_ = gridGeometry().elementMapper().index(element);
        makeElementGeometries(element);
    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

private:

    void makeElementGeometries(const Element& element)
    {
        hasBoundaryScvf_ = false;

        // get the element geometry
        auto elementGeometry = element.geometry();
        elemGeometryType_ = elementGeometry.type();
        const auto refElement = referenceElement(elementGeometry);

        // get the sub control volume geometries of this element
        GeometryHelper geometryHelper(elementGeometry);

        // construct the sub control volumes
        scvs_.resize(elementGeometry.corners());
        for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
        {
            // get asssociated dof index
            const auto dofIdxGlobal = gridGeometry().vertexMapper().subIndex(element, scvLocalIdx, dim);

            // add scv to the local container
            scvs_[scvLocalIdx] = SubControlVolume(geometryHelper,
                                                  scvLocalIdx,
                                                  eIdx_,
                                                  dofIdxGlobal);
        }

        // construct the sub control volume faces
        const auto numInnerScvf = (dim==1) ? 1 : element.subEntities(dim-1);
        scvfs_.resize(numInnerScvf);

        LocalIndexType scvfLocalIdx = 0;
        for (; scvfLocalIdx < numInnerScvf; ++scvfLocalIdx)
        {
            // find the local scv indices this scvf is connected to
            std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                         static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

            scvfs_[scvfLocalIdx] = SubControlVolumeFace(geometryHelper,
                                                        element,
                                                        elementGeometry,
                                                        scvfLocalIdx,
                                                        std::move(localScvIndices),
                                                        false);
        }

        // construct the sub control volume faces on the domain boundary
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (intersection.boundary() && !intersection.neighbor())
            {
                const auto isGeometry = intersection.geometry();
                hasBoundaryScvf_ = true;

                for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                {
                    // find the scv this scvf is connected to
                    const LocalIndexType insideScvIdx = static_cast<LocalIndexType>(refElement.subEntity(intersection.indexInInside(), 1, isScvfLocalIdx, dim));
                    std::vector<LocalIndexType> localScvIndices = {insideScvIdx, insideScvIdx};

                    scvfs_.emplace_back(geometryHelper,
                                        intersection,
                                        isGeometry,
                                        isScvfLocalIdx,
                                        scvfLocalIdx,
                                        std::move(localScvIndices),
                                        true);

                    // increment local counter
                    scvfLocalIdx++;
                }
            }
        }
    }

    //! The bound element
    Dune::GeometryType elemGeometryType_;
    GridIndexType eIdx_;

    //! The global geometry this is a restriction of
    const GridGeometry* gridGeometryPtr_;

    //! vectors to store the geometries locally after binding an element
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    bool hasBoundaryScvf_ = false;
};

} // end namespace Dumux

#endif
