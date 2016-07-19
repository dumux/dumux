// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for the local finite volume geometry for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for an element.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BOX_FV_ELEMENT_GEOMETRY_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/implicit/box/properties.hh>

namespace Dumux
{

//! forward declaration of the global finite volume geometry
template<class TypeTag, bool EnableGlobalFVGeometryCache>
class BoxGlobalFVGeometry;

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableGlobalFVGeometryCache>
class BoxFVElementGeometry
{};

//! specialization in case the FVElementGeometries are stored
template<class TypeTag>
class BoxFVElementGeometry<TypeTag, true>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalFVGeometry = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    using ScvIterator = ScvIterator<SubControlVolume, std::vector<IndexType>, ThisType>;
    using ScvfIterator = ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType>;

public:
    //! Constructor
    BoxFVElementGeometry(const GlobalFVGeometry& globalFvGeometry)
    : globalFvGeometryPtr_(&globalFvGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return globalFvGeometry().scvs(eIdx_)[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        return globalFvGeometry().scvfs(eIdx_)[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.globalFvGeometry();
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
        const auto& g = fvGeometry.globalFvGeometry();
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvfs(fvGeometry.eIdx_).begin(), g.scvfs(fvGeometry.eIdx_).end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return globalFvGeometry().feCache().get(elementPtr_->geometry().type()).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return globalFvGeometry().scvs(eIdx_).size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return globalFvGeometry().scvfs(eIdx_).size();
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
        elementPtr_ = &element;
        eIdx_ = globalFvGeometry().problem_().elementMapper().index(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GlobalFVGeometry& globalFvGeometry() const
    { return *globalFvGeometryPtr_; }

private:
    const Element* elementPtr_;
    const GlobalFVGeometry* globalFvGeometryPtr_;

    IndexType eIdx_;
};

//! specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class BoxFVElementGeometry<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GlobalFVGeometry = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    using GeometryHelper = BoxGeometryHelper<GridView, dim>;

public:
    //! Constructor
    BoxFVElementGeometry(const GlobalFVGeometry& globalFvGeometry)
    : globalFvGeometryPtr_(&globalFvGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
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
        return globalFvGeometry().feCache().get(elementPtr_->geometry().type()).localBasis();
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
        elementPtr_ = &element;
        eIdx_ = globalFvGeometry().problem_().elementMapper().index(element);
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GlobalFVGeometry& globalFvGeometry() const
    { return *globalFvGeometryPtr_; }

private:

    void makeElementGeometries(const Element& element)
    {
        auto eIdx = globalFvGeometry().problem_().elementMapper().index(element);

        // get the element geometry
        auto elementGeometry = element.geometry();
        const auto& referenceElement = ReferenceElements::general(elementGeometry.type());

        // get the sub control volume geometries of this element
        GeometryHelper geometryHelper(elementGeometry);

        // construct the sub control volumes
        for (unsigned int scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
        {
            // get asssociated dof index
            auto dofIdxGlobal = globalFvGeometry().problem_().vertexMapper().subIndex(element, scvLocalIdx, dim);

            // get the corners and the volume of the scv
            auto scvCorners = geometryHelper.getScvCorners(scvLocalIdx);
            auto volume = geometryHelper.volume(scvCorners);

            // add scv to the local container
            scvs_.emplace_back(std::move(scvCorners),
                               volume,
                               scvLocalIdx,
                               eIdx,
                               dofIdxGlobal);
        }

        // construct the sub control volume faces
        unsigned int scvfLocalIdx = 0;
        for (; scvfLocalIdx < element.subEntities(1); ++scvfLocalIdx)
        {
            // find the local scv indices this scvf is connsected to
            std::vector<IndexType> localScvIndices({static_cast<IndexType>(referenceElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                    static_cast<IndexType>(referenceElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

            // get the corner points, the area, and the unit normal vector of the scv face
            auto scvfCorners = geometryHelper.getScvfCorners(scvfLocalIdx);
            auto area = geometryHelper.area(scvfCorners);
            auto normal = geometryHelper.normal(elementGeometry, scvfCorners);
            const auto v = elementGeometry.corner(localScvIndices[1]) - elementGeometry.corner(localScvIndices[0]);
            const auto s = v*normal;
            if (std::signbit(s))
                normal *= -1;

            scvfs_.emplace_back(std::move(scvfCorners),
                                normal,
                                area,
                                scvfLocalIdx,
                                localScvIndices,
                                false);
        }

        // construct the sub control volume faces on the domain boundary
        for (const auto& intersection : intersections(globalFvGeometry().gridView(), element))
        {
            if (intersection.boundary())
            {
                auto isGeometry = intersection.geometry();

                for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                {
                    // find the scv this scvf is connected to
                    std::vector<IndexType> localScvIndices =
                       {static_cast<IndexType>(referenceElement.subEntity(intersection.indexInInside(), 1,
                                                                          isScvfLocalIdx, dim))};

                    // get the corner points and the area of the boundary scv face
                    auto scvfCorners = geometryHelper.getBoundaryScvfCorners(isGeometry, isScvfLocalIdx);
                    auto area = geometryHelper.area(scvfCorners);

                    scvfs_.emplace_back(std::move(scvfCorners),
                                        intersection.centerUnitOuterNormal(),
                                        area,
                                        scvfLocalIdx,
                                        localScvIndices,
                                        true);

                    // increment local counter
                    scvfLocalIdx++;
                }
            }
        }
    }

    //! The bound element
    const Element* elementPtr_;
    IndexType eIdx_;

    //! The global geometry this is a restriction of
    const GlobalFVGeometry* globalFvGeometryPtr_;

    //! vectors to store the geometries locally after binding an element
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
};

} // end namespace

#endif
